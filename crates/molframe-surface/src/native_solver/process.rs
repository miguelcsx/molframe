//! Provider-specific command protocols and bounded process lifetime.
use super::{NativeAtoms, NativeSurfaceError, NativeSurfaceOptions, NativeSurfaceProvider};
use molframe_core::execution::{ExecutionContext, TemporaryWorkspace};
use std::{
    fs::File,
    io::{BufWriter, Read, Seek, SeekFrom, Write},
    path::Path,
    process::{Child, Command, Stdio},
    time::{Duration, Instant},
};

pub(super) fn write_input(
    directory: &Path,
    provider: NativeSurfaceProvider,
    atoms: NativeAtoms<'_>,
    context: &ExecutionContext,
) -> Result<(), NativeSurfaceError> {
    let mut output = BufWriter::new(File::create(directory.join("atoms"))?);
    for (index, &row) in atoms.atoms.rows.iter().enumerate() {
        if index % 1024 == 0 {
            context.cancellation().check()?;
        }
        let [x, y, z] = atoms.position(row).ok_or(NativeSurfaceError::Invalid(
            "invalid selected surface coordinate",
        ))?;
        let radius = atoms.atoms.radii[row as usize];
        match provider {
            NativeSurfaceProvider::Surf => writeln!(output, "{index} {radius} {x} {y} {z}")?,
            NativeSurfaceProvider::Msms => writeln!(output, "{x} {y} {z} {radius}")?,
        }
    }
    output.flush()?;
    Ok(())
}

pub(super) fn run(
    executable: &Path,
    provider: NativeSurfaceProvider,
    options: NativeSurfaceOptions,
    workspace: &TemporaryWorkspace,
    context: &ExecutionContext,
) -> Result<(), NativeSurfaceError> {
    let log = File::create(workspace.path().join("provider.log"))?;
    let mut command = Command::new(executable);
    command
        .current_dir(workspace.path())
        .stdin(Stdio::null())
        .stdout(Stdio::from(log.try_clone()?))
        .stderr(Stdio::from(log));
    match provider {
        NativeSurfaceProvider::Surf => {
            command
                .args(["-W", "1", "-R"])
                .arg(options.probe_radius.to_string())
                .arg("atoms");
        }
        NativeSurfaceProvider::Msms => {
            command
                .args(["-if", "atoms", "-of", "surface", "-probe_radius"])
                .arg(options.probe_radius.to_string())
                .arg("-density")
                .arg(options.density.to_string())
                .args(["-no_area", "-no_header", "-no_rest"]);
            if options.all_components {
                command.arg("-all_components");
            }
        }
    }
    let mut child = Process(command.spawn()?);
    let start = Instant::now();
    loop {
        context.cancellation().check()?;
        workspace.check_size()?;
        if let Some(status) = child.0.try_wait()? {
            return if status.success() {
                Ok(())
            } else {
                Err(NativeSurfaceError::Failed {
                    status,
                    diagnostic: diagnostic_tail(workspace.path())?,
                })
            };
        }
        if start.elapsed() >= options.timeout {
            return Err(NativeSurfaceError::Timeout);
        }
        std::thread::sleep(Duration::from_millis(20));
    }
}

fn diagnostic_tail(directory: &Path) -> Result<String, std::io::Error> {
    // A failed solver can produce a large log. Keep its final explanation
    // without retaining its private files or reading the complete stream.
    let mut input = File::open(directory.join("provider.log"))?;
    let size = input.metadata()?.len();
    input.seek(SeekFrom::Start(size.saturating_sub(8192)))?;
    let mut bytes = Vec::with_capacity(8192);
    input.take(8192).read_to_end(&mut bytes)?;
    Ok(String::from_utf8_lossy(&bytes).trim().to_owned())
}

struct Process(Child);
impl Drop for Process {
    fn drop(&mut self) {
        // Every early return terminates and reaps before the workspace is removed.
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}
