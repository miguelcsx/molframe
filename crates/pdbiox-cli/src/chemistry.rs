//! Explicit chemistry-provider loading for operations that require CCD.

use crate::exit::Exit;
use crate::report::Context;
use std::path::Path;

pub(crate) struct CcdSource {
    pub(crate) path: &'static Path,
    pub(crate) version: &'static str,
    pub(crate) context: Context,
}

pub(crate) fn resolve_ccd(
    path: Option<&Path>,
    version: Option<&str>,
    mut context: Context,
) -> Result<CcdSource, Exit> {
    let (path, version) = match (path, version) {
        (Some(path), Some(version)) => (
            Box::leak(Box::new(path.to_path_buf())) as &'static Path,
            Box::leak(version.to_owned().into_boxed_str()) as &'static str,
        ),
        (None, None) => {
            if let (Some(path), Some(version)) = (context.ccd, context.ccd_version) {
                (path, version)
            } else {
                eprintln!(
                    "chemistry requires --ccd and --ccd-version, or chem.ccd_cache and \
                     chem.ccd_version in configuration"
                );
                return Err(Exit::Usage);
            }
        }
        _ => {
            eprintln!("--ccd and --ccd-version must be supplied together");
            return Err(Exit::Usage);
        }
    };
    context.ccd = Some(path);
    context.ccd_version = Some(version);
    Ok(CcdSource {
        path,
        version,
        context,
    })
}

pub(crate) fn with_ccd(
    path: Option<&Path>,
    version: Option<&str>,
    context: Context,
    run: impl FnOnce(&Path, &str, Context) -> Exit,
) -> Exit {
    match resolve_ccd(path, version, context) {
        Ok(source) => run(source.path, source.version, source.context),
        Err(exit) => exit,
    }
}

pub(crate) fn component(path: &Path, version: &str, id: &str, context: Context) -> Exit {
    use pdbiox::ComponentProvider as _;
    let provider = match load_ccd(path, version, context) {
        Ok(provider) => provider,
        Err(exit) => return exit,
    };
    let component = match provider.get(id) {
        Ok(Some(component)) => component,
        Ok(None) => {
            eprintln!("component {id:?} is absent from {}", path.display());
            return Exit::Indeterminate;
        }
        Err(finding) => {
            context.findings(&[finding], &path.display().to_string());
            return Exit::Consistency;
        }
    };
    if context.is_json() {
        let mut json = crate::report::Json::new();
        json.text("id", &component.id)
            .text("name", &component.name)
            .text("kind", component_kind(component.kind))
            .number("atoms", component.atoms.len())
            .number("bonds", component.bonds.len())
            .text("dictionary_version", version);
        if let Some(parent) = &component.parent {
            json.text("parent", parent);
        }
        context.result(&json.finish());
    } else {
        context.result(&format!(
            "{}\t{}\t{}\t{} atoms\t{} bonds\tCCD {}",
            component.id,
            component.name,
            component_kind(component.kind),
            component.atoms.len(),
            component.bonds.len(),
            version
        ));
    }
    Exit::Success
}

const fn component_kind(kind: pdbiox::ComponentKind) -> &'static str {
    match kind {
        pdbiox::ComponentKind::AminoAcid => "amino-acid",
        pdbiox::ComponentKind::Nucleotide => "nucleotide",
        pdbiox::ComponentKind::Saccharide => "saccharide",
        pdbiox::ComponentKind::Lipid => "lipid",
        pdbiox::ComponentKind::NonPolymer => "non-polymer",
        pdbiox::ComponentKind::Solvent => "solvent",
        pdbiox::ComponentKind::Ion => "ion",
        pdbiox::ComponentKind::Unknown => "unknown",
    }
}

pub(crate) fn load_ccd(
    path: &Path,
    version: &str,
    context: Context,
) -> Result<pdbiox::CifProvider, Exit> {
    match pdbiox::read_component_dictionary(path, pdbiox::DictionaryVersion::new(version)) {
        Ok((provider, findings)) => {
            context.findings(&findings, &path.display().to_string());
            Ok(provider)
        }
        Err(findings) => {
            context.findings(&findings, &path.display().to_string());
            Err(Exit::of(&findings))
        }
    }
}

pub(crate) fn annotate(
    structure: &pdbiox::Structure,
    path: &Path,
    version: &str,
    context: Context,
) -> Result<pdbiox::Structure, Exit> {
    let provider = load_ccd(path, version, context)?;
    match pdbiox::apply_component_chemistry(
        structure,
        &provider,
        pdbiox::PolymerLinkPolicy::Disabled,
    ) {
        Ok(report) => {
            context.findings(&report.findings, &path.display().to_string());
            Ok(report.structure)
        }
        Err(finding) => {
            context.findings(&[finding], &path.display().to_string());
            Err(Exit::Consistency)
        }
    }
}
