//! Component dictionary resolution for the commands that need chemistry.
//!
//! A user dictionary, when given, answers first; the dictionary the store's
//! profile names (`components:`) answers next; the structure's own definitions
//! complete the stack. Nothing is built in. An empty path stands for "no user
//! dictionary"; the layering itself is Rust's
//! [`molframe::chemistry::resolution_stack`].

use crate::exit::{Exit, or_exit};
use crate::report::Context;
use molframe::chemistry::{ComponentProvider, LayeredProvider, TableStore};
use std::path::Path;
use std::sync::Arc;

/// The path that stands for "no user dictionary".
const NO_USER_DICTIONARY: &str = "";

/// The version recorded when no dictionary is named: the profile decides.
const PROFILE_VERSION: &str = "profile";

/// The version recorded when a dictionary path is given without one.
const UNVERSIONED: &str = "unversioned";

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
    let leak_path = |path: &Path| Box::leak(Box::new(path.to_path_buf())) as &'static Path;
    let leak_version = |text: &str| Box::leak(text.to_owned().into_boxed_str()) as &'static str;
    let (path, version) = match (path, version) {
        (Some(path), Some(version)) => (leak_path(path), leak_version(version)),
        (Some(path), None) => (leak_path(path), UNVERSIONED),
        (None, Some(_)) => {
            eprintln!("--ccd-version needs --ccd");
            return Err(Exit::Usage);
        }
        (None, None) => match (context.ccd, context.ccd_version) {
            (Some(path), Some(version)) => (path, version),
            (Some(path), None) => (path, UNVERSIONED),
            _ => (Path::new(NO_USER_DICTIONARY), PROFILE_VERSION),
        },
    };
    if path.as_os_str().is_empty() {
        context.ccd = None;
    } else {
        context.ccd = Some(path);
    }
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
    let provider = or_exit!(load_ccd(path, version, context));
    let component = match provider.get(id) {
        Ok(Some(component)) => component,
        Ok(None) => {
            eprintln!("component {id:?} is absent from {}", describe(path));
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

const fn component_kind(kind: molframe::chemistry::ComponentKind) -> &'static str {
    match kind {
        molframe::chemistry::ComponentKind::AminoAcid => "amino-acid",
        molframe::chemistry::ComponentKind::Nucleotide => "nucleotide",
        molframe::chemistry::ComponentKind::Saccharide => "saccharide",
        molframe::chemistry::ComponentKind::Lipid => "lipid",
        molframe::chemistry::ComponentKind::NonPolymer => "non-polymer",
        molframe::chemistry::ComponentKind::Solvent => "solvent",
        molframe::chemistry::ComponentKind::Ion => "ion",
        molframe::chemistry::ComponentKind::Unknown => "unknown",
    }
}

fn describe(path: &Path) -> String {
    if path.as_os_str().is_empty() {
        "the profile's dictionary".to_owned()
    } else {
        path.display().to_string()
    }
}

pub(crate) fn load_ccd(
    path: &Path,
    version: &str,
    context: Context,
) -> Result<LayeredProvider, Exit> {
    let source = describe(path);
    let mut layers: Vec<Arc<dyn ComponentProvider>> = Vec::new();
    if !path.as_os_str().is_empty() {
        match molframe::read_component_dictionary(path, molframe::DictionaryVersion::new(version)) {
            Ok((provider, findings)) => {
                context.findings(&findings, &source);
                layers.push(Arc::new(provider));
            }
            Err(findings) => {
                context.findings(&findings, &source);
                return Err(Exit::of(&findings));
            }
        }
    }
    match molframe::chemistry::profile_dictionary(&TableStore::global()) {
        Ok(Some(cached)) => layers.push(cached),
        Ok(None) => {}
        Err(error) => {
            let finding = molframe::Diagnostic::from(error);
            context.findings(std::slice::from_ref(&finding), &source);
            return Err(Exit::of(std::slice::from_ref(&finding)));
        }
    }
    Ok(LayeredProvider::new(layers))
}

pub(crate) fn annotate(
    structure: &molframe::Structure,
    path: &Path,
    version: &str,
    context: Context,
) -> Result<molframe::Structure, Exit> {
    let source = describe(path);
    // Annotation touches only the components the structure uses, so a user
    // dictionary is indexed and lowered on demand instead of read whole.
    let user: Option<Arc<dyn ComponentProvider>> = if path.as_os_str().is_empty() {
        None
    } else {
        match molframe::open_component_dictionary(path, molframe::DictionaryVersion::new(version)) {
            Ok(provider) => Some(Arc::new(provider)),
            Err(findings) => {
                context.findings(&findings, &source);
                return Err(Exit::of(&findings));
            }
        }
    };
    let provider = match molframe::chemistry::resolution_stack(
        structure.engine(),
        user,
        &TableStore::global(),
    ) {
        Ok(provider) => provider,
        Err(error) => {
            let finding = molframe::Diagnostic::from(error);
            context.findings(std::slice::from_ref(&finding), &source);
            return Err(Exit::of(std::slice::from_ref(&finding)));
        }
    };
    match molframe::chemistry::apply_component_chemistry(
        structure.engine(),
        &provider,
        molframe::chemistry::PolymerLinkPolicy::Disabled,
    ) {
        Ok(report) => {
            context.findings(&report.findings, &source);
            Ok(report.structure.into())
        }
        Err(finding) => {
            context.findings(&[finding], &source);
            Err(Exit::Consistency)
        }
    }
}
