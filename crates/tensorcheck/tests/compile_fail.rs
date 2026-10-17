//! TC-F: mistakes must not compile, and the headline must be ours (DD-16).
//!
//! The cases live in the crate `ui/` (one binary per case). This test checks them all
//! with one `cargo check --message-format=json` and compares only the headline of each
//! error (`error[E0277]: …`), so changes elsewhere in the compiler output (source
//! excerpts of this library, notes) do not require updating the expectations.
//!
//! Each case states its expectation in comment lines:
//! - `// expect: error[E0277]: …` once per error, in the order the compiler reports them;
//! - `// expect: any error` when the headline is the compiler's own (ERR-18);
//! - `// expect: compiles` for code that must compile.

use std::collections::BTreeMap;
use std::path::Path;
use std::process::Command;

#[derive(Debug, PartialEq)]
enum Expect {
    Headlines(Vec<String>),
    AnyError,
    Compiles,
}

fn expectation(src: &str) -> Expect {
    let lines: Vec<&str> = src
        .lines()
        .filter_map(|l| l.strip_prefix("// expect: "))
        .collect();
    match lines.as_slice() {
        ["any error"] => Expect::AnyError,
        ["compiles"] => Expect::Compiles,
        [] => panic!("no `// expect:` line"),
        _ => Expect::Headlines(lines.iter().map(|l| l.to_string()).collect()),
    }
}

/// The JSON string value that follows the first `"key":"` in `line`, unescaped.
/// Enough for Cargo's messages, which are one JSON object per line.
fn json_str(line: &str, key: &str) -> Option<String> {
    let start = line.find(&format!("\"{key}\":\""))? + key.len() + 4;
    let mut out = String::new();
    let mut chars = line[start..].chars();
    while let Some(c) = chars.next() {
        match c {
            '"' => return Some(out),
            '\\' => match chars.next()? {
                'n' => out.push('\n'),
                't' => out.push('\t'),
                'r' => out.push('\r'),
                'u' => {
                    let hex: String = chars.by_ref().take(4).collect();
                    out.push(char::from_u32(u32::from_str_radix(&hex, 16).ok()?)?);
                }
                other => out.push(other), // `"`, `\`, `/`
            },
            c => out.push(c),
        }
    }
    None
}

#[test]
fn compile_fail() {
    let ui = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../ui");
    let target = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target/ui");
    let out = Command::new(env!("CARGO"))
        .args([
            "check",
            "--bins",
            "--keep-going",
            "--message-format=json",
            "--quiet",
        ])
        .current_dir(&ui)
        .env("CARGO_TARGET_DIR", &target)
        .output()
        .expect("run cargo check on ui/");
    let stdout = String::from_utf8(out.stdout).expect("cargo prints UTF-8");

    // Headlines and successful builds per case.
    let mut headlines: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut built: Vec<String> = Vec::new();
    for line in stdout.lines() {
        let Some(name) = json_str(line, "name") else {
            continue;
        };
        if line.contains("\"reason\":\"compiler-message\"") {
            let rendered = json_str(line, "rendered").unwrap_or_default();
            let head = rendered.lines().next().unwrap_or_default();
            if head.starts_with("error") && !head.starts_with("error: aborting") {
                headlines.entry(name).or_default().push(head.to_string());
            }
        } else if line.contains("\"reason\":\"compiler-artifact\"") {
            built.push(name);
        }
    }

    let mut cases = 0;
    let mut failures = Vec::new();
    for entry in std::fs::read_dir(ui.join("src/bin")).expect("ui/src/bin") {
        let path = entry.expect("dir entry").path();
        let name = path
            .file_stem()
            .expect("file name")
            .to_string_lossy()
            .into_owned();
        let src = std::fs::read_to_string(&path).expect("read case");
        let got = headlines.remove(&name).unwrap_or_default();
        cases += 1;
        let ok = match expectation(&src) {
            Expect::Headlines(want) => got == want,
            Expect::AnyError => !got.is_empty(),
            Expect::Compiles => got.is_empty() && built.contains(&name),
        };
        if !ok {
            failures.push(format!(
                "{name}:\n  expected {:?}\n  got      {got:?}",
                expectation(&src)
            ));
        }
    }
    assert!(cases > 0, "no cases found in ui/src/bin");
    assert!(
        headlines.is_empty(),
        "errors from unknown targets: {headlines:?}"
    );
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn json_strings_are_unescaped() {
    let line =
        r#"{"target":{"name":"err01"},"message":{"rendered":"error[E0277]: a \"b\" `c`\n -->é"}}"#;
    assert_eq!(json_str(line, "name").as_deref(), Some("err01"));
    assert_eq!(
        json_str(line, "rendered").as_deref(),
        Some("error[E0277]: a \"b\" `c`\n -->é")
    );
    assert_eq!(json_str(line, "missing"), None);
}
