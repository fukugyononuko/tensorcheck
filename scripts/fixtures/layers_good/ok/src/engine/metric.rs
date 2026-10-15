//! Valid forms: own module, lower layers, contracts, root re-exports of contracts.
//! A comment may mention crate::manager::tensor::Tensor without breaking a rule.
use crate::access::component::Store;
use crate::contracts::Point;
use crate::engine::metric::Metric as Same;
use crate::utility::{extent::Extent, scalar::Scalar};
use crate::Quantity;

const TEXT: &str = "crate::manager::tensor::Tensor";
const RAW: &str = r#"crate::engine::index::Slot"#;

pub fn first<'a>(x: &'a str) -> char {
    let _ = (TEXT, RAW);
    'x'
}

#[cfg(test)]
mod tests {
    use super::*;
    use self::super::first;
}
