//! Index slots, spaces, charts, labels and contraction (CMP-E02, DD-04, DD-06).

/// Marks a tensor: it obeys the tensor transformation law.
#[derive(Clone, Copy, Debug)]
pub struct Tensorial;
/// Marks a chart-bound quantity (not a tensor), such as Christoffel symbols.
#[derive(Clone, Copy, Debug)]
pub struct ChartBound;
