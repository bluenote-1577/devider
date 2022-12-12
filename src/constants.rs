use crate::types_structs::{GenotypeCount};
use ordered_float::OrderedFloat;
pub const NUM_ITER_OPTIMIZE:usize = 20;
pub const MIN_SHARED_READS_UNAMBIG: f64 = 3.;
pub const DIV_FACTOR: f64 = 0.05;
pub const PROB_CUTOFF: f64 = 0.0001;
pub const HAPQ_CUTOFF: u8 = 1;
pub const MERGE_CUTOFF: f64 = 0.95;
pub const SMALL_HAPLOGROUP_CUTOFF: usize = 20;
pub const SAME_SNP_DENSITY_CUTOFF: f64 = 1. / 10000.;
pub const DIST_COV_CUTOFF: GenotypeCount = OrderedFloat(5.);
pub const USE_QUAL_SCORES: bool = true;
pub const MERGE_SIMILAR_HAPLOGROUPS: bool = true;
pub const SEPARATE_BROKEN_HAPLOGROUPS: bool = true;
pub const WEIRD_SPLIT: bool = false;
