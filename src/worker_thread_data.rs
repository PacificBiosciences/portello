use std::sync::{Arc, Mutex};

use camino::Utf8Path;
use rust_htslib::bam;

/// For worker threads making indexed bam reads, this provides a persistent worker specific reader
/// for each bam file
pub struct BamReaderWorkerThreadData {
    pub bam_reader: bam::IndexedReader,
}

impl BamReaderWorkerThreadData {
    pub fn new(bam_filename: &Utf8Path) -> Self {
        let bam_reader = bam::IndexedReader::from_path(bam_filename).unwrap();
        Self { bam_reader }
    }
}

pub type BamReaderWorkerThreadDataSet = Arc<Vec<Mutex<BamReaderWorkerThreadData>>>;

pub fn get_bam_reader_worker_thread_data(
    thread_count: usize,
    bam_filename: &Utf8Path,
) -> BamReaderWorkerThreadDataSet {
    let mut worker_thread_data = Vec::new();
    for _ in 0..thread_count {
        worker_thread_data.push(Mutex::new(BamReaderWorkerThreadData::new(bam_filename)));
    }
    Arc::new(worker_thread_data)
}
