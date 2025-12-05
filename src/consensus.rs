use rust_htslib::bam::ext::BamRecordExtensions;
use fxhash::FxHashMap;
use rust_htslib::{bam, bam::Read as DUMMY_NAME1};
use crate::types_structs::*;
use crate::parse_cmd_line::*;
use bio::io::fasta::IndexedReader as FastaIndexedReader;
use std::io::{BufWriter};

pub fn simple_consensus(
    main_bam: &mut bam::IndexedReader,
    _chrom_seqs: &mut Option<FastaIndexedReader<std::fs::File>>,
    contig_range: (&str, Option<(usize,usize)>),
    partition: &Vec<HapFinalResultString>,
    options: &Options,
    _vcf_profile: &VcfProfile,
    first: bool,
){
    if partition.len() == 0{
        return;
    }

    let mut record_partition = vec![Vec::new(); partition.len()];

    // Get all positions... make N if cov < min-deptha and alsounder certain conditions
    let start_s;
    let end_s;
    if let Some((start, end)) = contig_range.1{
        main_bam.fetch((contig_range.0, start as i32, end as i32)).unwrap();
        start_s = format!("{}", start);
        end_s = format!("{}", end);
    }
    else{
        main_bam.fetch(contig_range.0).unwrap();
        start_s = String::from("ALL");
        end_s = String::from("ALL");
    }
    let mut inv_index = FxHashMap::default();
    for i in 0..partition.len(){
        for frag_name in partition[i].assigned_frags.iter(){
            inv_index.insert(frag_name, i);
        }
    }
    for record in main_bam.records(){
        //let id = record.unwrap().qname().to_string();
        let id = String::from_utf8_lossy(record.as_ref().unwrap().qname()).into_owned();
        if let Some(&i) = inv_index.get(&id){
            record_partition[i].push(record.unwrap());
        }
    }

    let mut consensus_strings = Vec::new();
    for part in record_partition.iter_mut(){
        let mut map_to_allele_count = FxHashMap::default();
        let mut min_pos = std::i64::MAX;
        let mut max_pos = std::i64::MIN;
        log::trace!("Processing partition with {} reads", part.len());
        if part.len() == 0{
            continue;
        }
        for record in part.iter(){
            if record.is_secondary(){
                continue;
            }
            for aligned_pairs in record.aligned_pairs(){
                let ref_pos = aligned_pairs[1];
                let pos = aligned_pairs[0];
                let base = BYTE_TO_SEQ[record.seq()[pos as usize] as usize];
                let counts = map_to_allele_count.entry(ref_pos).or_insert([0; 4]);
                counts[base as usize] += 1;
                if ref_pos < min_pos{
                    min_pos = ref_pos;
                }
                if ref_pos > max_pos{
                    max_pos = ref_pos;
                }
            }
        }

        log::trace!("Min pos: {}, Max pos: {}", min_pos, max_pos);
        let mut consensus_bases = vec![0; (max_pos - min_pos + 1) as usize];
        let ambiguity_threshold = options.n_fraction;
        for (pos, counts) in map_to_allele_count.iter(){
            let mut total_counts = 0;
            let pos = *pos - min_pos;
            let mut max_count = 0;
            let mut max_base = 0;
            for (base, &count) in counts.iter().enumerate(){
                total_counts += count;
                if count > max_count{
                    max_count = count;
                    max_base = base;
                }
            }
            if (max_count as f64 / total_counts as f64) < ambiguity_threshold{
                max_base = 4;
            }
            if (total_counts as f64) < options.min_cov{
                max_base = 4;
            }

            consensus_bases[pos as usize] = max_base;
        }
        log::trace!("finished consensus for partition");
        let mut new_string = Vec::new();
        for base in consensus_bases.iter(){
            if *base == 4{
                new_string.push(b'N');
            }else{
                new_string.push(SEQ_TO_ASCII[*base as usize]);
            }
        }
        consensus_strings.push(new_string);
    }
    //write consensus strings to file
    //consensus file goes to options.output_dir/consensus.fasta
    let consensus_file = format!("{}/majority_vote_haplotypes.fasta", options.output_dir);
    let bufwriter;
    
    if first{
        bufwriter = BufWriter::new(std::fs::File::create(consensus_file).unwrap());
    }
    else{
        bufwriter = BufWriter::new(std::fs::OpenOptions::new().append(true).create(true).open(consensus_file).unwrap());
    }
    let mut consensus_writer = bio::io::fasta::Writer::from_bufwriter(bufwriter);
    for (i, consensus_string) in consensus_strings.iter().enumerate(){
        let id = format!("Contig:{},Range:{}-{},Haplotype:{},Abundance:{},Depth:{} SimpleConsensus", contig_range.0, start_s, end_s, i, partition[i].relative_abundances, partition[i].depth);
        let seq = String::from_utf8(consensus_string.clone()).unwrap();
        consensus_writer.write(&id, None, seq.as_bytes()).unwrap();
    }
}
