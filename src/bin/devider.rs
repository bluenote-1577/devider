extern crate time;
use clap::Parser;
use devider::file_reader;
use devider::dbg;
use devider::consensus;
use devider::parse_cmd_line;
use devider::utils_frags;
use std::fs;
use std::path::Path;
use std::time::Instant;
use devider::types_structs::*;

//This makes statically compiled musl library
//much much faster. Set to default for x86 systems...
#[cfg(target_env = "musl")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[allow(deprecated)]
fn main() {
    #![allow(warnings)]
    let options = parse_cmd_line::Options::parse();
    //set threads
    rayon::ThreadPoolBuilder::new()
        .num_threads(options.num_threads)
        .build_global()
        .unwrap();
    if options.trace {
        simple_logger::SimpleLogger::new()
            .with_level(log::LevelFilter::Trace)
            .init()
            .unwrap();
    } else if options.debug{
        simple_logger::SimpleLogger::new()
            .with_level(log::LevelFilter::Debug)
            .init()
            .unwrap();
    }
    else{
        simple_logger::SimpleLogger::new()
            .with_level(log::LevelFilter::Info)
            .init()
            .unwrap();
    }

    let contig_out_dir = format!("{}", options.output_dir);

    if Path::new(&format!("{}/intermediate", contig_out_dir)).exists() && options.overwrite {

        let dir_path = &contig_out_dir;
        // Make an exception for this folder
        let pipeline_files = "pipeline_files";
        //join the path with the folder to keep
        
        let entries = fs::read_dir(dir_path).expect("Could not read output directory");
        for entry in entries {
            let entry = entry.unwrap();
            let file_name = entry.file_name();

            // If the current entry is not the folder you want to keep
            if file_name != pipeline_files {
                let path = entry.path();
                if path.is_file() {
                    fs::remove_file(path).expect("Could not remove file for devider run"); // Remove file
                } else if path.is_dir() {
                    fs::remove_dir_all(path).expect("Could not remove file for devider run") ; // Recursively remove directory
                }
            }
        }
    }
    else if Path::new(&format!("{}/intermediate",contig_out_dir)).exists() && !options.overwrite {
        log::error!(
            "Output directory {} already exists. Use --overwrite to overwrite.",
            &contig_out_dir
        );
        std::process::exit(1);
    }



    let start_t_initial = Instant::now();
    log::info!("Preprocessing VCF/Reference");
    let start_t = Instant::now();
    let mut all_contigs = file_reader::get_contigs_to_phase(&options.bam_file);
    let mut main_bam = file_reader::get_bam_readers(&options);
    log::debug!("Read BAM file successfully.");

    let mut chrom_seqs = None;

    let vcf_profile = file_reader::get_vcf_profile(&options.vcf_file, &all_contigs);
    log::debug!("Read VCF successfully.");
    if options.reference_fasta != "" {
        chrom_seqs = Some(file_reader::get_fasta_seqs(&options.reference_fasta));
        log::debug!("Read reference fasta successfully.");
    }
    log::debug!("Finished preprocessing in {:?}", Instant::now() - start_t);

    // Parse bed file and sequence ranges
    let mut bed_sequences = file_reader::get_bed_sequences(&options.bed_file);
    if let Some(seqs_to_phase) = &options.sequences_to_phase{
        for seq in seqs_to_phase{
            let seqs = seq.split(":").collect::<Vec<&str>>();
            if seqs.len() == 1{
                let seq_name = seqs[0];
                bed_sequences.push((seq_name.to_string(), None));
            }
            else if seqs.len() == 2{
                let range = seqs[1].split("-").collect::<Vec<&str>>();
                if range.len() != 2{
                    log::warn!("Invalid range for sequence {}. Ignoring ranges.", seq);
                    bed_sequences.push((seq.to_string(), None));
                }
                else{
                    let start = range[0].parse::<usize>();
                    let end = range[1].parse::<usize>();
                    if start.is_err() || end.is_err(){
                        let seq_name = seq;
                        log::warn!("Something went wrong when parsing the string {}. Should be of the form STRING:START-END. Using no range info.", seq);
                        bed_sequences.push((seq_name.to_string(), None));
                    }
                    let seq_name = seqs[0];
                    bed_sequences.push((seq_name.to_string(), Some((start.unwrap(), end.unwrap()))) );
                }
            }
            else{
                log::warn!("Something went wrong when parsing the string {}. Should be of the form STRING:START-END. Using no range info.", seq);
                bed_sequences.push((seq.to_string(), None));
            }
        }
    }

    let mut contigs_to_phase = vec![];
    if bed_sequences.len() != 0{
        contigs_to_phase = bed_sequences.iter().map(|x| (x.0.clone(), x.1)).collect();
    }
    else{
        contigs_to_phase = all_contigs.iter().map(|x| (x.clone(), None)).collect();
    }


    fs::create_dir_all(&format!("{}/intermediate", contig_out_dir)).unwrap();
    let mut warn_first_length = true;
    for (contig, range) in contigs_to_phase.iter() {
        if !vcf_profile.vcf_pos_allele_map.contains_key(contig.as_str())
            || vcf_profile.vcf_pos_allele_map[contig.as_str()].len() < options.snp_count_filter
        {
            if warn_first_length {
                log::warn!(
                    "A contig ({}) is not present, has invalid range, or has < {} variants. This warning will not be shown from now on.",
                    contig,
                    options.snp_count_filter,
                );
            }
            warn_first_length = false;
            continue;
        }

        let range_contig_str = if range != &None{
            format!("{}:{}-{}", contig, range.unwrap().0, range.unwrap().1)
        }
        else{
            format!("{}", contig)
        };

        let start_t = Instant::now();
        //log::info!("-----{}-----", contig);
        if range != &None{
            log::info!(
                "Phasing {}",
                range_contig_str
            );
        }

        let mut all_frags;
        let frags_without_snps;
        (all_frags, frags_without_snps) = file_reader::get_frags_from_bamvcf_rewrite(
            &mut main_bam,
            &vcf_profile,
            &options,
            &mut chrom_seqs,
            &contig,
            *range
        );

        log::debug!("Number of reads passing filtering: {}", all_frags.len());
        if all_frags.len() == 0 {
            log::debug!("Contig {} has no fragments", range_contig_str);
            continue;
        }

        if vcf_profile.vcf_snp_pos_to_gn_pos_map.contains_key(contig.as_str()) {
            
            let snp_to_genome_pos = &vcf_profile.vcf_snp_pos_to_gn_pos_map[contig.as_str()];

            all_frags.sort();
            for (i, frag) in all_frags.iter_mut().enumerate() {
                frag.counter_id = i;
            }

            //Get last SNP on the genome covered over all fragments.
            let length_gn = utils_frags::get_length_gn(&all_frags);
            log::info!("Contig {} has {} SNPs", range_contig_str, length_gn);

            let mut final_frags;
            final_frags = all_frags;
            let dbg_frags : Vec<FragDBG> = final_frags.iter().map(|x| dbg::frag_to_dbgfrag(x, &options)).collect();

            let dbg_frags = dbg_frags[0..options.max_frags.min(dbg_frags.len())].to_vec();
            log::debug!(
                "Reading inputs, realigning time taken {:?}",
                Instant::now() - start_t
            );

            let final_partitions = dbg::devider_run(dbg_frags, &options, &snp_to_genome_pos, &contig, *range, &vcf_profile);

            if let Some(final_partitions) = final_partitions {
                consensus::simple_consensus(
                    &mut main_bam,
                    &mut chrom_seqs, 
                    (&contig, *range),
                    &final_partitions,
                    &options,
                    &vcf_profile,
                );
            }
        }
    }
    log::info!("Total time taken is {:?}", Instant::now() - start_t_initial);
}
