use clap::Parser;

#[derive(clap::ValueEnum, Clone, Debug, Default)]
pub enum Preset{
    #[default]
    OldLongReads,
    NanoporeR9,
    NanoporeR10,
    HiFi,
}

#[derive(Parser, Debug)]
#[command(name = "devider", version, about = "Long-read haplotyping for diverse small sequences (e.g. viruses, genes).", long_about = None)]
pub struct Options{

    /// Number of threads to use.
    #[arg(short='t', long="threads", default_value_t = 10)]
    pub num_threads: usize,

    /// Presets for different technologies. More accurate technologies use more aggressive
    /// parameters.
    #[arg(short, long, value_enum, default_value_t = Preset::NanoporeR9, help_heading = "PRESETS")]
    pub preset: Preset,

    //heading = input?
    /// Indexed bam file to phase.
    #[arg(short, long, help_heading = "INPUT")]
    pub bam_file: String,

    /// VCF file with SNPs.
    #[arg(short, long, help_heading = "INPUT")]
    pub vcf_file: String,

    /// Reference fasta file.
    #[arg(short, long, help_heading = "INPUT")]
    pub reference_fasta: String,

    /// Output directory.
    #[arg(short, long, default_value = "devider_output", help_heading = "OUTPUT")]
    pub output_dir: String,

    /// Sequences to phase separated by commas (e.g. NC_001802.1:1-1000,NC_045512.2). See also
    /// --bed-file.
    #[arg(short='S', long,value_delimiter = ',', help_heading = "INPUT")]
    pub sequences_to_phase: Option<Vec<String>>,

    /// BED file with >= 3 columns (contig, start, end). Only these regions will be phased; reads
    /// must cover 50% of the region to be considered.
    #[arg(long, help_heading = "INPUT")]
    pub bed_file: Option<String>,

    /// Do not use supplementary alignments.
    #[arg(long, help_heading = "OPTIONS")]
    pub dont_use_supp_aln: bool,

    /// Overwrite output directory if it exists.
    #[arg(short='O',long, help_heading = "OUTPUT")]
    pub overwrite: bool,

    /// Output haplotype-tagged reads.
    #[arg(long, help_heading = "OUTPUT")]
    pub output_reads: bool,

    /// Output nucleotide alleles instead of 0-1 (ref,alt) alleles.
    #[arg(long, help_heading = "OUTPUT")]
    pub allele_output: bool,

    /// > this fraction supporting the major allele to call a non-N base. 
    #[arg(long, help_heading = "OUTPUT", default_value_t = 0.66)]
    pub n_fraction: f64,


    /// Value of "k". Set automatically if not provided.
    #[arg(short, help_heading = "ALGORITHM")]
    pub k: Option<usize>,

    /// Don't use primary mappings with < --mapq-cutoff.
    #[arg(long, default_value_t = 5, help_heading = "OPTIONS")]
    pub mapq_cutoff: u8,

    /// Don't use supp. mappings with < --supp-mapq-cutoff.
    #[arg(long, default_value_t = 30, help_heading = "OPTIONS")]
    pub supp_mapq_cutoff: u8,

    /// Require supplementary alignments to be within this distance of the primary alignment.
    #[arg(long, default_value_t = 5000, help_heading = "OPTIONS")]
    pub supp_aln_dist_cutoff: i64,

    
    /// Only phase contigs with > this # of SNPs.
    #[arg(long, default_value = "1", help_heading = "OPTIONS")]
    pub snp_count_filter: usize,

    /// Trace logging (VERY VERBOSE).
    #[arg(long)]
    pub trace: bool,

    /// Debug logging.
    #[arg(long)]
    pub debug: bool,

    /// Minimum abundance (in %) of a haplotype to be considered.
    #[arg(long, default_value_t = 0.25, help_heading = "ALGORITHM")]
    pub min_abund: f64,

    /// Minimum coverage (depth) of a haplotype to be considered.
    #[arg(long, default_value_t = 5., help_heading = "ALGORITHM")]
    pub min_cov: f64,

    /// Minimum base quality to consider for fastq.
    #[arg(long, default_value_t = 3, help_heading = "OPTIONS")]
    pub min_qual: u8,

    /// Maximum number of alignments per contig.
    #[arg(long, default_value_t = 1000000000000, help_heading = "OPTIONS")]
    pub max_frags: usize,

    /// No base realignment against SNPs.
    #[arg(long, help_heading = "OPTIONS")]
    pub no_realign: bool,

    /// Merge haplotypes with < this fractional difference. Value depends on preset. 
    #[arg(long, help_heading = "ALGORITHM")]
    pub resolution: Option<f64>,

    /// FDR for strand bias filtering.
    #[arg(long, default_value_t = 0.005, help_heading = "ALGORITHM")]
    pub strand_bias_fdr: f64
}


//pub fn parse_cmd_line(matches : ArgMatches) -> Options{
//    // Set up our logger if the user passed the debug flag
//    if matches.is_present("trace") {
//        simple_logger::SimpleLogger::new()
//            .with_level(log::LevelFilter::Trace)
//            .init()
//            .unwrap();
//    } else if matches.is_present("debug"){
//        simple_logger::SimpleLogger::new()
//            .with_level(log::LevelFilter::Debug)
//            .init()
//            .unwrap();
//    }
//    else{
//        simple_logger::::new()
//            .with_level(log::LevelFilter::Info)
//            .init()
//            .unwrap();
//    }
//
//
//    let overwrite = matches.is_present("overwrite");
//    //Parse command line args.
//    let max_number_solns_str = matches.value_of("max_number_solns").unwrap_or("10");
//    let supp_aln_dist_cutoff = matches.value_of("supp_aln_dist_cutoff").unwrap_or("40000").parse::<i64>().unwrap();
//    let max_number_solns = max_number_solns_str.parse::<usize>().unwrap();
//    let num_t_str = matches.value_of("threads").unwrap_or("10");
//    let num_threads = match num_t_str.parse::<usize>() {
//        Ok(num_threads) => num_threads,
//        Err(_) => panic!("Number of threads must be positive integer"),
//    };
//
//    let max_ploidy = matches.value_of("max_ploidy").unwrap_or("5").parse::<usize>().unwrap();
//    let hybrid = matches.is_present("hybrid");
//    let reassign_short = matches.is_present("reassign_short");
//    let do_binning = matches.is_present("do_binning");
//    let trim_reads = matches.is_present("trim");
//    let short_bam_file = matches.value_of("hybrid").unwrap_or("").to_string();
//    let list_to_phase: Vec<String>;
//    if let Some(values) = matches.values_of("list_to_phase") {
//        list_to_phase = values.map(|x| x.to_string()).collect();
//    } else {
//        list_to_phase = vec![];
//    }
//
//    let bam;
//    let bam_file = match matches.value_of("bam") {
//        None => {
//            bam = false;
//            "_"
//        }
//        Some(bam_file) => {
//            bam = true;
//            bam_file
//        }
//    };
//    let bam_file = bam_file.to_string();
//
//
//    let epsilon;
//    let block_length;
//    if !matches.is_present("epsilon") || !matches.is_present("bam_block_length"){
//        //let (est_block_len, est_epsilon) = file_reader::l_epsilon_auto_detect(&bam_file);
//        let (est_block_len, est_epsilon) = (1000, 0.04);
//        if !matches.is_present("epsilon"){
//            epsilon = est_epsilon;
//        }
//        else{
//            epsilon = matches.value_of("epsilon").unwrap().parse::<f64>().unwrap();
//        }
//        if !matches.is_present("bam_block_length"){
//            block_length = est_block_len;
//        }
//        else{
//            block_length = matches.value_of("bam_block_length").unwrap().parse::<usize>().unwrap();
//        }
//    }
//    else {
//        epsilon = matches.value_of("epsilon").unwrap().parse::<f64>().unwrap();
//        block_length = matches.value_of("bam_block_length").unwrap().parse::<usize>().unwrap();
//    }
//
////    let mut epsilon = 0.04;
////    if hybrid{
////        epsilon = 0.03
////    }
////    if matches.is_present("epsilon"){
////        epsilon = matches.value_of("epsilon").unwrap().parse::<f64>().unwrap();
////    }
////    let block_length = matches.value_of("bam_block_length").unwrap_or("15000");
////    let block_length = block_length.parse::<usize>().unwrap();
//    //    let use_mec = matches.is_present("use_mec");
//    let reference_fasta = matches.value_of("reference_fasta").unwrap_or("").to_string();
//    let dont_use_supp_aln = matches.is_present("dont use supplementary");
//    let gzip = matches.is_present("gzip-reads");
//    //If the user is splitting the bam file according to the output partition.
//    let out_dir = matches
//        .value_of("output dir")
//        .unwrap_or("floria_out_dir")
//        .to_string();
//    let snp_density = matches
//        .value_of("snp_density")
//        .unwrap_or("0.0005")
//        .parse::<f64>()
//        .unwrap();
//
//    if Path::new(&out_dir).exists() && !overwrite{
//        error!("Output directory exists; output directory must not be an existing directory. Use --overwrite to overwrite existing directory.");
//        std::process::exit(1);
//    }
//
//    std::fs::create_dir_all(&out_dir).unwrap();
//    let mut cmd_file =
//        File::create(format!("{}/cmd.log", out_dir)).expect("Can't create file");
//    for arg in env::args() {
//        write!(cmd_file, "{} ", arg).unwrap();
//    }
//
//    //Overwrite the contig_ploidy_info file. 
//    let mut all_ploidy_file = File::create(format!("{}/contig_ploidy_info.tsv", out_dir)).expect("Can't create file");
//    write!(
//        all_ploidy_file,
//        "{}",
//        CONTIG_PLOIDY_HEADER
//    )
//    .unwrap();
//
//
//    //If the user is getting frag files from BAM and VCF.
//
//    let vcf_file = matches.value_of("vcf").unwrap().to_string();
//
//    if !bam{
//        panic!("Must input a BAM file.")
//    }
//
//    let snp_count_filter = matches.value_of("snp_count_filter").unwrap_or("100").parse::<usize>().unwrap();
//    let use_qual_scores = matches.is_present("use_qual_scores");
//    let output_reads = matches.is_present("output reads");
//    let mapq_cutoff = matches.value_of("mapq_cutoff").unwrap_or("5").parse::<u8>().unwrap();
//    
//    
//
//    rayon::ThreadPoolBuilder::new()
//        .num_threads(num_threads)
//        .build_global()
//        .unwrap();
//
//    let stopping_heuristic = !matches.is_present("no stop heuristic");
//    let ignore_monomorphic = matches.is_present("ignore monomorphic");
//    let ploidy_sensitivity = matches.value_of("ploidy sensitivity").unwrap_or("2").parse::<u8>().unwrap();
//    if !(ploidy_sensitivity >= 1 && ploidy_sensitivity <= 3){
//        log::error!("Ploidy sensitivty option must be between 1 and 3");
//        std::process::exit(1);
//    }
//
//    let opt = Options{
//        bam_file,
//        vcf_file,
//        use_qual_scores,
//        gzip,
//        output_reads,
//        mapq_cutoff,
//        epsilon,
//        dont_use_supp_aln,
//        reassign_short,
//        do_binning,
//        max_number_solns,
//        snp_density,
//        max_ploidy,
//        out_dir,
//        hybrid,
//        list_to_phase,
//        block_length,
//        reference_fasta,
//        trim_reads,
//        short_bam_file,
//        snp_count_filter,
//        stopping_heuristic,
//        ignore_monomorphic,
//        num_threads,
//        overwrite,
//        ploidy_sensitivity,
//        supp_aln_dist_cutoff
//    };
//    opt
//}
