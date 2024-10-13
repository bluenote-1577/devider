use crate::parse_cmd_line;
use fishers_exact::fishers_exact;
use rust_lapper::*;
use std::sync::Mutex;
use rayon::prelude::*;
use crate::parse_cmd_line::Options;
use crate::types_structs::*;
use crate::utils_frags;
use disjoint_sets::UnionFind;
use fxhash::{FxHashMap, FxHashSet};
use ordered_float::*;
use statrs::distribution::{Binomial, DiscreteCDF};
use std::collections::VecDeque;
use std::fmt;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Arc;

pub fn devider_run(
    mut dbg_frags:Vec<FragDBG>,
    options: &Options,
    snp_pos_to_genome_pos: &Vec<usize>,
    contig_name: &str,
    range: Option<(usize, usize)>,
    vcf_profile: &VcfProfile,
) -> Option<Vec<HapFinalResultString>> {
    let k;
    let mut thirty = utils_frags::get_avg_length_dbgf(&dbg_frags, 0.33);
    let fifty = utils_frags::get_avg_length_dbgf(&dbg_frags, 0.5);
    let max_k_preset;
    let max_median;
    let mut resolution;
    match options.preset {
        parse_cmd_line::Preset::OldLongReads => {
            max_k_preset = 10;
            max_median = 50;
            resolution = 0.02;
        }

        parse_cmd_line::Preset::NanoporeR9 => {
            max_k_preset = 20;
            max_median = 150;
            resolution = 0.01;
        }
        parse_cmd_line::Preset::NanoporeR10 => {
            max_k_preset = 35;
            max_median = 250;
            resolution = 0.005;
        }
        parse_cmd_line::Preset::HiFi => {
            max_k_preset = 100;
            max_median = 500;
            resolution = 0.001;
        }
    }

    let coverage_divider = (100. / options.min_abund) as u64;

    if let Some(res) = options.resolution {
        resolution = res;
    }

    log::debug!("33% number of SNPs in a read is {}", thirty);
    log::debug!("50th perc. number of SNPs in a read is {}", fifty);

    //subsample if too many snps
    let division_factor = fifty as f64 / max_median as f64;
    let mut snp_pos_to_genome_pos_new = snp_pos_to_genome_pos.clone();
    let mut subsampled_positions = FxHashSet::default();
    let mut old_pos_to_index_map = FxHashMap::default();
    if division_factor > 1.{
        log::trace!("Subsampling SNPs; division factor is {}", division_factor);
        for i in 0..snp_pos_to_genome_pos.len() {
            let subsampled_pos = (i as f64 * division_factor) as usize;
            if subsampled_pos >= snp_pos_to_genome_pos.len() {
                break;
            }
            subsampled_positions.insert(subsampled_pos as u32 + 1);
            old_pos_to_index_map.insert(subsampled_pos as u32 + 1, i as u32 + 1);
        }
        let mut ss_pos_vec = subsampled_positions.iter().collect::<Vec<_>>();
        ss_pos_vec.sort();
        log::trace!("Subsampled indices: {:?}", ss_pos_vec);
    }

    if !subsampled_positions.is_empty(){
        log::debug!("Subsampling SNPs");
        snp_pos_to_genome_pos_new = subsample_positions_fragdbg(&mut dbg_frags, &subsampled_positions, &old_pos_to_index_map, snp_pos_to_genome_pos);
        thirty = utils_frags::get_avg_length_dbgf(&dbg_frags, 0.33);
    }

    let num_snps = snp_pos_to_genome_pos_new.len();
    let snp_pos_to_genome_pos_new = strand_bias_filter(&mut dbg_frags, options, num_snps, &snp_pos_to_genome_pos_new);
    if snp_pos_to_genome_pos_new.len() < num_snps / 10 && num_snps > 10{
        log::warn!("{} has > 90% of SNPs filtered out by strand bias. Maybe coverage is very high. ", contig_name);
    }
    let num_snps = snp_pos_to_genome_pos_new.len();

    if let Some(opt_k) = options.k {
        k = opt_k;
    } else {
        k = thirty.min(num_snps * 3 / 4).min(max_k_preset).max(1);
    }

    if k > num_snps {
        log::warn!("Not enough SNPs; exiting.");
        return None;
    }

    //disable this for now
    log::trace!("Start k: {}", k);
    let end = 0;

    let dbg = dbg_from_frags(&dbg_frags, k, None, None, None);
    log::debug!("Constructed DBG for k = {}", k);
    let (mut kmer_count, used_snp_positions) = count_kmers(&dbg_frags, k);
    let num_snps_range = used_snp_positions.len();
    let total_cov = kmer_count.iter().fold(0, |acc, (_varmer, cov)| acc + cov);
    let min_cov = u64::max(
        total_cov / (num_snps_range as u64 - k as u64 + 1) / coverage_divider,
        2,
    );
    log::debug!("Minimum coverage for global filter is : {:?}", min_cov);

    let mut dbg = filter_dbg(dbg, Some(min_cov), None, k, false, num_snps_range);
    print_dbg(&dbg, &format!("{}/intermediate/dbg.dot", options.output_dir));

    let mut uni = get_unitigs(&dbg, k, false);
    kmer_count.retain(|varmer, _cov| dbg.contains_key(varmer));
    let step = 1;

    for l in (step..end + 1).step_by(step) {
        dbg = dbg_from_frags(&dbg_frags, l + k, Some(dbg), Some(&uni), Some(step));
        uni = get_unitigs(&dbg, k + l, false);
    }
    print_dbg(&uni, format!("{}/intermediate/unitigs.dot", options.output_dir).as_str());

    //Remove tips
//    for _ in 0..2{
//        let tips = remove_tips(&uni, options, k + end);
//        uni = filter_dbg(uni, None, Some(tips), k + end, false, num_snps_range);
//        print_dbg(&uni, format!("{}/intermediate/tips_removed.dot", options.output_dir).as_str());
//        uni = get_unitigs(&uni, k + end, true);
//        print_dbg(&uni, format!("{}/intermediate/tips_removed_unitigs.dot", options.output_dir).as_str());
//    }

    //Unitigging
    let unitigs = uni;

    //Query and clean
    log::debug!("Cleaning unitigs");
    let mut final_unitigs = unitigs;
    for i in 1..3 {
        let bad_unitigs = query_unitigs(&final_unitigs, i);
        log::debug!("Number of bad unitigs {}", bad_unitigs.len());
        let filtered_unitigs = filter_dbg(final_unitigs, None, Some(bad_unitigs), k + end, false, num_snps_range);
        print_dbg(&filtered_unitigs, format!("{}/intermediate/clean_dbg_{}.dot", options.output_dir, i).as_str());
        final_unitigs = get_unitigs(&filtered_unitigs, k + end, true);
        print_dbg(&final_unitigs, format!("{}/intermediate/clean_unitigs_{}.dot",options.output_dir, i).as_str());
    }

    //Remove tips again
    let tips = remove_tips(&final_unitigs, options, k + end);
    final_unitigs = filter_dbg(final_unitigs, None, Some(tips), k + end, false, num_snps_range);
    print_dbg(&final_unitigs, format!("{}/intermediate/tips_removed_round2.dot", options.output_dir).as_str());
    final_unitigs = get_unitigs(&final_unitigs, k + end, true);

    
    //Remove small disconnected components that have length < 1.5 * k and coverage < mean_cov / 100
    let min_cov_small_disconnected = u64::max(total_cov / (num_snps_range as u64 - k as u64 + 1) / (coverage_divider) * 4, 2);
    log::trace!("Second round min cov :{}", min_cov_small_disconnected);
    final_unitigs = filter_dbg(final_unitigs, Some(min_cov_small_disconnected), None, k + end, true, num_snps_range);
    final_unitigs = get_unitigs(&final_unitigs, k + end, true);

    let final_unitigs = clean_hanging_kmers(final_unitigs, k + end - 1);
    print_dbg(&final_unitigs, format!("{}/intermediate/cleaned_unitigs.dot", options.output_dir).as_str());

    //Try aligning reads to graph
    log::debug!("Aligning reads to graph of size {}", final_unitigs.len());
    let vec_df_all: Vec<DictFrag> = dbg_frags
        .iter()
        .map(|frag| fragdbg_to_dictfrag(frag))
        .collect();
    let vec_df_graph: Vec<DictFrag> = final_unitigs
        .iter()
        .map(|(varmer, info)| DictFrag {
            seq: varmer.iter().cloned().collect(),
            seq_vec: varmer.clone(),
            first_position: varmer.first().unwrap().0,
            last_position: varmer.last().unwrap().0,
            cov: info.coverage,
        })
        .collect();

    let dict_frag_to_index = vec_df_graph
        .iter()
        .enumerate()
        .map(|(i, frag)| (frag, i))
        .collect::<FxHashMap<&DictFrag, usize>>();

    let path_dict = Mutex::new(FxHashMap::default());
    //for dict_frag in vec_df_all.iter() {
    log::trace!("ALIGN TO GRAPH");
    vec_df_all.into_par_iter().for_each(|dict_frag| {
        print_varmer_d(&dict_frag, true);
        let hits = get_hits(&dict_frag, &vec_df_graph, 100000, 10, true);
        let dp_res = dp_hits(
            &hits,
            &dict_frag,
            &final_unitigs,
            10000,
            -3.0,
            false,
            false,
            GraphConstraint::RequireDagRescue,
            1000,
            k,
            &FxHashSet::default()
        );
        let varmers = varmers_from_dp_res(&dp_res, 0.2);
        log::trace!("HITS: {:?}", hits.len());
        log::trace!("4P RES: {:?}", dp_res.score);
        for varmer in varmers.iter() {
            print_varmer_d(varmer, true);
        }
        log::trace!("FIN DP RES");
        *path_dict.lock().unwrap().
            entry(varmers).or_insert(0) += 1;
    });

    let mut unitig_paths = vec![];
    let path_dict = path_dict.into_inner().unwrap();
    if path_dict.is_empty() {
        log::error!("No paths found. Exiting.");
        return None;
    }
    let mut counts = path_dict.iter().map(|(_, count)| *count).collect::<Vec<_>>();
    counts.sort();
    let median_count = counts[counts.len() / 2];

    for (varmers, count) in path_dict.into_iter() {
        log::trace!("INITIAL PATH COUNT: {}", count);
        for varmer in varmers.iter() {
            print_varmer_d(varmer, true);
        }
        if varmers.len() == 0 {
            continue;
        }
        if count < (min_cov - 1).min(3).min(median_count) {
            continue;
        }
        let path = VarmerPath {
            first: varmers.first().unwrap().first_position,
            last: varmers.last().unwrap().last_position,
            varmers,
            total_avg_cov: count,
        };
        unitig_paths.push(path);
    }

    log::debug!("Number of total paths: {}", unitig_paths.len());
    let integer_paths = get_outside_paths_and_integers(&unitig_paths, &dict_frag_to_index);
    log::debug!("Number of candidate outside paths: {}", integer_paths.len());
    let assembly_graph = get_assembly_integer_graph(&integer_paths);

    print_dbg(&assembly_graph, format!("{}/intermediate/assembly_graph.dot", options.output_dir).as_str());

    let integer_unitigs = get_unitigs(&assembly_graph, 1, true);
    //let integer_unitigs = assembly_graph;
    let mut paths = vec![];
    for (int_unitig, info) in integer_unitigs.iter() {
        let path_as_df = int_unitig
            .iter()
            .map(|(_, x)| &vec_df_graph[(*x) as usize])
            .collect::<Vec<_>>();
        let min_cov_path = path_as_df.iter().map(|x| x.cov).min().unwrap();
        if min_cov_path / 3 > info.coverage {
            log::trace!(
                "PATH ALIGN VS UNITIG COV CUTOFF - FAILED -- MIN COV PATH: {} INFO COV: {}",
                min_cov_path,
                info.coverage
            );
            for df in path_as_df.iter() {
                print_varmer_d(df, true);
            }
        } else {
            log::trace!(
                "PATH ALIGN VS UNITIG COV CUTOFF - INTEGER UNITIG PASSED -- MIN COV PATH: {} INFO COV: {}",
                min_cov_path,
                info.coverage
            );
            for df in path_as_df.iter() {
                print_varmer_d(df, true);
            }
            let path_as_varmers = int_unitig
                .iter()
                .map(|(_, x)| (&vec_df_graph[(*x) as usize].seq_vec, info.coverage as usize))
                .collect::<Vec<_>>();
            paths.push(path_as_varmers);
        }
    }

    log::debug!("Number of candidate integer unitig paths passing filters: {}", paths.len());
    let mut hap_path_results = get_path_haps(&dbg_frags, &final_unitigs, paths, num_snps, options);
    print_final_hap_results(
        &hap_path_results,
        num_snps,
        options,
        "intermediate/hap_before.txt",
        "intermediate/id_before.txt",
        None, 
        (contig_name,range),
        None,
        vcf_profile,
        &snp_pos_to_genome_pos_new,
    );

    let mut j = 0;
    loop{
        j+=1;
        log::debug!("Consensus round {}", j);
        let final_results_consensus = consensus(
            &hap_path_results,
            num_snps,
            &snp_pos_to_genome_pos_new,
            &dbg_frags,
            &format!("intermediate/hap_info-{}.txt", j),
            options,
            (contig_name,range),
            false,
            0.0,
        ).0;

        //let final_results_consensus = filter_final_haplotypes(final_results_consensus, options);

        print_final_hap_results(
            &final_results_consensus,
            num_snps,
            options,
            format!("intermediate/haplotypes-{}.fasta", j).as_str(),
            format!("intermediate/ids-{}.txt", j).as_str(),
            None,
            (contig_name, range),
            None,
            vcf_profile,
            &snp_pos_to_genome_pos_new,
        );


        let mut same = false;
        if hap_path_results.len() == final_results_consensus.len() {
            same = true;
        }

        if same{
            hap_path_results = final_results_consensus;
            let mut unassigned;
            loop{
                log::debug!("Semifinal consensus");
                let (final_results, unassigned_loop) = consensus(
                    &hap_path_results,
                    num_snps,
                    &snp_pos_to_genome_pos_new,
                    &dbg_frags,
                    &format!("intermediate/hap_info-semi.txt"),
                    options,
                    (contig_name, range),
                    false,
                    resolution
                );
                unassigned = unassigned_loop;

                let final_results_filtered = filter_final_haplotypes(final_results, options);
                //let final_results_filtered = final_results;
                if final_results_filtered.len() == hap_path_results.len(){
                    break;
                }
                else{
                    hap_path_results = final_results_filtered;
                }
            }

            let final_results_filtered = hap_path_results;
            let output_reads = if options.output_reads {
                Some("reads.fq")
            } else {
                None
            };
            print_final_hap_results(
                &final_results_filtered,
                num_snps,
                options,
                "snp_haplotypes.fasta",
                "ids.txt",
                output_reads,
                (contig_name, range),
                Some(&unassigned),
                vcf_profile,
                &snp_pos_to_genome_pos_new,
            );

            log::debug!("Final consensus");
            let _only_for_printing = consensus(
                &final_results_filtered,
                num_snps,
                &snp_pos_to_genome_pos_new,
                &dbg_frags,
                &format!("hap_info.txt"),
                options,
                (contig_name, range),
                false,
                resolution,
            );
            hap_path_results = final_results_filtered.clone();
            break;
        }
        hap_path_results = final_results_consensus;
    }

    let mut final_results_strings = vec![];

    for res in hap_path_results.iter() {
        let abund = res.relative_abundances;
        let depth = res.depth;
        let mut read_ids = vec![];
        for frag in res.assigned_frags.iter() {
            read_ids.push(frag.id.clone());
        }

        let hap_res_str = HapFinalResultString {
            relative_abundances: abund,
            depth,
            assigned_frags: read_ids,
        };
        
        final_results_strings.push(hap_res_str);
    }

    return Some(final_results_strings);

}

fn filter_final_haplotypes<'a>(
    final_results: Vec<HapFinalResult<'a>>,
    options: &'a Options,
) -> Vec<HapFinalResult<'a>> {
    let mut filtered_results = vec![];
    for (i,res) in final_results.iter().enumerate() {
        if res.relative_abundances < options.min_abund {
            log::debug!("Haplotype {} has relative abundance of {} which is less than the minimum abundance of {}. Skipping", i, res.relative_abundances, options.min_abund);
            continue;
        }
        if res.depth < options.min_cov {
            log::debug!("Haplotype {} has average coverage depth of {}x which is less than the minimum depth of {}x. Skipping", i, res.depth, options.min_cov);
            continue;
        }
        filtered_results.push(res.clone());
    }
    return filtered_results;
}

fn consensus<'a>(
    hap_path_results: &Vec<HapFinalResult>,
    snps: usize,
    snp_pos_to_genome_pos: &Vec<usize>,
    dbg_frags: &'a Vec<FragDBG>,
    consensus_file_loc: &str,
    options: &Options,
    contig: (&str, Option<(usize,usize)>),
    only_print: bool,
    resolution: f64,
) -> (Vec<HapFinalResult<'a>>, Vec<&'a FragDBG>) {

    let contig_name = contig.0;
    let start;
    let end;
    if let Some(range) = contig.1 {
        start = format!("{}", range.0);
        end = format!("{}", range.1);
    } else {
        start = String::from("ALL");
        end = String::from("ALL");
    }
    let previous_total_depth = hap_path_results.iter().map(|x| x.depth).sum::<f64>();
    let dir = Path::new(&options.output_dir);
    let cons_file = dir.join(consensus_file_loc);
    let cons_file = cons_file.to_str().unwrap();

    let mut consensus_file;
    if Path::exists(Path::new(cons_file)) {
        consensus_file = BufWriter::new(std::fs::OpenOptions::new().append(true).open(cons_file).unwrap());
    }
    else{
        consensus_file = BufWriter::new(std::fs::File::create(cons_file).unwrap());
    }
    let mut haps = vec![];
    let mut vec_haps = vec![];
    for res in hap_path_results.iter() {
        let hap = utils_frags::fragdbg_to_seq_dict(&res.assigned_frags, true);

        let vec_form = hap
            .iter()
            .map(|(snp_pos, geno_dicts)| {
                let mut sorted_genos_by_cov = geno_dicts.iter().collect::<Vec<_>>();
                sorted_genos_by_cov.sort_by(|a, b| b.1.cmp(&a.1));
                let best_geno = sorted_genos_by_cov[0].0;
                let ratio =
                    sorted_genos_by_cov[0].1 / geno_dicts.values().sum::<OrderedFloat<f64>>();
                (*snp_pos, format!("{}:{:.2}", best_geno, ratio))
            })
            .collect::<FxHashMap<u32, String>>();

        haps.push(hap);
        vec_haps.push(vec_form);
    }
    consensus_file.write(&format!("Contig:{},Range:{}-{}", contig_name, start, end).as_bytes()).unwrap();
    for i in 0..vec_haps.len() {
        consensus_file
            .write(format!("\tHaplotype:{}", i).as_bytes())
            .unwrap();
    }
    consensus_file.write(b"\n").unwrap();
    for i in 1..snps + 1 {
        consensus_file
            .write(format!("{}", snp_pos_to_genome_pos[i - 1] + 1).as_bytes())
            .unwrap();
        for j in 0..vec_haps.len() {
            if vec_haps[j].contains_key(&(i as u32)) {
                let geno = vec_haps[j].get(&(i as u32)).unwrap();
                consensus_file
                    .write(format!("\t{}", geno).as_bytes())
                    .unwrap();
            } else {
                consensus_file.write(b"\t-").unwrap();
            }
        }
        consensus_file.write(b"\n").unwrap();
    }
    if only_print {
        return (vec![], vec![]);
    }

    let mut union_find = UnionFind::new(haps.len());
    let dist_cutoff = (previous_total_depth/400.).max(2.);
    for i in 0..haps.len() {
        for j in 0..i {
            log::trace!("Comparing haplotypes {} and {}", i, j);
            let (same, diff) =
                utils_frags::distance_between_haplotypes(&haps[i], &haps[j], &(0, u32::MAX), 0.75, dist_cutoff);
            log::trace!("Same: {}, Diff: {}", same, diff);
            if same < 2. {
                continue;
            }
            if diff == 0. || (diff as f64 / same as f64) < resolution{
                union_find.union(i,j);
            }
//            if (diff as f64 / same as f64) < error_ratio {
//                union_find.union(i, j);
//            }
        }
    }

    //deduplicate haps based on parameters
    union_find.force();

    //force representative to be the highest depth rep
    let mut best_depth_root = FxHashMap::default();
    for i in 0..haps.len() {
        let root = union_find.find(i);
        if let Some((_, depth)) = best_depth_root.get(&root) {
            if hap_path_results[i].depth > *depth {
                best_depth_root.insert(root, (i, hap_path_results[i].depth));
            }
        } else {
            best_depth_root.insert(root, (i, hap_path_results[i].depth));
        }
    }

    let mut new_final_results_map = FxHashMap::default();
    for i in 0..haps.len() {
        let root = best_depth_root[&union_find.find(i)].0;

        let res = new_final_results_map.entry(root).or_insert(HapFinalResult {
            relative_abundances: 0.,
            depth: 0.,
            assigned_frags: vec![],
            path_frag: DictFrag::default(),
        });

        if res.path_frag.seq.is_empty() {
            let consensus_seq_dict = utils_frags::get_consensus_seq_dict(&haps[root]);
            let mut seq_varmer = consensus_seq_dict
                .iter()
                .map(|(x, y)| (*x, *y))
                .collect::<Vec<_>>();
            seq_varmer.sort_by(|a, b| a.0.cmp(&b.0));
            if seq_varmer.is_empty() {
                continue;
            }
            let first = seq_varmer.first().unwrap().0;
            let last = seq_varmer.last().unwrap().0;
            let path_frag_cons = DictFrag {
                seq: consensus_seq_dict,
                seq_vec: seq_varmer,
                first_position: first,
                last_position: last,
                cov: 0,
            };
            res.path_frag = path_frag_cons;
        }
        //this is not correct, but not used
        res.relative_abundances += hap_path_results[i].relative_abundances;
        res.depth += hap_path_results[i].depth;
    }

    let mut final_results_consensus = new_final_results_map.into_values().collect::<Vec<_>>();
    let unassigned = reassign_frags(dbg_frags, &mut final_results_consensus, true);

    for res in final_results_consensus.iter_mut() {
        if res.assigned_frags.is_empty() {
            continue;
        }
        let seq_dict = utils_frags::fragdbg_to_seq_dict(&res.assigned_frags, true);
        let total_bases = seq_dict.values().map(|x| x.values().sum::<OrderedFloat<f64>>()).sum::<OrderedFloat<f64>>();
        let depth = total_bases/seq_dict.len() as f64;
        let consensus_seq_dict = utils_frags::get_consensus_seq_dict(&seq_dict);
        res.depth = depth.into_inner();
        let mut seq_varmer = consensus_seq_dict
            .iter()
            .map(|(x, y)| (*x, *y))
            .collect::<Vec<_>>();
        seq_varmer.sort_by(|a, b| a.0.cmp(&b.0));
        let first = seq_varmer.first().unwrap().0;
        let last = seq_varmer.last().unwrap().0;
        res.path_frag = DictFrag {
            seq: consensus_seq_dict,
            seq_vec: seq_varmer,
            first_position: first,
            last_position: last,
            cov: 0,
        };
    }

    let total_depth = final_results_consensus.iter().map(|x| x.depth).sum::<f64>();
    final_results_consensus.iter_mut().for_each(|res| { res.relative_abundances = 100. * res.depth / total_depth; });
    final_results_consensus
        .sort_by(|a, b| a.path_frag.first_position.cmp(&b.path_frag.first_position));
    return (final_results_consensus, unassigned);
}

pub fn frag_to_dbgfrag(frag: &Frag, options: &Options) -> FragDBG {
    let mut seq = vec![];
    let mut seq_dict = FxHashMap::default();
    let mut first_position = u32::MAX;
    let mut last_position = 0;
    for (pos, geno) in frag.seq_dict.iter() {
        if frag.qual_dict[pos] < options.min_qual {
            continue;
        }
        if *pos < first_position {
            first_position = *pos;
        }
        if *pos > last_position {
            last_position = *pos;
        }
        seq.push((*pos, *geno));
        seq_dict.insert(*pos, *geno);
    }
    seq.sort_by(|a, b| a.0.cmp(&b.0));
    let toret = FragDBG {
        id: frag.id.clone(),
        counter_id: frag.counter_id,
        seq,
        seq_dict: frag.seq_dict.clone(),
        seq_string: frag.seq_string.clone(),
        qual_string: frag.qual_string.clone(),
        first_position,
        last_position,
        snp_pos_to_seq_pos: frag.snp_pos_to_seq_pos.clone(),
        qual_dict: frag.qual_dict.clone(),
        forward_strand: frag.forward_strand,
    };
    toret
}

impl fmt::Display for VarmerPath {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let seq1 = self
            .varmers
            .iter()
            .map(|varmer| fmt_varmer_d(varmer))
            .collect::<Vec<String>>()
            .join("->");
        write!(f, "{}", seq1)
    }
}

fn get_edges_varmers(dbg: &mut FxHashMap<VarMer, DBGInfo>, k: usize) {
    let mut suffixes: FxHashMap<VarMer, Vec<Arc<VarMer>>> = FxHashMap::default();
    let mut prefixes: FxHashMap<VarMer, Vec<Arc<VarMer>>> = FxHashMap::default();

    for (varmer, _info) in dbg.iter() {
        //dbg!(varmer.len(),k);
        let suffix = varmer[varmer.len() - k + 1..].to_vec();
        let prefix = varmer[..k - 1].to_vec();
        let varmer_ref = Arc::new(varmer.clone());

        suffixes
            .entry(suffix)
            .or_insert_with(Vec::new)
            .push(varmer_ref.clone());
        prefixes
            .entry(prefix)
            .or_insert_with(Vec::new)
            .push(varmer_ref.clone());
    }

    for (prefix, out_varmers) in prefixes.iter() {
        if let Some(in_varmers) = suffixes.get(prefix) {
            for in_varmer in in_varmers.iter() {
                for out_varmer in out_varmers.iter() {
                    if let Some(in_info) = dbg.get_mut(&**in_varmer) {
                        in_info.out_varmers.push(Arc::clone(out_varmer));
                    }
                    if let Some(out_info) = dbg.get_mut(&**out_varmer) {
                        out_info.in_varmers.push(Arc::clone(in_varmer));
                    }
                }
            }
        }
    }
}

pub fn print_dbg(dbg: &FxHashMap<VarMer, DBGInfo>, file_name: &str) {
    let mut dot = String::from("digraph G {\n");
    for (node, info) in dbg.iter() {
        let out_edges = &info.out_varmers;
        if out_edges.len() == 0 && info.in_varmers.len() == 0 {
            let mut s1 = node
                .iter()
                .map(|(_, geno)| format!("{}", geno))
                .collect::<Vec<String>>()
                .join("");
            s1.push_str(&format!(
                "_{}-{}_COV:{}_LEN:{}",
                node.first().unwrap().0,
                node.last().unwrap().0,
                info.coverage,
                node.len()
            ));
            dot.push_str(&format!("    \"{}\" -> \"{}\";\n", s1, s1));
        }
        for out in out_edges {
            let mut s1 = node
                .iter()
                .map(|(_, geno)| format!("{}", geno))
                .collect::<Vec<String>>()
                .join("");
            let mut s2 = out
                .iter()
                .map(|(_, geno)| format!("{}", geno))
                .collect::<Vec<String>>()
                .join("");
            //append coverage to s1 and s2
            s1.push_str(&format!(
                "_{}-{}_COV:{}_LEN:{}",
                node.first().unwrap().0,
                node.last().unwrap().0,
                info.coverage,
                node.len()
            ));
            s2.push_str(&format!(
                "_{}-{}_COV:{}_LEN:{}",
                out.first().unwrap().0,
                out.last().unwrap().0,
                dbg.get(&**out).unwrap().coverage,
                out.len()
            ));
            dot.push_str(&format!("    \"{}\" -> \"{}\";\n", s1, s2));
        }
    }
    dot.push_str("}\n");

    //write to file called t.dot
    std::fs::write(file_name, dot).expect("Unable to write file");
}

fn dbg_from_frags(
    dbg_frags: &Vec<FragDBG>,
    k: usize,
    prev_dbg: Option<FxHashMap<VarMer, DBGInfo>>,
    unitigs: Option<&FxHashMap<VarMer, DBGInfo>>,
    base_k: Option<usize>,
) -> FxHashMap<VarMer, DBGInfo> {
    // DBG Format: key -> edges
    let mut dbg: FxHashMap<VarMer, DBGInfo> = FxHashMap::default();
    let mut sequences = vec![];

    if let Some(unitigs) = unitigs.as_ref() {
        for (unitig, _) in unitigs.iter() {
            sequences.push(unitig);
        }
    }
    for frag in dbg_frags.iter() {
        sequences.push(&frag.seq);
    }

    for seq in sequences.iter() {
        let mut varmer = VecDeque::new();
        for (curr_pos, curr_geno) in seq.into_iter() {
            varmer.push_back((*curr_pos, *curr_geno));
            if varmer.len() == k {
                let new_varmer = varmer.iter().cloned().collect::<VarMer>();
                if let Some(prev_dbg) = prev_dbg.as_ref() {
                    let base_k = base_k.unwrap();
                    let prev_kmer = &new_varmer[0..k - base_k];
                    let end_base_kmer = &new_varmer[base_k..];
                    let min_cov;
                    if let Some(prev_info) = prev_dbg.get(prev_kmer) {
                        if let Some(end_cov) = prev_dbg.get(end_base_kmer) {
                            min_cov = prev_info.coverage.min(end_cov.coverage);
                        } else {
                            min_cov = 0;
                        }
                    } else {
                        min_cov = 0;
                    }

                    if min_cov != 0 {
                        dbg.entry(new_varmer).or_insert(DBGInfo {
                            out_varmers: vec![],
                            in_varmers: vec![],
                            coverage: min_cov,
                        });
                    }
                } else {
                    let info = dbg.entry(new_varmer).or_insert(DBGInfo {
                        out_varmers: vec![],
                        in_varmers: vec![],
                        coverage: 0,
                    });
                    info.coverage += 1;
                }
                varmer.pop_front();
            }
        }
    }

    get_edges_varmers(&mut dbg, k);
    return dbg;
}

fn filter_dbg(
    dbg: FxHashMap<VarMer, DBGInfo>,
    min_cov: Option<u64>,
    bad_nodes: Option<Vec<VarMer>>,
    k: usize,
    only_disconnected_small: bool,
    num_snps : usize,
) -> FxHashMap<VarMer, DBGInfo> {
    let mut new_dbg = FxHashMap::default();
    for (node, info) in dbg.iter() {
        if let Some(bad_nodes) = bad_nodes.as_ref() {
            if bad_nodes.contains(node) {
                continue;
            }
        }
        if let Some(min_cov) = min_cov {
            if only_disconnected_small{
                if info.coverage < min_cov && info.in_varmers.len() == 0 && info.out_varmers.len() == 0 && node.len() < num_snps / 2{
                    continue;
                }
                
            }
            else if info.coverage <= min_cov {
                continue;
            }
        }
        new_dbg.insert(
            node.clone(),
            DBGInfo {
                out_varmers: vec![],
                in_varmers: vec![],
                coverage: info.coverage,
            },
        );
    }
    get_edges_varmers(&mut new_dbg, k);
    return new_dbg;
}

fn get_unitigs(
    dbg: &FxHashMap<VarMer, DBGInfo>,
    k: usize,
    unitig_graph: bool,
) -> FxHashMap<VarMer, DBGInfo> {
    let mut unitigs = FxHashMap::default();
    let mut visited_nodes = FxHashSet::default();
    for (node, info) in dbg.iter() {
        if visited_nodes.contains(node) {
            continue;
        }
        if info.in_varmers.len() != 1
            || (info.in_varmers.len() == 1
                && dbg.get(&**info.in_varmers[0]).unwrap().out_varmers.len() != 1)
        {
            let mut cov = info.coverage;
            let mut count = 1;
            visited_nodes.insert(node);
            let mut unitig = node.clone();

            if info.out_varmers.len() == 1 {
                let mut prev = &node.clone();
                let mut next = &info.out_varmers[0];
                while dbg.get(&**prev).unwrap().out_varmers.len() == 1
                    && dbg.get(&**next).unwrap().in_varmers.len() == 1
                {
                    visited_nodes.insert(next);
                    if unitig_graph {
                        for pg in next.iter() {
                            if unitig.contains(pg) {
                                continue;
                            } else {
                                unitig.push(*pg);
                            }
                        }
                    } else {
                        unitig.push(next.last().unwrap().clone());
                    }

                    cov += dbg.get(&**next).unwrap().coverage;
                    count += 1;

                    if dbg.get(&**next).unwrap().out_varmers.len() == 0 {
                        break;
                    }
                    prev = next;
                    next = &dbg.get(&**next).unwrap().out_varmers[0];
                }
            }
            unitigs.insert(
                unitig,
                DBGInfo {
                    out_varmers: vec![],
                    in_varmers: vec![],
                    coverage: cov / count,
                },
            );
        }
    }
    get_edges_varmers(&mut unitigs, k);
    return unitigs;
}

fn count_kmers(dbg_frags: &Vec<FragDBG>, k: usize) -> (FxHashMap<VarMer, u64>, FxHashSet<u32>) {
    let mut snps = FxHashSet::default();
    let mut kmers = FxHashMap::default();
    for frag in dbg_frags.iter() {
        let seq = &frag.seq;
        let mut varmer = VecDeque::new();
        for (curr_pos, curr_geno) in seq.iter() {
            snps.insert(*curr_pos);
            varmer.push_back((*curr_pos, *curr_geno));
            if varmer.len() == k {
                let new_varmer = varmer.iter().cloned().collect::<VarMer>();
                *kmers.entry(new_varmer).or_insert(0) += 1;
                varmer.pop_front();
            }
        }
    }
    return (kmers, snps);
}

pub fn get_paths(
    unitigs: &FxHashMap<VarMer, DBGInfo>,
    conservative: bool,
) -> Vec<Vec<(&VarMer, usize)>> {
    let mut visited_nodes = FxHashSet::default();
    let mut paths = vec![];
    let mut sinks_and_forks = vec![];
    for (node, info) in unitigs.iter() {
        if info.in_varmers.len() == 0 {
            sinks_and_forks.push((node, info));
        } else if info.out_varmers.len() > 1 {
            sinks_and_forks.push((node, info));
        }
    }
    sinks_and_forks.sort_by(|a, b| b.1.coverage.cmp(&a.1.coverage));
    for (node, info) in sinks_and_forks.iter() {
        if info.out_varmers.len() == 0 {
            paths.push(vec![(*node, info.coverage as usize)]);
        }
        for out in info.out_varmers.iter() {
            let mut path = vec![(*node, info.coverage as usize)];
            let mut curr_node = &**out;
            loop {
                if visited_nodes.contains(curr_node) && conservative {
                    break;
                }
                visited_nodes.insert(curr_node.clone());
                let info = unitigs.get(&**curr_node).unwrap();
                path.push((curr_node, info.coverage as usize));
                let num_in_edges = info.in_varmers.len();
                if num_in_edges == 1 {
                    curr_node = &unitigs.get(&**curr_node).unwrap().out_varmers[0];
                } else {
                    paths.push(path);
                    break;
                }
            }
        }
    }
    return paths;
}

fn get_hits<'a>(
    varmer_d: &'a DictFrag,
    dict_frags: &'a Vec<DictFrag>,
    hit_threshold: usize,
    penalty_allow: usize,
    allow_equal: bool,
) -> Vec<Hit<'a>> {
    let mut hits = vec![];
    for ref_varmer in dict_frags.iter() {
        if varmer_d == ref_varmer && !allow_equal {
            continue;
        }

        let mut same = FxHashSet::default();
        let mut r_to_a = FxHashSet::default();
        let mut a_to_r = FxHashSet::default();
        let mut del = FxHashSet::default();

        //check if unitigs overlap based on first and last position
        if (varmer_d.first_position <= ref_varmer.last_position
            && varmer_d.last_position >= ref_varmer.first_position)
            || (ref_varmer.first_position <= varmer_d.last_position
                && ref_varmer.last_position >= varmer_d.first_position)
        {
            for (pos, geno) in ref_varmer.seq.iter() {
                if varmer_d.seq.contains_key(pos) {
                    let gn2 = *geno;
                    let gn1 = *varmer_d.seq.get(pos).unwrap();
                    if gn1 == gn2 {
                        same.insert(*pos);
                    } else if gn1 == 0 && gn2 != 0 {
                        r_to_a.insert(*pos);
                    } else if gn1 != 0 && gn2 == 0 {
                        a_to_r.insert(*pos);
                    } else {
                        a_to_r.insert(*pos);
                    }
                } else if *pos <= varmer_d.last_position && *pos >= varmer_d.first_position {
                    del.insert(*pos);
                }
            }

            for (pos, _) in varmer_d.seq.iter() {
                if !ref_varmer.seq.contains_key(pos)
                    && *pos <= ref_varmer.last_position
                    && *pos >= ref_varmer.first_position
                {
                    del.insert(*pos);
                }
            }
            if r_to_a.len() + a_to_r.len() + del.len() <= hit_threshold {
                if same.len() + penalty_allow > (r_to_a.len() + a_to_r.len()){
                    hits.push(Hit {
                        varmer: &ref_varmer,
                        same,
                        r_to_a,
                        a_to_r,
                        del,
                    });
                }
            }
        }
    }

    hits.sort_by(|a, b| {
        (a.varmer.first_position, a.varmer.last_position)
            .cmp(&(b.varmer.first_position, b.varmer.last_position))
    });
    return hits;
}

fn varmers_from_dp_res<'a>(dp_res: &'a DpResult, cutoff_ratio: f64) -> Vec<DictFrag> {
    if dp_res.total_errs as f64 > dp_res.total_matches as f64 * cutoff_ratio {
        return vec![];
    }
    if dp_res.traceback.len() == 0 {
        return vec![];
    }
    let mut varmers = vec![];
    let mut traceback = dp_res.max_index;
    while traceback != dp_res.traceback[traceback] {
        varmers.push(dp_res.hits[traceback].varmer);
        traceback = dp_res.traceback[traceback];
    }
    varmers.push(dp_res.hits[traceback].varmer);
    let varmers = varmers.into_iter().rev().cloned().collect::<Vec<_>>();
    return varmers;
}

fn dp_hits<'a>(
    hits: &'a Vec<Hit>,
    varmer_d: &DictFrag,
    graph: &FxHashMap<VarMer, DBGInfo>,
    threshold: usize,
    mismatch_pen: f64,
    ambiguous_allowed: bool,
    del_penalty: bool,
    require_dag: GraphConstraint,
    band: usize,
    k: usize,
    forbidden_nodes: &FxHashSet<&VarMer>,
) -> DpResult<'a> {
    let empty_result = DpResult {
        score: (0., 0),
        traceback: vec![],
        max_index: 0,
        hits: &hits,
        dels_max: FxHashSet::default(),
        rtoa_max: FxHashSet::default(),
        ator_max: FxHashSet::default(),
        same_max: FxHashSet::default(),
        total_errs: 0,
        total_matches: 0,
    };

    if hits.len() == 0 {
        return empty_result;
    }

    let mut dp_vec = hits
        .iter()
        .map(|x| {
            if del_penalty {
                (
                    x.same.len() as f64
                        + mismatch_pen * (x.r_to_a.len() + x.a_to_r.len() + x.del.len()) as f64,
                    x.varmer.cov,
                )
            } else {
                (
                    x.same.len() as f64 + mismatch_pen * (x.r_to_a.len() + x.a_to_r.len()) as f64,
                    x.varmer.cov,
                )
            }
        })
        .collect::<Vec<(f64, u64)>>();
    let mut traceback_vec = (0..hits.len()).map(|x| x).collect::<Vec<usize>>();
    //        let mut same_vecs : Vec<FxHashSet<u32>> = vec![FxHashSet::default(); hits.len()];
    //        let mut rtoa_vecs : Vec<FxHashSet<u32>> = vec![FxHashSet::default(); hits.len()];
    //        let mut ator_vecs : Vec<FxHashSet<u32>> = vec![FxHashSet::default(); hits.len()];
    //        let mut del_vecs : Vec<FxHashSet<u32>> = vec![FxHashSet::default(); hits.len()];
    let mut same_vecs = hits
        .iter()
        .map(|x| x.same.clone())
        .collect::<Vec<FxHashSet<u32>>>();
    let mut rtoa_vecs = hits
        .iter()
        .map(|x| x.r_to_a.clone())
        .collect::<Vec<FxHashSet<u32>>>();
    let mut ator_vecs = hits
        .iter()
        .map(|x| x.a_to_r.clone())
        .collect::<Vec<FxHashSet<u32>>>();
    let mut del_vecs = hits
        .iter()
        .map(|x| x.del.clone())
        .collect::<Vec<FxHashSet<u32>>>();
    let mut chain_length = vec![1; hits.len()];

    for i in 0..hits.len() {
        if forbidden_nodes.contains(&hits[i].varmer.seq_vec) {
            continue;
        }
        let mut best_index = i;
        //let mut new_scores = vec![0.; i];
        let mut best_score = 0.;
        let mut last_tied = false;

        if varmer_d.last_position == 234 && varmer_d.first_position == 210 {
            let hit = &hits[i];
            log::trace!(
                "TEST HIT HIT: {:?}, score {}, bad {}, cov {}",
                hit.varmer.seq_vec,
                hit.same.len(),
                hit.r_to_a.len() + hit.a_to_r.len() + hit.del.len(),
                dp_vec[i].1
            );
        }

        let ind = if i >= band { i - band } else { 0 };

        for j in ind..i {
            let cov_j = hits[j].varmer.cov;
            let mut score = 0;
            let mut bad = 0;
            let mut inside = false;
            if forbidden_nodes.contains(&hits[j].varmer.seq_vec) {
                continue;
            }
            if require_dag == GraphConstraint::RequireDag{
                for varmer in graph[&hits[i].varmer.seq_vec].in_varmers.iter() {
                    if &**varmer == &hits[j].varmer.seq_vec {
                        inside = true;
                        break;
                    }
                }
            } else if require_dag == GraphConstraint::NoRequireDag {
                if hits[j].varmer.last_position  < hits[i].varmer.first_position + k as u32{
                    inside = true;
                }
            }
            else if require_dag == GraphConstraint::RequireDagRescue{
                if graph[&hits[j].varmer.seq_vec].out_varmers.len() == 0 &&
                 graph[&hits[i].varmer.seq_vec].in_varmers.len() == 0 {
                    if hits[j].varmer.last_position < hits[i].varmer.first_position  + k as u32 {
                        inside = true;
                    }
                 }
                else{
                    for varmer in graph[&hits[i].varmer.seq_vec].in_varmers.iter() {
                        if &**varmer == &hits[j].varmer.seq_vec {
                            inside = true;
                            break;
                        }
                    }
                }
            }
            else{
                panic!("Invalid graph constraint");
            }

            if !inside {
                continue;
            }

            for pos in hits[i].same.iter() {
                if !same_vecs[j].contains(pos) {
                    score += 1;
                }
            }
            for pos in hits[i].r_to_a.iter() {
                if !rtoa_vecs[j].contains(pos) {
                    bad += 1;
                }
            }
            for pos in hits[i].a_to_r.iter() {
                if !ator_vecs[j].contains(pos) {
                    bad += 1;
                }
            }
            if del_penalty {
                for pos in hits[i].del.iter() {
                    if !del_vecs[j].contains(pos) {
                        bad += 1;
                    }
                }
            }

//            if varmer_d.last_position == 142 && varmer_d.first_position == 1 {
//                log::trace!(
//                    "TEST HIT HIT j: {:?} best_index {} new scores {:?}, bad {}, score {}, INSIDE {}",
//                    hits[j].varmer.seq_vec,
//                    best_index,
//                    best_score,
//                    bad,
//                    score,
//                    inside
//                );
//            }
//
            if bad + rtoa_vecs[j].len() + ator_vecs[j].len() + del_vecs[j].len() <= threshold {
                let put_score = score as f64 + mismatch_pen * bad as f64 + dp_vec[j].0;
                if best_score < put_score || (best_score == put_score && cov_j > hits[best_index].varmer.cov) {
                    best_index = j;
                    best_score = put_score;
                    last_tied = false;
                } 
            }
        }
        if varmer_d.last_position == 60 && varmer_d.first_position == 38 {
            log::trace!("best_index {}, i {}, score,cov {} {} ", best_index, i, best_score, dp_vec[best_index].1);
        }

        let j = best_index;
        if j == i || last_tied {
            continue;
        }
        traceback_vec[i] = j;
        chain_length[i] = chain_length[j] + 1;
        rtoa_vecs[i] = rtoa_vecs[i].union(&rtoa_vecs[j]).cloned().collect();
        ator_vecs[i] = ator_vecs[i].union(&ator_vecs[j]).cloned().collect();
        del_vecs[i] = del_vecs[i].union(&del_vecs[j]).cloned().collect();
        same_vecs[i] = same_vecs[i].union(&same_vecs[j]).cloned().collect();
        //            ator_vecs[i] = ator_vecs[i].union(&hits[j].3).cloned().collect();
        //            del_vecs[i] = del_vecs[i].union(&hits[j].4).cloned().collect();
        //            same_vecs[i] = same_vecs[i].union(&hits[j].1).cloned().collect();
        dp_vec[i] = (best_score, ((dp_vec[i].1 * chain_length[i] - 1) + dp_vec[j].1) / chain_length[i]);
    }

    //get traceback with max score
    //sort the vector and remember the index
    let mut score_vec = dp_vec
        .iter()
        .enumerate()
        .map(|(i, x)| (x.0, x.1, i))
        .collect::<Vec<(f64, u64, usize)>>();
    score_vec.sort_by(|a, b| b.partial_cmp(&a).unwrap());
    let max_index = score_vec[0].2;
    let max_score = (score_vec[0].0, score_vec[0].1);

    //When mapping reads, don't allow ambiguous paths.
    if !ambiguous_allowed {
        if score_vec.len() > 1 && score_vec[1].0 == max_score.0 {
            log::trace!("AMBIGUOUS PATHS. MAX SCORES and INDICES ARE");
            for i in 0..score_vec.len().min(5){
                log::trace!("{}: {} {}", i, score_vec[i].0, score_vec[i].2);
                print_varmer_d(hits[score_vec[i].2].varmer, true);
            }
            return empty_result;
        }
    }

    return DpResult {
        score: max_score,
        traceback: traceback_vec,
        max_index,
        hits,
        dels_max: del_vecs[max_index].clone(),
        rtoa_max: rtoa_vecs[max_index].clone(),
        ator_max: ator_vecs[max_index].clone(),
        same_max: same_vecs[max_index].clone(),
        total_errs: del_vecs[max_index].len()
            + rtoa_vecs[max_index].len()
            + ator_vecs[max_index].len(),
        total_matches: same_vecs[max_index].len(),
    };
}

fn fmt_varmer_d(varmer_d: &DictFrag) -> String {
    let mut seq1 = varmer_d.seq.iter().collect::<Vec<(&u32, &u8)>>();
    seq1.sort_by(|a, b| a.0.cmp(&b.0));
    let mut seq1 = seq1
        .iter()
        .map(|(_, geno)| format!("{}", *geno))
        .collect::<Vec<String>>()
        .join("");
    seq1.push_str(&format!(
        "_{}-{}--{}-LEN{}",
        varmer_d.first_position,
        varmer_d.last_position,
        varmer_d.cov,
        varmer_d.seq.len()
    ));
    seq1
}

fn print_varmer(varmer: &VarMer, trace: bool) {
    let mut seq1 = varmer
        .iter()
        .map(|(_, geno)| format!("{}", geno))
        .collect::<Vec<String>>()
        .join("");
    seq1.push('_');
    seq1.push_str(&format!(
        "{}-{}",
        varmer.first().unwrap().0,
        varmer.last().unwrap().0
    ));
    if trace {
        log::trace!("{}", seq1);
    } else {
        log::debug!("{}", seq1);
    }
}

fn print_varmer_d(varmer_d: &DictFrag, trace: bool) {
    let seq1 = fmt_varmer_d(varmer_d);

    if trace {
        log::trace!("{}", seq1);
    } else {
        log::debug!("{}", seq1);
    }
}

pub fn query_unitigs(unitigs: &FxHashMap<VarMer, DBGInfo>, threshold: usize) -> Vec<VarMer> {
    let mut dict_unitigs = vec![];
    let mut bad_unitigs = vec![];
    for unitig in unitigs.iter() {
        let mut dict_seq = FxHashMap::default();
        let cov = unitig.1.coverage;
        let mut first_position = u32::MAX;
        let mut last_position = 0;
        for (pos, geno) in unitig.0.iter() {
            dict_seq.insert(*pos, *geno);
            if *pos < first_position {
                first_position = *pos;
            }
            if *pos > last_position {
                last_position = *pos;
            }
        }
        dict_unitigs.push(DictFrag {
            seq: dict_seq,
            seq_vec: unitig.0.clone(),
            first_position,
            last_position,
            cov,
        });
    }

    dict_unitigs.sort_by(|a, b| a.cov.cmp(&b.cov));  

    let mut failed_unitigs = FxHashSet::default();
    for unitig1 in dict_unitigs.iter() {
        let hits = get_hits(unitig1, &dict_unitigs, threshold, usize::MAX/2, false);

        log::trace!("QUERY");
        print_varmer_d(unitig1, true);

        if hits.len() == 0 {
            continue;
        }

        let dp_res = dp_hits(
            &hits,
            unitig1,
            unitigs,
            threshold,
            0.2,
            true,
            true,
            GraphConstraint::RequireDag,
            usize::MAX,
            0,
            &failed_unitigs
        );

        if dp_res.total_errs > threshold {
            log::trace!("FAILED ERRS: {}", dp_res.total_errs);
            continue;
        }

        if dp_res.total_errs + dp_res.same_max.len() < unitig1.seq.len() as usize {
            log::trace!("FAILED LEN: {}", dp_res.total_errs + dp_res.same_max.len());
            log::trace!(
                "VECS {:?} {:?} {:?} {:?}",
                dp_res.rtoa_max,
                dp_res.ator_max,
                dp_res.dels_max,
                dp_res.same_max
            );
            log::trace!("BEST HIT: {:?}", hits[dp_res.max_index].varmer.seq_vec);
            log::trace!("ITSELF {:?}", unitig1.seq_vec);
            continue;
        }

        log::trace!("PASSED");
        print_varmer_d(unitig1, true);
        log::trace!(
            "MAX SCORE: {:?}, RTOA {:?}, ATOR {:?}, DEL {:?}, max_index {}",
            dp_res.score,
            dp_res.rtoa_max,
            dp_res.ator_max,
            dp_res.dels_max,
            dp_res.max_index
        );

        log::trace!("HIT");
        print_varmer_d(hits[dp_res.max_index].varmer, true);
        let max_index = dp_res.max_index;
        let traceback_vec = &dp_res.traceback;

        let mut avg_cov = hits[max_index].varmer.cov;
        let mut count = 1;
        let mut final_cov = 0;

        if hits[max_index].varmer.first_position <= unitig1.first_position
            && hits[max_index].varmer.last_position >= unitig1.last_position
        {
            final_cov = hits[max_index].varmer.cov;
        }

        let mut trace_index = max_index;

        while trace_index != traceback_vec[trace_index] {
            trace_index = traceback_vec[trace_index];
            log::trace!("HIT");
            print_varmer_d(hits[trace_index].varmer, true);
            avg_cov += hits[trace_index].varmer.cov;
            count += 1;

            if hits[trace_index].varmer.first_position <= unitig1.first_position
                && hits[trace_index].varmer.last_position >= unitig1.last_position
            {
                if hits[trace_index].varmer.cov > final_cov{
                    final_cov = hits[trace_index].varmer.cov;
                }
            }
        }

        if final_cov == 0 {
            final_cov = avg_cov / count;
        }
        log::trace!("QUERY UNITIGS: FINAL COV: {}, contig_COV {}", final_cov, unitig1.cov);

        let mut failed = false;
        let num_errs = dp_res.total_errs as f64;

        let mult;
        if dp_res.rtoa_max.len() == 0 && dp_res.ator_max.len() == 0 {
            mult = 0.35.powi(dp_res.dels_max.len() as i32);
        } else {
            mult = 0.15.powi(dp_res.rtoa_max.len() as i32)
                * 0.10.powi(dp_res.ator_max.len() as i32);
        }

        if num_errs > 0.{
            if binomial_test(
                final_cov as u64,
                unitig1.cov as u64,
                mult,
            ) > 0.005
            {
                failed = true;
                log::trace!("BAD ERR: {}", final_cov);
                print_varmer_d(unitig1, true);
                print_varmer_d(hits[dp_res.max_index].varmer, true);
            }
        }

        if failed {
            failed_unitigs.insert(&unitig1.seq_vec);
            bad_unitigs.push(unitig1.seq_vec.clone());
        }
    }
    return bad_unitigs;
}

fn fragdbg_to_dictfrag(frag: &FragDBG) -> DictFrag {
    let mut seq = FxHashMap::default();
    for (pos, geno) in frag.seq.iter() {
        seq.insert(*pos, *geno);
    }
    DictFrag {
        seq,
        seq_vec: frag.seq.iter().cloned().collect(),
        first_position: frag.first_position,
        last_position: frag.last_position,
        cov: frag.seq.len() as u64,
    }
}

fn binomial_test(n: u64, k: u64, p: f64) -> f64 {
    // n: number of trials
    // k: number of successes
    // p: probability of success

    // Create a binomial distribution
    let binomial = Binomial::new(p, n).unwrap();

    // Calculate the probability of observing k or more successes
    let p_value = 1.0 - binomial.cdf(k);

    p_value
}

fn get_path_haps<'a>(
    dbg_frags: &'a Vec<FragDBG>,
    final_unitigs: &FxHashMap<VarMer, DBGInfo>,
    paths: Vec<Vec<(&VarMer, usize)>>,
    _snps: usize,
    options: &Options,
) -> Vec<HapFinalResult<'a>> {
    let mut path_frags: Vec<DictFrag> = vec![];
    for path in paths {
        let mut seq = vec![];
        let mut covs = vec![];
        for (varmer, cov) in path.iter() {
            for (pos, geno) in varmer.iter() {
                if !seq.contains(&(*pos, *geno)) {
                    seq.push((*pos, *geno));
                }
            }
            covs.push(*cov);
        }
        let final_cov = covs.iter().fold(0, |acc, cov| acc + cov) / covs.len();
        let mut s = seq
            .iter()
            .map(|(_, geno)| format!("{}", geno))
            .collect::<Vec<String>>()
            .join("");
        s.push_str(&format!(
            "_{}-{}:{}",
            seq.first().unwrap().0,
            seq.last().unwrap().0,
            final_cov
        ));
        log::trace!("GETTING PATH HAP: {}", s);
        for (varmer, _) in path.iter() {
            print_varmer(varmer, true);
            log::trace!("COV: {}", final_unitigs[*varmer].coverage);
        }
        let path_frag = DictFrag {
            seq: seq.iter().cloned().collect::<FxHashMap<u32, u8>>(),
            first_position: seq.first().unwrap().0,
            last_position: seq.last().unwrap().0,
            seq_vec: seq,
            cov: 0,
        };
        path_frags.push(path_frag);
    }
    path_frags.sort_by(|a, b| {
        (a.first_position, a.last_position).cmp(&(b.first_position, b.last_position))
    });
    let mut final_results: Vec<HapFinalResult> = vec![];
    for path_frag in path_frags.iter() {
        let final_res = HapFinalResult {
            relative_abundances: 0.,
            depth: 0.,
            assigned_frags: vec![],
            path_frag: path_frag.clone(),
        };
        final_results.push(final_res);
    }

    reassign_frags(dbg_frags, &mut final_results, false);
    //print relative percentage
    let total_cov = final_results
        .iter()
        .map(|x| &x.assigned_frags)
        .fold(0, |acc, assignment| {
            acc + assignment.iter().fold(0, |acc, frag| acc + frag.seq.len())
        });

    let mut ret_results = vec![];
    for (i,mut final_res) in final_results.into_iter().enumerate(){
        let mut relative_cov = 0;
        let mut depth = 0;
        for frag in final_res.assigned_frags.iter() {
            relative_cov += frag.seq.len();
            depth += frag.seq.len().min(final_res.path_frag.seq.len());
        }

        let abundance = relative_cov as f64 / total_cov as f64 * 100.;
        let depth = depth as f64 / final_res.path_frag.seq.len() as f64;

        log::debug!(
            "PATH: {}, COV: {}, RELATIVE COV: {}",
            i,
            depth,
            abundance
        );

        final_res.depth = depth;
        final_res.relative_abundances = abundance;

        if abundance > options.min_abund && depth > 1. {
            ret_results.push(final_res);
        }

    }

    return ret_results;
}

fn reassign_frags<'a>(dbg_frags: &'a Vec<FragDBG>, final_results: &mut Vec<HapFinalResult<'a>>, assign_ambiguous: bool) -> Vec<&'a FragDBG>{
    let mut assignments = vec![vec![]; final_results.len()];
    let mut unassignable = vec![];
    let seq_lens = final_results
        .iter()
        .map(|x| x.path_frag.seq.len())
        .collect::<Vec<usize>>();
    let seq_covs = final_results
        .iter()
        .map(|x| x.relative_abundances)
        .collect::<Vec<f64>>();


    for frag in dbg_frags {
        let mut ambig = false;
        let mut best_score = 0;
        let mut best_index = 0;
        for (i, res) in final_results.iter_mut().enumerate() {
            let mut score = 0;
            for (pos, geno) in frag.seq.iter() {
                if res.path_frag.seq.contains_key(pos) {
                    if *geno == *res.path_frag.seq.get(pos).unwrap() {
                        score += 2;
                    } else {
                        // if fragments genotype is 0, possible reference bias
                        if *geno == 0{
                            score -= 3
                        }
                        else{
                            score -= 5
                        }
                    }
                } else {
                    score -= 1;
                }
            }
            if score > best_score {
                best_score = score;
                best_index = i;
                ambig = false;
            } else if score == best_score {
                ambig = true;
                if res.path_frag.seq.len() as f64 * seq_covs[i] > seq_lens[best_index] as f64 * seq_covs[best_index] {
                    best_index = i;
                }
            }
        }
        if best_score == 0 {
            unassignable.push(frag);
            continue;
        }
        if ambig && !assign_ambiguous {
            unassignable.push(frag);
            continue;
        }

        assignments[best_index].push(frag);
    }
    for (i, assignment) in assignments.into_iter().enumerate() {
        final_results[i].assigned_frags = assignment;
    }
    return unassignable;
}

fn print_final_hap_results(
    final_results: &Vec<HapFinalResult>,
    snps: usize,
    options: &Options,
    hap_file: &str,
    id_file: &str,
    fastq_file: Option<&str>,
    contig_range: (&str, Option<(usize,usize)>),
    unassigned: Option<&Vec<&FragDBG>>,
    vcf_profile: &VcfProfile,
    snp_pos_to_genome_pos: &Vec<usize>
) {
    //prepend outdir_dir
    let contig_name = contig_range.0;
    let start;
    let end;
    if let Some((s, e)) = contig_range.1 {
        // start as string
        start = format!("{}", s);
        end = format!("{}", e);
    } else {
        start = String::from("ALL");
        end = String::from("ALL");
    }
    let dir = Path::new(&options.output_dir);
    let hap_file = dir.join(hap_file);
    let hap_file = hap_file.to_str().unwrap();
    let id_file = dir.join(id_file);
    let id_file = id_file.to_str().unwrap();
    let mut haplotype_writer;
    let mut id_writer;
    let mut fastq_writer = None;
    if Path::exists(Path::new(hap_file)) {
        haplotype_writer = BufWriter::new(
            std::fs::File::options()
                .append(true)
                .open(hap_file)
                .expect(&format!("Could not open haplotype file {}", hap_file)),
        );
        id_writer = BufWriter::new(
            std::fs::File::options()
                .append(true)
                .open(id_file)
                .expect(&format!("Could not open id file {}", id_file)),
        );
        if let Some(fastq_file) = fastq_file {
            let fastq_file = dir.join(fastq_file);
            let fastq_file = fastq_file.to_str().unwrap();

            fastq_writer = Some(BufWriter::new(
                std::fs::File::options()
                    .append(true)
                    .open(fastq_file)
                    .expect(&format!("Could not open fastq file {}", fastq_file)),
            ));
        }
    } else {
        haplotype_writer = BufWriter::new(
            std::fs::File::create(hap_file)
                .expect(&format!("Could not create haplotype file {}", hap_file)),
        );
        id_writer = BufWriter::new(
            std::fs::File::create(id_file).expect(&format!("Could not create id file {}", id_file)),
        );
        if let Some(fastq_file) = fastq_file {
            let fastq_file = dir.join(fastq_file);
            let fastq_file = fastq_file.to_str().unwrap();

            fastq_writer = Some(BufWriter::new(
                std::fs::File::create(fastq_file)
                    .expect(&format!("Could not create fastq file {}", fastq_file)),
            ));
        }
    }

    for (i, res) in final_results.iter().enumerate() {
        haplotype_writer
            .write_all(
                format!(
                    ">Contig:{},Range:{}-{},Haplotype:{},Abundance:{:.2},Depth:{:.2}\n",
                    contig_name, start,end, i, res.relative_abundances, res.depth
                )
                .as_bytes(),
            )
            .unwrap();
        let seq = &res.path_frag.seq;
        let mut printable_seq = vec![b'-'; snps];
        if options.allele_output{
            let pos_to_allele = &vcf_profile.vcf_pos_allele_map[contig_name];
            for (pos, geno) in seq.iter() {
                let seq_pos = snp_pos_to_genome_pos[*pos as usize - 1];
                let allele = pos_to_allele[&seq_pos][*geno as usize];
                printable_seq[*pos as usize - 1] = allele;
            }
        }
        else{
            for (pos, geno) in seq.iter() {
                printable_seq[*pos as usize - 1] = *geno + 48;
            }
        }
        //print wrapped lines of 80
        let mut j = 0;
        while j < snps {
            let end = std::cmp::min(j + 80, snps);
            haplotype_writer.write_all(&printable_seq[j..end]).unwrap();
            haplotype_writer.write_all(b"\n").unwrap();
            j += 80;
        }
        haplotype_writer.write_all(b"\n").unwrap();

        //print ids to a file where each row is a path and each column is a frag id, tab sep
        id_writer
            .write_all(format!("Contig:{}\tRange:{}-{}\tHaplotype:{}\t", contig_name,start,end, i).as_bytes())
            .unwrap();
        for frag in res.assigned_frags.iter() {
            id_writer.write_all(frag.id.as_bytes()).unwrap();
            id_writer.write_all(b"\t").unwrap();
        }
        id_writer.write_all(b"\n").unwrap();

        //Write these seq strings to a fastq file, diff identifier
        let reads = &res.assigned_frags;
        let seqs = reads  
            .iter()
            .map(|frag| frag.seq_string[0].to_ascii_vec())
            .collect::<Vec<Vec<u8>>>();
        if let Some(fastq_writer) = &mut fastq_writer {
            for (j, seq) in seqs.iter().enumerate() {
                if j != 0 && i != 0 {
                    fastq_writer.write_all(b"\n").unwrap();
                }
                let rec_str = format!("@Contig:{},Range:{}-{},Haplotype:{},Read:{},\n", contig_name, start,end, i, reads[j].id);
                fastq_writer.write_all(rec_str.as_bytes()).unwrap();
                fastq_writer.write_all(seq).unwrap();
                fastq_writer.write_all(b"\n+\n").unwrap();
                if reads[j].qual_string[0].len() != reads[j].seq_string[0].len() {
                    fastq_writer.write_all(&vec![b'I'; reads[j].seq_string[0].len()]).unwrap();
                } else {
                    fastq_writer.write_all(&reads[j].qual_string[0]).unwrap();
                }
                fastq_writer.write_all(b"\n").unwrap();
            }
        }
    }

    if let Some(unassigned) = unassigned {
        id_writer.write_all(format!("Contig:{}\tRange:{}-{}\tHaplotype:unassigned\t", contig_name, start, end).as_bytes()).unwrap();
        for frag in unassigned.iter() {
            id_writer.write_all(frag.id.as_bytes()).unwrap();
            id_writer.write_all(b"\t").unwrap();

            if let Some(fastq_writer) = &mut fastq_writer{
                let rec_str = format!("@Contig:{},Range:{}-{},Haplotype:unassigned,Read:{},\n", contig_name, start, end, frag.id);
                fastq_writer.write_all(rec_str.as_bytes()).unwrap();
                fastq_writer.write_all(&frag.seq_string[0].to_ascii_vec()).unwrap();
                fastq_writer.write_all(b"\n+\n").unwrap();
                if frag.qual_string[0].len() != frag.seq_string[0].len() {
                    fastq_writer.write_all(&vec![b'I'; frag.seq_string[0].len()]).unwrap();
                } else {
                    fastq_writer.write_all(&frag.qual_string[0]).unwrap();
                }
                fastq_writer.write_all(b"\n").unwrap();
            }
        }
    }
}

fn clean_hanging_kmers(
    unitigs: FxHashMap<VarMer, DBGInfo>,
    k: usize,
) -> FxHashMap<VarMer, DBGInfo> {
    let mut int_rep: FxHashMap<usize, (Vec<usize>, Vec<usize>, u64)> = FxHashMap::default();
    let unitigs_to_ints: FxHashMap<&VarMer, usize> = unitigs
        .iter()
        .enumerate()
        .map(|(i, (varmer, _))| (varmer, i))
        .collect();

    for (varmer, info) in unitigs.iter() {
        let mut in_ints = vec![];
        let mut out_ints = vec![];
        for in_varmer in info.in_varmers.iter() {
            in_ints.push(*unitigs_to_ints.get(&**in_varmer).unwrap());
        }
        for out_varmer in info.out_varmers.iter() {
            out_ints.push(*unitigs_to_ints.get(&**out_varmer).unwrap());
        }
        int_rep.insert(
            unitigs_to_ints.get(varmer).unwrap().clone(),
            (in_ints, out_ints, info.coverage),
        );
    }

    let mut new_unitigs_map = FxHashMap::default();

    let mut new_unitigs = FxHashMap::default();
    let k_l = k / 2;
    let k_r = k - k_l;

    for (varmer, info) in unitigs.iter() {
        let ind = unitigs_to_ints.get(varmer).unwrap();
        let cut_varmer;
        if info.in_varmers.len() == 0 && info.out_varmers.len() == 0 {
            cut_varmer = varmer.clone();
        } else if info.in_varmers.len() != 0 && info.out_varmers.len() != 0 {
            cut_varmer = varmer[k_l..varmer.len() - k_r]
                .iter()
                .cloned()
                .collect::<VarMer>();
        } else if info.in_varmers.len() == 0 && info.out_varmers.len() != 0 {
            cut_varmer = varmer[..varmer.len() - k_r]
                .iter()
                .cloned()
                .collect::<VarMer>();
        } else if info.in_varmers.len() != 0 && info.out_varmers.len() == 0 {
            cut_varmer = varmer[k_l..].iter().cloned().collect::<VarMer>();
        } else {
            cut_varmer = varmer.clone();
        }
        new_unitigs_map.insert(ind, cut_varmer.clone());
        new_unitigs.insert(cut_varmer, DBGInfo::default());
    }

    for (ind, (in_ints, out_ints, cov)) in int_rep.iter() {
        let mut new_in_ints = vec![];
        let mut new_out_ints = vec![];
        for in_int in in_ints.iter() {
            if new_unitigs_map.contains_key(in_int) {
                new_in_ints.push(*in_int);
            }
        }
        for out_int in out_ints.iter() {
            if new_unitigs_map.contains_key(out_int) {
                new_out_ints.push(*out_int);
            }
        }
        let varmer = new_unitigs_map.get(ind).unwrap().clone();
        new_unitigs.get_mut(&varmer).unwrap().in_varmers = new_in_ints
            .iter()
            .map(|x| Arc::new(new_unitigs_map.get(x).unwrap().clone()))
            .collect();
        new_unitigs.get_mut(&varmer).unwrap().out_varmers = new_out_ints
            .iter()
            .map(|x| Arc::new(new_unitigs_map.get(x).unwrap().clone()))
            .collect();
        new_unitigs.get_mut(&varmer).unwrap().coverage = *cov;
    }

    return new_unitigs;
}

fn get_outside_paths_and_integers(
    unitig_paths: &Vec<VarmerPath>,
    dict_frag_to_index: &FxHashMap<&DictFrag, usize>,
) -> Vec<VarmerPathInteger> {
    //1) get outside paths: paths not contained in any reads

    let unitig_paths = unitig_paths
        .iter()
        .filter(|x| x.varmers.len() > 0)
        .collect::<Vec<&VarmerPath>>();
    for path in unitig_paths.iter() {
        log::trace!("UNITIG PATH: COV {}", path.total_avg_cov);
        for varmer in path.varmers.iter() {
            print_varmer_d(varmer, true);
        }
    }
    type Iv<'a> = Interval<u32, usize>;
    let data = unitig_paths
        .iter().enumerate()
        .map(|(y,x)| Iv{start: x.first - 1,
            stop: x.last,
            val: y})
        .collect::<Vec<Iv>>();
    let lapper = Lapper::new(data);

    let mut graph_dict = FxHashMap::default();
    for path in unitig_paths.iter() {
        graph_dict.entry(path).or_insert((vec![], vec![]));
    }
    let hashsets_paths = unitig_paths
        .iter()
        .map(|x| x.varmers.iter().collect::<FxHashSet<_>>())
        .collect::<Vec<FxHashSet<_>>>();
    for i in 0..unitig_paths.len() {
        let overlaps = lapper.find(unitig_paths[i].first, unitig_paths[i].last);
        for ol in overlaps {
            let j = ol.val;
            if i == j {
                continue;
            }
            let path1 = &unitig_paths[i];
            let path2 = &unitig_paths[j];
            let mut contained = false;

            if hashsets_paths[i].is_subset(&hashsets_paths[j]) {
                contained = true;
            }

            if contained {
                let (out_edges, _) = graph_dict.get_mut(path1).unwrap();
                out_edges.push(path2);
                let (_, in_edges) = graph_dict.get_mut(path2).unwrap();
                in_edges.push(path1);
            }
        }
    }

    let mut outside_paths = FxHashMap::default();
    for path in unitig_paths.iter() {
        let (out_edges, _) = graph_dict.get(path).unwrap();
        if out_edges.len() == 0 {
            outside_paths.insert(path.varmers.clone(), path.total_avg_cov);
        }
    }

    // Distribute coverage for unambiguous contained paths
    for path in unitig_paths.iter() {
        let (out_edges, _) = graph_dict.get(&path).unwrap();
        let mut outside_out_edges = vec![];
        for path in out_edges.iter() {
            if outside_paths.contains_key(&path.varmers) {
                outside_out_edges.push(path);
            }
        }
        if outside_out_edges.len() != 1 {
            continue;
        }
        for outside_path in outside_out_edges.iter() {
            let len = outside_out_edges.len() as u64;
            let cov = outside_paths.get_mut(&outside_path.varmers).unwrap();
            *cov += path.total_avg_cov / len;
        }
    }

    for path in outside_paths.iter() {
        log::trace!("OUTSIDE PATH: COV {}", outside_paths[path.0]);
        for varmer in path.0.iter() {
            print_varmer_d(varmer, true);
        }
    }

    //Turn into "integer" paths with varmers as alleles and positions, where positions are
    //arbitrary but alleles are unitigs
    let outside_paths = outside_paths
        .into_iter()
        .map(|(varmers, cov)| VarmerPath {
            first: varmers.first().unwrap().first_position,
            last: varmers.last().unwrap().last_position,
            varmers,
            total_avg_cov: cov,
        })
        .collect::<Vec<VarmerPath>>();
    let mut integer_paths = vec![];
    for path in outside_paths.iter() {
        let integer_node = path
            .varmers
            .iter()
            .map(|varmer| (*dict_frag_to_index.get(varmer).unwrap() as u32, *dict_frag_to_index.get(varmer).unwrap() as u8))
            .collect::<Vec<(u32, u8)>>();
        let integer_path = VarmerPathInteger {
            first: path.first,
            last: path.last,
            intver: integer_node,
            total_avg_cov: path.total_avg_cov,
        };
        integer_paths.push(integer_path);
    }

    return integer_paths;
}

fn get_assembly_integer_graph(
    integer_paths: &Vec<VarmerPathInteger>,
) -> FxHashMap<VarMer, DBGInfo> {
    let mut assembly_graph: FxHashMap<VarMer, DBGInfo> = FxHashMap::default();
    for path1 in integer_paths.iter() {
        assembly_graph
            .entry(path1.intver.clone())
            .or_insert(DBGInfo {
                out_varmers: vec![],
                in_varmers: vec![],
                coverage: path1.total_avg_cov,
            });
        log::trace!("INTEGER PATH - FIRST {}, LAST {}, total_avg_cov {}", path1.first, path1.last, path1.total_avg_cov);
        let mut string = String::new();
        for integer in path1.intver.iter() {
            string.push_str(&format!(" -> {}", integer.0));
        }
        log::trace!("{}", string);
    }
    let data = integer_paths
        .iter().enumerate()
        .map(|(y,x)| Interval{start: x.first - 1,
            stop: x.last,
            val: y})
        .collect::<Vec<Interval<u32, usize>>>();
    let lapper = Lapper::new(data);
    for path1 in integer_paths.iter() {
        let overlaps = lapper.find(path1.first, path1.last);
        for ol in overlaps {
            let path2 = &integer_paths[ol.val];
            if path1.intver == path2.intver {
                continue;
            }
            log::trace!("PATH 1 {} {}", path1.first, path1.last);
            log::trace!("PATH 2 {} {}", path2.first, path2.last);
            //check for suffix prefix overlaps of VarmerPaths.varmers
            let mut overlap_len = 0;
            for k in 0..path1.intver.len().min(path2.intver.len()) {
                if path1.intver[k..] == path2.intver[..path2.intver.len() - k] {
                    overlap_len = k;
                    break;
                }
            }
            if overlap_len > 0 {
                log::trace!("OVERLAP LEN {}", overlap_len);
                let info1 = assembly_graph.get_mut(&path1.intver).unwrap();
                info1.out_varmers.push(Arc::new(path2.intver.clone()));
                let info2 = assembly_graph.get_mut(&path2.intver).unwrap();
                info2.in_varmers.push(Arc::new(path1.intver.clone()));
            }
        }
    }
    return assembly_graph;
}

fn subsample_positions_fragdbg(
    fragdbg: &mut Vec<FragDBG>,
    positions: &FxHashSet<u32>,
    old_pos_to_index_map: &FxHashMap<u32, u32>,
    old_snp_pos_to_gn_pos: &Vec<usize>
) -> Vec<usize>{
    for frag in fragdbg.iter_mut(){
        let mut new_seq = vec![];
        let mut new_seq_dict = FxHashMap::default();
        let mut new_snp_pos_to_seq_pos = FxHashMap::default();
        let mut new_qual_dict = FxHashMap::default();

        for (pos, geno) in frag.seq.iter() {
            if positions.contains(pos) {
                let new_pos = *old_pos_to_index_map.get(pos).unwrap();
                new_seq.push((new_pos, *geno));
                new_seq_dict.insert(new_pos, *geno);
                new_snp_pos_to_seq_pos.insert(new_pos, frag.snp_pos_to_seq_pos[pos]);
                new_qual_dict.insert(new_pos, frag.qual_dict[pos]);
            }
        }

        let new_first_pos;
        let new_last_pos;
        if new_seq.len() > 0 {
            new_first_pos = new_seq.first().unwrap().0;
            new_last_pos = new_seq.last().unwrap().0;
        }
        else{
            new_first_pos = 0;
            new_last_pos = 0;
        }

        frag.seq = new_seq;
        frag.seq_dict = new_seq_dict;
        frag.qual_dict = new_qual_dict;
        frag.snp_pos_to_seq_pos = new_snp_pos_to_seq_pos;
        frag.first_position = new_first_pos;
        frag.last_position = new_last_pos;
    }
    let mut new_snp_pos_to_gn_pos = vec![];
    for (i, pos) in old_snp_pos_to_gn_pos.iter().enumerate(){
        if positions.contains(&(i as u32 + 1)){
            new_snp_pos_to_gn_pos.push(*pos);
        }
    }

    return new_snp_pos_to_gn_pos
}

fn remove_tips(
    unitigs: &FxHashMap<VarMer, DBGInfo>,
    _options: &Options,
    k: usize,
) -> Vec<VarMer> {

    let mut bad_unitigs = vec![];
    //Prev-tip
    for varmer in unitigs.keys(){
        //Only pop tips for deletions
        if varmer.len() == varmer.last().unwrap().0 as usize - varmer[0].0 as usize + 1{
            continue
        }
        let into = &unitigs[varmer].in_varmers;
        let out = &unitigs[varmer].out_varmers;
        let cov = unitigs[varmer].coverage;
        if into.len() == 0 && out.len() == 1{
            let test_cov = unitigs[out[0].as_ref()].coverage;
            //only goes 1 k-mer back
            if varmer[0].0 + k as u32 > out[0][0].0{
                if binomial_test(test_cov, cov, 0.10) > 0.005{
                    bad_unitigs.push(varmer.clone());
                }
            }
        }
        if into.len() == 1 && out.len() == 0{
            let test_cov = unitigs[into[0].as_ref()].coverage;
            // only goes 1 k-mer forward
            if varmer.last().unwrap().0 < into[0].last().unwrap().0 + k as u32{
                if binomial_test(test_cov, cov, 0.10) > 0.005{
                    bad_unitigs.push(varmer.clone());
                }
            }
        }
    }

    for bad_unitig in bad_unitigs.iter(){
        log::trace!("REMOVED TIP");
        print_varmer(bad_unitig, true);
    }

    return bad_unitigs;
}
   
fn strand_bias_filter(dbg_frags: &mut Vec<FragDBG>, options: &Options, num_snps: usize, snp_pos_to_gn: &Vec<usize>) -> Vec<usize>{
    let mut pvalues = vec![];
    let mut snps_to_4_table: Vec<[u32;4]> = vec![[0; 4]; num_snps];
    for frag in dbg_frags.iter() {
        let ind;
        if frag.forward_strand{
            ind = 0
        }
        else{
            ind = 2
        }
        for (snp_pos,geno) in frag.seq.iter() {
            let mut geno = *geno as usize;
            if geno > 1{
                geno = 1;
            }
            snps_to_4_table[(*snp_pos - 1) as usize][geno + ind] += 1; 
        }
    }
    for (snp, table) in snps_to_4_table.iter().enumerate(){
        let p = fishers_exact(table).unwrap().two_tail_pvalue;
        pvalues.push((p, snp));
    }

    pvalues.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    let mut good_snps = FxHashSet::default();

    for i in 0..pvalues.len(){
        if pvalues[i].0 < options.strand_bias_fdr * (i as f64) / (pvalues.len() as f64){
            for j in i..pvalues.len(){
                log::trace!("THRESHOLD STRAND BIAS SNP {} : {}", pvalues[j].1 + 1, pvalues[j].0);
                log::trace!("TABLE {:?}", snps_to_4_table[pvalues[j].1]);
                let table = snps_to_4_table[pvalues[j].1];
                let ratio;
                if table[1] == 0 || table[2] == 0{
                    ratio = f64::MAX;
                }
                else{
                    ratio = (table[0] * table[3]) as f64 / (table[1] * table[2]) as f64;
                }
                let ratio = ratio.max(1./ratio);
                log::trace!("OR: {}", ratio);

                //Require high odds ratio filter for high coverage datasets
                if ratio < 1.5{
                    good_snps.insert((pvalues[j].1 + 1) as u32);
                }
            }
            break;
        }
        else{
            good_snps.insert((pvalues[i].1 + 1) as u32);
        }
    }

    let old_pos_to_new_pos_map = (1..=num_snps).filter(|x| good_snps.contains(&(*x as u32))).enumerate().map(|(i, x)| (x as u32, i as u32 + 1)).collect::<FxHashMap<u32, u32>>();
    log::trace!("GOOD SNPS: {:?}", good_snps);

    subsample_positions_fragdbg(dbg_frags, &good_snps, &old_pos_to_new_pos_map, snp_pos_to_gn)
}
