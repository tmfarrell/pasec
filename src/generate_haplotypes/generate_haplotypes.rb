#!/usr/bin/env ruby

require 'pp'
require 'ostruct'
require 'optparse'
require_relative 'utilities.rb'

# cmd-line parser 
class AmpliconSeqAnalysisParser
  def self.parse(args) 
    opt_parser = OptionParser.new do |opts| 
      opts.banner = "\nUsage: amplicon_seq_analysis.rb haplotypes [options]\n\n" +
                    "Computes haplotype coverage or allele frequencies over a given genomic region." 

      opts.separator ""
      opts.separator "Common args:" 
      opts.on_tail("-h", "--help", "Show this message") do
        puts opts
        exit
      end
    end
    opt_parser.parse!(args)
  end
end

class HaplotypesParser
  def self.parse(args)
    options = OpenStruct.new
    options.verbose = false
    options.print_to = 'stdout'
    options.print_reads = false
    options.mask_intervals_file = nil
    options.mask_base = "T"
    options.deletion_char = "-"
    options.min_haplotype_freq = 0.01
    options.min_haplotype_coverage = 2
    options.max_ref_mismatch = 100
    options.bwa_aligner_used = 'aln'
    options.min_mapping_quality = 0
    options.haplotype_clustering_edit_dist = 1
    options.run_haplotype_cluster_contraction = true
    options.haplotype_clustering_coverage_ratio = 0.15
    options.max_N_fraction = 2.0 / 300.0 # 2 N per 300 bases

    opt_parser = OptionParser.new do |opts|
      opts.banner = "\nUsage: amplicon_seq_analysis.rb haplotypes [options]\n\n" +
                    "Computes haplotype coverage over a given genomic region (specified by bed_file and bed_file_index) " +
                    "for a given alignment (bam_file)." 

      opts.separator ""
      opts.separator "Required args:" #required
      opts.on("--bam BAM_FILE", "Path to bam file.") do |bam_file|
        options.bam_file = bam_file
        if bam_file.include? "sorted"
          options.sample_id = File.basename(bam_file, ".sorted.bam")
        else
          options.sample_id = File.basename(bam_file, ".bam")
        end
      end
      opts.on("--bed BED_FILE", 
              "Path to bed file, containing genomic coordinates like ('chrom','start','end').") do |bed_file|
        options.bed_file = bed_file
      end
      opts.on("--bed_index BED_FILE_INDEX", 
              "Indices (1-based) of loci in bed file to compute haplotypes for. Either single index like '2', "+\
              "comma-separated list like '1,3,4' or range like '1-3'.") do |bed_file_index|
        options.bed_file_index = bed_file_index
      end

      opts.separator ""
      opts.separator "Optional args:" 
      opts.on("--min_haplotype_coverage MIN_HAPLOTYPE_COVERAGE", 
              Integer, "Filter haplotypes with coverage below this threshold. Default=0.") do |coverage| 
        options.min_haplotype_coverage = coverage
      end
      opts.on("--min_haplotype_freq MIN_HAPLOTYPE_FREQ", 
              Float, "Filter haplotypes with frequency below this threshold. Default=0.01.") do |freq| 
        options.min_haplotype_freq = freq
      end
      opts.on("--min_mapping_quality MIN_MAPPING_QUAL", 
              Integer, "Filter haplotypes with mapping quality below this threshold. Default=0.") do |qual| 
        options.min_mapping_quality = qual
      end
      opts.on("--max_ref_mismatch MAX_MISMATCH", 
              Integer, "Filter haplotypes with # reference mismatches over this threshold."+\
              " Default=20.") do |mismatch| 
        options.max_ref_mismatch = mismatch
      end
      opts.on("--max_N_fraction MAX_N_FRACTION", 
              Float, "Filter haplotypes with proportion of 'N' below this threshold. Must be 0.0 <= x < 1."+\
              " Default=0.006 (or 2 in 300 bases).") do |frac| 
        options.max_N_fraction = frac
      end
      opts.on("--bwa_aligner_used BWA_ALIGNER", 
              "Either 'aln' or 'mem'. Determines how SAM fields are parsed, since these differ slightly for each aligner. Default='aln'.") do |aligner|
        if ! ['aln','mem'].include?(aligner)
          puts "Invalid input for flag '--bwa_aligner_used': " + aligner
          exit
        else
          options.bwa_aligner_used = aligner
        end
      end
      opts.on("--haplotype_clustering-off", 
              "Turn haplotype cluster contraction off. If this is not passed, haplotype clustering contraction will be run with the below options.") do |run|
        options.run_haplotype_cluster_contraction = false
      end 
      opts.on("--haplotype_clustering-coverage_ratio MINOR_TO_MAJOR_COVERAGE_RATIO", 
                Float, "Threshold ratio between minor and major haplotypes to be grouped as one (i.e. clustered).",
                "Minor haplotypes less than or equal this threshold, which also satisfy --haplotype_clustering_edit_dist",
                "with a major haplotype, will be filtered and their counts included as major haplotype. Default=0.125 (or 1 to 8).") do |ratio| 
        options.haplotype_clustering_coverage_ratio = ratio
      end
      opts.on("--haplotype_clustering-edit_dist MISMATCH", 
                Integer, "Threshold edit distance between minor and major haplotypes to be grouped as one (i.e. clustered).",
                "Minor haplotypes less than or equal to this threshold, which also satisfy --haplotype_clustering_coverage_ratio",
                "with a major haplotype, will be filtered and their counts included as major haplotype. Default=1.") do |mismatch| 
        options.haplotype_clustering_edit_dist = mismatch
      end
      opts.on("--mask_bed_file MASK_FILE", 
              "List of coordinates to mask for each loci. Should be tab-delimited with columns 'chrom','start','end'.") do |mask_file| 
        options.mask_intervals_file = mask_file
      end
      opts.on("--mask_base MASK_BASE", 
              "Character to replace masked positions with. Default='T'.") do |mask_base| 
        options.mask_base = mask_base
      end
      opts.on("--deletion_char DELETION_CHAR", 
              "Character to fill in/ represent deletions. Default='-'.") do |deletion_char| 
        options.deletion_char = deletion_char
      end
      opts.on("--print_reads", 
              "Print reads along w/ haplotypes. If --print_to 'stdout', will be printed before haplotypes. If --print_to $dir, will be",
              "printed to its own file $dir/$sample_id.$amplicon.reads.tsv.") do |print_reads| 
        options.print_reads = true
      end

      opts.separator ""
      opts.separator "Common args:"
      opts.on("--print_to OUT_STREAM", "Out stream to print output to. Either 'stdout' or a valid directory, where if",
              "valid directory path will output files {sample_id}.read.tsv and {sample_id}.haplotypes.tsv. Default='stdout'") do |out|
        options.print_to = out
      end
      opts.on_tail("-h", "--help", "Show this message") do
        puts opts
        exit
      end
      opts.on("-v", "--verbose", "Run verbosely.") do |v| 
        options.verbose = v
      end 
    end # opt_parser
    opt_parser.parse!(args)
    options
  end # parse
end

subcommands = {
  'haplotypes' => HaplotypesParser
}

# main 
# parse opts
command = ARGV.shift
if ["--help", "-h"].include?(command) 
  AmpliconSeqAnalysisParser.parse(args)
else
  options = subcommands[command].parse(ARGV)
end

if ['haplotypes'].include?(command)
  # get loci 
  loci = parse_bed_file(options.bed_file, options.bed_file_index)
  puts "loci: #{loci.to_s}" if options.verbose
  # get mask intervals, which gives hash like {chrom => [(start0, end0),(start1,end0),...]}
  # for each chrom in loci 
  if options.mask_intervals_file
    mask_intervals = parse_mask_intervals_file(options.mask_intervals_file, loci)
    puts "mask_intervals: #{mask_intervals.to_s}" if options.verbose
  else
    mask_intervals = Hash[loci.map{|chrom,_,_,_| [chrom, []]}]
  end 
  loci.each{ |locus| 
    # get read set, filtering w/ cmd-line opts
    reads, haplotype_set = get_haplotype_set(options.bam_file, locus, options.min_haplotype_coverage, options.min_haplotype_freq,
                                             options.min_mapping_quality, options.max_ref_mismatch, options.max_N_fraction, 
                                             options.verbose, options.bwa_aligner_used, mask_intervals=mask_intervals[locus.first], 
                                             trim_to='locus', deletion_char=options.deletion_char, mask_base=options.mask_base)
    # if no reads, go to next locus
    haplotype_set_filtered = haplotype_set.select{|h| h[:passed_filters]}
    if haplotype_set_filtered == [] or haplotype_set.length == 0
      next
    end
    # sort haplotype set 
    haplotype_set.sort!{|x,y| x[:haplotype_idx] <=> y[:haplotype_idx]}
    # run clustering contraction, if enabled 
    if options.run_haplotype_cluster_contraction
      haplotype_set_clustered, clusters = apply_haplotype_cluster_contraction(haplotype_set_filtered, coverage_ratio_threshold=options.haplotype_clustering_coverage_ratio, 
                                                                              mismatch_threshold=options.haplotype_clustering_edit_dist,
                                                                              sum_cluster_coverage=true, collect_subclusters=false, verbose=options.verbose)
      puts "Those then clustered into #{haplotype_set.length} haplotypes, using distance #{options.haplotype_clustering_edit_dist}" +\
      " and minor to major coverage ratio of #{options.haplotype_clustering_coverage_ratio}." if options.verbose
      # add haplotype_cluster field to haplotype_set
      # first add it to those haplotypes that were clustered
      clusters.each{ |cluster, haplotype_indices|
                     ([[cluster, 0]] + haplotype_indices).each{ |haplotype_idx, _|  haplotype_set[haplotype_idx][:haplotype_cluster] = cluster }
      }
    end
    # then add it to those that were not clustered
    haplotype_set = haplotype_set.each{|h| (not h.keys.include?(:haplotype_cluster)) ? h.update({:haplotype_cluster => h[:haplotype_idx]}) : h }
    if options.print_reads
      print_reads(reads, options.sample_id, locus, to=options.print_to)
    end
    print_haplotypes(haplotype_set, options.sample_id, locus, to=options.print_to, include_read_idx=options.print_reads)
  }
end
