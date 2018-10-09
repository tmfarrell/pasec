#!/usr/bin/env ruby

require 'fileutils'
require_relative './cigar_utils.rb'
require_relative './haplotype_clustering_utils.rb'

## consts
BASES = ["A","T","C","G","-","N","*"]

## utilities
def bool_to_int(bool)
  return(bool ? 1 : 0)
end

def to_freq(count, total) 
  return (count.to_f / total.to_f).round(4)
end

def safe_div(num, denom)
  r = (num.to_f / denom.to_f)
  if r.nan?
    return 0
  else
    return r.round(4)
  end
end 

def mean(arr) 
  return((arr.reduce(:+).to_f / arr.length.to_f).round(4))
end 

def num_mismatches(mismatch_str)
  return BASES.map{|base| mismatch_str.count(base)}.reduce(:+) - mismatch_str.count('^')
end

def metrics_to_str(metrics) 
  return(["MM","D","I","E"].map{|op| 
         "#{op}:#{metrics[op]}"}.join('/'))
end 

def to_metrics(haplotype)
  ## give mismatch, deletion and insertion metrics for a given haplotype
  return({"MM" => mean(haplotype[:mismatches]), 
          "D" => haplotype[:haplotype].count('-'),
          "I" => mean(haplotype[:cigars].map{|cigar| cigar.count('I')}),
          "E" => mean(haplotype[:edit_dists])})
end 

def locus_to_str(locus, trimmed_to) 
  chrom, coords, _ = locus
  if trimmed_to == 'locus'
    return "#{chrom}:#{coords.first}-#{coords.last}"
  elsif trimmed_to == 'locus_start'
    return "#{chrom}:#{coords.first}"
  else
    return "#{chrom}"
  end
end 

def allele_counts_to_str(allele_counts, to_freqs=false)
  if not to_freqs
    return(BASES.map{|base| "#{base}:#{allele_counts[base]}"}.join('/'))
  else
    total = allele_freqs.values.reduce(:+).to_f 
    return(BASES.map{|base| "#{base}:#{to_freq(allele_freqs[base], total)}"}.join('/'))
  end
end 

## i/o functions
def print_run_stats(bam_file, loci, solexa)
  chroms = loci.map{|chrom, coords, id, len| chrom} + ["*"]
  run_stats = (`samtools idxstats #{bam_file}`).split("\n")
                                               .map{|line| line.split("\t")}
                                               .select{|chrom,_,mapped,unmapped| chroms.include?(chrom)} 
  #puts "#{run_stats}"
  total_aligned = run_stats.map{|chrom,_,mapped,unmapped| mapped.to_i + unmapped.to_i}.reduce(:+)
  run_stats_str = (solexa + "\t" + \
                   run_stats.select{|chrom,_,_,_| chrom == "*"}.map{|_,_,_,unmapped| "total/#{safe_div(unmapped.to_f, total_aligned.to_f)}/#{total_aligned}"}.join() + "\t" +\
                   run_stats.select{|chrom,_,_,_| chrom != "*"}
                            .map{|chrom,_,mapped,unmapped| "#{chrom}/#{safe_div(mapped.to_f, mapped.to_f + unmapped.to_f)}/#{(mapped.to_i + unmapped.to_i)}"}.join("\t"))
  puts run_stats_str
end

def parse_reads_file(read_file)
  reads = File.open(read_file).read.split("\n")[1..-1]
              .map{|line| line.split("\t")}
              .map{|_,read_index,coverage,_,_,_,_,_,read,cigar,mismatch_str,_,_|
                   {:read_index => read_index, :coverage => coverage.to_i, :read => read, :cigar => cigar, :mismatch_str => mismatch_str}}
  return(reads)
end 

def print_reads(read_set, solexa, locus, to='stdout', format='tsv', file_tag='')
  # prints read set to stdout
  chrom, coords, amplicon = locus
  if ! ['fasta','tsv'].include?(format)
    raise(Exception, "Invalid output format. Must be 'fasta' or 'tsv'.")
  end
  if to != 'stdout' and ! File.directory?(to)
    raise(Exception, "Invalid out stream. Must be either 'stdout' or valid directory.")
  elsif to != 'stdout'
    if file_tag.empty?
      #FileUtils.touch(to + "/#{solexa}.#{amplicon}.reads.#{format}")
      out = open(to + "/#{solexa}.#{amplicon}.reads.#{format}", 'w')
    else
      #FileUtils.touch(to + "/#{solexa}.#{locus_str}.#{file_tag}.reads.#{format}")
      out = open(to + "/#{solexa}.#{locus_str}.#{file_tag}.reads.#{format}", 'w')
    end
    out.write("")
  end
  if read_set.empty? 
    FileUtils.touch(to + "/#{solexa}.#{amplicon}.reads.#{format}")
    return
  end
  read_set.sort!{|x,y| x[:haplotype_idx] <=> y[:haplotype_idx]}
  if format == 'tsv'
    line = "Solexa_ID\tlocus\tamplicon\tread_idx\thaplotype_idx\tedit_dist\tlen\tcigar\tmismatch\tqual\tread\n"
    (to == 'stdout') ? (puts line) : (out.write(line))
    read_set.each{ |read| 
      if read[:mismatch_str] != nil
        line = "#{solexa}\t#{read[:chrom]}:#{read[:pos]}\t#{amplicon}\t#{read[:read_idx]}\t#{read[:haplotype_idx]}\t#{read[:edit_dist]}\t#{read[:seq].length}\t" +
          "#{to_cigar_str(read[:cigar])}\t#{read[:mismatch_str]}\t#{read[:qual_str]}\t#{read[:seq]}\n"
      else
        line = "#{solexa}\t#{read[:chrom]}\t#{read[:pos]}\t#{amplicon}\t#{read[:read_idx]}\t#{to_cigar_str(read[:cigar])}\t#{read[:seq].length}\t#{read[:seq]}\t#{read[:qual_str]}\n"
      end
      (to == 'stdout') ? (puts line) : (out.write(line))
    }
  end
end 
 
def print_haplotypes(haplotype_set, solexa, locus, to='stdout', include_read_idx=false, 
                     format='tsv', with_extras=false, file_tag='')
=begin
    Prints final haplotype set to either STDOUT or the path provided by to, in either 'fasta' or 'tsv' format. 
    input: 
      haplotype_set :: [{:haplotype => String, :freq => Integer}]
      locus :: [[chrom :: String, [start :: Integer, end :: Integer]]
=end
  chrom, coords, amplicon = locus
  if ! ['fasta','tsv'].include?(format)
    raise(Exception, "Invalid output format. Must be 'fasta' or 'tsv'.")
  end
  if to != 'stdout' and ! File.directory?(to)
    raise(Exception, "Invalid out stream. Must be either 'stdout' or valid directory.")
  elsif to != 'stdout'
    if file_tag.empty?
      out = open(to + "/#{solexa}.#{amplicon}.haplotypes.#{format}", 'w')
    else
      out = open(to + "/#{solexa}.#{amplicon}.#{file_tag}.haplotypes.#{format}", 'w')
    end
    out.write("")
  end
  # sort by coverage
  haplotype_set.sort!{|x, y| y[:coverage] <=> x[:coverage]}
  if format == 'tsv'
    if include_read_idx != false
      line = "Solexa_ID\tlocus\tamplicon\tcoverage\tfreq\tlen\tread_idx\thaplotype_cluster_idx\thaplotype_idx\tpass_coverage\tpass_freq\tpass_len\thaplotype\n"
    else
      line = "Solexa_ID\tlocus\tamplicon\tcoverage\tfreq\tlen\thaplotype_cluster_idx\thaplotype_idx\tpass_coverage\tpass_freq\tpass_len\thaplotype\n"
    end
    (to == 'stdout') ? (puts line) : (out.write(line))
  end
  if haplotype_set.empty?
    FileUtils.touch(to + "/#{solexa}.#{amplicon}.haplotypes.#{format}")
    return
  end
  haplotype_set.each{ |haplotype|
    #haplotype_locus_str = locus_to_str([haplotype[:chrom], [haplotype[:start], haplotype[:end]]])
    if format == 'fasta'
      line = ">#{solexa}__#{haplotype[:haplotype_locus]}__#{haplotype[:coverage]}__#{haplotype[:final_freq]}"+"\n"+haplotype[:haplotype]+'\n'
    elsif format == 'tsv'
      if ! with_extras
        if include_read_idx
          line = "#{solexa}\t#{haplotype[:haplotype_locus]}\t#{amplicon}\t#{haplotype[:coverage]}\t#{haplotype[:final_freq].round(8)}\t#{haplotype[:haplotype_len]}\t" + 
            "#{haplotype[:read_idx]}\t#{haplotype[:haplotype_cluster]}\t#{haplotype[:haplotype_idx]}\t#{bool_to_int(haplotype[:passed_coverage_filter])}\t" + 
            "#{bool_to_int(haplotype[:passed_freq_filter])}\t#{bool_to_int(haplotype[:passed_haplotype_len_filter])}\t#{haplotype[:haplotype]}\n"
        else 
          line = "#{solexa}\t#{haplotype[:haplotype_locus]}\t#{amplicon}\t#{haplotype[:coverage]}\t#{haplotype[:final_freq].round(8)}\t#{haplotype[:haplotype_len]}\t" + 
            "#{haplotype[:haplotype_cluster]}\t#{haplotype[:haplotype_idx]}\t#{bool_to_int(haplotype[:passed_coverage_filter])}\t" + 
            "#{bool_to_int(haplotype[:passed_freq_filter])}\t#{bool_to_int(haplotype[:passed_haplotype_len_filter])}\t#{haplotype[:haplotype]}\n"
        end
      else
        line = "#{solexa}\t#{haplotype_locus_str}\t#{amplicon}\t#{haplotype[:coverage]}\t#{to_freq(haplotype[:coverage],total)}\t"+\
               "#{haplotype[:cigars]}\t#{haplotype[:edit_dists]}\t#{haplotype[:mismatches]}\t#{haplotype[:mismatch_strs]}\t#{haplotype[:haplotype].length}\t#{haplotype[:haplotype]}\n"
      end
    end
    (to == 'stdout') ? (puts line) : (out.write(line))
  }
  out.close() if out
end

def print_alleles(allele_counts, locus, solexa, to='stdout')
  locus_str = locus_to_str(locus)
  if to != 'stdout' and ! File.directory?(to)
    raise(Exception, "Invalid out stream. Must be either 'stdout' or valid directory.")
  elsif to != 'stdout'
    out = open(to + "/#{solexa}.#{locus_str}.allele_counts.tsv", 'w')
    out.write("solexa\t#{solexa}\n")
  else
    puts "solexa\t#{solexa}\n"
  end
  chrom, coords, id = locus
  start = coords.first.to_i
  for i in 0...allele_counts.length
    pos = start + i
    counts = allele_counts[i]
    # make order explicit: A, T, C, G 
    line = "#{chrom}:#{pos}\t#{allele_counts_to_str(counts)}\n"
    (to == 'stdout') ? (puts line) : (out.write(line))
  end
  out.close() if out
end 

def parse_bed_file(filepath, idx=nil)
=begin
    Parses BED file.
    input: 
      filepath :: String : path to bed file
      idx :: nil | String : limit returned genomic coordinates to some index of the BED file (1-based), 
                            either some integer like '1', an interval of indices like '2-4' or a string of comma-separated indices like '2,4,5'
    output: 
      genomic_coords :: [[chrom :: String, [start :: Integer, end :: Integer]]]
=end
  # open bed, split lines and fields
  genomic_coords = File.open(filepath).read.split("\n")
                                           .map{|line| line.split("\t")}
                                           .map{|chrom,start,end_,id,len| [chrom, [start.to_i, end_.to_i], id, len]}
  # return coords in index, if given
  if !idx
      return genomic_coords
  elsif idx.length == 1
    return [genomic_coords[idx.to_i - 1]]
  elsif idx.include?('-')
    start_idx, end_idx = idx.split('-')
    return genomic_coords.values_at((start_idx.to_i - 1)..(end_idx.to_i - 1))
  elsif idx.include?(',')
    idx = idx.split(',').map{|i| i.to_i - 1}
    return genomic_coords.select{|gc| idx.include?(genomic_coords.index(gc))}
  end
end

def parse_mask_intervals_file(mask_file, loci)
  mask_intervals = {} 
  genomic_coords = parse_bed_file(mask_file) 
  loci.each{|chrom,_,_|
    mask_intervals[chrom] = genomic_coords.select{|chr,coords,_| chr == chrom}
                                          .map{|chr,coords,_| [coords.first, coords.last]}
  }
  return(mask_intervals)
end

def read_haplotype_file(path)
  haplotype_file_lines = File.open(path, 'r').read.split("\n").map{|line| line.split("\t")} 
  haplotype_set = haplotype_file_lines.map{|solexa,locus,coverage,freq,haplotype,metrics|
    {:solexa => solexa, :locus => locus, :coverage => coverage.to_i, :freq => freq.to_f, :haplotype => haplotype, :metrics => metrics}}
  return haplotype_set
end 

## main functions 
def get_haplotype_set(bam_file, locus, min_read_coverage, min_read_freq, min_mapping_quality, max_ref_mismatch, 
                      max_N_fraction, verbose=false, bwa_aligner_used=nil, mask_intervals=[], trim_to='locus')
=begin
   Gets unique set of reads (with metrics) from bam_file within a some locus, with some filtering. 

   input:
     bam_file :: String : filepath of input bam file 
     locus :: [chromosome :: String, [start :: Integer, end :: Integer]] : genomic region to fetch aligned reads over  
     min_read_coverage :: Integer : minimum threshold for read coverage 
     min_mapping_quality :: Integer : minimum threshold for mapping quality
     max_N_fraction :: Float (0.0 <= x < 1) : maximum threshold for fraction of read containing 'N'
   
   output: 
     read_set :: [{:chrom :: String, :pos => Integer, :cigar :: String, :seq :: String, :len :: Integer, :coverage :: Integer}]
=end
  chromosome, coords, id, min_haplotype_len = locus
  # get raw reads, filtering reads that do not overlap locus, have low quality or are reverse complements or secondary alignments (bit flags 4 and 256) 
  puts "" if verbose
  puts "Getting reads..." if verbose  
  reads = (`samtools view -F 4 -F 256 -q #{min_mapping_quality} #{bam_file} #{chromosome}:#{coords.first}-#{coords.last}`)
            .split("\n")
            .map{|line| line.split("\t")}
  # select relevant fields and put into list of hashes 
  if bwa_aligner_used == 'mem'
    reads = reads.map{|_,_,chrom,start,_,cigar,_,_,_,seq,qual_str,edit_dist,mismatch_str,_,_| 
      {:chrom => chrom, :pos => start.to_i, :cigar => to_cigar_array(cigar), :seq => seq, 
       :mismatches => num_mismatches(mismatch_str[5..-1]), :mismatch_str => mismatch_str[5..-1], 
       :len => seq.length, :edit_dist => edit_dist[5..-1].to_i, :qual_str => qual_str}}.select{|read| read[:mismatches] < max_ref_mismatch}
  elsif bwa_aligner_used == 'aln'
    reads = reads.map{|_,_,chrom,start,_,cigar,_,_,_,seq,qual_str,_,edit_dist,_,_,mismatches,_,_,mismatch_str| 
      {:chrom => chrom, :pos => start.to_i, :cigar => to_cigar_array(cigar), :seq => seq, 
       :mismatches => mismatches[5..-1].to_i, :mismatch_str => mismatch_str[5..-1], 
       :len => seq.length, :edit_dist => edit_dist[5..-1].to_i, :qual_str => qual_str}}.select{|read| read[:mismatches] < max_ref_mismatch}
  else
    reads = reads.map{|_,_,chrom,start,_,cigar,_,_,_,seq,qual_str| 
                      {:chrom => chrom, :pos => start.to_i, :cigar => to_cigar_array(cigar), :seq => seq, :qual_str => qual_str}}
  end
  puts "#{reads.length} reads with overlapping #{locus}, with minimum quality #{min_mapping_quality} and less than #{max_ref_mismatch} mismatches w/ reference." if verbose
  # get read set, with coverage metric
  puts "" if verbose
  puts "Compiling haplotypes..." if verbose
  reads = reads.map do |read|
    aligned_seq, aligned_start = get_aligned_seq(read[:seq], read[:cigar], locus, read[:pos], mask_intervals, trim_to=trim_to, verbose=verbose) 
    read.update({:haplotype => aligned_seq, :haplotype_locus => locus_to_str([read[:chrom], [aligned_start, aligned_start + aligned_seq.length], id, ''], trim_to)})
  end
  unique_read_seqs = reads.map{|read| read[:haplotype]}.uniq
  puts "Of those, #{unique_read_seqs.length} had a unique haplotype (after trimming and masking)." if verbose
  puts "" if verbose
  puts "Counting coverage..." if verbose
  # group reads by haplotype (unique read seq) index
  reads = reads.each_with_index.map{|read, i| read.update({:haplotype_idx => unique_read_seqs.index(read[:haplotype]), :read_idx => i})}
  # gives back a map like {unique_read_idx_0 => [read_idx_j,...], ..., unique_read_idx_M => [read_idx_j',...]} 
  read_seq_groups = (0...reads.length).group_by{|i| reads[i][:haplotype_idx]}
  # convert to array of 3-tuples, like (unique_read_idx_0, [read_idx_j,...], coverage) where coverage is # reads in read index
  read_seq_groups = read_seq_groups.to_a.map{|h_idx, read_seq_idx| [h_idx, read_seq_idx, read_seq_idx.length]}
  # build haplotype_set from read_seq_groups 
  haplotype_set = []
  read_seq_groups.each{|h_idx, read_seq_idx, coverage|
    haplotype_set << {:haplotype => unique_read_seqs[h_idx], :coverage => coverage, :passed_coverage_filter => (coverage >= min_read_coverage), 
                      :haplotype_idx => h_idx, :haplotype_locus => reads[read_seq_idx[0]][:haplotype_locus], :haplotype_len => unique_read_seqs[h_idx].length,
                      :read_idx => read_seq_idx[0]}
  }
  # remove haplotype fields from reads
  reads = reads.map{|read| read.reject{|k,v| [:haplotype,:haplotype_locus].include?(k)}}
  # select groups that pass coverage filter
  read_seq_groups_filtered = read_seq_groups.select{|_, _, coverage| coverage >= min_read_coverage}
  # then group filtered reads by start pos, sorting by most popular start pos
  read_pos_groups = read_seq_groups_filtered.map{|_, read_seq_idx, _| read_seq_idx}.flatten
                                            .group_by{|i| reads[i][:pos]}.to_a.sort{|x, y| y[1].length <=> x[1].length}
  # add best_read field to haplotype set, indicating which read is the best in that haplotype group 
  # for each haplotype w/ minimum coverage
  read_seq_groups_filtered.each{|h_idx, read_seq_idx, _|
    # loop through from most to least popular position
    for i in 0...read_pos_groups.length
      pos, read_pos_idx = read_pos_groups[i]
      # if reads w/ this haplotype are in this position read index
      x = (read_seq_idx & read_pos_idx)
      if x.length > 0
        # take the first one as the best read
        best_read = x[0]
        break
      end
    end
    # set that to haplotype's read_idx 
    haplotype_set[h_idx] = haplotype_set[h_idx].update({:read_idx => best_read})
  }
  puts "reads:" if verbose
  puts reads if verbose
  puts "before filtering:" if verbose
  puts haplotype_set if verbose
  passed_coverage_filter_idx = haplotype_set.select{|read| read[:passed_coverage_filter]}.map{|read| read[:haplotype_idx]}
  total_coverage = haplotype_set.select{|read| passed_coverage_filter_idx.include?(read[:haplotype_idx])}
                                .map{|read| read[:coverage]}.reduce(:+).to_f
  puts "Of those, #{passed_coverage_filter_idx.length} were >= the minimum haplotype coverage #{min_read_coverage}, with total coverage #{total_coverage}." if verbose
  # then filter reads w/ haplotype shorter than certain length
  haplotype_set = haplotype_set.map{|read| read.update({:passed_haplotype_len_filter => (read[:haplotype].length >= min_haplotype_len.to_i)})}
  passed_len_filter_idx = haplotype_set.select{|read| read[:passed_haplotype_len_filter]}
                                       .map{|read| read[:haplotype_idx]}
  passed_filters_idx = (passed_coverage_filter_idx & passed_len_filter_idx)
  puts "Of those, #{passed_filters_idx.length} were >= the minimum haplotype length #{min_haplotype_len}." if verbose
  # then filter reads w/ low frequency
  total_coverage = haplotype_set.select{|read| passed_filters_idx.include?(read[:haplotype_idx])}.map{|read| read[:coverage]}.reduce(:+).to_f
  # if passed filters, add freq; else set freq to 0.0
  haplotype_set = haplotype_set.map{|r| (passed_filters_idx.include?(r[:haplotype_idx])) ? 
                                        r.update({:freq => (r[:coverage].to_f / total_coverage)}) : 
                                        r.update({:freq => 0.0})} 
  haplotype_set.sort!{|r_i, r_j| r_j[:coverage] <=> r_i[:coverage]}
  #puts "unfiltered reads:" if verbose
  #print_reads(haplotype_set, "-----", locus) if verbose
  haplotype_set = haplotype_set.map{|r| r.update({:passed_freq_filter => (r[:freq] >= min_read_freq)})} 
  # update index of haplotypes that passed filters, then recompute freq
  passed_filters_idx = (passed_filters_idx & (haplotype_set.select{|read| read[:passed_freq_filter]}.map{|read| read[:haplotype_idx]}))
  total_coverage = haplotype_set.select{|read| passed_filters_idx.include?(read[:haplotype_idx])}.map{|read| read[:coverage]}.reduce(:+).to_f
  haplotype_set = haplotype_set.map{|r| (passed_filters_idx.include?(r[:haplotype_idx])) ? 
                                        r.update({:final_freq => (r[:coverage].to_f / total_coverage)}) :
                                        r.update({:final_freq => 0.0})} # add freq
  puts "Of those, #{passed_filters_idx.length} were above the minimum read freq #{min_read_freq}." if verbose
  haplotype_set = haplotype_set.map{|r| r.update({:passed_filters => [r[:passed_coverage_filter],r[:passed_freq_filter],r[:passed_haplotype_len_filter]].all?})}
  # and w/ high N fraction  
  #haplotype_set.reject!{|r| (Float(r[:seq].count('N')) / Float(r[:seq].length)) >= max_N_fraction}   # in-place
  #puts "Of those, #{haplotype_set.length} were below the max N fraction %0.5f." % max_N_fraction if verbose
  if haplotype_set == nil
    return nil
  end 
  # add padded version of each read
  #insertion_positions = get_insertion_position_set(haplotype_set)
  #haplotype_set = haplotype_set.map{|read| read.update({:padded_cigar => to_padded_cigar_arr(read[:cigar], insertion_positions)})}
  #haplotype_set = haplotype_set.map{|read| read.update({:padded_seq => get_aligned_seq(read[:seq], read[:padded_cigar], locus, read[:pos], trim_to=nil)})}
  # sort reads by coverage
  #haplotype_set.sort!{|r_i, r_j| r_j[:coverage] <=> r_i[:coverage]}
  puts "top 10 unclustered haplotypes:" if verbose
  print_haplotypes(haplotype_set[0...10], locus, '-----') if (verbose and not haplotype_set.empty?)
  #if verbose
  #  haplotype_set.each{|read| puts read} 
  #end
  return [reads, haplotype_set]
end 

def get_allele_counts(haplotype_set, count_coverage=true, verbose=false)
  allele_counts = []
  max_haplotype_len = haplotype_set.map{|haplotype| haplotype[:haplotype].length}.max
  puts "max haplotype length: #{max_haplotype_len}" if verbose
  allele_counts_map = {"A" => 0, "T" => 0, "C" => 0, "G" => 0, "-" => 0, "N" => 0, "*" => 0}
  if max_haplotype_len
    for i in 0...max_haplotype_len
      allele_counts.push(allele_counts_map.dup)
    end
    for i in 0...haplotype_set.length
      coverage = haplotype_set[i][:coverage]
      haplotype = haplotype_set[i][:haplotype]
      for j in 0...haplotype.length
        base = haplotype[j]
        puts "base: #{base}" if verbose
        puts "allele counts at #{j}: #{allele_counts[j]}" if verbose
        if count_coverage
          allele_counts[j][base] = allele_counts[j][base] + coverage
        else
          allele_counts[j][base] = allele_counts[j][base] + 1
        end
      end
    end
  end
  return allele_counts
end 
