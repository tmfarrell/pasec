#!/usr/bin/env ruby

def to_cigar_str(cigar_arr)
  return cigar_arr.flatten.join()
end 

def to_cigar_array(cigar_str)
=begin
    Returns cigar array from cigar string.
    input:   cigar :: String like \*|([0-9]+[MIDNSHPX=])+
    output:  cigar_arr :: [[count :: Integer, cigar_op :: String]] 
=end
  return cigar_str.scan(/[0-9]+[MIDNSHPX=]/)
                  .map{|c| [c.match(/[0-9]+/)[0].to_i, c.match(/[MIDNSHPX=]/)[0]]}
end

def cigar_arr_to_long_cigar(cigar_arr) 
    return(cigar_arr.map{|count, op| op * count}.join)
end

def long_cigar_to_cigar_arr(long_cigar) 
    arr = []
    count = 0
    curr = long_cigar[0]
    for i in 0...long_cigar.length
      if curr == long_cigar[i] 
        count = count + 1
      else 
        arr = arr + [[count, curr]]
        curr = long_cigar[i]
        count = 1
      end
      if i == (long_cigar.length - 1)  # if last 
        arr = arr + [[count, curr]]
      end
    end     
  return(arr)
end

def get_insertion_positions(cigar_arr)
  positions =[]
  for i in 0...cigar_arr.length 
    count, op = cigar_arr[i]
    if op == "I" 
      curr_pos = cigar_arr[0...i].select{|count, op| op != "I"}.map{|count, op| count}.reduce(:+)
      positions = positions + (0...count).map{|i| curr_pos + i}
    end
  end
  return(positions)
end 

def get_insertion_position_set(read_set, verbose=false)
  insertion_positions = []
  read_set.each do |read|
    puts "for read with cigar #{read[:cigar]} and locus #{read[:chrom]}:#{read[:pos]}-#{read[:pos] + read[:seq].length}" if verbose
    insert_pos = get_insertion_positions(read[:cigar])
    puts "its insertion positions are: #{insert_pos}" if verbose
    insertion_positions = insertion_positions + insert_pos
  end 
  return(insertion_positions.flatten.uniq.sort)
end

def to_padded_cigar_arr(cigar_arr, insertion_positions, verbose=false) 
  long_cigar = cigar_arr_to_long_cigar(cigar_arr)
  #puts long_cigar if verbose
  return(long_cigar_to_cigar_arr(to_padded_long_cigar(long_cigar, insertion_positions, verbose)))
end

def to_padded_long_cigar(long_cigar, insertion_positions, verbose=false)    
  padded = ''
  insertion_pos_idx, unpadded_idx, padded_idx = 0, 0, 0
  while true do 
    # if curr_pos (in padded space) is an insertion position (in unpadded space)
    if padded_idx == insertion_positions[insertion_pos_idx] 
      # if base at curr_pos is not insert
      if long_cigar[unpadded_idx] != "I"
        # add P for padding 
        padded = padded + "P"
        insertion_pos_idx = insertion_pos_idx + 1
      else # if curr_pos is an insert
        # add to padded 
        padded = padded + "I"
        # incr unpadded index  
        unpadded_idx = unpadded_idx + 1
        # move to next insertion position
        insertion_pos_idx = insertion_pos_idx + 1
      end
    else
      padded = padded + long_cigar[unpadded_idx]
      padded_idx = padded_idx + 1
      unpadded_idx = unpadded_idx + 1
    end       
    if (unpadded_idx - padded.count("I")) >= long_cigar.length
      break
    end
    if insertion_pos_idx >= insertion_positions.length
      padded = padded + long_cigar[unpadded_idx..-1]
      break
    end
    #puts "padded: #{padded}, padded_idx: #{padded_idx}, unpadded_idx: #{unpadded_idx}, insertion_pos_idx: #{insertion_pos_idx}" if verbose 
  end
  return(padded)
end

def get_seq_intervals(seq, intervals, verbose=false) 
  new_seq = ''
  puts "#{intervals}" if verbose
  intervals.each do |interval|
    if interval.class == Array
      start, end_ = interval
      puts "#{start}, #{end_}" if verbose
      section = seq[start...end_]
      if section != nil
        new_seq = new_seq + section
      end
    elsif interval.class == String
      new_seq = new_seq + interval
    end    
  end
  return(new_seq)
end

def get_mask_interval_complement(mask_intervals, start_pos, end_pos, replace_w_homopolymer=false, 
                                 verbose=false, homopolymer_base="T")
  complement = []
  curr_pos = start_pos
  to_nonnegative = Proc.new {|x| (x < 0) ? 0 : x }
  puts "mask_intervals: #{mask_intervals.to_s}" if verbose
  puts "end_pos: #{end_pos.to_s}" if verbose
  mask_intervals.each{ |start, end_|
    complement = complement + [[curr_pos, start]]
    if replace_w_homopolymer
      complement = complement + [homopolymer_base * (end_ + 1 - start)]
    end
    curr_pos = end_ + 1
    if curr_pos > end_pos
      return(complement.map{|arr| (arr.class == String) ? arr : arr.map{|x| to_nonnegative.call(x)}})
    end
  }
  return((complement + [[curr_pos, end_pos]]).map{|arr| (arr.class == String) ? arr : arr.map{|x| to_nonnegative.call(x)}})
end

# return haplotype, haplotype_start_pos from read and cigar array 
def get_aligned_seq(seq, cigar_arr, interval, aligned_start_pos, mask_intervals=[], 
                    trim_to=nil, verbose=false, deletion_char="-")
=begin
    Retrieves portion of read seq which aligns to specified genomic interval, with
    optional masking and trimming.

    input: 
      seq :: String : read sequence
      cigar_arr ::  [[count :: Integer, cigar_op :: String]] 
      interval :: [chrom :: String, [start :: Integer, end :: Integer]]
      aligned_start_pos :: Int : start position of alignment.

    output: 
      aligned_seq :: String : portion of read sequence aligning to desired interval
      start_pos :: Int : start pos of haplotype, accounting for trimming/etc. 
=end
  if ! [nil, 'locus', 'locus_start'].include?(trim_to)
    raise(Exception, "trim_to parameter must be nil, 'locus' or 'locus_start'.")
  end
  aligned_seq = ''
  if verbose 
    puts ""
    puts 'Getting aligned sequence for:' 
    puts 'cigar: ' + cigar_arr.to_s
    puts "seq: #{seq.length} #{seq}"
  end 
  cigar_arr.each_with_index do |cigar_pair, i|
    count, cigar_op = cigar_pair 
    if i == 0 and cigar_op == 'S'
      aligned_start_pos = aligned_start_pos + count
    end
    puts "#{[count, cigar_op]}" if verbose
    case cigar_op
    when 'M', 'N'  # if match or N, add base as usual
      aligned_seq += seq[0...count]
      seq = seq[count...seq.length] 
    when 'P'       # if padding position (i.e. where insertions occur in any another read), add "*"
      aligned_seq += (['*'] * count).join
    when 'D'       # if deletion, insert "-" 
      aligned_seq += ([deletion_char] * count).join
    when 'H', 'S', 'I'
      seq = seq[count...seq.length]
    end
    puts "aligned_seq so far: #{aligned_seq.length} #{aligned_seq}" if verbose
  end
  if ! mask_intervals.empty?
    mask_intervals = mask_intervals.map{|start, end_| [start - aligned_start_pos, end_ - aligned_start_pos]}
  end
  if trim_to
    _, coords, _, min_haplotype_len = interval  
    puts "aligned seq before trimming: #{aligned_seq} #{aligned_seq.length}" if verbose
    #min_haplotype_len = coords.last - coords.first - 25
    if aligned_seq.length < min_haplotype_len.to_i
      puts "aligned seq too short: #{aligned_seq.length}" if verbose
      return [aligned_seq, aligned_start_pos]
    end
    if trim_to == 'locus'
      start_pos = coords.first - aligned_start_pos
      end_pos = aligned_seq.length - ((aligned_start_pos + aligned_seq.length) - coords.last) + 1
      mask_interval_complement = get_mask_interval_complement(mask_intervals, start_pos, end_pos, true, verbose)
      puts "mask intervals: #{mask_intervals.to_s}" if verbose
      puts "mask intervals complement: #{mask_interval_complement}" if verbose
      puts "start_pos, end_pos: #{start_pos.to_s}, #{end_pos.to_s}" if verbose
      aligned_seq = get_seq_intervals(aligned_seq, mask_interval_complement) 
    elsif trim_to == 'locus_start'
      start_pos = coords.first - aligned_start_pos
      #mask_intervals = [[0, (coords.first - aligned_start_pos - 1)]] + mask_intervals
      aligned_seq = get_seq_intervals(aligned_seq, get_mask_interval_complement(mask_intervals, start_pos, aligned_seq.length))
    end
    puts "aligned seq after trimming: #{aligned_seq} #{aligned_seq.length}" if verbose 
    return([aligned_seq, aligned_start_pos])
  end
  puts "aligned seq: #{aligned_seq} #{aligned_seq.length}" if verbose 
  return([get_seq_intervals(aligned_seq, get_mask_interval_complement(mask_intervals, aligned_seq.length)), aligned_start_pos])
end 
