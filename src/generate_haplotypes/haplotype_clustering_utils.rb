#!/usr/bin/env ruby

require 'levenshtein'

def hamming_distance(s1, s2)
=begin
     For two strings of equal length, returns number of positions at which they differ,
     or equivalently the number of substitutions needed to turn one of the strings into the other
     input:    s1 :: String, s2 :: String
     output:   hamming_dist :: Integer 

     For strings of unequal length, use Levenshtein.distance(s1, s2) 
=end
  if s1.length != s2.length 
    raise(Exception, 'Strings not of equal length.')
  end 
  return (0...s1.length).map{|i| s1[i] != s2[i] ? 1 : 0}.reduce(:+)
end

## haplotype clustering functions
def collect_subclusters(haplotype_clusters, key, verbose=false)
=begin 
   Recursively consolidate subclusters    
=end
  puts "collecting subclusters for #{key}" if verbose
  if ! haplotype_clusters.keys.include?(key)
    return []
  end
  cluster_haplotypes = haplotype_clusters[key].map{|h, dist| h}
  puts "cluster haplotypes: #{cluster_haplotypes}" if verbose
  if cluster_haplotypes == []
    return []
  else 
    collected_haplotype_cluster = haplotype_clusters[key].dup
    puts "collected haplotype cluster before: #{collected_haplotype_cluster}" if verbose 
    for j in cluster_haplotypes
      # recursive call 
      collected_haplotype_cluster = collected_haplotype_cluster.concat(collect_subclusters(haplotype_clusters, key=j))
    end
    puts "collected haplotype cluster after: #{collected_haplotype_cluster}" if verbose 
    return collected_haplotype_cluster
  end
end 

def build_haplotype_cluster_graph(haplotype_set, coverage_ratio_threshold, mismatch_threshold, 
                                  collect_subclusters=false, min_clusterability_freq=0.005, 
                                  verbose=false, print_clustering_metrics=false)
=begin
  input: 
    haplotype_set :: [{:haplotype => String, :coverage => Integer}]
      which is sorted by descending coverage

  output:    
    haplotype_clusters :: { haplotype => (haplotype_cluster :: [haplotype]) }
      where haplotype_cluster is list of haplotypes with coverage_ratio <= coverage_ratio_threshold, 
      mismatch < mismatch_threshold in relation to haplotype key and which contains no haplotypes 
      of any other haplotype's cluster (i.e. haplotypes belong to at most one cluster) 
         
  This can be thought of as building a graph where: (a) haplotypes are nodes; (b) nodes that are 
  less than Math.sqrt(coverage_ratio_threshold**2 + mismatch_threshold**2) from each other in 
  coverage_ratio-mismatch 2D-space have a directed edge b/t them (minor to major); (c) all nodes 
  have at most 1 out-edge where its out-edge connects it to the closest node and, if that minimum
  distance is not unique, that edge points to the node with the higher coverage. During "contraction", 
  all nodes with out-edges are clustered to the major node they point to 
=end  
  haplotype_clusters = Hash.new
  total_coverage = haplotype_set.map{|haplotype| haplotype[:coverage]}.reduce(:+).to_f
  for i in 0...haplotype_set.length
    id1 = haplotype_set[i][:haplotype_idx] 
    haplotype_clusters[id1] = [] 
    coverage1 = haplotype_set[i][:coverage].to_f
    # don't assess "clusterability" for haplotypes w/ freq below a certain threshold
    if (coverage1 / total_coverage) > min_clusterability_freq
      haplotype1 = haplotype_set[i][:haplotype]
      for j in (i+1)...haplotype_set.length 
        # treat id1 as major haplotype
        coverage_ratio = (haplotype_set[j][:coverage].to_f / coverage1)
        haplotype2 = haplotype_set[j][:haplotype]
        if haplotype1.length == haplotype2.length
          mismatches = hamming_distance(haplotype1, haplotype2)
        else
          mismatches = Levenshtein.distance(haplotype1, haplotype2) 
        end
        puts "coverage_ratio(#{id1}, #{j.to_s}) = #{coverage_ratio}" if print_clustering_metrics
        puts "edit_distance(#{id1}, #{j.to_s}) = #{mismatches}" if print_clustering_metrics
        if mismatches <= mismatch_threshold and coverage_ratio <= coverage_ratio_threshold
          puts "#{j.to_s} had edit_distance <= #{mismatch_threshold} and " + \
               "coverage_ratio <= #{coverage_ratio_threshold} with #{id1}" if print_clustering_metrics
          dist = Math.sqrt(mismatches**2 + coverage_ratio**2)
          haplotype_clusters[id1].push([haplotype_set[j][:haplotype_idx], dist])
        end 
      end
    end
  end
  # drop those haplotypes with empty clusters
  haplotype_clusters.delete_if{|k, v| v == []}
  puts "Haplotype clusters:\n#{haplotype_clusters}" if verbose
  # resolve intersections between haplotype clusters, keeping edges with minimum distance 
  for i in 0...haplotype_clusters.length
    id1 = haplotype_clusters.keys[i]
    cluster1 = haplotype_clusters[id1]
    cluster1_haplotypes = cluster1.map{|h, dist| h}
    for j in (i+1)...haplotype_clusters.length
      id2 = haplotype_clusters.keys[j]
      cluster2 = haplotype_clusters[id2]
      cluster2_haplotypes = cluster2.map{|h, dist| h}
      intersection = cluster1_haplotypes & cluster2_haplotypes
      for x in intersection
        # get haplotype cluster id to filter, based on which has larger dist to x 
        _, dist1 = cluster1.select{|h, dist| h == x}[0]
        _, dist2 = cluster2.select{|h, dist| h == x}[0]
        filter_id = dist1 <= dist2 ? id2 : id1
        # filter the one with larger distance, giving ties to id1 (higher coverage) 
        haplotype_clusters[filter_id] = haplotype_clusters[filter_id].select{|h, dist| h != x}
      end
    end
  end
  # drop those haplotypes with empty clusters, again
  haplotype_clusters.delete_if{|k, v| v == []}
  if ! collect_subclusters
    puts "Haplotype clusters before filtering haplotypes that have their own cluster:\n#{haplotype_clusters}" if verbose
    # drop those haplotypes in haplotype clusters, if they themselves have clusters
    for haplotype in haplotype_clusters.keys
      haplotype_clusters[haplotype] = haplotype_clusters[haplotype].select{|h, dist| !haplotype_clusters.keys.include?(h)}
    end
    puts "Haplotype clusters after filtering haplotypes that have their own cluster:\n#{haplotype_clusters}" if verbose
    return haplotype_clusters
  else
    # if haplotypes identified to be contracted have clusters themselves
    # put their cluster contents into the higher-coverage cluster, recursively
    puts "Clusters before collecting subclusters:\n#{haplotype_clusters}" if print_clustering_metrics
    collected_haplotype_clusters = Hash.new
    for i in haplotype_clusters.keys
      collected_haplotype_clusters[i] = collect_subclusters(haplotype_clusters, key=i, 
                                                            verbose=print_clustering_metrics)
    end
    # drop those haplotypes with empty clusters, again
    collected_haplotype_clusters.delete_if{|k, v| v == []}
    puts "Clusters after collecting subclusters:\n#{collected_haplotype_clusters}" if verbose 
    return collected_haplotype_clusters
  end
end

def apply_haplotype_cluster_contraction(haplotype_set, coverage_ratio_threshold=0.125, mismatch_threshold=1,
                                        sum_cluster_coverage=true, collect_subclusters=false, verbose=false, 
                                        print_clustering_metrics=false, min_clusterability_freq=0.005)
=begin
    Merges/ contracts low-coverage haplotypes into highly-similar high-coverage haplotypes, 
    using a distance-based clustering, in an effort to correct for errors and reduce noise    
    
    input:   haplotype_set :: [{:haplotype => String, :coverage => Integer}]

    output:  contracted_haplotype_set :: [{:haplotype => String, :coverage => Integer}]
=end
  if haplotype_set == []
    return []
  end 
  # make sure haplotype_set is sorted by coverage, descending
  haplotype_set.sort!{|h1, h2| h2[:coverage] <=> h1[:coverage]}
  # build haplotype cluster graph 
  haplotype_clusters = build_haplotype_cluster_graph(haplotype_set, coverage_ratio_threshold, mismatch_threshold,
                                                     collect_subclusters=collect_subclusters, min_clusterability_freq=min_clusterability_freq, 
                                                     verbose=verbose, print_clustering_metrics=print_clustering_metrics)
  # build contracted haplotype set 
  contracted_haplotype_set = []
  haplotypes_to_be_contracted = haplotype_clusters.values.flatten(1).map{|h, dist| h}
  puts "Haplotypes to be contracted:\n#{haplotypes_to_be_contracted}" if print_clustering_metrics 
  for i in 0...haplotype_set.length
    id1 = i.to_s
    if ! haplotypes_to_be_contracted.include?(id1)
      haplotype_data = haplotype_set[i]
      coverage = haplotype_data[:coverage]
      haplotype = haplotype_data[:haplotype]
      if haplotype_clusters.keys.include?(id1) 
        for j, _ in haplotype_clusters[id1]
          if sum_cluster_coverage
            coverage = coverage + haplotype_set[j.to_i][:coverage]
          end
        end
      end
      contracted_haplotype_set.push(haplotype_data.merge({:coverage => coverage, :haplotype => haplotype}))
    end
  end
  return([contracted_haplotype_set, haplotype_clusters]) 
end 

## main 
if __FILE__ == $0 
  require_relative './utilities.rb'
  # parse read file as haplotypes
  reads = parse_reads_file(ARGV[0]).map{|read| read.update({:haplotype => read[:read]})}
  # cluster 
  clustered_haplotypes, clusters = apply_haplotype_cluster_contraction(reads, coverage_ratio_threshold=0.125, mismatch_threshold=1,
                                                                       sum_cluster_coverage=true, collect_subclusters=false) 
  # print reads 
  #reads.each{|read|
  #  puts "#{read[:read_index]}\t#{read[:coverage]}\t#{read[:cigar]}\n"
  #}
  #puts "# reads: #{reads.length}"
  # print clustered haplotypes
  clustered_haplotypes.each{|haplotype|
    if clusters[haplotype[:read_index][1..-1]]
      clustered = clusters[haplotype[:read_index][1..-1]].map{|id,_| "R" + id}.join('/')
    else
      clustered = ""
    end
    puts "#{haplotype[:read_index]}\t#{haplotype[:coverage]}\t#{haplotype[:cigar]}\t#{clustered}\n"
  }
  #puts "# clustered haplotypes: #{clustered_haplotypes.length}"
end
