#
# pasec.uger.wdl:
#  
#   A workflow to be used on the Broad Institute UGER cluster. 
#
#   Workflow summary: 
# 
#   For each bam_path in bam_paths_file, in parallel: 
#		
# 	a) sorts bam and converts bam to fastq ;
#    	b) merges paired-end reads using FLASH ;
#    	c) realigns with bwa and generates sam ;
#    	d) converts sam to bam, sorts and indexes ; and 
#    	e) computes haplotype coverage for each amplicon in bed file. 
# 
#
# Tim Farrell
# Broad Institute, IDMP, Malaria 
# tfarrell@broadinstitute.org
# 20170810
#

import "/seq/plasmodium/tfarrell/rtss/code/pasec/src/wdl/generate_haplotypes.wdl" as generate_haplotypes

## WORKFLOW DEFINITION
workflow pasec {
	## inputs
	# data
	File ref
	String output_dir
	File bam_paths_file
	File amplicon_bed_file
	Array[String] amplicons
	File run_seq_metadata
	File run_sample_metadata
	File? mask_file
	File? known_haplotypes
	Map[String, String] amplicon_bed_index_map

	# output params
	Boolean print_reads

	# filtering/ clustering configs
	# each a map like { amplicon: value } 
	Map[String, Int] min_haplotype_cov
	Map[String, Float] min_haplotype_freq
	Map[String, Int] haplotype_cluster_dist
	Map[String, Float] haplotype_cluster_cov_ratio

	# tools 
	String bwa_aligner
	String path_to_flash
	String generate_haplotypes_script_path
	String compute_read_metrics_script_path
	String run_analysis_script_path 

	# workflow control
	Boolean realign

	# parse bam_paths_file 
	Array[File] bam_paths = read_lines(bam_paths_file) 

	if (realign == true) { 
	   	# build bwa ref index 
		call build_bwa_index { 
	     	     input: 
	     	     ref = ref 
		}

		# run alignment pipeline on each bam_path in parallel 
		scatter (bam_path in bam_paths) { 
	 		String id = basename(bam_path, '.bam')
	 		# merge paired-end	       
		 	call merge_paired_end_reads { 
		 	     input: 
		 	     id = id,
		 	     bam = bam_path, 
		 	     flash = path_to_flash,
		 	     merged_dir = (output_dir + "/merged_fastqs")
		 	}  
		 	# alignment 
			call bwa_mem { 
			     input: 
			     id = id,
			     ref = ref,
			     sa = build_bwa_index.sa, 
			     amb = build_bwa_index.amb, 
			     ann = build_bwa_index.ann, 
			     pac = build_bwa_index.pac,
			     bwt = build_bwa_index.bwt, 
			     fastq = merge_paired_end_reads.merged_fastq,
			     aligned_dir = (output_dir + "/alignment")
			}
		}
	}

	Array[File] aligned_bams = select_first([bwa_mem.aligned_bam, bam_paths])

	# loop over amplicons 
	scatter (amplicon in amplicons) { 
		# generate haplotypes
		call generate_haplotypes.generate_amplicon_haplotypes { 
			input: 
			amplicon = amplicon, 
			bed = amplicon_bed_file,
			bam_paths = aligned_bams,
			bwa_aligner = bwa_aligner,  
			output_dir = (output_dir + "/haplotypes"),
			mask_file = mask_file, 
			print_reads = print_reads, 
			bed_index = amplicon_bed_index_map[amplicon],
			min_haplotype_cov = min_haplotype_cov[amplicon],
			min_haplotype_freq = min_haplotype_freq[amplicon],
			haplotype_cluster_dist = haplotype_cluster_dist[amplicon],
			haplotype_cluster_cov_ratio = haplotype_cluster_cov_ratio[amplicon],
			generate_haplotypes_script_path = generate_haplotypes_script_path, 
			compute_plot_read_metrics_path = compute_read_metrics_script_path
		}
		# analyze run 
		call analyze_run { 
			input: 
			amplicon = amplicon, 
			seq_metadata = run_seq_metadata,
			sample_metadata = run_sample_metadata, 
			output_dir = (output_dir + "/haplotypes"),
			known_haplotypes = known_haplotypes, 
			amplicon_haplotypes_files = generate_amplicon_haplotypes.haplotypes_files, 
			analysis_script_path = run_analysis_script_path
		}
	}
	# pool different amplicon analyses together
	call pool_amplicon_analyses { 
		input: 
		output_dir = output_dir,
		seq_index_files = analyze_run.seq_index, 
		seq_stats_files = analyze_run.seq_stats,
		filter_summary_files = analyze_run.filter_summary,
		haplotype_index_files = analyze_run.haplotype_index,
		read_metrics_files = flatten(generate_amplicon_haplotypes.read_metrics_files)
	}

	output { 
	       File seq_index = pool_amplicon_analyses.pooled_seq_index
	       File seq_stats = pool_amplicon_analyses.pooled_seq_stats
	       File filter_summary = pool_amplicon_analyses.pooled_filter_summary
	       File haplotype_index_file = pool_amplicon_analyses.pooled_haplotype_index
	       File read_metrics_file = pool_amplicon_analyses.read_metrics_file
	}	 
}


## TASK DEFINITIONS
# build index of reference file for bwa
task build_bwa_index { 
     File ref 

     command {
     	     source /broad/software/scripts/useuse
	     use BWA
	     bwa index ${ref}
     } 

     output { 
     	    File sa = "${ref}.sa"
	    File amb = "${ref}.amb"
	    File ann = "${ref}.ann"
	    File pac = "${ref}.pac"
	    File bwt = "${ref}.bwt"
     } 
} 

# flash task 
# merge paired end reads using FLASH (https://ccb.jhu.edu/software/FLASH/)
task merge_paired_end_reads { 
	File bam
	String id
	String flash
	String merged_dir
	
	command { 
		# load tools
		source /broad/software/scripts/useuse
		use Samtools 
		use .zlib-1.2.8 
		# mk dir if not there 
		if [ ! -d "${merged_dir}" ]; then 
			mkdir -p ${merged_dir}
		fi 
		# sort bam and bam2fq
		samtools sort -o ${id}.sorted.bam ${bam}
		samtools bam2fq ${id}.sorted.bam > ${id}.raw.fastq
		# merge
		${flash} -I ${id}.raw.fastq -o ${id} -d ${merged_dir} -M 200
	} 

	output { 
		File merged_fastq = "${merged_dir}/${id}.extendedFrags.fastq"
	} 
} 

# bwa tasks
task bwa_aln { 
	File ref
	File sa 
	File amb
	File ann 
	File pac 
	File bwt 
	String id
	File fastq
	String aligned_dir

	command {
		# load bwa
		source /broad/software/scripts/useuse
		use BWA
		use Samtools
		# mk dir if not there 
		if [ ! -d "${aligned_dir}" ]; then 
			mkdir -p ${aligned_dir}
		fi
		# bwa aln 
		bwa aln ${ref} ${fastq} -f ${id}.sai
		bwa samse ${ref} ${id}.sai ${fastq} > ${id}.sam
		# sam2bam, sort and index
		samtools view ${id}.sam -bS > ${id}.unsorted.bam
		samtools sort -o ${aligned_dir}/${id}.bam ${id}.unsorted.bam
		samtools index ${aligned_dir}/${id}.bam
	} 

	output { 
		File aligned_bam = "${aligned_dir}/${id}.bam"
		File aligned_bam_index = "${aligned_dir}/${id}.bam.bai"
	} 
}  

task bwa_mem { 
	File ref
	File sa 
	File amb
	File ann 
	File pac 
	File bwt
	String id
	File fastq
	String aligned_dir
	
	command {
		# load bwa/ samtools 
		source /broad/software/scripts/useuse
		use BWA
		use Samtools 
		# mk dir if not there 
		if [ ! -d "${aligned_dir}" ]; then 
			mkdir -p ${aligned_dir}
		fi
		# bwa mem 
		bwa mem -M ${ref} ${fastq} > ${id}.sam
		# sam2bam, sort and index
		samtools view ${id}.sam -bS > ${id}.unsorted.bam
		samtools sort -o ${aligned_dir}/${id}.bam ${id}.unsorted.bam
		samtools index ${aligned_dir}/${id}.bam
	} 

	output { 
		File aligned_bam = "${aligned_dir}/${id}.bam"
		File aligned_bam_index = "${aligned_dir}/${id}.bam.bai"
	} 
} 

task analyze_run { 
	String run_id
	String amplicon
	String output_dir
	String seq_metadata
	String sample_metadata
	File? known_haplotypes
	Float? min_population_freq
	String analysis_script_path
	Array[File] amplicon_haplotypes_files

	String raw_haplotypes_file = output_dir + "/run.haplotypes." + amplicon + ".raw.tsv"

	command { 
		# load python 
		source /broad/software/scripts/useuse
		use Anaconda
		source activate updated_packages_env
		# make dir, if doesn't exist
		if [ ! -d "${output_dir}" ]; then 
			mkdir -p ${output_dir}
		fi
		# pool raw seq haplotypes for amplicon 
		head -n 1 ${amplicon_haplotypes_files[5]} > ${raw_haplotypes_file}
		tail -q -n +2 ${sep=" " amplicon_haplotypes_files} >> ${raw_haplotypes_file}
		rm ${sep=" " amplicon_haplotypes_files}
		# run 
		python ${analysis_script_path} \
			--run_id ${run_id} \
			--amplicon ${amplicon} \
			--run_metadata ${seq_metadata} \
			--sample_metadata ${sample_metadata} \
			--known_haplotypes_file ${default='""' known_haplotypes} \
			--min_population_freq ${default="0.00001" min_population_freq} \
			--run_haplotype_coverage_file ${raw_haplotypes_file} \
			--output_dir ${output_dir} 
	}

	output { 
		File seq_index = "${output_dir}/run.${amplicon}.index.tsv"
		File? seq_stats = "${output_dir}/run.${amplicon}.stats.tsv"
		File filter_summary = "${output_dir}/${amplicon}.filter.cluster.summary.tsv"
		File haplotype_index = "${output_dir}/haplotypes.${amplicon}.index.tsv"
	}
}

task pool_amplicon_analyses {
	String output_dir  
	Array[File] seq_index_files
	Array[File?] seq_stats_files
	Array[File?] read_metrics_files
	Array[File] filter_summary_files
	Array[File] haplotype_index_files

	command { 
		# pool solexa indexes
		head -n 1 ${seq_index_files[0]} > ${output_dir}/run.index.tsv 
		tail -q -n +2 ${sep=" " seq_index_files} >> ${output_dir}/run.index.tsv
		# pool solexa stats
		head -n 1 ${seq_stats_files[0]} > ${output_dir}/run.stats.tsv 
		tail -q -n +2 ${sep=" " seq_stats_files} >> ${output_dir}/run.stats.tsv
		# pool filter summary 
		head -n 1 ${filter_summary_files[0]} > ${output_dir}/filter.cluster.summary.tsv 
		tail -q -n +2 ${sep=" " filter_summary_files} >> ${output_dir}/filter.cluster.summary.tsv
		# pool haplotype indexes
		head -n 1 ${haplotype_index_files[0]} > ${output_dir}/haplotypes.index.tsv 
		tail -q -n +2 ${sep=" " haplotype_index_files} >> ${output_dir}/haplotypes.index.tsv
		# pool read metrics 
		head -n 1 ${read_metrics_files[0]} > ${output_dir}/read_metrics.csv
		tail -q -n +2 ${sep=" " read_metrics_files} >> ${output_dir}/read_metrics.csv 
	}

	output { 
		File pooled_seq_index = "${output_dir}/run.index.tsv"
		File pooled_seq_stats = "${output_dir}/run.stats.tsv"
		File pooled_filter_summary = "${output_dir}/filter.cluster.summary.tsv"
		File pooled_haplotype_index = "${output_dir}/haplotypes.index.tsv"
		File read_metrics_file = "${output_dir}/read_metrics.csv"
	}
}
