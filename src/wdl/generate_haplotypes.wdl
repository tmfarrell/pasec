##
##	generate_haplotypes.wdl: 
## 	 - A subworkflow to be called by amplicon_seq.wdl 
##
## 	tfarrell@broadinstitute.org
##  	20180312 
## 

workflow generate_amplicon_haplotypes { 
	
	File bed
	String amplicon
	String bed_index
	String output_dir 
	Array[File] bam_paths
	String generate_haplotypes_script_path
	String compute_plot_read_metrics_path

	String bwa_aligner
	Boolean print_reads
	Int min_haplotype_cov
	Float min_haplotype_freq
	Int haplotype_cluster_dist
	Float haplotype_cluster_cov_ratio

	scatter (bam_path in bam_paths) { 

		call generate_haplotypes { 
			input:
			bed = bed,
			bam = bam_path, 
			amplicon = amplicon, 
			bed_index = bed_index,
			print_reads = print_reads, 
			output_dir = output_dir,
			bam_index = bam_path + ".bai",
			bwa_aligner = bwa_aligner, 
			min_haplotype_cov = min_haplotype_cov,
			min_haplotype_freq = min_haplotype_freq,
			haplotype_cluster_dist = haplotype_cluster_dist,
			haplotype_cluster_cov_ratio = haplotype_cluster_cov_ratio,
			script_path = generate_haplotypes_script_path
		}

		call compute_read_metrics { 
		     input:
		     amplicon = amplicon, 
		     output_dir = output_dir, 
		     reads_file = generate_haplotypes.reads_file, 
		     haplotypes_file = generate_haplotypes.haplotypes_file,
		     script_path = compute_plot_read_metrics_path
		} 
	}

	output { 
		Array[File] haplotypes_files = generate_haplotypes.haplotypes_file
		Array[File] reads_files = generate_haplotypes.reads_file
		Array[File] read_metrics_files = compute_read_metrics.read_metrics_file
	} 
}

task generate_haplotypes { 
	File bam
	File bed
	File bam_index
	String bed_index 
	String amplicon
	String output_dir
	String script_path

	String bwa_aligner
	Boolean print_reads 
	String min_haplotype_cov
	String min_haplotype_freq
	Int haplotype_cluster_dist
	Float haplotype_cluster_cov_ratio 

	String id = basename(bam, ".bam")

	command { 
		# load ruby and samtools
		source /broad/software/scripts/useuse
		use Ruby-1.9
		use Samtools
		# mk dir if not there
		if [ ! -d "${output_dir}" ]; then 
			mkdir -p ${output_dir} 
		fi
		if ${print_reads} ; then 
		   ruby ${script_path} haplotypes \
			--bam ${bam} \
			--print_reads --print_to ${output_dir} \
			--bwa_aligner_used ${bwa_aligner} \
			--bed ${bed} --bed_index ${bed_index} \
			--min_haplotype_coverage ${min_haplotype_cov} \
			--min_haplotype_freq ${min_haplotype_freq} \
			--haplotype_clustering-edit_dist ${haplotype_cluster_dist} \
			--haplotype_clustering-coverage_ratio ${haplotype_cluster_cov_ratio}
		else
		   ruby ${script_path} haplotypes \
			--bam ${bam} \
			--bwa_aligner_used ${bwa_aligner} \
			--bed ${bed} --bed_index ${bed_index} \
			--min_haplotype_coverage ${min_haplotype_cov} \
			--min_haplotype_freq ${min_haplotype_freq} \
			--haplotype_clustering-edit_dist ${haplotype_cluster_dist} \
			--haplotype_clustering-coverage_ratio ${haplotype_cluster_cov_ratio} \
			> ${output_dir}/${id}.${amplicon}.haplotypes.tsv	
		fi
		if [ ! -f ${output_dir}/${id}.${amplicon}.haplotypes.tsv ] ; then 
		   touch ${output_dir}/${id}.${amplicon}.haplotypes.tsv
		fi	 
		if [ ! -f ${output_dir}/${id}.${amplicon}.reads.tsv ] ; then
		   touch ${output_dir}/${id}.${amplicon}.reads.tsv
		fi 
	} 

	output { 
		File haplotypes_file = "${output_dir}/${id}.${amplicon}.haplotypes.tsv"
		File reads_file = "${output_dir}/${id}.${amplicon}.reads.tsv"
	}
} 

task compute_read_metrics { 
     String amplicon 
     String output_dir 
     String script_path     

     File reads_file
     File haplotypes_file

     String id = basename(reads_file, "." + amplicon + ".reads.tsv")
     String outfile = output_dir + "/" + id + "." + amplicon + ".read_metrics.csv"
     
     command { 
     	     # load python w/ updated packages
     	     source /broad/software/scripts/useuse
	     use Anaconda
	     source activate updated_packages_env
	     # run 
	     python ${script_path} \
	     	     --reads_file ${reads_file} \
		     --haplotypes_file ${haplotypes_file} \
		     --output_dir ${output_dir} \
		     --do_not_plot
	     if [ ! -f ${outfile} ] ; then 
	         touch ${outfile}
	     fi 	     
     } 

     output { 
     	    File read_metrics_file = "${outfile}"
     }    
} 
