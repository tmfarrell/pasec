#!/usr/bin/env python
##
##  compute_plot_read_metrics.py
##
##  computes and plots base, mismatch, indel and quality metrics
##  for reads/ haplotypes that were filtered, clustered and included 
##  in the end results. 
##
##  Tim Farrell, tfarrell@broadinstitute.org
##  IDMP, Malaria
##  20180717
## 

import os
import argparse
#import numpy as np
#import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
curr_dir = os.path.dirname(os.path.realpath(__file__))
exec(open(os.path.join(curr_dir, '../process_haplotypes/utilities.py')).read())

# methods 
def min_freq(freqs): 
    return(max([v for v in freqs if (v != max(freqs))]))

def quality_str_to_arr(qstr): 
    return([ord(c)-33 for c in qstr])

def parse_solexa_amplicon_from_fname(fname): 
	file_prefix = os.path.splitext(os.path.basename(fname))[0]
	solexa = file_prefix[:file_prefix.find('.')]
	amplicon = file_prefix[file_prefix.find('.')+1:file_prefix.rfind('.')]
	return(solexa, amplicon)

def assign_haplotype_filter_class(metrics): 
    if np.prod(metrics[['pass_coverage','pass_freq','pass_len']]) != 1: 
        return('filtered')
    elif metrics['haplotype_cluster_idx'] != metrics['haplotype_idx']: 
        return('clustered')
    else: 
        return('result')

def plot_quality(df, title=None, save_as=None, verbose=False,
                 filter_classes=['filtered','clustered','result']): 
    Q = pd.DataFrame()
    f, axes = plt.subplots(len(filter_classes), 1, sharex=True, sharey=True)
    for filter_class, ax in zip(filter_classes, axes): 
        qual = (pd.DataFrame(df[df.filter_class == filter_class]['qual']
                             .apply(quality_str_to_arr).values.tolist()))
        if not qual.empty: 
            avg = qual.mean()
            std = qual.std()
            ax.plot(qual.columns, avg, alpha=0.75)
            ax.fill_between(qual.columns, avg - std, avg + std, alpha=0.5)
            ax.set_title(filter_class)
            avg = avg.to_frame().transpose()
            avg['filter_class'] = [filter_class]
            avg['metric'] = ['quality']
            Q = Q.append(avg)
    #Q['filter_class'] = filter_classes
    #Q['metric'] = ['quality'] * len(filter_classes)
    plt.suptitle('quality' if not title else title)
    f.set_size_inches((10, 6))
    if save_as: 
        plt.savefig(save_as)
    plt.close(f)
    return(Q)

def plot_freqs(df, freq_type='mismatch', title=None, save_as=None, verbose=False,
               filter_classes=['filtered','clustered','result']): 
    ops_to_plot = {'mismatch': ['X'], 'cigar': ['I','D'],
                   'read': ['major','minor','alt']}[freq_type]
    F = pd.DataFrame()
    f, axes = plt.subplots(len(filter_classes), 1, sharex=True, sharey=True)
    for filter_class, ax in zip(filter_classes, axes): 
        cts = (pd.DataFrame(df[df.filter_class == filter_class]
                            ['long_' + freq_type if freq_type in ['cigar','mismatch'] else freq_type]
                            .apply(list).values.tolist())
               .apply(lambda col: col.value_counts()).fillna(0.0))
        freqs = cts.apply(lambda col: [float(v)/float(col.sum()) for v in col.values])
        if verbose:
            print(freq_type, filter_class)
            print('counts:')
            print(cts)
            print('freqs:')
            print(freqs)
        if not freqs.empty: 
            if freq_type == 'read': 
                freqs.loc['major',:] = freqs.apply(lambda col: col.max())
                freqs.loc['minor',:] = freqs.apply(lambda col: min_freq(col.values))
                freqs.loc['alt',:] = 1.0 - (freqs.loc['major',:] + freqs.loc['minor',:])
                freqs = freqs.loc[['major','minor','alt'],:]
            for op in ops_to_plot:
                if op in freqs.index: 
                    ax.plot(freqs.columns, freqs.loc[op,:].values, label=op, alpha=0.6)
                else: 
                    freqs.loc[op,:] = [0.0] * len(freqs.columns)
                    ax.plot(freqs.columns, freqs.loc[op,:].values, label=op, alpha=0.6)
            ax.set_title(filter_class, fontsize='medium')
            freqs['filter_class'] = [filter_class] * len(freqs.index)
            F = F.append(freqs.copy())
    axes[0].legend()
    plt.ylim(0.0, 1.0)
    plt.suptitle(freq_type + ' freqs' if not title else title)
    f.set_size_inches((10, 6))
    if save_as: 
        plt.savefig(save_as)
    plt.close(f)
    return(F)


# main 
def main(): 
	# parse args
	parser = argparse.ArgumentParser(description="computes and plots base, mismatch, indel and quality metrics " +\
						     "for reads/ haplotypes that were filtered, clustered or included in the end results.")
	parser.add_argument('--reads_files', nargs='+', required=True, 
						help="Path to reads files. Expected in format for each: /path/to/${solexa}.${amplicon}.reads.tsv.")
	parser.add_argument('--haplotypes_files', nargs='+', required=True, 
					    help="Path to haplotypes files, corresponding to each of the reads files. " +
					         "Expected in format for each: /path/to/${solexa}.${amplicon}.haplotypes.tsv.")
	parser.add_argument('--output_dir', default=None, help="Path to directory to save output. Default is directory of reads file.")
	parser.add_argument('--plot_dir', default=None, help="Path to directory to save plots. Default is output_dir.")
	parser.add_argument('--do_not_plot', action='store_true', help="Flag to pass if do not want to save plots.")
	args = parser.parse_args()
	output_dir = args.output_dir if (not args.output_dir == None) else os.path.dirname(args.reads_files[0])
	plot_dir = args.plot_dir if (not args.plot_dir == None) else output_dir
	# for each reads file 
	for reads_file, haplotypes_file in zip(args.reads_files, args.haplotypes_files):
		solexa_r, amplicon_r = parse_solexa_amplicon_from_fname(reads_file)
		solexa_h, amplicon_h = parse_solexa_amplicon_from_fname(haplotypes_file)
		assert solexa_r == solexa_h and amplicon_r == amplicon_h, \
			"solexa and amplicon of reads and haplotypes must be the same."
		solexa = solexa_r ; amplicon = amplicon_r
		# compute and plot
		#print(solexa, amplicon)
		try: 
			reads = pd.read_table(reads_file)
			haplotypes = pd.read_table(haplotypes_file)
		except: 
			continue
		if not reads[reads.len < 350].empty:
			read_metrics = pd.DataFrame()
			# label reads/haplotypes as either filtered/clustered/results
			haplotypes['filter_class'] = haplotypes.apply(lambda row: assign_haplotype_filter_class(row), axis=1)
			haplotype_read_idx = haplotypes[['amplicon','read_idx','haplotype_idx','filter_class']]
			haplotype_read_idx = haplotype_read_idx.rename(columns={'read_idx':'best_read_idx'})
			reads = reads[reads.len < 350]
			reads = reads.merge(haplotype_read_idx, on=['amplicon','haplotype_idx'], how='left')
			reads['read_filter_class'] = (reads.apply(lambda row: 'collapsed' 
			                                          if row['best_read_idx'] != row['read_idx']
			                                          else 'haplotype', axis=1))
			# add read cols
			reads['long_mismatch'] = reads['mismatch'].apply(lambda m: mismatch_str_to_long_mismatch(m))
			reads['long_cigar'] = reads['cigar'].apply(lambda c: cigar_arr_to_long_cigar(to_cigar_arr(c)))
			# for each metric
			for freq_type in ['read','mismatch','cigar']:
			    # plot/ compute metrics 
			    F = plot_freqs(reads, freq_type, title=(freq_type + ' freqs: ' + amplicon + " " + solexa), 
			                   save_as=(os.path.join(plot_dir, '.'.join([solexa,freq_type+'_freqs',amplicon,'png']))
			                   	    if not args.do_not_plot else None), verbose=True)
			    # add solexa id and append to read metrics
			    F['Solexa_ID'] = [solexa] * len(F.index)
			    F['amplicon'] = [amplicon] * len(F.index)
			    read_metrics = read_metrics.append(F.reset_index().rename(columns={'index':'metric'}))
			# plot/ compute quality 
			Q = plot_quality(reads, title=('quality: ' + amplicon + ' ' + solexa), 
			                 save_as=(os.path.join(plot_dir,'.'.join([solexa,'quality',amplicon,'png']))
			                 	  if not args.do_not_plot else None), verbose=True)
			Q['Solexa_ID'] = [solexa] * len(Q.index)
			Q['amplicon'] = [amplicon] * len(Q.index)
			read_metrics = read_metrics.append(Q)
			read_metrics = read_metrics[[c for c in read_metrics.columns.tolist() if type(c) == str] + \
                            			[c for c in read_metrics.columns.tolist() if type(c) != str]].fillna(0.0)
			read_metrics.to_csv(os.path.join(output_dir, '.'.join([solexa,amplicon,'read_metrics.csv'])), index=False)
	return


if __name__ == '__main__': 
	main()
