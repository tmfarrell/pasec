#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np
import pandas as pd
code_dir = os.path.dirname(__file__)
exec(open(os.path.join(code_dir, 'utilities.py')).read())

## stats computations
def compute_summary(haplotype_cov_df, step):
    return({'step': step, 
            'amplicon': haplotype_cov_df.iloc[0]['amplicon'],
            'total_coverage': haplotype_cov_df.coverage.sum(),
            'num_samples': haplotype_cov_df.sample_id.nunique(), 
            'total_num_haplotypes': haplotype_cov_df.haplotype.nunique(),
            'avg_sample_coverage': haplotype_cov_df.groupby('sample_id').coverage.sum().mean(),
            'avg_num_haplotypes': haplotype_cov_df.groupby('sample_id').haplotype.count().mean()})
           
def compute_haplotype_stats(run_haplotype_coverage, amplicon, min_population_freq=None):
    # groupby haplotype
    try: 
        haplotype_index = (run_haplotype_coverage.groupby('haplotype')
                    .agg({'coverage': np.sum, 
                          'locus': lambda l: str(set(l).pop()),
                          'amplicon': lambda l: str(set(l).pop())})
                    .reset_index()[['coverage','amplicon','locus','haplotype']])
    except: 
        haplotype_index = (run_haplotype_coverage.groupby('haplotype')
                    .agg({'coverage': np.sum})
                    .reset_index()[['coverage','haplotype']])
    #print run_haplotype_coverage
    #sys.stderr.write(haplotype_index.to_csv(sep='\t'))
    # sort and add index 
    haplotype_index = haplotype_index.sort_values('coverage', ascending=False).reset_index(drop=True)
    haplotype_index['haplotype_index'] = haplotype_index.index.values
    haplotype_index['amplicon_index'] = [amplicon + '-' + str(i) for i in haplotype_index.index.values]
    # add coverage frequencies
    total_coverage = float(haplotype_index.coverage.sum())
    haplotype_index['freq'] = haplotype_index.coverage.apply(lambda c: float(c) / total_coverage)
    if not min_population_freq is None: 
        # population-level filter
        haplotype_index['pass_population_freq'] = haplotype_index.freq.apply(lambda x: 1 if x > min_population_freq else 0)
    if not min_population_freq is None: 
        # index run_haplotypes with haplotype_index
        run_haplotype_coverage_indexed = run_haplotype_coverage.merge(haplotype_index[['haplotype','haplotype_index','amplicon_index','pass_population_freq']], on='haplotype', how='left')
    else: 
        run_haplotype_coverage_indexed = run_haplotype_coverage.merge(haplotype_index[['haplotype','haplotype_index','amplicon_index']], on='haplotype', how='left')
    del run_haplotype_coverage_indexed['haplotype']
    #haplotype_index['haplotype'] = haplotype_index['haplotype']
    #print run_haplotype_coverage
    try: 
        # add run/sample_indices with coverages
        haplotype_index = haplotype_index.merge((run_haplotype_coverage_indexed.groupby('haplotype_index')
                                                 .apply(lambda df: metrics_to_str(dict(zip(df['sample_index'], df['coverage']))))
                                                 .reset_index()
                                                 .rename(columns={0:'sample_index'})), on='haplotype_index', how='left')
        #print haplotype_index
        # add sample counts
        haplotype_index['sample_ct'] = haplotype_index.sample_index.apply(lambda s: s.count('/') + 1)
        # add sample percents
        num_samples = run_haplotype_coverage['sample_id'].nunique()
        haplotype_index['sample_pct'] = haplotype_index.sample_ct.apply(lambda c: float(c) / num_samples)
    except: pass
    # add start pos
    try: 
        haplotype_index['start_pos'] = (haplotype_index.locus
                                        .apply(lambda l: int(l.split(':')[1].split('-')[0]) if '-' in l 
                                               else min([int(s.split(':')[1]) for s in l.split('/')])))
    except: pass 
    return(haplotype_index, run_haplotype_coverage_indexed)

def compute_run_haplotype_stats(run_haplotype_coverage_indexed, haplotype_index, metadata, 
                                amplicon, known_haplotypes=None):
    # get data 
    run_haplotype_stats = (run_haplotype_coverage_indexed.groupby('sample_index')
                         .agg({'coverage': np.sum })
                         .reset_index()[['sample_index','coverage']])
    # add haplotype freqs 
    run_haplotype_stats = run_haplotype_stats.merge((run_haplotype_coverage_indexed.groupby('sample_index')
                                                     .apply(lambda df: metrics_to_str(cts_dict_to_freqs(dict(zip(df['haplotype_index'].values,
                                                                                                         df['coverage'])))))
                                                     .reset_index()
                                                     .rename(columns={0:'haplotype_freqs'})), 
                                                    on='sample_index', how='left')
    # add haplotype count 
    run_haplotype_stats['haplotype_ct'] = run_haplotype_stats.haplotype_freqs.apply(lambda s: len(metrics_to_dict(s).keys()))
    ## rewrite: add haplotype mismatch ## 
    # join w/ metadata
    run_haplotype_stats = run_haplotype_stats.merge(metadata, on='sample_index', how='outer')
    if not known_haplotypes is None and not known_haplotypes.empty:
        '''try: 
            # add mock sample metrics
            mock_samples = metadata[metadata.contents == 'p_falciparum_mock'].sample_index.values.tolist()
            haplotype_index['mock_coverage'] = haplotype_index.sample_index.apply(lambda s: sum(map(lambda (sample, x): x, 
                                                                                                    filter(lambda (sample, x): sample in mock_samples, 
                                                                                                           metrics_to_dict(s).items()))))
            haplotype_index['frac_mock'] = haplotype_index.apply(lambda row: safe_div(row['mock_coverage'], row['coverage']), axis=1)
            haplotype_index['in_mock'] = haplotype_index.apply(lambda row: 1 if (row['mock_coverage'] > 0) else 0, axis=1)
            haplotype_index['mock_haplotype'] = haplotype_index.apply(lambda row: 1 if (row['frac_mock'] > 0.01 or row['mock_coverage'] > 1500) else 0, axis=1)
        except: pass
        # compute strain strain assignments for haplotype index 
        haplotype_index, strain_assign_fields = assign_strains(haplotype_index, run_haplotype_stats, known_haplotypes=known_haplotypes)
        # add those assignments to run_haplotype_coverage_indexed 
        run_haplotype_coverage_indexed = run_haplotype_coverage_indexed.merge(haplotype_index[['haplotype_index'] + strain_assign_fields],
                                                                          on='haplotype_index', how='left')
        # add haplotype indices to run_haplotype_stats, w/ strain assign fields
        for strain_assign_field in strain_assign_fields: 
            run_haplotype_stats = run_haplotype_stats.merge((run_haplotype_coverage_indexed.groupby('sample_index')
                                                         .apply(lambda df: dict_to_str(dict(zip(df['haplotype_index'].values,
                                                                                                df[strain_assign_field]))))
                                                         .reset_index()
                                                         .rename(columns={0:'haplotype_to_' + strain_assign_field})), on='sample_index', how='left')
            run_haplotype_stats['haplotype_to_' + strain_assign_field] = run_haplotype_stats['haplotype_to_' + strain_assign_field].fillna('/')
            run_haplotype_stats['strains_' + strain_assign_field] = (run_haplotype_stats['haplotype_to_' + strain_assign_field]
                                                                   .apply(lambda s: ';'.join(set([rs.split(':')[1] for rs in s.split('/') if rs]))) # first get strains assigned
                                                                   .apply(lambda s: ';'.join(set(s.split(';')))))                                   # then make unique    
        try: 
            # compute strain resolution metrics
            for strain_assign_field in strain_assign_fields: 
                for metric_name, metric_fcn in [('jaccard', jaccard_index),('precision',precision_multilabel),
                                                ('recall', recall_multilabel), ('fdr',fdr),('f1_score',f1_score)]: 
                    run_haplotype_stats['_'.join([metric_name, strain_assign_field])] = (run_haplotype_stats
                                                                                       .apply(lambda row: metric_fcn(str(row['Strains']).split(';'), 
                                                                                                                     str(row['strains_' + strain_assign_field]).split(';')),
                                                                                              axis=1))
        except: pass
        '''
        # rewrite above to deal with list of provided known haplotypes (i.e. truth set) 
        pass 
    # fill NaNs appropriately
    run_haplotype_stats['coverage'] = run_haplotype_stats['coverage'].fillna(0)
    run_haplotype_stats['haplotype_ct'] = run_haplotype_stats['haplotype_ct'].fillna(0)
    run_haplotype_stats['haplotype_freqs'] = run_haplotype_stats['haplotype_freqs'].fillna("")
    try: run_haplotype_stats['haplotype_mismatch'] = run_haplotype_stats['haplotype_mismatch'].fillna("")
    except: pass
    # sort_values by run_index
    run_haplotype_stats['sample_index_#'] = run_haplotype_stats['sample_index'].apply(lambda s: int(s[1:]))
    run_haplotype_stats = run_haplotype_stats.sort_values('sample_index_#').reset_index(drop=True)
    del run_haplotype_stats['sample_index_#']
    run_haplotype_stats['amplicon'] = [amplicon] * len(run_haplotype_stats.index)
    # include haplotype_to_{strain_assign_field} and {metric}_{strain_assign_field}
    try: 
        return(run_haplotype_stats[['sample_index','amplicon','coverage'] +\
                                 filter(lambda col: any(map(lambda field: (field in col) and (not 'haplotype' in col), strain_assign_fields)), run_haplotype_stats.columns) +\
                                 filter(lambda col: 'haplotype' in col, run_haplotype_stats.columns)], 
               run_haplotype_coverage_indexed, haplotype_index)
    except: 
        return(run_haplotype_stats[['sample_index','amplicon','coverage','haplotype_ct','haplotype_freqs']],
               run_haplotype_coverage_indexed, haplotype_index)

## main 
def main():
    ## init some useful generic parsers before building subparsers
    verbose_parser = argparse.ArgumentParser(add_help=False)
    verbose_parser.add_argument('--verbose', action='store_true', help="Run verbosely.") 
    metadata_parser = argparse.ArgumentParser(add_help=False)
    metadata_parser.add_argument('--metadata_file', 
                                  help="Path to the metadata file. Should be .tsv that relates sample_id.bam to each sample's attributes.")
    haplotype_coverage_file_parser = argparse.ArgumentParser(add_help=False)
    haplotype_coverage_file_parser.add_argument('--input', required=True,
                                           help="Path to the main data file, a .tsv with 3 columns: sample_id, haplotype and coverage.")
    ## build main parser
    main_parser = argparse.ArgumentParser(parents=[haplotype_coverage_file_parser, metadata_parser, verbose_parser])
    main_parser.add_argument('--output_dir', required=True, help='Directory to save output files to.')
    main_parser.add_argument('--amplicon', required=True, help='Amplicon computing data for.')
    main_parser.add_argument('--known_haplotypes_file', default=None, help='A file listing (amplicon, strain, haplotype) for each known haplotype.')
    main_parser.add_argument('--filename_id', default='', help='Identifier to be added to all filenames.')
    main_parser.add_argument('--run_id', required=True, help='Prefix to all output files (e.g. "${run_id}.filtered.tsv").')
    main_parser.add_argument('--output_fastas', default=False, action='store_const', const=True, 
                             help='Whether to save ${amplicon}.haplotypes.index.fasta.')
    main_parser.add_argument('--min_population_freq', default=None, type=float, help="Minimum threshold frequency for valid haplotypes in the population.")
    main_parser.add_argument('--no_filter_cluster', action='store_true', 
                             help="Pass this flag to not apply filtering/ clustering.")
    main_parser.add_argument('--compute_run_stats', action='store_true', 
                             help="Compute statistics for each sample_id.")
    ## parse 
    args = main_parser.parse_args()
    ## get metadata 
    if args.metadata_file: 
        metadata = parse_metadata_file(args.metadata_file, args.verbose)
    else: 
        metadata = None 
    #metadata = pd.haplotype_csv(args.metadata)
    if args.known_haplotypes_file: 
        known_haplotypes = pd.read_table(args.known_haplotypes_file)
        if '_' in args.amplicon: 
            base_amplicon = (lambda s: s[:s.find('_')])(args.amplicon)
        else: 
            base_amplicon = args.amplicon
        #print(base_amplicon)
        known_haplotypes = known_haplotypes[known_haplotypes.amplicon == base_amplicon]
        #print(known_haplotypes)
    else: 
        known_haplotypes = None
    ## get run haplotype coverage data
    run_haplotype_coverage = parse_run_haplotype_coverage_file(args.input)
    if run_haplotype_coverage.empty: 
        if args.verbose: print("Input run haplotypes file is empty.")
        return
    if not 'amplicon' in run_haplotype_coverage.columns: run_haplotype_coverage['amplicon'] = args.amplicon
    if not args.no_filter_cluster:
        filter_summary = [compute_summary(run_haplotype_coverage, 'initial')]
        # apply filters first 
        run_haplotype_coverage = (run_haplotype_coverage[run_haplotype_coverage
                                .apply(lambda row: row[['pass_coverage','pass_freq','pass_len']].sum() == 3, axis=1)])
        filter_summary = filter_summary + [compute_summary(run_haplotype_coverage, 'after_filtering')]
        # cluster
        run_haplotype_coverage = (run_haplotype_coverage
                                .groupby(['sample_id','haplotype_cluster_idx'])
                                .agg({'freq': np.sum, 
                                      'coverage': np.sum, 
                                      'len': lambda df: df.iloc[0], 
                                      'locus': lambda df: df.iloc[0], 
                                      'haplotype': lambda df: df.iloc[0],
                                      'amplicon': lambda df: df.iloc[0]})
                                .reset_index())
        del run_haplotype_coverage['haplotype_cluster_idx']
        filter_summary = filter_summary + [compute_summary(run_haplotype_coverage, 'after_clustering')]
    run_haplotype_coverage = add_metadata(run_haplotype_coverage, metadata, args.verbose)
    # compute haplotype index 
    if args.verbose: print("\tcomputing haplotypes index...")
    haplotype_index, run_haplotype_coverage_indexed = compute_haplotype_stats(run_haplotype_coverage, args.amplicon, min_population_freq=args.min_population_freq)
    # add amplicon as column 
    haplotype_index['amplicon'] = [args.amplicon] * len(haplotype_index.index)
    run_haplotype_coverage_indexed['amplicon'] = [args.amplicon] * len(run_haplotype_coverage_indexed.index)
    # save filtered, filtered/indexed and haplotype index data to files 
    filtered_filename = os.path.join(args.output_dir, '.'.join(filter(lambda s: s != '', [args.run_id,args.amplicon,args.filename_id,'filtered.tsv'])))
    run_haplotype_coverage.to_csv(filtered_filename, sep='\t', index=False)
    filtered_idxed_filename = os.path.join(args.output_dir, '.'.join(filter(lambda s: s != '', [args.run_id,args.amplicon,args.filename_id,'filtered.indexed.tsv'])))
    run_haplotype_coverage_indexed.to_csv(filtered_idxed_filename, sep='\t', index=False)
    haplotype_index_filename = os.path.join(args.output_dir, '.'.join(filter(lambda s: s != '', [args.run_id,args.amplicon,'haplotypes',args.filename_id,'index.tsv']))) 
    haplotype_index.to_csv(haplotype_index_filename, sep='\t', index=False)
    # if output_fastas, save to haplotypes to fasta
    if args.output_fastas: 
        print_to_fasta(haplotype_index, args.output_dir, args.amplicon, args.filename_id)
    # if compute run stats, do so 
    if args.compute_run_stats: 
        if args.verbose: print("\tcomputing run haplotype stats...")
        run_haplotype_stats, run_haplotype_coverage_indexed, haplotype_index = \
            compute_run_haplotype_stats(run_haplotype_coverage_indexed, haplotype_index, metadata, 
                                        args.amplicon, known_haplotypes)
        run_haplotype_stats.to_csv(args.output_dir + '/' + '.'.join(filter(lambda s: s != '', ['run',args.amplicon,args.filename_id,'stats.tsv'])), sep='\t')
    if not args.no_filter_cluster: 
        pd.DataFrame(filter_summary).to_csv(os.path.join(args.output_dir, '%s.%s.filter_summary.tsv' % (args.run_id, args.amplicon)), sep='\t', index=False)
    else: 
        pd.DataFrame([compute_summary(run_haplotype_coverage, 'initial')]).to_csv(os.path.join(args.output_dir, '%s.filter.cluster.summary.tsv' % args.amplicon), 
                                                                                sep='\t', index=False)
    return

if __name__ == '__main__': 
    main()

