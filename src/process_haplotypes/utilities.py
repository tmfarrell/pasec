#!/usr/bin/env python

import re
import sys
import csv
import random
import Levenshtein
import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.linear_model import OrthogonalMatchingPursuit
bases = ['A','T','C','G']

## utils  
def safe_div(num, denom): 
    if denom == 0: return(0) 
    else: return(float(num)/float(denom))

def to_str(X): 
    return('/'.join(map(str, X)))

def from_str(X_str): 
    return(tuple(X_str.split('/')))

def to_str2(X, Y): 
    if X and Y and type(X) == tuple and type(Y) == tuple:  
        return('/'.join(map(lambda pair: '%s:%0.4f' % pair, zip(X, Y))))
    return('')

def dict_to_str(D): 
    return('/'.join(['%s:%s' % item for item in D.items()]))

def to_dict(str_): 
    try: return(dict([s.split(':') for s in str_.split('/')]))
    except: return({})

def replace_keys(dict_, key_to_new_key): 
    new_dict = {}
    for key in dict_.keys():
        new_dict[key_to_new_key[key]] = dict_[key]
    return(new_dict)

def metrics_to_str(metrics): 
    return('/'.join(['%s:%0.4f' % (str(k), float(v)) for k, v in metrics.items()]))

def metrics_to_dict(metrics_str): 
    try: return({k: float(v) for k, v in [s.split(':') for s in metrics_str.split('/')]})
    except: return({})

def filter_metrics_dict(metrics, lower_bound): 
    return({m: v for m, v in metrics.items() if v > lower_bound})

def cts_dict_to_freqs(cts): 
    total = float(sum(cts.values()))
    return({k: float(v)/total for k, v in cts.items()})

def metrics_dicts_mean(metrics_dicts): 
    total = {k: 0 for k in set([key for keys in map(lambda d: d.keys(), metrics_dicts) for key in keys])}
    num_dicts = len(metrics_dicts)
    for i in range(num_dicts): 
        d = metrics_dicts[i]
        for k in d.keys(): 
            total[k] = total[k] + d[k]
    return({k: float(v)/float(num_dicts) for k, v in total.items()})

def cts_to_freqs(position_counts_df): 
    return(position_counts_df.apply(lambda row: 
                                    [save_div(row[op], row.sum()) for op in position_counts_df.columns], 
                                    axis=1))

def print_to_fasta(haplotypes_df, output_dir, amplicon, filename_id, remove_deletion_char=True):
    f = open(output_dir + '/' + amplicon + '.' + (filename_id + '.' if filename_id else '') + 'index.fasta', 'w')
    for i in haplotypes_df.index:
        row = haplotypes_df.loc[i,:].to_dict()
        f.write('-'.join([">Pf3D7", row['amplicon'], (row['run'] if 'run' in row.keys() else ''), row['haplotype_index'], str(row['coverage'])]) + '\n')
        f.write(row['haplotype'].replace('-','' if remove_deletion_char else '-') + '\n\n')
    f.close()
    return
    

## false negative haplotype analysis utils
def lookup(df, queries, field=None):
    for query_field, query_val in queries: 
        df = df[df[query_field] == query_val]
    try:    
        if field: return(df.iloc[0][field])
        else:     return(df.iloc[0])
    except: return('')

def edit_distance(r1_index, r2_index, amplicon, haplotypes_index_df): 
    r1 = lookup(haplotypes_index_df, [('amplicon',amplicon), ('haplotype_index',r1_index)], 'haplotype')
    r2 = lookup(haplotypes_index_df, [('amplicon',amplicon), ('haplotype_index',r2_index)], 'haplotype')
    return(Levenshtein.distance(r1, r2))

def edit_ops(r1_index, r2_index, amplicon, haplotypes_index_df): 
    if r1_index != 'None' and r2_index != 'None':
        def edit_op_to_str(edit_op, s1, s2): 
            op, p1, p2 = edit_op
            if op == 'replace': return(str(p1) + s1[p1] + '->' + s2[p2])
            elif op == 'insert': return(str(p1) + 'I' + '->' + s2[p2])
            elif op == 'delete': return(str(p1) + '->D')
        r1 = lookup(haplotypes_index_df, [('amplicon',amplicon), ('haplotype_index',r1_index)], 'haplotype')
        r2 = lookup(haplotypes_index_df, [('amplicon',amplicon), ('haplotype_index',r2_index)], 'haplotype')
        if abs(len(r1) - len(r2)) < 20: 
            editops = Levenshtein.editops(r1, r2)
            return(map(lambda edit_op: edit_op_to_str(edit_op, r1, r2), editops))
        else: 
            return([])
    else: 
        return([])

def to_mismatch_cts(coverage, edit_pos_list, length=350): 
    cts = [0] * length
    for edit_pos in edit_pos_list: 
        try: cts[int(re.match('[1-9][0-9]*', edit_pos).group())] = coverage
        except: pass
    return(cts)

def dists_to_fn(strain_to_read, read_to_assign, amplicon, haplotypes_index_df, reduce_fcn=np.mean): 
    strain_to_read = to_dict(strain_to_read)
    read_to_assign = to_dict(read_to_assign)
    false_negatives = [r for r in strain_to_read.values() if not r in read_to_assign.keys()]
    return({fn: reduce_fcn([edit_distance(ri, fn, amplicon, haplotypes_index_df) for ri in read_to_assign.keys()
                             if not ri in strain_to_read.values()])
            for fn in false_negatives})

def edit_pos_to_fn(strain_to_read, read_to_assign, amplicon, read_coverages, haplotypes_index_df): 
    strain_to_read = to_dict(strain_to_read)
    #print('strains: ' + str(strain_to_read))
    read_to_assign = to_dict(read_to_assign)
    #print('reads: ' + str(read_to_assign))
    false_negatives = [r for r in strain_to_read.values() if not r in read_to_assign.keys()]
    edit_pos_to_fns = [(fn, [strain for strain in strain_to_read.keys() if strain_to_read[strain] == fn][0], 
                        [(ri, read_to_assign[ri], read_coverages[ri], edit_ops(fn, ri, amplicon, haplotypes_index_df)) 
                             for ri in read_to_assign.keys()
                             if not ri in strain_to_read.values()])
                       for fn in false_negatives]
    return(edit_pos_to_fns)

# for quantifying strain assignment performance
def f1_score(truth_set, prediction_set): 
    recall = recall_multilabel(truth_set, prediction_set)
    precision = precision_multilabel(truth_set, prediction_set)
    return(safe_div(2.0, (safe_div(1.0, precision) + safe_div(1.0, recall))))

def fdr(truth_set, prediction_set): 
    return(safe_div(len([p for p in prediction_set if not p in truth_set]), len(prediction_set)))

def jaccard_index(set1, set2): 
    # len(intersection(set1, set2)) / len(union(set1, set2))
    return(safe_div(float(len(set(set1).intersection(set2))), float(len(set(set1).union(set2)))))

def precision_multilabel(truth_set, prediction_set): 
    return(safe_div(float(len(set(prediction_set).intersection(truth_set))), float(len(set(prediction_set)))))

def recall_multilabel(truth_set, prediction_set): 
    return(safe_div(float(len(set(prediction_set).intersection(truth_set))), float(len(set(truth_set)))))

## cigar/ mismatch_str ops
def to_cigar_arr(cigar_str): 
    if "/" in cigar_str: 
        # just take first one
        cigar_str = cigar_str.split("/")[0]
    return(map(lambda match: (int(re.search('[0-9]+', match).group(0)), re.search('[MIDNSHPX=]', match).group(0)), 
               re.findall('[0-9]+[MIDNSHPX=]', cigar_str)))

def to_cigar_str(cigar_arr): 
    return(''.join(map(lambda arr: ''.join(map(str, arr)), cigar_arr)))

def cigar_arr_to_long_cigar(cigar_arr): 
    return(''.join(map(lambda count_op: count_op[0] * count_op[1], cigar_arr)))

def mismatch_str_to_long_mismatch(mismatch_str): 
    if "/" in mismatch_str: 
        # just take first one
        mismatch_str = mismatch_str.split("/")[0]
    long_mismatch = ''
    ops = [s for s in re.split('([0-9]+|\^?[ATCG])', mismatch_str) if s != '']
    for op in ops: 
        if '^' in op: #skip deletions 
            continue
        elif not op in ['A','T','C','G']: # then op is # of ref matches
            # insert Rs for that # of ref matches 
            long_mismatch = long_mismatch + (int(op) * 'R')
        else:                             # insert X for ref  
            long_mismatch = long_mismatch + 'X'
    return(long_mismatch)

def mismatch_str_to_cts(mismatch_str): 
    long_mismatch = mismatch_str_to_long_mismatch(mismatch_str)
    return({op: long_mismatch.count(op) for op in ['X','R']})

def long_cigar_to_cigar_arr(long_cigar): 
    arr = []
    count = 0
    curr = long_cigar[0]
    for i in range(len(long_cigar)): 
        if curr == long_cigar[i]: 
            count = count + 1
        else: 
            arr = arr + [[count, curr]]
            curr = long_cigar[i]
            count = 1
        if i == (len(long_cigar) - 1):  # if last 
            arr = arr + [[count, curr]]
    return(arr)

def get_snps(mismatch_str, read, pos_or_base='pos'): 
    #sys.stderr.write(read)
    pos = 0
    snps = []
    for token in re.split('(\^?[ATCG]+)', mismatch_str): 
        if token in bases: 
            pos = pos + 1
            # fix this later
            if pos_or_base == '': 
                snps = snps + [(pos, read[pos])]
            else: # read_or_base == 'pos': 
                snps = snps + [(pos, '')]
        elif '^' in token: 
            pos = pos + len(token) - 1
        else: 
            pos = pos + int(token)
    if not snps:
        return(';')
    elif pos_or_base == 'both':
        return(';'.join(map(lambda snp: ''.join(map(str, snp)), snps)))
    else: #pos_or_base == 'pos': 
        return(';'.join(map(lambda snp: snp[0], snps)))
    #elif pos_or_base == 'base': 
    #    return(';'.join(map(lambda (pos, base): str(base), snps)))

def get_insertion_positions(cigar_arr):
    positions = []
    for i in range(len(cigar_arr)):  
        count, op = cigar_arr[i]
        if op == "I":  
            curr_pos = sum(map(lambda count_op: count_op[0], filter(lambda count_op: count_op[1] != "I", cigar_arr[:i])))
            positions = positions + [(curr_pos, count)]
    return(positions)

def get_insertion_position_set(read_stats_cigars, verbose=False):
    if not type(read_stats_cigars) == list: 
        read_stats_cigars = read_stats_cigars.values.tolist()
    insert_pos = list(set([v for L in map(lambda c: get_insertion_positions(to_cigar_arr(c)), read_stats_cigars)
                           for v in L]))
    # make unique, taking the highest count at each position
    insert_pos = [(v, max(map(lambda t: t[1], filter(lambda t: t[0] == v, insert_pos)))) 
                  for v in list(set(map(lambda t: t[0], insert_pos)))]
    insert_pos.sort()
    return(insert_pos)

def get_insertion_pos_df(read_stats_file):
    rs_stats = pd.read_csv(read_stats_file, sep='\t')
    insertion_positions = get_insertion_position_set(rs_stats.cigar.values.tolist())
    return(pd.DataFrame(map(lambda pos_num: 
                            {'pos': pos_num[0], 'num_inserts': pos_num[1]}, insertion_positions),
                        columns=['pos', 'num_inserts']))

def parse_insertion_position_set_file(insertion_pos_file, amplicon=None):
    insertions = pd.read_csv(insertion_pos_file, sep='\t')
    if insertions.empty: 
        return([])
    elif amplicon: 
        return(insertions[insertions['amplicon'] == amplicon][['pos','num_inserts']].values.tolist())
    else: 
        return(insertions[['pos','num_inserts']].values.tolist())
 
def to_padded_cigar_arr(cigar_arr, insertion_positions, verbose=False): 
    long_cigar = cigar_arr_to_long_cigar(cigar_arr)
    return(long_cigar_to_cigar_arr(to_padded_long_cigar(long_cigar, insertion_positions, verbose)))

def to_padded_cigar(cigar, insertion_positions): 
    return(to_cigar_str(to_padded_cigar_arr(to_cigar_arr(cigar), insertion_positions)))

def to_padded_long_cigar(long_cigar, insertion_positions, start_pos_delta, verbose=False):    
    if not insertion_positions: return(long_cigar)
    num_soft_clip = len(re.match('S*', long_cigar).group())
    padded = 'D' * (start_pos_delta - num_soft_clip)
    insertion_pos_idx, unpadded_idx, padded_idx = 0, 0, (-1)*num_soft_clip
    while True: 
        curr_op = long_cigar[unpadded_idx]
        if verbose: 
            print("curr_op: %s, unpadded_idx: %d, padded_idx: %s, insertion_pos_idx: %d, padded: %s" % \
                  (curr_op, unpadded_idx, padded_idx, insertion_pos_idx, padded)) 
        # if curr_pos (in padded space) is an insertion position (in unpadded space)
        if padded_idx == (insertion_positions[insertion_pos_idx][0] - start_pos_delta):
            num_insertions = insertion_positions[insertion_pos_idx][1]
            # if base at curr_pos is not insert
            if curr_op != "I": 
                # add P for padding 
                padded = padded + ("P" * num_insertions)
                insertion_pos_idx = insertion_pos_idx + 1
            else: # if curr_pos is an insert
                insert_section = long_cigar[unpadded_idx:(unpadded_idx + num_insertions)]
                if verbose: print('insert_section: %s' % (insert_section))
                insert_section_arr = long_cigar_to_cigar_arr(insert_section)
                if verbose: print('insert_section_arr: %s' % (str(insert_section_arr)))
                n_I = insert_section_arr[0][0]  # since first cigar op must be I
                n_P = num_insertions - n_I
                # add to padded 
                padded = padded + ("I" * n_I) + ("P" * n_P)
                # incr unpadded index  
                unpadded_idx = unpadded_idx + 1
                # move to next insertion position
                insertion_pos_idx = insertion_pos_idx + 1
        # otherwise just add current op to padded cigar
        else: 
            if curr_op != 'I': 
                padded = padded + curr_op
                padded_idx = padded_idx + 1
                unpadded_idx = unpadded_idx + 1
            else: 
                unpadded_idx = unpadded_idx + 1
        # then break if either unpadded index greater than long cigar
        if (unpadded_idx - padded.count("I")) >= len(long_cigar): 
            break
        # or if ran through all the insertion positions
        if insertion_pos_idx >= len(insertion_positions): 
            padded = padded + long_cigar[unpadded_idx:]
            break
    return(padded)

def get_locus_ref_seq(ref_fasta_path, locus): 
    chrom, start, end = locus
    ref = SeqIO.parse(ref_fasta_path, 'fasta')
    for record in ref: 
        if record.id == chrom: 
            return(record.seq[(start-1):(end-1)])
    return('')

def get_aligned_seq(seq, cigar_arr, start_pos_delta, trim_to=None, 
                    remove_T_homopolymers_of_len=4, verbose=False): 
    assert trim_to in [None, 'locus', 'locus_start'],\
           "trim_to parameter must be nil, 'locus' or 'locus_start'."
    aligned_seq = ''
    if verbose:  
        print("")
        print('Getting aligned sequence for:') 
        print('cigar: ' + str(cigar_arr))
        print("seq: #{seq.length} #{seq}")
    # build aligned seq
    for count, cigar_op in cigar_arr: 
        if verbose: print("#{[count, cigar_op]}") 
        if cigar_op in ['M','N','I','S']: 
            aligned_seq += seq[:count]
            seq = seq[count:] 
        elif cigar_op == 'P': 
            aligned_seq += ''.join(['*'] * count)
        elif cigar_op == 'D':     # if deletion, insert "-" 
            aligned_seq += ''.join(['-'] * count)
        elif cigar_op in ['H']: 
            seq = seq[count:]
        if verbose: print("aligned_seq so far: #{aligned_seq.length} #{aligned_seq}")
    # trim it, if applicable
    if trim_to: 
        chrom, coords, id = interval  
        if verbose: print("aligned seq before trimming: #{aligned_seq} #{aligned_seq.length}") 
        if trim_to == 'locus':
            aligned_seq = aligned_seq[(coords.first - aligned_start_pos):(coords.last - aligned_start_pos)]
        elif trim_to == 'locus_start':
            aligned_seq = aligned_seq[(coords.first - aligned_start_pos):]
        if verbose: print("aligned seq after trimming: #{aligned_seq} #{aligned_seq.length}")
    else:
        if verbose: print("aligned seq: #{aligned_seq} #{aligned_seq.length}") 
    if remove_T_homopolymers_of_len: 
        aligned_seq = ''.join(re.split('[T]{%s,}' % str(remove_T_homopolymers_of_len), aligned_seq))
    return(aligned_seq)

## strain fcns
def encode_strain_label(strain_label, strains_list): 
    return(map(lambda s: 1 if s in strain_label else 0, strains_list))

def encode_strain_mix(row, strains_list): 
    if pd.notnull(row['Mixture']) and row['Mixture'] != '': 
        strain_mix = dict(zip(row['Strains'].split(';'), map(lambda x: float(x) / 100.0, row['Mixture'].split(':')))) 
    else: 
        strain_mix = {}
    return(map(lambda s: strain_mix[s] if s in strain_mix else 0, strains_list))

def to_full_haplotype_arr(read_freqs, num_haplotypes, read_to_idx_map, 
                          haplotype_mock_coverages=None, coverage=None):
    arr = [0] * num_haplotypes
    read_freqs = metrics_to_dict(read_freqs)
    for read in [r for r in read_freqs.keys() if r in read_to_idx_map.keys()]: 
        if not haplotype_mock_coverages: 
            if coverage: 
                arr[read_to_idx_map[read]] = read_freqs[read] * coverage
            else: 
                arr[read_to_idx_map[read]] = read_freqs[read] 
        else: 
            arr[read_to_idx_map[read]] = read_freqs[read] * haplotype_mock_coverages[read]
    return(arr)

def assign_strains_by_label(haplotype_index, solexa_read_stats): 
    try: 
        # get solexa read stats mock data 
        solexa_read_stats = solexa_read_stats[(solexa_read_stats.contents == 'p_falciparum_mock')]
        # get haplotypes that occur in mock samples
        mock_haplotype_index = haplotype_index[(haplotype_index.in_mock == 1)]
    except: 
        mock_haplotype_index = haplotype_index.copy()
    # let N be # solexas
    N = len(solexa_read_stats.solexa_index.unique())
    # let M be # haplotypes unique for that amplicon
    M = len(mock_haplotype_index.haplotype_index.unique())
    h_idx_to_m = dict(zip(mock_haplotype_index.haplotype_index.values.tolist(),
                          range(len(mock_haplotype_index.index))))
    try: h_idx_to_mock_coverage = dict(mock_haplotype_index[['haplotype_index','mock_coverage']].values.tolist())
    except: h_idx_to_mock_coverage = dict(mock_haplotype_index[['haplotype_index','coverage']].values.tolist())
    m_idx_to_h = dict(zip(range(len(mock_haplotype_index.index)),
                          mock_haplotype_index.haplotype_index.values.tolist()))
    # let S be list of strains, and K be # strains
    S = sorted(list(set([strain for strains in (solexa_read_stats.Strains.apply(lambda s: s.split(';')).values.tolist()) 
                         for strain in strains])))
    K = len(S)
    # build X, where x_nm in [0...inf] is coverage of haplotype m for solexa n
    #          for n in 0...N-1, m in 0...M-1
    X = pd.DataFrame(solexa_read_stats.apply(lambda row: to_full_haplotype_arr(row['haplotype_freqs'], M, h_idx_to_m), axis=1)
                     .values.tolist()).as_matrix()
    # build Y, where y_nk in {0,1} is an indicator of strain k in solexa n
    #          for n in 0...N-1, k in 0...K-1
    Y = pd.DataFrame(solexa_read_stats
                     .apply(lambda row: encode_strain_mix(row, S), axis=1)
                     .values.tolist()).as_matrix()
    # learn map from M' -> K, where M' subset of M
    m_to_k_map = OrthogonalMatchingPursuit(n_nonzero_coefs=2)
    m_to_k_map.fit(X, Y) 
    strain_label_assign = [(strain, [(m_idx_to_h[i], c) for i, c in [(i, x) for i, x in enumerate(coefs) if x > 0]])
                           for strain, coefs in zip(S, m_to_k_map.coef_)]
    # label those haplotypes w/ indices in M', as associated strain mapped to K
    for strain, haplotype_indices in strain_label_assign: 
        for hi, coef in haplotype_indices: 
            haplotype_index.loc[(haplotype_index.haplotype_index == hi),'label'] = strain
            haplotype_index.loc[(haplotype_index.haplotype_index == hi),'label_coef'] = coef
    haplotype_index['label'] = haplotype_index['label'].fillna('UNK')
    haplotype_index['label_coef'] = haplotype_index['label_coef'].fillna(0.0)
    return(haplotype_index)

def assign_strains(hs_df, solexa_read_stats, known_haplotypes=None): 
    # method to break ties
    def break_tie(best_two): 
        vals = map(lambda x: x[0], best_two) 
        if vals[0] == vals[1]: return(best_two[0][1] if random.random() > 0.5 else best_two[1][1])
        else:                  return(best_two[0][1])
    if not known_haplotypes is None: 
        # for strain in known haplotypes, compute distance to that strain for each haplotype  
        for strain in known_haplotypes.strain.unique(): 
            hs_df['dist_to_' + strain] = \
                (hs_df.apply(lambda row: Levenshtein.distance(row['haplotype'], 
                                                              known_haplotypes[(known_haplotypes.strain == strain)].iloc[0].haplotype), 
                             axis=1))
        # assign each haplotype its closest strain match  
        hs_df['closest'] = (hs_df.loc[:, filter(lambda c: 'dist_to' in c, hs_df.columns)]
                            .apply(lambda row: break_tie(zip(row.nsmallest(2).tolist(), 
                                                             row.nsmallest(2).index.tolist())), axis=1)
                            .apply(lambda s: s[s.rfind('_')+1:]))
        # assign one haplotype as best match to strain, taking the closest match w/ highest mock coverage
        hs_df['best'] = len(hs_df.index) * ["UNK"]
        for strain in hs_df.closest.unique(): 
            try: hs_df.loc[hs_df[hs_df.closest == strain].mock_coverage.idxmax(),'best'] = strain
            except: hs_df.loc[hs_df[hs_df.closest == strain].coverage.idxmax(),'best'] = strain
    # assign strains using strain composition labels 
    try: 
        hs_df = assign_strains_by_label(hs_df, solexa_read_stats)
        strain_assign_fields = ['label']
    except: 
        strain_assign_fields = [] 
    if not known_haplotypes is None: 
        strain_assign_fields = strain_assign_fields + ['best','closest']
    # get consensus strain assignment
    hs_df['consensus'] = hs_df.apply(lambda row: row[strain_assign_fields].value_counts().idxmax() 
                                                 if row[strain_assign_fields].value_counts().max() > 1 else 'UNK', axis=1)
    try: 
        # for those labels that are not in known_haplotypes, make that label the consensus
        for strain in [s for s in hs_df.label.unique() if not s in list(known_haplotypes.strain.unique()) and s != 'UNK']: 
            idx = hs_df[hs_df.label == strain].index
            hs_df.loc[idx, 'consensus'] = [strain] * len(idx)
    except: pass
    strain_assign_fields = strain_assign_fields + ['consensus']
    return(hs_df, strain_assign_fields)

def to_strain_cts(sample_coverage_cts_str, sample_metadata, count_coverage=True, scale_by_num_strains=False): 
    strain_cts = {}
    sample_coverage_cts = metrics_to_dict(sample_coverage_cts_str)
    is_mock_sample = lambda s: ('batch' not in s) and \
                               (not s in ['patient_blood_spot','blank','Human','plasmid','Water'])
    for sample in sample_coverage_cts.keys(): 
        strains = filter(is_mock_sample, (sample_metadata[sample_metadata['sample_index'] == sample]['Strains']
                                          .apply(lambda s: re.split(';', s)).values.tolist()[0]))
        if strains:
            scale_factor = (1.0 / float(len(strains))) if scale_by_num_strains else 1.0
            for strain in strains: 
                if not strain in strain_cts.keys(): 
                    if count_coverage: strain_cts[strain] = scale_factor * sample_coverage_cts[sample]
                    else: strain_cts[strain] = scale_factor
                else: 
                    if count_coverage: strain_cts[strain] = strain_cts[strain] + (scale_factor * sample_coverage_cts[sample])
                    else: strain_cts[strain] = strain_cts[strain] + scale_factor
    return(strain_cts)

## allele fcns
def add_ref_pos_index(alleles, ref_file, chrom, min_start_pos, max_len, insertion_positions): 
    ref, pos = [], []
    i, j, insert_pos_idx = 0, 0, 0
    ref_seq = get_locus_ref_seq(ref_file, (chrom, min_start_pos, (min_start_pos + max_len)))
    while i < max_len: 
        if insert_pos_idx < len(insertion_positions) and \
           j == insertion_positions[insert_pos_idx][0]: 
            num_insertions = insertion_positions[insert_pos_idx][1]
            ref = ref + (["*"] * num_insertions)
            pos = pos + ([min_start_pos + j] * num_insertions)
            i = i + num_insertions
            insert_pos_idx = insert_pos_idx + 1
        else:
            ref = ref + [ref_seq[j]]
            pos = pos + [min_start_pos + j]
            i = i + 1
            j = j + 1
    alleles['pos'] = pos
    alleles['ref'] = ref
    return(alleles)

def get_allele_cts(rs, amplicon, count_coverage=True, insertion_pos_file=None, ref_file=None): 
    read_field = 'read' if (not 'padded_read' in rs.columns) else 'padded_read' 
    max_len = rs[read_field].apply(len).max()
    #print('max_len: ' + str(max_len))
    if not 'start_pos' in rs.columns: 
        rs['start_pos'] = (rs.locus.apply(
                            lambda l: int(l.split(':')[1]) 
                            if not '/' in l 
                            else min([int(s.split(':')[1]) for s in l.split('/')])))
    min_start_pos = rs['start_pos'].min()
    #print('min_start_pos: ' + str(min_start_pos))
    allele_counts = []
    for i in range(max_len): 
        allele_counts = allele_counts + [{'A':0,'T':0,'C':0,'G':0,'*':0,'-':0,'N':0}]
    for i in rs.index:
        read = rs[read_field][i]
        #print('read start_pos: ' + str(rs['start_pos'][i]))
        delta = rs['start_pos'][i] - min_start_pos
        #print('delta: ' + str(delta))
        for j in range(len(read)): 
            #print j+delta
            if count_coverage: 
                allele_counts[j+delta][read[j]] += rs.loc[i, 'coverage']
            else: 
                allele_counts[j+delta][read[j]] += 1
    allele_cts = pd.DataFrame(allele_counts)
    if insertion_pos_file: 
        chrom = rs.loc[:,'locus'].iloc[0].split(':')[0]
        allele_cts = add_ref_pos_index(allele_cts, ref_file, chrom, min_start_pos, max_len,
                                       parse_insertion_position_set_file(insertion_pos_file, amplicon))
    return(allele_cts)
    
def get_allele_freqs(rs, count_coverage=True, padded=True, insertion_pos_file=None): 
    min_start_pos = rs['start_pos'].min()
    read_field = 'read' if (not padded) else 'padded_read' 
    max_len = rs[read_field].apply(len).max()
    allele_freqs = cts_to_freqs(get_allele_cts(rs, count_coverage, padded))
    if insertion_pos_file: 
        allele_freqs = add_ref_pos_index(allele_freqs, min_start_pos, max_len,
                                         parse_insertion_position_set_file(insertion_pos_file))
    return(allele_freqs)

## indel/ mismatch computation fcns
def alleles_to_mismatch(alleles): 
    return(pd.DataFrame(alleles.apply(lambda row: 
                                      {'X': row[['*','-','A','T','C','G']].sum() - row[row['ref']]}, 
                                      axis=1).values.tolist()))

def get_indel_cts(read_stats, amplicon=None):
    if amplicon: 
        read_stats = read_stats[read_stats['amplicon'] == amplicon]
    # get df of unique cigars, with desired metric
    cigar_coverage = read_stats.groupby('cigar').coverage.sum()
    # get max len of cigars
    long_cigars = map(lambda c: cigar_arr_to_long_cigar(to_cigar_arr(c)), cigar_coverage.index)
    max_len = max(map(len, long_cigars))
    # init empty counts struct, with len max_len
    position_cigar_op_cts = []
    for i in range(max_len): 
        position_cigar_op_cts = position_cigar_op_counts + [{'M':0,'D':0,'I':0,'S':0,'H':0}]
    # for each cigar, add cigar op coverage counts to counts at each position
    for cigar in cigar_coverage.index: 
        long_cigar = cigar_arr_to_long_cigar(to_cigar_arr(cigar))
        for i in range(len(long_cigar)): 
            position_cigar_op_cts[i][long_cigar[i]] += cigar_coverage[cigar]
    return(pd.DataFrame(position_cigar_op_cts)[['D','I']])

def get_mismatch_cts(read_stats, amplicon=None):
    if amplicon: 
        read_stats = read_stats[read_stats['amplicon'] == amplicon]
    # basically same as get_indel_cts
    mismatch_coverage = read_stats.groupby('mismatch_str').coverage.sum().reset_index()
    max_len = mismatch_coverage.mismatch_str.apply(lambda s: len(mismatch_str_to_long_mismatch(s))).max()
    position_mismatch_cts = []
    for i in range(max_len): 
        position_mismatch_cts = position_mismatch_cts + [{'R':0,'X':0}]
    for i in mismatch_coverage.index:
        long_mismatch = mismatch_str_to_long_mismatch(mismatch_coverage.loc[i,'mismatch_str'])
        for j in range(len(long_mismatch)): 
            coverage = mismatch_coverage.loc[i,'coverage']
            mismatch_counts[i][long_mismatch[i]] += coverage
    return(pd.DataFrame(mismatch_counts))

def get_indel_mismatch_freqs(read_stats, amplicon=None): 
    indel_freqs = cts_to_freqs(get_indel_cts(read_stats, amplicon))
    mismatch_freqs = cts_to_freqs(get_mismatch_cts(read_stats, amplicon))
    return(pd.concat([indel_freqs[['I','D']], mismatch_freqs['X']], axis=1))


## file parsing 
def  parse_run_haplotype_coverage_file(path): 
    t = pd.read_csv(path, sep='\t')
    assert 'sample_id' in t.columns, "sample_id must be a column in the haplotype coverage file."
    return(t)

def get_sample_metadata(sample_metadata_path, verbose=False):
    # get metadata
    sample_metadata = pd.read_csv(sample_metadata_path)
    assert 'Sample_ID' in sample_metadata.columns, "Sample_ID needs to be a column in sample metadata"
    sample_metadata['sample_index'] = map(lambda i: 'S' + str(i), sample_metadata.index.values)
    try: 
        sample_metadata['malaria'] = sample_metadata['Parasite_Density'].apply(lambda d: 1 if d > 0 else 0)
        sample_metadata['strain_mix'] = sample_metadata.apply(lambda row: row['Strains'] + '-' + str(row['Mixture']), axis=1)
    except: pass
    if verbose:
        print("sample_metadata:\n" + sample_metadata.to_csv(sep='\t'))
    return(sample_metadata)
    
def parse_metadata_file(metadata_file, verbose=False): 
    metadata = pd.read_csv(metadata_file, sep='\t')
    assert 'sample_id' in metadata.columns, "sample_id needs to be a valid column in metadata."
    metadata['sample_index'] = ['S' + str(i) for i in metadata.index.values] 
    return(metadata)

def get_solexa_metadata(solexa_metadata_path, verbose=False): 
    solexa_metadata = pd.read_csv(solexa_metadata_path)
    assert 'Sample_ID' in solexa_metadata.columns, "Sample_ID needs to be a column in solexa metadata"
    assert 'Solexa_ID' in solexa_metadata.columns, "Solexa_ID needs to be a column in solexa metadata"
    # add index
    solexa_metadata['solexa_index'] = map(lambda i: 'X' + str(i), solexa_metadata.index.values)
    try: 
        # add plate # for well_pos plotting
        well_pos_cts = solexa_metadata.well_pos.value_counts()
        for well_pos in well_pos_cts.index: 
            solexa_metadata.loc[solexa_metadata.well_pos == well_pos, 'plate_num'] = range(1, well_pos_cts[well_pos] + 1)
    except: pass
    if verbose: 
        print("solexa_metadata:\n" + solexa_metadata.to_csv(sep='\t'))
    return(solexa_metadata)

def add_sample_measures(solexa_metadata, sample_metadata, verbose=False):
    if 'malaria' in solexa_metadata.columns: 
        del solexa_metadata['malaria']
    if 'Solexa_ID' in sample_metadata.columns: 
        sample_metadata = sample_metadata.drop(columns=['Solexa_ID'])
    #print("solexa_metadata:\n" + solexa_metadata.to_csv(sep='\t'))
    solexa_metadata = solexa_metadata.merge(sample_metadata, on='Sample_ID', how='left')
    try: 
        solexa_metadata['Parasite_Density'] = solexa_metadata.apply(lambda row: row['contents'].split('_')[0] \
                                                                  if (row['contents'] != 'p_falciparum_mock') \
                                                                  else row['Parasite_Density'], axis=1)
        solexa_metadata['num_strains'] = solexa_metadata.apply(lambda row: row['contents'].split('_')[0] \
                                                               if (row['contents'] != 'p_falciparum_mock') \
                                                               else row['Strains'].count(';') + 1, axis=1)
    except: pass
    if verbose: 
        print('solexa_metadata after adding sample relations:\n' + str(solexa_metadata.head()))
    return(solexa_metadata)

def add_metadata(df, metadata=None, verbose=False): 
    if verbose: 
        print(df.head())
        if metadata: print(metadata.head())
    if not metadata: 
        metadata = pd.DataFrame()
        metadata['sample_id'] = df.sample_id.unique()
        metadata['sample_index'] = ['S' + str(i) for i in metadata.index.values]
    df = df.merge(metadata, on='sample_id', how='left')
    if verbose: 
        print('after adding metadata:\n' + str(df.head()))
    return(df)

def add_false_negative_haplotype_fields(solexa_hs, hs): 
    solexa_hs['strain_to_haplotype'] = (solexa_hs
                                   .apply(lambda row: '/'.join([s + ':' + str(lookup(hs, [('amplicon', row['amplicon']), ('consensus', s)], 
                                                                                     'haplotype_index'))
                                                                for s in row['Strains'].split(';') if (s != '' and s != 'UNK')]), axis=1))
    solexa_hs['min_dist_to_false_neg'] = (solexa_hs.apply(lambda row: 
                                                      metrics_to_str(dists_to_fn(row['strain_to_haplotype'], 
                                                                                 row['haplotype_to_consensus'], 
                                                                                 row['amplicon'], hs, 
                                                                                 lambda ls: min(ls) if ls else np.nan)), 
                                                      axis=1))
    solexa_hs['mean_dist_to_false_neg'] = (solexa_hs.apply(lambda row: 
                                                       metrics_to_str(dists_to_fn(row['strain_to_haplotype'], 
                                                                                  row['haplotype_to_consensus'], 
                                                                                  row['amplicon'], hs)), 
                                                       axis=1))
    solexa_hs['edit_pos_freq_to_false_neg'] = (solexa_hs.apply(lambda row: 
                                                           edit_pos_to_fn(row['strain_to_haplotype'], 
                                                                           row['haplotype_to_consensus'], 
                                                                           row['amplicon'],
                                                                           metrics_to_dict(row['haplotype_freqs']), hs), 
                                                           axis=1))
    return(solexa_hs)
