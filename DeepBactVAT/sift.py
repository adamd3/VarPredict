import subprocess
import tempfile
import os
import sys
import re
import collections
from itertools import chain
import numpy as np
import pandas as pd
from pyfaidx import Fasta
import math


## TODO:
##  - accept GFF format input (instead of / as well as) TSV
##  - rename column containing the reference aa as 'ref'

def make_sift_files(
    vars_file, # file containing table of variants in protein-coding regions
    ann_file, # genome annotation file
    fasta_file, # amino acid multi-fasta file
    counts_file, # normalised read counts file
    refname,
    outdir): # output directory

    # read input
    count_data = pd.read_csv(counts_file, sep = "\t")
    vars_data = pd.read_csv(vars_file, sep = "\t", low_memory=False)
    annotation = pd.read_csv(ann_file, sep = "\t")
    gene_seq = Fasta(fasta_file)

    # subset to strains present in count data
    strain_names = count_data.columns.values.tolist()
    keep_cols = ['pos', 'ref'] + strain_names
    vars_data = vars_data[keep_cols]

    # remove variants where all strains identical to ref:
    vars_data_alleles = vars_data[['ref'] + strain_names]
    rm_rows = vars_data_alleles.eq(vars_data_alleles.iloc[:, 0], axis=0).all(1)
    vars_data = vars_data[~rm_rows]

    # extract gene name + position
    pos_split = vars_data['pos'].str.split("_")
    var_pos = [int(p[0]) for p in pos_split]
    var_genes = [p[1] for p in pos_split]
    vars_data['pos'] = var_pos
    vars_data['gene'] = var_genes

    # prepare annotation data
    annot_sub = annotation[['Locus Tag', 'Feature Type', 'Start', 'End', 'Strand', 'Gene Name']]
    cds_annot = annot_sub[annot_sub['Feature Type']=='CDS']
    cds_strands = cds_annot['Strand'].tolist()
    vars_sub = vars_data.copy(deep=True)

    # get codon/amino acid positions in the reference strain
    codon_pos_list = []
    for index, row in vars_sub.iterrows():
        pos = int(row['pos'])
        gene = row['gene']
        ref_aa = row['ref']
        if gene in cds_annot['Locus Tag'].tolist():
            annot = (cds_annot[cds_annot['Locus Tag']==gene])
            gene_start = int((cds_annot[cds_annot['Locus Tag']==gene])['Start'])
            gene_stop = int((cds_annot[cds_annot['Locus Tag']==gene])['End'])
            gene_strand = str(cds_annot.loc[cds_annot['Locus Tag'] == gene, 'Strand'].item())
            gene_len = (gene_stop-gene_start)+1
            if gene in gene_seq.keys():
                if gene_strand == "+":
                    codon_pos = math.floor((pos-gene_start)/3) ## NB this is 0-based indexing!
                    actual_codon_pos = codon_pos+1
                    if codon_pos < len(gene_seq[gene])+1:
                        try:
                            ref_seq_aa = str(gene_seq[gene])[codon_pos]
                        except IndexError:
                            ref_seq_aa = '*' ## stop codon
                        if ref_seq_aa==ref_aa:
                            codon_pos_list.append(actual_codon_pos)
                        else:
                            codon_pos_list.append(np.nan)
                    else:
                        codon_pos_list.append(np.nan)
                elif gene_strand == "-":
                    codon_pos = math.floor((gene_stop-pos)/3) ## NB this is 0-based indexing!
                    actual_codon_pos = codon_pos+1
                    if codon_pos < len(gene_seq[gene])+1:
                        try:
                            ref_seq_aa = str(gene_seq[gene])[codon_pos]
                        except IndexError:
                            ref_seq_aa = '*' ## stop codon
                        if ref_seq_aa==ref_aa:
                            codon_pos_list.append(actual_codon_pos)
                        else:
                            codon_pos_list.append(np.nan)
                    else:
                        codon_pos_list.append(np.nan)
        else:
            codon_pos_list.append(np.nan)
    vars_sub['codon_pos'] = codon_pos_list

    ## remove variants where the given reference AA didn't match the sequence
    vars_sub = vars_sub[vars_sub['codon_pos'].notna()]
    vars_sub['codon_pos'] = pd.to_numeric(vars_sub['codon_pos'], downcast='integer')
    vars_sub['codon_pos'] = vars_sub['codon_pos'].astype(str)

    # lists of strains and genes
    colnames = vars_sub.columns.values.tolist()
    rm_colnames = ['pos', 'ref', 'gene', 'codon_pos']
    strain_names = [x for x in colnames if x not in rm_colnames]
    gene_names =  vars_sub['gene'].tolist()

    # dict of dicts containing variants per strain, per gene
    strain_var_dict = {}
    for gene_name in gene_names:
        gene_vars = vars_sub[vars_sub['gene']==gene_name]
        strain_var_dict[gene_name] = {}
        for strain in strain_names:
            (strain_var_dict[gene_name])[strain] = []
            strain_vars = gene_vars[['ref', 'codon_pos', strain, 'gene']]
            strain_vars = strain_vars[~(strain_vars['ref']==strain_vars[strain])]
            vars_list = (strain_vars['ref'] + strain_vars['codon_pos'] + strain_vars[strain]).tolist()
            if vars_list:
                (strain_var_dict[gene_name][strain]).append(vars_list)

    # save input files for SIFT
    for gene_name in strain_var_dict.keys():
        all_vars = list(chain.from_iterable(strain_var_dict[gene_name].values()))
        all_vars = list(set([x for xs in all_vars for x in xs]))
        outfile = os.path.join(outdir, 'subst_files', gene_name+'.subst')
        with open(outfile, 'a+') as f:
            for var in all_vars:
                f.write("%s\n" % var)


def run_sift(
    query, # query sequences
    subst, # aa substitutions directory
    database, # database to search
    outdir, # output directory
    kmer=5, # length of kmers used for database search (default = 5)
    maxcand=5000, # number of database sequences retained for local alignment (default = 5000)
    threads=32, # number of threads to use
    evalue=0.0001, # e-value threshold (default = 0.0001)
    quiet=False):

    cmd = "sift4g"
    cmd += " --query " + query
    cmd += " --subst " + subst
    cmd += " --database " + database
    cmd += " --out " + outdir
    cmd += " --kmer-length " + str(kmer)
    cmd += " --max-candidates " + str(maxcand)
    cmd += " --threads " + str(threads)
    cmd += " --evalue " + str(evalue)

    if not quiet:
        print("executing command: " + cmd)
    else:
        cmd += " > /dev/null"

    subprocess.run(cmd, shell=True, check=True)

    return
