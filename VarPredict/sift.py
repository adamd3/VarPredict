# Set up data and run SIFT

import os
import sys
import re
import numpy as np
import pandas as pd
import subprocess
import tempfile
import collections
import math
import ftplib
from itertools import chain
from pyfaidx import Fasta


## TODO:
##  - rename column containing the reference aa as 'ref'


def get_db(
    db_name = 'uniref90', # name of the db (uniref90 / uniref100)
    out_dir = './'): # output directory for storing db
    ftp_host='ftp.uniprot.org'
    if db_name=='uniref90':
        ftp_path='/pub/databases/uniprot/uniref/uniref90'
        filename='uniref90.fasta.gz'
    elif db_name=='uniref100':
        ftp_path='/pub/databases/uniprot/uniref/uniref100'
        filename='uniref100.fasta.gz'
    else:
        raise ValueError("`db_name` must be one of `uniref90` or `uniref100`")
    outf=os.path.join(out_dir, 'uniref_db.fasta.gz')
    ftp=ftplib.FTP(ftp_host)
    ftp.login()
    ftp.cwd(ftp_path)
    # ftp.retrlines("LIST")
    ftp.retrbinary(f"RETR {filename}", open(outf, 'wb').write)

def get_codon_pos(
    vars_data, # table of variants
    cds_annot, # annotations of protein-coding genes in reference strain
    gene_seq): # amino acid sequences of genes

    codon_pos_list = []
    for index, row in vars_data.iterrows():
        pos = int(row['pos'])
        gene = row['gene']
        ref_aa = row['ref']
        if gene in cds_annot['locus_tag'].tolist():
            annot = (cds_annot[cds_annot['locus_tag']==gene])
            gene_start = int((cds_annot[cds_annot['locus_tag']==gene])['start'])
            gene_stop = int((cds_annot[cds_annot['locus_tag']==gene])['end'])
            gene_strand = str(cds_annot.loc[cds_annot['locus_tag'] == gene, 'strand'].item())
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

    return codon_pos_list


def make_sift_files(
    fasta_file, # amino acid multi-fasta file
    vars_file, # file containing table of variants in protein-coding regions
    gff_file, # genome annotation file
    # counts_file, # normalised read counts file
    out_dir): # output directory

    # read input
    vars_data = pd.read_csv(vars_file, sep = "\t", low_memory=False)
    annotation = pd.read_csv(gff_file, sep = "\t")
    gene_seq = Fasta(fasta_file)

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
    cds_annot = gff_to_df(gff_file)

    vars_sub = vars_data.copy(deep=True)
    vars_sub['codon_pos'] = get_codon_pos(vars_sub, cds_annot, gene_seq)

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
    os.mkdir(os.path.join(out_dir, 'subst_files'))
    for gene_name in strain_var_dict.keys():
        all_vars = list(chain.from_iterable(strain_var_dict[gene_name].values()))
        all_vars = list(set([x for xs in all_vars for x in xs]))
        outfile = os.path.join(out_dir, 'subst_files', gene_name+'.subst')
        with open(outfile, 'a+') as f:
            for var in all_vars:
                f.write("%s\n" % var)

    return



def run_sift(fasta_file, vars_file, gff_file, out_dir):
    # download UNIPROT db
    get_db(out_dir = out_dir)

    # make SIFT input files
    make_sift_files(
        fasta_file = fasta_file,
        vars_file = vars_file,
        gff_file = gff_file,
        # counts_file = counts_file,
        out_dir = out_dir
    )

    cmd = "sift4g"
    cmd += " --query " + fasta_file
    cmd += " --subst " + os.path.join(out_dir, 'subst_files')
    cmd += " --database " + os.path.join(out_dir, 'uniref_db.fasta.gz')
    cmd += " --out " + os.path.join(out_dir, 'sift_out')
    # default params (not currently modifiable):
    cmd += " --kmer-length " + "5"
    cmd += " --max-candidates " + "5000"
    cmd += " --threads " + "1" # crashes when >1 core used
    cmd += " --evalue " + "0.0001"

    if not quiet:
        print("executing command: " + cmd)
    else:
        cmd += " > /dev/null"

    subprocess.run(cmd, shell=True, check=True)

    return