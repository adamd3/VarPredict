import gffpandas.gffpandas as gffpd

# Converts a GFF file to a df of gene annotations
def gff_to_df(gff_f):
    annot_dat = (gffpd.read_gff3(gff_f)).df
    cds_annot = annot_dat[annot_dat['type']=='CDS']
    annot_split = cds_annot['attributes'].str.split(";")
    locus_tags = (annot_split.str[1]).str.replace(r'Parent=gene-','')
    gene_names = ((annot_split.str[6]).str.split("=")).str[1]
    cds_annot['locus_tag'] = locus_tags
    cds_annot['gene_name'] = gene_names
    cds_annot = cds_annot[['locus_tag', 'gene_name', 'start', 'end', 'strand']]
    return(cds_annot)
