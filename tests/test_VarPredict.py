# Test pipeline
import sys, os
import tempfile
from BactVarMiner.BactVarMiner import main

def test_pipeline(data_dir):

    with tempfile.TemporaryDirectory() as tmpdirname:

        vars_f = os.path.join(data_dir[0], 'sample2variant_coding.head.txt')
        gff_f = os.path.join(data_dir[0], 'Pseudomonas_aeruginosa_PAO1.gff')
        fasta_f = os.path.join(data_dir[0], 'Pseudomonas_aeruginosa_PAO1_107.faa')

        ## run pipeline
        sys.argv = [
            "run", 
            "-o", tmpdirname, 
            "-v", vars_f, 
            "-g", gff_f, 
            "-f", fasta_f
        ]
        main()

        # check directory is present
        outdir = os.path.join(tmpdirname, 'subst_files')


    return
