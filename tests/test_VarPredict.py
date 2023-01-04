# Test pipeline
import sys, os
import tempfile
from VarPredict.rf_model import main

def test_rf(data_dir):

    with tempfile.TemporaryDirectory() as tmpdirname:

        geno_f = os.path.join(data_dir[0], 'gene_burden_scores.sub.txt')
        counts_f = os.path.join(data_dir[0], 'counts_matrix.sub.txt')
        meta_f = os.path.join(data_dir[0], 'metadata_merged.txt')

        ## run random-forest modelling
        sys.argv = [
            "random-forest", 
            "-o", tmpdirname, 
            "-g", geno_f, 
            "-c", counts_f, 
            "-m", meta_f
        ]
        main()


        # TODO: check that output matches expectation

    return
