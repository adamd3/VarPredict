# Test pipeline

import sys, os
import tempfile
from VarPredict.rf_model import main as rf_main
from VarPredict.en_model import main as en_main

def test_rf(data_dir):

    with tempfile.TemporaryDirectory() as tmpdirname:

        geno_f = os.path.join(data_dir[0], 'gene_burden_scores.sub.txt')
        counts_f = os.path.join(data_dir[0], 'counts_matrix.sub.txt')
        meta_f = os.path.join(data_dir[0], 'metadata_merged.txt')

        sys.argv = [
            'random-forest', 
            '-o', tmpdirname, 
            '-g', geno_f, 
            '-c', counts_f, 
            '-m', meta_f
        ]
        rf_main()

        assert os.path.isfile(tmpdirname + 'acc_df_rf.txt')
        assert os.path.isfile(tmpdirname + 'imp_df_rf.txt')

    return


def test_en(data_dir):

    with tempfile.TemporaryDirectory() as tmpdirname:

        geno_f = os.path.join(data_dir[0], 'gene_burden_scores.sub.txt')
        counts_f = os.path.join(data_dir[0], 'counts_matrix.sub.txt')
        meta_f = os.path.join(data_dir[0], 'metadata_merged.txt')

        sys.argv = [
            'elastic-net', 
            '-o', tmpdirname, 
            '-g', geno_f, 
            '-c', counts_f, 
            '-m', meta_f
        ]
        en_main()

        assert os.path.isfile(tmpdirname + 'acc_df_en.txt')
        assert os.path.isfile(tmpdirname + 'coef_df_en.txt')

    return

