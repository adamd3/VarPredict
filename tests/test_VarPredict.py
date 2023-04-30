# Test pipeline

import os
import sys
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

        outf1 = os.path.join(tmpdirname, 'acc_df_rf.txt')
        outf2 = os.path.join(tmpdirname, 'imp_df_rf.txt')

        assert os.path.isfile(outf1)
        assert os.path.isfile(outf2)

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

        outf1 = os.path.join(tmpdirname, 'acc_df_en.txt')
        outf2 = os.path.join(tmpdirname, 'coef_df_en.txt')

        assert os.path.isfile(outf1)
        assert os.path.isfile(outf2)

    return

