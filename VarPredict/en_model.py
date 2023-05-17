# Elastic Net-penalised Logistic Regression model

import argparse
import numpy as np
import os
import pandas as pd
import pickle
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import (
    cross_val_score,
    GridSearchCV,
    StratifiedKFold
)
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import balanced_accuracy_score
from .is_valid import *
from .__init__ import __version__


def en_parser(parser):

    parser.description = 'Elastic Net-penalised Logistic Regression model'

    # input/output options
    io_opts = parser.add_argument_group('Input/output options')
    io_opts.add_argument(
        '-o',
        '--out_dir',
        dest = 'output_dir',
        required = True,
        help = 'Location of output directory, which should already exist',
        type = lambda x: is_valid_dir(parser, x)
    )
    io_opts.add_argument(
        '-g',
        '--genotypes',
        dest = 'geno_f',
        required = True,
        help = 'TSV file containing genotypes across strains (can be ' +
        'integers or continuous values)',
        type = argparse.FileType('r'),
        default = None
    )
    io_opts.add_argument(
        '-c',
        '--counts',
        dest = 'counts_f',
        required = True,
        help = 'TSV file containing normalised expression values per feature',
        type = argparse.FileType('r'),
        default = None
    )
    io_opts.add_argument(
        '-m',
        '--metadata',
        dest = 'meta_f',
        required = False,
        help = 'TSV file containing metadata for samples. Columns contain ' + 
        'model covariates.',
        type = argparse.FileType('r'),
        default = None
    )
    
    # hyperparameter options
    hyper_opts = parser.add_argument_group('Hyperparameter options')
    hyper_opts.add_argument(
        '-v',
        '--c_vals',
        dest = 'c_vals',
        nargs = '+',
        default = list(np.power(10.0, np.arange(-2, 2))),
        help = 'C values to try for cross-validation. Note: smaller ' +
            'C = stronger penalisation.'
    )
    hyper_opts.add_argument(
        '-l',
        '--l1_ratios',
        dest = 'l1_ratios',
        nargs = '+',
        default = np.arange(0, 1.10, 0.1),
        help = 'l1 ratios to try for cross-validation. Note: l1 penalty ' +
            '= lasso; l2 penalty = ridge'
    )

    parser.set_defaults(func = en_model)

    return parser

def en_model(args):
    print('Running Elastic Net-penalised Logistic Regression modelling')

    genotypes_data = pd.read_csv(args.geno_f, sep = '\t')
    counts_data = pd.read_csv(args.counts_f, sep = '\t')
    meta_data = pd.read_csv(args.meta_f, sep = '\t')

    args.output_dir = os.path.join(args.output_dir, '')

    genotypes_data = genotypes_data[counts_data.columns.tolist()]

    strains_present = counts_data.columns.values.tolist()[1:len(
        counts_data.columns)]
    meta_data = meta_data[meta_data['sample_name'].isin(strains_present)]
    strain_list = meta_data['sample_name'].tolist()
    counts_data = counts_data[['feature_id'] + strain_list]
    st_encod = pd.get_dummies(meta_data['majority_ST'], prefix='ST')
    st_encod.index = strain_list

    genotypes_t = genotypes_data.iloc[:,1:len(
        genotypes_data.columns)].transpose()
    genotypes_t.columns = genotypes_data['feature_id'].tolist()
    counts_t = counts_data.iloc[:,1:len(counts_data.columns)].transpose()
    feature_list = counts_data['feature_id'].tolist()
    counts_t.columns = feature_list

    # drop genes with low variance
    counts_var = counts_t.var()
    cols_to_drop = counts_var[counts_var < 1e-5].index
    counts_t = counts_t.drop(cols_to_drop, axis=1)
    feature_list_sub = list(counts_t.columns)

    vars_st = genotypes_t.merge(st_encod, left_index=True, right_index=True)

    acc_df = en_nested_cv(feature_list_sub, counts_t, vars_st, args)
    outf1 = os.path.join(args.output_dir, 'acc_df_en.txt')
    acc_df.to_csv(outf1, index=False, sep='\t')

    coef_df_combined = en_coefs(feature_list_sub, counts_t, vars_st, args)
    outf2 = os.path.join(args.output_dir, 'coef_df_en.txt')
    coef_df_combined.to_csv(outf2, index=False, sep='\t')

    return

def en_nested_cv(feature_list, counts_t, vars_st, args):

    X = vars_st.values
    scaler = StandardScaler()
    X = scaler.fit_transform(X)  

    cv_outer = StratifiedKFold(n_splits=5, shuffle=True)
    cv_inner = StratifiedKFold(n_splits=3, shuffle=True)
    model = LogisticRegression(penalty = 'elasticnet', solver = 'saga')

    grid = {
        'C': args.c_vals, 
        'l1_ratio': args.l1_ratios
    }

    mean_outer_accs={}
    for feat in feature_list:
        try:
            y = pd.qcut(counts_t[feat].values, 2, labels = [0,1])
            search = GridSearchCV(
                model, grid, scoring='balanced_accuracy', cv=cv_inner, 
                refit=True, n_jobs=-1
            )
            scores = cross_val_score(
                search, X, y, cv=cv_outer, 
                scoring = 'balanced_accuracy', n_jobs=-1
            )
            mean_acc = np.mean(scores)
            mean_std = np.std(scores)
            mean_outer_accs[feat] = mean_acc
        except ValueError:
            mean_outer_accs[feat] = np.nan

    acc_df = pd.DataFrame.from_dict(mean_outer_accs, orient='index')
    acc_df = acc_df.rename_axis('feature').reset_index()
    acc_df = acc_df.rename(columns={0:'mean_acc'})
    acc_df = acc_df.dropna()
    acc_df = acc_df.sort_values('mean_acc', ascending=False)
    acc_df = acc_df.reset_index(drop=True)

    return acc_df

def en_coefs(feature_list, counts_t, vars_st, args):

    X = vars_st.values
    scaler = StandardScaler()
    X = scaler.fit_transform(X)  

    cv_split = StratifiedKFold(n_splits=5, shuffle=True)
    model = LogisticRegression(penalty = 'elasticnet', solver = 'saga')

    grid = {
        'C': args.c_vals, 
        'l1_ratio': args.l1_ratios
    }

    coefs_dict = {}
    best_models = {}
    for feat in feature_list:
        y = pd.qcut(counts_t[feat].values, 2, labels = [0,1])
        try:
            search = GridSearchCV(
                model, grid, scoring='balanced_accuracy', 
                cv=cv_split, refit=True
            )
            result = search.fit(X, y)
            best_model = result.best_estimator_
            best_models[feat] = best_model
            coefs_df = pd.DataFrame((best_model.coef_).transpose())
            coefs_df['var'] = vars_st.columns
            coefs_df = coefs_df.rename(columns={0: 'coef'})
            coefs_df = coefs_df.sort_values('coef', ascending=False)
            coefs_df = coefs_df.reset_index(drop=True)
            coefs_df['feat'] = feat 
            coefs_df = coefs_df[coefs_df['coef']!=0]
            coefs_dict[feat] = coefs_df
        except ValueError:
            coefs_dict[feat] = np.nan

    outf3 = os.path.join(args.output_dir, 'best_models_en.pkl')
    with open(outf3, 'wb') as file:
        pickle.dump(best_models, file)

    df_names = []
    for i in coefs_dict.keys():
        vali = coefs_dict[i]
        if isinstance(vali, pd.DataFrame):
            df_names.append(coefs_dict[i])
    coefs_df_combined = pd.concat(df_names)

    return coefs_df_combined

def main():
    parser = argparse.ArgumentParser()
    parser = en_parser(parser)
    args = parser.parse_args()

    args.func(args)  

    return

if __name__ == '__main__':
    main()