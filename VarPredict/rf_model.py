# Random Forest model

import os
import argparse
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score, KFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import balanced_accuracy_score


from .__init__ import __version__


def rf_parser(parser):

    parser.description = "Random Forest model"

    # input/output options
    io_opts = parser.add_argument_group("Input/output options")
    io_opts.add_argument(
        "-o",
        "--out_dir",
        dest = "output_dir",
        required = True,
        help = "Location of output directory, which should already exist",
        type = lambda x: is_valid_dir(parser, x)
    )
    io_opts.add_argument(
        "-g",
        "--genotypes",
        dest = "geno_f",
        required = True,
        help = "TSV file containing genotypes across strains (can be " +
        "integers or continuous values)",
        type = argparse.FileType("r"),
        default = None
    )
    io_opts.add_argument(
        "-c",
        "--counts",
        dest = "counts_f",
        required = True,
        help = "TSV file containing normalised expression values per feature",
        type = argparse.FileType("r"),
        default = None
    )
    io_opts.add_argument(
        "-m",
        "--metadata",
        dest = "meta_f",
        required = False,
        help = "TSV file containing metadata for samples. Columns contain " + 
        "model covariates.",
        type = argparse.FileType("r"),
        default = None
    )
    
    # hyperparameter options
    hyper_opts = parser.add_argument_group("Hyperparameter options")
    hyper_opts.add_argument(
            "-d",
            "--max_depth",
            dest = "max_depth",
            nargs = '+',
            default = [5,10, None],
            help = "The number of splits that each decision tree is allowed to make"
            )
    hyper_opts.add_argument(
            "-f",
            "--max_features",
            dest = "max_features",
            nargs = '+',
            default = ['sqrt', 'log2', None],
            help = "The number of features to consider when looking for the best split"
            )
    hyper_opts.add_argument(
            "-s",
            "--min_samples_split",
            dest = "min_samples_split",
            nargs = '+',
            default = [2, 4, 8],
            help = "The minimum number of samples required to split an internal node"
            )
    hyper_opts.add_argument(
            "-e",
            "--n_estimators",
            dest = "n_estimators",
            nargs = '+',
            default = [300],
            help = "The number of trees in the forest"
            )

    # add functions to be run (these will be called in `main()`)
    parser.set_defaults(func = rf_model)

    return parser

def rf_model(args):
    print("Running Random Forest modelling")

    # read input files
    genotypes_data = pd.read_csv(args.geno_f, sep = "\t")
    counts_data = pd.read_csv(args.counts_f, sep = "\t")
    meta_data = pd.read_csv(args.meta_f, sep = "\t")

    # ensure valid file path from output_dir
    args.output_dir = os.path.join(args.output_dir, "")

    # subset variants to strains present in expression data
    genotypes_data = genotypes_data[counts_data.columns.tolist()]

    ## Subset metadata to strains in counts
    strains_present = counts_data.columns.values.tolist()[1:len(counts_data.columns)]
    meta_data = meta_data[meta_data['sample_name'].isin(strains_present)]
    strain_list = meta_data['sample_name'].tolist()
    counts_data = counts_data[['feature_id'] + strain_list]
    st_encod = pd.get_dummies(meta_data['majority_ST'], prefix='ST')
    st_encod.index = strain_list

    # transpose the tables
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

    # merge variants + STs
    vars_st = genotypes_t.merge(st_encod, left_index=True, right_index=True)

    # estimate model accuracy via nested CV
    acc_df = rf_nested_cv(feature_list_sub, counts_t, vars_st)
    outf1 = os.path.join(args.output_dir, 'acc_df_rf.txt')
    acc_df.to_csv(outf1, index=False, sep='\t')

    # get importance scores by model for predictors
    imp_df_combined = rf_imp_scores(feature_list_sub, counts_t, vars_st)
    outf2 = os.path.join(args.output_dir, 'imp_df_rf.txt')
    imp_df_combined.to_csv(outf2, index=False, sep='\t')

    return

def rf_nested_cv(feature_list, counts_t, vars_st):
    print("Estimating Random Forest model accuracy via nested CV")

    X = vars_st.values
    cv_outer = StratifiedKFold(n_splits=5, shuffle=True)
    cv_inner = StratifiedKFold(n_splits=3, shuffle=True)
    model = RandomForestClassifier()
    grid = {
        'max_depth': args.max_depth, 
        'max_features': args.max_features,
        'min_samples_split': args.min_samples_split,
        'n_estimators': args.n_estimators
    }

    # get mean model accuracies from nested CV
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

    ## export table of results
    acc_df = pd.DataFrame.from_dict(mean_outer_accs, orient='index')
    acc_df = acc_df.rename_axis("feature").reset_index()
    acc_df = acc_df.rename(columns={0:'mean_acc'})
    acc_df = acc_df.dropna()
    acc_df = acc_df.sort_values('mean_acc', ascending=False)
    acc_df = acc_df.reset_index(drop=True)

    return acc_df

def rf_imp_scores(feature_list, counts_t, vars_st):
    print("Generating Random Forest Gini Importance Scores for predictors")

    imp_dict = {}
    for feat in feature_list:
        X = vars_st.values
        y = pd.qcut(counts_t[feat].values, 2, labels = [0,1])
        model = RandomForestClassifier()
        try:
            search = GridSearchCV(
                model, grid, scoring='balanced_accuracy', 
                cv=cv_outer, refit=True
            )
            result = search.fit(X, y)
            best_model = result.best_estimator_
            imp_df = pd.DataFrame((best_model.feature_importances_).transpose())
            imp_df['var'] = vars_st.columns
            imp_df = imp_df.rename(columns={0: 'importance'})
            imp_df = imp_df.sort_values('importance', ascending=False)
            imp_df = imp_df.reset_index(drop=True)
            imp_df['feat'] = feat 
            imp_df = imp_df[imp_df['importance']>0]
            imp_dict[feat] = imp_df
        except ValueError:
            imp_dict[feat] = np.nan

    # remove nans + combine dfs
    df_names = []
    for i in imp_dict.keys():
        vali = imp_dict[i]
        if isinstance(vali, pd.DataFrame):
            df_names.append(imp_dict[i])
    imp_df_combined = pd.concat(df_names)

    return imp_df_combined

def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser()
    parser = rf_parser(parser)
    args = parser.parse_args()

    # AD: the value for `func` is defined above in parser
    args.func(args)  

    return

if __name__ == "__main__":
    main()