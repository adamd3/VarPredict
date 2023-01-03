# Random Forest model

import os
import argparse

from .__init__ import __version__


def rf_parser(parser):

    parser.description = "Random Forest model"

    # hyperparameter options
    hyper_opts = parser.add_argument_group("Hyperparameter options")
    hyper_opts.add(
            "-d",
            "--max_depth",
            dest = "max_depth",
            nargs = '+',
            default = [5,10, None],
            help = "The number of splits that each decision tree is allowed to make"
            )
    hyper_opts.add(
            "-f",
            "--max_features",
            dest = "max_features",
            nargs = '+',
            default = ['sqrt', 'log2', None],
            help = "The number of features to consider when looking for the best split"
            )
    hyper_opts.add(
            "-s",
            "--min_samples_split",
            dest = "min_samples_split",
            nargs = '+',
            default = [2, 4, 8],
            help = "The minimum number of samples required to split an internal node"
            )
    hyper_opts.add(
            "-e",
            "--n_estimators",
            dest = "n_estimators",
            nargs = '+',
            default = [300],
            help = "The number of trees in the forest"
            )

    parser.set_defaults(func = rf_model)

    return parser

def rf_model(args):
    print("Running Random Forest modelling")

    genotypes_data = pd.read_csv(args.geno_f, sep = "\t")
    counts_data = pd.read_csv(args.counts_f, sep = "\t")
    meta_data = pd.read_csv(args.meta_f, sep = "\t")
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
    # drop genes with a single value (no variation)
    counts_var = counts_t.var()
    cols_to_drop = counts_var[counts_var < 1e-5].index
    counts_t = counts_t.drop(cols_to_drop, axis=1)
    feature_list_sub = list(counts_t.columns)
    # merge variants + STs
    vars_st = genotypes_t.merge(st_encod, left_index=True, right_index=True)
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

    ## Nested CV using `cross_val_score` convenience function:
    ## Inner CV to choose hyperparams + best model
    ## when refit=True, the entire training set will be used to refit the model
    mean_outer_accs={}
    for feat in feature_list_sub:
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
    ## Make table of results
    acc_df = pd.DataFrame.from_dict(mean_outer_accs, orient='index')
    acc_df = acc_df.rename_axis("feature").reset_index()
    acc_df = acc_df.rename(columns={0:'mean_acc'})
    acc_df = acc_df.dropna()
    acc_df = acc_df.sort_values('mean_acc', ascending=False)
    acc_df = acc_df.reset_index(drop=True)
    outf = os.path.join(args.output_dir, 'acc_df_rf.txt')
    acc_df.to_csv(outf, index=False, sep='\t')

    return

def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser()
    parser = rf_parser(parser)
    args = parser.parse_args()

    # run modelling
    args.func(args)

    return

if __name__ == "__main__":
    main()