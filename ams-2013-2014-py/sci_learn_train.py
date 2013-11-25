import warnings
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
import sklearn.metrics as skMetrics
import argparse
import pandas as pd
import timeit
import numpy as np
import ml_metrics
import ast
import pickle


def train_and_predict(train_data_file, test_data_file, target_col, test_pred_file,
                      model_type, model_file, fit_args, test_metric, na_fill_value,
                      silent):

    start = timeit.default_timer()

    train_x = pd.read_csv(train_data_file)

    mappings = dict()
    for col in train_x.columns:
        if train_x[col].dtype == np.dtype('object'):
            s = np.unique(train_x[col].fillna(na_fill_value).values)
            mappings[col] = pd.Series([x[0] for x in enumerate(s)], index=s)
            train_x[col] = train_x[col].map(mappings[col]).fillna(na_fill_value)
        else:
            train_x[col] = train_x[col].fillna(na_fill_value)
    train_y = train_x[target_col]
    del train_x[target_col]

    x_cols = train_x.columns
    feat_importance_fun = lambda (fitted_model): fitted_model.feature_importances_
    staged_predict = lambda (fitted_model, pred_x): [fitted_model.predict(pred_x)]
    predict = lambda (fitted_model, pred_x): fitted_model.predict(pred_x)

    model = None
    if model_type == "RandomForestRegressor":
        model = RandomForestRegressor(**fit_args)
        model.fit(X=train_x, y=train_y)

    elif model_type == "GradientBoostingRegressor":
        model = GradientBoostingRegressor(**fit_args)
        model.fit(X=train_x, y=train_y)
        staged_predict = lambda (fitted_model, pred_x): fitted_model.staged_predict(pred_x)

    save_model(model=model, model_file=model_file)

    del train_x, train_y

    test_x = pd.read_csv(test_data_file)
    for col in test_x.columns:
        if col in mappings:
            test_x[col] = test_x[col].map(mappings[col]).fillna(na_fill_value)
        else:
            test_x[col] = test_x[col].fillna(na_fill_value)

    test_y = None
    if target_col in test_x.columns:
        test_y = test_x[target_col][test_x[target_col] != na_fill_value]
        if len(test_y) != len(test_x):
            test_y = None

    test_x = test_x[x_cols]

    test_pred = pd.DataFrame({'pred': predict((model, test_x))})
    if not silent and test_y is not None:
        print_stages(test_y=test_y, stage_predictions=staged_predict((model, test_x)), test_metric=test_metric)

    if not silent:
        feat_importance = feat_importance_fun(model)
        if feat_importance is not None:
            feat_importance = pd.DataFrame({'Features': x_cols,
                                            'Importance': feat_importance})
            pd.set_printoptions(max_columns=len(test_x.columns), max_rows=len(test_x.columns))
            print("Feature importances:")
            feat_importance.sort(columns='Importance', ascending=False, inplace=True)
            feat_importance.index = range(1, len(feat_importance) + 1)
            print(feat_importance)

    test_pred.to_csv(test_pred_file, index=False)
    stop = timeit.default_timer()
    if not silent:
        print "Time: %d s" % (stop - start)


def print_stages(test_y, stage_predictions, test_metric):
    if hasattr(ml_metrics, test_metric):
        eval_metric = getattr(ml_metrics, test_metric)
    else:
        eval_metric = getattr(skMetrics, test_metric)
    count = 0
    iters = []
    loss = []
    for prediction in stage_predictions:
        count += 1
        if count in [1, 5, 10, 30] or count % 50 == 0:
            iters.append(count)
            loss.append(eval_metric(test_y, prediction))
    loss_df = pd.DataFrame({'Iteration': iters, 'Loss': loss})
    loss_df.rename(columns={'Loss': test_metric}, inplace=True)
    pd.set_printoptions(max_columns=len(loss_df.columns), max_rows=len(loss_df))
    print("Loss:")
    print(loss_df)


def load_model(model_file):
    model_file = open(model_file, 'rb')
    return pickle.load(model_file)


def save_model(model, model_file):
    model_file = open(model_file, 'wb')
    pickle.dump(model, model_file, pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Train and predict data using some sklearn algorithms.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-train_data_file',
                        # default='../data/output-py/rf2_gefs_0_k_1_tr.csv',
                        type=str,
                        required=True,
                        help='CSV with training data.')

    parser.add_argument('-test_data_file',
                        # default='../data/output-py/rf2_gefs_0_k_1_test.csv',
                        type=str,
                        required=True,
                        help='CSV with testing data.')

    parser.add_argument('-target_col',
                        # default='power',
                        type=str,
                        required=True,
                        help='Name of target variable.')

    parser.add_argument('-test_pred_file',
                        # default='../data/output-py/rf2_gefs_0_k_1_test_pred.csv',
                        type=str,
                        required=True,
                        help='Path to output testing predictions.')

    parser.add_argument('-test_metric',
                        # default='mae',
                        type=str,
                        required=True,
                        help='Metric to compute on test set. Any metric on ml_metrics or sklearn.metrics')

    parser.add_argument('-model_type',
                        # default='RandomForestRegressor',
                        type=str,
                        help='Type of model to fit.',
                        required=True,
                        choices=["RandomForestRegressor",
                                 "GradientBoostingRegressor"])

    parser.add_argument('-model_file',
                        #default=None,
                        type=str,
                        required=True,
                        help='File to save the model to.')

    parser.add_argument('-na_fill_value',
                        default=-1,
                        type=int,
                        help='Value to fill in NAs.')

    parser.add_argument('-fit_args',
                        # default='{"criterion": "mse", "n_estimators": 1, "oob_score": True, "max_features": 10,' +
                        #         '"random_state": 67899, "compute_importances": True, "verbose": 2}',
                        type=str,
                        required=True,
                        help='String in dictionary form of fit params.')

    parser.add_argument('-silent',
                        default=False,
                        action='store_true',
                        help="Don't print execution information.")

    args = vars(parser.parse_args())

    args['fit_args'] = ast.literal_eval(args['fit_args'])
    for key in args['fit_args']:
        if args['fit_args'][key] in args:
            args['fit_args'][key] = args[args['fit_args'][key]]

    if not args['silent']:
        print(args)

    # train_data_file = args['train_data_file']
    # test_data_file = args['test_data_file']
    # target_col = args['target_col']
    # test_pred_file = args['test_pred_file']
    # model_type = args['model_type']
    # fit_args = args['fit_args']
    # test_metric = args['test_metric']
    # silent = args['silent']

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=Warning)
        train_and_predict(**args)
