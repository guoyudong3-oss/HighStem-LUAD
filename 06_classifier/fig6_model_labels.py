"""Mapping from internal model keys to the scikit-learn classes defined in model_zoo(seed) of make_fig6_reference_like.py; used for figure and table labels."""
MODEL_LABELS = {
    "classif.log_reg":     {"display": "LogisticRegression",             "sklearn_class": "sklearn.pipeline.Pipeline[StandardScaler -> LogisticRegression(max_iter=5000)]"},
    "classif.xgboost":     {"display": "HistGradientBoosting",           "sklearn_class": "sklearn.ensemble.HistGradientBoostingClassifier(max_depth=4, learning_rate=0.05, max_iter=300)"},
    "classif.lda":         {"display": "LinearDiscriminantAnalysis",     "sklearn_class": "sklearn.pipeline.Pipeline[StandardScaler -> LinearDiscriminantAnalysis]"},
    "classif.svm":         {"display": "SVC",                            "sklearn_class": "sklearn.pipeline.Pipeline[StandardScaler -> SVC(probability=True)]"},
    "classif.naive_bayes": {"display": "GaussianNB",                     "sklearn_class": "sklearn.naive_bayes.GaussianNB"},
    "classif.kknn":        {"display": "KNeighbors",                     "sklearn_class": "sklearn.pipeline.Pipeline[StandardScaler -> KNeighborsClassifier(n_neighbors=31)]"},
    "classif.rpart":       {"display": "DecisionTree",                   "sklearn_class": "sklearn.tree.DecisionTreeClassifier(max_depth=4)"},
    "classif.ranger":      {"display": "RandomForest",                   "sklearn_class": "sklearn.ensemble.RandomForestClassifier(n_estimators=400)"},
}
def relabel(df, col="model"):
    out = df.copy()
    out.insert(out.columns.get_loc(col) + 1, "model_display", out[col].map(lambda k: MODEL_LABELS[k]["display"]))
    out.insert(out.columns.get_loc(col) + 2, "sklearn_class", out[col].map(lambda k: MODEL_LABELS[k]["sklearn_class"]))
    return out
