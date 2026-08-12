## Using EAGLES scores for phenotype predictions
This vignette demonstrates how to use EAGLES scores for downstream phenotype prediction. As an example, we show how we used EAGLES pathways scores to predict immunotherapy response across multiple melanoma cohorts. This guide asses you have already generated EAGLES scores (see README.md for installation and general instructions). Here we focus on loading EAGLES scores, running classifiers for predictions, and visualizing results.

## Loading scores
EAGLES scores are outputted as matrices with rowd corresponding sampled and columnd to a gene-level features.
Below is a minimal example of loading in EAGLES scores created from finemapped eQTLs.
```
def load_scores(scores_path):
    path = scores_path
    raw = pd.read_csv(path, sep="\t", header=None, low_memory=False)
    
    # gene IDs
    gene_ids = raw.iloc[1].tolist()
    
    # sample info
    df_scores = raw.iloc[3:, 1:].copy()
    
    # sample IDs
    samples = raw.iloc[3:, 0].tolist()
    
    df_scores.index = samples
    df_scores.columns = gene_ids[1:] 
    df_scores.index.name = "sample"
    
    return df_scores
```
Then to get scores, use `df_scores_finemap = load_scores("projects/eagles/icb_melanoma_scores/finemap/scores.tsv")` for example.
For our purposes, we split up the scores into 4 pathways (interferon alpha response pathway, interferon gamma response pathway, reactive oxygen species pathway, and DNA repair pathway) and used classifiers to test each of the pathway gene sets. We used Hallmark gene sets to define pathway membership.
```
# hallmark membership
df_hallmark = pd.read_csv("projects/eagles/icb_scores_blood/hallmark_members.tsv", sep="\t")
df_hallmark = df_hallmark.set_index("Unnamed: 0")

# extract genes from pathway
def genes_for(pathway):
    return df_hallmark.index[df_hallmark[pathway] == 1].tolist()

# gene sets
genes_ifna = genes_for("HALLMARK_INTERFERON_ALPHA_RESPONSE")
genes_ifng = genes_for("HALLMARK_INTERFERON_GAMMA_RESPONSE")
genes_ros = genes_for("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")
genes_dnar = genes_for("HALLMARK_DNA_REPAIR")

# subset genes
df_ifna = subset_genes(df_scores_finemap, genes_ifna)
df_ifng = subset_genes(df_scores_finemap, genes_ifng)
df_ros  = subset_genes(df_scores_finemap, genes_ros)
df_dnar = subset_genes(df_scores_finemap, genes_dnar)
```

## Running phenotype prediction
Below we have skeleton code for XGBoost classifier we used for immunotherapy response prediction. The outputs from this function allow us to plot ROC curves, compares AUC values, and look at SHAP summary plots.
```
# XGBoost model
def run_xgb_simple(df, model_label, clinical):
    df_norm = df.apply(pd.to_numeric, errors="coerce")
    
    # cleaning up sample IDs
    df_norm.index = (
        df_norm.index
        .str.replace(r"^0_", "", regex=True)
        .str.replace(r"_normal_sample$", "", regex=True)
    )

    # match with clinical data
    common = df_norm.index.intersection(clinical.index)
    responses = clinical.loc[common, "Response_standardized"]
    valid = responses.dropna().index

    X_all = df_norm.loc[valid]
    y_all = responses.loc[valid].astype(int)
    cohorts_all = clinical.loc[valid, "study"]

    # train/test split by cohort
    train_cohorts = []
    test_cohorts = []

    train_idx = cohorts_all[cohorts_all.isin(train_cohorts)].index
    test_idx  = cohorts_all[cohorts_all.isin(test_cohorts)].index

    X_train_full = X_all.loc[train_idx]
    y_train_full = y_all.loc[train_idx]

    X_test = X_all.loc[test_idx]
    y_test = y_all.loc[test_idx]
    test_cohorts_series = cohorts_all.loc[test_idx]

    # only evaluate cohorts with >=10 samples
    cohort_counts = test_cohorts_series.value_counts()
    valid_cohorts = cohort_counts[cohort_counts >= 10].index

    # feature alignment
    feature_order = X_train_full.columns.tolist()
    X_train_full = X_train_full[feature_order]
    X_test = X_test[feature_order]

    # parameters can be changed / tuned
    xgb_params = dict(
        n_estimators=250,
        max_depth=5,
        learning_rate=0.05,
        subsample=0.6,
        colsample_bytree=0.6,
        min_child_weight=4,
        gamma=2.0,
        reg_alpha=1.5,
        reg_lambda=5.0,
        objective="binary:logistic",
        eval_metric="logloss"
    ) 

    model = XGBClassifier(
        **xgb_params,
        random_state=0,
        n_jobs=-1
    )

    # fit on full training cohort
    model.fit(X_train_full, y_train_full)

    # training predictions
    y_train_pred = model.predict_proba(X_train_full)[:, 1]
    auc_train = roc_auc_score(y_train_full, y_train_pred)

    # compute training ROC curve
    fpr_train, tpr_train, _ = roc_curve(y_train_full, y_train_pred)

    # SHAP values
    explainer = shap.TreeExplainer(model)
    shap_train = explainer.shap_values(X_train_full)

    # test predictions
    y_pred_proba = model.predict_proba(X_test)[:, 1]
    y_pred_series = pd.Series(y_pred_proba, index=X_test.index)

    # ROC per cohort
    auc_summary = {}
    for cohort in valid_cohorts:
        idx = test_cohorts_series[test_cohorts_series == cohort].index
        y_true = y_test.loc[idx]
        y_scores = y_pred_series.loc[idx]

        if len(np.unique(y_true)) < 2:
            continue

        fpr, tpr, _ = roc_curve(y_true, y_scores)
        auc = roc_auc_score(y_true, y_scores)

        auc_summary[cohort] = {"fpr": fpr, "tpr": tpr, "auc": auc}

    return (
        auc_summary,
        model_label,
        model,
        X_test,
        y_test,
        X_all,
        X_train_full,
        y_train_full,
        y_train_pred,
        auc_train,
        shap_train,
        fpr_train,
        tpr_train
    )
```

We can then run these classifiers on EAGLES scores. Each call returns the trained models, train/test predictions, AUC values per cohort, and SHAP values.
```
# running xgb on gene sets
auc_ifng = run_xgb_simple(df_ifng, "IFNG", clinical)
auc_ifna = run_xgb_simple(df_ifna, "IFNA", clinical)
auc_ros  = run_xgb_simple(df_ros,  "ROS",  clinical)
auc_dnar = run_xgb_simple(df_dnar, "DNAR", clinical)
```

## Visualization
We can visualize ROC curves for each pathway using the returned summaries.
```
# plot 
results = {
    "IFN$\\alpha$": auc_ifna[0],
    "IFN$\\gamma$": auc_ifng[0],
    "ROS":  auc_ros[0],
    "DNAR": auc_dnar[0]
}

colors = {
    "IFN$\\alpha$": "tab:blue",
    "IFN$\\gamma$": "tab:red",
    "ROS":  "tab:green",
    "DNAR": "tab:purple",
}

plt.figure(figsize=(8, 8))

for label, auc_summary in results.items():
    d = auc_summary["pooled_test"]
    fpr = d["fpr"]
    tpr = d["tpr"]
    auc = d["auc"]

    plt.plot(
        fpr,
        tpr,
        color=colors[label],
        linewidth=2.5,
        label=f"{label} (AUC={auc:.3f})"
    )

plt.plot([0, 1], [0, 1], "k--", alpha=0.4)

plt.xlabel("False Positive Rate", fontsize=18)
plt.ylabel("True Positive Rate", fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14)
plt.tight_layout()
plt.show()
```

Finally, we can look at SHAP summary plots to help identify which genes contribute most to immunotherapy response prediction.
```
# pull model results
(
    auc_ifng_xgb_summary,
    _,
    model_ifng_xgb,
    X_test_ifng_xgb,
    y_test_ifng_xgb,
    X_all_ifng_xgb,
    X_train_ifng_xgb,
    y_train_ifng_xgb,
    y_train_pred_ifng_xgb,
    auc_train_ifng_xgb,
    shap_train_ifng_xgb,
    _,
    _
) = run_xgb_simple(df_ifng, "IFNG", clinical)

# compute SHAP
expl_ifng = shap.TreeExplainer(model_ifng_xgb)
shap_ifng_xgb = expl_ifng.shap_values(X_test_ifng_xgb)

plt.figure(figsize=(8, 6))

shap.summary_plot(
    shap_ifng_xgb,
    X_ifng_named_xgb,
    plot_type="dot",
    max_display=10,
    show=False
)

plt.title("IFNG SHAP", fontsize=18)
plt.tight_layout()
plt.show()
```
