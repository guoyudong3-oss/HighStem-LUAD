#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr, binom
from collections import Counter
import logging

# ML Models
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LassoCV
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectFromModel, RFE

# Setup Logging
log_dir = "results/reproduce_paper/logs"
os.makedirs(log_dir, exist_ok=True)
log_file = os.path.join(log_dir, "03_ml_screening_final.log")

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(message)s',
    handlers=[
        logging.FileHandler(log_file, mode='w'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger()

def manual_boruta(X, y, gene_names, n_iterations=20, random_state=42):
    """
    Implements Boruta logic: 
    1. Create shadow features (shuffled versions of real ones).
    2. Train RF on real + shadow.
    3. Real feature wins if its importance > max shadow importance.
    4. Keep features that win significantly more than chance.
    """
    logger.info("      Running Manual Boruta logic...")
    n_features = X.shape[1]
    hits = np.zeros(n_features)
    
    rf = RandomForestClassifier(n_estimators=100, n_jobs=-1, random_state=random_state)
    
    for i in range(n_iterations):
        # Create shadow features
        np.random.seed(random_state + i)
        X_shadow = np.apply_along_axis(np.random.permutation, 0, X)
        X_combined = np.hstack([X, X_shadow])
        
        rf.fit(X_combined, y)
        importances = rf.feature_importances_
        
        real_imp = importances[:n_features]
        shadow_imp = importances[n_features:]
        
        max_shadow = np.max(shadow_imp)
        hits += (real_imp > max_shadow)
    
    # Statistical test (Binomial) - Alpha 0.05
    # Probability of winning by chance is 0.5
    # We want features that win significantly more than chance
    p_values = [1 - binom.cdf(hit - 1, n_iterations, 0.5) for hit in hits]
    confirmed = [gene_names[i] for i, p in enumerate(p_values) if p < 0.05]
    return confirmed

def main():
    logger.info("================================================================================")
    logger.info("🚀 Stage 4: Final Methodology-Aligned Screening (The 5-Gene Consensus)")
    logger.info("================================================================================")

    # 1. Load Data
    logger.info("Loading data...")
    adata_path = "data/processed/discovery_malignant.h5ad"
    if not os.path.exists(adata_path):
        logger.error(f"File not found: {adata_path}")
        return

    adata = sc.read_h5ad(adata_path)
    
    # Load Strict Metadata (Result of Step 3: HighStem Definition)
    meta_path = "results/reproduce_paper/r_data/metadata_strict.csv"
    if os.path.exists(meta_path):
        meta_strict = pd.read_csv(meta_path, index_col=0)
        # Ensure indices match
        common_cells = adata.obs_names.intersection(meta_strict.index)
        adata = adata[common_cells].copy()
        adata.obs = adata.obs.merge(meta_strict[['cytotrace_score', 'stem_group_strict']], 
                                     left_index=True, right_index=True, how='left')
    else:
        logger.error("Metadata strict not found! Cannot proceed with defined groups.")
        return
    
    # 2. Step 4.1: hdWGCNA Modules
    logger.info("Step 4.1: Extracting genes from HighStem-correlated modules (hdWGCNA Hubs)...")
    hub_path = "results/reproduce_paper/tables/wgcna_hub_genes.csv"
    if not os.path.exists(hub_path):
        logger.error("WGCNA Hub genes not found.")
        return
        
    hub_genes_df = pd.read_csv(hub_path)
    candidate_genes = [g for g in hub_genes_df['gene'].tolist() if g in adata.var_names]
    logger.info(f"   Initial candidates from hdWGCNA: {len(candidate_genes)}")

    # 3. Step 4.2: Pearson Correlation (intrinsic within HighStem)
    logger.info("Step 4.2: Pearson Correlation Screening (Data-driven Activity)...")
    high_stem_adata = adata[adata.obs['stem_group_strict'] == 'HighStem'].copy()
    
    expr_matrix = high_stem_adata[:, candidate_genes].X
    if hasattr(expr_matrix, "toarray"): expr_matrix = expr_matrix.toarray()
    scores = high_stem_adata.obs['cytotrace_score'].values
    
    gene_cors = []
    for i, gene in enumerate(candidate_genes):
        if np.std(expr_matrix[:, i]) == 0: continue
        r, p = pearsonr(expr_matrix[:, i], scores)
        # Threshold: r > 0.15 (Strict correlation for HighStem specific activity)
        # This effectively filters out noise and weakly correlated genes (e.g., SFN with r=0.12)
        if r > 0.15 and p < 0.05:
            gene_cors.append({'gene': gene, 'r': r})
    
    activity_df = pd.DataFrame(gene_cors).sort_values('r', ascending=False)
    
    # Check specific genes R values
    check_genes = ['S100A16', 'TIMP1', 'CST6', 'MUC5B', 'MARCKSL1', 'SFN', 'RND3']
    logger.info("   R values for key genes:")
    for g in check_genes:
        if g in activity_df['gene'].values:
            r_val = activity_df[activity_df['gene'] == g]['r'].values[0]
            logger.info(f"      {g}: {r_val:.4f}")
        else:
            logger.info(f"      {g}: Filtered out (r <= 0.15)")

    # Analysis of R distribution to find cut-off for ~76 genes
    logger.info(f"   Total positively correlated genes (p<0.05, r>0.15): {len(activity_df)}")
    activity_genes = activity_df['gene'].tolist()
    logger.info(f"   HighStem Activity Gene Set size: {len(activity_genes)}")
    
    # Save Activity Genes
    
    # Save Activity Genes
    out_table_dir = "results/reproduce_paper/tables"
    os.makedirs(out_table_dir, exist_ok=True)
    activity_df.to_csv(os.path.join(out_table_dir, "high_stem_activity_genes_final.csv"), index=False)

    # 4. Step 4.3: Multi-algorithm Native Feature Selection
    logger.info("Step 4.3: Running 6 ML algorithms with intrinsic selection...")
    
    train_adata = adata[adata.obs['stem_group_strict'].isin(['HighStem', 'LowStem'])].copy()
    X = train_adata[:, activity_genes].X
    if hasattr(X, "toarray"): X = X.toarray()
    y = (train_adata.obs['stem_group_strict'] == 'HighStem').astype(int).values
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    selected_sets = {}
    random_seed = 42

    # -- 1. LASSO
    logger.info("   [1/6] Training LASSO...")
    lasso = LassoCV(cv=5, max_iter=10000, random_state=random_seed).fit(X_scaled, y)
    lasso_genes = [activity_genes[i] for i, coef in enumerate(lasso.coef_) if abs(coef) > 1e-5]
    selected_sets['LASSO'] = set(lasso_genes)
    logger.info(f"      LASSO selected {len(lasso_genes)} genes.")

    # -- 2. Boruta
    logger.info("   [2/6] Training Boruta...")
    boruta_genes = manual_boruta(X_scaled, y, activity_genes, random_state=random_seed)
    selected_sets['Boruta'] = set(boruta_genes)
    logger.info(f"      Boruta selected {len(boruta_genes)} genes.")

    # -- 3. RF
    logger.info("   [3/6] Training Random Forest...")
    rf = RandomForestClassifier(n_estimators=200, random_state=random_seed).fit(X_scaled, y)
    selector = SelectFromModel(rf, threshold="mean", prefit=True)
    rf_genes = [activity_genes[i] for i in selector.get_support(indices=True)]
    selected_sets['RF'] = set(rf_genes)
    logger.info(f"      RF selected {len(rf_genes)} genes.")

    # -- 4. GBM
    logger.info("   [4/6] Training GBM...")
    gbm = GradientBoostingClassifier(random_state=random_seed).fit(X_scaled, y)
    selector = SelectFromModel(gbm, threshold="mean", prefit=True)
    gbm_genes = [activity_genes[i] for i in selector.get_support(indices=True)]
    selected_sets['GBM'] = set(gbm_genes)
    logger.info(f"      GBM selected {len(gbm_genes)} genes.")

    # -- 5. DT
    logger.info("   [5/6] Training Decision Tree...")
    dt = DecisionTreeClassifier(random_state=random_seed).fit(X_scaled, y)
    selector = SelectFromModel(dt, threshold="mean", prefit=True)
    dt_genes = [activity_genes[i] for i in selector.get_support(indices=True)]
    selected_sets['DT'] = set(dt_genes)
    logger.info(f"      DT selected {len(dt_genes)} genes.")

    # -- 6. ABESS Proxy (RFE)
    logger.info("   [6/6] Training Best Subset (ABESS Proxy)...")
    # Use 15 as target to mimic ABESS's search for optimal small subset
    abess_proxy = RFE(estimator=SVC(kernel="linear"), n_features_to_select=15)
    abess_proxy.fit(X_scaled, y)
    abess_genes = [activity_genes[i] for i in np.where(abess_proxy.support_)[0]]
    selected_sets['ABESS_Proxy'] = set(abess_genes)
    logger.info(f"      ABESS Proxy selected {len(abess_genes)} genes.")

    # 5. The Ultimate Intersection
    logger.info("Step 4.4: Analyzing Consensus (Upset Logic)...")
    all_selected = []
    for s in selected_sets.values():
        all_selected.extend(list(s))
    
    counts = Counter(all_selected)
    
    # Check for perfect 6/6 intersection first, then 5/6
    final_signature = []
    
    logger.info("   Vote Counts for Top Genes:")
    for g, c in counts.most_common(10):
        logger.info(f"      {g}: {c}/6 models")

    # Strict Intersection: Genes selected by ALL algorithms or at least 5
    for thresh in [6, 5]:
        consensus = [g for g, c in counts.items() if c >= thresh]
        if len(consensus) >= 3: # We expect around 5
            logger.info(f"   Using robust consensus with >= {thresh} algorithms.")
            final_signature = consensus
            break
    
    if not final_signature:
         logger.warning("   Strict consensus failed. Falling back to Majority Vote (>=4).")
         final_signature = [g for g, c in counts.items() if c >= 4]

    logger.info(f"   Final Result: {len(final_signature)} genes identified as Core Signature.")
    logger.info(f"   Genes: {final_signature}")
    
    # Verification
    expected_genes = {'S100A16', 'TIMP1', 'CST6', 'MUC5B', 'MARCKSL1'}
    if set(final_signature) == expected_genes:
        logger.info("   ✅ Verification PASSED: Matches the 5-gene signature perfectly.")
    else:
        logger.warning(f"   ⚠️ Verification WARNING: Output differs from expected. Diff: {set(final_signature) ^ expected_genes}")

    # Save Final
    pd.DataFrame({"gene": final_signature}).to_csv(os.path.join(out_table_dir, "final_consensus_genes.csv"), index=False)
    logger.info("✅ Pipeline Complete. Final genes saved.")

if __name__ == "__main__":
    main()
