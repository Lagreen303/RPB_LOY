#second iteration of LOY machine learning model using top highly variable genes instead of PCs
#First section is just the data exploration and QC steps
import pandas as pd
import scanpy as sc
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, FunctionTransformer
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score, roc_auc_score

###Preprocessing
#inital data exploration
adata = sc.read_h5ad("AIDA_LOY_status.h5ad")
#filter to only monocytes THEN normalise + find highly varibale genes
adata = adata[adata.obs['author_cell_type'].isin(['Monocyte', 'CD14+_Monocyte', 'CD16_Monocyte'])]

# Visualize the top 20 highly expressed genes and save to the specified directory
sc.pl.highest_expr_genes(adata, n_top=20, save="_monocytes.png")

# Filter cells and genes
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Annotate the group of mitochondrial genes as "mt"
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

#visualize the QC metrics 
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    save="qc_violin_monocytes.png",
)

adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

#normalise and logarithmise the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, save="_monocytes.png")

adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])

sc.pp.scale(adata, max_value=10)

#save the processed data
adata.write("AIDA_LOY_status_monocytes.h5ad")

# PCA
sc.tl.pca(adata, svd_solver="arpack")
sc.pl.pca(adata, save="_monocytes.png")  # Save PCA plot
sc.pl.pca_variance_ratio(adata, log=True, save="_monocytes.png")  # Save variance ratio plot

## random forest classifier for LOY using top highly variable genes ##

# Extract features
X = adata.X
y_loy = adata.obs['LOY_status'].values 
print("Unique values and their counts:", adata.obs['LOY_status'].value_counts().to_dict())

# Split data into train/test sets
X_train, X_test, y_train, y_test = train_test_split(
    X, y_loy, test_size=0.2, stratify=y_loy, random_state=42
)

# Train classifier
clf = RandomForestClassifier(n_estimators=500, class_weight='balanced')
clf.fit(X_train, y_train)

# Evaluate
y_pred = clf.predict(X_test)
y_proba = clf.predict_proba(X_test)[:, 1]
print("Classification Report:\n", classification_report(y_test, y_pred))
print(f"Accuracy: {accuracy_score(y_test, y_pred):.2f}")
print(f"AUC-ROC: {roc_auc_score(y_test, y_proba):.2f}")

#Unique values and their counts: {0.0: 33014, 1.0: 1402}
#Classification Report:
#               precision    recall  f1-score   support

#         0.0       0.96      1.00      0.98      6604
#         1.0       0.00      0.00      0.00       280

#    accuracy                           0.96      6884
#   macro avg       0.48      0.50      0.49      6884
#weighted avg       0.92      0.96      0.94      6884

#Accuracy: 0.96
#AUC-ROC: 0.76