#test random forest classifier for LOY using 50 PCs and 100 trees
#results were First LOY classification report:
 #              precision    recall  f1-score   support

 #          0       0.99      1.00      0.99     50477
 #          1       0.73      0.13      0.22       850

  #  accuracy                           0.98     51327
  # macro avg       0.86      0.56      0.60     51327
#weighted avg       0.98      0.98      0.98     51327

import scanpy as sc
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, FunctionTransformer
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score, roc_auc_score

# Load .h5ad file
adata = sc.read_h5ad("AIDA_male_LOY.h5ad")

# Step 3: Extract data
X = adata.X  
y_loy = adata.obs['LOY_status'].values 
print("Unique values and their counts:", adata.obs['LOY_status'].value_counts().to_dict())

# Step 4: Preprocess
# Log-transform 
X_log = FunctionTransformer(np.log1p).fit_transform(X)

# Scale features (sparse-safe)
X_scaled = StandardScaler(with_mean=False).fit_transform(X_log)

# PCA (convert to dense)
X_pca = PCA(n_components=50).fit_transform(X_scaled.toarray())

# Step 5: Predict LOY
# Split data into train/test sets
X_train, X_test, y_train, y_test = train_test_split(
    X_pca, y_loy, test_size=0.2, stratify=y_loy, random_state=42
)

# Train classifier
clf = RandomForestClassifier(n_estimators=100, class_weight='balanced', random_state=42)
clf.fit(X_train, y_train)

gclf = HistGradientBoostingClassifier(n_estimators=100, random_state=42)
gclf.fit(X_train, y_train)

# Evaluate
y_pred = clf.predict(X_test)
y_proba = clf.predict_proba(X_test)[:, 1]
print("Classification Report:\n", classification_report(y_test, y_pred))
print(f"Accuracy: {accuracy_score(y_test, y_pred):.2f}")
print(f"AUC-ROC: {roc_auc_score(y_test, y_proba):.2f}")