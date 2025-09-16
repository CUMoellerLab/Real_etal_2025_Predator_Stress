#!/usr/bin/env python3
"""
predict_phenotypes.py

Unified script for predicting host phenotypes from microbiome or transcriptome data
using Random Forest regression.

Supports:
- Predictor types: microbiome, transcriptome
- Analyses: across (no residualizing), within (residualizing for treatment)
- Cross-validation with permutation testing
- Feature importance
- Predicted vs Actual and Residual plots

Author: Madalena VF Real
Date: 2025-09-16
"""

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for plotting
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold, cross_val_predict, permutation_test_score

# ------------------------------
# FUNCTIONS
# ------------------------------

def parse_args():
    """
    Parse command-line arguments to configure the analysis.
    """
    parser = argparse.ArgumentParser(description="Predict phenotypes using Random Forest regression")
    parser.add_argument("--predictor", required=True, choices=["microbiome", "transcriptome"], help="Predictor type")
    parser.add_argument("--analysis", required=True, choices=["across", "within"], help="Across or within-treatment analysis")
    parser.add_argument("--phenotype_file", required=True, help="Behavior or VAT metadata file")
    parser.add_argument("--abundance_file", help="Microbiome abundance file (required for microbiome)")
    parser.add_argument("--transcriptome_file", help="Transcriptome file (required for transcriptome)")
    parser.add_argument("--output_dir", required=True, help="Directory to save outputs")
    parser.add_argument("--phenotype_vars", type=str, nargs="+", help="Phenotype variables to predict")
    parser.add_argument("--timepoint", help="Timepoint for microbiome data (e.g., D0, D15)")
    return parser.parse_args()


def load_data(args):
    """
    Load phenotype metadata and predictor data, merge into a single dataframe.
    
    For microbiome:
        - Requires a timepoint to create unique Sample_IDs
        - Loads abundance table and merges with metadata
    For transcriptome:
        - Uses Mouse_ID as Sample_ID
        - Loads transcriptome table and merges with metadata
    """
    # Load phenotype metadata
    phenotype = pd.read_csv(args.phenotype_file, sep="\t")
    
    if args.predictor == "microbiome":
        # Ensure timepoint is provided
        if args.timepoint is None:
            raise ValueError("timepoint must be specified for microbiome analysis")
        
        # Create Sample_ID by combining Mouse_ID and timepoint
        phenotype["Sample_ID"] = phenotype["Mouse_ID"].astype(str) + "_" + args.timepoint
        
        # Load microbiome abundance table
        abundance = pd.read_csv(args.abundance_file, sep="\t", index_col=0).T.reset_index()
        abundance = abundance.rename(columns={"index": "Sample_ID"})
        # Ensure all abundance values are numeric
        abundance.iloc[:, 1:] = abundance.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')

        # Merge phenotype and abundance
        merged = pd.merge(phenotype, abundance, on="Sample_ID")
        
    else:
        # For transcriptome, Sample_ID = Mouse_ID
        phenotype["Sample_ID"] = phenotype["Mouse_ID"].astype(str)

        transcriptome = pd.read_csv(args.transcriptome_file, sep="\t", index_col=0).T.reset_index()
        transcriptome = transcriptome.rename(columns={"index": "Mouse_ID"})
        transcriptome.iloc[:, 1:] = transcriptome.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')
        transcriptome["Sample_ID"] = transcriptome["Mouse_ID"].astype(str)

        # Merge phenotype and transcriptome
        merged = pd.merge(phenotype, transcriptome, on="Sample_ID")
    
    return merged


def residualize_within_treatment(df, phenotype_vars):
    """
    Residualize phenotype values for treatment effects using linear regression.
    
    This is used for "within" analyses to remove treatment effects.
    """
    df_resid = df.copy()
    for var in phenotype_vars:
        y = df[var]
        # Create dummy variables for Treatment
        X = pd.get_dummies(df["Treatment"], drop_first=True)
        # Fit linear model using least squares
        model = np.linalg.lstsq(X, y, rcond=None)
        coef = model[0]
        predicted = X.values @ coef
        # Replace original phenotype with residuals
        df_resid[var] = y - predicted
    return df_resid


def run_random_forest(X, y, cv_splits=5, n_permutations=1000):
    """
    Run Random Forest regression with cross-validation and permutation test.
    
    Returns:
        - Cross-validated predictions
        - R2 score from permutation test
        - P-value from permutation test
        - Feature importances
    """
    rf_model = RandomForestRegressor(n_estimators=20, random_state=42, n_jobs=-1)
    
    # Cross-validated predictions
    y_pred_cv = cross_val_predict(rf_model, X, y, cv=KFold(n_splits=cv_splits, shuffle=True, random_state=42))
    
    # Permutation test to assess significance of model performance
    score, permutation_scores, pvalue = permutation_test_score(
        rf_model, X, y, cv=KFold(n_splits=cv_splits, shuffle=True, random_state=42),
        n_permutations=n_permutations, scoring="r2", random_state=42
    )
    
    # Fit on full data to get feature importances
    rf_model.fit(X, y)
    feature_importances = pd.DataFrame({
        "Feature": X.columns,
        "Importance": rf_model.feature_importances_
    }).sort_values(by="Importance", ascending=False)
    
    return y_pred_cv, score, pvalue, feature_importances


def plot_results(pred_df, feature_imp, output_dir, base_filename):
    """
    Generate plots for predicted vs actual values, residuals, and top features.
    
    Files are saved as PDF with descriptive filenames.
    """
    # Define treatment colors
    treatment_color_map = {
        "Pair_H2O": "#63B8FF",     # steelblue2
        "Pair_TMT": "#EEC900",     # gold2
        "Single_H2O": "#CD1076",   # deeppink3
        "Single_TMT": "#EE9A49"    # sienna2
    }

    # Predicted vs Actual
    plt.figure(figsize=(4,4))
    sns.scatterplot(data=pred_df, x="Predicted", y="Actual", hue="Treatment", palette=treatment_color_map, alpha=0.6)
    sns.regplot(data=pred_df, x="Predicted", y="Actual", scatter=False, color="red")
    slope, intercept, r_value, p_val, std_err = linregress(pred_df["Predicted"], pred_df["Actual"])
    plt.text(pred_df["Predicted"].min(), pred_df["Actual"].max()*0.95, f"p = {p_val:.3g}")
    plt.xlabel(f"Predicted {pred_df['Phenotype'].iloc[0]}")
    plt.ylabel(f"Actual {pred_df['Phenotype'].iloc[0]}")
    plt.title(f"{pred_df['Phenotype'].iloc[0]} Prediction")
    plt.savefig(os.path.join(output_dir, f"Predicted_vs_actual_{base_filename}.pdf"))
    plt.close()
    
    # Residual plot
    residuals = pred_df["Actual"] - pred_df["Predicted"]
    plt.figure(figsize=(4,4))
    plt.scatter(pred_df["Predicted"], residuals, alpha=0.6, color="green")
    plt.axhline(0, color="red", linestyle="--")
    plt.xlabel(f"Predicted {pred_df['Phenotype'].iloc[0]}")
    plt.ylabel("Residuals")
    plt.title(f"{pred_df['Phenotype'].iloc[0]} Residuals")
    plt.savefig(os.path.join(output_dir, f"Residuals_{base_filename}.pdf"))
    plt.close()
    
    # Feature importance
    plt.figure(figsize=(6,4))
    sns.barplot(data=feature_imp.head(10), x="Importance", y="Feature", palette="viridis")
    plt.xlabel("Feature Importance")
    plt.ylabel("Feature")
    plt.title(f"Top 10 Features")
    plt.savefig(os.path.join(output_dir, f"Feature_importance_{base_filename}.pdf"))
    plt.close()


# ------------------------------
# MAIN SCRIPT
# ------------------------------

def main():
    """
    Main function to orchestrate phenotype prediction analysis.
    
    Steps:
    1. Parse arguments and create output directory
    2. Load and merge phenotype and predictor data
    3. Residualize phenotypes if "within" analysis
    4. Loop through phenotype variables:
       - Prepare predictors and outcome
       - Run Random Forest regression
       - Store results, feature importances, and predicted vs actual values
       - Generate plots
    5. Save all results as tab-separated files with descriptive filenames
    """
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load and merge data
    merged = load_data(args)
    
    # Residualize if performing "within" analysis
    if args.analysis == "within":
        merged = residualize_within_treatment(merged, args.phenotype_vars)
    
    # Initialize storage for outputs
    results_list = []
    feature_list = []
    plot_records = []

    # Loop through each phenotype variable
    for phenotype_var in args.phenotype_vars:
        print(f"Processing {phenotype_var}...")
        
        # Select predictor columns (exclude metadata and phenotype columns)
        phenotype_cols = {"Mouse_ID", "Treatment", "Sample_ID"} | set(args.phenotype_vars)
        X = merged[[c for c in merged.columns if c not in phenotype_cols]]
        y = merged[phenotype_var]
        
        # Run Random Forest regression with cross-validation
        y_pred_cv, score, pval, feature_imp = run_random_forest(X, y)
        
        # Store summary results
        results_list.append({"Phenotype": phenotype_var, "Perm_test_R2": score, "P_value": pval})
        feature_imp["Phenotype"] = phenotype_var
        feature_list.append(feature_imp)
        
        # Store predicted vs actual dataframe
        plot_df = pd.DataFrame({
            "Phenotype": phenotype_var,
            "Treatment": merged["Treatment"],
            "Predicted": y_pred_cv,
            "Actual": y
        })
        plot_records.append(plot_df)
        
        # Construct descriptive base filename for plots
        base_filename = f"{args.predictor}_{args.analysis}_{phenotype_var}"
        if args.predictor == "microbiome":
            base_filename += f"_{args.timepoint}"
        
        # Generate plots
        plot_results(plot_df, feature_imp, args.output_dir, base_filename)
    
    # Save all results as tab-separated files
    results_df = pd.DataFrame(results_list)
    features_df = pd.concat(feature_list, ignore_index=True)
    plots_df = pd.concat(plot_records, ignore_index=True)
    
    results_file = os.path.join(args.output_dir, f"RF_permutation_test_results_{args.predictor}_{args.analysis}.tsv")
    features_file = os.path.join(args.output_dir, f"RF_feature_importances_{args.predictor}_{args.analysis}.tsv")
    plots_file = os.path.join(args.output_dir, f"RF_predicted_vs_actual_{args.predictor}_{args.analysis}.tsv")
    
    results_df.to_csv(results_file, sep="\t", index=False)
    features_df.to_csv(features_file, sep="\t", index=False)
    plots_df.to_csv(plots_file, sep="\t", index=False)
    
    print(f"Analysis complete. Results saved to {args.output_dir}")


if __name__ == "__main__":
    main()
