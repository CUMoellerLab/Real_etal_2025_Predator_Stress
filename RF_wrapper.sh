#!/bin/bash
# RF_run_all.sh
# Wrapper to run phenotype prediction for microbiome and transcriptome data
# Loops over predictors, analysis types, timepoints, and phenotype types
# Generates organized output directories

# ------------------------------
# CONFIGURATION
# ------------------------------

OUTPUT_BASE="output"

# Input files
ABUNDANCE_FILE="data/rep_SGB_raref_relative_abundance_MAG_ID.txt"
TRANSCRIPTOME_FILE="data/VAT_normalized_counts.txt"
BEHAVIOR_FILE="data/Behavior_data.txt"
VAT_FILE="data/VAT_data.txt"

# Ensure base output directory exists
mkdir -p "$OUTPUT_BASE"

# ------------------------------
# LOOP OVER PREDICTORS AND ANALYSES
# ------------------------------

for PREDICTOR in microbiome transcriptome; do
  for ANALYSIS in across within; do

    if [[ "$PREDICTOR" == "microbiome" ]]; then
      # Microbiome: loop over timepoints and phenotype types
      for TP in T1 T3 T4 T5; do
        for PHENOTYPE_TYPE in Behavior VAT; do

          # Set phenotype variables and file based on type
          if [[ "$PHENOTYPE_TYPE" == "Behavior" ]]; then
            PHENOTYPE_VARS="Center_occupancy Grooming_duration Social_preference"
            PHENOTYPE_FILE="$BEHAVIOR_FILE"
          else
            PHENOTYPE_VARS="PC1 PC2 PC3 PC4"
            PHENOTYPE_FILE="$VAT_FILE"
          fi

          # Define output directory
          OUTPUT_DIR="$OUTPUT_BASE/$PREDICTOR/$ANALYSIS/$PHENOTYPE_TYPE/$TP"
          mkdir -p "$OUTPUT_DIR"

          # Run Python script
          python RF_predict_phenotypes.py \
            --predictor "$PREDICTOR" \
            --analysis "$ANALYSIS" \
            --phenotype_file "$PHENOTYPE_FILE" \
            --abundance_file "$ABUNDANCE_FILE" \
            --output_dir "$OUTPUT_DIR" \
            --phenotype_vars $PHENOTYPE_VARS \
            --timepoint "$TP"

        done
      done

    else
      # Transcriptome: only Behavior phenotype
      PHENOTYPE_VARS="Center_occupancy Grooming_duration Social_preference"
      OUTPUT_DIR="$OUTPUT_BASE/$PREDICTOR/$ANALYSIS"
      mkdir -p "$OUTPUT_DIR"

      python RF_predict_phenotypes.py \
        --predictor "$PREDICTOR" \
        --analysis "$ANALYSIS" \
        --phenotype_file "$BEHAVIOR_FILE" \
        --transcriptome_file "$TRANSCRIPTOME_FILE" \
        --output_dir "$OUTPUT_DIR" \
        --phenotype_vars $PHENOTYPE_VARS

    fi

  done
done

echo "All analyses completed. Results saved in $OUTPUT_BASE."