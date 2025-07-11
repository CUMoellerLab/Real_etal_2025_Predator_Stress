## Bash script for the processing of metagenomic data from wild-derived mouse (Mus musculus domesticus) fecal samples

## Metagenomic data was processed using MAGmaker, our lab's custom made snakemake pipeline (https://github.com/CUMoellerLab/sn-mg-pipeline)
## This pipeline takes fastq files, quality filters them, removes host reads, and assembles reads into contiguous sequences (contigs) 
## These contigs will then be binned into metagenome-assembled genomes (MAGs)


## Please note that timepoints T1=D0, T3=D15, T4=D22, and T5=D30.


## Download FASTQ sequence data from NCBI’s Sequence Read Archive under accession no. PRJNA1271077.

# Path to directory with samples
cd /path/to/working/directory


## Set up

# Download MAGmaker pipeline
git clone https://github.com/CUMoellerLab/sn-mg-pipeline.git
cd sn-mg-pipeline 

# Download mamba to handle package installation
conda install -c conda-forge mamba

# Download snakemake with mamba and create a specific environment
mamba env create -n snakemake -f resources/env/snakemake.yaml

# Activate snakemake environment
conda activate snakemake

# Update samples.txt, units.txt, and config.yaml files in the ./resources/config directory
# See more information here: https://github.com/CUMoellerLab/sn-mg-pipeline

# Start pipeline and install necessary packages
conda install -n base -c conda-forge mamba
snakemake --cores 8 --use-conda

# The first time you run this, it may take longer to set up your conda environment. 
# Be sure to select the appropriate number of cores for your analysis.



## Quality Control

# Create a .txt file that is a list of all sample names, where each sample is a row. 
# Make sure to add an empty line and the end, or the last sample in the list may be missed.

# NOTE:
    # Sample names should start with a LETTER not a NUMBER.
    # Python does not accept numbers as the start of a variable name, and sample names starting with a number may cause issues later on.

# For each sample in your list, snakemake will run all steps needed to output host-filtered reads
cat sample_list.txt | while read -r LINE;
  do
    sample=$LINE;

    snakemake -c all --use-conda --conda-prefix ~/snakemake_envs \
    -k output/qc/host_filter/nonhost/$sample.{R1,R2}.fastq.gz \
    --rerun-incomplete --restart-times 5 -n # Remove dry-run marker (-n) when ready to run pipeline
done


## NOTE: 
    # We only assembled D0 and D30 samples into MAGs, and subsequently mapped 
    # quality controled and host-filtered D15 and D22 samples onto our custom built SGB database.
    # Thus, the following assembly and binning steps were run only on D0 and D30 samples.


## Assembly

# Assemble reads into contigs with MEGAHIT using assemble.smk
# Snakemake will run all steps needed to produce a multiqc.html file acessing assembly quality

snakemake -c all --use-conda --conda-prefix ~/snakemake_envs \
-k output/assemble/multiqc_metaquast/multiqc.html \
--rerun-incomplete --restart-times 5 -n # Remove dry-run marker (-n) when ready to run pipeline

# Alternatively, assembly step can be run using only QUAST for quality assessment, using multiqc_assemble.log as the end point.

snakemake -c all --use-conda --conda-prefix ~/snakemake_envs \
-k output/logs/assemble/multiqc_assemble/multiqc_assemble.log \
--rerun-incomplete --restart-times 5 -n # Remove dry-run marker (-n) when ready to run pipeline


## Bin assemblies into MAGs

# Prototype selection
    # Instead of binning 'all-by-all', where all samples are compared to each other, 
    # our pipeline improves computational efficiency by limiting the comparison to representative sequences (prototypes). 
    # It calculates the sourmash distances between samples, and provides a list of sample groupings (of increasing sizes) 
    # that are representative of the genetic diversity in the dataset.

snakemake -c all --use-conda --conda-prefix ~/snakemake_envs \
-k output/prototype_selection/prototype_selection/selected_prototypes.yaml \
--rerun-incomplete --restart-times 5 -n

# We selected the group containing 12 prototypes:
  12:
  - S_619_T5.R1.fastq.gz
  - S_632_T5.R1.fastq.gz
  - S_631_T5.R1.fastq.gz
  - S_606A_T1.R1.fastq.gz
  - S_594_T1.R1.fastq.gz
  - S_599_T5.R1.fastq.gz
  - S_643_T1.R1.fastq.gz
  - S_628_T1.R1.fastq.gz
  - S_617_T1.R1.fastq.gz
  - S_612_T1.R1.fastq.gz
  - S_642_T1.R1.fastq.gz

# Check prototype quality (number of reads and contigs)
	# We dropped sample 635_T1 from our prototype list because it had a lower than average contig number (~190K).
    # We ran subsequent binning steps with the remaining 11 prototypes.

# Update binning.txt
	# Put an "A" next to selected prototypes in the Read_Groups column
	# Put an "A" next to every other sample (including prototypes) in Contig_Groups colum


# Run binning Snakefile

snakemake -s Snakefile-bin -c all --use-conda --conda-prefix ~/snakemake_envs \
--rerun-incomplete --restart-times 5 -n

# At the end of this step, the pipeline should have produced bins using three different algorithms: CONCOCT, MetaBAT2, and MaxBin
# DASTool will have selected the optimal set of bins from among the three binners.
# In our case, DASTool selected 6,404 bins.


## MAGs Quality Control

# GUNC
  # Python package for detection of chimerism and contamination in prokaryotic genomes.
  # Detects whether there is mis-binning (contigs incorrectly assigned to a bin)

# Please refer to the GUNC GitHub (https://github.com/grp-bork/gunc) for installation information.

## Run GUNC w/ GTDB95 database

singularity exec --bind $PWD --pwd $PWD /programs/gunc-1.0.6/gunc.sif \
gunc run --input_dir output/selected_bins/minimap2/DAS_Tool_Fastas/ \
--db_file resources/db/gunc/gunc_db_gtdb95.dmnd \
--threads 56 --out_dir output/selected_bins/minimap2/GUNC_gtdb_bins/

# Out of the initial 6,404 bins, 4,850 passed GUNC (GTDB95)


# CheckM
    # MAG quality assessment (https://github.com/Ecogenomics/CheckM)
    # Check for completeness and contamination:
    ## Completeness: estimated completeness of genome as determined from the presence/absence of marker genes and the expected collocalization of these genes
    ## Contamination: estimated contamination of genome as determined by the presence of multi-copy marker genes and the expected collocalization of these genes
    ## Strain heterogeneity: estimated strain heterogeneity as determined from the number of multi-copy marker pairs which exceed a specified amino acid identity threshold (default = 90%). 
      ## High strain heterogeneity suggests the majority of reported contamination is from one or more closely related organisms (i.e. potentially the same species), 
      ## while low strain heterogeneity suggests the majority of contamination is from more phylogenetically diverse sources


# Create CheckM environment
conda create -n checkm
conda activate checkm

export PATH=/programs/hmmer/bin:$PATH
export PATH=/programs/prodigal-2.6.3:$PATH
export PATH=/programs/pplacer-Linux-v1.1.alpha19:$PATH
export PATH=/programs/checkm-1.2.2/bin:$PATH
export PYTHONPATH=/programs/checkm-1.2.2/lib/python3.9/site-packages

# Set database directory
export CHECKM_DATA_PATH=/local/workdir/dbs/checkm

# Run CheckM
  # lineage_wf: runs tree, lineage_set, analyze, qa
    # tree: place bins in the reference genome tree
    # lineage_set: infer lineage-specific marker sets for each bin
    # analyze: identify marker genes in bins
    # qa: assess bins for contamination and completeness

# Have CheckM internally call genes with Prodigal (same as GUNC)
checkm lineage_wf \
output/selected_bins/minimap2/DAS_Tool_Fastas \
output/selected_bins/minimap2/CheckM_bins  \
--threads 48 --tmpdir output/selected_bins/minimap2 -x fa

# Run CheckM quality assessment with -f qa.tsv -o 2 --tab_table parameters
checkm qa \
output/selected_bins/minimap2/CheckM_bins/lineage.ms \
output/selected_bins/minimap2/CheckM_bins/ \
--out_format 2 --tab_table --file output/selected_bins/minimap2/CheckM_bins/qa.tsv


# Produce a merged file combining the outputs of GUNC and checkM
  # Both should have been run on the same input files.
  # CheckM qa should be run with -f qa.tsv -o 2 --tab_table parameters.

singularity exec --bind $PWD --pwd $PWD /programs/gunc-1.0.6/gunc.sif \
gunc merge_checkm \
--gunc_file output/selected_bins/minimap2/GUNC_gtdb_bins/GUNC.gtdb_95.maxCSS_level.tsv \
--checkm_file output/selected_bins/minimap2/CheckM_bins/qa.tsv \
--out_dir output/selected_bins/minimap2/


# 3,542/6,404 bins passed both GUNC and high MIMAG (Completion: >90%; Contamination: <5%)

# Move high quality bins to a separate directory

cat High_MIMAG_GUNC_bins_list.txt | while read -r LINE;
  do
    file_name=$LINE;

    cp output/selected_bins/minimap2/DAS_Tool_Fastas/${file_name}.fa output/selected_bins/minimap2/High_MIMAG_GUNC_bins/
done

# Move MAGs to new work directory
cp -r /path/to/directory/sn-mg-pipeline/output/selected_bins/minimap2 /path/to/directory/MAGs

cd /path/to/directory/MAGs


## MAG Annotation

# Annotate high quality MAGs w/ Prokka (https://github.com/tseemann/prokka)

# Run following commands from the directory where the genome fasta file is located. 

cd High_MIMAG_GUNC_bins

mkdir Prokka

# We ran Prokka through a singularity image is built from Docker images staphb/prokka:1.14.5. The image has a build-in db
export PROKKA_CMD="singularity run -C -B $PWD --pwd $PWD /programs/prokka-1.14.5-r9/prokka.sif"

for mag in *.fa; do
  mag_name=$(basename -s .fa $mag)

  $PROKKA_CMD nice -n 10 prokka -out Prokka/$mag_name --prefix $mag_name $mag --cpus 48
done



## Make SGB database

# Dereplicate high quality MAGs into SGBs (95% ANI) w/ dRep (https://drep.readthedocs.io/en/latest/index.html)

cd /path/to/directory/MAGs

mkdir dRep
mkdir dRep/SGBs_95ANI

# Install dRep (Find more information here: https://drep.readthedocs.io/en/latest/installation.html)
conda config --add channels bioconda; conda install drep

# Run dRep
dRep dereplicate dRep/SGBs_95ANI -g High_MIMAG_GUNC_bins/*.fa -p 64 -pa 0.90 -sa 0.95 --genomeInfo dRep/checkM_genome_info.csv

    # -p PROCESSORS threads (default: 6)
    # -pa P_ANI: ANI threshold to form primary (MASH) clusters (default: 0.9)
    # -sa S_ANI: ANI threshold to form secondary clusters (default:0.95)
    # --genomeInfo GENOMEINFO: location of .csv file containing quality information on the genomes.

# dRep made 173 representative SBGs (95% ANI) out of ~3.5k high qual MAGs


# Concatenate representative SGBs into a single fasta file (our custom SGB database)

# Append MAG_ID to contig name
for file in dRep/SGBs_95ANI/dereplicated_genomes/*.fa; do
    filename=$(basename "$file" .fa)

    # Prepend filename to sequence identifiers using awk
    awk -v filename="$filename" '/^>/{print ">" filename "_" substr($0, 2); next} 1' "$file" > dRep/SGBs_95ANI/representative_SGBs/${filename}.fa
done

# Concatenate all files into a single file
cat dRep/SGBs_95ANI/representative_SGBs/*.fa > Rep_SGBs_95ANI.fasta

mv MAGs/dRep/SGBs_95ANI/dereplicated_genomes/Rep_SGBs_95ANI.fasta MAGs/Rep_SGBs_95ANI.fasta


## Calculate genome size of each representative SGBs
    # This information is needed to calculate the abundance profile of each representative SGB in our samples

# Write the header to the output file
echo "genome,genome_size" > "dRep/SGBs_95ANI/rep_SGBs_genome_sizes.csv"

for fasta in dRep/SGBs_95ANI/representative_SGBs/*.fa; do
  # Get the base name of the FASTA file (without path and extension)
  genome_name=$(basename "$fasta" .fa)

  # Calculate the total size by summing up the lengths of all sequence lines
  genome_size=$(grep -v '^>' "$fasta" | tr -d '\n' | wc -c)

  # Append the result to the output file
  echo "$genome_name,$genome_size" >> "dRep/SGBs_95ANI/rep_SGBs_genome_sizes.csv"
done



## Taxonomically classify SGB database 

# We will use the GTDB-Tk database (https://github.com/Ecogenomics/GTDBTk) 

# Set up
cd /path/to/directory/MAGs
mkdir GTDB
mkdir GTDB/SGBs_95ANI

export OMP_NUM_THREADS=64
export PYTHONPATH=/programs/gtdbtk-2.4.0/lib64/python3.9/site-packages:/programs/gtdbtk-2.4.0/lib/python3.9/site-packages
export PATH=/programs/gtdbtk-2.4.0/bin:/programs/hmmer/binaries:/programs/prodigal-2.6.3:/programs/FastTree-2.1.11:/programs/fastANI-1.32:/programs/pplacer-Linux-v1.1.alpha19:/programs/mash-Linux64-v2.3:/programs/skani-0.2.1:$PATH

# Download latest GTDB release (R220)
mkdir GTDB/release220
cd GTDB/release220

wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
tar xvfz gtdbtk_data.tar.gz

#Set GTDBTK_DATA_PATH to reference db directory
export GTDBTK_DATA_PATH=/path/to/directory/MAGs/GTDB/release220


# Run classify_wf 
    # (https://ecogenomics.github.io/GTDBTk/commands/classify_wf.html)
    # The classify workflow consists of four steps: ani_screen, identify, align, and classify.
    # ani_screen: compares user genomes against a Mash database composed of all GTDB representative genomes, then verify the best mash hits using skani. 
        # User genomes classified with skani are not run through the rest of the pipeline (identify, align, classify) and are reported in the summary file.
    # identify: calls genes using Prodigal, and uses HMM models and the HMMER package to identify the 120 bacterial and 53 archaeal marker genes used for phylogenetic inference (Parks et al., 2018). 
        # Multiple sequence alignments (MSA) are obtained by aligning marker genes to their respective HMM model.
    # align: concatenates the aligned marker genes and filters the concatenated MSA to approximately 5,000 amino acids.
    # classify: uses pplacer to find the maximum-likelihood placement of each genome in the GTDB-Tk reference tree. 
        # GTDB-Tk classifies each genome based on its placement in the reference tree, its relative evolutionary divergence, and/or average nucleotide identity (ANI) to reference genomes.

gtdbtk classify_wf \
    --genome_dir dRep/SGBs_95ANI/dereplicated_genomes \
    --out_dir GTDB/SGBs_95ANI \
    --extension fa \
    --cpus 64 \
    --skip_ani_screen



## Identify genes in the SGB database

# Gene calling w/ Prodigal (https://github.com/hyattpd/Prodigal) 
# using a wrapper for parallelization (https://github.com/sjaenick/pprodigal)

mkdir Prodigal
mkdir Prodigal/SGBs_95ANI

export PATH=/programs/prodigal-2.6.3:$PATH

# Install
conda create -n mypprodigal -c conda-forge -c bioconda pprodigal
conda activate mypprodigal

pprodigal \
    --input "Rep_SGBs_95ANI.fasta" \
    --nucl "Prodigal/SGBs_95ANI/Rep_SGBs_95ANI.fna" \
    --proteins "Prodigal/SGBs_95ANI/Rep_SGBs_95ANI.faa" \
    --output "Prodigal/SGBs_95ANI/Rep_SGBs_95ANI.gff" \
    --tasks 64




## Profile representative SGBs abundance in all samples

# Map all fasta files (D0, D15, D22, D30) to the Rep_SGBs_95ANI.fasta file with Bowtie2

export PATH=/programs/bowtie2-2.5.1-linux-x86_64:$PATH

# Make directory to save Bowtie2 outputs
mkdir Bowtie2
mkdir Bowtie2/index
mkdir Bowtie2/SAM
mkdir Bowtie2/BAM


# Create Bowtie2 index for the representative FASTA file
bowtie2-build Rep_SGBs_95ANI.fasta Bowtie2/index/Rep_SGBs_95ANI

# Map quality-controlled and host-filtered reads to Rep_SGBs_95ANI.fasta

dir_path=/path/to/directory/with/quality/controled/and/host/filtered/reads

# Loop through all FASTA files
for fasta_file in $dir_path/*.gz; do
    base_name=$(basename "$fasta_file")
    sample=$(echo "$base_name" | sed 's/\.R[12]\.fastq\.gz//')
    echo $(date +"%Y-%m-%d %H:%M:%S");
    echo "Processing Sample ${sample}";

    R1="${sample}.R1.fastq.gz";
    R2="${sample}.R2.fastq.gz";
    
    # Run Bowtie2 alignment
    nice -n 10 bowtie2 -x Bowtie2/index/Rep_SGBs_95ANI \
        -1 ${dir_path}/${R1} -2 ${dir_path}/${R2} \
        -S Bowtie2/SAM/${sample}.sam \
        --threads 64
done >> Bowtie2/mapping.log.txt 2>&1



# Convert SAM files to BAM, sort BAM files, index them, and get mapping statistics
    #!/bin/bash
    sam_file=$1
    sample=$(basename "$sam_file" .sam)
    wd="/path/to/directory/MAGs"

    # Convert SAM to BAM and sort it
    samtools view -bS "${sam_file}" | samtools sort -o $wd/Bowtie2/BAM/"${sample}".sorted.bam

    # Index the BAM file
    samtools index $wd/Bowtie2/BAM/"${sample}".sorted.bam

    # Get mapping statistics and save to a text file
    samtools flagstat $wd/Bowtie2/BAM/"${sample}".sorted.bam > $wd/Bowtie2/Mapping_stats/"${sample}".mapping_stats.txt

chmod +x Bowtie2/sam_to_bam.sh

find /path/to/directory/MAGs/Bowtie2/SAM/ -name "*.sam" | parallel -j 32 --verbose Bowtie2/sam_to_bam.sh


# Compile mapping statistics into a summary table
    # This provides information on the mappability of our samples to our custom SGB database

#!/bin/bash

# Output file
output="mapping_summary.tsv"

# Header for the output file
echo -e "Sample\tTotal_Reads\tMapped_Reads\tMapped_Percentage\tProperly_Paired\tPaired_Percentage\tSingletons" > "$output"

# Loop through all .mapping_stats.txt files
for file in Bowtie2/Mapping_stats/*.mapping_stats.txt; do
    # Extract the sample name
    sample=$(basename "$file" .mapping_stats.txt)

    # Use grep and awk to extract statistics, ensuring no duplication
    total_reads=$(grep "in total" "$file" | awk '{print $1}')
    mapped_reads=$(grep "mapped (" "$file" | head -n 1 | awk '{print $1}')
    mapped_percentage=$(grep "mapped (" "$file" | head -n 1 | awk -F '[()%]' '{print $2}')
    properly_paired=$(grep "properly paired" "$file" | awk '{print $1}')
    paired_percentage=$(grep "properly paired" "$file" | awk -F '[()%]' '{print $2}')
    singletons=$(grep "singletons" "$file" | awk '{print $1}')

    # Append the extracted data to the output file
    echo -e "${sample}\t${total_reads}\t${mapped_reads}\t${mapped_percentage}\t${properly_paired}\t${paired_percentage}\t${singletons}" >> "$output"
done

echo "Mapping summary saved to $output"

####

chmod +x Bowtie2/extract_mapping_stats.sh

Bowtie2/extract_mapping_stats.sh



# Create a scaffold_to_bin file 
source /programs/miniconda3/bin/activate drep

chmod +x dRep/parse_stb.py

dRep/parse_stb.py --reverse -f dRep/SGBs_95ANI/representative_SGBs/* -o dRep/SGBs_95ANI/scaffold_to_bin.stb


# Profile mapping with InStrain profile (https://instrain.readthedocs.io/en/latest/user_manual.html)

export PYTHONPATH=/programs/instrain-1.8.0/lib/python3.9/site-packages:/programs/instrain-1.8.0/lib64/python3.9/site-packages
export PATH=/programs/instrain-1.8.0/bin:$PATH

mkdir InStrain
mkdir InStrain/Profile

#!/bin/bash
    file=$1
    sample=$(basename "$file" .sorted.bam)

    nice -n 10 inStrain profile $file Rep_SGBs_95ANI.fasta \
        -o InStrain/Profile/${sample} \
        --stb dRep/SGBs_95ANI/scaffold_to_bin.stb \
        --gene_file Prodigal/SGBs_95ANI/Rep_SGBs_95ANI.fna \
        -p 6 \
        --detailed_mapping_info \
        --database_mode

chmod +x InStrain/inStrain_profile.sh

find Bowtie2/BAM/ -name "*.bam" | parallel -j 10 --verbose InStrain/inStrain_profile.sh



# Copy all *_genome_info.tsv outputs to common directory
    # These are the files that will allow us to calculate the abundance of representative SGBs in each sample

mkdir InStrain/SGB_profiles

for file in InStrain/Profile/*/output/*_genome_info.tsv; do
    if [ -f "$file" ]; then
        cp "$file" InStrain/SGB_profiles/
    fi
done


# We now have all necessary files to get an abundance table of each of our representative SGBs across all samples.
# Please download the following files and refer to Stress_analyses.Rmd for the next steps in the analyses.

# Files needed for subsequent analyses:
    # MAGs/Bowtie2/mapping_summary.tsv
    # output/qc/multiqc_data/multiqc_fastqc.txt
    # MAGs/dRep/SGBs_95ANI/rep_SGBs_genome_sizes.csv
    # MAGs/GTDB/SGBs_95ANI/classify/gtdbtk.bac120.summary.tsv
    # MAGs/InStrain/SGB_profiles

