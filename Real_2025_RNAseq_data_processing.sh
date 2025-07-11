# Bash script for the processing of RNAseq data from the visceral adipose tissue (VAT) 

## Download FASTQ sequence data from NCBI’s Sequence Read Archive under accession no. PRJNA1271077.

# Path to directory with samples
cd /path/to/directory

mkdir VAT/Raw_data


# Assess quality of raw fastq reads

cd VAT

mkdir qc
mkdir qc/fastqc_raw_fastq
mkdir logs
mkdir logs/qc

## QC raw fastq with FastQC

export PATH=/programs/FastQC-0.12.1:$PATH

for file in Raw_data/*.fastq.gz; do
    echo "Processing Sample ${file}";

    fastqc -o qc/fastqc_raw_fastq/ \
    --threads 4 \
    $file
done >> logs/qc/fastqc_raw_fastq_log.txt 2>&1


## MultiQC

export PYTHONPATH=/programs/multiqc-1.15/lib64/python3.9/site-packages:/programs/multiqc-1.15/lib/python3.9/site-packages
export PATH=/programs/multiqc-1.15/bin:$PATH

multiqc -n qc/multiqc_fastq_raw.html --interactive qc/fastqc_raw_fastq/


# Map reads to host genome with STAR

mkdir dbs
mkdir dbs/STAR_Mus_musculus_genome

## Note: We used a previously built STAR reference genome index, based on NCBI's GCF_000001635.27 mouse genome
  # (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/) Date: 2020-06-24
  # GCF_000001635.27_GRCm39_genomic.fna ~2,7 Gb
  # GCF_000001635.27_GRCm39_genomic.gtf ~1 Gb

## The STAR genome index can be built from scratch with the following code:

export PATH=/programs/STAR-2.7.10b:$PATH

STAR --runMode genomeGenerate \
     --genomeDir dbs/STAR_Mus_musculus_genome \
     --genomeFastaFiles /path/to/directory/with/all/GCF_000001635.27_GRCm39_genomic.fna \
     --sjdbGTFfile /path/to/directory/with/GCF_000001635.27_GRCm39_genomic.gtf \
     --runThreadN 20

## Map transcripts onto the reference genome

mkdir map
mkdir logs/map

for file in Raw_data/*.fastq; do
	basename=$(basename "$file" .fastq)
    echo "Processing Sample ${basename}";

    STAR --runThreadN 20 \
		--genomeDir dbs/STAR_Mus_musculus_genome \
		--readFilesIn $file \
		--quantMode GeneCounts \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix map/${basename}/${basename}.
done >> logs/map/STAR_map_log.txt 2>&1


# Assess quality of mapping results

multiqc -d map/ -f --title "Mapping QC" -n qc/multiqc_map.html

## Copy STAR's ReadsPerGene.out.tab files to a shared folder

mkdir map/STAR_ReadsPerGene_untrimmed

for dir in map/*; do
    find "${dir}" -name "*.ReadsPerGene.out.tab" -exec cp {} map/STAR_ReadsPerGene_untrimmed/ \;
done

## Download the STAR_ReadsPerGene_untrimmed folder for both Adr and VAT samples to continue the transcriptomic analyses using R.