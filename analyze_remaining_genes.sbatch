#!/bin/bash
#SBATCH --job-name=inf515-genes
#SBATCH --array=0-6
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err

set -euo pipefail

# Submit this script from the project root.
PROJ_DIR="$SLURM_SUBMIT_DIR"
genes=(foxp2 egfr rag1 rbp3 ghr vwf apob)
gene="${genes[$SLURM_ARRAY_TASK_ID]}"

module load mafft
# Load the course Python and Java modules if required.

mkdir -p "$PROJ_DIR/$gene"
cd "$PROJ_DIR/$gene"

echo "Analyzing $gene"

../datasets download gene symbol "$gene" \
    --ortholog primates,9615 \
    --include cds

unzip ncbi_dataset.zip
cp ncbi_dataset/data/cds.fna \
    "${gene}_mammals_ncbi_datasets.fna"

python ../filter_fasta.py \
    "${gene}_mammals_ncbi_datasets.fna" \
    "${gene}_mammals_filtered.fna"

grep -c '^>' "${gene}_mammals_filtered.fna"

if ! grep -q '^>Canis_lupus_familiaris\([[:space:]]\|$\)' \
    "${gene}_mammals_filtered.fna"; then
    echo "Expected dog outgroup missing: inspect filtered headers." >&2
    exit 1
fi

srun mafft "${gene}_mammals_filtered.fna" \
    > "${gene}_mammals_aligned.fna" \
    2> "${gene}_mafft.log"

srun ../iqtree-2.4.0-Linux-intel/bin/iqtree2 \
    -s "${gene}_mammals_aligned.fna" \
    -o Canis_lupus_familiaris \
    -B 1000 \
    -T "$SLURM_CPUS_PER_TASK" \
    --prefix "$gene"

echo "Finished $gene"
