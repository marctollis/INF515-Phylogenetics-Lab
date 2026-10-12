#!/usr/bin/env bash
# Run from the project root after completing the manual BRCA1 analysis.
set -euo pipefail

PROJ_DIR="$(pwd)"
genes=(foxp2 egfr rag1 rbp3 ghr vwf apob)

module load mafft

for gene in "${genes[@]}"; do
    echo "Analyzing $gene"
    mkdir -p "$PROJ_DIR/$gene"

    (
        cd "$PROJ_DIR/$gene"

        # Download primate and dog ortholog coding sequences.
        ../datasets download gene symbol "$gene" \
            --ortholog primates,9615 \
            --include cds

        unzip ncbi_dataset.zip
        cp ncbi_dataset/data/cds.fna \
            "${gene}_mammals_ncbi_datasets.fna"

        # Filter to one representative transcript per species.
        python ../filter_fasta.py \
            "${gene}_mammals_ncbi_datasets.fna" \
            "${gene}_mammals_filtered.fna"

        echo "Filtered sequence count:"
        grep -c '^>' "${gene}_mammals_filtered.fna"

        # Confirm that the expected dog outgroup is present.
        if ! grep -q '^>Canis_lupus_familiaris\([[:space:]]\|$\)' \
            "${gene}_mammals_filtered.fna"; then
            echo "Expected dog outgroup missing: inspect $gene headers." >&2
            exit 1
        fi

        # Align sequences.
        srun mafft "${gene}_mammals_filtered.fna" \
            > "${gene}_mammals_aligned.fna" \
            2> "${gene}_mafft.log"

        # Infer the gene tree with ultrafast bootstrap support.
        srun -c 2 ../iqtree-2.4.0-Linux-intel/bin/iqtree2 \
            -s "${gene}_mammals_aligned.fna" \
            -o Canis_lupus_familiaris \
            -B 1000 \
            -T 2 \
            --prefix "$gene"
    )

    echo "Finished $gene"
done

echo "All seven additional gene analyses are complete."
