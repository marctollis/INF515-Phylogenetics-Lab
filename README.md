# INF515 Comparative Genomics: Phylogenetics Lab

Infer gene trees for eight protein-coding loci, combine them into an ASTRAL species tree, and compare relationships and branch support. Sample primates and use dog (*Canis lupus familiaris*) as the outgroup.

## 1. Set up the project

Run these commands on the course Linux cluster. Use the course-provided interactive Slurm allocation for the manual analysis. Python, MAFFT, Java, `wget`, and `unzip` must be available.

```bash
export PROJ_DIR="/scratch/$USER/CompGenomicsCourse/fall26/phylogenetics"
mkdir -p "$PROJ_DIR"
cd "$PROJ_DIR"

mkdir -p genes species_tree

# NCBI Datasets
wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
chmod +x datasets

# Course scripts
wget https://raw.githubusercontent.com/marctollis/INF515-Phylogenetics-Lab/main/filter_fasta.py
wget https://raw.githubusercontent.com/marctollis/INF515-Phylogenetics-Lab/main/analyze_remaining_genes.sbatch

# IQ-TREE 2.4.0
wget https://github.com/iqtree/iqtree2/releases/download/v2.4.0/iqtree-2.4.0-Linux-intel.tar.gz
tar -xzvf iqtree-2.4.0-Linux-intel.tar.gz

# ASTRAL 5.7.8
wget https://github.com/smirarab/ASTRAL/raw/master/Astral.5.7.8.zip
unzip Astral.5.7.8.zip

module load mafft
java -version
```

Load the course Python and Java modules if required. Python must include any dependencies used by `filter_fasta.py`.

Keep shared software and scripts at the project root. Store gene analyses under `genes/` and the combined species-tree analysis under `species_tree/`:

```text
phylogenetics/
├── datasets
├── filter_fasta.py
├── analyze_remaining_genes.sbatch
├── iqtree-2.4.0-Linux-intel/
├── Astral/
├── genes/
│   ├── brca1/
│   ├── foxp2/
│   ├── egfr/
│   ├── rag1/
│   ├── rbp3/
│   ├── ghr/
│   ├── vwf/
│   └── apob/
└── species_tree/
```

The gene folders will be created as each analysis starts.

## 2. Analyze BRCA1 step by step

### Download coding sequences

Create the BRCA1 folder and work inside it so all outputs stay together.

```bash
cd "$PROJ_DIR"
mkdir -p genes/brca1
cd genes/brca1

"$PROJ_DIR/datasets" download gene symbol brca1 \
    --ortholog primates,9615 \
    --include cds

unzip ncbi_dataset.zip
cp ncbi_dataset/data/cds.fna brca1_mammals_ncbi_datasets.fna

# Inspect headers and count sequences.
grep '^>' brca1_mammals_ncbi_datasets.fna | head
grep -c '^>' brca1_mammals_ncbi_datasets.fna
```

Taxon 9615 identifies dog. Keep the downloaded archive and metadata as a record of the input data.

### Filter sequences

Use the course script to retain one representative transcript per species.

```bash
python "$PROJ_DIR/filter_fasta.py" \
    brca1_mammals_ncbi_datasets.fna \
    brca1_mammals_filtered.fna

grep '^>' brca1_mammals_filtered.fna
grep -c '^>' brca1_mammals_filtered.fna
```

Confirm that the filtered headers include `Canis_lupus_familiaris`. These species labels will also be used in the species-tree analysis. Stop and investigate if filtering fails or the dog outgroup is missing.

### Align sequences

```bash
srun mafft brca1_mammals_filtered.fna \
    > brca1_mammals_aligned.fna \
    2> brca1_mafft.log

less brca1_mammals_aligned.fna
```

Inspect the alignment for extensive gaps or incomplete sequences. Press `q` to exit `less`.

### Infer the gene tree

```bash
srun -c 2 "$PROJ_DIR/iqtree-2.4.0-Linux-intel/bin/iqtree2" \
    -s brca1_mammals_aligned.fna \
    -o Canis_lupus_familiaris \
    -B 1000 \
    -T 2 \
    --prefix brca1

less brca1.iqtree
cat brca1.treefile

cd "$PROJ_DIR"
```

IQ-TREE selects a substitution model by default. The options specify:

- `-o`: dog as the outgroup.
- `-B 1000`: 1,000 ultrafast bootstrap replicates.
- `-T 2`: two threads.
- `--prefix brca1`: output filenames beginning with `brca1`.

Find the selected model in `brca1.iqtree`. The inferred gene tree is stored in `brca1.treefile` in Newick format.

## 3. Analyze the remaining genes with a Slurm job array

The script downloaded in Section 1, `analyze_remaining_genes.sbatch`, analyzes **FOXP2, EGFR, RAG1, RBP3, GHR, VWF, and APOB** using the same workflow as BRCA1.

The script stores the gene names in a Bash array. Slurm creates seven tasks, and each task selects one gene, creates its folder under `genes/`, downloads and filters its sequences, runs MAFFT, and infers its IQ-TREE tree.

| Task index | Gene |
|---|---|
| 0 | FOXP2 |
| 1 | EGFR |
| 2 | RAG1 |
| 3 | RBP3 |
| 4 | GHR |
| 5 | VWF |
| 6 | APOB |

### Inspect and submit the downloaded script

```bash
cd "$PROJ_DIR"
less analyze_remaining_genes.sbatch
```

Identify the gene array, the Slurm array directive, and the commands you ran manually for BRCA1. Press `q` to exit.

Submit from the project root so the script can find the shared software:

```bash
sbatch analyze_remaining_genes.sbatch
squeue -u "$USER"
```

Use any account or partition options specified by your instructor when submitting. Record the job ID printed by `sbatch`. Slurm may run multiple tasks concurrently, depending on available resources.

### Check the results

After the tasks finish, inspect their Slurm logs at the project root. Filenames include the job ID and task index:

```text
slurm-JOBID_0.out
slurm-JOBID_0.err
```

A failed task does not stop the other tasks. Check each task’s completion status and logs.

Confirm that all eight gene trees exist:

```bash
cd "$PROJ_DIR"

for gene in brca1 foxp2 egfr rag1 rbp3 ghr vwf apob; do
    if test -s "genes/$gene/$gene.treefile"; then
        echo "$gene: tree present"
    else
        echo "$gene: tree missing or empty"
    fi
done
```

Inspect the filtered headers, alignment, and IQ-TREE report for each gene. For example:

```bash
grep '^>' genes/foxp2/foxp2_mammals_filtered.fna
less genes/foxp2/foxp2_mammals_aligned.fna
less genes/foxp2/foxp2.iqtree
```

NCBI records change, so sequence counts and species coverage may differ among genes. Confirm consistent labels for shared species and examine incomplete sequences or extensive gaps.

If a task fails, inspect the error and consult your instructor before resubmitting it. Proceed to ASTRAL only after all eight analyses complete successfully.

## 4. Infer the ASTRAL species tree

Combine one best gene tree per locus. Use `.treefile` outputs rather than bootstrap replicate trees.

Run this block within the course-provided Slurm allocation. It stops if any required tree is missing or empty.

```bash
cd "$PROJ_DIR"

(
    set -e

    for gene in brca1 foxp2 egfr rag1 rbp3 ghr vwf apob; do
        test -s "genes/$gene/$gene.treefile"
    done

    for gene in brca1 foxp2 egfr rag1 rbp3 ghr vwf apob; do
        cat "genes/$gene/$gene.treefile"
        printf '\n'
    done > species_tree/gene_trees.txt

    srun java -jar "$PROJ_DIR/Astral/astral.5.7.8.jar" \
        -i species_tree/gene_trees.txt \
        -o species_tree/species_tree.nwk \
        2> species_tree/astral.log
)

cat species_tree/species_tree.nwk
```

ASTRAL treats the input gene trees as unrooted and estimates an unrooted species tree. Its default branch-support values are **local posterior probabilities**, ranging from 0 to 1. IQ-TREE’s gene-tree support values are **ultrafast bootstrap percentages**, ranging from 0 to 100. These are different support measures.

## 5. View, root, and examine trees in Phylo.io

1. Print the species tree:

   ```bash
   cat "$PROJ_DIR/species_tree/species_tree.nwk"
   ```

2. Copy the entire Newick string, including the final semicolon.
3. Open [Phylo.io](https://www.phylo.io/viewer/). Choose **Add Tree**, paste the text, and select **Newick**. If offered, enable **Use internal label for branches (e.g. bootstrap values)** before loading.
4. Find `Canis_lupus_familiaris`. Select its terminal branch and use the rooting/rerooting action to place the root on that branch. Dog should appear as sister to the primate clade.
5. Toggle **Branch Labels/Support** in the display settings to show or hide internal support values. Menu wording may vary by viewer version.
6. Load individual gene trees using the same copy-and-paste process:

   ```bash
   cat "$PROJ_DIR/genes/brca1/brca1.treefile"
   ```

   Substitute other gene names to view their trees. Check the dog outgroup and display branch support. Use the comparison view to examine discordance.

Rooting in the browser changes the displayed tree; it does not overwrite the Newick file on the cluster. Save or export the rooted view if needed.

## 6. Interpret your results

- Which genes agree or disagree about human, chimpanzee, and gorilla relationships?
- Does the species tree recover the expected human–chimpanzee grouping? What support does it have?
- Which species-tree branches have the weakest support?
- How might incomplete lineage sorting, gene-tree estimation error, or missing data explain discordance?
- Why might adding genes improve inference without guaranteeing a particular topology?

Eight loci provide a classroom demonstration rather than a definitive phylogenomic analysis. Base your conclusions on the observed topology, support, and data quality.

## Resources

- [NCBI Datasets: download ortholog data](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genes/download-ortholog-data-package/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [IQ-TREE command reference](https://iqtree.github.io/doc/Command-Reference)
- [ASTRAL documentation](https://github.com/smirarab/ASTRAL/blob/master/README.md)
- [Phylo.io](https://www.phylo.io/)
