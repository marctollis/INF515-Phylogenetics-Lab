# INF515 Comparative Genomics: Phylogenetics Lab

Infer eight gene trees, combine them into an ASTRAL species tree, and examine gene-tree discordance and branch support. Sample primates and use dog (*Canis lupus familiaris*) as the outgroup.

## 1. Set up your project

Run these commands on the course Linux cluster. Use the course-provided Slurm allocation for computational steps. Python, MAFFT, Java, `wget`, and `unzip` must be available.

```bash
export PROJ_DIR="/scratch/$USER/CompGenomicsCourse/fall26/phylogenetics"
mkdir -p "$PROJ_DIR"
cd "$PROJ_DIR"

# Shared software and filtering script
wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
chmod +x datasets

wget https://raw.githubusercontent.com/marctollis/INF515-Phylogenetics-Lab/main/filter_fasta.py

wget https://github.com/iqtree/iqtree2/releases/download/v2.4.0/iqtree-2.4.0-Linux-intel.tar.gz
tar -xzvf iqtree-2.4.0-Linux-intel.tar.gz

wget https://github.com/smirarab/ASTRAL/raw/master/Astral.5.7.8.zip
unzip Astral.5.7.8.zip

module load mafft
java -version
mkdir -p species_tree
```

Keep shared software and scripts at the project root. Store each gene’s sequences, alignment, tree, and logs in its own folder:

```text
phylogenetics/
├── datasets
├── filter_fasta.py
├── iqtree-2.4.0-Linux-intel/
├── Astral/
├── brca1/
├── foxp2/
├── egfr/
├── rag1/
├── rbp3/
├── ghr/
├── vwf/
├── apob/
└── species_tree/
```

## 2. Analyze BRCA1

Run the following workflow from the project root:

```bash
cd "$PROJ_DIR"
gene=brca1
mkdir -p "$gene"
cd "$gene"

# Download primate and dog ortholog coding sequences.
../datasets download gene symbol "$gene" \
    --ortholog primates,9615 \
    --include cds

unzip ncbi_dataset.zip
cp ncbi_dataset/data/cds.fna "${gene}_mammals_ncbi_datasets.fna"

# Inspect sequence headers and count records.
grep '^>' "${gene}_mammals_ncbi_datasets.fna" | head
grep -c '^>' "${gene}_mammals_ncbi_datasets.fna"

# Retain one representative transcript per species.
python ../filter_fasta.py \
    "${gene}_mammals_ncbi_datasets.fna" \
    "${gene}_mammals_filtered.fna"

grep '^>' "${gene}_mammals_filtered.fna"
grep -c '^>' "${gene}_mammals_filtered.fna"
```

Check that the filtered sequences include `Canis_lupus_familiaris` and that headers use species names consistently. Stop and investigate if filtering fails or the dog outgroup is missing.

Align the sequences and infer the gene tree:

```bash
srun mafft "${gene}_mammals_filtered.fna" \
    > "${gene}_mammals_aligned.fna" \
    2> "${gene}_mafft.log"

less "${gene}_mammals_aligned.fna"

srun -c 2 ../iqtree-2.4.0-Linux-intel/bin/iqtree2 \
    -s "${gene}_mammals_aligned.fna" \
    -o Canis_lupus_familiaris \
    -B 1000 \
    -T 2 \
    --prefix "$gene"

less "${gene}.iqtree"
cat "${gene}.treefile"

cd "$PROJ_DIR"
```

Press `q` to exit `less`. Inspect the alignment for extensive gaps or incomplete sequences. In the IQ-TREE report, identify the selected substitution model and examine support.

IQ-TREE selects a substitution model by default. `-B 1000` requests 1,000 ultrafast bootstrap replicates, `-T 2` uses two threads, and `-o` specifies the dog outgroup. The `.treefile` contains the inferred tree in Newick format.

## 3. Analyze the remaining genes

Repeat **all commands in Section 2** for each gene below, replacing `gene=brca1` with the indicated assignment. Begin each analysis with `cd "$PROJ_DIR"`.

| Gene | Assignment | Analysis folder |
|---|---|---|
| FOXP2 | `gene=foxp2` | `foxp2/` |
| EGFR | `gene=egfr` | `egfr/` |
| RAG1 | `gene=rag1` | `rag1/` |
| RBP3 | `gene=rbp3` | `rbp3/` |
| GHR | `gene=ghr` | `ghr/` |
| VWF | `gene=vwf` | `vwf/` |
| APOB | `gene=apob` | `apob/` |

For every gene, download CDS sequences, filter, inspect the headers, align, and run IQ-TREE. Confirm successful completion before moving to the next gene.

NCBI records change, so species coverage and sequence counts may differ among genes. A species present in multiple gene trees must have the same label in each tree.

## 4. Infer the ASTRAL species tree

Combine one best gene tree per locus. Use `.treefile` outputs rather than bootstrap replicate trees.

The following block stops if any required tree is missing or empty:

```bash
cd "$PROJ_DIR"

(
    set -e

    for gene in brca1 foxp2 egfr rag1 rbp3 ghr vwf apob; do
        test -s "$gene/$gene.treefile"
    done

    for gene in brca1 foxp2 egfr rag1 rbp3 ghr vwf apob; do
        cat "$gene/$gene.treefile"
        printf '\n'
    done > species_tree/gene_trees.txt

    srun java -jar Astral/astral.5.7.8.jar \
        -i species_tree/gene_trees.txt \
        -o species_tree/species_tree.nwk \
        2> species_tree/astral.log
)

cat species_tree/species_tree.nwk
```

ASTRAL treats input gene trees as unrooted and estimates an unrooted species tree. Its default support labels are **local posterior probabilities**, ranging from 0 to 1. IQ-TREE’s gene-tree labels are **ultrafast bootstrap percentages**, ranging from 0 to 100. These are different support measures. See the [ASTRAL documentation](https://github.com/smirarab/ASTRAL/blob/master/README.md).

## 5. View and root trees in Phylo.io

1. Print the species tree:

   ```bash
   cat "$PROJ_DIR/species_tree/species_tree.nwk"
   ```

2. Copy the entire Newick string, including its final semicolon.
3. Open [Phylo.io](https://www.phylo.io/viewer/). Choose **Add Tree**, paste the text, and select **Newick**. If offered, enable **Use internal label for branches (e.g. bootstrap values)** before loading.
4. Find `Canis_lupus_familiaris`. Select its terminal branch and use the rooting/rerooting action to place the root on that branch. Dog should appear as sister to the primate clade.
5. Toggle **Branch Labels/Support** in the display settings to show or hide internal support values. Menu wording may vary by viewer version.
6. Load individual gene trees using the same copy-and-paste process:

   ```bash
   cat "$PROJ_DIR/brca1/brca1.treefile"
   ```

   Substitute the other gene names to view their trees. Check the dog outgroup and display branch support. Use the comparison view to examine discordance.

Rooting in the browser changes the displayed tree; it does not overwrite the cluster’s Newick file. Save or export the rooted view if needed.

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
