# SARS-CoV-2 Genome and Spike Protein Phylogenetic Comparison Pipeline - Standard Operating Procedure
## Matthew Nguyen
## Background and Rationale

### Dependencies
The following dependencies are required, and defined in the environment YAML file for easy installation via Conda.

* Python (version 3.6.11)
* Snakemake (version 5.28.0)
* MAFFT (version 7.471)
* RAxML-NG (version 1.0.1)
* ETE 3 (version 3.1.2)

### Directed Acyclic Graph for Workflow
![](images/dag.svg) 


## Usage
All datasets are provided in the Github repository, along with the Snakemake workflow and a YAML file for creating a conda environment. 

### Installation
First, clone this repository by running: 

`git clone git@github.com:matnguyen/BIOF501A-snakemake.git`

Then create a Conda environment named `sop` using the included `env.yml` file:

`conda create -f env.yml`

This will install all the required dependencies. Activate the environment by running:

`conda activate sop`

### Running the workflow
The snakemake pipeline can be run using:

`snakemake --cores 1`

### Steps of the workflow
#### decompress_data
This step decompresses the `msa_1119.tar.xz` and `spikeprot1119.fasta.xz` using the commands `tar -xf msa_1119.tar.xz` and `unxz spikeprot1119.fasta.xz` respectively to produce two uncompressed FASTA files: `msa_1119.fasta` and `spikeprot1119.fasta`

#### obtain_isolate_ids
This step produces two text files containing a list of IDs of Vietnamese SARS-CoV-2 isolates, one from the genome FASTA and one from the spike protein FASTA respectively. This is done using a simple `bash` pipeline comprising of `grep`, `cut`, `sort` and `sed`.

#### obtain_paired_ids
This step produces a single text file containing a list of IDs present in both the genome and spike protein FASTA files, to filter out all IDs that only occur in once, and don't have a paired spike protein and genome sequence. This is done by using the `comm` command in `bash`.

#### preprocess_spike_fasta
This step processes the spike protein FASTA to only include samples present in the paired IDs list. This is a simple `while` loop with `grep`. 

#### preprocess_genome_fasta
The same process outlined for `preprocess_spike_fasta` is done for the genome fasta, producing a processed FASTA file with only samples contained in the paired IDs list.

#### run_mafft_spike
Since the spike protein FASTA is unaligned, this step performs a multiple sequence alignment of the sequences in the FASTA using `mafft`.

#### rename_fasta_headers
This step renames the FASTA headers for both the spike protein and genome FASTA files in order to be compatible with RAxML. This is done using a simple `sed` command.

#### run_raxml_spike:
This step runs `RAxML` to produce a phylogenetic tree for the spike protein sequences from the FASTA file.

#### run_raxml_genome:
Similarly, this step runs `RAxML` to produce a phylogenetic tree for the genome sequences.

#### visualize_tree:
This step produces `PNG` files for visualizing the phylogenetic trees of the spike protein and genome sequences. It takes in the best derived tree from multiple bootstrap runs of `RAxML`. This uses `ete3 view` for visualization.

#### compare_trees
This step uses `ete3 compare` to quantitatively compare the two phylogenetic trees obtained. The reference tree is the genome tree, and the target tree is the spike protein tree. It outputs a text file containing different metrics, which will be explained in the Output section below.

#### all
This step just ensures all the correct output files are generated.


## Input
`config.yaml` is a configuration file passed in to the Snakemake workflow. Here we define a filter country, a prefix for the filenames, an output directory and the number of threads. For this example, our filter country is `Vietnam`, our prefix is `vn`, our output directory is `./outputs` and we define `4` threads. These settings can be changed to fit the user's goal and computational resources. 

The main input for the workflow are two FASTA files: one containing the multiple sequence alignment of 184,718 SARS-CoV-2 genomes from around the world, and one containing 196,674 SARS-CoV-2 spike protein sequences from around the world. These FASTA files will then be preprocessed in order to obtain paired genome and spike protein sequences of isolates from a country or location determined by the configuration file. Both FASTA files are single-line formatted, whereby the first line starting with the `>` character contains the metadata of the isolate, and the second line containing the entire genome or protein sequence. The metadata field contains multiple fields separated by the `|` character: the protein name (if applicable), the virus name/collection location/isolate ID/collection year, the accesion ID, the collection date, and the collection location. There may also be other metadata fields such as the collection institution or the name of the submitter.   

## Output
### Overview
By default, this workflow will create a directory `outputs` to hold the output files (asides from the initial decompressed data). This directory can be change in `config.yaml`. The three main final outputs are a phylogenetic tree for the genome sequences, a phylogenetic tree for the spike protein sequences, and a text file quantitatively comparing the two phylogenetic trees. `outputs/intermediate_files` contains other files that are generated from the workflow, but are not important for the final analysis. These include aligned and processed FASTA files, and phylogenetic trees and models generated by RAxML.

Below, we see the two phylogenetic trees generated by the workflow. 

### Phylogenetic tree obtained from genome sequences
![](/images/vn_msa_tree.png)

### Phylogenetic tree obtained from spike protein sequences
![](/images/vn_spikeprot_tree.png)

Qualitatively comparing the two phylogenetic trees, we can see that the two trees are quite different. Most sequences in the spike protein sequences are very genetically similar, with only two small clades that are significantly genetically distinct from the rest. Whereas when we look at the phylogenetic tree generated from genome sequences, we see a lot more distinct clades and a lot of genetic variation. This is an ideal outcome, to see that even genetically different SARS-CoV-2 isolates have genetically similar spike protein sequences, thus meaning the sequence is highly conserved for Vietnamese samples. This means that a vaccine targetting the spike protein may be very effective against SARS-CoV-2 strains found in Vietnam. 

### Quantitative tree comparison
```
 E.size  | nRF     | RF      | maxRF   | src-br+ | ref-br+ | subtre+ | treekoD
 ======+ | ======+ | ======+ | ======+ | ======+ | ======+ | ======+ | ======+
 86      | 0.96    | 160.00  | 166.00  | 0.52    | 0.52    | 1       | NA
 ```
 
Here we see the output of the text file containing the quantitative measures for phylogenetic tree difference between the tree generated by genome sequences and the tree generated by spike protein sequences. We arbitrarily define the target tree as the spike protein tree, and the reference tree as the genome sequence tree. The columns of the text file are:

* `E.size`: the number of isolates in the phylogenetic tree
* `nRF`: normalized Robinson-Foulds distance (RF/maxRF)
* `RF`: Robinson-Foulds symmetric distance
* `maxRF`: maximum Robinson-Foulds value for this comparison
* `src-br+`: frequency of edges in target tree found in the reference
* `ref_br+`: frequency of edges in the reference tree found in target
* `subtre+`: number of subtrees used for the comparison
* `treekoD`: average distance among all possible subtrees in the original target trees to the reference tree
