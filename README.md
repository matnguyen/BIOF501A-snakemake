# Standard Operating Procedure
## Background and Rationale

### Dependencies
The following dependencies are required, and defined in the environment YAML file for easy installation via Conda.

* Python (version 3.6.11)
* Snakemake (version 5.28.0)
* MAFFT (version 7.471)
* RAxML-NG (version 1.0.1)
* ETE 3 (version 3.1.2)


## Usage
All datasets are provided in the Github repository, along with the Snakemake workflow and a YAML file for creating a conda environment. 

### Installation
First, clone this repository by running: 

`git clone git@github.com:matnguyen/BIOF501A-snakemake.git`

Then create a Conda environment named `sop` using the included `env.yml` file:

`conda create -f env.yml`

This will install all the required dependencies. Activate the environment by running:

`conda activate sop`

### Running the pipeline
The snakemake pipeline can be run using:

`snakemake -s analysis.smk --cores 1`

## Input
The input are two FASTA files: one containing the multiple sequence alignment of 184,718 SARS-CoV-2 genomes from around the world, and one containing 196,674 SARS-CoV-2 spike protein sequences from around the world. These FASTA files will then be preprocessed in order to obtain paired genome and spike protein sequences of isolates from British Columbia, Canada. Both FASTA files are single-line formatted, whereby the first line starting with the `>` character contains the metadata of the isolate, and the second line containing the entire genome or protein sequence. The metadata field contains multiple fields separated by the `|` character: the protein name (if applicable), the virus name/collection location/isolate ID/collection year, the accesion ID, the collection date, and the collection location. There may also be other metadata fields such as the collection institution or the name of the submitter.   

## Output
