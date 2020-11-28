configfile: 'config.yaml'

COUNTRY = config['filter']['country']
PREFIX = config['filter']['prefix']

rule all:
    input:
        expand('{prefix}_msa_tree.png', prefix=PREFIX),
        expand('{prefix}_spikeprot_tree.png', prefix=PREFIX),
        expand('{prefix}_trees_comparison.txt', prefix=PREFIX)

# Decompress the dataset to obtain a multiple sequence alignment of genomes (fasta)
# and a fasta of unaligned spike protein sequences
# Output:
#   - decompressed multiple sequence alignment fasta
#   - decompressed unaligned spike protein fasta
rule decompress_data:
    input:
        genome = 'msa_1119.tar.xz',
        spike = 'spikeprot1119.fasta.xz'
    output:
        spike_fasta = 'spikeprot1119.fasta',
        genome_fasta = 'msa_1119.fasta'
    shell:
        'tar -xf {input.genome}; unxz {input.spike}'

# Obtain list of all British Columbia isolates with genomes, and all isolates
# with apike protein sequences
# Output:
#   - text file of list of Vietnam isolates with genomes
#   - text file of list of Vietnam isolates with spike protein sequences
rule obtain_isolate_ids:
    input: 
        genome_fasta = 'msa_1119.fasta',
        spike_fasta = 'spikeprot1119.fasta'
    output:
        genome_txt = expand('{prefix}_genome_ids.txt', prefix=PREFIX),
        spike_txt = expand('{prefix}_spike_ids.txt', prefix=PREFIX)
    shell:
        'grep {COUNTRY} {input.spike_fasta} | cut -f 2 -d "|" | sort >> {output.spike_txt}; '
        'grep {COUNTRY} {input.genome_fasta} | cut -f 1 -d "|" | sed "s/>//"  | '
        'sort >> {output.genome_txt}' 

# Obtain list of Vietnam isolates with both spike protein sequence and 
# genome sequence
# Output:
#   - text file of Vietnam isolates with paired protein and genome sequences
rule obtain_paired_ids:
    input:
        genome_txt = expand('{prefix}_genome_ids.txt', prefix=PREFIX),
        spike_txt = expand('{prefix}_spike_ids.txt', prefix=PREFIX)
    output:
        ids = expand('{prefix}_paired_ids.txt', prefix=PREFIX)
    shell:
        'comm -12 {input.genome_txt} {input.spike_txt} >> {output.ids}'

# Preprocess the spike protein sequences fasta to obtain only sequences of 
# isolates from British Columbia
# Output: 
#   - fasta of only Vietnam spike protein sequences
rule preprocess_spike_fasta:
    input:
        fasta = 'spikeprot1119.fasta',
        ids = expand('{prefix}_paired_ids.txt', prefix=PREFIX)
    output:
        fasta = expand('{prefix}_spikeprot.fasta', prefix=PREFIX)
    shell:
        'while read p; do grep -A1 -m 1 $p {input.fasta} >> {output.fasta}; '
        'done < {input.ids}'      

# Preprocessing multiple sequence alignment fasta to obtain only sequences of 
# isolates from British Columbia
# Output:
#   - fasta of only aligned Vietnam genome sequences
rule preprocess_genome_msa:
    input: 
        fasta = 'msa_1119.fasta',
        ids = expand('{prefix}_paired_ids.txt', prefix=PREFIX)
    output:
        fasta = expand('{prefix}_msa.fasta', prefix=PREFIX)
    shell:
        'while read p; do grep -A1 -m 1 $p {input.fasta} >> {output.fasta}; '
        'done < {input.ids}'

# Run MAFFT to obtain a multiple sequence alignment of the protein spike sequences
# Output:
#   - fasta of multiple sequence alignment of spike protein sequences of Vietnam isolates
rule run_mafft_spike:
    input:
        fasta = expand('{prefix}_spikeprot.fasta', prefix=PREFIX)
    output:
        fasta = expand('{prefix}_spikeprot_aligned.fasta', prefix=PREFIX)
    shell:
        'mafft {input.fasta} >> {output.fasta}'

# Trims the fasta headers so they are compatible with RAxML
# Output:
#   - fasta files with corrected headers
rule rename_fasta_headers:
    input:
        genome_fasta = expand('{prefix}_msa.fasta', prefix=PREFIX)
        spike_fasta = expand('{prefix}_spikeprot.fasta', prefix=PREFIX)
    output:
        genome_fasta = expand('{prefix}_msa_fixed.fasta', prefix=PREFIX),
        spike_fasta = expand('{prefix}_spikeprot_fixed.fasta', prefix=PREFIX)
    shell:
        'sed "s/.*V{COUNTRY}\/\(.*\)\/2020.*/>\\1/" {input.genome_fasta} >> {output.genome_fasta}; '
        'sed "s/.*{COUNTRY}\/\(.*\)\/2020.*/>\\1/" {input.spike_fasta} >> {output.spike_fasta}'

# Runs RAxML for the spike protein alignment file to produce a 
# phylogenetic tree
# Output:
#   - Best scoring maximum likelihood (ML) tree
#   - Multifurcating version of the best-scoring ML tree
#   - ML trees for each starting tree
#   - Optimized model parameters for best-scoring ML tree
rule run_raxml_spike:
    input:
        fasta = expand('{prefix}_spikeprot_fixed.fasta', prefix=PREFIX)
    output:
        bestTree = expand('{prefix}_spikeprot.raxml.bestTree', prefix=PREFIX),
        bestTree_collapse = expand('{prefix}_spikeprot.raxml.bestTreeCollapsed', prefix=PREFIX),
        trees = expand('{prefix}_spikeprot.raxml.mlTrees', prefix=PREFIX),
        model = expand('{prefix}_spikeprot.raxml.bestModel', prefix=PREFIX)
    shell:
        'raxml-ng --msa {input.fasta} --msa-format FASTA --model LG+G '
        '--prefix {PREFIX}_spikeprot --threads 4 --seed 667'

# Runs RAxML for the genome alignment file to produce a 
# phylogenetic tree
# Output:
#   - Best scoring maximum likelihood (ML) tree
#   - Multifurcating version of the best-scoring ML tree
#   - ML trees for each starting tree
#   - Optimized model parameters for best-scoring ML tree
rule run_raxml_genome:
    input:
        fasta = 'vn_msa_fixed.fasta'
    output:
        bestTree = expand('{prefix}_msa.raxml.bestTree', prefix=PREFIX),
        bestTree_collapse = expand('{prefix}_msa.raxml.bestTreeCollapsed', prefix=PREFIX),
        trees = expand('{prefix}_msa.raxml.mlTrees', prefix=PREFIX),
        model = expand('{prefix}_msa.raxml.bestModel', prefix=PREFIX)
    shell:
        'raxml-ng --msa {input.fasta} --msa-format FASTA --model GTR+G '
        '--prefix {PREFIX}_msa --threads 4 --seed 667'

# Visualize the RAxML trees using ETE
# Output:
#   - PNG files for each phylogenetic tree
rule visualize_trees:
    input:
        genome_tree = expand('{prefix}_msa.raxml.bestTree', prefix=PREFIX),
        spike_tree = expand('{prefix}_spikeprot.raxml.bestTree', prefix=PREFIX)
    output:
        genome_png = expand('{prefix}_msa_tree.png', prefix=PREFIX),
        spike_png = expand('{prefix}_spikeprot_tree.png', prefix=PREFIX)
    shell:
        'ete3 view -i {output.genome_png} -t {input.genome_tree} -m c; '
        'ete3 view -i {output.spike_png} -t {input.spike_tree} -m c'

# Compare the two phylogenetic trees quantitatively
# Output:
#   - text file with metrics for comparison
rule compare_trees:
    input:
        genome_tree = expand('{prefix}_msa.raxml.bestTree', prefix=PREFIX),
        spike_tree = expand('{prefix}_spikeprot.raxml.bestTree', prefix=PREFIX)
    output:
        txt = expand('{prefix}_trees_comparison.txt', prefix=PREFIX)
    shell:
        'ete3 compare -t {input.spike_tree} -r {input.genome_tree} --unrooted | '
        'cut -f 3-10 -d "|" >> {output.txt}'