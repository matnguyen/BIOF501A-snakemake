rule all:
    input:
        'bc_msa_tree.png',
        'bc_spikeprot_tree.png'

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
#   - text file of list of BC isolates with genomes
#   - text file of list of BC isolates with spike protein sequences
rule obtain_isolate_ids:
    input: 
        genome_fasta = 'msa_1119.fasta',
        spike_fasta = 'spikeprot1119.fasta'
    output:
        genome_txt = 'bc_genome_ids.txt',
        spike_txt = 'bc_spike_ids.txt'
    shell:
        'grep BCCDC {input.spike_fasta} | cut -f 2 -d "|" | '
        'cut -f 3 -d "/" | sort >> {output.spike_txt}; '
        'grep BC_ {input.genome_fasta} | cut -f 3 -d "/" | '
        'sort >> {output.genome_txt}' 

# Obtain list of BC isolates with both spike protein sequence and 
# genome sequence
# Output:
#   - text file of BC isolates with paired protein and genome sequences
rule obtain_paired_ids:
    input:
        genome_txt = 'bc_genome_ids.txt',
        spike_txt = 'bc_spike_ids.txt'
    output:
        ids = 'bc_paired_ids.txt'
    shell:
        'comm -12 {input.genome_txt} {input.spike_txt} >> {output.ids}'

# Preprocess the spike protein sequences fasta to obtain only sequences of 
# isolates from British Columbia
# Output: 
#   - fasta of only BC spike protein sequences
rule preprocess_spike_fasta:
    input:
        fasta = 'spikeprot1119.fasta',
        ids = 'bc_paired_ids.txt'
    output:
        fasta = 'bc_spikeprot.fasta',
    shell:
        'while read p; do grep -A1 -m 1 $p {input.fasta} >> {output.fasta}; '
        'done < {input.ids}'      

# Preprocessing multiple sequence alignment fasta to obtain only sequences of 
# isolates from British Columbia
# Output:
#   - fasta of only aligned BC genome sequences
rule preprocess_genome_msa:
    input: 
        fasta = 'msa_1119.fasta',
        ids = 'bc_paired_ids.txt'
    output:
        fasta = 'bc_msa.fasta'
    shell:
        'while read p; do grep -A1 -m 1 $p {input.fasta} >> {output.fasta}; '
        'done < {input.ids}'

# Run MAFFT to obtain a multiple sequence alignment of the protein spike sequences
# Output:
#   - fasta of multiple sequence alignment of spike protein sequences of BC isolates
rule run_mafft_spike:
    input:
        fasta = 'bc_spikeprot.fasta'
    output:
        fasta = 'bc_spikeprot_aligned.fasta'
    shell:
        'mafft {input.fasta} >> {output.fasta}'

# Trims the fasta headers so they are compatible with RAxML
# Output:
#   - fasta files with corrected headers
rule rename_fasta_headers:
    input:
        genome_fasta = 'bc_msa.fasta',
        spike_fasta = 'bc_spikeprot.fasta'
    output:
        genome_fasta = 'bc_msa_fixed.fasta',
        spike_fasta = 'bc_spikeprot_fixed.fasta'
    shell:
        'sed "s/.*\(BC.*\/2020\).*/>\\1/" {input.genome_fasta} >> {output.genome_fasta}; '
        'sed "s/.*\(BC.*\/2020\).*/>\\1/" {input.spike_fasta} >> {output.spike_fasta}'

# Runs RAxML for the spike protein alignment file to produce a 
# phylogenetic tree
# Output:
#   - Best scoring maximum likelihood (ML) tree
#   - Multifurcating version of the best-scoring ML tree
#   - ML trees for each starting tree
#   - Optimized model parameters for best-scoring ML tree
rule run_raxml_spike:
    input:
        fasta = 'bc_spikeprot_fixed.fasta'
    output:
        bestTree = 'bc_spikeprot.raxml.bestTree',
        bestTree_collapse = 'bc_spikeprot.raxml.bestTreeCollapsed',
        trees = 'bc_spikeprot.raxml.mlTrees',
        model = 'bc_spikeprot.raxml.bestModel'
    shell:
        'raxml-ng --msa {input.fasta} --msa-format FASTA --model LG+G '
        '--prefix bc_spikeprot --threads 4 --seed 667'

# Runs RAxML for the genome alignment file to produce a 
# phylogenetic tree
# Output:
#   - Best scoring maximum likelihood (ML) tree
#   - Multifurcating version of the best-scoring ML tree
#   - ML trees for each starting tree
#   - Optimized model parameters for best-scoring ML tree
rule run_raxml_genome:
    input:
        fasta = 'bc_msa_fixed.fasta'
    output:
        bestTree = 'bc_msa.raxml.bestTree',
        bestTree_collapse = 'bc_msa.raxml.bestTreeCollapsed',
        trees = 'bc_msa.raxml.mlTrees',
        model = 'bc_msa.raxml.bestModel'
    shell:
        'raxml-ng --msa {input.fasta} --msa-format FASTA --model GTR+G '
        '--prefix bc_msa --threads 4 --seed 667'

# Visualize the RAxML trees using ETE
# Output:
#   - PNG files for each phylogenetic tree
rule visualize_trees:
    input:
        genome_tree = 'bc_msa.raxml.bestTree',
        spike_tree = 'bc_spikeprot.raxml.bestTree'
    output:
        genome_png = 'bc_msa_tree.png',
        spike_png = 'bc_spikeprot_tree.png'
    shell:
        'ete3 view -i {output.genome_png} -t {input.genome_tree} -m c; '
        'ete3 view -i {output.spike_png} -t {input.spike_tree} -m c'