rule all:
    input:
        'vn_msa_tree.pdf',
        'vn_spikeprot_tree.pdf'

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
        genome_txt = 'vn_genome_ids.txt',
        spike_txt = 'vn_spike_ids.txt'
    shell:
        'grep Vietnam {input.spike_fasta} | cut -f 2 -d "|" | sort >> {output.spike_txt}; '
        'grep Vietnam {input.genome_fasta} | cut -f 1 -d "|" | sed "s/>//"  | '
        'sort >> {output.genome_txt}' 

# Obtain list of Vietnam isolates with both spike protein sequence and 
# genome sequence
# Output:
#   - text file of Vietnam isolates with paired protein and genome sequences
rule obtain_paired_ids:
    input:
        genome_txt = 'vn_genome_ids.txt',
        spike_txt = 'vn_spike_ids.txt'
    output:
        ids = 'vn_paired_ids.txt'
    shell:
        'comm -12 {input.genome_txt} {input.spike_txt} >> {output.ids}'

# Preprocess the spike protein sequences fasta to obtain only sequences of 
# isolates from British Columbia
# Output: 
#   - fasta of only Vietnam spike protein sequences
rule preprocess_spike_fasta:
    input:
        fasta = 'spikeprot1119.fasta',
        ids = 'vn_paired_ids.txt'
    output:
        fasta = 'vn_spikeprot.fasta',
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
        ids = 'vn_paired_ids.txt'
    output:
        fasta = 'vn_msa.fasta'
    shell:
        'while read p; do grep -A1 -m 1 $p {input.fasta} >> {output.fasta}; '
        'done < {input.ids}'

# Run MAFFT to obtain a multiple sequence alignment of the protein spike sequences
# Output:
#   - fasta of multiple sequence alignment of spike protein sequences of Vietnam isolates
rule run_mafft_spike:
    input:
        fasta = 'vn_spikeprot.fasta'
    output:
        fasta = 'vn_spikeprot_aligned.fasta'
    shell:
        'mafft {input.fasta} >> {output.fasta}'

# Trims the fasta headers so they are compatible with RAxML
# Output:
#   - fasta files with corrected headers
rule rename_fasta_headers:
    input:
        genome_fasta = 'vn_msa.fasta',
        spike_fasta = 'vn_spikeprot_aligned.fasta'
    output:
        genome_fasta = 'vn_msa_fixed.fasta',
        spike_fasta = 'vn_spikeprot_fixed.fasta'
    shell:
        'sed "s/.*Vietnam\/\(.*\)\/2020.*/>\\1/" {input.genome_fasta} >> {output.genome_fasta}; '
        'sed "s/.*Vietnam\/\(.*\)\/2020.*/>\\1/" {input.spike_fasta} >> {output.spike_fasta}'

# Runs RAxML for the spike protein alignment file to produce a 
# phylogenetic tree
# Output:
#   - Best scoring maximum likelihood (ML) tree
#   - Multifurcating version of the best-scoring ML tree
#   - ML trees for each starting tree
#   - Optimized model parameters for best-scoring ML tree
rule run_raxml_spike:
    input:
        fasta = 'vn_spikeprot_fixed.fasta'
    output:
        bestTree = 'vn_spikeprot.raxml.bestTree',
        bestTree_collapse = 'vn_spikeprot.raxml.bestTreeCollapsed',
        trees = 'vn_spikeprot.raxml.mlTrees',
        model = 'vn_spikeprot.raxml.bestModel'
    shell:
        'raxml-ng --msa {input.fasta} --msa-format FASTA --model LG+G '
        '--prefix vn_spikeprot --threads 4 --seed 667'

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
        bestTree = 'vn_msa.raxml.bestTree',
        bestTree_collapse = 'vn_msa.raxml.bestTreeCollapsed',
        trees = 'vn_msa.raxml.mlTrees',
        model = 'vn_msa.raxml.bestModel'
    shell:
        'raxml-ng --msa {input.fasta} --msa-format FASTA --model GTR+G '
        '--prefix vn_msa --threads 4 --seed 667'

# Visualize the RAxML trees using ETE
# Output:
#   - PNG files for each phylogenetic tree
rule visualize_trees:
    input:
        genome_tree = 'vn_msa.raxml.bestTree',
        spike_tree = 'vn_spikeprot.raxml.bestTree'
    output:
        genome_png = 'vn_msa_tree.pdf',
        spike_png = 'vn_spikeprot_tree.pdf'
    shell:
        'ete3 view -i {output.genome_png} -t {input.genome_tree} -m c; '
        'ete3 view -i {output.spike_png} -t {input.spike_tree} -m c'