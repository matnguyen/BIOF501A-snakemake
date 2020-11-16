# Decompress the dataset to obtain a multiple sequence alignment of genomes (fasta)
# and a fasta of unaligned spike protein sequences
# Output:
#   - decompressed multiple sequence alignment fasta
#   - decompressed unaligned spike protein fasta
rule decompress_data:
    input:
        genome = 'msa_1114.tar.xz',
        spike = 'spikeprot1116.fasta.xz'
    output:
        spike_fasta = 'spikeprot1114.fasta',
        genome_fasta = 'msa_1114.fasta'
    shell:
        'tar -zvf {input.genome}; unxz {input.spike}'

# Obtain list of all British Columbia isolates with genomes, and all isolates
# with apike protein sequences
# Output:
#   - text file of list of BC isolates with genomes
#   - text file of list of BC isolates with spike protein sequences
rule obtain_isolate_ids:
    input: 
        genome_fasta = 'msa_1114.fasta'
        spike_fasta = 'spikeprot1114.fasta'
    output:
        genome_txt = 'bc_genome_ids.txt'
        spike_txt = 'bc_spike_ids.txt'
    shell:
        'grep BCCDC {input.spike_fasta} | cut -f 2 -d "|" | '
        'cut -f 3 -d "/" | sort >> {output.spike_txt}; '
        'grep BC_ {input.genome_fasta} | cut -f 3 -d "/" | '
        'sort >> {output.spike_txt}' 

# Obtain list of BC isolates with both spike protein sequence and 
# genome sequence
# Output:
#   - text file of BC isolates with paired protein and genome sequences
rule obtain_paired_ids:
    input:
        genome_txt = 'bc_genome_ids.txt'
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
        fasta = 'spikeprot1114.fasta'
        ids = 'bc_paired_ids.txt'
    output:
        fasta = 'bc_spikeprot.fasta',
        id_list = 'bc_ids.txt'
    shell:
        'while read p; do grep -A1 -m 1 $p {input.fasta} >> {output.fasta}; '
        'done < {input.ids}'      

# Preprocessing multiple sequence alignment fasta to obtain only sequences of 
# isolates from British Columbia
# Output:
#   - fasta of only aligned BC genome sequences
rule preprocess_genome_msa:
    input: 
        fasta = 'msa_1114.fasta'
        ids = 'bc_ids.txt'
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
        mafft {input.fasta} >> {output.fasta}