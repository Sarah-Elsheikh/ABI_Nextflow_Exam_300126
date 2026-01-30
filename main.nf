params.accession = "M21012"
params.in = null 
params.out = "data_out"

// nextflow run main.nf --in "Nextflow/exam/hepatitis/*.fasta"

workflow {
// check that the user provides input Fasta files
    if (params.in == null) {
        println ("Error: provide path to directory containing Fasta files")
        exit 1
    }

//    println("$params.accession, $params.in, $params.out")

    def ch_reference = fetch_reference(params.accession) 
    
    def ch_input = Channel.fromPath("${params.in}/*.fasta", checkIfExists: true)
    
    ch_input
        .mix(ch_reference)
        .collect()
        .set {ch_combined_for_merge}
    
    def ch_combine_fastas = combine_fastas(ch_combined_for_merge)

    def ch_alignment = align_mafft(ch_combine_fastas)

    trimal_cleanup(ch_alignment).fasta
        // |view
}

// fetch the reference from the database
process fetch_reference {
    conda "bioconda::entrez-direct=24.0"

    input:
    val accession

    output:
    path "${accession}.fasta"  

    script:
    """
    esearch -db nucleotide -query "$accession" \\
        | efetch -format fasta > "${accession}.fasta"  
    """
}

// Merges all individual FASTA files into a single multi-FASTA [cite: 22]
process combine_fastas {
    conda 'conda-forge::sed'

    input:
    path fasta_files

    output:
    path "combined_sequences.fasta"

    script:
    """
    cat ${fasta_files} > combined_sequences.fasta
    """
}

process align_mafft {
    conda 'bioconda::mafft=7.525'

    input:
    path combined_fasta

    output:
    path "alignment.fasta"

    script:
    """
    mafft ${combined_fasta} > alignment.fasta
    """
}

process trimal_cleanup {
    publishDir "${params.out}", mode: 'copy'
    conda 'bioconda::trimal=1.5.0'

    input:
    path alignment

    output:
    path "alignment_trimmed.fasta", emit: fasta
    path "alignment_trimmed.html", emit: html

    script:
    """
    trimal -in ${alignment} \
           -out alignment_trimmed.fasta \
           -htmlout alignment_trimmed.html \
           -automated1
    """
}

