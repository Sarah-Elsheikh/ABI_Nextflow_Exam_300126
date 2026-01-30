params.accession = "M21012"
params.in = null 
params.out = "data_out"

// nextflow run main.nf --in "Nextflow/exam/hepatitis/*.fasta"

workflow {
// check that the user provides input Fasta files
    if (params.in == null) {
        println ("Error: provide input directory containing Fasta files")
        exit 1
    }

//    println("$params.accession, $params.in, $params.out")

    def ch_reference = fetch_reference(params.accession) 
        | view
    
    def ch_input = Channel.fromPath(params.in)
        // .mix (ch_reference)
        // .collect ()
    
    combine_fastas(ch_input)
        | view




    //def ch_input = channel.fromSRA(params.data, apiKey: "d5822ef54698cb072e0cf866736fd5f6ab08")
    //    | view

//     def ch_input = fetch_data(params.data) 
//     //    | view

//     // now this is how we connect


//     // we will use ch_input multiple times, so we define them here
//     // we output 3 different channels in the output
//     fastp(ch_input).sample_ID
//         | view

//     fastqc(ch_input)
//         | view
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

// process fetch_data {
//     conda "bioconda::sra-tools=3.2.1"

//     input:
//     val sra_num

//     output:
//     path "*.fastq.gz"

//     script:
//     """
//     prefetch "$sra_num"
//     """
// }
