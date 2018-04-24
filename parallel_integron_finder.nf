#!/usr/bin/env nextflow

replicons_file = file(params.replicons)
params.cpu = 1
params.gbk = false

println(params)

process split{

    input:
        file replicons from replicons_file

    output:
        file "*.fst" into chunk_files mode flatten

    """
    integron_split ${replicons}
    """
}


process integron_finder{
    input:
        file one_chunk from chunk_files

    output:
        file "Results_Integron_Finder_${one_chunk.baseName}" into all_chunk_results_dir

    script:
        """
        integron_finder --gbk --pdf --cpu ${params.cpu} --mute ${one_chunk}
        """
}


process merge{
    publishDir "Results_Integron_Finder_${replicons_file.simpleName}"

    input:
        file all_chunk_results from all_chunk_results_dir.toList()

    output:
        file "Results_Integron_Finder_${replicons_file.simpleName}"

    """
    integron_merge "Results_Integron_Finder_${replicons_file.baseName}" "${replicons_file.baseName}" ${all_chunk_results}
    """
}


workflow.onComplete {
    if ( workflow.success )
        println("\nDone!")
        println("Results are in --> ${workDir}/zzzzzz/xxxx/Results_Integron_Finder_${replicons_file.simpleName}")

}

workflow.onError {
    println "Oops .. something went wrong"
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}



