#!/usr/bin/env nextflow

replicons_file = Channel.fromPath(params.replicons)

/*************************
 * Default options
 *************************/

params.gbk = false
params.pdf = false
params['local-max']= false
params['func-annot'] = false
params['distance-threshold'] = false
params['union-integrases'] = false
params['path-func-annot'] = false
params['attc-model'] = false
params['evalue-attc'] = false
params['keep-palindrome'] = false
params['no-proteins'] = false
params['max-attc-size'] = false
params['min-attc-size'] = false
params.circ = false
params.linear = false
params['topology-file'] = false
params['keep-tmp'] = false
params['calin-threshold'] = false


gbk = params.gbk ? '--gbk' : ''
pdf = params.pdf ? '--pdf' : ''
local_max = params['local-max'] ? '--local-max' : ''
func_annot = params['func-annot'] ? '--func-annot' : ''
path_func_annot = params['path-func-annot'] ? "--path-func-annot ${params['path-func-annot']}" : ''
dist_thr = params['distance-threshold'] ? "--distance-thresh ${params['distance-threshold']}" : ''
union_integrases = params['union-integrases'] ? '--union-integrase' : ''
attc_model = params['attc-model'] ? "--attc-model ${params['attc-model']}" : ''
evalue_attc = params['evalue-attc'] ? "--evalue-attc ${params['evalue-attc']}" : ''
keep_palindrome = params['keep-palindrome'] ? '--keep-palindrome' : ''
no_proteins = params['no-proteins'] ? '--no-proteins' : ''
max_attc_size = params['max-attc-size'] ? "--max-attc-size ${params['max-attc-size']}" : ''
min_attc_size = params['min-attc-size'] ? "--min-attc-size ${params['min-attc-size']}" : ''
circ = params.circ ? '--circ' : ''
linear = params.linear ? '--linear' : ''
topology_file = params['topology-file'] ? "--topology-file ${params['topology-file']}" : ''
keep_tmp = params['keep-tmp'] ? '--keep-tmp' : ''
calin_threshold = params['calin-threshold'] ? "--calin-threshold ${params['calin-threshold']}" : ''

if (! params.replicons){
    throw new Exception("The option '--replicons' is mandatory.")
}

params.out = false
replicon_file = new File(params.replicons)
res_dir_suffix = params.out ? params.out : replicon_file.name.split("\\.", 2)[0]
result_dir = "Results_Integron_Finder_${res_dir_suffix}"


/****************************************
 *           The workflow               *
 ****************************************/

process split{

    input:
        file replicons from replicons_file

    output:
        file "*.fst" into chunk_files mode flatten

    script:
        """
        integron_split --mute ${replicons}
        """
}


process integron_finder{
    input:
        file one_chunk from chunk_files
        val gbk
        val pdf
        val local_max
        val func_annot
        val path_func_annot
        val circ
        val linear
        val topology_file
        val dist_thr
        val union_integrases
        val attc_model
        val evalue_attc
        val keep_palindrome
        val no_proteins
        val max_attc_size
        val min_attc_size
        val keep_tmp
        val calin_threshold
    output:
        file "Results_Integron_Finder_${one_chunk.baseName}" into all_chunk_results_dir

    script:
        """
        integron_finder ${local_max} ${func_annot} ${path_func_annot} ${dist_thr} ${union_integrases} ${attc_model} ${evalue_attc} ${keep_palindrome} ${no_proteins} ${max_attc_size} ${min_attc_size} ${calin_threshold} ${circ} ${linear} ${topology_file} ${gbk} ${pdf} ${keep_tmp} --cpu ${task.cpus} --mute ${one_chunk}
        """
}


process merge{
    publishDir "${result_dir}", mode:'copy'

    input:
        file all_chunk_results from all_chunk_results_dir.toList()

    output:
        file "${result_dir}"
        
    script:
        """
        integron_merge "${result_dir}" "${res_dir_suffix}" ${all_chunk_results}
        """
}


workflow.onComplete {
    if ( workflow.success )
        println("\nDone!")
        println("Results are in --> ${result_dir}")

}

workflow.onError {
    println "Oops .. something went wrong"
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}



