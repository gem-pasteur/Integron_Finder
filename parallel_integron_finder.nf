#!/usr/bin/env nextflow

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
params['promoter-attI'] = false
params['max-attc-size'] = false
params['min-attc-size'] = false
params.circ = false
params.linear = false
params['topology-file'] = false
params['keep-tmp'] = false
params['calin-threshold'] = false
params.out = false

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
promoter = params['promoter-attI'] ? '--promoter-attI' : ''
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
if (params.circ && params.linear){
    throw new Exception("The options '--linear' and '--circ' are mutually exclusive.")
}

if(params.replicons.tokenize(',').size() > 1 && params.out){
    throw new Exception("The options '--out' cannot be set when several replicon files are provides.")
}

replicons_file = Channel.fromPath(params.replicons.tokenize(',')).flatten()


/****************************************
 *           The workflow               *
 ****************************************/

process split{

    input:
        file(replicons) from replicons_file

    output:
        val("${replicons.baseName}") into input_id
        set val("${replicons.baseName}"), file("*.fst") into chunk_files mode flatten
        stdout nb_chunks
    script:
        """
        integron_split --mute ${replicons} | wc -w
        """
}

/**********************************************************
 need to emit as nb_chunks values as chunk_files
 otherwise only one integron_finder process is executed
 The first step is to create a channel containing [input_id, chunk_nb]
 The second step is to emit as many times there is chunk_files
************************************************************/
id_nb = input_id.merge(nb_chunks)
if_inputs = id_nb.cross(chunk_files)

process integron_finder{

    input:
        set val(input_id_nb), val(id_one_replicon) from if_inputs
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
        val promoter
        val max_attc_size
        val min_attc_size
        val keep_tmp
        val calin_threshold
    output:
        set val(id), file("Results_Integron_Finder_${one_replicon.baseName}") into all_chunk_results_dir

    script:
        (id, nb) = input_id_nb
        (_, one_replicon) = id_one_replicon

        if (params.circ){
            topo = '--circ'
        } else if (params.linear){
            topo = '--linear'
        } else if ( nb_chunks == 1) {
            topo = '--circ'
        } else {
            topo = '--linear'
        }

        """
        integron_finder ${local_max} ${func_annot} ${path_func_annot} ${dist_thr} ${union_integrases} ${attc_model} ${evalue_attc} ${keep_palindrome} ${no_proteins} ${promoter} ${max_attc_size} ${min_attc_size} ${calin_threshold} ${topo} ${topology_file} ${gbk} ${pdf} ${keep_tmp} --cpu ${task.cpus} --mute ${one_replicon}
        """
}


grouped_results = all_chunk_results_dir.groupTuple(by:0)


process merge{

    input:
        set val(input_id), file(all_chunk_results) from grouped_results

    output:
        set val(input_id), file ("${result_dir}/*") into final_res mode flatten
        
    script:
        res_dir_suffix = params.out ? params.out : input_id
        result_dir = "Results_Integron_Finder_${res_dir_suffix}"
        """
        integron_merge "${result_dir}" "${res_dir_suffix}" ${all_chunk_results}
        """
}


final_res.subscribe{
    input_id, result ->
        res_dir_suffix = params.out ? params.out : input_id
        result_dir = "Results_Integron_Finder_${res_dir_suffix}"
        result.copyTo("${result_dir}/" + result.name);
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
