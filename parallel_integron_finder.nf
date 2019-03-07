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
params.gembase = false
params['gembase-path'] = false
params.debug = false

/****************************************
*  parameter not really used in the wf  *
*  used to improve UI                   *
*****************************************/
params.help = false
params.profile = false

/****************************************
*             real parameters           *
*****************************************/
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
gembase_path = params['gembase-path'] ? "--gembase-path ${params['gembase-path']}" : ''
debug = params.debug ? '-vv' : ''


if (params.help){
   msg = '''
parallel_integron_finder available options:

 --gbk
 --pdf
 --local-max
 --func-annot
 --distance-threshold
 --union-integrases
 --path-func-annot
 --attc-model
 --evalue-attc
 --keep-palindrome
 --no-proteins
 --promoter-attI
 --max-attc-size
 --min-attc-size
 --circ
 --linear
 --topology-file
 --keep-tmp
 --calin-threshold
 --gembase-path
 --replicons

Please refer to the integron_finder documentation for the meaning of each options.
'''
    println(msg)
    System.exit(0)
}
if (params.profile){
    throw new Exception("The integron_finder option '--profile' does not exists. May be you want to use the nextflow option '-profile'.")
}
if (! params.replicons){
    throw new Exception("The option '--replicons' is mandatory.")
}
if (params.circ && params.linear){
    throw new Exception("The options '--linear' and '--circ' are mutually exclusive.")
}
if (params.gembase){
    throw new Exception("The options '--gembase' is not available for parallel_integron_finder. Use '--gembase-path' instead")
}
if (params['gembase-path']){
    // need to compute the abspath to gembase
    // the integron_finder step will be compute in the work subdirectory
    file = new File(params['gembase-path'])
    gembase_path = file.getCanonicalPath();
    gembase_path = "--gembase-path ${gembase_path}"
} else {
    gembase_path = ''
}
if (params.replicons.contains(',')){
    paths = params.replicons.tokenize(',')
    replicons_file = Channel.fromPath(paths)
} else {
    replicons_file = Channel.fromPath(params.replicons)
}

/****************************************
 *           The workflow               *
 ****************************************/

process split{

    input:
        file(replicons) from replicons_file

    output:
        set val("${replicons.baseName}"), stdout, file("*.fst") into chunk_files mode flatten
    script:
        """
        integron_split --mute ${replicons} | wc -w
        """
}


process integron_finder{

    input:
        set val(id_input), val(nb_chunks), file(one_replicon) from chunk_files
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
        val gembase_path
        val debug
    output:
        set val(id_input), file("Results_Integron_Finder_${one_replicon.baseName}") into all_chunk_results_dir

    script:
        /*******************************************************************************************************
        For sequential IF the default topology behavior is
        if there is only one replicon in the replicon file then the topology is circular
        if there are several replicons in the replicon file then the topology is linear
        but for parallelisation we split replicon files in as file as replicon so the default topology is always circular.
        To restore the same behavior we need to analyse the number of replicon per file and force the topology
        ********************************************************************************************************/
        nb_chunks = nb_chunks as int

        if (params.circ){
            topo = '--circ'
        } else if (params.linear){
            topo = '--linear'
        } else if (nb_chunks == 1) {
            topo = '--circ'
        } else {
            topo = '--linear'
        }
        
        """
        integron_finder ${local_max} ${func_annot} ${path_func_annot} ${dist_thr} ${union_integrases} ${attc_model} ${evalue_attc} ${keep_palindrome} ${no_proteins} ${promoter} ${max_attc_size} ${min_attc_size} ${calin_threshold} ${topo} ${topology_file} ${gbk} ${pdf} ${keep_tmp} --cpu ${task.cpus} ${gembase_path} --mute ${debug} ${one_replicon}
        """
}


grouped_results = all_chunk_results_dir.groupTuple(by:0)


process merge{

    input:
        set val(input_id), file(all_chunk_results) from grouped_results

    output:
        set val(input_id), file ("${result_dir}/*") into final_res mode flatten
        
    script:
        result_dir = "Results_Integron_Finder_${input_id}"
        """
        integron_merge "${result_dir}" "${input_id}" ${all_chunk_results}
        """
}


final_res.subscribe{
    input_id, result ->
        result_dir = "Results_Integron_Finder_${input_id}"
        result.copyTo("${result_dir}/" + result.name);
}

workflow.onComplete {
    if ( workflow.success )
        println("\nDone!")
}

workflow.onError {
    println "Oops .. something went wrong"
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
