#!/usr/bin/env nextflow

replicons_file = Channel.fromPath(params.replicons)

/*************************
 * Default options
 *************************/

params.gbk = false
gbk = ''
params.pdf = false
pdf = ''
params.localMax= false
local_max = ''
params.funcAnnot = false
func_annot = ''
params.distanceThresh = false
dist_thr = ''
params.unionIntegrases = false
union_integrases = ''
params.pathFuncAnnot = false
path_func_annot = ''
params.attcModel = false
attc_model = ''
params.evalueAttc = false
evalue_attc = ''
params.keepPalindrome = false
keep_palindrome = ''
params.noProteins = false
no_proteins = ''
params.maxAttcSize = false
max_attc_size = ''
params.minAttcSize = false
min_attc_size = ''
params.circ = false
circ = ''
params.linear = false
linear = ''
params.topologyFile = false
topology_file = ''
params.keepTmp = false
keep_tmp = ''
params.calin_threshold = false
calin_threshold = ''
promoter_attI = false

if (params.gbk){
    gbk = '--gbk'
}

if (params.pdf){
    pdf = '--pdf'
}

if (params.localMax){
    local_max = '--local-max'
}

if (params.funcAnnot){
    func_annot = '--func-annot'
}

if (params.pathFuncAnnot){
    path_func_annot = "--path-func-annot ${params.pathFuncAnnot}"
}

if (params.distanceThresh){
    dist_thr = "--distance-thresh ${params.distanceThreshold}"
}

if (params.unionIntegrases){
    union_integrases = '--union-integrase'
}

if (params.attcModel){
    attc_model = "--attc-model ${params.attcModel}"
}

if (params.evalueAttc){
    evalue_attc = "--evalue-attc ${params.evalueAttc}"
}

if (params.keepPalindrome){
    keep_palindrome = '--keep-palindrome'
}

if (params.noProteins){
    no_proteins = '--no-proteins'
}

if (params.maxAttcSize){
    max_attc_size = "--max-attc-size ${params.maxAttcSize}"
}

if (params.minAttcSize){
    min_attc_size = "--min-attc-size ${params.minAttcSize}"
}

if (params.circ){
    circ = '--circ'
}

if (params.linear){
    linear = '--linear'
}

if (params.topologyFile){
    topology_file = "--topology-file ${params.topologyFile}"
}

if (params.keepTmp){
    keep_tmp = '--keepTmp'
}
if (params.calin_threshold){
    calin_threshold = "--calin-threshold ${params.calin_threshold}"
}
if (params.promoter_attI){
    promoter_attI = "--promoter-attI ${params.promoter_attI}"
}

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
        integron_finder ${local_max} ${func_annot} ${path_func_annot} ${dist_thr} ${union_integrases} ${attc_model} ${evalue_attc} ${keep_palindrome} ${no_proteins} ${max_attc_size} ${min_attc_size} ${calin_threshold} ${params.promoter_attI} ${circ} ${linear} ${topology_file} ${gbk} ${pdf} ${keep_tmp} --cpu ${task.cpu} --mute ${one_chunk}
        """
}


process merge{
    publishDir "Results_Integron_Finder_${params.out}"

    input:
        file all_chunk_results from all_chunk_results_dir.toList()

    output:
        file "Results_Integron_Finder_${params.out}"
        
    script:
        """
        integron_merge "Results_Integron_Finder_${params.out}" "${params.out}" ${all_chunk_results}
        """
}


workflow.onComplete {
    if ( workflow.success )
        println("\nDone!")
        println("Results are in --> Results_Integron_Finder_${params.out}")

}

workflow.onError {
    println "Oops .. something went wrong"
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}



