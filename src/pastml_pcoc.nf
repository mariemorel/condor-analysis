//ancestral reconstruction and pcoc

path_file = "$baseDir"

//rhodopsin 
params.tree = path_file + "data/rhodopsin/tree_rerooted.nhx"
params.phenotype = path_file + "data/rhodopsin/spp_ismarine.tsv"
params.align = path_file + "data/rhodopsin/lessgappy.align_reroot.fa"
params.method = "MAP" // DOWNPASS, ACCTRAN
params.resdir= path_file + "results/" //do not forget to change the resdir

/*
//HIV real data
params.tree = path_file + "data/HIV_africa/pcoc_file/rooted_align.treefile" <---- to check
params.phenotype = path_file + "data/HIV_africa/pcoc_file/cvgt_phenotype.tsv"
params.align = path_file + "data/HIV_africa/pcoc_file/align.noCRF.jphmm.aa.fa"
params.method = "MAP" // DOWNPASS, ACCTRAN
params.resdir= path_file + "results/" //do not forget to change the resdir
*/

/*
//sedges c3c4
params.tree = path_file + "data/c3c4/cyp_coding.phy_phyml_tree.txt"
params.phenotype = path_file + "data/c3c4/spp_cvgt.tsv"
params.align = path_file + "data/c3c4/cyp_coding.aa.coor_mays.fa"
params.method = "MAP" // DOWNPASS, ACCTRAN
params.resdir=path_file + "results/" //do not forget to change the resdir
*/

tree = file(params.tree)
pheno = file(params.phenotype)
align = file(params.align)
method = params.method

resdir=file(params.resdir)
resdir.with {mkdirs()}

process acr_pastml{
    label 'pastml'

    publishDir "${resdir}", mode: 'copy'
    input:
    file tree
    file pheno
    val method

    output:
    tuple file("work_pastml/*nwk"), file("rhodo.*pastml.ML.out") into pastml_ML_out
    
    //rate sed numerotation from 1 ($line -1)
    //input numerotation from 0
    //create parameter file for pastml including rate per site (scaling factor) for each site and frequencies of amino acids
    //run pastml 
    //remove some temp files
    //retrive marginal proba for root (sed -n 2p)

    //rhodo: "cvgt" for fresh brackish, "ismarine" for marine
    //ID="cvgt"
    shell:
    '''
    ID="ismarine"  
    pastml --threads !{task.cpus} -t !{tree} -d !{pheno} --prediction_method !{method} -m F81 -c $ID -o rhodo.!{method}.pastml.ML.out --work_dir work_pastml
    '''
}

process scenario{
    label 'python'

    publishDir "${resdir}", mode: 'copy'
    input:
    file tree
    file pheno
    output:
    tuple file("pcoc_tree_rerooted.nw"), file("scenario_pcoc_noacr.txt") into PcocScenario
    shell:
    '''
    scenario_pcoc.py !{tree} !{pheno} ismarine
    ''' 
}

process pcoc_run{
    label 'pcoc_old'
    publishDir "${resdir}", mode: 'copy'
    input:
    tuple file(pcoc_tree), file(scenario) from PcocScenario
    file align
    file tree
    output:
    file 'output_pcoc'
    
    shell:
    '''
    S=`cat !{scenario}`
    pcoc_det.py -t !{tree} -aa !{align} -o output_pcoc -m $S -f 0.8 --max_gap_allowed 20 -CATX_est 10 --gamma 
    #pcoc_det.py -t !{tree} -aa !{align} -cpu 4 -o output_pcoc -m $S -f 0.8 -f_oc 0.5 --max_gap_allowed 100 -est_profiles C10 --gamma 
    ''' 
}
