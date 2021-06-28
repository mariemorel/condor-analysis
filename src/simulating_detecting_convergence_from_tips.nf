path_file = "/pasteur/zeus/projets/p01/Evolbioinfo/users/mamorel/Projet_Convergence/Data/"

params.tree = path_file+"/preprocessing/reroot.B.pol.aa.norecomb.noPRandRTdrm.treefile" 
params.root = path_file+"/preprocessing/5_RT_hxb2_naturepaper.txt"

params.ref_freq = "0.060490222, 0.066039665, 0.044127815, 0.042109048, 0.020075899, 0.053606488, 0.071567447, 0.072308239, 0.022293943, 0.069730629, 0.098851122, 0.056968211, 0.019768318, 0.028809447, 0.046025282, 0.05060433, 0.053636813, 0.033011601, 0.028350243, 0.061625237"
params.rates_init = path_file+"/preprocessing/5RT_hxb2_naturepaper.rate"
params.mode = 'HIVb' //never changes
params.simulation_model = path_file+"args/raw/HIVb_phyml.model" //changes to JTT_phyml.model or HIV HIVb_phyml.model
params.iqtreemode = 'HIVb+G6' //changes to JTT+G6 or HIVb+G6
//params.sdrms = path_file+"args/processed/naturepaper/sequences_with_most_common_RT_SDRMs.txt"
params.sdrms = path_file+"args/processed/naturepaper/sequences_with_DRMs_RT.txt"
params.outgroup = path_file+"trees/processed/naturepaper/B.pol.aa.woutgroup.norecomb.maskeddrm.outgroup"
params.matrix = path_file+"args/raw/HIV_phyml.pastml_matrix" //changes to JTT.pastml_matrix or HIV_phyml.pastml_matrix
params.reconstruction = "false" //do we reconstruct the tree from scratch once simulated or do we only reoptimize ?

params.resdir=path_file+"results/naturepaper/synthetic_HIVG6maskedDRMsparse_12/" //do not forget to change the resdir

params.nb_simu = 10000
params.length = 1250 // length of the root = 5*250
params.nb_tips = 3587

tree = file(params.tree)
root = file(params.root)
ref_freq = params.ref_freq
rates_init = file(params.rates_init)
mode = params.mode
simulation_model = file(params.simulation_model)
iqtreemode = params.iqtreemode
sdrm = file(params.sdrms)
outgroup = file(params.outgroup)
matrix = file(params.matrix)
reconstruction = params.reconstruction


nb_simu = params.nb_simu
length = params.length
nb_tips = params.nb_tips


resdir=file(params.resdir)
resdir.with {mkdirs()}


//create an alignment corresponding to this tree
//with rates following gamma shape 
process reference_align_simulation {
    conda '/pasteur/sonic/homes/mamorel/miniconda3/envs/jupyter-notebook'
    input:
    file root
    file tree
    file rates_init
    val ref_freq
    val mode  
    output:
    file("out_simulated_5RT*.fasta") into SimulatedChannel
    shell:
    '''
    R=`cat !{root}`
    for i in !{ref_freq} ; do echo ${i%,}  >> freq ; done
    simulator.py ${R} !{tree} !{rates_init} freq !{mode} 5RT_freerate_simu
    '''
}


// add convergence in alignement following amount of the true DRMs or not. 
process add_convergence{
    conda '/pasteur/sonic/homes/mamorel/miniconda3/envs/jupyter-notebook'
    
    publishDir "${resdir}", mode: 'copy'
    input :
    file(align) from SimulatedChannel
    file root
    file sdrm

    output:
    file "increase*" into ConvergentChannel, TipsChannel
    file 'positions*'

    shell:
    '''
    R=`cat !{root}`
    add_convergence_DRMtips.py ${R} !{align} !{sdrm} 
    '''
}


process construct_tree {
    publishDir "${resdir}", pattern: "*_rate", mode: 'copy'
    publishDir "${resdir}", pattern: "*.iqtree", mode: 'copy'
    input:
    file(align) from ConvergentChannel
    file tree
    val iqtreemode
    val length
    output:
    file "*.treefile" into TreeChannel 
    tuple val(length), file(align), file ("reestimated_rate") into RatesparamsChannel
    file "reestimated_rate" into  RatesChannel
    file "*.iqtree" 
    
    shell:  
     if( reconstruction == 'true' )

//starting tree 
    '''
    iqtree -m !{iqtreemode} -nt AUTO -s !{align} -t !{tree} -wsr -pre align 
    tail -n+10 align.rate | cut -f 2 > reestimated_rate
    len=`wc -l reestimated_rate | cut -d " " -f 1`
    if [ "$len" -eq "0" ]; then for i in {1..!{length}} ; do echo 1 >> reestimated_rate ;done ;  fi
    
    '''
    //gotree divide -i scale_rooted_tree.nw -o split_tree 

     else if (reconstruction == 'false')
//fixed tree no tree search performed
    '''
    iqtree -m !{iqtreemode} -nt AUTO -s !{align} -te !{tree} -wsr -pre align
    tail -n+10 align.rate | cut -f 2 > reestimated_rate
    len=`wc -l reestimated_rate | cut -d " " -f 1`
    if [ "$len" -eq "0" ] ; then for i in {1..!{length}} ; do echo 1 >> reestimated_rate ;done ;  fi
    '''
}


//rename it, and reroot it, make input file for seqgen
//export name tree for my simulator 
process tree_rename{
    publishDir "${resdir}", mode: 'copy'
    input:
    file tree from TreeChannel
    file outgroup
    output:
    file "named_tree" into SimulatorChannel, NamedtreeChannel
    file "rooted_*" into RootedtreeChannel
    //file "named_*" into MySimulatorChannel
    shell:
    ''' 
    gotree reroot outgroup -i !{tree} -l !{outgroup} -o rooted_!{tree}
    gotree rename --internal --tips=false --auto -i rooted_!{tree} -o named_tree
    '''
}

process frequencies_file{
    input:
    val ref_freq
    output:
    file "frequencies.txt" into FreqChannel, FreqparamsChannel, FrequenciesChannel
    shell:
    '''
    for i in !{ref_freq} ; do echo ${i%,}  >> freq ; done
    for i in A R N D C Q E G H I L K M F P S T W Y V ; do printf $i"\n" ; done > amino
	paste amino freq > frequencies.txt
    '''    
}

ParametersChannel = RatesparamsChannel.merge(FreqparamsChannel)

process pars_align_file{
    conda '/pasteur/sonic/homes/mamorel/miniconda3/envs/jupyter-notebook'
    input : 
    tuple val(length), file(align), file(rates), file(freq) from ParametersChannel

    output:
    tuple val(length), file(rates), file(freq), file('*pastml_input.tsv.gz') into Pastml_align mode flatten
  
    shell:
    '''
    pars_fasta.py !{align} !{length} refalign_
    '''
}


process acr_pastml{
    conda '/pasteur/sonic/homes/mamorel/miniconda3/envs/pastmlmatrix'
    publishDir "${resdir}", pattern: "work_pastml/named*.nw",  mode: 'copy'
    input:
    tuple val(length), file(rate), file(freq), file(input) from Pastml_align
    file tree from RootedtreeChannel
    file matrix

    output:
    tuple val(length), file(rate), file("work_pastml/named.*nwk"), file("*pastml.ML.out.gz"), file("marginal_root.txt") into pastml_ML_out

    shell:
    '''
    align="!{input}"
    for id in {1..!{length}} ; do 
    R=`sed -n "${id}"p !{rate}` 
    echo -e 'parameter\tvalue' > parameter_$((${id}-1)) 
    cat !{freq} >> parameter_$((${id}-1))
    echo -e "scaling_factor\t${R}" >> parameter_$((${id}-1)) ; done
    gunzip -c !{input} > ${align%.*.*}.input

    VAR=`for id in {1..!{length}}; do echo parameter_$((${id}-1)); done`
    ID=`for id in {1..!{length}}; do echo $((${id}-1)); done`
    MATRIX=`for id in {1..!{length}}; do echo !{matrix}; done`
    pastml -t !{tree}  -d ${align%.*.*}.input --prediction_method MAP -m CUSTOM_RATES --rate_matrix $MATRIX -c $ID -o ${align%.*.*}.pastml.ML.out --work_dir work_pastml --parameter $VAR
    gzip ${align%.*.*}.pastml.ML.out
    rm work_pastml/params*.tab
    rm parameter_*
    for i in `ls work_pastml/marginal_probabilities.character_*` ; do echo ${i//[^0-9]/} ; sed -n 2p $i ;  done > marginal_root.txt 

    '''
}

process pre_count{
    conda '/pasteur/sonic/homes/mamorel/miniconda3/envs/jupyter-notebook'
    publishDir "${resdir}", mode: 'copy'
    input : 
    tuple val(length), file(rate), file(tree), file(pastml_acr), file(marginal_root) from pastml_ML_out
  
    output : 
    tuple file(rate), file(tree), file("*pastml_acr.fasta") into python_count
    file "reconstructed_root" into Root_seq
    file "*marginal_posterior.txt" into Marginal_root

    shell:
    '''
    pastml_fasta.py !{pastml_acr} !{length} !{marginal_root} 5RT_
    '''
    //gunzip -c '1.out.gz' > pastml_output
    //for i in {2..!{length}};  do gunzip -c ${i}.out.gz | cut -f 2 > temp && paste pastml_output temp > test2 && mv test2 pastml_output ;  done
    //rm temp

}

//my simulator of sequence in python, using root, nb simulations and the ROOTED tree.
process simulator {
    //errorStrategy 'retry'
    //maxRetries 3
    conda '/pasteur/sonic/homes/mamorel/miniconda3/envs/jupyter-notebook'
    publishDir "${resdir}", pattern: "count*.tsv.gz", mode: 'copy'
    input : 
    each x from 1..length
    val nb_simu
    file simulation_model
    //val simufreq
    file (rates) from RatesChannel
    file (freq) from FreqChannel
    file (named_tree) from SimulatorChannel
    file (root) from Marginal_root

    output : 
    //tuple val(x), file("count*.tsv.gz") into MysimulationsChannel.collect() // We take all the counts into a single Channel
    file("count*npz") into MysimulationsChannel
    shell:
    '''
    sed -n "/^$((!{x}-1))\t/p" !{root} | cut -f 2 > root.txt
    rate=`sed -n "!{x}p" !{rates}` # sed numerotation from 1
    
    output="!{x}_!{named_tree}_"
    simulator_counting_rates_from_root.py root.txt !{named_tree} !{nb_simu} ${output} ${rate} !{freq} !{simulation_model} 
    '''
}

Collect_simulations = MysimulationsChannel.collect()

//fasta file into a tab separated table understandable by pastml
//only input needed is align 
//56 different alignments 

process count_apparitions{
    conda '/pasteur/sonic/homes/mamorel/miniconda3/envs/jupyter-notebook'
    input : 
    tuple file(rate), file(tree), file (align) from python_count
    output : 
    tuple file(rate), file(align), file ("*substitutions_even_root.tsv"), file("*substitutions_aa_tips_per_base.tsv") into Subscribe_matrices, Ref_couting
    
    shell:
    '''
    count_substitutions_from_tips.py !{align} !{tree} !{length}
    '''
}

Subscribe_matrices.subscribe{rate, align, freqs, substitutions ->  freqs.copyTo(file("${resdir}").resolve('ref_substitutions.txt'));}
 
// should 
process conclude_convergence{
    publishDir "${resdir}", mode: 'copy'
    conda '/pasteur/sonic/homes/mamorel/miniconda3/envs/jupyter-notebook'
    input: 
    file simulation_model
    file (tipsalign) from TipsChannel
    file (freq) from FrequenciesChannel
    tuple file (rate), file (acralign), file(ref_matrix), file(substitutions) from Ref_couting 
    file (counts) from Collect_simulations
    file(root) from Root_seq
    val nb_simu

    //named*.phy
     
    output:
    file("all_results_metrics.tsv") 

    shell:
    '''
    R=`cat !{root}`
    convergent_substitutions_pvalue_simu.py ${R} !{rate} !{tipsalign} !{ref_matrix} !{substitutions} !{nb_simu} !{freq} !{simulation_model}
    '''
}


