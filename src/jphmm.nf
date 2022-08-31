#!/usr/bin/env nextflow

path_local = "U:/users/mamorel/Projet_Convergence/condor-analysis/"


//change to your local directory with JPHMM
pathjphmm = "/pasteur/zeus/projets/p01/Evolbioinfo/users/mamorel/bin/jpHMM/"

align = path_local + "data/HIV_africa/preprocessing/noCRF.fa"

params.resdir="path_local/results/outputjphmm/"

resdir=file(params.resdir)
resdir.with {mkdirs()}

SplitChannel = Channel.fromPath(align).splitFasta( by: 5, file:true )

process jphmm {
    publishDir "${resdir}", mode: 'copy'
    label 'jphmm'
    input:
    file align from SplitChannel
    output: 
    file ("output/recombination_*.txt") 
    shell:
    '''
    pathjphmm="/pasteur/zeus/projets/p01/Evolbioinfo/users/mamorel/bin/jpHMM/"
    jpHMM -s !{align} -v HIV -a ${pathjphmm}priors/emissionPriors_HIV.txt -b ${pathjphmm}priors/transition_priors.txt -i ${pathjphmm}input/HIV_alignment.fas
    mv output/recombination.txt output/recombination_!{align}.txt
    '''
}


