# :warning: Repository moved to a [new location](https://github.com/evolbioinfo/condor-analysis/) :warning:

# Accurate Detection of Convergent Mutations in Large Protein Alignments with ConDor


### Marie MOREL, Frédéric LEMOINE, Anna ZHUKOVA and Olivier GASCUEL

In this repository you can find data and files used in the ConDor paper. \
ConDor is a workflow to detect convergent evolution at the scale of the mutation in protein alignments. 
It is available from a web service located at https://condor.pasteur.cloud/

## Help
* data folder: data described in the Materials and Methods section, including sequence alignments and phylogenetic trees.  
* results folder: csv files corresponding to the outputs provided by ConDor, as well as the results from FADE and PCOC.
    * __Note__: For Rhodopsin data, amino-acid positions presented in results files `results/rhodopsin` have a -12 offset compared to numbering presented in the manuscript.
* src folder: scripts used to analyze the data and obtain tables and figures as in the ConDor Paper. You can also find the commands used to run the various tools (FADE, PCOC, JPHMM...).



