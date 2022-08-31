
#### rhodospine
##conjunction
/pasteur/zeus/projets/p01/Evolbioinfo/users/mamorel/Projet_Convergence/Data/args/processed/rhodopsin/lessgappy.align_reroot.fa
/pasteur/zeus/projets/p01/Evolbioinfo/users/mamorel/Projet_Convergence/Data/args/processed/rhodopsin/FADE/tree_rerooted_foreground_conjunction.nhx
sbatch --qos bioevo -p bioevo --wrap="./hyphy fade --alignment rhodopsin/lessgappy.align_reroot.fa --tree rhodopsin/FADE/tree_rerooted_foreground_conjunction.nhx --output rhodopsin/FADE/foreground_aa.LG.FADE.json --branches Foreground --model LG --cache rhodopsin/FADE/foreground_aa.mtZOA.FADE.cache"



### c3c4
sbatch --qos bioevo -p bioevo --mem-per-cpu 10G --wrap="./hyphy fade --alignment C3C4/cyp_coding.aa.coor_mays.fa --tree C3C4/FADE/cyp_coding.phy_phyml_tree_foreground.txt --output C3C4/FADE/foreground_aa.FADE.json --branches Foreground"
sbatch --qos bioevo -p bioevo --mem-per-cpu 10G --wrap="./hyphy fade --alignment c3c4/new_homologues/align.subseq.alignment.aa.fasta --tree c3c4/new_homologues/align.subseq.alignment.aa.foregroundc3andc4.treefile --output c3c4/new_homologues/foreground_aa_c3c4.FADE.json --cache c3c4/new_homologues/foreground_aa_c3c4.FADE.cache --branches Foreground --model JTT"
sbatch --qos bioevo -p bioevo --mem-per-cpu 10G --wrap="./hyphy fade --alignment c3c4/no_c3c4_clade/cyp_coding_noc3c4.aa.coor_mays.fa --tree c3c4/no_c3c4_clade/foreground_aa_c4_pheno.FADE.cache --branches Foreground --model JTT"


### HIV africa
sbatch --qos bioevo -p bioevo --mem-per-cpu 10G --wrap="./hyphy fade --model HIVBm --alignment africa/recomb_jphmm/align.noCRF.jphmm_outgroup.aa.fa --tree africa/recomb_jphmm/foreground_tips_conjunction_aa.tree --output africa/recomb_jphmm/foreground_aa_conjunction.FADE.json --branches Foreground"