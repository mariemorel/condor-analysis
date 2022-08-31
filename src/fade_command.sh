#FADE version 0.2

### Sedges c3c4
sbatch --qos bioevo -p bioevo --mem-per-cpu 10G --wrap="./hyphy fade \\
--alignment data/c3c4/cyp_coding.aa.coor_mays.fa \\
--tree data/c3c4/cyp_coding.phy_phyml_tree_foreground.txt \\
--output results/c3c4/FADE/foreground_aaJTT.FADE.json \\
--cache results/c3c4/FADE/foreground_aaJTT.FADE.cache \\
--branches Foreground --model JTT"

sbatch --qos bioevo -p bioevo --mem-per-cpu 10G --wrap="./hyphy fade \\
--alignment data/c3c4/pheno_noc3c4/cyp_coding_noc3c4.aa.coor_mays.fa \\
--tree data/c3c4/pheno_noc3c4/cyp_coding_noc3c4_foreground.phy_phyml_tree.txt
--cache results/c3c4/FADE/foreground_aa_noc3c4_pheno.FADE.cache \\
--output results/c3c4/FADE/foreground_aa_noc3c4_pheno.FADE.json \\
--branches Foreground --model JTT"

### HIV africa
sbatch --qos bioevo -p bioevo --mem-per-cpu 10G --wrap="./hyphy fade \\
--model HIVBm \\
--alignment data/HIV_africa/align.noCRF.jphmm_outgroup.aa.fa \\
--tree data/HIV_africa/foreground_tips_conjunction_aa.tree \\
--output results/HIV_africa/FADE/foreground_aa_conjunction_HIVb.FADE.json \\
--branches Foreground"

#model misspecification JTT
sbatch --qos bioevo -p bioevo --mem-per-cpu 10G --wrap="./hyphy fade \\
--model JTT \\
--alignment data/HIV_africa/align.noCRF.jphmm_outgroup.aa.fa \\
--tree data/HIV_africa/foreground_tips_conjunction_aa.tree \\
--output results/HIV_africa/FADE/foreground_aa_conjunction_JTT.FADE.json \\
--branches Foreground"


#### rhodospine
##fresh/brackish
sbatch --qos bioevo -p bioevo --wrap="./hyphy fade \\
--alignment data/rhodopsin/lessgappy.align_reroot.fa \\
--tree data/rhodopsin/tree_rerooted_foreground_conjunction.nhx \\
--output results/rhodopsin/FADE/foreground_aa_conjunction_brackish.LG.FADE.json \\
--branches Foreground \\
--model LG \\
--cache results/rhodopsin/FADE/foreground_aa_conjunction_brackish.LG.FADE.cache"

###marine
sbatch --qos bioevo -p bioevo --wrap="./hyphy fade \\
--alignment data/rhodopsin/lessgappy.align_reroot.fa \\
--tree data/rhodopsin/tree_rerooted_foreground_marine_conjunction.nhx \\
--output results/rhodopsin/FADE/foreground_aa_conjunction_marine.LG.FADE.json \\
--branches Foreground \\
--model LG \\
--cache results/rhodopsin/FADE/foreground_aa_conjunction_marine.LG.FADE.cache"
