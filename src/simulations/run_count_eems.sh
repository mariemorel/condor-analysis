set -ue

PREFIX=../
MATRICES=$PREFIX/condor/assets/protein_model.txt
TREE=$PREFIX/data/cyp_coding.phy_phyml_tree.txt
REALALIGN=$PREFIX/data/cyp_coding.aa.coor_mays.fa
OUTGROUP=$PREFIX/data/outgroup.txt
PHENOTYPE=$PREFIX/data/besnard2009_convergent_species.txt
CONVPOS=$PREFIX/data/convmuts_12.txt
SIMU=$PREFIX/simulation.nf
SIMUCONF=$PREFIX/simulation.config
CONDOR=$PREFIX/condor/condor.nf
COMPUTETABLE=$PREFIX/bin/compute_tables_12.pl
OUTDIR=res_simu_12
RESDIR=res_cond_12

### Simulations + CONDOR for comparing Condor estimated EEMs with Simulated EEMs
################################################################################

for i in {1..10}
do
    # Simulation
    nextflow run -c $SIMUCONF -profile pasteurcluster $SIMU --matrices $MATRICES --outdir ${OUTDIR}_nb${i} --pastml true --snag false --tree $TREE --msa $REALALIGN --outgroup $OUTGROUP --convpos $CONVPOS
    # CONDOR
    nextflow run -profile pasteurcluster $CONDOR --align ${OUTDIR}_nb${i}/simulated_alignment_jtt.fasta --tree ${OUTDIR}_nb${i}/cyp_coding.aa.coor_mays_gaps.fasta_reroot_keep.treefile --outgroup $OUTGROUP --phenotype $PHENOTYPE --resdir ${RESDIR}_nb${i} --model 'JTT+R3' --optimize true --matrices $MATRICES --nb_simu 10 --min_seq 1 --min_eem 1 --freqmode 'Fmodel' --branches condor --correction holm --alpha 0.1 --bayes 2 -resume
    # Model Violation
    nextflow run -profile pasteurcluster $CONDOR --align ${OUTDIR}_nb${i}/simulated_alignment_jtt.fasta --tree ${OUTDIR}_nb${i}/cyp_coding.aa.coor_mays_gaps.fasta_reroot_keep.treefile --outgroup $OUTGROUP --phenotype $PHENOTYPE --resdir ${RESDIR}_LG_G4_nb${i} --model 'LG+G4' --optimize true --matrices $MATRICES --nb_simu 10 --min_seq 1 --min_eem 1 --freqmode 'Fmodel' --branches condor --correction holm --alpha 0.1 --bayes 2

done

# High diversity Tree scaling x 3.0
for i in {1..10}
do
    nextflow run -c $SIMUCONF -profile pasteurcluster $SIMU --matrices $MATRICES --outdir ${OUTDIR}_scale3.0_nb${i}  --pastml true --snag false --tree $TREE --msa $REALALIGN --outgroup $OUTGROUP --convpos $CONVPOS --scale 3
    nextflow run -profile pasteurcluster $CONDOR --align ${OUTDIR}_scale3.0_nb${i}/simulated_alignment_jtt.fasta --tree ${OUTDIR}_scale3.0_nb${i}/cyp_coding.aa.coor_mays_gaps.fasta_reroot_keep.treefile --outgroup $OUTGROUP --phenotype $PHENOTYPE --resdir ${RESDIR}_scale3.0_nb${i} --model 'JTT+R3' --optimize true --matrices $MATRICES --nb_simu 10 --min_seq 1 --min_eem 1 --freqmode 'Fmodel' --branches condor --correction holm --alpha 0.1 --bayes 2
done

# Low diversity Tree scaling x 0.33
for i in {1..10}
do
    nextflow run -c $SIMUCONF -profile pasteurcluster $SIMU --matrices $MATRICES --outdir ${OUTDIR}_scale0.33_nb${i}  --pastml true --snag false --tree $TREE --msa $REALALIGN --outgroup $OUTGROUP --convpos $CONVPOS --scale 0.33
    nextflow run -profile pasteurcluster $CONDOR --align ${OUTDIR}_scale0.33_nb${i}/simulated_alignment_jtt.fasta --tree ${OUTDIR}_scale0.33_nb${i}/cyp_coding.aa.coor_mays_gaps.fasta_reroot_keep.treefile --outgroup $OUTGROUP --phenotype $PHENOTYPE --resdir ${RESDIR}_scale0.33_nb${i} --model 'JTT+R3' --optimize true --matrices $MATRICES --nb_simu 10 --min_seq 1 --min_eem 1 --freqmode 'Fmodel' --branches condor --correction holm --alpha 0.1 --bayes 2
done

# Extract CONDOR and Simulated EEMs Counts
for i in {1..10}; do ../bin/compare_eems_simu_condor.pl res_simu_12_nb${i}/simulated_alignment_eems.tsv res_cond_12_nb${i}/tested_results.tsv | tail -n+2 | awk -v i=$i '{print i"\t" $0}'; done > comp_simu_cond.txt
for i in {1..10}; do ../bin/compare_eems_simu_condor.pl res_simu_12_scale0.33_nb${i}/simulated_alignment_eems.tsv res_cond_12_scale0.33_nb${i}/tested_results.tsv | tail -n+2 | awk -v i=$i '{print i"\t" $0}'; done > comp_simu_cond_scale0.33.txt
for i in {1..10}; do ../bin/compare_eems_simu_condor.pl res_simu_12_scale3.0_nb${i}/simulated_alignment_eems.tsv res_cond_12_scale3.0_nb${i}/tested_results.tsv | tail -n+2 | awk -v i=$i '{print i"\t" $0}'; done > comp_simu_cond_scale3.00.txt
for i in {1..10}; do ../bin/compare_eems_simu_condor.pl res_simu_12_nb${i}/simulated_alignment_eems.tsv res_cond_12_LG_G4_nb${i}/tested_results.tsv | tail -n+2 | awk -v i=$i '{print i"\t" $0}'; done > comp_simu_cond_model.txt

# Count number of different amino acids per sites
for i in {1..10}; do goalign stats char -i res_simu_12_nb${i}/simulated_alignment_jtt.fasta --per-sites | awk -v rep=$i '{if(NR==1){print "Replicate\t" $0 "\tNbAA"}else{if($1>452 && $1<908){naa=0;for(i=2;i<=NF;i++){if($i>0){naa++}};print rep "\t" $0 "\t" naa}}}' > count_diff_aa_nb${i}.txt ; done
for i in {1..10}; do goalign stats char -i res_simu_12_scale0.33_nb${i}/simulated_alignment_jtt.fasta --per-sites | awk -v rep=$i '{if(NR==1){print "Replicate\t" $0 "\tNbAA"}else{if($1>452 && $1<908){naa=0;for(i=2;i<=NF;i++){if($i>0){naa++}};print rep "\t" $0 "\t" naa}}}' > count_diff_aa_scale0.33_nb${i}.txt ; done
for i in {1..10}; do goalign stats char -i res_simu_12_scale3.0_nb${i}/simulated_alignment_jtt.fasta --per-sites | awk -v rep=$i '{if(NR==1){print "Replicate\t" $0 "\tNbAA"}else{if($1>452 && $1<908){naa=0;for(i=2;i<=NF;i++){if($i>0){naa++}};print rep "\t" $0 "\t" naa}}}' > count_diff_aa_scale3.00_nb${i}.txt ; done

# Sum of the nb amino acids per replicate
rm -f count_diff_aa.txt count_diff_aa_scale0.33.txt count_diff_aa_scale3.0.txt
for i in {1..10}; do goalign stats char -i res_simu_12_nb${i}/simulated_alignment_jtt.fasta --per-sites | awk -v rep=$i 'BEGIN{total=0}{if(NR>1){if($1>452 && $1<908){naa=0;for(i=2;i<=NF;i++){if($i>0){naa++}};total+=naa}}}END{print rep "\t" total}' >> count_diff_aa.txt ; done
for i in {1..10}; do goalign stats char -i res_simu_12_scale0.33_nb${i}/simulated_alignment_jtt.fasta --per-sites | awk -v rep=$i 'BEGIN{total=0}{if(NR>1){if($1>452 && $1<908){naa=0;for(i=2;i<=NF;i++){if($i>0){naa++}};total+=naa}}}END{print rep "\t" total}' >> count_diff_aa_scale0.33.txt ; done
for i in {1..10}; do goalign stats char -i res_simu_12_scale3.0_nb${i}/simulated_alignment_jtt.fasta --per-sites | awk -v rep=$i 'BEGIN{total=0}{if(NR>1){if($1>452 && $1<908){naa=0;for(i=2;i<=NF;i++){if($i>0){naa++}};total+=naa}}}END{print rep "\t" total}' >> count_diff_aa_scale3.0.txt ; done

# Run R script to compute metrics
# Rscript $PREFIX/bin/compare_EEMs.R
