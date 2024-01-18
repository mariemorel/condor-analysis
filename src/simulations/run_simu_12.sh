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

### Simulations + CONDOR
########################

## 10 times SIMU+CONDOR
for i in {1..10}
do
    # Simulation
    nextflow run -c $SIMUCONF -profile pasteurcluster $SIMU --matrices $MATRICES --outdir ${OUTDIR}_nb${i} --pastml true --snag false --tree $TREE --msa $REALALIGN --outgroup $OUTGROUP --convpos $CONVPOS
    # True run
    nextflow run -profile pasteurcluster $CONDOR --align ${OUTDIR}_nb${i}/simulated_alignment_jtt_convergent.fasta --tree ${OUTDIR}_nb${i}/cyp_coding.aa.coor_mays_gaps.fasta_reroot_keep.treefile --outgroup $OUTGROUP --phenotype $PHENOTYPE --resdir ${RESDIR}_nb${i} --model 'JTT+R3' --optimize true --matrices $MATRICES --nb_simu 10000 --min_seq 3 --min_eem 3 --freqmode 'Fmodel' --branches condor --correction holm --alpha 0.1 --bayes 2
    # Model Violation : LG+G4
    nextflow run -profile pasteurcluster $CONDOR --align ${OUTDIR}_nb${i}/simulated_alignment_jtt_convergent.fasta --tree ${OUTDIR}_nb${i}/cyp_coding.aa.coor_mays_gaps.fasta_reroot_keep.treefile --outgroup $OUTGROUP --phenotype $PHENOTYPE --resdir ${RESDIR}_LG_G4_nb${i} --model 'LG+G4' --optimize true --matrices $MATRICES --nb_simu 10000 --min_seq 3 --min_eem 3 --freqmode 'Fmodel' --branches condor --correction holm --alpha 0.1 --bayes 2
done

## High diversity Tree Scaling x3.0
for i in {1..10}
do
    nextflow run -c $SIMUCONF -profile pasteurcluster $SIMU --matrices $MATRICES --outdir ${OUTDIR}_scale3.0_nb${i}  --pastml true --snag false --tree $TREE --msa $REALALIGN --outgroup $OUTGROUP --convpos $CONVPOS --scale 3
    nextflow run -profile pasteurcluster $CONDOR --align ${OUTDIR}_scale3.0_nb${i}/simulated_alignment_jtt_convergent.fasta --tree ${OUTDIR}_scale3.0_nb${i}/cyp_coding.aa.coor_mays_gaps.fasta_reroot_keep.treefile --outgroup $OUTGROUP --phenotype $PHENOTYPE --resdir ${RESDIR}_scale3.0_nb${i} --model 'JTT+R3' --optimize true --matrices $MATRICES --nb_simu 10000 --min_seq 3 --min_eem 3 --freqmode 'Fmodel' --branches condor --correction holm --alpha 0.1 --bayes 2
done

## Low diversity : Scaling x0.33
for i in {1..10}
do
    nextflow run -c $SIMUCONF -profile pasteurcluster $SIMU --matrices $MATRICES --outdir ${OUTDIR}_scale0.33_nb${i}  --pastml true --snag false --tree $TREE --msa $REALALIGN --outgroup $OUTGROUP --convpos $CONVPOS --scale 0.33
    nextflow run -profile pasteurcluster $CONDOR --align ${OUTDIR}_scale0.33_nb${i}/simulated_alignment_jtt_convergent.fasta --tree ${OUTDIR}_scale0.33_nb${i}/cyp_coding.aa.coor_mays_gaps.fasta_reroot_keep.treefile --outgroup $OUTGROUP --phenotype $PHENOTYPE --resdir ${RESDIR}_scale0.33_nb${i} --model 'JTT+R3' --optimize true --matrices $MATRICES --nb_simu 10000 --min_seq 3 --min_eem 3 --freqmode 'Fmodel' --branches condor --correction holm --alpha 0.1 --bayes 2
done

## Phenotypic noise
for i in {1..10}
do
    mkdir -p ${RESDIR}_noisy8_nb${i}
    # convergent species -8  to keep
    shuf $PHENOTYPE| tail -n+9 > ${RESDIR}_noisy8_nb${i}/phenotype.txt
    # 8 random additional convergent species
    gotree labels -i $TREE | grep -v -f $PHENOTYPE | shuf | head -n 8 >> ${RESDIR}_noisy8_nb${i}/phenotype.txt
    # Now condor on this noisy convergent annotation dataset
    nextflow run -profile pasteurcluster $CONDOR --align ${OUTDIR}_nb${i}/simulated_alignment_jtt_convergent.fasta --tree ${OUTDIR}_nb${i}/cyp_coding.aa.coor_mays_gaps.fasta_reroot_keep.treefile --outgroup $OUTGROUP --phenotype ${RESDIR}_noisy8_nb${i}/phenotype.txt --resdir ${RESDIR}_noisy8_nb${i} --model 'JTT+R3' --optimize true --matrices $MATRICES --nb_simu 10000 --min_seq 3 --min_eem 3 --freqmode 'Fmodel' --branches condor --correction holm --alpha 0.1 --bayes 2
done

## More phenotypic noise
for i in {1..10}
do
    mkdir -p ${RESDIR}_noisy10_nb${i}
    # convergent species -10  to keep
    shuf $PHENOTYPE| tail -n+11 > ${RESDIR}_noisy10_nb${i}/phenotype.txt
    # 4 random additional convergent species
    gotree labels -i $TREE | grep -v -f $PHENOTYPE | shuf | head -n 10 >> ${RESDIR}_noisy10_nb${i}/phenotype.txt
    # Now condor on this noisy convergent annotation dataset
    nextflow run -profile pasteurcluster $CONDOR --align ${OUTDIR}_nb${i}/simulated_alignment_jtt_convergent.fasta --tree ${OUTDIR}_nb${i}/cyp_coding.aa.coor_mays_gaps.fasta_reroot_keep.treefile --outgroup $OUTGROUP --phenotype ${RESDIR}_noisy10_nb${i}/phenotype.txt --resdir ${RESDIR}_noisy10_nb${i} --model 'JTT+R3' --optimize true --matrices $MATRICES --nb_simu 10000 --min_seq 3 --min_eem 3 --freqmode 'Fmodel' --branches condor --correction holm --alpha 0.1 --bayes 2
done

$COMPUTETABLE ${RESDIR}_nb* > stats_12_replicates.tsv
$COMPUTETABLE ${RESDIR}_LG_G4_nb* > stats_12_lg_g4_replicates.tsv
$COMPUTETABLE ${RESDIR}_scale3.0_nb* > stats_12_scale3.0.tsv
$COMPUTETABLE ${RESDIR}_scale0.33_nb* > stats_12_scale0.33.tsv
$COMPUTETABLE ${RESDIR}_noisy8_nb* > stats_12_noisy8_replicates.tsv
$COMPUTETABLE ${RESDIR}_noisy10_nb* > stats_12_noisy10_replicates.tsv
