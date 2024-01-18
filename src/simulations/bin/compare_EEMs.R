setwd("/pasteur/zeus/projets/p01/Evolbioinfo/users/flemoine/projets/2023_10_13_Condor_Revisions/results_compare_eems")

eemanalysis("comp_simu_cond.txt")
eemanalysis("comp_simu_cond_model.txt")
eemanalysis("comp_simu_cond_scale0.33.txt")
eemanalysis("comp_simu_cond_scale3.00.txt")


eemanalysis=function(infile){
  simucondor=read.table(infile,header=F)
  colnames(simucondor) = c("rep","pos","mut","simu","condor")
  
  print("correlation simu + condor of all mutations")
  plot(simucondor$simu,simucondor$condor)
  print(cor(simucondor$simu,simucondor$condor))
  
  simucondorsup2simu=simucondor[simucondor$simu>2,]
  simucondorsup2condor=simucondor[simucondor$condor>2,]
  simucondorsup2all=simucondor[simucondor$simu>2 & simucondor$condor>2,]
  #plot(simucondorsup2simu$simu,simucondorsup2simu$condor)
  print("Distribution of SIMU EEMs - CONDOR EEMs among mutations that have >=3 Simulated & CONDOR EEMs")
  print(table(simucondorsup2all$simu-simucondorsup2all$condor))
  print("Number of mutations that have >=3 Simulated EEMs")
  print(dim(simucondorsup2simu))
  print("Number of mutations that have >=3 CONDOR EEMs")
  print(dim(simucondorsup2condor))
  print("Number of mutations that have >=3 Simulated EEMs and >=3 Condor EEMs")
  print(dim(simucondorsup2simu[simucondorsup2simu$condor>2,]))
  print("Number of mutations that have <3 Simulated EEMs and >=3 Condor EEMs")
  print(dim(simucondorsup2condor[simucondorsup2condor$simu<3,]))
  simucondorsup2all$diff=simucondorsup2all$simu-simucondorsup2all$condor
  dcast(simucondorsup2all, rep~diff, length) 
}
