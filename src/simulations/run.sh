git clone git@github.com:mariemorel/condor.git

mkdir results_simulations
pushd results_simulations
../run_simu_12.sh
popd

mkdir results_eems
pushd results_eems
../run_count_eems.sh
popd
