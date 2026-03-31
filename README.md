Scripts to simulate trees in R (FossilSim) and alignments in IQ-TREE \
1. Go to the `fossilWithSeqs` folder \
2. `Rscript Rscripts/simulateTreesSeqs.R`. Modify `Rscripts/simulateTreesSeqs.R` to define branching and clock parameters of your choice and IQ-TREE alignment simulation settings. Simulated parameters and trees will be stored in a `truth` folder. \
2. `bash submit_all_reps.sh` will submit BEAST 2 analysis of the simulated datasets. Modify `XML/run.xml` to define BEAST 2 analysis settings. \