#! /bin/sh 
# ./fromSingleJobResToRooDataSet.sh:
# it runs several scripts that bring from the results of the selection stored in the directory "results/"
# to the merged RooDataSets of the species used for the fit

lumi=$1

echo "COUNTING EVENTS PROCESSED IN MC..."
root -b -l <<EOF
.L AxigluonEventCounter.cc+
countEvents();
.q
EOF
echo "COUNTING EVENTS PROCESSED IN MC DONE."

echo "MERGING TREES STEP..."
./mergeTrees.sh
echo "MERGING TREES STEP DONE."

echo "WEIGHTING TREES STEP..."
chmod a+x weightTrees.sh
./weightTrees.sh $lumi $lumi 
echo "MERGING TREES STEP DONE."

echo "MERGING WEIGHTED TREES ACCORDING FIT SPECIES DEFINITION..."
./mergeMultiSamples.sh
echo "MERGING WEIGHTED TREES ACCORDING FIT SPECIES DEFINITION DONE."

echo "MERGING FINAL STATES..."
./mergeFinalStates.sh
echo "MERGING FINAL STATES DONE."

#echo "CREATING FIT ROODATASETS STEP..."
#chmod u+x ./createFitDatasets.sh
#./createFitDatasets.sh $Hmass
#echo "CREATING FIT ROODATASETS DONE."
