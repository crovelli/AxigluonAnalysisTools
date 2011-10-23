#! /bin/sh


lumiEE=$1
lumiMM=$2
echo "Adding weights for ele datasets for " $lumiEE " pb-1..."
root -l -b <<EOF
.L addWeightsToTree.cc+
addWeights("results/merged/axiW150_Ele.root", 0.00155066*$lumiEE, 1150 ,1);
addWeights("results/merged/axiW500_Ele.root", 8.05829e-05*$lumiEE, 1500 ,1);
addWeights("results/merged/axiW1000_Ele.root", 6.42645e-06*$lumiEE, 2000 ,1);
addWeights("results/merged/axiW1500_Ele.root", 8.38308e-07*$lumiEE, 2500 ,1);
addWeights("results/merged/Wjets_Ele.root", 0.000413442*$lumiEE, 80 ,1);
addWeights("results/merged/SingleT_sChannel_Ele.root", 9.00485e-06*$lumiEE, 15 ,1);
addWeights("results/merged/SingleTbar_sChannel_Ele.root", 2.04042e-05*$lumiEE, 16 ,1);
addWeights("results/merged/SingleT_tChannel_Ele.root", 1.03531e-06*$lumiEE, 13 ,1);
addWeights("results/merged/SingleTbar_tChannel_Ele.root", 9.9901e-07*$lumiEE, 14 ,1);
addWeights("results/merged/SingleT_tWChannel_Ele.root", 1.15222e-05*$lumiEE, 19 ,1);
addWeights("results/merged/SingleTbar_tWChannel_Ele.root", 1.12429e-05*$lumiEE, 20 ,1);
addWeights("results/merged/TTbar_Ele.root", 4.25452e-05*$lumiEE, 10 ,1);
.q

EOF

echo "Adding weights for muon datasets for " $lumiMM " pb-1..."
root -l -b <<EOF
.L addWeightsToTree.cc+
addWeights("results/merged/axiW150_Mu.root", 0.00155066*$lumiMM, 1150 ,0);
addWeights("results/merged/axiW500_Mu.root", 8.05829e-05*$lumiMM, 1500 ,0);
addWeights("results/merged/axiW1000_Mu.root", 6.42645e-06*$lumiMM, 2000 ,0);
addWeights("results/merged/axiW1500_Mu.root", 8.38308e-07*$lumiMM, 2500 ,0);
addWeights("results/merged/Wjets_Mu.root", 0.000413442*$lumiMM, 80 ,0);
addWeights("results/merged/SingleT_sChannel_Mu.root", 9.00485e-06*$lumiMM, 15 ,0);
addWeights("results/merged/SingleTbar_sChannel_Mu.root", 2.04042e-05*$lumiMM, 16 ,0);
addWeights("results/merged/SingleT_tChannel_Mu.root", 1.03531e-06*$lumiMM, 13 ,0);
addWeights("results/merged/SingleTbar_tChannel_Mu.root", 9.9901e-07*$lumiMM, 14 ,0);
addWeights("results/merged/SingleT_tWChannel_Mu.root", 1.15222e-05*$lumiMM, 19 ,0);
addWeights("results/merged/SingleTbar_tWChannel_Mu.root", 1.12429e-05*$lumiMM, 20 ,0);
addWeights("results/merged/TTbar_Mu.root", 4.25452e-05*$lumiMM, 10 ,0);
.q

EOF

echo "done weighting."
