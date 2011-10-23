#! /bin/sh 
# ./mergeTrees.sh expects results of the selection in subdirectories such as results/H120, results/bkg1, results/bkg2...
# it creates a merged root file in the self-created results/merged/

mkdir -p results/merged

echo "Now merging Ele datasets for mass $Hmass ..."
for i in 150 500 1000 1500
do
  hadd results/merged/axiW$i\_Ele.root results/Summer11_V1/AxigluonW_M-$i\_TuneZ2_7TeV-calchep-pythia/*datasetEle.root
done
hadd results/merged/TTbar_Ele.root results/Summer11_V1/TTJets_TuneZ2_7TeV-madgraph-tauola/*datasetEle.root
hadd results/merged/SingleT_sChannel_Ele.root results/Summer11_V1/T_TuneZ2_s-channel_7TeV-powheg-tauola/*datasetEle.root
hadd results/merged/SingleTbar_sChannel_Ele.root results/Summer11_V1/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/*datasetEle.root
hadd results/merged/SingleT_tChannel_Ele.root results/Summer11_V1/T_TuneZ2_t-channel_7TeV-powheg-tauola/*datasetEle.root
hadd results/merged/SingleTbar_tChannel_Ele.root results/Summer11_V1/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/*datasetEle.root
hadd results/merged/SingleT_tWChannel_Ele.root results/Summer11_V1/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*datasetEle.root
hadd results/merged/SingleTbar_tWChannel_Ele.root results/Summer11_V1/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*datasetEle.root
hadd results/merged/Wjets_Ele.root results/Summer11_V1/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/*datasetEle.root

echo "Now merging Mu datasets for mass $Hmass ..."
for i in 150 500 1000 1500
do
  hadd results/merged/axiW$i\_Mu.root results/Summer11_V1/AxigluonW_M-$i\_TuneZ2_7TeV-calchep-pythia/*datasetMu.root
done
hadd results/merged/TTbar_Mu.root results/Summer11_V1/TTJets_TuneZ2_7TeV-madgraph-tauola/*datasetMu.root
hadd results/merged/SingleT_sChannel_Mu.root results/Summer11_V1/T_TuneZ2_s-channel_7TeV-powheg-tauola/*datasetMu.root
hadd results/merged/SingleTbar_sChannel_Mu.root results/Summer11_V1/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/*datasetMu.root
hadd results/merged/SingleT_tChannel_Mu.root results/Summer11_V1/T_TuneZ2_t-channel_7TeV-powheg-tauola/*datasetMu.root
hadd results/merged/SingleTbar_tChannel_Mu.root results/Summer11_V1/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/*datasetMu.root
hadd results/merged/SingleT_tWChannel_Mu.root results/Summer11_V1/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*datasetMu.root
hadd results/merged/SingleTbar_tWChannel_Mu.root results/Summer11_V1/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*datasetMu.root
hadd results/merged/Wjets_Mu.root results/Summer11_V1/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/*datasetMu.root

