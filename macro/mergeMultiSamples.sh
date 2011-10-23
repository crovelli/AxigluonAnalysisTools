#! /bin/sh 
# ./mergeTrees.sh expects results of single sample files in results/merged/
# it creates a merged root file for each single specie of the fit in results/datasets

mkdir -p results/datasets_trees

echo "Now merging species for ele..."
# signal is always a species per se
for i in 150 500 1000 1500
do
  hadd results/datasets_trees/axiW$i\_Ele.root results/merged/axiW$i\_Ele.root 
done

# Wjets is a species per se
cp results/merged/Wjets_Ele.root results/datasets_trees/Wjets_Ele.root

# merging ttbar and single t in a species
hadd results/datasets_trees/top_Ele.root results/merged/TTbar_Ele.root results/merged/SingleT_tWChannel_Ele.root results/merged/SingleT_tChannel_Ele.root results/merged/SingleT_sChannel_Ele.root results/merged/SingleTbar_tWChannel_Ele.root results/merged/SingleTbar_tChannel_Ele.root results/merged/SingleTbar_sChannel_Ele.root 

# merging dibosons
hadd results/datasets_trees/dibosons_Ele.root results/merged/WW_Ele.root results/merged/WZ_Ele.root

# others (DY, gamma+jets)
hadd results/datasets_trees/others_Ele.root results/merged/GJ_40-100_Ele.root results/merged/GJ_100-200_Ele.root results/merged/DY_Ele.root 

echo "Now merging species for mu..."
# signal is always a species per se
for i in 150 500 1000 1500
do
  hadd results/datasets_trees/axiW$i\_Mu.root results/merged/axiW$i\_Mu.root 
done

# Wjets is a species per se
cp results/merged/Wjets_Mu.root results/datasets_trees/Wjets_Mu.root

# merging ttbar and single t in a species
hadd results/datasets_trees/top_Mu.root results/merged/TTbar_Mu.root results/merged/SingleT_tWChannel_Mu.root results/merged/SingleT_tChannel_Mu.root results/merged/SingleT_sChannel_Mu.root results/merged/SingleTbar_tWChannel_Mu.root results/merged/SingleTbar_tChannel_Mu.root results/merged/SingleTbar_sChannel_Mu.root 

# merging dibosons
hadd results/datasets_trees/dibosons_Mu.root results/merged/WW_Mu.root results/merged/WZ_Mu.root

# others (DY, gamma+jets)
hadd results/datasets_trees/others_Mu.root results/merged/GJ_40-100_Mu.root results/merged/GJ_100-200_Mu.root results/merged/DY_Mu.root 
