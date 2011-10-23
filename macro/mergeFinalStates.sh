#! /bin/sh
# ./mergeFinalStates.sh merges the trees for ee/emu/mm final states in one for a combined fit
# the variable "finalstate in the tree will distinguish the three"



# usage: ./mergeFinalStates.sh Hmass
for i in 150 500 1000 1500
do
  hadd -f results/datasets_trees/axiW$i\_ll.root results/datasets_trees/axiW$i\_Ele.root results/datasets_trees/axiW$i\_Mu.root 
done
hadd -f results/datasets_trees/top_ll.root results/datasets_trees/top_Ele.root results/datasets_trees/top_Mu.root 
hadd -f results/datasets_trees/Wjets_ll.root results/datasets_trees/Wjets_Ele.root results/datasets_trees/Wjets_Mu.root 

