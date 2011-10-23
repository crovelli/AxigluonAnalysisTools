#!/usr/bin/perl

use Time::Local;
use Getopt::Std;

getopts('p:');
if($opt_p) {$prefix = $opt_p;}
else { die "usage: ./submitMassDependentAnalysis.pl -p <prefix>";}

@masses = (150,500,1000,1500);

for($i=0; $i<($#masses+1); $i++) {
    $mass = $masses[$i];
    print "-------------------------->\n";
    print "SUBMITTING MASS SELECTION: mAxi = $mass ...\n\n";
    open(MASSFILE,">config/axi/axiMass.txt");
    print MASSFILE "axiMass\t$mass\n";

    print "submitting signals...\n";
    $axiList = "AxigluonW_M-".$mass."_TuneZ2_7TeV-calchep-pythia";
    system("python cmst3_submit_manyfilesperjob.py Summer11_V1 $axiList 15 AxigluonApp 8nh $prefix 1");
}   

print "done with signals. \n";

open(MASSFILE,">config/axi/axiMass.txt");
print MASSFILE "axiMass\t 150\n";

print  "submitting top...\n";
system("python cmst3_submit_manyfilesperjob.py Summer11_V1 T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola 15 AxigluonApp 8nh $prefix 1");
system("python cmst3_submit_manyfilesperjob.py Summer11_V1 T_TuneZ2_tW-channel-DS_7TeV-powheg-tauola 15 AxigluonApp 8nh $prefix 1");
system("python cmst3_submit_manyfilesperjob.py Summer11_V1 Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola 15 AxigluonApp 8nh $prefix 1");
system("python cmst3_submit_manyfilesperjob.py Summer11_V1 Tbar_TuneZ2_tW-channel-DS_7TeV-powheg-tauola 15 AxigluonApp 8nh $prefix 1");
system("python cmst3_submit_manyfilesperjob.py Summer11_V1 T_TuneZ2_t-channel_7TeV-powheg-tauola 15 AxigluonApp 8nh $prefix 1");
system("python cmst3_submit_manyfilesperjob.py Summer11_V1 Tbar_TuneZ2_t-channel_7TeV-powheg-tauola 15 AxigluonApp 8nh $prefix 1");
system("python cmst3_submit_manyfilesperjob.py Summer11_V1 T_TuneZ2_s-channel_7TeV-powheg-tauola 15 AxigluonApp 8nh $prefix 1");
system("python cmst3_submit_manyfilesperjob.py Summer11_V1 Tbar_TuneZ2_s-channel_7TeV-powheg-tauola 15 AxigluonApp 8nh $prefix 1");
system("python cmst3_submit_manyfilesperjob.py Summer11_V1 TTJets_TuneZ2_7TeV-madgraph-tauola 15 AxigluonApp 8nh $prefix 1");
print  "done with top.\n";

print  "submitting V+jets...\n";
system("python cmst3_submit_manyfilesperjob.py Summer11_V1 WJetsToLNu_TuneZ2_7TeV-madgraph-tauola 15  AxigluonApp 8nh $prefix 1");

sleep 800;

print "submitting dibosons...\n";
system("python cmst3_submit_manyfilesperjob.py Summer11_V1 WW_TuneZ2_7TeV_pythia6_tauola 15  AxigluonApp 8nh $prefix 1");
system("python cmst3_submit_manyfilesperjob.py Summer11_V1 WZ_TuneZ2_7TeV_pythia6_tauola 15  AxigluonApp 8nh $prefix 1");

print "submitting gamma+jets...\n";
system("python cmst3_submit_manyfilesperjob.py Summer11_V1 GJets_TuneZ2_100_HT_200_7TeV-madgraph 15  AxigluonApp 8nh $prefix 1");
system("python cmst3_submit_manyfilesperjob.py Summer11_V1 GJets_TuneZ2_40_HT_100_7TeV-madgraph 15  AxigluonApp 8nh $prefix 1");

print "submitting DY...\n";
system("python cmst3_submit_manyfilesperjob.py Summer11_V1 DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1 15  AxigluonApp 8nh $prefix 1"); 

print "\nDONE WITH MASS $mass GeV\n";
print "<--------------------------\n";
    


