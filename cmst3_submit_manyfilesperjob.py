#! /usr/bin/env python
import os
import sys
import time
import commands
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
if len(sys.argv) != 8:
    print "usage cmst3_submit_manyfilesperjob.py process dataset nfileperjob applicationName queue prefix isMC"
    sys.exit(1)
process = sys.argv[1]
dataset = sys.argv[2]
output = dataset
# choose among cmt3 8nm 1nh 8nh 1nd 1nw 
queue = sys.argv[5]
nfileperjob = int(sys.argv[3])
application = sys.argv[4]
prefix = sys.argv[6]
isMC = int(sys.argv[7])
if isMC != 0:
    inputlist = "../HiggsAnalysisTools/cmst3_42X/MC/"+process+"/"+dataset+".list"
else:
    inputlist = "../HiggsAnalysisTools/cmst3_42X/"+process+"/"+dataset+".list"
# to write on local disks
################################################
diskoutputdir = "/cmsrm/pc23_2/crovelli/data/Axi4.2.X/"
diskoutputmain = diskoutputdir+"/"+prefix+"/"+process+"/"+output
# prepare job to write on the cmst3 cluster disks
################################################
os.system("mkdir -p "+prefix+"/"+process+"/"+output)
os.system("mkdir -p "+prefix+"/"+process+"/"+output+"/log/")
os.system("mkdir -p "+prefix+"/"+process+"/"+output+"/input/")
os.system("mkdir -p "+prefix+"/"+process+"/"+output+"/src/")

if diskoutputdir != "none": 
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm23 mkdir -p "+diskoutputmain)

#look for the current directory
#######################################
pwd = os.environ['PWD']
#######################################
inputListfile=open(inputlist)
inputfiles = inputListfile.readlines()
ijob=0

#copy the configuration in the actual run directory
os.system("cp -r config "+prefix)

while (len(inputfiles) > 0):
    inputfilename = pwd+"/"+prefix+"/"+process+"/"+output+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename,'w')
    for line in range(min(nfileperjob,len(inputfiles))):
        ntpfile = inputfiles.pop()
        if ntpfile != '':
            inputfile.write(ntpfile)
            

    inputfile.close()

    # prepare the script to run
    outputname = prefix+"/"+process+"/"+output+"/src/submit_"+str(ijob)+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export STAGE_HOST=castorcms\n')
    outputfile.write('export STAGE_SVCCLASS=cmst3\n')
    outputfile.write('cp -r '+pwd+"/"+prefix+'/config $WORKDIR\n')
    outputfile.write('cp -r '+pwd+'/data $WORKDIR\n')
    outputfile.write('export SCRAM_ARCH=slc5_amd64_gcc434\n')
    outputfile.write('cd /afs/cern.ch/user/c/crovelli/scratch0/Vecbos2010/CMSSW_4_2_6/\n')
    outputfile.write('eval `scramv1 runtime -sh`\n')
    outputfile.write('cd $WORKDIR\n')
    outputfile.write(pwd+'/'+application+' '+inputfilename+" "+output+"_"+str(ijob)+" "+str(isMC)+" "+dataset+"\n")

    outputfile.write('ls *.root | grep -v Z_calibFall08 | xargs -i scp -o BatchMode=yes -o StrictHostKeyChecking=no {} pccmsrm23:'+diskoutputmain+'/{}\n') 
    outputfile.close
    os.system("echo bsub -q "+queue+" -o "+prefix+"/"+output+"/log/"+output+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    os.system("bsub -q "+queue+" -o "+prefix+"/"+process+"/"+output+"/log/"+output+"_"+str(ijob)+".log source "+pwd+"/"+outputname+" -copyInput="+process+"_"+str(ijob))
    ijob = ijob+1
    #time.sleep(2)
    continue
