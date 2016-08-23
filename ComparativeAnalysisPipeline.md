#**Positive Selection Pipeline for Comparative Analysis of *P. eremicus* transcriptome relative to *P. maniculatus* and *M. musculus* transcriptomes**


___________________________________________________________

***DATA SOURCES***
-
*Two eremicus files:* 

cds (nucl) and pep (prot)
/mnt/data3/lauren/NEWtranscriptome/dammit/good.reproductive_RERUNhighexpHALF.trinity.Trinity.fasta.dammit

62M Mar  8 18:10 
good.reproductive_RERUNhighexpHALF.trinity.Trinity.fasta.transdecoder.cds

27M Mar  8 18:09 
good.reproductive_RERUNhighexpHALF.trinity.Trinity.fasta.transdecoder.pep

These have been copied and placed in the new directory PAML4.8 (within two directories respectively as below)

Files were re-named: 

good.reproductive_RERUNhighexpHALF.trinity.Trinity.fasta.transdecoder.cds
/mnt/data3/lauren/PAML4.8/original_files/pero_erem.cds

good.reproductive_RERUNhighexpHALF.trinity.Trinity.fasta.transdecoder.pep
/mnt/data3/lauren/PAML4.8/original_files/pero_erem.pep

-
*Two mus files:* 

cds and pep

ftp://ftp.ensembl.org/pub/release-84/fasta/mus_musculus

Index of /pub/release-84/fasta/mus_musculus/cds/  	15.9 MB	2/25/16, 12:22:00 AM
Mus_musculus.GRCm38.cds.all.fa.gz

We renamed this file to the following:
/mnt/data3/lauren/PAML4.8/original_files/muss_musc.cds

Index of /pub/release-84/fasta/mus_musculus/pep/	10.5 MB	3/11/16, 11:34:00 AM
Mus_musculus.GRCm38.pep.all.fa.gz

We renamed this file to the following:
/mnt/data3/lauren/PAML4.8/original_files/muss_musc.pep

-
*Two maniculatus files:*

cds and pep

Peromyscus maniculatus (North American deer mouse)
Representative genome: Peromyscus maniculatus bairdii (assembly Pman_1.0)
http://www.ncbi.nlm.nih.gov/genome/?term=txid10042[orgn]

We confirmed that the "rna" file is a cds file and we downloaded it to the cds folder:

cds file was renamed from GCF_000500345.1_Pman_1.0_rna.fna: 
/mnt/data3/lauren/PAML4.8/original_files/pero_mani.cds

We confirmed that the "protein" file is a pep file and we downloaded it to the pep folder:

pep file was renamed from GCF_000500345.1_Pman_1.0_protein.faa :
/mnt/data3/lauren/PAML4.8/original_files/pero_mani.pep

__________________________________________________________

***RENAMING HEADERS***

*script renames all of the headers in a particular file to be more useable.*

The two renaming scripts were used and then placed in the paml_scripts folder.  
They are called:
/mnt/data3/lauren/PAML4.8/paml_scripts/renaming_pep_files.py 
/mnt/data3/lauren/PAML4.8/paml_scripts/renaming_cds_files.py

We ran the pep one on all three pep files, so that all three have good headers. The three re-named 
pep files were then afterwards relabeled as fasta files, and they are named:
/mnt/data3/lauren/PAML4.8/pep_files/peps/pero_mani_renamed.fasta
/mnt/data3/lauren/PAML4.8/pep_files/peps/muss_musc_renamed.fasta
/mnt/data3/lauren/PAML4.8/pep_files/peps/pero_erem_renamed.fasta

We ran cds one on all three cds files, so that all three have good headers. The three re-named 
cds files are named:
/mnt/data3/lauren/PAML4.8/cds_files/pero_mani_renamed.cds
/mnt/data3/lauren/PAML4.8/cds_files/muss_musc_renamed.cds
/mnt/data3/lauren/PAML4.8/cds_files/pero_erem_renamed.cds


***These are the two scripts***

**renaming_pep_files.py**

!/usr/bin/python3
A program for renaming fasta headers
USAGE: ./renaming_pep_headers_TA.py --input path_to_input_directory
Author: Taruna Aggarwal
Affiliation: University of New Hampshire, Durham, NH, USA
Date: 1/27/2016
Purpose is to replace headers with their corresponding file names
and add consecutive numbers to the new headers
This script assumes that your pep files are named in a specific manner:
SpeciesInitial_genusName.pep and an example is F_solani.pep
The script will generate new files in the same directory as itself.

-

import argparse
import os

parser = argparse.ArgumentParser(description="This script renames files and their headers in a directory.")
parser.add_argument('--input', help="PATH to the directory with input files.", required=True)
args = parser.parse_args()

for file in os.listdir(args.input):
    if file.endswith(".pep"):
        working_file = open(args.input + '/' + file, "r")
        new_file = open(file[:-4] + "_renamed.pep", "w")
        print("Renaming {0}".format(file))
        counter = 1
        for currentLine in working_file:
            currentLine = currentLine.rstrip()
            if currentLine.startswith(">"):
                new_file.write(">{0}_{1}\n".format((file[:-4]), counter))
                counter += 1
            else:
                new_file.write("{0}\n".format(currentLine))
 
-
               


**renaming_cds_files.py**

!/usr/bin/python3
A program for renaming fasta headers
USAGE: ./renaming_cds_headers_TA.py --input path_to_input_directory
Author: Taruna Aggarwal
Affiliation: University of New Hampshire, Durham, NH, USA
Date: 1/27/2016
Purpose is to replace headers with their corresponding file names
and add consecutive numbers to the new headers
This script assumes that your cds files are named in a specific manner:
SpeciesInitial_genusName.cds and an example is F_solani.cds
The script will generate new files in the same directory as itself.

-
import argparse
import os

parser = argparse.ArgumentParser(description="This script renames files and their headers in a directory.")
parser.add_argument('--input', help="PATH to the directory with input files.", required=True)
args = parser.parse_args()

for file in os.listdir(args.input):
    if file.endswith(".cds"):
        working_file = open(args.input + '/' + file, "r")
        new_file = open(file[:-4] + "_renamed.cds", "w")
        print("Renaming {0}".format(file))
        counter = 1
        for currentLine in working_file:
            currentLine = currentLine.rstrip()
            if currentLine.startswith(">"):
                new_file.write(">{0}_{1}\n".format((file[:-4]), counter))
                counter += 1
            else:
                new_file.write("{0}\n".format(currentLine))

-
________________________________________

***ORTHOFINDER: v0.6.1***

*This software finds orthologous groups of protein sequences using two scripts*

1) The first orthofinder command finds all of the orthologs (but does not trim down to SCOs):

You can retrieve the orthofinder.py script on the github website.


The command I used to run this script is:

python /mnt/data3/lauren/bin/OrthoFinder/orthofinder.py -f peps/ -t 16


*David Emms wrote orthofinder.py*

To site the program, use the citation below on the website: Emms & Kelly 2015

When publishing work that uses OrthoFinder please cite:

D.M. Emms & S. Kelly (2015), OrthoFinder: solving fundamental biases in whole genome comparisons
dramatically improves orthogroup inference accuracy, Genome Biology 16:157.
   
-

2) The second orthofinder command was slightly modified to make trees with FastTreeMP, because the
script on the website uses FastTree (a different version of FastTree)
That was the only change to the script.

You can retrieve the trees_for_orthogroups.py script on the github website.

The command I used to run this script is:

python /mnt/data3/lauren/bin/OrthoFinder/trees_for_orthogroups.py Results_Jun28_5 -t 16

To site the program, use the citation on the website: Emms & Kelly 2015

When publishing work that uses OrthoFinder please cite:

D.M. Emms & S. Kelly (2015), OrthoFinder: solving fundamental biases in whole genome comparisons
dramatically improves orthogroup inference accuracy, Genome Biology 16:157.


Fasta files for orthogroups have been written to:
   Results_Jun28_5/Sequences/

Multiple sequences alignments have been written to:
   Results_Jun28_5/Alignments/

Gene trees have been written to:
   Results_Jun28_5/Trees/

______________________________________

***Selecting Single Copy Orthologs (SCOs) from the orthologous groups generated by OthoFinder***


We parsed orthogroups from the alignments directory because we only want the SCOs:


This is the command to run the script:

bash ./sco.sh


***This is the script:***

**sco.sh**

!/bin/bash

To run, you must be in the `Results` folder produced by OrthoFinder

bash ./sco.sh


This script takes the output from OrthoFinder, contained in the OrthologousGroups.txt file, and pulls out the single copy orthologes
The only thing you have to change is the number of species in your OrthoFinder run, so, change NSPECIES from 2 to whatever umber is appropriate

The only little issue is that you have to have your species names in some cogent way, so that each species is uniquely recognizable
My protein names are like `Haliaeetus_albicilla_N`, where N is an integer.
To disambiguate `Haliaeetus_albicilla` from `Haliaeetus_leucocephalus` I need 12 characters, which is why `awk '{print substr($0,0,12)}'` is 12 and not something else.
look at `Results*/WorkingDirectory/SpeciesIDs.txt` to figure this out.


-
NSPECIES=3  
TRUNC=6

min=$(expr $NSPECIES / 2 + 2)
equal=$(expr $NSPECIES + 1)
input=OrthologousGroups.txt
END=$(wc -l $input | awk '{print $1}')
START=1
LIMIT=$(expr $NSPECIES + 2)

rm SCOs.txt 2> /dev/null

for i in $(eval echo "{$START..$END}") ; do
    sed -n ''$i'p' $input | awk '!$var' var="$LIMIT" | awk '!array[$0]++'  | tr -s ' ' \\n | awk '{print substr($0,0,tr)}' tr="$TRUNC" | sort -u | wc -l > 1.txt;
    sed -n ''$i'p' $input | awk '!$var' var="$LIMIT" | awk '!array[$0]++'  | tr -s ' ' \\n | awk '{print substr($0,0,tr)}' tr="$TRUNC" | wc -l > 2.txt;
    if [ $(cat 1.txt) -eq $(cat 2.txt) ] && [ $(cat 1.txt) -gt "$min" ]; then sed -n ''$i'p' $input >> SCOs.txt; fi ;
done

rm 1.txt 2.txt 

-


The SCO output file is here: 
/mnt/data3/lauren/PAML4.8/pep_files/peps/Results_Jun28_5/SCOs.txt

This txt file contains only the Orthogroup IDs which are SCOs and the corresponding transcript ID matches for each species for each of these SCOs.

This is an example line from this .txt file:

OG0010552: muss_musc_1368 pero_erem_53235 pero_mani_8257
_________________________________________________________________________________

***Selecting cds files for SCOs***

We selected the cds file transcripts for all three species corresponding to the SCOs

Recall that the renamed cds files are named:

/mnt/data3/lauren/PAML4.8/cds_files/pero_mani_renamed.cds

/mnt/data3/lauren/PAML4.8/cds_files/muss_musc_renamed.cds

/mnt/data3/lauren/PAML4.8/cds_files/pero_erem_renamed.cds


I copied all of those cds files into a new directory called cds_cat

I made a combined.fasta file with these concatenated cds files into cds_cat using the following command:

cat muss_musc_renamed.cds pero_erem_renamed.cds pero_mani_renamed.cds > combined.fasta

The concatenated cds file is located here:
/mnt/data3/lauren/PAML4.8/cds_files/cds_cat/combined.fasta


Next UNWRAP the sequence lines in the combined.fasta file that is in the cds_cat directory.
This was necessary because the script that I will execute needs unwrapped lines

The unwrapping command is as follows:
sed -i ':begin;$!N;/[ACTGNn-]\n[ACTGNn-]/s/\n//;tbegin;P;D' combined.fasta

I renamed this file as UNWRAPPEDcombined.fasta
It is located here:
/mnt/data3/lauren/PAML4.8/cds_files/cds_cat/UNWRAPPEDcombined.fasta


Next I ran a script which selects all of the cds transcripts from the UNWRAPPEDcombined.fasta file according to the transcript IDs for each SCO line from the SCOs.txt file

The command I used to run this script is:

bash ./sco2cds.sh


***This is the script:*** 


**sco2cds.sh**


!/bin/bash

To run, you must be in the `Results` folder produced by OrthoFinder

  bash ./sco2cds.sh


Give the name of the input file - this should be produced by sco.sh and filtered however you want to, e.g., to contain only certaint OrthoGroups


-
input=SCOs.txt
START=1
END=$(wc -l $input | awk '{print $1}')
mkdir OGs
for i in $(eval echo "{$START..$END}") ; do
    for j in $(sed -n ''$i'p' $input | awk '!array[$0]++'  | tr -s ' ' \\n) ; do
        grep --no-group-separator --max-count=1 -w -hA1 $j /mnt/data3/lauren/PAML4.8/cds_files/cds_cat/UNWRAPPEDcombined.fasta >> OGs/$(sed -n ''$i'p' $input | awk -F : '{print $1}').fasta;
    done
done

-

This script generates an orthogroups directory that only contains the SCO fasta files.
Specifically, each fasta file in this directory has the header and sequence corresponding to each of the three species for the SCO orthogroup.

This is one example of a fasta file in the directory:  OG0012588.fasta

>muss_musc_26677
ATGGCCCAGGCGAAGATCAGCGCCAAGGCCCACGAGGGCCGCTTCTGCCGCTCGTCTTCCATGGCCGACCGCTCCAGCCGCCTGCTGGAGAGTCTGGACCAGCTGGAGC
TCAGGGTGGAAGCTTTGAGAGACGCAGCTACTGCTGTTGAGCAGGAGAAAGAAATCCTTCTGGAGATGATCCACAGTATCCAAAACAGCCAGGACATGAGGCAGATTAG
CGATGGAGAAAGAGAGGAACTAAACCTGACTGCCAACCGTCTGATGGGCCGGACCCTCACGGTTGAGGTCTCGGTGGAAACGATCCGAAACCCCCAGCAGGAGGAATCC
TTGAAGCATGCCACGAGGATTATAGACGAGGTGGTCAGCAAGTTCCTAGATGACCTGGGGAATGCCAAGAGCCACTTAATGTCACTTTACAGTGCCTGCTCATCGGAGG
TGCCGCCTGGGCCGGTGGACCAGAAGTTTCAATCGATAGTCATCGGTTGCGCTCTTGAGGATCAGAAGAAAATCAAGAGGCGATTGGAGACTCTGCTGAGGAACATTGA
CAACTCCGACAAGGCCATTAAACTCCTAGAGCATGCTAAAGGAGCTGGTTCCAAAAGCCTGCAGAACACTGACGGCAAATTTAATTAG
>pero_erem_59032
ATGGCCCAGGCGAAGATCAGCGCCAAGGCCCACGAGGGCCGCTTCTGCCGCTCCTCGTCCATGGCCGACCGCTCCAGCCGCCTGCTGGAGAGCCTGGACCAGCTGGAGC
TCAGGGTTGAAGCTTTGAGAGACGCAGCTACTGCTGTTGAACAAGAGAAAGAAATCCTTCTGGAGATGATCCACAGCATCCAGAACAGCCAGGACATGAGGCAGATTAG
TGACGGAGAAAGAGAGGAATTAAACCTCACTGCAAACCGCCTCATGGGCCGAACTCTCACTGTCGAGGTCTCAGTGGAAACAATTCGAAGTCCCCAGCAGGAGGAATCC
TTGAAGCACGCCACGAAGATTATAGATGAGGTGGTCACTAAGTTCCTGGACGACCTGGGAAATGCCAAGAGTCACTTAATGTCACTCTACAGTGCGTGTTCATCCGAGG
TGCCCCCCGGGCCAGTTGACCAAAAATTTCAGTCGATAGTAATCGGTTGCGCTCTTGAAGATCAGAAGAAAATTAAGAGGCGACTAGAGACTCTGCTGAGGAACATCGA
CAACTCCGACAAGGCCATTAAACTATTAGAGCACTCTAAAGGAGGCGGTCCCAAAAGTCTGCAAAACACTGACGGCAAATTCAATTAG
>pero_mani_22364
ATGTCAGAAAAGCACTTAAAACTTCGAGAGGCTAAGCCACAAGGCAAAGCCACGATAAGAACCACATCAGAGGTTGAAGCTTTGAGAGATGCAGCTACTGCTGTTGAAC
AAGAGAAAGAAATCCTTCTGGAGATGATCCACAGCATCCAGAACAGCCAGGACATGAGGCAGATTAGTGACGGAGAAAGAGAAGAATTAAACCTCACTGCAAACCGCCT
GATGGGCCGCACGCTCACTGTTGAGGTCTCAGTGGAGACAATTCGAAGTCCCCAGCAGGAGGAATCCTTGAAGCACGCCACGAAGATTATAGATGAGGTGGTCAGTAAG
TTCCTGGACGACCTGGGAAATGCCAAGAGCCACTTAATGTCACTCTACAGTGCGTGTTCATCTGAGGTGCCGCCCGGGCCAGTTGACCAAAAATTTCAGTCGATAGTAA
TCGGTTGCGCTCTCGAAGATCAGAAGAAAATTAAGAGGCGACTAGAGACTTTGCTGAGGAACATTGACAACTCCGACAAGGCCATTAAACTATTAGAGCACTCTAAAGG
AGGCGGTCCCAAGAGTCTGCCAAACATTGACGGCAAGTTCAATTAGTCTTCCAACCCACAGGCCCTTACAAACGATGTAACAGAGAAAAACCACTATTTTAAATAACTA
GTTCTTTATGTTAGGCATAACCACTTATACCTTATACTGATAACTGTTTCCGATGAGGAAGTATTCCAGTGTGGATAAGGTACGCAAATCACACTTAAACCTAGAAGTA
TTTCAGTCGCCTGCTCACAGGTTTTGGACAGCTTTGCCGTCTTGAACTGGACTCTAGCCCTGGATATTCAACTTCTGAAGATGCACACATGGTCTGTGGGAAAATAAAT
AAGGGGCAAATGCATCACACTGCTCTTTCCTTAACCTCAGTACAGATGCACCCAGCCCCCGGAGCAGTGCACTCACATCACATTGTGAACAG

_____________________________________________

Now that we have all the cds transcripts for the SCOs, we want to run our comparative analysis for positive selection

However, this analysis requires tree files, so we need to modify the tree files generated several steps above by the second OrthoFinder script (trees_for_orthogroups.py), which are in this directory:
Gene trees have been written to:
   Results_Jun28_5/Trees/

However, we must modify the tree files so that the syntax of the species naming is compatible with our next command script. 


This command modifies the names in the trees files:

for i in $(ls *txt); do sed -i 's/muss_musc_renamed_muss_musc_/muss_musc_/;s/\pero_erem_renamed_pero_erem_/pero_erem_/;s/\pero_mani_renamed_pero_mani_/pero_mani_/' $i; done


This directory contains Trees for all possible orthogroups for our three species. We only be using the trees from the SCOs.  

However, the script we will use will only select the appropriate trees (those with matches in the OGs file). So there is actually no need to pre-select the SCOs tree files to put in the OGs directory.

___________________________________________

***PRANK v150803 alignment and analysis for positive selection with codeml in PAML 4.8 codeml***

The M2 Branch-site test for positive selection requires the user to run two different models:

M1a (the null model - where omega is set at 1)

M2a (the alternative model - for positive selection)


Therefore I made two new directories to run the models separately:

/mnt/data3/lauren/PAML4.8/pep_files/peps/Results_Jun28_5/OGs/M1a

/mnt/data3/lauren/PAML4.8/pep_files/peps/Results_Jun28_5/OGs/M2a


In order to do this analysis, I copied all the modified tree files and all of SCO fasta files into both of these directories.

However, in order to run this branch-site test, I needed to place *P. eremicus* in the foreground, therefore I had to add #1 to this species for every tree file in both model directories.

The command for executing this action is below:

for i in $(ls *tree.txt); do sed -i 's_,_\n,_g' $i && sed -i '/erem/s/$/ #1/' $i; done


Within both model directories (M1a and M2a) all of the tree files were edited appropriately.


*I also added three scripts to both of these directories:*

codeml.ctl

autoPAML.py

autoPAML.sh


The only file that is different for the two directories is the codeml.ctl file.  
This control file has different settings specific to the M1a and M2a models


***The codeml.ctl file for M1a is below:***

**control.ctl**

-
noisy = 1
verbose = 0
runmode = 0
seqtype = 1
CodonFreq = 2
model = 2
NSsites = 2
fix_omega = 1
omega = 1

-

***The codeml.ctl file for M2a is below:***

**control.ctl**

-
noisy = 1
verbose = 0
runmode = 0
seqtype = 1
CodonFreq = 2
model = 2
NSsites = 2
fix_omega = 0

-


The autoPAML.sh file is the file that you implement for running the PRANK aligner and the codeml analysis in the command line. In contrast, the python file and the control file are utilized within this shell script 

***This is the command I used to run the script***

bash ./autoPAML.sh -t 25

***This is the script for autoPAML.sh***

**autoPAML.sh**

! /bin/bash

usage=$(cat << EOF
    This script runs a pipeline that takes a fasta file and BAMfiles and tests for selection:

-
 
   autoPAML.sh [options]
   Options:
      -t <v> : *required* Numberof threads to use.
EOF
);

while getopts o:t: option
do
        case "${option}"
        in
		t) TC=${OPTARG};;
        esac
done


for inputaln in $(ls *fasta); do
    F=$(basename "$inputaln" .fasta)
    if [ $(ps -all | grep 'prank\|codeml' | wc -l | awk '{print $1}') -lt $TC ] ;
    then
        if [ ! -f $F.out ] ;
        then
            echo 'processing' $inputaln
            awk '{print $1}' $inputaln > tmp && mv tmp $inputaln
            prank -d="$inputaln" -translate -F -o=$F &&
            perl /share/pal2nal.v14/pal2nal.pl $F.best.pep.fas $F.best.nuc.fas -output fasta -nogap -nomismatch > $F.clean || true &&
            python autoPAML.py $F.clean $F_tree.txt $F.out &&
            python autoPAMLresults.py $F.out | tee -a paml.results &
        else
            echo 'next'
        fi
    else
        until [ $(ps -all | grep 'prank\|codeml' | wc -l | awk '{print $1}') -lt $TC ] ;
        do
            echo 'waiting for' $inputaln
            sleep 25s;
        done
        if [ ! -f $F.out ] ;
        then
            awk '{print $1}' $inputaln > tmp && mv tmp $inputaln
            prank -d="$inputaln" -translate -F -o=$F &&
            perl /share/pal2nal.v14/pal2nal.pl $F.best.pep.fas $F.best.nuc.fas -output fasta -nogap -nomismatch > $F.clean || true &&
            python autoPAML.py $F.clean $F_tree.txt $F.out &&
            python autoPAMLresults.py $F.out | tee -a paml.results &
        else
            echo 'moving on'
        fi
    fi
done

-



***Below is the autoPAML.py script utilized by the above shell script***


**autoPAML.py**

!/usr/bin/python

usage  python autoPAML.py alignment.fasta tree.nexus paml.output

-

from __future__ import division
from Bio.Phylo.PAML import codeml
import sys

if len(sys.argv) != 4:
        print "Error.  our missing an input file, see the usage"
        quit()

cml = codeml.Codeml()
cml.read_ctl_file("codeml.ctl")
cml.alignment = sys.argv[1]
cml.tree = sys.argv[2]
cml.out_file = sys.argv[3]
cml.run()

-


Once the M1a and M2a analyses were performed, each directory contained .out files to correspond to the tests completed for each SCO.

Within the M1a directory all output files were renamed .out -> .null.out 

rename -v 's/.out/.null.out/' *

Within the M2a directory all output files renamed .out -> .alt.out 

rename -v 's/.out/.alt.out/' *

-

NEXT I compared the results for the M2a and M1a runs of PAML:

So it is a matter of selecting  the appropriate numbers into a file from the corresponding .null.out and .alt.out files for each SCO, and performing a likelihood ratio test on these values to determine if there is evidence of positive selection for any of the *P. eremicus* genes within these SCOs. 

This was implemented with a python script that parsed the results for the output files and another python script which was implemented within this executed script.

The names of these two scripts are below: 

runPAML_results_TA.py

results_paml.py



I ran the following command to execute the runPAML_results_TA.py script:

python3 runPAML_results_TA.py --null /mnt/data3/lauren/PAML4.8/pep_files/peps/Results_Jun28_5/OGs/M1a --alt /mnt/data3/lauren/PAML4.8/pep_files/peps/Results_Jun28_5/OGs/M2a 


***The runPAML_results_TA.py script is below:***


**runPAML_results_TA.py**

!/usr/bin/python3
A program for calculating adjusted P-value for PAML's codemls outputs.
USAGE: ./runPAML_results_TA.py
Author: Taruna Aggarwal
Affiliation: University of New Hampshire, Durham, NH, USA
Date: 01/27/2016

-

import sys
import os
import subprocess
import argparse

parser = argparse.ArgumentParser(description="This script calculates p-value for PAML's output.")
parser.add_argument('--null', help="PATH to null files", required=True)
parser.add_argument('--alt', help="PATH to alt files", required=True)
args = parser.parse_args()


def analyzePAML(nullfile, altfile):
    analyzePAML_cmd = "python results_paml.py {0} {1} | tee -a PAML.results".format(nullfile, altfile)
    subprocess.call(analyzePAML_cmd, shell=True)

for file in os.listdir(args.null):
  if file.endswith(".null.out"):
    output = analyzePAML(args.null + "/" + file, args.alt + "/" + file[:-9] + ".alt.out")

-

***The results_paml.py script is below:***

**results_paml.py**

!/usr/bin/python

-

from __future__ import division
from Bio.Phylo.PAML import codeml
import sys
import os
import shutil
import time
import glob
from math import sqrt
from rpy2 import robjects

r = robjects.r
len = len(glob.glob1(".","*.out"))

def compare_models(null_lnl, alt_lnl, df):
    likelihood = 2*(abs(null_lnl-alt_lnl))
    p = 1 - robjects.r.pchisq(likelihood, df)[0]
    return p


null_results = codeml.read(sys.argv[1])
alt_results = codeml.read(sys.argv[2])

null_nssites = null_results.get("NSsites")
alt_nssites = alt_results.get("NSsites")

null_value = null_nssites.get(2)
null_lnl = null_value.get("lnL")

alt_value = alt_nssites.get(2)
alt_lnl = alt_value.get("lnL")

bs_p_pos = compare_models(null_lnl,alt_lnl,1)

r.assign('bs_p_pos', bs_p_pos)
r.assign('len', len)

paBS = robjects.r('p.adjust(bs_p_pos, "BH", len)')

print 'Ajusted p-val for null vs. alt of orthogroup {} is {}'.format(sys.argv[1], paBS[0])

-


The output file generated by running the runPAML_results_TA.py script is in:
/mnt/data3/lauren/PAML4.8/pep_files/peps/Results_Jun28_5/OGs/PAML.results


The likelihood ratio test (LRT) was run successfully for 2698 paired null and alt .out files for my analysis

I also pulled out the significant p-values (p < 0.05) for the LRT from the PAML.resuls output file with this command:

cat PAML.results | awk '$11<0.05' > significant_orthogroups.txt


There are 42 Significant results in this significant_orthogroups.txt file reflecting evidence for 42 genes under positive selection in the *P. eremicus* lineage according to the likelihood ratio test comparisons between the null and alternative models


-

For each of these 42 genes under positive selection in *P. eremicus*, I selected the *M. musculus* cds sequence from the corresponding SCO fasta file (in the OG directory).
 
AND then I performed BLASTn on these sequences to find matches on NCBI in order to determine the gene IDs for the *P. eremicus* genes under positive selection within the SCO groups tested in this analysis.

The results are in a text wrangler file called: 42SeqsM2aPosSel.txt

