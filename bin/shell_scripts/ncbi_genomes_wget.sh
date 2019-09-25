#!/bin/bash

#Motivation: download Genbank genomes from NCBI FTP in parallel.
#md5 hashes of downloaded genomes will be compared to those on NCBI FTP.

#Author: Ian Rambo
#Contact: ian.rambo@utexas.edu, imrambo@lbl.gov
#Thirteen... that's a mighty unlucky number... for somebody!

#Dependencies: wget

#Download FTP tables from: https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/
#==============================================================================
#Positional Parameters
genomeTable=$1 #file from NCBI containing FTP path(s)
dbDir=$2 #output directory
nJob=$3 #number of simultaneous parallel jobs
joblog=$4 #joblog from gnu parallel
sleeptime=${5:-0s} #sleep time between jobs

#Script is currently set for Genbank genomes
#==============================================================================
###---MAIN---###
#==============================================================================
#Create the output directory if not present
[[ -d $dbDir ]] || mkdir $dbDir

printf "\n################ Downloading genome files ###############\n"
cd $dbDir
grep -o "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/[0-9]\{3\}/[0-9]\{3\}/[0-9]\{3\}/GCA_[0-9]\{9\}.[0-9]*[\.\_A-Za-z0-9\-]*" $genomeTable | \
    sed  -r 's|(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/)(GCA/)([0-9]{3}/)([0-9]{3}/)([0-9]{3}/)(GCA_.+)|\1\2\3\4\5\6/\6_genomic.fna.gz|' | \
    sort | \
    uniq | \
    parallel --jobs $nJob --joblog $joblog "wget --quiet {} && sleep $sleeptime"

#Did wget jobs complete successfully?
exitStats=$( tail -n +2 $joblog | awk '$7 != 0' | wc -l )
if [ "$exitStats" -gt 0 ]
then
    #If any wget jobs have errors, re-attempt the download(s)
    echo "WARNING: $exitStats wget jobs exited with errors"
    echo "Re-attemping to download $exitStats files..." && sleep 3s
    parallel --retry-failed --joblog $joblog

    # contCommands=$( tail -n +2 $joblog | awk '$7 != 0' | awk 'BEGIN {FS="\t"} {print $9}' | sed 's/wget --quiet/wget -c/g' )
    # contJoblog="$(echo ${joblog} | cut -f1 -d'.')_continue.log"
    # parallel --jobs $nJob --joblog $contJoblog {1} ::: $contCommands && sleep $sleeptime
    # contStats=$( tail -n +2 $contJoblog | awk '$7 != 0' | wc -l )
    # if [ "$contStats" -gt 0 ]
    # then
    #     echo "WARNING: $contStats re-attempted wget jobs exited with errors"
    # else
    #     echo "No re-attempted wget jobs exited with errors"
    #     echo "Download re-attempts: DONE"
    # fi
else
    echo "No wget jobs exited with errors"
    echo "DONE"
fi

#Do md5 checksums match?
printf "\n################ Downloading md5sum files ###############\n"
prefix=$(basename "$genomeTable" | cut -f1 -d'.')
md5sum_file=${prefix}_md5sum.out
md5sum_check=${prefix}_md5sumchk.out
awk '{FS=","} {print $15}' $genomeTable | grep 'ftp' | sed 's/"//g' | \
    parallel --jobs $nJob --joblog md5sum.log "wget --quiet -O - {}/md5checksums.txt && sleep $sleeptime" | \
    grep '_genomic.fna.gz' | \
    grep -v 'rna_from\|cds_from' > $md5sum_file

printf "\n################ Download finished ###############\n"
printf "\n################ Verifying downloaded files using  md5sum ###############\n"

# verification of downloaded files
md5sum --quiet -c $md5sum_file > $md5sum_check

if [ -s "$md5sum_check" ]
then
    nfail=$(wc -l "$md5sum_check" | cut -f1 -d' ')
    printf "\n WARNING !! $nfail downloaded files were not OK. Check $md5sum_check for the list of filenames \n"
else
    printf "\nDownloaded files verified. All OK ! \n"
fi

exit 0
