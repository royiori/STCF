#!/bin/bash
########################################################
#
#
########################################################
set +x
current=`pwd`
temsh="./template.sh"
temcmd="./template.cmd"
filename="FileName"
jobs_path="${current}/jobs/"
exes_path="${current}/../RICH"
args_path="${current}/../sim/stcf-5mm-batch.txt"
if [ ! -d $jobs_path ] ; then
        mkdir -p $jobs_path
fi

for num_i in `seq 1 102` ; do
#num_i=RUN
        job_path=`echo $jobs_path | sed "s/\//@/g"`
        exe_path=`echo $exes_path | sed "s/\//@/g"`
        arg_path=`echo $args_path | sed "s/\//@/g"`
        #cat $temsh > $jobs_path"/"$filename$num_i.sh
        cat $temcmd | sed "s/EXECPATH/$exe_path/g" | sed "s/ARGU1/$arg_path/g" | sed "s/ARGU2/$num_i/g" | sed "s/OUT/$job_path$filename$num_i.out/g" | sed "s/ERR/$job_path$filename$num_i.err/g"| sed "s/LOG/$job_path$filename$num_i.log/g"| sed "s/@/\//g" > $jobs_path"/"$filename$num_i.cmd
        #condor_submit $jobs_path$filename$num_i.cmd
done
