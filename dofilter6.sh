#!/bin/bash
# 2016-10-18 isradelacon@gmail.com

contaminants="/home/barrantes/data/contaminant_list.fa"

# decompress if needed
echo "$0 and gzip start $(date)"
gzip -d --force *.fastq.gz

# filter PhiX
echo "filterPhageSE.sh start $(date)"
filterPhageSE.sh

# adapter removal
# parallel run
pids=""
RESULT=0

# scythe
echo "scythe start $(date)"
for file in *.nophix.fastq; do
 scythe --quiet -a $contaminants -o $file.scythe.fastq $file &
 pids="$pids $!"
done

# waits until each parallel process ends
# source http://stackoverflow.com/a/26240420/1588082
for pid in $pids; do
    wait $pid || let "RESULT=1"
done

if [ "$RESULT" == "1" ];
    then
       exit 1
fi

rename 's/\.nophix\.fastq//' *.scythe.fastq

# quality filtering
# parallel run
pids=""
RESULT=0

# sickle
echo "sickle start $(date)"
for file in *.scythe.fastq; do
 sickle se --qual-type sanger --fastq-file $file --output-file $file.sickle.fastq &
 pids="$pids $!"
done

# waits until each parallel process ends
# source http://stackoverflow.com/a/26240420/1588082
for pid in $pids; do 
    wait $pid || let "RESULT=1"
done

if [ "$RESULT" == "1" ];
    then
       exit 1
fi

rename 's/\.scythe\.fastq//' *.sickle.fastq
echo "$0 end $(date)" 
