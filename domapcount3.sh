#!/bin/bash
# mapping and counts
# USAGE: domapcount3.sh
# 2015-10-18 isradelacon@gmail.com
echo "$0 start $(date)"

# mapping variables for STAR
DEFAULT_STAR_PARAMS="--outReadsUnmapped Fastx --runThreadN 24"
MAPPING_STAR_PARAMS="--alignEndsType EndToEnd --alignIntronMax 1 --alignIntronMin 2 --scoreDelOpen -10000 --scoreInsOpen -10000 --outFilterMultimapNmax 100 --outFilterMismatchNoverLmax 0.1"
EXPRESS_MAPPING_TARGET="SPY2001.ffn"

# variables for running eXpress in parallel
pids=""
RESULT=0
blank=""

# mapping
if [ ! -d GenomeDir ]; then
  echo "STAR Genome file not found!"
  exit 1
else
  echo "mapping start $(date)"
  ls -1 *.sickle.fastq | sed -e 's/.sickle.fastq//;' | xargs -I{} STAR $DEFAULT_STAR_PARAMS $MAPPING_STAR_PARAMS --readFilesIn {}.sickle.fastq --outFileNamePrefix {}
  echo "mapping complete $(date)"
fi

# counting on SAMs
ls -1 *.sam > /dev/null 2>&1
if [ "$?" = "0" ]; then
 echo "count start $(date)"

 for file in *.sam; do 
  # replace string
  # source http://stackoverflow.com/a/13210909/1588082
  file="${file/Aligned.out.sam/$blank}"
  express --f-stranded --output-dir mix.$file $EXPRESS_MAPPING_TARGET $file\Aligned.out.sam &
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
 echo "count end $(date)"
 
 echo "merge XPRS outputs $(date)"
 find . -name 'results.xprs' | perl -lne '$orfile = $_; s{^\.\/}{}; s{\/results}{}; system "ln -s $orfile $_";'
 mergexprs.sh
 mv allxprs.tsv mix.allxprs.tsv
 mergetpmxprs.sh 
 mv alltpmxprs.tsv mix.alltpmxprs.tsv
 dosummarySTAR.sh mix-mapping-summary-$(date '+%Y-%m-%d').tsv
 tar czf map-$(date '+%Y-%m-%d').tar.gz *.tab *.out --remove-files
# rm -rf *.sam
 echo "all processes complete $(date)"
else
    echo "SAM files not found!"
    exit 1
fi

# report on Excel format
ln -s mix.allxprs.tsv readcounts
ln -s mix.alltpmxprs.tsv tpmnormalization
for file in readcounts tpmnormalization; do tab2xls2.pl $file $file.xls; done
xlsmerge.pl -s -o mapping.xls readcounts.xls tpmnormalization.xls
for file in readcounts tpmnormalization; do unlink $file; rm -rf $file.xls; done

echo "$0 end $(date)"
