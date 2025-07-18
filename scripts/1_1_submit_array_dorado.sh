#!/bin/bash

#variables 
PROJECT=$1
LIB=$2
POD5_DIR=$SCRATCH/projects/$PROJECT/raw_data/pod5/$PROJECT'_'$LIB

ls $POD5_DIR/*pod5 > $POD5_DIR/$PROJECT'_'$LIB'_'files.txt
N_files=$( wc -l < $POD5_DIR/$PROJECT'_'$LIB'_'files.txt )

split -d -l 500 $POD5_DIR/$PROJECT'_'$LIB'_'files.txt $POD5_DIR/$PROJECT'_'$LIB'_'part
rm $POD5_DIR/$PROJECT'_'$LIB'_'files.txt

for i in $POD5_DIR/$PROJECT'_'$LIB'_'part*
do
    echo "$i"
    sbatch dorado.sh $PROJECT $LIB $i
    sleep 10
done

## Define start (in most cases this should be 0, but if you, for example, have donw the first 1500 you can indicate that here )
#start_after=0
## Define end
#end_at=2203
## Define the batch size
#batch_size=250

## Calculate the total number of tasks
#total_tasks=$(( end_at - start_after )) 
## Calculate the number of batches needed
#num_batches=$(( total_tasks / batch_size ))
## Calculate size of last batch
#last_batch=$((  total_tasks-(num_batches*batch_size) ))

## Loop through the batches and submit sbatch jobs
#for ((i=0; i<num_batches; i++)); do
#    start_index=$((i * batch_size + start_after))
#    echo "submitted  $batch_size jobs, starting after $start_index"
#    sbatch --array=1-$batch_size 3_dorado.sh $PROJECT $LIB $start_index
#    sleep 10
#done

## Submit any remaining tasks (if the total_tasks is not evenly divisible by batch_size)
#if ((total_tasks % batch_size != 0)); then
#    start_index=$((num_batches * batch_size + start_after))
#    echo "submitted  $last_batch jobs, starting after $start_index"
#    sbatch --array=1-$last_batch 3_dorado.sh $PROJECT $LIB $start_index
#fi