#!/bin/bash



starttime=$(date +%Y-%m-%dT%X)

echo "starttime set to ${starttime}"
echo


while true; do
    
    echo

    raw=$(sacct -p --noheader --delimiter="" -S ${starttime} --format State) #| cut -d"," -f 6 #| grep -vE ""

    #echo ${raw}

    n_pending=$(echo $raw | grep -o PENDING | wc -l)
    n_running=$(echo $raw | grep -o RUNNING | wc -l)
    n_failed=$(echo $raw | grep -o FAILED | wc -l)
    n_completed=$(echo $raw | grep -o COMPLETED | wc -l)
    n_sum=$(echo $raw | wc -w)




    echo "pending:    ${n_pending}"
    echo "running:    ${n_running}"
    echo "failed:     ${n_failed}"
    echo "completed:  ${n_completed}"
    echo "sum: ${n_sum}"



    if [ $n_failed != "0" ]; then
        echo "One or more jobs has failed:"
        echo $(sacct -p --noheader --delimiter="" -S ${starttime} --format JobID,State | grep "FAILED")
        echo "exiting"
        break
    fi


    if [ $n_pending == "0" ] && [ $n_running == "0" ] && [ $n_failed == "0" ]; then 
        echo "All jobs are done. Exiting while loop"
        break

    fi



    sleep 8



done