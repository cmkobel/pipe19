#!/bin/bash

# Any subsequent(*) commands which fail will cause the shell script to exit immediately
set -e

# Ensure existence of the directory for logs.



default_waitingtime=8
log_dir="controller_log"
main_log="main_log.txt"


mkdir -p $log_dir

function get_time() {
    echo $(date +%Y-%m-%dT%X)
}

function enlog() {
    echo "[$(get_time)] ${1}" | tee -a $main_log
}




# This loop waits for a pipeline step to complete.
function wait_for_jobs_to_complete() {

    parameter_starttime=$1

    echo "function called with time $1"

    while true; do
        
        echo

        raw=$(sacct -p --noheader --delimiter="" -S ${parameter_starttime} --format State) #| cut -d"," -f 6 #| grep -vE ""

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


        # Exit erroneously if one or more jobs fail
        if [ $n_failed != "0" ]; then
            echo "One or more jobs has failed:"
            echo $(sacct -p --noheader --delimiter="" -S ${parameter_starttime} --format JobID,State | grep "FAILED")
            echo "exiting"
            return 1
        fi

        # Exit gracefully if all jobs complete
        if [ $n_pending == "0" ] && [ $n_running == "0" ] && [ $n_failed == "0" ]; then 
            echo "All jobs are done. Exiting while loop"
            return 0

        fi



        sleep $default_waitingtime



    done
}



# Master loop that runs forever
while true; do


    level="pipe19"
    enlog "Entering level: ${level}"

    gwf status --status shouldrun  > ${log_dir}/gwf_status_${level}_shouldrun.stdout 2> ${log_dir}/gwf_status_${level}_shouldrun.stderr

    # Check if stderr above is non-empty:
    if [ -s ${log_dir}/gwf_status_${level}_shouldrun.stderr ]; then

        #echo "Jobs present. Calling gwf ..."
        enlog "$(cat ${log_dir}/gwf_status_shouldrun.stderr | wc -l) jobs to submit"
        enlog "Calling gwf ..."
        gwf run > ${log_dir}/gwf_run.stdout 2> ${log_dir}/gwf_run.stderr

        wait_for_jobs_to_complete $(get_time) && echo "yes" || echo "no"

    else
        enlog "No jobs yet. Waiting ..."
        echo && sleep 60
    fi










    







done

