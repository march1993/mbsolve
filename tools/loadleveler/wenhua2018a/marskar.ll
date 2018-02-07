#!/bin/bash
#@ wall_clock_limit = 24:00:00
#@ job_name = mbsolve
#@ job_type = parallel
#@ class = fat
#@ node = 1
#@tasks_per_node = 1
#@ node_usage = not_shared
#@ initialdir = $(home)/simulations/
#@ output = wenhua2018a-marskar-$(schedd_host).$(jobid).$(stepid).out
#@ error = wenhua2018a-marskar-$(schedd_host).$(jobid).$(stepid).out
#@ notification=always
#@ queue
#@ notify_user = wenhua.shi@tum.de

. /etc/profile
. /etc/profile.d/modules.sh

module load boost/1.61_icc

thread_s=1
thread_i=39
thread_e=40
iterations=5

# gridpoints=32768
# endtime=4e-12
name=wenhua2018a-marskar
method=openmp-6lvl-os-red
device=marskar2011multilevel

# vary thread count
for threads in `seq $thread_s $thread_i $thread_e`; do

# reproducibility
for it in `seq 1 $iterations`; do

out_dir=$name-$LOADL_STEP_ID/$threads/$it/

mkdir -p $out_dir

echo "Thread count: " $threads

KMP_AFFINITY=granularity=fine,proclist=[`seq -s , 0 $(($threads - 1))`],explicit OMP_NUM_THREADS=$threads ../build-wenhua2018a/mbsolve-tool/mbsolve-tool -m $method -d $device -w matlab -o $out_dir/$name.mat

done

done

# extract performance data
cat $name-$LOADL_STEP_ID.out | ../mbsolve/tools/loadleveler/generate_plot.sh > $name-$LOADL_STEP_ID/perf-$name.csv

mv $name-$LOADL_STEP_ID.out $name-$LOADL_STEP_ID/
