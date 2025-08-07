#!/bin/bash
#PBS -N gbm_spheroid
#PBS -l select=1:ncpus=12:mem=10gb:scratch_local=10gb
#PBS -l walltime=24:00:00
# The 4 lines above are options for the scheduling system: the job will run 1 hour at maximum, 1 machine with 4 processors + 4gb RAM memory + 10gb scratch memory are requested

# define a DATADIR variable: directory where the input files are taken from and where the output will be copied to
DATADIR=/your/path # substitute username and path to your real username and path

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of the node it is run on, and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails, and you need to remove the scratch directory manually
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" >> $DATADIR/jobs_info.txt

SESSION_ID=$(echo $PBS_JOBID | cut -d. -f1)

mkdir -p $DATADIR/sim_backup/${PBS_JOBNAME}_${SESSION_ID} || { echo >&2 "Failed to create backup directory!"; exit 4; }

(
  while true; do
    sleep 60
    echo "Backing up sim_output to $DATADIR/sim_backup at $(date)"
    cp -r $SCRATCHDIR/sim_output $DATADIR/sim_backup/${PBS_JOBNAME}_${SESSION_ID}
  done
) &



module add julia/1.11.3

export JULIA_NUM_THREADS=12
export JULIA_PKGDIR=$DATADIR/.julia
export JULIA_DEPOT_PATH=$DATADIR/.julia
export JULIA_PROJECT=$DATADIR/spheroid_model_cluster

# test if the scratch directory is set
# if scratch directory is not set, issue error message and exit
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

cp -r $DATADIR/spheroid_model_cluster/scripts  $SCRATCHDIR || { echo >&2 "Error while copying input file(s)!"; exit 2; }
cp -r $DATADIR/spheroid_model_cluster/src  $SCRATCHDIR || { echo >&2 "Error while copying input file(s)!"; exit 2; }

cp -r $DATADIR/spheroid_model_cluster/Project.toml $SCRATCHDIR
cp -r $DATADIR/spheroid_model_cluster/Manifest.toml $SCRATCHDIR


# move into scratch directory
cd $SCRATCHDIR || { echo >&2 "Failed to cd into $SCRATCHDIR"; exit 2; }

julia scripts/simulations_runnable.jl $ARGS || { echo >&2 "Calculation ended up erroneously (with a code $?) !!"; exit 3; }

# move the output to user's DATADIR or exit in case of failure
if [ ! -d "$DATADIR"/results/ ]; then
        mkdir "$DATADIR"/results
fi

ls -l $SCRATCHDIR > output_text.txt
cp output_text.txt $DATADIR/results

# move the output to user's DATADIR or exit in case of failure
# create a subdirectory for this job's results
mkdir -p $DATADIR/results/${PBS_JOBNAME}_${SESSION_ID} || { echo >&2 "Failed to create result directory!"; exit 4; }

# copy the sim_output to that directory
cp -r $SCRATCHDIR/sim_output/* $DATADIR/results/${PBS_JOBNAME}_${SESSION_ID} || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }


# clean the SCRATCH directory
clean_scratch
