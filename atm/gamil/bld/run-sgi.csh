#! /usr/bin/csh -f

#-----------------------------------------------------------------------
## SGI
##------------
##
## This is an example script to build and run the default GAMIL configuration on
## an SGI.  The default configuration is T42L26, Eulerian dynamics, CLM2 land
## model, and CSIM4 ice model.
##
## Set NQS options 
## -q  is the queue in which to run the job
## -l  on NCAR's SGI this number should be 4 greater than the number of CPUs 
##       specified on the -q line
## -lT is the time limit in seconds
## -eo says to include both stdout and stderr in the output file

#QSUB -q ded_32
#QSUB -l mpp_p=36
#QSUB -lT 20000
#QSUB -eo

## Number of OMP threads.  Only need to set if requesting
## fewer than the maximum number for the queue.
#setenv OMP_NUM_THREADS 32

## Set runtime variables for optimal execution on a DSM
setenv OMP_DYNAMIC FALSE
setenv _DSM_PLACEMENT ROUND_ROBIN
setenv _DSM_WAIT SPIN
setenv MPC_GANG    OFF

## MP_SLAVE_STACKSIZE sets the size of the thread stack
setenv MP_SLAVE_STACKSIZE 40000000
## Do our best to get sufficient stack memory
limit stacksize unlimited

## SGI-specific stuff. MIPSpro is required for f90
source /opt/modules/modules/init/csh
module purge
module load MIPSpro modules nqe mpt

## ROOT OF GAMIL DISTRIBUTION - probably needs to be customized.
## Contains the source code for the GAMIL distribution.
## (the root directory contains the subdirectory "models")
set gamilroot      = /fs/cgd/data0/$LOGNAME/gamil1.0

## ROOT OF GAMIL DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the GAMIL distribution.
## (the root directory contains the subdirectories "atm" and "lnd")
setenv CSMDATA    /fs/cgd/csm/inputdata

## Default namelist settings:
## $case is the case identifier for this run. It will be placed in the namelist.
## $runtype is the run type: initial, restart, or branch.
## $nelapse is the number of timesteps to integrate, or number of days if negative.
set case         = run
set runtype      = initial
set nelapse      = -1

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the GAMIL configuration scripts.
set wrkdir       = /ptmp/$LOGNAME
set blddir       = $wrkdir/$case/bld
set rundir       = $wrkdir/$case
set cfgdir       = $gamilroot/models/atm/gamil/bld

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## If an executable doesn't exist, build one.
if ( ! -x $blddir/gamil ) then
    cd $blddir                  || echo "cd $blddir failed" && exit 1
    $cfgdir/configure           || echo "configure failed" && exit 1
    echo "building GAMIL in $blddir ..."
    rm -f Depends
    gmake -j4 >&! MAKE.out      || echo "GAMIL build failed: see $blddir/MAKE.out" && exit 1
endif

## Create the namelist
cd $blddir                      || echo "cd $blddir failed" && exit 1
$cfgdir/build-namelist -s -case $case -runtype $runtype -o $rundir/namelist \
 -namelist "&atmexp nelapse=$nelapse /"  || echo "build-namelist failed" && exit 1

## Run GAMIL
cd $rundir                      || echo "cd $rundir failed" && exit 1
echo "running GAMIL in $rundir"
$blddir/gamil < namelist          || echo "GAMIL run failed" && exit 1

exit 0
