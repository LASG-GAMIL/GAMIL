#! /bin/tcsh -f

#-----------------------------------------------------------------------
## PC-linux
##------------

## nthreads is the number of Open-MP threads.  On the PC we assume a 
## pure shared-memory configuration which assigns 1 thread per processor.
set nthreads = 2

## Do our best to get sufficient stack memory
limit stacksize unlimited

## ROOT OF GAMIL DISTRIBUTION - probably needs to be customized.
## Contains the source code for the GAMIL distribution.
## (the root directory contains the subdirectory "models")
set gamilroot      = /fs/cgd/data0/$LOGNAME/gamil1.0

## ROOT OF GAMIL DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the GAMIL distribution.
## (the root directory contains the subdirectories "atm" and "lnd")
setenv CSMDATA     /fs/cgd/csm/inputdata

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
    $cfgdir/configure -smp      || echo "configure failed" && exit 1
    echo "building GAMIL in $blddir ..."
    rm -f Depends
    gmake -j2 >&! MAKE.out      || echo "GAMIL build failed: see $blddir/MAKE.out" && exit 1
endif

## Create the namelist
cd $blddir                      || echo "cd $blddir failed" && exit 1
$cfgdir/build-namelist -s -case $case -runtype $runtype -o $rundir/namelist \
 -namelist "&atmexp nelapse=$nelapse /"  || echo "build-namelist failed" && exit 1

## Run GAMIL
cd $rundir                      || echo "cd $rundir failed" && exit 1
echo "running GAMIL in $rundir"
env OMP_NUM_TASKS=$nthreads MPSTKZ="128M" $blddir/gamil < namelist  || echo "GAMIL run failed" && exit 1

exit 0
