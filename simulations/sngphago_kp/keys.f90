MODULE KEYS
  ! keyword parameters that are globally used in many different places in the code
  IMPLICIT NONE

   ! -------- General program control ---------------
  CHARACTER*100 :: ACTION
  ! Random generator seed
  ! positive value uses that exact seed
  ! 0: uses current time in milliseconds
  ! -1:  last 5 characters of commandline argument
  ! -2: last 4 characters and time in milliseconds  
  INTEGER :: RNGSEED
  LOGICAL :: VERBOSE  
  
  ! ----------------------
  ! Output / input
  ! -----------------------
  CHARACTER*100 :: OUTFILE, SNAPSHOTFILE
  ! dumpsnapshot: whether we dump snapshots
  LOGICAL :: DUMPSNAPSHOTS, APPENDSNAPSHOTS
  ! how often to dump snapshots
  INTEGER :: SNAPSHOTEVERY, OUTPUTEVERY, PRINTEVERY
  ! which iteration to take snapshots of
  INTEGER :: WHICHSNAPSHOT

  ! ------------
  ! Parameters for cluster dynamics
  ! ----------------
  ! max allowed clusters
  INTEGER :: MAXNCLUST
  ! fission and fusion rates
  DOUBLE PRECISION :: PS, KW, AVGSM
  DOUBLE PRECISION :: KPROD ! production rate
  DOUBLE PRECISION :: UNITPROT !unit protein content
  DOUBLE PRECISION :: VEL ! mito velocity
  DOUBLE PRECISION :: KDECAY ! protein decay rate
  DOUBLE PRECISION :: L, SINKWIDTH ! length of domain and width of sinks
  DOUBLE PRECISION :: MITOSIZE
  DOUBLE PRECISION :: NSINKS ! Number of sinks
  INTEGER :: SM
  INTEGER :: NREGIONS ! Number of demand regions
  DOUBLE PRECISION :: TOTMITOINL !total mito in L, used to define velocity
  INTEGER :: MODELFLAG ! which shuttle model will run, 0 is as many sinks
  DOUBLE PRECISION :: DEATHPROTLEVEL ! minimum protein below which mito dies


  !-----------
  ! parameters for running dynamics
  ! ----------
  INTEGER :: NSTEP ! total steps to run
  DOUBLE PRECISION :: DELT ! time step
  INTEGER :: NITER ! number of iterations 

  ! fill up all regions with 1 sink before placing remainder randomly
  LOGICAL :: FILLUPREGIONS
END MODULE KEYS
