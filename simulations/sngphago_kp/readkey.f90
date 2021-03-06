SUBROUTINE READKEY
  ! this subroutine reads in keywords from a parameter file
  ! it sets the various global variables defined in KEYS module
  ! name of the parameter file is param.* where * is a keyword argument
  ! if no keyword argument is supplied, the default is just a file called param
  ! The EXTRAPARAMFILES keyword will allow extra parameter files to be 
  ! read in as well

  USE KEYS
  USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO
  USE GENUTIL

  IMPLICIT NONE

  ! ---- stuff for inputing the parameter file in free format --------
  CHARACTER*100 :: ARG ! command line argument
  INTEGER :: NUMARG ! number of command line arguments
  INTEGER :: NITEMS ! number of items on the line in the parameter file
  INTEGER :: PF ! input file unit
  LOGICAL :: FILEEND=.FALSE. ! done reading file?
  CHARACTER*100 :: WORD ! keyword
  ! -------------- for reading multiple parameter files --------  
  INTEGER, PARAMETER :: MAXNFILES = 10
  CHARACTER*100 :: PARAMFILES(MAXNFILES)
  INTEGER :: NPARAMFILES, NPARAMREAD
  ! ------ for initializing random number generator
  INTEGER :: TIMEVAL(8), SEED
  ! ---------------- temporary variables ---------------
  INTEGER :: DUMI, I, TMPI, DUMI1, DUMI2, DUMI3
  CHARACTER*100 :: DUMSTR
  LOGICAL :: LDUM

  ! ------------------------
  ! set variable defaults
  ! ------------------------
  ACTION = 'NONE'
  RNGSEED = 0
  VERBOSE = .FALSE.


  ! input/output  
  OUTFILE = '*.out'
  OUTPUTEVERY = 1
  DUMPSNAPSHOTS = .FALSE. ! periodically dump chain snapshots
  APPENDSNAPSHOTS = .FALSE. 
  SNAPSHOTEVERY = 1D4 ! how often to dump snapshots
  SNAPSHOTFILE = '*.snap.out' ! snapshot file
  PRINTEVERY = 1 ! how often to print output
  WHICHSNAPSHOT = 1 ! which iteration to snapshot

  ! kinetic rates
  PS = 0D0
  KW = -1D0
  KDECAY = 1D0
  KPROD = 1D0
  MODELFLAG = 0 !default model is as many sinks allowed 

  ! unit protein content
  UNITPROT = 1D0
  ! mito velocity
  VEL = 1D0
  ! length of domain
  L = 5D2
  ! number of sinks if provided
  NSINKS = 1
  ! avg SM if provided
  AVGSM = -1
  ! max number of stationary mito per region
  SM = -1
  ! width of sinks
  SINKWIDTH = 1
  ! number of demand regions
  NREGIONS = -1
  ! size of newly born mitochondria
  MITOSIZE = 1
  ! total mito in the domain
  TOTMITOINL = -1

  ! running dynamics
  NSTEP = 1000
  DELT = 1D-3
  NITER = 10

  ! max allowed clusters
  MAXNCLUST = 1D5

  ! kill mitochondria if protein gets below this level
  DEATHPROTLEVEL = 0D0

  ! fill up all the regions before spreading remaining sinks randomly
  FILLUPREGIONS = .FALSE.


  ! -------------------------
  ! Read in all parameter files, starting with the ones specified on command line
  ! --------------------------

  PF = 55 ! i/o unit number to be used for parameter files

  ! get input parameter files from command line
  NPARAMFILES = 0
  NUMARG = COMMAND_ARGUMENT_COUNT()  
  IF (NUMARG==0) THEN
     NPARAMFILES = 1
     PARAMFILES(1) = 'param'
     ARG = ''
  ELSE
     DO I = 1,NUMARG
        CALL GETARG(I, ARG)
        NPARAMFILES = NPARAMFILES + 1
        WRITE(DUMSTR,'(A)') 'param.' //TRIM(ADJUSTL(ARG))
        PARAMFILES(NPARAMFILES) = DUMSTR
     ENDDO
     ! reset arg to its original value
     IF (NUMARG.GT.1) CALL GETARG(1,ARG)
  ENDIF

  NPARAMREAD = 0 ! keep track of how many files have been read
  DO WHILE (NPARAMREAD.LT.NPARAMFILES)
     NPARAMREAD = NPARAMREAD + 1

     PRINT*, 'Reading parameter file: ', PARAMFILES(NPARAMREAD)
     INQUIRE(FILE=PARAMFILES(NPARAMREAD),EXIST=LDUM)
     IF (.NOT.LDUM) THEN
        PRINT*, 'ERROR in READKEY: Parameter file ', TRIM(ADJUSTL(PARAMFILES(NPARAMREAD))), ' does not exist.'
        STOP 1
     ENDIF
     OPEN(UNIT=PF, FILE=PARAMFILES(NPARAMREAD), STATUS='OLD')

     ! read in the keywords one line at a time
     DO 
        CALL READLINE(PF,FILEEND,NITEMS)
        IF (FILEEND.and.nitems.eq.0) EXIT

        ! skip empty lines
        IF (NITEMS.EQ.0) CYCLE

        ! Read in the keyword for this line
        CALL READA(WORD,CASESET=1)

        ! Skip any empty lines or any comment lines
        IF (WORD(1:1).EQ.'#') CYCLE

        SELECT CASE(WORD) ! pick which keyword
        CASE('ACTION')
           CALL READA(ACTION, CASESET=1) ! force upper case
        CASE('AVGSM')
           CALL READF(AVGSM)
        CASE('DEATHPROTLEVEL')
           CALL READF(DEATHPROTLEVEL)
        CASE('DELT')
           CALL READF(DELT)
        CASE('FILLUPREGIONS')
           CALL READO(FILLUPREGIONS)
        CASE('KDECAY')
           CALL READF(KDECAY)
        CASE('KPROD')
           CALL READF(KPROD)
        CASE('KW')
           CALL READF(KW)
        CASE('L')
           CALL READF(L)
        CASE('MAXNCLUST')
           CALL READI(MAXNCLUST)
        CASE('MITOSIZE')
           CALL READF(MITOSIZE)
        CASE('MODELFLAG')
           CALL READI(MODELFLAG)
        CASE('NITER')
           CALL READI(NITER)
        CASE('NREGIONS')
           CALL READI(NREGIONS)
        CASE('NSTEP')
           CALL READI(NSTEP)
        CASE('NSINKS')
           CALL READF(NSINKS)
        CASE('OUTFILE')
           CALL READA(OUTFILE)
           IF (NITEMS.GT.2) CALL READI(OUTPUTEVERY)
        CASE('OUTPUTEVERY')
           CALL READI(OUTPUTEVERY)
        CASE('PRINTEVERY')
           CALL READI(PRINTEVERY)
        CASE('PS')
           CALL READF(PS)
        CASE('RNGSEED')
           CALL READI(RNGSEED)
        CASE('SINKWIDTH')
           CALL READF(SINKWIDTH)
        CASE('SM')
           CALL READI(SM)
        CASE('SNAPSHOTFILE')
           CALL READA (SNAPSHOTFILE)
        CASE('SNAPSHOT')
           ! DUMPSNAPSHOTS = .TRUE.           
           IF (NITEMS.GT.1) CALL READI(SNAPSHOTEVERY)
           IF (NITEMS.GT.2) CALL READA(SNAPSHOTFILE)
           IF (NITEMS.GT.3) CALL READI(WHICHSNAPSHOT)
           IF (NITEMS.GT.4) CALL READO(APPENDSNAPSHOTS)
        CASE('SNAPSHOTEVERY')
           CALL READI(SNAPSHOTEVERY)
        CASE('TOTMITOINL')
           CALL READF(TOTMITOINL)
        CASE('UNITPROT')
           CALL READF(UNITPROT)
        CASE('VEL')
           CALL READF(VEL)
        CASE('VERBOSE')
           CALL READO(VERBOSE)
        CASE DEFAULT
           print*, 'ERROR: unidentified keyword ', TRIM(WORD), " Will ignore."
        END SELECT
     ENDDO
     CLOSE(PF)
  ENDDO


  ! -----------------
  ! check validity of some values, raise errors or adjust as necessary
  ! -----------------  

  IF (PS.LT.0 .or. PS.GT.1) THEN
     PRINT*, 'ERROR IN PS VALUE', PS
     STOP 1
  ENDIF

  IF (KDECAY.LT.0) THEN
     PRINT*, 'ERROR IN KDECAY VALUE', KDECAY
     STOP 1
  ENDIF



!!$  ! Assign production rate to give the desired total number of mitochondria
!!$  ! only assign it if preassigned value is negative 
!!$  IF (KPROD .LT. 0) THEN
!!$
!!$     IF (TOTMITOINL.GT.0) THEN
!!$        KPROD = (TOTMITOINL - REAL(NSINKS)) * VEL / (2*L)
!!$     ENDIF
!!$     
!!$  ENDIF

  ! Update demand regions
  IF (NREGIONS.LT.0) THEN
     NREGIONS = INT(NSINKS)
  ENDIF

  ! assign value of kw is AVGSM is provided
  ! otherwise set AVGSM to NSINKS; one of them has to be provided

  IF (AVGSM .LT. 0) THEN
     AVGSM = NSINKS/NREGIONS
  END IF




  IF (KPROD .GT. 0) THEN

     IF (TOTMITOINL.GT.0) THEN
        KPROD = (TOTMITOINL - REAL(NREGIONS*AVGSM)) * VEL / (2*L)
        IF (KPROD.LT.0) THEN
           PRINT*, 'ERROR IN SETTING KPROD:', KPROD
           STOP 1
        ENDIF
     ENDIF

  ENDIF

  IF (KW .LT. 0) THEN
     IF (TOTMITOINL .GT. 0) THEN
        KW = PS*VEL*(TOTMITOINL - NREGIONS*AVGSM)/ L / AVGSM
     ELSE
        KW = 2/AVGSM * PS * KPROD
     END IF

  END IF








  ! ----------- fix file names -----------
  CALL REPLACESUBSTR(OUTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(SNAPSHOTFILE,'*',TRIM(ADJUSTL(ARG)))
  ! ---------------------------

  ! Initiate random number generator 
  IF (RNGSEED.EQ.0) THEN
     ! use the current time of day in milliseconds
     CALL DATE_AND_TIME(VALUES=TIMEVAL)
     SEED = TIMEVAL(5)*3600*1000 + TIMEVAL(6)*60*1000 + TIMEVAL(7)*1000 + TIMEVAL(8)
  ELSEIF (RNGSEED.EQ.-1) THEN
     ! use the last 5 characters in the command-line argument
     SEED = STRING2NUM(TRIM(ADJUSTL(ARG)))    
  ELSEIF (RNGSEED.EQ.-2) THEN
     ! use the last 4 characters in the command-line argument 
     ! and additionally the millisecond time 
     CALL DATE_AND_TIME(VALUES=TIMEVAL)
     SEED = STRING2NUM(TRIM(ADJUSTL(ARG)),TIMEVAL(8))
  ELSE
     ! use this seed directly
     SEED = RNGSEED
  ENDIF

  print*, 'Initiating Mersenne twister random number generator with seed:', SEED
  CALL SGRND(SEED)

  print*, '------------Parameter values : -------------------'
  print*, 'ACTION: ', TRIM(ADJUSTL(ACTION))
  print*, 'Output file: ', TRIM(OUTFILE)  
  IF (DUMPSNAPSHOTS) THEN
     PRINT*, 'Dumping snapshot every', SNAPSHOTEVERY,'steps. In file:', TRIM(ADJUSTL(SNAPSHOTFILE))
     print*, 'Snapshot iteration(s):', WHICHSNAPSHOT
  ENDIF

  print*, 'Number of mito:', TOTMITOINL
  print*, 'Avg SM:', AVGSM
  print*, 'kw, ps,kdecay,vel,Nregions,kprod:', &
       & KW,PS,KDECAY,VEL,NREGIONS,KPROD
  print*, 'NITER:', niter
  print*, '----------------------------------------------------'


END SUBROUTINE READKEY
