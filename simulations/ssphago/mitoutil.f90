MODULE MITOUTIL
  IMPLICIT NONE

  TYPE CLUSTERSET
     INTEGER :: MAXNCLUST ! max allowed number of clusters
     INTEGER :: MAXCLUSTID !
     INTEGER :: MAXTAG ! maximum id tag that has ever been used in this set of clusters
     DOUBLE PRECISION, POINTER :: POS(:) ! position of cluster
     DOUBLE PRECISION, POINTER :: SIZE(:) ! size of cluster
     INTEGER, POINTER :: STATE(:) ! stopped or walking
     INTEGER, POiNTER :: IDTAG(:) ! permanent id tag associated with each cluster, to keep track of individual trajectories
     
! TODO: get rid of PFORWARD, implement states:
! -2) stopped and planning to move backward
! -1) moving backward
! 0) stopped
! 1) moving forward
! 2) stopped and planning to move forward
     DOUBLE PRECISION, POINTER :: PROT(:) ! protein content
     DOUBLE PRECISION, POINTER :: PFORWARD(:) ! probability of moving forward
     
     LOGICAL, POINTER :: EXISTS(:) ! does this cluster currently exist?
     ! first available spot for a new cluster
     INTEGER :: FIRSTEMPTY, NCLUST
     
     LOGICAL :: ARRAYSET=.FALSE.

     DOUBLE PRECISION :: MITOSIZE
     DOUBLE PRECISION :: KFMOV, KFSTAT, KFISSION, KPROD, UNITPROT, &
          & VEL, KDECAY, PAUSERATE, WALKRATE, L, SINKWIDTH
     DOUBLE PRECISION :: NSINKS

      ! dead mitochondria are no longer capable of fusion
     LOGICAL, POINTER :: DEAD(:) ! is this cluster currently dead (incapable of fusion)

  END type CLUSTERSET

  
  TYPE CLUSTER
     DOUBLE PRECISION :: POS ! position of cluster
     DOUBLE PRECISION :: SIZE ! size of cluster
     DOUBLE PRECISiON :: STATE ! stopped or walking
     DOUBLE PRECISION :: PROT ! protein content
     DOUBLE PRECISION :: PFORWARD ! probability of moving forward
  END type CLUSTER
  
CONTAINS
  
  SUBROUTINE MAKENEWCLUSTER(CLUSTP, CNEW)
    IMPLICIT NONE
    TYPE(CLUSTERSET), POINTER :: CLUSTP
    INTEGER, INTENT(OUT) :: CNEW
    INTEGER :: CC
    LOGICAL :: FIRSTSET
    
    CNEW = CLUSTP%FIRSTEMPTY

    IF (CNEW.GT.CLUSTP%MAXCLUSTID) THEN
       CLUSTP%MAXCLUSTID = CLUSTP%MAXCLUSTID + 1
    ENDIF
    
    CLUSTP%EXISTS(CLUSTP%FIRSTEMPTY) = .TRUE.
    CLUSTP%DEAD(CLUSTP%FIRSTEMPTY) = .FALSE.
    ! generate a unique tag for the new cluster
    CLUSTP%MAXTAG = CLUSTP%MAXTAG + 1
    CLUSTP%IDTAG(CLUSTP%FIRSTEMPTY) = CLUSTP%MAXTAG
    
    FIRSTSET = .FALSE.
    DO CC = CNEW+1,CLUSTP%MAXCLUSTID
       IF (.NOT.CLUSTP%EXISTS(CC)) THEN
          CLUSTP%FIRSTEMPTY = CC
          FIRSTSET = .TRUE.
          EXIT
       ENDIF
    ENDDO
    IF (.NOT.FIRSTSET) THEN
	   CLUSTP%FIRSTEMPTY = CLUSTP%MAXCLUSTID + 1

       IF (CLUSTP%MAXCLUSTID.GT.CLUSTP%MAXNCLUST) THEN
          PRINT*, 'ERROR IN MAKENEWCLUSTER: you have exceeded the max allowed number of clusters!'
          STOP 1
       ENDIF
    ENDIF
    
    !PRINT*, 'CNEW,FIRSTEMPTY,MAXCLUSTID', CNEW, CLUSTP%FIRSTEMPTY, CLUSTP%MAXCLUSTID
    if (clustp%firstempty.LE.cnew) then
       PRINT*, 'ERROR'
       STOP 1
    ENDIF
  END SUBROUTINE MAKENEWCLUSTER


  SUBROUTINE INITIALIZECLUSTERSET(CLUSTP)
    USE KEYS, ONLY : KFISSION, KFMOV, KFSTAT, KPROD, UNITPROT, &
         & VEL, KDECAY, L, NSINKS, SINKWIDTH, TOTMITOINL
    IMPLICIT NONE
    TYPE(CLUSTERSET), POINTER :: CLUSTP

    IF (.NOT.CLUSTP%ARRAYSET) THEN
       PRINT*, 'ERROR: must set arrays before initializing'
       STOP 1
    ENDIF

    CLUSTP%KFMOV = KFMOV
    CLUSTP%KFSTAT = KFSTAT
    CLUSTP%KPROD = KPROD
    CLUSTP%UNITPROT = UNITPROT
    CLUSTP%VEL = VEL
    CLUSTP%KDECAY = KDECAY
    CLUSTP%L = L
    CLUSTP%SINKWIDTH = SINKWIDTH
    CLUSTP%NSINKS = NSINKS
    CLUSTP%NCLUST = TOTMITOINL
    CLUSTP%MAXTAG = 0 ! maximum tag id that has been used


  END SUBROUTINE INITIALIZECLUSTERSET
  
  SUBROUTINE SETUPCLUSTERSET(CLUSTP, MAXNCLUST, MAXCLUSTID)    
    IMPLICIT NONE
    TYPE(CLUSTERSET), POINTER :: CLUSTP
    INTEGER, INTENT(IN) :: MAXNCLUST, MAXCLUSTID

    IF (MAXCLUSTID.GT.MAXNCLUST) THEN
       PRINT*, 'ERROR IN SETUPCLUSTERSET: too many clusters', MAXNCLUST, MAXCLUSTID
       STOP 1
    ENDIF
    
    CLUSTP%MAXNCLUST = MAXNCLUST
    CLUSTP%MAXCLUSTID = MAXCLUSTID

    ALLOCATE(CLUSTP%POS(MAXNCLUST), CLUSTP%SIZE(MAXNCLUST), &
         & CLUSTP%STATE(MAXNCLUST), CLUSTP%PROT(MAXNCLUST), &
         & CLUSTP%PFORWARD(MAXNCLUST), CLUSTP%EXISTS(MAXNCLUST),&
         CLUSTP%DEAD(MAXNCLUST), CLUSTP%IDTAG(MAXNCLUST))

    CLUSTP%DEAD = .FALSE.
    CLUSTP%EXISTS = .FALSE.
    CLUSTP%EXISTS(1:MAXCLUSTID) = .TRUE.
    CLUSTP%FIRSTEMPTY = MAXCLUSTID+1
    
    
    CLUSTP%ARRAYSET = .TRUE.

    CLUSTP%POS = 0
    CLUSTP%STATE = 0
    CLUSTP%PROT = 0
    
  END SUBROUTINE SETUPCLUSTERSET

  SUBROUTINE CLEANUPCLUSTERSET(CLUSTP)
    IMPLICIT NONE
    TYPE(CLUSTERSET), POINTER :: CLUSTP
    
    DEALLOCATE(CLUSTP%POS, CLUSTP%SIZE, CLUSTP%STATE, CLUSTP%PROT, &
         & CLUSTP%PFORWARD,CLUSTP%EXISTS, CLUSTP%DEAD,CLUSTP%IDTAG)
  END SUBROUTINE CLEANUPCLUSTERSET
  
  SUBROUTINE HELLOWORLD(X)

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: X

    PRINT*, 'HELLO ', X
    
  END SUBROUTINE HELLOWORLD

  SUBROUTINE OUTPUTSNAPSHOT(CLUSTP,FILENAME,INFO,APPEND)
    ! Output a snapshot of current particle configuration
    ! ------------------
    ! input parameters:
    ! CLUSTP: structure describing mito clusters
    ! FILENAME: output file
    ! INFO: additional informational float (eg: time)
    ! APPEND: append to file or rewrite?
    ! -------------
    ! output format:
    ! 1st line -> number of clusters
    ! Successive lines -> one line per particle
    ! for each particle, lists CLUST ID,CLUSTP%POS,CLUSTP%PROT
    ! 1 to NSINKS ind of mito are sink mito

    IMPLICIT NONE

    TYPE(CLUSTERSET), POINTER :: CLUSTP
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    DOUBLE PRECISION, INTENT(IN) :: INFO(:)
    LOGICAL, INTENT(IN) :: APPEND
    INTEGER :: CC,NINFO,NLIVE
    !DOUBLE PRECISION :: ROTMAT(3,3)
    CHARACTER(LEN=2) :: X1
    CHARACTER(LEN=19) :: FMT1

    IF (APPEND) THEN
       OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN',POSITION='APPEND')
    ELSE
       OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN')
    END IF
    !NINFO = SIZE(INFO)
    !WRITE (X1,'(I2)') NINFO
    FMT1 = '(3ES10.2,3I4)'
    !WRITE(99,FMT1) NTRIALS,PARTP%NP, PARTP%NWALK, PARTP%NTETH, INFO
    NLIVE = COUNT(CLUSTP%EXISTS)
    
    WRITE(99,*) INFO, NLIVE !, INT(CLUSTP%NSINKS), CLUSTP%MAXCLUSTID

    DO CC = 1,CLUSTP%MAXCLUSTID
       ! write positions and protein content of particles
       IF (CLUSTP%EXISTS(CC)) THEN
          WRITE(99,*) CC, CLUSTP%POS(CC), CLUSTP%PROT(CC), CLUSTP%STATE(CC),CLUSTP%DEAD(CC), CLUSTP%IDTAG(CC)
       ENDIF
       
    END DO

    CLOSE(99)

  END SUBROUTINE OUTPUTSNAPSHOT

    SUBROUTINE OUTPUTFULLSNAPSHOT(CLUSTERLIST,FILENAME,INFO,APPEND)
    ! different format for outputing full snapshot info for many iterations
    ! ------------------
    ! input parameters:
    ! CLUSTERLIST: list of structures describing mito clusters
    ! FILENAME: output file
    ! INFO: additional informational float (eg: time)
    ! APPEND: append to file or rewrite?
    ! -------------
   

    USE KEYS, ONLY : MAXNCLUST
    IMPLICIT NONE

    TYPE(CLUSTERSET), INTENT(IN), TARGET :: CLUSTERLIST(:)
    TYPE(CLUSTERSET), POINTER :: CLUSTP
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    DOUBLE PRECISION, INTENT(IN) :: INFO(:)
    LOGICAL, INTENT(IN) :: APPEND
    INTEGER :: CC,NINFO,NLIVE, NCLUST, NITER, SC, LC
    !DOUBLE PRECISION :: ROTMAT(3,3)
    CHARACTER(LEN=2) :: X1
    CHARACTER(LEN=19) :: FMT1
    CHARACTER(LEN=MAXNCLUST) :: DEADSTRING

    IF (APPEND) THEN
       OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN',POSITION='APPEND')
    ELSE
       OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN')
    END IF


    NITER = SIZE(CLUSTERLIST)

    ! info for this snapshot: number of iterations, info (sim time)
    WRITE(99,*) NITER, INFO
    
    DO CC = 1,NITER
       CLUSTP=>CLUSTERLIST(CC)
       NCLUST = COUNT(CLUSTP%EXISTS)
       NLIVE = COUNT(CLUSTP%EXISTS.AND..NOT.CLUSTP%DEAD)

       ! info line: which iteration, current number of clusters, current number live, extra info (ie: snapshot number or time)       
       WRITE(99,*) CC, NCLUST,NLIVE,INFO

       IF (NCLUST.EQ.0) THEN
          ! just put a 0, avoid empty lines
          DO LC= 1,4
             WRITE(99,*) 0
          ENDDO
       ELSE
          ! output positions of existing clusters
          WRITE(99,*) PACK(CLUSTP%POS,CLUSTP%EXISTS)
          ! write protein contents
          WRITE(99,*) PACK(CLUSTP%PROT,CLUSTP%EXISTS)
          ! write state of clusters
          WRITE(99,*) PACK(CLUSTP%STATE,CLUSTP%EXISTS)

          ! write whether clusters are dead
          ! converting logicals to numbers
          WRITE(DEADSTRING,*) PACK(CLUSTP%DEAD,CLUSTP%EXISTS)
          DO SC = 1,len(DEADSTRING)
             IF (DEADSTRING(SC:SC).EQ.'T') THEN
                DEADSTRING(SC:SC)='1'
             ELSEIF (DEADSTRING(SC:SC).EQ.'F') THEN
                DEADSTRING(SC:SC) = '0'
             END IF
          ENDDO

          WRITE(99,*) TRIM(DEADSTRING)
       ENDIF
    ENDDO       

    CLOSE(99)

  END SUBROUTINE OUTPUTFULLSNAPSHOT
  
  
END MODULE MITOUTIL
