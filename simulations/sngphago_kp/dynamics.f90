MODULE DYNAMICS  
  USE MITOUTIL, ONLY : CLUSTERSET,OUTPUTSNAPSHOT, OUTPUTFULLSNAPSHOT
  USE MT19937, ONLY : GRND
  IMPLICIT NONE

CONTAINS
  SUBROUTINE RUNDYNAMICS(CLUSTERLIST, NSTEP, NITER, NREGIONS, &
       & MODELFLAG, DELT, OUTFILE, OUTPUTEVERY,SNAPSHOTFILE,&
       & SNAPSHOTEVERY, WHICHSNAPSHOT)
    USE KEYS, ONLY : PRINTEVERY, APPENDSNAPSHOTS,FILLUPREGIONS,&
         DEATHPROTLEVEL, VERBOSE, SM, AVGSM
    USE MITOUTIL, ONLY : MAKENEWCLUSTER
    IMPLICIT NONE
    TYPE(CLUSTERSET), POINTER :: CLUSTP
    TYPE(CLUSTERSET), INTENT(IN), TARGET :: CLUSTERLIST(NITER)
    ! SNAPSHOTFILE: file for dumping cluster position snapshots
    ! SNAPSHOTEVERY: how often to dump snapshots
    ! APPENDSNAPSHOTS: append snapshots to existing files (otherwise, start from scratch and append as we go)

    CHARACTER(LEN=*), INTENT(IN) :: SNAPSHOTFILE, OUTFILE
    INTEGER, INTENT(IN) :: SNAPSHOTEVERY, OUTPUTEVERY
    !    LOGICAL, INTENT(IN) :: APPENDSNAPSHOTS
    INTEGER, INTENT(IN) :: NSTEP ! number of steps to run
    DOUBLE PRECISION, INTENT(IN) :: DELT ! timestep size
    INTEGER, INTENT(IN) :: NITER ! number of iterations to run code
    INTEGER, INTENT(IN) :: NREGIONS ! number of demand regions in the domain
    INTEGER, INTENT(IN) :: WHICHSNAPSHOT ! which iteration to take snapshots of

    DOUBLE PRECISION :: CHECKCROSS, POSPREV(CLUSTERLIST(1)%MAXNCLUST)

    !INTEGER, POINTER :: WILLDELETE    
    ! store integer counting variables
    INTEGER :: STEP,CNEW, MAXN, CC,RC,I,RANDINT, NSINKS,SC,ITER, &
         & REACHEND, CROSSEVENTS, REMOVECLUST,NC,NCROSS, NINREG, &
         ATTEMPT, IC 
    DOUBLE PRECISION :: U, MAKEPROB, DECAYMULT, CURTIME, &
         & WALKPROB,SIGNVAL, MEANPROT


    DOUBLE PRECISION :: INFO(2)
    DOUBLE PRECISION :: RANDVAL
    LOGICAL :: MITOINSINK(CLUSTERLIST(1)%MAXNCLUST),&
         & MITOSTOP(CLUSTERLIST(1)%MAXNCLUST)
    INTEGER, INTENT(IN) :: MODELFLAG
    DOUBLE PRECISION :: SINKPOS,REGPOS
    DOUBLE PRECISION, DIMENSION(NREGIONS) :: POSNREGIONS,PROTNREGIONS,&
         & AVGPROTSINREG, VARPROTSINREG
    DOUBLE PRECISION, DIMENSION(NREGIONS,NITER) :: STOPPEDPROTS, PROTSINREG
    DOUBLE PRECISION, DIMENSION(NREGIONS,NITER) :: NREGIONPROTS
    DOUBLE PRECISION, DIMENSION(NREGIONS) :: AVGSTOPPEDPROT,VARSTOPPEDPROT
    DOUBLE PRECISION :: AVGTOTPROT, AVGCLUST, LEAVEPROT,VARTOTPROT,AVGSTOPPED
    DOUBLE PRECISION, DIMENSION(NITER) :: TOTPROT, TOTCLUST,STOPPEDMITO
    INTEGER :: CTLEAVEPROT
    INTEGER :: SCCOUNT, SCCHOOSE
    INTEGER, DIMENSION(NITER, CLUSTERLIST(1)%MAXNCLUST) :: SINKREGNO
    !INTEGER, DIMENSION(NITER, NREGIONS) ::  REGIONNEEDSMITOS
    LOGICAL, DIMENSION(CLUSTERLIST(1)%MAXNCLUST) :: HAVERESOLVED
    INTEGER, DIMENSION(CLUSTERLIST(1)%MAXNCLUST) :: DIDCROSS,INTHISREGION
    DOUBLE PRECISION, DIMENSION(CLUSTERLIST(1)%MAXNCLUST) :: PREVDIST
    DOUBLE PRECISION :: TPROD, TMPTIME, VARTOTSTOPPROT


    STOPPEDPROTS = 0D0
    PROTSINREG = 0D0
    NREGIONPROTS = 0D0
    PROTNREGIONS = 0D0
    !REGIONNEEDSMITOS = 0
    AVGCLUST = 0
    AVGSTOPPEDPROT = 0
    AVGTOTPROT = 0

    CLUSTP=>CLUSTERLIST(1)
    SIGNVAL = 1.0
 
    MAKEPROB = 1D0-EXP(-CLUSTP%KPROD*DELT)
    !PRINT*, 'KPROD=',CLUSTP%KPROD, 'MAKEPROB', MAKEPROB
    !WRITE(*,*) "makeprob", MAKEPROB
    DECAYMULT = EXP(-CLUSTP%KDECAY*DELT)

    ! probability of restarting
    WALKPROB = 1D0-EXP(-CLUSTP%KW*DELT)

    print*, 'MODELFLAG = ', MODELFLAG

    CURTIME = 0D0

    ! declare array
    !TYPE(CLUSTERSET), TARGET :: CLUSTERLIST(:)

    ! Initial snapshot
    INFO =(/ CURTIME, NREGIONS*1d0 /)


    ! position nregions in the domain
    DO RC = 1,NREGIONS
       POSNREGIONS(RC) = (CLUSTP%L/(NREGIONS+1)) * RC
       ! PROBVALREG(RC) =  EXP(RC) / (EXP(1) * (1-EXP(RC))/(1 - EXP(1)))
    END DO


    ! position stationary mitochondria
    ! in sngphago, have an initial number of stationary mitochondria in the domain,
    ! similar to ssphago
    ! have similar implementation of -ve kp, to also have motile mito 
    ! not sure what to do with sinkregno 
    DO ITER = 1,NITER
       CLUSTP=>CLUSTERLIST(ITER)
       NSINKS = INT(AVGSM*NREGIONS) !set NSINKS to the expected stationary fraction
       DO CC = 1,NSINKS
          CALL MAKENEWCLUSTER(CLUSTP,CNEW)
          ! if Nregions is not provided, position sinks uniformly
          IF (NREGIONS .LT. 0) THEN
             CLUSTP%POS(CNEW) = (CLUSTP%L/(NSINKS+1))* CNEW
             SINKREGNO(ITER,CNEW) = CNEW 
          ELSE IF (FILLUPREGIONS) THEN
             ! fill up each region (in order, SEVERAL TIMES AROUND), 
             ! then distribute remaining sinks randomly
             IF (CC.LE.NREGIONS*(NSINKS/NREGIONS)) THEN
                CLUSTP%POS(CNEW) = POSNREGIONS(MOD(CC-1,NREGIONS)+1)
                SINKREGNO(ITER,CNEW) = MOD(CC-1,NREGIONS)+1 ! to a separate array, assign which region it is confined to
             ELSE
                RANDVAL = GRND()
                RANDINT = 1 + FLOOR(RANDVAL * NREGIONS) ! choose an index RC
                CLUSTP%POS(CNEW) = POSNREGIONS(RANDINT)
                SINKREGNO(ITER,CNEW) = RANDINT 
             ENDIF
          ELSE             
             ! position all sinks randomly across domains
             RANDVAL = GRND()
             RANDINT = 1 + FLOOR(RANDVAL * NREGIONS) ! choose an index RC
             CLUSTP%POS(CNEW) = POSNREGIONS(RANDINT)
             SINKREGNO(ITER,CNEW) = RANDINT 
          END IF

          ! position them all in one region
          ! CLUSTP%POS(CNEW) = POSNREGIONS(1) ! all at the start
          ! CLUSTP%POS(CNEW) = POSNREGIONS(NREGIONS) ! all at end
          !CLUSTP%POS(CNEW) = POSNREGIONS(INT(NREGIONS/2)) ! all in middle

          ! no doubling up
          ! CLUSTP%POS(CNEW) = POSNREGIONS(CC)



          ! position them in specific indices
          ! CLUSTP%POS(CNEW) = POSNREGIONS(CNEW * &
          !     & INT(NREGIONS/MINVAL(CLUSTP%NSINKS,NREGIONS)))
          ! CLUSTP%POS(CNEW) = (CLUSTP%L/(NSINKS+1))* CNEW !position when &
          ! & ! nregions is not defined
          !CLUSTP%POS(1:NSINKS) = (CLUSTP%L/(NSINKS+1) * (/ (I, I = 1,NSINKS,1/) )
          CLUSTP%STATE(CNEW) = 0 ! permanently stopped

          ! sample starting protein to generate uniformly distributed death times
          CLUSTP%PROT(CNEW) = DEATHPROTLEVEL*EXP(GRND()*EXP(-CLUSTP%KDECAY*DELT))
          CLUSTP%PROT(CNEW) = grnd()*clustp%unitprot
          !          CLUSTP%PROT(CNEW) = GRND()*(CLUSTP%UNITPROT-DEATHPROTLEVEL) + DEATHPROTLEVEL
          CLUSTP%SIZE(CNEW) = CLUSTP%MITOSIZE ! size is the same as the sink size, not relevant for now          
          !          PRINT*, 'CNEW,FIRSTEMPTY,MAXCLUSTID', CNEW, CLUSTP%FIRSTEMPTY, CLUSTP%MAXCLUSTID
       ENDDO

       ! initialize moving mito, if kp is set -ve
       IF (CLUSTP%KPROD .LT. 0) THEN 
          DO CC = NSINKS+1,CLUSTP%NCLUST
             CALL MAKENEWCLUSTER(CLUSTP,CNEW)
             CLUSTP%POS(CC) = CLUSTP%L*GRND()
             U = GRND()
             IF (U.LT.0.5D0) THEN
                CLUSTP%STATE(CC) = -1
             ELSE
                CLUSTP%STATE(CC) = 1
             ENDIF

             ! sample starting protein to generate uniformly distributed death times
             CLUSTP%PROT(CC) = DEATHPROTLEVEL*EXP(GRND()*EXP(-CLUSTP%KDECAY*DELT))
             clustp%prot(cc) = grnd()*clustp%unitprot
             ! CLUSTP%PROT(CC) =  GRND()*(CLUSTP%UNITPROT-DEATHPROTLEVEL) + DEATHPROTLEVEL
          ENDDO
       ENDIF

       !PRINT*, 'TESTX1:', ITER, CLUSTP%PROT(1:CLUSTP%NCLUST)
    END DO

    ! Counting cluster reaching ends and crossing events
    REACHEND = 0
    CROSSEVENTS = 0
    REMOVECLUST = 0

    IF (WHICHSNAPSHOT.LE.0) THEN
       CALL OUTPUTFULLSNAPSHOT(CLUSTERLIST,SNAPSHOTFILE,INFO,APPENDSNAPSHOTS)
    ELSE
       CLUSTP=> CLUSTERLIST(WHICHSNAPSHOT)
       CALL OUTPUTSNAPSHOT(CLUSTP,SNAPSHOTFILE,INFO,APPENDSNAPSHOTS)
    END IF

    ! output sink protein content
    ! will output average over all iterations
    ! they start identical, so just output the values for the first iter initially

    OPEN(UNIT=88,FILE=OUTFILE,STATUS='UNKNOWN')
    WRITE(88,*) NSINKS,CLUSTP%KPROD,CLUSTP%PS, &
         & CLUSTP%KW, CLUSTP%VEL, NREGIONS, CLUSTP%KDECAY
    STEP = 0
    WRITE(88,*) STEP, CURTIME, CLUSTP%PROT(1:NSINKS), CLUSTP%PROT(1:NSINKS)

    LEAVEPROT = 0; CTLEAVEPROT = 0

    IF (VERBOSE) THEN
       DO CC = 1,CLUSTP%NCLUST
          PRINT*, 'INIT POS, PROT:', CC, CLUSTP%POS(CC), CLUSTP%PROT(CC)
       ENDDO
    ENDIF

    DO STEP = 1,NSTEP       
       TOTPROT = 0D0
       TOTCLUST = 0
       DO ITER = 1,NITER

          CLUSTP=> CLUSTERLIST(ITER)

          ! IF (STEP == 1) THEN
          !    PRINT*, 'CLUSTP%FIRSTEMPTY',CLUSTP%FIRSTEMPTY
          ! END IF

          IF (CLUSTP%FIRSTEMPTY.GT.CLUSTP%MAXCLUSTID+1) THEN
             print*, 'ERROR IN DYNAMICS: firstempty is bigger than maxclustid+1'
             PRINT*, 'TESTX2:', CLUSTP%FIRSTEMPTY, CLUSTP%MAXCLUSTID, ITER, CLUSTERLIST(ITER)%FIRSTEMPTY
             PRINT*, 'TESTX3:', CLUSTP%FIRSTEMPTY, CLUSTP%MAXCLUSTID, CLUSTP%EXISTS(9:15), REMOVECLUST
             STOP 1
          ENDIF


          ! order of the operations set to ensure no cluster is
          ! outside the domain


          ! ------------- TURN OFF IF NOT USING KP ---------------
          ! decide whether to produce a mito.
          ! Can produce multiple ones in a timestep.
          ! TUNR OFF FOR -VE KPROD
          IF (CLUSTP%KPROD .GE. 0) THEN
             TMPTIME = 0D0 ! time into the step
             DO WHILE (TMPTIME.LT.DELT)
                ! time to produce next mito
                U = GRND()
                TPROD = -1/CLUSTP%KPROD*LOG(U)
                TMPTIME = TMPTIME+TPROD
                !PRINT*, 'TESTX1:', COUNT(CLUSTP%EXISTS), TPROD, TMPTIME
                IF (TMPTIME.LT.DELT) THEN ! next mito gets made before the step is over
                   CALL MAKENEWCLUSTER(CLUSTP,CNEW)
                   ! where does new mito start (will move v*delt by end of step)
                   CLUSTP%POS(CNEW) = -TMPTIME*CLUSTP%VEL
                   CLUSTP%STATE(CNEW) = 1
                   CLUSTP%SIZE(CNEW) = 1
                   ! starting protein content.
                   ! If produced close to end of step, will have content 1 by end of step
                   CLUSTP%PROT(CNEW) = CLUSTP%UNITPROT*EXP(CLUSTP%KDECAY*TMPTIME)
                ENDIF
             ENDDO
          ENDIF
          ! --------------------------------------------

          ! degrade proteins within the cluster
          CLUSTP%PROT(1:CLUSTP%MAXCLUSTID) = CLUSTP%PROT(1:CLUSTP%MAXCLUSTID) * DECAYMULT

          ! step all clusters forward in time
          ! current implementation : 5 states (-2,-1,0,1,2), -2 and 2 are not implemented currently

          ! first save previous position
          POSPREV = CLUSTP%POS

          CLUSTP%POS(1:CLUSTP%MAXCLUSTID) = CLUSTP%POS(1:CLUSTP%MAXCLUSTID) &
               & +CLUSTP%VEL*(MOD(CLUSTP%STATE(1:CLUSTP%MAXCLUSTID),2))*DELT


          ! reflect clusters off boundary

          DO CC = 1,CLUSTP%MAXCLUSTID

             IF (CLUSTP%EXISTS(CC) .AND. CLUSTP%POS(CC).GT. CLUSTP%L) THEN !CHECK IF RIGHT EDGE IS OUT OF DOMAIN
                IF (ITER.EQ.1) REACHEND = REACHEND + 1
                CLUSTP%POS(CC) = CLUSTP%L - (CLUSTP%POS(CC) - CLUSTP%L) ! REFLECT OFF CLUSTER
                CLUSTP%STATE(CC) = -1 ! PAUSE AND WALK BACKWARDS NEXT
             END IF

          END DO

          ! ----- CHANGE BELOW IF NOT USING KP --------------
          IF (CLUSTP%KPROD .GE. 0) THEN 

             ! absorb clusters which left the domain
             MAXN = CLUSTP%MAXCLUSTID 
             CLUSTP%MAXCLUSTID = 1
             DO CC = 1,MAXN

                IF (CLUSTP%EXISTS(CC) .AND. CLUSTP%POS(CC).LT.0) THEN !CHECK IF CENTER IS OUT OF DOMAIN
                   ! keep track of avg protein content leaving the domain
                   LEAVEPROT = LEAVEPROT + CLUSTP%PROT(CC)
                   CTLEAVEPROT = CTLEAVEPROT + 1
                   CLUSTP%EXISTS(CC) = .FALSE.  ! DELETE THE CLUSTER

                   ! IF a cluster gets deleted, move the index of first empty position to that cluster's index
                   IF (CLUSTP%FIRSTEMPTY .GT. CC) THEN
                      CLUSTP%FIRSTEMPTY = CC ! SET THE FIRST OF SUCH AS FIRSTEMPTY
                      ! ADD A STEP WHICH SHORTENS MAXCLUSTID?
                      ! OR HAVE A DIFF SUBROUTINE TO DELETE CLUSTER
                   END IF

                END IF

                IF (CLUSTP%EXISTS(CC)) CLUSTP%MAXCLUSTID = CC

             END DO

          ELSE ! implement constant Nmito case
             DO CC = 1,CLUSTP%NCLUST

                IF (CLUSTP%EXISTS(CC) .AND. CLUSTP%POS(CC).LT.0) THEN !CHECK IF CENTER IS OUT OF DOMAIN
                   ! keep track of avg protein content leaving the domain RETAINING THE OLD STRUCTURE FOR NOW
                   LEAVEPROT = LEAVEPROT + CLUSTP%PROT(CC)
                   CTLEAVEPROT = CTLEAVEPROT + 1

                   CLUSTP%POS(CC) =  - CLUSTP%POS(CC) ! REFLECT OFF CLUSTER
                   CLUSTP%STATE(CC) = 1 ! MOVE CLUSTER ANTEROGRADE AGAIN
                   CLUSTP%DEAD(CC) = .FALSE.  ! REVIVE THE CLUSTER IF IT WAS DEAD
                   CLUSTP%PROT(CC) = CLUSTP%UNITPROT ! PROTEIN REPLENISHES
                   IF (VERBOSE) print*, ITER, STEP, 'RESURRECTED: Cluster ', CC
                END IF

             END DO
          END IF

          ! -----------------------------------------------------
          ! PRINT*, 'TESTX1:', ITER, STEP, CLUSTP%PROT(8:9), CLUSTP%PROT(11), &
          !      & SQRT(SUM(CLUSTP%PROT(1:10)**2)/10 - (SUM(CLUSTP%PROT(1:10))/10)**2)


          ! check if any have passed below a cutoff and should be dead (phagocytosed)          
          DO CC = 1,CLUSTP%MAXCLUSTID

             ! IF (SINKREGNO(ITER,CC).EQ.NREGIONs.AND.ITER.EQ.1) THEN
             !    PRINT*, 'last reg prot levels:', STEP, CC, SINKREGNO(ITER,CC), CLUSTP%PROT(CC)
             ! ENDIF

             IF (CLUSTP%EXISTS(CC).AND.CLUSTP%PROT(CC).LT.DEATHPROTLEVEL.AND..NOT.CLUSTP%DEAD(CC)) THEN
                CLUSTP%DEAD(CC) = .TRUE.
                IF (VERBOSE) PRINT*, ITER, STEP, 'DIED: Cluster ', CC, 'in region ', SINKREGNO(ITER,CC), &
                     & clustp%prot(cc), CLUSTP%POS(CC)

                IF (CLUSTP%STATE(CC).EQ.0) THEN
                   ! IF (ITER.EQ.1) THEN
                   !    PRINT*, 'stationary mito just died: step, iter, cc, reg no, prot', STEP, ITER, CC, &
                   !    SINKREGNO(ITER,CC),CLUSTP%PROT(CC)
                   ! ENDIF

                   ! this region will need to get a replacement mitochondrion
                   ! REMOVE THIS LINE LATER, NOT NEEDED FOR SNGPHAGO
                   !REGIONNEEDSMITOS(ITER,SINKREGNO(ITER,CC)) = REGIONNEEDSMITOS(ITER,SINKREGNO(ITER,CC)) + 1
                ENDIF

                CLUSTP%STATE(CC) = -1
                SINKREGNO(ITER,CC) = 0
             ENDIF
          END DO





          ! --------------
          ! implement stopngo dynamics
          ! --------------
          ! cycle through all clusters; if stopped mito should start
          ! or if moving mito crossed region, they could stop
          DO CC = 1,CLUSTP%MAXCLUSTID
             IF (.NOT.CLUSTP%EXISTS(CC).OR.CLUSTP%DEAD(CC)) CYCLE
             IF (CLUSTP%STATE(CC) .EQ. 0) THEN
                ! does this stopped cluster restart
                U = GRND()
                IF (U .LT. WALKPROB) THEN
                   RANDVAL = GRND()
                   IF (RANDVAL .LT. 0.5) THEN
                      CLUSTP%STATE(CC) = -1
                   ELSE
                      CLUSTP%STATE(CC) = 1
                   END IF
                END IF
             ELSE
                ! if walking, resolve if it crossed a region
                DO RC = 1,NREGIONS
                   CHECKCROSS = 0
                   PREVDIST(CC) = (POSPREV(CC) - POSNREGIONS(RC))
                   CHECKCROSS = PREVDIST(CC) *&
                        & (CLUSTP%POS(CC) - POSNREGIONS(RC))
                   IF (CHECKCROSS .LT. 0) THEN
                      U = GRND()
                      IF (U .LT. CLUSTP%PS) THEN
                         CLUSTP%STATE(CC) = 0
                         CLUSTP%POS(CC) = POSNREGIONS(RC) ! update to position of region crossed
                      END IF
                   END IF
                   IF (CHECKCROSS .LT. 0) EXIT !exit loop since no more resolution required
                END DO
             END IF
          END DO






!!$             
!!$          DO RC = 1,NREGIONS
!!$             NCROSS = 0
!!$             NINREG = 0
!!$             DO CC = 1,CLUSTP%MAXCLUSTID                
!!$                IF (.NOT.CLUSTP%EXISTS(CC).OR.CLUSTP%DEAD(CC)) CYCLE
!!$
!!$                IF (CLUSTP%STATE(CC).NE.0) THEN
!!$                   ! did this cluster cross this region on this step?
!!$                   PREVDIST(CC) = (POSPREV(CC) - POSNREGIONS(RC))
!!$                   CHECKCROSS = PREVDIST(CC) *&
!!$                        & (CLUSTP%POS(CC) - POSNREGIONS(RC))
!!$                   IF (CHECKCROSS.LT.0) THEN
!!$                      NCROSS = NCROSS+1
!!$                      DIDCROSS(NCROSS) = CC                  
!!$                   ENDIF
!!$                ENDIF
!!$                
!!$                ! keep track of which clusters are in this region
!!$                IF (CLUSTP%STATE(CC).EQ.0.AND.SINKREGNO(ITER,CC).EQ.RC) THEN
!!$                   NINREG=NINREG+1
!!$                   INTHISREGION(NINREG) = CC
!!$                ENDIF
!!$             ENDDO
!!$
!!$             ! sort crossing mitochondria by position
!!$             ! WARNING: THIS PART IS CURRENTLY VERY INEFFICIENT. REPLACE WITH PROPER SORTING SOMETIME...
!!$             HAVERESOLVED = .FALSE.
!!$             DO ATTEMPT = 1,NCROSS ! resolve the next crossing mitochondria
!!$                ! find the next one to cross
!!$                IC = MINLOC(ABS(PREVDIST(DIDCROSS(1:NCROSS))),DIM=1,MASK = (.NOT.HAVERESOLVED(1:NCROSS)))
!!$                NC = DIDCROSS(IC)
!!$
!!$                ! decide whether to fuse with existing mitochondria
!!$                IF (MODELFLAG.EQ.0) THEN ! option to fuse with all mito in sink
!!$                   DO SC = 1,NINREG ! go through all mitos in this region
!!$                      RANDVAL = GRND()
!!$                      IF (RANDVAL.LT.FUSPROBSTAT) THEN                         
!!$                         !IF (ITER.EQ.1) CROSSEVENTS = CROSSEVENTS + 1
!!$                         IF (VERBOSE) print*, iter, step, 'FUSE: region ', rc, 'clusters ', nc, &
!!$                              & inthisregion(sc), CLUSTP%PROT(NC), CLUSTP%PROT(INTHISREGION(SC))!clustp%state(nc), clustp%state(inthisregion(sc)) 
!!$                         MEANPROT = (CLUSTP%PROT(NC) + CLUSTP%PROT(INTHISREGION(SC))) /2
!!$                         CLUSTP%PROT(NC) = MEANPROT
!!$                         CLUSTP%PROT(INTHISREGION(SC)) = MEANPROT                         
!!$                      END IF
!!$                   ENDDO
!!$                ELSEIF (MODELFLAG.EQ.1) THEN
!!$                   ! decide whether to fuse in this region
!!$                   RANDVAL = GRND()
!!$                   IF (RANDVAL.LT.FUSPROBSTAT .AND. NINREG .GT. 0) THEN                      
!!$                      ! pick which mito to fuse with
!!$                      SC = FLOOR(GRND()*NINREG)+1
!!$                      IF (SC.GT.NINREG) THEN
!!$                         PRINT*, 'ERROR: BAD CHOICE OF FUSION MITO', NINREG
!!$                         STOP 1
!!$                      ENDIF
!!$                      MEANPROT = (CLUSTP%PROT(NC) + CLUSTP%PROT(INTHISREGION(SC))) /2
!!$                      CLUSTP%PROT(NC) = MEANPROT
!!$                      CLUSTP%PROT(INTHISREGION(SC)) = MEANPROT
!!$                   ENDIF
!!$                ELSE
!!$                   PRINT*, 'ERROR: BAD MODEL FLAG'
!!$                   STOP 1
!!$                ENDIF
!!$
!!$                ! decide whether this mito should stop
!!$                IF (REGIONNEEDSMITOS(ITER,RC).GT.0) THEN
!!$                   CLUSTP%STATE(NC) = 0
!!$                   REGIONNEEDSMITOS(ITER,RC) = REGIONNEEDSMITOS(ITER,RC)-1
!!$                   SINKREGNO(ITER,NC) = RC ! set which sink this stationary mito is at
!!$                   IF (VERBOSE) PRINT*, ITER, STEP, 'STOPPED: cluster ', NC, ' region ', rc, CLUSTP%PROT(NC)
!!$                ENDIF
!!$
!!$                HAVERESOLVED(IC) = .TRUE.
!!$             ENDDO
!!$          ENDDO
!!$
!!$
!!$
          ! exchange proteins with sink mitochondria AND replace missing mitochondria
          ! if MODELFLAG = 0, then exchange possible with all mito in sink
          ! if MODELFLAG = 1, then only 1 exchange at a time possible
          !IF (MODELFLAG .EQ. 0) THEN ! go ahead with old algorithm
          !   DO CC = 1,CLUSTP%MAXCLUSTID
          !      IF (.NOT.CLUSTP%EXISTS(CC) .OR. CLUSTP%DEAD(CC)) CYCLE
          !      
          !      DO SC = 1,CLUSTP%MAXCLUSTID
          !         IF (CLUSTP%STATE(SC) .NE. 0 .OR. &
          !              & (.NOT. CLUSTP%EXISTS(SC)) .OR. &
          !              & CLUSTP%DEAD(SC) ) CYCLE !only check with stopped  
          !         CHECKCROSS = (POSPREV(CC) - POSPREV(SC)) *&
          !              & (CLUSTP%POS(CC) - CLUSTP%POS(SC))

          !         IF (CHECKCROSS .LT. 0) THEN
          ! if there has been an intersection, prot equilibration occurs
          !            RANDVAL = GRND()
          !            IF (RANDVAL.LT.FUSPROBSTAT) THEN
          !               IF (ITER.EQ.1) CROSSEVENTS = CROSSEVENTS + 1
          !               MEANPROT = (CLUSTP%PROT(CC) + CLUSTP%PROT(SC)) /2
          !               CLUSTP%PROT(CC) = MEANPROT
          !               CLUSTP%PROT(SC) = MEANPROT
          !            END IF
          !         END IF
          !print*,'cluster and sink no', CC,SC
          !      END DO
          !   END DO


          !ELSE IF (MODELFLAG .EQ. 1) THEN ! new algorithm             
          !  DO CC = 1,CLUSTP%MAXCLUSTID
          !    IF (.NOT.CLUSTP%EXISTS(CC) .OR. CLUSTP%DEAD(CC)) CYCLE
          !    DO RC = 1, NREGIONS
          !       REGPOS = POSNREGIONS(RC)
          !       CHECKCROSS = (POSPREV(CC) - POSNREGIONS(RC)) &
          !            & * (CLUSTP%POS(CC) - POSNREGIONS(RC))
          !       IF (CHECKCROSS .LT. 0) THEN ! check if the cluster has intersected with the region
          !          RANDVAL = GRND()
          !          IF (RANDVAL .LT. FUSPROBSTAT) THEN
          !             IF (ITER .EQ. 1) CROSSEVENTS = CROSSEVENTS + 1
          !             ! choose which of the mito in that region it fuses with
          !             SCCOUNT = 0
          !             DO SC = 1,CLUSTP%MAXCLUSTID
          !                IF (CLUSTP%STATE(SC) .NE. 0 .OR. CLUSTP%DEAD(SC) &
          !                     & .OR. .NOT.(CLUSTP%EXISTS(SC))) CYCLE
          !                IF (CLUSTP%POS(SC).GT.REGPOS-CLUSTP%SINKWIDTH/2 &
          !                     & .AND.CLUSTP%POS(SC).LT.REGPOS+CLUSTP%SINKWIDTH/2) THEN
          !                   SCCOUNT = SCCOUNT + 1
          !                   RANDVAL = GRND()
          !                   RANDINT = 1 + FLOOR(SCCOUNT * RANDVAL)
          !                   IF (RANDINT .EQ. SCCOUNT) SCCHOOSE = SC
          !                END IF
          !             END DO
          !             MEANPROT = (CLUSTP%PROT(CC) + CLUSTP%PROT(SCCHOOSE)) /2
          !             CLUSTP%PROT(CC) = MEANPROT
          !             CLUSTP%PROT(SCCHOOSE) = MEANPROT
          !          END IF

          !       END IF

          !    END DO
          ! END DO

          !END IF






          ! exchange with moving mito and moving mito
          !DO CC = NSINKS+1, CLUSTP%MAXCLUSTID

          !  DO SC = NSINKS+1, CC
          !     CHECKCROSS = (POSPREV(CC) - POSPREV(SC)) * (CLUSTP%POS(CC) - CLUSTP%POS(SC))

          !    IF (CHECKCROSS .LT. 0) THEN
          ! if there has been an intersection, prot equilibration occurs
          !      RANDVAL = GRND()
          !print*, 'Randval exchange',RANDVAL
          !     IF (RANDVAL.LT.FUSPROBMOV) THEN
          !       IF (ITER.EQ.1) CROSSEVENTS = CROSSEVENTS + 1
          !       MEANPROT = (CLUSTP%PROT(CC) + CLUSTP%PROT(SC)) /2
          !      CLUSTP%PROT(CC) = MEANPROT
          !     CLUSTP%PROT(SC) = MEANPROT
          !  END IF
          ! END IF

          !  END DO
          !  END DO

          IF (MOD(STEP,OUTPUTEVERY).EQ.0) THEN
             ! get average Nregion  content for this iteration

             ! compile total protein content in Nregions
             DO RC = 1,NREGIONS
                REGPOS = POSNREGIONS(RC)
                MITOINSINK = CLUSTP%EXISTS .AND. &
                     & CLUSTP%POS.GT.REGPOS-CLUSTP%SINKWIDTH/2 &
                     & .AND.CLUSTP%POS.LT.REGPOS+CLUSTP%SINKWIDTH/2
                ! this is total amount of proteins within a certain width
                ! of each region center (walking and not)
                NREGIONPROTS(RC,ITER) = SUM(CLUSTP%PROT,MASK=MITOINSINK)
                ! get avg stopped content
                MITOSTOP = CLUSTP%EXISTS .AND. (.NOT.CLUSTP%DEAD) &
                     & .AND.CLUSTP%POS.GT.REGPOS-CLUSTP%SINKWIDTH/2 &
                     & .AND.CLUSTP%POS.LT.REGPOS+CLUSTP%SINKWIDTH/2 &
                     & .AND. CLUSTP%STATE .EQ. 0
                STOPPEDPROTS(RC, ITER) = SUM(CLUSTP%PROT, MASK = MITOSTOP)
             END DO
          ENDIF


          ! Track total protein content
          TOTPROT(ITER) =  SUM(CLUSTP%PROT,MASK=CLUSTP%EXISTS.and..NOT.CLUSTP%DEAD)
          TOTCLUST(ITER) = COUNT(CLUSTP%EXISTS)
          STOPPEDMITO(ITER) = COUNT(CLUSTP%EXISTS .AND. CLUSTP%STATE .EQ. 0)
       ENDDO

       ! write to snapshot file (first iter only)
       IF (MOD(STEP,SNAPSHOTEVERY).EQ.0) THEN
          INFO(1) = CURTIME
          IF (WHICHSNAPSHOT.LE.0) THEN ! output all iterations
             CALL OUTPUTFULLSNAPSHOT(CLUSTERLIST,SNAPSHOTFILE,INFO,.TRUE.)
          ELSE
             ! output single iteration snapshot
             CLUSTP=>CLUSTERLIST(WHICHSNAPSHOT)
             INFO(1) = CURTIME
             CALL OUTPUTSNAPSHOT(CLUSTP,SNAPSHOTFILE,INFO,.TRUE.)
          ENDIF
       END IF

       IF (MOD(STEP,OUTPUTEVERY).EQ.0) THEN
          AVGSTOPPEDPROT =SUM(STOPPEDPROTS,2)/NITER
          AVGPROTSINREG = SUM(NREGIONPROTS,2)/NITER ! prots in each region
          VARSTOPPEDPROT = SUM(STOPPEDPROTS**2,2)/NITER - AVGSTOPPEDPROT**2 
          VARPROTSINREG = SUM(NREGIONPROTS**2,2)/NITER - AVGPROTSINREG**2

          AVGTOTPROT = SUM(TOTPROT)/NITER
          VARTOTPROT = SUM(TOTPROT**2)/NITER - AVGTOTPROT**2
          ! variance in total stopped protein 
          VARTOTSTOPPROT = SUM(SUM(STOPPEDPROTS,1)**2)/NITER - SUM(AVGSTOPPEDPROT)**2

          AVGCLUST = SUM(TOTCLUST)/NITER
          AVGSTOPPED = SUM(STOPPEDMITO)/NITER


          ! output data          
          WRITE(88,*) STEP, CURTIME, AVGPROTSINREG,VARPROTSINREG,&
               & AVGSTOPPEDPROT,VARSTOPPEDPROT, AVGTOTPROT,VARTOTPROT, &
               & AVGCLUST, AVGSTOPPED, LEAVEPROT/MAX(CTLEAVEPROT,1),VARTOTSTOPPROT

          !          PRINT*, 'TESTX1:', SUM(STOPPEDPROTS,1)
       ENDIF


       CURTIME = CURTIME + DELT

       IF (MOD(STEP,PRINTEVERY).EQ.0) THEN          
          PRINT*, 'STEP,time,VEL, number clusters, tot prot, stopped prot:', STEP,CURTIME, &
               & CLUSTP%VEL, AVGCLUST, AVGTOTPROT, SUM(AVGSTOPPEDPROT)!LEAVEPROT/MAX(CTLEAVEPROT,1)
          ! PRINT*, 'TESTX1:', CLUSTP%PROT(1:CLUSTP%NCLUST)
          !  PRINT*, 'TESTX2:', CLUSTP%Pos(1:CLUSTP%NCLUST)
       ENDIF


    END DO

    !PRINT*, 'Reflection events, crossing events:', REACHEND,CROSSEVENTS




  END SUBROUTINE RUNDYNAMICS
END MODULE DYNAMICS




! reset stopped or walking state
! OLD CODE: Continuous stopping and walking of mitochondria
!       MAXN = CLUSTP%MAXCLUSTID
!       DO CC = 1,CLUSTP%MAXCLUSTID
!          IF(CLUSTP%EXISTS(CC)) THEN
!             ! TEST IF CLUSTER IS STOPPED
!             IF (CLUSTP%STATE(CC) == 0) THEN
!                U = GRND()
!                IF (U.LT.WALKPROB) THEN
!                    CLUSTP%STATE(CC) = SIGN(SIGNVAL,GRND()-1) ! 1 IF MOVING FORWARD, -1 IF MOVING BACKWARDS
!!OLD CODE: USES PFORWARD (COMMENTED)
!!CLUSTP%STATE(CC) = SIGN(SIGNVAL,CLUSTP%PFORWARD(CC) - GRND()) ! SET TO WALKING
!!CLUSTP%PFORWARD(CC) = (CLUSTP%STATE(CC) + 1)/2

!ENDIF
!ELSE !IF WALKING STATE
!U = GRND()
!IF (U.LT.PAUSEPROB) THEN
!CLUSTP%STATE(CC) = 0 !SET TO STOPPED
!END IF
!END IF
!END IF
!END DO





