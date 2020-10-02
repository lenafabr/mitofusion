PROGRAM MAIN
  USE KEYS

  IMPLICIT NONE

  CALL READKEY

  SELECT CASE(ACTION)
  CASE('RUNDYNAMICS')
     CALL DYNAMICDRIVER
  CASE DEFAULT
     PRINT*, 'I do not understand this action: ', ACTION
     STOP 1
  END SELECT

CONTAINS
  SUBROUTINE DYNAMICDRIVER
    USE DYNAMICS, ONLY : RUNDYNAMICS
    USE MITOUTIL, ONLY : CLUSTERSET, SETUPCLUSTERSET, INITIALIZECLUSTERSET, CLEANUPCLUSTERSET

    IMPLICIT NONE

    TYPE(CLUSTERSET), TARGET :: OURCLUSTERS(NITER)
    TYPE(CLUSTERSET), POINTER :: CLUSTP
    INTEGER :: ITER

    DO ITER = 1,NITER
       CLUSTP=>OURCLUSTERS(ITER)    
       CALL SETUPCLUSTERSET(CLUSTP, MAXNCLUST, 0)
       CALL INITIALIZECLUSTERSET(CLUSTP)
    ENDDO

    CALL RUNDYNAMICS(OURCLUSTERS, NSTEP,NITER,NREGIONS, MODELFLAG, &
         & DELT, OUTFILE, OUTPUTEVERY, SNAPSHOTFILE, SNAPSHOTEVERY,&
         & whichsnapshot)

    CALL CLEANUPCLUSTERSET(CLUSTP)
  END SUBROUTINE DYNAMICDRIVER
END PROGRAM MAIN