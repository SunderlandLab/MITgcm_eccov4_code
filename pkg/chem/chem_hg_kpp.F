C $Header: /u/gcmpack/MITgcm/pkg/chem/chem_hg_kpp.F,v 1.10 2014/03/26
C 00:14:47 jmc Exp $
C $Name: checkpoint64g $

C To calculate the chemistry of Hg species in the ocean using KPP solver
C (finally :) ), yxzhang,
C 03/26/14
C

CStartOfInterFace
      SUBROUTINE CHEM_HG_KPP( DeltaT, CONC )

C     /==========================================================\
C     | SUBROUTINE CHEM_HG_KPP                                  |
C     |   Wrapper for other KPP files                            |
C     | o Calculate the changes to Hg species because            |
C     |   photo- and bio-chemical reactions                      |
C     |   Modified based on Anne Soerensen's code for GEOS-Chem  |
C     |==========================================================|
      IMPLICIT NONE
      INCLUDE 'chem_hg_Parameters.h'
      INCLUDE 'chem_hg_Global.h'

C     == Routine arguments ==
      LOGICAL, SAVE ::   FIRST
      DATA FIRST / .TRUE. /   
      REAL*8 T0, DeltaT
      REAL*8 CONC(NVAR)
      INTEGER i

C     time variables
      T0 = 0.0

C     Set tolerances
      RTOLS = 0.001 
      DO i=1,NVAR
        RTOL(i) = RTOLS  ! Relative tolerance
        ATOL(i) = 0.0001  ! Absolute tolerance, 0.0001 fM
      END DO
      ATOL(11) = 0.0001  !DMHG
      ATOL(12) = 0.0001  !MMHG
      ATOL(13) = 0.0001  !HG0
      ATOL(14) = 0.0001  !HG2

C     Initialize
      IF ( FIRST ) THEN
        CALL INITIALIZE
        FIRST = .FALSE.
      ENDIF
                                          
C     Map the vector CONC to the KPP concentration 
C      vector from chem_hg_Global.f.
      DO i=1,NVAR
        C(i) = CONC(i)
      ENDDO

C     Update the rate constants just prior to integration
      CALL Update_RCONST

C     Integrate for one step
      CALL INTEGRATE( T0, DeltaT )

C     Return results
      DO i=1,NVAR
        CONC(i) = C(i)
      ENDDO

      RETURN
      END

