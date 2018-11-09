C $Header: /u/gcmpack/MITgcm/pkg/gchem/gchem_calc_tendency.F,v 1.5 2013/06/10 02:52:57 jmc Exp $
C $Name:  $

#include "GCHEM_OPTIONS.h"

CBOP
C !ROUTINE: GCHEM_CALC_TENDENCY

C !INTERFACE: ==========================================================
      SUBROUTINE GCHEM_CALC_TENDENCY(
     I     myTime, myIter, myThid )

C !DESCRIPTION:
C In the case of GCHEM_SEPARATE_FORCING not defined,
C this subroutine computes the tracer tendencies due to a
C bio-geogchemistry or ecosystem model and stores them on an array
C gchemTendency, that will be incorporated into regular timestepping in
C in ptracers_intergrate.F
C The current example uses the CFC package, but his is meant to
C be replaced by anything that the user provides.
C Add Hg support, yxzhang, 5/16/13
C Add PFOS/PFC and PCB (bmy, 2/29/16)

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "GCHEM.h"
#include "GCHEM_FIELDS.h"
#ifdef ALLOW_HG
#include "HG_OPTIONS.h"
#include "HG_SIZE.h"
#endif
#ifdef ALLOW_PFC
#include "PFC_OPTIONS.h"
#include "PFC_SIZE.h"
#endif
#ifdef ALLOW_PCB
#include "PCB_OPTIONS.h"
#include "PCB_SIZE.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  myTime               :: current time
C  myIter               :: current timestep
C  myThid               :: thread number
      _RL myTime
      INTEGER myIter, myThid
CEOP

#ifdef ALLOW_GCHEM
# ifndef GCHEM_SEPARATE_FORCING
C !LOCAL VARIABLES: ====================================================
C  i,j                  :: loop indices
C  iTracer              :: ptracer number
C  bi,bj                :: tile indices
C  k                    :: vertical level
      INTEGER i,j
      INTEGER iTracer
      INTEGER bi,bj,iMin,iMax,jMin,jMax,k

C     gchemTendency is re-initialized here
      DO iTracer = 1, PTRACERS_numInUse
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO k = 1, Nr
          DO j = 1-OLy, sNy+OLy
           DO i = 1-OLx, sNx+OLx
            gchemTendency(i,j,k,bi,bj,iTracer) = 0. _d 0
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
C
C     Here is the place for code to compute bio-geochemical
C     tendency terms (sometimes referred to as source-minus-sink
C     terms). The tendencies are stored on gchemTendency, as show
C     in the CFC-example.
C
C     loop over tiles
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
C     define horizontal loop ranges
        iMin = 1
        iMax = sNx
        jMin = 1
        jMax = sNy
c       DO iTracer = 1, PTRACER_numInUse
c        DO k = 1, Nr
c         DO j = 1-OLy, sNy+OLy
c          DO i = 1-OLx, sNx+OLx
c           gchemTendency(i,j,k,bi,bj,iTracer) = your specific model
c          ENDDO
c         ENDDO
c        ENDDO
c       ENDDO


#ifdef ALLOW_CFC
        !==============================================================
        ! %%%%%       CHEMICAL FORCING BY CFC11 and CFC12        %%%%%
        !==============================================================
        IF ( useCFC ) THEN
           CALL CFC11_FORCING(
     I                          pTracer     (1-OLx,1-OLy,1,bi,bj,1),
     U                          gchemTendency(1-OLx,1-OLy,1,bi,bj,1),
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myTime, myIter, myThid )

           CALL CFC12_FORCING(
     I                          pTracer     (1-OLx,1-OLy,1,bi,bj,2),
     U                          gchemTendency(1-OLx,1-OLy,1,bi,bj,2),
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myTime, myIter, myThid )
        ENDIF
#endif /* ALLOW_CFC */

c yxzhang, 5/16/13
#ifdef ALLOW_HG
        !==============================================================
        ! %%%%%         CHEMICAL FORCING BY Hg SPECIES           %%%%%
        !==============================================================
        IF ( useHG ) THEN

           ! Partitioning, this should be put as the first step
           CALL HG_PARTITION(   pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )

           ! Bioaccumulation, this should be put as the second step
           CALL HG_BIOACCUM(    pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )

           ! Atmospheric forcings
           CALL HG_FORCING(     pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                          gchemTendency(1-OLx,1-OLy,1,1,1,1),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )

           ! Chmistry, redox + methylation + demethylation
           CALL HG_CHEMISTRY(   pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                          gchemTendency(1-OLx,1-OLy,1,1,1,1),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )

           ! chemistry, riverine particle dissolution
           CALL HG_DISSOLUTION( pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                          gchemTendency(1-OLx,1-OLy,1,1,1,1),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )

           ! Foodweb dynamics
           CALL HG_FOODWEB(     pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                          gchemTendency(1-OLx,1-OLy,1,1,1,1),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )

           ! Biological pump, particulate matter sinking
           CALL HG_BIOPUMP(     pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                          gchemTendency(1-OLx,1-OLy,1,1,1,1),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )

        ENDIF

#endif /* ALLOW_HG */


#ifdef ALLOW_PFC
        !==============================================================
        ! %%%%%         CHEMICAL FORCING BY PFC SPECIES          %%%%%
        !==============================================================
        IF ( usePFC ) THEN

!          CALL PFC_ARCTICBOUNDARY( pTracer     (1-OLx,1-OLy,1,1,1,1),
!     &                             bi, bj, iMin, iMax, jMin, jMax,
!     &                             myIter, myTime, myThid )

           ! Partitioning
           !write( 6, '(a)' ) 'xz checking call PFC_PARTITION'
           CALL PFC_PARTITION(  pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )

           ! Atmospheric forcings
           CALL PFC_FORCING(    pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                          gchemTendency(1-OLx,1-OLy,1,1,1,1),  
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )

!           ! Partitioning
!           CALL PFC_PARTITION(  pTracer     (1-OLx,1-OLy,1,1,1,1),
!     &                          bi, bj, iMin, iMax, jMin, jMax,
!     &                          myIter, myTime, myThid )

!           ! Chemistry, no chemistry considered for PFC -XZ
!           CALL PFC_CHEMISTRY(  pTracer     (1-OLx,1-OLy,1,1,1,1),
!     &                          gchemTendency(1-OLx,1-OLy,1,1,1,1),   
!     &                          bi, bj, iMin, iMax, jMin, jMax,
!     &                          myIter, myTime, myThid )

        ! Biological pump, particulate matter sinking
           CALL PFC_BIOPUMP(    pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                          gchemTendency(1-OLx,1-OLy,1,1,1,1), 
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )

!          CALL PFC_ARCTICBOUNDARY( pTracer     (1-OLx,1-OLy,1,1,1,1),
!     &                             bi, bj, iMin, iMax, jMin, jMax,
!     &                             myIter, myTime, myThid )

        ENDIF
#endif /* ALLOW_PFC */

#ifdef ALLOW_PCB
        !==================================================================== 
        ! %%%%%       CHEMICAL FORCING BY PCB SPECIES %%%%%
        !====================================================================
        IF ( usePCB ) THEN

           ! Partitioning
           !write( 6, '(a)' ) 'xz checking call PCB_PARTITION'
           CALL PCB_PARTITION(  pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )

           ! Atmospheric forcings
           CALL PCB_FORCING(    pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                          gchemTendency(1-OLx,1-OLy,1,1,1,1),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )

!           ! Partitioning 
!           CALL PCB_PARTITION(  pTracer     (1-OLx,1-OLy,1,1,1,1),
!     &                          bi, bj, iMin, iMax, jMin, jMax,
!     &                          myIter, myTime, myThid )

           ! Degradation
           CALL PCB_CHEMISTRY(  pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                          gchemTendency(1-OLx,1-OLy,1,1,1,1),   
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )

           ! Biological pump, particulate matter sinking
           CALL PCB_BIOPUMP(    pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                          gchemTendency(1-OLx,1-OLy,1,1,1,1),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )

        ENDIF
#endif /* ALLOW_PCB */

C     end of tile-loops
       ENDDO
      ENDDO
# endif /* GCHEM_SEPARATE_FORCING */
#endif /* ALLOW_GCHEM */

      RETURN
      END
