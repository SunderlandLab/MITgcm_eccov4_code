C $Header: /u/gcmpack/MITgcm/pkg/gchem/gchem_fields_load.F,v 1.9 2013/06/10 02:58:12 jmc Exp $
C $Name:  $

#include "GCHEM_OPTIONS.h"

CBOP
C !ROUTINE: GCHEM_FIELDS_LOAD

C !INTERFACE: ==========================================================
      SUBROUTINE GCHEM_FIELDS_LOAD (
     I           myTime, myIter, myThid)

C !DESCRIPTION:
C  calls routines which read in fields needed for any tracer experiment
C
C !MODIFICATION HISTORY
C 16 May 2013 - Yanxu Zhang    - Add Hg support
C 16 May 2014 - Xianming Zhang - Add PFCs
C 08 Jul 2014 - Helen Amos     - Add PCBs
C 29 Feb 2016 - R. Yantosca    - Update comment headers 
C 29 Feb 2016 - R. Yantosca    - Bracket HG_FIELDS_LOAD  with IF ( useHg  )
C 29 Feb 2016 - R. Yantosca    - Bracket PFC_FIELDS_LOAD with IF ( usePFC )
C 29 Feb 2016 - R. Yantosca    - Bracket PCB_FIELDS_LOAD with IF ( usePCB )
C
C !USES: ===============================================================
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "GCHEM.h"

C !INPUT PARAMETERS: ===================================================
C  myTime               :: current time
C  myIter               :: current iteration
C  myThid               :: thread number
      _RL myTime
      INTEGER myIter
      INTEGER myThid

C !LOCAL VARIABLES: ====================================================
CEOP

#ifdef ALLOW_GCHEM
      !=================================================================
      ! %%%%% LOAD EXTERNAL DATA FOR THE VARIOUS GCHEM PACKAGES  %%%%%
      !=================================================================

#ifdef ALLOW_DIC
      !================================================================
      ! %%%%%         Load fields for the DIC package             %%%%%
      !================================================================
      IF ( useDIC ) THEN
         CALL DIC_FIELDS_LOAD( myIter, myTime, myThid )
      ENDIF
#endif

#ifdef ALLOW_CFC
      !================================================================
      ! %%%%%         Load fields for the CFC package            %%%%%
      !================================================================
      IF ( useCFC ) THEN
         CALL CFC_FIELDS_LOAD( myTime, myIter, myThid )
      ENDIF
#endif

#ifdef ALLOW_HG
      !================================================================
      ! %%%%%          Load fields for the HG package            %%%%%
      ! %%%%%             (Yanxu Zhang, 5/16/13)                 %%%%%
      !================================================================
      IF ( useHG ) THEN
         CALL HG_FIELDS_LOAD( myIter, myTime, myThid )
      ENDIF
#endif

#ifdef ALLOW_PFC
      !================================================================
      ! %%%%%          Load fields for the PFC package           %%%%%
      ! %%%%%            (Xianming Zhang, 5/16/2014)             %%%%%
      !================================================================
      IF ( usePFC ) THEN
         CALL PFC_FIELDS_LOAD( myIter, myTime, myThid )
      ENDIF
#endif

#ifdef ALLOW_PCB
      !================================================================
      ! %%%%%          Load fields for the PCB package           %%%%%
      ! %%%%%               (Helen Amos, 7/8/14)                 %%%%%
      !================================================================
      IF ( usePCB ) THEN
         CALL PCB_FIELDS_LOAD( myIter, myTime, myThid )
      ENDIF
#endif

#ifdef ALLOW_DARWIN
      !================================================================
      ! %%%%%         Load fields for the DARWIN packag e        %%%%%
      !================================================================
      IF ( useDARWIN ) THEN
         CALL DARWIN_FIELDS_LOAD( myIter, myTime ,myThid )
      ENDIF
#endif

#endif /* ALLOW_GCHEM */

      RETURN
      END

