C $Header: /u/gcmpack/MITgcm/pkg/gchem/gchem_init_fixed.F,v 1.11 2009/06/30 16:43:00 jahn Exp $
C $Name:  $

#include "GCHEM_OPTIONS.h"

C !INTERFACE: ==========================================================
      SUBROUTINE GCHEM_INIT_FIXED(myThid )

C !DESCRIPTION:
C calls subroutines that initialized fixed variables for any tracer
C experiment
C
C !MODIFICATION HISTORY:
C 16 May 2013 - Yanxu Zhang    - Add Hg support
C        2014 - Xianming Zhang - Add PFCs
C 08 Jul 2014 - Helen Amos     - Add PCBs
C 29 Feb 2016 - R. Yantosca    - Cosmetic changes, update comments


C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GCHEM.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid
CEOP

#ifdef ALLOW_GCHEM
      !=================================================================
      ! %%%%%  INITIALIZE FIXED DATA FOR VARIOUS GCHEM PACKAGES   %%%%%
      !=================================================================

C jmc: simpler to put MASTER thread directive within individual S/R :
C      that leave the choice, for each S/R, to do either single (master)
C      thread tasks or Multi-threaded tasks.
c        _BEGIN_MASTER( mythid )


#ifdef ALLOW_DIC
      !================================================================
      ! %%%%%  Initialize fixed variables for the DIC package    %%%%%
      !================================================================
      IF ( useDIC ) THEN
          CALL DIC_INIT_FIXED( myThid )
      ENDIF
#endif /* ALLOW_DIC */

#ifdef ALLOW_CFC
      !================================================================
      ! %%%%%  Initialize fixed variables for the CFC package    %%%%%
      !================================================================
      IF ( useCFC ) THEN
          CALL CFC_PARAM( myThid )
          CALL CFC_ATMOS( myThid )
      ENDIF
#endif /* ALLOW_CFC */

#ifdef ALLOW_DARWIN
      !================================================================
      ! %%%%%  Initialize fixed variables for DARWIN package     %%%%%
      !================================================================
      IF ( useDARWIN ) THEN
         CALL DARWIN_INIT_FIXED( myThid )
      ENDIF
#endif

#ifdef ALLOW_HG
      !================================================================
      ! %%%%%  Initialize fixed variables for the HG package     %%%%%
      ! %%%%%             (Yanxu Zhang, 5/16/13)                 %%%%%
      !================================================================
      IF ( useHG ) THEN
         CALL HG_INIT_FIXED( myThid )
      ENDIF
#endif

#ifdef ALLOW_PFC
      !================================================================
      ! %%%%%  Initialize fixed variables for the PFC package    %%%%%
      ! %%%%%            (Xianming Zhang, 5/16/2014)             %%%%%
      !================================================================
      IF ( usePFC ) THEN
         CALL PFC_PARAM( myThid )
         CALL PFC_INIT_FIXED( myThid )
      ENDIF
#endif

#ifdef ALLOW_PCB
      !================================================================
      ! %%%%%  Initialize fixed variables for the PCB package    %%%%%
      ! %%%%%               (Helen Amos, 7/8/14)                 %%%%%
      !================================================================
      IF ( usePCB ) THEN
         call PCB_PARAM( myThid )
         call PCB_INIT_FIXED( myThid )
      ENDIF
#endif

#ifdef ALLOW_DIAGNOSTICS
      !================================================================
      ! %%%%%   Define diagnostics names for the GCHEM package   %%%%%
      !================================================================
      IF ( useDiagnostics ) THEN
        CALL GCHEM_DIAGNOSTICS_INIT( myThid )
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_GCHEM */

c         _END_MASTER( mythid )

      RETURN
      END
