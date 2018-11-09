!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pfc_init_fixed.F
!
! !DESCRIPTION: Sets up diagnostics for the PFC simulation.
!\\
!\\
! !INTERFACE:
!
#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "PFC_OPTIONS.h"

      SUBROUTINE PFC_INIT_FIXED( myThid )
!
! !USES:
!
      IMPLICIT NONE

      ! Global variables
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "GCHEM.h"
#include "PFC_SIZE.h"
#include "PFC.h"
!
! !INPUT PARAMETERS: 
!
      ! myThid -- Index of the thread we are on
      INTEGER myThid
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Feb 2016 - R. Yantosca - Added ProTex headers
!  29 Feb 2016 - R. Yantosca - Moved #ifdef ALLOW_PTRACERS inside routine
!  29 Feb 2016 - R. Yantosca - Moved #ifdef ALLOW_PFC      inside routine
!EOP
!------------------------------------------------------------------------------
!BOC
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PFC
#ifdef ALLOW_DIAGNOSTICS

      !================================================================
      ! %%%%%     Set up diagnostics for the PFC simulation      %%%%%
      !================================================================
      IF ( useDIAGNOSTICS ) THEN
         CALL PFC_DIAGNOSTICS_INIT( myThid )
      ENDIF

#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_PFC         */
#endif /* ALLOW_PTRACERS    */

      RETURN
      END
!EOC
