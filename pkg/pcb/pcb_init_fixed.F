!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pcb_init_fixed.F
!
! !DESCRIPTION: Sets up diagnostics for the PCB simulation.
!\\
!\\
! !INTERFACE:
!
#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "PCB_OPTIONS.h"

      SUBROUTINE PCB_INIT_FIXED( myThid )
!
! !USES:
!
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "GCHEM.h"
#include "PCB_SIZE.h"
#include "PCB.h"
!
! !INPUT PARAMETERS: 
!
      ! myThid - thread ID #
      INTEGER myThid
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Feb 2016 - R. Yantosca - Added ProTex headers
!  01 Mar 2016 - R. Yantosca - Moved #ifdef ALLOW_PTRACERS inside routine
!  01 Mar 2016 - R. Yantosca - Moved #ifdef ALLOW_PCB      inside routine
!EOP
!------------------------------------------------------------------------------
!BOC
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PCB
#ifdef ALLOW_DIAGNOSTICS

      IF ( useDIAGNOSTICS ) THEN
        CALL PCB_DIAGNOSTICS_INIT( myThid )
      ENDIF

#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_PCB         */
#endif /* ALLOW_PTRACERS    */

      RETURN
      END
!EOC
