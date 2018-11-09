!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: hg_init_fixed
!
! !DESCRIPTION: Sets up diagnostics for the Hg simulation.
!\\
!\\
! !INTERFACE:
!
#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "HG_OPTIONS.h"

      SUBROUTINE HG_INIT_FIXED( myThid )
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
#include "HG_SIZE.h"
#include "HG.h"
!
! !INPUT PARAMETERS: 
!
      INTEGER myThid
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Feb 2016 - R. Yantosca - Added ProTex headers
!EOP
!------------------------------------------------------------------------------
!BOC
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_HG
#ifdef ALLOW_DIAGNOSTICS

      !================================================================
      ! %%%%%      Set up diagnostics for the Hg simulation      %%%%%
      !================================================================
      IF ( useDIAGNOSTICS ) THEN
         CALL HG_DIAGNOSTICS_INIT( myThid )
      ENDIF

#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_HG        */
#endif /* ALLOW_PTRACERS */
      RETURN
      END
!EOC
