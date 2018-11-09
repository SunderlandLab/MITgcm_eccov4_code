!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pfc_diagnostics_init.F
!
! !DESCRIPTION: Sets up diagnostics for the PFC (aka PFOS) simulation.
!\\
!\\
! !INTERFACE:
!
#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "PFC_OPTIONS.h"

      SUBROUTINE PFC_DIAGNOSTICS_INIT( myThid )
!
! !USES:
!
      IMPLICIT NONE

      ! Global variables
#include "SIZE.h"
#include "PTRACERS_SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GCHEM.h"
#include "PFC_SIZE.h"
#include "PFC.h"
!
! !INPUT PARAMETERS: 
!
      ! myThid -  Number of this instance of HG_DIAGNOSTICS_INIT
      INTEGER myThid
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Feb 2016 - R. Yantosca - Added ProTex headers
!  29 Feb 2016 - R. Yantosca - Moved #ifdef ALLOW_PTRACERS inside routine
!  29 Feb 2016 - R. Yantosca - Moved #ifdef ALLOW_PFC      inside routine
!  02 Mar 2016 - R. Yantosca - Reorder diagnostics by tracer
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PFC
#ifdef ALLOW_DIAGNOSTICS
      INTEGER        diagNum
      CHARACTER*8    diagName
      CHARACTER*16   diagCode
      CHARACTER*16   diagUnits
      CHARACTER*(80) diagTitle

      INTEGER         np

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|


      !================================================================
      ! %%%%%      Add PFCn diagnostics to the (long) list       %%%%%
      !================================================================

      diagName = 'PFCnCONC'
      WRITE(diagTitle, '(A)') 'Atmospheric PFCn concentration'
      diagUnits = 'ng/m^3          '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PFCnSPRY'
      WRITE(diagTitle, '(A)') 'PFCn seaspray '
      diagUnits = 'mol/m^2/s          '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PFCnEVAS'
      WRITE(diagTitle, '(A)') 'PFCn evasion flux from ocean'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PFCnRIVR'
      WRITE(diagTitle, '(A)') 'PFCn riverine runoff'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PFCnPART'
      WRITE(diagTitle, '(A)') 'PFCn partitioned onto particles'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      !================================================================
      ! %%%%%      Add PFCi diagnostics to the (long) list       %%%%%
      !================================================================

      diagName = 'PFCiSPRY'
      WRITE(diagTitle, '(A)') 'PFCi seaspray '
      diagUnits = 'mol/m^2/s          '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PFCiDEPO'
      WRITE(diagTitle, '(A)') 'Atmospheric PFCi deposition flux'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )


      diagName = 'PFCiRIVR'
      WRITE(diagTitle, '(A)') 'PFCi riverine runoff'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PFCiPART'
      WRITE(diagTitle, '(A)') 'PFCi partitioned onto particles'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      !================================================================
      ! %%%%%      Add PFCp diagnostics to the (long) list       %%%%%
      !================================================================

      diagName = 'PFCpRIVR'
      WRITE(diagTitle, '(A)') 'PFCp riverine runoff'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PFCpPART'
      WRITE(diagTitle, '(A)') 'PFCp partitioned onto particles'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName = 'PFCpSINK'
      WRITE(diagTitle, '(A)') 'PFCp sinking flux at level bottom'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )  


#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_PFC         */
#endif /* ALLOW_PTRACERS    */
      RETURN
      END     
!EOC
