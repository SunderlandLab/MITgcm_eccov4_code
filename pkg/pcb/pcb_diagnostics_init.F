!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pcb_diagnostics_init.F
!
! !DESCRIPTION: Initializes diagnostics for the PCB simulation.
!\\
!\\
! !INTERFACE:
!
#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "PCB_OPTIONS.h"

      SUBROUTINE PCB_DIAGNOSTICS_INIT( myThid )
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
#include "PCB_SIZE.h"
#include "PCB.h"
!
! !INPUT PARAMETERS: 
!
      ! myThid -  Number of this instance of HG_DIAGNOSTICS_INIT
      INTEGER myThid
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jul 2014 - H. Amos     - Comment out PCBi tracer. Rename PCBn to PCBa.    
!  03 Feb 2015 - H. Amos     - Update rivers based on Yanxu Zhang's Hg code.
!  01 Mar 2016 - R. Yantosca - Added ProTex headers
!  01 Mar 2016 - R. Yantosca - Moved #ifdef ALLOW_TRACERS into routine
!  01 Mar 2016 - R. Yantosca - Moved #ifdef ALLOW_PCB     into routine
!  01 Mar 2016 - R. Yantosca - Removed commented-out code
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PCB
#ifdef ALLOW_DIAGNOSTICS
      INTEGER        diagNum
      CHARACTER*8    diagName
      CHARACTER*16   diagCode
      CHARACTER*16   diagUnits
      CHARACTER*(80) diagTitle
      INTEGER        np

      !================================================================
      ! %%%%%            Diagnostics for PCBa tracer             %%%%%
      !================================================================

      diagName = 'PCBaSPRY'
      WRITE(diagTitle, '(A)') 'PCBa seaspray '
      diagUnits = 'mol/m^2/s          '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )


      diagName = 'PCBaCONC'
      WRITE(diagTitle, '(A)') 'Atmospheric PCBa concentration'
      diagUnits = 'mol/m^3          '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PCBaEVAS'
      WRITE(diagTitle, '(A)') 'PCBa evasion flux from ocean'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
          
      diagName = 'PCBaDEPO'
      WRITE(diagTitle, '(A)') 'Atmospheric PCBa deposition flux'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PCBaRIVR'
      WRITE(diagTitle, '(A)') 'PCBa riverine runoff'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      ! Mimic Hg for PCB river inputs
      ! (H. Amos, 03 Feb 2015) 
#ifdef PCBAREFRACTORY
      diagName = 'PCBaRRIV'
      WRITE(diagTitle, '(A)') 'Refractory riverine PCBa riverine runoff'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )          

      diagName = 'PCBaRDIS'
      WRITE(diagTitle, '(A)') 'Refractory river PCBa dissolution flux'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#endif /* PCBAREFRACTORY */
     
      diagName = 'PCBaPART'
      WRITE(diagTitle, '(A)') 'PCBa partitioned onto particles'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PCBaDEGR'
      WRITE(diagTitle, '(A)') 'PCBa degradation flux'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )    

      !================================================================
      ! %%%%%            Diagnostics for PCBp tracer             %%%%%
      !================================================================

      diagName = 'PCBpDEPO'
      WRITE(diagTitle, '(A)') 'Atmospheric PCBp deposition flux'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PCBpRIVR'
      WRITE(diagTitle, '(A)') 'PCBp riverine runoff'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid ) 

#ifdef PCBAPREFRACTORY

      diagName = 'PCBPRSNK'
      WRITE(diagTitle, '(A)') 'Refrac PCBp sinking flux at level bottom'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PCBPRRIV'
      WRITE(diagTitle, '(A)') 'Refrac. PCBp riverine runoff'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )          
          
      diagName = 'PCBPRDIS'
      WRITE(diagTitle, '(A)') 'Refractory river PCBp dissolution flux'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#endif /* PCBAREFRACTORY */

      diagName = 'PCBpPART'
      WRITE(diagTitle, '(A)') 'PCBp partitioned onto particles'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName = 'PCBpSINK'
      WRITE(diagTitle, '(A)') 'PCBp sinking flux at level bottom'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )         

#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_PCB         */
#endif /* ALLOW_PTRACERS    */
      RETURN
      END     
!EOC
