!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: hg_diagnostics_init.F
!
! !DESCRIPTION: Routine to initialize Hg diagnostics
!\\
!\\
! !INTERFACE:
!
#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "HG_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_HG

      SUBROUTINE HG_DIAGNOSTICS_INIT( myThid )
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
#include "HG_SIZE.h"
#include "HG.h"
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

#ifdef ALLOW_DIAGNOSTICS
C     !LOCAL VARIABLES:
C     === Local variables ===
C     msgBuf      - Informational/error meesage buffer
c     CHARACTER*(MAX_LEN_MBUF) msgBuf

      INTEGER       diagNum
      CHARACTER*8   diagName
      CHARACTER*16  diagCode
      CHARACTER*16  diagUnits
      CHARACTER*(80) diagTitle

      INTEGER      np,nz

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C--   Add diagnostics to the (long) list

C--   Add diagnostic HG2DEPO, HG0CONC, HG0EVAS, HG2REDN, and HGPSINK
      diagName = 'HG2DEPO '
      WRITE(diagTitle, '(A)') 'Atmospheric Hg2 deposition flux'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'HG0CONC '
      WRITE(diagTitle, '(A)') 'Atmospheric Hg0 concentration'
      diagUnits = 'ng/m^3          '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'HG2RIVER'
      WRITE(diagTitle, '(A)') 'Hg2 riverine runoff'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'HGPRIVER'
      WRITE(diagTitle, '(A)') 'HgP riverine runoff'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'HG0EVAS '
      WRITE(diagTitle, '(A)') 'Hg0 evasion flux from ocean'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName = 'HG2REDN '
      WRITE(diagTitle, '(A)') 'Net Hg2 reduction rate'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName = 'HG2RBIO '
      WRITE(diagTitle, '(A)') 'Gross Hg2 bio-reduction rate'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'HG2RPHO '
      WRITE(diagTitle, '(A)') 'Gross Hg2 photo-reduction rate'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'HG0OBIO '
      WRITE(diagTitle, '(A)') 'Gross Hg0 bio-oxidation rate'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
                    
      diagName = 'HG0OPHO '
      WRITE(diagTitle, '(A)') 'Gross Hg0 photo-oxidation rate'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName = 'HG0ODAK '
      WRITE(diagTitle, '(A)') 'Gross Hg0 dark-oxidation rate'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName = 'HGPSINK '
      WRITE(diagTitle, '(A)') 'Hgp sinking flux at level bottom'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'HG2PART '
      WRITE(diagTitle, '(A)') 'Hg2 partitioned onto particles'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#ifdef HG2REFRACTORY
      diagName = 'HG2RRIVE'
      WRITE(diagTitle, '(A)') 'Refractory riverine Hg2 riverine runoff'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'HG2RDISS'
      WRITE(diagTitle, '(A)') 'Refractory riverine Hg2 dissolution flux'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* HG2REFRACTORY */
     
#ifdef HGPREFRACTORY
      diagName = 'HGPRSINK'
      WRITE(diagTitle, '(A)') 'Refrac Hgp sinking flux at level bottom'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'HGPRRIVE'
      WRITE(diagTitle, '(A)') 'Refrac. HgP riverine runoff'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
          
      diagName = 'HGPRDISS'
      WRITE(diagTitle, '(A)') 'Refractory river HgP dissolution flux'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* HGPREFRACTORY */
               
#ifdef METHY_MODE
      diagName = 'DMHGEVAS'
      WRITE(diagTitle, '(A)') 'DMHg evasion flux from ocean'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName = 'MHGPSINK'
      WRITE(diagTitle, '(A)') 'MMHgp sinking flux at level bottom'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName = 'MMHGPART'
      WRITE(diagTitle, '(A)') 'MMHg partitioned onto particles'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'HG2METB '
      WRITE(diagTitle, '(A)') 'Hg2 biological methylation'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'MMHGDMEB'
      WRITE(diagTitle, '(A)') 'MMHg biological demethylation'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'MMHGDMEP'
      WRITE(diagTitle, '(A)') 'MMHg photochemical demethylation'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'MMHGMETB'
      WRITE(diagTitle, '(A)') 'MMHg biological methylation to DMHG'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'DMHGDMEP'
      WRITE(diagTitle, '(A)') 'DMHg photochemical demethylation'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#ifdef FOODW_MODE
      DO np=1,npmax
       WRITE(diagName,'(A6,I2.2)') 'MMHGBC',np
       WRITE(diagTitle,'(A,I2.2)') 'Phytoplankton bioaccumulation', np
       diagUnits = 'mol/m^3/s       '
       diagCode  = 'SM P    MR      '
       CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      ENDDO
      DO np=1,npmax
      DO nz=1,nzmax
       WRITE(diagName,'(A4,I2.2,I2.2)') 'MHGZ',np,nz
       WRITE(diagTitle,'(A,I2.2,I2.2)') 'Zooplankton graze', np, nz
       diagUnits = 'mol/m^3/s       '
       diagCode  = 'SM P    MR      '
       CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      ENDDO
      ENDDO
      DO nz=1,nzmax
       WRITE(diagName,'(A6,I2.2)') 'MHMORT',nz
       WRITE(diagTitle,'(A,I2.2)') 'Release of dead zooplankton', nz
       diagUnits = 'mol/m^3/s       '
       diagCode  = 'SM P    MR      '
       CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      ENDDO
      diagName = 'MMHGFLOP'
      WRITE(diagTitle, '(A)') 'Release of MMHg via floopy feeding'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* FOODW_MODE        */

#endif /* METHY_MODE        */
     
#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END     
#endif /* ALLOW_HG          */
#endif /* ALLOW_PTRACERS    */
!EOC
