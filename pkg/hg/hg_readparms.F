!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: hg_readparms.F
!
! !DESCRIPTION: Reads parameters for the Hg simulation.
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "HG_OPTIONS.h"

      SUBROUTINE HG_READPARMS( myThid )
!
! !USES:
!
      IMPLICIT NONE

      ! Global variables
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "HG_SIZE.h"
#include "HG.h"
!
! !INPUT PARAMETERS: 
!
      ! myThid :: My Thread Id. number
      INTEGER myThid
!
! !REMARKS:
!     *==========================================================*
!     | S/R HG_READPARMS
!     | o Initialise and read hg package parameters
!     *==========================================================*
!
! !REVISION HISTORY:
!  29 Feb 2016 - R. Yantosca - Added ProTex headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#ifdef ALLOW_HG

C     === Local variables ===
C     msgBuf    :: Informational/error message buffer
C     iUnit     :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      NAMELIST /HG_FORCING/
     &          HG_windFile, HG_iceFile, radsw_File,
     &          HG_forcingPeriod, HG_forcingCycle,
     &          HG_chemPeriod, HG_chemCycle,
     &          HG0_concFile, HG2_depoFile,
     &          HG2_riverFile, HGP_riverFile,
#ifdef HG2REFRACTORY
     &          HG2R_riverFile,
#endif
#ifdef HGPREFRACTORY
     &          HGPR_riverFile,
#endif     
     &          chl_file, npp_file,
     &          doc_file, poc_file,
     &          Fdoc_file, Fpoc_file, Fremin_file,
     &          o2_file, dfe_file
#ifdef FOODW_MODE     
     &          ,NO3_file,phytoP_file, zooP_file 
     &          ,phygrow_file,zoograze_file,zoomort_file           
#endif     

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      _BEGIN_MASTER(myThid)

       HG_windFile  = ' '
       HG_iceFile   = ' '
       radsw_File   = ' '
       HG0_concFile = ' '
       HG2_depoFile = ' '
       HG2_riverFile= ' '
       HGP_riverFile= ' '
#ifdef HG2REFRACTORY
       HG2R_riverFile= ' '
#endif
#ifdef HGPREFRACTORY
       HGPR_riverFile= ' '
#endif       
       chl_file     = ' '
       npp_file     = ' '
       doc_file     = ' '
       poc_file     = ' '
       Fdoc_file    = ' '
       Fpoc_file    = ' '
       Fremin_file  = ' '
       o2_file      = ' '
       dfe_file     = ' '
#ifdef FOODW_MODE       
       NO3_file     = ' '
       phytoP_file  = ' '
       zooP_file    = ' '
       phygrow_file = ' '
       zoograze_file= ' '
       zoomort_file = ' '
#endif       
c air-sea exchange: default periodic forcing to same as for physics
       HG_forcingPeriod = externForcingPeriod
       HG_forcingCycle  = externForcingCycle
c hg chemistry: default periodic forcing to same as for physics
       HG_chemPeriod = externForcingPeriod
       HG_chemCycle  = externForcingCycle

      WRITE(msgBuf,'(A)') ' HG_READPARMS: opening data.hg'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )

      CALL OPEN_COPY_DATA_FILE( 'data.hg', 'HG_READPARMS',
     O                          iUnit, myThid )

C--   Read parameters from open data file:

C-    forcing filenames and parameters
      READ(UNIT=iUnit,NML=HG_FORCING)

      WRITE(msgBuf,'(A)')
     &   ' HG_READPARMS: finished reading data.hg'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )

C--   Close the open data file
      CLOSE(iUnit)

      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#endif /* ALLOW_HG */

      RETURN
      END
!EOC
