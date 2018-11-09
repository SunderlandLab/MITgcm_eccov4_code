!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pfc_readparms.F
!
! !DESCRIPTION: Reads input data for the PFC simulation.
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "PFC_OPTIONS.h"

      SUBROUTINE PFC_READPARMS( myThid )
!
! !USES:
!
      IMPLICIT NONE

      ! Global variables
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "PFC_SIZE.h"
#include "PFC.h"

#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h"
# include "PTRACERS_PARAMS.h"
#endif
!
! !INPUT PARAMETERS: 
!
      ! myThid :: My Thread Id. number
      INTEGER myThid
!
! !REMARKS:
!     *==========================================================*
!     | S/R PFC_READPARMS
!     | o Initialise and read pfc package parameters
!     *==========================================================*
!
! !REVISION HISTORY:
!  29 Feb 2016 - R. Yantosca - Added ProTex headers
!  03 Mar 2016 - R. Yantosca - Bug fix: set o2_file and dfe_file to ' '
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

#ifdef ALLOW_PFC

C     === Local variables ===
C     msgBuf    :: Informational/error message buffer
C     iUnit     :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      NAMELIST /PFC_FORCING/
     &          PFC_windFile, PFC_iceFile,
     &          PFC_forcingPeriod, PFC_forcingCycle,
c xz:the following lines are modified from hg code  
c PFCn-neutral PFC; PFCi-ionized PFC; PFCp-particle bounded PFC      
     &          PFC_chemPeriod, PFC_chemCycle,
     &          PFC_emisPeriod, PFC_emisCycle,
     &          PFCn_concFile, PFCi_depoFile,
     &          PFCi_riverFile, PFCp_riverFile,
     &          chl_file, npp_file,
     &          doc_file, poc_file,
     &          Fdoc_file, Fpoc_file, Fremin_file
     
     

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      _BEGIN_MASTER(myThid)

      ! Set all files in the common block /PFC_FILENAMES/ (located in
      ! include file PFC.h) to ' ', which will skip them from being 
      ! read in from disk.  Also initialize o2_file and dfe_file,
      ! which seemed to be missing here. (bmy, 3/3/16)
      PFC_windFile   = ' '
      PFC_iceFile    = ' '
      radsw_File     = ' '
      PFCn_concFile  = ' '
      PFCi_depoFile  = ' '
      PFCi_riverFile = ' '
      PFCp_riverFile = ' '      
      chl_file       = ' '
      npp_file       = ' '
      doc_file       = ' '
      poc_file       = ' '
      Fdoc_file      = ' '
      Fpoc_file      = ' '
      Fremin_file    = ' '
      o2_file        = ' '
      dfe_file       = ' '
      
c default periodic forcing to same as for physics
      PFC_forcingPeriod = externForcingPeriod
      PFC_forcingCycle  = externForcingCycle
c====== added based on hg code default periodic forcing to same as for physics
      PFC_chemPeriod = externForcingPeriod
      PFC_chemCycle  = externForcingCycle       
c=======
      PFC_emisPeriod = externForcingPeriod
      PFC_emisCycle  = externForcingCycle       



      WRITE(msgBuf,'(A)') ' PFC_READPARMS: opening data.pfc'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )

      CALL OPEN_COPY_DATA_FILE( 'data.pfc', 'PFC_READPARMS',
     O                          iUnit, myThid )

C--   Read parameters from open data file:

C-    forcing filenames and parameters
      READ(UNIT=iUnit,NML=PFC_FORCING)

      WRITE(msgBuf,'(A)')
     &   ' PFC_READPARMS: finished reading data.pfc'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )

C--   Close the open data file
      CLOSE(iUnit)

c==============the following used in cfc but not in hg
c      IF ( atmPFC_timeOffset .EQ. UNSET_RL ) THEN
C-    Assume that pfcTime=0 corresponds to the beginning of the 1rst record
C     time-period. This is consistent with 1rst record value = time-averaged
C     atmos-PFC over time period: pfcTime= 0 to pfcTime= 1 x atmPFC_recSepTime
C-    To reproduce previous time interpolation (formerly in pfc02_forcing.F),
C     set default timeOffset according to PTRACERS_Iter0 :
c        atmPFC_timeOffset = atmPFC_recSepTime
c     &                     - deltaTClock*PTRACERS_Iter0
c      ENDIF

c      PFC_monFreq = 0. _d 0
c      IF ( debugLevel.GE.debLevB ) PFC_monFreq = PTRACERS_monitorFreq
c==============


      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#endif /* ALLOW_PFC */

      RETURN
      END
!EOC
