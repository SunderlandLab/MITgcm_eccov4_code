!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pcb_readparms.F
!
! !DESCRIPTION: Initialize and read PCB package parameters.
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "PCB_OPTIONS.h"

      SUBROUTINE PCB_READPARMS( myThid )
!
! !USES:
!
      IMPLICIT NONE

      ! Global variables
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "PCB_SIZE.h"
#include "PCB.h"
!
! !INPUT PARAMETERS: 
!
      ! myThid :: My Thread Id. number
      INTEGER myThid
!
! !REMARKS:
!
! !REVISION HISTORY:
!  03 Feb 2015 - H. Amos     - Update river inputs of PCBs to
!                              mimic Yanxu Zhang's Hg code.  
!  01 Mar 2016 - R. Yantosca - Added ProTex headers
!  01 Mar 2016 - R. Yantosca - Removed commented-out code
!  01 Mar 2016 - R. Yantosca - Removed #ifdef FOODW_MODE
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#ifdef ALLOW_PCB

      ! Local variables
      !   msgBuf    :: Informational/error message buffer
      !   iUnit     :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit

      NAMELIST /PCB_FORCING/
          ! %%% FILENAMES %%%
     &          PCB_windFile,      
     &          PCB_iceFile,
     &          PCBp_depoFile,     
     &          PCBa_depoFile,
     &          PCBa_concFile,    
     &          PCBp_riverFile,
     &          PCBa_riverFile,
#ifdef PCBAREFRACTORY
     &          PCBAR_riverFile,
#endif /* PCBAREFRACTORY */
#ifdef PCBPREFRACTORY
     &          PCBPR_riverFile,
#endif /* PCBPREFRACTORY */    
     &          chl_file, 
     &          npp_file,
     &          doc_file, 
     &          poc_file,
     &          Fpoc_file, 
     &          Fremin_file,
          ! %%% TIMESTEPS %%%
     &          PCB_forcingPeriod, 
     &          PCB_forcingCycle,
     &          PCB_chemPeriod,    
     &          PCB_chemCycle,
     &          PCB_emisPeriod,    
     &          PCB_emisCycle

     
      _BEGIN_MASTER(myThid)

      ! Set all files in the common block /PCB_FILENAMES/ (located in
      ! include file PFC.h) to ' ', which will skip them from being 
      ! read in from disk, unless they are explicitly specified in the 
      ! data.pcb namelist file. (bmy, 3/4/16)
      PCB_windFile      = ' '
      PCB_iceFile       = ' '
      PCBp_depoFile     = ' '
      PCBa_depoFile     = ' '
      PCBa_concFile     = ' '
      PCBp_riverFile    = ' '      
      PCBa_riverFile    = ' '      
#ifdef PCBAREFRACTORY
      PCBaR_riverFile   = ' '
#endif /* PCBAREFRACTORY */
#ifdef PCBPREFRACTORY
      PCBpR_riverFile   = ' '
#endif /* PCBPREFRACTORY */    
      chl_file          = ' '
      npp_file          = ' '
      doc_file          = ' '
      poc_file          = ' '
      Fpoc_file         = ' '
      Fremin_file       = ' '
      radsw_File        = ' '
      
      ! Default periodic forcing to same as for physics
      PCB_forcingPeriod = externForcingPeriod
      PCB_forcingCycle  = externForcingCycle
      PCB_emisPeriod    = externForcingPeriod
      PCB_emisCycle     = externForcingCycle       

      ! Added based on hg code default periodic forcing 
      ! to same as for physics (Xianming)
      PCB_chemPeriod    = externForcingPeriod
      PCB_chemCycle     = externForcingCycle       

      WRITE(msgBuf,'(A)') ' PCB_READPARMS: opening data.pcb'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )

      CALL OPEN_COPY_DATA_FILE( 'data.pcb', 'PCB_READPARMS',
     O                          iUnit, myThid )

      !----------------------------------------------------------------
      ! Read parameters from open data file
      !----------------------------------------------------------------

      ! Forcing filenames and parameters
      READ(UNIT=iUnit,NML=PCB_FORCING)

      WRITE(msgBuf,'(A)')
     &   ' PCB_READPARMS: finished reading data.pcb'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )

      ! Close the open data file
      CLOSE(iUnit)

      _END_MASTER(myThid)

       ! Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_PCB */

      RETURN
      END
!EOC
