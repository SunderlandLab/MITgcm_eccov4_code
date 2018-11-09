!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pcb_partition.F 
!
! !DESCRIPTION: Calculate paritioning between dissolved and particulate 
!  PCB phases in the ocean. 
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "PCB_OPTIONS.h"

      SUBROUTINE PCB_PARTITION( Ptr,
     &                          bi    , bj    ,
     &                          imin  , imax  , jmin  , jmax ,
     &                          myIter, myTime, myThid        )
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
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PCB_SIZE.h"
#include "PCB.h"
!
! !INPUT PARAMETERS: 
!
      INTEGER myIter
      _RL     myTime
      INTEGER myThid
      INTEGER i,j,k,bi, bj, imin, imax, jmin, jmax
!
! !INPUT/OUTPUT PARAMETERS: 
!
      _RL  Ptr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
!
! !REVISION HISTORY:
!  From Xianming Zhang:
!    Based on the code of Hg species in the ocean by Yanxu Zhang.
!    To calculate the partition of PFOS species in the ocean
!    Local equlibrium between the anionic, neutral, and 
!    particulate species are assumed to be reached instantaneously
!
!  From Helen Amos (06 Aug 2014):
!    PFOS code adopted for the simulation of PCBs in the ocean.
!    Comments added and obsolete/commented-out code deleted to 
!    clean things up.
!
!  06 Aug 2015 - H. Amos    -  PFOS code adopted for the simulation of 
!                              PCBs in the ocean. Comments added and 
!                              obsolete/commented-out code deleted to 
!                              clean things up.
!  06 Jul 2015 - H. Amos     - Move Koc and Kow to pcb_param.F
!  05 Oct 2015 - H. Amos     - Updated Kow-Koc relationship
!  01 Mar 2016 - R. Yantosca - Added ProTex headers
!  01 Mar 2016 - R. Yantosca - Deleted commented-out code
!  01 Mar 2016 - R. Yantosca - Added missing #endif for PART_EQ
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#ifdef ALLOW_DIAGNOSTICS      
      _RL     pcbapart_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL     pcbppart_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)            
      CHARACTER*8 diagname
#endif /* ALLOW_DIAGNOSTICS */
      _RL     TOT_PCB
      _RL     FRAC_TotDissolved
      _RL     FRAC_Dissolved
      _RL     FRAC_POCBound
      _RL     poco, doco
      _RL     pcba_old, pcbp_old, pcba_new, pcbp_new

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PCB

      !--------------------------------------------------------
      ! Initialize diagnostics
      !--------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS
      DO k=1,Nr
        DO j=jMin,jMax
          DO i=iMin,iMax
            pcbapart_diag(i,j,k) = 0d0  
            pcbppart_diag(i,j,k) = 0d0         
          ENDDO
        ENDDO
      ENDDO
#endif /* ALLOW_DIAGNOSTICS */

      !--------------------------------------------------------
      ! Start calculation
      !--------------------------------------------------------
      DO j=jMin,jMax
       DO i=iMin,iMax
        DO k=1,Nr
         
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN

           !===========================================================
           ! %%%%% Instantaneous equilibrium %%%%%
           !===========================================================

           ! POC and DOC concentration (units?)
           poco = max( poc(i,j,k,bi,bj), 0d0 ) 
           doco = max( doc(i,j,k,bi,bj), 0d0 )

           !-----------------------------------------------------------
           ! PCB concentrations before partitioning 
           !-----------------------------------------------------------
           pcba_old = max(Ptr(i,j,k,bi,bj,ipcba),0d0)
           pcbp_old = max(Ptr(i,j,k,bi,bj,ipcbp),0d0)
           TOT_PCB = pcba_old + pcbp_old

           IF (TOT_PCB .gt. 0. _d 0) THEN

             !---------------------------------------------------------
             ! Calculate fraction in the dissolved phase
             ! convert mmolC/m3 -> kgC/L
             !---------------------------------------------------------
             FRAC_Dissolved = 1. _d 0 
     &             /(1. _d 0 + Koc_pcb * poco 
     &             * 1.2 _d - 8 )  ! mmolC/m3 -> kgC/L

             !---------------------------------------------------------
             ! pH of seawater
             !
             ! Note: Currently not considering the real pH
             !       (xz)
             !---------------------------------------------------------
             pH(i,j,k,bi,bj) = 8. _d 0  
           
             !---------------------------------------------------------
             ! Particulate fraction (i.e., PCB bound to POC)
             !---------------------------------------------------------
             FRAC_POCBound = 1.d0 - FRAC_Dissolved

             !---------------------------------------------------------
             ! New PCB concentrations (mol/m3) in the dissolved
             ! and particle phases after consdering partitioning
             !---------------------------------------------------------
             pcba_new = TOT_PCB * FRAC_Dissolved
             pcbp_new = TOT_PCB * FRAC_POCBound

             !---------------------------------------------------------
             ! Update tracer array
             !---------------------------------------------------------
             Ptr(i,j,k,bi,bj,ipcba) = pcba_new  
C     &                                - pcba_deg * deltaTclock
             Ptr(i,j,k,bi,bj,ipcbp) = pcbp_new

#ifdef ALLOW_DIAGNOSTICS
             pcbapart_diag(i,j,k) = (pcba_old - pcba_new) / deltaTclock
             pcbppart_diag(i,j,k) = (pcbp_old - pcbp_new) / deltaTclock
#endif /* ALLOW_DIAGNOSTICS */            
           ENDIF
         
         ENDIF !(HFacC(i,j,k,bi,bj).gt.0. _d 0) 
        ENDDO   !k
       ENDDO    !i
      ENDDO     !j

      !--------------------------------------------------------
      ! Write diagnostics
      !--------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS   
      WRITE(diagname,'(A8)') 'PCBaPART'
      CALL DIAGNOSTICS_FILL( pcbapart_diag(1-Olx,1-Oly,1), diagname,
     &                       0,Nr,2,bi,bj,myThid )
     
      WRITE(diagname,'(A8)') 'PCBpPART'
      CALL DIAGNOSTICS_FILL( pcbppart_diag(1-Olx,1-Oly,1), diagname,
     &                       0,Nr,2,bi,bj,myThid )   
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_PCB */
#endif /* ALLOW_PTRACERS */

      RETURN
      END
!EOC
