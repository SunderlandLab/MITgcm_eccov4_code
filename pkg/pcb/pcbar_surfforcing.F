!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pcbar_surfforcing.F
!
! !DESCRIPTION: External inputs of refractory PCBa from rivers and 
!   atmospheric deposition
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "PCB_OPTIONS.h"
#ifdef PCBAREFRACTORY

      SUBROUTINE PCBaR_SURFFORCING( Ptr, surfpcbar,
     I           bi,bj,iMin,iMax,jMin,jMax,
     I           myIter, myTime, myThid )
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

! Note Yanxu has his HG_SIZE.h in a difference place than my 
! PCB_SIZE.h, so that might create problems at compile or 
! run time. Something to look out for. 
! (H. Amos, 04 Feb 2015)
#include "PCB_SIZE.h"
#include "PCB.h"
!
! !INPUT PARAMETERS: 
!
      ! bi, bj     :: current tile indices
      ! iMin,iMax  :: computation domain, 1rst index bounds
      ! jMin,jMax  :: computation domain, 2nd  index bounds
      ! myTime     :: current time in simulation
      ! myIter     :: current iteration number
      ! myThid     :: my Thread Id number

      INTEGER myIter, myThid
      _RL myTime
      INTEGER iMin,iMax,jMin,jMax, bi, bj
!
! !INPUT/OUTPUT PARAMETERS: 
!
      ! pTr        :: ocean PCB01 concentration
      ! SurfPCBar  :: atmospheric PCBar field
      _RL  Ptr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      _RL  surfpcbar(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
!
! !REMARKS:
!
! !REVISION HISTORY:
!  04 Feb 2015 - H. Amos     - Copied from Yanxu Zhang's Hg
!                              code and adapted for PCBa.    
!  01 Mar 2016 - R. Yantosca - Added ProTex Headers; updated comments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PCB
      _RL  Fluxpcbar(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  pcbadepoo,pcbarrivero
#ifdef ALLOW_DIAGNOSTICS
      _RL  pcbarriver_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      CHARACTER*8 diagname

      INTEGER I,J    ! loop counters

      !-----------------------------------------------------------
      ! Initialize diagnostics
      !-----------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS
      DO j=jMin,jMax
        DO i=iMin,iMax
            pcbarriver_diag(i,j)=0d0
        ENDDO
      ENDDO
#endif /* ALLOW_DIAGNOSTICS */
      
      DO j=jMin,jMax
        DO i=iMin,iMax
          IF (HFacC(i,j,1,bi,bj).GT.0.5) THEN
             !----------------------------------------------------
             ! Localize variables
             !----------------------------------------------------

             ! Atmospheric deposition of PCBa 
             ! Convert kg/m2/s -> mol/m2/s
             pcbadepoo = PCBaDepo(i,j,bi,bj)*1000.d0/Mb

             ! Assume part of the deposition is refractory
     &                 * f_depo_refrac

             ! Inputs of PCBa from rivers
             ! Convert kg/s -> mol/m2/s
             pcbarrivero = pcbarriver(i,j,bi,bj)
     &           *recip_rA(i,j,bi,bj)*1000.d0/Mb
                           
             Fluxpcbar(i,j)=pcbadepoo+pcbarrivero
#ifdef ALLOW_DIAGNOSTICS             
             pcbarriver_diag(i,j)=pcbarriver_diag(i,j)+
     &           Fluxpcbar(i,j)
#endif /* ALLOW_DIAGNOSTICS */
          ENDIF
        ENDDO
      ENDDO

        !---------------------------------------------------------
        ! update surface tendencies
        !---------------------------------------------------------
        DO j=jMin,jMax
          DO i=iMin,iMax
            SURFpcbar(i,j)=
     &         HFacC(i,j,1,bi,bj)*Fluxpcbar(i,j)*recip_drF(1)
          ENDDO
         ENDDO
         
        !---------------------------------------------------------
        ! write diagnostics
        !---------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS
        WRITE(diagname,'(A8)') 'PCBaRRIV'
        CALL DIAGNOSTICS_FILL( pcbarriver_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )         
#endif /* ALLOW_DIAGNOSTICS */
         
#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_PCB         */
#endif /* ALLOW_PTRACERS    */
        RETURN
        END
#endif /* PCBAREFRACTORY    */
!EOC
