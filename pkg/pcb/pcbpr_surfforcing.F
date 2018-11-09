!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pcbpr_surfforcing.F
!
! !DESCRIPTION: External inputs of refractory PCBp from rivers. 
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "PCB_OPTIONS.h"
#ifdef PCBPREFRACTORY

      SUBROUTINE PCBpR_SURFFORCING( Ptr, surfpcbpr,
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
#include "PCB_SIZE.h"
#include "PCB.h"
!
! !INPUT PARAMETERS: 
!
      INTEGER myIter, myThid
      _RL myTime
      INTEGER iMin,iMax,jMin,jMax, bi, bj
!
! !INPUT/OUTPUT PARAMETERS: 
!
      _RL  Ptr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      _RL  surfpcbpr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nrmax)
!
! !REMARKS:
!
! !REVISION HISTORY:
!  04 Feb 2015 - H. Amos     - Copied from Yanxu Zhang's Hg
!                              code and adapted for PCBp. 
!  29 Feb 2016 - R. Yantosca - Added ProTex headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!      
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PCB
      _RL  Fluxpcbpr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nrmax)
      _RL  pcbprrivero(nrmax)
#ifdef ALLOW_DIAGNOSTICS
      _RL  pcbprriver_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      CHARACTER*8 diagname

      ! Local variables
      INTEGER I,J,nriver  ! I, J - loop counters

#ifdef ALLOW_DIAGNOSTICS
      !----------------------------------------------------------------
      ! Initialize diagnostics
      !----------------------------------------------------------------
      DO j=jMin,jMax
        DO i=iMin,iMax
           pcbprriver_diag(i,j)=0d0
        ENDDO
      ENDDO
#endif /* ALLOW_DIAGNOSTICS */
      
      DO j=jMin,jMax
        DO i=iMin,iMax
          IF (HFacC(i,j,1,bi,bj).GT.0.5) THEN
            !----------------------------------------------------------
            ! Localize variables
            !----------------------------------------------------------

            ! Inputs from rivers
            ! Convert kg/s -> mol/m2/s
            DO nriver=1,nrmax
             pcbprrivero(nriver) = pcbprriver(i,j,bi,bj,nriver)
     &           *recip_rA(i,j,bi,bj)*1000.d0/Mb
                           
             Fluxpcbpr(i,j,nriver)=pcbprrivero(nriver)
#ifdef ALLOW_DIAGNOSTICS             
             pcbprriver_diag(i,j)=pcbprriver_diag(i,j)+
     &           Fluxpcbpr(i,j,nriver)
#endif /* ALLOW_DIAGNOSTICS */             
            ENDDO
          ENDIF
        ENDDO
      ENDDO

      !----------------------------------------------------------------
      ! update surface tendencies
      !----------------------------------------------------------------
      DO j=jMin,jMax
        DO i=iMin,iMax
          DO nriver=1,nrmax
            SURFpcbpr(i,j,nriver)=
     &         HFacC(i,j,1,bi,bj)*Fluxpcbpr(i,j,nriver)*recip_drF(1)
          ENDDO
        ENDDO
      ENDDO
         
#ifdef ALLOW_DIAGNOSTICS
      !----------------------------------------------------------------
      ! write diagnostics
      !----------------------------------------------------------------
      WRITE(diagname,'(A8)') 'PCBPRRIV'
      CALL DIAGNOSTICS_FILL( pcbprriver_diag(1-Olx,1-Oly), diagname,
     &                       0,1,2,bi,bj,myThid )         
#endif /* ALLOW_DIAGNOSTICS */
         
#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_PCB         */
#endif /* ALLOW_PTRACERS    */
      RETURN
      END
#endif /* PCBPREFRACTORY    */
!EOC

