!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: hgpr_surfforcing.F 
!
! !DESCRIPTION: Computes surface forcing of the riverine HgP species.
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "HG_OPTIONS.h"
#ifdef HGPREFRACTORY

      SUBROUTINE hgpr_SURFFORCING( Ptr, surfhgpr,
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
#include "HG_SIZE.h"
#include "HG.h"
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
      _RL  surfhgpr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nrmax)
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
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_HG
      _RL  Fluxhgpr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nrmax)
      _RL  hgprrivero(nrmax)
#ifdef ALLOW_DIAGNOSTICS
      _RL  hgprriver_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      CHARACTER*8 diagname

C     == Local variables ==
C I, J - Loop counters
      INTEGER I,J,nriver

C Initialize diagnostics
#ifdef ALLOW_DIAGNOSTICS
      DO j=jMin,jMax
        DO i=iMin,iMax
            hgprriver_diag(i,j)=0d0
        ENDDO
      ENDDO
#endif  /* ALLOW_DIAGNOSTICS */
      
      DO j=jMin,jMax
        DO i=iMin,iMax
          IF (HFacC(i,j,1,bi,bj).GT.0.5) THEN
C Localize variabls
C            kg/s -> mol/m2/s
            DO nriver=1,nrmax
             hgprrivero(nriver) = hgprriver(i,j,bi,bj,nriver)
     &           *recip_rA(i,j,bi,bj)*4.985
                           
             Fluxhgpr(i,j,nriver)=hgprrivero(nriver)
#ifdef ALLOW_DIAGNOSTICS             
             hgprriver_diag(i,j)=hgprriver_diag(i,j)+
     &           Fluxhgpr(i,j,nriver)
#endif  /* ALLOW_DIAGNOSTICS */            
            ENDDO
          ENDIF
        ENDDO
      ENDDO

C update surface tendencies
        DO j=jMin,jMax
          DO i=iMin,iMax
           DO nriver=1,nrmax
            SURFhgpr(i,j,nriver)=
     &         HFacC(i,j,1,bi,bj)*Fluxhgpr(i,j,nriver)*recip_drF(1)
           ENDDO
          ENDDO
         ENDDO
         
c write diagnostics
#ifdef ALLOW_DIAGNOSTICS
        WRITE(diagname,'(A8)') 'HGPRRIVE'
        CALL DIAGNOSTICS_FILL( hgprriver_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )         
#endif /* ALLOW_DIAGNOSTICS */
         
#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_HG          */
#endif /* ALLOW_PTRACERS    */
        RETURN
        END
#endif /* HGPREFRACTORY     */
!EOC
