!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: hgp_surfforcing.F
!
! !DESCRIPTION: Handles surface forcing of the HgP species.
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "HG_OPTIONS.h"

      SUBROUTINE hgp_SURFFORCING( PTR, surfhgp,
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
      _RL  surfhgp(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
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
      _RL  Fluxhgp(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  hgprivero
#ifdef ALLOW_DIAGNOSTICS      
      _RL  hgpriver_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      CHARACTER*8 diagname
#endif /* ALLOW_DIAGNOSTICS  */     

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_HG

C     == Local variables ==
C I, J - Loop counters
      INTEGER I,J

C Initialize diagnostics
#ifdef ALLOW_DIAGNOSTICS
      DO j=jMin,jMax
        DO i=iMin,iMax
            hgpriver_diag(i,j)=0d0
        ENDDO
      ENDDO
#endif /* ALLOW_DIAGNOSTICS  */
      
      DO j=jMin,jMax
        DO i=iMin,iMax
          IF (HFacC(i,j,1,bi,bj).GT.0.5) THEN
C Localize variabls
C            kg/s -> mol/m2/s
             hgprivero = Hgpriver(i,j,bi,bj)*recip_rA(i,j,bi,bj)*4.985
                           
             Fluxhgp(i,j)=hgprivero
#ifdef ALLOW_DIAGNOSTICS             
             hgpriver_diag(i,j)=FluxHGp(i,j)
#endif /* ALLOW_DIAGNOSTICS  */            
          ENDIF
        ENDDO
      ENDDO

C update surface tendencies
        DO j=jMin,jMax
          DO i=iMin,iMax
           SURFHGP(i,j)=
     &         HFacC(i,j,1,bi,bj)*Fluxhgp(i,j)*recip_drF(1)
          ENDDO
         ENDDO
         
c write diagnostics
#ifdef ALLOW_DIAGNOSTICS
        WRITE(diagname,'(A8)') 'HGPRIVER'
        CALL DIAGNOSTICS_FILL( hgpriver_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )         
#endif /* ALLOW_DIAGNOSTICS  */
         
#endif /* ALLOW_HG           */   
#endif /* ALLOW_PTRACERS     */   

        RETURN
        END
!EOC
