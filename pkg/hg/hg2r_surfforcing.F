!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: hg2r_surfforcing.F
!
! !DESCRIPTION: Surface forcing of the riverine Hg2 tracer.
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "HG_OPTIONS.h"
#ifdef HG2REFRACTORY

      SUBROUTINE hg2r_SURFFORCING( Ptr, surfhg2r,
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
      _RL  surfhg2r(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
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
      
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_HG
      _RL  Fluxhg2r(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  hg2depoo,hg2rrivero
#ifdef ALLOW_DIAGNOSTICS
      _RL  hg2rriver_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      CHARACTER*8 diagname

C     == Local variables ==
C I, J - Loop counters
      INTEGER I,J

C Initialize diagnostics
#ifdef ALLOW_DIAGNOSTICS
      DO j=jMin,jMax
        DO i=iMin,iMax
            hg2rriver_diag(i,j)=0d0
        ENDDO
      ENDDO
#endif /* ALLOW_DIAGNOSTICS */
      
      DO j=jMin,jMax
        DO i=iMin,iMax
          IF (HFacC(i,j,1,bi,bj).GT.0.5) THEN
C Localize variabls
C            kg/s -> mol/m2/s
C            depo  : kg/m2/s -> mol/m2/s
             hg2depoo = Hg2Depo(i,j,bi,bj)*4.985
C            Assume part of the deposition is refractory
     &                 * f_depo_refrac
C            river : kg/s    -> mol/m2/s
             hg2rrivero = hg2rriver(i,j,bi,bj)
     &           *recip_rA(i,j,bi,bj)*4.985
                           
             Fluxhg2r(i,j)=hg2depoo+hg2rrivero
#ifdef ALLOW_DIAGNOSTICS             
             hg2rriver_diag(i,j)=hg2rriver_diag(i,j)+
     &           Fluxhg2r(i,j)
#endif /* ALLOW_DIAGNOSTICS */      
          ENDIF
        ENDDO
      ENDDO

C update surface tendencies
        DO j=jMin,jMax
          DO i=iMin,iMax
            SURFhg2r(i,j)=
     &         HFacC(i,j,1,bi,bj)*Fluxhg2r(i,j)*recip_drF(1)
          ENDDO
         ENDDO
         
c write diagnostics
#ifdef ALLOW_DIAGNOSTICS
        WRITE(diagname,'(A8)') 'HG2RRIVE'
        CALL DIAGNOSTICS_FILL( hg2rriver_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )         
#endif /* ALLOW_DIAGNOSTICS */
         
#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_HG          */
#endif /* ALLOW_PTRACERS    */
        RETURN
        END
#endif /* HG2REFRACTORY */

