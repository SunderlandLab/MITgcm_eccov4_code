!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: hg2_surfforcing.F
!
! !DESCRIPTION: Surface forcing of the Hg2 tracer.
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "HG_OPTIONS.h"

      SUBROUTINE hg2_SURFFORCING( PTR, surfhg2,
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
      _RL  surfhg2(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
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
      _RL  FluxHG2(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  hg2depoo,hg2rivero

#ifdef ALLOW_DIAGNOSTICS      
      _RL  hg2depo_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  hg2river_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      CHARACTER*8 diagname
#endif /* ALLOW_DIAGNOSTICS */   

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_HG

C     == Local variables ==
C I, J - Loop counters
      INTEGER I,J

C Initialize diagnostics
#ifdef ALLOW_DIAGNOSTICS
      DO j=jMin,jMax
        DO i=iMin,iMax
            hg2depo_diag(i,j)=0d0
            hg2river_diag(i,j)=0d0
        ENDDO
      ENDDO
#endif /* ALLOW_DIAGNOSTICS */
      DO j=jMin,jMax
        DO i=iMin,iMax
          IF (HFacC(i,j,1,bi,bj).GT.0.5) THEN
C Localize variabls
C            depo  : kg/m2/s -> mol/m2/s
             hg2depoo = Hg2Depo(i,j,bi,bj)*4.985
#ifdef HG2REFRACTORY
C            Assume part of the deposition is refractory
     &                 * (1.0 _d 0 - f_depo_refrac) 
#endif /* HG2REFRACTORY */
C            river : kg/s    -> mol/m2/s
             hg2rivero = Hg2River(i,j,bi,bj)*recip_rA(i,j,bi,bj)*4.985
                           
             FluxHG2(i,j)=hg2depoo+hg2rivero
#ifdef ALLOW_DIAGNOSTICS             
             hg2depo_diag(i,j)=hg2depoo
             hg2river_diag(i,j)=hg2rivero
#endif /* ALLOW_DIAGNOSTICS */          
          ENDIF
        ENDDO
      ENDDO

C update surface tendencies
        DO j=jMin,jMax
          DO i=iMin,iMax
           SURFHG2(i,j)=
     &         HFacC(i,j,1,bi,bj)*FluxHG2(i,j)*recip_drF(1)
          ENDDO
         ENDDO
         
c write diagnostics
#ifdef ALLOW_DIAGNOSTICS
        WRITE(diagname,'(A8)') 'HG2DEPO '
        CALL DIAGNOSTICS_FILL( hg2depo_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )         
        WRITE(diagname,'(A8)') 'HG2RIVER'
        CALL DIAGNOSTICS_FILL( hg2river_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )         
#endif /* ALLOW_DIAGNOSTICS */
         
#endif /* ALLOW_HG          */
#endif /* ALLOW_PTRACERS    */

        RETURN
        END
!EOC
