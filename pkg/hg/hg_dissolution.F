!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: hg_dissolution.F
!
! !DESCRIPTION: To calculate the chemistry of Hg species in the ocean
!  (Yanxu Zhang, 5/30/13)
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "HG_OPTIONS.h"

      SUBROUTINE HG_DISSOLUTION( Ptr,
     &                         Gchem,    
     &                         bi,bj,imin,imax,jmin,jmax,
     &                         myIter,myTime,myThid)
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
#include "GCHEM.h"
#include "HG_SIZE.h"
#include "HG.h"
!
! !INPUT PARAMETERS: 
!
      INTEGER myIter
      _RL myTime
      INTEGER myThid
      INTEGER i,j,k,l,bi, bj, imin, imax, jmin, jmax
!
! !INPUT/OUTPUT PARAMETERS: 
!
      _RL  Ptr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      _RL  Gchem(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)    
!
! !REMARKS:
!     /==========================================================\
!     | SUBROUTINE HG_DISSOLUTION                                |
!     | o Calculate the dissolution of HgR to Hg2                |
!     |==========================================================|
!
! !REVISION HISTORY:
!  29 Feb 2016 - R. Yantosca - Added ProTex headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#ifdef ALLOW_DIAGNOSTICS      
#ifdef HG2REFRACTORY       
      _RL  hg2rdiss_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif /*  HG2REFRACTORY    */
#ifdef HGPREFRACTORY       
      _RL  hgprdiss_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif /* HGPREFRACTORY     */               
      CHARACTER*8 diagname
#endif /* ALLOW_DIAGNOSTICS */
      
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_HG

C     -- Gross tendencies
#ifdef HG2REFRACTORY
      _RL hg2r_hg2
#endif /* HG2REFRACTORY */ 
#ifdef HGPREFRACTORY
      _RL hgpr_hg2
#endif /* HGPREFRACTORY */

C     == Parameters
#ifdef HG2REFRACTORY       
      _RL k_hg2r_hg2
      data k_hg2r_hg2 / 3.85 _d -8 /  !10%/ (30 day)
C      data k_hg2r_hg2 / 0.0 _d 0 /
#endif /* HG2REFRACTORY */
#ifdef HGPREFRACTORY       
      _RL k_hgpr_hg2
      data k_hgpr_hg2 / 3.85 _d -8 /  !10%/ (30 day)
      INTEGER nriver 
#endif /*  HGPREFRACTORY */          
      
C     == Start calculation
C Initialize diagnostics
#ifdef ALLOW_DIAGNOSTICS
      DO k=1,Nr
        DO j=jMin,jMax
          DO i=iMin,iMax  
#ifdef HG2REFRACTORY     
            hg2rdiss_diag(i,j,k)=0d0
#endif /* HG2REFRACTORY     */
#ifdef HGPREFRACTORY     
            hgprdiss_diag(i,j,k)=0d0
#endif /* HGPREFRACTORY     */
          ENDDO
        ENDDO
      ENDDO
#endif /* ALLOW_DIAGNOSTICS */

#if defined(HG2REFRACTORY) & defined(HGPREFRACTORY) & defined(HGPR_SLOWRELEASING)
      DO j=jMin,jMax
       DO i=iMin,iMax
        DO k=1,Nr
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN
C     == Tendency calcualtion                      
#ifdef HG2REFRACTORY
            hg2r_hg2 = 0d0
            hg2r_hg2 = hg2r_hg2 + k_hg2r_hg2
     &           * max( Ptr(i,j,k,bi,bj,ihg2r), 0d0)
#endif /* HG2REFRACTORY */
#ifdef HGPREFRACTORY
            hgpr_hg2 = 0d0
            DO nriver=1,nrmax
              hgpr_hg2 = hgpr_hg2 + k_hgpr_hg2 
     &           * max( Ptr(i,j,k,bi,bj,ihgpr+nriver-1), 0d0)
            ENDDO
#endif /* HGPREFRACTORY */
           
C Calculate diagnistics                                                  
#ifdef ALLOW_DIAGNOSTICS
#ifdef HG2REFRACTORY
            hg2rdiss_diag(i,j,k) = hg2r_hg2
#endif /* HG2REFRACTORY     */
#ifdef HGPREFRACTORY
            hgprdiss_diag(i,j,k) = hgpr_hg2
#endif /* HGPREFRACTORY     */                         
#endif /* ALLOW_DIAGNOSTICS */

C     == Accumulate tendency
#ifdef HG2REFRACTORY
c           Gchem(i,j,k,bi,bj,ihg2)=Gchem(i,j,k,bi,bj,ihg2)+hg2r_hg2
c           Gchem(i,j,k,bi,bj,ihg2r)=
c     &            Gchem(i,j,k,bi,bj,ihg2r) - k_hg2r_hg2 *
c     &            max(Ptr(i,j,k,bi,bj,ihg2r), 0d0)
            Ptr(i,j,k,bi,bj,ihg2) = Ptr(i,j,k,bi,bj,ihg2)
     &                            + hg2r_hg2 * deltaTclock
            Ptr(i,j,k,bi,bj,ihg2r) =
     &            Ptr(i,j,k,bi,bj,ihg2r) - k_hg2r_hg2 *
     &            max(Ptr(i,j,k,bi,bj,ihg2r), 0d0)
#endif /* HG2REFRACTORY */
#ifdef HGPREFRACTORY
c           Gchem(i,j,k,bi,bj,ihg2)=Gchem(i,j,k,bi,bj,ihg2)+hgpr_hg2
            Ptr(i,j,k,bi,bj,ihg2) = Ptr(i,j,k,bi,bj,ihg2)
     &                            + hgpr_hg2 * deltaTclock
           DO nriver=1,nrmax
c             Gchem(i,j,k,bi,bj,ihgpr+nriver-1)=
c     &            Gchem(i,j,k,bi,bj,ihgpr+nriver-1) - k_hgpr_hg2 * 
c     &            max(Ptr(i,j,k,bi,bj,ihgpr+nriver-1), 0d0)
             Ptr(i,j,k,bi,bj,ihgpr+nriver-1) =
     &            Ptr(i,j,k,bi,bj,ihgpr+nriver-1) - k_hgpr_hg2 *
     &            max(Ptr(i,j,k,bi,bj,ihgpr+nriver-1), 0d0)
           ENDDO
#endif /* HGPREFRACTORY */                  
         ENDIF
        ENDDO
       ENDDO
      ENDDO
#endif /* HG2REFRACTORY & HGPREFRACTORY & HGPR_SLOWRELEASING */
      
c write diagnostics
#ifdef ALLOW_DIAGNOSTICS
#ifdef HG2REFRACTORY
        WRITE(diagname,'(A8)') 'HG2RDISS'
        CALL DIAGNOSTICS_FILL( hg2rdiss_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#endif /* HG2REFRACTORY     */
#ifdef HGPREFRACTORY
        WRITE(diagname,'(A8)') 'HGPRDISS'
        CALL DIAGNOSTICS_FILL( hgprdiss_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#endif /* HGPREFRACTORY     */     
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_HG          */
#endif /* ALLOW_PTRACERS    */

       RETURN
       END
!EOC
