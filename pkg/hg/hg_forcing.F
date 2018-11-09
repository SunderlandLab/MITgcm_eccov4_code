!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: hg_forcing.F
!
! !DESCRIPTION: Driver routine which applies the various individual
!  forcings to Hg tracers through air-sea fluxes.
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "HG_OPTIONS.h"

      SUBROUTINE HG_FORCING( Ptr,
     &                       Gchem,
     &                       bi,bj,imin,imax,jmin,jmax,
     &                       myIter,myTime,myThid)

!
! !USES:
!
      IMPLICIT NONE

      ! Global variables
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "HG_SIZE.h"
#include "HG.h"
!
! !INPUT PARAMETERS: 
!
      INTEGER myIter
      _RL myTime
      INTEGER myThid
      INTEGER i,j,bi, bj, imin, imax, jmin, jmax
!
! !INPUT/OUTPUT PARAMETERS: 
!
      _RL  Ptr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      _RL  Gchem(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
!
! !REMARKS:
!     /==========================================================\
!     | SUBROUTINE HG_FORCING                                   |
!     | o Calculate the changes to Hg0 and Hg2                  |
!     |   through air-sea  fluxes                               |
!     |==========================================================|
!
! !REVISION HISTORY:
!  modified for external_forcing_DIC.F  August 1999
!  modified swd Oct 01 and Feb 02, for use as package for c40_patch1
!  modified to use with c44 and ptracers: swd May 2002
!  modified to have carbonate and biological influences: swd June 2002
!  modified for hg: swd Sep 2003
!  29 Feb 2016 - R. Yantosca - Added ProTex headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER nriver     

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_HG
C     == Local variables ==
      _RL  SURHg0(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  SURHg2(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  SURHgP(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#ifdef HG2REFRACTORY
      _RL  SURHg2R(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#endif /* HG2REFRACTORY */
#ifdef HGPREFRACTORY
      _RL  SURHgPR(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nrmax) 
#endif /* HGPREFRACTORY */           
#ifdef METHY_MODE      
      _RL  SURDMHg(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#endif /* METHY_MODE    */

         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           SURhg0(i,j)=0.d0
           SURhg2(i,j)=0.d0
           SURhgp(i,j)=0.d0
#ifdef HG2REFRACTORY
           SURhg2r(i,j)=0.d0
#endif /* HG2REFRACTORY */
#ifdef HGPREFRACTORY
           DO nriver=1,nrmax
            SURhgpr(i,j,nriver)=0.d0
           ENDDO 
#endif /* HGPREFRACTORY */          
#ifdef METHY_MODE
           SURDMHg(i,j)=0.d0
#endif /* METHY_MODE    */
          ENDDO
         ENDDO

C hg0 air-sea interaction
       CALL hg0_SURFFORCING( Ptr(1-OLx,1-OLy,1,1,1,1), SURhg0,
     &                    bi,bj,imin,imax,jmin,jmax,
     &                    myIter,myTime,myThid)

C hg2 air-sea interaction and river runoff
       CALL hg2_SURFFORCING( Ptr(1-OLx,1-OLy,1,1,1,1), SURhg2,
     &                    bi,bj,imin,imax,jmin,jmax,
     &                    myIter,myTime,myThid)

C hgp river runoff
       CALL hgp_SURFFORCING( Ptr(1-OLx,1-OLy,1,1,1,1), SURhgp,
     &                    bi,bj,imin,imax,jmin,jmax,
     &                    myIter,myTime,myThid)

C hgr river runoff
#ifdef HG2REFRACTORY     
       CALL hg2r_SURFFORCING( Ptr(1-OLx,1-OLy,1,1,1,1), SURhg2r,
     &                    bi,bj,imin,imax,jmin,jmax,
     &                    myIter,myTime,myThid)
#endif /* HG2REFRACTORY */
#ifdef HGPREFRACTORY     
       CALL hgpr_SURFFORCING( Ptr(1-OLx,1-OLy,1,1,1,1), SURhgpr,
     &                    bi,bj,imin,imax,jmin,jmax,
     &                    myIter,myTime,myThid) 
#endif /* HGPREFRACTORY */          

C DMHg air-sea exchange
#ifdef METHY_MODE
       CALL dmhg_SURFFORCING( Ptr(1-OLx,1-OLy,1,1,1,1), SURdmhg,
     &                    bi,bj,imin,imax,jmin,jmax,
     &                    myIter,myTime,myThid)
#endif /* METHY_MODE */
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
c           Gchem(i,j,1,bi,bj,ihg0)=Gchem(i,j,1,bi,bj,ihg0)+SURhg0(i,j)
c           Gchem(i,j,1,bi,bj,ihg2)=Gchem(i,j,1,bi,bj,ihg2)+SURhg2(i,j)
c           Gchem(i,j,1,bi,bj,ihgp)=Gchem(i,j,1,bi,bj,ihgp)+SURhgp(i,j)
            Ptr(i,j,1,bi,bj,ihg0) = Ptr(i,j,1,bi,bj,ihg0)
     &                           + SURhg0(i,j) * deltaTclock
            Ptr(i,j,1,bi,bj,ihg2) = Ptr(i,j,1,bi,bj,ihg2)
     &                           + SURhg2(i,j) * deltaTclock
            Ptr(i,j,1,bi,bj,ihgp) = Ptr(i,j,1,bi,bj,ihgp)
     &                           + SURhgp(i,j) * deltaTclock
#ifdef HG2REFRACTORY
c            Gchem(i,j,1,bi,bj,ihg2R)=
c     &      Gchem(i,j,1,bi,bj,ihg2R)+SURhg2r(i,j)
            Ptr(i,j,1,bi,bj,ihg2R) = Ptr(i,j,1,bi,bj,ihg2R)
     &                           + SURhg2r(i,j) * deltaTclock
#endif /* HG2REFRACTORY */
#ifdef HGPREFRACTORY
           DO nriver=1,nrmax
c            Gchem(i,j,1,bi,bj,ihgPR+nriver-1)=
c     &      Gchem(i,j,1,bi,bj,ihgPR+nriver-1)+SURhgpr(i,j,nriver)
            Ptr(i,j,1,bi,bj,ihgPR+nriver-1) = 
     &                           Ptr(i,j,1,bi,bj,ihgPR+nriver-1)
     &                           + SURhgpr(i,j,nriver) * deltaTclock
           ENDDO
#endif /* HGPREFRACTORY */         
#ifdef METHY_MODE
c           Gchem(i,j,1,bi,bj,idmhg)=Gchem(i,j,1,bi,bj,idmhg)
c     &          +SURdmhg(i,j)
            Ptr(i,j,1,bi,bj,idmhg) = Ptr(i,j,1,bi,bj,idmhg)
     &                           + SURdmhg(i,j) * deltaTclock
#endif /* METHY_MODE */
          ENDDO
         ENDDO

#endif /* ALLOW_HG       */
#endif /* ALLOW_PTRACERS */

       RETURN
       END
