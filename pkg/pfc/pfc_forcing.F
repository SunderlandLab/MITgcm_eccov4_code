!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pfc_forcing
!
! !DESCRIPTION: Applies forcings (from sea-air fluxes) to the PFC tracers.
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "PFC_OPTIONS.h"

      SUBROUTINE PFC_FORCING( Ptr,
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
#include "PFC_SIZE.h"
#include "PFC.h"
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
!     | SUBROUTINE PFC_FORCING                                   |
!     | o Calculate the changes of PFCi, PFCn, PFCP,             |
!     |   through air-sea fluxes                                 |
!     |==========================================================|
!
! !REVISION HISTORY:
!  modified for external_forcing_DIC.F  August 1999
!  modified swd Oct 01 and Feb 02, for use as package for c40_patch1
!  modified to use with c44 and ptracers: swd May 2002
!  modified to have carbonate and biological influences: swd June 2002
!  modified for hg: swd Sep 2003
!  modified for PFCs xz 
!  modified for PFCs and checkpoint65r Bob Y. (2/26/16)
!  29 Feb 2016 - R. Yantosca - Added ProTex headers; updated comments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PFC
      _RL  SURPFCn(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  SURPFCi(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  SURPFCP(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      DO j=1-OLy,sNy+OLy
      DO i=1-OLx,sNx+OLx
         SURPFCn(i,j)=0.d0
         SURPFCi(i,j)=0.d0
         SURPFCp(i,j)=0.d0
      ENDDO
      ENDDO

       !============================================================
       ! %%%%%            PFCn air-sea interaction            %%%%%
       !============================================================
      CALL pfcn_SURFFORCING( 
     I                   Ptr(1-OLx,1-OLy,1,1,1,1),
     O                   SURPFCn,
     I                   bi,bj,imin,imax,jmin,jmax,
     I                   myIter,myTime,myThid)

       !============================================================
       ! %%%%%   PFCi air-sea interaction and river runoff   %%%%%
       !============================================================
       CALL pfci_SURFFORCING(
     I                    Ptr(1-OLx,1-OLy,1,1,1,1),
     O                    SURPFCi,
     I                    bi,bj,imin,imax,jmin,jmax,
     I                    myIter,myTime,myThid)


       !============================================================
       ! %%%%%              PFCp river runoff                 %%%%%
       !============================================================
       CALL pfcp_SURFFORCING( 
     I                    Ptr(1-OLx,1-OLy,1,1,1,1), 
     O                    SURPFCp,
     I                    bi,bj,imin,imax,jmin,jmax,
     I                    myIter,myTime,myThid)

         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
!#############################################################################
!### COMMENT OUT ORIGINAL CODE TO MATCH WHAT WE DO IN pkg/hg/hg_forcing.F
!### (bmy, 2/26/16)
!          Gchem(i,j,1,bi,bj,ipfcn)=Gchem(i,j,1,bi,bj,ipfcn)+SURPFCn(i,j)
!          Gchem(i,j,1,bi,bj,ipfci)=Gchem(i,j,1,bi,bj,ipfci)+SURPFCi(i,j)
!          Gchem(i,j,1,bi,bj,ipfcp)=Gchem(i,j,1,bi,bj,ipfcp)+SURPFCp(i,j)
!#############################################################################

            ! Update PFCn  
            Ptr(i,j,1,bi,bj,iPFCn) = Ptr(i,j,1,bi,bj,iPFCn)
     &                             + SURPFCn(i,j) * deltaTclock

            ! Update PFCi
            Ptr(i,j,1,bi,bj,iPFCi) = Ptr(i,j,1,bi,bj,iPFCi)
     &                             + SURPFCi(i,j) * deltaTclock

            ! Update PFCp
            Ptr(i,j,1,bi,bj,iPFCp) = Ptr(i,j,1,bi,bj,iPFCp)
     &                             + SURPFCp(i,j) * deltaTclock

          ENDDO
         ENDDO

#endif /* ALLOW_PFC      */
#endif /* ALLOW_PTRACERS */

       RETURN
       END
!EOC
