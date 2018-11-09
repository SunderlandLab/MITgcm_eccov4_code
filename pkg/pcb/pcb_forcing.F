!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pcb_forcing.F
!
! !DESCRIPTION: Calculate the changes of PCBa and PCBP through air-sea  
!  fluxes and through rivers          
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "PCB_OPTIONS.h"

      SUBROUTINE PCB_FORCING( Ptr,
     &                       Gchem,
     &                       bi,bj,imin,imax,jmin,jmax,
     &                       myIter,myTime,myThid)

      IMPLICIT NONE
!
! !USES:
!
      ! Global variables
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PCB_SIZE.h"
#include "PCB.h"
!
! !INPUT PARAMETERS: 
!
      _RL myTime
      INTEGER myThid
      INTEGER myIter
      INTEGER i,j,bi, bj, imin, imax, jmin, jmax
!
! !INPUT/OUTPUT PARAMETERS: 
!
      _RL Ptr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      _RL Gchem(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
!
! !REMARKS:
!
! !REVISION HISTORY:
!  modified for external_forcing_DIC.F  August 1999
!  modified swd Oct 01 and Feb 02, for use as package for 
!    c40_patch1
!  modified to use with c44 and ptracers: swd May 2002
!  modified to have carbonate and biological influences: 
!    swd June 2002
!  modified for hg: swd Sep 2003
!  modifed for PFCs  xz
!  03 Feb 2015 - H. Amos - Update river inputs of PCBs to 
!                          mimic Hg river code.
!  29 Feb 2016 - R. Yantosca - Added ProTex headers
!  01 Mar 2016 - R. Yantosca - Removed commented-out code
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PCB
      INTEGER nriver     

      _RL  SURPCBa(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  SURPCBP(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
 
      ! Mimic Hg for PCB river inputs
      ! (H. Amos, 03 Feb 2015) 
#ifdef PCBAREFRACTORY
      _RL  SURPCBaR(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#endif /* PCBAREFRACTORY */
#ifdef PCBPREFRACTORY
      _RL  SURPCBpR(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nrmax) 
#endif /* PCBPREFRACTORY */

         ! Initialize
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           SURPCBa(i,j)=0.d0
           SURPCBp(i,j)=0.d0

           ! Mimic Hg for PCB river inputs
           ! (H. Amos, 03 Feb 2015) 
#ifdef PCBAREFRACTORY
           SURPCBar(i,j)=0.d0
#endif /* PCBAREFRACTORY */
#ifdef PCBPREFRACTORY
           DO nriver=1,nrmax
            SURPCBpR(i,j,nriver)=0.d0
           ENDDO 
#endif /* PCBPREFRACTORY */
          ENDDO
         ENDDO

       ! PCBa air-sea interaction
       CALL pcba_SURFFORCING( 
     &                   Ptr(1-OLx,1-OLy,1,1,1,1),
     &                   SURPCBa,
     &                    bi,bj,imin,imax,jmin,jmax,
     &                    myIter,myTime,myThid)

       ! PCBp river runoff
       CALL pcbp_SURFFORCING( 
     &                      Ptr(1-OLx,1-OLy,1,1,1,1), 
     &                      SURPCBp,
     &                    bi,bj,imin,imax,jmin,jmax,
     &                    myIter,myTime,myThid)

       ! Refractory river runoff
#ifdef PCBAREFRACTORY     
       CALL PCBaR_SURFFORCING( Ptr(1-OLx,1-OLy,1,1,1,1), SURPCBaR,
     &                    bi,bj,imin,imax,jmin,jmax,
     &                    myIter,myTime,myThid)
#endif /* PCBAREFRACTORY */
#ifdef PCBPREFRACTORY     
       CALL PCBpR_SURFFORCING( Ptr(1-OLx,1-OLy,1,1,1,1), SURPCBpR,
     &                    bi,bj,imin,imax,jmin,jmax,
     &                    myIter,myTime,myThid) 
#endif /* PCBPREFRACTORY */


       DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
          Gchem(i,j,1,bi,bj,ipcba)=Gchem(i,j,1,bi,bj,ipcba)+SURPCBa(i,j)
          Gchem(i,j,1,bi,bj,ipcbp)=Gchem(i,j,1,bi,bj,ipcbp)+SURPCBp(i,j)

          ! Mimic Hg for PCB river inputs
          ! (H. Amos, 03 Feb 2015)
#ifdef PCBAREFRACTORY
          Gchem(i,j,1,bi,bj,ipcbaR)=
     &      Gchem(i,j,1,bi,bj,ipcbaR)+SURPCBaR(i,j)
#endif /* PCBAREFRACTORY */
#ifdef PCBPREFRACTORY
          DO nriver=1,nrmax
             Gchem(i,j,1,bi,bj,ipcbPR+nriver-1)=
     &        Gchem(i,j,1,bi,bj,ipcbPR+nriver-1)+SURPCBpR(i,j,nriver)
          ENDDO
#endif /* PCBPREFRACTORY */          

       ENDDO
      ENDDO

#endif /* ALLOW_PCB      */
#endif /* ALLOW_PTRACERS */

      RETURN
       END
!EOC
