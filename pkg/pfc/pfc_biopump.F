!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pfc_biopump.F
!
! !DESCRIPTION: To calculate the sinking of PFC along POM.
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "PFC_OPTIONS.h"

      SUBROUTINE PFC_BIOPUMP(  Ptr,Gchem,      
     &                        bi,bj,imin,imax,jmin,jmax,
     &                        myIter,myTime,myThid)
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
      INTEGER i,j,k,bi, bj, imin, imax, jmin, jmax      
!
! !INPUT/OUTPUT PARAMETERS: 
!
      _RL  Ptr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      _RL  Gchem(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)    
!
! !REMARKS:
!     /==========================================================\
!     | SUBROUTINE PFC_BIOPUMP                                   |
!     | o Calculate the changes to PFCP                          |
!     |   through particulate sinking                            |
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
      _RL  pfcpsink_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      CHARACTER*8 diagname
#endif /* ALLOW_DIAGNOSTICS */    
      _RL  Fout_PFCp(Nr)
      _RL  poco, fpoco
      _RL  tiny
      data tiny /1d-32/

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PFC
C     == Local variables == 

C Initialize diagnostics
#ifdef ALLOW_DIAGNOSTICS
      DO k=1,Nr
        DO j=jMin,jMax
          DO i=iMin,iMax
            pfcpsink_diag(i,j,k)=0d0
          ENDDO
        ENDDO
      ENDDO
#endif /* ALLOW_DIAGNOSTICS */     
      
      DO j=jMin,jMax
       DO i=iMin,iMax

C       Zero flux
        DO k=1,Nr
         Fout_PFCp(k) = 0. _d 0
        ENDDO

C       Flux calculation
        DO k=1,Nr
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN
           ! localize variables
           poco =  max(poc(i,j,k,bi,bj),0d0) !mmolC/m3
           fpoco = max(Fpoc(i,j,k,bi,bj), 0d0) !mmolC/m2/d
             
           Fout_PFCp(k) = max(PTR(i,j,k,bi,bj,ipfcp),0d0)
     &          * fpoco / ( poco + tiny )
     &          / 8.64 d 4
c       mol/m2/s = mol/m3 * (mmolC/m2/d) / (mmolC/m3) /(s/d)     
     
C Calculate diagnistics
#ifdef ALLOW_DIAGNOSTICS
           pfcpsink_diag(i,j,k) = Fout_PFCp(k)  !mol/m^2/s   
#endif /* ALLOW_DIAGNOSTICS */

         ENDIF
        ENDDO

C----&---1---------2---------3---------4---------5---------6---------7--    
C       Tendency calculation
        IF (HFacC(i,j,1,bi,bj).gt.0. _d 0) THEN
         Gchem(i,j,1,bi,bj,ipfcp)=Gchem(i,j,1,bi,bj,ipfcp) 
     &                          - Fout_PFCp(1) * recip_drF(1)
        ENDIF
  
        DO k=2,Nr
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN
          Gchem(i,j,k,bi,bj,ipfcp)=Gchem(i,j,k,bi,bj,ipfcp) + 
     &             (Fout_PFCp(k-1) - Fout_PFCp(k)) * recip_drF(k)
         ENDIF
        ENDDO

       ENDDO
      ENDDO

c write diagnostics
#ifdef ALLOW_DIAGNOSTICS
        WRITE(diagname,'(A8)') 'PFCpSINK'
        CALL DIAGNOSTICS_FILL( PFCpsink_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_PFC         */
#endif /* ALLOW_PTRACERS    */

       RETURN
       END
!EOC
