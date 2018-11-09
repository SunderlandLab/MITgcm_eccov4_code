!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: 
!
! !DESCRIPTION: To calculate the sinking of Hg species along POM,
!   also include the methylation and demethylation
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "HG_OPTIONS.h"

      SUBROUTINE HG_BIOPUMP(  Ptr,Gchem,      
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
#include "HG_SIZE.h"
#include "HG.h"
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
!     | SUBROUTINE HG_BIOPUMP                                   |
!     | o Calculate the changes to Hgp and MMHgp                |
!     |   through particulate sinking                           |
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
C     == Routine arguments ==
#ifdef ALLOW_DIAGNOSTICS      
      _RL  hgpsink_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      CHARACTER*8 diagname
#endif /* ALLOW_DIAGNOSTICS */
      _RL  Fout_Hgp(Nr)
#ifdef HGPREFRACTORY
      INTEGER nriver
      _RL  whgr(nrmax)
C     data whgr /100.0 _d 0, 10.0 _d 0, 1.0 _d 0/
C     data whgr /3.0 _d 0, 3.0 _d 0, 3.0 _d 0/        
      data whgr /3.0 _d 0 /
C     data whgr /30.0 _d 0 /
#if nrmax > 3
      What is the sinking speed of riverine particles?
#endif /* nrmax > 3 */             
      _RL  Fout_HgPR(Nr,nrmax)
#ifdef ALLOW_DIAGNOSTICS      
      _RL  hgprsink_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif /* ALLOW_DIAGNOSTICS */
#endif /* HGPREFRACTORY     */
#ifdef METHY_MODE
      INTEGER np,nz       
      _RL  Fout_MMHgp(Nr)
      _RL  Fout_MMHgphy(Nr,npmax)
      _RL  wphy(npmax)  ! sinking speed of phy, m/d
      data wphy /0.5,0.5,0,0,0.5,0.5/
#ifdef ALLOW_DIAGNOSTICS      
      _RL  mmhgpsink_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif /* ALLOW_DIAGNOSTICS */    
#endif /* HGPREFRACTORY     */
      
      _RL  poco, fpoco
      _RL  tiny
      data tiny /1d-32/

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_HG

C Initialize diagnostics
#ifdef ALLOW_DIAGNOSTICS
      DO k=1,Nr
        DO j=jMin,jMax
          DO i=iMin,iMax
            hgpsink_diag(i,j,k)=0d0
#ifdef HGPREFRACTORY 
            hgprsink_diag(i,j,k)=0d0
#endif /* HGPREFRACTORY */                   
#ifdef METHY_MODE 
            mmhgpsink_diag(i,j,k)=0d0
#endif /* METHY_MODE */         
          ENDDO
        ENDDO
      ENDDO
#endif /* ALLOW_DIAGNOSTICS */   
      
      DO j=jMin,jMax
       DO i=iMin,iMax

C       Zero flux
        DO k=1,Nr
         Fout_Hgp(k) = 0. _d 0
#ifdef HGPREFRACTORY
         DO nriver=1,nrmax 
          Fout_Hgpr(k,nriver) = 0. _d 0
         ENDDO
#endif /* HGPREFRACTORY */                
#ifdef METHY_MODE 
         Fout_MMHgp(k) = 0. _d 0
#endif /* METHY_MODE */        
        ENDDO

C       Flux calculation
        DO k=1,Nr
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN
           ! localize variables
           poco =  max(poc(i,j,k,bi,bj),0d0)  ! real poc
           fpoco = poco * 10.0 _d 0           ! real poc sinking
           DO np=1,npmax
               fpoco = fpoco + max(phytoP(i,j,k,bi,bj,np),0. _d 0) !poc & phyto sink
     &                * wphy(np)
               poco = poco + max(phytoP(i,j,k,bi,bj,np),0. _d 0) !mmolC/m3, poc+phyto
           ENDDO
           DO nz=1,nzmax
               !note: zooplankton not sinking 
               poco = poco + max(zooP(i,j,k,bi,bj,nz),0. _d 0) !mmolC/m3
           ENDDO
           IF ( poco .gt. tiny .and. fpoco .gt. tiny ) THEN   
             Fout_Hgp(k) = max(PTR(i,j,k,bi,bj,ihgp),0d0)
     &          * fpoco / poco  ! include phytoplankton sinking
     &          / 8.64 _d 4
           ELSE
             Fout_Hgp(k) = 0. _d 0
           ENDIF 
     
C Calculate diagnistics
#ifdef ALLOW_DIAGNOSTICS
           hgpsink_diag(i,j,k) = Fout_Hgp(k)  !mol/m^2/s   
#endif /* ALLOW_DIAGNOSTICS */

#ifdef HGPREFRACTORY
           DO nriver=1,nrmax
            Fout_Hgpr(k,nriver) = max(PTR(i,j,k,bi,bj,
     &            ihgpr+nriver-1),0d0)
#ifdef HGPR_IRREDUCIBLE 
     &          * fpoco / poco  ! include phytoplankton sinking
#elif defined(HGPR_SLOWRELEASING)
     &          * whgr(nriver)
#endif
     &          / 8.64 d 4

#ifdef WALSH_FILTER
#ifdef walsh_depth_0
C           do nothing
            IF ( 1 .eq. 0 ) THEN
#endif                
#ifdef walsh_depth_1
c           the particle sinking at river mouth is handled by the Walsh Typology filter, 4/18/2014
            IF ( hgprriver(i,j,bi,bj,nriver) .gt. tiny ) THEN
#endif                
#ifdef walsh_depth_2
c           the particle sinking over region shallower than 55m depth is handled by the Walsh Typology filter, 4/19/2014
            IF ( HFacC(i,j,5,bi,bj) .lt. tiny ) THEN
#endif                
#ifdef walsh_depth_3
c           the particle sinking over region shallower than 185m depth is handled by the Walsh Typology filter, 4/19/2014
            IF ( HFacC(i,j,9,bi,bj) .lt. tiny ) THEN
#endif
#ifdef walsh_depth_4
c           the particle sinking over region shallower than 55m depth, or
C           at the river mouth is handled by the Walsh Typology filter, 4/19/2014
            IF ( HFacC(i,j,5,bi,bj) .lt. tiny .or. 
     &           hgprriver(i,j,bi,bj,nriver) .gt. tiny) THEN
#endif                  
             Fout_Hgpr(k,nriver) = 0. _d 0
            ENDIF
#endif /* WALSH_FILTER      */
           ENDDO
#ifdef ALLOW_DIAGNOSTICS
           DO nriver=1,nrmax
            hgprsink_diag(i,j,k) = hgprsink_diag(i,j,k)
     &           +Fout_Hgpr(k,nriver)  !mol/m^2/s
           ENDDO 
#endif /* ALLOW_DIAGNOSTICS */
#endif /* HGPREFRACTORY     */
                                                
#ifdef METHY_MODE
           Fout_MMHgp(k) = max(PTR(i,j,k,bi,bj,immhgp),0d0)
     &          * 10.0 _d 0    ! poc w/o phytoplankton
     &          / 8.64 _d 4
           DO np=1,npmax
             Fout_MMHgphy(k,np) = max(PTR(i,j,k,bi,bj,iPhy+np-1),0d0)
     &          * wphy(np)     ! phytoplankton
     &          / 8.64 _d 4
           ENDDO
#ifdef ALLOW_DIAGNOSTICS
           mmhgpsink_diag(i,j,k) = Fout_MMHgp(k)   !mol/m^2/s   
           DO np=1,npmax
             mmhgpsink_diag(i,j,k) = mmhgpsink_diag(i,j,k)
     &          + Fout_MMHgphy(k,np)
           ENDDO
#endif /* ALLOW_DIAGNOSTICS */
#endif /* METHY_MODE        */
         ENDIF
        ENDDO

C----&---1---------2---------3---------4---------5---------6---------7--    
C       Tendency calculation
        IF (HFacC(i,j,1,bi,bj).gt.0. _d 0) THEN
c         Gchem(i,j,1,bi,bj,ihgp)=Gchem(i,j,1,bi,bj,ihgp) 
c     &                          - Fout_Hgp(1) * recip_drF(1)
         Ptr(i,j,1,bi,bj,ihgp)=Ptr(i,j,1,bi,bj,ihgp)
     &                - Fout_Hgp(1) * recip_drF(1) * deltaTclock 
#ifdef HGPREFRACTORY
         DO nriver=1,nrmax
c          Gchem(i,j,1,bi,bj,ihgpr+nriver-1)=
c     &     Gchem(i,j,1,bi,bj,ihgpr+nriver-1) 
c     &     - Fout_Hgpr(1,nriver) * recip_drF(1)
          Ptr(i,j,1,bi,bj,ihgpr+nriver-1)=
     &     Ptr(i,j,1,bi,bj,ihgpr+nriver-1)
     &     - Fout_Hgpr(1,nriver) * recip_drF(1) * deltaTclock
         ENDDO
#endif /* HGPREFRACTORY */        
#ifdef METHY_MODE         
c         Gchem(i,j,1,bi,bj,immhgp)=Gchem(i,j,1,bi,bj,immhgp)
c     &     - Fout_MMHgp(1) * recip_drF(1)
         Ptr(i,j,1,bi,bj,immhgp)=Ptr(i,j,1,bi,bj,immhgp)
     &     - Fout_MMHgp(1) * recip_drF(1) * deltaTclock
         DO np=1,npmax
          Ptr(i,j,1,bi,bj,iphy+np-1)=Ptr(i,j,1,bi,bj,iphy+np-1)
     &     - Fout_MMHgphy(1,np) * recip_drF(1) * deltaTclock
         ENDDO
#endif /* METHY_MODE */     
        ENDIF
  
        DO k=2,Nr
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN
c          Gchem(i,j,k,bi,bj,ihgp)=Gchem(i,j,k,bi,bj,ihgp) + 
c     &             (Fout_Hgp(k-1) - Fout_Hgp(k)) * recip_drF(k)
          Ptr(i,j,k,bi,bj,ihgp)=Ptr(i,j,k,bi,bj,ihgp) + 
     &             (Fout_Hgp(k-1) - Fout_Hgp(k)) * recip_drF(k)
     &             * deltaTclock
#ifdef HGPREFRACTORY
          DO nriver=1,nrmax     
c           Gchem(i,j,k,bi,bj,ihgpr+nriver-1)=
c     &        Gchem(i,j,k,bi,bj,ihgpr+nriver-1) + 
c     &        (Fout_Hgpr(k-1,nriver) - Fout_Hgpr(k,nriver)) 
c     &        * recip_drF(k)
           Ptr(i,j,k,bi,bj,ihgpr+nriver-1)=
     &        Ptr(i,j,k,bi,bj,ihgpr+nriver-1) +
     &        (Fout_Hgpr(k-1,nriver) - Fout_Hgpr(k,nriver))
     &        * recip_drF(k) * deltaTclock
          ENDDO
#endif /* HGPREFRACTORY */
#ifdef METHY_MODE     
c          Gchem(i,j,k,bi,bj,immhgp)=Gchem(i,j,k,bi,bj,immhgp) + 
c     &             (Fout_MMHgp(k-1) - Fout_MMHgp(k)) * recip_drF(k)
          Ptr(i,j,k,bi,bj,immhgp)=Ptr(i,j,k,bi,bj,immhgp) +
     &             (Fout_MMHgp(k-1) - Fout_MMHgp(k)) * recip_drF(k)
     &             * deltaTclock
          DO np=1,npmax
           Ptr(i,j,k,bi,bj,iphy+np-1)=Ptr(i,j,k,bi,bj,iphy+np-1)
     &             + (Fout_MMHgphy(k-1,np) - Fout_MMHgphy(k,np)) 
     &             * recip_drF(k) * deltaTclock
          ENDDO
#endif /* METHY_MODE */   
         ENDIF
        ENDDO

       ENDDO
      ENDDO

c write diagnostics
#ifdef ALLOW_DIAGNOSTICS
        WRITE(diagname,'(A8)') 'HGPSINK '
        CALL DIAGNOSTICS_FILL( hgpsink_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#ifdef HGPREFRACTORY
        WRITE(diagname,'(A8)') 'HGPRSINK'
        CALL DIAGNOSTICS_FILL( hgprsink_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#endif /* HGPREFRACTORY     */
#ifdef METHY_MODE
        WRITE(diagname,'(A8)') 'MHGPSINK'
        CALL DIAGNOSTICS_FILL( mmhgpsink_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#endif /* METHY_MODE        */
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_HG          */
#endif /* ALLOW_PTRACERS    */

       RETURN
       END
!EOC
