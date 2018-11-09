!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: hg_foodweb.F
!
! !DESCRIPTION: To apply foodweb dynamics to Hg species.
!  (Yanxu Zhang, 5/30/13)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HG_FOODWEB(   Ptr,
     &                         Gchem,
     &                         bi,bj,imin,imax,jmin,jmax,
     &                         myIter,myTime,myThid)
!
! !USES:
!
#include "GCHEM_OPTIONS.h"
#include "HG_OPTIONS.h"

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
      INTEGER myIter
      _RL myTime
      INTEGER myThid
      INTEGER bi, bj, imin, imax, jmin, jmax
!
! !INPUT/OUTPUT PARAMETERS: 
!
      _RL  Ptr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      _RL  Gchem(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
!
! !REMARKS:
!     /==========================================================\
!     | SUBROUTINE HG_FOODWEB                                    |
!     | o Calculate the changes to Hg species because            |
!     |   of foodweb dynamics                                    |
!     |   based on Asif Quesi 2013 paper                         |
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
#ifdef FOODW_MODE
C     == Local variables ==            
      INTEGER i,j,k
      INTEGER np,nz
      _RL phytopo(npmax), zoopo(nzmax)
      _RL phygrowo(npmax),zoograzeo(npmax,nzmax),zoomorto(nzmax)
      _RL grazeo,deatho, phymmhgo(npmax), zoommhgo(nzmax)
      _RL tiny
      PARAMETER(tiny=1d-32)
      
C     == Plankton variables
      _RL e_mehg_assim(nzmax)   ! assimilation efficiency
      _RL k_mehg_excret(nzmax)   ! excretion rate     
      DATA e_mehg_assim / 7 _d -1, 7 _d -1 / ! Lee 2015
      DATA k_mehg_excret / 5.787 _d -7, 5.787 _d -7  / ! 0.05 d-1, Tsui and Wang ES&T 2004/

#ifdef ALLOW_DIAGNOSTICS      
      CHARACTER*8   diagName
      _RL  mmhggz_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &                  Nr,npmax,nzmax)  ! MMHG graze
      _RL  mmhgmt_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nzmax)   ! release of MMHg from zooplankton after death
      _RL  mmhgfp_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)         ! release of MMHg via floopy feeding
#endif /* ALLOW_DIAGNOSTICS */

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_HG
      
C     == Biomagnification
      DO j=jMin,jMax
       DO i=iMin,iMax
        DO k=1,Nr
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN
           ! Localize variables
           DO np=1,npmax
               phytopo(np) = max(phytoP(i,j,k,bi,bj,np),0. _d 0)   !mmol C/m3
               phymmhgo(np) = max(Ptr(i,j,k,bi,bj,iPhy+np-1),0. _d 0) !mol/m3
           ENDDO
           DO nz=1,nzmax
               zoommhgo(nz) = max(Ptr(i,j,k,bi,bj,iZoo+nz-1),0. _d 0) !mol/m3
               zoopo(nz) = max(zooP(i,j,k,bi,bj,nz),0. _d 0)       !mmolC/m3
               zoomorto(nz) = max(zoomort(i,j,k,bi,bj,nz),0. _d 0) !mmolC/m3/s 
               DO np=1,npmax
                 zoograzeo(np,nz)=max(zoograze(i,j,k,bi,bj,np,nz), !mmolC/m3/s
     &                                0. _d 0)
               ENDDO 
           ENDDO

           ! Graze
           mmhgfp_diag(i,j,k) = 0. _d 0
           DO np=1,npmax
           DO nz=1,nzmax
            IF ( zoograzeo(np,nz) .gt. tiny .and. 
     &           phytopo(np) .gt. tiny ) THEN            
             grazeo = zoograzeo(np,nz) / phytopo(np)                 ! /s
             grazeo = ( 1. _d 0 - exp( - grazeo * deltaTclock ) ) 
     &                * phymmhgo(np) / deltaTclock                   ! mol/m3/s
c             Gchem(i,j,k,bi,bj,iPhy+np-1) = 
c     &                Gchem(i,j,k,bi,bj,iPhy+np-1) -
c     &                grazeo
c             Gchem(i,j,k,bi,bj,iZoo+nz-1) = 
c     &                Gchem(i,j,k,bi,bj,iZoo+nz-1) +
c     &                grazeo * e_mehg_assim(nz)
c             Gchem(i,j,k,bi,bj,immhg) = 
c     &                Gchem(i,j,k,bi,bj,immhg)
c     &                + grazeo * (1.0 _d 0 - e_mehg_assim(nz))
              Ptr(i,j,k,bi,bj,iPhy+np-1) = 
     &                Ptr(i,j,k,bi,bj,iPhy+np-1) -
     &                grazeo * deltaTclock 
              Ptr(i,j,k,bi,bj,iZoo+nz-1) = 
     &                Ptr(i,j,k,bi,bj,iZoo+nz-1) +
     &                grazeo * e_mehg_assim(nz) * deltaTclock
              Ptr(i,j,k,bi,bj,immhg) = 
     &                Ptr(i,j,k,bi,bj,immhg)
     &                + grazeo * (1.0 _d 0 - e_mehg_assim(nz)) 
     &                * deltaTclock
#ifdef ALLOW_DIAGNOSTICS
             mmhggz_diag(i,j,k,np,nz) = grazeo
             mmhgfp_diag(i,j,k) = mmhgfp_diag(i,j,k) + grazeo * 
     &                (1.0 _d 0 - e_mehg_assim(nz))  
#endif /* ALLOW_DIAGNOSTICS */
            ELSE
#ifdef ALLOW_DIAGNOSTICS
             mmhggz_diag(i,j,k,np,nz) = 0. _d 0
#endif /* ALLOW_DIAGNOSTICS */
            ENDIF
           ENDDO
           ENDDO

           ! release via death and excretion
           DO nz=1,nzmax
            IF ( zoomorto(nz) .gt. tiny .and.
     &           zoopo(nz) .gt. tiny ) THEN
             deatho = zoomorto(nz) / zoopo(nz)     ! death
             deatho = deatho + k_mehg_excret(nz)   ! excretion
             deatho = ( 1. _d 0 - exp( - deatho * deltaTclock ) ) 
     &                * zoommhgo(nz) / deltaTclock                    ! mol/m3/s
c             Gchem(i,j,k,bi,bj,iZoo+nz-1)=Gchem(i,j,k,bi,bj,iZoo+nz-1)
c     &                - deatho
c             Gchem(i,j,k,bi,bj,immhg)=Gchem(i,j,k,bi,bj,immhg)
c     &                + deatho
             Ptr(i,j,k,bi,bj,iZoo+nz-1)=Ptr(i,j,k,bi,bj,iZoo+nz-1)
     &                - deatho * deltaTclock
             Ptr(i,j,k,bi,bj,immhg)=Ptr(i,j,k,bi,bj,immhg)
     &                + deatho * deltaTclock
#ifdef ALLOW_DIAGNOSTICS
             mmhgmt_diag(i,j,k,nz) = deatho
#endif /* ALLOW_DIAGNOSTICS */
           ELSE
#ifdef ALLOW_DIAGNOSTICS
             mmhgmt_diag(i,j,k,nz) = 0. _d 0
#endif /* ALLOW_DIAGNOSTICS */
           ENDIF
          ENDDO
         ENDIF
        ENDDO
       ENDDO
      ENDDO

c write diagnostics
#ifdef ALLOW_DIAGNOSTICS
        DO np=1,npmax
         DO nz=1,nzmax
          WRITE(diagName,'(A4,I2.2,I2.2)') 'MHGZ',np,nz
          CALL DIAGNOSTICS_FILL( mmhggz_diag(1-Olx,1-Oly,1,np,nz),
     &                          diagname,
     &                          0,Nr,2,bi,bj,myThid )
         ENDDO
        ENDDO
        DO nz=1,nzmax
         WRITE(diagName,'(A6,I2.2)') 'MHMORT',nz
         CALL DIAGNOSTICS_FILL( mmhgmt_diag(1-Olx,1-Oly,1,nz), diagname,
     &                          0,Nr,2,bi,bj,myThid )
        ENDDO
        WRITE(diagName,'(A8)') 'MMHGFLOP'
        CALL DIAGNOSTICS_FILL( mmhgfp_diag(1-Olx,1-Oly,1), diagname,
     &                          0,Nr,2,bi,bj,myThid ) 
#endif /* ALLOW_DIAGNOSTICS */

#endif /* FOODW_MODE        */
#endif /* ALLOW_HG          */
#endif /* ALLOW_PTRACERS    */
       RETURN
       END
