!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: 
!
! !DESCRIPTION: Handles the partitioning between Hg species.
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "HG_OPTIONS.h"

      SUBROUTINE HG_PARTITION( Ptr,
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
#include "HG_SIZE.h"
#include "HG.h"
!
! !INPUT PARAMETERS: 
!
      INTEGER myIter
      _RL myTime
      INTEGER myThid
!
! !INPUT/OUTPUT PARAMETERS: 
!
      _RL  Ptr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
!
! !REMARKS:
!     /==========================================================\
!     | SUBROUTINE HG_PARTITION                                  |
!     | o Calculate the changes to Hg species because            |
!     |   of partitioning onto suspended particulate matter      |
!     |   based on Anne Sorerensen 2010 paper                    |
!     |   We assume an instantaneous equilibrium                 |
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
      _RL  hg2part_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#ifdef METHY_MODE      
      _RL  mmhgpart_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif /* METHY_MODE        */
      CHARACTER*8 diagname
#endif /* ALLOW_DIAGNOSTICS */
      _RL  TOT_HG, FRAC_DISSOLVED, poco, doco
      INTEGER i,j,k,bi, bj, imin, imax, jmin, jmax
      _RL  hg2_old,hgp_old,hg2_new,hgp_new
#ifdef HG2REFRACTORY
      _RL hg2r_old,hg2r_new
#endif /* HG2REFRACTORY     */          
#if defined(HG2REFRACTORY) & defined(HGPR_IRREDUCIBLE)
      _RL hgpr_old,hgpr_new
#endif /* HG2REFRACTPRY & HGPR_IRREDUCIBLE */
#ifdef METHY_MODE
      INTEGER np,nz
      _RL phytopt,zoopt
      _RL  mmhg_old,mmhgp_old,mmhg_new,mmhgp_new
#endif /* METHY_MODE        */     

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_HG

C     == Parameters
C     -- Partition coefficient
C     Kd value based on NP and North Atlantic measurements (5.5+/-0.5)
C     if based on poc, the Kd should be multiplied by 10: (6.5+/-0.5)
C     Soerensen et al. 2010; Zhang et al. 2014
      _RL Kd_hg2_poc  ! L/kg
C      data  Kd_hg2_poc / 3.5d6 /  ! Kd=10**5.5, for DARWIN field 05.27.2014
C                                    which has ballast effect and very
C                                    fast POC remineralization. This is
C                                    now obselete.
      data Kd_hg2_poc   / 1.4d6 /  ! Kd=10**5.15, for DARWIN field with
C                                    updated temperature dependent of
C                                    remineralization based on Marsay
C                                    PNAS 2015 paper. No ballast effect
C                                    is introduced.
C      data  Kd_hg2_poc / 3.2d5 /  ! Kd=10**4.5 for Arctic, (Fisher et
C      al. 2012)
#ifdef METHY_MODE 
      _RL Kd_mmhg_poc  ! L/kg
      data Kd_mmhg_poc /6.3 _d 4/  !Note: this is Koc, i.e. = 10 * Kd, using logKd=3.8
#endif /* METHY_MODE        */
C     == Start calculation
C Initialize diagnostics
#ifdef ALLOW_DIAGNOSTICS
      DO k=1,Nr
        DO j=jMin,jMax
          DO i=iMin,iMax
            hg2part_diag(i,j,k)=0d0
#ifdef METHY_MODE            
            mmhgpart_diag(i,j,k)=0d0
#endif /* METHY_MODE        */            
          ENDDO
        ENDDO
      ENDDO
#endif /* ALLOW_DIAGNOSTICS */ 

      DO j=jMin,jMax
       DO i=iMin,iMax
        DO k=1,Nr

         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN
           !===========================================================
           ! %%%%% Instantaneous equilibrium  of Hg2, HgP%%%%%
           !===========================================================

           ! POC and DOC concentrations
           poco = max(poc(i,j,k,bi,bj),0d0)
           doco = max(doc(i,j,k,bi,bj),0d0)

           ! Phytoplankton (mmolC/m3)
           phytopt = 0. _d 0
           DO np=1,npmax
               phytopt = phytopt + max(phytoP(i,j,k,bi,bj,np),0. _d 0)
           ENDDO

           ! Zooplankton (mmolC/m3)
           zoopt = 0. _d 0
           DO nz=1,nzmax
               zoopt = zoopt + max(zooP(i,j,k,bi,bj,nz),0. _d 0)
           ENDDO
 
           !-----------------------------------------------------------
           ! Hg2, HgP concentrations before partitioning 
           !-----------------------------------------------------------
           hg2_old = max(Ptr(i,j,k,bi,bj,ihg2),0d0)
           hgp_old = max(Ptr(i,j,k,bi,bj,ihgp),0d0)
           TOT_HG  = hg2_old + hgp_old

           IF (TOT_HG .gt. 0. _d 0) THEN

             !---------------------------------------------------------
             ! Calculate fraction in the dissolved phase
             !---------------------------------------------------------
             FRAC_DISSOLVED = 1. _d 0 /
     &             (1. _d 0 + Kd_hg2_poc * ( poco + phytopt + zoopt)
     &             * 1.2 _d - 8 )  ! mmolC/m3 -> kgC/L

             !---------------------------------------------------------
             ! New Hg2, Hgp concentrations (units?) in the dissolved
             ! and particle phases after consdering partitioning
             !---------------------------------------------------------
             hg2_new = TOT_HG * FRAC_DISSOLVED
             hgp_new = TOT_HG * (1. _d 0 - FRAC_DISSOLVED)

             !---------------------------------------------------------
             ! Update tracer array
             !---------------------------------------------------------
             Ptr(i,j,k,bi,bj,ihg2) = hg2_new
             Ptr(i,j,k,bi,bj,ihgp) = hgp_new

!             Gchem(i,j,k,bi,bj,ihg2) = Gchem(i,j,k,bi,bj,ihg2)
!     &           + ( hg2_new - hg2_old ) / deltaTclock
!             Gchem(i,j,k,bi,bj,ihgp) = Gchem(i,j,k,bi,bj,ihgp)
!     &           + ( hgp_new - hgp_old ) / deltaTclock 

#ifdef ALLOW_DIAGNOSTICS
             hg2part_diag(i,j,k) = (hg2_old - hg2_new
     &                 ) / deltaTclock
#endif /* ALLOW_DIAGNOSTICS */            

           ENDIF                  ! TOT_HG > 0

#if defined(HG2REFRACTORY) & defined(HGPREFRACTORY) &  defined(HGPR_IRREDUCIBLE)
           !===========================================================
           ! %%%%% Instantaneous equilibrium  of Hg2R, HgPR %%%%%
           !===========================================================

           !-----------------------------------------------------------
           ! Hg2R, HgPR concentrations before partitioning 
           !-----------------------------------------------------------
           hg2r_old = max(Ptr(i,j,k,bi,bj,ihg2r),0d0)
           hgpr_old = max(Ptr(i,j,k,bi,bj,ihgpr),0d0)
           TOT_HG   = hg2r_old + hgpr_old

           IF (TOT_HG .gt. 0. _d 0) THEN

             FRAC_DISSOLVED = 1. _d 0 /
     &             (1. _d 0 + Kd_hg2_poc * ( poco + phytopt + zoopt)
     &             * 1.2 _d - 8 )  ! mmolC/m3 -> kgC/L

             !---------------------------------------------------------
             ! New Hg2R, HgpR concentrations (units?) in the dissolved
             ! and particle phases after consdering partitioning
             !---------------------------------------------------------
             hg2r_new = TOT_HG * FRAC_DISSOLVED
             hgpr_new = TOT_HG * (1. _d 0 - FRAC_DISSOLVED)

             !---------------------------------------------------------
             ! Update tracer array
             !---------------------------------------------------------
             Ptr(i,j,k,bi,bj,ihg2r) = hg2r_new
             Ptr(i,j,k,bi,bj,ihgpr) = hgpr_new

!             Gchem(i,j,k,bi,bj,ihg2r) = Gchem(i,j,k,bi,bj,ihg2r)
!     &           + ( hg2r_new - hg2r_old ) / deltaTclock
!             Gchem(i,j,k,bi,bj,ihgpr) = Gchem(i,j,k,bi,bj,ihgpr)
!     &           + ( hgpr_new - hgpr_old ) / deltaTclock

#ifdef ALLOW_DIAGNOSTICS
             hg2part_diag(i,j,k) = hg2part_diag(i,j,k)
     &                 + (hg2r_old - hg2r_new
     &                 ) / deltaTclock
#endif /* ALLOW_DIAGNOSTICS */             
 
           ENDIF                  ! TOT_HG > 0d0

#endif /* HG2REFRACTORY & HGPREFRACTORY & HGPRIRREDUCIBLE */

#ifdef METHY_MODE
           !===========================================================
           ! %%%%% Instantaneous equilibrium of MMHg, MMHgP %%%%%
           !===========================================================

           !-----------------------------------------------------------
           ! MMHg, MMHgp concentrations before partitioning 
           !-----------------------------------------------------------
           mmhg_old  = max(Ptr(i,j,k,bi,bj,immhg),0d0)
           mmhgp_old = max(Ptr(i,j,k,bi,bj,immhgp),0d0)
           TOT_HG    = mmhg_old + mmhgp_old

           IF (TOT_HG .gt. 0. _d 0) THEN

             !---------------------------------------------------------
             ! Calculate fraction in the dissolved phase
             !---------------------------------------------------------
             FRAC_DISSOLVED = 1. _d 0 /
     &             (1. _d 0 + Kd_mmhg_poc * poco   
             ! Note: phyto partition in bioaccum.F
     &             * 1.2 _d - 8 )  ! mmolC/m3 -> kgC/L

             !---------------------------------------------------------
             ! New MMHg, MMHgP concentrations (units?) in the dissolved
             ! and particle phases after consdering partitioning
             !---------------------------------------------------------
             mmhg_new  = TOT_HG * FRAC_DISSOLVED
             mmhgp_new = TOT_HG - mmhg_new

             !---------------------------------------------------------
             ! Update tracer array
             !---------------------------------------------------------
             Ptr(i,j,k,bi,bj,immhg)  = mmhg_new
             Ptr(i,j,k,bi,bj,immhgp) = mmhgp_new

!             Gchem(i,j,k,bi,bj,immhg) = Gchem(i,j,k,bi,bj,immhg)
!     &           + ( mmhg_new - mmhg_old ) / deltaTclock
!             Gchem(i,j,k,bi,bj,immhgp) = Gchem(i,j,k,bi,bj,immhgp)
!     &           + ( mmhgp_new - mmhgp_old ) / deltaTclock

#ifdef ALLOW_DIAGNOSTICS
             mmhgpart_diag(i,j,k) = (mmhg_old - mmhg_new) / deltaTclock
#endif /* ALLOW_DIAGNOSTICS */            

           ENDIF                  ! TOT_HG > 0

#endif /* METHY_MODE        */

         ENDIF                    ! HFacC(i,j,k,bi,bj) > 0
        ENDDO
       ENDDO
      ENDDO

c write diagnostics
#ifdef ALLOW_DIAGNOSTICS
        WRITE(diagname,'(A8)') 'HG2PART '
        CALL DIAGNOSTICS_FILL( hg2part_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#ifdef METHY_MODE     
        WRITE(diagname,'(A8)') 'MMHGPART'
        CALL DIAGNOSTICS_FILL( mmhgpart_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#endif /* METHY_MODE        */
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_HG          */
#endif /* ALLOW_PTRACERS    */

       RETURN
       END
!EOC
