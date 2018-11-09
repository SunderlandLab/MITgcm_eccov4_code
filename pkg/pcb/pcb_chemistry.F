!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pcb_chemistry.F
!
! !DESCRIPTION: Calculate the changes to dissolved-phased PCBs due to  
!   degradation.                                           
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "PCB_OPTIONS.h"

#ifdef USE_EXFIWR
#include "EXF_OPTIONS.h"
#include "SEAICE_OPTIONS.h"
#endif /* USE_EXFIWR */

      SUBROUTINE PCB_CHEMISTRY( Ptr,
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
#include "PCB_SIZE.h"
#include "PCB.h"

#ifdef USE_EXFIWR
#include "SEAICE_SIZE.h"
#include "SEAICE.h"
#include "EXF_FIELDS.h"
#include "FFIELDS.h"
#endif /* USE_EXFIWR */
!
! !INPUT PARAMETERS: 
!
      INTEGER myThid
      INTEGER myIter
      _RL myTime
      INTEGER i,j,k,l,m,bi, bj, imin, imax, jmin, jmax
!
! !INPUT/OUTPUT PARAMETERS: 
!
      _RL Ptr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      _RL Gchem(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)    
!
! !REMARKS:
!
! !REVISION HISTORY:
!  14 Jan 2015 - H. Amos     - Adopted copy of hg_chemistry.F
!                              for PCB degradation.
!  19 May 2015 - H. Amos     - Update USE_EXFIWR calls based on Hg
!                              so you can use online ice, radiation, 
!                              etc. when running with ECCOv4. Add
!                              call PCB_SUNCOS function.
!  01 Oct 2015 - H. Amos     - Add temperature dependence to degradation.
!  29 Feb 2016 - R. Yantosca - Added ProTex headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#ifdef ALLOW_DIAGNOSTICS
      _RL pcbaDegrad_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr) ! PCBa lost via degradation      
      CHARACTER*8 diagname
#endif /* ALLOW_DIAGNOSTICS */
      _RL,SAVE :: meansuncos(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      LOGICAL,SAVE :: FIRST
      DATA FIRST /.TRUE./
      INTEGER nl
      
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PCB
      !--------------------------------------------------------
      ! Local variables
      !--------------------------------------------------------
      _RL radz    , radz0    , ecz1, ecz2
      _RL suncos
      _RL pcba_old, pcba_new
      _RL Kdeg_pcb
      _RL T_CORREC ! HM Amos, added 01 Oct 2015

      !--------------------------------------------------------
      ! Parameters
      !--------------------------------------------------------

      ! Degradation
      _RL R                ! universal gas consntant [J/mol/K]
      _RL T0               ! standard temperature    [K]
      _RL C2K              ! conversion factor, degrees C to K
      _RL EACT             ! activation energy       [J/mol]
      data R    /8.3145d0/
      data T0   /298d0/
      data C2K  /273.15d0/
      data EACT /30000d0/  ! value selected by Friedman & Selin 
                           ! [2015] for consistency with Lamon
                           ! et al. [2009] and Gouin et al. [2013]
      

      ! Radiation extinction
      _RL EC_w    ! 1/m
      _RL EC_doc  ! 1/m
      _RL EC_chla ! 1/m
      data EC_w /0.0145 _d 0/
      data EC_doc /0.026 _d 0/
      data EC_chla /3.1 _d 1/

      !--------------------------------------------------------
      ! Local variables      
      !--------------------------------------------------------
      _RL iceo,nppo,chlo,fremino,doco,radswo,radswo0,meank     
      _RL myTimeo

      !--------------------------------------------------------
      ! Initialize diagnostics
      !--------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS
      DO k=1,Nr
        DO j=jMin,jMax
          DO i=iMin,iMax
            pcbaDegrad_diag(i,j,k)=0d0
          ENDDO
        ENDDO
      ENDDO
#endif /* ALLOW_DIAGNOSTICS */

      !--------------------------------------------------------
      ! Calculate monthly mean suncos
      !--------------------------------------------------------
#ifndef USE_EXFIWR
      IF ( mod( myTime, 2629800.0 ) .lt. 0.1 .or. FIRST) THEN
       FIRST = .FALSE.   
       nl = int(2629800 / deltaTclock)
       DO j=jMin,jMax
        DO i=iMin,iMax
         IF (HFacC(i,j,1,bi,bj).gt.0. _d 0) THEN
           meansuncos(i,j)=0. _d 0  
           DO l=1,nl
            myTimeo=myTime+dble(l-1)*deltaTclock  
            CALL PCB_SUNCOS(i,j,bi,bj,myTimeo,suncos)  
            meansuncos(i,j)=meansuncos(i,j)+suncos
           ENDDO
           meansuncos(i,j)=meansuncos(i,j)/dble(nl)   
         ENDIF                                   
        ENDDO
       ENDDO
      ENDIF
#endif /* USE_EXFIWR */
      
      DO j=jMin,jMax
       DO i=iMin,iMax
       IF (HFacC(i,j,1,bi,bj).gt.0. _d 0) THEN

#ifdef USE_EXFIWR
        iceo = AREA(i,j,bi,bj)

        ! radswo0 = swdown(i,j,bi,bj)*(1d0-iceo)
        ! Qsw   :: Net upward shortwave radiation in W/m^2
        !       Qsw = - ( downward - ice and snow absorption - reflected)
        !       > 0 for decrease in theta (ocean cooling)
        !       Typical range: -350 < Qsw < 0  
        ! If seaice package is enabled, ice fraction is already taken into
        ! account.     
        radswo0 = -Qsw(i,j,bi,bj)
        radswo = radswo0
#else
        iceo = fice(i,j,bi,bj)
        radswo0 = max(radsw(i,j,bi,bj),0d0)*(1d0-iceo)  ! wm-2

        ! Impose a diurnal cycle
        CALL PCB_SUNCOS(i,j,bi,bj,myTime,suncos)
        radswo = radswo0 * suncos / (meansuncos(i,j) + 1. _d -32)
#endif /* USE_EXFIWR */ 

 
        ! Radiation attenuation calculation
        ecz1 = 0d0
        ecz2 = 0d0

        !--------------------------------------------------------
        ! Localize variables
        !--------------------------------------------------------                          
        DO k=1,Nr
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN               
           
           !--------------------------------------------------------
           ! Local concentrations
           !--------------------------------------------------------
           pcba_old = max( PTR(i,j,k,bi,bj,ipcba), 0d0)

!Updated to be consistent with Carey Friedman's GEOS-Chem simulation
!and Matt McCleod's BETR-GLOBAL models.
!01 October 2015
!
!           !--------------------------------------------------------
!           ! Temperature dependence of degradation
!           ! Kdeg_pcb has units of 1/sec
!           !--------------------------------------------------------
!           Kdeg_pcb = a_deg * EXP( b_deg * theta(i,j,k,bi,bj) )

           !--------------------------------------------------------
           ! Temperature dependence of degradation
           !   Kdeg = a* ( ln(2)/tHalf )
           !   where a is a temperature correction
           !--------------------------------------------------------
           ! Temperature dependence
           T_CORREC = EXP((-EACT / R) 
     &                   * (1d0/(theta(i,j,k,bi,bj)+C2K) - 1/T0))

           ! Calculate degradation rate [1/sec]
           Kdeg_pcb = 0.693d0 / tHalf      ! [1/hr]
           Kdeg_pcb = Kdeg_pcb / 3600d0    ! convert [1/hr] --> [1/sec]
           Kdeg_pcb = T_CORREC * Kdeg_pcb  ! apply temperature dependence

           !--------------------------------------------------------
           ! Calculate new concentration of PCBa after 
           ! degradation (photolytic + biodegradation)
           !--------------------------------------------------------
           pcba_new = pcba_old - Kdeg_pcb * pcba_old / deltaTclock

           !--------------------------------------------------------
           ! Calculate diagnostics                                             
           !--------------------------------------------------------     
#ifdef ALLOW_DIAGNOSTICS
           ! bulk degradation (mol m^-3 sec^-1)
           pcbaDegrad_diag(i,j,k) = Kdeg_pcb * pcba_old / deltaTclock
#endif /* ALLOW_DIAGNOSTICS */

           !--------------------------------------------------------
           ! Accumulate tendency
           !--------------------------------------------------------
           Gchem(i,j,k,bi,bj,ipcba)=Gchem(i,j,k,bi,bj,ipcba)+
     &                             (pcba_new - pcba_old)
     &                             /deltaTclock            

         ENDIF
        ENDDO
        
       endif
       ENDDO
      ENDDO

      !-------------------------------------------------------------    
      !  Write diagnostics
      !-------------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS
      WRITE(diagname,'(A8)') 'PCBaDEGR'
      CALL DIAGNOSTICS_FILL( pcbaDegrad_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_PCB */
#endif /* ALLOW_PTRACERS */

      RETURN
      END
!EOC
