!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: hg0_surfforcing.F
!
! !DESCRIPTION: Surface forcing of the Hg0 tracer. 
!\\
!\\
! !INTERFACE:
!
#include "CPP_OPTIONS.h"
#include "GCHEM_OPTIONS.h"
#include "HG_OPTIONS.h"
#ifdef USE_EXFIWR
#include "SEAICE_OPTIONS.h"
#include "EXF_OPTIONS.h"
#endif

      SUBROUTINE Hg0_SURFFORCING( PTR, surfhg0,
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
#include "FFIELDS.h"

#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "GCHEM.h"
#include "HG_SIZE.h"
#include "HG.h"

#ifdef USE_EXFIWR
#include "SEAICE_SIZE.h"
#include "SEAICE.h"
#include "EXF_FIELDS.h"
#endif
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
      _RL  surfhg0(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
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

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_HG

C     == Local variables ==
C I, J, K - Loop counters
      INTEGER I,J
C Solubility relation coefficients
      _RL SchmidtNohg0, ScCO2
      _RL Kwexch
      _RL FluxHG0(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL tk, tc, vi, D, H
      _RL doco, iceo, pisvelo
      _RL total_river_1
      LOGICAL isrivermouth
      INTEGER i1,j1,month
#ifdef ALLOW_DIAGNOSTICS      
      _RL hg0conc_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL hg0evas_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL hg0conco
      CHARACTER*8 diagname
#endif /* ALLOW_DIAGNOSTICS */

C Initialize diagnostics
#ifdef ALLOW_DIAGNOSTICS
      DO j=jMin,jMax
        DO i=iMin,iMax
            hg0conc_diag(i,j)=0d0
            hg0evas_diag(i,j)=0d0
        ENDDO
      ENDDO
#endif /* ALLOW_DIAGNOSTICS */

C Get the current month
#ifdef USE_EXFIWR
        month=mod(int(myTime/2629800.0), 12) + 1
C        WRITE(standardmessageunit,'(A,1x,I2)')
C     &               'current month: ',month
#endif /* USE_EXFIWR */
                  
C calculate SCHMIDT NO.
        DO j=jMin,jMax
          DO i=iMin,iMax
C Locoalize variables              
            hg0conco = max(hg0conc(i,j,bi,bj),0d0)
#ifdef USE_EXFIWR
            iceo = AREA(i,j,bi,bj)
C           Soerensen et al. 2010: 4x5
C           pisvelo = 6.94 _d -7 * wspeed(i,j,bi,bj) ** 2
C           Nightingale et al. 2010: instantaneous wind
            pisvelo = (0.222 * wspeed(i,j,bi,bj) ** 2 + 
     &                 0.333 * wspeed(i,j,bi,bj)) *
     &                 2.7778e-06  ! cm/h->m/s
#else
            iceo = fice(i,j,bi,bj)
            pisvelo = pisvel(i,j,bi,bj)
#endif /* USE_EXFIWR */
            IF (HFacC(i,j,1,bi,bj).GT.0.5) THEN

C Is this a river mouth grid
#ifdef USE_EXFIWR
              total_river_1 = 0d0
              DO j1=j-4,j+4
              DO i1=i-4,i+4
               IF (i1 .GE. iMin .AND. i1 .LE. iMax .AND.
     &             j1 .GE. jMin .AND. j1 .LE. jMax .AND.
     &             HFacC(i1,j1,1,bi,bj) .GT. 0.5 _d 0 ) THEN
                total_river_1 = total_river_1
     &           + hg2river(i1,j1,bi,bj) + hgpriver(i1,j1,bi,bj)
#ifdef HG2REFRACTORY
     &           + hg2rriver(i1,j1,bi,bj)
#endif /* HG2REFRACTORY */
#ifdef HGPREFRACTORY
     &           + hgprriver(i1,j1,bi,bj,1)
#endif /*HGPREFRACTORY */ 
               ENDIF
              ENDDO
              ENDDO

              IF ( total_river_1 .gt. 1d-32 ) THEN
                 isrivermouth = .TRUE.
              ELSE
                 isrivermouth = .FALSE.
              ENDIF
#endif /* USE_EXFIWR */

C Locoalize variables
              doco     = max(doc(i,j,1,bi,bj),0d0)
                              
C calculate SCHMIDT NO. for Hg0
C yxzhang, theta in Celcuis degree
              tc = max(theta(i,j,1,bi,bj), 0 _d 0)
              !Note: this is a drawback of current settings. We don't
              !consider the temperature of runoff, and assume the runoff
              !temperature equals to that of SST. This is not the case.
              !According to the study by
              !Nghiem and colleagues (2014), the average surface temperature of
              !the open water rose by 6.5°C (11.7°Fahrenheit) after the
              !pulse from the Mackenzie. Temperatures were 10°C at 287
              !kilometers (178 miles) from the coast and 8°C as far as
              !350 kilometers (217 miles) away.
              !Improve this when we get better information for the
              !runoff temperature.
#ifdef USE_EXFIWR
              IF ( isrivermouth .AND. 
     &             month .GE. 6 .AND. month .LE. 7) THEN
C                disable this because setting runoff temperature in
C                exf_mapfields.F
                 tc = tc + 10.0
              ENDIF
#endif /* USE_EXFIWR */
              tk = tc + 273.15 _d 0
              vi = ( 10**( ( 1301.0d0 / ( 998.333d0 + 8.1855d0
     &              * ( tc - 20.0d0 )+ 0.00585d0 * (tc - 20.0d0 )**2 ) )
     &              - 3.30233d0 ) ) * 100.0d0
              D = 9.37 _d -8 * tk / vi
              SchmidtNohg0 =
     &            0.017 _d 0 * exp( - 0.025 _d 0 * tc ) / D
              ScCO2  = 644.7d0 + tc * ( -6.16d0 + tc * ( 0.11d0 ) )
c calculate Henry's law constant for Hg0, L air/L water
              H = exp( ( -2404.3d0 / tk ) + 6.92d0 )
C Determine surface flux (Fhg0), m/s

              ! Seasonal sea ice region, air-sea exchange is enhanced
              ! because of turbulance and mixing (Loose et al. Ocean
              ! Science 2014)
              IF ( (iceo .gt. 0.01 .and. iceo .lt. 0.99)
#ifdef USE_EXFIWR
     &             .OR. isrivermouth
#endif /* USE_EXFIWR */
     &             ) THEN
                  pisvelo = pisvelo * 2.0
              ENDIF
              ! increase another factor of 2 when fraction of open water
              ! is between 0.2 and 0.6, i.e. non linear effect in Fig.
              ! 10 of Loose et al. 2014.
              IF ( iceo .gt. 0.4 .and. iceo .lt. 0.8 ) THEN
                  pisvelo = pisvelo * 2.0
              ENDIF
              ! increase another factor of 2 when the salinity is
              ! largely decreased during June and July in river mouth.
              ! This is an assumption so far.
              IF ( month .gt. 5 .and. month .lt. 8 
     &              .and. isrivermouth ) THEN
                  pisvelo = pisvelo * 2.0
              ENDIF

              ! ice fraction is zero at river mouth (and nearby grids) if runoff exists               
#ifdef USE_EXFIWR
              IF ( isrivermouth .AND. 
     &             month .GT. 2 .AND. month .LT. 11 ) THEN
                  iceo = 0.0 _d 0 
              ENDIF
#endif /* USE_EXFIWR */
              Kwexch = (1. _d 0 - iceo)
     &                    * pisvelo
     &                    / sqrt( SchmidtNoHG0 / ScCO2 )

              FluxHG0(i,j) =
     &                    Kwexch * (hg0conco * 4.985 _d -12   !ng/m3->mol/m3
     &                    / H - max(PTR(i,j,1,bi,bj,ihg0),0d0))
            ELSE
              FluxHG0(i,j) = 0. _d 0
            ENDIF

C Calculate diagnistics
#ifdef ALLOW_DIAGNOSTICS
            hg0conc_diag(i,j) = HG0conco  !ng/m^3
            hg0evas_diag(i,j) = -FluxHG0(i,j)  !mol/m^2/s
#endif /* ALLOW_DIAGNOSTICS */        
          ENDDO
        ENDDO

C update surface tendencies      
        DO j=jMin,jMax
          DO i=iMin,iMax
           SURFHG0(i,j)=
     &         HFacC(i,j,1,bi,bj)*FluxHG0(i,j)*recip_drF(1)
          ENDDO
        ENDDO

C        write(*,*)'Evasion flux from Mackenzie River:',FluxHG0(152,116),
C     &         -FluxHG0(152,116) / 
C     &        (HG2RIVER(152,116,bi,bj)+HGPRIVER(152,116,bi,bj) * 
C     &        recip_rA(152,116,bi,bj)*4.985 )
C        write(*,*)'Shortwave radiation here:',swdown(152,116,bi,bj)
C        write(*,*)'Ice fraction here:',AREA(152,116,bi,bj)
       
c write diagnostics
#ifdef ALLOW_DIAGNOSTICS
        WRITE(diagname,'(A8)') 'HG0CONC '
        CALL DIAGNOSTICS_FILL( hg0conc_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )
        WRITE(diagname,'(A8)') 'HG0EVAS '
        CALL DIAGNOSTICS_FILL( hg0evas_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )     
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_HG          */
#endif /* ALLOW_PTRACERS    */
        RETURN
        END
!EOC
