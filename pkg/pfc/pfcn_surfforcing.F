!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pfcn_surfforcing.F
!
! !DESCRIPTION: Surface forcing of the PFCn tracer.
!\\
!\\
! !INTERFACE:
!
#include "CPP_OPTIONS.h"
#include "GCHEM_OPTIONS.h"
#include "PFC_OPTIONS.h"
#ifdef USE_EXFIWR
#include "SEAICE_OPTIONS.h"
#include "EXF_OPTIONS.h"
#endif

      SUBROUTINE PFCn_SURFFORCING( Ptr, SurfPFCn,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )
!
! !USES:
!
      IMPLICIT NONE

C     == GLobal variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"

#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "GCHEM.h"
#include "PFC_SIZE.h"
#include "PFC.h"

#ifdef USE_EXFIWR
#include "SEAICE_SIZE.h"
#include "SEAICE.h"
#include "EXF_FIELDS.h"
#endif
!
! !INPUT PARAMETERS: 
!
      _RL myTime
      INTEGER iMin, iMax, jMin, jMax, bi, bj
      INTEGER myIter, myThid
!
! !INPUT/OUTPUT PARAMETERS: 
!
      _RL  PTR (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr, 
     &             nSx, nSy, PTRACERS_num)
      _RL  SurfPFCn (1-OLx:sNx+OLx,1-OLy:sNy+OLy)

! !REMARKS:
!
! !REVISION HISTORY:
!  29 Feb 2016 - R. Yantosca - Added ProTex headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PFC

C     pTr_PFC01  :: ocean PFC01 concentration
C     AtmosPFC01 :: atmospheric PFC01 field
C     fluxPFC01  :: air-sea PFC01 fluxes
C     bi, bj     :: current tile indices
C     iMin,iMax  :: computation domain, 1rst index bounds
C     jMin,jMax  :: computation domain, 2nd  index bounds
C     myTime     :: current time in simulation
C     myIter     :: current iteration number
C     myThid     :: my Thread Id number
      _RL  FluxPFCn(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  PFCnRiver(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  pfcnrivero, pfcnevaso
#ifdef ALLOW_SEASPRAY
      _RL pfcnseasprayo
#endif  /* ALLOW_SEASPRAY     */    
CEOP

C     !LOCAL VARIABLES:
C i, j, k - Loop counters
      INTEGER i,j,k
C Solubility relation coefficients
      _RL etaA, rhoA
      _RL etaW
      _RL SchmidtAPFCn
      _RL SchmidtWPFCn
      _RL DiffusivityAir(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL DiffusivitySurfWater(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL ka(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL kw(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL Kaexch(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
c      _RL Kwexch(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      
      _RL KawPFCn(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL KawPFCeff(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

c     _RL pfc01sat(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

c      _RL Csat(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL ttemp, ttemp2
      _RL temp1, temp2, temp3, temp4, temp5, temp6
      _RL doco, iceo
      _RL tc
      _RL total_river_1
      LOGICAL isrivermouth
      INTEGER i1, j1, month
#ifdef ALLOW_DIAGNOSTICS      
      _RL pfcnconc_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL pfcnevas_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#ifdef ALLOW_SEASPRAY
      _RL pfcnseaspray_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#endif /* ALLOW_SEASPRAY     */    
      _RL pfcnriver_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL pfcnconco
      CHARACTER*8 diagname
#endif  /* ALLOW_DIAGNOSTICS */    

                                    
C Initialize diagnostics
#ifdef ALLOW_DIAGNOSTICS
      DO j=jMin,jMax
        DO i=iMin,iMax
            pfcnconc_diag(i,j)=0.d0
            pfcnevas_diag(i,j)=0.d0
#ifdef ALLOW_SEASPRAY
            pfcnseaspray_diag(i,j)=0.d0
#endif
            pfcnriver_diag(i,j)=0.d0
        ENDDO
      ENDDO
#endif      

C Get the current month
#ifdef USE_EXFIWR
        month=mod(int(myTime/2629800.0), 12) + 1
C        WRITE(standardmessageunit,'(A,1x,I2)')
C     &               'current month: ',month
#endif /* USE_EXFIWR */

      
c      write(*,*) 'pfcn air-sea exchange'
C calculate SCHMIDT NO. 
        DO j=jMin,jMax
         DO i=iMin,iMax

#ifdef USE_EXFIWR
            iceo = AREA(i,j,bi,bj)
C           Nightingale et al. 2010: instantaneous wind
C            pisvelo = (0.222 * wspeed(i,j,bi,bj) ** 2 +
C     &                 0.333 * wspeed(i,j,bi,bj)) *
C     &                 2.7778e-06  ! cm/h->m/s
#else
            iceo = fice(i,j,bi,bj)
C            pisvelo = pisvel(i,j,bi,bj)
#endif /* USE_EXFIWR */

C           IF ( maskC(i,j,k,bi,bj).NE.zeroRS ) THEN !CFC code
C     yanxu hg code  IF (HFacC(i,j,1,bi,bj).GT.0.5) THEN 
            IF (HFacC(i,j,1,bi,bj).GT.0.1) THEN 
C Locoalize variables
              pfcnconco = max(pfcnconc(i,j,bi,bj),0.d0)
C            river : kg/s    -> mol/m2/s
             pfcnrivero = PFCnRiver(i,j,bi,bj)*recip_rA(i,j,bi,bj)
     &           *1000.d0/(Mb-1.d0)
              doco     = max(doc(i,j,1,bi,bj),0.d0)
     
C Is this a river mouth grid
#ifdef USE_EXFIWR
              total_river_1 = 0d0
              DO j1=j-4,j+4
              DO i1=i-4,i+4
               IF (i1 .GE. iMin .AND. i1 .LE. iMax .AND.
     &             j1 .GE. jMin .AND. j1 .LE. jMax .AND.
     &             HFacC(i1,j1,1,bi,bj) .GT. 0.1 _d 0 ) THEN
                total_river_1 = total_river_1
     &           + PFCiRiver(i1,j1,bi,bj) + PFCpRiver(i1,j1,bi,bj)
     &          +PFCnRiver(i,j,bi,bj)
               ENDIF
              ENDDO
              ENDDO

              IF ( total_river_1 .gt. 1d-32 ) THEN
                 isrivermouth = .TRUE.
              ELSE
                 isrivermouth = .FALSE.
              ENDIF
#endif /* USE_EXFIWR */





             
C calculate SCHMIDT NO. for PFCn used by XZ following Tsilingiris 2008

            tc = max(theta(i,j,1,bi,bj), 0 _d 0)
#ifdef USE_EXFIWR
              IF ( isrivermouth .AND.
     &             month .GE. 6 .AND. month .LE. 7) THEN
C                disable this because setting runoff temperature in
C                exf_mapfields.F
                 tc = tc + 10.0
              ENDIF
#endif /* USE_EXFIWR */


            etaA = sv_0
     &          + sv_1 * tc
     &          + sv_2 * tc*tc
     &          + sv_3 * tc**3
     &          + sv_4 * tc**4
C       write(*,*) 'xmz check sv_1, etaA', sv_1, etaA
            rhoA = sd_0
     &          + sd_1 * tc
     &          + sd_2 * tc*tc
     &          + sd_3 * tc**3
     
            temp1 = (Ma + Mb) / (Ma * Mb)
            temp2 = Va**(1.0d0/3.0d0) +Vb**(1.0d0/3.0d0)
            DiffusivityAir (i,j) = 1.0d-3*(tc+273.15d0)
     &       **1.75d0*temp1**0.5d0/(Patm*temp2**2.0d0) !cm2/s
         
            SchmidtAPFCn = etaA/rhoA/DiffusivityAir(i,j) 
            
c            DiffusivitySurfWater (i,j) !cm2/s
c     &   = 1.25d-8*(theta(i,j,k,bi,bj)+273.15d0)**1.52d0
c     &     *1.219**(9.58d0/Vb-1.12d0)*(Vb**(-0.19d0)-0.292d0)
            
            DiffusivitySurfWater (i,j) !cm2/s
     &  = 7.4d-8*(tc+273.15d0)*(2.6d0*18)**0.5d0
     &    /1d0/(Vb**0.6d0) 
            temp3 = 0.05594d0*tc**2
     &    +5.2842d0*tc+137.37d0
            
            etaW = (tc+246d0)/(temp3)
            
            SchmidtWPFCn = etaW/1.03d3/DiffusivitySurfWater (i,j)

#ifdef USE_EXFIWR            
            temp4 =  wspeed(i,j,bi,bj)
     &           *sqrt(6.1d-4+6.3d-5*wind(i,j,bi,bj))
#else
            temp4 =  wind(i,j,bi,bj)
     &           *sqrt(6.1d-4+6.3d-5*wind(i,j,bi,bj))
#endif /* USE_EXFIWR */

            temp5 = 13.3d0*sqrt(SchmidtAPFCn)
     &           +1.3d-3**(-0.5d0)-5
     &           +log(SchmidtAPFCn)/2.d0/0.4d0
            
            ka (i,j) = (1.0d-3+temp4/temp5)
#ifdef USE_EXFIWR
            temp6 = 0.222d0*wspeed(i,j,bi,bj)*wspeed(i,j,bi,bj)
     &   +0.333d0*wspeed(i,j,bi,bj)
#else            
            temp6 = 0.222d0*wind(i,j,bi,bj)*wind(i,j,bi,bj)
     &   +0.333d0*wind(i,j,bi,bj)
#endif /* USE_EXFIWR */
     
       kw = temp6*(SchmidtWPFCn/600.d0)**(-0.5d0)
            
c            KawPFCn (i,j) = 1.0d-3 ! Li et al. for PFOA
            KawPFCn (i,j) = 10.0**1.65 !  for PFOS
c==================================================================
c COSMOtherm estimated properties for netural PFCAs and PFSAs
c 2011 Environ. Chem. 2011, 8, 389–398.doi:10.1071/EN10143_AC
c Abbreviation log KAW log KOW, dry log KOA log PL (Pa) log SW (mol L–1)
c Perfluorinated carboxylate acids (PFCAs)
c PFBA –3.23 2.82 6.04 3.59 0.42
c PFPeA –2.90 3.43 6.33 3.13 –0.37
c PFHxA –2.58 4.06 6.63 2.66 –1.16
c PFHpA –2.25 4.67 6.92 2.20 –1.94
c PFOA –1.93 (–3.0) 5.30 7.23 1.73 (0.62) –2.73
c PFNA –1.58 5.92 7.50 1.27 (0.10) –3.55
c PFDA –1.27 6.50 7.77 0.82 (–0.64) –4.31
c PFUnDA –0.92 7.15 8.08 0.34 (–0.98) –5.13
c PFDoDA –0.58 7.77 8.36 –0.13 (–2.29) –5.94
c PFTrDA –0.38 8.25 8.63 –0.57 –6.59
c PFTeDA 0.03 8.90 8.87 –0.99 –7.42 
c Perfluorinated sulfonic acids (PFSAs)
c PFBS –2.59 3.90 6.49 2.80 –1.00
c PFHxS –2.38 5.17 7.55 1.77 –2.24
c PFOS –1.65 6.43 8.07 0.83 –3.92
c PFDS –1.15 7.66 8.82 –0.15 –5.39   
c================================================================
          !currently as a constant 
c                               simplifications good for PFCs but not for others
c                       Kaw changes w/ Temperature, Salinity, etc

            KawPFCeff(i,j) = KawPFCn (i,j)
     &   /(1.d0+10.d0**(pH(i,j,1,bi,bj)-pKa)) 
                        
            Kaexch (i,j) = (1. _d 0 - iceo)
     &   /(1.d0/ka(i, j)+KawPFCeff(i,j)/kw(i,j)) !m/s
            
            pfcnevaso =      Kaexch(i,j)
     &        *(pfcnconco*1.0d-9/Mb !ng/m3 to mol/m3
     &          - PTR(i,j,1,bi,bj,ipfcn)*KawPFCeff(i,j))  
     & !!mol/m2/s          
C            write(*,*)'xz check pfcnevaso:',pfcnevaso                                                                                                                   
            fluxPFCn(i,j) = -pfcnevaso +pfcnrivero !mol/m2/s
            
     
            If (i .eq. 180 .and. j .eq. 80) then
            write(*,*)'fice     at (180,80):',fice(i,j,bi,bj)
            write(*,*)'pisvel   at (180,80):',pisvel(i,j,bi,bj)
            write(*,*)'Schmidt#air at (180,80):',SchmidtAPFCn
            write(*,*)'Schmidt#water at (180,80):',SchmidtWPFCn
            write(*,*)'PFCnSurfOcean at (180,80):',
     &                     max(PTR(i,j,1,bi,bj,ipfcn),0d0)
            write(*,*)'KawPFCeff        at (180,80):',KawPFCeff(i,j)
            write(*,*)'Kaexch   at (180,80):',kaexch (i,j)
            write(*,*)'FluxPFCn no seaspray at (180,80):',fluxPFCn(i,j)
            Endif

#ifdef ALLOW_SEASPRAY
        pfcnseasprayo =  (1. _d 0 - iceo)
     &          *seaspraywtrflux(i,j,bi,bj) !m3/m2/s
     &          *enrichfactor*PTR(i,j,1,bi,bj,ipfcn)!mol/m3
C            write(*,*)'xz check enrichfactor:',enrichfactor                                             
            fluxPFCn(i,j) = fluxPFCn(i,j)-pfcnseasprayo
     &       
        If (i .eq. 180 .and. j .eq. 80) then          
        write(*,*)'FluxPFCn w/ seaspray at (180,80):',fluxPFCn(i,j)
        Endif
#endif /* ALLOW_SEASPRAY    */       


           ELSE
             fluxPFCn(i,j) = 0. _d 0
           ENDIF
           
C Calculate diagnistics
#ifdef ALLOW_DIAGNOSTICS
            pfcnconc_diag(i,j) = pfcnconco  !air concentration ng/m^3
            pfcnriver_diag(i,j) = pfcnrivero !mol/m^2/s
            pfcnevas_diag(i,j) = pfcnevaso  !ocean to air evaporations mol/m^2/s
#ifdef ALLOW_SEASPRAY
            pfcnseaspray_diag(i,j) = pfcnseasprayo
#endif /* ALLOW_SEASPRAY    */        
#endif /* ALLOW_DIAGNOSTICS */     

         ENDDO
        ENDDO

C update surface tendencies      
        DO j=jMin,jMax
          DO i=iMin,iMax
           SURFPFCn(i,j)=
     &         HFacC(i,j,1,bi,bj)*FluxPFCn(i,j)*recip_drF(1)
C       xz: drF-Cell face separation along Z axis ( units of r ).     
          ENDDO
         ENDDO

c write diagnostics
#ifdef ALLOW_DIAGNOSTICS
        WRITE(diagname,'(A8)') 'PFCnCONC'
        CALL DIAGNOSTICS_FILL( pfcnconc_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )
        WRITE(diagname,'(A8)') 'PFCnRIVR'
        CALL DIAGNOSTICS_FILL( pfcnriver_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )
#ifdef ALLOW_SEASPRAY
        WRITE(diagname,'(A8)') 'PFCnSPRY'
        CALL DIAGNOSTICS_FILL( pfcnseaspray_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )
#endif /* ALLOW_SEASPRAY    */
        WRITE(diagname,'(A8)') 'PFCnEVAS'
        CALL DIAGNOSTICS_FILL( pfcnevas_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )     
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_PFC         */
#endif /* ALLOW_PTRACERS    */

      RETURN
      END
!EOC



        
        

        
        
        
