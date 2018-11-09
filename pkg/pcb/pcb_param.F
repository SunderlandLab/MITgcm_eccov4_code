!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pcb_param.F
!
! !DESCRIPTION: Set parameters needed for PCB runs. 
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "PCB_OPTIONS.h"

      SUBROUTINE PCB_PARAM( myThid )
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
#include "PCB.h"
!
! !INPUT PARAMETERS: 
!
      ! Number of this instance of PCB_PARAM
      INTEGER myThid   
!
! !REMARKS:
!
! !REVISION HISTORY:
!  02 Jul 2015 - H. Amos      - Update air-sea exchange equations.
!  04 Mar 2016 - R. Yantosca - Added ProTex headers
!  04 Mar 2016 - R. Yantosca - Deleted obsolete code, lined up code
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#ifdef ALLOW_PCB

      !---------------------------------------------------------
      ! Local variables
      !---------------------------------------------------------
      _RL logKow
      _RL logKoc

      _BEGIN_MASTER( myThid )
      !---------------------------------------------------------
      ! Coefficients for calculating the Schmidt number in air
      !   From Johnson (2010), Ocean Science:
      !
      !   Schmidt number in air:
      !   Sca = v_a / D_a = eta_a / ( rho_a * D_a )                [Eq. 23]
      !   
      !   dynamic viscosity:
      !   eta_a = sv_0 + sv_1*T + sv_2*T^2 + sv_3*T^3 + sv_4*T^4   [Eq. 24]
      !
      !   density:
      !   rho_a = sd_0 + sd_1*T + sd_2*T^2 + sd_3*T^3              [Eq. 25]
      !
      !   diffusion coefficient in air:
      !   Da = 0.001 * (T^1.75) * (Mr^0.5) / ( (P*Va^1/3) + Vb^1/3 )^2   [Eq. 26]
      !
      !   Mr = ( Ma + Mb ) / (Ma * Mb)                             [Eq. 27]
      !
      !   where t is temperature in degrees Celcius, the units on eta_a are
      !   kg m^-3, and the units on rho_a are kg m^-3. Values of sv_* and 
      !   sd_* are in Table A1 of Johnson (2010), which references 
      !   Tsilingiris (2008). P is the atmospheric pressure in atm. Va and
      !   Vb are the molar volumes of air and of the gas of interest, 
      !   respectively. Ma and Mb are the molecular masses of air and the
      !   gas of interest, respectively. 
      !  
      ! References:
      ! Johnson, M. T. A numerical scheme to calculate temperature and 
      !   salinity dependent air-water transfer velocities for any gas.
      !   Ocean Science, 6, 913-932, 2010.
      !
      ! Tsilingiris, P. Thermophysical and transport properties of 
      !   humid air at temperature ranges between 0 and 100 C. Energ.
      !   Convers. Manage, 49, 1098-1110, 2008.
      !
      ! Tucker, W. A. and Nelken, L. H. Diffusion coefficients in air
      !   and water, in: Handbook of Chemical Property Estimatation
      !   Methods, American Chemical Society, Washington D. C., 1990.
      !---------------------------------------------------------

      !---------------------------------------------------------
      ! Parameters from Table A1 in Johnson (2010)
      !---------------------------------------------------------
      sv_0 = 1.71575 _d -5
      sv_1 = 4.7224 _d -8
      sv_2 = -3.6630 _d -10
      sv_3 = 1.87324 _d -12
      sv_4 = -8.0502 _d -14
      sd_0 = 1.2934 _d 0 
      sd_1 = -5.5384 _d -3
      sd_2 = 3.8602 _d -5
      sd_3 = -5.2536 _d -7

      !---------------------------------------------------------     
      ! Molar volume of air (cm^3 mol^-1)
      ! Reference: Tucker & Nelken (1990) 
      !---------------------------------------------------------
      Va = 20.1 _d 0  
      
!Replace molar volumes with Schwarzenbach et al. 2003 for
!consistency w/ Freidman & Selin 2015 (in review)
!30 Sep 2015, HM Amos
!      !---------------------------------------------------------
!      ! Molar volume of gas of interest (cm^3 mol^-1)
!      !  PCB 28  : 247.3 cm3/mol
!      !  PCB 52  : 268.2 
!      !  PCB 101 : 289.1 
!      !  PCB 118 : 289.1
!      !  PCB 138 : 310
!      !  PCB 153 : 310
!      !  PCB 180 : 330.9
!      !
!      ! Reference: Shiu & Mackay (1986)
!      ! Notes: Vb at 25 degrees C. Vb for PCB 180 is not listed
!      !        in Shiu & Mackay (1986), so what's used here for 180
!      !        is the same Vb as for PCBs 171, 179, 183, 185, and 187.
!      !---------------------------------------------------------
!      !Vb = 247.3d0   ! for PCB 28
!      !Vb = 289.1d0   ! for PCB 118
!      Vb = 310.0d0   ! for PCB 153
!      !Vb = 330.9d0   ! for PCB 180

      !---------------------------------------------------------
      ! Molar volume of gas of interest (cm^3 mol^-1)
      !  PCB 28  : 169.14
      !  PCB 52  : 181.38 
      !  PCB 101 : 193.62 
      !  PCB 118 : 193.62
      !  PCB 138 : 205.86
      !  PCB 153 : 205.86
      !  PCB 180 : 218.1
      !
      ! Reference: Schwarzenbach, R. P. , Gschwend, P. M., and
      !            D. M. Imboden, Environmental Organic,
      !            Chemistry, 2nd Ed, 2003.
      !---------------------------------------------------------
      Vb = 169.14d0   ! for PCB 28
      !Vb = 193.62d0   ! for PCB 101
      !Vb = 205.86d0   ! for PCB 153
      !Vb = 218.1d0    ! for PCB 180

      !---------------------------------------------------------
      ! Molecular mass of air (g mol^-1)
      ! Reference: Tucker & Nelken (1990)
      !---------------------------------------------------------
      Ma = 28.97 _d 0 

      !---------------------------------------------------------
      ! Molecular mass of gas of interest (g mol^-1)
      !  PCB 28  : 257.54  g/mol
      !  PCB 52  : 291.99 
      !  PCB 101 : 326.43
      !  PCB 118 : 326.43
      !  PCB 138 : 360.88
      !  PCB 153 : 360.88
      !  PCB 180 : 395.32
      !
      ! Reference: Li et al. (2003), J. Phys. Chem. Ref. Data
      !---------------------------------------------------------
      Mb = 257.54d0  ! for PCB 28
      !Mb = 326.43d0  ! for PCB 101
      !Mb = 326.43d0  ! for PCB 118
      !Mb = 360.88d0  ! for PCB 153
      !Mb = 395.32d0  ! for PCB 180

      !---------------------------------------------------------
      ! Atmospheric pressure (atm). Johnson (2010) assumes P
      ! to be unity for all calculations.       
      !---------------------------------------------------------
      Patm = 1.0 _d 0 !atm
           
!Kow-Koc relationship updated to use Sobek et al. 2004 ES&T,
!and they use Kow from Hawker & Connell 1998
!05 October 2015, HM Amos
      !---------------------------------------------------------
      ! Octanol-water partition coefficient (L/kg)
      ! 
      ! *********************
      ! * OLD MODEL VERSION *
      ! *********************
      ! log Kow from Schenker et al. 2005:   
      !  PCB 28  : 5.92  
      !  PCB 52  : 6.26  
      !  PCB 101 : 6.76  
      !  PCB 118 : 7.08  
      !  PCB 138 : 7.70 
      !  PCB 153 : 7.31 
      !  PCB 180 : 7.66 
      !
      ! *********************
      ! *    USE THESE !    *
      ! *********************
      ! log Kow from Hawker & Connell 1988:
      !  PCB 28  : 5.67  
      !  PCB 52  : 5.84  
      !  PCB 101 : 6.38  
      !  PCB 118 : 6.74  
      !  PCB 138 : 6.83 
      !  PCB 153 : 6.92 
      !  PCB 180 : 7.36 
      !
      ! References: 
      ! Hawker & Connell (1988), Environ. Sci. Technol., Table 2
      ! Schenker et al. (2005), Environ. Sci. Technol., Table 2
      !---------------------------------------------------------
      logKow = 5.67d0   ! for PCB-28
      !logKow = 6.38d0   ! for PCB-101
      !logKow = 6.92d0   ! for PCB-153
      !logKow = 7.36d0   ! for PCB-180

!Update Koc-Kow relationship
!05 Oct 2015, HM Amos
      !---------------------------------------------------------
      ! Partition coefficient to organic carbon (L/kg)
      !
      ! *********************
      ! * OLD MODEL VERSION *
      ! *********************
      ! Koc from Green et al. for soil-sediment:   
      !    Koc = 0.33 * Kow
      !
      ! *********************
      ! *     USE THIS !    *
      ! *********************
      ! log Koc from Sobek et al. 2004:
      !    logKow = 0.88(+/-0.07)logKow + 0.90(+/-0.47)
      ! where the (+/-) are 95% confidence intervals      
      !
      ! References: 
      ! Greene et al. (2013), Environ Sci Technol
      ! Sobek et al. (2004), Particle-water partitioning of PCBs
      !   in the photic zone: a 25-month study in the open 
      !   Baltic Sea, Environ Sci Technol
      !---------------------------------------------------------
      !Koc_pcb = 0.33d0 * Kow
      logKoc  = 0.88d0 * logKow + 0.90d0
      Koc_pcb = 10d0**logKoc

!Degradation half lives replaced by Mackay et al. 2006, for
!consistency w/ Carey Friedman GEOS-Chem & Matt McCleod BETR-GLOBAL
!01 October 2015, HM Amos
!      !---------------------------------------------------------
!      ! Degradation half-lives (hr)
!      ! 
!      !  PCB 28  :   1,450 hrs
!      !  PCB 52  :  30,000  
!      !  PCB 101 :  60,000
!      !  PCB 118 :  60,000
!      !  PCB 138 : 120,000
!      !  PCB 153 : 120,000
!      !  PCB 180 : 240,000
!      !
!      ! Reference: Sinkkonen & Paasivirta (2000), Table 2, 
!      !            Chemosphere
!      !---------------------------------------------------------
!      !tHalf     = 1450d0          ! for PCB 28
!      !tHalf     = 60000d0         ! for PCB 118
!      tHalf     = 120000d0        ! for PCB 153
!      !tHalf     = 240000d0        ! for PCB 180

      !---------------------------------------------------------
      ! Degradation half-lives (hr)
      ! 
      !  PCB 28  : 17,000
      !  PCB 52  : 55,000  
      !  PCB 101 : 55,000
      !  PCB 118 : 55,000
      !  PCB 138 : 55,000
      !  PCB 153 : 55,000
      !  PCB 180 : 55,000
      !
      ! Reference: Mackay et al. (2006), CRC Handbook
      !---------------------------------------------------------
      tHalf     = 17d3          ! for PCB 28
      !tHalf     = 55d3          ! for PCB 101, 153, 180
      
!Adopt temperature dependence used by Carey Friedman in GEOS-Chem
!and Matt McCleod in BETR-GLOBAL. Put it in pcb_chemistry.F
!01 October 2015, HM Amos
!
!     !---------------------------------------------------------
!      ! Temperature dependent degradation
!      ! Kdeg(T) = a * exp(b*T)
!      !
!      ! Kdeg : degradation rate coefficient (1/sec)
!      ! a    : congener dependent coefficient (1/sec)
!      ! b    : constant coefficient (1/C)
!      ! T    : seawater temperature (degrees C)
!      !
!      ! The value of "a" depends on the congener:
!      !  PCB 28  :    8.174d-8 1/sec 
!      !  PCB 52  :   0.3951d-8   
!      !  PCB 101 :   0.1975d-8
!      !  PCB 118 :   0.1975d-8
!      !  PCB 138 :  0.09877d-8
!      !  PCB 153 :  0.09877d-8
!      !  PCB 180 :  0.04938d-8
!      !---------------------------------------------------------
!      b_deg = 0.06931d0
!
!      !a_deg = 8.174d-8 ! for PCB-28
!      !a_deg = 0.1975d-8 ! for PCB-118
!      a_deg = 0.09877d-8 ! for PCB-153
!      !a_deg = 0.04938d-8 ! for PCB-180
 
c Do not need a Kiw for PCBs (hma, 21 july 2014)
c      !---------------------------------------------------------
c      ! Interface/water partition coefficient (m)
c      !
c      ! At dilute condition of the Langmuir adsorption equation.
c      ! Reference: Vecitis et al. 2008 JPC-C
c      !---------------------------------------------------------
c      Kiw = 5.1d-6*1.97d0   !xz PFOS
 
c Disable sea spray (H. Amos, )
c#ifdef ALLOW_SEASPRAY
c      enrichfactor = 5.6d0 !for PFOA MCMURDO et al. Environ. Sci. Technol. 2008, 42, 3969–3974
cc        write(*,*)'xz check enrichfactor in pcb_param:',enrichfactor 
c#endif /* ALLOW_SEASPRAY */

      !---------------------------------------------------------
      ! Ehthalpy of air-water exchange (kJ mol^-1)
      !
      !  PCB 28  : 51.8
      !  PCB 52  : 53.8
      !  PCB 101 : 65.2 
      !  PCB 118 : 65.2
      !  PCB 138 : 64.7
      !  PCB 153 : 68.2
      !  PCB 180 : 69.0
      !
      ! Factor of 1000 is to convert from kJ/mol to J/mol
      !
      ! Reference: 
      ! Schenker et al. (2005), Environ. Sci. Technol., Table 2
      !---------------------------------------------------------
      dH   = 51.8d0 * 1000.0d0 ! for PCB 28
      !dH   = 65.2d0 * 1000.0d0 ! for PCB 101
      !dH   = 65.2d0 * 1000.0d0 ! for PCB 118 
      !dH   = 68.2d0 * 1000.0d0 ! for PCB 153
      !dH   = 69.0d0 * 1000.0d0 ! for PCB 180 

      !---------------------------------------------------------
      ! Dimensionless gas-over-liquid Henry's law constant
      !
      ! log10(KH)
      !           
      !  PCB 28  : -1.93
      !  PCB 52  : -1.96
      !  PCB 101 : -2.08  
      !  PCB 118 : -2.36
      !  PCB 138 : -1.97
      !  PCB 153 : -2.13
      !  PCB 180 : -2.15
      ! 
      ! Note: Transform KH from log10 by 10**(KH)
      !
      ! References:
      !  Schenker et al. (2005), Environ. Sci. Technol., Table 2
      !---------------------------------------------------------
      KH0 = 10.d0**-1.93d0      ! for PCB 28
      !KH0 = 10.d0**-2.08d0      ! for PCB 101
      !KH0 = 10.d0**-2.36d0      ! for PCB 118
      !KH0 = 10.d0**-2.13d0      ! for PCB 153
      !KH0 = 10.d0**-2.15d0      ! for PCB 180

      _END_MASTER( myThid )

       ! Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_PCB */

      RETURN
      END
!EOC
