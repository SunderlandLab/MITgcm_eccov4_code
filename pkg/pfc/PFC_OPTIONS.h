!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: PFC_OPTIONS.h
!
! !DESCRIPTION: Header file specifying #ifdefs for the pfc package
!  (i.e. pfos simulation).
!
! !REMARKS:
!  This should be placed in your verification/pfos/code directory.
!
! !REVISION HISTORY:
!  24 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

C use external forcing for seaice, windspeed and radiation of short wave
#define USE_EXFIWR

C chemical transport from seawater to the air via waterdroplet
#undef ALLOW_SEASPRAY

!EOC
