!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: PCB_SIZE.h
!
! !DESCRIPTION: Contains tracer indices for the PCB simulation.
!\\
!\\
! !REMARKS:
!  We are not using food web model for PCBs (H. Amos, 04 Feb 2015)
!
! !REVISION HISTORY:
!  04 Feb 2015 - H. Amos     - Copied from Yanxu Zhang's Hg model
!                              and adapted for PCBs.
!  29 Feb 2016 - R. Yantosca - Added ProTex headers; cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC

      ! PCBa is tracer #1
      INTEGER    iPCBa
      PARAMETER( iPCBa  = 1 )

      ! PCBp is tracer #2
      INTEGER    iPCBp
      PARAMETER( iPCBp  = 2 )

#ifdef PCBPREFRACTORY

      ! Number of rivers
      INTEGER    nrmax
      PARAMETER( nrmax  = 1 )

#ifdef PCBAREFRACTORY

      ! PCBaR is tracer #3 (when PCBAREFRACTORY is defined)
      INTEGER    iPCBaR
      PARAMETER( iPCBaR = 3 )

      ! PCBpR is tracer #4 (when PCBAREFRACTORY is defined)
      INTEGER    iPCBpR
      PARAMETER( iPCBpR = 4 )
#else

      ! PCBpR is tracer #3 (when PCBAREFRACTORY is NOT defined)
      INTEGER    iPCBpR
      PARAMETER( iPCBpR = 3 )

#endif

#endif
!EOC
