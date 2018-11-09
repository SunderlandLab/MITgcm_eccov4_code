!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: PCB_OPTIONS.h
!
! !DESCRIPTION: Contains switches to turn on/off various options
!  for the PCB simulation.
!\\
!\\
! !REMARKS:
!  To activate    an option, use #define OPTION
!  To de-activate an option, use #undef OPTION
!
! !REVISION HISTORY:
!  14 Jul 2014 - hma         - Copied from Xianming's PFC simulation
!                              for my PCB simulation.
!  02 Feb 2015 - hma         - Disable sea spray option.
!  04 Feb 2015 - hma         - Adopt options to mimic Yanxu's Hg river
!                              inputs for PCBs. 
!  29 Feb 2016 - R. Yantosca - Added ProTex headers
!  02 Mar 2016 - R. Yantosca - Remove partitioning options
!EOP
!------------------------------------------------------------------------------
!BOC

!==============================================================================
! %%%%%             ! A seperate tracer for riverine PCBp                %%%%%
!==============================================================================
#undef PCBPRIVER

!==============================================================================
! %%%%%             A separate tracer for refractory PCBa,               %%%%%
! %%%%%           either from river or atmospheric deposition            %%%%%
!==============================================================================
#undef PCBAREFRACTORY

#ifdef PCBAREFRACTORY
# define f_depo_refrac 0.6
#else
# define f_depo_refrac 0
#endif

!==============================================================================
! %%%%%             A separate tracer for refractory PCBp,               %%%%%
! %%%%%           this is only contributed from river source             %%%%%
!==============================================================================
#undef PCBPREFRACTORY

! what's the meaning of refractory? irreducible or slow-releasing
#ifdef PCBPREFRACTORY
#define PCBPR_IRREDUCIBLE
#undef PCBPR_SLOWRELEASING


#if defined(PCBPR_IRREDUCIBLE) & defined(PCBPR_SLOWRELEASING)
PCBPR_SLOWRELEASING and PCBPR_IRREDUCIBLE can not be turned on at the same time
#endif
#endif

! = Walsh & Nittrouer Typology filter for riverine particles
#ifdef PCBPREFRACTORY

#define WALSH_FILTER
! = disable sedimentation:
! = 0: no disable
! = 1: only river mouth
! = 2: over depth <= 55 m (top 4 levels)
! = 3: over depth <= 185 m (top 8 levels)
#define walsh_depth_0
#undef walsh_depth_1
#undef walsh_depth_2
#undef walsh_depth_3

#endif

! == HgII reducible pool is a function of DOC concentrations
#undef VARI_RED_POOL
!EOC
