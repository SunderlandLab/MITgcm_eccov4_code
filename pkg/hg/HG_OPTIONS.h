!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: HG_OPTIONS.h
!
! !DESCRIPTION: Header file specifying #ifdefs for the hg package
!  (i.e. pfos simulation).
!
! !REMARKS:
!  This should be placed in your verification/global_hg_llc90/code directory.
!
! !REVISION HISTORY:
!  24 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

C == choose one mode to build the code

C = inorganic mercury chemistry is the default mode

C = mercury chemistry with methylmercury 
#define METHY_MODE

C = with food web, following model of Asif Quesi
#define FOODW_MODE

C = FOODW_MODE needs METHY_MODE
#ifdef FOODW_MODE
#define METHY_MODE
#endif

C = use EXF fields, needs to turn on for online run
#define USE_EXFIWR 

C = mehg bioaccumulation
C Mason et al 1996, bioduluaiton?
#define BC_Mason1996
C Correct BCF by Hammerschmidt et al. 2013 data?
#define BC_Hammerschmidt2013

C = partition
C no partition
#undef PART_NO
C instantaneous equilibrium
#define PART_EQ
C release
#undef PART_DI

C = a seperate tracer for riverine Hg2, 
C = the refractory part of atmospheric deposition can also use this tracer
#undef HG2REFRACTORY

C - is the hg2 tracer excludely labile?
c for atm depo only case, enable this for the labil pool goto hg2, 
c      refractory portion goes to hg2r, need to enable HG2REFRACTORY
c for river case, we need hg2r as riverine hg-doc complexes, 
c      so we need to disable this triger
#undef ISHG2LABILE

C = also the fraction of refractory for atmospheric deposition?
#ifdef HG2REFRACTORY
#ifdef   ISHG2LABILE
#define    f_depo_refrac 0.6
#else
#define    f_depo_refrac 0
#endif
#else
#define   f_depo_refrac 0
#endif

C = a seperate tracer for riverine HgP
#undef HGPREFRACTORY

C what's the meaning of refractory? irreducible or slow-releasing
#ifdef HGPREFRACTORY
c this case is for the partitioning of refractory pool of atm depo
#undef HGPR_IRREDUCIBLE
c this case is for the riverine refractory hgp
#define HGPR_SLOWRELEASING
#endif

C = Walsh & Nittrouer Typology filter for riverine particles
#ifdef HGPREFRACTORY
#define WALSH_FILTER
C = disable sedimentation:
C = 0: no disable
C = 1: only river mouth
C = 2: over depth <= 55 m (top 4 levels)
C = 3: over depth <= 185 m (top 8 levels)
C = 4: over river mouth AND depth shallower than 55 m
#define walsh_depth_0
#undef walsh_depth_1
#undef walsh_depth_2
#undef walsh_depth_3
#undef walsh_depth_4
#endif

C == HgII reducible pool is a function of DOC concentrations
#undef VARI_RED_POOL

!EOC
