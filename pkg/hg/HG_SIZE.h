
!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HG_SIZE.h
!
! !DESCRIPTION: Contains size dimensions and tracer indices for
!  the Hg simulation.
!\\
!\\
! !REMARKS:

! !REVISION HISTORY:
!  29 Feb 2016 - R. Yantosca - Added ProTex headers
!  29 Feb 2016 - R. Yantosca - Consolidated some #ifdefs for clarity
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !===================================================================
      ! %%%%%                   SIZE PARAMETERS                     %%%%%
      !===================================================================
#ifdef FOODW_MODE      

      ! npmax = no of "functional groups" of phytoplankton
      INTEGER    npmax
      PARAMETER( npmax   = 6              )

      ! nzmax = no of "functional groups" of zooplankton
      INTEGER    nzmax
      PARAMETER( nzmax   = 2              )

#endif /* FOODW_MODE */

#ifdef HGPREFRACTORY
      INTEGER    nrmax
      PARAMETER( nrmax   = 1              )
#endif /* HGPREFRACTORY */

      !===================================================================
      ! %%%%%                 Hg SPECIES INDICES                    %%%%%
      !===================================================================

      ! Hg0 is tracer #1
      INTEGER     iHg0
      PARAMETER ( iHg0   = 1              )

      ! Hg2 is tracer #2
      INTEGER     iHg2
      PARAMETER ( iHg2   = 2              )

      ! HgP is tracer #3
      INTEGER     iHgP
      PARAMETER ( iHgP   = 3              )

#ifdef HG2REFRACTORY
      ! Hg2R is tracer #4 (if HG2REFRACTORY is defined)
      INTEGER     iHg2R
      PARAMETER ( iHg2R  = 4              )
#endif /* HG2REFRACTORY */

#ifdef HGPREFRACTORY
      ! iHgPR is tracer #5, if HG2REFRACTORY is defined
      ! or is tracer #4 otherwise
      INTEGER     iHgPR
#ifdef HG2REFRACTORY
      PARAMETER ( iHgPR  = 5              )
#else
      PARAMETER ( iHgPR  = 4              )
#endif /* HG2REFRACTORY */

#endif /* HGPREFRACTORY */
      
#ifdef METHY_MODE
      ! Compute tracer index for methyl mercury
      INTEGER     iMMHg
      
#ifdef HGPREFRACTORY
      PARAMETER ( iMMHg  = iHgPR + nrmax  )
#else
#ifdef HG2REFRACTORY
      PARAMETER ( iMMHg  = iHg2R + 1      )
#else
      PARAMETER ( iMMHg  = 4              )
#endif /* HG2REFRACTORY */
#endif /* HGPREFRACTORY */           

      INTEGER     iMMHgP
      PARAMETER ( iMMHgP = iMMHg  + 1     )

      INTEGER     iDMHg
      PARAMETER ( iDMHg  = iMMHgP + 1     ) 

#ifdef FOODW_MODE   
      ! remember to bring the fields in data.ptracers in the right order !

      ! iPhy  :: index of first phytoplankton
      INTEGER     iPhy
      PARAMETER ( iPhy   = iDMHg  + 1     )

      ! iZoo  :: index of first zooplankton
      INTEGER     iZoo
      PARAMETER ( iZoo   = iPhy   + npmax )

#endif /* FOODW_MODE */

#endif /* METHY_MODE */
!EOC
