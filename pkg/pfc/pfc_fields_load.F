!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pfc_fields_load.F
!
! !DESCRIPTION: Loads data arrays for the PFC simulation.
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "PFC_OPTIONS.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_OPTIONS.h"
#endif /* ALLOW_SEAICE */
#ifdef ALLOW_THSICE
# include "THSICE_OPTIONS.h"
#endif /* ALLOW_THSICE */

      SUBROUTINE PFC_FIELDS_LOAD (
     I          myIter, myTime , myThid )

!
! !USES:
!
      IMPLICIT NONE

      ! Global variables
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "PFC_SIZE.h"
#include "PFC.h"
#ifdef ALLOW_EXF
# include "EXF_FIELDS.h"
#endif
#ifdef ALLOW_SEAICE
# include "SEAICE_SIZE.h"
# include "SEAICE.h"
#endif
#ifdef ALLOW_THSICE
# include "THSICE_VARS.h"
#endif /* ALLOW_THSICE */
!
! !INPUT PARAMETERS: 
!
      ! myTime  :: current time in simulation
      ! myIter  :: current iteration number
      !CmyThid  :: my Thread Id number
      _RL myTime
      INTEGER myIter
      INTEGER myThid
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Feb 2016 - R. Yantosca - Added ProTex headers
!  29 Feb 2016 - R. Yantosca - Removed #ifdef FOODW_MODE blocks
!  29 Feb 2016 - R. Yantosca - Removed commented-out code
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

c===========from hg code=========
      INTEGER intimeP, intime0, intime1
      INTEGER bi, bj, i, j, k
      INTEGER prec
      LOGICAL nowRead

      _RL aWght, bWght
c===========from hg code=========



      IF ( PFC_forcingCycle .GT. 0. _d 0 ) THEN

C First call requires that we initialize everything to zero for safety
cQQQ need to check timing
       IF ( myIter .EQ. nIter0 ) THEN
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
            PFC_ldRec_forcing(bi,bj) = 0
c xz changed from original PFC_ldRec(bi,bj) in cfc to PFC_ldRec_forcing(bi,bj) = 0
c xz to accomodate PFC_ldRec_chem(bi,bj) later on      
          ENDDO
         ENDDO
       ENDIF

C--   Now calculate whether it is time to update the forcing arrays
       
       CALL GET_PERIODIC_INTERVAL(
     O                   intimeP, intime0, intime1, bWght, aWght,
     I                   PFC_forcingCycle, PFC_forcingPeriod,
     I                   deltaTclock, myTime, myThid )
       
       bi = myBxLo(myThid)
       bj = myByLo(myThid)

#ifdef ALLOW_DEBUG
       IF ( debugLevel.GE.debLevB ) THEN
        _BEGIN_MASTER(myThid)
        WRITE(standardMessageUnit,'(A,I10,A,4I5,A,2F14.10)')
     &   ' PFC_FIELDS_LOAD,', myIter,
     &   ' : iP,iLd,i0,i1=', intimeP,PFC_ldRec_forcing(bi,bj), 
     &                          intime0,intime1,
     &   ' ; Wght=', bWght, aWght
        _END_MASTER(myThid)
       ENDIF
#endif /* ALLOW_DEBUG */


c#ifdef ALLOW_CAL
c       IF ( nowRead ) THEN
c#else 


#ifdef ALLOW_AUTODIFF_TAMC
C-    assuming that we call S/R PFC_FIELDS_LOAD at each time-step and
C     with increasing time, this will catch when we need to load new records;
C     But with Adjoint run, this is not always the case => might end-up using
C     the wrong time-records
       IF ( intime0.NE.intimeP .OR. myIter.EQ.nIter0 ) THEN
#else /* ALLOW_AUTODIFF_TAMC */
C-    Make no assumption on sequence of calls to PFC_FIELDS_LOAD ;
C     This is the correct formulation (works in Adjoint run).
C     Unfortunatly, might produce many recomputations <== not used until it is fixed
       IF ( intime1.NE.PFC_ldRec_forcing(bi,bj) ) THEN         
#endif /* ALLOW_AUTODIFF_TAMC */

C--   If the above condition is met then we need to read in
C     data for the period ahead and the period behind myTime.
        IF ( debugLevel.GE.debLevZero ) THEN
         _BEGIN_MASTER(myThid)
         WRITE(standardMessageUnit,'(A,I10,A,2(2I5,A))')
     &    ' PFC_FIELDS_LOAD, it=', myIter,
     &    ' : Reading new data, i0,i1=', intime0, intime1,
     &    ' (prev=', intimeP, PFC_ldRec_forcing(bi,bj), ' )'
         _END_MASTER(myThid)
        ENDIF

        _BARRIER

        IF ( PFC_windFile .NE. ' '  .AND. .NOT.useEXF ) THEN
         CALL READ_REC_XY_RS( PFC_windFile,wind0,intime0,
     &        myIter,myThid )
         CALL READ_REC_XY_RS( PFC_windFile,wind1,intime1,
     &        myIter,myThid )
        ENDIF
        
         _EXCH_XY_RS(wind0, myThid )
         _EXCH_XY_RS(wind1, myThid )

C-    save newly loaded time-record
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           PFC_ldRec_forcing(bi,bj) = intime1
         ENDDO
        ENDDO

C--   end if-block for loading new time-records
       ENDIF


c==================from hg code===
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx

c           IF ( PFC_windFile .NE. ' '  ) THEN
c             wind(i,j,bi,bj) = bWght*wind0(i,j,bi,bj)
c     &                           +aWght*wind1(i,j,bi,bj)
c       xz: change for ecco4 
           IF ( PFC_windFile .NE. ' '  ) THEN
             wind = bWght*wind0(i,j,bi,bj)
     &                           +aWght*wind1(i,j,bi,bj)

           ELSE
             wind = 5. _d 0*maskC(i,j,1,bi,bj)
           ENDIF
c calculate piston velocity
c QQ: note - we should have wind speed variance in here
c following Wannikof (1992)
c m/s
c F = - D/Zfilm (Cg - Cl) = - D/Zfilm (KHPg - Cl) 
c Because D/Zfilm has velocity units, it has been called the Piston Velocity 
c piston velocity calculation is different for Hg0
c according to Soerensen et al. 2010
c yzh, 6/18/2013
c xz: may need to be changed for POPs if gasesou air-sea exchange is important
c xz: for PFC, this process is not important
           pisvel(i,j,bi,bj)=(0.31 _d 0*wind(i,j,bi,bj)**2)/3.6 _d 5
          ENDDO
         ENDDO
        ENDDO
       ENDDO
c===============================

c==================xz pfc seaspray start ===
#ifdef ALLOW_SEASPRAY
c  xz: numberical integration able to adjust steps 
C               get values of seaspraywtrflux in m3/m2/s
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx

c       write(*,*) 'myIter, myTime, myThid', myIter,myTime,myThid

           
            Call PFC_CalcWaterSpray( 
     I                          i, j, bi, bj) 

           
      
            if (i .eq. 180 .and. j .eq. 80) then
             write(*,*)'seaspraywtrflux at (180,80):'
     &          ,seaspraywtrflux (i,j,bi,bj)
            endif

          ENDDO
         ENDDO
        ENDDO
       ENDDO
#endif
c===============xz pfc seaspray end================

C endif for periodicForcing
      ENDIF

C=======================================================================
C read fields for PFC chemistry,based on Hg code by yanxu
C in yanxu's Hg code,period cycle for Hg and for Darwin are the same 
C So  PFC_chemForcing, PFC_chemCycle,
C Hg_ldRec_chem(bi, bj) were target for both Hg and Darwin input
C Here, for PFC input, I want to use PFC cycling period of 1 year and
C cycle of multiple years. So I added PFC_emisForcing,
C PFC_emisCycle, PFC_ldRec_emis(bi, bi) to stand for the  original _chem
C for PFC emission (riverine discharge, atmos deposition, etc)
C=======================================================================
      IF ( PFC_chemCycle .GT. 0. _d 0 ) THEN

C First call requires that we initialize everything to zero for safety
       IF ( myIter .EQ. nIter0 ) THEN
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
            PFC_ldRec_chem(bi,bj) = 0
          ENDDO
         ENDDO
       ENDIF

C--   Now calculate whether it is time to update the forcing arrays
c#ifdef ALLOW_CAL
c       CALL CAL_GET_PERIODIC_INTERVAL(
c     O                   intimeP, intime0, intime1, bWght, aWght,
c     I                   PFC_chemCycle, PFC_chemPeriod, nowRead,
c     I                   deltaTclock, myIter, myTime, myThid )
c        if ( myIter .eq. nIter0 ) THEN
c            nowRead=.TRUE.
c        endif    
c#else
       CALL GET_PERIODIC_INTERVAL(
     O                   intimeP, intime0, intime1, bWght, aWght,
     I                   PFC_chemCycle, PFC_chemPeriod,
     I                   deltaTclock, myTime, myThid )
c#endif

       bi = myBxLo(myThid)
       bj = myByLo(myThid)
#ifdef ALLOW_DEBUG
       IF ( debugLevel.GE.debLevB ) THEN
        _BEGIN_MASTER(myThid)
        WRITE(standardMessageUnit,'(A,I10,A,4I5,A,2F14.10)')
     &   ' PFC_FIELDS_LOAD,', myIter,
     &   ' : iP,iLd,i0,i1=', intimeP,PFC_ldRec_chem(bi,bj), 
     &                       intime0,intime1,
     &   ' ; Wght=', bWght, aWght
        _END_MASTER(myThid)
       ENDIF
#endif /* ALLOW_DEBUG */

c#ifdef ALLOW_CAL
c       IF ( nowRead ) THEN
c#else
#ifdef ALLOW_AUTODIFF_TAMC
C-    assuming that we call S/R HG_FIELDS_LOAD at each time-step and
C     with increasing time, this will catch when we need to load new records;
C     But with Adjoint run, this is not always the case => might end-up using
C     the wrong time-records
       IF ( intime0.NE.intimeP .OR. myIter.EQ.nIter0 ) THEN
#else /* ALLOW_AUTODIFF_TAMC */
C-    Make no assumption on sequence of calls to HG_FIELDS_LOAD ;
C     This is the correct formulation (works in Adjoint run).
C     Unfortunatly, might produce many recomputations <== not used until it is fixed
       IF ( intime1.NE.PFC_ldRec_chem(bi,bj) ) THEN
#endif /* ALLOW_AUTODIFF_TAMC */
c#endif
C--   If the above condition is met then we need to read in
C     data for the period ahead and the period behind myTime.
        IF ( debugLevel.GE.debLevZero ) THEN
         _BEGIN_MASTER(myThid)
         WRITE(standardMessageUnit,'(A,I10,A,2(2I5,A))')
     &    ' PFC_FIELDS_LOAD, it=', myIter,
     &    ' : Reading new data, i0,i1=', intime0, intime1,
     &    ' (prev=', intimeP, PFC_ldRec_chem(bi,bj), ' )'
         _END_MASTER(myThid)
        ENDIF

        _BARRIER

        IF ( radsw_File .NE. ' ' ) THEN
         CALL READ_REC_XY_RS( radsw_File,radsw0,intime0,
     &        myIter,myThid )
         CALL READ_REC_XY_RS( radsw_File,radsw1,intime1,
     &        myIter,myThid )
        ENDIF
        IF ( PFC_iceFile .NE. ' '  .AND. .NOT.useSEAICE ) THEN
         CALL READ_REC_XY_RS( PFC_iceFile,ice0,intime0,
     &       myIter,myThid )
         CALL READ_REC_XY_RS( PFC_iceFile,ice1,intime1,
     &       myIter,myThid )
        ENDIF
        IF ( PFCi_riverFile .NE. ' ' ) THEN
        CALL READ_REC_XY_RS( PFCi_riverFile,pfciriver0,intime0,
     &       myIter,myThid )
        CALL READ_REC_XY_RS( PFCi_riverFile,pfciriver1,intime1,
     &       myIter,myThid )
        ENDIF
                
        prec=32
        IF ( chl_File .NE. ' ' ) THEN
         CALL READ_REC_3D_RS( chl_File,prec,Nr,chl0,intime0,
     &        myIter,myThid )
         CALL READ_REC_3D_RS( chl_File,prec,Nr,chl1,intime1,
     &        myIter,myThid )
        ENDIF
        IF ( npp_File .NE. ' ' ) THEN
         CALL READ_REC_3D_RS( npp_File,prec,Nr,npp0,intime0,
     &        myIter,myThid )
         CALL READ_REC_3D_RS( npp_File,prec,Nr,npp1,intime1,
     &        myIter,myThid )
        ENDIF
        IF ( doc_file .NE. ' ' ) THEN
         CALL READ_REC_3D_RS( doc_file,prec,Nr,doc0,intime0,
     &        myIter,myThid )
         CALL READ_REC_3D_RS( doc_file,prec,Nr,doc1,intime1,
     &        myIter,myThid )
        ENDIF
        IF ( poc_file .NE. ' ' ) THEN
         CALL READ_REC_3D_RS( poc_file,prec,Nr,poc0,intime0,
     &        myIter,myThid )
         CALL READ_REC_3D_RS( poc_file,prec,Nr,poc1,intime1,
     &        myIter,myThid )
        ENDIF
        IF ( Fdoc_file .NE. ' ' ) THEN
         CALL READ_REC_3D_RS( Fdoc_file,prec,Nr,Fdoc0,intime0,
     &        myIter,myThid )
         CALL READ_REC_3D_RS( Fdoc_file,prec,Nr,Fdoc1,intime1,
     &        myIter,myThid )
        ENDIF
        IF ( Fpoc_file .NE. ' ' ) THEN
         CALL READ_REC_3D_RS( Fpoc_file,prec,Nr,Fpoc0,intime0,
     &        myIter,myThid )
         CALL READ_REC_3D_RS( Fpoc_file,prec,Nr,Fpoc1,intime1,
     &        myIter,myThid )
        ENDIF
        IF ( Fremin_file .NE. ' ' ) THEN
         CALL READ_REC_3D_RS( Fremin_file,prec,Nr,Fremin0,intime0,
     &        myIter,myThid )
         CALL READ_REC_3D_RS( Fremin_file,prec,Nr,Fremin1,intime1,
     &        myIter,myThid )
        ENDIF
        IF ( o2_file .NE. ' ' ) THEN
         CALL READ_REC_3D_RS( o2_file,prec,Nr,o20,intime0,
     &        myIter,myThid )
         CALL READ_REC_3D_RS( o2_file,prec,Nr,o21,intime1,
     &        myIter,myThid )
        ENDIF
        IF ( dfe_file .NE. ' ' ) THEN
         CALL READ_REC_XY_RS( dfe_file,dfe0,intime0,
     &        myIter,myThid )
         CALL READ_REC_XY_RS( dfe_file,dfe1,intime1,
     &        myIter,myThid )
        ENDIF

         _EXCH_XY_RS(radsw0, myThid )
         _EXCH_XY_RS(radsw1, myThid )
        IF (.NOT.useSEAICE) THEN
         _EXCH_XY_RS(ice0, myThid )
         _EXCH_XY_RS(ice1, myThid )
        ENDIF

         _EXCH_XYZ_RS(chl0, myThid )
         _EXCH_XYZ_RS(chl1, myThid )
         _EXCH_XYZ_RS(npp0, myThid )
         _EXCH_XYZ_RS(npp1, myThid )
         _EXCH_XYZ_RS(doc0, myThid )
         _EXCH_XYZ_RS(doc1, myThid )
         _EXCH_XYZ_RS(poc0, myThid )
         _EXCH_XYZ_RS(poc1, myThid )
         _EXCH_XYZ_RS(Fdoc0, myThid )
         _EXCH_XYZ_RS(Fdoc1, myThid )
         _EXCH_XYZ_RS(Fpoc0, myThid )
         _EXCH_XYZ_RS(Fpoc1, myThid )
         _EXCH_XYZ_RS(Fremin0, myThid )
         _EXCH_XYZ_RS(Fremin1, myThid )
C-    save newly loaded time-record
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           PFC_ldRec_chem(bi,bj) = intime1
         ENDDO
        ENDDO

C--   end if-block for loading new time-records
       ENDIF

C
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           IF ( radsw_File .NE. ' '  ) THEN
             radsw(i,j,bi,bj) = bWght*radsw0(i,j,bi,bj)
     &                          +aWght*radsw1(i,j,bi,bj)
           ELSE
             radsw(i,j,bi,bj) = 0. _d 0
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           IF ( PFC_iceFile .NE. ' '  ) THEN
            FIce(i,j,bi,bj)    = bWght*ice0(i,j,bi,bj)
     &                          +aWght*ice1(i,j,bi,bj)
           ELSE
            FIce(i,j,bi,bj) =0. _d 0
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           IF ( chl_File .NE. ' '  ) THEN
             chl(i,j,k,bi,bj) = bWght*chl0(i,j,k,bi,bj)
     &                          +aWght*chl1(i,j,k,bi,bj)
           ELSE
             chl(i,j,k,bi,bj) = 0. _d 0
           ENDIF
           IF ( npp_File .NE. ' '  ) THEN
             npp(i,j,k,bi,bj) = bWght*npp0(i,j,k,bi,bj)
     &                          +aWght*npp1(i,j,k,bi,bj)
           ELSE
             npp(i,j,k,bi,bj) = 0. _d 0
           ENDIF
           IF ( doc_File .NE. ' '  ) THEN
             doc(i,j,k,bi,bj) = bWght*doc0(i,j,k,bi,bj)
     &                          +aWght*doc1(i,j,k,bi,bj)
           ELSE
             doc(i,j,k,bi,bj) = 0. _d 0
           ENDIF
           IF ( poc_File .NE. ' '  ) THEN
             poc(i,j,k,bi,bj) = bWght*poc0(i,j,k,bi,bj)
     &                          +aWght*poc1(i,j,k,bi,bj)
           ELSE
             poc(i,j,k,bi,bj) = 0. _d 0
           ENDIF
           IF ( Fdoc_File .NE. ' '  ) THEN
             Fdoc(i,j,k,bi,bj) = bWght*Fdoc0(i,j,k,bi,bj)
     &                          +aWght*Fdoc1(i,j,k,bi,bj)
           ELSE
             Fdoc(i,j,k,bi,bj) = 0. _d 0
           ENDIF
           IF ( Fpoc_File .NE. ' '  ) THEN
             Fpoc(i,j,k,bi,bj) = bWght*Fpoc0(i,j,k,bi,bj)
     &                          +aWght*Fpoc1(i,j,k,bi,bj)
           ELSE
             Fpoc(i,j,k,bi,bj) = 0. _d 0
           ENDIF
           IF ( Fremin_File .NE. ' '  ) THEN
             Fremin(i,j,k,bi,bj) = bWght*Fremin0(i,j,k,bi,bj)
     &                          +aWght*Fremin1(i,j,k,bi,bj)
           ELSE
             Fremin(i,j,k,bi,bj) = 0. _d 0
           ENDIF
           IF ( o2_File .NE. ' '  ) THEN
             o2(i,j,k,bi,bj) = bWght*o20(i,j,k,bi,bj)
     &                          +aWght*o21(i,j,k,bi,bj)
           ELSE
             o2(i,j,k,bi,bj) = 0. _d 0
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C endif for periodicChem !for IF ( PFC_chemCycle .GT. 0. _d 0 ) THEN
      ENDIF 

C=========================================================================
C xz for read in PFC input (emission) files including riverine discharge 
C atmospheric deposition, atmospheric concentration...
C=========================================================================

      IF ( PFC_emisCycle .GT. 0. _d 0 ) THEN

C First call requires that we initialize everything to zero for safety
       IF ( myIter .EQ. nIter0 ) THEN
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
            PFC_ldRec_emis(bi,bj) = 0
          ENDDO
         ENDDO
       ENDIF

       CALL GET_PERIODIC_INTERVAL(
     O                   intimeP, intime0, intime1, bWght, aWght,
     I                   PFC_emisCycle, PFC_emisPeriod,
     I                   deltaTclock, myTime, myThid )

       bi = myBxLo(myThid)
       bj = myByLo(myThid)
#ifdef ALLOW_DEBUG
       IF ( debugLevel.GE.debLevB ) THEN
        _BEGIN_MASTER(myThid)
        WRITE(standardMessageUnit,'(A,I10,A,4I5,A,2F14.10)')
     &   ' PFC_FIELDS_LOAD,', myIter,
     &   ' : iP,iLd,i0,i1=', intimeP,PFC_ldRec_emis(bi,bj), 
     &                       intime0,intime1,
     &   ' ; Wght=', bWght, aWght
        _END_MASTER(myThid)
       ENDIF
#endif /* ALLOW_DEBUG */

#ifdef ALLOW_CAL
       IF ( nowRead ) THEN
#else
#ifdef ALLOW_AUTODIFF_TAMC
C-    assuming that we call S/R HG_FIELDS_LOAD at each time-step and
C     with increasing time, this will catch when we need to load new records;
C     But with Adjoint run, this is not always the case => might end-up using
C     the wrong time-records
       IF ( intime0.NE.intimeP .OR. myIter.EQ.nIter0 ) THEN
#else /* ALLOW_AUTODIFF_TAMC */
C-    Make no assumption on sequence of calls to HG_FIELDS_LOAD ;
C     This is the correct formulation (works in Adjoint run).
C     Unfortunatly, might produce many recomputations <== not used until it is fixed
       IF ( intime1.NE.PFC_ldRec_emis(bi,bj) ) THEN
#endif /* ALLOW_AUTODIFF_TAMC */
#endif
C--   If the above condition is met then we need to read in
C     data for the period ahead and the period behind myTime.
        IF ( debugLevel.GE.debLevZero ) THEN
         _BEGIN_MASTER(myThid)
         WRITE(standardMessageUnit,'(A,I10,A,2(2I5,A))')
     &    ' PFC_FIELDS_LOAD, it=', myIter,
     &    ' : Reading new data, i0,i1=', intime0, intime1,
     &    ' (prev=', intimeP, PFC_ldRec_chem(bi,bj), ' )'
         _END_MASTER(myThid)
        ENDIF

        _BARRIER

 
        IF ( PFCn_concFile .NE. ' ' ) THEN
         CALL READ_REC_XY_RS( PFCn_concFile,pfcnconc0,intime0,
     &       myIter,myThid )
        CALL READ_REC_XY_RS( PFCn_concFile,pfcnconc1,intime1,
     &       myIter,myThid )
        ENDIF
        IF ( PFCi_depoFile .NE. ' ' ) THEN
         CALL READ_REC_XY_RS( PFCi_depoFile,pfcidepo0,intime0,
     &       myIter,myThid )
        CALL READ_REC_XY_RS( PFCi_depoFile,pfcidepo1,intime1,
     &       myIter,myThid )
        ENDIF
        IF ( PFCi_riverFile .NE. ' ' ) THEN
         CALL READ_REC_XY_RS( PFCi_riverFile,pfciriver0,intime0,
     &       myIter,myThid )
        CALL READ_REC_XY_RS( PFCi_riverFile,pfciriver1,intime1,
     &       myIter,myThid )
        ENDIF
        IF ( PFCP_riverFile .NE. ' ' ) THEN
         CALL READ_REC_XY_RS( PFCP_riverFile,pfcpriver0,intime0,
     &       myIter,myThid )
        CALL READ_REC_XY_RS( PFCP_riverFile,pfcpriver1,intime1,
     &       myIter,myThid )
        ENDIF
     
                
         _EXCH_XY_RS(pfcnconc0, myThid )
         _EXCH_XY_RS(pfcnconc1, myThid )
         _EXCH_XY_RS(pfcidepo0, myThid )
         _EXCH_XY_RS(pfcidepo1, myThid )
         _EXCH_XY_RS(pfciriver0, myThid )
         _EXCH_XY_RS(pfciriver1, myThid )
         _EXCH_XY_RS(pfcpriver0, myThid )
         _EXCH_XY_RS(pfcpriver1, myThid ) 

C-    save newly loaded time-record
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           PFC_ldRec_emis(bi,bj) = intime1
         ENDDO
        ENDDO

C--   end if-block for loading new time-records
       ENDIF

C       
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           IF ( PFCn_concFile .NE. ' '  ) THEN
             PFCnConc(i,j,bi,bj) = bWght*pfcnconc0(i,j,bi,bj)
     &                          +aWght*pfcnconc1(i,j,bi,bj)
           ELSE
             PFCnConc(i,j,bi,bj) = 0. _d 0
           ENDIF
           IF ( PFCi_depoFile .NE. ' '  ) THEN
             PFCiDepo(i,j,bi,bj) = bWght*pfcidepo0(i,j,bi,bj)
     &                          +aWght*pfcidepo1(i,j,bi,bj)
           ELSE
             PFCiDepo(i,j,bi,bj) = 0. _d 0
           ENDIF
           IF ( PFCi_riverFile .NE. ' '  ) THEN
             PFCiRiver(i,j,bi,bj) = bWght*pfciriver0(i,j,bi,bj)
     &                          +aWght*pfciriver1(i,j,bi,bj)
           ELSE
             PFCiRiver(i,j,bi,bj) = 0. _d 0
           ENDIF
           IF ( PFCP_riverFile .NE. ' '  ) THEN
             PFCPRiver(i,j,bi,bj) = bWght*pfcpriver0(i,j,bi,bj)
     &                          +aWght*pfcpriver1(i,j,bi,bj)
           ELSE
             PFCPRiver(i,j,bi,bj) = 0. _d 0
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO
C


C endif for periodicEmis ! of IF ( PFC_emisCycle .GT. 0. _d 0 ) 
      ENDIF 

      RETURN
      END
!EOC
