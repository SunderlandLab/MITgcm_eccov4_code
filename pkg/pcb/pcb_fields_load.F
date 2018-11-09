!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pcb_fields_load.F
!
! !DESCRIPTION: Initializes data arrays for the PCB simulation.
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "PCB_OPTIONS.h"

      SUBROUTINE PCB_FIELDS_LOAD (
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
#include "PCB_SIZE.h"
#include "PCB.h"
!
! !INPUT PARAMETERS:
!
      _RL myTime
      INTEGER myIter
      INTEGER myThid
!
! !REMARKS:
!
! !REVISION HISTORY:
!  07 Jan 2015 - H. Amos     - Removed inactive CFC code to clean
!                              up the module a bit. 
!  03 Feb 2015 - H. Amos     - Update rivers based on Yanxu Zhang's Hg code.
!  29 Jun 2015 - H. Amos     - Move atmospheric concentration and
!                              deposition so they are controlled by the
!                              PCB_forcingCycle.
!  04 Feb 2016 - R. Yantosca - Added ProTex headers
!  04 Mar 2016 - R. Yantosca - Removed #ifdef FOODW_MODE
!  04 Mar 2016 - R. Yantosca - Removed commented-out code
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      _RL aWght, bWght
      INTEGER intimeP, intime0, intime1
      INTEGER bi, bj, i, j, k
      INTEGER prec

      ! Mimic Hg for PCB river inputs
      ! (H. Amos, 03 Feb 2015)
#ifdef PCBPREFRACTORY
      CHARACTER*(MAX_LEN_FNAM) pcbR_File_Act
      INTEGER nriver
      CHARACTER*(2) nrivers
#endif /* PCBPREFRACTORY */

      LOGICAL nowRead

      !--------------------------------------------------------------
      ! 2D forcing fields (e.g., atmospheric deposition)
      !--------------------------------------------------------------
      IF ( PCB_forcingCycle .GT. 0. _d 0 ) THEN

       ! First call requires that we initialize everything to zero 
       ! for safety
       IF ( myIter .EQ. nIter0 ) THEN
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
            PCB_ldRec_forcing(bi,bj) = 0
          ENDDO
         ENDDO
       ENDIF

       ! Now calculate whether it is time to update the forcing arrays
       CALL GET_PERIODIC_INTERVAL(
     O                   intimeP, intime0, intime1, bWght, aWght,
     I                   PCB_forcingCycle, PCB_forcingPeriod,
     I                   deltaTclock, myTime, myThid )

       bi = myBxLo(myThid)
       bj = myByLo(myThid)

#ifdef ALLOW_DEBUG
       IF ( debugLevel.GE.debLevB ) THEN
        _BEGIN_MASTER(myThid)
        WRITE(standardMessageUnit,'(A,I10,A,4I5,A,2F14.10)')
     &   ' PCB_FIELDS_LOAD,', myIter,
     &   ' : iP,iLd,i0,i1=', intimeP,PCB_ldRec_forcing(bi,bj), 
     &                          intime0,intime1,
     &   ' ; Wght=', bWght, aWght
        _END_MASTER(myThid)
       ENDIF
#endif /* ALLOW_DEBUG */

#ifdef ALLOW_AUTODIFF_TAMC
       ! Assuming that we call S/R PCB_FIELDS_LOAD at each time-step and
       ! with increasing time, this will catch when we need to load new records;
       ! But with Adjoint run, this is not always the case => might end-up using
       ! the wrong time-records
       IF ( intime0.NE.intimeP .OR. myIter.EQ.nIter0 ) THEN
#else /* ALLOW_AUTODIFF_TAMC */
       ! Make no assumption on sequence of calls to PCB_FIELDS_LOAD ;
       ! This is the correct formulation (works in Adjoint run).
       ! Unfortunatly, might produce many recomputations <== not used until it is fixed
       IF ( intime1.NE.PCB_ldRec_forcing(bi,bj) ) THEN         
#endif /* ALLOW_AUTODIFF_TAMC */

        !-----------------------------------------------------------
        ! If the above condition is met then we need to read in
        ! data for the period ahead and the period behind myTime.
        !-----------------------------------------------------------
        IF ( debugLevel.GE.debLevZero ) THEN
         _BEGIN_MASTER(myThid)
         WRITE(standardMessageUnit,'(A,I10,A,2(2I5,A))')
     &    ' PCB_FIELDS_LOAD, it=', myIter,
     &    ' : Reading new data, i0,i1=', intime0, intime1,
     &    ' (prev=', intimeP, PCB_ldRec_forcing(bi,bj), ' )'
         _END_MASTER(myThid)
        ENDIF

        _BARRIER

        ! Note: In ECCOv4, wind speed, solar radiation, and sea ice 
        !       are part of the online physical simulation. You do
        !       not need to provide an external forcing file.

        ! Wind speed (m/s)
        IF ( PCB_windFile .NE. ' '  .AND. .NOT.useEXF ) THEN
         CALL READ_REC_XY_RS( PCB_windFile,wind0,intime0,
     &        myIter,myThid )
         CALL READ_REC_XY_RS( PCB_windFile,wind1,intime1,
     &        myIter,myThid )
        ENDIF
        
        ! Solar radiation (units?)
        IF ( radsw_File .NE. ' ' ) THEN
         CALL READ_REC_XY_RS( radsw_File,radsw0,intime0,
     &        myIter,myThid )
         CALL READ_REC_XY_RS( radsw_File,radsw1,intime1,
     &        myIter,myThid )
        ENDIF

        ! Sea ice (units?)
        IF ( PCB_iceFile .NE. ' '  .AND. .NOT.useSEAICE ) THEN
         CALL READ_REC_XY_RS( PCB_iceFile,ice0,intime0,
     &       myIter,myThid )
         CALL READ_REC_XY_RS( PCB_iceFile,ice1,intime1,
     &       myIter,myThid )
        ENDIF

        ! atmospheric concentration (mol/m3)
        IF ( PCBa_concFile .NE. ' ' ) THEN
         CALL READ_REC_XY_RS( PCBa_concFile,pcbaconc0,intime0,
     &       myIter,myThid )
        CALL READ_REC_XY_RS( PCBa_concFile,pcbaconc1,intime1,
     &       myIter,myThid )
        ENDIF

        ! gas-phase PCB atmospheric depostiion (mol/m2/sec)
        IF ( PCBa_depoFile .NE. ' ' ) THEN
         CALL READ_REC_XY_RS( PCBa_depoFile,pcbadepo0,intime0,
     &       myIter,myThid )
         CALL READ_REC_XY_RS( PCBa_depoFile,pcbadepo1,intime1,
     &       myIter,myThid )
        ENDIF

        ! particulate PCB atmospheric deposition (mol/m2/sec)
        IF ( PCBp_depoFile .NE. ' ' ) THEN
         CALL READ_REC_XY_RS( PCBp_depoFile,pcbpdepo0,intime0,
     &       myIter,myThid )
         CALL READ_REC_XY_RS( PCBp_depoFile,pcbpdepo1,intime1,
     &       myIter,myThid )
 
        ENDIF

        !-----------------------------------------------------------
        ! Handle exchange of _RS 2D scalar arrays 
        !-----------------------------------------------------------
         _EXCH_XY_RS(wind0, myThid )     ! wind speed
         _EXCH_XY_RS(wind1, myThid )     ! " "
         _EXCH_XY_RS(radsw0, myThid )    ! solar radiation
         _EXCH_XY_RS(radsw1, myThid )    ! " "
        IF (.NOT.useSEAICE) THEN
         _EXCH_XY_RS(ice0, myThid )      ! sea ice
         _EXCH_XY_RS(ice1, myThid )      ! " "
        ENDIF
         _EXCH_XY_RS(pcbaconc0, myThid ) ! PCBa atmospheric concentration
         _EXCH_XY_RS(pcbaconc1, myThid ) ! " "
         _EXCH_XY_RS(pcbadepo0, myThid ) ! PCBa atmopsheric deposition
         _EXCH_XY_RS(pcbadepo1, myThid ) ! " "
         _EXCH_XY_RS(pcbpdepo0, myThid ) ! PCBp atmospheric deposition
         _EXCH_XY_RS(pcbpdepo1, myThid ) ! " "


        !-----------------------------------------------------------
        ! Save newly loaded time-record
        !-----------------------------------------------------------
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           PCB_ldRec_forcing(bi,bj) = intime1
         ENDDO
        ENDDO

       ! End if-block for loading new time-records
       ENDIF
       
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx

           ! Wind speed (m/s)
           IF ( PCB_windFile .NE. ' '  ) THEN
             wind(i,j,bi,bj) = bWght*wind0(i,j,bi,bj)
     &                           +aWght*wind1(i,j,bi,bj)
           ELSE
             wind = 5. _d 0*maskC(i,j,1,bi,bj)
           ENDIF

           ! Solar radiation (units?)
           IF ( radsw_File .NE. ' '  ) THEN
             radsw(i,j,bi,bj) = bWght*radsw0(i,j,bi,bj)
     &                          +aWght*radsw1(i,j,bi,bj)
           ELSE
             radsw(i,j,bi,bj) = 0. _d 0
           ENDIF

           ! sea ice (units?)
           IF ( PCB_iceFile .NE. ' '  ) THEN
            FIce(i,j,bi,bj)    = bWght*ice0(i,j,bi,bj)
     &                          +aWght*ice1(i,j,bi,bj)
           ELSE
            FIce(i,j,bi,bj) =0. _d 0
           ENDIF

           ! Gaseous atmospheric concentrations (mol/m3)
           IF ( PCBa_concFile .NE. ' '  ) THEN
             PCBaConc(i,j,bi,bj) = bWght*pcbaconc0(i,j,bi,bj)
     &                          +aWght*pcbaconc1(i,j,bi,bj)
           ELSE
             PCBaConc(i,j,bi,bj) = 0. _d 0
           ENDIF

           ! Gas-phase atmospheric deposition (mol/m2/sec)
           IF ( PCBa_depoFile .NE. ' '  ) THEN
             PCBaDEP(i,j,bi,bj) = bWght*pcbadepo0(i,j,bi,bj)
     &                          +aWght*pcbadepo1(i,j,bi,bj)
           ELSE
             PCBaDEP(i,j,bi,bj) = 0. _d 0
           ENDIF

           ! Particulate-phase atmospheric deposition (mol/m2/sec)
           IF ( PCBp_depoFile .NE. ' '  ) THEN
             PCBpDEP(i,j,bi,bj) = bWght*pcbpdepo0(i,j,bi,bj)
     &                          +aWght*pcbpdepo1(i,j,bi,bj)
           ELSE
             PCBpDEP(i,j,bi,bj) = 0. _d 0
           ENDIF

          ENDDO
         ENDDO
        ENDDO
       ENDDO


      ENDIF ! endif for periodicForcing

      !--------------------------------------------------------------
      ! Biogeochemistry & DARWIN
      !--------------------------------------------------------------
      IF ( PCB_chemCycle .GT. 0. _d 0 ) THEN

       ! First call requires that we initialize everything to zero 
       ! for safety
       IF ( myIter .EQ. nIter0 ) THEN
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
            PCB_ldRec_chem(bi,bj) = 0
          ENDDO
         ENDDO
       ENDIF

       ! Now calculate whether it is time to update the forcing arrays
       CALL GET_PERIODIC_INTERVAL(
     O                   intimeP, intime0, intime1, bWght, aWght,
     I                   PCB_chemCycle, PCB_chemPeriod,
     I                   deltaTclock, myTime, myThid )

       bi = myBxLo(myThid)
       bj = myByLo(myThid)

#ifdef ALLOW_DEBUG
       IF ( debugLevel.GE.debLevB ) THEN
        _BEGIN_MASTER(myThid)
        WRITE(standardMessageUnit,'(A,I10,A,4I5,A,2F14.10)')
     &   ' PCB_FIELDS_LOAD,', myIter,
     &   ' : iP,iLd,i0,i1=', intimeP,PCB_ldRec_chem(bi,bj), 
     &                       intime0,intime1,
     &   ' ; Wght=', bWght, aWght
        _END_MASTER(myThid)
       ENDIF
#endif /* ALLOW_DEBUG */

#ifdef ALLOW_AUTODIFF_TAMC
       ! Assuming that we call S/R PCB_FIELDS_LOAD at each time-step and
       ! with increasing time, this will catch when we need to load new records;
       ! But with Adjoint run, this is not always the case => might end-up using
       ! the wrong time-records
       IF ( intime0.NE.intimeP .OR. myIter.EQ.nIter0 ) THEN
#else /* ALLOW_AUTODIFF_TAMC */
       ! Make no assumption on sequence of calls to HG_FIELDS_LOAD ;
       ! This is the correct formulation (works in Adjoint run).
       ! Unfortunatly, might produce many recomputations <== not used until it is fixed
       IF ( intime1.NE.PCB_ldRec_chem(bi,bj) ) THEN
#endif /* ALLOW_AUTODIFF_TAMC */

        !-----------------------------------------------------------
        ! If the above condition is met then we need to read in
        ! data for the period ahead and the period behind myTime.
        !-----------------------------------------------------------
        IF ( debugLevel.GE.debLevZero ) THEN
         _BEGIN_MASTER(myThid)
         WRITE(standardMessageUnit,'(A,I10,A,2(2I5,A))')
     &    ' PCB_FIELDS_LOAD, it=', myIter,
     &    ' : Reading new data, i0,i1=', intime0, intime1,
     &    ' (prev=', intimeP, PCB_ldRec_chem(bi,bj), ' )'
         _END_MASTER(myThid)
        ENDIF

        _BARRIER


        ! Riverine inputs of dissolved PCB (kg/s)
        IF ( PCBa_riverFile .NE. ' ' ) THEN
         CALL READ_REC_XY_RS( PCBa_riverFile,pcbariver0,intime0,
     &       myIter,myThid )
        CALL READ_REC_XY_RS( PCBa_riverFile,pcbariver1,intime1,
     &       myIter,myThid )
        ENDIF

        ! Riverine inputs of particle-bound PCB (kg/s)
        IF ( PCBp_riverFile .NE. ' ' ) THEN
         CALL READ_REC_XY_RS( PCBp_riverFile,pcbpriver0,intime0,
     &       myIter,myThid )
        CALL READ_REC_XY_RS( PCBp_riverFile,pcbpriver1,intime1,
     &       myIter,myThid )
        ENDIF

        ! Refractory riverine inputs
        ! Reference: Zhang et al., 2015, Global Biogeochem. Cycles
#ifdef PCBAREFRACTORY        
        IF ( PCBaR_riverFile .NE. ' ' ) THEN
         CALL READ_REC_XY_RS( PCBaR_riverFile,pcbarriver0,intime0,
     &       myIter,myThid )
        CALL READ_REC_XY_RS( PCBaR_riverFile,pcbarriver1,intime1,
     &       myIter,myThid )
        ENDIF
#endif /* PCBAREFRACTORY */

#ifdef PCBPREFRACTORY        
        IF ( PCBpR_riverFile .NE. ' ' ) THEN
         DO nriver=1,nrmax
          IF (nrmax .EQ. 1) THEN
            PCBpR_file_Act = trim(PCBpR_riverFile)
     &                   // '.bin'
          ELSE
            WRITE (nrivers,'(I2.2)') nriver
            PCBpR_file_Act = trim(PCBpR_riverFile) // '.' 
     &                   // nrivers // '.bin'
          ENDIF
          CALL READ_REC_XY_RS( PCBpR_file_Act,
     &       pcbprriver0(1-OLx,1-OLy,1,1,nriver),intime0,
     &       myIter,myThid )
          CALL READ_REC_XY_RS( PCBpR_file_Act,
     &       pcbprriver1(1-OLx,1-OLy,1,1,nriver),intime1,
     &       myIter,myThid )
         ENDDO
        ENDIF
#endif /* PCBPREFRACTORY */              
      
        ! Precision
        prec=32

        ! Chlorophyll (mg chl/m3)
        IF ( chl_File .NE. ' ' ) THEN
         CALL READ_REC_3D_RS( chl_File,prec,Nr,chl0,intime0,
     &        myIter,myThid )
         CALL READ_REC_3D_RS( chl_File,prec,Nr,chl1,intime1,
     &        myIter,myThid )
        ENDIF

        ! Net primary productivity (mmol P/m3/s)
        IF ( npp_File .NE. ' ' ) THEN
         CALL READ_REC_3D_RS( npp_File,prec,Nr,npp0,intime0,
     &        myIter,myThid )
         CALL READ_REC_3D_RS( npp_File,prec,Nr,npp1,intime1,
     &        myIter,myThid )
        ENDIF

        ! Dissolved organic carbon (mmol C/m3)
        IF ( doc_file .NE. ' ' ) THEN
         CALL READ_REC_3D_RS( doc_file,prec,Nr,doc0,intime0,
     &        myIter,myThid )
         CALL READ_REC_3D_RS( doc_file,prec,Nr,doc1,intime1,
     &        myIter,myThid )
        ENDIF

        ! Particulate organic carbon (mmol C/m3)
        IF ( poc_file .NE. ' ' ) THEN
         CALL READ_REC_3D_RS( poc_file,prec,Nr,poc0,intime0,
     &        myIter,myThid )
         CALL READ_REC_3D_RS( poc_file,prec,Nr,poc1,intime1,
     &        myIter,myThid )
        ENDIF

        ! Sinking flux of POC (mmol C/m3 * m/d)
        IF ( Fpoc_file .NE. ' ' ) THEN
         CALL READ_REC_3D_RS( Fpoc_file,prec,Nr,Fpoc0,intime0,
     &        myIter,myThid )
         CALL READ_REC_3D_RS( Fpoc_file,prec,Nr,Fpoc1,intime1,
     &        myIter,myThid )
        ENDIF

        ! Remineralization flux of POC and DOC (mmol C/m3/d)
        IF ( Fremin_file .NE. ' ' ) THEN
         CALL READ_REC_3D_RS( Fremin_file,prec,Nr,Fremin0,intime0,
     &        myIter,myThid )
         CALL READ_REC_3D_RS( Fremin_file,prec,Nr,Fremin1,intime1,
     &        myIter,myThid )
        ENDIF

         !-------------------------------------------------------
         ! Handle exchange of _RS 2D scalar arrays              
         !-------------------------------------------------------
         _EXCH_XY_RS(pcbariver0, myThid )  ! PCBa river inputs
         _EXCH_XY_RS(pcbariver1, myThid )  ! " "
         _EXCH_XY_RS(pcbpriver0, myThid )  ! PCBp river inputs
         _EXCH_XY_RS(pcbpriver1, myThid )  ! " "
#ifdef PCBAREFRACTORY
          ! Refractory PCBa river inputs
          _EXCH_XY_RS(pcbarriver0(1-OLx,1-OLy,1,1),myThid)
          _EXCH_XY_RS(pcbarriver1(1-OLx,1-OLy,1,1),myThid)
#endif /* PCBAREFRACTORY */
#ifdef PCBPREFRACTORY
         DO nriver=1,nrmax
          ! Refractory PCBp river inputs
          _EXCH_XY_RS(pcbprriver0(1-OLx,1-OLy,1,1,nriver),myThid)
          _EXCH_XY_RS(pcbprriver1(1-OLx,1-OLy,1,1,nriver),myThid)
         ENDDO 
#endif /* PCBPREFRACTORY */                                 
         _EXCH_XYZ_RS(chl0, myThid )       ! Chlorophyll
         _EXCH_XYZ_RS(chl1, myThid )       ! " "
         _EXCH_XYZ_RS(npp0, myThid )       ! Net primary productivity
         _EXCH_XYZ_RS(npp1, myThid )       ! " "
         _EXCH_XYZ_RS(doc0, myThid )       ! Dissolved organic carbon
         _EXCH_XYZ_RS(doc1, myThid )       ! " "
         _EXCH_XYZ_RS(poc0, myThid )       ! Particulate organic carbon
         _EXCH_XYZ_RS(poc1, myThid )       ! " "
         _EXCH_XYZ_RS(Fpoc0, myThid )      ! Sinking flux of POC
         _EXCH_XYZ_RS(Fpoc1, myThid )      ! " "
         _EXCH_XYZ_RS(Fremin0, myThid )    ! Remineralization of DOC and POC
         _EXCH_XYZ_RS(Fremin1, myThid )    ! " "

        ! Save newly loaded time-record
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           PCB_ldRec_chem(bi,bj) = intime1
         ENDDO
        ENDDO

       ENDIF ! ENDIF for loading new time records

       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx

           ! Dissolved riverine inputs (kg/s)
           IF ( PCBa_riverFile .NE. ' '  ) THEN
             PCBaRiver(i,j,bi,bj) = bWght*pcbariver0(i,j,bi,bj)
     &                          +aWght*pcbariver1(i,j,bi,bj)
           ELSE
             PCBaRiver(i,j,bi,bj) = 0. _d 0
           ENDIF

           ! Particulate riverine inputs (kg/s)
           IF ( PCBP_riverFile .NE. ' '  ) THEN
             PCBPRiver(i,j,bi,bj) = bWght*pcbpriver0(i,j,bi,bj)
     &                          +aWght*pcbpriver1(i,j,bi,bj)
           ELSE
             PCBPRiver(i,j,bi,bj) = 0. _d 0
           ENDIF

#ifdef PCBAREFRACTORY
            ! Refractory dissolved riverine inputs (kg/s)
            IF ( PCBaR_riverFile .NE. ' '  ) THEN
              PCBaRRiver(i,j,bi,bj) =
     &                           bWght*pcbarriver0(i,j,bi,bj)
     &                          +aWght*pcbarriver1(i,j,bi,bj)
            ELSE
              PCBarRiver(i,j,bi,bj) = 0. _d 0
            ENDIF
#endif /* PCBAREFRACTORY */

#ifdef PCBPREFRACTORY
           ! Refractory particulate riverine inputs (kg/s)
           DO nriver=1,nrmax
            IF ( PCBPR_riverFile .NE. ' '  ) THEN
              PCBPRRiver(i,j,bi,bj,nriver) =
     &                           bWght*pcbaprriver0(i,j,bi,bj,nriver)
     &                          +aWght*pcbaprriver1(i,j,bi,bj,nriver)
            ELSE
              PCBprRiver(i,j,bi,bj,nriver) = 0. _d 0
            ENDIF 
           ENDDO
#endif /* PCBPREFRACTORY */           

          ENDDO
         ENDDO
        ENDDO
       ENDDO


      ! DARWIN & FOODWEB parameters
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx

           ! Chlorophyll (mg chl/m3)
           IF ( chl_File .NE. ' '  ) THEN
             chl(i,j,k,bi,bj) = bWght*chl0(i,j,k,bi,bj)
     &                          +aWght*chl1(i,j,k,bi,bj)
           ELSE
             chl(i,j,k,bi,bj) = 0. _d 0
           ENDIF

           ! Net primary productivity (mmol P/m3/s)
           IF ( npp_File .NE. ' '  ) THEN
             npp(i,j,k,bi,bj) = bWght*npp0(i,j,k,bi,bj)
     &                          +aWght*npp1(i,j,k,bi,bj)
           ELSE
             npp(i,j,k,bi,bj) = 0. _d 0
           ENDIF

           ! Dissolved organic carbon (mmol C/m3)
           IF ( doc_File .NE. ' '  ) THEN
             doc(i,j,k,bi,bj) = bWght*doc0(i,j,k,bi,bj)
     &                          +aWght*doc1(i,j,k,bi,bj)
           ELSE
             doc(i,j,k,bi,bj) = 0. _d 0
           ENDIF

           ! Particulate organic carbon (mmol C/m3)
           IF ( poc_File .NE. ' '  ) THEN
             poc(i,j,k,bi,bj) = bWght*poc0(i,j,k,bi,bj)
     &                          +aWght*poc1(i,j,k,bi,bj)
           ELSE
             poc(i,j,k,bi,bj) = 0. _d 0
           ENDIF

           ! Sinking flux of POC (mmol C/m3 * m/d)
           IF ( Fpoc_File .NE. ' '  ) THEN
             Fpoc(i,j,k,bi,bj) = bWght*Fpoc0(i,j,k,bi,bj)
     &                          +aWght*Fpoc1(i,j,k,bi,bj)
           ELSE
             Fpoc(i,j,k,bi,bj) = 0. _d 0
           ENDIF
   
           ! Remineralization flux of POC and DOC (mmol C/m3/d)
           IF ( Fremin_File .NE. ' '  ) THEN
             Fremin(i,j,k,bi,bj) = bWght*Fremin0(i,j,k,bi,bj)
     &                          +aWght*Fremin1(i,j,k,bi,bj)
           ELSE
             Fremin(i,j,k,bi,bj) = 0. _d 0
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      ENDIF ! ENDIF for periodicChem

      !-------------------------------------------------------------
      ! Emissions
      ! - Note: Xianming Zhang added this for PFOS, but it's not
      !         being used by PCBs or Hg. (hma, 29 Jun 2015)
      !-------------------------------------------------------------
      IF ( PCB_emisCycle .GT. 0. _d 0 ) THEN

       ! First call requires that we initialize everything to zero 
       ! for safety
       IF ( myIter .EQ. nIter0 ) THEN
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
            PCB_ldRec_emis(bi,bj) = 0
          ENDDO
         ENDDO
       ENDIF

       ! Now calculate whether it is time to update the forcing arrays
       CALL GET_PERIODIC_INTERVAL(
     O                   intimeP, intime0, intime1, bWght, aWght,
     I                   PCB_emisCycle, PCB_emisPeriod,
     I                   deltaTclock, myTime, myThid )

       bi = myBxLo(myThid)
       bj = myByLo(myThid)
#ifdef ALLOW_DEBUG
       IF ( debugLevel.GE.debLevB ) THEN
        _BEGIN_MASTER(myThid)
        WRITE(standardMessageUnit,'(A,I10,A,4I5,A,2F14.10)')
     &   ' PCB_FIELDS_LOAD,', myIter,
     &   ' : iP,iLd,i0,i1=', intimeP,PCB_ldRec_emis(bi,bj), 
     &                       intime0,intime1,
     &   ' ; Wght=', bWght, aWght
        _END_MASTER(myThid)
       ENDIF
#endif /* ALLOW_DEBUG */

#ifdef ALLOW_AUTODIFF_TAMC
       ! Assuming that we call S/R HG_FIELDS_LOAD at each time-step and
       ! with increasing time, this will catch when we need to load new records;
       ! But with Adjoint run, this is not always the case => might end-up using
       ! the wrong time-records
       IF ( intime0.NE.intimeP .OR. myIter.EQ.nIter0 ) THEN
#else /* ALLOW_AUTODIFF_TAMC */
       ! Make no assumption on sequence of calls to HG_FIELDS_LOAD ;
       ! This is the correct formulation (works in Adjoint run).
       ! Unfortunatly, might produce many recomputations <== not used until it is fixed
       IF ( intime1.NE.PCB_ldRec_emis(bi,bj) ) THEN
#endif /* ALLOW_AUTODIFF_TAMC */

        !-----------------------------------------------------------
        ! If the above condition is met then we need to read in
        ! data for the period ahead and the period behind myTime.
        !-----------------------------------------------------------
        IF ( debugLevel.GE.debLevZero ) THEN
         _BEGIN_MASTER(myThid)
         WRITE(standardMessageUnit,'(A,I10,A,2(2I5,A))')
     &    ' PCB_FIELDS_LOAD, it=', myIter,
     &    ' : Reading new data, i0,i1=', intime0, intime1,
     &    ' (prev=', intimeP, PCB_ldRec_chem(bi,bj), ' )'
         _END_MASTER(myThid)
        ENDIF

        _BARRIER

        !-----------------------------------------------------------
        ! Save newly loaded time-record
        !-----------------------------------------------------------
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           PCB_ldRec_emis(bi,bj) = intime1
         ENDDO
        ENDDO

       ENDIF ! endif for loading new time records

      ENDIF ! endif for periodicEmis

      RETURN
      END
!EOC
