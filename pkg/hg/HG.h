!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HG.h
!
! !DESCRIPTION: Global variables and common blocks for the Hg simulation.
!\\
!\\
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
!
       COMMON /HG_NEEDS/
     &              Hg0Conc, Hg2Depo, Hg2River, HgpRiver,
#ifdef HG2REFRACTORY
     &              Hg2rRiver,
#endif
#ifdef HGPREFRACTORY
     &              HgprRiver,
#endif      
     &              pisvel, fice, radsw,
     &              chl, npp, doc, poc,
     &              Fdoc, Fpoc, Fremin,
     &              o2, dfe
#ifdef FOODW_MODE     
     &              , NO3, phytoP, zooP
     &              , phygrow, zoograze, zoomort
#endif     
      _RL  Hg0Conc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  Hg2Depo(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  Hg2River(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  HgpRiver(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#ifdef HG2REFRACTORY
      _RL  Hg2rRiver(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#endif
#ifdef HGPREFRACTORY
      _RL  HgprRiver(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,nrmax)
#endif                  
      _RL  pisvel(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  fice  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  radsw (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  chl   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL  npp   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL  doc   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)   
      _RL  poc   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL  Fdoc  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL  Fpoc  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL  Fremin(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL  o2    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL  dfe   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#ifdef FOODW_MODE      
      _RL  NO3   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)      
      _RL  phytoP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,npmax)
      _RL  zooP  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,nzmax)
      _RL  phygrow(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,npmax)
      _RL  zoograze(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,
     &              nSx,nSy,npmax,nzmax)
      _RL  zoomort(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,nzmax)
#endif

C     time-record currently loaded (in temp arrays *[1])
      COMMON /HG_LOAD_I/ HG_ldRec_forcing, HG_ldRec_chem
      INTEGER HG_ldRec_forcing(nSx,nSy), HG_ldRec_chem(nSx,nSy)

      COMMON /HG_LOAD/
     &    wind0, wind1, ice0, ice1, radsw0, radsw1,
     &    hg0conc0, hg0conc1, hg2depo0, hg2depo1,
     &    hg2river0, hg2river1, hgpriver0, hgpriver1,
#ifdef HG2REFRACTORY
     &    hg2rriver0, hg2rriver1,
#endif
#ifdef HGPREFRACTORY
     &    hgprriver0, hgprriver1,
#endif     
     &    chl0, chl1, npp0, npp1, doc0, doc1,
     &    poc0, poc1, Fdoc0, Fdoc1, Fpoc0, Fpoc1,
     &    Fremin0, Fremin1,
     &    o20,o21,dfe0,dfe1
#ifdef FOODW_MODE     
     &    , NO30, NO31,phytoP0, phytoP1, zooP0, zoop1
     &    , phygrow0,phygrow1,zoograze0,zoograze1
     &    , zoomort0,zoomort1
#endif
      _RS wind0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS wind1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS ice0    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS ice1    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS radsw0   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS radsw1   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS hg0conc0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS hg0conc1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS hg2depo0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS hg2depo1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS hg2river0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS hg2river1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS hgpriver0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS hgpriver1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) 
#ifdef HG2REFRACTORY
      _RS hg2rriver0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS hg2rriver1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#endif
#ifdef HGPREFRACTORY
      _RS hgprriver0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,nrmax)
      _RS hgprriver1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,nrmax) 
#endif           
      _RS chl0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS chl1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS npp0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS npp1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS doc0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS doc1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS poc0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS poc1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS Fdoc0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS Fdoc1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS Fpoc0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS Fpoc1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS Fremin0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS Fremin1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS o20 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS o21 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS dfe0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS dfe1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#ifdef FOODW_MODE      
      _RS NO30    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS NO31    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)      
      _RS phytoP0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,npmax)
      _RS phytoP1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,npmax)
      _RS zooP0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,nzmax)
      _RS zooP1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,nzmax)
      _RS phygrow0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,npmax)
      _RS phygrow1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,npmax)
      _RS zoograze0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,
     &              nSx,nSy,npmax,nzmax)
      _RS zoograze1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,
     &              nSx,nSy,npmax,nzmax)
      _RS zoomort0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,nzmax)
      _RS zoomort1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,nzmax)
#endif
C--   COMMON /HG_FILENAMES/
C  HG_windFile    :: file name of wind speeds
C  HG_iceFile     :: file name of seaice fraction
C  HG_forcingPeriod :: periodic forcing parameter specific for hg (seconds)
C  HG_forcingCycle  :: periodic forcing parameter specific for hg (seconds)
C  HG0_concFile     :: file name of hg0 concentration in atmosphere
C  HG2_depoFile     :: file name of hg2 + hgp deposition flux from atmosphere
C  HG2_riverFile    :: file name of Hg2 riverine runoff
C  HGP_riverFile    :: file name of HgP riverine runoff
C  HG2R_riverFile    :: file name of Hg2R riverine runoff
C  HGPR_riverFile    :: file name of HgPR riverine runoff
C  radsw_file       :: file name of short-wave radiation
C  chl_file         :: file name of chl concentration
C  npp_file         :: file name of net primary production
C  doc_file         :: file name of dissolved organic carbon concentration
C  poc_file         :: file name of particualte organic carbon concentration
C  Fdoc_file        :: file name of sinking dissolved organic carbon
C  Fpoc_file        :: file name of sinking particualte organic carbon
C  Fremin_file      :: file name of particualte organic carbon remineralization rate
C  o2_file          :: file name of o2 concentration
C  dfe_file         :: file name of dissolved fe deposition flux
C  NO3_file         :: file name of NO3 concentration
C  phytoP_file      :: file name of phytoplankton concentration
C  zooP_file        :: file name of zooplankton concentration
C  phygrow_file     :: file name of phytoplankton growth
C  zoograze_file    :: file name of zooplankton graze
C  zoomort_file     :: file name of zooplankton mortality
      COMMON /HG_FILENAMES/
     &        HG_windFile, HG_iceFile, radsw_File,
     &        HG_forcingPeriod, HG_forcingCycle,
     &        HG_chemPeriod, HG_chemCycle,
     &        HG0_concFile, HG2_depoFile,
     &        HG2_riverFile, HGP_riverFile,
#ifdef HG2REFRACTORY
     &        HG2R_riverFile,
#endif
#ifdef HGPREFRACTORY
     &        HGPR_riverFile,
#endif       
     &        chl_file, npp_file,
     &        doc_file, poc_file,
     &        Fdoc_file, Fpoc_file, Fremin_file,
     &        o2_file, dfe_file 
#ifdef FOODW_MODE     
     &        , NO3_file,phytoP_file, zooP_file
     &        , phygrow_file,zoograze_file,zoomort_file
#endif
      CHARACTER*(MAX_LEN_FNAM) HG_windFile
      CHARACTER*(MAX_LEN_FNAM) HG_iceFile
      CHARACTER*(MAX_LEN_FNAM) HG0_concFile
      CHARACTER*(MAX_LEN_FNAM) HG2_depoFile
      CHARACTER*(MAX_LEN_FNAM) HG2_riverFile
      CHARACTER*(MAX_LEN_FNAM) HGP_riverFile  
#ifdef HG2REFRACTORY
      CHARACTER*(MAX_LEN_FNAM) HG2R_riverFile
#endif
#ifdef HGPREFRACTORY
      CHARACTER*(MAX_LEN_FNAM) HGPR_riverFile 
#endif             
      CHARACTER*(MAX_LEN_FNAM) radsw_File
      CHARACTER*(MAX_LEN_FNAM) chl_File
      CHARACTER*(MAX_LEN_FNAM) npp_File
      CHARACTER*(MAX_LEN_FNAM) doc_File
      CHARACTER*(MAX_LEN_FNAM) poc_File
      CHARACTER*(MAX_LEN_FNAM) Fdoc_File
      CHARACTER*(MAX_LEN_FNAM) Fpoc_File
      CHARACTER*(MAX_LEN_FNAM) Fremin_File
      CHARACTER*(MAX_LEN_FNAM) o2_File
      CHARACTER*(MAX_LEN_FNAM) dfe_File
#ifdef FOODW_MODE      
      CHARACTER*(MAX_LEN_FNAM) NO3_File
      CHARACTER*(MAX_LEN_FNAM) phytoP_File
      CHARACTER*(MAX_LEN_FNAM) zooP_File
      CHARACTER*(MAX_LEN_FNAM) phygrow_File
      CHARACTER*(MAX_LEN_FNAM) zoograze_File
      CHARACTER*(MAX_LEN_FNAM) zoomort_File
#endif
      _RL     HG_forcingPeriod   ! data related to air-sea exchange
      _RL     HG_forcingCycle
      _RL     HG_chemPeriod      ! data related to Hg chemistry
      _RL     HG_chemCycle


C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!EOC
