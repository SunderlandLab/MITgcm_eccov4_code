!------------------------------------------------------------------------------
!              Harvard Biogeochemistry of Global Pollutants Group             !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pfc_partition.F
!
! !DESCRIPTION: Partitions the PFC species between each other.
!\\
!\\
! !INTERFACE:
!
#include "GCHEM_OPTIONS.h"
#include "PFC_OPTIONS.h"

      SUBROUTINE PFC_PARTITION( Ptr,
     &                         bi,bj,imin,imax,jmin,jmax,
     &                         myIter,myTime,myThid)
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
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PFC_SIZE.h"
#include "PFC.h"
!
! !INPUT PARAMETERS: 
!
      INTEGER myIter
      _RL myTime
      INTEGER myThid
      INTEGER i,j,k,bi, bj, imin, imax, jmin, jmax
!
! !INPUT/OUTPUT PARAMETERS: 
!
      _RL  Ptr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
!
! !REMARKS:
!  based on the code of Hg species in the ocean by yanxu
!  To calculate the partition of PFC species in the ocean
!  Local equlibrium between the anionic, neutral, and particulate species are 
!  assumed to be reached instantaneously
!
!
! !REVISION HISTORY:
!  29 Feb 2016 - R. Yantosca - Added ProTex headers
!  01 Mar 2016 - R. Yantosca - Fixed program flow, updated comments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#ifdef ALLOW_DIAGNOSTICS      
      _RL  pfcipart_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  pfcnpart_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  pfcppart_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)            
      CHARACTER*8 diagname
#endif /* ALLOW_DIAGNOSTICS */
      _RL  TOT_PFC
      _RL FRAC_AnionDissolved, FRAC_NeutralDissolved, FRAC_TotDissolved,
      _RL FRAC_POCBound
      _RL poco, doco
      _RL  pfci_old,pfcn_old,pfcp_old,pfci_new,pfcn_new,pfcp_new
      _RL  FRAC_Dissolved

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PFC

#ifdef ALLOW_DIAGNOSTICS
      !----------------------------------------------------------------
      ! Initialize diagnostics
      !----------------------------------------------------------------
      DO k=1,Nr
        DO j=jMin,jMax
          DO i=iMin,iMax
            pfcipart_diag(i,j,k)=0d0
            pfcnpart_diag(i,j,k)=0d0  
            pfcppart_diag(i,j,k)=0d0         
          ENDDO
        ENDDO
      ENDDO
#endif  /* ALLOW_DIAGNOSTICS */

      !----------------------------------------------------------------
      ! Start calculation
      !----------------------------------------------------------------
      DO j=jMin,jMax
       DO i=iMin,iMax
        DO k=1,Nr

         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN

           !===========================================================
           ! %%%%% Instantaneous equilibrium %%%%%
           !===========================================================
          
           ! POC and DOC concentration (units?)
           poco = max(poc(i,j,k,bi,bj),0d0)
           doco = max(doc(i,j,k,bi,bj),0d0)
 
           !-----------------------------------------------------------
           ! PFC concentrations before partitioning 
           !-----------------------------------------------------------
           pfci_old = max(Ptr(i,j,k,bi,bj,ipfci),0d0)
           pfcn_old = max(Ptr(i,j,k,bi,bj,ipfcn),0d0)
           pfcp_old = max(Ptr(i,j,k,bi,bj,ipfcp),0d0)
           TOT_PFC  = pfci_old + pfcn_old + pfcp_old

           IF (TOT_PFC .gt. 0. _d 0) THEN

             !---------------------------------------------------------
             ! Calculate fraction in the dissolved phase
             !---------------------------------------------------------
             FRAC_Dissolved = 1. _d 0 
     &             /(1. _d 0 + Koc_pfc * poco 
     &             * 1.2 _d - 8 )  ! mmolC/m3 -> kgC/L

!              write(*,*) 
!     &          'xz check FRAC_Dissolved in PFC_Partition.f',
!     &           FRAC_Dissolved 

             !---------------------------------------------------------
             ! pH of seawater
             !
             ! Note: Currently not considering the real pH
             !       (xz)
             !---------------------------------------------------------
             pH(i,j,k,bi,bj) = 8. _d 0 
           
             !---------------------------------------------------------
             ! Particulate fraction (i.e., PFC bound to POC)
             !---------------------------------------------------------
             FRAC_NeutralDissolved =  FRAC_Dissolved 
     &            *1.d0/(1.d0+10.d0**(pH(i,j,k,bi,bj)-pKa))
             FRAC_AnionDissolved =  FRAC_Dissolved 
     &            *10.d0**(pH(i,j,k,bi,bj)-pKa)
     &                  /(1.d0+10.d0**(pH(i,j, k, bi,bj)-pKa))
             FRAC_POCBound = 1.d0 - FRAC_Dissolved
        
!              write(*,*)'xz check PFRC_AnionDissolved ', 
!     &            ' FRAC_NeutralDissolved ', 
!     &            ' FRAC_POCBound in PFC_Partition.f',
!     &            FRAC_AnionDissolved, 
!     &            FRAC_NeutralDissolved, 
!     &            FRAC_POCBound

             !---------------------------------------------------------
             ! New PFC concentrations (units?) in the dissolved
             ! and particle phases after consdering partitioning
             !---------------------------------------------------------
             pfci_new = TOT_PFC * FRAC_AnionDissolved
             pfcn_new = TOT_PFC * FRAC_NeutralDissolved
             pfcp_new = TOT_PFC * FRAC_POCBound

             !---------------------------------------------------------
             ! Update tracer array
             !---------------------------------------------------------
             Ptr(i,j,k,bi,bj,ipfci) = pfci_new
             Ptr(i,j,k,bi,bj,ipfcn) = pfcn_new
             Ptr(i,j,k,bi,bj,ipfcp) = pfcp_new

#ifdef ALLOW_DIAGNOSTICS
             pfcipart_diag(i,j,k) = (pfci_old - pfci_new) / deltaTclock
             pfcnpart_diag(i,j,k) = (pfcn_old - pfcn_new) / deltaTclock
             pfcppart_diag(i,j,k) = (pfcp_old - pfcp_new) / deltaTclock
#endif /* ALLOW_DIAGNOSTICS */

           ENDIF        ! TOT_PFC > 0
         ENDIF          ! (HFacC(i,j,k,bi,bj).gt.0. _d 0) 
        ENDDO           ! k
       ENDDO            ! i
      ENDDO             ! j

c write diagnostics
#ifdef ALLOW_DIAGNOSTICS   
      WRITE(diagname,'(A8)') 'PFCnPART'
      CALL DIAGNOSTICS_FILL( pfcnpart_diag(1-Olx,1-Oly,1), diagname,
     &                       0,Nr,2,bi,bj,myThid )
     
      WRITE(diagname,'(A8)') 'PFCiPART'
      CALL DIAGNOSTICS_FILL( pfcipart_diag(1-Olx,1-Oly,1), diagname,
     &                       0,Nr,2,bi,bj,myThid )

      WRITE(diagname,'(A8)') 'PFCpPART'
      CALL DIAGNOSTICS_FILL( pfcppart_diag(1-Olx,1-Oly,1), diagname,
     &                       0,Nr,2,bi,bj,myThid )   
#endif /* ALLOW_DIAGNOSTICS */
 
#endif /* ALLOW_PFC         */
#endif /* ALLOW_PGRACERS    */
C     write(*,*) 'xz checking call PFC_PARTITION end'

      RETURN
      END
!EOC
