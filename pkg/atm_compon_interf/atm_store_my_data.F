C $Header: /u/gcmpack/MITgcm/pkg/atm_compon_interf/atm_store_my_data.F,v 1.7 2015/12/25 04:37:19 jmc Exp $
C $Name:  $

#include "ATM_CPL_OPTIONS.h"

CBOP
C     !ROUTINE: ATM_STORE_MY_DATA
C     !INTERFACE:
      SUBROUTINE ATM_STORE_MY_DATA(
     I                              myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE ATM_STORE_MY_DATA
C     | o Routine for controlling storage of coupling data to
C     |   coupler layer.
C     *==========================================================*
C     | This version interfaces to the MITgcm AIMPHYS package.
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "CPL_PARAMS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myTime    :: Current internal time.
C     myIter    :: Current timestep number.
C     myThid    :: my Thread Id number
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef COMPONENT_MODULE
C     !LOCAL VARIABLES:
C     bi, bj    :: Tile indices
      INTEGER bi,bj

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

C     Store coupling data ready for export
        CALL ATM_STORE_DYNVARS(     bi, bj, myTime, myIter, myThid )
        CALL ATM_STORE_SURFFLUX(    bi, bj, myTime, myIter, myThid )
        IF ( useAIM ) THEN
         CALL ATM_STORE_AIM_WNDSTR( bi, bj, myTime, myIter, myThid )
        ENDIF
        IF ( useAIM .AND. atm_cplExch_DIC ) THEN
         CALL ATM_STORE_AIM_FIELDS( bi, bj, myTime, myIter, myThid )
        ENDIF
        IF ( useLand ) THEN
         CALL ATM_STORE_LAND(       bi, bj, myTime, myIter, myThid )
        ENDIF
C-    This call has been moved to CPL_EXPORT_IMPORT_DATA so that
C     thSIce state-vars are stored just before beeing exported to CPL
c       IF ( useThSIce ) THEN
c        CALL ATM_STORE_THSICE(     bi, bj, myTime, myIter, myThid )
c       ENDIF

C--   end bi,bj loops.
       ENDDO
      ENDDO

#endif /* COMPONENT_MODULE */

      RETURN
      END
