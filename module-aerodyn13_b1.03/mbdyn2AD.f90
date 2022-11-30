! the module store variables getting from mbdyn solver
MODULE mbdyn2aerodyn
	
	USE SharedTypes
     
    IMPLICIT NONE
	
	TYPE(MARKER),save,allocatable			:: Blades_components(:)
	TYPE(MARKER),save			:: Nacelle_components
	TYPE(MARKER),save			:: Hub_components
	TYPE(MARKER),save			:: RotorFurl_components
	TYPE(MARKER),save			:: Tower_components
	
	TYPE(AllAeroMarkers),save	:: BladeElem_components

END MODULE

! these subroutines are interface to setting variable into mbdyn2aerodyn module
!
SUBROUTINE set_nacelle_components(position, orientation, transVel, rotVel)
	
	USE SharedTypes
	USE mbdyn2aerodyn
	
	IMPLICIT NONE
	
	REAL(4),INTENT(IN)			::position(3)
	REAL(4),INTENT(IN)			::orientation(9)
	REAL(4),INTENT(IN)			::transVel(3)
	REAL(4),INTENT(IN)			::rotVel(3)
	INTEGER                     :: i, j
    

	do i = 1,3
		Nacelle_components%Position(i) 			= position(i)
		Nacelle_components%TranslationVel(i)	= transVel(i)
		Nacelle_components%RotationVel(i)		= rotVel(i)
		do j = 1,3
			Nacelle_components%Orientation(j,i) = orientation(j+(i-1)*3)
		end do
    end do
	
END SUBROUTINE set_nacelle_components

SUBROUTINE set_hub_components(position, orientation, transVel, rotVel)
	
	USE SharedTypes
	USE mbdyn2aerodyn
	
	IMPLICIT NONE
	
	REAL(4),INTENT(IN)			::position(3)
	REAL(4),INTENT(IN)			::orientation(9)
	REAL(4),INTENT(IN)			::transVel(3)
	REAL(4),INTENT(IN)			::rotVel(3)
	INTEGER                     :: i, j


	do i = 1,3
		Hub_components%Position(i) 			= position(i)
		Hub_components%TranslationVel(i)	= transVel(i)
		Hub_components%RotationVel(i)		= rotVel(i)
		
		do j = 1,3
			Hub_components%Orientation(j,i) = orientation(j+(i-1)*3)
		end do
	end do
            
END SUBROUTINE set_hub_components

SUBROUTINE set_rotorfurl_components(position, orientation, transVel, rotVel)
	
	USE SharedTypes
	USE mbdyn2aerodyn
	
	IMPLICIT NONE
	
	REAL(4),INTENT(IN)			::position(3)
	REAL(4),INTENT(IN)			::orientation(9)
	REAL(4),INTENT(IN)			::transVel(3)
	REAL(4),INTENT(IN)			::rotVel(3)
    INTEGER                     :: i, j
	
	do i = 1,3
		RotorFurl_components%Position(i) 		= position(i)
		RotorFurl_components%TranslationVel(i)	= transVel(i)
		RotorFurl_components%RotationVel(i)		= rotVel(i)
		
		do j = 1,3
			RotorFurl_components%Orientation(j,i) = orientation(j+(i-1)*3)
		end do
	end do
	
END SUBROUTINE set_rotorfurl_components

SUBROUTINE set_tower_components(position, orientation, transVel, rotVel)
	
	USE SharedTypes
	USE mbdyn2aerodyn
	
	IMPLICIT NONE
	
	REAL(4),INTENT(IN)			::position(3)
	REAL(4),INTENT(IN)			::orientation(9)
	REAL(4),INTENT(IN)			::transVel(3)
	REAL(4),INTENT(IN)			::rotVel(3)
    INTEGER                     :: i, j
	
	do i = 1,3
		Tower_components%Position(i) 		= position(i)
		Tower_components%TranslationVel(i)	= transVel(i)
		Tower_components%RotationVel(i)		= rotVel(i)
		
		do j = 1,3
			Tower_components%Orientation(j,i) = orientation(j+(i-1)*3)
		end do
	end do
	
END SUBROUTINE set_tower_components

SUBROUTINE set_blades_components(c_blade, position, orientation, transVel, rotVel)
	
	USE SharedTypes
	USE mbdyn2aerodyn
	
    IMPLICIT NONE
    
	INTEGER(4)							::c_blade
	REAL(4),INTENT(IN)				::position(3)
	REAL(4),INTENT(IN)				::orientation(9)
	REAL(4),INTENT(IN)				::transVel(3)
	REAL(4),INTENT(IN)				::rotVel(3)
    INTEGER                         :: i, j, k
    
    
	if (c_blade == 1) then
		do j= 1,3
			Blades_components(1)%Position(j) 		= position(j)
			Blades_components(1)%TranslationVel(j)	= transVel(j)
			Blades_components(1)%RotationVel(j)		= rotVel(j)
                do k = 1,3
					Blades_components(1)%Orientation(k ,j) = orientation(k+(j-1)*3)
				end do

		end do
	else if (c_blade == 2) then
		do j= 1,3
			Blades_components(2)%Position(j) 		= position(j)
			Blades_components(2)%TranslationVel(j)	= transVel(j)
			Blades_components(2)%RotationVel(j)		= rotVel(j)
                do k = 1,3
					Blades_components(2)%Orientation(k ,j) = orientation(k+(j-1)*3)
				end do

		end do
	else if (c_blade == 3) then
		do j= 1,3
			Blades_components(3)%Position(j) 		= position(j)
			Blades_components(3)%TranslationVel(j)	= transVel(j)
			Blades_components(3)%RotationVel(j)		= rotVel(j)
		
				do k = 1,3
                    Blades_components(3)%Orientation(k,j) = orientation(k+(j-1)*3)
				end do
        end do
    end if

    
END SUBROUTINE set_blades_components

	
SUBROUTINE set_bladeelem_component(c_blade, c_elem, BladElem_position, BladElem_orientation, BladElem_transVel , BladElem_rotateVel)
    
    USE SharedTypes
	USE mbdyn2aerodyn
     
    IMPLICIT NONE
    
    ! passed variables
    INTEGER(4)              :: c_blade
    INTEGER(4)              :: c_elem
    REAL(4)                 :: BladElem_position(3)          
    REAL(4)                 :: BladElem_orientation(9)     
    REAL(4)                 :: BladElem_transVel(3)
    REAL(4)                 :: BladElem_rotateVel(3)
    INTEGER                 :: i, j
    
    
    
	do i = 1,3
		BladeElem_components%Blade(c_elem,c_blade)%Position(i) = BladElem_position(i)
		BladeElem_components%Blade(c_elem,c_blade)%TranslationVel(i) = BladElem_transVel(i)
		BladeElem_components%Blade(c_elem,c_blade)%RotationVel(i) = BladElem_rotatevel(i)
		
		do j = 1,3
			BladeElem_components%Blade(c_elem,c_blade)%Orientation(j,i) = BladElem_orientation(j+(i-1)*3)
		end do
	end do
   
END SUBROUTINE set_bladeelem_component

SUBROUTINE mbdyn_init( Version, NBlades, nelems )

!USE Identify
    USE Blade
    USE Wind  ! Use Global variable VX,VY,VZ.
    USE Element
    USE Precision
    USE mbdyn2aerodyn

    IMPLICIT NONE

    CHARACTER(26) Version
    INTEGER(4) NBlades
    INTEGER(4) nelems
    REAL(ReKi) RotRadius

        !DynProg = 'MBDyn '
        !DynVer = Version
        !Prog = ' '

        !CALL SetProgName

        NB = NBlades

        !B = NB    ! Need this value when we swith on "WAKE"
                  ! or "SWIRL" options
        !R = RotRadius
		
		NELM = nelems
		
        ALLOCATE(Blades_components(NBlades))
        
        RETURN

    END SUBROUTINE mbdyn_init

SUBROUTINE mbdyn_ad_inputgate(Length ,AD_input_file, FileNameLen, Out_root_name,ElemFileNameLen, *)
    
    USE AeroDyn
	USE mbdyn2aerodyn
    USE Blade
    USE Element
    IMPLICIT NONE
    
    ! passed Variables
    REAL(4),INTENT(IN)                     :: Length                  
    INTEGER(4) FileNameLen,ElemFileNameLen

    CHARACTER(FileNameLen),INTENT(IN)                :: AD_input_file             ! Name of the AeroDyn input file
    CHARACTER(ElemFileNameLen),INTENT(IN)                :: Out_root_name             ! Root name of the AeroDyn summary and element files
    
    
    ! FUNCTION AD_Initに渡す変数
    TYPE(AD_InitOptions)        :: AD_inputdata
    TYPE(AeroConfig)            :: Initial_Components
    INTEGER(4)                  :: ErrStat = 0
    

    ! Function definition
    TYPE(AllAeroMarkers)             :: AD_initial
	REAL(4) :: HubRadius
	
	! AllAeroMarkersの配列大きさ割り当て
	ALLOCATE(BladeElem_components%Blade(NELM,NB))
   
    !AD_Initの引数ADOptionsのセット
    AD_inputdata%ADInputFile = AD_input_file
    AD_inputdata%OutRootName = Out_root_name
    AD_inputdata%WrSumFile = .true.

    Initial_Components%Blade        = Blades_components
    Initial_Components%BladeLength  = Length
    Initial_Components%Hub          = HUB_components
    Initial_Components%RotorFurl    = RotorFurl_components
    Initial_Components%Nacelle      = Nacelle_components
    Initial_Components%Tower        = Tower_components
      
    
    AD_initial = AD_Init(AD_inputdata, Initial_Components, ErrStat)
    !ADOptions               ! Options for AeroDyn
    !TurbineComponents       ! initial configuration of the turbine components
    !ErrStat                 ! Determines if an error was encountered
	return 0
	
END SUBROUTINE mbdyn_ad_inputgate
	
	
SUBROUTINE call_ad_calculateloads(Curr_time,  Length,  DFN, DFT, PMA)
    
    USE AeroDyn
	USE mbdyn2aerodyn
    USE Blade
    USE Element
    
    IMPLICIT NONE
    
    INTEGER                 :: i, j
    
    ! passed Variables
    REAL(4),INTENT(IN)                     :: Curr_time
    REAL(4),INTENT(IN)                     :: Length                   ! 翼長さ
    
    TYPE(AeroConfig)            :: Turbine_components
    TYPE(AeroLoadsOptions)     :: Curr_ADoption
    INTEGER(4)                  :: ErrStat = 0
    TYPE(AllAeroMarkers)		:: All_aero_markers
    
    ! Function definition
    TYPE( AllAeroLoads )                :: AD_CalcuLoads           
    REAL(4)        :: DFN(NB*NELM)
    REAL(4)        :: DFT(NB*NELM)
    REAL(4)        :: PMA(NB*NELM)
    
	 ALLOCATE(Curr_ADoption%SetMulTabLoc(NELM,NB))
    
	All_aero_markers = BladeElem_components
    
    Turbine_components%Blade     	= Blades_components
    Turbine_components%BladeLength  = Length
    Turbine_components%Hub          = Hub_components
    Turbine_components%RotorFurl    = RotorFurl_components
    Turbine_components%Nacelle      = Nacelle_components
    Turbine_components%Tower        = Tower_components
    
    !AD_Calculateloadsの引数CurrentADOptionsのセット
    Curr_ADoption%LinearizeFlag     =.true.
    
    do  i = 1,NB
        do j = 1,NELM
            Curr_ADoption%SetMulTabLoc(j,i) = .false.
        end do
    end do
    
    
    
    AD_CalcuLoads = AD_CalculateLoads( Curr_time, All_aero_markers, Turbine_components, Curr_ADoption, ErrStat )

    do i = 1,NB
        do j = 1,NELM
            DFN(j+(i-1)*NELM) = AD_CalcuLoads%Blade(j,i)%Force(1)*DR(j)
            DFT(j+(i-1)*NELM) = AD_CalcuLoads%Blade(j,i)%Force(2)*DR(j)
            PMA(j+(i-1)*NELM) = AD_CalcuLoads%Blade(j,i)%Moment(3)*DR(j)
        end do
    end do
            
END SUBROUTINE call_ad_calculateloads
    

SUBROUTINE get_horizonal_wind_speed(ZTime, HorWindV)

    USE AeroDyn
    USE NWTC_Library
	USE mbdyn2aerodyn

    IMPLICIT NONE

    REAL(4)     :: ZTime
    REAL(4)     :: HorWindV
    REAL(4)     :: HHWndVec  (3)     ! Hub-height wind vector in the AeroDyn coordinate system.
    INTEGER(4)  :: ErrStat = 0
    INTEGER(4)  :: Sttus
    REAL(4)     :: HubHeight

    HubHeight = AD_GetConstant('RefHt', ErrStat )  !Get the hub height and use it instead of getting values directly from AeroDyn in future.  This SHOULD be changed to the actual structural hub height.
    HHWndVec(:) = AD_GetUndisturbedWind( ZTime, (/ REAL(0.0, ReKi), REAL(0.0, ReKi), HubHeight /), Sttus )  !bjj use turbine hub height!

    HorWindV = SQRT(HHWndVec(1)*HHWndVec(1) + HHWndVec(2)*HHWndVec(2))


END SUBROUTINE get_horizonal_wind_speed


SUBROUTINE MBDyn_true( val )

IMPLICIT NONE

LOGICAL val

        val = .TRUE.

        RETURN

END SUBROUTINE MBDyn_true

SUBROUTINE MBDyn_false( val )

IMPLICIT NONE

LOGICAL val

        val = .FALSE.

        RETURN

END SUBROUTINE MBDyn_false

! This subroutine is to pass the current simulation time 
! of MBDyn to AeroDyn!
! c_time: current time

SUBROUTINE MBDyn_sim_time(c_time)

USE Aerodyn
USE Blade
USE Element
USE Precision
USE AeroTime

IMPLICIT NONE

REAL(DbKi) c_time

        TIME = c_time
        
        RETURN
END SUBROUTINE MBDyn_sim_time

! This subroutine is to pass the current simulation time step 
! of MBDyn to AeroDyn!
! dt: time step

SUBROUTINE MBDyn_time_step(time_step)

USE Aerodyn
USE Blade
USE Element
USE Precision
USE AeroTime

IMPLICIT NONE

REAL(ReKi) time_step

        DT = time_step
        
        RETURN
END SUBROUTINE MBDyn_time_step

!=========================================================================================================================================
!======================================== control subroutines ============================================================================
!=========================================================================================================================================
MODULE VALtoDISCON

REAL(4),SAVE        :: Time
INTEGER,SAVE        :: Blnum

REAL(4),SAVE        :: BlPitch(3)
REAL(4),SAVE        :: GenSpeed
REAL(4),SAVE        :: HorWindV

END MODULE VALtoDISCON

!========================================================================================================
MODULE BladedDLLParameters


   ! This MODULE stores various PARAMETER constants used in the
   !   Bladed-style master controller DLL interface.  Users will need
   !   to set these values as required by their models:


!USE                             Precision


INTEGER(4), PARAMETER        :: N                 = 0                           ! No. of points in torque-speed look-up table: 0 = none and use the optimal mode PARAMETERs instead, nonzero = ignore the optimal mode PARAMETERs by setting Record 16 to 0.0 (-)
INTEGER(4), PARAMETER        :: Ptch_Cntrl        = 0                           ! Pitch control: 0 = collective, 1 = individual (-)

REAL(4), PARAMETER           :: DTCntrl           = 0.0001                      ! Communication interval for controller (sec) (this works in conjunction with FAST's time step (DT) the same way DTAero in AeroDyn works with FAST--see the FAST User's Guide description of DTAero for more information)
REAL(4), PARAMETER           :: Gain_OM           = 0.0                         ! Optimal mode gain (Nm/(rad/s)^2)
REAL(4), PARAMETER           :: GenPwr_Dem        = 0.0                         ! Demanded power (W)
REAL(4), PARAMETER           :: GenSpd_Dem        = 0.0                         ! Demanded generator speed above rated (rad/s)
REAL(4), PARAMETER           :: GenSpd_MaxOM      = 0.0                         ! Optimal mode maximum speed (rad/s)
REAL(4), PARAMETER           :: GenSpd_MinOM      = 0.0                         ! Minimum generator speed (rad/s)
REAL(4), PARAMETER           :: GenSpd_TLU(N)     = 0.0                         ! Table (array) containing N generator speeds  for the torque-speed table look-up (TLU) (rad/s) -- this should be defined using an array constructor; for example, if N = 3: GenSpd_TLU(N)    = (/ 0.0, 99.9, 999.9 /)
REAL(4), PARAMETER           :: GenTrq_Dem        = 0.0                         ! Demanded generator torque (Nm)
REAL(4), PARAMETER           :: GenTrq_TLU(N)     = 0.0                         ! Table (array) containing N generator torques for the torque-speed table look-up (TLU) (Nm   ) -- this should be defined using an array constructor, for example, if N = 3: GenTrq_TLU(N)    = (/ 0.0, 10.0, 200.0 /)
REAL(4), PARAMETER           :: Ptch_Max          = 0.0                         ! Maximum pitch angle (rad)
REAL(4), PARAMETER           :: Ptch_Min          = 0.0                         ! Minimum pitch angle (rad)
REAL(4), PARAMETER           :: Ptch_SetPnt       = 0.0                         ! Below-rated pitch angle set-point (rad)
REAL(4), PARAMETER           :: PtchRate_Max      = 0.0                         ! Maximum pitch rate                               (rad/s)
REAL(4), PARAMETER           :: PtchRate_Min      = 0.0                         ! Minimum pitch rate (most negative value allowed) (rad/s)
REAL(4), PARAMETER           :: NacYaw_North      = 0.0                         ! Reference yaw angle of the nacelle when the upwind end points due North (rad)

CHARACTER(1024), PARAMETER   :: DLL_FileName      = 'DISCON.dll'                ! The name of the DLL file including the full path to the current working directory.
CHARACTER(1024), PARAMETER   :: DLL_ProcName      = 'DISCON'                    ! The name of the procedure in the DLL that will be called.
CHARACTER(1024), PARAMETER   :: DLL_InFile        = 'DISCON.IN'

END MODULE BladedDLLParameters


!==========================================================================================================================================
! this subroutine is used to calculate the generator torque for Variable Speed Controll by calling the subroutine DISCON
!==========================================================================================================================================
SUBROUTINE call_controller(HSS_Spd, BlPitch_in, NumBl, ZTime,  GenEff, GenTrq, ElecPwr, BlPitchCom)

    USE AeroDyn
    USE NWTC_Library
    USE VALtoDISCON

    IMPLICIT NONE

    ! Passed variables:
    INTEGER(4), INTENT(IN )      :: NumBl                                           ! Number of blades, (-).

    !REAL(4), INTENT(IN )     :: DT                                              ! Integration time step, sec.
    REAL(4), INTENT(OUT)      :: ElecPwr                                         ! Electrical power (account for losses), watts.
    !REAL(4), INTENT(IN )     :: GBRatio                                         ! Gearbox ratio, (-).
    REAL(4), INTENT(IN )      :: GenEff                                          ! Generator efficiency, (-).
    REAL(4), INTENT(OUT)      :: GenTrq                                          ! Electrical generator torque, N-m.
    REAL(4), INTENT(IN )      :: HSS_Spd                                         ! HSS speed, rad/s.
    REAL(4), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.
    REAL(4), INTENT(IN )      :: BlPitch_in   (NumBl) 

    !CHARACTER(1024), INTENT(IN ) :: DirRoot       

       ! Local variables:

    REAL(4), INTENT(OUT)      :: BlPitchCom(NumBl)                               ! Commanded blade pitch angles (demand pitch angles), rad. 
    REAL(4)                   :: HSSBrFrac                                       ! Fraction of full braking torque: 0 (off) <= HSSBrFrac <= 1 (full), (-). 
    REAL(4)                   :: YawRateCom                                      ! Commanded nacelle-yaw angular rate (demand yaw rate), rad/s. 
    INTEGER(4)                :: ErrStat = 0
    INTEGER(4)                :: Sttus

    REAL(4)                   :: HubHeight   
    REAL(4)                   :: HHWndVec  (3)                                   ! Hub-height wind vector in the AeroDyn coordinate system.

    ! setting the variable to module
    Time = ZTime
    Blnum = NumBl

    BlPitch = BlPitch_in
    GenSpeed = HSS_Spd

    ! Get horizonal wind speed at hub height from aerodyn interface
    HubHeight = AD_GetConstant('RefHt', ErrStat )  !Get the hub height and use it instead of getting values directly from AeroDyn in future.  This SHOULD be changed to the actual structural hub height.
    HHWndVec(:) = AD_GetUndisturbedWind( ZTime, (/ REAL(0.0, ReKi), REAL(0.0, ReKi), HubHeight /), Sttus )  !bjj use turbine hub height!

    HorWindV = SQRT(HHWndVec(1)*HHWndVec(1) + HHWndVec(2)*HHWndVec(2))

   

   ! Compute the electric generator power:
   ! The generator efficiency is either additive for motoring,
   !   or subtractive for generating power.

    ! call subroutine BladedDLLInterface(args) for getting Gentrq and BlPitchCom
    CALL BladedDLLInterface ( NumBl, BlPitchCom, HSSBrFrac, GenTrq, YawRateCom )

    IF ( GenTrq > 0.0 )  THEN
       ElecPwr = GenTrq*HSS_Spd*GenEff
    ELSE
       ElecPwr = GenTrq*HSS_Spd/GenEff
    ENDIF

END SUBROUTINE call_controller
!=========================================================================================================================================
SUBROUTINE BladedDLLInterface ( NumBl, BlPitchCom, HSSBrFrac, GenTrq, YawRateCom )
   ! This SUBROUTINE is used to call a master controller implemented as
   !   a dynamic-link-library (DLL) in the style of Garrad Hassan's
   !   Bladed wind turbine software package.

    USE                             BladedDLLParameters
    USE                             VALtoDISCON

    IMPLICIT                        NONE


       ! Passed Variables:

    INTEGER(4), INTENT(IN )      :: NumBl                                           ! Number of blades, (-).

    REAL(4), INTENT(OUT)      :: BlPitchCom(NumBl)                               ! Commanded blade pitch angles (demand pitch angles) (rad).
    REAL(4), INTENT(OUT)      :: GenTrq                                          ! Electrical generator torque (N-m).
    REAL(4), INTENT(OUT)      :: HSSBrFrac                                       ! Fraction of full braking torque: 0 (off) <= HSSBrFrac <= 1 (full), (-).
    REAL(4), INTENT(OUT)      :: YawRateCom                                      ! Commanded nacelle-yaw angular rate (demand yaw rate) (rad/s).

    !CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.


       ! Local Variables:

    INTEGER(4), PARAMETER        :: R                 = 85                          ! Start of below-rated torque-speed look-up table (record no.)

    REAL(4)                      :: avrSWAP   (R+(2*N)-1)                           ! The swap array, used to pass data to, and receive data from, the DLL controller.
    REAL(4), SAVE,ALLOCATABLE    :: BlPitchCom_SAVE(:)                              ! Commanded blade pitch angles (demand pitch angles) (rad).
    REAL(4), SAVE                :: GenTrq_SAVE       = 0.0                         ! Electrical generator torque (N-m).
    REAL(4), SAVE                :: LastTime          = 0.0                         ! The simulation time at the end of the last DLL call (sec).
    REAL(4), PARAMETER           :: OnePlusEps        = 1.0 + EPSILON(OnePlusEps)   ! The number slighty greater than unity in the precision of 4.
    REAL(4), SAVE                :: YawRateCom_SAVE   = 0.0                         ! Commanded nacelle-yaw angular rate (demand yaw rate) (rad/s).

    INTEGER(4)                   :: aviFAIL                                         ! A flag used to indicate the success of the DLL call as follows: 0 if the DLL call was successful, >0 if the DLL call was successful but cMessage should be issued as a warning messsage, <0 if the DLL call was unsuccessful or for any other reason the simulation is to be stopped at this point with cMessage as the error message
    INTEGER(4), SAVE             :: BrkState_SAVE     = 0                           ! Shaft brake status: 0 = off, 1 = on (full) (-).
    INTEGER(4), SAVE             :: GenState_SAVE     = 1                           ! Generator contactor: 0 = off, 1 = main (high speed) or variable speed generator (-).
    INTEGER(4)                   :: I                                               ! Generic index.
    INTEGER(4)                   :: K                                               ! Loops through blades.
    INTEGER(4)                   :: Sttus                                           ! Status returned by an attempted allocation.

    INTEGER(1)                   :: accINFILE (  256)                               ! The address of the first record of an array of 1-byte CHARACTERs giving the name of the parameter input file, 'DISCON.IN'.
    INTEGER(1)                   :: avcMSG    (  256)                               ! The address of the first record of an array of 1-byte CHARACTERS returning a message that will be displayed if aviFAIL <> 0.
    INTEGER(1)                   :: avcOUTNAME( 1024)                               ! The address of the first record of an array of 1-byte CHARACTERS giving the simulation run name without extension.
    INTEGER(1)                   :: iInFile   (  256)                               ! CHARACTER string cInFile  stored as a 1-byte array.
    INTEGER(1)                   :: iMessage  (  256)                               ! CHARACTER string cMessage stored as a 1-byte array.
    INTEGER(1)                   :: iOutName  ( 1024)                               ! CHARACTER string cOutName stored as a 1-byte array.

    CHARACTER( 256)              :: cInFile                                         ! CHARACTER string giving the name of the parameter input file, 'DISCON.IN'
    CHARACTER( 256)              :: cMessage                                        ! CHARACTER string giving a message that will be displayed by the calling program if aviFAIL <> 0.
    CHARACTER(1024)              :: cOutName                                        ! CHARACTER string giving the simulation run name without extension.

    LOGICAL,    SAVE             :: FirstPas          = .TRUE.                      ! When .TRUE., indicates we're on the first pass.

   ! Set EQUIVALENCE relationships between INTEGER(1) byte arrays and CHARACTER strings:

   ! ALLOCATE and initilize the BlPitchCom_SAVE first pass only:

    EQUIVALENCE (iInFile , cInFile )
    EQUIVALENCE (iMessage, cMessage)
    EQUIVALENCE (iOutName, cOutName)

    IF ( FirstPas .AND. ( .NOT. ALLOCATED(BlPitchCom_SAVE) ) )  THEN

       ALLOCATE ( BlPitchCom_SAVE(NumBl) , STAT=Sttus )

       BlPitchCom_SAVE = 0.0

       ! Don't set FirstPas = .FALSE. here yet, since FirstPas is also needed in
       !   the next block of code.

    ENDIF


   ! Since multiple control routines (one each for pitch, yaw, and torque
   !   control) are calling this routine, let's add a test condition to make
   !   sure we don't call the DLL master controller more than once per time
   !   step (all subsequent calls to this routine after the first one will
   !   use the SAVEd values from the first call of the time step):
   ! NOTE: AllOuts(Time) is scaled by OnePlusEps to ensure that the controller
   !       is called at every time step when DTCntrl = DT, even in the presence
   !       of numerical precision errors.

    IF ( ( Time  - LastTime ) >= DTCntrl )  THEN  ! Make sure time has incremented a controller time step.  

        avrSWAP     = 0.0               ! Initialize the avrSWAP by zero

        IF ( FirstPas )  THEN
        ! Set status flag:

            avrSWAP( 1) = 0.0                                  ! Status flag set as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation (-)

            FirstPas = .FALSE.                                 ! Don't enter here again!

        ELSE

            avrSWAP( 1) = 1.0                                  ! Status flag set as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation (-)

        ENDIF

        ! Setting the avrSWAP by module
        avrSWAP( 2) = Time                          ! Current time (sec)
        avrSWAP( 3) = Time - LastTime               ! Communication interval (sec)
        avrSWAP( 4) = BlPitch(1)                    ! Blade 1 pitch angle (rad)
        avrSWAP(20) = GenSpeed                      ! Measured generator speed (rad/s)
        avrSWAP(27) = HorWindV                      ! Hub wind speed (m/s)
        avrSWAP(33) = BlPitch(2)                    ! Blade 2 pitch angle (rad)
        avrSWAP(34) = BlPitch(3)                    ! Blade 3 pitch angle (rad)     
        avrSWAP(61) = NumBl                         ! No. of blades (-)

        avrSWAP(10) = 0.0                                     ! 0 = pitch position actuator, 1 = pitch rate actuator (-) -- must be 0 for FAST

        avrSWAP(49) = 256.0                                   ! Maximum no. of characters allowed in the returned message (-)
        avrSWAP(50) = REAL( LEN_TRIM(cInFile)  )              ! No. of characters in the "INFILE"  argument (-)
        avrSWAP(51) = REAL( LEN_TRIM(cOutName) )              ! No. of characters in the "OUTNAME" argument (-)  

        ! Create the input file and outname file arguments to the DLL (this requires
        !   the CHARACTER strings to be converted to byte arrays):

        cInFile = DLL_InFile
        !cOutName = TRIM( DirRoot )

        DO I = 1,MIN(  256, NINT( avrSWAP(50) ) )
           accINFILE (I) = iInFile (I)   ! Same as cInFile  by EQUIVALENCE
        ENDDO
        DO I = 1,MIN( 1024, NINT( avrSWAP(51) ) )
           avcOUTNAME(I) = iOutName(I)   ! Same as cOutName by EQUIVALENCE
        ENDDO       

        ! Call the Bladed-style DLL controller:

        CALL DISCON (  avrSWAP, aviFAIL, accINFILE, avcOUTNAME, avcMSG )

        ! Load control demands (commands) out of the avrSWAP array according to
        !   Appendix A of the Bladed User Manual:

        GenState_SAVE = NINT( avrSWAP(35) )                   ! Generator contactor (-)
        BrkState_SAVE = NINT( avrSWAP(36) )                   ! Shaft brake status (-)

        ! Record 41, demanded yaw actuator torque, is ignored since record 29 is set to 0 by FAST indicating yaw rate control

        IF (     Ptch_Cntrl == 1 )  THEN ! Individual pitch control

            DO K = 1,NumBl ! Loop through all blades avrSWAP(42), avrSWAP(43), and, if NumBl = 3, avrSWAP(44)
                BlPitchCom_SAVE(K) = avrSWAP( 41 + K )          ! Demanded individual pitch position of blade K (rad)
            ENDDO ! K - blades

        ELSEIF ( Ptch_Cntrl == 0 )  THEN ! Collective pitch control

            BlPitchCom_SAVE       = avrSWAP(45)                ! Demanded pitch angle (Collective pitch) (rad)
        ! Record 46, demanded pitch rate (Collective pitch), is ingored since record 10 is set to 0 by FAST indicating pitch position actuator

        ENDIF

        GenTrq_SAVE     = avrSWAP(47)                         ! Demanded generator torque (Nm)
        YawRateCom_SAVE = avrSWAP(48)                         ! Demanded nacelle yaw rate (rad/s)      

        LastTime = Time

ENDIF




    ! Return the SAVEd demand values of pitch, HSS brake fraction, generator
    !   torque, and yaw at every call:

    BlPitchCom = BlPitchCom_SAVE
    HSSBrFrac  = BrkState_SAVE
    GenTrq     = GenState_SAVE*GenTrq_SAVE ! Set GenTrq to 0.0 if GenState = 0 = off, or GenTrq_SAVE if GenState = 1 = main (high speed) or variable speed generator
    YawRateCom = YawRateCom_SAVE



RETURN
END SUBROUTINE BladedDLLInterface
