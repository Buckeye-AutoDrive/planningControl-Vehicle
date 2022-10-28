#include "ModelWithControllersOnly.h"

/* Storage class 'PageSwitching' */
ModelWithControllersOn_cal_type ModelWithControllersOn_cal_impl = {
  /* Mask Parameter: VehicleBody3DOFSingleTrack_Cs
   * Referenced by: '<S142>/Cs'
   */
  { 0.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3,
    0.32999999999999996, 0.36, 0.39, 0.42, 0.45, 0.48000000000000004, 0.51, 0.54,
    0.57000000000000006, 0.60000000000000009, 0.63, 0.66, 0.69000000000000006,
    0.72, 0.75, 0.78, 0.81, 0.84000000000000008, 0.87, 0.9 },

  /* Mask Parameter: VehicleBody3DOFSingleTrack_Cym
   * Referenced by: '<S142>/Cym'
   */
  { 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12,
    0.13, 0.14, 0.15, 0.15999999999999998, 0.16999999999999998, 0.18, 0.19,
    0.19999999999999998, 0.21, 0.21999999999999997, 0.22999999999999998, 0.24,
    0.25, 0.26, 0.27, 0.27999999999999997, 0.29, 0.3 },

  /* Mask Parameter: MappedSteering_StrgAngBpts
   * Referenced by:
   *   '<S219>/1-D Lookup Table'
   *   '<S219>/1-D Lookup Table1'
   */
  { -4.71238898038469, 4.71238898038469 },

  /* Mask Parameter: VehicleBody3DOFSingleTrack_beta
   * Referenced by:
   *   '<S142>/Cs'
   *   '<S142>/Cym'
   */
  { 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12,
    0.13, 0.14, 0.15, 0.15999999999999998, 0.16999999999999998, 0.18, 0.19,
    0.19999999999999998, 0.21, 0.21999999999999997, 0.22999999999999998, 0.24,
    0.25, 0.26, 0.27, 0.27999999999999997, 0.29, 0.3 },

  /* Mask Parameter: LookupGain_bpts
   * Referenced by: '<S218>/1-D Lookup Table'
   */
  { -5000.0, 5000.0 },

  /* Mask Parameter: LookupGain_tbl
   * Referenced by: '<S218>/1-D Lookup Table'
   */
  { 1.0, 1.0 },

  /* Expression: w
   * Referenced by: '<S126>/vehicle model'
   */
  { 1.4, 1.4 },

  /* Expression: [0.5 0.5 0.5]
   * Referenced by: '<Root>/Constant4'
   */
  { 0.5, 0.5, 0.5 },

  /* Expression: [0 0]
   * Referenced by: '<S142>/Crm'
   */
  { 0.0, 0.0 },

  /* Expression: [-1 1]
   * Referenced by: '<S142>/Crm'
   */
  { -1.0, 1.0 },

  /* Expression: [4.*ones(2,1); 0]
   * Referenced by: '<S142>/4'
   */
  { 4.0, 4.0, 0.0 },

  /* Expression: [0; 0; 1]
   * Referenced by: '<S142>/Constant4'
   */
  { 0.0, 0.0, 1.0 },

  /* Expression: WhlLftTbl
   * Referenced by: '<S219>/1-D Lookup Table1'
   */
  { 0.27435999882697704, -0.27435999882697704 },

  /* Expression: WhlRghtTbl
   * Referenced by: '<S219>/1-D Lookup Table'
   */
  { 0.27435999882697704, -0.27435999882697704 },

  /* Mask Parameter: VehicleBody3DOFSingleTrack_Af
   * Referenced by: '<S142>/.5.*A.*Pabs.//R.//T'
   */
  2.0,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_Cd
   * Referenced by: '<S142>/Constant'
   */
  0.3,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_Cl
   * Referenced by: '<S142>/Constant1'
   */
  0.1,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_Cpm
   * Referenced by: '<S142>/Constant2'
   */
  0.1,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_Cy_f
   * Referenced by: '<S141>/Cyf'
   */
  1.2E+6,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_Cy_r
   * Referenced by: '<S141>/Cyr'
   */
  1.1E+6,

  /* Mask Parameter: LateralControllerStanley_DelayG
   * Referenced by: '<S8>/Gain2'
   */
  0.2,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_Fzno
   * Referenced by: '<S126>/vehicle model'
   */
  5000.0,

  /* Mask Parameter: VehicleDynamics_InitSpeed
   * Referenced by: '<S126>/xdot_oConstant'
   */
  0.0,

  /* Mask Parameter: PIForward_InitialConditionForIn
   * Referenced by: '<S49>/Integrator'
   */
  0.0,

  /* Mask Parameter: PIReverse_InitialConditionForIn
   * Referenced by: '<S100>/Integrator'
   */
  0.0,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_Izz
   * Referenced by:
   *   '<S126>/vehicle model'
   *   '<S148>/Constant2'
   */
  2050.0,

  /* Mask Parameter: LongitudinalControllerStanley_K
   * Referenced by:
   *   '<S46>/Integral Gain'
   *   '<S97>/Integral Gain'
   */
  2.5,

  /* Mask Parameter: LongitudinalControllerStanley_g
   * Referenced by:
   *   '<S54>/Proportional Gain'
   *   '<S105>/Proportional Gain'
   */
  10.0,

  /* Mask Parameter: PIForward_LowerSaturationLimit
   * Referenced by:
   *   '<S56>/Saturation'
   *   '<S42>/DeadZone'
   */
  -40.0,

  /* Mask Parameter: PIReverse_LowerSaturationLimit
   * Referenced by:
   *   '<S107>/Saturation'
   *   '<S93>/DeadZone'
   */
  -40.0,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_NF
   * Referenced by: '<S126>/vehicle model'
   */
  2.0,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_NR
   * Referenced by: '<S126>/vehicle model'
   */
  2.0,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_Pabs
   * Referenced by: '<S142>/.5.*A.*Pabs.//R.//T'
   */
  101325.0,

  /* Mask Parameter: LateralControllerStanley_Positi
   * Referenced by: '<S1>/Kinematic'
   */
  2.5,

  /* Mask Parameter: LateralControllerStanley_Posi_i
   * Referenced by: '<S1>/Kinematic'
   */
  2.5,

  /* Mask Parameter: DragForce_R
   * Referenced by: '<S142>/.5.*A.*Pabs.//R.//T'
   */
  287.058,

  /* Mask Parameter: HardPointCoordinateTransformFro
   * Referenced by: '<S150>/R_T2'
   */
  0.0,

  /* Mask Parameter: HardPointCoordinateTransformRea
   * Referenced by: '<S152>/R_T2'
   */
  0.0,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_Tair
   * Referenced by: '<S126>/AirTempConstant'
   */
  273.0,

  /* Mask Parameter: PIForward_UpperSaturationLimit
   * Referenced by:
   *   '<S56>/Saturation'
   *   '<S42>/DeadZone'
   */
  40.0,

  /* Mask Parameter: PIReverse_UpperSaturationLimit
   * Referenced by:
   *   '<S107>/Saturation'
   *   '<S93>/DeadZone'
   */
  40.0,

  /* Mask Parameter: LateralControllerStanley_YawRat
   * Referenced by: '<S8>/Gain'
   */
  0.2,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_g
   * Referenced by:
   *   '<S126>/vehicle model'
   *   '<S148>/Constant'
   */
  9.81,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_h
   * Referenced by:
   *   '<S126>/vehicle model'
   *   '<S150>/R_T3'
   *   '<S152>/R_T3'
   *   '<S153>/Constant2'
   */
  0.35,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_latO
   * Referenced by: '<S151>/latOff'
   */
  0.0,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_long
   * Referenced by: '<S151>/longOff'
   */
  0.0,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_mu
   * Referenced by: '<S205>/Constant'
   */
  1.0,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_r_o
   * Referenced by: '<S126>/r_oConstant'
   */
  0.0,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_sigm
   * Referenced by: '<S211>/Constant1'
   */
  0.1,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_si_n
   * Referenced by: '<S211>/Constant2'
   */
  0.1,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_vert
   * Referenced by: '<S151>/vertOff'
   */
  0.0,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_xdot
   * Referenced by:
   *   '<S126>/vehicle model'
   *   '<S201>/Constant'
   *   '<S202>/Constant'
   *   '<S168>/Constant'
   *   '<S169>/Constant'
   */
  0.01,

  /* Mask Parameter: VehicleBody3DOFSingleTrack_ydot
   * Referenced by: '<S126>/ydot_oConstant'
   */
  0.0,

  /* Expression: MaxSteeringAngle
   * Referenced by: '<S1>/Kinematic'
   */
  80.0,

  /* Expression: Wheelbase
   * Referenced by: '<S1>/Kinematic'
   */
  3.0,

  /* Expression: 0
   * Referenced by: '<S40>/Constant1'
   */
  0.0,

  /* Computed Parameter: Integrator_gainval
   * Referenced by: '<S49>/Integrator'
   */
  0.05,

  /* Expression: 0
   * Referenced by: '<S40>/ZeroGain'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S91>/Constant1'
   */
  0.0,

  /* Computed Parameter: Integrator_gainval_e
   * Referenced by: '<S100>/Integrator'
   */
  0.05,

  /* Expression: 0
   * Referenced by: '<S91>/ZeroGain'
   */
  0.0,

  /* Expression: Fxtire_sat
   * Referenced by: '<S126>/vehicle model'
   */
  70000.0,

  /* Expression: Fytire_sat
   * Referenced by: '<S126>/vehicle model'
   */
  70000.0,

  /* Expression: d
   * Referenced by: '<S126>/vehicle model'
   */
  0.0,

  /* Expression: 0.001
   * Referenced by: '<S5>/Saturation'
   */
  0.001,

  /* Expression: 180/pi
   * Referenced by: '<S3>/Gain1'
   */
  57.295779513082323,

  /* Expression: 1
   * Referenced by: '<Root>/Constant1'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<Root>/Constant'
   */
  0.0,

  /* Expression: 180/pi
   * Referenced by: '<S9>/Gain'
   */
  57.295779513082323,

  /* Expression: VehicleMass/(2*TireStiffness*(1+Lf/Lr))
   * Referenced by: '<S8>/Gain1'
   */
  0.09265,

  /* Expression: 180/pi
   * Referenced by: '<Root>/Gain1'
   */
  57.295779513082323,

  /* Expression: 0
   * Referenced by: '<S8>/Unit Delay'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Unit Delay'
   */
  0.0,

  /* Expression: MaxSteeringAngle
   * Referenced by: '<S8>/Saturation'
   */
  80.0,

  /* Expression: -MaxSteeringAngle
   * Referenced by: '<S8>/Saturation'
   */
  -80.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S119>/Integrator'
   */
  0.0,

  /* Expression: 1.2*0.3*2/2
   * Referenced by: '<S119>/Gain1'
   */
  0.36,

  /* Expression: 0.7*9.81
   * Referenced by: '<S119>/Gain2'
   */
  6.867,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S4>/Constant'
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S118>/Integrator'
   */
  0.0,

  /* Expression: 1.2*0.3*2/2
   * Referenced by: '<S118>/Gain1'
   */
  0.36,

  /* Expression: 0.7*9.81
   * Referenced by: '<S118>/Gain2'
   */
  6.867,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S4>/Switch'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S4>/Switch1'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S5>/Constant9'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S149>/Constant2'
   */
  0.0,

  /* Expression: -180/pi
   * Referenced by: '<S120>/Gain'
   */
  -57.295779513082323,

  /* Expression: -1
   * Referenced by: '<S120>/Gain5'
   */
  -1.0,

  /* Expression: 0
   * Referenced by: '<S149>/Constant8'
   */
  0.0,

  /* Expression: -1
   * Referenced by: '<S120>/Gain6'
   */
  -1.0,

  /* Expression: 0
   * Referenced by: '<S149>/Constant7'
   */
  0.0,

  /* Expression: -180/pi
   * Referenced by: '<S120>/Gain1'
   */
  -57.295779513082323,

  /* Expression: -180/pi
   * Referenced by: '<S120>/Gain2'
   */
  -57.295779513082323,

  /* Expression: -1
   * Referenced by: '<S120>/Gain3'
   */
  -1.0,

  /* Expression: 0
   * Referenced by: '<S149>/Constant'
   */
  0.0,

  /* Expression: -1
   * Referenced by: '<S120>/Gain4'
   */
  -1.0,

  /* Expression: 0
   * Referenced by: '<S120>/Rate Transition'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S120>/Rate Transition1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S120>/Rate Transition2'
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by:
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S121>/Constant2'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S121>/First Order Hold'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S121>/First Order Hold1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S129>/Constant'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S149>/Constant5'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S149>/Constant4'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S149>/Constant6'
   */
  0.0,

  /* Expression: dh
   * Referenced by: '<S153>/Constant'
   */
  1.0,

  /* Expression: hh
   * Referenced by: '<S153>/Constant3'
   */
  0.2,

  /* Expression: 0
   * Referenced by: '<S153>/Constant4'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S149>/Constant1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S149>/Constant10'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S149>/Constant3'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S149>/Constant9'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S212>/lateral'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S213>/lateral'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S127>/second-order low-pass torque model'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S127>/second-order low-pass brake model'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S127>/Unit Delay'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S127>/second-order low-pass steering model'
   */
  0.0,

  /* Expression: pi/180
   * Referenced by: '<S127>/deg-to-rad1'
   */
  0.017453292519943295,

  /* Expression: 0.5
   * Referenced by: '<S127>/Gain2'
   */
  0.5,

  /* Expression: 1
   * Referenced by: '<S5>/Constant3'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S5>/Constant4'
   */
  1.0,

  /* Expression: 120
   * Referenced by: '<S5>/Constant'
   */
  120.0,

  /* Expression: 0
   * Referenced by: '<S5>/Constant1'
   */
  0.0,

  /* Expression: -1000
   * Referenced by: '<S5>/Constant2'
   */
  -1000.0,

  /* Expression: 1
   * Referenced by: '<S5>/Constant5'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<Root>/Constant2'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<Root>/Constant3'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S12>/Constant4'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S12>/Constant1'
   */
  0.0,

  /* Expression: -1
   * Referenced by: '<S12>/Constant3'
   */
  -1.0,

  /* Expression: 1
   * Referenced by: '<S12>/Constant2'
   */
  1.0,

  /* Computed Parameter: Merge_InitialOutput
   * Referenced by: '<S11>/Merge'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S10>/Constant'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S10>/Switch'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S10>/Switch1'
   */
  0.0,

  /* Expression: session.data.EgoCarId
   * Referenced by: '<S120>/Constant'
   */
  1.0,

  /* Computed Parameter: uDLookupTable_numYWorkElts
   * Referenced by: '<S218>/1-D Lookup Table'
   */
  { 1U, 0U },

  /* Computed Parameter: uDLookupTable_maxIndex
   * Referenced by: '<S218>/1-D Lookup Table'
   */
  1U,

  /* Computed Parameter: uDLookupTable_dimSizes
   * Referenced by: '<S218>/1-D Lookup Table'
   */
  1U
};

ModelWithControllersOn_cal_type *ModelWithControllersOnly_cal =
  &ModelWithControllersOn_cal_impl;
