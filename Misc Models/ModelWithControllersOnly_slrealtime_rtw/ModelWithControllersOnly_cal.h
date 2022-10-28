#ifndef RTW_HEADER_ModelWithControllersOnly_cal_h_
#define RTW_HEADER_ModelWithControllersOnly_cal_h_
#include "rtwtypes.h"
#include "ModelWithControllersOnly_types.h"

/* Storage class 'PageSwitching', for system '<Root>' */
struct ModelWithControllersOn_cal_type {
  real_T VehicleBody3DOFSingleTrack_Cs[31];
                                /* Mask Parameter: VehicleBody3DOFSingleTrack_Cs
                                 * Referenced by: '<S142>/Cs'
                                 */
  real_T VehicleBody3DOFSingleTrack_Cym[31];
                               /* Mask Parameter: VehicleBody3DOFSingleTrack_Cym
                                * Referenced by: '<S142>/Cym'
                                */
  real_T MappedSteering_StrgAngBpts[2];
                                   /* Mask Parameter: MappedSteering_StrgAngBpts
                                    * Referenced by:
                                    *   '<S219>/1-D Lookup Table'
                                    *   '<S219>/1-D Lookup Table1'
                                    */
  real_T VehicleBody3DOFSingleTrack_beta[31];
                              /* Mask Parameter: VehicleBody3DOFSingleTrack_beta
                               * Referenced by:
                               *   '<S142>/Cs'
                               *   '<S142>/Cym'
                               */
  real_T LookupGain_bpts[2];           /* Mask Parameter: LookupGain_bpts
                                        * Referenced by: '<S218>/1-D Lookup Table'
                                        */
  real_T LookupGain_tbl[2];            /* Mask Parameter: LookupGain_tbl
                                        * Referenced by: '<S218>/1-D Lookup Table'
                                        */
  real_T vehiclemodel_w[2];            /* Expression: w
                                        * Referenced by: '<S126>/vehicle model'
                                        */
  real_T Constant4_Value[3];           /* Expression: [0.5 0.5 0.5]
                                        * Referenced by: '<Root>/Constant4'
                                        */
  real_T Crm_tableData[2];             /* Expression: [0 0]
                                        * Referenced by: '<S142>/Crm'
                                        */
  real_T Crm_bp01Data[2];              /* Expression: [-1 1]
                                        * Referenced by: '<S142>/Crm'
                                        */
  real_T u_Gain[3];                    /* Expression: [4.*ones(2,1); 0]
                                        * Referenced by: '<S142>/4'
                                        */
  real_T Constant4_Value_a[3];         /* Expression: [0; 0; 1]
                                        * Referenced by: '<S142>/Constant4'
                                        */
  real_T uDLookupTable1_tableData[2];  /* Expression: WhlLftTbl
                                        * Referenced by: '<S219>/1-D Lookup Table1'
                                        */
  real_T uDLookupTable_tableData[2];   /* Expression: WhlRghtTbl
                                        * Referenced by: '<S219>/1-D Lookup Table'
                                        */
  real_T VehicleBody3DOFSingleTrack_Af;
                                /* Mask Parameter: VehicleBody3DOFSingleTrack_Af
                                 * Referenced by: '<S142>/.5.*A.*Pabs.//R.//T'
                                 */
  real_T VehicleBody3DOFSingleTrack_Cd;
                                /* Mask Parameter: VehicleBody3DOFSingleTrack_Cd
                                 * Referenced by: '<S142>/Constant'
                                 */
  real_T VehicleBody3DOFSingleTrack_Cl;
                                /* Mask Parameter: VehicleBody3DOFSingleTrack_Cl
                                 * Referenced by: '<S142>/Constant1'
                                 */
  real_T VehicleBody3DOFSingleTrack_Cpm;
                               /* Mask Parameter: VehicleBody3DOFSingleTrack_Cpm
                                * Referenced by: '<S142>/Constant2'
                                */
  real_T VehicleBody3DOFSingleTrack_Cy_f;
                              /* Mask Parameter: VehicleBody3DOFSingleTrack_Cy_f
                               * Referenced by: '<S141>/Cyf'
                               */
  real_T VehicleBody3DOFSingleTrack_Cy_r;
                              /* Mask Parameter: VehicleBody3DOFSingleTrack_Cy_r
                               * Referenced by: '<S141>/Cyr'
                               */
  real_T LateralControllerStanley_DelayG;
                              /* Mask Parameter: LateralControllerStanley_DelayG
                               * Referenced by: '<S8>/Gain2'
                               */
  real_T VehicleBody3DOFSingleTrack_Fzno;
                              /* Mask Parameter: VehicleBody3DOFSingleTrack_Fzno
                               * Referenced by: '<S126>/vehicle model'
                               */
  real_T VehicleDynamics_InitSpeed; /* Mask Parameter: VehicleDynamics_InitSpeed
                                     * Referenced by: '<S126>/xdot_oConstant'
                                     */
  real_T PIForward_InitialConditionForIn;
                              /* Mask Parameter: PIForward_InitialConditionForIn
                               * Referenced by: '<S49>/Integrator'
                               */
  real_T PIReverse_InitialConditionForIn;
                              /* Mask Parameter: PIReverse_InitialConditionForIn
                               * Referenced by: '<S100>/Integrator'
                               */
  real_T VehicleBody3DOFSingleTrack_Izz;
                               /* Mask Parameter: VehicleBody3DOFSingleTrack_Izz
                                * Referenced by:
                                *   '<S126>/vehicle model'
                                *   '<S148>/Constant2'
                                */
  real_T LongitudinalControllerStanley_K;
                              /* Mask Parameter: LongitudinalControllerStanley_K
                               * Referenced by:
                               *   '<S46>/Integral Gain'
                               *   '<S97>/Integral Gain'
                               */
  real_T LongitudinalControllerStanley_g;
                              /* Mask Parameter: LongitudinalControllerStanley_g
                               * Referenced by:
                               *   '<S54>/Proportional Gain'
                               *   '<S105>/Proportional Gain'
                               */
  real_T PIForward_LowerSaturationLimit;
                               /* Mask Parameter: PIForward_LowerSaturationLimit
                                * Referenced by:
                                *   '<S56>/Saturation'
                                *   '<S42>/DeadZone'
                                */
  real_T PIReverse_LowerSaturationLimit;
                               /* Mask Parameter: PIReverse_LowerSaturationLimit
                                * Referenced by:
                                *   '<S107>/Saturation'
                                *   '<S93>/DeadZone'
                                */
  real_T VehicleBody3DOFSingleTrack_NF;
                                /* Mask Parameter: VehicleBody3DOFSingleTrack_NF
                                 * Referenced by: '<S126>/vehicle model'
                                 */
  real_T VehicleBody3DOFSingleTrack_NR;
                                /* Mask Parameter: VehicleBody3DOFSingleTrack_NR
                                 * Referenced by: '<S126>/vehicle model'
                                 */
  real_T VehicleBody3DOFSingleTrack_Pabs;
                              /* Mask Parameter: VehicleBody3DOFSingleTrack_Pabs
                               * Referenced by: '<S142>/.5.*A.*Pabs.//R.//T'
                               */
  real_T LateralControllerStanley_Positi;
                              /* Mask Parameter: LateralControllerStanley_Positi
                               * Referenced by: '<S1>/Kinematic'
                               */
  real_T LateralControllerStanley_Posi_i;
                              /* Mask Parameter: LateralControllerStanley_Posi_i
                               * Referenced by: '<S1>/Kinematic'
                               */
  real_T DragForce_R;                  /* Mask Parameter: DragForce_R
                                        * Referenced by: '<S142>/.5.*A.*Pabs.//R.//T'
                                        */
  real_T HardPointCoordinateTransformFro;
                              /* Mask Parameter: HardPointCoordinateTransformFro
                               * Referenced by: '<S150>/R_T2'
                               */
  real_T HardPointCoordinateTransformRea;
                              /* Mask Parameter: HardPointCoordinateTransformRea
                               * Referenced by: '<S152>/R_T2'
                               */
  real_T VehicleBody3DOFSingleTrack_Tair;
                              /* Mask Parameter: VehicleBody3DOFSingleTrack_Tair
                               * Referenced by: '<S126>/AirTempConstant'
                               */
  real_T PIForward_UpperSaturationLimit;
                               /* Mask Parameter: PIForward_UpperSaturationLimit
                                * Referenced by:
                                *   '<S56>/Saturation'
                                *   '<S42>/DeadZone'
                                */
  real_T PIReverse_UpperSaturationLimit;
                               /* Mask Parameter: PIReverse_UpperSaturationLimit
                                * Referenced by:
                                *   '<S107>/Saturation'
                                *   '<S93>/DeadZone'
                                */
  real_T LateralControllerStanley_YawRat;
                              /* Mask Parameter: LateralControllerStanley_YawRat
                               * Referenced by: '<S8>/Gain'
                               */
  real_T VehicleBody3DOFSingleTrack_g;
                                 /* Mask Parameter: VehicleBody3DOFSingleTrack_g
                                  * Referenced by:
                                  *   '<S126>/vehicle model'
                                  *   '<S148>/Constant'
                                  */
  real_T VehicleBody3DOFSingleTrack_h;
                                 /* Mask Parameter: VehicleBody3DOFSingleTrack_h
                                  * Referenced by:
                                  *   '<S126>/vehicle model'
                                  *   '<S150>/R_T3'
                                  *   '<S152>/R_T3'
                                  *   '<S153>/Constant2'
                                  */
  real_T VehicleBody3DOFSingleTrack_latO;
                              /* Mask Parameter: VehicleBody3DOFSingleTrack_latO
                               * Referenced by: '<S151>/latOff'
                               */
  real_T VehicleBody3DOFSingleTrack_long;
                              /* Mask Parameter: VehicleBody3DOFSingleTrack_long
                               * Referenced by: '<S151>/longOff'
                               */
  real_T VehicleBody3DOFSingleTrack_mu;
                                /* Mask Parameter: VehicleBody3DOFSingleTrack_mu
                                 * Referenced by: '<S205>/Constant'
                                 */
  real_T VehicleBody3DOFSingleTrack_r_o;
                               /* Mask Parameter: VehicleBody3DOFSingleTrack_r_o
                                * Referenced by: '<S126>/r_oConstant'
                                */
  real_T VehicleBody3DOFSingleTrack_sigm;
                              /* Mask Parameter: VehicleBody3DOFSingleTrack_sigm
                               * Referenced by: '<S211>/Constant1'
                               */
  real_T VehicleBody3DOFSingleTrack_si_n;
                              /* Mask Parameter: VehicleBody3DOFSingleTrack_si_n
                               * Referenced by: '<S211>/Constant2'
                               */
  real_T VehicleBody3DOFSingleTrack_vert;
                              /* Mask Parameter: VehicleBody3DOFSingleTrack_vert
                               * Referenced by: '<S151>/vertOff'
                               */
  real_T VehicleBody3DOFSingleTrack_xdot;
                              /* Mask Parameter: VehicleBody3DOFSingleTrack_xdot
                               * Referenced by:
                               *   '<S126>/vehicle model'
                               *   '<S201>/Constant'
                               *   '<S202>/Constant'
                               *   '<S168>/Constant'
                               *   '<S169>/Constant'
                               */
  real_T VehicleBody3DOFSingleTrack_ydot;
                              /* Mask Parameter: VehicleBody3DOFSingleTrack_ydot
                               * Referenced by: '<S126>/ydot_oConstant'
                               */
  real_T Kinematic_MaxSteeringAngle;   /* Expression: MaxSteeringAngle
                                        * Referenced by: '<S1>/Kinematic'
                                        */
  real_T Kinematic_Wheelbase;          /* Expression: Wheelbase
                                        * Referenced by: '<S1>/Kinematic'
                                        */
  real_T Constant1_Value;              /* Expression: 0
                                        * Referenced by: '<S40>/Constant1'
                                        */
  real_T Integrator_gainval;           /* Computed Parameter: Integrator_gainval
                                        * Referenced by: '<S49>/Integrator'
                                        */
  real_T ZeroGain_Gain;                /* Expression: 0
                                        * Referenced by: '<S40>/ZeroGain'
                                        */
  real_T Constant1_Value_f;            /* Expression: 0
                                        * Referenced by: '<S91>/Constant1'
                                        */
  real_T Integrator_gainval_e;       /* Computed Parameter: Integrator_gainval_e
                                      * Referenced by: '<S100>/Integrator'
                                      */
  real_T ZeroGain_Gain_i;              /* Expression: 0
                                        * Referenced by: '<S91>/ZeroGain'
                                        */
  real_T vehiclemodel_Fxtire_sat;      /* Expression: Fxtire_sat
                                        * Referenced by: '<S126>/vehicle model'
                                        */
  real_T vehiclemodel_Fytire_sat;      /* Expression: Fytire_sat
                                        * Referenced by: '<S126>/vehicle model'
                                        */
  real_T vehiclemodel_d;               /* Expression: d
                                        * Referenced by: '<S126>/vehicle model'
                                        */
  real_T Saturation_LowerSat;          /* Expression: 0.001
                                        * Referenced by: '<S5>/Saturation'
                                        */
  real_T Gain1_Gain;                   /* Expression: 180/pi
                                        * Referenced by: '<S3>/Gain1'
                                        */
  real_T Constant1_Value_k;            /* Expression: 1
                                        * Referenced by: '<Root>/Constant1'
                                        */
  real_T Constant_Value;               /* Expression: 0
                                        * Referenced by: '<Root>/Constant'
                                        */
  real_T Gain_Gain;                    /* Expression: 180/pi
                                        * Referenced by: '<S9>/Gain'
                                        */
  real_T Gain1_Gain_o;    /* Expression: VehicleMass/(2*TireStiffness*(1+Lf/Lr))
                           * Referenced by: '<S8>/Gain1'
                           */
  real_T Gain1_Gain_a;                 /* Expression: 180/pi
                                        * Referenced by: '<Root>/Gain1'
                                        */
  real_T UnitDelay_InitialCondition;   /* Expression: 0
                                        * Referenced by: '<S8>/Unit Delay'
                                        */
  real_T UnitDelay_InitialCondition_g; /* Expression: 0
                                        * Referenced by: '<Root>/Unit Delay'
                                        */
  real_T Saturation_UpperSat_n;        /* Expression: MaxSteeringAngle
                                        * Referenced by: '<S8>/Saturation'
                                        */
  real_T Saturation_LowerSat_f;        /* Expression: -MaxSteeringAngle
                                        * Referenced by: '<S8>/Saturation'
                                        */
  real_T TmpRTBAtEqual1Inport1_InitialCo;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtEqual2Inport1_InitialCo;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtEqual4Inport2_InitialCo;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtEqual5Inport2_InitialCo;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T xdot_InitialCondition;        /* Expression: 0
                                        * Referenced by:
                                        */
  real_T xdot1_InitialCondition;       /* Expression: 0
                                        * Referenced by:
                                        */
  real_T TmpRTBAtMultiplyInport1_Initial;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtSwitchCaseInport1_Initi;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T xdot_InitialCondition_b;      /* Expression: 0
                                        * Referenced by:
                                        */
  real_T TmpRTBAtIntegratorInport2_Initi;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T Integrator_IC;                /* Expression: 0
                                        * Referenced by: '<S119>/Integrator'
                                        */
  real_T Gain1_Gain_d;                 /* Expression: 1.2*0.3*2/2
                                        * Referenced by: '<S119>/Gain1'
                                        */
  real_T Gain2_Gain;                   /* Expression: 0.7*9.81
                                        * Referenced by: '<S119>/Gain2'
                                        */
  real_T TmpRTBAtSumInport3_InitialCondi;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T Constant_Value_d;             /* Expression: 0
                                        * Referenced by: '<S4>/Constant'
                                        */
  real_T TmpRTBAtIntegratorInport2_Ini_j;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T Integrator_IC_g;              /* Expression: 0
                                        * Referenced by: '<S118>/Integrator'
                                        */
  real_T Gain1_Gain_p;                 /* Expression: 1.2*0.3*2/2
                                        * Referenced by: '<S118>/Gain1'
                                        */
  real_T Gain2_Gain_h;                 /* Expression: 0.7*9.81
                                        * Referenced by: '<S118>/Gain2'
                                        */
  real_T TmpRTBAtSumInport3_InitialCon_j;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtIntegratorInport1_Initi;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtIntegratorInport1_Ini_j;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T Switch_Threshold;             /* Expression: 0
                                        * Referenced by: '<S4>/Switch'
                                        */
  real_T Switch1_Threshold;            /* Expression: 0
                                        * Referenced by: '<S4>/Switch1'
                                        */
  real_T Constant9_Value;              /* Expression: 1
                                        * Referenced by: '<S5>/Constant9'
                                        */
  real_T Constant2_Value;              /* Expression: 0
                                        * Referenced by: '<S149>/Constant2'
                                        */
  real_T Gain_Gain_l;                  /* Expression: -180/pi
                                        * Referenced by: '<S120>/Gain'
                                        */
  real_T Gain5_Gain;                   /* Expression: -1
                                        * Referenced by: '<S120>/Gain5'
                                        */
  real_T Constant8_Value;              /* Expression: 0
                                        * Referenced by: '<S149>/Constant8'
                                        */
  real_T Gain6_Gain;                   /* Expression: -1
                                        * Referenced by: '<S120>/Gain6'
                                        */
  real_T Constant7_Value;              /* Expression: 0
                                        * Referenced by: '<S149>/Constant7'
                                        */
  real_T Gain1_Gain_b;                 /* Expression: -180/pi
                                        * Referenced by: '<S120>/Gain1'
                                        */
  real_T Gain2_Gain_c;                 /* Expression: -180/pi
                                        * Referenced by: '<S120>/Gain2'
                                        */
  real_T Gain3_Gain;                   /* Expression: -1
                                        * Referenced by: '<S120>/Gain3'
                                        */
  real_T Constant_Value_a;             /* Expression: 0
                                        * Referenced by: '<S149>/Constant'
                                        */
  real_T Gain4_Gain;                   /* Expression: -1
                                        * Referenced by: '<S120>/Gain4'
                                        */
  real_T RateTransition_InitialCondition;/* Expression: 0
                                          * Referenced by: '<S120>/Rate Transition'
                                          */
  real_T RateTransition1_InitialConditio;/* Expression: 0
                                          * Referenced by: '<S120>/Rate Transition1'
                                          */
  real_T RateTransition2_InitialConditio;/* Expression: 0
                                          * Referenced by: '<S120>/Rate Transition2'
                                          */
  real_T TmpRTBAtMATLABFunction1Inport1_;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtMATLABFunction1Inport2_;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtMATLABFunction1Inport3_;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtMATLABFunction2Inport1_;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtMATLABFunction2Inport2_;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtMATLABFunction2Inport3_;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtMATLABFunction3Inport1_;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtMATLABFunction3Inport2_;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtMATLABFunction3Inport3_;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtMATLABFunctionInport1_I;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtMATLABFunctionInport2_I;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtMATLABFunctionInport3_I;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T Constant2_Value_g;            /* Expression: 0
                                        * Referenced by: '<S121>/Constant2'
                                        */
  real_T FirstOrderHold_IniOut;        /* Expression: 0
                                        * Referenced by: '<S121>/First Order Hold'
                                        */
  real_T FirstOrderHold1_IniOut;       /* Expression: 0
                                        * Referenced by: '<S121>/First Order Hold1'
                                        */
  real_T Constant_Value_k;             /* Expression: 0
                                        * Referenced by: '<S129>/Constant'
                                        */
  real_T Constant5_Value;              /* Expression: 0
                                        * Referenced by: '<S149>/Constant5'
                                        */
  real_T Constant4_Value_g;            /* Expression: 0
                                        * Referenced by: '<S149>/Constant4'
                                        */
  real_T Constant6_Value;              /* Expression: 0
                                        * Referenced by: '<S149>/Constant6'
                                        */
  real_T Constant_Value_dz;            /* Expression: dh
                                        * Referenced by: '<S153>/Constant'
                                        */
  real_T Constant3_Value;              /* Expression: hh
                                        * Referenced by: '<S153>/Constant3'
                                        */
  real_T Constant4_Value_b;            /* Expression: 0
                                        * Referenced by: '<S153>/Constant4'
                                        */
  real_T Constant1_Value_o;            /* Expression: 0
                                        * Referenced by: '<S149>/Constant1'
                                        */
  real_T Constant10_Value;             /* Expression: 0
                                        * Referenced by: '<S149>/Constant10'
                                        */
  real_T Constant3_Value_j;            /* Expression: 0
                                        * Referenced by: '<S149>/Constant3'
                                        */
  real_T Constant9_Value_n;            /* Expression: 0
                                        * Referenced by: '<S149>/Constant9'
                                        */
  real_T lateral_IC;                   /* Expression: 0
                                        * Referenced by: '<S212>/lateral'
                                        */
  real_T lateral_IC_a;                 /* Expression: 0
                                        * Referenced by: '<S213>/lateral'
                                        */
  real_T secondorderlowpasstorquemodel_I;/* Expression: 0
                                          * Referenced by: '<S127>/second-order low-pass torque model'
                                          */
  real_T secondorderlowpassbrakemodel_In;/* Expression: 0
                                          * Referenced by: '<S127>/second-order low-pass brake model'
                                          */
  real_T UnitDelay_InitialCondition_f; /* Expression: 0
                                        * Referenced by: '<S127>/Unit Delay'
                                        */
  real_T secondorderlowpasssteeringmodel;/* Expression: 0
                                          * Referenced by: '<S127>/second-order low-pass steering model'
                                          */
  real_T degtorad1_Gain;               /* Expression: pi/180
                                        * Referenced by: '<S127>/deg-to-rad1'
                                        */
  real_T Gain2_Gain_h0;                /* Expression: 0.5
                                        * Referenced by: '<S127>/Gain2'
                                        */
  real_T Constant3_Value_g;            /* Expression: 1
                                        * Referenced by: '<S5>/Constant3'
                                        */
  real_T Constant4_Value_ad;           /* Expression: 1
                                        * Referenced by: '<S5>/Constant4'
                                        */
  real_T Constant_Value_h;             /* Expression: 120
                                        * Referenced by: '<S5>/Constant'
                                        */
  real_T Constant1_Value_b;            /* Expression: 0
                                        * Referenced by: '<S5>/Constant1'
                                        */
  real_T Constant2_Value_o;            /* Expression: -1000
                                        * Referenced by: '<S5>/Constant2'
                                        */
  real_T Constant5_Value_h;            /* Expression: 1
                                        * Referenced by: '<S5>/Constant5'
                                        */
  real_T Constant2_Value_e;            /* Expression: 1
                                        * Referenced by: '<Root>/Constant2'
                                        */
  real_T Constant3_Value_j2;           /* Expression: 0
                                        * Referenced by: '<Root>/Constant3'
                                        */
  real_T Constant4_Value_i;            /* Expression: 0
                                        * Referenced by: '<S12>/Constant4'
                                        */
  real_T Constant1_Value_ks;           /* Expression: 0
                                        * Referenced by: '<S12>/Constant1'
                                        */
  real_T Constant3_Value_n;            /* Expression: -1
                                        * Referenced by: '<S12>/Constant3'
                                        */
  real_T Constant2_Value_i;            /* Expression: 1
                                        * Referenced by: '<S12>/Constant2'
                                        */
  real_T Merge_InitialOutput;         /* Computed Parameter: Merge_InitialOutput
                                       * Referenced by: '<S11>/Merge'
                                       */
  real_T Constant_Value_c;             /* Expression: 0
                                        * Referenced by: '<S10>/Constant'
                                        */
  real_T Switch_Threshold_l;           /* Expression: 0
                                        * Referenced by: '<S10>/Switch'
                                        */
  real_T Switch1_Threshold_i;          /* Expression: 0
                                        * Referenced by: '<S10>/Switch1'
                                        */
  real_T Constant_Value_ke;            /* Expression: session.data.EgoCarId
                                        * Referenced by: '<S120>/Constant'
                                        */
  uint32_T uDLookupTable_numYWorkElts[2];
                               /* Computed Parameter: uDLookupTable_numYWorkElts
                                * Referenced by: '<S218>/1-D Lookup Table'
                                */
  uint32_T uDLookupTable_maxIndex; /* Computed Parameter: uDLookupTable_maxIndex
                                    * Referenced by: '<S218>/1-D Lookup Table'
                                    */
  uint32_T uDLookupTable_dimSizes; /* Computed Parameter: uDLookupTable_dimSizes
                                    * Referenced by: '<S218>/1-D Lookup Table'
                                    */
};

/* Storage class 'PageSwitching' */
extern ModelWithControllersOn_cal_type ModelWithControllersOn_cal_impl;
extern ModelWithControllersOn_cal_type *ModelWithControllersOnly_cal;

#endif                          /* RTW_HEADER_ModelWithControllersOnly_cal_h_ */
