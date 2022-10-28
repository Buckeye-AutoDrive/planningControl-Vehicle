#ifndef RTW_HEADER_speedgoat_target_model_2021b_cal_h_
#define RTW_HEADER_speedgoat_target_model_2021b_cal_h_
#include "rtwtypes.h"
#include "speedgoat_target_model_2021b_types.h"

/* Storage class 'PageSwitching', for system '<Root>' */
struct speedgoat_target_model_cal_type {
  real_T Constant_Value_kp[2];         /* Expression: [28.17 -141.5]
                                        * Referenced by: '<S1>/Constant'
                                        */
  real_T ym_zero_Value[2];             /* Expression: zeros(nym,1)
                                        * Referenced by: '<S23>/ym_zero'
                                        */
  real_T F_zero_Value[2];              /* Expression: zeros(1,2)
                                        * Referenced by: '<S22>/F_zero'
                                        */
  real_T ywt_zero_Value[2];            /* Expression: zeros(2,1)
                                        * Referenced by: '<S22>/y.wt_zero'
                                        */
  real_T ymin_scale1_Gain[2];       /* Expression: Yscale(:,ones(1,max(nCC,1)))'
                                     * Referenced by: '<S23>/ymin_scale1'
                                     */
  real_T last_x_InitialCondition[4];   /* Expression: lastx+xoff
                                        * Referenced by: '<S23>/last_x'
                                        */
  real_T umin_scale3_Gain[11];         /* Expression: MVscale(:,ones(1,p+1))'
                                        * Referenced by: '<S23>/umin_scale3'
                                        */
  real_T umin_scale5_Gain[22];         /* Expression: Yscale(:,ones(1,p+1))'
                                        * Referenced by: '<S23>/umin_scale5'
                                        */
  real_T ACCsystem_DefaultSpacing;   /* Mask Parameter: ACCsystem_DefaultSpacing
                                      * Referenced by: '<S4>/Default spacing constant'
                                      */
  real_T LateralControllerStanley_DelayG;
                              /* Mask Parameter: LateralControllerStanley_DelayG
                               * Referenced by: '<S47>/Gain2'
                               */
  real_T PIForward_InitialConditionForIn;
                              /* Mask Parameter: PIForward_InitialConditionForIn
                               * Referenced by: '<S88>/Integrator'
                               */
  real_T PIReverse_InitialConditionForIn;
                              /* Mask Parameter: PIReverse_InitialConditionForIn
                               * Referenced by: '<S139>/Integrator'
                               */
  real_T LongitudinalControllerStanley_K;
                              /* Mask Parameter: LongitudinalControllerStanley_K
                               * Referenced by:
                               *   '<S85>/Integral Gain'
                               *   '<S136>/Integral Gain'
                               */
  real_T LongitudinalControllerStanley_n;
                              /* Mask Parameter: LongitudinalControllerStanley_n
                               * Referenced by:
                               *   '<S93>/Proportional Gain'
                               *   '<S144>/Proportional Gain'
                               */
  real_T PIForward_LowerSaturationLimit;
                               /* Mask Parameter: PIForward_LowerSaturationLimit
                                * Referenced by:
                                *   '<S95>/Saturation'
                                *   '<S81>/DeadZone'
                                */
  real_T PIReverse_LowerSaturationLimit;
                               /* Mask Parameter: PIReverse_LowerSaturationLimit
                                * Referenced by:
                                *   '<S146>/Saturation'
                                *   '<S132>/DeadZone'
                                */
  real_T ACCsystem_MaxAcceleration; /* Mask Parameter: ACCsystem_MaxAcceleration
                                     * Referenced by: '<S4>/Maximum longitudinal acceleration constant'
                                     */
  real_T ACCsystem_MaxVelocity;        /* Mask Parameter: ACCsystem_MaxVelocity
                                        * Referenced by: '<S4>/Maximum velocity constant'
                                        */
  real_T ACCsystem_MinAcceleration; /* Mask Parameter: ACCsystem_MinAcceleration
                                     * Referenced by: '<S4>/Minimum longitudinal acceleration constant'
                                     */
  real_T LateralControllerStanley_Positi;
                              /* Mask Parameter: LateralControllerStanley_Positi
                               * Referenced by: '<S7>/Kinematic'
                               */
  real_T LateralControllerStanley_Posi_i;
                              /* Mask Parameter: LateralControllerStanley_Posi_i
                               * Referenced by: '<S7>/Kinematic'
                               */
  real_T PIForward_UpperSaturationLimit;
                               /* Mask Parameter: PIForward_UpperSaturationLimit
                                * Referenced by:
                                *   '<S95>/Saturation'
                                *   '<S81>/DeadZone'
                                */
  real_T PIReverse_UpperSaturationLimit;
                               /* Mask Parameter: PIReverse_UpperSaturationLimit
                                * Referenced by:
                                *   '<S146>/Saturation'
                                *   '<S132>/DeadZone'
                                */
  real_T LateralControllerStanley_YawRat;
                              /* Mask Parameter: LateralControllerStanley_YawRat
                               * Referenced by: '<S47>/Gain'
                               */
  real_T CompareToConstant_const;     /* Mask Parameter: CompareToConstant_const
                                       * Referenced by: '<S175>/Constant'
                                       */
  real_T CompareToConstant_const_e; /* Mask Parameter: CompareToConstant_const_e
                                     * Referenced by: '<S176>/Constant'
                                     */
  real_T CompareToConstant_const_l; /* Mask Parameter: CompareToConstant_const_l
                                     * Referenced by: '<S177>/Constant'
                                     */
  real_T CompareToConstant_const_k; /* Mask Parameter: CompareToConstant_const_k
                                     * Referenced by: '<S178>/Constant'
                                     */
  real_T CompareToConstant_const_g; /* Mask Parameter: CompareToConstant_const_g
                                     * Referenced by: '<S179>/Constant'
                                     */
  real_T Kinematic_MaxSteeringAngle;   /* Expression: MaxSteeringAngle
                                        * Referenced by: '<S7>/Kinematic'
                                        */
  real_T Kinematic_Wheelbase;          /* Expression: Wheelbase
                                        * Referenced by: '<S7>/Kinematic'
                                        */
  real_T Constant1_Value;              /* Expression: 0
                                        * Referenced by: '<S79>/Constant1'
                                        */
  real_T Integrator_gainval;           /* Computed Parameter: Integrator_gainval
                                        * Referenced by: '<S88>/Integrator'
                                        */
  real_T ZeroGain_Gain;                /* Expression: 0
                                        * Referenced by: '<S79>/ZeroGain'
                                        */
  real_T Constant1_Value_f;            /* Expression: 0
                                        * Referenced by: '<S130>/Constant1'
                                        */
  real_T Integrator_gainval_b;       /* Computed Parameter: Integrator_gainval_b
                                      * Referenced by: '<S139>/Integrator'
                                      */
  real_T ZeroGain_Gain_k;              /* Expression: 0
                                        * Referenced by: '<S130>/ZeroGain'
                                        */
  real_T DiscreteTimeIntegrator_gainval;
                           /* Computed Parameter: DiscreteTimeIntegrator_gainval
                            * Referenced by: '<S164>/Discrete-Time Integrator'
                            */
  real_T DiscreteTimeIntegrator_IC;    /* Expression: 0
                                        * Referenced by: '<S164>/Discrete-Time Integrator'
                                        */
  real_T Gain1_Gain;                   /* Expression: 1.2*0.3*2/2
                                        * Referenced by: '<S164>/Gain1'
                                        */
  real_T Constant1_Value_i;            /* Expression: 1853
                                        * Referenced by: '<S164>/Constant1'
                                        */
  real_T Gain2_Gain;                   /* Expression: 0.7*9.81
                                        * Referenced by: '<S164>/Gain2'
                                        */
  real_T Gain_Gain;                    /* Expression: 1.1*1853
                                        * Referenced by: '<S164>/Gain'
                                        */
  real_T Gain4_Gain;                   /* Expression: -0.2159/(7.05*0.8)
                                        * Referenced by: '<S164>/Gain4'
                                        */
  real_T BrakeSaturation_UpperSat;     /* Expression: 0
                                        * Referenced by: '<S157>/Brake Saturation'
                                        */
  real_T BrakeSaturation_LowerSat;     /* Expression: -65534
                                        * Referenced by: '<S157>/Brake Saturation'
                                        */
  real_T Constant_Value;               /* Expression: 0
                                        * Referenced by: '<S157>/Constant'
                                        */
  real_T Constant1_Value_d;            /* Expression: 1853
                                        * Referenced by: '<S163>/Constant1'
                                        */
  real_T DiscreteTimeIntegrator_gainva_j;
                          /* Computed Parameter: DiscreteTimeIntegrator_gainva_j
                           * Referenced by: '<S163>/Discrete-Time Integrator'
                           */
  real_T DiscreteTimeIntegrator_IC_a;  /* Expression: 0
                                        * Referenced by: '<S163>/Discrete-Time Integrator'
                                        */
  real_T Gain_Gain_n;                  /* Expression: 1.1*1853
                                        * Referenced by: '<S163>/Gain'
                                        */
  real_T Gain1_Gain_d;                 /* Expression: 1.2*0.3*2/2
                                        * Referenced by: '<S163>/Gain1'
                                        */
  real_T Gain2_Gain_n;                 /* Expression: 0.7*9.81
                                        * Referenced by: '<S163>/Gain2'
                                        */
  real_T Gain4_Gain_b;                 /* Expression: 0.2159/(7.05*0.8)
                                        * Referenced by: '<S163>/Gain4'
                                        */
  real_T ThrottleSaturation_UpperSat;  /* Expression: 22534
                                        * Referenced by: '<S157>/Throttle Saturation'
                                        */
  real_T ThrottleSaturation_LowerSat;  /* Expression: -22534
                                        * Referenced by: '<S157>/Throttle Saturation'
                                        */
  real_T Switch_Threshold;             /* Expression: 500
                                        * Referenced by: '<S157>/Switch'
                                        */
  real_T Switch1_Threshold;            /* Expression: -500
                                        * Referenced by: '<S157>/Switch1'
                                        */
  real_T DiscreteTimeIntegrator_gainva_c;
                          /* Computed Parameter: DiscreteTimeIntegrator_gainva_c
                           * Referenced by: '<S160>/Discrete-Time Integrator'
                           */
  real_T DiscreteTimeIntegrator_IC_i;  /* Expression: 0
                                        * Referenced by: '<S160>/Discrete-Time Integrator'
                                        */
  real_T Gain1_Gain_m;                 /* Expression: 1.2*0.3*2/2
                                        * Referenced by: '<S160>/Gain1'
                                        */
  real_T Constant1_Value_o;            /* Expression: 1853
                                        * Referenced by: '<S160>/Constant1'
                                        */
  real_T Gain2_Gain_l;                 /* Expression: 0.7*9.81
                                        * Referenced by: '<S160>/Gain2'
                                        */
  real_T Gain_Gain_m;                  /* Expression: 1.1*1853
                                        * Referenced by: '<S160>/Gain'
                                        */
  real_T Gain4_Gain_m;                 /* Expression: -0.2159/(7.05*0.8)
                                        * Referenced by: '<S160>/Gain4'
                                        */
  real_T BrakeSaturation_UpperSat_a;   /* Expression: 0
                                        * Referenced by: '<S156>/Brake Saturation'
                                        */
  real_T BrakeSaturation_LowerSat_m;   /* Expression: -65534
                                        * Referenced by: '<S156>/Brake Saturation'
                                        */
  real_T Constant_Value_d;             /* Expression: 0
                                        * Referenced by: '<S156>/Constant'
                                        */
  real_T Constant1_Value_fy;           /* Expression: 1853
                                        * Referenced by: '<S159>/Constant1'
                                        */
  real_T DiscreteTimeIntegrator_gainva_a;
                          /* Computed Parameter: DiscreteTimeIntegrator_gainva_a
                           * Referenced by: '<S159>/Discrete-Time Integrator'
                           */
  real_T DiscreteTimeIntegrator_IC_d;  /* Expression: 0
                                        * Referenced by: '<S159>/Discrete-Time Integrator'
                                        */
  real_T Gain_Gain_c;                  /* Expression: 1.1*1853
                                        * Referenced by: '<S159>/Gain'
                                        */
  real_T Gain1_Gain_j;                 /* Expression: 1.2*0.3*2/2
                                        * Referenced by: '<S159>/Gain1'
                                        */
  real_T Gain2_Gain_h;                 /* Expression: 0.7*9.81
                                        * Referenced by: '<S159>/Gain2'
                                        */
  real_T Gain4_Gain_n;                 /* Expression: 0.2159/(7.05*0.8)
                                        * Referenced by: '<S159>/Gain4'
                                        */
  real_T ThrottleSaturation_UpperSat_e;/* Expression: 22534
                                        * Referenced by: '<S156>/Throttle Saturation'
                                        */
  real_T ThrottleSaturation_LowerSat_g;/* Expression: -22534
                                        * Referenced by: '<S156>/Throttle Saturation'
                                        */
  real_T Switch_Threshold_f;           /* Expression: 500
                                        * Referenced by: '<S156>/Switch'
                                        */
  real_T Switch1_Threshold_m;          /* Expression: -500
                                        * Referenced by: '<S156>/Switch1'
                                        */
  real_T CurrPoses_Y0;                 /* Computed Parameter: CurrPoses_Y0
                                        * Referenced by: '<S165>/CurrPoses'
                                        */
  real_T ObjAhead_Y0;                  /* Computed Parameter: ObjAhead_Y0
                                        * Referenced by: '<S166>/Obj Ahead'
                                        */
  real_T Longitudinalvelocity_Y0; /* Computed Parameter: Longitudinalvelocity_Y0
                                   * Referenced by: '<S166>/Longitudinal velocity'
                                   */
  real_T SceneID_Y0;                   /* Computed Parameter: SceneID_Y0
                                        * Referenced by: '<S166>/SceneID'
                                        */
  real_T Relativevelocity_Y0;         /* Computed Parameter: Relativevelocity_Y0
                                       * Referenced by: '<S166>/Relative velocity'
                                       */
  real_T RefPoses_Y0;                  /* Computed Parameter: RefPoses_Y0
                                        * Referenced by: '<S167>/RefPoses'
                                        */
  real_T ObjectDist_Y0;                /* Computed Parameter: ObjectDist_Y0
                                        * Referenced by: '<S168>/ObjectDist'
                                        */
  real_T D2S_Y0;                       /* Computed Parameter: D2S_Y0
                                        * Referenced by: '<S168>/D2S'
                                        */
  real_T Relativedistance_Y0;         /* Computed Parameter: Relativedistance_Y0
                                       * Referenced by: '<S168>/Relative distance'
                                       */
  real_T TLstate_Y0;                   /* Computed Parameter: TLstate_Y0
                                        * Referenced by: '<S168>/TL state'
                                        */
  real_T Direction_Y0;                 /* Computed Parameter: Direction_Y0
                                        * Referenced by: '<S169>/Direction'
                                        */
  real_T Curvature_Y0;                 /* Computed Parameter: Curvature_Y0
                                        * Referenced by: '<S169>/Curvature'
                                        */
  real_T yawrate_Y0;                   /* Computed Parameter: yawrate_Y0
                                        * Referenced by: '<S169>/yaw rate'
                                        */
  real_T NumMsg_Y0;                    /* Computed Parameter: NumMsg_Y0
                                        * Referenced by: '<S170>/NumMsg'
                                        */
  real_T Constant_Value_k;             /* Expression: 1
                                        * Referenced by: '<S170>/Constant'
                                        */
  real_T NumMsg_Y0_b;                  /* Computed Parameter: NumMsg_Y0_b
                                        * Referenced by: '<S171>/NumMsg'
                                        */
  real_T Constant_Value_dy;            /* Expression: 1
                                        * Referenced by: '<S171>/Constant'
                                        */
  real_T NumMsg_Y0_bf;                 /* Computed Parameter: NumMsg_Y0_bf
                                        * Referenced by: '<S172>/NumMsg'
                                        */
  real_T Constant_Value_p;             /* Expression: 1
                                        * Referenced by: '<S172>/Constant'
                                        */
  real_T NumMsg_Y0_f;                  /* Computed Parameter: NumMsg_Y0_f
                                        * Referenced by: '<S173>/NumMsg'
                                        */
  real_T Constant_Value_g;             /* Expression: 1
                                        * Referenced by: '<S173>/Constant'
                                        */
  real_T NumMsg_Y0_c;                  /* Computed Parameter: NumMsg_Y0_c
                                        * Referenced by: '<S174>/NumMsg'
                                        */
  real_T Constant_Value_j;             /* Expression: 1
                                        * Referenced by: '<S174>/Constant'
                                        */
  real_T Gain_Gain_g;                  /* Expression: 180/pi
                                        * Referenced by: '<S48>/Gain'
                                        */
  real_T Gain1_Gain_j3;   /* Expression: VehicleMass/(2*TireStiffness*(1+Lf/Lr))
                           * Referenced by: '<S47>/Gain1'
                           */
  real_T Gain_Gain_k;                  /* Expression: 180/pi
                                        * Referenced by: '<S1>/Gain'
                                        */
  real_T UnitDelay_InitialCondition;   /* Expression: 0
                                        * Referenced by: '<S47>/Unit Delay'
                                        */
  real_T UnitDelay_InitialCondition_j; /* Expression: 0
                                        * Referenced by: '<S1>/Unit Delay'
                                        */
  real_T Saturation_UpperSat;          /* Expression: MaxSteeringAngle
                                        * Referenced by: '<S47>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: -MaxSteeringAngle
                                        * Referenced by: '<S47>/Saturation'
                                        */
  real_T ApproachingSpeed_Value;       /* Expression: 5
                                        * Referenced by: '<S1>/Approaching Speed'
                                        */
  real_T SceneSpeed_Value;             /* Expression: 5
                                        * Referenced by: '<S1>/Scene Speed'
                                        */
  real_T Constant1_Value_i4;           /* Expression: 25
                                        * Referenced by: '<S1>/Constant1'
                                        */
  real_T RRCrossing_Value;             /* Expression: 0
                                        * Referenced by: '<S1>/RR Crossing'
                                        */
  real_T NavigatingFlag_Value;         /* Expression: 1
                                        * Referenced by: '<S1>/Navigating Flag'
                                        */
  real_T Clear2Go_Value;               /* Expression: 1
                                        * Referenced by: '<S1>/Clear2Go'
                                        */
  real_T RoadSpeed_Value;              /* Expression: 15
                                        * Referenced by: '<S1>/RoadSpeed'
                                        */
  real_T WaitTime_Value;               /* Expression: 3
                                        * Referenced by: '<S1>/WaitTime'
                                        */
  real_T MinSpeed_Value;               /* Expression: 2
                                        * Referenced by: '<S1>/MinSpeed'
                                        */
  real_T LaneAvailable_Value;          /* Expression: 0
                                        * Referenced by: '<S1>/LaneAvailable'
                                        */
  real_T LaneChanged_Value;            /* Expression: 0
                                        * Referenced by: '<S1>/LaneChanged'
                                        */
  real_T acceleration_InitialCondition;/* Expression: 0
                                        * Referenced by:
                                        */
  real_T TmpRTBAtSwitchOutport1_InitialC;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtSwitch1Outport1_Initial;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T Merge1_InitialOutput;       /* Computed Parameter: Merge1_InitialOutput
                                      * Referenced by: '<S11>/Merge1'
                                      */
  real_T Merge2_InitialOutput;       /* Computed Parameter: Merge2_InitialOutput
                                      * Referenced by: '<S11>/Merge2'
                                      */
  real_T p_zero_Value;                 /* Expression: zeros(1,1)
                                        * Referenced by: '<S22>/p_zero'
                                        */
  real_T m_zero_Value;                 /* Expression: zeros(1,1)
                                        * Referenced by: '<S22>/m_zero'
                                        */
  real_T TmpRTBAtoptimizerInport3_Initia;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtoptimizerInport5_Initia;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T Minimumvelocityconstant_Value;/* Expression: MinVelocity
                                        * Referenced by: '<S4>/Minimum velocity constant'
                                        */
  real_T min_velocity_InitialCondition;/* Expression: 0
                                        * Referenced by:
                                        */
  real_T TmpRTBAtDataTypeConversion_vset;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T default_spacing_InitialConditio;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T Externalcontrolsignalconstant_V;/* Expression: 0
                                          * Referenced by: '<S4>/External control signal constant'
                                          */
  real_T TmpRTBAtEqual1Inport1_InitialCo;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtEqual2Inport1_InitialCo;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtEqual3Inport1_InitialCo;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtEqual4Inport2_InitialCo;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtEqual5Inport2_InitialCo;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtSign1Inport1_InitialCon;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T set_velocity_InitialCondition;/* Expression: 0
                                        * Referenced by:
                                        */
  real_T set_velocity1_InitialCondition;/* Expression: 0
                                         * Referenced by:
                                         */
  real_T TmpRTBAtMultiplyInport1_Initial;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtSwitchCaseInport1_Initi;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T TmpRTBAtMinusInport2_InitialCon;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T set_velocity_InitialCondition_h;/* Expression: 0
                                          * Referenced by:
                                          */
  real_T Enableoptimizationconstant_Valu;/* Expression: 0
                                          * Referenced by: '<S4>/Enable optimization constant'
                                          */
  real_T E_zero_Value;                 /* Expression: zeros(1,1)
                                        * Referenced by: '<S22>/E_zero'
                                        */
  real_T G_zero_Value;                 /* Expression: zeros(1,1)
                                        * Referenced by: '<S22>/G_zero'
                                        */
  real_T uwt_zero_Value;               /* Expression: zeros(1,1)
                                        * Referenced by: '<S22>/u.wt_zero'
                                        */
  real_T duwt_zero_Value;              /* Expression: zeros(1,1)
                                        * Referenced by: '<S22>/du.wt_zero'
                                        */
  real_T umin_scale4_Gain;         /* Expression: MVscale(:,ones(1,max(nCC,1)))'
                                    * Referenced by: '<S23>/umin_scale4'
                                    */
  real_T S_zero_Value;                 /* Expression: zeros(1,1)
                                        * Referenced by: '<S22>/S_zero'
                                        */
  real_T ymin_scale2_Gain;         /* Expression: MDscale(:,ones(1,max(nCC,1)))'
                                    * Referenced by: '<S23>/ymin_scale2'
                                    */
  real_T extmv_zero_Value;             /* Expression: zeros(1,1)
                                        * Referenced by: '<S22>/ext.mv_zero'
                                        */
  real_T extmv_scale_Gain;             /* Expression: RMVscale
                                        * Referenced by: '<S23>/ext.mv_scale'
                                        */
  real_T mvtarget_zero_Value;          /* Expression: zeros(1,1)
                                        * Referenced by: '<S22>/mv.target_zero'
                                        */
  real_T extmv_scale1_Gain;            /* Expression: RMVscale
                                        * Referenced by: '<S23>/ext.mv_scale1'
                                        */
  real_T last_mv_InitialCondition;     /* Expression: lastu+uoff
                                        * Referenced by: '<S23>/last_mv'
                                        */
  real_T Unconstrained_Value;          /* Expression: 0
                                        * Referenced by: '<S4>/Unconstrained'
                                        */
  real_T ecrwt_zero_Value;             /* Expression: zeros(1,1)
                                        * Referenced by: '<S22>/ecr.wt_zero'
                                        */
  real_T umin_scale1_Gain;             /* Expression: MVscale
                                        * Referenced by: '<S23>/umin_scale1'
                                        */
  real_T constant_Value;               /* Expression: lastu+uoff
                                        * Referenced by: '<S23>/constant'
                                        */
  real_T umin_scale2_Gain;             /* Expression: MVscale
                                        * Referenced by: '<S23>/umin_scale2'
                                        */
  real_T Constant4_Value;              /* Expression: 0
                                        * Referenced by: '<S51>/Constant4'
                                        */
  real_T Constant1_Value_id;           /* Expression: 0
                                        * Referenced by: '<S51>/Constant1'
                                        */
  real_T Constant3_Value;              /* Expression: -1
                                        * Referenced by: '<S51>/Constant3'
                                        */
  real_T Constant2_Value;              /* Expression: 1
                                        * Referenced by: '<S51>/Constant2'
                                        */
  real_T reset_Value;                  /* Expression: 0
                                        * Referenced by: '<S1>/reset'
                                        */
  real_T Merge_InitialOutput;         /* Computed Parameter: Merge_InitialOutput
                                       * Referenced by: '<S50>/Merge'
                                       */
  real_T Constant_Value_m;             /* Expression: 0
                                        * Referenced by: '<S49>/Constant'
                                        */
  real_T Switch_Threshold_k;           /* Expression: 0
                                        * Referenced by: '<S49>/Switch'
                                        */
  real_T Switch1_Threshold_g;          /* Expression: 0
                                        * Referenced by: '<S49>/Switch1'
                                        */
  uint16_T MatrixDimensionCheck_P1[64];
                                  /* Computed Parameter: MatrixDimensionCheck_P1
                                   * Referenced by: '<S28>/Matrix Dimension Check'
                                   */
  uint16_T MatrixDimensionCheck_P2[2];
                                  /* Computed Parameter: MatrixDimensionCheck_P2
                                   * Referenced by: '<S28>/Matrix Dimension Check'
                                   */
  uint16_T MatrixDimensionCheck_P1_j[64];
                                /* Computed Parameter: MatrixDimensionCheck_P1_j
                                 * Referenced by: '<S34>/Matrix Dimension Check'
                                 */
  uint16_T MatrixDimensionCheck_P2_j[4];
                                /* Computed Parameter: MatrixDimensionCheck_P2_j
                                 * Referenced by: '<S34>/Matrix Dimension Check'
                                 */
  uint16_T VectorDimensionCheck_P1[64];
                                  /* Computed Parameter: VectorDimensionCheck_P1
                                   * Referenced by: '<S38>/Vector Dimension Check'
                                   */
  uint16_T VectorDimensionCheck_P1_b[64];
                                /* Computed Parameter: VectorDimensionCheck_P1_b
                                 * Referenced by: '<S42>/Vector Dimension Check'
                                 */
  uint16_T VectorDimensionCheck_P2_l[2];
                                /* Computed Parameter: VectorDimensionCheck_P2_l
                                 * Referenced by: '<S42>/Vector Dimension Check'
                                 */
  uint16_T MatrixDimensionCheck_P1_p[64];
                                /* Computed Parameter: MatrixDimensionCheck_P1_p
                                 * Referenced by: '<S24>/Matrix Dimension Check'
                                 */
  uint16_T MatrixDimensionCheck_P2_c[18];
                                /* Computed Parameter: MatrixDimensionCheck_P2_c
                                 * Referenced by: '<S24>/Matrix Dimension Check'
                                 */
  uint16_T MatrixDimensionCheck_P1_c[64];
                                /* Computed Parameter: MatrixDimensionCheck_P1_c
                                 * Referenced by: '<S25>/Matrix Dimension Check'
                                 */
  uint16_T MatrixDimensionCheck_P2_i[18];
                                /* Computed Parameter: MatrixDimensionCheck_P2_i
                                 * Referenced by: '<S25>/Matrix Dimension Check'
                                 */
  uint16_T MatrixDimensionCheck_P1_i[64];
                                /* Computed Parameter: MatrixDimensionCheck_P1_i
                                 * Referenced by: '<S26>/Matrix Dimension Check'
                                 */
  uint16_T MatrixDimensionCheck_P2_e[18];
                                /* Computed Parameter: MatrixDimensionCheck_P2_e
                                 * Referenced by: '<S26>/Matrix Dimension Check'
                                 */
  uint16_T MatrixDimensionCheck_P1_cd[64];
                               /* Computed Parameter: MatrixDimensionCheck_P1_cd
                                * Referenced by: '<S27>/Matrix Dimension Check'
                                */
  uint16_T MatrixDimensionCheck_P2_k[3];
                                /* Computed Parameter: MatrixDimensionCheck_P2_k
                                 * Referenced by: '<S27>/Matrix Dimension Check'
                                 */
  uint16_T MatrixDimensionCheck_P1_d[64];
                                /* Computed Parameter: MatrixDimensionCheck_P1_d
                                 * Referenced by: '<S29>/Matrix Dimension Check'
                                 */
  uint16_T MatrixDimensionCheck_P2_j0[4];
                               /* Computed Parameter: MatrixDimensionCheck_P2_j0
                                * Referenced by: '<S29>/Matrix Dimension Check'
                                */
  uint16_T MatrixDimensionCheck_P1_dw[64];
                               /* Computed Parameter: MatrixDimensionCheck_P1_dw
                                * Referenced by: '<S30>/Matrix Dimension Check'
                                */
  uint16_T MatrixDimensionCheck_P2_d[4];
                                /* Computed Parameter: MatrixDimensionCheck_P2_d
                                 * Referenced by: '<S30>/Matrix Dimension Check'
                                 */
  uint16_T MatrixDimensionCheck_P1_jj[64];
                               /* Computed Parameter: MatrixDimensionCheck_P1_jj
                                * Referenced by: '<S31>/Matrix Dimension Check'
                                */
  uint16_T MatrixDimensionCheck_P2_ie[5];
                               /* Computed Parameter: MatrixDimensionCheck_P2_ie
                                * Referenced by: '<S31>/Matrix Dimension Check'
                                */
  uint16_T MatrixDimensionCheck_P1_o[64];
                                /* Computed Parameter: MatrixDimensionCheck_P1_o
                                 * Referenced by: '<S32>/Matrix Dimension Check'
                                 */
  uint16_T MatrixDimensionCheck_P2_n[4];
                                /* Computed Parameter: MatrixDimensionCheck_P2_n
                                 * Referenced by: '<S32>/Matrix Dimension Check'
                                 */
  uint16_T MatrixDimensionCheck_P1_l[64];
                                /* Computed Parameter: MatrixDimensionCheck_P1_l
                                 * Referenced by: '<S33>/Matrix Dimension Check'
                                 */
  uint16_T MatrixDimensionCheck_P2_d3[4];
                               /* Computed Parameter: MatrixDimensionCheck_P2_d3
                                * Referenced by: '<S33>/Matrix Dimension Check'
                                */
  uint16_T MatrixDimensionCheck_P1_jjz[64];
                              /* Computed Parameter: MatrixDimensionCheck_P1_jjz
                               * Referenced by: '<S35>/Matrix Dimension Check'
                               */
  uint16_T MatrixDimensionCheck_P2_g[4];
                                /* Computed Parameter: MatrixDimensionCheck_P2_g
                                 * Referenced by: '<S35>/Matrix Dimension Check'
                                 */
  uint16_T VectorDimensionCheck_P1_f[64];
                                /* Computed Parameter: VectorDimensionCheck_P1_f
                                 * Referenced by: '<S36>/Vector Dimension Check'
                                 */
  uint16_T VectorDimensionCheck_P2_b[6];
                                /* Computed Parameter: VectorDimensionCheck_P2_b
                                 * Referenced by: '<S36>/Vector Dimension Check'
                                 */
  uint16_T VectorDimensionCheck_P1_fq[64];
                               /* Computed Parameter: VectorDimensionCheck_P1_fq
                                * Referenced by: '<S37>/Vector Dimension Check'
                                */
  uint16_T VectorDimensionCheck_P2_n[6];
                                /* Computed Parameter: VectorDimensionCheck_P2_n
                                 * Referenced by: '<S37>/Vector Dimension Check'
                                 */
  uint16_T VectorDimensionCheck_P1_h[64];
                                /* Computed Parameter: VectorDimensionCheck_P1_h
                                 * Referenced by: '<S39>/Vector Dimension Check'
                                 */
  uint16_T VectorDimensionCheck_P2_h[6];
                                /* Computed Parameter: VectorDimensionCheck_P2_h
                                 * Referenced by: '<S39>/Vector Dimension Check'
                                 */
  uint16_T VectorDimensionCheck_P1_o[64];
                                /* Computed Parameter: VectorDimensionCheck_P1_o
                                 * Referenced by: '<S40>/Vector Dimension Check'
                                 */
  uint16_T VectorDimensionCheck_P2_p[18];
                                /* Computed Parameter: VectorDimensionCheck_P2_p
                                 * Referenced by: '<S40>/Vector Dimension Check'
                                 */
  uint16_T VectorDimensionCheck_P1_oo[64];
                               /* Computed Parameter: VectorDimensionCheck_P1_oo
                                * Referenced by: '<S41>/Vector Dimension Check'
                                */
  uint16_T VectorDimensionCheck_P2_f[9];
                                /* Computed Parameter: VectorDimensionCheck_P2_f
                                 * Referenced by: '<S41>/Vector Dimension Check'
                                 */
  uint16_T UDPSendControllerOutput_toPort;
                           /* Computed Parameter: UDPSendControllerOutput_toPort
                            * Referenced by: '<S3>/UDP Send: Controller Output'
                            */
  uint16_T UDPSendControllerIntermediateSi;
                          /* Computed Parameter: UDPSendControllerIntermediateSi
                           * Referenced by: '<S3>/UDP Send: Controller Intermediate Signals'
                           */
  uint16_T MatrixDimensionCheck_P3;    /* Expression: nrow
                                        * Referenced by: '<S28>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P4;    /* Expression: ncol
                                        * Referenced by: '<S28>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P3_a;  /* Expression: nrow
                                        * Referenced by: '<S34>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P4_m;  /* Expression: ncol
                                        * Referenced by: '<S34>/Matrix Dimension Check'
                                        */
  uint16_T VectorDimensionCheck_P2;
                                  /* Computed Parameter: VectorDimensionCheck_P2
                                   * Referenced by: '<S38>/Vector Dimension Check'
                                   */
  uint16_T VectorDimensionCheck_P3;    /* Expression: n
                                        * Referenced by: '<S42>/Vector Dimension Check'
                                        */
  uint16_T VectorDimensionCheck_P4;    /* Expression: option
                                        * Referenced by: '<S42>/Vector Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P3_d;  /* Expression: nrow
                                        * Referenced by: '<S24>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P4_e;  /* Expression: ncol
                                        * Referenced by: '<S24>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P5;    /* Expression: nsteps
                                        * Referenced by: '<S24>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P6;    /* Expression: isltv
                                        * Referenced by: '<S24>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P3_n;  /* Expression: nrow
                                        * Referenced by: '<S25>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P4_p;  /* Expression: ncol
                                        * Referenced by: '<S25>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P5_m;  /* Expression: nsteps
                                        * Referenced by: '<S25>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P6_d;  /* Expression: isltv
                                        * Referenced by: '<S25>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P3_m;  /* Expression: nrow
                                        * Referenced by: '<S26>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P4_ms; /* Expression: ncol
                                        * Referenced by: '<S26>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P5_h;  /* Expression: nsteps
                                        * Referenced by: '<S26>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P6_l;  /* Expression: isltv
                                        * Referenced by: '<S26>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P3_l;  /* Expression: nrow
                                        * Referenced by: '<S27>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P4_k;  /* Expression: ncol
                                        * Referenced by: '<S27>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P3_i;  /* Expression: nrow
                                        * Referenced by: '<S29>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P4_j;  /* Expression: ncol
                                        * Referenced by: '<S29>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P3_k;  /* Expression: nrow
                                        * Referenced by: '<S30>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P4_c;  /* Expression: ncol
                                        * Referenced by: '<S30>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P3_p;  /* Expression: nrow
                                        * Referenced by: '<S31>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P4_j2; /* Expression: ncol
                                        * Referenced by: '<S31>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P3_p0; /* Expression: nrow
                                        * Referenced by: '<S32>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P4_i;  /* Expression: ncol
                                        * Referenced by: '<S32>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P3_e;  /* Expression: nrow
                                        * Referenced by: '<S33>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P4_jh; /* Expression: ncol
                                        * Referenced by: '<S33>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P3_e2; /* Expression: nrow
                                        * Referenced by: '<S35>/Matrix Dimension Check'
                                        */
  uint16_T MatrixDimensionCheck_P4_jd; /* Expression: ncol
                                        * Referenced by: '<S35>/Matrix Dimension Check'
                                        */
  uint16_T VectorDimensionCheck_P3_i;  /* Expression: n
                                        * Referenced by: '<S39>/Vector Dimension Check'
                                        */
  uint16_T VectorDimensionCheck_P4_c;  /* Expression: option
                                        * Referenced by: '<S39>/Vector Dimension Check'
                                        */
  uint16_T VectorDimensionCheck_P3_d;  /* Expression: n
                                        * Referenced by: '<S40>/Vector Dimension Check'
                                        */
  uint16_T VectorDimensionCheck_P4_f;  /* Expression: option
                                        * Referenced by: '<S40>/Vector Dimension Check'
                                        */
  uint16_T VectorDimensionCheck_P3_dg; /* Expression: n
                                        * Referenced by: '<S41>/Vector Dimension Check'
                                        */
  uint16_T VectorDimensionCheck_P4_b;  /* Expression: option
                                        * Referenced by: '<S41>/Vector Dimension Check'
                                        */
  uint8_T UDPReceive2_fmAddress[4]; /* Computed Parameter: UDPReceive2_fmAddress
                                     * Referenced by: '<S170>/UDP Receive2'
                                     */
  uint8_T UDPReceive2_fmAddress_e[4];
                                  /* Computed Parameter: UDPReceive2_fmAddress_e
                                   * Referenced by: '<S171>/UDP Receive2'
                                   */
  uint8_T UDPReceive2_fmAddress_c[4];
                                  /* Computed Parameter: UDPReceive2_fmAddress_c
                                   * Referenced by: '<S172>/UDP Receive2'
                                   */
  uint8_T UDPReceive2_fmAddress_h[4];
                                  /* Computed Parameter: UDPReceive2_fmAddress_h
                                   * Referenced by: '<S173>/UDP Receive2'
                                   */
  uint8_T UDPReceive2_fmAddress_i[4];
                                  /* Computed Parameter: UDPReceive2_fmAddress_i
                                   * Referenced by: '<S174>/UDP Receive2'
                                   */
  uint8_T UDPSendControllerOutput_toAddre[4];
                          /* Computed Parameter: UDPSendControllerOutput_toAddre
                           * Referenced by: '<S3>/UDP Send: Controller Output'
                           */
  uint8_T UDPSendControllerIntermediate_k[4];
                          /* Computed Parameter: UDPSendControllerIntermediate_k
                           * Referenced by: '<S3>/UDP Send: Controller Intermediate Signals'
                           */
  uint8_T Data_Y0;                     /* Computed Parameter: Data_Y0
                                        * Referenced by: '<S170>/Data'
                                        */
  uint8_T Data_Y0_j;                   /* Computed Parameter: Data_Y0_j
                                        * Referenced by: '<S171>/Data'
                                        */
  uint8_T Data_Y0_k;                   /* Computed Parameter: Data_Y0_k
                                        * Referenced by: '<S172>/Data'
                                        */
  uint8_T Data_Y0_ja;                  /* Computed Parameter: Data_Y0_ja
                                        * Referenced by: '<S173>/Data'
                                        */
  uint8_T Data_Y0_l;                   /* Computed Parameter: Data_Y0_l
                                        * Referenced by: '<S174>/Data'
                                        */
  boolean_T Memory_InitialCondition[34];/* Expression: iA
                                         * Referenced by: '<S23>/Memory'
                                         */
};

/* Storage class 'PageSwitching' */
extern speedgoat_target_model_cal_type speedgoat_target_model_cal_impl;
extern speedgoat_target_model_cal_type *speedgoat_target_model_2021_cal;

#endif                      /* RTW_HEADER_speedgoat_target_model_2021b_cal_h_ */
