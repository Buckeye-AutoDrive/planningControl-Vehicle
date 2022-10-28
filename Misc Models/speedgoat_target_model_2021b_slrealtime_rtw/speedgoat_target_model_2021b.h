/*
 * speedgoat_target_model_2021b.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "speedgoat_target_model_2021b".
 *
 * Model version              : 1.14
 * Simulink Coder version : 9.6 (R2021b) 14-May-2021
 * C++ source code generated on : Tue Oct 25 18:51:47 2022
 *
 * Target selection: slrealtime.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_speedgoat_target_model_2021b_h_
#define RTW_HEADER_speedgoat_target_model_2021b_h_
#include <cmath>
#include <cstring>
#include <cstring>
#include "rtwtypes.h"
#include "rtw_extmode.h"
#include "sysran_types.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "verify/verifyIntrf.h"
#include "slrealtime/libsrc/IP/udp.hpp"
#include "sf_runtime/slTestResult.h"
#include "sf_runtime/slTestTypes.h"
#include "speedgoat_target_model_2021b_types.h"

/* Shared type includes */
#include "multiword_types.h"

/* Child system includes */
#include "speedgoat_target_model_2021b_cal.h"
#include "rtGetNaN.h"
#include "rt_nonfinite.h"
#include "crl_mutex.hpp"
#include "rt_assert.h"
#include "slrealtime_cblas.hpp"
#include "rtGetInf.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmCounterLimit
#define rtmCounterLimit(rtm, idx)      ((rtm)->Timing.TaskCounters.cLimit[(idx)])
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmStepTask
#define rtmStepTask(rtm, idx)          ((rtm)->Timing.TaskCounters.TID[(idx)] == 0)
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   ((rtm)->Timing.taskTime0)
#endif

#ifndef rtmTaskCounter
#define rtmTaskCounter(rtm, idx)       ((rtm)->Timing.TaskCounters.TID[(idx)])
#endif

/* Block signals for system '<S4>/DataTypeConversion_L0' */
struct B_DataTypeConversion_L0_speed_T {
  real_T y;                            /* '<S4>/DataTypeConversion_L0' */
};

/* Block signals for system '<S4>/DataTypeConversion_amax' */
struct B_DataTypeConversion_amax_spe_T {
  real_T y;                            /* '<S4>/DataTypeConversion_amax' */
};

/* Block signals (default storage) */
struct B_speedgoat_target_model_2021_T {
  real_T Product;                      /* '<S47>/Product' */
  real_T Gain;                         /* '<S48>/Gain' */
  real_T Multiply;                     /* '<S47>/Multiply' */
  real_T CurvedSteadyYaw;              /* '<S47>/Gain1' */
  real_T Gain_d;                       /* '<S1>/Gain' */
  real_T Minus;                        /* '<S47>/Minus' */
  real_T YawRateFeedback;              /* '<S47>/Gain' */
  real_T UnitDelay;                    /* '<S47>/Unit Delay' */
  real_T UnitDelay_n;                  /* '<S1>/Unit Delay' */
  real_T Minus1;                       /* '<S47>/Minus1' */
  real_T SteeringDelay;                /* '<S47>/Gain2' */
  real_T Add;                          /* '<S47>/Add' */
  real_T Saturation;                   /* '<S47>/Saturation' */
  real_T ApproachingSpeed;             /* '<S1>/Approaching Speed' */
  real_T SceneSpeed;                   /* '<S1>/Scene Speed' */
  real_T Constant1;                    /* '<S1>/Constant1' */
  real_T RRCrossing;                   /* '<S1>/RR Crossing' */
  real_T NavigatingFlag;               /* '<S1>/Navigating Flag' */
  real_T Clear2Go;                     /* '<S1>/Clear2Go' */
  real_T RoadSpeed;                    /* '<S1>/RoadSpeed' */
  real_T WaitTime;                     /* '<S1>/WaitTime' */
  real_T MinSpeed;                     /* '<S1>/MinSpeed' */
  real_T LaneAvailable;                /* '<S1>/LaneAvailable' */
  real_T LaneChanged;                  /* '<S1>/LaneChanged' */
  real_T acceleration;
  real_T TmpRTBAtSwitchOutport1;       /* '<S49>/Switch' */
  real_T TmpRTBAtSwitch1Outport1;      /* '<S49>/Switch1' */
  real_T Merge1;                       /* '<S11>/Merge1' */
  real_T Merge2;                       /* '<S11>/Merge2' */
  real_T Product2;                     /* '<S4>/Product2' */
  real_T safe_distance;                /* '<S4>/Sum1' */
  real_T lead_velocity;                /* '<S4>/Sum6' */
  real_T Floor;                        /* '<S23>/Floor' */
  real_T Floor1;                       /* '<S23>/Floor1' */
  real_T ym_zero[2];                   /* '<S23>/ym_zero' */
  real_T default_spacing;              /* '<S4>/Default spacing constant' */
  real_T MathFunction[2];              /* '<S23>/Math Function' */
  real_T MathFunction1;                /* '<S23>/Math Function1' */
  real_T MathFunction2;                /* '<S23>/Math Function2' */
  real_T umin_scale4;                  /* '<S23>/umin_scale4' */
  real_T Reshape;                      /* '<S23>/Reshape' */
  real_T ymin_scale1[2];               /* '<S23>/ymin_scale1' */
  real_T Reshape1[2];                  /* '<S23>/Reshape1' */
  real_T ymin_scale2;                  /* '<S23>/ymin_scale2' */
  real_T Reshape2;                     /* '<S23>/Reshape2' */
  real_T Reshape3[2];                  /* '<S23>/Reshape3' */
  real_T Reshape4;                     /* '<S23>/Reshape4' */
  real_T Reshape5;                     /* '<S23>/Reshape5' */
  real_T extmv_scale;                  /* '<S23>/ext.mv_scale' */
  real_T extmv_scale1;                 /* '<S23>/ext.mv_scale1' */
  real_T last_mv;                      /* '<S23>/last_mv' */
  real_T last_x[4];                    /* '<S23>/last_x' */
  real_T TmpRTBAtoptimizerInport3[2];
  real_T TmpRTBAtoptimizerInport5;
  real_T min_velocity[2];
  real_T TmpRTBAtDataTypeConversion_vset;/* '<S4>/DataTypeConversion_vset' */
  real_T umin_scale1;                  /* '<S23>/umin_scale1' */
  real_T umin_scale3[11];              /* '<S23>/umin_scale3' */
  real_T umin_scale5[22];              /* '<S23>/umin_scale5' */
  real_T umin_scale2;                  /* '<S23>/umin_scale2' */
  real_T TmpRTBAtEqual2Inport1;
  real_T set_velocity1;
  real_T Sign2;                        /* '<S51>/Sign2' */
  real_T set_velocity;
  real_T TmpRTBAtEqual1Inport1;
  real_T TmpRTBAtSign1Inport1;
  real_T Sign1;                        /* '<S51>/Sign1' */
  real_T TmpRTBAtEqual3Inport1;
  real_T TmpRTBAtEqual5Inport2;
  real_T TmpRTBAtEqual4Inport2;
  real_T TmpRTBAtMultiplyInport1;
  real_T TmpRTBAtSwitchCaseInport1;
  real_T set_velocity_c;
  real_T TmpRTBAtMinusInport2;
  real_T error;                        /* '<S8>/Minus' */
  real_T Merge;                        /* '<S50>/Merge' */
  real_T Multiply_k;                   /* '<S49>/Multiply' */
  real_T Switch;                       /* '<S49>/Switch' */
  real_T Switch1;                      /* '<S49>/Switch1' */
  real_T UDPReceive2_o2;               /* '<S174>/UDP Receive2' */
  real_T Subtract;                     /* '<S174>/Subtract' */
  real_T UDPReceive2_o2_o;             /* '<S173>/UDP Receive2' */
  real_T Subtract_p;                   /* '<S173>/Subtract' */
  real_T UDPReceive2_o2_h;             /* '<S172>/UDP Receive2' */
  real_T Subtract_d;                   /* '<S172>/Subtract' */
  real_T UDPReceive2_o2_j;             /* '<S171>/UDP Receive2' */
  real_T Subtract_b;                   /* '<S171>/Subtract' */
  real_T UDPReceive2_o2_c;             /* '<S170>/UDP Receive2' */
  real_T Subtract_bo;                  /* '<S170>/Subtract' */
  real_T Direction;                    /* '<S169>/Byte Unpacking ' */
  real_T Curvature;                    /* '<S169>/Byte Unpacking ' */
  real_T yawrate;                      /* '<S169>/Byte Unpacking ' */
  real_T ObjectDist;                   /* '<S168>/Byte Unpacking ' */
  real_T D2S;                          /* '<S168>/Byte Unpacking ' */
  real_T Relativedistance;             /* '<S168>/Byte Unpacking ' */
  real_T TLstate[90];                  /* '<S168>/Byte Unpacking ' */
  real_T RefPoses[3];                  /* '<S167>/Byte Unpacking ' */
  real_T ObjAhead;                     /* '<S166>/Byte Unpacking ' */
  real_T Longitudinalvelocity;         /* '<S166>/Byte Unpacking ' */
  real_T SceneID;                      /* '<S166>/Byte Unpacking ' */
  real_T Relativevelocity;             /* '<S166>/Byte Unpacking ' */
  real_T CurrPoses[3];                 /* '<S165>/Byte Unpacking 1' */
  real_T v;                            /* '<S160>/Discrete-Time Integrator' */
  real_T Square;                       /* '<S160>/Square' */
  real_T F_aero;                       /* '<S160>/Gain1' */
  real_T Gain2;                        /* '<S160>/Gain2' */
  real_T Gain_dx;                      /* '<S160>/Gain' */
  real_T Sum;                          /* '<S160>/Sum' */
  real_T T_eng;                        /* '<S160>/Gain4' */
  real_T BrakeSaturation;              /* '<S156>/Brake Saturation' */
  real_T v_g;                          /* '<S159>/Discrete-Time Integrator' */
  real_T Gain_h;                       /* '<S159>/Gain' */
  real_T Square_i;                     /* '<S159>/Square' */
  real_T F_aero_k;                     /* '<S159>/Gain1' */
  real_T Gain2_k;                      /* '<S159>/Gain2' */
  real_T Sum_o;                        /* '<S159>/Sum' */
  real_T T_eng_h;                      /* '<S159>/Gain4' */
  real_T ThrottleSaturation;           /* '<S156>/Throttle Saturation' */
  real_T v_p;                          /* '<S164>/Discrete-Time Integrator' */
  real_T Square_j;                     /* '<S164>/Square' */
  real_T F_aero_h;                     /* '<S164>/Gain1' */
  real_T Gain2_p;                      /* '<S164>/Gain2' */
  real_T Gain_g;                       /* '<S164>/Gain' */
  real_T Sum_d;                        /* '<S164>/Sum' */
  real_T T_eng_c;                      /* '<S164>/Gain4' */
  real_T BrakeSaturation_n;            /* '<S157>/Brake Saturation' */
  real_T v_o;                          /* '<S163>/Discrete-Time Integrator' */
  real_T Gain_b;                       /* '<S163>/Gain' */
  real_T Square_k;                     /* '<S163>/Square' */
  real_T F_aero_hy;                    /* '<S163>/Gain1' */
  real_T Gain2_b;                      /* '<S163>/Gain2' */
  real_T Sum_f;                        /* '<S163>/Sum' */
  real_T T_eng_c0;                     /* '<S163>/Gain4' */
  real_T ThrottleSaturation_m;         /* '<S157>/Throttle Saturation' */
  real_T accel;                        /* '<S157>/MATLAB Function' */
  real_T decel;                        /* '<S157>/MATLAB Function' */
  real_T RefVelCmd;                    /* '<S1>/MATLAB Function2' */
  real_T D2D;                          /* '<S1>/MATLAB Function' */
  real_T ProportionalGain;             /* '<S144>/Proportional Gain' */
  real_T Integrator;                   /* '<S139>/Integrator' */
  real_T Sum_fg;                       /* '<S148>/Sum' */
  real_T ZeroGain;                     /* '<S130>/ZeroGain' */
  real_T DeadZone;                     /* '<S132>/DeadZone' */
  real_T SignPreSat;                   /* '<S130>/SignPreSat' */
  real_T IntegralGain;                 /* '<S136>/Integral Gain' */
  real_T SignPreIntegrator;            /* '<S130>/SignPreIntegrator' */
  real_T Switch_o;                     /* '<S130>/Switch' */
  real_T ProportionalGain_f;           /* '<S93>/Proportional Gain' */
  real_T Integrator_b;                 /* '<S88>/Integrator' */
  real_T Sum_k;                        /* '<S97>/Sum' */
  real_T ZeroGain_l;                   /* '<S79>/ZeroGain' */
  real_T DeadZone_n;                   /* '<S81>/DeadZone' */
  real_T SignPreSat_m;                 /* '<S79>/SignPreSat' */
  real_T IntegralGain_a;               /* '<S85>/Integral Gain' */
  real_T SignPreIntegrator_d;          /* '<S79>/SignPreIntegrator' */
  real_T Switch_l;                     /* '<S79>/Switch' */
  real_T Abs;                          /* '<S49>/Abs' */
  real_T steerCmd;                     /* '<S7>/Kinematic' */
  real_T ACC;                          /* '<S1>/Collision Avoidance' */
  real_T LaneChangeCmd;                /* '<S1>/Collision Avoidance' */
  real_T VelCmd;                       /* '<S1>/Chart' */
  real_T ACC_i;                        /* '<S1>/Chart' */
  real_T Tgap;                         /* '<S1>/Chart' */
  real_T Trans;                        /* '<S1>/Chart' */
  real_T StopSim;                      /* '<S1>/Chart' */
  real_T TmpSignalConversionAtSFunctionI[2];/* '<S43>/optimizer' */
  real_T TmpSignalConversionAtSFunctio_m[2];/* '<S43>/optimizer' */
  real_T xk1[4];                       /* '<S43>/optimizer' */
  real_T u;                            /* '<S43>/optimizer' */
  real_T cost;                         /* '<S43>/optimizer' */
  real_T useq[11];                     /* '<S43>/optimizer' */
  real_T xseq[44];                     /* '<S43>/optimizer' */
  real_T yseq[22];                     /* '<S43>/optimizer' */
  real_T status;                       /* '<S43>/optimizer' */
  real_T xest[4];                      /* '<S43>/optimizer' */
  real_T y;                            /* '<S4>/DataTypeConversion_optsgn' */
  int32_T WhileIterator;               /* '<S174>/While Iterator' */
  int32_T WhileIterator_k;             /* '<S173>/While Iterator' */
  int32_T WhileIterator_c;             /* '<S172>/While Iterator' */
  int32_T WhileIterator_p;             /* '<S171>/While Iterator' */
  int32_T WhileIterator_e;             /* '<S170>/While Iterator' */
  int16_T OR;                          /* '<S1>/OR' */
  uint8_T BytePacking[24];             /* '<S181>/Byte Packing' */
  uint8_T BytePacking1[16];            /* '<S180>/Byte Packing1' */
  uint8_T UDPReceive2_o1[24];          /* '<S174>/UDP Receive2' */
  uint8_T UDPReceive2_o1_o[744];       /* '<S173>/UDP Receive2' */
  uint8_T UDPReceive2_o1_d[24];        /* '<S172>/UDP Receive2' */
  uint8_T UDPReceive2_o1_p[32];        /* '<S171>/UDP Receive2' */
  uint8_T UDPReceive2_o1_m[24];        /* '<S170>/UDP Receive2' */
  int8_T DataTypeConv1;                /* '<S130>/DataTypeConv1' */
  int8_T DataTypeConv2;                /* '<S130>/DataTypeConv2' */
  int8_T DataTypeConv1_p;              /* '<S79>/DataTypeConv1' */
  int8_T DataTypeConv2_l;              /* '<S79>/DataTypeConv2' */
  boolean_T DataTypeConversion4;       /* '<S2>/Data Type Conversion4' */
  boolean_T DataTypeConversion5;       /* '<S2>/Data Type Conversion5' */
  boolean_T DataTypeConversion2;       /* '<S2>/Data Type Conversion2' */
  boolean_T DataTypeConversion3;       /* '<S2>/Data Type Conversion3' */
  boolean_T DataTypeConversion1;       /* '<S2>/Data Type Conversion1' */
  boolean_T Memory[34];                /* '<S23>/Memory' */
  boolean_T Equal2;                    /* '<S51>/Equal2' */
  boolean_T Equal6;                    /* '<S51>/Equal6' */
  boolean_T AND3;                      /* '<S51>/AND3' */
  boolean_T Equal1;                    /* '<S51>/Equal1' */
  boolean_T Equal3;                    /* '<S51>/Equal3' */
  boolean_T AND1;                      /* '<S51>/AND1' */
  boolean_T OR_h;                      /* '<S51>/OR' */
  boolean_T NOT1;                      /* '<S51>/NOT1' */
  boolean_T Equal5;                    /* '<S51>/Equal5' */
  boolean_T Equal4;                    /* '<S51>/Equal4' */
  boolean_T AND2;                      /* '<S51>/AND2' */
  boolean_T NOT;                       /* '<S51>/NOT' */
  boolean_T DataTypeConversion;        /* '<S50>/Data Type Conversion' */
  boolean_T Compare;                   /* '<S179>/Compare' */
  boolean_T Compare_j;                 /* '<S178>/Compare' */
  boolean_T Compare_n;                 /* '<S177>/Compare' */
  boolean_T Compare_p;                 /* '<S176>/Compare' */
  boolean_T Compare_b;                 /* '<S175>/Compare' */
  boolean_T NotEqual;                  /* '<S130>/NotEqual' */
  boolean_T Equal1_c;                  /* '<S130>/Equal1' */
  boolean_T AND3_b;                    /* '<S130>/AND3' */
  boolean_T NotEqual_b;                /* '<S79>/NotEqual' */
  boolean_T Equal1_a;                  /* '<S79>/Equal1' */
  boolean_T AND3_j;                    /* '<S79>/AND3' */
  boolean_T iAout[34];                 /* '<S43>/optimizer' */
  B_DataTypeConversion_L0_speed_T sf_DataTypeConversion_vset;/* '<S4>/DataTypeConversion_vset' */
  B_DataTypeConversion_L0_speed_T sf_DataTypeConversion_vlead;/* '<S4>/DataTypeConversion_vlead' */
  B_DataTypeConversion_L0_speed_T sf_DataTypeConversion_vego;/* '<S4>/DataTypeConversion_vego' */
  B_DataTypeConversion_L0_speed_T sf_DataTypeConversion_reldist;/* '<S4>/DataTypeConversion_reldist' */
  B_DataTypeConversion_L0_speed_T sf_DataTypeConversion_dmin;/* '<S4>/DataTypeConversion_dmin' */
  B_DataTypeConversion_amax_spe_T sf_DataTypeConversion_atrack;/* '<S4>/DataTypeConversion_atrack' */
  B_DataTypeConversion_amax_spe_T sf_DataTypeConversion_amin;/* '<S4>/DataTypeConversion_amin' */
  B_DataTypeConversion_amax_spe_T sf_DataTypeConversion_amax;/* '<S4>/DataTypeConversion_amax' */
  B_DataTypeConversion_L0_speed_T sf_DataTypeConversion_L0;/* '<S4>/DataTypeConversion_L0' */
};

/* Block states (default storage) for system '<Root>' */
struct DW_speedgoat_target_model_202_T {
  slTestBlkInfo_t Assertion_sltestBlkInfo;/* '<S51>/Assertion' */
  slTestBlkInfo_t Assertion1_sltestBlkInfo;/* '<S51>/Assertion1' */
  real_T UnitDelay_DSTATE;             /* '<S47>/Unit Delay' */
  real_T UnitDelay_DSTATE_o;           /* '<S1>/Unit Delay' */
  real_T last_mv_DSTATE;               /* '<S23>/last_mv' */
  real_T DiscreteTimeIntegrator_DSTATE;/* '<S160>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator_DSTATE_a;/* '<S159>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator_DSTATE_f;/* '<S164>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator_DSTATE_e;/* '<S163>/Discrete-Time Integrator' */
  real_T Integrator_DSTATE;            /* '<S139>/Integrator' */
  real_T Integrator_DSTATE_k;          /* '<S88>/Integrator' */
  real_T acceleration_Buf[2];          /* synthesized block */
  real_T TmpRTBAtSwitchOutport1_Buf[2];/* synthesized block */
  real_T TmpRTBAtSwitch1Outport1_Buf[2];/* synthesized block */
  real_T TmpRTBAtoptimizerInport3_Buf[4];/* synthesized block */
  real_T TmpRTBAtoptimizerInport5_Buf[2];/* synthesized block */
  real_T min_velocity_Buf0[2];         /* synthesized block */
  real_T min_velocity_Buf1[2];         /* synthesized block */
  real_T min_velocity_Buf2[2];         /* synthesized block */
  real_T TmpRTBAtDataTypeConversion_vset[2];/* synthesized block */
  real_T default_spacing_Buf0;         /* synthesized block */
  real_T default_spacing_Buf1;         /* synthesized block */
  real_T default_spacing_Buf2;         /* synthesized block */
  real_T TmpRTBAtEqual1Inport1_Buf[2]; /* synthesized block */
  real_T TmpRTBAtEqual2Inport1_Buf[2]; /* synthesized block */
  real_T TmpRTBAtEqual3Inport1_Buf[2]; /* synthesized block */
  real_T TmpRTBAtEqual4Inport2_Buf[2]; /* synthesized block */
  real_T TmpRTBAtEqual5Inport2_Buf[2]; /* synthesized block */
  real_T TmpRTBAtSign1Inport1_Buf[2];  /* synthesized block */
  real_T set_velocity_Buf[2];          /* synthesized block */
  real_T set_velocity1_Buf[2];         /* synthesized block */
  real_T TmpRTBAtMultiplyInport1_Buf[2];/* synthesized block */
  real_T TmpRTBAtSwitchCaseInport1_Buf[2];/* synthesized block */
  real_T TmpRTBAtMinusInport2_Buf[2];  /* synthesized block */
  real_T set_velocity_Buf_n[2];        /* synthesized block */
  real_T last_x_PreviousInput[4];      /* '<S23>/last_x' */
  real_T Assertion_sltestLastResultTime;/* '<S51>/Assertion' */
  real_T Assertion1_sltestLastResultTime;/* '<S51>/Assertion1' */
  real_T Subtract_DWORK1;              /* '<S174>/Subtract' */
  real_T Stuck;                        /* '<S1>/Collision Avoidance' */
  real_T DestReached;                  /* '<S1>/Chart' */
  real_T ExitFlag;                     /* '<S1>/Chart' */
  real_T CurrentVelCmd;                /* '<S1>/Chart' */
  real_T Navi;                         /* '<S1>/Chart' */
  void *UDPSendControllerOutput_PWORK; /* '<S3>/UDP Send: Controller Output' */
  void *UDPSendControllerIntermediateSi;
                          /* '<S3>/UDP Send: Controller Intermediate Signals' */
  void* min_velocity_d0_SEMAPHORE;     /* synthesized block */
  void* default_spacing_d0_SEMAPHORE;  /* synthesized block */
  void *UDPReceive2_PWORK[2];          /* '<S174>/UDP Receive2' */
  void *UDPReceive2_PWORK_b[2];        /* '<S173>/UDP Receive2' */
  void *UDPReceive2_PWORK_f[2];        /* '<S172>/UDP Receive2' */
  void *UDPReceive2_PWORK_e[2];        /* '<S171>/UDP Receive2' */
  void *UDPReceive2_PWORK_e0[2];       /* '<S170>/UDP Receive2' */
  int32_T Assertion_sltestFinalResult; /* '<S51>/Assertion' */
  int32_T Assertion_sltestCurrentResult;/* '<S51>/Assertion' */
  int32_T Assertion1_sltestFinalResult;/* '<S51>/Assertion1' */
  int32_T Assertion1_sltestCurrentResult;/* '<S51>/Assertion1' */
  int32_T sfEvent;                     /* '<S1>/Collision Avoidance' */
  int32_T sfEvent_o;                   /* '<S1>/Chart' */
  uint32_T is_c9_speedgoat_target_model_20;/* '<S1>/Collision Avoidance' */
  uint32_T temporalCounter_i1;         /* '<S1>/Collision Avoidance' */
  uint32_T is_c21_speedgoat_target_model_2;/* '<S1>/Chart' */
  uint32_T is_Navigation;              /* '<S1>/Chart' */
  uint32_T is_RailCrossing;            /* '<S1>/Chart' */
  uint32_T is_RoundAbout;              /* '<S1>/Chart' */
  uint32_T is_TrafficLight;            /* '<S1>/Chart' */
  uint32_T is_StopSign;                /* '<S1>/Chart' */
  int_T BytePacking_IWORK[6];          /* '<S181>/Byte Packing' */
  int_T UDPSendControllerOutput_IWORK[2];/* '<S3>/UDP Send: Controller Output' */
  int_T BytePacking1_IWORK[4];         /* '<S180>/Byte Packing1' */
  int_T UDPSendControllerIntermediate_a[2];
                          /* '<S3>/UDP Send: Controller Intermediate Signals' */
  int_T UDPReceive2_IWORK;             /* '<S174>/UDP Receive2' */
  int_T UDPReceive2_IWORK_l;           /* '<S173>/UDP Receive2' */
  int_T UDPReceive2_IWORK_lp;          /* '<S172>/UDP Receive2' */
  int_T UDPReceive2_IWORK_h;           /* '<S171>/UDP Receive2' */
  int_T UDPReceive2_IWORK_a;           /* '<S170>/UDP Receive2' */
  int_T ByteUnpacking_IWORK[6];        /* '<S169>/Byte Unpacking ' */
  int_T ByteUnpacking_IWORK_l[8];      /* '<S168>/Byte Unpacking ' */
  int_T ByteUnpacking_IWORK_i[2];      /* '<S167>/Byte Unpacking ' */
  int_T ByteUnpacking_IWORK_d[8];      /* '<S166>/Byte Unpacking ' */
  int_T ByteUnpacking1_IWORK[2];       /* '<S165>/Byte Unpacking 1' */
  int8_T acceleration_RdBufIdx;        /* synthesized block */
  int8_T acceleration_WrBufIdx;        /* synthesized block */
  int8_T TmpRTBAtSwitchOutport1_RdBufIdx;/* synthesized block */
  int8_T TmpRTBAtSwitchOutport1_WrBufIdx;/* synthesized block */
  int8_T TmpRTBAtSwitch1Outport1_RdBufId;/* synthesized block */
  int8_T TmpRTBAtSwitch1Outport1_WrBufId;/* synthesized block */
  int8_T TmpRTBAtoptimizerInport3_RdBufI;/* synthesized block */
  int8_T TmpRTBAtoptimizerInport3_WrBufI;/* synthesized block */
  int8_T TmpRTBAtoptimizerInport5_RdBufI;/* synthesized block */
  int8_T TmpRTBAtoptimizerInport5_WrBufI;/* synthesized block */
  int8_T min_velocity_LstBufWR;        /* synthesized block */
  int8_T min_velocity_RDBuf;           /* synthesized block */
  int8_T TmpRTBAtDataTypeConversion_vs_n;/* synthesized block */
  int8_T TmpRTBAtDataTypeConversion_vs_m;/* synthesized block */
  int8_T default_spacing_LstBufWR;     /* synthesized block */
  int8_T default_spacing_RDBuf;        /* synthesized block */
  int8_T TmpRTBAtEqual1Inport1_RdBufIdx;/* synthesized block */
  int8_T TmpRTBAtEqual1Inport1_WrBufIdx;/* synthesized block */
  int8_T TmpRTBAtEqual2Inport1_RdBufIdx;/* synthesized block */
  int8_T TmpRTBAtEqual2Inport1_WrBufIdx;/* synthesized block */
  int8_T TmpRTBAtEqual3Inport1_RdBufIdx;/* synthesized block */
  int8_T TmpRTBAtEqual3Inport1_WrBufIdx;/* synthesized block */
  int8_T TmpRTBAtEqual4Inport2_RdBufIdx;/* synthesized block */
  int8_T TmpRTBAtEqual4Inport2_WrBufIdx;/* synthesized block */
  int8_T TmpRTBAtEqual5Inport2_RdBufIdx;/* synthesized block */
  int8_T TmpRTBAtEqual5Inport2_WrBufIdx;/* synthesized block */
  int8_T TmpRTBAtSign1Inport1_RdBufIdx;/* synthesized block */
  int8_T TmpRTBAtSign1Inport1_WrBufIdx;/* synthesized block */
  int8_T set_velocity_RdBufIdx;        /* synthesized block */
  int8_T set_velocity_WrBufIdx;        /* synthesized block */
  int8_T set_velocity1_RdBufIdx;       /* synthesized block */
  int8_T set_velocity1_WrBufIdx;       /* synthesized block */
  int8_T TmpRTBAtMultiplyInport1_RdBufId;/* synthesized block */
  int8_T TmpRTBAtMultiplyInport1_WrBufId;/* synthesized block */
  int8_T TmpRTBAtSwitchCaseInport1_RdBuf;/* synthesized block */
  int8_T TmpRTBAtSwitchCaseInport1_WrBuf;/* synthesized block */
  int8_T TmpRTBAtMinusInport2_RdBufIdx;/* synthesized block */
  int8_T TmpRTBAtMinusInport2_WrBufIdx;/* synthesized block */
  int8_T set_velocity_RdBufIdx_p;      /* synthesized block */
  int8_T set_velocity_WrBufIdx_k;      /* synthesized block */
  int8_T ByteUnpackThirdset_SubsysRanBC;/* '<S2>/Byte Unpack: Third set' */
  int8_T ByteUnpackSecondset_SubsysRanBC;/* '<S2>/Byte Unpack: Second set' */
  int8_T ByteUnpackFourthset_SubsysRanBC;/* '<S2>/Byte Unpack: Fourth set' */
  int8_T ByteUnpackFirstset_SubsysRanBC;/* '<S2>/Byte Unpack: First set' */
  int8_T ByteUnpackFifthset_SubsysRanBC;/* '<S2>/Byte Unpack: Fifth set' */
  int8_T Subsystem1_SubsysRanBC;       /* '<S11>/Subsystem1' */
  int8_T DiscreteTimeIntegrator_PrevRese;/* '<S160>/Discrete-Time Integrator' */
  int8_T DiscreteTimeIntegrator_PrevRe_a;/* '<S159>/Discrete-Time Integrator' */
  int8_T Subsystem2_SubsysRanBC;       /* '<S11>/Subsystem2' */
  int8_T DiscreteTimeIntegrator_PrevRe_b;/* '<S164>/Discrete-Time Integrator' */
  int8_T DiscreteTimeIntegrator_PrevR_bg;/* '<S163>/Discrete-Time Integrator' */
  int8_T Reverse_SubsysRanBC;          /* '<S50>/Reverse' */
  int8_T Integrator_PrevResetState;    /* '<S139>/Integrator' */
  int8_T Forward_SubsysRanBC;          /* '<S50>/Forward' */
  int8_T Integrator_PrevResetState_j;  /* '<S88>/Integrator' */
  uint8_T is_active_c9_speedgoat_target_m;/* '<S1>/Collision Avoidance' */
  uint8_T is_active_c21_speedgoat_target_;/* '<S1>/Chart' */
  uint8_T temporalCounter_i1_i;        /* '<S1>/Chart' */
  boolean_T Memory_PreviousInput[34];  /* '<S23>/Memory' */
};

/* Invariant block signals (default storage) */
struct ConstB_speedgoat_target_model_T {
  real_T Width;                        /* '<S3>/Width' */
  real_T Width1;                       /* '<S3>/Width1' */
};

/* Real-time Model Data Structure */
struct tag_RTM_speedgoat_target_mode_T {
  const char_T *errorStatus;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    time_T taskTime0;
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    struct {
      uint8_T TID[2];
      uint8_T cLimit[2];
    } TaskCounters;

    struct {
      uint8_T TID0_1;
    } RateInteraction;

    boolean_T stopRequestedFlag;
  } Timing;
};

/* Block signals (default storage) */
#ifdef __cplusplus

extern "C" {

#endif

  extern struct B_speedgoat_target_model_2021_T speedgoat_target_model_2021b_B;

#ifdef __cplusplus

}
#endif

/* Block states (default storage) */
extern struct DW_speedgoat_target_model_202_T speedgoat_target_model_2021b_DW;
extern const ConstB_speedgoat_target_model_T speedgoat_target_model_2_ConstB;/* constant block i/o */
extern void rate_scheduler(void);

#ifdef __cplusplus

extern "C" {

#endif

  /* Model entry point functions */
  extern void speedgoat_target_model_2021b_initialize(void);
  extern void speedgoat_target_model_2021b_step0(void);
  extern void speedgoat_target_model_2021b_step1(void);
  extern void speedgoat_target_model_2021b_terminate(void);

#ifdef __cplusplus

}
#endif

/* Real-time Model object */
#ifdef __cplusplus

extern "C" {

#endif

  extern RT_MODEL_speedgoat_target_mod_T *const speedgoat_target_model_2021b_M;

#ifdef __cplusplus

}
#endif

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'speedgoat_target_model_2021b'
 * '<S1>'   : 'speedgoat_target_model_2021b/Planning and Control'
 * '<S2>'   : 'speedgoat_target_model_2021b/UDP Data Receive'
 * '<S3>'   : 'speedgoat_target_model_2021b/UDP Data Send'
 * '<S4>'   : 'speedgoat_target_model_2021b/Planning and Control/ACC system'
 * '<S5>'   : 'speedgoat_target_model_2021b/Planning and Control/Chart'
 * '<S6>'   : 'speedgoat_target_model_2021b/Planning and Control/Collision Avoidance'
 * '<S7>'   : 'speedgoat_target_model_2021b/Planning and Control/Lateral Controller Stanley'
 * '<S8>'   : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley'
 * '<S9>'   : 'speedgoat_target_model_2021b/Planning and Control/MATLAB Function'
 * '<S10>'  : 'speedgoat_target_model_2021b/Planning and Control/MATLAB Function2'
 * '<S11>'  : 'speedgoat_target_model_2021b/Planning and Control/Subsystem3'
 * '<S12>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/DataTypeConversion_L0'
 * '<S13>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/DataTypeConversion_amax'
 * '<S14>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/DataTypeConversion_amin'
 * '<S15>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/DataTypeConversion_atrack'
 * '<S16>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/DataTypeConversion_dmin'
 * '<S17>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/DataTypeConversion_optsgn'
 * '<S18>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/DataTypeConversion_reldist'
 * '<S19>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/DataTypeConversion_vego'
 * '<S20>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/DataTypeConversion_vlead'
 * '<S21>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/DataTypeConversion_vset'
 * '<S22>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC'
 * '<S23>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC'
 * '<S24>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/MPC Matrix Signal Check'
 * '<S25>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/MPC Matrix Signal Check1'
 * '<S26>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/MPC Matrix Signal Check2'
 * '<S27>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/MPC Preview Signal Check'
 * '<S28>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/MPC Preview Signal Check1'
 * '<S29>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/MPC Preview Signal Check2'
 * '<S30>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/MPC Preview Signal Check3'
 * '<S31>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/MPC Preview Signal Check4'
 * '<S32>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/MPC Preview Signal Check5'
 * '<S33>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/MPC Preview Signal Check6'
 * '<S34>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/MPC Preview Signal Check7'
 * '<S35>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/MPC Preview Signal Check8'
 * '<S36>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/MPC Scalar Signal Check'
 * '<S37>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/MPC Scalar Signal Check1'
 * '<S38>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/MPC Scalar Signal Check2'
 * '<S39>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/MPC Vector Signal Check'
 * '<S40>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/MPC Vector Signal Check1'
 * '<S41>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/MPC Vector Signal Check6'
 * '<S42>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/moorx'
 * '<S43>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/optimizer'
 * '<S44>'  : 'speedgoat_target_model_2021b/Planning and Control/ACC system/MPC/MPC/optimizer/optimizer'
 * '<S45>'  : 'speedgoat_target_model_2021b/Planning and Control/Lateral Controller Stanley/Dynamic'
 * '<S46>'  : 'speedgoat_target_model_2021b/Planning and Control/Lateral Controller Stanley/Kinematic'
 * '<S47>'  : 'speedgoat_target_model_2021b/Planning and Control/Lateral Controller Stanley/Dynamic/Dynamic Enabled'
 * '<S48>'  : 'speedgoat_target_model_2021b/Planning and Control/Lateral Controller Stanley/Dynamic/Dynamic Enabled/Radians to Degrees'
 * '<S49>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/Convert Command '
 * '<S50>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller'
 * '<S51>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/Verify Direction'
 * '<S52>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward'
 * '<S53>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse'
 * '<S54>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward '
 * '<S55>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Anti-windup'
 * '<S56>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /D Gain'
 * '<S57>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Filter'
 * '<S58>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Filter ICs'
 * '<S59>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /I Gain'
 * '<S60>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Ideal P Gain'
 * '<S61>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Ideal P Gain Fdbk'
 * '<S62>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Integrator'
 * '<S63>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Integrator ICs'
 * '<S64>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /N Copy'
 * '<S65>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /N Gain'
 * '<S66>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /P Copy'
 * '<S67>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Parallel P Gain'
 * '<S68>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Reset Signal'
 * '<S69>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Saturation'
 * '<S70>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Saturation Fdbk'
 * '<S71>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Sum'
 * '<S72>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Sum Fdbk'
 * '<S73>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Tracking Mode'
 * '<S74>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Tracking Mode Sum'
 * '<S75>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Tsamp - Integral'
 * '<S76>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Tsamp - Ngain'
 * '<S77>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /postSat Signal'
 * '<S78>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /preSat Signal'
 * '<S79>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Anti-windup/Disc. Clamping Parallel'
 * '<S80>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Anti-windup/Disc. Clamping Parallel/Dead Zone'
 * '<S81>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Anti-windup/Disc. Clamping Parallel/Dead Zone/Enabled'
 * '<S82>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /D Gain/Disabled'
 * '<S83>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Filter/Disabled'
 * '<S84>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Filter ICs/Disabled'
 * '<S85>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /I Gain/Internal Parameters'
 * '<S86>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Ideal P Gain/Passthrough'
 * '<S87>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Ideal P Gain Fdbk/Disabled'
 * '<S88>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Integrator/Discrete'
 * '<S89>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Integrator ICs/Internal IC'
 * '<S90>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /N Copy/Disabled wSignal Specification'
 * '<S91>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /N Gain/Disabled'
 * '<S92>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /P Copy/Disabled'
 * '<S93>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Parallel P Gain/Internal Parameters'
 * '<S94>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Reset Signal/External Reset'
 * '<S95>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Saturation/Enabled'
 * '<S96>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Saturation Fdbk/Disabled'
 * '<S97>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Sum/Sum_PI'
 * '<S98>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Sum Fdbk/Disabled'
 * '<S99>'  : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Tracking Mode/Disabled'
 * '<S100>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Tracking Mode Sum/Passthrough'
 * '<S101>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Tsamp - Integral/Passthrough'
 * '<S102>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Tsamp - Ngain/Passthrough'
 * '<S103>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /postSat Signal/Forward_Path'
 * '<S104>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /preSat Signal/Forward_Path'
 * '<S105>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse'
 * '<S106>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Anti-windup'
 * '<S107>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/D Gain'
 * '<S108>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Filter'
 * '<S109>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Filter ICs'
 * '<S110>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/I Gain'
 * '<S111>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Ideal P Gain'
 * '<S112>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Ideal P Gain Fdbk'
 * '<S113>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Integrator'
 * '<S114>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Integrator ICs'
 * '<S115>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/N Copy'
 * '<S116>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/N Gain'
 * '<S117>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/P Copy'
 * '<S118>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Parallel P Gain'
 * '<S119>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Reset Signal'
 * '<S120>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Saturation'
 * '<S121>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Saturation Fdbk'
 * '<S122>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Sum'
 * '<S123>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Sum Fdbk'
 * '<S124>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Tracking Mode'
 * '<S125>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Tracking Mode Sum'
 * '<S126>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Tsamp - Integral'
 * '<S127>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Tsamp - Ngain'
 * '<S128>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/postSat Signal'
 * '<S129>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/preSat Signal'
 * '<S130>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Anti-windup/Disc. Clamping Parallel'
 * '<S131>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Anti-windup/Disc. Clamping Parallel/Dead Zone'
 * '<S132>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Anti-windup/Disc. Clamping Parallel/Dead Zone/Enabled'
 * '<S133>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/D Gain/Disabled'
 * '<S134>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Filter/Disabled'
 * '<S135>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Filter ICs/Disabled'
 * '<S136>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/I Gain/Internal Parameters'
 * '<S137>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Ideal P Gain/Passthrough'
 * '<S138>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Ideal P Gain Fdbk/Disabled'
 * '<S139>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Integrator/Discrete'
 * '<S140>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Integrator ICs/Internal IC'
 * '<S141>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/N Copy/Disabled wSignal Specification'
 * '<S142>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/N Gain/Disabled'
 * '<S143>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/P Copy/Disabled'
 * '<S144>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Parallel P Gain/Internal Parameters'
 * '<S145>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Reset Signal/External Reset'
 * '<S146>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Saturation/Enabled'
 * '<S147>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Saturation Fdbk/Disabled'
 * '<S148>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Sum/Sum_PI'
 * '<S149>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Sum Fdbk/Disabled'
 * '<S150>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Tracking Mode/Disabled'
 * '<S151>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Tracking Mode Sum/Passthrough'
 * '<S152>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Tsamp - Integral/Passthrough'
 * '<S153>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Tsamp - Ngain/Passthrough'
 * '<S154>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/postSat Signal/Forward_Path'
 * '<S155>' : 'speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/preSat Signal/Forward_Path'
 * '<S156>' : 'speedgoat_target_model_2021b/Planning and Control/Subsystem3/Subsystem1'
 * '<S157>' : 'speedgoat_target_model_2021b/Planning and Control/Subsystem3/Subsystem2'
 * '<S158>' : 'speedgoat_target_model_2021b/Planning and Control/Subsystem3/Subsystem1/Subsystem2'
 * '<S159>' : 'speedgoat_target_model_2021b/Planning and Control/Subsystem3/Subsystem1/Subsystem2/Subsystem'
 * '<S160>' : 'speedgoat_target_model_2021b/Planning and Control/Subsystem3/Subsystem1/Subsystem2/Subsystem1'
 * '<S161>' : 'speedgoat_target_model_2021b/Planning and Control/Subsystem3/Subsystem2/MATLAB Function'
 * '<S162>' : 'speedgoat_target_model_2021b/Planning and Control/Subsystem3/Subsystem2/Subsystem2'
 * '<S163>' : 'speedgoat_target_model_2021b/Planning and Control/Subsystem3/Subsystem2/Subsystem2/Subsystem'
 * '<S164>' : 'speedgoat_target_model_2021b/Planning and Control/Subsystem3/Subsystem2/Subsystem2/Subsystem1'
 * '<S165>' : 'speedgoat_target_model_2021b/UDP Data Receive/Byte Unpack: Fifth set'
 * '<S166>' : 'speedgoat_target_model_2021b/UDP Data Receive/Byte Unpack: First set'
 * '<S167>' : 'speedgoat_target_model_2021b/UDP Data Receive/Byte Unpack: Fourth set'
 * '<S168>' : 'speedgoat_target_model_2021b/UDP Data Receive/Byte Unpack: Second set'
 * '<S169>' : 'speedgoat_target_model_2021b/UDP Data Receive/Byte Unpack: Third set'
 * '<S170>' : 'speedgoat_target_model_2021b/UDP Data Receive/UDP Receive: Fifth set'
 * '<S171>' : 'speedgoat_target_model_2021b/UDP Data Receive/UDP Receive: First set'
 * '<S172>' : 'speedgoat_target_model_2021b/UDP Data Receive/UDP Receive: Fourth set'
 * '<S173>' : 'speedgoat_target_model_2021b/UDP Data Receive/UDP Receive: Second set'
 * '<S174>' : 'speedgoat_target_model_2021b/UDP Data Receive/UDP Receive: Third set'
 * '<S175>' : 'speedgoat_target_model_2021b/UDP Data Receive/UDP Receive: Fifth set/Compare To Constant'
 * '<S176>' : 'speedgoat_target_model_2021b/UDP Data Receive/UDP Receive: First set/Compare To Constant'
 * '<S177>' : 'speedgoat_target_model_2021b/UDP Data Receive/UDP Receive: Fourth set/Compare To Constant'
 * '<S178>' : 'speedgoat_target_model_2021b/UDP Data Receive/UDP Receive: Second set/Compare To Constant'
 * '<S179>' : 'speedgoat_target_model_2021b/UDP Data Receive/UDP Receive: Third set/Compare To Constant'
 * '<S180>' : 'speedgoat_target_model_2021b/UDP Data Send/Byte Pack: Controller Intermediate Signals'
 * '<S181>' : 'speedgoat_target_model_2021b/UDP Data Send/Byte Pack: Controller Output'
 */
#endif                          /* RTW_HEADER_speedgoat_target_model_2021b_h_ */
