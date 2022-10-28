/*
 * ModelWithControllersOnly.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "ModelWithControllersOnly".
 *
 * Model version              : 2.117
 * Simulink Coder version : 9.6 (R2021b) 14-May-2021
 * C++ source code generated on : Thu Sep 29 18:30:07 2022
 *
 * Target selection: slrealtime.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_ModelWithControllersOnly_h_
#define RTW_HEADER_ModelWithControllersOnly_h_
#include <cmath>
#include <cstring>
#include "rtwtypes.h"
#include "zero_crossing_types.h"
#include "rtw_extmode.h"
#include "sysran_types.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "verify/verifyIntrf.h"
#include "sf_runtime/slTestResult.h"
#include "sf_runtime/slTestTypes.h"
#include "ModelWithControllersOnly_types.h"

/* Shared type includes */
#include "multiword_types.h"

/* Child system includes */
#include "ModelWithControllersOnly_cal.h"
#include "rtGetNaN.h"
#include "rt_nonfinite.h"
#include "rtsplntypes.h"
#include "rt_zcfcn.h"
#include "crl_mutex.hpp"
#include "rt_assert.h"
#include "rt_defines.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
#define rtmGetContStateDisabled(rtm)   ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
#define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
#define rtmGetContStates(rtm)          ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
#define rtmSetContStates(rtm, val)     ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetIntgData
#define rtmGetIntgData(rtm)            ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
#define rtmSetIntgData(rtm, val)       ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
#define rtmGetOdeF(rtm)                ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
#define rtmSetOdeF(rtm, val)           ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
#define rtmGetOdeY(rtm)                ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
#define rtmSetOdeY(rtm, val)           ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
#define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
#define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
#define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
#define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
#define rtmGetdX(rtm)                  ((rtm)->derivs)
#endif

#ifndef rtmSetdX
#define rtmSetdX(rtm, val)             ((rtm)->derivs = (val))
#endif

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
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

#ifndef rtmTaskCounter
#define rtmTaskCounter(rtm, idx)       ((rtm)->Timing.TaskCounters.TID[(idx)])
#endif

/* Block signals (default storage) */
struct B_ModelWithControllersOnly_T {
  BusActorState BusCreator;            /* '<S120>/Bus Creator' */
  real_T VectorConcatenate[4];         /* '<S126>/Vector Concatenate' */
  real_T Integrator[4];                /* '<S214>/Integrator' */
  real_T Saturation;                   /* '<S5>/Saturation' */
  real_T VectorConcatenate3[2];        /* '<S126>/Vector Concatenate3' */
  real_T Integrator_m[2];              /* '<S149>/Integrator' */
  real_T Gain1;                        /* '<S3>/Gain1' */
  real_T TmpSignalConversionAtTransposeI[3];
  real_T Transpose[3];                 /* '<S3>/Transpose' */
  real_T Product;                      /* '<S8>/Product' */
  real_T Gain;                         /* '<S9>/Gain' */
  real_T Multiply;                     /* '<S8>/Multiply' */
  real_T CurvedSteadyYaw;              /* '<S8>/Gain1' */
  real_T Gain1_g;                      /* '<Root>/Gain1' */
  real_T Minus;                        /* '<S8>/Minus' */
  real_T YawRateFeedback;              /* '<S8>/Gain' */
  real_T UnitDelay;                    /* '<S8>/Unit Delay' */
  real_T UnitDelay_l;                  /* '<Root>/Unit Delay' */
  real_T Minus1;                       /* '<S8>/Minus1' */
  real_T SteeringDelay;                /* '<S8>/Gain2' */
  real_T Add;                          /* '<S8>/Add' */
  real_T Saturation_k;                 /* '<S8>/Saturation' */
  real_T Gain_n;                       /* '<Root>/Gain' */
  real_T TmpRTBAtIntegratorInport2;
  real_T v;                            /* '<S119>/Integrator' */
  real_T Square;                       /* '<S119>/Square' */
  real_T F_aero;                       /* '<S119>/Gain1' */
  real_T Gain2;                        /* '<S119>/Gain2' */
  real_T TmpRTBAtSumInport3;           /* '<S119>/Gain' */
  real_T Sum;                          /* '<S119>/Sum' */
  real_T T_eng;                        /* '<S119>/Gain4' */
  real_T BrakeSaturation;              /* '<S4>/Brake Saturation' */
  real_T TmpRTBAtIntegratorInport2_e;
  real_T v_g;                          /* '<S118>/Integrator' */
  real_T Square_g;                     /* '<S118>/Square' */
  real_T F_aero_i;                     /* '<S118>/Gain1' */
  real_T Gain2_d;                      /* '<S118>/Gain2' */
  real_T TmpRTBAtSumInport3_l;         /* '<S118>/Gain' */
  real_T Sum_n;                        /* '<S118>/Sum' */
  real_T T_eng_c;                      /* '<S118>/Gain4' */
  real_T TmpRTBAtIntegratorInport1;
  real_T TmpRTBAtIntegratorInport1_c;
  real_T ThrottleSaturation;           /* '<S4>/Throttle Saturation' */
  real_T Switch;                       /* '<S4>/Switch' */
  real_T Switch1;                      /* '<S4>/Switch1' */
  real_T gear_cmd;                     /* '<S5>/Constant9' */
  real_T Roll;                         /* '<S120>/Gain' */
  real_T Derivative;                   /* '<S120>/Derivative' */
  real_T Gain5;                        /* '<S120>/Gain5' */
  real_T Derivative1;                  /* '<S120>/Derivative1' */
  real_T Gain6;                        /* '<S120>/Gain6' */
  real_T Derivative2;                  /* '<S120>/Derivative2' */
  real_T Derivative3;                  /* '<S120>/Derivative3' */
  real_T Pitch;                        /* '<S120>/Gain1' */
  real_T Derivative4;                  /* '<S120>/Derivative4' */
  real_T Yawdeg;                       /* '<S120>/Gain2' */
  real_T Derivative5;                  /* '<S120>/Derivative5' */
  real_T Gain3;                        /* '<S120>/Gain3' */
  real_T Gain4;                        /* '<S120>/Gain4' */
  real_T EgoVehLongAcc;                /* '<S121>/Derivative' */
  real_T FirstOrderHold;               /* '<S121>/First Order Hold' */
  real_T TnetNm;                       /* '<S121>/First Order Hold1' */
  real_T wheelspeed;                   /* '<S121>/TireRadius' */
  real_T VectorConcatenate1[3];        /* '<S129>/Vector Concatenate1' */
  real_T TrigonometricFunction_o1;     /* '<S143>/Trigonometric Function' */
  real_T TrigonometricFunction_o2;     /* '<S143>/Trigonometric Function' */
  real_T VectorConcatenate1_b[3];      /* '<S215>/Vector Concatenate1' */
  real_T Product3;                     /* '<S143>/Product3' */
  real_T Product1;                     /* '<S143>/Product1' */
  real_T Product_k;                    /* '<S143>/Product' */
  real_T Product2;                     /* '<S143>/Product2' */
  real_T VectorConcatenate2[3];        /* '<S143>/Vector Concatenate2' */
  real_T Add1[3];                      /* '<S142>/Add1' */
  real_T Product_i[3];                 /* '<S142>/Product' */
  real_T SumofElements;                /* '<S142>/Sum of Elements' */
  real_T Sqrt;                         /* '<S142>/Sqrt' */
  real_T Product2_h;                   /* '<S142>/Product2' */
  real_T TrigonometricFunction;        /* '<S142>/Trigonometric Function' */
  real_T u[3];                         /* '<S142>/4' */
  real_T Tanh[3];                      /* '<S142>/Tanh' */
  real_T VectorConcatenate_b[6];       /* '<S142>/Vector Concatenate' */
  real_T Product1_a[6];                /* '<S142>/Product1' */
  real_T uAPabsRT[6];                  /* '<S142>/.5.*A.*Pabs.//R.//T' */
  real_T Product4[3];                  /* '<S142>/Product4' */
  real_T UnaryMinus1[3];               /* '<S129>/Unary Minus1' */
  real_T Add_d[3];                     /* '<S126>/Add' */
  real_T Add2[3];                      /* '<S142>/Add2' */
  real_T Product3_o[3];                /* '<S142>/Product3' */
  real_T UnaryMinus[3];                /* '<S129>/Unary Minus' */
  real_T Add1_d[3];                    /* '<S126>/Add1' */
  real_T UnitConversion1;              /* '<S207>/Unit Conversion1' */
  real_T TmpSignalConversionAtsincosInpo[3];
  real_T sincos_o1[3];                 /* '<S154>/sincos' */
  real_T sincos_o2[3];                 /* '<S154>/sincos' */
  real_T VectorConcatenate_f[9];       /* '<S158>/Vector Concatenate' */
  real_T Reshape9to3x3columnmajor[9];
                                /* '<S158>/Reshape (9) to [3x3] column-major' */
  real_T Transpose1[9];                /* '<S150>/Transpose1' */
  real_T VectorConcatenate_d[3];       /* '<S150>/Vector Concatenate' */
  real_T Reshape1[3];                  /* '<S156>/Reshape1' */
  real_T Product_f[3];                 /* '<S156>/Product' */
  real_T Reshape2[3];                  /* '<S156>/Reshape2' */
  real_T Add_e[3];                     /* '<S150>/Add' */
  real_T jxk;                          /* '<S159>/j x k' */
  real_T kxi;                          /* '<S159>/k x i' */
  real_T ixj;                          /* '<S159>/i x j' */
  real_T kxj;                          /* '<S160>/k x j' */
  real_T ixk;                          /* '<S160>/i x k' */
  real_T jxi;                          /* '<S160>/j x i' */
  real_T Sum_m[3];                     /* '<S157>/Sum' */
  real_T Add1_f[3];                    /* '<S150>/Add1' */
  real_T Reshape1_m[3];                /* '<S155>/Reshape1' */
  real_T Product_m[3];                 /* '<S155>/Product' */
  real_T Reshape2_p[3];                /* '<S155>/Reshape2' */
  real_T V_wb[3];                      /* '<S150>/Add4' */
  real_T sincos_o1_h[3];               /* '<S163>/sincos' */
  real_T sincos_o2_g[3];               /* '<S163>/sincos' */
  real_T VectorConcatenate_bj[9];      /* '<S170>/Vector Concatenate' */
  real_T Reshape9to3x3columnmajor_e[9];
                                /* '<S170>/Reshape (9) to [3x3] column-major' */
  real_T Transpose1_a[9];              /* '<S161>/Transpose1' */
  real_T VectorConcatenate_c[3];       /* '<S151>/Vector Concatenate' */
  real_T Reshape1_e[3];                /* '<S165>/Reshape1' */
  real_T Product_n[3];                 /* '<S165>/Product' */
  real_T Reshape2_m[3];                /* '<S165>/Reshape2' */
  real_T Add_p[3];                     /* '<S161>/Add' */
  real_T jxk_j;                        /* '<S171>/j x k' */
  real_T kxi_g;                        /* '<S171>/k x i' */
  real_T ixj_d;                        /* '<S171>/i x j' */
  real_T kxj_j;                        /* '<S172>/k x j' */
  real_T ixk_j;                        /* '<S172>/i x k' */
  real_T jxi_i;                        /* '<S172>/j x i' */
  real_T Sum_b[3];                     /* '<S166>/Sum' */
  real_T Add1_o[3];                    /* '<S161>/Add1' */
  real_T Reshape1_j[3];                /* '<S164>/Reshape1' */
  real_T Product_p[3];                 /* '<S164>/Product' */
  real_T Reshape2_c[3];                /* '<S164>/Reshape2' */
  real_T V_wb_p[3];                    /* '<S161>/Add4' */
  real_T Switch_d;                     /* '<S167>/Switch' */
  real_T Divide;                       /* '<S162>/Divide' */
  real_T Beta;                         /* '<S162>/Trigonometric Function' */
  real_T sincos_o1_n[3];               /* '<S173>/sincos' */
  real_T sincos_o2_k[3];               /* '<S173>/sincos' */
  real_T VectorConcatenate_f5[9];      /* '<S177>/Vector Concatenate' */
  real_T Reshape9to3x3columnmajor_n[9];
                                /* '<S177>/Reshape (9) to [3x3] column-major' */
  real_T Transpose1_k[9];              /* '<S152>/Transpose1' */
  real_T VectorConcatenate_o[3];       /* '<S152>/Vector Concatenate' */
  real_T Reshape1_n[3];                /* '<S175>/Reshape1' */
  real_T Product_g[3];                 /* '<S175>/Product' */
  real_T Reshape2_g[3];                /* '<S175>/Reshape2' */
  real_T Add_pa[3];                    /* '<S152>/Add' */
  real_T jxk_m;                        /* '<S178>/j x k' */
  real_T kxi_k;                        /* '<S178>/k x i' */
  real_T ixj_e;                        /* '<S178>/i x j' */
  real_T kxj_b;                        /* '<S179>/k x j' */
  real_T ixk_m;                        /* '<S179>/i x k' */
  real_T jxi_h;                        /* '<S179>/j x i' */
  real_T Sum_f[3];                     /* '<S176>/Sum' */
  real_T Add1_o0[3];                   /* '<S152>/Add1' */
  real_T Reshape1_k[3];                /* '<S174>/Reshape1' */
  real_T Product_ky[3];                /* '<S174>/Product' */
  real_T Reshape2_f[3];                /* '<S174>/Reshape2' */
  real_T V_wb_g[3];                    /* '<S152>/Add4' */
  real_T sincos_o1_b[3];               /* '<S181>/sincos' */
  real_T sincos_o2_gp[3];              /* '<S181>/sincos' */
  real_T VectorConcatenate_bx[9];      /* '<S185>/Vector Concatenate' */
  real_T Reshape9to3x3columnmajor_h[9];
                                /* '<S185>/Reshape (9) to [3x3] column-major' */
  real_T Transpose1_m[9];              /* '<S180>/Transpose1' */
  real_T VectorConcatenate_h[3];       /* '<S153>/Vector Concatenate' */
  real_T Reshape1_f[3];                /* '<S183>/Reshape1' */
  real_T Product_a[3];                 /* '<S183>/Product' */
  real_T Reshape2_mj[3];               /* '<S183>/Reshape2' */
  real_T Add_c[3];                     /* '<S180>/Add' */
  real_T jxk_b;                        /* '<S186>/j x k' */
  real_T kxi_e;                        /* '<S186>/k x i' */
  real_T ixj_c;                        /* '<S186>/i x j' */
  real_T kxj_f;                        /* '<S187>/k x j' */
  real_T ixk_p;                        /* '<S187>/i x k' */
  real_T jxi_h1;                       /* '<S187>/j x i' */
  real_T Sum_i[3];                     /* '<S184>/Sum' */
  real_T Add1_h[3];                    /* '<S180>/Add1' */
  real_T Reshape1_nq[3];               /* '<S182>/Reshape1' */
  real_T Product_iz[3];                /* '<S182>/Product' */
  real_T Reshape2_e[3];                /* '<S182>/Reshape2' */
  real_T V_wb_e[3];                    /* '<S180>/Add4' */
  real_T Product8[3];                  /* '<S148>/Product8' */
  real_T Product2_m[3];                /* '<S148>/Product2' */
  real_T Product14[3];                 /* '<S148>/Product14' */
  real_T Product3_b[3];                /* '<S148>/Product3' */
  real_T Product_po;                   /* '<S148>/Product' */
  real_T VectorConcatenate1_f[2];      /* '<S205>/Vector Concatenate1' */
  real_T VectorConcatenate1_b2[2];     /* '<S206>/Vector Concatenate1' */
  real_T VectorConcatenate1_c[2];      /* '<S209>/Vector Concatenate1' */
  real_T Product12[3];                 /* '<S148>/Product12' */
  real_T Product13[3];                 /* '<S148>/Product13' */
  real_T Divide1;                      /* '<S198>/Divide1' */
  real_T ax;                           /* '<S198>/Sum of Elements1' */
  real_T UnitConversion1_n;            /* '<S204>/Unit Conversion1' */
  real_T UnitConversion2;              /* '<S148>/Unit Conversion2' */
  real_T Product5;                     /* '<S148>/Product5' */
  real_T Divide_c;                     /* '<S198>/Divide' */
  real_T ay;                           /* '<S198>/Sum of Elements' */
  real_T UnitConversion;               /* '<S204>/Unit Conversion' */
  real_T UnitConversion1_f;            /* '<S148>/Unit Conversion1' */
  real_T Product6;                     /* '<S148>/Product6' */
  real_T Product7;                     /* '<S148>/Product7' */
  real_T VectorConcatenate1_n[6];      /* '<S148>/Vector Concatenate1' */
  real_T SumofElements1;               /* '<S148>/Sum of Elements1' */
  real_T VectorConcatenate2_o[6];      /* '<S148>/Vector Concatenate2' */
  real_T SumofElements2;               /* '<S148>/Sum of Elements2' */
  real_T VectorConcatenate3_l[6];      /* '<S148>/Vector Concatenate3' */
  real_T SumofElements3;               /* '<S148>/Sum of Elements3' */
  real_T UnaryMinus5;                  /* '<S148>/Unary Minus5' */
  real_T Sum1[2];                      /* '<S193>/Sum1' */
  real_T Sum_e;                        /* '<S193>/Sum' */
  real_T VectorConcatenate_j[2];       /* '<S144>/Vector Concatenate' */
  real_T VectorConcatenate1_i[2];      /* '<S144>/Vector Concatenate1' */
  real_T Switch_p;                     /* '<S200>/Switch' */
  real_T Divide_p;                     /* '<S195>/Divide' */
  real_T Beta_n;                       /* '<S195>/Trigonometric Function' */
  real_T az;                           /* '<S149>/Constant1' */
  real_T Constant10;                   /* '<S149>/Constant10' */
  real_T zddot;                        /* '<S149>/Constant3' */
  real_T Constant9;                    /* '<S149>/Constant9' */
  real_T sincos_o1_bs[3];              /* '<S197>/sincos' */
  real_T sincos_o2_f[3];               /* '<S197>/sincos' */
  real_T VectorConcatenate_jc[9];      /* '<S203>/Vector Concatenate' */
  real_T Reshape9to3x3columnmajor_l[9];
                                /* '<S203>/Reshape (9) to [3x3] column-major' */
  real_T lateral;                      /* '<S212>/lateral' */
  real_T Sum_c;                        /* '<S212>/Sum' */
  real_T Product1_h;                   /* '<S212>/Product1' */
  real_T lateral_g;                    /* '<S213>/lateral' */
  real_T Sum_g;                        /* '<S213>/Sum' */
  real_T Product1_af;                  /* '<S213>/Product1' */
  real_T ThrottleSaturation_p;         /* '<S127>/Throttle Saturation' */
  real_T TengineNm;            /* '<S127>/second-order low-pass torque model' */
  real_T DeadZone1;                    /* '<S127>/Dead Zone1' */
  real_T BrakeSaturation_m;            /* '<S127>/Brake Saturation' */
  real_T TbrakeNm;              /* '<S127>/second-order low-pass brake model' */
  real_T UnitDelay_p;                  /* '<S127>/Unit Delay' */
  real_T DeadZone;                     /* '<S127>/Dead Zone' */
  real_T SteeringSaturation;           /* '<S127>/Steering Saturation' */
  real_T steeringwheelangledeg;
                             /* '<S127>/second-order low-pass steering model' */
  real_T steeringwheelanglerad;        /* '<S127>/deg-to-rad1' */
  real_T uDLookupTable;                /* '<S218>/1-D Lookup Table' */
  real_T Product_kk;                   /* '<S218>/Product' */
  real_T uDLookupTable1;               /* '<S219>/1-D Lookup Table1' */
  real_T uDLookupTable_c;              /* '<S219>/1-D Lookup Table' */
  real_T Add_p3;                       /* '<S127>/Add' */
  real_T Gain2_o;                      /* '<S127>/Gain2' */
  real_T rollingresistance;            /* '<S127>/rolling resistance' */
  real_T TnetNm_e;                     /* '<S127>/Sum' */
  real_T torquetoaccelerationconversion1;
                               /* '<S127>/torque to acceleration conversion1' */
  real_T gear_cmd_e;                   /* '<S5>/Constant3' */
  real_T gear_cmd_i;                   /* '<S5>/Constant4' */
  real_T brake_cmdNm;                  /* '<S5>/From Workspace' */
  real_T steer_cmddeg;                 /* '<S5>/From Workspace1' */
  real_T torque_cmdNm;                 /* '<S5>/From Workspace4' */
  real_T steer_cmddeg_n;               /* '<S5>/Constant' */
  real_T torque_cmdNm_j;               /* '<S5>/Constant1' */
  real_T brake_cmdNm_j;                /* '<S5>/Constant2' */
  real_T gear_cmd_d;                   /* '<S5>/Constant5' */
  real_T TmpRTBAtEqual2Inport1;
  real_T Sign2;                        /* '<S12>/Sign2' */
  real_T TmpRTBAtEqual1Inport1;
  real_T xdot1;
  real_T Sign1;                        /* '<S12>/Sign1' */
  real_T xdot;
  real_T TmpRTBAtEqual5Inport2;
  real_T TmpRTBAtEqual4Inport2;
  real_T TmpRTBAtMultiplyInport1;
  real_T TmpRTBAtSwitchCaseInport1;
  real_T xdot_p;
  real_T error;                        /* '<S2>/Minus' */
  real_T Merge;                        /* '<S11>/Merge' */
  real_T Multiply_n;                   /* '<S10>/Multiply' */
  real_T Switch_l;                     /* '<S10>/Switch' */
  real_T Switch1_l;                    /* '<S10>/Switch1' */
  real_T Gain_g;                       /* '<S118>/Gain' */
  real_T Gain_ne;                      /* '<S119>/Gain' */
  real_T TmpRTBAtMATLABFunction3Inport1;/* '<S120>/Derivative3' */
  real_T TmpRTBAtMATLABFunction3Inport2;/* '<S120>/Derivative4' */
  real_T TmpRTBAtMATLABFunction3Inport3;/* '<S120>/Derivative5' */
  real_T X;
  real_T TmpRTBAtMATLABFunctionInport2;/* '<S120>/Gain3' */
  real_T TmpRTBAtMATLABFunctionInport3;/* '<S120>/Gain4' */
  real_T Xdot;
  real_T TmpRTBAtMATLABFunction1Inport2;/* '<S120>/Gain5' */
  real_T TmpRTBAtMATLABFunction1Inport3;/* '<S120>/Gain6' */
  real_T TmpRTBAtMATLABFunction2Inport1;/* '<S120>/Derivative' */
  real_T TmpRTBAtMATLABFunction2Inport2;/* '<S120>/Derivative1' */
  real_T TmpRTBAtMATLABFunction2Inport3;/* '<S120>/Derivative2' */
  real_T Tnet;                         /* '<S127>/compute net torque' */
  real_T TmpSignalConversionAtSFunctionI[3];/* '<S126>/vehicle model' */
  real_T TmpSignalConversionAtSFunctio_n[3];/* '<S126>/vehicle model' */
  real_T yOut[3];                      /* '<S126>/vehicle model' */
  real_T FBody[3];                     /* '<S126>/vehicle model' */
  real_T MBody[3];                     /* '<S126>/vehicle model' */
  real_T FOut[6];                      /* '<S126>/vehicle model' */
  real_T FTire[6];                     /* '<S126>/vehicle model' */
  real_T Fg[3];                        /* '<S126>/vehicle model' */
  real_T wheelInfo[4];                 /* '<S126>/vehicle model' */
  real_T stateDer[4];                  /* '<S126>/vehicle model' */
  real_T status;                       /* '<S126>/vehicle model' */
  real_T y[2];                         /* '<S149>/COMB2I' */
  real_T Abs;                          /* '<S200>/Abs' */
  real_T Fcn;                          /* '<S200>/Fcn' */
  real_T Abs_d;                        /* '<S167>/Abs' */
  real_T Fcn_g;                        /* '<S167>/Fcn' */
  real_T AngularVelocity[3];           /* '<S120>/MATLAB Function3' */
  real_T ProportionalGain;             /* '<S105>/Proportional Gain' */
  real_T Integrator_f;                 /* '<S100>/Integrator' */
  real_T Sum_mr;                       /* '<S109>/Sum' */
  real_T ZeroGain;                     /* '<S91>/ZeroGain' */
  real_T DeadZone_f;                   /* '<S93>/DeadZone' */
  real_T SignPreSat;                   /* '<S91>/SignPreSat' */
  real_T IntegralGain;                 /* '<S97>/Integral Gain' */
  real_T SignPreIntegrator;            /* '<S91>/SignPreIntegrator' */
  real_T Switch_m;                     /* '<S91>/Switch' */
  real_T ProportionalGain_o;           /* '<S54>/Proportional Gain' */
  real_T Integrator_c;                 /* '<S49>/Integrator' */
  real_T Sum_gu;                       /* '<S58>/Sum' */
  real_T ZeroGain_p;                   /* '<S40>/ZeroGain' */
  real_T DeadZone_o;                   /* '<S42>/DeadZone' */
  real_T SignPreSat_n;                 /* '<S40>/SignPreSat' */
  real_T IntegralGain_e;               /* '<S46>/Integral Gain' */
  real_T SignPreIntegrator_k;          /* '<S40>/SignPreIntegrator' */
  real_T Switch_a;                     /* '<S40>/Switch' */
  real_T Abs_k;                        /* '<S10>/Abs' */
  real_T steerCmd;                     /* '<S1>/Kinematic' */
  int8_T DataTypeConv1;                /* '<S91>/DataTypeConv1' */
  int8_T DataTypeConv2;                /* '<S91>/DataTypeConv2' */
  int8_T DataTypeConv1_j;              /* '<S40>/DataTypeConv1' */
  int8_T DataTypeConv2_g;              /* '<S40>/DataTypeConv2' */
  boolean_T Compare;                   /* '<S168>/Compare' */
  boolean_T Compare_f;                 /* '<S169>/Compare' */
  boolean_T LogicalOperator;           /* '<S167>/Logical Operator' */
  boolean_T Compare_i;                 /* '<S201>/Compare' */
  boolean_T Compare_p;                 /* '<S202>/Compare' */
  boolean_T LogicalOperator_e;         /* '<S200>/Logical Operator' */
  boolean_T Equal2;                    /* '<S12>/Equal2' */
  boolean_T Equal6;                    /* '<S12>/Equal6' */
  boolean_T AND3;                      /* '<S12>/AND3' */
  boolean_T Equal1;                    /* '<S12>/Equal1' */
  boolean_T Equal3;                    /* '<S12>/Equal3' */
  boolean_T AND1;                      /* '<S12>/AND1' */
  boolean_T OR;                        /* '<S12>/OR' */
  boolean_T NOT1;                      /* '<S12>/NOT1' */
  boolean_T Equal5;                    /* '<S12>/Equal5' */
  boolean_T Equal4;                    /* '<S12>/Equal4' */
  boolean_T AND2;                      /* '<S12>/AND2' */
  boolean_T NOT;                       /* '<S12>/NOT' */
  boolean_T DataTypeConversion;        /* '<S11>/Data Type Conversion' */
  boolean_T NotEqual;                  /* '<S91>/NotEqual' */
  boolean_T Equal1_e;                  /* '<S91>/Equal1' */
  boolean_T AND3_h;                    /* '<S91>/AND3' */
  boolean_T NotEqual_m;                /* '<S40>/NotEqual' */
  boolean_T Equal1_j;                  /* '<S40>/Equal1' */
  boolean_T AND3_k;                    /* '<S40>/AND3' */
};

/* Block states (default storage) for system '<Root>' */
struct DW_ModelWithControllersOnly_T {
  slTestBlkInfo_t Assertion_sltestBlkInfo;/* '<S12>/Assertion' */
  slTestBlkInfo_t Assertion1_sltestBlkInfo;/* '<S12>/Assertion1' */
  real_T UnitDelay_DSTATE;             /* '<S8>/Unit Delay' */
  real_T UnitDelay_DSTATE_e;           /* '<Root>/Unit Delay' */
  real_T secondorderlowpasstorquemodel_s[2];
                               /* '<S127>/second-order low-pass torque model' */
  real_T secondorderlowpassbrakemodel_st[2];
                                /* '<S127>/second-order low-pass brake model' */
  real_T UnitDelay_DSTATE_l;           /* '<S127>/Unit Delay' */
  real_T secondorderlowpasssteeringmodel[2];
                             /* '<S127>/second-order low-pass steering model' */
  real_T Integrator_DSTATE;            /* '<S100>/Integrator' */
  real_T Integrator_DSTATE_o;          /* '<S49>/Integrator' */
  real_T TmpRTBAtEqual1Inport1_Buf0;   /* synthesized block */
  real_T TmpRTBAtEqual1Inport1_Buf1;   /* synthesized block */
  real_T TmpRTBAtEqual1Inport1_Buf2;   /* synthesized block */
  real_T TmpRTBAtEqual2Inport1_Buf0;   /* synthesized block */
  real_T TmpRTBAtEqual2Inport1_Buf1;   /* synthesized block */
  real_T TmpRTBAtEqual2Inport1_Buf2;   /* synthesized block */
  real_T TmpRTBAtEqual4Inport2_Buf0;   /* synthesized block */
  real_T TmpRTBAtEqual4Inport2_Buf1;   /* synthesized block */
  real_T TmpRTBAtEqual4Inport2_Buf2;   /* synthesized block */
  real_T TmpRTBAtEqual5Inport2_Buf0;   /* synthesized block */
  real_T TmpRTBAtEqual5Inport2_Buf1;   /* synthesized block */
  real_T TmpRTBAtEqual5Inport2_Buf2;   /* synthesized block */
  real_T xdot_Buf[2];                  /* synthesized block */
  real_T xdot1_Buf[2];                 /* synthesized block */
  real_T TmpRTBAtMultiplyInport1_Buf0; /* synthesized block */
  real_T TmpRTBAtMultiplyInport1_Buf1; /* synthesized block */
  real_T TmpRTBAtMultiplyInport1_Buf2; /* synthesized block */
  real_T TmpRTBAtSwitchCaseInport1_Buf0;/* synthesized block */
  real_T TmpRTBAtSwitchCaseInport1_Buf1;/* synthesized block */
  real_T TmpRTBAtSwitchCaseInport1_Buf2;/* synthesized block */
  real_T xdot_Buf_d[2];                /* synthesized block */
  real_T TmpRTBAtIntegratorInport2_Buf[2];/* synthesized block */
  real_T TmpRTBAtSumInport3_Buf[2];    /* synthesized block */
  real_T TmpRTBAtIntegratorInport2_Buf_m[2];/* synthesized block */
  real_T TmpRTBAtSumInport3_Buf_f[2];  /* synthesized block */
  real_T TmpRTBAtIntegratorInport1_Buf[2];/* synthesized block */
  real_T TmpRTBAtIntegratorInport1_Buf_p[2];/* synthesized block */
  real_T TimeStampA;                   /* '<S120>/Derivative' */
  real_T LastUAtTimeA;                 /* '<S120>/Derivative' */
  real_T TimeStampB;                   /* '<S120>/Derivative' */
  real_T LastUAtTimeB;                 /* '<S120>/Derivative' */
  real_T TimeStampA_l;                 /* '<S120>/Derivative1' */
  real_T LastUAtTimeA_d;               /* '<S120>/Derivative1' */
  real_T TimeStampB_l;                 /* '<S120>/Derivative1' */
  real_T LastUAtTimeB_o;               /* '<S120>/Derivative1' */
  real_T TimeStampA_e;                 /* '<S120>/Derivative2' */
  real_T LastUAtTimeA_n;               /* '<S120>/Derivative2' */
  real_T TimeStampB_k;                 /* '<S120>/Derivative2' */
  real_T LastUAtTimeB_d;               /* '<S120>/Derivative2' */
  real_T TimeStampA_o;                 /* '<S120>/Derivative3' */
  real_T LastUAtTimeA_nv;              /* '<S120>/Derivative3' */
  real_T TimeStampB_lf;                /* '<S120>/Derivative3' */
  real_T LastUAtTimeB_g;               /* '<S120>/Derivative3' */
  real_T TimeStampA_l5;                /* '<S120>/Derivative4' */
  real_T LastUAtTimeA_k;               /* '<S120>/Derivative4' */
  real_T TimeStampB_e;                 /* '<S120>/Derivative4' */
  real_T LastUAtTimeB_p;               /* '<S120>/Derivative4' */
  real_T TimeStampA_k;                 /* '<S120>/Derivative5' */
  real_T LastUAtTimeA_h;               /* '<S120>/Derivative5' */
  real_T TimeStampB_l4;                /* '<S120>/Derivative5' */
  real_T LastUAtTimeB_f;               /* '<S120>/Derivative5' */
  real_T RateTransition_Buf[2];        /* '<S120>/Rate Transition' */
  real_T RateTransition1_Buf[2];       /* '<S120>/Rate Transition1' */
  real_T RateTransition2_Buf[2];       /* '<S120>/Rate Transition2' */
  real_T TmpRTBAtMATLABFunction1Inport1_[2];/* synthesized block */
  real_T TmpRTBAtMATLABFunction1Inport2_[2];/* synthesized block */
  real_T TmpRTBAtMATLABFunction1Inport3_;/* synthesized block */
  real_T TmpRTBAtMATLABFunction1Inport_f;/* synthesized block */
  real_T TmpRTBAtMATLABFunction1Inport_o;/* synthesized block */
  real_T TmpRTBAtMATLABFunction2Inport1_[2];/* synthesized block */
  real_T TmpRTBAtMATLABFunction2Inport2_[2];/* synthesized block */
  real_T TmpRTBAtMATLABFunction2Inport3_[2];/* synthesized block */
  real_T TmpRTBAtMATLABFunction3Inport1_[2];/* synthesized block */
  real_T TmpRTBAtMATLABFunction3Inport2_[2];/* synthesized block */
  real_T TmpRTBAtMATLABFunction3Inport3_[2];/* synthesized block */
  real_T TmpRTBAtMATLABFunctionInport1_B[2];/* synthesized block */
  real_T TmpRTBAtMATLABFunctionInport2_B[2];/* synthesized block */
  real_T TmpRTBAtMATLABFunctionInport3_B;/* synthesized block */
  real_T TmpRTBAtMATLABFunctionInport3_c;/* synthesized block */
  real_T TmpRTBAtMATLABFunctionInport3_f;/* synthesized block */
  real_T TimeStampA_h;                 /* '<S121>/Derivative' */
  real_T LastUAtTimeA_o;               /* '<S121>/Derivative' */
  real_T TimeStampB_h;                 /* '<S121>/Derivative' */
  real_T LastUAtTimeB_c;               /* '<S121>/Derivative' */
  real_T Tk;                           /* '<S121>/First Order Hold' */
  real_T Ck;                           /* '<S121>/First Order Hold' */
  real_T Mk;                           /* '<S121>/First Order Hold' */
  real_T Uk;                           /* '<S121>/First Order Hold' */
  real_T Tk_d;                         /* '<S121>/First Order Hold1' */
  real_T Ck_o;                         /* '<S121>/First Order Hold1' */
  real_T Mk_n;                         /* '<S121>/First Order Hold1' */
  real_T Uk_d;                         /* '<S121>/First Order Hold1' */
  real_T secondorderlowpasstorquemodel_t;
                               /* '<S127>/second-order low-pass torque model' */
  real_T secondorderlowpassbrakemodel_tm;
                                /* '<S127>/second-order low-pass brake model' */
  real_T secondorderlowpasssteeringmod_c;
                             /* '<S127>/second-order low-pass steering model' */
  real_T m_bpLambda;                   /* '<S218>/1-D Lookup Table' */
  real_T m_yyA;                        /* '<S218>/1-D Lookup Table' */
  real_T m_yyB;                        /* '<S218>/1-D Lookup Table' */
  real_T m_yy2;                        /* '<S218>/1-D Lookup Table' */
  real_T m_up[2];                      /* '<S218>/1-D Lookup Table' */
  real_T m_y2[2];                      /* '<S218>/1-D Lookup Table' */
  real_T prevBp0AndTableData[4];       /* '<S218>/1-D Lookup Table' */
  real_T Assertion_sltestLastResultTime;/* '<S12>/Assertion' */
  real_T Assertion1_sltestLastResultTime;/* '<S12>/Assertion1' */
  struct {
    void *LoggedData[2];
  } Scope_PWORK;                       /* '<Root>/Scope' */

  void* TmpRTBAtEqual1Inport1_d0_SEMAPH;/* synthesized block */
  void* TmpRTBAtEqual2Inport1_d0_SEMAPH;/* synthesized block */
  void* TmpRTBAtEqual4Inport2_d0_SEMAPH;/* synthesized block */
  void* TmpRTBAtEqual5Inport2_d0_SEMAPH;/* synthesized block */
  void* TmpRTBAtMultiplyInport1_d0_SEMA;/* synthesized block */
  void* TmpRTBAtSwitchCaseInport1_d0_SE;/* synthesized block */
  struct {
    void *LoggedData;
  } Scope1_PWORK;                      /* '<S120>/Scope1' */

  void* TmpRTBAtMATLABFunction1Inport_k;/* synthesized block */
  void* TmpRTBAtMATLABFunctionInport3_d;/* synthesized block */
  struct {
    void *LoggedData[2];
  } Scope1_PWORK_k;                    /* '<S127>/Scope1' */

  struct {
    void *LoggedData[2];
  } Scope2_PWORK;                      /* '<S127>/Scope2' */

  struct {
    void *LoggedData;
  } Scope3_PWORK;                      /* '<S127>/Scope3' */

  void* m_bpDataSet;                   /* '<S218>/1-D Lookup Table' */
  void* TWork[6];                      /* '<S218>/1-D Lookup Table' */
  void* SWork[9];                      /* '<S218>/1-D Lookup Table' */
  struct {
    void *TimePtr;
    void *DataPtr;
    void *RSimInfoPtr;
  } FromWorkspace_PWORK;               /* '<S5>/From Workspace' */

  struct {
    void *TimePtr;
    void *DataPtr;
    void *RSimInfoPtr;
  } FromWorkspace1_PWORK;              /* '<S5>/From Workspace1' */

  struct {
    void *TimePtr;
    void *DataPtr;
    void *RSimInfoPtr;
  } FromWorkspace4_PWORK;              /* '<S5>/From Workspace4' */

  struct {
    void *LoggedData;
  } Scope_PWORK_c;                     /* '<S120>/Scope' */

  int32_T Assertion_sltestFinalResult; /* '<S12>/Assertion' */
  int32_T Assertion_sltestCurrentResult;/* '<S12>/Assertion' */
  int32_T Assertion1_sltestFinalResult;/* '<S12>/Assertion1' */
  int32_T Assertion1_sltestCurrentResult;/* '<S12>/Assertion1' */
  uint32_T m_bpIndex;                  /* '<S218>/1-D Lookup Table' */
  int_T Integrator_IWORK;              /* '<S214>/Integrator' */
  int_T Integrator_IWORK_c;            /* '<S149>/Integrator' */
  struct {
    int_T PrevIndex;
  } FromWorkspace_IWORK;               /* '<S5>/From Workspace' */

  struct {
    int_T PrevIndex;
  } FromWorkspace1_IWORK;              /* '<S5>/From Workspace1' */

  struct {
    int_T PrevIndex;
  } FromWorkspace4_IWORK;              /* '<S5>/From Workspace4' */

  int8_T TmpRTBAtEqual1Inport1_LstBufWR;/* synthesized block */
  int8_T TmpRTBAtEqual1Inport1_RDBuf;  /* synthesized block */
  int8_T TmpRTBAtEqual2Inport1_LstBufWR;/* synthesized block */
  int8_T TmpRTBAtEqual2Inport1_RDBuf;  /* synthesized block */
  int8_T TmpRTBAtEqual4Inport2_LstBufWR;/* synthesized block */
  int8_T TmpRTBAtEqual4Inport2_RDBuf;  /* synthesized block */
  int8_T TmpRTBAtEqual5Inport2_LstBufWR;/* synthesized block */
  int8_T TmpRTBAtEqual5Inport2_RDBuf;  /* synthesized block */
  int8_T xdot_RdBufIdx;                /* synthesized block */
  int8_T xdot_WrBufIdx;                /* synthesized block */
  int8_T xdot1_RdBufIdx;               /* synthesized block */
  int8_T xdot1_WrBufIdx;               /* synthesized block */
  int8_T TmpRTBAtMultiplyInport1_LstBufW;/* synthesized block */
  int8_T TmpRTBAtMultiplyInport1_RDBuf;/* synthesized block */
  int8_T TmpRTBAtSwitchCaseInport1_LstBu;/* synthesized block */
  int8_T TmpRTBAtSwitchCaseInport1_RDBuf;/* synthesized block */
  int8_T xdot_RdBufIdx_m;              /* synthesized block */
  int8_T xdot_WrBufIdx_b;              /* synthesized block */
  int8_T TmpRTBAtIntegratorInport2_RdBuf;/* synthesized block */
  int8_T TmpRTBAtIntegratorInport2_WrBuf;/* synthesized block */
  int8_T TmpRTBAtSumInport3_RdBufIdx;  /* synthesized block */
  int8_T TmpRTBAtSumInport3_WrBufIdx;  /* synthesized block */
  int8_T TmpRTBAtIntegratorInport2_RdB_d;/* synthesized block */
  int8_T TmpRTBAtIntegratorInport2_WrB_d;/* synthesized block */
  int8_T TmpRTBAtSumInport3_RdBufIdx_l;/* synthesized block */
  int8_T TmpRTBAtSumInport3_WrBufIdx_o;/* synthesized block */
  int8_T TmpRTBAtIntegratorInport1_RdBuf;/* synthesized block */
  int8_T TmpRTBAtIntegratorInport1_WrBuf;/* synthesized block */
  int8_T TmpRTBAtIntegratorInport1_RdB_k;/* synthesized block */
  int8_T TmpRTBAtIntegratorInport1_WrB_k;/* synthesized block */
  int8_T RateTransition_RdBufIdx;      /* '<S120>/Rate Transition' */
  int8_T RateTransition_WrBufIdx;      /* '<S120>/Rate Transition' */
  int8_T RateTransition1_RdBufIdx;     /* '<S120>/Rate Transition1' */
  int8_T RateTransition1_WrBufIdx;     /* '<S120>/Rate Transition1' */
  int8_T RateTransition2_RdBufIdx;     /* '<S120>/Rate Transition2' */
  int8_T RateTransition2_WrBufIdx;     /* '<S120>/Rate Transition2' */
  int8_T TmpRTBAtMATLABFunction1Inport_h;/* synthesized block */
  int8_T TmpRTBAtMATLABFunction1Inport_j;/* synthesized block */
  int8_T TmpRTBAtMATLABFunction1Inport_p;/* synthesized block */
  int8_T TmpRTBAtMATLABFunction1Inport_m;/* synthesized block */
  int8_T TmpRTBAtMATLABFunction1Inpor_oh;/* synthesized block */
  int8_T TmpRTBAtMATLABFunction1Inport_b;/* synthesized block */
  int8_T TmpRTBAtMATLABFunction2Inport_e;/* synthesized block */
  int8_T TmpRTBAtMATLABFunction2Inport_d;/* synthesized block */
  int8_T TmpRTBAtMATLABFunction2Inport_m;/* synthesized block */
  int8_T TmpRTBAtMATLABFunction2Inport_k;/* synthesized block */
  int8_T TmpRTBAtMATLABFunction2Inport_a;/* synthesized block */
  int8_T TmpRTBAtMATLABFunction2Inport_p;/* synthesized block */
  int8_T TmpRTBAtMATLABFunction3Inport_a;/* synthesized block */
  int8_T TmpRTBAtMATLABFunction3Inport_g;/* synthesized block */
  int8_T TmpRTBAtMATLABFunction3Inport_f;/* synthesized block */
  int8_T TmpRTBAtMATLABFunction3Inport_b;/* synthesized block */
  int8_T TmpRTBAtMATLABFunction3Inpor_gh;/* synthesized block */
  int8_T TmpRTBAtMATLABFunction3Inpor_bx;/* synthesized block */
  int8_T TmpRTBAtMATLABFunctionInport1_R;/* synthesized block */
  int8_T TmpRTBAtMATLABFunctionInport1_W;/* synthesized block */
  int8_T TmpRTBAtMATLABFunctionInport2_R;/* synthesized block */
  int8_T TmpRTBAtMATLABFunctionInport2_W;/* synthesized block */
  int8_T TmpRTBAtMATLABFunctionInport3_L;/* synthesized block */
  int8_T TmpRTBAtMATLABFunctionInport3_R;/* synthesized block */
  int8_T SwitchCase_ActiveSubsystem;   /* '<S11>/Switch Case' */
  int8_T Reverse_SubsysRanBC;          /* '<S11>/Reverse' */
  int8_T Integrator_PrevResetState;    /* '<S100>/Integrator' */
  int8_T Forward_SubsysRanBC;          /* '<S11>/Forward' */
  int8_T Integrator_PrevResetState_n;  /* '<S49>/Integrator' */
  uint8_T reCalcSecDerivFirstDimCoeffs;/* '<S218>/1-D Lookup Table' */
};

/* Continuous states (default storage) */
struct X_ModelWithControllersOnly_T {
  real_T Integrator_CSTATE[4];         /* '<S214>/Integrator' */
  real_T Integrator_CSTATE_b[2];       /* '<S149>/Integrator' */
  real_T Integrator_CSTATE_o;          /* '<S119>/Integrator' */
  real_T Integrator_CSTATE_c;          /* '<S118>/Integrator' */
  real_T lateral_CSTATE;               /* '<S212>/lateral' */
  real_T lateral_CSTATE_k;             /* '<S213>/lateral' */
};

/* State derivatives (default storage) */
struct XDot_ModelWithControllersOnly_T {
  real_T Integrator_CSTATE[4];         /* '<S214>/Integrator' */
  real_T Integrator_CSTATE_b[2];       /* '<S149>/Integrator' */
  real_T Integrator_CSTATE_o;          /* '<S119>/Integrator' */
  real_T Integrator_CSTATE_c;          /* '<S118>/Integrator' */
  real_T lateral_CSTATE;               /* '<S212>/lateral' */
  real_T lateral_CSTATE_k;             /* '<S213>/lateral' */
};

/* State disabled  */
struct XDis_ModelWithControllersOnly_T {
  boolean_T Integrator_CSTATE[4];      /* '<S214>/Integrator' */
  boolean_T Integrator_CSTATE_b[2];    /* '<S149>/Integrator' */
  boolean_T Integrator_CSTATE_o;       /* '<S119>/Integrator' */
  boolean_T Integrator_CSTATE_c;       /* '<S118>/Integrator' */
  boolean_T lateral_CSTATE;            /* '<S212>/lateral' */
  boolean_T lateral_CSTATE_k;          /* '<S213>/lateral' */
};

/* Zero-crossing (trigger) state */
struct PrevZCX_ModelWithControllersO_T {
  ZCSigState Integrator_Reset_ZCE;     /* '<S119>/Integrator' */
  ZCSigState Integrator_Reset_ZCE_c;   /* '<S118>/Integrator' */
};

#ifndef ODE3_INTG
#define ODE3_INTG

/* ODE3 Integration Data */
struct ODE3_IntgData {
  real_T *y;                           /* output */
  real_T *f[3];                        /* derivatives */
};

#endif

/* Parameters (default storage) */
struct P_ModelWithControllersOnly_T_ {
  real_T Saturation_UpperSat;          /* Expression: inf
                                        * Referenced by: '<S5>/Saturation'
                                        */
  real_T FirstOrderHold_ErrTol;        /* Expression: inf
                                        * Referenced by: '<S121>/First Order Hold'
                                        */
  real_T FirstOrderHold1_ErrTol;       /* Expression: inf
                                        * Referenced by: '<S121>/First Order Hold1'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_ModelWithControllersO_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_ModelWithControllersOnly_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[10];
  real_T odeF[3][10];
  ODE3_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    uint32_T clockTick2;
    uint32_T clockTickH2;
    boolean_T firstInitCondFlag;
    struct {
      uint8_T TID[3];
      uint8_T cLimit[3];
    } TaskCounters;

    struct {
      uint8_T TID1_2;
    } RateInteraction;

    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[3];
  } Timing;
};

/* Block parameters (default storage) */
#ifdef __cplusplus

extern "C" {

#endif

  extern P_ModelWithControllersOnly_T ModelWithControllersOnly_P;

#ifdef __cplusplus

}
#endif

/* Block signals (default storage) */
#ifdef __cplusplus

extern "C" {

#endif

  extern struct B_ModelWithControllersOnly_T ModelWithControllersOnly_B;

#ifdef __cplusplus

}
#endif

/* Continuous states (default storage) */
extern X_ModelWithControllersOnly_T ModelWithControllersOnly_X;

/* Block states (default storage) */
extern struct DW_ModelWithControllersOnly_T ModelWithControllersOnly_DW;

/* Zero-crossing (trigger) state */
extern PrevZCX_ModelWithControllersO_T ModelWithControllersOnl_PrevZCX;

/* External data declarations for dependent source files */
extern const real_T ModelWithControllersOnly_RGND;/* real_T ground */
extern void rate_scheduler(void);

#ifdef __cplusplus

extern "C" {

#endif

  /* Model entry point functions */
  extern void ModelWithControllersOnly_initialize(void);
  extern void ModelWithControllersOnly_step0(void);
  extern void ModelWithControllersOnly_step2(void);
  extern void ModelWithControllersOnly_terminate(void);

#ifdef __cplusplus

}
#endif

/* Real-time Model object */
#ifdef __cplusplus

extern "C" {

#endif

  extern RT_MODEL_ModelWithControllers_T *const ModelWithControllersOnly_M;

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
 * '<Root>' : 'ModelWithControllersOnly'
 * '<S1>'   : 'ModelWithControllersOnly/Lateral Controller Stanley'
 * '<S2>'   : 'ModelWithControllersOnly/Longitudinal Controller Stanley'
 * '<S3>'   : 'ModelWithControllersOnly/Subsystem1'
 * '<S4>'   : 'ModelWithControllersOnly/Subsystem3'
 * '<S5>'   : 'ModelWithControllersOnly/Vehicle Dynamics1'
 * '<S6>'   : 'ModelWithControllersOnly/Lateral Controller Stanley/Dynamic'
 * '<S7>'   : 'ModelWithControllersOnly/Lateral Controller Stanley/Kinematic'
 * '<S8>'   : 'ModelWithControllersOnly/Lateral Controller Stanley/Dynamic/Dynamic Enabled'
 * '<S9>'   : 'ModelWithControllersOnly/Lateral Controller Stanley/Dynamic/Dynamic Enabled/Radians to Degrees'
 * '<S10>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/Convert Command '
 * '<S11>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller'
 * '<S12>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/Verify Direction'
 * '<S13>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward'
 * '<S14>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse'
 * '<S15>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward '
 * '<S16>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Anti-windup'
 * '<S17>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /D Gain'
 * '<S18>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Filter'
 * '<S19>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Filter ICs'
 * '<S20>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /I Gain'
 * '<S21>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Ideal P Gain'
 * '<S22>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Ideal P Gain Fdbk'
 * '<S23>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Integrator'
 * '<S24>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Integrator ICs'
 * '<S25>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /N Copy'
 * '<S26>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /N Gain'
 * '<S27>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /P Copy'
 * '<S28>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Parallel P Gain'
 * '<S29>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Reset Signal'
 * '<S30>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Saturation'
 * '<S31>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Saturation Fdbk'
 * '<S32>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Sum'
 * '<S33>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Sum Fdbk'
 * '<S34>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Tracking Mode'
 * '<S35>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Tracking Mode Sum'
 * '<S36>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Tsamp - Integral'
 * '<S37>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Tsamp - Ngain'
 * '<S38>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /postSat Signal'
 * '<S39>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /preSat Signal'
 * '<S40>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Anti-windup/Disc. Clamping Parallel'
 * '<S41>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Anti-windup/Disc. Clamping Parallel/Dead Zone'
 * '<S42>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Anti-windup/Disc. Clamping Parallel/Dead Zone/Enabled'
 * '<S43>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /D Gain/Disabled'
 * '<S44>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Filter/Disabled'
 * '<S45>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Filter ICs/Disabled'
 * '<S46>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /I Gain/Internal Parameters'
 * '<S47>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Ideal P Gain/Passthrough'
 * '<S48>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Ideal P Gain Fdbk/Disabled'
 * '<S49>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Integrator/Discrete'
 * '<S50>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Integrator ICs/Internal IC'
 * '<S51>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /N Copy/Disabled wSignal Specification'
 * '<S52>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /N Gain/Disabled'
 * '<S53>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /P Copy/Disabled'
 * '<S54>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Parallel P Gain/Internal Parameters'
 * '<S55>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Reset Signal/External Reset'
 * '<S56>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Saturation/Enabled'
 * '<S57>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Saturation Fdbk/Disabled'
 * '<S58>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Sum/Sum_PI'
 * '<S59>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Sum Fdbk/Disabled'
 * '<S60>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Tracking Mode/Disabled'
 * '<S61>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Tracking Mode Sum/Passthrough'
 * '<S62>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Tsamp - Integral/Passthrough'
 * '<S63>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /Tsamp - Ngain/Passthrough'
 * '<S64>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /postSat Signal/Forward_Path'
 * '<S65>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Forward/PI Forward /preSat Signal/Forward_Path'
 * '<S66>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse'
 * '<S67>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Anti-windup'
 * '<S68>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/D Gain'
 * '<S69>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Filter'
 * '<S70>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Filter ICs'
 * '<S71>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/I Gain'
 * '<S72>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Ideal P Gain'
 * '<S73>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Ideal P Gain Fdbk'
 * '<S74>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Integrator'
 * '<S75>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Integrator ICs'
 * '<S76>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/N Copy'
 * '<S77>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/N Gain'
 * '<S78>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/P Copy'
 * '<S79>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Parallel P Gain'
 * '<S80>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Reset Signal'
 * '<S81>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Saturation'
 * '<S82>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Saturation Fdbk'
 * '<S83>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Sum'
 * '<S84>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Sum Fdbk'
 * '<S85>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Tracking Mode'
 * '<S86>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Tracking Mode Sum'
 * '<S87>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Tsamp - Integral'
 * '<S88>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Tsamp - Ngain'
 * '<S89>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/postSat Signal'
 * '<S90>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/preSat Signal'
 * '<S91>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Anti-windup/Disc. Clamping Parallel'
 * '<S92>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Anti-windup/Disc. Clamping Parallel/Dead Zone'
 * '<S93>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Anti-windup/Disc. Clamping Parallel/Dead Zone/Enabled'
 * '<S94>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/D Gain/Disabled'
 * '<S95>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Filter/Disabled'
 * '<S96>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Filter ICs/Disabled'
 * '<S97>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/I Gain/Internal Parameters'
 * '<S98>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Ideal P Gain/Passthrough'
 * '<S99>'  : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Ideal P Gain Fdbk/Disabled'
 * '<S100>' : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Integrator/Discrete'
 * '<S101>' : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Integrator ICs/Internal IC'
 * '<S102>' : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/N Copy/Disabled wSignal Specification'
 * '<S103>' : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/N Gain/Disabled'
 * '<S104>' : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/P Copy/Disabled'
 * '<S105>' : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Parallel P Gain/Internal Parameters'
 * '<S106>' : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Reset Signal/External Reset'
 * '<S107>' : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Saturation/Enabled'
 * '<S108>' : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Saturation Fdbk/Disabled'
 * '<S109>' : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Sum/Sum_PI'
 * '<S110>' : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Sum Fdbk/Disabled'
 * '<S111>' : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Tracking Mode/Disabled'
 * '<S112>' : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Tracking Mode Sum/Passthrough'
 * '<S113>' : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Tsamp - Integral/Passthrough'
 * '<S114>' : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/Tsamp - Ngain/Passthrough'
 * '<S115>' : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/postSat Signal/Forward_Path'
 * '<S116>' : 'ModelWithControllersOnly/Longitudinal Controller Stanley/PI Controller/Reverse/PI Reverse/preSat Signal/Forward_Path'
 * '<S117>' : 'ModelWithControllersOnly/Subsystem3/Subsystem2'
 * '<S118>' : 'ModelWithControllersOnly/Subsystem3/Subsystem2/Subsystem'
 * '<S119>' : 'ModelWithControllersOnly/Subsystem3/Subsystem2/Subsystem1'
 * '<S120>' : 'ModelWithControllersOnly/Vehicle Dynamics1/Info2ActorState'
 * '<S121>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics'
 * '<S122>' : 'ModelWithControllersOnly/Vehicle Dynamics1/Info2ActorState/MATLAB Function'
 * '<S123>' : 'ModelWithControllersOnly/Vehicle Dynamics1/Info2ActorState/MATLAB Function1'
 * '<S124>' : 'ModelWithControllersOnly/Vehicle Dynamics1/Info2ActorState/MATLAB Function2'
 * '<S125>' : 'ModelWithControllersOnly/Vehicle Dynamics1/Info2ActorState/MATLAB Function3'
 * '<S126>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track'
 * '<S127>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/engine-steer-brake dynamics'
 * '<S128>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Cy'
 * '<S129>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Drag'
 * '<S130>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing'
 * '<S131>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/friction'
 * '<S132>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/front forces'
 * '<S133>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/front steer'
 * '<S134>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/hitch geometry parameters'
 * '<S135>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/rear forces'
 * '<S136>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/rear steer'
 * '<S137>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/sigma'
 * '<S138>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/state'
 * '<S139>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/vehicle model'
 * '<S140>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/wind'
 * '<S141>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Cy/Cy const'
 * '<S142>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Drag/Drag Force'
 * '<S143>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Drag/inertial2body'
 * '<S144>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing'
 * '<S145>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Forces 3DOF'
 * '<S146>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF'
 * '<S147>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Moments'
 * '<S148>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Power'
 * '<S149>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/state2bus'
 * '<S150>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front'
 * '<S151>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric'
 * '<S152>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear'
 * '<S153>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hitch Coordinate Transform'
 * '<S154>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/Rotation Angles to Direction Cosine Matrix'
 * '<S155>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/transform to Inertial axes'
 * '<S156>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/transform to Inertial axes1'
 * '<S157>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/wxR'
 * '<S158>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S159>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/wxR/Subsystem'
 * '<S160>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/wxR/Subsystem1'
 * '<S161>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta'
 * '<S162>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip'
 * '<S163>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Rotation Angles to Direction Cosine Matrix'
 * '<S164>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/transform to Inertial axes'
 * '<S165>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/transform to Inertial axes1'
 * '<S166>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/wxR'
 * '<S167>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip/div0protect - abs poly'
 * '<S168>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip/div0protect - abs poly/Compare To Constant'
 * '<S169>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip/div0protect - abs poly/Compare To Constant1'
 * '<S170>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S171>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/wxR/Subsystem'
 * '<S172>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/wxR/Subsystem1'
 * '<S173>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/Rotation Angles to Direction Cosine Matrix'
 * '<S174>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/transform to Inertial axes'
 * '<S175>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/transform to Inertial axes1'
 * '<S176>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/wxR'
 * '<S177>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S178>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/wxR/Subsystem'
 * '<S179>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/wxR/Subsystem1'
 * '<S180>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hitch Coordinate Transform/Hard Point Coordinate Transform External Displacement'
 * '<S181>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hitch Coordinate Transform/Hard Point Coordinate Transform External Displacement/Rotation Angles to Direction Cosine Matrix'
 * '<S182>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hitch Coordinate Transform/Hard Point Coordinate Transform External Displacement/transform to Inertial axes'
 * '<S183>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hitch Coordinate Transform/Hard Point Coordinate Transform External Displacement/transform to Inertial axes1'
 * '<S184>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hitch Coordinate Transform/Hard Point Coordinate Transform External Displacement/wxR'
 * '<S185>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hitch Coordinate Transform/Hard Point Coordinate Transform External Displacement/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S186>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hitch Coordinate Transform/Hard Point Coordinate Transform External Displacement/wxR/Subsystem'
 * '<S187>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Lateral 3DOF/Hitch Coordinate Transform/Hard Point Coordinate Transform External Displacement/wxR/Subsystem1'
 * '<S188>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Power/Power Accounting Bus Creator'
 * '<S189>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Power/xdot mode'
 * '<S190>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Power/Power Accounting Bus Creator/PwrNotTrnsfrd Input'
 * '<S191>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Power/Power Accounting Bus Creator/PwrStored Input'
 * '<S192>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Power/Power Accounting Bus Creator/PwrTrnsfrd Input'
 * '<S193>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/Power/xdot mode/FxIn'
 * '<S194>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/state2bus/Angle Wrap'
 * '<S195>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/state2bus/Body Slip'
 * '<S196>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/state2bus/COMB2I'
 * '<S197>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/state2bus/Rotation Angles to Direction Cosine Matrix'
 * '<S198>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/state2bus/xddot2ax'
 * '<S199>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/state2bus/Angle Wrap/None'
 * '<S200>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/state2bus/Body Slip/div0protect - abs poly'
 * '<S201>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/state2bus/Body Slip/div0protect - abs poly/Compare To Constant'
 * '<S202>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/state2bus/Body Slip/div0protect - abs poly/Compare To Constant1'
 * '<S203>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/state2bus/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S204>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/Signal Routing/Signal Routing/state2bus/xddot2ax/m^22gn'
 * '<S205>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/friction/mu int'
 * '<S206>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/front forces/ext long'
 * '<S207>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/front steer/delta ext'
 * '<S208>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/hitch geometry parameters/hitch inactive'
 * '<S209>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/rear forces/ext long'
 * '<S210>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/rear steer/delta int'
 * '<S211>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/sigma/sigma'
 * '<S212>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/sigma/sigma/relaxation front'
 * '<S213>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/sigma/sigma/relaxation rear'
 * '<S214>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/state/xdot int'
 * '<S215>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/Vehicle Body 3DOF Single Track/wind/wind int'
 * '<S216>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/engine-steer-brake dynamics/Mapped Steering'
 * '<S217>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/engine-steer-brake dynamics/compute net torque'
 * '<S218>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/engine-steer-brake dynamics/Mapped Steering/LookupGain'
 * '<S219>' : 'ModelWithControllersOnly/Vehicle Dynamics1/VehicleDynamics/engine-steer-brake dynamics/Mapped Steering/Mapped'
 */
#endif                              /* RTW_HEADER_ModelWithControllersOnly_h_ */
