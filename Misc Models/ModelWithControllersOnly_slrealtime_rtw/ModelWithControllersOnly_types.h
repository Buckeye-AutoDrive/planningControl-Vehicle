/*
 * ModelWithControllersOnly_types.h
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

#ifndef RTW_HEADER_ModelWithControllersOnly_types_h_
#define RTW_HEADER_ModelWithControllersOnly_types_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#include "zero_crossing_types.h"

/* Model Code Variants */
#ifndef DEFINED_TYPEDEF_FOR_BusActorState_
#define DEFINED_TYPEDEF_FOR_BusActorState_

struct BusActorState
{
  real_T ActorID;
  real_T Position[3];
  real_T Velocity[3];
  real_T Acceleration[3];
  real_T Roll;
  real_T Pitch;
  real_T Yaw;
  real_T AngularVelocity[3];
};

#endif

/* Parameters (default storage) */
typedef struct P_ModelWithControllersOnly_T_ P_ModelWithControllersOnly_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_ModelWithControllersO_T RT_MODEL_ModelWithControllers_T;

#endif                        /* RTW_HEADER_ModelWithControllersOnly_types_h_ */
