/*
 *  rtmodel.cpp:
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

#include "rtmodel.h"

/* Use this function only if you need to maintain compatibility with an existing static main program. */
void speedgoat_target_model_2021b_step(int_T tid)
{
  switch (tid) {
   case 0 :
    speedgoat_target_model_2021b_step0();
    break;

   case 1 :
    speedgoat_target_model_2021b_step1();
    break;

   default :
    /* do nothing */
    break;
  }
}
