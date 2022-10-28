/*
 * rt_assert.h
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

#ifndef RTW_HEADER_rt_assert_h_
#define RTW_HEADER_rt_assert_h_

/*=========*
 * Asserts *
 *=========*/
#ifndef utAssert
#if defined(DOASSERTS)
#if !defined(PRINT_ASSERTS)
#include <assert.h>
#define utAssert(exp)                  assert(exp)
#else
#include <stdio.h>

static void _assert(char *statement, char *file, int line)
{
  printf("%s in %s on line %d\n", statement, file, line);
}

#define utAssert(_EX)                  ((_EX) ? (void)0 : _assert(#_EX, __FILE__, __LINE__))
#endif

#else
#define utAssert(exp)                                            /* do nothing */
#endif
#endif
#endif                                 /* RTW_HEADER_rt_assert_h_ */
