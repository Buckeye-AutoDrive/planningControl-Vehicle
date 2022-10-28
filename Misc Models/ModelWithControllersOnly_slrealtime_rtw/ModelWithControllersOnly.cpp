/*
 * ModelWithControllersOnly.cpp
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

#include "rte_ModelWithControllersOnly_parameters.h"
#include "ModelWithControllersOnly.h"
#include "ModelWithControllersOnly_private.h"

const real_T ModelWithControllersOnly_RGND = 0.0;/* real_T ground */

/* Block signals (default storage) */
B_ModelWithControllersOnly_T ModelWithControllersOnly_B;

/* Continuous states */
X_ModelWithControllersOnly_T ModelWithControllersOnly_X;

/* Block states (default storage) */
DW_ModelWithControllersOnly_T ModelWithControllersOnly_DW;

/* Previous zero-crossings (trigger) states */
PrevZCX_ModelWithControllersO_T ModelWithControllersOnl_PrevZCX;

/* Real-time model */
RT_MODEL_ModelWithControllers_T ModelWithControllersOnly_M_ =
  RT_MODEL_ModelWithControllers_T();
RT_MODEL_ModelWithControllers_T *const ModelWithControllersOnly_M =
  &ModelWithControllersOnly_M_;

/* Forward declaration for local functions */
static void ModelW_angleUtilities_wrapTo2Pi(real_T *theta);
static void ModelWithC_automlvehdynftiresat(real_T Ftire_x, real_T Ftire_y,
  real_T b_Fxtire_sat, real_T b_Fytire_sat, real_T *Ftire_xs, real_T *Ftire_ys);

/* n-D Spline interpolation function */
real_T look_SplNBinXZcd(uint32_T numDims, const real_T* u, const
  rt_LUTSplineWork * const SWork)
{
  /*
   *   n-D column-major table lookup operating on real_T with:
   *       - Spline interpolation
   *       - Linear extrapolation
   *       - Binary breakpoint search
   *       - Index search starts at the same place each time
   */
  rt_LUTnWork * const TWork_look = SWork->m_TWork;
  real_T* const fraction = static_cast<real_T*>(TWork_look->m_bpLambda);
  uint32_T* const bpIdx = TWork_look->m_bpIndex;
  const uint32_T* const maxIndex = TWork_look->m_maxIndex;
  uint32_T k;
  for (k = 0U; k < numDims; k++) {
    const real_T* const bpData = ((const real_T * const *)
      TWork_look->m_bpDataSet)[k];
    bpIdx[k] = plook_binx(u[k], bpData, maxIndex[k], &fraction[k]);
  }

  return(intrp_NSplcd(numDims, SWork, 2U));
}

/*
 * Second derivative initialization function for spline
 * for last dimension.
 */
void rt_Spline2Derivd(const real_T *x, const real_T *y, uint32_T n, real_T *u,
                      real_T *y2)
{
  real_T p, qn, sig, un;
  uint32_T n1, i, k;
  n1 = n - 1U;
  y2[0U] = 0.0;
  u[0U] = 0.0;
  for (i = 1U; i < n1; i++) {
    real_T dxm1 = x[i] - x[i - 1U];
    real_T dxp1 = x[i + 1U] - x[i];
    real_T dxpm = dxp1 + dxm1;
    sig = dxm1 / dxpm;
    p = (sig * y2[i - 1U]) + 2.0;
    y2[i] = (sig - 1.0) / p;
    u[i] = ((y[i + 1U] - y[i]) / dxp1) - ((y[i] - y[i - 1U]) / dxm1);
    u[i] = (((6.0 * u[i]) / dxpm) - (sig * u[i - 1U])) / p;
  }

  qn = 0.0;
  un = 0.0;
  y2[n1] = (un - (qn * u[n1 - 1U])) / ((qn * y2[n1 - 1U]) + 1.0);
  for (k = n1; k > 0U; k--) {
    y2[k-1U] = (y2[k-1U] * y2[k]) + u[k-1U];
  }

  return;
}

/* n-D natural spline calculation function */
real_T intrp_NSplcd(uint32_T numDims, const rt_LUTSplineWork * const splWork,
                    uint32_T extrapMethod)
{
  uint32_T il;
  uint32_T iu, k, i;
  real_T h, s, p, smsq, pmsq;

  /* intermediate results work areas "this" and "next" */
  const rt_LUTnWork *TWork_interp = static_cast<const rt_LUTnWork *>
    (splWork->m_TWork);
  const real_T *fraction = static_cast<real_T *>(TWork_interp->m_bpLambda);
  const real_T *yp = static_cast<real_T *>(TWork_interp->m_tableData);
  real_T *yyA = static_cast<real_T *>(splWork->m_yyA);
  real_T *yyB = static_cast<real_T *>(splWork->m_yyB);
  real_T *yy2 = static_cast<real_T *>(splWork->m_yy2);
  real_T *up = static_cast<real_T *>(splWork->m_up);
  real_T *y2 = static_cast<real_T *>(splWork->m_y2);
  uint8_T* reCalc = splWork->m_reCalc;
  real_T *dp = static_cast<real_T *>(splWork->m_preBp0AndTable);
  const real_T **bpDataSet = (const real_T **) TWork_interp->m_bpDataSet;
  const real_T *xp = bpDataSet[0U];
  real_T *yy = yyA;
  uint32_T bufBank = 0U;
  uint32_T len = TWork_interp->m_maxIndex[0U] + 1U;

  /* compare bp0 and table to see whether they get changed */
  {
    /* compare the bp0 data */
    if (std::memcmp(dp, xp,
                    len * sizeof(real_T)) != 0) {
      *reCalc = 1;
      (void) std::memcpy(dp, xp,
                         len * sizeof(real_T));
    }

    /* compare the table data */
    dp = &(dp[len]);
    if (std::memcmp(dp, yp,
                    len * splWork->m_numYWorkElts[0U] * sizeof(real_T)) != 0) {
      *reCalc = 1;
      (void) std::memcpy(dp, yp,
                         len * splWork->m_numYWorkElts[0U] * sizeof(real_T));
    }
  }

  if (*reCalc == 1) {
    /* If table and bps are tunable calculate 1st dim 2nd deriv */
    /* Generate first dimension's second derivatives */
    for (i = 0U; i < splWork->m_numYWorkElts[0U]; i++) {
      rt_Spline2Derivd(xp, yp, len, up, y2);
      yp = &yp[len];
      y2 = &y2[len];
    }

    /* Set pointers back to beginning */
    yp = (const real_T *) TWork_interp->m_tableData;
    y2 = (real_T *) splWork->m_y2;
  }

  *reCalc = 0;

  /* Generate at-point splines in each dimension */
  for (k = 0U; k < numDims; k++ ) {
    /* this dimension's input setup */
    xp = bpDataSet[k];
    len = TWork_interp->m_maxIndex[k] + 1U;
    il = TWork_interp->m_bpIndex[k];
    iu = il + 1U;
    h = xp[iu] - xp[il];
    p = fraction[k];
    s = 1.0 - p;
    pmsq = p * ((p*p) - 1.0);
    smsq = s * ((s*s) - 1.0);

    /*
     * Calculate spline curves for input in this
     * dimension at each value of the higher
     * other dimensions\' points in the table.
     */
    if ((p > 1.0) && (extrapMethod == 2U) ) {
      real_T slope;
      for (i = 0U; i < splWork->m_numYWorkElts[k]; i++) {
        slope = (yp[iu] - yp[il]) + ((y2[il]*h*h)*(1.0/6.0));
        yy[i] = yp[iu] + (slope * (p-1.0));
        yp = &yp[len];
        y2 = &y2[len];
      }
    } else if ((p < 0.0) && (extrapMethod == 2U) ) {
      real_T slope;
      for (i = 0U; i < splWork->m_numYWorkElts[k]; i++) {
        slope = (yp[iu] - yp[il]) - ((y2[iu]*h*h)*(1.0/6.0));
        yy[i] = yp[il] + (slope * p);
        yp = &yp[len];
        y2 = &y2[len];
      }
    } else {
      for (i = 0U; i < splWork->m_numYWorkElts[k]; i++) {
        yy[i] = yp[il] + p * (yp[iu] - yp[il]) +
          ((smsq * y2[il] + pmsq * y2[iu])*h*h)*(1.0/6.0);
        yp = &yp[len];
        y2 = &y2[len];
      }
    }

    /* set pointers to new result and calculate second derivatives */
    yp = yy;
    y2 = yy2;
    if (splWork->m_numYWorkElts[k+1U] > 0U ) {
      uint32_T nextLen = TWork_interp->m_maxIndex[k+1U] + 1U;
      const real_T *nextXp = bpDataSet[k+1U];
      for (i = 0U; i < splWork->m_numYWorkElts[k+1U]; i++) {
        rt_Spline2Derivd(nextXp, yp, nextLen, up, y2);
        yp = &yp[nextLen];
        y2 = &y2[nextLen];
      }
    }

    /*
     * Set work vectors yp, y2 and yy for next iteration;
     * the yy just calculated becomes the yp in the
     * next iteration, y2 was just calculated for these
     * new points and the yy buffer is swapped to the space
     * for storing the next iteration\'s results.
     */
    yp = yy;
    y2 = yy2;

    /*
     * Swap buffers for next dimension and
     * toggle bufBank for next iteration.
     */
    if (bufBank == 0U) {
      yy = yyA;
      bufBank = 1U;
    } else {
      yy = yyB;
      bufBank = 0U;
    }
  }

  return( yp[0U] );
}

real_T look1_binlcpw(real_T u0, const real_T bp0[], const real_T table[],
                     uint32_T maxIndex)
{
  real_T frac;
  real_T yL_0d0;
  uint32_T iLeft;

  /* Column-major Lookup 1-D
     Search method: 'binary'
     Use previous index: 'off'
     Interpolation method: 'Linear point-slope'
     Extrapolation method: 'Clip'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Clip'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u0 <= bp0[0U]) {
    iLeft = 0U;
    frac = 0.0;
  } else if (u0 < bp0[maxIndex]) {
    uint32_T bpIdx;
    uint32_T iRght;

    /* Binary Search */
    bpIdx = maxIndex >> 1U;
    iLeft = 0U;
    iRght = maxIndex;
    while (iRght - iLeft > 1U) {
      if (u0 < bp0[bpIdx]) {
        iRght = bpIdx;
      } else {
        iLeft = bpIdx;
      }

      bpIdx = (iRght + iLeft) >> 1U;
    }

    frac = (u0 - bp0[iLeft]) / (bp0[iLeft + 1U] - bp0[iLeft]);
  } else {
    iLeft = maxIndex - 1U;
    frac = 1.0;
  }

  /* Column-major Interpolation 1-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  yL_0d0 = table[iLeft];
  return (table[iLeft + 1U] - yL_0d0) * frac + yL_0d0;
}

real_T look1_binlxpw(real_T u0, const real_T bp0[], const real_T table[],
                     uint32_T maxIndex)
{
  real_T frac;
  real_T yL_0d0;
  uint32_T iLeft;

  /* Column-major Lookup 1-D
     Search method: 'binary'
     Use previous index: 'off'
     Interpolation method: 'Linear point-slope'
     Extrapolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u0 <= bp0[0U]) {
    iLeft = 0U;
    frac = (u0 - bp0[0U]) / (bp0[1U] - bp0[0U]);
  } else if (u0 < bp0[maxIndex]) {
    uint32_T bpIdx;
    uint32_T iRght;

    /* Binary Search */
    bpIdx = maxIndex >> 1U;
    iLeft = 0U;
    iRght = maxIndex;
    while (iRght - iLeft > 1U) {
      if (u0 < bp0[bpIdx]) {
        iRght = bpIdx;
      } else {
        iLeft = bpIdx;
      }

      bpIdx = (iRght + iLeft) >> 1U;
    }

    frac = (u0 - bp0[iLeft]) / (bp0[iLeft + 1U] - bp0[iLeft]);
  } else {
    iLeft = maxIndex - 1U;
    frac = (u0 - bp0[maxIndex - 1U]) / (bp0[maxIndex] - bp0[maxIndex - 1U]);
  }

  /* Column-major Interpolation 1-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  yL_0d0 = table[iLeft];
  return (table[iLeft + 1U] - yL_0d0) * frac + yL_0d0;
}

uint32_T plook_binx(real_T u, const real_T bp[], uint32_T maxIndex, real_T
                    *fraction)
{
  uint32_T bpIndex;

  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u <= bp[0U]) {
    bpIndex = 0U;
    *fraction = (u - bp[0U]) / (bp[1U] - bp[0U]);
  } else if (u < bp[maxIndex]) {
    bpIndex = binsearch_u32d(u, bp, maxIndex >> 1U, maxIndex);
    *fraction = (u - bp[bpIndex]) / (bp[bpIndex + 1U] - bp[bpIndex]);
  } else {
    bpIndex = maxIndex - 1U;
    *fraction = (u - bp[maxIndex - 1U]) / (bp[maxIndex] - bp[maxIndex - 1U]);
  }

  return bpIndex;
}

uint32_T binsearch_u32d(real_T u, const real_T bp[], uint32_T startIndex,
  uint32_T maxIndex)
{
  uint32_T bpIdx;
  uint32_T bpIndex;
  uint32_T iRght;

  /* Binary Search */
  bpIdx = startIndex;
  bpIndex = 0U;
  iRght = maxIndex;
  while (iRght - bpIndex > 1U) {
    if (u < bp[bpIdx]) {
      iRght = bpIdx;
    } else {
      bpIndex = bpIdx;
    }

    bpIdx = (iRght + bpIndex) >> 1U;
  }

  return bpIndex;
}

/*
 *         This function updates active task flag for each subrate.
 *         The function is called in the model base rate function.
 *         It maintains SampleHit information to allow scheduling
 *         of the subrates from the base rate function.
 */
void rate_scheduler(void)
{
  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (ModelWithControllersOnly_M->Timing.TaskCounters.TID[2])++;
  if ((ModelWithControllersOnly_M->Timing.TaskCounters.TID[2]) > 4) {/* Sample time: [0.05s, 0.0s] */
    ModelWithControllersOnly_M->Timing.TaskCounters.TID[2] = 0;
  }
}

/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static const real_T rt_ODE3_A[3] = {
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3] = {
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE3_IntgData *id = static_cast<ODE3_IntgData *>(rtsiGetSolverData(si));
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T hB[3];
  int_T i;
  int_T nXc = 10;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) std::memcpy(y, x,
                     static_cast<uint_T>(nXc)*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  ModelWithControllersOnly_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  ModelWithControllersOnly_step0();
  ModelWithControllersOnly_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  ModelWithControllersOnly_step0();
  ModelWithControllersOnly_derivatives();

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Function for MATLAB Function: '<S1>/Kinematic' */
static void ModelW_angleUtilities_wrapTo2Pi(real_T *theta)
{
  real_T x;
  boolean_T positiveInput;
  positiveInput = (*theta > 0.0);
  x = *theta;
  if (rtIsNaN(*theta) || rtIsInf(*theta)) {
    *theta = (rtNaN);
  } else if (*theta == 0.0) {
    *theta = 0.0;
  } else {
    boolean_T rEQ0;
    *theta = std::fmod(*theta, 6.2831853071795862);
    rEQ0 = (*theta == 0.0);
    if (!rEQ0) {
      real_T q;
      q = std::abs(x / 6.2831853071795862);
      rEQ0 = !(std::abs(q - std::floor(q + 0.5)) > 2.2204460492503131E-16 * q);
    }

    if (rEQ0) {
      *theta = 0.0;
    } else if (x < 0.0) {
      *theta += 6.2831853071795862;
    }
  }

  *theta += static_cast<real_T>((*theta == 0.0) && positiveInput) *
    6.2831853071795862;
}

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    int32_T tmp;
    int32_T tmp_0;
    if (u1 > 0.0) {
      tmp = 1;
    } else {
      tmp = -1;
    }

    if (u0 > 0.0) {
      tmp_0 = 1;
    } else {
      tmp_0 = -1;
    }

    y = std::atan2(static_cast<real_T>(tmp_0), static_cast<real_T>(tmp));
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = std::atan2(u0, u1);
  }

  return y;
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    real_T tmp;
    real_T tmp_0;
    tmp = std::abs(u0);
    tmp_0 = std::abs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = 1.0;
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = (rtNaN);
    } else {
      y = std::pow(u0, u1);
    }
  }

  return y;
}

/* Function for MATLAB Function: '<S126>/vehicle model' */
static void ModelWithC_automlvehdynftiresat(real_T Ftire_x, real_T Ftire_y,
  real_T b_Fxtire_sat, real_T b_Fytire_sat, real_T *Ftire_xs, real_T *Ftire_ys)
{
  real_T Ftire_mag;
  real_T Ftire_x_max;
  real_T theta_Ftire;
  theta_Ftire = rt_atan2d_snf(Ftire_x, Ftire_y);
  Ftire_x_max = b_Fxtire_sat * std::cos(theta_Ftire);
  Ftire_mag = b_Fytire_sat * std::sin(theta_Ftire);
  Ftire_mag = b_Fxtire_sat * b_Fytire_sat / std::sqrt(Ftire_x_max * Ftire_x_max
    + Ftire_mag * Ftire_mag);
  Ftire_x_max = Ftire_mag * std::sin(theta_Ftire);
  theta_Ftire = Ftire_mag * std::cos(theta_Ftire);
  *Ftire_xs = Ftire_x;
  if (std::abs(Ftire_x) > std::abs(Ftire_x_max)) {
    *Ftire_xs = Ftire_x_max;
  }

  *Ftire_ys = Ftire_y;
  if (std::abs(Ftire_y) > std::abs(theta_Ftire)) {
    *Ftire_ys = theta_Ftire;
  }
}

/* Model step function for TID0 */
void ModelWithControllersOnly_step0(void) /* Sample time: [0.0s, 0.0s] */
{
  if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
    /* set solver stop time */
    if (!(ModelWithControllersOnly_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&ModelWithControllersOnly_M->solverInfo,
                            ((ModelWithControllersOnly_M->Timing.clockTickH0 + 1)
        * ModelWithControllersOnly_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&ModelWithControllersOnly_M->solverInfo,
                            ((ModelWithControllersOnly_M->Timing.clockTick0 + 1)
        * ModelWithControllersOnly_M->Timing.stepSize0 +
        ModelWithControllersOnly_M->Timing.clockTickH0 *
        ModelWithControllersOnly_M->Timing.stepSize0 * 4294967296.0));
    }

    /* Update the flag to indicate when data transfers from
     *  Sample time: [0.01s, 0.0s] to Sample time: [0.05s, 0.0s]  */
    (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2)++;
    if ((ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2) > 4) {
      ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 = 0;
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(ModelWithControllersOnly_M)) {
    ModelWithControllersOnly_M->Timing.t[0] = rtsiGetT
      (&ModelWithControllersOnly_M->solverInfo);
  }

  {
    real_T poseF[3];
    real_T refPose[3];
    real_T B1;
    real_T B2;
    real_T Cy_f;
    real_T Cy_r;
    real_T Fx_f;
    real_T Fx_ft;
    real_T Fx_r;
    real_T Fx_rt;
    real_T Fy_ft;
    real_T Fy_rt;
    real_T FzCalc_idx_0;
    real_T FzCalc_idx_1;
    real_T Nf;
    real_T Nr;
    real_T alfa_f;
    real_T alfa_r;
    real_T b_B1;
    real_T b_B2;
    real_T b_Fxtire_sat;
    real_T b_Fytire_sat;
    real_T b_Fznom;
    real_T b_Izz;
    real_T b_b;
    real_T b_g;
    real_T b_h;
    real_T b_m;
    real_T b_y;
    real_T c_a;
    real_T d_a;
    real_T d_idx_0;
    real_T d_idx_1;
    real_T maxFerr;
    real_T posError;
    real_T r;
    real_T rdot;
    real_T tmp_5;
    real_T tmp_6;
    real_T tmp_7;
    real_T tmp_8;
    real_T tmp_9;
    real_T tmp_a;
    real_T tmp_b;
    real_T tmp_c;
    real_T u;
    real_T u2;
    real_T x_data;
    real_T xddot;
    real_T xdot;
    real_T yddot;
    real_T ydot;
    real_T *lastU;
    real_T *tmp;
    real_T *tmp_0;
    real_T *tmp_1;
    real_T *tmp_2;
    real_T *tmp_3;
    real_T *tmp_4;
    int32_T i;
    int32_T loop_ub;
    int8_T wrBufIdx;
    boolean_T zcEvent;
    ZCEventType zcEvent_0;
    tmp_c = *get_MAX_TORQUE();
    tmp_b = *get_MIN_TORQUE();
    tmp_a = *get_MAX_BRAKE();
    tmp_9 = *get_MIN_BRAKE();
    tmp_8 = *get_LF();
    b_b = *get_LR();
    b_m = *get_MCAR();
    tmp_7 = *get_Tdead();
    tmp_6 = *get_Vdead();
    tmp_5 = *get_R();
    tmp_4 = get_NUMD_S();
    tmp_3 = get_DEND_S();
    tmp_2 = get_NUMD_B();
    tmp_1 = get_DEND_B();
    tmp_0 = get_NUMD_T();
    tmp = get_DEND_T();
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Constant: '<S126>/xdot_oConstant' incorporates:
       *  Concatenate: '<S126>/Vector Concatenate'
       */
      ModelWithControllersOnly_B.VectorConcatenate[0] =
        ModelWithControllersOnly_cal->VehicleDynamics_InitSpeed;

      /* Constant: '<S126>/ydot_oConstant' incorporates:
       *  Concatenate: '<S126>/Vector Concatenate'
       */
      ModelWithControllersOnly_B.VectorConcatenate[1] =
        ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_ydot;

      /* Constant: '<S121>/Constant4' incorporates:
       *  Concatenate: '<S126>/Vector Concatenate'
       */
      ModelWithControllersOnly_B.VectorConcatenate[2] = *get_psi_o();

      /* Constant: '<S126>/r_oConstant' incorporates:
       *  Concatenate: '<S126>/Vector Concatenate'
       */
      ModelWithControllersOnly_B.VectorConcatenate[3] =
        ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_r_o;
    }

    /* Integrator: '<S214>/Integrator' */
    if (ModelWithControllersOnly_DW.Integrator_IWORK != 0) {
      ModelWithControllersOnly_X.Integrator_CSTATE[0] =
        ModelWithControllersOnly_B.VectorConcatenate[0];
      ModelWithControllersOnly_X.Integrator_CSTATE[1] =
        ModelWithControllersOnly_B.VectorConcatenate[1];
      ModelWithControllersOnly_X.Integrator_CSTATE[2] =
        ModelWithControllersOnly_B.VectorConcatenate[2];
      ModelWithControllersOnly_X.Integrator_CSTATE[3] =
        ModelWithControllersOnly_B.VectorConcatenate[3];
    }

    /* Integrator: '<S214>/Integrator' */
    ModelWithControllersOnly_B.Integrator[0] =
      ModelWithControllersOnly_X.Integrator_CSTATE[0];
    ModelWithControllersOnly_B.Integrator[1] =
      ModelWithControllersOnly_X.Integrator_CSTATE[1];
    ModelWithControllersOnly_B.Integrator[2] =
      ModelWithControllersOnly_X.Integrator_CSTATE[2];
    ModelWithControllersOnly_B.Integrator[3] =
      ModelWithControllersOnly_X.Integrator_CSTATE[3];

    /* Saturate: '<S5>/Saturation' */
    u = ModelWithControllersOnly_B.Integrator[0];
    posError = ModelWithControllersOnly_cal->Saturation_LowerSat;
    u2 = ModelWithControllersOnly_P.Saturation_UpperSat;
    if (u > u2) {
      /* Saturate: '<S5>/Saturation' */
      ModelWithControllersOnly_B.Saturation = u2;
    } else if (u < posError) {
      /* Saturate: '<S5>/Saturation' */
      ModelWithControllersOnly_B.Saturation = posError;
    } else {
      /* Saturate: '<S5>/Saturation' */
      ModelWithControllersOnly_B.Saturation = u;
    }

    /* End of Saturate: '<S5>/Saturation' */
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Constant: '<S121>/Constant1' incorporates:
       *  Concatenate: '<S126>/Vector Concatenate3'
       */
      ModelWithControllersOnly_B.VectorConcatenate3[0] = *get_X_o();

      /* Constant: '<S121>/Constant3' incorporates:
       *  Concatenate: '<S126>/Vector Concatenate3'
       */
      ModelWithControllersOnly_B.VectorConcatenate3[1] = *get_Y_o();
    }

    /* Integrator: '<S149>/Integrator' */
    if (ModelWithControllersOnly_DW.Integrator_IWORK_c != 0) {
      ModelWithControllersOnly_X.Integrator_CSTATE_b[0] =
        ModelWithControllersOnly_B.VectorConcatenate3[0];
      ModelWithControllersOnly_X.Integrator_CSTATE_b[1] =
        ModelWithControllersOnly_B.VectorConcatenate3[1];
    }

    /* Integrator: '<S149>/Integrator' */
    ModelWithControllersOnly_B.Integrator_m[0] =
      ModelWithControllersOnly_X.Integrator_CSTATE_b[0];
    ModelWithControllersOnly_B.Integrator_m[1] =
      ModelWithControllersOnly_X.Integrator_CSTATE_b[1];

    /* Gain: '<S3>/Gain1' */
    ModelWithControllersOnly_B.Gain1 = ModelWithControllersOnly_cal->Gain1_Gain *
      ModelWithControllersOnly_B.Integrator[2];

    /* SignalConversion generated from: '<S3>/Transpose' */
    ModelWithControllersOnly_B.TmpSignalConversionAtTransposeI[0] =
      ModelWithControllersOnly_B.Integrator_m[0];
    ModelWithControllersOnly_B.TmpSignalConversionAtTransposeI[1] =
      ModelWithControllersOnly_B.Integrator_m[1];
    ModelWithControllersOnly_B.TmpSignalConversionAtTransposeI[2] =
      ModelWithControllersOnly_B.Gain1;

    /* Math: '<S3>/Transpose' */
    ModelWithControllersOnly_B.Transpose[0] =
      ModelWithControllersOnly_B.TmpSignalConversionAtTransposeI[0];

    /* MATLAB Function: '<S1>/Kinematic' incorporates:
     *  Constant: '<Root>/Constant4'
     *  Math: '<S3>/Transpose'
     */
    refPose[0] = ModelWithControllersOnly_cal->Constant4_Value[0];
    poseF[0] = ModelWithControllersOnly_B.Transpose[0];

    /* Math: '<S3>/Transpose' */
    ModelWithControllersOnly_B.Transpose[1] =
      ModelWithControllersOnly_B.TmpSignalConversionAtTransposeI[1];

    /* MATLAB Function: '<S1>/Kinematic' incorporates:
     *  Constant: '<Root>/Constant4'
     *  Math: '<S3>/Transpose'
     */
    refPose[1] = ModelWithControllersOnly_cal->Constant4_Value[1];
    poseF[1] = ModelWithControllersOnly_B.Transpose[1];

    /* Math: '<S3>/Transpose' */
    ModelWithControllersOnly_B.Transpose[2] =
      ModelWithControllersOnly_B.TmpSignalConversionAtTransposeI[2];

    /* MATLAB Function: '<S1>/Kinematic' incorporates:
     *  Constant: '<Root>/Constant1'
     *  Constant: '<Root>/Constant4'
     *  Math: '<S3>/Transpose'
     */
    if (ModelWithControllersOnly_cal->Constant1_Value_k == 1.0) {
      u = ModelWithControllersOnly_cal->LateralControllerStanley_Positi;
    } else {
      u = ModelWithControllersOnly_cal->LateralControllerStanley_Posi_i;
    }

    refPose[2] = 0.017453292519943295 *
      ModelWithControllersOnly_cal->Constant4_Value[2];
    ModelW_angleUtilities_wrapTo2Pi(&refPose[2]);
    poseF[2] = 0.017453292519943295 * ModelWithControllersOnly_B.Transpose[2];
    ModelW_angleUtilities_wrapTo2Pi(&poseF[2]);
    FzCalc_idx_0 = std::cos(refPose[2]);
    FzCalc_idx_1 = std::sin(refPose[2]);
    if (ModelWithControllersOnly_cal->Constant1_Value_k == 1.0) {
      poseF[0] = ModelWithControllersOnly_cal->Kinematic_Wheelbase * std::cos
        (poseF[2]) + ModelWithControllersOnly_B.Transpose[0];
      poseF[1] = ModelWithControllersOnly_cal->Kinematic_Wheelbase * std::sin
        (poseF[2]) + ModelWithControllersOnly_B.Transpose[1];
      d_idx_0 = poseF[0] - refPose[0];
      d_idx_1 = poseF[1] - refPose[1];
    } else {
      d_idx_0 = poseF[0] - refPose[0];
      d_idx_1 = poseF[1] - refPose[1];
    }

    posError = -(d_idx_0 * FzCalc_idx_1 - FzCalc_idx_0 * d_idx_1);
    u2 = (poseF[2] - refPose[2]) + 3.1415926535897931;
    ModelW_angleUtilities_wrapTo2Pi(&u2);
    if (ModelWithControllersOnly_cal->Constant1_Value_k == 1.0) {
      u = -(std::atan(u * posError / (ModelWithControllersOnly_B.Saturation +
              1.0)) + (u2 - 3.1415926535897931));
    } else {
      u = std::atan(u * posError / (ModelWithControllersOnly_B.Saturation + -1.0))
        + (u2 - 3.1415926535897931);
    }

    ModelWithControllersOnly_B.steerCmd = 57.295779513082323 * u;
    u = ModelWithControllersOnly_B.steerCmd;
    if (u < 0.0) {
      u2 = -1.0;
    } else if (u > 0.0) {
      u2 = 1.0;
    } else if (u == 0.0) {
      u2 = 0.0;
    } else {
      u2 = (rtNaN);
    }

    u = std::abs(ModelWithControllersOnly_B.steerCmd);
    posError = ModelWithControllersOnly_cal->Kinematic_MaxSteeringAngle;
    if ((u <= posError) || rtIsNaN(posError)) {
      posError = u;
    }

    ModelWithControllersOnly_B.steerCmd = u2 * posError;

    /* Product: '<S8>/Product' incorporates:
     *  Constant: '<Root>/Constant'
     */
    ModelWithControllersOnly_B.Product = ModelWithControllersOnly_B.Saturation *
      ModelWithControllersOnly_cal->Constant_Value;

    /* Gain: '<S9>/Gain' */
    ModelWithControllersOnly_B.Gain = ModelWithControllersOnly_cal->Gain_Gain *
      ModelWithControllersOnly_B.Product;

    /* Product: '<S8>/Multiply' */
    ModelWithControllersOnly_B.Multiply = ModelWithControllersOnly_B.Saturation *
      ModelWithControllersOnly_B.Gain;

    /* Gain: '<S8>/Gain1' */
    ModelWithControllersOnly_B.CurvedSteadyYaw =
      ModelWithControllersOnly_cal->Gain1_Gain_o *
      ModelWithControllersOnly_B.Multiply;

    /* Gain: '<Root>/Gain1' */
    ModelWithControllersOnly_B.Gain1_g =
      ModelWithControllersOnly_cal->Gain1_Gain_a *
      ModelWithControllersOnly_B.Integrator[3];

    /* Sum: '<S8>/Minus' */
    ModelWithControllersOnly_B.Minus = ModelWithControllersOnly_B.Gain -
      ModelWithControllersOnly_B.Gain1_g;

    /* Gain: '<S8>/Gain' */
    ModelWithControllersOnly_B.YawRateFeedback =
      ModelWithControllersOnly_cal->LateralControllerStanley_YawRat *
      ModelWithControllersOnly_B.Minus;
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* UnitDelay: '<S8>/Unit Delay' */
      ModelWithControllersOnly_B.UnitDelay =
        ModelWithControllersOnly_DW.UnitDelay_DSTATE;

      /* UnitDelay: '<Root>/Unit Delay' */
      ModelWithControllersOnly_B.UnitDelay_l =
        ModelWithControllersOnly_DW.UnitDelay_DSTATE_e;

      /* Sum: '<S8>/Minus1' */
      ModelWithControllersOnly_B.Minus1 = ModelWithControllersOnly_B.UnitDelay -
        ModelWithControllersOnly_B.UnitDelay_l;

      /* Gain: '<S8>/Gain2' */
      ModelWithControllersOnly_B.SteeringDelay =
        ModelWithControllersOnly_cal->LateralControllerStanley_DelayG *
        ModelWithControllersOnly_B.Minus1;
    }

    /* Sum: '<S8>/Add' */
    ModelWithControllersOnly_B.Add = ((ModelWithControllersOnly_B.steerCmd +
      ModelWithControllersOnly_B.CurvedSteadyYaw) +
      ModelWithControllersOnly_B.YawRateFeedback) +
      ModelWithControllersOnly_B.SteeringDelay;

    /* Saturate: '<S8>/Saturation' */
    u = ModelWithControllersOnly_B.Add;
    posError = ModelWithControllersOnly_cal->Saturation_LowerSat_f;
    u2 = ModelWithControllersOnly_cal->Saturation_UpperSat_n;
    if (u > u2) {
      /* Saturate: '<S8>/Saturation' */
      ModelWithControllersOnly_B.Saturation_k = u2;
    } else if (u < posError) {
      /* Saturate: '<S8>/Saturation' */
      ModelWithControllersOnly_B.Saturation_k = posError;
    } else {
      /* Saturate: '<S8>/Saturation' */
      ModelWithControllersOnly_B.Saturation_k = u;
    }

    /* End of Saturate: '<S8>/Saturation' */

    /* Gain: '<Root>/Gain' */
    posError = 0.0 / *get_STEER_RATIO();

    /* Gain: '<Root>/Gain' */
    ModelWithControllersOnly_B.Gain_n = posError *
      ModelWithControllersOnly_B.Saturation_k;
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* RateTransition generated from: '<S12>/Equal1' incorporates:
       *  Constant: '<Root>/Constant1'
       */
      rtw_slrealtime_mutex_lock
        (ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_d0_SEMAPH);
      wrBufIdx = static_cast<int8_T>
        (ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_LstBufWR + 1);
      if (wrBufIdx == 3) {
        wrBufIdx = 0;
      }

      if (wrBufIdx == ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_RDBuf) {
        wrBufIdx = static_cast<int8_T>(wrBufIdx + 1);
        if (wrBufIdx == 3) {
          wrBufIdx = 0;
        }
      }

      rtw_slrealtime_mutex_unlock
        (ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_d0_SEMAPH);
      switch (wrBufIdx) {
       case 0:
        ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_Buf0 =
          ModelWithControllersOnly_cal->Constant1_Value_k;
        break;

       case 1:
        ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_Buf1 =
          ModelWithControllersOnly_cal->Constant1_Value_k;
        break;

       case 2:
        ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_Buf2 =
          ModelWithControllersOnly_cal->Constant1_Value_k;
        break;
      }

      ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_LstBufWR = wrBufIdx;

      /* End of RateTransition generated from: '<S12>/Equal1' */

      /* RateTransition generated from: '<S12>/Equal2' incorporates:
       *  Constant: '<Root>/Constant1'
       */
      rtw_slrealtime_mutex_lock
        (ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_d0_SEMAPH);
      wrBufIdx = static_cast<int8_T>
        (ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_LstBufWR + 1);
      if (wrBufIdx == 3) {
        wrBufIdx = 0;
      }

      if (wrBufIdx == ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_RDBuf) {
        wrBufIdx = static_cast<int8_T>(wrBufIdx + 1);
        if (wrBufIdx == 3) {
          wrBufIdx = 0;
        }
      }

      rtw_slrealtime_mutex_unlock
        (ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_d0_SEMAPH);
      switch (wrBufIdx) {
       case 0:
        ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_Buf0 =
          ModelWithControllersOnly_cal->Constant1_Value_k;
        break;

       case 1:
        ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_Buf1 =
          ModelWithControllersOnly_cal->Constant1_Value_k;
        break;

       case 2:
        ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_Buf2 =
          ModelWithControllersOnly_cal->Constant1_Value_k;
        break;
      }

      ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_LstBufWR = wrBufIdx;

      /* End of RateTransition generated from: '<S12>/Equal2' */

      /* RateTransition generated from: '<S12>/Equal4' incorporates:
       *  Constant: '<Root>/Constant1'
       */
      rtw_slrealtime_mutex_lock
        (ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_d0_SEMAPH);
      wrBufIdx = static_cast<int8_T>
        (ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_LstBufWR + 1);
      if (wrBufIdx == 3) {
        wrBufIdx = 0;
      }

      if (wrBufIdx == ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_RDBuf) {
        wrBufIdx = static_cast<int8_T>(wrBufIdx + 1);
        if (wrBufIdx == 3) {
          wrBufIdx = 0;
        }
      }

      rtw_slrealtime_mutex_unlock
        (ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_d0_SEMAPH);
      switch (wrBufIdx) {
       case 0:
        ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_Buf0 =
          ModelWithControllersOnly_cal->Constant1_Value_k;
        break;

       case 1:
        ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_Buf1 =
          ModelWithControllersOnly_cal->Constant1_Value_k;
        break;

       case 2:
        ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_Buf2 =
          ModelWithControllersOnly_cal->Constant1_Value_k;
        break;
      }

      ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_LstBufWR = wrBufIdx;

      /* End of RateTransition generated from: '<S12>/Equal4' */

      /* RateTransition generated from: '<S12>/Equal5' incorporates:
       *  Constant: '<Root>/Constant1'
       */
      rtw_slrealtime_mutex_lock
        (ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_d0_SEMAPH);
      wrBufIdx = static_cast<int8_T>
        (ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_LstBufWR + 1);
      if (wrBufIdx == 3) {
        wrBufIdx = 0;
      }

      if (wrBufIdx == ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_RDBuf) {
        wrBufIdx = static_cast<int8_T>(wrBufIdx + 1);
        if (wrBufIdx == 3) {
          wrBufIdx = 0;
        }
      }

      rtw_slrealtime_mutex_unlock
        (ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_d0_SEMAPH);
      switch (wrBufIdx) {
       case 0:
        ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_Buf0 =
          ModelWithControllersOnly_cal->Constant1_Value_k;
        break;

       case 1:
        ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_Buf1 =
          ModelWithControllersOnly_cal->Constant1_Value_k;
        break;

       case 2:
        ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_Buf2 =
          ModelWithControllersOnly_cal->Constant1_Value_k;
        break;
      }

      ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_LstBufWR = wrBufIdx;

      /* End of RateTransition generated from: '<S12>/Equal5' */

      /* RateTransition generated from: '<S12>/Equal3' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.xdot_WrBufIdx = static_cast<int8_T>
          (ModelWithControllersOnly_DW.xdot_WrBufIdx == 0);
      }

      ModelWithControllersOnly_DW.xdot_Buf[ModelWithControllersOnly_DW.xdot_WrBufIdx]
        = ModelWithControllersOnly_B.Saturation;

      /* End of RateTransition generated from: '<S12>/Equal3' */

      /* RateTransition generated from: '<S12>/Sign1' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.xdot1_WrBufIdx = static_cast<int8_T>
          (ModelWithControllersOnly_DW.xdot1_WrBufIdx == 0);
      }

      ModelWithControllersOnly_DW.xdot1_Buf[ModelWithControllersOnly_DW.xdot1_WrBufIdx]
        = ModelWithControllersOnly_B.Saturation;

      /* End of RateTransition generated from: '<S12>/Sign1' */

      /* RateTransition generated from: '<S10>/Multiply' incorporates:
       *  Constant: '<Root>/Constant1'
       */
      rtw_slrealtime_mutex_lock
        (ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_d0_SEMA);
      wrBufIdx = static_cast<int8_T>
        (ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_LstBufW + 1);
      if (wrBufIdx == 3) {
        wrBufIdx = 0;
      }

      if (wrBufIdx == ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_RDBuf)
      {
        wrBufIdx = static_cast<int8_T>(wrBufIdx + 1);
        if (wrBufIdx == 3) {
          wrBufIdx = 0;
        }
      }

      rtw_slrealtime_mutex_unlock
        (ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_d0_SEMA);
      switch (wrBufIdx) {
       case 0:
        ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_Buf0 =
          ModelWithControllersOnly_cal->Constant1_Value_k;
        break;

       case 1:
        ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_Buf1 =
          ModelWithControllersOnly_cal->Constant1_Value_k;
        break;

       case 2:
        ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_Buf2 =
          ModelWithControllersOnly_cal->Constant1_Value_k;
        break;
      }

      ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_LstBufW = wrBufIdx;

      /* End of RateTransition generated from: '<S10>/Multiply' */

      /* RateTransition generated from: '<S11>/Switch Case' incorporates:
       *  Constant: '<Root>/Constant1'
       */
      rtw_slrealtime_mutex_lock
        (ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_d0_SE);
      wrBufIdx = static_cast<int8_T>
        (ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_LstBu + 1);
      if (wrBufIdx == 3) {
        wrBufIdx = 0;
      }

      if (wrBufIdx ==
          ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_RDBuf) {
        wrBufIdx = static_cast<int8_T>(wrBufIdx + 1);
        if (wrBufIdx == 3) {
          wrBufIdx = 0;
        }
      }

      rtw_slrealtime_mutex_unlock
        (ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_d0_SE);
      switch (wrBufIdx) {
       case 0:
        ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_Buf0 =
          ModelWithControllersOnly_cal->Constant1_Value_k;
        break;

       case 1:
        ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_Buf1 =
          ModelWithControllersOnly_cal->Constant1_Value_k;
        break;

       case 2:
        ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_Buf2 =
          ModelWithControllersOnly_cal->Constant1_Value_k;
        break;
      }

      ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_LstBu = wrBufIdx;

      /* End of RateTransition generated from: '<S11>/Switch Case' */

      /* RateTransition generated from: '<S2>/Minus' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.xdot_WrBufIdx_b = static_cast<int8_T>
          (ModelWithControllersOnly_DW.xdot_WrBufIdx_b == 0);
      }

      ModelWithControllersOnly_DW.xdot_Buf_d[ModelWithControllersOnly_DW.xdot_WrBufIdx_b]
        = ModelWithControllersOnly_B.Saturation;

      /* End of RateTransition generated from: '<S2>/Minus' */

      /* RateTransition generated from: '<S119>/Integrator' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_RdBuf =
          static_cast<int8_T>
          (ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_RdBuf == 0);
      }

      /* RateTransition generated from: '<S119>/Integrator' */
      ModelWithControllersOnly_B.TmpRTBAtIntegratorInport2 =
        ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_Buf[ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_RdBuf];
    }

    /* Integrator: '<S119>/Integrator' */
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      zcEvent_0 = rt_ZCFcn(FALLING_ZERO_CROSSING,
                           &ModelWithControllersOnl_PrevZCX.Integrator_Reset_ZCE,
                           (ModelWithControllersOnly_B.TmpRTBAtIntegratorInport2));
      zcEvent = (zcEvent_0 != NO_ZCEVENT);

      /* evaluate zero-crossings */
      if (zcEvent) {
        ModelWithControllersOnly_X.Integrator_CSTATE_o =
          ModelWithControllersOnly_cal->Integrator_IC;
      }
    }

    /* Integrator: '<S119>/Integrator' */
    ModelWithControllersOnly_B.v =
      ModelWithControllersOnly_X.Integrator_CSTATE_o;

    /* Math: '<S119>/Square' */
    ModelWithControllersOnly_B.Square = ModelWithControllersOnly_B.v *
      ModelWithControllersOnly_B.v;

    /* Gain: '<S119>/Gain1' */
    ModelWithControllersOnly_B.F_aero =
      ModelWithControllersOnly_cal->Gain1_Gain_d *
      ModelWithControllersOnly_B.Square;
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Gain: '<S119>/Gain2' incorporates:
       *  Constant: '<S119>/Constant1'
       */
      ModelWithControllersOnly_B.Gain2 =
        ModelWithControllersOnly_cal->Gain2_Gain * b_m;

      /* RateTransition generated from: '<S119>/Sum' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.TmpRTBAtSumInport3_RdBufIdx =
          static_cast<int8_T>
          (ModelWithControllersOnly_DW.TmpRTBAtSumInport3_RdBufIdx == 0);
      }

      /* RateTransition generated from: '<S119>/Sum' */
      ModelWithControllersOnly_B.TmpRTBAtSumInport3 =
        ModelWithControllersOnly_DW.TmpRTBAtSumInport3_Buf[ModelWithControllersOnly_DW.TmpRTBAtSumInport3_RdBufIdx];
    }

    /* Sum: '<S119>/Sum' */
    ModelWithControllersOnly_B.Sum = (ModelWithControllersOnly_B.F_aero +
      ModelWithControllersOnly_B.Gain2) +
      ModelWithControllersOnly_B.TmpRTBAtSumInport3;

    /* Gain: '<S119>/Gain4' */
    posError = -tmp_5 / 5.6400000000000006;

    /* Gain: '<S119>/Gain4' */
    ModelWithControllersOnly_B.T_eng = posError * ModelWithControllersOnly_B.Sum;

    /* Saturate: '<S4>/Brake Saturation' */
    u = ModelWithControllersOnly_B.T_eng;
    posError = tmp_9;
    u2 = tmp_a;
    if (u > u2) {
      /* Saturate: '<S4>/Brake Saturation' */
      ModelWithControllersOnly_B.BrakeSaturation = u2;
    } else if (u < posError) {
      /* Saturate: '<S4>/Brake Saturation' */
      ModelWithControllersOnly_B.BrakeSaturation = posError;
    } else {
      /* Saturate: '<S4>/Brake Saturation' */
      ModelWithControllersOnly_B.BrakeSaturation = u;
    }

    /* End of Saturate: '<S4>/Brake Saturation' */

    /* RateTransition generated from: '<S118>/Integrator' */
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_RdB_d =
          static_cast<int8_T>
          (ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_RdB_d == 0);
      }

      /* RateTransition generated from: '<S118>/Integrator' */
      ModelWithControllersOnly_B.TmpRTBAtIntegratorInport2_e =
        ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_Buf_m[ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_RdB_d];
    }

    /* End of RateTransition generated from: '<S118>/Integrator' */

    /* Integrator: '<S118>/Integrator' */
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      zcEvent_0 = rt_ZCFcn(FALLING_ZERO_CROSSING,
                           &ModelWithControllersOnl_PrevZCX.Integrator_Reset_ZCE_c,
                           (ModelWithControllersOnly_B.TmpRTBAtIntegratorInport2_e));
      zcEvent = (zcEvent_0 != NO_ZCEVENT);

      /* evaluate zero-crossings */
      if (zcEvent) {
        ModelWithControllersOnly_X.Integrator_CSTATE_c =
          ModelWithControllersOnly_cal->Integrator_IC_g;
      }
    }

    /* Integrator: '<S118>/Integrator' */
    ModelWithControllersOnly_B.v_g =
      ModelWithControllersOnly_X.Integrator_CSTATE_c;

    /* Math: '<S118>/Square' */
    ModelWithControllersOnly_B.Square_g = ModelWithControllersOnly_B.v_g *
      ModelWithControllersOnly_B.v_g;

    /* Gain: '<S118>/Gain1' */
    ModelWithControllersOnly_B.F_aero_i =
      ModelWithControllersOnly_cal->Gain1_Gain_p *
      ModelWithControllersOnly_B.Square_g;
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Gain: '<S118>/Gain2' incorporates:
       *  Constant: '<S118>/Constant1'
       */
      ModelWithControllersOnly_B.Gain2_d =
        ModelWithControllersOnly_cal->Gain2_Gain_h * b_m;

      /* RateTransition generated from: '<S118>/Sum' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.TmpRTBAtSumInport3_RdBufIdx_l = static_cast<
          int8_T>(ModelWithControllersOnly_DW.TmpRTBAtSumInport3_RdBufIdx_l == 0);
      }

      /* RateTransition generated from: '<S118>/Sum' */
      ModelWithControllersOnly_B.TmpRTBAtSumInport3_l =
        ModelWithControllersOnly_DW.TmpRTBAtSumInport3_Buf_f[ModelWithControllersOnly_DW.TmpRTBAtSumInport3_RdBufIdx_l];
    }

    /* Sum: '<S118>/Sum' */
    ModelWithControllersOnly_B.Sum_n = (ModelWithControllersOnly_B.F_aero_i +
      ModelWithControllersOnly_B.Gain2_d) +
      ModelWithControllersOnly_B.TmpRTBAtSumInport3_l;

    /* Gain: '<S118>/Gain4' */
    posError = tmp_5 / 5.6400000000000006;

    /* Gain: '<S118>/Gain4' */
    ModelWithControllersOnly_B.T_eng_c = posError *
      ModelWithControllersOnly_B.Sum_n;

    /* RateTransition generated from: '<S118>/Integrator' */
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_RdBuf =
          static_cast<int8_T>
          (ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_RdBuf == 0);
      }

      /* RateTransition generated from: '<S118>/Integrator' */
      ModelWithControllersOnly_B.TmpRTBAtIntegratorInport1 =
        ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_Buf[ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_RdBuf];

      /* RateTransition generated from: '<S119>/Integrator' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_RdB_k =
          static_cast<int8_T>
          (ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_RdB_k == 0);
      }

      /* RateTransition generated from: '<S119>/Integrator' */
      ModelWithControllersOnly_B.TmpRTBAtIntegratorInport1_c =
        ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_Buf_p[ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_RdB_k];
    }

    /* End of RateTransition generated from: '<S118>/Integrator' */

    /* Saturate: '<S4>/Throttle Saturation' */
    u = ModelWithControllersOnly_B.T_eng_c;
    posError = tmp_b;
    u2 = tmp_c;
    if (u > u2) {
      /* Saturate: '<S4>/Throttle Saturation' */
      ModelWithControllersOnly_B.ThrottleSaturation = u2;
    } else if (u < posError) {
      /* Saturate: '<S4>/Throttle Saturation' */
      ModelWithControllersOnly_B.ThrottleSaturation = posError;
    } else {
      /* Saturate: '<S4>/Throttle Saturation' */
      ModelWithControllersOnly_B.ThrottleSaturation = u;
    }

    /* End of Saturate: '<S4>/Throttle Saturation' */

    /* Switch: '<S4>/Switch' */
    if (ModelWithControllersOnly_B.ThrottleSaturation >
        ModelWithControllersOnly_cal->Switch_Threshold) {
      /* Switch: '<S4>/Switch' */
      ModelWithControllersOnly_B.Switch =
        ModelWithControllersOnly_B.ThrottleSaturation;
    } else {
      /* Switch: '<S4>/Switch' incorporates:
       *  Constant: '<S4>/Constant'
       */
      ModelWithControllersOnly_B.Switch =
        ModelWithControllersOnly_cal->Constant_Value_d;
    }

    /* End of Switch: '<S4>/Switch' */

    /* Switch: '<S4>/Switch1' */
    if (ModelWithControllersOnly_B.BrakeSaturation >
        ModelWithControllersOnly_cal->Switch1_Threshold) {
      /* Switch: '<S4>/Switch1' incorporates:
       *  Constant: '<S4>/Constant'
       */
      ModelWithControllersOnly_B.Switch1 =
        ModelWithControllersOnly_cal->Constant_Value_d;
    } else {
      /* Switch: '<S4>/Switch1' */
      ModelWithControllersOnly_B.Switch1 =
        ModelWithControllersOnly_B.BrakeSaturation;
    }

    /* End of Switch: '<S4>/Switch1' */
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Constant: '<S5>/Constant9' */
      ModelWithControllersOnly_B.gear_cmd =
        ModelWithControllersOnly_cal->Constant9_Value;

      /* Gain: '<S120>/Gain' incorporates:
       *  Constant: '<S149>/Constant2'
       */
      ModelWithControllersOnly_B.Roll =
        ModelWithControllersOnly_cal->Gain_Gain_l *
        ModelWithControllersOnly_cal->Constant2_Value;
    }

    /* MATLAB Function: '<S149>/COMB2I' */
    ModelWithControllersOnly_B.y[0] = ModelWithControllersOnly_B.Integrator[0] *
      std::cos(ModelWithControllersOnly_B.Integrator[2]) -
      ModelWithControllersOnly_B.Integrator[1] * std::sin
      (ModelWithControllersOnly_B.Integrator[2]);
    ModelWithControllersOnly_B.y[1] = ModelWithControllersOnly_B.Integrator[0] *
      std::sin(ModelWithControllersOnly_B.Integrator[2]) +
      ModelWithControllersOnly_B.Integrator[1] * std::cos
      (ModelWithControllersOnly_B.Integrator[2]);

    /* Derivative: '<S120>/Derivative' */
    if ((ModelWithControllersOnly_DW.TimeStampA >=
         ModelWithControllersOnly_M->Timing.t[0]) &&
        (ModelWithControllersOnly_DW.TimeStampB >=
         ModelWithControllersOnly_M->Timing.t[0])) {
      /* Derivative: '<S120>/Derivative' */
      ModelWithControllersOnly_B.Derivative = 0.0;
    } else {
      u = ModelWithControllersOnly_DW.TimeStampA;
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA;
      if (ModelWithControllersOnly_DW.TimeStampA <
          ModelWithControllersOnly_DW.TimeStampB) {
        if (ModelWithControllersOnly_DW.TimeStampB <
            ModelWithControllersOnly_M->Timing.t[0]) {
          u = ModelWithControllersOnly_DW.TimeStampB;
          lastU = &ModelWithControllersOnly_DW.LastUAtTimeB;
        }
      } else if (ModelWithControllersOnly_DW.TimeStampA >=
                 ModelWithControllersOnly_M->Timing.t[0]) {
        u = ModelWithControllersOnly_DW.TimeStampB;
        lastU = &ModelWithControllersOnly_DW.LastUAtTimeB;
      }

      u = ModelWithControllersOnly_M->Timing.t[0] - u;

      /* Derivative: '<S120>/Derivative' */
      ModelWithControllersOnly_B.Derivative = (ModelWithControllersOnly_B.y[0] -
        *lastU) / u;
    }

    /* End of Derivative: '<S120>/Derivative' */

    /* Gain: '<S120>/Gain5' */
    ModelWithControllersOnly_B.Gain5 = ModelWithControllersOnly_cal->Gain5_Gain *
      ModelWithControllersOnly_B.y[1];

    /* Derivative: '<S120>/Derivative1' */
    if ((ModelWithControllersOnly_DW.TimeStampA_l >=
         ModelWithControllersOnly_M->Timing.t[0]) &&
        (ModelWithControllersOnly_DW.TimeStampB_l >=
         ModelWithControllersOnly_M->Timing.t[0])) {
      /* Derivative: '<S120>/Derivative1' */
      ModelWithControllersOnly_B.Derivative1 = 0.0;
    } else {
      u = ModelWithControllersOnly_DW.TimeStampA_l;
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA_d;
      if (ModelWithControllersOnly_DW.TimeStampA_l <
          ModelWithControllersOnly_DW.TimeStampB_l) {
        if (ModelWithControllersOnly_DW.TimeStampB_l <
            ModelWithControllersOnly_M->Timing.t[0]) {
          u = ModelWithControllersOnly_DW.TimeStampB_l;
          lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_o;
        }
      } else if (ModelWithControllersOnly_DW.TimeStampA_l >=
                 ModelWithControllersOnly_M->Timing.t[0]) {
        u = ModelWithControllersOnly_DW.TimeStampB_l;
        lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_o;
      }

      u = ModelWithControllersOnly_M->Timing.t[0] - u;

      /* Derivative: '<S120>/Derivative1' */
      ModelWithControllersOnly_B.Derivative1 = (ModelWithControllersOnly_B.Gain5
        - *lastU) / u;
    }

    /* End of Derivative: '<S120>/Derivative1' */
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Gain: '<S120>/Gain6' incorporates:
       *  Constant: '<S149>/Constant8'
       */
      ModelWithControllersOnly_B.Gain6 =
        ModelWithControllersOnly_cal->Gain6_Gain *
        ModelWithControllersOnly_cal->Constant8_Value;
    }

    /* Derivative: '<S120>/Derivative2' */
    if ((ModelWithControllersOnly_DW.TimeStampA_e >=
         ModelWithControllersOnly_M->Timing.t[0]) &&
        (ModelWithControllersOnly_DW.TimeStampB_k >=
         ModelWithControllersOnly_M->Timing.t[0])) {
      /* Derivative: '<S120>/Derivative2' */
      ModelWithControllersOnly_B.Derivative2 = 0.0;
    } else {
      u = ModelWithControllersOnly_DW.TimeStampA_e;
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA_n;
      if (ModelWithControllersOnly_DW.TimeStampA_e <
          ModelWithControllersOnly_DW.TimeStampB_k) {
        if (ModelWithControllersOnly_DW.TimeStampB_k <
            ModelWithControllersOnly_M->Timing.t[0]) {
          u = ModelWithControllersOnly_DW.TimeStampB_k;
          lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_d;
        }
      } else if (ModelWithControllersOnly_DW.TimeStampA_e >=
                 ModelWithControllersOnly_M->Timing.t[0]) {
        u = ModelWithControllersOnly_DW.TimeStampB_k;
        lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_d;
      }

      u = ModelWithControllersOnly_M->Timing.t[0] - u;

      /* Derivative: '<S120>/Derivative2' */
      ModelWithControllersOnly_B.Derivative2 = (ModelWithControllersOnly_B.Gain6
        - *lastU) / u;
    }

    /* End of Derivative: '<S120>/Derivative2' */

    /* Derivative: '<S120>/Derivative3' */
    if ((ModelWithControllersOnly_DW.TimeStampA_o >=
         ModelWithControllersOnly_M->Timing.t[0]) &&
        (ModelWithControllersOnly_DW.TimeStampB_lf >=
         ModelWithControllersOnly_M->Timing.t[0])) {
      /* Derivative: '<S120>/Derivative3' */
      ModelWithControllersOnly_B.Derivative3 = 0.0;
    } else {
      u = ModelWithControllersOnly_DW.TimeStampA_o;
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA_nv;
      if (ModelWithControllersOnly_DW.TimeStampA_o <
          ModelWithControllersOnly_DW.TimeStampB_lf) {
        if (ModelWithControllersOnly_DW.TimeStampB_lf <
            ModelWithControllersOnly_M->Timing.t[0]) {
          u = ModelWithControllersOnly_DW.TimeStampB_lf;
          lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_g;
        }
      } else if (ModelWithControllersOnly_DW.TimeStampA_o >=
                 ModelWithControllersOnly_M->Timing.t[0]) {
        u = ModelWithControllersOnly_DW.TimeStampB_lf;
        lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_g;
      }

      u = ModelWithControllersOnly_M->Timing.t[0] - u;

      /* Derivative: '<S120>/Derivative3' */
      ModelWithControllersOnly_B.Derivative3 = (ModelWithControllersOnly_B.Roll
        - *lastU) / u;
    }

    /* End of Derivative: '<S120>/Derivative3' */
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Gain: '<S120>/Gain1' incorporates:
       *  Constant: '<S149>/Constant7'
       */
      ModelWithControllersOnly_B.Pitch =
        ModelWithControllersOnly_cal->Gain1_Gain_b *
        ModelWithControllersOnly_cal->Constant7_Value;
    }

    /* Derivative: '<S120>/Derivative4' */
    if ((ModelWithControllersOnly_DW.TimeStampA_l5 >=
         ModelWithControllersOnly_M->Timing.t[0]) &&
        (ModelWithControllersOnly_DW.TimeStampB_e >=
         ModelWithControllersOnly_M->Timing.t[0])) {
      /* Derivative: '<S120>/Derivative4' */
      ModelWithControllersOnly_B.Derivative4 = 0.0;
    } else {
      u = ModelWithControllersOnly_DW.TimeStampA_l5;
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA_k;
      if (ModelWithControllersOnly_DW.TimeStampA_l5 <
          ModelWithControllersOnly_DW.TimeStampB_e) {
        if (ModelWithControllersOnly_DW.TimeStampB_e <
            ModelWithControllersOnly_M->Timing.t[0]) {
          u = ModelWithControllersOnly_DW.TimeStampB_e;
          lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_p;
        }
      } else if (ModelWithControllersOnly_DW.TimeStampA_l5 >=
                 ModelWithControllersOnly_M->Timing.t[0]) {
        u = ModelWithControllersOnly_DW.TimeStampB_e;
        lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_p;
      }

      u = ModelWithControllersOnly_M->Timing.t[0] - u;

      /* Derivative: '<S120>/Derivative4' */
      ModelWithControllersOnly_B.Derivative4 = (ModelWithControllersOnly_B.Pitch
        - *lastU) / u;
    }

    /* End of Derivative: '<S120>/Derivative4' */

    /* Gain: '<S120>/Gain2' */
    ModelWithControllersOnly_B.Yawdeg =
      ModelWithControllersOnly_cal->Gain2_Gain_c *
      ModelWithControllersOnly_B.Integrator[2];

    /* Derivative: '<S120>/Derivative5' */
    if ((ModelWithControllersOnly_DW.TimeStampA_k >=
         ModelWithControllersOnly_M->Timing.t[0]) &&
        (ModelWithControllersOnly_DW.TimeStampB_l4 >=
         ModelWithControllersOnly_M->Timing.t[0])) {
      /* Derivative: '<S120>/Derivative5' */
      ModelWithControllersOnly_B.Derivative5 = 0.0;
    } else {
      u = ModelWithControllersOnly_DW.TimeStampA_k;
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA_h;
      if (ModelWithControllersOnly_DW.TimeStampA_k <
          ModelWithControllersOnly_DW.TimeStampB_l4) {
        if (ModelWithControllersOnly_DW.TimeStampB_l4 <
            ModelWithControllersOnly_M->Timing.t[0]) {
          u = ModelWithControllersOnly_DW.TimeStampB_l4;
          lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_f;
        }
      } else if (ModelWithControllersOnly_DW.TimeStampA_k >=
                 ModelWithControllersOnly_M->Timing.t[0]) {
        u = ModelWithControllersOnly_DW.TimeStampB_l4;
        lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_f;
      }

      u = ModelWithControllersOnly_M->Timing.t[0] - u;

      /* Derivative: '<S120>/Derivative5' */
      ModelWithControllersOnly_B.Derivative5 =
        (ModelWithControllersOnly_B.Yawdeg - *lastU) / u;
    }

    /* End of Derivative: '<S120>/Derivative5' */

    /* Gain: '<S120>/Gain3' */
    ModelWithControllersOnly_B.Gain3 = ModelWithControllersOnly_cal->Gain3_Gain *
      ModelWithControllersOnly_B.Integrator_m[1];
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Gain: '<S120>/Gain4' incorporates:
       *  Constant: '<S149>/Constant'
       */
      ModelWithControllersOnly_B.Gain4 =
        ModelWithControllersOnly_cal->Gain4_Gain *
        ModelWithControllersOnly_cal->Constant_Value_a;

      /* RateTransition: '<S120>/Rate Transition' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.RateTransition_WrBufIdx = static_cast<int8_T>
          (ModelWithControllersOnly_DW.RateTransition_WrBufIdx == 0);
      }

      ModelWithControllersOnly_DW.RateTransition_Buf[ModelWithControllersOnly_DW.RateTransition_WrBufIdx]
        = ModelWithControllersOnly_B.Roll;

      /* End of RateTransition: '<S120>/Rate Transition' */

      /* RateTransition: '<S120>/Rate Transition1' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.RateTransition1_WrBufIdx =
          static_cast<int8_T>
          (ModelWithControllersOnly_DW.RateTransition1_WrBufIdx == 0);
      }

      ModelWithControllersOnly_DW.RateTransition1_Buf[ModelWithControllersOnly_DW.RateTransition1_WrBufIdx]
        = ModelWithControllersOnly_B.Pitch;

      /* End of RateTransition: '<S120>/Rate Transition1' */

      /* RateTransition: '<S120>/Rate Transition2' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.RateTransition2_WrBufIdx =
          static_cast<int8_T>
          (ModelWithControllersOnly_DW.RateTransition2_WrBufIdx == 0);
      }

      ModelWithControllersOnly_DW.RateTransition2_Buf[ModelWithControllersOnly_DW.RateTransition2_WrBufIdx]
        = ModelWithControllersOnly_B.Yawdeg;

      /* End of RateTransition: '<S120>/Rate Transition2' */

      /* RateTransition generated from: '<S120>/MATLAB Function1' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_j =
          static_cast<int8_T>
          (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_j == 0);
      }

      ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport1_[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_j]
        = ModelWithControllersOnly_B.y[0];

      /* End of RateTransition generated from: '<S120>/MATLAB Function1' */

      /* RateTransition generated from: '<S120>/MATLAB Function1' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_m =
          static_cast<int8_T>
          (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_m == 0);
      }

      ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport2_[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_m]
        = ModelWithControllersOnly_B.Gain5;

      /* End of RateTransition generated from: '<S120>/MATLAB Function1' */

      /* RateTransition generated from: '<S120>/MATLAB Function1' */
      rtw_slrealtime_mutex_lock
        (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_k);
      wrBufIdx = static_cast<int8_T>
        (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inpor_oh + 1);
      if (wrBufIdx == 3) {
        wrBufIdx = 0;
      }

      if (wrBufIdx ==
          ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_b) {
        wrBufIdx = static_cast<int8_T>(wrBufIdx + 1);
        if (wrBufIdx == 3) {
          wrBufIdx = 0;
        }
      }

      rtw_slrealtime_mutex_unlock
        (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_k);
      switch (wrBufIdx) {
       case 0:
        ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport3_ =
          ModelWithControllersOnly_B.Gain6;
        break;

       case 1:
        ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_f =
          ModelWithControllersOnly_B.Gain6;
        break;

       case 2:
        ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_o =
          ModelWithControllersOnly_B.Gain6;
        break;
      }

      ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inpor_oh = wrBufIdx;

      /* End of RateTransition generated from: '<S120>/MATLAB Function1' */

      /* RateTransition generated from: '<S120>/MATLAB Function2' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_d =
          static_cast<int8_T>
          (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_d == 0);
      }

      ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport1_[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_d]
        = ModelWithControllersOnly_B.Derivative;

      /* End of RateTransition generated from: '<S120>/MATLAB Function2' */

      /* RateTransition generated from: '<S120>/MATLAB Function2' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_k =
          static_cast<int8_T>
          (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_k == 0);
      }

      ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport2_[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_k]
        = ModelWithControllersOnly_B.Derivative1;

      /* End of RateTransition generated from: '<S120>/MATLAB Function2' */

      /* RateTransition generated from: '<S120>/MATLAB Function2' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_p =
          static_cast<int8_T>
          (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_p == 0);
      }

      ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport3_[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_p]
        = ModelWithControllersOnly_B.Derivative2;

      /* End of RateTransition generated from: '<S120>/MATLAB Function2' */

      /* RateTransition generated from: '<S120>/MATLAB Function3' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport_g =
          static_cast<int8_T>
          (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport_g == 0);
      }

      ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport1_[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport_g]
        = ModelWithControllersOnly_B.Derivative3;

      /* End of RateTransition generated from: '<S120>/MATLAB Function3' */

      /* RateTransition generated from: '<S120>/MATLAB Function3' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport_b =
          static_cast<int8_T>
          (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport_b == 0);
      }

      ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport2_[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport_b]
        = ModelWithControllersOnly_B.Derivative4;

      /* End of RateTransition generated from: '<S120>/MATLAB Function3' */

      /* RateTransition generated from: '<S120>/MATLAB Function3' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inpor_bx =
          static_cast<int8_T>
          (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inpor_bx == 0);
      }

      ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport3_[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inpor_bx]
        = ModelWithControllersOnly_B.Derivative5;

      /* End of RateTransition generated from: '<S120>/MATLAB Function3' */

      /* RateTransition generated from: '<S120>/MATLAB Function' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport1_W =
          static_cast<int8_T>
          (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport1_W == 0);
      }

      ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport1_B[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport1_W]
        = ModelWithControllersOnly_B.Integrator_m[0];

      /* End of RateTransition generated from: '<S120>/MATLAB Function' */

      /* RateTransition generated from: '<S120>/MATLAB Function' */
      if (ModelWithControllersOnly_M->Timing.RateInteraction.TID1_2 == 1) {
        ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport2_W =
          static_cast<int8_T>
          (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport2_W == 0);
      }

      ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport2_B[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport2_W]
        = ModelWithControllersOnly_B.Gain3;

      /* End of RateTransition generated from: '<S120>/MATLAB Function' */

      /* RateTransition generated from: '<S120>/MATLAB Function' */
      rtw_slrealtime_mutex_lock
        (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_d);
      wrBufIdx = static_cast<int8_T>
        (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_L + 1);
      if (wrBufIdx == 3) {
        wrBufIdx = 0;
      }

      if (wrBufIdx ==
          ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_R) {
        wrBufIdx = static_cast<int8_T>(wrBufIdx + 1);
        if (wrBufIdx == 3) {
          wrBufIdx = 0;
        }
      }

      rtw_slrealtime_mutex_unlock
        (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_d);
      switch (wrBufIdx) {
       case 0:
        ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_B =
          ModelWithControllersOnly_B.Gain4;
        break;

       case 1:
        ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_c =
          ModelWithControllersOnly_B.Gain4;
        break;

       case 2:
        ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_f =
          ModelWithControllersOnly_B.Gain4;
        break;
      }

      ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_L = wrBufIdx;

      /* End of RateTransition generated from: '<S120>/MATLAB Function' */
    }

    /* Derivative: '<S121>/Derivative' */
    if ((ModelWithControllersOnly_DW.TimeStampA_h >=
         ModelWithControllersOnly_M->Timing.t[0]) &&
        (ModelWithControllersOnly_DW.TimeStampB_h >=
         ModelWithControllersOnly_M->Timing.t[0])) {
      /* Derivative: '<S121>/Derivative' */
      ModelWithControllersOnly_B.EgoVehLongAcc = 0.0;
    } else {
      u = ModelWithControllersOnly_DW.TimeStampA_h;
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA_o;
      if (ModelWithControllersOnly_DW.TimeStampA_h <
          ModelWithControllersOnly_DW.TimeStampB_h) {
        if (ModelWithControllersOnly_DW.TimeStampB_h <
            ModelWithControllersOnly_M->Timing.t[0]) {
          u = ModelWithControllersOnly_DW.TimeStampB_h;
          lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_c;
        }
      } else if (ModelWithControllersOnly_DW.TimeStampA_h >=
                 ModelWithControllersOnly_M->Timing.t[0]) {
        u = ModelWithControllersOnly_DW.TimeStampB_h;
        lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_c;
      }

      u = ModelWithControllersOnly_M->Timing.t[0] - u;

      /* Derivative: '<S121>/Derivative' */
      ModelWithControllersOnly_B.EgoVehLongAcc =
        (ModelWithControllersOnly_B.Integrator[0] - *lastU) / u;
    }

    /* End of Derivative: '<S121>/Derivative' */

    /* FirstOrderHold: '<S121>/First Order Hold' */
    ModelWithControllersOnly_B.FirstOrderHold = ModelWithControllersOnly_DW.Ck;

    /* FirstOrderHold: '<S121>/First Order Hold' */
    if (ModelWithControllersOnly_DW.Tk != (rtInf)) {
      u = ModelWithControllersOnly_M->Timing.t[0] -
        ModelWithControllersOnly_DW.Tk;

      /* FirstOrderHold: '<S121>/First Order Hold' */
      ModelWithControllersOnly_B.FirstOrderHold +=
        ModelWithControllersOnly_DW.Mk * u;
    }

    /* FirstOrderHold: '<S121>/First Order Hold1' */
    ModelWithControllersOnly_B.TnetNm = ModelWithControllersOnly_DW.Ck_o;

    /* FirstOrderHold: '<S121>/First Order Hold1' */
    if (ModelWithControllersOnly_DW.Tk_d != (rtInf)) {
      u = ModelWithControllersOnly_M->Timing.t[0] -
        ModelWithControllersOnly_DW.Tk_d;

      /* FirstOrderHold: '<S121>/First Order Hold1' */
      ModelWithControllersOnly_B.TnetNm += ModelWithControllersOnly_DW.Mk_n * u;
    }

    /* Gain: '<S121>/TireRadius' */
    ModelWithControllersOnly_B.wheelspeed = tmp_5 *
      ModelWithControllersOnly_B.Integrator[0];

    /* SignalConversion generated from: '<S129>/Vector Concatenate1' incorporates:
     *  Concatenate: '<S129>/Vector Concatenate1'
     */
    ModelWithControllersOnly_B.VectorConcatenate1[0] =
      ModelWithControllersOnly_B.Integrator[0];

    /* SignalConversion generated from: '<S129>/Vector Concatenate1' incorporates:
     *  Concatenate: '<S129>/Vector Concatenate1'
     */
    ModelWithControllersOnly_B.VectorConcatenate1[1] =
      ModelWithControllersOnly_B.Integrator[1];
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Constant: '<S129>/Constant' incorporates:
       *  Concatenate: '<S129>/Vector Concatenate1'
       */
      ModelWithControllersOnly_B.VectorConcatenate1[2] =
        ModelWithControllersOnly_cal->Constant_Value_k;

      /* SignalConversion generated from: '<S215>/Vector Concatenate1' incorporates:
       *  Concatenate: '<S215>/Vector Concatenate1'
       */
      ModelWithControllersOnly_B.VectorConcatenate1_b[0] = 0.0;

      /* SignalConversion generated from: '<S215>/Vector Concatenate1' incorporates:
       *  Concatenate: '<S215>/Vector Concatenate1'
       */
      ModelWithControllersOnly_B.VectorConcatenate1_b[1] = 0.0;

      /* SignalConversion generated from: '<S215>/Vector Concatenate1' incorporates:
       *  Concatenate: '<S215>/Vector Concatenate1'
       */
      ModelWithControllersOnly_B.VectorConcatenate1_b[2] = 0.0;
    }

    /* Trigonometry: '<S143>/Trigonometric Function' */
    posError = ModelWithControllersOnly_B.Integrator[2];
    u = std::sin(posError);
    posError = std::cos(posError);

    /* Trigonometry: '<S143>/Trigonometric Function' */
    ModelWithControllersOnly_B.TrigonometricFunction_o1 = u;

    /* Trigonometry: '<S143>/Trigonometric Function' */
    ModelWithControllersOnly_B.TrigonometricFunction_o2 = posError;

    /* Product: '<S143>/Product3' */
    ModelWithControllersOnly_B.Product3 =
      ModelWithControllersOnly_B.TrigonometricFunction_o2 *
      ModelWithControllersOnly_B.VectorConcatenate1_b[0];

    /* Product: '<S143>/Product1' */
    ModelWithControllersOnly_B.Product1 =
      ModelWithControllersOnly_B.TrigonometricFunction_o1 *
      ModelWithControllersOnly_B.VectorConcatenate1_b[1];

    /* Sum: '<S143>/Add' incorporates:
     *  Concatenate: '<S143>/Vector Concatenate2'
     */
    ModelWithControllersOnly_B.VectorConcatenate2[0] =
      ModelWithControllersOnly_B.Product3 + ModelWithControllersOnly_B.Product1;

    /* Product: '<S143>/Product' */
    ModelWithControllersOnly_B.Product_k =
      ModelWithControllersOnly_B.TrigonometricFunction_o1 *
      ModelWithControllersOnly_B.VectorConcatenate1_b[0];

    /* Product: '<S143>/Product2' */
    ModelWithControllersOnly_B.Product2 =
      ModelWithControllersOnly_B.TrigonometricFunction_o2 *
      ModelWithControllersOnly_B.VectorConcatenate1_b[1];

    /* Sum: '<S143>/Add1' incorporates:
     *  Concatenate: '<S143>/Vector Concatenate2'
     */
    ModelWithControllersOnly_B.VectorConcatenate2[1] =
      ModelWithControllersOnly_B.Product2 - ModelWithControllersOnly_B.Product_k;

    /* SignalConversion generated from: '<S143>/Vector Concatenate2' incorporates:
     *  Concatenate: '<S143>/Vector Concatenate2'
     */
    ModelWithControllersOnly_B.VectorConcatenate2[2] =
      ModelWithControllersOnly_B.VectorConcatenate1_b[2];

    /* Sum: '<S142>/Add1' */
    ModelWithControllersOnly_B.Add1[0] =
      ModelWithControllersOnly_B.VectorConcatenate1[0] -
      ModelWithControllersOnly_B.VectorConcatenate2[0];

    /* Product: '<S142>/Product' */
    ModelWithControllersOnly_B.Product_i[0] = ModelWithControllersOnly_B.Add1[0]
      * ModelWithControllersOnly_B.Add1[0];

    /* Sum: '<S142>/Add1' */
    ModelWithControllersOnly_B.Add1[1] =
      ModelWithControllersOnly_B.VectorConcatenate1[1] -
      ModelWithControllersOnly_B.VectorConcatenate2[1];

    /* Product: '<S142>/Product' */
    ModelWithControllersOnly_B.Product_i[1] = ModelWithControllersOnly_B.Add1[1]
      * ModelWithControllersOnly_B.Add1[1];

    /* Sum: '<S142>/Add1' */
    ModelWithControllersOnly_B.Add1[2] =
      ModelWithControllersOnly_B.VectorConcatenate1[2] -
      ModelWithControllersOnly_B.VectorConcatenate2[2];

    /* Product: '<S142>/Product' */
    ModelWithControllersOnly_B.Product_i[2] = ModelWithControllersOnly_B.Add1[2]
      * ModelWithControllersOnly_B.Add1[2];

    /* Sum: '<S142>/Sum of Elements' */
    posError = ModelWithControllersOnly_B.Product_i[0];
    posError += ModelWithControllersOnly_B.Product_i[1];
    posError += ModelWithControllersOnly_B.Product_i[2];

    /* Sum: '<S142>/Sum of Elements' */
    ModelWithControllersOnly_B.SumofElements = posError;

    /* Sqrt: '<S142>/Sqrt' */
    ModelWithControllersOnly_B.Sqrt = std::sqrt
      (ModelWithControllersOnly_B.SumofElements);

    /* Product: '<S142>/Product2' */
    ModelWithControllersOnly_B.Product2_h = ModelWithControllersOnly_B.Sqrt *
      ModelWithControllersOnly_B.Sqrt;
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Constant: '<S142>/Constant' incorporates:
       *  Concatenate: '<S142>/Vector Concatenate'
       */
      ModelWithControllersOnly_B.VectorConcatenate_b[0] =
        ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_Cd;
    }

    /* Trigonometry: '<S142>/Trigonometric Function' */
    ModelWithControllersOnly_B.TrigonometricFunction = rt_atan2d_snf
      (ModelWithControllersOnly_B.Add1[1], ModelWithControllersOnly_B.Add1[0]);

    /* Lookup_n-D: '<S142>/Cs' incorporates:
     *  Concatenate: '<S142>/Vector Concatenate'
     *  Trigonometry: '<S142>/Trigonometric Function'
     */
    ModelWithControllersOnly_B.VectorConcatenate_b[1] = look1_binlcpw
      (ModelWithControllersOnly_B.TrigonometricFunction,
       ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_beta,
       ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_Cs, 30U);
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Constant: '<S142>/Constant1' incorporates:
       *  Concatenate: '<S142>/Vector Concatenate'
       */
      ModelWithControllersOnly_B.VectorConcatenate_b[2] =
        ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_Cl;
    }

    /* Lookup_n-D: '<S142>/Crm' incorporates:
     *  Concatenate: '<S142>/Vector Concatenate'
     *  Trigonometry: '<S142>/Trigonometric Function'
     */
    ModelWithControllersOnly_B.VectorConcatenate_b[3] = look1_binlxpw
      (ModelWithControllersOnly_B.TrigonometricFunction,
       ModelWithControllersOnly_cal->Crm_bp01Data,
       ModelWithControllersOnly_cal->Crm_tableData, 1U);

    /* Gain: '<S142>/4' */
    ModelWithControllersOnly_B.u[0] = ModelWithControllersOnly_cal->u_Gain[0] *
      ModelWithControllersOnly_B.Add1[0];

    /* Trigonometry: '<S142>/Tanh' */
    ModelWithControllersOnly_B.Tanh[0] = std::tanh(ModelWithControllersOnly_B.u
      [0]);

    /* Gain: '<S142>/4' */
    ModelWithControllersOnly_B.u[1] = ModelWithControllersOnly_cal->u_Gain[1] *
      ModelWithControllersOnly_B.Add1[1];

    /* Trigonometry: '<S142>/Tanh' */
    ModelWithControllersOnly_B.Tanh[1] = std::tanh(ModelWithControllersOnly_B.u
      [1]);

    /* Gain: '<S142>/4' */
    ModelWithControllersOnly_B.u[2] = ModelWithControllersOnly_cal->u_Gain[2] *
      ModelWithControllersOnly_B.Add1[2];

    /* Trigonometry: '<S142>/Tanh' */
    ModelWithControllersOnly_B.Tanh[2] = std::tanh(ModelWithControllersOnly_B.u
      [2]);

    /* Product: '<S142>/Product5' incorporates:
     *  Concatenate: '<S142>/Vector Concatenate'
     *  Constant: '<S142>/Constant2'
     */
    ModelWithControllersOnly_B.VectorConcatenate_b[4] =
      ModelWithControllersOnly_B.Tanh[0] *
      ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_Cpm;

    /* Lookup_n-D: '<S142>/Cym' incorporates:
     *  Concatenate: '<S142>/Vector Concatenate'
     *  Trigonometry: '<S142>/Trigonometric Function'
     */
    ModelWithControllersOnly_B.VectorConcatenate_b[5] = look1_binlxpw
      (ModelWithControllersOnly_B.TrigonometricFunction,
       ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_beta,
       ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_Cym, 30U);

    /* Gain: '<S142>/.5.*A.*Pabs.//R.//T' */
    posError = 0.5 * ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_Af
      * ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_Pabs /
      ModelWithControllersOnly_cal->DragForce_R;
    for (i = 0; i < 6; i++) {
      /* Product: '<S142>/Product1' incorporates:
       *  Constant: '<S126>/AirTempConstant'
       */
      ModelWithControllersOnly_B.Product1_a[i] =
        ModelWithControllersOnly_B.Product2_h *
        ModelWithControllersOnly_B.VectorConcatenate_b[i] /
        ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_Tair;

      /* Gain: '<S142>/.5.*A.*Pabs.//R.//T' */
      ModelWithControllersOnly_B.uAPabsRT[i] = posError *
        ModelWithControllersOnly_B.Product1_a[i];
    }

    /* Product: '<S142>/Product4' incorporates:
     *  Constant: '<S142>/Constant3'
     */
    posError = tmp_8 + b_b;

    /* Product: '<S142>/Product4' */
    ModelWithControllersOnly_B.Product4[0] =
      ModelWithControllersOnly_B.uAPabsRT[3] * posError;

    /* UnaryMinus: '<S129>/Unary Minus1' */
    ModelWithControllersOnly_B.UnaryMinus1[0] =
      -ModelWithControllersOnly_B.Product4[0];

    /* Sum: '<S126>/Add' */
    ModelWithControllersOnly_B.Add_d[0] =
      ModelWithControllersOnly_B.UnaryMinus1[0];

    /* Sum: '<S142>/Add2' incorporates:
     *  Constant: '<S142>/Constant4'
     */
    ModelWithControllersOnly_B.Add2[0] = ModelWithControllersOnly_B.Tanh[0] -
      ModelWithControllersOnly_cal->Constant4_Value_a[0];

    /* Product: '<S142>/Product3' */
    ModelWithControllersOnly_B.Product3_o[0] = ModelWithControllersOnly_B.Add2[0]
      * ModelWithControllersOnly_B.uAPabsRT[0];

    /* UnaryMinus: '<S129>/Unary Minus' */
    ModelWithControllersOnly_B.UnaryMinus[0] =
      -ModelWithControllersOnly_B.Product3_o[0];

    /* Sum: '<S126>/Add1' */
    ModelWithControllersOnly_B.Add1_d[0] =
      ModelWithControllersOnly_B.UnaryMinus[0];

    /* Product: '<S142>/Product4' */
    ModelWithControllersOnly_B.Product4[1] =
      ModelWithControllersOnly_B.uAPabsRT[4] * posError;

    /* UnaryMinus: '<S129>/Unary Minus1' */
    ModelWithControllersOnly_B.UnaryMinus1[1] =
      -ModelWithControllersOnly_B.Product4[1];

    /* Sum: '<S126>/Add' */
    ModelWithControllersOnly_B.Add_d[1] =
      ModelWithControllersOnly_B.UnaryMinus1[1];

    /* Sum: '<S142>/Add2' incorporates:
     *  Constant: '<S142>/Constant4'
     */
    ModelWithControllersOnly_B.Add2[1] = ModelWithControllersOnly_B.Tanh[1] -
      ModelWithControllersOnly_cal->Constant4_Value_a[1];

    /* Product: '<S142>/Product3' */
    ModelWithControllersOnly_B.Product3_o[1] = ModelWithControllersOnly_B.Add2[1]
      * ModelWithControllersOnly_B.uAPabsRT[1];

    /* UnaryMinus: '<S129>/Unary Minus' */
    ModelWithControllersOnly_B.UnaryMinus[1] =
      -ModelWithControllersOnly_B.Product3_o[1];

    /* Sum: '<S126>/Add1' */
    ModelWithControllersOnly_B.Add1_d[1] =
      ModelWithControllersOnly_B.UnaryMinus[1];

    /* Product: '<S142>/Product4' */
    ModelWithControllersOnly_B.Product4[2] =
      ModelWithControllersOnly_B.uAPabsRT[5] * posError;

    /* UnaryMinus: '<S129>/Unary Minus1' */
    ModelWithControllersOnly_B.UnaryMinus1[2] =
      -ModelWithControllersOnly_B.Product4[2];

    /* Sum: '<S126>/Add' */
    ModelWithControllersOnly_B.Add_d[2] =
      ModelWithControllersOnly_B.UnaryMinus1[2];

    /* Sum: '<S142>/Add2' incorporates:
     *  Constant: '<S142>/Constant4'
     */
    ModelWithControllersOnly_B.Add2[2] = ModelWithControllersOnly_B.Tanh[2] -
      ModelWithControllersOnly_cal->Constant4_Value_a[2];

    /* Product: '<S142>/Product3' */
    ModelWithControllersOnly_B.Product3_o[2] = ModelWithControllersOnly_B.Add2[2]
      * ModelWithControllersOnly_B.uAPabsRT[2];

    /* UnaryMinus: '<S129>/Unary Minus' */
    ModelWithControllersOnly_B.UnaryMinus[2] =
      -ModelWithControllersOnly_B.Product3_o[2];

    /* Sum: '<S126>/Add1' */
    ModelWithControllersOnly_B.Add1_d[2] =
      ModelWithControllersOnly_B.UnaryMinus[2];
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* SignalConversion generated from: '<S144>/Vector Concatenate1' incorporates:
       *  Concatenate: '<S144>/Vector Concatenate1'
       */
      ModelWithControllersOnly_B.VectorConcatenate1_i[0] = 0.0;

      /* SignalConversion generated from: '<S144>/Vector Concatenate1' incorporates:
       *  Concatenate: '<S144>/Vector Concatenate1'
       */
      ModelWithControllersOnly_B.VectorConcatenate1_i[1] = 0.0;
    }

    /* UnitConversion: '<S207>/Unit Conversion1' */
    /* Unit Conversion - from: rad to: rad
       Expression: output = (1*input) + (0) */
    ModelWithControllersOnly_B.UnitConversion1 =
      ModelWithControllersOnly_B.FirstOrderHold;

    /* SignalConversion generated from: '<S144>/Vector Concatenate' incorporates:
     *  Concatenate: '<S144>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_j[0] =
      ModelWithControllersOnly_B.UnitConversion1;

    /* SignalConversion generated from: '<S144>/Vector Concatenate' incorporates:
     *  Concatenate: '<S144>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_j[1] =
      ModelWithControllersOnly_B.UnitConversion1;

    /* SignalConversion generated from: '<S154>/sincos' incorporates:
     *  Constant: '<S149>/Constant2'
     *  Constant: '<S149>/Constant7'
     */
    ModelWithControllersOnly_B.TmpSignalConversionAtsincosInpo[0] =
      ModelWithControllersOnly_B.Integrator[2];
    ModelWithControllersOnly_B.TmpSignalConversionAtsincosInpo[1] =
      ModelWithControllersOnly_cal->Constant7_Value;
    ModelWithControllersOnly_B.TmpSignalConversionAtsincosInpo[2] =
      ModelWithControllersOnly_cal->Constant2_Value;

    /* Trigonometry: '<S154>/sincos' */
    posError = ModelWithControllersOnly_B.TmpSignalConversionAtsincosInpo[0];
    u = std::sin(posError);
    posError = std::cos(posError);
    ModelWithControllersOnly_B.sincos_o1[0] = u;
    ModelWithControllersOnly_B.sincos_o2[0] = posError;
    posError = ModelWithControllersOnly_B.TmpSignalConversionAtsincosInpo[1];
    u = std::sin(posError);
    posError = std::cos(posError);
    ModelWithControllersOnly_B.sincos_o1[1] = u;
    ModelWithControllersOnly_B.sincos_o2[1] = posError;
    posError = ModelWithControllersOnly_B.TmpSignalConversionAtsincosInpo[2];
    u = std::sin(posError);
    posError = std::cos(posError);
    ModelWithControllersOnly_B.sincos_o1[2] = u;
    ModelWithControllersOnly_B.sincos_o2[2] = posError;

    /* Fcn: '<S154>/Fcn11' incorporates:
     *  Concatenate: '<S158>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_f[0] =
      ModelWithControllersOnly_B.sincos_o2[0] *
      ModelWithControllersOnly_B.sincos_o2[1];

    /* Fcn: '<S154>/Fcn21' incorporates:
     *  Concatenate: '<S158>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_f[1] =
      ModelWithControllersOnly_B.sincos_o1[1] *
      ModelWithControllersOnly_B.sincos_o1[2] *
      ModelWithControllersOnly_B.sincos_o2[0] -
      ModelWithControllersOnly_B.sincos_o1[0] *
      ModelWithControllersOnly_B.sincos_o2[2];

    /* Fcn: '<S154>/Fcn31' incorporates:
     *  Concatenate: '<S158>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_f[2] =
      ModelWithControllersOnly_B.sincos_o1[1] *
      ModelWithControllersOnly_B.sincos_o2[2] *
      ModelWithControllersOnly_B.sincos_o2[0] +
      ModelWithControllersOnly_B.sincos_o1[0] *
      ModelWithControllersOnly_B.sincos_o1[2];

    /* Fcn: '<S154>/Fcn12' incorporates:
     *  Concatenate: '<S158>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_f[3] =
      ModelWithControllersOnly_B.sincos_o1[0] *
      ModelWithControllersOnly_B.sincos_o2[1];

    /* Fcn: '<S154>/Fcn22' incorporates:
     *  Concatenate: '<S158>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_f[4] =
      ModelWithControllersOnly_B.sincos_o1[1] *
      ModelWithControllersOnly_B.sincos_o1[2] *
      ModelWithControllersOnly_B.sincos_o1[0] +
      ModelWithControllersOnly_B.sincos_o2[0] *
      ModelWithControllersOnly_B.sincos_o2[2];

    /* Fcn: '<S154>/Fcn32' incorporates:
     *  Concatenate: '<S158>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_f[5] =
      ModelWithControllersOnly_B.sincos_o1[1] *
      ModelWithControllersOnly_B.sincos_o2[2] *
      ModelWithControllersOnly_B.sincos_o1[0] -
      ModelWithControllersOnly_B.sincos_o2[0] *
      ModelWithControllersOnly_B.sincos_o1[2];

    /* Fcn: '<S154>/Fcn13' incorporates:
     *  Concatenate: '<S158>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_f[6] =
      -ModelWithControllersOnly_B.sincos_o1[1];

    /* Fcn: '<S154>/Fcn23' incorporates:
     *  Concatenate: '<S158>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_f[7] =
      ModelWithControllersOnly_B.sincos_o2[1] *
      ModelWithControllersOnly_B.sincos_o1[2];

    /* Fcn: '<S154>/Fcn33' incorporates:
     *  Concatenate: '<S158>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_f[8] =
      ModelWithControllersOnly_B.sincos_o2[1] *
      ModelWithControllersOnly_B.sincos_o2[2];

    /* Reshape: '<S158>/Reshape (9) to [3x3] column-major' */
    std::memcpy(&ModelWithControllersOnly_B.Reshape9to3x3columnmajor[0],
                &ModelWithControllersOnly_B.VectorConcatenate_f[0], 9U * sizeof
                (real_T));
    for (i = 0; i < 3; i++) {
      /* Math: '<S150>/Transpose1' incorporates:
       *  Reshape: '<S158>/Reshape (9) to [3x3] column-major'
       */
      ModelWithControllersOnly_B.Transpose1[3 * i] =
        ModelWithControllersOnly_B.Reshape9to3x3columnmajor[i];
      ModelWithControllersOnly_B.Transpose1[3 * i + 1] =
        ModelWithControllersOnly_B.Reshape9to3x3columnmajor[i + 3];
      ModelWithControllersOnly_B.Transpose1[3 * i + 2] =
        ModelWithControllersOnly_B.Reshape9to3x3columnmajor[i + 6];
    }

    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Constant: '<S150>/R_T1' incorporates:
       *  Concatenate: '<S150>/Vector Concatenate'
       */
      ModelWithControllersOnly_B.VectorConcatenate_d[0] = tmp_8;

      /* Constant: '<S150>/R_T2' incorporates:
       *  Concatenate: '<S150>/Vector Concatenate'
       */
      ModelWithControllersOnly_B.VectorConcatenate_d[1] =
        ModelWithControllersOnly_cal->HardPointCoordinateTransformFro;

      /* Constant: '<S150>/R_T3' incorporates:
       *  Concatenate: '<S150>/Vector Concatenate'
       */
      ModelWithControllersOnly_B.VectorConcatenate_d[2] =
        ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_h;

      /* Reshape: '<S156>/Reshape1' */
      ModelWithControllersOnly_B.Reshape1[0] =
        ModelWithControllersOnly_B.VectorConcatenate_d[0];
      ModelWithControllersOnly_B.Reshape1[1] =
        ModelWithControllersOnly_B.VectorConcatenate_d[1];
      ModelWithControllersOnly_B.Reshape1[2] =
        ModelWithControllersOnly_B.VectorConcatenate_d[2];
    }

    /* Product: '<S156>/Product' incorporates:
     *  Math: '<S150>/Transpose1'
     *  Product: '<S155>/Product'
     *  Product: '<S164>/Product'
     *  Product: '<S165>/Product'
     *  Product: '<S174>/Product'
     *  Product: '<S175>/Product'
     *  Product: '<S182>/Product'
     *  Product: '<S183>/Product'
     *  Reshape: '<S156>/Reshape1'
     */
    lastU = &ModelWithControllersOnly_B.Transpose1[0];
    refPose[0] = ModelWithControllersOnly_B.Reshape1[0];
    refPose[1] = ModelWithControllersOnly_B.Reshape1[1];
    refPose[2] = ModelWithControllersOnly_B.Reshape1[2];
    for (i = 0; i < 3; i++) {
      /* Product: '<S156>/Product' */
      ModelWithControllersOnly_B.Product_f[i] = 0.0;
      ModelWithControllersOnly_B.Product_f[i] += lastU[i] * refPose[0];
      ModelWithControllersOnly_B.Product_f[i] += lastU[i + 3] * refPose[1];
      ModelWithControllersOnly_B.Product_f[i] += lastU[i + 6] * refPose[2];

      /* Reshape: '<S156>/Reshape2' incorporates:
       *  Product: '<S156>/Product'
       */
      ModelWithControllersOnly_B.Reshape2[i] =
        ModelWithControllersOnly_B.Product_f[i];
    }

    /* Sum: '<S150>/Add' incorporates:
     *  Constant: '<S149>/Constant'
     */
    ModelWithControllersOnly_B.Add_e[0] =
      ModelWithControllersOnly_B.Integrator_m[0] +
      ModelWithControllersOnly_B.Reshape2[0];
    ModelWithControllersOnly_B.Add_e[1] =
      ModelWithControllersOnly_B.Integrator_m[1] +
      ModelWithControllersOnly_B.Reshape2[1];
    ModelWithControllersOnly_B.Add_e[2] =
      ModelWithControllersOnly_cal->Constant_Value_a +
      ModelWithControllersOnly_B.Reshape2[2];
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Product: '<S159>/j x k' incorporates:
       *  Constant: '<S149>/Constant5'
       */
      ModelWithControllersOnly_B.jxk =
        ModelWithControllersOnly_cal->Constant5_Value *
        ModelWithControllersOnly_B.VectorConcatenate_d[2];
    }

    /* Product: '<S159>/k x i' */
    ModelWithControllersOnly_B.kxi =
      ModelWithControllersOnly_B.VectorConcatenate_d[0] *
      ModelWithControllersOnly_B.Integrator[3];
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Product: '<S159>/i x j' incorporates:
       *  Constant: '<S149>/Constant4'
       */
      ModelWithControllersOnly_B.ixj =
        ModelWithControllersOnly_cal->Constant4_Value_g *
        ModelWithControllersOnly_B.VectorConcatenate_d[1];
    }

    /* Product: '<S160>/k x j' */
    ModelWithControllersOnly_B.kxj =
      ModelWithControllersOnly_B.VectorConcatenate_d[1] *
      ModelWithControllersOnly_B.Integrator[3];
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Product: '<S160>/i x k' incorporates:
       *  Constant: '<S149>/Constant4'
       */
      ModelWithControllersOnly_B.ixk =
        ModelWithControllersOnly_cal->Constant4_Value_g *
        ModelWithControllersOnly_B.VectorConcatenate_d[2];

      /* Product: '<S160>/j x i' incorporates:
       *  Constant: '<S149>/Constant5'
       */
      ModelWithControllersOnly_B.jxi =
        ModelWithControllersOnly_cal->Constant5_Value *
        ModelWithControllersOnly_B.VectorConcatenate_d[0];
    }

    /* Sum: '<S157>/Sum' */
    ModelWithControllersOnly_B.Sum_m[0] = ModelWithControllersOnly_B.jxk -
      ModelWithControllersOnly_B.kxj;
    ModelWithControllersOnly_B.Sum_m[1] = ModelWithControllersOnly_B.kxi -
      ModelWithControllersOnly_B.ixk;
    ModelWithControllersOnly_B.Sum_m[2] = ModelWithControllersOnly_B.ixj -
      ModelWithControllersOnly_B.jxi;

    /* Sum: '<S150>/Add1' incorporates:
     *  Constant: '<S149>/Constant6'
     */
    ModelWithControllersOnly_B.Add1_f[0] = ModelWithControllersOnly_B.Sum_m[0] +
      ModelWithControllersOnly_B.Integrator[0];
    ModelWithControllersOnly_B.Add1_f[1] = ModelWithControllersOnly_B.Sum_m[1] +
      ModelWithControllersOnly_B.Integrator[1];
    ModelWithControllersOnly_B.Add1_f[2] = ModelWithControllersOnly_B.Sum_m[2] +
      ModelWithControllersOnly_cal->Constant6_Value;

    /* Reshape: '<S155>/Reshape1' */
    ModelWithControllersOnly_B.Reshape1_m[0] = ModelWithControllersOnly_B.Sum_m
      [0];
    ModelWithControllersOnly_B.Reshape1_m[1] = ModelWithControllersOnly_B.Sum_m
      [1];
    ModelWithControllersOnly_B.Reshape1_m[2] = ModelWithControllersOnly_B.Sum_m
      [2];

    /* Product: '<S155>/Product' incorporates:
     *  Math: '<S150>/Transpose1'
     *  Product: '<S156>/Product'
     *  Product: '<S164>/Product'
     *  Product: '<S165>/Product'
     *  Product: '<S174>/Product'
     *  Product: '<S175>/Product'
     *  Product: '<S182>/Product'
     *  Product: '<S183>/Product'
     *  Reshape: '<S155>/Reshape1'
     */
    lastU = &ModelWithControllersOnly_B.Transpose1[0];
    refPose[0] = ModelWithControllersOnly_B.Reshape1_m[0];
    refPose[1] = ModelWithControllersOnly_B.Reshape1_m[1];
    refPose[2] = ModelWithControllersOnly_B.Reshape1_m[2];
    for (i = 0; i < 3; i++) {
      /* Product: '<S155>/Product' */
      ModelWithControllersOnly_B.Product_m[i] = 0.0;
      ModelWithControllersOnly_B.Product_m[i] += lastU[i] * refPose[0];
      ModelWithControllersOnly_B.Product_m[i] += lastU[i + 3] * refPose[1];
      ModelWithControllersOnly_B.Product_m[i] += lastU[i + 6] * refPose[2];

      /* Reshape: '<S155>/Reshape2' incorporates:
       *  Product: '<S155>/Product'
       */
      ModelWithControllersOnly_B.Reshape2_p[i] =
        ModelWithControllersOnly_B.Product_m[i];
    }

    /* Sum: '<S150>/Add4' incorporates:
     *  Constant: '<S149>/Constant8'
     */
    ModelWithControllersOnly_B.V_wb[0] = ModelWithControllersOnly_B.y[0] +
      ModelWithControllersOnly_B.Reshape2_p[0];
    ModelWithControllersOnly_B.V_wb[1] = ModelWithControllersOnly_B.y[1] +
      ModelWithControllersOnly_B.Reshape2_p[1];
    ModelWithControllersOnly_B.V_wb[2] =
      ModelWithControllersOnly_cal->Constant8_Value +
      ModelWithControllersOnly_B.Reshape2_p[2];

    /* Trigonometry: '<S163>/sincos' */
    posError = ModelWithControllersOnly_B.TmpSignalConversionAtsincosInpo[0];
    u = std::sin(posError);
    posError = std::cos(posError);
    ModelWithControllersOnly_B.sincos_o1_h[0] = u;
    ModelWithControllersOnly_B.sincos_o2_g[0] = posError;
    posError = ModelWithControllersOnly_B.TmpSignalConversionAtsincosInpo[1];
    u = std::sin(posError);
    posError = std::cos(posError);
    ModelWithControllersOnly_B.sincos_o1_h[1] = u;
    ModelWithControllersOnly_B.sincos_o2_g[1] = posError;
    posError = ModelWithControllersOnly_B.TmpSignalConversionAtsincosInpo[2];
    u = std::sin(posError);
    posError = std::cos(posError);
    ModelWithControllersOnly_B.sincos_o1_h[2] = u;
    ModelWithControllersOnly_B.sincos_o2_g[2] = posError;

    /* Fcn: '<S163>/Fcn11' incorporates:
     *  Concatenate: '<S170>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_bj[0] =
      ModelWithControllersOnly_B.sincos_o2_g[0] *
      ModelWithControllersOnly_B.sincos_o2_g[1];

    /* Fcn: '<S163>/Fcn21' incorporates:
     *  Concatenate: '<S170>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_bj[1] =
      ModelWithControllersOnly_B.sincos_o1_h[1] *
      ModelWithControllersOnly_B.sincos_o1_h[2] *
      ModelWithControllersOnly_B.sincos_o2_g[0] -
      ModelWithControllersOnly_B.sincos_o1_h[0] *
      ModelWithControllersOnly_B.sincos_o2_g[2];

    /* Fcn: '<S163>/Fcn31' incorporates:
     *  Concatenate: '<S170>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_bj[2] =
      ModelWithControllersOnly_B.sincos_o1_h[1] *
      ModelWithControllersOnly_B.sincos_o2_g[2] *
      ModelWithControllersOnly_B.sincos_o2_g[0] +
      ModelWithControllersOnly_B.sincos_o1_h[0] *
      ModelWithControllersOnly_B.sincos_o1_h[2];

    /* Fcn: '<S163>/Fcn12' incorporates:
     *  Concatenate: '<S170>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_bj[3] =
      ModelWithControllersOnly_B.sincos_o1_h[0] *
      ModelWithControllersOnly_B.sincos_o2_g[1];

    /* Fcn: '<S163>/Fcn22' incorporates:
     *  Concatenate: '<S170>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_bj[4] =
      ModelWithControllersOnly_B.sincos_o1_h[1] *
      ModelWithControllersOnly_B.sincos_o1_h[2] *
      ModelWithControllersOnly_B.sincos_o1_h[0] +
      ModelWithControllersOnly_B.sincos_o2_g[0] *
      ModelWithControllersOnly_B.sincos_o2_g[2];

    /* Fcn: '<S163>/Fcn32' incorporates:
     *  Concatenate: '<S170>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_bj[5] =
      ModelWithControllersOnly_B.sincos_o1_h[1] *
      ModelWithControllersOnly_B.sincos_o2_g[2] *
      ModelWithControllersOnly_B.sincos_o1_h[0] -
      ModelWithControllersOnly_B.sincos_o2_g[0] *
      ModelWithControllersOnly_B.sincos_o1_h[2];

    /* Fcn: '<S163>/Fcn13' incorporates:
     *  Concatenate: '<S170>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_bj[6] =
      -ModelWithControllersOnly_B.sincos_o1_h[1];

    /* Fcn: '<S163>/Fcn23' incorporates:
     *  Concatenate: '<S170>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_bj[7] =
      ModelWithControllersOnly_B.sincos_o2_g[1] *
      ModelWithControllersOnly_B.sincos_o1_h[2];

    /* Fcn: '<S163>/Fcn33' incorporates:
     *  Concatenate: '<S170>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_bj[8] =
      ModelWithControllersOnly_B.sincos_o2_g[1] *
      ModelWithControllersOnly_B.sincos_o2_g[2];

    /* Reshape: '<S170>/Reshape (9) to [3x3] column-major' */
    std::memcpy(&ModelWithControllersOnly_B.Reshape9to3x3columnmajor_e[0],
                &ModelWithControllersOnly_B.VectorConcatenate_bj[0], 9U * sizeof
                (real_T));
    for (i = 0; i < 3; i++) {
      /* Math: '<S161>/Transpose1' incorporates:
       *  Reshape: '<S170>/Reshape (9) to [3x3] column-major'
       */
      ModelWithControllersOnly_B.Transpose1_a[3 * i] =
        ModelWithControllersOnly_B.Reshape9to3x3columnmajor_e[i];
      ModelWithControllersOnly_B.Transpose1_a[3 * i + 1] =
        ModelWithControllersOnly_B.Reshape9to3x3columnmajor_e[i + 3];
      ModelWithControllersOnly_B.Transpose1_a[3 * i + 2] =
        ModelWithControllersOnly_B.Reshape9to3x3columnmajor_e[i + 6];
    }

    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Constant: '<S151>/longOff' incorporates:
       *  Concatenate: '<S151>/Vector Concatenate'
       */
      ModelWithControllersOnly_B.VectorConcatenate_c[0] =
        ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_long;

      /* Constant: '<S151>/latOff' incorporates:
       *  Concatenate: '<S151>/Vector Concatenate'
       */
      ModelWithControllersOnly_B.VectorConcatenate_c[1] =
        ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_latO;

      /* Constant: '<S151>/vertOff' incorporates:
       *  Concatenate: '<S151>/Vector Concatenate'
       */
      ModelWithControllersOnly_B.VectorConcatenate_c[2] =
        ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_vert;

      /* Reshape: '<S165>/Reshape1' */
      ModelWithControllersOnly_B.Reshape1_e[0] =
        ModelWithControllersOnly_B.VectorConcatenate_c[0];
      ModelWithControllersOnly_B.Reshape1_e[1] =
        ModelWithControllersOnly_B.VectorConcatenate_c[1];
      ModelWithControllersOnly_B.Reshape1_e[2] =
        ModelWithControllersOnly_B.VectorConcatenate_c[2];
    }

    /* Product: '<S165>/Product' incorporates:
     *  Math: '<S161>/Transpose1'
     *  Product: '<S155>/Product'
     *  Product: '<S156>/Product'
     *  Product: '<S164>/Product'
     *  Product: '<S174>/Product'
     *  Product: '<S175>/Product'
     *  Product: '<S182>/Product'
     *  Product: '<S183>/Product'
     *  Reshape: '<S165>/Reshape1'
     */
    lastU = &ModelWithControllersOnly_B.Transpose1_a[0];
    refPose[0] = ModelWithControllersOnly_B.Reshape1_e[0];
    refPose[1] = ModelWithControllersOnly_B.Reshape1_e[1];
    refPose[2] = ModelWithControllersOnly_B.Reshape1_e[2];
    for (i = 0; i < 3; i++) {
      /* Product: '<S165>/Product' */
      ModelWithControllersOnly_B.Product_n[i] = 0.0;
      ModelWithControllersOnly_B.Product_n[i] += lastU[i] * refPose[0];
      ModelWithControllersOnly_B.Product_n[i] += lastU[i + 3] * refPose[1];
      ModelWithControllersOnly_B.Product_n[i] += lastU[i + 6] * refPose[2];

      /* Reshape: '<S165>/Reshape2' incorporates:
       *  Product: '<S165>/Product'
       */
      ModelWithControllersOnly_B.Reshape2_m[i] =
        ModelWithControllersOnly_B.Product_n[i];
    }

    /* Sum: '<S161>/Add' incorporates:
     *  Constant: '<S149>/Constant'
     */
    ModelWithControllersOnly_B.Add_p[0] =
      ModelWithControllersOnly_B.Integrator_m[0] +
      ModelWithControllersOnly_B.Reshape2_m[0];
    ModelWithControllersOnly_B.Add_p[1] =
      ModelWithControllersOnly_B.Integrator_m[1] +
      ModelWithControllersOnly_B.Reshape2_m[1];
    ModelWithControllersOnly_B.Add_p[2] =
      ModelWithControllersOnly_cal->Constant_Value_a +
      ModelWithControllersOnly_B.Reshape2_m[2];
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Product: '<S171>/j x k' incorporates:
       *  Constant: '<S149>/Constant5'
       */
      ModelWithControllersOnly_B.jxk_j =
        ModelWithControllersOnly_cal->Constant5_Value *
        ModelWithControllersOnly_B.VectorConcatenate_c[2];
    }

    /* Product: '<S171>/k x i' */
    ModelWithControllersOnly_B.kxi_g =
      ModelWithControllersOnly_B.VectorConcatenate_c[0] *
      ModelWithControllersOnly_B.Integrator[3];
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Product: '<S171>/i x j' incorporates:
       *  Constant: '<S149>/Constant4'
       */
      ModelWithControllersOnly_B.ixj_d =
        ModelWithControllersOnly_cal->Constant4_Value_g *
        ModelWithControllersOnly_B.VectorConcatenate_c[1];
    }

    /* Product: '<S172>/k x j' */
    ModelWithControllersOnly_B.kxj_j =
      ModelWithControllersOnly_B.VectorConcatenate_c[1] *
      ModelWithControllersOnly_B.Integrator[3];
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Product: '<S172>/i x k' incorporates:
       *  Constant: '<S149>/Constant4'
       */
      ModelWithControllersOnly_B.ixk_j =
        ModelWithControllersOnly_cal->Constant4_Value_g *
        ModelWithControllersOnly_B.VectorConcatenate_c[2];

      /* Product: '<S172>/j x i' incorporates:
       *  Constant: '<S149>/Constant5'
       */
      ModelWithControllersOnly_B.jxi_i =
        ModelWithControllersOnly_cal->Constant5_Value *
        ModelWithControllersOnly_B.VectorConcatenate_c[0];
    }

    /* Sum: '<S166>/Sum' */
    ModelWithControllersOnly_B.Sum_b[0] = ModelWithControllersOnly_B.jxk_j -
      ModelWithControllersOnly_B.kxj_j;
    ModelWithControllersOnly_B.Sum_b[1] = ModelWithControllersOnly_B.kxi_g -
      ModelWithControllersOnly_B.ixk_j;
    ModelWithControllersOnly_B.Sum_b[2] = ModelWithControllersOnly_B.ixj_d -
      ModelWithControllersOnly_B.jxi_i;

    /* Sum: '<S161>/Add1' incorporates:
     *  Constant: '<S149>/Constant6'
     */
    ModelWithControllersOnly_B.Add1_o[0] = ModelWithControllersOnly_B.Sum_b[0] +
      ModelWithControllersOnly_B.Integrator[0];
    ModelWithControllersOnly_B.Add1_o[1] = ModelWithControllersOnly_B.Sum_b[1] +
      ModelWithControllersOnly_B.Integrator[1];
    ModelWithControllersOnly_B.Add1_o[2] = ModelWithControllersOnly_B.Sum_b[2] +
      ModelWithControllersOnly_cal->Constant6_Value;

    /* Reshape: '<S164>/Reshape1' */
    ModelWithControllersOnly_B.Reshape1_j[0] = ModelWithControllersOnly_B.Sum_b
      [0];
    ModelWithControllersOnly_B.Reshape1_j[1] = ModelWithControllersOnly_B.Sum_b
      [1];
    ModelWithControllersOnly_B.Reshape1_j[2] = ModelWithControllersOnly_B.Sum_b
      [2];

    /* Product: '<S164>/Product' incorporates:
     *  Math: '<S161>/Transpose1'
     *  Product: '<S155>/Product'
     *  Product: '<S156>/Product'
     *  Product: '<S165>/Product'
     *  Product: '<S174>/Product'
     *  Product: '<S175>/Product'
     *  Product: '<S182>/Product'
     *  Product: '<S183>/Product'
     *  Reshape: '<S164>/Reshape1'
     */
    lastU = &ModelWithControllersOnly_B.Transpose1_a[0];
    refPose[0] = ModelWithControllersOnly_B.Reshape1_j[0];
    refPose[1] = ModelWithControllersOnly_B.Reshape1_j[1];
    refPose[2] = ModelWithControllersOnly_B.Reshape1_j[2];
    for (i = 0; i < 3; i++) {
      /* Product: '<S164>/Product' */
      ModelWithControllersOnly_B.Product_p[i] = 0.0;
      ModelWithControllersOnly_B.Product_p[i] += lastU[i] * refPose[0];
      ModelWithControllersOnly_B.Product_p[i] += lastU[i + 3] * refPose[1];
      ModelWithControllersOnly_B.Product_p[i] += lastU[i + 6] * refPose[2];

      /* Reshape: '<S164>/Reshape2' incorporates:
       *  Product: '<S164>/Product'
       */
      ModelWithControllersOnly_B.Reshape2_c[i] =
        ModelWithControllersOnly_B.Product_p[i];
    }

    /* Sum: '<S161>/Add4' incorporates:
     *  Constant: '<S149>/Constant8'
     */
    ModelWithControllersOnly_B.V_wb_p[0] = ModelWithControllersOnly_B.y[0] +
      ModelWithControllersOnly_B.Reshape2_c[0];
    ModelWithControllersOnly_B.V_wb_p[1] = ModelWithControllersOnly_B.y[1] +
      ModelWithControllersOnly_B.Reshape2_c[1];
    ModelWithControllersOnly_B.V_wb_p[2] =
      ModelWithControllersOnly_cal->Constant8_Value +
      ModelWithControllersOnly_B.Reshape2_c[2];

    /* RelationalOperator: '<S168>/Compare' incorporates:
     *  Constant: '<S168>/Constant'
     */
    posError = -ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_xdot;

    /* RelationalOperator: '<S168>/Compare' */
    ModelWithControllersOnly_B.Compare = (ModelWithControllersOnly_B.Add1_o[0] >=
      posError);

    /* RelationalOperator: '<S169>/Compare' incorporates:
     *  Constant: '<S169>/Constant'
     */
    ModelWithControllersOnly_B.Compare_f = (ModelWithControllersOnly_B.Add1_o[0]
      <= ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_xdot);

    /* Logic: '<S167>/Logical Operator' */
    ModelWithControllersOnly_B.LogicalOperator =
      (ModelWithControllersOnly_B.Compare &&
       ModelWithControllersOnly_B.Compare_f);

    /* Switch: '<S167>/Switch' */
    if (ModelWithControllersOnly_B.LogicalOperator) {
      /* Fcn: '<S167>/Fcn' */
      u = ModelWithControllersOnly_B.Add1_o[0] / 0.01;
      posError = rt_powd_snf(u, 2.0);

      /* Fcn: '<S167>/Fcn' */
      ModelWithControllersOnly_B.Fcn_g = 0.02 / (3.0 - posError);

      /* Switch: '<S167>/Switch' */
      ModelWithControllersOnly_B.Switch_d = ModelWithControllersOnly_B.Fcn_g;
    } else {
      /* Abs: '<S167>/Abs' */
      ModelWithControllersOnly_B.Abs_d = std::abs
        (ModelWithControllersOnly_B.Add1_o[0]);

      /* Switch: '<S167>/Switch' */
      ModelWithControllersOnly_B.Switch_d = ModelWithControllersOnly_B.Abs_d;
    }

    /* End of Switch: '<S167>/Switch' */

    /* Product: '<S162>/Divide' */
    ModelWithControllersOnly_B.Divide = ModelWithControllersOnly_B.Add1_o[1] /
      ModelWithControllersOnly_B.Switch_d;

    /* Trigonometry: '<S162>/Trigonometric Function' */
    ModelWithControllersOnly_B.Beta = std::atan
      (ModelWithControllersOnly_B.Divide);

    /* Trigonometry: '<S173>/sincos' */
    posError = ModelWithControllersOnly_B.TmpSignalConversionAtsincosInpo[0];
    u = std::sin(posError);
    posError = std::cos(posError);
    ModelWithControllersOnly_B.sincos_o1_n[0] = u;
    ModelWithControllersOnly_B.sincos_o2_k[0] = posError;
    posError = ModelWithControllersOnly_B.TmpSignalConversionAtsincosInpo[1];
    u = std::sin(posError);
    posError = std::cos(posError);
    ModelWithControllersOnly_B.sincos_o1_n[1] = u;
    ModelWithControllersOnly_B.sincos_o2_k[1] = posError;
    posError = ModelWithControllersOnly_B.TmpSignalConversionAtsincosInpo[2];
    u = std::sin(posError);
    posError = std::cos(posError);
    ModelWithControllersOnly_B.sincos_o1_n[2] = u;
    ModelWithControllersOnly_B.sincos_o2_k[2] = posError;

    /* Fcn: '<S173>/Fcn11' incorporates:
     *  Concatenate: '<S177>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_f5[0] =
      ModelWithControllersOnly_B.sincos_o2_k[0] *
      ModelWithControllersOnly_B.sincos_o2_k[1];

    /* Fcn: '<S173>/Fcn21' incorporates:
     *  Concatenate: '<S177>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_f5[1] =
      ModelWithControllersOnly_B.sincos_o1_n[1] *
      ModelWithControllersOnly_B.sincos_o1_n[2] *
      ModelWithControllersOnly_B.sincos_o2_k[0] -
      ModelWithControllersOnly_B.sincos_o1_n[0] *
      ModelWithControllersOnly_B.sincos_o2_k[2];

    /* Fcn: '<S173>/Fcn31' incorporates:
     *  Concatenate: '<S177>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_f5[2] =
      ModelWithControllersOnly_B.sincos_o1_n[1] *
      ModelWithControllersOnly_B.sincos_o2_k[2] *
      ModelWithControllersOnly_B.sincos_o2_k[0] +
      ModelWithControllersOnly_B.sincos_o1_n[0] *
      ModelWithControllersOnly_B.sincos_o1_n[2];

    /* Fcn: '<S173>/Fcn12' incorporates:
     *  Concatenate: '<S177>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_f5[3] =
      ModelWithControllersOnly_B.sincos_o1_n[0] *
      ModelWithControllersOnly_B.sincos_o2_k[1];

    /* Fcn: '<S173>/Fcn22' incorporates:
     *  Concatenate: '<S177>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_f5[4] =
      ModelWithControllersOnly_B.sincos_o1_n[1] *
      ModelWithControllersOnly_B.sincos_o1_n[2] *
      ModelWithControllersOnly_B.sincos_o1_n[0] +
      ModelWithControllersOnly_B.sincos_o2_k[0] *
      ModelWithControllersOnly_B.sincos_o2_k[2];

    /* Fcn: '<S173>/Fcn32' incorporates:
     *  Concatenate: '<S177>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_f5[5] =
      ModelWithControllersOnly_B.sincos_o1_n[1] *
      ModelWithControllersOnly_B.sincos_o2_k[2] *
      ModelWithControllersOnly_B.sincos_o1_n[0] -
      ModelWithControllersOnly_B.sincos_o2_k[0] *
      ModelWithControllersOnly_B.sincos_o1_n[2];

    /* Fcn: '<S173>/Fcn13' incorporates:
     *  Concatenate: '<S177>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_f5[6] =
      -ModelWithControllersOnly_B.sincos_o1_n[1];

    /* Fcn: '<S173>/Fcn23' incorporates:
     *  Concatenate: '<S177>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_f5[7] =
      ModelWithControllersOnly_B.sincos_o2_k[1] *
      ModelWithControllersOnly_B.sincos_o1_n[2];

    /* Fcn: '<S173>/Fcn33' incorporates:
     *  Concatenate: '<S177>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_f5[8] =
      ModelWithControllersOnly_B.sincos_o2_k[1] *
      ModelWithControllersOnly_B.sincos_o2_k[2];

    /* Reshape: '<S177>/Reshape (9) to [3x3] column-major' */
    std::memcpy(&ModelWithControllersOnly_B.Reshape9to3x3columnmajor_n[0],
                &ModelWithControllersOnly_B.VectorConcatenate_f5[0], 9U * sizeof
                (real_T));
    for (i = 0; i < 3; i++) {
      /* Math: '<S152>/Transpose1' incorporates:
       *  Reshape: '<S177>/Reshape (9) to [3x3] column-major'
       */
      ModelWithControllersOnly_B.Transpose1_k[3 * i] =
        ModelWithControllersOnly_B.Reshape9to3x3columnmajor_n[i];
      ModelWithControllersOnly_B.Transpose1_k[3 * i + 1] =
        ModelWithControllersOnly_B.Reshape9to3x3columnmajor_n[i + 3];
      ModelWithControllersOnly_B.Transpose1_k[3 * i + 2] =
        ModelWithControllersOnly_B.Reshape9to3x3columnmajor_n[i + 6];
    }

    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Constant: '<S152>/R_T1' incorporates:
       *  Concatenate: '<S152>/Vector Concatenate'
       */
      u = -b_b;
      ModelWithControllersOnly_B.VectorConcatenate_o[0] = u;

      /* Constant: '<S152>/R_T2' incorporates:
       *  Concatenate: '<S152>/Vector Concatenate'
       */
      ModelWithControllersOnly_B.VectorConcatenate_o[1] =
        ModelWithControllersOnly_cal->HardPointCoordinateTransformRea;

      /* Constant: '<S152>/R_T3' incorporates:
       *  Concatenate: '<S152>/Vector Concatenate'
       */
      ModelWithControllersOnly_B.VectorConcatenate_o[2] =
        ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_h;

      /* Reshape: '<S175>/Reshape1' */
      ModelWithControllersOnly_B.Reshape1_n[0] =
        ModelWithControllersOnly_B.VectorConcatenate_o[0];
      ModelWithControllersOnly_B.Reshape1_n[1] =
        ModelWithControllersOnly_B.VectorConcatenate_o[1];
      ModelWithControllersOnly_B.Reshape1_n[2] =
        ModelWithControllersOnly_B.VectorConcatenate_o[2];
    }

    /* Product: '<S175>/Product' incorporates:
     *  Math: '<S152>/Transpose1'
     *  Product: '<S155>/Product'
     *  Product: '<S156>/Product'
     *  Product: '<S164>/Product'
     *  Product: '<S165>/Product'
     *  Product: '<S174>/Product'
     *  Product: '<S182>/Product'
     *  Product: '<S183>/Product'
     *  Reshape: '<S175>/Reshape1'
     */
    lastU = &ModelWithControllersOnly_B.Transpose1_k[0];
    refPose[0] = ModelWithControllersOnly_B.Reshape1_n[0];
    refPose[1] = ModelWithControllersOnly_B.Reshape1_n[1];
    refPose[2] = ModelWithControllersOnly_B.Reshape1_n[2];
    for (i = 0; i < 3; i++) {
      /* Product: '<S175>/Product' */
      ModelWithControllersOnly_B.Product_g[i] = 0.0;
      ModelWithControllersOnly_B.Product_g[i] += lastU[i] * refPose[0];
      ModelWithControllersOnly_B.Product_g[i] += lastU[i + 3] * refPose[1];
      ModelWithControllersOnly_B.Product_g[i] += lastU[i + 6] * refPose[2];

      /* Reshape: '<S175>/Reshape2' incorporates:
       *  Product: '<S175>/Product'
       */
      ModelWithControllersOnly_B.Reshape2_g[i] =
        ModelWithControllersOnly_B.Product_g[i];
    }

    /* Sum: '<S152>/Add' incorporates:
     *  Constant: '<S149>/Constant'
     */
    ModelWithControllersOnly_B.Add_pa[0] =
      ModelWithControllersOnly_B.Integrator_m[0] +
      ModelWithControllersOnly_B.Reshape2_g[0];
    ModelWithControllersOnly_B.Add_pa[1] =
      ModelWithControllersOnly_B.Integrator_m[1] +
      ModelWithControllersOnly_B.Reshape2_g[1];
    ModelWithControllersOnly_B.Add_pa[2] =
      ModelWithControllersOnly_cal->Constant_Value_a +
      ModelWithControllersOnly_B.Reshape2_g[2];
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Product: '<S178>/j x k' incorporates:
       *  Constant: '<S149>/Constant5'
       */
      ModelWithControllersOnly_B.jxk_m =
        ModelWithControllersOnly_cal->Constant5_Value *
        ModelWithControllersOnly_B.VectorConcatenate_o[2];
    }

    /* Product: '<S178>/k x i' */
    ModelWithControllersOnly_B.kxi_k =
      ModelWithControllersOnly_B.VectorConcatenate_o[0] *
      ModelWithControllersOnly_B.Integrator[3];
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Product: '<S178>/i x j' incorporates:
       *  Constant: '<S149>/Constant4'
       */
      ModelWithControllersOnly_B.ixj_e =
        ModelWithControllersOnly_cal->Constant4_Value_g *
        ModelWithControllersOnly_B.VectorConcatenate_o[1];
    }

    /* Product: '<S179>/k x j' */
    ModelWithControllersOnly_B.kxj_b =
      ModelWithControllersOnly_B.VectorConcatenate_o[1] *
      ModelWithControllersOnly_B.Integrator[3];
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Product: '<S179>/i x k' incorporates:
       *  Constant: '<S149>/Constant4'
       */
      ModelWithControllersOnly_B.ixk_m =
        ModelWithControllersOnly_cal->Constant4_Value_g *
        ModelWithControllersOnly_B.VectorConcatenate_o[2];

      /* Product: '<S179>/j x i' incorporates:
       *  Constant: '<S149>/Constant5'
       */
      ModelWithControllersOnly_B.jxi_h =
        ModelWithControllersOnly_cal->Constant5_Value *
        ModelWithControllersOnly_B.VectorConcatenate_o[0];
    }

    /* Sum: '<S176>/Sum' */
    ModelWithControllersOnly_B.Sum_f[0] = ModelWithControllersOnly_B.jxk_m -
      ModelWithControllersOnly_B.kxj_b;
    ModelWithControllersOnly_B.Sum_f[1] = ModelWithControllersOnly_B.kxi_k -
      ModelWithControllersOnly_B.ixk_m;
    ModelWithControllersOnly_B.Sum_f[2] = ModelWithControllersOnly_B.ixj_e -
      ModelWithControllersOnly_B.jxi_h;

    /* Sum: '<S152>/Add1' incorporates:
     *  Constant: '<S149>/Constant6'
     */
    ModelWithControllersOnly_B.Add1_o0[0] = ModelWithControllersOnly_B.Sum_f[0]
      + ModelWithControllersOnly_B.Integrator[0];
    ModelWithControllersOnly_B.Add1_o0[1] = ModelWithControllersOnly_B.Sum_f[1]
      + ModelWithControllersOnly_B.Integrator[1];
    ModelWithControllersOnly_B.Add1_o0[2] = ModelWithControllersOnly_B.Sum_f[2]
      + ModelWithControllersOnly_cal->Constant6_Value;

    /* Reshape: '<S174>/Reshape1' */
    ModelWithControllersOnly_B.Reshape1_k[0] = ModelWithControllersOnly_B.Sum_f
      [0];
    ModelWithControllersOnly_B.Reshape1_k[1] = ModelWithControllersOnly_B.Sum_f
      [1];
    ModelWithControllersOnly_B.Reshape1_k[2] = ModelWithControllersOnly_B.Sum_f
      [2];

    /* Product: '<S174>/Product' incorporates:
     *  Math: '<S152>/Transpose1'
     *  Product: '<S155>/Product'
     *  Product: '<S156>/Product'
     *  Product: '<S164>/Product'
     *  Product: '<S165>/Product'
     *  Product: '<S175>/Product'
     *  Product: '<S182>/Product'
     *  Product: '<S183>/Product'
     *  Reshape: '<S174>/Reshape1'
     */
    lastU = &ModelWithControllersOnly_B.Transpose1_k[0];
    refPose[0] = ModelWithControllersOnly_B.Reshape1_k[0];
    refPose[1] = ModelWithControllersOnly_B.Reshape1_k[1];
    refPose[2] = ModelWithControllersOnly_B.Reshape1_k[2];
    for (i = 0; i < 3; i++) {
      /* Product: '<S174>/Product' */
      ModelWithControllersOnly_B.Product_ky[i] = 0.0;
      ModelWithControllersOnly_B.Product_ky[i] += lastU[i] * refPose[0];
      ModelWithControllersOnly_B.Product_ky[i] += lastU[i + 3] * refPose[1];
      ModelWithControllersOnly_B.Product_ky[i] += lastU[i + 6] * refPose[2];

      /* Reshape: '<S174>/Reshape2' incorporates:
       *  Product: '<S174>/Product'
       */
      ModelWithControllersOnly_B.Reshape2_f[i] =
        ModelWithControllersOnly_B.Product_ky[i];
    }

    /* Sum: '<S152>/Add4' incorporates:
     *  Constant: '<S149>/Constant8'
     */
    ModelWithControllersOnly_B.V_wb_g[0] = ModelWithControllersOnly_B.y[0] +
      ModelWithControllersOnly_B.Reshape2_f[0];
    ModelWithControllersOnly_B.V_wb_g[1] = ModelWithControllersOnly_B.y[1] +
      ModelWithControllersOnly_B.Reshape2_f[1];
    ModelWithControllersOnly_B.V_wb_g[2] =
      ModelWithControllersOnly_cal->Constant8_Value +
      ModelWithControllersOnly_B.Reshape2_f[2];
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Constant: '<S153>/Constant4' incorporates:
       *  Concatenate: '<S153>/Vector Concatenate'
       */
      ModelWithControllersOnly_B.VectorConcatenate_h[1] =
        ModelWithControllersOnly_cal->Constant4_Value_b;
    }

    /* Trigonometry: '<S181>/sincos' */
    posError = ModelWithControllersOnly_B.TmpSignalConversionAtsincosInpo[0];
    u = std::sin(posError);
    posError = std::cos(posError);
    ModelWithControllersOnly_B.sincos_o1_b[0] = u;
    ModelWithControllersOnly_B.sincos_o2_gp[0] = posError;
    posError = ModelWithControllersOnly_B.TmpSignalConversionAtsincosInpo[1];
    u = std::sin(posError);
    posError = std::cos(posError);
    ModelWithControllersOnly_B.sincos_o1_b[1] = u;
    ModelWithControllersOnly_B.sincos_o2_gp[1] = posError;
    posError = ModelWithControllersOnly_B.TmpSignalConversionAtsincosInpo[2];
    u = std::sin(posError);
    posError = std::cos(posError);
    ModelWithControllersOnly_B.sincos_o1_b[2] = u;
    ModelWithControllersOnly_B.sincos_o2_gp[2] = posError;

    /* Fcn: '<S181>/Fcn11' incorporates:
     *  Concatenate: '<S185>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_bx[0] =
      ModelWithControllersOnly_B.sincos_o2_gp[0] *
      ModelWithControllersOnly_B.sincos_o2_gp[1];

    /* Fcn: '<S181>/Fcn21' incorporates:
     *  Concatenate: '<S185>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_bx[1] =
      ModelWithControllersOnly_B.sincos_o1_b[1] *
      ModelWithControllersOnly_B.sincos_o1_b[2] *
      ModelWithControllersOnly_B.sincos_o2_gp[0] -
      ModelWithControllersOnly_B.sincos_o1_b[0] *
      ModelWithControllersOnly_B.sincos_o2_gp[2];

    /* Fcn: '<S181>/Fcn31' incorporates:
     *  Concatenate: '<S185>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_bx[2] =
      ModelWithControllersOnly_B.sincos_o1_b[1] *
      ModelWithControllersOnly_B.sincos_o2_gp[2] *
      ModelWithControllersOnly_B.sincos_o2_gp[0] +
      ModelWithControllersOnly_B.sincos_o1_b[0] *
      ModelWithControllersOnly_B.sincos_o1_b[2];

    /* Fcn: '<S181>/Fcn12' incorporates:
     *  Concatenate: '<S185>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_bx[3] =
      ModelWithControllersOnly_B.sincos_o1_b[0] *
      ModelWithControllersOnly_B.sincos_o2_gp[1];

    /* Fcn: '<S181>/Fcn22' incorporates:
     *  Concatenate: '<S185>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_bx[4] =
      ModelWithControllersOnly_B.sincos_o1_b[1] *
      ModelWithControllersOnly_B.sincos_o1_b[2] *
      ModelWithControllersOnly_B.sincos_o1_b[0] +
      ModelWithControllersOnly_B.sincos_o2_gp[0] *
      ModelWithControllersOnly_B.sincos_o2_gp[2];

    /* Fcn: '<S181>/Fcn32' incorporates:
     *  Concatenate: '<S185>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_bx[5] =
      ModelWithControllersOnly_B.sincos_o1_b[1] *
      ModelWithControllersOnly_B.sincos_o2_gp[2] *
      ModelWithControllersOnly_B.sincos_o1_b[0] -
      ModelWithControllersOnly_B.sincos_o2_gp[0] *
      ModelWithControllersOnly_B.sincos_o1_b[2];

    /* Fcn: '<S181>/Fcn13' incorporates:
     *  Concatenate: '<S185>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_bx[6] =
      -ModelWithControllersOnly_B.sincos_o1_b[1];

    /* Fcn: '<S181>/Fcn23' incorporates:
     *  Concatenate: '<S185>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_bx[7] =
      ModelWithControllersOnly_B.sincos_o2_gp[1] *
      ModelWithControllersOnly_B.sincos_o1_b[2];

    /* Fcn: '<S181>/Fcn33' incorporates:
     *  Concatenate: '<S185>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_bx[8] =
      ModelWithControllersOnly_B.sincos_o2_gp[1] *
      ModelWithControllersOnly_B.sincos_o2_gp[2];

    /* Reshape: '<S185>/Reshape (9) to [3x3] column-major' */
    std::memcpy(&ModelWithControllersOnly_B.Reshape9to3x3columnmajor_h[0],
                &ModelWithControllersOnly_B.VectorConcatenate_bx[0], 9U * sizeof
                (real_T));
    for (i = 0; i < 3; i++) {
      /* Math: '<S180>/Transpose1' incorporates:
       *  Reshape: '<S185>/Reshape (9) to [3x3] column-major'
       */
      ModelWithControllersOnly_B.Transpose1_m[3 * i] =
        ModelWithControllersOnly_B.Reshape9to3x3columnmajor_h[i];
      ModelWithControllersOnly_B.Transpose1_m[3 * i + 1] =
        ModelWithControllersOnly_B.Reshape9to3x3columnmajor_h[i + 3];
      ModelWithControllersOnly_B.Transpose1_m[3 * i + 2] =
        ModelWithControllersOnly_B.Reshape9to3x3columnmajor_h[i + 6];
    }

    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* UnaryMinus: '<S153>/Unary Minus2' incorporates:
       *  Concatenate: '<S153>/Vector Concatenate'
       *  Constant: '<S153>/Constant'
       */
      ModelWithControllersOnly_B.VectorConcatenate_h[0] =
        -ModelWithControllersOnly_cal->Constant_Value_dz;

      /* Sum: '<S153>/Sum' incorporates:
       *  Concatenate: '<S153>/Vector Concatenate'
       *  Constant: '<S153>/Constant2'
       *  Constant: '<S153>/Constant3'
       */
      ModelWithControllersOnly_B.VectorConcatenate_h[2] =
        ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_h -
        ModelWithControllersOnly_cal->Constant3_Value;

      /* Reshape: '<S183>/Reshape1' */
      ModelWithControllersOnly_B.Reshape1_f[0] =
        ModelWithControllersOnly_B.VectorConcatenate_h[0];
      ModelWithControllersOnly_B.Reshape1_f[1] =
        ModelWithControllersOnly_B.VectorConcatenate_h[1];
      ModelWithControllersOnly_B.Reshape1_f[2] =
        ModelWithControllersOnly_B.VectorConcatenate_h[2];
    }

    /* Product: '<S183>/Product' incorporates:
     *  Math: '<S180>/Transpose1'
     *  Product: '<S155>/Product'
     *  Product: '<S156>/Product'
     *  Product: '<S164>/Product'
     *  Product: '<S165>/Product'
     *  Product: '<S174>/Product'
     *  Product: '<S175>/Product'
     *  Product: '<S182>/Product'
     *  Reshape: '<S183>/Reshape1'
     */
    lastU = &ModelWithControllersOnly_B.Transpose1_m[0];
    refPose[0] = ModelWithControllersOnly_B.Reshape1_f[0];
    refPose[1] = ModelWithControllersOnly_B.Reshape1_f[1];
    refPose[2] = ModelWithControllersOnly_B.Reshape1_f[2];
    for (i = 0; i < 3; i++) {
      /* Product: '<S183>/Product' */
      ModelWithControllersOnly_B.Product_a[i] = 0.0;
      ModelWithControllersOnly_B.Product_a[i] += lastU[i] * refPose[0];
      ModelWithControllersOnly_B.Product_a[i] += lastU[i + 3] * refPose[1];
      ModelWithControllersOnly_B.Product_a[i] += lastU[i + 6] * refPose[2];

      /* Reshape: '<S183>/Reshape2' incorporates:
       *  Product: '<S183>/Product'
       */
      ModelWithControllersOnly_B.Reshape2_mj[i] =
        ModelWithControllersOnly_B.Product_a[i];
    }

    /* Sum: '<S180>/Add' incorporates:
     *  Constant: '<S149>/Constant'
     */
    ModelWithControllersOnly_B.Add_c[0] =
      ModelWithControllersOnly_B.Integrator_m[0] +
      ModelWithControllersOnly_B.Reshape2_mj[0];
    ModelWithControllersOnly_B.Add_c[1] =
      ModelWithControllersOnly_B.Integrator_m[1] +
      ModelWithControllersOnly_B.Reshape2_mj[1];
    ModelWithControllersOnly_B.Add_c[2] =
      ModelWithControllersOnly_cal->Constant_Value_a +
      ModelWithControllersOnly_B.Reshape2_mj[2];
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Product: '<S186>/j x k' incorporates:
       *  Constant: '<S149>/Constant5'
       */
      ModelWithControllersOnly_B.jxk_b =
        ModelWithControllersOnly_cal->Constant5_Value *
        ModelWithControllersOnly_B.VectorConcatenate_h[2];
    }

    /* Product: '<S186>/k x i' */
    ModelWithControllersOnly_B.kxi_e =
      ModelWithControllersOnly_B.VectorConcatenate_h[0] *
      ModelWithControllersOnly_B.Integrator[3];
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Product: '<S186>/i x j' incorporates:
       *  Constant: '<S149>/Constant4'
       */
      ModelWithControllersOnly_B.ixj_c =
        ModelWithControllersOnly_cal->Constant4_Value_g *
        ModelWithControllersOnly_B.VectorConcatenate_h[1];
    }

    /* Product: '<S187>/k x j' */
    ModelWithControllersOnly_B.kxj_f =
      ModelWithControllersOnly_B.VectorConcatenate_h[1] *
      ModelWithControllersOnly_B.Integrator[3];
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Product: '<S187>/i x k' incorporates:
       *  Constant: '<S149>/Constant4'
       */
      ModelWithControllersOnly_B.ixk_p =
        ModelWithControllersOnly_cal->Constant4_Value_g *
        ModelWithControllersOnly_B.VectorConcatenate_h[2];

      /* Product: '<S187>/j x i' incorporates:
       *  Constant: '<S149>/Constant5'
       */
      ModelWithControllersOnly_B.jxi_h1 =
        ModelWithControllersOnly_cal->Constant5_Value *
        ModelWithControllersOnly_B.VectorConcatenate_h[0];
    }

    /* Sum: '<S184>/Sum' */
    ModelWithControllersOnly_B.Sum_i[0] = ModelWithControllersOnly_B.jxk_b -
      ModelWithControllersOnly_B.kxj_f;
    ModelWithControllersOnly_B.Sum_i[1] = ModelWithControllersOnly_B.kxi_e -
      ModelWithControllersOnly_B.ixk_p;
    ModelWithControllersOnly_B.Sum_i[2] = ModelWithControllersOnly_B.ixj_c -
      ModelWithControllersOnly_B.jxi_h1;

    /* Sum: '<S180>/Add1' incorporates:
     *  Constant: '<S149>/Constant6'
     */
    ModelWithControllersOnly_B.Add1_h[0] = ModelWithControllersOnly_B.Sum_i[0] +
      ModelWithControllersOnly_B.Integrator[0];
    ModelWithControllersOnly_B.Add1_h[1] = ModelWithControllersOnly_B.Sum_i[1] +
      ModelWithControllersOnly_B.Integrator[1];
    ModelWithControllersOnly_B.Add1_h[2] = ModelWithControllersOnly_B.Sum_i[2] +
      ModelWithControllersOnly_cal->Constant6_Value;

    /* Reshape: '<S182>/Reshape1' */
    ModelWithControllersOnly_B.Reshape1_nq[0] =
      ModelWithControllersOnly_B.Sum_i[0];
    ModelWithControllersOnly_B.Reshape1_nq[1] =
      ModelWithControllersOnly_B.Sum_i[1];
    ModelWithControllersOnly_B.Reshape1_nq[2] =
      ModelWithControllersOnly_B.Sum_i[2];

    /* Product: '<S182>/Product' incorporates:
     *  Math: '<S180>/Transpose1'
     *  Product: '<S155>/Product'
     *  Product: '<S156>/Product'
     *  Product: '<S164>/Product'
     *  Product: '<S165>/Product'
     *  Product: '<S174>/Product'
     *  Product: '<S175>/Product'
     *  Product: '<S183>/Product'
     *  Reshape: '<S182>/Reshape1'
     */
    lastU = &ModelWithControllersOnly_B.Transpose1_m[0];
    refPose[0] = ModelWithControllersOnly_B.Reshape1_nq[0];
    refPose[1] = ModelWithControllersOnly_B.Reshape1_nq[1];
    refPose[2] = ModelWithControllersOnly_B.Reshape1_nq[2];
    for (i = 0; i < 3; i++) {
      /* Product: '<S182>/Product' */
      ModelWithControllersOnly_B.Product_iz[i] = 0.0;
      ModelWithControllersOnly_B.Product_iz[i] += lastU[i] * refPose[0];
      ModelWithControllersOnly_B.Product_iz[i] += lastU[i + 3] * refPose[1];
      ModelWithControllersOnly_B.Product_iz[i] += lastU[i + 6] * refPose[2];

      /* Reshape: '<S182>/Reshape2' incorporates:
       *  Product: '<S182>/Product'
       */
      ModelWithControllersOnly_B.Reshape2_e[i] =
        ModelWithControllersOnly_B.Product_iz[i];
    }

    /* Sum: '<S180>/Add4' incorporates:
     *  Constant: '<S149>/Constant8'
     */
    ModelWithControllersOnly_B.V_wb_e[0] = ModelWithControllersOnly_B.y[0] +
      ModelWithControllersOnly_B.Reshape2_e[0];
    ModelWithControllersOnly_B.V_wb_e[1] = ModelWithControllersOnly_B.y[1] +
      ModelWithControllersOnly_B.Reshape2_e[1];
    ModelWithControllersOnly_B.V_wb_e[2] =
      ModelWithControllersOnly_cal->Constant8_Value +
      ModelWithControllersOnly_B.Reshape2_e[2];

    /* Product: '<S148>/Product8' */
    ModelWithControllersOnly_B.Product8[0] = 0.0 *
      ModelWithControllersOnly_B.Add1_h[0];

    /* SignalConversion generated from: '<S148>/Vector Concatenate1' incorporates:
     *  Concatenate: '<S148>/Vector Concatenate1'
     */
    ModelWithControllersOnly_B.VectorConcatenate1_n[0] =
      ModelWithControllersOnly_B.Product8[0];

    /* Product: '<S148>/Product8' */
    ModelWithControllersOnly_B.Product8[1] = 0.0 *
      ModelWithControllersOnly_B.Add1_h[1];

    /* SignalConversion generated from: '<S148>/Vector Concatenate1' incorporates:
     *  Concatenate: '<S148>/Vector Concatenate1'
     */
    ModelWithControllersOnly_B.VectorConcatenate1_n[1] =
      ModelWithControllersOnly_B.Product8[1];

    /* Product: '<S148>/Product8' */
    ModelWithControllersOnly_B.Product8[2] = 0.0 *
      ModelWithControllersOnly_B.Add1_h[2];

    /* SignalConversion generated from: '<S148>/Vector Concatenate1' incorporates:
     *  Concatenate: '<S148>/Vector Concatenate1'
     */
    ModelWithControllersOnly_B.VectorConcatenate1_n[2] =
      ModelWithControllersOnly_B.Product8[2];

    /* Product: '<S148>/Product2' incorporates:
     *  Constant: '<S149>/Constant4'
     *  Constant: '<S149>/Constant5'
     */
    ModelWithControllersOnly_B.Product2_m[0] =
      ModelWithControllersOnly_cal->Constant4_Value_g * 0.0;
    ModelWithControllersOnly_B.Product2_m[1] =
      ModelWithControllersOnly_cal->Constant5_Value * 0.0;
    ModelWithControllersOnly_B.Product2_m[2] =
      ModelWithControllersOnly_B.Integrator[3] * 0.0;

    /* SignalConversion generated from: '<S148>/Vector Concatenate1' incorporates:
     *  Concatenate: '<S148>/Vector Concatenate1'
     */
    ModelWithControllersOnly_B.VectorConcatenate1_n[3] =
      ModelWithControllersOnly_B.Product2_m[0];
    ModelWithControllersOnly_B.VectorConcatenate1_n[4] =
      ModelWithControllersOnly_B.Product2_m[1];
    ModelWithControllersOnly_B.VectorConcatenate1_n[5] =
      ModelWithControllersOnly_B.Product2_m[2];

    /* Product: '<S148>/Product14' incorporates:
     *  Constant: '<S149>/Constant6'
     */
    ModelWithControllersOnly_B.Product14[0] = 0.0 *
      ModelWithControllersOnly_B.Integrator[0];
    ModelWithControllersOnly_B.Product14[1] = 0.0 *
      ModelWithControllersOnly_B.Integrator[1];
    ModelWithControllersOnly_B.Product14[2] = 0.0 *
      ModelWithControllersOnly_cal->Constant6_Value;

    /* SignalConversion generated from: '<S148>/Vector Concatenate2' incorporates:
     *  Concatenate: '<S148>/Vector Concatenate2'
     */
    ModelWithControllersOnly_B.VectorConcatenate2_o[0] =
      ModelWithControllersOnly_B.Product14[0];
    ModelWithControllersOnly_B.VectorConcatenate2_o[1] =
      ModelWithControllersOnly_B.Product14[1];
    ModelWithControllersOnly_B.VectorConcatenate2_o[2] =
      ModelWithControllersOnly_B.Product14[2];

    /* Product: '<S148>/Product3' incorporates:
     *  Constant: '<S149>/Constant4'
     *  Constant: '<S149>/Constant5'
     */
    ModelWithControllersOnly_B.Product3_b[0] = 0.0 *
      ModelWithControllersOnly_cal->Constant4_Value_g;
    ModelWithControllersOnly_B.Product3_b[1] = 0.0 *
      ModelWithControllersOnly_cal->Constant5_Value;
    ModelWithControllersOnly_B.Product3_b[2] = 0.0 *
      ModelWithControllersOnly_B.Integrator[3];

    /* SignalConversion generated from: '<S148>/Vector Concatenate2' incorporates:
     *  Concatenate: '<S148>/Vector Concatenate2'
     */
    ModelWithControllersOnly_B.VectorConcatenate2_o[3] =
      ModelWithControllersOnly_B.Product3_b[0];
    ModelWithControllersOnly_B.VectorConcatenate2_o[4] =
      ModelWithControllersOnly_B.Product3_b[1];
    ModelWithControllersOnly_B.VectorConcatenate2_o[5] =
      ModelWithControllersOnly_B.Product3_b[2];
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Product: '<S148>/Product' incorporates:
       *  Constant: '<S148>/Constant'
       */
      posError = b_m *
        ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_g;

      /* Product: '<S148>/Product' incorporates:
       *  Constant: '<S149>/Constant8'
       */
      ModelWithControllersOnly_B.Product_po = posError *
        ModelWithControllersOnly_cal->Constant8_Value;

      /* SignalConversion generated from: '<S205>/Vector Concatenate1' incorporates:
       *  Concatenate: '<S205>/Vector Concatenate1'
       *  Constant: '<S205>/Constant'
       */
      ModelWithControllersOnly_B.VectorConcatenate1_f[0] =
        ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_mu;

      /* SignalConversion generated from: '<S205>/Vector Concatenate1' incorporates:
       *  Concatenate: '<S205>/Vector Concatenate1'
       *  Constant: '<S205>/Constant'
       */
      ModelWithControllersOnly_B.VectorConcatenate1_f[1] =
        ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_mu;

      /* UnitConversion: '<S209>/Unit Conversion1' incorporates:
       *  Concatenate: '<S209>/Vector Concatenate1'
       *  Constant: '<S121>/Constant2'
       */
      /* Unit Conversion - from: N to: N
         Expression: output = (1*input) + (0) */
      ModelWithControllersOnly_B.VectorConcatenate1_c[0] =
        ModelWithControllersOnly_cal->Constant2_Value_g;

      /* SignalConversion generated from: '<S209>/Vector Concatenate1' incorporates:
       *  Concatenate: '<S209>/Vector Concatenate1'
       */
      ModelWithControllersOnly_B.VectorConcatenate1_c[1] = 0.0;
    }

    /* UnitConversion: '<S206>/Unit Conversion1' incorporates:
     *  Concatenate: '<S206>/Vector Concatenate1'
     */
    /* Unit Conversion - from: N to: N
       Expression: output = (1*input) + (0) */
    ModelWithControllersOnly_B.VectorConcatenate1_b2[0] =
      ModelWithControllersOnly_B.TnetNm;

    /* SignalConversion generated from: '<S206>/Vector Concatenate1' incorporates:
     *  Concatenate: '<S206>/Vector Concatenate1'
     */
    ModelWithControllersOnly_B.VectorConcatenate1_b2[1] = 0.0;

    /* SignalConversion generated from: '<S139>/ SFunction ' incorporates:
     *  MATLAB Function: '<S126>/vehicle model'
     */
    ModelWithControllersOnly_B.TmpSignalConversionAtSFunctionI[0] = 0.0;

    /* SignalConversion generated from: '<S139>/ SFunction ' incorporates:
     *  MATLAB Function: '<S126>/vehicle model'
     */
    ModelWithControllersOnly_B.TmpSignalConversionAtSFunctio_n[0] = 0.0;

    /* SignalConversion generated from: '<S139>/ SFunction ' incorporates:
     *  MATLAB Function: '<S126>/vehicle model'
     */
    ModelWithControllersOnly_B.TmpSignalConversionAtSFunctionI[1] = 0.0;

    /* SignalConversion generated from: '<S139>/ SFunction ' incorporates:
     *  MATLAB Function: '<S126>/vehicle model'
     */
    ModelWithControllersOnly_B.TmpSignalConversionAtSFunctio_n[1] = 0.0;

    /* SignalConversion generated from: '<S139>/ SFunction ' incorporates:
     *  MATLAB Function: '<S126>/vehicle model'
     */
    ModelWithControllersOnly_B.TmpSignalConversionAtSFunctionI[2] = 0.0;

    /* SignalConversion generated from: '<S139>/ SFunction ' incorporates:
     *  MATLAB Function: '<S126>/vehicle model'
     */
    ModelWithControllersOnly_B.TmpSignalConversionAtSFunctio_n[2] = 0.0;

    /* MATLAB Function: '<S126>/vehicle model' incorporates:
     *  Constant: '<S141>/Cyf'
     *  Constant: '<S141>/Cyr'
     */
    posError = ModelWithControllersOnly_B.UnitConversion1;
    u2 = tmp_8;
    b_h = ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_h;
    Nf = ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_NF;
    Nr = ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_NR;
    b_Izz = ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_Izz;
    b_g = ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_g;
    b_Fxtire_sat = ModelWithControllersOnly_cal->vehiclemodel_Fxtire_sat;
    b_Fytire_sat = ModelWithControllersOnly_cal->vehiclemodel_Fytire_sat;
    b_Fznom = ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_Fzno;
    Cy_f = ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_Cy_f;
    Cy_r = ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_Cy_r;
    xdot = ModelWithControllersOnly_B.Integrator[0];
    ydot = ModelWithControllersOnly_B.Integrator[1];
    r = ModelWithControllersOnly_B.Integrator[3];
    c_a = tmp_8 * ModelWithControllersOnly_B.Integrator[3] +
      ModelWithControllersOnly_B.Integrator[1];
    u = ModelWithControllersOnly_B.Integrator[0] *
      ModelWithControllersOnly_B.Integrator[0];
    d_a = ModelWithControllersOnly_B.Integrator[1] - b_b *
      ModelWithControllersOnly_B.Integrator[3];
    b_y = ModelWithControllersOnly_B.Integrator[0] *
      ModelWithControllersOnly_B.Integrator[0];
    alfa_f = std::abs(ModelWithControllersOnly_B.Integrator[0]);
    i = 0;
    if (alfa_f < ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_xdot)
    {
      i = 1;
    }

    loop_ub = i - 1;
    if (0 <= loop_ub) {
      x_data = alfa_f /
        ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_xdot;
    }

    if (0 <= i - 1) {
      x_data *= x_data;
    }

    loop_ub = i - 1;
    for (i = 0; i <= loop_ub; i++) {
      alfa_r = x_data;
      alfa_r = 2.0 *
        ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_xdot / (3.0 -
        alfa_r);
      x_data = alfa_r;
    }

    alfa_r = alfa_f;
    if (alfa_f < ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_xdot)
    {
      alfa_r = x_data;
    }

    alfa_f = rt_atan2d_snf(tmp_8 * ModelWithControllersOnly_B.Integrator[3] +
      ModelWithControllersOnly_B.Integrator[1], alfa_r) - std::tanh(4.0 *
      ModelWithControllersOnly_B.Integrator[0]) *
      ModelWithControllersOnly_B.UnitConversion1;
    alfa_r = rt_atan2d_snf(ModelWithControllersOnly_B.Integrator[1] - b_b *
      ModelWithControllersOnly_B.Integrator[3], alfa_r) - std::tanh(4.0 *
      ModelWithControllersOnly_B.Integrator[0]) * 0.0;
    FzCalc_idx_0 = 0.0;
    FzCalc_idx_1 = 0.0;
    for (i = 0; i < 6; i++) {
      if (i == 0) {
        B1 = ((((((ModelWithControllersOnly_B.Add1_d[0] +
                   ModelWithControllersOnly_B.TmpSignalConversionAtSFunctionI[0])
                  - (0.0 - r * 0.0) * b_m) - b_m * b_g * 0.0) * b_h -
                ModelWithControllersOnly_B.TmpSignalConversionAtSFunctio_n[1]) -
               ModelWithControllersOnly_B.TmpSignalConversionAtSFunctionI[2] *
               0.0) -
              ModelWithControllersOnly_B.TmpSignalConversionAtSFunctionI[0] *
              b_h) - ModelWithControllersOnly_B.Add_d[1];
        B2 = (-ModelWithControllersOnly_B.Add1_d[2] -
              ModelWithControllersOnly_B.TmpSignalConversionAtSFunctionI[2]) -
          b_m * b_g;
        FzCalc_idx_0 = -(B1 - B2 * b_b) / (u2 + b_b);
        FzCalc_idx_1 = (B2 * u2 + B1) / (u2 + b_b);
        FzCalc_idx_0 = -FzCalc_idx_0;
        d_idx_0 = FzCalc_idx_0;
        if (FzCalc_idx_0 < 0.0) {
          d_idx_0 = 0.0;
        }

        FzCalc_idx_1 = -FzCalc_idx_1;
        d_idx_1 = FzCalc_idx_1;
        if (FzCalc_idx_1 < 0.0) {
          d_idx_1 = 0.0;
        }
      } else {
        d_idx_0 = FzCalc_idx_0;
        d_idx_1 = FzCalc_idx_1;
      }

      B1 = -Cy_f * alfa_f * ModelWithControllersOnly_B.VectorConcatenate1_f[0] *
        d_idx_0 / b_Fznom;
      B2 = -Cy_r * alfa_r * ModelWithControllersOnly_B.VectorConcatenate1_f[1] *
        d_idx_1 / b_Fznom;
      ModelWithC_automlvehdynftiresat
        (ModelWithControllersOnly_B.VectorConcatenate1_b2[0], B1, b_Fxtire_sat *
         d_idx_0 / b_Fznom, b_Fytire_sat * d_idx_0 / b_Fznom, &Fx_ft, &Fy_ft);
      ModelWithC_automlvehdynftiresat
        (ModelWithControllersOnly_B.VectorConcatenate1_c[0], B2, b_Fxtire_sat *
         d_idx_1 / b_Fznom, b_Fytire_sat * d_idx_1 / b_Fznom, &Fx_rt, &Fy_rt);
      Fx_f = Fx_ft * std::cos(posError) - Fy_ft * std::sin(posError);
      B1 = -Fx_f * std::sin(posError) + B1 * std::cos(posError);
      Fx_r = Fx_rt - Fy_rt * 0.0;
      B2 += -Fx_r * 0.0;
      xddot = ((((Fx_f + Fx_r) - b_m * b_g * 0.0) +
                ModelWithControllersOnly_B.Add1_d[0]) +
               ModelWithControllersOnly_B.TmpSignalConversionAtSFunctionI[0]) /
        b_m + ydot * r;
      yddot = (((B1 + B2) + ModelWithControllersOnly_B.Add1_d[1]) +
               ModelWithControllersOnly_B.TmpSignalConversionAtSFunctionI[1]) /
        b_m + -xdot * r;
      rdot = ((((u2 * B1 - b_b * B2) + ModelWithControllersOnly_B.Add_d[2]) +
               ModelWithControllersOnly_B.TmpSignalConversionAtSFunctio_n[2]) -
              ModelWithControllersOnly_B.TmpSignalConversionAtSFunctionI[1] *
              0.0) / b_Izz;
      b_B1 = ((((((ModelWithControllersOnly_B.Add1_d[0] +
                   ModelWithControllersOnly_B.TmpSignalConversionAtSFunctionI[0])
                  - (xddot - r * ydot) * b_m) - b_m * b_g * 0.0) * b_h -
                ModelWithControllersOnly_B.TmpSignalConversionAtSFunctio_n[1]) -
               ModelWithControllersOnly_B.TmpSignalConversionAtSFunctionI[2] *
               0.0) -
              ModelWithControllersOnly_B.TmpSignalConversionAtSFunctionI[0] *
              b_h) - ModelWithControllersOnly_B.Add_d[1];
      b_B2 = (-ModelWithControllersOnly_B.Add1_d[2] -
              ModelWithControllersOnly_B.TmpSignalConversionAtSFunctionI[2]) -
        b_m * b_g;
      FzCalc_idx_0 = -(b_B1 - b_B2 * b_b) / (u2 + b_b);
      FzCalc_idx_1 = (b_B2 * u2 + b_B1) / (u2 + b_b);
      FzCalc_idx_0 = -FzCalc_idx_0;
      if (FzCalc_idx_0 < 0.0) {
        FzCalc_idx_0 = 0.0;
      }

      d_idx_0 = FzCalc_idx_0 - d_idx_0;
      maxFerr = std::abs(d_idx_0);
      FzCalc_idx_1 = -FzCalc_idx_1;
      if (FzCalc_idx_1 < 0.0) {
        FzCalc_idx_1 = 0.0;
      }

      d_idx_1 = FzCalc_idx_1 - d_idx_1;
      d_idx_0 = std::abs(d_idx_1);
      b_B1 = FzCalc_idx_0;
      b_B2 = FzCalc_idx_1;
      if ((maxFerr < d_idx_0) || (rtIsNaN(maxFerr) && (!rtIsNaN(d_idx_0)))) {
        maxFerr = d_idx_0;
      }
    }

    ModelWithControllersOnly_B.stateDer[0] = xddot;
    ModelWithControllersOnly_B.stateDer[1] = yddot;
    ModelWithControllersOnly_B.stateDer[2] = r;
    ModelWithControllersOnly_B.stateDer[3] = rdot;
    ModelWithControllersOnly_B.wheelInfo[0] = alfa_f;
    ModelWithControllersOnly_B.wheelInfo[1] = std::sqrt(c_a * c_a + u);
    ModelWithControllersOnly_B.wheelInfo[2] = alfa_r;
    ModelWithControllersOnly_B.wheelInfo[3] = std::sqrt(d_a * d_a + b_y);
    ModelWithControllersOnly_B.yOut[0] = xddot;
    ModelWithControllersOnly_B.yOut[1] = yddot;
    ModelWithControllersOnly_B.yOut[2] = rdot;
    ModelWithControllersOnly_B.FBody[0] = (xddot - ydot * r) * b_m;
    ModelWithControllersOnly_B.FBody[1] = (xdot * r + yddot) * b_m;
    ModelWithControllersOnly_B.FBody[2] = 0.0;
    ModelWithControllersOnly_B.MBody[0] = 0.0;
    ModelWithControllersOnly_B.MBody[1] = rdot * b_Izz;
    ModelWithControllersOnly_B.MBody[2] = 0.0;
    ModelWithControllersOnly_B.FTire[0] = Fx_ft / Nf;
    ModelWithControllersOnly_B.FTire[1] = Fy_ft / Nf;
    ModelWithControllersOnly_B.FTire[2] = b_B1 / Nf;
    ModelWithControllersOnly_B.FTire[3] = Fx_rt / Nr;
    ModelWithControllersOnly_B.FTire[4] = Fy_rt / Nr;
    ModelWithControllersOnly_B.FTire[5] = b_B2 / Nr;
    ModelWithControllersOnly_B.FOut[0] = Fx_f;
    ModelWithControllersOnly_B.FOut[1] = B1;
    ModelWithControllersOnly_B.FOut[2] = b_B1;
    ModelWithControllersOnly_B.FOut[3] = Fx_r;
    ModelWithControllersOnly_B.FOut[4] = B2;
    ModelWithControllersOnly_B.FOut[5] = b_B2;
    ModelWithControllersOnly_B.Fg[0] = 0.0;
    ModelWithControllersOnly_B.Fg[1] = 0.0;
    ModelWithControllersOnly_B.Fg[2] = b_m * b_g;
    ModelWithControllersOnly_B.status = maxFerr;

    /* Product: '<S148>/Product12' */
    ModelWithControllersOnly_B.Product12[0] = ModelWithControllersOnly_B.FOut[0]
      * ModelWithControllersOnly_B.Add1_f[0];
    ModelWithControllersOnly_B.Product12[1] = ModelWithControllersOnly_B.FOut[1]
      * ModelWithControllersOnly_B.Add1_f[1];
    ModelWithControllersOnly_B.Product12[2] = ModelWithControllersOnly_B.FOut[2]
      * ModelWithControllersOnly_B.Add1_f[2];

    /* Product: '<S148>/Product13' */
    ModelWithControllersOnly_B.Product13[0] =
      ModelWithControllersOnly_B.Add1_o0[0] * ModelWithControllersOnly_B.FOut[3];
    ModelWithControllersOnly_B.Product13[1] =
      ModelWithControllersOnly_B.Add1_o0[1] * ModelWithControllersOnly_B.FOut[4];
    ModelWithControllersOnly_B.Product13[2] =
      ModelWithControllersOnly_B.Add1_o0[2] * ModelWithControllersOnly_B.FOut[5];

    /* Product: '<S148>/Product15' incorporates:
     *  Concatenate: '<S148>/Vector Concatenate3'
     *  Constant: '<S149>/Constant6'
     */
    ModelWithControllersOnly_B.VectorConcatenate3_l[0] =
      ModelWithControllersOnly_B.Integrator[0] *
      ModelWithControllersOnly_B.UnaryMinus[0];
    ModelWithControllersOnly_B.VectorConcatenate3_l[1] =
      ModelWithControllersOnly_B.Integrator[1] *
      ModelWithControllersOnly_B.UnaryMinus[1];
    ModelWithControllersOnly_B.VectorConcatenate3_l[2] =
      ModelWithControllersOnly_cal->Constant6_Value *
      ModelWithControllersOnly_B.UnaryMinus[2];

    /* Product: '<S148>/Product4' incorporates:
     *  Concatenate: '<S148>/Vector Concatenate3'
     *  Constant: '<S149>/Constant4'
     *  Constant: '<S149>/Constant5'
     */
    ModelWithControllersOnly_B.VectorConcatenate3_l[3] =
      ModelWithControllersOnly_cal->Constant4_Value_g *
      ModelWithControllersOnly_B.UnaryMinus1[0];
    ModelWithControllersOnly_B.VectorConcatenate3_l[4] =
      ModelWithControllersOnly_cal->Constant5_Value *
      ModelWithControllersOnly_B.UnaryMinus1[1];
    ModelWithControllersOnly_B.VectorConcatenate3_l[5] =
      ModelWithControllersOnly_B.UnaryMinus1[2] *
      ModelWithControllersOnly_B.Integrator[3];

    /* Product: '<S198>/Divide1' */
    ModelWithControllersOnly_B.Divide1 = ModelWithControllersOnly_B.Integrator[1]
      * ModelWithControllersOnly_B.Integrator[3];

    /* Sum: '<S198>/Sum of Elements1' */
    ModelWithControllersOnly_B.ax = ModelWithControllersOnly_B.yOut[0] -
      ModelWithControllersOnly_B.Divide1;

    /* UnitConversion: '<S204>/Unit Conversion1' */
    /* Unit Conversion - from: m/s^2 to: gn
       Expression: output = (0.101972*input) + (0) */
    ModelWithControllersOnly_B.UnitConversion1_n = 0.10197162129779282 *
      ModelWithControllersOnly_B.ax;

    /* UnitConversion: '<S148>/Unit Conversion2' */
    /* Unit Conversion - from: gn to: m/s^2
       Expression: output = (9.80665*input) + (0) */
    ModelWithControllersOnly_B.UnitConversion2 = 9.8066500000000012 *
      ModelWithControllersOnly_B.UnitConversion1_n;

    /* Product: '<S148>/Product5' incorporates:
     *  Constant: '<S148>/Constant1'
     */
    ModelWithControllersOnly_B.Product5 = b_m *
      ModelWithControllersOnly_B.Integrator[0] *
      ModelWithControllersOnly_B.UnitConversion2;

    /* Product: '<S198>/Divide' */
    ModelWithControllersOnly_B.Divide_c = ModelWithControllersOnly_B.Integrator
      [0] * ModelWithControllersOnly_B.Integrator[3];

    /* Sum: '<S198>/Sum of Elements' */
    ModelWithControllersOnly_B.ay = ModelWithControllersOnly_B.yOut[1] +
      ModelWithControllersOnly_B.Divide_c;

    /* UnitConversion: '<S204>/Unit Conversion' */
    /* Unit Conversion - from: m/s^2 to: gn
       Expression: output = (0.101972*input) + (0) */
    ModelWithControllersOnly_B.UnitConversion = 0.10197162129779282 *
      ModelWithControllersOnly_B.ay;

    /* UnitConversion: '<S148>/Unit Conversion1' */
    /* Unit Conversion - from: gn to: m/s^2
       Expression: output = (9.80665*input) + (0) */
    ModelWithControllersOnly_B.UnitConversion1_f = 9.8066500000000012 *
      ModelWithControllersOnly_B.UnitConversion;

    /* Product: '<S148>/Product6' incorporates:
     *  Constant: '<S148>/Constant1'
     */
    ModelWithControllersOnly_B.Product6 = b_m *
      ModelWithControllersOnly_B.Integrator[1] *
      ModelWithControllersOnly_B.UnitConversion1_f;

    /* Product: '<S148>/Product7' incorporates:
     *  Constant: '<S148>/Constant2'
     */
    ModelWithControllersOnly_B.Product7 =
      ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_Izz *
      ModelWithControllersOnly_B.Integrator[3] *
      ModelWithControllersOnly_B.yOut[2];

    /* Sum: '<S148>/Sum of Elements1' */
    posError = -0.0;
    for (i = 0; i < 6; i++) {
      posError += ModelWithControllersOnly_B.VectorConcatenate1_n[i];
    }

    /* Sum: '<S148>/Sum of Elements1' */
    ModelWithControllersOnly_B.SumofElements1 = posError;

    /* Sum: '<S148>/Sum of Elements2' */
    posError = -0.0;
    for (i = 0; i < 6; i++) {
      posError += ModelWithControllersOnly_B.VectorConcatenate2_o[i];
    }

    /* Sum: '<S148>/Sum of Elements2' */
    ModelWithControllersOnly_B.SumofElements2 = posError;

    /* Sum: '<S148>/Sum of Elements3' */
    posError = -0.0;
    for (i = 0; i < 6; i++) {
      posError += ModelWithControllersOnly_B.VectorConcatenate3_l[i];
    }

    /* Sum: '<S148>/Sum of Elements3' */
    ModelWithControllersOnly_B.SumofElements3 = posError;
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* UnaryMinus: '<S148>/Unary Minus5' */
      ModelWithControllersOnly_B.UnaryMinus5 =
        -ModelWithControllersOnly_B.Product_po;
    }

    /* Sum: '<S193>/Sum1' */
    ModelWithControllersOnly_B.Sum1[0] = ModelWithControllersOnly_B.Product14[0]
      + ModelWithControllersOnly_B.Product8[0];
    ModelWithControllersOnly_B.Sum1[1] = ModelWithControllersOnly_B.Product14[1]
      + ModelWithControllersOnly_B.Product8[1];

    /* Sum: '<S193>/Sum' */
    ModelWithControllersOnly_B.Sum_e = ModelWithControllersOnly_B.Product3_b[2]
      + ModelWithControllersOnly_B.Product2_m[2];

    /* RelationalOperator: '<S201>/Compare' incorporates:
     *  Constant: '<S201>/Constant'
     */
    posError = -ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_xdot;

    /* RelationalOperator: '<S201>/Compare' */
    ModelWithControllersOnly_B.Compare_i =
      (ModelWithControllersOnly_B.Integrator[0] >= posError);

    /* RelationalOperator: '<S202>/Compare' incorporates:
     *  Constant: '<S202>/Constant'
     */
    ModelWithControllersOnly_B.Compare_p =
      (ModelWithControllersOnly_B.Integrator[0] <=
       ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_xdot);

    /* Logic: '<S200>/Logical Operator' */
    ModelWithControllersOnly_B.LogicalOperator_e =
      (ModelWithControllersOnly_B.Compare_i &&
       ModelWithControllersOnly_B.Compare_p);

    /* Switch: '<S200>/Switch' */
    if (ModelWithControllersOnly_B.LogicalOperator_e) {
      /* Fcn: '<S200>/Fcn' */
      u = ModelWithControllersOnly_B.Integrator[0] / 0.01;
      posError = rt_powd_snf(u, 2.0);

      /* Fcn: '<S200>/Fcn' */
      ModelWithControllersOnly_B.Fcn = 0.02 / (3.0 - posError);

      /* Switch: '<S200>/Switch' */
      ModelWithControllersOnly_B.Switch_p = ModelWithControllersOnly_B.Fcn;
    } else {
      /* Abs: '<S200>/Abs' */
      ModelWithControllersOnly_B.Abs = std::abs
        (ModelWithControllersOnly_B.Integrator[0]);

      /* Switch: '<S200>/Switch' */
      ModelWithControllersOnly_B.Switch_p = ModelWithControllersOnly_B.Abs;
    }

    /* End of Switch: '<S200>/Switch' */

    /* Product: '<S195>/Divide' */
    ModelWithControllersOnly_B.Divide_p = ModelWithControllersOnly_B.Integrator
      [1] / ModelWithControllersOnly_B.Switch_p;

    /* Trigonometry: '<S195>/Trigonometric Function' */
    ModelWithControllersOnly_B.Beta_n = std::atan
      (ModelWithControllersOnly_B.Divide_p);
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Constant: '<S149>/Constant1' */
      ModelWithControllersOnly_B.az =
        ModelWithControllersOnly_cal->Constant1_Value_o;

      /* Constant: '<S149>/Constant10' */
      ModelWithControllersOnly_B.Constant10 =
        ModelWithControllersOnly_cal->Constant10_Value;

      /* Constant: '<S149>/Constant3' */
      ModelWithControllersOnly_B.zddot =
        ModelWithControllersOnly_cal->Constant3_Value_j;

      /* Constant: '<S149>/Constant9' */
      ModelWithControllersOnly_B.Constant9 =
        ModelWithControllersOnly_cal->Constant9_Value_n;
    }

    /* Trigonometry: '<S197>/sincos' */
    posError = ModelWithControllersOnly_B.TmpSignalConversionAtsincosInpo[0];
    u = std::sin(posError);
    posError = std::cos(posError);
    ModelWithControllersOnly_B.sincos_o1_bs[0] = u;
    ModelWithControllersOnly_B.sincos_o2_f[0] = posError;
    posError = ModelWithControllersOnly_B.TmpSignalConversionAtsincosInpo[1];
    u = std::sin(posError);
    posError = std::cos(posError);
    ModelWithControllersOnly_B.sincos_o1_bs[1] = u;
    ModelWithControllersOnly_B.sincos_o2_f[1] = posError;
    posError = ModelWithControllersOnly_B.TmpSignalConversionAtsincosInpo[2];
    u = std::sin(posError);
    posError = std::cos(posError);
    ModelWithControllersOnly_B.sincos_o1_bs[2] = u;
    ModelWithControllersOnly_B.sincos_o2_f[2] = posError;

    /* Fcn: '<S197>/Fcn11' incorporates:
     *  Concatenate: '<S203>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_jc[0] =
      ModelWithControllersOnly_B.sincos_o2_f[0] *
      ModelWithControllersOnly_B.sincos_o2_f[1];

    /* Fcn: '<S197>/Fcn21' incorporates:
     *  Concatenate: '<S203>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_jc[1] =
      ModelWithControllersOnly_B.sincos_o1_bs[1] *
      ModelWithControllersOnly_B.sincos_o1_bs[2] *
      ModelWithControllersOnly_B.sincos_o2_f[0] -
      ModelWithControllersOnly_B.sincos_o1_bs[0] *
      ModelWithControllersOnly_B.sincos_o2_f[2];

    /* Fcn: '<S197>/Fcn31' incorporates:
     *  Concatenate: '<S203>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_jc[2] =
      ModelWithControllersOnly_B.sincos_o1_bs[1] *
      ModelWithControllersOnly_B.sincos_o2_f[2] *
      ModelWithControllersOnly_B.sincos_o2_f[0] +
      ModelWithControllersOnly_B.sincos_o1_bs[0] *
      ModelWithControllersOnly_B.sincos_o1_bs[2];

    /* Fcn: '<S197>/Fcn12' incorporates:
     *  Concatenate: '<S203>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_jc[3] =
      ModelWithControllersOnly_B.sincos_o1_bs[0] *
      ModelWithControllersOnly_B.sincos_o2_f[1];

    /* Fcn: '<S197>/Fcn22' incorporates:
     *  Concatenate: '<S203>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_jc[4] =
      ModelWithControllersOnly_B.sincos_o1_bs[1] *
      ModelWithControllersOnly_B.sincos_o1_bs[2] *
      ModelWithControllersOnly_B.sincos_o1_bs[0] +
      ModelWithControllersOnly_B.sincos_o2_f[0] *
      ModelWithControllersOnly_B.sincos_o2_f[2];

    /* Fcn: '<S197>/Fcn32' incorporates:
     *  Concatenate: '<S203>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_jc[5] =
      ModelWithControllersOnly_B.sincos_o1_bs[1] *
      ModelWithControllersOnly_B.sincos_o2_f[2] *
      ModelWithControllersOnly_B.sincos_o1_bs[0] -
      ModelWithControllersOnly_B.sincos_o2_f[0] *
      ModelWithControllersOnly_B.sincos_o1_bs[2];

    /* Fcn: '<S197>/Fcn13' incorporates:
     *  Concatenate: '<S203>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_jc[6] =
      -ModelWithControllersOnly_B.sincos_o1_bs[1];

    /* Fcn: '<S197>/Fcn23' incorporates:
     *  Concatenate: '<S203>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_jc[7] =
      ModelWithControllersOnly_B.sincos_o2_f[1] *
      ModelWithControllersOnly_B.sincos_o1_bs[2];

    /* Fcn: '<S197>/Fcn33' incorporates:
     *  Concatenate: '<S203>/Vector Concatenate'
     */
    ModelWithControllersOnly_B.VectorConcatenate_jc[8] =
      ModelWithControllersOnly_B.sincos_o2_f[1] *
      ModelWithControllersOnly_B.sincos_o2_f[2];

    /* Reshape: '<S203>/Reshape (9) to [3x3] column-major' */
    std::memcpy(&ModelWithControllersOnly_B.Reshape9to3x3columnmajor_l[0],
                &ModelWithControllersOnly_B.VectorConcatenate_jc[0], 9U * sizeof
                (real_T));

    /* Integrator: '<S212>/lateral' */
    ModelWithControllersOnly_B.lateral =
      ModelWithControllersOnly_X.lateral_CSTATE;

    /* Sum: '<S212>/Sum' */
    ModelWithControllersOnly_B.Sum_c = ModelWithControllersOnly_B.wheelInfo[0] -
      ModelWithControllersOnly_B.lateral;

    /* Product: '<S212>/Product1' incorporates:
     *  Constant: '<S211>/Constant1'
     */
    ModelWithControllersOnly_B.Product1_h = ModelWithControllersOnly_B.Sum_c *
      ModelWithControllersOnly_B.wheelInfo[1] /
      ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_sigm;

    /* Integrator: '<S213>/lateral' */
    ModelWithControllersOnly_B.lateral_g =
      ModelWithControllersOnly_X.lateral_CSTATE_k;

    /* Sum: '<S213>/Sum' */
    ModelWithControllersOnly_B.Sum_g = ModelWithControllersOnly_B.wheelInfo[2] -
      ModelWithControllersOnly_B.lateral_g;

    /* Product: '<S213>/Product1' incorporates:
     *  Constant: '<S211>/Constant2'
     */
    ModelWithControllersOnly_B.Product1_af = ModelWithControllersOnly_B.Sum_g *
      ModelWithControllersOnly_B.wheelInfo[3] /
      ModelWithControllersOnly_cal->VehicleBody3DOFSingleTrack_si_n;

    /* Saturate: '<S127>/Throttle Saturation' */
    u = ModelWithControllersOnly_B.Switch;
    posError = tmp_b;
    u2 = tmp_c;
    if (u > u2) {
      /* Saturate: '<S127>/Throttle Saturation' */
      ModelWithControllersOnly_B.ThrottleSaturation_p = u2;
    } else if (u < posError) {
      /* Saturate: '<S127>/Throttle Saturation' */
      ModelWithControllersOnly_B.ThrottleSaturation_p = posError;
    } else {
      /* Saturate: '<S127>/Throttle Saturation' */
      ModelWithControllersOnly_B.ThrottleSaturation_p = u;
    }

    /* End of Saturate: '<S127>/Throttle Saturation' */
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* DiscreteTransferFcn: '<S127>/second-order low-pass torque model' */
      B1 = ModelWithControllersOnly_B.ThrottleSaturation_p;
      B1 -= ModelWithControllersOnly_DW.secondorderlowpasstorquemodel_s[0] *
        tmp[1];
      B1 -= ModelWithControllersOnly_DW.secondorderlowpasstorquemodel_s[1] *
        tmp[2];
      B1 /= tmp[0];
      ModelWithControllersOnly_DW.secondorderlowpasstorquemodel_t = B1;
      B1 = tmp_0[0] *
        ModelWithControllersOnly_DW.secondorderlowpasstorquemodel_t;
      B1 += ModelWithControllersOnly_DW.secondorderlowpasstorquemodel_s[0] *
        tmp_0[1];
      B1 += ModelWithControllersOnly_DW.secondorderlowpasstorquemodel_s[1] *
        tmp_0[2];

      /* DiscreteTransferFcn: '<S127>/second-order low-pass torque model' */
      ModelWithControllersOnly_B.TengineNm = B1;

      /* DeadZone: '<S127>/Dead Zone1' */
      if (ModelWithControllersOnly_B.TengineNm > tmp_7) {
        /* DeadZone: '<S127>/Dead Zone1' */
        ModelWithControllersOnly_B.DeadZone1 =
          ModelWithControllersOnly_B.TengineNm - tmp_7;
      } else {
        posError = -tmp_7;
        if (ModelWithControllersOnly_B.TengineNm >= posError) {
          /* DeadZone: '<S127>/Dead Zone1' */
          ModelWithControllersOnly_B.DeadZone1 = 0.0;
        } else {
          posError = -tmp_7;

          /* DeadZone: '<S127>/Dead Zone1' */
          ModelWithControllersOnly_B.DeadZone1 =
            ModelWithControllersOnly_B.TengineNm - posError;
        }
      }

      /* End of DeadZone: '<S127>/Dead Zone1' */
    }

    /* Saturate: '<S127>/Brake Saturation' */
    u = ModelWithControllersOnly_B.Switch1;
    posError = tmp_9;
    u2 = tmp_a;
    if (u > u2) {
      /* Saturate: '<S127>/Brake Saturation' */
      ModelWithControllersOnly_B.BrakeSaturation_m = u2;
    } else if (u < posError) {
      /* Saturate: '<S127>/Brake Saturation' */
      ModelWithControllersOnly_B.BrakeSaturation_m = posError;
    } else {
      /* Saturate: '<S127>/Brake Saturation' */
      ModelWithControllersOnly_B.BrakeSaturation_m = u;
    }

    /* End of Saturate: '<S127>/Brake Saturation' */
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* DiscreteTransferFcn: '<S127>/second-order low-pass brake model' */
      B1 = ModelWithControllersOnly_B.BrakeSaturation_m;
      B1 -= ModelWithControllersOnly_DW.secondorderlowpassbrakemodel_st[0] *
        tmp_1[1];
      B1 -= ModelWithControllersOnly_DW.secondorderlowpassbrakemodel_st[1] *
        tmp_1[2];
      B1 /= tmp_1[0];
      ModelWithControllersOnly_DW.secondorderlowpassbrakemodel_tm = B1;
      B1 = tmp_2[0] *
        ModelWithControllersOnly_DW.secondorderlowpassbrakemodel_tm;
      B1 += ModelWithControllersOnly_DW.secondorderlowpassbrakemodel_st[0] *
        tmp_2[1];
      B1 += ModelWithControllersOnly_DW.secondorderlowpassbrakemodel_st[1] *
        tmp_2[2];

      /* DiscreteTransferFcn: '<S127>/second-order low-pass brake model' */
      ModelWithControllersOnly_B.TbrakeNm = B1;

      /* UnitDelay: '<S127>/Unit Delay' */
      ModelWithControllersOnly_B.UnitDelay_p =
        ModelWithControllersOnly_DW.UnitDelay_DSTATE_l;

      /* DeadZone: '<S127>/Dead Zone' */
      if (ModelWithControllersOnly_B.UnitDelay_p > tmp_6) {
        /* DeadZone: '<S127>/Dead Zone' */
        ModelWithControllersOnly_B.DeadZone =
          ModelWithControllersOnly_B.UnitDelay_p - tmp_6;
      } else {
        posError = -tmp_6;
        if (ModelWithControllersOnly_B.UnitDelay_p >= posError) {
          /* DeadZone: '<S127>/Dead Zone' */
          ModelWithControllersOnly_B.DeadZone = 0.0;
        } else {
          posError = -tmp_6;

          /* DeadZone: '<S127>/Dead Zone' */
          ModelWithControllersOnly_B.DeadZone =
            ModelWithControllersOnly_B.UnitDelay_p - posError;
        }
      }

      /* End of DeadZone: '<S127>/Dead Zone' */

      /* MATLAB Function: '<S127>/compute net torque' */
      B1 = ModelWithControllersOnly_B.DeadZone1;
      if (ModelWithControllersOnly_B.gear_cmd == 0.0) {
        B1 = 0.0;
      }

      ModelWithControllersOnly_B.Tnet = std::tanh(5.0 *
        ModelWithControllersOnly_B.DeadZone) *
        ModelWithControllersOnly_B.TbrakeNm + B1;

      /* End of MATLAB Function: '<S127>/compute net torque' */
    }

    /* Saturate: '<S127>/Steering Saturation' */
    u = ModelWithControllersOnly_B.Gain_n;
    posError = *get_MIN_STEER();
    u2 = *get_MAX_STEER();
    if (u > u2) {
      /* Saturate: '<S127>/Steering Saturation' */
      ModelWithControllersOnly_B.SteeringSaturation = u2;
    } else if (u < posError) {
      /* Saturate: '<S127>/Steering Saturation' */
      ModelWithControllersOnly_B.SteeringSaturation = posError;
    } else {
      /* Saturate: '<S127>/Steering Saturation' */
      ModelWithControllersOnly_B.SteeringSaturation = u;
    }

    /* End of Saturate: '<S127>/Steering Saturation' */
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* DiscreteTransferFcn: '<S127>/second-order low-pass steering model' */
      B1 = ModelWithControllersOnly_B.SteeringSaturation;
      B1 -= ModelWithControllersOnly_DW.secondorderlowpasssteeringmodel[0] *
        tmp_3[1];
      B1 -= ModelWithControllersOnly_DW.secondorderlowpasssteeringmodel[1] *
        tmp_3[2];
      B1 /= tmp_3[0];
      ModelWithControllersOnly_DW.secondorderlowpasssteeringmod_c = B1;
      B1 = tmp_4[0] *
        ModelWithControllersOnly_DW.secondorderlowpasssteeringmod_c;
      B1 += ModelWithControllersOnly_DW.secondorderlowpasssteeringmodel[0] *
        tmp_4[1];
      B1 += ModelWithControllersOnly_DW.secondorderlowpasssteeringmodel[1] *
        tmp_4[2];

      /* DiscreteTransferFcn: '<S127>/second-order low-pass steering model' */
      ModelWithControllersOnly_B.steeringwheelangledeg = B1;

      /* Gain: '<S127>/deg-to-rad1' */
      ModelWithControllersOnly_B.steeringwheelanglerad =
        ModelWithControllersOnly_cal->degtorad1_Gain *
        ModelWithControllersOnly_B.steeringwheelangledeg;

      /* Lookup_n-D: '<S218>/1-D Lookup Table' */
      /*
       * About '<S218>/1-D Lookup Table':
       *       Table size:  2
       *    Interpolation:  Spline
       *    Extrapolation:  Linear
       *   Breakpt Search:  Binary
       *    Breakpt Cache:  OFF
       */
      ModelWithControllersOnly_B.uDLookupTable = look_SplNBinXZcd(1U, (
        const_cast<real_T*>(&ModelWithControllersOnly_RGND)), (rt_LUTSplineWork*)
        &ModelWithControllersOnly_DW.SWork[0]);

      /* Product: '<S218>/Product' */
      ModelWithControllersOnly_B.Product_kk =
        ModelWithControllersOnly_B.steeringwheelanglerad *
        ModelWithControllersOnly_B.uDLookupTable;

      /* Lookup_n-D: '<S219>/1-D Lookup Table1' incorporates:
       *  Product: '<S218>/Product'
       */
      ModelWithControllersOnly_B.uDLookupTable1 = look1_binlxpw
        (ModelWithControllersOnly_B.Product_kk,
         ModelWithControllersOnly_cal->MappedSteering_StrgAngBpts,
         ModelWithControllersOnly_cal->uDLookupTable1_tableData, 1U);

      /* Lookup_n-D: '<S219>/1-D Lookup Table' incorporates:
       *  Product: '<S218>/Product'
       */
      ModelWithControllersOnly_B.uDLookupTable_c = look1_binlxpw
        (ModelWithControllersOnly_B.Product_kk,
         ModelWithControllersOnly_cal->MappedSteering_StrgAngBpts,
         ModelWithControllersOnly_cal->uDLookupTable_tableData, 1U);

      /* Sum: '<S127>/Add' */
      ModelWithControllersOnly_B.Add_p3 =
        ModelWithControllersOnly_B.uDLookupTable1 +
        ModelWithControllersOnly_B.uDLookupTable_c;

      /* Gain: '<S127>/Gain2' */
      ModelWithControllersOnly_B.Gain2_o =
        ModelWithControllersOnly_cal->Gain2_Gain_h0 *
        ModelWithControllersOnly_B.Add_p3;

      /* Gain: '<S127>/rolling resistance' */
      ModelWithControllersOnly_B.rollingresistance = *get_RRdamp() *
        ModelWithControllersOnly_B.UnitDelay_p;

      /* Sum: '<S127>/Sum' */
      ModelWithControllersOnly_B.TnetNm_e = ModelWithControllersOnly_B.Tnet -
        ModelWithControllersOnly_B.rollingresistance;

      /* Gain: '<S127>/torque to acceleration conversion1' */
      posError = 1.0 / tmp_5;

      /* Gain: '<S127>/torque to acceleration conversion1' */
      ModelWithControllersOnly_B.torquetoaccelerationconversion1 = posError *
        ModelWithControllersOnly_B.TnetNm_e;

      /* Constant: '<S5>/Constant3' */
      ModelWithControllersOnly_B.gear_cmd_e =
        ModelWithControllersOnly_cal->Constant3_Value_g;

      /* Constant: '<S5>/Constant4' */
      ModelWithControllersOnly_B.gear_cmd_i =
        ModelWithControllersOnly_cal->Constant4_Value_ad;

      /* FromWorkspace: '<S5>/From Workspace' */
      {
        real_T *pDataValues = (real_T *)
          ModelWithControllersOnly_DW.FromWorkspace_PWORK.DataPtr;
        real_T *pTimeValues = (real_T *)
          ModelWithControllersOnly_DW.FromWorkspace_PWORK.TimePtr;
        int_T currTimeIndex =
          ModelWithControllersOnly_DW.FromWorkspace_IWORK.PrevIndex;
        real_T t = (((ModelWithControllersOnly_M->Timing.clockTick1+
                      ModelWithControllersOnly_M->Timing.clockTickH1*
                      4294967296.0)) * 0.01);

        /* Get index */
        if (t <= pTimeValues[0]) {
          currTimeIndex = 0;
        } else if (t >= pTimeValues[7999]) {
          currTimeIndex = 7998;
        } else {
          if (t < pTimeValues[currTimeIndex]) {
            while (t < pTimeValues[currTimeIndex]) {
              currTimeIndex--;
            }
          } else {
            while (t >= pTimeValues[currTimeIndex + 1]) {
              currTimeIndex++;
            }
          }
        }

        ModelWithControllersOnly_DW.FromWorkspace_IWORK.PrevIndex =
          currTimeIndex;

        /* Post output */
        {
          real_T t1 = pTimeValues[currTimeIndex];
          real_T t2 = pTimeValues[currTimeIndex + 1];
          if (t1 == t2) {
            if (t < t1) {
              ModelWithControllersOnly_B.brake_cmdNm = pDataValues[currTimeIndex];
            } else {
              ModelWithControllersOnly_B.brake_cmdNm = pDataValues[currTimeIndex
                + 1];
            }
          } else {
            real_T f1 = (t2 - t) / (t2 - t1);
            real_T f2 = 1.0 - f1;
            real_T d1;
            real_T d2;
            int_T TimeIndex = currTimeIndex;
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            ModelWithControllersOnly_B.brake_cmdNm = (real_T) rtInterpolate(d1,
              d2, f1, f2);
            pDataValues += 8000;
          }
        }
      }

      /* FromWorkspace: '<S5>/From Workspace1' */
      {
        real_T *pDataValues = (real_T *)
          ModelWithControllersOnly_DW.FromWorkspace1_PWORK.DataPtr;
        real_T *pTimeValues = (real_T *)
          ModelWithControllersOnly_DW.FromWorkspace1_PWORK.TimePtr;
        int_T currTimeIndex =
          ModelWithControllersOnly_DW.FromWorkspace1_IWORK.PrevIndex;
        real_T t = (((ModelWithControllersOnly_M->Timing.clockTick1+
                      ModelWithControllersOnly_M->Timing.clockTickH1*
                      4294967296.0)) * 0.01);

        /* Get index */
        if (t <= pTimeValues[0]) {
          currTimeIndex = 0;
        } else if (t >= pTimeValues[7999]) {
          currTimeIndex = 7998;
        } else {
          if (t < pTimeValues[currTimeIndex]) {
            while (t < pTimeValues[currTimeIndex]) {
              currTimeIndex--;
            }
          } else {
            while (t >= pTimeValues[currTimeIndex + 1]) {
              currTimeIndex++;
            }
          }
        }

        ModelWithControllersOnly_DW.FromWorkspace1_IWORK.PrevIndex =
          currTimeIndex;

        /* Post output */
        {
          real_T t1 = pTimeValues[currTimeIndex];
          real_T t2 = pTimeValues[currTimeIndex + 1];
          if (t1 == t2) {
            if (t < t1) {
              ModelWithControllersOnly_B.steer_cmddeg =
                pDataValues[currTimeIndex];
            } else {
              ModelWithControllersOnly_B.steer_cmddeg =
                pDataValues[currTimeIndex + 1];
            }
          } else {
            real_T f1 = (t2 - t) / (t2 - t1);
            real_T f2 = 1.0 - f1;
            real_T d1;
            real_T d2;
            int_T TimeIndex = currTimeIndex;
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            ModelWithControllersOnly_B.steer_cmddeg = (real_T) rtInterpolate(d1,
              d2, f1, f2);
            pDataValues += 8000;
          }
        }
      }

      /* FromWorkspace: '<S5>/From Workspace4' */
      {
        real_T *pDataValues = (real_T *)
          ModelWithControllersOnly_DW.FromWorkspace4_PWORK.DataPtr;
        real_T *pTimeValues = (real_T *)
          ModelWithControllersOnly_DW.FromWorkspace4_PWORK.TimePtr;
        int_T currTimeIndex =
          ModelWithControllersOnly_DW.FromWorkspace4_IWORK.PrevIndex;
        real_T t = (((ModelWithControllersOnly_M->Timing.clockTick1+
                      ModelWithControllersOnly_M->Timing.clockTickH1*
                      4294967296.0)) * 0.01);

        /* Get index */
        if (t <= pTimeValues[0]) {
          currTimeIndex = 0;
        } else if (t >= pTimeValues[7999]) {
          currTimeIndex = 7998;
        } else {
          if (t < pTimeValues[currTimeIndex]) {
            while (t < pTimeValues[currTimeIndex]) {
              currTimeIndex--;
            }
          } else {
            while (t >= pTimeValues[currTimeIndex + 1]) {
              currTimeIndex++;
            }
          }
        }

        ModelWithControllersOnly_DW.FromWorkspace4_IWORK.PrevIndex =
          currTimeIndex;

        /* Post output */
        {
          real_T t1 = pTimeValues[currTimeIndex];
          real_T t2 = pTimeValues[currTimeIndex + 1];
          if (t1 == t2) {
            if (t < t1) {
              ModelWithControllersOnly_B.torque_cmdNm =
                pDataValues[currTimeIndex];
            } else {
              ModelWithControllersOnly_B.torque_cmdNm =
                pDataValues[currTimeIndex + 1];
            }
          } else {
            real_T f1 = (t2 - t) / (t2 - t1);
            real_T f2 = 1.0 - f1;
            real_T d1;
            real_T d2;
            int_T TimeIndex = currTimeIndex;
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            ModelWithControllersOnly_B.torque_cmdNm = (real_T) rtInterpolate(d1,
              d2, f1, f2);
            pDataValues += 8000;
          }
        }
      }

      /* Constant: '<S5>/Constant' */
      ModelWithControllersOnly_B.steer_cmddeg_n =
        ModelWithControllersOnly_cal->Constant_Value_h;

      /* Constant: '<S5>/Constant1' */
      ModelWithControllersOnly_B.torque_cmdNm_j =
        ModelWithControllersOnly_cal->Constant1_Value_b;

      /* Constant: '<S5>/Constant2' */
      ModelWithControllersOnly_B.brake_cmdNm_j =
        ModelWithControllersOnly_cal->Constant2_Value_o;

      /* Constant: '<S5>/Constant5' */
      ModelWithControllersOnly_B.gear_cmd_d =
        ModelWithControllersOnly_cal->Constant5_Value_h;
    }
  }

  if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
    real_T *lastU;

    /* Update for Integrator: '<S214>/Integrator' */
    ModelWithControllersOnly_DW.Integrator_IWORK = 0;

    /* Update for Integrator: '<S149>/Integrator' */
    ModelWithControllersOnly_DW.Integrator_IWORK_c = 0;
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      /* Update for UnitDelay: '<S8>/Unit Delay' */
      ModelWithControllersOnly_DW.UnitDelay_DSTATE =
        ModelWithControllersOnly_B.UnitDelay_l;

      /* Update for UnitDelay: '<Root>/Unit Delay' */
      ModelWithControllersOnly_DW.UnitDelay_DSTATE_e =
        ModelWithControllersOnly_B.Saturation_k;
    }

    /* Update for Derivative: '<S120>/Derivative' */
    if (ModelWithControllersOnly_DW.TimeStampA == (rtInf)) {
      ModelWithControllersOnly_DW.TimeStampA =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA;
    } else if (ModelWithControllersOnly_DW.TimeStampB == (rtInf)) {
      ModelWithControllersOnly_DW.TimeStampB =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeB;
    } else if (ModelWithControllersOnly_DW.TimeStampA <
               ModelWithControllersOnly_DW.TimeStampB) {
      ModelWithControllersOnly_DW.TimeStampA =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA;
    } else {
      ModelWithControllersOnly_DW.TimeStampB =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeB;
    }

    *lastU = ModelWithControllersOnly_B.y[0];

    /* End of Update for Derivative: '<S120>/Derivative' */

    /* Update for Derivative: '<S120>/Derivative1' */
    if (ModelWithControllersOnly_DW.TimeStampA_l == (rtInf)) {
      ModelWithControllersOnly_DW.TimeStampA_l =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA_d;
    } else if (ModelWithControllersOnly_DW.TimeStampB_l == (rtInf)) {
      ModelWithControllersOnly_DW.TimeStampB_l =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_o;
    } else if (ModelWithControllersOnly_DW.TimeStampA_l <
               ModelWithControllersOnly_DW.TimeStampB_l) {
      ModelWithControllersOnly_DW.TimeStampA_l =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA_d;
    } else {
      ModelWithControllersOnly_DW.TimeStampB_l =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_o;
    }

    *lastU = ModelWithControllersOnly_B.Gain5;

    /* End of Update for Derivative: '<S120>/Derivative1' */

    /* Update for Derivative: '<S120>/Derivative2' */
    if (ModelWithControllersOnly_DW.TimeStampA_e == (rtInf)) {
      ModelWithControllersOnly_DW.TimeStampA_e =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA_n;
    } else if (ModelWithControllersOnly_DW.TimeStampB_k == (rtInf)) {
      ModelWithControllersOnly_DW.TimeStampB_k =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_d;
    } else if (ModelWithControllersOnly_DW.TimeStampA_e <
               ModelWithControllersOnly_DW.TimeStampB_k) {
      ModelWithControllersOnly_DW.TimeStampA_e =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA_n;
    } else {
      ModelWithControllersOnly_DW.TimeStampB_k =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_d;
    }

    *lastU = ModelWithControllersOnly_B.Gain6;

    /* End of Update for Derivative: '<S120>/Derivative2' */

    /* Update for Derivative: '<S120>/Derivative3' */
    if (ModelWithControllersOnly_DW.TimeStampA_o == (rtInf)) {
      ModelWithControllersOnly_DW.TimeStampA_o =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA_nv;
    } else if (ModelWithControllersOnly_DW.TimeStampB_lf == (rtInf)) {
      ModelWithControllersOnly_DW.TimeStampB_lf =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_g;
    } else if (ModelWithControllersOnly_DW.TimeStampA_o <
               ModelWithControllersOnly_DW.TimeStampB_lf) {
      ModelWithControllersOnly_DW.TimeStampA_o =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA_nv;
    } else {
      ModelWithControllersOnly_DW.TimeStampB_lf =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_g;
    }

    *lastU = ModelWithControllersOnly_B.Roll;

    /* End of Update for Derivative: '<S120>/Derivative3' */

    /* Update for Derivative: '<S120>/Derivative4' */
    if (ModelWithControllersOnly_DW.TimeStampA_l5 == (rtInf)) {
      ModelWithControllersOnly_DW.TimeStampA_l5 =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA_k;
    } else if (ModelWithControllersOnly_DW.TimeStampB_e == (rtInf)) {
      ModelWithControllersOnly_DW.TimeStampB_e =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_p;
    } else if (ModelWithControllersOnly_DW.TimeStampA_l5 <
               ModelWithControllersOnly_DW.TimeStampB_e) {
      ModelWithControllersOnly_DW.TimeStampA_l5 =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA_k;
    } else {
      ModelWithControllersOnly_DW.TimeStampB_e =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_p;
    }

    *lastU = ModelWithControllersOnly_B.Pitch;

    /* End of Update for Derivative: '<S120>/Derivative4' */

    /* Update for Derivative: '<S120>/Derivative5' */
    if (ModelWithControllersOnly_DW.TimeStampA_k == (rtInf)) {
      ModelWithControllersOnly_DW.TimeStampA_k =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA_h;
    } else if (ModelWithControllersOnly_DW.TimeStampB_l4 == (rtInf)) {
      ModelWithControllersOnly_DW.TimeStampB_l4 =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_f;
    } else if (ModelWithControllersOnly_DW.TimeStampA_k <
               ModelWithControllersOnly_DW.TimeStampB_l4) {
      ModelWithControllersOnly_DW.TimeStampA_k =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA_h;
    } else {
      ModelWithControllersOnly_DW.TimeStampB_l4 =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_f;
    }

    *lastU = ModelWithControllersOnly_B.Yawdeg;

    /* End of Update for Derivative: '<S120>/Derivative5' */

    /* Update for Derivative: '<S121>/Derivative' */
    if (ModelWithControllersOnly_DW.TimeStampA_h == (rtInf)) {
      ModelWithControllersOnly_DW.TimeStampA_h =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA_o;
    } else if (ModelWithControllersOnly_DW.TimeStampB_h == (rtInf)) {
      ModelWithControllersOnly_DW.TimeStampB_h =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_c;
    } else if (ModelWithControllersOnly_DW.TimeStampA_h <
               ModelWithControllersOnly_DW.TimeStampB_h) {
      ModelWithControllersOnly_DW.TimeStampA_h =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeA_o;
    } else {
      ModelWithControllersOnly_DW.TimeStampB_h =
        ModelWithControllersOnly_M->Timing.t[0];
      lastU = &ModelWithControllersOnly_DW.LastUAtTimeB_c;
    }

    *lastU = ModelWithControllersOnly_B.Integrator[0];

    /* End of Update for Derivative: '<S121>/Derivative' */

    /* Update for FirstOrderHold: '<S121>/First Order Hold' */
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      real_T err;
      real_T tol;
      boolean_T guard1 = false;
      boolean_T resetCoeff;
      resetCoeff = (ModelWithControllersOnly_DW.Tk == (rtInf));
      guard1 = false;
      if (!resetCoeff) {
        if ((ModelWithControllersOnly_B.Gain2_o >= -1.0) &&
            (ModelWithControllersOnly_B.Gain2_o <= 1.0)) {
          tol = ModelWithControllersOnly_P.FirstOrderHold_ErrTol;
        } else if (ModelWithControllersOnly_B.Gain2_o > 1.0) {
          tol = ModelWithControllersOnly_B.Gain2_o *
            ModelWithControllersOnly_P.FirstOrderHold_ErrTol;
        } else {
          tol = -(ModelWithControllersOnly_B.Gain2_o *
                  ModelWithControllersOnly_P.FirstOrderHold_ErrTol);
        }

        err = ModelWithControllersOnly_B.FirstOrderHold -
          ModelWithControllersOnly_B.Gain2_o;
        if ((err > tol) || (err < -tol)) {
          guard1 = true;
        } else {
          tol = ModelWithControllersOnly_M->Timing.t[0] -
            ModelWithControllersOnly_DW.Tk;
          ModelWithControllersOnly_DW.Mk = (ModelWithControllersOnly_B.Gain2_o -
            ModelWithControllersOnly_DW.Uk) / tol;
          ModelWithControllersOnly_DW.Ck =
            ModelWithControllersOnly_B.FirstOrderHold;
        }
      } else {
        guard1 = true;
      }

      if (guard1) {
        if (ModelWithControllersOnly_B.Gain2_o !=
            ModelWithControllersOnly_B.FirstOrderHold) {
          rtsiSetBlockStateForSolverChangedAtMajorStep
            (&ModelWithControllersOnly_M->solverInfo, true);
          rtsiSetContTimeOutputInconsistentWithStateAtMajorStep
            (&ModelWithControllersOnly_M->solverInfo, true);
        }

        ModelWithControllersOnly_DW.Ck = ModelWithControllersOnly_B.Gain2_o;
        ModelWithControllersOnly_DW.Mk = 0.0;
      }

      ModelWithControllersOnly_DW.Uk = ModelWithControllersOnly_B.Gain2_o;
      ModelWithControllersOnly_DW.Tk = ModelWithControllersOnly_M->Timing.t[0];

      /* Update for FirstOrderHold: '<S121>/First Order Hold1' incorporates:
       *  FirstOrderHold: '<S121>/First Order Hold'
       */
      resetCoeff = (ModelWithControllersOnly_DW.Tk_d == (rtInf));
      guard1 = false;
      if (!resetCoeff) {
        if ((ModelWithControllersOnly_B.torquetoaccelerationconversion1 >= -1.0)
            && (ModelWithControllersOnly_B.torquetoaccelerationconversion1 <=
                1.0)) {
          tol = ModelWithControllersOnly_P.FirstOrderHold1_ErrTol;
        } else if (ModelWithControllersOnly_B.torquetoaccelerationconversion1 >
                   1.0) {
          tol = ModelWithControllersOnly_B.torquetoaccelerationconversion1 *
            ModelWithControllersOnly_P.FirstOrderHold1_ErrTol;
        } else {
          tol = -(ModelWithControllersOnly_B.torquetoaccelerationconversion1 *
                  ModelWithControllersOnly_P.FirstOrderHold1_ErrTol);
        }

        err = ModelWithControllersOnly_B.TnetNm -
          ModelWithControllersOnly_B.torquetoaccelerationconversion1;
        if ((err > tol) || (err < -tol)) {
          guard1 = true;
        } else {
          tol = ModelWithControllersOnly_M->Timing.t[0] -
            ModelWithControllersOnly_DW.Tk_d;
          ModelWithControllersOnly_DW.Mk_n =
            (ModelWithControllersOnly_B.torquetoaccelerationconversion1 -
             ModelWithControllersOnly_DW.Uk_d) / tol;
          ModelWithControllersOnly_DW.Ck_o = ModelWithControllersOnly_B.TnetNm;
        }
      } else {
        guard1 = true;
      }

      if (guard1) {
        if (ModelWithControllersOnly_B.torquetoaccelerationconversion1 !=
            ModelWithControllersOnly_B.TnetNm) {
          rtsiSetBlockStateForSolverChangedAtMajorStep
            (&ModelWithControllersOnly_M->solverInfo, true);
          rtsiSetContTimeOutputInconsistentWithStateAtMajorStep
            (&ModelWithControllersOnly_M->solverInfo, true);
        }

        ModelWithControllersOnly_DW.Ck_o =
          ModelWithControllersOnly_B.torquetoaccelerationconversion1;
        ModelWithControllersOnly_DW.Mk_n = 0.0;
      }

      ModelWithControllersOnly_DW.Uk_d =
        ModelWithControllersOnly_B.torquetoaccelerationconversion1;
      ModelWithControllersOnly_DW.Tk_d = ModelWithControllersOnly_M->Timing.t[0];

      /* End of Update for FirstOrderHold: '<S121>/First Order Hold1' */

      /* Update for DiscreteTransferFcn: '<S127>/second-order low-pass torque model' */
      ModelWithControllersOnly_DW.secondorderlowpasstorquemodel_s[1] =
        ModelWithControllersOnly_DW.secondorderlowpasstorquemodel_s[0];
      ModelWithControllersOnly_DW.secondorderlowpasstorquemodel_s[0] =
        ModelWithControllersOnly_DW.secondorderlowpasstorquemodel_t;

      /* Update for DiscreteTransferFcn: '<S127>/second-order low-pass brake model' */
      ModelWithControllersOnly_DW.secondorderlowpassbrakemodel_st[1] =
        ModelWithControllersOnly_DW.secondorderlowpassbrakemodel_st[0];
      ModelWithControllersOnly_DW.secondorderlowpassbrakemodel_st[0] =
        ModelWithControllersOnly_DW.secondorderlowpassbrakemodel_tm;

      /* Update for UnitDelay: '<S127>/Unit Delay' */
      ModelWithControllersOnly_DW.UnitDelay_DSTATE_l =
        ModelWithControllersOnly_B.Integrator[0];

      /* Update for DiscreteTransferFcn: '<S127>/second-order low-pass steering model' */
      ModelWithControllersOnly_DW.secondorderlowpasssteeringmodel[1] =
        ModelWithControllersOnly_DW.secondorderlowpasssteeringmodel[0];
      ModelWithControllersOnly_DW.secondorderlowpasssteeringmodel[0] =
        ModelWithControllersOnly_DW.secondorderlowpasssteeringmod_c;
    }

    /* End of Update for FirstOrderHold: '<S121>/First Order Hold' */

    /* ContTimeOutputInconsistentWithStateAtMajorOutputFlag is set, need to run a minor output */
    if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
      if (rtsiGetContTimeOutputInconsistentWithStateAtMajorStep
          (&ModelWithControllersOnly_M->solverInfo)) {
        rtsiSetSimTimeStep(&ModelWithControllersOnly_M->solverInfo,
                           MINOR_TIME_STEP);
        rtsiSetContTimeOutputInconsistentWithStateAtMajorStep
          (&ModelWithControllersOnly_M->solverInfo, false);
        ModelWithControllersOnly_step0();
        rtsiSetSimTimeStep(&ModelWithControllersOnly_M->solverInfo,
                           MAJOR_TIME_STEP);
      }
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
    rt_ertODEUpdateContinuousStates(&ModelWithControllersOnly_M->solverInfo);

    /* Update absolute time */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++ModelWithControllersOnly_M->Timing.clockTick0)) {
      ++ModelWithControllersOnly_M->Timing.clockTickH0;
    }

    ModelWithControllersOnly_M->Timing.t[0] = rtsiGetSolverStopTime
      (&ModelWithControllersOnly_M->solverInfo);

    /* Update absolute time */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The resolution of this integer timer is 0.01, which is the step size
     * of the task. Size of "clockTick1" ensures timer will not overflow during the
     * application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    ModelWithControllersOnly_M->Timing.clockTick1++;
    if (!ModelWithControllersOnly_M->Timing.clockTick1) {
      ModelWithControllersOnly_M->Timing.clockTickH1++;
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void ModelWithControllersOnly_derivatives(void)
{
  XDot_ModelWithControllersOnly_T *_rtXdot;
  _rtXdot = ((XDot_ModelWithControllersOnly_T *)
             ModelWithControllersOnly_M->derivs);

  /* Derivatives for Integrator: '<S214>/Integrator' */
  _rtXdot->Integrator_CSTATE[0] = ModelWithControllersOnly_B.stateDer[0];
  _rtXdot->Integrator_CSTATE[1] = ModelWithControllersOnly_B.stateDer[1];
  _rtXdot->Integrator_CSTATE[2] = ModelWithControllersOnly_B.stateDer[2];
  _rtXdot->Integrator_CSTATE[3] = ModelWithControllersOnly_B.stateDer[3];

  /* Derivatives for Integrator: '<S149>/Integrator' */
  _rtXdot->Integrator_CSTATE_b[0] = ModelWithControllersOnly_B.y[0];
  _rtXdot->Integrator_CSTATE_b[1] = ModelWithControllersOnly_B.y[1];

  /* Derivatives for Integrator: '<S119>/Integrator' */
  _rtXdot->Integrator_CSTATE_o =
    ModelWithControllersOnly_B.TmpRTBAtIntegratorInport1_c;

  /* Derivatives for Integrator: '<S118>/Integrator' */
  _rtXdot->Integrator_CSTATE_c =
    ModelWithControllersOnly_B.TmpRTBAtIntegratorInport1;

  /* Derivatives for Integrator: '<S212>/lateral' */
  _rtXdot->lateral_CSTATE = ModelWithControllersOnly_B.Product1_h;

  /* Derivatives for Integrator: '<S213>/lateral' */
  _rtXdot->lateral_CSTATE_k = ModelWithControllersOnly_B.Product1_af;
}

/* Model step function for TID2 */
void ModelWithControllersOnly_step2(void) /* Sample time: [0.05s, 0.0s] */
{
  real_T tmp;
  real_T u;
  int8_T rtAction;
  tmp = *get_MCAR();

  /* Reset subsysRan breadcrumbs */
  srClearBC(ModelWithControllersOnly_DW.Forward_SubsysRanBC);

  /* Reset subsysRan breadcrumbs */
  srClearBC(ModelWithControllersOnly_DW.Reverse_SubsysRanBC);

  /* RateTransition generated from: '<S12>/Equal2' */
  rtw_slrealtime_mutex_lock
    (ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_d0_SEMAPH);
  ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_RDBuf =
    ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_LstBufWR;
  rtw_slrealtime_mutex_unlock
    (ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_d0_SEMAPH);
  switch (ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_RDBuf) {
   case 0:
    /* RateTransition generated from: '<S12>/Equal2' */
    ModelWithControllersOnly_B.TmpRTBAtEqual2Inport1 =
      ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_Buf0;
    break;

   case 1:
    /* RateTransition generated from: '<S12>/Equal2' */
    ModelWithControllersOnly_B.TmpRTBAtEqual2Inport1 =
      ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_Buf1;
    break;

   case 2:
    /* RateTransition generated from: '<S12>/Equal2' */
    ModelWithControllersOnly_B.TmpRTBAtEqual2Inport1 =
      ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_Buf2;
    break;
  }

  /* End of RateTransition generated from: '<S12>/Equal2' */

  /* Signum: '<S12>/Sign2' incorporates:
   *  Constant: '<Root>/Constant2'
   */
  u = ModelWithControllersOnly_cal->Constant2_Value_e;
  if (u < 0.0) {
    /* Signum: '<S12>/Sign2' */
    ModelWithControllersOnly_B.Sign2 = -1.0;
  } else if (u > 0.0) {
    /* Signum: '<S12>/Sign2' */
    ModelWithControllersOnly_B.Sign2 = 1.0;
  } else if (u == 0.0) {
    /* Signum: '<S12>/Sign2' */
    ModelWithControllersOnly_B.Sign2 = 0.0;
  } else {
    /* Signum: '<S12>/Sign2' */
    ModelWithControllersOnly_B.Sign2 = (rtNaN);
  }

  /* End of Signum: '<S12>/Sign2' */

  /* RelationalOperator: '<S12>/Equal2' */
  ModelWithControllersOnly_B.Equal2 =
    (ModelWithControllersOnly_B.TmpRTBAtEqual2Inport1 !=
     ModelWithControllersOnly_B.Sign2);

  /* RelationalOperator: '<S12>/Equal6' incorporates:
   *  Constant: '<Root>/Constant2'
   *  Constant: '<S12>/Constant4'
   */
  ModelWithControllersOnly_B.Equal6 =
    (ModelWithControllersOnly_cal->Constant2_Value_e !=
     ModelWithControllersOnly_cal->Constant4_Value_i);

  /* Logic: '<S12>/AND3' */
  ModelWithControllersOnly_B.AND3 = (ModelWithControllersOnly_B.Equal2 &&
    ModelWithControllersOnly_B.Equal6);

  /* RateTransition generated from: '<S12>/Equal1' */
  rtw_slrealtime_mutex_lock
    (ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_d0_SEMAPH);
  ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_RDBuf =
    ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_LstBufWR;
  rtw_slrealtime_mutex_unlock
    (ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_d0_SEMAPH);
  switch (ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_RDBuf) {
   case 0:
    /* RateTransition generated from: '<S12>/Equal1' */
    ModelWithControllersOnly_B.TmpRTBAtEqual1Inport1 =
      ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_Buf0;
    break;

   case 1:
    /* RateTransition generated from: '<S12>/Equal1' */
    ModelWithControllersOnly_B.TmpRTBAtEqual1Inport1 =
      ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_Buf1;
    break;

   case 2:
    /* RateTransition generated from: '<S12>/Equal1' */
    ModelWithControllersOnly_B.TmpRTBAtEqual1Inport1 =
      ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_Buf2;
    break;
  }

  /* End of RateTransition generated from: '<S12>/Equal1' */

  /* RateTransition generated from: '<S12>/Sign1' */
  ModelWithControllersOnly_DW.xdot1_RdBufIdx = static_cast<int8_T>
    (ModelWithControllersOnly_DW.xdot1_RdBufIdx == 0);

  /* RateTransition generated from: '<S12>/Sign1' */
  ModelWithControllersOnly_B.xdot1 =
    ModelWithControllersOnly_DW.xdot1_Buf[ModelWithControllersOnly_DW.xdot1_RdBufIdx];

  /* Signum: '<S12>/Sign1' */
  u = ModelWithControllersOnly_B.xdot1;
  if (u < 0.0) {
    /* Signum: '<S12>/Sign1' */
    ModelWithControllersOnly_B.Sign1 = -1.0;
  } else if (u > 0.0) {
    /* Signum: '<S12>/Sign1' */
    ModelWithControllersOnly_B.Sign1 = 1.0;
  } else if (u == 0.0) {
    /* Signum: '<S12>/Sign1' */
    ModelWithControllersOnly_B.Sign1 = 0.0;
  } else {
    /* Signum: '<S12>/Sign1' */
    ModelWithControllersOnly_B.Sign1 = (rtNaN);
  }

  /* End of Signum: '<S12>/Sign1' */

  /* RelationalOperator: '<S12>/Equal1' */
  ModelWithControllersOnly_B.Equal1 =
    (ModelWithControllersOnly_B.TmpRTBAtEqual1Inport1 !=
     ModelWithControllersOnly_B.Sign1);

  /* RateTransition generated from: '<S12>/Equal3' */
  ModelWithControllersOnly_DW.xdot_RdBufIdx = static_cast<int8_T>
    (ModelWithControllersOnly_DW.xdot_RdBufIdx == 0);

  /* RateTransition generated from: '<S12>/Equal3' */
  ModelWithControllersOnly_B.xdot =
    ModelWithControllersOnly_DW.xdot_Buf[ModelWithControllersOnly_DW.xdot_RdBufIdx];

  /* RelationalOperator: '<S12>/Equal3' incorporates:
   *  Constant: '<S12>/Constant1'
   */
  ModelWithControllersOnly_B.Equal3 = (ModelWithControllersOnly_B.xdot !=
    ModelWithControllersOnly_cal->Constant1_Value_ks);

  /* Logic: '<S12>/AND1' */
  ModelWithControllersOnly_B.AND1 = (ModelWithControllersOnly_B.Equal1 &&
    ModelWithControllersOnly_B.Equal3);

  /* Logic: '<S12>/OR' */
  ModelWithControllersOnly_B.OR = (ModelWithControllersOnly_B.AND3 ||
    ModelWithControllersOnly_B.AND1);

  /* Logic: '<S12>/NOT1' */
  ModelWithControllersOnly_B.NOT1 = !ModelWithControllersOnly_B.OR;

  /* Assertion: '<S12>/Assertion' */
  ModelWithControllersOnly_DW.Assertion_sltestCurrentResult =
    !ModelWithControllersOnly_B.NOT1;
  if (ModelWithControllersOnly_DW.Assertion_sltestFinalResult <
      ModelWithControllersOnly_DW.Assertion_sltestCurrentResult) {
    ModelWithControllersOnly_DW.Assertion_sltestFinalResult =
      ModelWithControllersOnly_DW.Assertion_sltestCurrentResult;
    ModelWithControllersOnly_DW.Assertion_sltestLastResultTime =
      ModelWithControllersOnly_M->Timing.t[0];
  }

  slTestLogAssessments(&ModelWithControllersOnly_DW.Assertion_sltestBlkInfo,
                       ModelWithControllersOnly_M->Timing.t[0]);
  utAssert(ModelWithControllersOnly_B.NOT1);

  /* End of Assertion: '<S12>/Assertion' */

  /* RateTransition generated from: '<S12>/Equal5' */
  rtw_slrealtime_mutex_lock
    (ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_d0_SEMAPH);
  ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_RDBuf =
    ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_LstBufWR;
  rtw_slrealtime_mutex_unlock
    (ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_d0_SEMAPH);
  switch (ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_RDBuf) {
   case 0:
    /* RateTransition generated from: '<S12>/Equal5' */
    ModelWithControllersOnly_B.TmpRTBAtEqual5Inport2 =
      ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_Buf0;
    break;

   case 1:
    /* RateTransition generated from: '<S12>/Equal5' */
    ModelWithControllersOnly_B.TmpRTBAtEqual5Inport2 =
      ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_Buf1;
    break;

   case 2:
    /* RateTransition generated from: '<S12>/Equal5' */
    ModelWithControllersOnly_B.TmpRTBAtEqual5Inport2 =
      ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_Buf2;
    break;
  }

  /* End of RateTransition generated from: '<S12>/Equal5' */

  /* RelationalOperator: '<S12>/Equal5' incorporates:
   *  Constant: '<S12>/Constant3'
   */
  ModelWithControllersOnly_B.Equal5 =
    (ModelWithControllersOnly_cal->Constant3_Value_n !=
     ModelWithControllersOnly_B.TmpRTBAtEqual5Inport2);

  /* RateTransition generated from: '<S12>/Equal4' */
  rtw_slrealtime_mutex_lock
    (ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_d0_SEMAPH);
  ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_RDBuf =
    ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_LstBufWR;
  rtw_slrealtime_mutex_unlock
    (ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_d0_SEMAPH);
  switch (ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_RDBuf) {
   case 0:
    /* RateTransition generated from: '<S12>/Equal4' */
    ModelWithControllersOnly_B.TmpRTBAtEqual4Inport2 =
      ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_Buf0;
    break;

   case 1:
    /* RateTransition generated from: '<S12>/Equal4' */
    ModelWithControllersOnly_B.TmpRTBAtEqual4Inport2 =
      ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_Buf1;
    break;

   case 2:
    /* RateTransition generated from: '<S12>/Equal4' */
    ModelWithControllersOnly_B.TmpRTBAtEqual4Inport2 =
      ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_Buf2;
    break;
  }

  /* End of RateTransition generated from: '<S12>/Equal4' */

  /* RelationalOperator: '<S12>/Equal4' incorporates:
   *  Constant: '<S12>/Constant2'
   */
  ModelWithControllersOnly_B.Equal4 =
    (ModelWithControllersOnly_cal->Constant2_Value_i !=
     ModelWithControllersOnly_B.TmpRTBAtEqual4Inport2);

  /* Logic: '<S12>/AND2' */
  ModelWithControllersOnly_B.AND2 = (ModelWithControllersOnly_B.Equal5 &&
    ModelWithControllersOnly_B.Equal4);

  /* Logic: '<S12>/NOT' */
  ModelWithControllersOnly_B.NOT = !ModelWithControllersOnly_B.AND2;

  /* Assertion: '<S12>/Assertion1' */
  ModelWithControllersOnly_DW.Assertion1_sltestCurrentResult =
    !ModelWithControllersOnly_B.NOT;
  if (ModelWithControllersOnly_DW.Assertion1_sltestFinalResult <
      ModelWithControllersOnly_DW.Assertion1_sltestCurrentResult) {
    ModelWithControllersOnly_DW.Assertion1_sltestFinalResult =
      ModelWithControllersOnly_DW.Assertion1_sltestCurrentResult;
    ModelWithControllersOnly_DW.Assertion1_sltestLastResultTime =
      ModelWithControllersOnly_M->Timing.t[0];
  }

  slTestLogAssessments(&ModelWithControllersOnly_DW.Assertion1_sltestBlkInfo,
                       ModelWithControllersOnly_M->Timing.t[0]);
  utAssert(ModelWithControllersOnly_B.NOT);

  /* End of Assertion: '<S12>/Assertion1' */

  /* RateTransition generated from: '<S10>/Multiply' */
  rtw_slrealtime_mutex_lock
    (ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_d0_SEMA);
  ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_RDBuf =
    ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_LstBufW;
  rtw_slrealtime_mutex_unlock
    (ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_d0_SEMA);
  switch (ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_RDBuf) {
   case 0:
    /* RateTransition generated from: '<S10>/Multiply' */
    ModelWithControllersOnly_B.TmpRTBAtMultiplyInport1 =
      ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_Buf0;
    break;

   case 1:
    /* RateTransition generated from: '<S10>/Multiply' */
    ModelWithControllersOnly_B.TmpRTBAtMultiplyInport1 =
      ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_Buf1;
    break;

   case 2:
    /* RateTransition generated from: '<S10>/Multiply' */
    ModelWithControllersOnly_B.TmpRTBAtMultiplyInport1 =
      ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_Buf2;
    break;
  }

  /* End of RateTransition generated from: '<S10>/Multiply' */

  /* RateTransition generated from: '<S11>/Switch Case' */
  rtw_slrealtime_mutex_lock
    (ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_d0_SE);
  ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_RDBuf =
    ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_LstBu;
  rtw_slrealtime_mutex_unlock
    (ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_d0_SE);
  switch (ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_RDBuf) {
   case 0:
    /* RateTransition generated from: '<S11>/Switch Case' */
    ModelWithControllersOnly_B.TmpRTBAtSwitchCaseInport1 =
      ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_Buf0;
    break;

   case 1:
    /* RateTransition generated from: '<S11>/Switch Case' */
    ModelWithControllersOnly_B.TmpRTBAtSwitchCaseInport1 =
      ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_Buf1;
    break;

   case 2:
    /* RateTransition generated from: '<S11>/Switch Case' */
    ModelWithControllersOnly_B.TmpRTBAtSwitchCaseInport1 =
      ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_Buf2;
    break;
  }

  /* End of RateTransition generated from: '<S11>/Switch Case' */

  /* RateTransition generated from: '<S2>/Minus' */
  ModelWithControllersOnly_DW.xdot_RdBufIdx_m = static_cast<int8_T>
    (ModelWithControllersOnly_DW.xdot_RdBufIdx_m == 0);

  /* RateTransition generated from: '<S2>/Minus' */
  ModelWithControllersOnly_B.xdot_p =
    ModelWithControllersOnly_DW.xdot_Buf_d[ModelWithControllersOnly_DW.xdot_RdBufIdx_m];

  /* Sum: '<S2>/Minus' incorporates:
   *  Constant: '<Root>/Constant2'
   */
  ModelWithControllersOnly_B.error =
    ModelWithControllersOnly_cal->Constant2_Value_e -
    ModelWithControllersOnly_B.xdot_p;

  /* DataTypeConversion: '<S11>/Data Type Conversion' incorporates:
   *  Constant: '<Root>/Constant3'
   */
  ModelWithControllersOnly_B.DataTypeConversion =
    (ModelWithControllersOnly_cal->Constant3_Value_j2 != 0.0);

  /* SwitchCase: '<S11>/Switch Case' */
  rtAction = -1;
  if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
    u = ModelWithControllersOnly_B.TmpRTBAtSwitchCaseInport1;
    if (u < 0.0) {
      u = std::ceil(u);
    } else {
      u = std::floor(u);
    }

    if (rtIsNaN(u) || rtIsInf(u)) {
      u = 0.0;
    } else {
      u = std::fmod(u, 4.294967296E+9);
    }

    switch (u < 0.0 ? -static_cast<int32_T>(static_cast<uint32_T>(-u)) :
            static_cast<int32_T>(static_cast<uint32_T>(u))) {
     case 1:
      rtAction = 0;
      break;

     case -1:
      rtAction = 1;
      break;
    }

    ModelWithControllersOnly_DW.SwitchCase_ActiveSubsystem = rtAction;
  } else {
    rtAction = ModelWithControllersOnly_DW.SwitchCase_ActiveSubsystem;
  }

  switch (rtAction) {
   case 0:
    {
      real_T u1;
      real_T u2;

      /* Outputs for IfAction SubSystem: '<S11>/Forward' incorporates:
       *  ActionPort: '<S13>/Action Port'
       */
      /* Gain: '<S54>/Proportional Gain' */
      ModelWithControllersOnly_B.ProportionalGain_o =
        ModelWithControllersOnly_cal->LongitudinalControllerStanley_g *
        ModelWithControllersOnly_B.error;

      /* DiscreteIntegrator: '<S49>/Integrator' */
      if (ModelWithControllersOnly_B.DataTypeConversion ||
          (ModelWithControllersOnly_DW.Integrator_PrevResetState_n != 0)) {
        ModelWithControllersOnly_DW.Integrator_DSTATE_o =
          ModelWithControllersOnly_cal->PIForward_InitialConditionForIn;
      }

      /* DiscreteIntegrator: '<S49>/Integrator' */
      ModelWithControllersOnly_B.Integrator_c =
        ModelWithControllersOnly_DW.Integrator_DSTATE_o;

      /* Sum: '<S58>/Sum' */
      ModelWithControllersOnly_B.Sum_gu =
        ModelWithControllersOnly_B.ProportionalGain_o +
        ModelWithControllersOnly_B.Integrator_c;

      /* Gain: '<S40>/ZeroGain' */
      ModelWithControllersOnly_B.ZeroGain_p =
        ModelWithControllersOnly_cal->ZeroGain_Gain *
        ModelWithControllersOnly_B.Sum_gu;

      /* DeadZone: '<S42>/DeadZone' */
      if (ModelWithControllersOnly_B.Sum_gu >
          ModelWithControllersOnly_cal->PIForward_UpperSaturationLimit) {
        /* DeadZone: '<S42>/DeadZone' */
        ModelWithControllersOnly_B.DeadZone_o =
          ModelWithControllersOnly_B.Sum_gu -
          ModelWithControllersOnly_cal->PIForward_UpperSaturationLimit;
      } else if (ModelWithControllersOnly_B.Sum_gu >=
                 ModelWithControllersOnly_cal->PIForward_LowerSaturationLimit) {
        /* DeadZone: '<S42>/DeadZone' */
        ModelWithControllersOnly_B.DeadZone_o = 0.0;
      } else {
        /* DeadZone: '<S42>/DeadZone' */
        ModelWithControllersOnly_B.DeadZone_o =
          ModelWithControllersOnly_B.Sum_gu -
          ModelWithControllersOnly_cal->PIForward_LowerSaturationLimit;
      }

      /* End of DeadZone: '<S42>/DeadZone' */

      /* RelationalOperator: '<S40>/NotEqual' */
      ModelWithControllersOnly_B.NotEqual_m =
        (ModelWithControllersOnly_B.ZeroGain_p !=
         ModelWithControllersOnly_B.DeadZone_o);

      /* Signum: '<S40>/SignPreSat' */
      u = ModelWithControllersOnly_B.DeadZone_o;
      if (u < 0.0) {
        /* Signum: '<S40>/SignPreSat' */
        ModelWithControllersOnly_B.SignPreSat_n = -1.0;
      } else if (u > 0.0) {
        /* Signum: '<S40>/SignPreSat' */
        ModelWithControllersOnly_B.SignPreSat_n = 1.0;
      } else if (u == 0.0) {
        /* Signum: '<S40>/SignPreSat' */
        ModelWithControllersOnly_B.SignPreSat_n = 0.0;
      } else {
        /* Signum: '<S40>/SignPreSat' */
        ModelWithControllersOnly_B.SignPreSat_n = (rtNaN);
      }

      /* End of Signum: '<S40>/SignPreSat' */

      /* DataTypeConversion: '<S40>/DataTypeConv1' */
      u = std::floor(ModelWithControllersOnly_B.SignPreSat_n);
      if (rtIsNaN(u) || rtIsInf(u)) {
        u = 0.0;
      } else {
        u = std::fmod(u, 256.0);
      }

      /* DataTypeConversion: '<S40>/DataTypeConv1' */
      ModelWithControllersOnly_B.DataTypeConv1_j = static_cast<int8_T>(u < 0.0 ?
        static_cast<int32_T>(static_cast<int8_T>(-static_cast<int8_T>(
        static_cast<uint8_T>(-u)))) : static_cast<int32_T>(static_cast<int8_T>(
        static_cast<uint8_T>(u))));

      /* Gain: '<S46>/Integral Gain' */
      ModelWithControllersOnly_B.IntegralGain_e =
        ModelWithControllersOnly_cal->LongitudinalControllerStanley_K *
        ModelWithControllersOnly_B.error;

      /* Signum: '<S40>/SignPreIntegrator' */
      u = ModelWithControllersOnly_B.IntegralGain_e;
      if (u < 0.0) {
        /* Signum: '<S40>/SignPreIntegrator' */
        ModelWithControllersOnly_B.SignPreIntegrator_k = -1.0;
      } else if (u > 0.0) {
        /* Signum: '<S40>/SignPreIntegrator' */
        ModelWithControllersOnly_B.SignPreIntegrator_k = 1.0;
      } else if (u == 0.0) {
        /* Signum: '<S40>/SignPreIntegrator' */
        ModelWithControllersOnly_B.SignPreIntegrator_k = 0.0;
      } else {
        /* Signum: '<S40>/SignPreIntegrator' */
        ModelWithControllersOnly_B.SignPreIntegrator_k = (rtNaN);
      }

      /* End of Signum: '<S40>/SignPreIntegrator' */

      /* DataTypeConversion: '<S40>/DataTypeConv2' */
      u = std::floor(ModelWithControllersOnly_B.SignPreIntegrator_k);
      if (rtIsNaN(u) || rtIsInf(u)) {
        u = 0.0;
      } else {
        u = std::fmod(u, 256.0);
      }

      /* DataTypeConversion: '<S40>/DataTypeConv2' */
      ModelWithControllersOnly_B.DataTypeConv2_g = static_cast<int8_T>(u < 0.0 ?
        static_cast<int32_T>(static_cast<int8_T>(-static_cast<int8_T>(
        static_cast<uint8_T>(-u)))) : static_cast<int32_T>(static_cast<int8_T>(
        static_cast<uint8_T>(u))));

      /* RelationalOperator: '<S40>/Equal1' */
      ModelWithControllersOnly_B.Equal1_j =
        (ModelWithControllersOnly_B.DataTypeConv1_j ==
         ModelWithControllersOnly_B.DataTypeConv2_g);

      /* Logic: '<S40>/AND3' */
      ModelWithControllersOnly_B.AND3_k = (ModelWithControllersOnly_B.NotEqual_m
        && ModelWithControllersOnly_B.Equal1_j);

      /* Switch: '<S40>/Switch' */
      if (ModelWithControllersOnly_B.AND3_k) {
        /* Switch: '<S40>/Switch' incorporates:
         *  Constant: '<S40>/Constant1'
         */
        ModelWithControllersOnly_B.Switch_a =
          ModelWithControllersOnly_cal->Constant1_Value;
      } else {
        /* Switch: '<S40>/Switch' */
        ModelWithControllersOnly_B.Switch_a =
          ModelWithControllersOnly_B.IntegralGain_e;
      }

      /* End of Switch: '<S40>/Switch' */

      /* Saturate: '<S56>/Saturation' */
      u = ModelWithControllersOnly_B.Sum_gu;
      u1 = ModelWithControllersOnly_cal->PIForward_LowerSaturationLimit;
      u2 = ModelWithControllersOnly_cal->PIForward_UpperSaturationLimit;
      if (u > u2) {
        /* Merge: '<S11>/Merge' */
        ModelWithControllersOnly_B.Merge = u2;
      } else if (u < u1) {
        /* Merge: '<S11>/Merge' */
        ModelWithControllersOnly_B.Merge = u1;
      } else {
        /* Merge: '<S11>/Merge' */
        ModelWithControllersOnly_B.Merge = u;
      }

      /* End of Saturate: '<S56>/Saturation' */

      /* Update for DiscreteIntegrator: '<S49>/Integrator' */
      ModelWithControllersOnly_DW.Integrator_DSTATE_o +=
        ModelWithControllersOnly_cal->Integrator_gainval *
        ModelWithControllersOnly_B.Switch_a;
      ModelWithControllersOnly_DW.Integrator_PrevResetState_n =
        static_cast<int8_T>(ModelWithControllersOnly_B.DataTypeConversion);
      if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
        srUpdateBC(ModelWithControllersOnly_DW.Forward_SubsysRanBC);
      }

      /* End of Outputs for SubSystem: '<S11>/Forward' */
    }
    break;

   case 1:
    {
      real_T u1;
      real_T u2;

      /* Outputs for IfAction SubSystem: '<S11>/Reverse' incorporates:
       *  ActionPort: '<S14>/Action Port'
       */
      /* Gain: '<S105>/Proportional Gain' */
      ModelWithControllersOnly_B.ProportionalGain =
        ModelWithControllersOnly_cal->LongitudinalControllerStanley_g *
        ModelWithControllersOnly_B.error;

      /* DiscreteIntegrator: '<S100>/Integrator' */
      if (ModelWithControllersOnly_B.DataTypeConversion ||
          (ModelWithControllersOnly_DW.Integrator_PrevResetState != 0)) {
        ModelWithControllersOnly_DW.Integrator_DSTATE =
          ModelWithControllersOnly_cal->PIReverse_InitialConditionForIn;
      }

      /* DiscreteIntegrator: '<S100>/Integrator' */
      ModelWithControllersOnly_B.Integrator_f =
        ModelWithControllersOnly_DW.Integrator_DSTATE;

      /* Sum: '<S109>/Sum' */
      ModelWithControllersOnly_B.Sum_mr =
        ModelWithControllersOnly_B.ProportionalGain +
        ModelWithControllersOnly_B.Integrator_f;

      /* Gain: '<S91>/ZeroGain' */
      ModelWithControllersOnly_B.ZeroGain =
        ModelWithControllersOnly_cal->ZeroGain_Gain_i *
        ModelWithControllersOnly_B.Sum_mr;

      /* DeadZone: '<S93>/DeadZone' */
      if (ModelWithControllersOnly_B.Sum_mr >
          ModelWithControllersOnly_cal->PIReverse_UpperSaturationLimit) {
        /* DeadZone: '<S93>/DeadZone' */
        ModelWithControllersOnly_B.DeadZone_f =
          ModelWithControllersOnly_B.Sum_mr -
          ModelWithControllersOnly_cal->PIReverse_UpperSaturationLimit;
      } else if (ModelWithControllersOnly_B.Sum_mr >=
                 ModelWithControllersOnly_cal->PIReverse_LowerSaturationLimit) {
        /* DeadZone: '<S93>/DeadZone' */
        ModelWithControllersOnly_B.DeadZone_f = 0.0;
      } else {
        /* DeadZone: '<S93>/DeadZone' */
        ModelWithControllersOnly_B.DeadZone_f =
          ModelWithControllersOnly_B.Sum_mr -
          ModelWithControllersOnly_cal->PIReverse_LowerSaturationLimit;
      }

      /* End of DeadZone: '<S93>/DeadZone' */

      /* RelationalOperator: '<S91>/NotEqual' */
      ModelWithControllersOnly_B.NotEqual = (ModelWithControllersOnly_B.ZeroGain
        != ModelWithControllersOnly_B.DeadZone_f);

      /* Signum: '<S91>/SignPreSat' */
      u = ModelWithControllersOnly_B.DeadZone_f;
      if (u < 0.0) {
        /* Signum: '<S91>/SignPreSat' */
        ModelWithControllersOnly_B.SignPreSat = -1.0;
      } else if (u > 0.0) {
        /* Signum: '<S91>/SignPreSat' */
        ModelWithControllersOnly_B.SignPreSat = 1.0;
      } else if (u == 0.0) {
        /* Signum: '<S91>/SignPreSat' */
        ModelWithControllersOnly_B.SignPreSat = 0.0;
      } else {
        /* Signum: '<S91>/SignPreSat' */
        ModelWithControllersOnly_B.SignPreSat = (rtNaN);
      }

      /* End of Signum: '<S91>/SignPreSat' */

      /* DataTypeConversion: '<S91>/DataTypeConv1' */
      u = std::floor(ModelWithControllersOnly_B.SignPreSat);
      if (rtIsNaN(u) || rtIsInf(u)) {
        u = 0.0;
      } else {
        u = std::fmod(u, 256.0);
      }

      /* DataTypeConversion: '<S91>/DataTypeConv1' */
      ModelWithControllersOnly_B.DataTypeConv1 = static_cast<int8_T>(u < 0.0 ?
        static_cast<int32_T>(static_cast<int8_T>(-static_cast<int8_T>(
        static_cast<uint8_T>(-u)))) : static_cast<int32_T>(static_cast<int8_T>(
        static_cast<uint8_T>(u))));

      /* Gain: '<S97>/Integral Gain' */
      ModelWithControllersOnly_B.IntegralGain =
        ModelWithControllersOnly_cal->LongitudinalControllerStanley_K *
        ModelWithControllersOnly_B.error;

      /* Signum: '<S91>/SignPreIntegrator' */
      u = ModelWithControllersOnly_B.IntegralGain;
      if (u < 0.0) {
        /* Signum: '<S91>/SignPreIntegrator' */
        ModelWithControllersOnly_B.SignPreIntegrator = -1.0;
      } else if (u > 0.0) {
        /* Signum: '<S91>/SignPreIntegrator' */
        ModelWithControllersOnly_B.SignPreIntegrator = 1.0;
      } else if (u == 0.0) {
        /* Signum: '<S91>/SignPreIntegrator' */
        ModelWithControllersOnly_B.SignPreIntegrator = 0.0;
      } else {
        /* Signum: '<S91>/SignPreIntegrator' */
        ModelWithControllersOnly_B.SignPreIntegrator = (rtNaN);
      }

      /* End of Signum: '<S91>/SignPreIntegrator' */

      /* DataTypeConversion: '<S91>/DataTypeConv2' */
      u = std::floor(ModelWithControllersOnly_B.SignPreIntegrator);
      if (rtIsNaN(u) || rtIsInf(u)) {
        u = 0.0;
      } else {
        u = std::fmod(u, 256.0);
      }

      /* DataTypeConversion: '<S91>/DataTypeConv2' */
      ModelWithControllersOnly_B.DataTypeConv2 = static_cast<int8_T>(u < 0.0 ?
        static_cast<int32_T>(static_cast<int8_T>(-static_cast<int8_T>(
        static_cast<uint8_T>(-u)))) : static_cast<int32_T>(static_cast<int8_T>(
        static_cast<uint8_T>(u))));

      /* RelationalOperator: '<S91>/Equal1' */
      ModelWithControllersOnly_B.Equal1_e =
        (ModelWithControllersOnly_B.DataTypeConv1 ==
         ModelWithControllersOnly_B.DataTypeConv2);

      /* Logic: '<S91>/AND3' */
      ModelWithControllersOnly_B.AND3_h = (ModelWithControllersOnly_B.NotEqual &&
        ModelWithControllersOnly_B.Equal1_e);

      /* Switch: '<S91>/Switch' */
      if (ModelWithControllersOnly_B.AND3_h) {
        /* Switch: '<S91>/Switch' incorporates:
         *  Constant: '<S91>/Constant1'
         */
        ModelWithControllersOnly_B.Switch_m =
          ModelWithControllersOnly_cal->Constant1_Value_f;
      } else {
        /* Switch: '<S91>/Switch' */
        ModelWithControllersOnly_B.Switch_m =
          ModelWithControllersOnly_B.IntegralGain;
      }

      /* End of Switch: '<S91>/Switch' */

      /* Saturate: '<S107>/Saturation' */
      u = ModelWithControllersOnly_B.Sum_mr;
      u1 = ModelWithControllersOnly_cal->PIReverse_LowerSaturationLimit;
      u2 = ModelWithControllersOnly_cal->PIReverse_UpperSaturationLimit;
      if (u > u2) {
        /* Merge: '<S11>/Merge' */
        ModelWithControllersOnly_B.Merge = u2;
      } else if (u < u1) {
        /* Merge: '<S11>/Merge' */
        ModelWithControllersOnly_B.Merge = u1;
      } else {
        /* Merge: '<S11>/Merge' */
        ModelWithControllersOnly_B.Merge = u;
      }

      /* End of Saturate: '<S107>/Saturation' */

      /* Update for DiscreteIntegrator: '<S100>/Integrator' */
      ModelWithControllersOnly_DW.Integrator_DSTATE +=
        ModelWithControllersOnly_cal->Integrator_gainval_e *
        ModelWithControllersOnly_B.Switch_m;
      ModelWithControllersOnly_DW.Integrator_PrevResetState = static_cast<int8_T>
        (ModelWithControllersOnly_B.DataTypeConversion);
      if (rtmIsMajorTimeStep(ModelWithControllersOnly_M)) {
        srUpdateBC(ModelWithControllersOnly_DW.Reverse_SubsysRanBC);
      }

      /* End of Outputs for SubSystem: '<S11>/Reverse' */
    }
    break;
  }

  /* End of SwitchCase: '<S11>/Switch Case' */

  /* Product: '<S10>/Multiply' */
  ModelWithControllersOnly_B.Multiply_n =
    ModelWithControllersOnly_B.TmpRTBAtMultiplyInport1 *
    ModelWithControllersOnly_B.Merge;

  /* Switch: '<S10>/Switch' */
  if (ModelWithControllersOnly_B.Multiply_n >
      ModelWithControllersOnly_cal->Switch_Threshold_l) {
    /* Switch: '<S10>/Switch' */
    ModelWithControllersOnly_B.Switch_l = ModelWithControllersOnly_B.Multiply_n;
  } else {
    /* Switch: '<S10>/Switch' incorporates:
     *  Constant: '<S10>/Constant'
     */
    ModelWithControllersOnly_B.Switch_l =
      ModelWithControllersOnly_cal->Constant_Value_c;
  }

  /* End of Switch: '<S10>/Switch' */

  /* Switch: '<S10>/Switch1' */
  if (ModelWithControllersOnly_B.Multiply_n >
      ModelWithControllersOnly_cal->Switch1_Threshold_i) {
    /* Switch: '<S10>/Switch1' incorporates:
     *  Constant: '<S10>/Constant'
     */
    ModelWithControllersOnly_B.Switch1_l =
      ModelWithControllersOnly_cal->Constant_Value_c;
  } else {
    /* Abs: '<S10>/Abs' */
    ModelWithControllersOnly_B.Abs_k = std::abs
      (ModelWithControllersOnly_B.Multiply_n);

    /* Switch: '<S10>/Switch1' */
    ModelWithControllersOnly_B.Switch1_l = ModelWithControllersOnly_B.Abs_k;
  }

  /* End of Switch: '<S10>/Switch1' */

  /* Gain: '<S118>/Gain' */
  u = 1.1 * tmp;

  /* Gain: '<S118>/Gain' */
  ModelWithControllersOnly_B.Gain_g = u * ModelWithControllersOnly_B.Switch_l;

  /* RateTransition generated from: '<S118>/Integrator' */
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_WrBuf =
    static_cast<int8_T>
    (ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_WrBuf == 0);
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_Buf[ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_WrBuf]
    = ModelWithControllersOnly_B.Switch_l;

  /* RateTransition generated from: '<S118>/Integrator' */
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_WrB_d =
    static_cast<int8_T>
    (ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_WrB_d == 0);
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_Buf_m[ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_WrB_d]
    = ModelWithControllersOnly_B.Switch_l;

  /* RateTransition generated from: '<S118>/Sum' */
  ModelWithControllersOnly_DW.TmpRTBAtSumInport3_WrBufIdx_o = static_cast<int8_T>
    (ModelWithControllersOnly_DW.TmpRTBAtSumInport3_WrBufIdx_o == 0);
  ModelWithControllersOnly_DW.TmpRTBAtSumInport3_Buf_f[ModelWithControllersOnly_DW.TmpRTBAtSumInport3_WrBufIdx_o]
    = ModelWithControllersOnly_B.Gain_g;

  /* Gain: '<S119>/Gain' */
  u = 1.1 * tmp;

  /* Gain: '<S119>/Gain' */
  ModelWithControllersOnly_B.Gain_ne = u * ModelWithControllersOnly_B.Switch1_l;

  /* RateTransition generated from: '<S119>/Integrator' */
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_WrB_k =
    static_cast<int8_T>
    (ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_WrB_k == 0);
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_Buf_p[ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_WrB_k]
    = ModelWithControllersOnly_B.Switch1_l;

  /* RateTransition generated from: '<S119>/Integrator' */
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_WrBuf =
    static_cast<int8_T>
    (ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_WrBuf == 0);
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_Buf[ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_WrBuf]
    = ModelWithControllersOnly_B.Switch1_l;

  /* RateTransition generated from: '<S119>/Sum' */
  ModelWithControllersOnly_DW.TmpRTBAtSumInport3_WrBufIdx = static_cast<int8_T>
    (ModelWithControllersOnly_DW.TmpRTBAtSumInport3_WrBufIdx == 0);
  ModelWithControllersOnly_DW.TmpRTBAtSumInport3_Buf[ModelWithControllersOnly_DW.TmpRTBAtSumInport3_WrBufIdx]
    = ModelWithControllersOnly_B.Gain_ne;

  /* RateTransition generated from: '<S120>/MATLAB Function3' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport_a =
    static_cast<int8_T>
    (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport_a == 0);

  /* RateTransition generated from: '<S120>/MATLAB Function3' */
  ModelWithControllersOnly_B.TmpRTBAtMATLABFunction3Inport1 =
    ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport1_[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport_a];

  /* RateTransition generated from: '<S120>/MATLAB Function3' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport_f =
    static_cast<int8_T>
    (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport_f == 0);

  /* RateTransition generated from: '<S120>/MATLAB Function3' */
  ModelWithControllersOnly_B.TmpRTBAtMATLABFunction3Inport2 =
    ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport2_[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport_f];

  /* RateTransition generated from: '<S120>/MATLAB Function3' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inpor_gh =
    static_cast<int8_T>
    (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inpor_gh == 0);

  /* RateTransition generated from: '<S120>/MATLAB Function3' */
  ModelWithControllersOnly_B.TmpRTBAtMATLABFunction3Inport3 =
    ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport3_[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inpor_gh];

  /* MATLAB Function: '<S120>/MATLAB Function3' */
  ModelWithControllersOnly_B.AngularVelocity[0] =
    ModelWithControllersOnly_B.TmpRTBAtMATLABFunction3Inport1;
  ModelWithControllersOnly_B.AngularVelocity[1] =
    ModelWithControllersOnly_B.TmpRTBAtMATLABFunction3Inport2;
  ModelWithControllersOnly_B.AngularVelocity[2] =
    ModelWithControllersOnly_B.TmpRTBAtMATLABFunction3Inport3;

  /* RateTransition generated from: '<S120>/MATLAB Function' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport1_R =
    static_cast<int8_T>
    (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport1_R == 0);

  /* RateTransition generated from: '<S120>/MATLAB Function' */
  ModelWithControllersOnly_B.X =
    ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport1_B[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport1_R];

  /* RateTransition generated from: '<S120>/MATLAB Function' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport2_R =
    static_cast<int8_T>
    (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport2_R == 0);

  /* RateTransition generated from: '<S120>/MATLAB Function' */
  ModelWithControllersOnly_B.TmpRTBAtMATLABFunctionInport2 =
    ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport2_B[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport2_R];

  /* RateTransition generated from: '<S120>/MATLAB Function' */
  rtw_slrealtime_mutex_lock
    (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_d);
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_R =
    ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_L;
  rtw_slrealtime_mutex_unlock
    (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_d);
  switch (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_R) {
   case 0:
    /* RateTransition generated from: '<S120>/MATLAB Function' */
    ModelWithControllersOnly_B.TmpRTBAtMATLABFunctionInport3 =
      ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_B;
    break;

   case 1:
    /* RateTransition generated from: '<S120>/MATLAB Function' */
    ModelWithControllersOnly_B.TmpRTBAtMATLABFunctionInport3 =
      ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_c;
    break;

   case 2:
    /* RateTransition generated from: '<S120>/MATLAB Function' */
    ModelWithControllersOnly_B.TmpRTBAtMATLABFunctionInport3 =
      ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_f;
    break;
  }

  /* End of RateTransition generated from: '<S120>/MATLAB Function' */

  /* MATLAB Function: '<S120>/MATLAB Function' */
  ModelWithControllersOnly_B.BusCreator.Position[0] =
    ModelWithControllersOnly_B.X;
  ModelWithControllersOnly_B.BusCreator.Position[1] =
    ModelWithControllersOnly_B.TmpRTBAtMATLABFunctionInport2;
  ModelWithControllersOnly_B.BusCreator.Position[2] =
    ModelWithControllersOnly_B.TmpRTBAtMATLABFunctionInport3;

  /* RateTransition generated from: '<S120>/MATLAB Function1' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_h =
    static_cast<int8_T>
    (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_h == 0);

  /* RateTransition generated from: '<S120>/MATLAB Function1' */
  ModelWithControllersOnly_B.Xdot =
    ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport1_[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_h];

  /* RateTransition generated from: '<S120>/MATLAB Function1' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_p =
    static_cast<int8_T>
    (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_p == 0);

  /* RateTransition generated from: '<S120>/MATLAB Function1' */
  ModelWithControllersOnly_B.TmpRTBAtMATLABFunction1Inport2 =
    ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport2_[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_p];

  /* RateTransition generated from: '<S120>/MATLAB Function1' */
  rtw_slrealtime_mutex_lock
    (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_k);
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_b =
    ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inpor_oh;
  rtw_slrealtime_mutex_unlock
    (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_k);
  switch (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_b) {
   case 0:
    /* RateTransition generated from: '<S120>/MATLAB Function1' */
    ModelWithControllersOnly_B.TmpRTBAtMATLABFunction1Inport3 =
      ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport3_;
    break;

   case 1:
    /* RateTransition generated from: '<S120>/MATLAB Function1' */
    ModelWithControllersOnly_B.TmpRTBAtMATLABFunction1Inport3 =
      ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_f;
    break;

   case 2:
    /* RateTransition generated from: '<S120>/MATLAB Function1' */
    ModelWithControllersOnly_B.TmpRTBAtMATLABFunction1Inport3 =
      ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_o;
    break;
  }

  /* End of RateTransition generated from: '<S120>/MATLAB Function1' */

  /* MATLAB Function: '<S120>/MATLAB Function1' */
  ModelWithControllersOnly_B.BusCreator.Velocity[0] =
    ModelWithControllersOnly_B.Xdot;
  ModelWithControllersOnly_B.BusCreator.Velocity[1] =
    ModelWithControllersOnly_B.TmpRTBAtMATLABFunction1Inport2;
  ModelWithControllersOnly_B.BusCreator.Velocity[2] =
    ModelWithControllersOnly_B.TmpRTBAtMATLABFunction1Inport3;

  /* RateTransition generated from: '<S120>/MATLAB Function2' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_e =
    static_cast<int8_T>
    (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_e == 0);

  /* RateTransition generated from: '<S120>/MATLAB Function2' */
  ModelWithControllersOnly_B.TmpRTBAtMATLABFunction2Inport1 =
    ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport1_[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_e];

  /* RateTransition generated from: '<S120>/MATLAB Function2' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_m =
    static_cast<int8_T>
    (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_m == 0);

  /* RateTransition generated from: '<S120>/MATLAB Function2' */
  ModelWithControllersOnly_B.TmpRTBAtMATLABFunction2Inport2 =
    ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport2_[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_m];

  /* RateTransition generated from: '<S120>/MATLAB Function2' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_a =
    static_cast<int8_T>
    (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_a == 0);

  /* RateTransition generated from: '<S120>/MATLAB Function2' */
  ModelWithControllersOnly_B.TmpRTBAtMATLABFunction2Inport3 =
    ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport3_[ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_a];

  /* MATLAB Function: '<S120>/MATLAB Function2' */
  ModelWithControllersOnly_B.BusCreator.Acceleration[0] =
    ModelWithControllersOnly_B.TmpRTBAtMATLABFunction2Inport1;
  ModelWithControllersOnly_B.BusCreator.Acceleration[1] =
    ModelWithControllersOnly_B.TmpRTBAtMATLABFunction2Inport2;
  ModelWithControllersOnly_B.BusCreator.Acceleration[2] =
    ModelWithControllersOnly_B.TmpRTBAtMATLABFunction2Inport3;

  /* RateTransition: '<S120>/Rate Transition' */
  ModelWithControllersOnly_DW.RateTransition_RdBufIdx = static_cast<int8_T>
    (ModelWithControllersOnly_DW.RateTransition_RdBufIdx == 0);
  ModelWithControllersOnly_B.BusCreator.Roll =
    ModelWithControllersOnly_DW.RateTransition_Buf[ModelWithControllersOnly_DW.RateTransition_RdBufIdx];

  /* RateTransition: '<S120>/Rate Transition1' */
  ModelWithControllersOnly_DW.RateTransition1_RdBufIdx = static_cast<int8_T>
    (ModelWithControllersOnly_DW.RateTransition1_RdBufIdx == 0);
  ModelWithControllersOnly_B.BusCreator.Pitch =
    ModelWithControllersOnly_DW.RateTransition1_Buf[ModelWithControllersOnly_DW.RateTransition1_RdBufIdx];

  /* RateTransition: '<S120>/Rate Transition2' */
  ModelWithControllersOnly_DW.RateTransition2_RdBufIdx = static_cast<int8_T>
    (ModelWithControllersOnly_DW.RateTransition2_RdBufIdx == 0);
  ModelWithControllersOnly_B.BusCreator.Yaw =
    ModelWithControllersOnly_DW.RateTransition2_Buf[ModelWithControllersOnly_DW.RateTransition2_RdBufIdx];

  /* BusCreator: '<S120>/Bus Creator' incorporates:
   *  Constant: '<S120>/Constant'
   */
  ModelWithControllersOnly_B.BusCreator.ActorID =
    ModelWithControllersOnly_cal->Constant_Value_ke;
  ModelWithControllersOnly_B.BusCreator.AngularVelocity[0] =
    ModelWithControllersOnly_B.AngularVelocity[0];
  ModelWithControllersOnly_B.BusCreator.AngularVelocity[1] =
    ModelWithControllersOnly_B.AngularVelocity[1];
  ModelWithControllersOnly_B.BusCreator.AngularVelocity[2] =
    ModelWithControllersOnly_B.AngularVelocity[2];

  /* Update absolute time */
  /* The "clockTick2" counts the number of times the code of this task has
   * been executed. The resolution of this integer timer is 0.05, which is the step size
   * of the task. Size of "clockTick2" ensures timer will not overflow during the
   * application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick2 and the high bits
   * Timing.clockTickH2. When the low bit overflows to 0, the high bits increment.
   */
  ModelWithControllersOnly_M->Timing.clockTick2++;
  if (!ModelWithControllersOnly_M->Timing.clockTick2) {
    ModelWithControllersOnly_M->Timing.clockTickH2++;
  }
}

/* Model initialize function */
void ModelWithControllersOnly_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  ModelWithControllersOnly_P.Saturation_UpperSat = rtInf;
  ModelWithControllersOnly_P.FirstOrderHold_ErrTol = rtInf;
  ModelWithControllersOnly_P.FirstOrderHold1_ErrTol = rtInf;
  (ModelWithControllersOnly_M)->Timing.TaskCounters.cLimit[0] = 1;
  (ModelWithControllersOnly_M)->Timing.TaskCounters.cLimit[1] = 1;
  (ModelWithControllersOnly_M)->Timing.TaskCounters.cLimit[2] = 5;

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&ModelWithControllersOnly_M->solverInfo,
                          &ModelWithControllersOnly_M->Timing.simTimeStep);
    rtsiSetTPtr(&ModelWithControllersOnly_M->solverInfo, &rtmGetTPtr
                (ModelWithControllersOnly_M));
    rtsiSetStepSizePtr(&ModelWithControllersOnly_M->solverInfo,
                       &ModelWithControllersOnly_M->Timing.stepSize0);
    rtsiSetdXPtr(&ModelWithControllersOnly_M->solverInfo,
                 &ModelWithControllersOnly_M->derivs);
    rtsiSetContStatesPtr(&ModelWithControllersOnly_M->solverInfo, (real_T **)
                         &ModelWithControllersOnly_M->contStates);
    rtsiSetNumContStatesPtr(&ModelWithControllersOnly_M->solverInfo,
      &ModelWithControllersOnly_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&ModelWithControllersOnly_M->solverInfo,
      &ModelWithControllersOnly_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&ModelWithControllersOnly_M->solverInfo,
      &ModelWithControllersOnly_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&ModelWithControllersOnly_M->solverInfo,
      &ModelWithControllersOnly_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&ModelWithControllersOnly_M->solverInfo,
                          (&rtmGetErrorStatus(ModelWithControllersOnly_M)));
    rtsiSetRTModelPtr(&ModelWithControllersOnly_M->solverInfo,
                      ModelWithControllersOnly_M);
  }

  rtsiSetSimTimeStep(&ModelWithControllersOnly_M->solverInfo, MAJOR_TIME_STEP);
  ModelWithControllersOnly_M->intgData.y = ModelWithControllersOnly_M->odeY;
  ModelWithControllersOnly_M->intgData.f[0] = ModelWithControllersOnly_M->odeF[0];
  ModelWithControllersOnly_M->intgData.f[1] = ModelWithControllersOnly_M->odeF[1];
  ModelWithControllersOnly_M->intgData.f[2] = ModelWithControllersOnly_M->odeF[2];
  ModelWithControllersOnly_M->contStates = ((X_ModelWithControllersOnly_T *)
    &ModelWithControllersOnly_X);
  rtsiSetSolverData(&ModelWithControllersOnly_M->solverInfo, static_cast<void *>
                    (&ModelWithControllersOnly_M->intgData));
  rtsiSetSolverName(&ModelWithControllersOnly_M->solverInfo,"ode3");
  rtmSetTPtr(ModelWithControllersOnly_M,
             &ModelWithControllersOnly_M->Timing.tArray[0]);
  ModelWithControllersOnly_M->Timing.stepSize0 = 0.01;
  rtmSetFirstInitCond(ModelWithControllersOnly_M, 1);

  /* block I/O */
  (void) std::memset((static_cast<void *>(&ModelWithControllersOnly_B)), 0,
                     sizeof(B_ModelWithControllersOnly_T));

  /* states (continuous) */
  {
    (void) std::memset(static_cast<void *>(&ModelWithControllersOnly_X), 0,
                       sizeof(X_ModelWithControllersOnly_T));
  }

  /* states (dwork) */
  (void) std::memset(static_cast<void *>(&ModelWithControllersOnly_DW), 0,
                     sizeof(DW_ModelWithControllersOnly_T));

  /* Start for RateTransition generated from: '<S12>/Equal1' */
  ModelWithControllersOnly_B.TmpRTBAtEqual1Inport1 =
    ModelWithControllersOnly_cal->TmpRTBAtEqual1Inport1_InitialCo;

  /* Start for RateTransition generated from: '<S12>/Equal1' */
  rtw_slrealtime_mutex_init
    (&ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_d0_SEMAPH);

  /* Start for RateTransition generated from: '<S12>/Equal2' */
  ModelWithControllersOnly_B.TmpRTBAtEqual2Inport1 =
    ModelWithControllersOnly_cal->TmpRTBAtEqual2Inport1_InitialCo;

  /* Start for RateTransition generated from: '<S12>/Equal2' */
  rtw_slrealtime_mutex_init
    (&ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_d0_SEMAPH);

  /* Start for RateTransition generated from: '<S12>/Equal4' */
  ModelWithControllersOnly_B.TmpRTBAtEqual4Inport2 =
    ModelWithControllersOnly_cal->TmpRTBAtEqual4Inport2_InitialCo;

  /* Start for RateTransition generated from: '<S12>/Equal4' */
  rtw_slrealtime_mutex_init
    (&ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_d0_SEMAPH);

  /* Start for RateTransition generated from: '<S12>/Equal5' */
  ModelWithControllersOnly_B.TmpRTBAtEqual5Inport2 =
    ModelWithControllersOnly_cal->TmpRTBAtEqual5Inport2_InitialCo;

  /* Start for RateTransition generated from: '<S12>/Equal5' */
  rtw_slrealtime_mutex_init
    (&ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_d0_SEMAPH);

  /* Start for RateTransition generated from: '<S10>/Multiply' */
  ModelWithControllersOnly_B.TmpRTBAtMultiplyInport1 =
    ModelWithControllersOnly_cal->TmpRTBAtMultiplyInport1_Initial;

  /* Start for RateTransition generated from: '<S10>/Multiply' */
  rtw_slrealtime_mutex_init
    (&ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_d0_SEMA);

  /* Start for RateTransition generated from: '<S11>/Switch Case' */
  ModelWithControllersOnly_B.TmpRTBAtSwitchCaseInport1 =
    ModelWithControllersOnly_cal->TmpRTBAtSwitchCaseInport1_Initi;

  /* Start for RateTransition generated from: '<S11>/Switch Case' */
  rtw_slrealtime_mutex_init
    (&ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_d0_SE);

  /* Start for RateTransition generated from: '<S119>/Integrator' */
  ModelWithControllersOnly_B.TmpRTBAtIntegratorInport2 =
    ModelWithControllersOnly_cal->TmpRTBAtIntegratorInport2_Initi;

  /* Start for RateTransition generated from: '<S119>/Sum' */
  ModelWithControllersOnly_B.TmpRTBAtSumInport3 =
    ModelWithControllersOnly_cal->TmpRTBAtSumInport3_InitialCondi;

  /* Start for RateTransition generated from: '<S118>/Integrator' */
  ModelWithControllersOnly_B.TmpRTBAtIntegratorInport2_e =
    ModelWithControllersOnly_cal->TmpRTBAtIntegratorInport2_Ini_j;

  /* Start for RateTransition generated from: '<S118>/Sum' */
  ModelWithControllersOnly_B.TmpRTBAtSumInport3_l =
    ModelWithControllersOnly_cal->TmpRTBAtSumInport3_InitialCon_j;

  /* Start for RateTransition generated from: '<S118>/Integrator' */
  ModelWithControllersOnly_B.TmpRTBAtIntegratorInport1 =
    ModelWithControllersOnly_cal->TmpRTBAtIntegratorInport1_Initi;

  /* Start for RateTransition generated from: '<S119>/Integrator' */
  ModelWithControllersOnly_B.TmpRTBAtIntegratorInport1_c =
    ModelWithControllersOnly_cal->TmpRTBAtIntegratorInport1_Ini_j;

  /* Start for Constant: '<S5>/Constant9' */
  ModelWithControllersOnly_B.gear_cmd =
    ModelWithControllersOnly_cal->Constant9_Value;

  /* Start for RateTransition generated from: '<S120>/MATLAB Function1' */
  ModelWithControllersOnly_B.TmpRTBAtMATLABFunction1Inport3 =
    ModelWithControllersOnly_cal->TmpRTBAtMATLABFunction1Inport3_;

  /* Start for RateTransition generated from: '<S120>/MATLAB Function1' */
  rtw_slrealtime_mutex_init
    (&ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_k);

  /* Start for RateTransition generated from: '<S120>/MATLAB Function' */
  ModelWithControllersOnly_B.TmpRTBAtMATLABFunctionInport3 =
    ModelWithControllersOnly_cal->TmpRTBAtMATLABFunctionInport3_I;

  /* Start for RateTransition generated from: '<S120>/MATLAB Function' */
  rtw_slrealtime_mutex_init
    (&ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_d);

  /* Start for Lookup_n-D: '<S218>/1-D Lookup Table' */
  {
    rt_LUTnWork *TWork_start = (rt_LUTnWork *)
      &ModelWithControllersOnly_DW.TWork[0];
    void **bpDataSet = static_cast<void **>
      (&ModelWithControllersOnly_DW.m_bpDataSet);
    TWork_start->m_dimSizes = static_cast<const uint32_T *>
      (&ModelWithControllersOnly_cal->uDLookupTable_dimSizes);
    TWork_start->m_tableData = (void *)
      ModelWithControllersOnly_cal->LookupGain_tbl;
    TWork_start->m_bpDataSet = bpDataSet;
    TWork_start->m_bpIndex = &ModelWithControllersOnly_DW.m_bpIndex;
    TWork_start->m_bpLambda = &ModelWithControllersOnly_DW.m_bpLambda;
    TWork_start->m_maxIndex = static_cast<const uint32_T *>
      (&ModelWithControllersOnly_cal->uDLookupTable_maxIndex);
    bpDataSet[0] = (void *) ModelWithControllersOnly_cal->LookupGain_bpts;
  }

  {
    const real_T **bpDataSet;
    const real_T *xp, *yp;
    real_T *dp;
    uint32_T len;
    const rt_LUTnWork *TWork_interp;
    rt_LUTSplineWork *rt_SplWk = (rt_LUTSplineWork*)
      &ModelWithControllersOnly_DW.SWork[0];
    rt_SplWk->m_TWork = (rt_LUTnWork*)&ModelWithControllersOnly_DW.TWork[0];
    rt_SplWk->m_yyA = &ModelWithControllersOnly_DW.m_yyA;
    rt_SplWk->m_yyB = &ModelWithControllersOnly_DW.m_yyB;
    rt_SplWk->m_yy2 = &ModelWithControllersOnly_DW.m_yy2;
    rt_SplWk->m_up = &ModelWithControllersOnly_DW.m_up[0];
    rt_SplWk->m_y2 = &ModelWithControllersOnly_DW.m_y2[0];
    rt_SplWk->m_numYWorkElts =
      ModelWithControllersOnly_cal->uDLookupTable_numYWorkElts;
    rt_SplWk->m_reCalc =
      &ModelWithControllersOnly_DW.reCalcSecDerivFirstDimCoeffs;
    rt_SplWk->m_preBp0AndTable =
      &ModelWithControllersOnly_DW.prevBp0AndTableData[0];
    *rt_SplWk->m_reCalc = 1;

    /* cache table data and first breakpoint data */
    TWork_interp = static_cast<const rt_LUTnWork *>(rt_SplWk->m_TWork);
    bpDataSet = (const real_T **) TWork_interp->m_bpDataSet;
    xp = bpDataSet[0U];
    len = TWork_interp->m_maxIndex[0U] + 1U;
    dp = (real_T *) rt_SplWk->m_preBp0AndTable;
    yp = (real_T *) TWork_interp->m_tableData;
    (void) std::memcpy(dp, xp,
                       len * sizeof(real_T));
    dp = &(dp[len]);

    /* save the table data */
    (void) std::memcpy(dp, yp,
                       len * rt_SplWk->m_numYWorkElts[0U] * sizeof(real_T));
  }

  /* Start for FromWorkspace: '<S5>/From Workspace' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06,
      0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18,
      0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31,
      0.32, 0.33, 0.34, 0.35000000000000003, 0.36, 0.37, 0.38, 0.39, 0.4,
      0.41000000000000003, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47000000000000003,
      0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57000000000000006,
      0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68,
      0.69000000000000006, 0.70000000000000007, 0.71, 0.72, 0.73, 0.74, 0.75,
      0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82000000000000006,
      0.83000000000000007, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92,
      0.93, 0.94000000000000006, 0.95000000000000007, 0.96, 0.97, 0.98, 0.99,
      1.0, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.12,
      1.1300000000000001, 1.1400000000000001, 1.1500000000000001, 1.16, 1.17,
      1.18, 1.19, 1.2, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3,
      1.31, 1.32, 1.33, 1.34, 1.35, 1.36, 1.37, 1.3800000000000001,
      1.3900000000000001, 1.4000000000000001, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46,
      1.47, 1.48, 1.49, 1.5, 1.51, 1.52, 1.53, 1.54, 1.55, 1.56, 1.57, 1.58,
      1.59, 1.6, 1.61, 1.62, 1.6300000000000001, 1.6400000000000001,
      1.6500000000000001, 1.6600000000000001, 1.67, 1.68, 1.69, 1.7, 1.71, 1.72,
      1.73, 1.74, 1.75, 1.76, 1.77, 1.78, 1.79, 1.8, 1.81, 1.82, 1.83, 1.84,
      1.85, 1.86, 1.87, 1.8800000000000001, 1.8900000000000001,
      1.9000000000000001, 1.9100000000000001, 1.92, 1.93, 1.94, 1.95, 1.96, 1.97,
      1.98, 1.99, 2.0, 2.0100000000000002, 2.02, 2.0300000000000002, 2.04, 2.05,
      2.06, 2.07, 2.08, 2.09, 2.1, 2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17,
      2.18, 2.19, 2.2, 2.21, 2.22, 2.23, 2.24, 2.25, 2.2600000000000002, 2.27,
      2.2800000000000002, 2.29, 2.3000000000000003, 2.31, 2.32, 2.33, 2.34, 2.35,
      2.36, 2.37, 2.38, 2.39, 2.4, 2.41, 2.42, 2.43, 2.44, 2.45, 2.46, 2.47,
      2.48, 2.49, 2.5, 2.5100000000000002, 2.52, 2.5300000000000002, 2.54,
      2.5500000000000003, 2.56, 2.57, 2.58, 2.59, 2.6, 2.61, 2.62, 2.63, 2.64,
      2.65, 2.66, 2.67, 2.68, 2.69, 2.7, 2.71, 2.72, 2.73, 2.74, 2.75,
      2.7600000000000002, 2.77, 2.7800000000000002, 2.79, 2.8000000000000003,
      2.81, 2.82, 2.83, 2.84, 2.85, 2.86, 2.87, 2.88, 2.89, 2.9, 2.91, 2.92,
      2.93, 2.94, 2.95, 2.96, 2.97, 2.98, 2.99, 3.0, 3.0100000000000002, 3.02,
      3.0300000000000002, 3.04, 3.0500000000000003, 3.06, 3.0700000000000003,
      3.08, 3.09, 3.1, 3.11, 3.12, 3.13, 3.14, 3.15, 3.16, 3.17, 3.18, 3.19, 3.2,
      3.21, 3.22, 3.23, 3.24, 3.25, 3.2600000000000002, 3.27, 3.2800000000000002,
      3.29, 3.3000000000000003, 3.31, 3.3200000000000003, 3.33, 3.34, 3.35, 3.36,
      3.37, 3.38, 3.39, 3.4, 3.41, 3.42, 3.43, 3.44, 3.45, 3.46, 3.47, 3.48,
      3.49, 3.5, 3.5100000000000002, 3.52, 3.5300000000000002, 3.54,
      3.5500000000000003, 3.56, 3.5700000000000003, 3.58, 3.59, 3.6, 3.61, 3.62,
      3.63, 3.64, 3.65, 3.66, 3.67, 3.68, 3.69, 3.7, 3.71, 3.72, 3.73, 3.74,
      3.75, 3.7600000000000002, 3.77, 3.7800000000000002, 3.79,
      3.8000000000000003, 3.81, 3.8200000000000003, 3.83, 3.84, 3.85, 3.86, 3.87,
      3.88, 3.89, 3.9, 3.91, 3.92, 3.93, 3.94, 3.95, 3.96, 3.97, 3.98, 3.99, 4.0,
      4.01, 4.0200000000000005, 4.03, 4.04, 4.05, 4.0600000000000005, 4.07, 4.08,
      4.09, 4.1, 4.11, 4.12, 4.13, 4.14, 4.15, 4.16, 4.17, 4.18, 4.19, 4.2, 4.21,
      4.22, 4.23, 4.24, 4.25, 4.26, 4.2700000000000005, 4.28, 4.29, 4.3,
      4.3100000000000005, 4.32, 4.33, 4.34, 4.3500000000000005, 4.36, 4.37, 4.38,
      4.39, 4.4, 4.41, 4.42, 4.43, 4.44, 4.45, 4.46, 4.47, 4.48, 4.49, 4.5, 4.51,
      4.5200000000000005, 4.53, 4.54, 4.55, 4.5600000000000005, 4.57, 4.58, 4.59,
      4.6000000000000005, 4.61, 4.62, 4.63, 4.64, 4.65, 4.66, 4.67, 4.68, 4.69,
      4.7, 4.71, 4.72, 4.73, 4.74, 4.75, 4.76, 4.7700000000000005, 4.78, 4.79,
      4.8, 4.8100000000000005, 4.82, 4.83, 4.84, 4.8500000000000005, 4.86, 4.87,
      4.88, 4.89, 4.9, 4.91, 4.92, 4.93, 4.94, 4.95, 4.96, 4.97, 4.98, 4.99, 5.0,
      5.01, 5.0200000000000005, 5.03, 5.04, 5.05, 5.0600000000000005, 5.07, 5.08,
      5.09, 5.1000000000000005, 5.11, 5.12, 5.13, 5.14, 5.15, 5.16, 5.17, 5.18,
      5.19, 5.2, 5.21, 5.22, 5.23, 5.24, 5.25, 5.26, 5.2700000000000005, 5.28,
      5.29, 5.3, 5.3100000000000005, 5.32, 5.33, 5.34, 5.3500000000000005, 5.36,
      5.37, 5.38, 5.39, 5.4, 5.41, 5.42, 5.43, 5.44, 5.45, 5.46, 5.47, 5.48,
      5.49, 5.5, 5.51, 5.5200000000000005, 5.53, 5.54, 5.55, 5.5600000000000005,
      5.57, 5.58, 5.59, 5.6000000000000005, 5.61, 5.62, 5.63, 5.64, 5.65, 5.66,
      5.67, 5.68, 5.69, 5.7, 5.71, 5.72, 5.73, 5.74, 5.75, 5.76,
      5.7700000000000005, 5.78, 5.79, 5.8, 5.8100000000000005, 5.82, 5.83, 5.84,
      5.8500000000000005, 5.86, 5.87, 5.88, 5.89, 5.9, 5.91, 5.92, 5.93, 5.94,
      5.95, 5.96, 5.97, 5.98, 5.99, 6.0, 6.01, 6.0200000000000005, 6.03, 6.04,
      6.05, 6.0600000000000005, 6.07, 6.08, 6.09, 6.1000000000000005, 6.11, 6.12,
      6.13, 6.1400000000000006, 6.15, 6.16, 6.17, 6.18, 6.19, 6.2, 6.21, 6.22,
      6.23, 6.24, 6.25, 6.26, 6.2700000000000005, 6.28, 6.29, 6.3,
      6.3100000000000005, 6.32, 6.33, 6.34, 6.3500000000000005, 6.36, 6.37, 6.38,
      6.3900000000000006, 6.4, 6.41, 6.42, 6.43, 6.44, 6.45, 6.46, 6.47, 6.48,
      6.49, 6.5, 6.51, 6.5200000000000005, 6.53, 6.54, 6.55, 6.5600000000000005,
      6.57, 6.58, 6.59, 6.6000000000000005, 6.61, 6.62, 6.63, 6.6400000000000006,
      6.65, 6.66, 6.67, 6.68, 6.69, 6.7, 6.71, 6.72, 6.73, 6.74, 6.75, 6.76,
      6.7700000000000005, 6.78, 6.79, 6.8, 6.8100000000000005, 6.82, 6.83, 6.84,
      6.8500000000000005, 6.86, 6.87, 6.88, 6.8900000000000006, 6.9, 6.91, 6.92,
      6.93, 6.94, 6.95, 6.96, 6.97, 6.98, 6.99, 7.0, 7.01, 7.0200000000000005,
      7.03, 7.04, 7.05, 7.0600000000000005, 7.07, 7.08, 7.09, 7.1000000000000005,
      7.11, 7.12, 7.13, 7.1400000000000006, 7.15, 7.16, 7.17, 7.18, 7.19, 7.2,
      7.21, 7.22, 7.23, 7.24, 7.25, 7.26, 7.2700000000000005, 7.28, 7.29, 7.3,
      7.3100000000000005, 7.32, 7.33, 7.34, 7.3500000000000005, 7.36, 7.37, 7.38,
      7.3900000000000006, 7.4, 7.41, 7.42, 7.43, 7.44, 7.45, 7.46, 7.47, 7.48,
      7.49, 7.5, 7.51, 7.5200000000000005, 7.53, 7.54, 7.55, 7.5600000000000005,
      7.57, 7.58, 7.59, 7.6000000000000005, 7.61, 7.62, 7.63, 7.6400000000000006,
      7.65, 7.66, 7.67, 7.68, 7.69, 7.7, 7.71, 7.72, 7.73, 7.74, 7.75, 7.76,
      7.7700000000000005, 7.78, 7.79, 7.8, 7.8100000000000005, 7.82, 7.83, 7.84,
      7.8500000000000005, 7.86, 7.87, 7.88, 7.8900000000000006, 7.9, 7.91, 7.92,
      7.9300000000000006, 7.94, 7.95, 7.96, 7.97, 7.98, 7.99, 8.0, 8.01, 8.02,
      8.03, 8.0400000000000009, 8.05, 8.06, 8.07, 8.08, 8.09, 8.1, 8.11,
      8.120000000000001, 8.13, 8.14, 8.15, 8.16, 8.17, 8.18, 8.19, 8.2, 8.21,
      8.22, 8.23, 8.24, 8.25, 8.26, 8.27, 8.28, 8.2900000000000009, 8.3, 8.31,
      8.32, 8.33, 8.34, 8.35, 8.36, 8.370000000000001, 8.38, 8.39, 8.4, 8.41,
      8.42, 8.43, 8.44, 8.45, 8.46, 8.47, 8.48, 8.49, 8.5, 8.51, 8.52, 8.53,
      8.5400000000000009, 8.55, 8.56, 8.57, 8.58, 8.59, 8.6, 8.61,
      8.620000000000001, 8.63, 8.64, 8.65, 8.66, 8.67, 8.68, 8.69,
      8.7000000000000011, 8.71, 8.72, 8.73, 8.74, 8.75, 8.76, 8.77, 8.78,
      8.7900000000000009, 8.8, 8.81, 8.82, 8.83, 8.84, 8.85, 8.86,
      8.870000000000001, 8.88, 8.89, 8.9, 8.91, 8.92, 8.93, 8.94,
      8.9500000000000011, 8.96, 8.97, 8.98, 8.99, 9.0, 9.01, 9.02, 9.03,
      9.0400000000000009, 9.05, 9.06, 9.07, 9.08, 9.09, 9.1, 9.11,
      9.120000000000001, 9.13, 9.14, 9.15, 9.16, 9.17, 9.18, 9.19,
      9.2000000000000011, 9.21, 9.22, 9.23, 9.24, 9.25, 9.26, 9.27, 9.28,
      9.2900000000000009, 9.3, 9.31, 9.32, 9.33, 9.34, 9.35, 9.36,
      9.370000000000001, 9.38, 9.39, 9.4, 9.41, 9.42, 9.43, 9.44,
      9.4500000000000011, 9.46, 9.47, 9.48, 9.49, 9.5, 9.51, 9.52, 9.53,
      9.5400000000000009, 9.55, 9.56, 9.57, 9.58, 9.59, 9.6, 9.61,
      9.620000000000001, 9.63, 9.64, 9.65, 9.66, 9.67, 9.68, 9.69,
      9.7000000000000011, 9.71, 9.72, 9.73, 9.74, 9.75, 9.76, 9.77, 9.78,
      9.7900000000000009, 9.8, 9.81, 9.82, 9.83, 9.84, 9.85, 9.86,
      9.870000000000001, 9.88, 9.89, 9.9, 9.91, 9.92, 9.93, 9.94,
      9.9500000000000011, 9.96, 9.97, 9.98, 9.99, 10.0, 10.01, 10.02, 10.03,
      10.040000000000001, 10.05, 10.06, 10.07, 10.08, 10.09, 10.1, 10.11,
      10.120000000000001, 10.13, 10.14, 10.15, 10.16, 10.17, 10.18, 10.19,
      10.200000000000001, 10.21, 10.22, 10.23, 10.24, 10.25, 10.26, 10.27, 10.28,
      10.290000000000001, 10.3, 10.31, 10.32, 10.33, 10.34, 10.35, 10.36,
      10.370000000000001, 10.38, 10.39, 10.4, 10.41, 10.42, 10.43, 10.44,
      10.450000000000001, 10.46, 10.47, 10.48, 10.49, 10.5, 10.51, 10.52, 10.53,
      10.540000000000001, 10.55, 10.56, 10.57, 10.58, 10.59, 10.6, 10.61,
      10.620000000000001, 10.63, 10.64, 10.65, 10.66, 10.67, 10.68, 10.69,
      10.700000000000001, 10.71, 10.72, 10.73, 10.74, 10.75, 10.76, 10.77, 10.78,
      10.790000000000001, 10.8, 10.81, 10.82, 10.83, 10.84, 10.85, 10.86,
      10.870000000000001, 10.88, 10.89, 10.9, 10.91, 10.92, 10.93, 10.94,
      10.950000000000001, 10.96, 10.97, 10.98, 10.99, 11.0, 11.01, 11.02, 11.03,
      11.040000000000001, 11.05, 11.06, 11.07, 11.08, 11.09, 11.1, 11.11,
      11.120000000000001, 11.13, 11.14, 11.15, 11.16, 11.17, 11.18, 11.19,
      11.200000000000001, 11.21, 11.22, 11.23, 11.24, 11.25, 11.26, 11.27, 11.28,
      11.290000000000001, 11.3, 11.31, 11.32, 11.33, 11.34, 11.35, 11.36,
      11.370000000000001, 11.38, 11.39, 11.4, 11.41, 11.42, 11.43, 11.44,
      11.450000000000001, 11.46, 11.47, 11.48, 11.49, 11.5, 11.51, 11.52, 11.53,
      11.540000000000001, 11.55, 11.56, 11.57, 11.58, 11.59, 11.6, 11.61,
      11.620000000000001, 11.63, 11.64, 11.65, 11.66, 11.67, 11.68, 11.69,
      11.700000000000001, 11.71, 11.72, 11.73, 11.74, 11.75, 11.76, 11.77, 11.78,
      11.790000000000001, 11.8, 11.81, 11.82, 11.83, 11.84, 11.85, 11.86,
      11.870000000000001, 11.88, 11.89, 11.9, 11.91, 11.92, 11.93, 11.94,
      11.950000000000001, 11.96, 11.97, 11.98, 11.99, 12.0, 12.01, 12.02,
      12.030000000000001, 12.040000000000001, 12.05, 12.06, 12.07, 12.08, 12.09,
      12.1, 12.11, 12.120000000000001, 12.13, 12.14, 12.15, 12.16, 12.17, 12.18,
      12.19, 12.200000000000001, 12.21, 12.22, 12.23, 12.24, 12.25, 12.26, 12.27,
      12.280000000000001, 12.290000000000001, 12.3, 12.31, 12.32, 12.33, 12.34,
      12.35, 12.36, 12.370000000000001, 12.38, 12.39, 12.4, 12.41, 12.42, 12.43,
      12.44, 12.450000000000001, 12.46, 12.47, 12.48, 12.49, 12.5, 12.51, 12.52,
      12.530000000000001, 12.540000000000001, 12.55, 12.56, 12.57, 12.58, 12.59,
      12.6, 12.61, 12.620000000000001, 12.63, 12.64, 12.65, 12.66, 12.67, 12.68,
      12.69, 12.700000000000001, 12.71, 12.72, 12.73, 12.74, 12.75, 12.76, 12.77,
      12.780000000000001, 12.790000000000001, 12.8, 12.81, 12.82, 12.83, 12.84,
      12.85, 12.86, 12.870000000000001, 12.88, 12.89, 12.9, 12.91, 12.92, 12.93,
      12.94, 12.950000000000001, 12.96, 12.97, 12.98, 12.99, 13.0, 13.01, 13.02,
      13.030000000000001, 13.040000000000001, 13.05, 13.06, 13.07, 13.08, 13.09,
      13.1, 13.11, 13.120000000000001, 13.13, 13.14, 13.15, 13.16, 13.17, 13.18,
      13.19, 13.200000000000001, 13.21, 13.22, 13.23, 13.24, 13.25, 13.26, 13.27,
      13.280000000000001, 13.290000000000001, 13.3, 13.31, 13.32, 13.33, 13.34,
      13.35, 13.36, 13.370000000000001, 13.38, 13.39, 13.4, 13.41, 13.42, 13.43,
      13.44, 13.450000000000001, 13.46, 13.47, 13.48, 13.49, 13.5, 13.51, 13.52,
      13.530000000000001, 13.540000000000001, 13.55, 13.56, 13.57, 13.58, 13.59,
      13.6, 13.61, 13.620000000000001, 13.63, 13.64, 13.65, 13.66, 13.67, 13.68,
      13.69, 13.700000000000001, 13.71, 13.72, 13.73, 13.74, 13.75, 13.76, 13.77,
      13.780000000000001, 13.790000000000001, 13.8, 13.81, 13.82, 13.83, 13.84,
      13.85, 13.86, 13.870000000000001, 13.88, 13.89, 13.9, 13.91, 13.92, 13.93,
      13.94, 13.950000000000001, 13.96, 13.97, 13.98, 13.99, 14.0, 14.01, 14.02,
      14.030000000000001, 14.040000000000001, 14.05, 14.06, 14.07, 14.08, 14.09,
      14.1, 14.11, 14.120000000000001, 14.13, 14.14, 14.15, 14.16, 14.17, 14.18,
      14.19, 14.200000000000001, 14.21, 14.22, 14.23, 14.24, 14.25, 14.26, 14.27,
      14.280000000000001, 14.290000000000001, 14.3, 14.31, 14.32, 14.33, 14.34,
      14.35, 14.36, 14.370000000000001, 14.38, 14.39, 14.4, 14.41, 14.42, 14.43,
      14.44, 14.450000000000001, 14.46, 14.47, 14.48, 14.49, 14.5, 14.51, 14.52,
      14.530000000000001, 14.540000000000001, 14.55, 14.56, 14.57, 14.58, 14.59,
      14.6, 14.61, 14.620000000000001, 14.63, 14.64, 14.65, 14.66, 14.67, 14.68,
      14.69, 14.700000000000001, 14.71, 14.72, 14.73, 14.74, 14.75, 14.76, 14.77,
      14.780000000000001, 14.790000000000001, 14.8, 14.81, 14.82, 14.83, 14.84,
      14.85, 14.86, 14.870000000000001, 14.88, 14.89, 14.9, 14.91, 14.92, 14.93,
      14.94, 14.950000000000001, 14.96, 14.97, 14.98, 14.99, 15.0, 15.01, 15.02,
      15.030000000000001, 15.040000000000001, 15.05, 15.06, 15.07, 15.08, 15.09,
      15.1, 15.11, 15.120000000000001, 15.13, 15.14, 15.15, 15.16, 15.17, 15.18,
      15.19, 15.200000000000001, 15.21, 15.22, 15.23, 15.24, 15.25, 15.26, 15.27,
      15.280000000000001, 15.290000000000001, 15.3, 15.31, 15.32, 15.33, 15.34,
      15.35, 15.36, 15.370000000000001, 15.38, 15.39, 15.4, 15.41, 15.42, 15.43,
      15.44, 15.450000000000001, 15.46, 15.47, 15.48, 15.49, 15.5, 15.51, 15.52,
      15.530000000000001, 15.540000000000001, 15.55, 15.56, 15.57, 15.58, 15.59,
      15.6, 15.610000000000001, 15.620000000000001, 15.63, 15.64, 15.65, 15.66,
      15.67, 15.68, 15.69, 15.700000000000001, 15.71, 15.72, 15.73, 15.74, 15.75,
      15.76, 15.77, 15.780000000000001, 15.790000000000001, 15.8, 15.81, 15.82,
      15.83, 15.84, 15.85, 15.860000000000001, 15.870000000000001, 15.88, 15.89,
      15.9, 15.91, 15.92, 15.93, 15.94, 15.950000000000001, 15.96, 15.97, 15.98,
      15.99, 16.0, 16.01, 16.02, 16.03, 16.04, 16.05, 16.06, 16.07,
      16.080000000000002, 16.09, 16.1, 16.11, 16.12, 16.13, 16.14, 16.15, 16.16,
      16.17, 16.18, 16.19, 16.2, 16.21, 16.22, 16.23, 16.240000000000002, 16.25,
      16.26, 16.27, 16.28, 16.29, 16.3, 16.31, 16.32, 16.330000000000002, 16.34,
      16.35, 16.36, 16.37, 16.38, 16.39, 16.4, 16.41, 16.42, 16.43, 16.44, 16.45,
      16.46, 16.47, 16.48, 16.490000000000002, 16.5, 16.51, 16.52, 16.53, 16.54,
      16.55, 16.56, 16.57, 16.580000000000002, 16.59, 16.6, 16.61, 16.62, 16.63,
      16.64, 16.65, 16.66, 16.67, 16.68, 16.69, 16.7, 16.71, 16.72, 16.73,
      16.740000000000002, 16.75, 16.76, 16.77, 16.78, 16.79, 16.8, 16.81, 16.82,
      16.830000000000002, 16.84, 16.85, 16.86, 16.87, 16.88, 16.89, 16.9, 16.91,
      16.92, 16.93, 16.94, 16.95, 16.96, 16.97, 16.98, 16.990000000000002, 17.0,
      17.01, 17.02, 17.03, 17.04, 17.05, 17.06, 17.07, 17.080000000000002, 17.09,
      17.1, 17.11, 17.12, 17.13, 17.14, 17.150000000000002, 17.16, 17.17, 17.18,
      17.19, 17.2, 17.21, 17.22, 17.23, 17.240000000000002, 17.25, 17.26, 17.27,
      17.28, 17.29, 17.3, 17.31, 17.32, 17.330000000000002, 17.34, 17.35, 17.36,
      17.37, 17.38, 17.39, 17.400000000000002, 17.41, 17.42, 17.43, 17.44, 17.45,
      17.46, 17.47, 17.48, 17.490000000000002, 17.5, 17.51, 17.52, 17.53, 17.54,
      17.55, 17.56, 17.57, 17.580000000000002, 17.59, 17.6, 17.61, 17.62, 17.63,
      17.64, 17.650000000000002, 17.66, 17.67, 17.68, 17.69, 17.7, 17.71, 17.72,
      17.73, 17.740000000000002, 17.75, 17.76, 17.77, 17.78, 17.79, 17.8, 17.81,
      17.82, 17.830000000000002, 17.84, 17.85, 17.86, 17.87, 17.88, 17.89,
      17.900000000000002, 17.91, 17.92, 17.93, 17.94, 17.95, 17.96, 17.97, 17.98,
      17.990000000000002, 18.0, 18.01, 18.02, 18.03, 18.04, 18.05, 18.06, 18.07,
      18.080000000000002, 18.09, 18.1, 18.11, 18.12, 18.13, 18.14,
      18.150000000000002, 18.16, 18.17, 18.18, 18.19, 18.2, 18.21, 18.22, 18.23,
      18.240000000000002, 18.25, 18.26, 18.27, 18.28, 18.29, 18.3, 18.31, 18.32,
      18.330000000000002, 18.34, 18.35, 18.36, 18.37, 18.38, 18.39,
      18.400000000000002, 18.41, 18.42, 18.43, 18.44, 18.45, 18.46, 18.47, 18.48,
      18.490000000000002, 18.5, 18.51, 18.52, 18.53, 18.54, 18.55, 18.56, 18.57,
      18.580000000000002, 18.59, 18.6, 18.61, 18.62, 18.63, 18.64,
      18.650000000000002, 18.66, 18.67, 18.68, 18.69, 18.7, 18.71, 18.72, 18.73,
      18.740000000000002, 18.75, 18.76, 18.77, 18.78, 18.79, 18.8, 18.81, 18.82,
      18.830000000000002, 18.84, 18.85, 18.86, 18.87, 18.88, 18.89,
      18.900000000000002, 18.91, 18.92, 18.93, 18.94, 18.95, 18.96, 18.97, 18.98,
      18.990000000000002, 19.0, 19.01, 19.02, 19.03, 19.04, 19.05, 19.06, 19.07,
      19.080000000000002, 19.09, 19.1, 19.11, 19.12, 19.13, 19.14,
      19.150000000000002, 19.16, 19.17, 19.18, 19.19, 19.2, 19.21, 19.22, 19.23,
      19.240000000000002, 19.25, 19.26, 19.27, 19.28, 19.29, 19.3, 19.31, 19.32,
      19.330000000000002, 19.34, 19.35, 19.36, 19.37, 19.38, 19.39,
      19.400000000000002, 19.41, 19.42, 19.43, 19.44, 19.45, 19.46, 19.47, 19.48,
      19.490000000000002, 19.5, 19.51, 19.52, 19.53, 19.54, 19.55, 19.56, 19.57,
      19.580000000000002, 19.59, 19.6, 19.61, 19.62, 19.63, 19.64,
      19.650000000000002, 19.66, 19.67, 19.68, 19.69, 19.7, 19.71, 19.72, 19.73,
      19.740000000000002, 19.75, 19.76, 19.77, 19.78, 19.79, 19.8, 19.81, 19.82,
      19.830000000000002, 19.84, 19.85, 19.86, 19.87, 19.88, 19.89,
      19.900000000000002, 19.91, 19.92, 19.93, 19.94, 19.95, 19.96, 19.97, 19.98,
      19.990000000000002, 20.0, 20.01, 20.02, 20.03, 20.04, 20.05, 20.06, 20.07,
      20.080000000000002, 20.09, 20.1, 20.11, 20.12, 20.13, 20.14,
      20.150000000000002, 20.16, 20.17, 20.18, 20.19, 20.2, 20.21, 20.22, 20.23,
      20.240000000000002, 20.25, 20.26, 20.27, 20.28, 20.29, 20.3, 20.31, 20.32,
      20.330000000000002, 20.34, 20.35, 20.36, 20.37, 20.38, 20.39,
      20.400000000000002, 20.41, 20.42, 20.43, 20.44, 20.45, 20.46, 20.47, 20.48,
      20.490000000000002, 20.5, 20.51, 20.52, 20.53, 20.54, 20.55, 20.56, 20.57,
      20.580000000000002, 20.59, 20.6, 20.61, 20.62, 20.63, 20.64,
      20.650000000000002, 20.66, 20.67, 20.68, 20.69, 20.7, 20.71, 20.72, 20.73,
      20.740000000000002, 20.75, 20.76, 20.77, 20.78, 20.79, 20.8, 20.81, 20.82,
      20.830000000000002, 20.84, 20.85, 20.86, 20.87, 20.88, 20.89,
      20.900000000000002, 20.91, 20.92, 20.93, 20.94, 20.95, 20.96, 20.97, 20.98,
      20.990000000000002, 21.0, 21.01, 21.02, 21.03, 21.04, 21.05, 21.06, 21.07,
      21.080000000000002, 21.09, 21.1, 21.11, 21.12, 21.13, 21.14,
      21.150000000000002, 21.16, 21.17, 21.18, 21.19, 21.2, 21.21, 21.22, 21.23,
      21.240000000000002, 21.25, 21.26, 21.27, 21.28, 21.29, 21.3, 21.31, 21.32,
      21.330000000000002, 21.34, 21.35, 21.36, 21.37, 21.38, 21.39,
      21.400000000000002, 21.41, 21.42, 21.43, 21.44, 21.45, 21.46, 21.47, 21.48,
      21.490000000000002, 21.5, 21.51, 21.52, 21.53, 21.54, 21.55, 21.56, 21.57,
      21.580000000000002, 21.59, 21.6, 21.61, 21.62, 21.63, 21.64,
      21.650000000000002, 21.66, 21.67, 21.68, 21.69, 21.7, 21.71, 21.72, 21.73,
      21.740000000000002, 21.75, 21.76, 21.77, 21.78, 21.79, 21.8, 21.81, 21.82,
      21.830000000000002, 21.84, 21.85, 21.86, 21.87, 21.88, 21.89,
      21.900000000000002, 21.91, 21.92, 21.93, 21.94, 21.95, 21.96, 21.97, 21.98,
      21.990000000000002, 22.0, 22.01, 22.02, 22.03, 22.04, 22.05, 22.06, 22.07,
      22.080000000000002, 22.09, 22.1, 22.11, 22.12, 22.13, 22.14,
      22.150000000000002, 22.16, 22.17, 22.18, 22.19, 22.2, 22.21, 22.22, 22.23,
      22.240000000000002, 22.25, 22.26, 22.27, 22.28, 22.29, 22.3, 22.31, 22.32,
      22.330000000000002, 22.34, 22.35, 22.36, 22.37, 22.38, 22.39,
      22.400000000000002, 22.41, 22.42, 22.43, 22.44, 22.45, 22.46, 22.47, 22.48,
      22.490000000000002, 22.5, 22.51, 22.52, 22.53, 22.54, 22.55, 22.56, 22.57,
      22.580000000000002, 22.59, 22.6, 22.61, 22.62, 22.63, 22.64,
      22.650000000000002, 22.66, 22.67, 22.68, 22.69, 22.7, 22.71, 22.72, 22.73,
      22.740000000000002, 22.75, 22.76, 22.77, 22.78, 22.79, 22.8, 22.81, 22.82,
      22.830000000000002, 22.84, 22.85, 22.86, 22.87, 22.88, 22.89,
      22.900000000000002, 22.91, 22.92, 22.93, 22.94, 22.95, 22.96, 22.97, 22.98,
      22.990000000000002, 23.0, 23.01, 23.02, 23.03, 23.04, 23.05, 23.06, 23.07,
      23.080000000000002, 23.09, 23.1, 23.11, 23.12, 23.13, 23.14,
      23.150000000000002, 23.16, 23.17, 23.18, 23.19, 23.2, 23.21, 23.22, 23.23,
      23.240000000000002, 23.25, 23.26, 23.27, 23.28, 23.29, 23.3, 23.31, 23.32,
      23.330000000000002, 23.34, 23.35, 23.36, 23.37, 23.38, 23.39,
      23.400000000000002, 23.41, 23.42, 23.43, 23.44, 23.45, 23.46, 23.47, 23.48,
      23.490000000000002, 23.5, 23.51, 23.52, 23.53, 23.54, 23.55, 23.56, 23.57,
      23.580000000000002, 23.59, 23.6, 23.61, 23.62, 23.63, 23.64,
      23.650000000000002, 23.66, 23.67, 23.68, 23.69, 23.7, 23.71, 23.72, 23.73,
      23.740000000000002, 23.75, 23.76, 23.77, 23.78, 23.79, 23.8, 23.81, 23.82,
      23.830000000000002, 23.84, 23.85, 23.86, 23.87, 23.88, 23.89,
      23.900000000000002, 23.91, 23.92, 23.93, 23.94, 23.95, 23.96, 23.97, 23.98,
      23.990000000000002, 24.0, 24.01, 24.02, 24.03, 24.04, 24.05,
      24.060000000000002, 24.07, 24.080000000000002, 24.09, 24.1, 24.11, 24.12,
      24.13, 24.14, 24.150000000000002, 24.16, 24.17, 24.18, 24.19, 24.2, 24.21,
      24.22, 24.23, 24.240000000000002, 24.25, 24.26, 24.27, 24.28, 24.29, 24.3,
      24.310000000000002, 24.32, 24.330000000000002, 24.34, 24.35, 24.36, 24.37,
      24.38, 24.39, 24.400000000000002, 24.41, 24.42, 24.43, 24.44, 24.45, 24.46,
      24.47, 24.48, 24.490000000000002, 24.5, 24.51, 24.52, 24.53, 24.54, 24.55,
      24.560000000000002, 24.57, 24.580000000000002, 24.59, 24.6, 24.61, 24.62,
      24.63, 24.64, 24.650000000000002, 24.66, 24.67, 24.68, 24.69, 24.7, 24.71,
      24.72, 24.73, 24.740000000000002, 24.75, 24.76, 24.77, 24.78, 24.79, 24.8,
      24.810000000000002, 24.82, 24.830000000000002, 24.84, 24.85, 24.86, 24.87,
      24.88, 24.89, 24.900000000000002, 24.91, 24.92, 24.93, 24.94, 24.95, 24.96,
      24.97, 24.98, 24.990000000000002, 25.0, 25.01, 25.02, 25.03, 25.04, 25.05,
      25.060000000000002, 25.07, 25.080000000000002, 25.09, 25.1, 25.11, 25.12,
      25.13, 25.14, 25.150000000000002, 25.16, 25.17, 25.18, 25.19, 25.2, 25.21,
      25.22, 25.23, 25.240000000000002, 25.25, 25.26, 25.27, 25.28, 25.29, 25.3,
      25.310000000000002, 25.32, 25.330000000000002, 25.34, 25.35, 25.36, 25.37,
      25.38, 25.39, 25.400000000000002, 25.41, 25.42, 25.43, 25.44, 25.45, 25.46,
      25.47, 25.48, 25.490000000000002, 25.5, 25.51, 25.52, 25.53, 25.54, 25.55,
      25.560000000000002, 25.57, 25.580000000000002, 25.59, 25.6, 25.61, 25.62,
      25.63, 25.64, 25.650000000000002, 25.66, 25.67, 25.68, 25.69, 25.7, 25.71,
      25.72, 25.73, 25.740000000000002, 25.75, 25.76, 25.77, 25.78, 25.79, 25.8,
      25.810000000000002, 25.82, 25.830000000000002, 25.84, 25.85, 25.86, 25.87,
      25.88, 25.89, 25.900000000000002, 25.91, 25.92, 25.93, 25.94, 25.95, 25.96,
      25.97, 25.98, 25.990000000000002, 26.0, 26.01, 26.02, 26.03, 26.04, 26.05,
      26.060000000000002, 26.07, 26.080000000000002, 26.09, 26.1, 26.11, 26.12,
      26.13, 26.14, 26.150000000000002, 26.16, 26.17, 26.18, 26.19, 26.2, 26.21,
      26.22, 26.23, 26.240000000000002, 26.25, 26.26, 26.27, 26.28, 26.29, 26.3,
      26.310000000000002, 26.32, 26.330000000000002, 26.34, 26.35, 26.36, 26.37,
      26.38, 26.39, 26.400000000000002, 26.41, 26.42, 26.43, 26.44, 26.45, 26.46,
      26.47, 26.48, 26.490000000000002, 26.5, 26.51, 26.52, 26.53, 26.54, 26.55,
      26.560000000000002, 26.57, 26.580000000000002, 26.59, 26.6, 26.61, 26.62,
      26.63, 26.64, 26.650000000000002, 26.66, 26.67, 26.68, 26.69, 26.7, 26.71,
      26.72, 26.73, 26.740000000000002, 26.75, 26.76, 26.77, 26.78, 26.79, 26.8,
      26.810000000000002, 26.82, 26.830000000000002, 26.84, 26.85, 26.86, 26.87,
      26.88, 26.89, 26.900000000000002, 26.91, 26.92, 26.93, 26.94, 26.95, 26.96,
      26.97, 26.98, 26.990000000000002, 27.0, 27.01, 27.02, 27.03, 27.04, 27.05,
      27.060000000000002, 27.07, 27.080000000000002, 27.09, 27.1, 27.11, 27.12,
      27.13, 27.14, 27.150000000000002, 27.16, 27.17, 27.18, 27.19, 27.2, 27.21,
      27.22, 27.23, 27.240000000000002, 27.25, 27.26, 27.27, 27.28, 27.29, 27.3,
      27.310000000000002, 27.32, 27.330000000000002, 27.34, 27.35, 27.36, 27.37,
      27.38, 27.39, 27.400000000000002, 27.41, 27.42, 27.43, 27.44, 27.45, 27.46,
      27.47, 27.48, 27.490000000000002, 27.5, 27.51, 27.52, 27.53, 27.54, 27.55,
      27.560000000000002, 27.57, 27.580000000000002, 27.59, 27.6, 27.61, 27.62,
      27.63, 27.64, 27.650000000000002, 27.66, 27.67, 27.68, 27.69, 27.7, 27.71,
      27.72, 27.73, 27.740000000000002, 27.75, 27.76, 27.77, 27.78, 27.79, 27.8,
      27.810000000000002, 27.82, 27.830000000000002, 27.84, 27.85, 27.86, 27.87,
      27.88, 27.89, 27.900000000000002, 27.91, 27.92, 27.93, 27.94, 27.95, 27.96,
      27.97, 27.98, 27.990000000000002, 28.0, 28.01, 28.02, 28.03, 28.04, 28.05,
      28.060000000000002, 28.07, 28.080000000000002, 28.09, 28.1, 28.11, 28.12,
      28.13, 28.14, 28.150000000000002, 28.16, 28.17, 28.18, 28.19, 28.2, 28.21,
      28.22, 28.23, 28.240000000000002, 28.25, 28.26, 28.27, 28.28, 28.29, 28.3,
      28.310000000000002, 28.32, 28.330000000000002, 28.34, 28.35, 28.36, 28.37,
      28.38, 28.39, 28.400000000000002, 28.41, 28.42, 28.43, 28.44, 28.45, 28.46,
      28.47, 28.48, 28.490000000000002, 28.5, 28.51, 28.52, 28.53, 28.54, 28.55,
      28.560000000000002, 28.57, 28.580000000000002, 28.59, 28.6, 28.61, 28.62,
      28.63, 28.64, 28.650000000000002, 28.66, 28.67, 28.68, 28.69, 28.7, 28.71,
      28.72, 28.73, 28.740000000000002, 28.75, 28.76, 28.77, 28.78, 28.79, 28.8,
      28.810000000000002, 28.82, 28.830000000000002, 28.84, 28.85, 28.86, 28.87,
      28.88, 28.89, 28.900000000000002, 28.91, 28.92, 28.93, 28.94, 28.95, 28.96,
      28.97, 28.98, 28.990000000000002, 29.0, 29.01, 29.02, 29.03, 29.04, 29.05,
      29.060000000000002, 29.07, 29.080000000000002, 29.09, 29.1, 29.11, 29.12,
      29.13, 29.14, 29.150000000000002, 29.16, 29.17, 29.18, 29.19, 29.2, 29.21,
      29.22, 29.23, 29.240000000000002, 29.25, 29.26, 29.27, 29.28, 29.29, 29.3,
      29.310000000000002, 29.32, 29.330000000000002, 29.34, 29.35, 29.36, 29.37,
      29.38, 29.39, 29.400000000000002, 29.41, 29.42, 29.43, 29.44, 29.45, 29.46,
      29.47, 29.48, 29.490000000000002, 29.5, 29.51, 29.52, 29.53, 29.54, 29.55,
      29.560000000000002, 29.57, 29.580000000000002, 29.59, 29.6, 29.61, 29.62,
      29.63, 29.64, 29.650000000000002, 29.66, 29.67, 29.68, 29.69, 29.7, 29.71,
      29.72, 29.73, 29.740000000000002, 29.75, 29.76, 29.77, 29.78, 29.79, 29.8,
      29.810000000000002, 29.82, 29.830000000000002, 29.84, 29.85, 29.86, 29.87,
      29.88, 29.89, 29.900000000000002, 29.91, 29.92, 29.93, 29.94, 29.95, 29.96,
      29.97, 29.98, 29.990000000000002, 30.0, 30.01, 30.02, 30.03, 30.04, 30.05,
      30.060000000000002, 30.07, 30.080000000000002, 30.09, 30.1, 30.11, 30.12,
      30.13, 30.14, 30.150000000000002, 30.16, 30.17, 30.18, 30.19, 30.2, 30.21,
      30.22, 30.23, 30.240000000000002, 30.25, 30.26, 30.27, 30.28, 30.29, 30.3,
      30.310000000000002, 30.32, 30.330000000000002, 30.34, 30.35, 30.36, 30.37,
      30.38, 30.39, 30.400000000000002, 30.41, 30.42, 30.43, 30.44, 30.45, 30.46,
      30.47, 30.48, 30.490000000000002, 30.5, 30.51, 30.52, 30.53, 30.54, 30.55,
      30.560000000000002, 30.57, 30.580000000000002, 30.59, 30.6, 30.61, 30.62,
      30.63, 30.64, 30.650000000000002, 30.66, 30.67, 30.68, 30.69, 30.7, 30.71,
      30.72, 30.73, 30.740000000000002, 30.75, 30.76, 30.77, 30.78, 30.79, 30.8,
      30.810000000000002, 30.82, 30.830000000000002, 30.84, 30.85, 30.86, 30.87,
      30.88, 30.89, 30.900000000000002, 30.91, 30.92, 30.93, 30.94, 30.95, 30.96,
      30.970000000000002, 30.98, 30.990000000000002, 31.0, 31.01, 31.02, 31.03,
      31.04, 31.05, 31.060000000000002, 31.07, 31.080000000000002, 31.09, 31.1,
      31.11, 31.12, 31.13, 31.14, 31.150000000000002, 31.16, 31.17, 31.18, 31.19,
      31.2, 31.21, 31.220000000000002, 31.23, 31.240000000000002, 31.25, 31.26,
      31.27, 31.28, 31.29, 31.3, 31.310000000000002, 31.32, 31.330000000000002,
      31.34, 31.35, 31.36, 31.37, 31.38, 31.39, 31.400000000000002, 31.41, 31.42,
      31.43, 31.44, 31.45, 31.46, 31.470000000000002, 31.48, 31.490000000000002,
      31.5, 31.51, 31.52, 31.53, 31.54, 31.55, 31.560000000000002, 31.57,
      31.580000000000002, 31.59, 31.6, 31.61, 31.62, 31.63, 31.64,
      31.650000000000002, 31.66, 31.67, 31.68, 31.69, 31.7, 31.71,
      31.720000000000002, 31.73, 31.740000000000002, 31.75, 31.76, 31.77, 31.78,
      31.79, 31.8, 31.810000000000002, 31.82, 31.830000000000002, 31.84, 31.85,
      31.86, 31.87, 31.88, 31.89, 31.900000000000002, 31.91, 31.92, 31.93, 31.94,
      31.95, 31.96, 31.970000000000002, 31.98, 31.990000000000002, 32.0, 32.01,
      32.02, 32.03, 32.04, 32.05, 32.06, 32.07, 32.08, 32.09, 32.1, 32.11, 32.12,
      32.13, 32.14, 32.15, 32.160000000000004, 32.17, 32.18, 32.19, 32.2, 32.21,
      32.22, 32.230000000000004, 32.24, 32.25, 32.26, 32.27, 32.28, 32.29, 32.3,
      32.31, 32.32, 32.33, 32.34, 32.35, 32.36, 32.37, 32.38, 32.39, 32.4,
      32.410000000000004, 32.42, 32.43, 32.44, 32.45, 32.46, 32.47,
      32.480000000000004, 32.49, 32.5, 32.51, 32.52, 32.53, 32.54, 32.55, 32.56,
      32.57, 32.58, 32.59, 32.6, 32.61, 32.62, 32.63, 32.64, 32.65,
      32.660000000000004, 32.67, 32.68, 32.69, 32.7, 32.71, 32.72,
      32.730000000000004, 32.74, 32.75, 32.76, 32.77, 32.78, 32.79, 32.8, 32.81,
      32.82, 32.83, 32.84, 32.85, 32.86, 32.87, 32.88, 32.89, 32.9,
      32.910000000000004, 32.92, 32.93, 32.94, 32.95, 32.96, 32.97,
      32.980000000000004, 32.99, 33.0, 33.01, 33.02, 33.03, 33.04, 33.05, 33.06,
      33.07, 33.08, 33.09, 33.1, 33.11, 33.12, 33.13, 33.14, 33.15,
      33.160000000000004, 33.17, 33.18, 33.19, 33.2, 33.21, 33.22,
      33.230000000000004, 33.24, 33.25, 33.26, 33.27, 33.28, 33.29, 33.3, 33.31,
      33.32, 33.33, 33.34, 33.35, 33.36, 33.37, 33.38, 33.39, 33.4,
      33.410000000000004, 33.42, 33.43, 33.44, 33.45, 33.46, 33.47,
      33.480000000000004, 33.49, 33.5, 33.51, 33.52, 33.53, 33.54, 33.55, 33.56,
      33.57, 33.58, 33.59, 33.6, 33.61, 33.62, 33.63, 33.64, 33.65,
      33.660000000000004, 33.67, 33.68, 33.69, 33.7, 33.71, 33.72,
      33.730000000000004, 33.74, 33.75, 33.76, 33.77, 33.78, 33.79, 33.8, 33.81,
      33.82, 33.83, 33.84, 33.85, 33.86, 33.87, 33.88, 33.89, 33.9,
      33.910000000000004, 33.92, 33.93, 33.94, 33.95, 33.96, 33.97,
      33.980000000000004, 33.99, 34.0, 34.01, 34.02, 34.03, 34.04, 34.05, 34.06,
      34.07, 34.08, 34.09, 34.1, 34.11, 34.12, 34.13, 34.14, 34.15,
      34.160000000000004, 34.17, 34.18, 34.19, 34.2, 34.21, 34.22,
      34.230000000000004, 34.24, 34.25, 34.26, 34.27, 34.28, 34.29,
      34.300000000000004, 34.31, 34.32, 34.33, 34.34, 34.35, 34.36, 34.37, 34.38,
      34.39, 34.4, 34.410000000000004, 34.42, 34.43, 34.44, 34.45, 34.46, 34.47,
      34.480000000000004, 34.49, 34.5, 34.51, 34.52, 34.53, 34.54,
      34.550000000000004, 34.56, 34.57, 34.58, 34.59, 34.6, 34.61, 34.62, 34.63,
      34.64, 34.65, 34.660000000000004, 34.67, 34.68, 34.69, 34.7, 34.71, 34.72,
      34.730000000000004, 34.74, 34.75, 34.76, 34.77, 34.78, 34.79,
      34.800000000000004, 34.81, 34.82, 34.83, 34.84, 34.85, 34.86, 34.87, 34.88,
      34.89, 34.9, 34.910000000000004, 34.92, 34.93, 34.94, 34.95, 34.96, 34.97,
      34.980000000000004, 34.99, 35.0, 35.01, 35.02, 35.03, 35.04,
      35.050000000000004, 35.06, 35.07, 35.08, 35.09, 35.1, 35.11, 35.12, 35.13,
      35.14, 35.15, 35.160000000000004, 35.17, 35.18, 35.19, 35.2, 35.21, 35.22,
      35.230000000000004, 35.24, 35.25, 35.26, 35.27, 35.28, 35.29,
      35.300000000000004, 35.31, 35.32, 35.33, 35.34, 35.35, 35.36, 35.37, 35.38,
      35.39, 35.4, 35.410000000000004, 35.42, 35.43, 35.44, 35.45, 35.46, 35.47,
      35.480000000000004, 35.49, 35.5, 35.51, 35.52, 35.53, 35.54,
      35.550000000000004, 35.56, 35.57, 35.58, 35.59, 35.6, 35.61, 35.62, 35.63,
      35.64, 35.65, 35.660000000000004, 35.67, 35.68, 35.69, 35.7, 35.71, 35.72,
      35.730000000000004, 35.74, 35.75, 35.76, 35.77, 35.78, 35.79,
      35.800000000000004, 35.81, 35.82, 35.83, 35.84, 35.85, 35.86, 35.87, 35.88,
      35.89, 35.9, 35.910000000000004, 35.92, 35.93, 35.94, 35.95, 35.96, 35.97,
      35.980000000000004, 35.99, 36.0, 36.01, 36.02, 36.03, 36.04,
      36.050000000000004, 36.06, 36.07, 36.08, 36.09, 36.1, 36.11, 36.12, 36.13,
      36.14, 36.15, 36.160000000000004, 36.17, 36.18, 36.19, 36.2, 36.21, 36.22,
      36.230000000000004, 36.24, 36.25, 36.26, 36.27, 36.28, 36.29,
      36.300000000000004, 36.31, 36.32, 36.33, 36.34, 36.35, 36.36, 36.37, 36.38,
      36.39, 36.4, 36.410000000000004, 36.42, 36.43, 36.44, 36.45, 36.46, 36.47,
      36.480000000000004, 36.49, 36.5, 36.51, 36.52, 36.53, 36.54,
      36.550000000000004, 36.56, 36.57, 36.58, 36.59, 36.6, 36.61, 36.62, 36.63,
      36.64, 36.65, 36.660000000000004, 36.67, 36.68, 36.69, 36.7, 36.71, 36.72,
      36.730000000000004, 36.74, 36.75, 36.76, 36.77, 36.78, 36.79,
      36.800000000000004, 36.81, 36.82, 36.83, 36.84, 36.85, 36.86, 36.87, 36.88,
      36.89, 36.9, 36.910000000000004, 36.92, 36.93, 36.94, 36.95, 36.96, 36.97,
      36.980000000000004, 36.99, 37.0, 37.01, 37.02, 37.03, 37.04,
      37.050000000000004, 37.06, 37.07, 37.08, 37.09, 37.1, 37.11, 37.12, 37.13,
      37.14, 37.15, 37.160000000000004, 37.17, 37.18, 37.19, 37.2, 37.21, 37.22,
      37.230000000000004, 37.24, 37.25, 37.26, 37.27, 37.28, 37.29,
      37.300000000000004, 37.31, 37.32, 37.33, 37.34, 37.35, 37.36, 37.37, 37.38,
      37.39, 37.4, 37.410000000000004, 37.42, 37.43, 37.44, 37.45, 37.46, 37.47,
      37.480000000000004, 37.49, 37.5, 37.51, 37.52, 37.53, 37.54,
      37.550000000000004, 37.56, 37.57, 37.58, 37.59, 37.6, 37.61, 37.62, 37.63,
      37.64, 37.65, 37.660000000000004, 37.67, 37.68, 37.69, 37.7, 37.71, 37.72,
      37.730000000000004, 37.74, 37.75, 37.76, 37.77, 37.78, 37.79,
      37.800000000000004, 37.81, 37.82, 37.83, 37.84, 37.85, 37.86, 37.87, 37.88,
      37.89, 37.9, 37.910000000000004, 37.92, 37.93, 37.94, 37.95, 37.96, 37.97,
      37.980000000000004, 37.99, 38.0, 38.01, 38.02, 38.03, 38.04,
      38.050000000000004, 38.06, 38.07, 38.08, 38.09, 38.1, 38.11, 38.12, 38.13,
      38.14, 38.15, 38.160000000000004, 38.17, 38.18, 38.19, 38.2, 38.21, 38.22,
      38.230000000000004, 38.24, 38.25, 38.26, 38.27, 38.28, 38.29,
      38.300000000000004, 38.31, 38.32, 38.33, 38.34, 38.35, 38.36, 38.37, 38.38,
      38.39, 38.4, 38.410000000000004, 38.42, 38.43, 38.44, 38.45, 38.46, 38.47,
      38.480000000000004, 38.49, 38.5, 38.51, 38.52, 38.53, 38.54,
      38.550000000000004, 38.56, 38.57, 38.58, 38.59, 38.6, 38.61, 38.62, 38.63,
      38.64, 38.65, 38.660000000000004, 38.67, 38.68, 38.69, 38.7, 38.71, 38.72,
      38.730000000000004, 38.74, 38.75, 38.76, 38.77, 38.78, 38.79,
      38.800000000000004, 38.81, 38.82, 38.83, 38.84, 38.85, 38.86, 38.87, 38.88,
      38.89, 38.9, 38.910000000000004, 38.92, 38.93, 38.94, 38.95, 38.96, 38.97,
      38.980000000000004, 38.99, 39.0, 39.01, 39.02, 39.03, 39.04,
      39.050000000000004, 39.06, 39.07, 39.08, 39.09, 39.1, 39.11, 39.12, 39.13,
      39.14, 39.15, 39.160000000000004, 39.17, 39.18, 39.19, 39.2, 39.21, 39.22,
      39.230000000000004, 39.24, 39.25, 39.26, 39.27, 39.28, 39.29,
      39.300000000000004, 39.31, 39.32, 39.33, 39.34, 39.35, 39.36, 39.37, 39.38,
      39.39, 39.4, 39.410000000000004, 39.42, 39.43, 39.44, 39.45, 39.46, 39.47,
      39.480000000000004, 39.49, 39.5, 39.51, 39.52, 39.53, 39.54,
      39.550000000000004, 39.56, 39.57, 39.58, 39.59, 39.6, 39.61, 39.62, 39.63,
      39.64, 39.65, 39.660000000000004, 39.67, 39.68, 39.69, 39.7, 39.71, 39.72,
      39.730000000000004, 39.74, 39.75, 39.76, 39.77, 39.78, 39.79,
      39.800000000000004, 39.81, 39.82, 39.83, 39.84, 39.85, 39.86, 39.87, 39.88,
      39.89, 39.9, 39.910000000000004, 39.92, 39.93, 39.94, 39.95, 39.96, 39.97,
      39.980000000000004, 39.99, 40.0, 40.01, 40.02, 40.03, 40.04,
      40.050000000000004, 40.06, 40.07, 40.08, 40.09, 40.1, 40.11, 40.12, 40.13,
      40.14, 40.15, 40.160000000000004, 40.17, 40.18, 40.19, 40.2, 40.21, 40.22,
      40.230000000000004, 40.24, 40.25, 40.26, 40.27, 40.28, 40.29,
      40.300000000000004, 40.31, 40.32, 40.33, 40.34, 40.35, 40.36, 40.37, 40.38,
      40.39, 40.4, 40.410000000000004, 40.42, 40.43, 40.44, 40.45, 40.46, 40.47,
      40.480000000000004, 40.49, 40.5, 40.51, 40.52, 40.53, 40.54,
      40.550000000000004, 40.56, 40.57, 40.58, 40.59, 40.6, 40.61, 40.62, 40.63,
      40.64, 40.65, 40.660000000000004, 40.67, 40.68, 40.69, 40.7, 40.71, 40.72,
      40.730000000000004, 40.74, 40.75, 40.76, 40.77, 40.78, 40.79,
      40.800000000000004, 40.81, 40.82, 40.83, 40.84, 40.85, 40.86, 40.87, 40.88,
      40.89, 40.9, 40.910000000000004, 40.92, 40.93, 40.94, 40.95, 40.96, 40.97,
      40.980000000000004, 40.99, 41.0, 41.01, 41.02, 41.03, 41.04,
      41.050000000000004, 41.06, 41.07, 41.08, 41.09, 41.1, 41.11, 41.12, 41.13,
      41.14, 41.15, 41.160000000000004, 41.17, 41.18, 41.19, 41.2, 41.21, 41.22,
      41.230000000000004, 41.24, 41.25, 41.26, 41.27, 41.28, 41.29,
      41.300000000000004, 41.31, 41.32, 41.33, 41.34, 41.35, 41.36, 41.37, 41.38,
      41.39, 41.4, 41.410000000000004, 41.42, 41.43, 41.44, 41.45, 41.46, 41.47,
      41.480000000000004, 41.49, 41.5, 41.51, 41.52, 41.53, 41.54,
      41.550000000000004, 41.56, 41.57, 41.58, 41.59, 41.6, 41.61, 41.62, 41.63,
      41.64, 41.65, 41.660000000000004, 41.67, 41.68, 41.69, 41.7, 41.71, 41.72,
      41.730000000000004, 41.74, 41.75, 41.76, 41.77, 41.78, 41.79,
      41.800000000000004, 41.81, 41.82, 41.83, 41.84, 41.85, 41.86, 41.87, 41.88,
      41.89, 41.9, 41.910000000000004, 41.92, 41.93, 41.94, 41.95, 41.96, 41.97,
      41.980000000000004, 41.99, 42.0, 42.01, 42.02, 42.03, 42.04,
      42.050000000000004, 42.06, 42.07, 42.08, 42.09, 42.1, 42.11, 42.12, 42.13,
      42.14, 42.15, 42.160000000000004, 42.17, 42.18, 42.19, 42.2, 42.21, 42.22,
      42.230000000000004, 42.24, 42.25, 42.26, 42.27, 42.28, 42.29,
      42.300000000000004, 42.31, 42.32, 42.33, 42.34, 42.35, 42.36, 42.37, 42.38,
      42.39, 42.4, 42.410000000000004, 42.42, 42.43, 42.44, 42.45, 42.46, 42.47,
      42.480000000000004, 42.49, 42.5, 42.51, 42.52, 42.53, 42.54,
      42.550000000000004, 42.56, 42.57, 42.58, 42.59, 42.6, 42.61, 42.62, 42.63,
      42.64, 42.65, 42.660000000000004, 42.67, 42.68, 42.69, 42.7, 42.71, 42.72,
      42.730000000000004, 42.74, 42.75, 42.76, 42.77, 42.78, 42.79,
      42.800000000000004, 42.81, 42.82, 42.83, 42.84, 42.85, 42.86, 42.87, 42.88,
      42.89, 42.9, 42.910000000000004, 42.92, 42.93, 42.94, 42.95, 42.96, 42.97,
      42.980000000000004, 42.99, 43.0, 43.01, 43.02, 43.03, 43.04,
      43.050000000000004, 43.06, 43.07, 43.08, 43.09, 43.1, 43.11, 43.12, 43.13,
      43.14, 43.15, 43.160000000000004, 43.17, 43.18, 43.19, 43.2, 43.21, 43.22,
      43.230000000000004, 43.24, 43.25, 43.26, 43.27, 43.28, 43.29,
      43.300000000000004, 43.31, 43.32, 43.33, 43.34, 43.35, 43.36, 43.37, 43.38,
      43.39, 43.4, 43.410000000000004, 43.42, 43.43, 43.44, 43.45, 43.46, 43.47,
      43.480000000000004, 43.49, 43.5, 43.51, 43.52, 43.53, 43.54,
      43.550000000000004, 43.56, 43.57, 43.58, 43.59, 43.6, 43.61, 43.62, 43.63,
      43.64, 43.65, 43.660000000000004, 43.67, 43.68, 43.69, 43.7, 43.71, 43.72,
      43.730000000000004, 43.74, 43.75, 43.76, 43.77, 43.78, 43.79,
      43.800000000000004, 43.81, 43.82, 43.83, 43.84, 43.85, 43.86, 43.87, 43.88,
      43.89, 43.9, 43.910000000000004, 43.92, 43.93, 43.94, 43.95, 43.96, 43.97,
      43.980000000000004, 43.99, 44.0, 44.01, 44.02, 44.03, 44.04,
      44.050000000000004, 44.06, 44.07, 44.08, 44.09, 44.1, 44.11, 44.12, 44.13,
      44.14, 44.15, 44.160000000000004, 44.17, 44.18, 44.19, 44.2, 44.21, 44.22,
      44.230000000000004, 44.24, 44.25, 44.26, 44.27, 44.28, 44.29,
      44.300000000000004, 44.31, 44.32, 44.33, 44.34, 44.35, 44.36, 44.37, 44.38,
      44.39, 44.4, 44.410000000000004, 44.42, 44.43, 44.44, 44.45, 44.46, 44.47,
      44.480000000000004, 44.49, 44.5, 44.51, 44.52, 44.53, 44.54,
      44.550000000000004, 44.56, 44.57, 44.58, 44.59, 44.6, 44.61, 44.62, 44.63,
      44.64, 44.65, 44.660000000000004, 44.67, 44.68, 44.69, 44.7, 44.71, 44.72,
      44.730000000000004, 44.74, 44.75, 44.76, 44.77, 44.78, 44.79,
      44.800000000000004, 44.81, 44.82, 44.83, 44.84, 44.85, 44.86, 44.87, 44.88,
      44.89, 44.9, 44.910000000000004, 44.92, 44.93, 44.94, 44.95, 44.96, 44.97,
      44.980000000000004, 44.99, 45.0, 45.01, 45.02, 45.03, 45.04,
      45.050000000000004, 45.06, 45.07, 45.08, 45.09, 45.1, 45.11, 45.12, 45.13,
      45.14, 45.15, 45.160000000000004, 45.17, 45.18, 45.19, 45.2, 45.21, 45.22,
      45.230000000000004, 45.24, 45.25, 45.26, 45.27, 45.28, 45.29,
      45.300000000000004, 45.31, 45.32, 45.33, 45.34, 45.35, 45.36, 45.37, 45.38,
      45.39, 45.4, 45.410000000000004, 45.42, 45.43, 45.44, 45.45, 45.46, 45.47,
      45.480000000000004, 45.49, 45.5, 45.51, 45.52, 45.53, 45.54,
      45.550000000000004, 45.56, 45.57, 45.58, 45.59, 45.6, 45.61, 45.62, 45.63,
      45.64, 45.65, 45.660000000000004, 45.67, 45.68, 45.69, 45.7, 45.71, 45.72,
      45.730000000000004, 45.74, 45.75, 45.76, 45.77, 45.78, 45.79,
      45.800000000000004, 45.81, 45.82, 45.83, 45.84, 45.85, 45.86, 45.87, 45.88,
      45.89, 45.9, 45.910000000000004, 45.92, 45.93, 45.94, 45.95, 45.96, 45.97,
      45.980000000000004, 45.99, 46.0, 46.01, 46.02, 46.03, 46.04,
      46.050000000000004, 46.06, 46.07, 46.08, 46.09, 46.1, 46.11, 46.12, 46.13,
      46.14, 46.15, 46.160000000000004, 46.17, 46.18, 46.19, 46.2, 46.21, 46.22,
      46.230000000000004, 46.24, 46.25, 46.26, 46.27, 46.28, 46.29,
      46.300000000000004, 46.31, 46.32, 46.33, 46.34, 46.35, 46.36, 46.37, 46.38,
      46.39, 46.4, 46.410000000000004, 46.42, 46.43, 46.44, 46.45, 46.46, 46.47,
      46.480000000000004, 46.49, 46.5, 46.51, 46.52, 46.53, 46.54,
      46.550000000000004, 46.56, 46.57, 46.58, 46.59, 46.6, 46.61, 46.62, 46.63,
      46.64, 46.65, 46.660000000000004, 46.67, 46.68, 46.69, 46.7, 46.71, 46.72,
      46.730000000000004, 46.74, 46.75, 46.76, 46.77, 46.78, 46.79,
      46.800000000000004, 46.81, 46.82, 46.83, 46.84, 46.85, 46.86, 46.87, 46.88,
      46.89, 46.9, 46.910000000000004, 46.92, 46.93, 46.94, 46.95, 46.96, 46.97,
      46.980000000000004, 46.99, 47.0, 47.01, 47.02, 47.03, 47.04,
      47.050000000000004, 47.06, 47.07, 47.08, 47.09, 47.1, 47.11, 47.12, 47.13,
      47.14, 47.15, 47.160000000000004, 47.17, 47.18, 47.19, 47.2, 47.21, 47.22,
      47.230000000000004, 47.24, 47.25, 47.26, 47.27, 47.28, 47.29,
      47.300000000000004, 47.31, 47.32, 47.33, 47.34, 47.35, 47.36, 47.37, 47.38,
      47.39, 47.4, 47.410000000000004, 47.42, 47.43, 47.44, 47.45, 47.46, 47.47,
      47.480000000000004, 47.49, 47.5, 47.51, 47.52, 47.53, 47.54,
      47.550000000000004, 47.56, 47.57, 47.58, 47.59, 47.6, 47.61, 47.62, 47.63,
      47.64, 47.65, 47.660000000000004, 47.67, 47.68, 47.69, 47.7, 47.71, 47.72,
      47.730000000000004, 47.74, 47.75, 47.76, 47.77, 47.78, 47.79,
      47.800000000000004, 47.81, 47.82, 47.83, 47.84, 47.85, 47.86,
      47.870000000000005, 47.88, 47.89, 47.9, 47.910000000000004, 47.92, 47.93,
      47.94, 47.95, 47.96, 47.97, 47.980000000000004, 47.99, 48.0, 48.01, 48.02,
      48.03, 48.04, 48.050000000000004, 48.06, 48.07, 48.08, 48.09, 48.1, 48.11,
      48.120000000000005, 48.13, 48.14, 48.15, 48.160000000000004, 48.17, 48.18,
      48.19, 48.2, 48.21, 48.22, 48.230000000000004, 48.24, 48.25, 48.26, 48.27,
      48.28, 48.29, 48.300000000000004, 48.31, 48.32, 48.33, 48.34, 48.35, 48.36,
      48.370000000000005, 48.38, 48.39, 48.4, 48.410000000000004, 48.42, 48.43,
      48.44, 48.45, 48.46, 48.47, 48.480000000000004, 48.49, 48.5, 48.51, 48.52,
      48.53, 48.54, 48.550000000000004, 48.56, 48.57, 48.58, 48.59, 48.6, 48.61,
      48.620000000000005, 48.63, 48.64, 48.65, 48.660000000000004, 48.67, 48.68,
      48.69, 48.7, 48.71, 48.72, 48.730000000000004, 48.74, 48.75, 48.76, 48.77,
      48.78, 48.79, 48.800000000000004, 48.81, 48.82, 48.83, 48.84, 48.85, 48.86,
      48.870000000000005, 48.88, 48.89, 48.9, 48.910000000000004, 48.92, 48.93,
      48.94, 48.95, 48.96, 48.97, 48.980000000000004, 48.99, 49.0, 49.01, 49.02,
      49.03, 49.04, 49.050000000000004, 49.06, 49.07, 49.08, 49.09, 49.1, 49.11,
      49.120000000000005, 49.13, 49.14, 49.15, 49.160000000000004, 49.17, 49.18,
      49.19, 49.2, 49.21, 49.22, 49.230000000000004, 49.24, 49.25, 49.26, 49.27,
      49.28, 49.29, 49.300000000000004, 49.31, 49.32, 49.33, 49.34, 49.35, 49.36,
      49.370000000000005, 49.38, 49.39, 49.4, 49.410000000000004, 49.42, 49.43,
      49.44, 49.45, 49.46, 49.47, 49.480000000000004, 49.49, 49.5, 49.51, 49.52,
      49.53, 49.54, 49.550000000000004, 49.56, 49.57, 49.58, 49.59, 49.6, 49.61,
      49.620000000000005, 49.63, 49.64, 49.65, 49.660000000000004, 49.67, 49.68,
      49.69, 49.7, 49.71, 49.72, 49.730000000000004, 49.74, 49.75, 49.76, 49.77,
      49.78, 49.79, 49.800000000000004, 49.81, 49.82, 49.83, 49.84, 49.85, 49.86,
      49.870000000000005, 49.88, 49.89, 49.9, 49.910000000000004, 49.92, 49.93,
      49.94, 49.95, 49.96, 49.97, 49.980000000000004, 49.99, 50.0, 50.01, 50.02,
      50.03, 50.04, 50.050000000000004, 50.06, 50.07, 50.08, 50.09, 50.1, 50.11,
      50.120000000000005, 50.13, 50.14, 50.15, 50.160000000000004, 50.17, 50.18,
      50.19, 50.2, 50.21, 50.22, 50.230000000000004, 50.24, 50.25, 50.26, 50.27,
      50.28, 50.29, 50.300000000000004, 50.31, 50.32, 50.33, 50.34, 50.35, 50.36,
      50.370000000000005, 50.38, 50.39, 50.4, 50.410000000000004, 50.42, 50.43,
      50.44, 50.45, 50.46, 50.47, 50.480000000000004, 50.49, 50.5, 50.51, 50.52,
      50.53, 50.54, 50.550000000000004, 50.56, 50.57, 50.58, 50.59, 50.6, 50.61,
      50.620000000000005, 50.63, 50.64, 50.65, 50.660000000000004, 50.67, 50.68,
      50.69, 50.7, 50.71, 50.72, 50.730000000000004, 50.74, 50.75, 50.76, 50.77,
      50.78, 50.79, 50.800000000000004, 50.81, 50.82, 50.83, 50.84, 50.85, 50.86,
      50.870000000000005, 50.88, 50.89, 50.9, 50.910000000000004, 50.92, 50.93,
      50.94, 50.95, 50.96, 50.97, 50.980000000000004, 50.99, 51.0, 51.01, 51.02,
      51.03, 51.04, 51.050000000000004, 51.06, 51.07, 51.08, 51.09, 51.1, 51.11,
      51.120000000000005, 51.13, 51.14, 51.15, 51.160000000000004, 51.17, 51.18,
      51.19, 51.2, 51.21, 51.22, 51.230000000000004, 51.24, 51.25, 51.26, 51.27,
      51.28, 51.29, 51.300000000000004, 51.31, 51.32, 51.33, 51.34, 51.35, 51.36,
      51.370000000000005, 51.38, 51.39, 51.4, 51.410000000000004, 51.42, 51.43,
      51.44, 51.45, 51.46, 51.47, 51.480000000000004, 51.49, 51.5, 51.51, 51.52,
      51.53, 51.54, 51.550000000000004, 51.56, 51.57, 51.58, 51.59, 51.6, 51.61,
      51.620000000000005, 51.63, 51.64, 51.65, 51.660000000000004, 51.67, 51.68,
      51.69, 51.7, 51.71, 51.72, 51.730000000000004, 51.74, 51.75, 51.76, 51.77,
      51.78, 51.79, 51.800000000000004, 51.81, 51.82, 51.83, 51.84, 51.85, 51.86,
      51.870000000000005, 51.88, 51.89, 51.9, 51.910000000000004, 51.92, 51.93,
      51.94, 51.95, 51.96, 51.97, 51.980000000000004, 51.99, 52.0, 52.01, 52.02,
      52.03, 52.04, 52.050000000000004, 52.06, 52.07, 52.08, 52.09, 52.1, 52.11,
      52.120000000000005, 52.13, 52.14, 52.15, 52.160000000000004, 52.17, 52.18,
      52.19, 52.2, 52.21, 52.22, 52.230000000000004, 52.24, 52.25, 52.26, 52.27,
      52.28, 52.29, 52.300000000000004, 52.31, 52.32, 52.33, 52.34, 52.35, 52.36,
      52.370000000000005, 52.38, 52.39, 52.4, 52.410000000000004, 52.42, 52.43,
      52.44, 52.45, 52.46, 52.47, 52.480000000000004, 52.49, 52.5, 52.51, 52.52,
      52.53, 52.54, 52.550000000000004, 52.56, 52.57, 52.58, 52.59, 52.6, 52.61,
      52.620000000000005, 52.63, 52.64, 52.65, 52.660000000000004, 52.67, 52.68,
      52.69, 52.7, 52.71, 52.72, 52.730000000000004, 52.74, 52.75, 52.76, 52.77,
      52.78, 52.79, 52.800000000000004, 52.81, 52.82, 52.83, 52.84, 52.85, 52.86,
      52.870000000000005, 52.88, 52.89, 52.9, 52.910000000000004, 52.92, 52.93,
      52.94, 52.95, 52.96, 52.97, 52.980000000000004, 52.99, 53.0, 53.01, 53.02,
      53.03, 53.04, 53.050000000000004, 53.06, 53.07, 53.08, 53.09, 53.1, 53.11,
      53.120000000000005, 53.13, 53.14, 53.15, 53.160000000000004, 53.17, 53.18,
      53.19, 53.2, 53.21, 53.22, 53.230000000000004, 53.24, 53.25, 53.26, 53.27,
      53.28, 53.29, 53.300000000000004, 53.31, 53.32, 53.33, 53.34, 53.35, 53.36,
      53.370000000000005, 53.38, 53.39, 53.4, 53.410000000000004, 53.42, 53.43,
      53.44, 53.45, 53.46, 53.47, 53.480000000000004, 53.49, 53.5, 53.51, 53.52,
      53.53, 53.54, 53.550000000000004, 53.56, 53.57, 53.58, 53.59, 53.6, 53.61,
      53.620000000000005, 53.63, 53.64, 53.65, 53.660000000000004, 53.67, 53.68,
      53.69, 53.7, 53.71, 53.72, 53.730000000000004, 53.74, 53.75, 53.76, 53.77,
      53.78, 53.79, 53.800000000000004, 53.81, 53.82, 53.83, 53.84, 53.85, 53.86,
      53.870000000000005, 53.88, 53.89, 53.9, 53.910000000000004, 53.92, 53.93,
      53.94, 53.95, 53.96, 53.97, 53.980000000000004, 53.99, 54.0, 54.01, 54.02,
      54.03, 54.04, 54.050000000000004, 54.06, 54.07, 54.08, 54.09, 54.1, 54.11,
      54.120000000000005, 54.13, 54.14, 54.15, 54.160000000000004, 54.17, 54.18,
      54.19, 54.2, 54.21, 54.22, 54.230000000000004, 54.24, 54.25, 54.26, 54.27,
      54.28, 54.29, 54.300000000000004, 54.31, 54.32, 54.33, 54.34, 54.35, 54.36,
      54.370000000000005, 54.38, 54.39, 54.4, 54.410000000000004, 54.42, 54.43,
      54.44, 54.45, 54.46, 54.47, 54.480000000000004, 54.49, 54.5, 54.51, 54.52,
      54.53, 54.54, 54.550000000000004, 54.56, 54.57, 54.58, 54.59, 54.6, 54.61,
      54.620000000000005, 54.63, 54.64, 54.65, 54.660000000000004, 54.67, 54.68,
      54.69, 54.7, 54.71, 54.72, 54.730000000000004, 54.74, 54.75, 54.76, 54.77,
      54.78, 54.79, 54.800000000000004, 54.81, 54.82, 54.83, 54.84, 54.85, 54.86,
      54.870000000000005, 54.88, 54.89, 54.9, 54.910000000000004, 54.92, 54.93,
      54.94, 54.95, 54.96, 54.97, 54.980000000000004, 54.99, 55.0, 55.01, 55.02,
      55.03, 55.04, 55.050000000000004, 55.06, 55.07, 55.08, 55.09, 55.1, 55.11,
      55.120000000000005, 55.13, 55.14, 55.15, 55.160000000000004, 55.17, 55.18,
      55.19, 55.2, 55.21, 55.22, 55.230000000000004, 55.24, 55.25, 55.26, 55.27,
      55.28, 55.29, 55.300000000000004, 55.31, 55.32, 55.33, 55.34, 55.35, 55.36,
      55.370000000000005, 55.38, 55.39, 55.4, 55.410000000000004, 55.42, 55.43,
      55.44, 55.45, 55.46, 55.47, 55.480000000000004, 55.49, 55.5, 55.51, 55.52,
      55.53, 55.54, 55.550000000000004, 55.56, 55.57, 55.58, 55.59, 55.6, 55.61,
      55.620000000000005, 55.63, 55.64, 55.65, 55.660000000000004, 55.67, 55.68,
      55.69, 55.7, 55.71, 55.72, 55.730000000000004, 55.74, 55.75, 55.76, 55.77,
      55.78, 55.79, 55.800000000000004, 55.81, 55.82, 55.83, 55.84, 55.85, 55.86,
      55.870000000000005, 55.88, 55.89, 55.9, 55.910000000000004, 55.92, 55.93,
      55.94, 55.95, 55.96, 55.97, 55.980000000000004, 55.99, 56.0, 56.01, 56.02,
      56.03, 56.04, 56.050000000000004, 56.06, 56.07, 56.08, 56.09, 56.1, 56.11,
      56.120000000000005, 56.13, 56.14, 56.15, 56.160000000000004, 56.17, 56.18,
      56.19, 56.2, 56.21, 56.22, 56.230000000000004, 56.24, 56.25, 56.26, 56.27,
      56.28, 56.29, 56.300000000000004, 56.31, 56.32, 56.33, 56.34, 56.35, 56.36,
      56.370000000000005, 56.38, 56.39, 56.4, 56.410000000000004, 56.42, 56.43,
      56.44, 56.45, 56.46, 56.47, 56.480000000000004, 56.49, 56.5, 56.51, 56.52,
      56.53, 56.54, 56.550000000000004, 56.56, 56.57, 56.58, 56.59, 56.6, 56.61,
      56.620000000000005, 56.63, 56.64, 56.65, 56.660000000000004, 56.67, 56.68,
      56.69, 56.7, 56.71, 56.72, 56.730000000000004, 56.74, 56.75, 56.76, 56.77,
      56.78, 56.79, 56.800000000000004, 56.81, 56.82, 56.83, 56.84, 56.85, 56.86,
      56.870000000000005, 56.88, 56.89, 56.9, 56.910000000000004, 56.92, 56.93,
      56.94, 56.95, 56.96, 56.97, 56.980000000000004, 56.99, 57.0, 57.01, 57.02,
      57.03, 57.04, 57.050000000000004, 57.06, 57.07, 57.08, 57.09, 57.1, 57.11,
      57.120000000000005, 57.13, 57.14, 57.15, 57.160000000000004, 57.17, 57.18,
      57.19, 57.2, 57.21, 57.22, 57.230000000000004, 57.24, 57.25, 57.26, 57.27,
      57.28, 57.29, 57.300000000000004, 57.31, 57.32, 57.33, 57.34, 57.35, 57.36,
      57.370000000000005, 57.38, 57.39, 57.4, 57.410000000000004, 57.42, 57.43,
      57.44, 57.45, 57.46, 57.47, 57.480000000000004, 57.49, 57.5, 57.51, 57.52,
      57.53, 57.54, 57.550000000000004, 57.56, 57.57, 57.58, 57.59, 57.6, 57.61,
      57.620000000000005, 57.63, 57.64, 57.65, 57.660000000000004, 57.67, 57.68,
      57.69, 57.7, 57.71, 57.72, 57.730000000000004, 57.74, 57.75, 57.76, 57.77,
      57.78, 57.79, 57.800000000000004, 57.81, 57.82, 57.83, 57.84, 57.85, 57.86,
      57.870000000000005, 57.88, 57.89, 57.9, 57.910000000000004, 57.92, 57.93,
      57.94, 57.95, 57.96, 57.97, 57.980000000000004, 57.99, 58.0, 58.01, 58.02,
      58.03, 58.04, 58.050000000000004, 58.06, 58.07, 58.08, 58.09, 58.1, 58.11,
      58.120000000000005, 58.13, 58.14, 58.15, 58.160000000000004, 58.17, 58.18,
      58.19, 58.2, 58.21, 58.22, 58.230000000000004, 58.24, 58.25, 58.26, 58.27,
      58.28, 58.29, 58.300000000000004, 58.31, 58.32, 58.33, 58.34, 58.35, 58.36,
      58.370000000000005, 58.38, 58.39, 58.4, 58.410000000000004, 58.42, 58.43,
      58.44, 58.45, 58.46, 58.47, 58.480000000000004, 58.49, 58.5, 58.51, 58.52,
      58.53, 58.54, 58.550000000000004, 58.56, 58.57, 58.58, 58.59, 58.6, 58.61,
      58.620000000000005, 58.63, 58.64, 58.65, 58.660000000000004, 58.67, 58.68,
      58.69, 58.7, 58.71, 58.72, 58.730000000000004, 58.74, 58.75, 58.76, 58.77,
      58.78, 58.79, 58.800000000000004, 58.81, 58.82, 58.83, 58.84, 58.85, 58.86,
      58.870000000000005, 58.88, 58.89, 58.9, 58.910000000000004, 58.92, 58.93,
      58.94, 58.95, 58.96, 58.97, 58.980000000000004, 58.99, 59.0, 59.01, 59.02,
      59.03, 59.04, 59.050000000000004, 59.06, 59.07, 59.08, 59.09, 59.1, 59.11,
      59.120000000000005, 59.13, 59.14, 59.15, 59.160000000000004, 59.17, 59.18,
      59.19, 59.2, 59.21, 59.22, 59.230000000000004, 59.24, 59.25, 59.26, 59.27,
      59.28, 59.29, 59.300000000000004, 59.31, 59.32, 59.33, 59.34, 59.35, 59.36,
      59.370000000000005, 59.38, 59.39, 59.4, 59.410000000000004, 59.42, 59.43,
      59.44, 59.45, 59.46, 59.47, 59.480000000000004, 59.49, 59.5, 59.51, 59.52,
      59.53, 59.54, 59.550000000000004, 59.56, 59.57, 59.58, 59.59, 59.6, 59.61,
      59.620000000000005, 59.63, 59.64, 59.65, 59.660000000000004, 59.67, 59.68,
      59.69, 59.7, 59.71, 59.72, 59.730000000000004, 59.74, 59.75, 59.76, 59.77,
      59.78, 59.79, 59.800000000000004, 59.81, 59.82, 59.83, 59.84, 59.85, 59.86,
      59.870000000000005, 59.88, 59.89, 59.9, 59.910000000000004, 59.92, 59.93,
      59.94, 59.95, 59.96, 59.97, 59.980000000000004, 59.99, 60.0, 60.01, 60.02,
      60.03, 60.04, 60.050000000000004, 60.06, 60.07, 60.08, 60.09, 60.1, 60.11,
      60.120000000000005, 60.13, 60.14, 60.15, 60.160000000000004, 60.17, 60.18,
      60.19, 60.2, 60.21, 60.22, 60.230000000000004, 60.24, 60.25, 60.26, 60.27,
      60.28, 60.29, 60.300000000000004, 60.31, 60.32, 60.33, 60.34, 60.35, 60.36,
      60.370000000000005, 60.38, 60.39, 60.4, 60.410000000000004, 60.42, 60.43,
      60.44, 60.45, 60.46, 60.47, 60.480000000000004, 60.49, 60.5, 60.51, 60.52,
      60.53, 60.54, 60.550000000000004, 60.56, 60.57, 60.58, 60.59, 60.6, 60.61,
      60.620000000000005, 60.63, 60.64, 60.65, 60.660000000000004, 60.67, 60.68,
      60.69, 60.7, 60.71, 60.72, 60.730000000000004, 60.74, 60.75, 60.76, 60.77,
      60.78, 60.79, 60.800000000000004, 60.81, 60.82, 60.83, 60.84, 60.85, 60.86,
      60.870000000000005, 60.88, 60.89, 60.9, 60.910000000000004, 60.92, 60.93,
      60.94, 60.95, 60.96, 60.97, 60.980000000000004, 60.99, 61.0, 61.01, 61.02,
      61.03, 61.04, 61.050000000000004, 61.06, 61.07, 61.08, 61.09, 61.1, 61.11,
      61.120000000000005, 61.13, 61.14, 61.15, 61.160000000000004, 61.17, 61.18,
      61.19, 61.2, 61.21, 61.22, 61.230000000000004, 61.24, 61.25, 61.26, 61.27,
      61.28, 61.29, 61.300000000000004, 61.31, 61.32, 61.33, 61.34, 61.35, 61.36,
      61.370000000000005, 61.38, 61.39, 61.4, 61.410000000000004, 61.42, 61.43,
      61.44, 61.45, 61.46, 61.47, 61.480000000000004, 61.49, 61.5, 61.51, 61.52,
      61.53, 61.54, 61.550000000000004, 61.56, 61.57, 61.58, 61.59, 61.6, 61.61,
      61.620000000000005, 61.63, 61.64, 61.65, 61.660000000000004, 61.67, 61.68,
      61.690000000000005, 61.7, 61.71, 61.72, 61.730000000000004, 61.74, 61.75,
      61.76, 61.77, 61.78, 61.79, 61.800000000000004, 61.81, 61.82, 61.83, 61.84,
      61.85, 61.86, 61.870000000000005, 61.88, 61.89, 61.9, 61.910000000000004,
      61.92, 61.93, 61.940000000000005, 61.95, 61.96, 61.97, 61.980000000000004,
      61.99, 62.0, 62.01, 62.02, 62.03, 62.04, 62.050000000000004, 62.06, 62.07,
      62.08, 62.09, 62.1, 62.11, 62.120000000000005, 62.13, 62.14, 62.15,
      62.160000000000004, 62.17, 62.18, 62.190000000000005, 62.2, 62.21, 62.22,
      62.230000000000004, 62.24, 62.25, 62.26, 62.27, 62.28, 62.29,
      62.300000000000004, 62.31, 62.32, 62.33, 62.34, 62.35, 62.36,
      62.370000000000005, 62.38, 62.39, 62.4, 62.410000000000004, 62.42, 62.43,
      62.440000000000005, 62.45, 62.46, 62.47, 62.480000000000004, 62.49, 62.5,
      62.51, 62.52, 62.53, 62.54, 62.550000000000004, 62.56, 62.57, 62.58, 62.59,
      62.6, 62.61, 62.620000000000005, 62.63, 62.64, 62.65, 62.660000000000004,
      62.67, 62.68, 62.690000000000005, 62.7, 62.71, 62.72, 62.730000000000004,
      62.74, 62.75, 62.76, 62.77, 62.78, 62.79, 62.800000000000004, 62.81, 62.82,
      62.83, 62.84, 62.85, 62.86, 62.870000000000005, 62.88, 62.89, 62.9,
      62.910000000000004, 62.92, 62.93, 62.940000000000005, 62.95, 62.96, 62.97,
      62.980000000000004, 62.99, 63.0, 63.01, 63.02, 63.03, 63.04,
      63.050000000000004, 63.06, 63.07, 63.08, 63.09, 63.1, 63.11,
      63.120000000000005, 63.13, 63.14, 63.15, 63.160000000000004, 63.17, 63.18,
      63.190000000000005, 63.2, 63.21, 63.22, 63.230000000000004, 63.24, 63.25,
      63.26, 63.27, 63.28, 63.29, 63.300000000000004, 63.31, 63.32, 63.33, 63.34,
      63.35, 63.36, 63.370000000000005, 63.38, 63.39, 63.4, 63.410000000000004,
      63.42, 63.43, 63.440000000000005, 63.45, 63.46, 63.47, 63.480000000000004,
      63.49, 63.5, 63.51, 63.52, 63.53, 63.54, 63.550000000000004, 63.56, 63.57,
      63.58, 63.59, 63.6, 63.61, 63.620000000000005, 63.63, 63.64, 63.65,
      63.660000000000004, 63.67, 63.68, 63.690000000000005, 63.7, 63.71, 63.72,
      63.730000000000004, 63.74, 63.75, 63.76, 63.77, 63.78, 63.79,
      63.800000000000004, 63.81, 63.82, 63.83, 63.84, 63.85, 63.86,
      63.870000000000005, 63.88, 63.89, 63.9, 63.910000000000004, 63.92, 63.93,
      63.940000000000005, 63.95, 63.96, 63.97, 63.980000000000004, 63.99, 64.0,
      64.01, 64.02, 64.03, 64.04, 64.05, 64.06, 64.070000000000007, 64.08, 64.09,
      64.1, 64.11, 64.12, 64.13, 64.14, 64.15, 64.16, 64.17, 64.18, 64.19, 64.2,
      64.210000000000008, 64.22, 64.23, 64.24, 64.25, 64.26, 64.27, 64.28, 64.29,
      64.3, 64.31, 64.320000000000007, 64.33, 64.34, 64.35, 64.36, 64.37, 64.38,
      64.39, 64.4, 64.41, 64.42, 64.43, 64.44, 64.45, 64.460000000000008, 64.47,
      64.48, 64.49, 64.5, 64.51, 64.52, 64.53, 64.54, 64.55, 64.56,
      64.570000000000007, 64.58, 64.59, 64.6, 64.61, 64.62, 64.63, 64.64, 64.65,
      64.66, 64.67, 64.68, 64.69, 64.7, 64.710000000000008, 64.72, 64.73, 64.74,
      64.75, 64.76, 64.77, 64.78, 64.79, 64.8, 64.81, 64.820000000000007, 64.83,
      64.84, 64.85, 64.86, 64.87, 64.88, 64.89, 64.9, 64.91, 64.92, 64.93, 64.94,
      64.95, 64.960000000000008, 64.97, 64.98, 64.99, 65.0, 65.01, 65.02, 65.03,
      65.04, 65.05, 65.06, 65.070000000000007, 65.08, 65.09, 65.1, 65.11, 65.12,
      65.13, 65.14, 65.15, 65.16, 65.17, 65.18, 65.19, 65.2, 65.210000000000008,
      65.22, 65.23, 65.24, 65.25, 65.26, 65.27, 65.28, 65.29, 65.3, 65.31,
      65.320000000000007, 65.33, 65.34, 65.35, 65.36, 65.37, 65.38, 65.39, 65.4,
      65.41, 65.42, 65.43, 65.44, 65.45, 65.460000000000008, 65.47, 65.48, 65.49,
      65.5, 65.51, 65.52, 65.53, 65.54, 65.55, 65.56, 65.570000000000007, 65.58,
      65.59, 65.6, 65.61, 65.62, 65.63, 65.64, 65.65, 65.66, 65.67, 65.68, 65.69,
      65.7, 65.710000000000008, 65.72, 65.73, 65.74, 65.75, 65.76, 65.77, 65.78,
      65.79, 65.8, 65.81, 65.820000000000007, 65.83, 65.84, 65.85, 65.86, 65.87,
      65.88, 65.89, 65.9, 65.91, 65.92, 65.93, 65.94, 65.95, 65.960000000000008,
      65.97, 65.98, 65.99, 66.0, 66.01, 66.02, 66.03, 66.04, 66.05, 66.06,
      66.070000000000007, 66.08, 66.09, 66.1, 66.11, 66.12, 66.13, 66.14, 66.15,
      66.16, 66.17, 66.18, 66.19, 66.2, 66.210000000000008, 66.22, 66.23, 66.24,
      66.25, 66.26, 66.27, 66.28, 66.29, 66.3, 66.31, 66.320000000000007, 66.33,
      66.34, 66.35, 66.36, 66.37, 66.38, 66.39, 66.4, 66.41, 66.42, 66.43, 66.44,
      66.45, 66.460000000000008, 66.47, 66.48, 66.49, 66.5, 66.51, 66.52, 66.53,
      66.54, 66.55, 66.56, 66.570000000000007, 66.58, 66.59, 66.6, 66.61, 66.62,
      66.63, 66.64, 66.65, 66.66, 66.67, 66.68, 66.69, 66.7, 66.710000000000008,
      66.72, 66.73, 66.74, 66.75, 66.76, 66.77, 66.78, 66.79, 66.8, 66.81,
      66.820000000000007, 66.83, 66.84, 66.85, 66.86, 66.87, 66.88, 66.89, 66.9,
      66.91, 66.92, 66.93, 66.94, 66.95, 66.960000000000008, 66.97, 66.98, 66.99,
      67.0, 67.01, 67.02, 67.03, 67.04, 67.05, 67.06, 67.070000000000007, 67.08,
      67.09, 67.1, 67.11, 67.12, 67.13, 67.14, 67.15, 67.16, 67.17, 67.18, 67.19,
      67.2, 67.210000000000008, 67.22, 67.23, 67.24, 67.25, 67.26, 67.27, 67.28,
      67.29, 67.3, 67.31, 67.320000000000007, 67.33, 67.34, 67.35, 67.36, 67.37,
      67.38, 67.39, 67.4, 67.41, 67.42, 67.43, 67.44, 67.45, 67.460000000000008,
      67.47, 67.48, 67.49, 67.5, 67.51, 67.52, 67.53, 67.54, 67.55, 67.56,
      67.570000000000007, 67.58, 67.59, 67.6, 67.61, 67.62, 67.63, 67.64, 67.65,
      67.66, 67.67, 67.68, 67.69, 67.7, 67.710000000000008, 67.72, 67.73, 67.74,
      67.75, 67.76, 67.77, 67.78, 67.79, 67.8, 67.81, 67.820000000000007, 67.83,
      67.84, 67.85, 67.86, 67.87, 67.88, 67.89, 67.9, 67.91, 67.92, 67.93, 67.94,
      67.95, 67.960000000000008, 67.97, 67.98, 67.99, 68.0, 68.01, 68.02, 68.03,
      68.04, 68.05, 68.06, 68.070000000000007, 68.08, 68.09, 68.1, 68.11, 68.12,
      68.13, 68.14, 68.15, 68.16, 68.17, 68.18, 68.19, 68.2, 68.210000000000008,
      68.22, 68.23, 68.24, 68.25, 68.26, 68.27, 68.28, 68.29, 68.3, 68.31,
      68.320000000000007, 68.33, 68.34, 68.350000000000009, 68.36, 68.37, 68.38,
      68.39, 68.4, 68.41, 68.42, 68.43, 68.44, 68.45, 68.460000000000008, 68.47,
      68.48, 68.49, 68.5, 68.51, 68.52, 68.53, 68.54, 68.55, 68.56,
      68.570000000000007, 68.58, 68.59, 68.600000000000009, 68.61, 68.62, 68.63,
      68.64, 68.65, 68.66, 68.67, 68.68, 68.69, 68.7, 68.710000000000008, 68.72,
      68.73, 68.74, 68.75, 68.76, 68.77, 68.78, 68.79, 68.8, 68.81,
      68.820000000000007, 68.83, 68.84, 68.850000000000009, 68.86, 68.87, 68.88,
      68.89, 68.9, 68.91, 68.92, 68.93, 68.94, 68.95, 68.960000000000008, 68.97,
      68.98, 68.99, 69.0, 69.01, 69.02, 69.03, 69.04, 69.05, 69.06,
      69.070000000000007, 69.08, 69.09, 69.100000000000009, 69.11, 69.12, 69.13,
      69.14, 69.15, 69.16, 69.17, 69.18, 69.19, 69.2, 69.210000000000008, 69.22,
      69.23, 69.24, 69.25, 69.26, 69.27, 69.28, 69.29, 69.3, 69.31,
      69.320000000000007, 69.33, 69.34, 69.350000000000009, 69.36, 69.37, 69.38,
      69.39, 69.4, 69.41, 69.42, 69.43, 69.44, 69.45, 69.460000000000008, 69.47,
      69.48, 69.49, 69.5, 69.51, 69.52, 69.53, 69.54, 69.55, 69.56,
      69.570000000000007, 69.58, 69.59, 69.600000000000009, 69.61, 69.62, 69.63,
      69.64, 69.65, 69.66, 69.67, 69.68, 69.69, 69.7, 69.710000000000008, 69.72,
      69.73, 69.74, 69.75, 69.76, 69.77, 69.78, 69.79, 69.8, 69.81,
      69.820000000000007, 69.83, 69.84, 69.850000000000009, 69.86, 69.87, 69.88,
      69.89, 69.9, 69.91, 69.92, 69.93, 69.94, 69.95, 69.960000000000008, 69.97,
      69.98, 69.99, 70.0, 70.01, 70.02, 70.03, 70.04, 70.05, 70.06,
      70.070000000000007, 70.08, 70.09, 70.100000000000009, 70.11, 70.12, 70.13,
      70.14, 70.15, 70.16, 70.17, 70.18, 70.19, 70.2, 70.210000000000008, 70.22,
      70.23, 70.24, 70.25, 70.26, 70.27, 70.28, 70.29, 70.3, 70.31,
      70.320000000000007, 70.33, 70.34, 70.350000000000009, 70.36, 70.37, 70.38,
      70.39, 70.4, 70.41, 70.42, 70.43, 70.44, 70.45, 70.460000000000008, 70.47,
      70.48, 70.49, 70.5, 70.51, 70.52, 70.53, 70.54, 70.55, 70.56,
      70.570000000000007, 70.58, 70.59, 70.600000000000009, 70.61, 70.62, 70.63,
      70.64, 70.65, 70.66, 70.67, 70.68, 70.69, 70.7, 70.710000000000008, 70.72,
      70.73, 70.74, 70.75, 70.76, 70.77, 70.78, 70.79, 70.8, 70.81,
      70.820000000000007, 70.83, 70.84, 70.850000000000009, 70.86, 70.87, 70.88,
      70.89, 70.9, 70.91, 70.92, 70.93, 70.94, 70.95, 70.960000000000008, 70.97,
      70.98, 70.99, 71.0, 71.01, 71.02, 71.03, 71.04, 71.05, 71.06,
      71.070000000000007, 71.08, 71.09, 71.100000000000009, 71.11, 71.12, 71.13,
      71.14, 71.15, 71.16, 71.17, 71.18, 71.19, 71.2, 71.210000000000008, 71.22,
      71.23, 71.24, 71.25, 71.26, 71.27, 71.28, 71.29, 71.3, 71.31,
      71.320000000000007, 71.33, 71.34, 71.350000000000009, 71.36, 71.37, 71.38,
      71.39, 71.4, 71.41, 71.42, 71.43, 71.44, 71.45, 71.460000000000008, 71.47,
      71.48, 71.49, 71.5, 71.51, 71.52, 71.53, 71.54, 71.55, 71.56,
      71.570000000000007, 71.58, 71.59, 71.600000000000009, 71.61, 71.62, 71.63,
      71.64, 71.65, 71.66, 71.67, 71.68, 71.69, 71.7, 71.710000000000008, 71.72,
      71.73, 71.74, 71.75, 71.76, 71.77, 71.78, 71.79, 71.8, 71.81,
      71.820000000000007, 71.83, 71.84, 71.850000000000009, 71.86, 71.87, 71.88,
      71.89, 71.9, 71.91, 71.92, 71.93, 71.94, 71.95, 71.960000000000008, 71.97,
      71.98, 71.99, 72.0, 72.01, 72.02, 72.03, 72.04, 72.05, 72.06,
      72.070000000000007, 72.08, 72.09, 72.100000000000009, 72.11, 72.12, 72.13,
      72.14, 72.15, 72.16, 72.17, 72.18, 72.19, 72.2, 72.210000000000008, 72.22,
      72.23, 72.24, 72.25, 72.26, 72.27, 72.28, 72.29, 72.3, 72.31,
      72.320000000000007, 72.33, 72.34, 72.350000000000009, 72.36, 72.37, 72.38,
      72.39, 72.4, 72.41, 72.42, 72.43, 72.44, 72.45, 72.460000000000008, 72.47,
      72.48, 72.49, 72.5, 72.51, 72.52, 72.53, 72.54, 72.55, 72.56,
      72.570000000000007, 72.58, 72.59, 72.600000000000009, 72.61, 72.62, 72.63,
      72.64, 72.65, 72.66, 72.67, 72.68, 72.69, 72.7, 72.710000000000008, 72.72,
      72.73, 72.74, 72.75, 72.76, 72.77, 72.78, 72.79, 72.8, 72.81,
      72.820000000000007, 72.83, 72.84, 72.850000000000009, 72.86, 72.87, 72.88,
      72.89, 72.9, 72.91, 72.92, 72.93, 72.94, 72.95, 72.960000000000008, 72.97,
      72.98, 72.99, 73.0, 73.01, 73.02, 73.03, 73.04, 73.05, 73.06,
      73.070000000000007, 73.08, 73.09, 73.100000000000009, 73.11, 73.12, 73.13,
      73.14, 73.15, 73.16, 73.17, 73.18, 73.19, 73.2, 73.210000000000008, 73.22,
      73.23, 73.24, 73.25, 73.26, 73.27, 73.28, 73.29, 73.3, 73.31,
      73.320000000000007, 73.33, 73.34, 73.350000000000009, 73.36, 73.37, 73.38,
      73.39, 73.4, 73.41, 73.42, 73.43, 73.44, 73.45, 73.460000000000008, 73.47,
      73.48, 73.49, 73.5, 73.51, 73.52, 73.53, 73.54, 73.55, 73.56,
      73.570000000000007, 73.58, 73.59, 73.600000000000009, 73.61, 73.62, 73.63,
      73.64, 73.65, 73.66, 73.67, 73.68, 73.69, 73.7, 73.710000000000008, 73.72,
      73.73, 73.74, 73.75, 73.76, 73.77, 73.78, 73.79, 73.8, 73.81,
      73.820000000000007, 73.83, 73.84, 73.850000000000009, 73.86, 73.87, 73.88,
      73.89, 73.9, 73.91, 73.92, 73.93, 73.94, 73.95, 73.960000000000008, 73.97,
      73.98, 73.99, 74.0, 74.01, 74.02, 74.03, 74.04, 74.05, 74.06,
      74.070000000000007, 74.08, 74.09, 74.100000000000009, 74.11, 74.12, 74.13,
      74.14, 74.15, 74.16, 74.17, 74.18, 74.19, 74.2, 74.210000000000008, 74.22,
      74.23, 74.24, 74.25, 74.26, 74.27, 74.28, 74.29, 74.3, 74.31,
      74.320000000000007, 74.33, 74.34, 74.350000000000009, 74.36, 74.37, 74.38,
      74.39, 74.4, 74.41, 74.42, 74.43, 74.44, 74.45, 74.460000000000008, 74.47,
      74.48, 74.49, 74.5, 74.51, 74.52, 74.53, 74.54, 74.55, 74.56,
      74.570000000000007, 74.58, 74.59, 74.600000000000009, 74.61, 74.62, 74.63,
      74.64, 74.65, 74.66, 74.67, 74.68, 74.69, 74.7, 74.710000000000008, 74.72,
      74.73, 74.74, 74.75, 74.76, 74.77, 74.78, 74.79, 74.8, 74.81,
      74.820000000000007, 74.83, 74.84, 74.850000000000009, 74.86, 74.87, 74.88,
      74.89, 74.9, 74.91, 74.92, 74.93, 74.94, 74.95, 74.960000000000008, 74.97,
      74.98, 74.99, 75.0, 75.01, 75.02, 75.03, 75.04, 75.05, 75.06,
      75.070000000000007, 75.08, 75.09, 75.100000000000009, 75.11, 75.12, 75.13,
      75.14, 75.15, 75.16, 75.17, 75.18, 75.19, 75.2, 75.210000000000008, 75.22,
      75.23, 75.24, 75.25, 75.26, 75.27, 75.28, 75.29, 75.3, 75.31,
      75.320000000000007, 75.33, 75.34, 75.350000000000009, 75.36, 75.37, 75.38,
      75.39, 75.4, 75.41, 75.42, 75.43, 75.44, 75.45, 75.460000000000008, 75.47,
      75.48, 75.49, 75.5, 75.51, 75.52, 75.53, 75.54, 75.55, 75.56,
      75.570000000000007, 75.58, 75.59, 75.600000000000009, 75.61, 75.62, 75.63,
      75.64, 75.65, 75.66, 75.67, 75.68, 75.69, 75.7, 75.710000000000008, 75.72,
      75.73, 75.74, 75.75, 75.76, 75.77, 75.78, 75.79, 75.8, 75.81,
      75.820000000000007, 75.83, 75.84, 75.850000000000009, 75.86, 75.87, 75.88,
      75.89, 75.9, 75.91, 75.92, 75.93, 75.94, 75.95, 75.960000000000008, 75.97,
      75.98, 75.99, 76.0, 76.01, 76.02, 76.03, 76.04, 76.05, 76.06,
      76.070000000000007, 76.08, 76.09, 76.100000000000009, 76.11, 76.12, 76.13,
      76.14, 76.15, 76.16, 76.17, 76.18, 76.19, 76.2, 76.210000000000008, 76.22,
      76.23, 76.24, 76.25, 76.26, 76.27, 76.28, 76.29, 76.3, 76.31,
      76.320000000000007, 76.33, 76.34, 76.350000000000009, 76.36, 76.37, 76.38,
      76.39, 76.4, 76.41, 76.42, 76.43, 76.44, 76.45, 76.460000000000008, 76.47,
      76.48, 76.49, 76.5, 76.51, 76.52, 76.53, 76.54, 76.55, 76.56,
      76.570000000000007, 76.58, 76.59, 76.600000000000009, 76.61, 76.62, 76.63,
      76.64, 76.65, 76.66, 76.67, 76.68, 76.69, 76.7, 76.710000000000008, 76.72,
      76.73, 76.74, 76.75, 76.76, 76.77, 76.78, 76.79, 76.8, 76.81,
      76.820000000000007, 76.83, 76.84, 76.850000000000009, 76.86, 76.87, 76.88,
      76.89, 76.9, 76.91, 76.92, 76.93, 76.94, 76.95, 76.960000000000008, 76.97,
      76.98, 76.99, 77.0, 77.01, 77.02, 77.03, 77.04, 77.05, 77.06,
      77.070000000000007, 77.08, 77.09, 77.100000000000009, 77.11, 77.12, 77.13,
      77.14, 77.15, 77.16, 77.17, 77.18, 77.19, 77.2, 77.210000000000008, 77.22,
      77.23, 77.24, 77.25, 77.26, 77.27, 77.28, 77.29, 77.3, 77.31,
      77.320000000000007, 77.33, 77.34, 77.350000000000009, 77.36, 77.37, 77.38,
      77.39, 77.4, 77.41, 77.42, 77.43, 77.44, 77.45, 77.460000000000008, 77.47,
      77.48, 77.49, 77.5, 77.51, 77.52, 77.53, 77.54, 77.55, 77.56,
      77.570000000000007, 77.58, 77.59, 77.600000000000009, 77.61, 77.62, 77.63,
      77.64, 77.65, 77.66, 77.67, 77.68, 77.69, 77.7, 77.710000000000008, 77.72,
      77.73, 77.74, 77.75, 77.76, 77.77, 77.78, 77.79, 77.8, 77.81,
      77.820000000000007, 77.83, 77.84, 77.850000000000009, 77.86, 77.87, 77.88,
      77.89, 77.9, 77.91, 77.92, 77.93, 77.94, 77.95, 77.960000000000008, 77.97,
      77.98, 77.99, 78.0, 78.01, 78.02, 78.03, 78.04, 78.05, 78.06,
      78.070000000000007, 78.08, 78.09, 78.100000000000009, 78.11, 78.12, 78.13,
      78.14, 78.15, 78.16, 78.17, 78.18, 78.19, 78.2, 78.210000000000008, 78.22,
      78.23, 78.24, 78.25, 78.26, 78.27, 78.28, 78.29, 78.3, 78.31,
      78.320000000000007, 78.33, 78.34, 78.350000000000009, 78.36, 78.37, 78.38,
      78.39, 78.4, 78.41, 78.42, 78.43, 78.44, 78.45, 78.460000000000008, 78.47,
      78.48, 78.49, 78.5, 78.51, 78.52, 78.53, 78.54, 78.55, 78.56,
      78.570000000000007, 78.58, 78.59, 78.600000000000009, 78.61, 78.62, 78.63,
      78.64, 78.65, 78.66, 78.67, 78.68, 78.69, 78.7, 78.710000000000008, 78.72,
      78.73, 78.74, 78.75, 78.76, 78.77, 78.78, 78.79, 78.8, 78.81,
      78.820000000000007, 78.83, 78.84, 78.850000000000009, 78.86, 78.87, 78.88,
      78.89, 78.9, 78.91, 78.92, 78.93, 78.94, 78.95, 78.960000000000008, 78.97,
      78.98, 78.99, 79.0, 79.01, 79.02, 79.03, 79.04, 79.05, 79.06,
      79.070000000000007, 79.08, 79.09, 79.100000000000009, 79.11, 79.12, 79.13,
      79.14, 79.15, 79.16, 79.17, 79.18, 79.19, 79.2, 79.210000000000008, 79.22,
      79.23, 79.24, 79.25, 79.26, 79.27, 79.28, 79.29, 79.3, 79.31,
      79.320000000000007, 79.33, 79.34, 79.350000000000009, 79.36, 79.37, 79.38,
      79.39, 79.4, 79.41, 79.42, 79.43, 79.44, 79.45, 79.460000000000008, 79.47,
      79.48, 79.49, 79.5, 79.51, 79.52, 79.53, 79.54, 79.55, 79.56,
      79.570000000000007, 79.58, 79.59, 79.600000000000009, 79.61, 79.62, 79.63,
      79.64, 79.65, 79.66, 79.67, 79.68, 79.69, 79.7, 79.710000000000008, 79.72,
      79.73, 79.74, 79.75, 79.76, 79.77, 79.78, 79.79, 79.8, 79.81,
      79.820000000000007, 79.83, 79.84, 79.850000000000009, 79.86, 79.87, 79.88,
      79.89, 79.9, 79.91, 79.92, 79.93, 79.94, 79.95, 79.960000000000008, 79.97,
      79.98, 79.99 } ;

    static real_T pDataValues0[] = { -1000.0, -1000.0, -1000.0, -1000.0, -1000.0,
      -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0,
      -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0,
      -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0,
      -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0,
      -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0,
      -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0,
      -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0,
      -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0,
      -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0,
      -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0,
      -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0,
      -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    ModelWithControllersOnly_DW.FromWorkspace_PWORK.TimePtr = static_cast<void *>
      (pTimeValues0);
    ModelWithControllersOnly_DW.FromWorkspace_PWORK.DataPtr = static_cast<void *>
      (pDataValues0);
    ModelWithControllersOnly_DW.FromWorkspace_IWORK.PrevIndex = 0;
  }

  /* Start for FromWorkspace: '<S5>/From Workspace1' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06,
      0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18,
      0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31,
      0.32, 0.33, 0.34, 0.35000000000000003, 0.36, 0.37, 0.38, 0.39, 0.4,
      0.41000000000000003, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47000000000000003,
      0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57000000000000006,
      0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68,
      0.69000000000000006, 0.70000000000000007, 0.71, 0.72, 0.73, 0.74, 0.75,
      0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82000000000000006,
      0.83000000000000007, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92,
      0.93, 0.94000000000000006, 0.95000000000000007, 0.96, 0.97, 0.98, 0.99,
      1.0, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.12,
      1.1300000000000001, 1.1400000000000001, 1.1500000000000001, 1.16, 1.17,
      1.18, 1.19, 1.2, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3,
      1.31, 1.32, 1.33, 1.34, 1.35, 1.36, 1.37, 1.3800000000000001,
      1.3900000000000001, 1.4000000000000001, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46,
      1.47, 1.48, 1.49, 1.5, 1.51, 1.52, 1.53, 1.54, 1.55, 1.56, 1.57, 1.58,
      1.59, 1.6, 1.61, 1.62, 1.6300000000000001, 1.6400000000000001,
      1.6500000000000001, 1.6600000000000001, 1.67, 1.68, 1.69, 1.7, 1.71, 1.72,
      1.73, 1.74, 1.75, 1.76, 1.77, 1.78, 1.79, 1.8, 1.81, 1.82, 1.83, 1.84,
      1.85, 1.86, 1.87, 1.8800000000000001, 1.8900000000000001,
      1.9000000000000001, 1.9100000000000001, 1.92, 1.93, 1.94, 1.95, 1.96, 1.97,
      1.98, 1.99, 2.0, 2.0100000000000002, 2.02, 2.0300000000000002, 2.04, 2.05,
      2.06, 2.07, 2.08, 2.09, 2.1, 2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17,
      2.18, 2.19, 2.2, 2.21, 2.22, 2.23, 2.24, 2.25, 2.2600000000000002, 2.27,
      2.2800000000000002, 2.29, 2.3000000000000003, 2.31, 2.32, 2.33, 2.34, 2.35,
      2.36, 2.37, 2.38, 2.39, 2.4, 2.41, 2.42, 2.43, 2.44, 2.45, 2.46, 2.47,
      2.48, 2.49, 2.5, 2.5100000000000002, 2.52, 2.5300000000000002, 2.54,
      2.5500000000000003, 2.56, 2.57, 2.58, 2.59, 2.6, 2.61, 2.62, 2.63, 2.64,
      2.65, 2.66, 2.67, 2.68, 2.69, 2.7, 2.71, 2.72, 2.73, 2.74, 2.75,
      2.7600000000000002, 2.77, 2.7800000000000002, 2.79, 2.8000000000000003,
      2.81, 2.82, 2.83, 2.84, 2.85, 2.86, 2.87, 2.88, 2.89, 2.9, 2.91, 2.92,
      2.93, 2.94, 2.95, 2.96, 2.97, 2.98, 2.99, 3.0, 3.0100000000000002, 3.02,
      3.0300000000000002, 3.04, 3.0500000000000003, 3.06, 3.0700000000000003,
      3.08, 3.09, 3.1, 3.11, 3.12, 3.13, 3.14, 3.15, 3.16, 3.17, 3.18, 3.19, 3.2,
      3.21, 3.22, 3.23, 3.24, 3.25, 3.2600000000000002, 3.27, 3.2800000000000002,
      3.29, 3.3000000000000003, 3.31, 3.3200000000000003, 3.33, 3.34, 3.35, 3.36,
      3.37, 3.38, 3.39, 3.4, 3.41, 3.42, 3.43, 3.44, 3.45, 3.46, 3.47, 3.48,
      3.49, 3.5, 3.5100000000000002, 3.52, 3.5300000000000002, 3.54,
      3.5500000000000003, 3.56, 3.5700000000000003, 3.58, 3.59, 3.6, 3.61, 3.62,
      3.63, 3.64, 3.65, 3.66, 3.67, 3.68, 3.69, 3.7, 3.71, 3.72, 3.73, 3.74,
      3.75, 3.7600000000000002, 3.77, 3.7800000000000002, 3.79,
      3.8000000000000003, 3.81, 3.8200000000000003, 3.83, 3.84, 3.85, 3.86, 3.87,
      3.88, 3.89, 3.9, 3.91, 3.92, 3.93, 3.94, 3.95, 3.96, 3.97, 3.98, 3.99, 4.0,
      4.01, 4.0200000000000005, 4.03, 4.04, 4.05, 4.0600000000000005, 4.07, 4.08,
      4.09, 4.1, 4.11, 4.12, 4.13, 4.14, 4.15, 4.16, 4.17, 4.18, 4.19, 4.2, 4.21,
      4.22, 4.23, 4.24, 4.25, 4.26, 4.2700000000000005, 4.28, 4.29, 4.3,
      4.3100000000000005, 4.32, 4.33, 4.34, 4.3500000000000005, 4.36, 4.37, 4.38,
      4.39, 4.4, 4.41, 4.42, 4.43, 4.44, 4.45, 4.46, 4.47, 4.48, 4.49, 4.5, 4.51,
      4.5200000000000005, 4.53, 4.54, 4.55, 4.5600000000000005, 4.57, 4.58, 4.59,
      4.6000000000000005, 4.61, 4.62, 4.63, 4.64, 4.65, 4.66, 4.67, 4.68, 4.69,
      4.7, 4.71, 4.72, 4.73, 4.74, 4.75, 4.76, 4.7700000000000005, 4.78, 4.79,
      4.8, 4.8100000000000005, 4.82, 4.83, 4.84, 4.8500000000000005, 4.86, 4.87,
      4.88, 4.89, 4.9, 4.91, 4.92, 4.93, 4.94, 4.95, 4.96, 4.97, 4.98, 4.99, 5.0,
      5.01, 5.0200000000000005, 5.03, 5.04, 5.05, 5.0600000000000005, 5.07, 5.08,
      5.09, 5.1000000000000005, 5.11, 5.12, 5.13, 5.14, 5.15, 5.16, 5.17, 5.18,
      5.19, 5.2, 5.21, 5.22, 5.23, 5.24, 5.25, 5.26, 5.2700000000000005, 5.28,
      5.29, 5.3, 5.3100000000000005, 5.32, 5.33, 5.34, 5.3500000000000005, 5.36,
      5.37, 5.38, 5.39, 5.4, 5.41, 5.42, 5.43, 5.44, 5.45, 5.46, 5.47, 5.48,
      5.49, 5.5, 5.51, 5.5200000000000005, 5.53, 5.54, 5.55, 5.5600000000000005,
      5.57, 5.58, 5.59, 5.6000000000000005, 5.61, 5.62, 5.63, 5.64, 5.65, 5.66,
      5.67, 5.68, 5.69, 5.7, 5.71, 5.72, 5.73, 5.74, 5.75, 5.76,
      5.7700000000000005, 5.78, 5.79, 5.8, 5.8100000000000005, 5.82, 5.83, 5.84,
      5.8500000000000005, 5.86, 5.87, 5.88, 5.89, 5.9, 5.91, 5.92, 5.93, 5.94,
      5.95, 5.96, 5.97, 5.98, 5.99, 6.0, 6.01, 6.0200000000000005, 6.03, 6.04,
      6.05, 6.0600000000000005, 6.07, 6.08, 6.09, 6.1000000000000005, 6.11, 6.12,
      6.13, 6.1400000000000006, 6.15, 6.16, 6.17, 6.18, 6.19, 6.2, 6.21, 6.22,
      6.23, 6.24, 6.25, 6.26, 6.2700000000000005, 6.28, 6.29, 6.3,
      6.3100000000000005, 6.32, 6.33, 6.34, 6.3500000000000005, 6.36, 6.37, 6.38,
      6.3900000000000006, 6.4, 6.41, 6.42, 6.43, 6.44, 6.45, 6.46, 6.47, 6.48,
      6.49, 6.5, 6.51, 6.5200000000000005, 6.53, 6.54, 6.55, 6.5600000000000005,
      6.57, 6.58, 6.59, 6.6000000000000005, 6.61, 6.62, 6.63, 6.6400000000000006,
      6.65, 6.66, 6.67, 6.68, 6.69, 6.7, 6.71, 6.72, 6.73, 6.74, 6.75, 6.76,
      6.7700000000000005, 6.78, 6.79, 6.8, 6.8100000000000005, 6.82, 6.83, 6.84,
      6.8500000000000005, 6.86, 6.87, 6.88, 6.8900000000000006, 6.9, 6.91, 6.92,
      6.93, 6.94, 6.95, 6.96, 6.97, 6.98, 6.99, 7.0, 7.01, 7.0200000000000005,
      7.03, 7.04, 7.05, 7.0600000000000005, 7.07, 7.08, 7.09, 7.1000000000000005,
      7.11, 7.12, 7.13, 7.1400000000000006, 7.15, 7.16, 7.17, 7.18, 7.19, 7.2,
      7.21, 7.22, 7.23, 7.24, 7.25, 7.26, 7.2700000000000005, 7.28, 7.29, 7.3,
      7.3100000000000005, 7.32, 7.33, 7.34, 7.3500000000000005, 7.36, 7.37, 7.38,
      7.3900000000000006, 7.4, 7.41, 7.42, 7.43, 7.44, 7.45, 7.46, 7.47, 7.48,
      7.49, 7.5, 7.51, 7.5200000000000005, 7.53, 7.54, 7.55, 7.5600000000000005,
      7.57, 7.58, 7.59, 7.6000000000000005, 7.61, 7.62, 7.63, 7.6400000000000006,
      7.65, 7.66, 7.67, 7.68, 7.69, 7.7, 7.71, 7.72, 7.73, 7.74, 7.75, 7.76,
      7.7700000000000005, 7.78, 7.79, 7.8, 7.8100000000000005, 7.82, 7.83, 7.84,
      7.8500000000000005, 7.86, 7.87, 7.88, 7.8900000000000006, 7.9, 7.91, 7.92,
      7.9300000000000006, 7.94, 7.95, 7.96, 7.97, 7.98, 7.99, 8.0, 8.01, 8.02,
      8.03, 8.0400000000000009, 8.05, 8.06, 8.07, 8.08, 8.09, 8.1, 8.11,
      8.120000000000001, 8.13, 8.14, 8.15, 8.16, 8.17, 8.18, 8.19, 8.2, 8.21,
      8.22, 8.23, 8.24, 8.25, 8.26, 8.27, 8.28, 8.2900000000000009, 8.3, 8.31,
      8.32, 8.33, 8.34, 8.35, 8.36, 8.370000000000001, 8.38, 8.39, 8.4, 8.41,
      8.42, 8.43, 8.44, 8.45, 8.46, 8.47, 8.48, 8.49, 8.5, 8.51, 8.52, 8.53,
      8.5400000000000009, 8.55, 8.56, 8.57, 8.58, 8.59, 8.6, 8.61,
      8.620000000000001, 8.63, 8.64, 8.65, 8.66, 8.67, 8.68, 8.69,
      8.7000000000000011, 8.71, 8.72, 8.73, 8.74, 8.75, 8.76, 8.77, 8.78,
      8.7900000000000009, 8.8, 8.81, 8.82, 8.83, 8.84, 8.85, 8.86,
      8.870000000000001, 8.88, 8.89, 8.9, 8.91, 8.92, 8.93, 8.94,
      8.9500000000000011, 8.96, 8.97, 8.98, 8.99, 9.0, 9.01, 9.02, 9.03,
      9.0400000000000009, 9.05, 9.06, 9.07, 9.08, 9.09, 9.1, 9.11,
      9.120000000000001, 9.13, 9.14, 9.15, 9.16, 9.17, 9.18, 9.19,
      9.2000000000000011, 9.21, 9.22, 9.23, 9.24, 9.25, 9.26, 9.27, 9.28,
      9.2900000000000009, 9.3, 9.31, 9.32, 9.33, 9.34, 9.35, 9.36,
      9.370000000000001, 9.38, 9.39, 9.4, 9.41, 9.42, 9.43, 9.44,
      9.4500000000000011, 9.46, 9.47, 9.48, 9.49, 9.5, 9.51, 9.52, 9.53,
      9.5400000000000009, 9.55, 9.56, 9.57, 9.58, 9.59, 9.6, 9.61,
      9.620000000000001, 9.63, 9.64, 9.65, 9.66, 9.67, 9.68, 9.69,
      9.7000000000000011, 9.71, 9.72, 9.73, 9.74, 9.75, 9.76, 9.77, 9.78,
      9.7900000000000009, 9.8, 9.81, 9.82, 9.83, 9.84, 9.85, 9.86,
      9.870000000000001, 9.88, 9.89, 9.9, 9.91, 9.92, 9.93, 9.94,
      9.9500000000000011, 9.96, 9.97, 9.98, 9.99, 10.0, 10.01, 10.02, 10.03,
      10.040000000000001, 10.05, 10.06, 10.07, 10.08, 10.09, 10.1, 10.11,
      10.120000000000001, 10.13, 10.14, 10.15, 10.16, 10.17, 10.18, 10.19,
      10.200000000000001, 10.21, 10.22, 10.23, 10.24, 10.25, 10.26, 10.27, 10.28,
      10.290000000000001, 10.3, 10.31, 10.32, 10.33, 10.34, 10.35, 10.36,
      10.370000000000001, 10.38, 10.39, 10.4, 10.41, 10.42, 10.43, 10.44,
      10.450000000000001, 10.46, 10.47, 10.48, 10.49, 10.5, 10.51, 10.52, 10.53,
      10.540000000000001, 10.55, 10.56, 10.57, 10.58, 10.59, 10.6, 10.61,
      10.620000000000001, 10.63, 10.64, 10.65, 10.66, 10.67, 10.68, 10.69,
      10.700000000000001, 10.71, 10.72, 10.73, 10.74, 10.75, 10.76, 10.77, 10.78,
      10.790000000000001, 10.8, 10.81, 10.82, 10.83, 10.84, 10.85, 10.86,
      10.870000000000001, 10.88, 10.89, 10.9, 10.91, 10.92, 10.93, 10.94,
      10.950000000000001, 10.96, 10.97, 10.98, 10.99, 11.0, 11.01, 11.02, 11.03,
      11.040000000000001, 11.05, 11.06, 11.07, 11.08, 11.09, 11.1, 11.11,
      11.120000000000001, 11.13, 11.14, 11.15, 11.16, 11.17, 11.18, 11.19,
      11.200000000000001, 11.21, 11.22, 11.23, 11.24, 11.25, 11.26, 11.27, 11.28,
      11.290000000000001, 11.3, 11.31, 11.32, 11.33, 11.34, 11.35, 11.36,
      11.370000000000001, 11.38, 11.39, 11.4, 11.41, 11.42, 11.43, 11.44,
      11.450000000000001, 11.46, 11.47, 11.48, 11.49, 11.5, 11.51, 11.52, 11.53,
      11.540000000000001, 11.55, 11.56, 11.57, 11.58, 11.59, 11.6, 11.61,
      11.620000000000001, 11.63, 11.64, 11.65, 11.66, 11.67, 11.68, 11.69,
      11.700000000000001, 11.71, 11.72, 11.73, 11.74, 11.75, 11.76, 11.77, 11.78,
      11.790000000000001, 11.8, 11.81, 11.82, 11.83, 11.84, 11.85, 11.86,
      11.870000000000001, 11.88, 11.89, 11.9, 11.91, 11.92, 11.93, 11.94,
      11.950000000000001, 11.96, 11.97, 11.98, 11.99, 12.0, 12.01, 12.02,
      12.030000000000001, 12.040000000000001, 12.05, 12.06, 12.07, 12.08, 12.09,
      12.1, 12.11, 12.120000000000001, 12.13, 12.14, 12.15, 12.16, 12.17, 12.18,
      12.19, 12.200000000000001, 12.21, 12.22, 12.23, 12.24, 12.25, 12.26, 12.27,
      12.280000000000001, 12.290000000000001, 12.3, 12.31, 12.32, 12.33, 12.34,
      12.35, 12.36, 12.370000000000001, 12.38, 12.39, 12.4, 12.41, 12.42, 12.43,
      12.44, 12.450000000000001, 12.46, 12.47, 12.48, 12.49, 12.5, 12.51, 12.52,
      12.530000000000001, 12.540000000000001, 12.55, 12.56, 12.57, 12.58, 12.59,
      12.6, 12.61, 12.620000000000001, 12.63, 12.64, 12.65, 12.66, 12.67, 12.68,
      12.69, 12.700000000000001, 12.71, 12.72, 12.73, 12.74, 12.75, 12.76, 12.77,
      12.780000000000001, 12.790000000000001, 12.8, 12.81, 12.82, 12.83, 12.84,
      12.85, 12.86, 12.870000000000001, 12.88, 12.89, 12.9, 12.91, 12.92, 12.93,
      12.94, 12.950000000000001, 12.96, 12.97, 12.98, 12.99, 13.0, 13.01, 13.02,
      13.030000000000001, 13.040000000000001, 13.05, 13.06, 13.07, 13.08, 13.09,
      13.1, 13.11, 13.120000000000001, 13.13, 13.14, 13.15, 13.16, 13.17, 13.18,
      13.19, 13.200000000000001, 13.21, 13.22, 13.23, 13.24, 13.25, 13.26, 13.27,
      13.280000000000001, 13.290000000000001, 13.3, 13.31, 13.32, 13.33, 13.34,
      13.35, 13.36, 13.370000000000001, 13.38, 13.39, 13.4, 13.41, 13.42, 13.43,
      13.44, 13.450000000000001, 13.46, 13.47, 13.48, 13.49, 13.5, 13.51, 13.52,
      13.530000000000001, 13.540000000000001, 13.55, 13.56, 13.57, 13.58, 13.59,
      13.6, 13.61, 13.620000000000001, 13.63, 13.64, 13.65, 13.66, 13.67, 13.68,
      13.69, 13.700000000000001, 13.71, 13.72, 13.73, 13.74, 13.75, 13.76, 13.77,
      13.780000000000001, 13.790000000000001, 13.8, 13.81, 13.82, 13.83, 13.84,
      13.85, 13.86, 13.870000000000001, 13.88, 13.89, 13.9, 13.91, 13.92, 13.93,
      13.94, 13.950000000000001, 13.96, 13.97, 13.98, 13.99, 14.0, 14.01, 14.02,
      14.030000000000001, 14.040000000000001, 14.05, 14.06, 14.07, 14.08, 14.09,
      14.1, 14.11, 14.120000000000001, 14.13, 14.14, 14.15, 14.16, 14.17, 14.18,
      14.19, 14.200000000000001, 14.21, 14.22, 14.23, 14.24, 14.25, 14.26, 14.27,
      14.280000000000001, 14.290000000000001, 14.3, 14.31, 14.32, 14.33, 14.34,
      14.35, 14.36, 14.370000000000001, 14.38, 14.39, 14.4, 14.41, 14.42, 14.43,
      14.44, 14.450000000000001, 14.46, 14.47, 14.48, 14.49, 14.5, 14.51, 14.52,
      14.530000000000001, 14.540000000000001, 14.55, 14.56, 14.57, 14.58, 14.59,
      14.6, 14.61, 14.620000000000001, 14.63, 14.64, 14.65, 14.66, 14.67, 14.68,
      14.69, 14.700000000000001, 14.71, 14.72, 14.73, 14.74, 14.75, 14.76, 14.77,
      14.780000000000001, 14.790000000000001, 14.8, 14.81, 14.82, 14.83, 14.84,
      14.85, 14.86, 14.870000000000001, 14.88, 14.89, 14.9, 14.91, 14.92, 14.93,
      14.94, 14.950000000000001, 14.96, 14.97, 14.98, 14.99, 15.0, 15.01, 15.02,
      15.030000000000001, 15.040000000000001, 15.05, 15.06, 15.07, 15.08, 15.09,
      15.1, 15.11, 15.120000000000001, 15.13, 15.14, 15.15, 15.16, 15.17, 15.18,
      15.19, 15.200000000000001, 15.21, 15.22, 15.23, 15.24, 15.25, 15.26, 15.27,
      15.280000000000001, 15.290000000000001, 15.3, 15.31, 15.32, 15.33, 15.34,
      15.35, 15.36, 15.370000000000001, 15.38, 15.39, 15.4, 15.41, 15.42, 15.43,
      15.44, 15.450000000000001, 15.46, 15.47, 15.48, 15.49, 15.5, 15.51, 15.52,
      15.530000000000001, 15.540000000000001, 15.55, 15.56, 15.57, 15.58, 15.59,
      15.6, 15.610000000000001, 15.620000000000001, 15.63, 15.64, 15.65, 15.66,
      15.67, 15.68, 15.69, 15.700000000000001, 15.71, 15.72, 15.73, 15.74, 15.75,
      15.76, 15.77, 15.780000000000001, 15.790000000000001, 15.8, 15.81, 15.82,
      15.83, 15.84, 15.85, 15.860000000000001, 15.870000000000001, 15.88, 15.89,
      15.9, 15.91, 15.92, 15.93, 15.94, 15.950000000000001, 15.96, 15.97, 15.98,
      15.99, 16.0, 16.01, 16.02, 16.03, 16.04, 16.05, 16.06, 16.07,
      16.080000000000002, 16.09, 16.1, 16.11, 16.12, 16.13, 16.14, 16.15, 16.16,
      16.17, 16.18, 16.19, 16.2, 16.21, 16.22, 16.23, 16.240000000000002, 16.25,
      16.26, 16.27, 16.28, 16.29, 16.3, 16.31, 16.32, 16.330000000000002, 16.34,
      16.35, 16.36, 16.37, 16.38, 16.39, 16.4, 16.41, 16.42, 16.43, 16.44, 16.45,
      16.46, 16.47, 16.48, 16.490000000000002, 16.5, 16.51, 16.52, 16.53, 16.54,
      16.55, 16.56, 16.57, 16.580000000000002, 16.59, 16.6, 16.61, 16.62, 16.63,
      16.64, 16.65, 16.66, 16.67, 16.68, 16.69, 16.7, 16.71, 16.72, 16.73,
      16.740000000000002, 16.75, 16.76, 16.77, 16.78, 16.79, 16.8, 16.81, 16.82,
      16.830000000000002, 16.84, 16.85, 16.86, 16.87, 16.88, 16.89, 16.9, 16.91,
      16.92, 16.93, 16.94, 16.95, 16.96, 16.97, 16.98, 16.990000000000002, 17.0,
      17.01, 17.02, 17.03, 17.04, 17.05, 17.06, 17.07, 17.080000000000002, 17.09,
      17.1, 17.11, 17.12, 17.13, 17.14, 17.150000000000002, 17.16, 17.17, 17.18,
      17.19, 17.2, 17.21, 17.22, 17.23, 17.240000000000002, 17.25, 17.26, 17.27,
      17.28, 17.29, 17.3, 17.31, 17.32, 17.330000000000002, 17.34, 17.35, 17.36,
      17.37, 17.38, 17.39, 17.400000000000002, 17.41, 17.42, 17.43, 17.44, 17.45,
      17.46, 17.47, 17.48, 17.490000000000002, 17.5, 17.51, 17.52, 17.53, 17.54,
      17.55, 17.56, 17.57, 17.580000000000002, 17.59, 17.6, 17.61, 17.62, 17.63,
      17.64, 17.650000000000002, 17.66, 17.67, 17.68, 17.69, 17.7, 17.71, 17.72,
      17.73, 17.740000000000002, 17.75, 17.76, 17.77, 17.78, 17.79, 17.8, 17.81,
      17.82, 17.830000000000002, 17.84, 17.85, 17.86, 17.87, 17.88, 17.89,
      17.900000000000002, 17.91, 17.92, 17.93, 17.94, 17.95, 17.96, 17.97, 17.98,
      17.990000000000002, 18.0, 18.01, 18.02, 18.03, 18.04, 18.05, 18.06, 18.07,
      18.080000000000002, 18.09, 18.1, 18.11, 18.12, 18.13, 18.14,
      18.150000000000002, 18.16, 18.17, 18.18, 18.19, 18.2, 18.21, 18.22, 18.23,
      18.240000000000002, 18.25, 18.26, 18.27, 18.28, 18.29, 18.3, 18.31, 18.32,
      18.330000000000002, 18.34, 18.35, 18.36, 18.37, 18.38, 18.39,
      18.400000000000002, 18.41, 18.42, 18.43, 18.44, 18.45, 18.46, 18.47, 18.48,
      18.490000000000002, 18.5, 18.51, 18.52, 18.53, 18.54, 18.55, 18.56, 18.57,
      18.580000000000002, 18.59, 18.6, 18.61, 18.62, 18.63, 18.64,
      18.650000000000002, 18.66, 18.67, 18.68, 18.69, 18.7, 18.71, 18.72, 18.73,
      18.740000000000002, 18.75, 18.76, 18.77, 18.78, 18.79, 18.8, 18.81, 18.82,
      18.830000000000002, 18.84, 18.85, 18.86, 18.87, 18.88, 18.89,
      18.900000000000002, 18.91, 18.92, 18.93, 18.94, 18.95, 18.96, 18.97, 18.98,
      18.990000000000002, 19.0, 19.01, 19.02, 19.03, 19.04, 19.05, 19.06, 19.07,
      19.080000000000002, 19.09, 19.1, 19.11, 19.12, 19.13, 19.14,
      19.150000000000002, 19.16, 19.17, 19.18, 19.19, 19.2, 19.21, 19.22, 19.23,
      19.240000000000002, 19.25, 19.26, 19.27, 19.28, 19.29, 19.3, 19.31, 19.32,
      19.330000000000002, 19.34, 19.35, 19.36, 19.37, 19.38, 19.39,
      19.400000000000002, 19.41, 19.42, 19.43, 19.44, 19.45, 19.46, 19.47, 19.48,
      19.490000000000002, 19.5, 19.51, 19.52, 19.53, 19.54, 19.55, 19.56, 19.57,
      19.580000000000002, 19.59, 19.6, 19.61, 19.62, 19.63, 19.64,
      19.650000000000002, 19.66, 19.67, 19.68, 19.69, 19.7, 19.71, 19.72, 19.73,
      19.740000000000002, 19.75, 19.76, 19.77, 19.78, 19.79, 19.8, 19.81, 19.82,
      19.830000000000002, 19.84, 19.85, 19.86, 19.87, 19.88, 19.89,
      19.900000000000002, 19.91, 19.92, 19.93, 19.94, 19.95, 19.96, 19.97, 19.98,
      19.990000000000002, 20.0, 20.01, 20.02, 20.03, 20.04, 20.05, 20.06, 20.07,
      20.080000000000002, 20.09, 20.1, 20.11, 20.12, 20.13, 20.14,
      20.150000000000002, 20.16, 20.17, 20.18, 20.19, 20.2, 20.21, 20.22, 20.23,
      20.240000000000002, 20.25, 20.26, 20.27, 20.28, 20.29, 20.3, 20.31, 20.32,
      20.330000000000002, 20.34, 20.35, 20.36, 20.37, 20.38, 20.39,
      20.400000000000002, 20.41, 20.42, 20.43, 20.44, 20.45, 20.46, 20.47, 20.48,
      20.490000000000002, 20.5, 20.51, 20.52, 20.53, 20.54, 20.55, 20.56, 20.57,
      20.580000000000002, 20.59, 20.6, 20.61, 20.62, 20.63, 20.64,
      20.650000000000002, 20.66, 20.67, 20.68, 20.69, 20.7, 20.71, 20.72, 20.73,
      20.740000000000002, 20.75, 20.76, 20.77, 20.78, 20.79, 20.8, 20.81, 20.82,
      20.830000000000002, 20.84, 20.85, 20.86, 20.87, 20.88, 20.89,
      20.900000000000002, 20.91, 20.92, 20.93, 20.94, 20.95, 20.96, 20.97, 20.98,
      20.990000000000002, 21.0, 21.01, 21.02, 21.03, 21.04, 21.05, 21.06, 21.07,
      21.080000000000002, 21.09, 21.1, 21.11, 21.12, 21.13, 21.14,
      21.150000000000002, 21.16, 21.17, 21.18, 21.19, 21.2, 21.21, 21.22, 21.23,
      21.240000000000002, 21.25, 21.26, 21.27, 21.28, 21.29, 21.3, 21.31, 21.32,
      21.330000000000002, 21.34, 21.35, 21.36, 21.37, 21.38, 21.39,
      21.400000000000002, 21.41, 21.42, 21.43, 21.44, 21.45, 21.46, 21.47, 21.48,
      21.490000000000002, 21.5, 21.51, 21.52, 21.53, 21.54, 21.55, 21.56, 21.57,
      21.580000000000002, 21.59, 21.6, 21.61, 21.62, 21.63, 21.64,
      21.650000000000002, 21.66, 21.67, 21.68, 21.69, 21.7, 21.71, 21.72, 21.73,
      21.740000000000002, 21.75, 21.76, 21.77, 21.78, 21.79, 21.8, 21.81, 21.82,
      21.830000000000002, 21.84, 21.85, 21.86, 21.87, 21.88, 21.89,
      21.900000000000002, 21.91, 21.92, 21.93, 21.94, 21.95, 21.96, 21.97, 21.98,
      21.990000000000002, 22.0, 22.01, 22.02, 22.03, 22.04, 22.05, 22.06, 22.07,
      22.080000000000002, 22.09, 22.1, 22.11, 22.12, 22.13, 22.14,
      22.150000000000002, 22.16, 22.17, 22.18, 22.19, 22.2, 22.21, 22.22, 22.23,
      22.240000000000002, 22.25, 22.26, 22.27, 22.28, 22.29, 22.3, 22.31, 22.32,
      22.330000000000002, 22.34, 22.35, 22.36, 22.37, 22.38, 22.39,
      22.400000000000002, 22.41, 22.42, 22.43, 22.44, 22.45, 22.46, 22.47, 22.48,
      22.490000000000002, 22.5, 22.51, 22.52, 22.53, 22.54, 22.55, 22.56, 22.57,
      22.580000000000002, 22.59, 22.6, 22.61, 22.62, 22.63, 22.64,
      22.650000000000002, 22.66, 22.67, 22.68, 22.69, 22.7, 22.71, 22.72, 22.73,
      22.740000000000002, 22.75, 22.76, 22.77, 22.78, 22.79, 22.8, 22.81, 22.82,
      22.830000000000002, 22.84, 22.85, 22.86, 22.87, 22.88, 22.89,
      22.900000000000002, 22.91, 22.92, 22.93, 22.94, 22.95, 22.96, 22.97, 22.98,
      22.990000000000002, 23.0, 23.01, 23.02, 23.03, 23.04, 23.05, 23.06, 23.07,
      23.080000000000002, 23.09, 23.1, 23.11, 23.12, 23.13, 23.14,
      23.150000000000002, 23.16, 23.17, 23.18, 23.19, 23.2, 23.21, 23.22, 23.23,
      23.240000000000002, 23.25, 23.26, 23.27, 23.28, 23.29, 23.3, 23.31, 23.32,
      23.330000000000002, 23.34, 23.35, 23.36, 23.37, 23.38, 23.39,
      23.400000000000002, 23.41, 23.42, 23.43, 23.44, 23.45, 23.46, 23.47, 23.48,
      23.490000000000002, 23.5, 23.51, 23.52, 23.53, 23.54, 23.55, 23.56, 23.57,
      23.580000000000002, 23.59, 23.6, 23.61, 23.62, 23.63, 23.64,
      23.650000000000002, 23.66, 23.67, 23.68, 23.69, 23.7, 23.71, 23.72, 23.73,
      23.740000000000002, 23.75, 23.76, 23.77, 23.78, 23.79, 23.8, 23.81, 23.82,
      23.830000000000002, 23.84, 23.85, 23.86, 23.87, 23.88, 23.89,
      23.900000000000002, 23.91, 23.92, 23.93, 23.94, 23.95, 23.96, 23.97, 23.98,
      23.990000000000002, 24.0, 24.01, 24.02, 24.03, 24.04, 24.05,
      24.060000000000002, 24.07, 24.080000000000002, 24.09, 24.1, 24.11, 24.12,
      24.13, 24.14, 24.150000000000002, 24.16, 24.17, 24.18, 24.19, 24.2, 24.21,
      24.22, 24.23, 24.240000000000002, 24.25, 24.26, 24.27, 24.28, 24.29, 24.3,
      24.310000000000002, 24.32, 24.330000000000002, 24.34, 24.35, 24.36, 24.37,
      24.38, 24.39, 24.400000000000002, 24.41, 24.42, 24.43, 24.44, 24.45, 24.46,
      24.47, 24.48, 24.490000000000002, 24.5, 24.51, 24.52, 24.53, 24.54, 24.55,
      24.560000000000002, 24.57, 24.580000000000002, 24.59, 24.6, 24.61, 24.62,
      24.63, 24.64, 24.650000000000002, 24.66, 24.67, 24.68, 24.69, 24.7, 24.71,
      24.72, 24.73, 24.740000000000002, 24.75, 24.76, 24.77, 24.78, 24.79, 24.8,
      24.810000000000002, 24.82, 24.830000000000002, 24.84, 24.85, 24.86, 24.87,
      24.88, 24.89, 24.900000000000002, 24.91, 24.92, 24.93, 24.94, 24.95, 24.96,
      24.97, 24.98, 24.990000000000002, 25.0, 25.01, 25.02, 25.03, 25.04, 25.05,
      25.060000000000002, 25.07, 25.080000000000002, 25.09, 25.1, 25.11, 25.12,
      25.13, 25.14, 25.150000000000002, 25.16, 25.17, 25.18, 25.19, 25.2, 25.21,
      25.22, 25.23, 25.240000000000002, 25.25, 25.26, 25.27, 25.28, 25.29, 25.3,
      25.310000000000002, 25.32, 25.330000000000002, 25.34, 25.35, 25.36, 25.37,
      25.38, 25.39, 25.400000000000002, 25.41, 25.42, 25.43, 25.44, 25.45, 25.46,
      25.47, 25.48, 25.490000000000002, 25.5, 25.51, 25.52, 25.53, 25.54, 25.55,
      25.560000000000002, 25.57, 25.580000000000002, 25.59, 25.6, 25.61, 25.62,
      25.63, 25.64, 25.650000000000002, 25.66, 25.67, 25.68, 25.69, 25.7, 25.71,
      25.72, 25.73, 25.740000000000002, 25.75, 25.76, 25.77, 25.78, 25.79, 25.8,
      25.810000000000002, 25.82, 25.830000000000002, 25.84, 25.85, 25.86, 25.87,
      25.88, 25.89, 25.900000000000002, 25.91, 25.92, 25.93, 25.94, 25.95, 25.96,
      25.97, 25.98, 25.990000000000002, 26.0, 26.01, 26.02, 26.03, 26.04, 26.05,
      26.060000000000002, 26.07, 26.080000000000002, 26.09, 26.1, 26.11, 26.12,
      26.13, 26.14, 26.150000000000002, 26.16, 26.17, 26.18, 26.19, 26.2, 26.21,
      26.22, 26.23, 26.240000000000002, 26.25, 26.26, 26.27, 26.28, 26.29, 26.3,
      26.310000000000002, 26.32, 26.330000000000002, 26.34, 26.35, 26.36, 26.37,
      26.38, 26.39, 26.400000000000002, 26.41, 26.42, 26.43, 26.44, 26.45, 26.46,
      26.47, 26.48, 26.490000000000002, 26.5, 26.51, 26.52, 26.53, 26.54, 26.55,
      26.560000000000002, 26.57, 26.580000000000002, 26.59, 26.6, 26.61, 26.62,
      26.63, 26.64, 26.650000000000002, 26.66, 26.67, 26.68, 26.69, 26.7, 26.71,
      26.72, 26.73, 26.740000000000002, 26.75, 26.76, 26.77, 26.78, 26.79, 26.8,
      26.810000000000002, 26.82, 26.830000000000002, 26.84, 26.85, 26.86, 26.87,
      26.88, 26.89, 26.900000000000002, 26.91, 26.92, 26.93, 26.94, 26.95, 26.96,
      26.97, 26.98, 26.990000000000002, 27.0, 27.01, 27.02, 27.03, 27.04, 27.05,
      27.060000000000002, 27.07, 27.080000000000002, 27.09, 27.1, 27.11, 27.12,
      27.13, 27.14, 27.150000000000002, 27.16, 27.17, 27.18, 27.19, 27.2, 27.21,
      27.22, 27.23, 27.240000000000002, 27.25, 27.26, 27.27, 27.28, 27.29, 27.3,
      27.310000000000002, 27.32, 27.330000000000002, 27.34, 27.35, 27.36, 27.37,
      27.38, 27.39, 27.400000000000002, 27.41, 27.42, 27.43, 27.44, 27.45, 27.46,
      27.47, 27.48, 27.490000000000002, 27.5, 27.51, 27.52, 27.53, 27.54, 27.55,
      27.560000000000002, 27.57, 27.580000000000002, 27.59, 27.6, 27.61, 27.62,
      27.63, 27.64, 27.650000000000002, 27.66, 27.67, 27.68, 27.69, 27.7, 27.71,
      27.72, 27.73, 27.740000000000002, 27.75, 27.76, 27.77, 27.78, 27.79, 27.8,
      27.810000000000002, 27.82, 27.830000000000002, 27.84, 27.85, 27.86, 27.87,
      27.88, 27.89, 27.900000000000002, 27.91, 27.92, 27.93, 27.94, 27.95, 27.96,
      27.97, 27.98, 27.990000000000002, 28.0, 28.01, 28.02, 28.03, 28.04, 28.05,
      28.060000000000002, 28.07, 28.080000000000002, 28.09, 28.1, 28.11, 28.12,
      28.13, 28.14, 28.150000000000002, 28.16, 28.17, 28.18, 28.19, 28.2, 28.21,
      28.22, 28.23, 28.240000000000002, 28.25, 28.26, 28.27, 28.28, 28.29, 28.3,
      28.310000000000002, 28.32, 28.330000000000002, 28.34, 28.35, 28.36, 28.37,
      28.38, 28.39, 28.400000000000002, 28.41, 28.42, 28.43, 28.44, 28.45, 28.46,
      28.47, 28.48, 28.490000000000002, 28.5, 28.51, 28.52, 28.53, 28.54, 28.55,
      28.560000000000002, 28.57, 28.580000000000002, 28.59, 28.6, 28.61, 28.62,
      28.63, 28.64, 28.650000000000002, 28.66, 28.67, 28.68, 28.69, 28.7, 28.71,
      28.72, 28.73, 28.740000000000002, 28.75, 28.76, 28.77, 28.78, 28.79, 28.8,
      28.810000000000002, 28.82, 28.830000000000002, 28.84, 28.85, 28.86, 28.87,
      28.88, 28.89, 28.900000000000002, 28.91, 28.92, 28.93, 28.94, 28.95, 28.96,
      28.97, 28.98, 28.990000000000002, 29.0, 29.01, 29.02, 29.03, 29.04, 29.05,
      29.060000000000002, 29.07, 29.080000000000002, 29.09, 29.1, 29.11, 29.12,
      29.13, 29.14, 29.150000000000002, 29.16, 29.17, 29.18, 29.19, 29.2, 29.21,
      29.22, 29.23, 29.240000000000002, 29.25, 29.26, 29.27, 29.28, 29.29, 29.3,
      29.310000000000002, 29.32, 29.330000000000002, 29.34, 29.35, 29.36, 29.37,
      29.38, 29.39, 29.400000000000002, 29.41, 29.42, 29.43, 29.44, 29.45, 29.46,
      29.47, 29.48, 29.490000000000002, 29.5, 29.51, 29.52, 29.53, 29.54, 29.55,
      29.560000000000002, 29.57, 29.580000000000002, 29.59, 29.6, 29.61, 29.62,
      29.63, 29.64, 29.650000000000002, 29.66, 29.67, 29.68, 29.69, 29.7, 29.71,
      29.72, 29.73, 29.740000000000002, 29.75, 29.76, 29.77, 29.78, 29.79, 29.8,
      29.810000000000002, 29.82, 29.830000000000002, 29.84, 29.85, 29.86, 29.87,
      29.88, 29.89, 29.900000000000002, 29.91, 29.92, 29.93, 29.94, 29.95, 29.96,
      29.97, 29.98, 29.990000000000002, 30.0, 30.01, 30.02, 30.03, 30.04, 30.05,
      30.060000000000002, 30.07, 30.080000000000002, 30.09, 30.1, 30.11, 30.12,
      30.13, 30.14, 30.150000000000002, 30.16, 30.17, 30.18, 30.19, 30.2, 30.21,
      30.22, 30.23, 30.240000000000002, 30.25, 30.26, 30.27, 30.28, 30.29, 30.3,
      30.310000000000002, 30.32, 30.330000000000002, 30.34, 30.35, 30.36, 30.37,
      30.38, 30.39, 30.400000000000002, 30.41, 30.42, 30.43, 30.44, 30.45, 30.46,
      30.47, 30.48, 30.490000000000002, 30.5, 30.51, 30.52, 30.53, 30.54, 30.55,
      30.560000000000002, 30.57, 30.580000000000002, 30.59, 30.6, 30.61, 30.62,
      30.63, 30.64, 30.650000000000002, 30.66, 30.67, 30.68, 30.69, 30.7, 30.71,
      30.72, 30.73, 30.740000000000002, 30.75, 30.76, 30.77, 30.78, 30.79, 30.8,
      30.810000000000002, 30.82, 30.830000000000002, 30.84, 30.85, 30.86, 30.87,
      30.88, 30.89, 30.900000000000002, 30.91, 30.92, 30.93, 30.94, 30.95, 30.96,
      30.970000000000002, 30.98, 30.990000000000002, 31.0, 31.01, 31.02, 31.03,
      31.04, 31.05, 31.060000000000002, 31.07, 31.080000000000002, 31.09, 31.1,
      31.11, 31.12, 31.13, 31.14, 31.150000000000002, 31.16, 31.17, 31.18, 31.19,
      31.2, 31.21, 31.220000000000002, 31.23, 31.240000000000002, 31.25, 31.26,
      31.27, 31.28, 31.29, 31.3, 31.310000000000002, 31.32, 31.330000000000002,
      31.34, 31.35, 31.36, 31.37, 31.38, 31.39, 31.400000000000002, 31.41, 31.42,
      31.43, 31.44, 31.45, 31.46, 31.470000000000002, 31.48, 31.490000000000002,
      31.5, 31.51, 31.52, 31.53, 31.54, 31.55, 31.560000000000002, 31.57,
      31.580000000000002, 31.59, 31.6, 31.61, 31.62, 31.63, 31.64,
      31.650000000000002, 31.66, 31.67, 31.68, 31.69, 31.7, 31.71,
      31.720000000000002, 31.73, 31.740000000000002, 31.75, 31.76, 31.77, 31.78,
      31.79, 31.8, 31.810000000000002, 31.82, 31.830000000000002, 31.84, 31.85,
      31.86, 31.87, 31.88, 31.89, 31.900000000000002, 31.91, 31.92, 31.93, 31.94,
      31.95, 31.96, 31.970000000000002, 31.98, 31.990000000000002, 32.0, 32.01,
      32.02, 32.03, 32.04, 32.05, 32.06, 32.07, 32.08, 32.09, 32.1, 32.11, 32.12,
      32.13, 32.14, 32.15, 32.160000000000004, 32.17, 32.18, 32.19, 32.2, 32.21,
      32.22, 32.230000000000004, 32.24, 32.25, 32.26, 32.27, 32.28, 32.29, 32.3,
      32.31, 32.32, 32.33, 32.34, 32.35, 32.36, 32.37, 32.38, 32.39, 32.4,
      32.410000000000004, 32.42, 32.43, 32.44, 32.45, 32.46, 32.47,
      32.480000000000004, 32.49, 32.5, 32.51, 32.52, 32.53, 32.54, 32.55, 32.56,
      32.57, 32.58, 32.59, 32.6, 32.61, 32.62, 32.63, 32.64, 32.65,
      32.660000000000004, 32.67, 32.68, 32.69, 32.7, 32.71, 32.72,
      32.730000000000004, 32.74, 32.75, 32.76, 32.77, 32.78, 32.79, 32.8, 32.81,
      32.82, 32.83, 32.84, 32.85, 32.86, 32.87, 32.88, 32.89, 32.9,
      32.910000000000004, 32.92, 32.93, 32.94, 32.95, 32.96, 32.97,
      32.980000000000004, 32.99, 33.0, 33.01, 33.02, 33.03, 33.04, 33.05, 33.06,
      33.07, 33.08, 33.09, 33.1, 33.11, 33.12, 33.13, 33.14, 33.15,
      33.160000000000004, 33.17, 33.18, 33.19, 33.2, 33.21, 33.22,
      33.230000000000004, 33.24, 33.25, 33.26, 33.27, 33.28, 33.29, 33.3, 33.31,
      33.32, 33.33, 33.34, 33.35, 33.36, 33.37, 33.38, 33.39, 33.4,
      33.410000000000004, 33.42, 33.43, 33.44, 33.45, 33.46, 33.47,
      33.480000000000004, 33.49, 33.5, 33.51, 33.52, 33.53, 33.54, 33.55, 33.56,
      33.57, 33.58, 33.59, 33.6, 33.61, 33.62, 33.63, 33.64, 33.65,
      33.660000000000004, 33.67, 33.68, 33.69, 33.7, 33.71, 33.72,
      33.730000000000004, 33.74, 33.75, 33.76, 33.77, 33.78, 33.79, 33.8, 33.81,
      33.82, 33.83, 33.84, 33.85, 33.86, 33.87, 33.88, 33.89, 33.9,
      33.910000000000004, 33.92, 33.93, 33.94, 33.95, 33.96, 33.97,
      33.980000000000004, 33.99, 34.0, 34.01, 34.02, 34.03, 34.04, 34.05, 34.06,
      34.07, 34.08, 34.09, 34.1, 34.11, 34.12, 34.13, 34.14, 34.15,
      34.160000000000004, 34.17, 34.18, 34.19, 34.2, 34.21, 34.22,
      34.230000000000004, 34.24, 34.25, 34.26, 34.27, 34.28, 34.29,
      34.300000000000004, 34.31, 34.32, 34.33, 34.34, 34.35, 34.36, 34.37, 34.38,
      34.39, 34.4, 34.410000000000004, 34.42, 34.43, 34.44, 34.45, 34.46, 34.47,
      34.480000000000004, 34.49, 34.5, 34.51, 34.52, 34.53, 34.54,
      34.550000000000004, 34.56, 34.57, 34.58, 34.59, 34.6, 34.61, 34.62, 34.63,
      34.64, 34.65, 34.660000000000004, 34.67, 34.68, 34.69, 34.7, 34.71, 34.72,
      34.730000000000004, 34.74, 34.75, 34.76, 34.77, 34.78, 34.79,
      34.800000000000004, 34.81, 34.82, 34.83, 34.84, 34.85, 34.86, 34.87, 34.88,
      34.89, 34.9, 34.910000000000004, 34.92, 34.93, 34.94, 34.95, 34.96, 34.97,
      34.980000000000004, 34.99, 35.0, 35.01, 35.02, 35.03, 35.04,
      35.050000000000004, 35.06, 35.07, 35.08, 35.09, 35.1, 35.11, 35.12, 35.13,
      35.14, 35.15, 35.160000000000004, 35.17, 35.18, 35.19, 35.2, 35.21, 35.22,
      35.230000000000004, 35.24, 35.25, 35.26, 35.27, 35.28, 35.29,
      35.300000000000004, 35.31, 35.32, 35.33, 35.34, 35.35, 35.36, 35.37, 35.38,
      35.39, 35.4, 35.410000000000004, 35.42, 35.43, 35.44, 35.45, 35.46, 35.47,
      35.480000000000004, 35.49, 35.5, 35.51, 35.52, 35.53, 35.54,
      35.550000000000004, 35.56, 35.57, 35.58, 35.59, 35.6, 35.61, 35.62, 35.63,
      35.64, 35.65, 35.660000000000004, 35.67, 35.68, 35.69, 35.7, 35.71, 35.72,
      35.730000000000004, 35.74, 35.75, 35.76, 35.77, 35.78, 35.79,
      35.800000000000004, 35.81, 35.82, 35.83, 35.84, 35.85, 35.86, 35.87, 35.88,
      35.89, 35.9, 35.910000000000004, 35.92, 35.93, 35.94, 35.95, 35.96, 35.97,
      35.980000000000004, 35.99, 36.0, 36.01, 36.02, 36.03, 36.04,
      36.050000000000004, 36.06, 36.07, 36.08, 36.09, 36.1, 36.11, 36.12, 36.13,
      36.14, 36.15, 36.160000000000004, 36.17, 36.18, 36.19, 36.2, 36.21, 36.22,
      36.230000000000004, 36.24, 36.25, 36.26, 36.27, 36.28, 36.29,
      36.300000000000004, 36.31, 36.32, 36.33, 36.34, 36.35, 36.36, 36.37, 36.38,
      36.39, 36.4, 36.410000000000004, 36.42, 36.43, 36.44, 36.45, 36.46, 36.47,
      36.480000000000004, 36.49, 36.5, 36.51, 36.52, 36.53, 36.54,
      36.550000000000004, 36.56, 36.57, 36.58, 36.59, 36.6, 36.61, 36.62, 36.63,
      36.64, 36.65, 36.660000000000004, 36.67, 36.68, 36.69, 36.7, 36.71, 36.72,
      36.730000000000004, 36.74, 36.75, 36.76, 36.77, 36.78, 36.79,
      36.800000000000004, 36.81, 36.82, 36.83, 36.84, 36.85, 36.86, 36.87, 36.88,
      36.89, 36.9, 36.910000000000004, 36.92, 36.93, 36.94, 36.95, 36.96, 36.97,
      36.980000000000004, 36.99, 37.0, 37.01, 37.02, 37.03, 37.04,
      37.050000000000004, 37.06, 37.07, 37.08, 37.09, 37.1, 37.11, 37.12, 37.13,
      37.14, 37.15, 37.160000000000004, 37.17, 37.18, 37.19, 37.2, 37.21, 37.22,
      37.230000000000004, 37.24, 37.25, 37.26, 37.27, 37.28, 37.29,
      37.300000000000004, 37.31, 37.32, 37.33, 37.34, 37.35, 37.36, 37.37, 37.38,
      37.39, 37.4, 37.410000000000004, 37.42, 37.43, 37.44, 37.45, 37.46, 37.47,
      37.480000000000004, 37.49, 37.5, 37.51, 37.52, 37.53, 37.54,
      37.550000000000004, 37.56, 37.57, 37.58, 37.59, 37.6, 37.61, 37.62, 37.63,
      37.64, 37.65, 37.660000000000004, 37.67, 37.68, 37.69, 37.7, 37.71, 37.72,
      37.730000000000004, 37.74, 37.75, 37.76, 37.77, 37.78, 37.79,
      37.800000000000004, 37.81, 37.82, 37.83, 37.84, 37.85, 37.86, 37.87, 37.88,
      37.89, 37.9, 37.910000000000004, 37.92, 37.93, 37.94, 37.95, 37.96, 37.97,
      37.980000000000004, 37.99, 38.0, 38.01, 38.02, 38.03, 38.04,
      38.050000000000004, 38.06, 38.07, 38.08, 38.09, 38.1, 38.11, 38.12, 38.13,
      38.14, 38.15, 38.160000000000004, 38.17, 38.18, 38.19, 38.2, 38.21, 38.22,
      38.230000000000004, 38.24, 38.25, 38.26, 38.27, 38.28, 38.29,
      38.300000000000004, 38.31, 38.32, 38.33, 38.34, 38.35, 38.36, 38.37, 38.38,
      38.39, 38.4, 38.410000000000004, 38.42, 38.43, 38.44, 38.45, 38.46, 38.47,
      38.480000000000004, 38.49, 38.5, 38.51, 38.52, 38.53, 38.54,
      38.550000000000004, 38.56, 38.57, 38.58, 38.59, 38.6, 38.61, 38.62, 38.63,
      38.64, 38.65, 38.660000000000004, 38.67, 38.68, 38.69, 38.7, 38.71, 38.72,
      38.730000000000004, 38.74, 38.75, 38.76, 38.77, 38.78, 38.79,
      38.800000000000004, 38.81, 38.82, 38.83, 38.84, 38.85, 38.86, 38.87, 38.88,
      38.89, 38.9, 38.910000000000004, 38.92, 38.93, 38.94, 38.95, 38.96, 38.97,
      38.980000000000004, 38.99, 39.0, 39.01, 39.02, 39.03, 39.04,
      39.050000000000004, 39.06, 39.07, 39.08, 39.09, 39.1, 39.11, 39.12, 39.13,
      39.14, 39.15, 39.160000000000004, 39.17, 39.18, 39.19, 39.2, 39.21, 39.22,
      39.230000000000004, 39.24, 39.25, 39.26, 39.27, 39.28, 39.29,
      39.300000000000004, 39.31, 39.32, 39.33, 39.34, 39.35, 39.36, 39.37, 39.38,
      39.39, 39.4, 39.410000000000004, 39.42, 39.43, 39.44, 39.45, 39.46, 39.47,
      39.480000000000004, 39.49, 39.5, 39.51, 39.52, 39.53, 39.54,
      39.550000000000004, 39.56, 39.57, 39.58, 39.59, 39.6, 39.61, 39.62, 39.63,
      39.64, 39.65, 39.660000000000004, 39.67, 39.68, 39.69, 39.7, 39.71, 39.72,
      39.730000000000004, 39.74, 39.75, 39.76, 39.77, 39.78, 39.79,
      39.800000000000004, 39.81, 39.82, 39.83, 39.84, 39.85, 39.86, 39.87, 39.88,
      39.89, 39.9, 39.910000000000004, 39.92, 39.93, 39.94, 39.95, 39.96, 39.97,
      39.980000000000004, 39.99, 40.0, 40.01, 40.02, 40.03, 40.04,
      40.050000000000004, 40.06, 40.07, 40.08, 40.09, 40.1, 40.11, 40.12, 40.13,
      40.14, 40.15, 40.160000000000004, 40.17, 40.18, 40.19, 40.2, 40.21, 40.22,
      40.230000000000004, 40.24, 40.25, 40.26, 40.27, 40.28, 40.29,
      40.300000000000004, 40.31, 40.32, 40.33, 40.34, 40.35, 40.36, 40.37, 40.38,
      40.39, 40.4, 40.410000000000004, 40.42, 40.43, 40.44, 40.45, 40.46, 40.47,
      40.480000000000004, 40.49, 40.5, 40.51, 40.52, 40.53, 40.54,
      40.550000000000004, 40.56, 40.57, 40.58, 40.59, 40.6, 40.61, 40.62, 40.63,
      40.64, 40.65, 40.660000000000004, 40.67, 40.68, 40.69, 40.7, 40.71, 40.72,
      40.730000000000004, 40.74, 40.75, 40.76, 40.77, 40.78, 40.79,
      40.800000000000004, 40.81, 40.82, 40.83, 40.84, 40.85, 40.86, 40.87, 40.88,
      40.89, 40.9, 40.910000000000004, 40.92, 40.93, 40.94, 40.95, 40.96, 40.97,
      40.980000000000004, 40.99, 41.0, 41.01, 41.02, 41.03, 41.04,
      41.050000000000004, 41.06, 41.07, 41.08, 41.09, 41.1, 41.11, 41.12, 41.13,
      41.14, 41.15, 41.160000000000004, 41.17, 41.18, 41.19, 41.2, 41.21, 41.22,
      41.230000000000004, 41.24, 41.25, 41.26, 41.27, 41.28, 41.29,
      41.300000000000004, 41.31, 41.32, 41.33, 41.34, 41.35, 41.36, 41.37, 41.38,
      41.39, 41.4, 41.410000000000004, 41.42, 41.43, 41.44, 41.45, 41.46, 41.47,
      41.480000000000004, 41.49, 41.5, 41.51, 41.52, 41.53, 41.54,
      41.550000000000004, 41.56, 41.57, 41.58, 41.59, 41.6, 41.61, 41.62, 41.63,
      41.64, 41.65, 41.660000000000004, 41.67, 41.68, 41.69, 41.7, 41.71, 41.72,
      41.730000000000004, 41.74, 41.75, 41.76, 41.77, 41.78, 41.79,
      41.800000000000004, 41.81, 41.82, 41.83, 41.84, 41.85, 41.86, 41.87, 41.88,
      41.89, 41.9, 41.910000000000004, 41.92, 41.93, 41.94, 41.95, 41.96, 41.97,
      41.980000000000004, 41.99, 42.0, 42.01, 42.02, 42.03, 42.04,
      42.050000000000004, 42.06, 42.07, 42.08, 42.09, 42.1, 42.11, 42.12, 42.13,
      42.14, 42.15, 42.160000000000004, 42.17, 42.18, 42.19, 42.2, 42.21, 42.22,
      42.230000000000004, 42.24, 42.25, 42.26, 42.27, 42.28, 42.29,
      42.300000000000004, 42.31, 42.32, 42.33, 42.34, 42.35, 42.36, 42.37, 42.38,
      42.39, 42.4, 42.410000000000004, 42.42, 42.43, 42.44, 42.45, 42.46, 42.47,
      42.480000000000004, 42.49, 42.5, 42.51, 42.52, 42.53, 42.54,
      42.550000000000004, 42.56, 42.57, 42.58, 42.59, 42.6, 42.61, 42.62, 42.63,
      42.64, 42.65, 42.660000000000004, 42.67, 42.68, 42.69, 42.7, 42.71, 42.72,
      42.730000000000004, 42.74, 42.75, 42.76, 42.77, 42.78, 42.79,
      42.800000000000004, 42.81, 42.82, 42.83, 42.84, 42.85, 42.86, 42.87, 42.88,
      42.89, 42.9, 42.910000000000004, 42.92, 42.93, 42.94, 42.95, 42.96, 42.97,
      42.980000000000004, 42.99, 43.0, 43.01, 43.02, 43.03, 43.04,
      43.050000000000004, 43.06, 43.07, 43.08, 43.09, 43.1, 43.11, 43.12, 43.13,
      43.14, 43.15, 43.160000000000004, 43.17, 43.18, 43.19, 43.2, 43.21, 43.22,
      43.230000000000004, 43.24, 43.25, 43.26, 43.27, 43.28, 43.29,
      43.300000000000004, 43.31, 43.32, 43.33, 43.34, 43.35, 43.36, 43.37, 43.38,
      43.39, 43.4, 43.410000000000004, 43.42, 43.43, 43.44, 43.45, 43.46, 43.47,
      43.480000000000004, 43.49, 43.5, 43.51, 43.52, 43.53, 43.54,
      43.550000000000004, 43.56, 43.57, 43.58, 43.59, 43.6, 43.61, 43.62, 43.63,
      43.64, 43.65, 43.660000000000004, 43.67, 43.68, 43.69, 43.7, 43.71, 43.72,
      43.730000000000004, 43.74, 43.75, 43.76, 43.77, 43.78, 43.79,
      43.800000000000004, 43.81, 43.82, 43.83, 43.84, 43.85, 43.86, 43.87, 43.88,
      43.89, 43.9, 43.910000000000004, 43.92, 43.93, 43.94, 43.95, 43.96, 43.97,
      43.980000000000004, 43.99, 44.0, 44.01, 44.02, 44.03, 44.04,
      44.050000000000004, 44.06, 44.07, 44.08, 44.09, 44.1, 44.11, 44.12, 44.13,
      44.14, 44.15, 44.160000000000004, 44.17, 44.18, 44.19, 44.2, 44.21, 44.22,
      44.230000000000004, 44.24, 44.25, 44.26, 44.27, 44.28, 44.29,
      44.300000000000004, 44.31, 44.32, 44.33, 44.34, 44.35, 44.36, 44.37, 44.38,
      44.39, 44.4, 44.410000000000004, 44.42, 44.43, 44.44, 44.45, 44.46, 44.47,
      44.480000000000004, 44.49, 44.5, 44.51, 44.52, 44.53, 44.54,
      44.550000000000004, 44.56, 44.57, 44.58, 44.59, 44.6, 44.61, 44.62, 44.63,
      44.64, 44.65, 44.660000000000004, 44.67, 44.68, 44.69, 44.7, 44.71, 44.72,
      44.730000000000004, 44.74, 44.75, 44.76, 44.77, 44.78, 44.79,
      44.800000000000004, 44.81, 44.82, 44.83, 44.84, 44.85, 44.86, 44.87, 44.88,
      44.89, 44.9, 44.910000000000004, 44.92, 44.93, 44.94, 44.95, 44.96, 44.97,
      44.980000000000004, 44.99, 45.0, 45.01, 45.02, 45.03, 45.04,
      45.050000000000004, 45.06, 45.07, 45.08, 45.09, 45.1, 45.11, 45.12, 45.13,
      45.14, 45.15, 45.160000000000004, 45.17, 45.18, 45.19, 45.2, 45.21, 45.22,
      45.230000000000004, 45.24, 45.25, 45.26, 45.27, 45.28, 45.29,
      45.300000000000004, 45.31, 45.32, 45.33, 45.34, 45.35, 45.36, 45.37, 45.38,
      45.39, 45.4, 45.410000000000004, 45.42, 45.43, 45.44, 45.45, 45.46, 45.47,
      45.480000000000004, 45.49, 45.5, 45.51, 45.52, 45.53, 45.54,
      45.550000000000004, 45.56, 45.57, 45.58, 45.59, 45.6, 45.61, 45.62, 45.63,
      45.64, 45.65, 45.660000000000004, 45.67, 45.68, 45.69, 45.7, 45.71, 45.72,
      45.730000000000004, 45.74, 45.75, 45.76, 45.77, 45.78, 45.79,
      45.800000000000004, 45.81, 45.82, 45.83, 45.84, 45.85, 45.86, 45.87, 45.88,
      45.89, 45.9, 45.910000000000004, 45.92, 45.93, 45.94, 45.95, 45.96, 45.97,
      45.980000000000004, 45.99, 46.0, 46.01, 46.02, 46.03, 46.04,
      46.050000000000004, 46.06, 46.07, 46.08, 46.09, 46.1, 46.11, 46.12, 46.13,
      46.14, 46.15, 46.160000000000004, 46.17, 46.18, 46.19, 46.2, 46.21, 46.22,
      46.230000000000004, 46.24, 46.25, 46.26, 46.27, 46.28, 46.29,
      46.300000000000004, 46.31, 46.32, 46.33, 46.34, 46.35, 46.36, 46.37, 46.38,
      46.39, 46.4, 46.410000000000004, 46.42, 46.43, 46.44, 46.45, 46.46, 46.47,
      46.480000000000004, 46.49, 46.5, 46.51, 46.52, 46.53, 46.54,
      46.550000000000004, 46.56, 46.57, 46.58, 46.59, 46.6, 46.61, 46.62, 46.63,
      46.64, 46.65, 46.660000000000004, 46.67, 46.68, 46.69, 46.7, 46.71, 46.72,
      46.730000000000004, 46.74, 46.75, 46.76, 46.77, 46.78, 46.79,
      46.800000000000004, 46.81, 46.82, 46.83, 46.84, 46.85, 46.86, 46.87, 46.88,
      46.89, 46.9, 46.910000000000004, 46.92, 46.93, 46.94, 46.95, 46.96, 46.97,
      46.980000000000004, 46.99, 47.0, 47.01, 47.02, 47.03, 47.04,
      47.050000000000004, 47.06, 47.07, 47.08, 47.09, 47.1, 47.11, 47.12, 47.13,
      47.14, 47.15, 47.160000000000004, 47.17, 47.18, 47.19, 47.2, 47.21, 47.22,
      47.230000000000004, 47.24, 47.25, 47.26, 47.27, 47.28, 47.29,
      47.300000000000004, 47.31, 47.32, 47.33, 47.34, 47.35, 47.36, 47.37, 47.38,
      47.39, 47.4, 47.410000000000004, 47.42, 47.43, 47.44, 47.45, 47.46, 47.47,
      47.480000000000004, 47.49, 47.5, 47.51, 47.52, 47.53, 47.54,
      47.550000000000004, 47.56, 47.57, 47.58, 47.59, 47.6, 47.61, 47.62, 47.63,
      47.64, 47.65, 47.660000000000004, 47.67, 47.68, 47.69, 47.7, 47.71, 47.72,
      47.730000000000004, 47.74, 47.75, 47.76, 47.77, 47.78, 47.79,
      47.800000000000004, 47.81, 47.82, 47.83, 47.84, 47.85, 47.86,
      47.870000000000005, 47.88, 47.89, 47.9, 47.910000000000004, 47.92, 47.93,
      47.94, 47.95, 47.96, 47.97, 47.980000000000004, 47.99, 48.0, 48.01, 48.02,
      48.03, 48.04, 48.050000000000004, 48.06, 48.07, 48.08, 48.09, 48.1, 48.11,
      48.120000000000005, 48.13, 48.14, 48.15, 48.160000000000004, 48.17, 48.18,
      48.19, 48.2, 48.21, 48.22, 48.230000000000004, 48.24, 48.25, 48.26, 48.27,
      48.28, 48.29, 48.300000000000004, 48.31, 48.32, 48.33, 48.34, 48.35, 48.36,
      48.370000000000005, 48.38, 48.39, 48.4, 48.410000000000004, 48.42, 48.43,
      48.44, 48.45, 48.46, 48.47, 48.480000000000004, 48.49, 48.5, 48.51, 48.52,
      48.53, 48.54, 48.550000000000004, 48.56, 48.57, 48.58, 48.59, 48.6, 48.61,
      48.620000000000005, 48.63, 48.64, 48.65, 48.660000000000004, 48.67, 48.68,
      48.69, 48.7, 48.71, 48.72, 48.730000000000004, 48.74, 48.75, 48.76, 48.77,
      48.78, 48.79, 48.800000000000004, 48.81, 48.82, 48.83, 48.84, 48.85, 48.86,
      48.870000000000005, 48.88, 48.89, 48.9, 48.910000000000004, 48.92, 48.93,
      48.94, 48.95, 48.96, 48.97, 48.980000000000004, 48.99, 49.0, 49.01, 49.02,
      49.03, 49.04, 49.050000000000004, 49.06, 49.07, 49.08, 49.09, 49.1, 49.11,
      49.120000000000005, 49.13, 49.14, 49.15, 49.160000000000004, 49.17, 49.18,
      49.19, 49.2, 49.21, 49.22, 49.230000000000004, 49.24, 49.25, 49.26, 49.27,
      49.28, 49.29, 49.300000000000004, 49.31, 49.32, 49.33, 49.34, 49.35, 49.36,
      49.370000000000005, 49.38, 49.39, 49.4, 49.410000000000004, 49.42, 49.43,
      49.44, 49.45, 49.46, 49.47, 49.480000000000004, 49.49, 49.5, 49.51, 49.52,
      49.53, 49.54, 49.550000000000004, 49.56, 49.57, 49.58, 49.59, 49.6, 49.61,
      49.620000000000005, 49.63, 49.64, 49.65, 49.660000000000004, 49.67, 49.68,
      49.69, 49.7, 49.71, 49.72, 49.730000000000004, 49.74, 49.75, 49.76, 49.77,
      49.78, 49.79, 49.800000000000004, 49.81, 49.82, 49.83, 49.84, 49.85, 49.86,
      49.870000000000005, 49.88, 49.89, 49.9, 49.910000000000004, 49.92, 49.93,
      49.94, 49.95, 49.96, 49.97, 49.980000000000004, 49.99, 50.0, 50.01, 50.02,
      50.03, 50.04, 50.050000000000004, 50.06, 50.07, 50.08, 50.09, 50.1, 50.11,
      50.120000000000005, 50.13, 50.14, 50.15, 50.160000000000004, 50.17, 50.18,
      50.19, 50.2, 50.21, 50.22, 50.230000000000004, 50.24, 50.25, 50.26, 50.27,
      50.28, 50.29, 50.300000000000004, 50.31, 50.32, 50.33, 50.34, 50.35, 50.36,
      50.370000000000005, 50.38, 50.39, 50.4, 50.410000000000004, 50.42, 50.43,
      50.44, 50.45, 50.46, 50.47, 50.480000000000004, 50.49, 50.5, 50.51, 50.52,
      50.53, 50.54, 50.550000000000004, 50.56, 50.57, 50.58, 50.59, 50.6, 50.61,
      50.620000000000005, 50.63, 50.64, 50.65, 50.660000000000004, 50.67, 50.68,
      50.69, 50.7, 50.71, 50.72, 50.730000000000004, 50.74, 50.75, 50.76, 50.77,
      50.78, 50.79, 50.800000000000004, 50.81, 50.82, 50.83, 50.84, 50.85, 50.86,
      50.870000000000005, 50.88, 50.89, 50.9, 50.910000000000004, 50.92, 50.93,
      50.94, 50.95, 50.96, 50.97, 50.980000000000004, 50.99, 51.0, 51.01, 51.02,
      51.03, 51.04, 51.050000000000004, 51.06, 51.07, 51.08, 51.09, 51.1, 51.11,
      51.120000000000005, 51.13, 51.14, 51.15, 51.160000000000004, 51.17, 51.18,
      51.19, 51.2, 51.21, 51.22, 51.230000000000004, 51.24, 51.25, 51.26, 51.27,
      51.28, 51.29, 51.300000000000004, 51.31, 51.32, 51.33, 51.34, 51.35, 51.36,
      51.370000000000005, 51.38, 51.39, 51.4, 51.410000000000004, 51.42, 51.43,
      51.44, 51.45, 51.46, 51.47, 51.480000000000004, 51.49, 51.5, 51.51, 51.52,
      51.53, 51.54, 51.550000000000004, 51.56, 51.57, 51.58, 51.59, 51.6, 51.61,
      51.620000000000005, 51.63, 51.64, 51.65, 51.660000000000004, 51.67, 51.68,
      51.69, 51.7, 51.71, 51.72, 51.730000000000004, 51.74, 51.75, 51.76, 51.77,
      51.78, 51.79, 51.800000000000004, 51.81, 51.82, 51.83, 51.84, 51.85, 51.86,
      51.870000000000005, 51.88, 51.89, 51.9, 51.910000000000004, 51.92, 51.93,
      51.94, 51.95, 51.96, 51.97, 51.980000000000004, 51.99, 52.0, 52.01, 52.02,
      52.03, 52.04, 52.050000000000004, 52.06, 52.07, 52.08, 52.09, 52.1, 52.11,
      52.120000000000005, 52.13, 52.14, 52.15, 52.160000000000004, 52.17, 52.18,
      52.19, 52.2, 52.21, 52.22, 52.230000000000004, 52.24, 52.25, 52.26, 52.27,
      52.28, 52.29, 52.300000000000004, 52.31, 52.32, 52.33, 52.34, 52.35, 52.36,
      52.370000000000005, 52.38, 52.39, 52.4, 52.410000000000004, 52.42, 52.43,
      52.44, 52.45, 52.46, 52.47, 52.480000000000004, 52.49, 52.5, 52.51, 52.52,
      52.53, 52.54, 52.550000000000004, 52.56, 52.57, 52.58, 52.59, 52.6, 52.61,
      52.620000000000005, 52.63, 52.64, 52.65, 52.660000000000004, 52.67, 52.68,
      52.69, 52.7, 52.71, 52.72, 52.730000000000004, 52.74, 52.75, 52.76, 52.77,
      52.78, 52.79, 52.800000000000004, 52.81, 52.82, 52.83, 52.84, 52.85, 52.86,
      52.870000000000005, 52.88, 52.89, 52.9, 52.910000000000004, 52.92, 52.93,
      52.94, 52.95, 52.96, 52.97, 52.980000000000004, 52.99, 53.0, 53.01, 53.02,
      53.03, 53.04, 53.050000000000004, 53.06, 53.07, 53.08, 53.09, 53.1, 53.11,
      53.120000000000005, 53.13, 53.14, 53.15, 53.160000000000004, 53.17, 53.18,
      53.19, 53.2, 53.21, 53.22, 53.230000000000004, 53.24, 53.25, 53.26, 53.27,
      53.28, 53.29, 53.300000000000004, 53.31, 53.32, 53.33, 53.34, 53.35, 53.36,
      53.370000000000005, 53.38, 53.39, 53.4, 53.410000000000004, 53.42, 53.43,
      53.44, 53.45, 53.46, 53.47, 53.480000000000004, 53.49, 53.5, 53.51, 53.52,
      53.53, 53.54, 53.550000000000004, 53.56, 53.57, 53.58, 53.59, 53.6, 53.61,
      53.620000000000005, 53.63, 53.64, 53.65, 53.660000000000004, 53.67, 53.68,
      53.69, 53.7, 53.71, 53.72, 53.730000000000004, 53.74, 53.75, 53.76, 53.77,
      53.78, 53.79, 53.800000000000004, 53.81, 53.82, 53.83, 53.84, 53.85, 53.86,
      53.870000000000005, 53.88, 53.89, 53.9, 53.910000000000004, 53.92, 53.93,
      53.94, 53.95, 53.96, 53.97, 53.980000000000004, 53.99, 54.0, 54.01, 54.02,
      54.03, 54.04, 54.050000000000004, 54.06, 54.07, 54.08, 54.09, 54.1, 54.11,
      54.120000000000005, 54.13, 54.14, 54.15, 54.160000000000004, 54.17, 54.18,
      54.19, 54.2, 54.21, 54.22, 54.230000000000004, 54.24, 54.25, 54.26, 54.27,
      54.28, 54.29, 54.300000000000004, 54.31, 54.32, 54.33, 54.34, 54.35, 54.36,
      54.370000000000005, 54.38, 54.39, 54.4, 54.410000000000004, 54.42, 54.43,
      54.44, 54.45, 54.46, 54.47, 54.480000000000004, 54.49, 54.5, 54.51, 54.52,
      54.53, 54.54, 54.550000000000004, 54.56, 54.57, 54.58, 54.59, 54.6, 54.61,
      54.620000000000005, 54.63, 54.64, 54.65, 54.660000000000004, 54.67, 54.68,
      54.69, 54.7, 54.71, 54.72, 54.730000000000004, 54.74, 54.75, 54.76, 54.77,
      54.78, 54.79, 54.800000000000004, 54.81, 54.82, 54.83, 54.84, 54.85, 54.86,
      54.870000000000005, 54.88, 54.89, 54.9, 54.910000000000004, 54.92, 54.93,
      54.94, 54.95, 54.96, 54.97, 54.980000000000004, 54.99, 55.0, 55.01, 55.02,
      55.03, 55.04, 55.050000000000004, 55.06, 55.07, 55.08, 55.09, 55.1, 55.11,
      55.120000000000005, 55.13, 55.14, 55.15, 55.160000000000004, 55.17, 55.18,
      55.19, 55.2, 55.21, 55.22, 55.230000000000004, 55.24, 55.25, 55.26, 55.27,
      55.28, 55.29, 55.300000000000004, 55.31, 55.32, 55.33, 55.34, 55.35, 55.36,
      55.370000000000005, 55.38, 55.39, 55.4, 55.410000000000004, 55.42, 55.43,
      55.44, 55.45, 55.46, 55.47, 55.480000000000004, 55.49, 55.5, 55.51, 55.52,
      55.53, 55.54, 55.550000000000004, 55.56, 55.57, 55.58, 55.59, 55.6, 55.61,
      55.620000000000005, 55.63, 55.64, 55.65, 55.660000000000004, 55.67, 55.68,
      55.69, 55.7, 55.71, 55.72, 55.730000000000004, 55.74, 55.75, 55.76, 55.77,
      55.78, 55.79, 55.800000000000004, 55.81, 55.82, 55.83, 55.84, 55.85, 55.86,
      55.870000000000005, 55.88, 55.89, 55.9, 55.910000000000004, 55.92, 55.93,
      55.94, 55.95, 55.96, 55.97, 55.980000000000004, 55.99, 56.0, 56.01, 56.02,
      56.03, 56.04, 56.050000000000004, 56.06, 56.07, 56.08, 56.09, 56.1, 56.11,
      56.120000000000005, 56.13, 56.14, 56.15, 56.160000000000004, 56.17, 56.18,
      56.19, 56.2, 56.21, 56.22, 56.230000000000004, 56.24, 56.25, 56.26, 56.27,
      56.28, 56.29, 56.300000000000004, 56.31, 56.32, 56.33, 56.34, 56.35, 56.36,
      56.370000000000005, 56.38, 56.39, 56.4, 56.410000000000004, 56.42, 56.43,
      56.44, 56.45, 56.46, 56.47, 56.480000000000004, 56.49, 56.5, 56.51, 56.52,
      56.53, 56.54, 56.550000000000004, 56.56, 56.57, 56.58, 56.59, 56.6, 56.61,
      56.620000000000005, 56.63, 56.64, 56.65, 56.660000000000004, 56.67, 56.68,
      56.69, 56.7, 56.71, 56.72, 56.730000000000004, 56.74, 56.75, 56.76, 56.77,
      56.78, 56.79, 56.800000000000004, 56.81, 56.82, 56.83, 56.84, 56.85, 56.86,
      56.870000000000005, 56.88, 56.89, 56.9, 56.910000000000004, 56.92, 56.93,
      56.94, 56.95, 56.96, 56.97, 56.980000000000004, 56.99, 57.0, 57.01, 57.02,
      57.03, 57.04, 57.050000000000004, 57.06, 57.07, 57.08, 57.09, 57.1, 57.11,
      57.120000000000005, 57.13, 57.14, 57.15, 57.160000000000004, 57.17, 57.18,
      57.19, 57.2, 57.21, 57.22, 57.230000000000004, 57.24, 57.25, 57.26, 57.27,
      57.28, 57.29, 57.300000000000004, 57.31, 57.32, 57.33, 57.34, 57.35, 57.36,
      57.370000000000005, 57.38, 57.39, 57.4, 57.410000000000004, 57.42, 57.43,
      57.44, 57.45, 57.46, 57.47, 57.480000000000004, 57.49, 57.5, 57.51, 57.52,
      57.53, 57.54, 57.550000000000004, 57.56, 57.57, 57.58, 57.59, 57.6, 57.61,
      57.620000000000005, 57.63, 57.64, 57.65, 57.660000000000004, 57.67, 57.68,
      57.69, 57.7, 57.71, 57.72, 57.730000000000004, 57.74, 57.75, 57.76, 57.77,
      57.78, 57.79, 57.800000000000004, 57.81, 57.82, 57.83, 57.84, 57.85, 57.86,
      57.870000000000005, 57.88, 57.89, 57.9, 57.910000000000004, 57.92, 57.93,
      57.94, 57.95, 57.96, 57.97, 57.980000000000004, 57.99, 58.0, 58.01, 58.02,
      58.03, 58.04, 58.050000000000004, 58.06, 58.07, 58.08, 58.09, 58.1, 58.11,
      58.120000000000005, 58.13, 58.14, 58.15, 58.160000000000004, 58.17, 58.18,
      58.19, 58.2, 58.21, 58.22, 58.230000000000004, 58.24, 58.25, 58.26, 58.27,
      58.28, 58.29, 58.300000000000004, 58.31, 58.32, 58.33, 58.34, 58.35, 58.36,
      58.370000000000005, 58.38, 58.39, 58.4, 58.410000000000004, 58.42, 58.43,
      58.44, 58.45, 58.46, 58.47, 58.480000000000004, 58.49, 58.5, 58.51, 58.52,
      58.53, 58.54, 58.550000000000004, 58.56, 58.57, 58.58, 58.59, 58.6, 58.61,
      58.620000000000005, 58.63, 58.64, 58.65, 58.660000000000004, 58.67, 58.68,
      58.69, 58.7, 58.71, 58.72, 58.730000000000004, 58.74, 58.75, 58.76, 58.77,
      58.78, 58.79, 58.800000000000004, 58.81, 58.82, 58.83, 58.84, 58.85, 58.86,
      58.870000000000005, 58.88, 58.89, 58.9, 58.910000000000004, 58.92, 58.93,
      58.94, 58.95, 58.96, 58.97, 58.980000000000004, 58.99, 59.0, 59.01, 59.02,
      59.03, 59.04, 59.050000000000004, 59.06, 59.07, 59.08, 59.09, 59.1, 59.11,
      59.120000000000005, 59.13, 59.14, 59.15, 59.160000000000004, 59.17, 59.18,
      59.19, 59.2, 59.21, 59.22, 59.230000000000004, 59.24, 59.25, 59.26, 59.27,
      59.28, 59.29, 59.300000000000004, 59.31, 59.32, 59.33, 59.34, 59.35, 59.36,
      59.370000000000005, 59.38, 59.39, 59.4, 59.410000000000004, 59.42, 59.43,
      59.44, 59.45, 59.46, 59.47, 59.480000000000004, 59.49, 59.5, 59.51, 59.52,
      59.53, 59.54, 59.550000000000004, 59.56, 59.57, 59.58, 59.59, 59.6, 59.61,
      59.620000000000005, 59.63, 59.64, 59.65, 59.660000000000004, 59.67, 59.68,
      59.69, 59.7, 59.71, 59.72, 59.730000000000004, 59.74, 59.75, 59.76, 59.77,
      59.78, 59.79, 59.800000000000004, 59.81, 59.82, 59.83, 59.84, 59.85, 59.86,
      59.870000000000005, 59.88, 59.89, 59.9, 59.910000000000004, 59.92, 59.93,
      59.94, 59.95, 59.96, 59.97, 59.980000000000004, 59.99, 60.0, 60.01, 60.02,
      60.03, 60.04, 60.050000000000004, 60.06, 60.07, 60.08, 60.09, 60.1, 60.11,
      60.120000000000005, 60.13, 60.14, 60.15, 60.160000000000004, 60.17, 60.18,
      60.19, 60.2, 60.21, 60.22, 60.230000000000004, 60.24, 60.25, 60.26, 60.27,
      60.28, 60.29, 60.300000000000004, 60.31, 60.32, 60.33, 60.34, 60.35, 60.36,
      60.370000000000005, 60.38, 60.39, 60.4, 60.410000000000004, 60.42, 60.43,
      60.44, 60.45, 60.46, 60.47, 60.480000000000004, 60.49, 60.5, 60.51, 60.52,
      60.53, 60.54, 60.550000000000004, 60.56, 60.57, 60.58, 60.59, 60.6, 60.61,
      60.620000000000005, 60.63, 60.64, 60.65, 60.660000000000004, 60.67, 60.68,
      60.69, 60.7, 60.71, 60.72, 60.730000000000004, 60.74, 60.75, 60.76, 60.77,
      60.78, 60.79, 60.800000000000004, 60.81, 60.82, 60.83, 60.84, 60.85, 60.86,
      60.870000000000005, 60.88, 60.89, 60.9, 60.910000000000004, 60.92, 60.93,
      60.94, 60.95, 60.96, 60.97, 60.980000000000004, 60.99, 61.0, 61.01, 61.02,
      61.03, 61.04, 61.050000000000004, 61.06, 61.07, 61.08, 61.09, 61.1, 61.11,
      61.120000000000005, 61.13, 61.14, 61.15, 61.160000000000004, 61.17, 61.18,
      61.19, 61.2, 61.21, 61.22, 61.230000000000004, 61.24, 61.25, 61.26, 61.27,
      61.28, 61.29, 61.300000000000004, 61.31, 61.32, 61.33, 61.34, 61.35, 61.36,
      61.370000000000005, 61.38, 61.39, 61.4, 61.410000000000004, 61.42, 61.43,
      61.44, 61.45, 61.46, 61.47, 61.480000000000004, 61.49, 61.5, 61.51, 61.52,
      61.53, 61.54, 61.550000000000004, 61.56, 61.57, 61.58, 61.59, 61.6, 61.61,
      61.620000000000005, 61.63, 61.64, 61.65, 61.660000000000004, 61.67, 61.68,
      61.690000000000005, 61.7, 61.71, 61.72, 61.730000000000004, 61.74, 61.75,
      61.76, 61.77, 61.78, 61.79, 61.800000000000004, 61.81, 61.82, 61.83, 61.84,
      61.85, 61.86, 61.870000000000005, 61.88, 61.89, 61.9, 61.910000000000004,
      61.92, 61.93, 61.940000000000005, 61.95, 61.96, 61.97, 61.980000000000004,
      61.99, 62.0, 62.01, 62.02, 62.03, 62.04, 62.050000000000004, 62.06, 62.07,
      62.08, 62.09, 62.1, 62.11, 62.120000000000005, 62.13, 62.14, 62.15,
      62.160000000000004, 62.17, 62.18, 62.190000000000005, 62.2, 62.21, 62.22,
      62.230000000000004, 62.24, 62.25, 62.26, 62.27, 62.28, 62.29,
      62.300000000000004, 62.31, 62.32, 62.33, 62.34, 62.35, 62.36,
      62.370000000000005, 62.38, 62.39, 62.4, 62.410000000000004, 62.42, 62.43,
      62.440000000000005, 62.45, 62.46, 62.47, 62.480000000000004, 62.49, 62.5,
      62.51, 62.52, 62.53, 62.54, 62.550000000000004, 62.56, 62.57, 62.58, 62.59,
      62.6, 62.61, 62.620000000000005, 62.63, 62.64, 62.65, 62.660000000000004,
      62.67, 62.68, 62.690000000000005, 62.7, 62.71, 62.72, 62.730000000000004,
      62.74, 62.75, 62.76, 62.77, 62.78, 62.79, 62.800000000000004, 62.81, 62.82,
      62.83, 62.84, 62.85, 62.86, 62.870000000000005, 62.88, 62.89, 62.9,
      62.910000000000004, 62.92, 62.93, 62.940000000000005, 62.95, 62.96, 62.97,
      62.980000000000004, 62.99, 63.0, 63.01, 63.02, 63.03, 63.04,
      63.050000000000004, 63.06, 63.07, 63.08, 63.09, 63.1, 63.11,
      63.120000000000005, 63.13, 63.14, 63.15, 63.160000000000004, 63.17, 63.18,
      63.190000000000005, 63.2, 63.21, 63.22, 63.230000000000004, 63.24, 63.25,
      63.26, 63.27, 63.28, 63.29, 63.300000000000004, 63.31, 63.32, 63.33, 63.34,
      63.35, 63.36, 63.370000000000005, 63.38, 63.39, 63.4, 63.410000000000004,
      63.42, 63.43, 63.440000000000005, 63.45, 63.46, 63.47, 63.480000000000004,
      63.49, 63.5, 63.51, 63.52, 63.53, 63.54, 63.550000000000004, 63.56, 63.57,
      63.58, 63.59, 63.6, 63.61, 63.620000000000005, 63.63, 63.64, 63.65,
      63.660000000000004, 63.67, 63.68, 63.690000000000005, 63.7, 63.71, 63.72,
      63.730000000000004, 63.74, 63.75, 63.76, 63.77, 63.78, 63.79,
      63.800000000000004, 63.81, 63.82, 63.83, 63.84, 63.85, 63.86,
      63.870000000000005, 63.88, 63.89, 63.9, 63.910000000000004, 63.92, 63.93,
      63.940000000000005, 63.95, 63.96, 63.97, 63.980000000000004, 63.99, 64.0,
      64.01, 64.02, 64.03, 64.04, 64.05, 64.06, 64.070000000000007, 64.08, 64.09,
      64.1, 64.11, 64.12, 64.13, 64.14, 64.15, 64.16, 64.17, 64.18, 64.19, 64.2,
      64.210000000000008, 64.22, 64.23, 64.24, 64.25, 64.26, 64.27, 64.28, 64.29,
      64.3, 64.31, 64.320000000000007, 64.33, 64.34, 64.35, 64.36, 64.37, 64.38,
      64.39, 64.4, 64.41, 64.42, 64.43, 64.44, 64.45, 64.460000000000008, 64.47,
      64.48, 64.49, 64.5, 64.51, 64.52, 64.53, 64.54, 64.55, 64.56,
      64.570000000000007, 64.58, 64.59, 64.6, 64.61, 64.62, 64.63, 64.64, 64.65,
      64.66, 64.67, 64.68, 64.69, 64.7, 64.710000000000008, 64.72, 64.73, 64.74,
      64.75, 64.76, 64.77, 64.78, 64.79, 64.8, 64.81, 64.820000000000007, 64.83,
      64.84, 64.85, 64.86, 64.87, 64.88, 64.89, 64.9, 64.91, 64.92, 64.93, 64.94,
      64.95, 64.960000000000008, 64.97, 64.98, 64.99, 65.0, 65.01, 65.02, 65.03,
      65.04, 65.05, 65.06, 65.070000000000007, 65.08, 65.09, 65.1, 65.11, 65.12,
      65.13, 65.14, 65.15, 65.16, 65.17, 65.18, 65.19, 65.2, 65.210000000000008,
      65.22, 65.23, 65.24, 65.25, 65.26, 65.27, 65.28, 65.29, 65.3, 65.31,
      65.320000000000007, 65.33, 65.34, 65.35, 65.36, 65.37, 65.38, 65.39, 65.4,
      65.41, 65.42, 65.43, 65.44, 65.45, 65.460000000000008, 65.47, 65.48, 65.49,
      65.5, 65.51, 65.52, 65.53, 65.54, 65.55, 65.56, 65.570000000000007, 65.58,
      65.59, 65.6, 65.61, 65.62, 65.63, 65.64, 65.65, 65.66, 65.67, 65.68, 65.69,
      65.7, 65.710000000000008, 65.72, 65.73, 65.74, 65.75, 65.76, 65.77, 65.78,
      65.79, 65.8, 65.81, 65.820000000000007, 65.83, 65.84, 65.85, 65.86, 65.87,
      65.88, 65.89, 65.9, 65.91, 65.92, 65.93, 65.94, 65.95, 65.960000000000008,
      65.97, 65.98, 65.99, 66.0, 66.01, 66.02, 66.03, 66.04, 66.05, 66.06,
      66.070000000000007, 66.08, 66.09, 66.1, 66.11, 66.12, 66.13, 66.14, 66.15,
      66.16, 66.17, 66.18, 66.19, 66.2, 66.210000000000008, 66.22, 66.23, 66.24,
      66.25, 66.26, 66.27, 66.28, 66.29, 66.3, 66.31, 66.320000000000007, 66.33,
      66.34, 66.35, 66.36, 66.37, 66.38, 66.39, 66.4, 66.41, 66.42, 66.43, 66.44,
      66.45, 66.460000000000008, 66.47, 66.48, 66.49, 66.5, 66.51, 66.52, 66.53,
      66.54, 66.55, 66.56, 66.570000000000007, 66.58, 66.59, 66.6, 66.61, 66.62,
      66.63, 66.64, 66.65, 66.66, 66.67, 66.68, 66.69, 66.7, 66.710000000000008,
      66.72, 66.73, 66.74, 66.75, 66.76, 66.77, 66.78, 66.79, 66.8, 66.81,
      66.820000000000007, 66.83, 66.84, 66.85, 66.86, 66.87, 66.88, 66.89, 66.9,
      66.91, 66.92, 66.93, 66.94, 66.95, 66.960000000000008, 66.97, 66.98, 66.99,
      67.0, 67.01, 67.02, 67.03, 67.04, 67.05, 67.06, 67.070000000000007, 67.08,
      67.09, 67.1, 67.11, 67.12, 67.13, 67.14, 67.15, 67.16, 67.17, 67.18, 67.19,
      67.2, 67.210000000000008, 67.22, 67.23, 67.24, 67.25, 67.26, 67.27, 67.28,
      67.29, 67.3, 67.31, 67.320000000000007, 67.33, 67.34, 67.35, 67.36, 67.37,
      67.38, 67.39, 67.4, 67.41, 67.42, 67.43, 67.44, 67.45, 67.460000000000008,
      67.47, 67.48, 67.49, 67.5, 67.51, 67.52, 67.53, 67.54, 67.55, 67.56,
      67.570000000000007, 67.58, 67.59, 67.6, 67.61, 67.62, 67.63, 67.64, 67.65,
      67.66, 67.67, 67.68, 67.69, 67.7, 67.710000000000008, 67.72, 67.73, 67.74,
      67.75, 67.76, 67.77, 67.78, 67.79, 67.8, 67.81, 67.820000000000007, 67.83,
      67.84, 67.85, 67.86, 67.87, 67.88, 67.89, 67.9, 67.91, 67.92, 67.93, 67.94,
      67.95, 67.960000000000008, 67.97, 67.98, 67.99, 68.0, 68.01, 68.02, 68.03,
      68.04, 68.05, 68.06, 68.070000000000007, 68.08, 68.09, 68.1, 68.11, 68.12,
      68.13, 68.14, 68.15, 68.16, 68.17, 68.18, 68.19, 68.2, 68.210000000000008,
      68.22, 68.23, 68.24, 68.25, 68.26, 68.27, 68.28, 68.29, 68.3, 68.31,
      68.320000000000007, 68.33, 68.34, 68.350000000000009, 68.36, 68.37, 68.38,
      68.39, 68.4, 68.41, 68.42, 68.43, 68.44, 68.45, 68.460000000000008, 68.47,
      68.48, 68.49, 68.5, 68.51, 68.52, 68.53, 68.54, 68.55, 68.56,
      68.570000000000007, 68.58, 68.59, 68.600000000000009, 68.61, 68.62, 68.63,
      68.64, 68.65, 68.66, 68.67, 68.68, 68.69, 68.7, 68.710000000000008, 68.72,
      68.73, 68.74, 68.75, 68.76, 68.77, 68.78, 68.79, 68.8, 68.81,
      68.820000000000007, 68.83, 68.84, 68.850000000000009, 68.86, 68.87, 68.88,
      68.89, 68.9, 68.91, 68.92, 68.93, 68.94, 68.95, 68.960000000000008, 68.97,
      68.98, 68.99, 69.0, 69.01, 69.02, 69.03, 69.04, 69.05, 69.06,
      69.070000000000007, 69.08, 69.09, 69.100000000000009, 69.11, 69.12, 69.13,
      69.14, 69.15, 69.16, 69.17, 69.18, 69.19, 69.2, 69.210000000000008, 69.22,
      69.23, 69.24, 69.25, 69.26, 69.27, 69.28, 69.29, 69.3, 69.31,
      69.320000000000007, 69.33, 69.34, 69.350000000000009, 69.36, 69.37, 69.38,
      69.39, 69.4, 69.41, 69.42, 69.43, 69.44, 69.45, 69.460000000000008, 69.47,
      69.48, 69.49, 69.5, 69.51, 69.52, 69.53, 69.54, 69.55, 69.56,
      69.570000000000007, 69.58, 69.59, 69.600000000000009, 69.61, 69.62, 69.63,
      69.64, 69.65, 69.66, 69.67, 69.68, 69.69, 69.7, 69.710000000000008, 69.72,
      69.73, 69.74, 69.75, 69.76, 69.77, 69.78, 69.79, 69.8, 69.81,
      69.820000000000007, 69.83, 69.84, 69.850000000000009, 69.86, 69.87, 69.88,
      69.89, 69.9, 69.91, 69.92, 69.93, 69.94, 69.95, 69.960000000000008, 69.97,
      69.98, 69.99, 70.0, 70.01, 70.02, 70.03, 70.04, 70.05, 70.06,
      70.070000000000007, 70.08, 70.09, 70.100000000000009, 70.11, 70.12, 70.13,
      70.14, 70.15, 70.16, 70.17, 70.18, 70.19, 70.2, 70.210000000000008, 70.22,
      70.23, 70.24, 70.25, 70.26, 70.27, 70.28, 70.29, 70.3, 70.31,
      70.320000000000007, 70.33, 70.34, 70.350000000000009, 70.36, 70.37, 70.38,
      70.39, 70.4, 70.41, 70.42, 70.43, 70.44, 70.45, 70.460000000000008, 70.47,
      70.48, 70.49, 70.5, 70.51, 70.52, 70.53, 70.54, 70.55, 70.56,
      70.570000000000007, 70.58, 70.59, 70.600000000000009, 70.61, 70.62, 70.63,
      70.64, 70.65, 70.66, 70.67, 70.68, 70.69, 70.7, 70.710000000000008, 70.72,
      70.73, 70.74, 70.75, 70.76, 70.77, 70.78, 70.79, 70.8, 70.81,
      70.820000000000007, 70.83, 70.84, 70.850000000000009, 70.86, 70.87, 70.88,
      70.89, 70.9, 70.91, 70.92, 70.93, 70.94, 70.95, 70.960000000000008, 70.97,
      70.98, 70.99, 71.0, 71.01, 71.02, 71.03, 71.04, 71.05, 71.06,
      71.070000000000007, 71.08, 71.09, 71.100000000000009, 71.11, 71.12, 71.13,
      71.14, 71.15, 71.16, 71.17, 71.18, 71.19, 71.2, 71.210000000000008, 71.22,
      71.23, 71.24, 71.25, 71.26, 71.27, 71.28, 71.29, 71.3, 71.31,
      71.320000000000007, 71.33, 71.34, 71.350000000000009, 71.36, 71.37, 71.38,
      71.39, 71.4, 71.41, 71.42, 71.43, 71.44, 71.45, 71.460000000000008, 71.47,
      71.48, 71.49, 71.5, 71.51, 71.52, 71.53, 71.54, 71.55, 71.56,
      71.570000000000007, 71.58, 71.59, 71.600000000000009, 71.61, 71.62, 71.63,
      71.64, 71.65, 71.66, 71.67, 71.68, 71.69, 71.7, 71.710000000000008, 71.72,
      71.73, 71.74, 71.75, 71.76, 71.77, 71.78, 71.79, 71.8, 71.81,
      71.820000000000007, 71.83, 71.84, 71.850000000000009, 71.86, 71.87, 71.88,
      71.89, 71.9, 71.91, 71.92, 71.93, 71.94, 71.95, 71.960000000000008, 71.97,
      71.98, 71.99, 72.0, 72.01, 72.02, 72.03, 72.04, 72.05, 72.06,
      72.070000000000007, 72.08, 72.09, 72.100000000000009, 72.11, 72.12, 72.13,
      72.14, 72.15, 72.16, 72.17, 72.18, 72.19, 72.2, 72.210000000000008, 72.22,
      72.23, 72.24, 72.25, 72.26, 72.27, 72.28, 72.29, 72.3, 72.31,
      72.320000000000007, 72.33, 72.34, 72.350000000000009, 72.36, 72.37, 72.38,
      72.39, 72.4, 72.41, 72.42, 72.43, 72.44, 72.45, 72.460000000000008, 72.47,
      72.48, 72.49, 72.5, 72.51, 72.52, 72.53, 72.54, 72.55, 72.56,
      72.570000000000007, 72.58, 72.59, 72.600000000000009, 72.61, 72.62, 72.63,
      72.64, 72.65, 72.66, 72.67, 72.68, 72.69, 72.7, 72.710000000000008, 72.72,
      72.73, 72.74, 72.75, 72.76, 72.77, 72.78, 72.79, 72.8, 72.81,
      72.820000000000007, 72.83, 72.84, 72.850000000000009, 72.86, 72.87, 72.88,
      72.89, 72.9, 72.91, 72.92, 72.93, 72.94, 72.95, 72.960000000000008, 72.97,
      72.98, 72.99, 73.0, 73.01, 73.02, 73.03, 73.04, 73.05, 73.06,
      73.070000000000007, 73.08, 73.09, 73.100000000000009, 73.11, 73.12, 73.13,
      73.14, 73.15, 73.16, 73.17, 73.18, 73.19, 73.2, 73.210000000000008, 73.22,
      73.23, 73.24, 73.25, 73.26, 73.27, 73.28, 73.29, 73.3, 73.31,
      73.320000000000007, 73.33, 73.34, 73.350000000000009, 73.36, 73.37, 73.38,
      73.39, 73.4, 73.41, 73.42, 73.43, 73.44, 73.45, 73.460000000000008, 73.47,
      73.48, 73.49, 73.5, 73.51, 73.52, 73.53, 73.54, 73.55, 73.56,
      73.570000000000007, 73.58, 73.59, 73.600000000000009, 73.61, 73.62, 73.63,
      73.64, 73.65, 73.66, 73.67, 73.68, 73.69, 73.7, 73.710000000000008, 73.72,
      73.73, 73.74, 73.75, 73.76, 73.77, 73.78, 73.79, 73.8, 73.81,
      73.820000000000007, 73.83, 73.84, 73.850000000000009, 73.86, 73.87, 73.88,
      73.89, 73.9, 73.91, 73.92, 73.93, 73.94, 73.95, 73.960000000000008, 73.97,
      73.98, 73.99, 74.0, 74.01, 74.02, 74.03, 74.04, 74.05, 74.06,
      74.070000000000007, 74.08, 74.09, 74.100000000000009, 74.11, 74.12, 74.13,
      74.14, 74.15, 74.16, 74.17, 74.18, 74.19, 74.2, 74.210000000000008, 74.22,
      74.23, 74.24, 74.25, 74.26, 74.27, 74.28, 74.29, 74.3, 74.31,
      74.320000000000007, 74.33, 74.34, 74.350000000000009, 74.36, 74.37, 74.38,
      74.39, 74.4, 74.41, 74.42, 74.43, 74.44, 74.45, 74.460000000000008, 74.47,
      74.48, 74.49, 74.5, 74.51, 74.52, 74.53, 74.54, 74.55, 74.56,
      74.570000000000007, 74.58, 74.59, 74.600000000000009, 74.61, 74.62, 74.63,
      74.64, 74.65, 74.66, 74.67, 74.68, 74.69, 74.7, 74.710000000000008, 74.72,
      74.73, 74.74, 74.75, 74.76, 74.77, 74.78, 74.79, 74.8, 74.81,
      74.820000000000007, 74.83, 74.84, 74.850000000000009, 74.86, 74.87, 74.88,
      74.89, 74.9, 74.91, 74.92, 74.93, 74.94, 74.95, 74.960000000000008, 74.97,
      74.98, 74.99, 75.0, 75.01, 75.02, 75.03, 75.04, 75.05, 75.06,
      75.070000000000007, 75.08, 75.09, 75.100000000000009, 75.11, 75.12, 75.13,
      75.14, 75.15, 75.16, 75.17, 75.18, 75.19, 75.2, 75.210000000000008, 75.22,
      75.23, 75.24, 75.25, 75.26, 75.27, 75.28, 75.29, 75.3, 75.31,
      75.320000000000007, 75.33, 75.34, 75.350000000000009, 75.36, 75.37, 75.38,
      75.39, 75.4, 75.41, 75.42, 75.43, 75.44, 75.45, 75.460000000000008, 75.47,
      75.48, 75.49, 75.5, 75.51, 75.52, 75.53, 75.54, 75.55, 75.56,
      75.570000000000007, 75.58, 75.59, 75.600000000000009, 75.61, 75.62, 75.63,
      75.64, 75.65, 75.66, 75.67, 75.68, 75.69, 75.7, 75.710000000000008, 75.72,
      75.73, 75.74, 75.75, 75.76, 75.77, 75.78, 75.79, 75.8, 75.81,
      75.820000000000007, 75.83, 75.84, 75.850000000000009, 75.86, 75.87, 75.88,
      75.89, 75.9, 75.91, 75.92, 75.93, 75.94, 75.95, 75.960000000000008, 75.97,
      75.98, 75.99, 76.0, 76.01, 76.02, 76.03, 76.04, 76.05, 76.06,
      76.070000000000007, 76.08, 76.09, 76.100000000000009, 76.11, 76.12, 76.13,
      76.14, 76.15, 76.16, 76.17, 76.18, 76.19, 76.2, 76.210000000000008, 76.22,
      76.23, 76.24, 76.25, 76.26, 76.27, 76.28, 76.29, 76.3, 76.31,
      76.320000000000007, 76.33, 76.34, 76.350000000000009, 76.36, 76.37, 76.38,
      76.39, 76.4, 76.41, 76.42, 76.43, 76.44, 76.45, 76.460000000000008, 76.47,
      76.48, 76.49, 76.5, 76.51, 76.52, 76.53, 76.54, 76.55, 76.56,
      76.570000000000007, 76.58, 76.59, 76.600000000000009, 76.61, 76.62, 76.63,
      76.64, 76.65, 76.66, 76.67, 76.68, 76.69, 76.7, 76.710000000000008, 76.72,
      76.73, 76.74, 76.75, 76.76, 76.77, 76.78, 76.79, 76.8, 76.81,
      76.820000000000007, 76.83, 76.84, 76.850000000000009, 76.86, 76.87, 76.88,
      76.89, 76.9, 76.91, 76.92, 76.93, 76.94, 76.95, 76.960000000000008, 76.97,
      76.98, 76.99, 77.0, 77.01, 77.02, 77.03, 77.04, 77.05, 77.06,
      77.070000000000007, 77.08, 77.09, 77.100000000000009, 77.11, 77.12, 77.13,
      77.14, 77.15, 77.16, 77.17, 77.18, 77.19, 77.2, 77.210000000000008, 77.22,
      77.23, 77.24, 77.25, 77.26, 77.27, 77.28, 77.29, 77.3, 77.31,
      77.320000000000007, 77.33, 77.34, 77.350000000000009, 77.36, 77.37, 77.38,
      77.39, 77.4, 77.41, 77.42, 77.43, 77.44, 77.45, 77.460000000000008, 77.47,
      77.48, 77.49, 77.5, 77.51, 77.52, 77.53, 77.54, 77.55, 77.56,
      77.570000000000007, 77.58, 77.59, 77.600000000000009, 77.61, 77.62, 77.63,
      77.64, 77.65, 77.66, 77.67, 77.68, 77.69, 77.7, 77.710000000000008, 77.72,
      77.73, 77.74, 77.75, 77.76, 77.77, 77.78, 77.79, 77.8, 77.81,
      77.820000000000007, 77.83, 77.84, 77.850000000000009, 77.86, 77.87, 77.88,
      77.89, 77.9, 77.91, 77.92, 77.93, 77.94, 77.95, 77.960000000000008, 77.97,
      77.98, 77.99, 78.0, 78.01, 78.02, 78.03, 78.04, 78.05, 78.06,
      78.070000000000007, 78.08, 78.09, 78.100000000000009, 78.11, 78.12, 78.13,
      78.14, 78.15, 78.16, 78.17, 78.18, 78.19, 78.2, 78.210000000000008, 78.22,
      78.23, 78.24, 78.25, 78.26, 78.27, 78.28, 78.29, 78.3, 78.31,
      78.320000000000007, 78.33, 78.34, 78.350000000000009, 78.36, 78.37, 78.38,
      78.39, 78.4, 78.41, 78.42, 78.43, 78.44, 78.45, 78.460000000000008, 78.47,
      78.48, 78.49, 78.5, 78.51, 78.52, 78.53, 78.54, 78.55, 78.56,
      78.570000000000007, 78.58, 78.59, 78.600000000000009, 78.61, 78.62, 78.63,
      78.64, 78.65, 78.66, 78.67, 78.68, 78.69, 78.7, 78.710000000000008, 78.72,
      78.73, 78.74, 78.75, 78.76, 78.77, 78.78, 78.79, 78.8, 78.81,
      78.820000000000007, 78.83, 78.84, 78.850000000000009, 78.86, 78.87, 78.88,
      78.89, 78.9, 78.91, 78.92, 78.93, 78.94, 78.95, 78.960000000000008, 78.97,
      78.98, 78.99, 79.0, 79.01, 79.02, 79.03, 79.04, 79.05, 79.06,
      79.070000000000007, 79.08, 79.09, 79.100000000000009, 79.11, 79.12, 79.13,
      79.14, 79.15, 79.16, 79.17, 79.18, 79.19, 79.2, 79.210000000000008, 79.22,
      79.23, 79.24, 79.25, 79.26, 79.27, 79.28, 79.29, 79.3, 79.31,
      79.320000000000007, 79.33, 79.34, 79.350000000000009, 79.36, 79.37, 79.38,
      79.39, 79.4, 79.41, 79.42, 79.43, 79.44, 79.45, 79.460000000000008, 79.47,
      79.48, 79.49, 79.5, 79.51, 79.52, 79.53, 79.54, 79.55, 79.56,
      79.570000000000007, 79.58, 79.59, 79.600000000000009, 79.61, 79.62, 79.63,
      79.64, 79.65, 79.66, 79.67, 79.68, 79.69, 79.7, 79.710000000000008, 79.72,
      79.73, 79.74, 79.75, 79.76, 79.77, 79.78, 79.79, 79.8, 79.81,
      79.820000000000007, 79.83, 79.84, 79.850000000000009, 79.86, 79.87, 79.88,
      79.89, 79.9, 79.91, 79.92, 79.93, 79.94, 79.95, 79.960000000000008, 79.97,
      79.98, 79.99 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.050050050050050053, 0.10010010010010011, 0.15015015015015015,
      0.20020020020020021, 0.25025025025025027, 0.3003003003003003,
      0.35035035035035034, 0.40040040040040042, 0.45045045045045046,
      0.50050050050050054, 0.55055055055055058, 0.60060060060060061,
      0.65065065065065064, 0.70070070070070067, 0.75075075075075071,
      0.80080080080080085, 0.85085085085085088, 0.90090090090090091,
      0.950950950950951, 1.0010010010010011, 1.0510510510510511,
      1.1011011011011012, 1.1511511511511512, 1.2012012012012012,
      1.2512512512512513, 1.3013013013013013, 1.3513513513513513,
      1.4014014014014013, 1.4514514514514514, 1.5015015015015014,
      1.5515515515515514, 1.6016016016016017, 1.6516516516516517,
      1.7017017017017018, 1.7517517517517518, 1.8018018018018018,
      1.8518518518518519, 1.9019019019019019, 1.9519519519519519,
      2.0020020020020022, 2.0520520520520522, 2.1021021021021022,
      2.1521521521521523, 2.2022022022022023, 2.2522522522522523,
      2.3023023023023024, 2.3523523523523524, 2.4024024024024024,
      2.4524524524524525, 2.5025025025025025, 2.5525525525525525,
      2.6026026026026026, 2.6526526526526526, 2.7027027027027026,
      2.7527527527527527, 2.8028028028028027, 2.8528528528528527,
      2.9029029029029028, 2.9529529529529528, 3.0030030030030028,
      3.0530530530530529, 3.1031031031031029, 3.1531531531531534,
      3.2032032032032034, 3.2532532532532534, 3.3033033033033035,
      3.3533533533533535, 3.4034034034034035, 3.4534534534534536,
      3.5035035035035036, 3.5535535535535536, 3.6036036036036037,
      3.6536536536536537, 3.7037037037037037, 3.7537537537537538,
      3.8038038038038038, 3.8538538538538538, 3.9039039039039038,
      3.9539539539539539, 4.0040040040040044, 4.0540540540540544,
      4.1041041041041044, 4.1541541541541545, 4.2042042042042045,
      4.2542542542542545, 4.3043043043043046, 4.3543543543543546,
      4.4044044044044046, 4.4544544544544546, 4.5045045045045047,
      4.5545545545545547, 4.6046046046046047, 4.6546546546546548,
      4.7047047047047048, 4.7547547547547548, 4.8048048048048049,
      4.8548548548548549, 4.9049049049049049, 4.954954954954955,
      5.005005005005005, 5.055055055055055, 5.1051051051051051,
      5.1551551551551551, 5.2052052052052051, 5.2552552552552552,
      5.3053053053053052, 5.3553553553553552, 5.4054054054054053,
      5.4554554554554553, 5.5055055055055053, 5.5555555555555554,
      5.6056056056056054, 5.6556556556556554, 5.7057057057057055,
      5.7557557557557555, 5.8058058058058055, 5.8558558558558556,
      5.9059059059059056, 5.9559559559559556, 6.0060060060060056,
      6.0560560560560557, 6.1061061061061057, 6.1561561561561557,
      6.2062062062062058, 6.2562562562562567, 6.3063063063063067,
      6.3563563563563568, 6.4064064064064068, 6.4564564564564568,
      6.5065065065065069, 6.5565565565565569, 6.6066066066066069,
      6.656656656656657, 6.706706706706707, 6.756756756756757,
      6.8068068068068071, 6.8568568568568571, 6.9069069069069071,
      6.9569569569569571, 7.0070070070070072, 7.0570570570570572,
      7.1071071071071072, 7.1571571571571573, 7.2072072072072073,
      7.2572572572572573, 7.3073073073073074, 7.3573573573573574,
      7.4074074074074074, 7.4574574574574575, 7.5075075075075075,
      7.5575575575575575, 7.6076076076076076, 7.6576576576576576,
      7.7077077077077076, 7.7577577577577577, 7.8078078078078077,
      7.8578578578578577, 7.9079079079079078, 7.9579579579579578,
      8.0080080080080087, 8.0580580580580587, 8.1081081081081088,
      8.1581581581581588, 8.2082082082082088, 8.2582582582582589,
      8.3083083083083089, 8.3583583583583589, 8.408408408408409,
      8.458458458458459, 8.508508508508509, 8.5585585585585591,
      8.6086086086086091, 8.65865865865866, 8.70870870870871, 8.75875875875876,
      8.80880880880881, 8.85885885885886, 8.90890890890891, 8.95895895895896,
      9.00900900900901, 9.05905905905906, 9.10910910910911, 9.15915915915916,
      9.20920920920921, 9.25925925925926, 9.30930930930931, 9.35935935935936,
      9.40940940940941, 9.45945945945946, 9.50950950950951, 9.55955955955956,
      9.60960960960961, 9.65965965965966, 9.70970970970971, 9.75975975975976,
      9.80980980980981, 9.85985985985986, 9.90990990990991, 9.95995995995996,
      10.01001001001001, 10.06006006006006, 10.11011011011011, 10.16016016016016,
      10.21021021021021, 10.26026026026026, 10.31031031031031, 10.36036036036036,
      10.41041041041041, 10.46046046046046, 10.51051051051051, 10.56056056056056,
      10.61061061061061, 10.66066066066066, 10.71071071071071, 10.76076076076076,
      10.810810810810811, 10.860860860860861, 10.910910910910911,
      10.960960960960961, 11.011011011011011, 11.061061061061061,
      11.111111111111111, 11.161161161161161, 11.211211211211211,
      11.261261261261261, 11.311311311311311, 11.361361361361361,
      11.411411411411411, 11.461461461461461, 11.511511511511511,
      11.561561561561561, 11.611611611611611, 11.661661661661661,
      11.711711711711711, 11.761761761761761, 11.811811811811811,
      11.861861861861861, 11.911911911911911, 11.961961961961961,
      12.012012012012011, 12.062062062062061, 12.112112112112111,
      12.162162162162161, 12.212212212212211, 12.262262262262261,
      12.312312312312311, 12.362362362362362, 12.412412412412412,
      12.462462462462462, 12.512512512512513, 12.562562562562563,
      12.612612612612613, 12.662662662662663, 12.712712712712714,
      12.762762762762764, 12.812812812812814, 12.862862862862864,
      12.912912912912914, 12.962962962962964, 13.013013013013014,
      13.063063063063064, 13.113113113113114, 13.163163163163164,
      13.213213213213214, 13.263263263263264, 13.313313313313314,
      13.363363363363364, 13.413413413413414, 13.463463463463464,
      13.513513513513514, 13.563563563563564, 13.613613613613614,
      13.663663663663664, 13.713713713713714, 13.763763763763764,
      13.813813813813814, 13.863863863863864, 13.913913913913914,
      13.963963963963964, 14.014014014014014, 14.064064064064064,
      14.114114114114114, 14.164164164164164, 14.214214214214214,
      14.264264264264265, 14.314314314314315, 14.364364364364365,
      14.414414414414415, 14.464464464464465, 14.514514514514515,
      14.564564564564565, 14.614614614614615, 14.664664664664665,
      14.714714714714715, 14.764764764764765, 14.814814814814815,
      14.864864864864865, 14.914914914914915, 14.964964964964965,
      15.015015015015015, 15.065065065065065, 15.115115115115115,
      15.165165165165165, 15.215215215215215, 15.265265265265265,
      15.315315315315315, 15.365365365365365, 15.415415415415415,
      15.465465465465465, 15.515515515515515, 15.565565565565565,
      15.615615615615615, 15.665665665665665, 15.715715715715715,
      15.765765765765765, 15.815815815815816, 15.865865865865866,
      15.915915915915916, 15.965965965965966, 16.016016016016017,
      16.066066066066067, 16.116116116116117, 16.166166166166168,
      16.216216216216218, 16.266266266266268, 16.316316316316318,
      16.366366366366368, 16.416416416416418, 16.466466466466468,
      16.516516516516518, 16.566566566566568, 16.616616616616618,
      16.666666666666668, 16.716716716716718, 16.766766766766768,
      16.816816816816818, 16.866866866866868, 16.916916916916918,
      16.966966966966968, 17.017017017017018, 17.067067067067068,
      17.117117117117118, 17.167167167167168, 17.217217217217218,
      17.267267267267268, 17.317317317317318, 17.367367367367368,
      17.417417417417418, 17.467467467467468, 17.517517517517518,
      17.567567567567568, 17.617617617617618, 17.667667667667668,
      17.717717717717719, 17.767767767767769, 17.817817817817819,
      17.867867867867869, 17.917917917917919, 17.967967967967969,
      18.018018018018019, 18.068068068068069, 18.118118118118119,
      18.168168168168169, 18.218218218218219, 18.268268268268269,
      18.318318318318319, 18.368368368368369, 18.418418418418419,
      18.468468468468469, 18.518518518518519, 18.568568568568569,
      18.618618618618619, 18.668668668668669, 18.718718718718719,
      18.768768768768769, 18.818818818818819, 18.868868868868869,
      18.918918918918919, 18.968968968968969, 19.019019019019019,
      19.069069069069069, 19.119119119119119, 19.169169169169169,
      19.219219219219219, 19.26926926926927, 19.31931931931932,
      19.36936936936937, 19.41941941941942, 19.46946946946947, 19.51951951951952,
      19.56956956956957, 19.61961961961962, 19.66966966966967, 19.71971971971972,
      19.76976976976977, 19.81981981981982, 19.86986986986987, 19.91991991991992,
      19.96996996996997, 20.02002002002002, 20.07007007007007, 20.12012012012012,
      20.17017017017017, 20.22022022022022, 20.27027027027027, 20.32032032032032,
      20.37037037037037, 20.42042042042042, 20.47047047047047, 20.52052052052052,
      20.57057057057057, 20.62062062062062, 20.67067067067067, 20.72072072072072,
      20.77077077077077, 20.820820820820821, 20.870870870870871,
      20.920920920920921, 20.970970970970971, 21.021021021021021,
      21.071071071071071, 21.121121121121121, 21.171171171171171,
      21.221221221221221, 21.271271271271271, 21.321321321321321,
      21.371371371371371, 21.421421421421421, 21.471471471471471,
      21.521521521521521, 21.571571571571571, 21.621621621621621,
      21.671671671671671, 21.721721721721721, 21.771771771771771,
      21.821821821821821, 21.871871871871871, 21.921921921921921,
      21.971971971971971, 22.022022022022021, 22.072072072072071,
      22.122122122122121, 22.172172172172171, 22.222222222222221,
      22.272272272272271, 22.322322322322321, 22.372372372372372,
      22.422422422422422, 22.472472472472472, 22.522522522522522,
      22.572572572572572, 22.622622622622622, 22.672672672672672,
      22.722722722722722, 22.772772772772772, 22.822822822822822,
      22.872872872872872, 22.922922922922922, 22.972972972972972,
      23.023023023023022, 23.073073073073072, 23.123123123123122,
      23.173173173173172, 23.223223223223222, 23.273273273273272,
      23.323323323323322, 23.373373373373372, 23.423423423423422,
      23.473473473473472, 23.523523523523522, 23.573573573573572,
      23.623623623623622, 23.673673673673672, 23.723723723723722,
      23.773773773773772, 23.823823823823822, 23.873873873873872,
      23.923923923923923, 23.973973973973973, 24.024024024024023,
      24.074074074074073, 24.124124124124123, 24.174174174174173,
      24.224224224224223, 24.274274274274273, 24.324324324324323,
      24.374374374374373, 24.424424424424423, 24.474474474474473,
      24.524524524524523, 24.574574574574573, 24.624624624624623,
      24.674674674674673, 24.724724724724723, 24.774774774774773,
      24.824824824824823, 24.874874874874873, 24.924924924924923,
      24.974974974974973, 25.025025025025027, 25.075075075075077,
      25.125125125125127, 25.175175175175177, 25.225225225225227,
      25.275275275275277, 25.325325325325327, 25.375375375375377,
      25.425425425425427, 25.475475475475477, 25.525525525525527,
      25.575575575575577, 25.625625625625627, 25.675675675675677,
      25.725725725725727, 25.775775775775777, 25.825825825825827,
      25.875875875875877, 25.925925925925927, 25.975975975975977,
      26.026026026026027, 26.076076076076077, 26.126126126126128,
      26.176176176176178, 26.226226226226228, 26.276276276276278,
      26.326326326326328, 26.376376376376378, 26.426426426426428,
      26.476476476476478, 26.526526526526528, 26.576576576576578,
      26.626626626626628, 26.676676676676678, 26.726726726726728,
      26.776776776776778, 26.826826826826828, 26.876876876876878,
      26.926926926926928, 26.976976976976978, 27.027027027027028,
      27.077077077077078, 27.127127127127128, 27.177177177177178,
      27.227227227227228, 27.277277277277278, 27.327327327327328,
      27.377377377377378, 27.427427427427428, 27.477477477477478,
      27.527527527527528, 27.577577577577578, 27.627627627627628,
      27.677677677677679, 27.727727727727729, 27.777777777777779,
      27.827827827827829, 27.877877877877879, 27.927927927927929,
      27.977977977977979, 28.028028028028029, 28.078078078078079,
      28.128128128128129, 28.178178178178179, 28.228228228228229,
      28.278278278278279, 28.328328328328329, 28.378378378378379,
      28.428428428428429, 28.478478478478479, 28.528528528528529,
      28.578578578578579, 28.628628628628629, 28.678678678678679,
      28.728728728728729, 28.778778778778779, 28.828828828828829,
      28.878878878878879, 28.928928928928929, 28.978978978978979,
      29.029029029029029, 29.079079079079079, 29.129129129129129,
      29.179179179179179, 29.22922922922923, 29.27927927927928,
      29.32932932932933, 29.37937937937938, 29.42942942942943, 29.47947947947948,
      29.52952952952953, 29.57957957957958, 29.62962962962963, 29.67967967967968,
      29.72972972972973, 29.77977977977978, 29.82982982982983, 29.87987987987988,
      29.92992992992993, 29.97997997997998, 30.03003003003003, 30.08008008008008,
      30.13013013013013, 30.18018018018018, 30.23023023023023, 30.28028028028028,
      30.33033033033033, 30.38038038038038, 30.43043043043043, 30.48048048048048,
      30.53053053053053, 30.58058058058058, 30.63063063063063, 30.68068068068068,
      30.73073073073073, 30.780780780780781, 30.830830830830831,
      30.880880880880881, 30.930930930930931, 30.980980980980981,
      31.031031031031031, 31.081081081081081, 31.131131131131131,
      31.181181181181181, 31.231231231231231, 31.281281281281281,
      31.331331331331331, 31.381381381381381, 31.431431431431431,
      31.481481481481481, 31.531531531531531, 31.581581581581581,
      31.631631631631631, 31.681681681681681, 31.731731731731731,
      31.781781781781781, 31.831831831831831, 31.881881881881881,
      31.931931931931931, 31.981981981981981, 32.032032032032035,
      32.082082082082081, 32.132132132132135, 32.182182182182181,
      32.232232232232235, 32.282282282282281, 32.332332332332335,
      32.382382382382382, 32.432432432432435, 32.482482482482482,
      32.532532532532535, 32.582582582582582, 32.632632632632635,
      32.682682682682682, 32.732732732732735, 32.782782782782782,
      32.832832832832835, 32.882882882882882, 32.932932932932935,
      32.982982982982982, 33.033033033033036, 33.083083083083082,
      33.133133133133136, 33.183183183183182, 33.233233233233236,
      33.283283283283282, 33.333333333333336, 33.383383383383382,
      33.433433433433436, 33.483483483483482, 33.533533533533536,
      33.583583583583582, 33.633633633633636, 33.683683683683682,
      33.733733733733736, 33.783783783783782, 33.833833833833836,
      33.883883883883883, 33.933933933933936, 33.983983983983983,
      34.034034034034036, 34.084084084084083, 34.134134134134136,
      34.184184184184183, 34.234234234234236, 34.284284284284283,
      34.334334334334336, 34.384384384384383, 34.434434434434436,
      34.484484484484483, 34.534534534534536, 34.584584584584583,
      34.634634634634637, 34.684684684684683, 34.734734734734737,
      34.784784784784783, 34.834834834834837, 34.884884884884883,
      34.934934934934937, 34.984984984984983, 35.035035035035037,
      35.085085085085083, 35.135135135135137, 35.185185185185183,
      35.235235235235237, 35.285285285285283, 35.335335335335337,
      35.385385385385383, 35.435435435435437, 35.485485485485484,
      35.535535535535537, 35.585585585585584, 35.635635635635637,
      35.685685685685684, 35.735735735735737, 35.785785785785784,
      35.835835835835837, 35.885885885885884, 35.935935935935937,
      35.985985985985984, 36.036036036036037, 36.086086086086084,
      36.136136136136138, 36.186186186186184, 36.236236236236238,
      36.286286286286284, 36.336336336336338, 36.386386386386384,
      36.436436436436438, 36.486486486486484, 36.536536536536538,
      36.586586586586584, 36.636636636636638, 36.686686686686684,
      36.736736736736738, 36.786786786786784, 36.836836836836838,
      36.886886886886884, 36.936936936936938, 36.986986986986985,
      37.037037037037038, 37.087087087087085, 37.137137137137138,
      37.187187187187185, 37.237237237237238, 37.287287287287285,
      37.337337337337338, 37.387387387387385, 37.437437437437438,
      37.487487487487485, 37.537537537537538, 37.587587587587585,
      37.637637637637638, 37.687687687687685, 37.737737737737739,
      37.787787787787785, 37.837837837837839, 37.887887887887885,
      37.937937937937939, 37.987987987987985, 38.038038038038039,
      38.088088088088085, 38.138138138138139, 38.188188188188185,
      38.238238238238239, 38.288288288288285, 38.338338338338339,
      38.388388388388385, 38.438438438438439, 38.488488488488485,
      38.538538538538539, 38.588588588588586, 38.638638638638639,
      38.688688688688686, 38.738738738738739, 38.788788788788786,
      38.838838838838839, 38.888888888888886, 38.938938938938939,
      38.988988988988986, 39.039039039039039, 39.089089089089086,
      39.139139139139139, 39.189189189189186, 39.23923923923924,
      39.289289289289286, 39.33933933933934, 39.389389389389386,
      39.43943943943944, 39.489489489489486, 39.53953953953954,
      39.589589589589586, 39.63963963963964, 39.689689689689686,
      39.73973973973974, 39.789789789789786, 39.83983983983984,
      39.889889889889886, 39.93993993993994, 39.989989989989986,
      40.04004004004004, 40.090090090090094, 40.14014014014014,
      40.190190190190194, 40.24024024024024, 40.290290290290294,
      40.34034034034034, 40.390390390390394, 40.44044044044044,
      40.490490490490494, 40.54054054054054, 40.590590590590594,
      40.64064064064064, 40.690690690690694, 40.74074074074074,
      40.790790790790794, 40.840840840840841, 40.890890890890894,
      40.940940940940941, 40.990990990990994, 41.041041041041041,
      41.091091091091094, 41.141141141141141, 41.191191191191194,
      41.241241241241241, 41.291291291291294, 41.341341341341341,
      41.391391391391394, 41.441441441441441, 41.491491491491495,
      41.541541541541541, 41.591591591591595, 41.641641641641641,
      41.691691691691695, 41.741741741741741, 41.791791791791795,
      41.841841841841841, 41.891891891891895, 41.941941941941941,
      41.991991991991995, 42.042042042042041, 42.092092092092095,
      42.142142142142141, 42.192192192192195, 42.242242242242241,
      42.292292292292295, 42.342342342342342, 42.392392392392395,
      42.442442442442442, 42.492492492492495, 42.542542542542542,
      42.592592592592595, 42.642642642642642, 42.692692692692695,
      42.742742742742742, 42.792792792792795, 42.842842842842842,
      42.892892892892895, 42.942942942942942, 42.992992992992995,
      43.043043043043042, 43.093093093093096, 43.143143143143142,
      43.193193193193196, 43.243243243243242, 43.293293293293296,
      43.343343343343342, 43.393393393393396, 43.443443443443442,
      43.493493493493496, 43.543543543543542, 43.593593593593596,
      43.643643643643642, 43.693693693693696, 43.743743743743742,
      43.793793793793796, 43.843843843843842, 43.893893893893896,
      43.943943943943943, 43.993993993993996, 44.044044044044043,
      44.094094094094096, 44.144144144144143, 44.194194194194196,
      44.244244244244243, 44.294294294294296, 44.344344344344343,
      44.394394394394396, 44.444444444444443, 44.4944944944945,
      44.544544544544543, 44.5945945945946, 44.644644644644643, 44.6946946946947,
      44.744744744744743, 44.7947947947948, 44.844844844844843, 44.8948948948949,
      44.944944944944943, 44.994994994995, 45.045045045045043, 45.0950950950951,
      45.145145145145143, 45.1951951951952, 45.245245245245243, 45.2952952952953,
      45.345345345345343, 45.3953953953954, 45.445445445445444, 45.4954954954955,
      45.545545545545544, 45.5955955955956, 45.645645645645644, 45.6956956956957,
      45.745745745745744, 45.7957957957958, 45.845845845845844, 45.8958958958959,
      45.945945945945944, 45.995995995996, 46.046046046046044, 46.0960960960961,
      46.146146146146144, 46.1961961961962, 46.246246246246244, 46.2962962962963,
      46.346346346346344, 46.3963963963964, 46.446446446446444, 46.4964964964965,
      46.546546546546544, 46.5965965965966, 46.646646646646644, 46.6966966966967,
      46.746746746746744, 46.7967967967968, 46.846846846846844, 46.8968968968969,
      46.946946946946944, 46.996996996997, 47.047047047047045, 47.0970970970971,
      47.147147147147145, 47.1971971971972, 47.247247247247245, 47.2972972972973,
      47.347347347347345, 47.3973973973974, 47.447447447447445, 47.4974974974975,
      47.547547547547545, 47.5975975975976, 47.647647647647645, 47.6976976976977,
      47.747747747747745, 47.7977977977978, 47.847847847847845, 47.8978978978979,
      47.947947947947945, 47.997997997998, 48.048048048048045, 48.0980980980981,
      48.148148148148145, 48.1981981981982, 48.248248248248245, 48.2982982982983,
      48.348348348348345, 48.3983983983984, 48.448448448448445, 48.4984984984985,
      48.548548548548546, 48.5985985985986, 48.648648648648646, 48.6986986986987,
      48.748748748748746, 48.7987987987988, 48.848848848848846, 48.8988988988989,
      48.948948948948946, 48.998998998999, 49.049049049049046, 49.0990990990991,
      49.149149149149146, 49.1991991991992, 49.249249249249246, 49.2992992992993,
      49.349349349349346, 49.3993993993994, 49.449449449449446, 49.4994994994995,
      49.549549549549546, 49.5995995995996, 49.649649649649646, 49.6996996996997,
      49.749749749749746, 49.7997997997998, 49.849849849849846, 49.8998998998999,
      49.949949949949946, 50.0, 50.0, 49.8998998998999, 49.7997997997998,
      49.6996996996997, 49.5995995995996, 49.4994994994995, 49.3993993993994,
      49.2992992992993, 49.1991991991992, 49.0990990990991, 48.998998998999,
      48.8988988988989, 48.7987987987988, 48.6986986986987, 48.5985985985986,
      48.4984984984985, 48.3983983983984, 48.2982982982983, 48.1981981981982,
      48.0980980980981, 47.997997997998, 47.8978978978979, 47.7977977977978,
      47.6976976976977, 47.5975975975976, 47.4974974974975, 47.3973973973974,
      47.2972972972973, 47.1971971971972, 47.0970970970971, 46.996996996997,
      46.8968968968969, 46.7967967967968, 46.6966966966967, 46.5965965965966,
      46.4964964964965, 46.3963963963964, 46.2962962962963, 46.1961961961962,
      46.0960960960961, 45.995995995996, 45.8958958958959, 45.7957957957958,
      45.6956956956957, 45.5955955955956, 45.4954954954955, 45.3953953953954,
      45.2952952952953, 45.1951951951952, 45.0950950950951, 44.994994994995,
      44.8948948948949, 44.7947947947948, 44.6946946946947, 44.5945945945946,
      44.4944944944945, 44.394394394394396, 44.294294294294296,
      44.194194194194196, 44.094094094094096, 43.993993993993996,
      43.893893893893896, 43.793793793793796, 43.693693693693696,
      43.593593593593596, 43.493493493493496, 43.393393393393396,
      43.293293293293296, 43.193193193193196, 43.093093093093096,
      42.992992992992995, 42.892892892892895, 42.792792792792795,
      42.692692692692695, 42.592592592592595, 42.492492492492495,
      42.392392392392395, 42.292292292292295, 42.192192192192195,
      42.092092092092095, 41.991991991991995, 41.891891891891895,
      41.791791791791795, 41.691691691691695, 41.591591591591595,
      41.491491491491495, 41.391391391391394, 41.291291291291294,
      41.191191191191194, 41.091091091091094, 40.990990990990994,
      40.890890890890894, 40.790790790790794, 40.690690690690694,
      40.590590590590594, 40.490490490490494, 40.390390390390394,
      40.290290290290294, 40.190190190190194, 40.090090090090094,
      39.989989989989994, 39.889889889889893, 39.789789789789793,
      39.689689689689693, 39.589589589589593, 39.489489489489493,
      39.389389389389393, 39.289289289289293, 39.189189189189193,
      39.089089089089093, 38.988988988988993, 38.888888888888893,
      38.788788788788793, 38.688688688688693, 38.588588588588593,
      38.488488488488493, 38.388388388388393, 38.288288288288292,
      38.188188188188192, 38.088088088088092, 37.987987987987992,
      37.887887887887892, 37.787787787787792, 37.687687687687692,
      37.587587587587592, 37.487487487487492, 37.387387387387392,
      37.287287287287292, 37.187187187187192, 37.087087087087092,
      36.986986986986992, 36.886886886886892, 36.786786786786791,
      36.686686686686691, 36.586586586586591, 36.486486486486491,
      36.386386386386391, 36.286286286286291, 36.186186186186191,
      36.086086086086091, 35.985985985985991, 35.885885885885891,
      35.785785785785791, 35.685685685685691, 35.585585585585591,
      35.485485485485491, 35.385385385385391, 35.285285285285291,
      35.18518518518519, 35.08508508508509, 34.98498498498499, 34.88488488488489,
      34.78478478478479, 34.68468468468469, 34.58458458458459, 34.48448448448449,
      34.38438438438439, 34.28428428428429, 34.18418418418419,
      34.084084084084083, 33.983983983983983, 33.883883883883883,
      33.783783783783782, 33.683683683683682, 33.583583583583582,
      33.483483483483482, 33.383383383383382, 33.283283283283282,
      33.183183183183182, 33.083083083083082, 32.982982982982982,
      32.882882882882882, 32.782782782782782, 32.682682682682682,
      32.582582582582582, 32.482482482482482, 32.382382382382382,
      32.282282282282281, 32.182182182182181, 32.082082082082081,
      31.981981981981985, 31.881881881881885, 31.781781781781785,
      31.681681681681685, 31.581581581581585, 31.481481481481485,
      31.381381381381384, 31.281281281281284, 31.181181181181184,
      31.081081081081084, 30.980980980980984, 30.880880880880884,
      30.780780780780784, 30.680680680680684, 30.580580580580584,
      30.480480480480484, 30.380380380380384, 30.280280280280284,
      30.180180180180184, 30.080080080080084, 29.979979979979984,
      29.87987987987988, 29.77977977977978, 29.67967967967968, 29.57957957957958,
      29.47947947947948, 29.37937937937938, 29.27927927927928,
      29.179179179179179, 29.079079079079079, 28.978978978978979,
      28.878878878878879, 28.778778778778779, 28.678678678678679,
      28.578578578578579, 28.478478478478479, 28.378378378378379,
      28.278278278278279, 28.178178178178179, 28.078078078078079,
      27.977977977977979, 27.877877877877879, 27.777777777777779,
      27.677677677677679, 27.577577577577578, 27.477477477477478,
      27.377377377377378, 27.277277277277278, 27.177177177177178,
      27.077077077077078, 26.976976976976978, 26.876876876876878,
      26.776776776776778, 26.676676676676678, 26.576576576576578,
      26.476476476476478, 26.376376376376378, 26.276276276276278,
      26.176176176176178, 26.076076076076077, 25.975975975975977,
      25.875875875875877, 25.775775775775777, 25.675675675675677,
      25.575575575575577, 25.475475475475477, 25.375375375375377,
      25.275275275275277, 25.175175175175177, 25.075075075075077,
      24.974974974974977, 24.874874874874877, 24.774774774774777,
      24.674674674674677, 24.574574574574577, 24.474474474474476,
      24.374374374374376, 24.274274274274276, 24.174174174174176,
      24.074074074074076, 23.973973973973976, 23.873873873873876,
      23.773773773773776, 23.673673673673676, 23.573573573573576,
      23.473473473473476, 23.373373373373376, 23.273273273273276,
      23.173173173173176, 23.073073073073076, 22.972972972972975,
      22.872872872872875, 22.772772772772775, 22.672672672672675,
      22.572572572572575, 22.472472472472475, 22.372372372372375,
      22.272272272272275, 22.172172172172175, 22.072072072072075,
      21.971971971971975, 21.871871871871875, 21.771771771771775,
      21.671671671671675, 21.571571571571575, 21.471471471471475,
      21.371371371371374, 21.271271271271271, 21.171171171171171,
      21.071071071071071, 20.970970970970971, 20.870870870870871,
      20.77077077077077, 20.67067067067067, 20.57057057057057, 20.47047047047047,
      20.37037037037037, 20.27027027027027, 20.17017017017017, 20.07007007007007,
      19.96996996996997, 19.86986986986987, 19.76976976976977, 19.66966966966967,
      19.56956956956957, 19.46946946946947, 19.36936936936937, 19.26926926926927,
      19.169169169169169, 19.069069069069069, 18.968968968968969,
      18.868868868868869, 18.768768768768769, 18.668668668668669,
      18.568568568568569, 18.468468468468469, 18.368368368368369,
      18.268268268268269, 18.168168168168169, 18.068068068068069,
      17.967967967967969, 17.867867867867869, 17.767767767767769,
      17.667667667667668, 17.567567567567568, 17.467467467467468,
      17.367367367367368, 17.267267267267268, 17.167167167167168,
      17.067067067067068, 16.966966966966968, 16.866866866866868,
      16.766766766766768, 16.666666666666668, 16.566566566566568,
      16.466466466466468, 16.366366366366368, 16.266266266266268,
      16.166166166166168, 16.066066066066067, 15.965965965965967,
      15.865865865865867, 15.765765765765767, 15.665665665665667,
      15.565565565565567, 15.465465465465467, 15.365365365365367,
      15.265265265265267, 15.165165165165167, 15.065065065065067,
      14.964964964964967, 14.864864864864865, 14.764764764764765,
      14.664664664664665, 14.564564564564565, 14.464464464464465,
      14.364364364364365, 14.264264264264265, 14.164164164164164,
      14.064064064064064, 13.963963963963964, 13.863863863863864,
      13.763763763763764, 13.663663663663664, 13.563563563563564,
      13.463463463463464, 13.363363363363364, 13.263263263263264,
      13.163163163163164, 13.063063063063064, 12.962962962962964,
      12.862862862862864, 12.762762762762764, 12.662662662662663,
      12.562562562562563, 12.462462462462463, 12.362362362362363,
      12.262262262262263, 12.162162162162163, 12.062062062062063,
      11.961961961961963, 11.861861861861863, 11.761761761761763,
      11.661661661661663, 11.561561561561563, 11.461461461461463,
      11.361361361361363, 11.261261261261263, 11.161161161161163,
      11.061061061061062, 10.960960960960962, 10.860860860860862,
      10.760760760760762, 10.66066066066066, 10.56056056056056,
      10.46046046046046, 10.36036036036036, 10.26026026026026, 10.16016016016016,
      10.06006006006006, 9.95995995995996, 9.85985985985986, 9.75975975975976,
      9.65965965965966, 9.55955955955956, 9.45945945945946, 9.35935935935936,
      9.25925925925926, 9.15915915915916, 9.05905905905906, 8.95895895895896,
      8.85885885885886, 8.75875875875876, 8.65865865865866, 8.5585585585585591,
      8.458458458458459, 8.3583583583583589, 8.2582582582582589,
      8.1581581581581588, 8.0580580580580587, 7.9579579579579587,
      7.8578578578578586, 7.7577577577577586, 7.6576576576576585,
      7.5575575575575584, 7.4574574574574575, 7.3573573573573574,
      7.2572572572572573, 7.1571571571571573, 7.0570570570570572,
      6.9569569569569571, 6.8568568568568571, 6.756756756756757,
      6.656656656656657, 6.5565565565565569, 6.4564564564564568,
      6.3563563563563568, 6.2562562562562567, 6.1561561561561566,
      6.0560560560560566, 5.9559559559559565, 5.8558558558558564,
      5.7557557557557564, 5.6556556556556563, 5.5555555555555562,
      5.4554554554554562, 5.3553553553553561, 5.2552552552552552,
      5.1551551551551551, 5.055055055055055, 4.954954954954955,
      4.8548548548548549, 4.7547547547547548, 4.6546546546546548,
      4.5545545545545547, 4.4544544544544546, 4.3543543543543546,
      4.2542542542542545, 4.1541541541541545, 4.0540540540540544,
      3.9539539539539543, 3.8538538538538543, 3.7537537537537542,
      3.6536536536536537, 3.5535535535535536, 3.4534534534534536,
      3.3533533533533535, 3.2532532532532534, 3.1531531531531534,
      3.0530530530530533, 2.9529529529529532, 2.8528528528528532,
      2.7527527527527531, 2.6526526526526526, 2.5525525525525525,
      2.4524524524524525, 2.3523523523523524, 2.2522522522522523,
      2.1521521521521523, 2.0520520520520522, 1.9519519519519521,
      1.8518518518518519, 1.7517517517517518, 1.6516516516516517,
      1.5515515515515517, 1.4514514514514516, 1.3513513513513515,
      1.2512512512512513, 1.1511511511511512, 1.0510510510510511,
      0.95095095095095106, 0.85085085085085088, 0.75075075075075082,
      0.65065065065065064, 0.55055055055055058, 0.45045045045045046,
      0.35035035035035039, 0.25025025025025027, 0.15015015015015015,
      0.050050050050050053, -0.050050050050050053, -0.15015015015015015,
      -0.25025025025025027, -0.35035035035035039, -0.45045045045045046,
      -0.55055055055055058, -0.65065065065065064, -0.75075075075075082,
      -0.85085085085085088, -0.95095095095095106, -1.0510510510510511,
      -1.1511511511511512, -1.2512512512512513, -1.3513513513513515,
      -1.4514514514514516, -1.5515515515515517, -1.6516516516516517,
      -1.7517517517517518, -1.8518518518518519, -1.9519519519519521,
      -2.0520520520520522, -2.1521521521521523, -2.2522522522522523,
      -2.3523523523523524, -2.4524524524524525, -2.5525525525525525,
      -2.6526526526526526, -2.7527527527527531, -2.8528528528528532,
      -2.9529529529529532, -3.0530530530530533, -3.1531531531531534,
      -3.2532532532532534, -3.3533533533533535, -3.4534534534534536,
      -3.5535535535535536, -3.6536536536536537, -3.7537537537537542,
      -3.8538538538538543, -3.9539539539539543, -4.0540540540540544,
      -4.1541541541541545, -4.2542542542542545, -4.3543543543543546,
      -4.4544544544544546, -4.5545545545545547, -4.6546546546546548,
      -4.7547547547547548, -4.8548548548548549, -4.954954954954955,
      -5.055055055055055, -5.1551551551551551, -5.2552552552552552,
      -5.3553553553553561, -5.4554554554554562, -5.5555555555555562,
      -5.6556556556556563, -5.7557557557557564, -5.8558558558558564,
      -5.9559559559559565, -6.0560560560560566, -6.1561561561561566,
      -6.2562562562562567, -6.3563563563563568, -6.4564564564564568,
      -6.5565565565565569, -6.656656656656657, -6.756756756756757,
      -6.8568568568568571, -6.9569569569569571, -7.0570570570570572,
      -7.1571571571571573, -7.2572572572572573, -7.3573573573573574,
      -7.4574574574574575, -7.5575575575575584, -7.6576576576576585,
      -7.7577577577577586, -7.8578578578578586, -7.9579579579579587,
      -8.0580580580580587, -8.1581581581581588, -8.2582582582582589,
      -8.3583583583583589, -8.458458458458459, -8.5585585585585591,
      -8.65865865865866, -8.75875875875876, -8.85885885885886, -8.95895895895896,
      -9.05905905905906, -9.15915915915916, -9.25925925925926, -9.35935935935936,
      -9.45945945945946, -9.55955955955956, -9.65965965965966, -9.75975975975976,
      -9.85985985985986, -9.95995995995996, -10.06006006006006,
      -10.16016016016016, -10.26026026026026, -10.36036036036036,
      -10.46046046046046, -10.56056056056056, -10.66066066066066,
      -10.760760760760762, -10.860860860860862, -10.960960960960962,
      -11.061061061061062, -11.161161161161163, -11.261261261261263,
      -11.361361361361363, -11.461461461461463, -11.561561561561563,
      -11.661661661661663, -11.761761761761763, -11.861861861861863,
      -11.961961961961963, -12.062062062062063, -12.162162162162163,
      -12.262262262262263, -12.362362362362363, -12.462462462462463,
      -12.562562562562563, -12.662662662662663, -12.762762762762764,
      -12.862862862862864, -12.962962962962964, -13.063063063063064,
      -13.163163163163164, -13.263263263263264, -13.363363363363364,
      -13.463463463463464, -13.563563563563564, -13.663663663663664,
      -13.763763763763764, -13.863863863863864, -13.963963963963964,
      -14.064064064064064, -14.164164164164164, -14.264264264264265,
      -14.364364364364365, -14.464464464464465, -14.564564564564565,
      -14.664664664664665, -14.764764764764765, -14.864864864864865,
      -14.964964964964967, -15.065065065065067, -15.165165165165167,
      -15.265265265265267, -15.365365365365367, -15.465465465465467,
      -15.565565565565567, -15.665665665665667, -15.765765765765767,
      -15.865865865865867, -15.965965965965967, -16.066066066066067,
      -16.166166166166168, -16.266266266266268, -16.366366366366368,
      -16.466466466466468, -16.566566566566568, -16.666666666666668,
      -16.766766766766768, -16.866866866866868, -16.966966966966968,
      -17.067067067067068, -17.167167167167168, -17.267267267267268,
      -17.367367367367368, -17.467467467467468, -17.567567567567568,
      -17.667667667667668, -17.767767767767769, -17.867867867867869,
      -17.967967967967969, -18.068068068068069, -18.168168168168169,
      -18.268268268268269, -18.368368368368369, -18.468468468468469,
      -18.568568568568569, -18.668668668668669, -18.768768768768769,
      -18.868868868868869, -18.968968968968969, -19.069069069069069,
      -19.169169169169169, -19.26926926926927, -19.36936936936937,
      -19.46946946946947, -19.56956956956957, -19.66966966966967,
      -19.76976976976977, -19.86986986986987, -19.96996996996997,
      -20.07007007007007, -20.17017017017017, -20.27027027027027,
      -20.37037037037037, -20.47047047047047, -20.57057057057057,
      -20.67067067067067, -20.77077077077077, -20.870870870870871,
      -20.970970970970971, -21.071071071071071, -21.171171171171171,
      -21.271271271271271, -21.371371371371374, -21.471471471471475,
      -21.571571571571575, -21.671671671671675, -21.771771771771775,
      -21.871871871871875, -21.971971971971975, -22.072072072072075,
      -22.172172172172175, -22.272272272272275, -22.372372372372375,
      -22.472472472472475, -22.572572572572575, -22.672672672672675,
      -22.772772772772775, -22.872872872872875, -22.972972972972975,
      -23.073073073073076, -23.173173173173176, -23.273273273273276,
      -23.373373373373376, -23.473473473473476, -23.573573573573576,
      -23.673673673673676, -23.773773773773776, -23.873873873873876,
      -23.973973973973976, -24.074074074074076, -24.174174174174176,
      -24.274274274274276, -24.374374374374376, -24.474474474474476,
      -24.574574574574577, -24.674674674674677, -24.774774774774777,
      -24.874874874874877, -24.974974974974977, -25.075075075075077,
      -25.175175175175177, -25.275275275275277, -25.375375375375377,
      -25.475475475475477, -25.575575575575577, -25.675675675675677,
      -25.775775775775777, -25.875875875875877, -25.975975975975977,
      -26.076076076076077, -26.176176176176178, -26.276276276276278,
      -26.376376376376378, -26.476476476476478, -26.576576576576578,
      -26.676676676676678, -26.776776776776778, -26.876876876876878,
      -26.976976976976978, -27.077077077077078, -27.177177177177178,
      -27.277277277277278, -27.377377377377378, -27.477477477477478,
      -27.577577577577578, -27.677677677677679, -27.777777777777779,
      -27.877877877877879, -27.977977977977979, -28.078078078078079,
      -28.178178178178179, -28.278278278278279, -28.378378378378379,
      -28.478478478478479, -28.578578578578579, -28.678678678678679,
      -28.778778778778779, -28.878878878878879, -28.978978978978979,
      -29.079079079079079, -29.179179179179179, -29.27927927927928,
      -29.37937937937938, -29.47947947947948, -29.57957957957958,
      -29.67967967967968, -29.77977977977978, -29.87987987987988,
      -29.979979979979984, -30.080080080080084, -30.180180180180184,
      -30.280280280280284, -30.380380380380384, -30.480480480480484,
      -30.580580580580584, -30.680680680680684, -30.780780780780784,
      -30.880880880880884, -30.980980980980984, -31.081081081081084,
      -31.181181181181184, -31.281281281281284, -31.381381381381384,
      -31.481481481481485, -31.581581581581585, -31.681681681681685,
      -31.781781781781785, -31.881881881881885, -31.981981981981985,
      -32.082082082082081, -32.182182182182181, -32.282282282282281,
      -32.382382382382382, -32.482482482482482, -32.582582582582582,
      -32.682682682682682, -32.782782782782782, -32.882882882882882,
      -32.982982982982982, -33.083083083083082, -33.183183183183182,
      -33.283283283283282, -33.383383383383382, -33.483483483483482,
      -33.583583583583582, -33.683683683683682, -33.783783783783782,
      -33.883883883883883, -33.983983983983983, -34.084084084084083,
      -34.18418418418419, -34.28428428428429, -34.38438438438439,
      -34.48448448448449, -34.58458458458459, -34.68468468468469,
      -34.78478478478479, -34.88488488488489, -34.98498498498499,
      -35.08508508508509, -35.18518518518519, -35.285285285285291,
      -35.385385385385391, -35.485485485485491, -35.585585585585591,
      -35.685685685685691, -35.785785785785791, -35.885885885885891,
      -35.985985985985991, -36.086086086086091, -36.186186186186191,
      -36.286286286286291, -36.386386386386391, -36.486486486486491,
      -36.586586586586591, -36.686686686686691, -36.786786786786791,
      -36.886886886886892, -36.986986986986992, -37.087087087087092,
      -37.187187187187192, -37.287287287287292, -37.387387387387392,
      -37.487487487487492, -37.587587587587592, -37.687687687687692,
      -37.787787787787792, -37.887887887887892, -37.987987987987992,
      -38.088088088088092, -38.188188188188192, -38.288288288288292,
      -38.388388388388393, -38.488488488488493, -38.588588588588593,
      -38.688688688688693, -38.788788788788793, -38.888888888888893,
      -38.988988988988993, -39.089089089089093, -39.189189189189193,
      -39.289289289289293, -39.389389389389393, -39.489489489489493,
      -39.589589589589593, -39.689689689689693, -39.789789789789793,
      -39.889889889889893, -39.989989989989994, -40.090090090090094,
      -40.190190190190194, -40.290290290290294, -40.390390390390394,
      -40.490490490490494, -40.590590590590594, -40.690690690690694,
      -40.790790790790794, -40.890890890890894, -40.990990990990994,
      -41.091091091091094, -41.191191191191194, -41.291291291291294,
      -41.391391391391394, -41.491491491491495, -41.591591591591595,
      -41.691691691691695, -41.791791791791795, -41.891891891891895,
      -41.991991991991995, -42.092092092092095, -42.192192192192195,
      -42.292292292292295, -42.392392392392395, -42.492492492492495,
      -42.592592592592595, -42.692692692692695, -42.792792792792795,
      -42.892892892892895, -42.992992992992995, -43.093093093093096,
      -43.193193193193196, -43.293293293293296, -43.393393393393396,
      -43.493493493493496, -43.593593593593596, -43.693693693693696,
      -43.793793793793796, -43.893893893893896, -43.993993993993996,
      -44.094094094094096, -44.194194194194196, -44.294294294294296,
      -44.394394394394396, -44.4944944944945, -44.5945945945946,
      -44.6946946946947, -44.7947947947948, -44.8948948948949, -44.994994994995,
      -45.0950950950951, -45.1951951951952, -45.2952952952953, -45.3953953953954,
      -45.4954954954955, -45.5955955955956, -45.6956956956957, -45.7957957957958,
      -45.8958958958959, -45.995995995996, -46.0960960960961, -46.1961961961962,
      -46.2962962962963, -46.3963963963964, -46.4964964964965, -46.5965965965966,
      -46.6966966966967, -46.7967967967968, -46.8968968968969, -46.996996996997,
      -47.0970970970971, -47.1971971971972, -47.2972972972973, -47.3973973973974,
      -47.4974974974975, -47.5975975975976, -47.6976976976977, -47.7977977977978,
      -47.8978978978979, -47.997997997998, -48.0980980980981, -48.1981981981982,
      -48.2982982982983, -48.3983983983984, -48.4984984984985, -48.5985985985986,
      -48.6986986986987, -48.7987987987988, -48.8988988988989, -48.998998998999,
      -49.0990990990991, -49.1991991991992, -49.2992992992993, -49.3993993993994,
      -49.4994994994995, -49.5995995995996, -49.6996996996997, -49.7997997997998,
      -49.8998998998999, -50.0, -50.0, -49.949949949949946, -49.8998998998999,
      -49.849849849849846, -49.7997997997998, -49.749749749749746,
      -49.6996996996997, -49.649649649649646, -49.5995995995996,
      -49.549549549549546, -49.4994994994995, -49.449449449449446,
      -49.3993993993994, -49.349349349349346, -49.2992992992993,
      -49.249249249249246, -49.1991991991992, -49.149149149149146,
      -49.0990990990991, -49.049049049049046, -48.998998998999,
      -48.948948948948946, -48.8988988988989, -48.848848848848846,
      -48.7987987987988, -48.748748748748746, -48.6986986986987,
      -48.648648648648646, -48.5985985985986, -48.548548548548546,
      -48.4984984984985, -48.448448448448445, -48.3983983983984,
      -48.348348348348345, -48.2982982982983, -48.248248248248245,
      -48.1981981981982, -48.148148148148145, -48.0980980980981,
      -48.048048048048045, -47.997997997998, -47.947947947947945,
      -47.8978978978979, -47.847847847847845, -47.7977977977978,
      -47.747747747747745, -47.6976976976977, -47.647647647647645,
      -47.5975975975976, -47.547547547547545, -47.4974974974975,
      -47.447447447447445, -47.3973973973974, -47.347347347347345,
      -47.2972972972973, -47.247247247247245, -47.1971971971972,
      -47.147147147147145, -47.0970970970971, -47.047047047047045,
      -46.996996996997, -46.946946946946944, -46.8968968968969,
      -46.846846846846844, -46.7967967967968, -46.746746746746744,
      -46.6966966966967, -46.646646646646644, -46.5965965965966,
      -46.546546546546544, -46.4964964964965, -46.446446446446444,
      -46.3963963963964, -46.346346346346344, -46.2962962962963,
      -46.246246246246244, -46.1961961961962, -46.146146146146144,
      -46.0960960960961, -46.046046046046044, -45.995995995996,
      -45.945945945945944, -45.8958958958959, -45.845845845845844,
      -45.7957957957958, -45.745745745745744, -45.6956956956957,
      -45.645645645645644, -45.5955955955956, -45.545545545545544,
      -45.4954954954955, -45.445445445445444, -45.3953953953954,
      -45.345345345345343, -45.2952952952953, -45.245245245245243,
      -45.1951951951952, -45.145145145145143, -45.0950950950951,
      -45.045045045045043, -44.994994994995, -44.944944944944943,
      -44.8948948948949, -44.844844844844843, -44.7947947947948,
      -44.744744744744743, -44.6946946946947, -44.644644644644643,
      -44.5945945945946, -44.544544544544543, -44.4944944944945,
      -44.444444444444443, -44.394394394394396, -44.344344344344343,
      -44.294294294294296, -44.244244244244243, -44.194194194194196,
      -44.144144144144143, -44.094094094094096, -44.044044044044043,
      -43.993993993993996, -43.943943943943943, -43.893893893893896,
      -43.843843843843842, -43.793793793793796, -43.743743743743742,
      -43.693693693693696, -43.643643643643642, -43.593593593593596,
      -43.543543543543542, -43.493493493493496, -43.443443443443442,
      -43.393393393393396, -43.343343343343342, -43.293293293293296,
      -43.243243243243242, -43.193193193193196, -43.143143143143142,
      -43.093093093093096, -43.043043043043042, -42.992992992992995,
      -42.942942942942942, -42.892892892892895, -42.842842842842842,
      -42.792792792792795, -42.742742742742742, -42.692692692692695,
      -42.642642642642642, -42.592592592592595, -42.542542542542542,
      -42.492492492492495, -42.442442442442442, -42.392392392392395,
      -42.342342342342342, -42.292292292292295, -42.242242242242241,
      -42.192192192192195, -42.142142142142141, -42.092092092092095,
      -42.042042042042041, -41.991991991991995, -41.941941941941941,
      -41.891891891891888, -41.841841841841841, -41.791791791791795,
      -41.741741741741741, -41.691691691691688, -41.641641641641641,
      -41.591591591591595, -41.541541541541541, -41.491491491491487,
      -41.441441441441441, -41.391391391391394, -41.341341341341341,
      -41.291291291291287, -41.241241241241241, -41.191191191191194,
      -41.141141141141141, -41.091091091091087, -41.041041041041041,
      -40.990990990990994, -40.940940940940941, -40.890890890890887,
      -40.840840840840841, -40.790790790790794, -40.74074074074074,
      -40.690690690690687, -40.64064064064064, -40.590590590590594,
      -40.54054054054054, -40.490490490490487, -40.44044044044044,
      -40.390390390390394, -40.34034034034034, -40.290290290290287,
      -40.24024024024024, -40.190190190190194, -40.14014014014014,
      -40.090090090090087, -40.04004004004004, -39.989989989989994,
      -39.93993993993994, -39.889889889889886, -39.83983983983984,
      -39.789789789789793, -39.73973973973974, -39.689689689689686,
      -39.63963963963964, -39.589589589589593, -39.53953953953954,
      -39.489489489489486, -39.43943943943944, -39.389389389389393,
      -39.33933933933934, -39.289289289289286, -39.23923923923924,
      -39.189189189189193, -39.139139139139139, -39.089089089089086,
      -39.039039039039039, -38.988988988988993, -38.938938938938939,
      -38.888888888888886, -38.838838838838839, -38.788788788788793,
      -38.738738738738739, -38.688688688688686, -38.638638638638639,
      -38.588588588588593, -38.538538538538539, -38.488488488488485,
      -38.438438438438439, -38.388388388388393, -38.338338338338339,
      -38.288288288288285, -38.238238238238239, -38.188188188188192,
      -38.138138138138139, -38.088088088088085, -38.038038038038039,
      -37.987987987987992, -37.937937937937939, -37.887887887887885,
      -37.837837837837839, -37.787787787787792, -37.737737737737739,
      -37.687687687687685, -37.637637637637638, -37.587587587587592,
      -37.537537537537538, -37.487487487487485, -37.437437437437438,
      -37.387387387387385, -37.337337337337338, -37.287287287287285,
      -37.237237237237238, -37.187187187187185, -37.137137137137138,
      -37.087087087087085, -37.037037037037038, -36.986986986986985,
      -36.936936936936938, -36.886886886886884, -36.836836836836838,
      -36.786786786786784, -36.736736736736738, -36.686686686686684,
      -36.636636636636638, -36.586586586586584, -36.536536536536538,
      -36.486486486486484, -36.436436436436438, -36.386386386386384,
      -36.336336336336338, -36.286286286286284, -36.236236236236238,
      -36.186186186186184, -36.136136136136138, -36.086086086086084,
      -36.036036036036037, -35.985985985985984, -35.935935935935937,
      -35.885885885885884, -35.835835835835837, -35.785785785785784,
      -35.735735735735737, -35.685685685685684, -35.635635635635637,
      -35.585585585585584, -35.535535535535537, -35.485485485485484,
      -35.435435435435437, -35.385385385385383, -35.335335335335337,
      -35.285285285285283, -35.235235235235237, -35.185185185185183,
      -35.135135135135137, -35.085085085085083, -35.035035035035037,
      -34.984984984984983, -34.934934934934937, -34.884884884884883,
      -34.834834834834837, -34.784784784784783, -34.734734734734737,
      -34.684684684684683, -34.634634634634637, -34.584584584584583,
      -34.534534534534536, -34.484484484484483, -34.434434434434436,
      -34.384384384384383, -34.334334334334336, -34.284284284284283,
      -34.234234234234236, -34.184184184184183, -34.134134134134136,
      -34.084084084084083, -34.034034034034036, -33.983983983983983,
      -33.933933933933929, -33.883883883883883, -33.833833833833836,
      -33.783783783783782, -33.733733733733729, -33.683683683683682,
      -33.633633633633636, -33.583583583583582, -33.533533533533529,
      -33.483483483483482, -33.433433433433436, -33.383383383383382,
      -33.333333333333329, -33.283283283283282, -33.233233233233236,
      -33.183183183183182, -33.133133133133128, -33.083083083083082,
      -33.033033033033036, -32.982982982982982, -32.932932932932928,
      -32.882882882882882, -32.832832832832835, -32.782782782782782,
      -32.732732732732728, -32.682682682682682, -32.632632632632635,
      -32.582582582582582, -32.532532532532528, -32.482482482482482,
      -32.432432432432435, -32.382382382382382, -32.332332332332328,
      -32.282282282282281, -32.232232232232235, -32.182182182182181,
      -32.132132132132128, -32.082082082082081, -32.032032032032035,
      -31.981981981981981, -31.931931931931931, -31.881881881881881,
      -31.831831831831831, -31.781781781781781, -31.731731731731731,
      -31.681681681681681, -31.631631631631631, -31.581581581581581,
      -31.531531531531531, -31.481481481481481, -31.431431431431431,
      -31.381381381381381, -31.331331331331331, -31.281281281281281,
      -31.231231231231231, -31.181181181181181, -31.131131131131131,
      -31.081081081081081, -31.031031031031031, -30.980980980980981,
      -30.930930930930931, -30.880880880880881, -30.830830830830831,
      -30.780780780780781, -30.73073073073073, -30.68068068068068,
      -30.63063063063063, -30.58058058058058, -30.53053053053053,
      -30.48048048048048, -30.43043043043043, -30.38038038038038,
      -30.33033033033033, -30.28028028028028, -30.23023023023023,
      -30.18018018018018, -30.13013013013013, -30.08008008008008,
      -30.03003003003003, -29.97997997997998, -29.92992992992993,
      -29.87987987987988, -29.82982982982983, -29.77977977977978,
      -29.72972972972973, -29.67967967967968, -29.62962962962963,
      -29.57957957957958, -29.52952952952953, -29.47947947947948,
      -29.42942942942943, -29.37937937937938, -29.32932932932933,
      -29.27927927927928, -29.22922922922923, -29.179179179179179,
      -29.129129129129129, -29.079079079079079, -29.029029029029029,
      -28.978978978978979, -28.928928928928929, -28.878878878878879,
      -28.828828828828829, -28.778778778778779, -28.728728728728729,
      -28.678678678678679, -28.628628628628629, -28.578578578578579,
      -28.528528528528529, -28.478478478478479, -28.428428428428429,
      -28.378378378378379, -28.328328328328329, -28.278278278278279,
      -28.228228228228229, -28.178178178178179, -28.128128128128129,
      -28.078078078078079, -28.028028028028029, -27.977977977977979,
      -27.927927927927929, -27.877877877877879, -27.827827827827829,
      -27.777777777777779, -27.727727727727729, -27.677677677677679,
      -27.627627627627628, -27.577577577577578, -27.527527527527528,
      -27.477477477477478, -27.427427427427428, -27.377377377377378,
      -27.327327327327328, -27.277277277277278, -27.227227227227228,
      -27.177177177177178, -27.127127127127128, -27.077077077077078,
      -27.027027027027028, -26.976976976976978, -26.926926926926928,
      -26.876876876876878, -26.826826826826828, -26.776776776776778,
      -26.726726726726728, -26.676676676676678, -26.626626626626628,
      -26.576576576576578, -26.526526526526528, -26.476476476476478,
      -26.426426426426428, -26.376376376376378, -26.326326326326328,
      -26.276276276276278, -26.226226226226228, -26.176176176176178,
      -26.126126126126128, -26.076076076076077, -26.026026026026027,
      -25.975975975975977, -25.925925925925927, -25.875875875875877,
      -25.825825825825827, -25.775775775775777, -25.725725725725727,
      -25.675675675675677, -25.625625625625627, -25.575575575575577,
      -25.525525525525527, -25.475475475475477, -25.425425425425427,
      -25.375375375375377, -25.325325325325327, -25.275275275275277,
      -25.225225225225227, -25.175175175175177, -25.125125125125127,
      -25.075075075075077, -25.025025025025027, -24.974974974974973,
      -24.924924924924923, -24.874874874874873, -24.824824824824823,
      -24.774774774774773, -24.724724724724723, -24.674674674674673,
      -24.624624624624623, -24.574574574574573, -24.524524524524523,
      -24.474474474474473, -24.424424424424423, -24.374374374374373,
      -24.324324324324323, -24.274274274274273, -24.224224224224223,
      -24.174174174174173, -24.124124124124123, -24.074074074074073,
      -24.024024024024023, -23.973973973973973, -23.923923923923923,
      -23.873873873873872, -23.823823823823822, -23.773773773773772,
      -23.723723723723722, -23.673673673673672, -23.623623623623622,
      -23.573573573573572, -23.523523523523522, -23.473473473473472,
      -23.423423423423422, -23.373373373373372, -23.323323323323322,
      -23.273273273273272, -23.223223223223222, -23.173173173173172,
      -23.123123123123122, -23.073073073073072, -23.023023023023022,
      -22.972972972972972, -22.922922922922922, -22.872872872872872,
      -22.822822822822822, -22.772772772772772, -22.722722722722722,
      -22.672672672672672, -22.622622622622622, -22.572572572572572,
      -22.522522522522522, -22.472472472472472, -22.422422422422422,
      -22.372372372372372, -22.322322322322321, -22.272272272272271,
      -22.222222222222221, -22.172172172172171, -22.122122122122121,
      -22.072072072072071, -22.022022022022021, -21.971971971971971,
      -21.921921921921921, -21.871871871871871, -21.821821821821821,
      -21.771771771771771, -21.721721721721721, -21.671671671671671,
      -21.621621621621621, -21.571571571571571, -21.521521521521521,
      -21.471471471471471, -21.421421421421421, -21.371371371371371,
      -21.321321321321321, -21.271271271271271, -21.221221221221221,
      -21.171171171171171, -21.121121121121121, -21.071071071071071,
      -21.021021021021021, -20.970970970970971, -20.920920920920921,
      -20.870870870870871, -20.820820820820821, -20.77077077077077,
      -20.72072072072072, -20.67067067067067, -20.62062062062062,
      -20.57057057057057, -20.52052052052052, -20.47047047047047,
      -20.42042042042042, -20.37037037037037, -20.32032032032032,
      -20.27027027027027, -20.22022022022022, -20.17017017017017,
      -20.12012012012012, -20.07007007007007, -20.02002002002002,
      -19.96996996996997, -19.91991991991992, -19.86986986986987,
      -19.81981981981982, -19.76976976976977, -19.71971971971972,
      -19.66966966966967, -19.61961961961962, -19.56956956956957,
      -19.51951951951952, -19.46946946946947, -19.41941941941942,
      -19.36936936936937, -19.31931931931932, -19.26926926926927,
      -19.219219219219219, -19.169169169169169, -19.119119119119119,
      -19.069069069069069, -19.019019019019019, -18.968968968968969,
      -18.918918918918919, -18.868868868868869, -18.818818818818819,
      -18.768768768768769, -18.718718718718719, -18.668668668668669,
      -18.618618618618619, -18.568568568568569, -18.518518518518519,
      -18.468468468468469, -18.418418418418419, -18.368368368368369,
      -18.318318318318319, -18.268268268268269, -18.218218218218219,
      -18.168168168168169, -18.118118118118119, -18.068068068068069,
      -18.018018018018019, -17.967967967967965, -17.917917917917919,
      -17.867867867867865, -17.817817817817819, -17.767767767767765,
      -17.717717717717719, -17.667667667667665, -17.617617617617618,
      -17.567567567567565, -17.517517517517518, -17.467467467467465,
      -17.417417417417418, -17.367367367367365, -17.317317317317318,
      -17.267267267267265, -17.217217217217218, -17.167167167167165,
      -17.117117117117118, -17.067067067067065, -17.017017017017018,
      -16.966966966966964, -16.916916916916918, -16.866866866866864,
      -16.816816816816818, -16.766766766766764, -16.716716716716718,
      -16.666666666666664, -16.616616616616618, -16.566566566566564,
      -16.516516516516518, -16.466466466466464, -16.416416416416418,
      -16.366366366366364, -16.316316316316318, -16.266266266266264,
      -16.216216216216218, -16.166166166166164, -16.116116116116117,
      -16.066066066066064, -16.016016016016017, -15.965965965965964,
      -15.915915915915917, -15.865865865865864, -15.815815815815817,
      -15.765765765765764, -15.715715715715717, -15.665665665665664,
      -15.615615615615617, -15.565565565565564, -15.515515515515517,
      -15.465465465465464, -15.415415415415417, -15.365365365365363,
      -15.315315315315317, -15.265265265265263, -15.215215215215217,
      -15.165165165165163, -15.115115115115117, -15.065065065065063,
      -15.015015015015017, -14.964964964964963, -14.914914914914917,
      -14.864864864864863, -14.814814814814817, -14.764764764764763,
      -14.714714714714717, -14.664664664664663, -14.614614614614617,
      -14.564564564564563, -14.514514514514516, -14.464464464464463,
      -14.414414414414416, -14.364364364364363, -14.314314314314316,
      -14.264264264264263, -14.214214214214216, -14.164164164164163,
      -14.114114114114116, -14.064064064064063, -14.014014014014016,
      -13.963963963963963, -13.913913913913916, -13.863863863863862,
      -13.813813813813816, -13.763763763763762, -13.713713713713716,
      -13.663663663663662, -13.613613613613616, -13.563563563563562,
      -13.513513513513516, -13.463463463463462, -13.413413413413416,
      -13.363363363363362, -13.313313313313316, -13.263263263263262,
      -13.213213213213216, -13.163163163163162, -13.113113113113116,
      -13.063063063063062, -13.013013013013015, -12.962962962962962,
      -12.912912912912915, -12.862862862862862, -12.812812812812815,
      -12.762762762762762, -12.712712712712715, -12.662662662662662,
      -12.612612612612615, -12.562562562562562, -12.512512512512515,
      -12.462462462462462, -12.412412412412415, -12.362362362362362,
      -12.312312312312315, -12.262262262262261, -12.212212212212215,
      -12.162162162162161, -12.112112112112115, -12.062062062062061,
      -12.012012012012015, -11.961961961961961, -11.911911911911915,
      -11.861861861861861, -11.811811811811815, -11.761761761761761,
      -11.711711711711715, -11.661661661661661, -11.611611611611615,
      -11.561561561561561, -11.511511511511515, -11.461461461461461,
      -11.411411411411414, -11.361361361361361, -11.311311311311314,
      -11.261261261261261, -11.211211211211214, -11.161161161161161,
      -11.111111111111114, -11.061061061061061, -11.011011011011014,
      -10.960960960960961, -10.910910910910914, -10.860860860860861,
      -10.810810810810814, -10.76076076076076, -10.710710710710714,
      -10.66066066066066, -10.610610610610614, -10.56056056056056,
      -10.510510510510514, -10.46046046046046, -10.410410410410414,
      -10.36036036036036, -10.310310310310314, -10.26026026026026,
      -10.210210210210214, -10.16016016016016, -10.110110110110114,
      -10.06006006006006, -10.010010010010014, -9.95995995995996,
      -9.9099099099099064, -9.85985985985986, -9.8098098098098063,
      -9.75975975975976, -9.7097097097097063, -9.65965965965966,
      -9.6096096096096062, -9.55955955955956, -9.5095095095095061,
      -9.45945945945946, -9.4094094094094061, -9.35935935935936,
      -9.309309309309306, -9.25925925925926, -9.2092092092092059,
      -9.15915915915916, -9.1091091091091059, -9.05905905905906,
      -9.0090090090090058, -8.95895895895896, -8.9089089089089057,
      -8.85885885885886, -8.8088088088088057, -8.75875875875876,
      -8.7087087087087056, -8.65865865865866, -8.6086086086086056,
      -8.5585585585585591, -8.5085085085085055, -8.458458458458459,
      -8.4084084084084054, -8.3583583583583589, -8.3083083083083054,
      -8.2582582582582589, -8.2082082082082053, -8.1581581581581588,
      -8.1081081081081052, -8.0580580580580587, -8.0080080080080052,
      -7.9579579579579587, -7.9079079079079051, -7.8578578578578586,
      -7.807807807807805, -7.7577577577577586, -7.707707707707705,
      -7.6576576576576585, -7.6076076076076049, -7.5575575575575584,
      -7.5075075075075048, -7.4574574574574584, -7.4074074074074048,
      -7.3573573573573583, -7.3073073073073047, -7.2572572572572582,
      -7.2072072072072046, -7.1571571571571582, -7.1071071071071046,
      -7.0570570570570581, -7.0070070070070045, -6.956956956956958,
      -6.9069069069069045, -6.856856856856858, -6.8068068068068044,
      -6.7567567567567579, -6.7067067067067043, -6.6566566566566578,
      -6.6066066066066043, -6.5565565565565578, -6.5065065065065042,
      -6.4564564564564577, -6.4064064064064041, -6.3563563563563577,
      -6.3063063063063041, -6.2562562562562576, -6.206206206206204,
      -6.1561561561561575, -6.1061061061061039, -6.0560560560560575,
      -6.0060060060060039, -5.9559559559559574, -5.9059059059059038,
      -5.8558558558558573, -5.8058058058058037, -5.7557557557557573,
      -5.7057057057057037, -5.6556556556556572, -5.6056056056056036,
      -5.5555555555555571, -5.5055055055055035, -5.4554554554554571,
      -5.4054054054054035, -5.355355355355357, -5.3053053053053034,
      -5.2552552552552569, -5.2052052052052034, -5.1551551551551569,
      -5.1051051051051033, -5.0550550550550568, -5.0050050050050032,
      -4.9549549549549567, -4.9049049049049032, -4.8548548548548567,
      -4.8048048048048031, -4.7547547547547566, -4.704704704704703,
      -4.6546546546546566, -4.604604604604603, -4.5545545545545565,
      -4.5045045045045029, -4.4544544544544564, -4.4044044044044028,
      -4.3543543543543564, -4.3043043043043028, -4.2542542542542563,
      -4.2042042042042027, -4.1541541541541562, -4.1041041041041026,
      -4.0540540540540562, -4.0040040040040026, -3.9539539539539561,
      -3.9039039039039025, -3.853853853853856, -3.8038038038038025,
      -3.753753753753756, -3.7037037037037024, -3.6536536536536559,
      -3.6036036036036023, -3.5535535535535558, -3.5035035035035023,
      -3.4534534534534558, -3.4034034034034022, -3.3533533533533557,
      -3.3033033033033021, -3.2532532532532557, -3.2032032032032021,
      -3.1531531531531556, -3.103103103103102, -3.0530530530530555,
      -3.0030030030030019, -2.9529529529529555, -2.9029029029029019,
      -2.8528528528528554, -2.8028028028028018, -2.7527527527527553,
      -2.7027027027027017, -2.6526526526526553, -2.6026026026026017,
      -2.5525525525525552, -2.5025025025025016, -2.4524524524524551,
      -2.4024024024024015, -2.3523523523523551, -2.3023023023023015,
      -2.252252252252255, -2.2022022022022014, -2.1521521521521549,
      -2.1021021021021014, -2.0520520520520549, -2.0020020020020013,
      -1.9519519519519548, -1.9019019019019012, -1.8518518518518547,
      -1.8018018018018012, -1.7517517517517547, -1.7017017017017011,
      -1.6516516516516546, -1.601601601601601, -1.5515515515515546,
      -1.501501501501501, -1.4514514514514545, -1.4014014014014009,
      -1.3513513513513544, -1.3013013013013008, -1.2512512512512544,
      -1.2012012012012008, -1.1511511511511543, -1.1011011011011007,
      -1.0510510510510542, -1.0010010010010006, -0.95095095095095417,
      -0.90090090090090058, -0.8508508508508541, -0.80080080080080052,
      -0.750750750750754, -0.70070070070070045, -0.650650650650654,
      -0.60060060060060039, -0.55055055055055391, -0.50050050050050032,
      -0.45045045045045384, -0.40040040040040026, -0.35035035035035378,
      -0.30030030030030019, -0.25025025025025371, -0.20020020020020013,
      -0.15015015015015365, -0.10010010010010006, -0.050050050050053585, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    } ;

    ModelWithControllersOnly_DW.FromWorkspace1_PWORK.TimePtr = static_cast<void *>
      (pTimeValues0);
    ModelWithControllersOnly_DW.FromWorkspace1_PWORK.DataPtr = static_cast<void *>
      (pDataValues0);
    ModelWithControllersOnly_DW.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* Start for FromWorkspace: '<S5>/From Workspace4' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06,
      0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18,
      0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31,
      0.32, 0.33, 0.34, 0.35000000000000003, 0.36, 0.37, 0.38, 0.39, 0.4,
      0.41000000000000003, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47000000000000003,
      0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57000000000000006,
      0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68,
      0.69000000000000006, 0.70000000000000007, 0.71, 0.72, 0.73, 0.74, 0.75,
      0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82000000000000006,
      0.83000000000000007, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92,
      0.93, 0.94000000000000006, 0.95000000000000007, 0.96, 0.97, 0.98, 0.99,
      1.0, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.12,
      1.1300000000000001, 1.1400000000000001, 1.1500000000000001, 1.16, 1.17,
      1.18, 1.19, 1.2, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3,
      1.31, 1.32, 1.33, 1.34, 1.35, 1.36, 1.37, 1.3800000000000001,
      1.3900000000000001, 1.4000000000000001, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46,
      1.47, 1.48, 1.49, 1.5, 1.51, 1.52, 1.53, 1.54, 1.55, 1.56, 1.57, 1.58,
      1.59, 1.6, 1.61, 1.62, 1.6300000000000001, 1.6400000000000001,
      1.6500000000000001, 1.6600000000000001, 1.67, 1.68, 1.69, 1.7, 1.71, 1.72,
      1.73, 1.74, 1.75, 1.76, 1.77, 1.78, 1.79, 1.8, 1.81, 1.82, 1.83, 1.84,
      1.85, 1.86, 1.87, 1.8800000000000001, 1.8900000000000001,
      1.9000000000000001, 1.9100000000000001, 1.92, 1.93, 1.94, 1.95, 1.96, 1.97,
      1.98, 1.99, 2.0, 2.0100000000000002, 2.02, 2.0300000000000002, 2.04, 2.05,
      2.06, 2.07, 2.08, 2.09, 2.1, 2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17,
      2.18, 2.19, 2.2, 2.21, 2.22, 2.23, 2.24, 2.25, 2.2600000000000002, 2.27,
      2.2800000000000002, 2.29, 2.3000000000000003, 2.31, 2.32, 2.33, 2.34, 2.35,
      2.36, 2.37, 2.38, 2.39, 2.4, 2.41, 2.42, 2.43, 2.44, 2.45, 2.46, 2.47,
      2.48, 2.49, 2.5, 2.5100000000000002, 2.52, 2.5300000000000002, 2.54,
      2.5500000000000003, 2.56, 2.57, 2.58, 2.59, 2.6, 2.61, 2.62, 2.63, 2.64,
      2.65, 2.66, 2.67, 2.68, 2.69, 2.7, 2.71, 2.72, 2.73, 2.74, 2.75,
      2.7600000000000002, 2.77, 2.7800000000000002, 2.79, 2.8000000000000003,
      2.81, 2.82, 2.83, 2.84, 2.85, 2.86, 2.87, 2.88, 2.89, 2.9, 2.91, 2.92,
      2.93, 2.94, 2.95, 2.96, 2.97, 2.98, 2.99, 3.0, 3.0100000000000002, 3.02,
      3.0300000000000002, 3.04, 3.0500000000000003, 3.06, 3.0700000000000003,
      3.08, 3.09, 3.1, 3.11, 3.12, 3.13, 3.14, 3.15, 3.16, 3.17, 3.18, 3.19, 3.2,
      3.21, 3.22, 3.23, 3.24, 3.25, 3.2600000000000002, 3.27, 3.2800000000000002,
      3.29, 3.3000000000000003, 3.31, 3.3200000000000003, 3.33, 3.34, 3.35, 3.36,
      3.37, 3.38, 3.39, 3.4, 3.41, 3.42, 3.43, 3.44, 3.45, 3.46, 3.47, 3.48,
      3.49, 3.5, 3.5100000000000002, 3.52, 3.5300000000000002, 3.54,
      3.5500000000000003, 3.56, 3.5700000000000003, 3.58, 3.59, 3.6, 3.61, 3.62,
      3.63, 3.64, 3.65, 3.66, 3.67, 3.68, 3.69, 3.7, 3.71, 3.72, 3.73, 3.74,
      3.75, 3.7600000000000002, 3.77, 3.7800000000000002, 3.79,
      3.8000000000000003, 3.81, 3.8200000000000003, 3.83, 3.84, 3.85, 3.86, 3.87,
      3.88, 3.89, 3.9, 3.91, 3.92, 3.93, 3.94, 3.95, 3.96, 3.97, 3.98, 3.99, 4.0,
      4.01, 4.0200000000000005, 4.03, 4.04, 4.05, 4.0600000000000005, 4.07, 4.08,
      4.09, 4.1, 4.11, 4.12, 4.13, 4.14, 4.15, 4.16, 4.17, 4.18, 4.19, 4.2, 4.21,
      4.22, 4.23, 4.24, 4.25, 4.26, 4.2700000000000005, 4.28, 4.29, 4.3,
      4.3100000000000005, 4.32, 4.33, 4.34, 4.3500000000000005, 4.36, 4.37, 4.38,
      4.39, 4.4, 4.41, 4.42, 4.43, 4.44, 4.45, 4.46, 4.47, 4.48, 4.49, 4.5, 4.51,
      4.5200000000000005, 4.53, 4.54, 4.55, 4.5600000000000005, 4.57, 4.58, 4.59,
      4.6000000000000005, 4.61, 4.62, 4.63, 4.64, 4.65, 4.66, 4.67, 4.68, 4.69,
      4.7, 4.71, 4.72, 4.73, 4.74, 4.75, 4.76, 4.7700000000000005, 4.78, 4.79,
      4.8, 4.8100000000000005, 4.82, 4.83, 4.84, 4.8500000000000005, 4.86, 4.87,
      4.88, 4.89, 4.9, 4.91, 4.92, 4.93, 4.94, 4.95, 4.96, 4.97, 4.98, 4.99, 5.0,
      5.01, 5.0200000000000005, 5.03, 5.04, 5.05, 5.0600000000000005, 5.07, 5.08,
      5.09, 5.1000000000000005, 5.11, 5.12, 5.13, 5.14, 5.15, 5.16, 5.17, 5.18,
      5.19, 5.2, 5.21, 5.22, 5.23, 5.24, 5.25, 5.26, 5.2700000000000005, 5.28,
      5.29, 5.3, 5.3100000000000005, 5.32, 5.33, 5.34, 5.3500000000000005, 5.36,
      5.37, 5.38, 5.39, 5.4, 5.41, 5.42, 5.43, 5.44, 5.45, 5.46, 5.47, 5.48,
      5.49, 5.5, 5.51, 5.5200000000000005, 5.53, 5.54, 5.55, 5.5600000000000005,
      5.57, 5.58, 5.59, 5.6000000000000005, 5.61, 5.62, 5.63, 5.64, 5.65, 5.66,
      5.67, 5.68, 5.69, 5.7, 5.71, 5.72, 5.73, 5.74, 5.75, 5.76,
      5.7700000000000005, 5.78, 5.79, 5.8, 5.8100000000000005, 5.82, 5.83, 5.84,
      5.8500000000000005, 5.86, 5.87, 5.88, 5.89, 5.9, 5.91, 5.92, 5.93, 5.94,
      5.95, 5.96, 5.97, 5.98, 5.99, 6.0, 6.01, 6.0200000000000005, 6.03, 6.04,
      6.05, 6.0600000000000005, 6.07, 6.08, 6.09, 6.1000000000000005, 6.11, 6.12,
      6.13, 6.1400000000000006, 6.15, 6.16, 6.17, 6.18, 6.19, 6.2, 6.21, 6.22,
      6.23, 6.24, 6.25, 6.26, 6.2700000000000005, 6.28, 6.29, 6.3,
      6.3100000000000005, 6.32, 6.33, 6.34, 6.3500000000000005, 6.36, 6.37, 6.38,
      6.3900000000000006, 6.4, 6.41, 6.42, 6.43, 6.44, 6.45, 6.46, 6.47, 6.48,
      6.49, 6.5, 6.51, 6.5200000000000005, 6.53, 6.54, 6.55, 6.5600000000000005,
      6.57, 6.58, 6.59, 6.6000000000000005, 6.61, 6.62, 6.63, 6.6400000000000006,
      6.65, 6.66, 6.67, 6.68, 6.69, 6.7, 6.71, 6.72, 6.73, 6.74, 6.75, 6.76,
      6.7700000000000005, 6.78, 6.79, 6.8, 6.8100000000000005, 6.82, 6.83, 6.84,
      6.8500000000000005, 6.86, 6.87, 6.88, 6.8900000000000006, 6.9, 6.91, 6.92,
      6.93, 6.94, 6.95, 6.96, 6.97, 6.98, 6.99, 7.0, 7.01, 7.0200000000000005,
      7.03, 7.04, 7.05, 7.0600000000000005, 7.07, 7.08, 7.09, 7.1000000000000005,
      7.11, 7.12, 7.13, 7.1400000000000006, 7.15, 7.16, 7.17, 7.18, 7.19, 7.2,
      7.21, 7.22, 7.23, 7.24, 7.25, 7.26, 7.2700000000000005, 7.28, 7.29, 7.3,
      7.3100000000000005, 7.32, 7.33, 7.34, 7.3500000000000005, 7.36, 7.37, 7.38,
      7.3900000000000006, 7.4, 7.41, 7.42, 7.43, 7.44, 7.45, 7.46, 7.47, 7.48,
      7.49, 7.5, 7.51, 7.5200000000000005, 7.53, 7.54, 7.55, 7.5600000000000005,
      7.57, 7.58, 7.59, 7.6000000000000005, 7.61, 7.62, 7.63, 7.6400000000000006,
      7.65, 7.66, 7.67, 7.68, 7.69, 7.7, 7.71, 7.72, 7.73, 7.74, 7.75, 7.76,
      7.7700000000000005, 7.78, 7.79, 7.8, 7.8100000000000005, 7.82, 7.83, 7.84,
      7.8500000000000005, 7.86, 7.87, 7.88, 7.8900000000000006, 7.9, 7.91, 7.92,
      7.9300000000000006, 7.94, 7.95, 7.96, 7.97, 7.98, 7.99, 8.0, 8.01, 8.02,
      8.03, 8.0400000000000009, 8.05, 8.06, 8.07, 8.08, 8.09, 8.1, 8.11,
      8.120000000000001, 8.13, 8.14, 8.15, 8.16, 8.17, 8.18, 8.19, 8.2, 8.21,
      8.22, 8.23, 8.24, 8.25, 8.26, 8.27, 8.28, 8.2900000000000009, 8.3, 8.31,
      8.32, 8.33, 8.34, 8.35, 8.36, 8.370000000000001, 8.38, 8.39, 8.4, 8.41,
      8.42, 8.43, 8.44, 8.45, 8.46, 8.47, 8.48, 8.49, 8.5, 8.51, 8.52, 8.53,
      8.5400000000000009, 8.55, 8.56, 8.57, 8.58, 8.59, 8.6, 8.61,
      8.620000000000001, 8.63, 8.64, 8.65, 8.66, 8.67, 8.68, 8.69,
      8.7000000000000011, 8.71, 8.72, 8.73, 8.74, 8.75, 8.76, 8.77, 8.78,
      8.7900000000000009, 8.8, 8.81, 8.82, 8.83, 8.84, 8.85, 8.86,
      8.870000000000001, 8.88, 8.89, 8.9, 8.91, 8.92, 8.93, 8.94,
      8.9500000000000011, 8.96, 8.97, 8.98, 8.99, 9.0, 9.01, 9.02, 9.03,
      9.0400000000000009, 9.05, 9.06, 9.07, 9.08, 9.09, 9.1, 9.11,
      9.120000000000001, 9.13, 9.14, 9.15, 9.16, 9.17, 9.18, 9.19,
      9.2000000000000011, 9.21, 9.22, 9.23, 9.24, 9.25, 9.26, 9.27, 9.28,
      9.2900000000000009, 9.3, 9.31, 9.32, 9.33, 9.34, 9.35, 9.36,
      9.370000000000001, 9.38, 9.39, 9.4, 9.41, 9.42, 9.43, 9.44,
      9.4500000000000011, 9.46, 9.47, 9.48, 9.49, 9.5, 9.51, 9.52, 9.53,
      9.5400000000000009, 9.55, 9.56, 9.57, 9.58, 9.59, 9.6, 9.61,
      9.620000000000001, 9.63, 9.64, 9.65, 9.66, 9.67, 9.68, 9.69,
      9.7000000000000011, 9.71, 9.72, 9.73, 9.74, 9.75, 9.76, 9.77, 9.78,
      9.7900000000000009, 9.8, 9.81, 9.82, 9.83, 9.84, 9.85, 9.86,
      9.870000000000001, 9.88, 9.89, 9.9, 9.91, 9.92, 9.93, 9.94,
      9.9500000000000011, 9.96, 9.97, 9.98, 9.99, 10.0, 10.01, 10.02, 10.03,
      10.040000000000001, 10.05, 10.06, 10.07, 10.08, 10.09, 10.1, 10.11,
      10.120000000000001, 10.13, 10.14, 10.15, 10.16, 10.17, 10.18, 10.19,
      10.200000000000001, 10.21, 10.22, 10.23, 10.24, 10.25, 10.26, 10.27, 10.28,
      10.290000000000001, 10.3, 10.31, 10.32, 10.33, 10.34, 10.35, 10.36,
      10.370000000000001, 10.38, 10.39, 10.4, 10.41, 10.42, 10.43, 10.44,
      10.450000000000001, 10.46, 10.47, 10.48, 10.49, 10.5, 10.51, 10.52, 10.53,
      10.540000000000001, 10.55, 10.56, 10.57, 10.58, 10.59, 10.6, 10.61,
      10.620000000000001, 10.63, 10.64, 10.65, 10.66, 10.67, 10.68, 10.69,
      10.700000000000001, 10.71, 10.72, 10.73, 10.74, 10.75, 10.76, 10.77, 10.78,
      10.790000000000001, 10.8, 10.81, 10.82, 10.83, 10.84, 10.85, 10.86,
      10.870000000000001, 10.88, 10.89, 10.9, 10.91, 10.92, 10.93, 10.94,
      10.950000000000001, 10.96, 10.97, 10.98, 10.99, 11.0, 11.01, 11.02, 11.03,
      11.040000000000001, 11.05, 11.06, 11.07, 11.08, 11.09, 11.1, 11.11,
      11.120000000000001, 11.13, 11.14, 11.15, 11.16, 11.17, 11.18, 11.19,
      11.200000000000001, 11.21, 11.22, 11.23, 11.24, 11.25, 11.26, 11.27, 11.28,
      11.290000000000001, 11.3, 11.31, 11.32, 11.33, 11.34, 11.35, 11.36,
      11.370000000000001, 11.38, 11.39, 11.4, 11.41, 11.42, 11.43, 11.44,
      11.450000000000001, 11.46, 11.47, 11.48, 11.49, 11.5, 11.51, 11.52, 11.53,
      11.540000000000001, 11.55, 11.56, 11.57, 11.58, 11.59, 11.6, 11.61,
      11.620000000000001, 11.63, 11.64, 11.65, 11.66, 11.67, 11.68, 11.69,
      11.700000000000001, 11.71, 11.72, 11.73, 11.74, 11.75, 11.76, 11.77, 11.78,
      11.790000000000001, 11.8, 11.81, 11.82, 11.83, 11.84, 11.85, 11.86,
      11.870000000000001, 11.88, 11.89, 11.9, 11.91, 11.92, 11.93, 11.94,
      11.950000000000001, 11.96, 11.97, 11.98, 11.99, 12.0, 12.01, 12.02,
      12.030000000000001, 12.040000000000001, 12.05, 12.06, 12.07, 12.08, 12.09,
      12.1, 12.11, 12.120000000000001, 12.13, 12.14, 12.15, 12.16, 12.17, 12.18,
      12.19, 12.200000000000001, 12.21, 12.22, 12.23, 12.24, 12.25, 12.26, 12.27,
      12.280000000000001, 12.290000000000001, 12.3, 12.31, 12.32, 12.33, 12.34,
      12.35, 12.36, 12.370000000000001, 12.38, 12.39, 12.4, 12.41, 12.42, 12.43,
      12.44, 12.450000000000001, 12.46, 12.47, 12.48, 12.49, 12.5, 12.51, 12.52,
      12.530000000000001, 12.540000000000001, 12.55, 12.56, 12.57, 12.58, 12.59,
      12.6, 12.61, 12.620000000000001, 12.63, 12.64, 12.65, 12.66, 12.67, 12.68,
      12.69, 12.700000000000001, 12.71, 12.72, 12.73, 12.74, 12.75, 12.76, 12.77,
      12.780000000000001, 12.790000000000001, 12.8, 12.81, 12.82, 12.83, 12.84,
      12.85, 12.86, 12.870000000000001, 12.88, 12.89, 12.9, 12.91, 12.92, 12.93,
      12.94, 12.950000000000001, 12.96, 12.97, 12.98, 12.99, 13.0, 13.01, 13.02,
      13.030000000000001, 13.040000000000001, 13.05, 13.06, 13.07, 13.08, 13.09,
      13.1, 13.11, 13.120000000000001, 13.13, 13.14, 13.15, 13.16, 13.17, 13.18,
      13.19, 13.200000000000001, 13.21, 13.22, 13.23, 13.24, 13.25, 13.26, 13.27,
      13.280000000000001, 13.290000000000001, 13.3, 13.31, 13.32, 13.33, 13.34,
      13.35, 13.36, 13.370000000000001, 13.38, 13.39, 13.4, 13.41, 13.42, 13.43,
      13.44, 13.450000000000001, 13.46, 13.47, 13.48, 13.49, 13.5, 13.51, 13.52,
      13.530000000000001, 13.540000000000001, 13.55, 13.56, 13.57, 13.58, 13.59,
      13.6, 13.61, 13.620000000000001, 13.63, 13.64, 13.65, 13.66, 13.67, 13.68,
      13.69, 13.700000000000001, 13.71, 13.72, 13.73, 13.74, 13.75, 13.76, 13.77,
      13.780000000000001, 13.790000000000001, 13.8, 13.81, 13.82, 13.83, 13.84,
      13.85, 13.86, 13.870000000000001, 13.88, 13.89, 13.9, 13.91, 13.92, 13.93,
      13.94, 13.950000000000001, 13.96, 13.97, 13.98, 13.99, 14.0, 14.01, 14.02,
      14.030000000000001, 14.040000000000001, 14.05, 14.06, 14.07, 14.08, 14.09,
      14.1, 14.11, 14.120000000000001, 14.13, 14.14, 14.15, 14.16, 14.17, 14.18,
      14.19, 14.200000000000001, 14.21, 14.22, 14.23, 14.24, 14.25, 14.26, 14.27,
      14.280000000000001, 14.290000000000001, 14.3, 14.31, 14.32, 14.33, 14.34,
      14.35, 14.36, 14.370000000000001, 14.38, 14.39, 14.4, 14.41, 14.42, 14.43,
      14.44, 14.450000000000001, 14.46, 14.47, 14.48, 14.49, 14.5, 14.51, 14.52,
      14.530000000000001, 14.540000000000001, 14.55, 14.56, 14.57, 14.58, 14.59,
      14.6, 14.61, 14.620000000000001, 14.63, 14.64, 14.65, 14.66, 14.67, 14.68,
      14.69, 14.700000000000001, 14.71, 14.72, 14.73, 14.74, 14.75, 14.76, 14.77,
      14.780000000000001, 14.790000000000001, 14.8, 14.81, 14.82, 14.83, 14.84,
      14.85, 14.86, 14.870000000000001, 14.88, 14.89, 14.9, 14.91, 14.92, 14.93,
      14.94, 14.950000000000001, 14.96, 14.97, 14.98, 14.99, 15.0, 15.01, 15.02,
      15.030000000000001, 15.040000000000001, 15.05, 15.06, 15.07, 15.08, 15.09,
      15.1, 15.11, 15.120000000000001, 15.13, 15.14, 15.15, 15.16, 15.17, 15.18,
      15.19, 15.200000000000001, 15.21, 15.22, 15.23, 15.24, 15.25, 15.26, 15.27,
      15.280000000000001, 15.290000000000001, 15.3, 15.31, 15.32, 15.33, 15.34,
      15.35, 15.36, 15.370000000000001, 15.38, 15.39, 15.4, 15.41, 15.42, 15.43,
      15.44, 15.450000000000001, 15.46, 15.47, 15.48, 15.49, 15.5, 15.51, 15.52,
      15.530000000000001, 15.540000000000001, 15.55, 15.56, 15.57, 15.58, 15.59,
      15.6, 15.610000000000001, 15.620000000000001, 15.63, 15.64, 15.65, 15.66,
      15.67, 15.68, 15.69, 15.700000000000001, 15.71, 15.72, 15.73, 15.74, 15.75,
      15.76, 15.77, 15.780000000000001, 15.790000000000001, 15.8, 15.81, 15.82,
      15.83, 15.84, 15.85, 15.860000000000001, 15.870000000000001, 15.88, 15.89,
      15.9, 15.91, 15.92, 15.93, 15.94, 15.950000000000001, 15.96, 15.97, 15.98,
      15.99, 16.0, 16.01, 16.02, 16.03, 16.04, 16.05, 16.06, 16.07,
      16.080000000000002, 16.09, 16.1, 16.11, 16.12, 16.13, 16.14, 16.15, 16.16,
      16.17, 16.18, 16.19, 16.2, 16.21, 16.22, 16.23, 16.240000000000002, 16.25,
      16.26, 16.27, 16.28, 16.29, 16.3, 16.31, 16.32, 16.330000000000002, 16.34,
      16.35, 16.36, 16.37, 16.38, 16.39, 16.4, 16.41, 16.42, 16.43, 16.44, 16.45,
      16.46, 16.47, 16.48, 16.490000000000002, 16.5, 16.51, 16.52, 16.53, 16.54,
      16.55, 16.56, 16.57, 16.580000000000002, 16.59, 16.6, 16.61, 16.62, 16.63,
      16.64, 16.65, 16.66, 16.67, 16.68, 16.69, 16.7, 16.71, 16.72, 16.73,
      16.740000000000002, 16.75, 16.76, 16.77, 16.78, 16.79, 16.8, 16.81, 16.82,
      16.830000000000002, 16.84, 16.85, 16.86, 16.87, 16.88, 16.89, 16.9, 16.91,
      16.92, 16.93, 16.94, 16.95, 16.96, 16.97, 16.98, 16.990000000000002, 17.0,
      17.01, 17.02, 17.03, 17.04, 17.05, 17.06, 17.07, 17.080000000000002, 17.09,
      17.1, 17.11, 17.12, 17.13, 17.14, 17.150000000000002, 17.16, 17.17, 17.18,
      17.19, 17.2, 17.21, 17.22, 17.23, 17.240000000000002, 17.25, 17.26, 17.27,
      17.28, 17.29, 17.3, 17.31, 17.32, 17.330000000000002, 17.34, 17.35, 17.36,
      17.37, 17.38, 17.39, 17.400000000000002, 17.41, 17.42, 17.43, 17.44, 17.45,
      17.46, 17.47, 17.48, 17.490000000000002, 17.5, 17.51, 17.52, 17.53, 17.54,
      17.55, 17.56, 17.57, 17.580000000000002, 17.59, 17.6, 17.61, 17.62, 17.63,
      17.64, 17.650000000000002, 17.66, 17.67, 17.68, 17.69, 17.7, 17.71, 17.72,
      17.73, 17.740000000000002, 17.75, 17.76, 17.77, 17.78, 17.79, 17.8, 17.81,
      17.82, 17.830000000000002, 17.84, 17.85, 17.86, 17.87, 17.88, 17.89,
      17.900000000000002, 17.91, 17.92, 17.93, 17.94, 17.95, 17.96, 17.97, 17.98,
      17.990000000000002, 18.0, 18.01, 18.02, 18.03, 18.04, 18.05, 18.06, 18.07,
      18.080000000000002, 18.09, 18.1, 18.11, 18.12, 18.13, 18.14,
      18.150000000000002, 18.16, 18.17, 18.18, 18.19, 18.2, 18.21, 18.22, 18.23,
      18.240000000000002, 18.25, 18.26, 18.27, 18.28, 18.29, 18.3, 18.31, 18.32,
      18.330000000000002, 18.34, 18.35, 18.36, 18.37, 18.38, 18.39,
      18.400000000000002, 18.41, 18.42, 18.43, 18.44, 18.45, 18.46, 18.47, 18.48,
      18.490000000000002, 18.5, 18.51, 18.52, 18.53, 18.54, 18.55, 18.56, 18.57,
      18.580000000000002, 18.59, 18.6, 18.61, 18.62, 18.63, 18.64,
      18.650000000000002, 18.66, 18.67, 18.68, 18.69, 18.7, 18.71, 18.72, 18.73,
      18.740000000000002, 18.75, 18.76, 18.77, 18.78, 18.79, 18.8, 18.81, 18.82,
      18.830000000000002, 18.84, 18.85, 18.86, 18.87, 18.88, 18.89,
      18.900000000000002, 18.91, 18.92, 18.93, 18.94, 18.95, 18.96, 18.97, 18.98,
      18.990000000000002, 19.0, 19.01, 19.02, 19.03, 19.04, 19.05, 19.06, 19.07,
      19.080000000000002, 19.09, 19.1, 19.11, 19.12, 19.13, 19.14,
      19.150000000000002, 19.16, 19.17, 19.18, 19.19, 19.2, 19.21, 19.22, 19.23,
      19.240000000000002, 19.25, 19.26, 19.27, 19.28, 19.29, 19.3, 19.31, 19.32,
      19.330000000000002, 19.34, 19.35, 19.36, 19.37, 19.38, 19.39,
      19.400000000000002, 19.41, 19.42, 19.43, 19.44, 19.45, 19.46, 19.47, 19.48,
      19.490000000000002, 19.5, 19.51, 19.52, 19.53, 19.54, 19.55, 19.56, 19.57,
      19.580000000000002, 19.59, 19.6, 19.61, 19.62, 19.63, 19.64,
      19.650000000000002, 19.66, 19.67, 19.68, 19.69, 19.7, 19.71, 19.72, 19.73,
      19.740000000000002, 19.75, 19.76, 19.77, 19.78, 19.79, 19.8, 19.81, 19.82,
      19.830000000000002, 19.84, 19.85, 19.86, 19.87, 19.88, 19.89,
      19.900000000000002, 19.91, 19.92, 19.93, 19.94, 19.95, 19.96, 19.97, 19.98,
      19.990000000000002, 20.0, 20.01, 20.02, 20.03, 20.04, 20.05, 20.06, 20.07,
      20.080000000000002, 20.09, 20.1, 20.11, 20.12, 20.13, 20.14,
      20.150000000000002, 20.16, 20.17, 20.18, 20.19, 20.2, 20.21, 20.22, 20.23,
      20.240000000000002, 20.25, 20.26, 20.27, 20.28, 20.29, 20.3, 20.31, 20.32,
      20.330000000000002, 20.34, 20.35, 20.36, 20.37, 20.38, 20.39,
      20.400000000000002, 20.41, 20.42, 20.43, 20.44, 20.45, 20.46, 20.47, 20.48,
      20.490000000000002, 20.5, 20.51, 20.52, 20.53, 20.54, 20.55, 20.56, 20.57,
      20.580000000000002, 20.59, 20.6, 20.61, 20.62, 20.63, 20.64,
      20.650000000000002, 20.66, 20.67, 20.68, 20.69, 20.7, 20.71, 20.72, 20.73,
      20.740000000000002, 20.75, 20.76, 20.77, 20.78, 20.79, 20.8, 20.81, 20.82,
      20.830000000000002, 20.84, 20.85, 20.86, 20.87, 20.88, 20.89,
      20.900000000000002, 20.91, 20.92, 20.93, 20.94, 20.95, 20.96, 20.97, 20.98,
      20.990000000000002, 21.0, 21.01, 21.02, 21.03, 21.04, 21.05, 21.06, 21.07,
      21.080000000000002, 21.09, 21.1, 21.11, 21.12, 21.13, 21.14,
      21.150000000000002, 21.16, 21.17, 21.18, 21.19, 21.2, 21.21, 21.22, 21.23,
      21.240000000000002, 21.25, 21.26, 21.27, 21.28, 21.29, 21.3, 21.31, 21.32,
      21.330000000000002, 21.34, 21.35, 21.36, 21.37, 21.38, 21.39,
      21.400000000000002, 21.41, 21.42, 21.43, 21.44, 21.45, 21.46, 21.47, 21.48,
      21.490000000000002, 21.5, 21.51, 21.52, 21.53, 21.54, 21.55, 21.56, 21.57,
      21.580000000000002, 21.59, 21.6, 21.61, 21.62, 21.63, 21.64,
      21.650000000000002, 21.66, 21.67, 21.68, 21.69, 21.7, 21.71, 21.72, 21.73,
      21.740000000000002, 21.75, 21.76, 21.77, 21.78, 21.79, 21.8, 21.81, 21.82,
      21.830000000000002, 21.84, 21.85, 21.86, 21.87, 21.88, 21.89,
      21.900000000000002, 21.91, 21.92, 21.93, 21.94, 21.95, 21.96, 21.97, 21.98,
      21.990000000000002, 22.0, 22.01, 22.02, 22.03, 22.04, 22.05, 22.06, 22.07,
      22.080000000000002, 22.09, 22.1, 22.11, 22.12, 22.13, 22.14,
      22.150000000000002, 22.16, 22.17, 22.18, 22.19, 22.2, 22.21, 22.22, 22.23,
      22.240000000000002, 22.25, 22.26, 22.27, 22.28, 22.29, 22.3, 22.31, 22.32,
      22.330000000000002, 22.34, 22.35, 22.36, 22.37, 22.38, 22.39,
      22.400000000000002, 22.41, 22.42, 22.43, 22.44, 22.45, 22.46, 22.47, 22.48,
      22.490000000000002, 22.5, 22.51, 22.52, 22.53, 22.54, 22.55, 22.56, 22.57,
      22.580000000000002, 22.59, 22.6, 22.61, 22.62, 22.63, 22.64,
      22.650000000000002, 22.66, 22.67, 22.68, 22.69, 22.7, 22.71, 22.72, 22.73,
      22.740000000000002, 22.75, 22.76, 22.77, 22.78, 22.79, 22.8, 22.81, 22.82,
      22.830000000000002, 22.84, 22.85, 22.86, 22.87, 22.88, 22.89,
      22.900000000000002, 22.91, 22.92, 22.93, 22.94, 22.95, 22.96, 22.97, 22.98,
      22.990000000000002, 23.0, 23.01, 23.02, 23.03, 23.04, 23.05, 23.06, 23.07,
      23.080000000000002, 23.09, 23.1, 23.11, 23.12, 23.13, 23.14,
      23.150000000000002, 23.16, 23.17, 23.18, 23.19, 23.2, 23.21, 23.22, 23.23,
      23.240000000000002, 23.25, 23.26, 23.27, 23.28, 23.29, 23.3, 23.31, 23.32,
      23.330000000000002, 23.34, 23.35, 23.36, 23.37, 23.38, 23.39,
      23.400000000000002, 23.41, 23.42, 23.43, 23.44, 23.45, 23.46, 23.47, 23.48,
      23.490000000000002, 23.5, 23.51, 23.52, 23.53, 23.54, 23.55, 23.56, 23.57,
      23.580000000000002, 23.59, 23.6, 23.61, 23.62, 23.63, 23.64,
      23.650000000000002, 23.66, 23.67, 23.68, 23.69, 23.7, 23.71, 23.72, 23.73,
      23.740000000000002, 23.75, 23.76, 23.77, 23.78, 23.79, 23.8, 23.81, 23.82,
      23.830000000000002, 23.84, 23.85, 23.86, 23.87, 23.88, 23.89,
      23.900000000000002, 23.91, 23.92, 23.93, 23.94, 23.95, 23.96, 23.97, 23.98,
      23.990000000000002, 24.0, 24.01, 24.02, 24.03, 24.04, 24.05,
      24.060000000000002, 24.07, 24.080000000000002, 24.09, 24.1, 24.11, 24.12,
      24.13, 24.14, 24.150000000000002, 24.16, 24.17, 24.18, 24.19, 24.2, 24.21,
      24.22, 24.23, 24.240000000000002, 24.25, 24.26, 24.27, 24.28, 24.29, 24.3,
      24.310000000000002, 24.32, 24.330000000000002, 24.34, 24.35, 24.36, 24.37,
      24.38, 24.39, 24.400000000000002, 24.41, 24.42, 24.43, 24.44, 24.45, 24.46,
      24.47, 24.48, 24.490000000000002, 24.5, 24.51, 24.52, 24.53, 24.54, 24.55,
      24.560000000000002, 24.57, 24.580000000000002, 24.59, 24.6, 24.61, 24.62,
      24.63, 24.64, 24.650000000000002, 24.66, 24.67, 24.68, 24.69, 24.7, 24.71,
      24.72, 24.73, 24.740000000000002, 24.75, 24.76, 24.77, 24.78, 24.79, 24.8,
      24.810000000000002, 24.82, 24.830000000000002, 24.84, 24.85, 24.86, 24.87,
      24.88, 24.89, 24.900000000000002, 24.91, 24.92, 24.93, 24.94, 24.95, 24.96,
      24.97, 24.98, 24.990000000000002, 25.0, 25.01, 25.02, 25.03, 25.04, 25.05,
      25.060000000000002, 25.07, 25.080000000000002, 25.09, 25.1, 25.11, 25.12,
      25.13, 25.14, 25.150000000000002, 25.16, 25.17, 25.18, 25.19, 25.2, 25.21,
      25.22, 25.23, 25.240000000000002, 25.25, 25.26, 25.27, 25.28, 25.29, 25.3,
      25.310000000000002, 25.32, 25.330000000000002, 25.34, 25.35, 25.36, 25.37,
      25.38, 25.39, 25.400000000000002, 25.41, 25.42, 25.43, 25.44, 25.45, 25.46,
      25.47, 25.48, 25.490000000000002, 25.5, 25.51, 25.52, 25.53, 25.54, 25.55,
      25.560000000000002, 25.57, 25.580000000000002, 25.59, 25.6, 25.61, 25.62,
      25.63, 25.64, 25.650000000000002, 25.66, 25.67, 25.68, 25.69, 25.7, 25.71,
      25.72, 25.73, 25.740000000000002, 25.75, 25.76, 25.77, 25.78, 25.79, 25.8,
      25.810000000000002, 25.82, 25.830000000000002, 25.84, 25.85, 25.86, 25.87,
      25.88, 25.89, 25.900000000000002, 25.91, 25.92, 25.93, 25.94, 25.95, 25.96,
      25.97, 25.98, 25.990000000000002, 26.0, 26.01, 26.02, 26.03, 26.04, 26.05,
      26.060000000000002, 26.07, 26.080000000000002, 26.09, 26.1, 26.11, 26.12,
      26.13, 26.14, 26.150000000000002, 26.16, 26.17, 26.18, 26.19, 26.2, 26.21,
      26.22, 26.23, 26.240000000000002, 26.25, 26.26, 26.27, 26.28, 26.29, 26.3,
      26.310000000000002, 26.32, 26.330000000000002, 26.34, 26.35, 26.36, 26.37,
      26.38, 26.39, 26.400000000000002, 26.41, 26.42, 26.43, 26.44, 26.45, 26.46,
      26.47, 26.48, 26.490000000000002, 26.5, 26.51, 26.52, 26.53, 26.54, 26.55,
      26.560000000000002, 26.57, 26.580000000000002, 26.59, 26.6, 26.61, 26.62,
      26.63, 26.64, 26.650000000000002, 26.66, 26.67, 26.68, 26.69, 26.7, 26.71,
      26.72, 26.73, 26.740000000000002, 26.75, 26.76, 26.77, 26.78, 26.79, 26.8,
      26.810000000000002, 26.82, 26.830000000000002, 26.84, 26.85, 26.86, 26.87,
      26.88, 26.89, 26.900000000000002, 26.91, 26.92, 26.93, 26.94, 26.95, 26.96,
      26.97, 26.98, 26.990000000000002, 27.0, 27.01, 27.02, 27.03, 27.04, 27.05,
      27.060000000000002, 27.07, 27.080000000000002, 27.09, 27.1, 27.11, 27.12,
      27.13, 27.14, 27.150000000000002, 27.16, 27.17, 27.18, 27.19, 27.2, 27.21,
      27.22, 27.23, 27.240000000000002, 27.25, 27.26, 27.27, 27.28, 27.29, 27.3,
      27.310000000000002, 27.32, 27.330000000000002, 27.34, 27.35, 27.36, 27.37,
      27.38, 27.39, 27.400000000000002, 27.41, 27.42, 27.43, 27.44, 27.45, 27.46,
      27.47, 27.48, 27.490000000000002, 27.5, 27.51, 27.52, 27.53, 27.54, 27.55,
      27.560000000000002, 27.57, 27.580000000000002, 27.59, 27.6, 27.61, 27.62,
      27.63, 27.64, 27.650000000000002, 27.66, 27.67, 27.68, 27.69, 27.7, 27.71,
      27.72, 27.73, 27.740000000000002, 27.75, 27.76, 27.77, 27.78, 27.79, 27.8,
      27.810000000000002, 27.82, 27.830000000000002, 27.84, 27.85, 27.86, 27.87,
      27.88, 27.89, 27.900000000000002, 27.91, 27.92, 27.93, 27.94, 27.95, 27.96,
      27.97, 27.98, 27.990000000000002, 28.0, 28.01, 28.02, 28.03, 28.04, 28.05,
      28.060000000000002, 28.07, 28.080000000000002, 28.09, 28.1, 28.11, 28.12,
      28.13, 28.14, 28.150000000000002, 28.16, 28.17, 28.18, 28.19, 28.2, 28.21,
      28.22, 28.23, 28.240000000000002, 28.25, 28.26, 28.27, 28.28, 28.29, 28.3,
      28.310000000000002, 28.32, 28.330000000000002, 28.34, 28.35, 28.36, 28.37,
      28.38, 28.39, 28.400000000000002, 28.41, 28.42, 28.43, 28.44, 28.45, 28.46,
      28.47, 28.48, 28.490000000000002, 28.5, 28.51, 28.52, 28.53, 28.54, 28.55,
      28.560000000000002, 28.57, 28.580000000000002, 28.59, 28.6, 28.61, 28.62,
      28.63, 28.64, 28.650000000000002, 28.66, 28.67, 28.68, 28.69, 28.7, 28.71,
      28.72, 28.73, 28.740000000000002, 28.75, 28.76, 28.77, 28.78, 28.79, 28.8,
      28.810000000000002, 28.82, 28.830000000000002, 28.84, 28.85, 28.86, 28.87,
      28.88, 28.89, 28.900000000000002, 28.91, 28.92, 28.93, 28.94, 28.95, 28.96,
      28.97, 28.98, 28.990000000000002, 29.0, 29.01, 29.02, 29.03, 29.04, 29.05,
      29.060000000000002, 29.07, 29.080000000000002, 29.09, 29.1, 29.11, 29.12,
      29.13, 29.14, 29.150000000000002, 29.16, 29.17, 29.18, 29.19, 29.2, 29.21,
      29.22, 29.23, 29.240000000000002, 29.25, 29.26, 29.27, 29.28, 29.29, 29.3,
      29.310000000000002, 29.32, 29.330000000000002, 29.34, 29.35, 29.36, 29.37,
      29.38, 29.39, 29.400000000000002, 29.41, 29.42, 29.43, 29.44, 29.45, 29.46,
      29.47, 29.48, 29.490000000000002, 29.5, 29.51, 29.52, 29.53, 29.54, 29.55,
      29.560000000000002, 29.57, 29.580000000000002, 29.59, 29.6, 29.61, 29.62,
      29.63, 29.64, 29.650000000000002, 29.66, 29.67, 29.68, 29.69, 29.7, 29.71,
      29.72, 29.73, 29.740000000000002, 29.75, 29.76, 29.77, 29.78, 29.79, 29.8,
      29.810000000000002, 29.82, 29.830000000000002, 29.84, 29.85, 29.86, 29.87,
      29.88, 29.89, 29.900000000000002, 29.91, 29.92, 29.93, 29.94, 29.95, 29.96,
      29.97, 29.98, 29.990000000000002, 30.0, 30.01, 30.02, 30.03, 30.04, 30.05,
      30.060000000000002, 30.07, 30.080000000000002, 30.09, 30.1, 30.11, 30.12,
      30.13, 30.14, 30.150000000000002, 30.16, 30.17, 30.18, 30.19, 30.2, 30.21,
      30.22, 30.23, 30.240000000000002, 30.25, 30.26, 30.27, 30.28, 30.29, 30.3,
      30.310000000000002, 30.32, 30.330000000000002, 30.34, 30.35, 30.36, 30.37,
      30.38, 30.39, 30.400000000000002, 30.41, 30.42, 30.43, 30.44, 30.45, 30.46,
      30.47, 30.48, 30.490000000000002, 30.5, 30.51, 30.52, 30.53, 30.54, 30.55,
      30.560000000000002, 30.57, 30.580000000000002, 30.59, 30.6, 30.61, 30.62,
      30.63, 30.64, 30.650000000000002, 30.66, 30.67, 30.68, 30.69, 30.7, 30.71,
      30.72, 30.73, 30.740000000000002, 30.75, 30.76, 30.77, 30.78, 30.79, 30.8,
      30.810000000000002, 30.82, 30.830000000000002, 30.84, 30.85, 30.86, 30.87,
      30.88, 30.89, 30.900000000000002, 30.91, 30.92, 30.93, 30.94, 30.95, 30.96,
      30.970000000000002, 30.98, 30.990000000000002, 31.0, 31.01, 31.02, 31.03,
      31.04, 31.05, 31.060000000000002, 31.07, 31.080000000000002, 31.09, 31.1,
      31.11, 31.12, 31.13, 31.14, 31.150000000000002, 31.16, 31.17, 31.18, 31.19,
      31.2, 31.21, 31.220000000000002, 31.23, 31.240000000000002, 31.25, 31.26,
      31.27, 31.28, 31.29, 31.3, 31.310000000000002, 31.32, 31.330000000000002,
      31.34, 31.35, 31.36, 31.37, 31.38, 31.39, 31.400000000000002, 31.41, 31.42,
      31.43, 31.44, 31.45, 31.46, 31.470000000000002, 31.48, 31.490000000000002,
      31.5, 31.51, 31.52, 31.53, 31.54, 31.55, 31.560000000000002, 31.57,
      31.580000000000002, 31.59, 31.6, 31.61, 31.62, 31.63, 31.64,
      31.650000000000002, 31.66, 31.67, 31.68, 31.69, 31.7, 31.71,
      31.720000000000002, 31.73, 31.740000000000002, 31.75, 31.76, 31.77, 31.78,
      31.79, 31.8, 31.810000000000002, 31.82, 31.830000000000002, 31.84, 31.85,
      31.86, 31.87, 31.88, 31.89, 31.900000000000002, 31.91, 31.92, 31.93, 31.94,
      31.95, 31.96, 31.970000000000002, 31.98, 31.990000000000002, 32.0, 32.01,
      32.02, 32.03, 32.04, 32.05, 32.06, 32.07, 32.08, 32.09, 32.1, 32.11, 32.12,
      32.13, 32.14, 32.15, 32.160000000000004, 32.17, 32.18, 32.19, 32.2, 32.21,
      32.22, 32.230000000000004, 32.24, 32.25, 32.26, 32.27, 32.28, 32.29, 32.3,
      32.31, 32.32, 32.33, 32.34, 32.35, 32.36, 32.37, 32.38, 32.39, 32.4,
      32.410000000000004, 32.42, 32.43, 32.44, 32.45, 32.46, 32.47,
      32.480000000000004, 32.49, 32.5, 32.51, 32.52, 32.53, 32.54, 32.55, 32.56,
      32.57, 32.58, 32.59, 32.6, 32.61, 32.62, 32.63, 32.64, 32.65,
      32.660000000000004, 32.67, 32.68, 32.69, 32.7, 32.71, 32.72,
      32.730000000000004, 32.74, 32.75, 32.76, 32.77, 32.78, 32.79, 32.8, 32.81,
      32.82, 32.83, 32.84, 32.85, 32.86, 32.87, 32.88, 32.89, 32.9,
      32.910000000000004, 32.92, 32.93, 32.94, 32.95, 32.96, 32.97,
      32.980000000000004, 32.99, 33.0, 33.01, 33.02, 33.03, 33.04, 33.05, 33.06,
      33.07, 33.08, 33.09, 33.1, 33.11, 33.12, 33.13, 33.14, 33.15,
      33.160000000000004, 33.17, 33.18, 33.19, 33.2, 33.21, 33.22,
      33.230000000000004, 33.24, 33.25, 33.26, 33.27, 33.28, 33.29, 33.3, 33.31,
      33.32, 33.33, 33.34, 33.35, 33.36, 33.37, 33.38, 33.39, 33.4,
      33.410000000000004, 33.42, 33.43, 33.44, 33.45, 33.46, 33.47,
      33.480000000000004, 33.49, 33.5, 33.51, 33.52, 33.53, 33.54, 33.55, 33.56,
      33.57, 33.58, 33.59, 33.6, 33.61, 33.62, 33.63, 33.64, 33.65,
      33.660000000000004, 33.67, 33.68, 33.69, 33.7, 33.71, 33.72,
      33.730000000000004, 33.74, 33.75, 33.76, 33.77, 33.78, 33.79, 33.8, 33.81,
      33.82, 33.83, 33.84, 33.85, 33.86, 33.87, 33.88, 33.89, 33.9,
      33.910000000000004, 33.92, 33.93, 33.94, 33.95, 33.96, 33.97,
      33.980000000000004, 33.99, 34.0, 34.01, 34.02, 34.03, 34.04, 34.05, 34.06,
      34.07, 34.08, 34.09, 34.1, 34.11, 34.12, 34.13, 34.14, 34.15,
      34.160000000000004, 34.17, 34.18, 34.19, 34.2, 34.21, 34.22,
      34.230000000000004, 34.24, 34.25, 34.26, 34.27, 34.28, 34.29,
      34.300000000000004, 34.31, 34.32, 34.33, 34.34, 34.35, 34.36, 34.37, 34.38,
      34.39, 34.4, 34.410000000000004, 34.42, 34.43, 34.44, 34.45, 34.46, 34.47,
      34.480000000000004, 34.49, 34.5, 34.51, 34.52, 34.53, 34.54,
      34.550000000000004, 34.56, 34.57, 34.58, 34.59, 34.6, 34.61, 34.62, 34.63,
      34.64, 34.65, 34.660000000000004, 34.67, 34.68, 34.69, 34.7, 34.71, 34.72,
      34.730000000000004, 34.74, 34.75, 34.76, 34.77, 34.78, 34.79,
      34.800000000000004, 34.81, 34.82, 34.83, 34.84, 34.85, 34.86, 34.87, 34.88,
      34.89, 34.9, 34.910000000000004, 34.92, 34.93, 34.94, 34.95, 34.96, 34.97,
      34.980000000000004, 34.99, 35.0, 35.01, 35.02, 35.03, 35.04,
      35.050000000000004, 35.06, 35.07, 35.08, 35.09, 35.1, 35.11, 35.12, 35.13,
      35.14, 35.15, 35.160000000000004, 35.17, 35.18, 35.19, 35.2, 35.21, 35.22,
      35.230000000000004, 35.24, 35.25, 35.26, 35.27, 35.28, 35.29,
      35.300000000000004, 35.31, 35.32, 35.33, 35.34, 35.35, 35.36, 35.37, 35.38,
      35.39, 35.4, 35.410000000000004, 35.42, 35.43, 35.44, 35.45, 35.46, 35.47,
      35.480000000000004, 35.49, 35.5, 35.51, 35.52, 35.53, 35.54,
      35.550000000000004, 35.56, 35.57, 35.58, 35.59, 35.6, 35.61, 35.62, 35.63,
      35.64, 35.65, 35.660000000000004, 35.67, 35.68, 35.69, 35.7, 35.71, 35.72,
      35.730000000000004, 35.74, 35.75, 35.76, 35.77, 35.78, 35.79,
      35.800000000000004, 35.81, 35.82, 35.83, 35.84, 35.85, 35.86, 35.87, 35.88,
      35.89, 35.9, 35.910000000000004, 35.92, 35.93, 35.94, 35.95, 35.96, 35.97,
      35.980000000000004, 35.99, 36.0, 36.01, 36.02, 36.03, 36.04,
      36.050000000000004, 36.06, 36.07, 36.08, 36.09, 36.1, 36.11, 36.12, 36.13,
      36.14, 36.15, 36.160000000000004, 36.17, 36.18, 36.19, 36.2, 36.21, 36.22,
      36.230000000000004, 36.24, 36.25, 36.26, 36.27, 36.28, 36.29,
      36.300000000000004, 36.31, 36.32, 36.33, 36.34, 36.35, 36.36, 36.37, 36.38,
      36.39, 36.4, 36.410000000000004, 36.42, 36.43, 36.44, 36.45, 36.46, 36.47,
      36.480000000000004, 36.49, 36.5, 36.51, 36.52, 36.53, 36.54,
      36.550000000000004, 36.56, 36.57, 36.58, 36.59, 36.6, 36.61, 36.62, 36.63,
      36.64, 36.65, 36.660000000000004, 36.67, 36.68, 36.69, 36.7, 36.71, 36.72,
      36.730000000000004, 36.74, 36.75, 36.76, 36.77, 36.78, 36.79,
      36.800000000000004, 36.81, 36.82, 36.83, 36.84, 36.85, 36.86, 36.87, 36.88,
      36.89, 36.9, 36.910000000000004, 36.92, 36.93, 36.94, 36.95, 36.96, 36.97,
      36.980000000000004, 36.99, 37.0, 37.01, 37.02, 37.03, 37.04,
      37.050000000000004, 37.06, 37.07, 37.08, 37.09, 37.1, 37.11, 37.12, 37.13,
      37.14, 37.15, 37.160000000000004, 37.17, 37.18, 37.19, 37.2, 37.21, 37.22,
      37.230000000000004, 37.24, 37.25, 37.26, 37.27, 37.28, 37.29,
      37.300000000000004, 37.31, 37.32, 37.33, 37.34, 37.35, 37.36, 37.37, 37.38,
      37.39, 37.4, 37.410000000000004, 37.42, 37.43, 37.44, 37.45, 37.46, 37.47,
      37.480000000000004, 37.49, 37.5, 37.51, 37.52, 37.53, 37.54,
      37.550000000000004, 37.56, 37.57, 37.58, 37.59, 37.6, 37.61, 37.62, 37.63,
      37.64, 37.65, 37.660000000000004, 37.67, 37.68, 37.69, 37.7, 37.71, 37.72,
      37.730000000000004, 37.74, 37.75, 37.76, 37.77, 37.78, 37.79,
      37.800000000000004, 37.81, 37.82, 37.83, 37.84, 37.85, 37.86, 37.87, 37.88,
      37.89, 37.9, 37.910000000000004, 37.92, 37.93, 37.94, 37.95, 37.96, 37.97,
      37.980000000000004, 37.99, 38.0, 38.01, 38.02, 38.03, 38.04,
      38.050000000000004, 38.06, 38.07, 38.08, 38.09, 38.1, 38.11, 38.12, 38.13,
      38.14, 38.15, 38.160000000000004, 38.17, 38.18, 38.19, 38.2, 38.21, 38.22,
      38.230000000000004, 38.24, 38.25, 38.26, 38.27, 38.28, 38.29,
      38.300000000000004, 38.31, 38.32, 38.33, 38.34, 38.35, 38.36, 38.37, 38.38,
      38.39, 38.4, 38.410000000000004, 38.42, 38.43, 38.44, 38.45, 38.46, 38.47,
      38.480000000000004, 38.49, 38.5, 38.51, 38.52, 38.53, 38.54,
      38.550000000000004, 38.56, 38.57, 38.58, 38.59, 38.6, 38.61, 38.62, 38.63,
      38.64, 38.65, 38.660000000000004, 38.67, 38.68, 38.69, 38.7, 38.71, 38.72,
      38.730000000000004, 38.74, 38.75, 38.76, 38.77, 38.78, 38.79,
      38.800000000000004, 38.81, 38.82, 38.83, 38.84, 38.85, 38.86, 38.87, 38.88,
      38.89, 38.9, 38.910000000000004, 38.92, 38.93, 38.94, 38.95, 38.96, 38.97,
      38.980000000000004, 38.99, 39.0, 39.01, 39.02, 39.03, 39.04,
      39.050000000000004, 39.06, 39.07, 39.08, 39.09, 39.1, 39.11, 39.12, 39.13,
      39.14, 39.15, 39.160000000000004, 39.17, 39.18, 39.19, 39.2, 39.21, 39.22,
      39.230000000000004, 39.24, 39.25, 39.26, 39.27, 39.28, 39.29,
      39.300000000000004, 39.31, 39.32, 39.33, 39.34, 39.35, 39.36, 39.37, 39.38,
      39.39, 39.4, 39.410000000000004, 39.42, 39.43, 39.44, 39.45, 39.46, 39.47,
      39.480000000000004, 39.49, 39.5, 39.51, 39.52, 39.53, 39.54,
      39.550000000000004, 39.56, 39.57, 39.58, 39.59, 39.6, 39.61, 39.62, 39.63,
      39.64, 39.65, 39.660000000000004, 39.67, 39.68, 39.69, 39.7, 39.71, 39.72,
      39.730000000000004, 39.74, 39.75, 39.76, 39.77, 39.78, 39.79,
      39.800000000000004, 39.81, 39.82, 39.83, 39.84, 39.85, 39.86, 39.87, 39.88,
      39.89, 39.9, 39.910000000000004, 39.92, 39.93, 39.94, 39.95, 39.96, 39.97,
      39.980000000000004, 39.99, 40.0, 40.01, 40.02, 40.03, 40.04,
      40.050000000000004, 40.06, 40.07, 40.08, 40.09, 40.1, 40.11, 40.12, 40.13,
      40.14, 40.15, 40.160000000000004, 40.17, 40.18, 40.19, 40.2, 40.21, 40.22,
      40.230000000000004, 40.24, 40.25, 40.26, 40.27, 40.28, 40.29,
      40.300000000000004, 40.31, 40.32, 40.33, 40.34, 40.35, 40.36, 40.37, 40.38,
      40.39, 40.4, 40.410000000000004, 40.42, 40.43, 40.44, 40.45, 40.46, 40.47,
      40.480000000000004, 40.49, 40.5, 40.51, 40.52, 40.53, 40.54,
      40.550000000000004, 40.56, 40.57, 40.58, 40.59, 40.6, 40.61, 40.62, 40.63,
      40.64, 40.65, 40.660000000000004, 40.67, 40.68, 40.69, 40.7, 40.71, 40.72,
      40.730000000000004, 40.74, 40.75, 40.76, 40.77, 40.78, 40.79,
      40.800000000000004, 40.81, 40.82, 40.83, 40.84, 40.85, 40.86, 40.87, 40.88,
      40.89, 40.9, 40.910000000000004, 40.92, 40.93, 40.94, 40.95, 40.96, 40.97,
      40.980000000000004, 40.99, 41.0, 41.01, 41.02, 41.03, 41.04,
      41.050000000000004, 41.06, 41.07, 41.08, 41.09, 41.1, 41.11, 41.12, 41.13,
      41.14, 41.15, 41.160000000000004, 41.17, 41.18, 41.19, 41.2, 41.21, 41.22,
      41.230000000000004, 41.24, 41.25, 41.26, 41.27, 41.28, 41.29,
      41.300000000000004, 41.31, 41.32, 41.33, 41.34, 41.35, 41.36, 41.37, 41.38,
      41.39, 41.4, 41.410000000000004, 41.42, 41.43, 41.44, 41.45, 41.46, 41.47,
      41.480000000000004, 41.49, 41.5, 41.51, 41.52, 41.53, 41.54,
      41.550000000000004, 41.56, 41.57, 41.58, 41.59, 41.6, 41.61, 41.62, 41.63,
      41.64, 41.65, 41.660000000000004, 41.67, 41.68, 41.69, 41.7, 41.71, 41.72,
      41.730000000000004, 41.74, 41.75, 41.76, 41.77, 41.78, 41.79,
      41.800000000000004, 41.81, 41.82, 41.83, 41.84, 41.85, 41.86, 41.87, 41.88,
      41.89, 41.9, 41.910000000000004, 41.92, 41.93, 41.94, 41.95, 41.96, 41.97,
      41.980000000000004, 41.99, 42.0, 42.01, 42.02, 42.03, 42.04,
      42.050000000000004, 42.06, 42.07, 42.08, 42.09, 42.1, 42.11, 42.12, 42.13,
      42.14, 42.15, 42.160000000000004, 42.17, 42.18, 42.19, 42.2, 42.21, 42.22,
      42.230000000000004, 42.24, 42.25, 42.26, 42.27, 42.28, 42.29,
      42.300000000000004, 42.31, 42.32, 42.33, 42.34, 42.35, 42.36, 42.37, 42.38,
      42.39, 42.4, 42.410000000000004, 42.42, 42.43, 42.44, 42.45, 42.46, 42.47,
      42.480000000000004, 42.49, 42.5, 42.51, 42.52, 42.53, 42.54,
      42.550000000000004, 42.56, 42.57, 42.58, 42.59, 42.6, 42.61, 42.62, 42.63,
      42.64, 42.65, 42.660000000000004, 42.67, 42.68, 42.69, 42.7, 42.71, 42.72,
      42.730000000000004, 42.74, 42.75, 42.76, 42.77, 42.78, 42.79,
      42.800000000000004, 42.81, 42.82, 42.83, 42.84, 42.85, 42.86, 42.87, 42.88,
      42.89, 42.9, 42.910000000000004, 42.92, 42.93, 42.94, 42.95, 42.96, 42.97,
      42.980000000000004, 42.99, 43.0, 43.01, 43.02, 43.03, 43.04,
      43.050000000000004, 43.06, 43.07, 43.08, 43.09, 43.1, 43.11, 43.12, 43.13,
      43.14, 43.15, 43.160000000000004, 43.17, 43.18, 43.19, 43.2, 43.21, 43.22,
      43.230000000000004, 43.24, 43.25, 43.26, 43.27, 43.28, 43.29,
      43.300000000000004, 43.31, 43.32, 43.33, 43.34, 43.35, 43.36, 43.37, 43.38,
      43.39, 43.4, 43.410000000000004, 43.42, 43.43, 43.44, 43.45, 43.46, 43.47,
      43.480000000000004, 43.49, 43.5, 43.51, 43.52, 43.53, 43.54,
      43.550000000000004, 43.56, 43.57, 43.58, 43.59, 43.6, 43.61, 43.62, 43.63,
      43.64, 43.65, 43.660000000000004, 43.67, 43.68, 43.69, 43.7, 43.71, 43.72,
      43.730000000000004, 43.74, 43.75, 43.76, 43.77, 43.78, 43.79,
      43.800000000000004, 43.81, 43.82, 43.83, 43.84, 43.85, 43.86, 43.87, 43.88,
      43.89, 43.9, 43.910000000000004, 43.92, 43.93, 43.94, 43.95, 43.96, 43.97,
      43.980000000000004, 43.99, 44.0, 44.01, 44.02, 44.03, 44.04,
      44.050000000000004, 44.06, 44.07, 44.08, 44.09, 44.1, 44.11, 44.12, 44.13,
      44.14, 44.15, 44.160000000000004, 44.17, 44.18, 44.19, 44.2, 44.21, 44.22,
      44.230000000000004, 44.24, 44.25, 44.26, 44.27, 44.28, 44.29,
      44.300000000000004, 44.31, 44.32, 44.33, 44.34, 44.35, 44.36, 44.37, 44.38,
      44.39, 44.4, 44.410000000000004, 44.42, 44.43, 44.44, 44.45, 44.46, 44.47,
      44.480000000000004, 44.49, 44.5, 44.51, 44.52, 44.53, 44.54,
      44.550000000000004, 44.56, 44.57, 44.58, 44.59, 44.6, 44.61, 44.62, 44.63,
      44.64, 44.65, 44.660000000000004, 44.67, 44.68, 44.69, 44.7, 44.71, 44.72,
      44.730000000000004, 44.74, 44.75, 44.76, 44.77, 44.78, 44.79,
      44.800000000000004, 44.81, 44.82, 44.83, 44.84, 44.85, 44.86, 44.87, 44.88,
      44.89, 44.9, 44.910000000000004, 44.92, 44.93, 44.94, 44.95, 44.96, 44.97,
      44.980000000000004, 44.99, 45.0, 45.01, 45.02, 45.03, 45.04,
      45.050000000000004, 45.06, 45.07, 45.08, 45.09, 45.1, 45.11, 45.12, 45.13,
      45.14, 45.15, 45.160000000000004, 45.17, 45.18, 45.19, 45.2, 45.21, 45.22,
      45.230000000000004, 45.24, 45.25, 45.26, 45.27, 45.28, 45.29,
      45.300000000000004, 45.31, 45.32, 45.33, 45.34, 45.35, 45.36, 45.37, 45.38,
      45.39, 45.4, 45.410000000000004, 45.42, 45.43, 45.44, 45.45, 45.46, 45.47,
      45.480000000000004, 45.49, 45.5, 45.51, 45.52, 45.53, 45.54,
      45.550000000000004, 45.56, 45.57, 45.58, 45.59, 45.6, 45.61, 45.62, 45.63,
      45.64, 45.65, 45.660000000000004, 45.67, 45.68, 45.69, 45.7, 45.71, 45.72,
      45.730000000000004, 45.74, 45.75, 45.76, 45.77, 45.78, 45.79,
      45.800000000000004, 45.81, 45.82, 45.83, 45.84, 45.85, 45.86, 45.87, 45.88,
      45.89, 45.9, 45.910000000000004, 45.92, 45.93, 45.94, 45.95, 45.96, 45.97,
      45.980000000000004, 45.99, 46.0, 46.01, 46.02, 46.03, 46.04,
      46.050000000000004, 46.06, 46.07, 46.08, 46.09, 46.1, 46.11, 46.12, 46.13,
      46.14, 46.15, 46.160000000000004, 46.17, 46.18, 46.19, 46.2, 46.21, 46.22,
      46.230000000000004, 46.24, 46.25, 46.26, 46.27, 46.28, 46.29,
      46.300000000000004, 46.31, 46.32, 46.33, 46.34, 46.35, 46.36, 46.37, 46.38,
      46.39, 46.4, 46.410000000000004, 46.42, 46.43, 46.44, 46.45, 46.46, 46.47,
      46.480000000000004, 46.49, 46.5, 46.51, 46.52, 46.53, 46.54,
      46.550000000000004, 46.56, 46.57, 46.58, 46.59, 46.6, 46.61, 46.62, 46.63,
      46.64, 46.65, 46.660000000000004, 46.67, 46.68, 46.69, 46.7, 46.71, 46.72,
      46.730000000000004, 46.74, 46.75, 46.76, 46.77, 46.78, 46.79,
      46.800000000000004, 46.81, 46.82, 46.83, 46.84, 46.85, 46.86, 46.87, 46.88,
      46.89, 46.9, 46.910000000000004, 46.92, 46.93, 46.94, 46.95, 46.96, 46.97,
      46.980000000000004, 46.99, 47.0, 47.01, 47.02, 47.03, 47.04,
      47.050000000000004, 47.06, 47.07, 47.08, 47.09, 47.1, 47.11, 47.12, 47.13,
      47.14, 47.15, 47.160000000000004, 47.17, 47.18, 47.19, 47.2, 47.21, 47.22,
      47.230000000000004, 47.24, 47.25, 47.26, 47.27, 47.28, 47.29,
      47.300000000000004, 47.31, 47.32, 47.33, 47.34, 47.35, 47.36, 47.37, 47.38,
      47.39, 47.4, 47.410000000000004, 47.42, 47.43, 47.44, 47.45, 47.46, 47.47,
      47.480000000000004, 47.49, 47.5, 47.51, 47.52, 47.53, 47.54,
      47.550000000000004, 47.56, 47.57, 47.58, 47.59, 47.6, 47.61, 47.62, 47.63,
      47.64, 47.65, 47.660000000000004, 47.67, 47.68, 47.69, 47.7, 47.71, 47.72,
      47.730000000000004, 47.74, 47.75, 47.76, 47.77, 47.78, 47.79,
      47.800000000000004, 47.81, 47.82, 47.83, 47.84, 47.85, 47.86,
      47.870000000000005, 47.88, 47.89, 47.9, 47.910000000000004, 47.92, 47.93,
      47.94, 47.95, 47.96, 47.97, 47.980000000000004, 47.99, 48.0, 48.01, 48.02,
      48.03, 48.04, 48.050000000000004, 48.06, 48.07, 48.08, 48.09, 48.1, 48.11,
      48.120000000000005, 48.13, 48.14, 48.15, 48.160000000000004, 48.17, 48.18,
      48.19, 48.2, 48.21, 48.22, 48.230000000000004, 48.24, 48.25, 48.26, 48.27,
      48.28, 48.29, 48.300000000000004, 48.31, 48.32, 48.33, 48.34, 48.35, 48.36,
      48.370000000000005, 48.38, 48.39, 48.4, 48.410000000000004, 48.42, 48.43,
      48.44, 48.45, 48.46, 48.47, 48.480000000000004, 48.49, 48.5, 48.51, 48.52,
      48.53, 48.54, 48.550000000000004, 48.56, 48.57, 48.58, 48.59, 48.6, 48.61,
      48.620000000000005, 48.63, 48.64, 48.65, 48.660000000000004, 48.67, 48.68,
      48.69, 48.7, 48.71, 48.72, 48.730000000000004, 48.74, 48.75, 48.76, 48.77,
      48.78, 48.79, 48.800000000000004, 48.81, 48.82, 48.83, 48.84, 48.85, 48.86,
      48.870000000000005, 48.88, 48.89, 48.9, 48.910000000000004, 48.92, 48.93,
      48.94, 48.95, 48.96, 48.97, 48.980000000000004, 48.99, 49.0, 49.01, 49.02,
      49.03, 49.04, 49.050000000000004, 49.06, 49.07, 49.08, 49.09, 49.1, 49.11,
      49.120000000000005, 49.13, 49.14, 49.15, 49.160000000000004, 49.17, 49.18,
      49.19, 49.2, 49.21, 49.22, 49.230000000000004, 49.24, 49.25, 49.26, 49.27,
      49.28, 49.29, 49.300000000000004, 49.31, 49.32, 49.33, 49.34, 49.35, 49.36,
      49.370000000000005, 49.38, 49.39, 49.4, 49.410000000000004, 49.42, 49.43,
      49.44, 49.45, 49.46, 49.47, 49.480000000000004, 49.49, 49.5, 49.51, 49.52,
      49.53, 49.54, 49.550000000000004, 49.56, 49.57, 49.58, 49.59, 49.6, 49.61,
      49.620000000000005, 49.63, 49.64, 49.65, 49.660000000000004, 49.67, 49.68,
      49.69, 49.7, 49.71, 49.72, 49.730000000000004, 49.74, 49.75, 49.76, 49.77,
      49.78, 49.79, 49.800000000000004, 49.81, 49.82, 49.83, 49.84, 49.85, 49.86,
      49.870000000000005, 49.88, 49.89, 49.9, 49.910000000000004, 49.92, 49.93,
      49.94, 49.95, 49.96, 49.97, 49.980000000000004, 49.99, 50.0, 50.01, 50.02,
      50.03, 50.04, 50.050000000000004, 50.06, 50.07, 50.08, 50.09, 50.1, 50.11,
      50.120000000000005, 50.13, 50.14, 50.15, 50.160000000000004, 50.17, 50.18,
      50.19, 50.2, 50.21, 50.22, 50.230000000000004, 50.24, 50.25, 50.26, 50.27,
      50.28, 50.29, 50.300000000000004, 50.31, 50.32, 50.33, 50.34, 50.35, 50.36,
      50.370000000000005, 50.38, 50.39, 50.4, 50.410000000000004, 50.42, 50.43,
      50.44, 50.45, 50.46, 50.47, 50.480000000000004, 50.49, 50.5, 50.51, 50.52,
      50.53, 50.54, 50.550000000000004, 50.56, 50.57, 50.58, 50.59, 50.6, 50.61,
      50.620000000000005, 50.63, 50.64, 50.65, 50.660000000000004, 50.67, 50.68,
      50.69, 50.7, 50.71, 50.72, 50.730000000000004, 50.74, 50.75, 50.76, 50.77,
      50.78, 50.79, 50.800000000000004, 50.81, 50.82, 50.83, 50.84, 50.85, 50.86,
      50.870000000000005, 50.88, 50.89, 50.9, 50.910000000000004, 50.92, 50.93,
      50.94, 50.95, 50.96, 50.97, 50.980000000000004, 50.99, 51.0, 51.01, 51.02,
      51.03, 51.04, 51.050000000000004, 51.06, 51.07, 51.08, 51.09, 51.1, 51.11,
      51.120000000000005, 51.13, 51.14, 51.15, 51.160000000000004, 51.17, 51.18,
      51.19, 51.2, 51.21, 51.22, 51.230000000000004, 51.24, 51.25, 51.26, 51.27,
      51.28, 51.29, 51.300000000000004, 51.31, 51.32, 51.33, 51.34, 51.35, 51.36,
      51.370000000000005, 51.38, 51.39, 51.4, 51.410000000000004, 51.42, 51.43,
      51.44, 51.45, 51.46, 51.47, 51.480000000000004, 51.49, 51.5, 51.51, 51.52,
      51.53, 51.54, 51.550000000000004, 51.56, 51.57, 51.58, 51.59, 51.6, 51.61,
      51.620000000000005, 51.63, 51.64, 51.65, 51.660000000000004, 51.67, 51.68,
      51.69, 51.7, 51.71, 51.72, 51.730000000000004, 51.74, 51.75, 51.76, 51.77,
      51.78, 51.79, 51.800000000000004, 51.81, 51.82, 51.83, 51.84, 51.85, 51.86,
      51.870000000000005, 51.88, 51.89, 51.9, 51.910000000000004, 51.92, 51.93,
      51.94, 51.95, 51.96, 51.97, 51.980000000000004, 51.99, 52.0, 52.01, 52.02,
      52.03, 52.04, 52.050000000000004, 52.06, 52.07, 52.08, 52.09, 52.1, 52.11,
      52.120000000000005, 52.13, 52.14, 52.15, 52.160000000000004, 52.17, 52.18,
      52.19, 52.2, 52.21, 52.22, 52.230000000000004, 52.24, 52.25, 52.26, 52.27,
      52.28, 52.29, 52.300000000000004, 52.31, 52.32, 52.33, 52.34, 52.35, 52.36,
      52.370000000000005, 52.38, 52.39, 52.4, 52.410000000000004, 52.42, 52.43,
      52.44, 52.45, 52.46, 52.47, 52.480000000000004, 52.49, 52.5, 52.51, 52.52,
      52.53, 52.54, 52.550000000000004, 52.56, 52.57, 52.58, 52.59, 52.6, 52.61,
      52.620000000000005, 52.63, 52.64, 52.65, 52.660000000000004, 52.67, 52.68,
      52.69, 52.7, 52.71, 52.72, 52.730000000000004, 52.74, 52.75, 52.76, 52.77,
      52.78, 52.79, 52.800000000000004, 52.81, 52.82, 52.83, 52.84, 52.85, 52.86,
      52.870000000000005, 52.88, 52.89, 52.9, 52.910000000000004, 52.92, 52.93,
      52.94, 52.95, 52.96, 52.97, 52.980000000000004, 52.99, 53.0, 53.01, 53.02,
      53.03, 53.04, 53.050000000000004, 53.06, 53.07, 53.08, 53.09, 53.1, 53.11,
      53.120000000000005, 53.13, 53.14, 53.15, 53.160000000000004, 53.17, 53.18,
      53.19, 53.2, 53.21, 53.22, 53.230000000000004, 53.24, 53.25, 53.26, 53.27,
      53.28, 53.29, 53.300000000000004, 53.31, 53.32, 53.33, 53.34, 53.35, 53.36,
      53.370000000000005, 53.38, 53.39, 53.4, 53.410000000000004, 53.42, 53.43,
      53.44, 53.45, 53.46, 53.47, 53.480000000000004, 53.49, 53.5, 53.51, 53.52,
      53.53, 53.54, 53.550000000000004, 53.56, 53.57, 53.58, 53.59, 53.6, 53.61,
      53.620000000000005, 53.63, 53.64, 53.65, 53.660000000000004, 53.67, 53.68,
      53.69, 53.7, 53.71, 53.72, 53.730000000000004, 53.74, 53.75, 53.76, 53.77,
      53.78, 53.79, 53.800000000000004, 53.81, 53.82, 53.83, 53.84, 53.85, 53.86,
      53.870000000000005, 53.88, 53.89, 53.9, 53.910000000000004, 53.92, 53.93,
      53.94, 53.95, 53.96, 53.97, 53.980000000000004, 53.99, 54.0, 54.01, 54.02,
      54.03, 54.04, 54.050000000000004, 54.06, 54.07, 54.08, 54.09, 54.1, 54.11,
      54.120000000000005, 54.13, 54.14, 54.15, 54.160000000000004, 54.17, 54.18,
      54.19, 54.2, 54.21, 54.22, 54.230000000000004, 54.24, 54.25, 54.26, 54.27,
      54.28, 54.29, 54.300000000000004, 54.31, 54.32, 54.33, 54.34, 54.35, 54.36,
      54.370000000000005, 54.38, 54.39, 54.4, 54.410000000000004, 54.42, 54.43,
      54.44, 54.45, 54.46, 54.47, 54.480000000000004, 54.49, 54.5, 54.51, 54.52,
      54.53, 54.54, 54.550000000000004, 54.56, 54.57, 54.58, 54.59, 54.6, 54.61,
      54.620000000000005, 54.63, 54.64, 54.65, 54.660000000000004, 54.67, 54.68,
      54.69, 54.7, 54.71, 54.72, 54.730000000000004, 54.74, 54.75, 54.76, 54.77,
      54.78, 54.79, 54.800000000000004, 54.81, 54.82, 54.83, 54.84, 54.85, 54.86,
      54.870000000000005, 54.88, 54.89, 54.9, 54.910000000000004, 54.92, 54.93,
      54.94, 54.95, 54.96, 54.97, 54.980000000000004, 54.99, 55.0, 55.01, 55.02,
      55.03, 55.04, 55.050000000000004, 55.06, 55.07, 55.08, 55.09, 55.1, 55.11,
      55.120000000000005, 55.13, 55.14, 55.15, 55.160000000000004, 55.17, 55.18,
      55.19, 55.2, 55.21, 55.22, 55.230000000000004, 55.24, 55.25, 55.26, 55.27,
      55.28, 55.29, 55.300000000000004, 55.31, 55.32, 55.33, 55.34, 55.35, 55.36,
      55.370000000000005, 55.38, 55.39, 55.4, 55.410000000000004, 55.42, 55.43,
      55.44, 55.45, 55.46, 55.47, 55.480000000000004, 55.49, 55.5, 55.51, 55.52,
      55.53, 55.54, 55.550000000000004, 55.56, 55.57, 55.58, 55.59, 55.6, 55.61,
      55.620000000000005, 55.63, 55.64, 55.65, 55.660000000000004, 55.67, 55.68,
      55.69, 55.7, 55.71, 55.72, 55.730000000000004, 55.74, 55.75, 55.76, 55.77,
      55.78, 55.79, 55.800000000000004, 55.81, 55.82, 55.83, 55.84, 55.85, 55.86,
      55.870000000000005, 55.88, 55.89, 55.9, 55.910000000000004, 55.92, 55.93,
      55.94, 55.95, 55.96, 55.97, 55.980000000000004, 55.99, 56.0, 56.01, 56.02,
      56.03, 56.04, 56.050000000000004, 56.06, 56.07, 56.08, 56.09, 56.1, 56.11,
      56.120000000000005, 56.13, 56.14, 56.15, 56.160000000000004, 56.17, 56.18,
      56.19, 56.2, 56.21, 56.22, 56.230000000000004, 56.24, 56.25, 56.26, 56.27,
      56.28, 56.29, 56.300000000000004, 56.31, 56.32, 56.33, 56.34, 56.35, 56.36,
      56.370000000000005, 56.38, 56.39, 56.4, 56.410000000000004, 56.42, 56.43,
      56.44, 56.45, 56.46, 56.47, 56.480000000000004, 56.49, 56.5, 56.51, 56.52,
      56.53, 56.54, 56.550000000000004, 56.56, 56.57, 56.58, 56.59, 56.6, 56.61,
      56.620000000000005, 56.63, 56.64, 56.65, 56.660000000000004, 56.67, 56.68,
      56.69, 56.7, 56.71, 56.72, 56.730000000000004, 56.74, 56.75, 56.76, 56.77,
      56.78, 56.79, 56.800000000000004, 56.81, 56.82, 56.83, 56.84, 56.85, 56.86,
      56.870000000000005, 56.88, 56.89, 56.9, 56.910000000000004, 56.92, 56.93,
      56.94, 56.95, 56.96, 56.97, 56.980000000000004, 56.99, 57.0, 57.01, 57.02,
      57.03, 57.04, 57.050000000000004, 57.06, 57.07, 57.08, 57.09, 57.1, 57.11,
      57.120000000000005, 57.13, 57.14, 57.15, 57.160000000000004, 57.17, 57.18,
      57.19, 57.2, 57.21, 57.22, 57.230000000000004, 57.24, 57.25, 57.26, 57.27,
      57.28, 57.29, 57.300000000000004, 57.31, 57.32, 57.33, 57.34, 57.35, 57.36,
      57.370000000000005, 57.38, 57.39, 57.4, 57.410000000000004, 57.42, 57.43,
      57.44, 57.45, 57.46, 57.47, 57.480000000000004, 57.49, 57.5, 57.51, 57.52,
      57.53, 57.54, 57.550000000000004, 57.56, 57.57, 57.58, 57.59, 57.6, 57.61,
      57.620000000000005, 57.63, 57.64, 57.65, 57.660000000000004, 57.67, 57.68,
      57.69, 57.7, 57.71, 57.72, 57.730000000000004, 57.74, 57.75, 57.76, 57.77,
      57.78, 57.79, 57.800000000000004, 57.81, 57.82, 57.83, 57.84, 57.85, 57.86,
      57.870000000000005, 57.88, 57.89, 57.9, 57.910000000000004, 57.92, 57.93,
      57.94, 57.95, 57.96, 57.97, 57.980000000000004, 57.99, 58.0, 58.01, 58.02,
      58.03, 58.04, 58.050000000000004, 58.06, 58.07, 58.08, 58.09, 58.1, 58.11,
      58.120000000000005, 58.13, 58.14, 58.15, 58.160000000000004, 58.17, 58.18,
      58.19, 58.2, 58.21, 58.22, 58.230000000000004, 58.24, 58.25, 58.26, 58.27,
      58.28, 58.29, 58.300000000000004, 58.31, 58.32, 58.33, 58.34, 58.35, 58.36,
      58.370000000000005, 58.38, 58.39, 58.4, 58.410000000000004, 58.42, 58.43,
      58.44, 58.45, 58.46, 58.47, 58.480000000000004, 58.49, 58.5, 58.51, 58.52,
      58.53, 58.54, 58.550000000000004, 58.56, 58.57, 58.58, 58.59, 58.6, 58.61,
      58.620000000000005, 58.63, 58.64, 58.65, 58.660000000000004, 58.67, 58.68,
      58.69, 58.7, 58.71, 58.72, 58.730000000000004, 58.74, 58.75, 58.76, 58.77,
      58.78, 58.79, 58.800000000000004, 58.81, 58.82, 58.83, 58.84, 58.85, 58.86,
      58.870000000000005, 58.88, 58.89, 58.9, 58.910000000000004, 58.92, 58.93,
      58.94, 58.95, 58.96, 58.97, 58.980000000000004, 58.99, 59.0, 59.01, 59.02,
      59.03, 59.04, 59.050000000000004, 59.06, 59.07, 59.08, 59.09, 59.1, 59.11,
      59.120000000000005, 59.13, 59.14, 59.15, 59.160000000000004, 59.17, 59.18,
      59.19, 59.2, 59.21, 59.22, 59.230000000000004, 59.24, 59.25, 59.26, 59.27,
      59.28, 59.29, 59.300000000000004, 59.31, 59.32, 59.33, 59.34, 59.35, 59.36,
      59.370000000000005, 59.38, 59.39, 59.4, 59.410000000000004, 59.42, 59.43,
      59.44, 59.45, 59.46, 59.47, 59.480000000000004, 59.49, 59.5, 59.51, 59.52,
      59.53, 59.54, 59.550000000000004, 59.56, 59.57, 59.58, 59.59, 59.6, 59.61,
      59.620000000000005, 59.63, 59.64, 59.65, 59.660000000000004, 59.67, 59.68,
      59.69, 59.7, 59.71, 59.72, 59.730000000000004, 59.74, 59.75, 59.76, 59.77,
      59.78, 59.79, 59.800000000000004, 59.81, 59.82, 59.83, 59.84, 59.85, 59.86,
      59.870000000000005, 59.88, 59.89, 59.9, 59.910000000000004, 59.92, 59.93,
      59.94, 59.95, 59.96, 59.97, 59.980000000000004, 59.99, 60.0, 60.01, 60.02,
      60.03, 60.04, 60.050000000000004, 60.06, 60.07, 60.08, 60.09, 60.1, 60.11,
      60.120000000000005, 60.13, 60.14, 60.15, 60.160000000000004, 60.17, 60.18,
      60.19, 60.2, 60.21, 60.22, 60.230000000000004, 60.24, 60.25, 60.26, 60.27,
      60.28, 60.29, 60.300000000000004, 60.31, 60.32, 60.33, 60.34, 60.35, 60.36,
      60.370000000000005, 60.38, 60.39, 60.4, 60.410000000000004, 60.42, 60.43,
      60.44, 60.45, 60.46, 60.47, 60.480000000000004, 60.49, 60.5, 60.51, 60.52,
      60.53, 60.54, 60.550000000000004, 60.56, 60.57, 60.58, 60.59, 60.6, 60.61,
      60.620000000000005, 60.63, 60.64, 60.65, 60.660000000000004, 60.67, 60.68,
      60.69, 60.7, 60.71, 60.72, 60.730000000000004, 60.74, 60.75, 60.76, 60.77,
      60.78, 60.79, 60.800000000000004, 60.81, 60.82, 60.83, 60.84, 60.85, 60.86,
      60.870000000000005, 60.88, 60.89, 60.9, 60.910000000000004, 60.92, 60.93,
      60.94, 60.95, 60.96, 60.97, 60.980000000000004, 60.99, 61.0, 61.01, 61.02,
      61.03, 61.04, 61.050000000000004, 61.06, 61.07, 61.08, 61.09, 61.1, 61.11,
      61.120000000000005, 61.13, 61.14, 61.15, 61.160000000000004, 61.17, 61.18,
      61.19, 61.2, 61.21, 61.22, 61.230000000000004, 61.24, 61.25, 61.26, 61.27,
      61.28, 61.29, 61.300000000000004, 61.31, 61.32, 61.33, 61.34, 61.35, 61.36,
      61.370000000000005, 61.38, 61.39, 61.4, 61.410000000000004, 61.42, 61.43,
      61.44, 61.45, 61.46, 61.47, 61.480000000000004, 61.49, 61.5, 61.51, 61.52,
      61.53, 61.54, 61.550000000000004, 61.56, 61.57, 61.58, 61.59, 61.6, 61.61,
      61.620000000000005, 61.63, 61.64, 61.65, 61.660000000000004, 61.67, 61.68,
      61.690000000000005, 61.7, 61.71, 61.72, 61.730000000000004, 61.74, 61.75,
      61.76, 61.77, 61.78, 61.79, 61.800000000000004, 61.81, 61.82, 61.83, 61.84,
      61.85, 61.86, 61.870000000000005, 61.88, 61.89, 61.9, 61.910000000000004,
      61.92, 61.93, 61.940000000000005, 61.95, 61.96, 61.97, 61.980000000000004,
      61.99, 62.0, 62.01, 62.02, 62.03, 62.04, 62.050000000000004, 62.06, 62.07,
      62.08, 62.09, 62.1, 62.11, 62.120000000000005, 62.13, 62.14, 62.15,
      62.160000000000004, 62.17, 62.18, 62.190000000000005, 62.2, 62.21, 62.22,
      62.230000000000004, 62.24, 62.25, 62.26, 62.27, 62.28, 62.29,
      62.300000000000004, 62.31, 62.32, 62.33, 62.34, 62.35, 62.36,
      62.370000000000005, 62.38, 62.39, 62.4, 62.410000000000004, 62.42, 62.43,
      62.440000000000005, 62.45, 62.46, 62.47, 62.480000000000004, 62.49, 62.5,
      62.51, 62.52, 62.53, 62.54, 62.550000000000004, 62.56, 62.57, 62.58, 62.59,
      62.6, 62.61, 62.620000000000005, 62.63, 62.64, 62.65, 62.660000000000004,
      62.67, 62.68, 62.690000000000005, 62.7, 62.71, 62.72, 62.730000000000004,
      62.74, 62.75, 62.76, 62.77, 62.78, 62.79, 62.800000000000004, 62.81, 62.82,
      62.83, 62.84, 62.85, 62.86, 62.870000000000005, 62.88, 62.89, 62.9,
      62.910000000000004, 62.92, 62.93, 62.940000000000005, 62.95, 62.96, 62.97,
      62.980000000000004, 62.99, 63.0, 63.01, 63.02, 63.03, 63.04,
      63.050000000000004, 63.06, 63.07, 63.08, 63.09, 63.1, 63.11,
      63.120000000000005, 63.13, 63.14, 63.15, 63.160000000000004, 63.17, 63.18,
      63.190000000000005, 63.2, 63.21, 63.22, 63.230000000000004, 63.24, 63.25,
      63.26, 63.27, 63.28, 63.29, 63.300000000000004, 63.31, 63.32, 63.33, 63.34,
      63.35, 63.36, 63.370000000000005, 63.38, 63.39, 63.4, 63.410000000000004,
      63.42, 63.43, 63.440000000000005, 63.45, 63.46, 63.47, 63.480000000000004,
      63.49, 63.5, 63.51, 63.52, 63.53, 63.54, 63.550000000000004, 63.56, 63.57,
      63.58, 63.59, 63.6, 63.61, 63.620000000000005, 63.63, 63.64, 63.65,
      63.660000000000004, 63.67, 63.68, 63.690000000000005, 63.7, 63.71, 63.72,
      63.730000000000004, 63.74, 63.75, 63.76, 63.77, 63.78, 63.79,
      63.800000000000004, 63.81, 63.82, 63.83, 63.84, 63.85, 63.86,
      63.870000000000005, 63.88, 63.89, 63.9, 63.910000000000004, 63.92, 63.93,
      63.940000000000005, 63.95, 63.96, 63.97, 63.980000000000004, 63.99, 64.0,
      64.01, 64.02, 64.03, 64.04, 64.05, 64.06, 64.070000000000007, 64.08, 64.09,
      64.1, 64.11, 64.12, 64.13, 64.14, 64.15, 64.16, 64.17, 64.18, 64.19, 64.2,
      64.210000000000008, 64.22, 64.23, 64.24, 64.25, 64.26, 64.27, 64.28, 64.29,
      64.3, 64.31, 64.320000000000007, 64.33, 64.34, 64.35, 64.36, 64.37, 64.38,
      64.39, 64.4, 64.41, 64.42, 64.43, 64.44, 64.45, 64.460000000000008, 64.47,
      64.48, 64.49, 64.5, 64.51, 64.52, 64.53, 64.54, 64.55, 64.56,
      64.570000000000007, 64.58, 64.59, 64.6, 64.61, 64.62, 64.63, 64.64, 64.65,
      64.66, 64.67, 64.68, 64.69, 64.7, 64.710000000000008, 64.72, 64.73, 64.74,
      64.75, 64.76, 64.77, 64.78, 64.79, 64.8, 64.81, 64.820000000000007, 64.83,
      64.84, 64.85, 64.86, 64.87, 64.88, 64.89, 64.9, 64.91, 64.92, 64.93, 64.94,
      64.95, 64.960000000000008, 64.97, 64.98, 64.99, 65.0, 65.01, 65.02, 65.03,
      65.04, 65.05, 65.06, 65.070000000000007, 65.08, 65.09, 65.1, 65.11, 65.12,
      65.13, 65.14, 65.15, 65.16, 65.17, 65.18, 65.19, 65.2, 65.210000000000008,
      65.22, 65.23, 65.24, 65.25, 65.26, 65.27, 65.28, 65.29, 65.3, 65.31,
      65.320000000000007, 65.33, 65.34, 65.35, 65.36, 65.37, 65.38, 65.39, 65.4,
      65.41, 65.42, 65.43, 65.44, 65.45, 65.460000000000008, 65.47, 65.48, 65.49,
      65.5, 65.51, 65.52, 65.53, 65.54, 65.55, 65.56, 65.570000000000007, 65.58,
      65.59, 65.6, 65.61, 65.62, 65.63, 65.64, 65.65, 65.66, 65.67, 65.68, 65.69,
      65.7, 65.710000000000008, 65.72, 65.73, 65.74, 65.75, 65.76, 65.77, 65.78,
      65.79, 65.8, 65.81, 65.820000000000007, 65.83, 65.84, 65.85, 65.86, 65.87,
      65.88, 65.89, 65.9, 65.91, 65.92, 65.93, 65.94, 65.95, 65.960000000000008,
      65.97, 65.98, 65.99, 66.0, 66.01, 66.02, 66.03, 66.04, 66.05, 66.06,
      66.070000000000007, 66.08, 66.09, 66.1, 66.11, 66.12, 66.13, 66.14, 66.15,
      66.16, 66.17, 66.18, 66.19, 66.2, 66.210000000000008, 66.22, 66.23, 66.24,
      66.25, 66.26, 66.27, 66.28, 66.29, 66.3, 66.31, 66.320000000000007, 66.33,
      66.34, 66.35, 66.36, 66.37, 66.38, 66.39, 66.4, 66.41, 66.42, 66.43, 66.44,
      66.45, 66.460000000000008, 66.47, 66.48, 66.49, 66.5, 66.51, 66.52, 66.53,
      66.54, 66.55, 66.56, 66.570000000000007, 66.58, 66.59, 66.6, 66.61, 66.62,
      66.63, 66.64, 66.65, 66.66, 66.67, 66.68, 66.69, 66.7, 66.710000000000008,
      66.72, 66.73, 66.74, 66.75, 66.76, 66.77, 66.78, 66.79, 66.8, 66.81,
      66.820000000000007, 66.83, 66.84, 66.85, 66.86, 66.87, 66.88, 66.89, 66.9,
      66.91, 66.92, 66.93, 66.94, 66.95, 66.960000000000008, 66.97, 66.98, 66.99,
      67.0, 67.01, 67.02, 67.03, 67.04, 67.05, 67.06, 67.070000000000007, 67.08,
      67.09, 67.1, 67.11, 67.12, 67.13, 67.14, 67.15, 67.16, 67.17, 67.18, 67.19,
      67.2, 67.210000000000008, 67.22, 67.23, 67.24, 67.25, 67.26, 67.27, 67.28,
      67.29, 67.3, 67.31, 67.320000000000007, 67.33, 67.34, 67.35, 67.36, 67.37,
      67.38, 67.39, 67.4, 67.41, 67.42, 67.43, 67.44, 67.45, 67.460000000000008,
      67.47, 67.48, 67.49, 67.5, 67.51, 67.52, 67.53, 67.54, 67.55, 67.56,
      67.570000000000007, 67.58, 67.59, 67.6, 67.61, 67.62, 67.63, 67.64, 67.65,
      67.66, 67.67, 67.68, 67.69, 67.7, 67.710000000000008, 67.72, 67.73, 67.74,
      67.75, 67.76, 67.77, 67.78, 67.79, 67.8, 67.81, 67.820000000000007, 67.83,
      67.84, 67.85, 67.86, 67.87, 67.88, 67.89, 67.9, 67.91, 67.92, 67.93, 67.94,
      67.95, 67.960000000000008, 67.97, 67.98, 67.99, 68.0, 68.01, 68.02, 68.03,
      68.04, 68.05, 68.06, 68.070000000000007, 68.08, 68.09, 68.1, 68.11, 68.12,
      68.13, 68.14, 68.15, 68.16, 68.17, 68.18, 68.19, 68.2, 68.210000000000008,
      68.22, 68.23, 68.24, 68.25, 68.26, 68.27, 68.28, 68.29, 68.3, 68.31,
      68.320000000000007, 68.33, 68.34, 68.350000000000009, 68.36, 68.37, 68.38,
      68.39, 68.4, 68.41, 68.42, 68.43, 68.44, 68.45, 68.460000000000008, 68.47,
      68.48, 68.49, 68.5, 68.51, 68.52, 68.53, 68.54, 68.55, 68.56,
      68.570000000000007, 68.58, 68.59, 68.600000000000009, 68.61, 68.62, 68.63,
      68.64, 68.65, 68.66, 68.67, 68.68, 68.69, 68.7, 68.710000000000008, 68.72,
      68.73, 68.74, 68.75, 68.76, 68.77, 68.78, 68.79, 68.8, 68.81,
      68.820000000000007, 68.83, 68.84, 68.850000000000009, 68.86, 68.87, 68.88,
      68.89, 68.9, 68.91, 68.92, 68.93, 68.94, 68.95, 68.960000000000008, 68.97,
      68.98, 68.99, 69.0, 69.01, 69.02, 69.03, 69.04, 69.05, 69.06,
      69.070000000000007, 69.08, 69.09, 69.100000000000009, 69.11, 69.12, 69.13,
      69.14, 69.15, 69.16, 69.17, 69.18, 69.19, 69.2, 69.210000000000008, 69.22,
      69.23, 69.24, 69.25, 69.26, 69.27, 69.28, 69.29, 69.3, 69.31,
      69.320000000000007, 69.33, 69.34, 69.350000000000009, 69.36, 69.37, 69.38,
      69.39, 69.4, 69.41, 69.42, 69.43, 69.44, 69.45, 69.460000000000008, 69.47,
      69.48, 69.49, 69.5, 69.51, 69.52, 69.53, 69.54, 69.55, 69.56,
      69.570000000000007, 69.58, 69.59, 69.600000000000009, 69.61, 69.62, 69.63,
      69.64, 69.65, 69.66, 69.67, 69.68, 69.69, 69.7, 69.710000000000008, 69.72,
      69.73, 69.74, 69.75, 69.76, 69.77, 69.78, 69.79, 69.8, 69.81,
      69.820000000000007, 69.83, 69.84, 69.850000000000009, 69.86, 69.87, 69.88,
      69.89, 69.9, 69.91, 69.92, 69.93, 69.94, 69.95, 69.960000000000008, 69.97,
      69.98, 69.99, 70.0, 70.01, 70.02, 70.03, 70.04, 70.05, 70.06,
      70.070000000000007, 70.08, 70.09, 70.100000000000009, 70.11, 70.12, 70.13,
      70.14, 70.15, 70.16, 70.17, 70.18, 70.19, 70.2, 70.210000000000008, 70.22,
      70.23, 70.24, 70.25, 70.26, 70.27, 70.28, 70.29, 70.3, 70.31,
      70.320000000000007, 70.33, 70.34, 70.350000000000009, 70.36, 70.37, 70.38,
      70.39, 70.4, 70.41, 70.42, 70.43, 70.44, 70.45, 70.460000000000008, 70.47,
      70.48, 70.49, 70.5, 70.51, 70.52, 70.53, 70.54, 70.55, 70.56,
      70.570000000000007, 70.58, 70.59, 70.600000000000009, 70.61, 70.62, 70.63,
      70.64, 70.65, 70.66, 70.67, 70.68, 70.69, 70.7, 70.710000000000008, 70.72,
      70.73, 70.74, 70.75, 70.76, 70.77, 70.78, 70.79, 70.8, 70.81,
      70.820000000000007, 70.83, 70.84, 70.850000000000009, 70.86, 70.87, 70.88,
      70.89, 70.9, 70.91, 70.92, 70.93, 70.94, 70.95, 70.960000000000008, 70.97,
      70.98, 70.99, 71.0, 71.01, 71.02, 71.03, 71.04, 71.05, 71.06,
      71.070000000000007, 71.08, 71.09, 71.100000000000009, 71.11, 71.12, 71.13,
      71.14, 71.15, 71.16, 71.17, 71.18, 71.19, 71.2, 71.210000000000008, 71.22,
      71.23, 71.24, 71.25, 71.26, 71.27, 71.28, 71.29, 71.3, 71.31,
      71.320000000000007, 71.33, 71.34, 71.350000000000009, 71.36, 71.37, 71.38,
      71.39, 71.4, 71.41, 71.42, 71.43, 71.44, 71.45, 71.460000000000008, 71.47,
      71.48, 71.49, 71.5, 71.51, 71.52, 71.53, 71.54, 71.55, 71.56,
      71.570000000000007, 71.58, 71.59, 71.600000000000009, 71.61, 71.62, 71.63,
      71.64, 71.65, 71.66, 71.67, 71.68, 71.69, 71.7, 71.710000000000008, 71.72,
      71.73, 71.74, 71.75, 71.76, 71.77, 71.78, 71.79, 71.8, 71.81,
      71.820000000000007, 71.83, 71.84, 71.850000000000009, 71.86, 71.87, 71.88,
      71.89, 71.9, 71.91, 71.92, 71.93, 71.94, 71.95, 71.960000000000008, 71.97,
      71.98, 71.99, 72.0, 72.01, 72.02, 72.03, 72.04, 72.05, 72.06,
      72.070000000000007, 72.08, 72.09, 72.100000000000009, 72.11, 72.12, 72.13,
      72.14, 72.15, 72.16, 72.17, 72.18, 72.19, 72.2, 72.210000000000008, 72.22,
      72.23, 72.24, 72.25, 72.26, 72.27, 72.28, 72.29, 72.3, 72.31,
      72.320000000000007, 72.33, 72.34, 72.350000000000009, 72.36, 72.37, 72.38,
      72.39, 72.4, 72.41, 72.42, 72.43, 72.44, 72.45, 72.460000000000008, 72.47,
      72.48, 72.49, 72.5, 72.51, 72.52, 72.53, 72.54, 72.55, 72.56,
      72.570000000000007, 72.58, 72.59, 72.600000000000009, 72.61, 72.62, 72.63,
      72.64, 72.65, 72.66, 72.67, 72.68, 72.69, 72.7, 72.710000000000008, 72.72,
      72.73, 72.74, 72.75, 72.76, 72.77, 72.78, 72.79, 72.8, 72.81,
      72.820000000000007, 72.83, 72.84, 72.850000000000009, 72.86, 72.87, 72.88,
      72.89, 72.9, 72.91, 72.92, 72.93, 72.94, 72.95, 72.960000000000008, 72.97,
      72.98, 72.99, 73.0, 73.01, 73.02, 73.03, 73.04, 73.05, 73.06,
      73.070000000000007, 73.08, 73.09, 73.100000000000009, 73.11, 73.12, 73.13,
      73.14, 73.15, 73.16, 73.17, 73.18, 73.19, 73.2, 73.210000000000008, 73.22,
      73.23, 73.24, 73.25, 73.26, 73.27, 73.28, 73.29, 73.3, 73.31,
      73.320000000000007, 73.33, 73.34, 73.350000000000009, 73.36, 73.37, 73.38,
      73.39, 73.4, 73.41, 73.42, 73.43, 73.44, 73.45, 73.460000000000008, 73.47,
      73.48, 73.49, 73.5, 73.51, 73.52, 73.53, 73.54, 73.55, 73.56,
      73.570000000000007, 73.58, 73.59, 73.600000000000009, 73.61, 73.62, 73.63,
      73.64, 73.65, 73.66, 73.67, 73.68, 73.69, 73.7, 73.710000000000008, 73.72,
      73.73, 73.74, 73.75, 73.76, 73.77, 73.78, 73.79, 73.8, 73.81,
      73.820000000000007, 73.83, 73.84, 73.850000000000009, 73.86, 73.87, 73.88,
      73.89, 73.9, 73.91, 73.92, 73.93, 73.94, 73.95, 73.960000000000008, 73.97,
      73.98, 73.99, 74.0, 74.01, 74.02, 74.03, 74.04, 74.05, 74.06,
      74.070000000000007, 74.08, 74.09, 74.100000000000009, 74.11, 74.12, 74.13,
      74.14, 74.15, 74.16, 74.17, 74.18, 74.19, 74.2, 74.210000000000008, 74.22,
      74.23, 74.24, 74.25, 74.26, 74.27, 74.28, 74.29, 74.3, 74.31,
      74.320000000000007, 74.33, 74.34, 74.350000000000009, 74.36, 74.37, 74.38,
      74.39, 74.4, 74.41, 74.42, 74.43, 74.44, 74.45, 74.460000000000008, 74.47,
      74.48, 74.49, 74.5, 74.51, 74.52, 74.53, 74.54, 74.55, 74.56,
      74.570000000000007, 74.58, 74.59, 74.600000000000009, 74.61, 74.62, 74.63,
      74.64, 74.65, 74.66, 74.67, 74.68, 74.69, 74.7, 74.710000000000008, 74.72,
      74.73, 74.74, 74.75, 74.76, 74.77, 74.78, 74.79, 74.8, 74.81,
      74.820000000000007, 74.83, 74.84, 74.850000000000009, 74.86, 74.87, 74.88,
      74.89, 74.9, 74.91, 74.92, 74.93, 74.94, 74.95, 74.960000000000008, 74.97,
      74.98, 74.99, 75.0, 75.01, 75.02, 75.03, 75.04, 75.05, 75.06,
      75.070000000000007, 75.08, 75.09, 75.100000000000009, 75.11, 75.12, 75.13,
      75.14, 75.15, 75.16, 75.17, 75.18, 75.19, 75.2, 75.210000000000008, 75.22,
      75.23, 75.24, 75.25, 75.26, 75.27, 75.28, 75.29, 75.3, 75.31,
      75.320000000000007, 75.33, 75.34, 75.350000000000009, 75.36, 75.37, 75.38,
      75.39, 75.4, 75.41, 75.42, 75.43, 75.44, 75.45, 75.460000000000008, 75.47,
      75.48, 75.49, 75.5, 75.51, 75.52, 75.53, 75.54, 75.55, 75.56,
      75.570000000000007, 75.58, 75.59, 75.600000000000009, 75.61, 75.62, 75.63,
      75.64, 75.65, 75.66, 75.67, 75.68, 75.69, 75.7, 75.710000000000008, 75.72,
      75.73, 75.74, 75.75, 75.76, 75.77, 75.78, 75.79, 75.8, 75.81,
      75.820000000000007, 75.83, 75.84, 75.850000000000009, 75.86, 75.87, 75.88,
      75.89, 75.9, 75.91, 75.92, 75.93, 75.94, 75.95, 75.960000000000008, 75.97,
      75.98, 75.99, 76.0, 76.01, 76.02, 76.03, 76.04, 76.05, 76.06,
      76.070000000000007, 76.08, 76.09, 76.100000000000009, 76.11, 76.12, 76.13,
      76.14, 76.15, 76.16, 76.17, 76.18, 76.19, 76.2, 76.210000000000008, 76.22,
      76.23, 76.24, 76.25, 76.26, 76.27, 76.28, 76.29, 76.3, 76.31,
      76.320000000000007, 76.33, 76.34, 76.350000000000009, 76.36, 76.37, 76.38,
      76.39, 76.4, 76.41, 76.42, 76.43, 76.44, 76.45, 76.460000000000008, 76.47,
      76.48, 76.49, 76.5, 76.51, 76.52, 76.53, 76.54, 76.55, 76.56,
      76.570000000000007, 76.58, 76.59, 76.600000000000009, 76.61, 76.62, 76.63,
      76.64, 76.65, 76.66, 76.67, 76.68, 76.69, 76.7, 76.710000000000008, 76.72,
      76.73, 76.74, 76.75, 76.76, 76.77, 76.78, 76.79, 76.8, 76.81,
      76.820000000000007, 76.83, 76.84, 76.850000000000009, 76.86, 76.87, 76.88,
      76.89, 76.9, 76.91, 76.92, 76.93, 76.94, 76.95, 76.960000000000008, 76.97,
      76.98, 76.99, 77.0, 77.01, 77.02, 77.03, 77.04, 77.05, 77.06,
      77.070000000000007, 77.08, 77.09, 77.100000000000009, 77.11, 77.12, 77.13,
      77.14, 77.15, 77.16, 77.17, 77.18, 77.19, 77.2, 77.210000000000008, 77.22,
      77.23, 77.24, 77.25, 77.26, 77.27, 77.28, 77.29, 77.3, 77.31,
      77.320000000000007, 77.33, 77.34, 77.350000000000009, 77.36, 77.37, 77.38,
      77.39, 77.4, 77.41, 77.42, 77.43, 77.44, 77.45, 77.460000000000008, 77.47,
      77.48, 77.49, 77.5, 77.51, 77.52, 77.53, 77.54, 77.55, 77.56,
      77.570000000000007, 77.58, 77.59, 77.600000000000009, 77.61, 77.62, 77.63,
      77.64, 77.65, 77.66, 77.67, 77.68, 77.69, 77.7, 77.710000000000008, 77.72,
      77.73, 77.74, 77.75, 77.76, 77.77, 77.78, 77.79, 77.8, 77.81,
      77.820000000000007, 77.83, 77.84, 77.850000000000009, 77.86, 77.87, 77.88,
      77.89, 77.9, 77.91, 77.92, 77.93, 77.94, 77.95, 77.960000000000008, 77.97,
      77.98, 77.99, 78.0, 78.01, 78.02, 78.03, 78.04, 78.05, 78.06,
      78.070000000000007, 78.08, 78.09, 78.100000000000009, 78.11, 78.12, 78.13,
      78.14, 78.15, 78.16, 78.17, 78.18, 78.19, 78.2, 78.210000000000008, 78.22,
      78.23, 78.24, 78.25, 78.26, 78.27, 78.28, 78.29, 78.3, 78.31,
      78.320000000000007, 78.33, 78.34, 78.350000000000009, 78.36, 78.37, 78.38,
      78.39, 78.4, 78.41, 78.42, 78.43, 78.44, 78.45, 78.460000000000008, 78.47,
      78.48, 78.49, 78.5, 78.51, 78.52, 78.53, 78.54, 78.55, 78.56,
      78.570000000000007, 78.58, 78.59, 78.600000000000009, 78.61, 78.62, 78.63,
      78.64, 78.65, 78.66, 78.67, 78.68, 78.69, 78.7, 78.710000000000008, 78.72,
      78.73, 78.74, 78.75, 78.76, 78.77, 78.78, 78.79, 78.8, 78.81,
      78.820000000000007, 78.83, 78.84, 78.850000000000009, 78.86, 78.87, 78.88,
      78.89, 78.9, 78.91, 78.92, 78.93, 78.94, 78.95, 78.960000000000008, 78.97,
      78.98, 78.99, 79.0, 79.01, 79.02, 79.03, 79.04, 79.05, 79.06,
      79.070000000000007, 79.08, 79.09, 79.100000000000009, 79.11, 79.12, 79.13,
      79.14, 79.15, 79.16, 79.17, 79.18, 79.19, 79.2, 79.210000000000008, 79.22,
      79.23, 79.24, 79.25, 79.26, 79.27, 79.28, 79.29, 79.3, 79.31,
      79.320000000000007, 79.33, 79.34, 79.350000000000009, 79.36, 79.37, 79.38,
      79.39, 79.4, 79.41, 79.42, 79.43, 79.44, 79.45, 79.460000000000008, 79.47,
      79.48, 79.49, 79.5, 79.51, 79.52, 79.53, 79.54, 79.55, 79.56,
      79.570000000000007, 79.58, 79.59, 79.600000000000009, 79.61, 79.62, 79.63,
      79.64, 79.65, 79.66, 79.67, 79.68, 79.69, 79.7, 79.710000000000008, 79.72,
      79.73, 79.74, 79.75, 79.76, 79.77, 79.78, 79.79, 79.8, 79.81,
      79.820000000000007, 79.83, 79.84, 79.850000000000009, 79.86, 79.87, 79.88,
      79.89, 79.9, 79.91, 79.92, 79.93, 79.94, 79.95, 79.960000000000008, 79.97,
      79.98, 79.99 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0,
      250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0 } ;

    ModelWithControllersOnly_DW.FromWorkspace4_PWORK.TimePtr = static_cast<void *>
      (pTimeValues0);
    ModelWithControllersOnly_DW.FromWorkspace4_PWORK.DataPtr = static_cast<void *>
      (pDataValues0);
    ModelWithControllersOnly_DW.FromWorkspace4_IWORK.PrevIndex = 0;
  }

  /* Start for Assertion: '<S12>/Assertion' */
  ModelWithControllersOnly_DW.Assertion_sltestFinalResult =
    slTestResult_Untested;
  slTestInitialize(&ModelWithControllersOnly_DW.Assertion_sltestBlkInfo,
                   &ModelWithControllersOnly_DW.Assertion_sltestCurrentResult,
                   &ModelWithControllersOnly_DW.Assertion_sltestFinalResult,
                   &ModelWithControllersOnly_DW.Assertion_sltestLastResultTime,
                   1);
  slTestRegAssessment(&ModelWithControllersOnly_DW.Assertion_sltestBlkInfo, 0,
                      "", 0, 0, 0,
                      "ModelWithControllersOnly/Longitudinal Controller Stanley/Verify Direction/Assertion",
                      "", 0);

  /* Start for Assertion: '<S12>/Assertion1' */
  ModelWithControllersOnly_DW.Assertion1_sltestFinalResult =
    slTestResult_Untested;
  slTestInitialize(&ModelWithControllersOnly_DW.Assertion1_sltestBlkInfo,
                   &ModelWithControllersOnly_DW.Assertion1_sltestCurrentResult,
                   &ModelWithControllersOnly_DW.Assertion1_sltestFinalResult,
                   &ModelWithControllersOnly_DW.Assertion1_sltestLastResultTime,
                   1);
  slTestRegAssessment(&ModelWithControllersOnly_DW.Assertion1_sltestBlkInfo, 0,
                      "", 0, 0, 0,
                      "ModelWithControllersOnly/Longitudinal Controller Stanley/Verify Direction/Assertion1",
                      "", 0);

  /* Start for SwitchCase: '<S11>/Switch Case' */
  ModelWithControllersOnly_DW.SwitchCase_ActiveSubsystem = -1;
  ModelWithControllersOnl_PrevZCX.Integrator_Reset_ZCE = UNINITIALIZED_ZCSIG;
  ModelWithControllersOnl_PrevZCX.Integrator_Reset_ZCE_c = UNINITIALIZED_ZCSIG;

  /* InitializeConditions for Integrator: '<S214>/Integrator' incorporates:
   *  Integrator: '<S149>/Integrator'
   */
  if (rtmIsFirstInitCond(ModelWithControllersOnly_M)) {
    ModelWithControllersOnly_X.Integrator_CSTATE[0] = 0.0;
    ModelWithControllersOnly_X.Integrator_CSTATE[1] = 0.0;
    ModelWithControllersOnly_X.Integrator_CSTATE[2] = 0.0;
    ModelWithControllersOnly_X.Integrator_CSTATE[3] = 0.0;
    ModelWithControllersOnly_X.Integrator_CSTATE_b[0] = 0.0;
    ModelWithControllersOnly_X.Integrator_CSTATE_b[1] = 0.0;
  }

  ModelWithControllersOnly_DW.Integrator_IWORK = 1;

  /* End of InitializeConditions for Integrator: '<S214>/Integrator' */

  /* InitializeConditions for Integrator: '<S149>/Integrator' */
  ModelWithControllersOnly_DW.Integrator_IWORK_c = 1;

  /* InitializeConditions for UnitDelay: '<S8>/Unit Delay' */
  ModelWithControllersOnly_DW.UnitDelay_DSTATE =
    ModelWithControllersOnly_cal->UnitDelay_InitialCondition;

  /* InitializeConditions for UnitDelay: '<Root>/Unit Delay' */
  ModelWithControllersOnly_DW.UnitDelay_DSTATE_e =
    ModelWithControllersOnly_cal->UnitDelay_InitialCondition_g;

  /* InitializeConditions for RateTransition generated from: '<S12>/Equal1' */
  ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_Buf0 =
    ModelWithControllersOnly_cal->TmpRTBAtEqual1Inport1_InitialCo;

  /* InitializeConditions for RateTransition generated from: '<S12>/Equal2' */
  ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_Buf0 =
    ModelWithControllersOnly_cal->TmpRTBAtEqual2Inport1_InitialCo;

  /* InitializeConditions for RateTransition generated from: '<S12>/Equal4' */
  ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_Buf0 =
    ModelWithControllersOnly_cal->TmpRTBAtEqual4Inport2_InitialCo;

  /* InitializeConditions for RateTransition generated from: '<S12>/Equal5' */
  ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_Buf0 =
    ModelWithControllersOnly_cal->TmpRTBAtEqual5Inport2_InitialCo;

  /* InitializeConditions for RateTransition generated from: '<S12>/Equal3' */
  ModelWithControllersOnly_DW.xdot_Buf[0] =
    ModelWithControllersOnly_cal->xdot_InitialCondition;
  ModelWithControllersOnly_DW.xdot_WrBufIdx = 0;
  ModelWithControllersOnly_DW.xdot_RdBufIdx = 1;

  /* InitializeConditions for RateTransition generated from: '<S12>/Sign1' */
  ModelWithControllersOnly_DW.xdot1_Buf[0] =
    ModelWithControllersOnly_cal->xdot1_InitialCondition;
  ModelWithControllersOnly_DW.xdot1_WrBufIdx = 0;
  ModelWithControllersOnly_DW.xdot1_RdBufIdx = 1;

  /* InitializeConditions for RateTransition generated from: '<S10>/Multiply' */
  ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_Buf0 =
    ModelWithControllersOnly_cal->TmpRTBAtMultiplyInport1_Initial;

  /* InitializeConditions for RateTransition generated from: '<S11>/Switch Case' */
  ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_Buf0 =
    ModelWithControllersOnly_cal->TmpRTBAtSwitchCaseInport1_Initi;

  /* InitializeConditions for RateTransition generated from: '<S2>/Minus' */
  ModelWithControllersOnly_DW.xdot_Buf_d[0] =
    ModelWithControllersOnly_cal->xdot_InitialCondition_b;
  ModelWithControllersOnly_DW.xdot_WrBufIdx_b = 0;
  ModelWithControllersOnly_DW.xdot_RdBufIdx_m = 1;

  /* InitializeConditions for RateTransition generated from: '<S119>/Integrator' */
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_Buf[0] =
    ModelWithControllersOnly_cal->TmpRTBAtIntegratorInport2_Initi;
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_WrBuf = 0;
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_RdBuf = 1;

  /* InitializeConditions for Integrator: '<S119>/Integrator' */
  ModelWithControllersOnly_X.Integrator_CSTATE_o =
    ModelWithControllersOnly_cal->Integrator_IC;

  /* InitializeConditions for RateTransition generated from: '<S119>/Sum' */
  ModelWithControllersOnly_DW.TmpRTBAtSumInport3_Buf[0] =
    ModelWithControllersOnly_cal->TmpRTBAtSumInport3_InitialCondi;
  ModelWithControllersOnly_DW.TmpRTBAtSumInport3_WrBufIdx = 0;
  ModelWithControllersOnly_DW.TmpRTBAtSumInport3_RdBufIdx = 1;

  /* InitializeConditions for RateTransition generated from: '<S118>/Integrator' */
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_Buf_m[0] =
    ModelWithControllersOnly_cal->TmpRTBAtIntegratorInport2_Ini_j;
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_WrB_d = 0;
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport2_RdB_d = 1;

  /* InitializeConditions for Integrator: '<S118>/Integrator' */
  ModelWithControllersOnly_X.Integrator_CSTATE_c =
    ModelWithControllersOnly_cal->Integrator_IC_g;

  /* InitializeConditions for RateTransition generated from: '<S118>/Sum' */
  ModelWithControllersOnly_DW.TmpRTBAtSumInport3_Buf_f[0] =
    ModelWithControllersOnly_cal->TmpRTBAtSumInport3_InitialCon_j;
  ModelWithControllersOnly_DW.TmpRTBAtSumInport3_WrBufIdx_o = 0;
  ModelWithControllersOnly_DW.TmpRTBAtSumInport3_RdBufIdx_l = 1;

  /* InitializeConditions for RateTransition generated from: '<S118>/Integrator' */
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_Buf[0] =
    ModelWithControllersOnly_cal->TmpRTBAtIntegratorInport1_Initi;
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_WrBuf = 0;
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_RdBuf = 1;

  /* InitializeConditions for RateTransition generated from: '<S119>/Integrator' */
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_Buf_p[0] =
    ModelWithControllersOnly_cal->TmpRTBAtIntegratorInport1_Ini_j;
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_WrB_k = 0;
  ModelWithControllersOnly_DW.TmpRTBAtIntegratorInport1_RdB_k = 1;

  /* InitializeConditions for Derivative: '<S120>/Derivative' */
  ModelWithControllersOnly_DW.TimeStampA = (rtInf);
  ModelWithControllersOnly_DW.TimeStampB = (rtInf);

  /* InitializeConditions for Derivative: '<S120>/Derivative1' */
  ModelWithControllersOnly_DW.TimeStampA_l = (rtInf);
  ModelWithControllersOnly_DW.TimeStampB_l = (rtInf);

  /* InitializeConditions for Derivative: '<S120>/Derivative2' */
  ModelWithControllersOnly_DW.TimeStampA_e = (rtInf);
  ModelWithControllersOnly_DW.TimeStampB_k = (rtInf);

  /* InitializeConditions for Derivative: '<S120>/Derivative3' */
  ModelWithControllersOnly_DW.TimeStampA_o = (rtInf);
  ModelWithControllersOnly_DW.TimeStampB_lf = (rtInf);

  /* InitializeConditions for Derivative: '<S120>/Derivative4' */
  ModelWithControllersOnly_DW.TimeStampA_l5 = (rtInf);
  ModelWithControllersOnly_DW.TimeStampB_e = (rtInf);

  /* InitializeConditions for Derivative: '<S120>/Derivative5' */
  ModelWithControllersOnly_DW.TimeStampA_k = (rtInf);
  ModelWithControllersOnly_DW.TimeStampB_l4 = (rtInf);

  /* InitializeConditions for RateTransition: '<S120>/Rate Transition' */
  ModelWithControllersOnly_DW.RateTransition_Buf[0] =
    ModelWithControllersOnly_cal->RateTransition_InitialCondition;
  ModelWithControllersOnly_DW.RateTransition_WrBufIdx = 0;
  ModelWithControllersOnly_DW.RateTransition_RdBufIdx = 1;

  /* InitializeConditions for RateTransition: '<S120>/Rate Transition1' */
  ModelWithControllersOnly_DW.RateTransition1_Buf[0] =
    ModelWithControllersOnly_cal->RateTransition1_InitialConditio;
  ModelWithControllersOnly_DW.RateTransition1_WrBufIdx = 0;
  ModelWithControllersOnly_DW.RateTransition1_RdBufIdx = 1;

  /* InitializeConditions for RateTransition: '<S120>/Rate Transition2' */
  ModelWithControllersOnly_DW.RateTransition2_Buf[0] =
    ModelWithControllersOnly_cal->RateTransition2_InitialConditio;
  ModelWithControllersOnly_DW.RateTransition2_WrBufIdx = 0;
  ModelWithControllersOnly_DW.RateTransition2_RdBufIdx = 1;

  /* InitializeConditions for RateTransition generated from: '<S120>/MATLAB Function1' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport1_[0] =
    ModelWithControllersOnly_cal->TmpRTBAtMATLABFunction1Inport1_;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_j = 0;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_h = 1;

  /* InitializeConditions for RateTransition generated from: '<S120>/MATLAB Function1' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport2_[0] =
    ModelWithControllersOnly_cal->TmpRTBAtMATLABFunction1Inport2_;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_m = 0;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_p = 1;

  /* InitializeConditions for RateTransition generated from: '<S120>/MATLAB Function1' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport3_ =
    ModelWithControllersOnly_cal->TmpRTBAtMATLABFunction1Inport3_;

  /* InitializeConditions for RateTransition generated from: '<S120>/MATLAB Function2' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport1_[0] =
    ModelWithControllersOnly_cal->TmpRTBAtMATLABFunction2Inport1_;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_d = 0;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_e = 1;

  /* InitializeConditions for RateTransition generated from: '<S120>/MATLAB Function2' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport2_[0] =
    ModelWithControllersOnly_cal->TmpRTBAtMATLABFunction2Inport2_;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_k = 0;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_m = 1;

  /* InitializeConditions for RateTransition generated from: '<S120>/MATLAB Function2' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport3_[0] =
    ModelWithControllersOnly_cal->TmpRTBAtMATLABFunction2Inport3_;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_p = 0;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction2Inport_a = 1;

  /* InitializeConditions for RateTransition generated from: '<S120>/MATLAB Function3' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport1_[0] =
    ModelWithControllersOnly_cal->TmpRTBAtMATLABFunction3Inport1_;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport_g = 0;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport_a = 1;

  /* InitializeConditions for RateTransition generated from: '<S120>/MATLAB Function3' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport2_[0] =
    ModelWithControllersOnly_cal->TmpRTBAtMATLABFunction3Inport2_;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport_b = 0;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport_f = 1;

  /* InitializeConditions for RateTransition generated from: '<S120>/MATLAB Function3' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inport3_[0] =
    ModelWithControllersOnly_cal->TmpRTBAtMATLABFunction3Inport3_;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inpor_bx = 0;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction3Inpor_gh = 1;

  /* InitializeConditions for RateTransition generated from: '<S120>/MATLAB Function' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport1_B[0] =
    ModelWithControllersOnly_cal->TmpRTBAtMATLABFunctionInport1_I;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport1_W = 0;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport1_R = 1;

  /* InitializeConditions for RateTransition generated from: '<S120>/MATLAB Function' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport2_B[0] =
    ModelWithControllersOnly_cal->TmpRTBAtMATLABFunctionInport2_I;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport2_W = 0;
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport2_R = 1;

  /* InitializeConditions for RateTransition generated from: '<S120>/MATLAB Function' */
  ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_B =
    ModelWithControllersOnly_cal->TmpRTBAtMATLABFunctionInport3_I;

  /* InitializeConditions for Derivative: '<S121>/Derivative' */
  ModelWithControllersOnly_DW.TimeStampA_h = (rtInf);
  ModelWithControllersOnly_DW.TimeStampB_h = (rtInf);

  /* InitializeConditions for FirstOrderHold: '<S121>/First Order Hold' */
  ModelWithControllersOnly_DW.Tk = (rtInf);
  ModelWithControllersOnly_DW.Ck =
    ModelWithControllersOnly_cal->FirstOrderHold_IniOut;
  ModelWithControllersOnly_DW.Uk = (rtInf);
  ModelWithControllersOnly_DW.Mk = 0.0;

  /* InitializeConditions for FirstOrderHold: '<S121>/First Order Hold1' */
  ModelWithControllersOnly_DW.Tk_d = (rtInf);
  ModelWithControllersOnly_DW.Ck_o =
    ModelWithControllersOnly_cal->FirstOrderHold1_IniOut;
  ModelWithControllersOnly_DW.Uk_d = (rtInf);
  ModelWithControllersOnly_DW.Mk_n = 0.0;

  /* InitializeConditions for Integrator: '<S212>/lateral' */
  ModelWithControllersOnly_X.lateral_CSTATE =
    ModelWithControllersOnly_cal->lateral_IC;

  /* InitializeConditions for Integrator: '<S213>/lateral' */
  ModelWithControllersOnly_X.lateral_CSTATE_k =
    ModelWithControllersOnly_cal->lateral_IC_a;

  /* InitializeConditions for UnitDelay: '<S127>/Unit Delay' */
  ModelWithControllersOnly_DW.UnitDelay_DSTATE_l =
    ModelWithControllersOnly_cal->UnitDelay_InitialCondition_f;

  /* InitializeConditions for DiscreteTransferFcn: '<S127>/second-order low-pass torque model' */
  ModelWithControllersOnly_DW.secondorderlowpasstorquemodel_s[0] =
    ModelWithControllersOnly_cal->secondorderlowpasstorquemodel_I;

  /* InitializeConditions for DiscreteTransferFcn: '<S127>/second-order low-pass brake model' */
  ModelWithControllersOnly_DW.secondorderlowpassbrakemodel_st[0] =
    ModelWithControllersOnly_cal->secondorderlowpassbrakemodel_In;

  /* InitializeConditions for DiscreteTransferFcn: '<S127>/second-order low-pass steering model' */
  ModelWithControllersOnly_DW.secondorderlowpasssteeringmodel[0] =
    ModelWithControllersOnly_cal->secondorderlowpasssteeringmodel;

  /* InitializeConditions for DiscreteTransferFcn: '<S127>/second-order low-pass torque model' */
  ModelWithControllersOnly_DW.secondorderlowpasstorquemodel_s[1] =
    ModelWithControllersOnly_cal->secondorderlowpasstorquemodel_I;

  /* InitializeConditions for DiscreteTransferFcn: '<S127>/second-order low-pass brake model' */
  ModelWithControllersOnly_DW.secondorderlowpassbrakemodel_st[1] =
    ModelWithControllersOnly_cal->secondorderlowpassbrakemodel_In;

  /* InitializeConditions for DiscreteTransferFcn: '<S127>/second-order low-pass steering model' */
  ModelWithControllersOnly_DW.secondorderlowpasssteeringmodel[1] =
    ModelWithControllersOnly_cal->secondorderlowpasssteeringmodel;

  /* SystemInitialize for IfAction SubSystem: '<S11>/Forward' */
  /* InitializeConditions for DiscreteIntegrator: '<S49>/Integrator' */
  ModelWithControllersOnly_DW.Integrator_DSTATE_o =
    ModelWithControllersOnly_cal->PIForward_InitialConditionForIn;
  ModelWithControllersOnly_DW.Integrator_PrevResetState_n = 0;

  /* End of SystemInitialize for SubSystem: '<S11>/Forward' */

  /* SystemInitialize for IfAction SubSystem: '<S11>/Reverse' */
  /* InitializeConditions for DiscreteIntegrator: '<S100>/Integrator' */
  ModelWithControllersOnly_DW.Integrator_DSTATE =
    ModelWithControllersOnly_cal->PIReverse_InitialConditionForIn;
  ModelWithControllersOnly_DW.Integrator_PrevResetState = 0;

  /* End of SystemInitialize for SubSystem: '<S11>/Reverse' */

  /* SystemInitialize for Merge: '<S11>/Merge' */
  ModelWithControllersOnly_B.Merge =
    ModelWithControllersOnly_cal->Merge_InitialOutput;

  /* set "at time zero" to false */
  if (rtmIsFirstInitCond(ModelWithControllersOnly_M)) {
    rtmSetFirstInitCond(ModelWithControllersOnly_M, 0);
  }
}

/* Model terminate function */
void ModelWithControllersOnly_terminate(void)
{
  /* Terminate for RateTransition generated from: '<S12>/Equal1' */
  rtw_slrealtime_mutex_destroy
    (ModelWithControllersOnly_DW.TmpRTBAtEqual1Inport1_d0_SEMAPH);

  /* Terminate for RateTransition generated from: '<S12>/Equal2' */
  rtw_slrealtime_mutex_destroy
    (ModelWithControllersOnly_DW.TmpRTBAtEqual2Inport1_d0_SEMAPH);

  /* Terminate for RateTransition generated from: '<S12>/Equal4' */
  rtw_slrealtime_mutex_destroy
    (ModelWithControllersOnly_DW.TmpRTBAtEqual4Inport2_d0_SEMAPH);

  /* Terminate for RateTransition generated from: '<S12>/Equal5' */
  rtw_slrealtime_mutex_destroy
    (ModelWithControllersOnly_DW.TmpRTBAtEqual5Inport2_d0_SEMAPH);

  /* Terminate for RateTransition generated from: '<S10>/Multiply' */
  rtw_slrealtime_mutex_destroy
    (ModelWithControllersOnly_DW.TmpRTBAtMultiplyInport1_d0_SEMA);

  /* Terminate for RateTransition generated from: '<S11>/Switch Case' */
  rtw_slrealtime_mutex_destroy
    (ModelWithControllersOnly_DW.TmpRTBAtSwitchCaseInport1_d0_SE);

  /* Terminate for RateTransition generated from: '<S120>/MATLAB Function1' */
  rtw_slrealtime_mutex_destroy
    (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunction1Inport_k);

  /* Terminate for RateTransition generated from: '<S120>/MATLAB Function' */
  rtw_slrealtime_mutex_destroy
    (ModelWithControllersOnly_DW.TmpRTBAtMATLABFunctionInport3_d);

  /* Terminate for Assertion: '<S12>/Assertion' */
  slTestTerminate(&ModelWithControllersOnly_DW.Assertion_sltestBlkInfo);

  /* Terminate for Assertion: '<S12>/Assertion1' */
  slTestTerminate(&ModelWithControllersOnly_DW.Assertion1_sltestBlkInfo);
}
