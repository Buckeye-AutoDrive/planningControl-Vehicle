/*
 * speedgoat_target_model_2021b.cpp
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

#include "speedgoat_target_model_2021b.h"
#include "speedgoat_target_model_2021b_private.h"

/* Named constants for MATLAB Function: '<S43>/optimizer' */
const real_T speedgoat_target_model_RMDscale = 0.02;
const real_T speedgoat_target_model_RMVscale = 0.0125;

/* Named constants for Chart: '<S1>/Chart' */
const uint32_T speedgoat_t_IN_SceneRecognition = 5U;
const uint32_T speedgoat_targe_IN_RailCrossing = 3U;
const uint32_T speedgoat_targe_IN_TrafficLight = 7U;
const uint32_T speedgoat_target_IN_Approaching = 1U;
const uint32_T speedgoat_target__IN_Navigation = 1U;
const uint32_T speedgoat_target__IN_RoundAbout = 4U;
const uint32_T speedgoat_target__IN_Stationary = 2U;
const uint32_T speedgoat_target__IN_Transition = 4U;
const uint32_T speedgoat_target_m_IN_DestCheck = 1U;
const uint32_T speedgoat_target_m_IN_ExitScene = 2U;
const uint32_T speedgoat_target_mo_IN_StopSign = 6U;
const uint32_T speedgoat_target_model_20_IN_Go = 2U;
const uint32_T speedgoat_target_model__IN_Stop = 3U;

/* Named constants for Chart: '<S1>/Collision Avoidance' */
const uint32_T speedgoat_target__IN_LaneChange = 2U;
const uint32_T speedgoat_target_mode_IN_Follow = 1U;
const uint32_T speedgoat_target_model_IN_Start = 3U;

/* Block signals (default storage) */
B_speedgoat_target_model_2021_T speedgoat_target_model_2021b_B;

/* Block states (default storage) */
DW_speedgoat_target_model_202_T speedgoat_target_model_2021b_DW;

/* Real-time model */
RT_MODEL_speedgoat_target_mod_T speedgoat_target_model_2021b_M_ =
  RT_MODEL_speedgoat_target_mod_T();
RT_MODEL_speedgoat_target_mod_T *const speedgoat_target_model_2021b_M =
  &speedgoat_target_model_2021b_M_;

/* Forward declaration for local functions */
static void speedg_angleUtilities_wrapTo2Pi(real_T *theta);
static real_T speedgoat_target_model__xnrm2_g(const real_T x[4]);
static void speedgoat_target_model_20_xrotg(real_T *a, real_T *b, real_T *c,
  real_T *s);
static void speedgoat_target_model_2021_svd(const real_T A[4], real_T U[2]);
static void speedgoat_target_m_TrafficLight(void);
static void speedgoat_target_model_StopSign(void);
static void speedgoat_target_mod_Navigation(void);
static void speedgoat_target__Unconstrained(const real_T b_Hinv[9], const real_T
  f[3], real_T x[3], int16_T n);
static real_T speedgoat_target_model_202_norm(const real_T x[3]);
static void speedgoat_target_model_2021_abs(const real_T x[3], real_T y[3]);
static real_T speedgoat_target_model__maximum(const real_T x[3]);
static void speedgoat_target_model_20_abs_l(const real_T x[34], real_T y[34]);
static void speedgoat_target_model_maximum2(const real_T x[34], real_T y, real_T
  ex[34]);
static real_T speedgoat_target_model_20_xnrm2(int32_T n, const real_T x[9],
  int32_T ix0);
static void speedgoat_target_model_20_xgemv(int32_T b_m, int32_T n, const real_T
  b_A[9], int32_T ia0, const real_T x[9], int32_T ix0, real_T y[3]);
static void speedgoat_target_model_20_xgerc(int32_T b_m, int32_T n, real_T
  alpha1, int32_T ix0, const real_T y[3], real_T b_A[9], int32_T ia0);
static void speedgoat_target_model_2021b_qr(const real_T b_A[9], real_T Q[9],
  real_T R[9]);
static void speedgoat_target_mod_KWIKfactor(const real_T b_Ac[102], const
  int16_T iC[34], int16_T nA, const real_T b_Linv[9], real_T D[9], real_T b_H[9],
  int16_T n, real_T RLinv[9], real_T *Status);
static real_T speedgoat_target_model_2_mtimes(const real_T b_A[3], const real_T
  B[3]);
static void speedgoat_target_DropConstraint(int16_T kDrop, int16_T iA[34],
  int16_T *nA, int16_T iC[34]);
static void speedgoat_target_model_2_qpkwik(const real_T b_Linv[9], const real_T
  b_Hinv[9], const real_T f[3], const real_T b_Ac[102], const real_T b[34],
  int16_T iA[34], int16_T maxiter, real_T FeasTol, real_T x[3], real_T lambda[34],
  real_T *status);

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
  (speedgoat_target_model_2021b_M->Timing.TaskCounters.TID[1])++;
  if ((speedgoat_target_model_2021b_M->Timing.TaskCounters.TID[1]) > 9) {/* Sample time: [1.0s, 0.0s] */
    speedgoat_target_model_2021b_M->Timing.TaskCounters.TID[1] = 0;
  }
}

/*
 * Output and update for atomic system:
 *    '<S4>/DataTypeConversion_L0'
 *    '<S4>/DataTypeConversion_dmin'
 *    '<S4>/DataTypeConversion_reldist'
 *    '<S4>/DataTypeConversion_vego'
 *    '<S4>/DataTypeConversion_vlead'
 *    '<S4>/DataTypeConversion_vset'
 */
void speed_DataTypeConversion_L0(real_T rtu_u, B_DataTypeConversion_L0_speed_T
  *localB)
{
  localB->y = rtu_u;
}

/*
 * Output and update for atomic system:
 *    '<S4>/DataTypeConversion_amax'
 *    '<S4>/DataTypeConversion_amin'
 *    '<S4>/DataTypeConversion_atrack'
 */
void spe_DataTypeConversion_amax(real_T rtu_u, B_DataTypeConversion_amax_spe_T
  *localB)
{
  localB->y = rtu_u;
}

/* Function for MATLAB Function: '<S7>/Kinematic' */
static void speedg_angleUtilities_wrapTo2Pi(real_T *theta)
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

/* Function for MATLAB Function: '<S1>/MATLAB Function' */
static real_T speedgoat_target_model__xnrm2_g(const real_T x[4])
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  scale = 3.3121686421112381E-170;
  absxk = std::abs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }

  absxk = std::abs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  return scale * std::sqrt(y);
}

/* Function for MATLAB Function: '<S1>/MATLAB Function' */
static void speedgoat_target_model_20_xrotg(real_T *a, real_T *b, real_T *c,
  real_T *s)
{
  real_T absa;
  real_T absb;
  real_T roe;
  real_T scale;
  roe = *b;
  absa = std::abs(*a);
  absb = std::abs(*b);
  if (absa > absb) {
    roe = *a;
  }

  scale = absa + absb;
  if (scale == 0.0) {
    *s = 0.0;
    *c = 1.0;
    *a = 0.0;
    *b = 0.0;
  } else {
    real_T ads;
    real_T bds;
    ads = absa / scale;
    bds = absb / scale;
    scale *= std::sqrt(ads * ads + bds * bds);
    if (roe < 0.0) {
      scale = -scale;
    }

    *c = *a / scale;
    *s = *b / scale;
    if (absa > absb) {
      *b = *s;
    } else if (*c != 0.0) {
      *b = 1.0 / *c;
    } else {
      *b = 1.0;
    }

    *a = scale;
  }
}

/* Function for MATLAB Function: '<S1>/MATLAB Function' */
static void speedgoat_target_model_2021_svd(const real_T A[4], real_T U[2])
{
  real_T e[2];
  real_T s[2];
  real_T nrm;
  real_T rt;
  real_T smm1;
  real_T sqds;
  real_T ztest;
  int32_T b_ix;
  int32_T b_iy;
  rt = A[0];
  ztest = A[1];
  sqds = A[2];
  smm1 = A[3];
  nrm = speedgoat_target_model__xnrm2_g(A);
  if (nrm > 0.0) {
    if (A[0] < 0.0) {
      s[0] = -nrm;
    } else {
      s[0] = nrm;
    }

    if (std::abs(s[0]) >= 1.0020841800044864E-292) {
      nrm = 1.0 / s[0];
      rt *= nrm;
      ztest *= nrm;
    } else {
      rt /= s[0];
      ztest /= s[0];
    }

    rt++;
    s[0] = -s[0];
    nrm = rt * sqds;
    nrm += ztest * smm1;
    nrm = -(nrm / rt);
    if (!(nrm == 0.0)) {
      sqds += nrm * rt;
      smm1 += nrm * ztest;
    }
  } else {
    s[0] = 0.0;
  }

  b_iy = 0;
  s[1] = smm1;
  e[0] = sqds;
  e[1] = 0.0;
  ztest = e[0];
  if (s[0] != 0.0) {
    rt = std::abs(s[0]);
    nrm = s[0] / rt;
    s[0] = rt;
    ztest /= nrm;
  }

  if (ztest != 0.0) {
    rt = std::abs(ztest);
    nrm = rt / ztest;
    ztest = rt;
    s[1] *= nrm;
  }

  e[0] = ztest;
  ztest = e[0];
  if (s[1] != 0.0) {
    rt = std::abs(s[1]);
    s[1] = rt;
  }

  e[0] = ztest;
  b_ix = 0;
  ztest = s[0];
  nrm = e[0];
  if ((ztest >= nrm) || rtIsNaN(nrm)) {
    nrm = ztest;
  }

  ztest = s[1];
  if (!(ztest >= 0.0)) {
    ztest = 0.0;
  }

  if ((!(nrm >= ztest)) && (!rtIsNaN(ztest))) {
    nrm = ztest;
  }

  while ((b_iy + 2 > 0) && (b_ix < 75)) {
    int32_T kase;
    int32_T q;
    int32_T qs;
    boolean_T exitg1;
    q = b_iy + 1;
    exitg1 = false;
    while (!(exitg1 || (q == 0))) {
      rt = std::abs(e[0]);
      if ((rt <= (std::abs(s[0]) + std::abs(s[1])) * 2.2204460492503131E-16) ||
          ((rt <= 1.0020841800044864E-292) || ((b_ix > 20) && (rt <=
             2.2204460492503131E-16 * nrm)))) {
        e[0] = 0.0;
        exitg1 = true;
      } else {
        q = 0;
      }
    }

    if (b_iy + 1 == q) {
      kase = 4;
    } else {
      qs = b_iy + 2;
      kase = b_iy + 2;
      exitg1 = false;
      while ((!exitg1) && (kase >= q)) {
        qs = kase;
        if (kase == q) {
          exitg1 = true;
        } else {
          rt = 0.0;
          if (kase < b_iy + 2) {
            rt = std::abs(e[kase - 1]);
          }

          if (kase > q + 1) {
            rt += std::abs(e[0]);
          }

          ztest = std::abs(s[kase - 1]);
          if ((ztest <= 2.2204460492503131E-16 * rt) || (ztest <=
               1.0020841800044864E-292)) {
            s[kase - 1] = 0.0;
            exitg1 = true;
          } else {
            kase--;
          }
        }
      }

      if (qs == q) {
        kase = 3;
      } else if (b_iy + 2 == qs) {
        kase = 1;
      } else {
        kase = 2;
        q = qs;
      }
    }

    switch (kase) {
     case 1:
      rt = e[b_iy];
      e[b_iy] = 0.0;
      qs = b_iy + 1;
      while (qs >= q + 1) {
        speedgoat_target_model_20_xrotg(&s[0], &rt, &ztest, &sqds);
        qs = 0;
      }
      break;

     case 2:
      rt = e[q - 1];
      e[q - 1] = 0.0;
      while (q + 1 <= b_iy + 2) {
        speedgoat_target_model_20_xrotg(&s[q], &rt, &ztest, &sqds);
        rt = -sqds * e[q];
        e[q] *= ztest;
        q++;
      }
      break;

     case 3:
      {
        real_T emm1;
        real_T shift;
        ztest = std::abs(s[b_iy + 1]);
        rt = std::abs(s[b_iy]);
        if ((ztest >= rt) || rtIsNaN(rt)) {
          rt = ztest;
        }

        ztest = std::abs(e[b_iy]);
        if ((rt >= ztest) || rtIsNaN(ztest)) {
          ztest = rt;
        }

        rt = std::abs(s[q]);
        if ((ztest >= rt) || rtIsNaN(rt)) {
          rt = ztest;
        }

        ztest = std::abs(e[q]);
        if ((rt >= ztest) || rtIsNaN(ztest)) {
          ztest = rt;
        }

        rt = s[b_iy + 1] / ztest;
        smm1 = s[b_iy] / ztest;
        emm1 = e[b_iy] / ztest;
        sqds = s[q] / ztest;
        smm1 = ((smm1 + rt) * (smm1 - rt) + emm1 * emm1) / 2.0;
        emm1 *= rt;
        emm1 *= emm1;
        if ((smm1 != 0.0) || (emm1 != 0.0)) {
          shift = std::sqrt(smm1 * smm1 + emm1);
          if (smm1 < 0.0) {
            shift = -shift;
          }

          shift = emm1 / (smm1 + shift);
        } else {
          shift = 0.0;
        }

        rt = (sqds + rt) * (sqds - rt) + shift;
        ztest = e[q] / ztest * sqds;
        while (q + 1 <= b_iy + 1) {
          speedgoat_target_model_20_xrotg(&rt, &ztest, &sqds, &smm1);
          rt = sqds * s[0] + smm1 * e[0];
          e[0] = sqds * e[0] - smm1 * s[0];
          ztest = smm1 * s[1];
          s[1] *= sqds;
          s[0] = rt;
          speedgoat_target_model_20_xrotg(&s[0], &ztest, &sqds, &smm1);
          rt = sqds * e[0] + smm1 * s[1];
          s[1] = -smm1 * e[0] + sqds * s[1];
          ztest = smm1 * e[1];
          e[1] *= sqds;
          q = 1;
        }

        e[b_iy] = rt;
        b_ix++;
      }
      break;

     default:
      if (s[q] < 0.0) {
        s[q] = -s[q];
      }

      b_ix = q + 1;
      while ((q + 1 < 2) && (s[q] < s[b_ix])) {
        rt = s[q];
        s[q] = s[b_ix];
        s[b_ix] = rt;
        q = b_ix;
        b_ix++;
      }

      b_ix = 0;
      b_iy--;
      break;
    }
  }

  U[0] = s[0];
  U[1] = s[1];
}

/* Function for Chart: '<S1>/Chart' */
static void speedgoat_target_m_TrafficLight(void)
{
  if (speedgoat_target_model_2021b_DW.ExitFlag == 1.0) {
    speedgoat_target_model_2021b_DW.is_TrafficLight = 0U;
    speedgoat_target_model_2021b_DW.is_Navigation =
      speedgoat_target_model_20_IN_Go;
    speedgoat_target_model_2021b_B.VelCmd =
      speedgoat_target_model_2021b_B.RoadSpeed;
    speedgoat_target_model_2021b_DW.CurrentVelCmd =
      speedgoat_target_model_2021b_B.VelCmd;
    speedgoat_target_model_2021b_B.Trans = 0.0;
    speedgoat_target_model_2021b_DW.ExitFlag = 0.0;
    speedgoat_target_model_2021b_B.ACC_i = 0.0;
  } else {
    switch (speedgoat_target_model_2021b_DW.is_TrafficLight) {
     case speedgoat_target_IN_Approaching:
      speedgoat_target_model_2021b_B.Tgap = 0.0;
      speedgoat_target_model_2021b_B.ACC_i = 1.0;
      if ((speedgoat_target_model_2021b_B.D2S <= 4.0) &&
          (speedgoat_target_model_2021b_B.TLstate[54] < 3.0)) {
        speedgoat_target_model_2021b_DW.is_TrafficLight =
          speedgoat_target_model__IN_Stop;
        speedgoat_target_model_2021b_B.VelCmd = 0.0;

        /* SystemInitialize for S-Function (slrealtimebytepacking): '<S168>/Byte Unpacking ' */
      } else if (speedgoat_target_model_2021b_B.TLstate[54] == 3.0) {
        speedgoat_target_model_2021b_DW.is_TrafficLight =
          speedgoat_target__IN_Transition;
        speedgoat_target_model_2021b_B.Trans = 1.0;
        speedgoat_target_model_2021b_B.ACC_i = 1.0;
        speedgoat_target_model_2021b_B.Tgap = 1.0;
        speedgoat_target_model_2021b_B.VelCmd =
          speedgoat_target_model_2021b_B.SceneSpeed;
        speedgoat_target_model_2021b_DW.CurrentVelCmd =
          speedgoat_target_model_2021b_B.VelCmd;
      } else if (speedgoat_target_model_2021b_B.SceneID != 8.0) {
        speedgoat_target_model_2021b_DW.is_TrafficLight =
          speedgoat_target_m_IN_ExitScene;
        speedgoat_target_model_2021b_DW.ExitFlag = 1.0;
      }
      break;

     case speedgoat_target_m_IN_ExitScene:
      break;

     case speedgoat_target_model__IN_Stop:
      /* SystemInitialize for S-Function (slrealtimebytepacking): '<S168>/Byte Unpacking ' */
      if (speedgoat_target_model_2021b_B.TLstate[54] == 3.0) {
        speedgoat_target_model_2021b_DW.is_TrafficLight =
          speedgoat_target__IN_Transition;
        speedgoat_target_model_2021b_B.Trans = 1.0;
        speedgoat_target_model_2021b_B.ACC_i = 1.0;
        speedgoat_target_model_2021b_B.Tgap = 1.0;
        speedgoat_target_model_2021b_B.VelCmd =
          speedgoat_target_model_2021b_B.SceneSpeed;
        speedgoat_target_model_2021b_DW.CurrentVelCmd =
          speedgoat_target_model_2021b_B.VelCmd;
      }
      break;

     default:
      /* case IN_Transition: */
      speedgoat_target_model_2021b_B.ACC_i = 1.0;
      speedgoat_target_model_2021b_B.Tgap = 1.0;
      if (speedgoat_target_model_2021b_B.SceneID == 0.0) {
        speedgoat_target_model_2021b_DW.is_TrafficLight =
          speedgoat_target_m_IN_ExitScene;
        speedgoat_target_model_2021b_DW.ExitFlag = 1.0;
      } else if (speedgoat_target_model_2021b_B.TLstate[54] == 1.0) {
        speedgoat_target_model_2021b_DW.is_TrafficLight =
          speedgoat_target_IN_Approaching;
        speedgoat_target_model_2021b_B.Tgap = 0.0;
        speedgoat_target_model_2021b_B.Trans = 0.0;
        speedgoat_target_model_2021b_B.ACC_i = 1.0;
      }
      break;
    }
  }
}

/* Function for Chart: '<S1>/Chart' */
static void speedgoat_target_model_StopSign(void)
{
  if (speedgoat_target_model_2021b_DW.ExitFlag == 1.0) {
    speedgoat_target_model_2021b_DW.is_StopSign = 0U;
    speedgoat_target_model_2021b_DW.is_Navigation =
      speedgoat_target_model_20_IN_Go;
    speedgoat_target_model_2021b_B.VelCmd =
      speedgoat_target_model_2021b_B.RoadSpeed;
    speedgoat_target_model_2021b_DW.CurrentVelCmd =
      speedgoat_target_model_2021b_B.VelCmd;
    speedgoat_target_model_2021b_B.Trans = 0.0;
    speedgoat_target_model_2021b_DW.ExitFlag = 0.0;
    speedgoat_target_model_2021b_B.ACC_i = 0.0;
  } else {
    switch (speedgoat_target_model_2021b_DW.is_StopSign) {
     case speedgoat_target_IN_Approaching:
      speedgoat_target_model_2021b_B.ACC_i = 1.0;
      speedgoat_target_model_2021b_B.Tgap = 0.0;
      if (speedgoat_target_model_2021b_B.D2S <= 7.0) {
        speedgoat_target_model_2021b_DW.is_StopSign =
          speedgoat_target_model__IN_Stop;
        speedgoat_target_model_2021b_DW.temporalCounter_i1_i = 0U;
        speedgoat_target_model_2021b_B.VelCmd = 0.0;
        speedgoat_target_model_2021b_B.ACC_i = 0.0;
        speedgoat_target_model_2021b_DW.CurrentVelCmd =
          speedgoat_target_model_2021b_B.VelCmd;
      } else if (speedgoat_target_model_2021b_B.SceneID != 7.0) {
        speedgoat_target_model_2021b_DW.is_StopSign =
          speedgoat_target_m_IN_ExitScene;
        speedgoat_target_model_2021b_DW.ExitFlag = 1.0;
      }
      break;

     case speedgoat_target_m_IN_ExitScene:
      break;

     case speedgoat_target_model__IN_Stop:
      {
        boolean_T out;
        speedgoat_target_model_2021b_B.ACC_i = 0.0;
        out = ((speedgoat_target_model_2021b_DW.temporalCounter_i1_i >= 30U) &&
               (speedgoat_target_model_2021b_B.Clear2Go == 1.0));
        if (out) {
          speedgoat_target_model_2021b_DW.is_StopSign =
            speedgoat_target__IN_Transition;
          speedgoat_target_model_2021b_B.Trans = 1.0;
          speedgoat_target_model_2021b_B.ACC_i = 0.0;
          speedgoat_target_model_2021b_B.Tgap = 1.0;
          speedgoat_target_model_2021b_B.VelCmd =
            speedgoat_target_model_2021b_B.SceneSpeed;
          speedgoat_target_model_2021b_DW.CurrentVelCmd =
            speedgoat_target_model_2021b_B.VelCmd;
        }
      }
      break;

     default:
      /* case IN_Transition: */
      speedgoat_target_model_2021b_B.ACC_i = 0.0;
      speedgoat_target_model_2021b_B.Tgap = 1.0;
      if (speedgoat_target_model_2021b_B.SceneID == 0.0) {
        speedgoat_target_model_2021b_DW.is_StopSign =
          speedgoat_target_m_IN_ExitScene;
        speedgoat_target_model_2021b_DW.ExitFlag = 1.0;
      }
      break;
    }
  }
}

/* Function for Chart: '<S1>/Chart' */
static void speedgoat_target_mod_Navigation(void)
{
  boolean_T out;
  out = ((speedgoat_target_model_2021b_B.NavigatingFlag == 0.0) ||
         (speedgoat_target_model_2021b_DW.DestReached == 1.0));
  if (out) {
    speedgoat_target_model_2021b_DW.is_RailCrossing = 0U;
    speedgoat_target_model_2021b_DW.is_RoundAbout = 0U;
    speedgoat_target_model_2021b_DW.is_StopSign = 0U;
    speedgoat_target_model_2021b_DW.is_TrafficLight = 0U;
    speedgoat_target_model_2021b_DW.is_Navigation = 0U;
    speedgoat_target_model_2021b_DW.is_c21_speedgoat_target_model_2 =
      speedgoat_target__IN_Stationary;
  } else {
    switch (speedgoat_target_model_2021b_DW.is_Navigation) {
     case speedgoat_target_m_IN_DestCheck:
      if (speedgoat_target_model_2021b_DW.DestReached == 0.0) {
        speedgoat_target_model_2021b_DW.is_Navigation =
          speedgoat_t_IN_SceneRecognition;
        speedgoat_target_model_2021b_B.VelCmd =
          speedgoat_target_model_2021b_B.RoadSpeed;
        speedgoat_target_model_2021b_DW.CurrentVelCmd =
          speedgoat_target_model_2021b_B.VelCmd;
        speedgoat_target_model_2021b_DW.ExitFlag = 0.0;
      }
      break;

     case speedgoat_target_model_20_IN_Go:
      speedgoat_target_model_2021b_B.ACC_i = 0.0;
      if (speedgoat_target_model_2021b_B.SceneID == 0.0) {
        speedgoat_target_model_2021b_DW.is_Navigation =
          speedgoat_target_m_IN_DestCheck;
        if (speedgoat_target_model_2021b_B.D2D < 2.0) {
          speedgoat_target_model_2021b_DW.DestReached = 1.0;
          speedgoat_target_model_2021b_B.StopSim = 1.0;
        } else {
          speedgoat_target_model_2021b_DW.DestReached = 0.0;
        }
      } else if (speedgoat_target_model_2021b_B.SceneID == 10.0) {
        speedgoat_target_model_2021b_DW.is_Navigation =
          speedgoat_targe_IN_RailCrossing;
        speedgoat_target_model_2021b_DW.is_RailCrossing =
          speedgoat_target_IN_Approaching;
        speedgoat_target_model_2021b_B.ACC_i = 1.0;
        speedgoat_target_model_2021b_B.Tgap = 0.0;
        speedgoat_target_model_2021b_B.Trans = 0.0;
      } else if (speedgoat_target_model_2021b_B.SceneID == 8.0) {
        speedgoat_target_model_2021b_DW.is_Navigation =
          speedgoat_targe_IN_TrafficLight;
        speedgoat_target_model_2021b_DW.is_TrafficLight =
          speedgoat_target_IN_Approaching;
        speedgoat_target_model_2021b_B.Tgap = 0.0;
        speedgoat_target_model_2021b_B.Trans = 0.0;
        speedgoat_target_model_2021b_B.ACC_i = 1.0;
      } else if (speedgoat_target_model_2021b_B.SceneID == 9.0) {
        speedgoat_target_model_2021b_DW.is_Navigation =
          speedgoat_target__IN_RoundAbout;
        speedgoat_target_model_2021b_DW.is_RoundAbout =
          speedgoat_target_IN_Approaching;
        speedgoat_target_model_2021b_B.ACC_i = 1.0;
        speedgoat_target_model_2021b_B.Trans = 0.0;
        speedgoat_target_model_2021b_B.Tgap = 0.0;
      } else if (speedgoat_target_model_2021b_B.SceneID == 7.0) {
        speedgoat_target_model_2021b_DW.is_Navigation =
          speedgoat_target_mo_IN_StopSign;
        speedgoat_target_model_2021b_DW.is_StopSign =
          speedgoat_target_IN_Approaching;
        speedgoat_target_model_2021b_B.ACC_i = 1.0;
        speedgoat_target_model_2021b_B.Trans = 0.0;
        speedgoat_target_model_2021b_B.Tgap = 0.0;
      }
      break;

     case speedgoat_targe_IN_RailCrossing:
      if (speedgoat_target_model_2021b_DW.ExitFlag == 1.0) {
        speedgoat_target_model_2021b_DW.is_RailCrossing = 0U;
        speedgoat_target_model_2021b_DW.is_Navigation =
          speedgoat_target_model_20_IN_Go;
        speedgoat_target_model_2021b_B.VelCmd =
          speedgoat_target_model_2021b_B.RoadSpeed;
        speedgoat_target_model_2021b_DW.CurrentVelCmd =
          speedgoat_target_model_2021b_B.VelCmd;
        speedgoat_target_model_2021b_B.Trans = 0.0;
        speedgoat_target_model_2021b_DW.ExitFlag = 0.0;
        speedgoat_target_model_2021b_B.ACC_i = 0.0;
      } else {
        switch (speedgoat_target_model_2021b_DW.is_RailCrossing) {
         case speedgoat_target_IN_Approaching:
          speedgoat_target_model_2021b_B.ACC_i = 1.0;
          speedgoat_target_model_2021b_B.Tgap = 0.0;
          if (speedgoat_target_model_2021b_B.D2S <= 4.0) {
            speedgoat_target_model_2021b_DW.is_RailCrossing =
              speedgoat_target_model__IN_Stop;
            speedgoat_target_model_2021b_B.VelCmd = 0.0;
            speedgoat_target_model_2021b_B.ACC_i = 0.0;
            speedgoat_target_model_2021b_DW.CurrentVelCmd =
              speedgoat_target_model_2021b_B.VelCmd;
          } else if (speedgoat_target_model_2021b_B.SceneID != 10.0) {
            speedgoat_target_model_2021b_DW.is_RailCrossing =
              speedgoat_target_m_IN_ExitScene;
            speedgoat_target_model_2021b_DW.ExitFlag = 1.0;
          }
          break;

         case speedgoat_target_m_IN_ExitScene:
          break;

         case speedgoat_target_model__IN_Stop:
          speedgoat_target_model_2021b_B.ACC_i = 0.0;
          out = ((speedgoat_target_model_2021b_B.RRCrossing == 0.0) &&
                 (speedgoat_target_model_2021b_B.ObjectDist > 14.0));
          if (out) {
            speedgoat_target_model_2021b_DW.is_RailCrossing =
              speedgoat_target__IN_Transition;
            speedgoat_target_model_2021b_B.ACC_i = 1.0;
            speedgoat_target_model_2021b_B.Tgap = 1.0;
            speedgoat_target_model_2021b_B.Trans = 1.0;
            speedgoat_target_model_2021b_B.VelCmd =
              speedgoat_target_model_2021b_B.SceneSpeed;
            speedgoat_target_model_2021b_DW.CurrentVelCmd =
              speedgoat_target_model_2021b_B.VelCmd;
          }
          break;

         default:
          /* case IN_Transition: */
          speedgoat_target_model_2021b_B.ACC_i = 1.0;
          speedgoat_target_model_2021b_B.Tgap = 1.0;
          if (speedgoat_target_model_2021b_B.SceneID == 0.0) {
            speedgoat_target_model_2021b_DW.is_RailCrossing =
              speedgoat_target_m_IN_ExitScene;
            speedgoat_target_model_2021b_DW.ExitFlag = 1.0;
          }
          break;
        }
      }
      break;

     case speedgoat_target__IN_RoundAbout:
      if (speedgoat_target_model_2021b_DW.ExitFlag == 1.0) {
        speedgoat_target_model_2021b_DW.is_RoundAbout = 0U;
        speedgoat_target_model_2021b_DW.is_Navigation =
          speedgoat_target_model_20_IN_Go;
        speedgoat_target_model_2021b_B.VelCmd =
          speedgoat_target_model_2021b_B.RoadSpeed;
        speedgoat_target_model_2021b_DW.CurrentVelCmd =
          speedgoat_target_model_2021b_B.VelCmd;
        speedgoat_target_model_2021b_B.Trans = 0.0;
        speedgoat_target_model_2021b_DW.ExitFlag = 0.0;
        speedgoat_target_model_2021b_B.ACC_i = 0.0;
      } else {
        switch (speedgoat_target_model_2021b_DW.is_RoundAbout) {
         case speedgoat_target_IN_Approaching:
          speedgoat_target_model_2021b_B.ACC_i = 1.0;
          speedgoat_target_model_2021b_B.Tgap = 0.0;
          if (speedgoat_target_model_2021b_B.D2S <= 4.0) {
            speedgoat_target_model_2021b_DW.is_RoundAbout =
              speedgoat_target_model__IN_Stop;
            speedgoat_target_model_2021b_B.VelCmd = 0.0;
            speedgoat_target_model_2021b_B.ACC_i = 0.0;
            speedgoat_target_model_2021b_DW.CurrentVelCmd =
              speedgoat_target_model_2021b_B.VelCmd;
          } else if (speedgoat_target_model_2021b_B.SceneID != 9.0) {
            speedgoat_target_model_2021b_DW.is_RoundAbout =
              speedgoat_target_m_IN_ExitScene;
            speedgoat_target_model_2021b_DW.ExitFlag = 1.0;
          }
          break;

         case speedgoat_target_m_IN_ExitScene:
          break;

         case speedgoat_target_model__IN_Stop:
          speedgoat_target_model_2021b_B.ACC_i = 0.0;
          if (speedgoat_target_model_2021b_B.Clear2Go == 1.0) {
            speedgoat_target_model_2021b_DW.is_RoundAbout =
              speedgoat_target__IN_Transition;
            speedgoat_target_model_2021b_B.ACC_i = 1.0;
            speedgoat_target_model_2021b_B.Tgap = 1.0;
            speedgoat_target_model_2021b_B.Trans = 1.0;
            speedgoat_target_model_2021b_B.VelCmd =
              speedgoat_target_model_2021b_B.SceneSpeed;
            speedgoat_target_model_2021b_DW.CurrentVelCmd =
              speedgoat_target_model_2021b_B.VelCmd;
          }
          break;

         default:
          /* case IN_Transition: */
          speedgoat_target_model_2021b_B.ACC_i = 1.0;
          speedgoat_target_model_2021b_B.Tgap = 1.0;
          if (speedgoat_target_model_2021b_B.SceneID == 0.0) {
            speedgoat_target_model_2021b_DW.is_RoundAbout =
              speedgoat_target_m_IN_ExitScene;
            speedgoat_target_model_2021b_DW.ExitFlag = 1.0;
          }
          break;
        }
      }
      break;

     case speedgoat_t_IN_SceneRecognition:
      if (speedgoat_target_model_2021b_B.SceneID == 10.0) {
        speedgoat_target_model_2021b_DW.is_Navigation =
          speedgoat_targe_IN_RailCrossing;
        speedgoat_target_model_2021b_DW.is_RailCrossing =
          speedgoat_target_IN_Approaching;
        speedgoat_target_model_2021b_B.ACC_i = 1.0;
        speedgoat_target_model_2021b_B.Tgap = 0.0;
        speedgoat_target_model_2021b_B.Trans = 0.0;
      } else if (speedgoat_target_model_2021b_B.SceneID == 9.0) {
        speedgoat_target_model_2021b_DW.is_Navigation =
          speedgoat_target__IN_RoundAbout;
        speedgoat_target_model_2021b_DW.is_RoundAbout =
          speedgoat_target_IN_Approaching;
        speedgoat_target_model_2021b_B.ACC_i = 1.0;
        speedgoat_target_model_2021b_B.Trans = 0.0;
        speedgoat_target_model_2021b_B.Tgap = 0.0;
      } else if (speedgoat_target_model_2021b_B.SceneID == 7.0) {
        speedgoat_target_model_2021b_DW.is_Navigation =
          speedgoat_target_mo_IN_StopSign;
        speedgoat_target_model_2021b_DW.is_StopSign =
          speedgoat_target_IN_Approaching;
        speedgoat_target_model_2021b_B.ACC_i = 1.0;
        speedgoat_target_model_2021b_B.Trans = 0.0;
        speedgoat_target_model_2021b_B.Tgap = 0.0;
      } else if (speedgoat_target_model_2021b_B.SceneID == 8.0) {
        speedgoat_target_model_2021b_DW.is_Navigation =
          speedgoat_targe_IN_TrafficLight;
        speedgoat_target_model_2021b_DW.is_TrafficLight =
          speedgoat_target_IN_Approaching;
        speedgoat_target_model_2021b_B.Tgap = 0.0;
        speedgoat_target_model_2021b_B.Trans = 0.0;
        speedgoat_target_model_2021b_B.ACC_i = 1.0;
      } else if (speedgoat_target_model_2021b_B.D2D < 2.0) {
        speedgoat_target_model_2021b_DW.is_Navigation =
          speedgoat_target_m_IN_DestCheck;
        if (speedgoat_target_model_2021b_B.D2D < 2.0) {
          speedgoat_target_model_2021b_DW.DestReached = 1.0;
          speedgoat_target_model_2021b_B.StopSim = 1.0;
        } else {
          speedgoat_target_model_2021b_DW.DestReached = 0.0;
        }
      }
      break;

     case speedgoat_target_mo_IN_StopSign:
      speedgoat_target_model_StopSign();
      break;

     default:
      /* case IN_TrafficLight: */
      speedgoat_target_m_TrafficLight();
      break;
    }
  }
}

real_T rt_roundd_snf(real_T u)
{
  real_T y;
  if (std::abs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = std::floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = std::ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

/* Function for MATLAB Function: '<S43>/optimizer' */
static void speedgoat_target__Unconstrained(const real_T b_Hinv[9], const real_T
  f[3], real_T x[3], int16_T n)
{
  for (int32_T i = 1; i - 1 < n; i++) {
    real_T b_Hinv_0;
    int32_T i_0;
    i_0 = static_cast<int16_T>(i);
    b_Hinv_0 = -b_Hinv[i_0 - 1] * f[0];
    b_Hinv_0 += -b_Hinv[i_0 + 2] * f[1];
    b_Hinv_0 += -b_Hinv[i_0 + 5] * f[2];
    x[static_cast<int16_T>(i) - 1] = b_Hinv_0;
  }
}

/* Function for MATLAB Function: '<S43>/optimizer' */
static real_T speedgoat_target_model_202_norm(const real_T x[3])
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  scale = 3.3121686421112381E-170;
  absxk = std::abs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }

  absxk = std::abs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  absxk = std::abs(x[2]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  return scale * std::sqrt(y);
}

/* Function for MATLAB Function: '<S43>/optimizer' */
static void speedgoat_target_model_2021_abs(const real_T x[3], real_T y[3])
{
  y[0] = std::abs(x[0]);
  y[1] = std::abs(x[1]);
  y[2] = std::abs(x[2]);
}

/* Function for MATLAB Function: '<S43>/optimizer' */
static real_T speedgoat_target_model__maximum(const real_T x[3])
{
  real_T ex;
  int32_T idx;
  if (!rtIsNaN(x[0])) {
    idx = 1;
  } else {
    int32_T k;
    boolean_T exitg1;
    idx = 0;
    k = 2;
    exitg1 = false;
    while ((!exitg1) && (k < 4)) {
      if (!rtIsNaN(x[k - 1])) {
        idx = k;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  if (idx == 0) {
    ex = x[0];
  } else {
    ex = x[idx - 1];
    while (idx + 1 <= 3) {
      if (ex < x[idx]) {
        ex = x[idx];
      }

      idx++;
    }
  }

  return ex;
}

/* Function for MATLAB Function: '<S43>/optimizer' */
static void speedgoat_target_model_20_abs_l(const real_T x[34], real_T y[34])
{
  for (int32_T k = 0; k < 34; k++) {
    y[k] = std::abs(x[k]);
  }
}

/* Function for MATLAB Function: '<S43>/optimizer' */
static void speedgoat_target_model_maximum2(const real_T x[34], real_T y, real_T
  ex[34])
{
  for (int32_T k = 0; k < 34; k++) {
    real_T u0;
    u0 = x[k];
    if ((!(u0 >= y)) && (!rtIsNaN(y))) {
      u0 = y;
    }

    ex[k] = u0;
  }
}

/* Function for MATLAB Function: '<S43>/optimizer' */
static real_T speedgoat_target_model_20_xnrm2(int32_T n, const real_T x[9],
  int32_T ix0)
{
  real_T y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      real_T scale;
      int32_T kend;
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (int32_T k = ix0; k <= kend; k++) {
        real_T absxk;
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          real_T t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          real_T t;
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T a;
  real_T y;
  a = std::abs(u0);
  y = std::abs(u1);
  if (a < y) {
    a /= y;
    y *= std::sqrt(a * a + 1.0);
  } else if (a > y) {
    y /= a;
    y = std::sqrt(y * y + 1.0) * a;
  } else if (!rtIsNaN(y)) {
    y = a * 1.4142135623730951;
  }

  return y;
}

/* Function for MATLAB Function: '<S43>/optimizer' */
static void speedgoat_target_model_20_xgemv(int32_T b_m, int32_T n, const real_T
  b_A[9], int32_T ia0, const real_T x[9], int32_T ix0, real_T y[3])
{
  if ((b_m != 0) && (n != 0)) {
    int32_T b;
    int32_T b_iy;
    for (b_iy = 0; b_iy < n; b_iy++) {
      y[b_iy] = 0.0;
    }

    b_iy = 0;
    b = (n - 1) * 3 + ia0;
    for (int32_T iac = ia0; iac <= b; iac += 3) {
      real_T c;
      int32_T d;
      int32_T ix;
      ix = ix0;
      c = 0.0;
      d = (iac + b_m) - 1;
      for (int32_T ia = iac; ia <= d; ia++) {
        c += b_A[ia - 1] * x[ix - 1];
        ix++;
      }

      y[b_iy] += c;
      b_iy++;
    }
  }
}

/* Function for MATLAB Function: '<S43>/optimizer' */
static void speedgoat_target_model_20_xgerc(int32_T b_m, int32_T n, real_T
  alpha1, int32_T ix0, const real_T y[3], real_T b_A[9], int32_T ia0)
{
  if (!(alpha1 == 0.0)) {
    int32_T jA;
    int32_T jy;
    jA = ia0 - 1;
    jy = 0;
    for (int32_T j = 0; j < n; j++) {
      if (y[jy] != 0.0) {
        real_T temp;
        int32_T b;
        int32_T ijA;
        int32_T ix;
        temp = y[jy] * alpha1;
        ix = ix0;
        ijA = jA;
        b = b_m + jA;
        while (ijA + 1 <= b) {
          b_A[ijA] += b_A[ix - 1] * temp;
          ix++;
          ijA++;
        }
      }

      jy++;
      jA += 3;
    }
  }
}

/* Function for MATLAB Function: '<S43>/optimizer' */
static void speedgoat_target_model_2021b_qr(const real_T b_A[9], real_T Q[9],
  real_T R[9])
{
  real_T c_A[9];
  real_T work[3];
  real_T atmp;
  real_T tau_idx_0;
  real_T tau_idx_1;
  real_T xnorm;
  int32_T b_coltop;
  int32_T c_lastc;
  int32_T coltop;
  int32_T exitg1;
  int32_T knt;
  boolean_T exitg2;
  std::memcpy(&c_A[0], &b_A[0], 9U * sizeof(real_T));
  tau_idx_0 = 0.0;
  work[0] = 0.0;
  tau_idx_1 = 0.0;
  work[1] = 0.0;
  work[2] = 0.0;
  atmp = c_A[0];
  xnorm = speedgoat_target_model_20_xnrm2(2, c_A, 2);
  if (xnorm != 0.0) {
    xnorm = rt_hypotd_snf(c_A[0], xnorm);
    if (c_A[0] >= 0.0) {
      xnorm = -xnorm;
    }

    if (std::abs(xnorm) < 1.0020841800044864E-292) {
      knt = 0;
      do {
        knt++;
        for (b_coltop = 1; b_coltop < 3; b_coltop++) {
          c_A[b_coltop] *= 9.9792015476736E+291;
        }

        xnorm *= 9.9792015476736E+291;
        atmp *= 9.9792015476736E+291;
      } while ((std::abs(xnorm) < 1.0020841800044864E-292) && (knt < 20));

      xnorm = rt_hypotd_snf(atmp, speedgoat_target_model_20_xnrm2(2, c_A, 2));
      if (atmp >= 0.0) {
        xnorm = -xnorm;
      }

      tau_idx_0 = (xnorm - atmp) / xnorm;
      atmp = 1.0 / (atmp - xnorm);
      for (b_coltop = 1; b_coltop < 3; b_coltop++) {
        c_A[b_coltop] *= atmp;
      }

      for (c_lastc = 0; c_lastc < knt; c_lastc++) {
        xnorm *= 1.0020841800044864E-292;
      }

      atmp = xnorm;
    } else {
      tau_idx_0 = (xnorm - c_A[0]) / xnorm;
      atmp = 1.0 / (c_A[0] - xnorm);
      for (c_lastc = 1; c_lastc < 3; c_lastc++) {
        c_A[c_lastc] *= atmp;
      }

      atmp = xnorm;
    }
  }

  c_A[0] = atmp;
  xnorm = c_A[0];
  c_A[0] = 1.0;
  if (tau_idx_0 != 0.0) {
    knt = 3;
    c_lastc = 2;
    while ((knt > 0) && (c_A[c_lastc] == 0.0)) {
      knt--;
      c_lastc--;
    }

    c_lastc = 2;
    exitg2 = false;
    while ((!exitg2) && (c_lastc > 0)) {
      b_coltop = (c_lastc - 1) * 3 + 3;
      coltop = b_coltop;
      do {
        exitg1 = 0;
        if (coltop + 1 <= b_coltop + knt) {
          if (c_A[coltop] != 0.0) {
            exitg1 = 1;
          } else {
            coltop++;
          }
        } else {
          c_lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    knt = 0;
    c_lastc = 0;
  }

  if (knt > 0) {
    speedgoat_target_model_20_xgemv(knt, c_lastc, c_A, 4, c_A, 1, work);
    speedgoat_target_model_20_xgerc(knt, c_lastc, -tau_idx_0, 1, work, c_A, 4);
  }

  c_A[0] = xnorm;
  atmp = c_A[4];
  xnorm = speedgoat_target_model_20_xnrm2(1, c_A, 6);
  if (xnorm != 0.0) {
    xnorm = rt_hypotd_snf(c_A[4], xnorm);
    if (c_A[4] >= 0.0) {
      xnorm = -xnorm;
    }

    if (std::abs(xnorm) < 1.0020841800044864E-292) {
      knt = 0;
      do {
        knt++;
        for (b_coltop = 5; b_coltop < 6; b_coltop++) {
          c_A[b_coltop] *= 9.9792015476736E+291;
        }

        xnorm *= 9.9792015476736E+291;
        atmp *= 9.9792015476736E+291;
      } while ((std::abs(xnorm) < 1.0020841800044864E-292) && (knt < 20));

      xnorm = rt_hypotd_snf(atmp, speedgoat_target_model_20_xnrm2(1, c_A, 6));
      if (atmp >= 0.0) {
        xnorm = -xnorm;
      }

      tau_idx_1 = (xnorm - atmp) / xnorm;
      atmp = 1.0 / (atmp - xnorm);
      for (b_coltop = 5; b_coltop < 6; b_coltop++) {
        c_A[b_coltop] *= atmp;
      }

      for (c_lastc = 0; c_lastc < knt; c_lastc++) {
        xnorm *= 1.0020841800044864E-292;
      }

      atmp = xnorm;
    } else {
      tau_idx_1 = (xnorm - c_A[4]) / xnorm;
      atmp = 1.0 / (c_A[4] - xnorm);
      for (c_lastc = 5; c_lastc < 6; c_lastc++) {
        c_A[c_lastc] *= atmp;
      }

      atmp = xnorm;
    }
  }

  c_A[4] = atmp;
  xnorm = c_A[4];
  c_A[4] = 1.0;
  if (tau_idx_1 != 0.0) {
    knt = 2;
    c_lastc = 5;
    while ((knt > 0) && (c_A[c_lastc] == 0.0)) {
      knt--;
      c_lastc--;
    }

    c_lastc = 1;
    coltop = 7;
    do {
      exitg1 = 0;
      if (coltop + 1 <= 7 + knt) {
        if (c_A[coltop] != 0.0) {
          exitg1 = 1;
        } else {
          coltop++;
        }
      } else {
        c_lastc = 0;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    knt = 0;
    c_lastc = 0;
  }

  if (knt > 0) {
    speedgoat_target_model_20_xgemv(knt, c_lastc, c_A, 8, c_A, 5, work);
    speedgoat_target_model_20_xgerc(knt, c_lastc, -tau_idx_1, 5, work, c_A, 8);
  }

  c_A[4] = xnorm;
  R[0] = c_A[0];
  for (c_lastc = 1; c_lastc + 1 < 4; c_lastc++) {
    R[c_lastc] = 0.0;
  }

  work[0] = 0.0;
  for (c_lastc = 0; c_lastc < 2; c_lastc++) {
    R[c_lastc + 3] = c_A[c_lastc + 3];
  }

  while (c_lastc + 1 < 4) {
    R[c_lastc + 3] = 0.0;
    c_lastc++;
  }

  work[1] = 0.0;
  for (c_lastc = 0; c_lastc < 3; c_lastc++) {
    R[c_lastc + 6] = c_A[c_lastc + 6];
  }

  work[2] = 0.0;
  c_A[8] = 1.0;
  for (c_lastc = 0; c_lastc < 2; c_lastc++) {
    c_A[7 - c_lastc] = 0.0;
  }

  c_A[4] = 1.0;
  if (tau_idx_1 != 0.0) {
    b_coltop = 7;
    while ((c_lastc > 0) && (c_A[b_coltop - 2] == 0.0)) {
      c_lastc--;
      b_coltop--;
    }

    b_coltop = 1;
    knt = 8;
    do {
      exitg1 = 0;
      if (knt <= c_lastc + 7) {
        if (c_A[knt - 1] != 0.0) {
          exitg1 = 1;
        } else {
          knt++;
        }
      } else {
        b_coltop = 0;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    c_lastc = 0;
    b_coltop = 0;
  }

  if (c_lastc > 0) {
    speedgoat_target_model_20_xgemv(c_lastc, b_coltop, c_A, 8, c_A, 5, work);
    speedgoat_target_model_20_xgerc(c_lastc, b_coltop, -tau_idx_1, 5, work, c_A,
      8);
  }

  for (b_coltop = 5; b_coltop < 6; b_coltop++) {
    c_A[b_coltop] *= -tau_idx_1;
  }

  c_A[4] = 1.0 - tau_idx_1;
  c_A[3] = 0.0;
  c_A[0] = 1.0;
  if (tau_idx_0 != 0.0) {
    c_lastc = 3;
    b_coltop = 4;
    while ((c_lastc > 0) && (c_A[b_coltop - 2] == 0.0)) {
      c_lastc--;
      b_coltop--;
    }

    b_coltop = 2;
    exitg2 = false;
    while ((!exitg2) && (b_coltop > 0)) {
      coltop = (b_coltop - 1) * 3 + 4;
      knt = coltop;
      do {
        exitg1 = 0;
        if (knt <= (coltop + c_lastc) - 1) {
          if (c_A[knt - 1] != 0.0) {
            exitg1 = 1;
          } else {
            knt++;
          }
        } else {
          b_coltop--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    c_lastc = 0;
    b_coltop = 0;
  }

  if (c_lastc > 0) {
    speedgoat_target_model_20_xgemv(c_lastc, b_coltop, c_A, 4, c_A, 1, work);
    speedgoat_target_model_20_xgerc(c_lastc, b_coltop, -tau_idx_0, 1, work, c_A,
      4);
  }

  for (b_coltop = 1; b_coltop < 3; b_coltop++) {
    c_A[b_coltop] *= -tau_idx_0;
  }

  c_A[0] = 1.0 - tau_idx_0;
  for (c_lastc = 0; c_lastc < 3; c_lastc++) {
    Q[3 * c_lastc] = c_A[3 * c_lastc];
    Q[3 * c_lastc + 1] = c_A[3 * c_lastc + 1];
    Q[3 * c_lastc + 2] = c_A[3 * c_lastc + 2];
  }
}

/* Function for MATLAB Function: '<S43>/optimizer' */
static void speedgoat_target_mod_KWIKfactor(const real_T b_Ac[102], const
  int16_T iC[34], int16_T nA, const real_T b_Linv[9], real_T D[9], real_T b_H[9],
  int16_T n, real_T RLinv[9], real_T *Status)
{
  real_T QQ[9];
  real_T RR[9];
  real_T TL[9];
  int32_T f_i;
  int32_T f_i_0;
  int32_T i;
  int32_T i_0;
  *Status = 1.0;
  std::memset(&RLinv[0], 0, 9U * sizeof(real_T));
  for (i = 1; i - 1 < nA; i++) {
    f_i_0 = iC[static_cast<int16_T>(i) - 1];
    i_0 = static_cast<int16_T>(i) - 1;
    for (f_i = 0; f_i < 3; f_i++) {
      RLinv[f_i + 3 * i_0] = 0.0;
      RLinv[f_i + 3 * i_0] += b_Ac[f_i_0 - 1] * b_Linv[f_i];
      RLinv[f_i + 3 * i_0] += b_Linv[f_i + 3] * b_Ac[f_i_0 + 33];
      RLinv[f_i + 3 * i_0] += b_Linv[f_i + 6] * b_Ac[f_i_0 + 67];
    }
  }

  speedgoat_target_model_2021b_qr(RLinv, QQ, RR);
  i = 1;
  int32_T exitg1;
  do {
    exitg1 = 0;
    if (i - 1 <= nA - 1) {
      if (std::abs(RR[((static_cast<int16_T>(i) - 1) * 3 + static_cast<int16_T>
                       (i)) - 1]) < 1.0E-12) {
        *Status = -2.0;
        exitg1 = 1;
      } else {
        i++;
      }
    } else {
      int16_T b_j;
      int16_T c_k;
      for (i = 1; i - 1 < n; i++) {
        for (f_i = 1; f_i - 1 < n; f_i++) {
          real_T b_Linv_0;
          i_0 = static_cast<int16_T>(i);
          f_i_0 = static_cast<int16_T>(f_i);
          b_Linv_0 = b_Linv[(i_0 - 1) * 3] * QQ[(f_i_0 - 1) * 3];
          b_Linv_0 += b_Linv[(i_0 - 1) * 3 + 1] * QQ[(f_i_0 - 1) * 3 + 1];
          b_Linv_0 += b_Linv[(i_0 - 1) * 3 + 2] * QQ[(f_i_0 - 1) * 3 + 2];
          TL[(static_cast<int16_T>(i) + 3 * (static_cast<int16_T>(f_i) - 1)) - 1]
            = b_Linv_0;
        }
      }

      std::memset(&RLinv[0], 0, 9U * sizeof(real_T));
      for (b_j = nA; b_j > 0; b_j = static_cast<int16_T>(b_j - 1)) {
        RLinv[(b_j + 3 * (b_j - 1)) - 1] = 1.0;
        for (c_k = b_j; c_k <= nA; c_k = static_cast<int16_T>(c_k + 1)) {
          RLinv[(b_j + 3 * (c_k - 1)) - 1] /= RR[((b_j - 1) * 3 + b_j) - 1];
        }

        if (b_j > 1) {
          for (i = 1; i - 1 <= b_j - 2; i++) {
            for (c_k = b_j; c_k <= nA; c_k = static_cast<int16_T>(c_k + 1)) {
              RLinv[(static_cast<int16_T>(i) + 3 * (c_k - 1)) - 1] -= RR[((b_j -
                1) * 3 + static_cast<int16_T>(i)) - 1] * RLinv[((c_k - 1) * 3 +
                b_j) - 1];
            }
          }
        }
      }

      for (i = 1; i - 1 < n; i++) {
        for (b_j = static_cast<int16_T>(i); b_j <= n; b_j = static_cast<int16_T>
             (b_j + 1)) {
          b_H[(static_cast<int16_T>(i) + 3 * (b_j - 1)) - 1] = 0.0;
          f_i = nA + 1;
          if (f_i > 32767) {
            f_i = 32767;
          }

          for (c_k = static_cast<int16_T>(f_i); c_k <= n; c_k =
               static_cast<int16_T>(c_k + 1)) {
            b_H[(static_cast<int16_T>(i) + 3 * (b_j - 1)) - 1] -= TL[((c_k - 1) *
              3 + static_cast<int16_T>(i)) - 1] * TL[((c_k - 1) * 3 + b_j) - 1];
          }

          b_H[(b_j + 3 * (static_cast<int16_T>(i) - 1)) - 1] = b_H[((b_j - 1) *
            3 + static_cast<int16_T>(i)) - 1];
        }
      }

      for (i = 1; i - 1 < nA; i++) {
        for (f_i = 1; f_i - 1 < n; f_i++) {
          D[(static_cast<int16_T>(f_i) + 3 * (static_cast<int16_T>(i) - 1)) - 1]
            = 0.0;
          for (b_j = static_cast<int16_T>(i); b_j <= nA; b_j =
               static_cast<int16_T>(b_j + 1)) {
            D[(static_cast<int16_T>(f_i) + 3 * (static_cast<int16_T>(i) - 1)) -
              1] += TL[((b_j - 1) * 3 + static_cast<int16_T>(f_i)) - 1] * RLinv
              [((b_j - 1) * 3 + static_cast<int16_T>(i)) - 1];
          }
        }
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);
}

/* Function for MATLAB Function: '<S43>/optimizer' */
static real_T speedgoat_target_model_2_mtimes(const real_T b_A[3], const real_T
  B[3])
{
  real_T b_C;
  b_C = b_A[0] * B[0];
  b_C += b_A[1] * B[1];
  b_C += b_A[2] * B[2];
  return b_C;
}

/* Function for MATLAB Function: '<S43>/optimizer' */
static void speedgoat_target_DropConstraint(int16_T kDrop, int16_T iA[34],
  int16_T *nA, int16_T iC[34])
{
  int32_T tmp;
  iA[iC[kDrop - 1] - 1] = 0;
  if (kDrop < *nA) {
    int16_T b;
    tmp = *nA - 1;
    if (tmp < -32768) {
      tmp = -32768;
    }

    b = static_cast<int16_T>(tmp);
    for (int16_T i = kDrop; i <= b; i = static_cast<int16_T>(i + 1)) {
      iC[i - 1] = iC[i];
    }
  }

  iC[*nA - 1] = 0;
  tmp = *nA - 1;
  if (tmp < -32768) {
    tmp = -32768;
  }

  *nA = static_cast<int16_T>(tmp);
}

/* Function for MATLAB Function: '<S43>/optimizer' */
static void speedgoat_target_model_2_qpkwik(const real_T b_Linv[9], const real_T
  b_Hinv[9], const real_T f[3], const real_T b_Ac[102], const real_T b[34],
  int16_T iA[34], int16_T maxiter, real_T FeasTol, real_T x[3], real_T lambda[34],
  real_T *status)
{
  real_T cTol[34];
  real_T tmp[34];
  real_T D[9];
  real_T RLinv[9];
  real_T U[9];
  real_T b_H[9];
  real_T Opt[6];
  real_T Rhs[6];
  real_T b_Ac_0[3];
  real_T r[3];
  real_T z[3];
  real_T Xnorm0;
  real_T cMin;
  real_T rMin;
  int32_T i;
  int32_T i_0;
  int32_T kNext_0;
  int16_T iC[34];
  int16_T kDrop;
  int16_T kNext;
  int16_T nA;
  boolean_T ColdReset;
  boolean_T DualFeasible;
  boolean_T cTolComputed;
  boolean_T guard1 = false;
  *status = 1.0;
  x[0] = 0.0;
  r[0] = 0.0;
  x[1] = 0.0;
  r[1] = 0.0;
  x[2] = 0.0;
  r[2] = 0.0;
  rMin = 0.0;
  cTolComputed = false;
  for (i = 0; i < 34; i++) {
    lambda[i] = 0.0;
    cTol[i] = 1.0;
    iC[i] = 0;
  }

  nA = 0;
  for (i = 0; i < 34; i++) {
    if (iA[i] == 1) {
      i_0 = nA + 1;
      if (i_0 > 32767) {
        i_0 = 32767;
      }

      nA = static_cast<int16_T>(i_0);
      iC[nA - 1] = static_cast<int16_T>(i + 1);
    }
  }

  guard1 = false;
  if (nA > 0) {
    int32_T exitg3;
    uint16_T c_x;
    uint16_T q;
    for (i = 0; i < 6; i++) {
      Opt[i] = 0.0;
    }

    Rhs[0] = f[0];
    Rhs[3] = 0.0;
    Rhs[1] = f[1];
    Rhs[4] = 0.0;
    Rhs[2] = f[2];
    Rhs[5] = 0.0;
    DualFeasible = false;
    i_0 = 3 * nA;
    if (i_0 > 32767) {
      i_0 = 32767;
    }

    kNext = static_cast<int16_T>(i_0);
    if (kNext < 50) {
      kNext = 50;
    }

    q = static_cast<uint16_T>(kNext / 10U);
    c_x = static_cast<uint16_T>(static_cast<uint32_T>(kNext) - q * 10);
    if ((c_x > 0) && (c_x >= 5)) {
      q = static_cast<uint16_T>(q + 1);
    }

    ColdReset = false;
    do {
      exitg3 = 0;
      if ((!DualFeasible) && (nA > 0) && (static_cast<int32_T>(*status) <=
           maxiter)) {
        speedgoat_target_mod_KWIKfactor(b_Ac, iC, nA, b_Linv, D, b_H, 3, RLinv,
          &Xnorm0);
        if (Xnorm0 < 0.0) {
          if (ColdReset) {
            *status = -2.0;
            exitg3 = 2;
          } else {
            nA = 0;
            std::memset(&iA[0], 0, 34U * sizeof(int16_T));
            std::memset(&iC[0], 0, 34U * sizeof(int16_T));
            ColdReset = true;
          }
        } else {
          for (i = 1; i - 1 < nA; i++) {
            i_0 = static_cast<int16_T>(i) + 3;
            if (i_0 > 32767) {
              i_0 = 32767;
            }

            Rhs[i_0 - 1] = b[iC[static_cast<int16_T>(i) - 1] - 1];
            for (kNext = static_cast<int16_T>(i); kNext <= nA; kNext =
                 static_cast<int16_T>(kNext + 1)) {
              U[(kNext + 3 * (static_cast<int16_T>(i) - 1)) - 1] = 0.0;
              for (kNext_0 = 1; kNext_0 - 1 < nA; kNext_0++) {
                U[(kNext + 3 * (static_cast<int16_T>(i) - 1)) - 1] += RLinv[((
                  static_cast<int16_T>(kNext_0) - 1) * 3 + kNext) - 1] * RLinv
                  [((static_cast<int16_T>(kNext_0) - 1) * 3 +
                    static_cast<int16_T>(i)) - 1];
              }

              U[(static_cast<int16_T>(i) + 3 * (kNext - 1)) - 1] = U[((
                static_cast<int16_T>(i) - 1) * 3 + kNext) - 1];
            }
          }

          for (i = 0; i < 3; i++) {
            i_0 = i + 1;
            Xnorm0 = b_H[i_0 - 1] * Rhs[0];
            Xnorm0 += b_H[i_0 + 2] * Rhs[1];
            Xnorm0 += b_H[i_0 + 5] * Rhs[2];
            Opt[i] = Xnorm0;
            for (kNext_0 = 1; kNext_0 - 1 < nA; kNext_0++) {
              i_0 = static_cast<int16_T>(kNext_0) + 3;
              if (i_0 > 32767) {
                i_0 = 32767;
              }

              Opt[i] += D[(static_cast<int16_T>(kNext_0) - 1) * 3 + i] * Rhs[i_0
                - 1];
            }
          }

          for (i = 1; i - 1 < nA; i++) {
            i_0 = static_cast<int16_T>(i);
            Xnorm0 = D[(i_0 - 1) * 3] * Rhs[0];
            Xnorm0 += D[(i_0 - 1) * 3 + 1] * Rhs[1];
            Xnorm0 += D[(i_0 - 1) * 3 + 2] * Rhs[2];
            i_0 = static_cast<int16_T>(i) + 3;
            if (i_0 > 32767) {
              i_0 = 32767;
            }

            Opt[i_0 - 1] = Xnorm0;
            for (kNext_0 = 1; kNext_0 - 1 < nA; kNext_0++) {
              int32_T tmp_0;
              int32_T tmp_1;
              i_0 = static_cast<int16_T>(i) + 3;
              if (i_0 > 32767) {
                i_0 = 32767;
              }

              tmp_0 = static_cast<int16_T>(i) + 3;
              if (tmp_0 > 32767) {
                tmp_0 = 32767;
              }

              tmp_1 = static_cast<int16_T>(kNext_0) + 3;
              if (tmp_1 > 32767) {
                tmp_1 = 32767;
              }

              Opt[i_0 - 1] = U[((static_cast<int16_T>(kNext_0) - 1) * 3 +
                                static_cast<int16_T>(i)) - 1] * Rhs[tmp_1 - 1] +
                Opt[tmp_0 - 1];
            }
          }

          Xnorm0 = -1.0E-12;
          kDrop = 0;
          for (i = 1; i - 1 < nA; i++) {
            i_0 = static_cast<int16_T>(i) + 3;
            if (i_0 > 32767) {
              i_0 = 32767;
            }

            lambda[iC[static_cast<int16_T>(i) - 1] - 1] = Opt[i_0 - 1];
            i_0 = static_cast<int16_T>(i) + 3;
            if (i_0 > 32767) {
              i_0 = 32767;
            }

            if ((Opt[i_0 - 1] < Xnorm0) && (static_cast<int16_T>(i) <= nA)) {
              kDrop = static_cast<int16_T>(i);
              i_0 = static_cast<int16_T>(i) + 3;
              if (i_0 > 32767) {
                i_0 = 32767;
              }

              Xnorm0 = Opt[i_0 - 1];
            }
          }

          if (kDrop <= 0) {
            DualFeasible = true;
            x[0] = Opt[0];
            x[1] = Opt[1];
            x[2] = Opt[2];
          } else {
            (*status)++;
            if (static_cast<int32_T>(*status) > q) {
              nA = 0;
              std::memset(&iA[0], 0, 34U * sizeof(int16_T));
              std::memset(&iC[0], 0, 34U * sizeof(int16_T));
              ColdReset = true;
            } else {
              lambda[iC[kDrop - 1] - 1] = 0.0;
              speedgoat_target_DropConstraint(kDrop, iA, &nA, iC);
            }
          }
        }
      } else {
        if (nA <= 0) {
          std::memset(&lambda[0], 0, 34U * sizeof(real_T));
          speedgoat_target__Unconstrained(b_Hinv, f, x, 3);
        }

        exitg3 = 1;
      }
    } while (exitg3 == 0);

    if (exitg3 == 1) {
      guard1 = true;
    }
  } else {
    speedgoat_target__Unconstrained(b_Hinv, f, x, 3);
    guard1 = true;
  }

  if (guard1) {
    boolean_T exitg2;
    Xnorm0 = speedgoat_target_model_202_norm(x);
    exitg2 = false;
    while ((!exitg2) && (static_cast<int32_T>(*status) <= maxiter)) {
      real_T b_Ac_1;
      real_T cVal;
      real_T t;
      cMin = -FeasTol;
      kNext = 0;
      for (i = 0; i < 34; i++) {
        t = cTol[i];
        if (!cTolComputed) {
          i_0 = i + 1;
          b_Ac_0[0] = b_Ac[i_0 - 1] * x[0];
          b_Ac_0[1] = b_Ac[i_0 + 33] * x[1];
          b_Ac_0[2] = b_Ac[i_0 + 67] * x[2];
          speedgoat_target_model_2021_abs(b_Ac_0, z);
          cVal = speedgoat_target_model__maximum(z);
          if ((t >= cVal) || rtIsNaN(cVal)) {
            cVal = t;
          }

          t = cVal;
        }

        if (iA[i] == 0) {
          i_0 = i + 1;
          b_Ac_1 = b_Ac[i_0 - 1] * x[0];
          b_Ac_1 += b_Ac[i_0 + 33] * x[1];
          b_Ac_1 += b_Ac[i_0 + 67] * x[2];
          cVal = (b_Ac_1 - b[i]) / t;
          if (cVal < cMin) {
            cMin = cVal;
            kNext = static_cast<int16_T>(i + 1);
          }
        }

        cTol[i] = t;
      }

      cTolComputed = true;
      if (kNext <= 0) {
        exitg2 = true;
      } else if (static_cast<int32_T>(*status) == maxiter) {
        *status = 0.0;
        exitg2 = true;
      } else {
        int32_T exitg1;
        do {
          exitg1 = 0;
          if ((kNext > 0) && (static_cast<int32_T>(*status) <= maxiter)) {
            boolean_T guard2 = false;
            guard2 = false;
            if (nA == 0) {
              kNext_0 = kNext;
              for (i_0 = 0; i_0 < 3; i_0++) {
                cMin = b_Ac[kNext_0 - 1] * b_Hinv[i_0];
                cMin += b_Hinv[i_0 + 3] * b_Ac[kNext_0 + 33];
                cMin += b_Hinv[i_0 + 6] * b_Ac[kNext_0 + 67];
                z[i_0] = cMin;
              }

              guard2 = true;
            } else {
              speedgoat_target_mod_KWIKfactor(b_Ac, iC, nA, b_Linv, D, b_H, 3,
                RLinv, &cMin);
              if (cMin <= 0.0) {
                *status = -2.0;
                exitg1 = 1;
              } else {
                kNext_0 = kNext;
                for (i_0 = 0; i_0 < 9; i_0++) {
                  U[i_0] = -b_H[i_0];
                }

                for (i_0 = 0; i_0 < 3; i_0++) {
                  cMin = b_Ac[kNext_0 - 1] * U[i_0];
                  cMin += U[i_0 + 3] * b_Ac[kNext_0 + 33];
                  cMin += U[i_0 + 6] * b_Ac[kNext_0 + 67];
                  z[i_0] = cMin;
                }

                for (i = 1; i - 1 < nA; i++) {
                  kNext_0 = kNext;
                  i_0 = static_cast<int16_T>(i);
                  b_Ac_1 = D[(i_0 - 1) * 3] * b_Ac[kNext_0 - 1];
                  b_Ac_1 += D[(i_0 - 1) * 3 + 1] * b_Ac[kNext_0 + 33];
                  b_Ac_1 += D[(i_0 - 1) * 3 + 2] * b_Ac[kNext_0 + 67];
                  r[static_cast<int16_T>(i) - 1] = b_Ac_1;
                }

                guard2 = true;
              }
            }

            if (guard2) {
              kDrop = 0;
              cMin = 0.0;
              DualFeasible = true;
              ColdReset = true;
              if (nA > 0) {
                boolean_T exitg4;
                i = 0;
                exitg4 = false;
                while ((!exitg4) && (i <= nA - 1)) {
                  if (r[i] >= 1.0E-12) {
                    ColdReset = false;
                    exitg4 = true;
                  } else {
                    i++;
                  }
                }
              }

              ColdReset = ((nA == 0) || ColdReset);
              if (!ColdReset) {
                for (i = 1; i - 1 < nA; i++) {
                  if (r[static_cast<int16_T>(i) - 1] > 1.0E-12) {
                    cVal = lambda[iC[static_cast<int16_T>(i) - 1] - 1] / r[
                      static_cast<int16_T>(i) - 1];
                    if ((kDrop == 0) || (cVal < rMin)) {
                      rMin = cVal;
                      kDrop = static_cast<int16_T>(i);
                    }
                  }
                }

                if (kDrop > 0) {
                  cMin = rMin;
                  DualFeasible = false;
                }
              }

              kNext_0 = kNext;
              b_Ac_0[0] = b_Ac[kNext_0 - 1];
              b_Ac_0[1] = b_Ac[kNext_0 + 33];
              b_Ac_0[2] = b_Ac[kNext_0 + 67];
              cVal = speedgoat_target_model_2_mtimes(z, b_Ac_0);
              if (cVal <= 0.0) {
                cVal = 0.0;
                ColdReset = true;
              } else {
                kNext_0 = kNext;
                b_Ac_1 = b_Ac[kNext_0 - 1] * x[0];
                b_Ac_1 += b_Ac[kNext_0 + 33] * x[1];
                b_Ac_1 += b_Ac[kNext_0 + 67] * x[2];
                cVal = (b[kNext - 1] - b_Ac_1) / cVal;
                ColdReset = false;
              }

              if (DualFeasible && ColdReset) {
                *status = -1.0;
                exitg1 = 1;
              } else {
                if (ColdReset) {
                  t = cMin;
                } else if (DualFeasible) {
                  t = cVal;
                } else if ((cMin <= cVal) || rtIsNaN(cVal)) {
                  t = cMin;
                } else {
                  t = cVal;
                }

                for (i = 1; i - 1 < nA; i++) {
                  lambda[iC[static_cast<int16_T>(i) - 1] - 1] -= r
                    [static_cast<int16_T>(i) - 1] * t;
                  if ((iC[static_cast<int16_T>(i) - 1] <= 34) && (lambda[iC[
                       static_cast<int16_T>(i) - 1] - 1] < 0.0)) {
                    lambda[iC[static_cast<int16_T>(i) - 1] - 1] = 0.0;
                  }
                }

                lambda[kNext - 1] += t;
                if (t == cMin) {
                  speedgoat_target_DropConstraint(kDrop, iA, &nA, iC);
                }

                if (!ColdReset) {
                  cMin = x[0];
                  cMin += t * z[0];
                  x[0] = cMin;
                  cMin = x[1];
                  cMin += t * z[1];
                  x[1] = cMin;
                  cMin = x[2];
                  cMin += t * z[2];
                  x[2] = cMin;
                  if (t == cVal) {
                    if (nA == 3) {
                      *status = -1.0;
                      exitg1 = 1;
                    } else {
                      i_0 = nA + 1;
                      if (i_0 > 32767) {
                        i_0 = 32767;
                      }

                      nA = static_cast<int16_T>(i_0);
                      iC[nA - 1] = kNext;
                      kDrop = nA;
                      while ((kDrop > 1) && (!(iC[kDrop - 1] > iC[kDrop - 2])))
                      {
                        int16_T iSave;
                        iSave = iC[kDrop - 1];
                        iC[kDrop - 1] = iC[kDrop - 2];
                        iC[kDrop - 2] = iSave;
                        kDrop = static_cast<int16_T>(kDrop - 1);
                      }

                      iA[kNext - 1] = 1;
                      kNext = 0;
                      (*status)++;
                    }
                  } else {
                    (*status)++;
                  }
                } else {
                  (*status)++;
                }
              }
            }
          } else {
            cMin = speedgoat_target_model_202_norm(x);
            if (std::abs(cMin - Xnorm0) > 0.001) {
              Xnorm0 = cMin;
              speedgoat_target_model_20_abs_l(b, tmp);
              speedgoat_target_model_maximum2(tmp, 1.0, cTol);
              cTolComputed = false;
            }

            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    }
  }
}

/* Model step function for TID0 */
void speedgoat_target_model_2021b_step0(void) /* Sample time: [0.1s, 0.0s] */
{
  real_T x[4];
  real_T poseF[3];
  real_T refPose[3];
  real_T tmp[2];
  real_T absx;
  real_T d_idx_0;
  real_T d_idx_1;
  real_T tHat_idx_0;
  real_T tHat_idx_1;
  int32_T s172_iter;
  int8_T wrBufIdx;

  /* Update the flag to indicate when data transfers from
   *  Sample time: [0.1s, 0.0s] to Sample time: [1.0s, 0.0s]  */
  (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1)++;
  if ((speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1) > 9) {
    speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 = 0;
  }

  /* Reset subsysRan breadcrumbs */
  srClearBC(speedgoat_target_model_2021b_DW.Subsystem2_SubsysRanBC);

  /* Reset subsysRan breadcrumbs */
  srClearBC(speedgoat_target_model_2021b_DW.Subsystem1_SubsysRanBC);

  /* Reset subsysRan breadcrumbs */
  srClearBC(speedgoat_target_model_2021b_DW.ByteUnpackFifthset_SubsysRanBC);

  /* Reset subsysRan breadcrumbs */
  srClearBC(speedgoat_target_model_2021b_DW.ByteUnpackFirstset_SubsysRanBC);

  /* Reset subsysRan breadcrumbs */
  srClearBC(speedgoat_target_model_2021b_DW.ByteUnpackFourthset_SubsysRanBC);

  /* Reset subsysRan breadcrumbs */
  srClearBC(speedgoat_target_model_2021b_DW.ByteUnpackSecondset_SubsysRanBC);

  /* Reset subsysRan breadcrumbs */
  srClearBC(speedgoat_target_model_2021b_DW.ByteUnpackThirdset_SubsysRanBC);

  /* Outputs for Iterator SubSystem: '<S2>/UDP Receive: Fourth set' incorporates:
   *  WhileIterator: '<S172>/While Iterator'
   */
  s172_iter = 1;
  do {
    speedgoat_target_model_2021b_B.WhileIterator_c = s172_iter;

    {
      try {
        slrealtime::ip::udp::Socket *udpSock = reinterpret_cast<slrealtime::ip::
          udp::Socket*>(speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_f[0]);
        char *buffer = (char *)
          speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_f[1];
        memset(buffer,0,24);
        void *dataPort = &speedgoat_target_model_2021b_B.UDPReceive2_o1_d[0];
        int_T numBytesAvail = (int_T)(udpSock->bytesToRead());
        if (numBytesAvail > 0) {
          uint8_t* fmAddArg = (uint8_t *)
            speedgoat_target_model_2021_cal->UDPReceive2_fmAddress_c;
          size_t num_bytesRcvd = (size_t)(udpSock->receive(buffer,
            ( numBytesAvail<65504 )? numBytesAvail:65504, !0,fmAddArg));
          if (num_bytesRcvd == 0) {
            speedgoat_target_model_2021b_B.UDPReceive2_o2_h = 0;
          } else {
            speedgoat_target_model_2021b_B.UDPReceive2_o2_h = (double)
              num_bytesRcvd;
            memcpy(dataPort,(void*)buffer,24);
          }
        } else {
          speedgoat_target_model_2021b_B.UDPReceive2_o2_h = 0;
        }
      } catch (std::exception& e) {
        std::string tmp = std::string(e.what());
        static std::string eMsg = tmp;
        rtmSetErrorStatus(speedgoat_target_model_2021b_M, eMsg.c_str());
        rtmSetStopRequested(speedgoat_target_model_2021b_M, 1);
        ;
      }
    }

    speedgoat_target_model_2021b_B.Compare_n =
      (speedgoat_target_model_2021b_B.UDPReceive2_o2_h ==
       speedgoat_target_model_2021_cal->CompareToConstant_const_l);
    speedgoat_target_model_2021b_B.Subtract_d = static_cast<real_T>
      (speedgoat_target_model_2021b_B.WhileIterator_c) -
      speedgoat_target_model_2021_cal->Constant_Value_p;
    s172_iter++;
  } while (speedgoat_target_model_2021b_B.Compare_n && (s172_iter <= 5));

  /* End of Outputs for SubSystem: '<S2>/UDP Receive: Fourth set' */

  /* DataTypeConversion: '<S2>/Data Type Conversion4' */
  speedgoat_target_model_2021b_B.DataTypeConversion4 =
    (speedgoat_target_model_2021b_B.Subtract_d != 0.0);

  /* Outputs for Enabled SubSystem: '<S2>/Byte Unpack: Fourth set' incorporates:
   *  EnablePort: '<S167>/Enable'
   */
  if (speedgoat_target_model_2021b_B.DataTypeConversion4) {
    /* S-Function (slrealtimebytepacking): '<S167>/Byte Unpacking ' */

    /* Byte Unpacking: <S167>/Byte Unpacking  */
    (void)memcpy((uint8_T*)&speedgoat_target_model_2021b_B.RefPoses[0], (uint8_T*)
                 &speedgoat_target_model_2021b_B.UDPReceive2_o1_d[0] + 0, 24);
    srUpdateBC(speedgoat_target_model_2021b_DW.ByteUnpackFourthset_SubsysRanBC);
  }

  /* End of Outputs for SubSystem: '<S2>/Byte Unpack: Fourth set' */

  /* Outputs for Iterator SubSystem: '<S2>/UDP Receive: Fifth set' incorporates:
   *  WhileIterator: '<S170>/While Iterator'
   */
  s172_iter = 1;
  do {
    speedgoat_target_model_2021b_B.WhileIterator_e = s172_iter;

    {
      try {
        slrealtime::ip::udp::Socket *udpSock = reinterpret_cast<slrealtime::ip::
          udp::Socket*>(speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_e0[0]);
        char *buffer = (char *)
          speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_e0[1];
        memset(buffer,0,24);
        void *dataPort = &speedgoat_target_model_2021b_B.UDPReceive2_o1_m[0];
        int_T numBytesAvail = (int_T)(udpSock->bytesToRead());
        if (numBytesAvail > 0) {
          uint8_t* fmAddArg = (uint8_t *)
            speedgoat_target_model_2021_cal->UDPReceive2_fmAddress;
          size_t num_bytesRcvd = (size_t)(udpSock->receive(buffer,
            ( numBytesAvail<65504 )? numBytesAvail:65504, !0,fmAddArg));
          if (num_bytesRcvd == 0) {
            speedgoat_target_model_2021b_B.UDPReceive2_o2_c = 0;
          } else {
            speedgoat_target_model_2021b_B.UDPReceive2_o2_c = (double)
              num_bytesRcvd;
            memcpy(dataPort,(void*)buffer,24);
          }
        } else {
          speedgoat_target_model_2021b_B.UDPReceive2_o2_c = 0;
        }
      } catch (std::exception& e) {
        std::string tmp = std::string(e.what());
        static std::string eMsg = tmp;
        rtmSetErrorStatus(speedgoat_target_model_2021b_M, eMsg.c_str());
        rtmSetStopRequested(speedgoat_target_model_2021b_M, 1);
        ;
      }
    }

    speedgoat_target_model_2021b_B.Compare_b =
      (speedgoat_target_model_2021b_B.UDPReceive2_o2_c ==
       speedgoat_target_model_2021_cal->CompareToConstant_const);
    speedgoat_target_model_2021b_B.Subtract_bo = static_cast<real_T>
      (speedgoat_target_model_2021b_B.WhileIterator_e) -
      speedgoat_target_model_2021_cal->Constant_Value_k;
    s172_iter++;
  } while (speedgoat_target_model_2021b_B.Compare_b && (s172_iter <= 5));

  /* End of Outputs for SubSystem: '<S2>/UDP Receive: Fifth set' */

  /* DataTypeConversion: '<S2>/Data Type Conversion5' */
  speedgoat_target_model_2021b_B.DataTypeConversion5 =
    (speedgoat_target_model_2021b_B.Subtract_bo != 0.0);

  /* Outputs for Enabled SubSystem: '<S2>/Byte Unpack: Fifth set' incorporates:
   *  EnablePort: '<S165>/Enable'
   */
  if (speedgoat_target_model_2021b_B.DataTypeConversion5) {
    /* S-Function (slrealtimebytepacking): '<S165>/Byte Unpacking 1' */

    /* Byte Unpacking: <S165>/Byte Unpacking 1 */
    (void)memcpy((uint8_T*)&speedgoat_target_model_2021b_B.CurrPoses[0],
                 (uint8_T*)&speedgoat_target_model_2021b_B.UDPReceive2_o1_m[0] +
                 0, 24);
    srUpdateBC(speedgoat_target_model_2021b_DW.ByteUnpackFifthset_SubsysRanBC);
  }

  /* End of Outputs for SubSystem: '<S2>/Byte Unpack: Fifth set' */

  /* Outputs for Iterator SubSystem: '<S2>/UDP Receive: First set' incorporates:
   *  WhileIterator: '<S171>/While Iterator'
   */
  s172_iter = 1;
  do {
    speedgoat_target_model_2021b_B.WhileIterator_p = s172_iter;

    {
      try {
        slrealtime::ip::udp::Socket *udpSock = reinterpret_cast<slrealtime::ip::
          udp::Socket*>(speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_e[0]);
        char *buffer = (char *)
          speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_e[1];
        memset(buffer,0,32);
        void *dataPort = &speedgoat_target_model_2021b_B.UDPReceive2_o1_p[0];
        int_T numBytesAvail = (int_T)(udpSock->bytesToRead());
        if (numBytesAvail > 0) {
          uint8_t* fmAddArg = (uint8_t *)
            speedgoat_target_model_2021_cal->UDPReceive2_fmAddress_e;
          size_t num_bytesRcvd = (size_t)(udpSock->receive(buffer,
            ( numBytesAvail<65504 )? numBytesAvail:65504, !0,fmAddArg));
          if (num_bytesRcvd == 0) {
            speedgoat_target_model_2021b_B.UDPReceive2_o2_j = 0;
          } else {
            speedgoat_target_model_2021b_B.UDPReceive2_o2_j = (double)
              num_bytesRcvd;
            memcpy(dataPort,(void*)buffer,32);
          }
        } else {
          speedgoat_target_model_2021b_B.UDPReceive2_o2_j = 0;
        }
      } catch (std::exception& e) {
        std::string tmp = std::string(e.what());
        static std::string eMsg = tmp;
        rtmSetErrorStatus(speedgoat_target_model_2021b_M, eMsg.c_str());
        rtmSetStopRequested(speedgoat_target_model_2021b_M, 1);
        ;
      }
    }

    speedgoat_target_model_2021b_B.Compare_p =
      (speedgoat_target_model_2021b_B.UDPReceive2_o2_j ==
       speedgoat_target_model_2021_cal->CompareToConstant_const_e);
    speedgoat_target_model_2021b_B.Subtract_b = static_cast<real_T>
      (speedgoat_target_model_2021b_B.WhileIterator_p) -
      speedgoat_target_model_2021_cal->Constant_Value_dy;
    s172_iter++;
  } while (speedgoat_target_model_2021b_B.Compare_p && (s172_iter <= 5));

  /* End of Outputs for SubSystem: '<S2>/UDP Receive: First set' */

  /* DataTypeConversion: '<S2>/Data Type Conversion2' */
  speedgoat_target_model_2021b_B.DataTypeConversion2 =
    (speedgoat_target_model_2021b_B.Subtract_b != 0.0);

  /* Outputs for Enabled SubSystem: '<S2>/Byte Unpack: First set' incorporates:
   *  EnablePort: '<S166>/Enable'
   */
  if (speedgoat_target_model_2021b_B.DataTypeConversion2) {
    /* S-Function (slrealtimebytepacking): '<S166>/Byte Unpacking ' */

    /* Byte Unpacking: <S166>/Byte Unpacking  */
    (void)memcpy((uint8_T*)&speedgoat_target_model_2021b_B.ObjAhead, (uint8_T*)
                 &speedgoat_target_model_2021b_B.UDPReceive2_o1_p[0] + 0, 8);
    (void)memcpy((uint8_T*)&speedgoat_target_model_2021b_B.Longitudinalvelocity,
                 (uint8_T*)&speedgoat_target_model_2021b_B.UDPReceive2_o1_p[0] +
                 8, 8);
    (void)memcpy((uint8_T*)&speedgoat_target_model_2021b_B.SceneID, (uint8_T*)
                 &speedgoat_target_model_2021b_B.UDPReceive2_o1_p[0] + 16, 8);
    (void)memcpy((uint8_T*)&speedgoat_target_model_2021b_B.Relativevelocity,
                 (uint8_T*)&speedgoat_target_model_2021b_B.UDPReceive2_o1_p[0] +
                 24, 8);
    srUpdateBC(speedgoat_target_model_2021b_DW.ByteUnpackFirstset_SubsysRanBC);
  }

  /* End of Outputs for SubSystem: '<S2>/Byte Unpack: First set' */

  /* Outputs for Iterator SubSystem: '<S2>/UDP Receive: Third set' incorporates:
   *  WhileIterator: '<S174>/While Iterator'
   */
  s172_iter = 1;
  do {
    speedgoat_target_model_2021b_B.WhileIterator = s172_iter;

    {
      try {
        slrealtime::ip::udp::Socket *udpSock = reinterpret_cast<slrealtime::ip::
          udp::Socket*>(speedgoat_target_model_2021b_DW.UDPReceive2_PWORK[0]);
        char *buffer = (char *)
          speedgoat_target_model_2021b_DW.UDPReceive2_PWORK[1];
        memset(buffer,0,24);
        void *dataPort = &speedgoat_target_model_2021b_B.UDPReceive2_o1[0];
        int_T numBytesAvail = (int_T)(udpSock->bytesToRead());
        if (numBytesAvail > 0) {
          uint8_t* fmAddArg = (uint8_t *)
            speedgoat_target_model_2021_cal->UDPReceive2_fmAddress_i;
          size_t num_bytesRcvd = (size_t)(udpSock->receive(buffer,
            ( numBytesAvail<65504 )? numBytesAvail:65504, !0,fmAddArg));
          if (num_bytesRcvd == 0) {
            speedgoat_target_model_2021b_B.UDPReceive2_o2 = 0;
          } else {
            speedgoat_target_model_2021b_B.UDPReceive2_o2 = (double)
              num_bytesRcvd;
            memcpy(dataPort,(void*)buffer,24);
          }
        } else {
          speedgoat_target_model_2021b_B.UDPReceive2_o2 = 0;
        }
      } catch (std::exception& e) {
        std::string tmp = std::string(e.what());
        static std::string eMsg = tmp;
        rtmSetErrorStatus(speedgoat_target_model_2021b_M, eMsg.c_str());
        rtmSetStopRequested(speedgoat_target_model_2021b_M, 1);
        ;
      }
    }

    speedgoat_target_model_2021b_B.Compare =
      (speedgoat_target_model_2021b_B.UDPReceive2_o2 ==
       speedgoat_target_model_2021_cal->CompareToConstant_const_g);
    speedgoat_target_model_2021b_B.Subtract = static_cast<real_T>
      (speedgoat_target_model_2021b_B.WhileIterator) -
      speedgoat_target_model_2021_cal->Constant_Value_j;
    s172_iter++;
  } while (speedgoat_target_model_2021b_B.Compare && (s172_iter <= 5));

  /* End of Outputs for SubSystem: '<S2>/UDP Receive: Third set' */

  /* DataTypeConversion: '<S2>/Data Type Conversion3' */
  speedgoat_target_model_2021b_B.DataTypeConversion3 =
    (speedgoat_target_model_2021b_B.Subtract != 0.0);

  /* Outputs for Enabled SubSystem: '<S2>/Byte Unpack: Third set' incorporates:
   *  EnablePort: '<S169>/Enable'
   */
  if (speedgoat_target_model_2021b_B.DataTypeConversion3) {
    /* S-Function (slrealtimebytepacking): '<S169>/Byte Unpacking ' */

    /* Byte Unpacking: <S169>/Byte Unpacking  */
    (void)memcpy((uint8_T*)&speedgoat_target_model_2021b_B.Direction, (uint8_T*)
                 &speedgoat_target_model_2021b_B.UDPReceive2_o1[0] + 0, 8);
    (void)memcpy((uint8_T*)&speedgoat_target_model_2021b_B.Curvature, (uint8_T*)
                 &speedgoat_target_model_2021b_B.UDPReceive2_o1[0] + 8, 8);
    (void)memcpy((uint8_T*)&speedgoat_target_model_2021b_B.yawrate, (uint8_T*)
                 &speedgoat_target_model_2021b_B.UDPReceive2_o1[0] + 16, 8);
    srUpdateBC(speedgoat_target_model_2021b_DW.ByteUnpackThirdset_SubsysRanBC);
  }

  /* End of Outputs for SubSystem: '<S2>/Byte Unpack: Third set' */

  /* MATLAB Function: '<S7>/Kinematic' incorporates:
   *  S-Function (slrealtimebytepacking): '<S165>/Byte Unpacking 1'
   *  S-Function (slrealtimebytepacking): '<S167>/Byte Unpacking '
   */
  refPose[0] = speedgoat_target_model_2021b_B.RefPoses[0];
  refPose[1] = speedgoat_target_model_2021b_B.RefPoses[1];
  poseF[0] = speedgoat_target_model_2021b_B.CurrPoses[0];
  poseF[1] = speedgoat_target_model_2021b_B.CurrPoses[1];
  if (speedgoat_target_model_2021b_B.Direction == 1.0) {
    absx = speedgoat_target_model_2021_cal->LateralControllerStanley_Positi;
  } else {
    absx = speedgoat_target_model_2021_cal->LateralControllerStanley_Posi_i;
  }

  refPose[2] = 0.017453292519943295 * speedgoat_target_model_2021b_B.RefPoses[2];
  speedg_angleUtilities_wrapTo2Pi(&refPose[2]);
  poseF[2] = 0.017453292519943295 * speedgoat_target_model_2021b_B.CurrPoses[2];
  speedg_angleUtilities_wrapTo2Pi(&poseF[2]);
  tHat_idx_0 = std::cos(refPose[2]);
  tHat_idx_1 = std::sin(refPose[2]);
  if (speedgoat_target_model_2021b_B.Direction == 1.0) {
    poseF[0] = speedgoat_target_model_2021_cal->Kinematic_Wheelbase * std::cos
      (poseF[2]) + speedgoat_target_model_2021b_B.CurrPoses[0];
    poseF[1] = speedgoat_target_model_2021_cal->Kinematic_Wheelbase * std::sin
      (poseF[2]) + speedgoat_target_model_2021b_B.CurrPoses[1];
    d_idx_0 = poseF[0] - refPose[0];
    d_idx_1 = poseF[1] - refPose[1];
  } else {
    d_idx_0 = poseF[0] - refPose[0];
    d_idx_1 = poseF[1] - refPose[1];
  }

  tHat_idx_0 = -(d_idx_0 * tHat_idx_1 - tHat_idx_0 * d_idx_1);
  tHat_idx_1 = (poseF[2] - refPose[2]) + 3.1415926535897931;
  speedg_angleUtilities_wrapTo2Pi(&tHat_idx_1);
  if (speedgoat_target_model_2021b_B.Direction == 1.0) {
    absx = -(std::atan(absx * tHat_idx_0 /
                       (speedgoat_target_model_2021b_B.Longitudinalvelocity +
                        1.0)) + (tHat_idx_1 - 3.1415926535897931));
  } else {
    absx = std::atan(absx * tHat_idx_0 /
                     (speedgoat_target_model_2021b_B.Longitudinalvelocity + -1.0))
      + (tHat_idx_1 - 3.1415926535897931);
  }

  speedgoat_target_model_2021b_B.steerCmd = 57.295779513082323 * absx;
  absx = speedgoat_target_model_2021b_B.steerCmd;
  if (absx < 0.0) {
    tHat_idx_1 = -1.0;
  } else if (absx > 0.0) {
    tHat_idx_1 = 1.0;
  } else if (absx == 0.0) {
    tHat_idx_1 = 0.0;
  } else {
    tHat_idx_1 = (rtNaN);
  }

  absx = std::abs(speedgoat_target_model_2021b_B.steerCmd);
  tHat_idx_0 = speedgoat_target_model_2021_cal->Kinematic_MaxSteeringAngle;
  if ((absx <= tHat_idx_0) || rtIsNaN(tHat_idx_0)) {
    tHat_idx_0 = absx;
  }

  speedgoat_target_model_2021b_B.steerCmd = tHat_idx_1 * tHat_idx_0;

  /* End of MATLAB Function: '<S7>/Kinematic' */

  /* Product: '<S47>/Product' */
  speedgoat_target_model_2021b_B.Product =
    speedgoat_target_model_2021b_B.Longitudinalvelocity *
    speedgoat_target_model_2021b_B.Curvature;

  /* Gain: '<S48>/Gain' */
  speedgoat_target_model_2021b_B.Gain =
    speedgoat_target_model_2021_cal->Gain_Gain_g *
    speedgoat_target_model_2021b_B.Product;

  /* Product: '<S47>/Multiply' */
  speedgoat_target_model_2021b_B.Multiply =
    speedgoat_target_model_2021b_B.Longitudinalvelocity *
    speedgoat_target_model_2021b_B.Gain;

  /* Gain: '<S47>/Gain1' */
  speedgoat_target_model_2021b_B.CurvedSteadyYaw =
    speedgoat_target_model_2021_cal->Gain1_Gain_j3 *
    speedgoat_target_model_2021b_B.Multiply;

  /* Gain: '<S1>/Gain' */
  speedgoat_target_model_2021b_B.Gain_d =
    speedgoat_target_model_2021_cal->Gain_Gain_k *
    speedgoat_target_model_2021b_B.yawrate;

  /* Sum: '<S47>/Minus' */
  speedgoat_target_model_2021b_B.Minus = speedgoat_target_model_2021b_B.Gain -
    speedgoat_target_model_2021b_B.Gain_d;

  /* Gain: '<S47>/Gain' */
  speedgoat_target_model_2021b_B.YawRateFeedback =
    speedgoat_target_model_2021_cal->LateralControllerStanley_YawRat *
    speedgoat_target_model_2021b_B.Minus;

  /* UnitDelay: '<S47>/Unit Delay' */
  speedgoat_target_model_2021b_B.UnitDelay =
    speedgoat_target_model_2021b_DW.UnitDelay_DSTATE;

  /* UnitDelay: '<S1>/Unit Delay' */
  speedgoat_target_model_2021b_B.UnitDelay_n =
    speedgoat_target_model_2021b_DW.UnitDelay_DSTATE_o;

  /* Sum: '<S47>/Minus1' */
  speedgoat_target_model_2021b_B.Minus1 =
    speedgoat_target_model_2021b_B.UnitDelay -
    speedgoat_target_model_2021b_B.UnitDelay_n;

  /* Gain: '<S47>/Gain2' */
  speedgoat_target_model_2021b_B.SteeringDelay =
    speedgoat_target_model_2021_cal->LateralControllerStanley_DelayG *
    speedgoat_target_model_2021b_B.Minus1;

  /* Sum: '<S47>/Add' */
  speedgoat_target_model_2021b_B.Add = ((speedgoat_target_model_2021b_B.steerCmd
    + speedgoat_target_model_2021b_B.CurvedSteadyYaw) +
    speedgoat_target_model_2021b_B.YawRateFeedback) +
    speedgoat_target_model_2021b_B.SteeringDelay;

  /* Saturate: '<S47>/Saturation' */
  absx = speedgoat_target_model_2021b_B.Add;
  tHat_idx_0 = speedgoat_target_model_2021_cal->Saturation_LowerSat;
  tHat_idx_1 = speedgoat_target_model_2021_cal->Saturation_UpperSat;
  if (absx > tHat_idx_1) {
    /* Saturate: '<S47>/Saturation' */
    speedgoat_target_model_2021b_B.Saturation = tHat_idx_1;
  } else if (absx < tHat_idx_0) {
    /* Saturate: '<S47>/Saturation' */
    speedgoat_target_model_2021b_B.Saturation = tHat_idx_0;
  } else {
    /* Saturate: '<S47>/Saturation' */
    speedgoat_target_model_2021b_B.Saturation = absx;
  }

  /* End of Saturate: '<S47>/Saturation' */

  /* Outputs for Iterator SubSystem: '<S2>/UDP Receive: Second set' incorporates:
   *  WhileIterator: '<S173>/While Iterator'
   */
  s172_iter = 1;
  do {
    speedgoat_target_model_2021b_B.WhileIterator_k = s172_iter;

    {
      try {
        slrealtime::ip::udp::Socket *udpSock = reinterpret_cast<slrealtime::ip::
          udp::Socket*>(speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_b[0]);
        char *buffer = (char *)
          speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_b[1];
        memset(buffer,0,744);
        void *dataPort = &speedgoat_target_model_2021b_B.UDPReceive2_o1_o[0];
        int_T numBytesAvail = (int_T)(udpSock->bytesToRead());
        if (numBytesAvail > 0) {
          uint8_t* fmAddArg = (uint8_t *)
            speedgoat_target_model_2021_cal->UDPReceive2_fmAddress_h;
          size_t num_bytesRcvd = (size_t)(udpSock->receive(buffer,
            ( numBytesAvail<65504 )? numBytesAvail:65504, !0,fmAddArg));
          if (num_bytesRcvd == 0) {
            speedgoat_target_model_2021b_B.UDPReceive2_o2_o = 0;
          } else {
            speedgoat_target_model_2021b_B.UDPReceive2_o2_o = (double)
              num_bytesRcvd;
            memcpy(dataPort,(void*)buffer,744);
          }
        } else {
          speedgoat_target_model_2021b_B.UDPReceive2_o2_o = 0;
        }
      } catch (std::exception& e) {
        std::string tmp = std::string(e.what());
        static std::string eMsg = tmp;
        rtmSetErrorStatus(speedgoat_target_model_2021b_M, eMsg.c_str());
        rtmSetStopRequested(speedgoat_target_model_2021b_M, 1);
        ;
      }
    }

    speedgoat_target_model_2021b_B.Compare_j =
      (speedgoat_target_model_2021b_B.UDPReceive2_o2_o ==
       speedgoat_target_model_2021_cal->CompareToConstant_const_k);
    speedgoat_target_model_2021b_B.Subtract_p = static_cast<real_T>
      (speedgoat_target_model_2021b_B.WhileIterator_k) -
      speedgoat_target_model_2021_cal->Constant_Value_g;
    s172_iter++;
  } while (speedgoat_target_model_2021b_B.Compare_j && (s172_iter <= 5));

  /* End of Outputs for SubSystem: '<S2>/UDP Receive: Second set' */

  /* DataTypeConversion: '<S2>/Data Type Conversion1' */
  speedgoat_target_model_2021b_B.DataTypeConversion1 =
    (speedgoat_target_model_2021b_B.Subtract_p != 0.0);

  /* Outputs for Enabled SubSystem: '<S2>/Byte Unpack: Second set' incorporates:
   *  EnablePort: '<S168>/Enable'
   */
  if (speedgoat_target_model_2021b_B.DataTypeConversion1) {
    /* S-Function (slrealtimebytepacking): '<S168>/Byte Unpacking ' */

    /* Byte Unpacking: <S168>/Byte Unpacking  */
    (void)memcpy((uint8_T*)&speedgoat_target_model_2021b_B.ObjectDist, (uint8_T*)
                 &speedgoat_target_model_2021b_B.UDPReceive2_o1_o[0] + 0, 8);
    (void)memcpy((uint8_T*)&speedgoat_target_model_2021b_B.D2S, (uint8_T*)
                 &speedgoat_target_model_2021b_B.UDPReceive2_o1_o[0] + 8, 8);
    (void)memcpy((uint8_T*)&speedgoat_target_model_2021b_B.Relativedistance,
                 (uint8_T*)&speedgoat_target_model_2021b_B.UDPReceive2_o1_o[0] +
                 16, 8);
    (void)memcpy((uint8_T*)&speedgoat_target_model_2021b_B.TLstate[0], (uint8_T*)
                 &speedgoat_target_model_2021b_B.UDPReceive2_o1_o[0] + 24, 720);
    srUpdateBC(speedgoat_target_model_2021b_DW.ByteUnpackSecondset_SubsysRanBC);
  }

  /* End of Outputs for SubSystem: '<S2>/Byte Unpack: Second set' */

  /* Constant: '<S1>/Approaching Speed' */
  speedgoat_target_model_2021b_B.ApproachingSpeed =
    speedgoat_target_model_2021_cal->ApproachingSpeed_Value;

  /* Constant: '<S1>/Scene Speed' */
  speedgoat_target_model_2021b_B.SceneSpeed =
    speedgoat_target_model_2021_cal->SceneSpeed_Value;

  /* Constant: '<S1>/Constant1' */
  speedgoat_target_model_2021b_B.Constant1 =
    speedgoat_target_model_2021_cal->Constant1_Value_i4;

  /* Constant: '<S1>/RR Crossing' */
  speedgoat_target_model_2021b_B.RRCrossing =
    speedgoat_target_model_2021_cal->RRCrossing_Value;

  /* Constant: '<S1>/Navigating Flag' */
  speedgoat_target_model_2021b_B.NavigatingFlag =
    speedgoat_target_model_2021_cal->NavigatingFlag_Value;

  /* MATLAB Function: '<S1>/MATLAB Function' incorporates:
   *  Constant: '<S1>/Constant'
   */
  refPose[0] = speedgoat_target_model_2021b_B.CurrPoses[0];
  refPose[1] = -speedgoat_target_model_2021b_B.CurrPoses[1];
  speedgoat_target_model_2021b_B.D2D = 0.0;
  tHat_idx_0 = speedgoat_target_model_2021_cal->Constant_Value_kp[0] - refPose[0];
  absx = std::abs(tHat_idx_0);
  if (rtIsNaN(absx) || (absx > speedgoat_target_model_2021b_B.D2D)) {
    speedgoat_target_model_2021b_B.D2D = absx;
  }

  x[0] = tHat_idx_0;
  tHat_idx_0 = speedgoat_target_model_2021_cal->Constant_Value_kp[1] - refPose[0];
  absx = std::abs(tHat_idx_0);
  if (rtIsNaN(absx) || (absx > speedgoat_target_model_2021b_B.D2D)) {
    speedgoat_target_model_2021b_B.D2D = absx;
  }

  x[1] = tHat_idx_0;
  tHat_idx_0 = speedgoat_target_model_2021_cal->Constant_Value_kp[0] - refPose[1];
  absx = std::abs(tHat_idx_0);
  if (rtIsNaN(absx) || (absx > speedgoat_target_model_2021b_B.D2D)) {
    speedgoat_target_model_2021b_B.D2D = absx;
  }

  x[2] = tHat_idx_0;
  tHat_idx_0 = speedgoat_target_model_2021_cal->Constant_Value_kp[1] - refPose[1];
  absx = std::abs(tHat_idx_0);
  if (rtIsNaN(absx) || (absx > speedgoat_target_model_2021b_B.D2D)) {
    speedgoat_target_model_2021b_B.D2D = absx;
  }

  x[3] = tHat_idx_0;
  if ((!rtIsInf(speedgoat_target_model_2021b_B.D2D)) && (!rtIsNaN
       (speedgoat_target_model_2021b_B.D2D))) {
    speedgoat_target_model_2021_svd(x, tmp);
    speedgoat_target_model_2021b_B.D2D = tmp[0];
  }

  /* End of MATLAB Function: '<S1>/MATLAB Function' */

  /* Constant: '<S1>/Clear2Go' */
  speedgoat_target_model_2021b_B.Clear2Go =
    speedgoat_target_model_2021_cal->Clear2Go_Value;

  /* Constant: '<S1>/RoadSpeed' */
  speedgoat_target_model_2021b_B.RoadSpeed =
    speedgoat_target_model_2021_cal->RoadSpeed_Value;

  /* Chart: '<S1>/Chart' */
  if (speedgoat_target_model_2021b_DW.temporalCounter_i1_i < 31U) {
    speedgoat_target_model_2021b_DW.temporalCounter_i1_i = static_cast<uint8_T>
      (speedgoat_target_model_2021b_DW.temporalCounter_i1_i + 1U);
  }

  speedgoat_target_model_2021b_DW.sfEvent_o = -1;
  if (speedgoat_target_model_2021b_DW.is_active_c21_speedgoat_target_ == 0U) {
    speedgoat_target_model_2021b_DW.is_active_c21_speedgoat_target_ = 1U;
    speedgoat_target_model_2021b_DW.is_c21_speedgoat_target_model_2 =
      speedgoat_target__IN_Stationary;
  } else if (speedgoat_target_model_2021b_DW.is_c21_speedgoat_target_model_2 ==
             speedgoat_target__IN_Navigation) {
    speedgoat_target_mod_Navigation();

    /* case IN_Stationary: */
  } else if (speedgoat_target_model_2021b_B.NavigatingFlag == 1.0) {
    speedgoat_target_model_2021b_DW.is_c21_speedgoat_target_model_2 =
      speedgoat_target__IN_Navigation;
    speedgoat_target_model_2021b_DW.is_Navigation =
      speedgoat_target_m_IN_DestCheck;
    if (speedgoat_target_model_2021b_B.D2D < 2.0) {
      speedgoat_target_model_2021b_DW.DestReached = 1.0;
      speedgoat_target_model_2021b_B.StopSim = 1.0;
    } else {
      speedgoat_target_model_2021b_DW.DestReached = 0.0;
    }
  } else {
    speedgoat_target_model_2021b_B.VelCmd = 0.0;
    speedgoat_target_model_2021b_DW.CurrentVelCmd =
      speedgoat_target_model_2021b_B.VelCmd;
    speedgoat_target_model_2021b_DW.Navi = 0.0;
    speedgoat_target_model_2021b_B.Trans = 0.0;
    speedgoat_target_model_2021b_B.StopSim = 0.0;
  }

  /* End of Chart: '<S1>/Chart' */

  /* Constant: '<S1>/WaitTime' */
  speedgoat_target_model_2021b_B.WaitTime =
    speedgoat_target_model_2021_cal->WaitTime_Value;

  /* Constant: '<S1>/MinSpeed' */
  speedgoat_target_model_2021b_B.MinSpeed =
    speedgoat_target_model_2021_cal->MinSpeed_Value;

  /* Constant: '<S1>/LaneAvailable' */
  speedgoat_target_model_2021b_B.LaneAvailable =
    speedgoat_target_model_2021_cal->LaneAvailable_Value;

  /* Constant: '<S1>/LaneChanged' */
  speedgoat_target_model_2021b_B.LaneChanged =
    speedgoat_target_model_2021_cal->LaneChanged_Value;

  /* Chart: '<S1>/Collision Avoidance' */
  if (speedgoat_target_model_2021b_DW.temporalCounter_i1 < MAX_uint32_T) {
    speedgoat_target_model_2021b_DW.temporalCounter_i1++;
  }

  speedgoat_target_model_2021b_DW.sfEvent = -1;
  if (speedgoat_target_model_2021b_DW.is_active_c9_speedgoat_target_m == 0U) {
    speedgoat_target_model_2021b_DW.is_active_c9_speedgoat_target_m = 1U;
    speedgoat_target_model_2021b_DW.is_c9_speedgoat_target_model_20 =
      speedgoat_target_model_IN_Start;
    speedgoat_target_model_2021b_DW.Stuck = 0.0;
    speedgoat_target_model_2021b_B.ACC = 0.0;
    speedgoat_target_model_2021b_B.LaneChangeCmd = 0.0;
  } else {
    switch (speedgoat_target_model_2021b_DW.is_c9_speedgoat_target_model_20) {
     case speedgoat_target_mode_IN_Follow:
      {
        speedgoat_target_model_2021b_B.ACC = 1.0;
        if ((speedgoat_target_model_2021b_B.Longitudinalvelocity <=
             speedgoat_target_model_2021b_B.MinSpeed) &&
            (speedgoat_target_model_2021b_DW.temporalCounter_i1 >= static_cast<
             uint32_T>(std::ceil(speedgoat_target_model_2021b_B.WaitTime * 10.0))))
        {
          speedgoat_target_model_2021b_DW.Stuck = 1.0;
          speedgoat_target_model_2021b_DW.is_c9_speedgoat_target_model_20 =
            speedgoat_target_mode_IN_Follow;
          speedgoat_target_model_2021b_DW.temporalCounter_i1 = 0U;
          speedgoat_target_model_2021b_B.ACC = 1.0;
        } else {
          boolean_T out;
          out = ((speedgoat_target_model_2021b_B.LaneAvailable == 1.0) &&
                 (speedgoat_target_model_2021b_DW.Stuck == 1.0));
          if (out) {
            speedgoat_target_model_2021b_DW.is_c9_speedgoat_target_model_20 =
              speedgoat_target__IN_LaneChange;
            speedgoat_target_model_2021b_B.LaneChangeCmd = 1.0;
          } else if (speedgoat_target_model_2021b_B.ObjAhead == 0.0) {
            speedgoat_target_model_2021b_DW.is_c9_speedgoat_target_model_20 =
              speedgoat_target_model_IN_Start;
            speedgoat_target_model_2021b_DW.Stuck = 0.0;
            speedgoat_target_model_2021b_B.ACC = 0.0;
            speedgoat_target_model_2021b_B.LaneChangeCmd = 0.0;
          }
        }
      }
      break;

     case speedgoat_target__IN_LaneChange:
      speedgoat_target_model_2021b_B.LaneChangeCmd = 1.0;
      if (speedgoat_target_model_2021b_B.LaneChanged == 1.0) {
        speedgoat_target_model_2021b_DW.is_c9_speedgoat_target_model_20 =
          speedgoat_target_model_IN_Start;
        speedgoat_target_model_2021b_DW.Stuck = 0.0;
        speedgoat_target_model_2021b_B.ACC = 0.0;
        speedgoat_target_model_2021b_B.LaneChangeCmd = 0.0;
      }
      break;

     default:
      /* case IN_Start: */
      speedgoat_target_model_2021b_B.ACC = 0.0;
      speedgoat_target_model_2021b_B.LaneChangeCmd = 0.0;
      if (speedgoat_target_model_2021b_B.ObjAhead == 1.0) {
        speedgoat_target_model_2021b_DW.is_c9_speedgoat_target_model_20 =
          speedgoat_target_mode_IN_Follow;
        speedgoat_target_model_2021b_DW.temporalCounter_i1 = 0U;
        speedgoat_target_model_2021b_B.ACC = 1.0;
      }
      break;
    }
  }

  /* End of Chart: '<S1>/Collision Avoidance' */

  /* Logic: '<S1>/OR' */
  speedgoat_target_model_2021b_B.OR = ((speedgoat_target_model_2021b_B.ACC_i !=
    0.0) || (speedgoat_target_model_2021b_B.ACC != 0.0));

  /* RateTransition generated from: '<S11>/Subsystem2' */
  if (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 == 1) {
    speedgoat_target_model_2021b_DW.acceleration_RdBufIdx = static_cast<int8_T>
      (speedgoat_target_model_2021b_DW.acceleration_RdBufIdx == 0);
  }

  /* RateTransition generated from: '<S11>/Subsystem2' */
  speedgoat_target_model_2021b_B.acceleration =
    speedgoat_target_model_2021b_DW.acceleration_Buf[speedgoat_target_model_2021b_DW.acceleration_RdBufIdx];

  /* RateTransition generated from: '<S49>/Switch' */
  if (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 == 1) {
    speedgoat_target_model_2021b_DW.TmpRTBAtSwitchOutport1_RdBufIdx =
      static_cast<int8_T>
      (speedgoat_target_model_2021b_DW.TmpRTBAtSwitchOutport1_RdBufIdx == 0);
  }

  /* RateTransition generated from: '<S49>/Switch' */
  speedgoat_target_model_2021b_B.TmpRTBAtSwitchOutport1 =
    speedgoat_target_model_2021b_DW.TmpRTBAtSwitchOutport1_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtSwitchOutport1_RdBufIdx];

  /* RateTransition generated from: '<S49>/Switch1' */
  if (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 == 1) {
    speedgoat_target_model_2021b_DW.TmpRTBAtSwitch1Outport1_RdBufId =
      static_cast<int8_T>
      (speedgoat_target_model_2021b_DW.TmpRTBAtSwitch1Outport1_RdBufId == 0);
  }

  /* RateTransition generated from: '<S49>/Switch1' */
  speedgoat_target_model_2021b_B.TmpRTBAtSwitch1Outport1 =
    speedgoat_target_model_2021b_DW.TmpRTBAtSwitch1Outport1_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtSwitch1Outport1_RdBufId];

  /* SwitchCase: '<S11>/Switch Case1' */
  switch (speedgoat_target_model_2021b_B.OR) {
   case 1:
    /* Outputs for IfAction SubSystem: '<S11>/Subsystem2' incorporates:
     *  ActionPort: '<S157>/Action Port'
     */
    /* MATLAB Function: '<S157>/MATLAB Function' */
    if (speedgoat_target_model_2021b_B.acceleration > 0.0) {
      speedgoat_target_model_2021b_B.accel =
        speedgoat_target_model_2021b_B.acceleration;
      speedgoat_target_model_2021b_B.decel = 0.0;
    } else if (speedgoat_target_model_2021b_B.acceleration < 0.0) {
      speedgoat_target_model_2021b_B.accel = 0.0;
      speedgoat_target_model_2021b_B.decel =
        -speedgoat_target_model_2021b_B.acceleration;
    } else {
      speedgoat_target_model_2021b_B.accel = 0.0;
      speedgoat_target_model_2021b_B.decel = 0.0;
    }

    /* End of MATLAB Function: '<S157>/MATLAB Function' */

    /* DiscreteIntegrator: '<S164>/Discrete-Time Integrator' */
    if ((speedgoat_target_model_2021b_B.decel <= 0.0) &&
        (speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevRe_b == 1))
    {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_DSTATE_f =
        speedgoat_target_model_2021_cal->DiscreteTimeIntegrator_IC;
    }

    /* DiscreteIntegrator: '<S164>/Discrete-Time Integrator' */
    speedgoat_target_model_2021b_B.v_p =
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_DSTATE_f;

    /* Math: '<S164>/Square' */
    speedgoat_target_model_2021b_B.Square_j = speedgoat_target_model_2021b_B.v_p
      * speedgoat_target_model_2021b_B.v_p;

    /* Gain: '<S164>/Gain1' */
    speedgoat_target_model_2021b_B.F_aero_h =
      speedgoat_target_model_2021_cal->Gain1_Gain *
      speedgoat_target_model_2021b_B.Square_j;

    /* Gain: '<S164>/Gain2' incorporates:
     *  Constant: '<S164>/Constant1'
     */
    speedgoat_target_model_2021b_B.Gain2_p =
      speedgoat_target_model_2021_cal->Gain2_Gain *
      speedgoat_target_model_2021_cal->Constant1_Value_i;

    /* Gain: '<S164>/Gain' */
    speedgoat_target_model_2021b_B.Gain_g =
      speedgoat_target_model_2021_cal->Gain_Gain *
      speedgoat_target_model_2021b_B.decel;

    /* Sum: '<S164>/Sum' */
    speedgoat_target_model_2021b_B.Sum_d =
      (speedgoat_target_model_2021b_B.F_aero_h +
       speedgoat_target_model_2021b_B.Gain2_p) +
      speedgoat_target_model_2021b_B.Gain_g;

    /* Gain: '<S164>/Gain4' */
    speedgoat_target_model_2021b_B.T_eng_c =
      speedgoat_target_model_2021_cal->Gain4_Gain *
      speedgoat_target_model_2021b_B.Sum_d;

    /* Saturate: '<S157>/Brake Saturation' */
    absx = speedgoat_target_model_2021b_B.T_eng_c;
    tHat_idx_0 = speedgoat_target_model_2021_cal->BrakeSaturation_LowerSat;
    tHat_idx_1 = speedgoat_target_model_2021_cal->BrakeSaturation_UpperSat;
    if (absx > tHat_idx_1) {
      /* Saturate: '<S157>/Brake Saturation' */
      speedgoat_target_model_2021b_B.BrakeSaturation_n = tHat_idx_1;
    } else if (absx < tHat_idx_0) {
      /* Saturate: '<S157>/Brake Saturation' */
      speedgoat_target_model_2021b_B.BrakeSaturation_n = tHat_idx_0;
    } else {
      /* Saturate: '<S157>/Brake Saturation' */
      speedgoat_target_model_2021b_B.BrakeSaturation_n = absx;
    }

    /* End of Saturate: '<S157>/Brake Saturation' */

    /* DiscreteIntegrator: '<S163>/Discrete-Time Integrator' */
    if ((speedgoat_target_model_2021b_B.accel <= 0.0) &&
        (speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevR_bg == 1))
    {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_DSTATE_e =
        speedgoat_target_model_2021_cal->DiscreteTimeIntegrator_IC_a;
    }

    /* DiscreteIntegrator: '<S163>/Discrete-Time Integrator' */
    speedgoat_target_model_2021b_B.v_o =
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_DSTATE_e;

    /* Gain: '<S163>/Gain' */
    speedgoat_target_model_2021b_B.Gain_b =
      speedgoat_target_model_2021_cal->Gain_Gain_n *
      speedgoat_target_model_2021b_B.accel;

    /* Math: '<S163>/Square' */
    speedgoat_target_model_2021b_B.Square_k = speedgoat_target_model_2021b_B.v_o
      * speedgoat_target_model_2021b_B.v_o;

    /* Gain: '<S163>/Gain1' */
    speedgoat_target_model_2021b_B.F_aero_hy =
      speedgoat_target_model_2021_cal->Gain1_Gain_d *
      speedgoat_target_model_2021b_B.Square_k;

    /* Gain: '<S163>/Gain2' incorporates:
     *  Constant: '<S163>/Constant1'
     */
    speedgoat_target_model_2021b_B.Gain2_b =
      speedgoat_target_model_2021_cal->Gain2_Gain_n *
      speedgoat_target_model_2021_cal->Constant1_Value_d;

    /* Sum: '<S163>/Sum' */
    speedgoat_target_model_2021b_B.Sum_f =
      (speedgoat_target_model_2021b_B.F_aero_hy +
       speedgoat_target_model_2021b_B.Gain2_b) +
      speedgoat_target_model_2021b_B.Gain_b;

    /* Gain: '<S163>/Gain4' */
    speedgoat_target_model_2021b_B.T_eng_c0 =
      speedgoat_target_model_2021_cal->Gain4_Gain_b *
      speedgoat_target_model_2021b_B.Sum_f;

    /* Saturate: '<S157>/Throttle Saturation' */
    absx = speedgoat_target_model_2021b_B.T_eng_c0;
    tHat_idx_0 = speedgoat_target_model_2021_cal->ThrottleSaturation_LowerSat;
    tHat_idx_1 = speedgoat_target_model_2021_cal->ThrottleSaturation_UpperSat;
    if (absx > tHat_idx_1) {
      /* Saturate: '<S157>/Throttle Saturation' */
      speedgoat_target_model_2021b_B.ThrottleSaturation_m = tHat_idx_1;
    } else if (absx < tHat_idx_0) {
      /* Saturate: '<S157>/Throttle Saturation' */
      speedgoat_target_model_2021b_B.ThrottleSaturation_m = tHat_idx_0;
    } else {
      /* Saturate: '<S157>/Throttle Saturation' */
      speedgoat_target_model_2021b_B.ThrottleSaturation_m = absx;
    }

    /* End of Saturate: '<S157>/Throttle Saturation' */

    /* Switch: '<S157>/Switch' */
    if (speedgoat_target_model_2021b_B.ThrottleSaturation_m >
        speedgoat_target_model_2021_cal->Switch_Threshold) {
      /* Merge: '<S11>/Merge1' */
      speedgoat_target_model_2021b_B.Merge1 =
        speedgoat_target_model_2021b_B.ThrottleSaturation_m;
    } else {
      /* Merge: '<S11>/Merge1' incorporates:
       *  Constant: '<S157>/Constant'
       */
      speedgoat_target_model_2021b_B.Merge1 =
        speedgoat_target_model_2021_cal->Constant_Value;
    }

    /* End of Switch: '<S157>/Switch' */

    /* Switch: '<S157>/Switch1' */
    if (speedgoat_target_model_2021b_B.BrakeSaturation_n >
        speedgoat_target_model_2021_cal->Switch1_Threshold) {
      /* Merge: '<S11>/Merge2' incorporates:
       *  Constant: '<S157>/Constant'
       */
      speedgoat_target_model_2021b_B.Merge2 =
        speedgoat_target_model_2021_cal->Constant_Value;
    } else {
      /* Merge: '<S11>/Merge2' */
      speedgoat_target_model_2021b_B.Merge2 =
        speedgoat_target_model_2021b_B.BrakeSaturation_n;
    }

    /* End of Switch: '<S157>/Switch1' */

    /* Update for DiscreteIntegrator: '<S164>/Discrete-Time Integrator' */
    speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_DSTATE_f +=
      speedgoat_target_model_2021_cal->DiscreteTimeIntegrator_gainval *
      speedgoat_target_model_2021b_B.decel;
    if (speedgoat_target_model_2021b_B.decel > 0.0) {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevRe_b = 1;
    } else if (speedgoat_target_model_2021b_B.decel < 0.0) {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevRe_b = -1;
    } else if (speedgoat_target_model_2021b_B.decel == 0.0) {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevRe_b = 0;
    } else {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevRe_b = 2;
    }

    /* End of Update for DiscreteIntegrator: '<S164>/Discrete-Time Integrator' */

    /* Update for DiscreteIntegrator: '<S163>/Discrete-Time Integrator' */
    speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_DSTATE_e +=
      speedgoat_target_model_2021_cal->DiscreteTimeIntegrator_gainva_j *
      speedgoat_target_model_2021b_B.accel;
    if (speedgoat_target_model_2021b_B.accel > 0.0) {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevR_bg = 1;
    } else if (speedgoat_target_model_2021b_B.accel < 0.0) {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevR_bg = -1;
    } else if (speedgoat_target_model_2021b_B.accel == 0.0) {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevR_bg = 0;
    } else {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevR_bg = 2;
    }

    /* End of Update for DiscreteIntegrator: '<S163>/Discrete-Time Integrator' */
    /* End of Outputs for SubSystem: '<S11>/Subsystem2' */

    /* Update for IfAction SubSystem: '<S11>/Subsystem2' incorporates:
     *  ActionPort: '<S157>/Action Port'
     */
    /* Update for SwitchCase: '<S11>/Switch Case1' */
    srUpdateBC(speedgoat_target_model_2021b_DW.Subsystem2_SubsysRanBC);

    /* End of Update for SubSystem: '<S11>/Subsystem2' */
    break;

   case 0:
    /* Outputs for IfAction SubSystem: '<S11>/Subsystem1' incorporates:
     *  ActionPort: '<S156>/Action Port'
     */
    /* DiscreteIntegrator: '<S160>/Discrete-Time Integrator' */
    if ((speedgoat_target_model_2021b_B.TmpRTBAtSwitch1Outport1 <= 0.0) &&
        (speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevRese == 1))
    {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_DSTATE =
        speedgoat_target_model_2021_cal->DiscreteTimeIntegrator_IC_i;
    }

    /* DiscreteIntegrator: '<S160>/Discrete-Time Integrator' */
    speedgoat_target_model_2021b_B.v =
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_DSTATE;

    /* Math: '<S160>/Square' */
    speedgoat_target_model_2021b_B.Square = speedgoat_target_model_2021b_B.v *
      speedgoat_target_model_2021b_B.v;

    /* Gain: '<S160>/Gain1' */
    speedgoat_target_model_2021b_B.F_aero =
      speedgoat_target_model_2021_cal->Gain1_Gain_m *
      speedgoat_target_model_2021b_B.Square;

    /* Gain: '<S160>/Gain2' incorporates:
     *  Constant: '<S160>/Constant1'
     */
    speedgoat_target_model_2021b_B.Gain2 =
      speedgoat_target_model_2021_cal->Gain2_Gain_l *
      speedgoat_target_model_2021_cal->Constant1_Value_o;

    /* Gain: '<S160>/Gain' */
    speedgoat_target_model_2021b_B.Gain_dx =
      speedgoat_target_model_2021_cal->Gain_Gain_m *
      speedgoat_target_model_2021b_B.TmpRTBAtSwitch1Outport1;

    /* Sum: '<S160>/Sum' */
    speedgoat_target_model_2021b_B.Sum = (speedgoat_target_model_2021b_B.F_aero
      + speedgoat_target_model_2021b_B.Gain2) +
      speedgoat_target_model_2021b_B.Gain_dx;

    /* Gain: '<S160>/Gain4' */
    speedgoat_target_model_2021b_B.T_eng =
      speedgoat_target_model_2021_cal->Gain4_Gain_m *
      speedgoat_target_model_2021b_B.Sum;

    /* Saturate: '<S156>/Brake Saturation' */
    absx = speedgoat_target_model_2021b_B.T_eng;
    tHat_idx_0 = speedgoat_target_model_2021_cal->BrakeSaturation_LowerSat_m;
    tHat_idx_1 = speedgoat_target_model_2021_cal->BrakeSaturation_UpperSat_a;
    if (absx > tHat_idx_1) {
      /* Saturate: '<S156>/Brake Saturation' */
      speedgoat_target_model_2021b_B.BrakeSaturation = tHat_idx_1;
    } else if (absx < tHat_idx_0) {
      /* Saturate: '<S156>/Brake Saturation' */
      speedgoat_target_model_2021b_B.BrakeSaturation = tHat_idx_0;
    } else {
      /* Saturate: '<S156>/Brake Saturation' */
      speedgoat_target_model_2021b_B.BrakeSaturation = absx;
    }

    /* End of Saturate: '<S156>/Brake Saturation' */

    /* DiscreteIntegrator: '<S159>/Discrete-Time Integrator' */
    if ((speedgoat_target_model_2021b_B.TmpRTBAtSwitchOutport1 <= 0.0) &&
        (speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevRe_a == 1))
    {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_DSTATE_a =
        speedgoat_target_model_2021_cal->DiscreteTimeIntegrator_IC_d;
    }

    /* DiscreteIntegrator: '<S159>/Discrete-Time Integrator' */
    speedgoat_target_model_2021b_B.v_g =
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_DSTATE_a;

    /* Gain: '<S159>/Gain' */
    speedgoat_target_model_2021b_B.Gain_h =
      speedgoat_target_model_2021_cal->Gain_Gain_c *
      speedgoat_target_model_2021b_B.TmpRTBAtSwitchOutport1;

    /* Math: '<S159>/Square' */
    speedgoat_target_model_2021b_B.Square_i = speedgoat_target_model_2021b_B.v_g
      * speedgoat_target_model_2021b_B.v_g;

    /* Gain: '<S159>/Gain1' */
    speedgoat_target_model_2021b_B.F_aero_k =
      speedgoat_target_model_2021_cal->Gain1_Gain_j *
      speedgoat_target_model_2021b_B.Square_i;

    /* Gain: '<S159>/Gain2' incorporates:
     *  Constant: '<S159>/Constant1'
     */
    speedgoat_target_model_2021b_B.Gain2_k =
      speedgoat_target_model_2021_cal->Gain2_Gain_h *
      speedgoat_target_model_2021_cal->Constant1_Value_fy;

    /* Sum: '<S159>/Sum' */
    speedgoat_target_model_2021b_B.Sum_o =
      (speedgoat_target_model_2021b_B.F_aero_k +
       speedgoat_target_model_2021b_B.Gain2_k) +
      speedgoat_target_model_2021b_B.Gain_h;

    /* Gain: '<S159>/Gain4' */
    speedgoat_target_model_2021b_B.T_eng_h =
      speedgoat_target_model_2021_cal->Gain4_Gain_n *
      speedgoat_target_model_2021b_B.Sum_o;

    /* Saturate: '<S156>/Throttle Saturation' */
    absx = speedgoat_target_model_2021b_B.T_eng_h;
    tHat_idx_0 = speedgoat_target_model_2021_cal->ThrottleSaturation_LowerSat_g;
    tHat_idx_1 = speedgoat_target_model_2021_cal->ThrottleSaturation_UpperSat_e;
    if (absx > tHat_idx_1) {
      /* Saturate: '<S156>/Throttle Saturation' */
      speedgoat_target_model_2021b_B.ThrottleSaturation = tHat_idx_1;
    } else if (absx < tHat_idx_0) {
      /* Saturate: '<S156>/Throttle Saturation' */
      speedgoat_target_model_2021b_B.ThrottleSaturation = tHat_idx_0;
    } else {
      /* Saturate: '<S156>/Throttle Saturation' */
      speedgoat_target_model_2021b_B.ThrottleSaturation = absx;
    }

    /* End of Saturate: '<S156>/Throttle Saturation' */

    /* Switch: '<S156>/Switch' */
    if (speedgoat_target_model_2021b_B.ThrottleSaturation >
        speedgoat_target_model_2021_cal->Switch_Threshold_f) {
      /* Merge: '<S11>/Merge1' */
      speedgoat_target_model_2021b_B.Merge1 =
        speedgoat_target_model_2021b_B.ThrottleSaturation;
    } else {
      /* Merge: '<S11>/Merge1' incorporates:
       *  Constant: '<S156>/Constant'
       */
      speedgoat_target_model_2021b_B.Merge1 =
        speedgoat_target_model_2021_cal->Constant_Value_d;
    }

    /* End of Switch: '<S156>/Switch' */

    /* Switch: '<S156>/Switch1' */
    if (speedgoat_target_model_2021b_B.BrakeSaturation >
        speedgoat_target_model_2021_cal->Switch1_Threshold_m) {
      /* Merge: '<S11>/Merge2' incorporates:
       *  Constant: '<S156>/Constant'
       */
      speedgoat_target_model_2021b_B.Merge2 =
        speedgoat_target_model_2021_cal->Constant_Value_d;
    } else {
      /* Merge: '<S11>/Merge2' */
      speedgoat_target_model_2021b_B.Merge2 =
        speedgoat_target_model_2021b_B.BrakeSaturation;
    }

    /* End of Switch: '<S156>/Switch1' */

    /* Update for DiscreteIntegrator: '<S160>/Discrete-Time Integrator' */
    speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_DSTATE +=
      speedgoat_target_model_2021_cal->DiscreteTimeIntegrator_gainva_c *
      speedgoat_target_model_2021b_B.TmpRTBAtSwitch1Outport1;
    if (speedgoat_target_model_2021b_B.TmpRTBAtSwitch1Outport1 > 0.0) {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevRese = 1;
    } else if (speedgoat_target_model_2021b_B.TmpRTBAtSwitch1Outport1 < 0.0) {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevRese = -1;
    } else if (speedgoat_target_model_2021b_B.TmpRTBAtSwitch1Outport1 == 0.0) {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevRese = 0;
    } else {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevRese = 2;
    }

    /* End of Update for DiscreteIntegrator: '<S160>/Discrete-Time Integrator' */

    /* Update for DiscreteIntegrator: '<S159>/Discrete-Time Integrator' */
    speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_DSTATE_a +=
      speedgoat_target_model_2021_cal->DiscreteTimeIntegrator_gainva_a *
      speedgoat_target_model_2021b_B.TmpRTBAtSwitchOutport1;
    if (speedgoat_target_model_2021b_B.TmpRTBAtSwitchOutport1 > 0.0) {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevRe_a = 1;
    } else if (speedgoat_target_model_2021b_B.TmpRTBAtSwitchOutport1 < 0.0) {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevRe_a = -1;
    } else if (speedgoat_target_model_2021b_B.TmpRTBAtSwitchOutport1 == 0.0) {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevRe_a = 0;
    } else {
      speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevRe_a = 2;
    }

    /* End of Update for DiscreteIntegrator: '<S159>/Discrete-Time Integrator' */
    /* End of Outputs for SubSystem: '<S11>/Subsystem1' */

    /* Update for IfAction SubSystem: '<S11>/Subsystem1' incorporates:
     *  ActionPort: '<S156>/Action Port'
     */
    /* Update for SwitchCase: '<S11>/Switch Case1' */
    srUpdateBC(speedgoat_target_model_2021b_DW.Subsystem1_SubsysRanBC);

    /* End of Update for SubSystem: '<S11>/Subsystem1' */
    break;
  }

  /* End of SwitchCase: '<S11>/Switch Case1' */

  /* S-Function (slrealtimebytepacking): '<S181>/Byte Packing' */

  /* Byte Packing: <S181>/Byte Packing */
  (void)memcpy((uint8_T*)&speedgoat_target_model_2021b_B.BytePacking[0] + 0,
               (uint8_T*)&speedgoat_target_model_2021b_B.Saturation, 8);
  (void)memcpy((uint8_T*)&speedgoat_target_model_2021b_B.BytePacking[0] + 8,
               (uint8_T*)&speedgoat_target_model_2021b_B.Merge1, 8);
  (void)memcpy((uint8_T*)&speedgoat_target_model_2021b_B.BytePacking[0] + 16,
               (uint8_T*)&speedgoat_target_model_2021b_B.Merge2, 8);

  /* S-Function (slrealtimeUDPSend): '<S3>/UDP Send: Controller Output' */
  {
    try {
      slrealtime::ip::udp::Socket *udpSock = reinterpret_cast<slrealtime::ip::
        udp::Socket*>
        (speedgoat_target_model_2021b_DW.UDPSendControllerOutput_PWORK);
      uint8_t *remoteAddress = (uint8_t *)
        speedgoat_target_model_2021_cal->UDPSendControllerOutput_toAddre;
      uint16_t *remotePort = (uint16_t *)
        &speedgoat_target_model_2021_cal->UDPSendControllerOutput_toPort;
      udpSock->resetRemoteEndpoint(remoteAddress, remotePort);
      int dataLen = speedgoat_target_model_2_ConstB.Width;
      dataLen = (dataLen <=
                 speedgoat_target_model_2021b_DW.UDPSendControllerOutput_IWORK[0])
        ? dataLen :
        speedgoat_target_model_2021b_DW.UDPSendControllerOutput_IWORK[0];
      dataLen = (dataLen <= -1) ? 0 : dataLen;
      void *dataPort = &speedgoat_target_model_2021b_B.BytePacking[0];
      size_t numBytesSend = udpSock->send((const char*)dataPort,dataLen);
    } catch (std::exception& e) {
      std::string tmp = std::string(e.what());
      static std::string eMsg = tmp;
      rtmSetErrorStatus(speedgoat_target_model_2021b_M, eMsg.c_str());
      rtmSetStopRequested(speedgoat_target_model_2021b_M, 1);
      ;
    }
  }

  /* S-Function (slrealtimebytepacking): '<S180>/Byte Packing1' */

  /* Byte Packing: <S180>/Byte Packing1 */
  (void)memcpy((uint8_T*)&speedgoat_target_model_2021b_B.BytePacking1[0] + 0,
               (uint8_T*)&speedgoat_target_model_2021b_B.StopSim, 8);
  (void)memcpy((uint8_T*)&speedgoat_target_model_2021b_B.BytePacking1[0] + 8,
               (uint8_T*)&speedgoat_target_model_2021b_B.Trans, 8);

  /* S-Function (slrealtimeUDPSend): '<S3>/UDP Send: Controller Intermediate Signals' */
  {
    try {
      slrealtime::ip::udp::Socket *udpSock = reinterpret_cast<slrealtime::ip::
        udp::Socket*>
        (speedgoat_target_model_2021b_DW.UDPSendControllerIntermediateSi);
      uint8_t *remoteAddress = (uint8_t *)
        speedgoat_target_model_2021_cal->UDPSendControllerIntermediate_k;
      uint16_t *remotePort = (uint16_t *)
        &speedgoat_target_model_2021_cal->UDPSendControllerIntermediateSi;
      udpSock->resetRemoteEndpoint(remoteAddress, remotePort);
      int dataLen = speedgoat_target_model_2_ConstB.Width1;
      dataLen = (dataLen <=
                 speedgoat_target_model_2021b_DW.UDPSendControllerIntermediate_a[
                 0]) ? dataLen :
        speedgoat_target_model_2021b_DW.UDPSendControllerIntermediate_a[0];
      dataLen = (dataLen <= -1) ? 0 : dataLen;
      void *dataPort = &speedgoat_target_model_2021b_B.BytePacking1[0];
      size_t numBytesSend = udpSock->send((const char*)dataPort,dataLen);
    } catch (std::exception& e) {
      std::string tmp = std::string(e.what());
      static std::string eMsg = tmp;
      rtmSetErrorStatus(speedgoat_target_model_2021b_M, eMsg.c_str());
      rtmSetStopRequested(speedgoat_target_model_2021b_M, 1);
      ;
    }
  }

  /* Product: '<S4>/Product2' */
  speedgoat_target_model_2021b_B.Product2 =
    speedgoat_target_model_2021b_B.Longitudinalvelocity *
    speedgoat_target_model_2021b_B.Tgap;

  /* Sum: '<S4>/Sum1' incorporates:
   *  Constant: '<S4>/Default spacing constant'
   */
  speedgoat_target_model_2021b_B.safe_distance =
    speedgoat_target_model_2021_cal->ACCsystem_DefaultSpacing +
    speedgoat_target_model_2021b_B.Product2;

  /* MATLAB Function: '<S4>/DataTypeConversion_dmin' */
  speed_DataTypeConversion_L0(speedgoat_target_model_2021b_B.safe_distance,
    &speedgoat_target_model_2021b_B.sf_DataTypeConversion_dmin);

  /* MATLAB Function: '<S4>/DataTypeConversion_reldist' */
  speed_DataTypeConversion_L0(speedgoat_target_model_2021b_B.Relativedistance,
    &speedgoat_target_model_2021b_B.sf_DataTypeConversion_reldist);

  /* MATLAB Function: '<S4>/DataTypeConversion_vego' */
  speed_DataTypeConversion_L0
    (speedgoat_target_model_2021b_B.Longitudinalvelocity,
     &speedgoat_target_model_2021b_B.sf_DataTypeConversion_vego);

  /* Sum: '<S4>/Sum6' */
  speedgoat_target_model_2021b_B.lead_velocity =
    speedgoat_target_model_2021b_B.Longitudinalvelocity +
    speedgoat_target_model_2021b_B.Relativevelocity;

  /* MATLAB Function: '<S4>/DataTypeConversion_vlead' */
  speed_DataTypeConversion_L0(speedgoat_target_model_2021b_B.lead_velocity,
    &speedgoat_target_model_2021b_B.sf_DataTypeConversion_vlead);

  /* MATLAB Function: '<S1>/MATLAB Function2' */
  if ((std::abs(speedgoat_target_model_2021b_B.Curvature) >= 0.0) && (std::abs
       (speedgoat_target_model_2021b_B.Curvature) <= 0.01)) {
    speedgoat_target_model_2021b_B.RefVelCmd =
      speedgoat_target_model_2021b_B.VelCmd;
  } else {
    absx = speedgoat_target_model_2021b_B.VelCmd;
    tHat_idx_0 = std::sqrt(1.0 / std::abs
      (speedgoat_target_model_2021b_B.Curvature) * 9.81 * 0.7);
    if ((absx <= tHat_idx_0) || rtIsNaN(tHat_idx_0)) {
      speedgoat_target_model_2021b_B.RefVelCmd = absx;
    } else {
      speedgoat_target_model_2021b_B.RefVelCmd = tHat_idx_0;
    }
  }

  /* End of MATLAB Function: '<S1>/MATLAB Function2' */

  /* MATLAB Function: '<S4>/DataTypeConversion_vset' */
  speed_DataTypeConversion_L0(speedgoat_target_model_2021b_B.RefVelCmd,
    &speedgoat_target_model_2021b_B.sf_DataTypeConversion_vset);

  /* Rounding: '<S23>/Floor' incorporates:
   *  Constant: '<S22>/p_zero'
   */
  speedgoat_target_model_2021b_B.Floor = std::floor
    (speedgoat_target_model_2021_cal->p_zero_Value);

  /* Rounding: '<S23>/Floor1' incorporates:
   *  Constant: '<S22>/m_zero'
   */
  speedgoat_target_model_2021b_B.Floor1 = std::floor
    (speedgoat_target_model_2021_cal->m_zero_Value);

  /* RateTransition generated from: '<S43>/optimizer' */
  if (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 == 1) {
    speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport3_WrBufI =
      static_cast<int8_T>
      (speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport3_WrBufI == 0);
  }

  speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport3_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport3_WrBufI
    << 1] = speedgoat_target_model_2021b_B.sf_DataTypeConversion_reldist.y;
  speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport3_Buf
    [(speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport3_WrBufI << 1) + 1]
    = speedgoat_target_model_2021b_B.sf_DataTypeConversion_vego.y;

  /* End of RateTransition generated from: '<S43>/optimizer' */

  /* RateTransition generated from: '<S43>/optimizer' */
  if (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 == 1) {
    speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport5_WrBufI =
      static_cast<int8_T>
      (speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport5_WrBufI == 0);
  }

  speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport5_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport5_WrBufI]
    = speedgoat_target_model_2021b_B.sf_DataTypeConversion_vlead.y;

  /* End of RateTransition generated from: '<S43>/optimizer' */

  /* RateTransition generated from: '<S43>/optimizer' incorporates:
   *  Constant: '<S4>/Minimum velocity constant'
   */
  rtw_slrealtime_mutex_lock
    (speedgoat_target_model_2021b_DW.min_velocity_d0_SEMAPHORE);
  wrBufIdx = static_cast<int8_T>
    (speedgoat_target_model_2021b_DW.min_velocity_LstBufWR + 1);
  if (wrBufIdx == 3) {
    wrBufIdx = 0;
  }

  if (wrBufIdx == speedgoat_target_model_2021b_DW.min_velocity_RDBuf) {
    wrBufIdx = static_cast<int8_T>(wrBufIdx + 1);
    if (wrBufIdx == 3) {
      wrBufIdx = 0;
    }
  }

  rtw_slrealtime_mutex_unlock
    (speedgoat_target_model_2021b_DW.min_velocity_d0_SEMAPHORE);
  switch (wrBufIdx) {
   case 0:
    speedgoat_target_model_2021b_DW.min_velocity_Buf0[0] =
      speedgoat_target_model_2021b_B.sf_DataTypeConversion_dmin.y;
    speedgoat_target_model_2021b_DW.min_velocity_Buf0[1] =
      speedgoat_target_model_2021_cal->Minimumvelocityconstant_Value;
    break;

   case 1:
    speedgoat_target_model_2021b_DW.min_velocity_Buf1[0] =
      speedgoat_target_model_2021b_B.sf_DataTypeConversion_dmin.y;
    speedgoat_target_model_2021b_DW.min_velocity_Buf1[1] =
      speedgoat_target_model_2021_cal->Minimumvelocityconstant_Value;
    break;

   case 2:
    speedgoat_target_model_2021b_DW.min_velocity_Buf2[0] =
      speedgoat_target_model_2021b_B.sf_DataTypeConversion_dmin.y;
    speedgoat_target_model_2021b_DW.min_velocity_Buf2[1] =
      speedgoat_target_model_2021_cal->Minimumvelocityconstant_Value;
    break;
  }

  speedgoat_target_model_2021b_DW.min_velocity_LstBufWR = wrBufIdx;

  /* End of RateTransition generated from: '<S43>/optimizer' */

  /* Constant: '<S23>/ym_zero' */
  speedgoat_target_model_2021b_B.ym_zero[0] =
    speedgoat_target_model_2021_cal->ym_zero_Value[0];
  speedgoat_target_model_2021b_B.ym_zero[1] =
    speedgoat_target_model_2021_cal->ym_zero_Value[1];

  /* RateTransition generated from: '<S4>/DataTypeConversion_vset' */
  if (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 == 1) {
    speedgoat_target_model_2021b_DW.TmpRTBAtDataTypeConversion_vs_m =
      static_cast<int8_T>
      (speedgoat_target_model_2021b_DW.TmpRTBAtDataTypeConversion_vs_m == 0);
  }

  speedgoat_target_model_2021b_DW.TmpRTBAtDataTypeConversion_vset[speedgoat_target_model_2021b_DW.TmpRTBAtDataTypeConversion_vs_m]
    = speedgoat_target_model_2021b_B.sf_DataTypeConversion_vset.y;

  /* End of RateTransition generated from: '<S4>/DataTypeConversion_vset' */

  /* RateTransition generated from: '<S4>/DataTypeConversion_L0' incorporates:
   *  Constant: '<S4>/Default spacing constant'
   */
  rtw_slrealtime_mutex_lock
    (speedgoat_target_model_2021b_DW.default_spacing_d0_SEMAPHORE);
  wrBufIdx = static_cast<int8_T>
    (speedgoat_target_model_2021b_DW.default_spacing_LstBufWR + 1);
  if (wrBufIdx == 3) {
    wrBufIdx = 0;
  }

  if (wrBufIdx == speedgoat_target_model_2021b_DW.default_spacing_RDBuf) {
    wrBufIdx = static_cast<int8_T>(wrBufIdx + 1);
    if (wrBufIdx == 3) {
      wrBufIdx = 0;
    }
  }

  rtw_slrealtime_mutex_unlock
    (speedgoat_target_model_2021b_DW.default_spacing_d0_SEMAPHORE);
  switch (wrBufIdx) {
   case 0:
    speedgoat_target_model_2021b_DW.default_spacing_Buf0 =
      speedgoat_target_model_2021_cal->ACCsystem_DefaultSpacing;
    break;

   case 1:
    speedgoat_target_model_2021b_DW.default_spacing_Buf1 =
      speedgoat_target_model_2021_cal->ACCsystem_DefaultSpacing;
    break;

   case 2:
    speedgoat_target_model_2021b_DW.default_spacing_Buf2 =
      speedgoat_target_model_2021_cal->ACCsystem_DefaultSpacing;
    break;
  }

  speedgoat_target_model_2021b_DW.default_spacing_LstBufWR = wrBufIdx;

  /* End of RateTransition generated from: '<S4>/DataTypeConversion_L0' */

  /* MATLAB Function: '<S4>/DataTypeConversion_atrack' incorporates:
   *  Constant: '<S4>/External control signal constant'
   */
  spe_DataTypeConversion_amax
    (speedgoat_target_model_2021_cal->Externalcontrolsignalconstant_V,
     &speedgoat_target_model_2021b_B.sf_DataTypeConversion_atrack);

  /* RateTransition generated from: '<S51>/Equal1' */
  if (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 == 1) {
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual1Inport1_WrBufIdx =
      static_cast<int8_T>
      (speedgoat_target_model_2021b_DW.TmpRTBAtEqual1Inport1_WrBufIdx == 0);
  }

  speedgoat_target_model_2021b_DW.TmpRTBAtEqual1Inport1_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtEqual1Inport1_WrBufIdx]
    = speedgoat_target_model_2021b_B.Direction;

  /* End of RateTransition generated from: '<S51>/Equal1' */

  /* RateTransition generated from: '<S51>/Equal2' */
  if (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 == 1) {
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual2Inport1_WrBufIdx =
      static_cast<int8_T>
      (speedgoat_target_model_2021b_DW.TmpRTBAtEqual2Inport1_WrBufIdx == 0);
  }

  speedgoat_target_model_2021b_DW.TmpRTBAtEqual2Inport1_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtEqual2Inport1_WrBufIdx]
    = speedgoat_target_model_2021b_B.Direction;

  /* End of RateTransition generated from: '<S51>/Equal2' */

  /* RateTransition generated from: '<S51>/Equal3' */
  if (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 == 1) {
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual3Inport1_WrBufIdx =
      static_cast<int8_T>
      (speedgoat_target_model_2021b_DW.TmpRTBAtEqual3Inport1_WrBufIdx == 0);
  }

  speedgoat_target_model_2021b_DW.TmpRTBAtEqual3Inport1_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtEqual3Inport1_WrBufIdx]
    = speedgoat_target_model_2021b_B.Longitudinalvelocity;

  /* End of RateTransition generated from: '<S51>/Equal3' */

  /* RateTransition generated from: '<S51>/Equal4' */
  if (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 == 1) {
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual4Inport2_WrBufIdx =
      static_cast<int8_T>
      (speedgoat_target_model_2021b_DW.TmpRTBAtEqual4Inport2_WrBufIdx == 0);
  }

  speedgoat_target_model_2021b_DW.TmpRTBAtEqual4Inport2_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtEqual4Inport2_WrBufIdx]
    = speedgoat_target_model_2021b_B.Direction;

  /* End of RateTransition generated from: '<S51>/Equal4' */

  /* RateTransition generated from: '<S51>/Equal5' */
  if (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 == 1) {
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual5Inport2_WrBufIdx =
      static_cast<int8_T>
      (speedgoat_target_model_2021b_DW.TmpRTBAtEqual5Inport2_WrBufIdx == 0);
  }

  speedgoat_target_model_2021b_DW.TmpRTBAtEqual5Inport2_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtEqual5Inport2_WrBufIdx]
    = speedgoat_target_model_2021b_B.Direction;

  /* End of RateTransition generated from: '<S51>/Equal5' */

  /* RateTransition generated from: '<S51>/Sign1' */
  if (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 == 1) {
    speedgoat_target_model_2021b_DW.TmpRTBAtSign1Inport1_WrBufIdx =
      static_cast<int8_T>
      (speedgoat_target_model_2021b_DW.TmpRTBAtSign1Inport1_WrBufIdx == 0);
  }

  speedgoat_target_model_2021b_DW.TmpRTBAtSign1Inport1_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtSign1Inport1_WrBufIdx]
    = speedgoat_target_model_2021b_B.Longitudinalvelocity;

  /* End of RateTransition generated from: '<S51>/Sign1' */

  /* RateTransition generated from: '<S51>/Equal6' */
  if (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 == 1) {
    speedgoat_target_model_2021b_DW.set_velocity_WrBufIdx = static_cast<int8_T>
      (speedgoat_target_model_2021b_DW.set_velocity_WrBufIdx == 0);
  }

  speedgoat_target_model_2021b_DW.set_velocity_Buf[speedgoat_target_model_2021b_DW.set_velocity_WrBufIdx]
    = speedgoat_target_model_2021b_B.RefVelCmd;

  /* End of RateTransition generated from: '<S51>/Equal6' */

  /* RateTransition generated from: '<S51>/Sign2' */
  if (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 == 1) {
    speedgoat_target_model_2021b_DW.set_velocity1_WrBufIdx = static_cast<int8_T>
      (speedgoat_target_model_2021b_DW.set_velocity1_WrBufIdx == 0);
  }

  speedgoat_target_model_2021b_DW.set_velocity1_Buf[speedgoat_target_model_2021b_DW.set_velocity1_WrBufIdx]
    = speedgoat_target_model_2021b_B.RefVelCmd;

  /* End of RateTransition generated from: '<S51>/Sign2' */

  /* RateTransition generated from: '<S49>/Multiply' */
  if (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 == 1) {
    speedgoat_target_model_2021b_DW.TmpRTBAtMultiplyInport1_WrBufId =
      static_cast<int8_T>
      (speedgoat_target_model_2021b_DW.TmpRTBAtMultiplyInport1_WrBufId == 0);
  }

  speedgoat_target_model_2021b_DW.TmpRTBAtMultiplyInport1_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtMultiplyInport1_WrBufId]
    = speedgoat_target_model_2021b_B.Direction;

  /* End of RateTransition generated from: '<S49>/Multiply' */

  /* RateTransition generated from: '<S50>/Switch Case' */
  if (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 == 1) {
    speedgoat_target_model_2021b_DW.TmpRTBAtSwitchCaseInport1_WrBuf =
      static_cast<int8_T>
      (speedgoat_target_model_2021b_DW.TmpRTBAtSwitchCaseInport1_WrBuf == 0);
  }

  speedgoat_target_model_2021b_DW.TmpRTBAtSwitchCaseInport1_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtSwitchCaseInport1_WrBuf]
    = speedgoat_target_model_2021b_B.Direction;

  /* End of RateTransition generated from: '<S50>/Switch Case' */

  /* RateTransition generated from: '<S8>/Minus' */
  if (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 == 1) {
    speedgoat_target_model_2021b_DW.TmpRTBAtMinusInport2_WrBufIdx =
      static_cast<int8_T>
      (speedgoat_target_model_2021b_DW.TmpRTBAtMinusInport2_WrBufIdx == 0);
  }

  speedgoat_target_model_2021b_DW.TmpRTBAtMinusInport2_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtMinusInport2_WrBufIdx]
    = speedgoat_target_model_2021b_B.Longitudinalvelocity;

  /* End of RateTransition generated from: '<S8>/Minus' */

  /* RateTransition generated from: '<S8>/Minus' */
  if (speedgoat_target_model_2021b_M->Timing.RateInteraction.TID0_1 == 1) {
    speedgoat_target_model_2021b_DW.set_velocity_WrBufIdx_k = static_cast<int8_T>
      (speedgoat_target_model_2021b_DW.set_velocity_WrBufIdx_k == 0);
  }

  speedgoat_target_model_2021b_DW.set_velocity_Buf_n[speedgoat_target_model_2021b_DW.set_velocity_WrBufIdx_k]
    = speedgoat_target_model_2021b_B.RefVelCmd;

  /* End of RateTransition generated from: '<S8>/Minus' */

  /* Update for UnitDelay: '<S47>/Unit Delay' */
  speedgoat_target_model_2021b_DW.UnitDelay_DSTATE =
    speedgoat_target_model_2021b_B.UnitDelay_n;

  /* Update for UnitDelay: '<S1>/Unit Delay' */
  speedgoat_target_model_2021b_DW.UnitDelay_DSTATE_o =
    speedgoat_target_model_2021b_B.Saturation;

  /* Update absolute time */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++speedgoat_target_model_2021b_M->Timing.clockTick0)) {
    ++speedgoat_target_model_2021b_M->Timing.clockTickH0;
  }

  speedgoat_target_model_2021b_M->Timing.taskTime0 =
    speedgoat_target_model_2021b_M->Timing.clockTick0 *
    speedgoat_target_model_2021b_M->Timing.stepSize0 +
    speedgoat_target_model_2021b_M->Timing.clockTickH0 *
    speedgoat_target_model_2021b_M->Timing.stepSize0 * 4294967296.0;
}

/* Model step function for TID1 */
void speedgoat_target_model_2021b_step1(void) /* Sample time: [1.0s, 0.0s] */
{
  real_T Bc[34];
  real_T a__1[34];
  real_T vseq[22];
  real_T rseq[20];
  real_T xk[4];
  real_T f[3];
  real_T zopt[3];
  real_T y_innov[2];
  real_T ymax_incr[2];
  real_T ymin_incr[2];
  real_T old_u;
  real_T umax_incr;
  real_T umin_incr;
  int32_T ii;
  int16_T iAnew[34];
  boolean_T ymax_incr_flag[2];
  boolean_T ymin_incr_flag[2];
  boolean_T umax_incr_flag;
  static const real_T b_Mv[748] = { 1.3552527156068805E-19,
    1.0842021724855044E-19, 8.1315162936412833E-20, 1.0842021724855044E-19,
    1.0842021724855044E-19, 1.0842021724855044E-19, 1.3552527156068805E-19,
    1.3552527156068805E-19, 1.3552527156068805E-19, 8.1315162936412833E-20,
    0.99999999999999989, -1.3552527156068805E-19, 0.99999999999999989,
    -1.0842021724855044E-19, 0.99999999999999967, -8.1315162936412833E-20,
    0.99999999999999967, -1.0842021724855044E-19, 0.99999999999999967,
    -1.0842021724855044E-19, 0.99999999999999956, -1.0842021724855044E-19,
    0.99999999999999944, -1.3552527156068805E-19, 0.99999999999999944,
    -1.3552527156068805E-19, 0.99999999999999967, -1.3552527156068805E-19,
    0.99999999999999944, -8.1315162936412833E-20, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -0.0, 1.3552527156068805E-19, 1.0842021724855044E-19,
    8.1315162936412833E-20, 1.0842021724855044E-19, 1.0842021724855044E-19,
    1.0842021724855044E-19, 1.3552527156068805E-19, 1.3552527156068805E-19,
    1.3552527156068805E-19, 0.0, 0.0, 0.99999999999999989,
    -1.3552527156068805E-19, 0.99999999999999989, -1.0842021724855044E-19,
    0.99999999999999967, -8.1315162936412833E-20, 0.99999999999999967,
    -1.0842021724855044E-19, 0.99999999999999967, -1.0842021724855044E-19,
    0.99999999999999956, -1.0842021724855044E-19, 0.99999999999999944,
    -1.3552527156068805E-19, 0.99999999999999944, -1.3552527156068805E-19,
    0.99999999999999967, -1.3552527156068805E-19, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -0.0, -0.0, 1.3552527156068805E-19, 1.0842021724855044E-19,
    8.1315162936412833E-20, 1.0842021724855044E-19, 1.0842021724855044E-19,
    1.0842021724855044E-19, 1.3552527156068805E-19, 1.3552527156068805E-19, 0.0,
    0.0, 0.0, 0.0, 0.99999999999999989, -1.3552527156068805E-19,
    0.99999999999999989, -1.0842021724855044E-19, 0.99999999999999967,
    -8.1315162936412833E-20, 0.99999999999999967, -1.0842021724855044E-19,
    0.99999999999999967, -1.0842021724855044E-19, 0.99999999999999956,
    -1.0842021724855044E-19, 0.99999999999999944, -1.3552527156068805E-19,
    0.99999999999999944, -1.3552527156068805E-19, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -0.0, -0.0, -0.0, 1.3552527156068805E-19,
    1.0842021724855044E-19, 8.1315162936412833E-20, 1.0842021724855044E-19,
    1.0842021724855044E-19, 1.0842021724855044E-19, 1.3552527156068805E-19, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.99999999999999989, -1.3552527156068805E-19,
    0.99999999999999989, -1.0842021724855044E-19, 0.99999999999999967,
    -8.1315162936412833E-20, 0.99999999999999967, -1.0842021724855044E-19,
    0.99999999999999967, -1.0842021724855044E-19, 0.99999999999999956,
    -1.0842021724855044E-19, 0.99999999999999944, -1.3552527156068805E-19, 0.0,
    0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0,
    1.3552527156068805E-19, 1.0842021724855044E-19, 8.1315162936412833E-20,
    1.0842021724855044E-19, 1.0842021724855044E-19, 1.0842021724855044E-19, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99999999999999989,
    -1.3552527156068805E-19, 0.99999999999999989, -1.0842021724855044E-19,
    0.99999999999999967, -8.1315162936412833E-20, 0.99999999999999967,
    -1.0842021724855044E-19, 0.99999999999999967, -1.0842021724855044E-19,
    0.99999999999999956, -1.0842021724855044E-19, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.3552527156068805E-19,
    1.0842021724855044E-19, 8.1315162936412833E-20, 1.0842021724855044E-19,
    1.0842021724855044E-19, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.99999999999999989, -1.3552527156068805E-19, 0.99999999999999989,
    -1.0842021724855044E-19, 0.99999999999999967, -8.1315162936412833E-20,
    0.99999999999999967, -1.0842021724855044E-19, 0.99999999999999967,
    -1.0842021724855044E-19, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.3552527156068805E-19,
    1.0842021724855044E-19, 8.1315162936412833E-20, 1.0842021724855044E-19, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99999999999999989,
    -1.3552527156068805E-19, 0.99999999999999989, -1.0842021724855044E-19,
    0.99999999999999967, -8.1315162936412833E-20, 0.99999999999999967,
    -1.0842021724855044E-19, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.3552527156068805E-19,
    1.0842021724855044E-19, 8.1315162936412833E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99999999999999989,
    -1.3552527156068805E-19, 0.99999999999999989, -1.0842021724855044E-19,
    0.99999999999999967, -8.1315162936412833E-20, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    1.3552527156068805E-19, 1.0842021724855044E-19, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99999999999999989,
    -1.3552527156068805E-19, 0.99999999999999989, -1.0842021724855044E-19, 0.0,
    0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 1.3552527156068805E-19, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.99999999999999989, -1.3552527156068805E-19, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  real_T Bc_0;
  real_T v_idx_0;
  real_T v_idx_1;
  int32_T i;
  static const real_T f_a[8] = { -2.7529251322709752E-20, 0.0063245553203367215,
    -1.1159380919009811E-5, -0.012649109656163445, 0.028284269046030115,
    -4.99062686434405E-6, 0.0, 1.0 };

  static const real_T e_a[8] = { -10.538526224682116, 27.168712167896391,
    28.509166107776466, 0.3489966970615615, -10.197127455558066,
    -30.907649512684081, -2.1850309517854951, 0.41117101504513459 };

  static const real_T b_Mx[136] = { -0.0008559354856233582,
    -0.00011583827137906727, -1.567700526668816E-5, -2.1216519480315983E-6,
    -2.8713436727885496E-7, -3.8859410885137038E-8, -5.2590493410355552E-9,
    -7.1173489461200677E-10, -9.63228060394177E-11, -1.3035836725042097E-11,
    -0.00273430991735664, 0.0008559354856233582, -0.0031043585244787423,
    0.00011583827137906727, -0.0031544391575348881, 1.567700526668816E-5,
    -0.0031612168341941726, 2.1216519480315983E-6, -0.0031621340929845056,
    2.8713436727885496E-7, -0.0031622582304626588, 3.8859410885137038E-8,
    -0.0031622750306433871, 5.2590493410355552E-9, -0.0031622773043005672,
    7.1173489461200677E-10, -0.0031622776120065685, 9.63228060394177E-11,
    -0.0031622776536500093, 1.3035836725042097E-11, 0.0, 0.0, 0.0, 0.0,
    0.015383419360702295, 0.015753467939022679, 0.015803548568180972,
    0.015810326244312774, 0.015811243503031758, 0.015811367640500291,
    0.015811384440679754, 0.015811386714336798, 0.015811387022042822,
    0.015811387063686302, 0.014433072837015538, -0.015383419360702295,
    0.030059435618059321, -0.015753467939022679, 0.045845782373684141,
    -0.015803548568180972, 0.0616537806058222, -0.015810326244312774,
    0.077464709046666677, -0.015811243503031758, 0.093276034048136391,
    -0.015811367640500291, 0.10908741271825063, -0.015811384440679754,
    0.12489879865162608, -0.015811386714336798, 0.14071018556797704,
    -0.015811387022042822, 0.15652157261735927, -0.015811387063686302, 0.0, 0.0,
    0.0, 0.0, 6.0694316053764884E-6, 6.2154319505611135E-6,
    6.2351909486293173E-6, 6.237865038229349E-6, 6.2382269369027687E-6,
    6.238275914562239E-6, 6.2382825429676551E-6, 6.2382834400247792E-6,
    6.238283561428259E-6, 6.2382835778584324E-6, 0.028289967927240027,
    -6.0694316053764884E-6, 0.028296133210647861, -6.2154319505611135E-6,
    0.028302361614729255, -6.2351909486293173E-6, 0.028308598561264883,
    -6.237865038229349E-6, 0.028314836663895976, -6.2382269369027687E-6,
    0.028321074922987573, -6.238275914562239E-6, 0.028327313203253802,
    -6.2382825429676551E-6, 0.028333551486385705, -6.2382834400247792E-6,
    0.028339789769905432, -6.238283561428259E-6, 0.028346028053477648,
    -6.2382835778584324E-6, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0,
    -1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T b_Mlim_0[34] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.5, 0.5, 0.5, 0.5 };

  static const real_T b_Mu1[34] = { -0.90826822658930306, -2.4146525111110129,
    -4.0019830017413716, -5.6002683701023734, -7.2000363199438731,
    -8.8000049153699571, -10.400000665223061, -12.000000090028237,
    -13.600000012184093, -15.200000001649045, -0.34586588670534729,
    0.90826822658930306, -1.9926737444444755, 2.4146525111110129,
    -5.1990084991292642, 4.0019830017413716, -9.9998658149487163,
    5.6002683701023734, -16.399981840027902, 7.2000363199438731,
    -24.399997542314779, 8.8000049153699571, -33.999999667388131,
    10.400000665223061, -45.199999954985429, 12.000000090028237,
    -57.999999993907373, 13.600000012184093, -72.399999999174753,
    15.200000001649045, -1.0, -1.0, 1.0, 1.0 };

  static const real_T d_a[16] = { 0.13533528323661148, 3.8141507400383685E-5,
    -0.09667244669717949, 0.0, -0.43233232473224048, 0.99979851413824794,
    0.51068068776924647, 0.0, -0.00017057400661253083, -7.949498279156195E-8,
    1.0002014858617532, 0.0, 0.0, 0.0, 0.0, 1.0 };

  static const real_T c_a[4] = { -58.775954131573386, -101.18805268783535,
    -12.268129756687021, 0.0 };

  static const real_T b_a[8] = { -4.620254134992216E-17, -0.013949226148762242,
    35.355336307537648, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const int8_T b_Mrows_0[34] = { 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 21,
    22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
    41, 42, 51, 52 };

  static const real_T a[8] = { -13.177009845729751, 27.16283383490083,
    43.408232032105914, 0.34899669706156022, 11.982717542914219,
    -30.901810818400556, -16.983629654686244, 0.41117101504513554 };

  static const real_T b_Kx[8] = { 0.00011342231913128645, -0.12663236628870159,
    -4.9962005728040526E-5, 8.0125214113942338, 1.5350039003145413E-5,
    -0.10264809407563126, -4.0499161584700704E-5, 6.4925214112293288 };

  static const real_T b_Linv[9] = { 0.10837060546229786, -1.265397718601619, 0.0,
    0.0, 1.4878212318678712, 0.0, 0.0, 0.0, 0.0005623413251903491 };

  static const real_T b_Hinv[9] = { 1.612975574370447, -1.8826855924926547, 0.0,
    -1.8826855924926547, 2.2136120179968297, 0.0, 0.0, 0.0,
    3.1622776601683797E-7 };

  static const real_T b_Ac[102] = { -0.90826822658930306, -2.4146525111110129,
    -4.0019830017413716, -5.6002683701023734, -7.2000363199438731,
    -8.8000049153699571, -10.400000665223061, -12.000000090028237,
    -13.600000012184093, -15.200000001649045, -0.34586588670534729,
    0.90826822658930306, -1.9926737444444755, 2.4146525111110129,
    -5.1990084991292642, 4.0019830017413716, -9.9998658149487163,
    5.6002683701023734, -16.399981840027902, 7.2000363199438731,
    -24.399997542314779, 8.8000049153699571, -33.999999667388131,
    10.400000665223061, -45.199999954985429, 12.000000090028237,
    -57.999999993907373, 13.600000012184093, -72.399999999174753,
    15.200000001649045, -1.0, -1.0, 1.0, 1.0, -0.0, -0.90826822658930306,
    -2.4146525111110129, -4.0019830017413716, -5.6002683701023734,
    -7.2000363199438731, -8.8000049153699571, -10.400000665223061,
    -12.000000090028237, -13.600000012184093, 0.0, 0.0, -0.34586588670534729,
    0.90826822658930306, -1.9926737444444755, 2.4146525111110129,
    -5.1990084991292642, 4.0019830017413716, -9.9998658149487163,
    5.6002683701023734, -16.399981840027902, 7.2000363199438731,
    -24.399997542314779, 8.8000049153699571, -33.999999667388131,
    10.400000665223061, -45.199999954985429, 12.000000090028237,
    -57.999999993907373, 13.600000012184093, -0.0, -1.0, 0.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0 };

  static const real_T b_Kr[40] = { 0.0, -0.090826822658930315, 0.0,
    -0.24146525111110129, 0.0, -0.40019830017413716, 0.0, -0.56002683701023737,
    0.0, -0.72000363199438733, 0.0, -0.88000049153699578, 0.0,
    -1.0400000665223061, 0.0, -1.2000000090028238, 0.0, -1.3600000012184095, 0.0,
    -1.5200000001649046, -0.0, -0.0, 0.0, -0.090826822658930315, 0.0,
    -0.24146525111110129, 0.0, -0.40019830017413716, 0.0, -0.56002683701023737,
    0.0, -0.72000363199438733, 0.0, -0.88000049153699578, 0.0,
    -1.0400000665223061, 0.0, -1.2000000090028238, 0.0, -1.3600000012184095 };

  static const real_T b_Kv[44] = { -9.16712296476263E-19, 0.0,
    -9.6082584107482955E-19, 0.0, -9.0208634887205852E-19, 0.0,
    -8.21829325931044E-19, 0.0, -7.1991079458528791E-19, 0.0,
    -6.4184776082151158E-19, 0.0, -5.4643790504119014E-19, 0.0,
    -4.336808703598214E-19, 0.0, -3.4911309972333667E-19, 0.0,
    -2.0599841279459456E-19, 0.0, 0.0, 0.0, -7.450805640770756E-19, 0.0,
    -7.9311324879950636E-19, 0.0, -7.54827428280235E-19, 0.0,
    -6.9608793607746405E-19, 0.0, -6.1583091313644939E-19, 0.0,
    -5.5511206434961229E-19, 0.0, -4.77049030585836E-19, 0.0,
    -3.8163917480551443E-19, 0.0, -3.1008182268306464E-19, 0.0,
    -1.8431436948766103E-19, 0.0, 0.0, 0.0 };

  /* Reset subsysRan breadcrumbs */
  srClearBC(speedgoat_target_model_2021b_DW.Forward_SubsysRanBC);

  /* Reset subsysRan breadcrumbs */
  srClearBC(speedgoat_target_model_2021b_DW.Reverse_SubsysRanBC);

  /* RateTransition generated from: '<S4>/DataTypeConversion_L0' */
  rtw_slrealtime_mutex_lock
    (speedgoat_target_model_2021b_DW.default_spacing_d0_SEMAPHORE);
  speedgoat_target_model_2021b_DW.default_spacing_RDBuf =
    speedgoat_target_model_2021b_DW.default_spacing_LstBufWR;
  rtw_slrealtime_mutex_unlock
    (speedgoat_target_model_2021b_DW.default_spacing_d0_SEMAPHORE);
  switch (speedgoat_target_model_2021b_DW.default_spacing_RDBuf) {
   case 0:
    /* RateTransition generated from: '<S4>/DataTypeConversion_L0' */
    speedgoat_target_model_2021b_B.default_spacing =
      speedgoat_target_model_2021b_DW.default_spacing_Buf0;
    break;

   case 1:
    /* RateTransition generated from: '<S4>/DataTypeConversion_L0' */
    speedgoat_target_model_2021b_B.default_spacing =
      speedgoat_target_model_2021b_DW.default_spacing_Buf1;
    break;

   case 2:
    /* RateTransition generated from: '<S4>/DataTypeConversion_L0' */
    speedgoat_target_model_2021b_B.default_spacing =
      speedgoat_target_model_2021b_DW.default_spacing_Buf2;
    break;
  }

  /* End of RateTransition generated from: '<S4>/DataTypeConversion_L0' */

  /* MATLAB Function: '<S4>/DataTypeConversion_L0' */
  speed_DataTypeConversion_L0(speedgoat_target_model_2021b_B.default_spacing,
    &speedgoat_target_model_2021b_B.sf_DataTypeConversion_L0);

  /* MATLAB Function: '<S4>/DataTypeConversion_amax' incorporates:
   *  Constant: '<S4>/Maximum longitudinal acceleration constant'
   */
  spe_DataTypeConversion_amax
    (speedgoat_target_model_2021_cal->ACCsystem_MaxAcceleration,
     &speedgoat_target_model_2021b_B.sf_DataTypeConversion_amax);

  /* MATLAB Function: '<S4>/DataTypeConversion_amin' incorporates:
   *  Constant: '<S4>/Minimum longitudinal acceleration constant'
   */
  spe_DataTypeConversion_amax
    (speedgoat_target_model_2021_cal->ACCsystem_MinAcceleration,
     &speedgoat_target_model_2021b_B.sf_DataTypeConversion_amin);

  /* MATLAB Function: '<S4>/DataTypeConversion_optsgn' incorporates:
   *  Constant: '<S4>/Enable optimization constant'
   */
  speedgoat_target_model_2021b_B.y =
    ((speedgoat_target_model_2021_cal->Enableoptimizationconstant_Valu > 0.0) ||
     (speedgoat_target_model_2021_cal->Enableoptimizationconstant_Valu < 0.0));

  /* Math: '<S23>/Math Function' incorporates:
   *  Constant: '<S22>/y.wt_zero'
   */
  speedgoat_target_model_2021b_B.MathFunction[0] =
    speedgoat_target_model_2021_cal->ywt_zero_Value[0];
  speedgoat_target_model_2021b_B.MathFunction[1] =
    speedgoat_target_model_2021_cal->ywt_zero_Value[1];

  /* Math: '<S23>/Math Function1' incorporates:
   *  Constant: '<S22>/u.wt_zero'
   */
  speedgoat_target_model_2021b_B.MathFunction1 =
    speedgoat_target_model_2021_cal->uwt_zero_Value;

  /* Math: '<S23>/Math Function2' incorporates:
   *  Constant: '<S22>/du.wt_zero'
   */
  speedgoat_target_model_2021b_B.MathFunction2 =
    speedgoat_target_model_2021_cal->duwt_zero_Value;
  for (i = 0; i < 34; i++) {
    /* Memory: '<S23>/Memory' */
    speedgoat_target_model_2021b_B.Memory[i] =
      speedgoat_target_model_2021b_DW.Memory_PreviousInput[i];
  }

  /* Gain: '<S23>/umin_scale4' incorporates:
   *  Constant: '<S22>/E_zero'
   */
  speedgoat_target_model_2021b_B.umin_scale4 =
    speedgoat_target_model_2021_cal->umin_scale4_Gain *
    speedgoat_target_model_2021_cal->E_zero_Value;

  /* Reshape: '<S23>/Reshape' */
  speedgoat_target_model_2021b_B.Reshape =
    speedgoat_target_model_2021b_B.umin_scale4;

  /* Gain: '<S23>/ymin_scale1' incorporates:
   *  Constant: '<S22>/F_zero'
   */
  speedgoat_target_model_2021b_B.ymin_scale1[0] =
    speedgoat_target_model_2021_cal->ymin_scale1_Gain[0] *
    speedgoat_target_model_2021_cal->F_zero_Value[0];
  speedgoat_target_model_2021b_B.ymin_scale1[1] =
    speedgoat_target_model_2021_cal->ymin_scale1_Gain[1] *
    speedgoat_target_model_2021_cal->F_zero_Value[1];

  /* Gain: '<S23>/ymin_scale2' incorporates:
   *  Constant: '<S22>/S_zero'
   */
  speedgoat_target_model_2021b_B.ymin_scale2 =
    speedgoat_target_model_2021_cal->ymin_scale2_Gain *
    speedgoat_target_model_2021_cal->S_zero_Value;

  /* Reshape: '<S23>/Reshape2' */
  speedgoat_target_model_2021b_B.Reshape2 =
    speedgoat_target_model_2021b_B.ymin_scale2;

  /* Reshape: '<S23>/Reshape1' incorporates:
   *  Gain: '<S23>/ymin_scale1'
   */
  speedgoat_target_model_2021b_B.Reshape1[0] =
    speedgoat_target_model_2021b_B.ymin_scale1[0];

  /* Reshape: '<S23>/Reshape3' incorporates:
   *  Math: '<S23>/Math Function'
   */
  speedgoat_target_model_2021b_B.Reshape3[0] =
    speedgoat_target_model_2021b_B.MathFunction[0];

  /* Reshape: '<S23>/Reshape1' incorporates:
   *  Gain: '<S23>/ymin_scale1'
   */
  speedgoat_target_model_2021b_B.Reshape1[1] =
    speedgoat_target_model_2021b_B.ymin_scale1[1];

  /* Reshape: '<S23>/Reshape3' incorporates:
   *  Math: '<S23>/Math Function'
   */
  speedgoat_target_model_2021b_B.Reshape3[1] =
    speedgoat_target_model_2021b_B.MathFunction[1];

  /* Reshape: '<S23>/Reshape4' */
  speedgoat_target_model_2021b_B.Reshape4 =
    speedgoat_target_model_2021b_B.MathFunction1;

  /* Reshape: '<S23>/Reshape5' */
  speedgoat_target_model_2021b_B.Reshape5 =
    speedgoat_target_model_2021b_B.MathFunction2;

  /* Gain: '<S23>/ext.mv_scale' incorporates:
   *  Constant: '<S22>/ext.mv_zero'
   */
  speedgoat_target_model_2021b_B.extmv_scale =
    speedgoat_target_model_2021_cal->extmv_scale_Gain *
    speedgoat_target_model_2021_cal->extmv_zero_Value;

  /* Gain: '<S23>/ext.mv_scale1' incorporates:
   *  Constant: '<S22>/mv.target_zero'
   */
  speedgoat_target_model_2021b_B.extmv_scale1 =
    speedgoat_target_model_2021_cal->extmv_scale1_Gain *
    speedgoat_target_model_2021_cal->mvtarget_zero_Value;

  /* UnitDelay: '<S23>/last_mv' */
  speedgoat_target_model_2021b_B.last_mv =
    speedgoat_target_model_2021b_DW.last_mv_DSTATE;

  /* Memory: '<S23>/last_x' */
  speedgoat_target_model_2021b_B.last_x[0] =
    speedgoat_target_model_2021b_DW.last_x_PreviousInput[0];
  speedgoat_target_model_2021b_B.last_x[1] =
    speedgoat_target_model_2021b_DW.last_x_PreviousInput[1];
  speedgoat_target_model_2021b_B.last_x[2] =
    speedgoat_target_model_2021b_DW.last_x_PreviousInput[2];
  speedgoat_target_model_2021b_B.last_x[3] =
    speedgoat_target_model_2021b_DW.last_x_PreviousInput[3];

  /* RateTransition generated from: '<S43>/optimizer' */
  speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport3_RdBufI = static_cast<
    int8_T>(speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport3_RdBufI == 0);
  ii = speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport3_RdBufI << 1;

  /* RateTransition generated from: '<S43>/optimizer' */
  speedgoat_target_model_2021b_B.TmpRTBAtoptimizerInport3[0] =
    speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport3_Buf[ii];
  speedgoat_target_model_2021b_B.TmpRTBAtoptimizerInport3[1] =
    speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport3_Buf[ii + 1];

  /* RateTransition generated from: '<S43>/optimizer' */
  speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport5_RdBufI = static_cast<
    int8_T>(speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport5_RdBufI == 0);

  /* RateTransition generated from: '<S43>/optimizer' */
  speedgoat_target_model_2021b_B.TmpRTBAtoptimizerInport5 =
    speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport5_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport5_RdBufI];

  /* RateTransition generated from: '<S43>/optimizer' */
  rtw_slrealtime_mutex_lock
    (speedgoat_target_model_2021b_DW.min_velocity_d0_SEMAPHORE);
  speedgoat_target_model_2021b_DW.min_velocity_RDBuf =
    speedgoat_target_model_2021b_DW.min_velocity_LstBufWR;
  rtw_slrealtime_mutex_unlock
    (speedgoat_target_model_2021b_DW.min_velocity_d0_SEMAPHORE);
  switch (speedgoat_target_model_2021b_DW.min_velocity_RDBuf) {
   case 0:
    /* RateTransition generated from: '<S43>/optimizer' */
    speedgoat_target_model_2021b_B.min_velocity[0] =
      speedgoat_target_model_2021b_DW.min_velocity_Buf0[0];
    speedgoat_target_model_2021b_B.min_velocity[1] =
      speedgoat_target_model_2021b_DW.min_velocity_Buf0[1];
    break;

   case 1:
    /* RateTransition generated from: '<S43>/optimizer' */
    speedgoat_target_model_2021b_B.min_velocity[0] =
      speedgoat_target_model_2021b_DW.min_velocity_Buf1[0];
    speedgoat_target_model_2021b_B.min_velocity[1] =
      speedgoat_target_model_2021b_DW.min_velocity_Buf1[1];
    break;

   case 2:
    /* RateTransition generated from: '<S43>/optimizer' */
    speedgoat_target_model_2021b_B.min_velocity[0] =
      speedgoat_target_model_2021b_DW.min_velocity_Buf2[0];
    speedgoat_target_model_2021b_B.min_velocity[1] =
      speedgoat_target_model_2021b_DW.min_velocity_Buf2[1];
    break;
  }

  /* End of RateTransition generated from: '<S43>/optimizer' */

  /* RateTransition generated from: '<S4>/DataTypeConversion_vset' */
  speedgoat_target_model_2021b_DW.TmpRTBAtDataTypeConversion_vs_n = static_cast<
    int8_T>(speedgoat_target_model_2021b_DW.TmpRTBAtDataTypeConversion_vs_n == 0);

  /* RateTransition generated from: '<S4>/DataTypeConversion_vset' */
  speedgoat_target_model_2021b_B.TmpRTBAtDataTypeConversion_vset =
    speedgoat_target_model_2021b_DW.TmpRTBAtDataTypeConversion_vset[speedgoat_target_model_2021b_DW.TmpRTBAtDataTypeConversion_vs_n];

  /* SignalConversion generated from: '<S44>/ SFunction ' incorporates:
   *  MATLAB Function: '<S43>/optimizer'
   */
  speedgoat_target_model_2021b_B.TmpSignalConversionAtSFunctionI[0] =
    speedgoat_target_model_2021b_B.sf_DataTypeConversion_L0.y;
  speedgoat_target_model_2021b_B.TmpSignalConversionAtSFunctionI[1] =
    speedgoat_target_model_2021b_B.TmpRTBAtDataTypeConversion_vset;

  /* SignalConversion generated from: '<S44>/ SFunction ' incorporates:
   *  Constant: '<S4>/Maximum velocity constant'
   *  Constant: '<S4>/Unconstrained'
   *  MATLAB Function: '<S43>/optimizer'
   */
  speedgoat_target_model_2021b_B.TmpSignalConversionAtSFunctio_m[0] =
    speedgoat_target_model_2021_cal->Unconstrained_Value;
  speedgoat_target_model_2021b_B.TmpSignalConversionAtSFunctio_m[1] =
    speedgoat_target_model_2021_cal->ACCsystem_MaxVelocity;

  /* MATLAB Function: '<S43>/optimizer' */
  std::memset(&vseq[0], 0, 22U * sizeof(real_T));
  for (i = 0; i < 11; i++) {
    vseq[(i << 1) + 1] = 1.0;
  }

  for (i = 0; i < 10; i++) {
    rseq[i << 1] =
      speedgoat_target_model_2021b_B.TmpSignalConversionAtSFunctionI[0] * 0.02 -
      0.093699999999999992;
    rseq[(i << 1) + 1] =
      speedgoat_target_model_2021b_B.TmpSignalConversionAtSFunctionI[1] * 0.02;
  }

  for (i = 0; i < 11; i++) {
    vseq[i << 1] = speedgoat_target_model_RMDscale *
      speedgoat_target_model_2021b_B.TmpRTBAtoptimizerInport5;
  }

  v_idx_0 = vseq[0];
  v_idx_1 = vseq[1];
  old_u = speedgoat_target_model_2021b_B.last_mv;
  xk[0] = speedgoat_target_model_2021b_B.last_x[0];
  xk[1] = speedgoat_target_model_2021b_B.last_x[1];
  xk[2] = speedgoat_target_model_2021b_B.last_x[2];
  xk[3] = speedgoat_target_model_2021b_B.last_x[3];
  for (ii = 0; ii < 2; ii++) {
    umax_incr = f_a[ii] * xk[0];
    umax_incr += f_a[ii + 2] * xk[1];
    umax_incr += f_a[ii + 4] * xk[2];
    umax_incr += f_a[ii + 6] * xk[3];
    umin_incr = 0.0 * v_idx_0;
    umin_incr += 0.0 * v_idx_1;
    y_innov[ii] = (speedgoat_target_model_2021b_B.TmpRTBAtoptimizerInport3[ii] *
                   0.02 - (-0.093699999999999992 * static_cast<real_T>(ii) +
      0.093699999999999992)) - (umax_incr + umin_incr);
  }

  for (ii = 0; ii < 4; ii++) {
    umin_incr = e_a[ii] * y_innov[0];
    umin_incr += e_a[ii + 4] * y_innov[1];
    speedgoat_target_model_2021b_B.xest[ii] = xk[ii] + umin_incr;
  }

  speedgoat_target_model_2021b_B.status = 1.0;
  std::memset(&speedgoat_target_model_2021b_B.useq[0], 0, 11U * sizeof(real_T));
  for (i = 0; i < 34; i++) {
    speedgoat_target_model_2021b_B.iAout[i] = false;
  }

  umin_incr = rt_roundd_snf(speedgoat_target_model_2021b_B.y);
  if (umin_incr < 2.147483648E+9) {
    if (umin_incr >= -2.147483648E+9) {
      ii = static_cast<int32_T>(umin_incr);
    } else {
      ii = MIN_int32_T;
    }
  } else {
    ii = MAX_int32_T;
  }

  umax_incr_flag = (ii != 0);
  if (umax_incr_flag) {
    old_u = speedgoat_target_model_2021b_B.last_mv;
    for (i = 0; i < 11; i++) {
      speedgoat_target_model_2021b_B.useq[i] =
        speedgoat_target_model_2021b_B.last_mv;
    }
  } else {
    boolean_T umin_incr_flag;
    slrealtime::cblas_dgemv(102, 111, 34, 22, 1.0, (real_T *)&b_Mv[0], 34,
      &vseq[0], 1, 0.0, &a__1[0], 1);
    for (ii = 0; ii < 34; ii++) {
      umin_incr = b_Mx[ii] * speedgoat_target_model_2021b_B.xest[0];
      umin_incr += b_Mx[ii + 34] * speedgoat_target_model_2021b_B.xest[1];
      umin_incr += b_Mx[ii + 68] * speedgoat_target_model_2021b_B.xest[2];
      umin_incr += b_Mx[ii + 102] * speedgoat_target_model_2021b_B.xest[3];
      Bc[ii] = -(((b_Mlim_0[ii] + umin_incr) + b_Mu1[ii] *
                  speedgoat_target_model_2021b_B.last_mv) + a__1[ii]);
    }

    ymax_incr_flag[0] = false;
    ymax_incr[0] = 0.0;
    ymin_incr_flag[0] = false;
    ymin_incr[0] = 0.0;
    ymax_incr_flag[1] = false;
    ymax_incr[1] = 0.0;
    ymin_incr_flag[1] = false;
    ymin_incr[1] = 0.0;
    umax_incr = 0.0;
    umin_incr_flag = false;
    umin_incr = 0.0;
    for (i = 0; i < 34; i++) {
      real_T b_Mlim;
      int8_T b_Mrows;
      Bc_0 = Bc[i];
      b_Mlim = b_Mlim_0[i];
      b_Mrows = b_Mrows_0[i];
      if (b_Mrows <= 20) {
        boolean_T b_Del_Save_Flag0;
        ii = (b_Mrows - (((b_Mrows - 1) >> 1) << 1)) - 1;
        b_Del_Save_Flag0 = ymax_incr_flag[ii];
        if (!ymax_incr_flag[ii]) {
          b_Mlim = -(0.02 *
                     speedgoat_target_model_2021b_B.TmpSignalConversionAtSFunctio_m
                     [ii] - (-0.093699999999999992 * static_cast<real_T>(ii) +
                             0.093699999999999992)) - (-b_Mlim);
          b_Del_Save_Flag0 = true;
        } else {
          b_Mlim = ymax_incr[ii];
        }

        ymax_incr[ii] = b_Mlim;
        ymax_incr_flag[ii] = b_Del_Save_Flag0;
        Bc_0 += b_Mlim;
      } else if (b_Mrows <= 40) {
        boolean_T b_Del_Save_Flag0;
        ii = (b_Mrows - (((b_Mrows - 21) >> 1) << 1)) - 21;
        b_Del_Save_Flag0 = ymin_incr_flag[ii];
        if (!ymin_incr_flag[ii]) {
          b_Mlim = (0.02 * speedgoat_target_model_2021b_B.min_velocity[ii] -
                    (-0.093699999999999992 * static_cast<real_T>(ii) +
                     0.093699999999999992)) - (-b_Mlim);
          b_Del_Save_Flag0 = true;
        } else {
          b_Mlim = ymin_incr[ii];
        }

        ymin_incr[ii] = b_Mlim;
        ymin_incr_flag[ii] = b_Del_Save_Flag0;
        Bc_0 += b_Mlim;
      } else if (b_Mrows <= 50) {
        if (!umax_incr_flag) {
          umax_incr = -(speedgoat_target_model_RMVscale *
                        speedgoat_target_model_2021b_B.sf_DataTypeConversion_amax.y)
            - (-b_Mlim);
          umax_incr_flag = true;
        }

        Bc_0 += umax_incr;
      } else {
        if (!umin_incr_flag) {
          umin_incr = speedgoat_target_model_RMVscale *
            speedgoat_target_model_2021b_B.sf_DataTypeConversion_amin.y -
            (-b_Mlim);
          umin_incr_flag = true;
        }

        Bc_0 += umin_incr;
      }

      Bc[i] = Bc_0;
    }

    f[0] = 0.0;
    f[1] = 0.0;
    f[2] = 0.0;
    for (i = 0; i < 2; i++) {
      umax_incr = b_Kx[i << 2] * speedgoat_target_model_2021b_B.xest[0];
      umax_incr += b_Kx[(i << 2) + 1] * speedgoat_target_model_2021b_B.xest[1];
      umax_incr += b_Kx[(i << 2) + 2] * speedgoat_target_model_2021b_B.xest[2];
      umax_incr += b_Kx[(i << 2) + 3] * speedgoat_target_model_2021b_B.xest[3];
      umin_incr = 0.0;
      for (ii = 0; ii < 20; ii++) {
        umin_incr += b_Kr[20 * i + ii] * rseq[ii];
      }

      Bc_0 = 0.0;
      for (ii = 0; ii < 22; ii++) {
        Bc_0 += b_Kv[22 * i + ii] * vseq[ii];
      }

      f[i] = ((-12.72837101576097 * static_cast<real_T>(i) + 85.147499758214479)
              * speedgoat_target_model_2021b_B.last_mv + (umax_incr + umin_incr))
        + Bc_0;
    }

    for (i = 0; i < 34; i++) {
      iAnew[i] = speedgoat_target_model_2021b_B.Memory[i];
    }

    speedgoat_target_model_2_qpkwik(b_Linv, b_Hinv, f, b_Ac, Bc, iAnew, 148,
      1.0E-6, zopt, a__1, &umax_incr);
    speedgoat_target_model_2021b_B.status = umax_incr;
    for (i = 0; i < 34; i++) {
      speedgoat_target_model_2021b_B.iAout[i] = (iAnew[i] != 0);
    }

    umin_incr = rt_roundd_snf(umax_incr);
    if (umin_incr < 2.147483648E+9) {
      if (umin_incr >= -2.147483648E+9) {
        ii = static_cast<int32_T>(umin_incr);
      } else {
        ii = MIN_int32_T;
      }
    } else {
      ii = MAX_int32_T;
    }

    umin_incr = rt_roundd_snf(umax_incr);
    if (umin_incr < 2.147483648E+9) {
      if (umin_incr >= -2.147483648E+9) {
        i = static_cast<int32_T>(umin_incr);
      } else {
        i = MIN_int32_T;
      }
    } else {
      i = MAX_int32_T;
    }

    if ((ii < 0) || (i == 0)) {
      zopt[0] = 0.0;
    }

    old_u += zopt[0];
  }

  speedgoat_target_model_2021b_B.cost = 0.0;
  speedgoat_target_model_2021b_B.u = old_u;
  std::memset(&speedgoat_target_model_2021b_B.yseq[0], 0, 22U * sizeof(real_T));
  std::memset(&speedgoat_target_model_2021b_B.xseq[0], 0, 44U * sizeof(real_T));
  for (ii = 0; ii < 4; ii++) {
    umin_incr = d_a[ii] * xk[0];
    umin_incr += d_a[ii + 4] * xk[1];
    umin_incr += d_a[ii + 8] * xk[2];
    umin_incr += d_a[ii + 12] * xk[3];
    umax_incr = c_a[ii] * old_u + umin_incr;
    Bc_0 = b_a[ii] * v_idx_0;
    umin_incr = a[ii] * y_innov[0];
    Bc_0 += 0.0 * v_idx_1;
    umin_incr += a[ii + 4] * y_innov[1];
    speedgoat_target_model_2021b_B.xk1[ii] = (umax_incr + Bc_0) + umin_incr;
  }

  /* Gain: '<S23>/umin_scale1' */
  speedgoat_target_model_2021b_B.umin_scale1 =
    speedgoat_target_model_2021_cal->umin_scale1_Gain *
    speedgoat_target_model_2021b_B.u;
  for (i = 0; i < 11; i++) {
    /* Gain: '<S23>/umin_scale3' */
    speedgoat_target_model_2021b_B.umin_scale3[i] =
      speedgoat_target_model_2021_cal->umin_scale3_Gain[i] *
      speedgoat_target_model_2021b_B.useq[i];
  }

  for (i = 0; i < 22; i++) {
    /* Gain: '<S23>/umin_scale5' */
    speedgoat_target_model_2021b_B.umin_scale5[i] =
      speedgoat_target_model_2021_cal->umin_scale5_Gain[i] *
      speedgoat_target_model_2021b_B.yseq[i];
  }

  /* Gain: '<S23>/umin_scale2' incorporates:
   *  Constant: '<S23>/constant'
   */
  speedgoat_target_model_2021b_B.umin_scale2 =
    speedgoat_target_model_2021_cal->umin_scale2_Gain *
    speedgoat_target_model_2021_cal->constant_Value;

  /* RateTransition generated from: '<S51>/Equal2' */
  speedgoat_target_model_2021b_DW.TmpRTBAtEqual2Inport1_RdBufIdx =
    static_cast<int8_T>
    (speedgoat_target_model_2021b_DW.TmpRTBAtEqual2Inport1_RdBufIdx == 0);

  /* RateTransition generated from: '<S51>/Equal2' */
  speedgoat_target_model_2021b_B.TmpRTBAtEqual2Inport1 =
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual2Inport1_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtEqual2Inport1_RdBufIdx];

  /* RateTransition generated from: '<S51>/Sign2' */
  speedgoat_target_model_2021b_DW.set_velocity1_RdBufIdx = static_cast<int8_T>
    (speedgoat_target_model_2021b_DW.set_velocity1_RdBufIdx == 0);

  /* RateTransition generated from: '<S51>/Sign2' */
  speedgoat_target_model_2021b_B.set_velocity1 =
    speedgoat_target_model_2021b_DW.set_velocity1_Buf[speedgoat_target_model_2021b_DW.set_velocity1_RdBufIdx];

  /* Signum: '<S51>/Sign2' */
  v_idx_0 = speedgoat_target_model_2021b_B.set_velocity1;
  if (v_idx_0 < 0.0) {
    /* Signum: '<S51>/Sign2' */
    speedgoat_target_model_2021b_B.Sign2 = -1.0;
  } else if (v_idx_0 > 0.0) {
    /* Signum: '<S51>/Sign2' */
    speedgoat_target_model_2021b_B.Sign2 = 1.0;
  } else if (v_idx_0 == 0.0) {
    /* Signum: '<S51>/Sign2' */
    speedgoat_target_model_2021b_B.Sign2 = 0.0;
  } else {
    /* Signum: '<S51>/Sign2' */
    speedgoat_target_model_2021b_B.Sign2 = (rtNaN);
  }

  /* End of Signum: '<S51>/Sign2' */

  /* RelationalOperator: '<S51>/Equal2' */
  speedgoat_target_model_2021b_B.Equal2 =
    (speedgoat_target_model_2021b_B.TmpRTBAtEqual2Inport1 !=
     speedgoat_target_model_2021b_B.Sign2);

  /* RateTransition generated from: '<S51>/Equal6' */
  speedgoat_target_model_2021b_DW.set_velocity_RdBufIdx = static_cast<int8_T>
    (speedgoat_target_model_2021b_DW.set_velocity_RdBufIdx == 0);

  /* RateTransition generated from: '<S51>/Equal6' */
  speedgoat_target_model_2021b_B.set_velocity =
    speedgoat_target_model_2021b_DW.set_velocity_Buf[speedgoat_target_model_2021b_DW.set_velocity_RdBufIdx];

  /* RelationalOperator: '<S51>/Equal6' incorporates:
   *  Constant: '<S51>/Constant4'
   */
  speedgoat_target_model_2021b_B.Equal6 =
    (speedgoat_target_model_2021b_B.set_velocity !=
     speedgoat_target_model_2021_cal->Constant4_Value);

  /* Logic: '<S51>/AND3' */
  speedgoat_target_model_2021b_B.AND3 = (speedgoat_target_model_2021b_B.Equal2 &&
    speedgoat_target_model_2021b_B.Equal6);

  /* RateTransition generated from: '<S51>/Equal1' */
  speedgoat_target_model_2021b_DW.TmpRTBAtEqual1Inport1_RdBufIdx = static_cast<
    int8_T>(speedgoat_target_model_2021b_DW.TmpRTBAtEqual1Inport1_RdBufIdx == 0);

  /* RateTransition generated from: '<S51>/Equal1' */
  speedgoat_target_model_2021b_B.TmpRTBAtEqual1Inport1 =
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual1Inport1_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtEqual1Inport1_RdBufIdx];

  /* RateTransition generated from: '<S51>/Sign1' */
  speedgoat_target_model_2021b_DW.TmpRTBAtSign1Inport1_RdBufIdx =
    static_cast<int8_T>
    (speedgoat_target_model_2021b_DW.TmpRTBAtSign1Inport1_RdBufIdx == 0);

  /* RateTransition generated from: '<S51>/Sign1' */
  speedgoat_target_model_2021b_B.TmpRTBAtSign1Inport1 =
    speedgoat_target_model_2021b_DW.TmpRTBAtSign1Inport1_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtSign1Inport1_RdBufIdx];

  /* Signum: '<S51>/Sign1' */
  v_idx_0 = speedgoat_target_model_2021b_B.TmpRTBAtSign1Inport1;
  if (v_idx_0 < 0.0) {
    /* Signum: '<S51>/Sign1' */
    speedgoat_target_model_2021b_B.Sign1 = -1.0;
  } else if (v_idx_0 > 0.0) {
    /* Signum: '<S51>/Sign1' */
    speedgoat_target_model_2021b_B.Sign1 = 1.0;
  } else if (v_idx_0 == 0.0) {
    /* Signum: '<S51>/Sign1' */
    speedgoat_target_model_2021b_B.Sign1 = 0.0;
  } else {
    /* Signum: '<S51>/Sign1' */
    speedgoat_target_model_2021b_B.Sign1 = (rtNaN);
  }

  /* End of Signum: '<S51>/Sign1' */

  /* RelationalOperator: '<S51>/Equal1' */
  speedgoat_target_model_2021b_B.Equal1 =
    (speedgoat_target_model_2021b_B.TmpRTBAtEqual1Inport1 !=
     speedgoat_target_model_2021b_B.Sign1);

  /* RateTransition generated from: '<S51>/Equal3' */
  speedgoat_target_model_2021b_DW.TmpRTBAtEqual3Inport1_RdBufIdx =
    static_cast<int8_T>
    (speedgoat_target_model_2021b_DW.TmpRTBAtEqual3Inport1_RdBufIdx == 0);

  /* RateTransition generated from: '<S51>/Equal3' */
  speedgoat_target_model_2021b_B.TmpRTBAtEqual3Inport1 =
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual3Inport1_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtEqual3Inport1_RdBufIdx];

  /* RelationalOperator: '<S51>/Equal3' incorporates:
   *  Constant: '<S51>/Constant1'
   */
  speedgoat_target_model_2021b_B.Equal3 =
    (speedgoat_target_model_2021b_B.TmpRTBAtEqual3Inport1 !=
     speedgoat_target_model_2021_cal->Constant1_Value_id);

  /* Logic: '<S51>/AND1' */
  speedgoat_target_model_2021b_B.AND1 = (speedgoat_target_model_2021b_B.Equal1 &&
    speedgoat_target_model_2021b_B.Equal3);

  /* Logic: '<S51>/OR' */
  speedgoat_target_model_2021b_B.OR_h = (speedgoat_target_model_2021b_B.AND3 ||
    speedgoat_target_model_2021b_B.AND1);

  /* Logic: '<S51>/NOT1' */
  speedgoat_target_model_2021b_B.NOT1 = !speedgoat_target_model_2021b_B.OR_h;

  /* Assertion: '<S51>/Assertion' */
  speedgoat_target_model_2021b_DW.Assertion_sltestCurrentResult =
    !speedgoat_target_model_2021b_B.NOT1;
  if (speedgoat_target_model_2021b_DW.Assertion_sltestFinalResult <
      speedgoat_target_model_2021b_DW.Assertion_sltestCurrentResult) {
    speedgoat_target_model_2021b_DW.Assertion_sltestFinalResult =
      speedgoat_target_model_2021b_DW.Assertion_sltestCurrentResult;
    speedgoat_target_model_2021b_DW.Assertion_sltestLastResultTime =
      speedgoat_target_model_2021b_M->Timing.taskTime0;
  }

  slTestLogAssessments(&speedgoat_target_model_2021b_DW.Assertion_sltestBlkInfo,
                       speedgoat_target_model_2021b_M->Timing.taskTime0);
  utAssert(speedgoat_target_model_2021b_B.NOT1);

  /* End of Assertion: '<S51>/Assertion' */

  /* RateTransition generated from: '<S51>/Equal5' */
  speedgoat_target_model_2021b_DW.TmpRTBAtEqual5Inport2_RdBufIdx =
    static_cast<int8_T>
    (speedgoat_target_model_2021b_DW.TmpRTBAtEqual5Inport2_RdBufIdx == 0);

  /* RateTransition generated from: '<S51>/Equal5' */
  speedgoat_target_model_2021b_B.TmpRTBAtEqual5Inport2 =
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual5Inport2_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtEqual5Inport2_RdBufIdx];

  /* RelationalOperator: '<S51>/Equal5' incorporates:
   *  Constant: '<S51>/Constant3'
   */
  speedgoat_target_model_2021b_B.Equal5 =
    (speedgoat_target_model_2021_cal->Constant3_Value !=
     speedgoat_target_model_2021b_B.TmpRTBAtEqual5Inport2);

  /* RateTransition generated from: '<S51>/Equal4' */
  speedgoat_target_model_2021b_DW.TmpRTBAtEqual4Inport2_RdBufIdx =
    static_cast<int8_T>
    (speedgoat_target_model_2021b_DW.TmpRTBAtEqual4Inport2_RdBufIdx == 0);

  /* RateTransition generated from: '<S51>/Equal4' */
  speedgoat_target_model_2021b_B.TmpRTBAtEqual4Inport2 =
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual4Inport2_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtEqual4Inport2_RdBufIdx];

  /* RelationalOperator: '<S51>/Equal4' incorporates:
   *  Constant: '<S51>/Constant2'
   */
  speedgoat_target_model_2021b_B.Equal4 =
    (speedgoat_target_model_2021_cal->Constant2_Value !=
     speedgoat_target_model_2021b_B.TmpRTBAtEqual4Inport2);

  /* Logic: '<S51>/AND2' */
  speedgoat_target_model_2021b_B.AND2 = (speedgoat_target_model_2021b_B.Equal5 &&
    speedgoat_target_model_2021b_B.Equal4);

  /* Logic: '<S51>/NOT' */
  speedgoat_target_model_2021b_B.NOT = !speedgoat_target_model_2021b_B.AND2;

  /* Assertion: '<S51>/Assertion1' */
  speedgoat_target_model_2021b_DW.Assertion1_sltestCurrentResult =
    !speedgoat_target_model_2021b_B.NOT;
  if (speedgoat_target_model_2021b_DW.Assertion1_sltestFinalResult <
      speedgoat_target_model_2021b_DW.Assertion1_sltestCurrentResult) {
    speedgoat_target_model_2021b_DW.Assertion1_sltestFinalResult =
      speedgoat_target_model_2021b_DW.Assertion1_sltestCurrentResult;
    speedgoat_target_model_2021b_DW.Assertion1_sltestLastResultTime =
      speedgoat_target_model_2021b_M->Timing.taskTime0;
  }

  slTestLogAssessments(&speedgoat_target_model_2021b_DW.Assertion1_sltestBlkInfo,
                       speedgoat_target_model_2021b_M->Timing.taskTime0);
  utAssert(speedgoat_target_model_2021b_B.NOT);

  /* End of Assertion: '<S51>/Assertion1' */

  /* RateTransition generated from: '<S49>/Multiply' */
  speedgoat_target_model_2021b_DW.TmpRTBAtMultiplyInport1_RdBufId = static_cast<
    int8_T>(speedgoat_target_model_2021b_DW.TmpRTBAtMultiplyInport1_RdBufId == 0);

  /* RateTransition generated from: '<S49>/Multiply' */
  speedgoat_target_model_2021b_B.TmpRTBAtMultiplyInport1 =
    speedgoat_target_model_2021b_DW.TmpRTBAtMultiplyInport1_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtMultiplyInport1_RdBufId];

  /* RateTransition generated from: '<S50>/Switch Case' */
  speedgoat_target_model_2021b_DW.TmpRTBAtSwitchCaseInport1_RdBuf = static_cast<
    int8_T>(speedgoat_target_model_2021b_DW.TmpRTBAtSwitchCaseInport1_RdBuf == 0);

  /* RateTransition generated from: '<S50>/Switch Case' */
  speedgoat_target_model_2021b_B.TmpRTBAtSwitchCaseInport1 =
    speedgoat_target_model_2021b_DW.TmpRTBAtSwitchCaseInport1_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtSwitchCaseInport1_RdBuf];

  /* RateTransition generated from: '<S8>/Minus' */
  speedgoat_target_model_2021b_DW.set_velocity_RdBufIdx_p = static_cast<int8_T>
    (speedgoat_target_model_2021b_DW.set_velocity_RdBufIdx_p == 0);

  /* RateTransition generated from: '<S8>/Minus' */
  speedgoat_target_model_2021b_B.set_velocity_c =
    speedgoat_target_model_2021b_DW.set_velocity_Buf_n[speedgoat_target_model_2021b_DW.set_velocity_RdBufIdx_p];

  /* RateTransition generated from: '<S8>/Minus' */
  speedgoat_target_model_2021b_DW.TmpRTBAtMinusInport2_RdBufIdx =
    static_cast<int8_T>
    (speedgoat_target_model_2021b_DW.TmpRTBAtMinusInport2_RdBufIdx == 0);

  /* RateTransition generated from: '<S8>/Minus' */
  speedgoat_target_model_2021b_B.TmpRTBAtMinusInport2 =
    speedgoat_target_model_2021b_DW.TmpRTBAtMinusInport2_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtMinusInport2_RdBufIdx];

  /* Sum: '<S8>/Minus' */
  speedgoat_target_model_2021b_B.error =
    speedgoat_target_model_2021b_B.set_velocity_c -
    speedgoat_target_model_2021b_B.TmpRTBAtMinusInport2;

  /* DataTypeConversion: '<S50>/Data Type Conversion' incorporates:
   *  Constant: '<S1>/reset'
   */
  speedgoat_target_model_2021b_B.DataTypeConversion =
    (speedgoat_target_model_2021_cal->reset_Value != 0.0);

  /* SwitchCase: '<S50>/Switch Case' */
  v_idx_0 = speedgoat_target_model_2021b_B.TmpRTBAtSwitchCaseInport1;
  if (v_idx_0 < 0.0) {
    umin_incr = std::ceil(v_idx_0);
  } else {
    umin_incr = std::floor(v_idx_0);
  }

  if (rtIsNaN(umin_incr) || rtIsInf(umin_incr)) {
    umin_incr = 0.0;
  } else {
    umin_incr = std::fmod(umin_incr, 4.294967296E+9);
  }

  switch (umin_incr < 0.0 ? -static_cast<int32_T>(static_cast<uint32_T>
           (-umin_incr)) : static_cast<int32_T>(static_cast<uint32_T>(umin_incr)))
  {
   case 1:
    /* Outputs for IfAction SubSystem: '<S50>/Forward' incorporates:
     *  ActionPort: '<S52>/Action Port'
     */
    /* Gain: '<S93>/Proportional Gain' */
    speedgoat_target_model_2021b_B.ProportionalGain_f =
      speedgoat_target_model_2021_cal->LongitudinalControllerStanley_n *
      speedgoat_target_model_2021b_B.error;

    /* DiscreteIntegrator: '<S88>/Integrator' */
    if (speedgoat_target_model_2021b_B.DataTypeConversion ||
        (speedgoat_target_model_2021b_DW.Integrator_PrevResetState_j != 0)) {
      speedgoat_target_model_2021b_DW.Integrator_DSTATE_k =
        speedgoat_target_model_2021_cal->PIForward_InitialConditionForIn;
    }

    /* DiscreteIntegrator: '<S88>/Integrator' */
    speedgoat_target_model_2021b_B.Integrator_b =
      speedgoat_target_model_2021b_DW.Integrator_DSTATE_k;

    /* Sum: '<S97>/Sum' */
    speedgoat_target_model_2021b_B.Sum_k =
      speedgoat_target_model_2021b_B.ProportionalGain_f +
      speedgoat_target_model_2021b_B.Integrator_b;

    /* Gain: '<S79>/ZeroGain' */
    speedgoat_target_model_2021b_B.ZeroGain_l =
      speedgoat_target_model_2021_cal->ZeroGain_Gain *
      speedgoat_target_model_2021b_B.Sum_k;

    /* DeadZone: '<S81>/DeadZone' */
    if (speedgoat_target_model_2021b_B.Sum_k >
        speedgoat_target_model_2021_cal->PIForward_UpperSaturationLimit) {
      /* DeadZone: '<S81>/DeadZone' */
      speedgoat_target_model_2021b_B.DeadZone_n =
        speedgoat_target_model_2021b_B.Sum_k -
        speedgoat_target_model_2021_cal->PIForward_UpperSaturationLimit;
    } else if (speedgoat_target_model_2021b_B.Sum_k >=
               speedgoat_target_model_2021_cal->PIForward_LowerSaturationLimit)
    {
      /* DeadZone: '<S81>/DeadZone' */
      speedgoat_target_model_2021b_B.DeadZone_n = 0.0;
    } else {
      /* DeadZone: '<S81>/DeadZone' */
      speedgoat_target_model_2021b_B.DeadZone_n =
        speedgoat_target_model_2021b_B.Sum_k -
        speedgoat_target_model_2021_cal->PIForward_LowerSaturationLimit;
    }

    /* End of DeadZone: '<S81>/DeadZone' */

    /* RelationalOperator: '<S79>/NotEqual' */
    speedgoat_target_model_2021b_B.NotEqual_b =
      (speedgoat_target_model_2021b_B.ZeroGain_l !=
       speedgoat_target_model_2021b_B.DeadZone_n);

    /* Signum: '<S79>/SignPreSat' */
    v_idx_0 = speedgoat_target_model_2021b_B.DeadZone_n;
    if (v_idx_0 < 0.0) {
      /* Signum: '<S79>/SignPreSat' */
      speedgoat_target_model_2021b_B.SignPreSat_m = -1.0;
    } else if (v_idx_0 > 0.0) {
      /* Signum: '<S79>/SignPreSat' */
      speedgoat_target_model_2021b_B.SignPreSat_m = 1.0;
    } else if (v_idx_0 == 0.0) {
      /* Signum: '<S79>/SignPreSat' */
      speedgoat_target_model_2021b_B.SignPreSat_m = 0.0;
    } else {
      /* Signum: '<S79>/SignPreSat' */
      speedgoat_target_model_2021b_B.SignPreSat_m = (rtNaN);
    }

    /* End of Signum: '<S79>/SignPreSat' */

    /* DataTypeConversion: '<S79>/DataTypeConv1' */
    umin_incr = std::floor(speedgoat_target_model_2021b_B.SignPreSat_m);
    if (rtIsNaN(umin_incr) || rtIsInf(umin_incr)) {
      umin_incr = 0.0;
    } else {
      umin_incr = std::fmod(umin_incr, 256.0);
    }

    /* DataTypeConversion: '<S79>/DataTypeConv1' */
    speedgoat_target_model_2021b_B.DataTypeConv1_p = static_cast<int8_T>
      (umin_incr < 0.0 ? static_cast<int32_T>(static_cast<int8_T>
        (-static_cast<int8_T>(static_cast<uint8_T>(-umin_incr)))) : static_cast<
       int32_T>(static_cast<int8_T>(static_cast<uint8_T>(umin_incr))));

    /* Gain: '<S85>/Integral Gain' */
    speedgoat_target_model_2021b_B.IntegralGain_a =
      speedgoat_target_model_2021_cal->LongitudinalControllerStanley_K *
      speedgoat_target_model_2021b_B.error;

    /* Signum: '<S79>/SignPreIntegrator' */
    v_idx_0 = speedgoat_target_model_2021b_B.IntegralGain_a;
    if (v_idx_0 < 0.0) {
      /* Signum: '<S79>/SignPreIntegrator' */
      speedgoat_target_model_2021b_B.SignPreIntegrator_d = -1.0;
    } else if (v_idx_0 > 0.0) {
      /* Signum: '<S79>/SignPreIntegrator' */
      speedgoat_target_model_2021b_B.SignPreIntegrator_d = 1.0;
    } else if (v_idx_0 == 0.0) {
      /* Signum: '<S79>/SignPreIntegrator' */
      speedgoat_target_model_2021b_B.SignPreIntegrator_d = 0.0;
    } else {
      /* Signum: '<S79>/SignPreIntegrator' */
      speedgoat_target_model_2021b_B.SignPreIntegrator_d = (rtNaN);
    }

    /* End of Signum: '<S79>/SignPreIntegrator' */

    /* DataTypeConversion: '<S79>/DataTypeConv2' */
    umin_incr = std::floor(speedgoat_target_model_2021b_B.SignPreIntegrator_d);
    if (rtIsNaN(umin_incr) || rtIsInf(umin_incr)) {
      umin_incr = 0.0;
    } else {
      umin_incr = std::fmod(umin_incr, 256.0);
    }

    /* DataTypeConversion: '<S79>/DataTypeConv2' */
    speedgoat_target_model_2021b_B.DataTypeConv2_l = static_cast<int8_T>
      (umin_incr < 0.0 ? static_cast<int32_T>(static_cast<int8_T>
        (-static_cast<int8_T>(static_cast<uint8_T>(-umin_incr)))) : static_cast<
       int32_T>(static_cast<int8_T>(static_cast<uint8_T>(umin_incr))));

    /* RelationalOperator: '<S79>/Equal1' */
    speedgoat_target_model_2021b_B.Equal1_a =
      (speedgoat_target_model_2021b_B.DataTypeConv1_p ==
       speedgoat_target_model_2021b_B.DataTypeConv2_l);

    /* Logic: '<S79>/AND3' */
    speedgoat_target_model_2021b_B.AND3_j =
      (speedgoat_target_model_2021b_B.NotEqual_b &&
       speedgoat_target_model_2021b_B.Equal1_a);

    /* Switch: '<S79>/Switch' */
    if (speedgoat_target_model_2021b_B.AND3_j) {
      /* Switch: '<S79>/Switch' incorporates:
       *  Constant: '<S79>/Constant1'
       */
      speedgoat_target_model_2021b_B.Switch_l =
        speedgoat_target_model_2021_cal->Constant1_Value;
    } else {
      /* Switch: '<S79>/Switch' */
      speedgoat_target_model_2021b_B.Switch_l =
        speedgoat_target_model_2021b_B.IntegralGain_a;
    }

    /* End of Switch: '<S79>/Switch' */

    /* Saturate: '<S95>/Saturation' */
    v_idx_0 = speedgoat_target_model_2021b_B.Sum_k;
    v_idx_1 = speedgoat_target_model_2021_cal->PIForward_LowerSaturationLimit;
    old_u = speedgoat_target_model_2021_cal->PIForward_UpperSaturationLimit;
    if (v_idx_0 > old_u) {
      /* Merge: '<S50>/Merge' */
      speedgoat_target_model_2021b_B.Merge = old_u;
    } else if (v_idx_0 < v_idx_1) {
      /* Merge: '<S50>/Merge' */
      speedgoat_target_model_2021b_B.Merge = v_idx_1;
    } else {
      /* Merge: '<S50>/Merge' */
      speedgoat_target_model_2021b_B.Merge = v_idx_0;
    }

    /* End of Saturate: '<S95>/Saturation' */

    /* Update for DiscreteIntegrator: '<S88>/Integrator' */
    speedgoat_target_model_2021b_DW.Integrator_DSTATE_k +=
      speedgoat_target_model_2021_cal->Integrator_gainval *
      speedgoat_target_model_2021b_B.Switch_l;
    speedgoat_target_model_2021b_DW.Integrator_PrevResetState_j =
      static_cast<int8_T>(speedgoat_target_model_2021b_B.DataTypeConversion);

    /* End of Outputs for SubSystem: '<S50>/Forward' */

    /* Update for IfAction SubSystem: '<S50>/Forward' incorporates:
     *  ActionPort: '<S52>/Action Port'
     */
    /* Update for SwitchCase: '<S50>/Switch Case' */
    srUpdateBC(speedgoat_target_model_2021b_DW.Forward_SubsysRanBC);

    /* End of Update for SubSystem: '<S50>/Forward' */
    break;

   case -1:
    /* Outputs for IfAction SubSystem: '<S50>/Reverse' incorporates:
     *  ActionPort: '<S53>/Action Port'
     */
    /* Gain: '<S144>/Proportional Gain' */
    speedgoat_target_model_2021b_B.ProportionalGain =
      speedgoat_target_model_2021_cal->LongitudinalControllerStanley_n *
      speedgoat_target_model_2021b_B.error;

    /* DiscreteIntegrator: '<S139>/Integrator' */
    if (speedgoat_target_model_2021b_B.DataTypeConversion ||
        (speedgoat_target_model_2021b_DW.Integrator_PrevResetState != 0)) {
      speedgoat_target_model_2021b_DW.Integrator_DSTATE =
        speedgoat_target_model_2021_cal->PIReverse_InitialConditionForIn;
    }

    /* DiscreteIntegrator: '<S139>/Integrator' */
    speedgoat_target_model_2021b_B.Integrator =
      speedgoat_target_model_2021b_DW.Integrator_DSTATE;

    /* Sum: '<S148>/Sum' */
    speedgoat_target_model_2021b_B.Sum_fg =
      speedgoat_target_model_2021b_B.ProportionalGain +
      speedgoat_target_model_2021b_B.Integrator;

    /* Gain: '<S130>/ZeroGain' */
    speedgoat_target_model_2021b_B.ZeroGain =
      speedgoat_target_model_2021_cal->ZeroGain_Gain_k *
      speedgoat_target_model_2021b_B.Sum_fg;

    /* DeadZone: '<S132>/DeadZone' */
    if (speedgoat_target_model_2021b_B.Sum_fg >
        speedgoat_target_model_2021_cal->PIReverse_UpperSaturationLimit) {
      /* DeadZone: '<S132>/DeadZone' */
      speedgoat_target_model_2021b_B.DeadZone =
        speedgoat_target_model_2021b_B.Sum_fg -
        speedgoat_target_model_2021_cal->PIReverse_UpperSaturationLimit;
    } else if (speedgoat_target_model_2021b_B.Sum_fg >=
               speedgoat_target_model_2021_cal->PIReverse_LowerSaturationLimit)
    {
      /* DeadZone: '<S132>/DeadZone' */
      speedgoat_target_model_2021b_B.DeadZone = 0.0;
    } else {
      /* DeadZone: '<S132>/DeadZone' */
      speedgoat_target_model_2021b_B.DeadZone =
        speedgoat_target_model_2021b_B.Sum_fg -
        speedgoat_target_model_2021_cal->PIReverse_LowerSaturationLimit;
    }

    /* End of DeadZone: '<S132>/DeadZone' */

    /* RelationalOperator: '<S130>/NotEqual' */
    speedgoat_target_model_2021b_B.NotEqual =
      (speedgoat_target_model_2021b_B.ZeroGain !=
       speedgoat_target_model_2021b_B.DeadZone);

    /* Signum: '<S130>/SignPreSat' */
    v_idx_0 = speedgoat_target_model_2021b_B.DeadZone;
    if (v_idx_0 < 0.0) {
      /* Signum: '<S130>/SignPreSat' */
      speedgoat_target_model_2021b_B.SignPreSat = -1.0;
    } else if (v_idx_0 > 0.0) {
      /* Signum: '<S130>/SignPreSat' */
      speedgoat_target_model_2021b_B.SignPreSat = 1.0;
    } else if (v_idx_0 == 0.0) {
      /* Signum: '<S130>/SignPreSat' */
      speedgoat_target_model_2021b_B.SignPreSat = 0.0;
    } else {
      /* Signum: '<S130>/SignPreSat' */
      speedgoat_target_model_2021b_B.SignPreSat = (rtNaN);
    }

    /* End of Signum: '<S130>/SignPreSat' */

    /* DataTypeConversion: '<S130>/DataTypeConv1' */
    umin_incr = std::floor(speedgoat_target_model_2021b_B.SignPreSat);
    if (rtIsNaN(umin_incr) || rtIsInf(umin_incr)) {
      umin_incr = 0.0;
    } else {
      umin_incr = std::fmod(umin_incr, 256.0);
    }

    /* DataTypeConversion: '<S130>/DataTypeConv1' */
    speedgoat_target_model_2021b_B.DataTypeConv1 = static_cast<int8_T>(umin_incr
      < 0.0 ? static_cast<int32_T>(static_cast<int8_T>(-static_cast<int8_T>(
      static_cast<uint8_T>(-umin_incr)))) : static_cast<int32_T>(static_cast<
      int8_T>(static_cast<uint8_T>(umin_incr))));

    /* Gain: '<S136>/Integral Gain' */
    speedgoat_target_model_2021b_B.IntegralGain =
      speedgoat_target_model_2021_cal->LongitudinalControllerStanley_K *
      speedgoat_target_model_2021b_B.error;

    /* Signum: '<S130>/SignPreIntegrator' */
    v_idx_0 = speedgoat_target_model_2021b_B.IntegralGain;
    if (v_idx_0 < 0.0) {
      /* Signum: '<S130>/SignPreIntegrator' */
      speedgoat_target_model_2021b_B.SignPreIntegrator = -1.0;
    } else if (v_idx_0 > 0.0) {
      /* Signum: '<S130>/SignPreIntegrator' */
      speedgoat_target_model_2021b_B.SignPreIntegrator = 1.0;
    } else if (v_idx_0 == 0.0) {
      /* Signum: '<S130>/SignPreIntegrator' */
      speedgoat_target_model_2021b_B.SignPreIntegrator = 0.0;
    } else {
      /* Signum: '<S130>/SignPreIntegrator' */
      speedgoat_target_model_2021b_B.SignPreIntegrator = (rtNaN);
    }

    /* End of Signum: '<S130>/SignPreIntegrator' */

    /* DataTypeConversion: '<S130>/DataTypeConv2' */
    umin_incr = std::floor(speedgoat_target_model_2021b_B.SignPreIntegrator);
    if (rtIsNaN(umin_incr) || rtIsInf(umin_incr)) {
      umin_incr = 0.0;
    } else {
      umin_incr = std::fmod(umin_incr, 256.0);
    }

    /* DataTypeConversion: '<S130>/DataTypeConv2' */
    speedgoat_target_model_2021b_B.DataTypeConv2 = static_cast<int8_T>(umin_incr
      < 0.0 ? static_cast<int32_T>(static_cast<int8_T>(-static_cast<int8_T>(
      static_cast<uint8_T>(-umin_incr)))) : static_cast<int32_T>(static_cast<
      int8_T>(static_cast<uint8_T>(umin_incr))));

    /* RelationalOperator: '<S130>/Equal1' */
    speedgoat_target_model_2021b_B.Equal1_c =
      (speedgoat_target_model_2021b_B.DataTypeConv1 ==
       speedgoat_target_model_2021b_B.DataTypeConv2);

    /* Logic: '<S130>/AND3' */
    speedgoat_target_model_2021b_B.AND3_b =
      (speedgoat_target_model_2021b_B.NotEqual &&
       speedgoat_target_model_2021b_B.Equal1_c);

    /* Switch: '<S130>/Switch' */
    if (speedgoat_target_model_2021b_B.AND3_b) {
      /* Switch: '<S130>/Switch' incorporates:
       *  Constant: '<S130>/Constant1'
       */
      speedgoat_target_model_2021b_B.Switch_o =
        speedgoat_target_model_2021_cal->Constant1_Value_f;
    } else {
      /* Switch: '<S130>/Switch' */
      speedgoat_target_model_2021b_B.Switch_o =
        speedgoat_target_model_2021b_B.IntegralGain;
    }

    /* End of Switch: '<S130>/Switch' */

    /* Saturate: '<S146>/Saturation' */
    v_idx_0 = speedgoat_target_model_2021b_B.Sum_fg;
    v_idx_1 = speedgoat_target_model_2021_cal->PIReverse_LowerSaturationLimit;
    old_u = speedgoat_target_model_2021_cal->PIReverse_UpperSaturationLimit;
    if (v_idx_0 > old_u) {
      /* Merge: '<S50>/Merge' */
      speedgoat_target_model_2021b_B.Merge = old_u;
    } else if (v_idx_0 < v_idx_1) {
      /* Merge: '<S50>/Merge' */
      speedgoat_target_model_2021b_B.Merge = v_idx_1;
    } else {
      /* Merge: '<S50>/Merge' */
      speedgoat_target_model_2021b_B.Merge = v_idx_0;
    }

    /* End of Saturate: '<S146>/Saturation' */

    /* Update for DiscreteIntegrator: '<S139>/Integrator' */
    speedgoat_target_model_2021b_DW.Integrator_DSTATE +=
      speedgoat_target_model_2021_cal->Integrator_gainval_b *
      speedgoat_target_model_2021b_B.Switch_o;
    speedgoat_target_model_2021b_DW.Integrator_PrevResetState =
      static_cast<int8_T>(speedgoat_target_model_2021b_B.DataTypeConversion);

    /* End of Outputs for SubSystem: '<S50>/Reverse' */

    /* Update for IfAction SubSystem: '<S50>/Reverse' incorporates:
     *  ActionPort: '<S53>/Action Port'
     */
    /* Update for SwitchCase: '<S50>/Switch Case' */
    srUpdateBC(speedgoat_target_model_2021b_DW.Reverse_SubsysRanBC);

    /* End of Update for SubSystem: '<S50>/Reverse' */
    break;
  }

  /* End of SwitchCase: '<S50>/Switch Case' */

  /* Product: '<S49>/Multiply' */
  speedgoat_target_model_2021b_B.Multiply_k =
    speedgoat_target_model_2021b_B.TmpRTBAtMultiplyInport1 *
    speedgoat_target_model_2021b_B.Merge;

  /* Switch: '<S49>/Switch' */
  if (speedgoat_target_model_2021b_B.Multiply_k >
      speedgoat_target_model_2021_cal->Switch_Threshold_k) {
    /* Switch: '<S49>/Switch' */
    speedgoat_target_model_2021b_B.Switch =
      speedgoat_target_model_2021b_B.Multiply_k;
  } else {
    /* Switch: '<S49>/Switch' incorporates:
     *  Constant: '<S49>/Constant'
     */
    speedgoat_target_model_2021b_B.Switch =
      speedgoat_target_model_2021_cal->Constant_Value_m;
  }

  /* End of Switch: '<S49>/Switch' */

  /* Switch: '<S49>/Switch1' */
  if (speedgoat_target_model_2021b_B.Multiply_k >
      speedgoat_target_model_2021_cal->Switch1_Threshold_g) {
    /* Switch: '<S49>/Switch1' incorporates:
     *  Constant: '<S49>/Constant'
     */
    speedgoat_target_model_2021b_B.Switch1 =
      speedgoat_target_model_2021_cal->Constant_Value_m;
  } else {
    /* Abs: '<S49>/Abs' */
    speedgoat_target_model_2021b_B.Abs = std::abs
      (speedgoat_target_model_2021b_B.Multiply_k);

    /* Switch: '<S49>/Switch1' */
    speedgoat_target_model_2021b_B.Switch1 = speedgoat_target_model_2021b_B.Abs;
  }

  /* End of Switch: '<S49>/Switch1' */

  /* RateTransition generated from: '<S49>/Switch1' */
  speedgoat_target_model_2021b_DW.TmpRTBAtSwitch1Outport1_WrBufId =
    static_cast<int8_T>
    (speedgoat_target_model_2021b_DW.TmpRTBAtSwitch1Outport1_WrBufId == 0);
  speedgoat_target_model_2021b_DW.TmpRTBAtSwitch1Outport1_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtSwitch1Outport1_WrBufId]
    = speedgoat_target_model_2021b_B.Switch1;

  /* RateTransition generated from: '<S49>/Switch' */
  speedgoat_target_model_2021b_DW.TmpRTBAtSwitchOutport1_WrBufIdx =
    static_cast<int8_T>
    (speedgoat_target_model_2021b_DW.TmpRTBAtSwitchOutport1_WrBufIdx == 0);
  speedgoat_target_model_2021b_DW.TmpRTBAtSwitchOutport1_Buf[speedgoat_target_model_2021b_DW.TmpRTBAtSwitchOutport1_WrBufIdx]
    = speedgoat_target_model_2021b_B.Switch;

  /* RateTransition generated from: '<S11>/Subsystem2' */
  speedgoat_target_model_2021b_DW.acceleration_WrBufIdx = static_cast<int8_T>
    (speedgoat_target_model_2021b_DW.acceleration_WrBufIdx == 0);
  speedgoat_target_model_2021b_DW.acceleration_Buf[speedgoat_target_model_2021b_DW.acceleration_WrBufIdx]
    = speedgoat_target_model_2021b_B.umin_scale1;

  /* Update for Memory: '<S23>/Memory' */
  for (i = 0; i < 34; i++) {
    speedgoat_target_model_2021b_DW.Memory_PreviousInput[i] =
      speedgoat_target_model_2021b_B.iAout[i];
  }

  /* End of Update for Memory: '<S23>/Memory' */

  /* Update for UnitDelay: '<S23>/last_mv' */
  speedgoat_target_model_2021b_DW.last_mv_DSTATE =
    speedgoat_target_model_2021b_B.u;

  /* Update for Memory: '<S23>/last_x' */
  speedgoat_target_model_2021b_DW.last_x_PreviousInput[0] =
    speedgoat_target_model_2021b_B.xk1[0];
  speedgoat_target_model_2021b_DW.last_x_PreviousInput[1] =
    speedgoat_target_model_2021b_B.xk1[1];
  speedgoat_target_model_2021b_DW.last_x_PreviousInput[2] =
    speedgoat_target_model_2021b_B.xk1[2];
  speedgoat_target_model_2021b_DW.last_x_PreviousInput[3] =
    speedgoat_target_model_2021b_B.xk1[3];

  /* Update absolute time */
  /* The "clockTick1" counts the number of times the code of this task has
   * been executed. The resolution of this integer timer is 1.0, which is the step size
   * of the task. Size of "clockTick1" ensures timer will not overflow during the
   * application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick1 and the high bits
   * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
   */
  speedgoat_target_model_2021b_M->Timing.clockTick1++;
  if (!speedgoat_target_model_2021b_M->Timing.clockTick1) {
    speedgoat_target_model_2021b_M->Timing.clockTickH1++;
  }
}

/* Model initialize function */
void speedgoat_target_model_2021b_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));
  (speedgoat_target_model_2021b_M)->Timing.TaskCounters.cLimit[0] = 1;
  (speedgoat_target_model_2021b_M)->Timing.TaskCounters.cLimit[1] = 10;
  speedgoat_target_model_2021b_M->Timing.stepSize0 = 0.1;

  /* block I/O */
  (void) std::memset((static_cast<void *>(&speedgoat_target_model_2021b_B)), 0,
                     sizeof(B_speedgoat_target_model_2021_T));

  /* states (dwork) */
  (void) std::memset(static_cast<void *>(&speedgoat_target_model_2021b_DW), 0,
                     sizeof(DW_speedgoat_target_model_202_T));

  /* Start for Constant: '<S1>/Approaching Speed' */
  speedgoat_target_model_2021b_B.ApproachingSpeed =
    speedgoat_target_model_2021_cal->ApproachingSpeed_Value;

  /* Start for Constant: '<S1>/Scene Speed' */
  speedgoat_target_model_2021b_B.SceneSpeed =
    speedgoat_target_model_2021_cal->SceneSpeed_Value;

  /* Start for Constant: '<S1>/Constant1' */
  speedgoat_target_model_2021b_B.Constant1 =
    speedgoat_target_model_2021_cal->Constant1_Value_i4;

  /* Start for Constant: '<S1>/RR Crossing' */
  speedgoat_target_model_2021b_B.RRCrossing =
    speedgoat_target_model_2021_cal->RRCrossing_Value;

  /* Start for Constant: '<S1>/Navigating Flag' */
  speedgoat_target_model_2021b_B.NavigatingFlag =
    speedgoat_target_model_2021_cal->NavigatingFlag_Value;

  /* Start for Constant: '<S1>/Clear2Go' */
  speedgoat_target_model_2021b_B.Clear2Go =
    speedgoat_target_model_2021_cal->Clear2Go_Value;

  /* Start for Constant: '<S1>/RoadSpeed' */
  speedgoat_target_model_2021b_B.RoadSpeed =
    speedgoat_target_model_2021_cal->RoadSpeed_Value;

  /* Start for Constant: '<S1>/WaitTime' */
  speedgoat_target_model_2021b_B.WaitTime =
    speedgoat_target_model_2021_cal->WaitTime_Value;

  /* Start for Constant: '<S1>/MinSpeed' */
  speedgoat_target_model_2021b_B.MinSpeed =
    speedgoat_target_model_2021_cal->MinSpeed_Value;

  /* Start for Constant: '<S1>/LaneAvailable' */
  speedgoat_target_model_2021b_B.LaneAvailable =
    speedgoat_target_model_2021_cal->LaneAvailable_Value;

  /* Start for Constant: '<S1>/LaneChanged' */
  speedgoat_target_model_2021b_B.LaneChanged =
    speedgoat_target_model_2021_cal->LaneChanged_Value;

  /* Start for S-Function (slrealtimeUDPSend): '<S3>/UDP Send: Controller Output' */
  {
    try {
      slrealtime::ip::udp::Socket *udpSock = getUDPSocket("0.0.0.0",0U);
      uint8_t *remoteAddress = (uint8_t *)
        speedgoat_target_model_2021_cal->UDPSendControllerOutput_toAddre;
      uint16_t *remotePort = (uint16_t *)
        &speedgoat_target_model_2021_cal->UDPSendControllerOutput_toPort;
      udpSock->setRemoteEndpoint(remoteAddress, remotePort[0]);
      speedgoat_target_model_2021b_DW.UDPSendControllerOutput_IWORK[0] = 24;
      speedgoat_target_model_2021b_DW.UDPSendControllerOutput_IWORK[1] = 0U;
      speedgoat_target_model_2021b_DW.UDPSendControllerOutput_PWORK =
        reinterpret_cast<void*>(udpSock);
    } catch (std::exception& e) {
      std::string tmp = std::string(e.what());
      static std::string eMsg = tmp;
      rtmSetErrorStatus(speedgoat_target_model_2021b_M, eMsg.c_str());
      rtmSetStopRequested(speedgoat_target_model_2021b_M, 1);
      ;
    }
  }

  /* Start for S-Function (slrealtimeUDPSend): '<S3>/UDP Send: Controller Intermediate Signals' */
  {
    try {
      slrealtime::ip::udp::Socket *udpSock = getUDPSocket("0.0.0.0",0U);
      uint8_t *remoteAddress = (uint8_t *)
        speedgoat_target_model_2021_cal->UDPSendControllerIntermediate_k;
      uint16_t *remotePort = (uint16_t *)
        &speedgoat_target_model_2021_cal->UDPSendControllerIntermediateSi;
      udpSock->setRemoteEndpoint(remoteAddress, remotePort[0]);
      speedgoat_target_model_2021b_DW.UDPSendControllerIntermediate_a[0] = 16;
      speedgoat_target_model_2021b_DW.UDPSendControllerIntermediate_a[1] = 0U;
      speedgoat_target_model_2021b_DW.UDPSendControllerIntermediateSi =
        reinterpret_cast<void*>(udpSock);
    } catch (std::exception& e) {
      std::string tmp = std::string(e.what());
      static std::string eMsg = tmp;
      rtmSetErrorStatus(speedgoat_target_model_2021b_M, eMsg.c_str());
      rtmSetStopRequested(speedgoat_target_model_2021b_M, 1);
      ;
    }
  }

  /* Start for RateTransition generated from: '<S43>/optimizer' */
  speedgoat_target_model_2021b_B.min_velocity[0] =
    speedgoat_target_model_2021_cal->min_velocity_InitialCondition;
  speedgoat_target_model_2021b_B.min_velocity[1] =
    speedgoat_target_model_2021_cal->min_velocity_InitialCondition;

  /* Start for RateTransition generated from: '<S43>/optimizer' */
  rtw_slrealtime_mutex_init
    (&speedgoat_target_model_2021b_DW.min_velocity_d0_SEMAPHORE);

  /* Start for RateTransition generated from: '<S4>/DataTypeConversion_L0' */
  speedgoat_target_model_2021b_B.default_spacing =
    speedgoat_target_model_2021_cal->default_spacing_InitialConditio;

  /* Start for RateTransition generated from: '<S4>/DataTypeConversion_L0' */
  rtw_slrealtime_mutex_init
    (&speedgoat_target_model_2021b_DW.default_spacing_d0_SEMAPHORE);

  /* Start for Assertion: '<S51>/Assertion' */
  speedgoat_target_model_2021b_DW.Assertion_sltestFinalResult =
    slTestResult_Untested;
  slTestInitialize(&speedgoat_target_model_2021b_DW.Assertion_sltestBlkInfo,
                   &speedgoat_target_model_2021b_DW.Assertion_sltestCurrentResult,
                   &speedgoat_target_model_2021b_DW.Assertion_sltestFinalResult,
                   &speedgoat_target_model_2021b_DW.Assertion_sltestLastResultTime,
                   1);
  slTestRegAssessment(&speedgoat_target_model_2021b_DW.Assertion_sltestBlkInfo,
                      0, "", 0, 0, 0,
                      "speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/Verify Direction/Assertion",
                      "", 0);

  /* Start for Assertion: '<S51>/Assertion1' */
  speedgoat_target_model_2021b_DW.Assertion1_sltestFinalResult =
    slTestResult_Untested;
  slTestInitialize(&speedgoat_target_model_2021b_DW.Assertion1_sltestBlkInfo,
                   &speedgoat_target_model_2021b_DW.Assertion1_sltestCurrentResult,
                   &speedgoat_target_model_2021b_DW.Assertion1_sltestFinalResult,
                   &speedgoat_target_model_2021b_DW.Assertion1_sltestLastResultTime,
                   1);
  slTestRegAssessment(&speedgoat_target_model_2021b_DW.Assertion1_sltestBlkInfo,
                      0, "", 0, 0, 0,
                      "speedgoat_target_model_2021b/Planning and Control/Longitudinal Controller Stanley/Verify Direction/Assertion1",
                      "", 0);

  {
    int32_T i;

    /* InitializeConditions for UnitDelay: '<S47>/Unit Delay' */
    speedgoat_target_model_2021b_DW.UnitDelay_DSTATE =
      speedgoat_target_model_2021_cal->UnitDelay_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S1>/Unit Delay' */
    speedgoat_target_model_2021b_DW.UnitDelay_DSTATE_o =
      speedgoat_target_model_2021_cal->UnitDelay_InitialCondition_j;

    /* InitializeConditions for RateTransition generated from: '<S11>/Subsystem2' */
    speedgoat_target_model_2021b_DW.acceleration_Buf[0] =
      speedgoat_target_model_2021_cal->acceleration_InitialCondition;
    speedgoat_target_model_2021b_DW.acceleration_WrBufIdx = 0;
    speedgoat_target_model_2021b_DW.acceleration_RdBufIdx = 1;

    /* InitializeConditions for RateTransition generated from: '<S49>/Switch' */
    speedgoat_target_model_2021b_DW.TmpRTBAtSwitchOutport1_Buf[0] =
      speedgoat_target_model_2021_cal->TmpRTBAtSwitchOutport1_InitialC;
    speedgoat_target_model_2021b_DW.TmpRTBAtSwitchOutport1_WrBufIdx = 0;
    speedgoat_target_model_2021b_DW.TmpRTBAtSwitchOutport1_RdBufIdx = 1;

    /* InitializeConditions for RateTransition generated from: '<S49>/Switch1' */
    speedgoat_target_model_2021b_DW.TmpRTBAtSwitch1Outport1_Buf[0] =
      speedgoat_target_model_2021_cal->TmpRTBAtSwitch1Outport1_Initial;
    speedgoat_target_model_2021b_DW.TmpRTBAtSwitch1Outport1_WrBufId = 0;
    speedgoat_target_model_2021b_DW.TmpRTBAtSwitch1Outport1_RdBufId = 1;

    /* InitializeConditions for RateTransition generated from: '<S43>/optimizer' */
    speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport3_Buf[0] =
      speedgoat_target_model_2021_cal->TmpRTBAtoptimizerInport3_Initia;
    speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport3_Buf[1] =
      speedgoat_target_model_2021_cal->TmpRTBAtoptimizerInport3_Initia;
    speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport3_WrBufI = 0;
    speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport3_RdBufI = 1;

    /* InitializeConditions for RateTransition generated from: '<S43>/optimizer' */
    speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport5_Buf[0] =
      speedgoat_target_model_2021_cal->TmpRTBAtoptimizerInport5_Initia;
    speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport5_WrBufI = 0;
    speedgoat_target_model_2021b_DW.TmpRTBAtoptimizerInport5_RdBufI = 1;

    /* InitializeConditions for RateTransition generated from: '<S43>/optimizer' */
    speedgoat_target_model_2021b_DW.min_velocity_Buf0[0] =
      speedgoat_target_model_2021_cal->min_velocity_InitialCondition;
    speedgoat_target_model_2021b_DW.min_velocity_Buf0[1] =
      speedgoat_target_model_2021_cal->min_velocity_InitialCondition;

    /* InitializeConditions for RateTransition generated from: '<S4>/DataTypeConversion_vset' */
    speedgoat_target_model_2021b_DW.TmpRTBAtDataTypeConversion_vset[0] =
      speedgoat_target_model_2021_cal->TmpRTBAtDataTypeConversion_vset;
    speedgoat_target_model_2021b_DW.TmpRTBAtDataTypeConversion_vs_m = 0;
    speedgoat_target_model_2021b_DW.TmpRTBAtDataTypeConversion_vs_n = 1;

    /* InitializeConditions for RateTransition generated from: '<S4>/DataTypeConversion_L0' */
    speedgoat_target_model_2021b_DW.default_spacing_Buf0 =
      speedgoat_target_model_2021_cal->default_spacing_InitialConditio;

    /* InitializeConditions for RateTransition generated from: '<S51>/Equal1' */
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual1Inport1_Buf[0] =
      speedgoat_target_model_2021_cal->TmpRTBAtEqual1Inport1_InitialCo;
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual1Inport1_WrBufIdx = 0;
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual1Inport1_RdBufIdx = 1;

    /* InitializeConditions for RateTransition generated from: '<S51>/Equal2' */
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual2Inport1_Buf[0] =
      speedgoat_target_model_2021_cal->TmpRTBAtEqual2Inport1_InitialCo;
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual2Inport1_WrBufIdx = 0;
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual2Inport1_RdBufIdx = 1;

    /* InitializeConditions for RateTransition generated from: '<S51>/Equal3' */
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual3Inport1_Buf[0] =
      speedgoat_target_model_2021_cal->TmpRTBAtEqual3Inport1_InitialCo;
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual3Inport1_WrBufIdx = 0;
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual3Inport1_RdBufIdx = 1;

    /* InitializeConditions for RateTransition generated from: '<S51>/Equal4' */
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual4Inport2_Buf[0] =
      speedgoat_target_model_2021_cal->TmpRTBAtEqual4Inport2_InitialCo;
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual4Inport2_WrBufIdx = 0;
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual4Inport2_RdBufIdx = 1;

    /* InitializeConditions for RateTransition generated from: '<S51>/Equal5' */
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual5Inport2_Buf[0] =
      speedgoat_target_model_2021_cal->TmpRTBAtEqual5Inport2_InitialCo;
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual5Inport2_WrBufIdx = 0;
    speedgoat_target_model_2021b_DW.TmpRTBAtEqual5Inport2_RdBufIdx = 1;

    /* InitializeConditions for RateTransition generated from: '<S51>/Sign1' */
    speedgoat_target_model_2021b_DW.TmpRTBAtSign1Inport1_Buf[0] =
      speedgoat_target_model_2021_cal->TmpRTBAtSign1Inport1_InitialCon;
    speedgoat_target_model_2021b_DW.TmpRTBAtSign1Inport1_WrBufIdx = 0;
    speedgoat_target_model_2021b_DW.TmpRTBAtSign1Inport1_RdBufIdx = 1;

    /* InitializeConditions for RateTransition generated from: '<S51>/Equal6' */
    speedgoat_target_model_2021b_DW.set_velocity_Buf[0] =
      speedgoat_target_model_2021_cal->set_velocity_InitialCondition;
    speedgoat_target_model_2021b_DW.set_velocity_WrBufIdx = 0;
    speedgoat_target_model_2021b_DW.set_velocity_RdBufIdx = 1;

    /* InitializeConditions for RateTransition generated from: '<S51>/Sign2' */
    speedgoat_target_model_2021b_DW.set_velocity1_Buf[0] =
      speedgoat_target_model_2021_cal->set_velocity1_InitialCondition;
    speedgoat_target_model_2021b_DW.set_velocity1_WrBufIdx = 0;
    speedgoat_target_model_2021b_DW.set_velocity1_RdBufIdx = 1;

    /* InitializeConditions for RateTransition generated from: '<S49>/Multiply' */
    speedgoat_target_model_2021b_DW.TmpRTBAtMultiplyInport1_Buf[0] =
      speedgoat_target_model_2021_cal->TmpRTBAtMultiplyInport1_Initial;
    speedgoat_target_model_2021b_DW.TmpRTBAtMultiplyInport1_WrBufId = 0;
    speedgoat_target_model_2021b_DW.TmpRTBAtMultiplyInport1_RdBufId = 1;

    /* InitializeConditions for RateTransition generated from: '<S50>/Switch Case' */
    speedgoat_target_model_2021b_DW.TmpRTBAtSwitchCaseInport1_Buf[0] =
      speedgoat_target_model_2021_cal->TmpRTBAtSwitchCaseInport1_Initi;
    speedgoat_target_model_2021b_DW.TmpRTBAtSwitchCaseInport1_WrBuf = 0;
    speedgoat_target_model_2021b_DW.TmpRTBAtSwitchCaseInport1_RdBuf = 1;

    /* InitializeConditions for RateTransition generated from: '<S8>/Minus' */
    speedgoat_target_model_2021b_DW.TmpRTBAtMinusInport2_Buf[0] =
      speedgoat_target_model_2021_cal->TmpRTBAtMinusInport2_InitialCon;
    speedgoat_target_model_2021b_DW.TmpRTBAtMinusInport2_WrBufIdx = 0;
    speedgoat_target_model_2021b_DW.TmpRTBAtMinusInport2_RdBufIdx = 1;

    /* InitializeConditions for RateTransition generated from: '<S8>/Minus' */
    speedgoat_target_model_2021b_DW.set_velocity_Buf_n[0] =
      speedgoat_target_model_2021_cal->set_velocity_InitialCondition_h;
    speedgoat_target_model_2021b_DW.set_velocity_WrBufIdx_k = 0;
    speedgoat_target_model_2021b_DW.set_velocity_RdBufIdx_p = 1;

    /* InitializeConditions for Memory: '<S23>/Memory' */
    for (i = 0; i < 34; i++) {
      speedgoat_target_model_2021b_DW.Memory_PreviousInput[i] =
        speedgoat_target_model_2021_cal->Memory_InitialCondition[i];
    }

    /* End of InitializeConditions for Memory: '<S23>/Memory' */

    /* InitializeConditions for UnitDelay: '<S23>/last_mv' */
    speedgoat_target_model_2021b_DW.last_mv_DSTATE =
      speedgoat_target_model_2021_cal->last_mv_InitialCondition;

    /* InitializeConditions for Memory: '<S23>/last_x' */
    speedgoat_target_model_2021b_DW.last_x_PreviousInput[0] =
      speedgoat_target_model_2021_cal->last_x_InitialCondition[0];
    speedgoat_target_model_2021b_DW.last_x_PreviousInput[1] =
      speedgoat_target_model_2021_cal->last_x_InitialCondition[1];
    speedgoat_target_model_2021b_DW.last_x_PreviousInput[2] =
      speedgoat_target_model_2021_cal->last_x_InitialCondition[2];
    speedgoat_target_model_2021b_DW.last_x_PreviousInput[3] =
      speedgoat_target_model_2021_cal->last_x_InitialCondition[3];

    /* SystemInitialize for Chart: '<S1>/Chart' */
    speedgoat_target_model_2021b_DW.sfEvent_o = -1;
    speedgoat_target_model_2021b_DW.is_Navigation = 0U;
    speedgoat_target_model_2021b_DW.is_RailCrossing = 0U;
    speedgoat_target_model_2021b_DW.is_RoundAbout = 0U;
    speedgoat_target_model_2021b_DW.is_StopSign = 0U;
    speedgoat_target_model_2021b_DW.temporalCounter_i1_i = 0U;
    speedgoat_target_model_2021b_DW.is_TrafficLight = 0U;
    speedgoat_target_model_2021b_DW.is_active_c21_speedgoat_target_ = 0U;
    speedgoat_target_model_2021b_DW.is_c21_speedgoat_target_model_2 = 0U;
    speedgoat_target_model_2021b_DW.DestReached = 0.0;
    speedgoat_target_model_2021b_DW.ExitFlag = 0.0;
    speedgoat_target_model_2021b_DW.CurrentVelCmd = 0.0;
    speedgoat_target_model_2021b_DW.Navi = 0.0;
    speedgoat_target_model_2021b_B.VelCmd = 0.0;
    speedgoat_target_model_2021b_B.ACC_i = 0.0;
    speedgoat_target_model_2021b_B.Tgap = 0.0;
    speedgoat_target_model_2021b_B.Trans = 0.0;
    speedgoat_target_model_2021b_B.StopSim = 0.0;

    /* SystemInitialize for Chart: '<S1>/Collision Avoidance' */
    speedgoat_target_model_2021b_DW.sfEvent = -1;
    speedgoat_target_model_2021b_DW.temporalCounter_i1 = 0U;
    speedgoat_target_model_2021b_DW.is_active_c9_speedgoat_target_m = 0U;
    speedgoat_target_model_2021b_DW.is_c9_speedgoat_target_model_20 = 0U;
    speedgoat_target_model_2021b_DW.Stuck = 0.0;
    speedgoat_target_model_2021b_B.ACC = 0.0;
    speedgoat_target_model_2021b_B.LaneChangeCmd = 0.0;

    /* SystemInitialize for IfAction SubSystem: '<S50>/Forward' */
    /* InitializeConditions for DiscreteIntegrator: '<S88>/Integrator' */
    speedgoat_target_model_2021b_DW.Integrator_DSTATE_k =
      speedgoat_target_model_2021_cal->PIForward_InitialConditionForIn;
    speedgoat_target_model_2021b_DW.Integrator_PrevResetState_j = 0;

    /* End of SystemInitialize for SubSystem: '<S50>/Forward' */

    /* SystemInitialize for IfAction SubSystem: '<S50>/Reverse' */
    /* InitializeConditions for DiscreteIntegrator: '<S139>/Integrator' */
    speedgoat_target_model_2021b_DW.Integrator_DSTATE =
      speedgoat_target_model_2021_cal->PIReverse_InitialConditionForIn;
    speedgoat_target_model_2021b_DW.Integrator_PrevResetState = 0;

    /* End of SystemInitialize for SubSystem: '<S50>/Reverse' */

    /* SystemInitialize for Merge: '<S50>/Merge' */
    speedgoat_target_model_2021b_B.Merge =
      speedgoat_target_model_2021_cal->Merge_InitialOutput;

    /* SystemInitialize for IfAction SubSystem: '<S11>/Subsystem2' */
    /* InitializeConditions for DiscreteIntegrator: '<S164>/Discrete-Time Integrator' */
    speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_DSTATE_f =
      speedgoat_target_model_2021_cal->DiscreteTimeIntegrator_IC;
    speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevRe_b = 2;

    /* InitializeConditions for DiscreteIntegrator: '<S163>/Discrete-Time Integrator' */
    speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_DSTATE_e =
      speedgoat_target_model_2021_cal->DiscreteTimeIntegrator_IC_a;
    speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevR_bg = 2;

    /* End of SystemInitialize for SubSystem: '<S11>/Subsystem2' */

    /* SystemInitialize for IfAction SubSystem: '<S11>/Subsystem1' */
    /* InitializeConditions for DiscreteIntegrator: '<S160>/Discrete-Time Integrator' */
    speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_DSTATE =
      speedgoat_target_model_2021_cal->DiscreteTimeIntegrator_IC_i;
    speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevRese = 2;

    /* InitializeConditions for DiscreteIntegrator: '<S159>/Discrete-Time Integrator' */
    speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_DSTATE_a =
      speedgoat_target_model_2021_cal->DiscreteTimeIntegrator_IC_d;
    speedgoat_target_model_2021b_DW.DiscreteTimeIntegrator_PrevRe_a = 2;

    /* End of SystemInitialize for SubSystem: '<S11>/Subsystem1' */

    /* SystemInitialize for Merge: '<S11>/Merge1' */
    speedgoat_target_model_2021b_B.Merge1 =
      speedgoat_target_model_2021_cal->Merge1_InitialOutput;

    /* SystemInitialize for Merge: '<S11>/Merge2' */
    speedgoat_target_model_2021b_B.Merge2 =
      speedgoat_target_model_2021_cal->Merge2_InitialOutput;

    /* SystemInitialize for Enabled SubSystem: '<S2>/Byte Unpack: First set' */
    /* SystemInitialize for S-Function (slrealtimebytepacking): '<S166>/Byte Unpacking ' incorporates:
     *  Outport: '<S166>/Obj Ahead'
     */
    speedgoat_target_model_2021b_B.ObjAhead =
      speedgoat_target_model_2021_cal->ObjAhead_Y0;

    /* SystemInitialize for S-Function (slrealtimebytepacking): '<S166>/Byte Unpacking ' incorporates:
     *  Outport: '<S166>/Longitudinal velocity'
     */
    speedgoat_target_model_2021b_B.Longitudinalvelocity =
      speedgoat_target_model_2021_cal->Longitudinalvelocity_Y0;

    /* SystemInitialize for S-Function (slrealtimebytepacking): '<S166>/Byte Unpacking ' incorporates:
     *  Outport: '<S166>/SceneID'
     */
    speedgoat_target_model_2021b_B.SceneID =
      speedgoat_target_model_2021_cal->SceneID_Y0;

    /* SystemInitialize for S-Function (slrealtimebytepacking): '<S166>/Byte Unpacking ' incorporates:
     *  Outport: '<S166>/Relative velocity'
     */
    speedgoat_target_model_2021b_B.Relativevelocity =
      speedgoat_target_model_2021_cal->Relativevelocity_Y0;

    /* End of SystemInitialize for SubSystem: '<S2>/Byte Unpack: First set' */

    /* SystemInitialize for Iterator SubSystem: '<S2>/UDP Receive: First set' */
    /* Start for S-Function (slrealtimeUDPReceive): '<S171>/UDP Receive2' */
    {
      try {
        uint8_t *tempAddr = nullptr;
        uint8_t *tempInterface = nullptr;
        slrealtime::ip::udp::Socket *udpSock = getUDPSocket("0.0.0.0",8001U);
        if (tempAddr)
          delete tempAddr;
        if (tempInterface)
          delete tempInterface;
        speedgoat_target_model_2021b_DW.UDPReceive2_IWORK_h = 32;
        speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_e[0] =
          reinterpret_cast<void*>(udpSock);
        char *buffer = (char *)calloc(65504,sizeof(char));
        speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_e[1] = (void*)buffer;
      } catch (std::exception& e) {
        std::string tmp = std::string(e.what());
        static std::string eMsg = tmp;
        rtmSetErrorStatus(speedgoat_target_model_2021b_M, eMsg.c_str());
        rtmSetStopRequested(speedgoat_target_model_2021b_M, 1);
        ;
      }
    }

    /* SystemInitialize for Sum: '<S171>/Subtract' incorporates:
     *  Outport: '<S171>/NumMsg'
     */
    speedgoat_target_model_2021b_B.Subtract_b =
      speedgoat_target_model_2021_cal->NumMsg_Y0_b;
    for (i = 0; i < 32; i++) {
      /* SystemInitialize for S-Function (slrealtimeUDPReceive): '<S171>/UDP Receive2' incorporates:
       *  Outport: '<S171>/Data'
       */
      speedgoat_target_model_2021b_B.UDPReceive2_o1_p[i] =
        speedgoat_target_model_2021_cal->Data_Y0_j;
    }

    /* End of SystemInitialize for SubSystem: '<S2>/UDP Receive: First set' */

    /* SystemInitialize for Enabled SubSystem: '<S2>/Byte Unpack: Second set' */
    /* SystemInitialize for S-Function (slrealtimebytepacking): '<S168>/Byte Unpacking ' incorporates:
     *  Outport: '<S168>/ObjectDist'
     */
    speedgoat_target_model_2021b_B.ObjectDist =
      speedgoat_target_model_2021_cal->ObjectDist_Y0;

    /* SystemInitialize for S-Function (slrealtimebytepacking): '<S168>/Byte Unpacking ' incorporates:
     *  Outport: '<S168>/D2S'
     */
    speedgoat_target_model_2021b_B.D2S = speedgoat_target_model_2021_cal->D2S_Y0;

    /* SystemInitialize for S-Function (slrealtimebytepacking): '<S168>/Byte Unpacking ' incorporates:
     *  Outport: '<S168>/Relative distance'
     */
    speedgoat_target_model_2021b_B.Relativedistance =
      speedgoat_target_model_2021_cal->Relativedistance_Y0;
    for (i = 0; i < 90; i++) {
      /* SystemInitialize for S-Function (slrealtimebytepacking): '<S168>/Byte Unpacking ' incorporates:
       *  Outport: '<S168>/TL state'
       */
      speedgoat_target_model_2021b_B.TLstate[i] =
        speedgoat_target_model_2021_cal->TLstate_Y0;
    }

    /* End of SystemInitialize for SubSystem: '<S2>/Byte Unpack: Second set' */

    /* SystemInitialize for Iterator SubSystem: '<S2>/UDP Receive: Second set' */
    /* Start for S-Function (slrealtimeUDPReceive): '<S173>/UDP Receive2' */
    {
      try {
        uint8_t *tempAddr = nullptr;
        uint8_t *tempInterface = nullptr;
        slrealtime::ip::udp::Socket *udpSock = getUDPSocket("0.0.0.0",8002U);
        if (tempAddr)
          delete tempAddr;
        if (tempInterface)
          delete tempInterface;
        speedgoat_target_model_2021b_DW.UDPReceive2_IWORK_l = 744;
        speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_b[0] =
          reinterpret_cast<void*>(udpSock);
        char *buffer = (char *)calloc(65504,sizeof(char));
        speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_b[1] = (void*)buffer;
      } catch (std::exception& e) {
        std::string tmp = std::string(e.what());
        static std::string eMsg = tmp;
        rtmSetErrorStatus(speedgoat_target_model_2021b_M, eMsg.c_str());
        rtmSetStopRequested(speedgoat_target_model_2021b_M, 1);
        ;
      }
    }

    /* SystemInitialize for Sum: '<S173>/Subtract' incorporates:
     *  Outport: '<S173>/NumMsg'
     */
    speedgoat_target_model_2021b_B.Subtract_p =
      speedgoat_target_model_2021_cal->NumMsg_Y0_f;
    for (i = 0; i < 744; i++) {
      /* SystemInitialize for S-Function (slrealtimeUDPReceive): '<S173>/UDP Receive2' incorporates:
       *  Outport: '<S173>/Data'
       */
      speedgoat_target_model_2021b_B.UDPReceive2_o1_o[i] =
        speedgoat_target_model_2021_cal->Data_Y0_ja;
    }

    /* End of SystemInitialize for SubSystem: '<S2>/UDP Receive: Second set' */

    /* SystemInitialize for Enabled SubSystem: '<S2>/Byte Unpack: Third set' */
    /* SystemInitialize for S-Function (slrealtimebytepacking): '<S169>/Byte Unpacking ' incorporates:
     *  Outport: '<S169>/Direction'
     */
    speedgoat_target_model_2021b_B.Direction =
      speedgoat_target_model_2021_cal->Direction_Y0;

    /* SystemInitialize for S-Function (slrealtimebytepacking): '<S169>/Byte Unpacking ' incorporates:
     *  Outport: '<S169>/Curvature'
     */
    speedgoat_target_model_2021b_B.Curvature =
      speedgoat_target_model_2021_cal->Curvature_Y0;

    /* SystemInitialize for S-Function (slrealtimebytepacking): '<S169>/Byte Unpacking ' incorporates:
     *  Outport: '<S169>/yaw rate'
     */
    speedgoat_target_model_2021b_B.yawrate =
      speedgoat_target_model_2021_cal->yawrate_Y0;

    /* End of SystemInitialize for SubSystem: '<S2>/Byte Unpack: Third set' */

    /* SystemInitialize for Iterator SubSystem: '<S2>/UDP Receive: Third set' */
    /* Start for S-Function (slrealtimeUDPReceive): '<S174>/UDP Receive2' */
    {
      try {
        uint8_t *tempAddr = nullptr;
        uint8_t *tempInterface = nullptr;
        slrealtime::ip::udp::Socket *udpSock = getUDPSocket("0.0.0.0",8003U);
        if (tempAddr)
          delete tempAddr;
        if (tempInterface)
          delete tempInterface;
        speedgoat_target_model_2021b_DW.UDPReceive2_IWORK = 24;
        speedgoat_target_model_2021b_DW.UDPReceive2_PWORK[0] = reinterpret_cast<
          void*>(udpSock);
        char *buffer = (char *)calloc(65504,sizeof(char));
        speedgoat_target_model_2021b_DW.UDPReceive2_PWORK[1] = (void*)buffer;
      } catch (std::exception& e) {
        std::string tmp = std::string(e.what());
        static std::string eMsg = tmp;
        rtmSetErrorStatus(speedgoat_target_model_2021b_M, eMsg.c_str());
        rtmSetStopRequested(speedgoat_target_model_2021b_M, 1);
        ;
      }
    }

    /* SystemInitialize for Sum: '<S174>/Subtract' incorporates:
     *  Outport: '<S174>/NumMsg'
     */
    speedgoat_target_model_2021b_B.Subtract =
      speedgoat_target_model_2021_cal->NumMsg_Y0_c;
    for (i = 0; i < 24; i++) {
      /* SystemInitialize for S-Function (slrealtimeUDPReceive): '<S174>/UDP Receive2' incorporates:
       *  Outport: '<S174>/Data'
       */
      speedgoat_target_model_2021b_B.UDPReceive2_o1[i] =
        speedgoat_target_model_2021_cal->Data_Y0_l;
    }

    /* End of SystemInitialize for SubSystem: '<S2>/UDP Receive: Third set' */

    /* SystemInitialize for Enabled SubSystem: '<S2>/Byte Unpack: Fourth set' */
    /* SystemInitialize for S-Function (slrealtimebytepacking): '<S167>/Byte Unpacking ' incorporates:
     *  Outport: '<S167>/RefPoses'
     */
    speedgoat_target_model_2021b_B.RefPoses[0] =
      speedgoat_target_model_2021_cal->RefPoses_Y0;
    speedgoat_target_model_2021b_B.RefPoses[1] =
      speedgoat_target_model_2021_cal->RefPoses_Y0;
    speedgoat_target_model_2021b_B.RefPoses[2] =
      speedgoat_target_model_2021_cal->RefPoses_Y0;

    /* End of SystemInitialize for SubSystem: '<S2>/Byte Unpack: Fourth set' */

    /* SystemInitialize for Iterator SubSystem: '<S2>/UDP Receive: Fourth set' */
    /* Start for S-Function (slrealtimeUDPReceive): '<S172>/UDP Receive2' */
    {
      try {
        uint8_t *tempAddr = nullptr;
        uint8_t *tempInterface = nullptr;
        slrealtime::ip::udp::Socket *udpSock = getUDPSocket("0.0.0.0",8004U);
        if (tempAddr)
          delete tempAddr;
        if (tempInterface)
          delete tempInterface;
        speedgoat_target_model_2021b_DW.UDPReceive2_IWORK_lp = 24;
        speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_f[0] =
          reinterpret_cast<void*>(udpSock);
        char *buffer = (char *)calloc(65504,sizeof(char));
        speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_f[1] = (void*)buffer;
      } catch (std::exception& e) {
        std::string tmp = std::string(e.what());
        static std::string eMsg = tmp;
        rtmSetErrorStatus(speedgoat_target_model_2021b_M, eMsg.c_str());
        rtmSetStopRequested(speedgoat_target_model_2021b_M, 1);
        ;
      }
    }

    /* SystemInitialize for Sum: '<S172>/Subtract' incorporates:
     *  Outport: '<S172>/NumMsg'
     */
    speedgoat_target_model_2021b_B.Subtract_d =
      speedgoat_target_model_2021_cal->NumMsg_Y0_bf;
    for (i = 0; i < 24; i++) {
      /* SystemInitialize for S-Function (slrealtimeUDPReceive): '<S172>/UDP Receive2' incorporates:
       *  Outport: '<S172>/Data'
       */
      speedgoat_target_model_2021b_B.UDPReceive2_o1_d[i] =
        speedgoat_target_model_2021_cal->Data_Y0_k;
    }

    /* End of SystemInitialize for SubSystem: '<S2>/UDP Receive: Fourth set' */

    /* SystemInitialize for Enabled SubSystem: '<S2>/Byte Unpack: Fifth set' */
    /* SystemInitialize for S-Function (slrealtimebytepacking): '<S165>/Byte Unpacking 1' incorporates:
     *  Outport: '<S165>/CurrPoses'
     */
    speedgoat_target_model_2021b_B.CurrPoses[0] =
      speedgoat_target_model_2021_cal->CurrPoses_Y0;
    speedgoat_target_model_2021b_B.CurrPoses[1] =
      speedgoat_target_model_2021_cal->CurrPoses_Y0;
    speedgoat_target_model_2021b_B.CurrPoses[2] =
      speedgoat_target_model_2021_cal->CurrPoses_Y0;

    /* End of SystemInitialize for SubSystem: '<S2>/Byte Unpack: Fifth set' */

    /* SystemInitialize for Iterator SubSystem: '<S2>/UDP Receive: Fifth set' */
    /* Start for S-Function (slrealtimeUDPReceive): '<S170>/UDP Receive2' */
    {
      try {
        uint8_t *tempAddr = nullptr;
        uint8_t *tempInterface = nullptr;
        slrealtime::ip::udp::Socket *udpSock = getUDPSocket("0.0.0.0",8005U);
        if (tempAddr)
          delete tempAddr;
        if (tempInterface)
          delete tempInterface;
        speedgoat_target_model_2021b_DW.UDPReceive2_IWORK_a = 24;
        speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_e0[0] =
          reinterpret_cast<void*>(udpSock);
        char *buffer = (char *)calloc(65504,sizeof(char));
        speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_e0[1] = (void*)buffer;
      } catch (std::exception& e) {
        std::string tmp = std::string(e.what());
        static std::string eMsg = tmp;
        rtmSetErrorStatus(speedgoat_target_model_2021b_M, eMsg.c_str());
        rtmSetStopRequested(speedgoat_target_model_2021b_M, 1);
        ;
      }
    }

    /* SystemInitialize for Sum: '<S170>/Subtract' incorporates:
     *  Outport: '<S170>/NumMsg'
     */
    speedgoat_target_model_2021b_B.Subtract_bo =
      speedgoat_target_model_2021_cal->NumMsg_Y0;
    for (i = 0; i < 24; i++) {
      /* SystemInitialize for S-Function (slrealtimeUDPReceive): '<S170>/UDP Receive2' incorporates:
       *  Outport: '<S170>/Data'
       */
      speedgoat_target_model_2021b_B.UDPReceive2_o1_m[i] =
        speedgoat_target_model_2021_cal->Data_Y0;
    }

    /* End of SystemInitialize for SubSystem: '<S2>/UDP Receive: Fifth set' */
  }
}

/* Model terminate function */
void speedgoat_target_model_2021b_terminate(void)
{
  /* Terminate for Iterator SubSystem: '<S2>/UDP Receive: Fourth set' */

  /* Terminate for S-Function (slrealtimeUDPReceive): '<S172>/UDP Receive2' */
  {
    slrealtime::ip::udp::Socket *udpSock = getUDPSocket("0.0.0.0",8004U);
    if (udpSock)
      delete udpSock;
    char *buffer = (char *)speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_f[1];
    if (buffer)
      free(buffer);
  }

  /* End of Terminate for SubSystem: '<S2>/UDP Receive: Fourth set' */

  /* Terminate for Iterator SubSystem: '<S2>/UDP Receive: Fifth set' */

  /* Terminate for S-Function (slrealtimeUDPReceive): '<S170>/UDP Receive2' */
  {
    slrealtime::ip::udp::Socket *udpSock = getUDPSocket("0.0.0.0",8005U);
    if (udpSock)
      delete udpSock;
    char *buffer = (char *)speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_e0
      [1];
    if (buffer)
      free(buffer);
  }

  /* End of Terminate for SubSystem: '<S2>/UDP Receive: Fifth set' */

  /* Terminate for Iterator SubSystem: '<S2>/UDP Receive: First set' */

  /* Terminate for S-Function (slrealtimeUDPReceive): '<S171>/UDP Receive2' */
  {
    slrealtime::ip::udp::Socket *udpSock = getUDPSocket("0.0.0.0",8001U);
    if (udpSock)
      delete udpSock;
    char *buffer = (char *)speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_e[1];
    if (buffer)
      free(buffer);
  }

  /* End of Terminate for SubSystem: '<S2>/UDP Receive: First set' */

  /* Terminate for Iterator SubSystem: '<S2>/UDP Receive: Third set' */

  /* Terminate for S-Function (slrealtimeUDPReceive): '<S174>/UDP Receive2' */
  {
    slrealtime::ip::udp::Socket *udpSock = getUDPSocket("0.0.0.0",8003U);
    if (udpSock)
      delete udpSock;
    char *buffer = (char *)speedgoat_target_model_2021b_DW.UDPReceive2_PWORK[1];
    if (buffer)
      free(buffer);
  }

  /* End of Terminate for SubSystem: '<S2>/UDP Receive: Third set' */

  /* Terminate for Iterator SubSystem: '<S2>/UDP Receive: Second set' */

  /* Terminate for S-Function (slrealtimeUDPReceive): '<S173>/UDP Receive2' */
  {
    slrealtime::ip::udp::Socket *udpSock = getUDPSocket("0.0.0.0",8002U);
    if (udpSock)
      delete udpSock;
    char *buffer = (char *)speedgoat_target_model_2021b_DW.UDPReceive2_PWORK_b[1];
    if (buffer)
      free(buffer);
  }

  /* End of Terminate for SubSystem: '<S2>/UDP Receive: Second set' */

  /* Terminate for S-Function (slrealtimeUDPSend): '<S3>/UDP Send: Controller Output' */
  {
    slrealtime::ip::udp::Socket *udpSock = getUDPSocket("0.0.0.0",0U);
    if (udpSock)
      delete udpSock;
  }

  /* Terminate for S-Function (slrealtimeUDPSend): '<S3>/UDP Send: Controller Intermediate Signals' */
  {
    slrealtime::ip::udp::Socket *udpSock = getUDPSocket("0.0.0.0",0U);
    if (udpSock)
      delete udpSock;
  }

  /* Terminate for RateTransition generated from: '<S43>/optimizer' */
  rtw_slrealtime_mutex_destroy
    (speedgoat_target_model_2021b_DW.min_velocity_d0_SEMAPHORE);

  /* Terminate for RateTransition generated from: '<S4>/DataTypeConversion_L0' */
  rtw_slrealtime_mutex_destroy
    (speedgoat_target_model_2021b_DW.default_spacing_d0_SEMAPHORE);

  /* Terminate for Assertion: '<S51>/Assertion' */
  slTestTerminate(&speedgoat_target_model_2021b_DW.Assertion_sltestBlkInfo);

  /* Terminate for Assertion: '<S51>/Assertion1' */
  slTestTerminate(&speedgoat_target_model_2021b_DW.Assertion1_sltestBlkInfo);
}
