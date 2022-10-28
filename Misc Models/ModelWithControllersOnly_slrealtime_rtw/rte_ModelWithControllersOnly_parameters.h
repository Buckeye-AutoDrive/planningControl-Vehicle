#ifndef _RTE_MODELWITHCONTROLLERSONLY_PARAMETERS_H
#define _RTE_MODELWITHCONTROLLERSONLY_PARAMETERS_H
#include "rtwtypes.h"
#include "SegmentInfo.hpp"
#include "multiword_types.h"
#include "zero_crossing_types.h"
#include "ModelWithControllersOnly_types.h"

struct RTE_Param_Service_T {
  real_T DEND_B[3];
  real_T DEND_S[3];
  real_T DEND_T[3];
  real_T LF;
  real_T LR;
  real_T MAX_BRAKE;
  real_T MAX_STEER;
  real_T MAX_TORQUE;
  real_T MCAR;
  real_T MIN_BRAKE;
  real_T MIN_STEER;
  real_T MIN_TORQUE;
  real_T NUMD_B[3];
  real_T NUMD_S[3];
  real_T NUMD_T[3];
  real_T R;
  real_T RRdamp;
  real_T STEER_RATIO;
  real_T Tdead;
  real_T Vdead;
  real_T X_o;
  real_T Y_o;
  real_T psi_o;
};

extern RTE_Param_Service_T RTE_Param_Service;
extern RTE_Param_Service_T *RTE_Param_Service_ptr;
real_T* get_DEND_B(void);
real_T* get_DEND_S(void);
real_T* get_DEND_T(void);
real_T* get_LF(void);
real_T* get_LR(void);
real_T* get_MAX_BRAKE(void);
real_T* get_MAX_STEER(void);
real_T* get_MAX_TORQUE(void);
real_T* get_MCAR(void);
real_T* get_MIN_BRAKE(void);
real_T* get_MIN_STEER(void);
real_T* get_MIN_TORQUE(void);
real_T* get_NUMD_B(void);
real_T* get_NUMD_S(void);
real_T* get_NUMD_T(void);
real_T* get_R(void);
real_T* get_RRdamp(void);
real_T* get_STEER_RATIO(void);
real_T* get_Tdead(void);
real_T* get_Vdead(void);
real_T* get_X_o(void);
real_T* get_Y_o(void);
real_T* get_psi_o(void);
namespace slrealtime
{
  SegmentVector &getSegmentVector(void);
}                                      // slrealtime

#endif
