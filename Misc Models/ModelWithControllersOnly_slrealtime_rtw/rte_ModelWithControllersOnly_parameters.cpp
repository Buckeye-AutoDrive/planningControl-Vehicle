#include "rte_ModelWithControllersOnly_parameters.h"
#include "ModelWithControllersOnly.h"
#include "ModelWithControllersOnly_cal.h"

RTE_Param_Service_T RTE_Param_Service = {
  { 1.0, -1.803005404306504, 0.81115986469468793 },

  { 1.0, -1.8952164009111621, 0.901290812452544 },

  { 1.0, -1.6954106486625078, 0.72310058969318891 },
  1.3425,
  1.3425,
  0.0,
  540.0,
  22534.0,
  1853.0,
  -65534.0,
  -540.0,
  -22534.0,

  { 0.0020386150970459851, 0.00407723019409197, 0.0020386150970459851 },

  { 0.0015186028853454824, 0.0030372057706909649, 0.0015186028853454824 },

  { 0.0069224852576702861, 0.013844970515340572, 0.0069224852576702861 },
  0.21590043180086363,
  54.5,
  -0.058221,
  10.0,
  0.02,
  -4.1603370851845,
  1.3672916521692537,
  7.8974120497198439E-7
};

RTE_Param_Service_T *RTE_Param_Service_ptr = &RTE_Param_Service;
real_T* get_DEND_B(void)
{
  return RTE_Param_Service_ptr->DEND_B;
}

real_T* get_DEND_S(void)
{
  return RTE_Param_Service_ptr->DEND_S;
}

real_T* get_DEND_T(void)
{
  return RTE_Param_Service_ptr->DEND_T;
}

real_T* get_LF(void)
{
  return &RTE_Param_Service_ptr->LF;
}

real_T* get_LR(void)
{
  return &RTE_Param_Service_ptr->LR;
}

real_T* get_MAX_BRAKE(void)
{
  return &RTE_Param_Service_ptr->MAX_BRAKE;
}

real_T* get_MAX_STEER(void)
{
  return &RTE_Param_Service_ptr->MAX_STEER;
}

real_T* get_MAX_TORQUE(void)
{
  return &RTE_Param_Service_ptr->MAX_TORQUE;
}

real_T* get_MCAR(void)
{
  return &RTE_Param_Service_ptr->MCAR;
}

real_T* get_MIN_BRAKE(void)
{
  return &RTE_Param_Service_ptr->MIN_BRAKE;
}

real_T* get_MIN_STEER(void)
{
  return &RTE_Param_Service_ptr->MIN_STEER;
}

real_T* get_MIN_TORQUE(void)
{
  return &RTE_Param_Service_ptr->MIN_TORQUE;
}

real_T* get_NUMD_B(void)
{
  return RTE_Param_Service_ptr->NUMD_B;
}

real_T* get_NUMD_S(void)
{
  return RTE_Param_Service_ptr->NUMD_S;
}

real_T* get_NUMD_T(void)
{
  return RTE_Param_Service_ptr->NUMD_T;
}

real_T* get_R(void)
{
  return &RTE_Param_Service_ptr->R;
}

real_T* get_RRdamp(void)
{
  return &RTE_Param_Service_ptr->RRdamp;
}

real_T* get_STEER_RATIO(void)
{
  return &RTE_Param_Service_ptr->STEER_RATIO;
}

real_T* get_Tdead(void)
{
  return &RTE_Param_Service_ptr->Tdead;
}

real_T* get_Vdead(void)
{
  return &RTE_Param_Service_ptr->Vdead;
}

real_T* get_X_o(void)
{
  return &RTE_Param_Service_ptr->X_o;
}

real_T* get_Y_o(void)
{
  return &RTE_Param_Service_ptr->Y_o;
}

real_T* get_psi_o(void)
{
  return &RTE_Param_Service_ptr->psi_o;
}

extern ModelWithControllersOn_cal_type ModelWithControllersOn_cal_impl;
extern RTE_Param_Service_T RTE_Param_Service;
namespace slrealtime
{
  /* Description of SEGMENTS */
  SegmentVector segmentInfo {
    { (void*)&RTE_Param_Service, (void**)&RTE_Param_Service_ptr, sizeof
      (RTE_Param_Service_T), 2 },

    { (void*)&ModelWithControllersOn_cal_impl, (void**)
      &ModelWithControllersOnly_cal, sizeof(ModelWithControllersOn_cal_type), 2
    }
  };

  SegmentVector &getSegmentVector(void)
  {
    return segmentInfo;
  }
}                                      // slrealtime
