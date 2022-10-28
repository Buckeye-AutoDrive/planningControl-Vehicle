#include "rte_speedgoat_target_model_2021b_parameters.h"
#include "speedgoat_target_model_2021b.h"
#include "speedgoat_target_model_2021b_cal.h"

extern speedgoat_target_model_cal_type speedgoat_target_model_cal_impl;
namespace slrealtime
{
  /* Description of SEGMENTS */
  SegmentVector segmentInfo {
    { (void*)&speedgoat_target_model_cal_impl, (void**)
      &speedgoat_target_model_2021_cal, sizeof(speedgoat_target_model_cal_type),
      2 }
  };

  SegmentVector &getSegmentVector(void)
  {
    return segmentInfo;
  }
}                                      // slrealtime
