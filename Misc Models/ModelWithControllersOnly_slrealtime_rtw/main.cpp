/* Main generated for Simulink Real-Time model ModelWithControllersOnly */
#include <ModelInfo.hpp>
#include <utilities.hpp>
#include "ModelWithControllersOnly.h"
#include "rte_ModelWithControllersOnly_parameters.h"

/* Task descriptors */
slrealtime::TaskInfo task_1( 0u, std::bind(ModelWithControllersOnly_step0), slrealtime::TaskInfo::PERIODIC, 0.01, 0, 40);
slrealtime::TaskInfo task_2( 1u, std::bind(ModelWithControllersOnly_step2), slrealtime::TaskInfo::PERIODIC, 0.05, 0, 39);

/* Executable base address for XCP */
#ifdef __linux__
extern char __executable_start;
static uintptr_t const base_address = reinterpret_cast<uintptr_t>(&__executable_start);
#else
/* Set 0 as placeholder, to be parsed later from /proc filesystem */
static uintptr_t const base_address = 0;
#endif

/* Model descriptor */
slrealtime::ModelInfo ModelWithControllersOnly_Info =
{
    "ModelWithControllersOnly",
    ModelWithControllersOnly_initialize,
    ModelWithControllersOnly_terminate,
    []()->char const*& { return ModelWithControllersOnly_M->errorStatus; },
    []()->unsigned char& { return ModelWithControllersOnly_M->Timing.stopRequestedFlag; },
    { task_1, task_2 },
    slrealtime::getSegmentVector()
};

int main(int argc, char *argv[]) {
    slrealtime::BaseAddress::set(base_address);
    return slrealtime::runModel(argc, argv, ModelWithControllersOnly_Info);
}
