//=============================================================================
// FILENAME : TestImprove5BadLaunch.cu
//
// DESCRIPTION:
// Improve-1 I5 negative compile check (multi-GPU-improve1.md 3.3 gate): the
// raw <<<>>> backend must KEEP compile-time type checking of the kernel
// arguments. This file intentionally launches a kernel with a mismatched
// argument list; it must FAIL to compile under _CLG_LAUNCH_KERNEL=0.
//
// NOT part of any normal build target -- compiled only by the opt-in
// CheckRawLaunchTypeFail target (Tools/check_raw_launch_typefail.sh), which
// EXPECTS the failure. Under _CLG_LAUNCH_KERNEL=1 (launchKernel, type
// erasure) the bad call compiles; the script reports that as "not a raw
// build".
//
// REVISION:
//  [08/09/2026 Improve-1 I5 nbale]
//=============================================================================
#include "CLGLib.h"

__BEGIN_NAMESPACE

__global__ void _kernelImprove5TypeCheck(INT* pValue)
{
    if (NULL != pValue)
    {
        pValue[0] = 0;
    }
}

void appImprove5BadLaunch(INT* pValue)
{
    //Type error on purpose: an extra const char* the kernel does not take.
    _LAUNCH_KERNEL(_kernelImprove5TypeCheck, 1, 1, pValue, "type-mismatch");
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
