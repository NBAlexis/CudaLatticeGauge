//=============================================================================
// FILENAME : CudaHelper.h
// 
// DESCRIPTION:
// This is the file for some common CUDA usage
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CUDAHELPER_H_
#define _CUDAHELPER_H_

__BEGIN_NAMESPACE

extern "C" bool runCudaTest(const int argc, const char **argv,
    char *data, int2 *data_int2, unsigned int len);

extern "C" int GetCudaGPUCount();

class CLGAPI CCudaHelper
{
public:
    static void PrintGPUCount();
    static void DeviceQuery();
};

__END_NAMESPACE


#endif //#ifndef _CUDAHELPER_H_

//=============================================================================
// END OF FILE
//=============================================================================