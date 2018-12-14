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


enum { kContentLength = 1024,};

extern __constant__ UINT _constIntegers[kContentLength];
extern __constant__ FLOAT _constFloats[kContentLength];

extern __constant__ class CRandom* __r;
extern __constant__ class CRandomSchrage* __rs;


enum EConstIntId
{
    ECI_Dim,
    ECI_Dir,
    ECI_Lx,
    ECI_Ly,
    ECI_Lz,
    ECI_Lt,
    ECI_Volumn,
    ECI_MultX,
    ECI_MultY,
    ECI_MultZ,
    ECI_DecompX, //number of blocks
    ECI_DecompY,
    ECI_DecompZ,
    ECI_DecompLx, //threads per block
    ECI_DecompLy,
    ECI_DecompLz,
    ECI_RandomSeed,
    ECI_ExponentPrecision,
    ECI_IntegratorStepCount,
    ECI_UsingSchrageRandom,

    ECI_ForceDWORD = 0x7fffffff,
};

enum EConstFloatId
{
    ECF_PlaqutteBeta,
};

class CLGAPI CCudaHelper
{
public:
    CCudaHelper()
    {
        memset(m_ConstIntegers, 0, sizeof(UINT) * kContentLength);
        memset(m_ConstFloats, 0, sizeof(FLOAT) * kContentLength);
    }

    static void DeviceQuery();
    void CopyConstants() const;
    void CopyRandomPointer(const class CRandom* r, const class CRandomSchrage* rs);
    /**ret[0] = max thread count, ret[1,2,3] = max thread for x,y,z per block*/
    static TArray<UINT> GetMaxThreadCountAndThreadPerblock();


    UINT m_ConstIntegers[kContentLength];
    FLOAT m_ConstFloats[kContentLength];
};

__END_NAMESPACE


#endif //#ifndef _CUDAHELPER_H_

//=============================================================================
// END OF FILE
//=============================================================================