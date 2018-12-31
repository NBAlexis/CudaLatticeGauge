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


enum { kSharedLength = 1024, };

enum { kContentLength = 1024,};

extern __constant__ UINT _constIntegers[kContentLength];
extern __constant__ Real _constFloats[kContentLength];

extern __constant__ class CRandom* __r;
extern __constant__ class CRandomSchrage* __rs;
extern __constant__ class CIndex* __idx;

extern __constant__ class gammaMatrixSet* __diracGamma;
extern __constant__ class gammaMatrixSet* __chiralGamma;

extern __constant__ class deviceSU3* __SU3Generators[8];

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
    ECI_UsingSchrageRandom,
    ECI_ActionListLength,

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
        memset(m_ConstFloats, 0, sizeof(Real) * kContentLength);
    }

    static void DeviceQuery();
    void CopyConstants() const;
    void CopyRandomPointer(const class CRandom* r, const class CRandomSchrage* rs) const;
    void SetDeviceIndex(class CIndex** ppIdx) const;
    //we never need gamma matrix on host, so this is purely hiden in device
    void CreateGammaMatrix() const;

    /**ret[0] = max thread count, ret[1,2,3] = max thread for x,y,z per block*/
    static TArray<UINT> GetMaxThreadCountAndThreadPerblock();


    UINT m_ConstIntegers[kContentLength];
    Real m_ConstFloats[kContentLength];
};

__END_NAMESPACE


#endif //#ifndef _CUDAHELPER_H_

//=============================================================================
// END OF FILE
//=============================================================================