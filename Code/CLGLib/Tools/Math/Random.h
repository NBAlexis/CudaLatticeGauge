//=============================================================================
// FILENAME : Random.h
// 
// DESCRIPTION:
// This is random number for parallel
//
//
// REVISION:
//  [12/6/2018 nbale]
//=============================================================================

#ifndef _RANDOM_H_
#define _RANDOM_H_

//save some constant memory of cuda?
#define PI ((Real)3.141592653589)
// - 1/4294967296UL
#define AM ((Real)0.00000000023283064365386963)
// - _sqrt(2)
#define SQRT2 ((Real)1.4142135623730951)
// - 1 / _sqrt(2), or _sqrt(2)/2
#define InvSqrt2 ((Real)0.7071067811865475)
// - 2.0f * PI
#define PI2 ((Real)6.283185307179586)

#define __DefineRandomFuncion(rettype, funcname) __device__ __inline__ static rettype _deviceRandom##funcname(UINT uiFatIndex) \
{ \
    return (1 == _constIntegers[ECI_UsingSchrageRandom]) ? __rs->_deviceRandom##funcname(uiFatIndex) : __r->_deviceRandom##funcname(uiFatIndex); \
}


__BEGIN_NAMESPACE

__DEFINE_ENUM (ERandom,
    ER_Schrage,
    ER_ForceDWORD = 0x7fffffff,
    )

class CLGAPI CRandom
{
public:
    CRandom(UINT uiSeed) {}
    virtual __device__ Real _deviceRandomF(UINT fatIndex) = 0;
    __device__ __inline__ Real _deviceRandomGaussF(UINT fatIndex)
    {
        Real f1 = _deviceRandomF(fatIndex);
        Real f2 = _deviceRandomF(fatIndex) * PI2;

        Real amplitude = _sqrt(-_log(1.0f - f1) * 2.0f) * InvSqrt2;
        return cos(f2) * amplitude;
    }

    __device__ __inline__ _Complex _deviceRandomGaussC(UINT fatIndex)
    {
        Real f1 = _deviceRandomF(fatIndex);
        Real f2 = _deviceRandomF(fatIndex) * PI2;

        Real amplitude = _sqrt(-_log(1.0f - f1) * 2.0f) * InvSqrt2;
        return _make_cuComplex(_cos(f2) * amplitude, _sin(f2) * amplitude);
    }

    __host__ virtual Real GetRandomF() = 0;

protected:

    virtual __device__ UINT _deviceRandomUI(UINT fatIndex) = 0;
    virtual __host__ UINT GetRandomUI() = 0;
};

/**
* we cannot make it virtual inline
* so we have a random for testing usage only
* this random is a special random, the others should inherent from CRandom
* for real usage, (if that is ok), one do not need a table. The table only ensures, same seed same result.
*/
class CLGAPI CRandomSchrage
{
public:
    CRandomSchrage(UINT uiSeed);
    ~CRandomSchrage();

    UINT* m_pDeviceSeedTable;

    __device__ __inline__ Real _deviceRandomF(UINT fatIndex)
    {
        return AM * _deviceRandomUI(fatIndex);
    }

    /**
    * Although in bridge++, it says the deviation is 1/_sqrt(2)
    * In fact, the standard deviation of it is 1
    */
    __device__ __inline__ Real _deviceRandomGaussF(UINT fatIndex)
    {
        Real f1 = _deviceRandomF(fatIndex);
        Real f2 = _deviceRandomF(fatIndex) * PI2;

        Real amplitude = _sqrt(-_log(1.0f - f1) * 2.0f) * InvSqrt2;
        return cos(f2) * amplitude;
    }

    __device__ __inline__ _Complex _deviceRandomGaussC(UINT fatIndex)
    {
        Real f1 = _deviceRandomF(fatIndex);
        Real f2 = _deviceRandomF(fatIndex) * PI2;

        Real amplitude = _sqrt(-_log(1.0f - f1) * 2.0f) * InvSqrt2;
        return _make_cuComplex(_cos(f2) * amplitude, _sin(f2) * amplitude);
    }

    /**
    * run on device, parally set the table
    */
    __device__ __inline__ static void _deviceAsignSeeds(UINT* devicePtr, UINT uiSeed, UINT uiFatIndex)
    {
        devicePtr[uiFatIndex] = (1664525UL * (uiFatIndex+uiSeed) + 1013904223UL) & 0xffffffff;
    }

    __host__ __inline__ Real GetRandomF()
    {
        return AM * GetRandomUI();
    }

protected:

    __device__ __inline__ UINT _deviceRandomUI(UINT fatIndex)
    {
        m_pDeviceSeedTable[fatIndex] = ((1664525UL * m_pDeviceSeedTable[fatIndex] + 1013904223UL) & 0xffffffff);
        return m_pDeviceSeedTable[fatIndex];
    }

    __host__ __inline__ UINT GetRandomUI()
    {
        m_uiHostSeed = (1664525UL * m_uiHostSeed + 1013904223UL) & 0xffffffff;
        return m_uiHostSeed;
    }

    UINT m_uiHostSeed;
};

__DefineRandomFuncion(Real, F)

__DefineRandomFuncion(Real, GaussF)

__DefineRandomFuncion(_Complex, GaussC)

extern CLGAPI Real GetRandomReal();

//==========================
//functions for test
extern Real CLGAPI CalculatePi(const TArray<UINT> & decompose);

extern Real CLGAPI CalculateE(const TArray<UINT> & decompose);

__END_NAMESPACE

#endif //#ifndef _RANDOM_H_

//=============================================================================
// END OF FILE
//=============================================================================
