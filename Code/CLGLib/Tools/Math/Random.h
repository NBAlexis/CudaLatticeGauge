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
#define PI (3.141592653589f)
// - 1/4294967296UL
#define AM (0.00000000023283064365386963f)
// - sqrt(2)
#define SQRT2 (1.4142135623730951f)
// - 1 / sqrt(2), or sqrt(2)/2
#define InvSqrt2 (0.7071067811865475f)
// - 2.0f * PI
#define PI2 (6.283185307179586f)

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
    virtual __device__ FLOAT _deviceRandomF(UINT fatIndex) = 0;
    virtual __device__ UINT _deviceRandomUI(UINT fatIndex) = 0;
    __device__ __inline__ FLOAT _deviceRandomGaussF(UINT fatIndex)
    {
        FLOAT f1 = _deviceRandomF(fatIndex);
        FLOAT f2 = _deviceRandomF(fatIndex) * PI2;

        FLOAT amplitude = sqrt(-log(1.0f - f1) * 2.0f) * InvSqrt2;
        return cos(f2) * amplitude;
    }

    __device__ __inline__ cuComplex _deviceRandomGaussC(UINT fatIndex)
    {
        FLOAT f1 = _deviceRandomF(fatIndex);
        FLOAT f2 = _deviceRandomF(fatIndex) * PI2;

        FLOAT amplitude = sqrt(-log(1.0f - f1) * 2.0f) * InvSqrt2;
        return make_cuComplex(cos(f2) * amplitude, sin(f2) * amplitude);
    }
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

    __device__ __forceinline__ FLOAT _deviceRandomF(UINT fatIndex)
    {
        return AM * _deviceRandomUI(fatIndex);
    }

    __device__ __forceinline__ UINT _deviceRandomUI(UINT fatIndex)
    {
        m_pDeviceSeedTable[fatIndex] = (UINT)((1664525UL * m_pDeviceSeedTable[fatIndex] + 1013904223UL) & 0xffffffff);
        return m_pDeviceSeedTable[fatIndex];
    }

    /**
    * Although in bridge++, it says the deviation is 1/sqrt(2)
    * In fact, the standard deviation of it is 1
    */
    __device__ __forceinline__ FLOAT _deviceRandomGaussF(UINT fatIndex)
    {
        FLOAT f1 = _deviceRandomF(fatIndex);
        FLOAT f2 = _deviceRandomF(fatIndex) * PI2;

        FLOAT amplitude = sqrt(-log(1.0f - f1) * 2.0f) * InvSqrt2;
        return cos(f2) * amplitude;
    }

    __device__ __forceinline__ cuComplex _deviceRandomGaussC(UINT fatIndex)
    {
        FLOAT f1 = _deviceRandomF(fatIndex);
        FLOAT f2 = _deviceRandomF(fatIndex) * PI2;

        FLOAT amplitude = sqrt(-log(1.0f - f1) * 2.0f) * InvSqrt2;
        return make_cuComplex(cos(f2) * amplitude, sin(f2) * amplitude);
    }

    /**
    * run on device, parally set the table
    */
    __device__ __forceinline__ static void _deviceAsignSeeds(UINT* devicePtr, UINT uiSeed, UINT uiFatIndex)
    {
        devicePtr[uiFatIndex] = (1664525UL * (uiFatIndex^uiSeed) + 1013904223UL) & 0xffffffff;
    }
};

__DefineRandomFuncion(FLOAT, F)

__DefineRandomFuncion(UINT, UI)

__DefineRandomFuncion(FLOAT, GaussF)

__DefineRandomFuncion(cuComplex, GaussC)

__END_NAMESPACE

#endif //#ifndef _RANDOM_H_

//=============================================================================
// END OF FILE
//=============================================================================
