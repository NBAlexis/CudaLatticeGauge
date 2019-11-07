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


#define __SOBEL_OFFSET_MAX (4096)

#define __DefineRandomFuncion(rettype, funcname) __device__ __inline__ static rettype _deviceRandom##funcname(UINT uiFatIndex) \
{ \
    return __r->_deviceRandom##funcname(uiFatIndex); \
}


#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
    return;}} while(0)

__BEGIN_NAMESPACE

__DEFINE_ENUM (ERandom,
    ER_Schrage,

    ER_XORWOW,
    ER_MRG32K3A,
    //ER_MTGP32, //see the document, this has lots of constraints.
    ER_PHILOX4_32_10,
    ER_QUASI_SOBOL32,
    ER_SCRAMBLED_SOBOL32,

    ER_ForceDWORD = 0x7fffffff,
    )

__DEFINE_ENUM (ERandomSeedType,
    ERST_Number,
    ERST_Timestamp,

    ERST_ForceDWORD = 0x7fffffff,
    )

class CLGAPI CRandom
{
public:

    /**
    * There are nine types of random number generators in cuRAND, that fall into two categories. 
    *  CURAND_RNG_PSEUDO_XORWOW, CURAND_RNG_PSEUDO_MRG32K3A, CURAND_RNG_PSEUDO_MTGP32, CURAND_RNG_PSEUDO_PHILOX4_32_10 and CURAND_RNG_PSEUDO_MT19937 are pseudorandom number generators. 
    * CURAND_RNG_PSEUDO_XORWOW is implemented using the XORWOW algorithm, a member of the xor-shift family of pseudorandom number generators. 
    * CURAND_RNG_PSEUDO_MRG32K3A is a member of the Combined Multiple Recursive family of pseudorandom number generators. 
    *  CURAND_RNG_PSEUDO_MT19937 and CURAND_RNG_PSEUDO_MTGP32 are members of the Mersenne Twister family of pseudorandom number generators. 
    * CURAND_RNG_PSEUDO_MTGP32 has parameters customized for operation on the GPU. 
    * CURAND_RNG_PSEUDO_MT19937 has the same parameters as CPU version, but ordering is different. 
    * CURNAD_RNG_PSEUDO_MT19937 supports only HOST API and can be used only on architecture sm_35 or higher. 
    * CURAND_RNG_PHILOX4_32_10 is a member of Philox family, which is one of the three non-cryptographic Counter Based Random Number Generators presented on SC11 conference by D E Shaw Research. 
    *
    *  There are 4 variants of the basic SOBOL¡¯ quasi random number generator. All of the variants generate sequences in up to 20,000 dimensions. CURAND_RNG_QUASI_SOBOL32, CURAND_RNG_QUASI_SCRAMBLED_SOBOL32, CURAND_RNG_QUASI_SOBOL64, and CURAND_RNG_QUASI_SCRAMBLED_SOBOL64 are quasirandom number generator types. 
    * CURAND_RNG_QUASI_SOBOL32 is a Sobol¡¯ generator of 32-bit sequences. 
    * CURAND_RNG_QUASI_SCRAMBLED_SOBOL32 is a scrambled Sobol¡¯ generator of 32-bit sequences. 
    * CURAND_RNG_QUASI_SOBOL64 is a Sobol¡¯ generator of 64-bit sequences. 
    * CURAND_RNG_QUASI_SCRAMBLED_SOBOL64 is a scrambled Sobol¡¯ generator of 64-bit sequences.
    */
    CRandom(UINT uiSeed, ERandom er) 
        : m_eRandomType(er)
        , m_uiFatIdDivide(1)
        , m_uiHostSeed(uiSeed)
    { 
        switch (er)
        {
            case ER_Schrage:
                {
                    InitialTableSchrage(uiSeed);
                }
                break;
            case ER_MRG32K3A:
                {
                    CURAND_CALL(curandCreateGenerator(&m_HGen, CURAND_RNG_PSEUDO_MRG32K3A));
                    CURAND_CALL(curandSetPseudoRandomGeneratorSeed(m_HGen, uiSeed));
                    checkCudaErrors(cudaMalloc((void**)&m_deviceBuffer, sizeof(FLOAT)));
                    InitialStatesMRG(uiSeed);
                }
                break;
            case ER_PHILOX4_32_10:
                {
                    CURAND_CALL(curandCreateGenerator(&m_HGen, CURAND_RNG_PSEUDO_PHILOX4_32_10));
                    CURAND_CALL(curandSetPseudoRandomGeneratorSeed(m_HGen, uiSeed));
                    checkCudaErrors(cudaMalloc((void**)&m_deviceBuffer, sizeof(FLOAT)));
                    InitialStatesPhilox(uiSeed);
                }
                break;
            case ER_QUASI_SOBOL32:
                {
                    //for sobol, on the host, we use XORWOW
                    CURAND_CALL(curandCreateGenerator(&m_HGen, CURAND_RNG_QUASI_SOBOL32));
                    checkCudaErrors(cudaMalloc((void**)&m_deviceBuffer, sizeof(FLOAT)));
                    InitialStatesSobol32(uiSeed);
                }
                break;
            case ER_SCRAMBLED_SOBOL32:
                {
                    //for sobol, on the host, we use XORWOW
                    CURAND_CALL(curandCreateGenerator(&m_HGen, CURAND_RNG_QUASI_SCRAMBLED_SOBOL32));
                    checkCudaErrors(cudaMalloc((void**)&m_deviceBuffer, sizeof(FLOAT)));
                    InitialStatesScrambledSobol32(uiSeed);
                }
                break;
            case ER_XORWOW:
            default:
                {
                    CURAND_CALL(curandCreateGenerator(&m_HGen, CURAND_RNG_PSEUDO_XORWOW));
                    CURAND_CALL(curandSetPseudoRandomGeneratorSeed(m_HGen, uiSeed));
                    checkCudaErrors(cudaMalloc((void**)&m_deviceBuffer, sizeof(FLOAT)));
                    InitialStatesXORWOW(uiSeed);
                }
                break;
        }

        checkCudaErrors(cudaGetLastError());
    }

    ~CRandom();

    /**
    * Note that this gives [0, 1), and curand_uniform gives (0, 1]
    */
    __device__ __inline__ Real _deviceRandomF(UINT fatIndex) const
    {
        switch (m_eRandomType)
        {
            case ER_Schrage:
                return AM * _deviceRandomUISchrage(fatIndex);
            case ER_MRG32K3A:
                return 1 - curand_uniform(&(m_pDeviceRandStatesMRG[fatIndex]));
            case ER_PHILOX4_32_10:
                return 1 - curand_uniform(&(m_pDeviceRandStatesPhilox[fatIndex]));
            case ER_QUASI_SOBOL32:
                return 1 - curand_uniform(&(m_pDeviceRandStatesSobol32[(fatIndex / m_uiFatIdDivide)]));
            case ER_SCRAMBLED_SOBOL32:
                return 1 - curand_uniform(&(m_pDeviceRandStatesScrambledSobol32[fatIndex / m_uiFatIdDivide]));
            case ER_XORWOW:
            default:
                return 1 - curand_uniform(&(m_pDeviceRandStatesXORWOW[fatIndex / m_uiFatIdDivide]));
        }

        //return 0;
    }

    /**
    * Although in bridge++, it says the deviation is 1/_sqrt(2)
    * In fact, the standard deviation of it is 1
    */
    __device__ __inline__ Real _deviceRandomGaussF(UINT fatIndex) const
    {
        const Real f1 = _deviceRandomF(fatIndex);
        const Real f2 = _deviceRandomF(fatIndex) * PI2;

        const Real oneMinusf1 = F(1.0) - f1;
        const Real inSqrt = -F(2.0) * _log(oneMinusf1 > F(0.0) ? oneMinusf1 : (_CLG_FLT_MIN));
        const Real amplitude = (inSqrt > F(0.0) ? _sqrt(inSqrt) : F(0.0)) * InvSqrt2;
        return _cos(f2) * amplitude;
    }

    __device__ __inline__ Real _deviceRandomGaussFSqrt2(UINT fatIndex) const
    {
        const Real f1 = _deviceRandomF(fatIndex);
        const Real f2 = _deviceRandomF(fatIndex) * PI2;

        const Real oneMinusf1 = F(1.0) - f1;
        const Real inSqrt = -F(2.0) * _log(oneMinusf1 > F(0.0) ? oneMinusf1 : (_CLG_FLT_MIN));
        const Real amplitude = (inSqrt > F(0.0) ? _sqrt(inSqrt) : F(0.0)) * F(0.5);
        return _cos(f2) * amplitude;
    }

    __device__ __inline__ CLGComplex _deviceRandomGaussC(UINT fatIndex) const
    {
        const Real f1 = _deviceRandomF(fatIndex);
        const Real f2 = _deviceRandomF(fatIndex) * PI2;

        const Real oneMinusf1 = F(1.0) - f1;
        const Real inSqrt = -F(2.0) * _log(oneMinusf1 > F(0.0) ? oneMinusf1 : (_CLG_FLT_MIN));
        const Real amplitude = (inSqrt > F(0.0) ? _sqrt(inSqrt) : F(0.0)) * InvSqrt2;
        return _make_cuComplex(_cos(f2) * amplitude, _sin(f2) * amplitude);
    }

    __device__ __inline__ CLGComplex _deviceRandomZ4(UINT fatIndex) const
    {
        const INT byRandom = _floor2int(F(4.0) * _deviceRandomF(fatIndex));

        if (0 == byRandom)
        {
            return _make_cuComplex(F(1.0), F(0.0));
        }
        else if (1 == byRandom)
        {
            return _make_cuComplex(F(0.0), F(1.0));
        }
        else if (2 == byRandom)
        {
            return _make_cuComplex(-F(1.0), F(0.0));
        }
        return _make_cuComplex(F(0.0), -F(1.0));
    }

    __host__ __inline__ Real GetRandomF()
    {
        if (ER_Schrage == m_eRandomType)
        {
            return AM * GetRandomUISchrage();
        }

        curandGenerateUniform(m_HGen, m_deviceBuffer, 1);
        checkCudaErrors(cudaMemcpy(m_hostBuffer, m_deviceBuffer, sizeof(FLOAT), cudaMemcpyDeviceToHost));
        return (Real)m_hostBuffer[0];
    }

    FLOAT* m_deviceBuffer;
    FLOAT m_hostBuffer[1];
    curandGenerator_t m_HGen;
    ERandom m_eRandomType;
    UINT m_uiFatIdDivide;

protected:

    void InitialStatesXORWOW(UINT uiSeed);
    void InitialStatesPhilox(UINT uiSeed);
    void InitialStatesMRG(UINT uiSeed);
    void InitialStatesSobol32(UINT uiSeed);
    void InitialStatesScrambledSobol32(UINT uiSeed);

    curandState* m_pDeviceRandStatesXORWOW;
    curandStatePhilox4_32_10_t* m_pDeviceRandStatesPhilox;
    curandStateMRG32k3a* m_pDeviceRandStatesMRG;

    curandStateSobol32* m_pDeviceRandStatesSobol32;
    curandDirectionVectors32_t* m_pDeviceSobolDirVec;
    UINT* m_pDeviceSobelConsts;
    curandStateScrambledSobol32* m_pDeviceRandStatesScrambledSobol32;

#pragma region Schrage

public:

    UINT* m_pDeviceSeedTable;

    /**
    * run on device, parally set the table
    */
    __device__ __inline__ static void _deviceAsignSeeds(UINT* devicePtr, UINT uiSeed, UINT uiFatIndex)
    {
        devicePtr[uiFatIndex] = (1664525UL * (uiFatIndex + uiSeed) + 1013904223UL) & 0xffffffff;
    }

protected:

    void InitialTableSchrage(UINT uiSeed);

    __device__ __inline__ UINT _deviceRandomUISchrage(UINT fatIndex) const
    {
        m_pDeviceSeedTable[fatIndex] = ((1664525UL * m_pDeviceSeedTable[fatIndex] + 1013904223UL) & 0xffffffff);
        return m_pDeviceSeedTable[fatIndex];
    }

    __host__ __inline__ UINT GetRandomUISchrage()
    {
        m_uiHostSeed = (1664525UL * m_uiHostSeed + 1013904223UL) & 0xffffffff;
        return m_uiHostSeed;
    }

    UINT m_uiHostSeed;

#pragma endregion

};

__DefineRandomFuncion(Real, F)

__DefineRandomFuncion(Real, GaussF)

__DefineRandomFuncion(Real, GaussFSqrt2)

__DefineRandomFuncion(CLGComplex, GaussC)

__DefineRandomFuncion(CLGComplex, Z4)

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
