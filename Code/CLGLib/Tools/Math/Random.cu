//=============================================================================
// FILENAME : Random.h
// 
// DESCRIPTION:
//
//
// REVISION:
//  [12/6/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

///* include MTGP host helper functions */
//#include <curand_mtgp32_host.h>
///* include MTGP pre-computed parameter sets */
//#include <curand_mtgp32dc_p_11213.h>

__BEGIN_NAMESPACE

__global__ 
void _kernalAllocateSeedTable(UINT* pDevicePtr)
{
    intokernal;

    UINT uiSeed = _DC_Seed;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir + 1; ++idir)
        {
            UINT fatIndex = _deviceGetFatIndex(coord, idir);
            CRandom::_deviceAsignSeeds(pDevicePtr, uiSeed, fatIndex);
        }
    }
}

__global__
void _kernalInitialXORWOW(curandState * states)
{
    intokernal;

    UINT uiSeed = _DC_Seed;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir + 1; ++idir)
        {
            UINT fatIndex = _deviceGetFatIndex(coord, idir);
            curand_init(uiSeed, fatIndex, 0, &states[fatIndex]);
        }
    }
}

__global__
void _kernalInitialPhilox(curandStatePhilox4_32_10_t * states)
{
    intokernal;

    UINT uiSeed = _DC_Seed;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir + 1; ++idir)
        {
            UINT fatIndex = _deviceGetFatIndex(coord, idir);
            curand_init(uiSeed, fatIndex, 0, &states[fatIndex]);
        }
    }
}

__global__
void _kernalInitialMRG(curandStateMRG32k3a  * states)
{
    intokernal;

    UINT uiSeed = _DC_Seed;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir + 1; ++idir)
        {
            UINT fatIndex = _deviceGetFatIndex(coord, idir);
            curand_init(uiSeed, fatIndex, 0, &states[fatIndex]);
        }
    }
}

__global__
void _kernalInitialSobel32(curandStateSobol32* states, curandDirectionVectors32_t* dirs)
{
    intokernal;

    UINT uiSeed = _DC_Seed;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        UINT siteIndex = _deviceGetSiteIndex(coord);
        curand_init(dirs[siteIndex], uiSeed % 16, &states[siteIndex]);
    }
}

__global__
void _kernalInitialScrambledSobel32(curandStateScrambledSobol32* states, UINT* consts, curandDirectionVectors32_t* dirs)
{
    intokernal;

    UINT uiSeed = _DC_Seed;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        UINT siteIndex = _deviceGetSiteIndex(coord);
        curand_init(dirs[siteIndex], consts[siteIndex], uiSeed % __SOBEL_OFFSET_MAX, &states[siteIndex]);
    }
}

CRandom::~CRandom()
{

    switch (m_eRandomType)
    {
    case ER_Schrage:
        {
            checkCudaErrors(cudaFree(m_pDeviceSeedTable));
        }
        break;
    case ER_MRG32K3A:
        {
            CURAND_CALL(curandDestroyGenerator(m_HGen));
            checkCudaErrors(cudaFree(m_deviceBuffer));
            checkCudaErrors(cudaFree(m_pDeviceRandStatesMRG));
        }
        break;
    case ER_PHILOX4_32_10:
        {
            CURAND_CALL(curandDestroyGenerator(m_HGen));
            checkCudaErrors(cudaFree(m_deviceBuffer));
            checkCudaErrors(cudaFree(m_pDeviceRandStatesPhilox));
        }
        break;
    case ER_QUASI_SOBOL32:
        {
            CURAND_CALL(curandDestroyGenerator(m_HGen));
            checkCudaErrors(cudaFree(m_deviceBuffer));
            checkCudaErrors(cudaFree(m_pDeviceRandStatesSobol32));
            checkCudaErrors(cudaFree(m_pDeviceSobolDirVec));
        }
        break;
    case ER_SCRAMBLED_SOBOL32:
        {
            CURAND_CALL(curandDestroyGenerator(m_HGen));
            checkCudaErrors(cudaFree(m_deviceBuffer));
            checkCudaErrors(cudaFree(m_pDeviceRandStatesScrambledSobol32));
            checkCudaErrors(cudaFree(m_pDeviceSobolDirVec));
            checkCudaErrors(cudaFree(m_pDeviceSobelConsts));
        }
        break;
    case ER_XORWOW:
        default:
        {
            CURAND_CALL(curandDestroyGenerator(m_HGen));
            checkCudaErrors(cudaFree(m_deviceBuffer));
            checkCudaErrors(cudaFree(m_pDeviceRandStatesXORWOW));
        }
        break;
    }
}

void CRandom::InitialStatesXORWOW(UINT )
{
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceRandStatesXORWOW, sizeof(curandState) * _HC_Volumn * (_HC_Dir + 1)));
    preparethread;
    _kernalInitialXORWOW << <block, threads >> > (m_pDeviceRandStatesXORWOW);
}

void CRandom::InitialStatesPhilox(UINT )
{
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceRandStatesPhilox, sizeof(curandStatePhilox4_32_10_t) * _HC_Volumn * (_HC_Dir + 1)));
    preparethread;
    _kernalInitialPhilox << <block, threads >> > (m_pDeviceRandStatesPhilox);
}

void CRandom::InitialStatesMRG(UINT )
{
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceRandStatesMRG, sizeof(curandStateMRG32k3a) * _HC_Volumn * (_HC_Dir + 1)));
    preparethread;
    _kernalInitialMRG << <block, threads >> > (m_pDeviceRandStatesMRG);
}

void CRandom::InitialStatesSobol32(UINT )
{
    //support only 20000 dimensions, so using _HC_Volumn instead
    m_uiFatIdDivide = _HC_Dir + 1;
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceRandStatesSobol32, 
        sizeof(curandStateSobol32) * _HC_Volumn));
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceSobolDirVec, 
        sizeof(curandDirectionVectors32_t) * _HC_Volumn));

    //int[32]
    curandDirectionVectors32_t *hostVectors32;
    CURAND_CALL(curandGetDirectionVectors32(&hostVectors32, CURAND_DIRECTION_VECTORS_32_JOEKUO6));
    checkCudaErrors(cudaMemcpy(m_pDeviceSobolDirVec, hostVectors32, 
        _HC_Volumn * sizeof(curandDirectionVectors32_t),
        cudaMemcpyHostToDevice));

    preparethread;
    _kernalInitialSobel32 << <block, threads >> > (m_pDeviceRandStatesSobol32, m_pDeviceSobolDirVec);
}

void CRandom::InitialStatesScrambledSobol32(UINT uiSeed)
{
    m_uiFatIdDivide = _HC_Dir + 1;
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceRandStatesScrambledSobol32, 
        sizeof(curandStateScrambledSobol32) * _HC_Volumn));
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceSobolDirVec, 
        sizeof(curandDirectionVectors32_t) * _HC_Volumn));
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceSobelConsts, 
        sizeof(UINT) * _HC_Volumn));

    curandDirectionVectors32_t *hostVectors32;
    CURAND_CALL(curandGetDirectionVectors32(&hostVectors32, CURAND_SCRAMBLED_DIRECTION_VECTORS_32_JOEKUO6));
    checkCudaErrors(cudaMemcpy(
        m_pDeviceSobolDirVec, 
        hostVectors32, 
        _HC_Volumn * sizeof(curandDirectionVectors32_t),
        cudaMemcpyHostToDevice));

    UINT * hostScrambleConstants32;
    CURAND_CALL(curandGetScrambleConstants32(&hostScrambleConstants32));
    checkCudaErrors(cudaMemcpy(
        m_pDeviceSobelConsts, 
        hostScrambleConstants32, 
        _HC_Volumn * sizeof(UINT), 
        cudaMemcpyHostToDevice));

    preparethread;
    _kernalInitialScrambledSobel32 << <block, threads >> > (m_pDeviceRandStatesScrambledSobol32, m_pDeviceSobelConsts, m_pDeviceSobolDirVec);
}

void CRandom::InitialTableSchrage(UINT uiSeed)
{
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceSeedTable, sizeof(UINT) * _HC_Volumn * (_HC_Dir + 1)));
    preparethread;
    _kernalAllocateSeedTable << <block, threads >> > (m_pDeviceSeedTable);
}

Real GetRandomReal()
{
    return appGetLattice()->m_pRandom->GetRandomF();
}

#pragma region Test

__global__ void _kernelMCPi(UINT* output, UINT lengthyz, UINT lengthz, UINT uiLoop, UINT uithreadCount)
{
    __shared__ UINT sData1[1024];
    //__shared__ UINT sData2[1024];
    UINT uiToAdd = 0;
    //UINT uiToAdd2 = 0;
    UINT fatIndex = threadIdx.x * lengthyz + threadIdx.y * lengthz + threadIdx.z;
    for (UINT i = 0; i < uiLoop; ++i)
    {
        Real x = _deviceRandomF(fatIndex) * 2.0f - 1.0f;
        Real y = _deviceRandomF(fatIndex) * 2.0f - 1.0f;
        if (x * x + y * y < 1.0f)
        {
            ++uiToAdd;
        }
        //++uiToAdd2;
    }
    sData1[fatIndex] = uiToAdd;
    //sData2[fatIndex] = uiToAdd2;
    //if (0 == threadIdx.x)
    //{
    //    sData1[0] = 0;
    //    sData1[1] = 0;
    //}

    __syncthreads();
    if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {
        UINT all1 = 0;
        //UINT all2 = 0;
        for (UINT i = 0; i < uithreadCount; ++i)
        {
            all1 += sData1[i];
            //all2 += sData2[i];
        }
        //printf("how many?= %d\n", all1);
        atomicAdd(output, all1);
        //atomicAdd(output + 1, all2);
    }
}

__global__ void _kernelMCE(Real* output, UINT lengthyz, UINT lengthz, UINT uiLoop, UINT uithreadCount)
{
    __shared__ Real sData1[1024];
    __shared__ Real sData2[1024];
    Real fToAdd = 0;
    Real fToAdd2 = 0;
    UINT fatIndex = threadIdx.x * lengthyz + threadIdx.y * lengthz + threadIdx.z;
    for (UINT i = 0; i < uiLoop; ++i)
    {
        _Complex c = _deviceRandomGaussC(fatIndex);
        fToAdd += (c.x + c.y);
        fToAdd2 += (c.x * c.x + c.y * c.y);
    }
    sData1[fatIndex] = fToAdd;
    sData2[fatIndex] = fToAdd2;

    __syncthreads();
    if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {
        Real all1 = 0;
        Real all2 = 0;
        for (UINT i = 0; i < uithreadCount; ++i)
        {
            all1 += sData1[i];
            all2 += sData2[i];
        }
        //printf("how many?= %d\n", all1);
        atomicAdd(output, all1);
        atomicAdd(output + 1, all2);
    }
}

Real CLGAPI CalculatePi(const TArray<UINT> & decompose)
{
    dim3 blocknumber(decompose[0], decompose[1], decompose[2]);
    dim3 threadnumber(decompose[3], decompose[4], decompose[5]);
    UINT threadCount = decompose[3] * decompose[4] * decompose[5];
    UINT lengthyz = decompose[4] * decompose[5];
    UINT lengthz = decompose[5];
    UINT total = decompose[0] * decompose[1] * decompose[2] * decompose[3] * decompose[4] * decompose[5] * decompose[6];
    UINT uiLoop = decompose[6];

    UINT outPutHost[2];
    outPutHost[0] = 0;
    outPutHost[1] = 0;

    UINT *outPut;
    checkCudaErrors(cudaMalloc((void**)&outPut, sizeof(UINT) * 2));
    checkCudaErrors(cudaMemcpy(outPut, outPutHost, sizeof(UINT) * 2, cudaMemcpyHostToDevice));

    _kernelMCPi << <blocknumber, threadnumber >> > (outPut, lengthyz, lengthz, uiLoop, threadCount);
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaMemcpy(outPutHost, outPut, sizeof(UINT) * 2, cudaMemcpyDeviceToHost));

    appParanoiac(_T("==== results: %d / %d \n"), outPutHost[0], outPutHost[1]);

    return 4.0f * outPutHost[0] / (Real)(total);
}

Real CLGAPI CalculateE(const TArray<UINT> & decompose)
{
    dim3 blocknumber(decompose[0], decompose[1], decompose[2]);
    dim3 threadnumber(decompose[3], decompose[4], decompose[5]);
    UINT threadCount = decompose[3] * decompose[4] * decompose[5];
    UINT lengthyz = decompose[4] * decompose[5];
    UINT lengthz = decompose[5];
    UINT total = decompose[0] * decompose[1] * decompose[2] * decompose[3] * decompose[4] * decompose[5] * decompose[6];
    UINT uiLoop = decompose[6];

    Real outPutHost[2];
    outPutHost[0] = 0.0F;
    outPutHost[1] = 0.0F;

    Real *outPut;
    checkCudaErrors(cudaMalloc((void**)&outPut, sizeof(Real) * 2));
    checkCudaErrors(cudaMemcpy(outPut, outPutHost, sizeof(Real) * 2, cudaMemcpyHostToDevice));

    _kernelMCE << <blocknumber, threadnumber >> > (outPut, lengthyz, lengthz, uiLoop, threadCount);
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaMemcpy(outPutHost, outPut, sizeof(Real) * 2, cudaMemcpyDeviceToHost));

    Real fAv = outPutHost[0] / (2.0f * total);
    Real fBv = outPutHost[1] / (2.0f * total) - fAv * fAv;

    return _sqrt(fBv);
}

#pragma endregion

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
