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

__BEGIN_NAMESPACE

__global__ void _CLG_LAUNCH_BOUND
_kernalAllocateSeedTable(UINT* pDevicePtr)
{
    intokernaldir;

    const UINT uiSeed = _DC_Seed;

    for (UINT idir = 0; idir < uiDir + 1; ++idir)
    {
        UINT fatIndex = _deviceGetFatIndex(uiSiteIndex, idir);
        CRandom::_deviceAsignSeeds(pDevicePtr, uiSeed, fatIndex);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernalInitialXORWOW(curandState * states)
{
    UINT uiSiteIndex = ((threadIdx.x + blockIdx.x * blockDim.x) 
        * blockDim.y * gridDim.y * blockDim.z * gridDim.z 
        + (threadIdx.y + blockIdx.y * blockDim.y) 
        * blockDim.z * gridDim.z 
        + (threadIdx.z + blockIdx.z * blockDim.z));
    curand_init(_DC_Seed, uiSiteIndex, uiSiteIndex, &states[uiSiteIndex]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernalInitialPhilox(curandStatePhilox4_32_10_t * states)
{
    const UINT uiSiteIndex = ((threadIdx.x + blockIdx.x * blockDim.x) * blockDim.y * gridDim.y * blockDim.z * gridDim.z + (threadIdx.y + blockIdx.y * blockDim.y) * blockDim.z * gridDim.z + (threadIdx.z + blockIdx.z * blockDim.z));
    const UINT uiSeed = _DC_Seed;
    const UINT uiDir = _DC_Dir;
    for (UINT idir = 0; idir < uiDir + 1; ++idir)
    {
        UINT fatIndex = _deviceGetFatIndex(uiSiteIndex, idir);
        curand_init(uiSeed, fatIndex, 0, &states[fatIndex]);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernalInitialMRG(curandStateMRG32k3a  * states)
{
    const UINT uiSiteIndex = ((threadIdx.x + blockIdx.x * blockDim.x) * blockDim.y * gridDim.y * blockDim.z * gridDim.z + (threadIdx.y + blockIdx.y * blockDim.y) * blockDim.z * gridDim.z + (threadIdx.z + blockIdx.z * blockDim.z));
    const UINT uiSeed = _DC_Seed;
    const UINT uiDir = _DC_Dir;
    for (UINT idir = 0; idir < uiDir + 1; ++idir)
    {
        UINT fatIndex = _deviceGetFatIndex(uiSiteIndex, idir);
        curand_init(uiSeed, fatIndex, 0, &states[fatIndex]);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernalInitialSobel32(curandStateSobol32* states, curandDirectionVectors32_t* dirs)
{
    intokernal;
    curand_init(dirs[uiSiteIndex], _DC_Seed % 16, &states[uiSiteIndex]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernalInitialScrambledSobel32(curandStateScrambledSobol32* states, UINT* consts, curandDirectionVectors32_t* dirs)
{
    intokernal;
    curand_init(dirs[uiSiteIndex], consts[uiSiteIndex], _DC_Seed % __SOBEL_OFFSET_MAX, &states[uiSiteIndex]);
}

CRandom::~CRandom()
{

    switch (m_eRandomType)
    {
    case ER_Schrage:
        {
            checkCudaErrors(__cudaFree(m_pDeviceSeedTable));
        }
        break;
    case ER_MRG32K3A:
        {
            CURAND_CALL(curandDestroyGenerator(m_HGen));
            checkCudaErrors(__cudaFree(m_deviceBuffer));
            checkCudaErrors(__cudaFree(m_pDeviceRandStatesMRG));
        }
        break;
    case ER_PHILOX4_32_10:
        {
            CURAND_CALL(curandDestroyGenerator(m_HGen));
            checkCudaErrors(__cudaFree(m_deviceBuffer));
            checkCudaErrors(__cudaFree(m_pDeviceRandStatesPhilox));
        }
        break;
    case ER_QUASI_SOBOL32:
        {
            CURAND_CALL(curandDestroyGenerator(m_HGen));
            checkCudaErrors(__cudaFree(m_deviceBuffer));
            checkCudaErrors(__cudaFree(m_pDeviceRandStatesSobol32));
            checkCudaErrors(__cudaFree(m_pDeviceSobolDirVec));
        }
        break;
    case ER_SCRAMBLED_SOBOL32:
        {
            CURAND_CALL(curandDestroyGenerator(m_HGen));
            checkCudaErrors(__cudaFree(m_deviceBuffer));
            checkCudaErrors(__cudaFree(m_pDeviceRandStatesScrambledSobol32));
            checkCudaErrors(__cudaFree(m_pDeviceSobolDirVec));
            checkCudaErrors(__cudaFree(m_pDeviceSobelConsts));
        }
        break;
    case ER_XORWOW:
        default:
        {
            CURAND_CALL(curandDestroyGenerator(m_HGen));
            checkCudaErrors(__cudaFree(m_deviceBuffer));
            checkCudaErrors(__cudaFree(m_pDeviceRandStatesXORWOW));
        }
        break;
    }
}

//Initial XORWOW only support 512 threads per block
void CRandom::InitialStatesXORWOW(UINT )
{
    m_uiFatIdDivide = _HC_Dir + 1;
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceRandStatesXORWOW, sizeof(curandState) * _HC_Volume));
    TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock();
    deviceConstraints[0] = 512;
    TArray<UINT> latticeDim;
    latticeDim.AddItem(_HC_Lx * _HC_Ly);
    latticeDim.AddItem(_HC_Lz);
    latticeDim.AddItem(_HC_Lt);
    TArray <UINT> decomp = _getDecompose(deviceConstraints, latticeDim);
    dim3 block(decomp[0], decomp[1], decomp[2]);
    dim3 threads(decomp[3], decomp[4], decomp[5]);
    _kernalInitialXORWOW << <block, threads >> > (m_pDeviceRandStatesXORWOW);
}

//Initial Philox only support 256 threads per block
void CRandom::InitialStatesPhilox(UINT )
{
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceRandStatesPhilox, sizeof(curandStatePhilox4_32_10_t) * _HC_Volume * (_HC_Dir + 1)));

    TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock();
    deviceConstraints[0] = 256;
    TArray<UINT> latticeDim;
    latticeDim.AddItem(_HC_Lx * _HC_Ly);
    latticeDim.AddItem(_HC_Lz);
    latticeDim.AddItem(_HC_Lt);
    TArray <UINT> decomp = _getDecompose(deviceConstraints, latticeDim);
    dim3 block(decomp[0], decomp[1], decomp[2]);
    dim3 threads(decomp[3], decomp[4], decomp[5]);

    _kernalInitialPhilox << <block, threads >> > (m_pDeviceRandStatesPhilox);
}

//Initial MRG only support 256 threads per block
void CRandom::InitialStatesMRG(UINT )
{
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceRandStatesMRG, sizeof(curandStateMRG32k3a) * _HC_Volume * (_HC_Dir + 1)));
    TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock();
    deviceConstraints[0] = 256;
    TArray<UINT> latticeDim;
    latticeDim.AddItem(_HC_Lx * _HC_Ly);
    latticeDim.AddItem(_HC_Lz);
    latticeDim.AddItem(_HC_Lt);
    TArray <UINT> decomp = _getDecompose(deviceConstraints, latticeDim);
    dim3 block(decomp[0], decomp[1], decomp[2]);
    dim3 threads(decomp[3], decomp[4], decomp[5]);
    _kernalInitialMRG << <block, threads >> > (m_pDeviceRandStatesMRG);
}

void CRandom::InitialStatesSobol32(UINT )
{
    //support only 20000 dimensions, so using _HC_Volumn instead
    m_uiFatIdDivide = _HC_Dir + 1;
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceRandStatesSobol32,
        sizeof(curandStateSobol32) * _HC_Volume));
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceSobolDirVec,
        sizeof(curandDirectionVectors32_t) * _HC_Volume));

    //int[32]
    curandDirectionVectors32_t *hostVectors32;
    CURAND_CALL(curandGetDirectionVectors32(&hostVectors32, CURAND_DIRECTION_VECTORS_32_JOEKUO6));
    checkCudaErrors(cudaMemcpy(m_pDeviceSobolDirVec, hostVectors32, 
        _HC_Volume * sizeof(curandDirectionVectors32_t),
        cudaMemcpyHostToDevice));

    preparethread;
    _kernalInitialSobel32 << <block, threads >> > (m_pDeviceRandStatesSobol32, m_pDeviceSobolDirVec);
}

void CRandom::InitialStatesScrambledSobol32(UINT )
{
    m_uiFatIdDivide = _HC_Dir + 1;
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceRandStatesScrambledSobol32,
        sizeof(curandStateScrambledSobol32) * _HC_Volume));
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceSobolDirVec,
        sizeof(curandDirectionVectors32_t) * _HC_Volume));
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceSobelConsts,
        sizeof(UINT) * _HC_Volume));

    curandDirectionVectors32_t *hostVectors32;
    CURAND_CALL(curandGetDirectionVectors32(&hostVectors32, CURAND_SCRAMBLED_DIRECTION_VECTORS_32_JOEKUO6));
    checkCudaErrors(cudaMemcpy(
        m_pDeviceSobolDirVec, 
        hostVectors32, 
        _HC_Volume * sizeof(curandDirectionVectors32_t),
        cudaMemcpyHostToDevice));

    UINT * hostScrambleConstants32;
    CURAND_CALL(curandGetScrambleConstants32(&hostScrambleConstants32));
    checkCudaErrors(cudaMemcpy(
        m_pDeviceSobelConsts, 
        hostScrambleConstants32, 
        _HC_Volume * sizeof(UINT), 
        cudaMemcpyHostToDevice));

    preparethread;
    _kernalInitialScrambledSobel32 << <block, threads >> > (m_pDeviceRandStatesScrambledSobol32, m_pDeviceSobelConsts, m_pDeviceSobolDirVec);
}

void CRandom::InitialTableSchrage(UINT )
{
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceSeedTable, sizeof(UINT) * _HC_Volume * (_HC_Dir + 1)));
    preparethread;
    _kernalAllocateSeedTable << <block, threads >> > (m_pDeviceSeedTable);
}

UINT CRandom::DebugSeedTable() const
{
    UINT* hostseedtable = (UINT*)(malloc(sizeof(UINT) * _HC_Volume * (_HC_Dir + 1)));
    checkCudaErrors(cudaMemcpy(hostseedtable, m_pDeviceSeedTable, sizeof(UINT) * _HC_Volume * (_HC_Dir + 1), cudaMemcpyDeviceToHost));

    appGeneral(_T("%d, %d, %d, %d\n"), hostseedtable[0], hostseedtable[1], hostseedtable[2], hostseedtable[3]);

    UINT uiError = 0;
    if (static_cast<INT>(hostseedtable[0]) != 1566952946)
    {
        ++uiError;
    }
    if (static_cast<INT>(hostseedtable[1]) != -436264799)
    {
        ++uiError;
    }
    if (static_cast<INT>(hostseedtable[2]) != 846144630)
    {
        ++uiError;
    }
    if (static_cast<INT>(hostseedtable[3]) != 2128554059)
    {
        ++uiError;
    }

    appSafeFree(hostseedtable);
    return uiError;
}



Real GetRandomReal()
{
    return appGetLattice()->m_pRandom->GetRandomF();
}

#pragma region Test

__global__ void _CLG_LAUNCH_BOUND
_kernelMCPi(UINT* output, UINT lengthyz, UINT lengthz, UINT uiLoop, UINT uithreadCount)
{
    __shared__ UINT sData1[1024];
    __shared__ UINT sData2[1024];
    UINT uiToAdd = 0;
    UINT uiToAdd2 = 0;
    //We have a very large grid, but for a block, it is always smaller (or equval to volumn)
    const UINT fatIndex = threadIdx.x * lengthyz + threadIdx.y * lengthz + threadIdx.z;
    for (UINT i = 0; i < uiLoop; ++i)
    {
        const Real x = _deviceRandomF(fatIndex) * F(2.0) - F(1.0);
        const Real y = _deviceRandomF(fatIndex) * F(2.0) - F(1.0);
        if (x * x + y * y < F(1.0))
        {
            ++uiToAdd;
        }
        ++uiToAdd2;
    }
    sData1[fatIndex] = uiToAdd;
    sData2[fatIndex] = uiToAdd2;

    __syncthreads();
    if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {
        UINT all1 = 0;
        UINT all2 = 0;
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

__global__ void _CLG_LAUNCH_BOUND
_kernelMCE(Real* output, UINT lengthyz, UINT lengthz, UINT uiLoop, UINT uithreadCount)
{
    __shared__ Real sData1[1024];
    __shared__ Real sData2[1024];
    Real fToAdd = 0;
    Real fToAdd2 = 0;
    const UINT fatIndex = threadIdx.x * lengthyz + threadIdx.y * lengthz + threadIdx.z;
    for (UINT i = 0; i < uiLoop; ++i)
    {
        const CLGComplex c = _deviceRandomGaussC(fatIndex);
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
    const UINT threadCount = decompose[3] * decompose[4] * decompose[5];
    const UINT lengthyz = decompose[4] * decompose[5];
    const UINT lengthz = decompose[5];
    const UINT total = decompose[0] * decompose[1] * decompose[2] * decompose[3] * decompose[4] * decompose[5] * decompose[6];
    const UINT uiLoop = decompose[6];

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

    return F(4.0) * outPutHost[0] / (Real)(total);
}

Real CLGAPI CalculateE(const TArray<UINT> & decompose)
{
    dim3 blocknumber(decompose[0], decompose[1], decompose[2]);
    dim3 threadnumber(decompose[3], decompose[4], decompose[5]);
    const UINT threadCount = decompose[3] * decompose[4] * decompose[5];
    const UINT lengthyz = decompose[4] * decompose[5];
    const UINT lengthz = decompose[5];
    const UINT total = decompose[0] * decompose[1] * decompose[2] * decompose[3] * decompose[4] * decompose[5] * decompose[6];
    const UINT uiLoop = decompose[6];

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

    const Real fAv = outPutHost[0] / (2.0f * total);
    const Real fBv = outPutHost[1] / (2.0f * total) - fAv * fAv;

    return _hostsqrt(fBv);
}

#pragma endregion

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
