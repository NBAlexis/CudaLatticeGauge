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

            CRandomSchrage::_deviceAsignSeeds(pDevicePtr, uiSeed, fatIndex);
        }
    }
}

extern "C" {
    void _callKernelInitialRandomTable(UINT* devicePtr)
    {
        preparethread;
        _kernalAllocateSeedTable << <block, threads >> > (devicePtr);
    }
}

CRandomSchrage::CRandomSchrage(UINT uiSeed)
{
    m_uiHostSeed = uiSeed;
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceSeedTable, sizeof(UINT) * _HC_Volumn * (_HC_Dir + 1)));
    _callKernelInitialRandomTable(m_pDeviceSeedTable);
}

CRandomSchrage::~CRandomSchrage()
{
    checkCudaErrors(cudaFree(m_pDeviceSeedTable));
}

Real GetRandomReal()
{
    return (1 == appGetCudaHelper()->m_ConstIntegers[ECI_UsingSchrageRandom])
        ? appGetLattice()->m_pRandomSchrage->GetRandomF()
        : appGetLattice()->m_pRandom->GetRandomF();
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
