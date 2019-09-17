#include "CudaHelper.h"

#define _FFT_MAXDIMDIV 32
#define _FFT_MAXDIMDIVFACTORCOUNT 127
#define _FFT_MAXDIMDIVFACTOR 1021
#define _FFT_MAXDIM 4

#define _FFT_TESTVECTOR 10
#define _FFT_TESTMATRIXX 3
#define _FFT_TESTMATRIXY 4
#define _FFT_TESTMATRIXZ 4
#define _FFT_TESTMATRIXW 3
#define _FFT_TESTMATRIXV (_FFT_TESTMATRIXX * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ * _FFT_TESTMATRIXW)
#define _FFT_TESTMATRIXV3D (_FFT_TESTMATRIXX * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ)

#define _pi_ (3.141592653589f)

#pragma region divice functions

/**
 * DFT of the array
 * diviceSource[iStartS + n * iStrideS]
 * Write to diviceArrayRes[iStartT + n * iStrideT]
 * Tested
 */
__device__ void _deviceSmallDFT(
    bool bInverse, 
    cuComplex* diviceArrayRes,
    const cuComplex* __restrict__ diviceSource,
    int iN, 
    int iStartT, 
    int iStrideT, 
    int iStartS, 
    int iStrideS)
{
    for (int i = 0; i < iN; ++i)
    {
        cuComplex res = make_cuComplex(0.0f, 0.0f);
        for (int j = 0; j < iN; ++j)
        {
            float fAngle = 2.0f * _pi_ * i * j / iN;
            res = cuCaddf(res, cuCmulf(
                diviceSource[iStartS + iStrideS * j],
                make_cuComplex(cosf(fAngle), bInverse ? -sinf(fAngle) : sinf(fAngle))));
        }
        diviceArrayRes[iStartT + iStrideT * i] = res;
    }
}

#pragma endregion

#pragma region kernels

/**
 * test Small DFT
 */
__global__ void _kernelTestSmallDFT(
    cuComplex* diviceArrayRes, 
    const cuComplex* __restrict__ diviceSource,
    int iN)
{
    _deviceSmallDFT(false, diviceArrayRes, diviceSource, iN, 0, 1, 0, 1);
}

/**
 * 
 */
__global__ void CTMappingForward(
    const cuComplex* __restrict__ diviceSource, 
    const int* __restrict__ subDim,
    cuComplex* diviceArrayRes)
{
    
}

//__global__ void CopyData(
//    const cuComplex* __restrict__ diviceSource, 
//    cuComplex* res)
//{
//    
//}

#pragma endregion

/**
 * If iSize = 230
 * res = [4,  2, 3, 5, 7], with 4 factors, 2x3x5x7 = 230
 */
int FindDecomp(int* res, int iSize)
{
    bool bFound = true;
    int iFactorCount = 0;
    int iMaxFactor = 0;
    while (bFound)
    {
        bFound = false;
        for (int i = 2; i < iSize; ++i)
        {
            if ((iSize / i) * i == iSize)
            {
                ++iFactorCount;
                res[iFactorCount] = i;
                iSize = iSize / i;
                if (i > iMaxFactor)
                {
                    iMaxFactor = i;
                }
                bFound = true;
                break;
            }
        }

        if (1 == iSize)
        {
            break;
        }

        if (!bFound || iFactorCount >= _FFT_MAXDIMDIVFACTORCOUNT - 1)
        {
            ++iFactorCount;
            res[iFactorCount] = iSize;
            if (iSize > iMaxFactor)
            {
                iMaxFactor = iSize;
            }
            break;
        }
    }
    res[0] = iFactorCount;

    return iMaxFactor;
}

void GenerateTestArray(cuComplex* hostArray, int iSize)
{
    for (int i = 0; i < iSize; ++i)
    {
        hostArray[i] = make_cuComplex((rand() % 101 - 50) / 50.0f, (rand() % 101 - 50) / 50.0f);
    }
}

void GenerateTestArray4D(cuComplex* hostArray)
{
    for (int i = 0; i < _FFT_TESTMATRIXV; ++i)
    {
        hostArray[i] = make_cuComplex((rand() % 101 - 50) / 50.0f, (rand() % 101 - 50) / 50.0f);
    }
}

void PrintTestArray1D(cuComplex* hostArray)
{
    SaveLog("\n{");
    for (int i = 0; i < _FFT_TESTMATRIXX; ++i)
    {
        const int iIndex = i;
        SaveLog("%1.10f %s %1.10fi",
            hostArray[i].x,
            hostArray[i].y < 0.0f ? "-" : "+",
            abs(hostArray[i].y));

        if (i == _FFT_TESTMATRIXX - 1)
        {
            SaveLog("}\n");
        }
        else
        {
            SaveLog(",");
        }
    }
}

void PrintTestArray2D(cuComplex* hostArray)
{
    SaveLog("\n{");
    for (int i = 0; i < _FFT_TESTMATRIXX; ++i)
    {
        SaveLog("{");

        for (int j = 0; j < _FFT_TESTMATRIXY; ++j)
        {
            const int iIndex = i * _FFT_TESTMATRIXY + j;
            SaveLog("%1.10f %s %1.10f I",
                hostArray[iIndex].x,
                hostArray[iIndex].y < 0.0f ? "" : "+",
                hostArray[iIndex].y);

            if (j == _FFT_TESTMATRIXY - 1)
            {
                SaveLog("}");
            }
            else
            {
                SaveLog(",");
            }
        }

        if (i == _FFT_TESTMATRIXX - 1)
        {
            SaveLog("}\n");
        }
        else
        {
            SaveLog(",\n");
        }
    }
}

void PrintTestArray3D(cuComplex* hostArray)
{
    SaveLog("\n{");
    for (int i = 0; i < _FFT_TESTMATRIXX; ++i)
    {
        SaveLog("{");

        for (int j = 0; j < _FFT_TESTMATRIXY; ++j)
        {
            SaveLog("{");
            for (int k = 0; k < _FFT_TESTMATRIXZ; ++k)
            {
                const int iIndex = i * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ
                    + j * _FFT_TESTMATRIXZ
                    + k;
                SaveLog("%1.10f %s %1.10f I",
                    hostArray[iIndex].x,
                    hostArray[iIndex].y < 0.0f ? "" : "+",
                    hostArray[iIndex].y);

                if (k == _FFT_TESTMATRIXZ - 1)
                {
                    SaveLog("}");
                }
                else
                {
                    SaveLog(", ");
                }
            }
            if (j == _FFT_TESTMATRIXY - 1)
            {
                SaveLog("}\n");
            }
            else
            {
                SaveLog(",\n");
            }
        }

        if (i == _FFT_TESTMATRIXX - 1)
        {
            SaveLog("}\n");
        }
        else
        {
            SaveLog(",\n");
        }
    }
}

void PrintTestArray4D(cuComplex* hostArray)
{
    SaveLog("\n{\n");
    for (int i = 0; i < _FFT_TESTMATRIXX; ++i)
    {
        SaveLog("{");
        for (int j = 0; j < _FFT_TESTMATRIXY; ++j)
        {
            if (0 == j)
            {
                SaveLog("{");
            }
            else
            {
                SaveLog(" {");
            }
            for (int k = 0; k < _FFT_TESTMATRIXZ; ++k)
            {
                SaveLog("{");
                for (int l = 0; l < _FFT_TESTMATRIXW; ++l)
                {
                    const int iIndex = i * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ * _FFT_TESTMATRIXW
                                     + j * _FFT_TESTMATRIXZ * _FFT_TESTMATRIXW
                                     + k * _FFT_TESTMATRIXW
                                     + l;
                    SaveLog("%1.10f %s %1.10f I",
                        hostArray[iIndex].x,
                        hostArray[iIndex].y < 0.0f ? "-" : "+",
                        abs(hostArray[iIndex].y));
                    if (l == _FFT_TESTMATRIXW - 1)
                    {
                        SaveLog("}");
                    }
                    else
                    {
                        SaveLog(", ");
                    }
                }

                if (k == _FFT_TESTMATRIXZ - 1)
                {
                    SaveLog("}");
                }
                else
                {
                    SaveLog(", ");
                }
            }
            if (j == _FFT_TESTMATRIXY - 1)
            {
                SaveLog("\n}");
            }
            else
            {
                SaveLog(",\n");
            }
        }

        if (i == _FFT_TESTMATRIXX - 1)
        {
            SaveLog("}\n");
        }
        else
        {
            SaveLog(",\n");
        }
    }
}

cuComplex* GenerateTestMatrix(void)
{
    return nullptr;
}

/**
* 1D
* input [b * idist + x * istride]
* output[b * odist + x * ostride]
* 2D
* input [b * idist + (x * inembed[1] + y) * istride]
* output[b * odist + (x * onembed[1] + y) * ostride]
* 3D
* input [b * idist + ((x * inembed[1] + y) * inembed[2] + z) * istride]
* output[b * odist + ((x * onembed[1] + y) * onembed[2] + z) * ostride]
*/

int main()
{
    ClearLog();

    printf("decomp:");
    int decomp[128];
    FindDecomp(decomp, 210);
    for (int i = 0; i < decomp[0]; ++i)
    {
        printf("%d,", decomp[i + 1]);
    }
    printf("\n");


    cuComplex* dD1Res;
    cuComplex* dD1Source;
    cuComplex* dD2Res;
    cuComplex* dD2Source;
    cuComplex* dD3Res;
    cuComplex* dD3Source;
    cuComplex* dD4Res;
    cuComplex* dD4Source;
    cuComplex* hD1Res = (cuComplex*)malloc(_FFT_TESTMATRIXX * sizeof(cuComplex));
    cuComplex* hD1Source = (cuComplex*)malloc(_FFT_TESTMATRIXX * sizeof(cuComplex));
    cuComplex* hD2Res = (cuComplex*)malloc(_FFT_TESTMATRIXX * _FFT_TESTMATRIXY * sizeof(cuComplex));
    cuComplex* hD2Source = (cuComplex*)malloc(_FFT_TESTMATRIXX * _FFT_TESTMATRIXY * sizeof(cuComplex));
    cuComplex* hD3Res = (cuComplex*)malloc(_FFT_TESTMATRIXX * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ * sizeof(cuComplex));
    cuComplex* hD3Source = (cuComplex*)malloc(_FFT_TESTMATRIXX * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ * sizeof(cuComplex));
    cuComplex* hD4Res = (cuComplex*)malloc(_FFT_TESTMATRIXX * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ * _FFT_TESTMATRIXW * sizeof(cuComplex));
    cuComplex* hD4Source = (cuComplex*)malloc(_FFT_TESTMATRIXX * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ * _FFT_TESTMATRIXW * sizeof(cuComplex));

    checkCudaErrors(cudaMalloc((void**)&dD1Res, _FFT_TESTMATRIXX * sizeof(cuComplex)));
    checkCudaErrors(cudaMalloc((void**)&dD1Source, _FFT_TESTMATRIXX * sizeof(cuComplex)));
    checkCudaErrors(cudaMalloc((void**)&dD2Res, _FFT_TESTMATRIXX * _FFT_TESTMATRIXY * sizeof(cuComplex)));
    checkCudaErrors(cudaMalloc((void**)&dD2Source, _FFT_TESTMATRIXX * _FFT_TESTMATRIXY * sizeof(cuComplex)));
    checkCudaErrors(cudaMalloc((void**)&dD3Res, _FFT_TESTMATRIXX * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ * sizeof(cuComplex)));
    checkCudaErrors(cudaMalloc((void**)&dD3Source, _FFT_TESTMATRIXX * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ * sizeof(cuComplex)));
    checkCudaErrors(cudaMalloc((void**)&dD4Res, _FFT_TESTMATRIXX * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ * _FFT_TESTMATRIXW * sizeof(cuComplex)));
    checkCudaErrors(cudaMalloc((void**)&dD4Source, _FFT_TESTMATRIXX * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ * _FFT_TESTMATRIXW * sizeof(cuComplex)));

    GenerateTestArray(hD1Source, _FFT_TESTMATRIXX);
    GenerateTestArray(hD2Source, _FFT_TESTMATRIXX * _FFT_TESTMATRIXY);
    GenerateTestArray(hD3Source, _FFT_TESTMATRIXX * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ);
    GenerateTestArray(hD4Source, _FFT_TESTMATRIXX * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ * _FFT_TESTMATRIXW);

#pragma region 1D

    PrintTestArray1D(hD1Source);
    checkCudaErrors(cudaMemcpy(dD1Source, hD1Source, _FFT_TESTMATRIXX * sizeof(cuComplex), cudaMemcpyHostToDevice));
    cufftHandle plan1d;
    cufftPlan1d(&plan1d, _FFT_TESTMATRIXX, CUFFT_C2C, 1);
    cufftResult res1D = cufftExecC2C(plan1d, dD1Source, dD1Res, CUFFT_FORWARD);
    printf("1D res = %d\n", res1D);
    checkCudaErrors(cudaMemcpy(hD1Res, dD1Res, _FFT_TESTMATRIXX * sizeof(cuComplex), cudaMemcpyDeviceToHost));
    PrintTestArray1D(hD1Res);

#pragma endregion

#pragma region 2D

    PrintTestArray2D(hD2Source);
    checkCudaErrors(cudaMemcpy(dD2Source, hD2Source, _FFT_TESTMATRIXX * _FFT_TESTMATRIXY * sizeof(cuComplex), cudaMemcpyHostToDevice));
    cufftHandle plan2d;
    cufftPlan2d(&plan2d, _FFT_TESTMATRIXX, _FFT_TESTMATRIXY, CUFFT_C2C);
    cufftResult res2D = cufftExecC2C(plan2d, dD2Source, dD2Res, CUFFT_FORWARD);
    printf("2D res = %d\n", res2D);
    checkCudaErrors(cudaMemcpy(hD2Res, dD2Res, _FFT_TESTMATRIXX * _FFT_TESTMATRIXY * sizeof(cuComplex), cudaMemcpyDeviceToHost));
    PrintTestArray2D(hD2Res);

#pragma endregion

#pragma region 3D

    PrintTestArray3D(hD3Source);
    checkCudaErrors(cudaMemcpy(dD3Source, hD3Source, _FFT_TESTMATRIXX * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ * sizeof(cuComplex), cudaMemcpyHostToDevice));
    cufftHandle plan3d;
    cufftPlan3d(&plan3d, _FFT_TESTMATRIXX, _FFT_TESTMATRIXY, _FFT_TESTMATRIXZ, CUFFT_C2C);
    cufftResult res3D = cufftExecC2C(plan3d, dD3Source, dD3Res, CUFFT_FORWARD);
    printf("3D res = %d\n", res3D);
    checkCudaErrors(cudaMemcpy(hD3Res, dD3Res, _FFT_TESTMATRIXX * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ * sizeof(cuComplex), cudaMemcpyDeviceToHost));
    PrintTestArray3D(hD3Res);

#pragma endregion

#pragma region 4D

    PrintTestArray4D(hD4Source);
    checkCudaErrors(cudaMemcpy(dD4Source, hD4Source, _FFT_TESTMATRIXX * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ * _FFT_TESTMATRIXW * sizeof(cuComplex), cudaMemcpyHostToDevice));
    cufftHandle plan4d1;
    int n[3] = { _FFT_TESTMATRIXY, _FFT_TESTMATRIXZ, _FFT_TESTMATRIXW };
    int inembed[3] = { _FFT_TESTMATRIXY, _FFT_TESTMATRIXZ, _FFT_TESTMATRIXW };
    int dist = _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ * _FFT_TESTMATRIXW;
    cufftPlanMany(&plan4d1, 3, n,
        inembed, 1, dist,
        inembed, 1, dist,
        CUFFT_C2C, _FFT_TESTMATRIXX);
    
    cufftResult res4D1 = cufftExecC2C(plan4d1, dD4Source, dD4Res, CUFFT_FORWARD);
    printf("4D res 1 = %d\n", res4D1);

    checkCudaErrors(cudaMemcpy(hD4Res, dD4Res, _FFT_TESTMATRIXX * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ * _FFT_TESTMATRIXW * sizeof(cuComplex), cudaMemcpyDeviceToHost));
    PrintTestArray4D(hD4Res);

    cufftHandle plan4d2;
    int n2[1] = { _FFT_TESTMATRIXX };
    //note that if it was null, it will ignore the stride
    cufftPlanMany(&plan4d2, 1, n2,
        n2, dist, 1,
        n2, dist, 1,
        CUFFT_C2C, dist);

    //in out can be the same
    cufftResult res4D2 = cufftExecC2C(plan4d2, dD4Res, dD4Res, CUFFT_FORWARD);
    printf("4D res 2 = %d\n", res4D2);

    checkCudaErrors(cudaMemcpy(hD4Res, dD4Res, _FFT_TESTMATRIXX * _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ * _FFT_TESTMATRIXW * sizeof(cuComplex), cudaMemcpyDeviceToHost));
    PrintTestArray4D(hD4Res);

#pragma endregion

    //we have to do yzw 3D first
    //cufftHandle plan;
    //cufftPlan3d(&plan, _FFT_TESTMATRIXX, _FFT_TESTMATRIXY, _FFT_TESTMATRIXZ, CUFFT_C2C);
    //int n[3] = { _FFT_TESTMATRIXW, _FFT_TESTMATRIXZ, _FFT_TESTMATRIXY };
    //int inembed[3] = { _FFT_TESTMATRIXY, _FFT_TESTMATRIXZ, _FFT_TESTMATRIXW };
    //int dist = _FFT_TESTMATRIXY * _FFT_TESTMATRIXZ * _FFT_TESTMATRIXW;
    //cufftPlanMany(&plan, 3, n, 
    //    inembed, 1, dist,
    //    inembed, 1, dist,
    //    CUFFT_C2C, _FFT_TESTMATRIXX);

    


    

    //=======================
    checkCudaErrors(cudaFree(dD1Res));
    checkCudaErrors(cudaFree(dD1Source));
    checkCudaErrors(cudaFree(dD2Res));
    checkCudaErrors(cudaFree(dD2Source));
    checkCudaErrors(cudaFree(dD3Res));
    checkCudaErrors(cudaFree(dD3Source));
    checkCudaErrors(cudaFree(dD4Res));
    checkCudaErrors(cudaFree(dD4Source));
    free(hD1Res);
    free(hD1Source);
    free(hD2Res);
    free(hD2Source);
    free(hD3Res);
    free(hD3Source);
    free(hD4Res);
    free(hD4Source);

    return 0;
}
