//=============================================================================
// FILENAME : CLGFFT.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [09/16/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include <cufft.h>

#define _tfftMX 4
#define _tfftMY 5
#define _tfftMZ 6
#define _tfftMT 7

__BEGIN_NAMESPACE

#pragma region kernel

__global__ void _CLG_LAUNCH_BOUND
_kernelScale(CLGComplex* res, Real fScale)
{
    intokernal;
    res[uiSiteIndex] = cuCmulf_cr(res[uiSiteIndex], fScale);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCopyElementSU3ToC(deviceSU3* su3, CLGComplex* res, BYTE idx)
{
    intokernal;
    res[uiSiteIndex] = su3[uiSiteIndex].m_me[idx];
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCopyElementCToSU3(deviceSU3* su3, CLGComplex* res, BYTE idx)
{
    intokernal;
    su3[uiSiteIndex].m_me[idx] = res[uiSiteIndex];
}

#pragma endregion

UBOOL CCLGFFTHelper::FFT3DWithXYZ(CLGComplex* copied, TArray<INT> dims, UBOOL bForward)
{
    cufftHandle plan3d;
    INT n[3] = { dims[0], dims[1], dims[2] };

#if _CLG_DOUBLEFLOAT

    const cufftResult planRes = cufftPlanMany(&plan3d, 3, n,
        n, 1, 1,
        n, 1, 1,
        CUFFT_Z2Z, 1);

    if (CUFFT_SUCCESS != planRes)
    {
        appCrucial(_T("cufftPlanMany failed! %d\n"), planRes);
        return FALSE;
    }

    const cufftResult res3D = cufftExecZ2Z(plan3d, copied, copied, bForward ? CUFFT_FORWARD : CUFFT_INVERSE);
    if (CUFFT_SUCCESS != res3D)
    {
        appCrucial(_T("cufftResult failed! %d\n"), res3D);
        return FALSE;
    }

#else

    const cufftResult planRes = cufftPlanMany(&plan3d, 3, n,
        n, 1, 1,
        n, 1, 1,
        CUFFT_C2C, 1);

    if (CUFFT_SUCCESS != planRes)
    {
        appCrucial(_T("cufftPlanMany failed! %d\n"), planRes);
        return FALSE;
    }

    const cufftResult res3D = cufftExecC2C(plan3d, copied, copied, bForward ? CUFFT_FORWARD : CUFFT_INVERSE);
    if (CUFFT_SUCCESS != res3D)
    {
        appCrucial(_T("cufftResult failed! %d\n"), res3D);
        return FALSE;
    }

#endif

    return TRUE;
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
UBOOL CCLGFFTHelper::FFT3DWithXYZW(CLGComplex* copied, TArray<INT> dims, UBOOL bForward)
{
    cufftHandle plan4d1;
    INT n[3] = { dims[0], dims[1], dims[2] };
    INT inembed[3] = { dims[0], dims[1], dims[2] };
    const INT stride = dims[3];

#if _CLG_DOUBLEFLOAT

    const cufftResult planRes = cufftPlanMany(&plan4d1, 3, n,
        inembed, stride, 1,
        inembed, stride, 1,
        CUFFT_Z2Z, stride);

    if (CUFFT_SUCCESS != planRes)
    {
        appCrucial(_T("cufftPlanMany failed! %d\n"), planRes);
        return FALSE;
    }

    const cufftResult res4D1 = cufftExecZ2Z(plan4d1, copied, copied, bForward ? CUFFT_FORWARD : CUFFT_INVERSE);
    if (CUFFT_SUCCESS != res4D1)
    {
        appCrucial(_T("cufftResult failed! %d\n"), res4D1);
        return FALSE;
    }

#else

    const cufftResult planRes = cufftPlanMany(&plan4d1, 3, n,
        inembed, stride, 1,
        inembed, stride, 1,
        CUFFT_C2C, stride);

    if (CUFFT_SUCCESS != planRes)
    {
        appCrucial(_T("cufftPlanMany failed! %d\n"), planRes);
        return FALSE;
    }

    const cufftResult res4D1 = cufftExecC2C(plan4d1, copied, copied, bForward ? CUFFT_FORWARD : CUFFT_INVERSE);
    if (CUFFT_SUCCESS != res4D1)
    {
        appCrucial(_T("cufftResult failed! %d\n"), res4D1);
        return FALSE;
    }

#endif
    return TRUE;
}

UBOOL CCLGFFTHelper::FFT4DWithXYZW(CLGComplex* copied, TArray<INT> dims, UBOOL bForward)
{
    cufftHandle plan4d1;
    INT n[3] = { dims[1], dims[2], dims[3] };
    INT inembed[3] = { dims[1], dims[2], dims[3] };
    const INT dist = dims[1] * dims[2] * dims[3];
    cufftHandle plan4d2;
    INT n2[1] = { dims[0] };

#if _CLG_DOUBLEFLOAT

    const cufftResult planRes = cufftPlanMany(&plan4d1, 3, n,
        inembed, 1, dist,
        inembed, 1, dist,
        CUFFT_Z2Z, dims[0]);

    if (CUFFT_SUCCESS != planRes)
    {
        appCrucial(_T("cufftPlanMany failed! %d\n"), planRes);
        return FALSE;
    }


    const cufftResult res4D1 = cufftExecZ2Z(plan4d1, copied, copied, bForward ? CUFFT_FORWARD : CUFFT_INVERSE);
    if (CUFFT_SUCCESS != res4D1)
    {
        appCrucial(_T("cufftResult failed! %d\n"), res4D1);
        return FALSE;
    }

    //note that if it was null, it will ignore the stride
    const cufftResult planRes1 = cufftPlanMany(&plan4d2, 1, n2,
        n2, dist, 1,
        n2, dist, 1,
        CUFFT_Z2Z, dist);

    if (CUFFT_SUCCESS != planRes1)
    {
        appCrucial(_T("cufftPlanMany failed! %d\n"), planRes);
        return FALSE;
    }

    //in out can be the same
    const cufftResult res4D2 = cufftExecZ2Z(plan4d2, copied, copied, bForward ? CUFFT_FORWARD : CUFFT_INVERSE);
    if (CUFFT_SUCCESS != res4D2)
    {
        appCrucial(_T("cufftResult failed! %d\n"), res4D2);
        return FALSE;
    }

#else

    const cufftResult planRes = cufftPlanMany(&plan4d1, 3, n,
        inembed, 1, dist,
        inembed, 1, dist,
        CUFFT_C2C, dims[0]);

    if (CUFFT_SUCCESS != planRes)
    {
        appCrucial(_T("cufftPlanMany failed! %d\n"), planRes);
        return FALSE;
    }


    const cufftResult res4D1 = cufftExecC2C(plan4d1, copied, copied, bForward ? CUFFT_FORWARD : CUFFT_INVERSE);
    if (CUFFT_SUCCESS != res4D1)
    {
        appCrucial(_T("cufftResult failed! %d\n"), res4D1);
        return FALSE;
    }

    //note that if it was null, it will ignore the stride
    const cufftResult planRes2 = cufftPlanMany(&plan4d2, 1, n2,
        n2, dist, 1,
        n2, dist, 1,
        CUFFT_C2C, dist);

    if (CUFFT_SUCCESS != planRes2)
    {
        appCrucial(_T("cufftPlanMany failed! %d\n"), planRes2);
        return FALSE;
    }

    //in out can be the same
    const cufftResult res4D2 = cufftExecC2C(plan4d2, copied, copied, bForward ? CUFFT_FORWARD : CUFFT_INVERSE);
    if (CUFFT_SUCCESS != res4D2)
    {
        appCrucial(_T("cufftResult failed! %d\n"), res4D2);
        return FALSE;
    }

#endif

    return TRUE;
}

UBOOL CCLGFFTHelper::FFT3D(CLGComplex* copied, UBOOL bForward, EFFT_Scale eScale)
{
    TArray<INT> dims;
    dims.AddItem(static_cast<INT>(_HC_Lx));
    dims.AddItem(static_cast<INT>(_HC_Ly));
    dims.AddItem(static_cast<INT>(_HC_Lz));
    dims.AddItem(static_cast<INT>(_HC_Lt));
    if (!FFT3DWithXYZW(copied, dims, bForward))
    {
        return FALSE;
    }

    if (bForward)
    {
        if (ES_1OverNForward == eScale)
        {
            preparethread;
            _kernelScale << <block, threads >> > (copied, F(1.0) / _HC_Volume_xyz);
        }
        else if (ES_1OverSqrtNBoth == eScale)
        {
            preparethread;
            _kernelScale << <block, threads >> > (copied, F(1.0) / _hostsqrt(static_cast<Real>(_HC_Volume_xyz)));
        }
    }
    else
    {
        if (ES_1OverNInverse == eScale)
        {
            preparethread;
            _kernelScale << <block, threads >> > (copied, F(1.0) / _HC_Volume_xyz);
        }
        else if (ES_1OverSqrtNBoth == eScale)
        {
            preparethread;
            _kernelScale << <block, threads >> > (copied, F(1.0) / _hostsqrt(static_cast<Real>(_HC_Volume_xyz)));
        }
    }
    return TRUE;
}

UBOOL CCLGFFTHelper::FFT4D(CLGComplex* copied, UBOOL bForward, EFFT_Scale eScale)
{
    TArray<INT> dims;
    dims.AddItem(static_cast<INT>(_HC_Lx));
    dims.AddItem(static_cast<INT>(_HC_Ly));
    dims.AddItem(static_cast<INT>(_HC_Lz));
    dims.AddItem(static_cast<INT>(_HC_Lt));
    if (!FFT4DWithXYZW(copied, dims, bForward))
    {
        return FALSE;
    }

    if (bForward)
    {
        if (ES_1OverNForward == eScale)
        {
            preparethread;
            _kernelScale << <block, threads >> > (copied, F(1.0) / _HC_Volume);
        }
        else if (ES_1OverSqrtNBoth == eScale)
        {
            preparethread;
            _kernelScale << <block, threads >> > (copied, F(1.0) / _hostsqrt(static_cast<Real>(_HC_Volume)));
        }
    }
    else
    {
        if (ES_1OverNInverse == eScale)
        {
            preparethread;
            _kernelScale << <block, threads >> > (copied, F(1.0) / _HC_Volume);
        }
        else if (ES_1OverSqrtNBoth == eScale)
        {
            preparethread;
            _kernelScale << <block, threads >> > (copied, F(1.0) / _hostsqrt(static_cast<Real>(_HC_Volume)));
        }
    }

    return TRUE;
}

void CCLGFFTHelper::CheckBuffer()
{
    if (NULL == m_pDeviceBuffer)
    {
        checkCudaErrors(cudaMalloc((void**)& m_pDeviceBuffer, _HC_Volume * sizeof(CLGComplex)));
    }
}

/**
 * To save memory, we do the FFT element by element
 */
UBOOL CCLGFFTHelper::FFT3DSU3(deviceSU3* res, UBOOL bForward, EFFT_Scale eScale)
{
    CheckBuffer();
    preparethread;
    for (BYTE i = 0; i < 9; ++i)
    {
        _kernelCopyElementSU3ToC << <block, threads>> > (res, m_pDeviceBuffer, i);
        if (!FFT3D(m_pDeviceBuffer, bForward, eScale))
        {
            return FALSE;
        }
        _kernelCopyElementCToSU3 << <block, threads >> > (res, m_pDeviceBuffer, i);
    }
    return TRUE;
}

UBOOL CCLGFFTHelper::FFT4DSU3(deviceSU3* res, UBOOL bForward, EFFT_Scale eScale)
{
    CheckBuffer();
    preparethread;
    for (BYTE i = 0; i < 9; ++i)
    {
        _kernelCopyElementSU3ToC << <block, threads >> > (res, m_pDeviceBuffer, i);
        if (!FFT4D(m_pDeviceBuffer, bForward, eScale))
        {
            return FALSE;
        }
        _kernelCopyElementCToSU3 << <block, threads >> > (res, m_pDeviceBuffer, i);
    }
    return TRUE;
}

void CCLGFFTHelper::GenerateTestArray(CLGComplex * hostArray, INT iSize)
{
    for (INT i = 0; i < iSize; ++i)
    {
        hostArray[i] = _make_cuComplex((rand() % 101 - 50) / F(50.0), (rand() % 101 - 50) / F(50.0));
    }
}

void CCLGFFTHelper::PrintTestArray4D(CLGComplex* hostArray)
{
    appGeneral(_T("{\n"));
    for (INT i = 0; i < _tfftMX; ++i)
    {
        appGeneral(_T("{"));
        for (INT j = 0; j < _tfftMY; ++j)
        {
            if (0 == j)
            {
                appGeneral(_T("{"));
            }
            else
            {
                appGeneral(_T(" {"));
            }
            for (INT k = 0; k < _tfftMZ; ++k)
            {
                appGeneral("{");
                for (INT l = 0; l < _tfftMT; ++l)
                {
                    const int iIndex = i * _tfftMY * _tfftMZ * _tfftMT
                        + j * _tfftMZ * _tfftMT
                        + k * _tfftMT
                        + l;
                    appGeneral(_T("%1.10f %s %1.10f I"),
                        hostArray[iIndex].x,
                        hostArray[iIndex].y < 0.0f ? _T("-") : _T("+"),
                        abs(hostArray[iIndex].y));
                    if (l == _tfftMT - 1)
                    {
                        appGeneral(_T("}"));
                    }
                    else
                    {
                        appGeneral(_T(", "));
                    }
                }

                if (k == _tfftMZ - 1)
                {
                    appGeneral(_T("}"));
                }
                else
                {
                    appGeneral(_T(", "));
                }
            }
            if (j == _tfftMY - 1)
            {
                appGeneral(_T("\n}"));
            }
            else
            {
                appGeneral(_T(",\n"));
            }
        }

        if (i == _tfftMX - 1)
        {
            appGeneral(_T("};\n"));
        }
        else
        {
            appGeneral(_T(",\n"));
        }
    }
}

void CCLGFFTHelper::PrintTestArray3D(CLGComplex* hostArray)
{
    appGeneral(_T("{\n"));
    for (INT i = 0; i < _tfftMX; ++i)
    {
        appGeneral(_T("{"));
        for (INT j = 0; j < _tfftMY; ++j)
        {
            if (0 == j)
            {
                appGeneral(_T("{"));
            }
            else
            {
                appGeneral(_T(" {"));
            }
            for (INT k = 0; k < _tfftMZ; ++k)
            {
                const int iIndex = i * _tfftMY * _tfftMZ
                                 + j * _tfftMZ
                                 + k;
                appGeneral(_T("%1.10f %s %1.10f I"),
                    hostArray[iIndex].x,
                    hostArray[iIndex].y < 0.0f ? _T("-") : _T("+"),
                    abs(hostArray[iIndex].y));

                if (k == _tfftMZ - 1)
                {
                    appGeneral(_T("}"));
                }
                else
                {
                    appGeneral(_T(", "));
                }
            }
            if (j == _tfftMY - 1)
            {
                appGeneral(_T("\n}"));
            }
            else
            {
                appGeneral(_T(",\n"));
            }
        }

        if (i == _tfftMX - 1)
        {
            appGeneral(_T("};\n"));
        }
        else
        {
            appGeneral(_T(",\n"));
        }
    }
}


void CCLGFFTHelper::TestFFT()
{
    appSetLogDate(FALSE);

    CLGComplex* dD3Res;
    CLGComplex* hD3Res = (CLGComplex*)malloc(_tfftMX * _tfftMY * _tfftMZ * sizeof(CLGComplex));
    CLGComplex* hD3Source = (CLGComplex*)malloc(_tfftMX * _tfftMY * _tfftMZ * sizeof(CLGComplex));
    checkCudaErrors(cudaMalloc((void**)& dD3Res, _tfftMX * _tfftMY * _tfftMZ * sizeof(CLGComplex)));
    appGeneral(_T("(* Copy those to Mathematica to test *)\n"));
    GenerateTestArray(hD3Source, _tfftMX * _tfftMY * _tfftMZ);
    appGeneral(_T("sour3d="));
    PrintTestArray3D(hD3Source);
    checkCudaErrors(cudaMemcpy(dD3Res, hD3Source, _tfftMX * _tfftMY * _tfftMZ * sizeof(CLGComplex), cudaMemcpyHostToDevice));
    TArray<INT> dims3d;
    dims3d.AddItem(_tfftMX);
    dims3d.AddItem(_tfftMY);
    dims3d.AddItem(_tfftMZ);
    CCLGFFTHelper::FFT3DWithXYZ(dD3Res, dims3d, TRUE);
    checkCudaErrors(cudaMemcpy(hD3Res, dD3Res, _tfftMX * _tfftMY * _tfftMZ * sizeof(CLGComplex), cudaMemcpyDeviceToHost));
    appGeneral(_T("res3d="));
    PrintTestArray3D(hD3Res);


    CLGComplex* dD4Res;
    CLGComplex* hD4Res = (CLGComplex*)malloc(_tfftMX * _tfftMY * _tfftMZ * _tfftMT * sizeof(CLGComplex));
    CLGComplex* hD4Source = (CLGComplex*)malloc(_tfftMX * _tfftMY * _tfftMZ * _tfftMT * sizeof(CLGComplex));

    checkCudaErrors(cudaMalloc((void**)& dD4Res, _tfftMX * _tfftMY * _tfftMZ * _tfftMT * sizeof(CLGComplex)));

    GenerateTestArray(hD4Source, _tfftMX * _tfftMY * _tfftMZ * _tfftMT);
    appGeneral(_T("sour="));
    PrintTestArray4D(hD4Source);
    checkCudaErrors(cudaMemcpy(dD4Res, hD4Source, _tfftMX * _tfftMY * _tfftMZ * _tfftMT * sizeof(CLGComplex), cudaMemcpyHostToDevice));
    TArray<INT> dims;
    dims.AddItem(_tfftMX);
    dims.AddItem(_tfftMY);
    dims.AddItem(_tfftMZ);
    dims.AddItem(_tfftMT);
    CCLGFFTHelper::FFT3DWithXYZW(dD4Res, dims, TRUE);
    checkCudaErrors(cudaMemcpy(hD4Res, dD4Res, _tfftMX * _tfftMY * _tfftMZ * _tfftMT * sizeof(CLGComplex), cudaMemcpyDeviceToHost));
    appGeneral(_T("res1="));
    PrintTestArray4D(hD4Res);

    checkCudaErrors(cudaMemcpy(dD4Res, hD4Source, _tfftMX * _tfftMY * _tfftMZ * _tfftMT * sizeof(CLGComplex), cudaMemcpyHostToDevice));
    CCLGFFTHelper::FFT4DWithXYZW(dD4Res, dims, TRUE);
    checkCudaErrors(cudaMemcpy(hD4Res, dD4Res, _tfftMX * _tfftMY * _tfftMZ * _tfftMT * sizeof(CLGComplex), cudaMemcpyDeviceToHost));
    appGeneral(_T("res2="));
    PrintTestArray4D(hD4Res);

    CCLGFFTHelper::FFT4DWithXYZW(dD4Res, dims, FALSE);
    checkCudaErrors(cudaMemcpy(hD4Res, dD4Res, _tfftMX * _tfftMY * _tfftMZ * _tfftMT * sizeof(CLGComplex), cudaMemcpyDeviceToHost));
    appGeneral(_T("res3="));
    PrintTestArray4D(hD4Res);

    appGeneral(_T("InverseFourier[sour3d] Sqrt[%d] - res3d //Abs//Max\n"), _tfftMX * _tfftMY * _tfftMZ);

    for (INT i = 0; i < _tfftMT; ++i)
    {
        appGeneral(_T("sour%d=Table[sour[[i]][[j]][[k]][[%d]], {i, 1, %d}, {j, 1, %d}, {k, 1, %d}];\n"),
            i + 1,
            i + 1,
            _tfftMX,
            _tfftMY,
            _tfftMZ);
        appGeneral(_T("subres%d=Table[res1[[i]][[j]][[k]][[%d]], {i, 1, %d}, {j, 1, %d}, {k, 1, %d}];\n"),
            i + 1,
            i + 1,
            _tfftMX,
            _tfftMY,
            _tfftMZ);
    }

    for (INT i = 0; i < _tfftMX; ++i)
    {
        appGeneral(_T("InverseFourier[sour%d] Sqrt[%d] - subres%d //Abs//Max\n"),
            i + 1,
            _tfftMX * _tfftMY * _tfftMZ,
            i + 1);
    }

    appGeneral(_T("InverseFourier[sour] Sqrt[%d] - res2 //Abs//Max\n"), _tfftMX * _tfftMY * _tfftMZ * _tfftMT);
    appGeneral(_T("sour %d - res3 //Abs//Max\n"), _tfftMX * _tfftMY * _tfftMZ * _tfftMT);

    appSetLogDate(TRUE);

    free(hD3Res);
    free(hD3Source);
    free(hD4Res);
    free(hD4Source);
    checkCudaErrors(cudaFree(dD3Res));
    checkCudaErrors(cudaFree(dD4Res));
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
