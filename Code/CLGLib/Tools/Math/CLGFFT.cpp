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

#define _tfftMX 15
#define _tfftMY 16
#define _tfftMZ 17
#define _tfftMT 18

__BEGIN_NAMESPACE

UBOOL CCLGFFTHelper::FFT3DWithXYZW(CLGComplex* copied, TArray<INT> dims, UBOOL bForward)
{
    cufftHandle plan4d1;
    INT n[3] = { dims[1], dims[2], dims[3] };
    INT inembed[3] = { dims[1], dims[2], dims[3] };
    const INT dist = dims[1] * dims[2] * dims[3];

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
#endif
    return TRUE;
}

UBOOL CCLGFFTHelper::FFT4DWithXYZW(CLGComplex* copied, TArray<INT> dims, UBOOL bForward)
{
    if (!FFT3DWithXYZW(copied, dims, bForward))
    {
        return FALSE;
    }

    cufftHandle plan4d2;
    INT n2[1] = { dims[0] };
    const INT dist = dims[1] * dims[2] * dims[3];

#if _CLG_DOUBLEFLOAT

    //note that if it was null, it will ignore the stride
    const cufftResult planRes = cufftPlanMany(&plan4d2, 1, n2,
        n2, dist, 1,
        n2, dist, 1,
        CUFFT_Z2Z, dist);

    if (CUFFT_SUCCESS != planRes)
    {
        appCrucial(_T("cufftPlanMany failed! %d\n"), planRes);
        return FALSE;
    }

    //in out can be the same
    const cufftResult res4D2 = cufftExecZ2Z(plan4d2, copied, copied, CUFFT_FORWARD);
    if (CUFFT_SUCCESS != res4D2)
    {
        appCrucial(_T("cufftResult failed! %d\n"), res4D2);
        return FALSE;
    }

#else
    //note that if it was null, it will ignore the stride
    const cufftResult planRes = cufftPlanMany(&plan4d2, 1, n2,
        n2, dist, 1,
        n2, dist, 1,
        CUFFT_C2C, dist);

    if (CUFFT_SUCCESS != planRes)
    {
        appCrucial(_T("cufftPlanMany failed! %d\n"), planRes);
        return FALSE;
    }

    //in out can be the same
    const cufftResult res4D2 = cufftExecC2C(plan4d2, copied, copied, CUFFT_FORWARD);
    if (CUFFT_SUCCESS != res4D2)
    {
        appCrucial(_T("cufftResult failed! %d\n"), res4D2);
        return FALSE;
    }
#endif

    return TRUE;
}

UBOOL CCLGFFTHelper::FFT3D(CLGComplex* copied, UBOOL bForward)
{
    TArray<INT> dims;
    dims.AddItem(static_cast<INT>(_HC_Lx));
    dims.AddItem(static_cast<INT>(_HC_Ly));
    dims.AddItem(static_cast<INT>(_HC_Lz));
    dims.AddItem(static_cast<INT>(_HC_Lt));
    return FFT3DWithXYZW(copied, dims, bForward);
}

UBOOL CCLGFFTHelper::FFT4D(CLGComplex* copied, UBOOL bForward)
{
    TArray<INT> dims;
    dims.AddItem(static_cast<INT>(_HC_Lx));
    dims.AddItem(static_cast<INT>(_HC_Ly));
    dims.AddItem(static_cast<INT>(_HC_Lz));
    dims.AddItem(static_cast<INT>(_HC_Lt));
    return FFT4DWithXYZW(copied, dims, bForward);
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


void CCLGFFTHelper::TestFFT()
{
    appSetLogDate(FALSE);

    CLGComplex* dD4Res;
    CLGComplex* hD4Res = (CLGComplex*)malloc(_tfftMX * _tfftMY * _tfftMZ * _tfftMT * sizeof(CLGComplex));
    CLGComplex* hD4Source = (CLGComplex*)malloc(_tfftMX * _tfftMY * _tfftMZ * _tfftMT * sizeof(CLGComplex));

    checkCudaErrors(cudaMalloc((void**)& dD4Res, _tfftMX * _tfftMY * _tfftMZ * _tfftMT * sizeof(CLGComplex)));

    GenerateTestArray(hD4Source, _tfftMX * _tfftMY * _tfftMZ * _tfftMT);
    appGeneral(_T("(* Copy those to Mathematica to test *)\n"));
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

    for (INT i = 0; i < _tfftMX; ++i)
    {
        appGeneral(_T("InverseFourier[sour[[%d]]] Sqrt[%d] - (res1[[%d]]) //Abs//Max\n"),
            i + 1,
            _tfftMY * _tfftMZ * _tfftMT,
            i + 1);
    }

    appGeneral(_T("InverseFourier[sour] Sqrt[%d] - res2 //Abs//Max\n"), _tfftMX * _tfftMY * _tfftMZ * _tfftMT);

    appSetLogDate(TRUE);
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
