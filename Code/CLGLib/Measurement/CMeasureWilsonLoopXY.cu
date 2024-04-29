//=============================================================================
// FILENAME : CMeasureWilsonLoopXY.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [09/24/2021 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureWilsonLoopXY)

#pragma region kernles 

/**
 * block.x * thread.x = Lx * Ly
 * block.y * thread.y = Lz
 * block.z = thread.z = 1
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelWilsonLoopCalculatePXY(
    const deviceSU3* __restrict__ pDeviceData, deviceSU3* res, 
    const BYTE byFieldId, const INT tstart)
{
    const UINT uiXYZ = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_Lz + (threadIdx.y + blockIdx.y * blockDim.y);
    const UINT uiXYZ_Lt = uiXYZ * _DC_Lt;
    SSmallInt4 startSite = __deviceSiteIndexToInt4(uiXYZ_Lt + tstart);
    UINT uiSaveSiteIndex = 0;
    for (UINT t = 0; t < _DC_Lt; ++t)
    {
        if (0 == t)
        {
            const UINT uiSiteIndex = uiXYZ_Lt + tstart;
            uiSaveSiteIndex = uiXYZ_Lt;
            const UINT uiLinkIdx = _deviceGetLinkIndex(uiSiteIndex, 3);
            res[uiSaveSiteIndex] = pDeviceData[uiLinkIdx];
        }
        else
        {
            startSite.w = startSite.w + 1;
            const UINT uiSiteIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(startSite)].m_uiSiteIndex;
            const UINT uiLinkIdx = _deviceGetLinkIndex(uiSiteIndex, 3);
            res[uiSaveSiteIndex + 1] = res[uiSaveSiteIndex];
            res[uiSaveSiteIndex + 1].Mul(pDeviceData[uiLinkIdx]);

            startSite = __deviceSiteIndexToInt4(uiSiteIndex);
            uiSaveSiteIndex = uiSaveSiteIndex + 1;
        }
    }
}

/**
 * Block.x = Lz
 * Block.y = Lt
 *
 * For example, when linkLength is 3, links is [1, 1, 2] (x, x, y)
 * max product is 3, then we calculate:
 * W1=n-> n+(2,1)*i at t=0
 * W2=n-> n+(2,1)*i at t=t
 * for i from 1 to product
 * Then, tr(pP[n,t]*W2*pP[n+(2,1),t]^+*W1^+) is the loop for (2,1)*i
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelWilsonLoopsXY(
    const BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pP,
    const SSmallInt4 shift, 
    const SSmallInt4 center, 
    const INT product, const INT linkLength, const UINT shiftLength,
    const INT link1, const INT link2, const INT link3,
    const UINT maxLength, const BYTE tstart, UINT* counter, CLGComplex* correlator)
{
    INT iZ = threadIdx.x + blockIdx.x * blockDim.x;
    INT iT = threadIdx.y + blockIdx.y * blockDim.y + 1;
    if (iT > (_DC_Lt - 1))
    {
        return;
    }
    deviceSU3 w1 = deviceSU3::makeSU3Id();
    deviceSU3 w2 = deviceSU3::makeSU3Id();

    SSmallInt4 sSite4(center.x, center.y, iZ, iT);
    UINT uiSiteIndex = _deviceGetSiteIndex(sSite4);
    //with t - 1
    const INT uiStartPolyaIndex = uiSiteIndex - 1;
    SSmallInt4 point1(
        sSite4.x,
        sSite4.y,
        sSite4.z,
        tstart);
    SSmallInt4 point2 = sSite4;
    for (BYTE i = 0; i < tstart; ++i)
    {
        point2.w = point2.w + 1;
        point2 = __deviceSiteIndexToInt4(
            __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(point2)].m_uiSiteIndex);
    }

    INT links[3] = { link1, link2, link3};
    for (INT prod = 1; prod <= product; ++prod)
    {
        w1.Mul(_deviceLink(pDeviceData, point1, linkLength, byFieldId, links));
        w2.Mul(_deviceLink(pDeviceData, point2, linkLength, byFieldId, links));

        //======= shift the coordinate =====
        point1.Add(shift);
        point1 = __deviceSiteIndexToInt4(
            __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(point1)].m_uiSiteIndex);
        point2.Add(shift);
        point2 = __deviceSiteIndexToInt4(
            __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(point2)].m_uiSiteIndex);

        //(n0,0) to (n0,t)
        deviceSU3 loop = pP[uiStartPolyaIndex];
        //(n0,t) to (n',t)
        loop.Mul(w2);
        //(n',t) to (n',0)
        SSmallInt4 point3 = point2;
        point3.w = sSite4.w - 1;
        loop.MulDagger(pP[_deviceGetSiteIndex(point3)]);
        //(n',0) to (n0,0)
        loop.MulDagger(w1);
        CLGComplex res = loop.Tr();
                
        //======= record result ============
        INT c = prod * prod * shiftLength;
        if (c < static_cast<INT>(maxLength))
        {
            if (1 == iT)
            {
                atomicAdd(&counter[c], 1);
            }
            c = c * (_DC_Lt - 1) + static_cast<INT>(iT - 1);
            atomicAdd(&correlator[c].x, res.x);
            atomicAdd(&correlator[c].y, res.y);
        }
        else
        {
            break;
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelInitialWilsonCorrelatorOfSiteXY(UINT* counter, CLGComplex* correlator)
{
    UINT idx = threadIdx.x * (_DC_Lt - 1) + blockIdx.x;
    if (0 == blockIdx.x)
    {
        counter[threadIdx.x] = 0;
    }
    correlator[idx] = _make_cuComplex(F(0.0), F(0.0));
}

//block.x is Lt - 1
__global__ void _CLG_LAUNCH_BOUND
_kernelAverageWilsonLoopXY(UINT* counter, CLGComplex* correlator)
{
    const UINT uiIdx = threadIdx.x * (_DC_Lt - 1) + blockIdx.x;
    if (counter[threadIdx.x] > 0)
    {
        const Real fcounter = static_cast<Real>(counter[threadIdx.x]);
        correlator[uiIdx].x = correlator[uiIdx].x / fcounter;
        correlator[uiIdx].y = correlator[uiIdx].y / fcounter;
    }
}


#pragma endregion

CMeasureWilsonLoopXY::~CMeasureWilsonLoopXY()
{
    if (NULL != m_pHostCorrelator)
    {
        free(m_pHostCorrelator);
    }

    if (NULL != m_pHostCorrelatorCounter)
    {
        free(m_pHostCorrelatorCounter);
    }

    if (NULL != m_pTmpLoop)
    {
        checkCudaErrors(cudaFree(m_pTmpLoop));
    }

    if (NULL != m_pCorrelator)
    {
        checkCudaErrors(cudaFree(m_pCorrelator));
    }

    if (NULL != m_pCorrelatorCounter)
    {
        checkCudaErrors(cudaFree(m_pCorrelatorCounter));
    }
}

void CMeasureWilsonLoopXY::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    checkCudaErrors(cudaMalloc((void**)&m_pTmpLoop, sizeof(deviceSU3) * _HC_Lx * _HC_Ly * _HC_Lz * _HC_Lt));

    m_uiMaxLengthSq = ((_HC_Lx + 1) / 2) * ((_HC_Lx + 1) / 2)
        + ((_HC_Ly + 1) / 2) * ((_HC_Ly + 1) / 2)
        + 1;

    if (m_uiMaxLengthSq > _HC_ThreadConstraint)
    {
        m_uiMaxLengthSq = _HC_ThreadConstraint;
    }
    if (m_uiMaxLengthSq > _HC_ThreadConstraintX)
    {
        m_uiMaxLengthSq = _HC_ThreadConstraintX;
    }

    checkCudaErrors(cudaMalloc((void**)&m_pCorrelator, sizeof(CLGComplex) * m_uiMaxLengthSq * (_HC_Lt - 1)));
    checkCudaErrors(cudaMalloc((void**)&m_pCorrelatorCounter, sizeof(UINT) * m_uiMaxLengthSq));

    m_pHostCorrelator = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiMaxLengthSq * (_HC_Lt - 1));
    m_pHostCorrelatorCounter = (UINT*)malloc(sizeof(UINT) * m_uiMaxLengthSq);

    Reset();
}

//=================================
//(1) We calculate P(n,t) = prod P(n,t=0,1,...)
//(2) For each (n1,t=0), (n2 = n1 + R, t=1...Lt), we calculate the trace of the loop and store it in (R, t)
//This is W(R,t)
//(3) For each (n1,t from 2), we calculate log(W(R,t) / W(R,t+1)), store it in V(R)=(R, 0)
//(4) Then we calculate V(|R|) = average(V(R))
void CMeasureWilsonLoopXY::OnConfigurationAcceptedSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple)
{
    if (NULL == pAcceptGauge || EFT_GaugeSU3 != pAcceptGauge->GetFieldType())
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pAcceptGauge);

    const UINT tm1 = _HC_Lt - 1;
    dim3 block2(tm1, 1, 1);
    dim3 threads2(m_uiMaxLengthSq, 1, 1);

    _kernelInitialWilsonCorrelatorOfSiteXY << <block2, threads2 >> > (m_pCorrelatorCounter, m_pCorrelator);

    SSmallInt4 sOffsets[8] = 
    {
        SSmallInt4(1, 1, 0, 0),
        SSmallInt4(1, 1, 0, 0),
        SSmallInt4(-1, 1, 0, 0),
        SSmallInt4(-1, 1, 0, 0),
        SSmallInt4(1, -1, 0, 0),
        SSmallInt4(1, -1, 0, 0),
        SSmallInt4(-1, -1, 0, 0),
        SSmallInt4(-1, -1, 0, 0)
    };

    SBYTE sCenterX = static_cast<SBYTE>(_HC_Centerx);
    SBYTE sCenterY = static_cast<SBYTE>(_HC_Centery);
    SSmallInt4 sCenter[8] =
    {
        SSmallInt4(sCenterX, sCenterY, 0, 0),
        SSmallInt4(sCenterX, sCenterY, 0, 0),
        SSmallInt4(sCenterX - 1, sCenterY, 0, 0),
        SSmallInt4(sCenterX - 1, sCenterY, 0, 0),
        SSmallInt4(sCenterX, sCenterY - 1, 0, 0),
        SSmallInt4(sCenterX, sCenterY - 1, 0, 0),
        SSmallInt4(sCenterX - 1, sCenterY - 1, 0, 0),
        SSmallInt4(sCenterX - 1, sCenterY - 1, 0, 0)
    };

    INT products[8] =
    {
        _HC_Lxi / 2, _HC_Lxi / 2, _HC_Lxi / 2, _HC_Lxi / 2,
        _HC_Lxi / 2, _HC_Lxi / 2, _HC_Lxi / 2, _HC_Lxi / 2
    };

    INT iPathLengths[8] = 
    {
        2, 2, 2, 2,
        2, 2, 2, 2
    };

    INT iPaths[8][3] = 
    {
        {1, 2, 0},
        {2, 1, 0},
        {-1, 2, 0},
        {2, -1, 0},
        {1, -2, 0},
        {-2, 1, 0},
        {-1, -2, 0},
        {-2, -1, 0}
    };

    const dim3 block0(_HC_DecompY, _HC_DecompZ, 1);
    const dim3 threads0(_HC_DecompLy, _HC_DecompLz, 1);
    for (INT t = 0; t < _HC_Lti; ++t)
    {
        //_HC_DecompX * _HC_DecompLx =  Lx * Ly
        //_HC_DecompY * _HC_DecompLy = Lz
        dim3 block1(_HC_DecompX, _HC_DecompY, 1);
        dim3 threads1(_HC_DecompLx, _HC_DecompLy, 1);
        _kernelWilsonLoopCalculatePXY << <block1, threads1 >> > (
            pGaugeSU3->m_pDeviceData, m_pTmpLoop, pGaugeSU3->m_byFieldId, t);

        for (INT i = 0; i < 8; ++i)
        {
            UINT uiShiftLength = static_cast<UINT>(
                static_cast<INT>(sOffsets[i].x)* static_cast<INT>(sOffsets[i].x)
                + static_cast<INT>(sOffsets[i].y)* static_cast<INT>(sOffsets[i].y)
                );
            _kernelWilsonLoopsXY << <block0, threads0 >> > (
                pGaugeSU3->m_byFieldId,
                pGaugeSU3->m_pDeviceData,
                m_pTmpLoop,
                sOffsets[i],
                sCenter[i],
                products[i],
                iPathLengths[i],
                uiShiftLength,
                iPaths[i][0], iPaths[i][1], iPaths[i][2],
                m_uiMaxLengthSq,
                static_cast<BYTE>(t),
                m_pCorrelatorCounter,
                m_pCorrelator
                );
        }
    }

    _kernelAverageWilsonLoopXY << <block2, threads2 >> > (m_pCorrelatorCounter, m_pCorrelator);

    //extract res
    checkCudaErrors(cudaMemcpy(m_pHostCorrelatorCounter, m_pCorrelatorCounter, sizeof(UINT) * m_uiMaxLengthSq, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostCorrelator, m_pCorrelator, sizeof(CLGComplex) * m_uiMaxLengthSq * tm1, cudaMemcpyDeviceToHost));
    if (m_bShowResult)
    {
        appPushLogDate(FALSE);
    }
    TArray<TArray<CLGComplex>> thisConf;
    if (0 == m_uiConfigurationCount)
    {
        assert(0 == m_lstR.Num());
        assert(0 == m_lstC.Num());

        //we do not have L^2 < 1 Wilson loop
        for (UINT uiL = 1; uiL < m_uiMaxLengthSq; ++uiL)
        {
            if (m_pHostCorrelatorCounter[uiL] > 0)
            {
                m_lstR.AddItem(uiL);
                TArray<CLGComplex> thisR;
                for (UINT t = 0; t < tm1; ++t)
                {
                    thisR.AddItem(m_pHostCorrelator[uiL * tm1 + t]);
                }
                thisConf.AddItem(thisR);

                if (m_bShowResult)
                {
                    appDetailed(_T("C(%f): "), _hostsqrt(static_cast<Real>(uiL)));
                    for (UINT t = 0; t < tm1; ++t)
                    {
                        appDetailed(_T("t(%d)=%f + %f I"),
                            t + 1,
                            m_pHostCorrelator[uiL * tm1 + t].x,
                            m_pHostCorrelator[uiL * tm1 + t].y);
                    }
                    appDetailed(_T("\n"));
                }
            }
        }
    }
    else
    {
        for (INT i = 0; i < m_lstR.Num(); ++i)
        {
            assert(m_pHostCorrelatorCounter[m_lstR[i]] > 0);
            TArray<CLGComplex> thisR;
            for (UINT t = 0; t < tm1; ++t)
            {
                thisR.AddItem(m_pHostCorrelator[m_lstR[i] * tm1 + t]);
            }
            thisConf.AddItem(thisR);

            if (m_bShowResult)
            {
                appDetailed(_T("C(%f): "), _hostsqrt(static_cast<Real>(m_lstR[i])));
                for (UINT t = 0; t < tm1; ++t)
                {
                    appDetailed(_T("t(%d)=%f + %f I"),
                        t + 1,
                        m_pHostCorrelator[m_lstR[i] * tm1 + t].x,
                        m_pHostCorrelator[m_lstR[i] * tm1 + t].y);
                }
                appDetailed(_T("\n"));
            }
        }
    }
    m_lstC.AddItem(thisConf);

    if (m_bShowResult)
    {
        appPopLogDate();
    }
    ++m_uiConfigurationCount;
}

void CMeasureWilsonLoopXY::Report()
{
    assert(m_uiConfigurationCount > 0);
    assert(static_cast<UINT>(m_uiConfigurationCount)
        == static_cast<UINT>(m_lstC.Num()));
    assert(static_cast<UINT>(m_lstR.Num())
        == static_cast<UINT>(m_lstC[0].Num()));
    assert((_HC_Lt - 1)
        == static_cast<UINT>(m_lstC[0][0].Num()));

    appPushLogDate(FALSE);
    const UINT tm1 = _HC_Lt - 1;
    TArray<CLGComplex> tmpLoop;
    TArray<CLGComplex> tmpCorrelator;

    appGeneral(_T("\n\n==========================================================================\n"));
    appGeneral(_T("==================== V(r) (%d con)============================\n"), m_uiConfigurationCount);

    m_lstAverageC.RemoveAll();

    appGeneral(_T("\n ----------- Static Quark Correlator ------------- \n"));

    appGeneral(_T("lengthR = {"));
    for (INT i = 0; i < m_lstR.Num(); ++i)
    {
        appGeneral(_T("%2.12f%s "),
            _hostsqrt(static_cast<Real>(m_lstR[i])),
            (i == m_lstR.Num() - 1) ? _T("") : _T(","));
    }
    appGeneral(_T("}\n"));

    if (m_bShowResult)
    {
        appGeneral(_T("correlator = {\n"));
    }

    for (INT k = 0; k < m_lstR.Num(); ++k)
    {
        TArray<CLGComplex> finalAverage;
        appGeneral(_T("{\n"));
        for (UINT t = 0; t < tm1; ++t)
        {
            if (m_bShowResult)
            {
                appGeneral(_T("{"));
            }
            CLGComplex averageOfC_R = _make_cuComplex(F(0.0), F(0.0));
            for (UINT i = 0; i < m_uiConfigurationCount; ++i)
            {
                if (m_bShowResult)
                {
                    LogGeneralComplex(m_lstC[i][k][t]);
                }
                averageOfC_R = _cuCaddf(averageOfC_R, m_lstC[i][k][t]);
            }
            if (m_bShowResult)
            {
                appGeneral(_T("},\n"));
            }

            averageOfC_R.x = averageOfC_R.x / m_uiConfigurationCount;
            averageOfC_R.y = averageOfC_R.y / m_uiConfigurationCount;
            finalAverage.AddItem(averageOfC_R);
        }
        appGeneral(_T("},\n"));
        m_lstAverageC.AddItem(finalAverage);
    }
    if (m_bShowResult)
    {
        appGeneral(_T("}\n"));
    }

    appGeneral(_T("averagecorrelator = {\n"));
    for (INT k = 0; k < m_lstR.Num(); ++k)
    {
        appGeneral(_T("{\n"));
        for (UINT t = 0; t < tm1; ++t)
        {
            LogGeneralComplex(m_lstAverageC[k][t], t != tm1 - 1);
        }
        appGeneral(_T("}%s\n"), (k == m_lstR.Num() - 1) ? _T("") : _T(","));
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("\n==========================================================================\n"));
    appGeneral(_T("==========================================================================\n\n"));
    appPopLogDate();
}

void CMeasureWilsonLoopXY::Reset()
{
    CMeasure::Reset();
    m_lstR.RemoveAll();
    m_lstC.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================