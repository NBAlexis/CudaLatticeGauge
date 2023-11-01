//=============================================================================
// FILENAME : CMeasureWilsonLoop.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [09/26/2020 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureWilsonLoop)

#pragma region kernles 

/**
 * block.x * thread.x = Lx * Ly
 * block.y * thread.y = Lz
 * block.z = thread.z = 1
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelWilsonLoopCalculateP(
    const deviceSU3* __restrict__ pDeviceData, deviceSU3* res, 
    const BYTE byFieldId, const BYTE tstart)
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
 * Block.x = x * Ly + y, Thread.x = Lz
 * Block.y = Lt - 1, 
 *
 * For example, when linkLength is 3, links is [1, 1, 2] (x, x, y)
 * max product is 3, then we calculate:
 * W1=n-> n+(2,1)*i at t=0
 * W2=n-> n+(2,1)*i at t=t
 * for i from 1 to product
 * Then, tr(pP[n,t]*W2*pP[n+(2,1),t]^+*W1^+) is the loop for (2,1)*i
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelWilsonLoops(
    const BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pP,
    const SSmallInt4 shift, const INT product, const INT linkLength, const UINT shiftLength,
    const INT link1, const INT link2, const INT link3, const INT link4, const INT link5,
    const UINT maxLength, const BYTE tstart, UINT* counter, CLGComplex* correlator)
{
    intokernalInt4;
    if (0 == sSite4.w || sSite4.w > (_DC_Lt / 2))
    {
        return;
    }

    deviceSU3 w1 = deviceSU3::makeSU3Id();
    deviceSU3 w2 = deviceSU3::makeSU3Id();

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

    INT links[5] = { link1, link2, link3, link4, link5 };
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
            if (1 == sSite4.w)
            {
                atomicAdd(&counter[c], 1);
            }
            c = c * (_DC_Lt / 2) + static_cast<INT>(sSite4.w - 1);
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
_kernelInitialWilsonCorrelatorOfSite(UINT* counter, CLGComplex* correlator)
{
    UINT idx = threadIdx.x * (_DC_Lt / 2) + blockIdx.x;
    if (0 == blockIdx.x)
    {
        counter[threadIdx.x] = 0;
    }
    correlator[idx] = _make_cuComplex(F(0.0), F(0.0));
}

//block.x is Lt - 1
__global__ void _CLG_LAUNCH_BOUND
_kernelAverageWilsonLoop(UINT* counter, CLGComplex* correlator)
{
    const UINT uiIdx = threadIdx.x * (_DC_Lt / 2) + blockIdx.x;
    if (counter[threadIdx.x] > 0)
    {
        const Real fcounter = static_cast<Real>(counter[threadIdx.x]);
        correlator[uiIdx].x = correlator[uiIdx].x / fcounter;
        correlator[uiIdx].y = correlator[uiIdx].y / fcounter;
    }
}


#pragma endregion

CMeasureWilsonLoop::~CMeasureWilsonLoop()
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

void CMeasureWilsonLoop::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    checkCudaErrors(cudaMalloc((void**)&m_pTmpLoop, sizeof(deviceSU3) * _HC_Lx * _HC_Ly * _HC_Lz * _HC_Lt));

    m_uiMaxLengthSq = ((_HC_Lx + 1) / 2) * ((_HC_Lx + 1) / 2)
        + ((_HC_Ly + 1) / 2) * ((_HC_Ly + 1) / 2)
        + ((_HC_Lz + 1) / 2) * ((_HC_Lz + 1) / 2)
        + 1;

    if (m_uiMaxLengthSq > _HC_ThreadConstraint)
    {
        m_uiMaxLengthSq = _HC_ThreadConstraint;
    }
    if (m_uiMaxLengthSq > _HC_ThreadConstraintX)
    {
        m_uiMaxLengthSq = _HC_ThreadConstraintX;
    }

    checkCudaErrors(cudaMalloc((void**)&m_pCorrelator, sizeof(CLGComplex) * m_uiMaxLengthSq * (_HC_Lt / 2)));
    checkCudaErrors(cudaMalloc((void**)&m_pCorrelatorCounter, sizeof(UINT) * m_uiMaxLengthSq));

    m_pHostCorrelator = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiMaxLengthSq * (_HC_Lt / 2));
    m_pHostCorrelatorCounter = (UINT*)malloc(sizeof(UINT) * m_uiMaxLengthSq);

    Reset();

    INT iValue = 1;
    param.FetchValueINT(_T("FieldId"), iValue);
    m_byFieldId = static_cast<BYTE>(iValue);

    iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;
}

//=================================
//(1) We calculate P(n,t) = prod P(n,t=0,1,...)
//(2) For each (n1,t=0), (n2 = n1 + R, t=1...Lt), we calculate the trace of the loop and store it in (R, t)
//This is W(R,t)
//(3) For each (n1,t from 2), we calculate log(W(R,t) / W(R,t+1)), store it in V(R)=(R, 0)
//(4) Then we calculate V(|R|) = average(V(R))
void CMeasureWilsonLoop::OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple)
{
    if (NULL == pAcceptGauge || EFT_GaugeSU3 != pAcceptGauge->GetFieldType())
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pAcceptGauge);

    preparethread;
    const UINT halfT = _HC_Lt / 2;
    dim3 block2(halfT, 1, 1);
    dim3 threads2(m_uiMaxLengthSq, 1, 1);

    _kernelInitialWilsonCorrelatorOfSite << <block2, threads2 >> > (m_pCorrelatorCounter, m_pCorrelator);

    SSmallInt4 sOffsets[19] = 
    {
        SSmallInt4(1, 0, 0, 0),
        SSmallInt4(0, 1, 0, 0),
        SSmallInt4(0, 0, 1, 0),

        SSmallInt4(1, 1, 0, 0),
        SSmallInt4(1, 0, 1, 0),
        SSmallInt4(0, 1, 1, 0),

        SSmallInt4(2, 1, 0, 0),
        SSmallInt4(1, 2, 0, 0),
        SSmallInt4(0, 2, 1, 0),
        SSmallInt4(0, 1, 2, 0),
        SSmallInt4(2, 0, 1, 0),
        SSmallInt4(1, 0, 2, 0),

        SSmallInt4(1, 1, 1, 0),

        SSmallInt4(2, 1, 1, 0),
        SSmallInt4(1, 2, 1, 0),
        SSmallInt4(1, 1, 2, 0),

        SSmallInt4(2, 2, 1, 0),
        SSmallInt4(1, 2, 2, 0),
        SSmallInt4(2, 1, 2, 0)
    };

    INT products[19] =
    {
        _HC_Lxi / 2, _HC_Lyi / 2, _HC_Lzi / 2,
        _HC_Lxi / 2, _HC_Lyi / 2, _HC_Lzi / 2,

        _HC_Lxi / 4, _HC_Lyi / 4, _HC_Lyi / 4,
        _HC_Lzi / 4, _HC_Lxi / 4, _HC_Lzi / 4,

        _HC_Lxi / 2,

        _HC_Lxi / 4, _HC_Lyi / 4, _HC_Lzi / 4,
        _HC_Lxi / 4, _HC_Lyi / 4, _HC_Lzi / 4
    };

    INT iPathLengths[19] = 
    {
        1, 1, 1,
        2, 2, 2,
        3, 3, 3, 3, 3, 3,
        3,
        4, 4, 4,
        5, 5, 5
    };

    INT iPaths[19][5] = 
    {
        {1, 0, 0, 0, 0},
        {2, 0, 0, 0, 0},
        {3, 0, 0, 0, 0},

        {1, 2, 0, 0, 0},
        {1, 3, 0, 0, 0},
        {2, 3, 0, 0, 0},

        {1, 2, 1, 0, 0},
        {2, 1, 2, 0, 0},
        {2, 3, 2, 0, 0},
        {3, 2, 3, 0, 0},
        {1, 3, 1, 0, 0},
        {3, 1, 3, 0, 0},

        {1, 2, 3, 0, 0},

        {1, 2, 3, 1, 0},
        {2, 3, 1, 2, 0},
        {3, 1, 2, 3, 0},

        {1, 2, 3, 1, 2},
        {2, 3, 1, 2, 3},
        {3, 1, 2, 3, 1}
    };

    for (INT t = 0; t < _HC_Lti; ++t)
    {
        //_HC_DecompX * _HC_DecompLx =  Lx * Ly
        //_HC_DecompY * _HC_DecompLy = Lz
        dim3 block1(_HC_DecompX, _HC_DecompY, 1);
        dim3 threads1(_HC_DecompLx, _HC_DecompLy, 1);
        _kernelWilsonLoopCalculateP << <block1, threads1 >> > (
            pGaugeSU3->m_pDeviceData, m_pTmpLoop, pGaugeSU3->m_byFieldId, static_cast<BYTE>(t));

        for (INT i = 0; i < 19; ++i)
        {
            UINT uiShiftLength = static_cast<UINT>(
                static_cast<INT>(sOffsets[i].x)* static_cast<INT>(sOffsets[i].x)
                + static_cast<INT>(sOffsets[i].y)* static_cast<INT>(sOffsets[i].y)
                + static_cast<INT>(sOffsets[i].z)* static_cast<INT>(sOffsets[i].z)
                );
            _kernelWilsonLoops << <block, threads >> > (
                pGaugeSU3->m_byFieldId,
                pGaugeSU3->m_pDeviceData,
                m_pTmpLoop,
                sOffsets[i],
                products[i],
                iPathLengths[i],
                uiShiftLength,
                iPaths[i][0], iPaths[i][1], iPaths[i][2], iPaths[i][3], iPaths[i][4],
                m_uiMaxLengthSq,
                static_cast<BYTE>(t),
                m_pCorrelatorCounter,
                m_pCorrelator
                );
        }
    }

    _kernelAverageWilsonLoop << <block2, threads2 >> > (m_pCorrelatorCounter, m_pCorrelator);

    //extract res
    checkCudaErrors(cudaMemcpy(m_pHostCorrelatorCounter, m_pCorrelatorCounter, sizeof(UINT) * m_uiMaxLengthSq, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostCorrelator, m_pCorrelator, sizeof(CLGComplex) * m_uiMaxLengthSq * halfT, cudaMemcpyDeviceToHost));
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
                for (UINT t = 0; t < halfT; ++t)
                {
                    thisR.AddItem(m_pHostCorrelator[uiL * halfT + t]);
                }
                thisConf.AddItem(thisR);

                if (m_bShowResult)
                {
                    appDetailed(_T("C(%f): "), _hostsqrt(static_cast<Real>(uiL)));
                    for (UINT t = 0; t < halfT; ++t)
                    {
                        appDetailed(_T("t(%d)=%f + %f I"),
                            t + 1,
                            m_pHostCorrelator[uiL * halfT + t].x,
                            m_pHostCorrelator[uiL * halfT + t].y);
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
            for (UINT t = 0; t < halfT; ++t)
            {
                thisR.AddItem(m_pHostCorrelator[m_lstR[i] * halfT + t]);
            }
            thisConf.AddItem(thisR);

            if (m_bShowResult)
            {
                appDetailed(_T("C(%f): "), _hostsqrt(static_cast<Real>(m_lstR[i])));
                for (UINT t = 0; t < halfT; ++t)
                {
                    appDetailed(_T("t(%d)=%f + %f I"),
                        t + 1,
                        m_pHostCorrelator[m_lstR[i] * halfT + t].x,
                        m_pHostCorrelator[m_lstR[i] * halfT + t].y);
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

void CMeasureWilsonLoop::Report()
{
    assert(m_uiConfigurationCount > 0);
    assert(static_cast<UINT>(m_uiConfigurationCount)
        == static_cast<UINT>(m_lstC.Num()));
    assert(static_cast<UINT>(m_lstR.Num())
        == static_cast<UINT>(m_lstC[0].Num()));
    assert(_HC_Lt / 2
        == static_cast<UINT>(m_lstC[0][0].Num()));

    appPushLogDate(FALSE);
    const UINT halfT = _HC_Lt / 2;
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
        for (UINT t = 0; t < halfT; ++t)
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
        for (UINT t = 0; t < halfT; ++t)
        {
            LogGeneralComplex(m_lstAverageC[k][t], t != halfT - 1);
        }
        appGeneral(_T("}%s\n"), (k == m_lstR.Num() - 1) ? _T("") : _T(","));
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("\n==========================================================================\n"));
    appGeneral(_T("==========================================================================\n\n"));
    appPopLogDate();
}

void CMeasureWilsonLoop::Reset()
{
    CMeasure::Reset();
    m_lstR.RemoveAll();
    m_lstC.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================