//=============================================================================
// FILENAME : CGaugeSmearingHISQ.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [mm/dd/yy]
//  [11/14/2024 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "Data/Field/Gauge/CFieldGaugeU1Real.h"
#include "CGaugeSmearingASQTAD.h"
#include "CGaugeSmearingHISQ.h"

__BEGIN_NAMESPACE

#pragma region kernels

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelNaikLink(
    deviceGauge* pResult,
    const deviceGauge* __restrict__ pEffectiveGaugeL1,
    BYTE byGaugeFieldId)
{
    intokernalDirInt4;
    const SCHAR fwddir = dir + 1;
    SCHAR pathBuffer[3] = { fwddir, fwddir, fwddir };
    pResult[uiLinkIndex] = _deviceLinkT(pEffectiveGaugeL1, sSite4, 3, byGaugeFieldId, pathBuffer);
}

//__global__ void _CLG_LAUNCH_BOUND
//_kernelNaikLinkU1Real(
//    Real* pResult,
//    const Real* __restrict__ field,
//    BYTE byGaugeFieldId)
//{
//    intokernalDirInt4;
//    const SCHAR fwddir = dir + 1;
//    SCHAR pathBuffer[3] = { fwddir, fwddir, fwddir };
//    pResult[uiLinkIndex] = _deviceLinkT(field, sSite4, 3, byGaugeFieldId, pathBuffer);
//}

#pragma endregion

template<typename gaugetype, INT matrixN>
CGaugeSmearingHISQ<gaugetype, matrixN>::~CGaugeSmearingHISQ()
{
    //no matter whether it is returned, it will be destroyed
    //if (NULL != m_pOriginalGauge)
    //{
    //    m_pOriginalGauge->Return();
    //    m_pOriginalGauge = NULL;
    //}
    if (NULL != m_pDeviceRationalApproximation)
    {
        checkCudaErrors(__cudaFree(m_pDeviceRationalApproximation));
        m_pDeviceRationalApproximation = NULL;
    }
    if (NULL != m_pProjForcePtr)
    {
        checkCudaErrors(__cudaFree(m_pProjForcePtr));
        m_pProjForcePtr = NULL;
    }
    if (NULL != m_pDeviceDet)
    {
        checkCudaErrors(__cudaFree(m_pDeviceDet));
        m_pDeviceDet = NULL;
    }
    if (NULL != m_pDeviceDetPhase)
    {
        checkCudaErrors(__cudaFree(m_pDeviceDetPhase));
        m_pDeviceDetPhase = NULL;
    }
    if (NULL != m_pCaylayHamiltonCache)
    {
        checkCudaErrors(__cudaFree(m_pCaylayHamiltonCache));
        m_pCaylayHamiltonCache = NULL;
    }
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingHISQ<gaugetype, matrixN>::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    CGaugeSmearing::Initial(pOwner, params);

    //m_bProj = TRUE;
    //m_bProjDet = TRUE;
    m_bCalledWhenUpdate = TRUE;



    INT iVaule = 1;
    params.FetchValueINT(_T("Proj"), iVaule);
    m_bProj = (0 != iVaule);

    iVaule = 1;
    params.FetchValueINT(_T("ProjSUN"), iVaule);
    m_bProjDet = (0 != iVaule);

    iVaule = 1;
    params.FetchValueINT(_T("UseCaylayHamilton"), iVaule);
    m_bUseCaylayHamilton = (0 != iVaule);

    if (!UseCaylayHamilton())
    {
        Real fra[11] = { F(0.28017824602695207), F(0.019466137466992547), F(0.03561358805516765),
                  F(0.08215727655939538), F(0.21113622525403777), F(0.7946025292571568), F(0.00020653817365416855),
                  F(0.00302707751043529), F(0.020073267806506062), F(0.1251758627154607), F(1.0029328744648631) };
        checkCudaErrors(__cudaMalloc((void**)&m_pDeviceRationalApproximation, sizeof(Real) * 11));
        checkCudaErrors(cudaMemcpy(m_pDeviceRationalApproximation, fra, sizeof(Real) * 11, cudaMemcpyHostToDevice));
        m_uiRationalApproximationOrder = 5;

        //Real fra[9] = { F(0.0954113732807719), F(0.20212594413375087), F(0.30595220754840236),
        //            F(0.6336901399796908), F(2.269185448392556), F(0.02335218756717348),
        //            F(0.2862905885952218), F(1.4973064554346185), F(9.531358573058725) };
        //checkCudaErrors(__cudaMalloc((void**)&m_pDeviceRationalApproximation, sizeof(Real) * 9));
        //checkCudaErrors(cudaMemcpy(m_pDeviceRationalApproximation, fra, sizeof(Real) * 9, cudaMemcpyHostToDevice));
        //m_uiRationalApproximationOrder = 4;
    }

    Real fU0 = F(1.0);
    if (params.FetchValueReal(_T("U0"), fU0))
    {
        params.FetchValueReal(_T("OriginL1"), m_fOriginalL1);
        params.FetchValueReal(_T("Fat3L1"), m_fFat3L1);
        params.FetchValueReal(_T("Fat5L1"), m_fFat5L1);
        params.FetchValueReal(_T("Fat7L1"), m_fFat7L1);
        params.FetchValueReal(_T("LepageL1"), m_fLepageL1);

        Real invseu0 = F(1.0) / fU0;
        Real u2 = invseu0 * invseu0;
        Real u4 = u2 * u2;
        Real u6 = u2 * u4;
        m_fFat3L1 = m_fFat3L1 * u2;
        m_fFat5L1 = m_fFat5L1 * u4;
        m_fLepageL1 = m_fLepageL1 * u4;
        m_fFat7L1 = m_fFat7L1 * u6;
    }
    else
    {
        params.FetchValueReal(_T("OriginL1"), m_fOriginalL1);
        params.FetchValueReal(_T("Fat3L1"), m_fFat3L1);
        params.FetchValueReal(_T("Fat5L1"), m_fFat5L1);
        params.FetchValueReal(_T("Fat7L1"), m_fFat7L1);
        params.FetchValueReal(_T("LepageL1"), m_fLepageL1);
    }

    params.FetchValueReal(_T("OriginL2"), m_fOriginalL2);
    params.FetchValueReal(_T("Fat3L2"), m_fFat3L2);
    params.FetchValueReal(_T("Fat5L2"), m_fFat5L2);
    params.FetchValueReal(_T("Fat7L2"), m_fFat7L2);
    params.FetchValueReal(_T("LepageL2"), m_fLepageL2);

    params.FetchValueINT(_T("PhaseField"), m_iPhaseFieldId);

    CCString sCacheValue = _T("EHLC_Full");
    if (params.FetchStringValue(_T("Cache"), sCacheValue))
    {
        m_eCache = __STRING_TO_ENUM(EHISQLinkCache, sCacheValue);
    }
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingHISQ<gaugetype, matrixN>::GaugeSmearingFull(CFieldGauge* pGauge, const class CFieldGauge* pOrignalGauge, CFieldGauge* pStaple, UBOOL bProject)
{
    _RECORD(CGaugeSmearingHISQ::GaugeSmearingFull);
    appParanoiac(_T("CGaugeSmearingHISQ::GaugeSmearingFull\n"));

    if (NULL == m_pEffectiveGaugeL1)
    {
        m_pEffectiveGaugeL1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

        m_pP3_1_L1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP3_2_L1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP3_3_L1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

        m_pP5_1_1_L1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP5_1_2_L1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP5_2_1_L1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP5_2_2_L1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP5_3_1_L1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP5_3_2_L1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

        m_pP3_1_L2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP3_2_L2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP3_3_L2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

        m_pP5_1_1_L2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP5_1_2_L2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP5_2_1_L2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP5_2_2_L2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP5_3_1_L2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP5_3_2_L2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

        if (m_bProj)
        {
            m_pGaugeNotProjected = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

            if (UseCaylayHamilton())
            {
                m_pQ = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
                m_pQ2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
                m_pInverseSqrtQ = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
                checkCudaErrors(__cudaMalloc((void**)&m_pCaylayHamiltonCache, sizeof(DOUBLE) * _HC_Volume * _HC_Dir * 5));

                if (m_bProjDet)
                {
                    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceDetPhase, sizeof(DOUBLE) * _HC_Dir * _HC_Volume));
                }
            }
            else
            {
                TArray<gaugetype*> ptrs;
                for (UINT i = 0; i < m_uiRationalApproximationOrder + 1; ++i)
                {
                    m_pProjForce.AddItem(dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__)));
                    CFieldGaugeLink<gaugetype, matrixN>* pProjForce = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(m_pProjForce[i]);
                    ptrs.AddItem(pProjForce->m_pDeviceData);
                }
                checkCudaErrors(__cudaMalloc((void**)&m_pProjForcePtr, sizeof(gaugetype*) * (m_uiRationalApproximationOrder + 1)));
                checkCudaErrors(cudaMemcpy(m_pProjForcePtr, ptrs.GetData(), sizeof(gaugetype*) * (m_uiRationalApproximationOrder + 1), cudaMemcpyHostToDevice));

                if (m_bProjDet)
                {
                    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceDet, sizeof(CLGComplex) * _HC_Dir * _HC_Volume));
                }
            }
        }
    }

    //Level 1
    CGaugeSmearingASQTAD<gaugetype, matrixN>::Fat357Lepage(
        pGauge, //result
        NULL, //preserve
        m_pP3_1_L1,
        m_pP3_2_L1,
        m_pP3_3_L1,
        m_pP5_1_1_L1,
        m_pP5_1_2_L1,
        m_pP5_2_1_L1,
        m_pP5_2_2_L1,
        m_pP5_3_1_L1,
        m_pP5_3_2_L1,
        m_fOriginalL1,
        m_fFat3L1,
        m_fFat5L1,
        m_fFat7L1,
        m_fLepageL1,
        pOrignalGauge
    );

    if (m_bProj)
    {
        if (UseCaylayHamilton())
        {
            CGaugeSmearingASQTAD<gaugetype, matrixN>::ProjectCaylayHamilton(pGauge,
                m_pGaugeNotProjected,
                m_pQ,
                m_pQ2,
                m_pInverseSqrtQ,
                m_pCaylayHamiltonCache,
                m_bProjDet,
                m_pDeviceDetPhase);
        }
        else
        {
            CGaugeSmearingASQTAD<gaugetype, matrixN>::ProjectRationalApproximation(pGauge,
                m_pGaugeNotProjected,
                m_pProjForcePtr,
                m_pDeviceRationalApproximation,
                m_uiRationalApproximationOrder,
                m_bProjDet,
                m_pDeviceDet);
        }
    }

    //Level 2
    CGaugeSmearingASQTAD<gaugetype, matrixN>::Fat357Lepage(pGauge,
        m_pEffectiveGaugeL1, //preserve
        m_pP3_1_L2,
        m_pP3_2_L2,
        m_pP3_3_L2,
        m_pP5_1_1_L2,
        m_pP5_1_2_L2,
        m_pP5_2_1_L2,
        m_pP5_2_2_L2,
        m_pP5_3_1_L2,
        m_pP5_3_2_L2,
        m_fOriginalL2,
        m_fFat3L2,
        m_fFat5L2,
        m_fFat7L2,
        m_fLepageL2,
        NULL //if preserve field is provided, use preserve gauge for orignal gauge
    );

    if (NULL == m_pNaik)
    {
        m_pNaik = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_pEffectiveGaugeL1->m_byFieldId, _T(__FILE__), __LINE__));
    }
    CalculateNaikLink((gaugetype*)m_pNaik->GetData(), (const gaugetype*)m_pEffectiveGaugeL1->GetData(), m_pEffectiveGaugeL1->m_byFieldId);

    if (m_iPhaseFieldId > 0)
    {
        if (NULL == m_pNaikLinkPhase)
        {
            m_pNaikLinkPhase = dynamic_cast<CFieldGaugeU1Real*>(appGetLattice()->GetPooledFieldById(static_cast<BYTE>(m_iPhaseFieldId), _T(__FILE__), __LINE__));
            CalculateU1NaikLink(m_pNaikLinkPhase, dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(static_cast<BYTE>(m_iPhaseFieldId))));
        }
    }
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingHISQ<gaugetype, matrixN>::GaugeSmearingMedian(CFieldGauge* pGauge, const class CFieldGauge* pOrignalGauge, CFieldGauge* pStaple, UBOOL bProject)
{
    _RECORD(CGaugeSmearingHISQ::GaugeSmearingMedian);
    appParanoiac(_T("CGaugeSmearingHISQ::GaugeSmearingMedian\n"));

    if (NULL == m_pEffectiveGaugeL1)
    {
        m_pEffectiveGaugeL1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

        m_pP3_1_L1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP3_2_L1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP3_3_L1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

        m_pP3_1_L2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP3_2_L2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP3_3_L2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

        if (m_bProj)
        {
            m_pGaugeNotProjected = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

            if (UseCaylayHamilton())
            {
                m_pQ = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
                m_pQ2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
                m_pInverseSqrtQ = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
                checkCudaErrors(__cudaMalloc((void**)&m_pCaylayHamiltonCache, sizeof(DOUBLE) * _HC_Volume * _HC_Dir * 5));

                if (m_bProjDet)
                {
                    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceDetPhase, sizeof(DOUBLE) * _HC_Dir * _HC_Volume));
                }
            }
            else
            {
                TArray<gaugetype*> ptrs;
                for (UINT i = 0; i < m_uiRationalApproximationOrder + 1; ++i)
                {
                    m_pProjForce.AddItem(dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__)));
                    CFieldGaugeLink<gaugetype, matrixN>* pProjForce = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(m_pProjForce[i]);
                    ptrs.AddItem(pProjForce->m_pDeviceData);
                }
                checkCudaErrors(__cudaMalloc((void**)&m_pProjForcePtr, sizeof(gaugetype*) * (m_uiRationalApproximationOrder + 1)));
                checkCudaErrors(cudaMemcpy(m_pProjForcePtr, ptrs.GetData(), sizeof(gaugetype*) * (m_uiRationalApproximationOrder + 1), cudaMemcpyHostToDevice));

                if (m_bProjDet)
                {
                    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceDet, sizeof(CLGComplex) * _HC_Dir * _HC_Volume));
                }
            }
        }
    }

    //Level 1: p3 is cached, p5 is temporary
    {
        CFieldGauge* pP5_1_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_1_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_2_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_2_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_3_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_3_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

        CGaugeSmearingASQTAD<gaugetype, matrixN>::Fat357Lepage(
            pGauge,
            NULL,
            m_pP3_1_L1,
            m_pP3_2_L1,
            m_pP3_3_L1,
            pP5_1_1,
            pP5_1_2,
            pP5_2_1,
            pP5_2_2,
            pP5_3_1,
            pP5_3_2,
            m_fOriginalL1,
            m_fFat3L1,
            m_fFat5L1,
            m_fFat7L1,
            m_fLepageL1,
            pOrignalGauge
        );

        pP5_1_1->Return();
        pP5_1_2->Return();
        pP5_2_1->Return();
        pP5_2_2->Return();
        pP5_3_1->Return();
        pP5_3_2->Return();
    }

    if (m_bProj)
    {
        if (UseCaylayHamilton())
        {
            CGaugeSmearingASQTAD<gaugetype, matrixN>::ProjectCaylayHamilton(pGauge,
                m_pGaugeNotProjected,
                m_pQ,
                m_pQ2,
                m_pInverseSqrtQ,
                m_pCaylayHamiltonCache,
                m_bProjDet,
                m_pDeviceDetPhase);
        }
        else
        {
            CGaugeSmearingASQTAD<gaugetype, matrixN>::ProjectRationalApproximation(pGauge,
                m_pGaugeNotProjected,
                m_pProjForcePtr,
                m_pDeviceRationalApproximation,
                m_uiRationalApproximationOrder,
                m_bProjDet,
                m_pDeviceDet);
        }
    }

    //Level 2: p3 is cached, p5 is temporary
    {
        CFieldGauge* pP5_1_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_1_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_2_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_2_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_3_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_3_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

        CGaugeSmearingASQTAD<gaugetype, matrixN>::Fat357Lepage(pGauge,
            m_pEffectiveGaugeL1,
            m_pP3_1_L2,
            m_pP3_2_L2,
            m_pP3_3_L2,
            pP5_1_1,
            pP5_1_2,
            pP5_2_1,
            pP5_2_2,
            pP5_3_1,
            pP5_3_2,
            m_fOriginalL2,
            m_fFat3L2,
            m_fFat5L2,
            m_fFat7L2,
            m_fLepageL2,
            NULL
        );

        pP5_1_1->Return();
        pP5_1_2->Return();
        pP5_2_1->Return();
        pP5_2_2->Return();
        pP5_3_1->Return();
        pP5_3_2->Return();
    }

    if (NULL == m_pNaik)
    {
        m_pNaik = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_pEffectiveGaugeL1->m_byFieldId, _T(__FILE__), __LINE__));
    }
    CalculateNaikLink((gaugetype*)m_pNaik->GetData(), (const gaugetype*)m_pEffectiveGaugeL1->GetData(), m_pEffectiveGaugeL1->m_byFieldId);

    if (m_iPhaseFieldId > 0)
    {
        if (NULL == m_pNaikLinkPhase)
        {
            m_pNaikLinkPhase = dynamic_cast<CFieldGaugeU1Real*>(appGetLattice()->GetPooledFieldById(static_cast<BYTE>(m_iPhaseFieldId), _T(__FILE__), __LINE__));
            CalculateU1NaikLink(m_pNaikLinkPhase, dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(static_cast<BYTE>(m_iPhaseFieldId))));
        }
    }
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingHISQ<gaugetype, matrixN>::GaugeSmearingNone(CFieldGauge* pGauge, const class CFieldGauge* pOrignalGauge, CFieldGauge* pStaple, UBOOL bProject)
{
    _RECORD(CGaugeSmearingHISQ::GaugeSmearingNone);
    appParanoiac(_T("CGaugeSmearingHISQ::GaugeSmearingNone\n"));

    if (NULL == m_pEffectiveGaugeL1)
    {
        m_pEffectiveGaugeL1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

        if (m_bProj)
        {
            m_pGaugeNotProjected = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

            if (UseCaylayHamilton())
            {
                m_pQ = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
                m_pQ2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
                m_pInverseSqrtQ = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
                checkCudaErrors(__cudaMalloc((void**)&m_pCaylayHamiltonCache, sizeof(DOUBLE) * _HC_Volume * _HC_Dir * 5));

                if (m_bProjDet)
                {
                    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceDetPhase, sizeof(DOUBLE) * _HC_Dir * _HC_Volume));
                }
            }
            else
            {
                TArray<gaugetype*> ptrs;
                for (UINT i = 0; i < m_uiRationalApproximationOrder + 1; ++i)
                {
                    m_pProjForce.AddItem(dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__)));
                    CFieldGaugeLink<gaugetype, matrixN>* pProjForce = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(m_pProjForce[i]);
                    ptrs.AddItem(pProjForce->m_pDeviceData);
                }
                checkCudaErrors(__cudaMalloc((void**)&m_pProjForcePtr, sizeof(gaugetype*) * (m_uiRationalApproximationOrder + 1)));
                checkCudaErrors(cudaMemcpy(m_pProjForcePtr, ptrs.GetData(), sizeof(gaugetype*) * (m_uiRationalApproximationOrder + 1), cudaMemcpyHostToDevice));

                if (m_bProjDet)
                {
                    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceDet, sizeof(CLGComplex) * _HC_Dir * _HC_Volume));
                }
            }
        }
    }

    //Level 1: all staples are temporary
    {
        CFieldGauge* pP3_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP3_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP3_3 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_1_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_1_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_2_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_2_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_3_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_3_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

        CGaugeSmearingASQTAD<gaugetype, matrixN>::Fat357Lepage(
            pGauge,
            NULL,
            pP3_1,
            pP3_2,
            pP3_3,
            pP5_1_1,
            pP5_1_2,
            pP5_2_1,
            pP5_2_2,
            pP5_3_1,
            pP5_3_2,
            m_fOriginalL1,
            m_fFat3L1,
            m_fFat5L1,
            m_fFat7L1,
            m_fLepageL1,
            pOrignalGauge
        );

        pP3_1->Return();
        pP3_2->Return();
        pP3_3->Return();
        pP5_1_1->Return();
        pP5_1_2->Return();
        pP5_2_1->Return();
        pP5_2_2->Return();
        pP5_3_1->Return();
        pP5_3_2->Return();
    }

    if (m_bProj)
    {
        if (UseCaylayHamilton())
        {
            CGaugeSmearingASQTAD<gaugetype, matrixN>::ProjectCaylayHamilton(pGauge,
                m_pGaugeNotProjected,
                m_pQ,
                m_pQ2,
                m_pInverseSqrtQ,
                m_pCaylayHamiltonCache,
                m_bProjDet,
                m_pDeviceDetPhase);
        }
        else
        {
            CGaugeSmearingASQTAD<gaugetype, matrixN>::ProjectRationalApproximation(pGauge,
                m_pGaugeNotProjected,
                m_pProjForcePtr,
                m_pDeviceRationalApproximation,
                m_uiRationalApproximationOrder,
                m_bProjDet,
                m_pDeviceDet);
        }
    }

    //Level 2: all staples are temporary
    {
        CFieldGauge* pP3_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP3_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP3_3 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_1_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_1_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_2_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_2_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_3_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_3_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

        CGaugeSmearingASQTAD<gaugetype, matrixN>::Fat357Lepage(pGauge,
            m_pEffectiveGaugeL1,
            pP3_1,
            pP3_2,
            pP3_3,
            pP5_1_1,
            pP5_1_2,
            pP5_2_1,
            pP5_2_2,
            pP5_3_1,
            pP5_3_2,
            m_fOriginalL2,
            m_fFat3L2,
            m_fFat5L2,
            m_fFat7L2,
            m_fLepageL2,
            NULL
        );

        pP3_1->Return();
        pP3_2->Return();
        pP3_3->Return();
        pP5_1_1->Return();
        pP5_1_2->Return();
        pP5_2_1->Return();
        pP5_2_2->Return();
        pP5_3_1->Return();
        pP5_3_2->Return();
    }

    if (NULL == m_pNaik)
    {
        m_pNaik = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_pEffectiveGaugeL1->m_byFieldId, _T(__FILE__), __LINE__));
    }
    CalculateNaikLink((gaugetype*)m_pNaik->GetData(), (const gaugetype*)m_pEffectiveGaugeL1->GetData(), m_pEffectiveGaugeL1->m_byFieldId);

    if (m_iPhaseFieldId > 0)
    {
        if (NULL == m_pNaikLinkPhase)
        {
            m_pNaikLinkPhase = dynamic_cast<CFieldGaugeU1Real*>(appGetLattice()->GetPooledFieldById(static_cast<BYTE>(m_iPhaseFieldId), _T(__FILE__), __LINE__));
            CalculateU1NaikLink(m_pNaikLinkPhase, dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(static_cast<BYTE>(m_iPhaseFieldId))));
        }
    }
}

//m_pEffectiveGaugeL1 is unchaged
template<typename gaugetype, INT matrixN>
void CGaugeSmearingHISQ<gaugetype, matrixN>::CalculateNaikLink(gaugetype* naik, const gaugetype* effectivel1, BYTE byFieldId)
{
    _RECORD(CGaugeSmearingHISQ::CalculateNaikLink);
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelNaikLink<gaugetype>, block, threads, naik, effectivel1, byFieldId);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingHISQ<gaugetype, matrixN>::CalculateU1NaikLink(CFieldGaugeU1Real* naik, const CFieldGaugeU1Real* externalfield)
{
    _RECORD(CGaugeSmearingHISQ::CalculateU1NaikLink);
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelNaikLink<Real>, block, threads, (Real*)naik->GetData(), (const Real*)externalfield->GetData(), naik->m_byFieldId);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingHISQ<gaugetype, matrixN>::DerivateOnUFull(const CFieldGauge* pEffectiveGauge, const class CFieldGauge* pOrignalGauge, const class CFieldGauge* pNaikForce, CFieldGauge* pf0) const
{
    _RECORD(CGaugeSmearingHISQ::DerivateOnUFull);
    {
        _RECORD2(CGaugeSmearingHISQ::DerivateOnUFull::SmearingForce, a);
        CGaugeSmearingASQTAD<gaugetype, matrixN>::SmearingForce(m_pEffectiveGaugeL1,
            pf0,
            m_pP3_1_L2,
            m_pP3_2_L2,
            m_pP3_3_L2,
            m_pP5_1_1_L2,
            m_pP5_1_2_L2,
            m_pP5_2_1_L2,
            m_pP5_2_2_L2,
            m_pP5_3_1_L2,
            m_pP5_3_2_L2,
            m_fOriginalL2,
            m_fFat3L2,
            m_fFat5L2,
            m_fFat7L2,
            m_fLepageL2
        );
    }

    if (NULL != pNaikForce)
    {
        pf0->AxpyPlus(pNaikForce);
    }

    {
        _RECORD2(CGaugeSmearingHISQ::DerivateOnUFull::ProjectRationalApproximationForce, a);
        if (m_bProj)
        {
            if (UseCaylayHamilton())
            {
                CGaugeSmearingASQTAD<gaugetype, matrixN>::ProjectCaylayHamiltonForce(m_pEffectiveGaugeL1,
                    m_pGaugeNotProjected,
                    pf0,
                    m_pQ,
                    m_pQ2,
                    m_pInverseSqrtQ,
                    m_pCaylayHamiltonCache,
                    m_bProjDet,
                    m_pDeviceDetPhase);
            }
            else
            {
                CGaugeSmearingASQTAD<gaugetype, matrixN>::ProjectRationalApproximationForce(m_pEffectiveGaugeL1,
                    m_pGaugeNotProjected,
                    pf0,
                    m_pProjForcePtr,
                    m_pDeviceRationalApproximation,
                    m_uiRationalApproximationOrder,
                    m_bProjDet,
                    m_pDeviceDet);
            }
        }
    }

    {
        _RECORD2(CGaugeSmearingHISQ::DerivateOnUFull::SmearingForce, a);
        CGaugeSmearingASQTAD<gaugetype, matrixN>::SmearingForce(pOrignalGauge,
            pf0,
            m_pP3_1_L1,
            m_pP3_2_L1,
            m_pP3_3_L1,
            m_pP5_1_1_L1,
            m_pP5_1_2_L1,
            m_pP5_2_1_L1,
            m_pP5_2_2_L1,
            m_pP5_3_1_L1,
            m_pP5_3_2_L1,
            m_fOriginalL1,
            m_fFat3L1,
            m_fFat5L1,
            m_fFat7L1,
            m_fLepageL1
        );
    }
    _CHECKCUDA;
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingHISQ<gaugetype, matrixN>::DerivateOnUMedian(const CFieldGauge* pEffectiveGauge, const class CFieldGauge* pOrignalGauge, const class CFieldGauge* pNaikForce, CFieldGauge* pf0) const
{
    _RECORD(CGaugeSmearingHISQ::DerivateOnUMedian);

    //L2 force: p3_L2 cached, recompute p5_L2
    {
        CFieldGauge* pP5_1_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(m_pEffectiveGaugeL1, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_1_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(m_pEffectiveGaugeL1, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_2_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(m_pEffectiveGaugeL1, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_2_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(m_pEffectiveGaugeL1, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_3_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(m_pEffectiveGaugeL1, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_3_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(m_pEffectiveGaugeL1, _T(__FILE__), __LINE__));

        CGaugeSmearingASQTAD<gaugetype, matrixN>::CalculateFat5Only(
            m_pEffectiveGaugeL1,
            m_pP3_1_L2, m_pP3_2_L2, m_pP3_3_L2,
            pP5_1_1, pP5_1_2, pP5_2_1,
            pP5_2_2, pP5_3_1, pP5_3_2
        );

        CGaugeSmearingASQTAD<gaugetype, matrixN>::SmearingForce(m_pEffectiveGaugeL1,
            pf0,
            m_pP3_1_L2,
            m_pP3_2_L2,
            m_pP3_3_L2,
            pP5_1_1,
            pP5_1_2,
            pP5_2_1,
            pP5_2_2,
            pP5_3_1,
            pP5_3_2,
            m_fOriginalL2,
            m_fFat3L2,
            m_fFat5L2,
            m_fFat7L2,
            m_fLepageL2
        );

        pP5_1_1->Return();
        pP5_1_2->Return();
        pP5_2_1->Return();
        pP5_2_2->Return();
        pP5_3_1->Return();
        pP5_3_2->Return();
    }

    if (NULL != pNaikForce)
    {
        pf0->AxpyPlus(pNaikForce);
    }

    if (m_bProj)
    {
        if (UseCaylayHamilton())
        {
            CGaugeSmearingASQTAD<gaugetype, matrixN>::ProjectCaylayHamiltonForce(m_pEffectiveGaugeL1,
                m_pGaugeNotProjected,
                pf0,
                m_pQ,
                m_pQ2,
                m_pInverseSqrtQ,
                m_pCaylayHamiltonCache,
                m_bProjDet,
                m_pDeviceDetPhase);
        }
        else
        {
            CGaugeSmearingASQTAD<gaugetype, matrixN>::ProjectRationalApproximationForce(m_pEffectiveGaugeL1,
                m_pGaugeNotProjected,
                pf0,
                m_pProjForcePtr,
                m_pDeviceRationalApproximation,
                m_uiRationalApproximationOrder,
                m_bProjDet,
                m_pDeviceDet);
        }
    }

    //L1 force: p3_L1 cached, recompute p5_L1
    {
        CFieldGauge* pP5_1_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pOrignalGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_1_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pOrignalGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_2_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pOrignalGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_2_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pOrignalGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_3_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pOrignalGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_3_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pOrignalGauge, _T(__FILE__), __LINE__));

        CGaugeSmearingASQTAD<gaugetype, matrixN>::CalculateFat5Only(
            pOrignalGauge,
            m_pP3_1_L1, m_pP3_2_L1, m_pP3_3_L1,
            pP5_1_1, pP5_1_2, pP5_2_1,
            pP5_2_2, pP5_3_1, pP5_3_2
        );

        CGaugeSmearingASQTAD<gaugetype, matrixN>::SmearingForce(pOrignalGauge,
            pf0,
            m_pP3_1_L1,
            m_pP3_2_L1,
            m_pP3_3_L1,
            pP5_1_1,
            pP5_1_2,
            pP5_2_1,
            pP5_2_2,
            pP5_3_1,
            pP5_3_2,
            m_fOriginalL1,
            m_fFat3L1,
            m_fFat5L1,
            m_fFat7L1,
            m_fLepageL1
        );

        pP5_1_1->Return();
        pP5_1_2->Return();
        pP5_2_1->Return();
        pP5_2_2->Return();
        pP5_3_1->Return();
        pP5_3_2->Return();
    }
    _CHECKCUDA;
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingHISQ<gaugetype, matrixN>::DerivateOnUNone(const CFieldGauge* pEffectiveGauge, const class CFieldGauge* pOrignalGauge, const class CFieldGauge* pNaikForce, CFieldGauge* pf0) const
{
    _RECORD(CGaugeSmearingHISQ::DerivateOnUNone);

    //L2 force: recompute both p3_L2 and p5_L2
    {
        CFieldGauge* pP3_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(m_pEffectiveGaugeL1, _T(__FILE__), __LINE__));
        CFieldGauge* pP3_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(m_pEffectiveGaugeL1, _T(__FILE__), __LINE__));
        CFieldGauge* pP3_3 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(m_pEffectiveGaugeL1, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_1_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(m_pEffectiveGaugeL1, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_1_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(m_pEffectiveGaugeL1, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_2_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(m_pEffectiveGaugeL1, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_2_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(m_pEffectiveGaugeL1, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_3_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(m_pEffectiveGaugeL1, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_3_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(m_pEffectiveGaugeL1, _T(__FILE__), __LINE__));

        CGaugeSmearingASQTAD<gaugetype, matrixN>::CalculateFat3Only(
            m_pEffectiveGaugeL1,
            pP3_1, pP3_2, pP3_3
        );
        CGaugeSmearingASQTAD<gaugetype, matrixN>::CalculateFat5Only(
            m_pEffectiveGaugeL1,
            pP3_1, pP3_2, pP3_3,
            pP5_1_1, pP5_1_2, pP5_2_1,
            pP5_2_2, pP5_3_1, pP5_3_2
        );

        CGaugeSmearingASQTAD<gaugetype, matrixN>::SmearingForce(m_pEffectiveGaugeL1,
            pf0,
            pP3_1,
            pP3_2,
            pP3_3,
            pP5_1_1,
            pP5_1_2,
            pP5_2_1,
            pP5_2_2,
            pP5_3_1,
            pP5_3_2,
            m_fOriginalL2,
            m_fFat3L2,
            m_fFat5L2,
            m_fFat7L2,
            m_fLepageL2
        );

        pP3_1->Return();
        pP3_2->Return();
        pP3_3->Return();
        pP5_1_1->Return();
        pP5_1_2->Return();
        pP5_2_1->Return();
        pP5_2_2->Return();
        pP5_3_1->Return();
        pP5_3_2->Return();
    }

    if (NULL != pNaikForce)
    {
        pf0->AxpyPlus(pNaikForce);
    }

    if (m_bProj)
    {
        if (UseCaylayHamilton())
        {
            CGaugeSmearingASQTAD<gaugetype, matrixN>::ProjectCaylayHamiltonForce(m_pEffectiveGaugeL1,
                m_pGaugeNotProjected,
                pf0,
                m_pQ,
                m_pQ2,
                m_pInverseSqrtQ,
                m_pCaylayHamiltonCache,
                m_bProjDet,
                m_pDeviceDetPhase);
        }
        else
        {
            CGaugeSmearingASQTAD<gaugetype, matrixN>::ProjectRationalApproximationForce(m_pEffectiveGaugeL1,
                m_pGaugeNotProjected,
                pf0,
                m_pProjForcePtr,
                m_pDeviceRationalApproximation,
                m_uiRationalApproximationOrder,
                m_bProjDet,
                m_pDeviceDet);
        }
    }

    //L1 force: recompute both p3_L1 and p5_L1
    {
        CFieldGauge* pP3_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pOrignalGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP3_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pOrignalGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP3_3 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pOrignalGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_1_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pOrignalGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_1_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pOrignalGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_2_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pOrignalGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_2_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pOrignalGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_3_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pOrignalGauge, _T(__FILE__), __LINE__));
        CFieldGauge* pP5_3_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pOrignalGauge, _T(__FILE__), __LINE__));

        CGaugeSmearingASQTAD<gaugetype, matrixN>::CalculateFat3Only(
            pOrignalGauge,
            pP3_1, pP3_2, pP3_3
        );
        CGaugeSmearingASQTAD<gaugetype, matrixN>::CalculateFat5Only(
            pOrignalGauge,
            pP3_1, pP3_2, pP3_3,
            pP5_1_1, pP5_1_2, pP5_2_1,
            pP5_2_2, pP5_3_1, pP5_3_2
        );

        CGaugeSmearingASQTAD<gaugetype, matrixN>::SmearingForce(pOrignalGauge,
            pf0,
            pP3_1,
            pP3_2,
            pP3_3,
            pP5_1_1,
            pP5_1_2,
            pP5_2_1,
            pP5_2_2,
            pP5_3_1,
            pP5_3_2,
            m_fOriginalL1,
            m_fFat3L1,
            m_fFat5L1,
            m_fFat7L1,
            m_fLepageL1
        );

        pP3_1->Return();
        pP3_2->Return();
        pP3_3->Return();
        pP5_1_1->Return();
        pP5_1_2->Return();
        pP5_2_1->Return();
        pP5_2_2->Return();
        pP5_3_1->Return();
        pP5_3_2->Return();
    }
    _CHECKCUDA;
}

template<typename gaugetype, INT matrixN>
CCString CGaugeSmearingHISQ<gaugetype, matrixN>::GetInfos(const CCString &tab) const
{
    CCString sRet = CGaugeSmearing::GetInfos(tab);
    sRet = sRet + tab + _T("Original L1 : ") + appToString(m_fOriginalL1) + _T("\n");
    sRet = sRet + tab + _T("Fat3 L1     : ") + appToString(m_fFat3L1) + _T("\n");
    sRet = sRet + tab + _T("Fat5 L1     : ") + appToString(m_fFat5L1) + _T("\n");
    sRet = sRet + tab + _T("Fat7 L1     : ") + appToString(m_fFat7L1) + _T("\n");
    sRet = sRet + tab + _T("Lepage L1   : ") + appToString(m_fLepageL1) + _T("\n");
    
    sRet = sRet + tab + _T("Original L2 : ") + appToString(m_fOriginalL2) + _T("\n");
    sRet = sRet + tab + _T("Fat3 L2     : ") + appToString(m_fFat3L2) + _T("\n");
    sRet = sRet + tab + _T("Fat5 L2     : ") + appToString(m_fFat5L2) + _T("\n");
    sRet = sRet + tab + _T("Fat7 L2     : ") + appToString(m_fFat7L2) + _T("\n");
    sRet = sRet + tab + _T("Lepage L2   : ") + appToString(m_fLepageL2) + _T("\n");

    sRet = sRet + tab + _T("ProjectU    : ") + appToString(m_bProj) + _T("\n");
    sRet = sRet + tab + _T("ProjectSU3    : ") + appToString(m_bProjDet) + _T("\n");
    sRet = sRet + tab + _T("Using Caylay-Hamilton    : ") + appToString(m_bUseCaylayHamilton) + _T("\n");
    sRet = sRet + tab + _T("Cache       : ") + __ENUM_TO_STRING(EHISQLinkCache, m_eCache).c_str() + _T("\n");

    return sRet;
}

template class CGaugeSmearingHISQ<deviceSU3, 3>;
__CLGIMPLEMENT_CLASS(CGaugeSmearingHISQSU3)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================