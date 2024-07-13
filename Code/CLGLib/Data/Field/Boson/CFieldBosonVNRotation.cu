//=============================================================================
// FILENAME : CFieldBosonVNRotation.cpp
// 
// DESCRIPTION:
// This is the class for the spin fields
//
// REVISION:
//  [07/13/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "CFieldBosonVNRotation.h"

__BEGIN_NAMESPACE

__device__ Real testbosonlink(BYTE byFieldId, SSmallInt4 site, const SIndex& uiSiteBI)
{
    return F(1.0);
}

__device__ _deviceCoeffFunctionPointer pfttest = testbosonlink;

template<typename deviceDataBoson, typename deviceDataGauge>
CFieldBosonVNRotation<deviceDataBoson, deviceDataGauge>::CFieldBosonVNRotation()
    : CFieldBosonVN<deviceDataBoson, deviceDataGauge>()
    , m_fOmega(F(0.0))
    , m_pDevicePath(NULL)
    , m_pfCx(NULL)
{
    checkCudaErrors(cudaMalloc((void**)&m_pDevicePath, sizeof(INT) * _kLinkMaxLength));
    checkCudaErrors(cudaMemcpyFromSymbol(&m_pfCx, pfttest, sizeof(_deviceCoeffFunctionPointer)));
}

template<typename deviceDataBoson, typename deviceDataGauge>
CFieldBosonVNRotation<deviceDataBoson, deviceDataGauge>::~CFieldBosonVNRotation()
{
    checkCudaErrors(cudaFree(m_pDevicePath));
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNRotation<deviceDataBoson, deviceDataGauge>::CopyTo(CField* U) const
{
    CFieldBosonVN<deviceDataBoson, deviceDataGauge>::CopyTo(U);
    CFieldBosonVNRotation<deviceDataBoson, deviceDataGauge>* pU = static_cast<CFieldBosonVNRotation<deviceDataBoson, deviceDataGauge>*>(U);

    pU->m_fOmega = m_fOmega;
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNRotation<deviceDataBoson, deviceDataGauge>::DFromSource(const CFieldBoson* pSource, INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGaugeFields, const CFieldBoson* const* pBoson,
    EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    //SSmallInt4 x;
    //x.x = -1;
    //x.y = -1;
    //x.z = -1;
    //x.w = -1;
    //CIndexData::DebugEdgeMapping(m_byFieldId, x);


    CFieldBosonVN<deviceDataBoson, deviceDataGauge>::DFromSource(pSource, gaugeNum, bosonNum, pGaugeFields, pBoson, eCoeffType, fCoeffReal, fCoeffImg);

    const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, pGaugeFields);

    if (NULL != pGauge && VectorN() != pGauge->MatrixN())
    {
        appCrucial(_T("CFieldBosonVU can only play with gauge UN!"));
        return;
    }

    if (NULL == pSource || GetFieldType() != pSource->GetFieldType())
    {
        appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
        return;
    }

    const CFieldBosonVN<deviceDataBoson, deviceDataGauge>* pSourceVN = dynamic_cast<const CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(pSource);
    if (NULL == pSourceVN)
    {
        appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
        return;
    }

    //try external gauge field
    if (NULL == pGauge && m_byGaugeFieldIds.Num() > 0)
    {
        const CField* externelfield = appGetLattice()->GetFieldById(m_byGaugeFieldIds[0]);
        if (NULL != externelfield)
        {
            pGauge = dynamic_cast<const CFieldGauge*>(appGetLattice()->GetFieldById(m_byGaugeFieldIds[0]));
            if (pGauge->IsDynamic())
            {
                appCrucial(_T("CFieldBosonUN: A dynamic field is configured for this UN, but not for the action!\n"));
                pGauge = NULL;
            }
        }
    }

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    INT hostPath[3] = { 1, 2, -3 };
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 3, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    OneLinkS(
        pSourceVN->m_pDeviceData,
        NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
        NULL == pGauge ? 0 : pGauge->m_byFieldId,
        m_fOmega,
        m_pfCx,
        m_pDevicePath,
        3,
        eCoeffType,
        fRealCoeff,
        cCompCoeff);

    checkCudaErrors(cudaDeviceSynchronize());
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNRotation<deviceDataBoson, deviceDataGauge>::ForceOnGauge(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, CFieldGauge* const* pGaugeForce, const CFieldBoson* const* pBoson) const
{
    CFieldBosonVN<deviceDataBoson, deviceDataGauge>::ForceOnGauge(gaugeNum, bosonNum, pGauge, pGaugeForce, pBoson);

    if (m_byGaugeFieldIds.Num() < 1)
    {
        appCrucial(_T("CFieldBosonUN ForceOnGauge: there is no gauge!"));
        return;
    }
    INT gaugeidx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, pGauge, m_byGaugeFieldIds[0]);
    if (gaugeidx < 0 || gaugeidx >= m_byGaugeFieldIds.Num())
    {
        appCrucial(_T("CFieldBosonUN ForceOnGauge: there is no gauge!"));
        return;
    }

    const CFieldGauge* gauge = pGauge[gaugeidx];
    CFieldGauge* gaugeforce = pGaugeForce[gaugeidx];

    if (NULL == gauge || VectorN() != gauge->MatrixN())
    {
        appCrucial(_T("CFieldBosonUN can only play with gauge UN!"));
        return;
    }

    if (NULL == gaugeforce || VectorN() != gaugeforce->MatrixN())
    {
        appCrucial(_T("CFieldBosonUN can only play with gauge UN!"));
        return;
    }

    INT hostPath[3] = { 1, 2, -3 };
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 3, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());

    OneLinkForceGauge(
        (const deviceDataGauge*)gauge->GetData(),
        gauge->m_byFieldId,
        (deviceDataGauge*)gaugeforce->GetData(),
        m_fOmega,
        m_pfCx,
        m_pDevicePath,
        3
    );

    //gaugeu1->DebugPrintMe();
    checkCudaErrors(cudaDeviceSynchronize());
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNRotation<deviceDataBoson, deviceDataGauge>::InitialOtherParameters(CParameters& param)
{
    CFieldBosonVN<deviceDataBoson, deviceDataGauge>::InitialOtherParameters(param);

    Real fValue = F(0.1);
    if (param.FetchValueReal(_T("Omega"), fValue))
    {
        m_fOmega = fValue;
    }
}

template<typename deviceDataBoson, typename deviceDataGauge>
CCString CFieldBosonVNRotation<deviceDataBoson, deviceDataGauge>::GetInfos(const CCString& tab) const
{
    CCString sRet = CFieldBosonVN<deviceDataBoson, deviceDataGauge>::GetInfos(tab);
    sRet = sRet + appToString(m_fOmega) + _T("\n");
    return sRet;
}

__CLGIMPLEMENT_CLASS(CFieldBosonVNRotationU1)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================