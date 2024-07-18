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

__device__ Real _deviceBosonRotationCx(BYTE byFieldId, SSmallInt4 site1, SSmallInt4 site2, const SIndex& uiSiteBI1, const SIndex& uiSiteBI2)
{
    Real fX1, fX2;
    if (uiSiteBI1.IsDirichlet())
    {
        fX1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fX1 = static_cast<Real>(site1.x) - _DC_Centerx;
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fX2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fX2 = static_cast<Real>(site2.x) - _DC_Centerx;
    }

    return F(0.5) * (fX1 + fX2);
}

__device__ Real _deviceBosonRotationCy(BYTE byFieldId, SSmallInt4 site1, SSmallInt4 site2, const SIndex& uiSiteBI1, const SIndex& uiSiteBI2)
{
    Real fY1, fY2;
    if (uiSiteBI1.IsDirichlet())
    {
        fY1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fY1 = static_cast<Real>(site1.y) - _DC_Centery;
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fY2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fY2 = static_cast<Real>(site2.y) - _DC_Centery;
    }

    return F(0.5) * (fY1 + fY2);
}

__device__ Real _deviceBosonRotationCxy(BYTE byFieldId, SSmallInt4 site1, SSmallInt4 site2, const SIndex& uiSiteBI1, const SIndex& uiSiteBI2)
{
    Real fX1, fX2;
    if (uiSiteBI1.IsDirichlet())
    {
        fX1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fX1 = static_cast<Real>(site1.x) - _DC_Centerx;
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fX2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fX2 = static_cast<Real>(site2.x) - _DC_Centerx;
    }

    Real fY1, fY2;
    if (uiSiteBI1.IsDirichlet())
    {
        fY1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fY1 = static_cast<Real>(site1.y) - _DC_Centery;
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fY2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fY2 = static_cast<Real>(site2.y) - _DC_Centery;
    }

    return F(0.25) * ((fX1 + fX2) * (fY1 + fY2));
}

__device__ Real _deviceBosonRotationCxShift(BYTE byFieldId, SSmallInt4 site1, SSmallInt4 site2, const SIndex& uiSiteBI1, const SIndex& uiSiteBI2)
{
    Real fX1, fX2;
    if (uiSiteBI1.IsDirichlet())
    {
        fX1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fX1 = static_cast<Real>(site1.x) - _DC_Centerx + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fX2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fX2 = static_cast<Real>(site2.x) - _DC_Centerx + F(0.5);
    }

    return F(0.5) * (fX1 + fX2);
}

__device__ Real _deviceBosonRotationCyShift(BYTE byFieldId, SSmallInt4 site1, SSmallInt4 site2, const SIndex& uiSiteBI1, const SIndex& uiSiteBI2)
{
    Real fY1, fY2;
    if (uiSiteBI1.IsDirichlet())
    {
        fY1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fY1 = static_cast<Real>(site1.y) - _DC_Centery + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fY2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fY2 = static_cast<Real>(site2.y) - _DC_Centery + F(0.5);
    }

    return F(0.5) * (fY1 + fY2);
}

/**
* first x then y
*/
__device__ Real _deviceBosonRotationCxyShiftXYPP(BYTE byFieldId, SSmallInt4 site1, SSmallInt4 site2, const SIndex& uiSiteBI1, const SIndex& uiSiteBI2)
{
    UBOOL bOppsite = FALSE;
    if (_DC_Lxi == site2.x || -1 == site2.y || _DC_Lxi == site1.x || -1 == site1.y)
    {
        bOppsite = TRUE;
    }

    Real fX1, fX2;
    if (uiSiteBI1.IsDirichlet())
    {
        fX1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fX1 = static_cast<Real>(site1.x) - _DC_Centerx + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fX2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fX2 = static_cast<Real>(site2.x) - _DC_Centerx + F(0.5);
    }

    Real fY1, fY2;
    if (uiSiteBI1.IsDirichlet())
    {
        fY1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fY1 = static_cast<Real>(site1.y) - _DC_Centery + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fY2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fY2 = static_cast<Real>(site2.y) - _DC_Centery + F(0.5);
    }

    if (bOppsite)
    {
        return F(-0.25) * ((fX1 + fX2) * (fY1 + fY2));
    }
    return F(0.25)* ((fX1 + fX2) * (fY1 + fY2));
}

__device__ Real _deviceBosonRotationCxyShiftYXPP(BYTE byFieldId, SSmallInt4 site1, SSmallInt4 site2, const SIndex& uiSiteBI1, const SIndex& uiSiteBI2)
{
    UBOOL bOppsite = FALSE;
    if (_DC_Lyi == site2.y || -1 == site2.x || _DC_Lyi == site1.y || -1 == site1.x)
    {
        bOppsite = TRUE;
    }

    Real fX1, fX2;
    if (uiSiteBI1.IsDirichlet())
    {
        fX1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fX1 = static_cast<Real>(site1.x) - _DC_Centerx + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fX2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fX2 = static_cast<Real>(site2.x) - _DC_Centerx + F(0.5);
    }

    Real fY1, fY2;
    if (uiSiteBI1.IsDirichlet())
    {
        fY1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fY1 = static_cast<Real>(site1.y) - _DC_Centery + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fY2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fY2 = static_cast<Real>(site2.y) - _DC_Centery + F(0.5);
    }

    if (bOppsite)
    {
        return F(-0.25) * ((fX1 + fX2) * (fY1 + fY2));
    }
    return F(0.25) * ((fX1 + fX2) * (fY1 + fY2));
}

__device__ Real _deviceBosonRotationCxyShiftYXPM(BYTE byFieldId, SSmallInt4 site1, SSmallInt4 site2, const SIndex& uiSiteBI1, const SIndex& uiSiteBI2)
{
    UBOOL bOppsite = FALSE;
    if (_DC_Lyi == site2.y || _DC_Lxi == site2.x || _DC_Lyi == site1.y || _DC_Lxi == site1.x)
    {
        bOppsite = TRUE;
    }

    Real fX1, fX2;
    if (uiSiteBI1.IsDirichlet())
    {
        fX1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fX1 = static_cast<Real>(site1.x) - _DC_Centerx + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fX2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fX2 = static_cast<Real>(site2.x) - _DC_Centerx + F(0.5);
    }

    Real fY1, fY2;
    if (uiSiteBI1.IsDirichlet())
    {
        fY1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fY1 = static_cast<Real>(site1.y) - _DC_Centery + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fY2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fY2 = static_cast<Real>(site2.y) - _DC_Centery + F(0.5);
    }

    if (bOppsite)
    {
        return F(-0.25) * ((fX1 + fX2) * (fY1 + fY2));
    }
    return F(0.25) * ((fX1 + fX2) * (fY1 + fY2));
}

__device__ Real _deviceBosonRotationCxyShiftXYMP(BYTE byFieldId, SSmallInt4 site1, SSmallInt4 site2, const SIndex& uiSiteBI1, const SIndex& uiSiteBI2)
{
    UBOOL bOppsite = FALSE;
    if (-1 == site2.y || -1 == site2.x || -1 == site1.y || -1 == site1.x)
    {
        bOppsite = TRUE;
    }

    Real fX1, fX2;
    if (uiSiteBI1.IsDirichlet())
    {
        fX1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fX1 = static_cast<Real>(site1.x) - _DC_Centerx + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fX2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fX2 = static_cast<Real>(site2.x) - _DC_Centerx + F(0.5);
    }

    Real fY1, fY2;
    if (uiSiteBI1.IsDirichlet())
    {
        fY1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fY1 = static_cast<Real>(site1.y) - _DC_Centery + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fY2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fY2 = static_cast<Real>(site2.y) - _DC_Centery + F(0.5);
    }

    if (bOppsite)
    {
        return F(-0.25) * ((fX1 + fX2) * (fY1 + fY2));
    }
    return F(0.25) * ((fX1 + fX2) * (fY1 + fY2));
}

__device__ Real _deviceBosonRotationCxCySq(BYTE byFieldId, SSmallInt4 site1, const SIndex& uiSiteBI1)
{
    const Real x = static_cast<Real>(site1.x) - _DC_Centerx;
    const Real y = static_cast<Real>(site1.y) - _DC_Centery;
    return x * x + y * y;
}

__device__ Real _deviceBosonRotationCxCySqShift(BYTE byFieldId, SSmallInt4 site1, const SIndex& uiSiteBI1)
{
    const Real x = static_cast<Real>(site1.x) - _DC_Centerx + F(0.5);
    const Real y = static_cast<Real>(site1.y) - _DC_Centery + F(0.5);
    return x * x + y * y;
}


__device__ Real _deviceBosonRotationCxSq(BYTE byFieldId, SSmallInt4 site1, const SIndex& uiSiteBI1)
{
    if (site1.Out())
    {
        site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
    }
    const Real x = static_cast<Real>(site1.x) - _DC_Centerx;
    return x * x;
}

__device__ Real _deviceBosonRotationCySq(BYTE byFieldId, SSmallInt4 site1, const SIndex& uiSiteBI1)
{
    if (site1.Out())
    {
        site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
    }
    const Real y = static_cast<Real>(site1.y) - _DC_Centery;
    return y * y;
}

__device__ Real _deviceBosonRotationCxSqShift(BYTE byFieldId, SSmallInt4 site1, const SIndex& uiSiteBI1)
{
    if (site1.Out())
    {
        site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
    }
    const Real x = static_cast<Real>(site1.x) - _DC_Centerx + F(0.5);
    return x * x;
}

__device__ Real _deviceBosonRotationCySqShift(BYTE byFieldId, SSmallInt4 site1, const SIndex& uiSiteBI1)
{
    if (site1.Out())
    {
        site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
    }
    const Real y = static_cast<Real>(site1.y) - _DC_Centery + F(0.5);
    return y * y;
}

__device__ _deviceCoeffFunctionPointerTwoSites _deviceBosonRotationCxPf = _deviceBosonRotationCx;
__device__ _deviceCoeffFunctionPointerTwoSites _deviceBosonRotationCyPf = _deviceBosonRotationCy;
__device__ _deviceCoeffFunctionPointerTwoSites _deviceBosonRotationCxyPf = _deviceBosonRotationCxy;
__device__ _deviceCoeffFunctionPointerTwoSites _deviceBosonRotationCxPfShift = _deviceBosonRotationCxShift;
__device__ _deviceCoeffFunctionPointerTwoSites _deviceBosonRotationCyPfShift = _deviceBosonRotationCyShift;
__device__ _deviceCoeffFunctionPointerTwoSites _deviceBosonRotationCxyPfShiftXYPP = _deviceBosonRotationCxyShiftXYPP;
__device__ _deviceCoeffFunctionPointerTwoSites _deviceBosonRotationCxyPfShiftYXPP = _deviceBosonRotationCxyShiftYXPP;
__device__ _deviceCoeffFunctionPointerTwoSites _deviceBosonRotationCxyPfShiftYXPM = _deviceBosonRotationCxyShiftYXPM;
__device__ _deviceCoeffFunctionPointerTwoSites _deviceBosonRotationCxyPfShiftXYMP = _deviceBosonRotationCxyShiftXYMP;

__device__ _deviceCoeffFunctionPointer _deviceBosonRotationCxCySqPf = _deviceBosonRotationCxCySq;
__device__ _deviceCoeffFunctionPointer _deviceBosonRotationCxCySqShiftPf = _deviceBosonRotationCxCySqShift;
__device__ _deviceCoeffFunctionPointer _deviceBosonRotationCxSqPf = _deviceBosonRotationCxSq;
__device__ _deviceCoeffFunctionPointer _deviceBosonRotationCySqPf = _deviceBosonRotationCySq;
__device__ _deviceCoeffFunctionPointer _deviceBosonRotationCxSqShiftPf = _deviceBosonRotationCxSqShift;
__device__ _deviceCoeffFunctionPointer _deviceBosonRotationCySqShiftPf = _deviceBosonRotationCySqShift;

template<typename deviceDataBoson, typename deviceDataGauge>
CFieldBosonVNRotation<deviceDataBoson, deviceDataGauge>::CFieldBosonVNRotation()
    : CFieldBosonVN<deviceDataBoson, deviceDataGauge>()
    , m_fOmega(F(0.0))
    , m_pDevicePath(NULL)
    , m_bShiftCenter(FALSE)
    , m_pfCx(NULL)
    , m_pfCy(NULL)
    , m_pfCxy(NULL)
    , m_pfCxShift(NULL)
    , m_pfCyShift(NULL)
    , m_pfCxyShiftXYPP(NULL)
    , m_pfCxyShiftYXPP(NULL)
    , m_pfCxyShiftYXPM(NULL)
    , m_pfCxyShiftXYMP(NULL)
    , m_pfCxCySq(NULL)
    , m_pfCxCySqShift(NULL)
    , m_pfCxSq(NULL)
    , m_pfCySq(NULL)
    , m_pfCxSqShift(NULL)
    , m_pfCySqShift(NULL)
{
    checkCudaErrors(cudaMalloc((void**)&m_pDevicePath, sizeof(INT) * _kLinkMaxLength));

    checkCudaErrors(cudaMemcpyFromSymbol(&m_pfCx, _deviceBosonRotationCxPf, sizeof(_deviceCoeffFunctionPointerTwoSites)));
    checkCudaErrors(cudaMemcpyFromSymbol(&m_pfCy, _deviceBosonRotationCyPf, sizeof(_deviceCoeffFunctionPointerTwoSites)));
    checkCudaErrors(cudaMemcpyFromSymbol(&m_pfCxy, _deviceBosonRotationCxyPf, sizeof(_deviceCoeffFunctionPointerTwoSites)));
    checkCudaErrors(cudaMemcpyFromSymbol(&m_pfCxShift, _deviceBosonRotationCxPfShift, sizeof(_deviceCoeffFunctionPointerTwoSites)));
    checkCudaErrors(cudaMemcpyFromSymbol(&m_pfCyShift, _deviceBosonRotationCyPfShift, sizeof(_deviceCoeffFunctionPointerTwoSites)));
    checkCudaErrors(cudaMemcpyFromSymbol(&m_pfCxyShiftXYPP, _deviceBosonRotationCxyPfShiftXYPP, sizeof(_deviceCoeffFunctionPointerTwoSites)));
    checkCudaErrors(cudaMemcpyFromSymbol(&m_pfCxyShiftYXPP, _deviceBosonRotationCxyPfShiftYXPP, sizeof(_deviceCoeffFunctionPointerTwoSites)));
    checkCudaErrors(cudaMemcpyFromSymbol(&m_pfCxyShiftYXPM, _deviceBosonRotationCxyPfShiftYXPM, sizeof(_deviceCoeffFunctionPointerTwoSites)));
    checkCudaErrors(cudaMemcpyFromSymbol(&m_pfCxyShiftXYMP, _deviceBosonRotationCxyPfShiftXYMP, sizeof(_deviceCoeffFunctionPointerTwoSites)));

    checkCudaErrors(cudaMemcpyFromSymbol(&m_pfCxCySq, _deviceBosonRotationCxCySqPf, sizeof(_deviceCoeffFunctionPointer)));
    checkCudaErrors(cudaMemcpyFromSymbol(&m_pfCxCySqShift, _deviceBosonRotationCxCySqShiftPf, sizeof(_deviceCoeffFunctionPointer)));

    checkCudaErrors(cudaMemcpyFromSymbol(&m_pfCxSq, _deviceBosonRotationCxSqPf, sizeof(_deviceCoeffFunctionPointer)));
    checkCudaErrors(cudaMemcpyFromSymbol(&m_pfCySq, _deviceBosonRotationCySqPf, sizeof(_deviceCoeffFunctionPointer)));
    checkCudaErrors(cudaMemcpyFromSymbol(&m_pfCxSqShift, _deviceBosonRotationCxSqShiftPf, sizeof(_deviceCoeffFunctionPointer)));
    checkCudaErrors(cudaMemcpyFromSymbol(&m_pfCySqShift, _deviceBosonRotationCySqShiftPf, sizeof(_deviceCoeffFunctionPointer)));
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
    if (NULL != pU)
    {
        pU->m_fOmega = m_fOmega;
        pU->m_bShiftCenter = m_bShiftCenter;
    }
}

/**
* C [partial_{mu}phi partial_{nu}phi + partial_{nu}phi partial_{mu}phi]
* = C phi^2/a^2  - C (1/4a^2) 2Re[phi(n+mu) to phi(n+nu) - phi(n+mu) to phi(n-nu) - phi(n-mu) to phi(n+nu) + phi(n-mu) to phi(n-nu)]
* the 2 Re is in fact sum of two terms partial_{mu}phi partial_{nu}phi and partial_{nu}phi partial_{mu}phi which are Hermitian of each other
* 
* partial _x phi partial _x phi
* = + (2/a^2)phi^2 - phi(n)[U_{mu} phi(n+mu) + U_{-mu} phi(n-mu)]
* 
*/
template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNRotation<deviceDataBoson, deviceDataGauge>::DFromSource(const CFieldBoson* pSource, INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGaugeFields, const CFieldBoson* const* pBoson,
    EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    CFieldBosonVN<deviceDataBoson, deviceDataGauge>::DFromSource(pSource, gaugeNum, bosonNum, pGaugeFields, pBoson, eCoeffType, fCoeffReal, fCoeffImg);

    const CFieldGauge* pGauge = this->GetDefaultGauge(gaugeNum, pGaugeFields);

    if (NULL != pGauge && this->VectorN() != pGauge->MatrixN())
    {
        appCrucial(_T("CFieldBosonVU can only play with gauge UN!"));
        return;
    }

    if (NULL == pSource || this->GetFieldType() != pSource->GetFieldType())
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
    if (NULL == pGauge && this->m_byGaugeFieldIds.Num() > 0)
    {
        const CField* externelfield = appGetLattice()->GetFieldById(this->m_byGaugeFieldIds[0]);
        if (NULL != externelfield)
        {
            pGauge = dynamic_cast<const CFieldGauge*>(appGetLattice()->GetFieldById(this->m_byGaugeFieldIds[0]));
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

#pragma region +2 Re[Omega y pt px ]

    //====== +2 Re[Omega y pt px ]
    INT hostPath[2] = { 1, 4 };
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    OneLink(
        pSourceVN->m_pDeviceData,
        NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
        NULL == pGauge ? 0 : pGauge->m_byFieldId,
        -m_fOmega,
        m_bShiftCenter ? m_pfCyShift : m_pfCy,
        m_pDevicePath,
        2,
        eCoeffType,
        fRealCoeff,
        cCompCoeff);

    hostPath[0] = 4; hostPath[1] = 1;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    OneLink(
        pSourceVN->m_pDeviceData,
        NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
        NULL == pGauge ? 0 : pGauge->m_byFieldId,
        -m_fOmega,
        m_bShiftCenter ? m_pfCyShift : m_pfCy,
        m_pDevicePath,
        2,
        eCoeffType,
        fRealCoeff,
        cCompCoeff);

    hostPath[0] = 4; hostPath[1] = -1;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    OneLink(
        pSourceVN->m_pDeviceData,
        NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
        NULL == pGauge ? 0 : pGauge->m_byFieldId,
        m_fOmega,
        m_bShiftCenter ? m_pfCyShift : m_pfCy,
        m_pDevicePath,
        2,
        eCoeffType,
        fRealCoeff,
        cCompCoeff);

    hostPath[0] = -1; hostPath[1] = 4;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    OneLink(
        pSourceVN->m_pDeviceData,
        NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
        NULL == pGauge ? 0 : pGauge->m_byFieldId,
        m_fOmega,
        m_bShiftCenter ? m_pfCyShift : m_pfCy,
        m_pDevicePath,
        2,
        eCoeffType,
        fRealCoeff,
        cCompCoeff);

#pragma endregion

#pragma region -2 Re[Omega x pt py ]

    //====== -2 Re[Omega x pt py ]
    hostPath[0] = 2; hostPath[1] = 4;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    OneLink(
        pSourceVN->m_pDeviceData,
        NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
        NULL == pGauge ? 0 : pGauge->m_byFieldId,
        m_fOmega,
        m_bShiftCenter ? m_pfCxShift : m_pfCx,
        m_pDevicePath,
        2,
        eCoeffType,
        fRealCoeff,
        cCompCoeff);

    hostPath[0] = 4; hostPath[1] = 2;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    OneLink(
        pSourceVN->m_pDeviceData,
        NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
        NULL == pGauge ? 0 : pGauge->m_byFieldId,
        m_fOmega,
        m_bShiftCenter ? m_pfCxShift : m_pfCx,
        m_pDevicePath,
        2,
        eCoeffType,
        fRealCoeff,
        cCompCoeff);

    hostPath[0] = 4; hostPath[1] = -2;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    OneLink(
        pSourceVN->m_pDeviceData,
        NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
        NULL == pGauge ? 0 : pGauge->m_byFieldId,
        -m_fOmega,
        m_bShiftCenter ? m_pfCxShift : m_pfCx,
        m_pDevicePath,
        2,
        eCoeffType,
        fRealCoeff,
        cCompCoeff);

    hostPath[0] = -2; hostPath[1] = 4;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    OneLink(
        pSourceVN->m_pDeviceData,
        NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
        NULL == pGauge ? 0 : pGauge->m_byFieldId,
        -m_fOmega,
        m_bShiftCenter ? m_pfCxShift : m_pfCx,
        m_pDevicePath,
        2,
        eCoeffType,
        fRealCoeff,
        cCompCoeff);

#pragma endregion

#pragma region -2 Re[Omega x y px py ]

    //====== -2 Re[Omega x y px py ]
    hostPath[0] = 1; hostPath[1] = 2;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    OneLink(
        pSourceVN->m_pDeviceData,
        NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
        NULL == pGauge ? 0 : pGauge->m_byFieldId,
        m_fOmega * m_fOmega,
        m_bShiftCenter ? m_pfCxyShiftXYPP : m_pfCxy,
        m_pDevicePath,
        2,
        eCoeffType,
        fRealCoeff,
        cCompCoeff);

    hostPath[0] = 2; hostPath[1] = 1;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    OneLink(
        pSourceVN->m_pDeviceData,
        NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
        NULL == pGauge ? 0 : pGauge->m_byFieldId,
        m_fOmega * m_fOmega,
        m_bShiftCenter ? m_pfCxyShiftYXPP : m_pfCxy,
        m_pDevicePath,
        2,
        eCoeffType,
        fRealCoeff,
        cCompCoeff);

    hostPath[0] = 2; hostPath[1] = -1;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    OneLink(
        pSourceVN->m_pDeviceData,
        NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
        NULL == pGauge ? 0 : pGauge->m_byFieldId,
        -m_fOmega * m_fOmega,
        m_bShiftCenter ? m_pfCxyShiftYXPM : m_pfCxy,
        m_pDevicePath,
        2,
        eCoeffType,
        fRealCoeff,
        cCompCoeff);

    hostPath[0] = -1; hostPath[1] = 2;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    OneLink(
        pSourceVN->m_pDeviceData,
        NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
        NULL == pGauge ? 0 : pGauge->m_byFieldId,
        -m_fOmega * m_fOmega,
        m_bShiftCenter ? m_pfCxyShiftXYMP : m_pfCxy,
        m_pDevicePath,
        2,
        eCoeffType,
        fRealCoeff,
        cCompCoeff);

#pragma endregion

#pragma region +Omega^2 y^2 px^2 + Omega^2 x^2 py^2, square term

    DiagnalTerm(
        pSourceVN->m_pDeviceData, 
        -F(2.0) * m_fOmega * m_fOmega,
        m_bShiftCenter ? m_pfCxCySqShift : m_pfCxCySq,
        eCoeffType,
        fRealCoeff,
        cCompCoeff);

#pragma endregion

#pragma region +Omega^2 y^2 px^2 + Omega^2 x^2 py^2, shift term

    PartialSq(
        pSourceVN->m_pDeviceData,
        NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
        NULL == pGauge ? 0 : pGauge->m_byFieldId,
        m_fOmega * m_fOmega,
        m_bShiftCenter ? m_pfCySqShift : m_pfCySq,
        0,
        eCoeffType,
        fRealCoeff,
        cCompCoeff
    );

    PartialSq(
        pSourceVN->m_pDeviceData,
        NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
        NULL == pGauge ? 0 : pGauge->m_byFieldId,
        m_fOmega* m_fOmega,
        m_bShiftCenter ? m_pfCxSqShift : m_pfCxSq,
        1,
        eCoeffType,
        fRealCoeff,
        cCompCoeff
    );

#pragma endregion

    checkCudaErrors(cudaDeviceSynchronize());
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNRotation<deviceDataBoson, deviceDataGauge>::ForceOnGauge(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, CFieldGauge* const* pGaugeForce, const CFieldBoson* const* pBoson) const
{
    CFieldBosonVN<deviceDataBoson, deviceDataGauge>::ForceOnGauge(gaugeNum, bosonNum, pGauge, pGaugeForce, pBoson);

    if (this->m_byGaugeFieldIds.Num() < 1)
    {
        appCrucial(_T("CFieldBosonUN ForceOnGauge: there is no gauge!"));
        return;
    }
    INT gaugeidx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, pGauge, this->m_byGaugeFieldIds[0]);
    if (gaugeidx < 0 || gaugeidx >= this->m_byGaugeFieldIds.Num())
    {
        appCrucial(_T("CFieldBosonUN ForceOnGauge: there is no gauge!"));
        return;
    }

    const CFieldGauge* gauge = pGauge[gaugeidx];
    CFieldGauge* gaugeforce = pGaugeForce[gaugeidx];

    if (NULL == gauge || this->VectorN() != gauge->MatrixN())
    {
        appCrucial(_T("CFieldBosonUN can only play with gauge UN!"));
        return;
    }

    if (NULL == gaugeforce || this->VectorN() != gaugeforce->MatrixN())
    {
        appCrucial(_T("CFieldBosonUN can only play with gauge UN!"));
        return;
    }

#pragma region +2 Re[Omega y pt px ]

    INT hostPath[2] = { 1, 4 };
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    OneLinkForceGauge(
        (const deviceDataGauge*)gauge->GetData(),
        gauge->m_byFieldId,
        (deviceDataGauge*)gaugeforce->GetData(),
        -m_fOmega,
        m_bShiftCenter ? m_pfCyShift : m_pfCy,
        m_pDevicePath,
        2
    );

    hostPath[0] = 4; hostPath[1] = 1;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    OneLinkForceGauge(
        (const deviceDataGauge*)gauge->GetData(),
        gauge->m_byFieldId,
        (deviceDataGauge*)gaugeforce->GetData(),
        -m_fOmega,
        m_bShiftCenter ? m_pfCyShift : m_pfCy,
        m_pDevicePath,
        2
    );

    hostPath[0] = 4; hostPath[1] = -1;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    OneLinkForceGauge(
        (const deviceDataGauge*)gauge->GetData(),
        gauge->m_byFieldId,
        (deviceDataGauge*)gaugeforce->GetData(),
        m_fOmega,
        m_bShiftCenter ? m_pfCyShift : m_pfCy,
        m_pDevicePath,
        2
    );

    hostPath[0] = -1; hostPath[1] = 4;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    OneLinkForceGauge(
        (const deviceDataGauge*)gauge->GetData(),
        gauge->m_byFieldId,
        (deviceDataGauge*)gaugeforce->GetData(),
        m_fOmega,
        m_bShiftCenter ? m_pfCyShift : m_pfCy,
        m_pDevicePath,
        2
    );

#pragma endregion

#pragma region -2 Re[Omega x pt py ]

    hostPath[0] = 2; hostPath[1] = 4;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    OneLinkForceGauge(
        (const deviceDataGauge*)gauge->GetData(),
        gauge->m_byFieldId,
        (deviceDataGauge*)gaugeforce->GetData(),
        m_fOmega,
        m_bShiftCenter ? m_pfCxShift : m_pfCx,
        m_pDevicePath,
        2
    );

    hostPath[0] = 4; hostPath[1] = 2;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    OneLinkForceGauge(
        (const deviceDataGauge*)gauge->GetData(),
        gauge->m_byFieldId,
        (deviceDataGauge*)gaugeforce->GetData(),
        m_fOmega,
        m_bShiftCenter ? m_pfCxShift : m_pfCx,
        m_pDevicePath,
        2
    );

    hostPath[0] = 4; hostPath[1] = -2;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    OneLinkForceGauge(
        (const deviceDataGauge*)gauge->GetData(),
        gauge->m_byFieldId,
        (deviceDataGauge*)gaugeforce->GetData(),
        -m_fOmega,
        m_bShiftCenter ? m_pfCxShift : m_pfCx,
        m_pDevicePath,
        2
    );

    hostPath[0] = -2; hostPath[1] = 4;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    OneLinkForceGauge(
        (const deviceDataGauge*)gauge->GetData(),
        gauge->m_byFieldId,
        (deviceDataGauge*)gaugeforce->GetData(),
        -m_fOmega,
        m_bShiftCenter ? m_pfCxShift : m_pfCx,
        m_pDevicePath,
        2
    );

#pragma endregion

#pragma region -2 Re[Omega x y px py ]

    hostPath[0] = 1; hostPath[1] = 2;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    OneLinkForceGauge(
        (const deviceDataGauge*)gauge->GetData(),
        gauge->m_byFieldId,
        (deviceDataGauge*)gaugeforce->GetData(),
        m_fOmega * m_fOmega,
        m_bShiftCenter ? m_pfCxyShiftXYPP : m_pfCxy,
        m_pDevicePath,
        2
    );

    hostPath[0] = 2; hostPath[1] = 1;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    OneLinkForceGauge(
        (const deviceDataGauge*)gauge->GetData(),
        gauge->m_byFieldId,
        (deviceDataGauge*)gaugeforce->GetData(),
        m_fOmega * m_fOmega,
        m_bShiftCenter ? m_pfCxyShiftYXPP : m_pfCxy,
        m_pDevicePath,
        2
    );

    hostPath[0] = 2; hostPath[1] = -1;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    OneLinkForceGauge(
        (const deviceDataGauge*)gauge->GetData(),
        gauge->m_byFieldId,
        (deviceDataGauge*)gaugeforce->GetData(),
        -m_fOmega * m_fOmega,
        m_bShiftCenter ? m_pfCxyShiftYXPM : m_pfCxy,
        m_pDevicePath,
        2
    );

    hostPath[0] = -1; hostPath[1] = 2;
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    OneLinkForceGauge(
        (const deviceDataGauge*)gauge->GetData(),
        gauge->m_byFieldId,
        (deviceDataGauge*)gaugeforce->GetData(),
        -m_fOmega * m_fOmega,
        m_bShiftCenter ? m_pfCxyShiftXYMP : m_pfCxy,
        m_pDevicePath,
        2
    );

#pragma endregion

#pragma region +Omega^2 y^2 px^2 + Omega^2 x^2 py^2, shift term

    PartialSqForceGauge(
        (const deviceDataGauge*)gauge->GetData(),
        gauge->m_byFieldId,
        (deviceDataGauge*)gaugeforce->GetData(),
        m_fOmega * m_fOmega,
        m_bShiftCenter ? m_pfCySqShift : m_pfCySq,
        0
    );

    PartialSqForceGauge(
        (const deviceDataGauge*)gauge->GetData(),
        gauge->m_byFieldId,
        (deviceDataGauge*)gaugeforce->GetData(),
        m_fOmega* m_fOmega,
        m_bShiftCenter ? m_pfCxSqShift : m_pfCxSq,
        1
    );

#pragma endregion

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

    INT iValue = 0;
    if (param.FetchValueINT(_T("ShiftCenter"), iValue))
    {
        m_bShiftCenter = (0 != iValue);
    }
}

template<typename deviceDataBoson, typename deviceDataGauge>
CCString CFieldBosonVNRotation<deviceDataBoson, deviceDataGauge>::GetInfos(const CCString& tab) const
{
    CCString sRet = CFieldBosonVN<deviceDataBoson, deviceDataGauge>::GetInfos(tab);
    sRet = sRet + tab + _T("Omega : ") + appToString(m_fOmega) + _T("\n");
    sRet = sRet + tab + _T("ShiftCenter : ") + appToString(m_bShiftCenter) + _T("\n");
    return sRet;
}

__CLGIMPLEMENT_CLASS(CFieldBosonVNRotationU1)
__CLGIMPLEMENT_CLASS(CFieldBosonVNRotationSU2)
__CLGIMPLEMENT_CLASS(CFieldBosonVNRotationSU3)
__CLGIMPLEMENT_CLASS(CFieldBosonVNRotationSU4)
__CLGIMPLEMENT_CLASS(CFieldBosonVNRotationSU5)
__CLGIMPLEMENT_CLASS(CFieldBosonVNRotationSU6)
__CLGIMPLEMENT_CLASS(CFieldBosonVNRotationSU7)
__CLGIMPLEMENT_CLASS(CFieldBosonVNRotationSU8)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================