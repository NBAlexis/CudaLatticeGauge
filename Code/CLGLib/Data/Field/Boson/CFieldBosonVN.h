//=============================================================================
// FILENAME : CFieldBosonVN.h
// 
// DESCRIPTION:
// This is the class for all boson fields
//
// REVISION:
//  [mm/dd/yy]
//  [3/31/2024 nbale]
//=============================================================================
#pragma once

#include "Tools/Math/DeviceTemplates/DeviceInlineUseNoTemplateFunction.h"
#include "CFieldBosonVNKernel.h"
#include "CFieldBosonT.h"

#ifndef _CFIELDBOSONVN_H_
#define _CFIELDBOSONVN_H_

#define __DEFINE_BOSON_FIELD(FIELD_NAME, TYPE_BOSON, TYPE_GAUGE, VECTOR_N, FLOAT_N, ELEMENT_TYPE) \
__CLG_REGISTER_HELPER_HEADER(FIELD_NAME) \
class CLGAPI FIELD_NAME : public CFieldBosonVN<TYPE_BOSON, TYPE_GAUGE> \
{ \
    __CLGDECLARE_FIELDWITHOUTCOPYTO(FIELD_NAME) \
public: \
    EFieldType GetFieldType() const override { return ELEMENT_TYPE; } \
    UINT VectorN() const override { return VECTOR_N; } \
    UINT FloatN() const override { return FLOAT_N; } \
};


__BEGIN_NAMESPACE

template<typename deviceDataBoson, typename deviceDataGauge>
class __DLL_EXPORT CFieldBosonVN : public CFieldBosonT<deviceDataBoson>
{
public:

    void ForceOnGauge(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, CFieldGauge* const* pGaugeForce, const CFieldBoson* const* pBoson) const override
    {
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

        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::ForceOnGauge(this->m_pDeviceData, this->m_byFieldId, gauge->m_byFieldId, (const deviceDataGauge*)gauge->GetData(), (deviceDataGauge*)gaugeforce->GetData());
    }

    UINT CheckHermitian(INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields) const override
    {
        return CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CheckHermitian(this, this->m_uiSiteCount, gaugeNum, bosonNum, tensor2Num, gaugeFields, pBoson, tensor2Fields);
    }

protected:

    void DFromSource(const CFieldBoson* pSource, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override
    {
        const CFieldGauge* pGauge = this->GetDefaultGauge(gaugeNum, gaugeFields);

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

        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::DFromSource(pSourceVN->m_pDeviceData, this->m_pDeviceData, this->m_byFieldId,
            NULL == pGauge ? 0 : pGauge->m_byFieldId, 
            NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
            eCoeffType,
            fRealCoeff,
            cCompCoeff);
    }

    /**
    * NOTE: If put to D operator, this is minus
    * S = - (Dphi)^2 + c1 phi^2 + c2 phi^4 + ...
    * This is not D, but the minus terms in D
    */
    virtual void OneLink(const deviceDataBoson* pSource, const deviceDataGauge* pGauge, BYTE byGaugeFieldId,
        DOUBLE fCoeffiecient, _deviceCoeffFunctionPointerTwoSites fpCoeff,
        const SCHAR* pDevicePath, BYTE pathLength, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff)
    {
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::OneLink(this->m_pDeviceData, this->m_byFieldId, pSource, pGauge, byGaugeFieldId, fCoeffiecient, fpCoeff, pDevicePath, pathLength, eOCT, fRealCoeff, cCmpCoeff);
    }

    virtual void OneLinkForceGauge(const deviceDataGauge* pGuage, BYTE byGaugeFieldId, deviceDataGauge* pForce, DOUBLE fCoeffiecient, _deviceCoeffFunctionPointerTwoSites fpCoeff, const SCHAR* pDevicePath, BYTE pathLength) const
    {
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::OneLinkForceGauge(this->m_pDeviceData, this->m_byFieldId, pGuage, byGaugeFieldId, pForce, fCoeffiecient, fpCoeff, pDevicePath, pathLength);
    }

    /**
    * off-diagnal terms of (partial _x phi)^2
    * f(n) U_mu(n)phi(n+mu) + f(n-mu) U_{-mu}(n)phi(n-mu)
    */
    virtual void PartialSq(const deviceDataBoson* pSource, const deviceDataGauge* pGuage, BYTE byGaugeFieldId,
        DOUBLE fCoeffiecient, _deviceCoeffFunctionPointer fpCoeff, BYTE idir,
        EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff)
    {
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::PartialSq(this->m_pDeviceData, this->m_byFieldId, pSource, pGuage, byGaugeFieldId, fCoeffiecient, fpCoeff, idir, eOCT, fRealCoeff, cCmpCoeff);
    }

    virtual void PartialSqForceGauge(const deviceDataGauge* pGuage, BYTE byGaugeFieldId, deviceDataGauge* pForce, DOUBLE fCoeffiecient, _deviceCoeffFunctionPointer fpCoeff, BYTE idir) const
    {
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::PartialSqForceGauge(this->m_pDeviceData, this->m_byFieldId, pGuage, byGaugeFieldId, pForce, fCoeffiecient, fpCoeff, idir);
    }

};

__DEFINE_BOSON_FIELD(CFieldBosonU1, CLGComplex, CLGComplex, 1, 2, EFT_BosonComplex)
__DEFINE_BOSON_FIELD(CFieldBosonSU2, deviceSU2Vector, deviceSU2, 2, 4, EFT_BosonComplexVector2)
__DEFINE_BOSON_FIELD(CFieldBosonSU3, deviceSU3Vector, deviceSU3, 3, 6, EFT_BosonComplexVector3)

#if _CLG_SU4_BOSON
__DEFINE_BOSON_FIELD(CFieldBosonSU4, deviceSU4Vector, deviceSU4, 4, 8, EFT_BosonComplexVector4)
#endif
#if _CLG_SU5_BOSON
__DEFINE_BOSON_FIELD(CFieldBosonSU5, deviceSU5Vector, deviceSU5, 5, 10, EFT_BosonComplexVector5)
#endif
#if _CLG_SU6_BOSON
__DEFINE_BOSON_FIELD(CFieldBosonSU6, deviceSU6Vector, deviceSU6, 6, 12, EFT_BosonComplexVector6)
#endif
#if _CLG_SU7_BOSON
__DEFINE_BOSON_FIELD(CFieldBosonSU7, deviceSU7Vector, deviceSU7, 7, 14, EFT_BosonComplexVector7)
#endif
#if _CLG_SU8_BOSON
__DEFINE_BOSON_FIELD(CFieldBosonSU8, deviceSU8Vector, deviceSU8, 8, 16, EFT_BosonComplexVector8)
#endif

__END_NAMESPACE

#endif //#ifndef _CFIELDBOSONVN_H_

//=============================================================================
// END OF FILE
//=============================================================================