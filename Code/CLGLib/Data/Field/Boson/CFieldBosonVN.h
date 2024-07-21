//=============================================================================
// FILENAME : CFieldBosonVN.h
// 
// DESCRIPTION:
// This is the class for all boson fields
//
// REVISION:
//  [3/31/2024 nbale]
//=============================================================================
#include "Tools/Math/DeviceTemplates/DeviceInlineUseNoTemplateFunction.h"
#include "CFieldBosonVNKernel.h"

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
class __DLL_EXPORT CFieldBosonVN : public CFieldBoson
{
public:
    CFieldBosonVN() 
        : CFieldBoson()
        , m_pDeviceData(NULL)
    {
        CCommonKernel<deviceDataBoson>::AllocateBuffer(&m_pDeviceData, m_uiSiteCount);
    }

    ~CFieldBosonVN()
    {
        CCommonKernel<deviceDataBoson>::FreeBuffer(&m_pDeviceData);
    }

    void CopyTo(CField* U) const override
    {
        if (NULL == U || GetFieldType() != U->GetFieldType())
        {
            appCrucial(_T("EFT_BosonU1 can only copy to EFT_BosonU1!"));
            return;
        }

        CFieldBoson::CopyTo(U);

        CFieldBosonVN<deviceDataBoson, deviceDataGauge>* pField = dynamic_cast<CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(U);
        CCommonKernel<deviceDataBoson>::CopyBuffer(pField->m_pDeviceData, m_pDeviceData, m_uiSiteCount);
    }

    /**
    * This should be momentum field
    */
    void MakeRandomMomentum() override
    {
        if (m_bConstant)
        {
            Zero();
            return;
        }

        CCommonKernelSite<deviceDataBoson>::InitialBuffer(m_pDeviceData, m_byFieldId, EFIT_RandomGaussian);
    }

    void ForceOnGauge(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, CFieldGauge* const* pGaugeForce, const CFieldBoson* const* pBoson) const override
    {
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

        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::ForceOnGauge(m_pDeviceData, m_byFieldId, gauge->m_byFieldId, (const deviceDataGauge*)gauge->GetData(), (deviceDataGauge*)gaugeforce->GetData());
    }

    void InitialField(EFieldInitialType eInitialType) override
    {
        CCommonKernelSite<deviceDataBoson>::InitialBuffer(m_pDeviceData, m_byFieldId, eInitialType);
    }

    void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFieldType) override
    {
        if (eFieldType != EFFT_CLGBin
         && eFieldType != EFFT_CLGBinFloat
         && eFieldType != EFFT_CLGBinDouble)
        {
            appCrucial(_T("CFieldBosonU1::InitialFieldWithFile: Not support %s File\n"), __ENUM_TO_STRING(EFieldFileType, eFieldType).c_str());
            return;
        }

        UINT uiSize = static_cast<UINT>(sizeof(Real) * FloatN() * m_uiSiteCount);
        if (eFieldType == EFFT_CLGBinFloat)
        {
            uiSize = static_cast<UINT>(sizeof(FLOAT) * FloatN() * m_uiSiteCount);
        }
        else if (eFieldType == EFFT_CLGBinDouble)
        {
            uiSize = static_cast<UINT>(sizeof(DOUBLE) * FloatN() * m_uiSiteCount);
        }
        UINT uiReadSize = uiSize;
        BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiReadSize);
        if (NULL == data)
        {
            appCrucial(_T("File not found: %s\n"), sFileName.c_str());
            _FAIL_EXIT;
        }

        if (uiSize != uiReadSize)
        {
            appCrucial(_T("File size not correct: expecting: %d, found: %d\n"), uiSize, uiReadSize);
            _FAIL_EXIT;
        }

#if _CLG_DOUBLEFLOAT
        if (eFieldType == EFFT_CLGBinFloat)
        {
            FLOAT* data1 = (FLOAT*)data;
            DOUBLE* data2 = (DOUBLE*)(malloc(sizeof(DOUBLE) * FloatN() * m_uiSiteCount));
            for (UINT i = 0; i < m_uiSiteCount * FloatN(); ++i)
            {
                data2[i] = static_cast<DOUBLE>(data1[i]);
            }
            InitialWithByte((BYTE*)data2);
            free(data);
            free(data2);
        }
#else
        if (eFieldType == EFFT_CLGBinDouble)
        {
            DOUBLE* data1 = (DOUBLE*)data;
            FLOAT* data2 = (FLOAT*)(malloc(sizeof(FLOAT) * FloatN() * m_uiSiteCount));
            for (UINT i = 0; i < m_uiSiteCount * FloatN(); ++i)
            {
                data2[i] = static_cast<FLOAT>(data1[i]);
            }
            InitialWithByte((BYTE*)data2);
            free(data);
            free(data2);
        }
#endif
        else
        {
            InitialWithByte(data);
            free(data);
        }

        appCrucial(_T("CFieldBosonU1::InitialFieldWithFile: Not support %s File\n"), __ENUM_TO_STRING(EFieldFileType, eFieldType).c_str());
    }

    void InitialWithByte(BYTE* byData) override
    {
        CCommonKernelField<deviceDataBoson>::InitialWithByte(m_pDeviceData, m_uiSiteCount, byData);
    }

    //void InitialOtherParameters(CParameters& params) override;
    void DebugPrintMe() const override
    {
        CCommonKernelSite<deviceDataBoson>::DebugPrint(m_pDeviceData, m_uiSiteCount);
    }

    void Dagger() override
    {
        CCommonKernelSite<deviceDataBoson>::Dagger(m_pDeviceData, m_byFieldId);
    }

    //This is Axpy(1.0f, x)
    void FixBoundary() override
    {
        CCommonKernelSite<deviceDataBoson>::FixBoundary(m_pDeviceData, m_byFieldId);
    }

    void AxpyPlus(const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonVN!"));
            return;
        }
        const CFieldBosonVN<deviceDataBoson, deviceDataGauge>* pField = dynamic_cast<const CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(x);
        CCommonKernelSite<deviceDataBoson>::AxpyPlus(m_pDeviceData, m_byFieldId, pField->m_pDeviceData);
    }

    void AxpyMinus(const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonVN!"));
            return;
        }
        const CFieldBosonVN<deviceDataBoson, deviceDataGauge>* pField = dynamic_cast<const CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(x);
        CCommonKernelSite<deviceDataBoson>::AxpyMinus(m_pDeviceData, m_byFieldId, pField->m_pDeviceData);
    }

    void Axpy(Real a, const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonVN!"));
            return;
        }
        const CFieldBosonVN<deviceDataBoson, deviceDataGauge>* pField = dynamic_cast<const CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(x);
        CCommonKernelSite<deviceDataBoson>::Axpy(m_pDeviceData, m_byFieldId, a, pField->m_pDeviceData);
    }

    void Axpy(const CLGComplex& a, const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonVN!"));
            return;
        }
        const CFieldBosonVN<deviceDataBoson, deviceDataGauge>* pField = dynamic_cast<const CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(x);
        CCommonKernelSite<deviceDataBoson>::Axpy(m_pDeviceData, m_byFieldId, a, pField->m_pDeviceData);
    }

    void Mul(const CField* x, UBOOL bDagger = TRUE) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonVN!"));
            return;
        }
        const CFieldBosonVN<deviceDataBoson, deviceDataGauge>* pField = dynamic_cast<const CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(x);
        CCommonKernelSite<deviceDataBoson>::Mul(m_pDeviceData, m_byFieldId, pField->m_pDeviceData, bDagger);
    }

    void ScalarMultply(const CLGComplex& a) override
    {
        CCommonKernelSite<deviceDataBoson>::ScalarMultply(m_pDeviceData, m_byFieldId, a);
    }

    void ScalarMultply(Real a) override
    {
        CCommonKernelSite<deviceDataBoson>::ScalarMultply(m_pDeviceData, m_byFieldId, a);
    }

    cuDoubleComplex Dot(const CField* x) const override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonVN!"));
            return make_cuDoubleComplex(0.0, 0.0);
        }
        const CFieldBosonVN<deviceDataBoson, deviceDataGauge>* pField = dynamic_cast<const CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(x);
        return CCommonKernelSite<deviceDataBoson>::Dot(m_pDeviceData, m_byFieldId, pField->m_pDeviceData);
    }

    TArray<DOUBLE> Sum() const override
    {
        return CCommonKernelSite<deviceDataBoson>::Sum(m_pDeviceData, m_byFieldId);
    }

    BYTE* CopyDataOut(UINT& uiSize) const override
    {
        return CCommonKernelField<deviceDataBoson>::CopyDataOut(m_pDeviceData, m_uiSiteCount, uiSize);
    }

    BYTE* CopyDataOutFloat(UINT& uiSize) const override
    {
        return CCommonKernelField<deviceDataBoson>::CopyDataOutFloat(m_pDeviceData, m_uiSiteCount, uiSize);
    }

    BYTE* CopyDataOutDouble(UINT& uiSize) const override
    {
        return CCommonKernelField<deviceDataBoson>::CopyDataOutDouble(m_pDeviceData, m_uiSiteCount, uiSize);
    }

    UINT CheckHermitian(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson) const override
    {
        return CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CheckHermitian(this, m_uiSiteCount, gaugeNum, bosonNum, gaugeFields, pBoson);
    }

    void InitialAsSource(const SFermionBosonSource& sourceData) override
    {
        CCommonKernelSite<deviceDataBoson>::InitialSource(m_pDeviceData, m_byFieldId, sourceData);
    }

protected:

    void DFromSource(const CFieldBoson* pSource, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override
    {
        const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);

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

        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::DFromSource(pSourceVN->m_pDeviceData, m_pDeviceData, m_byFieldId, 
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
        const INT* pDevicePath, BYTE pathLength, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff)
    {
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::OneLink(m_pDeviceData, m_byFieldId, pSource, pGauge, byGaugeFieldId, fCoeffiecient, fpCoeff, pDevicePath, pathLength, eOCT, fRealCoeff, cCmpCoeff);
    }

    /**
    * f(x) phi^*(x)phi(x)
    */
    virtual void DiagnalTerm(const deviceDataBoson* pSource, DOUBLE fCoeffiecient, _deviceCoeffFunctionPointer fpCoeff, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff)
    {
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::DiagnalTerm(m_pDeviceData, m_byFieldId, pSource, fCoeffiecient, fpCoeff, eOCT, fRealCoeff, cCmpCoeff);
    }

    virtual void OneLinkForceGauge(const deviceDataGauge* pGuage, BYTE byGaugeFieldId, deviceDataGauge* pForce, DOUBLE fCoeffiecient, _deviceCoeffFunctionPointerTwoSites fpCoeff, const INT* pDevicePath, BYTE pathLength) const
    {
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::OneLinkForceGauge(m_pDeviceData, m_byFieldId, pGuage, byGaugeFieldId, pForce, fCoeffiecient, fpCoeff, pDevicePath, pathLength);
    }

    /**
    * off-diagnal terms of (partial _x phi)^2
    * f(n) U_mu(n)phi(n+mu) + f(n-mu) U_{-mu}(n)phi(n-mu)
    */
    virtual void PartialSq(const deviceDataBoson* pSource, const deviceDataGauge* pGuage, BYTE byGaugeFieldId,
        DOUBLE fCoeffiecient, _deviceCoeffFunctionPointer fpCoeff, BYTE idir,
        EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff)
    {
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::PartialSq(m_pDeviceData, m_byFieldId, pSource, pGuage, byGaugeFieldId, fCoeffiecient, fpCoeff, idir, eOCT, fRealCoeff, cCmpCoeff);
    }

    virtual void PartialSqForceGauge(const deviceDataGauge* pGuage, BYTE byGaugeFieldId, deviceDataGauge* pForce, DOUBLE fCoeffiecient, _deviceCoeffFunctionPointer fpCoeff, BYTE idir) const
    {
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::PartialSqForceGauge(m_pDeviceData, m_byFieldId, pGuage, byGaugeFieldId, pForce, fCoeffiecient, fpCoeff, idir);
    }

public:

    deviceDataBoson* m_pDeviceData;

    _GetData

};

__DEFINE_BOSON_FIELD(CFieldBosonU1, CLGComplex, CLGComplex, 1, 2, EFT_BosonComplex)
__DEFINE_BOSON_FIELD(CFieldBosonSU2, deviceSU2Vector, deviceSU2, 2, 4, EFT_BosonComplexVector2)
__DEFINE_BOSON_FIELD(CFieldBosonSU3, deviceSU3Vector, deviceSU3, 3, 6, EFT_BosonComplexVector3)
__DEFINE_BOSON_FIELD(CFieldBosonSU4, deviceSU4Vector, deviceSU4, 4, 8, EFT_BosonComplexVector4)
__DEFINE_BOSON_FIELD(CFieldBosonSU5, deviceSU5Vector, deviceSU5, 5, 10, EFT_BosonComplexVector5)
__DEFINE_BOSON_FIELD(CFieldBosonSU6, deviceSU6Vector, deviceSU6, 6, 12, EFT_BosonComplexVector6)
__DEFINE_BOSON_FIELD(CFieldBosonSU7, deviceSU7Vector, deviceSU7, 7, 14, EFT_BosonComplexVector7)
__DEFINE_BOSON_FIELD(CFieldBosonSU8, deviceSU8Vector, deviceSU8, 8, 16, EFT_BosonComplexVector8)

__END_NAMESPACE

#endif //#ifndef _CFIELDBOSONVN_H_

//=============================================================================
// END OF FILE
//=============================================================================