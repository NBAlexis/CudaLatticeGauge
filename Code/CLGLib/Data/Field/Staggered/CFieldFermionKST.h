//=============================================================================
// FILENAME : CFieldFermionKST.h
// 
// DESCRIPTION:
// This is the class for Kogut-Susskind staggered fermions
// For pseudo fermion, this is in fact a boson field phi.
//
// Current implementation, assumes square lattice
//
// REVISION:
//  [07/06/2024 nbale]
//=============================================================================
#include "Tools/Math/DeviceInlineTemplate.h"
#include "../Gauge/CFieldGaugeLink.h"
#include "CFieldFermionKSTKernel.h"

#ifndef _CFIELDFERMIONKST_H_
#define _CFIELDFERMIONKST_H_

#define _DEFINE_KS_FIELDSUN(n) \
__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSSU##n) \
class CLGAPI CFieldFermionKSSU##n : public CFieldFermionKST<deviceSU##n##Vector, deviceSU##n, n> \
{ \
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldFermionKSSU##n) \
public: \
    EFieldType GetFieldType() const override { return EFT_FermionStaggeredSU##n; } \
};

__BEGIN_NAMESPACE

template<typename deviceVector, typename deviceGauge, INT vectorN>
class __DLL_EXPORT CFieldFermionKST : public CFieldFermionKS
{
public:

    friend class CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>;

    CFieldFermionKST()
        : CFieldFermionKS()
        , m_pRationalFieldPointers(NULL)
    {
        CCommonKernel<deviceVector>::AllocateBuffer(&m_pDeviceData, m_uiSiteCount);
    }

    ~CFieldFermionKST()
    {
        CCommonKernel<deviceVector>::FreeBuffer(&m_pDeviceData);
        if (NULL != m_pRationalFieldPointers)
        {
            appSimpleFree(m_pRationalFieldPointers);
            m_pRationalFieldPointers = NULL;
        }
    }

    void InitialField(EFieldInitialType eInitialType) override
    {
        CCommonKernelSite<deviceVector>::InitialBuffer(m_pDeviceData, m_byFieldId, eInitialType);
    }

    void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFieldType) override
    {
        if (eFieldType != EFFT_CLGBin)
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN>::InitialFieldWithFile: Only support CLG Bin File\n"));
            return;
        }
        const UINT vectorn = _elementdim<deviceVector>();
        UINT uiSize = static_cast<UINT>(sizeof(Real) * vectorn * m_uiSiteCount);
        BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        if (NULL == data)
        {
            appCrucial(_T("File not found: %s\n"), sFileName.c_str());
            _FAIL_EXIT;
        }
        if (uiSize != (sizeof(Real) * vectorn * m_uiSiteCount))
        {
            appCrucial(_T("File size not correct: expecting: %d, found: %d\n"), static_cast<UINT>(sizeof(Real) * vectorn * m_uiSiteCount), uiSize);
            _FAIL_EXIT;
        }
        InitialWithByte(data);
        free(data);
    }

    void InitialWithByte(BYTE* byData) override
    {
        CCommonKernelField<deviceVector>::InitialWithByte(m_pDeviceData, m_uiSiteCount, byData);
    }

    void InitialOtherParameters(CParameters& params) override
    {
        CFieldFermionKS::InitialOtherParameters(params);

        if (NULL != m_pRationalFieldPointers)
        {
            appSimpleFree(m_pRationalFieldPointers);
        }
        appSimpleMalloc((void**)&m_pRationalFieldPointers, sizeof(deviceVector*) * 2 * m_rMD.m_uiDegree);
    }

    void DebugPrintMe() const override
    {
        CCommonKernelSite<deviceVector>::DebugPrint(m_pDeviceData, m_uiSiteCount);
    }

    void CopyTo(CField* U) const override
    {
        if (NULL == U || GetFieldType() != U->GetFieldType())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
            return;
        }

        CFieldFermionKS::CopyTo(U);

        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(U);

        CCommonKernel<deviceVector>::CopyBuffer(pField->m_pDeviceData, m_pDeviceData, m_uiSiteCount);

        if (NULL != pField->m_pRationalFieldPointers)
        {
            appSimpleFree(pField->m_pRationalFieldPointers);
        }
        appSimpleMalloc((void**)&pField->m_pRationalFieldPointers, sizeof(deviceVector*) * 2 * m_rMD.m_uiDegree);
    }

    void Dagger() override
    {
        CCommonKernelSite<deviceVector>::Dagger(m_pDeviceData, m_byFieldId);
    }

    //This is Axpy(1.0f, x)
    void AxpyPlus(const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
            return;
        }
        const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);

        CCommonKernelSite<deviceVector>::AxpyPlus(m_pDeviceData, m_byFieldId, pField->m_pDeviceData);
    }

    void AxpyMinus(const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
            return;
        }
        const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);

        CCommonKernelSite<deviceVector>::AxpyMinus(m_pDeviceData, m_byFieldId, pField->m_pDeviceData);
    }

    void Axpy(Real a, const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
            return;
        }
        const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);

        CCommonKernelSite<deviceVector>::Axpy(m_pDeviceData, m_byFieldId, a, pField->m_pDeviceData);
    }

    void Axpy(const CLGComplex& a, const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
            return;
        }
        const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);

        CCommonKernelSite<deviceVector>::Axpy(m_pDeviceData, m_byFieldId, a, pField->m_pDeviceData);
    }

    void Mul(const CField* other, UBOOL bDagger = TRUE) override
    {
        if (NULL == other || GetFieldType() != other->GetFieldType())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
            return;
        }
        const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(other);

        CCommonKernelSite<deviceVector>::Mul(m_pDeviceData, m_byFieldId, pField->m_pDeviceData, bDagger);
    }

    void ScalarMultply(const CLGComplex& a) override
    {
        CCommonKernelSite<deviceVector>::ScalarMultply(m_pDeviceData, m_byFieldId, a);
    }

    void ScalarMultply(Real a) override
    {
        CCommonKernelSite<deviceVector>::ScalarMultply(m_pDeviceData, m_byFieldId, a);
    }

    cuDoubleComplex Dot(const CField* x) const override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
            return make_cuDoubleComplex(0, 0);
        }
        const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);
        return CCommonKernelSite<deviceVector>::Dot(m_pDeviceData, m_byFieldId, pField->m_pDeviceData);
    }

protected:

    //pGauge must be gauge SU3
    //These are for Sparse linear algebra
    void DS(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override
    {
        const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
        if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
            return;
        }
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

        CCommonKernel<deviceVector>::CopyBuffer(pPooled->m_pDeviceData, m_pDeviceData, m_uiSiteCount);

        Real fRealCoeff = fCoeffReal;
        const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
        if (EOCT_Minus == eCoeffType)
        {
            eCoeffType = EOCT_Real;
            fRealCoeff = F(-1.0);
        }

        DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
            FALSE, eCoeffType, fRealCoeff, cCompCoeff);

        pPooled->Return();
    }

    void DdaggerS(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override
    {
        const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
        if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
            return;
        }
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
        CCommonKernel<deviceVector>::CopyBuffer(pPooled->m_pDeviceData, m_pDeviceData, m_uiSiteCount);

        Real fRealCoeff = fCoeffReal;
        const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
        if (EOCT_Minus == eCoeffType)
        {
            eCoeffType = EOCT_Real;
            fRealCoeff = F(-1.0);
        }

        DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
            TRUE, eCoeffType, fRealCoeff, cCompCoeff);


        pPooled->Return();
    }

    void DDS(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override
    {
        const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
        if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
            return;
        }

        Real fRealCoeff = fCoeffReal;
        const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
        if (EOCT_Minus == eCoeffType)
        {
            eCoeffType = EOCT_Real;
            fRealCoeff = F(-1.0);
        }
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

        DOperator(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
            FALSE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
        //why only apply coeff in the next step?
        DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
            FALSE, eCoeffType, fRealCoeff, cCompCoeff);

        pPooled->Return();
    }

    void DDdaggerS(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override
    {
        const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
        if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
            return;
        }

        Real fRealCoeff = fCoeffReal;
        const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
        if (EOCT_Minus == eCoeffType)
        {
            eCoeffType = EOCT_Real;
            fRealCoeff = F(-1.0);
        }
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

        DOperator(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
            TRUE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
        //why only apply coeff in the next step?
        DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
            FALSE, eCoeffType, fRealCoeff, cCompCoeff);

        pPooled->Return();
    }

    void DWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override
    {
        const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
        if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
            return;
        }
        if (m_bDiagonalMass)
        {
            appCrucial(_T("In the cass mass is not a number, should not in here!\n"));
        }

        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

        CCommonKernel<deviceVector>::CopyBuffer(pPooled->m_pDeviceData, m_pDeviceData, m_uiSiteCount);

        Real fRealCoeff = fCoeffReal;
        const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
        if (EOCT_Minus == eCoeffType)
        {
            eCoeffType = EOCT_Real;
            fRealCoeff = F(-1.0);
        }

        DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
            FALSE, eCoeffType, fRealCoeff, cCompCoeff);

        pPooled->Return();
    }

    void DdaggerWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override
    {
        const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
        if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
            return;
        }
        if (m_bDiagonalMass)
        {
            appCrucial(_T("In the cass mass is not a number, should not in here!\n"));
        }

        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
        CCommonKernel<deviceVector>::CopyBuffer(pPooled->m_pDeviceData, m_pDeviceData, m_uiSiteCount);

        Real fRealCoeff = fCoeffReal;
        const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
        if (EOCT_Minus == eCoeffType)
        {
            eCoeffType = EOCT_Real;
            fRealCoeff = F(-1.0);
        }

        DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
            TRUE, eCoeffType, fRealCoeff, cCompCoeff);


        pPooled->Return();
    }

    void DDWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override
    {
        const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
        if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
            return;
        }

        Real fRealCoeff = fCoeffReal;
        const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
        if (EOCT_Minus == eCoeffType)
        {
            eCoeffType = EOCT_Real;
            fRealCoeff = F(-1.0);
        }
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

        DOperatorKS(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
            FALSE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
        //why only apply coeff in the next step?
        DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
            FALSE, eCoeffType, fRealCoeff, cCompCoeff);

        pPooled->Return();
    }

    void DDdaggerWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override
    {
        const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
        if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
            return;
        }
        if (m_bDiagonalMass)
        {
            appCrucial(_T("In the cass mass is not a number, should not in here!\n"));
        }

        Real fRealCoeff = fCoeffReal;
        const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
        if (EOCT_Minus == eCoeffType)
        {
            eCoeffType = EOCT_Real;
            fRealCoeff = F(-1.0);
        }
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

        DOperatorKS(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
            TRUE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
        //why only apply coeff in the next step?
        DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
            FALSE, eCoeffType, fRealCoeff, cCompCoeff);

        pPooled->Return();
    }

public:

    void ApplyGamma(EGammaMatrix eGamma) override
    {
        appCrucial(_T("Not supported yet...\n"));
    }

    TArray<CFieldFermion*> GetSourcesAtSiteFromPool(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const SSmallInt4& site) const override
    {
        TArray<CFieldFermion*> ret;
        for (UINT j = 0; j < _dim<deviceVector>(); ++j)
        {
            ret.AddItem(dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId)));
            if (NULL == ret[j])
            {
                appCrucial(_T("GetSourcesAtSiteFromPool failed!\n"));
                _FAIL_EXIT;
            }
        }

        for (BYTE c = 0; c < _dim<deviceVector>(); ++c)
        {
            SFermionBosonSource sourceData;
            sourceData.m_eSourceType = EFS_Point;
            sourceData.m_sSourcePoint = site;
            sourceData.m_byColorIndex = c;
            sourceData.m_bySpinIndex = 0;

            ret[c]->InitialAsSource(sourceData);

            if (NULL != appGetFermionSolver(m_byFieldId) && !appGetFermionSolver(m_byFieldId)->IsAbsoluteAccuracy())
            {
                ret[c]->m_fLength = ret[c]->Dot(ret[c]).x;
            }

            ret[c]->InverseD(gaugeNum, bosonNum, gaugeFields, pBoson);
        }
        return ret;
    }

    /**
    * generate phi by gaussian random.
    * phi = (D^+D)^{1/4} phi
    */
    void PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson) override
    {
        CCommonKernelSite<deviceVector>::InitialBuffer(m_pDeviceData, m_byFieldId, EFIT_RandomGaussian);

        D_MC(gaugeNum, bosonNum, gaugeFields, pBoson);

        if (NULL != appGetFermionSolver(m_byFieldId) && !appGetFermionSolver(m_byFieldId)->IsAbsoluteAccuracy())
        {
            m_fLength = Dot(this).x;
        }

        //For KS, we generally use shifted solver, so do NOT cache the last result
    }

protected:

    void ApplyGammaKSS(const CFieldGauge* pGauge, EGammaMatrix eGamma) override
    {
        if (NULL == pGauge || vectorN != pGauge->MatrixN())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
            return;
        }
        const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
        CCommonKernel<deviceVector>::CopyBuffer(pPooled->m_pDeviceData, m_pDeviceData, m_uiSiteCount);
        InitialField(EFIT_Zero);

        //If it was gamma_mu or gamma_5 or sigmaij, it is i gamma mu and i gamma 5, therefore multiply -i
        UBOOL bImag = (GAMMA1 == eGamma) || (GAMMA2 == eGamma) || (GAMMA3 == eGamma) || (GAMMA4 == eGamma) || (GAMMA5 == eGamma)
            || (SIGMA12 == eGamma) || (SIGMA31 == eGamma) || (SIGMA41 == eGamma) || (SIGMA23 == eGamma) || (SIGMA42 == eGamma) || (SIGMA43 == eGamma);

        CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::appApplyGammaKS(
            m_pDeviceData,
            pPooled->m_pDeviceData,
            pFieldSU3->m_pDeviceData,
            eGamma,
            m_bEachSiteEta,
            FALSE,
            F(0.5),
            bImag ? EOCT_Complex : EOCT_None,
            F(1.0),
            bImag ? _make_cuComplex(F(0.0), -F(1.0)) : _onec,
            m_byFieldId,
            pGauge->m_byFieldId
        );

        pPooled->Return();
    }

    //================= test anti-hermitian =========
    UINT TestAntiHermitianS(const CFieldGauge* pGauge) const override
    {
        return CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::TestAntiHermitianS(m_byFieldId, pGauge);
    }

    //These are truely D or InverseD etc.

    /**
     * Use to calculate action, it is (D^+D)^{-1/4}
     */
    void D0S(const CField* pGauge) override
    {
        const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
        if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
            return;
        }
        if (m_bDiagonalMass)
        {
            appCrucial(_T("In the cass mass is not a number, should not in here except for check anti-hermiticity!\n"));
        }

        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

        CCommonKernel<deviceVector>::CopyBuffer(pPooled->m_pDeviceData, m_pDeviceData, m_uiSiteCount);

        DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, F(0.0),
            FALSE, EOCT_None, F(1.0), _onec);

        pPooled->Return();
    }

    UBOOL CalculateForceS(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const override
    {
        if (NULL == pGauge || vectorN != pGauge->MatrixN())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
            return FALSE;
        }
        if (NULL == pForce || vectorN != pForce->MatrixN())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
            return FALSE;
        }

        TArray<CField*> phii;
        TArray<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*> phiid;
        for (UINT i = 0; i < m_rMD.m_uiDegree; ++i)
        {
            CField* pPhi_i = dynamic_cast<CField*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
            phii.AddItem(pPhi_i);
            CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPhi_id = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
            phiid.AddItem(pPhi_id);
        }

        CMultiShiftSolver* solver = appGetMultiShiftSolver(m_byFieldId);
        if (NULL == solver)
        {
            appCrucial(_T("No multi solver found!"));
            _FAIL_EXIT;
        }

        TArray<CLGComplex> shifts;
        for (UINT i = 0; i < m_rMD.m_uiDegree; ++i)
        {
            shifts.AddItem(_make_cuComplex(m_rMD.m_lstB[i], F(0.0)));
        }
        TArray<const CFieldGauge*> gaues;
        gaues.AddItem(pGauge);
        solver->Solve(phii, shifts, this, 1, 0, gaues.GetData(), NULL, EFO_F_DDdagger);

        const UINT uiBufferSize = sizeof(deviceVector*) * 2 * m_rMD.m_uiDegree;
        deviceVector** hostPointers = (deviceVector**)appAlloca(uiBufferSize);
        for (UINT i = 0; i < m_rMD.m_uiDegree; ++i)
        {
            CFieldFermionKST<deviceVector, deviceGauge, vectorN>* phi_ks = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(phii[i]);
            phi_ks->CopyTo(phiid[i]);
            if (m_bDiagonalMass)
            {
                phiid[i]->DS(pGauge);
            }
            else
            {
                phiid[i]->D0S(pGauge);
            }

            hostPointers[i] = phi_ks->m_pDeviceData;
            hostPointers[i + m_rMD.m_uiDegree] = phiid[i]->m_pDeviceData;
        }

        appSimpleCopyHD(m_pRationalFieldPointers, hostPointers, uiBufferSize);

        const CFieldGaugeLink<deviceGauge, vectorN>* pGaugeSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
        CFieldGaugeLink<deviceGauge, vectorN>* pForceSU3 = dynamic_cast<CFieldGaugeLink<deviceGauge, vectorN>*>(pForce);

        DerivateD0(pForceSU3->m_pDeviceData, pGaugeSU3->m_pDeviceData, pGaugeSU3->m_byFieldId);

        for (UINT i = 0; i < m_rMD.m_uiDegree; ++i)
        {
            phii[i]->Return();
            phiid[i]->Return();
        }

        return TRUE;
    }

public:

    //For test only
    void PrepareForHMCOnlyRandomize() override
    {
        CCommonKernelSite<deviceVector>::InitialBuffer(m_pDeviceData, m_byFieldId, EFIT_RandomGaussian);
    }

    void PrepareForHMCNotRandomize(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson) override
    {
        D_MC(gaugeNum, bosonNum, gaugeFields, pBoson);
    }

    void InitialAsSource(const SFermionBosonSource& sourceData) override
    {
        CCommonKernelSite<deviceVector>::InitialSource(m_pDeviceData, m_byFieldId, sourceData);
    }

    BYTE* CopyDataOut(UINT& uiSize) const override
    {
        return CCommonKernelField<deviceVector>::CopyDataOut(m_pDeviceData, m_uiSiteCount, uiSize);
    }

    BYTE* CopyDataOutFloat(UINT& uiSize) const override
    {
        return CCommonKernelField<deviceVector>::CopyDataOutFloat(m_pDeviceData, m_uiSiteCount, uiSize);
    }

    BYTE* CopyDataOutDouble(UINT& uiSize) const override
    {
        return CCommonKernelField<deviceVector>::CopyDataOutDouble(m_pDeviceData, m_uiSiteCount, uiSize);
    }

protected:

    //============================
    //Override these two functions for KS
    void DerivateD0(void* pForce, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override
    {
        CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DerivateD0(m_pDeviceData, m_byFieldId, (deviceGauge*)pForce, (const deviceGauge*)pGaugeBuffer, byGaugeFieldId, m_pRationalFieldPointers, m_pMDNumerator, m_rMD.m_uiDegree);
    }

    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override
    {
        CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DOperatorKS(m_bEachSiteEta, (deviceVector*)pTargetBuffer, (const deviceVector*)pBuffer, (const deviceGauge*)pGaugeBuffer, m_byFieldId, byGaugeFieldId, f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);
    }

public:

    #pragma region Help functions to implement higher orders

    void OnlyMass(void* pTarget, Real f2am, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override
    {
        CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::OnlyMass(m_pDeviceData, (deviceVector*)pTarget, f2am, eOCT, fRealCoeff, cCmpCoeff);
    }

    #pragma endregion

    deviceVector* m_pDeviceData;

    _GetData

protected:

    void OneLinkS(const void* pGauge, BYTE byGaugeFieldId, void* pTarget, Real fCoefficient,
        const INT* pDevicePath, BYTE pathLength, BYTE byEtaIdx,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override
    {
        CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::OneLinkS(m_pDeviceData, m_byFieldId, (const deviceGauge*)pGauge, byGaugeFieldId, (deviceVector*)pTarget, fCoefficient, pDevicePath, pathLength, byEtaIdx, bDagger, eOCT, fRealCoeff, cCmpCoeff);
    }

    void OneLinkForceS(const void* pGauge, BYTE byGaugeFieldId, void* pForce, Real fCoefficient,
        const INT* pDevicePath, BYTE pathLength, BYTE byEtaIdx) const override
    {
        CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::OneLinkForceS(m_pDeviceData, m_byFieldId, (const deviceGauge*)pGauge, byGaugeFieldId, (deviceGauge*)pForce, fCoefficient, pDevicePath, pathLength, byEtaIdx, m_pRationalFieldPointers, m_pMDNumerator, m_rMD.m_uiDegree);
    }


    //phi _i and Dst0 phi _i
    deviceVector** m_pRationalFieldPointers;

};

template<typename deviceVector, typename deviceGauge, INT vectorN>
class __DLL_EXPORT CFieldMatrixOperationKST : public CFieldMatrixOperation
{
public:
    CFieldMatrixOperationKST()
        : CFieldMatrixOperation()
        , m_pResBuffer(NULL)
        , m_pLeftBuffer(NULL)
        , m_pHostResBuffer(NULL)
        , m_pHostLeftBuffer(NULL)
    {
        m_pHostResBuffer = (deviceVector**)malloc(sizeof(deviceVector*) * _kFieldMatrixMaxDim);
        m_pHostLeftBuffer = (deviceVector**)malloc(sizeof(deviceVector*) * _kFieldMatrixMaxDim);
        appSimpleMalloc((void**)&m_pResBuffer, sizeof(deviceVector*) * _kFieldMatrixMaxDim);
        appSimpleMalloc((void**)&m_pLeftBuffer, sizeof(deviceVector*) * _kFieldMatrixMaxDim);
    }
    ~CFieldMatrixOperationKST()
    {
        free(m_pHostResBuffer);
        free(m_pHostLeftBuffer);
        appSimpleFree(m_pResBuffer);
        appSimpleFree(m_pLeftBuffer);
    }

    //real left = (res,left)
    void VectorMultiplyMatrix(TArray<CField*>& res, const TArray<CField*>& left, const CLGComplex* deviceMatrix, UINT uiDimX, UINT uiDimY) override
    {
        CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::VectorMultiplyMatrix(m_pHostResBuffer, m_pHostLeftBuffer, m_pResBuffer, m_pLeftBuffer, res, left, deviceMatrix, uiDimX, uiDimY);
    }

    deviceVector** m_pResBuffer;
    deviceVector** m_pLeftBuffer;
    deviceVector** m_pHostResBuffer;
    deviceVector** m_pHostLeftBuffer;
};

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSU1)
class CLGAPI CFieldFermionKSU1 : public CFieldFermionKST<CLGComplex, CLGComplex, 1>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldFermionKSU1)
public:
    EFieldType GetFieldType() const override { return EFT_FermionStaggeredU1; }
};

_DEFINE_KS_FIELDSUN(2)
_DEFINE_KS_FIELDSUN(4)
_DEFINE_KS_FIELDSUN(5)
_DEFINE_KS_FIELDSUN(6)
_DEFINE_KS_FIELDSUN(7)
_DEFINE_KS_FIELDSUN(8)

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKST_H_

//=============================================================================
// END OF FILE
//=============================================================================