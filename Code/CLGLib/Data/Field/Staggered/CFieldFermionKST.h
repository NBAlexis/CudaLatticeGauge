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
//  [mm/dd/yy]
//  [07/06/2024 nbale]
//=============================================================================
#include "Tools/Math/DeviceInlineTemplate.h"
#include "../Gauge/CFieldGaugeLink.h"
#include "CFieldFermionKSTKernel.h"
#include "CFieldFermionKSTKernelR.h"
#include "CFieldFermionKSTKernelGamma.h"

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

class CLGAPI CRationalFieldPointer : public CRegisteredBufferCache
{
public:

    CRationalFieldPointer()
        : CRegisteredBufferCache()
        , m_pRationalFieldPointers(NULL)
        , m_uiRationFieldPointerBufferSize(0)
    {

    }
    ~CRationalFieldPointer()
    {
        if (NULL != m_pRationalFieldPointers)
        {
            //Improve-1 (I6): deregister the set before the outer extent dies.
            m_HaloBufferSet.Unbind();
            checkCudaErrors(__cudaFree(m_pRationalFieldPointers));
            m_pRationalFieldPointers = NULL;
            m_uiRationFieldPointerBufferSize = 0;
        }
        // The instance is managed by CLGManager as a registered cache and is
        // deleted on appQuitCLG().  Reset the singleton pointer so the next
        // process/test creates a fresh instance.
        m_pPointer = NULL;
    }

    template<typename T> 
    T** GetRationPoint(BYTE count)
    {
        UINT size = sizeof(T*) * count;
        if (NULL != m_pRationalFieldPointers && m_uiRationFieldPointerBufferSize < size)
        {
            //Improve-1 (I6): the outer extent is about to move; drop the set
            //first (the owner re-binds after the next fill).
            m_HaloBufferSet.Unbind();
            checkCudaErrors(__cudaFree(m_pRationalFieldPointers));
            checkCudaErrors(__cudaMalloc((void**)&m_pRationalFieldPointers, size));
            m_uiRationFieldPointerBufferSize = size;
        }
        else if (NULL == m_pRationalFieldPointers)
        {
            checkCudaErrors(__cudaMalloc((void**)&m_pRationalFieldPointers, size));
            m_uiRationFieldPointerBufferSize = size;
        }
        return reinterpret_cast<T**>(m_pRationalFieldPointers);
    }

    /**
    * Improve-1 (multi-GPU-improve1.md 3.4.1, I6): the device pointer array is
    * re-filled with (possibly different) pooled member fields on every force
    * evaluation. Re-bind the buffer-set handle so a launch argument hitting
    * the array expands to the CURRENT member handles. The owner calls this
    * right after each host-to-device fill of the array.
    */
    void RebindMemberHandles(const TArray<const CHaloBufferHandle*>& members)
    {
        if (m_HaloBufferSet.IsBound())
        {
            m_HaloBufferSet.Unbind();
        }
        if (NULL != m_pRationalFieldPointers && members.Num() > 0)
        {
            m_HaloBufferSet.Bind(reinterpret_cast<const BYTE*>(m_pRationalFieldPointers),
                static_cast<size_t>(m_uiRationFieldPointerBufferSize), members);
        }
    }

    void** m_pRationalFieldPointers;
    UINT m_uiRationFieldPointerBufferSize;

    //Improve-1 (I6): buffer-set handle over the pointer array (3.4.1).
    CHaloBufferSetHandle m_HaloBufferSet;

    static CRationalFieldPointer* GetInstance();
    static CRationalFieldPointer* m_pPointer;
};

template<typename deviceVector, typename deviceGauge, INT vectorN>
class __DLL_EXPORT CFieldFermionKST : public CFieldFermionKS
{
public:

    typedef deviceVector _Vector;
    typedef deviceGauge _Gauge;

    friend class CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>;
    friend class CFieldFermionKSTKernelR<deviceVector, deviceGauge, vectorN>;
    friend class CFieldFermionKSTKernelGamma<deviceVector, deviceGauge, vectorN>;

    typedef CCommonKernelSite<deviceVector> _SiteKernel;
    typedef CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN> _FermionKernel;
    typedef CFieldFermionKSTKernelR<deviceVector, deviceGauge, vectorN> _FermionKernelR;
    typedef CFieldFermionKSTKernelGamma<deviceVector, deviceGauge, vectorN> _FermionKernelGamma;

    CFieldFermionKST()
        : CFieldFermionKS()
        , m_byRationFieldPointerBufferLength(0)
    {
        //Design B halo storage: append halo slots after the local volume in the
        //same device buffer so the move-cache neighbour SIndex can point straight
        //into them. Single-GPU: _HC_HaloSiteCount()==0, so the layout is unchanged.
        m_uiHaloSiteCount = _HC_HaloSiteCount();
        checkCudaErrors(__cudaMalloc((void**)&m_pDeviceData,
            (m_uiSiteCount + m_uiHaloSiteCount) * sizeof(deviceVector)));

        //Improve-1 (multi-GPU-improve1.md 3.1/3.2): bind the halo handle to
        //this exact extent. Site field: one vector per site. The field id is
        //assigned after construction and synced through SetFieldId.
        SHaloBufferInfo sInfo;
        sInfo.m_pDeviceData = reinterpret_cast<BYTE*>(m_pDeviceData);
        sInfo.m_uiCapacityBytes = (m_uiSiteCount + m_uiHaloSiteCount) * sizeof(deviceVector);
        sInfo.m_uiBytesPerSite = sizeof(deviceVector);
        sInfo.m_uiLocalSiteCount = m_uiSiteCount;
        sInfo.m_uiHaloSiteCount = m_uiHaloSiteCount;
        sInfo.m_ullLayoutGeneration = appGetLayoutGeneration();
        sInfo.m_byFieldId = 0;
        sInfo.m_bHaloCapable = TRUE;
        m_HaloBuffer.Bind(sInfo);
    }

    ~CFieldFermionKST()
    {
        checkCudaErrors(__cudaFree(m_pDeviceData));
        m_pDeviceData = NULL;
    }

    void InitialField(EFieldInitialType eInitialType) override
    {
        if (EFIT_RandomZ4 == eInitialType)
        {
            CCommonKernelSite<deviceVector>::InitialBuffer(m_pDeviceData, m_byFieldId, eInitialType);
            NotifyWritten();
            return;
        }
        CCommonKernelField<deviceVector>::Initial(m_pDeviceData, m_uiSiteCount, eInitialType);
        NotifyWritten();
    }

    void ZeroOnEvenOdd(UBOOL bEven) override
    {
        _RECORD(CFieldFermionKST::ZeroOnEvenOdd);
        //CCommonKernelField<deviceVector>::InitialEvenOdd(m_pDeviceData, m_uiSiteCount, EFIT_Zero, TRUE, bEven, FALSE);
        CCommonKernelField<deviceVector>::ZeroEvenOddSite(m_pDeviceData, m_uiSiteCount, bEven);
        NotifyWritten();
    }

    void InitialWithByte(BYTE* byData) override
    {
        CCommonKernelField<deviceVector>::InitialWithByte(m_pDeviceData, m_uiSiteCount, byData);
        NotifyWritten();
    }

    void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFieldType) override
    {
        if (eFieldType != EFFT_CLGBin)
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN>::InitialFieldWithFile: Only support CLG Bin File\n"));
            return;
        }
        const UINT vectorn = _elementdim<deviceVector>();
        const UINT uiBytesPerSite = static_cast<UINT>(sizeof(Real) * vectorn);
        //P5-1.2: on disk the file is the whole GLOBAL lattice in global-site
        //order (CField::SaveToFile gathers to rank 0); read it fully, then
        //scatter to the per-rank sub-lattice on multi-GPU. Previously each rank
        //read by its LOCAL site count -> misaligned/incomplete data, silent
        //wrong (and the local-size check failed on the global file).
        UINT uiSize = 0;
        BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        if (NULL == data)
        {
            appCrucial(_T("File not found: %s\n"), sFileName.c_str());
            _FAIL_EXIT;
        }
        const UINT uiPerRankBytes = uiBytesPerSite * m_uiSiteCount;
#if _CLG_MULTI_GPU
        if (NULL != appGetComm() && appGetComm()->Size() > 1)
        {
            if (uiSize != uiPerRankBytes * appGetComm()->Size())
            {
                appCrucial(_T("File size not correct (MG): expecting global %d, found: %d\n"),
                    static_cast<UINT>(uiPerRankBytes * appGetComm()->Size()), uiSize);
                _FAIL_EXIT;
            }
            BYTE* byLocal = (BYTE*)malloc(uiPerRankBytes);
            appGetComm()->ScatterFieldFromRoot(data, uiBytesPerSite, byLocal);
            free(data);
            data = byLocal;
        }
        else
#endif
        if (uiSize != uiPerRankBytes)
        {
            appCrucial(_T("File size not correct: expecting: %d, found: %d\n"), uiPerRankBytes, uiSize);
            _FAIL_EXIT;
        }
        InitialWithByte(data);
        free(data);
        NotifyWritten();
    }

    void InitialOtherParameters(CParameters& params) override
    {
        CFieldFermionKS::InitialOtherParameters(params);
        m_byRationFieldPointerBufferLength = static_cast<BYTE>(2 * GRASet.m_pRASet[m_iMDIndex]->m_uiDegree);
    }

    void DebugPrintMe() const override
    {
        CCommonKernelSite<deviceVector>::DebugPrint(m_pDeviceData, m_uiSiteCount);
    }

    void CopyParamTo(CField* U) const override
    {
        CFieldFermionKS::CopyParamTo(U);
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(U);
        pField->m_byRationFieldPointerBufferLength = static_cast<BYTE>(2 * GRASet.m_pRASet[m_iMDIndex]->m_uiDegree);
    }

    void CopyBufferTo(CField* U) const override
    {
        if (NULL == U || GetFieldType() != U->GetFieldType())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
            return;
        }

        //CFieldFermionKS::CopyBufferTo(U);
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(U);
        CCommonKernel<deviceVector>::CopyBuffer(pField->m_pDeviceData, m_pDeviceData, m_uiSiteCount);
        U->NotifyWritten();
    }

    void Dagger() override
    {
        _RECORD(CFieldFermionKST::Dagger);
        CCommonKernelField<deviceVector>::Dagger(m_pDeviceData, m_uiSiteCount);
        NotifyWritten();
    }

    //This is Axpy(1.0f, x)
    void AxpyPlus(const CField* x) override
    {
        _RECORD(CFieldFermionKST::AxpyP);
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
            return;
        }
        const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);

        //For even pseudo-fermion, odd sites are always zero, it is safe to update even sites only
        if (m_bEvenPseudofermion)
        {
            CCommonKernelField<deviceVector>::AxpyPlusEvenOdd(m_pDeviceData, m_uiSiteCount, pField->m_pDeviceData, TRUE);
            NotifyWritten();
            return;
        }
        CCommonKernelField<deviceVector>::AxpyPlus(m_pDeviceData, m_uiSiteCount, pField->m_pDeviceData);
        NotifyWritten();
    }

    void AxpyMinus(const CField* x) override
    {
        _RECORD(CFieldFermionKST::AxpyM);
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
            return;
        }
        const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);

        //For even pseudo-fermion, odd sites are always zero, it is safe to update even sites only
        if (m_bEvenPseudofermion)
        {
            CCommonKernelField<deviceVector>::AxpyMinusEvenOdd(m_pDeviceData, m_uiSiteCount, pField->m_pDeviceData, TRUE);
            NotifyWritten();
            return;
        }
        CCommonKernelField<deviceVector>::AxpyMinus(m_pDeviceData, m_uiSiteCount, pField->m_pDeviceData);
        NotifyWritten();
    }

    void Axpy(Real a, const CField* x) override
    {
        _RECORD(CFieldFermionKST::AxpyR);
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
            return;
        }
        const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);

        //For even pseudo-fermion, odd sites are always zero, it is safe to update even sites only
        if (m_bEvenPseudofermion)
        {
            CCommonKernelField<deviceVector>::AxpyEvenOdd(m_pDeviceData, m_uiSiteCount, a, pField->m_pDeviceData, TRUE);
            NotifyWritten();
            return;
        }
        CCommonKernelField<deviceVector>::Axpy(m_pDeviceData, m_uiSiteCount, a, pField->m_pDeviceData);
        NotifyWritten();
    }

    void Axpy(const CLGComplex& a, const CField* x) override
    {
        _RECORD(CFieldFermionKST::AxpyC);
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
            return;
        }
        const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);

        //For even pseudo-fermion, odd sites are always zero, it is safe to update even sites only
        if (m_bEvenPseudofermion)
        {
            CCommonKernelField<deviceVector>::AxpyEvenOdd(m_pDeviceData, m_uiSiteCount, a, pField->m_pDeviceData, TRUE);
            NotifyWritten();
            return;
        }
        CCommonKernelField<deviceVector>::Axpy(m_pDeviceData, m_uiSiteCount, a, pField->m_pDeviceData);
        NotifyWritten();
    }

    void Mul(const CField* other, UBOOL bDaggerLeft = TRUE, UBOOL bDaggerRight = FALSE) override
    {
        if (NULL == other || GetFieldType() != other->GetFieldType())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
            return;
        }
        const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(other);

        CCommonKernelField<deviceVector>::Mul(m_pDeviceData, m_uiSiteCount, pField->m_pDeviceData, bDaggerLeft, bDaggerRight);
        NotifyWritten();
    }

    void LeftMul(const CField* other, UBOOL bDaggerLeft = FALSE, UBOOL bDaggerRight = FALSE) override
    {
        if (NULL == other || GetFieldType() != other->GetFieldType())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
            return;
        }
        const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(other);

        CCommonKernelField<deviceVector>::LeftMul(m_pDeviceData, m_uiSiteCount, pField->m_pDeviceData, bDaggerLeft, bDaggerRight);
        NotifyWritten();
    }

    void ScalarMultply(const CLGComplex& a) override
    {
        CCommonKernelField<deviceVector>::ScalarMultply(m_pDeviceData, m_uiSiteCount, a);
        NotifyWritten();
    }

    void ScalarMultply(Real a) override
    {
        CCommonKernelField<deviceVector>::ScalarMultply(m_pDeviceData, m_uiSiteCount, a);
        NotifyWritten();
    }

    cuDoubleComplex Dot(const CField* x) const override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
            return make_cuDoubleComplex(0, 0);
        }
        _RECORD(CFieldFermionKST::Dot);
        const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);
        return CCommonKernelField<deviceVector>::Dot(m_pDeviceData, m_uiSiteCount, pField->m_pDeviceData);
    }

    DOUBLE GetLength() const override
    {
        _RECORD(CFieldFermionKST::GetLength);
        return CCommonKernelField<deviceVector>::LengthSq(m_pDeviceData, m_uiSiteCount);
    }

    //Multi-GPU (Phase 2): per-site element size used by the halo manager to stage
    //this field's boundary planes through host buffers (Design B refill).
    UINT GetSiteElementBytes() const override { return static_cast<UINT>(sizeof(deviceVector)); }

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
        if (m_bEvenPseudofermion)
        {
            appCrucial(_T("DS not implemented for even pseudo-fermion, use DEven, note that this is only valid for test!\n"));
            D0OnEvenOrOddS(pGauge, TRUE);
            ZeroOnEvenOdd(TRUE);
            return;
        }
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));

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
        if (m_bEvenPseudofermion)
        {
            appCrucial(_T("DdaggerS not implemented\n"));
        }

        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));
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
        if (m_bEvenPseudofermion)
        {
            appCrucial(_T("DDS not implemented\n"));
        }

        Real fRealCoeff = fCoeffReal;
        const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
        if (EOCT_Minus == eCoeffType)
        {
            eCoeffType = EOCT_Real;
            fRealCoeff = F(-1.0);
        }
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));

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
        _RECORD(CFieldFermionKST::DDdaggerS);
        Real fRealCoeff = fCoeffReal;
        const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
        if (EOCT_Minus == eCoeffType)
        {
            eCoeffType = EOCT_Real;
            fRealCoeff = F(-1.0);
        }

        if (m_bEvenPseudofermion)
        {
            //appParanoiac(_T("DDdaggerS on even\n"));
            {
                _RECORD2(CFieldFermionKST::DDdaggerS::DOperatorKSOnEvenOrOdd, a);
                DOperatorKSOnEvenOrOdd(m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, TRUE, F(0.0),
                    TRUE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
            }
            {
                _RECORD2(CFieldFermionKST::DDdaggerS::DOperatorKSOnEvenOrOdd, a);
                DOperatorKSOnEvenOrOdd(m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, FALSE, m_f2am * m_f2am,
                    FALSE, eCoeffType, fRealCoeff, cCompCoeff);
            }
            {
                _RECORD2(CFieldFermionKST::DDdaggerS::ZeroOnEvenOdd, a);
                ZeroOnEvenOdd(FALSE);
            }
        }
        else
        {
            CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));

            DOperator(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
                TRUE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
            //why only apply coeff in the next step?
            DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
                FALSE, eCoeffType, fRealCoeff, cCompCoeff);

            pPooled->Return();
        }
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
            appCrucial(_T("DWithMassS does not support diagonal mass.\n"));
        }
        if (m_bEvenPseudofermion)
        {
            appCrucial(_T("DWithMassS not implemented\n"));
        }

        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));

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
            appCrucial(_T("DdaggerWithMassS does not support diagonal mass.\n"));
        }
        if (m_bEvenPseudofermion)
        {
            appCrucial(_T("DdaggerWithMassS not implemented\n"));
        }

        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));
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
        if (m_bDiagonalMass)
        {
            appCrucial(_T("DDWithMassS does not support diagonal mass.\n"));
        }
        if (m_bEvenPseudofermion)
        {
            appCrucial(_T("DDWithMassS not implemented\n"));
        }

        Real fRealCoeff = fCoeffReal;
        const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
        if (EOCT_Minus == eCoeffType)
        {
            eCoeffType = EOCT_Real;
            fRealCoeff = F(-1.0);
        }
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));

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
            appCrucial(_T("DDdaggerWithMassS does not support diagonal mass.\n"));
        }
        if (m_bEvenPseudofermion)
        {
            appCrucial(_T("DDdaggerWithMassS not implemented\n"));
        }

        Real fRealCoeff = fCoeffReal;
        const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
        if (EOCT_Minus == eCoeffType)
        {
            eCoeffType = EOCT_Real;
            fRealCoeff = F(-1.0);
        }
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));

        DOperatorKS(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
            TRUE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
        //why only apply coeff in the next step?
        DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
            FALSE, eCoeffType, fRealCoeff, cCompCoeff);

        pPooled->Return();
    }

public:

    TArray<CFieldFermion*> GetSourcesAtSiteFromPool(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const SSmallInt4& site) const override
    {
        TArray<CFieldFermion*> ret;
        for (UINT j = 0; j < _dim<deviceVector>(); ++j)
        {
            ret.AddItem(dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__)));
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

            ret[c]->InverseD(gaugeNum, bosonNum, 0, gaugeFields, pBoson, NULL);
        }
        return ret;
    }

    /**
    * generate phi by gaussian random.
    * phi = (D^+D)^{1/4} phi
    */
    void PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson) override
    {
        if (m_bEvenPseudofermion)
        {
            appParanoiac(_T("PrepareForHMC with Even field\n"));
            CCommonKernelSite<deviceVector>::InitialBufferEvenOdd(m_pDeviceData, m_byFieldId, TRUE, TRUE, EFIT_RandomGaussian);

            //this->DebugPrintMe();
            //Note that, using D0_MC, A=D0D0 + m^2, A^{-1/2} should be expand with m^2 considered.
            D_MC(gaugeNum, bosonNum, 0, gaugeFields, pBoson, NULL);
        }
        else
        {
            CCommonKernelSite<deviceVector>::InitialBuffer(m_pDeviceData, m_byFieldId, EFIT_RandomGaussian);
            D_MC(gaugeNum, bosonNum, 0, gaugeFields, pBoson, NULL);
        }

        //if (NULL != appGetFermionSolver(m_byFieldId) && !appGetFermionSolver(m_byFieldId)->IsAbsoluteAccuracy())
        //{
        //    m_fLength = Dot(this).x;
        //}

        //For KS, we generally use shifted solver, so do NOT cache the last result
    }

protected:

    void ApplyGammaS(const CFieldGauge* pGauge, EGammaMatrix eGamma) override
    {
        if (NULL == pGauge || vectorN != pGauge->MatrixN())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
            return;
        }
        const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));
        CCommonKernel<deviceVector>::CopyBuffer(pPooled->m_pDeviceData, m_pDeviceData, m_uiSiteCount);
        InitialField(EFIT_Zero);

        //If it was gamma_mu or gamma_5 or sigmaij, it is i gamma mu and i gamma 5, therefore multiply -i
        UBOOL bImag = (GAMMA1 == eGamma) || (GAMMA2 == eGamma) || (GAMMA3 == eGamma) || (GAMMA4 == eGamma) || (GAMMA5 == eGamma)
            || (SIGMA12 == eGamma) || (SIGMA31 == eGamma) || (SIGMA41 == eGamma) || (SIGMA23 == eGamma) || (SIGMA42 == eGamma) || (SIGMA43 == eGamma);

        CFieldFermionKSTKernelGamma<deviceVector, deviceGauge, vectorN>::appApplyGammaKS(
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
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));

        CCommonKernel<deviceVector>::CopyBuffer(pPooled->m_pDeviceData, m_pDeviceData, m_uiSiteCount);

        DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, F(0.0),
            FALSE, EOCT_None, F(1.0), _onec);

        pPooled->Return();
    }

    void D0OnEvenOrOddS(const CField* pGauge, UBOOL bEven) override
    {
        _RECORD(CFieldFermionKST::D0OnEvenOrOddS);
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
        if (m_bEachSiteEta && NULL != dynamic_cast<const CBoundaryConditionProjectivePlaneSquare*>(
            appGetLattice()->m_pIndex->GetBoudanryCondition()))
        {
            appCrucial(_T("Even odd does not support projective plane boundary.\n"));
        }
        appParanoiac(_T("D0OnEvenOrOddS\n"));
        DOperatorKSOnEvenOrOdd(m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, bEven, F(0.0),
            FALSE, EOCT_None, F(1.0), _onec);
    }

    UBOOL CalculateForceS(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const override
    {
        if (m_bEvenPseudofermion)
        {
            //CFieldGauge* thisfieldforce = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(pGauge->m_byFieldId));
            CalculateForceEvenOddS(pGauge, pForce, ePhase);
            //thisfieldforce->LeftMul(pGauge, FALSE, TRUE);
            //thisfieldforce->TA();
            //pForce->AxpyPlus(thisfieldforce);
            //thisfieldforce->Return();
            return TRUE;
        }

        //=================================== test ================================
        //CFieldGauge* thisfieldforce = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(pGauge->m_byFieldId));
        //CalculateForceTest(pGauge, thisfieldforce, ePhase);
        //thisfieldforce->LeftMul(pGauge, FALSE, TRUE);
        //thisfieldforce->TA();
        //pForce->AxpyPlus(thisfieldforce);
        //thisfieldforce->Return();
        //return TRUE;
        //=========================================================================

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
        for (UINT i = 0; i < GRASet.m_pRASet[m_iMDIndex]->m_uiDegree; ++i)
        {
            CField* pPhi_i = dynamic_cast<CField*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));
            phii.AddItem(pPhi_i);
            CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPhi_id = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));
            phiid.AddItem(pPhi_id);
        }

        CMultiShiftSolver* solver = appGetMultiShiftSolver(m_byFieldId);
        if (NULL == solver)
        {
            appCrucial(_T("No multi solver found!"));
            _FAIL_EXIT;
        }

        TArray<CLGComplex> shifts;
        for (UINT i = 0; i < GRASet.m_pRASet[m_iMDIndex]->m_uiDegree; ++i)
        {
            shifts.AddItem(_make_cuComplex(GRASet.m_pRASet[m_iMDIndex]->m_lstB[i], F(0.0)));
        }
        TArray<const CFieldGauge*> gaues;
        gaues.AddItem(pGauge);
        solver->Solve(phii, shifts, this, 1, 0, 0, gaues.GetData(), NULL, NULL, EFO_F_DDdagger);

        const UINT uiBufferSize = sizeof(deviceVector*) * 2 * GRASet.m_pRASet[m_iMDIndex]->m_uiDegree;
        deviceVector** hostPointers = (deviceVector**)appAlloca(uiBufferSize);
        for (UINT i = 0; i < GRASet.m_pRASet[m_iMDIndex]->m_uiDegree; ++i)
        {
            CFieldFermionKST<deviceVector, deviceGauge, vectorN>* phi_ks = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(phii[i]);
            phi_ks->FixBoundary(EFB_Field);
            phi_ks->CopyTo(phiid[i]);
            if (m_bDiagonalMass)
            {
                //phiid[i]->DS(pGauge);
                phiid[i]->DdaggerS(pGauge, EOCT_Minus);
            }
            else
            {
                phiid[i]->D0S(pGauge);
            }
            phiid[i]->FixBoundary(EFB_Field);
            hostPointers[i] = phi_ks->m_pDeviceData;
            hostPointers[i + GRASet.m_pRASet[m_iMDIndex]->m_uiDegree] = phiid[i]->m_pDeviceData;
        }

        appSimpleCopyHD(CRationalFieldPointer::GetInstance()->GetRationPoint<deviceVector>(m_byRationFieldPointerBufferLength), hostPointers, uiBufferSize);

        //Improve-1 (I6, 3.4.1): publish the current pooled membership to the
        //halo registry; a launch argument hitting the array expands to these
        //member handles (Ensure on reads, invalidate on writes).
        {
            TArray<const CHaloBufferHandle*> memberHandles;
            for (UINT i = 0; i < GRASet.m_pRASet[m_iMDIndex]->m_uiDegree; ++i)
            {
                if (NULL != phii[i]->GetHaloBufferHandle())
                {
                    memberHandles.AddItem(phii[i]->GetHaloBufferHandle());
                }
                if (NULL != phiid[i]->GetHaloBufferHandle())
                {
                    memberHandles.AddItem(phiid[i]->GetHaloBufferHandle());
                }
            }
            CRationalFieldPointer::GetInstance()->RebindMemberHandles(memberHandles);
        }

        const CFieldGaugeLink<deviceGauge, vectorN>* pGaugeSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
        CFieldGaugeLink<deviceGauge, vectorN>* pForceSU3 = dynamic_cast<CFieldGaugeLink<deviceGauge, vectorN>*>(pForce);

        DerivateD0(pForceSU3->m_pDeviceData, pGaugeSU3->m_pDeviceData, pGaugeSU3->m_byFieldId);

        for (UINT i = 0; i < GRASet.m_pRASet[m_iMDIndex]->m_uiDegree; ++i)
        {
            phii[i]->Return();
            phiid[i]->Return();
        }
        return TRUE;
    }

    void CalculateForceEvenOddS(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const override
    {
        appParanoiac(_T("CalculateForceEvenOddS\n"));
        TArray<CField*> shiftsolutions;
        RationalApproximationPooled(EFO_F_DDdagger, 1, 0, 0, &pGauge, NULL, NULL, m_iMDIndex, shiftsolutions);
        for (INT i = 0; i < shiftsolutions.Num(); ++i)
        {
            CFieldFermionKST<deviceVector, deviceGauge, vectorN>* fks = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(shiftsolutions[i]);
            fks->FixBoundary(EFB_Field);
            fks->D0OnEvenOrOddS(pGauge, TRUE);
            fks->FixBoundary(EFB_Field);
            //pForce->AddConnection(fks, GRASet.m_pRASet[m_iMDIndex]->m_lstA[i]);
            CalculateForceEvenOddS_SingleTermOfRational(pGauge, fks, pForce, GRASet.m_pRASet[m_iMDIndex]->m_lstA[i], i);
            fks->Return();
        }
    }

    void CalculateForceEvenOddS_SingleTermOfRational(const CFieldGauge* pGauge, const CFieldFermionKS* phi, CFieldGauge* pForce, Real fCoef, INT iRationalTermIndex) const override
    {
        pForce->AddConnection(phi, fCoef);
    }

    // right now, we use DerivateD0, which stores all shift solver results to m_pRationalFieldPointers
    // may change to 'CalculateForceS_SingleTermOfRational' in the future
    //void CalculateForceS_SingleTermOfRational(const CFieldFermionKS* phi, const CFieldFermionKS* phid, CFieldGauge* pForce, Real fCoef) const override
    //{
    //    appCrucial(_T("CalculateForceS_SingleTermOfRational not implemented\n"));
    //}

    TArray<CField*> CalculateRationalFields(const CFieldGauge* pGauge) const override
    {
        if (NULL == pGauge || vectorN != pGauge->MatrixN())
        {
            appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
            return TArray<CField*>();
        }

        TArray<CField*> phii;
        TArray<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*> phiid;
        for (UINT i = 0; i < GRASet.m_pRASet[m_iMDIndex]->m_uiDegree; ++i)
        {
            CField* pPhi_i = dynamic_cast<CField*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));
            phii.AddItem(pPhi_i);
            CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPhi_id = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));
            phiid.AddItem(pPhi_id);
        }

        CMultiShiftSolver* solver = appGetMultiShiftSolver(m_byFieldId);
        if (NULL == solver)
        {
            appCrucial(_T("No multi solver found!"));
            _FAIL_EXIT;
        }

        TArray<CLGComplex> shifts;
        for (UINT i = 0; i < GRASet.m_pRASet[m_iMDIndex]->m_uiDegree; ++i)
        {
            shifts.AddItem(_make_cuComplex(GRASet.m_pRASet[m_iMDIndex]->m_lstB[i], F(0.0)));
        }
        TArray<const CFieldGauge*> gaues;
        gaues.AddItem(pGauge);
        solver->Solve(phii, shifts, this, 1, 0, 0, gaues.GetData(), NULL, NULL, EFO_F_DDdagger);

        const UINT uiBufferSize = sizeof(deviceVector*) * 2 * GRASet.m_pRASet[m_iMDIndex]->m_uiDegree;
        deviceVector** hostPointers = (deviceVector**)appAlloca(uiBufferSize);
        for (UINT i = 0; i < GRASet.m_pRASet[m_iMDIndex]->m_uiDegree; ++i)
        {
            CFieldFermionKST<deviceVector, deviceGauge, vectorN>* phi_ks = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(phii[i]);
            phi_ks->CopyTo(phiid[i]);
            if (m_bDiagonalMass)
            {
                //phiid[i]->DS(pGauge);
                phiid[i]->DdaggerS(pGauge, EOCT_Minus);
            }
            else
            {
                phiid[i]->D0S(pGauge);
            }

            hostPointers[i] = phi_ks->m_pDeviceData;
            hostPointers[i + GRASet.m_pRASet[m_iMDIndex]->m_uiDegree] = phiid[i]->m_pDeviceData;
        }

        appSimpleCopyHD(CRationalFieldPointer::GetInstance()->GetRationPoint<deviceVector>(m_byRationFieldPointerBufferLength), hostPointers, uiBufferSize);

        //Improve-1 (I6, 3.4.1): publish the current pooled membership to the
        //halo registry; a launch argument hitting the array expands to these
        //member handles (Ensure on reads, invalidate on writes).
        {
            TArray<const CHaloBufferHandle*> memberHandles;
            for (UINT i = 0; i < GRASet.m_pRASet[m_iMDIndex]->m_uiDegree; ++i)
            {
                if (NULL != phii[i]->GetHaloBufferHandle())
                {
                    memberHandles.AddItem(phii[i]->GetHaloBufferHandle());
                }
                if (NULL != phiid[i]->GetHaloBufferHandle())
                {
                    memberHandles.AddItem(phiid[i]->GetHaloBufferHandle());
                }
            }
            CRationalFieldPointer::GetInstance()->RebindMemberHandles(memberHandles);
        }

        //const CFieldGaugeLink<deviceGauge, vectorN>* pGaugeSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
        //CFieldGaugeLink<deviceGauge, vectorN>* pForceSU3 = dynamic_cast<CFieldGaugeLink<deviceGauge, vectorN>*>(pForce);

        //DerivateD0(pForceSU3->m_pDeviceData, pGaugeSU3->m_pDeviceData, pGaugeSU3->m_byFieldId);

        TArray<CField*> rets;
        for (UINT i = 0; i < GRASet.m_pRASet[m_iMDIndex]->m_uiDegree; ++i)
        {
            rets.AddItem(phii[i]);
            rets.AddItem(phiid[i]);
        }

        return rets;
    }

public:

    //For test only
    void PrepareForHMCOnlyRandomize() override
    {
        CCommonKernelSite<deviceVector>::InitialBuffer(m_pDeviceData, m_byFieldId, EFIT_RandomGaussian);
        NotifyWritten();
    }

    void PrepareForHMCNotRandomize(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson) override
    {
        D_MC(gaugeNum, bosonNum, 0, gaugeFields, pBoson, NULL);
    }

    void InitialAsSource(const SFermionBosonSource& sourceData) override
    {
        CCommonKernelSite<deviceVector>::InitialSource(m_pDeviceData, m_byFieldId, sourceData);
        NotifyWritten();
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
        CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DerivateD0(m_pDeviceData, m_byFieldId, (deviceGauge*)pForce, (const deviceGauge*)pGaugeBuffer, byGaugeFieldId, 
            CRationalFieldPointer::GetInstance()->GetRationPoint<deviceVector>(m_byRationFieldPointerBufferLength),
            GRASet.m_pRASet[m_iMDIndex]->m_pDeviceData, GRASet.m_pRASet[m_iMDIndex]->m_uiDegree);
    }

    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override
    {
        CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DOperatorKS(m_bEachSiteEta, (deviceVector*)pTargetBuffer, (const deviceVector*)pBuffer, (const deviceGauge*)pGaugeBuffer, m_byFieldId, byGaugeFieldId, f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);
    }

    //even of target buffer is kept
    //odd of target buffer is changed to D0 buffer
    void DOperatorKSOnEvenOrOdd(void* pTargetBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, UBOOL bEven, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override
    {
        CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DOperatorKSOnEvenOrOdd((deviceVector*)pTargetBuffer, (const deviceGauge*)pGaugeBuffer, m_byFieldId, byGaugeFieldId, bEven, f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);
    }

public:

    #pragma region Help functions to implement higher orders

    void OnlyMass(void* pTarget, Real f2am, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override
    {
        if (m_bDiagonalMass)
        {
            appCrucial(_T("OnlyMass does not support diagonal mass.\n"));
        }
        CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::OnlyMass(m_pDeviceData, (deviceVector*)pTarget, f2am, eOCT, fRealCoeff, cCmpCoeff);
    }

    #pragma endregion

    void Connection(const CField* other, void* res) const override
    {
        CCommonKernelMV<deviceVector, deviceGauge, vectorN>::Connection(m_pDeviceData, (const deviceVector*)other->GetData(), (deviceGauge*)res, m_byFieldId);
    }

    void ConnectionSelf(void* res, Real fCoeff) const override
    {
        CCommonKernelMV<deviceVector, deviceGauge, vectorN>::ConnectionOneFieldStaggered(m_pDeviceData, (deviceGauge*)res, m_byFieldId, fCoeff);
    }

    void AddConnectionSelf(void* res, Real fCoeff) const override
    {
        CCommonKernelMV<deviceVector, deviceGauge, vectorN>::AddConnectionOneFieldStaggered(m_pDeviceData, (deviceGauge*)res, m_byFieldId, fCoeff);
    }

    void AddConnectionTwoField(void* res, const CField* other, Real fCoeff) const override
    {
        CCommonKernelMV<deviceVector, deviceGauge, vectorN>::AddConnection(m_pDeviceData, (const deviceVector*)other->GetData(), (deviceGauge*)res, fCoeff, m_byFieldId);
    }

    deviceVector* m_pDeviceData;

    //Improve-1: this field's own buffer handle (pool copies each bind theirs).
    CHaloBufferHandle m_HaloBuffer;
    CHaloBufferHandle* GetHaloBufferHandle() override { return &m_HaloBuffer; }
    const CHaloBufferHandle* GetHaloBufferHandle() const override { return &m_HaloBuffer; }

    _GetData

//protected:

    void OneLinkS(const void* pGauge, BYTE byGaugeFieldId, void* pTarget, Real fCoefficient,
        const SCHAR* pDevicePath, BYTE pathLength, BYTE byEtaIdx,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override
    {
        CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::OneLinkS(m_pDeviceData, m_byFieldId, (const deviceGauge*)pGauge, byGaugeFieldId, (deviceVector*)pTarget, fCoefficient, pDevicePath, pathLength, byEtaIdx, bDagger, eOCT, fRealCoeff, cCmpCoeff);
    }

    void OneLinkForceS(const void* pGauge, BYTE byGaugeFieldId, void* pForce, Real fCoefficient,
        const SCHAR* pDevicePath, BYTE pathLength, BYTE byEtaIdx) const override
    {
        CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::OneLinkForceS(m_pDeviceData, m_byFieldId, (const deviceGauge*)pGauge, byGaugeFieldId, (deviceGauge*)pForce, fCoefficient, pDevicePath, pathLength, byEtaIdx, 
            CRationalFieldPointer::GetInstance()->GetRationPoint<deviceVector>(m_byRationFieldPointerBufferLength),
            GRASet.m_pRASet[m_iMDIndex]->m_pDeviceData, GRASet.m_pRASet[m_iMDIndex]->m_uiDegree);
    }

protected:

    //phi _i and Dst0 phi _i
    //deviceVector** m_pRationalFieldPointers;
    BYTE m_byRationFieldPointerBufferLength;

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
        checkCudaErrors(__cudaMalloc((void**)&m_pResBuffer, sizeof(deviceVector*) * _kFieldMatrixMaxDim));
        checkCudaErrors(__cudaMalloc((void**)&m_pLeftBuffer, sizeof(deviceVector*) * _kFieldMatrixMaxDim));
    }
    ~CFieldMatrixOperationKST()
    {
        free(m_pHostResBuffer);
        free(m_pHostLeftBuffer);
        checkCudaErrors(__cudaFree(m_pResBuffer));
        checkCudaErrors(__cudaFree(m_pLeftBuffer));
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
_DEFINE_KS_FIELDSUN(3)
_DEFINE_KS_FIELDSUN(4)
//_DEFINE_KS_FIELDSUN(5)
//_DEFINE_KS_FIELDSUN(6)
//_DEFINE_KS_FIELDSUN(7)
//_DEFINE_KS_FIELDSUN(8)

typedef CFieldMatrixOperationKST<deviceSU3Vector, deviceSU3, 3> CFieldMatrixOperationKSSU3;

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKST_H_

//=============================================================================
// END OF FILE
//=============================================================================
