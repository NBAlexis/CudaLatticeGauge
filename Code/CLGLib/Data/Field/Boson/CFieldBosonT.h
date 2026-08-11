//=============================================================================
// FILENAME : CFieldBosonT.h
// 
// DESCRIPTION:
// This is the class for all boson fields
//
// REVISION:
//  [mm/dd/yy]
//  [3/31/2024 nbale]
//=============================================================================

#ifndef _CFIELDBOSONT_H_
#define _CFIELDBOSONT_H_


__BEGIN_NAMESPACE

template<typename deviceDataBoson>
class __DLL_EXPORT CFieldBosonT : public CFieldBoson
{
public:
    CFieldBosonT()
        : CFieldBoson()
        , m_pDeviceData(NULL)
    {
        //P4-4.1: append the halo site capacity (Design B flat buffer, same as
        //CFieldGaugeLink): stencil kernels read split-direction out-of-lattice
        //neighbours from the appended slots. m_uiSiteCount stays the physics
        //volume so CopyDataOut / gather-scatter are unchanged; only the
        //allocation grows. Single-GPU: _HC_HaloSiteCount() == 0 -> identical
        //size (bit-identical behaviour).
        const UINT uiHaloSiteCount = _HC_HaloSiteCount();
        checkCudaErrors(__cudaMalloc((void**)& m_pDeviceData, (m_uiSiteCount + uiHaloSiteCount) * sizeof(deviceDataBoson)));

        //Improve-1 (multi-GPU-improve1.md 3.1/3.2): bind the halo handle to
        //this exact extent. Site field: one element per site. The field id is
        //assigned after construction and synced through SetFieldId.
        SHaloBufferInfo sInfo;
        sInfo.m_pDeviceData = reinterpret_cast<BYTE*>(m_pDeviceData);
        sInfo.m_uiCapacityBytes = (m_uiSiteCount + uiHaloSiteCount) * sizeof(deviceDataBoson);
        sInfo.m_uiBytesPerSite = sizeof(deviceDataBoson);
        sInfo.m_uiLocalSiteCount = m_uiSiteCount;
        sInfo.m_uiHaloSiteCount = uiHaloSiteCount;
        sInfo.m_ullLayoutGeneration = appGetLayoutGeneration();
        sInfo.m_byFieldId = 0;
        sInfo.m_bHaloCapable = TRUE;
        m_HaloBuffer.Bind(sInfo);
    }

    ~CFieldBosonT()
    {
        checkCudaErrors(__cudaFree(m_pDeviceData));
        m_pDeviceData = NULL;
    }

    void CopyBufferTo(CField* U) const override
    {
        if (NULL == U || GetFieldType() != U->GetFieldType())
        {
            appCrucial(_T("EFT_BosonU1 can only copy to EFT_BosonU1!"));
            return;
        }
        CFieldBosonT<deviceDataBoson>* pField = dynamic_cast<CFieldBosonT<deviceDataBoson>*>(U);
        CCommonKernel<deviceDataBoson>::CopyBuffer(pField->m_pDeviceData, m_pDeviceData, m_uiSiteCount);
        U->NotifyWritten();
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
        NotifyWritten();
    }

    void InitialField(EFieldInitialType eInitialType) override
    {
        CCommonKernelSite<deviceDataBoson>::InitialBuffer(m_pDeviceData, m_byFieldId, eInitialType);
        NotifyWritten();
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
#if _CLG_MULTI_GPU
        //Whole-global-lattice file: read it fully on every rank so root can
        //scatter (non-root's copy is discarded by the scatter below).
        if (NULL != appGetComm() && appGetComm()->Size() > 1)
        {
            uiReadSize = uiSize * appGetComm()->Size();
        }
#endif
        BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiReadSize);
        if (NULL == data)
        {
            appCrucial(_T("File not found: %s\n"), sFileName.c_str());
            _FAIL_EXIT;
        }

#if _CLG_MULTI_GPU
        if (NULL != appGetComm() && appGetComm()->Size() > 1)
        {
            //Whole-global-lattice file: expect the GLOBAL size here (the local
            //comparison below would reject the 2x (Size x) global read).
            if (uiReadSize != uiSize * appGetComm()->Size())
            {
                appCrucial(_T("File size not correct (MG): expecting: %d, found: %d\n"), uiSize * appGetComm()->Size(), uiReadSize);
                _FAIL_EXIT;
            }
        }
        else
#endif
        if (uiSize != uiReadSize)
        {
            appCrucial(_T("File size not correct: expecting: %d, found: %d\n"), uiSize, uiReadSize);
            _FAIL_EXIT;
        }
#if _CLG_MULTI_GPU
        //Whole-global-lattice file (SaveToFile gather); scatter per rank.
        if (NULL != appGetComm() && appGetComm()->Size() > 1)
        {
            const UINT uiPerRank = static_cast<UINT>(sizeof(Real) * FloatN() * _HC_Volume);
            if (uiReadSize != uiPerRank * appGetComm()->Size())
            {
                appCrucial(_T("Loading file size not match (MG): %s, %d, expecting global %d\n"),
                    sFileName.c_str(), uiReadSize, static_cast<INT>(uiPerRank * appGetComm()->Size()));
                _FAIL_EXIT;
            }
            const UINT uiBytesPerSite = static_cast<UINT>(sizeof(Real) * FloatN());
            BYTE* byLocal = (BYTE*)malloc(uiPerRank);
            appGetComm()->ScatterFieldFromRoot(data, uiBytesPerSite, byLocal);
            free(data);
            data = byLocal;
        }
#endif

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
        NotifyWritten();
    }

    void InitialWithByte(BYTE* byData) override
    {
        CCommonKernelField<deviceDataBoson>::InitialWithByte(m_pDeviceData, m_uiSiteCount, byData);
        NotifyWritten();
    }

    void DebugPrintMe() const override
    {
        CCommonKernelSite<deviceDataBoson>::DebugPrint(m_pDeviceData, m_uiSiteCount);
    }

    void Dagger() override
    {
        CCommonKernelField<deviceDataBoson>::Dagger(m_pDeviceData, m_uiSiteCount);
        NotifyWritten();
    }

    void FixBoundary(EFixBoundary eType) override
    {
        CCommonKernelSite<deviceDataBoson>::FixBoundary(m_pDeviceData, m_byFieldId);
        NotifyWritten();
    }

    void AxpyPlus(const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonVN!"));
            return;
        }
        const CFieldBosonT<deviceDataBoson>* pField = dynamic_cast<const CFieldBosonT<deviceDataBoson>*>(x);
        CCommonKernelField<deviceDataBoson>::AxpyPlus(m_pDeviceData, m_uiSiteCount, pField->m_pDeviceData);
        NotifyWritten();
    }

    void AxpyMinus(const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonVN!"));
            return;
        }
        const CFieldBosonT<deviceDataBoson>* pField = dynamic_cast<const CFieldBosonT<deviceDataBoson>*>(x);
        CCommonKernelField<deviceDataBoson>::AxpyMinus(m_pDeviceData, m_uiSiteCount, pField->m_pDeviceData);
        NotifyWritten();
    }

    void Axpy(Real a, const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonVN!"));
            return;
        }
        const CFieldBosonT<deviceDataBoson>* pField = dynamic_cast<const CFieldBosonT<deviceDataBoson>*>(x);
        CCommonKernelField<deviceDataBoson>::Axpy(m_pDeviceData, m_uiSiteCount, a, pField->m_pDeviceData);
        NotifyWritten();
    }

    void Axpy(const CLGComplex& a, const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonVN!"));
            return;
        }
        const CFieldBosonT<deviceDataBoson>* pField = dynamic_cast<const CFieldBosonT<deviceDataBoson>*>(x);
        CCommonKernelField<deviceDataBoson>::Axpy(m_pDeviceData, m_uiSiteCount, a, pField->m_pDeviceData);
        NotifyWritten();
    }


    void Mul(const CField* x, UBOOL bDaggerLeft = TRUE, UBOOL bDaggerRight = FALSE) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonVN!"));
            return;
        }
        const CFieldBosonT<deviceDataBoson>* pField = dynamic_cast<const CFieldBosonT<deviceDataBoson>*>(x);
        CCommonKernelField<deviceDataBoson>::Mul(m_pDeviceData, m_uiSiteCount, pField->m_pDeviceData, bDaggerLeft, bDaggerRight);
        NotifyWritten();
    }

    void LeftMul(const CField* x, UBOOL bDaggerLeft = FALSE, UBOOL bDaggerRight = FALSE) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonVN!"));
            return;
        }
        const CFieldBosonT<deviceDataBoson>* pField = dynamic_cast<const CFieldBosonT<deviceDataBoson>*>(x);
        CCommonKernelField<deviceDataBoson>::LeftMul(m_pDeviceData, m_uiSiteCount, pField->m_pDeviceData, bDaggerLeft, bDaggerRight);
        NotifyWritten();
    }

    void ScalarMultply(const CLGComplex& a) override
    {
        CCommonKernelField<deviceDataBoson>::ScalarMultply(m_pDeviceData, m_uiSiteCount, a);
        NotifyWritten();
    }

    void ScalarMultply(Real a) override
    {
        CCommonKernelField<deviceDataBoson>::ScalarMultply(m_pDeviceData, m_uiSiteCount, a);
        NotifyWritten();
    }

    cuDoubleComplex Dot(const CField* x) const override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonVN!"));
            return make_cuDoubleComplex(0.0, 0.0);
        }
        const CFieldBosonT<deviceDataBoson>* pField = dynamic_cast<const CFieldBosonT<deviceDataBoson>*>(x);
        return CCommonKernelField<deviceDataBoson>::Dot(m_pDeviceData, m_uiSiteCount, pField->m_pDeviceData);
    }

    DOUBLE GetLength() const override
    {
        return CCommonKernelField<deviceDataBoson>::LengthSq(m_pDeviceData, m_uiSiteCount);
    }

    TArray<DOUBLE> Sum() const override
    {
        return CCommonKernelField<deviceDataBoson>::Sum(m_pDeviceData, m_uiSiteCount);
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

    void InitialAsSource(const SFermionBosonSource& sourceData) override
    {
        CCommonKernelSite<deviceDataBoson>::InitialSource(m_pDeviceData, m_byFieldId, sourceData);
        NotifyWritten();
    }

    deviceDataBoson* m_pDeviceData;

    //Improve-1: this field's own buffer handle (pool copies each bind theirs).
    CHaloBufferHandle m_HaloBuffer;
    CHaloBufferHandle* GetHaloBufferHandle() override { return &m_HaloBuffer; }
    const CHaloBufferHandle* GetHaloBufferHandle() const override { return &m_HaloBuffer; }

    _GetData

protected:


    /**
    * f(x) phi^*(x)phi(x)
    */
    virtual void DiagnalTerm(const deviceDataBoson* pSource, DOUBLE fCoeffiecient, _deviceCoeffFunctionPointer fpCoeff, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff)
    {
        CCommonKernelSite<deviceDataBoson>::DiagnalTerm(m_pDeviceData, m_byFieldId, pSource, fCoeffiecient, fpCoeff, eOCT, fRealCoeff, cCmpCoeff);
    }

};


__END_NAMESPACE

#endif //#ifndef _CFIELDBOSONVN_H_

//=============================================================================
// END OF FILE
//=============================================================================