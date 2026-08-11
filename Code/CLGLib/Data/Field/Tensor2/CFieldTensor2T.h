//=============================================================================
// FILENAME : CFieldTensor2T.h
//
// DESCRIPTION:
// This is the template class for tensor2 fields.
// The element count is 6 x volume, stored as data[plaqutteIndex * siteCount + siteIndex]
// (component-major, the same layout as the Fmunu buffers in CStapleCache),
// all the BLAS operators are the common kernels with count = 6 x volume.
//
// REVISION:
//  [mm/dd/yy]
//  [07/24/2026 nbale]
//=============================================================================

#ifndef _CFIELDTENSOR2T_H_
#define _CFIELDTENSOR2T_H_

__BEGIN_NAMESPACE

template<typename deviceDataTensor2>
class __DLL_EXPORT CFieldTensor2T : public CFieldTensor2
{
public:
    typedef CCommonKernelField<deviceDataTensor2> _FieldKernel;
    typedef CFieldTensor2Kernel<deviceDataTensor2> _Tensor2Kernel;

    CFieldTensor2T()
        : CFieldTensor2()
        , m_pDeviceData(NULL)
    {
        checkCudaErrors(__cudaMalloc((void**)& m_pDeviceData, GetElementCount() * sizeof(deviceDataTensor2)));

        //Improve-1 (multi-GPU-improve1.md 3.1/3.2, I3): component-major
        //(plaquette-index-major) storage has no per-site halo slot numbering
        //and no stencil kernel reads neighbours from this buffer -- bind the
        //handle as LocalOnly (halo site count 0, not halo-capable) so the
        //extent is still registry-tracked and versioned.
        SHaloBufferInfo sInfo;
        sInfo.m_pDeviceData = reinterpret_cast<BYTE*>(m_pDeviceData);
        sInfo.m_uiCapacityBytes = GetElementCount() * sizeof(deviceDataTensor2);
        sInfo.m_uiBytesPerSite = PlaqutteCountPerSite() * sizeof(deviceDataTensor2);
        sInfo.m_uiLocalSiteCount = m_uiSiteCount;
        sInfo.m_uiHaloSiteCount = 0;
        sInfo.m_ullLayoutGeneration = appGetLayoutGeneration();
        sInfo.m_byFieldId = 0;
        sInfo.m_bHaloCapable = FALSE;
        m_HaloBuffer.Bind(sInfo);
    }

    ~CFieldTensor2T()
    {
        checkCudaErrors(__cudaFree(m_pDeviceData));
        m_pDeviceData = NULL;
    }

    void CopyBufferTo(CField* U) const override
    {
        if (NULL == U || GetFieldType() != U->GetFieldType())
        {
            appCrucial(_T("CFieldTensor2 can only copy to the same type!"));
            return;
        }
        CFieldTensor2T<deviceDataTensor2>* pField = dynamic_cast<CFieldTensor2T<deviceDataTensor2>*>(U);
        CCommonKernel<deviceDataTensor2>::CopyBuffer(pField->m_pDeviceData, m_pDeviceData, GetElementCount());
        U->NotifyWritten();
    }

    void InitialField(EFieldInitialType eInitialType) override
    {
        _Tensor2Kernel::Initial(m_pDeviceData, m_uiSiteCount, PlaqutteCountPerSite(), eInitialType);
        NotifyWritten();
    }

    void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFieldType) override
    {
        if (eFieldType != EFFT_CLGBin
            && eFieldType != EFFT_CLGBinFloat
            && eFieldType != EFFT_CLGBinDouble)
        {
            appCrucial(_T("CFieldTensor2::InitialFieldWithFile: Not support %s File\n"), __ENUM_TO_STRING(EFieldFileType, eFieldType).c_str());
            return;
        }

        const UINT uiFloatCount = FloatN() * m_uiSiteCount;
        UINT uiSize = static_cast<UINT>(sizeof(Real) * uiFloatCount);
        if (eFieldType == EFFT_CLGBinFloat)
        {
            uiSize = static_cast<UINT>(sizeof(FLOAT) * uiFloatCount);
        }
        else if (eFieldType == EFFT_CLGBinDouble)
        {
            uiSize = static_cast<UINT>(sizeof(DOUBLE) * uiFloatCount);
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
            DOUBLE* data2 = (DOUBLE*)(malloc(sizeof(DOUBLE) * uiFloatCount));
            for (UINT i = 0; i < uiFloatCount; ++i)
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
            FLOAT* data2 = (FLOAT*)(malloc(sizeof(FLOAT) * uiFloatCount));
            for (UINT i = 0; i < uiFloatCount; ++i)
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
        _FieldKernel::InitialWithByte(m_pDeviceData, GetElementCount(), byData);
        NotifyWritten();
    }

    void DebugPrintMe() const override
    {
        UINT uiSize = 0;
        BYTE* pHostData = CopyDataOut(uiSize);
        const Real* pRealData = (const Real*)pHostData;
        CCString sOut;
        const UINT uiFloatCount = uiSize / sizeof(Real);
        const UINT uiPrintCount = uiFloatCount < 16 ? uiFloatCount : 16;
        for (UINT i = 0; i < uiPrintCount; ++i)
        {
            sOut = sOut + appToString(pRealData[i]) + _T(" ");
        }
        sOut = sOut + _T("\n");
        appGeneral(_T("%s"), sOut.c_str());
        free(pHostData);
    }

    void Dagger() override
    {
        _FieldKernel::Dagger(m_pDeviceData, GetElementCount());
        NotifyWritten();
    }

    void AxpyPlus(const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldTensor2 can only work with the same type!"));
            return;
        }
        const CFieldTensor2T<deviceDataTensor2>* pField = dynamic_cast<const CFieldTensor2T<deviceDataTensor2>*>(x);
        _FieldKernel::AxpyPlus(m_pDeviceData, GetElementCount(), pField->m_pDeviceData);
        NotifyWritten();
    }

    void AxpyMinus(const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldTensor2 can only work with the same type!"));
            return;
        }
        const CFieldTensor2T<deviceDataTensor2>* pField = dynamic_cast<const CFieldTensor2T<deviceDataTensor2>*>(x);
        _FieldKernel::AxpyMinus(m_pDeviceData, GetElementCount(), pField->m_pDeviceData);
        NotifyWritten();
    }

    void Axpy(Real a, const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldTensor2 can only work with the same type!"));
            return;
        }
        const CFieldTensor2T<deviceDataTensor2>* pField = dynamic_cast<const CFieldTensor2T<deviceDataTensor2>*>(x);
        _FieldKernel::Axpy(m_pDeviceData, GetElementCount(), a, pField->m_pDeviceData);
        NotifyWritten();
    }

    void Axpy(const CLGComplex& a, const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldTensor2 can only work with the same type!"));
            return;
        }
        const CFieldTensor2T<deviceDataTensor2>* pField = dynamic_cast<const CFieldTensor2T<deviceDataTensor2>*>(x);
        _FieldKernel::Axpy(m_pDeviceData, GetElementCount(), a, pField->m_pDeviceData);
        NotifyWritten();
    }

    void Mul(const CField* x, UBOOL bDaggerLeft = TRUE, UBOOL bDaggerRight = FALSE) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldTensor2 can only work with the same type!"));
            return;
        }
        const CFieldTensor2T<deviceDataTensor2>* pField = dynamic_cast<const CFieldTensor2T<deviceDataTensor2>*>(x);
        _FieldKernel::Mul(m_pDeviceData, GetElementCount(), pField->m_pDeviceData, bDaggerLeft, bDaggerRight);
        NotifyWritten();
    }

    void LeftMul(const CField* x, UBOOL bDaggerLeft = FALSE, UBOOL bDaggerRight = FALSE) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldTensor2 can only work with the same type!"));
            return;
        }
        const CFieldTensor2T<deviceDataTensor2>* pField = dynamic_cast<const CFieldTensor2T<deviceDataTensor2>*>(x);
        _FieldKernel::LeftMul(m_pDeviceData, GetElementCount(), pField->m_pDeviceData, bDaggerLeft, bDaggerRight);
        NotifyWritten();
    }

    void ScalarMultply(const CLGComplex& a) override
    {
        _FieldKernel::ScalarMultply(m_pDeviceData, GetElementCount(), a);
        NotifyWritten();
    }

    void ScalarMultply(Real a) override
    {
        _FieldKernel::ScalarMultply(m_pDeviceData, GetElementCount(), a);
        NotifyWritten();
    }

    cuDoubleComplex Dot(const CField* x) const override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial(_T("CFieldTensor2 can only work with the same type!"));
            return make_cuDoubleComplex(0.0, 0.0);
        }
        const CFieldTensor2T<deviceDataTensor2>* pField = dynamic_cast<const CFieldTensor2T<deviceDataTensor2>*>(x);
        //the layout is component-major, data[plaqutteIndex * siteCount + siteIndex],
        //so every segment is a common site field of count Volume
        cuDoubleComplex res = make_cuDoubleComplex(0.0, 0.0);
        for (UINT i = 0; i < PlaqutteCountPerSite(); ++i)
        {
            const cuDoubleComplex cSeg = _FieldKernel::Dot(m_pDeviceData + i * m_uiSiteCount, m_uiSiteCount, pField->m_pDeviceData + i * m_uiSiteCount);
            res = make_cuDoubleComplex(res.x + cSeg.x, res.y + cSeg.y);
        }
        return res;
    }

    DOUBLE GetLength() const override
    {
        DOUBLE res = 0.0;
        for (UINT i = 0; i < PlaqutteCountPerSite(); ++i)
        {
            res += _FieldKernel::LengthSq(m_pDeviceData + i * m_uiSiteCount, m_uiSiteCount);
        }
        return res;
    }

    TArray<DOUBLE> Sum() const
    {
        TArray<DOUBLE> ret = _FieldKernel::Sum(m_pDeviceData, m_uiSiteCount);
        for (UINT i = 1; i < PlaqutteCountPerSite(); ++i)
        {
            const TArray<DOUBLE> seg = _FieldKernel::Sum(m_pDeviceData + i * m_uiSiteCount, m_uiSiteCount);
            for (INT j = 0; j < ret.Num(); ++j)
            {
                ret[j] += seg[j];
            }
        }
        return ret;
    }

    BYTE* CopyDataOut(UINT& uiSize) const override
    {
        return _FieldKernel::CopyDataOut(m_pDeviceData, GetElementCount(), uiSize);
    }

    BYTE* CopyDataOutFloat(UINT& uiSize) const override
    {
        return _FieldKernel::CopyDataOutFloat(m_pDeviceData, GetElementCount(), uiSize);
    }

    BYTE* CopyDataOutDouble(UINT& uiSize) const override
    {
        return _FieldKernel::CopyDataOutDouble(m_pDeviceData, GetElementCount(), uiSize);
    }

    deviceDataTensor2* m_pDeviceData;

    //Improve-1: this field's own buffer handle (LocalOnly, pool copies each
    //bind theirs).
    CHaloBufferHandle m_HaloBuffer;
    CHaloBufferHandle* GetHaloBufferHandle() override { return &m_HaloBuffer; }
    const CHaloBufferHandle* GetHaloBufferHandle() const override { return &m_HaloBuffer; }

    _GetData

};


__END_NAMESPACE

#endif //#ifndef _CFIELDTENSOR2T_H_

//=============================================================================
// END OF FILE
//=============================================================================
