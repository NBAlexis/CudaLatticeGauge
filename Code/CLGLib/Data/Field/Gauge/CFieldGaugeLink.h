//=============================================================================
// FILENAME : CFieldGaugeLink.h
// 
// DESCRIPTION:
// This is the common class for all gauge fields
//
// REVISION:
//  [07/04/2018 nbale]
//=============================================================================
#include "CFieldGaugeKernel.h"

#ifndef _CFIELDGAUGE_LINK_H_
#define _CFIELDGAUGE_LINK_H_

#define __DEFINE_GAUGE_LINK(CLASSNAME, DEVICEDATA, N, FIELDTYPE) \
__CLG_REGISTER_HELPER_HEADER(CLASSNAME) \
class CLGAPI CLASSNAME : public CFieldGaugeLink<DEVICEDATA, N> \
{ \
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CLASSNAME) \
public: \
    EFieldType GetFieldType() const override { return FIELDTYPE; } \
};


__BEGIN_NAMESPACE

template<typename deviceGauge, INT matrixN>
class __DLL_EXPORT CFieldGaugeLink : public CFieldGauge
{

public:
    CFieldGaugeLink() : CFieldGauge()
    {
        CCommonKernel<deviceGauge>::AllocateBuffer(&m_pDeviceData, m_uiLinkeCount);
    }

    ~CFieldGaugeLink()
    {
        CCommonKernel<deviceGauge>::FreeBuffer(&m_pDeviceData);
    }

    void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFileType) override
    {
        if (!CFileSystem::IsFileExist(sFileName))
        {
            appCrucial(_T("File not exist!!! %s \n"), sFileName.c_str());
            _FAIL_EXIT;
        }

        switch (eFileType)
        {
        case EFFT_CLGBin:
#if _CLG_DOUBLEFLOAT
        case EFFT_CLGBinDouble:
#else
        case EFFT_CLGBinFloat:
#endif
        {
            UINT uiSize = static_cast<UINT>(sizeof(Real) * 2 * MatrixN() * MatrixN() * m_uiLinkeCount);
            BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
            if (uiSize != sizeof(Real) * 2 * MatrixN() * MatrixN() * _HC_LinkCount)
            {
                appCrucial(_T("Loading file size not match: %s, %d, expecting %d"), sFileName.c_str(), uiSize, static_cast<INT>(sizeof(Real) * 2 * MatrixN() * MatrixN() * _HC_LinkCount));
            }
            InitialWithByte(data);
            free(data);
            FixBoundary();
        }
        break;
#if _CLG_DOUBLEFLOAT
        case EFFT_CLGBinFloat:
        {
            UINT uiSize = static_cast<UINT>(sizeof(Real) * 2 * MatrixN() * MatrixN() * m_uiLinkeCount);
            BYTE* data = (BYTE*)malloc(uiSize);
            Real* rdata = (Real*)data;
            FLOAT* fdata = (FLOAT*)appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
            if (uiSize != sizeof(FLOAT) * 2 * MatrixN() * MatrixN() * _HC_LinkCount)
            {
                appCrucial(_T("Loading file size not match: %s, %d, expecting %d"), sFileName.c_str(), uiSize, static_cast<UINT>(sizeof(FLOAT) * 2 * MatrixN() * MatrixN() * _HC_LinkCount));
            }
            for (UINT i = 0; i < 2 * MatrixN() * MatrixN() * m_uiLinkeCount; ++i)
            {
                rdata[i] = static_cast<Real>(fdata[i]);
            }
            InitialWithByte(data);
            free(fdata);
            free(data);
            FixBoundary();
        }
        break;
#else
        case EFFT_CLGBinDouble:
        {
            UINT uiSize = static_cast<UINT>(sizeof(Real) * 2 * MatrixN() * MatrixN() * m_uiLinkeCount);
            BYTE* data = (BYTE*)malloc(uiSize);
            Real* rdata = (Real*)data;
            DOUBLE* ddata = (DOUBLE*)appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
            if (uiSize != sizeof(DOUBLE) * 2 * MatrixN() * MatrixN() * _HC_LinkCount)
            {
                appCrucial(_T("Loading file size not match: %s, %d, expecting %d"), sFileName.c_str(), uiSize, static_cast<INT>(sizeof(DOUBLE) * 2 * MatrixN() * MatrixN() * _HC_LinkCount));
            }
            for (UINT i = 0; i < 2 * MatrixN() * MatrixN() * m_uiLinkeCount; ++i)
            {
                rdata[i] = static_cast<Real>(ddata[i]);
            }
            InitialWithByte(data);
            free(ddata);
            free(data);
            FixBoundary();
        }
        break;
#endif
        case EFFT_CLGBinCompressed:
        {
            InitialWithByteCompressed(sFileName);
            FixBoundary();
        }
        break;
        default:
            appCrucial(_T("Not supported input file type %s\n"), __ENUM_TO_STRING(EFieldFileType, eFileType).c_str());
            break;

        }
    }

    void InitialWithByte(BYTE* byData) override
    {
        CCommonKernelField<deviceGauge>::InitialWithByte(m_pDeviceData, m_uiLinkeCount, byData);
    }

    void InitialWithByteCompressed(const CCString& fileName) override { appCrucial(_T("CFieldGaugeLink: InitialWithByteCompressed not supoorted!\n")); }

    void InitialField(EFieldInitialType eInitialType) override
    {
        CCommonKernelLink<deviceGauge>::InitialBuffer(m_pDeviceData, m_byFieldId, eInitialType);
    }

    void DebugPrintMe() const override
    {
        CCommonKernelLink<deviceGauge>::DebugPrint(m_pDeviceData, m_uiLinkeCount);
    }

#pragma region HMC

    void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override
    {
        if (NULL == pForce || GetFieldType() != pForce->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: force field is not SU3");
            return;
        }
        if (NULL != pStaple && GetFieldType() != pStaple->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: stape field is not SU3");
            return;
        }

        CFieldGaugeLink<deviceGauge, matrixN>* pForceSU3 = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pForce);
        CFieldGaugeLink<deviceGauge, matrixN>* pStableSU3 = NULL == pStaple ? NULL : dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pStaple);

        CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceAndStaple(
            m_pDeviceData,
            m_byFieldId,
            pForceSU3->m_pDeviceData,
            NULL == pStableSU3 ? NULL : pStableSU3->m_pDeviceData,
            betaOverN);
    }

    void CalculateOnlyStaple(CFieldGauge* pStaple) const override
    {
        if (NULL == pStaple || GetFieldType() != pStaple->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: stable field is not SU3");
            return;
        }
        CFieldGaugeLink<deviceGauge, matrixN>* pStableSU3 = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pStaple);
        CFieldGaugeKernel<deviceGauge, matrixN>::CalculateOnlyStaple(m_pDeviceData, m_byFieldId, pStableSU3->m_pDeviceData);
    }

    void MakeRandomGenerator() override
    {
        InitialField(EFIT_RandomGenerator);
    }

    DOUBLE CalculatePlaqutteEnergy(DOUBLE betaOverN) const override
    {
        return CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergy(m_pDeviceData, m_byFieldId, betaOverN);
    }

    DOUBLE CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const override
    {
        return CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergyUseClover(m_pDeviceData, m_byFieldId, betaOverN);
    }

    DOUBLE CalculatePlaqutteEnergyUsingStable(DOUBLE betaOverN, const CFieldGauge* pStaple) const override
    {
        if (NULL == pStaple || GetFieldType() != pStaple->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: stape field is not SU3");
            return F(0.0);
        }
        const CFieldGaugeLink<deviceGauge, matrixN>* pStableSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(pStaple);

        return CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergyUsingStable(m_pDeviceData, m_byFieldId, betaOverN, pStableSU3->m_pDeviceData);
    }

    DOUBLE CalculateKinematicEnergy() const override
    {
        return CCommonKernelField<deviceGauge>::LengthSq(m_pDeviceData, m_uiLinkeCount);
    }

#pragma endregion

#pragma region BLAS

    void Dagger() override
    {
        CCommonKernelField<deviceGauge>::Dagger(m_pDeviceData, m_uiLinkeCount);
    }

    void AxpyPlus(const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: axpy failed because the otherfield is not SU3");
            return;
        }
        const CFieldGaugeLink<deviceGauge, matrixN>* pSU3x = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(x);
        CCommonKernelField<deviceGauge>::AxpyPlus(m_pDeviceData, m_uiLinkeCount, pSU3x->m_pDeviceData);
    }

    void AxpyMinus(const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: axpy failed because the otherfield is not SU3");
            return;
        }
        const CFieldGaugeLink<deviceGauge, matrixN>* pSU3x = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(x);
        CCommonKernelField<deviceGauge>::AxpyMinus(m_pDeviceData, m_uiLinkeCount, pSU3x->m_pDeviceData);
    }

    void Axpy(Real a, const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: axpy failed because the otherfield is not SU3");
            return;
        }
        const CFieldGaugeLink<deviceGauge, matrixN>* pSU3x = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(x);
        CCommonKernelField<deviceGauge>::Axpy(m_pDeviceData, m_uiLinkeCount, a, pSU3x->m_pDeviceData);
    }

    void Axpy(const CLGComplex& a, const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: axpy failed because the otherfield is not SU3");
            return;
        }
        const CFieldGaugeLink<deviceGauge, matrixN>* pSU3x = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(x);
        CCommonKernelField<deviceGauge>::Axpy(m_pDeviceData, m_uiLinkeCount, a, pSU3x->m_pDeviceData);
    }

    void Mul(const CField* x, UBOOL bDagger = TRUE) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: axpy failed because the otherfield is not SU3");
            return;
        }
        const CFieldGaugeLink<deviceGauge, matrixN>* pSU3x = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(x);
        CCommonKernelField<deviceGauge>::Mul(m_pDeviceData, m_uiLinkeCount, pSU3x->m_pDeviceData, bDagger);
    }

    void ScalarMultply(const CLGComplex& a) override
    {
        CCommonKernelField<deviceGauge>::ScalarMultply(m_pDeviceData, m_uiLinkeCount, a);
    }

    void ScalarMultply(Real a) override
    {
        CCommonKernelField<deviceGauge>::ScalarMultply(m_pDeviceData, m_uiLinkeCount, a);
    }

    void SetOneDirectionUnity(BYTE byDir) override
    {
        if (0 == (byDir & 15))
        {
            return;
        }
        CCommonKernelLink<deviceGauge>::SetOneDirectionUnity(m_pDeviceData, m_byFieldId, byDir);
    }

    void SetOneDirectionZero(BYTE byDir) override
    {
        if (0 == (byDir & 15))
        {
            return;
        }
        CCommonKernelLink<deviceGauge>::SetOneDirectionZero(m_pDeviceData, m_byFieldId, byDir);
    }

#pragma endregion

#pragma region Test Functions to test gauge invarience of angular momentum

    /**
     * iA = U.TA() / 2
     */
    void TransformToIA() override
    {
        if (0 == _HC_ALog)
        {
            CCommonKernelLink<deviceGauge>::StrictLog(m_pDeviceData, m_byFieldId);
        }
        else
        {
            CCommonKernelLink<deviceGauge>::QuickLog(m_pDeviceData, m_byFieldId);
        }
    }

    /**
     * U=exp(iA)
     */
    void TransformToU() override
    {
        if (0 == _HC_ALog)
        {
            CCommonKernelLink<deviceGauge>::StrictExp(m_pDeviceData, m_byFieldId);
        }
        else
        {
            CCommonKernelLink<deviceGauge>::QuickExp(m_pDeviceData, m_byFieldId);
        }
    }

    void CalculateE_Using_U(CFieldGauge* pResoult) const override
    {
        if (NULL == pResoult || GetFieldType() != pResoult->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: U field is not SU3");
            return;
        }

        CFieldGaugeLink<deviceGauge, matrixN>* pUField = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pResoult);
        CCommonKernelLink<deviceGauge>::CalculateE_Using_U(m_pDeviceData, m_byFieldId, pUField->m_pDeviceData);
    }

    void CalculateNablaE_Using_U(CFieldGauge* pResoult, UBOOL bNaive = FALSE) const override
    {
        if (NULL == pResoult || GetFieldType() != pResoult->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: U field is not SU3");
            return;
        }

        CFieldGaugeLink<deviceGauge, matrixN>* pUField = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pResoult);
        CCommonKernelLink<deviceGauge>::CalculateNablaE_Using_U(m_pDeviceData, m_byFieldId, pUField->m_pDeviceData, bNaive);
    }

#pragma endregion

    void ExpMult(Real a, CField* U) const override
    {
        if (NULL == U || GetFieldType() != U->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: U field is not SU3");
            return;
        }

        CFieldGaugeLink<deviceGauge, matrixN>* pUField = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(U);
        CCommonKernelLink<deviceGauge>::ExpMul(pUField->m_pDeviceData, m_byFieldId, m_pDeviceData, a);
    }

    void ElementNormalize() override
    {
        CCommonKernelField<deviceGauge>::Norm(m_pDeviceData, m_uiLinkeCount);
    }

    cuDoubleComplex Dot(const CField* other) const override
    {
        if (NULL == other || GetFieldType() != other->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: U field is not SUN");
            return make_cuDoubleComplex(0, 0);
        }

        const CFieldGaugeLink<deviceGauge, matrixN>* pUField = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(other);
        return CCommonKernelField<deviceGauge>::Dot(m_pDeviceData, m_uiLinkeCount, pUField->m_pDeviceData);
    }

    CCString SaveToCompressedFile(const CCString& fileName) const override { appCrucial(_T("CFieldGaugeLink: SaveToCompressedFile not supoorted!\n")); return _T(""); };

    BYTE* CopyDataOut(UINT& uiSize) const override
    {
        return CCommonKernelField<deviceGauge>::CopyDataOut(m_pDeviceData, m_uiLinkeCount, uiSize);
    }

    BYTE* CopyDataOutFloat(UINT& uiSize) const override
    {
        return CCommonKernelField<deviceGauge>::CopyDataOutFloat(m_pDeviceData, m_uiLinkeCount, uiSize);
    }

    BYTE* CopyDataOutDouble(UINT& uiSize) const override
    {
        return CCommonKernelField<deviceGauge>::CopyDataOutDouble(m_pDeviceData, m_uiLinkeCount, uiSize);
    }

    void CopyTo(CField* pTarget) const override
    {
        if (NULL == pTarget || GetFieldType() != pTarget->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: target field is not SUN");
            return;
        }

        CFieldGauge::CopyTo(pTarget);

        CFieldGaugeLink<deviceGauge, matrixN>* pTargetField = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pTarget);
        CCommonKernel<deviceGauge>::CopyBuffer(pTargetField->m_pDeviceData, m_pDeviceData, m_uiLinkeCount);
    }

    void PolyakovOnSpatialSite(cuDoubleComplex* buffer) const override
    {
        CCommonKernelLink<deviceGauge>::PolyakovOnSpatialSite(m_pDeviceData, m_byFieldId, buffer);
    }

    UINT MatrixN() const override { return matrixN; }

    deviceGauge* m_pDeviceData;

    _GetData

};

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeU1)

class CLGAPI CFieldGaugeU1 : public CFieldGaugeLink<CLGComplex, 1>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldGaugeU1)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeU1; }

    void InitialWithByteCompressed(const CCString& sFileName) override;
    CCString SaveToCompressedFile(const CCString& fileName) const override;
};

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU2)

class CLGAPI CFieldGaugeSU2 : public CFieldGaugeLink<deviceSU2, 2>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldGaugeSU2)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU2; }

    void InitialWithByteCompressed(const CCString& sFileName) override;
    CCString SaveToCompressedFile(const CCString& fileName) const override;
};

__DEFINE_GAUGE_LINK(CFieldGaugeSU4, deviceSU4, 4, EFT_GaugeSU4)
__DEFINE_GAUGE_LINK(CFieldGaugeSU5, deviceSU5, 5, EFT_GaugeSU5)
__DEFINE_GAUGE_LINK(CFieldGaugeSU6, deviceSU6, 6, EFT_GaugeSU6)
__DEFINE_GAUGE_LINK(CFieldGaugeSU7, deviceSU7, 7, EFT_GaugeSU7)
__DEFINE_GAUGE_LINK(CFieldGaugeSU8, deviceSU8, 8, EFT_GaugeSU8)

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_LINK_H_

//=============================================================================
// END OF FILE
//=============================================================================