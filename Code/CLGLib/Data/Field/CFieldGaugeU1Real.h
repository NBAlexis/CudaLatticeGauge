//=============================================================================
// FILENAME : CFieldGaugeU1.h
// 
// DESCRIPTION:
// This is the class for the gauge fields
//
// REVISION:
//  [10/13/2020 nbale]
//=============================================================================

#ifndef _CFIELDGAUGE_U1_REAL_H_
#define _CFIELDGAUGE_U1_REAL_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeU1Real)

class CLGAPI CFieldGaugeU1Real : public CFieldGauge
{
    __CLGDECLARE_FIELD(CFieldGaugeU1Real)

public:
    CFieldGaugeU1Real();
    ~CFieldGaugeU1Real();

    void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFileType) override;
    void InitialWithByte(BYTE* byData) override;
    void InitialWithByteCompressed(BYTE*) override;
    void InitialField(EFieldInitialType eInitialType) override;
    EFieldType GetFieldType() const override { return EFT_GaugeReal; }
    void DebugPrintMe() const override;

#pragma region HMC

    void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override;
    void CalculateOnlyStaple(CFieldGauge* pStaple) const override;
    void MakeRandomGenerator() override;
#if !_CLG_DOUBLEFLOAT
    DOUBLE CalculatePlaqutteEnergy(DOUBLE betaOverN) const override;
    DOUBLE CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const override;
    DOUBLE CalculatePlaqutteEnergyUsingStable(DOUBLE betaOverN, const CFieldGauge* pStaple) const override;
    DOUBLE CalculateKinematicEnergy() const override;
#else
    Real CalculatePlaqutteEnergy(Real betaOverN) const override;
    Real CalculatePlaqutteEnergyUseClover(Real betaOverN) const override;
    Real CalculatePlaqutteEnergyUsingStable(Real betaOverN, const CFieldGauge *pStaple) const override;
    Real CalculateKinematicEnergy() const override;
#endif

#pragma endregion

#pragma region BLAS

    void Zero() override;
    void Identity() override;
    void Dagger() override;

    void AxpyPlus(const CField* x) override;
    void AxpyMinus(const CField* x) override;
    void Axpy(Real a, const CField* x) override;
    void Axpy(const CLGComplex& a, const CField* x) override;
    void ScalarMultply(const CLGComplex& a) override;
    void ScalarMultply(Real a) override;

    void SetOneDirectionUnity(BYTE byDir) override;
    void SetOneDirectionZero(BYTE byDir) override;

#pragma endregion

#pragma region Test Functions to test gauge invarience of angular momentum

    /**
     * iA = U.TA() / 2
     */
    void TransformToIA() override;

    /**
     * U=exp(iA)
     */
    void TransformToU() override;

    void CalculateE_Using_U(CFieldGauge* pResoult) const override;

    void CalculateNablaE_Using_U(CFieldGauge* pResoult, UBOOL bNaive = FALSE) const override;

#pragma endregion

    void ExpMult(Real a, CField* U) const override;

    //No need to normalize
    void ElementNormalize() override { ; }
#if !_CLG_DOUBLEFLOAT
    cuDoubleComplex Dot(const CField* other) const override;
#else
    CLGComplex Dot(const CField* other) const override;
#endif
    BYTE* CopyDataOut(UINT &uiSize) const override;
    BYTE* CopyDataOutFloat(UINT& uiSize) const override;
    BYTE* CopyDataOutDouble(UINT& uiSize) const override;
    CCString GetInfos(const CCString &tab) const override;

    Real* m_pDeviceData;
    EFieldInitialType m_eInitialType;

protected:

    void SetByArray(Real* array);
};

#pragma region Helper device functions

//static __device__ __inline__ CLGComplex _deviceLinkU1Real(
//    const Real* __restrict__ pDeviceData,
//    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
//    const INT* __restrict__ pDir)
//{
//    //length can be 0
//    CLGComplex sRet = _onec;
//    for (BYTE i = 0; i < byLength; ++i)
//    {
//        if (0 == pDir[i])
//        {
//            continue;
//        }
//        UBOOL bDagger = FALSE;
//        const BYTE byDir = pDir[i] > 0 ?
//            static_cast<BYTE>(pDir[i] - 1) : static_cast<BYTE>(-pDir[i] - 1);
//
//        if (pDir[i] < 0) //Move
//        {
//            bDagger = TRUE;
//            _deviceSmallInt4Offset(sStartSite, pDir[i]);
//        }
//        const SIndex& newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byDir];
//
//        if (0 == i)
//        {
//            if (!newLink.IsDirichlet())
//            {
//                sRet = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
//                if ((newLink.NeedToDagger() && !bDagger)
//                    || (!newLink.NeedToDagger() && bDagger)
//                    )
//                {
//                    sRet.y = -sRet.y;
//                }
//            }
//        }
//        else
//        {
//            if (!newLink.IsDirichlet())
//            {
//                if ((newLink.NeedToDagger() && !bDagger)
//                 || (!newLink.NeedToDagger() && bDagger)
//                    )
//                {
//                    const CLGComplex& toMul = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
//                    sRet = _make_cuComplex(sRet.x * toMul.x + sRet.y * toMul.y, sRet.y * toMul.x - sRet.x * toMul.y);
//                }
//                else
//                {
//                    const CLGComplex& toMul = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
//                    sRet = _make_cuComplex(sRet.x * toMul.x - sRet.y * toMul.y, sRet.y * toMul.x + sRet.x * toMul.y);
//                }
//            }
//        }
//
//        if (pDir[i] > 0) //Move
//        {
//            _deviceSmallInt4Offset(sStartSite, pDir[i]);
//        }
//    }
//
//    return sRet;
//}
//
//static __device__ __inline__ Real _deviceLinkU1RealArgSum(
//    const Real* __restrict__ pDeviceData,
//    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
//    const INT* __restrict__ pDir)
//{
//    //length can be 0
//    Real sRet = F(0.0);
//    for (BYTE i = 0; i < byLength; ++i)
//    {
//        if (0 == pDir[i])
//        {
//            continue;
//        }
//        UBOOL bDagger = FALSE;
//        const BYTE byDir = pDir[i] > 0 ?
//            static_cast<BYTE>(pDir[i] - 1) : static_cast<BYTE>(-pDir[i] - 1);
//
//        if (pDir[i] < 0) //Move
//        {
//            bDagger = TRUE;
//            _deviceSmallInt4Offset(sStartSite, pDir[i]);
//        }
//        const SIndex& newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byDir];
//
//        if (0 == i)
//        {
//            if (!newLink.IsDirichlet())
//            {
//                sRet = __cuCargf(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
//                if ((newLink.NeedToDagger() && !bDagger)
//                 || (!newLink.NeedToDagger() && bDagger)
//                    )
//                {
//                    sRet = -sRet;
//                }
//            }
//        }
//        else
//        {
//            if (!newLink.IsDirichlet())
//            {
//                if ((newLink.NeedToDagger() && !bDagger)
//                 || (!newLink.NeedToDagger() && bDagger)
//                    )
//                {
//                    sRet -= __cuCargf(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
//                }
//                else
//                {
//                    sRet += __cuCargf(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
//                }
//            }
//        }
//
//        if (pDir[i] > 0) //Move
//        {
//            _deviceSmallInt4Offset(sStartSite, pDir[i]);
//        }
//    }
//
//    return sRet;
//}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_U1_H_

//=============================================================================
// END OF FILE
//=============================================================================