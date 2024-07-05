//=============================================================================
// FILENAME : CFieldGaugeLink.h
// 
// DESCRIPTION:
// This is the common class for all gauge fields
//
// REVISION:
//  [07/04/2018 nbale]
//=============================================================================

#ifndef _CFIELDGAUGE_LINK_H_
#define _CFIELDGAUGE_LINK_H_

#define gaugeLinkKernelFuncionStart \
    intokernaldir; \
    for (UINT idir = 0; idir < uiDir; ++idir) \
    { \
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, idir); 


#define gaugeLinkKernelFuncionEnd \
    } 



__BEGIN_NAMESPACE

template<typename deviceGauge, INT matrixN>
class __DLL_EXPORT CFieldGaugeLink : public CFieldGauge
{

public:
    CFieldGaugeLink();
    ~CFieldGaugeLink();

    void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFileType) override;
    void InitialWithByte(BYTE* byData) override;
    void InitialWithByteCompressed(const CCString& fileName) override { appCrucial(_T("CFieldGaugeLink: InitialWithByteCompressed not supoorted!\n")); }
    void InitialField(EFieldInitialType eInitialType) override;

    void DebugPrintMe() const override;

#pragma region HMC

    void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override;
    void CalculateOnlyStaple(CFieldGauge* pStaple) const override;
    void MakeRandomGenerator() override;
    DOUBLE CalculatePlaqutteEnergy(DOUBLE betaOverN) const override;
    DOUBLE CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const override;
    DOUBLE CalculatePlaqutteEnergyUsingStable(DOUBLE betaOverN, const CFieldGauge* pStaple) const override;
    DOUBLE CalculateKinematicEnergy() const override;

#pragma endregion

#pragma region BLAS

    void Zero() override;
    void Identity() override;
    void Dagger() override;

    void AxpyPlus(const CField* x) override;
    void AxpyMinus(const CField* x) override;
    void Axpy(Real a, const CField* x) override;
    void Axpy(const CLGComplex& a, const CField* x) override;
    void Mul(const CField* other, UBOOL bDagger = TRUE) override;
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

    void ElementNormalize() override;
    cuDoubleComplex Dot(const CField* other) const override;
    CCString SaveToCompressedFile(const CCString& fileName) const override { appCrucial(_T("CFieldGaugeLink: SaveToCompressedFile not supoorted!\n")); return _T(""); };
    BYTE* CopyDataOut(UINT& uiSize) const override;
    BYTE* CopyDataOutFloat(UINT& uiSize) const override;
    BYTE* CopyDataOutDouble(UINT& uiSize) const override;

    void CopyTo(CField* pTarget) const override;

    void PolyakovOnSpatialSite(cuDoubleComplex* buffer) const override;

    UINT MatrixN() const override { return matrixN; }

    deviceGauge* m_pDeviceData;

    _GetData

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

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_LINK_H_

//=============================================================================
// END OF FILE
//=============================================================================