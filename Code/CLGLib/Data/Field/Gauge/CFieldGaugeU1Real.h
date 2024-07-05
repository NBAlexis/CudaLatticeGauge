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

__DEFINE_ENUM(EU1RealType,
    EURT_None,
    EURT_ImagineChemical,
    EURT_E_t,
    EURT_E_z,
    EURT_Bp_x,
    EURT_Bp_y,
    EURT_Bp_xy,
    EURT_Bp_x_notwist,
    EURT_Bp_y_notwist,
    EURT_Bp_xy_notwist,
    );

class CLGAPI CFieldGaugeU1Real : public CFieldGauge
{
    __CLGDECLARE_FIELD(CFieldGaugeU1Real)

public:
    CFieldGaugeU1Real();
    ~CFieldGaugeU1Real();

    void InitialOtherParameters(CParameters& param) override;

    void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFileType) override;
    void InitialWithByte(BYTE* byData) override;
    void InitialWithByteCompressed(const CCString& sFileName) override;
    void InitialField(EFieldInitialType eInitialType) override;
    EFieldType GetFieldType() const override { return EFT_GaugeReal; }
    UINT MatrixN() const override  { return 1; }
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

    void PolyakovOnSpatialSite(cuDoubleComplex* buffer) const override;

    Real* m_pDeviceData;
    EFieldInitialType m_eInitialType;

    _GetData

    /**
     * Note!!!
     * Do not check the conflicts, assume you know
     */
    void InitialU1Real(EU1RealType eChemicalType, EU1RealType eEType, EU1RealType eBType, Real fChemical, Real feEz, Real feBz, UBOOL bXYShiftCenter);

    Real CheckSliceSame(BYTE dir1, BYTE dir2) const;

    /**
    * check r[d1,d2,0,0].link[linkdirs] = 0, for all index of d1 and d2
    */
    Real CheckZero(BYTE dir1, BYTE dir2, const TArray<BYTE>& linkdirs) const;
    void DebugPrintSlice(BYTE dir1, BYTE dir2, const TArray<BYTE>& linkdirs) const;

    EU1RealType m_eChemical;
    EU1RealType m_eE;
    EU1RealType m_eB;
    Real m_fChemical;
    Real m_feEz;
    Real m_feBz;
    UBOOL m_bXYShiftCenter;

protected:

    void SetByArray(Real* array);
};


__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_U1_H_

//=============================================================================
// END OF FILE
//=============================================================================