//=============================================================================
// FILENAME : CUpdator.h
// 
// DESCRIPTION:
// This is the class for update the field
// It can be continous for U(1), SU(2), SU(3) gauge, using hybrid Monte Carlo
// It can also be descrete for Z2, Zn, tetrahydraul, octahydraul, icosahydraul, using Heatbath
//
// REVISION:
//  [mm/dd/yy]
//  [12/8/2018 nbale]
//=============================================================================

#ifndef _CUPDATOR_H_
#define _CUPDATOR_H_

__BEGIN_NAMESPACE

__DEFINE_ENUM(EUpdatorType,
    EUT_HMC,
    EUT_Heatbath,
    EUT_Max,

    EUT_ForceDWORD = 0x7fffffff,
)

class CLGAPI CUpdator : public CBase
{
public:

    CUpdator() 
        : m_pOwner(NULL)
        , m_uiUpdateCall(0)
        , m_iAcceptedConfigurationCount(0)
        , m_iSaveIndexStart(0)
        , m_bSaveConfigurations(FALSE)
        , m_eSaveFieldType(EFFT_CLGBin)
        , m_bTestHDiff(FALSE) 
        , m_bReport(TRUE)
        , m_bAdaptiveUpdate(FALSE)
        , m_uiMinStep(5)
        , m_uiMaxStep(100)
        , m_fGrowStep(F(-0.3))
        , m_fReduceStep(F(0.03))
        , m_fLastHDiff(0.0)
    {
    }

    virtual UINT Update(UINT iSteps, UBOOL bMeasure) = 0;
    virtual void UpdateUntileAccept(UINT iSteps, UBOOL bMeasure);
    virtual Real CalculateEnergy() = 0;

    virtual EUpdatorType GetUpdatorType() const { return EUT_Max; }
    virtual void Initial(class CLatticeData* pOwner, const CParameters& params) = 0;
    virtual CCString GetInfos(const CCString& sTab) const = 0;

    class CLatticeData* m_pOwner;

    void SetConfigurationCount(UINT uiCount) { m_iAcceptedConfigurationCount = uiCount; }
    UINT GetConfigurationCount() const { return m_iAcceptedConfigurationCount; }
    //void OnConfigurationAccepted();
    void SaveConfiguration(UINT uiUpdateStep) const;

    /**
    * Load the tensor2 companion files of a saved configuration into the
    * lattice tensor2 fields (in place).
    *
    * For every dynamic lattice tensor2 field with id F, the companion file
    *   <prefix>_<N>_t<F>.con
    * (written by SaveConfiguration) is loaded with eLoadType. Fields without a
    * matching file are left untouched. This is the counterpart of
    * SaveConfiguration's tensor2 saving and is meant to be called by
    * measurement applications right after loading the gauge .con file.
    */
    void LoadTensor2Configuration(const CCString& sPrefix, UINT uiN, EFieldFileType eLoadType) const;

    void SetTestHdiff(UBOOL bTestHDiff) 
    {
        m_bTestHDiff = bTestHDiff;
    }

    void ClearHDiffHistory()
    {
        m_lstHDiff.RemoveAll();
    }

    virtual void SetAutoCorrection(UBOOL bAutoCorrection) = 0;

#if !_CLG_DOUBLEFLOAT

    Real GetHDiff() const
    {
        if (0 == m_lstHDiff.Num())
        {
            return F(0.0);
        }

        if (1 == m_lstHDiff.Num())
        {
            return appAbs(static_cast<Real>(m_lstHDiff[0]));
        }

        Real fAdd = F(0.0);
        for (INT i = 0; i < m_lstHDiff.Num(); ++i)
        {
            fAdd += static_cast<Real>(m_lstHDiff[i] * m_lstHDiff[i]);
        }
        return _hostsqrt(fAdd / (m_lstHDiff.Num() - 1));
    }

    Real GetHValue() const
    {
        if (0 == m_lstH.Num())
        {
            return F(0.0);
        }

        Real fAdd = F(0.0);
        for (INT i = 0; i < m_lstH.Num(); ++i)
        {
            fAdd += static_cast<Real>(m_lstH[i]);
        }
        return fAdd / m_lstH.Num();
    }

#else
    Real GetHDiff() const
    {
        if (0 == m_lstHDiff.Num())
        {
            return F(0.0);
        }

        if (1 == m_lstHDiff.Num())
        {
            return abs(m_lstHDiff[0]);
        }

        Real fAdd = F(0.0);
        for (INT i = 0; i < m_lstHDiff.Num(); ++i)
        {
            fAdd += (m_lstHDiff[i] * m_lstHDiff[i]);
        }
        return _hostsqrt(fAdd / (m_lstHDiff.Num() - 1));
    }

    Real GetHValue() const
    {
        if (0 == m_lstH.Num())
        {
            return F(0.0);
        }

        Real fAdd = F(0.0);
        for (INT i = 0; i < m_lstH.Num(); ++i)
        {
            fAdd += m_lstH[i];
        }
        return fAdd / m_lstH.Num();
    }
#endif

    Real GetLastHDiff() const
    {
        return static_cast<Real>(m_fLastHDiff);
    }

    Real GetMaxHDiff() const
    {
        DOUBLE fMax = appAbs(m_lstHDiff[0]);
        INT idx = 0;
        for (INT i = 1; i < m_lstHDiff.Num(); ++i)
        {
            if (fMax < appAbs(m_lstHDiff[i]))
            {
                fMax = appAbs(m_lstHDiff[i]);
                idx = i;
            }
        }
        return static_cast<Real>(m_lstHDiff[idx]);
    }

    void SetSaveConfiguration(UBOOL bSave, const CCString& sPrefix, UINT saveIndexStart = 0, EFieldFileType efft = EFFT_CLGBin)
    {
        m_bSaveConfigurations = bSave;
        m_sConfigurationPrefix = sPrefix;
        m_iSaveIndexStart = saveIndexStart;
        m_eSaveFieldType = efft;
    }

protected:

    UINT m_uiUpdateCall;
    UINT m_iAcceptedConfigurationCount;
    UINT m_iSaveIndexStart;
    UBOOL m_bSaveConfigurations;
    EFieldFileType m_eSaveFieldType;
    UBOOL m_bTestHDiff;
    UBOOL m_bReport;
    UBOOL m_bAdaptiveUpdate;
    UINT m_uiMinStep;
    UINT m_uiMaxStep;
    Real m_fGrowStep;
    Real m_fReduceStep;
    CCString m_sConfigurationPrefix;
#if !_CLG_DOUBLEFLOAT
    TArray<DOUBLE> m_lstHDiff;
    TArray<DOUBLE> m_lstH;
    DOUBLE m_fLastHDiff;
#else
    TArray<Real> m_lstHDiff;
    TArray<Real> m_lstH;
    Real m_fLastHDiff;
#endif
};

__END_NAMESPACE

#endif //#ifndef _CUPDATOR_H_

//=============================================================================
// END OF FILE
//=============================================================================