//=============================================================================
// FILENAME : CMeasurePolyakovSUN.h
// 
// DESCRIPTION:
// This is measurement for Polyakov loop and static quark potential
// Suitable for translational invarient case
//
// REVISION:
//  [07/02/2024 nbale]
//=============================================================================

#ifndef _CMEASUREPOLYAKOVSUN_H_
#define _CMEASUREPOLYAKOVSUN_H_

__BEGIN_NAMESPACE


template <INT N, INT NoE>
class __DLL_EXPORT CMeasurePolyakovSUN : public CMeasure
{
public:

    CMeasurePolyakovSUN()
        : CMeasure()
        , m_pTmpLoop(NULL)
        , m_pTraceRes(NULL)
        , m_pTmpDeviceSum(NULL)

        , m_pCorrelator(NULL)
        , m_pCorrelatorCounter(NULL)
        , m_pHostCorrelator(NULL)
        , m_pHostCorrelatorCounter(NULL)
        , m_uiMaxLengthSq(1)
    {

    }

    ~CMeasurePolyakovSUN();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAcceptedSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    deviceSUN<N, NoE>* m_pTmpLoop;
    CLGComplex* m_pTraceRes;
    CLGComplex* m_pTmpDeviceSum;
    
    CLGComplex* m_pCorrelator;
    UINT* m_pCorrelatorCounter;
    CLGComplex* m_pHostCorrelator;
    UINT* m_pHostCorrelatorCounter;

    UINT m_uiMaxLengthSq;

public:

    TArray<UINT> m_lstR;
    TArray<CLGComplex> m_lstC;
    TArray<CLGComplex> m_lstAverageC;
};

__CLG_REGISTER_HELPER_HEADER(CMeasurePolyakovSU4)
__CLG_REGISTER_HELPER_HEADER(CMeasurePolyakovSU5)
__CLG_REGISTER_HELPER_HEADER(CMeasurePolyakovSU6)
__CLG_REGISTER_HELPER_HEADER(CMeasurePolyakovSU7)
__CLG_REGISTER_HELPER_HEADER(CMeasurePolyakovSU8)

class CLGAPI CMeasurePolyakovSU4 : public CMeasurePolyakovSUN<4, 16>
{
    __CLGDECLARE_CLASS(CMeasurePolyakovSU4)
};
class CLGAPI CMeasurePolyakovSU5 : public CMeasurePolyakovSUN<5, 32>
{
    __CLGDECLARE_CLASS(CMeasurePolyakovSU5)
};
class CLGAPI CMeasurePolyakovSU6 : public CMeasurePolyakovSUN<6, 64>
{
    __CLGDECLARE_CLASS(CMeasurePolyakovSU6)
};
class CLGAPI CMeasurePolyakovSU7 : public CMeasurePolyakovSUN<7, 64>
{
    __CLGDECLARE_CLASS(CMeasurePolyakovSU7)
};
class CLGAPI CMeasurePolyakovSU8 : public CMeasurePolyakovSUN<8, 64>
{
    __CLGDECLARE_CLASS(CMeasurePolyakovSU8)
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREPOLYAKOVSUN_H_

//=============================================================================
// END OF FILE
//=============================================================================