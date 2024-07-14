//=============================================================================
// FILENAME : CMeasureBerryPhase.h
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [10/13/2020 nbale]
//=============================================================================
#include "Data/Field/Gauge/CFieldGaugeLink.h"

#ifndef _CMEASUREBERRYPHASE_H_
#define _CMEASUREBERRYPHASE_H_

__BEGIN_NAMESPACE

static __device__ __inline__ Real _deviceLinkU1ArgSum(
    const CLGComplex* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    Real sRet = F(0.0);
    for (BYTE i = 0; i < byLength; ++i)
    {
        if (0 == pDir[i])
        {
            continue;
        }
        UBOOL bDagger = FALSE;
        const BYTE byDir = pDir[i] > 0 ?
            static_cast<BYTE>(pDir[i] - 1) : static_cast<BYTE>(-pDir[i] - 1);

        if (pDir[i] < 0) //Move
        {
            bDagger = TRUE;
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
        const SIndex& newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byDir];

        if (0 == i)
        {
            if (!newLink.IsDirichlet())
            {
                sRet = __cuCargf(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    sRet = -sRet;
                }
            }
        }
        else
        {
            if (!newLink.IsDirichlet())
            {
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    sRet -= __cuCargf(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
                else
                {
                    sRet += __cuCargf(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
            }
        }

        if (pDir[i] > 0) //Move
        {
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
    }

    return sRet;
}


__CLG_REGISTER_HELPER_HEADER(CMeasureBerryPhase)

class CLGAPI CMeasureBerryPhase : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureBerryPhase)
public:
    CMeasureBerryPhase()
        : CMeasure()
        , m_bWilsonDirac(TRUE)
        , m_bGuageFixing(FALSE)
        , m_pMomentumField(NULL)
        , m_pU1Field(NULL)
    {
        
    }

    ~CMeasureBerryPhase();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAccepted(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const class CFieldGauge* const* pCorrespondingStaple = NULL) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

    void CalculateMomentumSpacePhiWilsonDiracForPoint(const SSmallInt4& xprime, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields);
    void CalculateMomentumSpacePhiWilsonDirac(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields);
    void CalculateU1FieldWilsonDirac();
    void CalculateMomentumSpacePhiKSForPoint(const SSmallInt4& xprime, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields);
    void CalculateMomentumSpacePhiKS(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields);
    void CalculateU1FieldKS();

    void CalculateBerryPhase();
    void AllocateBuffers();

    TArray<TArray<DOUBLE>> m_lstData;
    TArray<TArray<DOUBLE>> m_lstDataXY;
    TArray<TArray<DOUBLE>> m_lstDataXZ;
    TArray<TArray<DOUBLE>> m_lstDataXT;
    TArray<TArray<DOUBLE>> m_lstDataYZ;
    TArray<TArray<DOUBLE>> m_lstDataYT;
    TArray<TArray<DOUBLE>> m_lstDataZT;

    UBOOL m_bWilsonDirac;
    UBOOL m_bGuageFixing;

    CFieldFermion* m_pMomentumField;
    CFieldGaugeU1* m_pU1Field;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREBERRYPHASE_H_

//=============================================================================
// END OF FILE
//=============================================================================