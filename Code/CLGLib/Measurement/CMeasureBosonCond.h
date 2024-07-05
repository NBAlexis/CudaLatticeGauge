//=============================================================================
// FILENAME : CMeasurePolyakov.h
// 
// DESCRIPTION:
// This is measurement for Polyakov loop and static quark potential
// Suitable for translational invarient case
//
// REVISION:
//  [07/07/2019 nbale]
//=============================================================================

#ifndef _CMEASUREBOSONCOND_H_
#define _CMEASUREBOSONCOND_H_

__BEGIN_NAMESPACE


__CLG_REGISTER_HELPER_HEADER(CMeasureBosonCond)

class CLGAPI CMeasureBosonCond : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureBosonCond)

public:

    void OnConfigurationAccepted(INT gaugeNum, INT bosonNum, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldGauge* const* pCorrespondingStaple) override;


    /**
    * NOTE: sources will be passed to multiple measures, do NOT change the content!
    * NOTE: site.x start from 1 to Lx - 1, 0 is not included
    */
    void SourceSanning(INT gaugeNum, INT bosonNum, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldGauge* const* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override
    {

    }

    /**
    * Z4 Source
    */
    void OnConfigurationAcceptedZ4(INT gaugeNum, INT bosonNum, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldGauge* const* pCorrespondingStaple, const class CFieldFermion* pZ4, const class CFieldFermion* pInverseZ4, UBOOL bStart, UBOOL bEnd) override
    {

    }

    void Report() override;
    void Reset() override;
    void Average() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

    TArray<TArray<DOUBLE>> m_lstElement;
    TArray<DOUBLE> m_lstAverageElement;

};

__END_NAMESPACE

#endif //#ifndef _CMEASUREBOSONCOND_H_

//=============================================================================
// END OF FILE
//=============================================================================