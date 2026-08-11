//=============================================================================
// FILENAME : CMeasureBosonValue.h
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [mm/dd/yy]
//  [11/07/2024 nbale]
//=============================================================================

#ifndef _CMEASUREBOSONVALUE_H_
#define _CMEASUREBOSONVALUE_H_

#define __DEFINE_BOSON_VALUE_MEASURE(MEASURE_NAME, TYPE_BOSON, TYPE_GAUGE) \
__CLG_REGISTER_HELPER_HEADER(MEASURE_NAME) \
class CLGAPI MEASURE_NAME : public TMeasureBosonValue<TYPE_BOSON, TYPE_GAUGE> \
{ \
    __CLGDECLARE_CLASS(MEASURE_NAME) \
};

__BEGIN_NAMESPACE

template<typename deviceDataBoson, typename deviceDataGauge>
class __DLL_EXPORT TMeasureBosonValue : public CMeasure
{

public:

    void OnConfigurationAccepted(INT gaugeNum, INT bosonNum, INT tensor2Num, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldTensor2* const* tensor2Fields, const class CFieldGauge* const* pCorrespondingStaple) override;


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
    void OnConfigurationAcceptedZ4(INT gaugeNum, INT bosonNum, INT tensor2Num, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldTensor2* const* tensor2Fields, const class CFieldGauge* const* pCorrespondingStaple, const class CFieldFermion* pZ4, const class CFieldFermion* pInverseZ4, UBOOL bStart, UBOOL bEnd) override
    {

    }

    void Reset() override;
    void Report() override {}

    UBOOL IsGaugeOrBosonMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

    TArray<cuDoubleComplex> m_lstEveryConfigurationC;
    TArray<DOUBLE> m_lstEveryConfigurationR;

    TArray<TArray<cuDoubleComplex>> m_lstEveryConfigurationZsliceC;
    TArray<TArray<DOUBLE>> m_lstEveryConfigurationZsliceR;

};

__DEFINE_BOSON_VALUE_MEASURE(CMeasureBosonValueReal, Real, Real)
__DEFINE_BOSON_VALUE_MEASURE(CMeasureBosonValueU1, CLGComplex, CLGComplex)
__DEFINE_BOSON_VALUE_MEASURE(CMeasureBosonValueSU2, deviceSU2Vector, deviceSU2)
__DEFINE_BOSON_VALUE_MEASURE(CMeasureBosonValueSU3, deviceSU3Vector, deviceSU3)

#if _CLG_SU4_BOSON
__DEFINE_BOSON_VALUE_MEASURE(CMeasureBosonValueSU4, deviceSU4Vector, deviceSU4)
#endif
#if _CLG_SU5_BOSON
__DEFINE_BOSON_VALUE_MEASURE(CMeasureBosonValueSU5, deviceSU5Vector, deviceSU5)
#endif
#if _CLG_SU6_BOSON
__DEFINE_BOSON_VALUE_MEASURE(CMeasureBosonValueSU6, deviceSU6Vector, deviceSU6)
#endif
#if _CLG_SU7_BOSON
__DEFINE_BOSON_VALUE_MEASURE(CMeasureBosonValueSU7, deviceSU7Vector, deviceSU7)
#endif
#if _CLG_SU8_BOSON
__DEFINE_BOSON_VALUE_MEASURE(CMeasureBosonValueSU8, deviceSU8Vector, deviceSU8)
#endif

__END_NAMESPACE

#endif //#ifndef _CMEASUREBOSONVALUE_H_

//=============================================================================
// END OF FILE
//=============================================================================