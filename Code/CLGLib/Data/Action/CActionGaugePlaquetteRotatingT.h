//=============================================================================
// FILENAME : CActionGaugePlaquetteRotatingT.h
// 
// DESCRIPTION:
// This is the class for rotating guage action
// Open boundary condition (identity Dirichlet boundary condition) is assumed 
// 
//
// REVISION:
//  [mm/dd/yy]
//  [07/08/2024 nbale]
//=============================================================================
#pragma once

#include "Data/Field/Gauge/CFieldGaugeLink.h"

#ifndef _CACTIONGAUGEPLAQUETTE_ROTATINGT_H_
#define _CACTIONGAUGEPLAQUETTE_ROTATINGT_H_

#define __DEFINE_ROTATIONACTION(n) \
__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteRotatingSU##n) \
class CLGAPI CActionGaugePlaquetteRotatingSU##n : public CActionGaugePlaquetteRotatingT<deviceSU##n, n> \
{ \
    __CLGDECLARE_CLASS(CActionGaugePlaquetteRotatingSU##n) \
};

__BEGIN_NAMESPACE

// The rotating (S1/S2) terms are derived assuming the Wilson normalization
// beta = 6/g^2. The improved gauge fields instead take beta = 10/g^2 (the
// plaquette coefficient is normalized to 1, see CFieldGaugeSU3TreeImproved.h).
// The scale factor below converts the incoming beta back to the 6/g^2 that the
// rotating kernels expect, and is applied ONLY to the rotating terms; the S0
// plaquette/rectangle term is delegated to the gauge field and keeps the raw beta.
//   ERBSF_Naive          -> 1.0   (plain Wilson beta = 6/g^2)
//   ERBSF_TreeImprove    -> 0.6   (6/10, tree level improved beta = 10/g^2)
//   ERBSF_OneLoopImprove -> 0.6   (6/10, one-loop improved also uses beta = 10/g^2)
//   ERBSF_Custom         -> read "ScaleFactor" (default 1.0)
__DEFINE_ENUM(ERotationBetaScaleFactor,
    ERBSF_Naive,
    ERBSF_TreeImprove,
    ERBSF_OneLoopImprove,
    ERBSF_Custom,

    ERBSF_ForceDWORD = 0x7fffffff,
    )

template<typename deviceGauge, INT matrixN>
class __DLL_EXPORT CActionGaugePlaquetteRotatingT : public CAction
{
public:

    CActionGaugePlaquetteRotatingT()
        : CAction()
        , m_fOmega(0.0)
        , m_fXi(1.0)
        , m_bCloverEnergy(FALSE)
        , m_bShiftHalfCoord(FALSE)
        , m_bTorus(FALSE)
        , m_fS0Energy(0.0)
        , m_fS1Energy(0.0)
        , m_fS2Energy(0.0)
        , m_eBetaScaleType(ERBSF_Naive)
        , m_fRotationScaleFactor(1.0)
    {
        //SetGaugeOmega(0.0);
    }
    ~CActionGaugePlaquetteRotatingT() {}

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;
    CCString GetInfos(const CCString& tab) const override;

    void SetGaugeOmega(DOUBLE fOmega)
    {
        m_fOmega = fOmega;
    }
    DOUBLE GetOmega() const { return m_fOmega; }
    DOUBLE GetAnisotropy() const { return m_fXi; }
    UBOOL IsCloverEnergy() const { return m_bCloverEnergy; }

    UINT GetDefaultMatrixN() const override { return matrixN; }

    DOUBLE GetS0Energy() const { return m_fS0Energy; }
    DOUBLE GetS1Energy() const { return m_fS1Energy; }
    DOUBLE GetS2Energy() const { return m_fS2Energy; }

    // beta/N to feed the rotating (S1/S2) terms, rescaled to the Wilson 6/g^2 convention
    DOUBLE GetRotationBetaOverN() const { return m_fBetaOverN * m_fRotationScaleFactor; }

protected:

    DOUBLE EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStaple = NULL) override;
    UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;

    void EnergyDirichlet(const class CFieldGaugeLink<deviceGauge, matrixN>* pGauge);
    void EnergyProjectivePlane(const class CFieldGaugeLink<deviceGauge, matrixN>* pGauge);
    void EnergyTorus(const class CFieldGaugeLink<deviceGauge, matrixN>* pGauge);

    void CalculateForceOnGaugeDirichlet(const class CFieldGaugeLink<deviceGauge, matrixN>* pGauge, class CFieldGaugeLink<deviceGauge, matrixN>* pForce) const;
    void CalculateForceOnGaugeProjectivePlane(const class CFieldGaugeLink<deviceGauge, matrixN>* pGauge, class CFieldGaugeLink<deviceGauge, matrixN>* pForce) const;
    void CalculateForceOnGaugeTorus(const class CFieldGaugeLink<deviceGauge, matrixN>* pGauge, class CFieldGaugeLink<deviceGauge, matrixN>* pForce) const;


    DOUBLE m_fOmega;
    DOUBLE m_fXi;
    UBOOL m_bCloverEnergy;
    UBOOL m_bShiftHalfCoord;
    UBOOL m_bTorus;

    DOUBLE m_fS0Energy;
    DOUBLE m_fS1Energy;
    DOUBLE m_fS2Energy;

    ERotationBetaScaleFactor m_eBetaScaleType;
    DOUBLE m_fRotationScaleFactor;

};

#if !_CLG_WIN
extern template class CActionGaugePlaquetteRotatingT<CLGComplex, 1>;
extern template class CActionGaugePlaquetteRotatingT<deviceSU2, 2>;
extern template class CActionGaugePlaquetteRotatingT<deviceSU3, 3>;
extern template class CActionGaugePlaquetteRotatingT<deviceSU4, 4>;
#endif

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteRotatingU1) 
class CLGAPI CActionGaugePlaquetteRotatingU1 : public CActionGaugePlaquetteRotatingT<CLGComplex, 1> 
{ 
    __CLGDECLARE_CLASS(CActionGaugePlaquetteRotatingU1)
};

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteRotating)
class CLGAPI CActionGaugePlaquetteRotating : public CActionGaugePlaquetteRotatingT<deviceSU3, 3>
{
    __CLGDECLARE_CLASS(CActionGaugePlaquetteRotating)
};

__DEFINE_ROTATIONACTION(2)
__DEFINE_ROTATIONACTION(4)

//cost too much time to build in release time, define them only when neccessary
//__DEFINE_ROTATIONACTION(5)
//__DEFINE_ROTATIONACTION(6)
//__DEFINE_ROTATIONACTION(7)
//__DEFINE_ROTATIONACTION(8)

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_ROTATINGT_H_

//=============================================================================
// END OF FILE
//=============================================================================