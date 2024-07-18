//=============================================================================
// FILENAME : CActionGaugePlaquetteRotatingT.h
// 
// DESCRIPTION:
// This is the class for rotating guage action
// Open boundary condition (identity Dirichlet boundary condition) is assumed 
// 
//
// REVISION:
//  [07/08/2024 nbale]
//=============================================================================
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

template<typename deviceGauge, INT matrixN>
class __DLL_EXPORT CActionGaugePlaquetteRotatingT : public CAction
{
public:

    CActionGaugePlaquetteRotatingT();
    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;
    CCString GetInfos(const CCString& tab) const override;

    void SetBeta(DOUBLE fBeta);
    void SetOmega(DOUBLE fOmega);
    DOUBLE GetOmega() const { return m_fOmega; }
    UBOOL IsCloverEnergy() const { return m_bCloverEnergy; }

protected:

    DOUBLE EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable = NULL) override;
    UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;

    void EnergyDirichlet(const class CFieldGaugeLink<deviceGauge, matrixN>* pGauge);
    void EnergyProjectivePlane(const class CFieldGaugeLink<deviceGauge, matrixN>* pGauge);
    void EnergyTorus(const class CFieldGaugeLink<deviceGauge, matrixN>* pGauge);

    void CalculateForceOnGaugeDirichlet(const class CFieldGaugeLink<deviceGauge, matrixN>* pGauge, class CFieldGaugeLink<deviceGauge, matrixN>* pForce) const;
    void CalculateForceOnGaugeProjectivePlane(const class CFieldGaugeLink<deviceGauge, matrixN>* pGauge, class CFieldGaugeLink<deviceGauge, matrixN>* pForce) const;
    void CalculateForceOnGaugeTorus(const class CFieldGaugeLink<deviceGauge, matrixN>* pGauge, class CFieldGaugeLink<deviceGauge, matrixN>* pForce) const;


    DOUBLE m_fOmega;
    UBOOL m_bCloverEnergy;
    UBOOL m_bShiftHalfCoord;
    UBOOL m_bTorus;

};

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteRotatingU1) 
class CLGAPI CActionGaugePlaquetteRotatingU1 : public CActionGaugePlaquetteRotatingT<CLGComplex, 1> 
{ 
    __CLGDECLARE_CLASS(CActionGaugePlaquetteRotatingU1)
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