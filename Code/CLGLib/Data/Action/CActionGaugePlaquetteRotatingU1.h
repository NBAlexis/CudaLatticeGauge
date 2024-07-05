//=============================================================================
// FILENAME : CActionGaugePlaquetteRotatingU1.h
// 
// DESCRIPTION:
// This is the class for rotating guage action
// Open boundary condition (identity Dirichlet boundary condition) is assumed 
// 
//
// REVISION:
//  [10/01/2021 nbale]
//=============================================================================

#ifndef _CACTIONGAUGEPLAQUETTE_ROTATINGU1_H_
#define _CACTIONGAUGEPLAQUETTE_ROTATINGU1_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteRotatingU1)

class CLGAPI CActionGaugePlaquetteRotatingU1 : public CAction
{
    __CLGDECLARE_CLASS(CActionGaugePlaquetteRotatingU1)
public:

    CActionGaugePlaquetteRotatingU1();

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;
    CCString GetInfos(const CCString &tab) const override;

    void SetBeta(DOUBLE fBeta);
    void SetOmega(DOUBLE fOmega);
    //void SetCenter(const SSmallInt4 &newCenter);
    //Real GetEnergyPerPlaqutte() const;

#if !_CLG_DOUBLEFLOAT
    DOUBLE m_fOmega;
#else
    Real m_fOmega;
#endif
    UBOOL m_bCloverEnergy;
    UBOOL m_bShiftHalfCoord;

    //===== test functions ======
    //DOUBLE XYTerm1(const class CFieldGauge* pGauge);
    //DOUBLE XYTerm2(const class CFieldGauge* pGauge);

protected:

    DOUBLE EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable = NULL) override;
    UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge* pGauge, class CFieldGauge* pForce, class CFieldGauge* pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;

    UINT m_uiPlaqutteCount;
};


__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_ROTATINGU1_H_

//=============================================================================
// END OF FILE
//=============================================================================