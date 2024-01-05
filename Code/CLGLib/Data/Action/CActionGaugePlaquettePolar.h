//=============================================================================
// FILENAME : CActionGaugePlaquetteRigidAcc.h
// 
// DESCRIPTION:
// 
// 
// it corresponds to:
// 
// p = phi
// S = r(frt^2 + frz^2 + fzt^2) + (1/r)(fpr^2+fpt^2+fpz^2)
// 
// It can be viewed as polar coordinate
// where:
// x = r
// y = a phi
// z = z
// t = t
//
// REVISION:
//  [05/01/2024 nbale]
//=============================================================================

#ifndef _CACTIONGAUGEPLAQUETTE_POLAR_H_
#define _CACTIONGAUGEPLAQUETTE_POLAR_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquettePolar)

class CLGAPI CActionGaugePlaquettePolar : public CAction
{
    __CLGDECLARE_CLASS(CActionGaugePlaquettePolar)
public:

    CActionGaugePlaquettePolar();
    ~CActionGaugePlaquettePolar();

#if !_CLG_DOUBLEFLOAT
    DOUBLE Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable = NULL) override;
#else
    Real Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable = NULL) override;
#endif
    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;

    UBOOL CalculateForceOnGauge(const class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMC(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;
    void OnFinishTrajectory(UBOOL bAccepted) override;
    CCString GetInfos(const CCString &tab) const override;

    void SetBetaList(const TArray<Real>& betalst);
    void SetR(Real fIn, Real fOut);

protected:

#if !_CLG_DOUBLEFLOAT
    DOUBLE m_fLastEnergy;
    DOUBLE m_fNewEnergy;
#else
    Real m_fLastEnergy;
    Real m_fNewEnergy;
#endif
    UINT m_uiPlaqutteCount;

    UBOOL m_bDirichlet;

    Real m_fRIn;
    Real m_fROut;
    Real m_fDeltaR;
    TArray<Real> m_lstBeta;
    Real* m_deviceBetaList;
};

#pragma region Device functions

//================= Put those device functions to header file because we will use them ==============

/**
 * We assume Z direction is Dirichlet, and Nz=0 is not included
 * g(n) = f(n) + f(n+mu) + f(n+nu) + f(n+mu+nu)
 * f(n) = beta * r
 */
static __device__ __inline__ Real _deviceGnPolarSpatialLeft(
    const SSmallInt4& sSite,
    const Real* __restrict__ betaList,
    Real fRIn,
    Real fRDeltaR,
    BYTE mu, BYTE nu, UBOOL bDirichlet)
{
    const Real betaR = betaList[sSite.x] * (fRIn + fRDeltaR * sSite.x);

    if (mu == 0 || nu == 0)
    {
        const SBYTE rp1 = (static_cast<SBYTE>(_DC_Lx - 1) == sSite.x) ? 0 : (sSite.x + 1);
        const Real betaR2 = betaList[rp1] * (fRIn + fRDeltaR * rp1);
        if (bDirichlet)
        {
            if (0 == sSite.x)
            {
                return betaR2 * F(0.5);
            }
            else if (0 == rp1)
            {
                return betaR * F(0.5);
            }
        }
        return F(0.5) * (betaR + betaR2);
    }
    return betaR;
}

/**
 * g(n-nu) = f(n-nu) + f(n+mu-nu) + f(n) + f(n+mu)
 */
static __device__ __inline__ Real _deviceGnPolarSpatialRight(
    const SSmallInt4& sSite,
    const Real* __restrict__ betaList,
    Real fRIn,
    Real fRDeltaR,
    BYTE mu, BYTE nu, UBOOL bDirichlet)
{
    const Real betaR = betaList[sSite.x] * (fRIn + fRDeltaR * sSite.x);
    if (0 == mu)
    {
        const SBYTE rp1 = (static_cast<SBYTE>(_DC_Lx - 1) == sSite.x) ? 0 : (sSite.x + 1);
        const Real betaR2 = betaList[rp1] * (fRIn + fRDeltaR * rp1);

        if (bDirichlet)
        {
            if (0 == sSite.x)
            {
                return betaR2 * F(0.5);
            }
            else if (0 == rp1)
            {
                return betaR * F(0.5);
            }
        }

        return F(0.5) * (betaR + betaR2);
    }

    if (0 == nu)
    {
        const SBYTE rm1 = (0 == sSite.x) ? (static_cast<SBYTE>(_DC_Lx) - 1) : (sSite.x - 1);
        const Real betaR2 = betaList[rm1] * (fRIn + fRDeltaR * rm1);

        if (bDirichlet)
        {
            if (0 == sSite.x)
            {
                printf("should never be here!\n");
                return betaR * F(0.5);
            }
            else if (1 == sSite.x)
            {
                return betaR * F(0.5);
            }
        }
        return F(0.5) * (betaR + betaR2);
    }

    return betaR;
}


/**
 * We assume Z direction is Dirichlet, and Nz=0 is not included
 * g(n) = f(n) + f(n+mu) + f(n+nu) + f(n+mu+nu)
 * f(n) = b(r) 1/r
 */
static __device__ __inline__ Real _deviceGnPolarPhiLeft(const SSmallInt4& sSite,
    const Real* __restrict__ betaList,
    Real fRIn,
    Real fRDeltaR,
    BYTE mu, BYTE nu, UBOOL bDirichlet)
{
    const Real betaOverR = betaList[sSite.x] / (fRIn + fRDeltaR * sSite.x);

    if (mu == 0 || nu == 0)
    {
        const SBYTE rp1 = (static_cast<SBYTE>(_DC_Lx - 1) == sSite.x) ? 0 : (sSite.x + 1);

        const Real betaOverR2 = betaList[rp1] / (fRIn + fRDeltaR * rp1);

        if (bDirichlet)
        {
            if (0 == sSite.x)
            {
                return betaOverR2 * F(0.5);
            }
            else if (0 == rp1)
            {
                return betaOverR * F(0.5);
            }
        }
        return F(0.5) * (betaOverR + betaOverR2);
    }
    return betaOverR;
}

/**
 * g(n-nu) = f(n-nu) + f(n+mu-nu) + f(n) + f(n+mu)
 */
static __device__ __inline__ Real _deviceGnPolarPhiRight(
    const SSmallInt4& sSite,
    const Real* __restrict__ betaList,
    Real fRIn,
    Real fRDeltaR,
    BYTE mu, BYTE nu, UBOOL bDirichlet)
{
    const Real betaOverR = betaList[sSite.x] / (fRIn + fRDeltaR * sSite.x);
    if (0 == mu)
    {
        const SBYTE rp1 = (static_cast<SBYTE>(_DC_Lx - 1) == sSite.x) ? 0 : (sSite.x + 1);
        const Real betaOverR2 = betaList[rp1] / (fRIn + fRDeltaR * rp1);

        if (bDirichlet)
        {
            if (0 == sSite.x)
            {
                return betaOverR2 * F(0.5);
            }
            else if (0 == rp1)
            {
                return betaOverR * F(0.5);
            }
        }

        return F(0.5) * (betaOverR + betaOverR2);
    }

    if (0 == nu)
    {
        const SBYTE rm1 = (0 == sSite.x) ? (static_cast<SBYTE>(_DC_Lx) - 1) : (sSite.x - 1);
        const Real betaOverR2 = betaList[rm1] / (fRIn + fRDeltaR * rm1);

        if (bDirichlet)
        {
            if (0 == sSite.x)
            {
                printf("should never be here!\n");
                return betaOverR * F(0.5);
            }
            else if (1 == sSite.x)
            {
                return betaOverR * F(0.5);
            }
        }
        return F(0.5) * (betaOverR + betaOverR2);
    }

    return betaOverR;
}


#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_POLAR_H_

//=============================================================================
// END OF FILE
//=============================================================================