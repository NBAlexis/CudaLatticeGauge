//=============================================================================
// FILENAME : CActionGaugePlaquetteRigidAcc.h
// 
// DESCRIPTION:
// 
// Dirichlet boundary is assumed
// We assume Z direction is Dirichlet, and Nz=0 is not included, for simplicity
//
// REVISION:
//  [07/31/2020 nbale]
//=============================================================================

#ifndef _CACTIONGAUGEPLAQUETTE_RIGIDACC_H_
#define _CACTIONGAUGEPLAQUETTE_RIGIDACC_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteRigidAcc)

class CLGAPI CActionGaugePlaquetteRigidAcc : public CAction
{
    __CLGDECLARE_CLASS(CActionGaugePlaquetteRigidAcc)
public:

    CActionGaugePlaquetteRigidAcc();

    Real Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable) override;
    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;

    UBOOL CalculateForceOnGauge(const class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMC(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;
    void OnFinishTrajectory(UBOOL bAccepted) override;
    CCString GetInfos(const CCString &tab) const override;

    void SetBeta(Real fBeta);
    static void SetG(Real fOmega);
    //void SetCenter(const SSmallInt4 &newCenter);
    //Real GetEnergyPerPlaqutte() const;

protected:

    Real m_fLastEnergy;
    Real m_fNewEnergy;
    Real m_fBetaOverN;
    UINT m_uiPlaqutteCount;
};

#pragma region Device functions

//================= Put those device functions to header file because we will use them ==============

/**
 * We assume Z direction is Dirichlet, and Nz=0 is not included
 * g(n) = f(n) + f(n+mu) + f(n+nu) + f(n+mu+nu)
 */
static __device__ __inline__ Real _deviceGnRigidAccLeft(const SSmallInt4& sSite, Real fG, BYTE mu, BYTE nu)
{
    if (mu == 2 || nu == 2)
    {
        if (sSite.z == static_cast<SBYTE>(_DC_Lz) - 1)
        {
            return (fG * sSite.z + F(1.0)) * F(0.5);
        }
        if (0 == sSite.z)
        {
            return (fG + F(1.0)) * F(0.5);
        }
        return fG * F(0.5) * (sSite.z + sSite.z + 1) + F(1.0);
    }
    return fG * sSite.z + F(1.0);
}

/**
 * g(n-nu) = f(n-nu) + f(n+mu-nu) + f(n) + f(n+mu)
 */
static __device__ __inline__ Real _deviceGnRigidAccRight(const SSmallInt4& sSite, Real fG, BYTE mu, BYTE nu)
{
    if (mu == 2)
    {
        if (sSite.z == static_cast<SBYTE>(_DC_Lz) - 1)
        {
            return (fG * sSite.z + F(1.0)) * F(0.5);
        }
        if (0 == sSite.z)
        {
            return (fG + F(1.0)) * F(0.5);
        }
        return fG * F(0.5) * (sSite.z + sSite.z + 1) + F(1.0);
    }
    if (nu == 2) //never here with z = 0
    {
        if (sSite.z <= 1)
        {
            return (fG + F(1.0)) * F(0.5);
        }
        return fG * F(0.5) * (sSite.z + sSite.z - 1) + F(1.0);
    }

    return fG * sSite.z + F(1.0);
}


/**
 * We assume Z direction is Dirichlet, and Nz=0 is not included
 * g(n) = f(n) + f(n+mu) + f(n+nu) + f(n+mu+nu)
 * Note: f(n) = (1 + g z)^3 so it is (1 + g z)^3 + (1 + g (z+1))^3
 * NOT (1 + g z + 1 + g(z+1) )^3
 */
static __device__ __inline__ Real _deviceGnRigidAccLeftTri(const SSmallInt4& sSite, Real fG, BYTE mu, BYTE nu)
{
    const Real ret = fG * sSite.z + F(1.0);
    if (mu == 2 || nu == 2)
    {
        if (sSite.z == static_cast<SBYTE>(_DC_Lz) - 1)
        {
            return ret * ret * ret * F(0.5);
        }
        const Real ret2 = fG * (sSite.z + 1) + F(1.0);
        if (0 == sSite.z)
        {
            return ret2 * ret2 * ret2 * F(0.5);
        }
        return F(0.5) * (ret * ret * ret + ret2 * ret2 * ret2);
    }
    return ret * ret * ret;
}

/**
 * g(n-nu) = f(n-nu) + f(n+mu-nu) + f(n) + f(n+mu)
 */
static __device__ __inline__ Real _deviceGnRigidAccRightTri(const SSmallInt4& sSite, Real fG, BYTE mu, BYTE nu)
{
    const Real ret = fG * sSite.z + F(1.0);
    if (mu == 2)
    {
        if (sSite.z == static_cast<SBYTE>(_DC_Lz) - 1)
        {
            return ret * ret * ret * F(0.5);
        }
        const Real ret2 = fG * (sSite.z + 1) + F(1.0);
        if (0 == sSite.z)
        {
            return ret2 * ret2 * ret2 * F(0.5);
        }
        return F(0.5) * (ret * ret * ret + ret2 * ret2 * ret2);
    }
    if (nu == 2) //never here with z = 0
    {
        if (sSite.z <= 1)
        {
            return ret * ret * ret * F(0.5);
        }
        const Real ret2 = fG * (sSite.z - 1) + F(1.0);
        return F(0.5) * (ret * ret * ret + ret2 * ret2 * ret2);
    }

    return ret * ret * ret;
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_RIGIDACC_H_

//=============================================================================
// END OF FILE
//=============================================================================