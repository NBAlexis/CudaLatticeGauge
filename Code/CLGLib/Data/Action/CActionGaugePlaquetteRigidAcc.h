//=============================================================================
// FILENAME : CActionGaugePlaquetteRigidAcc.h
// 
// DESCRIPTION:
// 
// Dirichlet boundary is assumed
// We assume Z direction is Dirichlet, and Nz=0 is not included, for simplicity
// 
// it corresponds to:
// 
// (1+gz)^2   0      0      0
//    0      -1      0      0
//    0       0     -1      0
//    0       0      0     -1
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

    void SetBeta(Real fBeta);
    static void SetG(Real fOmega);
    //void SetCenter(const SSmallInt4 &newCenter);
    //Real GetEnergyPerPlaqutte() const;

protected:

#if !_CLG_DOUBLEFLOAT
    DOUBLE m_fLastEnergy;
    DOUBLE m_fNewEnergy;
#else
    Real m_fLastEnergy;
    Real m_fNewEnergy;
#endif
    Real m_fBetaOverN;
    UINT m_uiPlaqutteCount;

    UBOOL m_bDirichlet;
};

#pragma region Device functions

//================= Put those device functions to header file because we will use them ==============

/**
 * We assume Z direction is Dirichlet, and Nz=0 is not included
 * g(n) = f(n) + f(n+mu) + f(n+nu) + f(n+mu+nu)
 */
static __device__ __inline__ Real _deviceGnRigidAccLeft(const SSmallInt4& sSite, Real fG, BYTE mu, BYTE nu, UBOOL bDirichlet)
{
    if (mu == 2 || nu == 2)
    {
        if (sSite.z == static_cast<SBYTE>(_DC_Lz) - 1)
        {
            if (bDirichlet)
            {
                return (fG * sSite.z + F(1.0)) * F(0.5);
            }
            else
            {
                //(gz+1) + (g(z+1)+1) = (gz+1)+ (g*0 + 1)
                return (fG * sSite.z + F(2.0)) * F(0.5);
            }
        }
        if ((0 == sSite.z) && bDirichlet)
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
static __device__ __inline__ Real _deviceGnRigidAccRight(const SSmallInt4& sSite, Real fG, BYTE mu, BYTE nu, UBOOL bDirichlet)
{
    if (mu == 2)
    {
        if (sSite.z == static_cast<SBYTE>(_DC_Lz) - 1)
        {
            if (bDirichlet)
            {
                return (fG * sSite.z + F(1.0)) * F(0.5);
            }
            else
            {
                //(gz+1) + (g(z+1)+1) = (gz+1)+ (g*0 + 1)
                return (fG * sSite.z + F(2.0)) * F(0.5);
            }
        }
        if ((0 == sSite.z) && bDirichlet)
        {
            return (fG + F(1.0)) * F(0.5);
        }
        return fG * F(0.5) * (sSite.z + sSite.z + 1) + F(1.0);
    }
    if (nu == 2) 
    {
        if (0 == sSite.z)
        {
            if (bDirichlet)
            {
                printf("should never be here!\n");
                return (fG + F(1.0)) * F(0.5);
            }
            // g*0 + 1 + g*(Lz-1) + 1
            return (fG * (_DC_Lz - 1) + F(2.0)) * F(0.5);
        }
        if (1 == sSite.z && bDirichlet)
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
 * Note: f(n) = (1 + g z)^3 so it is [(1 + g z)^3 + (1 + g (z+1))^3]/2
 * NOT (1 + g z + 1 + g(z+1) )^3
 */
static __device__ __inline__ Real _deviceGnRigidAccLeftTri(const SSmallInt4& sSite, Real fG, BYTE mu, BYTE nu, UBOOL bDirichlet)
{
    const Real ret = fG * sSite.z + F(1.0);
    if (mu == 2 || nu == 2)
    {
        Real ret2 = fG * (sSite.z + 1) + F(1.0);
        if (0 == sSite.z && bDirichlet)
        {
            return ret2 * ret2 * ret2 * F(0.5);
        }

        if (sSite.z == static_cast<SBYTE>(_DC_Lz) - 1)
        {
            if (bDirichlet)
            {
                return ret * ret * ret * F(0.5);
            }
            else
            {
                // sSite.z + 1 =0
                ret2 = F(1.0);
            }
        }

        return F(0.5) * (ret * ret * ret + ret2 * ret2 * ret2);
    }
    return ret * ret * ret;
}

/**
 * g(n-nu) = f(n-nu) + f(n+mu-nu) + f(n) + f(n+mu)
 */
static __device__ __inline__ Real _deviceGnRigidAccRightTri(const SSmallInt4& sSite, Real fG, BYTE mu, BYTE nu, UBOOL bDirichlet)
{
    Real ret = fG * sSite.z + F(1.0);
    if (2 == mu)
    {
        Real ret2 = fG * (sSite.z + 1) + F(1.0);
        if (0 == sSite.z && bDirichlet)
        {
            return ret2 * ret2 * ret2 * F(0.5);
        }

        if (sSite.z == static_cast<SBYTE>(_DC_Lz) - 1)
        {
            if (bDirichlet)
            {
                return ret * ret * ret * F(0.5);
            }
            else
            {
                // sSite.z + 1 =0
                ret2 = F(1.0);
            }
        }

        return F(0.5) * (ret * ret * ret + ret2 * ret2 * ret2);
    }

    if (2 == nu)
    {
        if (0 == sSite.z)
        {
            if (bDirichlet)
            {
                printf("should never be here!\n");
                return ret * ret * ret * F(0.5);
            }
            const Real ret2 = fG * (_DC_Lz - 1) + F(1.0);
            return F(0.5) * (ret * ret * ret + ret2 * ret2 * ret2);
        }
        if (1 == sSite.z && bDirichlet)
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