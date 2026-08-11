//=============================================================================
// FILENAME : CActionGaugePlaquetteRigidAcc.h
// 
// DESCRIPTION:
// 
// 
// it corresponds to:
// 
// (1+gz)^2   0      0      0
//    0      -1      0      0
//    0       0     -1      0
//    0       0      0     -1
// 
// S = (1+gz)(f12^2 + f13^2 + f23^2) + (1/(1+gz))(f01^2+f02^2+f03^2)
//
// REVISION:
//  [mm/dd/yy]
//  [07/31/2020 nbale]
//=============================================================================
#pragma once

#ifndef _CACTIONGAUGEPLAQUETTE_RIGIDACC_H_
#define _CACTIONGAUGEPLAQUETTE_RIGIDACC_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteRigidAcc)

class CLGAPI CActionGaugePlaquetteRigidAcc : public CAction
{
    __CLGDECLARE_CLASS(CActionGaugePlaquetteRigidAcc)
public:

    CActionGaugePlaquetteRigidAcc();

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;
    CCString GetInfos(const CCString &tab) const override;

    void SetBeta(Real fBeta);
    static void SetG(Real fOmega);
    //void SetCenter(const SSmallInt4 &newCenter);
    //Real GetEnergyPerPlaqutte() const;

protected:

    DOUBLE EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStaple = NULL) override;
    UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge* pGauge, class CFieldGauge* pForce, class CFieldGauge* pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;

    UINT m_uiPlaqutteCount;

    UBOOL m_bDirichlet;
};

#pragma region Device functions

//================= Put those device functions to header file because we will use them ==============

/**
 * We assume Z direction is Dirichlet, and Nz=0 is not included
 * g(n) = f(n) + f(n+mu) + f(n+nu) + f(n+mu+nu)
 * f(n) = 1 + g z
 */
static __device__ __inline__ Real _deviceGnRigidAccSpatialLeft(const SSmallInt4& sSite, Real fG, BYTE mu, BYTE nu, UBOOL bDirichlet)
{
    //Multi-GPU: f(n) = 1 + g z needs the GLOBAL z, and the z+1 wrap plus the
    //Dirichlet surface tests must use the global extent/coordinate
    //(P4-1.4/1.4-R1). Identity on single-GPU (offset 0, GlobalLz == Lz).
    const INT iGz = static_cast<INT>(sSite.z) + static_cast<INT>(_DC_OffsetZ);
    const Real ret = fG * (iGz - _DC_Centerz) + F(1.0);

    if (mu == 2 || nu == 2)
    {
        const INT iGzp1 = (static_cast<INT>(_DC_GlobalLz) - 1 == iGz) ? 0 : (iGz + 1);
        const Real ret2 = fG * (iGzp1 - _DC_Centerz) + F(1.0);
        if (bDirichlet)
        {
            if (0 == iGz)
            {
                return ret2 * F(0.5);
            }
            else if (0 == iGzp1)
            {
                return ret * F(0.5);
            }
        }
        return F(0.5) * (ret + ret2);
    }
    return ret;
}

/**
 * g(n-nu) = f(n-nu) + f(n+mu-nu) + f(n) + f(n+mu)
 */
static __device__ __inline__ Real _deviceGnRigidAccSpatialRight(const SSmallInt4& sSite, Real fG, BYTE mu, BYTE nu, UBOOL bDirichlet)
{
    //Multi-GPU: global z/wrap/surface tests (P4-1.4); identity on single-GPU.
    const INT iGz = static_cast<INT>(sSite.z) + static_cast<INT>(_DC_OffsetZ);
    const Real ret = fG * (iGz - _DC_Centerz) + F(1.0);
    if (2 == mu)
    {
        const INT iGzp1 = (static_cast<INT>(_DC_GlobalLz) - 1 == iGz) ? 0 : (iGz + 1);
        const Real ret2 = fG * (iGzp1 - _DC_Centerz) + F(1.0);

        if (bDirichlet)
        {
            if (0 == iGz)
            {
                return ret2 * F(0.5);
            }
            else if (0 == iGzp1)
            {
                return ret * F(0.5);
            }
        }

        return F(0.5) * (ret + ret2);
    }

    if (2 == nu)
    {
        const INT iGzm1 = (0 == iGz) ? (static_cast<INT>(_DC_GlobalLz) - 1) : (iGz - 1);
        const Real ret2 = fG * (iGzm1 - _DC_Centerz) + F(1.0);

        if (bDirichlet)
        {
            if (0 == iGz)
            {
                printf("should never be here!\n");
                return ret * F(0.5);
            }
            else if (1 == iGz)
            {
                return ret * F(0.5);
            }
        }
        return F(0.5) * (ret + ret2);
    }

    return ret;
}


/**
 * We assume Z direction is Dirichlet, and Nz=0 is not included
 * g(n) = f(n) + f(n+mu) + f(n+nu) + f(n+mu+nu)
 * f(n) = 1/(1 + g z)
 */
static __device__ __inline__ Real _deviceGnRigidAccTimeLeft(const SSmallInt4& sSite, Real fG, BYTE mu, BYTE nu, UBOOL bDirichlet)
{
    //Multi-GPU: global z/wrap/surface tests (P4-1.4); identity on single-GPU.
    const INT iGz = static_cast<INT>(sSite.z) + static_cast<INT>(_DC_OffsetZ);
    const Real ret = F(1.0) / (fG * (iGz - _DC_Centerz) + F(1.0));

    if (mu == 2 || nu == 2)
    {
        const INT iGzp1 = (static_cast<INT>(_DC_GlobalLz) - 1 == iGz) ? 0 : (iGz + 1);
        const Real ret2 = F(1.0) / (fG * (iGzp1 - _DC_Centerz) + F(1.0));
        if (bDirichlet)
        {
            if (0 == iGz)
            {
                return ret2 * F(0.5);
            }
            else if (0 == iGzp1)
            {
                return ret * F(0.5);
            }
        }
        return F(0.5) * (ret + ret2);
    }
    return ret;
}

/**
 * g(n-nu) = f(n-nu) + f(n+mu-nu) + f(n) + f(n+mu)
 */
static __device__ __inline__ Real _deviceGnRigidAccTimeRight(const SSmallInt4& sSite, Real fG, BYTE mu, BYTE nu, UBOOL bDirichlet)
{
    //Multi-GPU: global z/wrap/surface tests (P4-1.4); identity on single-GPU.
    const INT iGz = static_cast<INT>(sSite.z) + static_cast<INT>(_DC_OffsetZ);
    const Real ret = F(1.0) / (fG * (iGz - _DC_Centerz) + F(1.0));
    if (2 == mu)
    {
        const INT iGzp1 = (static_cast<INT>(_DC_GlobalLz) - 1 == iGz) ? 0 : (iGz + 1);
        const Real ret2 = F(1.0) / (fG * (iGzp1 - _DC_Centerz) + F(1.0));

        if (bDirichlet)
        {
            if (0 == iGz)
            {
                return ret2 * F(0.5);
            }
            else if (0 == iGzp1)
            {
                return ret * F(0.5);
            }
        }

        return F(0.5) * (ret + ret2);
    }

    if (2 == nu)
    {
        const INT iGzm1 = (0 == iGz) ? (static_cast<INT>(_DC_GlobalLz) - 1) : (iGz - 1);
        const Real ret2 = F(1.0) / (fG * (iGzm1 - _DC_Centerz) + F(1.0));

        if (bDirichlet)
        {
            if (0 == iGz)
            {
                printf("should never be here!\n");
                return ret * F(0.5);
            }
            else if (1 == iGz)
            {
                return ret * F(0.5);
            }
        }
        return F(0.5) * (ret + ret2);
    }

    return ret;
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_RIGIDACC_H_

//=============================================================================
// END OF FILE
//=============================================================================