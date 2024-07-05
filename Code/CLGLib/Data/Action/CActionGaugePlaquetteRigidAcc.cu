//=============================================================================
// FILENAME : CActionGaugePlaquetteRigidAcc.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [07/31/2020 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CActionGaugePlaquetteRigidAcc.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteRigidAcc)


#pragma region kernels

/**
* 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelEnergy_RigidAcc(
    const deviceSU3 * __restrict__ pDeviceData,
    BYTE byFieldId,
    Real betaOverN, Real fG,
    SSmallInt4 sCenter,
    DOUBLE* results)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];

    if (sIdx.IsDirichlet())
    {
        results[uiSiteIndex] = F(0.0);
        return;
    }

    const Real OnePlusGZ = fG * (sSite4.z - sCenter.z) + F(1.0);
    const Real OneOverOnePlusGZ = F(1.0) / (OnePlusGZ * OnePlusGZ);
    deviceSU3 toSub(_deviceClover(pDeviceData, sSite4, uiBigIdx, 3, 0, byFieldId));
    toSub.Add(_deviceClover(pDeviceData, sSite4, uiBigIdx, 3, 1, byFieldId));
    toSub.Add(_deviceClover(pDeviceData, sSite4, uiBigIdx, 3, 2, byFieldId));


    toSub.MulReal(OneOverOnePlusGZ); //Now this is (1/(1+gz)^2)(U14 + U24 + U34)
    toSub.Add(_deviceClover(pDeviceData, sSite4, uiBigIdx, 0, 1, byFieldId));
    toSub.Add(_deviceClover(pDeviceData, sSite4, uiBigIdx, 0, 2, byFieldId));
    toSub.Add(_deviceClover(pDeviceData, sSite4, uiBigIdx, 1, 2, byFieldId));
    //Now this is U12 + U13 + U23 + (1/(1+gz)^2)(U14 + U24 + U34)

    const Real toAdd = F(9.0) * OnePlusGZ + F(9.0) / OnePlusGZ - toSub.ReTr() * F(0.25) * OnePlusGZ;
    //Now this is ReTr[(1+gz)(3 - U12 - U13 - U23) + (1/(1+gz)) (3 - U14 - U24 - U34)]

    results[uiSiteIndex] = betaOverN * toAdd;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelEnergy_RigidAcc_Simplified(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
    BYTE plaqLength, BYTE plaqCount,
    Real betaOverN, Real fG,
    UBOOL bDirichlet,
#if !_CLG_DOUBLEFLOAT
    DOUBLE* results
#else
    Real* results
#endif
)
{
    intokernalInt4;

#if !_CLG_DOUBLEFLOAT
    DOUBLE resThisThread = 0.0;
#else
    Real resThisThread = F(0.0);
#endif
    const UINT plaqCountAll = plaqCount * plaqLength;
    for (BYTE i = 0; i < plaqCount; ++i)
    {
        SIndex first = pCachedIndex[i * plaqLength + uiSiteIndex * plaqCountAll];
        const BYTE mu = first.m_byDir;
        deviceSU3 toAdd(_deviceGetGaugeBCSU3(byFieldId, pDeviceData, first));

        if (first.NeedToDagger())
        {
            toAdd.Dagger();
        }

        first = pCachedIndex[i * plaqLength + 1 + uiSiteIndex * plaqCountAll];
        const BYTE nu = first.m_byDir;
        deviceSU3 toMul(_deviceGetGaugeBCSU3(byFieldId, pDeviceData, first));
        if (first.NeedToDagger())
        {
            toAdd.MulDagger(toMul);
        }
        else
        {
            toAdd.Mul(toMul);
        }

        for (BYTE j = 2; j < plaqLength; ++j)
        {
            first = pCachedIndex[i * plaqLength + j + uiSiteIndex * plaqCountAll];
            toMul = _deviceGetGaugeBCSU3(byFieldId, pDeviceData, first);
            if (first.NeedToDagger())
            {
                toAdd.MulDagger(toMul);
            }
            else
            {
                toAdd.Mul(toMul);
            }
        }

        //i=0: 12
        //  1: 13
        //  2: 14
        //  3: 23
        //  4: 24
        //  5: 34
        if (2 == i || 4 == i || 5 == i)
        {
            resThisThread += (F(3.0) - toAdd.ReTr()) * _deviceGnRigidAccTimeLeft(sSite4, fG, mu, nu, bDirichlet);
        }
        else
        {
            resThisThread += (F(3.0) - toAdd.ReTr()) * _deviceGnRigidAccSpatialLeft(sSite4, fG, mu, nu, bDirichlet);
        }
    }

    results[uiSiteIndex] = resThisThread * betaOverN;
}

/**
* Assuming x-y-t direction is NOT Dirichlet
* If \bar{U}_{mu nu} does not have Z direction, g(n) = f(n)
* If \bar{U}_{mu nu} does have Z direction, g(n) = [f(n) + f(n+Z)]/2
* This can be Dirichlet, because we have links like U_z(nz=0)
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForce4PlaqutteTermSU3_RigidAcc(
    const deviceSU3* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
    BYTE plaqLength, BYTE plaqCount,
    deviceSU3* pForceData,
    Real betaOverN, 
    Real fG,
    UBOOL bDirichlet,
    BYTE byFieldId)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    //const SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];

    //Real test_force = F(0.0);
    betaOverN = betaOverN * F(-0.5);
    const UINT plaqLengthm1 = plaqLength - 1;
    const UINT plaqCountAll = plaqCount * plaqLengthm1;

    #pragma unroll
    for (UINT idir = 0; idir < 4; ++idir)
    {
        if (__idx->_deviceIsBondOnSurface(uiBigIdx, idir))
        {
            continue;
        }

        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceSU3 res = deviceSU3::makeSU3Zero();
        const BYTE mu = idir;
        //there are 6 staples,
        //3 'other directions' each is sum of two plaquttes
        for (int i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            const BYTE nu = first.m_byDir;

            //if (mu != 2 && nu != 2 && 0 == sSite4.z)
            //{
            //    //All U on surface, the dynamic links with z=0 are U_z links
            //    continue;
            //}

            deviceSU3 toAdd(_deviceGetGaugeBCSU3(byFieldId, pDeviceData, first));
            Real fFactorG = F(0.0);
            if (first.NeedToDagger())
            {
                toAdd.Dagger();
                //if (0 == sSite4.z)
                //{
                //    continue; //fFactorG = 0
                //}

                if (mu == 3 || nu == 3) //one of mu nu is t
                {
                    fFactorG = _deviceGnRigidAccTimeRight(sSite4, fG, mu, nu, bDirichlet);
                }
                else
                {
                    fFactorG = _deviceGnRigidAccSpatialRight(sSite4, fG, mu, nu, bDirichlet);
                }
            }
            else
            {
                if (mu == 3 || nu == 3) //one of mu nu is t
                {
                    fFactorG = _deviceGnRigidAccTimeLeft(sSite4, fG, mu, nu, bDirichlet);
                }
                else
                {
                    fFactorG = _deviceGnRigidAccSpatialLeft(sSite4, fG, mu, nu, bDirichlet);
                }
            }

            for (BYTE j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                deviceSU3 toMul(_deviceGetGaugeBCSU3(byFieldId, pDeviceData, nextlink));

                if (nextlink.NeedToDagger())
                {
                    toAdd.MulDagger(toMul);
                }
                else
                {
                    toAdd.Mul(toMul);
                }
            }

            toAdd.MulReal(fFactorG);
            res.Add(toAdd);
        }

        //staple calculated
        deviceSU3 force(pDeviceData[linkIndex]);
        force.MulDagger(res);
        force.Ta();
        force.MulReal(betaOverN);

        //force is additive
        pForceData[linkIndex].Add(force);
    }
}


#pragma endregion


CActionGaugePlaquetteRigidAcc::CActionGaugePlaquetteRigidAcc()
    : CAction()
    , m_uiPlaqutteCount(0)
    , m_bDirichlet(FALSE)
{
}

void CActionGaugePlaquetteRigidAcc::PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate)
{
    if (0 == uiUpdateIterate)
    {
        m_fLastEnergy = EnergySingleField(FALSE, pGauge, NULL);
    }
}

void CActionGaugePlaquetteRigidAcc::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CAction::Initial(pOwner, param, byId);

    m_fBetaOverN = CCommonData::m_fBeta / static_cast<Real>(GetDefaultMatrixN());
    m_uiPlaqutteCount = _HC_Volume * (_HC_Dir - 1) * (_HC_Dir - 2);

    Real fG = 0.1f;
    param.FetchValueReal(_T("AccG"), fG);
    CCommonData::m_fG = fG;

    //TArray<INT> centerArray;
    //param.FetchValueArrayINT(_T("Center"), centerArray);
    //if (centerArray.Num() > 3)
    //{
    //    SSmallInt4 sCenter;
    //    sCenter.x = static_cast<SBYTE>(centerArray[0]);
    //    sCenter.y = static_cast<SBYTE>(centerArray[1]);
    //    sCenter.z = static_cast<SBYTE>(centerArray[2]);
    //    sCenter.w = static_cast<SBYTE>(centerArray[3]);
    //    CCommonData::m_sCenter = sCenter;
    //}
    //else
    //{
    //    CCommonData::m_sCenter.z = 0;
    //}

    INT iVaule = 1;
    param.FetchValueINT(_T("Dirichlet"), iVaule);
    m_bDirichlet = (0 != iVaule);
}

void CActionGaugePlaquetteRigidAcc::SetBeta(Real fBeta)
{
    CCommonData::m_fBeta = static_cast<DOUBLE>(fBeta);
    m_fBetaOverN = fBeta / static_cast<Real>(GetDefaultMatrixN());
}

UBOOL CActionGaugePlaquetteRigidAcc::CalculateForceOnGaugeSingleField(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
{
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);
    if (NULL == pGaugeSU3 || NULL == pForceSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteRigidAcc only work with SU3D now.\n"));
        return TRUE;
    }

    preparethread;

    _kernelAddForce4PlaqutteTermSU3_RigidAcc << <block, threads >> >(
        pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pForceSU3->m_pDeviceData, 
        m_fBetaOverNR,
        CCommonData::m_fG,
        m_bDirichlet,
        pGaugeSU3->m_byFieldId);

    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

DOUBLE CActionGaugePlaquetteRigidAcc::EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable)
{
    if (bBeforeEvolution)
    {
        return m_fLastEnergy;
    }

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    if (NULL == pGaugeSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteRotating only work with SU3-Dirichlet now.\n"));
        return m_fLastEnergy;
    }

    preparethread;
    //_kernelEnergy_RigidAcc << <block, threads >> > (
    //        pGaugeSU3->m_pDeviceData,
    //        pGaugeSU3->m_byFieldId,
    //        m_fBetaOverN,
    //        CCommonData::m_fG,
    //        _D_RealThreadBuffer);

    //Real fEnergy2 = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    _kernelEnergy_RigidAcc_Simplified << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        m_fBetaOverNR,
        CCommonData::m_fG,
        m_bDirichlet,
        _D_RealThreadBuffer);

    m_fNewEnergy = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    //appParanoiac(_T("E1 = %2.18f, E2 = %2.18f\n"), m_fNewEnergy, fEnergy2);

    return m_fNewEnergy;
}

void CActionGaugePlaquetteRigidAcc::SetG(Real fG)
{ 
    CCommonData::m_fG = fG;
}

CCString CActionGaugePlaquetteRigidAcc::GetInfos(const CCString &tab) const
{
    CCString sRet = CAction::GetInfos(tab);
    sRet = sRet + tab + _T("Beta : ") + appToString(CCommonData::m_fBeta) + _T("\n");
    sRet = sRet + tab + _T("Acc : ") + appToString(CCommonData::m_fG) + _T("\n");
    sRet = sRet + tab + _T("Dirichlet : ") + appToString(m_bDirichlet) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================