//=============================================================================
// FILENAME : CActionGaugePlaquettePolar.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [05/01/2024 nbale]
//=============================================================================
#include "CLGLib_Private.h"


__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquettePolar)


#pragma region kernels

/**
* 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelEnergy_Polar(
    const deviceSU3 * __restrict__ pDeviceData,
    BYTE byFieldId,
    const Real* __restrict__ betalist,
    Real fRIn,
    Real fDeltaR,
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

    const Real r = fRIn + fDeltaR * sSite4.x;
    const Real oneOverRSq = F(1.0) / (r * r);
    const Real betaOverN = betalist[sSite4.x];

    deviceSU3 toSub(_deviceClover(pDeviceData, sSite4, uiBigIdx, 1, 0, byFieldId));
    toSub.Add(_deviceClover(pDeviceData, sSite4, uiBigIdx, 1, 2, byFieldId));
    toSub.Add(_deviceClover(pDeviceData, sSite4, uiBigIdx, 1, 3, byFieldId));
    toSub.MulReal(oneOverRSq); //Now this is (1/r^2)(Upx + Upz + Upt)

    toSub.Add(_deviceClover(pDeviceData, sSite4, uiBigIdx, 0, 2, byFieldId));
    toSub.Add(_deviceClover(pDeviceData, sSite4, uiBigIdx, 0, 3, byFieldId));
    toSub.Add(_deviceClover(pDeviceData, sSite4, uiBigIdx, 2, 3, byFieldId));
    //Now this is Urz + Urt + Uzt + (1/r^2)(Upx + Upz + Upt)

    const Real toAdd = F(9.0) * r + F(9.0) / r - toSub.ReTr() * F(0.25) * r;
    //Now this is ReTr[r(3 - Urz - Urt - Uzt) + (1/r) (3 - U14 - U24 - U34)]

    results[uiSiteIndex] = betaOverN * toAdd;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelEnergy_Polar_Simplified(
    const deviceSU3* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
    BYTE plaqLength, BYTE plaqCount,
    const Real* __restrict__ betalist,
    Real fRIn,
    Real fDeltaR,
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
        deviceSU3 toAdd(_deviceGetGaugeBCSU3(pDeviceData, first));

        if (first.NeedToDagger())
        {
            toAdd.Dagger();
        }

        first = pCachedIndex[i * plaqLength + 1 + uiSiteIndex * plaqCountAll];
        const BYTE nu = first.m_byDir;
        deviceSU3 toMul(_deviceGetGaugeBCSU3(pDeviceData, first));
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
            toMul = _deviceGetGaugeBCSU3(pDeviceData, first);
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
        if (0 == i || 3 == i || 4 == i)
        {
            resThisThread += (F(3.0) - toAdd.ReTr()) * _deviceGnPolarPhiLeft(sSite4, betalist, fRIn, fDeltaR, mu, nu, bDirichlet);
        }
        else
        {
            resThisThread += (F(3.0) - toAdd.ReTr()) * _deviceGnPolarSpatialLeft(sSite4, betalist, fRIn, fDeltaR, mu, nu, bDirichlet);
        }
    }

    results[uiSiteIndex] = resThisThread;
}

/**
* Assuming x-y-t direction is NOT Dirichlet
* If \bar{U}_{mu nu} does not have Z direction, g(n) = f(n)
* If \bar{U}_{mu nu} does have Z direction, g(n) = [f(n) + f(n+Z)]/2
* This can be Dirichlet, because we have links like U_z(nz=0)
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForce4PlaqutteTermSU3_Polar(
    const deviceSU3* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
    BYTE plaqLength, BYTE plaqCount,
    deviceSU3* pForceData,
    const Real* __restrict__ betalist,
    Real fRIn,
    Real fDeltaR,
    UBOOL bDirichlet,
    BYTE byFieldId)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    //const SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];

    //Real test_force = F(0.0);
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

            deviceSU3 toAdd(_deviceGetGaugeBCSU3(pDeviceData, first));
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
                    fFactorG = _deviceGnPolarPhiRight(sSite4, betalist, fRIn, fDeltaR, mu, nu, bDirichlet);
                }
                else
                {
                    fFactorG = _deviceGnPolarSpatialRight(sSite4, betalist, fRIn, fDeltaR, mu, nu, bDirichlet);
                }
            }
            else
            {
                if (mu == 3 || nu == 3) //one of mu nu is t
                {
                    fFactorG = _deviceGnPolarPhiLeft(sSite4, betalist, fRIn, fDeltaR, mu, nu, bDirichlet);
                }
                else
                {
                    fFactorG = _deviceGnPolarSpatialLeft(sSite4, betalist, fRIn, fDeltaR, mu, nu, bDirichlet);
                }
            }

            for (BYTE j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                deviceSU3 toMul(_deviceGetGaugeBCSU3(pDeviceData, nextlink));

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
        force.MulReal(F(-0.5));

        //force is additive
        pForceData[linkIndex].Add(force);
    }
}


#pragma endregion


CActionGaugePlaquettePolar::CActionGaugePlaquettePolar()
    : CAction()
    , m_uiPlaqutteCount(0)
    , m_bDirichlet(FALSE)
    , m_fRIn(F(1.0))
    , m_fROut(F(2.0))
    , m_fDeltaR(F(0.1))
    , m_deviceBetaList(NULL)
{
    checkCudaErrors(cudaMalloc((void**)&m_deviceBetaList, sizeof(Real) * _HC_Lx));

    for (UINT i = 0; i < _HC_Lx; ++i)
    {
        m_lstBeta.AddItem(F(1.5));
    }

    checkCudaErrors(cudaMemcpy(m_deviceBetaList, m_lstBeta.GetData(), sizeof(Real) * _HC_Lx, cudaMemcpyHostToDevice));
}

CActionGaugePlaquettePolar::~CActionGaugePlaquettePolar()
{
    checkCudaErrors(cudaFree(m_deviceBetaList));
}

void CActionGaugePlaquettePolar::PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate)
{
    if (0 == uiUpdateIterate)
    {
        m_fLastEnergy = EnergySingleField(FALSE, pGauge, NULL);
    }
}

void CActionGaugePlaquettePolar::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CAction::Initial(pOwner, param, byId);

    m_lstBeta.RemoveAll();
    param.FetchValueArrayReal(_T("BetaList"), m_lstBeta);
    if (0 == m_lstBeta.Num())
    {
        m_lstBeta.AddItem(F(1.5));
    }
    for (INT i = 1; i < _HC_Lxi; ++i)
    {
        if (m_lstBeta.Num() < i + 1)
        {
            m_lstBeta.AddItem(m_lstBeta[i - 1]);
        }
    }
    assert(m_lstBeta.Num() == _HC_Lxi);
    checkCudaErrors(cudaMemcpy(m_deviceBetaList, m_lstBeta.GetData(), sizeof(Real) * _HC_Lx, cudaMemcpyHostToDevice));
    
    m_uiPlaqutteCount = _HC_Volume * (_HC_Dir - 1) * (_HC_Dir - 2);

    param.FetchValueReal(_T("RIn"), m_fRIn);
    param.FetchValueReal(_T("ROut"), m_fROut);

    assert(m_fRIn > _CLG_FLT_EPSILON);
    assert(m_fROut - m_fRIn > _CLG_FLT_EPSILON);
    assert(_HC_Lx > 1);

    if (m_fRIn < _CLG_FLT_EPSILON)
    {
        m_fRIn = F(1.0);
    }

    if (m_fROut - m_fRIn < _CLG_FLT_EPSILON)
    {
        m_fROut = m_fRIn + F(1.0);
    }
    m_fDeltaR = (m_fROut - m_fRIn) / (_HC_Lx - 1);

    INT iVaule = 1;
    param.FetchValueINT(_T("Dirichlet"), iVaule);
    m_bDirichlet = (0 != iVaule);
}

UBOOL CActionGaugePlaquettePolar::CalculateForceOnGaugeSingleField(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
{
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);
    if (NULL == pGaugeSU3 || NULL == pForceSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteRigidAcc only work with SU3D now.\n"));
        return TRUE;
    }

    preparethread;

    _kernelAddForce4PlaqutteTermSU3_Polar << <block, threads >> >(
        pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pForceSU3->m_pDeviceData, 
        m_deviceBetaList,
        m_fRIn,
        m_fDeltaR,
        m_bDirichlet,
        pGaugeSU3->m_byFieldId);

    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

DOUBLE CActionGaugePlaquettePolar::EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable)
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
    _kernelEnergy_Polar << <block, threads >> > (
            pGaugeSU3->m_pDeviceData,
            pGaugeSU3->m_byFieldId,
            m_deviceBetaList,
            m_fRIn,
            m_fDeltaR,
            _D_RealThreadBuffer);

    Real fEnergy2 = static_cast<Real>(appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer));

    _kernelEnergy_Polar_Simplified << <block, threads >> > (
        pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        m_deviceBetaList,
        m_fRIn,
        m_fDeltaR,
        m_bDirichlet,
        _D_RealThreadBuffer);

    m_fNewEnergy = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    appParanoiac(_T("E1 = %2.18f, E2 = %2.18f\n"), m_fNewEnergy, fEnergy2);

    return m_fNewEnergy;
}

void CActionGaugePlaquettePolar::SetBetaList(const TArray<Real>& betalst)
{ 
    m_lstBeta = betalst;
    for (INT i = 0; i < m_lstBeta.Num(); ++i)
    {
        m_lstBeta[i] = m_lstBeta[i] / static_cast<Real>(_HC_SUN);
    }
    if (0 == m_lstBeta.Num())
    {
        m_lstBeta.AddItem(F(1.5));
    }
    for (INT i = 1; i < _HC_Lxi; ++i)
    {
        if (m_lstBeta.Num() < i + 1)
        {
            m_lstBeta.AddItem(m_lstBeta[i - 1]);
        }
    }
    assert(m_lstBeta.Num() == _HC_Lxi);
    checkCudaErrors(cudaMemcpy(m_deviceBetaList, m_lstBeta.GetData(), sizeof(Real) * _HC_Lx, cudaMemcpyHostToDevice));
}

void CActionGaugePlaquettePolar::SetR(Real fIn, Real fOut)
{
    assert(m_fRIn > _CLG_FLT_EPSILON);
    assert(m_fROut - m_fRIn > _CLG_FLT_EPSILON);
    assert(_HC_Lx > 1);

    if (m_fRIn < _CLG_FLT_EPSILON)
    {
        m_fRIn = F(1.0);
    }

    if (m_fROut - m_fRIn < _CLG_FLT_EPSILON)
    {
        m_fROut = m_fRIn + F(1.0);
    }
    m_fDeltaR = (m_fROut - m_fRIn) / (_HC_Lx - 1);
}

CCString CActionGaugePlaquettePolar::GetInfos(const CCString &tab) const
{
    CCString sRet = tab + _T("Name : CActionGaugePlaquetteRigidAcc\n");
    sRet = sRet + CAction::GetInfos(tab);

    sRet = sRet + tab + _T("Beta : [") +  + _T("]\n");
    for (UINT i = 0; i < _HC_Lx; ++i)
    {
        sRet = sRet + appAnyToString(m_lstBeta[i] * _HC_SUN) + _T(", ");
    }
    sRet = sRet + _T("]\n");
    sRet = sRet + tab + _T("R : [");
    for (UINT i = 0; i < _HC_Lx; ++i)
    {
        sRet = sRet + appAnyToString(m_fRIn + i * m_fDeltaR) + _T(", ");
    }
    sRet = sRet + _T("]\n");
    sRet = sRet + tab + _T("Dirichlet : ") + appAnyToString(m_bDirichlet) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================