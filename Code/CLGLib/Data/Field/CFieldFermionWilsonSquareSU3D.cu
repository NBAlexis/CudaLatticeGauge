//=============================================================================
// FILENAME : CFieldFermionWilsonSquareSU3D.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [05/18/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"


__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionWilsonSquareSU3D)

#pragma region DOperator

#pragma region kernel

/**
* Dw phi(x) = phi(x) - kai sum _mu (1-gamma _mu) U(x,mu) phi(x+ mu) + (1+gamma _mu) U^{dagger}(x-mu) phi(x-mu)
* U act on su3
* gamma act on spinor
*
* If bDagger, it is gamma5, D, gamma5
*
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_D(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;
    BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex & sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        pResultData[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3(); 
        //((CFieldBoundaryWilsonSquareSU3*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[__idx->_devcieExchangeBoundaryFieldSiteIndex(sIdx)];
        return;
    }

    const gammaMatrix& gamma5 = __chiralGamma[GAMMA5];
    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];
    if (bDDagger)
    {
        pResultData[uiSiteIndex] = gamma5.MulWilsonC(pResultData[uiSiteIndex]);
    }

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        const gammaMatrix& gammaMu = __chiralGamma[GAMMA1 + idir];

        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        //deviceSU3 x_Gauge_element = pGauge[linkIndex];
        deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir);
        //deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
        if (x_m_mu_Gauge.NeedToDagger())
        {
            x_m_mu_Gauge_element.Dagger();
        }

        //deviceWilsonVectorSU3 x_p_mu_Fermion_element = pDeviceData[x_p_mu_Fermion.m_uiSiteIndex];
        //deviceWilsonVectorSU3 x_m_mu_Fermion_element = pDeviceData[x_m_mu_Fermion.m_uiSiteIndex];
        deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_p_mu_Fermion, byFieldId);
        deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_m_mu_Fermion, byFieldId);

        if (bDDagger)
        {
            x_p_mu_Fermion_element = gamma5.MulWilsonC(x_p_mu_Fermion_element);
            x_m_mu_Fermion_element = gamma5.MulWilsonC(x_m_mu_Fermion_element);
        }

        //hopping terms

        //U(x,mu) phi(x+ mu)
        deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            //printf("Opposite x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
            result.Sub(u_phi_x_p_m);

            //- gammamu U(x,mu) phi(x+ mu)
            result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
        }
        else
        {
            result.Add(u_phi_x_p_m);

            //- gammamu U(x,mu) phi(x+ mu)
            result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            result.Sub(u_dagger_phi_x_m_m);

            //gammamu U^{dagger}(x-mu) phi(x-mu)
            result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        }
        else
        {
            result.Add(u_dagger_phi_x_m_m);

            //gammamu U^{dagger}(x-mu) phi(x-mu)
            result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        }
    }

    //result = phi(x) - kai sum _mu result
    result.MulReal(kai);
    pResultData[uiSiteIndex].Sub(result);

    if (bDDagger)
    {
        pResultData[uiSiteIndex] = gamma5.MulWilsonC(pResultData[uiSiteIndex]);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        pResultData[uiSiteIndex].MulReal(fCoeff);
        break;
    case EOCT_Complex:
        pResultData[uiSiteIndex].MulComp(cCoeff);
        break;
    }
}

/**
* The output is on a gauge field
* Therefor cannot make together with _kernelDWilson
*
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDWilsonForceSU3_D(
    const deviceWilsonVectorSU3* __restrict__ pInverseD,
    const deviceWilsonVectorSU3* __restrict__ pInverseDDdagger,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    deviceSU3* pForce,
    Real fKai,
    BYTE byFieldId)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sSite = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    const deviceWilsonVectorSU3 x_Left(_deviceGetFermionBCWilsonSU3(pInverseDDdagger, sSite, byFieldId));
    const deviceWilsonVectorSU3 x_Right(_deviceGetFermionBCWilsonSU3(pInverseD, sSite, byFieldId));

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        const SIndex& x_p_mu_Fermion = pFermionMove[linkIndex * 2];

        //If one of the sites is on surface, it has no contribution.
        //Note that, the bond on surface is equivelant to both sites on surface.

        if (//__idx->_deviceIsBondOnSurface(uiBigIdx, idir) ||
            x_p_mu_Fermion.IsDirichlet()
         || sSite.IsDirichlet())
        {
            //continue;
        }
        else
        {
            //Get Gamma mu
            const gammaMatrix& gammaMu = __chiralGamma[GAMMA1 + idir];

            //SIndex x_m_mu_Gauge = __idx->_deviceGaugeIndexWalk(uiSiteIndex, -(idir + 1));
             // __idx->_deviceFermionIndexWalk(byFieldId, uiSiteIndex, (idir + 1));

            //all not on surface
            const deviceWilsonVectorSU3& x_p_mu_Right = pInverseD[x_p_mu_Fermion.m_uiSiteIndex];
            deviceWilsonVectorSU3 x_p_mu_Left(pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex]);
            //deviceWilsonVectorSU3 x_p_mu_Right = _deviceGetFermionBCWilsonSU3(pInverseD, x_p_mu_Fermion, byFieldId);
            //deviceWilsonVectorSU3 x_p_mu_Left = _deviceGetFermionBCWilsonSU3(pInverseDDdagger, x_p_mu_Fermion, byFieldId);

            deviceSU3 x_Gauge_element = pGauge[linkIndex]; // _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir); //pGauge[linkIndex];

            deviceWilsonVectorSU3 right1(x_p_mu_Right);
            right1.Sub(gammaMu.MulWilsonC(right1));
            deviceSU3 mid = deviceSU3::makeSU3Contract(x_Left, right1);

            deviceWilsonVectorSU3 right2(x_Right);
            right2.Add(gammaMu.MulWilsonC(right2));
            mid.Add(deviceSU3::makeSU3Contract(right2, x_p_mu_Left));

            deviceSU3 forceOfThisLink = x_Gauge_element.MulC(mid);
            forceOfThisLink.Ta();
            if (x_p_mu_Fermion.NeedToOpposite())
            {
                forceOfThisLink.MulReal(-fKai);
            }
            else
            {
                forceOfThisLink.MulReal(fKai);
            }

            pForce[linkIndex].Add(forceOfThisLink);
        }
    }
}

#pragma endregion

void CFieldFermionWilsonSquareSU3D::DOperator(void* pTargetBuffer, const void* pBuffer, 
    const void* pGaugeBuffer, 
    UBOOL bDagger, EOperatorCoefficientType eOCT, 
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    deviceWilsonVectorSU3* pTarget = (deviceWilsonVectorSU3*)pTargetBuffer;
    const deviceWilsonVectorSU3* pSource = (deviceWilsonVectorSU3*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;

    preparethread;
    _kernelDFermionWilsonSquareSU3_D << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        pTarget, 
        m_fKai, 
        m_byFieldId, 
        bDagger,
        eOCT,
        fRealCoeff, 
        cCmpCoeff);

}

void CFieldFermionWilsonSquareSU3D::DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const
{
    deviceSU3* pForceSU3 = (deviceSU3*)pForce;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;
    const deviceWilsonVectorSU3* pDphiBuffer = (deviceWilsonVectorSU3*)pDphi;
    const deviceWilsonVectorSU3* pDDphiBuffer = (deviceWilsonVectorSU3*)pDDphi;

    preparethread;
    _kernelDWilsonForceSU3_D << <block, threads >> > (
        pDphiBuffer,
        pDDphiBuffer,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        pForceSU3,
        m_fKai, m_byFieldId);
}

#pragma endregion

#pragma region Kernel

/**
* Initialize
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialFermionWilsonSquareSU3ForHMC(
    deviceWilsonVectorSU3 *pDevicePtr, 
    BYTE byFieldId)
{
    intokernalInt4;

    const UINT bigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][bigIdx];
    if (sIdx.IsDirichlet())
    {
        pDevicePtr[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3(); //((CFieldBoundaryWilsonSquareSU3*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[__idx->_devcieExchangeBoundaryFieldSiteIndex(sIdx)];
    }
    else
    {
        pDevicePtr[uiSiteIndex] = deviceWilsonVectorSU3::makeRandomGaussian(_deviceGetFatIndex(uiSiteIndex, 0));
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFixBoundaryWilsonSU3(deviceWilsonVectorSU3* pDeviceData, BYTE byFieldId)
{
    intokernalInt4;

    const SIndex idx = __idx->_deviceGetMappingIndex(sSite4, byFieldId);
    if (idx.IsDirichlet())
    {
        pDeviceData[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3(); //((CFieldBoundaryWilsonSquareSU3*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[__idx->_devcieExchangeBoundaryFieldSiteIndex(idx)];
    }
}

#pragma endregion


/**
* generate phi by gaussian random.
* phi = D phi
*/
void CFieldFermionWilsonSquareSU3D::PrepareForHMC(const CFieldGauge* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3 * pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldFermionWilsonSquareSU3* pPooled = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    preparethread;
    _kernelInitialFermionWilsonSquareSU3ForHMC << <block, threads >> > (
        pPooled->m_pDeviceData,
        m_byFieldId);

    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, 
        FALSE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));

    pPooled->Return();

    if (NULL != appGetFermionSolver(m_byFieldId) && !appGetFermionSolver(m_byFieldId)->IsAbsoluteAccuracy())
    {
        m_fLength = Dot(this).x;
    }
    //cache a inverse DDdagger field
    if (CCommonData::m_bStoreLastSolution)
    {
        CFieldCache* pCache = appGetLattice()->m_pFieldCache;
        CField* pField = pCache->GetCachedField(CFieldCache::CachedInverseDDdaggerField);
        if (NULL == pField)
        {
            pField = GetCopy();
            pCache->CacheField(CFieldCache::CachedInverseDDdaggerField, pField);
        }
        else
        {
            CopyTo(pField);
        }
    }
}

void CFieldFermionWilsonSquareSU3D::PrepareForHMCOnlyRandomize()
{
    preparethread;
    _kernelInitialFermionWilsonSquareSU3ForHMC << <block, threads >> > (
        m_pDeviceData,
        m_byFieldId);
}

void CFieldFermionWilsonSquareSU3D::PrepareForHMCNotRandomize(const CFieldGauge* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldFermionWilsonSquareSU3* pPooled = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    preparethread;
    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceWilsonVectorSU3) * _HC_Volume, cudaMemcpyDeviceToDevice));

    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData,
        FALSE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));

    pPooled->Return();

    if (NULL != appGetFermionSolver(m_byFieldId) && !appGetFermionSolver(m_byFieldId)->IsAbsoluteAccuracy())
    {
        m_fLength = Dot(this).x;
    }
    //cache a inverse DDdagger field
    if (CCommonData::m_bStoreLastSolution)
    {
        CFieldCache* pCache = appGetLattice()->m_pFieldCache;
        CField* pField = pCache->GetCachedField(CFieldCache::CachedInverseDDdaggerField);
        if (NULL == pField)
        {
            pField = GetCopy();
            pCache->CacheField(CFieldCache::CachedInverseDDdaggerField, pField);
        }
        else
        {
            CopyTo(pField);
        }
    }
}

void CFieldFermionWilsonSquareSU3D::FixBoundary()
{
    appDetailed(_T("CFieldFermionWilsonSquareSU3D::FixBoundary()\n"));

    preparethread;
    _kernelFixBoundaryWilsonSU3 << <block, threads >> >(m_pDeviceData, m_byFieldId);
}

void CFieldFermionWilsonSquareSU3D::CopyTo(CField* U) const
{
    CFieldFermionWilsonSquareSU3::CopyTo(U);
}

CCString CFieldFermionWilsonSquareSU3D::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = tab + _T("Name : CFieldFermionWilsonSquareSU3D\n");
    sRet = sRet + tab + _T("Hopping : ") + appFloatToString(CCommonData::m_fKai) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================