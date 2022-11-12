//=============================================================================
// FILENAME : CMeasureBerryPhase.cpp
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [10/13/2020 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region kernels

/**
 * The block is usual block * Lx, Ly, Lz
 * The thread is usual thread
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelMomentumFieldWilsonDiracSU3(
    const deviceWilsonVectorSU3* __restrict__ pInverseD, 
    deviceWilsonVectorSU3* pRes,
    SSmallInt4 xprime)
{
    const UINT blockIdxX1 = blockIdx.x / _DC_Lx;
    const UINT blockIdxX2 = blockIdx.x % _DC_Lx;
    const UINT blockIdxY1 = blockIdx.y / _DC_Ly;
    const UINT blockIdxY2 = blockIdx.y % _DC_Ly;
    const UINT blockIdxZ1 = blockIdx.z / _DC_Lz;
    const UINT blockIdxZ2 = blockIdx.z % _DC_Lz;

    const UINT uiSiteIndexP = ((threadIdx.x + blockIdxX1 * blockDim.x) * _DC_GridDimZT 
                             + (threadIdx.y + blockIdxY1 * blockDim.y) * _DC_Lt 
                             + (threadIdx.z + blockIdxZ1 * blockDim.z));

    const UINT uiSiteIndexX = (blockIdxX2 * _DC_Ly + blockIdxY2) * _DC_GridDimZT
        + blockIdxZ2 * _DC_Lt + (threadIdx.z + blockIdxZ1 * blockDim.z);

    SSmallInt4 p = __deviceSiteIndexToInt4(uiSiteIndexP);
    SSmallInt4 x = __deviceSiteIndexToInt4(uiSiteIndexX);

    //printf("p:(%d,%d,%d,%d) - x:(%d,%d,%d,%d)\n",
    //    p.x, p.y, p.z, p.w,
    //    x.x, x.y, x.z, x.w
    //    );

    //result will be add into uiSiteIndexP
    //We will calculate exp(I p.(x - xprime)) * D^-1(x,t|xprime,0)
    deviceWilsonVectorSU3 thisSite = pInverseD[uiSiteIndexX];
    //Fox x, y and z, there is no anti-periodic boundary condition, so 
    const Real px = PI2 * (p.x + F(1.0) - _DC_Lx / F(2.0)) / static_cast<Real>(_DC_Lx);
    const Real py = PI2 * (p.y + F(1.0) - _DC_Ly / F(2.0)) / static_cast<Real>(_DC_Ly);
    const Real pz = PI2 * (p.z + F(1.0) - _DC_Lz / F(2.0)) / static_cast<Real>(_DC_Lz);
    const Real fDot = px * static_cast<Real>(x.x - xprime.x)
                    + py * static_cast<Real>(x.y - xprime.y)
                    + pz * static_cast<Real>(x.z - xprime.z);
    thisSite.MulComp(_make_cuComplex(_cos(fDot), _sin(fDot)));

    #pragma unroll
    for (BYTE byElement = 0; byElement < 12; ++byElement)
    {
        atomicAdd(&pRes[uiSiteIndexP].m_me[byElement].x, thisSite.m_me[byElement].x);
        atomicAdd(&pRes[uiSiteIndexP].m_me[byElement].y, thisSite.m_me[byElement].y);
    }
}


__global__ void _CLG_LAUNCH_BOUND
_kernelBerryConnectWilsonDiracSU3(
    BYTE byFieldId,
    const deviceWilsonVectorSU3* __restrict__ pMomentumField,
    CLGComplex* pRes)
{
    intokernalInt4;

    #pragma unroll
    for (BYTE dir = 0; dir < 3; ++dir)
    {
        const UINT n_p_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(_deviceSmallInt4OffsetC(sSite4, dir + 1))].m_uiSiteIndex;
        CLGComplex dot = pMomentumField[uiSiteIndex].ConjugateDotC(pMomentumField[n_p_mu]);
        Real fArg = __cuCargf(dot);
        pRes[_deviceGetLinkIndex(uiSiteIndex, dir)] = _make_cuComplex(_cos(fArg), _sin(fArg));
    }
}


__global__ void _CLG_LAUNCH_BOUND
_kernelMomentumFieldKSSU3(
    const deviceSU3Vector* __restrict__ pInverseD,
    deviceSU3Vector* pRes,
    SSmallInt4 xprime)
{
    const UINT blockIdxX1 = blockIdx.x / _DC_Lx;
    const UINT blockIdxX2 = blockIdx.x % _DC_Lx;
    const UINT blockIdxY1 = blockIdx.y / _DC_Ly;
    const UINT blockIdxY2 = blockIdx.y % _DC_Ly;
    const UINT blockIdxZ1 = blockIdx.z / _DC_Lz;
    const UINT blockIdxZ2 = blockIdx.z % _DC_Lz;

    const UINT uiSiteIndexP = ((threadIdx.x + blockIdxX1 * blockDim.x) * _DC_GridDimZT
        + (threadIdx.y + blockIdxY1 * blockDim.y) * _DC_Lt
        + (threadIdx.z + blockIdxZ1 * blockDim.z));

    const UINT uiSiteIndexX = (blockIdxX2 * _DC_Ly + blockIdxY2) * _DC_GridDimZT
        + blockIdxZ2 * _DC_Lt + (threadIdx.z + blockIdxZ1 * blockDim.z);

    SSmallInt4 p = __deviceSiteIndexToInt4(uiSiteIndexP);
    SSmallInt4 x = __deviceSiteIndexToInt4(uiSiteIndexX);

    //printf("p:(%d,%d,%d,%d) - x:(%d,%d,%d,%d)\n",
    //    p.x, p.y, p.z, p.w,
    //    x.x, x.y, x.z, x.w
    //    );

    //result will be add into uiSiteIndexP
    //We will calculate exp(I p.(x - xprime)) * D^-1(x,t|xprime,0)
    deviceSU3Vector thisSite = pInverseD[uiSiteIndexX];
    //Fox x, y and z, there is no anti-periodic boundary condition, so 
    const Real px = PI2 * (p.x + F(1.0) - _DC_Lx / F(2.0)) / static_cast<Real>(_DC_Lx);
    const Real py = PI2 * (p.y + F(1.0) - _DC_Ly / F(2.0)) / static_cast<Real>(_DC_Ly);
    const Real pz = PI2 * (p.z + F(1.0) - _DC_Lz / F(2.0)) / static_cast<Real>(_DC_Lz);
    const Real fDot = px * static_cast<Real>(x.x - xprime.x)
                    + py * static_cast<Real>(x.y - xprime.y)
                    + pz * static_cast<Real>(x.z - xprime.z);
    thisSite.MulComp(_make_cuComplex(_cos(fDot), _sin(fDot)));

    #pragma unroll
    for (BYTE byElement = 0; byElement < 3; ++byElement)
    {
        atomicAdd(&pRes[uiSiteIndexP].m_ve[byElement].x, thisSite.m_ve[byElement].x);
        atomicAdd(&pRes[uiSiteIndexP].m_ve[byElement].y, thisSite.m_ve[byElement].y);
    }
}

/**
 * A_mu = phi^+(p) phi(p+mu) / |phi^+(p) phi(p+mu)|
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelBerryConnectKSSU3(
    BYTE byFieldId,
    const deviceSU3Vector* __restrict__ pMomentumField,
    CLGComplex* pRes)
{
    intokernalInt4;

    #pragma unroll
    for (BYTE dir = 0; dir < 3; ++dir)
    {
        const UINT n_p_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(_deviceSmallInt4OffsetC(sSite4, dir + 1))].m_uiSiteIndex;
        CLGComplex dot = pMomentumField[uiSiteIndex].ConjugateDotC(pMomentumField[n_p_mu]);
        Real fArg = __cuCargf(dot);
        pRes[_deviceGetLinkIndex(uiSiteIndex, dir)] = _make_cuComplex(_cos(fArg), _sin(fArg));
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelBerryCurvatureU1(
    BYTE byGaugeFieldId,
    const CLGComplex* __restrict__ pU1Field,
    BYTE byT,
#if !_CLG_DOUBLEFLOAT
    DOUBLE* pRes
#else
    Real* pRes
#endif
)
{
    SSmallInt4 sSite4; 
    const UINT _ixy = (threadIdx.x + blockIdx.x * blockDim.x); 
    sSite4.x = static_cast<SBYTE> (_ixy / _DC_Ly); 
    sSite4.y = static_cast<SBYTE> (_ixy % _DC_Ly); 
    sSite4.z = static_cast<SBYTE>(threadIdx.y + blockIdx.y * blockDim.y); 
    sSite4.w = byT;
    const UINT uiSiteSpatial = _ixy * _DC_Lz + sSite4.z;

    Real fArgSum = F(0.0);
    #pragma unroll
    for (BYTE dir = 0; dir < 3; ++dir)
    {
        //it is xy, yz, zt and tz plaqutte
        const INT iDir1 = ((dir + 1) % 3) + 1;
        const INT iDir2 = ((dir + 2) % 3) + 1;
        const INT path[4] = { iDir1, iDir2, -iDir1, -iDir2 };
        const Real curve = _deviceLinkU1ArgSum(pU1Field, sSite4, 4, byGaugeFieldId, path);
        fArgSum += curve;// *curve; // __cuCargf(pU1Field[_deviceGetSiteIndex(sSite4) * 4 + dir])* curve;
    }

    pRes[uiSiteSpatial] = fArgSum;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelBerryCurvatureU1XY(
    BYTE byGaugeFieldId,
    const CLGComplex* __restrict__ pU1Field,
    BYTE byT,
    INT dir1, INT dir2,
#if !_CLG_DOUBLEFLOAT
    DOUBLE* pRes
#else
    Real* pRes
#endif
)
{
    SSmallInt4 sSite4;
    const UINT _ixy = (threadIdx.x + blockIdx.x * blockDim.x);
    sSite4.x = static_cast<SBYTE> (_ixy / _DC_Ly);
    sSite4.y = static_cast<SBYTE> (_ixy % _DC_Ly);
    sSite4.z = static_cast<SBYTE>(threadIdx.y + blockIdx.y * blockDim.y);
    sSite4.w = byT;
    const UINT uiSiteSpatial = _ixy * _DC_Lz + sSite4.z;

    const INT path[4] = { dir1, dir2, -dir1, -dir2 };
    pRes[uiSiteSpatial] = _deviceLinkU1ArgSum(pU1Field, sSite4, 4, byGaugeFieldId, path);
}

#pragma endregion

__CLGIMPLEMENT_CLASS(CMeasureBerryPhase)

void CMeasureBerryPhase::CalculateMomentumSpacePhiWilsonDiracForPoint(const SSmallInt4& xprime, const CFieldGaugeSU3* pGauge)
{
    SFermionSource source;
    source.m_eSourceType = EFS_Point;
    source.m_byColorIndex = 4;
    source.m_sSourcePoint = xprime;

    const dim3 block(_HC_DecompX * _HC_Lx, _HC_DecompY * _HC_Ly, _HC_DecompZ * _HC_Lz);
    const dim3 threads(_HC_DecompLx, _HC_DecompLy, _HC_DecompLz);
    CFieldFermionWilsonSquareSU3* sourcefield = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(m_pMomentumField->m_byFieldId));
    CFieldFermionWilsonSquareSU3* momentumfield = dynamic_cast<CFieldFermionWilsonSquareSU3*>(m_pMomentumField);
    sourcefield->InitialAsSource(source);
    sourcefield->InverseD(pGauge);

    _kernelMomentumFieldWilsonDiracSU3 << <block, threads >> > (
        sourcefield->m_pDeviceData,
        momentumfield->m_pDeviceData,
        xprime
        );

    sourcefield->Return();
}

void CMeasureBerryPhase::CalculateMomentumSpacePhiWilsonDirac(const CFieldGaugeSU3* pGauge)
{
    m_pMomentumField->InitialField(EFIT_Zero);

    //Use transilational invarience, we only calculate xprime = 0 (otherwise, we will have run over x')
    SSmallInt4 xprime(0, 0, 0, 0);
    CalculateMomentumSpacePhiWilsonDiracForPoint(xprime, pGauge);
}

void CMeasureBerryPhase::CalculateU1FieldWilsonDirac()
{
    CFieldFermionWilsonSquareSU3* momentumfield = dynamic_cast<CFieldFermionWilsonSquareSU3*>(m_pMomentumField);
    m_pU1Field->InitialField(EFIT_Identity);
    preparethread;
    _kernelBerryConnectWilsonDiracSU3 << <block, threads>> > (
        momentumfield->m_byFieldId, 
        momentumfield->m_pDeviceData, 
        m_pU1Field->m_pDeviceData);

    //m_pU1Field->DebugPrintMe();
}


void CMeasureBerryPhase::CalculateMomentumSpacePhiKSForPoint(const SSmallInt4& xprime, const CFieldGaugeSU3* pGauge)
{
    SFermionSource source;
    source.m_eSourceType = EFS_Point;
    source.m_byColorIndex = 4;
    source.m_sSourcePoint = xprime;

    const dim3 block(_HC_DecompX * _HC_Lx, _HC_DecompY * _HC_Ly, _HC_DecompZ * _HC_Lz);
    const dim3 threads(_HC_DecompLx, _HC_DecompLy, _HC_DecompLz);
    CFieldFermionKSSU3* sourcefield = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(m_pMomentumField->m_byFieldId));
    CFieldFermionKSSU3* momentumfield = dynamic_cast<CFieldFermionKSSU3*>(m_pMomentumField);
    sourcefield->InitialAsSource(source);
    sourcefield->InverseD(pGauge);

    _kernelMomentumFieldKSSU3 << <block, threads >> > (
        sourcefield->m_pDeviceData,
        momentumfield->m_pDeviceData,
        xprime
        );

    sourcefield->Return();
}

void CMeasureBerryPhase::CalculateMomentumSpacePhiKS(const CFieldGaugeSU3* pGauge)
{
    m_pMomentumField->InitialField(EFIT_Zero);

    //Use transilational invarience, we only calculate xprime = 0 (otherwise, we will have run over x')
    SSmallInt4 xprime(0, 0, 0, 0);
    CalculateMomentumSpacePhiKSForPoint(xprime, pGauge);
}

void CMeasureBerryPhase::CalculateU1FieldKS()
{
    CFieldFermionKSSU3* momentumfield = dynamic_cast<CFieldFermionKSSU3*>(m_pMomentumField);
    m_pU1Field->InitialField(EFIT_Identity);
    preparethread;
    _kernelBerryConnectKSSU3 << <block, threads >> > (
        momentumfield->m_byFieldId,
        momentumfield->m_pDeviceData,
        m_pU1Field->m_pDeviceData);
}

void CMeasureBerryPhase::CalculateBerryPhase(BYTE byGaugeFieldId)
{
    const dim3 block(_HC_DecompX, _HC_DecompY, 1); 
    const dim3 threads(_HC_DecompLx, _HC_DecompLy, 1);

    TArray<DOUBLE> res;
    TArray<DOUBLE> resXY;
    TArray<DOUBLE> resXZ;
    TArray<DOUBLE> resXT;
    TArray<DOUBLE> resYZ;
    TArray<DOUBLE> resYT;
    TArray<DOUBLE> resZT;

    for (BYTE byT = 0; byT < _HC_Lt; ++byT)
    {
        _kernelBerryCurvatureU1 << <block, threads >> > (
            byGaugeFieldId,
            m_pU1Field->m_pDeviceData,
            byT,
            _D_RealThreadBuffer
            );

        res.AddItem(appGetCudaHelper()->ReduceReal(_D_RealThreadBuffer, _HC_Volume_xyz));

        _kernelBerryCurvatureU1XY << <block, threads >> > (
            byGaugeFieldId,
            m_pU1Field->m_pDeviceData,
            byT,
            1, 2,
            _D_RealThreadBuffer
            );

        resXY.AddItem(appGetCudaHelper()->ReduceReal(_D_RealThreadBuffer, _HC_Volume_xyz));

        _kernelBerryCurvatureU1XY << <block, threads >> > (
            byGaugeFieldId,
            m_pU1Field->m_pDeviceData,
            byT,
            1, 3,
            _D_RealThreadBuffer
            );

        resXZ.AddItem(appGetCudaHelper()->ReduceReal(_D_RealThreadBuffer, _HC_Volume_xyz));

        _kernelBerryCurvatureU1XY << <block, threads >> > (
            byGaugeFieldId,
            m_pU1Field->m_pDeviceData,
            byT,
            1, 4,
            _D_RealThreadBuffer
            );

        resXT.AddItem(appGetCudaHelper()->ReduceReal(_D_RealThreadBuffer, _HC_Volume_xyz));

        _kernelBerryCurvatureU1XY << <block, threads >> > (
            byGaugeFieldId,
            m_pU1Field->m_pDeviceData,
            byT,
            2, 3,
            _D_RealThreadBuffer
            );

        resYZ.AddItem(appGetCudaHelper()->ReduceReal(_D_RealThreadBuffer, _HC_Volume_xyz));

        _kernelBerryCurvatureU1XY << <block, threads >> > (
            byGaugeFieldId,
            m_pU1Field->m_pDeviceData,
            byT,
            2, 4,
            _D_RealThreadBuffer
            );

        resYT.AddItem(appGetCudaHelper()->ReduceReal(_D_RealThreadBuffer, _HC_Volume_xyz));

        _kernelBerryCurvatureU1XY << <block, threads >> > (
            byGaugeFieldId,
            m_pU1Field->m_pDeviceData,
            byT,
            3, 4,
            _D_RealThreadBuffer
            );

        resZT.AddItem(appGetCudaHelper()->ReduceReal(_D_RealThreadBuffer, _HC_Volume_xyz));
    }

    m_lstData.AddItem(res);
    m_lstDataXY.AddItem(resXY);
    m_lstDataXZ.AddItem(resXZ);
    m_lstDataXT.AddItem(resXT);
    m_lstDataYZ.AddItem(resYZ);
    m_lstDataYT.AddItem(resYT);
    m_lstDataZT.AddItem(resZT);
}

void CMeasureBerryPhase::AllocateBuffers()
{
    m_pU1Field = new CFieldGaugeU1();
    m_pMomentumField = dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
}

CMeasureBerryPhase::~CMeasureBerryPhase()
{
    appSafeDelete(m_pU1Field);
    appSafeDelete(m_pGaugeFixing);
    m_pMomentumField->Return();
}

void CMeasureBerryPhase::Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);
    INT iValue = 1;
    param.FetchValueINT(_T("WilsonDirac"), iValue);
    m_bWilsonDirac = (0 != iValue);

    iValue = 0;
    param.FetchValueINT(_T("DoGaugeFixing"), iValue);
    m_bGuageFixing = (0 != iValue);

    iValue = 0;
    param.FetchValueINT(_T("ShowRes"), iValue);
    m_bShowRes = (0 != iValue);

    iValue = 2;
    param.FetchValueINT(_T("FieldId"), iValue);
    m_byFieldId = static_cast<BYTE>(iValue);

    AllocateBuffers();
}

void CMeasureBerryPhase::OnConfigurationAccepted(const CFieldGauge* pAcceptGauge, const CFieldGauge* )
{
    if (m_bGuageFixing)
    {
        if (NULL == m_pGaugeFixing)
        {
            m_pGaugeFixing = dynamic_cast<CFieldGaugeSU3*>(pAcceptGauge->GetCopy());
        }
        else
        {
            pAcceptGauge->CopyTo(m_pGaugeFixing);
        }
        if (NULL != appGetLattice()->m_pGaugeFixing)
        {
            appGetLattice()->m_pGaugeFixing->GaugeFixing(m_pGaugeFixing);
        }

        if (m_bWilsonDirac)
        {
            CalculateMomentumSpacePhiWilsonDirac(m_pGaugeFixing);
            CalculateU1FieldWilsonDirac();
        }
        else
        {
            CalculateMomentumSpacePhiKS(m_pGaugeFixing);
            CalculateU1FieldKS();
        }
    }
    else
    {
        if (m_bWilsonDirac)
        {
            CalculateMomentumSpacePhiWilsonDirac(dynamic_cast<const CFieldGaugeSU3*>(pAcceptGauge));
            CalculateU1FieldWilsonDirac();
        }
        else
        {
            CalculateMomentumSpacePhiKS(dynamic_cast<const CFieldGaugeSU3*>(pAcceptGauge));
            CalculateU1FieldKS();
        }
    }

    CalculateBerryPhase(pAcceptGauge->m_byFieldId);

    //Gather result
    if (m_bShowRes)
    {
        appSetLogDate(FALSE);
        appGeneral(_T("Berry phase: {"));

        for (INT t = 0; t < _HC_Lti; ++t)
        {
            if (m_bShowRes)
            {
                appGeneral(_T("%2.18f%s"),
                    m_lstData[m_uiConfigurationCount][t],
                    t == (_HC_Lti - 1) ? _T("}") : _T(", ")
                );
            }
        }

        appGeneral(_T("\n"));

        appGeneral(_T("Berry phaseXY: {"));

        for (INT t = 0; t < _HC_Lti; ++t)
        {
            if (m_bShowRes)
            {
                appGeneral(_T("%2.18f%s"),
                    m_lstDataXY[m_uiConfigurationCount][t],
                    t == (_HC_Lti - 1) ? _T("}") : _T(", ")
                );
            }
        }

        appGeneral(_T("\n"));

        appGeneral(_T("Berry phaseZT: {"));

        for (INT t = 0; t < _HC_Lti; ++t)
        {
            if (m_bShowRes)
            {
                appGeneral(_T("%2.18f%s"),
                    m_lstDataZT[m_uiConfigurationCount][t],
                    t == (_HC_Lti - 1) ? _T("}") : _T(", ")
                );
            }
        }

        appGeneral(_T("\n"));
        appSetLogDate(TRUE);
    }

    ++m_uiConfigurationCount;
}

void CMeasureBerryPhase::Average(UINT )
{

}

void CMeasureBerryPhase::Report()
{

}

void CMeasureBerryPhase::Reset()
{
    m_lstData.Reset();
    m_lstDataXY.Reset();
    m_lstDataXZ.Reset();
    m_lstDataXT.Reset();
    m_lstDataYZ.Reset();
    m_lstDataYT.Reset();
    m_lstDataZT.Reset();
    m_uiConfigurationCount = 0;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================