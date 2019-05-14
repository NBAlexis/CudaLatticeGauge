//=============================================================================
// FILENAME : CFieldGaugeSU3.cu
// 
// DESCRIPTION:
// This is the device implementations of gauge SU3
//
// The SU3 Matrix is
// 0 1 2
// 3 4 5
// 6 7 8
//
// Number of threads: < 1024
// Number of blocks: V / 1024
//
// threadIdx.xyz = xyz, and we loop for t and dir
//
// REVISION:
//  [12/4/2018 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldGaugeSU3)

#pragma region Kernels

/**
* Initial SU3 Field with a value
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialSU3Feield(deviceSU3 *pDevicePtr, EFieldInitialType eInitialType)
{
    deviceSU3 id = deviceSU3::makeSU3Id();
    deviceSU3 zero = deviceSU3::makeSU3Zero();

    intokernalInt4;
    BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        switch (eInitialType)
        {
        case EFIT_Zero:
        {
            pDevicePtr[uiLinkIndex] = zero;
        }
        break;
        case EFIT_Identity:
        {
            pDevicePtr[uiLinkIndex] = id;
        }
        break;
        case EFIT_Random:
        {
            pDevicePtr[uiLinkIndex] = deviceSU3::makeSU3Random(_deviceGetFatIndex(uiSiteIndex, idir + 1));
        }
        break;
        case EFIT_RandomGenerator:
        {
            if (__idx->m_pDeviceIndexPositionToSIndex[1]
                [__idx->_deviceGetBigIndex(sSite4)].IsDirichlet())
            {
                pDevicePtr[uiLinkIndex] = zero;
            }
            else
            {
                pDevicePtr[uiLinkIndex] = deviceSU3::makeSU3RandomGenerator(_deviceGetFatIndex(uiSiteIndex, idir + 1));
            }
        }
        break;
        case EFIT_SumGenerator:
        {
            pDevicePtr[uiLinkIndex] = deviceSU3::makeSU3SumGenerator(F(1.0));
        }
        break;
        default:
        {
            printf("SU3 Field cannot be initialized with this type!");
        }
        break;
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpySU3A(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x, CLGComplex a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex].MulCompC(a));

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpySU3Real(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x, Real a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex].MulRealC(a));

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyPlusSU3(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex]);

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyMinusSU3(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].Sub(x[uiLinkIndex]);

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplySU3Complex(deviceSU3 *pDevicePtr, CLGComplex a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].MulComp(a);

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplySU3Real(deviceSU3 *pDevicePtr, Real a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].MulReal(a);

    gaugeSU3KernelFuncionEnd
}

/**
* debug kernel
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelPrintSU3(const deviceSU3 * __restrict__ pDeviceData)
{
    gaugeSU3KernelFuncionStart;

    printf("link at %d: %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i\n",
        uiLinkIndex,
        pDeviceData[uiLinkIndex].m_me[0].x, pDeviceData[uiLinkIndex].m_me[0].y,
        pDeviceData[uiLinkIndex].m_me[1].x, pDeviceData[uiLinkIndex].m_me[1].y,
        pDeviceData[uiLinkIndex].m_me[2].x, pDeviceData[uiLinkIndex].m_me[2].y,
        pDeviceData[uiLinkIndex].m_me[3].x, pDeviceData[uiLinkIndex].m_me[3].y,
        pDeviceData[uiLinkIndex].m_me[4].x, pDeviceData[uiLinkIndex].m_me[4].y,
        pDeviceData[uiLinkIndex].m_me[5].x, pDeviceData[uiLinkIndex].m_me[5].y,
        pDeviceData[uiLinkIndex].m_me[6].x, pDeviceData[uiLinkIndex].m_me[6].y,
        pDeviceData[uiLinkIndex].m_me[7].x, pDeviceData[uiLinkIndex].m_me[7].y,
        pDeviceData[uiLinkIndex].m_me[8].x, pDeviceData[uiLinkIndex].m_me[8].y
    );

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelStapleAtSiteSU3CacheIndex(
    const deviceSU3 * __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceSU3 *pStapleData, //can be NULL
    deviceSU3 *pForceData,
    Real betaOverN)
{
    intokernaldir;

    //Real test_force = F(0.0);
    betaOverN = betaOverN * F(-0.5);
    UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAll = plaqCount * plaqLengthm1;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceSU3 res = deviceSU3::makeSU3Zero();

        //there are 6 staples, each is sum of two plaquttes
        for (int i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            deviceSU3 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);

            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (int j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                deviceSU3 toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);

                if (nextlink.NeedToDagger())
                {
                    toAdd.MulDagger(toMul);
                }
                else
                {
                    toAdd.Mul(toMul);
                }
            }
            res.Add(toAdd);
        }
        if (NULL != pStapleData)
        {
            pStapleData[linkIndex] = res;
        }

        //staple calculated
        deviceSU3 force(pDeviceData[linkIndex]);
        force.MulDagger(res);
        //test_force += F(-2.0) * betaOverN * __SU3Generators[8].MulC(force).ImTr();
        force.Ta();
        force.MulReal(betaOverN);

        //force is additive
        pForceData[linkIndex].Add(force);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateOnlyStaple(
    const deviceSU3 * __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceSU3 *pStapleData)
{
    intokernaldir;

    //Real test_force = F(0.0);
    UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAll = plaqCount * plaqLengthm1;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceSU3 res = deviceSU3::makeSU3Zero();

        //there are 6 staples, each is sum of two plaquttes
        for (int i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            deviceSU3 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);

            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (int j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                deviceSU3 toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);
                if (nextlink.NeedToDagger())
                {
                    toAdd.MulDagger(toMul);
                }
                else
                {
                    toAdd.Mul(toMul);
                }
            }
            res.Add(toAdd);
        }
        pStapleData[linkIndex] = res;
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergySU3CacheIndex(
    const deviceSU3 * __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    Real betaOverN,
    Real* results)
{
    intokernal;

    Real resThisThread = F(0.0);
    UINT plaqCountAll = plaqCount * plaqLength;
    for (BYTE i = 0; i < plaqCount; ++i)
    {
        SIndex first = pCachedIndex[i * plaqLength + uiSiteIndex * plaqCountAll];
        deviceSU3 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
        if (first.NeedToDagger())
        {
            toAdd.Dagger();
        }

        for (BYTE j = 1; j < plaqLength; ++j)
        {
            first = pCachedIndex[i * plaqLength + j + uiSiteIndex * plaqCountAll];
            deviceSU3 toMul(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
            if (first.NeedToDagger())
            {
                toAdd.MulDagger(toMul);
            }
            else
            {
                toAdd.Mul(toMul);
            }
        }

#if _CLG_DEBUG
        Real reTr = toAdd.ReTr();
        assert(reTr > -F(1.50001));
        assert(reTr < F(3.00001));
#endif
        resThisThread += (F(3.0) - toAdd.ReTr());
    }

    results[uiSiteIndex] = resThisThread * betaOverN;

    //printf("  ---- energy: thread=%d, res=%f\n", __thread_id, results[__thread_id]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyUsingStableSU3(
    const deviceSU3 * __restrict__ pDeviceData,
    const deviceSU3 * __restrict__ pStableData,
    Real betaOverN,
    Real* results)
{
    intokernaldir;

    Real resThisThread = F(0.0);
    SIndex plaquttes[kMaxPlaqutteCache];
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        //For each link, there are 6 staples
        resThisThread += (F(18.0) - pDeviceData[linkIndex].MulDaggerC(pStableData[linkIndex]).ReTr());
    }

    results[uiSiteIndex] = resThisThread * betaOverN * F(0.25);

    //printf("  ---- energy: thread=%d, res=%f\n", __thread_id, results[__thread_id]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelExpMultSU3RealQ(
    const deviceSU3 * __restrict__ pMyDeviceData,
    Real a,
    deviceSU3 *pU)
{
    intokernaldir;

    for (BYTE idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        //deviceSU3 expP = pMyDeviceData[linkIndex].ExpReal(a, _DC_ExpPrecision);
        deviceSU3 expP = pMyDeviceData[linkIndex].QuickExp(a);
        expP.Mul(pU[linkIndex]);
        pU[linkIndex] = expP;
    }
}
__global__ void _CLG_LAUNCH_BOUND
_kernelExpMultSU3Real(
    const deviceSU3 * __restrict__ pMyDeviceData,
    Real a,
    deviceSU3 *pU,
    BYTE prec)
{
    intokernaldir;

    for (BYTE idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        //deviceSU3 expP = pMyDeviceData[linkIndex].ExpReal(a, _DC_ExpPrecision);
        deviceSU3 expP = pMyDeviceData[linkIndex].ExpReal(a, prec);
        expP.Mul(pU[linkIndex]);
        pU[linkIndex] = expP;
    }
}

/**
* Trace (P^2)
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateKinematicEnergySU3(const deviceSU3 * __restrict__ pDeviceData, Real* results)
{
    intokernaldir;

    Real resThisThread = F(0.0);
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        resThisThread += pDeviceData[linkIndex].DaggerMulC(pDeviceData[linkIndex]).ReTr();
    }

    results[uiSiteIndex] = resThisThread;
}


__global__ void _CLG_LAUNCH_BOUND
_kernelNormalizeSU3(deviceSU3 * pMyDeviceData)
{
    intokernaldir;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        pMyDeviceData[linkIndex].Norm();
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDotSU3(
    const deviceSU3 * __restrict__ pMyDeviceData, 
    const deviceSU3 * __restrict__ pOtherDeviceData,
    CLGComplex* result)
{
    intokernaldir;

    CLGComplex resThisThread = _make_cuComplex(0,0);
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        resThisThread = _cuCaddf(resThisThread, pMyDeviceData[linkIndex].DaggerMulC(pOtherDeviceData[linkIndex]).Tr());
    }

    result[uiSiteIndex] = resThisThread;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelSetConfigurationSU3(
    deviceSU3* pDeviceData,
    const Real* __restrict__ pRealData)
{
    gaugeSU3KernelFuncionStart

        //In Bridge, it is t,z,y,x
        //x + y * nx + z * nx * ny + t * nx * ny * nz
    SSmallInt4 xyzt = __deviceSiteIndexToInt4(uiSiteIndex);
    UINT uiBridgeSiteIndex = xyzt.w * _DC_Lx * _DC_Ly * _DC_Lz + xyzt.z * _DC_Lx * _DC_Ly + xyzt.y * _DC_Lx + xyzt.x;
    UINT uiBridgeLinkIndex = uiBridgeSiteIndex * _DC_Dir + idir;

    //0 1 2
    //3 4 5
    //6 7 8
    pDeviceData[uiLinkIndex].m_me[0] = _make_cuComplex(pRealData[18 * uiBridgeLinkIndex +  0], pRealData[18 * uiBridgeLinkIndex +  1]);
    pDeviceData[uiLinkIndex].m_me[1] = _make_cuComplex(pRealData[18 * uiBridgeLinkIndex +  2], pRealData[18 * uiBridgeLinkIndex +  3]);
    pDeviceData[uiLinkIndex].m_me[2] = _make_cuComplex(pRealData[18 * uiBridgeLinkIndex +  4], pRealData[18 * uiBridgeLinkIndex +  5]);
    pDeviceData[uiLinkIndex].m_me[3] = _make_cuComplex(pRealData[18 * uiBridgeLinkIndex +  6], pRealData[18 * uiBridgeLinkIndex +  7]);
    pDeviceData[uiLinkIndex].m_me[4] = _make_cuComplex(pRealData[18 * uiBridgeLinkIndex +  8], pRealData[18 * uiBridgeLinkIndex +  9]);
    pDeviceData[uiLinkIndex].m_me[5] = _make_cuComplex(pRealData[18 * uiBridgeLinkIndex + 10], pRealData[18 * uiBridgeLinkIndex + 11]);
    pDeviceData[uiLinkIndex].m_me[6] = _make_cuComplex(pRealData[18 * uiBridgeLinkIndex + 12], pRealData[18 * uiBridgeLinkIndex + 13]);
    pDeviceData[uiLinkIndex].m_me[7] = _make_cuComplex(pRealData[18 * uiBridgeLinkIndex + 14], pRealData[18 * uiBridgeLinkIndex + 15]);
    pDeviceData[uiLinkIndex].m_me[8] = _make_cuComplex(pRealData[18 * uiBridgeLinkIndex + 16], pRealData[18 * uiBridgeLinkIndex + 17]);

    //pDeviceData[uiLinkIndex].DebugPrint();
    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFixBoundarySU3(deviceSU3 * pDeviceData)
{
    intokernalInt4;

    SIndex idx = __idx->_deviceGetMappingIndex(sSite4, 1);
    if (idx.IsDirichlet())
    {
        UINT uiDir = _DC_Dir;

        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            pDeviceData[_deviceGetLinkIndex(uiSiteIndex, idir)] = ((CFieldBoundaryGaugeSU3*)__boundaryFieldPointers[1])->m_pDeviceData
            [
                __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idx.m_byDir
            ];
        }

        //printf("%d, %d, %d, %d\n", sSite4.x, sSite4.y, sSite4.z, sSite4.w);
    }
}

#pragma endregion

void CFieldGaugeSU3::AxpyPlus(const CField* x)
{
    if (NULL == x || EFT_GaugeSU3 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeSU3* pSU3x = dynamic_cast<const CFieldGaugeSU3*>(x);
    preparethread;
    _kernelAxpyPlusSU3 << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData);
}

void CFieldGaugeSU3::AxpyMinus(const CField* x)
{
    if (NULL == x || EFT_GaugeSU3 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeSU3* pSU3x = dynamic_cast<const CFieldGaugeSU3*>(x);
    preparethread;
    _kernelAxpyMinusSU3 << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData);

}

void CFieldGaugeSU3::ScalarMultply(const CLGComplex& a)
{
    preparethread;
    _kernelScalarMultiplySU3Complex << <block, threads >> > (m_pDeviceData, a);
}

void CFieldGaugeSU3::ScalarMultply(Real a)
{
    preparethread;
    _kernelScalarMultiplySU3Real << <block, threads >> > (m_pDeviceData, a);
}

void CFieldGaugeSU3::Axpy(Real a, const CField* x)
{
    if (NULL == x || EFT_GaugeSU3 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeSU3* pSU3x = dynamic_cast<const CFieldGaugeSU3*>(x);
    preparethread;
    _kernelAxpySU3Real << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData, a);
}

void CFieldGaugeSU3::Axpy(const CLGComplex& a, const CField* x)
{
    if (NULL == x || EFT_GaugeSU3 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeSU3* pSU3x = dynamic_cast<const CFieldGaugeSU3*>(x);
    preparethread;
    _kernelAxpySU3A << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData, a);
}


void CFieldGaugeSU3::Zero()
{
    preparethread;
    _kernelInitialSU3Feield << <block, threads >> > (m_pDeviceData, EFIT_Zero);
}

void CFieldGaugeSU3::Indentity()
{
    preparethread;
    _kernelInitialSU3Feield << <block, threads >> > (m_pDeviceData, EFIT_Identity);
}

void CFieldGaugeSU3::MakeRandomGenerator()
{
    preparethread;
    _kernelInitialSU3Feield << <block, threads >> > (m_pDeviceData, EFIT_RandomGenerator);
}

/**
*
*/
void CFieldGaugeSU3::InitialField(EFieldInitialType eInitialType)
{
    preparethread;
    _kernelInitialSU3Feield << <block, threads >> > (m_pDeviceData, eInitialType);
}

void CFieldGaugeSU3::InitialFieldWithFile(const CCString& sFileName, EFieldFileType eType)
{
    switch (eType)
    {
    case EFFT_BridgePPTXT:
    {
        CCString sContent = appGetFileSystem()->ReadAllText(sFileName);
        TArray<INT> seps;
        seps.AddItem(_T('\n'));
        seps.AddItem(_T('\r'));
        TArray<CCString> sStringlist = appGetStringList(sContent, seps, 0x7fffffff);
        assert(static_cast<UINT>(sStringlist.Num()) == _HC_LinkCount * 18);

        Real* pData = (Real*)malloc(sizeof(Real) * sStringlist.Num());
        for (INT i = 0; i < sStringlist.Num(); ++i)
        {
            pData[i] = static_cast<Real>(appStrToDOUBLE(sStringlist[i]));
        }

        SetByArray(pData);
    }
    break;
    case EFFT_BridgePPBin:
    {
        UINT uiSize = 0;
        BYTE* allBytes = appGetFileSystem()->ReadAllBytes(sFileName, uiSize);
        assert(uiSize == 8 * 18 * _HC_LinkCount);
        Real* pData = (Real*)malloc(sizeof(Real) * 18 * _HC_LinkCount);
        for (UINT i = 0; i < 18 * _HC_LinkCount; ++i)
        {
            BYTE data[8];
            for (INT j = 0; j < 8; ++j)
            {
                data[j] = allBytes[i * 8 + (7 - j)];
            }
            DOUBLE * dbData = (DOUBLE*)data;
            pData[i] = static_cast<Real>(*dbData);
        }
        free(allBytes);

        SetByArray(pData);
    }
    break;
    case EFFT_CLGBin:
    {
        UINT uiSize = static_cast<UINT>(sizeof(FLOAT) * 18 * m_uiLinkeCount);
        BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        InitialWithByte(data);
        free(data);
    }
    break;
    default:
        appCrucial(_T("Not supported input file type %s\n"), __ENUM_TO_STRING(EFieldFileType, eType).c_str());
        break;

    }
}

void CFieldGaugeSU3::InitialWithByte(BYTE* byData)
{
    deviceSU3* readData = (deviceSU3*)malloc(sizeof(deviceSU3) * m_uiLinkeCount);
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        for (UINT j = 0; j < 16; ++j)
        {
            FLOAT oneLink[18];
            memcpy(oneLink, byData + sizeof(FLOAT) * 18 * i, sizeof(FLOAT) * 18);
            if (j < 9)
            {
                readData[i].m_me[j] =
                    _make_cuComplex(
                        static_cast<Real>(oneLink[2 * j]),
                        static_cast<Real>(oneLink[2 * j + 1]));
            }
            else
            {
                readData[i].m_me[j] = _make_cuComplex(F(0.0), F(0.0));
            }
        }
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(deviceSU3) * m_uiLinkeCount, cudaMemcpyHostToDevice));
    free(readData);
}

void CFieldGaugeSU3::SetByArray(Real* array)
{
    assert(NULL != array);
    //we algin the su3 now
    //assert(sizeof(deviceSU3) == 32 * sizeof(Real));

    //checkCudaErrors(cudaMemcpy(m_pDeviceData, array, sizeof(Real) * _HC_LinkCount * 18, cudaMemcpyHostToDevice));
    
    Real* pDeviceArray;
    checkCudaErrors(__cudaMalloc((void**)&pDeviceArray, sizeof(Real) * _HC_LinkCount * 18));
    checkCudaErrors(cudaMemcpy(pDeviceArray, array, sizeof(Real) * _HC_LinkCount * 18, cudaMemcpyHostToDevice));
    preparethread;
    _kernelSetConfigurationSU3 << <block, threads >> > (m_pDeviceData, pDeviceArray);
    checkCudaErrors(__cudaFree(pDeviceArray));

    free(array);

    ElementNormalize();
}

/**
* (1) calculate staples
* (2) calculate force(additive)
* (3) calculate energy
*/
void CFieldGaugeSU3::CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStable, Real betaOverN) const
{
    if (NULL == pForce || EFT_GaugeSU3 != pForce->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: force field is not SU3");
        return;
    }
    if (NULL != pStable && EFT_GaugeSU3 != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: stape field is not SU3");
        return;
    }

    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);
    CFieldGaugeSU3* pStableSU3 = NULL == pStable ? NULL : dynamic_cast<CFieldGaugeSU3*>(pStable);

    preparethread;

    assert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache);

    _kernelStapleAtSiteSU3CacheIndex << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        NULL == pStableSU3 ? NULL : pStableSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData,
        betaOverN);
}

void CFieldGaugeSU3::CalculateOnlyStaple(CFieldGauge* pStable) const
{
    if (NULL == pStable || EFT_GaugeSU3 != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: stable field is not SU3");
        return;
    }
    CFieldGaugeSU3* pStableSU3 = dynamic_cast<CFieldGaugeSU3*>(pStable);

    preparethread;
    _kernelCalculateOnlyStaple << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStableSU3->m_pDeviceData);
}

Real CFieldGaugeSU3::CalculatePlaqutteEnergy(Real betaOverN) const
{
    assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);

    preparethread;
    _kernelPlaqutteEnergySU3CacheIndex << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

Real CFieldGaugeSU3::CalculatePlaqutteEnergyUsingStable(Real betaOverN, const CFieldGauge *pStable) const
{
    if (NULL == pStable || EFT_GaugeSU3 != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: stape field is not SU3");
        return F(0.0);
    }
    const CFieldGaugeSU3* pStableSU3 = dynamic_cast<const CFieldGaugeSU3*>(pStable);

    preparethread;
    _kernelPlaqutteEnergyUsingStableSU3 << <block, threads >> > (
        m_pDeviceData, 
        pStableSU3->m_pDeviceData, 
        betaOverN, 
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

Real CFieldGaugeSU3::CalculateKinematicEnergy() const
{
    preparethread;
    _kernelCalculateKinematicEnergySU3 << <block, threads >> > (m_pDeviceData, _D_RealThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

CFieldGaugeSU3::CFieldGaugeSU3() : CFieldGauge()
{
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceData, sizeof(deviceSU3) * m_uiLinkeCount));
}

CFieldGaugeSU3::~CFieldGaugeSU3()
{
    checkCudaErrors(__cudaFree(m_pDeviceData));
}

void CFieldGaugeSU3::ExpMult(Real a, CField* U) const
{
    if (NULL == U || EFT_GaugeSU3 != U->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: U field is not SU3");
        return;
    }

    CFieldGaugeSU3* pUField = dynamic_cast<CFieldGaugeSU3*>(U);

    preparethread;
    if (0 == _HC_ExpPrecision)
    {
        _kernelExpMultSU3RealQ << < block, threads >> > (m_pDeviceData, a, pUField->m_pDeviceData);
    }
    else
    {
        _kernelExpMultSU3Real << < block, threads >> > (m_pDeviceData, a, pUField->m_pDeviceData, static_cast<BYTE>(_HC_ExpPrecision));
    }
    
}

void CFieldGaugeSU3::ElementNormalize()
{
    preparethread;
    _kernelNormalizeSU3 << < block, threads >> > (m_pDeviceData);
}

CLGComplex CFieldGaugeSU3::Dot(const CField* other) const
{
    if (NULL == other || EFT_GaugeSU3 != other->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: U field is not SU3");
        return _make_cuComplex(0,0);
    }

    const CFieldGaugeSU3* pUField = dynamic_cast<const CFieldGaugeSU3*>(other);

    preparethread;
    _kernelDotSU3 << < block, threads >> > (m_pDeviceData, pUField->m_pDeviceData, _D_ComplexThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
}

void CFieldGaugeSU3::FixBoundary()
{
    appDetailed(_T("CFieldGaugeSU3::FixBoundary()\n"));

    preparethread;
    _kernelFixBoundarySU3 << <block, threads >> > (m_pDeviceData);
}

void CFieldGaugeSU3::CopyTo(CField* pTarget) const
{
    if (NULL == pTarget || EFT_GaugeSU3 != pTarget->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: target field is not SU3");
        return;
    }

    CFieldGauge::CopyTo(pTarget);

    CFieldGaugeSU3* pTargetField = dynamic_cast<CFieldGaugeSU3*>(pTarget);
    checkCudaErrors(cudaMemcpy(pTargetField->m_pDeviceData, m_pDeviceData, sizeof(deviceSU3) * m_uiLinkeCount, cudaMemcpyDeviceToDevice));
}

void CFieldGaugeSU3::DebugPrintMe() const
{
    preparethread;
    _kernelPrintSU3 << < block, threads >> > (m_pDeviceData);
}

void CFieldGaugeSU3::SaveToFile(const CCString &fileName) const
{
    UINT uiSize = 0;
    BYTE* byToSave = CopyDataOut(uiSize);
    appGetFileSystem()->WriteAllBytes(fileName.c_str(), byToSave, uiSize);
    free(byToSave);
}

BYTE* CFieldGaugeSU3::CopyDataOut(UINT &uiSize) const
{
    deviceSU3* toSave = (deviceSU3*)malloc(sizeof(deviceSU3) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceSU3) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    //fuck ofstream
    uiSize = static_cast<UINT>(sizeof(FLOAT) * m_uiLinkeCount * 18);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        FLOAT oneLink[18];
        for (UINT j = 0; j < 9; ++j)
        {
            oneLink[2 * j] = static_cast<FLOAT>(toSave[i].m_me[j].x);
            oneLink[2 * j + 1] = static_cast<FLOAT>(toSave[i].m_me[j].y);
        }
        memcpy(byToSave + i * sizeof(FLOAT) * 18, oneLink, sizeof(FLOAT) * 18);
    }
    free(toSave);

    return byToSave;
}

CCString CFieldGaugeSU3::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = tab + _T("Name : CFieldGaugeSU3\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================