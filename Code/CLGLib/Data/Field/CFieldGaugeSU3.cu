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
__global__
void _kernelInitialSU3Feield(deviceSU3 *pDevicePtr, EFieldInitialType eInitialType)
{
    deviceSU3 id = deviceSU3::makeSU3Id();
    deviceSU3 zero = deviceSU3::makeSU3Zero();

    gaugeSU3KernelFuncionStart

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
            pDevicePtr[uiLinkIndex] = deviceSU3::makeSU3Random(_deviceGetFatIndex(coord, idir + 1));
        }
        break;
        case EFIT_RandomGenerator:
        {
            pDevicePtr[uiLinkIndex] = deviceSU3::makeSU3RandomGenerator(_deviceGetFatIndex(coord, idir + 1));
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

    gaugeSU3KernelFuncionEnd
}

__global__
void _kernelAxpySU3A(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x, _Complex a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex].MulCompC(a));

    gaugeSU3KernelFuncionEnd
}

__global__
void _kernelAxpySU3Real(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x, Real a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex].MulRealC(a));

    gaugeSU3KernelFuncionEnd
}

__global__
void _kernelAxpyPlusSU3(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex]);

    gaugeSU3KernelFuncionEnd
}

__global__
void _kernelAxpyMinusSU3(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].Sub(x[uiLinkIndex]);

    gaugeSU3KernelFuncionEnd
}

__global__
void _kernelScalarMultiplySU3Complex(deviceSU3 *pDevicePtr, _Complex a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].MulComp(a);

    gaugeSU3KernelFuncionEnd
}

__global__
void _kernelScalarMultiplySU3Real(deviceSU3 *pDevicePtr, Real a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].MulReal(a);

    gaugeSU3KernelFuncionEnd
}

/**
* debug kernel
*/
__global__ void _kernelPrintSU3(const deviceSU3 * __restrict__ pDeviceData)
{
    intokernaldir;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            UINT linkIndex = _deviceGetLinkIndex(coord, idir);
            printf("link at %d: %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i\n",
                linkIndex,
                pDeviceData[linkIndex].m_me[0].x, pDeviceData[linkIndex].m_me[0].y,
                pDeviceData[linkIndex].m_me[1].x, pDeviceData[linkIndex].m_me[1].y,
                pDeviceData[linkIndex].m_me[2].x, pDeviceData[linkIndex].m_me[2].y,
                pDeviceData[linkIndex].m_me[3].x, pDeviceData[linkIndex].m_me[3].y,
                pDeviceData[linkIndex].m_me[4].x, pDeviceData[linkIndex].m_me[4].y,
                pDeviceData[linkIndex].m_me[5].x, pDeviceData[linkIndex].m_me[5].y,
                pDeviceData[linkIndex].m_me[6].x, pDeviceData[linkIndex].m_me[6].y,
                pDeviceData[linkIndex].m_me[7].x, pDeviceData[linkIndex].m_me[7].y,
                pDeviceData[linkIndex].m_me[8].x, pDeviceData[linkIndex].m_me[8].y
            );
        }
    }
}

/**
* calculate Staple and Force At Site
*/
__global__
void _kernelStapleAtSiteSU3(
    const deviceSU3 * __restrict__ pDeviceData,
    deviceSU3 *pStapleData, //can be NULL
    deviceSU3 *pForceData,
    Real betaOverN)
{
    intokernaldir;

    betaOverN = betaOverN * F(0.5);
    SIndex plaquttes[kMaxPlaqutteCache];
    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            UINT linkIndex = _deviceGetLinkIndex(coord, idir);
            UINT uiPlaqutteCount = 0;
            UINT uiPlaqutteLength = 0;

            //int2.x is linkIndex
            //int2.y is fieldIndex (may on bounday)
            //sign of int2.y is whether inverse
            __idx->_deviceGetPlaquttesAtLink(plaquttes, uiPlaqutteCount, uiPlaqutteLength, linkIndex);
            //printf("plaqutte count = %d, length = %d\n", uiPlaqutteCount, uiPlaqutteLength);
            deviceSU3 res = deviceSU3::makeSU3Zero();

            //there are 6 staples, each is sum of two plaquttes
            for (int i = 0; i < uiPlaqutteCount; ++i)
            {
                SIndex first = plaquttes[i * (uiPlaqutteLength - 1)];
                deviceSU3 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
                if (first.NeedToDagger())
                {
                    toAdd.Dagger();
                }

                for (int j = 1; j < uiPlaqutteLength - 1; ++j)
                {
                    SIndex nextlink = plaquttes[i * (uiPlaqutteLength - 1) + j];
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
            force.Ta();
            force.MulReal(betaOverN);

            //force is additive
            pForceData[linkIndex].Add(force);
        }
    }
}

__global__
void _kernelStapleAtSiteSU3CacheIndex(
    const deviceSU3 * __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceSU3 *pStapleData, //can be NULL
    deviceSU3 *pForceData,
    Real betaOverN)
{
    intokernaldir;

    betaOverN = betaOverN * F(0.5);
    UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAll = plaqCount * plaqLengthm1;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            UINT linkIndex = _deviceGetLinkIndex(coord, idir);
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
            force.Ta();
            force.MulReal(betaOverN);

            //force is additive
            pForceData[linkIndex].Add(force);
        }
    }
}

/**
* calculate Staple and eneregy At Site
*/
__global__
void _kernelPlaqutteEnergySU3(
    const deviceSU3 * __restrict__ pDeviceData,
    Real betaOverN,
    Real* results)
{
    intokernal;

    Real resThisThread = F(0.0);
    SIndex plaquttes[kMaxPlaqutteCache];
    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;

        UINT siteIndex = _deviceGetSiteIndex(coord);
        UINT uiPlaqutteCount = 0;
        UINT uiPlaqutteLength = 0;
        __idx->_deviceGetPlaquttesAtSite(plaquttes, uiPlaqutteCount, uiPlaqutteLength, siteIndex);

        for (int i = 0; i < uiPlaqutteCount; ++i)
        {
            SIndex first = plaquttes[i * uiPlaqutteLength];

            deviceSU3 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (int j = 1; j < uiPlaqutteLength; ++j)
            {
                SIndex nextlink = plaquttes[i * uiPlaqutteLength + j];
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

#if _CLG_DEBUG
            Real reTr = toAdd.ReTr();
            assert(reTr > -F(1.50001));
            assert(reTr < F(3.00001));
#endif
            resThisThread += (F(3.0)-toAdd.ReTr());
        }
    }

    results[__thread_id] = resThisThread * betaOverN;

    //printf("  ---- energy: thread=%d, res=%f\n", __thread_id, results[__thread_id]);
}

__global__
void _kernelPlaqutteEnergySU3CacheIndex(
    const deviceSU3 * __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    Real betaOverN,
    Real* results)
{
    intokernal;

    Real resThisThread = F(0.0);
    UINT plaqCountAll = plaqCount * plaqLength;
    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;

        UINT siteIndex = _deviceGetSiteIndex(coord);

        for (int i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLength + siteIndex * plaqCountAll];
            deviceSU3 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (int j = 1; j < plaqLength; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLength + j + siteIndex * plaqCountAll];
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

#if _CLG_DEBUG
            Real reTr = toAdd.ReTr();
            assert(reTr > -F(1.50001));
            assert(reTr < F(3.00001));
#endif
            resThisThread += (F(3.0) - toAdd.ReTr());
        }
    }

    results[__thread_id] = resThisThread * betaOverN;

    //printf("  ---- energy: thread=%d, res=%f\n", __thread_id, results[__thread_id]);
}

__global__
void _kernelPlaqutteEnergyUsingStableSU3(
    const deviceSU3 * __restrict__ pDeviceData,
    const deviceSU3 * __restrict__ pStableData,
    Real betaOverN,
    Real* results)
{
    intokernaldir;

    Real resThisThread = F(0.0);
    SIndex plaquttes[kMaxPlaqutteCache];
    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            UINT linkIndex = _deviceGetLinkIndex(coord, idir);
            //For each link, there are 6 staples
            resThisThread += (F(18.0) - pDeviceData[linkIndex].MulDaggerC(pStableData[linkIndex]).ReTr());
        }
    }

    results[__thread_id] = resThisThread * betaOverN * F(0.25);

    //printf("  ---- energy: thread=%d, res=%f\n", __thread_id, results[__thread_id]);
}

/**
*
*/
__global__
void _kernelExpMultSU3(
    const deviceSU3 * __restrict__ pMyDeviceData,
    _Complex a,
    deviceSU3 *pU)
{
    intokernaldir;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            UINT linkIndex = _deviceGetLinkIndex(coord, idir);
            deviceSU3 expP = pMyDeviceData[linkIndex].Exp(a, _DC_ExpPrecision);

            expP.Mul(pU[linkIndex]);
            //expP.Norm();
            pU[linkIndex] = expP;
        }
    }
}


/**
* Trace (P^2)
*/
__global__ 
void _kernelCalculateKinematicEnergySU3(const deviceSU3 * __restrict__ pDeviceData, Real* results)
{
    intokernaldir;

    Real resThisThread = 0;
    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            UINT linkIndex = _deviceGetLinkIndex(coord, idir);
            resThisThread += pDeviceData[linkIndex].DaggerMulC(pDeviceData[linkIndex]).ReTr();
        }
    }

    results[__thread_id] = resThisThread;
}


__global__
void _kernelNormalizeSU3(deviceSU3 * pMyDeviceData)
{
    intokernaldir;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            UINT linkIndex = _deviceGetLinkIndex(coord, idir);
            pMyDeviceData[linkIndex].Norm();
        }
    }
}

__global__
void _kernelDotSU3(
    const deviceSU3 * __restrict__ pMyDeviceData, 
    const deviceSU3 * __restrict__ pOtherDeviceData,
    _Complex* result)
{
    intokernaldir;

    _Complex resThisThread = _make_cuComplex(0,0);
    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            UINT linkIndex = _deviceGetLinkIndex(coord, idir);
            resThisThread = _cuCaddf(resThisThread, pMyDeviceData[linkIndex].DaggerMulC(pOtherDeviceData[linkIndex]).Tr());
        }
    }

    result[__thread_id] = resThisThread;
}

__global__
void _kernelSetConfigurationSU3(
    deviceSU3* pDeviceData,
    const Real* __restrict__ pRealData)
{
    gaugeSU3KernelFuncionStart

    //In Bridge, it is t,z,y,x
    //x + y * nx + z * nx * ny + t * nx * ny * nz
    UINT uiBridgeSiteIndex = coord[3] * _DC_Lx * _DC_Ly * _DC_Lz + coord[2] * _DC_Lx * _DC_Ly + coord[1] * _DC_Lx + coord[0];
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

void CFieldGaugeSU3::ScalarMultply(const _Complex& a)
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

void CFieldGaugeSU3::Axpy(const _Complex& a, const CField* x)
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
    default:
        appCrucial(_T("Not supported input file type %s\n"), __ENUM_TO_STRING(EFieldFileType, eType).c_str());
        break;

    }
}

void CFieldGaugeSU3::SetByArray(Real* array)
{
    assert(NULL != array);
    //we algin the su3 now
    //assert(sizeof(deviceSU3) == 32 * sizeof(Real));

    //checkCudaErrors(cudaMemcpy(m_pDeviceData, array, sizeof(Real) * _HC_LinkCount * 18, cudaMemcpyHostToDevice));
    
    Real* pDeviceArray;
    checkCudaErrors(cudaMalloc((void**)&pDeviceArray, sizeof(Real) * _HC_LinkCount * 18));
    checkCudaErrors(cudaMemcpy(pDeviceArray, array, sizeof(Real) * _HC_LinkCount * 18, cudaMemcpyHostToDevice));
    preparethread;
    _kernelSetConfigurationSU3 << <block, threads >> > (m_pDeviceData, pDeviceArray);
    checkCudaErrors(cudaFree(pDeviceArray));

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

    if (!m_bPlaqutteIndexCached)
    {
        _kernelStapleAtSiteSU3 << <block, threads >> > (
            m_pDeviceData, 
            NULL == pStableSU3 ? NULL : pStableSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, 
            betaOverN);
    }
    else 
    {
        _kernelStapleAtSiteSU3CacheIndex << <block, threads >> > (
            m_pDeviceData, 
            m_pPlaquttesPerLink,
            m_uiPlaqutteLength,
            m_uiPlaqutteCountPerLink,
            NULL == pStableSU3 ? NULL : pStableSU3->m_pDeviceData, 
            pForceSU3->m_pDeviceData, 
            betaOverN);
    }
}

Real CFieldGaugeSU3::CalculatePlaqutteEnergy(Real betaOverN) const
{
    preparethread;
    if (!m_bPlaqutteIndexCached)
    {
        _kernelPlaqutteEnergySU3 << <block, threads >> > (
            m_pDeviceData, 
            betaOverN, 
            _D_RealThreadBuffer);
    }
    else 
    {
        _kernelPlaqutteEnergySU3CacheIndex << <block, threads >> > (
            m_pDeviceData, 
            m_pPlaquttesPerSite,
            m_uiPlaqutteLength,
            m_uiPlaqutteCountPerSite,
            betaOverN, 
            _D_RealThreadBuffer);
    }
    

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
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceData, sizeof(deviceSU3) * m_uiLinkeCount));
}

CFieldGaugeSU3::~CFieldGaugeSU3()
{
    checkCudaErrors(cudaFree(m_pDeviceData));
}

void CFieldGaugeSU3::ExpMult(const _Complex& a, CField* U) const
{
    if (NULL == U || EFT_GaugeSU3 != U->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: U field is not SU3");
        return;
    }

    CFieldGaugeSU3* pUField = dynamic_cast<CFieldGaugeSU3*>(U);

    preparethread;
    _kernelExpMultSU3 << < block, threads >> > (m_pDeviceData, a, pUField->m_pDeviceData);
}

void CFieldGaugeSU3::ElementNormalize()
{
    preparethread;
    _kernelNormalizeSU3 << < block, threads >> > (m_pDeviceData);
}

_Complex CFieldGaugeSU3::Dot(const CField* other) const
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



__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================