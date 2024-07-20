//=============================================================================
// FILENAME : CFieldGaugeLinkDirichlet.cu
// 
// DESCRIPTION:
//
// REVISION:
//  [07/06/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CFieldGaugeLinkDirichlet.h"

__BEGIN_NAMESPACE

#pragma region Kernels

/**
* Initial SU3 Field with a value
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialGeneratorT_D(deviceGauge *pDevicePtr, BYTE byFieldId)
{
    deviceGauge zero = _makeZero<deviceGauge>();

    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        if (__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, idir))
        {
            pDevicePtr[uiLinkIndex] = zero;
        }
        else
        {
            pDevicePtr[uiLinkIndex] = _makeGaussian<deviceGauge>(_deviceGetFatIndex(uiSiteIndex, idir + 1));
        }
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelStapleAtSiteCacheIndexT_D(
    BYTE byFieldId,
    const deviceGauge * __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceGauge *pStapleData, //can be NULL
    deviceGauge *pForceData,
    Real betaOverN)
{
    intokernaldir;

    //Real test_force = F(0.0);
    betaOverN = betaOverN * F(-0.5);
    const UINT plaqLengthm1 = plaqLength - 1;
    const UINT plaqCountAll = plaqCount * plaqLengthm1;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceGauge res = _makeZero<deviceGauge>();

        //there are 6 staples, each is sum of two plaquttes
        for (BYTE i = 0; i < plaqCount; ++i)
        {
            BYTE diricCount = 0;
            const SIndex& first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            if (first.IsDirichlet())
            {
                ++diricCount;
            }
            //deviceGauge toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
            deviceGauge toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

            if (first.NeedToDagger())
            {
                _dagger(toAdd);
            }

            for (BYTE j = 1; j < plaqLengthm1; ++j)
            {
                const SIndex& nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                if (nextlink.IsDirichlet())
                {
                    ++diricCount;
                }
                //deviceGauge toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);
                deviceGauge toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink));

                if (nextlink.NeedToDagger())
                {
                    _muldag(toAdd, toMul);
                }
                else
                {
                    _mul(toAdd, toMul);
                }
            }
            if (diricCount < plaqLength - 1)
            {
                // If more than 3(including 3) of the edges are Dirichlet, 
                // the plaqutte dose NOT exist.
                _add(res, toAdd);
            }
        }
        if (NULL != pStapleData)
        {
            pStapleData[linkIndex] = res;
        }

        //staple calculated
        deviceGauge force(pDeviceData[linkIndex]);
        _muldag(force, res);
        //test_force += F(-2.0) * betaOverN * __SU3Generators[8].MulC(force).ImTr();
        _ta(force);
        _mul(force, betaOverN);

        //force is additive
        _add(pForceData[linkIndex], force);
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyCacheIndexT_D(
    BYTE byFieldId,
    const deviceGauge * __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    DOUBLE betaOverN,
    DOUBLE* results
)
{
    intokernal;

    DOUBLE resThisThread = 0.0;
    UINT plaqCountAll = plaqCount * plaqLength;
    for (BYTE i = 0; i < plaqCount; ++i)
    {
        SIndex first = pCachedIndex[i * plaqLength + uiSiteIndex * plaqCountAll];
        //deviceGauge toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
        deviceGauge toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

        if (first.NeedToDagger())
        {
            _dagger(toAdd);
        }

        for (BYTE j = 1; j < plaqLength; ++j)
        {
            first = pCachedIndex[i * plaqLength + j + uiSiteIndex * plaqCountAll];
            //deviceGauge toMul(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
            deviceGauge toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
            if (first.NeedToDagger())
            {
                _muldag(toAdd, toMul);
            }
            else
            {
                _mul(toAdd, toMul);
            }
        }

        resThisThread += (static_cast<DOUBLE>(_dim<deviceGauge>()) - _retr(toAdd));
    }

    results[uiSiteIndex] = resThisThread * betaOverN;
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateOnlyStapleT_D(
    BYTE byFieldId,
    const deviceGauge * __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceGauge *pStapleData)
{
    intokernaldir;

    //Real test_force = F(0.0);
    const UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAll = plaqCount * plaqLengthm1;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceGauge res = _makeZero<deviceGauge>();

        //there are 6 staples, each is sum of two plaquttes
        for (int i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            deviceGauge toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

            if (first.NeedToDagger())
            {
                _dagger(toAdd);
            }

            for (int j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                deviceGauge toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink));

                if (nextlink.NeedToDagger())
                {
                    _muldag(toAdd, toMul);
                }
                else
                {
                    _mul(toAdd, toMul);
                }
            }
            _add(res, toAdd);
        }
        pStapleData[linkIndex] = res;
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelExpMultRealT_D(
    BYTE byFieldId,
    const deviceGauge * __restrict__ pMyDeviceData,
    Real a,
    deviceGauge *pU)
{
    intokernalInt4;
    const UINT uiDir = _DC_Dir;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    for (BYTE idir = 0; idir < uiDir; ++idir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, idir))
        {
            UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
            deviceGauge expP = _expreal(pMyDeviceData[linkIndex], a);
            _mul(expP, pU[linkIndex]);
            pU[linkIndex] = expP;
        }
        else
        {
            UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
            pU[linkIndex] = _makeId<deviceGauge>();
        }
    }
}

/**
* Trace (P^2)
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateKinematicEnergyT_D(
    BYTE byFieldId,
    const deviceGauge * __restrict__ pDeviceData,
    DOUBLE* results
)
{
    intokernalInt4;
    const UINT uiDir = _DC_Dir;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    Real resThisThread = F(0.0);
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, idir))
        {
            UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
            resThisThread += _retr(_dagmulC(pDeviceData[linkIndex], pDeviceData[linkIndex]));
        }
    }
    results[uiSiteIndex] = resThisThread;
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelFixBoundaryT_D(deviceGauge * pDeviceData, BYTE byFieldId)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const UINT uiDir = _DC_Dir;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        if (__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, idir))
        {
            SIndex idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
            pDeviceData[_deviceGetLinkIndex(uiSiteIndex, idir)] =
                ((CFieldBoundary<deviceGauge>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData
                [
                    __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idir
                ];
        }
    }
}


/**
 * iA = U.TA() 
 */
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformUToIAT_D(
    BYTE byFieldId,
    deviceGauge* pDeviceData)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            _ta(pDeviceData[uiLinkIndex]);
        }
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformUToIALogT_D(
    BYTE byFieldId,
    deviceGauge* pDeviceData)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = _strictlog(pDeviceData[uiLinkIndex]);
        }
    }
}

/**
 * U = exp(A)
 */
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformIAToUT_D(
    BYTE byFieldId,
    deviceGauge* pDeviceData)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = _expreal(pDeviceData[uiLinkIndex], F(1.0));
        }
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformIAToULogT_D(
    BYTE byFieldId,
    deviceGauge* pDeviceData)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = _strictexp(pDeviceData[uiLinkIndex]);
        }
    }
}

#pragma endregion

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLinkD<deviceGauge, matrixN>::MakeRandomGenerator()
{
    preparethread;
    _kernelInitialGeneratorT_D << <block, threads >> > (this->m_pDeviceData, this->m_byFieldId);
}

/**
* (1) calculate staples
* (2) calculate force(additive)
* (3) calculate energy
*/
template<typename deviceGauge, INT matrixN>
void CFieldGaugeLinkD<deviceGauge, matrixN>::CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStable, Real betaOverN) const
{
    if (NULL == pForce || this->GetFieldType() != pForce->GetFieldType())
    {
        appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: force field is not SU3");
        return;
    }
    if (NULL != pStable && this->GetFieldType() != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: stape field is not SU3");
        return;
    }

    CFieldGaugeLink<deviceGauge, matrixN>* pForceSU3 = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pForce);
    CFieldGaugeLink<deviceGauge, matrixN>* pStableSU3 = NULL == pStable ? NULL : dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pStable);

    preparethread;

    assert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache);

    _kernelStapleAtSiteCacheIndexT_D << <block, threads >> > (
        this->m_byFieldId,
        this->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        NULL == pStableSU3 ? NULL : pStableSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData,
        betaOverN);
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeLinkD<deviceGauge, matrixN>::CalculatePlaqutteEnergy(DOUBLE betaOverN) const
{
    assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);

    preparethread;
    _kernelPlaqutteEnergyCacheIndexT_D << <block, threads >> > (
        this->m_byFieldId,
        this->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeLinkD<deviceGauge, matrixN>::CalculateKinematicEnergy() const
{
    preparethread;
    _kernelCalculateKinematicEnergyT_D << <block, threads >> > (this->m_byFieldId, this->m_pDeviceData, _D_RealThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLinkD<deviceGauge, matrixN>::CalculateOnlyStaple(CFieldGauge* pStaple) const
{
    if (NULL == pStaple || this->GetFieldType() != pStaple->GetFieldType())
    {
        appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: stable field is not SU3");
        return;
    }
    CFieldGaugeLink<deviceGauge, matrixN>* pStapleSU3 = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pStaple);

    preparethread;
    _kernelCalculateOnlyStapleT_D << <block, threads >> > (
        this->m_byFieldId,
        this->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStapleSU3->m_pDeviceData);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLinkD<deviceGauge, matrixN>::ExpMult(Real a, CField* U) const
{
    if (NULL == U || this->GetFieldType() != U->GetFieldType())
    {
        appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: U field is not SU3");
        return;
    }

    CFieldGaugeLink<deviceGauge, matrixN>* pUField = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(U);

    preparethread;
    _kernelExpMultRealT_D << < block, threads >> > (this->m_byFieldId, this->m_pDeviceData, a, pUField->m_pDeviceData);
    
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLinkD<deviceGauge, matrixN>::FixBoundary()
{
    appDetailed(_T("CFieldGaugeLinkD<deviceGauge, matrixN>::FixBoundary()\n"));

    preparethread;
    _kernelFixBoundaryT_D << <block, threads >> > (this->m_pDeviceData, this->m_byFieldId);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLinkD<deviceGauge, matrixN>::TransformToIA()
{
    preparethread;
    _kernelTransformUToIAT_D << <block, threads >> > (this->m_byFieldId, this->m_pDeviceData);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLinkD<deviceGauge, matrixN>::TransformToU()
{
    preparethread;
    _kernelTransformIAToUT_D << <block, threads >> > (this->m_byFieldId, this->m_pDeviceData);
}

template<typename deviceGauge, INT matrixN>
CCString CFieldGaugeLinkD<deviceGauge, matrixN>::GetInfos(const CCString &tab) const
{
    CCString sRet = CFieldGaugeLink<deviceGauge, matrixN>::GetInfos(tab);
    SSmallInt4 boundary = appGetLattice()->m_pIndex->GetBoudanryCondition()->GetFieldBC(this->m_byFieldId);
    sRet = sRet + tab + appToString(boundary) + _T("\n");

    return sRet;
}

void CFieldGaugeU1D::TransformToIA()
{
    preparethread;
    if (0 == _HC_ALog)
    {
        _kernelTransformUToIAT_D << <block, threads >> > (this->m_byFieldId, this->m_pDeviceData);
    }
    else
    {
        _kernelTransformUToIALogT_D << <block, threads >> > (this->m_byFieldId, this->m_pDeviceData);
    }
}

void CFieldGaugeU1D::InitialWithByteCompressed(const CCString& sFileName)
{
    UINT uiSize = static_cast<UINT>(sizeof(Real) * m_uiLinkeCount);
    BYTE* byData = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
    Real* fRead = (Real*)byData;
    CLGComplex* readData = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiLinkeCount);
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        readData[i].x = F(0.0);
        readData[i].y = fRead[i];
    }
    checkCudaErrors(cudaMemcpy(this->m_pDeviceData, readData, sizeof(CLGComplex) * m_uiLinkeCount, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    free(readData);

    preparethread;
    _kernelTransformIAToULogT_D << <block, threads >> > (this->m_byFieldId, this->m_pDeviceData);
    checkCudaErrors(cudaDeviceSynchronize());

    free(byData);
}

CCString CFieldGaugeU1D::SaveToCompressedFile(const CCString& fileName) const
{
    CFieldGaugeU1* pPooledGauge = dynamic_cast<CFieldGaugeU1*>(GetCopy());

    preparethread;
    _kernelTransformUToIALogT_D << <block, threads >> > (m_byFieldId, pPooledGauge->m_pDeviceData);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    CLGComplex* toSave = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, pPooledGauge->m_pDeviceData, sizeof(CLGComplex) * m_uiLinkeCount, cudaMemcpyDeviceToHost));

    //This is a traceless anti-Hermitian now, so we only save part of them
    const UINT uiSize = static_cast<UINT>(sizeof(Real) * m_uiLinkeCount);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    Real* fToSave = (Real*)byToSave;
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        fToSave[i] = static_cast<Real>(toSave[i].y);
    }

    appGetFileSystem()->WriteAllBytes(fileName.c_str(), byToSave, uiSize);
    //pPooledGauge->DebugPrintMe();
    free(toSave);
    CCString MD5 = CLGMD5Hash(byToSave, uiSize);
    free(byToSave);
    appSafeDelete(pPooledGauge);
    return MD5;
}

void CFieldGaugeSU2D::TransformToIA()
{
    preparethread;
    if (0 == _HC_ALog)
    {
        _kernelTransformUToIAT_D << <block, threads >> > (this->m_byFieldId, this->m_pDeviceData);
    }
    else
    {
        _kernelTransformUToIALogT_D << <block, threads >> > (this->m_byFieldId, this->m_pDeviceData);
    }
}

void CFieldGaugeSU2D::InitialWithByteCompressed(const CCString& sFileName)
{
    UINT uiSize = static_cast<UINT>(sizeof(Real) * 3 * m_uiLinkeCount);
    BYTE* byData = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);

    deviceSU2* readData = (deviceSU2*)malloc(sizeof(deviceSU2) * m_uiLinkeCount);
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[3];
        memcpy(oneLink, byData + sizeof(Real) * 3 * i, sizeof(Real) * 3);

        readData[i].m_me[1] = _make_cuComplex(oneLink[0], oneLink[1]);
        readData[i].m_me[2] = _make_cuComplex(-oneLink[0], oneLink[1]);

        readData[i].m_me[0] = _make_cuComplex(F(0.0), oneLink[2]);
        readData[i].m_me[3] = _make_cuComplex(F(0.0), -oneLink[2]);
    }
    checkCudaErrors(cudaMemcpy(this->m_pDeviceData, readData, sizeof(deviceSU2) * m_uiLinkeCount, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    free(readData);

    //DebugPrintMe();

    preparethread;
    _kernelTransformIAToULogT_D << <block, threads >> > (this->m_byFieldId, this->m_pDeviceData);
    checkCudaErrors(cudaDeviceSynchronize());

    //DebugPrintMe();
    free(byData);
}

CCString CFieldGaugeSU2D::SaveToCompressedFile(const CCString& fileName) const
{
    CFieldGaugeSU2* pPooledGauge = dynamic_cast<CFieldGaugeSU2*>(GetCopy());

    preparethread;
    _kernelTransformUToIALogT_D << <block, threads >> > (pPooledGauge->m_byFieldId, pPooledGauge->m_pDeviceData);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    deviceSU2* toSave = (deviceSU2*)malloc(sizeof(deviceSU2) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, pPooledGauge->m_pDeviceData, sizeof(deviceSU2) * m_uiLinkeCount, cudaMemcpyDeviceToHost));

    //This is a traceless anti-Hermitian now, so we only save part of them
    const UINT uiSize = static_cast<UINT>(sizeof(Real) * m_uiLinkeCount * 3);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[3];
        oneLink[0] = static_cast<Real>(toSave[i].m_me[1].x);
        oneLink[1] = static_cast<Real>(toSave[i].m_me[1].y);
        oneLink[2] = static_cast<Real>(toSave[i].m_me[0].x);

        memcpy(byToSave + i * sizeof(Real) * 3, oneLink, sizeof(Real) * 3);
    }

    appGetFileSystem()->WriteAllBytes(fileName.c_str(), byToSave, uiSize);
    //pPooledGauge->DebugPrintMe();
    free(toSave);
    CCString MD5 = CLGMD5Hash(byToSave, uiSize);
    free(byToSave);
    appSafeDelete(pPooledGauge);
    return MD5;
}

__CLG_FORCETEMPLATE_CONSTRUCTOR(CFieldGaugeLinkD, U1D, CLGComplex, 1)
__CLG_FORCETEMPLATE_CONSTRUCTOR(CFieldGaugeLinkD, SU2D, deviceSU2, 2)
__CLG_FORCETEMPLATE_CONSTRUCTOR(CFieldGaugeLinkD, SU4D, deviceSU4, 4)
__CLG_FORCETEMPLATE_CONSTRUCTOR(CFieldGaugeLinkD, SU5D, deviceSU5, 5)
__CLG_FORCETEMPLATE_CONSTRUCTOR(CFieldGaugeLinkD, SU6D, deviceSU6, 6)
__CLG_FORCETEMPLATE_CONSTRUCTOR(CFieldGaugeLinkD, SU7D, deviceSU7, 7)
__CLG_FORCETEMPLATE_CONSTRUCTOR(CFieldGaugeLinkD, SU8D, deviceSU8, 8)

__CLGIMPLEMENT_CLASS(CFieldGaugeU1D)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU2D)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU4D)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU5D)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU6D)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU7D)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU8D)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================