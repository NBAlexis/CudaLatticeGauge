//=============================================================================
// FILENAME : CFieldGaugeSU3_12.cu
//
// DESCRIPTION:
//   Compact SU(3) gauge field using deviceSU3_12 (12 reals instead of 18).
//   Used as backup buffer (m_pUPrime) in force-gradient integrators.
//   Supports cross-type CopyTo with CFieldGaugeSU3.
//
// REVISION:
//  [05/10/2026 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "CFieldGaugeSU3_12.h"
#include "CFieldGaugeLink.h"
#include "Tools/Math/DeviceInlineTemplate.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldGaugeSU3_12)

//=============================================================================
// CUDA Kernels for SU3 <-> SU3_12 conversion
//=============================================================================

__global__ void _CLG_LAUNCH_BOUND
_kernelSU3ToSU3_12(deviceSU3_12* dest, const deviceSU3* src, UINT count)
{
    __simplekernel(
        dest[idx] = deviceSU3_12(src[idx]);
    )
}

__global__ void _CLG_LAUNCH_BOUND
_kernelSU3_12ToSU3(deviceSU3* dest, const deviceSU3_12* src, UINT count)
{
    __simplekernel(
        dest[idx] = src[idx].toSU3();
    )
}

//=============================================================================
// Static member functions for SU3 <-> SU3_12 conversion
//=============================================================================

void CFieldGaugeSU3_12::CopySU3ToSU3_12(deviceSU3_12* pDest, const deviceSU3* pSrc, UINT uiCount)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelSU3ToSU3_12, block, threads, pDest, pSrc, uiCount);
}

void CFieldGaugeSU3_12::CopySU3_12ToSU3(deviceSU3* pDest, const deviceSU3_12* pSrc, UINT uiCount)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelSU3_12ToSU3, block, threads, pDest, pSrc, uiCount);
}

//=============================================================================
// MSE kernel: ||SU3_12-expanded - SU3||^2 per link, reduce-summed
//=============================================================================

__global__ void _CLG_LAUNCH_BOUND
_kernelSU3_12MSE(const deviceSU3_12* pCompact, const deviceSU3* pRef, DOUBLE* result, UINT count)
{
    __simplekernel(
        deviceSU3 expanded = pCompact[idx].toSU3();
        expanded.Sub(pRef[idx]);
        result[idx] = static_cast<DOUBLE>(_lensq(expanded));
    )
}

DOUBLE CFieldGaugeSU3_12::SU3_12MSE(const deviceSU3_12* pCompact, const deviceSU3* pRef, UINT uiLinkCount)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelSU3_12MSE, block, threads, pCompact, pRef, _D_RealThreadBuffer, uiLinkCount);
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

//=============================================================================
// Kernels for Dot and GetLength on deviceSU3_12
//=============================================================================

__global__ void _CLG_LAUNCH_BOUND
_kernelSU3_12Dot(const deviceSU3_12* pMe, const deviceSU3_12* pOther, cuDoubleComplex* result, UINT count)
{
    __simplekernel(
        result[idx] = _cToDouble(pMe[idx].DaggerMulC(pOther[idx]).Tr());
    )
}

__global__ void _CLG_LAUNCH_BOUND
_kernelSU3_12LengthSq(const deviceSU3_12* pMe, DOUBLE* result, UINT count)
{
    __simplekernel(
        result[idx] = static_cast<DOUBLE>(pMe[idx].DaggerMulC(pMe[idx]).ReTr());
    )
}

//=============================================================================
// Kernel for InitialField (Identity)
//=============================================================================

__global__ void _CLG_LAUNCH_BOUND
_kernelSU3_12InitialIdentity(deviceSU3_12* pData, UINT count)
{
    __simplekernel(
        pData[idx].Id();
    )
}

//=============================================================================
// Kernel for Dagger on deviceSU3_12
//=============================================================================

__global__ void _CLG_LAUNCH_BOUND
_kernelSU3_12Dagger(deviceSU3_12* pData, UINT count)
{
    __simplekernel(
        pData[idx].Dagger();
    )
}

//=============================================================================
// Kernels for Mul and LeftMul on deviceSU3_12
//=============================================================================

__global__ void _CLG_LAUNCH_BOUND
_kernelSU3_12Mul(deviceSU3_12* pMe, const deviceSU3_12* pOther, UINT count, UBOOL bDaggerLeft, UBOOL bDaggerRight)
{
    __simplekernel(
        if (bDaggerLeft)
        {
            pMe[idx].Dagger();
        }
        if (bDaggerRight)
        {
            pMe[idx].MulDagger(pOther[idx]);
        }
        else
        {
            pMe[idx].Mul(pOther[idx]);
        }
    )
}

__global__ void _CLG_LAUNCH_BOUND
_kernelSU3_12LeftMul(deviceSU3_12* pMe, const deviceSU3_12* pOther, UINT count, UBOOL bDaggerLeft, UBOOL bDaggerRight)
{
    __simplekernel(
        if (bDaggerLeft)
        {
            if (bDaggerRight)
            {
                pMe[idx] = pOther[idx].DaggerC().MulDaggerC(pMe[idx]);
            }
            else
            {
                pMe[idx] = pOther[idx].DaggerMulC(pMe[idx]);
            }
        }
        else
        {
            if (bDaggerRight)
            {
                pMe[idx] = pOther[idx].MulDaggerC(pMe[idx]);
            }
            else
            {
                pMe[idx] = pOther[idx].MulC(pMe[idx]);
            }
        }
    )
}

//=============================================================================
// CFieldGaugeSU3_12 Implementation
//=============================================================================

CFieldGaugeSU3_12::CFieldGaugeSU3_12()
    : CFieldGauge()
    , m_pDeviceSU3_12Data(NULL)
{
    //Multi-GPU (Improve-1 I3): append halo link slots after the local links,
    //same Design-B layout as CFieldGaugeLink; m_uiLinkeCount stays the physics
    //volume. Single-GPU: _HC_HaloLinkCount() == 0 -> identical allocation.
    m_uiHaloLinkCount = _HC_HaloLinkCount();
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceSU3_12Data,
        (m_uiLinkeCount + m_uiHaloLinkCount) * sizeof(deviceSU3_12)));

    //Improve-1 (multi-GPU-improve1.md 3.1/3.2): bind the halo handle to this
    //exact extent. Dir compact matrices per site; site counts are lattice
    //sites. The field id is assigned after construction and synced through
    //SetFieldId.
    SHaloBufferInfo sInfo;
    sInfo.m_pDeviceData = reinterpret_cast<BYTE*>(m_pDeviceSU3_12Data);
    sInfo.m_uiCapacityBytes = (m_uiLinkeCount + m_uiHaloLinkCount) * sizeof(deviceSU3_12);
    sInfo.m_uiBytesPerSite = _HC_Dir * sizeof(deviceSU3_12);
    sInfo.m_uiLocalSiteCount = _HC_Volume;
    sInfo.m_uiHaloSiteCount = _HC_HaloSiteCount();
    sInfo.m_ullLayoutGeneration = appGetLayoutGeneration();
    sInfo.m_byFieldId = 0;
    sInfo.m_bHaloCapable = TRUE;
    m_HaloBuffer.Bind(sInfo);
}

CFieldGaugeSU3_12::~CFieldGaugeSU3_12()
{
    checkCudaErrors(__cudaFree(m_pDeviceSU3_12Data));
    m_pDeviceSU3_12Data = NULL;
}

void CFieldGaugeSU3_12::InitialField(EFieldInitialType eInitialType)
{
    if (EFIT_Identity == eInitialType)
    {
        preparethreadDir;
        _LAUNCH_KERNEL(_kernelSU3_12InitialIdentity, block, threads, m_pDeviceSU3_12Data, m_uiLinkeCount);
    }
    else if (EFIT_Zero == eInitialType)
    {
        checkCudaErrors(cudaMemset(m_pDeviceSU3_12Data, 0, m_uiLinkeCount * sizeof(deviceSU3_12)));
    }
    else
    {
        appCrucial(_T("CFieldGaugeSU3_12::InitialField: only EFIT_Identity and EFIT_Zero are supported\n"));
    }
    NotifyWritten();
}

void CFieldGaugeSU3_12::InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFileType)
{
    // Read SU3_12 format directly (12 reals per link)
    if (!CFileSystem::IsFileExist(sFileName))
    {
        appCrucial(_T("File not exist!!! %s \n"), sFileName.c_str());
        _FAIL_EXIT;
    }

    UINT uiSize = static_cast<UINT>(sizeof(Real) * 12 * m_uiLinkeCount);
#if _CLG_MULTI_GPU
    //On disk the file is the whole global lattice; read it fully on every rank
    //so root can scatter (non-root's copy is discarded by the scatter).
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        uiSize = static_cast<UINT>(sizeof(Real) * 12 * _HC_LinkCount * appGetComm()->Size());
    }
#endif
    BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
#if _CLG_MULTI_GPU
    //On disk the file is the whole global lattice in global-site order (written
    //by SaveToFile's gather). Each rank read the full file; scatter so every
    //rank keeps only its own sub-lattice in local-site order (mirrors the
    //CFieldGaugeLink<float/double> branches; P5-1.x). No-op on a lone rank.
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        const UINT uiPerRank = static_cast<UINT>(sizeof(Real) * 12 * _HC_LinkCount);
        if (uiSize != uiPerRank * appGetComm()->Size())
        {
            appCrucial(_T("Loading file size not match (MG): %s, %d, expecting global %d\n"),
                sFileName.c_str(), uiSize, static_cast<INT>(uiPerRank * appGetComm()->Size()));
            _FAIL_EXIT;
        }
        const UINT uiBytesPerSite = static_cast<UINT>(sizeof(Real) * 12 * _HC_Dir);
        BYTE* byLocal = (BYTE*)malloc(uiPerRank);
        appGetComm()->ScatterFieldFromRoot(data, uiBytesPerSite, byLocal);
        free(data);
        data = byLocal;
    }
#endif
    InitialWithByte(data);
    free(data);
    FixBoundary(EFB_Field);
    NotifyWritten();
}

void CFieldGaugeSU3_12::InitialWithByte(BYTE* byData)
{
    // Read SU3_12 format: 6 complex = 12 reals per link
    deviceSU3_12* pHost = (deviceSU3_12*)malloc(m_uiLinkeCount * sizeof(deviceSU3_12));
    Real* freadData = (Real*)byData;
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        for (UINT j = 0; j < 6; ++j)
        {
            pHost[i].m_me[j].x = freadData[12 * i + 2 * j];
            pHost[i].m_me[j].y = freadData[12 * i + 2 * j + 1];
        }
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceSU3_12Data, pHost, m_uiLinkeCount * sizeof(deviceSU3_12), cudaMemcpyHostToDevice));
    free(pHost);
    NotifyWritten();
}

void CFieldGaugeSU3_12::CopyBufferTo(CField* pTarget) const
{
    if (NULL == pTarget) return;

    if (EFT_GaugeSU3_12 == pTarget->GetFieldType())
    {
        CFieldGaugeSU3_12* pTargetField = dynamic_cast<CFieldGaugeSU3_12*>(pTarget);
        checkCudaErrors(cudaMemcpy(pTargetField->m_pDeviceSU3_12Data, m_pDeviceSU3_12Data,
            m_uiLinkeCount * sizeof(deviceSU3_12), cudaMemcpyDeviceToDevice));
    }
    else if (EFT_GaugeSU3 == pTarget->GetFieldType())
    {
        CFieldGaugeSU3* pTargetField = dynamic_cast<CFieldGaugeSU3*>(pTarget);
        CopySU3_12ToSU3(pTargetField->m_pDeviceData, m_pDeviceSU3_12Data, m_uiLinkeCount);
    }
    else
    {
        appCrucial(_T("CFieldGaugeSU3_12::CopyBufferTo: unsupported target field type\n"));
    }
    pTarget->NotifyWritten();
}

void CFieldGaugeSU3_12::CopyTo(CField* U) const
{
    if (NULL == U) return;
    CopyParamTo(U);
    CopyBufferTo(U);
}

UINT CFieldGaugeSU3_12::GetDeviceMemorySize() const
{
    return static_cast<UINT>(m_uiLinkeCount * sizeof(deviceSU3_12));
}

void CFieldGaugeSU3_12::DebugPrintMe() const
{
    // Expand SU3_12 -> SU3 on device, then copy to host and print
    deviceSU3* pDeviceSU3 = NULL;
    checkCudaErrors(__cudaMalloc((void**)&pDeviceSU3, m_uiLinkeCount * sizeof(deviceSU3)));
    CopySU3_12ToSU3(pDeviceSU3, m_pDeviceSU3_12Data, m_uiLinkeCount);

    deviceSU3* pHostSU3 = (deviceSU3*)malloc(m_uiLinkeCount * sizeof(deviceSU3));
    checkCudaErrors(cudaMemcpy(pHostSU3, pDeviceSU3, m_uiLinkeCount * sizeof(deviceSU3), cudaMemcpyDeviceToHost));

    for (UINT uiLink = 0; uiLink < m_uiLinkeCount; ++uiLink)
    {
        UINT uiSite = uiLink / _HC_Dir;
        UINT uiDir = uiLink % _HC_Dir;
        SSmallInt4 site = __hostSiteIndexToInt4(uiSite);
        appGeneral(_T(" --- %d(%d, %d, %d, %d)_%d ---\n %s\n"),
            uiLink,
            static_cast<INT>(site.x), static_cast<INT>(site.y),
            static_cast<INT>(site.z), static_cast<INT>(site.w),
            uiDir,
            appToString(pHostSU3[uiLink]).c_str());
    }

    free(pHostSU3);
    checkCudaErrors(__cudaFree(pDeviceSU3));
}

const void* CFieldGaugeSU3_12::GetData() const { return m_pDeviceSU3_12Data; }
void* CFieldGaugeSU3_12::GetData() { return m_pDeviceSU3_12Data; }

//=============================================================================
// BLAS operations
//=============================================================================

void CFieldGaugeSU3_12::AxpyPlus(const CField* x)
{
    appCrucial(_T("CFieldGaugeSU3_12: AxpyPlus not supported\n"));
}

void CFieldGaugeSU3_12::AxpyMinus(const CField* x)
{
    appCrucial(_T("CFieldGaugeSU3_12: AxpyMinus not supported\n"));
}

void CFieldGaugeSU3_12::Axpy(Real a, const CField* x)
{
    appCrucial(_T("CFieldGaugeSU3_12: Axpy not supported\n"));
}

void CFieldGaugeSU3_12::Axpy(const CLGComplex& a, const CField* x)
{
    appCrucial(_T("CFieldGaugeSU3_12: Axpy(complex) not supported\n"));
}

void CFieldGaugeSU3_12::Mul(const CField* other, UBOOL bDaggerLeft, UBOOL bDaggerRight)
{
    if (EFT_GaugeSU3_12 == other->GetFieldType())
    {
        const CFieldGaugeSU3_12* pOther = dynamic_cast<const CFieldGaugeSU3_12*>(other);
        preparethreadDir;
        _LAUNCH_KERNEL(_kernelSU3_12Mul, block, threads,
            m_pDeviceSU3_12Data, pOther->m_pDeviceSU3_12Data, m_uiLinkeCount, bDaggerLeft, bDaggerRight);
    }
    else if (EFT_GaugeSU3 == other->GetFieldType())
    {
        // Mixed: expand self to SU3, do SU3 Mul, compress back
        deviceSU3* pTempSU3 = NULL;
        checkCudaErrors(__cudaMalloc((void**)&pTempSU3, m_uiLinkeCount * sizeof(deviceSU3)));
        CopySU3_12ToSU3(pTempSU3, m_pDeviceSU3_12Data, m_uiLinkeCount);

        CFieldGaugeSU3* pTempField = new CFieldGaugeSU3();
        checkCudaErrors(cudaMemcpy(pTempField->m_pDeviceData, pTempSU3, m_uiLinkeCount * sizeof(deviceSU3), cudaMemcpyDeviceToDevice));
        pTempField->Mul(other, bDaggerLeft, bDaggerRight);
        CopySU3ToSU3_12(m_pDeviceSU3_12Data, pTempField->m_pDeviceData, m_uiLinkeCount);

        appSafeDelete(pTempField);
        checkCudaErrors(__cudaFree(pTempSU3));
    }
    else
    {
        appCrucial(_T("CFieldGaugeSU3_12::Mul: unsupported other type\n"));
    }
    NotifyWritten();
}

void CFieldGaugeSU3_12::LeftMul(const CField* other, UBOOL bDaggerLeft, UBOOL bDaggerRight)
{
    if (EFT_GaugeSU3_12 == other->GetFieldType())
    {
        const CFieldGaugeSU3_12* pOther = dynamic_cast<const CFieldGaugeSU3_12*>(other);
        preparethreadDir;
        _LAUNCH_KERNEL(_kernelSU3_12LeftMul, block, threads,
            m_pDeviceSU3_12Data, pOther->m_pDeviceSU3_12Data, m_uiLinkeCount, bDaggerLeft, bDaggerRight);
    }
    else if (EFT_GaugeSU3 == other->GetFieldType())
    {
        // Mixed: expand self to SU3, do SU3 LeftMul, compress back
        deviceSU3* pTempSU3 = NULL;
        checkCudaErrors(__cudaMalloc((void**)&pTempSU3, m_uiLinkeCount * sizeof(deviceSU3)));
        CopySU3_12ToSU3(pTempSU3, m_pDeviceSU3_12Data, m_uiLinkeCount);

        CFieldGaugeSU3* pTempField = new CFieldGaugeSU3();
        checkCudaErrors(cudaMemcpy(pTempField->m_pDeviceData, pTempSU3, m_uiLinkeCount * sizeof(deviceSU3), cudaMemcpyDeviceToDevice));
        pTempField->LeftMul(other, bDaggerLeft, bDaggerRight);
        CopySU3ToSU3_12(m_pDeviceSU3_12Data, pTempField->m_pDeviceData, m_uiLinkeCount);

        appSafeDelete(pTempField);
        checkCudaErrors(__cudaFree(pTempSU3));
    }
    else
    {
        appCrucial(_T("CFieldGaugeSU3_12::LeftMul: unsupported other type\n"));
    }
    NotifyWritten();
}

void CFieldGaugeSU3_12::Dagger()
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelSU3_12Dagger, block, threads, m_pDeviceSU3_12Data, m_uiLinkeCount);
    NotifyWritten();
}

void CFieldGaugeSU3_12::ScalarMultply(Real a)
{
    appCrucial(_T("CFieldGaugeSU3_12: ScalarMultply not supported\n"));
}

void CFieldGaugeSU3_12::ScalarMultply(const CLGComplex& a)
{
    appCrucial(_T("CFieldGaugeSU3_12: ScalarMultply(complex) not supported\n"));
}

BYTE* CFieldGaugeSU3_12::CopyDataOut(UINT& uiSize) const
{
    // Save in SU3_12 format: 6 complex = 12 reals per link
    deviceSU3_12* pHost = (deviceSU3_12*)malloc(m_uiLinkeCount * sizeof(deviceSU3_12));
    checkCudaErrors(cudaMemcpy(pHost, m_pDeviceSU3_12Data, m_uiLinkeCount * sizeof(deviceSU3_12), cudaMemcpyDeviceToHost));
    uiSize = static_cast<UINT>(sizeof(Real) * 12 * m_uiLinkeCount);
    BYTE* saveData = (BYTE*)malloc(uiSize);
    Real* fsaveData = (Real*)saveData;
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        for (UINT j = 0; j < 6; ++j)
        {
            fsaveData[12 * i + 2 * j]     = pHost[i].m_me[j].x;
            fsaveData[12 * i + 2 * j + 1] = pHost[i].m_me[j].y;
        }
    }
    free(pHost);
    return saveData;
}

BYTE* CFieldGaugeSU3_12::CopyDataOutFloat(UINT& uiSize) const
{
    deviceSU3_12* pHost = (deviceSU3_12*)malloc(m_uiLinkeCount * sizeof(deviceSU3_12));
    checkCudaErrors(cudaMemcpy(pHost, m_pDeviceSU3_12Data, m_uiLinkeCount * sizeof(deviceSU3_12), cudaMemcpyDeviceToHost));
    uiSize = static_cast<UINT>(sizeof(FLOAT) * 12 * m_uiLinkeCount);
    BYTE* saveData = (BYTE*)malloc(uiSize);
    FLOAT* fsaveData = (FLOAT*)saveData;
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        for (UINT j = 0; j < 6; ++j)
        {
            fsaveData[12 * i + 2 * j]     = static_cast<FLOAT>(pHost[i].m_me[j].x);
            fsaveData[12 * i + 2 * j + 1] = static_cast<FLOAT>(pHost[i].m_me[j].y);
        }
    }
    free(pHost);
    return saveData;
}

BYTE* CFieldGaugeSU3_12::CopyDataOutDouble(UINT& uiSize) const
{
    deviceSU3_12* pHost = (deviceSU3_12*)malloc(m_uiLinkeCount * sizeof(deviceSU3_12));
    checkCudaErrors(cudaMemcpy(pHost, m_pDeviceSU3_12Data, m_uiLinkeCount * sizeof(deviceSU3_12), cudaMemcpyDeviceToHost));
    uiSize = static_cast<UINT>(sizeof(DOUBLE) * 12 * m_uiLinkeCount);
    BYTE* saveData = (BYTE*)malloc(uiSize);
    DOUBLE* fsaveData = (DOUBLE*)saveData;
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        for (UINT j = 0; j < 6; ++j)
        {
            fsaveData[12 * i + 2 * j]     = static_cast<DOUBLE>(pHost[i].m_me[j].x);
            fsaveData[12 * i + 2 * j + 1] = static_cast<DOUBLE>(pHost[i].m_me[j].y);
        }
    }
    free(pHost);
    return saveData;
}

cuDoubleComplex CFieldGaugeSU3_12::Dot(const CField* other) const
{
    if (EFT_GaugeSU3_12 != other->GetFieldType())
    {
        appCrucial(_T("CFieldGaugeSU3_12::Dot: other must be CFieldGaugeSU3_12\n"));
        cuDoubleComplex res = {0.0, 0.0};
        return res;
    }
    const CFieldGaugeSU3_12* pOther = dynamic_cast<const CFieldGaugeSU3_12*>(other);
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelSU3_12Dot, block, threads, m_pDeviceSU3_12Data, pOther->m_pDeviceSU3_12Data, _D_ComplexThreadBuffer, m_uiLinkeCount);
    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
}

DOUBLE CFieldGaugeSU3_12::GetLength() const
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelSU3_12LengthSq, block, threads, m_pDeviceSU3_12Data, _D_RealThreadBuffer, m_uiLinkeCount);
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

UBOOL CFieldGaugeSU3_12::ApplyOperator(EFieldOperator op, INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* pGauge, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields, EOperatorCoefficientType uiCoeffType, Real fCoeffReal, Real fCoeffImg, void* otherParameter)
{
    appCrucial(_T("CFieldGaugeSU3_12: ApplyOperator not supported\n"));
    return FALSE;
}

//=============================================================================
// CFieldGauge pure virtual stubs
//=============================================================================

void CFieldGaugeSU3_12::CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const
{
    appCrucial(_T("CFieldGaugeSU3_12: CalculateForceAndStaple not supported\n"));
}

void CFieldGaugeSU3_12::CalculateOnlyStaple(CFieldGauge* pStaple) const
{
    appCrucial(_T("CFieldGaugeSU3_12: CalculateOnlyStaple not supported\n"));
}

void CFieldGaugeSU3_12::MakeRandomGenerator()
{
    appCrucial(_T("CFieldGaugeSU3_12: MakeRandomGenerator not supported\n"));
}

DOUBLE CFieldGaugeSU3_12::CalculatePlaqutteEnergy(DOUBLE betaOverN) const
{
    appCrucial(_T("CFieldGaugeSU3_12: CalculatePlaqutteEnergy not supported\n"));
    return 0.0;
}

DOUBLE CFieldGaugeSU3_12::CalculatePlaqutteEnergyOriginal(DOUBLE betaOverN) const
{
    appCrucial(_T("CFieldGaugeSU3_12: CalculatePlaqutteEnergyOriginal not supported\n"));
    return 0.0;
}

DOUBLE CFieldGaugeSU3_12::CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const
{
    appCrucial(_T("CFieldGaugeSU3_12: CalculatePlaqutteEnergyUseClover not supported\n"));
    return 0.0;
}

DOUBLE CFieldGaugeSU3_12::CalculatePlaqutteEnergyUsingStaple(DOUBLE betaOverN, const CFieldGauge* pStaple) const
{
    appCrucial(_T("CFieldGaugeSU3_12: CalculatePlaqutteEnergyUsingStaple not supported\n"));
    return 0.0;
}

DOUBLE CFieldGaugeSU3_12::CalculateKinematicEnergy() const
{
    appCrucial(_T("CFieldGaugeSU3_12: CalculateKinematicEnergy not supported\n"));
    return 0.0;
}

void CFieldGaugeSU3_12::ExpMult(Real a, CField* U) const
{
    appCrucial(_T("CFieldGaugeSU3_12: ExpMult not supported\n"));
}

void CFieldGaugeSU3_12::ElementNormalize()
{
    appCrucial(_T("CFieldGaugeSU3_12: ElementNormalize not supported\n"));
}

void CFieldGaugeSU3_12::SetOneDirectionUnity(BYTE byDir)
{
    appCrucial(_T("CFieldGaugeSU3_12: SetOneDirectionUnity not supported\n"));
}

void CFieldGaugeSU3_12::SetOneDirectionZero(BYTE byDir)
{
    appCrucial(_T("CFieldGaugeSU3_12: SetOneDirectionZero not supported\n"));
}

void CFieldGaugeSU3_12::TransformToIA()
{
    appCrucial(_T("CFieldGaugeSU3_12: TransformToIA not supported\n"));
}

void CFieldGaugeSU3_12::TA()
{
    appCrucial(_T("CFieldGaugeSU3_12: TA not supported\n"));
}

void CFieldGaugeSU3_12::TransformToU()
{
    appCrucial(_T("CFieldGaugeSU3_12: TransformToU not supported\n"));
}

void CFieldGaugeSU3_12::CalculateE_Using_U(CFieldGauge* pResoult) const
{
    appCrucial(_T("CFieldGaugeSU3_12: CalculateE_Using_U not supported\n"));
}

void CFieldGaugeSU3_12::CalculateNablaE_Using_U(CFieldGauge* pResoult, UBOOL bNaive) const
{
    appCrucial(_T("CFieldGaugeSU3_12: CalculateNablaE_Using_U not supported\n"));
}

void CFieldGaugeSU3_12::PolyakovOnSpatialSite(cuDoubleComplex* buffer, BYTE byDir) const
{
    appCrucial(_T("CFieldGaugeSU3_12: PolyakovOnSpatialSite not supported\n"));
}

__END_NAMESPACE
