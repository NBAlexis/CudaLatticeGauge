//=============================================================================
// FILENAME : CFieldGaugeLink<deviceGauge, matrixN>.cpp
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
//  [mm/dd/yy]
//  [12/4/2018 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CFieldGaugeLink.h"
#include "CFieldGaugeSU3_12.h"

__BEGIN_NAMESPACE

#pragma region SU2 functions

void CFieldGaugeU1::InitialWithByteCompressed(const CCString& sFileName)
{
    //P5-1.3: the file is the whole GLOBAL lattice in global-site order
    //(SaveToCompressedFile gathers to rank 0); read it fully, then scatter to
    //the per-rank sub-lattice on multi-GPU. Previously each rank read by its
    //LOCAL link count -> misaligned/incomplete data, silent wrong.
    const UINT uiBytesPerSite = static_cast<UINT>(sizeof(Real) * _HC_Dir);
    const UINT uiPerRankBytes = static_cast<UINT>(sizeof(Real) * m_uiLinkeCount);
    UINT uiSize = 0;
    BYTE* byData = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
    if (NULL == byData)
    {
        appCrucial(_T("File not found: %s\n"), sFileName.c_str());
        _FAIL_EXIT;
    }
#if _CLG_MULTI_GPU
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        if (uiSize != uiPerRankBytes * appGetComm()->Size())
        {
            appCrucial(_T("File size not correct (MG): expecting global %d, found: %d\n"),
                static_cast<UINT>(uiPerRankBytes * appGetComm()->Size()), uiSize);
            _FAIL_EXIT;
        }
        BYTE* byLocal = (BYTE*)malloc(uiPerRankBytes);
        appGetComm()->ScatterFieldFromRoot(byData, uiBytesPerSite, byLocal);
        free(byData);
        byData = byLocal;
    }
    else
#endif
    if (uiSize != uiPerRankBytes)
    {
        appCrucial(_T("File size not correct: expecting: %d, found: %d\n"), uiPerRankBytes, uiSize);
        _FAIL_EXIT;
    }
    Real* fRead = (Real*)byData;
    CLGComplex* readData = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiLinkeCount);
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        readData[i].x = F(0.0);
        readData[i].y = fRead[i];
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(CLGComplex) * m_uiLinkeCount, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    free(readData);

    CCommonKernelLink<CLGComplex>::StrictExp(m_pDeviceData, m_byFieldId);
    checkCudaErrors(cudaDeviceSynchronize());

    free(byData);
}

CCString CFieldGaugeU1::SaveToCompressedFile(const CCString& fileName) const
{
    CFieldGaugeU1* pPooledGauge = dynamic_cast<CFieldGaugeU1*>(GetCopy());

    CCommonKernelLink<CLGComplex>::StrictLog(pPooledGauge->m_pDeviceData, m_byFieldId);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    CLGComplex* toSave = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, pPooledGauge->m_pDeviceData, sizeof(CLGComplex) * m_uiLinkeCount, cudaMemcpyDeviceToHost));

    //This is a traceless anti-Hermitian now, so we only save part of them
    const UINT uiPerRankBytes = static_cast<UINT>(sizeof(Real) * m_uiLinkeCount);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiPerRankBytes));
    Real* fToSave = (Real*)byToSave;
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        fToSave[i] = static_cast<Real>(toSave[i].y);
    }
    free(toSave);

    //P5-1.3: gather the compressed data to rank 0 in global-site order so the
    //file on disk is identical to the single-GPU layout (was appCrucial-rejected
    //under multi-GPU). Non-root ranks write nothing.
    UINT uiSize = uiPerRankBytes;
#if _CLG_MULTI_GPU
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        const UINT uiBytesPerSite = static_cast<UINT>(sizeof(Real) * _HC_Dir);
        UINT uiGlobalSize = 0;
        BYTE* byGlobal = appGetComm()->GatherFieldToRoot(byToSave, uiBytesPerSite, uiGlobalSize);
        free(byToSave);
        if (!appGetComm()->IsRoot())
        {
            appSafeDelete(pPooledGauge);
            return _T("");
        }
        byToSave = byGlobal;
        uiSize = uiGlobalSize;
    }
#endif

    appGetFileSystem()->WriteAllBytes(fileName.c_str(), byToSave, uiSize);
    //pPooledGauge->DebugPrintMe();
    CCString MD5 = CLGMD5Hash(byToSave, uiSize);
    free(byToSave);
    appSafeDelete(pPooledGauge);
    return MD5;
}

void CFieldGaugeSU2::InitialWithByteCompressed(const CCString& sFileName)
{
    //P5-1.3: the file is the whole GLOBAL lattice in global-site order
    //(SaveToCompressedFile gathers to rank 0); read it fully, then scatter to
    //the per-rank sub-lattice on multi-GPU. Previously each rank read by its
    //LOCAL link count -> misaligned/incomplete data, silent wrong.
    const UINT uiBytesPerSite = static_cast<UINT>(sizeof(Real) * 3 * _HC_Dir);
    const UINT uiPerRankBytes = static_cast<UINT>(sizeof(Real) * 3 * m_uiLinkeCount);
    UINT uiSize = 0;
    BYTE* byData = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
    if (NULL == byData)
    {
        appCrucial(_T("File not found: %s\n"), sFileName.c_str());
        _FAIL_EXIT;
    }
#if _CLG_MULTI_GPU
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        if (uiSize != uiPerRankBytes * appGetComm()->Size())
        {
            appCrucial(_T("File size not correct (MG): expecting global %d, found: %d\n"),
                static_cast<UINT>(uiPerRankBytes * appGetComm()->Size()), uiSize);
            _FAIL_EXIT;
        }
        BYTE* byLocal = (BYTE*)malloc(uiPerRankBytes);
        appGetComm()->ScatterFieldFromRoot(byData, uiBytesPerSite, byLocal);
        free(byData);
        byData = byLocal;
    }
    else
#endif
    if (uiSize != uiPerRankBytes)
    {
        appCrucial(_T("File size not correct: expecting: %d, found: %d\n"), uiPerRankBytes, uiSize);
        _FAIL_EXIT;
    }

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
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(deviceSU2) * m_uiLinkeCount, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    free(readData);

    //DebugPrintMe();

    CCommonKernelLink<deviceSU2>::StrictExp(m_pDeviceData, m_byFieldId);
    checkCudaErrors(cudaDeviceSynchronize());
    //DebugPrintMe();
    free(byData);
}

CCString CFieldGaugeSU2::SaveToCompressedFile(const CCString& fileName) const
{
    CFieldGaugeSU2* pPooledGauge = dynamic_cast<CFieldGaugeSU2*>(GetCopy());
    CCommonKernelLink<deviceSU2>::StrictLog(pPooledGauge->m_pDeviceData, m_byFieldId);

    deviceSU2* toSave = (deviceSU2*)malloc(sizeof(deviceSU2) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, pPooledGauge->m_pDeviceData, sizeof(deviceSU2) * m_uiLinkeCount, cudaMemcpyDeviceToHost));

    //This is a traceless anti-Hermitian now, so we only save part of them
    const UINT uiPerRankBytes = static_cast<UINT>(sizeof(Real) * m_uiLinkeCount * 3);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiPerRankBytes));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[3];
        oneLink[0] = static_cast<Real>(toSave[i].m_me[1].x);
        oneLink[1] = static_cast<Real>(toSave[i].m_me[1].y);
        oneLink[2] = static_cast<Real>(toSave[i].m_me[0].x);

        memcpy(byToSave + i * sizeof(Real) * 3, oneLink, sizeof(Real) * 3);
    }
    free(toSave);

    //P5-1.3: gather the compressed data to rank 0 in global-site order so the
    //file on disk is identical to the single-GPU layout (was appCrucial-rejected
    //under multi-GPU). Non-root ranks write nothing.
    UINT uiSize = uiPerRankBytes;
#if _CLG_MULTI_GPU
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        const UINT uiBytesPerSite = static_cast<UINT>(sizeof(Real) * 3 * _HC_Dir);
        UINT uiGlobalSize = 0;
        BYTE* byGlobal = appGetComm()->GatherFieldToRoot(byToSave, uiBytesPerSite, uiGlobalSize);
        free(byToSave);
        if (!appGetComm()->IsRoot())
        {
            appSafeDelete(pPooledGauge);
            return _T("");
        }
        byToSave = byGlobal;
        uiSize = uiGlobalSize;
    }
#endif

    appGetFileSystem()->WriteAllBytes(fileName.c_str(), byToSave, uiSize);
    //pPooledGauge->DebugPrintMe();
    CCString MD5 = CLGMD5Hash(byToSave, uiSize);
    free(byToSave);
    appSafeDelete(pPooledGauge);
    return MD5;
}

void CFieldGaugeSU3::CopyBufferTo(CField* pTarget) const
{
    if (NULL == pTarget) return;

    if (EFT_GaugeSU3 == pTarget->GetFieldType())
    {
        CFieldGaugeSU3* pTargetField = dynamic_cast<CFieldGaugeSU3*>(pTarget);
        checkCudaErrors(cudaMemcpy(pTargetField->m_pDeviceData, m_pDeviceData,
            m_uiLinkeCount * sizeof(deviceSU3), cudaMemcpyDeviceToDevice));
    }
    else if (EFT_GaugeSU3_12 == pTarget->GetFieldType())
    {
        CFieldGaugeSU3_12* pTargetField = dynamic_cast<CFieldGaugeSU3_12*>(pTarget);
        CFieldGaugeSU3_12::CopySU3ToSU3_12(pTargetField->m_pDeviceSU3_12Data, m_pDeviceData, m_uiLinkeCount);
    }
    else
    {
        appCrucial(_T("CFieldGaugeSU3::CopyBufferTo: unsupported target field type\n"));
    }
    pTarget->NotifyWritten();
}

void CFieldGaugeSU3::InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFileType)
{
    if (EFFT_CLGBinSU3_12 == eFileType)
    {
        if (!CFileSystem::IsFileExist(sFileName))
        {
            appCrucial(_T("File not exist!!! %s \n"), sFileName.c_str());
            _FAIL_EXIT;
        }
        CFieldGaugeSU3_12* pTemp = new CFieldGaugeSU3_12();
        pTemp->InitialFieldWithFile(sFileName, EFFT_CLGBin);
        pTemp->CopyBufferTo(this);
        appSafeDelete(pTemp);
        FixBoundary(EFB_Field);
    }
    else
    {
        CFieldGaugeLink<deviceSU3, 3>::InitialFieldWithFile(sFileName, eFileType);
    }
}

CCString CFieldGaugeSU3::SaveToFile(const CCString& fileName, EFieldFileType eType) const
{
    if (EFFT_CLGBinSU3_12 == eType)
    {
        CFieldGaugeSU3_12* pTemp = new CFieldGaugeSU3_12();
        CopyBufferTo(pTemp);
        CCString result = pTemp->SaveToFile(fileName, EFFT_CLGBin);
        appSafeDelete(pTemp);
        return result;
    }
    return CField::SaveToFile(fileName, eType);
}

void CFieldGaugeSU3::InitialWithByteCompressed(const CCString& sFileName)
{
    //P5-1.3: the file is the whole GLOBAL lattice in global-site order
    //(SaveToCompressedFile gathers to rank 0); read it fully, then scatter to
    //the per-rank sub-lattice on multi-GPU. Previously each rank read by its
    //LOCAL link count -> misaligned/incomplete data, silent wrong.
    const UINT uiBytesPerSite = static_cast<UINT>(sizeof(Real) * 9 * _HC_Dir);
    const UINT uiPerRankBytes = static_cast<UINT>(sizeof(Real) * 9 * m_uiLinkeCount);
    UINT uiSize = 0;
    BYTE* byData = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
    if (NULL == byData)
    {
        appCrucial(_T("File not found: %s\n"), sFileName.c_str());
        _FAIL_EXIT;
    }
#if _CLG_MULTI_GPU
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        if (uiSize != uiPerRankBytes * appGetComm()->Size())
        {
            appCrucial(_T("File size not correct (MG): expecting global %d, found: %d\n"),
                static_cast<UINT>(uiPerRankBytes * appGetComm()->Size()), uiSize);
            _FAIL_EXIT;
        }
        BYTE* byLocal = (BYTE*)malloc(uiPerRankBytes);
        appGetComm()->ScatterFieldFromRoot(byData, uiBytesPerSite, byLocal);
        free(byData);
        byData = byLocal;
    }
    else
#endif
    if (uiSize != uiPerRankBytes)
    {
        appCrucial(_T("File size not correct: expecting: %d, found: %d\n"), uiPerRankBytes, uiSize);
        _FAIL_EXIT;
    }

    deviceSU3* readData = (deviceSU3*)malloc(sizeof(deviceSU3) * m_uiLinkeCount);
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[9];
        memcpy(oneLink, byData + sizeof(Real) * 9 * i, sizeof(Real) * 9);

        readData[i].m_me[1] = _make_cuComplex(oneLink[0], oneLink[1]);
        readData[i].m_me[3] = _make_cuComplex(-oneLink[0], oneLink[1]);

        readData[i].m_me[2] = _make_cuComplex(oneLink[2], oneLink[3]);
        readData[i].m_me[6] = _make_cuComplex(-oneLink[2], oneLink[3]);

        readData[i].m_me[5] = _make_cuComplex(oneLink[4], oneLink[5]);
        readData[i].m_me[7] = _make_cuComplex(-oneLink[4], oneLink[5]);

        readData[i].m_me[0] = _make_cuComplex(F(0.0), oneLink[6]);
        readData[i].m_me[4] = _make_cuComplex(F(0.0), oneLink[7]);
        readData[i].m_me[8] = _make_cuComplex(F(0.0), oneLink[8]);

#if _CLG_PADDING
        //Only zero the padding slots when the matrix actually carries them;
        //without padding m_me has exactly 9 elements and writing [9..15] is a
        //heap overflow (fatal on the last link of large lattices).
        for (UINT j = 9; j < 16; ++j)
        {
            readData[i].m_me[j] = _zeroc;
        }
#endif
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(deviceSU3) * m_uiLinkeCount, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    free(readData);

    //DebugPrintMe();

    //preparethread;
    CCommonKernelLink<deviceSU3>::StrictExp(m_pDeviceData, m_byFieldId);
    checkCudaErrors(cudaDeviceSynchronize());

    //DebugPrintMe();
    free(byData);
    NotifyWritten();
}

CCString CFieldGaugeSU3::SaveToCompressedFile(const CCString& fileName) const
{
    CFieldGaugeSU3* pPooledGauge = dynamic_cast<CFieldGaugeSU3*>(GetCopy());
    //pPooledGauge->DebugPrintMe();

    //preparethread;
    CCommonKernelLink<deviceSU3>::StrictLog(pPooledGauge->m_pDeviceData, m_byFieldId);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    deviceSU3* toSave = (deviceSU3*)malloc(sizeof(deviceSU3) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, pPooledGauge->m_pDeviceData, sizeof(deviceSU3) * m_uiLinkeCount, cudaMemcpyDeviceToHost));

    //This is a traceless anti-Hermitian now, so we only save part of them
    const UINT uiPerRankBytes = static_cast<UINT>(sizeof(Real) * m_uiLinkeCount * 9);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiPerRankBytes));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[9];
        oneLink[0] = static_cast<Real>(toSave[i].m_me[1].x);
        oneLink[1] = static_cast<Real>(toSave[i].m_me[1].y);
        oneLink[2] = static_cast<Real>(toSave[i].m_me[2].x);
        oneLink[3] = static_cast<Real>(toSave[i].m_me[2].y);
        oneLink[4] = static_cast<Real>(toSave[i].m_me[5].x);
        oneLink[5] = static_cast<Real>(toSave[i].m_me[5].y);
        oneLink[6] = static_cast<Real>(toSave[i].m_me[0].y);
        oneLink[7] = static_cast<Real>(toSave[i].m_me[4].y);
        //The element is in fact can be NOT traceless!!!!, the trace can be 2 Pi or -2 Pi !!!
        oneLink[8] = static_cast<Real>(toSave[i].m_me[8].y);

        memcpy(byToSave + i * sizeof(Real) * 9, oneLink, sizeof(Real) * 9);
    }
    free(toSave);

    //P5-1.3: gather the compressed data to rank 0 in global-site order so the
    //file on disk is identical to the single-GPU layout (was appCrucial-rejected
    //under multi-GPU). Non-root ranks write nothing.
    UINT uiSize = uiPerRankBytes;
#if _CLG_MULTI_GPU
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        const UINT uiBytesPerSite = static_cast<UINT>(sizeof(Real) * 9 * _HC_Dir);
        UINT uiGlobalSize = 0;
        BYTE* byGlobal = appGetComm()->GatherFieldToRoot(byToSave, uiBytesPerSite, uiGlobalSize);
        free(byToSave);
        if (!appGetComm()->IsRoot())
        {
            appSafeDelete(pPooledGauge);
            return _T("");
        }
        byToSave = byGlobal;
        uiSize = uiGlobalSize;
    }
#endif

    appGetFileSystem()->WriteAllBytes(fileName.c_str(), byToSave, uiSize);
    //pPooledGauge->DebugPrintMe();
    CCString MD5 = CLGMD5Hash(byToSave, uiSize);
    free(byToSave);
    appSafeDelete(pPooledGauge);
    return MD5;
}

#pragma endregion

__CLGIMPLEMENT_CLASS(CFieldGaugeU1)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU2)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU3)

#if _CLG_SU4_GAUGE
__CLGIMPLEMENT_CLASS(CFieldGaugeSU4)
#endif
#if _CLG_SU5_GAUGE
__CLGIMPLEMENT_CLASS(CFieldGaugeSU5)
#endif
#if _CLG_SU6_GAUGE
__CLGIMPLEMENT_CLASS(CFieldGaugeSU6)
#endif
#if _CLG_SU7_GAUGE
__CLGIMPLEMENT_CLASS(CFieldGaugeSU7)
#endif
#if _CLG_SU8_GAUGE
__CLGIMPLEMENT_CLASS(CFieldGaugeSU8)
#endif
#if _CLG_Z2_GAUGE
__CLGIMPLEMENT_CLASS(CFieldGaugeZ2)
#endif
#if _CLG_Z3_GAUGE
__CLGIMPLEMENT_CLASS(CFieldGaugeZ3)
#endif
#if _CLG_Z4_GAUGE
__CLGIMPLEMENT_CLASS(CFieldGaugeZ4)
#endif
#if _CLG_Z5_GAUGE
__CLGIMPLEMENT_CLASS(CFieldGaugeZ5)
#endif
#if _CLG_Z6_GAUGE
__CLGIMPLEMENT_CLASS(CFieldGaugeZ6)
#endif
#if _CLG_D3_GAUGE
__CLGIMPLEMENT_CLASS(CFieldGaugeD3)
#endif
#if _CLG_D4_GAUGE
__CLGIMPLEMENT_CLASS(CFieldGaugeD4)
#endif
#if _CLG_D8_GAUGE
__CLGIMPLEMENT_CLASS(CFieldGaugeD8)
#endif

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================