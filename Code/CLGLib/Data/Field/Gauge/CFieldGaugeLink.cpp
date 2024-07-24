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
//  [12/4/2018 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CFieldGaugeLink.h"

__BEGIN_NAMESPACE

#pragma region SU2 functions

void CFieldGaugeU1::InitialWithByteCompressed(const CCString& sFileName)
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

void CFieldGaugeSU2::InitialWithByteCompressed(const CCString& sFileName)
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


#pragma endregion

__CLGIMPLEMENT_CLASS(CFieldGaugeU1)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU2)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU4)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU5)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU6)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU7)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU8)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================