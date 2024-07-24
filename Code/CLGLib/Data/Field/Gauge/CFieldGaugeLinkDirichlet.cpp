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

void CFieldGaugeU1D::InitialWithByteCompressed(const CCString& sFileName)
{
    UINT uiSize = static_cast<UINT>(sizeof(Real) * this->m_uiLinkeCount);
    BYTE* byData = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
    Real* fRead = (Real*)byData;
    CLGComplex* readData = (CLGComplex*)malloc(sizeof(CLGComplex) * this->m_uiLinkeCount);
    for (UINT i = 0; i < this->m_uiLinkeCount; ++i)
    {
        readData[i].x = F(0.0);
        readData[i].y = fRead[i];
    }
    checkCudaErrors(cudaMemcpy(this->m_pDeviceData, readData, sizeof(CLGComplex) * this->m_uiLinkeCount, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    free(readData);

    CCommonKernelLink<CLGComplex>::StrictExp(this->m_pDeviceData, this->m_byFieldId);
    checkCudaErrors(cudaDeviceSynchronize());

    free(byData);
}

CCString CFieldGaugeU1D::SaveToCompressedFile(const CCString& fileName) const
{
    CFieldGaugeU1* pPooledGauge = dynamic_cast<CFieldGaugeU1*>(GetCopy());

    CCommonKernelLink<CLGComplex>::StrictLog(pPooledGauge->m_pDeviceData, this->m_byFieldId);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    CLGComplex* toSave = (CLGComplex*)malloc(sizeof(CLGComplex) * this->m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, pPooledGauge->m_pDeviceData, sizeof(CLGComplex) * this->m_uiLinkeCount, cudaMemcpyDeviceToHost));

    //This is a traceless anti-Hermitian now, so we only save part of them
    const UINT uiSize = static_cast<UINT>(sizeof(Real) * this->m_uiLinkeCount);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    Real* fToSave = (Real*)byToSave;
    for (UINT i = 0; i < this->m_uiLinkeCount; ++i)
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

void CFieldGaugeSU2D::InitialWithByteCompressed(const CCString& sFileName)
{
    UINT uiSize = static_cast<UINT>(sizeof(Real) * 3 * this->m_uiLinkeCount);
    BYTE* byData = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);

    deviceSU2* readData = (deviceSU2*)malloc(sizeof(deviceSU2) * this->m_uiLinkeCount);
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[3];
        memcpy(oneLink, byData + sizeof(Real) * 3 * i, sizeof(Real) * 3);

        readData[i].m_me[1] = _make_cuComplex(oneLink[0], oneLink[1]);
        readData[i].m_me[2] = _make_cuComplex(-oneLink[0], oneLink[1]);

        readData[i].m_me[0] = _make_cuComplex(F(0.0), oneLink[2]);
        readData[i].m_me[3] = _make_cuComplex(F(0.0), -oneLink[2]);
    }
    checkCudaErrors(cudaMemcpy(this->m_pDeviceData, readData, sizeof(deviceSU2) * this->m_uiLinkeCount, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    free(readData);

    //DebugPrintMe();

    CCommonKernelLink<deviceSU2>::StrictExp(this->m_pDeviceData, this->m_byFieldId);
    checkCudaErrors(cudaDeviceSynchronize());

    //DebugPrintMe();
    free(byData);
}

CCString CFieldGaugeSU2D::SaveToCompressedFile(const CCString& fileName) const
{
    CFieldGaugeSU2* pPooledGauge = dynamic_cast<CFieldGaugeSU2*>(GetCopy());

    CCommonKernelLink<deviceSU2>::StrictLog(pPooledGauge->m_pDeviceData, this->m_byFieldId);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    deviceSU2* toSave = (deviceSU2*)malloc(sizeof(deviceSU2) * this->m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, pPooledGauge->m_pDeviceData, sizeof(deviceSU2) * this->m_uiLinkeCount, cudaMemcpyDeviceToHost));

    //This is a traceless anti-Hermitian now, so we only save part of them
    const UINT uiSize = static_cast<UINT>(sizeof(Real) * this->m_uiLinkeCount * 3);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < this->m_uiLinkeCount; ++i)
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