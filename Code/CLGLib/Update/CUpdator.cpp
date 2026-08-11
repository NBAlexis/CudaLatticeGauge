//=============================================================================
// FILENAME : CUpdator.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [mm/dd/yy]
//  [02/17/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

void CUpdator::SaveConfiguration(UINT uiUpdateStep) const
{
    _RECORD(CUpdator::SaveConfiguration);
    TCHAR buff1[256];
    TCHAR buff2[256];
    appGetTimeNow(buff1, 256);
    appGetTimeUtc(buff2, 256);
    //also save the infomations
    CCString sConf;
    sConf.Format(_T("%s_%d.txt"), m_sConfigurationPrefix.c_str(), m_iAcceptedConfigurationCount + m_iSaveIndexStart);
    CCString sInfo;
    sInfo.Format(_T("TimeStamp : %d\nTime : %s\nTimeUTC : %s\n"),
        appGetTimeStamp(),
        buff1,
        buff2);
    sInfo = sInfo + appGetLattice()->GetInfos(_T(""));
    appGetFileSystem()->WriteAllText(sConf, sInfo);

    if (appGetLattice()->m_pGaugeField.Num() > 0)
    {
        CCString sConf2;
        sConf2.Format(_T("%s_%d.con"), m_sConfigurationPrefix.c_str(), m_iAcceptedConfigurationCount + m_iSaveIndexStart);
        appGetLattice()->m_pGaugeField[0]->SaveToFile(sConf2, m_eSaveFieldType);
    }

    // save the dynamic tensor2 fields as companion files:
    //   <prefix>_<N>_t<fieldId>.con
    // This keeps lattice fields such as the PSU(3) Z3 2-form boundary field B
    // consistent with the gauge configuration. Non-dynamic tensor2 fields are
    // auxiliary and are not part of the configuration.
    const UINT uiN = m_iAcceptedConfigurationCount + m_iSaveIndexStart;
    for (INT i = 0; i < appGetLattice()->m_pTensor2Field.Num(); ++i)
    {
        CFieldTensor2* pTensor2 = appGetLattice()->m_pTensor2Field[i];
        if (NULL == pTensor2 || !pTensor2->IsDynamic())
        {
            continue;
        }
        CCString sTensor2;
        sTensor2.Format(_T("%s_%u_t%d.con"), m_sConfigurationPrefix.c_str(), uiN, pTensor2->m_byFieldId);
        pTensor2->SaveToFile(sTensor2, m_eSaveFieldType);
    }
}

void CUpdator::LoadTensor2Configuration(const CCString& sPrefix, UINT uiN, EFieldFileType eLoadType) const
{
    _RECORD(CUpdator::LoadTensor2Configuration);
    for (INT i = 0; i < appGetLattice()->m_pTensor2Field.Num(); ++i)
    {
        CFieldTensor2* pTensor2 = appGetLattice()->m_pTensor2Field[i];
        if (NULL == pTensor2 || !pTensor2->IsDynamic())
        {
            continue;
        }
        CCString sTensor2;
        sTensor2.Format(_T("%s_%u_t%d.con"), sPrefix.c_str(), uiN, pTensor2->m_byFieldId);
        if (CFileSystem::IsFileExist(sTensor2))
        {
            pTensor2->InitialFieldWithFile(sTensor2, eLoadType);
            appGeneral(_T("CUpdator::LoadTensor2Configuration: loaded %s\n"), sTensor2.c_str());
        }
    }
}

void CUpdator::UpdateUntileAccept(UINT iSteps, UBOOL bMeasure)
{
    m_iAcceptedConfigurationCount = 0;
    while (m_iAcceptedConfigurationCount < iSteps)
    {
        Update(1, bMeasure);
    }
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================