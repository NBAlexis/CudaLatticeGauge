//=============================================================================
// FILENAME : CUpdator.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [02/17/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

void CUpdator::SaveConfiguration(UINT uiUpdateStep) const
{
    TCHAR buff1[256];
    TCHAR buff2[256];
    appGetTimeNow(buff1, 256);
    appGetTimeUtc(buff2, 256);
    //also save the infomations
    CCString sConf;
    sConf.Format(_T("%s_%d.txt"), m_sConfigurationPrefix.c_str(), m_iAcceptedConfigurationCount);
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
        sConf2.Format(_T("%s_%d.con"), m_sConfigurationPrefix.c_str(), m_iAcceptedConfigurationCount);
        appGetLattice()->m_pGaugeField[0]->SaveToFile(sConf2);
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