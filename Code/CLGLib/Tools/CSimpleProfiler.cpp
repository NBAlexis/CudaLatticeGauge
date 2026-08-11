//=============================================================================
// FILENAME : CSimpleProfiler.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [mm/dd/yy]
//  [01/07/2025 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

CSimpleProfiler GSimpleProfiler;

CSimpleProfileRecord::CSimpleProfileRecord(const TCHAR* name)
    : m_sName(name)
    , m_uiStart(0)
{
    checkCudaErrors(cudaDeviceSynchronize());
    m_uiStart = std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
}

CSimpleProfileRecord::~CSimpleProfileRecord()
{
    if (!_HC_Profiler)
    {
        return;
    }
    checkCudaErrors(cudaDeviceSynchronize());
    const LONGLONG uiEclipse = std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count() - m_uiStart;
    GSimpleProfiler.AddRecord(m_sName, uiEclipse);
}

void CSimpleProfiler::Dump() const
{
    if (!_HC_Profiler)
    {
        return;
    }

    TArray<CCString> keys = m_NameHash.GetAllKeys();
    if (0 == keys.Num())
    {
        appGeneral(_T("It seems profiler is not turned on!\n"));
        return;
    }
    TArray<CCString> output;
    TArray<DOUBLE> totalduraction;
    TCHAR temp[4096];
    for (INT i = 0; i < keys.Num(); ++i)
    {
        INT idx = m_NameHash.GetAt(keys[i]);
        DOUBLE duration = static_cast<DOUBLE>(m_lstRecoreds[idx].m_uiTimeNanoseconds) * 0.000000001;
        totalduraction.AddItem(duration);

        appSprintf(temp, 4095, _T("%s:\nCalled %lld\tTotal %.6g (s)\tAverage %.6g (ms)\n"),
            keys[i].c_str(),
            m_lstRecoreds[idx].m_uiCalled,
            duration,
            static_cast<DOUBLE>(m_lstRecoreds[idx].m_uiTimeNanoseconds) * 0.000001 / m_lstRecoreds[idx].m_uiCalled
            );
        temp[4095] = 0;
        output.AddItem(CCString(temp));
    }

    for (INT i = 0; i < totalduraction.Num(); ++i)
    {
        for (INT j = i + 1; j < totalduraction.Num(); ++j)
        {
            if (totalduraction[i] < totalduraction[j])
            {
                DOUBLE tempduration = totalduraction[i];
                totalduraction[i] = totalduraction[j];
                totalduraction[j] = tempduration;

                CCString tempoutput = output[i];
                output[i] = output[j];
                output[j] = tempoutput;
            }
        }
    }

    appPushLogDate(FALSE);
    appGeneral(_T("\n=========\n\n"));
    for (INT i = 0; i < output.Num(); ++i)
    {
        appGeneral(output[i]);
    }
    appGeneral(_T("\n=========\n\n"));
    appPopLogDate();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
