//=============================================================================
// FILENAME : Timer.h
// 
// DESCRIPTION:
// This is the timer for tester use
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================

#ifndef _TIMER_H_
#define _TIMER_H_

__BEGIN_NAMESPACE

class CLGAPI CTimer
{
public:
    CTimer(const CCString& sId = _T(""), const bool bReport = false)
        : m_bStarted(FALSE)
        , m_uiStart(0)
        , m_fElapsed(0)
        , m_dwCounter(0)
        , m_sId(sId)
        , m_bReportOnExit(bReport)
    {}

    ~CTimer()
    {
        if (m_bReportOnExit) Report();
    }

public:

    void Start()
    {
        m_bStarted = TRUE;
        m_uiStart = (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    }

    void Stop()
    {
        if (m_bStarted)
        {
            m_bStarted = FALSE;
            ++m_dwCounter;
            m_fElapsed += ((std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count()) - m_uiStart) * 0.001f;
        }
    }

    void Reset()
    {
        m_bStarted = FALSE;
        m_fElapsed = 0;
        m_dwCounter = 0;
    }

    FLOAT Elapsed(void) const { return m_fElapsed; }
    DWORD GetCounter() const { return m_dwCounter; }

    void Report(const EVerboseLevel vl = GENERAL)
    {
        Stop();

        const DWORD dwCount = GetCounter();
        const FLOAT fElapsed = Elapsed();
        const FLOAT fAverage = (0 != dwCount) ? (fElapsed / dwCount) : 0;

        appVOut(vl, _T("Elapsed time: %s: total %12.2f sec, count %4d, average %12.2f sec\n")
            , m_sId.c_str()
            , fElapsed
            , dwCount
            , fAverage);
    }

private:

    UBOOL m_bStarted;
    ULONGLONG m_uiStart;
    FLOAT m_fElapsed;
    DWORD m_dwCounter;
    CCString m_sId;
    UBOOL m_bReportOnExit;
};

__END_NAMESPACE

#endif //#ifndef _TIMER_H_


//=============================================================================
// END OF FILE
//=============================================================================
