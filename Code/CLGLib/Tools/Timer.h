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
        , m_fStart(0.0F)
        , m_fElapsed(0.0F)
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
        m_fStart = appGetTime();
    }

    void Stop()
    {
        if (m_bStarted)
        {
            m_bStarted = FALSE;
            ++m_dwCounter;
            m_fElapsed += (appGetTime() - m_fStart);
        }
    }

    void Reset()
    {
        m_bStarted = FALSE;
        m_fElapsed = 0.0f;
        m_dwCounter = 0;
    }

    FLOAT Elapsed(void) const { return m_fElapsed; }
    DWORD GetCounter() const { return m_dwCounter; }

    void Report(const EVerboseLevel vl = GENERAL)
    {
        Stop();

        DWORD dwCount = GetCounter();
        FLOAT fElapsed = Elapsed();
        FLOAT fAverage = (0 != dwCount) ? (fElapsed / dwCount) : 0.0f;

        appVOut(vl, _T("Elapsed time: %s: total %12.2f sec, count %4d, average %12.2f sec\n")
            , m_sId.c_str()
            , fElapsed
            , dwCount
            , fAverage);
    }

private:

    BOOL m_bStarted;
    FLOAT m_fStart;
    FLOAT m_fElapsed;
    DWORD m_dwCounter;
    CCString m_sId;
    BOOL m_bReportOnExit;
};

__END_NAMESPACE

#endif //#ifndef _TIMER_H_


//=============================================================================
// END OF FILE
//=============================================================================
