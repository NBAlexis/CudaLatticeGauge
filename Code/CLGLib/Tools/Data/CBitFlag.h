//=============================================================================
// FILENAME : CBitFlag.h
// 
// DESCRIPTION:
//
// REVISION:
//  [3/13/2018 nbale]
//=============================================================================

#ifndef _CBITFLAG_H_
#define _CBITFLAG_H_

__BEGIN_NAMESPACE

class CLGAPI CCBitFlag
{
protected:
    QWORD m_nBits;
public:
    CCBitFlag(QWORD bits=0):m_nBits(bits){}
    QWORD GetValue() const { return m_nBits; }
    void SetValue(const QWORD &inValue) { m_nBits = inValue; }

    void SetFlag(QWORD flag) { m_nBits |= flag; }
    void ClearFlag(QWORD flag) { m_nBits &= (~flag); }
    UBOOL HasFlag(QWORD flag) const { return (m_nBits&flag) ? TRUE : FALSE; }
    void ToggleFlagBy(QWORD flag, UBOOL bSetFlag) {bSetFlag ? SetFlag(flag) : ClearFlag(flag);}
    CCBitFlag& operator = (const QWORD& inValue) { m_nBits = inValue; return *this; }

    operator QWORD() const { return m_nBits; }
};

__END_NAMESPACE

#endif //#ifndef _CBITFLAG_H_

//=============================================================================
// END OF FILE
//=============================================================================