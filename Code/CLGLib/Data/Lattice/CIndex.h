//=============================================================================
// FILENAME : CIndex.h
// 
// DESCRIPTION:
// This is the class for index on lattice
//
// REVISION:
//  [12/5/2018 nbale]
//=============================================================================

#ifndef _CINDEX_H_
#define _CINDEX_H_

__BEGIN_NAMESPACE

class CLGAPI CIndex
{

public:
    enum { kSpace = 0x01, kTime = 0x10, kSpaceTime = 0x11, };

    CIndex(class CLatticeData * pOwner) : m_pOwner(pOwner), m_pBoundaryCondition(NULL) { ; }

    class CLatticeData * GetOwner() const { return m_pOwner; }

    class CBoundaryCondition * GetBoundaryCondition() const { return m_pBoundaryCondition; }
    void SetBoundaryCondition(class CBoundaryCondition * pBoundaryCondition) { m_pBoundaryCondition = pBoundaryCondition; }

    /**
    * WOW, we can use virtual device functions!
    * This function to get the plaquttes indexes.
    * one can specify to get only spatial plaquttes or spacetime
    */
    virtual __device__ int2* GetPlaquttesAtLink(UINT& count, UINT& plaqutteLength, UINT uiDim, UINT uiLinkIndex, const UINT* length, const UINT* mult, UINT st = kSpaceTime) = 0;

protected:

    class CLatticeData * m_pOwner;
    class CBoundaryCondition * m_pBoundaryCondition;
};

__END_NAMESPACE

#endif //#ifndef _CINDEX_H_

//=============================================================================
// END OF FILE
//=============================================================================