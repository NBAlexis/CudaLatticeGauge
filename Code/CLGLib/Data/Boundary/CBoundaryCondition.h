//=============================================================================
// FILENAME : CBoundaryCondition.h
// 
// DESCRIPTION:
// This is the class for boundary conditions
// Note that, the boundary conditions should only make sense together with lattice!!
//
// REVISION:
//  [12/5/2018 nbale]
//=============================================================================

#ifndef _CBOUNDARYCONDITION_H_
#define _CBOUNDARYCONDITION_H_

__BEGIN_NAMESPACE

class CLGAPI CBoundaryCondition
{
public:

    CBoundaryCondition(class CIndex * pOwner) : m_pOwner(pOwner) { ; }
    class CIndex * GetOwner() const { return m_pOwner; }

    /**
    * The first index is the site index, the second index is index of field, it is 0 if it is not on boundary
    */
    virtual __device__ uint2 GetMappedIndex(const int4 &site, const int4 &fromsite, const UINT* length, const UINT* mult) = 0;

protected:
    class CIndex * m_pOwner;
};

__END_NAMESPACE

#endif //#ifndef _CBOUNDARYCONDITION_H_

//=============================================================================
// END OF FILE
//=============================================================================