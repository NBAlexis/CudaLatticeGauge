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

__DEFINE_ENUM(EBoundaryCondition,
    EBC_TorusSquare,
    EBC_Max,
    EBC_ForceDWORD = 0x7fffffff,
    )

class CLGAPI CBoundaryCondition : public CBase
{
public:

    CBoundaryCondition() : m_pOwner(NULL)
    {
        ;
    }

    /**
    * The first index is the site index, the second index is index of field, it is 0 if it is not on boundary
    */
    __device__ virtual uint2 _devcieGetMappedIndex(const int4 &site, const int4 &fromsite) const = 0;

    class CLatticeData* m_pOwner;
};

__END_NAMESPACE

#endif //#ifndef _CBOUNDARYCONDITION_H_

//=============================================================================
// END OF FILE
//=============================================================================