//=============================================================================
// FILENAME : CBoundaryConditionDirichlet.h
// 
// DESCRIPTION:
// This is the class for Torus(usually spatial) and Dirichlet(usually temporal) boundary condition
//
// REVISION:
//  [12/4/2018 nbale]
//=============================================================================

#ifndef _CLATTICEDATA_H_
#define _CLATTICEDATA_H_

__BEGIN_NAMESPACE

class CLGAPI CLatticecData
{
public:
    enum { kMaxDim = 4, kMaxFieldCount = 8,};

    CLatticecData(UINT dim, UINT lattice[], UINT dir, STRING fields[]);

    UINT GetDim() const
    {
        return m_uiDim;
    }

    UINT GetVolumn() const
    {
        return m_uiVolumn;
    }

    UINT GetDir() const
    {
        return m_uiDir;
    }

    UINT GetLatticeLength(UINT dim) const
    {
        appAssert(dim < m_uiDim);
        return m_uiLattice[dim];
    }

private:

    UINT m_uiDim;
    UINT m_uiVolumn;
    UINT m_uiDir;
    UINT m_uiLattice[kMaxDim];
    STRING m_sFields[kMaxFieldCount];
    CField* m_pFields[kMaxFieldCount];
};

__END_NAMESPACE

#endif //#ifndef _CLATTICEDATA_H_

//=============================================================================
// END OF FILE
//=============================================================================