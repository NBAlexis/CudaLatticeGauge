//=============================================================================
// FILENAME : CField.h
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CFIELD_H_
#define _CFIELD_H_

__BEGIN_NAMESPACE

class CLGAPI CField
{
public:

    CField(CLatticeData* pLattice)
        : m_pOwner(pLattice) 
        , m_uiDir(pLattice->m_uiDir)
    { 
        memcpy(m_uiLatticeLength, pLattice->m_uiLatticeLength, sizeof(UINT) * CCommonData::kMaxDim);
        memcpy(m_uiLatticeDecompose, pLattice->m_uiLatticeDecompose, sizeof(UINT) * CCommonData::kLatticeDecompose * 2);
        memcpy(m_uiLatticeMultipy, pLattice->m_uiLatticeMultipy, sizeof(UINT) * (CCommonData::kMaxDim + 1));
    }

    virtual EFieldType GetFieldType() const = 0;

#pragma region BLAS
    //what is BLAS? see: https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms

    virtual void axpy(FLOAT a, const CField *x) = 0;
    virtual void axpy(const cuComplex& a, const CField *x) = 0;

#pragma endregion BLAS

protected:

    CLatticeData * m_pOwner;
    UINT m_uiDir;
    UINT m_uiLatticeLength[CCommonData::kMaxDim];
    UINT m_uiLatticeDecompose[CCommonData::kLatticeDecompose * 2];
    UINT m_uiLatticeMultipy[CCommonData::kMaxDim + 1];
};

__END_NAMESPACE

#endif //#ifndef _CFIELD_H_

//=============================================================================
// END OF FILE
//=============================================================================