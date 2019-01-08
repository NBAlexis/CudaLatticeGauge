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

__DEFINE_ENUM(EFieldInitialType,

    EFIT_Zero,
    EFIT_Identity,
    EFIT_Random,
    EFIT_RandomGenerator,
    EFIT_ReadFromFile,

    EFIT_ForceDWORD = 0x7fffffff,
    )


class CLGAPI CField : public CBase
{
public:

    CField();

    virtual EFieldType GetFieldType() const = 0;
    virtual void InitialField(EFieldInitialType eInitialType) = 0;
    
    virtual void DebugPrintMe() const = 0;

#pragma region BLAS
    //what is BLAS? see: https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms

    virtual void Zero() = 0;
    virtual void Indentity() = 0;

    //This is Axpy(1.0f, x)
    virtual void Axpy(const CField* x) = 0;
    virtual void Axpy(Real a, const CField* x) = 0;
    virtual void Axpy(const _Complex& a, const CField* x) = 0;
    
#pragma endregion BLAS

#pragma region HMC

    /**
    * U = exp(a this)U
    */
    virtual void ExpMult(const _Complex& a, CField* U) const = 0;

    virtual void CopyTo(CField* U) const = 0;

#pragma endregion HMC

    class CLatticeData* m_pOwner;
    UINT m_uiThreadCount;
};

__END_NAMESPACE

#endif //#ifndef _CFIELD_H_

//=============================================================================
// END OF FILE
//=============================================================================