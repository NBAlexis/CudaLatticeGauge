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
    EFIT_RandomGaussian,
    EFIT_ReadFromFile,

    EFIT_ForceDWORD = 0x7fffffff,
    )



__DEFINE_ENUM(EFieldOperator,

    EFO_F_D,
    EFO_F_Ddagger,
    EFO_F_DDdagger,
    EFO_F_InverseDDdagger,

    EFO_ForceDWORD = 0x7fffffff,
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
    virtual void AxpyPlus(const CField* x) = 0;
    //This is Axpy(1.0f, x)
    virtual void AxpyMinus(const CField* x) = 0;

    virtual void Axpy(Real a, const CField* x) = 0;
    virtual void Axpy(const _Complex& a, const CField* x) = 0;

    //This is a * me
    virtual void ScalarMultply(const _Complex& a) = 0;
    virtual void ScalarMultply(Real a) = 0;
    
#pragma endregion BLAS

#pragma region Other useful operators

    /**
    * pDeviceBuffer is a Real array, with length of [thread count]
    * Using pDeviceBuffer, we make sure Dot function is a constant function as it should be.
    * The final result of dot, should be sum of pDeviceBuffer
    */
    virtual _Complex Dot(const CField* other) const = 0;

    /**
    * U = exp(a this)U
    */
    virtual void ExpMult(const _Complex& a, CField* U) const = 0;

    virtual void CopyTo(CField* U) const = 0;
    virtual CField* GetCopy() const = 0;
    virtual CField* GetZero() const = 0;

    virtual void ApplyOperator(EFieldOperator op, const CField* otherfield) = 0;

#pragma endregion HMC

    class CLatticeData* m_pOwner;
    BYTE m_byFieldId;
};

__END_NAMESPACE

#endif //#ifndef _CFIELD_H_

//=============================================================================
// END OF FILE
//=============================================================================