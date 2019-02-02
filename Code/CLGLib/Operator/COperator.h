//=============================================================================
// FILENAME : COperator.h
// 
// DESCRIPTION:
// 
// For better readability, we use an operator class to do all operations on fields.
// NOTE: 
// Only operators with physical meaning should be implemented with COperator, for example, 
// things like Axpy should NOT.
// things like normalize, make zero, etc, should NOT.
// things like D should, should it?
//
// REVISION:
//  [01/27/2019 nbale]
//=============================================================================

#ifndef _COPERATOR_H_
#define _COPERATOR_H_

__BEGIN_NAMESPACE

__DEFINE_ENUM(EOperatorType,

    Single,
    SingleWithComplex,
    SingleWithReal,

    Binary,
    BinaryWithComplex,
    BinaryWithReal,

    EOT_ForceDWORD = 0x7fffffff,
    )


__DEFINE_ENUM(EOperator,

        EO_MakeZero,
        EO_MakeId,
        EO_MakeGaussianRandom,
        EO_Add,
        EO_Sub,
        EO_Axpy,

        EO_ForceDWORD = 0x7fffffff,
        )


struct CLGAPI SOperatorCalculateInfo
{
    CField* m_pResult;
    CField* m_pFirst;
    CField* m_pSecond;

    _Complex m_c;
    Real m_r;
};

class CLGAPI COperator : public CBase
{
public:

    COperator();
    
    virtual void Calculate(const SOperatorCalculateInfo& sinfo) = 0;
};

__END_NAMESPACE

#endif //#ifndef _COPERATOR_H_

//=============================================================================
// END OF FILE
//=============================================================================