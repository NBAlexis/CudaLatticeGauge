//=============================================================================
// FILENAME : CField.h
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _COPERATOR_H_
#define _COPERATOR_H_

__BEGIN_NAMESPACE

__DEFINE_ENUM(EOperatorType,

    Single,
    SingleWithComplex,
    SingleWithReal,

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

class CLGAPI COperator : public CBase
{
public:

    COperator();
    
};

__END_NAMESPACE

#endif //#ifndef _COPERATOR_H_

//=============================================================================
// END OF FILE
//=============================================================================