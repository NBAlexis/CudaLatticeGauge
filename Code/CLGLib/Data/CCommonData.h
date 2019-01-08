//=============================================================================
// FILENAME : CCommonData.h
// 
// DESCRIPTION:
// This is the class for the common data
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CCOMMONDATA_H_
#define _CCOMMONDATA_H_

#pragma region constants

//DC for device constant
//HC for host constant

#define _DC_Dim (_constIntegers[ECI_Dim])
#define _HC_Dim (appGetCudaHelper()->m_ConstIntegers[ECI_Dim])

#define _DC_Dir (_constIntegers[ECI_Dir])
#define _HC_Dir (appGetCudaHelper()->m_ConstIntegers[ECI_Dir])

#define _DC_Lx (_constIntegers[ECI_Lx])
#define _HC_Lx (appGetCudaHelper()->m_ConstIntegers[ECI_Lx])
#define _DC_Ly (_constIntegers[ECI_Ly])
#define _HC_Ly (appGetCudaHelper()->m_ConstIntegers[ECI_Ly])
#define _DC_Lz (_constIntegers[ECI_Lz])
#define _HC_Lz (appGetCudaHelper()->m_ConstIntegers[ECI_Lz])
#define _DC_Lt (_constIntegers[ECI_Lt])
#define _HC_Lt (appGetCudaHelper()->m_ConstIntegers[ECI_Lt])

#define _DC_Volumn (_constIntegers[ECI_Volumn])
#define _HC_Volumn (appGetCudaHelper()->m_ConstIntegers[ECI_Volumn])

#define _DC_MultX (_constIntegers[ECI_MultX])
#define _HC_MultX (appGetCudaHelper()->m_ConstIntegers[ECI_MultX])
#define _DC_MultY (_constIntegers[ECI_MultY])
#define _HC_MultY (appGetCudaHelper()->m_ConstIntegers[ECI_MultY])
#define _DC_MultZ (_constIntegers[ECI_MultZ])
#define _HC_MultZ (appGetCudaHelper()->m_ConstIntegers[ECI_MultZ])

#define _DC_DecompX (_constIntegers[ECI_DecompX])
#define _HC_DecompX (appGetCudaHelper()->m_ConstIntegers[ECI_DecompX])
#define _DC_DecompY (_constIntegers[ECI_DecompY])
#define _HC_DecompY (appGetCudaHelper()->m_ConstIntegers[ECI_DecompY])
#define _DC_DecompZ (_constIntegers[ECI_DecompZ])
#define _HC_DecompZ (appGetCudaHelper()->m_ConstIntegers[ECI_DecompZ])
#define _DC_DecompLx (_constIntegers[ECI_DecompLx])
#define _HC_DecompLx (appGetCudaHelper()->m_ConstIntegers[ECI_DecompLx])
#define _DC_DecompLy (_constIntegers[ECI_DecompLy])
#define _HC_DecompLy (appGetCudaHelper()->m_ConstIntegers[ECI_DecompLy])
#define _DC_DecompLz (_constIntegers[ECI_DecompLz])
#define _HC_DecompLz (appGetCudaHelper()->m_ConstIntegers[ECI_DecompLz])
#define _DC_ThreadCount (_constIntegers[ECI_ThreadCount])
#define _HC_ThreadCount (appGetCudaHelper()->m_ConstIntegers[ECI_ThreadCount])

#define _DC_Seed (_constIntegers[ECI_RandomSeed])
#define _HC_Seed (appGetCudaHelper()->m_ConstIntegers[ECI_RandomSeed])
#define _DC_ExpPrecision (_constIntegers[ECI_ExponentPrecision])
#define _HC_ExpPrecision (appGetCudaHelper()->m_ConstIntegers[ECI_ExponentPrecision])

#define _DC_ActionListL (_constIntegers[ECI_ActionListLength])
#define _HC_ActionListL (appGetCudaHelper()->m_ConstIntegers[ECI_ActionListLength])

#define _D_ComplexThreadBuffer (appGetCudaHelper()->m_pComplexBufferThreadCount)
#define _D_RealThreadBuffer (appGetCudaHelper()->m_pRealBufferThreadCount)

#pragma endregion

__BEGIN_NAMESPACE

enum { kMaxFieldCount = 8, kMaxActionCount = 8 };

__DEFINE_ENUM(EFieldType,

    EFT_GaugeSU3,
    EFT_FermionWilsonSquareSU3,
    EFT_Max,
    EFT_ForceDword = 0x7fffffff,

    )



//====================================================
// Some common structures
//====================================================

struct CLGAPI SIndex
{
    enum
    {
        kDagger = 0x01,
    };

    __device__ SIndex()
        : m_uiSiteIndex(0)
        , m_byDir(0)
        , m_byTag(0)
        , m_byBoundaryFieldId(0)
    {

    }

    __device__ SIndex(UINT uiIndex, BYTE dir = 0, BYTE indexTag = 0, BYTE bcField = 0)
        : m_uiSiteIndex(uiIndex)
        , m_byDir(dir)
        , m_byTag(indexTag)
        , m_byBoundaryFieldId(bcField)
    {

    }

    __device__ __inline__ UBOOL NeedToDagger() const { return 0 != (kDagger & m_byTag); }
    __device__ __inline__ UBOOL NeedBoundaryField() const { return 0 != m_byBoundaryFieldId; }

    UINT m_uiSiteIndex;
    BYTE m_byDir;
    BYTE m_byTag;
    BYTE m_byBoundaryFieldId;
};


__END_NAMESPACE

#endif //#ifndef _CCOMMONDATA_H_

//=============================================================================
// END OF FILE
//=============================================================================