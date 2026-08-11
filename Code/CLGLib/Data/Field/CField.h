//=============================================================================
// FILENAME : CField.h
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
//
// REVISION:
//  [mm/dd/yy]
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CFIELD_H_
#define _CFIELD_H_

//Improve-1 (multi-GPU-improve1.md 3.1/3.2): buffer-extent halo identity; every
//field MAY expose its managed device buffer through a CHaloBufferHandle.
#include "Core/Distributed/CHaloBufferHandle.h"

#define _GetData \
const void* GetData() const override \
{ \
    return (const void*)m_pDeviceData; \
} \
void* GetData() override \
{ \
    return (void*)m_pDeviceData; \
}

__BEGIN_NAMESPACE

__DEFINE_ENUM(EFieldFileType,
    EFFT_BridgePPTXT,
    EFFT_BridgePPBin,
    EFFT_CLGBin,
    EFFT_CLGBinCompressed,
    EFFT_CLGBinFloat,
    EFFT_CLGBinDouble,
    EFFT_CLGBinSU3_12,

    EFFT_ForceDWORD = 0x7fffffff,
    )

__DEFINE_ENUM(EFieldOperator,

    EFO_F_D,
    EFO_F_Ddagger,
    EFO_F_DD,
    EFO_F_DDdagger,
    EFO_F_InverseD,
    EFO_F_InverseDdagger,
    EFO_F_InverseDD,
    EFO_F_InverseDDdagger,

    EFO_F_RationalD,

    EFO_F_D_WithMass,
    EFO_F_Ddagger_WithMass,
    EFO_F_DD_WithMass,
    EFO_F_DDdagger_WithMass,

    EFO_ForceDWORD = 0x7fffffff,
    )


/**
* When apply operator, it is convenient to apply a coefficient at the same time
*/
enum EOperatorCoefficientType
{
    EOCT_None, // (1.0)
    EOCT_Minus, // (-1.0)
    EOCT_Real, // Real Number
    EOCT_Complex, //Complex Number
};

__DEFINE_ENUM(EFermionBosonSource,
    EFS_Point,
    EFS_Wall,
    EFS_StaggeredWall,
    EFS_StaggeredZWall,
    EFS_MomentumWall,
    )

enum EFixBoundary
{
    EFB_Field,
    EFB_Momentum,
    EFB_Force,
};

struct SFermionBosonSource
{
    EFermionBosonSource m_eSourceType;
    SSmallInt4 m_sSourcePoint;
    //For staggered fermion, spin index is 0-7, for shift of x,y,z
    BYTE m_bySpinIndex;
    BYTE m_byColorIndex;
    CLGComplex m_cOtherParameters1;
    CLGComplex m_cOtherParameters2;
    CLGComplex m_cOtherParameters3;
    CLGComplex m_cOtherParameters4;
};

class CLGAPI CField : public CBase
{
public:

    CField();

    virtual EFieldType GetFieldType() const = 0;
    virtual void InitialField(EFieldInitialType eInitialType) = 0;
    virtual void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFile = EFFT_CLGBin) = 0;
    virtual void InitialWithByte(BYTE* byData) = 0;
    virtual void InitialWithByteCompressed(const CCString& sFileName) { appCrucial(_T("Not implemented compressed file format!\n")); }
    virtual void InitialOtherParameters(CParameters& param) 
    {
        m_pClass = GetClass();
        param.FetchValueArrayBYTE(_T("GaugeFields"), m_byGaugeFieldIds);
        param.FetchValueArrayBYTE(_T("BosonFields"), m_byBosonFieldIds);

        if (0 == m_byGaugeFieldIds.Num() && 0 == m_byBosonFieldIds.Num())
        {
            m_byGaugeFieldIds.AddItem(1);
        }

        INT iDynamic = 1;
        param.FetchValueINT(_T("Dynamic"), iDynamic);
        m_bDynamic = (0 != iDynamic);
    }

    UBOOL SingleField() const
    {
        return 1 == m_byGaugeFieldIds.Num() && 0 == m_byBosonFieldIds.Num();
    }

    virtual void DebugPrintMe() const = 0;
    virtual const void* GetData() const = 0;
    virtual void* GetData() = 0;

    BYTE GetFieldId() const { return m_byFieldId; }

    #pragma region BLAS
    //what is BLAS? see: https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms

    virtual void Zero()
    {
        InitialField(EFIT_Zero);
    }

    virtual void ZeroOnEvenOdd(UBOOL bEven)
    {
        appCrucial(_T("ZeroOnEvenOdd not implemented!\n"));
    }

    virtual void Identity()
    {
        InitialField(EFIT_Identity);
    }

    virtual void FixBoundary(EFixBoundary eType) {}

    //This is Axpy(1.0f, x)
    virtual void AxpyPlus(const CField* x) = 0;
    //This is Axpy(-1.0f, x)
    virtual void AxpyMinus(const CField* x) = 0;

    //me = me + a * x
    virtual void Axpy(Real a, const CField* x) = 0;
    virtual void Axpy(const CLGComplex& a, const CField* x) = 0;
    virtual void Mul(const CField* other, UBOOL bDaggerLeft = TRUE, UBOOL bDaggerRight = FALSE) = 0;
    virtual void LeftMul(const CField* other, UBOOL bDaggerLeft = FALSE, UBOOL bDaggerRight = FALSE) = 0;
    virtual void Dagger() = 0;

    virtual void ApplyPhaseC(const CField* other)
    {
        appCrucial(_T("ApplyPhaseC not implemented"));
    }
    virtual void ApplyPhaseR(const CField* other, Real fCharge)
    {
        appCrucial(_T("ApplyPhaseR not implemented"));
    }

    //This is a * me
    virtual void ScalarMultply(const CLGComplex& a) = 0;
    virtual void ScalarMultply(Real a) = 0;
    virtual CCString SaveToFile(const CCString &fileName, EFieldFileType eType = EFFT_CLGBin) const;
    virtual CCString SaveToCompressedFile(const CCString&) const { appCrucial(_T("Not supported compressed file format for this field!\n")); return _T("Not Supported"); }
    //Why we need this? because the data structure are aligned.
    //Remember to free the buffer
    virtual BYTE* CopyDataOut(UINT &uiSize) const = 0;
    virtual BYTE* CopyDataOutFloat(UINT& uiSize) const = 0;
    virtual BYTE* CopyDataOutDouble(UINT& uiSize) const = 0;
    CCString GetInfos(const CCString &tab) const override;

    virtual UBOOL IsDynamic() const
    {
        return m_bDynamic;
    }

    virtual UBOOL IsDirichlet() const
    {
        return FALSE;
    }

#pragma endregion

#pragma region Other useful operators

    /**
    * return me^* \cdot other
    * pDeviceBuffer is a Real array, with length of [thread count]
    * Using pDeviceBuffer, we make sure Dot function is a constant function as it should be.
    * The final result of dot, should be sum of pDeviceBuffer
    */
    virtual cuDoubleComplex Dot(const CField* other) const = 0;

    /**
    * Squared Length! this is just self dot self
    */
    virtual DOUBLE GetLength() const = 0;

    virtual CLGComplex DotReal(const CField* other) const
    {
#if _CLG_DOUBLEFLOAT
        return Dot(other);
#else
        return _cToFloat(Dot(other));
#endif
    }

    virtual CField* GetCopy() const = 0;

    virtual UBOOL ApplyOperator(EFieldOperator op, INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* pGauge, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields, EOperatorCoefficientType uiCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0), void* otherParameter = NULL) = 0;

    virtual UBOOL IsGaugeField() const { return FALSE; }
    virtual UBOOL IsFermionField() const { return FALSE; }
    virtual UBOOL IsBosonField() const { return FALSE; }
    virtual UBOOL IsSpinField() const { return FALSE; }
    virtual UBOOL IsTensor2Field() const { return FALSE; }

    /**
    * Improve-1 (multi-GPU-improve1.md 3.2): the handle describing this field's
    * managed device buffer extent. NULL means the object is NOT under the
    * automatic halo protocol (and must not be treated as a supported
    * LocalOnly buffer in a plain launch). Halo-capable template bases hold
    * and Bind one; pool copies each Bind their own extent.
    */
    virtual CHaloBufferHandle* GetHaloBufferHandle() { return NULL; }
    virtual const CHaloBufferHandle* GetHaloBufferHandle() const { return NULL; }

    /**
    * Improve-1 (multi-GPU-improve1.md 3.5, I4): every host-side or BLAS write
    * entry (InitialField / InitialWithByte / file load / CopyTo / CopyBufferTo
    * target / BLAS mutators / pointer rebind) MUST end with this, so the next
    * Ensure sees the write. Over-marking only costs an extra halo exchange;
    * under-marking is a correctness bug. No-op for objects without a bound
    * handle, so single-GPU and not-yet-wired fields are unaffected.
    */
    void NotifyWritten()
    {
        CHaloBufferHandle* pHandle = GetHaloBufferHandle();
        if (NULL != pHandle && pHandle->IsBound())
        {
            pHandle->NotifyWritten();
        }
    }

#pragma endregion

    class CLatticeData* m_pOwner;
    BYTE m_byFieldId;
    UBOOL m_bDynamic;
    DOUBLE m_fLength;

    void Return();

    friend class CFieldPool;

    //void UpdatePooledParamters() const;

    virtual void Connection(const CField* other, void* res) const
    {
        appCrucial(_T("Connect not implemented\n"));
    }

    virtual void ConnectionSelf(void* res, Real fCoeff) const
    {
        appCrucial(_T("ConnectionSelf not implemented\n"));
    }

    virtual void AddConnectionSelf(void* res, Real fCoeff) const
    {
        appCrucial(_T("AddConnectionSelf not implemented\n"));
    }

    virtual void AddConnectionTwoField(void* res, const CField* other, Real fCoeff) const
    {
        appCrucial(_T("AddConnection not implemented\n"));
    }

protected:

    const CFieldGauge* GetDefaultGauge(INT gaugeNum, const CFieldGauge* const* gaugeFields) const
    {
        if (0 == gaugeNum || m_byGaugeFieldIds.Num() < 1)
        {
            return NULL;
        }
        if (NULL == gaugeFields)
        {
            return NULL;
        }
        INT gaugeIdx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, gaugeFields, m_byGaugeFieldIds[0]);
        if (gaugeIdx < 0 || gaugeIdx >= gaugeNum)
        {
            return NULL;
        }
        return gaugeFields[gaugeIdx];
    }

    //class CFieldPool* m_pPool;

    TArray<BYTE> m_byGaugeFieldIds;
    TArray<BYTE> m_byBosonFieldIds;

#pragma region pool

public:

    virtual void CopyParamTo(CField* U) const
    {
        appAssert(NULL != U);
        U->m_pClass = m_pClass;
        U->m_pOwner = m_pOwner;
        U->m_byFieldId = m_byFieldId;
        U->m_bDynamic = m_bDynamic;
        U->m_fLength = m_fLength;

        U->m_byGaugeFieldIds = m_byGaugeFieldIds;
        U->m_byBosonFieldIds = m_byBosonFieldIds;

        //Improve-1: the pool copy's own handle stays bound to ITS buffer;
        //only the tag id follows the field id (identity is the buffer).
        if (NULL != U->GetHaloBufferHandle())
        {
            U->GetHaloBufferHandle()->SetFieldId(m_byFieldId);
        }
    }

    virtual void CopyBufferTo(CField* U) const = 0;

    virtual void CopyTo(CField* U) const
    {
        CopyParamTo(U);
        CopyBufferTo(U);
    }

    const CClass* m_pClass;

    //only use for debug the pool
    //const CFieldPool* ShowPool() const { return m_pPool; }

#pragma endregion

};

class CLGAPI CFieldPool
{

public:

    struct CLGAPI CPooledFields
    {
        UBOOL m_bInUse;
        CField* m_pField;

        UBOOL operator==(const CPooledFields& Other) const
        {
            return m_pField == Other.m_pField;
        }
    };

    CFieldPool()
    {

    }

    virtual ~CFieldPool()
    {
        FreeFields();
    }

    void FreeFields()
    {
        TArray<const CClass*> keys = m_pPool.GetAllKeys();
        for (INT i = 0; i < keys.Num(); ++i)
        {
            TArray<CPooledFields>& fields = m_pPool[keys[i]];
            for (INT j = 0; j < fields.Num(); ++j)
            {
                appSafeDelete(fields[j].m_pField);
            }
        }
        m_pPool.RemoveAll();
    }

    CField* GetOne(const CField* pOrignal);
    void Return(CField* pField);

    THashMap<const CClass*, TArray<CPooledFields>> m_pPool;
};

/**
* NOTE: CFieldCache is not tested for a long time
* To be removed
*/
class CLGAPI CFieldCache
{
public:
    enum ECacheReason
    {
        CachedInverseDDdaggerField = 1,

        //the real ID = 100 + action ID
        CachedForceFieldStart = 100,
    };

    CFieldCache()
    {

    }
    ~CFieldCache()
    {
        for (INT i = 0; i < m_pCachedFields.Num(); ++i)
        {
            appSafeDelete(m_pCachedFields[i]);
        }
    }

    UBOOL CacheField(UINT uiID, CField* pField)
    {
        if (!m_pCachedFieldMaps.Exist(uiID))
        {
            m_pCachedFields.AddItem(pField);
            m_pCachedFieldMaps.SetAt(uiID, pField);
            return TRUE;
        }
        return FALSE;
    }

    CField* GetCachedField(UINT uiID)
    {
        if (m_pCachedFieldMaps.Exist(uiID))
        {
            return m_pCachedFieldMaps[uiID];
        }
        return NULL;
    }

    TArray<CField*> m_pCachedFields;
    THashMap<UINT, CField*> m_pCachedFieldMaps;
};

/**
 * Maybe better to implement as a template?
 * No ..., for different inherent, the main work is just the kernel.
 * 
 * NOTE: CFieldMatrixOperation is implimented for deflation-restarted solver, not tested for a long time
 */
 class CLGAPI CFieldMatrixOperation
{
public:
    virtual ~CFieldMatrixOperation() {}
    enum {_kFieldMatrixMaxDim = 32 };
    static CFieldMatrixOperation* Create(EFieldType ef);

    virtual void VectorMultiplyMatrix(TArray<CField*>& res, const TArray<CField*>& left, const CLGComplex* deviceMatrix, UINT uiDimX, UINT uiDimY) = 0;
};

__END_NAMESPACE

#endif //#ifndef _CFIELD_H_

//=============================================================================
// END OF FILE
//=============================================================================