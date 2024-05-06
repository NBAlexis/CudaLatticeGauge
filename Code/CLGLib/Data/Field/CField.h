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

__DEFINE_ENUM(EFieldInitialType,

    EFIT_Zero,
    EFIT_Identity,
    EFIT_Random,
    EFIT_RandomGenerator,
    EFIT_SumGenerator, //for testing use only
    EFIT_RandomGaussian,
    EFIT_RandomZ4,

    EFIT_ReadFromFile,

    EFIT_U1Real,

    EFIT_ForceDWORD = 0x7fffffff,
    )


__DEFINE_ENUM(EFieldFileType,
    EFFT_BridgePPTXT,
    EFFT_BridgePPBin,
    EFFT_CLGBin,
    EFFT_CLGBinCompressed,
    EFFT_CLGBinFloat,
    EFFT_CLGBinDouble,
    
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

class CLGAPI CField : public CBase
{
public:

    CField();

    virtual EFieldType GetFieldType() const = 0;
    virtual void InitialField(EFieldInitialType eInitialType) = 0;
    virtual void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFile) = 0;
    virtual void InitialWithByte(BYTE* byData) = 0;
    virtual void InitialWithByteCompressed(BYTE* ) { appCrucial(_T("Not implemented compressed file format!\n")); }
    virtual void InitialOtherParameters(CParameters& param) 
    {
        param.FetchValueArrayBYTE(_T("GaugeFields"), m_byGaugeFieldIds);
        param.FetchValueArrayBYTE(_T("BosonFields"), m_byBosonFieldIds);

        if (0 == m_byGaugeFieldIds.Num() && 0 == m_byBosonFieldIds.Num())
        {
            m_byGaugeFieldIds.AddItem(1);
        }
    }

    UBOOL SingleField() const
    {
        return 1 == m_byGaugeFieldIds.Num() && 0 == m_byBosonFieldIds.Num();
    }

    virtual void DebugPrintMe() const = 0;
    virtual const void* GetData() const = 0;
    virtual void* GetData() = 0;

#pragma region BLAS
    //what is BLAS? see: https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms

    virtual void Zero()
    {
        InitialField(EFIT_Zero);
    }

    virtual void Identity()
    {
        InitialField(EFIT_Identity);
    }

    virtual void FixBoundary() {}

    //This is Axpy(1.0f, x)
    virtual void AxpyPlus(const CField* x) = 0;
    //This is Axpy(1.0f, x)
    virtual void AxpyMinus(const CField* x) = 0;

    virtual void Axpy(Real a, const CField* x) = 0;
    virtual void Axpy(const CLGComplex& a, const CField* x) = 0;
    virtual void Dagger() = 0;

    DOUBLE GetLength() const { return m_fLength; }

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
    virtual CCString GetInfos(const CCString &tab) const;

#pragma endregion

#pragma region Other useful operators

    /**
    * return me^* \cdot other
    * pDeviceBuffer is a Real array, with length of [thread count]
    * Using pDeviceBuffer, we make sure Dot function is a constant function as it should be.
    * The final result of dot, should be sum of pDeviceBuffer
    */
    virtual cuDoubleComplex Dot(const CField* other) const = 0;

    virtual CLGComplex DotReal(const CField* other) const
    {
#if _CLG_DOUBLEFLOAT
        return Dot(other);
#else
        return _cToFloat(Dot(other));
#endif
    }

    virtual void CopyTo(CField* U) const
    {
        assert(NULL != U);
        U->m_pOwner = m_pOwner;
        U->m_byFieldId = m_byFieldId;
        U->m_fLength = m_fLength;

        U->m_byGaugeFieldIds = m_byGaugeFieldIds;
        U->m_byBosonFieldIds = m_byBosonFieldIds;
    }

    virtual CField* GetCopy() const = 0;

    virtual UBOOL ApplyOperator(EFieldOperator op, INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, const CFieldBoson* const* pBoson, EOperatorCoefficientType uiCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0), void* otherParameter = NULL) = 0;

    virtual UBOOL IsGaugeField() const { return FALSE; }
    virtual UBOOL IsFermionField() const { return FALSE; }
    virtual UBOOL IsBosonField() const { return FALSE; }
    virtual UBOOL IsSpinField() const { return FALSE; }

#pragma endregion

    class CLatticeData* m_pOwner;
    BYTE m_byFieldId;
    DOUBLE m_fLength;

    void Return();

    friend class CFieldPool;

protected:

    const CFieldGauge* GetDefaultGauge(INT gaugeNum, const CFieldGauge* const* gaugeFields) const
    {
        return gaugeFields[CLatticeData::GetGaugeFieldIndexById(gaugeNum, gaugeFields, m_byGaugeFieldIds[0])];
    }

    class CFieldPool* m_pPool;

    TArray<BYTE> m_byGaugeFieldIds;
    TArray<BYTE> m_byBosonFieldIds;
};

class CLGAPI CFieldPool
{
public:
    CFieldPool(CField* pOrignal, UINT uiCount)
        : m_pOrignal(pOrignal)
    {
        for (UINT i = 0; i < uiCount; ++i)
        {
            m_pPool.AddItem(CreateNew());
        }
    }

    virtual ~CFieldPool()
    {
        for (INT i = 0; i < m_pPool.Num(); ++i)
        {
            appSafeDelete(m_pPool[i]);
        }
    }

    CField* GetOne()
    {
        if (m_pPool.Num() > 0)
        {
            return m_pPool.Pop();
        }
        appGeneral(_T("Warning: Field Pool Out Number!!!\n"));
        CField* newOne = CreateNew();
        return newOne;
    }

    void Return(CField* pField)
    {
        assert(NULL != pField 
            && pField->m_byFieldId == m_pOrignal->m_byFieldId 
            && pField->m_pPool == this
            && pField != m_pOrignal);
        m_pPool.PushBack(pField);
    }

    void ClearAll()
    {
        for (INT i = 0; i < m_pPool.Num(); ++i)
        {
            appSafeDelete(m_pPool[i]);
        }
        m_pPool.RemoveAll();
    }

    void ReCopyAll() const
    {
        for (INT i = 0; i < m_pAll.Num(); ++i)
        {
            m_pOrignal->CopyTo(m_pAll[i]);
        }
    }

    CField* m_pOrignal;
    TArray<CField*> m_pPool;
    TArray<CField*> m_pAll;

protected:

    CField * CreateNew()
    {
        CField* pNew = m_pOrignal->GetCopy();
        pNew->m_pPool = this;
        m_pAll.AddItem(pNew);
        return pNew;
    }
};

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