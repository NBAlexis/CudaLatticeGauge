//=============================================================================
// FILENAME : CField.cpp
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
//
// REVISION:
//  [mm/dd/yy]
//  [12/7/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "WilsonDirac/CFieldFermionWilsonSquareSU3.h"
#include "Staggered/CFieldFermionKST.h"

__BEGIN_NAMESPACE

CField::CField()
    : CBase()
    , m_pOwner(NULL)
    // Default to field id 1, the canonical (always-registered) field. An
    // ad-hoc field created via appCreate is never assigned an id, and kernels
    // that index the per-field index tables (e.g. StrictExp) would otherwise
    // dereference m_pDeviceIndexLinkToSIndex[garbage] -> illegal address.
    , m_byFieldId(1)
    , m_bDynamic(TRUE)
    , m_fLength(F(1.0))
    // , m_pPool(NULL)
    , m_pClass(NULL)
{
    //I cannot call virtual function in this construction function ...
    //m_pClass = GetClass();
    //printf("%s\n", GetClass()->GetName());
}

void CField::Return()
{
    //appAssert(NULL != m_pPool);
    //m_pPool->Return(this);
    appGetFieldPool()->Return(this);
}

CCString CField::SaveToFile(const CCString& fileName, EFieldFileType eType) const
{
    if (EFFT_CLGBinCompressed == eType)
    {
        return SaveToCompressedFile(fileName);
    }

    UINT uiSize = 0;
    BYTE* byToSave = NULL;
    switch (eType)
    {
    case EFFT_CLGBin:
        byToSave = CopyDataOut(uiSize);
        break;
    case EFFT_CLGBinFloat:
        byToSave = CopyDataOutFloat(uiSize);
        break;
    case EFFT_CLGBinDouble:
        byToSave = CopyDataOutDouble(uiSize);
        break;
    default:
        appCrucial(_T("Save for this type not implemented: CFieldGaugeSU3 : %s\n"), __ENUM_TO_STRING(EFieldFileType, eType).c_str());
        return _T("Not supported");
    }

#if _CLG_MULTI_GPU
    //Multi-GPU: each rank holds a sub-lattice. Gather to rank 0 in global-site
    //order so the file on disk is identical to the single-GPU layout (this is the
    //1-vs-N validation harness, Docs/MultiGPU-Plan.md section 8.5). Non-root ranks
    //write nothing.
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        const UINT uiLocalVolume = _HC_Volume;
        const UINT uiBytesPerSite = (0 == uiLocalVolume) ? 0 : (uiSize / uiLocalVolume);
        UINT uiGlobalSize = 0;
        BYTE* byGlobal = appGetComm()->GatherFieldToRoot(byToSave, uiBytesPerSite, uiGlobalSize);
        free(byToSave);

        if (!appGetComm()->IsRoot())
        {
            //Freed on non-root (GatherFieldToRoot returns NULL there). Wait for
            //root to finish writing so a later load by this rank cannot read a
            //missing / half-written file (cross-rank read-ahead race).
            appGetComm()->Barrier();
            return _T("");
        }
        byToSave = byGlobal;
        uiSize = uiGlobalSize;
    }
#endif

    appGetFileSystem()->WriteAllBytes(fileName.c_str(), byToSave, uiSize);
    CCString MD5 = CLGMD5Hash(byToSave, uiSize);
    free(byToSave);
#if _CLG_MULTI_GPU
    //Pair with the non-root barrier above: release them once the file is on disk.
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        appGetComm()->Barrier();
    }
#endif
    return MD5;
}

CCString CField::GetInfos(const CCString& tab) const
{
    CCString sRet = CBase::GetInfos(tab);
    sRet = sRet + tab + _T("FieldId : ") + appToString(m_byFieldId) + _T("\n");
    sRet = sRet + tab + _T("Dynamic : ") + appToString(m_bDynamic) + _T("\n");
    sRet = sRet + tab + _T("GaugeFields : ") + appToString(m_byGaugeFieldIds) + _T("\n");
    sRet = sRet + tab + _T("BosonFields : ") + appToString(m_byBosonFieldIds) + _T("\n");
    return sRet;
}

//void CField::UpdatePooledParamters() const
//{
//    appGetLattice()->ReCopyPooled(m_byFieldId);
//}


CField* CFieldPool::GetOne(const CField* pOrignal)
{
    appAssert(NULL != pOrignal && NULL != pOrignal->m_pClass);

    if (m_pPool.Exist(pOrignal->m_pClass))
    {
        TArray<CPooledFields>& fields = m_pPool[pOrignal->m_pClass];
        for (INT i = 0; i < fields.Num(); ++i)
        {
            if (!fields[i].m_bInUse)
            {
                fields[i].m_bInUse = TRUE;
                pOrignal->CopyParamTo(fields[i].m_pField);
                return fields[i].m_pField;
            }
        }

        CField* newOne = dynamic_cast<CField*>(pOrignal->m_pClass->Create());
        pOrignal->CopyParamTo(newOne);
        CPooledFields newOneRecord;
        newOneRecord.m_bInUse = TRUE;
        newOneRecord.m_pField = newOne;
        fields.AddItem(newOneRecord);

        return newOne;
    }

    TArray<CPooledFields> fields;
    CField* newOne = dynamic_cast<CField*>(pOrignal->m_pClass->Create());
    pOrignal->CopyParamTo(newOne);
    CPooledFields newOneRecord;
    newOneRecord.m_bInUse = TRUE;
    newOneRecord.m_pField = newOne;
    fields.AddItem(newOneRecord);
    m_pPool[pOrignal->m_pClass] = fields;

    return newOne;
}

void CFieldPool::Return(CField* pField)
{
    appAssert(NULL != pField && NULL != pField->m_pClass && m_pPool.Exist(pField->m_pClass));
    TArray<CPooledFields>& fields = m_pPool[pField->m_pClass];
    for (INT i = 0; i < fields.Num(); ++i)
    {
        if (fields[i].m_pField == pField)
        {
            appAssert(fields[i].m_bInUse);
            fields[i].m_bInUse = FALSE;
            return;
        }
    }
    appAssert(FALSE);
}


__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================