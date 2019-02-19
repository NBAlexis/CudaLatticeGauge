//=============================================================================
// FILENAME : CLatticeData.h
// 
// DESCRIPTION:
// This is the class for the lattce data
// NOTE:: We only have 4D case, 3D = 1xLxLxL, and 2D= 1x1xLxL
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CLATTICEDATA_H_
#define _CLATTICEDATA_H_

__BEGIN_NAMESPACE

class CLGAPI CLatticeData
{
public:
    /**
    * Initial with CommonData
    */
    CLatticeData();
    ~CLatticeData();

    void CreateFermionSolver(const CCString& sSolver, const CParameters& param, const class CField* pFermionField);
    void OnUpdatorConfigurationAccepted();
    void OnUpdatorFinished(UBOOL bMeasured);
    void GetPlaquetteLengthCount(BYTE& plaqLength, BYTE& countPerSite, BYTE& countPerLink);
    void CreateFieldPool(BYTE byFieldId, UINT uiCount);
    void SetFieldBoundaryCondition(BYTE byFieldId, const SBoundCondition& bc);
    CCString GetInfos(const CCString& sTab) const;

    class CRandom* m_pRandom;
    UINT m_uiRandomType;
    UINT m_uiRandomSeed;

    class CFieldGauge* m_pGaugeField;
    TArray<class CField*> m_pOtherFields; //for delete
    THashMap<BYTE, class CField*> m_pFieldMap;
    TArray<class CFieldPool*> m_pFieldPools;
    THashMap<BYTE, class CFieldPool*> m_pFieldPoolMap;
    class CFieldCache* m_pFieldCache;

    /**
    * \note we assume the first action is gauge action
    */
    TArray<class CAction*> m_pActionList;
    THashMap<BYTE, class CAction*> m_pActionMap;

    class CUpdator* m_pUpdator;
    class CMeasurementManager* m_pMeasurements;

    //this is a device copy
    //see:
    //https://stackoverflow.com/questions/53781421/cuda-the-member-field-with-device-ptr-and-device-member-function-to-visit-it-i
    class CRandom* m_pDeviceRandom;
    class CIndex* m_pDeviceIndex;
    class CIndexCache* m_pIndexCache;
    UINT m_uiIndexType;
    UINT m_uiBoundaryConditionType;
    class CSLASolver* m_pFermionSolver;

    class CField* GetFieldById(BYTE byId) const { return m_pFieldMap.GetAt(byId); }
    class CAction* GetActionById(BYTE byId) const { return m_pActionMap.GetAt(byId); }
    class CField* GetPooledFieldById(BYTE byId);
};

inline class CSLASolver* appGetFermionSolver();

__END_NAMESPACE

#endif //#ifndef _CLATTICEDATA_H_

//=============================================================================
// END OF FILE
//=============================================================================