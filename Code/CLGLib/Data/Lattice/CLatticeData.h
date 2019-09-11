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
    void OnUpdatorConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) const;
    void OnUpdatorFinished(UBOOL bMeasured, UBOOL bReport) const;
    //void GetPlaquetteLengthCount(BYTE& plaqLength, BYTE& countPerSite, BYTE& countPerLink);
    void CreateFieldPool(BYTE byFieldId, UINT uiCount);
    void SetFieldBoundaryCondition(BYTE byFieldId, const SBoundCondition& bc) const;
    void FixAllFieldBoundary() const;

    CCString GetInfos(const CCString& sTab) const;

    class CRandom* m_pRandom;
    UINT m_uiRandomType;
    UINT m_uiRandomSeed;

    class CFieldGauge* m_pGaugeField;
    /**
    * Now we support only 2x2x2x2 x 4 field. 
    * U_mu(nu=0) is the nu=0 boundary, U_mu(nu=1) is the nu=N_nu boundary
    */
    class CFieldGauge* m_pGaugeBoundary;

    TArray<class CField*> m_pOtherFields; //only for delete
    THashMap<BYTE, class CField*> m_pFieldMap;
    TArray<class CFieldBoundary*> m_pAllBoundaryFields; //only for delete
    THashMap<BYTE, class CFieldBoundary*> m_pBoundaryFieldMap;

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
    class CIndex* m_pIndex;
    class CIndexData* m_pIndexCache;

    class CSLASolver* m_pFermionSolver;

    class CGaugeSmearing* m_pGaugeSmearing;

    class CField* GetFieldById(BYTE byId) const { return m_pFieldMap.Exist(byId) ? m_pFieldMap.GetAt(byId) : NULL; }
    class CFieldBoundary* GetBoundaryFieldById(BYTE byId) const { return m_pBoundaryFieldMap.Exist(byId) ? m_pBoundaryFieldMap.GetAt(byId) : NULL; }
    class CAction* GetActionById(BYTE byId) const { return m_pActionMap.Exist(byId) ? m_pActionMap.GetAt(byId) : NULL; }
    class CField* GetPooledFieldById(BYTE byId);
};

inline class CSLASolver* appGetFermionSolver();

inline class CGaugeSmearing* appGetGaugeSmearing();


__END_NAMESPACE

#endif //#ifndef _CLATTICEDATA_H_

//=============================================================================
// END OF FILE
//=============================================================================