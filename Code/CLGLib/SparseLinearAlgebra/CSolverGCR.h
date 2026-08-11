//=============================================================================
// FILENAME : CSolverGCR.h
// 
// DESCRIPTION:
//
// Both GCR and ORTHOMIN has oscillation problem, ORTHODIR is good
// So the implementation is in fact ORTHODIR
//
// REVISION:
//  [mm/dd/yy]
//  [02/16/2019 nbale]
//=============================================================================

#ifndef _CSOLVERGCR_H_
#define _CSOLVERGCR_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CSLASolverGCR)

class CLGAPI CSLASolverGCR : public CSLASolver
{
    __CLGDECLARE_CLASS(CSLASolverGCR)

public:

    enum { _kMaxStep = 100, };

    CSLASolverGCR();
    ~CSLASolverGCR();

    void Configurate(const CParameters& param) override;
    void AllocateBuffers(const CField* pField) override;
    virtual void ReleaseBuffers();
    UBOOL Solve(CField* pFieldX, const CField* pFieldB, 
        INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
        EFieldOperator uiM, ESolverPhase ePhase = ESP_Once, const CField* pStart = NULL) override;

protected:

    void GenerateNextP(CField* pNextP, CField* pAAP, CField* pR, CField* pAPJ,
        const TArray<CField*>& pP, const TArray<CField*>& pAP, const TArray<Real>& lengthAP,
        UINT uiCurrentStep, UINT uiCurrentIndex, UINT uiNextIndex,
        INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
        EFieldOperator uiM) const;

    UINT m_uiReStart;
    UINT m_uiMaxDim;
    UINT m_uiIterateNumber;
    UINT m_uiCheckError;
    Real m_fAccuracy;
};

__END_NAMESPACE

#endif //#ifndef _CSOLVERGCR_H_

//=============================================================================
// END OF FILE
//=============================================================================
