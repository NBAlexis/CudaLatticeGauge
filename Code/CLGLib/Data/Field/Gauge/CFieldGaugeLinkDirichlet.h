//=============================================================================
// FILENAME : CFieldGaugeLinkDirichlet.h
// 
// DESCRIPTION:
// There is simplifications for periodic boundary condition
// which is invalid for Dirichlet.
// This is only for Dirichlet.
//
// REVISION:
//  [07/06/2019 nbale]
//=============================================================================
#include "CFieldGaugeLink.h"

#define __DEFINE_GAUGE_LINKD(N) \
__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU##N##D) \
class CLGAPI CFieldGaugeSU##N##D : public CFieldGaugeLinkD<deviceSU##N, N> \
{ \
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldGaugeSU##N##D) \
public: \
    EFieldType GetFieldType() const override { return EFT_GaugeSU##N; } \
};

#ifndef _CFIELDGAUGELINKDIRICHLET_H_
#define _CFIELDGAUGELINKDIRICHLET_H_

__BEGIN_NAMESPACE

template<typename deviceGauge, INT matrixN>
class __DLL_EXPORT CFieldGaugeLinkD : public CFieldGaugeLink<deviceGauge, matrixN>
{
public:

#pragma region HMC

    void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override
    {
        if (NULL == pForce || this->GetFieldType() != pForce->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: force field is not SU3");
            return;
        }
        if (NULL != pStaple && this->GetFieldType() != pStaple->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: stape field is not SU3");
            return;
        }

        CFieldGaugeLink<deviceGauge, matrixN>* pForceSU3 = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pForce);
        CFieldGaugeLink<deviceGauge, matrixN>* pStableSU3 = NULL == pStaple ? NULL : dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pStaple);

        CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceAndStaple_D(
            this->m_pDeviceData,
            this->m_byFieldId,
            pForceSU3->m_pDeviceData,
            NULL == pStableSU3 ? NULL : pStableSU3->m_pDeviceData,
            betaOverN);
    }

    void CalculateOnlyStaple(CFieldGauge* pStaple) const override
    {
        if (NULL == pStaple || this->GetFieldType() != pStaple->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: stable field is not SU3");
            return;
        }
        CFieldGaugeLink<deviceGauge, matrixN>* pStableSU3 = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pStaple);
        CFieldGaugeKernel<deviceGauge, matrixN>::CalculateOnlyStaple_D(this->m_pDeviceData, this->m_byFieldId, pStableSU3->m_pDeviceData);
    }

    DOUBLE CalculatePlaqutteEnergy(DOUBLE betaOverN) const override
    {
        return CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergy_D(this->m_pDeviceData, this->m_byFieldId, betaOverN);
    }

    DOUBLE CalculateKinematicEnergy() const override
    {
        return CCommonKernelLink<deviceGauge>::CalcKineticEnery(this->m_pDeviceData, this->m_byFieldId);
    }

    DOUBLE CalculatePlaqutteEnergyUsingStable(DOUBLE betaOverN, const CFieldGauge* pStaple) const override
    {
        return CalculatePlaqutteEnergy(betaOverN);
    }

#pragma endregion

#pragma region BLAS

    void FixBoundary() override
    {
        appDetailed(_T("CFieldGaugeLinkD<deviceGauge, matrixN>::FixBoundary()\n"));
        CCommonKernelLink<deviceGauge>::FixBoundary(this->m_pDeviceData, this->m_byFieldId);
    }

#pragma endregion

    CCString GetInfos(const CCString& tab) const override
    {
        CCString sRet = CFieldGaugeLink<deviceGauge, matrixN>::GetInfos(tab);
        SSmallInt4 boundary = appGetLattice()->m_pIndex->GetBoudanryCondition()->GetFieldBC(this->m_byFieldId);
        sRet = sRet + tab + appToString(boundary) + _T("\n");

        return sRet;
    }
};

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeU1D)

class CLGAPI CFieldGaugeU1D : public CFieldGaugeLinkD<CLGComplex, 1>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldGaugeU1D)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeU1; }

    void InitialWithByteCompressed(const CCString& sFileName) override;
    CCString SaveToCompressedFile(const CCString& fileName) const override;
};

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU2D)

class CLGAPI CFieldGaugeSU2D : public CFieldGaugeLinkD<deviceSU2, 2>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldGaugeSU2D)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU2; }

    void InitialWithByteCompressed(const CCString& sFileName) override;
    CCString SaveToCompressedFile(const CCString& fileName) const override;
};

__DEFINE_GAUGE_LINKD(4)
__DEFINE_GAUGE_LINKD(5)
__DEFINE_GAUGE_LINKD(6)
__DEFINE_GAUGE_LINKD(7)
__DEFINE_GAUGE_LINKD(8)


__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGELINKDIRICHLET_H_

//=============================================================================
// END OF FILE
//=============================================================================