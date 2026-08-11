//=============================================================================
// FILENAME : CFieldGaugeLinkDirichlet.h
// 
// DESCRIPTION:
// There is simplifications for periodic boundary condition
// which is invalid for Dirichlet.
// This is only for Dirichlet.
//
// REVISION:
//  [mm/dd/yy]
//  [07/06/2019 nbale]
//=============================================================================
#include "CFieldGaugeLink.h"

#ifndef _CFIELDGAUGELINKDIRICHLET_H_
#define _CFIELDGAUGELINKDIRICHLET_H_

__BEGIN_NAMESPACE

template<class CFeildG>
class __DLL_EXPORT CFieldGaugeLinkD : public CFeildG
{
public:
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

        CFeildG* pForceSU3 = dynamic_cast<CFeildG*>(pForce);
        CFeildG* pStapleSU3 = NULL == pStaple ? NULL : dynamic_cast<CFeildG*>(pStaple);

        CFeildG::_GaugeKernel::CalculateForceAndStaple_D(
            this->m_pDeviceData,
            this->m_byFieldId,
            pForceSU3->m_pDeviceData,
            NULL == pStapleSU3 ? NULL : pStapleSU3->m_pDeviceData,
            betaOverN);
    }

    void CalculateForceAnisotropy(CFieldGauge* pForce, DOUBLE betaOverN, DOUBLE xi) const override
    {
        if (NULL == pForce || this->GetFieldType() != pForce->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: force field is not SU3");
            return;
        }
        CFeildG* pForceSU3 = dynamic_cast<CFeildG*>(pForce);
        CFeildG::_GaugeKernel::CalculateForceAnisotropy_D(
            this->m_pDeviceData,
            this->m_byFieldId,
            pForceSU3->m_pDeviceData,
            betaOverN,
            xi);
    }

    void CalculateForceAndStapleClover(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override
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

        CFeildG* pForceSU3 = dynamic_cast<CFeildG*>(pForce);
        CFeildG* pStapleSU3 = NULL == pStaple ? NULL : dynamic_cast<CFeildG*>(pStaple);

        CFeildG::_GaugeKernel::CalculateForceAndStapleClover_D(
            this->m_pDeviceData,
            this->m_byFieldId,
            pForceSU3->m_pDeviceData,
            NULL == pStapleSU3 ? NULL : pStapleSU3->m_pDeviceData,
            betaOverN);
    }

    void CalculateForceCloverAnisotropy(CFieldGauge* pForce, DOUBLE betaOverN, DOUBLE xi) const override
    {
        if (NULL == pForce || this->GetFieldType() != pForce->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: force field is not SU3");
            return;
        }
        CFeildG* pForceSU3 = dynamic_cast<CFeildG*>(pForce);
        CFeildG::_GaugeKernel::CalculateForceCloverAnisotropy_D(
            this->m_pDeviceData,
            this->m_byFieldId,
            pForceSU3->m_pDeviceData,
            betaOverN,
            xi);
    }

    void CalculateOnlyStaple(CFieldGauge* pStaple) const override
    {
        if (NULL == pStaple || this->GetFieldType() != pStaple->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: staple field is not SU3");
            return;
        }
        CFeildG* pStapleSU3 = dynamic_cast<CFeildG*>(pStaple);
        CFeildG::_GaugeKernel::CalculateOnlyStaple_D(this->m_pDeviceData, this->m_byFieldId, pStapleSU3->m_pDeviceData);
    }

    DOUBLE CalculatePlaqutteEnergy(DOUBLE betaOverN) const override
    {
        return CFeildG::_GaugeKernel::CalculatePlaqutteEnergy_D(this->m_pDeviceData, this->m_byFieldId, betaOverN);
    }

    DOUBLE CalculatePlaqutteEnergyOriginal(DOUBLE betaOverN) const override
    {
        _RECORD(CFieldGaugeLink::CalculatePlaqutteEnergyOriginal);
        appWarning("Dirichlet boundary gauge not fully support u0 yet! the result will be larger than this.\n");
        return CFeildG::_GaugeKernel::CalculatePlaqutteEnergy_D(this->m_pDeviceData, this->m_byFieldId, betaOverN);
    }

    DOUBLE CalculatePlaqutteEnergyAnisotropy(DOUBLE betaOverN, DOUBLE xi) const override
    {
        return CFeildG::_GaugeKernel::CalculatePlaqutteEnergyAnisotropy_D(this->m_pDeviceData, this->m_byFieldId, betaOverN, xi);
    }

    DOUBLE CalculatePlaqutteEnergyUsingStaple(DOUBLE betaOverN, const CFieldGauge* pStaple) const override
    {
        return CalculatePlaqutteEnergy(betaOverN);
    }

    /**
    * for test
    */
    //DOUBLE CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const override
    //{
    //    return CalculatePlaqutteEnergy(betaOverN);
    //}

    DOUBLE CalculateKinematicEnergy() const override 
    { 
        if (abs(_HC_GaugeMomentumFactor - F(1.0)) > _CLG_FLT_EPSILON)
        {
            return CFeildG::_LinkKernel::CalcKineticEnery(this->m_pDeviceData, this->m_byFieldId) / _HC_GaugeMomentumFactor;
        }
        return CFeildG::_LinkKernel::CalcKineticEnery(this->m_pDeviceData, this->m_byFieldId);
    } 

    void InitialField(EFieldInitialType eInitialType) override
    {
        CFeildG::_LinkKernel::InitialBufferD(this->m_pDeviceData, this->m_byFieldId, eInitialType);
        if ((EFIT_RandomGenerator == eInitialType || EFIT_RandomGaussian == eInitialType) && abs(_HC_GaugeMomentumFactor - F(1.0)) > _CLG_FLT_EPSILON)
        {
            CFeildG::ScalarMultply(static_cast<Real>(sqrt(_HC_GaugeMomentumFactor)));
        }
    }

    void FixBoundary(EFixBoundary eType) override
    {
        //appDetailed(_T("CFieldGaugeLinkD<deviceGauge, matrixN>::FixBoundary()\n"));
        if (EFB_Field == eType)
        {
            CFeildG::_LinkKernel::FixBoundary(this->m_pDeviceData, this->m_byFieldId);
        }
        else
        {
            CFeildG::_LinkKernel::FixBoundaryZero(this->m_pDeviceData, this->m_byFieldId);
        }
    }

    CCString GetInfos(const CCString& tab) const override
    {
        CCString sRet = CFeildG::GetInfos(tab);
        SSmallInt4 boundary = appGetLattice()->m_pIndex->GetBoudanryCondition()->GetFieldBC(this->m_byFieldId);
        sRet = sRet + tab + appToString(boundary) + _T("\n");
        return sRet;
    }

    UBOOL IsDirichlet() const override
    {
        return TRUE;
    }
};

#define _DEFINE_Gauge_Dirichlet(N) \
__CLG_REGISTER_HELPER_HEADER(CFieldGauge##N##D) \
class CLGAPI CFieldGauge##N##D : public CFieldGaugeLinkD<CFieldGauge##N> \
{ \
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldGauge##N##D) \
};

_DEFINE_Gauge_Dirichlet(U1)
_DEFINE_Gauge_Dirichlet(SU2)
_DEFINE_Gauge_Dirichlet(SU3)

#if _CLG_SU4_GAUGE
_DEFINE_Gauge_Dirichlet(SU4)
#endif
#if _CLG_SU5_GAUGE
_DEFINE_Gauge_Dirichlet(SU5)
#endif
#if _CLG_SU6_GAUGE
_DEFINE_Gauge_Dirichlet(SU6)
#endif
#if _CLG_SU7_GAUGE
_DEFINE_Gauge_Dirichlet(SU7)
#endif
#if _CLG_SU8_GAUGE
_DEFINE_Gauge_Dirichlet(SU8)
#endif

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGELINKDIRICHLET_H_

//=============================================================================
// END OF FILE
//=============================================================================