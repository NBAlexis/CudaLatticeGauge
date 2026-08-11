//=============================================================================
// FILENAME : CFieldGaugeSU3TreeImproved.h
// 
// DESCRIPTION:
// Note:
// This was beta * [(5/3) plaq - (1/12) rect] with beta=6/g^2
// Use beta = 10/g^2 instead of 6/g^2,
// this is then beta * [plaq + (-1/20) rect]
//
// Note:
// Whenever staple is mentioned, we mean the staple in the sense of plaquettes only
// This is also used to define a 'fat' link
// So the rectangle terms are not included
//
// Note2: Eq.~(12) of hep-lat/0203008
// We use 1/12 x sum_{mu>nu}(6-Retr[rect1]-Retr[rect2]) = 1/12 x [36-sum_{mu>nu} (Retr[rect1]+Retr[rect2]) ]
// In some litertures, it  is 1/6 x sum_{mu>nu}[ 3 - (Retr[rect1]+Retr[rect2])/2 ]
//
// REVISION:
//  [mm/dd/yy]
//  [10/03/2020 nbale]
//=============================================================================
#pragma once

#ifndef _CFIELDGAUGE_SU3_TREEIMPROVED_H_
#define _CFIELDGAUGE_SU3_TREEIMPROVED_H_

#include "CFieldGaugeLink.h"
#include "CFieldGaugeLinkDirichlet.h"


__BEGIN_NAMESPACE

template<class CFeildG>
class __DLL_EXPORT CFieldGaugeTreeImproved : public CFeildG
{
public:
    CFieldGaugeTreeImproved()
        : CFeildG()
        , m_fRectOverPlaq(-0.05)
    {

    }

    void InitialOtherParameters(CParameters& param) override
    {
        CFeildG::InitialOtherParameters(param);
        DOUBLE fValue = -0.05;
        if (param.FetchValueDOUBLE(_T("RectOverPlaq"), fValue))
        {
            m_fRectOverPlaq = fValue;
        }
    }

    void CopyParamTo(CField* U) const override
    {
        CFeildG::CopyParamTo(U);
        CFieldGaugeTreeImproved<CFeildG>* pFieldGauge = dynamic_cast<CFieldGaugeTreeImproved<CFeildG>*>(U);
        if (NULL == pFieldGauge)
        {
            return;
        }
        pFieldGauge->m_fRectOverPlaq = m_fRectOverPlaq;
    }

    void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override
    {
        _RECORD(CFieldGaugeTreeImproved::CalculateForceAndStaple);
        CFeildG::CalculateForceAndStaple(pForce, pStaple, betaOverN);
        CFeildG* pForceSU3 = dynamic_cast<CFeildG*>(pForce);

        DOUBLE coef = m_fRectOverPlaq * betaOverN * -0.5;

        if (this->IsDirichlet())
        {
            CFeildG::_GaugeKernel::CalculateForceRectangular_D(
                this->m_pDeviceData,
                this->m_byFieldId,
                pForceSU3->m_pDeviceData,
                coef,
                appGetLattice()->m_pIndex->GetBoudanryCondition()->GetFieldBC(this->m_byFieldId)
                );
        }
        else
        {
            CFeildG::_GaugeKernel::CalculateForceRectangular(
                this->m_pDeviceData,
                this->m_byFieldId,
                pForceSU3->m_pDeviceData,
                coef);
        }
    }

    void CalculateForceAndStapleClover(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override
    {
        _RECORD(CFieldGaugeTreeImproved::CalculateForceAndStapleClover);
        CFeildG::CalculateForceAndStapleClover(pForce, pStaple, betaOverN);
        CFeildG* pForceSU3 = dynamic_cast<CFeildG*>(pForce);

        DOUBLE coef = m_fRectOverPlaq * betaOverN * -0.5;

        if (this->IsDirichlet())
        {
            CFeildG::_GaugeKernel::CalculateForceRectangularClover_D(
                this->m_pDeviceData,
                this->m_byFieldId,
                pForceSU3->m_pDeviceData,
                coef,
                appGetLattice()->m_pIndex->GetBoudanryCondition()->GetFieldBC(this->m_byFieldId));
        }
        else
        {
            CFeildG::_GaugeKernel::CalculateForceRectangular(
                this->m_pDeviceData,
                this->m_byFieldId,
                pForceSU3->m_pDeviceData,
                coef);
        }
    }

    DOUBLE CalculatePlaqutteEnergy(DOUBLE betaOverN) const override
    {
        _RECORD(CFieldGaugeTreeImproved::CalculatePlaqutteEnergy);
        DOUBLE fPlaqTerm = CFeildG::CalculatePlaqutteEnergy(betaOverN);
        //appParanoiac(_T("fPlaqTerm=%f\n"), fPlaqTerm);
        fPlaqTerm += CFeildG::_GaugeKernel::CalculateRectangularEnergy(
            this->m_pDeviceData,
            this->m_byFieldId,
            m_fRectOverPlaq * betaOverN);

        return fPlaqTerm;
    }

    DOUBLE CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const override
    {
        _RECORD(CFieldGaugeTreeImproved::CalculatePlaqutteEnergyUseClover);
        DOUBLE fPlaqTerm = CFeildG::CalculatePlaqutteEnergyUseClover(betaOverN);

        fPlaqTerm += CFeildG::_GaugeKernel::CalculateRectangularEnergyUseClover(
            this->m_pDeviceData,
            this->m_byFieldId,
            m_fRectOverPlaq * betaOverN);

        return fPlaqTerm;
    }

    DOUBLE CalculatePlaqutteEnergyUsingStaple(DOUBLE betaOverN, const CFieldGauge* pStaple) const override
    {
        //appGeneral(_T("Calculate energy using stable is not supported...\n"));
        //We don't known whether to use clover or not, so use normal
        DOUBLE fPlaqTerm = CFeildG::CalculatePlaqutteEnergyUsingStaple(betaOverN, pStaple);
        fPlaqTerm += CFeildG::_GaugeKernel::CalculateRectangularEnergy(
            this->m_pDeviceData,
            this->m_byFieldId,
            m_fRectOverPlaq * betaOverN);

        return fPlaqTerm;
    }

    void CalculateForceAnisotropy(CFieldGauge* pForce, DOUBLE betaOverN, DOUBLE xi) const override
    {
        _RECORD(CFieldGaugeTreeImproved::CalculateForceAnisotropy);
        CFeildG::CalculateForceAnisotropy(pForce, betaOverN, xi);
        CFeildG* pForceSU3 = dynamic_cast<CFeildG*>(pForce);

        DOUBLE coef = m_fRectOverPlaq * betaOverN * -0.5;

        if (this->IsDirichlet())
        {
            CFeildG::_GaugeKernel::CalculateForceRectangularAnisotropy_D(
                this->m_pDeviceData,
                this->m_byFieldId,
                pForceSU3->m_pDeviceData,
                coef,
                xi,
                appGetLattice()->m_pIndex->GetBoudanryCondition()->GetFieldBC(this->m_byFieldId)
                );
        }
        else
        {
            CFeildG::_GaugeKernel::CalculateForceRectangularAnisotropy(
                this->m_pDeviceData,
                this->m_byFieldId,
                pForceSU3->m_pDeviceData,
                coef,
                xi);
        }
    }

    void CalculateForceCloverAnisotropy(CFieldGauge* pForce, DOUBLE betaOverN, DOUBLE xi) const override
    {
        _RECORD(CFieldGaugeTreeImproved::CalculateForceCloverAnisotropy);
        CFeildG::CalculateForceCloverAnisotropy(pForce, betaOverN, xi);
        CFeildG* pForceSU3 = dynamic_cast<CFeildG*>(pForce);

        DOUBLE coef = m_fRectOverPlaq * betaOverN * -0.5;

        if (this->IsDirichlet())
        {
            CFeildG::_GaugeKernel::CalculateForceRectangularCloverAnisotropy_D(
                this->m_pDeviceData,
                this->m_byFieldId,
                pForceSU3->m_pDeviceData,
                coef,
                xi,
                appGetLattice()->m_pIndex->GetBoudanryCondition()->GetFieldBC(this->m_byFieldId));
        }
        else
        {
            CFeildG::_GaugeKernel::CalculateForceRectangularAnisotropy(
                this->m_pDeviceData,
                this->m_byFieldId,
                pForceSU3->m_pDeviceData,
                coef,
                xi);
        }
    }

    DOUBLE CalculatePlaqutteEnergyAnisotropy(DOUBLE betaOverN, DOUBLE xi) const override
    {
        _RECORD(CFieldGaugeTreeImproved::CalculatePlaqutteEnergyAnisotropy);
        DOUBLE fPlaqTerm = CFeildG::CalculatePlaqutteEnergyAnisotropy(betaOverN, xi);

        fPlaqTerm += CFeildG::_GaugeKernel::CalculateRectangularEnergyAnisotropy(
            this->m_pDeviceData,
            this->m_byFieldId,
            m_fRectOverPlaq * betaOverN,
            xi);

        return fPlaqTerm;
    }

    DOUBLE CalculatePlaqutteEnergyUseCloverAnisotropy(DOUBLE betaOverN, DOUBLE xi) const override
    {
        _RECORD(CFieldGaugeTreeImproved::CalculatePlaqutteEnergyUseCloverAnisotropy);
        DOUBLE fPlaqTerm = CFeildG::CalculatePlaqutteEnergyUseCloverAnisotropy(betaOverN, xi);

        fPlaqTerm += CFeildG::_GaugeKernel::CalculateRectangularEnergyUseCloverAnisotropy(
            this->m_pDeviceData,
            this->m_byFieldId,
            m_fRectOverPlaq * betaOverN,
            xi);

        return fPlaqTerm;
    }

    CCString GetInfos(const CCString& tab) const override
    {
        CCString sRet = CFeildG::GetInfos(tab);
        sRet = sRet + _T("Cr:") + appToString(m_fRectOverPlaq) + _T("\n");
        return sRet;
    }

    DOUBLE m_fRectOverPlaq;

};

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU3TreeImproved)
class CLGAPI CFieldGaugeSU3TreeImproved : public CFieldGaugeTreeImproved<CFieldGaugeSU3>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldGaugeSU3TreeImproved)
};

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU3TreeImprovedD)
class CLGAPI CFieldGaugeSU3TreeImprovedD : public CFieldGaugeTreeImproved<CFieldGaugeSU3D>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldGaugeSU3TreeImprovedD)
};

__END_NAMESPACE


#endif //#ifndef _CFIELDGAUGE_SU3_TREEIMPROVED_H_

//=============================================================================
// END OF FILE
//=============================================================================