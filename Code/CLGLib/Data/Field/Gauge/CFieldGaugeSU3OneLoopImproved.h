//=============================================================================
// FILENAME : CFieldGaugeSU3OneLoopImproved.h
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

#ifndef _CFIELDGAUGE_SU3_OneLoopIMPROVED_H_
#define _CFIELDGAUGE_SU3_OneLoopIMPROVED_H_

#include "CFieldGaugeLink.h"
#include "CFieldGaugeLinkDirichlet.h"


__BEGIN_NAMESPACE

template<class CFeildG>
class __DLL_EXPORT CFieldGaugeOneLoopImproved : public CFeildG
{
public:
    CFieldGaugeOneLoopImproved()
        : CFeildG()
        , m_fCr(-0.05)
        , m_fCt(0.0)
        , m_fU0(1.0)
        , m_iNf(2)
        , m_bDirectSetCoefficients(FALSE)
    {

    }

    void InitialOtherParameters(CParameters& param) override
    {
        CFeildG::InitialOtherParameters(param);
        if (param.FetchValueDOUBLE(_T("Cr"), m_fCr) && param.FetchValueDOUBLE(_T("Ct"), m_fCt))
        {
            m_bDirectSetCoefficients = TRUE;
        }
        else
        {
            m_bDirectSetCoefficients = FALSE;
            param.FetchValueDOUBLE(_T("U0"), m_fU0);
            param.FetchValueINT(_T("Nf"), m_iNf);

            //arXiv:1004.0342 Eq.~(A2)
            m_fCr = (- 1.0 / (20.0 * m_fU0 * m_fU0)) * (1.0 - (0.6264 - 1.1746*m_iNf) * log(m_fU0));
            // in pyquda, it was beta_1 = -beta_0 / (20 * u_0**2) * (1 + (0.4805 - 0.899 * n_flavor) * alpha_s), which is (- 1.0 / (20.0 * m_fU0 * m_fU0)) * (1.0 - (0.626387 - 1.17195*m_iNf) * log(m_fU0))
            // note the difference between 1.1746 and 1.17195

            m_fCt = (0.0433 - 0.0156*m_iNf) * log(m_fU0) / (m_fU0 * m_fU0);
            // in pyquda, it was beta_2 = -beta_0 / u_0**2 * (0.03325 - 0.0121 * n_flavor) * alpha_s, which is (0.0433452 - 0.0157737*m_iNf) * log(m_fU0) / (m_fU0 * m_fU0)
            // also note the difference in the second term
        }
    }

    void CopyParamTo(CField* U) const override
    {
        CFeildG::CopyParamTo(U);
        CFieldGaugeOneLoopImproved<CFeildG>* pFieldGauge = dynamic_cast<CFieldGaugeOneLoopImproved<CFeildG>*>(U);
        if (NULL == pFieldGauge)
        {
            return;
        }
        pFieldGauge->m_fCr = m_fCr;
        pFieldGauge->m_fCt = m_fCt;
        pFieldGauge->m_fU0 = m_fU0;
        pFieldGauge->m_iNf = m_iNf;
        pFieldGauge->m_bDirectSetCoefficients = m_bDirectSetCoefficients;
    }

    void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override
    {
        _RECORD(CFieldGaugeOneLoopImproved::CalculateForceAndStaple);
        CFeildG::CalculateForceAndStaple(pForce, pStaple, betaOverN);
        CFeildG* pForceSU3 = dynamic_cast<CFeildG*>(pForce);

        const DOUBLE coefCr = m_fCr * betaOverN * -0.5;
        const DOUBLE coefCt = m_fCt * betaOverN * -0.5;

        CFeildG::_GaugeKernel::CalculateForceRectangular(
            this->m_pDeviceData,
            this->m_byFieldId,
            pForceSU3->m_pDeviceData,
            coefCr);

        CFeildG::_GaugeKernel::CalculateForceTwistedLoop(
            this->m_pDeviceData,
            this->m_byFieldId,
            pForceSU3->m_pDeviceData,
            coefCt);
    }

    void CalculateForceAndStapleClover(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override
    {
        _RECORD(CFieldGaugeOneLoopImproved::CalculateForceAndStapleClover);
        CFeildG::CalculateForceAndStapleClover(pForce, pStaple, betaOverN);
        CFeildG* pForceSU3 = dynamic_cast<CFeildG*>(pForce);

        const DOUBLE coefCr = m_fCr * betaOverN * -0.5;
        const DOUBLE coefCt = m_fCt * betaOverN * -0.5;

        CFeildG::_GaugeKernel::CalculateForceRectangular(
            this->m_pDeviceData,
            this->m_byFieldId,
            pForceSU3->m_pDeviceData,
            coefCr);

        CFeildG::_GaugeKernel::CalculateForceTwistedLoop(
            this->m_pDeviceData,
            this->m_byFieldId,
            pForceSU3->m_pDeviceData,
            coefCt);
    }

    DOUBLE CalculatePlaqutteEnergy(DOUBLE betaOverN) const override
    {
        _RECORD(CFieldGaugeOneLoopImproved::CalculatePlaqutteEnergy);
        DOUBLE fPlaqTerm = CFeildG::CalculatePlaqutteEnergy(betaOverN);

        fPlaqTerm += CFeildG::_GaugeKernel::CalculateRectangularEnergy(
            this->m_pDeviceData,
            this->m_byFieldId,
            m_fCr * betaOverN);

        fPlaqTerm += CFeildG::_GaugeKernel::CalculateTwistedLoopEnergy(
            this->m_pDeviceData,
            this->m_byFieldId,
            m_fCt * betaOverN);

        return fPlaqTerm;
    }

    DOUBLE CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const override
    {
        _RECORD(CFieldGaugeOneLoopImproved::CalculatePlaqutteEnergyUseClover);
        DOUBLE fPlaqTerm = CFeildG::CalculatePlaqutteEnergyUseClover(betaOverN);

        fPlaqTerm += CFeildG::_GaugeKernel::CalculateRectangularEnergyUseClover(
            this->m_pDeviceData,
            this->m_byFieldId,
            m_fCr * betaOverN);

        fPlaqTerm += CFeildG::_GaugeKernel::CalculateTwistedLoopEnergy(
            this->m_pDeviceData,
            this->m_byFieldId,
            m_fCt * betaOverN);

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
            m_fCr * betaOverN);

        fPlaqTerm += CFeildG::_GaugeKernel::CalculateTwistedLoopEnergy(
            this->m_pDeviceData,
            this->m_byFieldId,
            m_fCt * betaOverN);

        return fPlaqTerm;
    }

    void CalculateForceAnisotropy(CFieldGauge* pForce, DOUBLE betaOverN, DOUBLE xi) const override
    {
        _RECORD(CFieldGaugeOneLoopImproved::CalculateForceAnisotropy);
        CFeildG::CalculateForceAnisotropy(pForce, betaOverN, xi);
        CFeildG* pForceSU3 = dynamic_cast<CFeildG*>(pForce);

        const DOUBLE coefCr = m_fCr * betaOverN * -0.5;
        const DOUBLE coefCt = m_fCt * betaOverN * -0.5;

        CFeildG::_GaugeKernel::CalculateForceRectangularAnisotropy(
            this->m_pDeviceData,
            this->m_byFieldId,
            pForceSU3->m_pDeviceData,
            coefCr,
            xi);

        CFeildG::_GaugeKernel::CalculateForceTwistedLoopAnisotropy(
            this->m_pDeviceData,
            this->m_byFieldId,
            pForceSU3->m_pDeviceData,
            coefCt,
            xi);
    }

    void CalculateForceCloverAnisotropy(CFieldGauge* pForce, DOUBLE betaOverN, DOUBLE xi) const override
    {
        _RECORD(CFieldGaugeOneLoopImproved::CalculateForceCloverAnisotropy);
        CFeildG::CalculateForceCloverAnisotropy(pForce, betaOverN, xi);
        CFeildG* pForceSU3 = dynamic_cast<CFeildG*>(pForce);

        const DOUBLE coefCr = m_fCr * betaOverN * -0.5;
        const DOUBLE coefCt = m_fCt * betaOverN * -0.5;

        CFeildG::_GaugeKernel::CalculateForceRectangularAnisotropy(
            this->m_pDeviceData,
            this->m_byFieldId,
            pForceSU3->m_pDeviceData,
            coefCr,
            xi);

        CFeildG::_GaugeKernel::CalculateForceTwistedLoopAnisotropy(
            this->m_pDeviceData,
            this->m_byFieldId,
            pForceSU3->m_pDeviceData,
            coefCt,
            xi);
    }

    DOUBLE CalculatePlaqutteEnergyAnisotropy(DOUBLE betaOverN, DOUBLE xi) const override
    {
        _RECORD(CFieldGaugeOneLoopImproved::CalculatePlaqutteEnergyAnisotropy);
        DOUBLE fPlaqTerm = CFeildG::CalculatePlaqutteEnergyAnisotropy(betaOverN, xi);

        fPlaqTerm += CFeildG::_GaugeKernel::CalculateRectangularEnergyAnisotropy(
            this->m_pDeviceData,
            this->m_byFieldId,
            m_fCr * betaOverN,
            xi);

        fPlaqTerm += CFeildG::_GaugeKernel::CalculateTwistedLoopEnergyAnisotropy(
            this->m_pDeviceData,
            this->m_byFieldId,
            m_fCt * betaOverN,
            xi);

        return fPlaqTerm;
    }

    DOUBLE CalculatePlaqutteEnergyUseCloverAnisotropy(DOUBLE betaOverN, DOUBLE xi) const override
    {
        _RECORD(CFieldGaugeOneLoopImproved::CalculatePlaqutteEnergyUseCloverAnisotropy);
        DOUBLE fPlaqTerm = CFeildG::CalculatePlaqutteEnergyUseCloverAnisotropy(betaOverN, xi);

        fPlaqTerm += CFeildG::_GaugeKernel::CalculateRectangularEnergyUseCloverAnisotropy(
            this->m_pDeviceData,
            this->m_byFieldId,
            m_fCr * betaOverN,
            xi);

        fPlaqTerm += CFeildG::_GaugeKernel::CalculateTwistedLoopEnergyAnisotropy(
            this->m_pDeviceData,
            this->m_byFieldId,
            m_fCt * betaOverN,
            xi);

        return fPlaqTerm;
    }

    CCString GetInfos(const CCString& tab) const override
    {
        CCString sRet = CFeildG::GetInfos(tab);
        if (m_bDirectSetCoefficients)
        {
            sRet = sRet + _T("Cr:") + appToString(m_fCr) + _T("\n");
            sRet = sRet + _T("Ct:") + appToString(m_fCt) + _T("\n");
        }
        else
        {
            sRet = sRet + _T("Cr:") + appToString(m_fCr) + _T("\n");
            sRet = sRet + _T("Ct:") + appToString(m_fCt) + _T("\n");
            sRet = sRet + _T("U0:") + appToString(m_fU0) + _T("\n");
            sRet = sRet + _T("Nf:") + appToString(m_iNf) + _T("\n");
        }
        
        return sRet;
    }

    DOUBLE m_fCr;
    DOUBLE m_fCt;
    DOUBLE m_fU0;
    INT m_iNf;
    UBOOL m_bDirectSetCoefficients;
};

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU3OneLoopImproved)
class CLGAPI CFieldGaugeSU3OneLoopImproved : public CFieldGaugeOneLoopImproved<CFieldGaugeSU3>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldGaugeSU3OneLoopImproved)
};

__END_NAMESPACE


#endif //#ifndef _CFIELDGAUGE_SU3_TREEIMPROVED_H_

//=============================================================================
// END OF FILE
//=============================================================================