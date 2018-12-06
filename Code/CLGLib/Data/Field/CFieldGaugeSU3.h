//=============================================================================
// FILENAME : CFieldGaugeSU3.h
// 
// DESCRIPTION:
// This is the class for the gauge fields
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CFIELDGAUGE_SU3_H_
#define _CFIELDGAUGE_SU3_H_

__BEGIN_NAMESPACE

class CLGAPI CFieldGaugeSU3 : public CFieldGauge
{
public:

    CFieldGaugeSU3(CLatticeData* pLattice, UBOOL bHot = FALSE);
    ~CFieldGaugeSU3();

    virtual EFieldType GetFieldType() const { return EFT_GaugeSU3; }
    virtual void CalculateStaple(void);
    //virtual void axpy(FLOAT a, const CField *x) { axpy(make_cuComplex(a, 0.0f), x); }
    //virtual void axpy(const cuComplex& a, const CField *x);

protected:

    /**
    * SU3(x=(x,y,z,t))_{n=a*3+b}=
    * m_pData[( 
        (x*m_uiDim[1]*m_uiDim[2]*m_uiDim[3] + y*m_uiDim[2]*m_uiDim[3] + z*m_uiDim[3] + t)
            * m_uiDir + dir) * 9 + n]
    */
    deviceSU3* m_pDeviceData;
    deviceSU3* m_pDeviceStaple;

    /**
    * For each dir, there exist a plaquette start from that link of that site
    * To calculate that, one need a action to index the plaquettes
    * Also, one need a boundary condition mapping
    */
    FLOAT* m_pDevicePlaquetteEnergy;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_SU3_H_

//=============================================================================
// END OF FILE
//=============================================================================