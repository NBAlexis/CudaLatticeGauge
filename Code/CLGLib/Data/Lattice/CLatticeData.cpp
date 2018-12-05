//=============================================================================
// FILENAME : CLatticeData.cpp
// 
// DESCRIPTION:
// This is the class for the lattce data
// NOTE:: We only have 4D case, 3D = 1xLxLxL, and 2D= 1x1xLxL
// REVISION:
//  [12/4/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

/**
* m_uiLatticeDecompose[0,1,2] is the blocks
* m_uiLatticeDecompose[3,4,5] is the threads in blocks
*/
CLatticeData::CLatticeData()
    : m_pIndex(NULL)
{
    m_uiDim = CCommonData::m_uiDim;
    m_uiDir = CCommonData::m_uiDir;
    memcpy(m_uiLatticeLength, CCommonData::m_uiLatticeLength, sizeof(UINT) * CCommonData::kMaxDim);

    m_uiLatticeMultipy[0] = m_uiLatticeLength[1] * m_uiLatticeLength[2] * m_uiLatticeLength[3];
    m_uiLatticeMultipy[1] = m_uiLatticeLength[2] * m_uiLatticeLength[3];
    m_uiLatticeMultipy[2] = m_uiLatticeLength[3];
    m_uiLatticeMultipy[3] = 1;
    m_uiLatticeMultipy[4] = 1;

    if (m_uiLatticeLength[0] * m_uiLatticeLength[1] * m_uiLatticeLength[2] <= CCommonData::m_uiMaxThread)
    {
        m_uiLatticeDecompose[0] = 1;
        m_uiLatticeDecompose[1] = 1;
        m_uiLatticeDecompose[2] = 1;

        m_uiLatticeDecompose[3] = m_uiLatticeLength[0];
        m_uiLatticeDecompose[4] = m_uiLatticeLength[1];
        m_uiLatticeDecompose[5] = m_uiLatticeLength[2];
    }
    else if (m_uiLatticeLength[1] * m_uiLatticeLength[2] <= CCommonData::m_uiMaxThread)
    {
        m_uiLatticeDecompose[0] = m_uiLatticeLength[0];
        m_uiLatticeDecompose[1] = 1;
        m_uiLatticeDecompose[2] = 1;

        m_uiLatticeDecompose[3] = 1;
        m_uiLatticeDecompose[4] = m_uiLatticeLength[1];
        m_uiLatticeDecompose[5] = m_uiLatticeLength[2];
    }
    else if (m_uiLatticeLength[2] <= CCommonData::m_uiMaxThread)
    {
        m_uiLatticeDecompose[0] = m_uiLatticeLength[0];
        m_uiLatticeDecompose[1] = m_uiLatticeLength[1];
        m_uiLatticeDecompose[2] = 1;

        m_uiLatticeDecompose[3] = 1;
        m_uiLatticeDecompose[4] = 1;
        m_uiLatticeDecompose[5] = m_uiLatticeLength[2];
    }
    else
    {
        appCrucial("CLatticeData:: Fail to divide the blocks!");
        m_uiLatticeDecompose[0] = m_uiLatticeLength[0];
        m_uiLatticeDecompose[1] = m_uiLatticeLength[1];
        m_uiLatticeDecompose[2] = m_uiLatticeLength[2];

        m_uiLatticeDecompose[3] = 1;
        m_uiLatticeDecompose[4] = 1;
        m_uiLatticeDecompose[5] = 1;
    }
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================

