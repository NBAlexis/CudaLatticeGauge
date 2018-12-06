//=============================================================================
// FILENAME : Random.h
// 
// DESCRIPTION:
// This is random number for parallel
//
//
// REVISION:
//  [12/6/2018 nbale]
//=============================================================================

#ifndef _RANDOM_H_
#define _RANDOM_H_

__BEGIN_NAMESPACE

class CLGAPI CRandom
{
public:
    CRandom(class CLatticeData* pOwner, UINT uiSeed) : m_pOwner(pOwner) {}
    virtual __device__ FLOAT _deviceRandomF(UINT fatIndex) = 0;
    virtual __device__ UINT _deviceRandomUI(UINT fatIndex) = 0;
protected:
    class CLatticeData* m_pOwner;
};

/**
* we cannot make it virtual inline
* so we have a random for testing usage only
* for real usage, (if that is ok), one do not need a table. The table only ensures, same seed same result.
*/
class CLGAPI CRandomSchrage
{
public:

    const FLOAT AM = (1.0 / 4294967296UL);
    
    __device__ __inline__ FLOAT _deviceRandomF(UINT fatIndex)
    {
        return AM * _deviceRandomUI(fatIndex);
    }

    __device__ __inline__ UINT _deviceRandomUI(UINT fatIndex)
    {
        m_pDeviceSeedTable[fatIndex] = (1664525UL * m_pDeviceSeedTable[fatIndex] + 1013904223UL) & 0xffffffff;
        return m_pDeviceSeedTable[fatIndex];
    }

    /**
    * run on device, parally set the table
    */
    __device__ __inline__ void _deviceAsignSeeds(UINT uiSeed, UINT uiFatIndex)
    {
        m_pDeviceSeedTable[uiFatIndex] = (1664525UL * (uiFatIndex^uiSeed) + 1013904223UL) & 0xffffffff;
    }

    friend class CLatticeData;
    friend class CDeviceLattice;

private:

    class CDeviceLattice* m_pOwner;
    CRandomSchrage(UINT uiSeed, class CDeviceLattice* pOwner);
    UINT* m_pDeviceSeedTable;
};

__END_NAMESPACE

#endif //#ifndef _RANDOM_H_

//=============================================================================
// END OF FILE
//=============================================================================
