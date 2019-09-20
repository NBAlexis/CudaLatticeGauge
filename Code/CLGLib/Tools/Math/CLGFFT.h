//=============================================================================
// FILENAME : CLGFFT.h
// 
// DESCRIPTION:
// This is helper to calculate FFT using cufft
//
//
// REVISION:
//  [09/16/2019 nbale]
//=============================================================================

#ifndef _CLGFFT_H_
#define _CLGFFT_H_

__BEGIN_NAMESPACE

/**
 * We prefer ES_None as default, because we will have to multiply other quantities
 */
enum EFFT_Scale
{
    ES_None,
    ES_1OverNForward,
    ES_1OverNInverse,
    ES_1OverSqrtNBoth,
};

class CLGAPI CCLGFFTHelper
{
public:

    CCLGFFTHelper()
    : m_pDeviceBuffer(NULL)
    {
        
    }

    ~CCLGFFTHelper()
    {
        if (NULL != m_pDeviceBuffer)
        {
            checkCudaErrors(cudaFree(m_pDeviceBuffer));
        }
    }

    /**
     * copied is the copy of the source, will be changed
     */
    static UBOOL FFT3DWithXYZ(CLGComplex* copied, TArray<INT> dims, UBOOL bForward);
    static UBOOL FFT3DWithXYZW(CLGComplex* copied, TArray<INT> dims, UBOOL bForward);
    static UBOOL FFT4DWithXYZW(CLGComplex* copied, TArray<INT> dims, UBOOL bForward);
    static UBOOL FFT3D(CLGComplex* res, UBOOL bForward, EFFT_Scale eScale = ES_None);
    static UBOOL FFT4D(CLGComplex* res, UBOOL bForward, EFFT_Scale eScale = ES_None);

    UBOOL FFT3DSU3(deviceSU3* res, UBOOL bForward, EFFT_Scale eScale = ES_None);
    UBOOL FFT4DSU3(deviceSU3* res, UBOOL bForward, EFFT_Scale eScale = ES_None);

    /**
     * Test function
     */
    static void TestFFT();

private:

    static void GenerateTestArray(CLGComplex* hostArray, INT iSize);
    static void PrintTestArray3D(CLGComplex* hostArray);
    static void PrintTestArray4D(CLGComplex* hostArray);
    void CheckBuffer();

    CLGComplex* m_pDeviceBuffer;
};

__END_NAMESPACE

#endif //#ifndef _CLGFFT_H_

//=============================================================================
// END OF FILE
//=============================================================================
