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

class CLGAPI CCLGFFTHelper
{
public:

    /**
     * copied is the copy of the source, will be changed
     */
    static UBOOL FFT3DWithXYZW(CLGComplex* copied, TArray<INT> dims, UBOOL bForward);
    static UBOOL FFT4DWithXYZW(CLGComplex* copied, TArray<INT> dims, UBOOL bForward);
    static UBOOL FFT3D(CLGComplex* res, UBOOL bForward);
    static UBOOL FFT4D(CLGComplex* res, UBOOL bForward);

    /**
     * Test function
     */
    static void TestFFT();

private:

    static void GenerateTestArray(CLGComplex* hostArray, INT iSize);
    static void PrintTestArray4D(CLGComplex* hostArray);
};

__END_NAMESPACE

#endif //#ifndef _CLGFFT_H_

//=============================================================================
// END OF FILE
//=============================================================================
