//=============================================================================
// FILENAME : TemplateFunction.h
// 
// DESCRIPTION:
//
// REVISION:
//  [3/13/2018 nbale]
//=============================================================================
#ifndef _TEMPLATEFUNCTION_H_
#define _TEMPLATEFUNCTION_H_

__BEGIN_NAMESPACE

//-----------------------------------------------------------------------
//template function
template< class T > FORCEINLINE T appAbs( const T A )
{
    return ( A >= (T)(0) ) ? A : -A;
}
template< class T > FORCEINLINE T appSgn( const T A )
{
    return ( A > 0 ) ? 1 : ( ( A < 0 ) ? -1 : 0);
}
template< class T > FORCEINLINE T appMax( const T A, const T B )
{
    return ( A >= B ) ? A : B;
}
template< class T > FORCEINLINE T appMin( const T A, const T B )
{
    return ( A <= B ) ? A : B;
}
template< class T > FORCEINLINE T appMax3( const T A, const T B, const T C )
{
    return appMax ( appMax( A, B ), C );
}
template< class T > FORCEINLINE T appMin3( const T A, const T B, const T C )
{
    return appMin ( appMin( A, B ), C );
}
template< class T > FORCEINLINE T appSquare( const T A )
{
    return A*A;
}
template< class T > FORCEINLINE T appClamp( const T X, const T Min, const T Max )
{
    return X<Min ? Min : X<Max ? X : Max;
}
template< class T > FORCEINLINE T appAlign( const T Ptr, INT Alignment )
{
    return (T)( ( ((PTRINT)Ptr) + Alignment - 1) & ~(Alignment-1) );
}
template< class T > FORCEINLINE void appExchange( T& A, T& B )
{
    const T Temp = A;
    A = B;
    B = Temp;
}
template< class T > FORCEINLINE T appLerp( T& A, T& B, FLOAT Alpha )
{
    return (T)(A + Alpha * (B-A));
}

//-----------------------------------------------------------------------
//define
#define appExchangeB(A,B) {UBOOL T=A; A=B; B=T;}

template<typename TYPE>
FORCEINLINE void appCopyElements(TYPE* pDest, const TYPE* pSrc, INT nCount)
{
    assert(nCount>=0);
    // default is element-copy using assignment
    while (nCount--)
        *pDest++ = *pSrc++;
}

__END_NAMESPACE

#endif //#ifndef _TEMPLATEFUNCTION_H_

//=============================================================================
// END OF FILE
//=============================================================================