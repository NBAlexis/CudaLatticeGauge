//=============================================================================
// FILENAME : CLGMD5.h
// 
// DESCRIPTION:
//  try not to depend on other libs
//
// REVISION:
//  [1/21/2021 nbale]
//=============================================================================
#pragma once

#ifndef _CLGMD5_H_
#define _CLGMD5_H_

__BEGIN_NAMESPACE

enum { md5_S11 = 7 };
enum { md5_S12 = 12 };
enum { md5_S13 = 17 };
enum { md5_S14 = 22 };

enum { md5_S21 = 5 };
enum { md5_S22 = 9 };
enum { md5_S23 = 14 };
enum { md5_S24 = 20 };

enum { md5_S31 = 4 };
enum { md5_S32 = 11 };
enum { md5_S33 = 16 };
enum { md5_S34 = 23 };

enum { md5_S41 = 6 };
enum { md5_S42 = 10 };
enum { md5_S43 = 15 };
enum { md5_S44 = 21 };

#define MD5_F(xx, yy, zz) (((xx) & (yy)) | ((~xx) & (zz)))
#define MD5_G(xx, yy, zz) (((xx) & (zz)) | ((yy) & (~zz)))
#define MD5_H(xx, yy, zz) ((xx) ^ (yy) ^ (zz))
#define MD5_I(xx, yy, zz) ((yy) ^ ((xx) | (~zz)))

#define ROTATE_LEFT(xx, nn) (((xx) << (nn)) | ((xx) >> (32 - (nn))))

#define MD5_FF(aa, bb, cc, dd, xx, ss, ac) {\
    (aa) += MD5_F ((bb), (cc), (dd)) + (xx) + (QWORD)(ac ); \
    (aa) = ROTATE_LEFT ((aa), (ss)); \
    (aa) += (bb);  }

#define MD5_GG(aa, bb, cc, dd, xx, ss, ac) { \
    (aa) += MD5_G ((bb), (cc), (dd)) + (xx) + (QWORD)(ac); \
    (aa) = ROTATE_LEFT ((aa), (ss)); \
    (aa) += (bb); }

#define MD5_HH(aa, bb, cc, dd, xx, ss, ac) { \
    (aa) += MD5_H ((bb), (cc), (dd)) + (xx) + (QWORD)(ac); \
    (aa) = ROTATE_LEFT ((aa), (ss)); \
    (aa) += (bb); }

#define MD5_II(aa, bb, cc, dd, xx, ss, ac) { \
    (aa) += MD5_I ((bb), (cc), (dd)) + (xx) + (QWORD)(ac); \
    (aa) = ROTATE_LEFT ((aa), (ss)); \
    (aa) += (bb); }



/**
*
*/
static QWORD* FoldDataMD5_OLD(UINT& iBlockCount, const BYTE* pInData, UINT iDataCount)
{
    const UINT iLength = (UINT)((INT)(iDataCount >> 6) + 1) * 16; //blocks is parted into every 64 chars, and so 16 DWORDs as a block
    QWORD* retVal = (QWORD*)malloc(sizeof(QWORD) * iLength);
    assert(retVal);
    memset(retVal, 0, sizeof(QWORD) * iLength);
    iBlockCount = iLength;

    for (UINT i = 0; i < iDataCount; ++i)
    {
        //every 4 char combine one 32-bit DWORD, and count of relVal is 16s
        const UINT iMove = i & 0x00000003;
        retVal[i >> 2] |= (QWORD)pInData[i] << (iMove * 8);
    }
    //AT LAST I REALIZED THAT THIS IS WHAT memcpy JUST DO! - Aleixs

    UINT iFinalCount = (((iDataCount + 8) >> 6) << 4) + 15;
    iFinalCount = iFinalCount < iLength ? iLength - 1 : iFinalCount;
    const UINT iNewLength = ((iFinalCount >> 4) + 1) * 16;
    QWORD* pNewsValue = (QWORD*)malloc(sizeof(QWORD) * iNewLength);
    memset(pNewsValue, 0, sizeof(QWORD) * iNewLength);
    memcpy(pNewsValue, retVal, sizeof(QWORD) * iLength);
    pNewsValue[iDataCount >> 2] |= 128 << (8 * (iDataCount & 0x00000003)); //Add length of string to end of array
    pNewsValue[(((iDataCount + 8) >> 6) << 4) + 14] = iDataCount * 8; //Add length of string to end of array why?

    free(retVal);

    iBlockCount = iNewLength;
    return pNewsValue;
}

static QWORD* MD5HashFoldeded_OLD(QWORD* pFoldedData, UINT iFoldCount)
{
    /*
    appAssert(1732584193 == 0x67452301);
    appAssert(-271733879 == 0xefcdab89);
    appAssert(-1732584194 == 0x98badcfe);
    appAssert(271733878 == 0x10325476);

    appAssert(-680876936 == 0xd76aa478);
    appAssert(-389564586 == 0xe8c7b756);
    appAssert(606105819 == 0x242070db);
    appAssert(-1044525330 == 0xc1bdceee);
    appAssert(-176418897 == 0xf57c0faf);
    appAssert(1200080426 == 0x4787c62a);
    appAssert(-1473231341 == 0xa8304613);
    appAssert(-45705983 == 0xfd469501);
    appAssert(1770035416 == 0x698098d8);
    appAssert(-1958414417 == 0x8b44f7af);
    appAssert(-42063 == 0xffff5bb1);
    appAssert(-1990404162 == 0x895cd7be);
    appAssert(1804603682 == 0x6b901122);
    appAssert(-40341101 == 0xfd987193);
    appAssert(-1502002290 == 0xa679438e);
    appAssert(1236535329 == 0x49b40821);

    appAssert(-165796510 == 0xf61e2562);
    appAssert(-1069501632 == 0xc040b340);
    appAssert(643717713 == 0x265e5a51);
    appAssert(-373897302 == 0xe9b6c7aa);
    appAssert(-701558691 == 0xd62f105d);
    appAssert(38016083 == 0x2441453);
    appAssert(-660478335 == 0xd8a1e681);
    appAssert(-405537848 == 0xe7d3fbc8);
    appAssert(568446438 == 0x21e1cde6);
    appAssert(-1019803690 == 0xc33707d6);
    appAssert(-187363961 == 0xf4d50d87);
    appAssert(1163531501 == 0x455a14ed);
    appAssert(-1444681467 == 0xa9e3e905);
    appAssert(-51403784 == 0xfcefa3f8);
    appAssert(1735328473 == 0x676f02d9);
    appAssert(-1926607734 == 0x8d2a4c8a);

    appAssert(-378558 == 0xfffa3942);
    appAssert(-2022574463 == 0x8771f681);
    appAssert(1839030562 == 0x6d9d6122);
    appAssert(-35309556 == 0xfde5380c);
    appAssert(-1530992060 == 0xa4beea44);
    appAssert(1272893353 == 0x4bdecfa9);
    appAssert(-155497632 == 0xf6bb4b60);
    appAssert(-1094730640 == 0xbebfbc70);
    appAssert(681279174 == 0x289b7ec6);
    appAssert(-358537222 == 0xeaa127fa);
    appAssert(-722521979 == 0xd4ef3085);
    appAssert(76029189 == 0x4881d05);
    appAssert(-640364487 == 0xd9d4d039);
    appAssert(-421815835 == 0xe6db99e5);
    appAssert(530742520 == 0x1fa27cf8);
    appAssert(-995338651 == 0xc4ac5665);

    appAssert(-198630844 == 0xf4292244);
    appAssert(1126891415 == 0x432aff97);
    appAssert(-1416354905 == 0xab9423a7);
    appAssert(-57434055 == 0xfc93a039);
    appAssert(1700485571 == 0x655b59c3);
    appAssert(-1894986606 == 0x8f0ccc92);
    appAssert(-1051523 == 0xffeff47d);
    appAssert(-2054922799 == 0x85845dd1);
    appAssert(1873313359 == 0x6fa87e4f);
    appAssert(-30611744 == 0xfe2ce6e0);
    appAssert(-1560198380 == 0xa3014314);
    appAssert(1309151649 == 0x4e0811a1);
    appAssert(-145523070 == 0xf7537e82);
    appAssert(-1120210379 == 0xbd3af235);
    appAssert(718787259 == 0x2ad7d2bb);
    appAssert(-343485551 == 0xeb86d391);
    */

    QWORD* pOutValue = (QWORD*)malloc(sizeof(QWORD) * 4);
    assert(pOutValue);

    pOutValue[0] = 0x67452301;
    pOutValue[1] = 0xefcdab89;
    pOutValue[2] = 0x98badcfe;
    pOutValue[3] = 0x10325476;

    assert((iFoldCount & 0x0000000f) == 0); //iLength must be 16s

    for (UINT i = 0; i < (iFoldCount >> 4); ++i)
    {
        QWORD a = pOutValue[0];
        QWORD b = pOutValue[1];
        QWORD c = pOutValue[2];
        QWORD d = pOutValue[3];
        QWORD x[16];

        memcpy(x, pFoldedData + 16 * i, sizeof(QWORD) * 16);

        /* Round 1 */
        MD5_FF(a, b, c, d, x[0], md5_S11, 0xd76aa478); /* 1 */
        MD5_FF(d, a, b, c, x[1], md5_S12, 0xe8c7b756); /* 2 */
        MD5_FF(c, d, a, b, x[2], md5_S13, 0x242070db); /* 3 */
        MD5_FF(b, c, d, a, x[3], md5_S14, 0xc1bdceee); /* 4 */
        MD5_FF(a, b, c, d, x[4], md5_S11, 0xf57c0faf); /* 5 */
        MD5_FF(d, a, b, c, x[5], md5_S12, 0x4787c62a); /* 6 */
        MD5_FF(c, d, a, b, x[6], md5_S13, 0xa8304613); /* 7 */
        MD5_FF(b, c, d, a, x[7], md5_S14, 0xfd469501); /* 8 */
        MD5_FF(a, b, c, d, x[8], md5_S11, 0x698098d8); /* 9 */
        MD5_FF(d, a, b, c, x[9], md5_S12, 0x8b44f7af); /* 10 */
        MD5_FF(c, d, a, b, x[10], md5_S13, 0xffff5bb1); /* 11 */
        MD5_FF(b, c, d, a, x[11], md5_S14, 0x895cd7be); /* 12 */
        MD5_FF(a, b, c, d, x[12], md5_S11, 0x6b901122); /* 13 */
        MD5_FF(d, a, b, c, x[13], md5_S12, 0xfd987193); /* 14 */
        MD5_FF(c, d, a, b, x[14], md5_S13, 0xa679438e); /* 15 */
        MD5_FF(b, c, d, a, x[15], md5_S14, 0x49b40821); /* 16 */

        /* Round 2 */
        MD5_GG(a, b, c, d, x[1], md5_S21, 0xf61e2562); /* 17 */
        MD5_GG(d, a, b, c, x[6], md5_S22, 0xc040b340); /* 18 */
        MD5_GG(c, d, a, b, x[11], md5_S23, 0x265e5a51); /* 19 */
        MD5_GG(b, c, d, a, x[0], md5_S24, 0xe9b6c7aa); /* 20 */
        MD5_GG(a, b, c, d, x[5], md5_S21, 0xd62f105d); /* 21 */
        MD5_GG(d, a, b, c, x[10], md5_S22, 0x2441453); /* 22 */
        MD5_GG(c, d, a, b, x[15], md5_S23, 0xd8a1e681); /* 23 */
        MD5_GG(b, c, d, a, x[4], md5_S24, 0xe7d3fbc8); /* 24 */
        MD5_GG(a, b, c, d, x[9], md5_S21, 0x21e1cde6); /* 25 */
        MD5_GG(d, a, b, c, x[14], md5_S22, 0xc33707d6); /* 26 */
        MD5_GG(c, d, a, b, x[3], md5_S23, 0xf4d50d87); /* 27 */
        MD5_GG(b, c, d, a, x[8], md5_S24, 0x455a14ed); /* 28 */
        MD5_GG(a, b, c, d, x[13], md5_S21, 0xa9e3e905); /* 29 */
        MD5_GG(d, a, b, c, x[2], md5_S22, 0xfcefa3f8); /* 30 */
        MD5_GG(c, d, a, b, x[7], md5_S23, 0x676f02d9); /* 31 */
        MD5_GG(b, c, d, a, x[12], md5_S24, 0x8d2a4c8a); /* 32 */

        /* Round 3 */
        MD5_HH(a, b, c, d, x[5], md5_S31, 0xfffa3942); /* 33 */
        MD5_HH(d, a, b, c, x[8], md5_S32, 0x8771f681); /* 34 */
        MD5_HH(c, d, a, b, x[11], md5_S33, 0x6d9d6122); /* 35 */
        MD5_HH(b, c, d, a, x[14], md5_S34, 0xfde5380c); /* 36 */
        MD5_HH(a, b, c, d, x[1], md5_S31, 0xa4beea44); /* 37 */
        MD5_HH(d, a, b, c, x[4], md5_S32, 0x4bdecfa9); /* 38 */
        MD5_HH(c, d, a, b, x[7], md5_S33, 0xf6bb4b60); /* 39 */
        MD5_HH(b, c, d, a, x[10], md5_S34, 0xbebfbc70); /* 40 */
        MD5_HH(a, b, c, d, x[13], md5_S31, 0x289b7ec6); /* 41 */
        MD5_HH(d, a, b, c, x[0], md5_S32, 0xeaa127fa); /* 42 */
        MD5_HH(c, d, a, b, x[3], md5_S33, 0xd4ef3085); /* 43 */
        MD5_HH(b, c, d, a, x[6], md5_S34, 0x4881d05); /* 44 */
        MD5_HH(a, b, c, d, x[9], md5_S31, 0xd9d4d039); /* 45 */
        MD5_HH(d, a, b, c, x[12], md5_S32, 0xe6db99e5); /* 46 */
        MD5_HH(c, d, a, b, x[15], md5_S33, 0x1fa27cf8); /* 47 */
        MD5_HH(b, c, d, a, x[2], md5_S34, 0xc4ac5665); /* 48 */

        /* Round 4 */
        MD5_II(a, b, c, d, x[0], md5_S41, 0xf4292244); /* 49 */
        MD5_II(d, a, b, c, x[7], md5_S42, 0x432aff97); /* 50 */
        MD5_II(c, d, a, b, x[14], md5_S43, 0xab9423a7); /* 51 */
        MD5_II(b, c, d, a, x[5], md5_S44, 0xfc93a039); /* 52 */
        MD5_II(a, b, c, d, x[12], md5_S41, 0x655b59c3); /* 53 */
        MD5_II(d, a, b, c, x[3], md5_S42, 0x8f0ccc92); /* 54 */
        MD5_II(c, d, a, b, x[10], md5_S43, 0xffeff47d); /* 55 */
        MD5_II(b, c, d, a, x[1], md5_S44, 0x85845dd1); /* 56 */
        MD5_II(a, b, c, d, x[8], md5_S41, 0x6fa87e4f); /* 57 */
        MD5_II(d, a, b, c, x[15], md5_S42, 0xfe2ce6e0); /* 58 */
        MD5_II(c, d, a, b, x[6], md5_S43, 0xa3014314); /* 59 */
        MD5_II(b, c, d, a, x[13], md5_S44, 0x4e0811a1); /* 60 */
        MD5_II(a, b, c, d, x[4], md5_S41, 0xf7537e82); /* 61 */
        MD5_II(d, a, b, c, x[11], md5_S42, 0xbd3af235); /* 62 */
        MD5_II(c, d, a, b, x[2], md5_S43, 0x2ad7d2bb); /* 63 */
        MD5_II(b, c, d, a, x[9], md5_S44, 0xeb86d391); /* 64 */

        pOutValue[0] += a;
        pOutValue[1] += b;
        pOutValue[2] += c;
        pOutValue[3] += d;

    }

    return pOutValue;
}

/**
 * Assuming the length of pData is 4 x bytes
 */
static CCString CLGMD5Hash_OLD(const BYTE * pData, UINT uiDataCount)
{
    UINT iBlockCount = 0;
    QWORD* folded = FoldDataMD5_OLD(iBlockCount, pData, uiDataCount);
    QWORD* pResoult = MD5HashFoldeded_OLD(folded, iBlockCount);

    const UINT iLength = 4;
    CCString sMid;
    for (UINT i = 0; i < (iLength << 3); ++i)
    {
        const BYTE byCurrent = (pResoult[i >> 3] >> ((i & 0x00000007) * 4)) & 0x0000000f;
        const TCHAR sThisChar = byCurrent > 9 ? byCurrent - 10 + _T('A') : byCurrent + _T('0');
        sMid += (TCHAR)sThisChar;
    }

    CCString sRet;
    for (UINT i = 0; i < (iLength << 2); ++i)
    {
        const TCHAR char1 = sMid.GetAt(2 * i);
        const TCHAR char2 = sMid.GetAt(2 * i + 1);
        sRet += char2;
        sRet += char1;
    }

    free(pResoult);
    free(folded);

    return sRet;
}




//=========================================================================
// The above is wrong on Ubuntu, preserved to check the old files
//=========================================================================
#undef MD5_FF
#undef MD5_GG
#undef MD5_HH
#undef MD5_II

#define MD5_FF(aa, bb, cc, dd, xx, ss, ac) {\
    (aa) += MD5_F ((bb), (cc), (dd)) + (xx) + (UINT)(ac ); \
    (aa) = ROTATE_LEFT ((aa), (ss)); \
    (aa) += (bb);  }

#define MD5_GG(aa, bb, cc, dd, xx, ss, ac) { \
    (aa) += MD5_G ((bb), (cc), (dd)) + (xx) + (UINT)(ac); \
    (aa) = ROTATE_LEFT ((aa), (ss)); \
    (aa) += (bb); }

#define MD5_HH(aa, bb, cc, dd, xx, ss, ac) { \
    (aa) += MD5_H ((bb), (cc), (dd)) + (xx) + (UINT)(ac); \
    (aa) = ROTATE_LEFT ((aa), (ss)); \
    (aa) += (bb); }

#define MD5_II(aa, bb, cc, dd, xx, ss, ac) { \
    (aa) += MD5_I ((bb), (cc), (dd)) + (xx) + (UINT)(ac); \
    (aa) = ROTATE_LEFT ((aa), (ss)); \
    (aa) += (bb); }




/**
*
*/
static UINT* FoldDataMD5(UINT& iBlockCount, const BYTE* pInData, UINT iDataCount)
{
    const UINT iLength = (UINT)((INT)(iDataCount >> 6) + 1) * 16; //blocks is parted into every 64 chars, and so 16 DWORDs as a block
    UINT* retVal = (UINT*)malloc(sizeof(UINT) * iLength);
    assert(retVal);
    memset(retVal, 0, sizeof(UINT) * iLength);
    iBlockCount = iLength;

    for (UINT i = 0; i < iDataCount; ++i)
    {
        //every 4 char combine one 32-bit DWORD, and count of relVal is 16s
        const UINT iMove = i & 0x00000003;
        retVal[i >> 2] |= (UINT)pInData[i] << (iMove * 8);
    }
    //AT LAST I REALIZED THAT THIS IS WHAT memcpy JUST DO! - Aleixs

    UINT iFinalCount = (((iDataCount + 8) >> 6) << 4) + 15;
    iFinalCount = iFinalCount < iLength ? iLength - 1 : iFinalCount;
    const UINT iNewLength = ((iFinalCount >> 4) + 1) * 16;
    UINT* pNewsValue = (UINT*)malloc(sizeof(UINT) * iNewLength);
    memset(pNewsValue, 0, sizeof(UINT) * iNewLength);
    memcpy(pNewsValue, retVal, sizeof(UINT) * iLength);
    pNewsValue[iDataCount >> 2] |= 128 << (8 * (iDataCount & 0x00000003)); //Add length of string to end of array
    pNewsValue[(((iDataCount + 8) >> 6) << 4) + 14] = iDataCount * 8; //Add length of string to end of array why?

    free(retVal);

    iBlockCount = iNewLength;
    return pNewsValue;
}

static UINT* MD5HashFoldeded(UINT* pFoldedData, UINT iFoldCount)
{
    /*
    appAssert(1732584193 == 0x67452301);
    appAssert(-271733879 == 0xefcdab89);
    appAssert(-1732584194 == 0x98badcfe);
    appAssert(271733878 == 0x10325476);

    appAssert(-680876936 == 0xd76aa478);
    appAssert(-389564586 == 0xe8c7b756);
    appAssert(606105819 == 0x242070db);
    appAssert(-1044525330 == 0xc1bdceee);
    appAssert(-176418897 == 0xf57c0faf);
    appAssert(1200080426 == 0x4787c62a);
    appAssert(-1473231341 == 0xa8304613);
    appAssert(-45705983 == 0xfd469501);
    appAssert(1770035416 == 0x698098d8);
    appAssert(-1958414417 == 0x8b44f7af);
    appAssert(-42063 == 0xffff5bb1);
    appAssert(-1990404162 == 0x895cd7be);
    appAssert(1804603682 == 0x6b901122);
    appAssert(-40341101 == 0xfd987193);
    appAssert(-1502002290 == 0xa679438e);
    appAssert(1236535329 == 0x49b40821);

    appAssert(-165796510 == 0xf61e2562);
    appAssert(-1069501632 == 0xc040b340);
    appAssert(643717713 == 0x265e5a51);
    appAssert(-373897302 == 0xe9b6c7aa);
    appAssert(-701558691 == 0xd62f105d);
    appAssert(38016083 == 0x2441453);
    appAssert(-660478335 == 0xd8a1e681);
    appAssert(-405537848 == 0xe7d3fbc8);
    appAssert(568446438 == 0x21e1cde6);
    appAssert(-1019803690 == 0xc33707d6);
    appAssert(-187363961 == 0xf4d50d87);
    appAssert(1163531501 == 0x455a14ed);
    appAssert(-1444681467 == 0xa9e3e905);
    appAssert(-51403784 == 0xfcefa3f8);
    appAssert(1735328473 == 0x676f02d9);
    appAssert(-1926607734 == 0x8d2a4c8a);

    appAssert(-378558 == 0xfffa3942);
    appAssert(-2022574463 == 0x8771f681);
    appAssert(1839030562 == 0x6d9d6122);
    appAssert(-35309556 == 0xfde5380c);
    appAssert(-1530992060 == 0xa4beea44);
    appAssert(1272893353 == 0x4bdecfa9);
    appAssert(-155497632 == 0xf6bb4b60);
    appAssert(-1094730640 == 0xbebfbc70);
    appAssert(681279174 == 0x289b7ec6);
    appAssert(-358537222 == 0xeaa127fa);
    appAssert(-722521979 == 0xd4ef3085);
    appAssert(76029189 == 0x4881d05);
    appAssert(-640364487 == 0xd9d4d039);
    appAssert(-421815835 == 0xe6db99e5);
    appAssert(530742520 == 0x1fa27cf8);
    appAssert(-995338651 == 0xc4ac5665);

    appAssert(-198630844 == 0xf4292244);
    appAssert(1126891415 == 0x432aff97);
    appAssert(-1416354905 == 0xab9423a7);
    appAssert(-57434055 == 0xfc93a039);
    appAssert(1700485571 == 0x655b59c3);
    appAssert(-1894986606 == 0x8f0ccc92);
    appAssert(-1051523 == 0xffeff47d);
    appAssert(-2054922799 == 0x85845dd1);
    appAssert(1873313359 == 0x6fa87e4f);
    appAssert(-30611744 == 0xfe2ce6e0);
    appAssert(-1560198380 == 0xa3014314);
    appAssert(1309151649 == 0x4e0811a1);
    appAssert(-145523070 == 0xf7537e82);
    appAssert(-1120210379 == 0xbd3af235);
    appAssert(718787259 == 0x2ad7d2bb);
    appAssert(-343485551 == 0xeb86d391);
    */

    UINT* pOutValue = (UINT*)malloc(sizeof(UINT) * 4);
    assert(pOutValue);

    pOutValue[0] = 0x67452301;
    pOutValue[1] = 0xefcdab89;
    pOutValue[2] = 0x98badcfe;
    pOutValue[3] = 0x10325476;

    assert((iFoldCount & 0x0000000f) == 0); //iLength must be 16s

    for (UINT i = 0; i < (iFoldCount >> 4); ++i)
    {
        UINT a = pOutValue[0];
        UINT b = pOutValue[1];
        UINT c = pOutValue[2];
        UINT d = pOutValue[3];
        UINT x[16];

        memcpy(x, pFoldedData + 16 * i, sizeof(UINT) * 16);

        /* Round 1 */
        MD5_FF(a, b, c, d, x[0], md5_S11, 0xd76aa478); /* 1 */
        MD5_FF(d, a, b, c, x[1], md5_S12, 0xe8c7b756); /* 2 */
        MD5_FF(c, d, a, b, x[2], md5_S13, 0x242070db); /* 3 */
        MD5_FF(b, c, d, a, x[3], md5_S14, 0xc1bdceee); /* 4 */
        MD5_FF(a, b, c, d, x[4], md5_S11, 0xf57c0faf); /* 5 */
        MD5_FF(d, a, b, c, x[5], md5_S12, 0x4787c62a); /* 6 */
        MD5_FF(c, d, a, b, x[6], md5_S13, 0xa8304613); /* 7 */
        MD5_FF(b, c, d, a, x[7], md5_S14, 0xfd469501); /* 8 */
        MD5_FF(a, b, c, d, x[8], md5_S11, 0x698098d8); /* 9 */
        MD5_FF(d, a, b, c, x[9], md5_S12, 0x8b44f7af); /* 10 */
        MD5_FF(c, d, a, b, x[10], md5_S13, 0xffff5bb1); /* 11 */
        MD5_FF(b, c, d, a, x[11], md5_S14, 0x895cd7be); /* 12 */
        MD5_FF(a, b, c, d, x[12], md5_S11, 0x6b901122); /* 13 */
        MD5_FF(d, a, b, c, x[13], md5_S12, 0xfd987193); /* 14 */
        MD5_FF(c, d, a, b, x[14], md5_S13, 0xa679438e); /* 15 */
        MD5_FF(b, c, d, a, x[15], md5_S14, 0x49b40821); /* 16 */

        /* Round 2 */
        MD5_GG(a, b, c, d, x[1], md5_S21, 0xf61e2562); /* 17 */
        MD5_GG(d, a, b, c, x[6], md5_S22, 0xc040b340); /* 18 */
        MD5_GG(c, d, a, b, x[11], md5_S23, 0x265e5a51); /* 19 */
        MD5_GG(b, c, d, a, x[0], md5_S24, 0xe9b6c7aa); /* 20 */
        MD5_GG(a, b, c, d, x[5], md5_S21, 0xd62f105d); /* 21 */
        MD5_GG(d, a, b, c, x[10], md5_S22, 0x2441453); /* 22 */
        MD5_GG(c, d, a, b, x[15], md5_S23, 0xd8a1e681); /* 23 */
        MD5_GG(b, c, d, a, x[4], md5_S24, 0xe7d3fbc8); /* 24 */
        MD5_GG(a, b, c, d, x[9], md5_S21, 0x21e1cde6); /* 25 */
        MD5_GG(d, a, b, c, x[14], md5_S22, 0xc33707d6); /* 26 */
        MD5_GG(c, d, a, b, x[3], md5_S23, 0xf4d50d87); /* 27 */
        MD5_GG(b, c, d, a, x[8], md5_S24, 0x455a14ed); /* 28 */
        MD5_GG(a, b, c, d, x[13], md5_S21, 0xa9e3e905); /* 29 */
        MD5_GG(d, a, b, c, x[2], md5_S22, 0xfcefa3f8); /* 30 */
        MD5_GG(c, d, a, b, x[7], md5_S23, 0x676f02d9); /* 31 */
        MD5_GG(b, c, d, a, x[12], md5_S24, 0x8d2a4c8a); /* 32 */

        /* Round 3 */
        MD5_HH(a, b, c, d, x[5], md5_S31, 0xfffa3942); /* 33 */
        MD5_HH(d, a, b, c, x[8], md5_S32, 0x8771f681); /* 34 */
        MD5_HH(c, d, a, b, x[11], md5_S33, 0x6d9d6122); /* 35 */
        MD5_HH(b, c, d, a, x[14], md5_S34, 0xfde5380c); /* 36 */
        MD5_HH(a, b, c, d, x[1], md5_S31, 0xa4beea44); /* 37 */
        MD5_HH(d, a, b, c, x[4], md5_S32, 0x4bdecfa9); /* 38 */
        MD5_HH(c, d, a, b, x[7], md5_S33, 0xf6bb4b60); /* 39 */
        MD5_HH(b, c, d, a, x[10], md5_S34, 0xbebfbc70); /* 40 */
        MD5_HH(a, b, c, d, x[13], md5_S31, 0x289b7ec6); /* 41 */
        MD5_HH(d, a, b, c, x[0], md5_S32, 0xeaa127fa); /* 42 */
        MD5_HH(c, d, a, b, x[3], md5_S33, 0xd4ef3085); /* 43 */
        MD5_HH(b, c, d, a, x[6], md5_S34, 0x4881d05); /* 44 */
        MD5_HH(a, b, c, d, x[9], md5_S31, 0xd9d4d039); /* 45 */
        MD5_HH(d, a, b, c, x[12], md5_S32, 0xe6db99e5); /* 46 */
        MD5_HH(c, d, a, b, x[15], md5_S33, 0x1fa27cf8); /* 47 */
        MD5_HH(b, c, d, a, x[2], md5_S34, 0xc4ac5665); /* 48 */

        /* Round 4 */
        MD5_II(a, b, c, d, x[0], md5_S41, 0xf4292244); /* 49 */
        MD5_II(d, a, b, c, x[7], md5_S42, 0x432aff97); /* 50 */
        MD5_II(c, d, a, b, x[14], md5_S43, 0xab9423a7); /* 51 */
        MD5_II(b, c, d, a, x[5], md5_S44, 0xfc93a039); /* 52 */
        MD5_II(a, b, c, d, x[12], md5_S41, 0x655b59c3); /* 53 */
        MD5_II(d, a, b, c, x[3], md5_S42, 0x8f0ccc92); /* 54 */
        MD5_II(c, d, a, b, x[10], md5_S43, 0xffeff47d); /* 55 */
        MD5_II(b, c, d, a, x[1], md5_S44, 0x85845dd1); /* 56 */
        MD5_II(a, b, c, d, x[8], md5_S41, 0x6fa87e4f); /* 57 */
        MD5_II(d, a, b, c, x[15], md5_S42, 0xfe2ce6e0); /* 58 */
        MD5_II(c, d, a, b, x[6], md5_S43, 0xa3014314); /* 59 */
        MD5_II(b, c, d, a, x[13], md5_S44, 0x4e0811a1); /* 60 */
        MD5_II(a, b, c, d, x[4], md5_S41, 0xf7537e82); /* 61 */
        MD5_II(d, a, b, c, x[11], md5_S42, 0xbd3af235); /* 62 */
        MD5_II(c, d, a, b, x[2], md5_S43, 0x2ad7d2bb); /* 63 */
        MD5_II(b, c, d, a, x[9], md5_S44, 0xeb86d391); /* 64 */

        pOutValue[0] += a;
        pOutValue[1] += b;
        pOutValue[2] += c;
        pOutValue[3] += d;

    }

    return pOutValue;
}

/**
 * Assuming the length of pData is 4 x bytes
 */
static CCString CLGMD5Hash(const BYTE* pData, UINT uiDataCount)
{
    UINT iBlockCount = 0;
    UINT* folded = FoldDataMD5(iBlockCount, pData, uiDataCount);
    UINT* pResoult = MD5HashFoldeded(folded, iBlockCount);

    const UINT iLength = 4;
    CCString sMid;
    for (UINT i = 0; i < (iLength << 3); ++i)
    {
        const BYTE byCurrent = (pResoult[i >> 3] >> ((i & 0x00000007) * 4)) & 0x0000000f;
        const TCHAR sThisChar = byCurrent > 9 ? byCurrent - 10 + _T('A') : byCurrent + _T('0');
        sMid += (TCHAR)sThisChar;
    }

    CCString sRet;
    for (UINT i = 0; i < (iLength << 2); ++i)
    {
        const TCHAR char1 = sMid.GetAt(2 * i);
        const TCHAR char2 = sMid.GetAt(2 * i + 1);
        sRet += char2;
        sRet += char1;
    }

    free(pResoult);
    free(folded);

    return sRet;
}


__END_NAMESPACE

#endif//#ifndef _CLGMD5_H_

//=============================================================================
// END OF FILE
//=============================================================================