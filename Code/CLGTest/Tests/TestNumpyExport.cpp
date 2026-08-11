//=============================================================================
// FILENAME : TestNumpyExport.cpp
//
// DESCRIPTION:
//     Test numpy export for nested (multi-rank) measurement data
//
// REVISION:
//  [mm/dd/yy]
//  [07/29/2026 nbale]
//=============================================================================

#include "CLGTest.h"

/**
 * Validate the NPY header of a file against the format spec.
 *
 * TestNumpyExportMultiRank below reads the header length and fseek()s straight
 * past the header dict, so it cannot see anything wrong *inside* the header.
 * That blind spot let a real bug ship: _SaveAsNumpyFileImpl appended the
 * terminating '\n' before the space padding, producing "{...}\n      ". The
 * header length stayed correct, the 64-byte alignment stayed correct, and the
 * data body started at the right offset -- so every value-based check passed
 * while numpy itself refused the file with "Cannot parse header".
 *
 * This checks the properties numpy actually enforces:
 *   - magic + version 1.0
 *   - (10 + headerLen) is a multiple of 64
 *   - the last header byte is '\n'
 *   - '\n' appears exactly once, at the end
 *   - the dict reports the expected descr and shape
 *   - the body holds exactly the expected number of elements, no trailing bytes
 *
 * Returns the number of errors found.
 */
static UINT _ValidateNpyHeader(
    const CCString& sFileName,
    const CCString& sExpectedDescr,
    const CCString& sExpectedShape,
    INT nExpectedElements,
    size_t uiElementBytes)
{
    UINT uiErrors = 0;

#if _CLG_WIN
    FILE* pFile = NULL;
    fopen_s(&pFile, sFileName.c_str(), "rb");
#else
    FILE* pFile = fopen(sFileName.c_str(), "rb");
#endif
    if (NULL == pFile)
    {
        appCrucial(_T("  %s: cannot open for header validation\n"), sFileName.c_str());
        return 1;
    }

    unsigned char magic[6];
    if (fread(magic, 1, 6, pFile) != 6 || memcmp(magic, "\x93NUMPY", 6) != 0)
    {
        appCrucial(_T("  %s: bad magic\n"), sFileName.c_str());
        fclose(pFile);
        return 1;
    }

    unsigned char major = 0, minor = 0;
    if (fread(&major, 1, 1, pFile) != 1 || fread(&minor, 1, 1, pFile) != 1)
    {
        appCrucial(_T("  %s: cannot read version\n"), sFileName.c_str());
        fclose(pFile);
        return 1;
    }
    if (1 != major || 0 != minor)
    {
        appCrucial(_T("  %s: expected version 1.0, got %d.%d\n"),
            sFileName.c_str(), major, minor);
        ++uiErrors;
    }

    uint16_t len16 = 0;
    if (fread(&len16, 2, 1, pFile) != 1)
    {
        appCrucial(_T("  %s: cannot read header length\n"), sFileName.c_str());
        fclose(pFile);
        return 1;
    }
    const size_t uiHeaderLen = len16;

    // numpy pads so that magic(6) + version(2) + len(2) + header is 64-aligned.
    if (0 != (10 + uiHeaderLen) % 64)
    {
        appCrucial(_T("  %s: header not 64-byte aligned (10 + %d = %d)\n"),
            sFileName.c_str(), static_cast<INT>(uiHeaderLen),
            static_cast<INT>(10 + uiHeaderLen));
        ++uiErrors;
    }

    std::string sHeader;
    sHeader.resize(uiHeaderLen);
    if (fread(&sHeader[0], 1, uiHeaderLen, pFile) != uiHeaderLen)
    {
        appCrucial(_T("  %s: header truncated\n"), sFileName.c_str());
        fclose(pFile);
        return uiErrors + 1;
    }

    // The spec requires the final byte to be '\n'. This is the check that the
    // padding-order bug violated.
    if (sHeader.empty() || '\n' != sHeader[uiHeaderLen - 1])
    {
        appCrucial(_T("  %s: last header byte is not '\\n' (got 0x%02x)\n"),
            sFileName.c_str(),
            sHeader.empty() ? 0 : static_cast<INT>(
                static_cast<unsigned char>(sHeader[uiHeaderLen - 1])));
        ++uiErrors;
    }

    // A '\n' anywhere else means padding was written after the terminator.
    const size_t uiFirstNewline = sHeader.find('\n');
    if (std::string::npos != uiFirstNewline && uiFirstNewline != uiHeaderLen - 1)
    {
        appCrucial(_T("  %s: '\\n' at offset %d, expected only at %d")
                   _T(" (padding written after the terminator?)\n"),
            sFileName.c_str(), static_cast<INT>(uiFirstNewline),
            static_cast<INT>(uiHeaderLen - 1));
        ++uiErrors;
    }

    // Everything between the dict and the terminator must be spaces.
    const size_t uiBrace = sHeader.rfind('}');
    if (std::string::npos == uiBrace)
    {
        appCrucial(_T("  %s: header has no closing '}'\n"), sFileName.c_str());
        ++uiErrors;
    }
    else
    {
        for (size_t i = uiBrace + 1; i + 1 < uiHeaderLen; ++i)
        {
            if (' ' != sHeader[i])
            {
                appCrucial(_T("  %s: non-space padding byte 0x%02x at offset %d\n"),
                    sFileName.c_str(),
                    static_cast<INT>(static_cast<unsigned char>(sHeader[i])),
                    static_cast<INT>(i));
                ++uiErrors;
                break;
            }
        }
    }

    // The dict must describe the array we asked for.
    CCString sDescrKey;
    sDescrKey.Format(_T("'descr': '%s'"), sExpectedDescr.c_str());
    if (std::string::npos == sHeader.find(sDescrKey.c_str()))
    {
        appCrucial(_T("  %s: expected %s in header\n"),
            sFileName.c_str(), sDescrKey.c_str());
        ++uiErrors;
    }
    if (std::string::npos == sHeader.find("'fortran_order': False"))
    {
        appCrucial(_T("  %s: expected 'fortran_order': False in header\n"),
            sFileName.c_str());
        ++uiErrors;
    }
    CCString sShapeKey;
    sShapeKey.Format(_T("'shape': (%s)"), sExpectedShape.c_str());
    if (std::string::npos == sHeader.find(sShapeKey.c_str()))
    {
        appCrucial(_T("  %s: expected %s in header\n"),
            sFileName.c_str(), sShapeKey.c_str());
        ++uiErrors;
    }

    // The body must hold exactly nExpectedElements, with nothing left over.
    const long lBodyStart = ftell(pFile);
    fseek(pFile, 0, SEEK_END);
    const long lFileSize = ftell(pFile);
    fclose(pFile);

    const long lBodyBytes = lFileSize - lBodyStart;
    const long lExpectedBytes =
        static_cast<long>(nExpectedElements) * static_cast<long>(uiElementBytes);
    if (lBodyBytes != lExpectedBytes)
    {
        appCrucial(_T("  %s: body is %d bytes, expected %d\n"),
            sFileName.c_str(), static_cast<INT>(lBodyBytes),
            static_cast<INT>(lExpectedBytes));
        ++uiErrors;
    }

    return uiErrors;
}

/**
 * Test nested TArray export to .npy files.
 *
 * Bug: SaveAsNumpyFile(TArray<TArray<TArray<T>>>) was writing TArray object
 * headers instead of leaf elements, producing garbage files. This test creates
 * a 3-level nested array (mimicking CMeasureWilsonLoop's [conf][r][t] structure),
 * exports it, reads it back, and validates the values.
 */
UINT TestNumpyExportMultiRank(CParameters& sParam)
{
    UINT uiErrors = 0;

    // Create a 3-level nested array: [Nconf][Nr][Nt]
    const INT Nconf = 5;
    const INT Nr = 10;
    const INT Nt = 8;

    TArray<TArray<TArray<CLGComplex>>> data;
    data.SetSize(Nconf);

    // Fill with known test pattern: data[c][r][t] = complex(c*100 + r*10 + t, -t)
    for (INT c = 0; c < Nconf; ++c)
    {
        data[c].SetSize(Nr);
        for (INT r = 0; r < Nr; ++r)
        {
            data[c][r].SetSize(Nt);
            for (INT t = 0; t < Nt; ++t)
            {
                Real re = static_cast<Real>(c * 100 + r * 10 + t);
                Real im = static_cast<Real>(-t);
                data[c][r][t] = _make_cuComplex(re, im);
            }
        }
    }

    // Export to .npy
    CCString sFileName = _T("test_multirank_export.npy");
    SaveAsNumpyFile(sFileName, data);

    // Read back and validate
    // The file should be: complex64 array of shape (5, 10, 8) = 400 elements
    // Each element: 8 bytes (float32 re + float32 im)
    // Header is ~128 bytes, body is 400*8 = 3200 bytes

#if _CLG_WIN
    FILE* pFile = NULL;
    fopen_s(&pFile, "test_multirank_export.npy", "rb");
#else
    FILE* pFile = fopen("test_multirank_export.npy", "rb");
#endif
    if (NULL == pFile)
    {
        appCrucial(_T("Failed to open %s for validation\n"), sFileName.c_str());
        return 1;
    }

    // Skip the numpy header (magic + version + header_len + header dict)
    // Format: 6 bytes magic "\x93NUMPY" + 1 byte major + 1 byte minor + 2/4 bytes header_len + header + data
    unsigned char magic[6];
    if (fread(magic, 1, 6, pFile) != 6 || memcmp(magic, "\x93NUMPY", 6) != 0)
    {
        appCrucial(_T("Invalid numpy magic or read error\n"));
        fclose(pFile);
        return 1;
    }

    unsigned char major, minor;
    if (fread(&major, 1, 1, pFile) != 1 || fread(&minor, 1, 1, pFile) != 1)
    {
        appCrucial(_T("Failed to read numpy version\n"));
        fclose(pFile);
        return 1;
    }

    UINT headerLen = 0;
    if (major == 1)
    {
        uint16_t len16;
        if (fread(&len16, 2, 1, pFile) != 1)
        {
            appCrucial(_T("Failed to read header length\n"));
            fclose(pFile);
            return 1;
        }
        headerLen = len16;
    }
    else if (major == 2 || major == 3)
    {
        uint32_t len32;
        if (fread(&len32, 4, 1, pFile) != 1)
        {
            appCrucial(_T("Failed to read header length\n"));
            fclose(pFile);
            return 1;
        }
        headerLen = len32;
    }
    else
    {
        appCrucial(_T("Unsupported numpy version %d.%d\n"), major, minor);
        fclose(pFile);
        return 1;
    }

    // Skip header dict
    fseek(pFile, headerLen, SEEK_CUR);

    // Read the data body
    const INT nTotal = Nconf * Nr * Nt;
    TArray<CLGComplex> readBack;
    readBack.SetSize(nTotal);
    size_t nRead = fread(readBack.GetData(), sizeof(CLGComplex), nTotal, pFile);
    fclose(pFile);

    if (static_cast<INT>(nRead) != nTotal)
    {
        appCrucial(_T("Read %d elements, expected %d\n"), static_cast<INT>(nRead), nTotal);
        ++uiErrors;
    }

    // Validate: readBack is row-major flattened [c][r][t]
    INT idx = 0;
    for (INT c = 0; c < Nconf; ++c)
    {
        for (INT r = 0; r < Nr; ++r)
        {
            for (INT t = 0; t < Nt; ++t)
            {
                Real expectedRe = static_cast<Real>(c * 100 + r * 10 + t);
                Real expectedIm = static_cast<Real>(-t);
                Real actualRe = readBack[idx].x;
                Real actualIm = readBack[idx].y;

                if (appAbs(actualRe - expectedRe) > F(1e-5) || appAbs(actualIm - expectedIm) > F(1e-5))
                {
                    if (uiErrors == 0)
                    {
                        appCrucial(_T("Mismatch at [%d][%d][%d]: expected (%f,%f) got (%f,%f)\n"),
                            c, r, t, expectedRe, expectedIm, actualRe, actualIm);
                    }
                    ++uiErrors;
                }
                ++idx;
            }
        }
    }

    if (uiErrors == 0)
    {
        appGeneral(_T("Multi-rank numpy export test PASSED: all %d elements validated\n"), nTotal);
    }
    else
    {
        appCrucial(_T("Multi-rank numpy export test FAILED: %d mismatches out of %d elements\n"),
            uiErrors, nTotal);
    }

    return uiErrors > 0 ? 1 : 0;
}

__REGIST_TEST(TestNumpyExportMultiRank, Common, TestNumpyExportMultiRank, NumpyExportMultiRank);

/**
 * Validate the NPY header for ranks 1 through 4.
 *
 * All SaveAsNumpyFile overloads funnel into _SaveAsNumpyFileImpl, so one bad
 * header affects every rank -- but the shape string differs per rank (rank 1
 * needs the trailing comma of a Python 1-tuple), and the padding lands at a
 * different offset for each, so each rank is worth exercising separately.
 *
 * Both real element types are covered: DOUBLE ('<f8') and CLGComplex ('<c8'),
 * since the descr and element size come from separate code paths.
 */
UINT TestNumpyExportHeader(CParameters& sParam)
{
    UINT uiErrors = 0;

    //---------------------------------------------------------------- rank 1
    {
        const INT N1 = 7;
        TArray<DOUBLE> a1;
        a1.SetSize(N1);
        for (INT i = 0; i < N1; ++i)
        {
            a1[i] = static_cast<DOUBLE>(i) * 1.5;
        }

        const CCString sFile = _T("test_npy_rank1.npy");
        SaveAsNumpyFile(sFile, a1);
        // A rank-1 shape must be a Python 1-tuple: "(7,)".
        const UINT uiRankErrors =
            _ValidateNpyHeader(sFile, _T("<f8"), _T("7,"), N1, sizeof(DOUBLE));
        appGeneral(_T("  rank 1 (7,) <f8    : %s\n"),
            0 == uiRankErrors ? _T("OK") : _T("FAILED"));
        uiErrors += uiRankErrors;
    }

    //---------------------------------------------------------------- rank 2
    {
        const INT N1 = 4, N2 = 6;
        TArray<TArray<CLGComplex>> a2;
        a2.SetSize(N1);
        for (INT i = 0; i < N1; ++i)
        {
            a2[i].SetSize(N2);
            for (INT j = 0; j < N2; ++j)
            {
                a2[i][j] = _make_cuComplex(
                    static_cast<Real>(i * 10 + j), static_cast<Real>(-j));
            }
        }

        const CCString sFile = _T("test_npy_rank2.npy");
        SaveAsNumpyFile(sFile, a2);
        const UINT uiRankErrors = _ValidateNpyHeader(
            sFile, _T("<c8"), _T("4, 6"), N1 * N2, sizeof(CLGComplex));
        appGeneral(_T("  rank 2 (4, 6) <c8  : %s\n"),
            0 == uiRankErrors ? _T("OK") : _T("FAILED"));
        uiErrors += uiRankErrors;
    }

    //---------------------------------------------------------------- rank 3
    {
        const INT N1 = 3, N2 = 5, N3 = 4;
        TArray<TArray<TArray<DOUBLE>>> a3;
        a3.SetSize(N1);
        for (INT i = 0; i < N1; ++i)
        {
            a3[i].SetSize(N2);
            for (INT j = 0; j < N2; ++j)
            {
                a3[i][j].SetSize(N3);
                for (INT k = 0; k < N3; ++k)
                {
                    a3[i][j][k] = static_cast<DOUBLE>(i * 100 + j * 10 + k);
                }
            }
        }

        const CCString sFile = _T("test_npy_rank3.npy");
        SaveAsNumpyFile(sFile, a3);
        const UINT uiRankErrors = _ValidateNpyHeader(
            sFile, _T("<f8"), _T("3, 5, 4"), N1 * N2 * N3, sizeof(DOUBLE));
        appGeneral(_T("  rank 3 (3, 5, 4) <f8: %s\n"),
            0 == uiRankErrors ? _T("OK") : _T("FAILED"));
        uiErrors += uiRankErrors;
    }

    //---------------------------------------------------------------- rank 4
    {
        const INT N1 = 2, N2 = 3, N3 = 4, N4 = 5;
        TArray<TArray<TArray<TArray<CLGComplex>>>> a4;
        a4.SetSize(N1);
        for (INT i = 0; i < N1; ++i)
        {
            a4[i].SetSize(N2);
            for (INT j = 0; j < N2; ++j)
            {
                a4[i][j].SetSize(N3);
                for (INT k = 0; k < N3; ++k)
                {
                    a4[i][j][k].SetSize(N4);
                    for (INT l = 0; l < N4; ++l)
                    {
                        a4[i][j][k][l] = _make_cuComplex(
                            static_cast<Real>(i * 1000 + j * 100 + k * 10 + l),
                            static_cast<Real>(l));
                    }
                }
            }
        }

        const CCString sFile = _T("test_npy_rank4.npy");
        SaveAsNumpyFile(sFile, a4);
        const UINT uiRankErrors = _ValidateNpyHeader(
            sFile, _T("<c8"), _T("2, 3, 4, 5"),
            N1 * N2 * N3 * N4, sizeof(CLGComplex));
        appGeneral(_T("  rank 4 (2, 3, 4, 5) <c8: %s\n"),
            0 == uiRankErrors ? _T("OK") : _T("FAILED"));
        uiErrors += uiRankErrors;
    }

    if (0 == uiErrors)
    {
        appGeneral(_T("Numpy header test PASSED: ranks 1-4 conform to the NPY spec\n"));
    }
    else
    {
        appCrucial(_T("Numpy header test FAILED: %d header problems\n"), uiErrors);
    }

    return uiErrors > 0 ? 1 : 0;
}

__REGIST_TEST(TestNumpyExportHeader, Common, TestNumpyExportHeader, NumpyExportHeader);

//=============================================================================
// END OF FILE
//=============================================================================
