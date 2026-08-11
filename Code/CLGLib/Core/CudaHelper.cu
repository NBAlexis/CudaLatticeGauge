//=============================================================================
// FILENAME : CudaHelper.h
// 
// DESCRIPTION:
// This is the file for CUDA Testing usage
//
// REVISION:
//  [mm/dd/yy]
//  [12/3/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

//below are added for test
#include "Measurement/CMeasurePlaqutteEnergy.h"
#include "Data/Field/Staggered/CFieldFermionKSHISQ.h"
#include "GaugeSmearing/CGaugeSmearingHISQ.h"
#include "Tools/Math/SU3_12.h"

__BEGIN_NAMESPACE

__device__ __constant__ UINT _constIntegers[kContentLength];
__device__ __constant__ INT _constSignedIntegers[kContentLength];
__device__ __constant__ Real _constFloats[kContentLength];
__device__ __constant__ CRandom* __r;
__device__ __constant__ CIndexData* __idx;
//why we create Dirac gamma matrix?
//__constant__ gammaMatrix __diracGamma[EGM_MAX];
__device__ __constant__ gammaMatrix __chiralGamma[EGM_MAX];
__device__ __constant__ deviceSU3 __SU3Generators[9];

__device__ __constant__ CField* __fieldPointers[kMaxFieldCount];
__device__ __constant__ CFieldBoundaryParent* __boundaryFieldPointers[kMaxFieldCount];

//__device__ __constant__ constexpr SCHAR _plaq_idx[6][2] = {
//    {1, 2},
//    {1, 3},
//    {1, 4},
//    {2, 3},
//    {2, 4},
//    {3, 4},
//};

#pragma region Kernels

/**
* The construction is on device
*/
__global__ void 
_CLG_LAUNCH_BOUND_SINGLE
_kernelCreateMatrix(/*gammaMatrix* pDirac,*/ gammaMatrix* pChiral, deviceSU3* pGenerator)
{
    gammaMatrixSet::CreateGammaMatrix(EGMS_Chiral, pChiral);
    #pragma unroll
    for (int i = 0; i < 9; ++i)
    {
        pGenerator[i] = deviceSU3::makeSU3Generator(i);
    }
}

__global__ void 
_CLG_LAUNCH_BOUND_SINGLE
_kernelDebugFunction()
{
    //ULONGLONG test[4];
    //test[0] = 0;
    //test[1] = 0;
    //test[2] = 0;
    //test[3] = 0;
    //((SIndex*)&test[2])->m_uiSiteIndex = 1;
    //((SIndex*)&test[2])->m_byDir = 1;

    //printf("testres %lld\n", test[0]);
    //printf("testres %lld\n", test[1]);
    //printf("testres %lld\n", test[2]);
    //printf("testres %lld\n", test[3]);

    //printf("testres %d\n", ((SIndex*)&test[2])->m_uiSiteIndex);
    //printf("testres %lld\n", ((SIndex*)&test[2])->m_ullData);
    //__idx->DebugPrintWalkingTable();
    
    //for (UINT i = 0; i < EGM_MAX; ++i)
    //{
    //    printf("gammamatrix is %d\n", i);
    //    __chiralGamma[i].Print();
    //}
    //printf("here1?\n");
    //for (INT i = 0; i < 15; ++i)
    //{
    //    deviceSU4 g = deviceSU4::makeSUNGenerator(i);
    //    g.DebugPrint("g");
    //}

    //deviceSU2 a = deviceSU2::makeSU2Random(0);
    //a.DebugPrint("a");
    //printf("sizeof(SU4) = %d\n", static_cast<INT>(sizeof(deviceSU4)));
    //printf("sizeof(SU8) = %d\n", static_cast<INT>(sizeof(deviceSU8)));

    printf("\n===== deviceSU3_12 Comprehensive Tests =====\n");
    printf("sizeof(deviceSU3) = %d\n", static_cast<INT>(sizeof(deviceSU3)));
    printf("sizeof(deviceSU3_12) = %d\n", static_cast<INT>(sizeof(deviceSU3_12)));

    deviceSU3 su3a = deviceSU3::makeSU3Random(0);
    deviceSU3 su3b = deviceSU3::makeSU3Random(1);
    deviceSU3_12 a12(su3a);
    deviceSU3_12 b12(su3b);

    // ====================================================================
    // SECTION 1: Constructors and Factory
    // ====================================================================
    printf("\n=== SECTION 1: Constructors and Factory ===\n");

    // T01: default constructor + copy constructor
    printf("\n--- T01: default + copy constructor ---\n");
    deviceSU3_12 copy12(a12);
    copy12.toSU3().DebugPrint("copy_of_a12");

    // T02: makeSU3_12Zero / Zero
    printf("\n--- T02: makeSU3_12Zero / Zero ---\n");
    deviceSU3_12 z12 = deviceSU3_12::makeSU3_12Zero();
    z12.toSU3().DebugPrint("makeZero");
    deviceSU3_12 z12b; z12b.Zero();
    z12b.toSU3().DebugPrint("Zero_mut");

    // T03: makeSU3_12Id / Id
    printf("\n--- T03: makeSU3_12Id / Id ---\n");
    deviceSU3_12 id12 = deviceSU3_12::makeSU3_12Id();
    id12.toSU3().DebugPrint("makeId");
    deviceSU3_12 id12b; id12b.Id();
    id12b.toSU3().DebugPrint("Id_mut");

    // T04: makeSU3_12Random
    printf("\n--- T04: makeSU3_12Random ---\n");
    deviceSU3_12 r12 = deviceSU3_12::makeSU3_12Random(3);
    deviceSU3 r3 = r12.toSU3();
    r3.DebugPrint("Random_12");
    printf("det=%f\n", r3.Det().x);

    // T05: makeSU3_12RandomGenerator
    printf("\n--- T05: makeSU3_12RandomGenerator ---\n");
    deviceSU3_12 rg12 = deviceSU3_12::makeSU3_12RandomGenerator(4);
    rg12.toSU3().DebugPrint("RandomGen_12");

    // T06: makeSU3_12Generator (all 8)
    printf("\n--- T06: makeSU3_12Generator (0..7) ---\n");
    for (UINT ig = 0; ig < 8; ++ig)
    {
        deviceSU3_12 g12 = deviceSU3_12::makeSU3_12Generator(ig);
        deviceSU3 g3 = deviceSU3::makeSU3Generator(ig);
        printf("gen%u: SU3[0]=(%f,%f) SU3_12[0]=(%f,%f)\n", ig,
            g3.m_me[0].x, g3.m_me[0].y, g12.m_me[0].x, g12.m_me[0].y);
    }

    // T07: makeSU3_12TA
    printf("\n--- T07: makeSU3_12TA ---\n");
    CLGComplex m12t = _make_cuComplex(F(0.1), F(0.2));
    CLGComplex m13t = _make_cuComplex(F(0.3), F(0.4));
    CLGComplex m23t = _make_cuComplex(F(0.5), F(0.6));
    deviceSU3_12 ta12 = deviceSU3_12::makeSU3_12TA(m12t, m13t, m23t, F(0.7), F(0.8));
    deviceSU3 ta3 = deviceSU3::makeSU3TA(m12t, m13t, m23t, F(0.7), F(0.8));
    printf("TA SU3[0]=(%f,%f) SU3_12[0]=(%f,%f)\n", ta3.m_me[0].x, ta3.m_me[0].y, ta12.m_me[0].x, ta12.m_me[0].y);

    // T08: makeSU3_12ContractV
    printf("\n--- T08: makeSU3_12ContractV ---\n");
    deviceSU3Vector lv, rv;
    lv.m_ve[0] = _make_cuComplex(F(0.1), F(0.2));
    lv.m_ve[1] = _make_cuComplex(F(0.3), F(0.4));
    lv.m_ve[2] = _make_cuComplex(F(0.5), F(0.6));
    rv.m_ve[0] = _make_cuComplex(F(0.7), F(0.8));
    rv.m_ve[1] = _make_cuComplex(F(0.9), F(0.1));
    rv.m_ve[2] = _make_cuComplex(F(0.2), F(0.3));
    deviceSU3_12 cv12 = deviceSU3_12::makeSU3_12ContractV(lv, rv);
    deviceSU3 cv3 = deviceSU3::makeSU3ContractV(lv, rv);
    printf("ContractV SU3[0]=(%f,%f) SU3_12[0]=(%f,%f)\n", cv3.m_me[0].x, cv3.m_me[0].y, cv12.m_me[0].x, cv12.m_me[0].y);

    // T08b: makeSU3_12Contract (WilsonVector version)
    printf("\n--- T08b: makeSU3_12Contract (WilsonVector) ---\n");
    deviceWilsonVectorSU3 wlv, wrv;
    for (int s = 0; s < 4; ++s) {
        wlv.m_d[s].m_ve[0] = _make_cuComplex(F(0.1) * (s + 1), F(0.2) * (s + 1));
        wlv.m_d[s].m_ve[1] = _make_cuComplex(F(0.3) * (s + 1), F(0.4) * (s + 1));
        wlv.m_d[s].m_ve[2] = _make_cuComplex(F(0.5) * (s + 1), F(0.6) * (s + 1));
        wrv.m_d[s].m_ve[0] = _make_cuComplex(F(0.7) * (s + 1), F(0.8) * (s + 1));
        wrv.m_d[s].m_ve[1] = _make_cuComplex(F(0.9) * (s + 1), F(0.1) * (s + 1));
        wrv.m_d[s].m_ve[2] = _make_cuComplex(F(0.2) * (s + 1), F(0.3) * (s + 1));
    }
    deviceSU3_12 wc12 = deviceSU3_12::makeSU3_12Contract(wlv, wrv);
    deviceSU3 wc3 = deviceSU3::makeSU3Contract(wlv, wrv);
    printf("Contract SU3[0]=(%f,%f) SU3_12[0]=(%f,%f)\n", wc3.m_me[0].x, wc3.m_me[0].y, wc12.m_me[0].x, wc12.m_me[0].y);

    // T09: makeSU3_12SumGenerator
    printf("\n--- T09: makeSU3_12SumGenerator ---\n");
    deviceSU3_12 sg12 = deviceSU3_12::makeSU3_12SumGenerator(F(0.5));
    deviceSU3 sg3 = deviceSU3::makeSU3SumGenerator(F(0.5));
    printf("SumGen SU3[0]=(%f,%f) SU3_12[0]=(%f,%f)\n", sg3.m_me[0].x, sg3.m_me[0].y, sg12.m_me[0].x, sg12.m_me[0].y);

    // ====================================================================
    // SECTION 2: Conversion
    // ====================================================================
    printf("\n=== SECTION 2: Conversion ===\n");

    // T10: compress/expand round-trip
    printf("\n--- T10: compress/expand round-trip ---\n");
    deviceSU3 su3a_back = a12.toSU3();
    su3a.DebugPrint("original");
    su3a_back.DebugPrint("roundtrip");

    // T11: cross product reconstruction
    printf("\n--- T11: cross product reconstruction ---\n");
    CLGComplex c6, c7, c8;
    a12._reconstructCol2(c6, c7, c8);
    printf("col2 orig: (%f,%f) (%f,%f) (%f,%f)\n",
        su3a.m_me[6].x, su3a.m_me[6].y, su3a.m_me[7].x, su3a.m_me[7].y, su3a.m_me[8].x, su3a.m_me[8].y);
    printf("col2 reconst: (%f,%f) (%f,%f) (%f,%f)\n",
        c6.x, c6.y, c7.x, c7.y, c8.x, c8.y);

    // T12: compress / expand static methods
    printf("\n--- T12: compress/expand static ---\n");
    deviceSU3_12 comp12 = deviceSU3_12::compress(su3b);
    deviceSU3 exp3 = deviceSU3_12::expand(comp12);
    su3b.DebugPrint("original_b");
    exp3.DebugPrint("expand_compress_b");

    // T13: free function conversions
    printf("\n--- T13: free function conversions ---\n");
    deviceSU3_12 conv12 = deviceSU3_12::deviceSU3_to_deviceSU3_12(su3a);
    deviceSU3 conv3 = deviceSU3_12::deviceSU3_12_to_deviceSU3(conv12);
    su3a.DebugPrint("original");
    conv3.DebugPrint("converted_back");

    // ====================================================================
    // SECTION 3: SU3-preserving operations
    // ====================================================================
    printf("\n=== SECTION 3: SU3-preserving operations ===\n");

    // T14: MulC (return version)
    printf("\n--- T14: MulC ---\n");
    deviceSU3_12 ab12 = a12.MulC(b12);
    deviceSU3 su3ab = su3a; su3ab.Mul(su3b);
    su3ab.DebugPrint("SU3_Mul");
    ab12.toSU3().DebugPrint("SU3_12_MulC");

    // T15: Mul (in-place)
    printf("\n--- T15: Mul (in-place) ---\n");
    deviceSU3_12 a12ip(a12);
    a12ip.Mul(b12);
    a12ip.toSU3().DebugPrint("SU3_12_Mul_ip");

    // T16: MulDaggerC
    printf("\n--- T16: MulDaggerC ---\n");
    deviceSU3_12 abd12 = a12.MulDaggerC(b12);
    deviceSU3 su3abd = su3a; su3abd.MulDagger(su3b);
    su3abd.DebugPrint("SU3_MulDagger");
    abd12.toSU3().DebugPrint("SU3_12_MulDaggerC");

    // T17: MulDagger (in-place)
    printf("\n--- T17: MulDagger (in-place) ---\n");
    deviceSU3_12 a12ip2(a12);
    a12ip2.MulDagger(b12);
    a12ip2.toSU3().DebugPrint("SU3_12_MulDagger_ip");

    // T18: DaggerMulC
    printf("\n--- T18: DaggerMulC ---\n");
    deviceSU3_12 dab12 = a12.DaggerMulC(b12);
    deviceSU3 su3dab = su3a; su3dab.DaggerMul(su3b);
    su3dab.DebugPrint("SU3_DaggerMul");
    dab12.toSU3().DebugPrint("SU3_12_DaggerMulC");

    // T19: DaggerMul (in-place)
    printf("\n--- T19: DaggerMul (in-place) ---\n");
    deviceSU3_12 a12ip3(a12);
    a12ip3.DaggerMul(b12);
    a12ip3.toSU3().DebugPrint("SU3_12_DaggerMul_ip");

    // T20: DaggerC
    printf("\n--- T20: DaggerC ---\n");
    deviceSU3_12 ad12 = a12.DaggerC();
    deviceSU3 su3ad = su3a; su3ad.Dagger();
    su3ad.DebugPrint("SU3_Dagger");
    ad12.toSU3().DebugPrint("SU3_12_DaggerC");

    // T21: Dagger (in-place)
    printf("\n--- T21: Dagger (in-place) ---\n");
    deviceSU3_12 a12ip4(a12);
    a12ip4.Dagger();
    a12ip4.toSU3().DebugPrint("SU3_12_Dagger_ip");

    // T22: MulVector
    printf("\n--- T22: MulVector ---\n");
    deviceSU3Vector v;
    v.m_ve[0] = _make_cuComplex(F(0.1), F(0.2));
    v.m_ve[1] = _make_cuComplex(F(0.3), F(0.4));
    v.m_ve[2] = _make_cuComplex(F(0.5), F(0.6));
    deviceSU3Vector vres3 = su3a.MulVector(v);
    deviceSU3Vector vres12 = a12.MulVector(v);
    printf("SU3: (%f,%f) (%f,%f) (%f,%f)\n",
        vres3.m_ve[0].x, vres3.m_ve[0].y, vres3.m_ve[1].x, vres3.m_ve[1].y, vres3.m_ve[2].x, vres3.m_ve[2].y);
    printf("SU3_12: (%f,%f) (%f,%f) (%f,%f)\n",
        vres12.m_ve[0].x, vres12.m_ve[0].y, vres12.m_ve[1].x, vres12.m_ve[1].y, vres12.m_ve[2].x, vres12.m_ve[2].y);

    // T23: DagMulVector
    printf("\n--- T23: DagMulVector ---\n");
    deviceSU3Vector dvres3 = su3a.DagMulVector(v);
    deviceSU3Vector dvres12 = a12.DagMulVector(v);
    printf("SU3: (%f,%f) (%f,%f) (%f,%f)\n",
        dvres3.m_ve[0].x, dvres3.m_ve[0].y, dvres3.m_ve[1].x, dvres3.m_ve[1].y, dvres3.m_ve[2].x, dvres3.m_ve[2].y);
    printf("SU3_12: (%f,%f) (%f,%f) (%f,%f)\n",
        dvres12.m_ve[0].x, dvres12.m_ve[0].y, dvres12.m_ve[1].x, dvres12.m_ve[1].y, dvres12.m_ve[2].x, dvres12.m_ve[2].y);

    // T24: MulWilsonVector / DagMulWilsonVector
    printf("\n--- T24: MulWilsonVector / DagMulWilsonVector ---\n");
    deviceWilsonVectorSU3 wv;
    for (INT d = 0; d < 4; ++d)
    {
        wv.m_d[d].m_ve[0] = _make_cuComplex(F(0.1) * (d+1), F(0.2) * (d+1));
        wv.m_d[d].m_ve[1] = _make_cuComplex(F(0.3) * (d+1), F(0.4) * (d+1));
        wv.m_d[d].m_ve[2] = _make_cuComplex(F(0.5) * (d+1), F(0.6) * (d+1));
    }
    deviceWilsonVectorSU3 wvres3 = su3a.MulWilsonVector(wv);
    deviceWilsonVectorSU3 wvres12 = a12.MulWilsonVector(wv);
    printf("SU3 W[0]: (%f,%f) (%f,%f) (%f,%f)\n",
        wvres3.m_d[0].m_ve[0].x, wvres3.m_d[0].m_ve[0].y, wvres3.m_d[0].m_ve[1].x, wvres3.m_d[0].m_ve[1].y, wvres3.m_d[0].m_ve[2].x, wvres3.m_d[0].m_ve[2].y);
    printf("SU3_12 W[0]: (%f,%f) (%f,%f) (%f,%f)\n",
        wvres12.m_d[0].m_ve[0].x, wvres12.m_d[0].m_ve[0].y, wvres12.m_d[0].m_ve[1].x, wvres12.m_d[0].m_ve[1].y, wvres12.m_d[0].m_ve[2].x, wvres12.m_d[0].m_ve[2].y);
    deviceWilsonVectorSU3 dwvres3 = su3a.DagMulWilsonVector(wv);
    deviceWilsonVectorSU3 dwvres12 = a12.DagMulWilsonVector(wv);
    printf("SU3 DW[0]: (%f,%f) (%f,%f) (%f,%f)\n",
        dwvres3.m_d[0].m_ve[0].x, dwvres3.m_d[0].m_ve[0].y, dwvres3.m_d[0].m_ve[1].x, dwvres3.m_d[0].m_ve[1].y, dwvres3.m_d[0].m_ve[2].x, dwvres3.m_d[0].m_ve[2].y);
    printf("SU3_12 DW[0]: (%f,%f) (%f,%f) (%f,%f)\n",
        dwvres12.m_d[0].m_ve[0].x, dwvres12.m_d[0].m_ve[0].y, dwvres12.m_d[0].m_ve[1].x, dwvres12.m_d[0].m_ve[1].y, dwvres12.m_d[0].m_ve[2].x, dwvres12.m_d[0].m_ve[2].y);

    // T25: Transpose
    printf("\n--- T25: Transpose ---\n");
    deviceSU3 tr12 = a12.Transpose();
    deviceSU3 manual_tr;
    manual_tr.m_me[0] = su3a.m_me[0]; manual_tr.m_me[1] = su3a.m_me[3]; manual_tr.m_me[2] = su3a.m_me[6];
    manual_tr.m_me[3] = su3a.m_me[1]; manual_tr.m_me[4] = su3a.m_me[4]; manual_tr.m_me[5] = su3a.m_me[7];
    manual_tr.m_me[6] = su3a.m_me[2]; manual_tr.m_me[7] = su3a.m_me[5]; manual_tr.m_me[8] = su3a.m_me[8];
    manual_tr.DebugPrint("manual_transpose");
    tr12.DebugPrint("SU3_12_Transpose");

    // T26: Inverse (= Dagger for SU(3))
    printf("\n--- T26: Inverse ---\n");
    deviceSU3_12 ai12 = a12.Inverse();
    deviceSU3 su3ai = su3a.Inverse();
    su3ai.DebugPrint("SU3_Inverse");
    ai12.toSU3().DebugPrint("SU3_12_Inverse");

    // T27: Projection no-ops
    printf("\n--- T27: Projection no-ops ---\n");
    deviceSU3_12 p12(a12);
    p12.Norm();
    p12.SU3Proj();
    p12.U3Proj();
    p12.Proj();
    p12.CabbiboMarinariProj();
    p12.toSU3().DebugPrint("after_projections");
    a12.toSU3().DebugPrint("original_a12");

    // ====================================================================
    // SECTION 4: Scalar queries
    // ====================================================================
    printf("\n=== SECTION 4: Scalar queries ===\n");

    // T28: Tr, ReTr, ImTr
    printf("\n--- T28: Tr, ReTr, ImTr ---\n");
    CLGComplex tr3 = su3a.Tr();
    CLGComplex tr12v = a12.Tr();
    printf("Tr SU3=(%f,%f) SU3_12=(%f,%f)\n", tr3.x, tr3.y, tr12v.x, tr12v.y);
    printf("ReTr SU3=%f SU3_12=%f\n", su3a.ReTr(), a12.ReTr());
    printf("ImTr SU3=%f SU3_12=%f\n", su3a.ImTr(), a12.ImTr());

    // T29: Det, DoubleDet
    printf("\n--- T29: Det, DoubleDet ---\n");
    printf("Det SU3=(%f,%f) SU3_12=(%f,%f)\n", su3a.Det().x, su3a.Det().y, a12.Det().x, a12.Det().y);
    cuDoubleComplex ddet3 = su3a.DoubleDet();
    cuDoubleComplex ddet12 = a12.DoubleDet();
    printf("DoubleDet SU3=(%f,%f) SU3_12=(%f,%f)\n", ddet3.x, ddet3.y, ddet12.x, ddet12.y);

    // T30: TrIm
    printf("\n--- T30: TrIm ---\n");
    Real trim3 = deviceSU3::TrIm(su3a, su3b);
    Real trim12 = deviceSU3_12::TrIm(a12, b12);
    printf("TrIm SU3=%f SU3_12=%f\n", trim3, trim12);

    // ====================================================================
    // SECTION 5: Non-SU(3) ops returning deviceSU3
    // ====================================================================
    printf("\n=== SECTION 5: Non-SU(3) ops returning deviceSU3 ===\n");

    // T31: AddC
    printf("\n--- T31: AddC ---\n");
    deviceSU3 add3 = su3a; add3.Add(su3b);
    deviceSU3 add12 = a12.AddC(b12);
    add3.DebugPrint("SU3_Add");
    add12.DebugPrint("SU3_12_AddC");

    // T32: AddDaggerC
    printf("\n--- T32: AddDaggerC ---\n");
    deviceSU3 addd3 = su3a; addd3.AddDagger(su3b);
    deviceSU3 addd12 = a12.AddDaggerC(b12);
    addd3.DebugPrint("SU3_AddDagger");
    addd12.DebugPrint("SU3_12_AddDaggerC");

    // T33: SubC
    printf("\n--- T33: SubC ---\n");
    deviceSU3 sub3 = su3a; sub3.Sub(su3b);
    deviceSU3 sub12 = a12.SubC(b12);
    sub3.DebugPrint("SU3_Sub");
    sub12.DebugPrint("SU3_12_SubC");

    // T34: SubDaggerC
    printf("\n--- T34: SubDaggerC ---\n");
    deviceSU3 subd3 = su3a; subd3.SubDagger(su3b);
    deviceSU3 subd12 = a12.SubDaggerC(b12);
    subd3.DebugPrint("SU3_SubDagger");
    subd12.DebugPrint("SU3_12_SubDaggerC");

    // T35: AddRealC
    printf("\n--- T35: AddRealC ---\n");
    deviceSU3 adr3 = su3a; adr3.AddReal(F(1.5));
    deviceSU3 adr12 = a12.AddRealC(F(1.5));
    adr3.DebugPrint("SU3_AddReal");
    adr12.DebugPrint("SU3_12_AddRealC");

    // T36: AddCompC
    printf("\n--- T36: AddCompC ---\n");
    CLGComplex ac = _make_cuComplex(F(0.5), F(0.3));
    deviceSU3 adc3 = su3a; adc3.AddComp(ac);
    deviceSU3 adc12 = a12.AddCompC(ac);
    adc3.DebugPrint("SU3_AddComp");
    adc12.DebugPrint("SU3_12_AddCompC");

    // T37: SubRealC
    printf("\n--- T37: SubRealC ---\n");
    deviceSU3 sur3 = su3a; sur3.SubReal(F(0.7));
    deviceSU3 sur12 = a12.SubRealC(F(0.7));
    sur3.DebugPrint("SU3_SubReal");
    sur12.DebugPrint("SU3_12_SubRealC");

    // T38: SubCompC
    printf("\n--- T38: SubCompC ---\n");
    deviceSU3 suc3 = su3a; suc3.SubComp(ac);
    deviceSU3 suc12 = a12.SubCompC(ac);
    suc3.DebugPrint("SU3_SubComp");
    suc12.DebugPrint("SU3_12_SubCompC");

    // T39: MulRealC
    printf("\n--- T39: MulRealC ---\n");
    deviceSU3 mr3 = su3a; mr3.MulReal(F(2.0));
    deviceSU3 mr12 = a12.MulRealC(F(2.0));
    mr3.DebugPrint("SU3_MulReal");
    mr12.DebugPrint("SU3_12_MulRealC");

    // T40: MulCompC
    printf("\n--- T40: MulCompC ---\n");
    CLGComplex mc = _make_cuComplex(F(0.5), F(0.3));
    deviceSU3 mco3 = su3a; mco3.MulComp(mc);
    deviceSU3 mco12 = a12.MulCompC(mc);
    mco3.DebugPrint("SU3_MulComp");
    mco12.DebugPrint("SU3_12_MulCompC");

    // T41: DivCompC
    printf("\n--- T41: DivCompC ---\n");
    CLGComplex dc = _make_cuComplex(F(2.0), F(1.0));
    deviceSU3 div3 = su3a; div3.DivComp(dc);
    deviceSU3 div12 = a12.DivCompC(dc);
    div3.DebugPrint("SU3_DivComp");
    div12.DebugPrint("SU3_12_DivCompC");

    // T42: OppositeC
    printf("\n--- T42: OppositeC ---\n");
    deviceSU3 opp3 = su3a; opp3.Opposite();
    deviceSU3 opp12 = a12.OppositeC();
    opp3.DebugPrint("SU3_Opposite");
    opp12.DebugPrint("SU3_12_OppositeC");

    // T43: ReC
    printf("\n--- T43: ReC ---\n");
    deviceSU3 re3 = su3a; re3.Re();
    deviceSU3 re12 = a12.ReC();
    re3.DebugPrint("SU3_Re");
    re12.DebugPrint("SU3_12_ReC");

    // T44: ImC
    printf("\n--- T44: ImC ---\n");
    deviceSU3 im3 = su3a; im3.Im();
    deviceSU3 im12 = a12.ImC();
    im3.DebugPrint("SU3_Im");
    im12.DebugPrint("SU3_12_ImC");

    // T45: TaC
    printf("\n--- T45: TaC ---\n");
    deviceSU3 tac3 = su3a; tac3.Ta();
    deviceSU3 tac12 = a12.TaC();
    tac3.DebugPrint("SU3_Ta");
    tac12.DebugPrint("SU3_12_TaC");

    // T46: ThC
    printf("\n--- T46: ThC ---\n");
    deviceSU3 th3 = su3a; th3.Th();
    deviceSU3 th12 = a12.ThC();
    th3.DebugPrint("SU3_Th");
    th12.DebugPrint("SU3_12_ThC");

    // T47: iIm2C
    printf("\n--- T47: iIm2C ---\n");
    deviceSU3 iim3 = su3a; iim3.iIm2();
    deviceSU3 iim12 = a12.iIm2C();
    iim3.DebugPrint("SU3_iIm2");
    iim12.DebugPrint("SU3_12_iIm2C");

    // T48: Re2C
    printf("\n--- T48: Re2C ---\n");
    deviceSU3 re23 = su3a; re23.Re2();
    deviceSU3 re212 = a12.Re2C();
    re23.DebugPrint("SU3_Re2");
    re212.DebugPrint("SU3_12_Re2C");

    // T49: Im2C
    printf("\n--- T49: Im2C ---\n");
    deviceSU3 im23 = su3a.Im2C();
    deviceSU3 im212 = a12.Im2C();
    im23.DebugPrint("SU3_Im2C");
    im212.DebugPrint("SU3_12_Im2C");

    // ====================================================================
    // SECTION 6: Exp/Log/Power (return deviceSU3)
    // ====================================================================
    printf("\n=== SECTION 6: Exp/Log/Power ===\n");

    // T50: QuickExp round-trip
    printf("\n--- T50: QuickExp round-trip ---\n");
    deviceSU3 gen = deviceSU3::makeSU3RandomGenerator(2);
    deviceSU3 exp_su3 = gen.QuickExp(F(0.5));
    deviceSU3_12 exp12(exp_su3);
    deviceSU3 exp_back = exp12.toSU3();
    exp_su3.DebugPrint("SU3_QuickExp");
    exp_back.DebugPrint("SU3_12_roundtrip");

    // T51: StrictExpTA round-trip
    printf("\n--- T51: StrictExpTA round-trip ---\n");
    deviceSU3 sexp3 = gen.StrictExpTA(F(0.5));
    deviceSU3_12 sexp12(sexp3);
    deviceSU3 sexp_back = sexp12.toSU3();
    sexp3.DebugPrint("SU3_StrictExpTA");
    sexp_back.DebugPrint("SU3_12_roundtrip");

    // T52: Exp
    printf("\n--- T52: Exp ---\n");
    CLGComplex ea = _make_cuComplex(F(0.5), F(0.1));
    deviceSU3 e3 = su3a.Exp(ea, 10);
    deviceSU3 e12 = a12.Exp(ea, 10);
    e3.DebugPrint("SU3_Exp");
    e12.DebugPrint("SU3_12_Exp");

    // T53: ExpReal
    printf("\n--- T53: ExpReal ---\n");
    deviceSU3 er3 = su3a.ExpReal(F(0.5), 10);
    deviceSU3 er12 = a12.ExpReal(F(0.5), 10);
    er3.DebugPrint("SU3_ExpReal");
    er12.DebugPrint("SU3_12_ExpReal");

    // T54: StrictExp
    printf("\n--- T54: StrictExp ---\n");
    deviceSU3 se3 = su3a.StrictExp();
    deviceSU3 se12 = a12.StrictExp();
    se3.DebugPrint("SU3_StrictExp");
    se12.DebugPrint("SU3_12_StrictExp");

    // T55: Log
    printf("\n--- T55: Log ---\n");
    deviceSU3 l3 = su3a.Log();
    deviceSU3 l12 = a12.Log();
    l3.DebugPrint("SU3_Log");
    l12.DebugPrint("SU3_12_Log");

    // T56: Power
    printf("\n--- T56: Power ---\n");
    deviceSU3 pw3 = su3a.Power(F(0.5));
    deviceSU3 pw12 = a12.Power(F(0.5));
    pw3.DebugPrint("SU3_Power");
    pw12.DebugPrint("SU3_12_Power");

    // ====================================================================
    // SECTION 7: EigenValues/Vectors
    // ====================================================================
    printf("\n=== SECTION 7: EigenValues/Vectors ===\n");

    // T57: CalculateEigenValues
    printf("\n--- T57: CalculateEigenValues ---\n");
    CLGComplex ev1_3, ev2_3, ev3_3;
    su3a.CalculateEigenValues(ev1_3, ev2_3, ev3_3);
    CLGComplex ev1_12, ev2_12, ev3_12;
    a12.CalculateEigenValues(ev1_12, ev2_12, ev3_12);
    printf("SU3 eigs: (%f,%f) (%f,%f) (%f,%f)\n", ev1_3.x, ev1_3.y, ev2_3.x, ev2_3.y, ev3_3.x, ev3_3.y);
    printf("SU3_12 eigs: (%f,%f) (%f,%f) (%f,%f)\n", ev1_12.x, ev1_12.y, ev2_12.x, ev2_12.y, ev3_12.x, ev3_12.y);

    // T58: EigenVectors
    printf("\n--- T58: EigenVectors ---\n");
    deviceSU3 evv3 = su3a.EigenVectors(ev1_3, ev2_3, ev3_3);
    deviceSU3 evv12 = a12.EigenVectors(ev1_12, ev2_12, ev3_12);
    evv3.DebugPrint("SU3_EigenVec");
    evv12.DebugPrint("SU3_12_EigenVec");

    // ====================================================================
    // SECTION 8: Static utility
    // ====================================================================
    printf("\n=== SECTION 8: Static utility ===\n");

    // T59: Determinent (static, takes CLGComplex*)
    printf("\n--- T59: Determinent (static) ---\n");
    CLGComplex det_static = deviceSU3_12::Determinent(su3a.m_me);
    printf("Determinent(static) = (%f,%f)\n", det_static.x, det_static.y);

    // T60: deviceUVW / deviceSqrtf012
    printf("\n--- T60: deviceUVW / deviceSqrtf012 ---\n");
    deviceSU3_12 q12 = a12.MulC(a12);
    deviceSU3_12 q212 = q12.MulC(q12);
    DOUBLE u3, v3, w3;
    deviceSU3 q3 = su3a; q3.Mul(su3a);
    deviceSU3 q23 = q3; q23.Mul(q3);
    deviceSU3::deviceUVW(q3, q23, u3, v3, w3);
    DOUBLE u12, v12, w12;
    deviceSU3_12::deviceUVW(q12, q212, u12, v12, w12);
    printf("SU3 UVW: %f %f %f\n", u3, v3, w3);
    printf("SU3_12 UVW: %f %f %f\n", u12, v12, w12);

    DOUBLE f0_3, f1_3, f2_3;
    deviceSU3::deviceSqrtf012(f0_3, f1_3, f2_3, u3, v3, w3);
    DOUBLE f0_12, f1_12, f2_12;
    deviceSU3_12::deviceSqrtf012(f0_12, f1_12, f2_12, u12, v12, w12);
    printf("SU3 sqrt012: %f %f %f\n", f0_3, f1_3, f2_3);
    printf("SU3_12 sqrt012: %f %f %f\n", f0_12, f1_12, f2_12);

    printf("\n===== deviceSU3_12 Tests Done =====\n");
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelDebugFunctionForSites()
{
    //intokernal;
    //printf("%d\n", uiSiteIndex);
    //for (UINT i = 0; i < EGM_MAX; ++i)
    //{
    //    printf("gammamatrix is %d\n", i);
    //    __chiralGamma[i].Print();
    //}
}

__global__ void _CLG_LAUNCH_BOUND
_kernelThreadBufferZeroReal(DOUBLE* arr, DOUBLE initial)
{
    intokernal;
    arr[uiSiteIndex] = initial;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelThreadBufferZeroComplex(cuDoubleComplex* arr, cuDoubleComplex initial)
{
    intokernal;
    arr[uiSiteIndex] = initial;
}

#if _CLG_DTK

__global__ void
_CLG_LAUNCH_BOUND
_kernelReduceRealOld(DOUBLE* arr, UINT uiJump, UINT uiMax)
{
    //for length 16 array
    //for jump = 1, this is 1->0, 3->2, 5->4, 7->6, 9->10, 11->10, 13->12, 15->14 
    //for jump = 2, this is 2->0, 6->4, 10->8, 14->12 
    //for jump = 4, this is 4->0, 12->8 
    //for jump = 8, this is 8->0, and is finished.

    //id target = idx * (jump << 1)
    //id from = target + jump
    UINT uiIdFrom = (threadIdx.x + blockIdx.x * blockDim.x) * (uiJump << 1) + uiJump;
    if (uiIdFrom < uiMax)
    {
        arr[uiIdFrom - uiJump] += arr[uiIdFrom];
    }
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelReduceCompOld(cuDoubleComplex* arr, UINT uiJump, UINT uiMax)
{
    UINT uiIdFrom = (threadIdx.x + blockIdx.x * blockDim.x) * (uiJump << 1) + uiJump;
    if (uiIdFrom < uiMax)
    {
        arr[uiIdFrom - uiJump] = cuCadd(arr[uiIdFrom - uiJump], arr[uiIdFrom]);
    }
}

#else

__global__ void
_CLG_LAUNCH_BOUND
_kernelReduceReal(DOUBLE* arr, UINT uiMax)
{
    extern __shared__ DOUBLE dsharedData[];
    UINT tid = threadIdx.x;
    UINT globalIdx = blockIdx.x * blockDim.x + tid;
    dsharedData[tid] = (globalIdx < uiMax) ? arr[globalIdx] : 0.0;

    __syncthreads();

    // Perform reduction in shared memory
    for (UINT stride = (blockDim.x >> 1); stride > 0; stride >>= 1)
    {
        if (tid < stride)
        {
            dsharedData[tid] += dsharedData[tid + stride];
        }
        __syncthreads();
    }

    // Write the result of this block to global memory
    if (0 == tid)
    {
        arr[blockIdx.x] = dsharedData[0];
    }
}

__global__ void 
_CLG_LAUNCH_BOUND
_kernelReduceComp(cuDoubleComplex* arr, UINT uiMax)
{
    extern __shared__ cuDoubleComplex sharedData[];
    UINT tid = threadIdx.x;
    UINT globalIdx = blockIdx.x * blockDim.x + tid;
    sharedData[tid] = (globalIdx < uiMax) ? arr[globalIdx] : make_cuDoubleComplex(0.0, 0.0);

    __syncthreads();

    // Perform reduction in shared memory
    for (UINT stride = (blockDim.x >> 1); stride > 0; stride >>= 1) 
    {
        if (tid < stride) 
        {
            sharedData[tid] = cuCadd(sharedData[tid], sharedData[tid + stride]);
        }
        __syncthreads();
    }

    // Write the result of this block to global memory
    if (0 == tid) 
    {
        arr[blockIdx.x] = sharedData[0];
    }
}

#endif

#pragma endregion

extern CLGAPI const char* _CLG_cudaGetErrorName(cudaError_t error)
{
    return cudaGetErrorName(error);
}

CCudaHelper::~CCudaHelper()
{
    ReleaseTemeraryBuffers();
    appSafeFree(m_pFunctionStackSpace);
}

void CCudaHelper::DeviceQuery()
{
    appGeneral(" CUDA Device Query (Runtime API) version (CUDART static linking)\n\n");

    INT deviceCount = 0;
    cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

    if (error_id != cudaSuccess) 
    {
        appGeneral("cudaGetDeviceCount returned %d\n-> %s\n",
            static_cast<INT>(error_id), cudaGetErrorString(error_id));
        appCrucial("Result = FAIL\n");
        _FAIL_EXIT;
    }

    // This function call returns 0 if there are no CUDA capable devices.
    if (deviceCount == 0) 
    {
        appGeneral("There are no available device(s) that support CUDA\n");
    }
    else 
    {
        appGeneral("Detected %d CUDA Capable device(s)\n", deviceCount);
    }

    INT dev, driverVersion = 0, runtimeVersion = 0;

    for (dev = 0; dev < deviceCount; ++dev) 
    {
        cudaSetDevice(dev);
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);

        appGeneral("\nDevice %d: \"%s\"\n", dev, deviceProp.name);

        // Console log
        cudaDriverGetVersion(&driverVersion);
        cudaRuntimeGetVersion(&runtimeVersion);
        appGeneral("  CUDA Driver Version / Runtime Version          %d.%d / %d.%d\n",
            driverVersion / 1000, (driverVersion % 100) / 10,
            runtimeVersion / 1000, (runtimeVersion % 100) / 10);
        appGeneral("  CUDA Capability Major/Minor version number:    %d.%d\n",
            deviceProp.major, deviceProp.minor);

        char msg[256];
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
        appSprintf(msg, sizeof(msg),
            "  Total amount of global memory:                 %.0f MBytes "
            "(%llu bytes)\n",
            static_cast<float>(deviceProp.totalGlobalMem / 1048576.0f),
            (unsigned long long)deviceProp.totalGlobalMem);
#else
        appSprintf(msg, sizeof(msg),
            "  Total amount of global memory:                 %.0f MBytes "
            "(%llu bytes)\n",
            static_cast<float>(deviceProp.totalGlobalMem / 1048576.0f),
            (unsigned long long)deviceProp.totalGlobalMem);
#endif
        appGeneral("%s", msg);

        appGeneral("  (%2d) Multiprocessors, (%3d) CUDA Cores/MP:     %d CUDA Cores\n",
            deviceProp.multiProcessorCount,
            _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor),
            _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) *
            deviceProp.multiProcessorCount);
        //appGeneral(
        //    "  GPU Max Clock rate:                            %.0f MHz (%0.2f "
        //    "GHz)\n",
        //    deviceProp.clockRate * 1e-3f, deviceProp.clockRate * 1e-6f);

#if CUDART_VERSION >= 5000
        // This is supported in CUDA 5.0 (runtime API device properties)
        //appGeneral("  Memory Clock rate:                             %.0f Mhz\n",
        //    deviceProp.memoryClockRate * 1e-3f);
        appGeneral("  Memory Bus Width:                              %d-bit\n",
            deviceProp.memoryBusWidth);

        if (deviceProp.l2CacheSize) 
        {
            appGeneral("  L2 Cache Size:                                 %d bytes\n",
                deviceProp.l2CacheSize);
        }

#else
        // This only available in CUDA 4.0-4.2 (but these were only exposed in the
        // CUDA Driver API)
        int memoryClock;
        getCudaAttribute<int>(&memoryClock, CU_DEVICE_ATTRIBUTE_MEMORY_CLOCK_RATE,
            dev);
        appGeneral("  Memory Clock rate:                             %.0f Mhz\n",
            memoryClock * 1e-3f);
        int memBusWidth;
        getCudaAttribute<int>(&memBusWidth,
            CU_DEVICE_ATTRIBUTE_GLOBAL_MEMORY_BUS_WIDTH, dev);
        appGeneral("  Memory Bus Width:                              %d-bit\n",
            memBusWidth);
        int L2CacheSize;
        getCudaAttribute<int>(&L2CacheSize, CU_DEVICE_ATTRIBUTE_L2_CACHE_SIZE, dev);

        if (L2CacheSize) 
        {
            appGeneral("  L2 Cache Size:                                 %d bytes\n",
                L2CacheSize);
        }

#endif

        appGeneral(
            "  Maximum Texture Dimension Size (x,y,z)         1D=(%d), 2D=(%d, "
            "%d), 3D=(%d, %d, %d)\n",
            deviceProp.maxTexture1D, deviceProp.maxTexture2D[0],
            deviceProp.maxTexture2D[1], deviceProp.maxTexture3D[0],
            deviceProp.maxTexture3D[1], deviceProp.maxTexture3D[2]);
        appGeneral(
            "  Maximum Layered 1D Texture Size, (num) layers  1D=(%d), %d layers\n",
            deviceProp.maxTexture1DLayered[0], deviceProp.maxTexture1DLayered[1]);
        appGeneral(
            "  Maximum Layered 2D Texture Size, (num) layers  2D=(%d, %d), %d "
            "layers\n",
            deviceProp.maxTexture2DLayered[0], deviceProp.maxTexture2DLayered[1],
            deviceProp.maxTexture2DLayered[2]);

        appGeneral("  Total amount of constant memory:               %lu bytes\n",
            deviceProp.totalConstMem);
        appGeneral("  Total amount of shared memory per block:       %lu bytes\n",
            deviceProp.sharedMemPerBlock);
        appGeneral("  Total number of registers available per block: %d\n",
            deviceProp.regsPerBlock);
        appGeneral("  Warp size:                                     %d\n",
            deviceProp.warpSize);
        appGeneral("  Maximum number of threads per multiprocessor:  %d\n",
            deviceProp.maxThreadsPerMultiProcessor);
        appGeneral("  Maximum number of threads per block:           %d\n",
            deviceProp.maxThreadsPerBlock);
        appGeneral("  Max dimension size of a thread block (x,y,z): (%d, %d, %d)\n",
            deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1],
            deviceProp.maxThreadsDim[2]);
        appGeneral("  Max dimension size of a grid size    (x,y,z): (%d, %d, %d)\n",
            deviceProp.maxGridSize[0], deviceProp.maxGridSize[1],
            deviceProp.maxGridSize[2]);
        appGeneral("  Maximum memory pitch:                          %lu bytes\n",
            deviceProp.memPitch);
        appGeneral("  Texture alignment:                             %lu bytes\n",
            deviceProp.textureAlignment);
#if CUDART_VERSION >= 13000
        //appGeneral(
        //    "  Concurrent copy and kernel execution:          %s with %d copy "
        //    "engine(s)\n",
        //    (deviceProp.deviceOverlap ? "Yes" : "No"), deviceProp.asyncEngineCount);
        //appGeneral("  Run time limit on kernels:                     %s\n",
        //    deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No");
#else
        appGeneral(
            "  Concurrent copy and kernel execution:          %s with %d copy "
            "engine(s)\n",
            (deviceProp.deviceOverlap ? "Yes" : "No"), deviceProp.asyncEngineCount);
        appGeneral("  Run time limit on kernels:                     %s\n",
            deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No");
#endif
        appGeneral("  Integrated GPU sharing Host Memory:            %s\n",
            deviceProp.integrated ? "Yes" : "No");
        appGeneral("  Support host page-locked memory mapping:       %s\n",
            deviceProp.canMapHostMemory ? "Yes" : "No");
        appGeneral("  Alignment requirement for Surfaces:            %s\n",
            deviceProp.surfaceAlignment ? "Yes" : "No");
        appGeneral("  Device has ECC support:                        %s\n",
            deviceProp.ECCEnabled ? "Enabled" : "Disabled");
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
        appGeneral("  CUDA Device Driver Mode (TCC or WDDM):         %s\n",
            deviceProp.tccDriver ? "TCC (Tesla Compute Cluster Driver)"
            : "WDDM (Windows Display Driver Model)");
#endif
        appGeneral("  Device supports Unified Addressing (UVA):      %s\n",
            deviceProp.unifiedAddressing ? "Yes" : "No");
        appGeneral("  Device supports Compute Preemption:            %s\n",
            deviceProp.computePreemptionSupported ? "Yes" : "No");
        appGeneral("  Supports Cooperative Kernel Launch:            %s\n",
            deviceProp.cooperativeLaunch ? "Yes" : "No");
#if CUDART_VERSION >= 13000
        //appGeneral("  Supports MultiDevice Co-op Kernel Launch:      %s\n",
        //    deviceProp.cooperativeMultiDeviceLaunch ? "Yes" : "No");
#else
        appGeneral("  Supports MultiDevice Co-op Kernel Launch:      %s\n",
            deviceProp.cooperativeMultiDeviceLaunch ? "Yes" : "No");
#endif
        appGeneral("  Device PCI Domain ID / Bus ID / location ID:   %d / %d / %d\n",
            deviceProp.pciDomainID, deviceProp.pciBusID, deviceProp.pciDeviceID);

        const char *sComputeMode[] = {
            "Default (multiple host threads can use ::cudaSetDevice() with device "
            "simultaneously)",
            "Exclusive (only one host thread in one process is able to use "
            "::cudaSetDevice() with this device)",
            "Prohibited (no host thread can use ::cudaSetDevice() with this "
            "device)",
            "Exclusive Process (many threads in one process is able to use "
            "::cudaSetDevice() with this device)",
            "Unknown",
            NULL };
        appGeneral("  Compute Mode:\n");
#if CUDART_VERSION >= 13000
        int computeMode = 0;
        cudaDeviceGetAttribute(&computeMode, cudaDevAttrComputeMode, dev);
        if (computeMode != cudaComputeModeProhibited) 
        {
            appGeneral("     < %s >\n", sComputeMode[computeMode]);
        }
#else
        appGeneral("     < %s >\n", sComputeMode[deviceProp.computeMode]);
#endif
    }

    // If there are 2 or more GPUs, query to determine whether RDMA is supported
    if (deviceCount >= 2) {
        cudaDeviceProp prop[64];
        int gpuid[64];  // we want to find the first two GPUs that can support P2P
        int gpu_p2p_count = 0;

        for (int i = 0; i < deviceCount; i++) 
        {
            checkCudaErrors(cudaGetDeviceProperties(&prop[i], i));

            // Only boards based on Fermi or later can support P2P
            if ((prop[i].major >= 2)
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
                // on Windows (64-bit), the Tesla Compute Cluster driver for windows
                // must be enabled to support this
                && prop[i].tccDriver
#endif
                ) {
                // This is an array of P2P capable GPUs
                gpuid[gpu_p2p_count++] = i;
            }
        }

        // Show all the combinations of support P2P GPUs
        int can_access_peer;

        if (gpu_p2p_count >= 2) {
            for (int i = 0; i < gpu_p2p_count; i++) 
            {
                for (int j = 0; j < gpu_p2p_count; j++) 
                {
                    if (gpuid[i] == gpuid[j]) 
                    {
                        continue;
                    }
                    checkCudaErrors(
                        cudaDeviceCanAccessPeer(&can_access_peer, gpuid[i], gpuid[j]));
                    appGeneral("> Peer access from %s (GPU%d) -> %s (GPU%d) : %s\n",
                        prop[gpuid[i]].name, gpuid[i], prop[gpuid[j]].name, gpuid[j],
                        can_access_peer ? "Yes" : "No");
                }
            }
        }
    }

    // csv masterlog info
    // *****************************
    // exe and CUDA driver name
    appGeneral("\n");
    std::string sProfileString = "deviceQuery, CUDA Driver = CUDART";
    char cTemp[16];

    // driver version
    sProfileString += ", CUDA Driver Version = ";
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    appSprintf(cTemp, 10, "%d.%d", driverVersion / 1000, (driverVersion % 100) / 10);
#else
    appSprintf(cTemp, sizeof(cTemp), "%d.%d", driverVersion / 1000,
        (driverVersion % 100) / 10);
#endif
    sProfileString += cTemp;

    // Runtime version
    sProfileString += ", CUDA Runtime Version = ";
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    appSprintf(cTemp, 10, "%d.%d", runtimeVersion / 1000, (runtimeVersion % 100) / 10);
#else
    appSprintf(cTemp, sizeof(cTemp), "%d.%d", runtimeVersion / 1000,
        (runtimeVersion % 100) / 10);
#endif
    sProfileString += cTemp;

    // Device count
    sProfileString += ", NumDevs = ";
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    appSprintf(cTemp, 10, "%d", deviceCount);
#else
    appSprintf(cTemp, sizeof(cTemp), "%d", deviceCount);
#endif
    sProfileString += cTemp;
    sProfileString += "\n";
    appGeneral("%s", sProfileString.c_str());

    appGeneral("Result = PASS\n");

}

void CCudaHelper::MemoryQuery()
{
    size_t availableMemory, totalMemory, usedMemory;
    cudaMemGetInfo(&availableMemory, &totalMemory);
    usedMemory = totalMemory - availableMemory;
    appGeneral(_T("Device Memory: used %llu, available %llu, total %llu\n"), usedMemory, availableMemory, totalMemory);
}

void CCudaHelper::DebugFunction()
{
    //appGeneral(_T("mult: %d, %d, %d\n"), _HC_MultX, _HC_MultY, _HC_MultZ);
    //appGeneral(_T("l: %d, %d, %d, %d\n"), _HC_Lx, _HC_Ly, _HC_Lz, _HC_Lt);

    _LAUNCH_KERNEL0(_kernelDebugFunction, 1, 1);
    checkCudaErrors(cudaDeviceSynchronize());

    preparethread;
    _LAUNCH_KERNEL0(_kernelDebugFunctionForSites, block, threads);
    checkCudaErrors(cudaDeviceSynchronize());

    //gammaMatrix testGamma[static_cast<INT>(EGM_MAX)];
    ////default is from device to host
    //checkCudaErrors(cudaMemcpyFromSymbol(testGamma, __chiralGamma, sizeof(gammaMatrix) * static_cast<INT>(EGM_MAX)));

    //for (INT i = 0; i < static_cast<INT>(EGM_MAX); ++i)
    //{
    //    appGeneral(_T("%s=\n"), __ENUM_TO_STRING(EGammaMatrix, (EGammaMatrix)i));
    //    testGamma[i].Print();
    //}

    //CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(1));
    //CFieldFermionHISQSU3* pFermion1 = dynamic_cast<CFieldFermionHISQSU3*>(appGetLattice()->GetFieldById(2));
    //CGaugeSmearingHISQSU3* pHISQSmearing = dynamic_cast<CGaugeSmearingHISQSU3*>(appGetLattice()->m_pGaugeSmearing[1]);

    //pHISQSmearing->GaugeSmearingC(pGauge);
    //TArray<const CFieldGauge*> gauges;
    //gauges.AddItem(pHISQSmearing->GetEffectiveGauge());

    //SFermionBosonSource source;
    //source.m_byColorIndex = 0;
    //source.m_eSourceType = EFS_Wall;
    //source.m_sSourcePoint = SSmallInt4(0, 0, 0, 0);

    //pFermion1->InitialAsSource(source);
    //CFieldFermionHISQSU3* pFermion2 = dynamic_cast<CFieldFermionHISQSU3*>(appGetLattice()->GetPooledCopy(pFermion1, __FILE__, __LINE__));
    ////pFermion->SaveToFile(_T("clg_wall_source_double.con"));
    //pFermion1->D0OnEvenOrOdd(TRUE, 1, 0, gauges.GetData(), NULL);
    //pFermion1->ZeroOnEvenOdd(TRUE);
    //pFermion2->D0OnEvenOrOdd(FALSE, 1, 0, gauges.GetData(), NULL);
    //pFermion2->ZeroOnEvenOdd(FALSE);
    //pFermion1->AxpyPlus(pFermion2);
    //pFermion1->SaveToFile(_T("clg_d_wall_2_cfg_100_r_double.con"));



    //pFermion2->Return();
    //CFieldGaugeU1Real testcomplete;

    //pFermion1->InverseD(1, 0, gauges.GetData(), NULL);
    //pFermion1->SaveToFile(_T("clg_wallsource_cfg_100_r_double.con"));

    //source.m_byColorIndex = 1;
    //pFermion1->InitialAsSource(source);
    //pFermion1->InverseD(1, 0, gauges.GetData(), NULL);
    //pFermion1->SaveToFile(_T("clg_wallsource_cfg_100_g_double.con"));

    //source.m_byColorIndex = 2;
    //pFermion1->InitialAsSource(source);
    //pFermion1->InverseD(1, 0, gauges.GetData(), NULL);
    //pFermion1->SaveToFile(_T("clg_wallsource_cfg_100_b_double.con"));

    //TArray<const CFieldGauge*> gaugeList;
    //gaugeList.AddItem(pGauge);
    //CAction* pAction = appGetLattice()->GetActionById(1);
    //appGeneral(_T("Action energy: %f\n"), pAction->Energy(FALSE, 1, 0, gaugeList.GetData(), NULL, NULL));

    //pGauge->DebugPrintMe();

    //CGaugeSmearing* pSmearing = appGetLattice()->m_pGaugeSmearing[1];
    //pSmearing->GaugeSmearingC(pGauge);

    //CFieldGaugeSU3* pGaugeCopy = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetPooledFieldById(1, __FILE__, __LINE__));

    //pSmearing->GetEffectiveGaugeLevel1()->CopyTo(pGaugeCopy);
    //pGaugeCopy->ApplyStaggeredPhase();
    //pGaugeCopy->SaveToFile(_T("cfg_100_clg_level1.con"));

    //pSmearing->GetEffectiveGauge()->CopyTo(pGaugeCopy);
    //pGaugeCopy->ApplyStaggeredPhase();
    //pGaugeCopy->SaveToFile(_T("cfg_100_clg_level2.con"));

    //pSmearing->GetNaikLink()->CopyTo(pGaugeCopy);
    //pGaugeCopy->ApplyStaggeredPhase();
    //pGaugeCopy->SaveToFile(_T("cfg_100_clg_naik.con"));
    //pGaugeCopy->Return();

    //CMeasurePlaqutteEnergy* pMeasure = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    //TArray<CFieldGauge*> gauges;
    //gauges.AddItem(pGauge);
    //pMeasure->OnConfigurationAccepted(1, 0, gauges.GetData(), NULL, NULL);
    
}

void CCudaHelper::CopyConstants() const
{
    checkCudaErrors(cudaMemcpyToSymbol(_constIntegers, m_ConstIntegers, sizeof(UINT) * kContentLength));
    checkCudaErrors(cudaMemcpyToSymbol(_constSignedIntegers, m_ConstSignedIntegers, sizeof(INT) * kContentLength));
    checkCudaErrors(cudaMemcpyToSymbol(_constFloats, m_ConstFloats, sizeof(Real) * kContentLength));
}

void CCudaHelper::CopyRandomPointer(const CRandom* r) const
{
    checkCudaErrors(cudaMemcpyToSymbol(__r, &r, sizeof(CRandom*)));
}

void CCudaHelper::CreateGammaMatrix() const
{
    gammaMatrix* pChiralGamma;
    deviceSU3* pSU3;

    //create pointer
    checkCudaErrors(__cudaMalloc((void**)&pChiralGamma, sizeof(gammaMatrix) * EGM_MAX));
    checkCudaErrors(__cudaMalloc((void**)&pSU3, sizeof(deviceSU3) * 9));

    //craete content
    _LAUNCH_KERNEL(_kernelCreateMatrix, 1, 1, pChiralGamma, pSU3);

    //copy to constant
    checkCudaErrors(cudaMemcpyToSymbol(__chiralGamma, pChiralGamma, sizeof(gammaMatrix) * EGM_MAX));
    checkCudaErrors(cudaMemcpyToSymbol(__SU3Generators, pSU3, sizeof(deviceSU3) * 9));

    //free pointers (already copy to constant, no need)
    checkCudaErrors(__cudaFree(pChiralGamma));
    checkCudaErrors(__cudaFree(pSU3));
}

void CCudaHelper::SetDeviceIndex(class CIndexData* pIdx) const
{
    checkCudaErrors(__cudaMalloc((void**)&m_pDevicePtrIndexData, sizeof(CIndexData)));
    checkCudaErrors(cudaMemcpy(m_pDevicePtrIndexData, pIdx, sizeof(CIndexData), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpyToSymbol(__idx, &m_pDevicePtrIndexData, sizeof(CIndexData*)));
}

void CCudaHelper::SetFieldPointers()
{
    for (BYTE i = 0; i < kMaxFieldCount; ++i)
    {
        CField * pField = appGetLattice()->GetFieldById(i);
        if (NULL == pField)
        {
            m_deviceFieldPointers[i] = NULL;
        }
        else
        {
            const UINT uiSize = pField->GetClass()->GetSize();
            CField* pDeviceField = NULL;
            checkCudaErrors(__cudaMalloc((void**)&pDeviceField, static_cast<size_t>(uiSize)));
            checkCudaErrors(cudaMemcpy(pDeviceField, pField, static_cast<size_t>(uiSize), cudaMemcpyHostToDevice));
            m_deviceFieldPointers[i] = pDeviceField;
        }

        CFieldBoundaryParent * pBoundaryField = appGetLattice()->GetBoundaryFieldById(i);
        if (NULL == pBoundaryField)
        {
            m_deviceBoundaryFieldPointers[i] = NULL;
        }
        else
        {
            const UINT uiSize = pBoundaryField->GetClass()->GetSize();
            CFieldBoundaryParent* pDeviceBoundaryField = NULL;
            checkCudaErrors(__cudaMalloc((void**)&pDeviceBoundaryField, static_cast<size_t>(uiSize)));
            checkCudaErrors(cudaMemcpy(pDeviceBoundaryField, pBoundaryField, static_cast<size_t>(uiSize), cudaMemcpyHostToDevice));
            m_deviceBoundaryFieldPointers[i] = pDeviceBoundaryField;
        }
    }

    checkCudaErrors(cudaMemcpyToSymbol(__fieldPointers, m_deviceFieldPointers, sizeof(CField*) * kMaxFieldCount));
    checkCudaErrors(cudaMemcpyToSymbol(__boundaryFieldPointers, m_deviceBoundaryFieldPointers, sizeof(CFieldBoundaryParent*) * kMaxFieldCount));
}

TArray<UINT> CCudaHelper::GetMaxThreadCountAndThreadPerblock(INT deviceId)
{
    TArray<UINT> ret;

    //INT deviceCount = 0;
    //const cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

    //if (error_id != cudaSuccess)
    //{
    //    appCrucial("cudaGetDeviceCount returned %d\n-> %s\n",
    //        static_cast<INT>(error_id), cudaGetErrorString(error_id));
    //    appCrucial("Result = FAIL\n");
    //    _FAIL_EXIT;
    //}

    //if (0 == deviceCount)
    //{
    //    appCrucial(_T("This program need GPU but you do NOT have a GPU.\n"));
    //    _FAIL_EXIT;
    //}

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, deviceId);

    //We need to constrain it further for shared memeory per block
    ret.AddItem(deviceProp.maxThreadsPerBlock);
    if (CCommonData::m_uiMaxThreadPerBlock > 0 && ret[0] > CCommonData::m_uiMaxThreadPerBlock)
    {
        ret[0] = CCommonData::m_uiMaxThreadPerBlock;
    }
    ret.AddItem(deviceProp.maxThreadsDim[0]);
    ret.AddItem(deviceProp.maxThreadsDim[1]);
    ret.AddItem(deviceProp.maxThreadsDim[2]);

    return ret;
}

/**
* The buffer size is HC_Volume
*/
void CCudaHelper::AllocateTemeraryBuffers(UINT uiThreadCount)
{
    m_uiThreadCount = uiThreadCount;
    m_uiReducePower = GetReduceDim((uiThreadCount + 1) >> 1);
    checkCudaErrors(__cudaMalloc((void**)&m_pRealBufferThreadCount, sizeof(DOUBLE) * uiThreadCount));
    checkCudaErrors(__cudaMalloc((void**)&m_pComplexBufferThreadCount, sizeof(cuDoubleComplex) * uiThreadCount));
}

cuDoubleComplex CCudaHelper::ThreadBufferSum(cuDoubleComplex* pDeviceBuffer)
{
    return ReduceComplexWithThreadCount(pDeviceBuffer);
}

DOUBLE CCudaHelper::ThreadBufferSum(DOUBLE* pDeviceBuffer)
{
    return ReduceRealWithThreadCount(pDeviceBuffer);
}

void CCudaHelper::ThreadBufferZero(cuDoubleComplex* pDeviceBuffer, cuDoubleComplex cInitial) const
{
    preparethread;
    _LAUN_KERNEL(_kernelThreadBufferZeroComplex, block, threads, pDeviceBuffer, cInitial);
    //_LAUNCH_KERNEL(_kernelThreadBufferZeroComplex, block, threads, pDeviceBuffer, cInitial);
}

void CCudaHelper::ThreadBufferZero(DOUBLE* pDeviceBuffer, DOUBLE fInitial) const
{
    preparethread;
    _LAUNCH_KERNEL(_kernelThreadBufferZeroReal, block, threads, pDeviceBuffer, fInitial);
}

#if _CLG_DTK

DOUBLE CCudaHelper::ReduceReal(DOUBLE* deviceBuffer, UINT uiLength)
{
    const UINT iRequiredDim = (uiLength + 1) >> 1;
    const UINT iPower = GetReduceDim(iRequiredDim);
    for (UINT i = 0; i <= iPower; ++i)
    {
        UINT iJump = 1 << i;
        UINT iThreadNeeded = 1 << (iPower - i);
        UINT iBlock = iThreadNeeded > _HC_ThreadConstraint ? iThreadNeeded / _HC_ThreadConstraint : 1;
        UINT iThread = iThreadNeeded > _HC_ThreadConstraint ? _HC_ThreadConstraint : iThreadNeeded;
        _LAUNCH_KERNEL(_kernelReduceRealOld, iBlock, iThread, deviceBuffer, iJump, uiLength);
    }
    DOUBLE result[1];
    cudaMemcpy(result, deviceBuffer, sizeof(DOUBLE), cudaMemcpyDeviceToHost);
    return result[0];
}

DOUBLE CCudaHelper::ReduceRealWithThreadCount(DOUBLE* deviceBuffer)
{
    for (UINT i = 0; i <= m_uiReducePower; ++i)
    {
        UINT iJump = 1 << i;
        UINT iThreadNeeded = 1 << (m_uiReducePower - i);
        UINT iBlock = iThreadNeeded > _HC_ThreadConstraint ? iThreadNeeded / _HC_ThreadConstraint : 1;
        UINT iThread = iThreadNeeded > _HC_ThreadConstraint ? _HC_ThreadConstraint : iThreadNeeded;
        _LAUNCH_KERNEL(_kernelReduceRealOld, iBlock, iThread, deviceBuffer, iJump, m_uiThreadCount);
    }
    DOUBLE result[1];
    cudaMemcpy(result, deviceBuffer, sizeof(DOUBLE), cudaMemcpyDeviceToHost);
    return result[0];
}

cuDoubleComplex CCudaHelper::ReduceComplex(cuDoubleComplex* deviceBuffer, UINT uiLength)
{
    const UINT iRequiredDim = (uiLength + 1) >> 1;
    const UINT iPower = GetReduceDim(iRequiredDim);
    for (UINT i = 0; i <= iPower; ++i)
    {
        UINT iJump = 1 << i;
        UINT iThreadNeeded = 1 << (iPower - i);
        UINT iBlock = iThreadNeeded > _HC_ThreadConstraint ? iThreadNeeded / _HC_ThreadConstraint : 1;
        UINT iThread = iThreadNeeded > _HC_ThreadConstraint ? _HC_ThreadConstraint : iThreadNeeded;
        _LAUNCH_KERNEL(_kernelReduceCompOld, iBlock, iThread, deviceBuffer, iJump, uiLength);
    }
    cuDoubleComplex result[1];
    cudaMemcpy(result, deviceBuffer, sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    return result[0];
}

cuDoubleComplex CCudaHelper::ReduceComplexWithThreadCount(cuDoubleComplex* deviceBuffer)
{
    for (UINT i = 0; i <= m_uiReducePower; ++i)
    {
        UINT iJump = 1 << i;
        UINT iThreadNeeded = 1 << (m_uiReducePower - i);
        UINT iBlock = iThreadNeeded > _HC_ThreadConstraint ? iThreadNeeded / _HC_ThreadConstraint : 1;
        UINT iThread = iThreadNeeded > _HC_ThreadConstraint ? _HC_ThreadConstraint : iThreadNeeded;
        _LAUNCH_KERNEL(_kernelReduceCompOld, iBlock, iThread, deviceBuffer, iJump, m_uiThreadCount);
    }
    cuDoubleComplex result[1];
    cudaMemcpy(result, deviceBuffer, sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    return result[0];
}

#else

DOUBLE CCudaHelper::ReduceReal(DOUBLE* deviceBuffer, UINT uiLength)
{
    UINT iThread = _HC_ThreadConstraint;
    UINT iBlock = (uiLength + iThread - 1) / iThread;

    // Allocate shared memory size
    size_t sharedMemSize = iThread * sizeof(DOUBLE);

    // Perform reduction in multiple steps
    while (uiLength > 1)
    {
        _LAUNCH_KERNELS(_kernelReduceReal, iBlock, iThread, sharedMemSize, deviceBuffer, uiLength);
        cudaDeviceSynchronize(); // Ensure the kernel completes before proceeding

        // Update the input buffer and dimensions for the next iteration
        uiLength = iBlock;
        iBlock = (uiLength + iThread - 1) / iThread;
    }

    DOUBLE result;
    cudaMemcpy(&result, deviceBuffer, sizeof(DOUBLE), cudaMemcpyDeviceToHost);
    return result;
}

DOUBLE CCudaHelper::ReduceRealWithThreadCount(DOUBLE* deviceBuffer)
{
    _RECORD(CCudaHelper::ReduceRealWithThreadCount);
    UINT uiLength = m_uiThreadCount;
    UINT iThread = _HC_ThreadConstraint;
    UINT iBlock = (uiLength + iThread - 1) / iThread;

    // Allocate shared memory size
    size_t sharedMemSize = iThread * sizeof(DOUBLE);

    // Perform reduction in multiple steps
    while (uiLength > 1)
    {
        _LAUNCH_KERNELS(_kernelReduceReal, iBlock, iThread, sharedMemSize, deviceBuffer, uiLength);
        cudaDeviceSynchronize(); // Ensure the kernel completes before proceeding

        // Update the input buffer and dimensions for the next iteration
        uiLength = iBlock;
        iBlock = (uiLength + iThread - 1) / iThread;
    }

    DOUBLE result;
    cudaMemcpy(&result, deviceBuffer, sizeof(DOUBLE), cudaMemcpyDeviceToHost);
    return result;
}

cuDoubleComplex CCudaHelper::ReduceComplex(cuDoubleComplex* deviceBuffer, UINT uiLength)
{
    UINT iThread = _HC_ThreadConstraint;
    UINT iBlock = (uiLength + iThread - 1) / iThread;

    // Allocate shared memory size
    size_t sharedMemSize = iThread * sizeof(cuDoubleComplex);

    // Perform reduction in multiple steps
    while (uiLength > 1)
    {
        _LAUNCH_KERNELS(_kernelReduceComp, iBlock, iThread, sharedMemSize, deviceBuffer, uiLength);
        cudaDeviceSynchronize(); // Ensure the kernel completes before proceeding

        // Update the input buffer and dimensions for the next iteration
        uiLength = iBlock;
        iBlock = (uiLength + iThread - 1) / iThread;
    }

    cuDoubleComplex result;
    cudaMemcpy(&result, deviceBuffer, sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    return result;
}

cuDoubleComplex CCudaHelper::ReduceComplexWithThreadCount(cuDoubleComplex* deviceBuffer)
{
    _RECORD(CCudaHelper::ReduceComplexWithThreadCount);
    UINT uiLength = m_uiThreadCount;
    UINT iThread = _HC_ThreadConstraint;
    UINT iBlock = (uiLength + iThread - 1) / iThread;

    // Allocate shared memory size
    size_t sharedMemSize = iThread * sizeof(cuDoubleComplex);

    // Perform reduction in multiple steps
    while (uiLength > 1)
    {
        _LAUNCH_KERNELS(_kernelReduceComp, iBlock, iThread, sharedMemSize, deviceBuffer, uiLength);
        cudaDeviceSynchronize(); // Ensure the kernel completes before proceeding

        // Update the input buffer and dimensions for the next iteration
        uiLength = iBlock;
        iBlock = (uiLength + iThread - 1) / iThread;
    }

    cuDoubleComplex result;
    cudaMemcpy(&result, deviceBuffer, sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    return result;
}

#endif

void CLGAPI appExistCuda()
{
    cudaDeviceReset();
}

void CLGAPI appSynchronize()
{
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
}

extern CLGAPI void launchKernel(const void* function, dim3 block, dim3 thread, void** args)
{
    checkCudaErrors(cudaLaunchKernel(function, block, thread, args));
    _CHECKCUDA;
}

extern CLGAPI void launchKernel(const void* function, UINT block, UINT thread, void** args)
{
    dim3 blockdim3(block, 1, 1);
    dim3 threaddim3(thread, 1, 1);
    checkCudaErrors(cudaLaunchKernel(function, blockdim3, threaddim3, args));
    _CHECKCUDA;
}

extern CLGAPI void launchKernel(const void* function, dim3 block, dim3 thread, size_t shareMem, void** args)
{
    checkCudaErrors(cudaLaunchKernel(function, block, thread, args, shareMem));
    _CHECKCUDA;
}

extern CLGAPI void launchKernel(const void* function, UINT block, UINT thread, size_t shareMem, void** args)
{
    dim3 blockdim3(block, 1, 1);
    dim3 threaddim3(thread, 1, 1);
    checkCudaErrors(cudaLaunchKernel(function, blockdim3, threaddim3, args, shareMem));
    _CHECKCUDA;
}

#pragma region Improve-1 single-evaluation probe

#if _CLG_MULTI_GPU

__global__ void _kernelImprove5Probe(UINT* pMarker, UINT uiA, UINT uiB)
{
    if (0 == threadIdx.x && 0 == blockIdx.x)
    {
        pMarker[0] = uiA + uiB;
    }
}

static UINT _improve5ProbeSideEffect(UINT* pCounts, UINT uiSlot, UINT uiRet)
{
    ++pCounts[uiSlot];
    return uiRet;
}

extern CLGAPI void appLaunchGuardSingleEvalProbe(UINT* puiCounts, UINT* puiDeviceMarker)
{
    //I5 gate: block, thread and every kernel-argument expression carries a
    //side effect and must be evaluated EXACTLY once by the macro (either
    //backend); the marker proves the kernel actually ran with 7 and 8.
    _LAUNCH_KERNEL(_kernelImprove5Probe,
        _improve5ProbeSideEffect(puiCounts, 0, 1),
        _improve5ProbeSideEffect(puiCounts, 1, 1),
        puiDeviceMarker,
        _improve5ProbeSideEffect(puiCounts, 2, 7),
        _improve5ProbeSideEffect(puiCounts, 3, 8));
    _CHECKCUDA;
}

#endif

#pragma endregion

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================