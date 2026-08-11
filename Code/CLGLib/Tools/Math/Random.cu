//=============================================================================
// FILENAME : Random.h
// 
// DESCRIPTION:
//
//
// REVISION:
//  [12/6/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

/**
* Multi-GPU (Phase 3, plan section 1.4-R2): RNG seeding by GLOBAL site index.
*
* curand states are stored at the LOCAL site/link index (that is this rank's
* storage), but the seed subsequence/offset must be derived from the GLOBAL
* coordinate so that the site living at global (gx,gy,gz,gt) draws the SAME
* stream regardless of how the lattice is decomposed. Without this, an n1 run
* and an nN run start HMC from different per-site momenta and their trajectories
* diverge, so the 1-vs-N energy comparison is meaningless.
*
* On single-GPU builds every offset is 0 and Global* == local L, so this is the
* identity map local-linear -> same-value global-linear, and behaviour is
* bit-identical to before.
*/
__device__ __inline__ static UINT _deviceGlobalSiteSeedIndex(UINT uiLocalSiteIndex)
{
    //NOTE: the RNG init kernels run BEFORE the device site-mapping table
    //(__idx->m_pSiteMappingTable) is populated, so we must NOT call
    //__deviceSiteIndexToInt4 here (it dereferences that table -> illegal address).
    //Constant memory (_constIntegers: _DC_Lx.., _DC_Offset.., _DC_GlobalL..) IS set
    //by this point (the kernels already read _DC_Seed), so decode the raw local
    //linear site index arithmetically. The RNG launch layout [Lx*Ly, Lz, Lt] yields
    //exactly the standard ordering ((x*Ly+y)*Lz+z)*Lt+w, matching this decode.
    const INT lx = _DC_Lxi;
    const INT ly = _DC_Lyi;
    const INT lz = _DC_Lzi;
    const INT lt = static_cast<INT>(_DC_Lt);
    const INT idx = static_cast<INT>(uiLocalSiteIndex);
    const INT w = idx % lt;
    const INT z = (idx / lt) % lz;
    const INT y = (idx / (lt * lz)) % ly;
    const INT x = idx / (lt * lz * ly);
    (void)lx;
    const INT gx = x + static_cast<INT>(_DC_OffsetX);
    const INT gy = y + static_cast<INT>(_DC_OffsetY);
    const INT gz = z + static_cast<INT>(_DC_OffsetZ);
    const INT gw = w + static_cast<INT>(_DC_OffsetT);
    return static_cast<UINT>(
        ((gx * static_cast<INT>(_DC_GlobalLy) + gy)
            * static_cast<INT>(_DC_GlobalLz) + gz)
            * static_cast<INT>(_DC_GlobalLt) + gw);
}

__global__ void _CLG_LAUNCH_BOUND
_kernalAllocateSeedTable(UINT* pDevicePtr)
{
    intokernaldir;

    const UINT uiSeed = _DC_Seed;
    const UINT uiGlobalSite = _deviceGlobalSiteSeedIndex(uiSiteIndex);

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //store at LOCAL link index, but derive the seed value from the GLOBAL
        //link index so the stream is decomposition-invariant (plan 1.4-R2).
        const UINT uiLocalFat = _deviceGetLinkIndex(uiSiteIndex, idir);
        const UINT uiGlobalFat = uiGlobalSite * uiDir + idir;
        pDevicePtr[uiLocalFat] = (1664525UL * (uiGlobalFat + uiSeed) + 1013904223UL) & 0xffffffffUL;
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernalInitialXORWOW(curandState * states)
{
    UINT uiSiteIndex = ((threadIdx.x + blockIdx.x * blockDim.x) 
        * blockDim.y * gridDim.y * blockDim.z * gridDim.z 
        + (threadIdx.y + blockIdx.y * blockDim.y) 
        * blockDim.z * gridDim.z 
        + (threadIdx.z + blockIdx.z * blockDim.z));
    const UINT uiSeedIndex = _deviceGlobalSiteSeedIndex(uiSiteIndex);
    curand_init(_DC_Seed, uiSeedIndex, uiSeedIndex, &states[uiSiteIndex]);
}

#if _CLG_USE_PHILOX4
__global__ void _CLG_LAUNCH_BOUND
_kernalInitialPhilox(curandStatePhilox4_32_10_t * states)
{
    const UINT uiSiteIndex = ((threadIdx.x + blockIdx.x * blockDim.x) * blockDim.y * gridDim.y * blockDim.z * gridDim.z + (threadIdx.y + blockIdx.y * blockDim.y) * blockDim.z * gridDim.z + (threadIdx.z + blockIdx.z * blockDim.z));
    const UINT uiSeed = _DC_Seed;
    const UINT uiDir = _DC_Dir;
    const UINT uiGlobalSite = _deviceGlobalSiteSeedIndex(uiSiteIndex);
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT fatIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        //store at local link index, seed by global link subsequence
        curand_init(uiSeed, uiGlobalSite * uiDir + idir, 0, &states[fatIndex]);
    }
}
#endif

#if _CLG_USE_MRG32K3A
__global__ void _CLG_LAUNCH_BOUND
_kernalInitialMRG(curandStateMRG32k3a  * states)
{
    const UINT uiSiteIndex = ((threadIdx.x + blockIdx.x * blockDim.x) * blockDim.y * gridDim.y * blockDim.z * gridDim.z + (threadIdx.y + blockIdx.y * blockDim.y) * blockDim.z * gridDim.z + (threadIdx.z + blockIdx.z * blockDim.z));
    const UINT uiSeed = _DC_Seed;
    const UINT uiDir = _DC_Dir;
    const UINT uiGlobalSite = _deviceGlobalSiteSeedIndex(uiSiteIndex);
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT fatIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        //store at local link index, seed by global link subsequence
        curand_init(uiSeed, uiGlobalSite * uiDir + idir, 0, &states[fatIndex]);
    }
}
#endif

#if _CLG_USE_QUASI_SOBOL32
__global__ void _CLG_LAUNCH_BOUND
_kernalInitialSobel32(curandStateSobol32* states, curandDirectionVectors32_t* dirs)
{
    intokernal;
    curand_init(dirs[uiSiteIndex], _DC_Seed % 16, &states[uiSiteIndex]);
}
#endif

#if _CLG_USE_SCRAMBLED_SOBOL32
__global__ void _CLG_LAUNCH_BOUND
_kernalInitialScrambledSobel32(curandStateScrambledSobol32* states, UINT* consts, curandDirectionVectors32_t* dirs)
{
    intokernal;
    curand_init(dirs[uiSiteIndex], consts[uiSiteIndex], _DC_Seed % __SOBEL_OFFSET_MAX, &states[uiSiteIndex]);
}
#endif

CRandom::~CRandom()
{

    switch (m_eRandomType)
    {
    case ER_Schrage:
        {
            checkCudaErrors(__cudaFree(m_pDeviceSeedTable));
        }
        break;
#if _CLG_USE_MRG32K3A
    case ER_MRG32K3A:
        {
            CURAND_CALL(curandDestroyGenerator(m_HGen));
            checkCudaErrors(__cudaFree(m_deviceBuffer));
            checkCudaErrors(__cudaFree(m_pDeviceRandStatesMRG));
        }
        break;
#endif
#if _CLG_USE_PHILOX4
    case ER_PHILOX4_32_10:
        {
            CURAND_CALL(curandDestroyGenerator(m_HGen));
            checkCudaErrors(__cudaFree(m_deviceBuffer));
            checkCudaErrors(__cudaFree(m_pDeviceRandStatesPhilox));
        }
        break;
#endif
#if _CLG_USE_QUASI_SOBOL32
    case ER_QUASI_SOBOL32:
        {
            CURAND_CALL(curandDestroyGenerator(m_HGen));
            checkCudaErrors(__cudaFree(m_deviceBuffer));
            checkCudaErrors(__cudaFree(m_pDeviceRandStatesSobol32));
            checkCudaErrors(__cudaFree(m_pDeviceSobolDirVec));
        }
        break;
#endif
#if _CLG_USE_SCRAMBLED_SOBOL32
    case ER_SCRAMBLED_SOBOL32:
        {
            CURAND_CALL(curandDestroyGenerator(m_HGen));
            checkCudaErrors(__cudaFree(m_deviceBuffer));
            checkCudaErrors(__cudaFree(m_pDeviceRandStatesScrambledSobol32));
            checkCudaErrors(__cudaFree(m_pDeviceSobolDirVec));
            checkCudaErrors(__cudaFree(m_pDeviceSobelConsts));
        }
        break;
#endif
    case ER_XORWOW:
        default:
        {
            CURAND_CALL(curandDestroyGenerator(m_HGen));
            checkCudaErrors(__cudaFree(m_deviceBuffer));
            checkCudaErrors(__cudaFree(m_pDeviceRandStatesXORWOW));
        }
        break;
    }
}

//Initial XORWOW only support 512 threads per block
void CRandom::InitialStatesXORWOW(UINT )
{
    //m_uiFatIdDivide = _HC_Dir + 1;
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceRandStatesXORWOW, sizeof(curandState) * _HC_Volume));
    TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock(appGetDeviceId());
    deviceConstraints[0] = min(512, _HC_ThreadCountPerBlock);
    TArray<UINT> latticeDim;
    latticeDim.AddItem(_HC_Lx * _HC_Ly);
    latticeDim.AddItem(_HC_Lz);
    latticeDim.AddItem(_HC_Lt);
    TArray <UINT> decomp = _getDecompose(deviceConstraints, latticeDim);
    dim3 block(decomp[0], decomp[1], decomp[2]);
    dim3 threads(decomp[3], decomp[4], decomp[5]);
    _LAUNCH_KERNEL(_kernalInitialXORWOW, block, threads, m_pDeviceRandStatesXORWOW);
}

#if _CLG_USE_PHILOX4
//Initial Philox only support 256 threads per block
void CRandom::InitialStatesPhilox(UINT )
{
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceRandStatesPhilox, sizeof(curandStatePhilox4_32_10_t) * _HC_Volume * _HC_Dir));

    TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock(appGetDeviceId());
    deviceConstraints[0] = min(256, _HC_ThreadCountPerBlock);
    TArray<UINT> latticeDim;
    latticeDim.AddItem(_HC_Lx * _HC_Ly);
    latticeDim.AddItem(_HC_Lz);
    latticeDim.AddItem(_HC_Lt);
    TArray <UINT> decomp = _getDecompose(deviceConstraints, latticeDim);
    dim3 block(decomp[0], decomp[1], decomp[2]);
    dim3 threads(decomp[3], decomp[4], decomp[5]);

    _LAUNCH_KERNEL(_kernalInitialPhilox, block, threads, m_pDeviceRandStatesPhilox);
}
#endif

#if _CLG_USE_MRG32K3A
//Initial MRG only support 256 threads per block
void CRandom::InitialStatesMRG(UINT )
{
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceRandStatesMRG, sizeof(curandStateMRG32k3a) * _HC_Volume * _HC_Dir));
    TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock(appGetDeviceId());
    deviceConstraints[0] = min(256, _HC_ThreadCountPerBlock);
    TArray<UINT> latticeDim;
    latticeDim.AddItem(_HC_Lx * _HC_Ly);
    latticeDim.AddItem(_HC_Lz);
    latticeDim.AddItem(_HC_Lt);
    TArray <UINT> decomp = _getDecompose(deviceConstraints, latticeDim);
    dim3 block(decomp[0], decomp[1], decomp[2]);
    dim3 threads(decomp[3], decomp[4], decomp[5]);
    _LAUNCH_KERNEL(_kernalInitialMRG, block, threads, m_pDeviceRandStatesMRG);
}
#endif

#if _CLG_USE_QUASI_SOBOL32
void CRandom::InitialStatesSobol32(UINT )
{
    //support only 20000 dimensions, so using _HC_Volumn instead
    //m_uiFatIdDivide = _HC_Dir + 1;
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceRandStatesSobol32,
        sizeof(curandStateSobol32) * _HC_Volume));
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceSobolDirVec,
        sizeof(curandDirectionVectors32_t) * _HC_Volume));

    //int[32]
    curandDirectionVectors32_t *hostVectors32;
    CURAND_CALL(curandGetDirectionVectors32(&hostVectors32, CURAND_DIRECTION_VECTORS_32_JOEKUO6));
    checkCudaErrors(cudaMemcpy(m_pDeviceSobolDirVec, hostVectors32, 
        _HC_Volume * sizeof(curandDirectionVectors32_t),
        cudaMemcpyHostToDevice));

    preparethread;
    _LAUNCH_KERNEL(_kernalInitialSobel32, block, threads, m_pDeviceRandStatesSobol32, m_pDeviceSobolDirVec);
}
#endif

#if _CLG_USE_SCRAMBLED_SOBOL32
void CRandom::InitialStatesScrambledSobol32(UINT )
{
    //m_uiFatIdDivide = _HC_Dir + 1;
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceRandStatesScrambledSobol32,
        sizeof(curandStateScrambledSobol32) * _HC_Volume));
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceSobolDirVec,
        sizeof(curandDirectionVectors32_t) * _HC_Volume));
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceSobelConsts,
        sizeof(UINT) * _HC_Volume));

    curandDirectionVectors32_t *hostVectors32;
    CURAND_CALL(curandGetDirectionVectors32(&hostVectors32, CURAND_SCRAMBLED_DIRECTION_VECTORS_32_JOEKUO6));
    checkCudaErrors(cudaMemcpy(
        m_pDeviceSobolDirVec, 
        hostVectors32, 
        _HC_Volume * sizeof(curandDirectionVectors32_t),
        cudaMemcpyHostToDevice));

    UINT * hostScrambleConstants32;
    CURAND_CALL(curandGetScrambleConstants32(&hostScrambleConstants32));
    checkCudaErrors(cudaMemcpy(
        m_pDeviceSobelConsts, 
        hostScrambleConstants32, 
        _HC_Volume * sizeof(UINT), 
        cudaMemcpyHostToDevice));

    preparethread;
    _LAUNCH_KERNEL(_kernalInitialScrambledSobel32, block, threads, m_pDeviceRandStatesScrambledSobol32, m_pDeviceSobelConsts, m_pDeviceSobolDirVec);
}
#endif

#if _CLG_MULTI_GPU
void CRandom::EnterGlobalContext()
{
    //P4-2.5 fix: rebuild the curand state arrays for the GLOBAL lattice. Save the
    //LOCAL (decomposed) arrays first; the rebuild allocates at the global volume
    //(_HC_Volume is global under the fixer context) and reseeds with
    //_deviceGlobalSiteSeedIndex, which degenerates to identity there (MGEnter
    //zeroes the offsets and sets GlobalL to the full lattice) -> the reseeded
    //states are globally consistent, identical to a single-GPU run.
    m_pSavedRandStatesXORWOW = m_pDeviceRandStatesXORWOW;
    InitialStatesXORWOW(_DC_Seed);
#if _CLG_USE_PHILOX4
    m_pSavedRandStatesPhilox = m_pDeviceRandStatesPhilox;
    InitialStatesPhilox(_DC_Seed);
#endif
#if _CLG_USE_MRG32K3A
    m_pSavedRandStatesMRG = m_pDeviceRandStatesMRG;
    InitialStatesMRG(_DC_Seed);
#endif
#if _CLG_USE_QUASI_SOBOL32
    m_pSavedRandStatesSobol32 = m_pDeviceRandStatesSobol32;
    m_pSavedSobolDirVec = m_pDeviceSobolDirVec;
    InitialStatesSobol32(_DC_Seed);
#endif
#if _CLG_USE_SCRAMBLED_SOBOL32
    m_pSavedRandStatesScrambledSobol32 = m_pDeviceRandStatesScrambledSobol32;
    m_pSavedSobolDirVec = m_pDeviceSobolDirVec;
    m_pSavedSobelConsts = m_pDeviceSobelConsts;
    InitialStatesScrambledSobol32(_DC_Seed);
#endif
    //P4-2.5 fix: __r (device constant) points at CLatticeData::m_pDeviceRandom,
    //the DEVICE COPY of this object. The rebuild above only updated the HOST
    //members (m_pDeviceRandStatesXORWOW & friends); kernels launched under the
    //temporary GLOBAL context read the copy via __r->_deviceRandomF and would
    //still see the LOCAL-volume arrays -> out-of-bounds on global site indices
    //(verified: cudaErrorIllegalAddress in _kernelRandomGauge on -n2). Re-sync
    //the whole object so the copy carries the rebuilt GLOBAL pointers.
    if (NULL != appGetLattice() && NULL != appGetLattice()->m_pDeviceRandom)
    {
        checkCudaErrors(cudaMemcpy(appGetLattice()->m_pDeviceRandom, this, sizeof(CRandom), cudaMemcpyHostToDevice));
    }
}

void CRandom::ExitGlobalContext()
{
    //Free the temporary GLOBAL arrays and restore the LOCAL ones. The destructor
    //frees whatever is current, so it must always see the local pointers after.
    checkCudaErrors(__cudaFree(m_pDeviceRandStatesXORWOW));
    m_pDeviceRandStatesXORWOW = m_pSavedRandStatesXORWOW;
    m_pSavedRandStatesXORWOW = NULL;
#if _CLG_USE_PHILOX4
    checkCudaErrors(__cudaFree(m_pDeviceRandStatesPhilox));
    m_pDeviceRandStatesPhilox = m_pSavedRandStatesPhilox;
    m_pSavedRandStatesPhilox = NULL;
#endif
#if _CLG_USE_MRG32K3A
    checkCudaErrors(__cudaFree(m_pDeviceRandStatesMRG));
    m_pDeviceRandStatesMRG = m_pSavedRandStatesMRG;
    m_pSavedRandStatesMRG = NULL;
#endif
#if _CLG_USE_QUASI_SOBOL32
    checkCudaErrors(__cudaFree(m_pDeviceRandStatesSobol32));
    m_pDeviceRandStatesSobol32 = m_pSavedRandStatesSobol32;
    m_pSavedRandStatesSobol32 = NULL;
#endif
#if _CLG_USE_SCRAMBLED_SOBOL32
    checkCudaErrors(__cudaFree(m_pDeviceRandStatesScrambledSobol32));
    m_pDeviceRandStatesScrambledSobol32 = m_pSavedRandStatesScrambledSobol32;
    m_pSavedRandStatesScrambledSobol32 = NULL;
    checkCudaErrors(__cudaFree(m_pDeviceSobelConsts));
    m_pDeviceSobelConsts = m_pSavedSobelConsts;
    m_pSavedSobelConsts = NULL;
#endif
#if _CLG_USE_QUASI_SOBOL32 || _CLG_USE_SCRAMBLED_SOBOL32
    //m_pDeviceSobolDirVec is shared by the Sobol variants and re-allocated by the
    //last InitialStates* call; free the rebuilt one and restore the original.
    checkCudaErrors(__cudaFree(m_pDeviceSobolDirVec));
    m_pDeviceSobolDirVec = m_pSavedSobolDirVec;
    m_pSavedSobolDirVec = NULL;
#endif
    //P4-2.5 fix: mirror Exit to the device copy (see EnterGlobalContext). The
    //destructor later frees whatever m_pDeviceRandStatesXORWOW & friends point
    //to on the HOST object, and the device copy must see the restored LOCAL
    //pointers again for any subsequent single-GPU (local-context) kernels.
    if (NULL != appGetLattice() && NULL != appGetLattice()->m_pDeviceRandom)
    {
        checkCudaErrors(cudaMemcpy(appGetLattice()->m_pDeviceRandom, this, sizeof(CRandom), cudaMemcpyHostToDevice));
    }
}
#endif

void CRandom::InitialTableSchrage(UINT )
{
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceSeedTable, sizeof(UINT) * _HC_Volume * _HC_Dir));
    preparethread;
    _LAUNCH_KERNEL(_kernalAllocateSeedTable, block, threads, m_pDeviceSeedTable);
}

UINT CRandom::DebugSeedTable() const
{
    UINT* hostseedtable = (UINT*)(malloc(sizeof(UINT) * _HC_Volume * _HC_Dir));
    checkCudaErrors(cudaMemcpy(hostseedtable, m_pDeviceSeedTable, sizeof(UINT) * _HC_Volume * _HC_Dir, cudaMemcpyDeviceToHost));

    appGeneral(_T("Cuda12.9 debug should be %d, %d, %d, %d\n"), 1187256132, -436264799, 846144630, 2128554059);
    appGeneral(_T("Vaules: %d, %d, %d, %d\n"), hostseedtable[0], hostseedtable[1], hostseedtable[2], hostseedtable[3]);

    
    UINT uiError = 0;
    if (static_cast<INT>(hostseedtable[0]) != 1566952946)
    {
        ++uiError;
    }
    if (static_cast<INT>(hostseedtable[1]) != -436264799)
    {
        ++uiError;
    }
    if (static_cast<INT>(hostseedtable[2]) != 846144630)
    {
        ++uiError;
    }
    if (static_cast<INT>(hostseedtable[3]) != 2128554059)
    {
        ++uiError;
    }

    appSafeFree(hostseedtable);
    return uiError;
}



Real GetRandomReal()
{
    return appGetLattice()->m_pRandom->GetRandomF();
}

#pragma region Test

__global__ void _CLG_LAUNCH_BOUND
_kernelMCPi(UINT* output, UINT lengthyz, UINT lengthz, UINT uiLoop, UINT uithreadCount)
{
    __shared__ UINT sData1[1024];
    __shared__ UINT sData2[1024];
    UINT uiToAdd = 0;
    UINT uiToAdd2 = 0;
    //We have a very large grid, but for a block, it is always smaller (or equval to volumn)
    const UINT fatIndex = threadIdx.x * lengthyz + threadIdx.y * lengthz + threadIdx.z;
    for (UINT i = 0; i < uiLoop; ++i)
    {
        const Real x = _deviceRandomF(fatIndex) * F(2.0) - F(1.0);
        const Real y = _deviceRandomF(fatIndex) * F(2.0) - F(1.0);
        if (x * x + y * y < F(1.0))
        {
            ++uiToAdd;
        }
        ++uiToAdd2;
    }
    sData1[fatIndex] = uiToAdd;
    sData2[fatIndex] = uiToAdd2;

    __syncthreads();
    if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {
        UINT all1 = 0;
        UINT all2 = 0;
        for (UINT i = 0; i < uithreadCount; ++i)
        {
            all1 += sData1[i];
            all2 += sData2[i];
        }
        //printf("how many?= %d\n", all1);
        atomicAdd(output, all1);
        atomicAdd(output + 1, all2);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelMCE(Real* output, UINT lengthyz, UINT lengthz, UINT uiLoop, UINT uithreadCount)
{
    __shared__ Real sData1[1024];
    __shared__ Real sData2[1024];
    Real fToAdd = 0;
    Real fToAdd2 = 0;
    const UINT fatIndex = threadIdx.x * lengthyz + threadIdx.y * lengthz + threadIdx.z;
    for (UINT i = 0; i < uiLoop; ++i)
    {
        const CLGComplex c = _deviceRandomGaussC(fatIndex);
        fToAdd += (c.x + c.y);
        fToAdd2 += (c.x * c.x + c.y * c.y);
    }
    sData1[fatIndex] = fToAdd;
    sData2[fatIndex] = fToAdd2;

    __syncthreads();
    if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {
        Real all1 = 0;
        Real all2 = 0;
        for (UINT i = 0; i < uithreadCount; ++i)
        {
            all1 += sData1[i];
            all2 += sData2[i];
        }
        //printf("how many?= %d\n", all1);
        atomicAdd(output, all1);
        atomicAdd(output + 1, all2);
    }
}

Real CLGAPI CalculatePi(const TArray<UINT> & decompose)
{
    dim3 blocknumber(decompose[0], decompose[1], decompose[2]);
    dim3 threadnumber(decompose[3], decompose[4], decompose[5]);
    const UINT threadCount = decompose[3] * decompose[4] * decompose[5];
    const UINT lengthyz = decompose[4] * decompose[5];
    const UINT lengthz = decompose[5];
    const UINT total = decompose[0] * decompose[1] * decompose[2] * decompose[3] * decompose[4] * decompose[5] * decompose[6];
    const UINT uiLoop = decompose[6];

    UINT outPutHost[2];
    outPutHost[0] = 0;
    outPutHost[1] = 0;

    UINT *outPut;
    checkCudaErrors(__cudaMalloc((void**)&outPut, sizeof(UINT) * 2));
    checkCudaErrors(cudaMemcpy(outPut, outPutHost, sizeof(UINT) * 2, cudaMemcpyHostToDevice));

    _LAUNCH_KERNEL(_kernelMCPi, blocknumber, threadnumber, outPut, lengthyz, lengthz, uiLoop, threadCount);
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaMemcpy(outPutHost, outPut, sizeof(UINT) * 2, cudaMemcpyDeviceToHost));

    appParanoiac(_T("==== results: %d / %d \n"), outPutHost[0], outPutHost[1]);

    return F(4.0) * outPutHost[0] / (Real)(total);
}

Real CLGAPI CalculateE(const TArray<UINT> & decompose)
{
    dim3 blocknumber(decompose[0], decompose[1], decompose[2]);
    dim3 threadnumber(decompose[3], decompose[4], decompose[5]);
    const UINT threadCount = decompose[3] * decompose[4] * decompose[5];
    const UINT lengthyz = decompose[4] * decompose[5];
    const UINT lengthz = decompose[5];
    const UINT total = decompose[0] * decompose[1] * decompose[2] * decompose[3] * decompose[4] * decompose[5] * decompose[6];
    const UINT uiLoop = decompose[6];

    Real outPutHost[2];
    outPutHost[0] = 0.0F;
    outPutHost[1] = 0.0F;

    Real *outPut;
    checkCudaErrors(__cudaMalloc((void**)&outPut, sizeof(Real) * 2));
    checkCudaErrors(cudaMemcpy(outPut, outPutHost, sizeof(Real) * 2, cudaMemcpyHostToDevice));

    _LAUNCH_KERNEL(_kernelMCE, blocknumber, threadnumber, outPut, lengthyz, lengthz, uiLoop, threadCount);
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaMemcpy(outPutHost, outPut, sizeof(Real) * 2, cudaMemcpyDeviceToHost));

    const Real fAv = outPutHost[0] / (2.0f * total);
    const Real fBv = outPutHost[1] / (2.0f * total) - fAv * fAv;

    return _hostsqrt(fBv);
}

#pragma endregion

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
