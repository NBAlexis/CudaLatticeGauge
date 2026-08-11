//=============================================================================
// FILENAME : SLNC.h
//
// DESCRIPTION:
// This is helper functions for calculate SL(N,C) complex matrix.
// SL(N,C) = N x N complex matrices with determinant = 1, no unitarity constraint.
//
// Layout (same as SUN.h):
// for N=4
// 0, 1, 2, 3
// 4, 5, 6, 7
// ...
//
// REVISION:
//  [05/29/2026 nbale]
//=============================================================================

#ifndef _SLNC_H_
#define _SLNC_H_

__BEGIN_NAMESPACE

template<INT N, INT NoE>
struct deviceSLNC
{

public:
    __device__ deviceSLNC()
    {

    }

    __device__ deviceSLNC(const deviceSLNC<N, NoE>& other)
    {
        memcpy(m_me, other.m_me, sizeof(CLGComplex) * N * N);
    }

    __device__ deviceSLNC(const CLGComplex* other)
    {
        memcpy(m_me, other, sizeof(CLGComplex) * N * N);
    }

    __device__ void DebugPrint(const char* header = NULL) const
    {
        printf("%s%s{{", NULL == header ? "" : header, NULL == header ? "" : "=");
        for (INT i = 0; i < N * N; ++i)
        {
            if ((N - 1) == (i % N))
            {
                printf(_CLGCMPFMT,
                    m_me[i].x,
                    m_me[i].y < 0 ? "" : "+",
                    m_me[i].y);
                if (i != (N * N - 1))
                {
                    printf("},\n{");
                }
            }
            else
            {
                printf(_CLGCMPFMT ",",
                    m_me[i].x,
                    m_me[i].y < 0 ? "" : "+",
                    m_me[i].y);
            }
        }
        printf("}};\n");
    }

#pragma region creation

    /**
    * ret = 0
    */
    __device__ __inline__ static deviceSLNC<N, NoE> makeSLNCZero()
    {
        deviceSLNC<N, NoE> ret;
        for (INT i = 0; i < N * N; ++i)
        {
            ret.m_me[i] = _make_cuComplex(F(0.0), F(0.0));
        }
        return ret;
    }

    __device__ __inline__ void Zero()
    {
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i] = _make_cuComplex(F(0.0), F(0.0));
        }
    }

    /**
    * ret = I
    */
    __device__ __inline__ static deviceSLNC<N, NoE> makeSLNCId()
    {
        deviceSLNC<N, NoE> ret;
        for (INT i = 0; i < N; ++i)
        {
            for (INT j = 0; j < N; ++j)
            {
                INT n = i * N + j;
                if (i == j)
                {
                    ret.m_me[n] = _make_cuComplex(F(1.0), F(0.0));
                }
                else
                {
                    ret.m_me[n] = _make_cuComplex(F(0.0), F(0.0));
                }
            }
        }
        return ret;
    }

    __device__ __inline__ void Id()
    {
        for (INT i = 0; i < N; ++i)
        {
            for (INT j = 0; j < N; ++j)
            {
                INT n = i * N + j;
                if (i == j)
                {
                    m_me[n] = _make_cuComplex(F(1.0), F(0.0));
                }
                else
                {
                    m_me[n] = _make_cuComplex(F(0.0), F(0.0));
                }
            }
        }
    }

    /**
    * The SU(N) generators form a basis of sl(N,C) over C.
    * For x < y: symmetric generator E_{yx} + E_{xy}
    * For x > y: antisymmetric generator i(E_{xy} - E_{yx})
    * For x == y: diagonal traceless generator
    */
    __device__ __inline__ static deviceSLNC<N, NoE> makeSLNCGenerator(UINT uiGenerator)
    {
        deviceSLNC<N, NoE> ret = makeSLNCZero();
        const Real half = F(1.0);
        if (uiGenerator < N * N - 1)
        {
            UINT x = uiGenerator % N;
            UINT y = uiGenerator / N;
            if (x < y)
            {
                ret.m_me[y * N + x] = _make_cuComplex(half, F(0.0));
                ret.m_me[x * N + y] = _make_cuComplex(half, F(0.0));
            }
            else if (x > y)
            {
                ret.m_me[y * N + x] = _make_cuComplex(F(0.0), -half);
                ret.m_me[x * N + y] = _make_cuComplex(F(0.0), half);
            }
            else
            {
                INT m = x + 1;
                Real fac = half * _sqrt(F(2.0) / (m * (m + F(1.0))));
                for (INT k = 0; k < m; ++k)
                {
                    ret.m_me[k * N + k] = _make_cuComplex(fac, F(0.0));
                }
                ret.m_me[m * N + m] = _make_cuComplex(fac * (-m), F(0.0));
            }
        }
        return ret;
    }

    /**
    * ret = c_a T_a, c_a are random complex Gaussians, T_a are generators.
    * This gives a random element of sl(N,C).
    */
    __device__ __inline__ static deviceSLNC<N, NoE> makeSLNCRandomGenerator(UINT fatIndex)
    {
        deviceSLNC<N, NoE> ret = makeSLNCZero();
        for (UINT i = 0; i < N * N - 1; ++i)
        {
            CLGComplex c;
            c.x = _deviceRandomGaussFSqrt2(fatIndex);
            c.y = _deviceRandomGaussFSqrt2(fatIndex);
            ret.Add(makeSLNCGenerator(i).MulCompC(c));
        }
        return ret;
    }

    __device__ __inline__ static deviceSLNC<N, NoE> makeSLNCSumGenerator(Real fDivide)
    {
        deviceSLNC<N, NoE> ret = makeSLNCZero();
        for (UINT i = 0; i < N * N - 1; ++i)
        {
            ret.Add(makeSLNCGenerator(i).MulRealC(fDivide));
        }
        return ret;
    }

    /**
    * can be called only after CLatticeData is created
    * ret = random complex entries, then normalized to det = 1
    */
    __device__ __inline__ static deviceSLNC<N, NoE> makeSLNCRandom(UINT fatIndex)
    {
        deviceSLNC<N, NoE> ret;
        for (UINT i = 0; i < N * N; ++i)
        {
            ret.m_me[i] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
        }
        ret.NormalizeDet();
        return ret;
    }

    /**
    * a random matrix which does not have to have det = 1
    */
    __device__ __inline__ static deviceSLNC<N, NoE> makeSLNCRandomAny(UINT fatIndex)
    {
        deviceSLNC<N, NoE> ret;
        for (UINT i = 0; i < N * N; ++i)
        {
            ret.m_me[i] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
        }
        return ret;
    }

    template<INT NoVE>
    __device__ __inline__ static deviceSLNC<N, NoE> makeSLNCContractV(const deviceSUNVector<N, NoVE>& left, const deviceSUNVector<N, NoVE>& right)
    {
        deviceSLNC<N, NoE> ret;
        for (INT y = 0; y < N; ++y)
        {
            for (INT x = 0; x < N; ++x)
            {
                ret.m_me[y * N + x] = _cuCmulf(_cuConjf(left.m_ve[x]), right.m_ve[y]);
            }
        }
        return ret;
    }

    //device SU3 SU2 is not implemented as device SUN Vector, so we need to specialize for them
    __device__ __inline__ static deviceSLNC<N, NoE> makeSLNCContractV(const deviceSU3Vector& left, const deviceSU3Vector& right)
    {
        if constexpr (3 == N)
        {
            deviceSLNC<N, NoE> ret;
            for (INT y = 0; y < N; ++y)
            {
                for (INT x = 0; x < N; ++x)
                {
                    ret.m_me[y * N + x] = _cuCmulf(_cuConjf(left.m_ve[x]), right.m_ve[y]);
                }
            }
            return ret;
        }
        return makeSLNCZero();
    }

    __device__ __inline__ static deviceSLNC<N, NoE> makeSLNCContractV(const deviceSU2Vector& left, const deviceSU2Vector& right)
    {
        if constexpr (2 == N)
        {
            deviceSLNC<N, NoE> ret;
            for (INT y = 0; y < N; ++y)
            {
                for (INT x = 0; x < N; ++x)
                {
                    ret.m_me[y * N + x] = _cuCmulf(_cuConjf(left.m_ve[x]), right.m_ve[y]);
                }
            }
            return ret;
        }
        return makeSLNCZero();
    }

#pragma endregion

#pragma region operators

    __device__ __inline__ void Add(const deviceSLNC<N, NoE>& right)
    {
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i] = _cuCaddf(m_me[i], right.m_me[i]);
        }
    }

    __device__ __inline__ void AddDagger(const deviceSLNC<N, NoE>& right)
    {
        for (INT y = 0; y < N; ++y)
        {
            for (INT x = 0; x < N; ++x)
            {
                m_me[y * N + x] = _cuCaddf(m_me[y * N + x], _cuConjf(right.m_me[x * N + y]));
            }
        }
    }

    __device__ __inline__ void Sub(const deviceSLNC<N, NoE>& right)
    {
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i] = _cuCsubf(m_me[i], right.m_me[i]);
        }
    }

    __device__ __inline__ void SubDagger(const deviceSLNC<N, NoE>& right)
    {
        for (INT y = 0; y < N; ++y)
        {
            for (INT x = 0; x < N; ++x)
            {
                m_me[y * N + x] = _cuCsubf(m_me[y * N + x], _cuConjf(right.m_me[x * N + y]));
            }
        }
    }

    __device__ __inline__ void AddReal(Real right)
    {
        for (INT i = 0; i < N; ++i)
        {
            m_me[i * N + i] = cuCaddf_cr(m_me[i * N + i], right);
        }
    }

    __device__ __inline__ void AddComp(const CLGComplex& right)
    {
        for (INT i = 0; i < N; ++i)
        {
            m_me[i * N + i] = _cuCaddf(m_me[i * N + i], right);
        }
    }

    __device__ __inline__ void SubReal(Real right)
    {
        for (INT i = 0; i < N; ++i)
        {
            m_me[i * N + i] = cuCsubf_cr(m_me[i * N + i], right);
        }
    }

    __device__ __inline__ void SubComp(const CLGComplex& right)
    {
        for (INT i = 0; i < N; ++i)
        {
            m_me[i * N + i] = _cuCsubf(m_me[i * N + i], right);
        }
    }

    __device__ __inline__ void Mul(const deviceSLNC<N, NoE>& right)
    {
        CLGComplex temp[N];
        for (INT y = 0; y < N; ++y)
        {
            for (INT x = 0; x < N; ++x)
            {
                temp[x] = _cuCmulf(m_me[y * N], right.m_me[x]);
                for (INT i = 1; i < N; ++i)
                {
                    temp[x] = _cuCaddf(temp[x], _cuCmulf(m_me[y * N + i], right.m_me[i * N + x]));
                }
            }
            memcpy(m_me + y * N, temp, sizeof(CLGComplex) * N);
        }
    }

    __device__ __inline__ void Mul(const deviceSLNC<N, NoE>& right, CLGComplex* temp)
    {
        for (INT y = 0; y < N; ++y)
        {
            for (INT x = 0; x < N; ++x)
            {
                temp[x] = _cuCmulf(m_me[y * N], right.m_me[x]);
                for (INT i = 1; i < N; ++i)
                {
                    temp[x] = _cuCaddf(temp[x], _cuCmulf(m_me[y * N + i], right.m_me[i * N + x]));
                }
            }
            memcpy(m_me + y * N, temp, sizeof(CLGComplex) * N);
        }
    }

    __device__ __inline__ void DaggerMul(const deviceSLNC<N, NoE>& right)
    {
        Dagger();
        Mul(right);
    }

    __device__ __inline__ void MulDagger(const deviceSLNC<N, NoE>& right)
    {
        CLGComplex temp[N];
        for (INT y = 0; y < N; ++y)
        {
            for (INT x = 0; x < N; ++x)
            {
                temp[x] = _cuCmulf(m_me[y * N], _cuConjf(right.m_me[x * N]));
                for (INT i = 1; i < N; ++i)
                {
                    temp[x] = _cuCaddf(temp[x], _cuCmulf(m_me[y * N + i], _cuConjf(right.m_me[x * N + i])));
                }
            }
            memcpy(m_me + y * N, temp, sizeof(CLGComplex) * N);
        }
    }

    __device__ __inline__ void MulDagger(const deviceSLNC<N, NoE>& right, CLGComplex* temp)
    {
        for (INT y = 0; y < N; ++y)
        {
            for (INT x = 0; x < N; ++x)
            {
                temp[x] = _cuCmulf(m_me[y * N], _cuConjf(right.m_me[x * N]));
                for (INT i = 1; i < N; ++i)
                {
                    temp[x] = _cuCaddf(temp[x], _cuCmulf(m_me[y * N + i], _cuConjf(right.m_me[x * N + i])));
                }
            }
            memcpy(m_me + y * N, temp, sizeof(CLGComplex) * N);
        }
    }

    __device__ __inline__ void MulOnMe(const deviceSLNC<N, NoE>& left)
    {
        CLGComplex temp[N];
        for (INT x = 0; x < N; ++x)
        {
            for (INT y = 0; y < N; ++y)
            {
                temp[y] = _cuCmulf(left.m_me[y * N], m_me[x]);
                for (INT i = 1; i < N; ++i)
                {
                    temp[y] = _cuCaddf(temp[y], _cuCmulf(left.m_me[y * N + i], m_me[i * N + x]));
                }
            }

            for (INT y = 0; y < N; ++y)
            {
                m_me[y * N + x] = temp[y];
            }
        }
    }

    __device__ __inline__ void MulOnMeDN(const deviceSLNC<N, NoE>& left)
    {
        CLGComplex temp[N];
        for (INT x = 0; x < N; ++x)
        {
            for (INT y = 0; y < N; ++y)
            {
                temp[y] = _cuCmulf(_cuConjf(left.m_me[y]), m_me[x]);
                for (INT i = 1; i < N; ++i)
                {
                    temp[y] = _cuCaddf(temp[y], _cuCmulf(_cuConjf(left.m_me[i * N + y]), m_me[i * N + x]));
                }
            }

            for (INT y = 0; y < N; ++y)
            {
                m_me[y * N + x] = temp[y];
            }
        }
    }

    __device__ __inline__ void MulOnMeND(const deviceSLNC<N, NoE>& left)
    {
        Dagger();
        MulOnMe(left);
    }

    __device__ __inline__ deviceSLNC<N, NoE> MulC(const deviceSLNC<N, NoE>& right) const
    {
        deviceSLNC<N, NoE> ret;
        for (INT x = 0; x < N; ++x)
        {
            for (INT y = 0; y < N; ++y)
            {
                ret.m_me[y * N + x] = _cuCmulf(m_me[y * N], right.m_me[x]);
                for (INT i = 1; i < N; ++i)
                {
                    ret.m_me[y * N + x] = _cuCaddf(ret.m_me[y * N + x], _cuCmulf(m_me[y * N + i], right.m_me[i * N + x]));
                }
            }
        }

        return ret;
    }

    __device__ __inline__ void MulReal(Real right)
    {
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i] = cuCmulf_cr(m_me[i], right);
        }
    }

    __device__ __inline__ void MulComp(const CLGComplex& right)
    {
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i] = _cuCmulf(m_me[i], right);
        }
    }

    __device__ __inline__ void DivComp(const CLGComplex& right)
    {
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i] = _cuCdivf(m_me[i], right);
        }
    }

    __device__ __inline__ void DivReal(const Real& right)
    {
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i] = cuCdivf_cr(m_me[i], right);
        }
    }

    __device__ __inline__ deviceSLNC<N, NoE> MulDaggerC(const deviceSLNC<N, NoE>& right) const
    {
        deviceSLNC<N, NoE> ret;
        for (INT x = 0; x < N; ++x)
        {
            for (INT y = 0; y < N; ++y)
            {
                ret.m_me[y * N + x] = _cuCmulf(m_me[y * N], _cuConjf(right.m_me[x * N]));
                for (INT i = 1; i < N; ++i)
                {
                    ret.m_me[y * N + x] = _cuCaddf(ret.m_me[y * N + x], _cuCmulf(m_me[y * N + i], _cuConjf(right.m_me[x * N + i])));
                }
            }
        }

        return ret;
    }

    __device__ __inline__ deviceSLNC<N, NoE> DaggerMulC(const deviceSLNC<N, NoE>& right) const
    {
        deviceSLNC<N, NoE> ret;
        for (INT x = 0; x < N; ++x)
        {
            for (INT y = 0; y < N; ++y)
            {
                ret.m_me[y * N + x] = _cuCmulf(_cuConjf(m_me[y]), right.m_me[x]);
                for (INT i = 1; i < N; ++i)
                {
                    ret.m_me[y * N + x] = _cuCaddf(ret.m_me[y * N + x], _cuCmulf(_cuConjf(m_me[i * N + y]), right.m_me[i * N + x]));
                }
            }
        }

        return ret;
    }

    __device__ __inline__ deviceSLNC<N, NoE> AddC(const deviceSLNC<N, NoE>& right) const { deviceSLNC<N, NoE> ret(*this); ret.Add(right); return ret; }
    __device__ __inline__ deviceSLNC<N, NoE> AddCompC(const CLGComplex& right) const { deviceSLNC<N, NoE> ret(*this); ret.AddComp(right); return ret; }
    __device__ __inline__ deviceSLNC<N, NoE> AddRealC(const Real& right) const { deviceSLNC<N, NoE> ret(*this); ret.AddReal(right); return ret; }

    __device__ __inline__ deviceSLNC<N, NoE> SubC(const deviceSLNC<N, NoE>& right) const { deviceSLNC<N, NoE> ret(*this); ret.Sub(right); return ret; }
    __device__ __inline__ deviceSLNC<N, NoE> SubCompC(const CLGComplex& right) const { deviceSLNC<N, NoE> ret(*this); ret.SubComp(right); return ret; }
    __device__ __inline__ deviceSLNC<N, NoE> SubRealC(const Real& right) const { deviceSLNC<N, NoE> ret(*this); ret.SubReal(right); return ret; }

    __device__ __inline__ deviceSLNC<N, NoE> MulCompC(const CLGComplex& right) const { deviceSLNC<N, NoE> ret(*this); ret.MulComp(right); return ret; }
    __device__ __inline__ deviceSLNC<N, NoE> MulRealC(const Real& right) const { deviceSLNC<N, NoE> ret(*this); ret.MulReal(right); return ret; }
    __device__ __inline__ deviceSLNC<N, NoE> DivCompC(const CLGComplex& right) const { deviceSLNC<N, NoE> ret(*this); ret.DivComp(right); return ret; }
    __device__ __inline__ deviceSLNC<N, NoE> DivRealC(const Real& right) const { deviceSLNC<N, NoE> ret(*this); ret.DivReal(right); return ret; }

    /**
    * self block mult must be square
    */
    __device__ __inline__ void BlockMul(
        const deviceSLNC<N, NoE>& right,
        UINT iWidth,
        UINT iLeftStartX, UINT iLeftStartY,
        UINT iRightStartX, UINT iRightStartY)
    {
        CLGComplex temp[N];
        for (INT y = 0; y < iWidth; ++y)
        {
            for (INT x = 0; x < iWidth; ++x)
            {
                temp[x] = _cuCmulf(m_me[(y + iLeftStartY) * N + iLeftStartX], right.m_me[iRightStartY * N + x + iRightStartX]);
                for (INT i = 1; i < iWidth; ++i)
                {
                    temp[x] = _cuCaddf(temp[x], _cuCmulf(m_me[(y + iLeftStartY) * N + (i + iLeftStartX)], right.m_me[(i + iRightStartY) * N + x + iRightStartX]));
                }
            }
            memcpy(m_me + (y + iLeftStartY) * N + iLeftStartX, temp, sizeof(CLGComplex) * iWidth);
        }
    }

    __device__ __inline__ void BlockCopy(const deviceSLNC<N, NoE>& right,
        UINT iLeftStartX, UINT iLeftStartY,
        UINT iRightStartX, UINT iRightStartY,
        UINT iWidth, UINT iHeight)
    {
        for (INT x = 0; x < iWidth; ++x)
        {
            for (INT y = 0; y < iHeight; ++y)
            {
                m_me[(y + iLeftStartY) * N + x + iLeftStartX] = right.m_me[(y + iRightStartY) * N + x + iRightStartX];
            }
        }
    }

    template<INT NofVE>
    __device__ __inline__ deviceSUNVector<N, NofVE> MulVector(const deviceSUNVector<N, NofVE>& v) const
    {
        deviceSUNVector<N, NofVE> ret;
        for (INT y = 0; y < N; ++y)
        {
            ret.m_ve[y] = _cuCmulf(m_me[y * N], v.m_ve[0]);
            for (INT x = 1; x < N; ++x)
            {
                ret.m_ve[y] = _cuCaddf(ret.m_ve[y], _cuCmulf(m_me[y * N + x], v.m_ve[x]));
            }
        }
        return ret;
    }

    //device SU3 SU2 is not implemented as device SUN Vector, so we need to specialize for them
    __device__ __inline__ deviceSU3Vector MulVector(const deviceSU3Vector& v) const
    {
        if constexpr (3 == N)
        {
            deviceSU3Vector ret;
            for (INT y = 0; y < N; ++y)
            {
                ret.m_ve[y] = _cuCmulf(m_me[y * N], v.m_ve[0]);
                for (INT x = 1; x < N; ++x)
                {
                    ret.m_ve[y] = _cuCaddf(ret.m_ve[y], _cuCmulf(m_me[y * N + x], v.m_ve[x]));
                }
            }
            return ret;
        }
        return deviceSU3Vector::makeZeroSU3Vector();
    }

    __device__ __inline__ deviceSU2Vector MulVector(const deviceSU2Vector& v) const
    {
        if constexpr (2 == N)
        {
            deviceSU2Vector ret;
            for (INT y = 0; y < N; ++y)
            {
                ret.m_ve[y] = _cuCmulf(m_me[y * N], v.m_ve[0]);
                for (INT x = 1; x < N; ++x)
                {
                    ret.m_ve[y] = _cuCaddf(ret.m_ve[y], _cuCmulf(m_me[y * N + x], v.m_ve[x]));
                }
            }
            return ret;
        }
        return deviceSU2Vector::makeZeroSU2Vector();
    }

    template<INT NofVE>
    __device__ __inline__ deviceSUNVector<N, NofVE> DagMulVector(const deviceSUNVector<N, NofVE>& v) const
    {
        deviceSUNVector<N, NofVE> ret;
        for (INT y = 0; y < N; ++y)
        {
            ret.m_ve[y] = _cuCmulf(_cuConjf(m_me[y]), v.m_ve[0]);
            for (INT x = 1; x < N; ++x)
            {
                ret.m_ve[y] = _cuCaddf(ret.m_ve[y], _cuCmulf(_cuConjf(m_me[x * N + y]), v.m_ve[x]));
            }
        }
        return ret;
    }

    //device SU3 SU2 is not implemented as device SUN Vector, so we need to specialize for them
    __device__ __inline__ deviceSU3Vector DagMulVector(const deviceSU3Vector& v) const
    {
        if constexpr (3 == N)
        {
            deviceSU3Vector ret;
            for (INT y = 0; y < N; ++y)
            {
                ret.m_ve[y] = _cuCmulf(_cuConjf(m_me[y]), v.m_ve[0]);
                for (INT x = 1; x < N; ++x)
                {
                    ret.m_ve[y] = _cuCaddf(ret.m_ve[y], _cuCmulf(_cuConjf(m_me[x * N + y]), v.m_ve[x]));
                }
            }
            return ret;
        }
        return deviceSU3Vector::makeZeroSU3Vector();
    }

    __device__ __inline__ deviceSU2Vector DagMulVector(const deviceSU2Vector& v) const
    {
        if constexpr (2 == N)
        {
            deviceSU2Vector ret;
            for (INT y = 0; y < N; ++y)
            {
                ret.m_ve[y] = _cuCmulf(_cuConjf(m_me[y]), v.m_ve[0]);
                for (INT x = 1; x < N; ++x)
                {
                    ret.m_ve[y] = _cuCaddf(ret.m_ve[y], _cuCmulf(_cuConjf(m_me[x * N + y]), v.m_ve[x]));
                }
            }
            return ret;
        }
        return deviceSU2Vector::makeZeroSU2Vector();
    }

    __device__ __inline__ void MulVectorElement(CLGComplex* res, const CLGComplex* v) const
    {
        for (INT y = 0; y < N; ++y)
        {
            res[y] = _cuCmulf(m_me[y * N], v[0]);
            for (INT x = 1; x < N; ++x)
            {
                res[y] = _cuCaddf(res[y], _cuCmulf(m_me[y * N + x], v[x]));
            }
        }
    }

#pragma endregion

#pragma region useful functions

public:

    __device__ __inline__ CLGComplex Determinent() const
    {
        deviceSLNC<N, NoE> tmp(*this);
        tmp.LUNoReturn();
        for (INT i = 1; i < N; ++i)
        {
            tmp.m_me[0] = _cuCmulf(tmp.m_me[0], tmp.m_me[i * N + i]);
        }
        return tmp.m_me[0];
    }

    /**
    * Normalize determinant to 1 by dividing each element by det^(1/N)
    */
    __device__ __inline__ void NormalizeDet()
    {
        CLGComplex det = Determinent();
        CLGComplex factor = __cuCpowerf(det, F(-1.0) / N);
        MulComp(factor);
    }

    __device__ __inline__ void Traceless()
    {
        CLGComplex tr = Tr();
        tr = cuCdivf_cr(tr, static_cast<Real>(N));
        for (INT i = 0; i < N; ++i)
        {
            m_me[i * N + i] = _cuCsubf(m_me[i * N + i], tr);
        }
    }

    __device__ __inline__ deviceSLNC<N, NoE> TracelessC() const
    {
        deviceSLNC<N, NoE> ret(*this);
        ret.Traceless();
        return ret;
    }

    __device__ __inline__ void Ta()
    {
        Real trace = F(0.0);
        for (INT y = 0; y < N; ++y)
        {
            for (INT x = y; x < N; ++x)
            {
                if (x == y)
                {
                    m_me[y * N + x].x = F(0.0);
                    trace += m_me[y * N + x].y;
                }
                else
                {
                    m_me[y * N + x].x = F(0.5) * (m_me[y * N + x].x - m_me[x * N + y].x);
                    m_me[y * N + x].y = F(0.5) * (m_me[y * N + x].y + m_me[x * N + y].y);
                    m_me[x * N + y].x = -m_me[y * N + x].x;
                    m_me[x * N + y].y = m_me[y * N + x].y;
                }
            }
        }

        trace = trace / N;
        for (INT x = 0; x < N; ++x)
        {
            m_me[x * N + x].y = m_me[x * N + x].y - trace;
        }
    }

    __device__ __inline__ CLGComplex Tr() const
    {
        CLGComplex ret = m_me[0];
        for (INT i = 1; i < N; ++i)
        {
            ret = _cuCaddf(ret, m_me[i * N + i]);
        }
        return ret;
    }

    /**
    * Re[Tr[U]]
    */
    __device__ __inline__ Real ReTr() const
    {
        Real ret = m_me[0].x;
        for (INT i = 1; i < N; ++i)
        {
            ret += m_me[i * N + i].x;
        }
        return ret;
    }

    __device__ __inline__ Real ImTr() const
    {
        Real ret = m_me[0].y;
        for (INT i = 1; i < N; ++i)
        {
            ret += m_me[i * N + i].y;
        }
        return ret;
    }

    __device__ __inline__ void Re()
    {
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i].y = F(0.0);
        }
    }

    __device__ __inline__ deviceSLNC<N, NoE> ReC() const { deviceSLNC<N, NoE> ret(*this); ret.Re(); return ret; }

    __device__ __inline__ void Im()
    {
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i].x = m_me[i].y;
            m_me[i].y = F(0.0);
        }
    }

    __device__ __inline__ deviceSLNC<N, NoE> ImC() const { deviceSLNC<N, NoE> ret(*this); ret.Im(); return ret; }

    /**
    * res = Conjugate[Transpose[U]]
    */
    __device__ __inline__ void Dagger()
    {
        for (INT x = 0; x < N; ++x)
        {
            for (INT y = x; y < N; ++y)
            {
                if (x == y)
                {
                    m_me[x * N + x].y = -m_me[x * N + x].y;
                }
                else
                {
                    CLGComplex tmp = m_me[y * N + x];
                    m_me[y * N + x] = _cuConjf(m_me[x * N + y]);
                    m_me[x * N + y] = _cuConjf(tmp);
                }
            }
        }
    }

    __device__ __inline__ deviceSLNC<N, NoE> DaggerC() const
    {
        deviceSLNC<N, NoE> ret;
        for (INT x = 0; x < N; ++x)
        {
            for (INT y = 0; y < N; ++y)
            {
                ret.m_me[y * N + x] = _cuConjf(m_me[x * N + y]);
            }
        }
        return ret;
    }

    __device__ __inline__ void Opposite()
    {
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i].x = -m_me[i].x;
            m_me[i].y = -m_me[i].y;
        }
    }

    __device__ __inline__ deviceSLNC<N, NoE> OppositeC() const
    {
        deviceSLNC<N, NoE> ret(*this);
        ret.Opposite();
        return ret;
    }

    __device__ __inline__ deviceSLNC<N, NoE> Transpose() const
    {
        deviceSLNC<N, NoE> ret;
        for (INT x = 0; x < N; ++x)
        {
            for (INT y = 0; y < N; ++y)
            {
                ret.m_me[y * N + x] = m_me[x * N + y];
            }
        }
        return ret;
    }

    /**
    * U = U + 1
    */
    __device__ __inline__ void AddId()
    {
        for (INT i = 0; i < N; ++i)
        {
            m_me[i * N + i].x += F(1.0);
        }
    }

    /**
    * U' = exp(aU) = (1 + a U + a^2 U^2/2 +  ... + a^N U^N/N!)
    *    = 1 + a U (1 + a U /2 (1 + a U/3 ...))
    */
    __device__ __inline__ deviceSLNC<N, NoE> Exp(const CLGComplex& a, BYTE uiPrecision = N + 1) const
    {
        deviceSLNC<N, NoE> tmp;

        for (BYTE i = 0; i < uiPrecision; ++i)
        {
            const Real exp_factor = __rcp(uiPrecision - i);
            CLGComplex alpha = cuCmulf_cr(a, exp_factor);
            deviceSLNC<N, NoE> aUoN(*this);
            aUoN.MulComp(alpha);
            if (0 == i)
            {
                tmp = aUoN;
            }
            else
            {
                tmp.Mul(aUoN);
            }
            tmp.AddId();
        }

        tmp.NormalizeDet();
        return tmp;
    }

    __device__ __inline__ deviceSLNC<N, NoE> ExpReal(Real a, BYTE uiPrecision = N + 1) const
    {
        deviceSLNC<N, NoE> tmp;
        for (BYTE i = 0; i < uiPrecision; ++i)
        {
            deviceSLNC<N, NoE> aUoN(*this);
            aUoN.MulReal(a * __rcp(uiPrecision - i));
            if (0 == i)
            {
                tmp = aUoN;
            }
            else
            {
                tmp.Mul(aUoN);
            }
            tmp.AddId();
        }
        tmp.NormalizeDet();
        return tmp;
    }

private:

    __device__ __inline__ void OneStepHouseHolderQR(deviceSLNC<N, NoE>& A, UINT i)
    {
        Real fLength = F(0.0);
        for (INT y = i; y < N; ++y)
        {
            fLength += __cuCabsSqf(A.m_me[y * N + i]);
        }

        if (fLength > _CLG_FLT_MIN_)
        {
            fLength = _sqrt(fLength);
            CLGComplex vlst[N];
            Real fULength = F(0.0);
            for (INT y = i; y < N; ++y)
            {
                vlst[y] = A.m_me[y * N + i];
                if (y == i)
                {
                    const Real arg = __cuCargf(vlst[y]);
                    vlst[y] = _cuCaddf(vlst[y], _make_cuComplex(_cos(arg) * fLength, _sin(arg) * fLength));
                }
                fULength += __cuCabsSqf(vlst[y]);
            }

            if (fULength > _CLG_FLT_MIN_)
            {
                const Real fBeta = F(2.0) / fULength;

                if (0 == i)
                {
                    for (INT px = i; px < N; ++px)
                    {
                        for (INT py = i; py < N; ++py)
                        {
                            m_me[py * N + px] = _cuCmulf(vlst[py], _cuConjf(vlst[px]));
                            if (px == py)
                            {
                                m_me[py * N + px].x = F(1.0) - fBeta * m_me[py * N + px].x;
                            }
                            else
                            {
                                m_me[py * N + px].x = -fBeta * m_me[py * N + px].x;
                            }
                            m_me[py * N + px].y = -fBeta * m_me[py * N + px].y;
                        }
                    }

                    A.MulOnMe(*this);
                }
                else
                {
                    deviceSLNC<N, NoE> tmp = makeSLNCId();
                    for (INT px = i; px < N; ++px)
                    {
                        for (INT py = i; py < N; ++py)
                        {
                            tmp.m_me[py * N + px] = _cuCmulf(vlst[py], _cuConjf(vlst[px]));
                            if (px == py)
                            {
                                tmp.m_me[py * N + px].x = F(1.0) - fBeta * tmp.m_me[py * N + px].x;
                            }
                            else
                            {
                                tmp.m_me[py * N + px].x = -fBeta * tmp.m_me[py * N + px].x;
                            }
                            tmp.m_me[py * N + px].y = -fBeta * tmp.m_me[py * N + px].y;
                        }
                    }

                    A.MulOnMe(tmp);
                    Mul(tmp);
                }
            }
        }
    }

public:

    /**
    * Let me be Q and return R
    */
    __device__ __inline__ deviceSLNC<N, NoE> QR()
    {
        deviceSLNC<N, NoE> r(*this);
        for (UINT i = 0; i < (N - 1); ++i)
        {
            OneStepHouseHolderQR(r, i);
        }
        return r;
    }

private:

    __device__ __inline__ void OneStepGaussianLU(deviceSLNC<N, NoE>& L, UINT i)
    {
        CLGComplex dobedivide = m_me[i * N + i];
        if (__cuCabsSqf(dobedivide) > _CLG_FLT_MIN_)
        {
            for (INT y = i + 1; y < N; ++y)
            {
                L.m_me[y * N + i] = _cuCdivf(m_me[y * N + i], dobedivide);
                m_me[y * N + i] = _zeroc;
            }

            for (INT x = i + 1; x < N; ++x)
            {
                for (INT y = i + 1; y < N; ++y)
                {
                    CLGComplex tobesub = _cuCmulf(L.m_me[y * N + i], m_me[i * N + x]);
                    m_me[y * N + x] = _cuCsubf(m_me[y * N + x], tobesub);
                }
            }
        }
    }

    __device__ __inline__ void OneStepGaussianLUNoReturn(UINT i)
    {
        CLGComplex dobedivide = m_me[i * N + i];
        if (__cuCabsSqf(dobedivide) > _CLG_FLT_MIN_)
        {
            for (INT y = i + 1; y < N; ++y)
            {
                m_me[y * N + i] = _cuCdivf(m_me[y * N + i], dobedivide);
            }

            for (INT x = i + 1; x < N; ++x)
            {
                for (INT y = i + 1; y < N; ++y)
                {
                    CLGComplex tobesub = _cuCmulf(m_me[y * N + i], m_me[i * N + x]);
                    m_me[y * N + x] = _cuCsubf(m_me[y * N + x], tobesub);
                }
            }
        }
    }

public:

    /**
    * return L, and let me be U
    */
    __device__ __inline__ deviceSLNC<N, NoE> LU()
    {
        deviceSLNC<N, NoE> L = makeSLNCId();
        for (UINT i = 0; i < (N - 1); ++i)
        {
            OneStepGaussianLU(L, i);
        }
        return L;
    }

    __device__ __inline__ void LUNoReturn()
    {
        for (UINT i = 0; i < (N - 1); ++i)
        {
            OneStepGaussianLUNoReturn(i);
        }
    }

    /**
    * let me be R matrix, set v=R^{-1}.v
    */
    template<INT NoVE>
    __device__ __inline__ void BackwardSubstitution(deviceSUNVector<N, NoVE>& v) const
    {
        for (INT i = N - 1; i >= 0; --i)
        {
            for (INT j = i + 1; j < N; ++j)
            {
                v.m_ve[i] = _cuCsubf(v.m_ve[i], _cuCmulf(m_me[i * N + j], v.m_ve[j]));
            }
            if (__cuCabsSqf(m_me[i * N + i]) > _CLG_FLT_MIN_)
            {
                v.m_ve[i] = _cuCdivf(v.m_ve[i], m_me[i * N + i]);
            }
        }
    }

    __device__ __inline__ deviceSLNC<N, NoE> InverseC() const
    {
        deviceSLNC<N, NoE> a(*this);
        deviceSLNC<N, NoE> inv = makeSLNCId();

        for (INT i = 0; i < N; ++i)
        {
            INT pivot = i;
            Real best = __cuCabsSqf(a.m_me[i * N + i]);
            for (INT y = i + 1; y < N; ++y)
            {
                const Real candidate = __cuCabsSqf(a.m_me[y * N + i]);
                if (candidate > best)
                {
                    best = candidate;
                    pivot = y;
                }
            }

            if (pivot != i)
            {
                for (INT x = 0; x < N; ++x)
                {
                    CLGComplex tmp = a.m_me[i * N + x];
                    a.m_me[i * N + x] = a.m_me[pivot * N + x];
                    a.m_me[pivot * N + x] = tmp;

                    tmp = inv.m_me[i * N + x];
                    inv.m_me[i * N + x] = inv.m_me[pivot * N + x];
                    inv.m_me[pivot * N + x] = tmp;
                }
            }

            const CLGComplex divider = a.m_me[i * N + i];
            for (INT x = 0; x < N; ++x)
            {
                a.m_me[i * N + x] = _cuCdivf(a.m_me[i * N + x], divider);
                inv.m_me[i * N + x] = _cuCdivf(inv.m_me[i * N + x], divider);
            }

            for (INT y = 0; y < N; ++y)
            {
                if (y == i)
                {
                    continue;
                }

                const CLGComplex factor = a.m_me[y * N + i];
                if (__cuCabsSqf(factor) <= _CLG_FLT_MIN_)
                {
                    continue;
                }

                for (INT x = 0; x < N; ++x)
                {
                    a.m_me[y * N + x] = _cuCsubf(a.m_me[y * N + x], _cuCmulf(factor, a.m_me[i * N + x]));
                    inv.m_me[y * N + x] = _cuCsubf(inv.m_me[y * N + x], _cuCmulf(factor, inv.m_me[i * N + x]));
                }
            }
        }

        return inv;
    }

    __device__ __inline__ void BackwardSubstitutionElement(CLGComplex* v) const
    {
        for (INT i = N - 1; i >= 0; --i)
        {
            for (INT j = i + 1; j < N; ++j)
            {
                v[i] = _cuCsubf(v[i], _cuCmulf(m_me[i * N + j], v[j]));
            }
            if (__cuCabsSqf(m_me[i * N + i]) > _CLG_FLT_MIN_)
            {
                v[i] = _cuCdivf(v[i], m_me[i * N + i]);
            }
        }
    }

private:

    __device__ __inline__ void OneStepHouseHolderHessenberg(deviceSLNC<N, NoE>& A, UINT i)
    {
        Real fLength = F(0.0);
        for (INT y = i + 1; y < N; ++y)
        {
            fLength += __cuCabsSqf(A.m_me[y * N + i]);
        }

        if (fLength > _CLG_FLT_MIN_)
        {
            fLength = _sqrt(fLength);
            CLGComplex vlst[N];
            Real fULength = F(0.0);
            for (INT y = i + 1; y < N; ++y)
            {
                vlst[y] = A.m_me[y * N + i];
                if (y == i + 1)
                {
                    const Real arg = __cuCargf(vlst[y]);
                    vlst[y] = _cuCaddf(vlst[y], _make_cuComplex(_cos(arg) * fLength, _sin(arg) * fLength));
                }
                fULength += __cuCabsSqf(vlst[y]);
            }

            if (fULength > _CLG_FLT_MIN_)
            {
                const Real fBeta = F(2.0) / fULength;

                if (0 == i)
                {
                    for (INT px = 0; px < N; ++px)
                    {
                        for (INT py = 0; py < N; ++py)
                        {
                            if (0 == px && 0 == py)
                            {
                                m_me[py * N + px] = _onec;
                            }
                            else if (px > 0 && py > 0)
                            {
                                m_me[py * N + px] = _cuCmulf(vlst[py], _cuConjf(vlst[px]));
                                if (px == py)
                                {
                                    m_me[py * N + px].x = F(1.0) - fBeta * m_me[py * N + px].x;
                                }
                                else
                                {
                                    m_me[py * N + px].x = -fBeta * m_me[py * N + px].x;
                                }
                                m_me[py * N + px].y = -fBeta * m_me[py * N + px].y;
                            }
                            else
                            {
                                m_me[py * N + px] = _zeroc;
                            }
                        }
                    }

                    A.MulOnMe(*this);
                    A.Mul(*this);
                }
                else
                {
                    deviceSLNC<N, NoE> tmp = makeSLNCId();
                    for (INT px = i + 1; px < N; ++px)
                    {
                        for (INT py = i + 1; py < N; ++py)
                        {
                            tmp.m_me[py * N + px] = _cuCmulf(vlst[py], _cuConjf(vlst[px]));
                            if (px == py)
                            {
                                tmp.m_me[py * N + px].x = F(1.0) - fBeta * tmp.m_me[py * N + px].x;
                            }
                            else
                            {
                                tmp.m_me[py * N + px].x = -fBeta * tmp.m_me[py * N + px].x;
                            }
                            tmp.m_me[py * N + px].y = -fBeta * tmp.m_me[py * N + px].y;
                        }
                    }

                    A.MulOnMe(tmp);
                    A.Mul(tmp);
                    Mul(tmp);
                }
            }
        }
    }

    __device__ __inline__ void OneStepHouseHolderHessenbergC(deviceSLNC<N, NoE>& A, UINT i) const
    {
        Real fLength = F(0.0);
        for (INT y = i + 1; y < N; ++y)
        {
            fLength += __cuCabsSqf(A.m_me[y * N + i]);
        }

        if (fLength > _CLG_FLT_MIN_)
        {
            fLength = _sqrt(fLength);
            CLGComplex vlst[N];
            Real fULength = F(0.0);
            for (INT y = i + 1; y < N; ++y)
            {
                vlst[y] = A.m_me[y * N + i];
                if (y == i + 1)
                {
                    const Real arg = __cuCargf(vlst[y]);
                    vlst[y] = _cuCaddf(vlst[y], _make_cuComplex(_cos(arg) * fLength, _sin(arg) * fLength));
                }
                fULength += __cuCabsSqf(vlst[y]);
            }

            if (fULength > _CLG_FLT_MIN_)
            {
                const Real fBeta = F(2.0) / fULength;
                deviceSLNC<N, NoE> tmp = makeSLNCId();
                for (INT px = i + 1; px < N; ++px)
                {
                    for (INT py = i + 1; py < N; ++py)
                    {
                        tmp.m_me[py * N + px] = _cuCmulf(vlst[py], _cuConjf(vlst[px]));
                        if (px == py)
                        {
                            tmp.m_me[py * N + px].x = F(1.0) - fBeta * tmp.m_me[py * N + px].x;
                        }
                        else
                        {
                            tmp.m_me[py * N + px].x = -fBeta * tmp.m_me[py * N + px].x;
                        }
                        tmp.m_me[py * N + px].y = -fBeta * tmp.m_me[py * N + px].y;
                    }
                }

                A.MulOnMe(tmp);
                A.Mul(tmp);
            }
        }
    }

public:

    __device__ __inline__ deviceSLNC<N, NoE> Hessenberg()
    {
        deviceSLNC<N, NoE> r(*this);
        for (UINT i = 0; i < (N - 2); ++i)
        {
            OneStepHouseHolderHessenberg(r, i);
        }
        return r;
    }

    __device__ __inline__ deviceSLNC<N, NoE> HessenbergC() const
    {
        deviceSLNC<N, NoE> r(*this);
        for (UINT i = 0; i < (N - 2); ++i)
        {
            OneStepHouseHolderHessenbergC(r, i);
        }
        return r;
    }

    #pragma region QR Iteration

    static __device__ __inline__ void checkMatrixIndexDoubleShift(CLGComplex* mtr, INT* decomp, INT dx)
    {
        decomp[0] = 0;
        decomp[1] = dx;

        for (INT i = dx - 2; i >= 0; --i)
        {
            if (__cuCabsSqf(mtr[(i + 1) * dx + i]) < _CLG_FLT_MIN_)
            {
                mtr[(i + 1) * dx + i] = _make_cuComplex(F(0.0), F(0.0));

                if (decomp[1] == i + 2)
                {
                    decomp[1] = i + 1;
                }

                if (i + 1 > decomp[0] && i + 1 < decomp[1])
                {
                    decomp[0] = i + 1;
                }
            }
        }
    }

    static __device__ __inline__ void doubleEigen(CLGComplex& s, CLGComplex& t,
        const CLGComplex& h00, const CLGComplex& h01,
        const CLGComplex& h10, const CLGComplex& h11)
    {
        s = _cuCaddf(h11, h00);
        t = _cuCsubf(_cuCmulf(h00, h11), _cuCmulf(h10, h01));
    }

    static __device__ __inline__  void doubleShift(
        CLGComplex& a, CLGComplex& b, CLGComplex& c,
        const CLGComplex& s, const CLGComplex& t,
        const CLGComplex& h00, const CLGComplex& h01,
        const CLGComplex& h10, const CLGComplex& h11,
        const CLGComplex& h21)
    {
        a = _cuCaddf(
            _cuCmulf(h00, _cuCsubf(h00, s)),
            _cuCaddf(_cuCmulf(h10, h01), t));
        b = _cuCmulf(h10,
            _cuCsubf(_cuCaddf(h00, h11), s)
        );
        c = _cuCmulf(h10, h21);
    }

    static __device__ __inline__ void threeHouseHolder(CLGComplex& a, CLGComplex& b, CLGComplex& c)
    {
        Real len = __cuCabsSqf(a) + __cuCabsSqf(b) + __cuCabsSqf(c);
        if (len < _CLG_FLT_MIN_)
        {
            a = _make_cuComplex(F(0.0), F(0.0));
            b = _make_cuComplex(F(0.0), F(0.0));
            c = _make_cuComplex(F(0.0), F(0.0));
            return;
        }
        len = _sqrt(len);
        Real lena = a.x * a.x + a.y * a.y;
        Real fCos = F(0.0);
        Real fSin = F(0.0);
        if (lena < _CLG_FLT_MIN_)
        {
            const Real fArg = _atan2(a.y, a.x);
            fCos = _cos(fArg);
            fSin = _sin(fArg);
        }
        else
        {
            lena = __div(F(1.0), _sqrt(lena));
            fCos = a.x * lena;
            fSin = a.y * lena;
        }
        a = _cuCaddf(a, _make_cuComplex(len * fCos, len * fSin));

        Real len2 = F(0.5) * (a.x * a.x + a.y * a.y + b.x * b.x + b.y * b.y + c.x * c.x + c.y * c.y);
        if (len2 < _CLG_FLT_MIN_)
        {
            a = _make_cuComplex(F(0.0), F(0.0));
            b = _make_cuComplex(F(0.0), F(0.0));
            c = _make_cuComplex(F(0.0), F(0.0));
            return;
        }
        len2 = __div(F(1.0), _sqrt(len2));
        a.x = a.x * len2;
        a.y = a.y * len2;
        b.x = b.x * len2;
        b.y = b.y * len2;
        c.x = c.x * len2;
        c.y = c.y * len2;
    }

    static __device__ __inline__ void twoHouseHolder(CLGComplex& a, CLGComplex& b)
    {
        Real len = __cuCabsSqf(a) + __cuCabsSqf(b);
        if (len < _CLG_FLT_MIN_)
        {
            a = _make_cuComplex(F(0.0), F(0.0));
            b = _make_cuComplex(F(0.0), F(0.0));
            return;
        }
        len = _sqrt(len);
        Real lena = a.x * a.x + a.y * a.y;
        Real fCos = F(0.0);
        Real fSin = F(0.0);
        if (lena < _CLG_FLT_MIN_)
        {
            const Real fArg = _atan2(a.y, a.x);
            fCos = _cos(fArg);
            fSin = _sin(fArg);
        }
        else
        {
            lena = __div(F(1.0), _sqrt(lena));
            fCos = a.x * lena;
            fSin = a.y * lena;
        }

        a = _cuCaddf(a, _make_cuComplex(len * fCos, len * fSin));

        Real len2 = F(0.5) * (a.x * a.x + a.y * a.y + b.x * b.x + b.y * b.y);
        if (len2 < _CLG_FLT_MIN_)
        {
            a = _make_cuComplex(F(0.0), F(0.0));
            b = _make_cuComplex(F(0.0), F(0.0));
            return;
        }
        len2 = __div(F(1.0), _sqrt(len2));
        a.x = a.x * len2;
        a.y = a.y * len2;
        b.x = b.x * len2;
        b.y = b.y * len2;
    }

    static __device__ __inline__ void FrancisQRStep(CLGComplex* H, INT dm)
    {
        CLGComplex s, t, x, y, z;
        doubleEigen(s, t,
            H[dm * dm - dm - 2],
            H[dm * dm - dm - 1],
            H[dm * dm - 2],
            H[dm * dm - 1]);

        doubleShift(x, y, z, s, t,
            H[0],
            H[1],
            H[dm],
            H[dm + 1],
            H[2 * dm + 1]);

        CLGComplex u[9];
        for (INT k = 0; k <= dm - 3; ++k)
        {
            threeHouseHolder(x, y, z);
            u[0] = _cuCmulf(_cuConjf(x), x);
            u[1] = _cuCmulf(_cuConjf(y), x);
            u[2] = _cuCmulf(_cuConjf(z), x);
            u[3] = _cuCmulf(_cuConjf(x), y);
            u[4] = _cuCmulf(_cuConjf(y), y);
            u[5] = _cuCmulf(_cuConjf(z), y);
            u[6] = _cuCmulf(_cuConjf(x), z);
            u[7] = _cuCmulf(_cuConjf(y), z);
            u[8] = _cuCmulf(_cuConjf(z), z);

            u[0].x = u[0].x - F(1.0);
            u[4].x = u[4].x - F(1.0);
            u[8].x = u[8].x - F(1.0);

            INT q = (k < 1) ? 1 : k;
            for (INT hx = (q - 1); hx < dm; ++hx)
            {
                const CLGComplex newhy1 = _cuCaddf(_cuCmulf(u[0], H[k * dm + hx]), _cuCaddf(_cuCmulf(u[1], H[(k + 1) * dm + hx]), _cuCmulf(u[2], H[(k + 2) * dm + hx])));
                const CLGComplex newhy2 = _cuCaddf(_cuCmulf(u[3], H[k * dm + hx]), _cuCaddf(_cuCmulf(u[4], H[(k + 1) * dm + hx]), _cuCmulf(u[5], H[(k + 2) * dm + hx])));
                const CLGComplex newhy3 = _cuCaddf(_cuCmulf(u[6], H[k * dm + hx]), _cuCaddf(_cuCmulf(u[7], H[(k + 1) * dm + hx]), _cuCmulf(u[8], H[(k + 2) * dm + hx])));

                H[k * dm + hx] = _make_cuComplex(-newhy1.x, -newhy1.y);
                H[(k + 1) * dm + hx] = _make_cuComplex(-newhy2.x, -newhy2.y);
                H[(k + 2) * dm + hx] = _make_cuComplex(-newhy3.x, -newhy3.y);
            }

            q = k + 3;
            if (q >= dm)
            {
                q = dm - 1;
            }
            for (INT hy = 0; hy <= q; ++hy)
            {
                const CLGComplex newhx1 = _cuCaddf(_cuCmulf(u[0], H[hy * dm + k]), _cuCaddf(_cuCmulf(u[3], H[hy * dm + k + 1]), _cuCmulf(u[6], H[hy * dm + k + 2])));
                const CLGComplex newhx2 = _cuCaddf(_cuCmulf(u[1], H[hy * dm + k]), _cuCaddf(_cuCmulf(u[4], H[hy * dm + k + 1]), _cuCmulf(u[7], H[hy * dm + k + 2])));
                const CLGComplex newhx3 = _cuCaddf(_cuCmulf(u[2], H[hy * dm + k]), _cuCaddf(_cuCmulf(u[5], H[hy * dm + k + 1]), _cuCmulf(u[8], H[hy * dm + k + 2])));

                H[hy * dm + k] = _make_cuComplex(-newhx1.x, -newhx1.y);
                H[hy * dm + k + 1] = _make_cuComplex(-newhx2.x, -newhx2.y);
                H[hy * dm + k + 2] = _make_cuComplex(-newhx3.x, -newhx3.y);
            }

            x = H[(k + 1) * dm + k];
            y = H[(k + 2) * dm + k];
            if (k < dm - 3)
            {
                z = H[(k + 3) * dm + k];
            }
        }

        twoHouseHolder(x, y);
        u[0] = _cuCmulf(_cuConjf(x), x);
        u[1] = _cuCmulf(_cuConjf(y), x);
        u[2] = _cuCmulf(_cuConjf(x), y);
        u[3] = _cuCmulf(_cuConjf(y), y);

        u[0].x = u[0].x - F(1.0);
        u[3].x = u[3].x - F(1.0);

        INT startX = (2 == dm) ? 0 : (dm - 3);
        for (INT hx = startX; hx < dm; ++hx)
        {
            const CLGComplex newhy1 = _cuCaddf(_cuCmulf(u[0], H[(dm - 2) * dm + hx]), _cuCmulf(u[1], H[(dm - 1) * dm + hx]));
            const CLGComplex newhy2 = _cuCaddf(_cuCmulf(u[2], H[(dm - 2) * dm + hx]), _cuCmulf(u[3], H[(dm - 1) * dm + hx]));

            H[(dm - 2) * dm + hx] = _make_cuComplex(-newhy1.x, -newhy1.y);
            H[(dm - 1) * dm + hx] = _make_cuComplex(-newhy2.x, -newhy2.y);
        }
        for (INT hy = 0; hy < dm; ++hy)
        {
            const CLGComplex newhx1 = _cuCaddf(_cuCmulf(u[0], H[hy * dm + dm - 2]), _cuCmulf(u[2], H[hy * dm + dm - 1]));
            const CLGComplex newhx2 = _cuCaddf(_cuCmulf(u[1], H[hy * dm + dm - 2]), _cuCmulf(u[3], H[hy * dm + dm - 1]));

            H[hy * dm + dm - 2] = _make_cuComplex(-newhx1.x, -newhx1.y);
            H[hy * dm + dm - 1] = _make_cuComplex(-newhx2.x, -newhx2.y);
        }
    }

    static __device__ __inline__ void CalculateEigen2x2(
        CLGComplex& h00, CLGComplex& h01,
        CLGComplex& h10, CLGComplex& h11)
    {
        const CLGComplex omega = _cuCmulf(h10, h01);
        const Real fOmegaSq = omega.x * omega.x + omega.y * omega.y;
        if (fOmegaSq > _CLG_FLT_MIN_)
        {
            const CLGComplex xi = _make_cuComplex(
                F(0.5) * (h11.x - h00.x),
                F(0.5) * (h11.y - h00.y));
            CLGComplex eta = _cuCaddf(_cuCmulf(xi, xi), omega);
            if (__cuCabsSqf(eta) > _CLG_FLT_MIN_)
            {
                eta = __cuCsqrtf(eta);
            }
            else
            {
                eta = _make_cuComplex(F(0.0), F(0.0));
            }
            const CLGComplex divider1 = _cuCaddf(eta, xi);
            const CLGComplex divider2 = _cuCsubf(eta, xi);
            if (__cuCabsSqf(divider1) > _CLG_FLT_MIN_ && __cuCabsSqf(divider2) > _CLG_FLT_MIN_)
            {
                if (xi.x * eta.x + xi.y * eta.y < F(0.0))
                {
                    h00 = _cuCaddf(h11, _cuCdivf(omega, _cuCaddf(eta, xi)));
                    h11 = _cuCsubf(h11, _cuCdivf(omega, _cuCsubf(eta, xi)));
                }
                else
                {
                    h00 = _cuCsubf(h11, _cuCdivf(omega, _cuCsubf(eta, xi)));
                    h11 = _cuCaddf(h11, _cuCdivf(omega, _cuCaddf(eta, xi)));
                }
            }
        }
        h10 = _make_cuComplex(F(0.0), F(0.0));
    }

    __device__ __inline__ void FrancisQRIteration()
    {
        CLGComplex submatrix[N * N];
        INT decomp[2];
        for (INT i = 0; i < N * N; ++i)
        {
            checkMatrixIndexDoubleShift(m_me, decomp, N);
            const INT dm = decomp[1] - decomp[0];
            if (dm < 2)
            {
                return;
            }

            if (2 == dm)
            {
                CalculateEigen2x2(m_me[decomp[0] * N + decomp[0]],       m_me[decomp[0] * N + decomp[0] + 1],
                                  m_me[(decomp[0] + 1) * N + decomp[0]], m_me[(decomp[0] + 1) * N + decomp[0] + 1]);
            }
            else
            {
                for (INT y = 0; y < dm; ++y)
                {
                    memcpy(submatrix + y * dm, m_me + (decomp[0] + y) * N + decomp[0], sizeof(CLGComplex) * dm);
                }

                FrancisQRStep(submatrix, dm);

                for (INT y = 0; y < dm; ++y)
                {
                    memcpy(m_me + (decomp[0] + y) * N + decomp[0], submatrix + y * dm, sizeof(CLGComplex) * dm);
                }
            }
        }
    }

    #pragma endregion

private:

    __device__ __inline__ void InversePower(const CLGComplex& ev, CLGComplex* startV, Real eps = F(1.0e-10)) const
    {
        INT iMaxIte = N;
        CLGComplex tmp[N];
        const deviceSLNC<N, NoE> m = SubCompC(ev);
        deviceSLNC<N, NoE> q(m);
        const deviceSLNC<N, NoE> r = q.QR();
        q.Dagger();
        CLGComplex mv[N];
        m.MulVectorElement(mv, startV);
        Real s = __cuCabsSqf(mv[0]);
        for (INT i = 1; i < N; ++i)
        {
            s += __cuCabsSqf(mv[i]);
        }
        while (iMaxIte > 0 && s > eps)
        {
            --iMaxIte;
            memcpy(tmp, startV, sizeof(CLGComplex) * N);
            q.MulVectorElement(startV, tmp);
            r.BackwardSubstitutionElement(startV);
            s = __cuCabsSqf(startV[0]);
            for (INT i = 1; i < N; ++i)
            {
                s += __cuCabsSqf(startV[i]);
            }

            if (s > _CLG_FLT_MIN_)
            {
                const Real fFactor = __rcp(_sqrt(s));
                for (INT i = 0; i < N; ++i)
                {
                    startV[i] = cuCmulf_cr(startV[i], fFactor);
                }
            }

            m.MulVectorElement(mv, startV);
            s = __cuCabsSqf(mv[0]);
            for (INT i = 1; i < N; ++i)
            {
                s += __cuCabsSqf(mv[i]);
            }
        }
        if (0 == iMaxIte && s > F(0.000001))
        {
            printf("warning: Inverse Power Max Iteration Reached, last s = %.20f\n", s);
        }
    }

public:

    __device__ __inline__ void EigenValues(CLGComplex* evs) const
    {
        deviceSLNC<N, NoE> tmp = HessenbergC();
        tmp.FrancisQRIteration();
        for (INT i = 0; i < N; ++i)
        {
            evs[i] = tmp.m_me[i * N + i];
        }
    }

    __device__ __inline__ deviceSLNC<N, NoE> EigenSystem(CLGComplex* evs) const
    {
        EigenValues(evs);
        CLGComplex tmpv[N];

        deviceSLNC<N, NoE> eigenVectors;
        for (INT i = 0; i < N; ++i)
        {
            for (INT j = 0; j < N; ++j)
            {
                tmpv[j] = (i == j) ? _onec : _zeroc;
            }

            InversePower(evs[i], tmpv);
            memcpy(eigenVectors.m_me + i * N, tmpv, sizeof(CLGComplex) * N);
        }
        return eigenVectors;
    }

    __device__ __inline__ deviceSLNC<N, NoE> Log() const
    {
        CLGComplex eign[N];
        deviceSLNC<N, NoE> egv = EigenSystem(eign).Transpose();
        deviceSLNC<N, NoE> egvInv = egv.InverseC();

        deviceSLNC<N, NoE> ret = makeSLNCZero();
        for (INT i = 0; i < N; ++i)
        {
            ret.m_me[i * N + i] = __cuClogf(eign[i]);
        }
        ret.Mul(egvInv);
        return egv.MulC(ret);
    }

    __device__ __inline__ deviceSLNC<N, NoE> StrictExp() const
    {
        CLGComplex eign[N];
        deviceSLNC<N, NoE> egv = EigenSystem(eign).Transpose();
        deviceSLNC<N, NoE> egvInv = egv.InverseC();

        deviceSLNC<N, NoE> ret = makeSLNCZero();
        for (INT i = 0; i < N; ++i)
        {
            ret.m_me[i * N + i] = __cuCexpf(eign[i]);
        }
        ret.Mul(egvInv);
        ret = egv.MulC(ret);
        ret.NormalizeDet();
        return ret;
    }

    __device__ __inline__ deviceSLNC<N, NoE> Power(Real fpower) const
    {
        CLGComplex eign[N];
        deviceSLNC<N, NoE> egv = EigenSystem(eign).Transpose();
        deviceSLNC<N, NoE> egvInv = egv.InverseC();

        deviceSLNC<N, NoE> ret = makeSLNCZero();
        for (INT i = 0; i < N; ++i)
        {
            ret.m_me[i * N + i] = __cuCpowerf(eign[i], fpower);
        }
        ret.Mul(egvInv);
        return egv.MulC(ret);
    }

    __device__ __inline__ void Inverse()
    {
        *this = InverseC();
    }

#pragma endregion

    CLGComplex m_me[NoE];
};

#define _TYPEDEFSLNC(n, moe) typedef deviceSLNC<n, moe> deviceSL##n##C;
#define _DEF_F2_TO_SLNC(n, imp) _DEF_F2SLNC_N(n, imp, 256, 256, 256, 256, 256, 128, 128, 128, 64, 64, 64, 32, 16, 16, 4)

_DEF_F2_TO_SLNC(_MAX_SLNC, _TYPEDEFSLNC)

__END_NAMESPACE

#endif //#ifndef _SLNC_H_

//=============================================================================
// END OF FILE
//=============================================================================
