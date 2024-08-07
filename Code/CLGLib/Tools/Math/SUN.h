//=============================================================================
// FILENAME : SUN.h
// 
// DESCRIPTION:
// This is helper functions for calculate D=N complex matrix
// NOTE: It is not nessary a SUN
//
// 
// for SU4
// 0, 1, 2, 3
// 4, 5, 6, 7
// ...
//
// REVISION:
//  [06/30/2024 nbale]
//=============================================================================

#ifndef _SUN_H_
#define _SUN_H_

__BEGIN_NAMESPACE

template<INT N, INT NoE>
struct deviceSUN
{

public:
    __device__ deviceSUN()
    {

    }

    __device__ deviceSUN(const deviceSUN<N, NoE>& other)
    {
        memcpy(m_me, other.m_me, sizeof(CLGComplex) * N * N);
    }

    __device__ deviceSUN(const CLGComplex* other)
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
                printf("%1.7f%s%1.7f I",
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
                printf("%1.7f%s%1.7f I,",
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
    __device__ __inline__ static deviceSUN<N, NoE> makeSUNZero()
    {
        deviceSUN<N, NoE> ret;
        for (INT i = 0; i < N * N; ++i)
        {
            ret.m_me[i] = _make_cuComplex(F(0.0), F(0.0));
        }
        return ret;
    }

    /**
    * ret = I
    */
    __device__ __inline__ static deviceSUN<N, NoE> makeSUNId()
    {
        deviceSUN<N, NoE> ret;
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

    __device__ __inline__ static deviceSUN<N, NoE> makeSUNGenerator(UINT uiGenerator)
    {
        deviceSUN<N, NoE> ret = makeSUNZero();
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
                ret.m_me[m * N + m] = _make_cuComplex(fac * (- m), F(0.0));
            }
        }

        return ret;
    }

    /**
    * Here, we keep same with Bridge++, that H(P)/D.O.F. = 0.5
    * can be called only after CLatticeData is created
    * ret = i r_a T_a, r_a is random real number, T_a are generators
    * r_a T_a = r1 T1 + r2 T2 + ...
    *
    *     r3 +r8/sqrt3     r1-ir2        r4-ir5
    * =   r1+ir2       -r3+r8/sqrt3      r6-ir7
    *     r4+ir5           r6+ir7      -2r8/sqrt3
    *
    */
    __device__ __inline__ static deviceSUN<N, NoE> makeSUNRandomGenerator(UINT fatIndex)
    {
        deviceSUN<N, NoE> ret = makeSUNZero();
        for (UINT i = 0; i < N * N - 1; ++i)
        {
            const Real r = _deviceRandomGaussFSqrt2(fatIndex);
            ret.Add(makeSUNGenerator(i).MulCompC(_make_cuComplex(F(0.0), r)));
        }
        return ret;
    }

    __device__ __inline__ static deviceSUN<N, NoE> makeSUNSumGenerator(Real fDivide)
    {
        deviceSUN<N, NoE> ret = makeSUNZero();
        for (UINT i = 0; i < N * N - 1; ++i)
        {
            ret.Add(makeSUNGenerator(i).MulRealC(fDivide));
        }
        return ret;
    }

    /**
    * can be called only after CLatticeData is created
    * ret = random
    */
    __device__ __inline__ static deviceSUN<N, NoE> makeSUNRandom(UINT fatIndex)
    {
        deviceSUN<N, NoE> ret;
        for (UINT i = 0; i < N * N; ++i)
        {
            ret.m_me[i] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));;
        }
        ret.Norm();
        return ret;
    }

    /**
    * a random matrix which does not has to be SU3
    */
    __device__ __inline__ static deviceSUN<N, NoE> makeSUNRandomAny(UINT fatIndex)
    {
        deviceSUN<N, NoE> ret;
        for (UINT i = 0; i < N * N; ++i)
        {
            ret.m_me[i] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));;
        }
        return ret;
    }

    template<INT NoVE>
    __device__ __inline__ static deviceSUN<N, NoE> makeSUNContractV(const deviceSUNVector<N, NoVE>& left, const deviceSUNVector<N, NoVE>& right)
    {
        deviceSUN<N, NoE> ret;
        for (INT y = 0; y < N; ++y)
        {
            for (INT x = 0; x < N; ++x)
            {
                ret.m_me[y * N + x] = _cuCmulf(_cuConjf(left.m_ve[x]), right.m_ve[y]);
            }
        }
        return ret;
    }


#pragma endregion

#pragma region operators

    __device__ __inline__ void Add(const deviceSUN<N, NoE>& right)
    {
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i] = _cuCaddf(m_me[i], right.m_me[i]);
        }
    }

    __device__ __inline__ void AddDagger(const deviceSUN<N, NoE>& right)
    {
        for (INT y = 0; y < N; ++y)
        {
            for (INT x = 0; x < N; ++x)
            {
                m_me[y * N + x] = _cuCaddf(m_me[y * N + x], _cuConjf(right.m_me[x * N + y]));
            }
        }
    }

    __device__ __inline__ void Sub(const deviceSUN<N, NoE>& right)
    {
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i] = _cuCsubf(m_me[i], right.m_me[i]);
        }
    }

    __device__ __inline__ void SubDagger(const deviceSUN<N, NoE>& right)
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

    __device__ __inline__ void Mul(const deviceSUN<N, NoE>& right)
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

    __device__ __inline__ void DaggerMul(const deviceSUN<N, NoE>& right)
    {
        Dagger();
        Mul(right);
    }

    __device__ __inline__ void MulDagger(const deviceSUN<N, NoE>& right)
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

    __device__ __inline__ void MulOnMe(const deviceSUN<N, NoE>& left)
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

    __device__ __inline__ void MulOnMeDN(const deviceSUN<N, NoE>& left)
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

    __device__ __inline__ void MulOnMeND(const deviceSUN<N, NoE>& left)
    {
        Dagger();
        MulOnMe(left);
    }

    __device__ __inline__ deviceSUN<N, NoE> MulC(const deviceSUN<N, NoE>& right) const
    {
        deviceSUN<N, NoE> ret;
        for (INT x = 0; x < N; ++x)
        {
            for (INT y = 0; y < N; ++y)
            {
                ret.m_me[y * N + x] = _cuCmulf(m_me[y * N], right.m_me[x]);
                for (INT i = 1; i < N; ++i)
                {
                    //ret.m_me[0] = __LINE_MUL(0, 1, 2, 0, 3, 6);
                    //ret.m_me[1] = __LINE_MUL(0, 1, 2, 1, 4, 7);
                    //ret.m_me[3] = __LINE_MUL(3, 4, 5, 0, 3, 6);
                    //for x=0,y=0, n=0, left is 0,1,2  right is 0,3,6
                    //for x=1,y=0, n=1, left is 0,1,2  right is 1,4,7
                    //for x=0,y=1, n=3, left is 3,4,5  right is 0,3,6
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

    __device__ __inline__ deviceSUN<N, NoE> MulDaggerC(const deviceSUN<N, NoE>& right) const
    {
        deviceSUN<N, NoE> ret;
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

    __device__ __inline__ deviceSUN<N, NoE> DaggerMulC(const deviceSUN<N, NoE>& right) const
    {
        deviceSUN<N, NoE> ret;
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

    __device__ __inline__ deviceSUN<N, NoE> AddC(const deviceSUN<N, NoE>& right) const { deviceSUN<N, NoE> ret(*this); ret.Add(right); return ret; }
    __device__ __inline__ deviceSUN<N, NoE> AddCompC(const CLGComplex& right) const { deviceSUN<N, NoE> ret(*this); ret.AddComp(right); return ret; }
    __device__ __inline__ deviceSUN<N, NoE> AddRealC(const Real& right) const { deviceSUN<N, NoE> ret(*this); ret.AddReal(right); return ret; }

    __device__ __inline__ deviceSUN<N, NoE> SubC(const deviceSUN<N, NoE>& right) const { deviceSUN<N, NoE> ret(*this); ret.Sub(right); return ret; }
    __device__ __inline__ deviceSUN<N, NoE> SubCompC(const CLGComplex& right) const { deviceSUN<N, NoE> ret(*this); ret.SubComp(right); return ret; }
    __device__ __inline__ deviceSUN<N, NoE> SubRealC(const Real& right) const { deviceSUN<N, NoE> ret(*this); ret.SubReal(right); return ret; }

    __device__ __inline__ deviceSUN<N, NoE> MulCompC(const CLGComplex& right) const { deviceSUN<N, NoE> ret(*this); ret.MulComp(right); return ret; }
    __device__ __inline__ deviceSUN<N, NoE> MulRealC(const Real& right) const { deviceSUN<N, NoE> ret(*this); ret.MulReal(right); return ret; }
    __device__ __inline__ deviceSUN<N, NoE> DivCompC(const CLGComplex& right) const { deviceSUN<N, NoE> ret(*this); ret.DivComp(right); return ret; }

    /**
    * self block mult must be square
    */
    __device__ __inline__ void BlockMul(
        const deviceSUN<N, NoE>& right,
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

    __device__ __inline__ void BlockCopy(const deviceSUN<N, NoE>& right,
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

    __device__ __inline__ void MulVectorElement(CLGComplex* res, const CLGComplex * v) const
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

#if 0

    //Avoid use of recursion on device, use LU decompose to compute determinet instead.

protected:

    /**
    * assume len(includeX) = len(includeY) = n
    */
    __device__ __inline__ CLGComplex SubDeterminent(UINT n, UINT* includeX, UINT* includeY) const
    {
        if (1 == n)
        {
            return m_me[includeY[0] * N + includeX[0]];
        }

        if (2 == n)
        {
            return _cuCaddf(_cuCmulf(m_me[includeY[0] * N + includeX[0]], m_me[includeY[1] * N + includeX[1]]), _cuCmulf(m_me[includeY[1] * N + includeX[0]], m_me[includeY[0] * N + includeX[1]]));
        }

        if (3 == n)
        {
            return _cuCsubf(
                _cuCaddf(
                    _cuCaddf(
                        _cuCmulf(_cuCmulf(m_me[includeY[0] * N + includeX[0]], m_me[includeY[1] * N + includeX[1]]), m_me[includeY[2] * N + includeX[2]]),
                        _cuCmulf(_cuCmulf(m_me[includeY[0] * N + includeX[1]], m_me[includeY[1] * N + includeX[2]]), m_me[includeY[2] * N + includeX[0]])
                    ),
                    _cuCmulf(_cuCmulf(m_me[includeY[0] * N + includeX[2]], m_me[includeY[1] * N + includeX[0]]), m_me[includeY[2] * N + includeX[1]])
                ),
                _cuCaddf(
                    _cuCaddf(
                        _cuCmulf(_cuCmulf(m_me[includeY[0] * N + includeX[0]], m_me[includeY[1] * N + includeX[2]]), m_me[includeY[2] * N + includeX[1]]),
                        _cuCmulf(_cuCmulf(m_me[includeY[0] * N + includeX[2]], m_me[includeY[1] * N + includeX[1]]), m_me[includeY[2] * N + includeX[0]])
                    ),
                    _cuCmulf(_cuCmulf(m_me[includeY[0] * N + includeX[1]], m_me[includeY[1] * N + includeX[0]]), m_me[includeY[2] * N + includeX[2]])
                )
            );
        }

        UINT tmpX[N];
        CLGComplex res = _zeroc;
        for (INT i = 0; i < n; ++i)
        {
            INT k = 0;
            for (INT j = 0; j < n; ++j)
            {
                if (i != j)
                {
                    tmpX[k] = includeX[j];
                    ++k;
                }
            }

            if ((i + n) & 1)
            {
                res = _cuCsubf(res, _cuCmulf(m_me[includeY[n - 1] * N + includeX[i]], SubDeterminent(n - 1, tmpX, includeY)));
            }
            else
            {
                res = _cuCaddf(res, _cuCmulf(m_me[includeY[n - 1] * N + includeX[i]], SubDeterminent(n - 1, tmpX, includeY)));
            }
        }
        return res;
    }

    __device__ __inline__ CLGComplex Determinent() const
    {
        UINT x[N];
        UINT y[N];
        for (INT i = 0; i < N; ++i)
        {
            x[i] = i;
            y[i] = i;
        }

        return SubDeterminent(N, x, y);
    }

#endif

public:

    __device__ __inline__ CLGComplex Determinent() const
    {
        deviceSUN<N, NoE> tmp(*this);
        tmp.LUNoReturn();
        for (INT i = 1; i < N; ++i)
        {
            tmp.m_me[0] = _cuCmulf(tmp.m_me[0], tmp.m_me[i * N + i]);
        }
        return tmp.m_me[0];
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

    __device__ __inline__ deviceSUN<N, NoE> ReC() const { deviceSUN<N, NoE> ret(*this); ret.Re(); return ret; }

    __device__ __inline__ void Im()
    {
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i].x = m_me[i].y;
            m_me[i].y = F(0.0);
        }
    }

    __device__ __inline__ deviceSUN<N, NoE> ImC() const { deviceSUN<N, NoE> ret(*this); ret.Im(); return ret; }

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

    __device__ __inline__ deviceSUN<N, NoE> DaggerC() const
    {
        deviceSUN<N, NoE> ret;
        for (INT x = 0; x < N; ++x)
        {
            for (INT y = 0; y < N; ++y)
            {
                ret.m_me[y * N + x] = _cuConjf(ret.m_me[x * N + y]);
            }
        }
        return ret;
    }


    __device__ __inline deviceSUN<N, NoE> Transpose() const
    {
        deviceSUN<N, NoE> ret;
        for (INT x = 0; x < N; ++x)
        {
            for (INT y = 0; y < N; ++y)
            {
                ret.m_me[y * N + x] = ret.m_me[x * N + y];
            }
        }
        return ret;
    }

    /**
    * A=(U-U^dagger)/2
    * res = A - tr(A)/3
    */
    __device__ __inline__ void Ta()
    {
        Real trace = F(0.0);
        for (INT y = 0; y < N; ++y)
        {
            for (INT x = y; x < N; ++x)
            {
                if (x == y)
                {
                    m_me[y * N + x].x = 0;
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

    /**
    * Return Tr[Im[a].Im[b]]
    * It is used in topological charge
    */
    ////__device__ __inline__ static Real TrIm(const deviceSU3& a, const deviceSU3& b)
    ////{
    ////    //0,1,2 * 0,3,6
    ////    Real ret = a.m_me[0].y * b.m_me[0].y + a.m_me[1].y * b.m_me[3].y + a.m_me[2].y * b.m_me[6].y;
    ////    //3,4,5 * 1,4,7
    ////    ret += a.m_me[3].y * b.m_me[1].y + a.m_me[4].y * b.m_me[4].y + a.m_me[5].y * b.m_me[7].y;
    ////    //6,7,8 * 2,5,8
    ////    ret += a.m_me[6].y * b.m_me[2].y + a.m_me[7].y * b.m_me[5].y + a.m_me[8].y * b.m_me[8].y;

    ////    return ret;
    ////}

    /**
    * res = U - U^dagger
    * It is like a matrix Im(M)
    * res = 2i Im(M), so is called iIm2
    */
    //__device__ __inline__ void iIm2()
    //{
    //    //0 1 2
    //    //3 4 5
    //    //6 7 8

    //    //new [1] = [1] - conj([3])
    //    //new [3] = [3] - conj([1]) = -conj(new [1])
    //    const CLGComplex new1 = _cuCsubf(m_me[1], _cuConjf(m_me[3]));
    //    const CLGComplex new2 = _cuCsubf(m_me[2], _cuConjf(m_me[6]));
    //    const CLGComplex new5 = _cuCsubf(m_me[5], _cuConjf(m_me[7]));
    //    m_me[1] = _make_cuComplex(_cuCrealf(new1), _cuCimagf(new1));
    //    m_me[3] = _make_cuComplex(-_cuCrealf(m_me[1]), _cuCimagf(m_me[1]));
    //    m_me[2] = _make_cuComplex(_cuCrealf(new2), _cuCimagf(new2));
    //    m_me[6] = _make_cuComplex(-_cuCrealf(m_me[2]), _cuCimagf(m_me[2]));
    //    m_me[5] = _make_cuComplex(_cuCrealf(new5), _cuCimagf(new5));
    //    m_me[7] = _make_cuComplex(-_cuCrealf(m_me[5]), _cuCimagf(m_me[5]));

    //    m_me[0] = _make_cuComplex(F(0.0), F(2.0) * m_me[0].y);
    //    m_me[4] = _make_cuComplex(F(0.0), F(2.0) * m_me[4].y);
    //    m_me[8] = _make_cuComplex(F(0.0), F(2.0) * m_me[8].y);
    //}

    /**
    * return -i(U-U^dagger) = ((-iU)+(-iU)dagger)
    */
    //__device__ __inline__ deviceSU3 Im2C() const
    //{
    //    deviceSU3 ret;
    //    //0 1 2
    //    //3 4 5
    //    //6 7 8
    //    // -i(x-y*)=i(y*-x)=i(yr-iyi-xr-ixi)=(yi+xi)+i(yr-xr)
    //    // new 1 = -i(1-3*)
    //    // new 2 = -i(2-6*)
    //    // new 5 = -i(5-7*)
    //    // (new 1)*= i(1*-3) = -i(3-1*) = new 3
    //    // new 0 = -i(0-0*) = 2Im(0)

    //    ret.m_me[1] = _make_cuComplex(m_me[1].y + m_me[3].y, m_me[3].x - m_me[1].x);
    //    ret.m_me[2] = _make_cuComplex(m_me[2].y + m_me[6].y, m_me[6].x - m_me[2].x);
    //    ret.m_me[5] = _make_cuComplex(m_me[5].y + m_me[7].y, m_me[7].x - m_me[5].x);
    //    ret.m_me[3] = _cuConjf(ret.m_me[1]);
    //    ret.m_me[6] = _cuConjf(ret.m_me[2]);
    //    ret.m_me[7] = _cuConjf(ret.m_me[5]);

    //    ret.m_me[0] = _make_cuComplex(F(2.0) * m_me[0].y, F(0.0));
    //    ret.m_me[4] = _make_cuComplex(F(2.0) * m_me[4].y, F(0.0));
    //    ret.m_me[8] = _make_cuComplex(F(2.0) * m_me[8].y, F(0.0));
    //    return ret;
    //}

    /**
    * make any matrix to SUN
    * 
    * This is from bridge++ but there is a problem
    * 
    * let M = {{0.9879009 + 0.0051302 I, -0.0682856 - 0.0641515 I, -0.0821236 - 0.0152904 I, -0.0859499 - 0.0296949 I}, 
               {0.0625330 - 0.0694527 I, 0.9924969 - 0.0488579 I, -0.0542777 - 0.0216560 I, 0.0198289 - 0.0045309 I}, 
               {0.0796094 - 0.0245828 I, 0.0462495 - 0.0277931 I, 0.9932243 - 0.0171367 I, -0.0316254 - 0.0485630 I}, 
               {0.0872017 - 0.0273023 I, -0.0270464 - 0.0118957 I, 0.0270664 - 0.0456755 I, 0.9920595 + 0.0614532 I}}

      Q.R = M
      with
      Q = {{-0.9879142 + 0.0000000 I, 0.0650371 + 0.0674427 I, 0.0818422 + 0.0167314 I, -0.0876184 - 0.0243346 I}, 
           {-0.0621715 + 0.0697765 I, -0.9936988 - 0.0001736 I, 0.0538887 + 0.0226065 I, 0.0195115 - 0.0057458 I}, 
           {-0.0794807 + 0.0249959 I, -0.0475645 + 0.0254772 I, -0.9933721 - 0.0003223 I, -0.0345618 - 0.0465189 I}, 
           {-0.0870588 + 0.0277547 I, 0.0264265 + 0.0132157 I, -0.0278649 + 0.0451927 I, 0.9939610 + 0.0001197 I}}

      but, Q/Exp[-I Arg[Det[Q]]/4], although was an SUN matrix, but not M (which is already an SUN matrix)
    */
    __device__ __inline__ void Norm()
    {
        //QR();
        //CLGComplex d = Determinent();
        //const Real fArg = -__cuCargf(d) / N;
        //MulComp(_make_cuComplex(_cos(fArg), _sin(fArg)));
        Proj();
    }

    /**
     * HYP projection
     * When Det[U] is large, it need large ite...
     */
    __device__ __inline__ void Proj(BYTE ite = N + 1)
    {
        deviceSUN<N, NoE> tmp;

        for (BYTE byIt = 0; byIt < ite; ++byIt)
        {
            memcpy(tmp.m_me, m_me, sizeof(CLGComplex) * N * N);
            tmp.DaggerMul(*this);

            const Real fDiv = __rcp(tmp.ReTr() / N);
            tmp.MulReal(fDiv);
            MulReal(_sqrt(fDiv));

            //x = (-1/2) me^+ me
            tmp.MulReal(F(-0.5));

            //x += (3/2)
            tmp.AddReal(F(1.5));

            //me = me.x
            Mul(tmp);

            //coef = det(me)
            // If the arg is more than pi/2, this will fall into -1
            //MulComp(_make_cuComplex(F(1.0), -Determinent().y / N));
            const Real fArg = -__cuCargf(Determinent()) / N;
            MulComp(_make_cuComplex(_cos(fArg), _sin(fArg)));
        }
    }

    /**
     * Cabbibo-Marinari Projection
     * Only one iteration can make it SU3, however, 6 iteration makes it large trace
     * NOTE: It has been tested using Mathematica, the trace is not increasing...
     *
     * Note for iteration:
     * u1 = m.cmproj()
     * m1 = u1.m
     * u2 = m1.cmproj().u1
     * m3 = u2.m
     * u3 = m3.cmproj().u2
     * ...
     * retr(u3.m)>retr(u2.m)>retr(u1.m), where u1,u2,u3 are su3, m is random matrix
     * finally, u = dagger(u3) so that retr(u3^+ . m) or retr(u3.m^+) is maximized
     */
    //__device__ __inline__ void CabbiboMarinariProj(/*BYTE ite = 1*/)
    //{
    //    //for (BYTE iteration = 0; iteration < ite; ++iteration)
    //    //{
    //    CLGComplex a11 = _cuCaddf(_cuConjf(m_me[0]), m_me[4]);
    //    CLGComplex b11 = _cuCaddf(_cuConjf(m_me[0]), m_me[8]);
    //    CLGComplex c22 = _cuCaddf(_cuConjf(m_me[4]), m_me[8]);
    //    CLGComplex a12 = _cuCsubf(_cuConjf(m_me[3]), m_me[1]);
    //    CLGComplex b13 = _cuCsubf(_cuConjf(m_me[6]), m_me[2]);
    //    CLGComplex c23 = _cuCsubf(_cuConjf(m_me[7]), m_me[5]);
    //    //CLGComplex a12 = _cuCsubf(_cuConjf(m_me[1]), m_me[3]);
    //    //CLGComplex b13 = _cuCsubf(_cuConjf(m_me[2]), m_me[6]);
    //    //CLGComplex c23 = _cuCsubf(_cuConjf(m_me[5]), m_me[7]);
    //    Real fNorm = __rcp(_sqrt(__cuCabsSqf(a11) + __cuCabsSqf(a12)));
    //    a11 = cuCmulf_cr(a11, fNorm);
    //    a12 = cuCmulf_cr(a12, fNorm);
    //    fNorm = __rcp(_sqrt(__cuCabsSqf(b11) + __cuCabsSqf(b13)));
    //    b11 = cuCmulf_cr(b11, fNorm);
    //    b13 = cuCmulf_cr(b13, fNorm);
    //    fNorm = __rcp(_sqrt(__cuCabsSqf(c22) + __cuCabsSqf(c23)));
    //    c22 = cuCmulf_cr(c22, fNorm);
    //    c23 = cuCmulf_cr(c23, fNorm);

    //    /**
    //     * ({
    //    {a11 b11,
    //    a12 c22 - a11 b13 Conjugate[c23],
    //    a11 b13 Conjugate[c22] + a12 c23},
    //    {-b11 Conjugate[a12],
    //    c22 Conjugate[a11] + b13 Conjugate[a12] Conjugate[c23],
    //    c23 Conjugate[a11] - b13 Conjugate[a12] Conjugate[c22]},
    //    {-Conjugate[b13],
    //    -Conjugate[b11] Conjugate[c23],
    //    Conjugate[b11] Conjugate[c22]}
    //    })
    //     */
    //    m_me[0] = _cuCmulf(a11, b11);
    //    m_me[1] = _cuCsubf(_cuCmulf(a12, c22), _cuCmulf(_cuCmulf(a11, b13), _cuConjf(c23)));
    //    m_me[2] = _cuCaddf(_cuCmulf(_cuCmulf(a11, b13), _cuConjf(c22)), _cuCmulf(a12, c23));

    //    m_me[3].x = -b11.x * a12.x - b11.y * a12.y;
    //    m_me[3].y = a12.y * b11.x - a12.x * b11.y;
    //    m_me[4] = _cuCaddf(_cuCmulf(c22, _cuConjf(a11)), _cuCmulf(_cuCmulf(b13, _cuConjf(a12)), _cuConjf(c23)));
    //    m_me[5] = _cuCsubf(_cuCmulf(c23, _cuConjf(a11)), _cuCmulf(_cuCmulf(b13, _cuConjf(a12)), _cuConjf(c22)));

    //    m_me[6].x = -b13.x; m_me[6].y = b13.y;
    //    m_me[7].x = b11.y * c23.y - b11.x * c23.x;
    //    m_me[7].y = b11.x * c23.y + b11.y * c23.x;
    //    m_me[8].x = b11.x * c22.x - b11.y * c22.y;
    //    m_me[8].y = -b11.x * c22.y - b11.y * c22.x;
    //    //}
    //}

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
    __device__ __inline__ deviceSUN<N, NoE> Exp(const CLGComplex& a, BYTE uiPrecision = N + 1) const
    {
        deviceSUN<N, NoE> tmp;

        /**
        * tmp = U
        * loop1: tmp = 1+aU/N
        * loop2: tmp = 1+(aU/(N-1))(1+aU/N)
        * loop3: tmp = 1+(aU/(N-2))(1+aU/(N-1)(1+aU/N))
        * ...
        * loopN: tmp = 1+aU (1 + aU/2 (1+...))
        */
        for (BYTE i = 0; i < uiPrecision; ++i)
        {
            const Real exp_factor = __rcp(uiPrecision - i);
            CLGComplex alpha = cuCmulf_cr(a, exp_factor);
            //aU/(N-i) = this x alpha
            deviceSUN<N, NoE> aUoN(*this); 
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
 
        return tmp;
    }

    __device__ __inline__ deviceSUN<N, NoE> ExpReal(Real a, BYTE uiPrecision = N + 1) const
    {
        deviceSUN<N, NoE> tmp;
        for (BYTE i = 0; i < uiPrecision; ++i)
        {
            deviceSUN<N, NoE> aUoN(*this); 
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
        tmp.Norm();
        return tmp;
    }

private:

    /**
    * when i = 0
    * R=1, A=A
    * 
    * left me be:
    * p
    * 
    * And
    * R be
    * X x x x
    * 0 x x x
    * 0 x x x
    * 0 x x x
    * so that Q.A = R
    * 
    * when i = 1
    * 
    * left me be:
    * me times
    * 1 0 0 0
    * 0 
    * 0   P
    * 0
    * 
    *
    * And A be
    * x x x x
    * 0 x x x
    * 0 x x x
    * 0 x x x
    * 
    * R be
    * x x x x
    * 0 x x x
    * 0 0 x x
    * 0 0 x x
    * so that Q.A = R
    * 
    */
    __device__ __inline__ void OneStepHouseHolderQR(deviceSUN<N, NoE>& A, UINT i)
    {

        //Calculate beta and v
        Real fLength = F(0.0);
        for (INT y = i; y < N; ++y)
        {
            fLength += __cuCabsSqf(A.m_me[y * N + i]);
        }

        if (fLength > _CLG_FLT_MIN_)
        {
            //Calculate P = 1 - beta (v^+ v)
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

                //Me = Me.P
                if (0 == i)
                {
                    //P = 1 - beta (v^+ v)
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
                                m_me[py * N + px].x = - fBeta * m_me[py * N + px].x;
                            }
                            m_me[py * N + px].y = - fBeta * m_me[py * N + px].y;
                        }
                    }

                    //Update A
                    A.MulOnMe(*this);
                }
                else
                {
                    //P = 1 - beta (v^+ v)
                    deviceSUN<N, NoE> tmp = makeSUNId();
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

                    //Update A
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
    __device__ __inline__ deviceSUN<N, NoE> QR()
    {
        deviceSUN<N, NoE> r(*this);
        for (UINT i = 0; i < (N - 1); ++i)
        {
            OneStepHouseHolderQR(r, i);
        }
        return r;
    }

private:

    __device__ __inline__ void OneStepGaussianLU(deviceSUN<N, NoE>& L, UINT i)
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
    __device__ __inline__ deviceSUN<N, NoE> LU()
    {
        deviceSUN<N, NoE> L = makeSUNId();
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

    __device__ __inline__ void OneStepHouseHolderHessenberg(deviceSUN<N, NoE>& A, UINT i)
    {
        //Calculate beta and v
        Real fLength = F(0.0);
        for (INT y = i + 1; y < N; ++y)
        {
            fLength += __cuCabsSqf(A.m_me[y * N + i]);
        }

        if (fLength > _CLG_FLT_MIN_)
        {
            //Calculate P = 1 - beta (v^+ v)
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

                //Me = Me.P
                if (0 == i)
                {
                    //P = 1 - beta (v^+ v)
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

                    //Update A
                    A.MulOnMe(*this);
                    A.Mul(*this);
                }
                else
                {
                    //P = 1 - beta (v^+ v)
                    deviceSUN<N, NoE> tmp = makeSUNId();
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

                    //Update A
                    A.MulOnMe(tmp);
                    A.Mul(tmp);
                    Mul(tmp);
                }
            }
        }
    }

    __device__ __inline__ void OneStepHouseHolderHessenbergC(deviceSUN<N, NoE>& A, UINT i) const
    {
        //Calculate beta and v
        Real fLength = F(0.0);
        for (INT y = i + 1; y < N; ++y)
        {
            fLength += __cuCabsSqf(A.m_me[y * N + i]);
        }

        if (fLength > _CLG_FLT_MIN_)
        {
            //Calculate P = 1 - beta (v^+ v)
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
                //P = 1 - beta (v^+ v)
                deviceSUN<N, NoE> tmp = makeSUNId();
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

                //Update A
                A.MulOnMe(tmp);
                A.Mul(tmp);
            }
        }
    }

public:

    __device__ __inline__ deviceSUN<N, NoE> Hessenberg()
    {
        deviceSUN<N, NoE> r(*this);
        for (UINT i = 0; i < (N - 2); ++i)
        {
            OneStepHouseHolderHessenberg(r, i);
        }
        return r;
    }

    __device__ __inline__ deviceSUN<N, NoE> HessenbergC() const
    {
        deviceSUN<N, NoE> r(*this);
        for (UINT i = 0; i < (N - 2); ++i)
        {
            OneStepHouseHolderHessenbergC(r, i);
        }
        return r;
    }

    #pragma region QR Iteration

    /**
    * for every 2x2 matrix of mtr which is dx x dx
    * from the right buttom corner to the left top
    * Denote it as
    * i  : m11 m12
    * i+1: m21 m22
    * 
    * if m21 ~ 0 (|m21|^2 << (|m11|^2 + |m22|^2))
    * set m21 = 0
    * 
    * decomp[0] is the start index of none-zero Henssenberg matrix
    * decomp[1] is the end   index of none-zero Henssenberg matrix
    * 
    * if (decomp[1] == i + 2)
    * it means, at each step of for, decomp[1] was updated (such that the right bottom of it are always 0)
    * 
    * if (i + 1 > decomp[0] && i + 1 < decomp[1])
    *    i + 1 > decomp[0] means decomp[0] was never updated
    *    i + 1 < decomp[1] means there was at least a two dimension matrix
    * 
    */
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


    /**
    * H = h00 h01
    *     h10 h11
    * a1, a2 be the eigenvalues
    * s = a1 + a2
    * t = a1.a2
    */
    static __device__ __inline__ void doubleEigen(CLGComplex& s, CLGComplex& t,
        const CLGComplex& h00, const CLGComplex& h01,
        const CLGComplex& h10, const CLGComplex& h11)
    {
        s = _cuCaddf(h11, h00);
        t = _cuCsubf(_cuCmulf(h00, h11), _cuCmulf(h10, h01));
    }

    /**
    * H = h00 h01
    *     h10 h11
    * a1, a2 be the eigenvalues
    * s = a1 + a2
    * t = a1.a2
    */
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


    /**
    * calculate v=[a, b, c] for house holder to transform input [a, b, c] to [x, 0, 0]
    */
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

    /**
    * calculate v=[a, b] for house holder to transform input [a, b] to [x, 0]
    */
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

    /**
    * The hessenberg is decomp[0] - decomp[1] square submatrix of dx x dx matrix mtr
    * It does not support when dm=2
    */
    static __device__ __inline__ void FrancisQRStep(CLGComplex* H, INT dm)
    {
        //========= start, calculate s, t, x, y, z
        CLGComplex s, t, x, y, z;
        doubleEigen(s, t,
            H[dm * dm - dm - 2/*(dm - 2) * dm + dm - 2*/],
            H[dm * dm - dm - 1/*(dm - 2) * dm + dm - 1*/],
            H[dm * dm - 2/*(dm - 1) * dm + dm - 2*/],
            H[dm * dm - 1/*(dm - 1) * dm + dm - 1*/]);

        doubleShift(x, y, z, s, t,
            H[0],
            H[1],
            H[dm],
            H[dm + 1],
            H[2 * dm + 1]);

        //3x3 matrix um = 1 - v^T v
        CLGComplex u[9]; 
        for (INT k = 0; k <= dm - 3; ++k)
        {
            //====== calculate um =====
            threeHouseHolder(x, y, z);
            //u [i * 3 + j] = h[j]* h[i]
            u[0] = _cuCmulf(_cuConjf(x), x);
            u[1] = _cuCmulf(_cuConjf(y), x);
            u[2] = _cuCmulf(_cuConjf(z), x);
            u[3] = _cuCmulf(_cuConjf(x), y);
            u[4] = _cuCmulf(_cuConjf(y), y);
            u[5] = _cuCmulf(_cuConjf(z), y);
            u[6] = _cuCmulf(_cuConjf(x), z);
            u[7] = _cuCmulf(_cuConjf(y), z);
            u[8] = _cuCmulf(_cuConjf(z), z);

            //u[i * 3 + i] -= 1
            u[0].x = u[0].x - F(1.0);
            u[4].x = u[4].x - F(1.0);
            u[8].x = u[8].x - F(1.0);

            // H(y=k:k+2; x=q-1:n-1) = um.H(y=k:k+2; x=q-1:n-1)
            INT q = (k < 1) ? 1 : k;
            for (INT hx = (q - 1); hx < dm; ++hx)
            {
                //um 0,1,2, h  k*dm+hx, (k+1)*dm+hx, (k+2)*dm+hx
                const CLGComplex newhy1 = _cuCaddf(_cuCmulf(u[0], H[k * dm + hx]), _cuCaddf(_cuCmulf(u[1], H[(k + 1) * dm + hx]), _cuCmulf(u[2], H[(k + 2) * dm + hx])));
                //um 3,4,5
                const CLGComplex newhy2 = _cuCaddf(_cuCmulf(u[3], H[k * dm + hx]), _cuCaddf(_cuCmulf(u[4], H[(k + 1) * dm + hx]), _cuCmulf(u[5], H[(k + 2) * dm + hx])));
                //um 6,7,8
                const CLGComplex newhy3 = _cuCaddf(_cuCmulf(u[6], H[k * dm + hx]), _cuCaddf(_cuCmulf(u[7], H[(k + 1) * dm + hx]), _cuCmulf(u[8], H[(k + 2) * dm + hx])));

                H[k * dm + hx] = _make_cuComplex(-newhy1.x, -newhy1.y);
                H[(k + 1) * dm + hx] = _make_cuComplex(-newhy2.x, -newhy2.y);
                H[(k + 2) * dm + hx] = _make_cuComplex(-newhy3.x, -newhy3.y);
            }

            // H[y=0:q; x=k:k+2] = H[y=0:q; x=k:k+2] u
            q = k + 3;
            if (q >= dm)
            {
                q = dm - 1;
            }
            for (INT hy = 0; hy <= q; ++hy)
            {
                //um 0,1,2, h  k*dm+hx, (k+1)*dm+hx, (k+2)*dm+hx
                const CLGComplex newhx1 = _cuCaddf(_cuCmulf(u[0], H[hy * dm + k]), _cuCaddf(_cuCmulf(u[3], H[hy * dm + k + 1]), _cuCmulf(u[6], H[hy * dm + k + 2])));
                //um 3,4,5
                const CLGComplex newhx2 = _cuCaddf(_cuCmulf(u[1], H[hy * dm + k]), _cuCaddf(_cuCmulf(u[4], H[hy * dm + k + 1]), _cuCmulf(u[7], H[hy * dm + k + 2])));
                //um 6,7,8
                const CLGComplex newhx3 = _cuCaddf(_cuCmulf(u[2], H[hy * dm + k]), _cuCaddf(_cuCmulf(u[5], H[hy * dm + k + 1]), _cuCmulf(u[8], H[hy * dm + k + 2])));

                H[hy * dm + k] = _make_cuComplex(-newhx1.x, -newhx1.y);
                H[hy * dm + k + 1] = _make_cuComplex(-newhx2.x, -newhx2.y);
                H[hy * dm + k + 2] = _make_cuComplex(-newhx3.x, -newhx3.y);
            }

            //update x,y,z
            x = H[(k + 1) * dm + k];
            y = H[(k + 2) * dm + k];
            if (k < dm - 3)
            {
                z = H[(k + 3) * dm + k];
            }
        }

        //Final step
        twoHouseHolder(x, y);
        u[0] = _cuCmulf(_cuConjf(x), x);
        u[1] = _cuCmulf(_cuConjf(y), x);
        u[2] = _cuCmulf(_cuConjf(x), y);
        u[3] = _cuCmulf(_cuConjf(y), y);

        //u[i * 2 + i] -= 1
        u[0].x = u[0].x - F(1.0);
        u[3].x = u[3].x - F(1.0);

        //H[y=n-2:n-1, x=n-3:n-1] = u.H[y=n-2:n-1, x=n-3:n-1]
        INT startX = (2 == dm) ? 0 : (dm - 3);
        for (INT hx = startX; hx < dm; ++hx)
        {
            const CLGComplex newhy1 = _cuCaddf(_cuCmulf(u[0], H[(dm - 2) * dm + hx]), _cuCmulf(u[1], H[(dm - 1) * dm + hx]));
            const CLGComplex newhy2 = _cuCaddf(_cuCmulf(u[2], H[(dm - 2) * dm + hx]), _cuCmulf(u[3], H[(dm - 1) * dm + hx]));

            H[(dm - 2) * dm + hx] = _make_cuComplex(-newhy1.x, -newhy1.y);
            H[(dm - 1) * dm + hx] = _make_cuComplex(-newhy2.x, -newhy2.y);
        }
        //H[y=0:n-1, x=n-2:n-1] = H[y=0:n-1, x=n-2:n-1].u
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
        //bc
        const CLGComplex omega = _cuCmulf(h10, h01);
        const Real fOmegaSq = omega.x * omega.x + omega.y * omega.y;
        if (fOmegaSq > _CLG_FLT_MIN_)
        {
            //(d-a)/2
            const CLGComplex xi = _make_cuComplex(
                F(0.5) * (h11.x - h00.x),
                F(0.5) * (h11.y - h00.y));
            //sqrt(((d-a)/2)^2 + bc)
            //when bc > fCrit, eta=sqrt(xi^2+bc) cannot be close to xi
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

    /**
    * Assume me is Hessenberg matrix
    * transform me to up traingular with diagonal the eigen values
    * too complicated to preserve all transformation matrices
    */
    __device__ __inline__ void FrancisQRIteration()
    {
        CLGComplex submatrix[N * N];
        INT decomp[2];
        //iteration for at most N times.
        for (INT i = 0; i < N * N; ++i)
        {
            //find hesenberg
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
                //block copy to submatrix
                for (INT y = 0; y < dm; ++y)
                {
                    memcpy(submatrix + y * dm, m_me + (decomp[0] + y) * N + decomp[0], sizeof(CLGComplex) * dm);
                }

                FrancisQRStep(submatrix, dm);

                //copy back
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
    //__device__ __inline__ void InversePower(const CLGComplex& ev, CLGComplex* startV) const
    {
        //for (INT i = 0; i < N; ++i)
        //{
        //    printf("at first startV = %f, %f\n", startV[i].x, startV[i].y);
        //}
        INT iMaxIte = N;
        CLGComplex tmp[N];
        const deviceSUN<N, NoE> m = SubCompC(ev);
        //m.DebugPrint("m");
        deviceSUN<N, NoE> q(m);
        const deviceSUN<N, NoE> r = q.QR();
        q.Dagger();
        CLGComplex mv[N];
        m.MulVectorElement(mv, startV);
        //q.DebugPrint("q");
        //r.DebugPrint("r");
        Real s = __cuCabsSqf(mv[0]);
        for (INT i = 1; i < N; ++i)
        {
            s += __cuCabsSqf(mv[i]);
        }
        //for (INT i = 0; i < N; ++i)
        //{
        //    printf("at first mv = %f, %f\n", mv[i].x, mv[i].y);
        //}
        while (iMaxIte > 0 && s > eps)
        {
            --iMaxIte;
            memcpy(tmp, startV, sizeof(CLGComplex) * N);
            //if (2 * N - 1 == iMaxIte)
            //{
            //    for (INT i = 0; i < N; ++i)
            //    {
            //        printf("at first iteration tmp = %f, %f\n", tmp[i].x, tmp[i].y);
            //        printf("at first iteration startV = %f, %f\n", startV[i].x, startV[i].y);
            //    }
            //}
            q.MulVectorElement(startV, tmp);
            //if (2 * N - 1 == iMaxIte)
            //{
            //    q.DebugPrint("q");
            //    for (INT i = 0; i < N; ++i)
            //    {
            //        printf("at first iteration q.startV = %f, %f\n", startV[i].x, startV[i].y);
            //    }
            //}
            r.BackwardSubstitutionElement(startV);
            //if (2 * N - 1 == iMaxIte)
            //{
            //    r.DebugPrint("r");
            //    for (INT i = 0; i < N; ++i)
            //    {
            //        printf("at first iteration r^{-1}.startV = %f, %f\n", startV[i].x, startV[i].y);
            //    }
            //}
            s = __cuCabsSqf(startV[0]);
            for (INT i = 1; i < N; ++i)
            {
                s += __cuCabsSqf(startV[i]);
                //printf("res startv = %f, %f\n", startV[i].x, startV[i].y);
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
            //if (2 * N - 1 == iMaxIte)
            //{
            //    m.DebugPrint("m");
            //    for (INT i = 0; i < N; ++i)
            //    {
            //        printf("at first iteration res mv = %f, %f\n", mv[i].x, mv[i].y);
            //        printf("at first iteration res startV = %f, %f\n", startV[i].x, startV[i].y);
            //    }
            //}
            s = __cuCabsSqf(mv[0]);
            for (INT i = 1; i < N; ++i)
            {
                s += __cuCabsSqf(mv[i]);
                //printf("res mv = %f, %f\n", mv[i].x, mv[i].y);
            }
            //printf("inverse power iteration: %d error: %f\n", iMaxIte, s);
        }
        if (0 == iMaxIte && s > F(0.000001))
        {
            printf("warning: Inverse Power Max Iteration Reached, last s = %.20f\n", s);
        }
    }

public:

    __device__ __inline__ void EigenValues(CLGComplex* evs) const
    {
        deviceSUN<N, NoE> tmp = HessenbergC();
        tmp.FrancisQRIteration();
        for (INT i = 0; i < N; ++i)
        {
            evs[i] = tmp.m_me[i * N + i];
        }
    }

    __device__ __inline__ deviceSUN<N, NoE> EigenSystem(CLGComplex* evs) const
    {
        EigenValues(evs);
        CLGComplex tmpv[N];

        deviceSUN<N, NoE> eigenVectors;
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

    __device__ __inline__ deviceSUN<N, NoE> Log() const
    {
        CLGComplex eign[N];
        deviceSUN<N, NoE> egv = EigenSystem(eign);

        deviceSUN<N, NoE> ret = makeSUNZero();
        for (INT i = 0; i < N; ++i)
        {
            ret.m_me[i * N + i] = __cuClogf(eign[i]);
        }
        ret.MulDagger(egv);
        return egv.MulC(ret);
    }

    __device__ __inline__ deviceSUN<N, NoE> StrictExp() const
    {
        CLGComplex eign[N];
        deviceSUN<N, NoE> egv = EigenSystem(eign);

        deviceSUN<N, NoE> ret = makeSUNZero();
        for (INT i = 0; i < N; ++i)
        {
            ret.m_me[i * N + i] = __cuCexpf(eign[i]);
        }
        ret.MulDagger(egv);
        return egv.MulC(ret);
    }

#pragma endregion


    CLGComplex m_me[NoE];
};

typedef deviceSUN<4, 16> deviceSU4;
typedef deviceSUN<5, 32> deviceSU5;
typedef deviceSUN<6, 64> deviceSU6;
typedef deviceSUN<7, 64> deviceSU7;
typedef deviceSUN<8, 64> deviceSU8;


__END_NAMESPACE

#endif //#ifndef _SUN_H_

//=============================================================================
// END OF FILE
//=============================================================================
