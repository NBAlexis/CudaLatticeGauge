//=============================================================================
// FILENAME : SON.h
//
// DESCRIPTION:
// This is helper functions for calculate SO(N) real special orthogonal matrix.
// SO(N) = N x N real matrices with R^T R = I, det = +1.
//
// Layout:
// for N=4
// 0, 1, 2, 3
// 4, 5, 6, 7
// ...
//
// REVISION:
//  [05/29/2026 nbale]
//=============================================================================

#ifndef _SON_H_
#define _SON_H_

__BEGIN_NAMESPACE

template<INT N, INT NoE>
struct deviceSON
{

public:
    __device__ deviceSON()
    {

    }

    __device__ deviceSON(const deviceSON<N, NoE>& other)
    {
        memcpy(m_me, other.m_me, sizeof(Real) * N * N);
    }

    __device__ deviceSON(const Real* other)
    {
        memcpy(m_me, other, sizeof(Real) * N * N);
    }

    __device__ void DebugPrint(const char* header = NULL) const
    {
        printf("%s%s{{", NULL == header ? "" : header, NULL == header ? "" : "=");
        for (INT i = 0; i < N * N; ++i)
        {
            if ((N - 1) == (i % N))
            {
                printf(" %.6f },\n{", m_me[i]);
            }
            else
            {
                printf(" %.6f ,", m_me[i]);
            }
        }
        printf("}};\n");
    }

#pragma region creation

    __device__ __inline__ static deviceSON<N, NoE> makeSONZero()
    {
        deviceSON<N, NoE> ret;
        for (INT i = 0; i < N * N; ++i)
        {
            ret.m_me[i] = F(0.0);
        }
        return ret;
    }

    __device__ __inline__ void Zero()
    {
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i] = F(0.0);
        }
    }

    __device__ __inline__ static deviceSON<N, NoE> makeSONId()
    {
        deviceSON<N, NoE> ret;
        for (INT i = 0; i < N; ++i)
        {
            for (INT j = 0; j < N; ++j)
            {
                INT n = i * N + j;
                if (i == j)
                {
                    ret.m_me[n] = F(1.0);
                }
                else
                {
                    ret.m_me[n] = F(0.0);
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
                    m_me[n] = F(1.0);
                }
                else
                {
                    m_me[n] = F(0.0);
                }
            }
        }
    }

    /**
    * Generators of o(N) = so(N): real antisymmetric matrices.
    * N(N-1)/2 generators, indexed by pairs (i,j) with i < j.
    */
    __device__ __inline__ static deviceSON<N, NoE> makeSONGenerator(UINT uiGenerator)
    {
        deviceSON<N, NoE> ret = makeSONZero();
        UINT count = 0;
        for (UINT i = 0; i < N; ++i)
        {
            for (UINT j = i + 1; j < N; ++j)
            {
                if (count == uiGenerator)
                {
                    ret.m_me[i * N + j] = F(1.0);
                    ret.m_me[j * N + i] = F(-1.0);
                    return ret;
                }
                ++count;
            }
        }
        return ret;
    }

    /**
    * ret = r_a T_a, r_a are random real Gaussians, T_a are generators.
    */
    __device__ __inline__ static deviceSON<N, NoE> makeSONRandomGenerator(UINT fatIndex)
    {
        deviceSON<N, NoE> ret = makeSONZero();
        for (UINT i = 0; i < N * (N - 1) / 2; ++i)
        {
            const Real r = _deviceRandomGaussFSqrt2(fatIndex);
            ret.Add(makeSONGenerator(i).MulRealC(r));
        }
        return ret;
    }

    __device__ __inline__ static deviceSON<N, NoE> makeSONSumGenerator(Real fDivide)
    {
        deviceSON<N, NoE> ret = makeSONZero();
        for (UINT i = 0; i < N * (N - 1) / 2; ++i)
        {
            ret.Add(makeSONGenerator(i).MulRealC(fDivide));
        }
        return ret;
    }

    __device__ __inline__ static deviceSON<N, NoE> makeSONRandom(UINT fatIndex)
    {
        deviceSON<N, NoE> ret;
        for (UINT i = 0; i < N * N; ++i)
        {
            ret.m_me[i] = _deviceRandomF(fatIndex) - F(0.5);
        }
        ret.Norm();
        return ret;
    }

    __device__ __inline__ static deviceSON<N, NoE> makeSONRandomAny(UINT fatIndex)
    {
        deviceSON<N, NoE> ret;
        for (UINT i = 0; i < N * N; ++i)
        {
            ret.m_me[i] = _deviceRandomF(fatIndex) - F(0.5);
        }
        return ret;
    }

#pragma endregion

#pragma region operators

    __device__ __inline__ void Add(const deviceSON<N, NoE>& right)
    {
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i] += right.m_me[i];
        }
    }

    __device__ __inline__ void Sub(const deviceSON<N, NoE>& right)
    {
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i] -= right.m_me[i];
        }
    }

    __device__ __inline__ void AddReal(Real right)
    {
        for (INT i = 0; i < N; ++i)
        {
            m_me[i * N + i] += right;
        }
    }

    __device__ __inline__ void SubReal(Real right)
    {
        for (INT i = 0; i < N; ++i)
        {
            m_me[i * N + i] -= right;
        }
    }

    __device__ __inline__ void Mul(const deviceSON<N, NoE>& right)
    {
        Real temp[N];
        for (INT y = 0; y < N; ++y)
        {
            for (INT x = 0; x < N; ++x)
            {
                temp[x] = m_me[y * N] * right.m_me[x];
                for (INT i = 1; i < N; ++i)
                {
                    temp[x] += m_me[y * N + i] * right.m_me[i * N + x];
                }
            }
            memcpy(m_me + y * N, temp, sizeof(Real) * N);
        }
    }

    __device__ __inline__ void Mul(const deviceSON<N, NoE>& right, Real* temp)
    {
        for (INT y = 0; y < N; ++y)
        {
            for (INT x = 0; x < N; ++x)
            {
                temp[x] = m_me[y * N] * right.m_me[x];
                for (INT i = 1; i < N; ++i)
                {
                    temp[x] += m_me[y * N + i] * right.m_me[i * N + x];
                }
            }
            memcpy(m_me + y * N, temp, sizeof(Real) * N);
        }
    }

    __device__ __inline__ void TransposeMul(const deviceSON<N, NoE>& right)
    {
        Transpose();
        Mul(right);
    }

    __device__ __inline__ void MulTranspose(const deviceSON<N, NoE>& right)
    {
        Real temp[N];
        for (INT y = 0; y < N; ++y)
        {
            for (INT x = 0; x < N; ++x)
            {
                temp[x] = m_me[y * N] * right.m_me[x * N];
                for (INT i = 1; i < N; ++i)
                {
                    temp[x] += m_me[y * N + i] * right.m_me[x * N + i];
                }
            }
            memcpy(m_me + y * N, temp, sizeof(Real) * N);
        }
    }

    __device__ __inline__ void MulTranspose(const deviceSON<N, NoE>& right, Real* temp)
    {
        for (INT y = 0; y < N; ++y)
        {
            for (INT x = 0; x < N; ++x)
            {
                temp[x] = m_me[y * N] * right.m_me[x * N];
                for (INT i = 1; i < N; ++i)
                {
                    temp[x] += m_me[y * N + i] * right.m_me[x * N + i];
                }
            }
            memcpy(m_me + y * N, temp, sizeof(Real) * N);
        }
    }

    __device__ __inline__ void MulOnMe(const deviceSON<N, NoE>& left)
    {
        Real temp[N];
        for (INT x = 0; x < N; ++x)
        {
            for (INT y = 0; y < N; ++y)
            {
                temp[y] = left.m_me[y * N] * m_me[x];
                for (INT i = 1; i < N; ++i)
                {
                    temp[y] += left.m_me[y * N + i] * m_me[i * N + x];
                }
            }
            for (INT y = 0; y < N; ++y)
            {
                m_me[y * N + x] = temp[y];
            }
        }
    }

    __device__ __inline__ deviceSON<N, NoE> MulC(const deviceSON<N, NoE>& right) const
    {
        deviceSON<N, NoE> ret;
        for (INT x = 0; x < N; ++x)
        {
            for (INT y = 0; y < N; ++y)
            {
                ret.m_me[y * N + x] = m_me[y * N] * right.m_me[x];
                for (INT i = 1; i < N; ++i)
                {
                    ret.m_me[y * N + x] += m_me[y * N + i] * right.m_me[i * N + x];
                }
            }
        }
        return ret;
    }

    __device__ __inline__ void MulReal(Real right)
    {
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i] *= right;
        }
    }

    __device__ __inline__ void DivReal(const Real& right)
    {
        Real inv = __rcp(right);
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i] *= inv;
        }
    }

    __device__ __inline__ deviceSON<N, NoE> MulTransposeC(const deviceSON<N, NoE>& right) const
    {
        deviceSON<N, NoE> ret;
        for (INT x = 0; x < N; ++x)
        {
            for (INT y = 0; y < N; ++y)
            {
                ret.m_me[y * N + x] = m_me[y * N] * right.m_me[x * N];
                for (INT i = 1; i < N; ++i)
                {
                    ret.m_me[y * N + x] += m_me[y * N + i] * right.m_me[x * N + i];
                }
            }
        }
        return ret;
    }

    __device__ __inline__ deviceSON<N, NoE> TransposeMulC(const deviceSON<N, NoE>& right) const
    {
        deviceSON<N, NoE> ret;
        for (INT x = 0; x < N; ++x)
        {
            for (INT y = 0; y < N; ++y)
            {
                ret.m_me[y * N + x] = m_me[y] * right.m_me[x];
                for (INT i = 1; i < N; ++i)
                {
                    ret.m_me[y * N + x] += m_me[i * N + y] * right.m_me[i * N + x];
                }
            }
        }
        return ret;
    }

    __device__ __inline__ deviceSON<N, NoE> AddC(const deviceSON<N, NoE>& right) const { deviceSON<N, NoE> ret(*this); ret.Add(right); return ret; }
    __device__ __inline__ deviceSON<N, NoE> AddRealC(const Real& right) const { deviceSON<N, NoE> ret(*this); ret.AddReal(right); return ret; }

    __device__ __inline__ deviceSON<N, NoE> SubC(const deviceSON<N, NoE>& right) const { deviceSON<N, NoE> ret(*this); ret.Sub(right); return ret; }
    __device__ __inline__ deviceSON<N, NoE> SubRealC(const Real& right) const { deviceSON<N, NoE> ret(*this); ret.SubReal(right); return ret; }

    __device__ __inline__ deviceSON<N, NoE> MulRealC(const Real& right) const { deviceSON<N, NoE> ret(*this); ret.MulReal(right); return ret; }
    __device__ __inline__ deviceSON<N, NoE> DivRealC(const Real& right) const { deviceSON<N, NoE> ret(*this); ret.DivReal(right); return ret; }

    __device__ __inline__ void BlockMul(
        const deviceSON<N, NoE>& right,
        UINT iWidth,
        UINT iLeftStartX, UINT iLeftStartY,
        UINT iRightStartX, UINT iRightStartY)
    {
        Real temp[N];
        for (INT y = 0; y < iWidth; ++y)
        {
            for (INT x = 0; x < iWidth; ++x)
            {
                temp[x] = m_me[(y + iLeftStartY) * N + iLeftStartX] * right.m_me[iRightStartY * N + x + iRightStartX];
                for (INT i = 1; i < iWidth; ++i)
                {
                    temp[x] += m_me[(y + iLeftStartY) * N + (i + iLeftStartX)] * right.m_me[(i + iRightStartY) * N + x + iRightStartX];
                }
            }
            memcpy(m_me + (y + iLeftStartY) * N + iLeftStartX, temp, sizeof(Real) * iWidth);
        }
    }

    __device__ __inline__ void BlockCopy(const deviceSON<N, NoE>& right,
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
            ret.m_ve[y] = _make_cuComplex(m_me[y * N] * v.m_ve[0].x, m_me[y * N] * v.m_ve[0].y);
            for (INT x = 1; x < N; ++x)
            {
                ret.m_ve[y].x += m_me[y * N + x] * v.m_ve[x].x;
                ret.m_ve[y].y += m_me[y * N + x] * v.m_ve[x].y;
            }
        }
        return ret;
    }

    __device__ __inline__ deviceSU3Vector MulVector(const deviceSU3Vector& v) const
    {
        if constexpr (3 == N)
        {
            deviceSU3Vector ret;
            for (INT y = 0; y < N; ++y)
            {
                ret.m_ve[y] = _make_cuComplex(m_me[y * N] * v.m_ve[0].x, m_me[y * N] * v.m_ve[0].y);
                for (INT x = 1; x < N; ++x)
                {
                    ret.m_ve[y].x += m_me[y * N + x] * v.m_ve[x].x;
                    ret.m_ve[y].y += m_me[y * N + x] * v.m_ve[x].y;
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
                ret.m_ve[y] = _make_cuComplex(m_me[y * N] * v.m_ve[0].x, m_me[y * N] * v.m_ve[0].y);
                for (INT x = 1; x < N; ++x)
                {
                    ret.m_ve[y].x += m_me[y * N + x] * v.m_ve[x].x;
                    ret.m_ve[y].y += m_me[y * N + x] * v.m_ve[x].y;
                }
            }
            return ret;
        }
        return deviceSU2Vector::makeZeroSU2Vector();
    }

    template<INT NofVE>
    __device__ __inline__ deviceSUNVector<N, NofVE> TransposeMulVector(const deviceSUNVector<N, NofVE>& v) const
    {
        deviceSUNVector<N, NofVE> ret;
        for (INT y = 0; y < N; ++y)
        {
            ret.m_ve[y] = _make_cuComplex(m_me[y] * v.m_ve[0].x, m_me[y] * v.m_ve[0].y);
            for (INT x = 1; x < N; ++x)
            {
                ret.m_ve[y].x += m_me[x * N + y] * v.m_ve[x].x;
                ret.m_ve[y].y += m_me[x * N + y] * v.m_ve[x].y;
            }
        }
        return ret;
    }

    __device__ __inline__ deviceSU3Vector TransposeMulVector(const deviceSU3Vector& v) const
    {
        if constexpr (3 == N)
        {
            deviceSU3Vector ret;
            for (INT y = 0; y < N; ++y)
            {
                ret.m_ve[y] = _make_cuComplex(m_me[y] * v.m_ve[0].x, m_me[y] * v.m_ve[0].y);
                for (INT x = 1; x < N; ++x)
                {
                    ret.m_ve[y].x += m_me[x * N + y] * v.m_ve[x].x;
                    ret.m_ve[y].y += m_me[x * N + y] * v.m_ve[x].y;
                }
            }
            return ret;
        }
        return deviceSU3Vector::makeZeroSU3Vector();
    }

    __device__ __inline__ deviceSU2Vector TransposeMulVector(const deviceSU2Vector& v) const
    {
        if constexpr (2 == N)
        {
            deviceSU2Vector ret;
            for (INT y = 0; y < N; ++y)
            {
                ret.m_ve[y] = _make_cuComplex(m_me[y] * v.m_ve[0].x, m_me[y] * v.m_ve[0].y);
                for (INT x = 1; x < N; ++x)
                {
                    ret.m_ve[y].x += m_me[x * N + y] * v.m_ve[x].x;
                    ret.m_ve[y].y += m_me[x * N + y] * v.m_ve[x].y;
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
            res[y] = _make_cuComplex(m_me[y * N] * v[0].x, m_me[y * N] * v[0].y);
            for (INT x = 1; x < N; ++x)
            {
                res[y].x += m_me[y * N + x] * v[x].x;
                res[y].y += m_me[y * N + x] * v[x].y;
            }
        }
    }

    template<INT NoVE>
    __device__ __inline__ static deviceSON<N, NoE> makeSONContractV(const deviceSUNVector<N, NoVE>& left, const deviceSUNVector<N, NoVE>& right)
    {
        deviceSON<N, NoE> ret;
        for (INT y = 0; y < N; ++y)
        {
            for (INT x = 0; x < N; ++x)
            {
                ret.m_me[y * N + x] = left.m_ve[x].x * right.m_ve[y].x + left.m_ve[x].y * right.m_ve[y].y;
            }
        }
        return ret;
    }

    __device__ __inline__ static deviceSON<N, NoE> makeSONContractV(const deviceSU3Vector& left, const deviceSU3Vector& right)
    {
        if constexpr (3 == N)
        {
            deviceSON<N, NoE> ret;
            for (INT y = 0; y < N; ++y)
            {
                for (INT x = 0; x < N; ++x)
                {
                    ret.m_me[y * N + x] = left.m_ve[x].x * right.m_ve[y].x + left.m_ve[x].y * right.m_ve[y].y;
                }
            }
            return ret;
        }
        return deviceSON<N, NoE>::makeSONZero();
    }

    __device__ __inline__ static deviceSON<N, NoE> makeSONContractV(const deviceSU2Vector& left, const deviceSU2Vector& right)
    {
        if constexpr (2 == N)
        {
            deviceSON<N, NoE> ret;
            for (INT y = 0; y < N; ++y)
            {
                for (INT x = 0; x < N; ++x)
                {
                    ret.m_me[y * N + x] = left.m_ve[x].x * right.m_ve[y].x + left.m_ve[x].y * right.m_ve[y].y;
                }
            }
            return ret;
        }
        return deviceSON<N, NoE>::makeSONZero();
    }

#pragma endregion

#pragma region useful functions

public:

    __device__ __inline__ Real Determinent() const
    {
        deviceSON<N, NoE> tmp(*this);
        tmp.LUNoReturn();
        Real det = tmp.m_me[0];
        for (INT i = 1; i < N; ++i)
        {
            det *= tmp.m_me[i * N + i];
        }
        return det;
    }

    __device__ __inline__ Real Tr() const
    {
        Real ret = m_me[0];
        for (INT i = 1; i < N; ++i)
        {
            ret += m_me[i * N + i];
        }
        return ret;
    }

    __device__ __inline__ void Transpose()
    {
        for (INT x = 0; x < N; ++x)
        {
            for (INT y = x + 1; y < N; ++y)
            {
                Real tmp = m_me[y * N + x];
                m_me[y * N + x] = m_me[x * N + y];
                m_me[x * N + y] = tmp;
            }
        }
    }

    __device__ __inline__ deviceSON<N, NoE> TransposeC() const
    {
        deviceSON<N, NoE> ret;
        for (INT x = 0; x < N; ++x)
        {
            for (INT y = 0; y < N; ++y)
            {
                ret.m_me[y * N + x] = m_me[x * N + y];
            }
        }
        return ret;
    }

    __device__ __inline__ void Opposite()
    {
        for (INT i = 0; i < N * N; ++i)
        {
            m_me[i] = -m_me[i];
        }
    }

    __device__ __inline__ deviceSON<N, NoE> OppositeC() const
    {
        deviceSON<N, NoE> ret(*this);
        ret.Opposite();
        return ret;
    }

    __device__ __inline__ void AddId()
    {
        for (INT i = 0; i < N; ++i)
        {
            m_me[i * N + i] += F(1.0);
        }
    }

    __device__ __inline__ void Norm()
    {
        QR();
        // Ensure det = +1 for SO(N)
        Real det = Determinent();
        if (det < 0)
        {
            for (INT i = 0; i < N; ++i)
            {
                m_me[i * N + (N - 1)] = -m_me[i * N + (N - 1)];
            }
        }
    }

    __device__ __inline__ void Proj()
    {
        Norm();
    }

    /**
    * U' = exp(aU) = (1 + a U + a^2 U^2/2 + ... + a^N U^N/N!)
    *    = 1 + a U (1 + a U /2 (1 + a U/3 ...))
    */
    __device__ __inline__ deviceSON<N, NoE> ExpReal(Real a, BYTE uiPrecision = N + 1) const
    {
        deviceSON<N, NoE> tmp;
        for (BYTE i = 0; i < uiPrecision; ++i)
        {
            deviceSON<N, NoE> aUoN(*this);
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

    __device__ __inline__ void Ta()
    {
        for (INT y = 0; y < N; ++y)
        {
            for (INT x = y + 1; x < N; ++x)
            {
                Real val = F(0.5) * (m_me[y * N + x] - m_me[x * N + y]);
                m_me[y * N + x] = val;
                m_me[x * N + y] = -val;
            }
            m_me[y * N + y] = F(0.0);
        }
    }

    /**
     * Approximate log for special orthogonal matrices: log(O) ~= (O - O^T)/2
     * Exact for O near identity.
     */
    __device__ __inline__ deviceSON<N, NoE> Log() const
    {
        deviceSON<N, NoE> res(*this);
        res.Ta();
        return res;
    }

private:

    __device__ __inline__ void OneStepHouseHolderQR(deviceSON<N, NoE>& A, UINT i)
    {
        Real fLength = F(0.0);
        for (INT y = i; y < N; ++y)
        {
            fLength += A.m_me[y * N + i] * A.m_me[y * N + i];
        }

        if (fLength > _CLG_FLT_MIN_)
        {
            fLength = _sqrt(fLength);
            Real vlst[N];
            Real fULength = F(0.0);
            for (UINT y = i; y < N; ++y)
            {
                vlst[y] = A.m_me[y * N + i];
                if (y == i)
                {
                    vlst[y] += (vlst[y] >= 0 ? fLength : -fLength);
                }
                fULength += vlst[y] * vlst[y];
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
                            m_me[py * N + px] = vlst[py] * vlst[px];
                            if (px == py)
                            {
                                m_me[py * N + px] = F(1.0) - fBeta * m_me[py * N + px];
                            }
                            else
                            {
                                m_me[py * N + px] = -fBeta * m_me[py * N + px];
                            }
                        }
                    }

                    A.MulOnMe(*this);
                }
                else
                {
                    deviceSON<N, NoE> tmp = makeSONId();
                    for (INT px = i; px < N; ++px)
                    {
                        for (INT py = i; py < N; ++py)
                        {
                            tmp.m_me[py * N + px] = vlst[py] * vlst[px];
                            if (px == py)
                            {
                                tmp.m_me[py * N + px] = F(1.0) - fBeta * tmp.m_me[py * N + px];
                            }
                            else
                            {
                                tmp.m_me[py * N + px] = -fBeta * tmp.m_me[py * N + px];
                            }
                        }
                    }

                    A.MulOnMe(tmp);
                    Mul(tmp);
                }
            }
        }
    }

public:

    __device__ __inline__ deviceSON<N, NoE> QR()
    {
        deviceSON<N, NoE> r(*this);
        for (UINT i = 0; i < (N - 1); ++i)
        {
            OneStepHouseHolderQR(r, i);
        }
        return r;
    }

private:

    __device__ __inline__ void OneStepGaussianLU(deviceSON<N, NoE>& L, UINT i)
    {
        Real dobedivide = m_me[i * N + i];
        if (dobedivide * dobedivide > _CLG_FLT_MIN_)
        {
            for (INT y = i + 1; y < N; ++y)
            {
                L.m_me[y * N + i] = m_me[y * N + i] / dobedivide;
                m_me[y * N + i] = F(0.0);
            }

            for (INT x = i + 1; x < N; ++x)
            {
                for (INT y = i + 1; y < N; ++y)
                {
                    m_me[y * N + x] -= L.m_me[y * N + i] * m_me[i * N + x];
                }
            }
        }
    }

    __device__ __inline__ void OneStepGaussianLUNoReturn(UINT i)
    {
        Real dobedivide = m_me[i * N + i];
        if (dobedivide * dobedivide > _CLG_FLT_MIN_)
        {
            for (INT y = i + 1; y < N; ++y)
            {
                m_me[y * N + i] /= dobedivide;
            }

            for (INT x = i + 1; x < N; ++x)
            {
                for (INT y = i + 1; y < N; ++y)
                {
                    m_me[y * N + x] -= m_me[y * N + i] * m_me[i * N + x];
                }
            }
        }
    }

public:

    __device__ __inline__ deviceSON<N, NoE> LU()
    {
        deviceSON<N, NoE> L = makeSONId();
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

    __device__ __inline__ void BackwardSubstitutionElement(Real* v) const
    {
        for (INT i = N - 1; i >= 0; --i)
        {
            for (INT j = i + 1; j < N; ++j)
            {
                v[i] -= m_me[i * N + j] * v[j];
            }
            if (m_me[i * N + i] * m_me[i * N + i] > _CLG_FLT_MIN_)
            {
                v[i] /= m_me[i * N + i];
            }
        }
    }

    template<INT NoVE>
    __device__ __inline__ void BackwardSubstitution(deviceSUNVector<N, NoVE>& v) const
    {
        for (INT i = N - 1; i >= 0; --i)
        {
            for (INT j = i + 1; j < N; ++j)
            {
                v.m_ve[i].x -= m_me[i * N + j] * v.m_ve[j].x;
                v.m_ve[i].y -= m_me[i * N + j] * v.m_ve[j].y;
            }
            if (m_me[i * N + i] * m_me[i * N + i] > _CLG_FLT_MIN_)
            {
                v.m_ve[i].x /= m_me[i * N + i];
                v.m_ve[i].y /= m_me[i * N + i];
            }
        }
    }

private:

    __device__ __inline__ void OneStepHouseHolderHessenberg(deviceSON<N, NoE>& A, UINT i)
    {
        Real fLength = F(0.0);
        for (INT y = i + 1; y < N; ++y)
        {
            fLength += A.m_me[y * N + i] * A.m_me[y * N + i];
        }

        if (fLength > _CLG_FLT_MIN_)
        {
            fLength = _sqrt(fLength);
            Real vlst[N];
            Real fULength = F(0.0);
            for (INT y = i + 1; y < N; ++y)
            {
                vlst[y] = A.m_me[y * N + i];
                if (y == i + 1)
                {
                    vlst[y] += (vlst[y] >= 0 ? fLength : -fLength);
                }
                fULength += vlst[y] * vlst[y];
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
                                m_me[py * N + px] = F(1.0);
                            }
                            else if (px > 0 && py > 0)
                            {
                                m_me[py * N + px] = vlst[py] * vlst[px];
                                if (px == py)
                                {
                                    m_me[py * N + px] = F(1.0) - fBeta * m_me[py * N + px];
                                }
                                else
                                {
                                    m_me[py * N + px] = -fBeta * m_me[py * N + px];
                                }
                            }
                            else
                            {
                                m_me[py * N + px] = F(0.0);
                            }
                        }
                    }

                    A.MulOnMe(*this);
                    A.Mul(*this);
                }
                else
                {
                    deviceSON<N, NoE> tmp = makeSONId();
                    for (INT px = i + 1; px < N; ++px)
                    {
                        for (INT py = i + 1; py < N; ++py)
                        {
                            tmp.m_me[py * N + px] = vlst[py] * vlst[px];
                            if (px == py)
                            {
                                tmp.m_me[py * N + px] = F(1.0) - fBeta * tmp.m_me[py * N + px];
                            }
                            else
                            {
                                tmp.m_me[py * N + px] = -fBeta * tmp.m_me[py * N + px];
                            }
                        }
                    }

                    A.MulOnMe(tmp);
                    A.Mul(tmp);
                    Mul(tmp);
                }
            }
        }
    }

    __device__ __inline__ void OneStepHouseHolderHessenbergC(deviceSON<N, NoE>& A, UINT i) const
    {
        Real fLength = F(0.0);
        for (INT y = i + 1; y < N; ++y)
        {
            fLength += A.m_me[y * N + i] * A.m_me[y * N + i];
        }

        if (fLength > _CLG_FLT_MIN_)
        {
            fLength = _sqrt(fLength);
            Real vlst[N];
            Real fULength = F(0.0);
            for (INT y = i + 1; y < N; ++y)
            {
                vlst[y] = A.m_me[y * N + i];
                if (y == i + 1)
                {
                    vlst[y] += (vlst[y] >= 0 ? fLength : -fLength);
                }
                fULength += vlst[y] * vlst[y];
            }

            if (fULength > _CLG_FLT_MIN_)
            {
                const Real fBeta = F(2.0) / fULength;
                deviceSON<N, NoE> tmp = makeSONId();
                for (INT px = i + 1; px < N; ++px)
                {
                    for (INT py = i + 1; py < N; ++py)
                    {
                        tmp.m_me[py * N + px] = vlst[py] * vlst[px];
                        if (px == py)
                        {
                            tmp.m_me[py * N + px] = F(1.0) - fBeta * tmp.m_me[py * N + px];
                        }
                        else
                        {
                            tmp.m_me[py * N + px] = -fBeta * tmp.m_me[py * N + px];
                        }
                    }
                }

                A.MulOnMe(tmp);
                A.Mul(tmp);
            }
        }
    }

public:

    __device__ __inline__ deviceSON<N, NoE> Hessenberg()
    {
        deviceSON<N, NoE> r(*this);
        for (UINT i = 0; i < (N - 2); ++i)
        {
            OneStepHouseHolderHessenberg(r, i);
        }
        return r;
    }

    __device__ __inline__ deviceSON<N, NoE> HessenbergC() const
    {
        deviceSON<N, NoE> r(*this);
        for (UINT i = 0; i < (N - 2); ++i)
        {
            OneStepHouseHolderHessenbergC(r, i);
        }
        return r;
    }

    #pragma region QR Iteration

    static __device__ __inline__ void checkMatrixIndexDoubleShift(Real* mtr, INT* decomp, INT dx)
    {
        decomp[0] = 0;
        decomp[1] = dx;

        for (INT i = dx - 2; i >= 0; --i)
        {
            if (mtr[(i + 1) * dx + i] * mtr[(i + 1) * dx + i] < _CLG_FLT_MIN_)
            {
                mtr[(i + 1) * dx + i] = F(0.0);

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

    static __device__ __inline__ void threeHouseHolder(Real& a, Real& b, Real& c)
    {
        Real len = a * a + b * b + c * c;
        if (len < _CLG_FLT_MIN_)
        {
            a = F(0.0);
            b = F(0.0);
            c = F(0.0);
            return;
        }
        len = _sqrt(len);
        a += len;

        Real len2 = F(0.5) * (a * a + b * b + c * c);
        if (len2 < _CLG_FLT_MIN_)
        {
            a = F(0.0);
            b = F(0.0);
            c = F(0.0);
            return;
        }
        len2 = __div(F(1.0), _sqrt(len2));
        a *= len2;
        b *= len2;
        c *= len2;
    }

    static __device__ __inline__ void twoHouseHolder(Real& a, Real& b)
    {
        Real len = a * a + b * b;
        if (len < _CLG_FLT_MIN_)
        {
            a = F(0.0);
            b = F(0.0);
            return;
        }
        len = _sqrt(len);
        a += len;

        Real len2 = F(0.5) * (a * a + b * b);
        if (len2 < _CLG_FLT_MIN_)
        {
            a = F(0.0);
            b = F(0.0);
            return;
        }
        len2 = __div(F(1.0), _sqrt(len2));
        a *= len2;
        b *= len2;
    }

    static __device__ __inline__ void FrancisQRStep(Real* H, INT dm)
    {
        Real s = H[dm * dm - dm - 2] + H[dm * dm - 1];
        Real t = H[dm * dm - dm - 2] * H[dm * dm - 1] - H[dm * dm - dm - 1] * H[dm * dm - 2];

        Real x = H[0] * H[0] + H[1] * H[dm] - s * H[0] + t;
        Real y = H[dm] * (H[0] + H[dm + 1] - s);
        Real z = H[dm] * H[2 * dm + 1];

        Real u[9];
        for (INT k = 0; k <= dm - 3; ++k)
        {
            threeHouseHolder(x, y, z);

            u[0] = x * x; u[1] = y * x; u[2] = z * x;
            u[3] = x * y; u[4] = y * y; u[5] = z * y;
            u[6] = x * z; u[7] = y * z; u[8] = z * z;

            u[0] -= F(1.0);
            u[4] -= F(1.0);
            u[8] -= F(1.0);

            INT q = (k < 1) ? 1 : k;
            for (INT hx = (q - 1); hx < dm; ++hx)
            {
                const Real newhy1 = u[0] * H[k * dm + hx] + u[1] * H[(k + 1) * dm + hx] + u[2] * H[(k + 2) * dm + hx];
                const Real newhy2 = u[3] * H[k * dm + hx] + u[4] * H[(k + 1) * dm + hx] + u[5] * H[(k + 2) * dm + hx];
                const Real newhy3 = u[6] * H[k * dm + hx] + u[7] * H[(k + 1) * dm + hx] + u[8] * H[(k + 2) * dm + hx];

                H[k * dm + hx] = -newhy1;
                H[(k + 1) * dm + hx] = -newhy2;
                H[(k + 2) * dm + hx] = -newhy3;
            }

            q = k + 3;
            if (q >= dm)
            {
                q = dm - 1;
            }
            for (INT hy = 0; hy <= q; ++hy)
            {
                const Real newhx1 = u[0] * H[hy * dm + k] + u[3] * H[hy * dm + k + 1] + u[6] * H[hy * dm + k + 2];
                const Real newhx2 = u[1] * H[hy * dm + k] + u[4] * H[hy * dm + k + 1] + u[7] * H[hy * dm + k + 2];
                const Real newhx3 = u[2] * H[hy * dm + k] + u[5] * H[hy * dm + k + 1] + u[8] * H[hy * dm + k + 2];

                H[hy * dm + k] = -newhx1;
                H[hy * dm + k + 1] = -newhx2;
                H[hy * dm + k + 2] = -newhx3;
            }

            x = H[(k + 1) * dm + k];
            y = H[(k + 2) * dm + k];
            if (k < dm - 3)
            {
                z = H[(k + 3) * dm + k];
            }
        }

        twoHouseHolder(x, y);
        u[0] = x * x; u[1] = y * x; u[2] = x * y; u[3] = y * y;

        u[0] -= F(1.0);
        u[3] -= F(1.0);

        INT startX = (2 == dm) ? 0 : (dm - 3);
        for (INT hx = startX; hx < dm; ++hx)
        {
            const Real newhy1 = u[0] * H[(dm - 2) * dm + hx] + u[1] * H[(dm - 1) * dm + hx];
            const Real newhy2 = u[2] * H[(dm - 2) * dm + hx] + u[3] * H[(dm - 1) * dm + hx];

            H[(dm - 2) * dm + hx] = -newhy1;
            H[(dm - 1) * dm + hx] = -newhy2;
        }
        for (INT hy = 0; hy < dm; ++hy)
        {
            const Real newhx1 = u[0] * H[hy * dm + dm - 2] + u[2] * H[hy * dm + dm - 1];
            const Real newhx2 = u[1] * H[hy * dm + dm - 2] + u[3] * H[hy * dm + dm - 1];

            H[hy * dm + dm - 2] = -newhx1;
            H[hy * dm + dm - 1] = -newhx2;
        }
    }

    static __device__ __inline__ void CalculateEigen2x2(
        Real& h00, Real& h01,
        Real& h10, Real& h11,
        CLGComplex* evs, INT& evCount)
    {
        Real trace = h00 + h11;
        Real det2 = h00 * h11 - h01 * h10;
        Real disc = trace * trace - F(4.0) * det2;
        if (disc >= F(0.0))
        {
            Real s = _sqrt(disc);
            evs[evCount++] = _make_cuComplex(F(0.5) * (trace + s), F(0.0));
            evs[evCount++] = _make_cuComplex(F(0.5) * (trace - s), F(0.0));
        }
        else
        {
            Real s = _sqrt(-disc);
            evs[evCount++] = _make_cuComplex(F(0.5) * trace, F(0.5) * s);
            evs[evCount++] = _make_cuComplex(F(0.5) * trace, -F(0.5) * s);
        }
        h10 = F(0.0);
    }

    __device__ __inline__ void FrancisQRIteration()
    {
        Real submatrix[N * N];
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
                INT offset = decomp[0];
                CLGComplex dummy[2];
                INT dummyCount = 0;
                CalculateEigen2x2(
                    m_me[offset * N + offset],       m_me[offset * N + offset + 1],
                    m_me[(offset + 1) * N + offset], m_me[(offset + 1) * N + offset + 1],
                    dummy, dummyCount);
            }
            else
            {
                for (INT y = 0; y < dm; ++y)
                {
                    memcpy(submatrix + y * dm, m_me + (decomp[0] + y) * N + decomp[0], sizeof(Real) * dm);
                }

                FrancisQRStep(submatrix, dm);

                for (INT y = 0; y < dm; ++y)
                {
                    memcpy(m_me + (decomp[0] + y) * N + decomp[0], submatrix + y * dm, sizeof(Real) * dm);
                }
            }
        }
    }

    #pragma endregion

public:

    __device__ __inline__ void EigenValues(CLGComplex* evs) const
    {
        deviceSON<N, NoE> tmp = HessenbergC();
        tmp.FrancisQRIteration();
        INT k = 0;
        for (INT i = 0; i < N; ++i)
        {
            if (i < N - 1 && tmp.m_me[(i + 1) * N + i] * tmp.m_me[(i + 1) * N + i] > _CLG_FLT_MIN_)
            {
                Real a = tmp.m_me[i * N + i];
                Real b = tmp.m_me[i * N + i + 1];
                Real c = tmp.m_me[(i + 1) * N + i];
                Real d = tmp.m_me[(i + 1) * N + i + 1];
                Real trace = a + d;
                Real det2 = a * d - b * c;
                Real disc = trace * trace - F(4.0) * det2;
                if (disc >= F(0.0))
                {
                    Real s = _sqrt(disc);
                    evs[k++] = _make_cuComplex(F(0.5) * (trace + s), F(0.0));
                    evs[k++] = _make_cuComplex(F(0.5) * (trace - s), F(0.0));
                }
                else
                {
                    Real s = _sqrt(-disc);
                    evs[k++] = _make_cuComplex(F(0.5) * trace, F(0.5) * s);
                    evs[k++] = _make_cuComplex(F(0.5) * trace, -F(0.5) * s);
                }
                ++i;
            }
            else
            {
                evs[k++] = _make_cuComplex(tmp.m_me[i * N + i], F(0.0));
            }
        }
    }

#pragma endregion

    Real m_me[NoE];
};

#define _TYPEDEFSON(n, moe) typedef deviceSON<n, moe> deviceSO##n;
#define _DEF_F2_TO_SON(n, imp) _DEF_F2SON_N(n, imp, 256, 256, 256, 256, 256, 128, 128, 128, 64, 64, 64, 32, 16, 16, 4)

_DEF_F2_TO_SON(_MAX_SON, _TYPEDEFSON)

__END_NAMESPACE

#endif //#ifndef _SON_H_

//=============================================================================
// END OF FILE
//=============================================================================
