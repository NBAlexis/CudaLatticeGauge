//=============================================================================
// FILENAME : ZN.h
//
// DESCRIPTION:
// deviceZN<N>: Z_N gauge group element (N-th roots of unity)
// Stored as 1x1 complex matrix, always projected to nearest root
//
// REVISION:
//  [05/16/2026 nbale]
//=============================================================================

#ifndef _ZN_H_
#define _ZN_H_

__BEGIN_NAMESPACE

template<INT N>
struct deviceZN
{

public:

    static constexpr INT matrixN = 1;
    static constexpr INT elementN = N;

    // Performance: default constructor leaves m_me uninitialized.
    // Caller must explicitly initialize via Id(), Zero(), or makeZN*() before use.
    __device__ deviceZN()
    {
    }

    __device__ deviceZN(const deviceZN<N>& other)
    {
        m_me = other.m_me;
    }

    __device__ deviceZN(const CLGComplex& other)
    {
        m_me = other;
    }

    __device__ void DebugPrint(const char* header = NULL) const
    {
        printf("%s%s{%f%s%fi};\n", NULL == header ? "" : header, NULL == header ? "" : "=",
            m_me.x, m_me.y < 0 ? "" : "+", m_me.y);
    }

#pragma region creation

    __device__ __inline__ static deviceZN<N> makeZNId()
    {
        deviceZN<N> ret;
        ret.m_me = _make_cuComplex(F(1.0), F(0.0));
        return ret;
    }

    __device__ __inline__ static deviceZN<N> makeAsK(UINT k)
    {
        deviceZN<N> ret;
        ret.SetFromIndex(k);
        return ret;
    }

    __device__ __inline__ void SetFromIndex(UINT idx)
    {
        Real theta = F(2.0) * PI * idx / N;
        m_me = _make_cuComplex(_cos(theta), _sin(theta));
    }

    __device__ __inline__ static deviceZN<N> makeZNRandom(UINT fatIndex)
    {
        deviceZN<N> ret;
        UINT k = static_cast<UINT>(_deviceRandomF(fatIndex) * N) % N;
        ret.SetFromIndex(k);
        return ret;
    }

    __device__ __inline__ static deviceZN<N> makeZNZero()
    {
        deviceZN<N> ret;
        ret.m_me = _make_cuComplex(F(0.0), F(0.0));
        return ret;
    }

    __device__ __inline__ void Id()
    {
        m_me = _make_cuComplex(F(1.0), F(0.0));
    }

    __device__ __inline__ void Zero()
    {
        m_me = _make_cuComplex(F(0.0), F(0.0));
    }

#pragma endregion

#pragma region operators

    __device__ __inline__ void Add(const deviceZN<N>& right)
    {
        m_me = _cuCaddf(m_me, right.m_me);
    }

    __device__ __inline__ void Sub(const deviceZN<N>& right)
    {
        m_me = _cuCsubf(m_me, right.m_me);
    }

    __device__ __inline__ void Mul(const deviceZN<N>& right)
    {
        m_me = _cuCmulf(m_me, right.m_me);
    }

    __device__ __inline__ void MulComp(const CLGComplex& right)
    {
        m_me = _cuCmulf(m_me, right);
    }

    __device__ __inline__ void MulReal(Real right)
    {
        m_me = cuCmulf_cr(m_me, right);
    }

    __device__ __inline__ void Dagger()
    {
        m_me.y = -m_me.y;
    }

    __device__ __inline__ void DaggerMul(const deviceZN<N>& right)
    {
        Dagger();
        Mul(right);
    }

    __device__ __inline__ void MulDagger(const deviceZN<N>& right)
    {
        m_me = _cuCmulf(m_me, _cuConjf(right.m_me));
    }

    __device__ __inline__ void Opposite()
    {
        m_me.x = -m_me.x;
        m_me.y = -m_me.y;
    }

    __device__ __inline__ deviceZN<N> AddC(const deviceZN<N>& right) const { deviceZN<N> ret(*this); ret.Add(right); return ret; }
    __device__ __inline__ deviceZN<N> SubC(const deviceZN<N>& right) const { deviceZN<N> ret(*this); ret.Sub(right); return ret; }
    __device__ __inline__ deviceZN<N> MulC(const deviceZN<N>& right) const { deviceZN<N> ret(*this); ret.Mul(right); return ret; }
    __device__ __inline__ deviceZN<N> MulCompC(const CLGComplex& right) const { deviceZN<N> ret(*this); ret.MulComp(right); return ret; }
    __device__ __inline__ deviceZN<N> MulRealC(const Real& right) const { deviceZN<N> ret(*this); ret.MulReal(right); return ret; }
    __device__ __inline__ deviceZN<N> DaggerC() const { deviceZN<N> ret(*this); ret.Dagger(); return ret; }
    __device__ __inline__ deviceZN<N> OppositeC() const { deviceZN<N> ret(*this); ret.Opposite(); return ret; }
    __device__ __inline__ deviceZN<N> DaggerMulC(const deviceZN<N>& right) const { deviceZN<N> ret(*this); ret.DaggerMul(right); return ret; }
    __device__ __inline__ deviceZN<N> MulDaggerC(const deviceZN<N>& right) const { deviceZN<N> ret(*this); ret.MulDagger(right); return ret; }

    __device__ __inline__ void AddDagger(const deviceZN<N>& right)
    {
        m_me = _cuCaddf(m_me, _cuConjf(right.m_me));
    }

    __device__ __inline__ void SubDagger(const deviceZN<N>& right)
    {
        m_me = _cuCsubf(m_me, _cuConjf(right.m_me));
    }

    __device__ __inline__ void AddReal(Real right)
    {
        m_me.x += right;
    }

    __device__ __inline__ void AddComp(const CLGComplex& right)
    {
        m_me = _cuCaddf(m_me, right);
    }

    __device__ __inline__ void SubReal(Real right)
    {
        m_me.x -= right;
    }

    __device__ __inline__ void SubComp(const CLGComplex& right)
    {
        m_me = _cuCsubf(m_me, right);
    }

    __device__ __inline__ void DivComp(const CLGComplex& right)
    {
        m_me = _cuCdivf(m_me, right);
    }

    __device__ __inline__ void DivReal(const Real& right)
    {
        m_me = cuCdivf_cr(m_me, right);
    }

    __device__ __inline__ deviceZN<N> AddCompC(const CLGComplex& right) const { deviceZN<N> ret(*this); ret.AddComp(right); return ret; }
    __device__ __inline__ deviceZN<N> AddRealC(const Real& right) const { deviceZN<N> ret(*this); ret.AddReal(right); return ret; }
    __device__ __inline__ deviceZN<N> SubCompC(const CLGComplex& right) const { deviceZN<N> ret(*this); ret.SubComp(right); return ret; }
    __device__ __inline__ deviceZN<N> SubRealC(const Real& right) const { deviceZN<N> ret(*this); ret.SubReal(right); return ret; }
    __device__ __inline__ deviceZN<N> DivCompC(const CLGComplex& right) const { deviceZN<N> ret(*this); ret.DivComp(right); return ret; }
    __device__ __inline__ deviceZN<N> DivRealC(const Real& right) const { deviceZN<N> ret(*this); ret.DivReal(right); return ret; }

    __device__ __inline__ void MulOnMe(const deviceZN<N>& left)
    {
        m_me = _cuCmulf(left.m_me, m_me);
    }

    __device__ __inline__ void MulOnMeDN(const deviceZN<N>& left)
    {
        m_me = _cuCmulf(_cuConjf(left.m_me), m_me);
    }

    __device__ __inline__ void MulOnMeND(const deviceZN<N>& left)
    {
        Dagger();
        MulOnMe(left);
    }

#pragma endregion

#pragma region useful functions

    __device__ __inline__ CLGComplex Tr() const
    {
        return m_me;
    }

    __device__ __inline__ Real ReTr() const
    {
        return m_me.x;
    }

    __device__ __inline__ Real ImTr() const
    {
        return m_me.y;
    }

    __device__ __inline__ void Re()
    {
        m_me.y = F(0.0);
    }

    __device__ __inline__ deviceZN<N> ReC() const { deviceZN<N> ret(*this); ret.Re(); return ret; }

    __device__ __inline__ void Im()
    {
        m_me.x = m_me.y;
        m_me.y = F(0.0);
    }

    __device__ __inline__ deviceZN<N> ImC() const { deviceZN<N> ret(*this); ret.Im(); return ret; }

    __device__ __inline__ void AddId()
    {
        m_me.x += F(1.0);
    }

    /**
     * Project to nearest N-th root of unity.
     * z = exp(2pi*i*round(N*arg(z)/(2pi))/N)
     */
    __device__ __inline__ void Proj(INT ite = 0)
    {
        UN_USE(ite);
        Real arg = _atan2(m_me.y, m_me.x);
        INT k = _round2int(arg * N / (F(2.0) * PI));
        k = ((k % N) + N) % N;
        Real theta = F(2.0) * PI * k / N;
        m_me = _make_cuComplex(_cos(theta), _sin(theta));
    }

    __device__ __inline__ void Norm()
    {
        Proj();
    }

    __device__ __inline__ deviceZN<N> Exp(const CLGComplex& a, BYTE uiPrecision = N + 1) const
    {
        deviceZN<N> tmp;
        for (BYTE i = 0; i < uiPrecision; ++i)
        {
            const Real exp_factor = __rcp(uiPrecision - i);
            CLGComplex alpha = cuCmulf_cr(a, exp_factor);
            deviceZN<N> aUoN(*this);
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

    __device__ __inline__ deviceZN<N> ExpReal(Real a, BYTE uiPrecision = N + 1) const
    {
        deviceZN<N> tmp;
        for (BYTE i = 0; i < uiPrecision; ++i)
        {
            deviceZN<N> aUoN(*this);
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

    __device__ __inline__ deviceZN<N> Log() const
    {
        deviceZN<N> ret;
        ret.m_me = _make_cuComplex(F(0.0), _atan2(m_me.y, m_me.x));
        return ret;
    }

    __device__ __inline__ deviceZN<N> StrictExp() const
    {
        deviceZN<N> ret;
        const Real ex = _exp(m_me.x);
        ret.m_me = _make_cuComplex(ex * _cos(m_me.y), ex * _sin(m_me.y));
        return ret;
    }

    // Ta/Th/Re2/iIm2 - for force computation compatibility
    // A = (U - U^dagger) / 2, then traceless anti-hermitian
    // For 1x1, this is just i*Im[U]
    __device__ __inline__ void Ta()
    {
        m_me = _make_cuComplex(F(0.0), m_me.y);
    }

    __device__ __inline__ void Th()
    {
        m_me = _make_cuComplex(m_me.x, F(0.0));
    }

    __device__ __inline__ void Re2()
    {
        m_me = _make_cuComplex(F(2.0) * m_me.x, F(0.0));
    }

    __device__ __inline__ void iIm2()
    {
        m_me = _make_cuComplex(F(0.0), F(2.0) * m_me.y);
    }

    __device__ __inline__ static Real TrIm(const deviceZN<N>& a, const deviceZN<N>& b)
    {
        return a.m_me.y * b.m_me.y;
    }

#pragma endregion

    CLGComplex m_me;
};

__END_NAMESPACE

#endif //#ifndef _ZN_H_

//=============================================================================
// END OF FILE
//=============================================================================
