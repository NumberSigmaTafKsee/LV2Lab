#pragma once

template<typename T>
// r = frac
// x = [i]
// y = [i+1]
T linear_interpolate(T x, T y, T r)
{        
    return x + r*(y-x);
}
template<typename T>
T cubic_interpolate(T finpos, T xm1, T x0, T x1, T x2)
{
    //T xm1 = x [inpos - 1];
    //T x0  = x [inpos + 0];
    //T x1  = x [inpos + 1];
    //T x2  = x [inpos + 2];
    T a = (3 * (x0-x1) - xm1 + x2) / 2;
    T b = 2*x1 + xm1 - (5*x0 + x2) / 2;
    T c = (x1 - xm1) / 2;
    return (((a * finpos) + b) * finpos + c) * finpos + x0;
}
// original
template<typename T>
// x = frac
// y0 = [i-1]
// y1 = [i]
// y2 = [i+1]
// y3 = [i+2]
T hermite1(T x, T y0, T y1, T y2, T y3)
{
    // 4-point, 3rd-order Hermite (x-form)
    T c0 = y1;
    T c1 = 0.5f * (y2 - y0);
    T c2 = y0 - 2.5f * y1 + 2.f * y2 - 0.5f * y3;
    T c3 = 1.5f * (y1 - y2) + 0.5f * (y3 - y0);
    return ((c3 * x + c2) * x + c1) * x + c0;
}

// james mccartney
template<typename T>
// x = frac
// y0 = [i-1]
// y1 = [i]
// y2 = [i+1]
// y3 = [i+2]
T hermite2(T x, T y0, T y1, T y2, T y3)
{
    // 4-point, 3rd-order Hermite (x-form)
    T c0 = y1;
    T c1 = 0.5f * (y2 - y0);
    T c3 = 1.5f * (y1 - y2) + 0.5f * (y3 - y0);
    T c2 = y0 - y1 + c1 - c3;
    return ((c3 * x + c2) * x + c1) * x + c0;
}

// james mccartney
template<typename T>
// x = frac
// y0 = [i-1]
// y1 = [i]
// y2 = [i+1]
// y3 = [i+2]
T hermite3(T x, T y0, T y1, T y2, T y3)
{
        // 4-point, 3rd-order Hermite (x-form)
        T c0 = y1;
        T c1 = 0.5f * (y2 - y0);
        T y0my1 = y0 - y1;
        T c3 = (y1 - y2) + 0.5f * (y3 - y0my1 - y2);
        T c2 = y0my1 + c1 - c3;

        return ((c3 * x + c2) * x + c1) * x + c0;
}

// laurent de soras
template<typename T>
// x[i-1]
// x[i]
// x[i+1]
// x[i+2]    
inline T hermite4(T frac_pos, T xm1, T x0, T x1, T x2)
{
    const T    c     = (x1 - xm1) * 0.5f;
    const T    v     = x0 - x1;
    const T    w     = c + v;
    const T    a     = w + v + (x2 - x0) * 0.5f;
    const T    b_neg = w + a;

    return ((((a * frac_pos) - b_neg) * frac_pos + c) * frac_pos + x0);
}

