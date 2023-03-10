%{
#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include <random> 

#include <kfr/kfr.h>
#include <kfr/dft.hpp>
#include <kfr/io.hpp>
#include <kfr/math.hpp>
#include "kfrcore.hpp"
%}

namespace kfr {
    using b8   = bool;
    using f32  = float;
    using f64  = double;
    using i8   = int8_t;
    using i16  = int16_t;
    using i32  = int32_t;
    using i64  = int64_t;
    using u8   = uint8_t;
    using u16  = uint16_t;
    using u32  = uint32_t;
    using u64  = uint64_t;
    using umax = uint64_t;
    using imax = int64_t;
    using fmax = double;
    using f80  = long double;
    using fbase = SampleType;
}

%inline %{
    typedef int8_t   i8;
    typedef uint8_t   u8;
    typedef int16_t  i16;
    typedef uint16_t  u16;
    typedef int32_t  i32;
    typedef uint32_t u32;
    typedef signed long   i64;
    typedef unsigned long  u64;
    typedef float    f32;
    typedef double   f64;
%};

namespace kfr {

    enum class audio_sample_type
    {
        unknown,
        i8,
        i16,
        i24,
        i32,
        i64,
        f32,
        f64,
        first_float = f32
    };

    enum class sample_rate_conversion_quality : int
    {
        draft   = 4,
        low     = 6,
        normal  = 8,
        high    = 10,
        perfect = 12,
    };

    enum class biquad_type
    {
        lowpass,
        highpass,
        bandpass,
        bandstop,
        peak,
        notch,
        lowshelf,
        highshelf
    };

    enum class Speaker : int
    {
        None          = -1,
        Mono          = 0,
        M             = static_cast<int>(Mono),
        Left          = 1,
        L             = static_cast<int>(Left),
        Right         = 2,
        R             = static_cast<int>(Right),
        Center        = 3,
        C             = static_cast<int>(Center),
        Lfe           = 4,
        Ls            = 5,
        LeftSurround  = static_cast<int>(Ls),
        Rs            = 6,
        RightSurround = static_cast<int>(Rs),
        Lc            = 7,
        Rc            = 8,
        S             = 9,
        Cs            = static_cast<int>(S),
        Sl            = 10,
        Sr            = 11,
        Tm            = 12,
        Tfl           = 13,
        Tfc           = 14,
        Tfr           = 15,
        Trl           = 16,
        Trc           = 17,
        Trr           = 18,
        Lfe2          = 19
    };

    enum class SpeakerArrangement : int
    {
        None           = -1,
        Mono           = 0,
        Stereo         = 1,
        StereoSurround = 2,
        StereoCenter   = 3,
        StereoSide     = 4,
        StereoCLfe     = 5,
        Cine30         = 6,
        Music30        = 7,
        Cine31         = 8,
        Music31        = 9,
        Cine40         = 10,
        Music40        = 11,
        Cine41         = 12,
        Music41        = 13,
        Arr50          = 14,
        Arr51          = 15,
        Cine60         = 16,
        Music60        = 17,
        Cine61         = 18,
        Music61        = 19,
        Cine70         = 20,
        Music70        = 21,
        Cine71         = 22,
        Music71        = 23,
        Cine80         = 24,
        Music80        = 25,
        Cine81         = 26,
        Music81        = 27,
        Arr102         = 28
    };

    /// @brief Seek origin
    enum class seek_origin : int
    {
        current = SEEK_CUR, ///< From the current position
        begin   = SEEK_SET, ///< From the beginning
        end     = SEEK_END, ///< From the end
    };

    struct audio_format
    {
        size_t channels        = 2;
        audio_sample_type type = audio_sample_type::i16;
        kfr::fmax samplerate        = 44100;
        bool use_w64           = false;
    };

    struct audio_format_and_length : audio_format
    {        
        constexpr audio_format_and_length();
        constexpr audio_format_and_length(const audio_format& fmt);

        imax length = 0; // in samples
    };

    constexpr size_t audio_sample_sizeof(audio_sample_type type);
    constexpr size_t audio_sample_bit_depth(audio_sample_type type);

    
    /*
    void CMT_ARCH_NAME::deinterleave(SampleType* out[], const SampleType *in, size_t channels, size_t size);
    void CMT_ARCH_NAME::interleave(SampleType* out, const SampleType* in[], size_t channels, size_t size);
    void CMT_ARCH_NAME::convert(SampleType* out, const SampleTYpe* in, size_t size);
    */
   
    struct fraction
    {
        fraction(i64 num = 0, i64 den = 1);
        void normalize();
        
        i64 numerator;
        i64 denominator;

        fraction operator+() const;
        fraction operator-() const;

        //explicit operator bool() const;
        //explicit operator double() const;
        //explicit operator float() const;
        //explicit operator kfr::signed long long() const;
    };    

    template <typename T>
    struct complex
    {        
        constexpr complex()  = default;
        constexpr complex(T re)  : re(re), im(0) {}
        constexpr complex(T re, T im)  : re(re), im(im) {}
        constexpr complex(const complex&)  = default;
        
        constexpr complex& operator=(const complex&)  = default;
        constexpr complex& operator=(complex&&)  = default;
        constexpr const T& real() const  { return re; }
        constexpr const T& imag() const  { return im; }
        constexpr void real(T value)  { re = value; }
        constexpr void imag(T value)  { im = value; }    
    };
}

%inline %{
    /*
    template<typename T> kfr::univector<T> csin(const kfr::univector<kfr::complex<T>> & v) { return kfr::csin(v); }
    template<typename T> kfr::univector<T> ccos(const kfr::univector<kfr::complex<T>> & v) { return kfr::ccos(v); }
    //template<typename T> kfr::univector<T> ctan(const kfr::univector<kfr::complex<T>> & v) { return kfr::ctan(v); }

    template<typename T> kfr::univector<T> csinh(const kfr::univector<kfr::complex<T>> & v) { return kfr::csinh(v); }
    template<typename T> kfr::univector<T> ccosh(const kfr::univector<kfr::complex<T>> & v) { return kfr::ccosh(v); }
    //template<typename T> kfr::univector<T> ctanh(const kfr::univector<kfr::complex<T>> & v) { return kfr::ctanh(v); }

    template<typename T> kfr::univector<T> cabssqr(const kfr::univector<kfr::complex<T>> & v) { return kfr::cabssqr(v); }
    template<typename T> kfr::univector<T> cabs(const kfr::univector<kfr::complex<T>> & v) { return kfr::cabs(v); }
    template<typename T> kfr::univector<T> carg(const kfr::univector<kfr::complex<T>> & v) { return kfr::carg(v); }
    
    template<typename T> kfr::univector<T> clog(const kfr::univector<kfr::complex<T>> & v) { return kfr::clog(v); }
    template<typename T> kfr::univector<T> clog2(const kfr::univector<kfr::complex<T>> & v) { return kfr::clog2(v); }
    template<typename T> kfr::univector<T> clog10(const kfr::univector<kfr::complex<T>> & v) { return kfr::clog10(v); }

    template<typename T> kfr::univector<T> cexp(const kfr::univector<kfr::complex<T>> & v) { return kfr::cexp(v); }
    template<typename T> kfr::univector<T> cexp2(const kfr::univector<kfr::complex<T>> & v) { return kfr::cexp2(v); }
    template<typename T> kfr::univector<T> cexp10(const kfr::univector<kfr::complex<T>> & v) { return kfr::cexp10(v); }

    template<typename T> kfr::univector<T> polar(const kfr::univector<kfr::complex<T>> & v) { return kfr::polar(v); }
    template<typename T> kfr::univector<T> cartesian(const kfr::univector<kfr::complex<T>> & v) { return kfr::cartesian(v); }
    //template<typename T> kfr::univector<T> cabsdup(const kfr::univector<kfr::complex<T>> & v) { return kfr::cabsdup(v); }

    template<typename T> kfr::univector<T> csqrt(const kfr::univector<kfr::complex<T>> & v) { return kfr::csqrt(v); }
    template<typename T> kfr::univector<T> csqr(const kfr::univector<kfr::complex<T>> & v) { return kfr::csqr(v); }
    */
    template<typename T> kfr::complex<T> csin(const kfr::complex<T> & v) { return kfr::csin(v); }
    template<typename T> kfr::complex<T> ccos(const kfr::complex<T> & v) { return kfr::ccos(v); }
    //template<typename T> kfr::univector<T> ctan(const kfr::univector<kfr::complex<T>> & v) { return kfr::ctan(v); }

    template<typename T> kfr::complex<T> csinh(const kfr::complex<T> & v) { return kfr::csinh(v); }
    template<typename T> kfr::complex<T> ccosh(const kfr::complex<T> & v) { return kfr::ccosh(v); }
    //template<typename T> kfr::univector<T> ctanh(const kfr::univector<kfr::complex<T>> & v) { return kfr::ctanh(v); }

    template<typename T> kfr::complex<T> cabssqr(const kfr::complex<T> & v) { return kfr::cabssqr(v); }
    template<typename T> kfr::complex<T> cabs(const kfr::complex<T> & v) { return kfr::cabs(v); }
    template<typename T> kfr::complex<T> carg(const kfr::complex<T> & v) { return kfr::carg(v); }
    
    template<typename T> kfr::complex<T> clog(const kfr::complex<T> & v) { return kfr::clog(v); }
    template<typename T> kfr::complex<T> clog2(const kfr::complex<T> & v) { return kfr::clog2(v); }
    template<typename T> kfr::complex<T> clog10(const kfr::complex<T> & v) { return kfr::clog10(v); }

    template<typename T> kfr::complex<T> cexp(const kfr::complex<T> & v) { return kfr::cexp(v); }
    template<typename T> kfr::complex<T> cexp2(const kfr::complex<T> & v) { return kfr::cexp2(v); }
    template<typename T> kfr::complex<T> cexp10(const kfr::complex<T> & v) { return kfr::cexp10(v); }

    template<typename T> kfr::complex<T> polar(const kfr::complex<T> & v) { return kfr::polar(v); }
    template<typename T> kfr::complex<T> cartesian(const kfr::complex<T> & v) { return kfr::cartesian(v); }
    //template<typename T> kfr::univector<T> cabsdup(const kfr::univector<kfr::complex<T>> & v) { return kfr::cabsdup(v); }

    template<typename T> kfr::complex<T> csqrt(const kfr::complex<T> & v) { return kfr::csqrt(v); }
    template<typename T> kfr::complex<T> csqr(const kfr::complex<T> & v) { return kfr::csqr(v); }

%}




%template(vectorf32)  std::vector<f32>;
%template(vectorf64)  std::vector<f64>;
%template(vectori8)   std::vector<i8>;
%template(vectorui8)  std::vector<u8>;
%template(vectori16)  std::vector<i16>;
%template(vectorui16) std::vector<u16>;
%template(vectori32)  std::vector<i32>;
%template(vectorui32) std::vector<u32>;
%template(vectori64)  std::vector<i64>;
%template(vectorui64)  std::vector<u64>;

%template(complex_vectorf32) std::vector<kfr::complex<float>>;
%template(complex_vectorf64) std::vector<kfr::complex<double>>;

%template(complexf32) kfr::complex<float>;
%template(complexf64) kfr::complex<double>;


namespace kfr
{
    // single channel vector (SampleVector)
    template<typename T> 
    struct univector
    {        
        univector() {}
        univector(size_t s);
        //univector(T* data, size_t size);
        //univector(std::vector<T> & data);
        univector(const univector<T> & u);        
        

        size_t size() const;
        void resize(size_t s);
        T& operator[](size_t i);
        T  operator[](size_t i);

        %extend {
            // lua is 1 based like fortran
            T  __getitem(size_t i) { assert(i > 0) ; return (*$self)[i-1]; }
            void __setitem(size_t i, const T & val) { assert(i > 0);(*$self)[i-1] = val; }

            univector<T> __add__(const univector& b) { return *$self + b; }
            univector<T> __sub__(const univector& b) { return *$self - b; }
            univector<T> __mul__(const univector& b) { return *$self * b; }
            univector<T> __div__(const univector& b) { return *$self / b; }
            univector<T> __unm__() { return -*$self; }
            univector<T> __pow__(const T& b) { return pow(*$self,b); }
            //bool         __eq__(const univector& b) { return (bool)*$self == b; }
            //bool         __lt__(const univector& b) { return *$self < b; }
            //bool         __le__(const univector& b) { return *$self <= b; }

            void fill(const T& val) { for(size_t i = 0; i < $self->size(); i++) (*$self)[i] = val; }
            void print() const { kfr::println(*$self); }
        }
                
        T& at(size_t pos);
        T& front();
        T& back();
        T* data();
            
        univector<T>& operator = (const univector<T>& u);
        
    };    
}
%inline %{

    kfr::univector<SampleType> to_univector(const std::vector<SampleType> & v) {
        kfr::univector<SampleType> r(v.size());
        std::copy(v.begin(),v.end(),r.begin());
        return r;
    }

    std::vector<SampleType> to_vector(const kfr::univector<SampleType> & v) {
        std::vector<SampleType> r(v.size());
        std::copy(v.begin(),v.end(),r.begin());
        return r;
    }

    template<typename X> X rol(X x, X y) { return kfr::rol(x,y); }
    template<typename X> X ror(X x, X y) { return kfr::ror(x,y); }
    template<typename X> X shl(X x, X y) { return kfr::shl(x,y); }
    template<typename X> X shr(X x, X y) { return kfr::rol(x,y); }

    template<typename T> kfr::univector<T> bitwiseand(const kfr::univector<T> & a, const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::bitwiseand(a,b); return r; }
    template<typename T> kfr::univector<T> bitwiseandnot(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::bitwiseandnot(a,b); return r; }
    template<typename T> kfr::univector<T> bitwisenot(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::bitwisenot(a); return r; }
    template<typename T> kfr::univector<T> bitwiseor(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::bitwiseor(a,b); return r; }
    template<typename T> kfr::univector<T> bitwisexor(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::bitwisexor(a,b); return r; }

    template<typename T>    
    kfr::univector<T> linspace(T start, T stop, size_t size, bool endpoint=false,bool trunc=false)
    {   
        kfr::univector<T> r; 
        r = kfr::linspace(start,stop,size,endpoint,trunc); 
        return r; 
    }

    template<typename T>    
    kfr::univector<T> pad(const kfr::univector<T> & in, const T & fill_value = T(0))
    { 
        kfr::univector<T> r; 
        r = kfr::padded(in,fill_value); 
        return r; 
    }

    template<typename T>    
    kfr::univector<T> slice(const kfr::univector<T> & v, size_t start, size_t end=kfr::max_size_t)
    {   
        kfr::univector<T> r;        
        r = v.slice(start,end);
        return r;
    }

    template<typename T>    
    kfr::univector<T> truncate(const kfr::univector<T> & v, size_t size)
    {   
        kfr::univector<T> r; 
        r = v.truncate();
        return r;
    }

    template<typename T>    
    kfr::univector<T> reverse(const kfr::univector<T> & v)
    {   
        kfr::univector<T> r;         
        r = kfr::reverse(v);
        return r;
    }


    template<typename T>    
    T& ringbuf_read(kfr::univector<T> &v,size_t & cursor, T& value) { v.ringbuf_read(cursor,value); return value; }

    template<typename T>    
    void ringbuf_write(kfr::univector<T> &v, size_t & cursor, T& value) { v.ringbuf_write(cursor,value); }
    
    template<typename T> kfr::univector<T> abs(const kfr::univector<T>& v) { return kfr::abs(v); }
    template<typename T> kfr::univector<T> add(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::add(a,b); return r; }
    template<typename T> kfr::univector<T> add(const kfr::univector<T> & a,const T & b) { kfr::univector<T> r; r = kfr::add(a,b); return r; }
    template<typename T> kfr::univector<T> absmax(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::absmax(a,b); return r; }
    template<typename T> kfr::univector<T> absmax(const kfr::univector<T> & a,const T & b) { kfr::univector<T> r; r = kfr::absmax(a,b); return r; }
    template<typename T> kfr::univector<T> absmin(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::absmin(a,b); return r; }
    template<typename T> kfr::univector<T> absmin(const kfr::univector<T> & a,const T & b) { kfr::univector<T> r; r = kfr::absmin(a,b); return r; }    
    template<typename T> kfr::univector<T> clamp(const kfr::univector<T> & a,const kfr::univector<T> & lo, const kfr::univector<T> &hi) { kfr::univector<T> r; r = kfr::clamp(a,lo,hi); return r; }
    template<typename T> kfr::univector<T> clamp(const kfr::univector<T> & a,const T& lo, const T &hi) { kfr::univector<T> r; r = kfr::clamp(a,lo,hi); return r; }
    template<typename T> kfr::univector<T> cube(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::cub(a); return r; }
    template<typename T> kfr::univector<T> div(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::div(a,b); return r; }    
    template<typename T> kfr::univector<T> fmadd(const kfr::univector<T> & a,const kfr::univector<T> & y, const kfr::univector<T> & z) { kfr::univector<T> r; r = kfr::fmadd(a,y,z); return r; }
    template<typename T> kfr::univector<T> fmsub(const kfr::univector<T> & a,const kfr::univector<T> & y, const kfr::univector<T> & z) { kfr::univector<T> r; r = kfr::fmsub(a,y,z); return r; }    
    template<typename T> kfr::univector<T> max(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::max(a,b); return r; }
    template<typename T> kfr::univector<T> max(const kfr::univector<T> & a, const T & b) { kfr::univector<T> r; r = kfr::max(a,b); return r; }
    template<typename T> kfr::univector<T> min(const kfr::univector<T> & a, const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::min(a,b); return r; }
    template<typename T> kfr::univector<T> min(const kfr::univector<T> & a, const T & b) { kfr::univector<T> r; r = kfr::min(a,b); return r; }
    template<typename T> kfr::univector<T> mix(const kfr::univector<T> & a, const T& c, const kfr::univector<T> & y) { kfr::univector<T> r; r = kfr::mix(c,a,y); return r; }
    template<typename T> kfr::univector<T> mixs(const kfr::univector<T> & a, const T& c, const kfr::univector<T> & y) { kfr::univector<T> r; r = kfr::mixs(c,a,y); return r; }
    template<typename T> kfr::univector<T> mul(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::mul(a,b); return r; }
    template<typename T> kfr::univector<T> mul(const kfr::univector<T> & a, const T & b) { kfr::univector<T> r; r = kfr::mul(a,b); return r; }
    template<typename T> kfr::univector<T> neg(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::neg(a); return r; }        
    template<typename T> kfr::univector<T> sqr(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::sqr(a); return r; }
    template<typename T> kfr::univector<T> sqrt(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::sqrt(a); return r; }
    template<typename T> kfr::univector<T> exp(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::exp(a); return r; }
    template<typename T> kfr::univector<T> exp10(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::exp10(a); return r; }
    template<typename T> kfr::univector<T> exp2(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::exp2(a); return r; }
    template<typename T> kfr::univector<T> exp_fmadd(const kfr::univector<T> & a,const kfr::univector<T> & y, const kfr::univector<T> & z) { kfr::univector<T> r; r = kfr::exp_fmadd(a,y,z); return r; }
    template<typename T> kfr::univector<T> log(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::log(a); return r; }
    template<typename T> kfr::univector<T> log10(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::log10(a); return r; }
    template<typename T> kfr::univector<T> log2(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::log2(a); return r; }
    template<typename T> kfr::univector<T> log_fmadd(const kfr::univector<T> & a,const kfr::univector<T> & y, const kfr::univector<T> & z) { kfr::univector<T> r; r = kfr::log_fmadd(a,y,z); return r; }
    template<typename T> kfr::univector<T> logb(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::logb(a); return r; }
    template<typename T> kfr::univector<T> logm(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::logm(a,b); return r; }
    template<typename T> kfr::univector<T> logn(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::logn(a,b); return r; }
    template<typename T> kfr::univector<T> pow(const kfr::univector<T> & a,const T & y) { kfr::univector<T> r; r = kfr::pow(a,y); return r; }
    template<typename T> kfr::univector<T> pow(const kfr::univector<T> & a,const kfr::univector<T> & y) { kfr::univector<T> r; r = kfr::pow(a,y); return r; }
    template<typename T> kfr::univector<T> root(const kfr::univector<T> & a,const kfr::univector<T> & y) { kfr::univector<T> r; r = kfr::root(a,y); return r; }
    template<typename T> kfr::univector<T> floor(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::floor(a); return r; }        
    template<typename T> kfr::univector<T> acos(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::acos(a); return r; }
    template<typename T> kfr::univector<T> asin(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::asin(a); return r; }
    template<typename T> kfr::univector<T> atan(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::atan(a); return r; }
    template<typename T> kfr::univector<T> atan2(const kfr::univector<T> & a,const T & b) { kfr::univector<T> r; r = kfr::atan2(a,b); return r; }
    template<typename T> kfr::univector<T> atan2(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::atan2(a,b); return r; }
    template<typename T> kfr::univector<T> atan2deg(const kfr::univector<T> & a,const T & b) { kfr::univector<T> r; r = kfr::atan2deg(a,b); return r; }
    template<typename T> kfr::univector<T> atan2deg(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::atan2deg(a,b); return r; }
    template<typename T> kfr::univector<T> atandeg(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::atandeg(a); return r; }
    template<typename T> kfr::univector<T> cos(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::cos(a); return r; }
    template<typename T> kfr::univector<T> sin(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::sin(a); return r; }    
    template<typename T> kfr::univector<T> cosdeg(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::cosdeg(a); return r; }        
    template<typename T> kfr::univector<T> sindeg(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::cosdeg(a); return r; }    
    template<typename T> kfr::univector<T> sinc(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::sinc(a); return r; }
    template<typename T> kfr::univector<T> tan(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::tan(a); return r; }        
    template<typename T> kfr::univector<T> cosh(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::cosh(a); return r; }
    template<typename T> kfr::univector<T> coth(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::coth(a); return r; }    
    template<typename T> kfr::univector<T> sinh(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::sinh(a); return r; }    
    template<typename T> kfr::univector<T> tanh(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::tanh(a); return r; }
    template<typename T> kfr::univector<T> gamma(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::gamma(a); return r; }

    template<typename T> T absmaxof(const kfr::univector<T> & a) { return kfr::absmaxof(a); }
    template<typename T> T absminof(const kfr::univector<T> & a) { return kfr::absminof(a); }
    template<typename T> T dot(const kfr::univector<T> & a,const kfr::univector<T> & b) { return kfr::dotproduct(a,b); }
    template<typename T> T maxof(const kfr::univector<T> & a) { return kfr::maxof(a); }
    template<typename T> T minof(const kfr::univector<T> & a) { return kfr::minof(a); }
    template<typename T> T mean(const kfr::univector<T> & a) { return kfr::mean(a); }
    template<typename T> T product(const kfr::univector<T> & a) { return kfr::product(a); }
    template<typename T> T rms(const kfr::univector<T> & a) { return kfr::rms(a); }
    template<typename T> T sum(const kfr::univector<T> & a) { return kfr::sum(a); }
    template<typename T> T sumsqr(const kfr::univector<T> & a) { return kfr::sumsqr(a); }

    // doesn't compile
    //template<typename T>    
    //kfr::univector<T> div(const kfr::univector<T> & a,const T b) { kfr::univector<T> r; r = kfr::div<T>(a,b); return r; }

    template<typename T>    
    kfr::univector<T> ipow(const kfr::univector<T> & v, int base) { kfr::univector<T> r; r = kfr::ipow(v,base); return r; }

    template<typename T>    
    T kcos2x(const T s, const T c) {return kfr::cos2x<SampleType>(s,c); }

    template<typename T>    
    T kcos3x(const T & s, const T & c) {return kfr::cos3x(s,c); }

    template<typename T>    
    T ksin2x(const T & s, const T & c) {return kfr::sin2x(s,c); }

    template<typename T>    
    T ksin3x(const T & s, const T & c) {return kfr::sin3x(s,c); }

    template<typename T>    
    kfr::univector<T> cossin(const kfr::univector<T> & v) { kfr::univector<T> r; r = kfr::cossin(v); return r; }

    template<typename T>    
    kfr::univector<T> sincos(const kfr::univector<T> & v) { kfr::univector<T> r; r = kfr::sincos(v); return r; }

    template<typename T>    
    T kreciprocal(const T & v) { return kfr::reciprocal(v); }

    template<typename T>    
    T rem(const T v,const T b) { return kfr::rem(v,b); }    

    template<typename T>    
    T satadd(const T v,const T y) { return kfr::satadd(v,y); }

    template<typename T>    
    T satsub(const T v,const T  y) { return kfr::satsub(v,y); }

    //? dont know how to make these work yet.
    template<typename T>    
    kfr::univector<T> fastcos(const kfr::univector<T> & v) { kfr::univector<T> r; r = kfr::fastcos(v); return r; }

    template<typename T>    
    kfr::univector<T> fastcosdeg(const kfr::univector<T> & v) { kfr::univector<T> r; r = kfr::fastcosdeg(v); return r; }

    template<typename T>    
    kfr::univector<T> fastsin(const kfr::univector<T> & v) { kfr::univector<T> r; r = kfr::fastsin(v); return r; }

    template<typename T>    
    kfr::univector<T> fastsindeg(const kfr::univector<T> & v) { kfr::univector<T> r; r = kfr::fastsindeg(v); return r; }        

    template<typename T>    
    kfr::univector<T> coshsinh(const kfr::univector<T> & v) { kfr::univector<T> r; r = kfr::coshsinh(v); return r; }

    template<typename T>    
    kfr::univector<T> sinhcosh(const kfr::univector<T> & v) { kfr::univector<T> r; r = kfr::sinhcosh(v); return r; }

    template<typename T>    
    kfr::univector<T> cossindeg(const kfr::univector<T> & v) { kfr::univector<T> r; r = kfr::cossindeg(v); return r; }    

    template<typename T>    
    kfr::univector<T> sincosdeg(const kfr::univector<T> & v) { kfr::univector<T> r; r = kfr::sincosdeg(v); return r; }    

    // I dont understand the kfr random at all yet
    template<typename T>    
    kfr::univector<T> random(size_t s) 
    {
        std::default_random_engine generator;
        std::uniform_real_distribution<T> distrib(0.0,1.0);
        kfr::univector<T> r(s);    
        for(size_t i = 0; i < s; i++)
            r[i] = distrib(generator);
        return r;
    }   

    template<typename T>    
    kfr::univector<T> random(size_t s, T min, T max) 
    {
        std::default_random_engine generator;
        std::uniform_real_distribution<T> distrib(min,max);
        kfr::univector<T> r(s);    
        for(size_t i = 0; i < s; i++)
            r[i] = distrib(generator);
        return r;
    }  

    template<typename T>
    void plot_save(const kfr::univector<T> & v, const std::string& name="", const std::string& options="") {
            kfr::plot_save(name,v,options);
        }
    template<typename T>    
    void plot_show(const kfr::univector<T> & v, const std::string& name="", const std::string&  options="") {
        kfr::plot_show(name,v,options);
    }

    //template<typename T> kfr::univector<T> make_univec(const T * data, size_t s) { return kfr::univector<T>(kfr::make_univector<T>(data,s));  }


%}

%template (univector) kfr::univector<SampleType>;

%template(csin) csin<SampleType>;
%template(ccos) ccos<SampleType>;
%template(csinh) csinh<SampleType>;
%template(ccosh) ccosh<SampleType>;
%template(cabssqr) cabssqr<SampleType>;
%template(cabs) cabs<SampleType>;
%template(carg) carg<SampleType>;
%template(clog) clog<SampleType>;
%template(clog10) clog10<SampleType>;
%template(clog2) clog2<SampleType>;
%template(cexp) cexp<SampleType>;
%template(cexp2) cexp2<SampleType>;
%template(cexp10) cexp10<SampleType>;
%template(cpolar) polar<SampleType>;
%template(ccartesian) cartesian<SampleType>;
%template(csqrt) csqrt<SampleType>;
%template(csqr) csqr<SampleType>;

%template(abs) abs<SampleType>;
%template(add) add<SampleType>;
%template(absmax) absmax<SampleType>;
%template(absmin) absmin<SampleType>;
%template(clamp) clamp<SampleType>;
%template(cube) cube<SampleType>;
%template(div) div<SampleType>;
%template(fmadd) fmadd<SampleType>;
%template(fmsub) fmsub<SampleType>;
%template(max) max<SampleType>;
%template(min) min<SampleType>;
%template(mix) mix<SampleType>;
%template(mixs) mixs<SampleType>;
%template(mul) mul<SampleType>;
%template(neg) neg<SampleType>;
%template(sqr) sqr<SampleType>;
%template(sqrt) sqrt<SampleType>;
%template(exp) exp<SampleType>;
%template(exp10) exp10<SampleType>;
%template(exp2) exp2<SampleType>;
%template(exp_fmadd) exp_fmadd<SampleType>;
%template(log) log<SampleType>;
%template(log10) log10<SampleType>;
%template(log2) log2<SampleType>;
%template(log_fmadd) log_fmadd<SampleType>;
%template(logb) logb<SampleType>;
%template(logm) logm<SampleType>;
%template(logn) logn<SampleType>;
%template(pow) pow<SampleType>;
%template(root) root<SampleType>;
%template(floor) floor<SampleType>;
%template(acos) acos<SampleType>;
%template(asin) asin<SampleType>;
%template(atan) atan<SampleType>;
%template(atan2) atan2<SampleType>;
%template(cos) cos<SampleType>;
%template(sin) sin<SampleType>;
%template(tan) tan<SampleType>;
%template(cosh) cosh<SampleType>;
%template(coth) coth<SampleType>;
%template(sinh) sinh<SampleType>;
%template(tanh) tanh<SampleType>;
%template(atan2deg) atan2deg<SampleType>;
%template(cosdeg) cosdeg<SampleType>;
%template(sindeg) sindeg<SampleType>;
%template(sinc) sinc<SampleType>;
%template(gamma) gamma<SampleType>;
%template(absmaxo) absmaxof<SampleType>;
%template(dot) dot<SampleType>;
%template(maxo) maxof<SampleType>;
%template(mino) minof<SampleType>;
%template(mean) mean<SampleType>;
%template(prdocut) product<SampleType>;
%template(rms) rms<SampleType>;
%template(sum) sum<SampleType>;
%template(sumsqr) sumsqr<SampleType>;
%template(ipow) ipow<SampleType>;
%template(cos2x) kcos2x<SampleType>;
%template(sin2x) ksin2x<SampleType>;
%template(cos3x) kcos3x<SampleType>;
%template(sin3x) ksin3x<SampleType>;
%template(reciprocal) kreciprocal<SampleType>;

%template(linspace) linspace<SampleType>;
%template(pad) pad<SampleType>;
%template(slice) slice<SampleType>;
%template(truncate) truncate<SampleType>;
%template(reverse) reverse<SampleType>;

// univector
%template(slice) slice<SampleType>;
%template(truncate) truncate<SampleType>;
%template(reverse) reverse<SampleType>;
%template(ringbuf_read) ringbuf_read<SampleType>;
%template(ringbuf_write) ringbuf_write<SampleType>;
%template (random) random<SampleType>;

// Plot
%template (plot_save) plot_save<SampleType>;
%template (plot_show) plot_show<SampleType>;

//%template(make_univector) make_univec<SampleType>;

%template(acosh) acosh<SampleType>;
%template(asinh) asinh<SampleType>;
%template(atanh) atanh<SampleType>;
%template(cbrt) cbrt<SampleType>;
%template(ceil) ceil<SampleType>;
%template(copysign) copysign<SampleType>;
%template(er) erf<SampleType>;
%template(erfc) erfc<SampleType>;
%template(expm1) expm1<SampleType>;
%template(fdim) fdim<SampleType>;
%template(fma) fma<SampleType>;
//%template(fmax) fmax<SampleType>;
//%template(fmin) fmin<SampleType>;
%template(fmod) fmod<SampleType>;
%template(fpclassify) fpclassify<SampleType>;
%template(hypot) hypot<SampleType>;
%template(ilogb) ilogb<SampleType>;
%template(isfinite) isfinite<SampleType>;
%template(isgreater) isgreater<SampleType>;
%template(isgreaterequal) isgreaterequal<SampleType>;
%template(isin) isinf<SampleType>;
%template(isless) isless<SampleType>;
%template(islessequal) islessequal<SampleType>;
%template(isnan) isnan<SampleType>;
%template(isnormal) isnormal<SampleType>;
%template(isunordered) isunordered<SampleType>;
%template(ldexp) ldexp<SampleType>;
%template(lgamma) lgamma<SampleType>;
%template(llrint) llrint<SampleType>;
%template(llround) llround<SampleType>;
%template(log1p) log1p<SampleType>;
%template(lrint) lrint<SampleType>;
%template(lround) lround<SampleType>;
%template(nan) nan<SampleType>;
%template(nanf) nanf<SampleType>;
%template(nanl) nanl<SampleType>;
%template(nearbyint) nearbyint<SampleType>;
%template(nextafter) nextafter<SampleType>;
%template(nexttoward) nexttoward<SampleType>;
%template(remainder) remainder<SampleType>;
%template(rint) rint<SampleType>;
%template(round) round<SampleType>;
%template(scalbln) scalbln<SampleType>;
%template(scalbn) scalbn<SampleType>;
%template(square) square<SampleType>;
%template(tgamma) tgamma<SampleType>;
%template(trunc) trunc<SampleType>;

/*
%template (fastcos) fastcos<f32>;
%template (fastcosdeg) fastcosdeg<f32>;
%template (fastsin) fastsin<f32>;
%template (fastsindeg) fastsindeg<f32>;
%template (coshsing) coshsinh<f32>;
%template (sinhcosh) sinhcosh<f32>;
%template (cossindeg) cossindeg<f32>;
%template (sincosdeg) fastsindeg<f32>;
*/
/*
%template(sataddi64) satadd<i64>;
%template(roli32) rol<i32>;
%template(ro4i32) ror<i32>;
%template(shli32) shl<i32>;
%template(shri32) shr<i32>;
%template(remi32)    rem<i32>;sh 
%template(sataddi32) satadd<i32>;
%template(satsubi32) satsub<i32>;
%template(bitwiseandi32) bitwiseand<i32>;
%template(bitwiseori32) bitwiseor<i32>;
%template(bitwisexori32) bitwisexor<i32>;
%template(bitwiseandnoti32) bitwiseandnot<i32>;
%template(bitwisenoti32) bitwisenot<i32>;
*/

%inline %{

    template<typename T>
    struct SampleVector
    {
        std::vector<kfr::univector<T>> samples;
        size_t                         channels;

        SampleVector(size_t channels) {
            samples.resize(channels);
            this->channels = channels;
        }
        
        T& operator()(size_t ch, size_t i) { return samples[ch][i]; }
        
        size_t num_channels() const { return channels; }
        size_t size() const { return samples[0].size(); }
        
        kfr::univector<T> get_channel(size_t channel) { return samples[channel]; }
        void set_channel(size_t channel, kfr::univector<T> & v) { samples[channel] = v; }

        kfr::univector<T> __getitem(size_t i ) { return samples[i]; }
        void __setitem(size_t i, kfr::univector<T> & v) { samples[i] = v; }

    };

    /*
    template<typename T>
    SampleVector<T> deinterleave(size_t channels, kfr::univector<T> & v) {
        SampleVector<T> r(channels);        
        for(size_t i = 0; i < channels; i++) {
            r.samples[i].resize(v.size()/channels)
            for(size_t j = i; j < v.size(); j += channels)
                r[channels][j] = v[j];
        }
        return r;
    }
    template<typename T>
    void interleave(SampleVector<T> & samples, kfr::univector<T> & out) {
        out.resize(samples.channels * samples[0].size());        
        for(size_t i = 0; i < samples.channels; i++)            
            for(size_t j = i; j < samples[i].size(); i+=samples.channels)
                out[j*channels + i] = samples[i][j];
        }
    */
    template<typename T>
    void copy(kfr::univector<T> & dst, std::vector<T> & src) {
        std::copy(src.begin(),src.end(),dst.begin());
    }
    template<typename T>
    void copy(std::vector<T> & dst, kfr::univector<T> & src) {
        std::copy(src.begin(),src.end(),dst.begin());
    }

%}



%inline %{

        template <typename T> T f_note_to_hertz(const T& input) {
            return kfr::note_to_hertz<T>(input);
        }    
        template <typename T> T f_hertz_to_note(const T& input) {
            return kfr::hertz_to_note<T>(input);
        }    
        template <typename T> T f_amp_to_dB(const T& input) {
            return kfr::amp_to_dB<T>(input);
        }    
        template <typename T> T f_dB_to_amp(const T& input) {
            return kfr::dB_to_amp<T>(input);
        }    
        template <typename T> T f_power_to_dB(const T& input) {
            return kfr::power_to_dB<T>(input);
        }    
        template <typename T> T f_dB_to_power(const T& input) {
            return kfr::dB_to_power<T>(input);
        }    
        
        /*
        template<typename T> kfr::complex<T> goertzal(kfr::complex<T> & result, T  omega) {
            kfr::complex<T> r(result);
            kfr::goertzal(r,omega);
            return r;
        }
        */
        

        template <typename T> T waveshaper_hardclip(T & input, double clip_level) 
        {            
            return kfr::waveshaper_hardclip(input,clip_level);
        }
        template <typename T> kfr::univector<T> waveshaper_hardclip(kfr::univector<T> & input, double clip_level) 
        {            
            kfr::univector r(input.size());
            for(size_t i = 0; i < input.size(); i++)
                r[i] = kfr::waveshaper_hardclip(input[i],clip_level);
            return r;
        }

        template <typename T> T waveshaper_tanh(T & input, double sat) 
        {            
            return kfr::waveshaper_tanh(input,sat);
        }
        template <typename T> kfr::univector<T> waveshaper_tanh(kfr::univector<T> & input, double sat) 
        {            
            kfr::univector r(input.size());
            for(size_t i = 0; i < input.size(); i++)
                r[i] = kfr::waveshaper_tanh(input[i],sat);
            return r;
        }

        template <typename T> T waveshaper_saturate_I(T & input, double sat) 
        {            
            return kfr::waveshaper_saturate_I(input,sat);
        }
        template <typename T> kfr::univector<T> waveshaper_saturate_I(kfr::univector<T> & input, double sat) 
        {            
            kfr::univector r(input.size());
            for(size_t i = 0; i < input.size(); i++)
                r[i] = kfr::waveshaper_saturate_I(input[i],sat);
            return r;
        }

        template <typename T> T waveshaper_saturate_II(T & input, double sat) 
        {            
            return kfr::waveshaper_saturate_II(input, sat);
        }
        template <typename T> kfr::univector<T> waveshaper_saturate_II(kfr::univector<T> & input, double sat) 
        {            
            kfr::univector r(input.size());
            for(size_t i = 0; i < input.size(); i++)
                r[i] = kfr::waveshaper_saturate_II(input[i],sat);
            return r;
        }
        /*
        template <typename T> T waveshaper_poly(T & input) 
        {            
            return kfr::waveshaper_poly(input);
        }
        template <typename T> kfr::univector<T> waveshaper_poly(kfr::univector<T> & input) 
        {            
            kfr::univector r(input.size());
            for(size_t i = 0; i < input.size(); i++)
                r[i] = kfr::waveshaper_poly(input[i]);
            return r;
        }
        */

%}

%template(waveshaper_hardclip) waveshaper_hardclip<SampleType>;
%template(waveshaper_tanh) waveshaper_tanh<SampleType>;
%template(waveshaper_saturateI) waveshaper_saturate_I<SampleType>;
%template(waveshaper_saturateII) waveshaper_saturate_II<SampleType>;
//%template(waveshaper_poly) waveshaper_poly<SampleType>;


%template(note_to_hertz) f_note_to_hertz<SampleType>;
%template(hertz_to_note) f_hertz_to_note<SampleType>;
%template(amp_to_dB) f_amp_to_dB<SampleType>;
%template(dB_to_amp) f_dB_to_amp<SampleType>;
%template(power_to_dB) f_power_to_dB<SampleType>;
%template(dB_to_power) f_dB_to_power<SampleType>;

%template(dcremove) DSP::dcremove<SampleType>;

%template(window_hann) DSP::make_window_hann<SampleType>;
%template(window_hamming) DSP::make_window_hamming<SampleType>;
%template(window_blackman) DSP::make_window_blackman<SampleType>;
%template(window_blackman_harris) DSP::make_window_blackman_harris<SampleType>;
%template(window_gaussian) DSP::make_window_gaussian<SampleType>;
%template(window_triangular) DSP::make_window_triangular<SampleType>;
%template(window_bartlett) DSP::make_window_bartlett<SampleType>;
%template(window_cosine) DSP::make_window_cosine<SampleType>;
%template(window_bartlett_hann) DSP::make_window_bartlett_hann<SampleType>;
%template(window_bohman) DSP::make_window_bohman<SampleType>;
%template(window_lanczos) DSP::make_window_lanczos<SampleType>;
%template(window_flattop) DSP::make_window_flattop<SampleType>;
%template(window_kaiser) DSP::make_window_kaiser<SampleType>;

%template(window_hann_ptr) DSP::make_window_hann_ptr<SampleType>;
%template(window_hamming_ptr) DSP::make_window_hamming_ptr<SampleType>;
%template(window_blackman_ptr) DSP::make_window_blackman_ptr<SampleType>;
%template(window_blackman_harris_ptr) DSP::make_window_blackman_harris_ptr<SampleType>;
%template(window_gaussian_ptr) DSP::make_window_gaussian_ptr<SampleType>;
%template(window_triangular_ptr) DSP::make_window_triangular_ptr<SampleType>;
%template(window_bartlett_ptr) DSP::make_window_bartlett_ptr<SampleType>;
%template(window_cosine_ptr) DSP::make_window_cosine_ptr<SampleType>;
%template(window_bartlett_hann_ptr) DSP::make_window_bartlett_hann_ptr<SampleType>;
%template(window_bohman_ptr) DSP::make_window_bohman_ptr<SampleType>;
%template(window_lanczos_ptr) DSP::make_window_lanczos_ptr<SampleType>;
%template(window_flattop_ptr) DSP::make_window_flattop_ptr<SampleType>;
%template(window_kaiser_ptr) DSP::make_window_kaiser_ptr<SampleType>;

%template(energy_to_loudness) DSP::energy_to_loudness<SampleType>;
%template(loudness_to_energy) DSP::loudness_to_energy<SampleType>;

%template(float_wavreader) DSP::WavReader<SampleType>;
%template(float_wavwriter) DSP::WavWriter<SampleType>;
%template(float_mp3_reader) DSP::MP3Reader<SampleType>;
%template(float_flac_reader) DSP::FlacReader<SampleType>;

%template (sinewave) DSP::sinewave<SampleType>;
%template (squarewave) DSP::squarewave<SampleType>;
%template (trianglewave) DSP::trianglewave<SampleType>;
%template (sawtoothwave) DSP::sawtoothwave<SampleType>;

%template (generate_sine) DSP::generate_sin<SampleType>;
%template (generate_linear) DSP::generate_linear<SampleType>;
%template (generate_exp) DSP::generate_exp<SampleType>;
%template (generate_exp2) DSP::generate_exp2<SampleType>;
%template (generate_cossin) DSP::generate_cossin<SampleType>;


%template (resample) DSP::resample<SampleType>;
%template (convert_sample) DSP::convert_sample<SampleType>;
//%template (amp_to_dB) DSP::amp_to_dB<SampleType>;

%template(interleave) DSP::do_interleave<SampleType>;
%template(deinterleave) DSP::do_deinterleave<SampleType>; 

%template(load_wav) DSP::load_wav<SampleType>;
%template(save_wav) DSP::write_wav<SampleType>;
%template(load_mp3) DSP::load_mp3<SampleType>;
%template(load_flac) DSP::load_flac<SampleType>;
