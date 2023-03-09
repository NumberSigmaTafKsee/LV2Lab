#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>
#include <complex>
#include <vector>
#include <new>
#include <chrono>
#include <random>
#include <cassert>

#include "SoundObject.hpp"
#include "stdsamples_allocator.hpp"
#include "stdsamples_sndfile.hpp"

#include "stdsamples_vector.hpp"
#include "stdsamples_matrix.hpp"
#include "stdsamples_complex_vector.hpp"
#include "stdsamples_complex_matrix.hpp"


namespace AudioDSP
{    
    struct wav_data {
        int64_t frames;
        size_t  size;
        int     samplerate;
        int     channels;
        int     format;
        int     sections;
    };
    
    void load_wave(const std::string& file, wav_data & info, sample_vector<float> & wav)
    {
        SndFileReaderFloat r(file.c_str());
        wav.resize(r.size());
        r.read(r.size(),wav.data());
        info.frames = r.frames();
        info.size   = r.size();
        info.samplerate = r.samplerate();
        info.channels = r.channels();
        info.format   = r.format();
        info.sections = r.sections();        
    }    
    void save_wave(const std::string& file, sample_vector<float> & samples, wav_data& info)
    {
        SndFileWriterFloat r(file.c_str(),info.format,info.channels,info.samplerate);
        r.write(samples.size(),samples.data());
    }
    void load_wave(const std::string& file, wav_data & info,sample_vector<double> & wav)
    {
        SndFileReaderDouble r(file.c_str());
        wav.resize(r.size());
        r.read(r.size(),wav.data());
        info.frames = r.frames();
        info.size   = r.size();
        info.samplerate = r.samplerate();
        info.channels = r.channels();
        info.format   = r.format();
        info.sections = r.sections();        
    }
    void save_wave(const std::string& file, sample_vector<double> & samples, wav_data& info)
    {
        SndFileWriterDouble r(file.c_str(),info.format,info.channels,info.samplerate);
        r.write(samples.size(),samples.data());
    }

    template<typename T>
    std::ostream& operator << (std::ostream & o, const sample_matrix<T> & m )
    {
        for(size_t i = 0; i < m.rows(); i++)
        {
            for(size_t j = 0; j < m.cols(); j++)
                o << m(i,j) << ",";
            o << std::endl;
        }
        return o;
    }
    
    template<typename T>
	void print_vector(const sample_vector<T> & v) {
		for(size_t i = 0; i < v.size(); i++) std::cout << v[i] << ",";
		std::cout << std::endl;
	}
	
    template<typename T>
    T get_stride(size_t ch, size_t num_channels, size_t pos, sample_vector<T> & samples)
    {
        return samples[pos*num_channels + ch];
    }
        
    template<typename T>
    void set_stride(size_t ch, size_t num_channels, size_t pos, sample_vector<T> & samples, T sample)
    {
        samples[pos*num_channels + ch] = sample;
    }

    template<typename T>
    sample_vector<T> get_left_channel(const sample_vector<T> & in) {
        sample_vector<T> r(in.size()/2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < in.size(); i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    sample_vector<T> get_right_channel(const sample_vector<T> & in) {
        sample_vector<T> r(in.size()/2);
        size_t x = 0;
        const T * input = in.data();
        T * output = r.data();
        #pragma omp simd aligned(input,output)
        for(size_t i = 1; i < in.size(); i+=2) output[x++] = input[i];
        return r;
    }
    template<typename T>
    sample_vector<T> get_channel(size_t ch, const sample_vector<T> & in) {
        sample_vector<T> r(in.size()/2);
        size_t x = 0;
        const T * input = in.data();
        T * output = r.data();
        #pragma omp simd aligned(input,output)
        for(size_t i = ch; i < in.size(); i+=2) output[x++] = input[i];
        return r;
    }
    template<typename T>
    void set_left_channel(const sample_vector<T> & left, sample_vector<T> & out) {
        size_t x = 0;
        const T * input = left.data();
        T * output = out.data();
        #pragma omp simd aligned(input,output)        
        for(size_t i = 0; i < out.size(); i+=2) output[i] = input[x++];
    }
    template<typename T>
    void set_right_channel(const sample_vector<T> & right, sample_vector<T> & out) {
        size_t x = 0;
        const T * input = right.data();
        T * output = out.data();
        #pragma omp simd aligned(input,output)
        for(size_t i = 1; i < out.size(); i+=2) output[i] = input[x++];
    }
    template<typename T>
    void set_channel(size_t ch, const sample_vector<T> & in, sample_vector<T> & out) {
        size_t x = 0;
        const T * input = in.data();
        T * output = out.data();
        #pragma omp simd aligned(input,output)
        for(size_t i = ch; i < out.size(); i+=2) output[i] = input[x++];
    }
    template<typename T>
    sample_vector<T> interleave(size_t n, size_t channels, const sample_vector<sample_vector<T>> & in) {
        sample_vector<T> r(n*channels);
        
        T * output = r.data();        
        for(size_t i = 0; i < channels; i++) {
			const T * input  = in[i].data();
			#pragma omp simd aligned(input,output)
            for(size_t j = 0; j < n; j++)
                output[j*channels + i] = input[j];
		}
        return r;
    }
    template<typename T>
    sample_vector<T> interleave(size_t n, size_t channels, const sample_vector<T*> & in) {
        sample_vector<T> r(n*channels);
        
        T * output = r.data();        
        for(size_t i = 0; i < channels; i++)
        {
			T * const input = in[i];
			#pragma omp simd aligned(input,output)
            for(size_t j = 0; j < n; j++)
                output[j*channels + i] = input[j];
		}
        return r;
    }
    template<typename T>
    sample_vector<sample_vector<T>> deinterleave(size_t n, size_t channels, const sample_vector<T> & in) {
        sample_vector<sample_vector<T>> r(n);
        const T * input  = in.data();        
        for(size_t i = 0; i < channels; i++)
        {			
            r[i].resize(n);            
            T * output = r[i].data();                    
            #pragma omp simd aligned(input,output)
            for(size_t j = 0; j < n; j++)
                r[j] = in[j*channels + i];
        }
        return r;
    }
    template<typename T>
    void fill_left_channel(sample_vector<T> & out, const T v) {
        size_t x = 0;
        T * output = out.data();
        #pragma omp simd aligned(output)
        for(size_t i = 0; i < out.size(); i+=2) output[i] = v;
    }
    template<typename T>
    void fill_right_channel(sample_vector<T> & out, const T v) {
        size_t x = 0;
        T * output = out.data();
        #pragma omp simd aligned(output)
        for(size_t i = 1; i < out.size(); i+=2) output[i] = v;
    }
    template<typename T>
    void fill_channel(size_t ch, sample_vector<T> & out, const T v) {
        size_t x = 0;
        T * output = out.data();
        #pragma omp simd aligned(output)
        for(size_t i = ch; i < out.size(); i+=2) output[i] = v;
    }


    template<typename T>
    sample_vector<T> interleave(size_t n, size_t channels, const T** & in) {
        sample_vector<T> r(n*channels);
        T * output = r.data();
        for(size_t i = 0; i < channels; i++)        
			#pragma omp simd aligned(output,in)
            for(size_t j = 0; j < n; j++)
                output[j*channels + i] = in[i][j];
        return r;
    }
    template<typename T>
    sample_vector<sample_vector<T>> deinterleave(size_t n, size_t channels, const T* & in) {
        sample_vector<sample_vector<T>> r(n);        
        for(size_t i = 0; i < channels; i++)
        {
            r[i].resize(n);
            T * output = r[i].data();
            #pragma omp simd aligned(output,in)
            for(size_t j = 0; j < n; j++)
                output[j] = in[j*channels + i];
        }
        return r;
    }

    template<typename T>
    bool equal_vector (sample_vector<T> & a, sample_vector<T> & b) {        
        return std::equal(a.begin(),a.end(),b.begin());
    }

    template<typename T>
    void copy_vector(sample_vector<T> & dst, sample_vector<T> & src) {        
        std::copy(src.begin(),src.end(),dst.begin());
    }
    template<typename T>
    void copy_vector(sample_vector<T> & dst, size_t n, T * src) {        
        std::copy(&src[0],&src[n-1],dst.begin());
    }
    template<typename T>
    sample_vector<T> slice_vector(size_t start, size_t end, sample_vector<T> & src) {
        sample_vector<T> r(end-start);        
        std::copy(src.begin()+start,src.begin()+end,r.begin());
        return r;
    }

    template<typename T>
    sample_vector<T> slice_buffer(size_t start, size_t end, T * ptr) {
        sample_vector<T> r(end-start);
        std::copy(&ptr[start],&ptr[end],r.begin());
        return r;
    }
	
	template<typename T>
    T get_stride(size_t ch, size_t num_channels, size_t pos, T * samples)
    {
        return samples[pos*num_channels + ch];
    }
    
    template<typename T>
    void set_stride(size_t ch, size_t num_channels, size_t pos, T * samples, T sample)
    {
        samples[pos*num_channels + ch] = sample;
    }

    template<typename T>
    sample_vector<T> get_left_channel(size_t n, const T* in) {
        sample_vector<T> r(n/2);
        size_t x = 0;
        T * output = r.data();
        #pragma omp simd aligned(in,output)
        for(size_t i = 0; i < n; i+=2) output[x++] = in[i];
        return r;
    }
    template<typename T>
    sample_vector<T> get_right_channel(size_t n, const T* & in) {
        sample_vector<T> r(n/2);
        size_t x = 0;
        T * output = r.data();
        #pragma omp simd aligned(in,output)
        for(size_t i = 1; i < n; i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    sample_vector<T> get_channel(size_t ch, size_t n, T* in) {
        sample_vector<T> r(n/2);
        size_t x = 0;
        T * output = r.data();
        #pragma omp simd aligned(in,output)
        for(size_t i = ch; i < n; i+=2) r[x++] = in[i];
        return r;
    }

    template<typename T>
    void set_left_channel(size_t n, const T* left, T* out) {
        size_t x = 0;
        #pragma omp simd aligned(left,out)        
        for(size_t i = 0; i < n; i+=2) out[i] = left[x++];
    }
    template<typename T>
    void set_right_channel(size_t n, const T* right, T* out) {
        size_t x = 0;
        #pragma omp simd aligned(right,out)
        for(size_t i = 1; i < n; i+=2) out[i] = right[x++];
    }
    template<typename T>
    void set_channel(size_t ch, size_t n, const T* in, T* out) {
        size_t x = 0;
        #pragma omp simd aligned(in,out)
        for(size_t i = ch; i < n; i+=2) out[i] = in[x++];
    }

    template<typename T>
    void fill_left_channel(size_t n, const T v, T* out) {
        size_t x = 0;
        #pragma omp simd aligned(out)
        for(size_t i = 0; i < n; i+=2) out[i] = v;
    }
    template<typename T>
    void fill_right_channel(size_t n, const T v, T* out) {
        size_t x = 0;
        #pragma omp simd aligned(out)
        for(size_t i = 1; i < n; i+=2) out[i] = v;
    }
    template<typename T>
    void fill_channel(size_t ch, size_t n, const T v, T* out) {
        size_t x = 0;
        #pragma omp simd aligned(out)
        for(size_t i = ch; i < n; i+=2) out[i] = v;
    }

    template<typename T>
    void copy_buffer(size_t n, T * dst, T * src) {                        
        memcpy(dst,src,n*sizeof(T));
    }
    template<typename T>
    void split_stereo(size_t n, const T* input, T * left, T * right)
    {
        size_t x=0;
        #pragma omp simd aligned(input,left,right)
        for(size_t i = 0; i < n; i+=2)
        {
            left[x] = input[i];
            right[x++] = input[i+1];
        }
    }
         
    template<typename T>
    void split_stereo(const sample_vector<T> & input, sample_vector<T> & left, sample_vector<T> & right) {
        size_t x = input.size();
        left.resize(x/2);
        right.resize(x/2);
        split_stereo(x,input.data(),left.data(),right.data());
    }

    template<typename T>
    void swap(sample_vector<T> & left, sample_vector<T> & right) {
        std::swap(left,right);
    }

    template<typename T>
    bool contains(const sample_vector<T> & v, const T val) {
        return std::find(v.begin(),v.end(),val) != v.end();
    }

    template<typename T>
    sample_vector<T> mix(const sample_vector<T> & a, const sample_vector<T> & b, T m = 0.5)
    {
        assert(a.size() == b.size());
        sample_vector<T> r(a.size());
        T max = -99999;
        T * output = r.data();
        T * ina    = a.data();
        T * inb    = b.data();
        #pragma omp simd aligned(output,ina,inb)
        for(size_t i = 0; i < r.size(); i++) 
        {            
            output[i] = m*(ina[i]+inb[i]);
        }                
        return r;
    }
    template<typename T>
    sample_vector<T> normalize(const sample_vector<T> & a) {
        sample_vector<T> r(a);        
        auto max = std::max_element(r.begin(),r.end());
        T * output = r.data();
        if(*max > 0) 
            #pragma omp simd aligned(output)
            for(size_t i = 0; i < r.size(); i++) output[i] /= *max;
        return r;
    }
    template<class A, class B>
    sample_vector<B> convert(const sample_vector<A> & v) {
        sample_vector<B> r(v.size());
        B * output = r.data();
        A * input  = v.data();
        #pragma omp simd aligned(output,input)
        for(size_t i = 0; i < v.size(); i++) output[i] = B(input[i]);
        return r;
    }
    template<class T>
    sample_vector<T> kernel(const sample_vector<T> & v, T (*f)(T value)) {
        sample_vector<T> r(v.size());
        T * output = r.data();
        T * input  = v.data();
        #pragma omp simd aligned(input,output)
        for(size_t i = 0; i < v.size(); i++) output[i] = f(input[i]);
        return r;
    }
    template<class T>
    sample_vector<T> kernel(const sample_vector<T> & v, std::function<T (T)> func) {
        sample_vector<T> r(v.size());
        T * output = r.data();
        T * input  = v.data();
        #pragma omp simd aligned(input,output)
        for(size_t i = 0; i < v.size(); i++) output[i] = f(input[i]);
        return r;
    }
    template<class T>
    void inplace_add(const sample_vector<T> & a, sample_vector<T> & r, std::function<T (T)> func) {        
        T * output = r.data();
        T * input  = a.data();
        #pragma omp simd aligned(input,output)
        for(size_t i = 0; i < a.size(); i++) output[i] += func(input[i]);        
    }
    template<class T>
    void inplace_sub(const sample_vector<T> & a, sample_vector<T> & r, std::function<T (T)> func) {        
        T * output = r.data();
        T * input  = a.data();
        #pragma omp simd aligned(input,output)
        for(size_t i = 0; i < a.size(); i++) output[i] -= func(input[i]);        
    }
    template<class T>
    void inplace_mul(const sample_vector<T> & a, sample_vector<T> & r, std::function<T (T)> func) {        
        T * output = r.data();
        T * input  = a.data();
        #pragma omp simd aligned(input,output)
        for(size_t i = 0; i < a.size(); i++) output[i] *= func(input[i]);        
    }
    template<class T>
    void inplace_div(const sample_vector<T> & a, sample_vector<T> & r, std::function<T (T)> func) {        
        T * output = r.data();
        T * input  = a.data();
        #pragma omp simd aligned(input,output)
        for(size_t i = 0; i < a.size(); i++) output[i] /= func(input[i]);
    }


    template<class T>
    void fill(sample_vector<T> & in, T x)
    {
        std::fill(in.begin(),in.end(),x);        
    }
    template<class T>
    void zeros(sample_vector<T> & in)
    {
        fill(in,T(0));
    }
    template<class T>
    void ones(sample_vector<T> & in)
    {
        fill(in,T(1));
    }
    
    template<typename T>
    T min(sample_vector<T> & v) {
		T x = *std::min_element(v.begin(),v.end());
		return x;
	}
	template<typename T>
    T max(sample_vector<T> & v) {
		T x = *std::max_element(v.begin(),v.end());
		return x;
	}
	template<typename T>
    T mean(sample_vector<T> & v) {
		T x = 0;		
        T * input  = v.data();
        #pragma omp simd aligned(input)
		for(size_t i = 0; i < v.size(); i++) x += input[i];
		return x/v.size();
	}
    
    template<typename T>
    sample_vector<T> linspace(int start, int end)
    {        
        T delta = start >= end? -1.0 : 1.0;
        int m = std::floor((end-start)/delta);
        if(start > end) m++;
        sample_vector<T> r(m);
        T * output  = r.data();
        #pragma omp simd aligned(output)
        for(size_t i = 0; i < m; i++) output[i] = start + i*delta;
        return r;
    }

    template<typename T>
    sample_vector<T> regspace(int start, int end, T delta=(T)0.0)
    {
        if(delta==0.0) delta = start <= end? 1.0 : -1.0;                
        int m = std::floor((end-start)/delta);
        if(start > end) m++;
        sample_vector<T> r(m);
        T * output  = r.data();
        #pragma omp simd aligned(output)
        for(size_t i = 0; i < m; i++) output[i] = start + i*delta;
        return r;
    }
    

    // linear or acyclic
    template<typename T>
    sample_vector<T> linear_convolution(sample_vector<T> & h, sample_vector<T> & x)
    {
        int M=h.size();
        int N=x.size();
        sample_vector<T> y(M+N-1);        
        T * output = y.data();		
        #pragma omp parallel for
        for(int i=0;i<M+N-1;i++)
        {  
            output[i]=0;                    
            #pragma omp simd aligned(h,x,output)
            for(int j=0;j<=i;j++)
                output[i]+=x[j]*h[i-j];
        }
        return y;
    }
    
    template<typename T>
    void cshiftr(sample_vector<T> & v)
    {
        T t = v[0];
        T * data = v.data();
        #pragma omp simd aligned(data)
        for(size_t i = 1; i < v.size(); i++)
            v[i-1] = v[i];
        v[v.size()-1] = t;
    }
    template<typename T>
    void cshiftl(sample_vector<T> & v)
    {
        T t = v[v.size()-1];
        T * data = v.data();
        #pragma omp simd aligned(data)
        for(size_t i = v.size()-1; i > 0; i--)
            v[i] = v[i-1];
        v[0] = t;
    }
    template<typename T>
    void cshiftleft(sample_vector<T> & v, size_t steps)
    {        
        for(size_t i = 0; i < steps; i++) cshiftl(v);
    }
    template<typename T>
    void cshiftright(sample_vector<T> & v, size_t steps)
    {
        for(size_t i = 0; i < steps; i++) cshiftr(v);
    }    
    template<typename T>
    void cshift(sample_vector<T> & v, int steps)
    {
        if(steps > 0) cshiftleft(v,steps);
        else if(steps < 0) cshiftright(v,std::abs(steps));
    }
    template<typename T>
    sample_vector<T> lag(sample_vector<T> & v, int lag)
    {
        sample_vector<T> r(v.size());
        if(lag < 0)
        {
            int n = (v.size() + lag) % v.size();
            int s = std::abs(lag) % v.size();
            memcpy(r.data(), v.data() + n,s);
            memcpy(r.data()+ n, v.data(), v.size() - s);
        }
        else
        {
            int n = (v.size() - lag) % v.size();
            int s = std::abs(lag) % v.size();
            memcpy(r.data(), v.data() + n, s);
            memcpy(r.data() + n, v.data(), v.size() - s);
        }
    }
    
    
    // convolution
    // will be enormously slow         
    
    template<typename T>
    sample_vector<T> linear_convolution(T * h, T * x, int M, int N)
    {
        sample_vector<T> y(M+N-1);        
        T * output = y.data();
        // nothing can really improve this but why not                
        #pragma omp parallel for        
        for(int i=0;i<M+N-1;i++)
        {  
            output[i]=0;        
            #pragma omp simd aligned(h,x,output)         
            for(int j=0;j<=i;j++)
                output[i]+=x[j]*h[i-j];
        }
        return y;
    }
            
    template<typename T>
    sample_vector<T> circular_convolution(T * h, T * x, int M, int N)
    {        
        sample_vector<T> y(M);  
        T * output = y.data();		
        #pragma omp parallel for              
        for(int n=0;n < M;n++)
        { 
            output[n]=0;                             
            #pragma omp simd aligned(h,x,output)     
            for(int k=0;k < N;k++)                
            {
                int j = (k>n)? n-k+N : n-k;
                output[n]=output[n]+x[k]*h[j];
            }
        }    
        return y;    
    }   
    template<typename T>
    sample_vector<T> circular_convolution(sample_vector<T>& h, sample_vector<T>& x)
    {
        int M = h.size();
        int N = x.size();
        sample_vector<T> y(M);  
        T * output = y.data();		
        #pragma omp parallel for              
        for(int n=0;n < M;n++)
        { 
            output[n]=0;                                  
            #pragma omp simd aligned(h,x,output)     
            for(int k=0;k < N;k++)                
            {
                int j = (k>n)? n-k+N : n-k;
                output[n]=output[n]+x[k]*h[j];
            }
        }    
        return y;    
    }       
    
       template<typename T>
    void cshiftr(size_t n, T * v)
    {
        T t = v[0];
        T * data = v.data();
        #pragma omp simd aligned(data)
        for(size_t i = 1; i < n; i++)
            v[i-1] = v[i];
        v[n-1] = t;
    }
    template<typename T>
    void cshiftl(size_t n, T* v)
    {
        T t = v[n-1];
        #pragma omp simd
        for(size_t i = n-1; i > 0; i--)
            v[i] = v[i-1];
        v[0] = t;
    }
    template<typename T>
    void cshiftleft(size_t n, T * v, size_t steps)
    {        
        for(size_t i = 0; i < steps; i++) cshiftl(n,v);
    }
    template<typename T>
    void cshiftright(size_t n, T * v, size_t steps)
    {
        for(size_t i = 0; i < steps; i++) cshiftr(n,v);
    }    
    template<typename T>
    void cshift(size_t n, T * v, int steps)
    {
        if(steps > 0) cshiftleft(n,v,steps);
        else if(steps < 0) cshiftright(n,v,std::abs(steps));
    }

    template<typename T>
    void lag(size_t n, T * v, int lag)
    {
        sample_vector<T> r(n);
        if(lag < 0)
        {
            int n = (n + lag) % n;
            int s = std::abs(lag) % n;
            memcpy(r.data(), v + n,s);
            memcpy(r.data()+ n, v, n - s);
        }
        else
        {
            int n = (n - lag) % n;
            int s = std::abs(lag) % n;
            memcpy(r.data(), v + n, s);
            memcpy(r.data() + n, v, n - s);
        }
        memcpy(v,r.data(),n*sizeof(T));
    }
}
