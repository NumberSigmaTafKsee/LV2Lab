#pragma once

template<typename T>
struct RingBuffer : public sample_vector<T>
{
    size_t r=0;
    size_t w=0;

    RingBuffer(size_t n) {
        sample_vector<T>::resize(n);
    }

    void set_write_position(size_t n) {
        w = (n % sample_vector<T>::size());
    }  
    T    get() {
        return (*this)[r++];
    }
    void push(T x) {
        (*this)[w++] = x;
        w = (w % sample_vector<T>::size());
    }
    T linear() {
        T x = (*this)[r];
        T x1= (*this)[r++];
        T f = x - floor(x);
        r = r % sample_vector<T>::size();
        return linear_interpolate(x,x1,f);        
    }
    T cubic() {
        T xm1= (*this)[(r-1) % sample_vector<T>::size()];
        T x = (*this)[r];
        T x1= (*this)[(r+1) % sample_vector<T>::size()];
        T x2= (*this)[(r+2) % sample_vector<T>::size()];
        T f = x - floor(x);
        r++;
        r = r % sample_vector<T>::size();
        return cubic_interpolate(f,xm1,x,x1,x2);        
    }
};
