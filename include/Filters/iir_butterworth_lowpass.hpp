#pragma once
#include "iir_filters.hpp"
#include "iir_butterworth_filters.hpp"

namespace iir_filters
{
	struct ButterworthLowpassFilter
	{
		int order;
		DspFloatType sampleRate,frequency,q;    
		Filters::BiquadFilterCascade filter;

		ButterworthLowpassFilter(int order, DspFloatType sr)
		{
			this->order = order;
			this->sampleRate = sr;        
			setFilter(1000.0,0.5);
		}
		void setFilter(DspFloatType f, DspFloatType Q) {
			auto x = Filters::butterlp(order,Q);
			frequency = f;
			if(f <= 30.0) f = 30.0f;
			if(f >= sampleRate/2.0) f = sampleRate/2.0;
			q = Q;
			auto c = AnalogBiquadCascade(x,f,sampleRate);
			filter.setCoefficients(c);
		}
		void setCutoff(DspFloatType f) {
			if(f <= 30) return;
			if(f >= sampleRate/2.0) return;
			setFilter(f,q);
		}
		void setQ(DspFloatType Q) {
			if(Q < 0.5) Q = 0.5;
			if(Q > 999) Q = 999;
			setFilter(frequency,Q);
		}
		DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
			return filter.Tick(I,A,X,Y);
		}
		void ProcessBlock(size_t n, float * in, float * out) {
			for(size_t i = 0; i < n; i++) in[i] = Tick(out[i]);
		}
		std::vector<DspFloatType> impulse_response(size_t n)
		{
			std::vector<DspFloatType> r(n);
			r[0] = Tick(1.0);
			for(size_t i = 1; i < n; i++) r[i] = Tick(0);
			return r;
		}
	};
}
