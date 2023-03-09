#include "stdsamples_casino.hpp"

using namespace std;
using namespace Casino;
using namespace Casino::IPP;

using ComplexFloat = std::complex<float>;
using ComplexDouble= std::complex<double>;

template<typename T>
struct generator {
	int start = 0;
	
	T operator()() {
		T r = start*start;
		start++;
		return r;
	}
};

int main()
{
	try {
		IPPArray<Ipp32f> af(100);
		IPPArray<Ipp64f> ad(100);
		IPPArray<Ipp32fc> acf(100);
		IPPArray<Ipp64fc> acd(100);
		
		std::generate(af.begin(),af.end(),generator<float>());
		std::cout << af << std::endl;
	}
	catch(std::runtime_error & what) {
		std::cout << what.what() << std::endl;
	}
}

