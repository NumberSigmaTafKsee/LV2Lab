%module floatbuffer
%{
#include <vector>
%}

%inline %{
	class FloatBuffer
	{
	private:
		float * buffer;
		size_t  size;
	public:
	
		FloatBuffer(size_t n, void * ptr) {
			buffer = (float*)ptr;
			size   = n;
		}
		float  __getitem__(size_t i) { return buffer[i]; }
		void   __setitem__(size_t i, float v) { buffer[i] = v; }
	};
%}
