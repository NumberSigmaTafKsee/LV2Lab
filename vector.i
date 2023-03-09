%module vector
%{
#include <vector>
%}

%include "std_vector.i"
 
%template(VectorXf) std::vector<float>;
%template(VectorXd) std::vector<double>;
%template(VectorXi8) std::vector<char>;
%template(VectorXu8) std::vector<unsigned char>;
%template(VectorXi16) std::vector<short>;
%template(VectorXu16) std::vector<unsigned short>;
%template(VectorXi32) std::vector<int>;
%template(VectorXu32) std::vector<unsigned int>;
%template(VectorXi64) std::vector<long long int>;
%template(VectorXu64) std::vector<unsigned long long int>;

%inline %{

	template<typename T>
	struct List : public std::vector<T>
	{
		List() = default;
		List(size_t n) : std::vector<T>(n) {
		
		}
		using std::vector<T>::resize;
		using std::vector<T>::size;
		using std::vector<T>::front;
		using std::vector<T>::back;
		using std::vector<T>::pop_back;
		using std::vector<T>::push_back;
		
		void pushBack(const T & value) {
			push_back(value); 
		}
		T popBack() {
			T r = back();
			pop_back();
			return r;
		}
		void pushFront(const T& value) {
			if(size() == 0) {
				push_back(value);
			}
			else {
				this->push_back(0);
				memcpy(this->data()+1,this->data(),(this->size()-2)*sizeof(T));
				(*this)[0] = value;
			}
		}
		
		T popFront() {
			if(size() == 0) throw std::runtime_error("pop empty list");
			T r = (*this)[0];
			memcpy(this->data(),this->data()+1,(this->size()-2)*sizeof(T));
			this->pop_back();
			return r;
		}
		
		T __getitem__(size_t i) { return (*this)[i]; }
		void __setitem__(size_t i, const T & value) { (*this)[i] = value; }
		
		void insertItem(size_t i, const T & value)
		{
			if(i >= size()) return;
			this->insert(this->begin()+i,value);
		}
		T removeItem(size_t i) {
			T r = (*this)[i];
			this->erase(this->begin()+i);
			return r;
		}
		void removeRange(size_t start, size_t end) {
			this->erase(this->begin()+start,this->begin()+end);
		}
			
	};
%}

%template (ListXf) List<float>;
