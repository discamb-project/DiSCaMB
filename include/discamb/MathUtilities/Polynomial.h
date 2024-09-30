#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include <vector>
#include <iostream>
#include <cmath>


namespace MChLib{

template<typename T1,typename T2>
class Polynomial;


template<typename U1,typename U2>
const Polynomial<U1,U2> operator*(const U1 &,const Polynomial<U1,U2> &);		

template<typename T1,typename T2>
const Polynomial<T1,T2> operator*(const Polynomial<T1,T2> &,const T1 &);	

template<typename T1,typename T2>
const Polynomial<T1,T2> operator*(const Polynomial<T1,T2> &,const Polynomial<T1,T2> &);

template<typename U1,typename U2>
const Polynomial<U1,U2> operator+(const Polynomial<U1,U2> &,const Polynomial<U1,U2> &);	

template<typename U1,typename U2>
const Polynomial<U1,U2> operator-(const Polynomial<U1,U2> &,const Polynomial<U1,U2> &);	

template<typename U1,typename U2>
const Polynomial<U1,U2> operator/(const Polynomial<U1,U2> &,const U1 &);	

template<typename U1,typename U2>
std::ostream &operator<<(std::ostream &,const Polynomial<U1,U2> &);

template<typename T1,typename T2>
class Polynomial
{
public:
	Polynomial();
	virtual ~Polynomial();
	Polynomial(const Polynomial<T1,T2> &);	
	Polynomial(const T1 &,unsigned int dimension);
	Polynomial<T1,T2> &operator=(const Polynomial<T1,T2> &);
	Polynomial<T1,T2> &operator=(const T1 &);
	void addTerm(const T1 &coefficient,const std::vector<T2> &powers);
	void addTerms(const std::vector<T1> &coefficients,
	              const std::vector<std::vector<T2> > &powers); 
	void get(std::vector<T1> &coefficients,
	         std::vector<std::vector<T2> > &powers) const;
	unsigned int getDimension() const;    
	T1 calculate(const std::vector<T1> &) const;  
	void clear();   	
	
	
	Polynomial &operator+=(const Polynomial &);
	Polynomial &operator-=(const Polynomial &);	
	Polynomial &operator*=(const Polynomial &);	
	Polynomial &operator*=(const T1 &);		
	Polynomial &operator/=(const T1 &);		
	
	void simplify(const T1 &treshold=0);


	friend const Polynomial<T1,T2> operator*<T1,T2>(const Polynomial<T1,T2> &,const Polynomial<T1,T2> &);

	friend const Polynomial<T1,T2> operator+<T1,T2>(const Polynomial<T1,T2> &,const Polynomial<T1,T2> &);	

	friend const Polynomial<T1,T2> operator-<T1,T2>(const Polynomial<T1,T2> &,const Polynomial<T1,T2> &);	

	friend const Polynomial<T1,T2> operator*<T1,T2>(const Polynomial<T1,T2> &,const T1 &);	

	friend const Polynomial<T1,T2> operator*<T1,T2>(const T1 &,const Polynomial<T1,T2> &);		


	friend const Polynomial<T1,T2> operator/<T1,T2>(const Polynomial<T1,T2> &,const T1 &);	


	friend std::ostream &operator<<<T1,T2>(std::ostream &,const Polynomial<T1,T2> &);
protected:
	void checkAndAdjustDimension(const Polynomial<T1,T2> &q,Polynomial<T1,T2> &p);
	void setDimension(unsigned int);
private:
	std::vector<T1> mCoefficients;
	std::vector<std::vector<T2> > mPowers;	
	unsigned int mDimension;
	//unsigned int mNTerms;
	bool reduce();
};

#include "Polynomial.hpp"

}

#endif /*POLYNOMIAL_H_*/
