//#include "Polynomial.h"
//#include <iostream>
//#include <cmath>
//#include <MChMathMyUtilities/MyUtilities.h>

using namespace std;

//-----------------------------------------------------------------------------

template <typename T1,typename T2>
Polynomial<T1,T2>::Polynomial()
{
	mDimension=0;
}

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
Polynomial<T1,T2>::~Polynomial()
{
}

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
Polynomial<T1,T2>::Polynomial(
const Polynomial<T1,T2> &p)
{
	mCoefficients=p.mCoefficients;
	mPowers=p.mPowers;
	mDimension=p.mDimension;
}

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
Polynomial<T1,T2>::Polynomial(
const T1 &d,
unsigned int dimension)
{
	mDimension=dimension;
	mCoefficients.assign(1,d);
	mPowers.assign(1,vector<T2>(dimension,0));
}
	
//-----------------------------------------------------------------------------

template<typename T1,typename T2>
Polynomial<T1,T2> &Polynomial<T1,T2>::operator=(
const Polynomial<T1,T2> &p)
{
	mCoefficients=p.mCoefficients;
	mPowers=p.mPowers;	
	mDimension=p.mDimension;
	return *this;
}

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
void Polynomial<T1,T2>::addTerm(
const T1 &coefficient,
const std::vector<T2> &powers)
{
	std::vector<T2> variable_powers=powers;
	unsigned int i;

	if(mPowers.size()>0)
	   {
	   if(variable_powers.size()>mDimension)
	      setDimension(variable_powers.size());
	   if(variable_powers.size()<mDimension)
	      variable_powers.resize(mDimension,0);
	   }
	else
	   {
	   mDimension=powers.size();	
	   }      
	mCoefficients.push_back(coefficient);
	mPowers.push_back(variable_powers);  
}

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
void Polynomial<T1,T2>::addTerms(
const std::vector<T1> &coefficients,
const std::vector<std::vector<T2> > &powers)
{
	unsigned int i;
	
	if(powers.size()!=coefficients.size())
	   {
	   	cerr<< "ERROR: coefficients and powers have different size"
	   	    << " in Polynomial::addTerms" << endl;
	   exit(1);	    
	   }
	
	for(i=0;i<powers.size();i++)
	   addTerm(coefficients[i],powers[i]);
}


//-----------------------------------------------------------------------------

template<typename T1,typename T2>
void Polynomial<T1,T2>::clear()
{
	mCoefficients.clear();
	mPowers.clear();
	mDimension=0;
}

//-----------------------------------------------------------------------------
	
template<typename T1,typename T2>
Polynomial<T1,T2> &Polynomial<T1,T2>::operator+=(
const Polynomial<T1,T2> &q)
{
	unsigned int i;
	Polynomial<T1,T2> p;
	
	checkAndAdjustDimension(q,p);
	   
	for(i=0;i<p.mPowers.size();i++)
	   {
	   mCoefficients.push_back(p.mCoefficients[i]);
	   mPowers.push_back(p.mPowers[i]);   	
	   }
	return *this;
}

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
Polynomial<T1,T2> &Polynomial<T1,T2>::operator-=(
const Polynomial<T1,T2> &q)
{
	unsigned int i;
	
	Polynomial<T1,T2> p;
	
	checkAndAdjustDimension(q,p);

		   
	for(i=0;i<p.mPowers.size();i++)
	   {
	   mCoefficients.push_back(-p.mCoefficients[i]);
	   mPowers.push_back(p.mPowers[i]);  	
	   }
	return *this;
}

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
Polynomial<T1,T2> &Polynomial<T1,T2>::operator*=(
const Polynomial<T1,T2> &q)
{
	unsigned int i,j,k,m;
	vector<vector<T2> > newPowers;
	vector<T1> newCoefficients;

	Polynomial<T1,T2> p;
	
	checkAndAdjustDimension(q,p);

	newPowers.resize(p.mPowers.size()*mPowers.size(),vector<T2>(mDimension));
	newCoefficients.resize(p.mPowers.size()*mPowers.size());
	      		
	m=0;
	for(i=0;i<mPowers.size();i++)
	   for(j=0;j<p.mPowers.size();j++)
	      {
	      newCoefficients[m]=mCoefficients[i]*p.mCoefficients[j];
	      for(k=0;k<mDimension;k++)
	         newPowers[m][k]=mPowers[i][k]+p.mPowers[j][k];
	      m++;   
	      }

	mPowers=newPowers;
	mCoefficients=newCoefficients;
	return *this;
}

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
Polynomial<T1,T2> &Polynomial<T1,T2>::operator*=(
const T1 &d)
{
	unsigned int i;
	for(i=0;i<mPowers.size();i++)
	   mCoefficients[i]*=d;
	return *this;   
}		

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
Polynomial<T1,T2> &Polynomial<T1,T2>::operator/=(
const T1 &d)
{
	unsigned int i;
	for(i=0;i<mPowers.size();i++)
	   mCoefficients[i]/=d;
	return *this;
}
	

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
void Polynomial<T1,T2>::simplify(
const T1 &treshold)
{
	unsigned int i;
	while(reduce());
	for(i=0;i<mPowers.size();i++)
	  if(fabs(mCoefficients[i])<=treshold)
	     {
	     mCoefficients.erase(mCoefficients.begin()+i);
	     mPowers.erase(mPowers.begin()+i);
	     i--;
	     }
	  
}

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
bool Polynomial<T1,T2>::reduce()
{
	unsigned int i,j;
	bool found=false;
	
	for(i=0;i<mPowers.size();i++)
	   for(j=i+1;j<mPowers.size();j++)
	      {
	      if(mPowers[i]==mPowers[j])
	         {
	         mCoefficients[i]+=mCoefficients[j];
	         found=true;
	         mCoefficients.erase(mCoefficients.begin()+j);
	         mPowers.erase(mPowers.begin()+j);
	         i=mPowers.size();
	         j=mPowers.size();	         
	         }	
	      }
	return found;
}
	
//-----------------------------------------------------------------------------

template<typename T1,typename T2>
void Polynomial<T1,T2>::get(
vector<T1> &coefficients,
vector<vector<T2> > &powers)
 const
{
 	coefficients=mCoefficients;
 	powers=mPowers;
}

//-----------------------------------------------------------------------------	

template<typename T1,typename T2>
const Polynomial<T1,T2> operator*(
const Polynomial<T1,T2> &p1,
const Polynomial<T1,T2> &p2)
{
	Polynomial<T1,T2> p(p1);
	p*=p2;
	return p;
}

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
const Polynomial<T1,T2> operator+(
const Polynomial<T1,T2> &p1,
const Polynomial<T1,T2> &p2)
{
	Polynomial<T1,T2> p(p1);
	p+=p2;
	return p;
}

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
const Polynomial<T1,T2> operator-(
const Polynomial<T1,T2> &p1,
const Polynomial<T1,T2> &p2)
{
	Polynomial<T1,T2> p(p1);
	p-=p2;
	return p;
}	

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
const Polynomial<T1,T2> operator*(
const Polynomial<T1,T2> &p,
const T1 &d)
{
	Polynomial<T1,T2> result(p);
	result*=d;
	return result;
}

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
const Polynomial<T1,T2> operator*(
const T1 &d,
const Polynomial<T1,T2> &p)
{
	Polynomial<T1,T2> result(p);
	result*=d;
	return result;	
}

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
const Polynomial<T1,T2> operator/(
const Polynomial<T1,T2> &p,
const T1 & d)
{
	Polynomial<T1,T2> result(p);
	result/=d;
	return result;		
}

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
ostream  &operator<<(
ostream &out,
const Polynomial<T1,T2> &p)
{
	unsigned int i,k,n;
	bool madeNewLine=false;
	
	
	for(i=0;i<p.mPowers.size();i++)
	   {
	   out<< p.mCoefficients[i] << "    ";
	   n=p.mPowers[i].size();
	   for(k=0;k<n;k++)
	      {
	      out<< p.mPowers[i][k] << " ";
	      if((k+1)%12==0)
	         {
	         madeNewLine=true;
	         out<< endl;	
	         }
	      else
	         madeNewLine=false;   
	      }      
	   if(!madeNewLine)
	      out<< endl;   
	   }  
	return out;     
}

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
void Polynomial<T1,T2>::checkAndAdjustDimension(
const Polynomial<T1,T2> &q,
Polynomial<T1,T2> &p)
{
	p=q;

	if(mPowers.size()>0&&p.mPowers.size()>0)
	   {
	   if(p.mDimension!=mDimension)
	      {
	      if(p.mDimension>mDimension)
			  setDimension(p.mDimension);
		  if(p.mDimension<mDimension)
			  p.setDimension(mDimension);
	      }
	   }
	else
	   {
	   if(mDimension==0)
	      mDimension=p.mDimension;
	   }      
	
}

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
unsigned int Polynomial<T1,T2>::getDimension()
 const
{
	return mDimension;
}
	
//-----------------------------------------------------------------------------

template<typename T1,typename T2>
Polynomial<T1,T2> &Polynomial<T1,T2>::operator=(
const T1 &d)
{
	if(d!=0)
	   {
	   mPowers.assign(1,vector<T2>(mDimension,0));
	   mCoefficients.assign(1,d);
	   }
	else
	   {
	   mPowers.clear();	
	   mCoefficients.clear();
	   }   
	return *this;   
}

//-----------------------------------------------------------------------------

template<typename T1,typename T2>
T1 Polynomial<T1,T2>::calculate(
const std::vector<T1> &variables)
 const
{
	T1 result=0,d;
	unsigned int i,j;
	
	for(i=0;i< mCoefficients.size();i++)
	   {
	   d=mCoefficients[i];
	   for(j=0;j<mDimension;j++)
	      d*=power(variables[j],mPowers[i][j]);
	   result+=d;   
	   }
	return result;
}


//-----------------------------------------------------------------------------

template<typename T1,typename T2>
void Polynomial<T1,T2>::setDimension(
unsigned int dimension)
{
	unsigned int i;
	mDimension=dimension;
	for(i=0;i<mPowers.size();i++)
		mPowers[i].resize(mDimension,0);
}


