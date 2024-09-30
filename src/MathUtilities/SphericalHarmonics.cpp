#include "discamb/MathUtilities/SphericalHarmonics.h"
#include "discamb/MathUtilities/MathUtilities.h"
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <limits>


#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif



using namespace std;

namespace MChLib{

namespace{
void makePolynomial(const vector<double> &coefficients1,
					const vector<vector<unsigned int> >&powers1,
					Polynomial<double,unsigned int> &polynomial,
					const vector<vector<double> > &coordinates);

double integrate(const Polynomial<double,unsigned int> &);

const double extremumAbsoluteValue[]=
{0.0795775,0.31831,0.31831,0.31831,0.375,0.375,0.413497,0.375,
0.375,0.424413,0.3849,0.441118,0.489707,0.441118,0.3849,0.424413,
0.46875,0.405949,0.425046,0.500559,0.55534,0.500559,0.425046,
0.405949,0.46875,0.509296,0.429325,0.43018,0.466054,0.554257,
0.613916,0.554257,0.466054,0.43018,0.429325,0.509296,0.546875,
0.452907,0.442045,0.459698,0.505212,0.603432,0.667334,0.603432,
0.505212,0.459698,0.442045,0.452907,0.546875,0.582052,0.476037,
0.456558,0.463212,0.489751,0.542215,0.648999,0.716766,0.648999,
0.542215,0.489751,0.463212,0.456558,0.476037,0.582052};


}

//-----------------------------------------------------------------------------
//   STATIC COMPONENTS INITIALIZATION
//-----------------------------------------------------------------------------

vector<vector<double> > SphericalHarmonics::mFlm;
unsigned int SphericalHarmonics::mMaxL=0;

vector<double> SphericalHarmonics::mPowerXValue(1);
vector<double> SphericalHarmonics::mPowerYValue(1);
vector<double> SphericalHarmonics::mPowerZValue(1);		
vector<double> SphericalHarmonics::mPowerRValue(1);

// sph

vector<vector<vector<vector<unsigned int> > > > SphericalHarmonics::mXyzPowers(1,
             vector<vector<vector<unsigned int> > >(1,
                vector<vector<unsigned int> >(1,
                   vector<unsigned int>(3,0)))); 

vector<vector<vector<double> > > SphericalHarmonics::mCoefficients
                                     (1,vector<vector<double> >
                                     (1,vector<double>
                                     (1,sqrt(0.25/M_PI) ))); 

// sph derivatives

vector<vector<vector<vector<double> > > > SphericalHarmonics::
                   mDerivativeCoefficients(1,vector<vector<vector<double> > >
                                          (1,vector<vector<double> >
                                          (3)));

vector<vector<vector<vector<vector<unsigned int> > > > > SphericalHarmonics::
                   mDerivativeXyzPowers(1,vector<vector<vector<vector<unsigned int> > > >
                                       (1,vector<vector<vector<unsigned int> > >
                                       (3)));//,vector<vector<unsigned int> >)));

vector<vector<vector<vector<double> > > > SphericalHarmonics::
                   mSolidDerivativeCoefficients(1,vector<vector<vector<double> > >
                                          (1,vector<vector<double> >
                                          (3)));

vector<vector<vector<vector<vector<unsigned int> > > > > SphericalHarmonics::
                   mSolidDerivativeXyzPowers(1,vector<vector<vector<vector<unsigned int> > > >
                                       (1,vector<vector<vector<unsigned int> > >
                                       (3)));//,vector<vector<unsigned int> >)));
                                       
// sph second derivatives

vector<vector<vector<vector<vector<double> > > > > SphericalHarmonics::
                   mSecondDerivativeCoefficients(1,vector<vector<vector<vector<double> > > >
                                                (1,vector<vector<vector<double> > >
                                                (3,vector<vector<double> >
                                                (3))));
                   	
vector<vector<vector<vector<vector<vector<unsigned int> > > > > > 
              SphericalHarmonics::mSecondDerivativeXyzPowers
                   (1,vector<vector<vector<vector<vector<unsigned int> > > > >
                   (1,vector<vector<vector<vector<unsigned int> > > >
                   (3,vector<vector<vector<unsigned int> > >
                   (3))));

vector<vector<vector<vector<vector<double> > > > > SphericalHarmonics::
                   mSolidSecondDerivativeCoefficients(1,vector<vector<vector<vector<double> > > >
                                                (1,vector<vector<vector<double> > >
                                                (3,vector<vector<double> >
                                                (3))));
                   	
vector<vector<vector<vector<vector<vector<unsigned int> > > > > > 
              SphericalHarmonics::mSolidSecondDerivativeXyzPowers
                   (1,vector<vector<vector<vector<vector<unsigned int> > > > >
                   (1,vector<vector<vector<vector<unsigned int> > > >
                   (3,vector<vector<vector<unsigned int> > >
                   (3))));

bool SphericalHarmonics::initialized=false;


//-----------------------------------------------------------------------------
//                        METHOD DEFINITIONS
//-----------------------------------------------------------------------------


SphericalHarmonics::SphericalHarmonics()
{
}

//-----------------------------------------------------------------------------

SphericalHarmonics::~SphericalHarmonics()
{
}

//-----------------------------------------------------------------------------

double SphericalHarmonics::calculate(
unsigned int l,
int m,
double x,
double y,
double z,
SphericalHarmonicForm form)
{
	double result=0;
	double r;
	int absM;
	if(!initialized)
	   initialize();

	r=sqrt(x*x+y*y+z*z);
	absM=abs(m);
	if(form==DensityNormalized||form==SolidDensityNormalized)
	   {
	   result=mFlm[l][absM]*numerator(l,m,x,y,z);
	   if(form==DensityNormalized)
	      result*=rPowerMinusL(l,r);
	   }
	else
	   {   
	   result=numerator(l,m,x,y,z);
	   if(form==WaveFunctionNormalized)
	      result*=rPowerMinusL(l,r);
	   }   	
                     
	return result;
}


//-----------------------------------------------------------------------------

double SphericalHarmonics::calculate(
unsigned int l,
int m,
double x,
double y,
double z,
double r,
SphericalHarmonicForm form)
{
	double result=0;
	int absM;
	
	if(!initialized)
	   initialize();
	absM=abs(m);
	if(form==DensityNormalized||form==SolidDensityNormalized)
	   {
	   result=mFlm[l][absM]*numerator(l,m,x,y,z);
	   if(form==DensityNormalized)
	      result*=rPowerMinusL(l,r);
	   }
	else
	   {   
	   result=numerator(l,m,x,y,z);
	   if(form==WaveFunctionNormalized)
	      result*=rPowerMinusL(l,r);
	   }   	
	                       
	return result;
}

//-----------------------------------------------------------------------------

void SphericalHarmonics::calculateGradient(
unsigned int l,
int m,
double x,
double y,
double z,
double r,
double (&gradient)[3],
SphericalHarmonicForm form)
{
	unsigned int i;
	double a,da,b,db;
	double coefficient;
	int absM;
	
	if(!initialized)
	   initialize();
	absM=abs(m);
	
	if(form==DensityNormalized||form==SolidDensityNormalized)
	   coefficient=mFlm[l][absM];
	else
	   coefficient=1.0;   
	
	for(i=0;i<3;i++)
	   {
	   a=numerator(l,m,x,y,z);
	   b=rPowerMinusL(l,r);
	   da=numeratorDerivative(i,l,m,x,y,z,r);
	   db=rPowerMinusLDerivative(i,l,x,y,z,r);
	   if(form==DensityNormalized||form==WaveFunctionNormalized)
	      gradient[i]=coefficient*(da*b+db*a);
	   else
	      gradient[i]=coefficient*da;   
	   }           
}

//-----------------------------------------------------------------------------

void SphericalHarmonics::calculateHessian(
unsigned int l,
int m,
double x,
double y,
double z,
double r,
double (&hessian)[3][3],
SphericalHarmonicForm form)
{
	unsigned int i,j;
	double num,d_num_i,d_num_j,dd_num,rpml,d_rpml_i,d_rpml_j,dd_rpml;
	double coefficient;
	int absM;
	if(!initialized)
	   initialize();

	absM=abs(m);
	if(form==DensityNormalized||form==SolidDensityNormalized)
	   coefficient=mFlm[l][absM];
	else
	   coefficient=1.0;   
	
	for(i=0;i<3;i++)
	   for(j=0;j<3;j++)
	      {
	      num=numerator(l,m,x,y,z);
	      d_num_i=numeratorDerivative(i,l,m,x,y,z,r);
	      d_num_j=numeratorDerivative(j,l,m,x,y,z,r);
	      dd_num=numeratorSecondDerivative(i,j,l,m,x,y,z,r);
	      rpml=rPowerMinusL(l,r);
	      d_rpml_i=rPowerMinusLDerivative(i,l,x,y,z,r);
	      d_rpml_j=rPowerMinusLDerivative(j,l,x,y,z,r);
	      dd_rpml=rPowerMinusLSecondDerivative(i,j,l,x,y,z,r);
	      if(form==DensityNormalized||form==WaveFunctionNormalized)	      
	         hessian[i][j]=coefficient*(dd_num*rpml+d_num_i*d_rpml_j+
	                                 d_num_j*d_rpml_i+num*dd_rpml);
	      else
	         hessian[i][j]=coefficient*dd_num;
	      }
}

//-----------------------------------------------------------------------------

void SphericalHarmonics::initialize()
{
	unsigned int l;
	mFlm.resize(8);
	for(l=0;l<=7;l++)
	   mFlm[l].resize(2*l+1);
	   
	mFlm[0][0]=0.28209479177388;
	mFlm[1][0]=0.65147001587056;
	mFlm[1][1]=0.65147001587056;
	mFlm[2][0]=0.65552905835525;
	mFlm[2][1]=0.68646842464783;
	mFlm[2][2]=0.68646842464783;
	mFlm[3][0]=0.65613421114746;
	mFlm[3][1]=0.70087743932709;
	mFlm[3][2]=0.69189513695864;
	mFlm[3][3]=0.71929123343438;
	mFlm[4][0]=0.65620991154888;
	mFlm[4][1]=0.70847465461627;
	mFlm[4][2]=0.6987955686655;
	mFlm[4][3]=0.70616251710906;
	mFlm[4][4]=0.74899845670137;
	mFlm[5][0]=0.65617176926179;
	mFlm[5][1]=0.71306675774595;
	mFlm[5][2]=0.70407303654785;
	mFlm[5][3]=0.7054752262822;
	mFlm[5][4]=0.72266090165344;
	mFlm[5][5]=0.7759136537823;
	mFlm[6][0]=0.65611016832343;
	mFlm[6][1]=0.71609595329005;
	mFlm[6][2]=0.70800876706179;
	mFlm[6][3]=0.70703232886769;
	mFlm[6][4]=0.71554818711837;
	mFlm[6][5]=0.73945146220726;
	mFlm[6][6]=0.8004796889213;
	mFlm[7][0]=0.65605000193133;
	mFlm[7][1]=0.71822036312964;
	mFlm[7][2]=0.71099470247414;
	mFlm[7][3]=0.70896558772175;
	mFlm[7][4]=0.71344649338342;
	mFlm[7][5]=0.72696198663051;
	mFlm[7][6]=0.75586913625466;
	mFlm[7][7]=0.82308112861385;

	initialized=true;
}

//-----------------------------------------------------------------------------

double SphericalHarmonics::rPowerMinusL(
unsigned int l,
double r)
{
	//return 1.0/power(r,l);
    return 1.0 / pow(r, l);
}

//-----------------------------------------------------------------------------  

double SphericalHarmonics::rPowerMinusLDerivative(
unsigned int var,
unsigned int l,
double x,
double y,
double z,
double r)
{
	double v=0,result=0,a,b,c;
	
	if(var==0) v=x;
	if(var==1) v=y;
	if(var==2) v=z;
	
	a=-static_cast<double>(l);
	b=v;
	c=pow(r,l+2);
	result=a*b/c;
	
	return result;
}

//-----------------------------------------------------------------------------	                          

double SphericalHarmonics::rPowerMinusLSecondDerivative(
unsigned int var1,
unsigned int var2,
unsigned int l,
double x,
double y,
double z,
double r)
{
	double v1=0,v2=0,result;
	
	if(var1==0) v1=x;
	if(var1==1) v1=y;
	if(var1==2) v1=z;

	if(var2==0) v2=x;
	if(var2==1) v2=y;
	if(var2==2) v2=z;
	
	
	result=static_cast<double>((2+l)*l)*v1*v2/pow(r,l+4);
	
	if(var1==var2)
	   result-=static_cast<double>(l)/pow(r,l+2);

	return result;   
}


//-----------------------------------------------------------------------------

void SphericalHarmonics::updateSph(
unsigned int l)
{
	unsigned int i;
	int intI;

	mXyzPowers.resize(l+1);
	mCoefficients.resize(l+1);
	for(i=mMaxL+1;i<=l;i++)
	   {
	   mXyzPowers[i].resize(2*l+1);
	   mCoefficients[i].resize(2*l+1);	
	   intI=i;
	   getRealSolidHarmonics(i,mXyzPowers[i],mCoefficients[i]);  
	   }

	mPowerXValue.resize(l+1);
	mPowerYValue.resize(l+1);
	mPowerZValue.resize(l+1);
	mPowerRValue.resize(l+1);

	updateSphDerivatives(l);   
}

//-----------------------------------------------------------------------------

void SphericalHarmonics::updateSphDerivatives(
unsigned int l)
{
	int i,m,c,intL;
	unsigned int j;
	intL=l;
	vector<vector<unsigned int> > powers;
	// without r^-l factor
		
	mSolidDerivativeXyzPowers.resize(l+1);
	mSolidDerivativeCoefficients.resize(l+1);
	for(i=mMaxL+1;i<=intL;i++)
	   {
	   mSolidDerivativeXyzPowers[i].resize(2*l+1);
	   mSolidDerivativeCoefficients[i].resize(2*l+1);	
	   for(m=-i;m<=i;m++)
	      {
	      mSolidDerivativeXyzPowers[i][i+m].resize(3);
	      mSolidDerivativeCoefficients[i][i+m].resize(3);		
	      for(c=0;c<3;c++)
	         differentiate(c,mXyzPowers[i][m+i],mCoefficients[i][m+i],
	            mSolidDerivativeXyzPowers[i][i+m][c],mSolidDerivativeCoefficients[i][i+m][c]);
	      }
	   }

	// with r^-l factor
	
	mDerivativeXyzPowers.resize(l+1);
	mDerivativeCoefficients.resize(l+1);
	for(i=mMaxL+1;i<=intL;i++)
	   {
	   	
	   mDerivativeXyzPowers[i].resize(2*l+1);
	   mDerivativeCoefficients[i].resize(2*l+1);	
	   for(m=-i;m<=i;m++)
	      {
	      powers=mXyzPowers[i][i+m];
	      
	      for(j=0;j<powers.size();j++)
	         {
	         powers[j].resize(4);
	         powers[j][3]=i;
	         }
	      	
	      mDerivativeXyzPowers[i][i+m].resize(3);
	      mDerivativeCoefficients[i][i+m].resize(3);		
	      for(c=0;c<3;c++)
	         rDifferentiate(c,powers,mCoefficients[i][m+i],
	            mDerivativeXyzPowers[i][i+m][c],mDerivativeCoefficients[i][i+m][c]);
	      }
	   }
	
		   
	updateSphSecondDerivatives(l);   	
}

//-----------------------------------------------------------------------------

void SphericalHarmonics::updateSphSecondDerivatives(
unsigned int l)
{
	int i,m,c1,c2,intL;
	intL=l;
	
	// Ylm
	
	mSolidSecondDerivativeXyzPowers.resize(l+1);
	mSolidSecondDerivativeCoefficients.resize(l+1);
	for(i=mMaxL+1;i<=intL;i++)
	   {
	   mSolidSecondDerivativeXyzPowers[i].resize(2*l+1);
	   mSolidSecondDerivativeCoefficients[i].resize(2*l+1);	
	   for(m=-i;m<=i;m++)
	      {
	      mSolidSecondDerivativeXyzPowers[i][i+m].resize(3);
	      mSolidSecondDerivativeCoefficients[i][i+m].resize(3);		
	      for(c1=0;c1<3;c1++)
	         {
	         mSolidSecondDerivativeXyzPowers[i][i+m][c1].resize(3);
	         mSolidSecondDerivativeCoefficients[i][i+m][c1].resize(3);		
	         for(c2=0;c2<=c1;c2++)
	            {
	            differentiate(c1,mSolidDerivativeXyzPowers[i][i+m][c2],
	                          mSolidDerivativeCoefficients[i][i+m][c2],
	                          mSolidSecondDerivativeXyzPowers[i][i+m][c1][c2],
	                          mSolidSecondDerivativeCoefficients[i][i+m][c1][c2]);
	            mSolidSecondDerivativeXyzPowers[i][i+m][c2][c1]=
	               mSolidSecondDerivativeXyzPowers[i][i+m][c1][c2];
	            mSolidSecondDerivativeCoefficients[i][i+m][c2][c1]=
	               mSolidSecondDerivativeCoefficients[i][i+m][c1][c2];
	                
	            }              
	         }   
	      }
	   }

// r^l * Ylm

	mSecondDerivativeXyzPowers.resize(l+1);
	mSecondDerivativeCoefficients.resize(l+1);
	for(i=mMaxL+1;i<=intL;i++)
	   {
	   mSecondDerivativeXyzPowers[i].resize(2*l+1);
	   mSecondDerivativeCoefficients[i].resize(2*l+1);	
	   for(m=-i;m<=i;m++)
	      {
	      mSecondDerivativeXyzPowers[i][i+m].resize(3);
	      mSecondDerivativeCoefficients[i][i+m].resize(3);		
	      for(c1=0;c1<3;c1++)
	         {
	         mSecondDerivativeXyzPowers[i][i+m][c1].resize(3);
	         mSecondDerivativeCoefficients[i][i+m][c1].resize(3);		
	         for(c2=0;c2<=c1;c2++)
	            {
	            rDifferentiate(c1,mDerivativeXyzPowers[i][i+m][c2],
	                       mDerivativeCoefficients[i][i+m][c2],
	                       mSecondDerivativeXyzPowers[i][i+m][c1][c2],
	                       mSecondDerivativeCoefficients[i][i+m][c1][c2]);
	            mSecondDerivativeXyzPowers[i][i+m][c2][c1]=
	               mSecondDerivativeXyzPowers[i][i+m][c1][c2];
	            mSecondDerivativeCoefficients[i][i+m][c2][c1]=
	               mSecondDerivativeCoefficients[i][i+m][c1][c2];
	            }              
	         }   
	      }
	   }
		   
	mMaxL=l;   	
}

//-----------------------------------------------------------------------------

void SphericalHarmonics::differentiate(
unsigned int component,
const vector<vector<unsigned int> > &powers,
const vector<double> &coefficients,
vector<vector<unsigned int> > &derivativePowers,
vector<double> &derivativeCoefficients)
{
	unsigned int i,n;
	derivativePowers.clear();
	derivativeCoefficients.clear();

	n=0;
	for(i=0;i<powers.size();i++)
	   {
	   if(powers[i][component]>0)
	      {
	      derivativePowers.push_back(powers[i]);
	      derivativePowers[n][component]--;
	      derivativeCoefficients.push_back(coefficients[i]*
	                                       powers[i][component]);
	      n++;	
	      }	
	   }
}

//-----------------------------------------------------------------------------

double SphericalHarmonics::numerator(
unsigned int l,
int m,
double x,
double y,
double z)
{
	if(l>mMaxL)
	   updateSph(l);

	return calculate(x,y,z,mXyzPowers[l][l+m],
	                 mCoefficients[l][l+m]);	
}

//-----------------------------------------------------------------------------

double SphericalHarmonics::numeratorDerivative(
unsigned int var,
unsigned int l,
int m,
double x,
double y,
double z,
double r)
{
	return calculate(x,y,z,mSolidDerivativeXyzPowers[l][l+m][var],
	                 mSolidDerivativeCoefficients[l][l+m][var]);	
}
	                          
//-----------------------------------------------------------------------------

double SphericalHarmonics::numeratorSecondDerivative(
unsigned int var1,
unsigned int var2,
unsigned int l,
int m,
double x,
double y,
double z,
double r)
{
	unsigned int c1,c2;
	
	c1=max(var1,var2);
	c2=min(var1,var2);	
	return calculate(x,y,z,mSolidSecondDerivativeXyzPowers[l][l+m][c1][c2],
	                 mSolidSecondDerivativeCoefficients[l][l+m][c1][c2]);
}

//-----------------------------------------------------------------------------

double SphericalHarmonics::calculate(
double x,
double y,
double z,
const vector<vector<unsigned int> > &xyzPowers,
const vector<double> &coefficients)
{
	unsigned int i,n;
	double result;
	n=coefficients.size();
	
	result=0.0;
	for(i=0;i<n;i++)
	   result+=coefficients[i]*pow(x,xyzPowers[i][0])*pow(y,xyzPowers[i][1])*
	           pow(z,xyzPowers[i][2]);
	return result;
}

//-----------------------------------------------------------------------------

void SphericalHarmonics::calculate(
unsigned int maxL,
double x,
double y,
double z,
std::vector<std::vector<double> > &values,
SphericalHarmonicForm form)
{
	unsigned int m_index,n,i,modM;
	unsigned int l;
	double r;

	r=sqrt(x*x+y*y+z*z);

	if(!initialized)
	   initialize();
	if(maxL>mMaxL)
	   updateSph(maxL);
	
	mPowerXValue[0]=1;
	mPowerYValue[0]=1;
	mPowerZValue[0]=1;
	mPowerRValue[0]=1;	
	
	if(form==DensityNormalized||form==SolidDensityNormalized)
	   values[0][0]=1/(4.0*M_PI);
	else
	   values[0][0]=sqrt(1/(4.0*M_PI));		
	
	for(l=1;l<=maxL;l++)
	   {	
	   mPowerXValue[l]=mPowerXValue[l-1]*x;
	   mPowerYValue[l]=mPowerYValue[l-1]*y;
	   mPowerZValue[l]=mPowerZValue[l-1]*z;		   
	   mPowerRValue[l]=mPowerRValue[l-1]*r;
	   	
	   for(m_index=0;m_index<2*l+1;m_index++)
	      {
	      n=mCoefficients[l][m_index].size();
	      values[l][m_index]=0;
	       
	      for(i=0;i<n;i++)
	         values[l][m_index]+=mCoefficients[l][m_index][i]*
	                             mPowerXValue[mXyzPowers[l][m_index][i][0]]*
	                             mPowerYValue[mXyzPowers[l][m_index][i][1]]*
	                             mPowerZValue[mXyzPowers[l][m_index][i][2]]; 
	                             
	      if(form==DensityNormalized||form==SolidDensityNormalized)
	         {
	         m_index>l ? modM=m_index-l : modM=l-m_index;
	         values[l][m_index]*=mFlm[l][modM];
	         }
	      if(form==WaveFunctionNormalized||form==DensityNormalized)
	         if(mPowerRValue[l]>numeric_limits<double>::min())
	            values[l][m_index]/=mPowerRValue[l];
	         else
	            values[l][m_index]=0.0;
	      }
	   }
}

//-----------------------------------------------------------------------------

void SphericalHarmonics::get(
unsigned int l,
int m,
vector<double> &coefficients,
vector<vector<unsigned int> > &powers,
SphericalHarmonicForm form)
{
	unsigned int i;
	if(!initialized)
	   initialize();

	if(l>mMaxL)
	   updateSph(l);   
	   
	coefficients=mCoefficients[l][l+m];
	powers=mXyzPowers[l][l+m];
	
	if(form==DensityNormalized||form==SolidDensityNormalized)
	   for(i=0;i<coefficients.size();i++)
	      coefficients[i]*=mFlm[l][abs(m)];
}

//-----------------------------------------------------------------------------

void SphericalHarmonics::get(
unsigned int maxL,
vector<vector<vector<double> > > &coefficients,
vector<vector<vector<vector<unsigned int> > > > &powers,
SphericalHarmonicForm form)
{
	unsigned int l;
	int m,intL;

	vector<vector<unsigned int> > auxPowers;
	if(!initialized)
	   initialize();
	   

	if(maxL>mMaxL)
	   updateSph(maxL);   
	   
	coefficients.resize(maxL+1);
	powers.resize(maxL+1);
	for(l=0;l<=maxL;l++)
	   {
	   coefficients[l].resize(2*l+1);
	   powers[l].resize(2*l+1);
	   intL=l;
	   for(m=-intL;m<=intL;m++)
	      get(l,m,coefficients[l][l+m],powers[l][l+m],form);
	   }
	
}

//-----------------------------------------------------------------------------

void SphericalHarmonics::getRealSolidHarmonics(
unsigned int l,
vector<vector<vector<unsigned int> > > &powers,
vector<vector<double> > &coefficients)  
{
	unsigned int t,u,v,vm,modM;	
	double coefficient,multiplier,sign;
	int m,iL;
	
	vector<unsigned int> xyzPowers(3);
	double fac1,fac2,fac3,bin1,bin2,bin3,bin4,powerOfTwo,powerOfFour;
	double n,x;
	
	multiplier=sqrt((2*l+1.0)/(4.0*M_PI));
	
		
	
	iL=l;	

	
	powers.resize(2*l+1); 

	coefficients.resize(2*l+1); 
	
			
	for(m=-iL;m<=iL;m++)
	   {
	   modM=abs(m);	

	   fac1=discamb::math_utilities::factorial<double>(l+modM);
	   fac2= discamb::math_utilities::factorial<double>(l-modM);
	   fac3= discamb::math_utilities::factorial<double>(l);
	   powerOfTwo=pow(2.0,modM);
	   
	   n=fac1*fac2*(2-(0==m?1:0));
	   n/=fac3*powerOfTwo*fac3*powerOfTwo;
	
	   for(t=0;t<=(l-modM)/2;t++)
	      for(u=0;u<=t;u++)
	         {
	         m>=0?vm=0:vm=1;	
	         for(v=0;v<=(modM-vm)/2;v++)
	            {
	            xyzPowers[0]=2*t+modM-2*(u+v)-vm;
	            xyzPowers[1]=2*(u+v)+vm;
	            xyzPowers[2]=l-2*t-modM;
	            powerOfFour=pow(4.0,t);
	            bin1=discamb::math_utilities::binomialCoefficient<double>(l,t);
	            bin2= discamb::math_utilities::binomialCoefficient<double>(l-t,modM+t);
	            bin3= discamb::math_utilities::binomialCoefficient<double>(t,u);
	            bin4= discamb::math_utilities::binomialCoefficient<double>(modM,2*v+vm);
	            
	            x=(0==(t+v)%2?1:-1)*bin1*bin2*bin3*bin4;
	            x/=powerOfFour;
	            
				x>0 ? sign=1 : sign=-1;
	            x*=x;
	            x*=n;
	            coefficient=x;
	            coefficient=sign*sqrt(coefficient);	            	            	            
	            powers[l+m].push_back(xyzPowers);
	            coefficients[l+m].push_back(coefficient*multiplier);      
	            }	
	         }   
	   }
}

//-----------------------------------------------------------------------------

void SphericalHarmonics::getDerivative(
unsigned int component,
unsigned int l,
int m,
vector<double> &coefficients,
vector<vector<unsigned int> > &powers,
SphericalHarmonicForm form)
{
	vector<double>::size_type i;
	if(!initialized)
	   initialize();

	if(l>mMaxL)
	   updateSph(l);   	
	
	if(form==WaveFunctionNormalized||form==DensityNormalized)
	   coefficients=mDerivativeCoefficients[l][l+m][component];
	else
	   coefficients=mSolidDerivativeCoefficients[l][l+m][component];
	
	if(form==DensityNormalized||form==SolidDensityNormalized)
	   {
	   for(i=0;i<coefficients.size();i++)
	      coefficients[i]*=mFlm[l][abs(m)];
	   }

	if(form==WaveFunctionNormalized||form==DensityNormalized)       
	   powers=mDerivativeXyzPowers[l][l+m][component];
	else
	   powers=mSolidDerivativeXyzPowers[l][l+m][component];
}
	                	                
//-----------------------------------------------------------------------------

void SphericalHarmonics::getSecondDerivative(
unsigned int component1,
unsigned int component2,
unsigned int l,
int m,
vector<double> &coefficients,
vector<vector<unsigned int> > &powers,
SphericalHarmonicForm form)
{
	unsigned int i;
	if(!initialized)
	   initialize();

	if(l>mMaxL)
	   updateSph(l);   
	
	if(form==WaveFunctionNormalized||form==DensityNormalized)
	   coefficients=mSecondDerivativeCoefficients[l][l+m][component1][component2];
	else
	   coefficients=mSolidSecondDerivativeCoefficients[l][l+m][component1][component2];
	
	if(form==DensityNormalized||form==SolidDensityNormalized)
	   {
	   for(i=0;i<coefficients.size();i++)
	      coefficients[i]*=mFlm[l][abs(m)];
	   }

	if(form==WaveFunctionNormalized||form==DensityNormalized)       
	   powers=mSecondDerivativeXyzPowers[l][l+m][component1][component2];
	else
	   powers=mSolidSecondDerivativeXyzPowers[l][l+m][component1][component2];
}

//-----------------------------------------------------------------------------

void SphericalHarmonics::rDifferentiate(
unsigned int component,
const std::vector<std::vector<unsigned int> > &powers,
const std::vector<double> &coefficients,
std::vector<std::vector<unsigned int> > &derivativePowers,
std::vector<double> &derivativeCoefficients)
{
	unsigned int i,n;
	derivativePowers.clear();
	derivativeCoefficients.clear();


	n=0;
	for(i=0;i<powers.size();i++)
	   {
	   if(powers[i][component]>0)
	      {
	      derivativePowers.push_back(powers[i]);
	      derivativePowers[n][component]--;
	      derivativeCoefficients.push_back(coefficients[i]*powers[i][component]);
	      n++;	
	      }
	   if(powers[i][3]>0)
	      {
	      derivativePowers.push_back(powers[i]);
	      derivativePowers[n][component]++;
	      derivativePowers[n][3]+=2;
	      derivativeCoefficients.push_back(-coefficients[i]*powers[i][3]);
	      n++;
	      }
	   }
}

//-----------------------------------------------------------------------------

double SphericalHarmonics::getConversionFactor(
unsigned int l,
int m)
{
	if(!initialized)
	   initialize();
	return mFlm[l][abs(m)];
}

//-----------------------------------------------------------------------------
/*
void SphericalHarmonics::getConversionMatrix(
unsigned int l,
const vector<vector<double> > &oldSystemOfCoordinates,
const vector<vector<double> > &newSystemOfCoordinates,
vector<vector<double> > &conversionMatrix)
{
	vector<double> coefficients1,coefficients2;
	vector<vector<unsigned int> > powers1,powers2;
	int m1,m2,intL;
	Polynomial<double,unsigned int> polynomial1,polynomial2,p;

	intL=l;

	conversionMatrix.assign(2*l+1,vector<double>(2*l+1,0.0));
	
	
	for(m1=-intL;m1<=intL;m1++)
		for(m2=-intL;m2<=intL;m2++)
		   {
	       get(l,m1,coefficients1,powers1,WaveFunctionNormalized);
		   get(l,m2,coefficients2,powers2,WaveFunctionNormalized);

	       makePolynomial(coefficients1,powers1,polynomial1,newSystemOfCoordinates);
	       makePolynomial(coefficients2,powers2,polynomial2,oldSystemOfCoordinates);
		   p=polynomial1*polynomial2;
	       conversionMatrix[m1+intL][m2+intL]=integrate(p);
		   }
}*/

void SphericalHarmonics::getConversionMatrix(
unsigned int l,
const vector<vector<double> > &oldSystemOfCoordinates,
const vector<vector<double> > &newSystemOfCoordinates,
vector<vector<double> > &conversionMatrix)
{
	vector<double> coefficients1,coefficients2;
	vector<vector<unsigned int> > powers1,powers2;
	int m1,m2,intL;
	Polynomial<double,unsigned int> polynomial1,polynomial2,p;

	intL=l;

	conversionMatrix.assign(2*l+1,vector<double>(2*l+1,0.0));
	
	
	for(m1=-intL;m1<=intL;m1++)
		for(m2=-intL;m2<=intL;m2++)
		   {
	       get(l,m1,coefficients1,powers1,WaveFunctionNormalized);
		   get(l,m2,coefficients2,powers2,WaveFunctionNormalized);

	       makePolynomial(coefficients1,powers1,polynomial1,newSystemOfCoordinates);
	       makePolynomial(coefficients2,powers2,polynomial2,oldSystemOfCoordinates);
		   p=polynomial1*polynomial2;
	       conversionMatrix[m1+intL][m2+intL]=integrate(p);
		   }
}

void SphericalHarmonics::getConversionMatrices(
unsigned int maxL,
const vector<vector<double> > &oldSystemOfCoordinates,
const vector<vector<double> > &newSystemOfCoordinates,
vector<vector<vector<double> > > &conversionMatrices)
{
	vector<double> coefficients1,coefficients2;
	vector<vector<unsigned int> > powers1,powers2;
	unsigned int l,m1,m2;
	// [l][m]
	vector<vector<Polynomial<double,unsigned int> > > oldPolynomials,newPolynomials;
	
	getPolynomials(maxL,oldSystemOfCoordinates,newSystemOfCoordinates,
				   oldPolynomials,newPolynomials);

	conversionMatrices.resize(maxL+1);
	for(l=0;l<=maxL;l++)
		conversionMatrices[l].assign(2*l+1,vector<double>(2*l+1,0.0));

	
	
	for(l=0;l<=maxL;l++)
	{
		for(m1=0;m1<2*l+1;m1++)
			for(m2=0;m2<2*l+1;m2++)
				conversionMatrices[l][m1][m2]=
					    integrate(newPolynomials[l][m1]*oldPolynomials[l][m2]);
	}
}


double SphericalHarmonics::maxAbs(
unsigned int l,
int m,
SphericalHarmonicForm form)
{
	int intL=l,n;
	n=(intL+1)*(intL+1)-1+intL+m;
	if(form==WaveFunctionNormalized)
		return extremumAbsoluteValue[n]/getConversionFactor(l,m);
	if(form==DensityNormalized)
		return extremumAbsoluteValue[n];
	return 0.0;
}

void SphericalHarmonics::getPolynomials(
unsigned int maxL,
const vector<vector<double> > &oldSystemOfCoordinates,
const vector<vector<double> > &newSystemOfCoordinates,
vector<vector<Polynomial<double,unsigned int> > > &oldPolynomials,
vector<vector<Polynomial<double,unsigned int> > > &newPolynomials)
{
	vector<Polynomial<double,unsigned int> > powersX,powersY,powersZ;
	Polynomial<double,unsigned int> one;
	vector<vector<double> > coordinates(3,vector<double>(3));
	vector<vector<unsigned int> > xyzPowers(3,vector<unsigned int>(3,0));
	unsigned int i,j,l,n,k;
	int intL,m;

	if(mMaxL<maxL)
		updateSph(maxL);

	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
		{
			coordinates[i][j]=0.0;
			for(k=0;k<3;k++)
				coordinates[i][j]+=newSystemOfCoordinates[i][k]*
				                  oldSystemOfCoordinates[j][k];
		}


	xyzPowers[0][0]=1;
	xyzPowers[1][1]=1;
	xyzPowers[2][2]=1;

	powersX.resize(max(maxL+1,2U));
	powersY.resize(max(maxL+1,2U));
	powersZ.resize(max(maxL+1,2U));
	powersX[0].addTerm(1.0,vector<unsigned int>(3,0));
	powersY[0].addTerm(1.0,vector<unsigned int>(3,0));
	powersZ[0].addTerm(1.0,vector<unsigned int>(3,0));
	powersX[1].addTerms(coordinates[0],xyzPowers);
	powersY[1].addTerms(coordinates[1],xyzPowers);
	powersZ[1].addTerms(coordinates[2],xyzPowers);

	for(l=2;l<=maxL;l++)
	{
		powersX[l]=powersX[l-1]*powersX[1];
		powersY[l]=powersY[l-1]*powersY[1];
		powersZ[l]=powersZ[l-1]*powersZ[1];
	}

	oldPolynomials.resize(maxL+1);
	newPolynomials.resize(maxL+1);
	for(l=0;l<=maxL;l++)
	{
		oldPolynomials[l].resize(2*l+1);
		newPolynomials[l].resize(2*l+1);
		intL=l;
		for(m=-intL;m<=intL;m++)
		{
			oldPolynomials[l][intL+m].addTerms(mCoefficients[l][intL+m],mXyzPowers[l][intL+m]);
			n=mCoefficients[l][intL+m].size();
			for(i=0;i<n;i++)
				newPolynomials[l][intL+m]+=mCoefficients[l][intL+m][i]*
					                       powersX[mXyzPowers[l][intL+m][i][0]]*
										   powersY[mXyzPowers[l][intL+m][i][1]]*
										   powersZ[mXyzPowers[l][intL+m][i][2]];
		}
	}


}




namespace{
void makePolynomial(
const vector<double> &coefficients,
const vector<vector<unsigned int> > &powers,
Polynomial<double,unsigned int> &polynomial,
const vector<vector<double> > &coordinates)
{
	Polynomial<double,unsigned int> x,y,z,one;
	unsigned int i;
	vector<vector<unsigned int> > xyzPowers(3,vector<unsigned int>(3,0));
	xyzPowers[0][0]=1;
	xyzPowers[1][1]=1;
	xyzPowers[2][2]=1;

	x.addTerms(coordinates[0],xyzPowers);
	y.addTerms(coordinates[1],xyzPowers);
	z.addTerms(coordinates[2],xyzPowers);
	one.addTerm(1.0,vector<unsigned int>(3,0));

	polynomial=0;
	for(i=0;i<powers.size();i++)
	   polynomial+=coefficients[i]*discamb::math_utilities::power(x,powers[i][0])* discamb::math_utilities::power(y,powers[i][1])
	                             * discamb::math_utilities::power(z,powers[i][2])*one;

}

double integrate(
const Polynomial<double,unsigned int> &p)
{
	vector<double> coefficients;
	vector<vector<unsigned int> > powers;
	double result;
	unsigned int i,l;
	p.get(coefficients,powers);
	l=0;
	if(powers.size()>0)
		l=powers[0][0]+powers[0][1]+powers[0][2];

	result=0;
	for(i=0;i<powers.size();i++)
		if(powers[i][0]%2==0&&powers[i][1]%2==0&&powers[i][2]%2==0)
			result+=4.0*M_PI*coefficients[i]*
			        discamb::math_utilities::doubleFactorial<double>(powers[i][0]-1)*
            discamb::math_utilities::doubleFactorial<double>(powers[i][1]-1)*
            discamb::math_utilities::doubleFactorial<double>(powers[i][2]-1)/
            discamb::math_utilities::doubleFactorial<double>(l+1);
	return result;
}

}

}

