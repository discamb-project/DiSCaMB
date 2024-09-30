#ifndef SPHERICALHARMONICS_H_
#define SPHERICALHARMONICS_H_

#include <vector>
#include "Polynomial.h"

namespace MChLib{

// internally represented by Clm*numerator*rPowerMinusL (denominator=1/r^l)


/**
\brief A class for evaluating values and formulas of real spherical harmonics and 
solid real spherical harmonics expressed in cartesian coordinates.

<b>Definitions</b>

The real spherical harmonics are expressed with use of complex
 spherical harmonics in the following way:

For m=0

\f[
R_{lm}=Y_{lm}
\f]
For m>0

\f[
R_{lm}=\frac{(-1)^m}{\sqrt{2}}(Y_{l,m}+Y_{l,-m})
\f]

\f[
R_{l,-m}=\frac{(-1)^m}{\mathrm{i}\sqrt{2}}(Y_{l,m}-Y_{l,-m})
\f]

 
 Where the complex spherical harmonics are given by:
 \f$Y_{lm}(\theta,\phi)\f$:


\f[
Y_{lm}(\theta,\phi)=(-1)^m\sqrt{\frac{(2l+1)(l-|m|)!}{4\pi(l+|m|)!}}
P^{|m|}_l(\cos\theta)e^{im\phi}
\f]

and Legendre polynomials \f$P^{m}_l\f$ as:

\f[
P^{m}_l(x)=\frac{1}{2^ll!}(1-x^2)^{m/2} \frac{\mathrm{d}^{l+m}}{\mathrm{d}x^{l+m}}
(x^2-1)^l
\f]

The solid real spherical harmonics \f$S_{lm}\f$ are given by:
\f[
S_{lm}=r^lR_{lm}
\f]

 */
class SphericalHarmonics
{
public:
	SphericalHarmonics();
	virtual ~SphericalHarmonics();
	enum SphericalHarmonicForm {WaveFunctionNormalized,DensityNormalized,
	                  SolidWaveFunctionNormalized,SolidDensityNormalized};
	 
	/** returns value of real spherical harmonics
	 possible normalization codes
	 <DL>
	 <DT><Strong> WaveFunctionNormalized </Strong>
	 <DD> returns \f$R_{lm}\f$	 
	 <DT><Strong> DensityNormalized</Strong>
	 <DD> returns \f$NR_{lm}\f$
	      with \f$N=(2-\delta_{l0})/\int |R_{lm}|\mathrm{d}\Omega \f$
	 <DT><Strong> SolidWaveFunctionNormalized</Strong>
	 <DD> returns \f$r^lR_{lm}\f$
	 <DT><Strong> SolidDensityNormalized</Strong>
	 <DD> returns \f$r^lNR_{lm}\f$
	 </DL> 
	 */
	 
	static double calculate(unsigned int l,int m, double x,double y,
	                        double z,SphericalHarmonicForm=WaveFunctionNormalized);	 
	/** Like the function above, but takes also r as an additional
	argument (instead of calculating it from x,y,z).
	*/ 
	static double calculate(unsigned int l,int m, double x,double y,
	                        double z, double r,
	                        SphericalHarmonicForm=WaveFunctionNormalized);

	static void calculate(unsigned int maxL,double x,double y,double z,
	                        std::vector<std::vector<double> > &values,
	                        SphericalHarmonicForm=WaveFunctionNormalized);    

	static void calculateGradient(unsigned int l,int m, double x,double y,
	                        double z, double r,double (&gradient)[3],
	                        SphericalHarmonicForm=WaveFunctionNormalized);                        
	                        
	static void calculateHessian(unsigned int l,int m, double x,double y,
	                        double z, double r,double (&hessian)[3][3],
	                        SphericalHarmonicForm=WaveFunctionNormalized); 
	
	static void get(unsigned int l,int m,std::vector<double> &coefficients,
	                std::vector<std::vector<unsigned int> > &powers,
	                SphericalHarmonicForm=WaveFunctionNormalized);
	                
	static void getDerivative(unsigned int component,unsigned int l,int m,
	                std::vector<double> &coefficients,
	                std::vector<std::vector<unsigned int> > &powers,
	                SphericalHarmonicForm=WaveFunctionNormalized);
	                
	static void getSecondDerivative(unsigned int component1,
	               unsigned int component2,unsigned int l,int m,
	               std::vector<double> &coefficients,
	               std::vector<std::vector<unsigned int> > &powers,
	               SphericalHarmonicForm=WaveFunctionNormalized);

	static void get(unsigned int maxL,
	    std::vector<std::vector<std::vector<double> > > &coefficients,
	    std::vector<std::vector<std::vector<std::vector<unsigned int> > > > &powers,
	    SphericalHarmonicForm=WaveFunctionNormalized);
	/// returns c such that 
	/// DensityNormalizedSph=c*WaveFunctionNormalizedSph
	static double getConversionFactor(unsigned int l,int m);

	static void getConversionMatrix(unsigned int l,
		     const std::vector<std::vector<double> > &oldSystemOfCoordinates,
			 const std::vector<std::vector<double> > &newSystemOfCoordinates,
			 std::vector<std::vector<double> > &conversionMatrix);
	static double maxAbs(unsigned int l,int m,
		                 SphericalHarmonicForm=WaveFunctionNormalized);
	static void getConversionMatrices(unsigned int maxL,
				const std::vector<std::vector<double> > &oldSystemOfCoordinates,
				const std::vector<std::vector<double> > &newSystemOfCoordinates,
				std::vector<std::vector<std::vector<double> > > &conversionMatrices);
	
private:


	//      DATA

	
	static std::vector<std::vector<double> > mFlm;	
	
	static std::vector<std::vector<std::vector<double> > > mCoefficients;	
	
	static std::vector<std::vector<std::vector<std::vector<unsigned int> > > >
	                                                                mXyzPowers;
	                                                                   
	static std::vector<std::vector<std::vector<std::vector<double> > > > 
	                                                   mDerivativeCoefficients;
	                                                   
	static std::vector<std::vector<std::vector<std::vector
	                   <std::vector<unsigned int> > > > > mDerivativeXyzPowers;
	                   
	static std::vector<std::vector<std::vector<std::vector
	                 <std::vector<double> > > > > mSecondDerivativeCoefficients;	
	                 
	static std::vector<std::vector<std::vector<std::vector<std::vector<
	             std::vector<unsigned int> > > > > > mSecondDerivativeXyzPowers;

	static std::vector<std::vector<std::vector<std::vector<double> > > > 
	                                               mSolidDerivativeCoefficients;
	                                               
	static std::vector<std::vector<std::vector<std::vector<
	                std::vector<unsigned int> > > > > mSolidDerivativeXyzPowers;
	                
	static std::vector<std::vector<std::vector<std::vector<
	             std::vector<double> > > > > mSolidSecondDerivativeCoefficients;
	             
	static std::vector<std::vector<std::vector<std::vector<std::vector<
	        std::vector<unsigned int> > > > > > mSolidSecondDerivativeXyzPowers;
	        
	static std::vector<double> mPowerXValue;
	
	static std::vector<double> mPowerYValue;
	
	static std::vector<double> mPowerZValue;
	
	static std::vector<double> mPowerRValue;	
	
	static bool initialized;
	
	static unsigned int mMaxL;
	
	
	//         METHODS
		
	
	static void initialize();

	static void initializeSph();


	static double higherSpH(unsigned int l,int m,double x,double y,
	                       double z,double r);

	static double numerator(unsigned int l,int m, double x,double y,
	                        double z);
		      
	                          
	static double rPowerMinusL(unsigned int l,double r);

	static double rPowerMinusLDerivative(unsigned int var,unsigned int l,
	                          double x,double y,double z, double r);	   

	static double rPowerMinusLSecondDerivative(unsigned int var1,unsigned int var2,
	                          unsigned int l,double x,double y,double z,
	                          double r);	   

	static double numeratorDerivative(unsigned int var,unsigned int l,
	                          int m, double x,double y,double z, double r);

	static double numeratorSecondDerivative(unsigned int var1,
	                          unsigned int var2,unsigned int l,int m,
	                          double x,double y,double z, double r);

	static void updateSph(unsigned int l);

	static void updateSphDerivatives(unsigned int l);

	static void updateSphSecondDerivatives(unsigned int l);	

	static double calculate(double x,double y,double z,
	       const std::vector<std::vector<unsigned int> > &xyzPowers,
	       const std::vector<double> &coefficients);

	static void differentiate(
	               unsigned int component,
	               const std::vector<std::vector<unsigned int> > &powers,
	               const std::vector<double> &coefficients,
                   std::vector<std::vector<unsigned int> > &derivativePowers,
                   std::vector<double> &derivativeCoefficients);

	static void rDifferentiate(
	               unsigned int component,
	               const std::vector<std::vector<unsigned int> > &powers,
	               const std::vector<double> &coefficients,
                   std::vector<std::vector<unsigned int> > &derivativePowers,
                   std::vector<double> &derivativeCoefficients);

	static void getRealSolidHarmonics(unsigned int l,
	          std::vector<std::vector<std::vector<unsigned int> > > &xyzPowers,
	          std::vector<std::vector<double> > &coefficients);  
	
	static void getPolynomials(unsigned int maxL,
		const std::vector<std::vector<double> > &oldSystemOfCoordinates,
	    const std::vector<std::vector<double> > &newSystemOfCoordinates,
		std::vector<std::vector<Polynomial<double,unsigned int> > > &oldPolynomials,
		std::vector<std::vector<Polynomial<double,unsigned int> > > &newPolynomials);
							   

};

}

#endif /*REALSPHERICALHARMONIC_H_*/
