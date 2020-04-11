#include "FirstOrderDerivedDensityMomentumCorrLine.hpp"
#include "../Integration/Integral_OSC_1D.hpp"
#include <iostream>
#include <gsl/gsl_math.h>


namespace KFT
{

	// ***********************************************************
	// ***                protected members                    ***
	// ***********************************************************
    
	double FirstOrderDerivedDensityMomentumCorrLine::Kernel_n1(const double K_ij_, const double q_, const bool Direction_) const
	{	
		double L2 = fabs(funcArgsRef.mu_ij * L2_ij());
		double Q = Q_ij(q_);
		double A4 = 6* funcArgsRef.L_k*funcArgsRef.L_j*M_kj()*iniCorrRef.a_02(q_)*iniCorrRef.get_sigma_p2();
		double A1 = (Q_kj(q_)*funcArgsRef.L_j*funcArgsRef.mu_j-funcArgsRef.mu_k*funcArgsRef.L_k);
		double A2 =6*iniCorrRef.get_sigma_p2()*funcArgsRef.L_j*funcArgsRef.L_j*funcArgsRef.L_k*M_kj()*iniCorrRef.a_02(q_)*funcArgsRef.mu_j;
		if(Direction_ == true)
		{
			if((fabs(L2) < 1.0e-3))
			{
				return  (1-L2)*exp(-Q)*q_ * q_ * iniCorrRef.a_11(q_) * funcArgsRef.mu_j * funcArgsRef.L_j  * 
						(Q_kj(q_)- 21/(funcArgsRef.K_ij*q_)/(funcArgsRef.K_ij*q_) * A4*B_ij(q_)*B_ij(q_)
						- 2*B_ij(q_)*A4/(funcArgsRef.K_ij*q_) - B_ij(q_)*B_ij(q_)*Q_kj(q_));
			}
			
			else if( (exp(-L2) > 0.0) && (fabs(Q) < 1.0e-4))
			{
				return   (1-Q)*exp(-L2)*q_ * q_ * iniCorrRef.a_11(q_) * funcArgsRef.mu_j * funcArgsRef.L_j  * 
						(Q_kj(q_)- 21/(funcArgsRef.K_ij*q_)/(funcArgsRef.K_ij*q_) * A4*B_ij(q_)*B_ij(q_)
						- 2*B_ij(q_)*A4/(funcArgsRef.K_ij*q_) - B_ij(q_)*B_ij(q_)*Q_kj(q_));
			}
		
			else
			{
				return  exp(-L2-Q)*q_ * q_ * iniCorrRef.a_11(q_) * funcArgsRef.mu_j * funcArgsRef.L_j  * 
						(Q_kj(q_)- 21/(funcArgsRef.K_ij*q_)/(funcArgsRef.K_ij*q_) * A4*B_ij(q_)*B_ij(q_)
						- 2*B_ij(q_)*A4/(funcArgsRef.K_ij*q_) - B_ij(q_)*B_ij(q_)*Q_kj(q_));
			}
	   }
	   
	   if(Direction_ == false)
	   {	
		if((fabs(L2) < 1.0e-3))
		{
			return  (1-L2)*q_ * q_ * iniCorrRef.a_11(q_) * 
					(funcArgsRef.mu_k*funcArgsRef.L_k+ exp(-Q)*(A1- 2*A2*B_ij(q_)/(funcArgsRef.K_ij*q_) - B_ij(q_)*B_ij(q_)*iniCorrRef.a_11(q_)
					- 21*A2*B_ij(q_)*B_ij(q_)/(funcArgsRef.K_ij*q_)/(funcArgsRef.K_ij*q_)));
		}
			
		else if( (exp(-L2) > 0.0) && (fabs(Q) < 1.0e-4))
		{
			return  exp(-L2)*q_ * q_ * iniCorrRef.a_11(q_) * 
					(funcArgsRef.mu_k*funcArgsRef.L_k+ (1-Q)*(A1- 2*A2*B_ij(q_)/(funcArgsRef.K_ij*q_) - B_ij(q_)*B_ij(q_)*iniCorrRef.a_11(q_)
					- 21*A2*B_ij(q_)*B_ij(q_)/(funcArgsRef.K_ij*q_)/(funcArgsRef.K_ij*q_)));
		}
	
		else
		{
			return  exp(-L2)*q_ * q_ * iniCorrRef.a_11(q_) * 
					(funcArgsRef.mu_k*funcArgsRef.L_k+ exp(-Q)*(A1- 2*A2*B_ij(q_)/(funcArgsRef.K_ij*q_) - B_ij(q_)*B_ij(q_)*iniCorrRef.a_11(q_)
					- 21*A2*B_ij(q_)*B_ij(q_)/(funcArgsRef.K_ij*q_)/(funcArgsRef.K_ij*q_)));
		}

	   }
	   return 0.0;
	}

	double FirstOrderDerivedDensityMomentumCorrLine::Kernel_n2(const double K_ij_, const double q_, const bool Direction_) const
	{	
		
		double L2 = fabs(funcArgsRef.mu_ij * L2_ij());
		double Q = Q_ij(q_);
		double A4 = 6* funcArgsRef.L_k*funcArgsRef.L_j*M_kj()*iniCorrRef.a_02(q_)*iniCorrRef.get_sigma_p2();
		double A1 = (Q_kj(q_)*funcArgsRef.L_j*funcArgsRef.mu_j-funcArgsRef.mu_k*funcArgsRef.L_k);
		double A2 =6*iniCorrRef.get_sigma_p2()*funcArgsRef.L_j*funcArgsRef.L_j*funcArgsRef.L_k*M_kj()*iniCorrRef.a_02(q_)*funcArgsRef.mu_j;

		if(Direction_== true)
		{
			if((fabs(L2) < 1.0e-3))
			{
				return q_* q_ *  (1-L2)*exp(-Q) * iniCorrRef.a_11(q_) * funcArgsRef.mu_j * funcArgsRef.L_j 
				     	*(A4/(funcArgsRef.K_ij*q_) + B_ij(q_)*Q_kj(q_) + 10*B_ij(q_)*A4/(funcArgsRef.K_ij*q_)/(funcArgsRef.K_ij*q_)
						 + B_ij(q_)*B_ij(q_)*Q_kj(q_)*5/(funcArgsRef.K_ij*q_)+3*B_ij(q_)*B_ij(q_)*
						 	A4/(funcArgsRef.K_ij*q_)*(35/(funcArgsRef.K_ij*q_)/(funcArgsRef.K_ij*q_) -1)  );
			}
			
			else if( (exp(-L2) > 0.0) && (fabs(Q) < 1.0e-4))
			{ 
				return q_* q_ * (1-Q)* exp(-L2) * iniCorrRef.a_11(q_) * funcArgsRef.mu_j * funcArgsRef.L_j 
				     	*(A4/(funcArgsRef.K_ij*q_) + B_ij(q_)*Q_kj(q_) + 10*B_ij(q_)*A4/(funcArgsRef.K_ij*q_)/(funcArgsRef.K_ij*q_)
						 + B_ij(q_)*B_ij(q_)*Q_kj(q_)*5/(funcArgsRef.K_ij*q_)+3*B_ij(q_)*B_ij(q_)*
						 	A4/(funcArgsRef.K_ij*q_)*(35/(funcArgsRef.K_ij*q_)/(funcArgsRef.K_ij*q_) -1)  );
			}
		
			else
			{
				return q_* q_ *  exp(-L2-Q) * iniCorrRef.a_11(q_) * funcArgsRef.mu_j * funcArgsRef.L_j 
				     	*(A4/(funcArgsRef.K_ij*q_) + B_ij(q_)*Q_kj(q_) + 10*B_ij(q_)*A4/(funcArgsRef.K_ij*q_)/(funcArgsRef.K_ij*q_)
						 + B_ij(q_)*B_ij(q_)*Q_kj(q_)*5/(funcArgsRef.K_ij*q_)+3*B_ij(q_)*B_ij(q_)*
						 	A4/(funcArgsRef.K_ij*q_)*(35/(funcArgsRef.K_ij*q_)/(funcArgsRef.K_ij*q_) -1)  );
			}
		}
		if(Direction_ == false)
		{
			if((fabs(L2) < 1.0e-3))
			{
				return q_* q_ *  (1-L2)*exp(-Q) * iniCorrRef.a_11(q_) *
						(A2/(funcArgsRef.K_ij*q_)+B_ij(q_)*A1+5/(funcArgsRef.K_ij*q_)*(21*A2*B_ij(q_)*B_ij(q_)/(funcArgsRef.K_ij*q_)/(funcArgsRef.K_ij*q_)
						+2*B_ij(q_)*A2/(funcArgsRef.K_ij*q_)+A1*B_ij(q_)*B_ij(q_)- 3 *A2 *B_ij(q_)*B_ij(q_)/(funcArgsRef.K_ij*q_) ));
			}
			
			else if( (exp(-L2) > 0.0) && (fabs(Q) < 1.0e-4))
			{ 
				return q_* q_ * (1-Q)* exp(-L2) * iniCorrRef.a_11(q_) *
						(A2/(funcArgsRef.K_ij*q_)+B_ij(q_)*A1+5/(funcArgsRef.K_ij*q_)*(21*A2*B_ij(q_)*B_ij(q_)/(funcArgsRef.K_ij*q_)/(funcArgsRef.K_ij*q_)
						+2*B_ij(q_)*A2/(funcArgsRef.K_ij*q_)+A1*B_ij(q_)*B_ij(q_)- 3 *A2 *B_ij(q_)*B_ij(q_)/(funcArgsRef.K_ij*q_) ));
			}
		
			else
			{
				return q_* q_ *  exp(-L2-Q) * iniCorrRef.a_11(q_) *
						(A2/(funcArgsRef.K_ij*q_)+B_ij(q_)*A1+5/(funcArgsRef.K_ij*q_)*(21*A2*B_ij(q_)*B_ij(q_)/(funcArgsRef.K_ij*q_)/(funcArgsRef.K_ij*q_)
						+2*B_ij(q_)*A2/(funcArgsRef.K_ij*q_)+A1*B_ij(q_)*B_ij(q_)- 3 *A2 *B_ij(q_)*B_ij(q_)/(funcArgsRef.K_ij*q_) ));
			}
		}
		return 0.0;
	   
	}	
	double FirstOrderDerivedDensityMomentumCorrLine::Kernel_n_plus_1(const double K_ij_, const double q_, const std::size_t n_, const bool Direction_) const
	{
		double L2 = fabs(funcArgsRef.mu_ij * L2_ij());
		double Q = Q_ij(q_);
		double A4 = 6* funcArgsRef.L_k*funcArgsRef.L_j*M_kj()*iniCorrRef.a_02(q_)*iniCorrRef.get_sigma_p2();
		double A1 = (Q_kj(q_)*funcArgsRef.L_j*funcArgsRef.mu_j-funcArgsRef.mu_k*funcArgsRef.L_k);
		double A2 =6*iniCorrRef.get_sigma_p2()*funcArgsRef.L_j*funcArgsRef.L_j*funcArgsRef.L_k*M_kj()*iniCorrRef.a_02(q_)*funcArgsRef.mu_j;
		if(Direction_ == true)
		{
			if((fabs(L2) < 1.0e-3))
			{
				return q_ * q_ * iniCorrRef.a_11(q_) * pow(B_ij(q_), double(n_)) * 
						(1.0-L2)  * exp(-Q) * funcArgsRef.mu_j * funcArgsRef.L_j * Q_kj(q_);
			}
			
			else if( (exp(-L2) > 0.0) && (fabs(Q) < 1.0e-4))
			{
				return  q_ * q_ * iniCorrRef.a_11(q_) * pow(B_ij(q_), double(n_)) * 
						exp(-L2) * funcArgsRef.mu_j * funcArgsRef.L_j * (1-Q) * Q_kj(q_);
			}
		
			else
			{
				return q_ * q_ * iniCorrRef.a_11(q_) * pow(B_ij(q_), double(n_)) * 
						funcArgsRef.mu_j * funcArgsRef.L_j * exp(-L2-Q) * Q_kj(q_);
			}
	   }
	   
	   if(Direction_ == false)
	   {	
			if((fabs(L2) < 1.0e-3))
			{
				return q_ * q_ * iniCorrRef.a_11(q_) * pow(B_ij(q_), double(n_)) 
						*(1-L2)* exp(-Q) * (A1);
			}
		
			else if( (exp(-L2) > 0.0) && (fabs(Q) < 1.0e-4))
			{
				return q_ * q_ * iniCorrRef.a_11(q_) * pow(B_ij(q_), double(n_)) 
						*(1-Q)* exp(-L2) * (A1);
			}
	
			else
			{
				return q_ * q_ * iniCorrRef.a_11(q_) * pow(B_ij(q_), double(n_)) 
						* exp(-L2-Q) * (A1);
			}

		}
	return 0.0;
	}

	double FirstOrderDerivedDensityMomentumCorrLine::Kernel_n_plus_2(const double K_ij_, const double q_, const std::size_t n_, const bool Direction_) const
	{
		double L2 = fabs(funcArgsRef.mu_ij * L2_ij());
		double Q = Q_ij(q_);
		double A4 = 6* funcArgsRef.L_k*funcArgsRef.L_j*M_kj()*iniCorrRef.a_02(q_)*iniCorrRef.get_sigma_p2();
		double A1 = (Q_kj(q_)*funcArgsRef.L_j*funcArgsRef.mu_j-funcArgsRef.mu_k*funcArgsRef.L_k);
		double A2 =6*iniCorrRef.get_sigma_p2()*funcArgsRef.L_j*funcArgsRef.L_j*funcArgsRef.L_k*M_kj()*iniCorrRef.a_02(q_)*funcArgsRef.mu_j;

		if(Direction_==true)
		{
			if((fabs(L2) < 1.0e-3))
			{
				return q_ * iniCorrRef.a_11(q_) * pow(B_ij(q_), double(n_)) * (1-L2) * exp(-Q)
						* 6.0 * funcArgsRef.mu_j * funcArgsRef.L_j  * funcArgsRef.L_j*iniCorrRef.get_sigma_p2()
						* funcArgsRef.L_k * M_kj() * iniCorrRef.a_02(q_) * (n_+1)/funcArgsRef.K_ij;
			}
			
			else if( (exp(-L2) > 0.0) && (fabs(Q) < 1.0e-4))
			{
				return  q_ * iniCorrRef.a_11(q_) * pow(B_ij(q_), double(n_)) * exp(-L2) * (1-Q) 
						 * 6.0 * funcArgsRef.mu_j * funcArgsRef.L_j  * funcArgsRef.L_j*iniCorrRef.get_sigma_p2()
						* funcArgsRef.L_k * M_kj() * iniCorrRef.a_02(q_) * (n_+1)/funcArgsRef.K_ij;
			}
		
			else
			{
				return q_ * iniCorrRef.a_11(q_) * pow(B_ij(q_), double(n_)) * exp(-L2-Q)
						 * 6.0 * funcArgsRef.mu_j * funcArgsRef.L_j  *funcArgsRef.L_j*iniCorrRef.get_sigma_p2()
						 * funcArgsRef.L_k * M_kj() * iniCorrRef.a_02(q_)*(n_+1)/funcArgsRef.K_ij;
			
			}
		}
		if(Direction_==false)
		{
			if((fabs(L2) < 1.0e-3))
			{
				return q_ * iniCorrRef.a_11(q_) * pow(B_ij(q_), double(n_)) * (1-L2)*exp(-Q)
						*(n_+1)*A2/funcArgsRef.K_ij;
			}
			
			else if( (exp(-L2) > 0.0) && (fabs(Q) < 1.0e-4))
			{
				return q_ * iniCorrRef.a_11(q_) * pow(B_ij(q_), double(n_)) * exp(-L2)*(1-Q)
						*(n_+1)*A2/funcArgsRef.K_ij;
			}
		
			else
			{
				return q_ * iniCorrRef.a_11(q_) * pow(B_ij(q_), double(n_)) * exp(-L2-Q)
						*(n_+1)*A2/funcArgsRef.K_ij;
			
			}
		}
		return 0.0;
	}

	// ***********************************************************
	// ***                   public members                    ***
	// ***********************************************************


    FirstOrderDerivedDensityMomentumCorrLine::FirstOrderDerivedDensityMomentumCorrLine(const iniCorr & iniCorrRef_, const funcArgs & funcArgsRef_,
                                                              const std::size_t BinsLevin_, const double q_min_, const double q_max_):
	CorrLineThimo(iniCorrRef_, funcArgsRef_)
	{
		BinsLevin = BinsLevin_;
		q_min = q_min_;
		q_max = q_max_;
	}

    
	
	double FirstOrderDerivedDensityMomentumCorrLine::LinearTerm(const bool Direction_) const
	{
		
		double damping = exp(-fabs(funcArgsRef.mu_ij * L2_ij()));
		if(Direction_==true)
		{
			return 0.0;
		}
		if(Direction_==false)
		{
			return -3*damping * funcArgsRef.mu_k * funcArgsRef.L_k/ funcArgsRef.K_ij * 
					iniCorrRef.get_sigma_p2() * iniCorrRef.RenormalisedInitialPowerSpectrum(funcArgsRef.K_ij);

		}
		return 0.0;
	}


	double FirstOrderDerivedDensityMomentumCorrLine::nth_NonLinearTerm(const std::size_t n_, const bool Direction_) const
	{
		
		if(n_ == 0)
		{
			//throw std::runtime_error("In DensityMomentumCorrLine::nth_NonLinearTerm: Bessel order n must be larger than zero!");
			return 0.0;
		}
		else if(n_ == 1)
		{
			Integral_OSC_1D::Integrand kernel_n1 = std::bind(&FirstOrderDerivedDensityMomentumCorrLine::Kernel_n1, this, std::placeholders::_1, std::placeholders::_2, Direction_);
			Integral_OSC_1D integrate_n1(kernel_n1, LevinInt::OSC_SPH_BESSEL, 1, q_min, q_max, KFT::BinningType::LOG_BINNING, BinsLevin);
			return 36.0 * M_PI *iniCorrRef.get_sigma_p2()*integrate_n1(funcArgsRef.K_ij);
		}
		else if(n_ == 2)
		{
			Integral_OSC_1D::Integrand kernel_n2 = std::bind(&FirstOrderDerivedDensityMomentumCorrLine::Kernel_n2, this, std::placeholders::_1, std::placeholders::_2, Direction_);
			Integral_OSC_1D integrate_n2(kernel_n2, LevinInt::OSC_SPH_BESSEL, 2, q_min, q_max, KFT::BinningType::LOG_BINNING, BinsLevin);
			return 36.0 * M_PI * iniCorrRef.get_sigma_p2()*integrate_n2(funcArgsRef.K_ij);
		}
		else
		{
			Integral_OSC_1D::Integrand Kernel_n_plus_1 = std::bind(&FirstOrderDerivedDensityMomentumCorrLine::Kernel_n, this, std::placeholders::_1, std::placeholders::_2, n_, Direction_);
			Integral_OSC_1D integrate_n_plus_1(Kernel_n_plus_1, LevinInt::OSC_SPH_BESSEL, n_+1, q_min, q_max, KFT::BinningType::LOG_BINNING, BinsLevin);
			Integral_OSC_1D::Integrand Kernel_n_plus_2 = std::bind(&FirstOrderDerivedDensityMomentumCorrLine::Kernel_n, this, std::placeholders::_1, std::placeholders::_2, n_, Direction_);
			Integral_OSC_1D integrate_n_plus_2(Kernel_n_plus_2, LevinInt::OSC_SPH_BESSEL, n_+2, q_min, q_max, KFT::BinningType::LOG_BINNING, BinsLevin);
			return 36.0 * M_PI *iniCorrRef.get_sigma_p2()*(integrate_n_plus_1(funcArgsRef.K_ij) + integrate_n_plus_2(funcArgsRef.K_ij));
		}
		return 0.0;
	}

	double FirstOrderDerivedDensityMomentumCorrLine::NonLinearTerms_FixedOrderSum(const std::size_t MaxBesselOrder_, const bool Direction_) const
	{
		double Result = 0.0;
		for(std::size_t i=0; i < MaxBesselOrder_+1; i++)
		{
			Result += nth_NonLinearTerm(i, Direction_);
		}
		return Result;
	}

	double FirstOrderDerivedDensityMomentumCorrLine::operator()(const double ConvergenceThreshold_, const bool Direction_) const
	{
		if(ConvergenceThreshold_ < 0.0)
		{
			throw std::runtime_error("In FirstOrderDerivedDensityMomentumCorrLine::operator(): convergence threshold must be positive!");
		}
		/*
		double Result = nth_NonLinearTerm(0, Direction_);
		bool converged = false;

		for(std::size_t n=1; converged == false; n++)
		{
			double Result_n = 0.0;
			Result_n = nth_NonLinearTerm(n,Direction_);
			if( (fabs(Result) == 0.0) || fabs(Result_n/Result) < ConvergenceThreshold_ )
			{
				converged = true;
			}
			Result += Result_n;
		}
		*/
		;
		double Result=nth_NonLinearTerm(1, Direction_)+nth_NonLinearTerm(2, Direction_)+nth_NonLinearTerm(3, Direction_)
		+nth_NonLinearTerm(4, Direction_) +nth_NonLinearTerm(5, Direction_)
		+nth_NonLinearTerm(0, Direction_);
		return LinearTerm(Direction_) + Result;
	}
}