#include "FreeCumulant_TwoPoint.hpp"
#include "../CorrLine/allCorrLines.hpp"


////////////////////////////////////////////////////////////////////////////////////////////////////
// public
////////////////////////////////////////////////////////////////////////////////////////////////////

KFT::FreeCumulant_TwoPoint::FreeCumulant_TwoPoint(const iniCorr& iniCorr_Ref_, const ParticleDynamics& ParticleDynamics_Ref_) :
FreeCumulant(2, iniCorr_Ref_, ParticleDynamics_Ref_)
{}

double KFT::FreeCumulant_TwoPoint::IntegratedInteractingCumulant(const double k1, const double t1) const
{		
	
		astro::integrator TimeIntegral([this, &t1, &k1](double t) {return this->TimeKernelInteractingCumulant(t, k1, t1);});
		return  TimeIntegral(0,t1);
}
double KFT::FreeCumulant_TwoPoint::IntegratedInteractingCumulantLinear(const double k1, const double t1) const
{		
	
		astro::integrator TimeIntegral([this, &t1, &k1](double t) {return this->TimeKernelInteractingCumulantLinear(t, k1, t1);});
		return  TimeIntegral(0,t1);
}
double KFT::FreeCumulant_TwoPoint::IntegratedInteractingCumulantNonLinear(const double k1, const double t1) const
{		
	
		astro::integrator TimeIntegral([this, &t1, &k1](double t) {return this->TimeKernelInteractingCumulantNonLinear(t, k1, t1);});
		return  TimeIntegral(0,t1);
}

double KFT::FreeCumulant_TwoPoint::TimeKernelInteractingCumulant(  double t ,const double k1,const double t1) const
{	
	const double gqp1 = g_qp(t1);
	const double gqp = g_qp(t);

	CorrLine::funcArgs CorrLineArgs;

	CorrLineArgs.K_ij = k1;
	CorrLineArgs.L_i = k1*gqp1;
	CorrLineArgs.L_j = k1*gqp;
	CorrLineArgs.mu_i = 1.0;
	CorrLineArgs.mu_j = -1.0;
	CorrLineArgs.mu_ij = -1.0;

	const double DampingAmplitude = 0.5*iniCorr_Ref.get_sigma_p2() * (std::pow(CorrLineArgs.L_i, 2.0) + std::pow(CorrLineArgs.L_j, 2.0)
																		- 2.0*std::fabs(CorrLineArgs.mu_ij * CorrLineArgs.L_i * CorrLineArgs.L_j));
	double damping = iniCorr_Ref.get_sigma_p2()*k1*k1*(t*t-t*t1);
	const std::size_t BinsLevin = 256;
	const double q_IntegrationMin = iniCorr_Ref.get_MinimumValueOf_q();
	const double q_IntegrationMax = 1.0e10;
	const double ConvergenceThreshold = 1.0e-3;

	DensityCorrLine DCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	DensityMomentumCorrLine DMCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	MomentumCorrLine MCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	DensityMomentumLoopCorrLine DMLCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);

	return 2*std::pow(rhoBar(),2.0) * std::exp(-DampingAmplitude) *std::exp(-damping)
		   * (t1-t)*1.5/(1+t)/(1+t)*
		   (DCL(ConvergenceThreshold) +DMCL(ConvergenceThreshold,true) +DMCL(ConvergenceThreshold,false) +
			MCL(ConvergenceThreshold) +DMLCL(ConvergenceThreshold));
}

double KFT::FreeCumulant_TwoPoint::TimeKernelInteractingCumulantNonLinear(  double t ,const double k1,const double t1) const
{	
	const double gqp1 = g_qp(t1);
	const double gqp = g_qp(t);

	CorrLine::funcArgs CorrLineArgs;

	CorrLineArgs.K_ij = k1;
	CorrLineArgs.L_i = k1*gqp1;
	CorrLineArgs.L_j = k1*gqp;
	CorrLineArgs.mu_i = 1.0;
	CorrLineArgs.mu_j = -1.0;
	CorrLineArgs.mu_ij = -1.0;

	const double DampingAmplitude = 0.5*iniCorr_Ref.get_sigma_p2() * (std::pow(CorrLineArgs.L_i, 2.0) + std::pow(CorrLineArgs.L_j, 2.0)
																		- 2.0*std::fabs(CorrLineArgs.mu_ij * CorrLineArgs.L_i * CorrLineArgs.L_j));
	double damping = iniCorr_Ref.get_sigma_p2()*k1*k1*(t*t-t*t1);
	const std::size_t BinsLevin = 256;
	const double q_IntegrationMin = iniCorr_Ref.get_MinimumValueOf_q();
	const double q_IntegrationMax = 1.0e10;
	const double ConvergenceThreshold = 1.0e-3;

	DensityCorrLine DCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	DensityMomentumCorrLine DMCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	MomentumCorrLine MCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	DensityMomentumLoopCorrLine DMLCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);

	return 2*std::pow(rhoBar(),2.0) * std::exp(-DampingAmplitude) *std::exp(-damping)
		   * (t1-t)*1.5/(1+t)/(1+t)*
		   (DCL.NonLinearTerms_FixedOrderSum(10,true) +DMCL.NonLinearTerms_FixedOrderSum(10,true) +DMCL.NonLinearTerms_FixedOrderSum(10,false) +
			MCL.NonLinearTerms_FixedOrderSum(10,true) +DMLCL.NonLinearTerms_FixedOrderSum(10,true));

}
double KFT::FreeCumulant_TwoPoint::TimeKernelInteractingCumulantLinear(  double t ,const double k1,const double t1) const
{	
	const double gqp1 = g_qp(t1);
	const double gqp = g_qp(t);

	CorrLine::funcArgs CorrLineArgs;

	CorrLineArgs.K_ij = k1;
	CorrLineArgs.L_i = k1*gqp1;
	CorrLineArgs.L_j = k1*gqp;
	CorrLineArgs.mu_i = 1.0;
	CorrLineArgs.mu_j = -1.0;
	CorrLineArgs.mu_ij = -1.0;

	const double DampingAmplitude = 0.5*iniCorr_Ref.get_sigma_p2() * (std::pow(CorrLineArgs.L_i, 2.0) + std::pow(CorrLineArgs.L_j, 2.0)
																		- 2.0*std::fabs(CorrLineArgs.mu_ij * CorrLineArgs.L_i * CorrLineArgs.L_j));
	//const double damping = 0.5*iniCorr_Ref.get_sigma_p2()*k1*k1*(gqp1-gqp)*(gqp1-gqp);
	double damping = iniCorr_Ref.get_sigma_p2()*k1*k1*(t*t-t*t1);
	const std::size_t BinsLevin = 256;
	const double q_IntegrationMin = iniCorr_Ref.get_MinimumValueOf_q();
	const double q_IntegrationMax = 1.0e10;
	const double ConvergenceThreshold = 1.0e-3;

	DensityCorrLine DCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	DensityMomentumCorrLine DMCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	MomentumCorrLine MCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	DensityMomentumLoopCorrLine DMLCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);

	return 2*std::pow(rhoBar(),2.0) * std::exp(-DampingAmplitude) *std::exp(-damping)
		   * (t1-t)*1.5/(1+t)/(1+t)*
		   (DCL.LinearTerm() +DMCL.LinearTerm(true) +DMCL.LinearTerm(false) +
			MCL.LinearTerm() +DMLCL.LinearTerm());
}


double KFT::FreeCumulant_TwoPoint::fF(const PhaseSpaceArguments& Args_) const
{
	double k1 = Args_[0];
	double l1 = Args_[1];
	double mu_kl1 = Args_[3];
	double t1 = Args_[6];
	double t2 = Args_[7];

	return -pot(k1, t2) * rhoBar() * (k1*k1*g_qp(t1,t2) + k1 * l1 * mu_kl1 * g_pp(t1, t2)) *
	       exp(- 0.5*sigma_v2() * ( pow(k1*(g_qp(t1) - g_qp(t2)), 2.0) + pow(l1* g_pp(t1),2.0) + 2.0*k1*l1*mu_kl1 * (g_qp(t1) - g_qp(t2)) * g_pp(t1) ));
}

double KFT::FreeCumulant_TwoPoint::rF(const ConfigurationSpaceArguments& Args_) const
{
	return fF(PhaseSpaceArguments({Args_[0], 0.0, 0.0, 0.0, 0.0, 0.0, Args_[1], Args_[2]}));
}

double KFT::FreeCumulant_TwoPoint::ff_1p(const PhaseSpaceArguments& Args_) const
{
	double k1 = Args_[0];
	double l1 = Args_[1];
	double l2 = Args_[2];
	double mu_kl1 = Args_[3];
	double mu_kl2 = Args_[4];
	double mu_l1l2 = Args_[5];
	double t1 = Args_[6];
	double t2 = Args_[7];
	
	double g_qp_Difference = g_qp(t1) - g_qp(t2);
	double gpp1 = g_pp(t1);
	double gpp2 = g_pp(t2);

	double DampingAmplitude = 0.5 * iniCorr_Ref.get_sigma_p2() * ( (k1*k1 * g_qp_Difference * g_qp_Difference) + (l1*l1*gpp1*gpp1 + l2*l2*gpp2*gpp2) +
																	(2.0 * k1 * l1 * mu_kl1 * g_qp_Difference * gpp1) + (2.0 * k1*l2*mu_kl2 * g_qp_Difference * gpp2) +
																	(2.0 * l1 * l2 * mu_l1l2 * gpp1 * gpp2) );

	return rhoBar() * exp(-DampingAmplitude);
}

double KFT::FreeCumulant_TwoPoint::rr_1p(const ConfigurationSpaceArguments& Args_) const
{
	return ff_1p(PhaseSpaceArguments({Args_[0], 0.0, 0.0, 0.0, 0.0, 0.0, Args_[1], Args_[2]}));
}

double KFT::FreeCumulant_TwoPoint::ff_2p(const PhaseSpaceArguments& Args_) const
{
	const double k1 = Args_[0];
	const double l1 = Args_[1];
	const double l2 = Args_[2];
	const double mu_kl1 = Args_[3];
	const double mu_kl2 = Args_[4];
	const double mu_l1l2 = Args_[5];
	const double t1 = Args_[6];
	const double t2 = Args_[7];

	const double gqp1 = g_qp(t1);
	const double gqp2 = g_qp(t2);
	const double gpp1 = g_pp(t1);
	const double gpp2 = g_pp(t2);

	CorrLine::funcArgs CorrLineArgs;

	CorrLineArgs.K_ij = k1;
	CorrLineArgs.L_i = std::sqrt( (k1*k1*gqp1*gqp1) + (l1*l1*gqp1*gqp1) + (2.0*k1*l1*mu_kl1*gqp1*gpp1) );
	CorrLineArgs.L_j = std::sqrt( (k1*k1*gqp1*gqp1) + (l2*l2*gqp2*gqp2) + (2.0*k1*l2*mu_kl2*gqp1*gpp2) );
	CorrLineArgs.mu_i = (k1*gqp1 + l1*mu_kl1*gpp1) / CorrLineArgs.L_i;
	CorrLineArgs.mu_j = (-k1*gqp2 + l2*mu_kl2*gpp2) / CorrLineArgs.L_j;
	CorrLineArgs.mu_ij = ( (-k1*k1*gqp1*gqp2) + k1*(l2*mu_kl2*gqp1*gpp2 - l1*mu_kl1*gqp2*gpp1) + (l1*l2*mu_l1l2*gpp1*gpp2) ) / CorrLineArgs.L_i / CorrLineArgs.L_j;

	const double DampingAmplitude = 0.5*iniCorr_Ref.get_sigma_p2() * (std::pow(CorrLineArgs.L_i, 2.0) + std::pow(CorrLineArgs.L_j, 2.0)
																		- 2.0*std::fabs(CorrLineArgs.mu_ij * CorrLineArgs.L_i * CorrLineArgs.L_j));

	const std::size_t BinsLevin = 256;
	const double q_IntegrationMin = iniCorr_Ref.get_MinimumValueOf_q();
	const double q_IntegrationMax = 1.0e10;
	const double ConvergenceThreshold = 1.0e-3;
	
	DensityCorrLine DCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	DensityMomentumCorrLine DMCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	MomentumCorrLine MCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	DensityMomentumLoopCorrLine DMLCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	
	return std::pow(rhoBar(),2.0) * std::exp(-DampingAmplitude) * (DCL(ConvergenceThreshold) +
																   DMCL(ConvergenceThreshold,true) +
																   DMCL(ConvergenceThreshold,false) +
																   MCL(ConvergenceThreshold) +
																   DMLCL(ConvergenceThreshold));
}

double KFT::FreeCumulant_TwoPoint::Linearff_2p(const PhaseSpaceArguments& Args_) const
{
	const double k1 = Args_[0];
	const double l1 = Args_[1];
	const double l2 = Args_[2];
	const double mu_kl1 = Args_[3];
	const double mu_kl2 = Args_[4];
	const double mu_l1l2 = Args_[5];
	const double t1 = Args_[6];
	const double t2 = Args_[7];

	const double gqp1 = g_qp(t1);
	const double gqp2 = g_qp(t2);
	const double gpp1 = g_pp(t1);
	const double gpp2 = g_pp(t2);

	CorrLine::funcArgs CorrLineArgs;

	CorrLineArgs.K_ij = k1;
	CorrLineArgs.L_i = std::sqrt( (k1*k1*gqp1*gqp1) + (l1*l1*gqp1*gqp1) + (2.0*k1*l1*mu_kl1*gqp1*gpp1) );
	CorrLineArgs.L_j = std::sqrt( (k1*k1*gqp1*gqp1) + (l2*l2*gqp2*gqp2) + (2.0*k1*l2*mu_kl2*gqp1*gpp2) );
	CorrLineArgs.mu_i = (k1*gqp1 + l1*mu_kl1*gpp1) / CorrLineArgs.L_i;
	CorrLineArgs.mu_j = (-k1*gqp2 + l2*mu_kl2*gpp2) / CorrLineArgs.L_j;
	CorrLineArgs.mu_ij = ( (-k1*k1*gqp1*gqp2) + k1*(l2*mu_kl2*gqp1*gpp2 - l1*mu_kl1*gqp2*gpp1) + (l1*l2*mu_l1l2*gpp1*gpp2) ) / CorrLineArgs.L_i / CorrLineArgs.L_j;

	const double DampingAmplitude = 0.5*iniCorr_Ref.get_sigma_p2() * (std::pow(CorrLineArgs.L_i, 2.0) + std::pow(CorrLineArgs.L_j, 2.0)
																		- 2.0*std::fabs(CorrLineArgs.mu_ij * CorrLineArgs.L_i * CorrLineArgs.L_j));

	const std::size_t BinsLevin = 256;
	const double q_IntegrationMin = iniCorr_Ref.get_MinimumValueOf_q();
	const double q_IntegrationMax = 1.0e10;
	const double ConvergenceThreshold = 1.0e-3;
	
	DensityCorrLine DCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	DensityMomentumCorrLine DMCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	MomentumCorrLine MCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	DensityMomentumLoopCorrLine DMLCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	
	return std::pow(rhoBar(),2.0) * std::exp(-DampingAmplitude) * (DCL.LinearTerm()+
																   DMCL.LinearTerm(true)+
																   DMCL.LinearTerm(false)+
																   MCL.LinearTerm() +
																   DMLCL.LinearTerm());

}
double KFT::FreeCumulant_TwoPoint::NonLinearff_2p(const PhaseSpaceArguments& Args_) const
{
	const double k1 = Args_[0];
	const double l1 = Args_[1];
	const double l2 = Args_[2];
	const double mu_kl1 = Args_[3];
	const double mu_kl2 = Args_[4];
	const double mu_l1l2 = Args_[5];
	const double t1 = Args_[6];
	const double t2 = Args_[7];

	const double gqp1 = g_qp(t1);
	const double gqp2 = g_qp(t2);
	const double gpp1 = g_pp(t1);
	const double gpp2 = g_pp(t2);

	CorrLine::funcArgs CorrLineArgs;

	CorrLineArgs.K_ij = k1;
	CorrLineArgs.L_i = std::sqrt( (k1*k1*gqp1*gqp1) + (l1*l1*gqp1*gqp1) + (2.0*k1*l1*mu_kl1*gqp1*gpp1) );
	CorrLineArgs.L_j = std::sqrt( (k1*k1*gqp1*gqp1) + (l2*l2*gqp2*gqp2) + (2.0*k1*l2*mu_kl2*gqp1*gpp2) );
	CorrLineArgs.mu_i = (k1*gqp1 + l1*mu_kl1*gpp1) / CorrLineArgs.L_i;
	CorrLineArgs.mu_j = (-k1*gqp2 + l2*mu_kl2*gpp2) / CorrLineArgs.L_j;
	CorrLineArgs.mu_ij = ( (-k1*k1*gqp1*gqp2) + k1*(l2*mu_kl2*gqp1*gpp2 - l1*mu_kl1*gqp2*gpp1) + (l1*l2*mu_l1l2*gpp1*gpp2) ) / CorrLineArgs.L_i / CorrLineArgs.L_j;

	const double DampingAmplitude = 0.5*iniCorr_Ref.get_sigma_p2() * (std::pow(CorrLineArgs.L_i, 2.0) + std::pow(CorrLineArgs.L_j, 2.0)
																		- 2.0*std::fabs(CorrLineArgs.mu_ij * CorrLineArgs.L_i * CorrLineArgs.L_j));

	const std::size_t BinsLevin = 256;
	const double q_IntegrationMin = iniCorr_Ref.get_MinimumValueOf_q();
	const double q_IntegrationMax = 1.0e10;
	const double ConvergenceThreshold = 1.0e-3;
	
	DensityCorrLine DCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	DensityMomentumCorrLine DMCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	MomentumCorrLine MCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	DensityMomentumLoopCorrLine DMLCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	
	return std::pow(rhoBar(),2.0) * std::exp(-DampingAmplitude) * (DCL.NonLinearTerms_FixedOrderSum(6,true)+
																   DMCL.NonLinearTerms_FixedOrderSum(6,true)+
																   DMCL.NonLinearTerms_FixedOrderSum(6,false)+
																   MCL.NonLinearTerms_FixedOrderSum(6,true)+
																   DMLCL.NonLinearTerms_FixedOrderSum(6,true));

}

double KFT::FreeCumulant_TwoPoint::rr_2p(const ConfigurationSpaceArguments& Args_) const
{
	return ff_2p(PhaseSpaceArguments({Args_[0], 0.0, 0.0, 0.0, 0.0, 0.0, Args_[1], Args_[2]}));
}

double KFT::FreeCumulant_TwoPoint::TimekernelfreePS(double t2, const double k1, const double t1) const
{	
	const double l1 = 0;
	const double l2 = 0;
	const double mu_kl1 = 0;
	const double mu_kl2 = 0;
	const double mu_l1l2 = 0;


	const double gqp1 = g_qp(t1);
	const double gqp2 = g_qp(t2);
	const double gpp1 = g_pp(t1);
	const double gpp2 = g_pp(t2);

	CorrLine::funcArgs CorrLineArgs;

	CorrLineArgs.K_ij = k1;
	CorrLineArgs.L_i = std::sqrt( (k1*k1*gqp1*gqp1) + (l1*l1*gqp1*gqp1) + (2.0*k1*l1*mu_kl1*gqp1*gpp1) );
	CorrLineArgs.L_j = std::sqrt( (k1*k1*gqp1*gqp1) + (l2*l2*gqp2*gqp2) + (2.0*k1*l2*mu_kl2*gqp1*gpp2) );
	CorrLineArgs.mu_i = (k1*gqp1 + l1*mu_kl1*gpp1) / CorrLineArgs.L_i;
	CorrLineArgs.mu_j = (-k1*gqp2 + l2*mu_kl2*gpp2) / CorrLineArgs.L_j;
	CorrLineArgs.mu_ij = ( (-k1*k1*gqp1*gqp2) + k1*(l2*mu_kl2*gqp1*gpp2 - l1*mu_kl1*gqp2*gpp1) + (l1*l2*mu_l1l2*gpp1*gpp2) ) / CorrLineArgs.L_i / CorrLineArgs.L_j;

	const double DampingAmplitude = 0.5*iniCorr_Ref.get_sigma_p2() * (std::pow(CorrLineArgs.L_i, 2.0) + std::pow(CorrLineArgs.L_j, 2.0)
																		- 2.0*std::fabs(CorrLineArgs.mu_ij * CorrLineArgs.L_i * CorrLineArgs.L_j));

	const std::size_t BinsLevin = 256;
	const double q_IntegrationMin = iniCorr_Ref.get_MinimumValueOf_q();
	const double q_IntegrationMax = 1.0e10;
	const double ConvergenceThreshold = 1.0e-3;
	
	DensityCorrLine DCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	DensityMomentumCorrLine DMCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	MomentumCorrLine MCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	DensityMomentumLoopCorrLine DMLCL(iniCorr_Ref, CorrLineArgs, BinsLevin, q_IntegrationMin, q_IntegrationMax);
	
	return std::pow(rhoBar(),2.0) * std::exp(-DampingAmplitude) * (DCL(ConvergenceThreshold) +
																   DMCL(ConvergenceThreshold,true) +
																   DMCL(ConvergenceThreshold,false) +
																   MCL(ConvergenceThreshold) +
																   DMLCL(ConvergenceThreshold));


}
double KFT::FreeCumulant_TwoPoint::rr_2pintegrated(const double k1, const double t1) const
{	
	astro::integrator TimeIntegral3([this, &t1, &k1](double t) {return this->TimekernelfreePS(t,k1, t1 );});

    
    return TimeIntegral3(0,t1);
}


double KFT::FreeCumulant_TwoPoint::TimeKernelFelixLinear(double t2, const double k1, const double t1) const
{
	

	return 2*fF(PhaseSpaceArguments({k1, 0.0, 0.0, 0.0, 0.0, 0.0, t1, t2 }))*
			Linearff_2p(PhaseSpaceArguments({k1, 0.0, 0.0, 0.0, 0.0, 0.0, t1 , t2}));
}
double KFT::FreeCumulant_TwoPoint::TimeKernelFelixNonLinear(double t2, const double k1, const double t1) const
{
	

	return 2*fF(PhaseSpaceArguments({k1, 0.0, 0.0, 0.0, 0.0, 0.0, t1, t2 }))*
			NonLinearff_2p(PhaseSpaceArguments({k1, 0.0, 0.0, 0.0, 0.0, 0.0, t1 , t2}));
}
double KFT::FreeCumulant_TwoPoint::TimeKernelFelix(double t2, const double k1, const double t1) const
{
	

	return 2*fF(PhaseSpaceArguments({k1, 0.0, 0.0, 0.0, 0.0, 0.0, t1, t2 }))*
			ff_2p(PhaseSpaceArguments({k1, 0.0, 0.0, 0.0, 0.0, 0.0, t1 , t2}));
}
double KFT::FreeCumulant_TwoPoint::InteractingCumulantFelix(const double k1, const double t1) const
{	


	astro::integrator TimeIntegral1([this, &k1, &t1](double t) {return this->TimeKernelFelix(t, k1, t1);});

    
    return TimeIntegral1(0,t1);
}

double KFT::FreeCumulant_TwoPoint::InteractingCumulantFelixLinear(const double k1, const double t1) const
{	
	astro::integrator TimeIntegral3([this, &k1, &t1](double t) {return this->TimeKernelFelixLinear(t, k1, t1);});

    
    return TimeIntegral3(0,t1);
}
double KFT::FreeCumulant_TwoPoint::InteractingCumulantFelixNonLinear(const double k1, const double t1) const
{	


	astro::integrator TimeIntegral3([this, &k1, &t1](double t) {return this->TimeKernelFelixNonLinear(t, k1, t1);});

    
    return TimeIntegral3(0,t1);
}