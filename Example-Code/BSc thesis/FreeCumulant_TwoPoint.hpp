#ifndef _FREE_CUMULANT_TWOPOINT_
#define _FREE_CUMULANT_TWOPOINT_

#include <array>

#include "FreeCumulant.hpp"
#include "astro/utilities/integrator.h"

namespace KFT
{
	class FreeCumulant_TwoPoint : public FreeCumulant
	{
		public:

			/**
			 * Defines a fixed-length array holding the full phase-space arguments in the format
			 * \f[ \mathrm{PhaseSpaceArguments} = \{ k_1, l_1, l_2, \mu_{k l_1}, \mu_{k l_2}, \mu_{l_1 l_2}, t_1, t_2 \} \;, \f]
			 * where \f$ \mu_{k l_{1,2}} \f$ are the cosines of the angles between \f$ \vec{k}_1 \f$ and \f$ \vec{l}_{1,2} \f$
			 * respectively. This also uses the fact that \f$ \vec{k}_2 = - \vec{k}_1 \f$ by virtue of the leading Dirac delta fucntion.
			 */
			using PhaseSpaceArguments = std::array<double,8>;

			/**
			 * Defines a fixed-length array holding the configuration arguments in the format
			 * \f[ \mathrm{PhaseSpaceArguments} = \{ k_1, t_1, t_2 \} \;, \f]
			 * This also uses the fact that \f$ \vec{k}_2 = - \vec{k}_1 \f$ by virtue of the leading Dirac delta fucntion.
			 */
			using ConfigurationSpaceArguments = std::array<double,3>;

			FreeCumulant_TwoPoint(const iniCorr& iniCorr_Ref_, const ParticleDynamics& ParticleDynamics_Ref_);
	
	
			/**
			 * Calculates the interacting cumulant 
			 * \f[ 2i \int_0^{t_1} \mathrm{d}t_{\bar{2}} \, G_{\rho_1 \mathcal{F}_{\hat{2}}}^{(0,1)} \, G_{\rho_{\hat{2}}\rho_2}^{(0,2)}\f]
			*/

			double IntegratedInteractingCumulant(const double k1, const double t1) const;


			double IntegratedInteractingCumulantLinear(const double k1, const double t1) const;

			double IntegratedInteractingCumulantNonLinear(const double k1, const double t1) const;
			//TimeKernels used for the calculation of IntegratedInteractingCumulant Full

			double TimeKernelInteractingCumulant(  double t,const double k1,const double t1) const;
			//TimeKernels used for the calculation of IntegratedInteractingCumulant Linear

			double TimeKernelInteractingCumulantLinear(  double t,const double k1,const double t1) const;


			double TimeKernelInteractingCumulantNonLinear(  double t,const double k1,const double t1) const;

			/**
			 * By virtue of causality this two-point cumulant is identically zero and only listed for completeness.
			 */
			inline double FF(const PhaseSpaceArguments& Args_) const
			{
				return 0.0;
			}

			/**
			 * Calculates the following expression:
			 * \f$ \tilde{G}^{(0,1)}_{f_1 \mathcal{F}_2} = \bar{\rho}\f$
			 * which defines the full cumulant through
			 * \f[ \mathrm{i} G^{0}_{f_1 \mathcal{F}_2} = (2\pi)^3 \delta_{\mathsc{D}}\left( \vec{k}_1 + \vec{k}_2 \right)(2\pi)^3 \delta_{\mathsc{D}}\left( \vec{l}_2 \right) \tilde{G}^{(0,1)}_{f_1 \mathcal{F}_2} \f]
			 * The arguments \f$ l_2, \mu_{k l_2} \f$ of the passed arguments array are ignored because \f$ \vec{l}_2 = 0 \f$.
			 */
			double fF(const PhaseSpaceArguments& Args_) const;

			/**
			 * Calculates the derived cumulant
			 * \f[ \tilde{G}^{(0,1)}_{\rho_1 \mathcal{F}_2} = \tilde{G}^{(0,1)}_{f_1 \mathcal{F}_2} \big|_{\vec{l}_1 = 0} \;.\f]
			 */
			double rF(const ConfigurationSpaceArguments& Args_) const;

			double ff_1p(const PhaseSpaceArguments& Args_) const;
			double rr_1p(const ConfigurationSpaceArguments& Args_) const;
			
			virtual double ff_2p(const PhaseSpaceArguments& Args_) const;
			double Linearff_2p(const PhaseSpaceArguments& Args_) const;
			double TimekernelfreePS(double t2, const double k1, const double t1) const;
			double NonLinearff_2p(const PhaseSpaceArguments& Args_) const;
			double rr_2p(const ConfigurationSpaceArguments& Args_) const;
			double rr_2pintegrated(const double k1, const double t1) const;
			double TimeKernelFelix(double t2, const double k1, const double t1) const;
			double TimeKernelFelixLinear(double t2, const double k1, const double t1) const;
			double TimeKernelFelixNonLinear(double t2, const double k1, const double t1) const;

			double InteractingCumulantFelix(const double k1, const double t1) const;
			double InteractingCumulantFelixLinear(const double k1, const double t1) const;
			double InteractingCumulantFelixNonLinear(const double k1, const double t1) const;
	};
}



#endif