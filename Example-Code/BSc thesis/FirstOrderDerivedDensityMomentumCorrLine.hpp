#ifndef FirstOrderDerived_DENSITY_MOMENTUM_CORRLINE_HPP_
#define FirstOrderDerived_DENSITY_MOMENTUM_CORRLINE_HPP_

#include "CorrLineThimo.hpp"
#include <array>
#include <functional>

namespace KFT
{

	/**
	* Implementation of the CorrLine base class as  the first order derivative of the momentum density correlation line,
	* in the form of \f$ \vec{L}_k \cdot \vec{\nabla}_{L_i} \mathcal{P}_{p_i \delta_j} (\vec{k}_{ij}) \f$.
	* For definitions of quantities appearing in the kernels see documentation of the base class.
	* Bool direction true =  \f$ \vec{L}_k \cdot \vec{\nabla}_{L_i} \mathcal{P}_{p_j \delta_i} (\vec{k}_{ij}) \f$
	* Bool direction false =  \f$ \vec{L}_k \cdot \vec{\nabla}_{L_i} \mathcal{P}_{p_i \delta_j} (\vec{k}_{ij}) \f$
	* 
	* n_bins is the number of bins the LevinIntegrator will split the integration interval
	* \f$[q_\mathrm{min},q_\mathrm{max}]\f$ into.
	* 
	* Thimo Preis, Heidelberg University, 2018
	*/
	class FirstOrderDerivedDensityMomentumCorrLine : public CorrLineThimo
	{
		protected:
			
			std::size_t BinsLevin;
			double q_min, q_max;
			
		public:
			
			/**
			* Zero-order Bessel kernel function for the non_linear_correction term. Defined as:
			* \f[\mathrm{Kernel}_0(K_{ij},q) := q^2 \exp\left(- \sigma_p^2 \left|\vec{L}_i \cdot \vec{L}_j\right| \right) \left( \exp\left(-Q_{ij}(q)\right) - 1 \right)\times (-Q_{kj}) \;.\f]
			* Note that the argument \c K_ij_ is has no effect here, since all wavevectors are supplied by the \c funcArgs struct
			* of the base class.
			*/
			double Kernel_n1(const double K_ij_, const double q_, const bool Direction_ = true) const;
			
			/**
			* First-order Bessel kernel function for the non_linear_correction term. Defined as:
			* \f[\mathrm{Kernel}_1(K_{ij},q) := q^2 \exp\left(- \sigma_p^2 \left|\vec{L}_i \cdot \vec{L}_j\right| \right) \left( \exp\left(-Q_{ij}(q)\right) - 1 \right) B_{ij}(q) \left[ \frac{L_k M_{kj}}{L_i M_{ij}} - Q_{kj}\right]\;.\f]
			* Note that the argument \c K_ij is has no effect here, since all wavevectors are supplied by the \c funcArgs struct
			* of the base class.
			*/
			double Kernel_n2(const double K_ij_, const double q_, const bool Direction_ = true) const;

			/**
			* nth-order \f$(n>1)\f$ Bessel kernel function for the non_linear_correction term. Defined as:
			* \f[\mathrm{Kernel}_n(K_{ij},q) := q^2 \exp\left(- \sigma_p^2 \left|\vec{L}_i \cdot \vec{L}_j\right| \right) \exp\left(-Q_{ij}(q)\right) \left(B_{ij}(q)\right)^n \left[ \frac{L_k M_{kj}}{L_i M_{ij}} \cdot n - Q_{kj}\right]\\;.\f]
			* Note that the argument \c K_ij is has no effect here, since all wavevectors are supplied by the \c funcArgs struct
			* of the base class.
			*/
			double Kernel_n_plus_1(const double K_ij_, const double q_, const std::size_t n_, const bool Direction_ = true) const;
			

			/**
			* nth-order \f$(n>1)\f$ Bessel kernel function for the non_linear_correction term. Defined as:
			* \f[\mathrm{Kernel}_n(K_{ij},q) := q^2 \exp\left(- \sigma_p^2 \left|\vec{L}_i \cdot \vec{L}_j\right| \right) \exp\left(-Q_{ij}(q)\right) \left(B_{ij}(q)\right)^n \left[ \frac{L_k M_{kj}}{L_i M_{ij}} \cdot n - Q_{kj}\right]\\;.\f]
			* Note that the argument \c K_ij is has no effect here, since all wavevectors are supplied by the \c funcArgs struct
			* of the base class.
			*/
			double Kernel_n_plus_2(const double K_ij_, const double q_, const std::size_t n_, const bool Direction_ = true) const;


			/**
			* Constructor will call the constructor of the base class to set the references.
			*/
			FirstOrderDerivedDensityMomentumCorrLine(const iniCorr & iniCorrRef_, const funcArgs & funcArgsRef_, std::size_t BinsLevin_, double q_min_, double q_max_);
			
			/**
			* Damped contribution from terms linear in the initial powerspectrum. The definition depends on the Direction.
			* Positive Direction (true) is defined by
			* \f[ \mathcal{P}_{\delta_i p_j}^{\mathrm{L}} = 0.0 .\f]
			* Negative Direction (false) is defined by
			* \f[ \mathcal{P}_{p_i \delta_j}^{\mathrm{L}} = -\frac{\mu_k \left|\vec{L}_k\right|}{K_{ij}} \exp\left(- \left|\vec{L}_i \cdot \vec{L}_j\right| \right) \sigma_p^2 \bar{P}_\delta(K_{ij}) \;.\f]
			*/
			double LinearTerm(const bool Direction_ = true) const;

			/**
			 * Returns the radial integral of the n-th order kernel times the Bessel function of the same order.
			 * \f[ \mathcal{P}_{p_i p_j}^{\mathrm{NL},(n)} = 4\pi \int_0^\infty \mathrm{d}q \mathrm{Kernel}_n(q) j_n\left(K_ij q\right) \f]
			 */
			double nth_NonLinearTerm(const std::size_t n_, const bool Direction_ = true) const;
			
			/**
			 * Returns the sum of the non-linear terms up to some specfied order \f$N\f$.
			 * \f[ \mathcal{P}_{p_i p_j}^{\mathrm{NL},N} = \sum_{n=0}^N \mathcal{P}_{p_i p_j}^{\mathrm{NL},(n)}\f]
			 */
			double NonLinearTerms_FixedOrderSum(std::size_t MaxBesselOrder_, bool Direction_ = true) const;

			/**
			 * Returns the full result as
			 * \f[ \mathcal{P}_{p_i p_j} = \mathcal{P}_{p_i p_j}^{\mathrm{L}} + \sum_{n=0}^N \mathcal{P}_{p_i p_j}^{\mathrm{NL},(n)} \;. \f]
			 * The truncation \f$N\f$ is chosen adaptively depending on the specified convergence threshold.
			 */
			double operator()(const double ConvergenceThreshold_ = 1.0e-3, const bool Direction_ = true) const;
	};

}
#endif