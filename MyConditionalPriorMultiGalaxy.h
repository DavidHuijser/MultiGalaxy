#ifndef DNest4_GalaxyField_MyConditionalPriorMultiGalaxy
#define DNest4_GalaxyField_MyConditionalPriorMultiGalaxy

#include "DNest4/code/DNest4.h"


#include "DH_Util.h"
#include <string> 
#include <sstream>
#include <fstream>         
#include <iostream>

#include <vector>
#include "MyConditionalPriorMultiGalaxy.h"


// Hyperparameters setting interim prior for galaxy properties
class MyConditionalPriorMultiGalaxy:public DNest4::ConditionalPrior
{
	private:
		// Limits
               std::vector <double> upper_limit;
               std::vector<double> lower_limit;                
               // Hyper-parameters
                double hyper_location;
                double hyper_scale;

//		double x_min, x_max, y_min, y_max;
//		double fluxlim_min, fluxlim_max;
//		double radiuslim_min, radiuslim_max;

		// Lower limit and 1/slope for Pareto interim prior
		// for masses
//		double fluxlim;
//		double gamma;

		// Lower limit and 1/slope for Pareto interim prior
		// for radii
//		double radiuslim;
//		double gamma_radius;

		// Uniform interim prior for radius ratio and mass ratio
          
//		double a1, b1;
//		double a2, b2;


		double perturb_hyperparameters(DNest4::RNG& rng);

	public:
		MyConditionalPriorMultiGalaxy(std::vector<double> lower_limit, std::vector<double> upper_limit);

		void from_prior(DNest4::RNG& rng);

		double log_pdf(const std::vector<double>& vec) const;
		void from_uniform(std::vector<double>& vec) const;
		void to_uniform(std::vector<double>& vec) const;

		void print(std::ostream& out) const;
};

#endif

