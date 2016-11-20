#ifndef DNest4_GalaxyField_MultiGalaxy
#define DNest4_GalaxyField_MultiGalaxy

#include "DNest4/code/DNest4.h"
#include "MyConditionalPriorMultiGalaxy.h"

#include "DH_Util.h"
#include <valarray>
#include <ostream>
#include <string> 
#include <vector>
#include <sstream>
#include <fstream>                 // includes
#include <iostream>                // includes:  cerr

#define ARMA_DONT_USE_CXX11 
#include <armadillo>


class MultiGalaxy
{
	private:
		DNest4::RJObject<MyConditionalPriorMultiGalaxy> objects;

		// The model image
		std::vector< std::vector<long double> > image;
		std::vector< std::vector<long double> > psf_image;
		void calculate_image();
	
		// How many steps since image was computed from scratch
		int staleness;

                void writefile(const std::vector< std::vector<long double> >  &image);

		// Noise standard deviation
		double sigma;

                double magzp;

               // FFT of the PSF
                arma::cx_mat fft_of_psf;
                bool fft_ready;
                int ni, nj;
                int Ni, Nj;

	public:
		MultiGalaxy();
                 int global_index;
		// Generate the point from the prior
		void from_prior(DNest4::RNG& rng);

		// Metropolis-Hastings proposals
		double perturb(DNest4::RNG& rng);

		// Likelihood function
		double log_likelihood() const;

		// Print to stream
		void print(std::ostream& out) const;

		// Return string with column information
		std::string description() const;

               // Set Limits
                void set_limit(); 

//                // Limits
                std::vector<double> upper_limit;
                std::vector<double> lower_limit;

                bool check_in_range();

	        std::string IntToStr(int n);
           
                // PSF
                void calculate_fft(int Ni, int Nj);
                void blur_image2(std::vector< std::vector<long double> >  &img);

};

#endif

