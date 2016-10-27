#ifndef DNest4_GalaxyField_Data
#define DNest4_GalaxyField_Data

#include <vector>

#include <assert.h>
#include <CCfits/CCfits>

using namespace CCfits;

class Data
{
	private:
		// Number of pixels
		int ni, nj;

		// Coordinates of image boundaries
		double x_min, x_max, y_min, y_max;

		// Pixel widths
		double dx, dy;

		// Coordinates of pixel centers
		std::vector< std::vector<double> > x_rays;
		std::vector< std::vector<double> > y_rays;

		// The pixels
		std::vector< std::vector<double> > image;

		// Sigma map
		std::vector< std::vector<double> > sig;

		// PSF 
		std::vector< std::vector<double> > psf_image;

                double max_mag;
                double min_mag; 
                double magzp;

                void set_upper_limit();
                void set_lower_limit();
   
                std::vector<double>  lower_limit;
                std::vector<double>   upper_limit;
//		void compute_ray_grid();

	public:
		Data();
		void load(const char* metadata_file, const char* image_file, const char* sigma_file, const char* psf_file);

		// Getters
		int get_ni() const { return ni; }
		int get_nj() const { return nj; }

		std::vector<double> get_upper_limit() const { return upper_limit; }
		std::vector<double> get_lower_limit() const { return lower_limit; }

		double get_x_min() const { return x_min; }
		double get_x_max() const { return x_max; }
		double get_y_min() const { return y_min; }
		double get_y_max() const { return y_max; }
		double get_dx() const { return dx; }
		double get_dy() const { return dy; }
		const std::vector< std::vector<double> >& get_x_rays() const
			{ return x_rays; }
		const std::vector< std::vector<double> >& get_y_rays() const
			{ return y_rays; }
		const std::vector< std::vector<double> >& get_image() const
			{ return image; }
		const std::vector< std::vector<double> >& get_sigma() const
			{ return sig; }
		const std::vector< std::vector<double> >& get_psf() const
			{ return psf_image; }

            	int get_magzp() const { return magzp; }

	// Singleton
	private:
		static Data instance;
	public:
		static Data& get_instance() { return instance; }
};

#endif

