#include "Data.h"
#include <iostream>
#include <fstream>
#include <cmath>

#include <assert.h>
#include <CCfits/CCfits>

using namespace std;
using namespace CCfits;

Data Data::instance;

Data::Data()
{
 set_lower_limit();
 set_upper_limit();
}



void Data::load(const char* metadata_file, const char* image_file, const char* sigma_file , const char* psf_file)
{

     std::auto_ptr<FITS> pInfile(new FITS(image_file,Read,true));
     PHDU& image_=pInfile->pHDU(); 

     std::valarray<double>  contents;      
     std::vector<String> name;

     name.push_back("MAGZP");
//     name.push_back("M01_MODMAG");

//     // read all user-specifed, coordinate, and checksum keys in the image
     image_.readAllKeys();
     image_.read(contents);

//     // this doesn't print the data, just header info.
      std::cout << image_ << std::endl;

//     // create an vector of vectors
     ni = image_.axis(0);
     nj = image_.axis(1);

      image_.keyWord(name[0]).value(magzp);
  //   image_.keyWord(name[1]).value(magzero);  
//      
      x_min = 0;
      y_min = 0;
      y_max = nj;
      x_max = ni;
     
      int counter =0; 
      image.assign(ni, vector<double>(nj));
      for (long j = 0;j< nj ; j++)
      {
           for (long i = 0;i < ni ; i++)
           {     
	         image[i][j] = double(contents[counter]);                         // pixels (x,y)
                 counter++;
            }
     }

     max_mag = contents.max();
     min_mag = contents.min();
     std::cout << "Max flux"  << contents.max() << std::endl;
     std::cout << "Min flux"  << contents.min() << std::endl;
     std::cout << "magzp"     << magzp << std::endl;
 
//	fstream fin(metadata_file, ios::in);
//	if(!fin)
//		cerr<<"# ERROR: couldn't open file "<<metadata_file<<"."<<endl;
//	fin>>ni>>nj;
//	fin>>x_min>>x_max>>y_min>>y_max;
//	fin.close();

	// Make sure maximum > minimum
	if(x_max <= x_min || y_max <= y_min)
		cerr<<"# ERROR: strange input in "<<metadata_file<<"."<<endl;

	// Compute pixel widths
	dx = (x_max - x_min)/nj;
	dy = (y_max - y_min)/ni;

	// Check that pixels are square
	if(abs(log(dx/dy)) >= 1E-3)
		cerr<<"# ERROR: pixels aren't square."<<endl;

	//compute_ray_grid();



     std::unique_ptr<FITS> pInfile2(new FITS(sigma_file,Read,true));
     PHDU& sigma_=pInfile2->pHDU(); 

     std::valarray<double>  contents_sigma;      
     std::vector<String> name_sigma;

//     // read all user-specifed, coordinate, and checksum keys in the image
     sigma_.readAllKeys();
     sigma_.read(contents_sigma);

//     // this doesn't print the data, just header info.
      std::cout << sigma_ << std::endl;

//     // create an vector of vectors
     ni = sigma_.axis(0);
     nj = sigma_.axis(1);

      x_min = 0;
      y_min = 0;
      y_max = nj;
      x_max = ni;
     
    
      counter =0; 
      sig.assign(ni, vector<double>(nj));
      for (long j = 0;j< nj ; j++)
      {
           for (long i = 0;i < ni ; i++)
           {     
	         sig[i][j] = double(contents_sigma[counter]);                         // pixels (x,y)
                 counter++;
            }
     }
    std::cout << "Finished imaged starting on psf" << std::endl;

    std::unique_ptr<FITS> pInfile3(new FITS(psf_file,Read,true));
    PHDU& psf_=pInfile3->pHDU(); 
    std::valarray<double> contents_psf;      
    std::vector<String> name_psf;
//     name.push_back("EXPTIME");
//     name.push_back("MAGZP");

     // read all user-specifed, coordinate, and checksum keys in the image
     psf_.readAllKeys();
     psf_.read(contents_psf);

     // this doesn't print the data, just header info.
     std::cout << psf_ << std::endl;

     // create an vector of vectors
     int N = psf_.axis(0);
     int M = psf_.axis(1);

//     std::vector< std::vector<double> > pixels(N, std::vector<double>(M));   // Create a 2D array of size N x M  
//      psf_image.assign(N, std::vector<double>(M,0.00));   // Create a 2D array of size N x M
      psf_image.assign(N, vector<double>(M));
      counter =0; 
      for (int j = 0;j< M; j++)
      {
           for (int i = 0;i < N ; i++)
           {                
                 psf_image[i][j] = double(contents_psf[counter]);                         // pixels (x,y)
                 counter++;
            }
      }

	
} //end-load



// A function to set the limit for the flat prior
void Data::set_lower_limit()
{ 
    double delta = 1e-2;
    lower_limit.assign(15,0);
    lower_limit[0] = 0; //x
    lower_limit[1] = 0; //y
    lower_limit[2] = Data::get_instance().get_magzp() - 18; // I
    lower_limit[3] = 1;  // R
    lower_limit[4] = 0.5; // n 
    lower_limit[5] = 0;  // q
    lower_limit[6] = 0;    // theta
    lower_limit[7] = -1.0;  // boxi

    lower_limit[8] = Data::get_instance().get_magzp() - 18; // I
    lower_limit[9] = delta;  // rout
    lower_limit[10] =  delta; // a                      // alpha [ 0, 10]     (Tayuen, Kim paper 2014) 
    lower_limit[11] =  delta;   // b                    // beta [0 ,2] 
    lower_limit[12] = 0;    // q
    lower_limit[13] = 0 ;  // theta  
    lower_limit[14] = -1.0;  // boxi
}


// A function to set the limit for the flat prior
void Data::set_upper_limit()
{ 
    double delta = 1e-2;
    upper_limit.assign(15,0);
    upper_limit[0] = 100; //x
    upper_limit[1] = 100; //y
    upper_limit[2] = Data::get_instance().get_magzp(); // I
    upper_limit[3] = 50;  // R
    upper_limit[4] = 10; // n 
    upper_limit[5] = 1;  // q
    upper_limit[6] = 180;    // theta
    upper_limit[7] = 1.0;  // boxi

    upper_limit[8] = Data::get_instance().get_magzp(); // I
    upper_limit[9] = 3;  // rout
    upper_limit[10] = 10; // a
    upper_limit[11] = 2.0 - delta;   // b
    upper_limit[12] = 1;    // q
    upper_limit[13] = 180 ;  // theta  
    upper_limit[14] = 1.0;  // boxi
}


//void Data::compute_ray_grid()
//{
//	// Make vectors of the correct size
//	x_rays.assign(ni, vector<double>(nj));
//	y_rays.assign(ni, vector<double>(nj));

//	// Distance between adjacent rays
//	double L = dx;

//	for(size_t i=0; i<x_rays.size(); i++)
//	{
//		for(size_t j=0; j<x_rays[i].size(); j++)
//		{
//			x_rays[i][j] = x_min + (j + 0.5)*L;
//			y_rays[i][j] = y_max - (i + 0.5)*L;
//		}
//	}
//}

