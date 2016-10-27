#include "MultiGalaxy.h"
#include "DNest4/code/DNest4.h"
#include "Data.h"
#include "DH_Util.h"
#include <string> 

#include <vector> 
#include <sstream>
#include <fstream>                 // includes
#include <iostream>


#define ARMA_DONT_USE_CXX11 
#include <armadillo>

/* Mandatory */
 #include <profit/profit.h>
//#include <"/home/dhui890/Documents/Research/Profit/libprofit/profit/profit.h>

/* Depending on which profiles you will use... */
#include <profit/sersic.h>
#include <profit/ferrer.h>
#include <profit/moffat.h>
#include <profit/sky.h>
#include <profit/psf.h>


// use:   more sample_info.txt 

#include <mutex>          // std::mutex

using namespace std;
using namespace DNest4;
using namespace arma;

// create maximum object with 15 parameters
MultiGalaxy::MultiGalaxy() 
:objects(15, 100, false, MyConditionalPriorMultiGalaxy(
      Data::get_instance().get_lower_limit(),Data::get_instance().get_upper_limit() ),  PriorType::log_uniform)
,image(Data::get_instance().get_ni(),vector<long double>(Data::get_instance().get_nj() ))
{
//    set_limit();
    global_index = 0;
    calculate_fft(Data::get_instance().get_ni(), Data::get_instance().get_nj());
    
}


//// A function to set the limit for the flat prior
//void MultiGalaxy::set_limit()
//{ 
//    double delta = 1e-1;
//    lower_limit.assign(15,0);
//    upper_limit.assign(15,0);
//    lower_limit[0] = 0; //x
//    lower_limit[1] = 0; //y
//    lower_limit[2] = Data::get_instance().get_magzp() - 18; // I
//    lower_limit[3] = 1;  // R
//    lower_limit[4] = 0.5; // n 
//    lower_limit[5] = 0;  // q
//    lower_limit[6] = 0;    // theta
//    lower_limit[7] = -1.0;  // boxi

//    lower_limit[8] = Data::get_instance().get_magzp() - 18; // I
//    lower_limit[9] = delta;  // rout
//    lower_limit[10] =  delta; // a                      // alpha [ 0, 10]     (Tayuen, Kim paper 2014) 
//    lower_limit[11] =  delta;   // b                    // beta [0 ,2] 
//    lower_limit[12] = 0;    // q
//    lower_limit[13] = 0 ;  // theta  
//    lower_limit[14] = -1.0;  // boxi

//    upper_limit[0] = 100; //x
//    upper_limit[1] = 100; //y
//    upper_limit[2] = Data::get_instance().get_magzp(); // I
//    upper_limit[3] = 50;  // R
//    upper_limit[4] = 10; // n 
//    upper_limit[5] = 1;  // q
//    upper_limit[6] = 180;    // theta
//    upper_limit[7] = 1.0;  // boxi

//    upper_limit[8] = Data::get_instance().get_magzp(); // I
//    upper_limit[9] = 3;  // rout
//    upper_limit[10] = 10; // a
//    upper_limit[11] = 2.0 - delta;   // b
//    upper_limit[12] = 1;    // q
//    upper_limit[13] = 180 ;  // theta  
//    upper_limit[14] = 1.0;  // boxi

//}

void MultiGalaxy::from_prior(RNG& rng)
{
   objects.from_prior(rng);
   calculate_image();
   sigma = exp(log(1.) + log(1E6)*rng.rand());
   calculate_image();
   //     std::cout << " Finished image " << std::endl; 
} // End-from prior


double MultiGalaxy::perturb(RNG& rng)
{
	double logH = 0.;

	if(rng.rand() <= 0.75)
	{
		logH += objects.perturb(rng);
		calculate_image();
	}
	else
	{
		sigma = log(sigma);
		sigma += log(1E6)*rng.randh();
		sigma = mod(sigma - log(1.), log(1E6)) + log(1.);
		sigma = exp(sigma);
	}

	return logH;
}





double MultiGalaxy::log_likelihood() const
{
  //      std::cout << std::endl << " Started eval loglikelihood " << std::endl;   
	const vector< vector<double> >& data = Data::get_instance().get_image();
	const vector< vector<double> >& sig =  Data::get_instance().get_sigma();

	double logL = 0.;
	double var;
	for(size_t i=0; i<data.size(); i++)
	{
		for(size_t j=0; j<data[i].size(); j++)
		{
			var = sigma*sigma + sig[i][j]*sig[i][j];
			logL += -0.5*log(2.*M_PI*var)
				-0.5*pow(data[i][j] - image[i][j], 2)/var;
		}
	}

	if(std::isnan(logL) || std::isinf(logL))
		logL = -1E300;

        //std::cout << std::endl << " Finished eval loglikelihood sigma:" <<  var << std::endl;    
        //std::cout << data.size() << " " << data[0].size() << std::endl;    
	return logL;

}


string MultiGalaxy::description() const
{
	return string("objects");
        return string("Each column is one of the 20 parameters.");
}



void MultiGalaxy::print(std::ostream& out) const
{	
         out<<setprecision(6);
         for (long i=0; i < long(image.size()); i++)
	  {
	     for (long j=0;  j < long(image[0].size()); j++)
	     {
	       out << image[i][j] <<" "; 
	     } //end-for-j	   	     
	  } //end for -i    
	out<<setprecision(10);
	objects.print(out); out<<' ';
	out<<sigma<<' ';   
} // 

void MultiGalaxy::calculate_image()
{
        global_index++;
        double x, y;
	// Get coordinatesfrom data and create an empty image
        int ni= Data::get_instance().get_ni();
        int nj= Data::get_instance().get_nj();
	image.assign(Data::get_instance().get_ni(), vector<long double>(Data::get_instance().get_nj(), 0.));
	


        bool update = objects.get_removed().size() == 0;
	// Components
	const vector< vector<double> >& components = (update)?(objects.get_added())
                                                 :(objects.get_components());
	if(!update)
	{
		// Zero the image
		image.assign(Data::get_instance().get_ni(),
			vector<long double>(Data::get_instance().get_nj(), 0.));
	}

        for(size_t k=0; k<components.size(); ++k)
	{

           profit::Model *model = new profit::Model();

           model->width = ni; 
           model->height = nj;
           model->magzero = Data::get_instance().get_magzp();
           double half_xbin = model->scale_x/2.;
           double half_ybin = model->scale_y/2.;

           profit::Profile *sersic_profile = model->add_profile("sersic");
           profit::SersicProfile *sp = static_cast<profit::SersicProfile *>(sersic_profile);
           sp->xcen = components[k][0];
           sp->ycen = components[k][1];
           sp->mag = components[k][2];
           sp->re =  components[k][3];
           sp->nser = components[k][4];
           sp->axrat = components[k][5];
           sp->ang = components[k][6];
           sp->box = components[k][7];


           profit::Profile *ferrer_profile = model->add_profile("ferrer");
           profit::FerrerProfile *fp = static_cast<profit::FerrerProfile*>(ferrer_profile);
           fp->xcen = components[k][0];
           fp->ycen = components[k][1];
           fp->mag = components[k][8];
           fp->rout = components[k][9]*components[k][3];
           fp->a = components[k][10];
           fp->b = components[k][11];
           fp->axrat = components[k][12]*components[k][5];  // ax_bar < ax_disc
           fp->ang = components[k][13];
           fp->box = components[k][14];

           // std::cout << "rout= " << params[9]*params[3] << " a= " << params[10] << " b= " << params[11] << std::endl;

            model->evaluate();

            for (int j=0; j < nj; j++) 
            {
	  	x += half_xbin;
	 	y = 0;
	 	for (int k=0; k < ni; k++) 
                {
		       y += half_ybin;                       
                       image[j][k] = image[j][k]  +   model->image[j + k*nj];
         	       y += half_ybin;
		} // end-for
		x += half_xbin;
             } //end-for
             delete model;


       } // end-for
        blur_image2(image);
        if (global_index % 500 == 0)
        {
             writefile(image); 
        } 
} // End-calculate image 


void MultiGalaxy::writefile(const std::vector< std::vector<long double> >  &image)
{

     int sub2 = int(global_index / 500); 
 //   std::stringstream result;
  //  result << sub2;

//     std::string sub = result.str();
//      std::string name="file_" + sub  +".txt";        

      std::string name="file_" +      std::to_string(sub2)  +".txt";        
      std::ofstream myout; // This is the object that sends data to the file. It is used like cout, but now writes in file. 
      myout.open(name.c_str());  // ofstream means Output File Stream.
    //  std::cout << "This is an " << image.size()  << " by " << image[0].size() << " image" << std::endl;
	  for (long i=0; i < long(image.size()); i++)
	  {
	     for (long j=0;  j < long(image[0].size()); j++)
	     {
	       myout << image[i][j] <<" "; 
	     } //end-for-j	   
	     myout << std::endl;      // preferable for c++  (but gives readin problems in read-in in IDL)     	     
	  }; //end for -i  
// create a new starting array, and filling it with random values
myout.close();  
}



void MultiGalaxy::calculate_fft(int Ni, int Nj)
{
   const vector< vector<double> >& psf_image =  Data::get_instance().get_psf();

// Make the psf the same size as the image
mat psf(Ni, Nj);
psf.zeros();

int ni = psf_image.size();
int nj = psf_image[0].size();

int m, n;
for(int i=0; i<ni; i++)
{
m = mod(i - ni/2, Ni);
for(int j=0; j<nj; j++)
{
n =  mod(j - nj/2, Nj);
psf(m, n) = psf_image[i][j];
}
}

fft_of_psf = fft2(psf);
fft_ready = true;
}

void MultiGalaxy::blur_image2(std::vector< std::vector<long double> >  &img) 
{
 
if(!fft_ready)
cerr<<"# Blurring failed."<<endl;

// Copy the image into an Armadillo matrix
mat A(img.size(), img[0].size());
for(size_t i=0; i<img.size(); i++)	
for(size_t j=0; j<img[0].size(); j++)
A(i, j) = img[i][j];

// Do the fft of it
cx_mat B = fft2(A);

// Multiply the two ffts
for(size_t i=0; i<img.size(); i++)
for(size_t j=0; j<img[0].size(); j++)
B(i, j) *= fft_of_psf(i, j);

// Do the inverse fft
B = ifft2(B);

// Put back in img
for(size_t i=0; i<img.size(); i++)
for(size_t j=0; j<img[0].size(); j++)
img[i][j] = real(B(i, j));
}





