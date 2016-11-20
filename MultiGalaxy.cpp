#include "MultiGalaxy.h"
#include "DNest4/code/DNest4.h"
#include "Data.h"
#include "DH_Util.h"
#include <string> 

#include <vector> 
#include <sstream>
#include <fstream>                 // includes
#include <iostream>




#include <gsl/gsl_errno.h>



#define ARMA_DONT_USE_CXX11 
#include <armadillo>

/* Mandatory */
#include <profit/profit.h>

//#include <"/home/dhui890/Documents/Research/Profit/libprofit/profit/profit.h>

/* Depending on which profiles you will use... */
//#include <profit/sersic.h>
//#include <profit/ferrer.h>
//#include <profit/moffat.h>
//#include <profit/sky.h>
//#include <profit/psf.h>


#define max_comp 3
#define max_param 15  
// use:   more sample_info.txt 

#include <mutex>          // std::mutex

using namespace std;
using namespace DNest4;
using namespace arma;

// create maximum object with 15 parameters
MultiGalaxy::MultiGalaxy() 
:objects(max_param, max_comp
, false, MyConditionalPriorMultiGalaxy(Data::get_instance().get_ni(), Data::get_instance().get_nj()),  PriorType::log_uniform)
,image(Data::get_instance().get_ni(),vector<long double>(Data::get_instance().get_nj() ))
{
//    set_limit();
    global_index = 0;
    calculate_fft(Data::get_instance().get_ni(), Data::get_instance().get_nj());
    ni  = Data::get_instance().get_ni();
    nj  = Data::get_instance().get_nj();
}


//// A function to set the limit for the flat prior
void MultiGalaxy::set_limit()
{ 
    double delta = 1e-2;
    lower_limit.assign(15,0);
    upper_limit.assign(15,0);
    lower_limit[0] = 0; //x
    lower_limit[1] = 0; //y
    lower_limit[2] = Data::get_instance().get_magzp() - 18; // I
    lower_limit[3] = 1;  // R
    lower_limit[4] = 0.5; // n 
    lower_limit[5] = delta;  // q
    lower_limit[6] = 0;    // theta
    lower_limit[7] = -1.0;  // boxi

    lower_limit[8] = Data::get_instance().get_magzp() - 18; // I
    lower_limit[9] = delta;  // rout
    lower_limit[10] =  delta; // a                      // alpha [ 0, 10]     (Tayuen, Kim paper 2014) 
    lower_limit[11] =  delta;   // b                    // beta [0 ,2] 
    lower_limit[12] = delta;    // q
    lower_limit[13] = 0 ;  // theta  
    lower_limit[14] = -1.0;  // boxi

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
    upper_limit[11] = 2.0 - 10*delta;   // b
    upper_limit[12] = 1;    // q
    upper_limit[13] = 180 ;  // theta  
    upper_limit[14] = 1.0;  // boxi

}

void MultiGalaxy::from_prior(RNG& rng)
{
  std::cout << std::endl << "+++Start from_prior+++" <<   std::endl; 
  do
  {
      objects.from_prior(rng);                           // initiate from_prior in MyConditionalPrior class using the hyper-priors 
                                                         // and from_unif in MyConditionalPrior class using the borders provided by limits. 
  }while (objects.get_components().size() == 0);

   std::cout << " From Prior number of objects:  " << objects.get_components().size() ;
   sigma = exp(log(1.) + log(1E6)*rng.rand());
   calculate_image();         
   std::cout << "+++Finished from prior+++" << std::endl; 

} // End-from prior




double MultiGalaxy::perturb(RNG& rng)
{
        std::cout << " Starting perturb  " << objects.get_components().size() ;
	double logH = 0.;
	if(rng.rand() <= 0.75)
	{
		logH += objects.perturb(rng);
                if (objects.get_components().size() == 0)
                {
                    return -std::numeric_limits<double>::max();
                }
                else
                {
		  calculate_image();
                }
	}
	else
	{
		sigma = log(sigma);
		sigma += log(1E6)*rng.randh();
		sigma = mod(sigma - log(1.), log(1E6)) + log(1.);
		sigma = exp(sigma);
	}
        std::cout << " End perturb  " << objects.get_components().size() ;
	return logH;
}


//double MultiGalaxy::perturb(RNG& rng)
//{
//         std::cout << " Starting perturb  " << objects.get_components().size() ;
//	double logH = 0.;
//        logH += objects.perturb(rng);
//        if (objects.get_components().size() == 0)
//        {
//                return -std::numeric_limits<double>::max();
//        }
////                else
////                {
////		  calculate_image();
////                }


////	if(rng.rand() <= 0.75)
////	{
////		logH += objects.perturb(rng);
////                if (objects.get_components().size() == 0)
////                {
////                    return -std::numeric_limits<double>::max();
////                }
////                else
////                {
////		  calculate_image();
////                }
////	}
////	else
////	{
////		sigma = log(sigma);
////		sigma += log(1E6)*rng.randh();
////		sigma = mod(sigma - log(1.), log(1E6)) + log(1.);
////		sigma = exp(sigma);
////	}
//        calculate_image(); 
//        std::cout << " End perturb  " << objects.get_components().size() ;
//	return logH;
//}





double MultiGalaxy::log_likelihood() const
{
        std::cout << std::endl << "           Started eval loglikelihood; "; 
	const vector< vector<double> >& data = Data::get_instance().get_image();
	const vector< vector<double> >& sig =  Data::get_instance().get_sigma();

	double logL = 0.;
	double var;
        std::cout << "sigma " << sigma << std::endl; 
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

        std::cout << std::endl << "           Finished eval loglikelihood; var:" <<  var << std::endl;    
	logL =0;  // testing the prior
        //std::cout << data.size() << " " << data[0].size() << std::endl;    
	return logL;

}



void MultiGalaxy::calculate_image()
{
        std::cout << std::endl << "           Started calculate image" <<   std::endl; 
        global_index++;
        double x, y;
	// Get coordinatesfrom data and create an empty image
        int ni= Data::get_instance().get_ni();
        int nj= Data::get_instance().get_nj();
	image.assign(Data::get_instance().get_ni(), vector<long double>(Data::get_instance().get_nj(), 0.));

        std::cout << std::endl << "           So far so good1" <<   std::endl; 
        vector<string> titels; 
        titels.assign(15, " ");
        titels[0] = "x"; titels[1] = "y";  titels[2] = "mag"; titels[3] = "Re"; titels[4] = "n"; titels[5] = "q";  titels[6] = "theta"; titels[7] = "boxi";
        titels[8] = "mag-bar"; titels[9] = "Rout";  titels[10] = "a"; titels[11] = "b"; titels[12] = "q-bar"; titels[13] = "theta-bar";  titels[14] = "box-bar"; 
//	titels = ["x", "y", "mag", "re", "n", "q", "theta", "boxi", "mag-bar", "Rour","a","b" "q-bar", "theta-bar", "box-bar"]

        bool update = objects.get_removed().size() == 0;
	// Components
	//const vector< vector<double> >& components = (update)?(objects.get_added())
         //                                        :(objects.get_components());
        std::cout << std::endl << "           So far so good2" <<   std::endl; 
        //bool test = check_in_range();
        //if(!update)
//	{
		// Zero the image
//	image.assign(Data::get_instance().get_ni(),vector<long double>(Data::get_instance().get_nj(), 0.));
//	}

	const vector< vector<double> >& components = objects.get_components();

        std::cout << std::endl << "           Print values in calculate image " << std::endl; 
        std::cout << std::endl << "           Dimensions " << std::endl << "Objects: " << components.size() << " Parameters: " << components[0].size() << std::endl; 

        for(size_t o=0; o<components.size(); ++o)
	{
            std::cout <<            "           ";
            for(size_t i=0; i<components[0].size(); i++)
            {  
                     std::cout << titels[i] << "=" << components[o][i] << " ; ";
            } // end for i-loop
            std::cout << std::endl;

           //profit::Model *model = new profit::Model();
           profit::Model model;
   		
           model.width = ni; 
           model.height = nj;
           model.magzero = Data::get_instance().get_magzp();
           double half_xbin = model.scale_x/2.;
           double half_ybin = model.scale_y/2.;

//           profit::Profile *sersic_profile = model->add_profile("sersic");
            profit::Profile &sersic_profile = model.add_profile("sersic");
           //profit::SersicProfile *sp = static_cast<profit::SersicProfile *>(sersic_profile);
  //         sp->xcen = components[o][0];
    //       sp->ycen = components[o][1];
      //     sp->mag = components[o][2];
        //   sp->re =  components[o][3];
         //  sp->nser = components[o][4];
         //  sp->axrat = components[o][5];
         //  sp->ang = components[o][6];
         //  sp->box = components[o][7];
           sersic_profile.parameter("xcen", components[o][0]);
           sersic_profile.parameter("ycen", components[o][1]);
           sersic_profile.parameter("mag", components[o][2]);
           sersic_profile.parameter("re", components[o][3]);
           sersic_profile.parameter("nser", components[o][4]);
           sersic_profile.parameter("axrat", components[o][5]);
           sersic_profile.parameter("ang", components[o][6]);
           sersic_profile.parameter("box", components[o][7]);


           profit::Profile &ferrer_profile = model.add_profile("ferrer");
           ferrer_profile.parameter("xcen", components[o][0]);
           ferrer_profile.parameter("ycen", components[o][1]);
           ferrer_profile.parameter("mag", components[o][8]);
           ferrer_profile.parameter("rout", components[o][9]*components[o][3]);
           ferrer_profile.parameter("a", components[o][10]);
           ferrer_profile.parameter("b", components[o][11]);
           ferrer_profile.parameter("axrat", components[o][12]*components[o][5]);
           ferrer_profile.parameter("ang", components[o][13]);
           ferrer_profile.parameter("box", components[o][14]);

 //          profit::Profile *ferrer_profile = model->add_profile("ferrer");
 //          profit::FerrerProfile *fp = static_cast<profit::FerrerProfile*>(ferrer_profile);
  //         fp->xcen = components[o][0];
  //        fp->ycen = components[o][1];
  //         fp->mag = components[o][8];
  //         fp->rout = components[o][9]*components[o][3];
  //         fp->a = components[o][10];
  //         fp->b = components[o][11];
   //        fp->axrat = components[o][12]*components[o][5];  // ax_bar < ax_disc
   //        fp->ang = components[o][13];
   //        fp->box = components[o][14];

     
      // std::cout << "rout= " << params[9]*params[3] << " a= " << params[10] << " b= " << params[11] << std::endl;
            std::cout << std::endl << "           Starting to evaluate image " <<   std::endl; 
            std::vector<double> new_image = model.evaluate();	
//	  try {
//                        std::vector<double> image_new = model.evaluate();	
//         } catch (invalid_parameter &e) {
//               cerr << "Oops! There was an error evaluating the model: " << e.what() << endl;
//          }
            std::cout << std::endl << "           Finished image " <<   std::endl; 
            std::cout << std::endl << " x-size / ni " << ni  << " " <<  image.size() << std::endl;   // x-size 
            std::cout << std::endl << " y-size / nj " << nj   << " " <<  image[0].size() << std::endl;   // y-size          
             
        //    std::cout << image.size() << " " << image[0].size() << std::endl;             
        //    std::cout << ni << " " << nj << std::endl;    


            for (int j=0; j < nj; j++) 

            {
	  	x += half_xbin;
	 	y = 0;
             //   std::cout << j;
	 	for (int k=0; k < ni; k++) 
                {
		       y += half_ybin;                       
                       image[k][j] = image[k][j]  +  new_image[k + j*ni];
         	       y += half_ybin;
		} // end-for k-loop
		x += half_xbin;
             } //end-for j-loop
         //    delete model;

           
       } // end-for o-loop
        std::cout << std::endl << "           End calculate image " <<   std::endl; 
     
        std::cout << std::endl << "           Starting to blur " <<   std::endl; 
        blur_image2(image);
        std::cout << std::endl << "           Finished  blur" <<   std::endl; 

        if (global_index % 500 == 0)
        {
             writefile(image); 
        } 
        std::cout << std::endl << "           End calculate image " <<   std::endl; 
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


bool MultiGalaxy::check_in_range()
{
   std::cout << "check in range";
   bool in_range = TRUE;
   const vector< vector<double> >& components =  objects.get_components();
   std::vector<string> titels; 
   titels.assign(15, " ");
   titels[0] = "x"; titels[1] = "y";  titels[2] = "mag"; titels[3] = "Re"; titels[4] = "n"; titels[5] = "q";  titels[6] = "theta"; titels[7] = "boxi";
   titels[8] = "mag-bar"; titels[9] = "Rout";  titels[10] = "a"; titels[11] = "b"; titels[12] = "q-bar"; titels[13] = "theta-bar";  titels[14] = "box-bar"; 

   for(size_t o=0; o<components.size(); ++o)
   {
            std::cout <<            "           ";
            for(size_t i=0; i<components[0].size(); i++)
            {  
                //        std::cout << titels[i] << "=" << components[o][i] << " ; ";
                if ((( i != 2) || (i != 8)) &&  (components[o][i]  < lower_limit[i] || components[o][i] > upper_limit[i] ))
                {
                   in_range = FALSE;
                   std::cout << "OUT OF RANGE " << titels[i] << "=" << components[o][i] << " ; ";
                }//end-if
           } //end for
   } // end-for
return in_range;
} // in_range 



void MultiGalaxy::print(std::ostream& out) const
{	
//         // print image 
//         for (long i=0; i < long(image.size()); i++)
//	  {
//	     for (long j=0;  j < long(image[0].size()); j++)
//	     {
//	       out << image[i][j] <<' '; 
//	     } //end-for-j	   	     
//          } //end for -i    
       out<<setprecision(10);
       out<<sigma<<' ';   
       out<<setprecision(6);
       objects.print(out); out<<' ';
       // store image
        for (long i=0; i < long(image.size()); i++)
	  {
	     for (long j=0;  j < long(image[0].size()); j++)
	     {
	       out << image[i][j] <<' '; 
	     } //end-for-j	   	     
	  } //end for -i  


} // 


string MultiGalaxy::description() const
{
    stringstream s;
//    // Anything printed by MyModel::print (except the last line)   
     s<<" sigma , ";

    // The rest is all what happens when you call .print on an RJObject
    s<<"dim_components, max_num_components, ";

    // Then the hyperparameters (i.e. whatever MyConditionalPrior::print prints)
    s<<"hyper_location, hyper_scale, ";

    // Then the actual number of components
    s<<"num_components, ";

    // Then it's all the components, padded with zeros
    // max_num_components is 100 in this model, so that's how far the
    // zero padding goes.

    const vector< vector<double> >& components =  objects.get_components();
    std::vector<string> titels; 
    titels.assign(15, " ");
    titels[0] = "x"; titels[1] = "y";  titels[2] = "mag"; titels[3] = "Re"; titels[4] = "n"; titels[5] = "q";  titels[6] = "theta"; titels[7] = "boxi";
    titels[8] = "mag-bar"; titels[9] = "Rout";  titels[10] = "a"; titels[11] = "b"; titels[12] = "q-bar"; titels[13] = "theta-bar";  titels[14] = "box-bar"; 


    for(size_t i=0; i<components[0].size(); i++)
    {  
         for(size_t o=0; o<max_comp; ++o)
         {
             s<< titels[i] << ",";
         }
    }  
     
     for (long i=0; i < long(image.size()); i++)
	  {
	     for (long j=0;  j < long(image[0].size()); j++)
             {
	        s<<" img ," ;
	     } //end-for-j	   	     
          } //end for -i 


//    for(int i=0; i<max_comp; ++i)
//        s<<"log_period["<<i<<"], ";
//    for(int i=0; i<max_comp; ++i)
//        s<<"log_amplitude["<<i<<"], ";
//    for(int i=0; i<max_comp; ++i)
//        s<<"phase["<<i<<"], ";
     return s.str();
}


//string MultiGalaxy::description() const
//{
//	return string("objects");
//        return string("Each column is one of the 20 parameters.");
//}



