#include "MyConditionalPriorMultiGalaxy.h"
#include "DNest4/code/DNest4.h"
#include <cmath>
#include <vector>

using namespace DNest4;

//MyConditionalPriorMultiGalaxy::MyConditionalPriorMultiGalaxy(
// std::vector<double> lower_limit,  std::vector<double> upper_limit)
MyConditionalPriorMultiGalaxy::MyConditionalPriorMultiGalaxy(double ni, double nj)
: ni(ni), nj(nj)
//: lower_limit[0](lower_limit[0]),upper_limit[0](upper_limit[0]),  
// lower_limit[1](lower_limit[1]),upper_limit[1](upper_limit[1]),  
// lower_limit[2](lower_limit[2]),upper_limit[2](upper_limit[2]),  
// lower_limit[3](lower_limit[3]),upper_limit[3](upper_limit[3]),  
// lower_limit[4](lower_limit[4]),upper_limit[4](upper_limit[4]),  
// lower_limit[5](lower_limit[5]),upper_limit[5](upper_limit[5]),  
// lower_limit[6](lower_limit[6]),upper_limit[6](upper_limit[6]),  
// lower_limit[7](lower_limit[7]),upper_limit[7](upper_limit[7]),  
// lower_limit[8](lower_limit[8]),upper_limit[8](upper_limit[8]),  
// lower_limit[9](lower_limit[9]),upper_limit[9](upper_limit[9]),  
// lower_limit[10](lower_limit[10]),upper_limit[10](upper_limit[10]),  
// lower_limit[11](lower_limit[11]),upper_limit[11](upper_limit[11]),  
// lower_limit[12](lower_limit[12]),upper_limit[12](upper_limit[12]),  
// lower_limit[13](lower_limit[13]),upper_limit[13](upper_limit[13]),  
// lower_limit[14](lower_limit[14]),upper_limit[14](upper_limit[14]))  
{
// std::cout << std::endl << "Make Conditional Prior object" <<   std::endl; 

//  lower_limit = new vector<int>(lower_limit_.size(),0);
  //uppper_limit = new vector<int>(lower_limit_.size(),0);

 //lower_limit.assign(lower_limit_.size(),0);
 //upper_limit.assign(upper_limit_.size(),0);
 set_upper_limit();
 set_lower_limit();

 for(size_t i=0; i<lower_limit.size(); i++)
 {  
     std::cout <<  lower_limit[i] << " " << upper_limit[i];
  // // lower_limit[i] =  lower_limit_[i];
   //// upper_limit[i] =  upper_limit_[i];
 }
 std::cout << std::endl << "Finished init Conditional Prior object" <<   std::endl; 
}


// Draw Hyper-Parameters from prior
void MyConditionalPriorMultiGalaxy::from_prior(RNG& rng)
{
        std::cout << std::endl << "---Start from_prior---" <<   std::endl; 
        hyper_location = 25 + 10*rng.randn();                  // location centre is at about 25, 20-30 
        hyper_scale =  5*rng.rand();                           // scale centred  

     //   wrap(hyper_scale, 1e-2, 10);  
//	fluxlim = exp(log(fluxlim_min) + log(fluxlim_max/fluxlim_min)*rng.rand());
//	gamma = rng.rand();

//	radiuslim = exp(log(radiuslim_min) + log(radiuslim_max/radiuslim_min)*rng.rand());
//	gamma_radius = rng.rand();

//	b1 = rng.rand();
//	a1 = rng.rand()*b1;

//	b2 = rng.rand();
//	a2 = rng.rand()*b2;
        std::cout << std::endl << "Just generated hyper priors: " <<  hyper_scale <<  " and  "  << hyper_location <<   std::endl; 
        std::cout << "---End from_prior---" <<   std::endl; 
}

double MyConditionalPriorMultiGalaxy::perturb_hyperparameters(RNG& rng)
{
        std::cout << std::endl << "Start perturb" <<   std::endl; 
	double logH = 0.;
        // choose which hyper parameter to change 
	int which = rng.rand_int(2);

	if(which == 0)
	{
                hyper_location += 3*rng.randn();
//		fluxlim = log(fluxlim);
//		fluxlim += log(fluxlim_max/fluxlim_min)*rng.randh();
//		fluxlim = mod(fluxlim - log(fluxlim_min),
//			log(fluxlim_max/fluxlim_min)) + log(fluxlim_min);
//		fluxlim = exp(fluxlim);
	}
	else if(which == 1)
	{
                hyper_scale += rng.randh();
		hyper_scale = mod(hyper_scale, 10.);
	}
//	else if(which == 2)
//	{
//		radiuslim = log(radiuslim);
//		radiuslim += log(radiuslim_max/radiuslim_min)*rng.randh();
//		radiuslim = mod(radiuslim - log(radiuslim_min),
//			log(radiuslim_max/radiuslim_min)) + log(radiuslim_min);
//		radiuslim = exp(radiuslim);
//	}
//	else if(which == 3)
//	{
//		gamma_radius += rng.randh();
//		gamma_radius = mod(gamma_radius, 1.);
//	}
//	else if(which == 4)
//	{
//		a1 /= b1;
//		b1 += rng.randh();
//		b1 = mod(b1, 1.);
//		a1 *= b1;
//	}
//	else if(which == 5)
//	{
//		a1 += b1*rng.randh();
//		a1 = mod(a1, b1);
//	}
//	else if(which == 6)
//	{
//		a2 /= b2;
//		b2 += rng.randh();
//		b2 = mod(b2, 1.);
//		a2 *= b2;
//	}
//	else
//	{
//		a2 += b2*rng.randh();
//		a2 = mod(a2, b2);
//	}
        std::cout << std::endl << "End perturb" <<   std::endl; 
	return logH;
}

// x, y, flux, radius, q, theta, rinner/router, Minner/Mouter

double MyConditionalPriorMultiGalaxy::log_pdf(const std::vector<double>& vec) const
{
        std::cout << std::endl << "Start  logpdf" <<   std::endl; 
        std::cout << "Vecsize " << vec.size(); 

        // if outside of prior range log prior = -infty 
	if(
             vec[0] < lower_limit[0] || vec[0] > upper_limit[0] ||   // x
             vec[1] < lower_limit[1] || vec[1] > upper_limit[1] ||   // y
             // vec[2] < lower_limit[2] || vec[2] > upper_limit[2] ||   // I
             vec[3] < lower_limit[3] || vec[3] > upper_limit[3] ||   // R 
             vec[4] < lower_limit[4] || vec[4] > upper_limit[4] ||   // n 
             vec[5] < lower_limit[5] || vec[5] > upper_limit[5] ||    //q 
             vec[6] < lower_limit[6] || vec[6] > upper_limit[6] ||   // theta
             vec[7] < lower_limit[7] || vec[7] > upper_limit[7] ||   // boxi
             //vec[8] < lower_limit[8] || vec[8] > upper_limit[8] ||   // I-bar
             vec[9] < lower_limit[9] || vec[9] > upper_limit[9] ||    // Rout
             vec[10] < lower_limit[10] || vec[10] > upper_limit[10] ||  // a
             vec[11] < lower_limit[11] || vec[11] > upper_limit[11] ||  // b
             vec[12] < lower_limit[12] || vec[12] > upper_limit[12] ||   // q-bar
             vec[13] < lower_limit[13] || vec[13] > upper_limit[13] ||    // theta-bar
             vec[14] < lower_limit[14] || vec[14] > upper_limit[14] )  // boxi-bar
        return -1E300;
             
       // if inside prior, eval
	double logp = 0.;

        std::cout << std::endl << "logp" << logp << std::endl; 
        Laplace l1(hyper_location, hyper_scale);
        logp +=   l1.log_pdf(vec[2]) ;
        std::cout << std::endl << "logp" << logp << std::endl; 
        logp +=   l1.log_pdf(vec[8]) ;

        std::cout << std::endl << "logp" << logp << std::endl; 

//        logp +=   dt((vec[2]-hyper_location)/hyper_scale, 2,TRUE);
 //       logp +=   dt((vec[8]-hyper_location)/hyper_scale, 2,TRUE);

//	logp += log(alpha) + alpha*log(fluxlim) - (alpha + 1.)*log(vec[2]);
//	logp += log(alpha_radius) + alpha_radius*log(radiuslim)
//			- (alpha_radius + 1.)*log(vec[3]);
//	logp += -log(b1 - a1) - log(b2 - a2);

      std::cout << std::endl << "Print values in log_pdf" << std::endl; 
      for(size_t i=0; i<vec.size()-1; i++)
      {  
        std::cout << vec[i] << " ";
      }
      std::cout << std::endl << "End logpdf,  return logp"  << std::endl; 
	return logp;
}

void MyConditionalPriorMultiGalaxy::from_uniform(std::vector<double>& vec) const
{  
    std::cout << std::endl << "---Start from_unif---" <<   std::endl; 
  // print parameters
    std::vector<string> titels; 
    titels.assign(15, " ");
    titels[0] = "x"; titels[1] = "y";  titels[2] = "mag"; titels[3] = "Re"; titels[4] = "n"; titels[5] = "q";  titels[6] = "theta"; titels[7] = "boxi";
    titels[8] = "mag-bar"; titels[9] = "Rout";  titels[10] = "a"; titels[11] = "b"; titels[12] = "q-bar"; titels[13] = "theta-bar";  titels[14] = "box-bar"; 

//    std::cout << "Vecsize " << vec.size() <<   std::endl;  
    // Laplace distribution
//    std::cout << "Hyper parameters "  <<  hyper_location << " " <<  hyper_scale <<   std::endl; 
    Laplace l1(hyper_location, hyper_scale);         
    std::cout << "Done hyper  parameters "  <<  hyper_location << " " <<  hyper_scale <<   std::endl;
    std::cout << "Check Vec parameters "  <<  hyper_location << " " <<  hyper_scale <<   std::endl;  

  for(size_t i=0; i<vec.size(); i++)
    {  
                     std::cout << titels[i] << "="  << "  " <<  vec[i] << "   "  << std::endl;
    }

   // Draw parameters 
   for(size_t i=0; i<vec.size(); i++)
   {  
       if(i == 2 || i == 8)
       {
            vec[i] = l1.cdf_inverse(vec[i]);
       }
       else
       {
           vec[i] = lower_limit[i]  +  (upper_limit[i]  - lower_limit[i])*vec[i];
       }

   } // end-for



  //   std::cout << "Print limit arrays" ; 
//     std::cout << lower_limit[0] ; 

     for(size_t i=0; i<vec.size(); i++)
    {  
               //     std::cout << titels[i] << "="  << vec[i] << std::endl;
                     std::cout << titels[i] << "=" << lower_limit[i]  << " < " <<  vec[i] << " <  " << upper_limit[i] << std::endl;
    }
  
  std::cout <<  "---End from_unif---" <<   std::endl; 
}


//	vec[0] = x_min + (x_max - x_min)*vec[0];
//	vec[1] = y_min + (y_max - y_min)*vec[1];
//	vec[2] = fluxlim*pow(1. - vec[2], -gamma);
//	vec[3] = radiuslim*pow(1. - vec[3], -gamma_radius);
//	vec[4] = 0.2 + 0.8*vec[4];
//	vec[5] = M_PI*vec[5];
//	vec[6] = a1 + (b1 - a1)*vec[6];
//	vec[7] = a2 + (b2 - a2)*vec[7];
//}

void MyConditionalPriorMultiGalaxy::to_uniform(std::vector<double>& vec) const
{
    std::cout << std::endl << "Start to_unif" <<   std::endl; 
   // Laplace distribution
    Laplace l1(hyper_location, hyper_scale);
   for(size_t i=0; i<vec.size(); i++)
   {  
       if (i != 2 ||  i != 8)
       {
            vec[i] = l1.cdf(vec[i]);
           
       }
       else
       {
            vec[i] = (vec[i] - lower_limit[i]) / (upper_limit[i]  - lower_limit[i]);
       }
     if ((vec[i] <0 ) ||  (vec[i] > 1 ) )
     {
           std::cout << std::endl << "Error in to_unif" <<   std::endl; 
           exit(1);
      }
   } // end-for
    std::cout << std::endl << "End to_unif" <<   std::endl; 
}




//{
//	double alpha = 1./gamma;
//	double alpha_radius = 1./gamma_radius;

//	vec[0] = (vec[0] - x_min)/(x_max - x_min);
//	vec[1] = (vec[1] - y_min)/(y_max - y_min);
//	vec[2] = 1. - pow(fluxlim/vec[2], alpha);
//	vec[3] = 1. - pow(radiuslim/vec[3], alpha_radius);
//	vec[4] = (vec[4] - 0.2)/0.8;
//	vec[5] = vec[5]/M_PI;
//	vec[6] = (vec[6] - a1)/(b1 - a1);
//	vec[7] = (vec[7] - a2)/(b2 - a2);
//}





void MyConditionalPriorMultiGalaxy::print(std::ostream& out) const
{
       out <<  hyper_location << ' ' <<  hyper_scale << ' ';
//	out<<vec[0]<<' '<<vec[1]<<' '<<vec[2]<<' '<<vec[3]<<' ';
//	out<<vec[4]<<' '<<vec[5]<<' '<<vec[6]<<' '<<vec[7]<<' ';
//	out<<vec[8]<<' '<<vec[9]<<' '<<vec[10]<<' '<<vec[11]<<' ';
//	out<<vec[12]<<' '<<vec[13]<<' '<<vec[14]<<' ';
}


// A function to set the limit for the flat prior
void MyConditionalPriorMultiGalaxy::set_lower_limit()
{ 
    double delta = 1e-2;
    lower_limit.assign(15,0);
    lower_limit[0] = 0; //x
    lower_limit[1] = 0; //y
    lower_limit[2] =  - 18; // I
    lower_limit[3] = 1;  // R
    lower_limit[4] = 0.5; // n 
    lower_limit[5] = delta;  // q
    lower_limit[6] = 0;    // theta
    lower_limit[7] = -1.0;  // boxi

    lower_limit[8] = - 18; // I
    lower_limit[9] = delta;  // rout
    lower_limit[10] =  delta; // a                      // alpha [ 0, 10]     (Tayuen, Kim paper 2014) 
    lower_limit[11] =  delta;   // b                    // beta [0 ,2] 
    lower_limit[12] = delta;    // q
    lower_limit[13] = 0 ;  // theta  
    lower_limit[14] = -1.0;  // boxi
}


// A function to set the limit for the flat prior
void MyConditionalPriorMultiGalaxy::set_upper_limit()
{ 
    double delta = 1e-2;
    upper_limit.assign(15,0);
    upper_limit[0] = ni; //x
    upper_limit[1] = nj; //y
    upper_limit[2] = 20; // I
    upper_limit[3] = 50;  // R
    upper_limit[4] = 10; // n 
    upper_limit[5] = 1;  // q
    upper_limit[6] = 180;    // theta
    upper_limit[7] = 1.0;  // boxi

    upper_limit[8] = 20; // I
    upper_limit[9] = 3;  // rout
    upper_limit[10] = 10; // a
    upper_limit[11] = 2.0 - 10*delta;   // b
    upper_limit[12] = 1;    // q
    upper_limit[13] = 180 ;  // theta  
    upper_limit[14] = 1.0;  // boxi

}


