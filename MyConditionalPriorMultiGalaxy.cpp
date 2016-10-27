#include "MyConditionalPriorMultiGalaxy.h"
#include "DNest4/code/DNest4.h"
#include <cmath>

using namespace DNest4;

MyConditionalPriorMultiGalaxy::MyConditionalPriorMultiGalaxy(
std::vector<double> lower_limit_, std::vector<double> upper_limit_)
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
 lower_limit.assign(lower_limit_.size(),0);
 upper_limit.assign(upper_limit_.size(),0);

 for(size_t i=0; i<lower_limit_.size(); i++)
 {  
    lower_limit[i] =  lower_limit_[i];
    upper_limit[i] =  upper_limit_[i];
 }
}


// Draw Hyper-Parameters from prior
void MyConditionalPriorMultiGalaxy::from_prior(RNG& rng)
{
        hyper_scale =  3 + 6*rng.rand();                    // scale centre is at about 5, 3-9 
        hyper_location = 20 + 10*rng.rand();                  // location centre is at about 25, 20-30 

//	fluxlim = exp(log(fluxlim_min) + log(fluxlim_max/fluxlim_min)*rng.rand());
//	gamma = rng.rand();

//	radiuslim = exp(log(radiuslim_min) + log(radiuslim_max/radiuslim_min)*rng.rand());
//	gamma_radius = rng.rand();

//	b1 = rng.rand();
//	a1 = rng.rand()*b1;

//	b2 = rng.rand();
//	a2 = rng.rand()*b2;
}

double MyConditionalPriorMultiGalaxy::perturb_hyperparameters(RNG& rng)
{
	double logH = 0.;
        // choose which hyper parameter to change 
	int which = rng.rand_int(2);

	if(which == 0)
	{
                hyper_location += 5*rng.randh();
//		fluxlim = log(fluxlim);
//		fluxlim += log(fluxlim_max/fluxlim_min)*rng.randh();
//		fluxlim = mod(fluxlim - log(fluxlim_min),
//			log(fluxlim_max/fluxlim_min)) + log(fluxlim_min);
//		fluxlim = exp(fluxlim);
	}
	else if(which == 1)
	{
                hyper_scale += 5*rng.randh();
//		gamma += rng.randh();
//		gamma = mod(gamma, 1.);
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

	return logH;
}

// x, y, flux, radius, q, theta, rinner/router, Minner/Mouter

double MyConditionalPriorMultiGalaxy::log_pdf(const std::vector<double>& vec) const
{
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
             vec[14] < lower_limit[14] || vec[15] > upper_limit[14]  )  // boxi-bar

       return -1E300;


       // if inside prior, eval
	double logp = 0.;

        Laplace l1(hyper_location, hyper_scale);
        logp +=   l1.log_pdf(vec[2]) ;
        logp +=   l1.log_pdf(vec[8]) ;

//        logp +=   dt((vec[2]-hyper_location)/hyper_scale, 2,TRUE);
 //       logp +=   dt((vec[8]-hyper_location)/hyper_scale, 2,TRUE);

//	logp += log(alpha) + alpha*log(fluxlim) - (alpha + 1.)*log(vec[2]);
//	logp += log(alpha_radius) + alpha_radius*log(radiuslim)
//			- (alpha_radius + 1.)*log(vec[3]);
//	logp += -log(b1 - a1) - log(b2 - a2);

	return logp;
}

void MyConditionalPriorMultiGalaxy::from_uniform(std::vector<double>& vec) const
{  
   // Laplace distribution
    Laplace l1(hyper_location, hyper_scale);

   for(size_t i=0; i<vec.size(); i++)
   {  
       if (i != 2 ||  i != 8)
       {
            vec[i] = lower_limit[i]  +  (upper_limit[i]  - lower_limit[i])*vec[i];
       }
       else
       {
            vec[i] = l1.cdf_inverse(vec[i]);
       }
   } // end-for
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
   // Laplace distribution
    Laplace l1(hyper_location, hyper_scale);
   for(size_t i=0; i<vec.size(); i++)
   {  
       if (i != 2 ||  i != 8)
       {
            vec[i] = (vec[i] - lower_limit[i]) / (upper_limit[i]  - lower_limit[i]);
       }
       else
       {
            vec[i] = l1.cdf(vec[i]);
       }

   } // end-for
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

