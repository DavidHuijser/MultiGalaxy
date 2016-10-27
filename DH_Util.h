
#ifndef DH_Util5_H
#define DH_Util5_H


// gcc -Wall -c DH_Util.cpp
// #include "RJ_Util5.h"
#include <sstream>
#include <string> 
#include <fstream>                 // includes
#include <iostream>                // includes:  cerr
#include <CCfits/CCfits>
#define PI 3.14159265352982
#define MATHLIB_STANDALONE 1

#include <Rmath.h>

struct OnePoint
{                                                     // Structure to store all info on one typical data-point
                                                      // each point is defined by the parameters,     
    std::vector< std::vector<double> >  parameters;   // parameters;    (first line for the first component)
                                                      // each line has several components:     [model, .............]
                                                      //
    std::vector<double> object_types;                 // object[(0= sersic, ................, -1=empty] 
    std::vector<double> object_parameters;             // number of parameters per object [8=sersic,                        , 0= empty]
    double logprior, loglikelihood;                   // logprior and loglikelihood 
    double sigma;
    int objects;
    bool in_prior_range; 
//    unsigned int comp;                                // number of components
//    std::vector<double>  global;                      // parameters;    (first line for the first component)
 //   unsigned int n;                                 // global [sigma, objects]             
}; // Onepoint


// FUNCTIONS PROTO_TYPE
double logsubstCpp(double x, double y); 

double DH_pTruncNormal(double x, double mean, double sigma, double a, double b, bool logp);

double DH_drandom(double lower_limit, double upper_limit);  //uniform
double DH_prandom(double val, double lower_limit, double upper_limit, bool logp);      //uniform
double drandom(double lower_limit, double upper_limit);
// Normal distribution 
double DH_dnorm(void);
double DH_pnorm(double x, double mean, double sd);
// double dnormal(void);

// void writefile(string filename,const  std::vector< std::vector<double> > &image);
//int mod(int y, int x);
//double mod(double y, double x);
std::string IntToStr(int n); 
double DH_dcauchy(double location, double scale);
double DH_pcauchy(double x, double location, double scale );
double DH_dTruncNormal(double mean, double sd, double a, double b);               // draw from a truncate normal

double DH_pt(double x, double gamma, bool logp);
double DH_dt(double gamma, bool  logp);

  // T-distribution
//   double draw_tdistribution(double gamma, double location, double scale);
 //  double prob_tdistribution(double x, double gamma, double location, double scale);


#endif
