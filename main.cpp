#include <iostream>
#include "DNest4/code/DNest4.h"
#include "MultiGalaxy.h"
#include "Data.h"


using namespace std;
using namespace DNest4;

int main(int argc, char** argv)
{
	Data::get_instance().load("Data/test_metadata.txt", "Data/fitim.fits", "Data/sigma.fits","Data/psfim.fits");

	Sampler<MultiGalaxy> sampler = setup<MultiGalaxy>(argc, argv);

	sampler.run();

	return 0;
}






