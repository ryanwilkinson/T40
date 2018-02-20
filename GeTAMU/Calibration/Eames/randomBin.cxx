#include <TRandom.h>

Double_t randomBin()
{
	return gRandom->Uniform(-0.5, 0.5);
}

class randomBinDummy{};
