#pragma once
#include "itensor/util/print_macro.h"

using std::vector;
using namespace itensor;

vector<BondGate>
makeGates(auto input, const SiteSet sites)
{
 
    auto N = input.getInt("N");
    //auto t = input.getInt("t");
    auto U = input.getInt("U");
    auto V = input.getInt("V");

    auto tstep = input.getReal("tau_tebd",0.001);

   vector<BondGate> gates;
   //forward gate
   for (int b = 1; b<=N-1; ++b)
   {
	   auto hterm = -1.*sites.op("Adagup*F",b)*sites.op("Aup",b+1);
	   hterm += -1.*sites.op("Adagdn",b)*sites.op("F*Adn",b+1);
	   hterm += +1.*sites.op("Aup*F",b)*sites.op("Adagup",b+1);
	   hterm += +1.*sites.op("Adn",b)*sites.op("F*Adagdn",b+1);
	   hterm += V*sites.op("Ntot",b)*sites.op("Ntot",b+1);

	   //hterm += U * sites.op("Nupdn",b)*sites.op("Id",b+1);
           //if(b==N-1) hterm += U*sites.op("Id",b)*sites.op("Nupdn",b+1);

	   auto g = BondGate(sites, b, b+1, BondGate::tImag, tstep/2., hterm);
           gates.push_back(g);

   }
   //reverse gate
   for (int b = N - 1; b >= 1; --b)
   {
	   auto hterm = -1.*sites.op("Adagup*F",b)*sites.op("Aup",b+1);
           hterm += -1.*sites.op("Adagdn",b)*sites.op("F*Adn",b+1);
           hterm += +1.*sites.op("Aup*F",b)*sites.op("Adagup",b+1);
           hterm += +1.*sites.op("Adn",b)*sites.op("F*Adagdn",b+1);
	   hterm += V*sites.op("Ntot",b)*sites.op("Ntot",b+1);

	   //if(b==N-1) hterm += U*sites.op("Id",b)*sites.op("Nupdn",b+1);
	   //hterm += U * sites.op("Nupdn",b)*sites.op("Id",b+1);

           auto g = BondGate(sites, b, b+1, BondGate::tImag, tstep/2., hterm);
           gates.push_back(g);
   }

   //onsite term
   for (int b = 1; b <= N - 1; ++b)
    {

        auto hterm = U * op(sites, "Nupdn", b) * op(sites, "Id", b + 1);
        hterm += U * op(sites, "Nupdn", b + 1) * op(sites, "Id", b);

        if (b == 1) //for boundary,as for boundary, gates will be applied only once per site
        {
            hterm += U * op(sites, "Id", b + 1) * op(sites, "Nupdn", b);
        }
        if (b == N - 1)
        {
            hterm += U * op(sites, "Id", b) * op(sites, "Nupdn", b + 1);
        }
	
        auto g = BondGate(sites, b, b+1, BondGate::tImag, tstep / 2., hterm);
        gates.push_back(g);
    }

   return gates;
}

