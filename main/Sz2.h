#pragma once
#include "itensor/util/print_macro.h"

using std::vector;
using namespace itensor;

MPO
makeSz2(auto input, const SiteSet sites)
{
    auto N = input.getInt("N");
    auto sz = AutoMPO(sites);
    for (auto j : range1(N))
    {
            sz += 1.0, "Sz", j;
    }
    auto Sz = toMPO(sz);
    auto Sz2 = MPO{};
    nmultMPO(prime(Sz),Sz,Sz2,{"MaxDim",1000,"Cutoff",1E-14});
    Sz2.mapPrime(2,1);
   return Sz2;
}
