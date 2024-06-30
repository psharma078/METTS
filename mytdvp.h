#pragma once
#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include "tdvp.h"
#include "basisextension.h"
#include <fstream>

using namespace itensor;

void tdvpEvol(MPS& psi, MPO& H, const Args& param) 
{
    // Start TDVP, either one site or two-site algorithm can be used by adjusting the "NumCenter" argument
    println("-----------------------GSE-TDVP-----------------------------");

    auto sweeps = Sweeps(1);
    sweeps.maxdim() = 500;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 20;

    auto t = param.getReal("tstep", 0.0);
    auto tend = param.getReal("ttotal", 0.0);
    auto n_gse = param.getReal("n_gse",1);
    auto n_center = param.getReal("n_cen",1);
    int nsw = tend / t;

    Real Energy; // Energy from tdvp steps
    for (int n = 1; n <= nsw; ++n)
    {
        if (n <= n_gse) 
	{   
	    for (int k:range(4)){	
            // Global subspace expansion
            std::vector<Real> epsilonK = {1E-12, 1E-12};
            addBasis(psi, H, epsilonK, {"Cutoff", 1E-12,
                                        "Method", "DensityMatrix",
                                        "KrylovOrd", 3,
                                        "DoNormalize", true,
                                        "Quiet", true});
	    } 
        }

        // TDVP sweep, real time evolution
	int numcen;
        if (n<=n_center) numcen=2;
        else numcen=1;
        Energy = tdvp(psi, H, -t, sweeps, {"Truncate", true,
                                                    "DoNormalize", true,
                                                    "Quiet", true,
                                                    "NumCenter", numcen,
                                                    "ErrGoal", 1E-7});

        printfln("\n Imaginary time t = %.10f, tdvp energy = %.10f", t*n , Energy);
    }
}
