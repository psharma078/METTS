#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include "mytdvp.h"
#include "measurement.h"
#include "tebdGates_Hubb.h"
#include "collapse.h"
#include "runningStats.h"

using namespace std;
using namespace itensor;

int 
main(int argc, char* argv[])
{
    if(argc < 2) 
        { 
        printfln("Usage: %s <input_file>", argv[0]); 
        return 0; 
        }
    auto infile = InputFile(argv[1]);
    auto input = InputGroup(infile,"input");

    auto N = input.getInt("N");
    auto t = input.getInt("t");
    auto U = input.getInt("U");
    auto V = input.getInt("V");

    auto beta = input.getReal("beta");
    Real tau_tebd = input.getReal("tau_tebd",0.1);
    Real tau_tdvp = input.getReal("tau_tdvp",0.1);

    auto cutoff = input.getReal("cutoff");
    auto maxdim = input.getInt("maxdim",5000);

    auto nmetts = input.getInt("nmetts",5000);
    auto nwarm = input.getInt("nwarm",5);
    
    auto sites = Electron(N,{"ConserveSz=",true});
     
    Args args;
    args.add("MaxDim",maxdim);
    args.add("Cutoff",cutoff);
    args.add("Method","DensityMatrix");
    args.add("verbose", false);

    auto ampo = AutoMPO(sites);
    for(auto j : range1(N-1))
        {
	ampo += -t,"Cdagup",j,"Cup",j+1;
	ampo += -t,"Cdagup",j+1,"Cup",j;
        ampo += -t,"Cdagdn",j,"Cdn",j+1;
        ampo += -t,"Cdagdn",j+1,"Cdn",j;
	ampo += V, "Ntot",j,"Ntot",j+1;
        }
   for(auto j : range1(N))
    {
        ampo += U,"Nupdn",j;
    } 
    
    auto H = toMPO(ampo);
    //auto H2 = MPO{};
    //nmultMPO(prime(H),H,H2,{"MaxDim",100,"Cutoff",1E-14});
    //H2.mapPrime(2,1);
 
    auto state = InitState(sites,"Emp");
    for(int j : range1(N))
        {
        state.set(j,j%2==0?"Up":"Dn");
        }
    auto psi0 = MPS(state); 
    
    /*auto sweeps = Sweeps(25);
    sweeps.maxdim() = 8,20,50,50,100,100,200,200,200,500,500,500,1000;
    sweeps.cutoff() = 1E-8;
    sweeps.niter() = 2;//change to 4
    sweeps.noise() = 1E-8,1E-9,0.0,0.0;
    //println(sweeps);

    auto [Energy, psi] = dmrg(H,psi0,sweeps,"Quiet");

    printfln("Ground state energy ", Energy);
    printfln("Ground state correlation ");
    auto SS = correlation_SoSi(psi,sites);
    for (int i : range(N))
    {
        printfln("",i+1," ",SS[i]);
    }*/
    
    //initial state for METTS cps construction
    auto sweeps = Sweeps(2);
    auto [E,psi1] = dmrg(H,psi0,sweeps,"Quiet");

    //observables
    bool verbose = true;
    Stats en_stat;//, en2_stat;
    
    size_t count = 0;
    std::vector<double> sum(N, 0.0);
    std::vector<double> sum_of_squares(N, 0.0);

  
    //printfln("total QN before time evolution = ", totalQN(psi)); 
    
    //auto gates = makeGates(input, sites);//making tebd gates
    //gateTEvol(gates, ttotal, tau_tebd, psi0, args);//tebd
    //printfln("TEBD Energy ", inner(psi0,H,psi0));
    //printfln("TEBD maxdim ", maxLinkDim(psi0));

    Args para;
    para.add("tstep", tau_tdvp);
    para.add("ttotal", 0.6);
    para.add("n_gse",0);
    para.add("n_cen",3);
    //printfln("maxdim of initial state for tdvp, ", maxLinkDim(psi));
    tdvpEvol(psi1, H, para);
     
    auto cps = collapse(psi0,{"Direction=","Z"});//starting cps
    
    //-------------------------------------------
    auto ttotal = beta/2.;
    const int nt = int(ttotal/tau_tdvp+(1e-9*(ttotal/tau_tdvp)));
    //cout << "ntot = " << nt << "  " << 1e-9*(ttotal/tau) << endl;
    if(fabs(nt*tau_tdvp-ttotal) > 1E-9)
        {
        Error("Timestep not commensurate with total time");
        }

    para.add("ttotal", ttotal);
    para.add("n_gse",1);
    para.add("n_cen",45);

    //--------------metts steps starts here--------------
    for(int step = 1; step <= (nwarm+nmetts); ++step)
    {
        for(int k : range1(N))
        {
            if (cps.at(k)==1) {state.set(k,"Emp");}
            else
            if (cps.at(k)==2) {state.set(k,"Up");}
	    else 
	    if (cps.at(k)==3) {state.set(k,"Dn");}
            else
	    if (cps.at(k)==4) {state.set(k,"UpDn");}
	    else Error("collapse doesn't executed properly");
	    //printfln("cps: ",cps.at(k));
        }
    	    auto phi = MPS(state);

	    auto ntot_i = expect(phi,sites,"Ntot");
	    auto ntot = std::accumulate(ntot_i.begin(), ntot_i.end(), 0.0);
	    printfln("total particle ", ntot);
	    if (fabs(N-ntot) > 10E-3) 
	    {Error("Total paticle number do not conserve, Ntot = "+std::to_string(ntot));}

        if(verbose)
        {
            if (step <= nwarm) 
                printfln("\nStarting step %d (warmup %d/%d)",step,step,nwarm);
            else
                printfln("\nMaking METTS number %d/%d",step-nwarm,nmetts);
        }
	printfln("total QN before time evolving cps = ", totalQN(phi));
        println("Doing time evolution");

	//gateTEvol(gates, tend_tebd, tau_tebd, phi, args);//first apply tebd
        //printfln("Energy after tebd upto ttotal= ",tend_tebd," is ", inner(phi,H,phi));
	//printfln("Max bond dim after tebd evolution ", maxLinkDim(phi));
	
	tdvpEvol(phi, H, para);

	printfln("total QN after time evolution = ", totalQN(phi));
	//printfln("overlap with initial state: ",inner(psi,phi));

        if(step > nwarm) println("Done making METTS ",step-nwarm);
            
        if(step > nwarm)
        {
            
            //Measure Energy 
            
            const auto en = inner(phi,H,phi);
            en_stat.putin(en);
            auto avgEn = en_stat.avg();
	    auto err_En = en_stat.err();
	    collectdata(step-nwarm,avgEn,err_En,"Average_Energy_N=80_beta=29.6.csv");
            printfln("Energy of METTS %d = %.14f",step-nwarm,en);
            printfln("Average energy = %.14f %.3E",avgEn, err_En);
            printfln("Average energy per site = %.14f %.3E",avgEn/N,err_En/N);

            std::cout<<std::endl;

	    //measure spin correlation
	    int metts_step = step-nwarm;
	    auto corr = correlation_SoSi(phi,sites);
            updateStatistics(corr, sum, sum_of_squares, count);
	    auto mean_corr = calculateMean(sum, count);
            auto std_err_corr = calculateStandardError(sum, sum_of_squares, count);
	    printfln("So.Si at metts step ", count);
	    for (int i:range(N))
	    {
		collectdata(count,i+1,mean_corr[i],std_err_corr[i], "runningAve_corr_N=80_beta=29.6.csv");
		printfln("", i+1, " ", corr[i], " ",std_err_corr[i]);
	    }

        }
        
        // Collapse into product state

        auto dir = (step%2==1) ? "X" : "Z";

        cps = collapse(phi,{"Direction=",dir});

        //Display new product state:
        printf("%s:",dir);
        for(int j = 1; j <= N; ++j)
            {
            //auto state_str = (cps.at(j)==1 ? "+" : "-");
            auto state_str = "NULL";
	    if (cps.at(j)==1)  state_str = "Em";
	    else if (cps.at(j)==2) state_str = "Up";
	    else if (cps.at(j)==3) state_str = "Dn";
	    else if (cps.at(j)==4) state_str = "UpDn";
            print(state_str," ");
            }
        println();	 

    }//metts step ends 
   
    return 0;
}
