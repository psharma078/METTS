#include "itensor/all.h"
#include "tebdGates_Hubb.h"
#include "Sz2.h"
#include "collapse.h"
#include "itensor/util/print_macro.h"

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

    auto beta = input.getReal("beta");
    Real tau = input.getReal("tau",0.1);

    auto cutoff = input.getReal("cutoff");
    auto maxdim = input.getInt("maxdim",5000);

    auto nmetts = input.getInt("nmetts",50000);
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
        }
   for(auto j : range1(N))
    {
        ampo += U,"Nupdn",j;
    }

    auto H = toMPO(ampo);
    auto H2 = MPO{};
    nmultMPO(prime(H),H,H2,{"MaxDim",100,"Cutoff",1E-14});
    H2.mapPrime(2,1);

    auto state = InitState(sites,"Emp");
    for(int j : range1(N))
        {
        auto st = (j%2==1 ? "Up" : "Up");
        state.set(j,st);
        }
    //state.set(N,"UpDn");
    auto psi0 = MPS(state);

    auto Sz2 = makeSz2(input,sites);
    auto sweeps = Sweeps(2);
    sweeps.maxdim() = 4,8;
    sweeps.cutoff() = 1E-7;
    sweeps.niter() = 2;//change to 4
    sweeps.noise() = 1E-4;
    //println(sweeps);

    auto [ENergy, psi] = dmrg(H,psi0,sweeps,"Quiet");

   //observables
    bool verbose = true;
    Stats en_stat, en2_stat, sz2_stat, cv_stat;

    auto ttotal = beta/2.;
    const int nt = int(ttotal/tau+(1e-9*(ttotal/tau)));
    //cout << "ntot = " << nt << "  " << 1e-9*(ttotal/tau) << endl;
    if(fabs(nt*tau-ttotal) > 1E-9)
        {
        Error("Timestep not commensurate with total time");
        }

    printfln("total QN before tebd = ", totalQN(psi));
    auto gates = makeGates(input, sites);//making tebd gates
    gateTEvol(gates, ttotal, tau, psi, args);//tebd
    printfln("total QN after tebd = ", totalQN(psi));
    auto cps = collapse(psi,{"Direction=","Z"});

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

     if(verbose)
            {
            if (step <= nwarm)
                printfln("\nStarting step %d (warmup %d/%d)",step,step,nwarm);
            else
                printfln("\nMaking METTS number %d/%d",step-nwarm,nmetts);
            }
        printfln("total QN before tebd = ", totalQN(phi));
        println("Doing time evolution");
        //for(int tt = 1; tt <= nt; ++tt)
        //    {
        //          gateTEvol(gates, tau, tau, psi, args);//tebd
        //    }

        gateTEvol(gates, ttotal, tau, phi, args);//tebd
        printfln("total QN after tebd = ", totalQN(phi));
        printfln("overlap with initial state: ",inner(psi,phi));
        printfln("Max bond dim after time evolution: %d",maxLinkDim(phi));

        if(step > nwarm) println("Done making METTS ",step-nwarm);

    if(step > nwarm)
            {
            //
            //Measure Energy
            //
            const auto en = inner(phi,H,phi);
            en_stat.putin(en);
            auto avgEn = en_stat.avg();
            printfln("Energy of METTS %d = %.14f",step-nwarm,en);
            printfln("Average energy = %.14f %.3E",avgEn,en_stat.err());
            printfln("Average energy per site = %.14f %.3E",avgEn/N,en_stat.err()/N);
            //
            //Measure CV
            //
            const auto en2 = inner(phi,H2,phi);
            en2_stat.putin(en2);
            auto avgE2 = en2_stat.avg();
            printfln("Avg E2  = %.14f %.3E", avgE2,en2_stat.err());
            //auto cv = (en2 - en*en) * beta*beta;
            //cv_stat.putin(cv);
            //auto aCv = cv_stat.avg();
            //auto eCv = cv_stat.err();
            //printfln("CV/N = %.14f %.3E (%.5f,%.5f)",aCv/N, eCv/N,(aCv-eCv)/N,(aCv+eCv)/N);
            auto Cv = (avgE2 - (avgEn*avgEn))*beta*beta;
            //auto eCv = (en2_stat.err() - en_stat.err()*en_stat.err())*beta*beta;
            cv_stat.putin(Cv);
            auto aCv = cv_stat.avg();
            auto eCv = cv_stat.err();

            printfln("CV/N = %.14f %.3E (%.5f,%.5f)",aCv/N, eCv/N,(aCv-eCv)/N,(aCv+eCv)/N);
            //
            //Measure Susceptibility
            //
            auto sz2val = inner(phi,Sz2,phi);
            sz2_stat.putin(sz2val);
            auto asus = (sz2_stat.avg()*beta);
            auto esus = (sz2_stat.err()*beta);
            printfln("<Sz^2> for METTS %d = %.14f",step-nwarm,sz2val);
            printfln("Average susceptibility per site = %.14f %.3E (%.5f,%.5f)",
                     asus/N,esus/N,(asus-esus)/N,(asus+esus)/N);

            }

        //
        // Collapse into product state
        //

        auto dir = (step%2==1) ? "X" : "X";

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

        /*auto ntotal = 0.0;
        for (int i : range1(N)){
        auto ntot = sites.op("Ntot",i);
        auto sz = sites.op("Sz",i);
        psi.position(i);
        auto C = psi(i)*ntot;
        auto Cz = psi(i)*sz;
        C *= dag(prime(psi(i),"Site"));
        Cz *= dag(prime(psi(i),"Site"));
        ntotal += elt(C);
        printfln("site ",i," <ntot>= ",elt(C), " <Sz> = ", elt(Cz));
        }
        printfln("Ntotal = ", ntotal);*/

        }//metts step ends

    return 0;
    }
                              
