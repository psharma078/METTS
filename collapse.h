#pragma once

#include "itensor/util/print_macro.h"

using namespace std;
using namespace itensor;

vector<int>
collapse(MPS& phi,
         Args const& args = Args::global())
    {
    auto sites = siteInds(phi);
    auto N = length(phi);
    //phi = removeQNs(phi);

    auto direction = args.getString("Direction");

    auto cps = vector<int>(N+1);//site indext staets from 1 whereas vector index starts from 0

    phi.position(1);
    for(int j = 1; j <= N; ++j)
        {
	//Index s = removeQNs(sites(j));
        Index sj = removeQNs(sites(j));

	auto PEm = ITensor(sj,prime(sj));
        auto PUp = ITensor(sj,prime(sj));
        auto PDn = ITensor(sj,prime(sj));
	auto PUD = ITensor(sj,prime(sj));
        if(direction == "Z")//this is a projection |z><z|
            {
	    PEm.set(1,1,1.0);
            PUp.set(2,2,1.0);
            PDn.set(3,3,1.0);
	    PUD.set(4,4,1.0);
	    }
        else if(direction == "X")//projection |x><x|
            {
	    PEm.set(1,1,1.0);

            PUp.set(2,2,0.5);
            PUp.set(2,3,0.5);
            PUp.set(3,2,0.5);
	    PUp.set(3,3,0.5);

            PDn.set(2,2,0.5);
            PDn.set(2,3,-0.5);
            PDn.set(3,2,-0.5);
            PDn.set(3,3,0.5);

	    PUD.set(4,4,1.0);
            }
        else Error("Direction '" + direction + "' not recognized");

	Real prob_up = elt(dag(prime(phi(j),"Site"))*PUp*phi(j));
        Real prob_dn = elt(dag(prime(phi(j),"Site"))*PDn*phi(j));
	Real prob_em = elt(dag(prime(phi(j),"Site"))*PEm*phi(j));
	Real prob_ud = elt(dag(prime(phi(j),"Site"))*PUD*phi(j));
        //printfln("probability = ",prob_em," ",prob_up," ", prob_dn," ",prob_ud," ",prob_up+prob_dn+prob_em+prob_ud);

	int st = 0;
        //if(Global::random() > prob_up) st = 2;
        //cps.at(j) = st;
        double rand = Global::random();
	if (rand <= prob_em) st = 1;
	else if (prob_em < rand and rand <= (prob_em + prob_up)) st = 2;
	else if (prob_em + prob_up < rand and rand <= (prob_em + prob_up \
				+ prob_dn)) st = 3;
	else if (prob_em + prob_up + prob_dn < rand and rand <= 1.0) st = 4;
	else Error("probability not conserved");

        cps.at(j) = st;	
 	
	//cout << "rand = " << rand << " st = " << st << endl;

        auto emState = ITensor(sj);
        auto upState = ITensor(sj);
        auto dnState = ITensor(sj);
        auto updnState = ITensor(sj);
        if(direction == "Z")
            {
            emState.set(1,1.0);
            upState.set(2,1.0);
            dnState.set(3,1.0);
            updnState.set(4,1.0);
            }
	else if(direction == "X")
            {
            emState.set(1,1.0);

            upState.set(2,1.0/sqrt(2.0));
            upState.set(3,1.0/sqrt(2.0));

            dnState.set(2,1.0/sqrt(2.0));
            dnState.set(3,-1.0/sqrt(2.0));

            updnState.set(4,1.0);
            }

	//Project state
        //ITensor jstate = (st==1) ? upState : downState;
        ITensor jstate;
        if (st==1) jstate = emState;
        else if (st==2) jstate = upState;
        else if (st==3) jstate = dnState;
        else if (st==4) jstate = updnState;
	else Error("st = "+str(st)+" can not be a cps");
        if(j < N)
            {
            auto newA = phi(j+1)*(dag(jstate)*phi(j));//calculating new A for site j+1
	    //printfln("norm at site ",j, " = ",norm(newA));
            newA /= norm(newA);//dividing by sqrt prob. (or divide by norm is the same thing)
            phi.set(j+1,newA);//replacing site j+1 MPS element
            }
        //Set site j tensor
        phi.set(j,jstate);//site j has cps 1x1 matrix
        }

    return cps;
    }

