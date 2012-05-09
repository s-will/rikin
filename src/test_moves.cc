#include  <stdlib.h>
#include  <string.h>
#include  <stdio.h>
#include  <math.h>
#include <assert.h>

extern "C" {
#include "ViennaRNA/fold_vars.h" // defines global variables
}


// #ifdef _OPENMP
// #include <omp.h>
// #endif


#include <LocARNA/matrices.hh>

#include "unpaired_pf.hh"

#include "hybrid_pf.hh"

#include "hybrid_ensemble_model.hh"

int
main()
{
    // set some global variables for Vienna libRNA
    dangles=2;

    //                      0        1         2         3         4         5         6         7         8         9         0         1         2
    //                      12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
    // const std::string seqA="GUUGGGAACUAGACCGAUCGCCAAUCCGUUUAUCUUUCAUAGAAGCCGGGAUUUAUCAGCUAUGUCGAAGAAUUUUAACUUGCUAUUGGGCACCCUGGUGGGGGUUAGUUUAGUUUUUCCCCAGG";
    //const std::string seqA="CCAAUCCGUUUAUCUUUCAUAGAAGCCGGGAUUUAUCAGCUAUGUCGAAGAAUUUUAACUUGCUA";
    
    
    //                                        UAUCUU...           ...UACA
    //const std::string seqA="CCCCGGGG";
    
    //                      1234567890123456789012345678901
    //const std::string seqA="CGCUACGUACGGAGCUAGCUGAAAC";
    const std::string seqA="ACGUACGGAGCGACGCUACGUACGGAGCUAGCUGAAACGUAGCAGACGUU";
    
    //                      12345678901234567890123456789012
    // const std::string seqB="GAUGCAUGCCUCGAU";
    const std::string seqB="GCUGCGAUGCAUGCCUCGAU";
    
    //const std::string seqB="GGGGGCCCCC";
        
    // ------------------------------------------------------------
    // enumerate moves

    std::cout << "Generate model ..." << std::endl;
    HybEnsModel model(seqA,seqB);
    std::cout << "    DONE." << std::endl;
    
    //    HybEnsModel::StateDescription state(8,9,16,17);
    // HybEnsModel::StateDescription state(3,1,5,3);
     //HybEnsModel::StateDescription state(1,1,3,13);
    //HybEnsModel::StateDescription state(3,1,14,9);
    HybEnsModel::StateDescription state(3,3,4,9,12,17,50,19);
    
    //check state energy
    std::cout << "Energy of state "<<state<<" = "<< model.energy(state) << std::endl;
    
    HybEnsModel::MoveIterator mi(state,model);
    
    size_t count_moves=0;
    size_t count_neighbors=0;

    HybEnsModel::Move *move = mi.firstMove();
    if (move==NULL) {
	std::cout << "No moves!" << std::endl;
    } else {
	
	do {
	    count_moves++;
	    //std::cout << "Try move "; move->print(std::cout); std::cout<<std::endl;
	    
	    HybEnsModel::energy_t tE=move->transitionEnergy();
	    if (tE < 1e6) {
		count_neighbors++;
		std::cout <<count_neighbors << " ("<<count_moves<<") ";
		std::cout << state<<" >>> ";
		move->print(std::cout);
		std::cout << " " << tE << " >>> ";
		HybEnsModel::StateDescription state2=state;
		move->apply(state2);
		std::cout << state2 << "( E="<<model.energy(state2)<<" )"<< std::endl;
		assert(state2.is_valid(model));
	    }
	} while ( (move = mi.nextMove(move))!=NULL );
	std::cout << std::endl;
    }
    
    exit(0);
}
