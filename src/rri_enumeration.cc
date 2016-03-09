#include "hybrid_ensemble_model.hh"
#include "rri_enumeration.hh"

void
RRIEnumeration::write_state(std::ostream &out,
			    double energy, 
			    const HybEnsModel::StateDescription &state, bool binary) {
    if (binary) {
	out.precision(6);
	out.setf(std::ios_base::fixed);
	out << energy << '\t';
	
	state.write_binary(out);
	out.put(0);
    } else {
	out.setf(std::ios_base::fixed);
	out.precision(6);
	out << energy;
	
	for (size_t i=0; i<state.size(); i++) {
	    out << '\t' << state[i].i1
		 << '\t' << state[i].i2
		 << '\t'  << state[i].j1
		 << '\t'  << state[i].j2;
	}
	out << '\n';
    }
}

bool 
RRIEnumeration::read_state(std::istream &in, 
			    HybEnsModel::StateDescription &state,
			    double &energy, 
			    size_t lineno,
			    bool binary) {
    if(binary) {
	in >> energy;
	if (in.eof()) return false;
	
	if(in.get()!=' ') {
	    std::cerr << "expected blank after energy at line "<< lineno <<std::endl;
	    return false; 
	}
	
	if ( !state.read_binary(in) || in.get()!=0) { 
	    std::cerr << "Error while parsing state at "<< lineno <<"."<<std::endl;
	    return false; 
	}
	
    }else{
	std::string line;
	
	if (!getline(in,line)) return false;
	
	std::istringstream linein(line);
	linein >> energy;
	
	std::vector<size_t> state_vec;
	for (size_t i; linein >> i;) {
	    state_vec.push_back(i);
	}
	
	switch ( state_vec.size() ) {
	case 0:
	    state = HybEnsModel::StateDescription();
	    break;
	    
	case 4:
	    state = HybEnsModel::StateDescription(state_vec[0],
						  state_vec[1],
						  state_vec[2],
						  state_vec[3]);
	    break;
	case 8:
	    state = HybEnsModel::StateDescription(state_vec[0],
						  state_vec[1],
						  state_vec[2],
						  state_vec[3],
						  state_vec[4],
						  state_vec[5],
						  state_vec[6],
						  state_vec[7]);
	    break;
	default:
	    std::cerr << "ERROR: invalid input line "<< lineno<<"." << std::endl;
	    exit(-1);
	}
	
    }

    return true;
}


// void
// RRIEnumeration::register_state(double energy, 
// 			       const HybEnsModel::StateDescription &state) {
//     states_.push_back(State(energy,state));
// }
    
size_t
RRIEnumeration::enumerate_single_sites(bool in_memory,std::ostream &out) const {
    size_t count_single_states=0;
    
    //cout << "Single Hybridisations:" << endl;
    for (size_t i1=model_.region_startA(); i1<=model_.region_endA(); i1++) {
	for (size_t i2=model_.region_startB(); i2<=model_.region_endB(); i2++) {
	    if ( (model_.pair_type(i1,i2)==0) ) {
		continue;
	    }
	    // enumerate j1 s.t. site length is between minimum site size and maximum hybridization site length
	    for (size_t j1=i1+model_.minsitesize()-1; j1<=std::min(model_.maxsitesize()+i1-1,model_.region_endA()); j1++) {

		// enumerate j2 s.t. site length is between minimum site size and maximum hybridization site length
		// and, furthermore, the maximum hybridization site length difference is not exceeded
		size_t from_j2 = std::max(i2+model_.minsitesize()-1+model_.maxsitesize_diff(),j1-i1+i2)-model_.maxsitesize_diff();
		size_t to_j2   = std::min(model_.maxsitesize()+i2-1,model_.region_endB());
		to_j2 = std::min(to_j2,j1-i1+i2+model_.maxsitesize_diff());
		for (size_t j2=from_j2; j2<=to_j2; j2++) {
		    if ( (model_.pair_type(j1,j2))==0 ) {
			continue;
		    }
		    
		    HybEnsModel::StateDescription state(i1,i2,j1,j2);
		    		    
		    double energy_hyb = 
			model_.energy_hybrid(state[0]);
		    
		    if ( energy_hyb > max_hyb_energy_ ) continue;
		    
		    // double energy_unpair =
		    // 	model_.energy_unpair(state[0]);
		    
		    // double total_energy = 
		    // 	energy_hyb
		    // 	+ energy_unpair;
		    
		    double total_energy=model_.energy(state);
		    
		    if (total_energy <= max_total_energy_) {
			write_state(out,total_energy,state,binary_);
			count_single_states++;
		    }
		}
	    }
	}
    }
    return count_single_states;
}


size_t
RRIEnumeration::enumerate_double_sites(bool in_memory,std::ostream &out) const {
    
    // Indexing for double hybridization 
    // ----\        /--------\       /---------
    //     i1------j1        k1-----l1
    //     i2------j2        k2-----l2
    //  ---/        \-------/         \------------
    
    size_t count_double_states=0;

    const std::string &seqA = model_.seqA();
    const std::string &seqB = model_.seqB();
	
    //cout << "Double Hybridizations:" << endl;
    for (size_t i1=model_.region_startA(); i1<=model_.region_endA(); i1++) {
	
	//double progress = (i1/(double)(model_.region_endA()-model_.region_startA()+1));
	//if (verbose) std::cerr << "\r" << (int(progress*10000)/100.0) << " %   ";

	for (size_t i2=model_.region_startB(); i2<=model_.region_endB(); i2++) {
	    if ( (model_.pair_type(i1,i2)==0) ) {
		continue;
	    }
	    
	    for (size_t j1=i1+model_.minsitesize()-1; j1<=std::min(model_.maxsitesize()+i1-1,model_.region_endA()); j1++) {

		size_t from_j2 = std::max(i2+model_.minsitesize()-1+model_.maxsitesize_diff(),j1-i1+i2)-model_.maxsitesize_diff();
		size_t to_j2   = std::min(model_.maxsitesize()+i2-1,model_.region_endB());
		to_j2 = std::min(to_j2,j1-i1+i2+model_.maxsitesize_diff());
		
		for (size_t j2=from_j2; j2<=to_j2; j2++) {
		    
		    if ( (model_.pair_type(j1,j2)==0) ) {
			continue;
		    }
		    
		    HybEnsModel::StateDescription::ISite is1(i1,i2,j1,j2);
		    
		    double energy_hyb1 = model_.energy_hybrid(is1);
		    
		    if ( energy_hyb1 > max_hyb_energy_ ) continue;
		    
		    for (size_t k1=j1+model_.minsitedist()+1; k1<=model_.region_endA(); k1++) {
			for (size_t k2=j2+model_.minsitedist()+1; k2<=model_.region_endB(); k2++) {
			    if ( (model_.pair_type(k1,k2)==0) ) {
				continue;
			    }
			    
			    for (size_t l1=k1+model_.minsitesize()-1; l1<=std::min(model_.maxsitesize()+k1-1,model_.region_endA()); l1++) {
				
				size_t from_l2 = std::max(k2+model_.minsitesize()-1+model_.maxsitesize_diff(),l1-k1+k2)-model_.maxsitesize_diff();
				size_t to_l2   = std::min(model_.maxsitesize()+k2-1,model_.region_endB());
				to_l2 = std::min(to_l2,l1-k1+k2+model_.maxsitesize_diff());
		
				for (size_t l2=from_l2; l2<=to_l2; l2++) {
				    if ( (model_.pair_type(l1,l2)==0) ) {
					continue;
				    }
				    
				    HybEnsModel::StateDescription::ISite is2(k1,k2,l1,l2);

				    HybEnsModel::StateDescription state(i1,i2,j1,j2,k1,k2,l1,l2);

				    double energy_hyb2 = model_.energy_hybrid(is2);
				    
				    if (energy_hyb2 > max_hyb_energy_) continue;
				    
				    //double energy_unpair=model_.energy_unpair(is1,is2);
				    
				    //double total_energy = energy_hyb1 + energy_hyb2 + energy_unpair;
				    
				    double total_energy=model_.energy(state);
				    
				    if (total_energy <= max_total_energy_) {
					write_state(out,total_energy,state,binary_);
					count_double_states++;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
    
    return count_double_states;
}
