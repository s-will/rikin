#include "basin_transition.hh"
#include "barrier_graph.hh"

#include <stack>
#include <iostream>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <cmath>
#include <set>
#include <sstream>

#include <algorithm>

#include <cassert>


BarrierGraph::BarrierGraph(std::istream &in,
			   double min_rate,
			   bool special_first_state,
			   bool verbose,
			   bool debug_out
			   )
    : special_first_state_(special_first_state),
      verbose_(verbose),
      debug_out_(debug_out)
{
    size_t stopper = std::numeric_limits<size_t>::max();
    size_t num_rates=0;
    
    // // get matrix dimension / number of basins
    // size_t dim;
    // in.read(reinterpret_cast<char *>(&dim),sizeof(dim));
    
    while( in ) {
	size_t i=0; // row index
	in.read(reinterpret_cast<char *>(&i),sizeof(i));
	if (i==stopper) { break; }

	double pf;
	in.read(reinterpret_cast<char *>(&pf),sizeof(pf));
	while (basins_.size()<i) {
	    Basin b = Basin(basins_.size(),0);
	    b.mark_merged();
	    basins_.push_back(b);
	}
	basins_.push_back(Basin(i,pf));
	
	while( in ) {
	    size_t j; // column index
	    in.read(reinterpret_cast<char *>(&j),sizeof(j));
	    if (j==stopper) { break; }
	        
	    double tpf;
	    in.read(reinterpret_cast<char *>(&tpf),sizeof(tpf));
	    
	    num_rates++;
	    transitions_[i][j]=BasinTransition(tpf);
	}
    }
    std::cerr << "Read "<<num_rates<<" rates."<<std::endl;
}


size_t
BarrierGraph::num_basins() const {
    size_t n=0;
    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    n++;
	}
    }
    return n;
}


size_t
BarrierGraph::num_transitions() const {
    size_t n=0;
    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    const transitions_map_t::const_iterator &it = transitions_.find(i);
	    if(it != transitions_.end()) {
		n += it->second.size();
	    }
	}
    }
    return n;
}

void
BarrierGraph::filter_basin_transitions(const Basin &x0, double min_rate) {
    transitions_map_row_t &trs_x0=transitions_[x0.idx()];
    for (transitions_map_row_t::iterator it=trs_x0.begin();
	 trs_x0.end()!=it; ) {
	
	//if (basins_[it->first].merged()) continue;
	assert(!basins_[it->first].merged());
	    
	if (it->second.get_Z() / std::min(x0.get_Z(), basins_[it->first].get_Z())
	    < min_rate) {
	    
	    //remove transition
	    //std::cerr << "Remove transition between "<<x0.idx()<<" and "<<it->first<<std::endl; 
	    transitions_[it->first].erase(x0.idx());
	    it=trs_x0.erase(it);
	    
	} else {
	    ++it;
	}
    }
    //std::cerr << "DONE"<<std::endl;
}

void
BarrierGraph::filter_transitions(double min_rate) {
    // run through basins
    for (size_t i=0; i<basins_.size(); ++i) {
	
	Basin &x0 = basins_[i];
	
	filter_basin_transitions(x0,min_rate);
    }    
}


void
BarrierGraph::dissolve_basin(Basin &x0) {

    // compute total outflow
    double total_out = outflow_pf(x0);
    
    transitions_map_row_t &trs_x0=transitions_[x0.idx()];
    for (transitions_map_row_t::iterator it=trs_x0.begin();
	 trs_x0.end()!=it; ++it) {
	
	Basin &y = basins_[it->first];
	if (y.merged()) continue;
	
	if (x0.idx()==y.idx()) continue;
		
	// 1) distribute the basin's partition function to its neighbors
	double Z_yx0 =  it->second.get_Z();
	double fraction = Z_yx0/total_out;
		
	y.merge_in(x0,fraction);
		
	if (debug_out_) {
	    std::cerr << "Transfer a fraction of " << fraction << " of " << x0.idx()
		      << "'s partfunc to " << y.idx()<<std::endl;
	}
		
	// 2) distribute the partition function of the transitions to this basin to its neighbors
		
	for (transitions_map_row_t::iterator it2=trs_x0.begin();
	     trs_x0.end()!=it2; ++it2) {
		    
	    Basin &x = basins_[it2->first];
	    if (x.merged()) continue;
				    
	    if (x.idx()==y.idx()) continue;
	    if (x0.idx()==x.idx()) continue;
		    
	    // update the transition from x to y (via x0)
		    
	    double Z_xx0 = transitions_[x.idx()][x0.idx()].get_Z();
		    
	    if (debug_out_) {
		std::cerr << "Transfer a fraction of " << fraction << " of transition "
			  << x.idx() <<"->"<< x0.idx() << " to " << x.idx() 
			  << "->"<< y.idx()<<std::endl;
	    }
		    
	    transitions_[x.idx()][y.idx()].update(Z_xx0 * fraction);		    
	} // end for it2 (over neighbors of x0)
		
    } // end for it (over neighbors of x0)
	    
    // finally, mark as merged
    if (debug_out_) {
	std::cerr << "Mark "<<x0.idx()<<" as merged."<<std::endl;
    }
    x0.mark_merged();
	    
    // and release all transitions from and to the merged basin
    for (transitions_map_row_t::iterator it=trs_x0.begin();
	 trs_x0.end()!=it; ++it) {
	transitions_[it->first].erase(x0.idx());
    }
    transitions_.erase(x0.idx());
    
}


bool
BarrierGraph::is_to_be_merged(const Basin &x0,double max_outflow,double min_p_equ) const {
    double total_Z = compute_Z();
    
    // compute total outflow
    double total_out = outflow_pf(x0);
    
    // probability in equilibrium
    double p_equ = x0.get_Z() / total_Z;
    
    bool to_be_merged = 
	(p_equ < min_p_equ) 
	|| (total_out/x0.get_Z() > max_outflow)
	;
    
    if (debug_out_) {
	if (to_be_merged) {
	    std::cerr << "Select basin "<<x0.idx()<<" for merge because ";
	    if (total_out/x0.get_Z()) {
		std::cerr << "total outflow is "<<total_out/x0.get_Z();
	    }
	    if (p_equ < min_p_equ) {
		std::cerr << "equilibrium probability is "<<p_equ;
	    }
	    std::cerr<<std::endl;
	}
    }
    
    if ( special_first_state_ &&  x0.idx()==0 ) {
	if (verbose_) {
	    std::cerr << "Keep (special) first basin ( p="<< p_equ << ", outflow="<< total_out/x0.get_Z() <<" )."<<std::endl;
	}
	to_be_merged=false;
    }
    
    return to_be_merged;
}

void
BarrierGraph::prune(double max_outflow,double min_p_equ, double min_rate) {
    // sort basins increasing by their partition function
    std::vector<size_t> sorted_basin_idxs;

    for (size_t i=0; i<basins_.size(); ++i) sorted_basin_idxs.push_back(i);
    sort(sorted_basin_idxs.begin(),sorted_basin_idxs.end(),compBasinIdxs(*this));
        
    // run through sorted basins and merge
    for (size_t i=0; i<basins_.size(); ++i) {
	
	Basin &x0 = basins_[sorted_basin_idxs[i]];
	
	//if (x0.merged()) continue;
	assert(!x0.merged());
	
	
	transitions_map_row_t &trs_x0=transitions_[x0.idx()];
	
	// remove transitions from state x0 to itself
	if (trs_x0.end()!=trs_x0.find(x0.idx())) {
	    trs_x0.erase(x0.idx());
	}
	
	// remove transitions wiht rates lower than min_rate
	filter_basin_transitions(x0,min_rate);
	            
        if (is_to_be_merged(x0,max_outflow,min_p_equ)) {
	    
	    if (debug_out_) {
		std::cerr << "Dissolve basin "<< x0.idx() << " with outflow "<< outflow_pf(x0)/x0.get_Z() << std::endl;
	    }
	    
	    dissolve_basin(x0);
	
	}
    }
}


void
BarrierGraph::reduce_basin_set(const std::set<size_t> &to_keep, double min_rate) {
    // sort basins increasing by their partition function
    std::vector<size_t> sorted_basin_idxs;

    for (size_t i=0; i<basins_.size(); ++i) sorted_basin_idxs.push_back(i);
    sort(sorted_basin_idxs.begin(),sorted_basin_idxs.end(),compBasinIdxs(*this));
        
    // run through sorted basins and merge
    for (size_t i=0; i<basins_.size(); ++i) {
	
	Basin &x0 = basins_[sorted_basin_idxs[i]];
	
	//if (x0.merged()) continue;
	assert(!x0.merged());
	
	
	transitions_map_row_t &trs_x0=transitions_[x0.idx()];
	
	// remove transitions from state x0 to itself
	if (trs_x0.end()!=trs_x0.find(x0.idx())) {
	    trs_x0.erase(x0.idx());
	}
	
	// remove transitions wiht rates lower than min_rate
	filter_basin_transitions(x0,min_rate);
	
        if (to_keep.count(x0.idx()) == 0) {
	    
	    if (verbose_) {
		std::cerr << "Dissolve basin "<< x0.idx() << std::endl;
	    }
	    
	    dissolve_basin(x0);
	
	}
    }
}

void
BarrierGraph::print_edges(std::ostream &out,const Basin &b) const {
    out << b.idx()
	// << " (" << b.get_Z() << ") "
	<< "  -> ";
	
    transitions_map_t::const_iterator
	ts = transitions_.find(b.idx());
	
    if (ts == transitions_.end()) return;
	
    for(transitions_map_row_t::const_iterator it=ts->second.begin();
	it != ts->second.end();
	++it) {
	    
	if ((it->first!=b.idx()) && (!basins_[it->first].merged())) {
	    out << " " << it->first
		// << " (" << it->second.get_Z() << ")"
		<< " " << (it->second.get_Z() / b.get_Z());
	}
    }
    out << std::endl;
}

void 
BarrierGraph::print_barrier_graph(std::ostream &out) const {
    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    print_edges(out,basins_[i]);
	}
    }	
}

double BarrierGraph::compute_Z() const {
    double total_Z=0.0;
    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    total_Z += basins_[i].get_Z();
	}
    }
    return total_Z;
}


void
BarrierGraph::print_basins(std::ostream &out) const {
    basins_[0].print_header(out);
    out<< " \tmax_out \ttotal_out \tp_equ";
    out<<std::endl;

    double total_Z = compute_Z();

    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    basins_[i].print(out);
	    
	    double max_out=max_outflow(basins_[i]);
	    
	    double p_equ = basins_[i].get_Z() / total_Z;
	    
	    out << " \t" << max_out <<" \t" <<outflow_rate(basins_[i])
		<<" \t" << p_equ;
	    
	    out<<std::endl;
	} else {
	}
    }
}

std::string
BarrierGraph::format_rate_for_treekin(double x)  {
    // using sprintf here is somwehat ugly, but there seem to be no
    // easy workarounds (i.e. without boost)
    const size_t bufsiz=256;
    char buf[bufsiz];
    snprintf(buf,bufsiz,"%10.4g",x);
    buf[bufsiz-1]=0;
    return std::string(buf);
}

void
BarrierGraph::print_treekin_ratesmatrix(std::ostream &out) const
{
    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    
	    double total=0.0; // row total rate (i!=j)
	    
	    std::vector<double> row(basins_.size());
	    
	    const transitions_map_t::const_iterator &ts = transitions_.find(i);
	    if (ts != transitions_.end()) {
		// "bucket sort"
		for(transitions_map_row_t::const_iterator it=ts->second.begin();
		    it != ts->second.end(); ++it) {
		    double rate = (it->second.get_Z() / basins_[i].get_Z());
		    row[it->first] = rate;
		    total += rate;
		}
	    }
	    
	    row[i] = -total;

	    for(size_t j=0; j<row.size();++j) {
		if (!basins_[j].merged()) {
		    out << format_rate_for_treekin(row[j]) << " ";
		}
	    }
	    out << "\n";
	}
    }
}


void
BarrierGraph::print_pfs(std::ostream &out,bool binary) const {
    if (binary) {
	// count unmerged basins => dim 
	size_t dim=0;
	for(size_t i=0; i<basins_.size(); ++i) {
	    if (!basins_[i].merged()) {
		dim++;
	    }
	}
	
	//write maxtrix dimension dim
	out.write(reinterpret_cast<const char *>(&dim),sizeof(dim));
    }
    
    
    // write matrix of transition pfs, where the diagonal consists of basin pfs
    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    
	    std::vector<double> row(basins_.size());
	    
	    const transitions_map_t::const_iterator &ts = transitions_.find(i);
	    if (ts != transitions_.end()) {
		// "bucket sort"
		for(transitions_map_row_t::const_iterator it=ts->second.begin();
		    it != ts->second.end(); ++it) {
		    row[it->first] = it->second.get_Z();
		}
	    }
	    
	    // set diagonal
	    row[i] = basins_[i].get_Z();
	    
	    for(size_t j=0; j<row.size();++j) {
		if (!basins_[j].merged()) {
		    if (binary) {
			out.write(reinterpret_cast<const char *>(&row[j]),sizeof(row[j]));
		    } else {
			out << row[j] << " ";
		    }
		}
	    }
	    if (!binary) {out << "\n";}
	}
    }
}

std::ostream &
BarrierGraph::write_binary(std::ostream &out,  double min_rate) const {
    size_t stopper = std::numeric_limits<size_t>::max();

    // size_t dim=basins_.size();
    // out.write(reinterpret_cast<const char *>(&dim),sizeof(dim));
    
    // write matrix of transition pfs, where the diagonal consists of basin pfs
    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    
	    out.write(reinterpret_cast<const char *>(&i),sizeof(i));
	    double pf = basins_[i].get_Z();
	    out.write(reinterpret_cast<const char *>(&pf),sizeof(pf));
	    
	    const transitions_map_t::const_iterator &ts = transitions_.find(i);
	    if (ts != transitions_.end()) {
		for(transitions_map_row_t::const_iterator it=ts->second.begin();
		    it != ts->second.end(); ++it) {
		    
		    size_t j=it->first;
		    double tpf = it->second.get_Z();
		    
		    // symmetric min-rate criterion!
		    if ( tpf / std::min(pf, basins_[j].get_Z()) >= min_rate ) {
			// write index
			out.write(reinterpret_cast<const char *>(&j),sizeof(j));
			// write pf
			out.write(reinterpret_cast<const char *>(&tpf),sizeof(tpf));
		    }
		}
	    }
	    out.write(reinterpret_cast<const char *>(&stopper),sizeof(stopper));
	}
    }
    out.write(reinterpret_cast<const char *>(&stopper),sizeof(stopper));
    
    return out;
}

void
BarrierGraph::print_rxns(std::ostream &out,
			 const std::string &nameA,
			 const std::string &nameB,
			 const std::string &nameAB) const {
    const size_t openstate_index=1;
   
    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    const transitions_map_t::const_iterator &ts = transitions_.find(i);
	    if (ts != transitions_.end()) {
		for(transitions_map_row_t::const_iterator it=ts->second.begin();
		    it != ts->second.end(); ++it) {
		    size_t j=it->first;
		    double rate = it->second.get_Z()/basins_[i].get_Z();
		    
		    std::ostringstream subs1;
		    std::ostringstream subs2;
		    std::ostringstream prods1;
		    std::ostringstream prods2;
		    
		    if (i!=openstate_index) {
			subs1<<nameAB<<i+1;
			subs2<<nameAB<<i+1;
		    } else {
			subs1<<nameA<<"+"<<nameB;
			subs2<<nameA<<" "<<nameB;
		    }
		    
		    if (j!=openstate_index) {
			prods1<<nameAB<<j+1;
			prods2<<nameAB<<j+1;
		    } else {
			prods1<<nameA<<"+"<<nameB;
			prods2<<nameA<<" "<<nameB;
		    }

		    out << "Rxn " << subs1.str() << "->" << prods1.str() <<std::endl
			<< "Subs " << subs2.str() << std::endl
			<< "Prods " << prods2.str() <<std::endl
			<< "Rate " << rate <<std::endl;    
		    
		    out << std::endl;
		}
	    }
	}
    }
}

void
BarrierGraph::print_spcs(std::ostream &out,
			 const std::string &nameA,
			 const std::string &nameB,
			 const std::string &nameAB) const {
    const size_t openstate_index=1;

    if (nameA!=nameB) {
	out << nameA << " " << 1 << std::endl;
	out << nameB << " " << 1 << std::endl;
    }

    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged() && (i!=openstate_index)) {
	    out << nameAB << (i+1) << " " << 0<< std::endl;
	}
    }
}

void
BarrierGraph::reindex() {
    std::vector<bool> keep(basins_.size());
    for (std::vector<Basin>::iterator it=basins_.begin();
	 basins_.end()!=it; ++it) {
	keep[it->idx()] = !it->merged();
    }
    reindex(keep);
}

void
BarrierGraph::keep_single_component(size_t c,
				    const std::vector<size_t> &components) {
    std::vector<bool> keep(basins_.size());
    for (std::vector<Basin>::iterator it=basins_.begin();
	 basins_.end()!=it; ++it) {
	keep[it->idx()] = components[it->idx()]==c;
    }
    reindex(keep);
}


void
BarrierGraph::reindex(const std::vector<bool> &keep) {
    
    std::vector<size_t> old2new(basins_.size());
    
    
    // phase 0: remove transitions of merged basins_
    for (std::vector<Basin>::iterator it=basins_.begin();
	 basins_.end()!=it; ++it) {
	if (!keep[it->idx()]) {
	    transitions_.erase(it->idx());
	}
    }

    for(transitions_map_t::iterator it=transitions_.begin();
	it != transitions_.end(); ++it) {
	for (std::vector<Basin>::iterator it2=basins_.begin();
	     basins_.end()!=it2; ++it2) {
	    if (!keep[it2->idx()]) {
		it->second.erase(it2->idx());
	    }
	}
    }

    // phase 1: reindex basins_; copy in place
    {
	size_t i=0;
	for (std::vector<Basin>::iterator it=basins_.begin();
	     basins_.end()!=it; ++it) {
	    if (keep[it->idx()]) {
		old2new[it->idx()]=i;
		if (it->idx()!=i) {
		    it->set_idx(i);
		    basins_[i]=*it;
		}
		i++;
	    }
	}
	basins_.resize(i);
    }
    
    // phase 2; remap transitions
    transitions_map_t new_transitions;

    for(transitions_map_t::const_iterator it=transitions_.begin();
	it != transitions_.end(); ++it) {
	for(transitions_map_row_t::const_iterator it2=it->second.begin();
	    it2 != it->second.end(); ++it2) {
	    new_transitions[old2new[it->first]][old2new[it2->first]] = it2->second;
	}
    }

    
    transitions_.clear();
    transitions_ = new_transitions;

    for (size_t i=0; i<basins_.size(); i++) {
	if (transitions_.find(i) == transitions_.end()) {
	    //if (verbose_) std::cerr << "WARNING: no transitions for state "<<i<<std::endl;
	    transitions_[i] = transitions_map_row_t();
	}
    }
}

double
BarrierGraph::check_rates() const
{
    // check whether transition rate matrix is symmetric 
    // (i.e. whether process is in detailed balance)

    double max_diff=0.0;

    for(size_t i=0; i<basins_.size(); ++i) {
     	std::vector<double> row(basins_.size());
	
	const transitions_map_t::const_iterator &ts = transitions_.find(i);
	if (ts != transitions_.end()) {
	    for(transitions_map_row_t::const_iterator it=ts->second.begin();
		it != ts->second.end(); ++it) {
		size_t j=it->first;
		
		//warning: the following finds can fail if data structures are not symmetric
		double diff = abs(it->second.get_Z() -
				  transitions_.find(j)->second.find(i)->second.get_Z());
		
		diff /= it->second.get_Z();
		
		max_diff = std::max(max_diff,diff);
	    }
	}
    }
    
    return max_diff;
}

std::vector<size_t>
BarrierGraph::connected_components(std::vector<size_t> &components) const
{
    // components[i] is the number of the component of basin i
    // or 0 if the basin is still unassigned
    components.clear();
    components.assign(basins_.size(),0);
    
    std::vector<size_t> component_size;
    
    // current component number
    size_t c=0;
    
    typedef std::pair<size_t,transitions_map_row_t::const_iterator> stack_entry;
    std::stack<stack_entry> stack;
    
    for (size_t i=0; i<basins_.size(); i++) {
	
	if (components[i] > 0) continue;
	c++;

	stack.push(stack_entry(i,transitions_.find(i)->second.begin()));
	components[i]=c;
	component_size.push_back(1);

	while (!stack.empty()) {
	    stack_entry &top = stack.top();
	    
	    while (top.second != transitions_.find(top.first)->second.end()
		   &&
		   components[top.second->first]>0) 
		top.second++;
	    
	    if (top.second == transitions_.find(top.first)->second.end()) {
		stack.pop();
		continue;
	    }
	    
	    size_t j=top.second->first;
	    stack.push(stack_entry(j,transitions_.find(j)->second.begin()));
	    components[j]=c;
	    component_size[c-1]++;
	}
    }
    

    return component_size;
}

void
BarrierGraph::swap_indices(size_t x, size_t y) {
	
    // first swap basins
    std::swap(basins_[x],basins_[y]);
    // (don't forget to swap indices of basins)
    size_t tmp = basins_[x].idx();
    basins_[x].set_idx(basins_[y].idx());
    basins_[y].set_idx(tmp);
    
    // then swap transitions:
    // -- first swap the rows
    std::swap(transitions_[x],transitions_[y]);
    
    // -- then columns
    for (transitions_map_t::iterator it=transitions_.begin(); transitions_.end()!=it; ++it) {
	transitions_map_row_t &row = it->second;
	
	transitions_map_row_t::iterator itx = row.find(x);
	transitions_map_row_t::iterator ity = row.find(y);
	
	if (itx!=row.end() && ity!=row.end()) {
	    std::swap(row[x],row[y]);
	} else if (itx==row.end() && ity!=row.end()) {
	    row[x] = ity->second;
	    row.erase(y);
	} else if (itx!=row.end() && ity==row.end()) {
	    row[y] = itx->second;
	    row.erase(x);
	}
    }
}

void
BarrierGraph::print_treekin_barriers(std::ostream &out, 
				     std::string header,
				     double RT) const {
    //out << "     " << model_.seqA() << "&" << model_.seqB() << std::endl;
    out << "     " << header << std::endl;
    
    size_t count=0;
    for(size_t i=0; i<basins_.size(); ++i) {
	if (basins_[i].merged()) continue;
	count++;
	size_t ow = out.width(4);

	std::ostringstream energy;
	energy.setf(std::ios_base::fixed, std::ios_base::floatfield);
	energy << std::setprecision(2)
	       << std::setw(6)
	       << (- RT * log(basins_[i].get_Z()));

	out << count << std::setw(ow) 
	    << "-" //<< " " << to_dotbracket(basins_[i].get_local_minimum())
	    << " " << energy.str() << std::endl;
    }
}


// determine outflow distribution
void
BarrierGraph::outflow_distribution(std::vector<double> &outflows) const {
    outflows.clear();
    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    outflows.push_back( outflow_rate(basins_[i]) );
	}
    }
    std::sort(outflows.begin(),outflows.end());
}

// determine equilibrium probability distribution
void
BarrierGraph::pequ_distribution(std::vector<double> &pequs) const {
    double total_Z = compute_Z();
    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    double p_equ = basins_[i].get_Z() / total_Z;
	    pequs.push_back( p_equ );
	}
    }
    std::sort(pequs.begin(),pequs.end());
}

double
BarrierGraph::max_outflow_by_quantile(double q) const {
    assert(0<=q && q<=100);
    assert(num_basins()>0);
    std::vector<double> outflows;
    outflow_distribution(outflows);
    return outflows[q*(num_basins()-1)/100];
}

double
BarrierGraph::max_outflow_by_number(size_t n) const {
    std::vector<double> outflows;
    outflow_distribution(outflows);
    if (n==outflows.size()) return outflows[outflows.size()-1];
    return outflows[n-1];
}

double
BarrierGraph::min_pequ_by_quantile(double q) const {
    assert(0<=q && q<=100);
    std::vector<double> pequs;
    pequ_distribution(pequs);
    if (q==0) return 1.0;
    return pequs[pequs.size() - q*pequs.size()/100];
}

double
BarrierGraph::min_pequ_by_number(size_t n) const {
    std::vector<double> pequs;
    pequ_distribution(pequs);
    if (n==0) return 1.0;
    return pequs[pequs.size() - n];
}


std::ostream &
BarrierGraph::print_stats(std::ostream &out) const {
    out << "basins:            "<<num_basins()<<std::endl;
    out << "mean transitions:  "<<num_transitions()/num_basins()<<std::endl;
    
    size_t nbasins=num_basins();
    
    // double total_outflow=0.0;
    // for(size_t i=0; i<basins_.size(); ++i) {
    // 	if (!basins_[i].merged()) {
    // 	    total_outflow += outflow_rate(basins_[i]);
    // 	}
    // }
    // out << "mean outflow (rate):      "<<total_outflow/nbasins<<std::endl;


    // determine outflow distribution
    std::vector<double> outflows;
    outflow_distribution(outflows);

    // determine p_equ distribution
    std::vector<double> pequs;
    pequ_distribution(pequs);
    
    out << "#basins\t";
    for(size_t n=nbasins; n>1; n/=2) {
	out << n << "\t";
    }
    out << "\n";
    
    out << "max-out\t";
    for(size_t n=nbasins; n>1; n/=2) {
	std::ostringstream fmtd;
	fmtd.precision(2);
	fmtd << outflows[n-1];
	out << fmtd.str() << "\t";
    }
    out << std::endl;

    out << "min-p\t";
    for(size_t n=nbasins; n>1; n/=2) {
	std::ostringstream fmtd;
	fmtd.precision(2);
	fmtd << pequs[nbasins-n];
	out << fmtd.str() << "\t";
    }
    out << std::endl;
    

    return out;
}

