#include <iostream>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <cmath>
#include <set>
#include <sstream>

#include <algorithm>

#include <cassert>

#include <zlib.h>

#include "basin_transition.hh"
#include "barrier_graph.hh"
#include "pair_pfs.hh"

BarrierGraph::BarrierGraph( bool special_first_state,
			    bool verbose,
			    bool debug_out
			    )
    : special_first_state_(special_first_state),
      verbose_(verbose),
      debug_out_(debug_out),
      track_pruning_(false),
      num_input_states_(0),
      min_contribution_(0)
{
}


void
BarrierGraph::push_back_basin(size_t index, double pf) {
    basins_.push_back(Basin(index,pf));
    if (track_pruning_) {
	basin_pruning_infos_.push_back(BasinPruningInfo(index));
    }
}


BarrierGraph::BarrierGraph(std::istream &in,
			   double min_rate,
			   bool special_first_state,
			   bool verbose,
			   bool debug_out,
			   bool track_pruning
			   )
    : special_first_state_(special_first_state),
      verbose_(verbose),
      debug_out_(debug_out),
      track_pruning_(track_pruning),
      num_input_states_(0),
      min_contribution_(1e-3)
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

	// pad with "merged" basins if index is not continuous in input
	while (basins_.size()<i) {
	    push_back_basin(basins_.size(),0);
	    basins_[basins_.size()-1].mark_merged();
	}
	// register input basin
	push_back_basin(i,pf);

	while( in ) {
	    size_t j; // column index
	    in.read(reinterpret_cast<char *>(&j),sizeof(j));
	    if (j==stopper) { break; }

	    double tpf;
	    in.read(reinterpret_cast<char *>(&tpf),sizeof(tpf));

	    num_rates++;
	    if (i!=j) { // ignore any i==j pfs in the input
                transitions_.set(i,j,tpf);
            }
	}
    }

    std::cerr << "Read "<<num_rates<<" rates."<<std::endl;

    num_input_states_ = basins_.size();
}

void
BarrierGraph::multiply_transitions_from_to(int basin_idx,
                                           double factor) {
    for(auto &x: transitions_.neighbors(basin_idx)) {
        transitions_.multiply(basin_idx, x.first, factor);
    }
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
            const auto neighbors = transitions_.neighbors(i);
     	    n += neighbors.size();
	}
    }
    return n;
}

void
BarrierGraph::filter_basin_transitions(const Basin &x0, double min_rate) {
    removed_transitions_ +=
        transitions_.filter(x0.idx(),
            [&](size_t j, double Ztr) { return Ztr / std::min(x0.Z(), basins_[j].Z()) >= min_rate; } );
}


double
BarrierGraph::outflow_pf(const Basin &x) const {
    double total_out=0;

    const auto &trs_x = transitions_.neighbors(x.idx());

    for (auto it=trs_x.begin(); trs_x.end()!=it; ++it) {
	if (basins_[it->first].merged()) continue;
	if (it->first==x.idx()) continue; // never include self transitions

	total_out += it->second;
    }
    return total_out;
}

double
BarrierGraph::max_outflow(const Basin &x) const {
    double max_outflow = 0.0;

    const auto &trs_x = transitions_.neighbors(x.idx());

    for (auto it=trs_x.begin(); trs_x.end()!=it; ++it) {
	if (basins_[it->first].merged()) continue;
	if (it->first==x.idx()) continue; // never include self transitions

	max_outflow = std::max(max_outflow, it->second);
    }
    return max_outflow/x.Z();
}

void
BarrierGraph::merge_in_basin(Basin &x0, Basin &y, double fraction) {
    y.merge_in(x0,fraction);

    if (track_pruning_) {
	size_t orig_index_x0 = basin_pruning_infos_[x0.idx()].orig_index();
	basin_pruning_infos_[y.idx()].merge_in_basin(orig_index_x0,
						     fraction,
						     basin_pruning_infos_[x0.idx()],
						     min_contribution_ / num_input_states_);
    }
}


void
BarrierGraph::dissolve_basin(Basin &x0, double min_rate) {

    transitions_.assert_consistency();

    // compute total outflow
    double total_out_x0 = outflow_pf(x0);

    double Z_x0x0 = transitions_.get(x0.idx(),x0.idx());

    // traverse neighbor basins of x0
    auto &trs_x0 = transitions_.neighbors(x0.idx());
    for (auto it=trs_x0.begin(); trs_x0.end()!=it; ++it) {

	Basin &y = basins_[it->first];

	// skip merged basins, and x0 itself
	if (y.merged()) continue;
	if (x0.idx()==y.idx()) continue;

	// 1) distribute the pf of x0 to its neighbors y
	double Z_yx0 =  it->second; // transition pf between x0 and y
	double fraction_yx0 = Z_yx0/total_out_x0; // fraction of the x0 outflow that flows to y

        if (x0.idx() != y.idx()) {
	    merge_in_basin(x0,y,fraction_yx0); // merge fraction of x0's pf into y's pf

	    if (debug_out_) {
	        std::cerr << "STATE Transfer " << fraction_yx0*100 << "% of " << x0.idx()
		          << "'s pf to " << y.idx()<<std::endl;
	    }
        }

	// 2) distribute the partition function of the
	//    transitions between x0 and its neighbors to transitions between neighbors of x0
	for (auto it2=trs_x0.begin(); trs_x0.end()!=it2; ++it2) {

	    Basin &x = basins_[it2->first];
	    if (x.merged()) continue;

	    if (x0.idx()==x.idx()) continue;

	    // make sure we see only one of the symmetric cases x,y and y,x
            // note that we cover x==y !
	    if (x.idx()<y.idx()) continue;

	    // update the transition from x to y (via x0) and vice versa

	    double Z_xx0 = transitions_.get(x.idx(), x0.idx());

	    double fraction_xx0 = Z_xx0/total_out_x0;


            double Z_add_x = Z_xx0 * fraction_yx0;
            double Z_add_y = Z_yx0 * fraction_xx0;
            double Z_add_x0 = Z_x0x0 * fraction_yx0 * fraction_xx0;

	    double Z_add = Z_add_x + Z_add_y + Z_add_x0;

	    if (debug_out_) {
		std::cerr << "Transfer " << fraction_yx0*100 << "% of "
			  << x.idx() <<"-"<< x0.idx() << " to " << y.idx()
			  << "-"<< x.idx()<<std::endl;
		std::cerr << "Transfer " << fraction_xx0*100 << "% of "
			  << y.idx() <<"-"<< x0.idx() << " to " << y.idx()
			  << "-"<< x.idx()<<std::endl;
		std::cerr << "Transfer " << fraction_yx0*fraction_xx0*100 << "% of "
			  << x0.idx() <<"-"<< x0.idx() << " to " << y.idx()
			  << "-"<< x.idx()<<std::endl;
	    }

	    // HEURISITC SPARSIFICATION; IMPORTANT FOR EFFICIENCY:
	    // update transitions pf only if it exceeds minrate in one of
            // the directions
	    //
            // Without this heuristic, the pruning procedure
	    // produces a lot of very small transitions.
	    if ( (transitions_.get(x.idx(), y.idx()) + Z_add) / std::min(x.Z(),y.Z()) < min_rate )  {
	    	if (debug_out_) {std::cerr << "Transition pf not updated due to sparsification"<<std::endl;}
                continue;
            }

            // Increase the transition pf between x and y due to dissolved
            // basin x0

            // add same term to both 'directions' to keep symmetry
            transitions_.add(x.idx(), y.idx(), Z_add);

	} // end for it2 (over neighbors of x0)

    } // end for it (over neighbors of x0)


    // finally, mark x0 as merged
    if (debug_out_) {
	std::cerr << "Mark "<<x0.idx()<<" as merged."<<std::endl;
    }
    x0.mark_merged();

    transitions_.erase_basin(x0.idx());
}

bool
BarrierGraph::is_to_be_merged(const Basin &x0,double max_outflow,double min_p_equ) const {
    double total_Z = compute_Z();

    // compute total outflow
    double total_out = outflow_pf(x0);

    // probability in equilibrium
    double p_equ = x0.Z() / total_Z;

    bool to_be_merged =
	(p_equ < min_p_equ)
	|| (total_out/x0.Z() > max_outflow)
	;

    if (debug_out_) {
	if (to_be_merged) {
	    std::cerr << "Select basin "<<x0.idx()<<" for merge because ";
	    if (total_out/x0.Z()) {
		std::cerr << "total outflow is "<<total_out/x0.Z();
	    }
	    if (p_equ < min_p_equ) {
		std::cerr << "equilibrium probability is "<<p_equ;
	    }
	    std::cerr<<std::endl;
	}
    }

    if ( special_first_state_ &&  x0.idx()==0 ) {
	if (verbose_) {
	    std::cerr << "Keep (special) first basin ( p="<< p_equ << ", outflow="<< total_out/x0.Z() <<" )."<<std::endl;
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

    removed_transitions_=0;

    // run through sorted basins and merge
    for (size_t i=0; i<basins_.size(); ++i) {

	Basin &x0 = basins_[sorted_basin_idxs[i]];

	//if (x0.merged()) continue;
	assert(!x0.merged());

	// remove transitions with rates lower than min_rate
	filter_basin_transitions(x0,min_rate);

        if (is_to_be_merged(x0,max_outflow,min_p_equ)) {

	    if (debug_out_) {
		std::cerr << "Dissolve basin "<< x0.idx() << " with outflow "<< outflow_pf(x0)/x0.Z() << std::endl;
	    }

	    dissolve_basin(x0, min_rate);

            if (track_pruning_) {
                // this saves a lot of memory!
                // remove pruning information of dissolved basin x0
                basin_pruning_infos_[x0.idx()].clear();
            }
	} else {
            if (track_pruning_) {
                //! @todo check the effect of this sparsification; seems to be small
                basin_pruning_infos_[x0.idx()].sparsify(min_contribution_);
            }
        }
    }

    std::cerr << "Removed "<<removed_transitions_<<" transitions due to min-rate filter;"
	      << "finally keeping "<<num_transitions()<<" tranisitions."
	      << std::endl;
}

// reduce to keep set
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

	// remove transitions with rates lower than min_rate
	filter_basin_transitions(x0,min_rate);

        // dissolve x0 unless it is in the to_keep set
        if (to_keep.count(x0.idx()) == 0) {
	    if (verbose_) {
		std::cerr << "Dissolve basin "<< x0.idx() << std::endl;
	    }

	    dissolve_basin(x0,min_rate);
	}
    }
}

void
BarrierGraph::print_edges(std::ostream &out,const Basin &b) const {
    out << b.idx()
	// << " (" << b.Z() << ") "
	<< "  -> ";

    auto ts = transitions_.neighbors(b.idx());

    for(auto it=ts.begin(); it != ts.end(); ++it) {
	if ((it->first!=b.idx()) && (!basins_[it->first].merged())) {
	    out << " " << it->first
		<< " " << (it->second / b.Z());
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
	    total_Z += basins_[i].Z();
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

	    double p_equ = basins_[i].Z() / total_Z;

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

bool
BarrierGraph::write_treekin_ratesmatrix_row(std::ostream &out, size_t i) const {
    if (!basins_[i].merged()) {

	double total=0.0; // row total rate (i!=j)

	std::vector<double> row(basins_.size());

    	auto &ts = transitions_.neighbors(i);
        // "bucket sort"
	for(auto it=ts.begin(); it != ts.end(); ++it) {
	    double rate = (it->second / basins_[i].Z());
	    row[it->first] = rate;
	    total += rate;
	}

	row[i] = -total;

	for(size_t j=0; j<row.size();++j) {
	    if (!basins_[j].merged()) {
		out << format_rate_for_treekin(row[j]) << " ";
	    }
	}
	out << "\n";
	return true;
    }
    return false;
}

void
BarrierGraph::write_treekin_ratesmatrix(std::ostream &out) const
{
    for(size_t i=0; i<basins_.size(); ++i) {
	write_treekin_ratesmatrix_row(out, i);
    }
}

gzFile
BarrierGraph::gzwrite_treekin_ratesmatrix(gzFile fh) const
{
    for(size_t i=0; i<basins_.size(); ++i) {
	std::stringstream out;
	if (write_treekin_ratesmatrix_row(out, i)) {
	    gzputs(fh,out.str().c_str());
	}
    }
    return fh;
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

	    auto &ts = transitions_.neighbors(i);
	    // "bucket sort"
	    for(auto it=ts.begin(); it != ts.end(); ++it) {
	        row[it->first] = it->second;
	    }

	    // set diagonal
	    row[i] = basins_[i].Z();

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

    // write transition pfs; for each i, write list of pairs j tpf_ij
    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {

	    out.write(reinterpret_cast<const char *>(&i),sizeof(i));
	    double pf = basins_[i].Z();
	    out.write(reinterpret_cast<const char *>(&pf),sizeof(pf));

	    const auto trs = transitions_.neighbors(i);
	    for(auto it=trs.begin(); it != trs.end(); ++it) {
		size_t j=it->first;
		double tpf = it->second;

		// symmetric min-rate criterion!
		if ( tpf / std::min(pf, basins_[j].Z()) >= min_rate ) {
		    // write index
		    out.write(reinterpret_cast<const char *>(&j),sizeof(j));
		    // write pf
		    out.write(reinterpret_cast<const char *>(&tpf),sizeof(tpf));
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
	    auto &ts = transitions_.neighbors(i);
	    for(auto it=ts.begin(); it != ts.end(); ++it) {
		size_t j=it->first;
		double rate = it->second/basins_[i].Z();

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
BarrierGraph::connect_components(size_t component_num,
                                 const std::vector<size_t> &components,
                                 double rate) {

    std::vector<size_t> max_idx;
    std::vector<double> max_Z;

    max_idx.assign(component_num,0);
    max_Z.assign(component_num,0.0);

    //get component maxima
    for (auto &basin: basins_) {
        size_t c=components[basin.idx()]-1;
        if ( basin.Z() > max_Z[c] ) {
            max_Z[c] = basin.Z();
            max_idx[c] = basin.idx();
        }
    }

    // connect all component maxima with rate
    for ( size_t x: max_idx ) {
        for ( size_t y: max_idx) {
            if (x<y) {
                double minZ = std::min(basins_[x].Z(),basins_[y].Z());
                double tZ = rate * minZ; // rate = tZ/Z, i.e. tZ=rate*Z
                transitions_.set(x, y, tZ);
            }
        }
    }
}



void
BarrierGraph::reindex_basin(Basin &x, size_t new_index) {
    size_t old_index=x.idx();

    x.set_idx(new_index);
    basins_[new_index]=x;

    if (track_pruning_) {

	// copy pruning information
	basin_pruning_infos_[new_index]=basin_pruning_infos_[old_index];

    }
}

void
BarrierGraph::reindex(const std::vector<bool> &keep) {

    std::vector<size_t> old2new(basins_.size());


    // phase 0: remove entries in the transitions matrix that correspond to merged basins

    //   for all basins not in keep: erase the corresponding transitions
    for (auto it=basins_.begin(); basins_.end()!=it; ++it) {
	if (!keep[it->idx()]) {
	    transitions_.erase_basin(it->idx());
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
		    reindex_basin(*it,i);
		}
		i++;
	    }
	}
	basins_.resize(i);
	if (track_pruning_) {
	    basin_pruning_infos_.resize(i);
	}
    }

    // phase 2; remap transitions
    BasinTransitions new_transitions;

    for(size_t i=0; i<basins_.size(); ++i) {
        if (!keep[i]) continue;
        auto &trs = transitions_.neighbors(i);
	for(auto it=trs.begin(); it != trs.end(); ++it) {
	    new_transitions.set( old2new[i], old2new[it->first], it->second );
	}
    }

    transitions_ = new_transitions;

}


void
BarrierGraph::connected_components(std::vector<size_t> &components,
                                   std::vector<size_t> &component_sizes,
                                   std::vector<double> &component_pfs
                                   ) const
{
    transitions_.connected_components(basins_,
        components, component_sizes, component_pfs);
}

void
BarrierGraph::swap_indices(size_t x, size_t y) {

    // first swap basins
    std::swap(basins_[x],basins_[y]);
    // (don't forget to swap indices of basins)
    size_t tmp = basins_[x].idx();
    basins_[x].set_idx(basins_[y].idx());
    basins_[y].set_idx(tmp);

    // swap pruning information
    if (track_pruning_) {
	// swap pruning information records x and y
	std::swap(basin_pruning_infos_[x],basin_pruning_infos_[y]);
    }




    // then swap transitions:
    transitions_.swap(x,y);
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
	       << (- RT * log(basins_[i].Z()));

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
	    double p_equ = basins_[i].Z() / total_Z;
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

std::ostream &
BarrierGraph::write_pruning_track(std::ostream &out, bool sparse) const {
    assert(track_pruning_);

    for (size_t i=0; i<basins_.size(); i++) {
	// out << i
	//     << "\t";
	basin_pruning_infos_[i].write(out, min_contribution_, num_input_states_, sparse);
	out << '\n';
    }

    return out;
}

gzFile
BarrierGraph::gzwrite_pruning_track(gzFile fh, bool sparse) const {
    assert(track_pruning_);

    for (size_t i=0; i<basins_.size(); i++) {
	std::stringstream out;
	basin_pruning_infos_[i].write(out, min_contribution_, num_input_states_, sparse);
	out << '\n';
	gzputs(fh,out.str().c_str());
    }

    return fh;
}


bool
BarrierGraph::compute_pruning_pps(size_t ms_idx,
                                  const PairPfs &ppfs,
                                  double &Z,
                                  LocARNA::SparseMatrix<double> &ms_ppfs) const {
    assert(ms_idx<basins_.size());

    // current macrostate
    const Basin &ms = basins_[ms_idx];
    if (ms.merged()) {return false;} // skip merged macrostates

    Z=ms.Z();

    //std::cerr << "BarrierGraph::compute_pruning_pps "<<ms_idx<<" "<<Z<<std::endl;


    // pruning info of macrostate
    const BasinPruningInfo &pi = basin_pruning_infos_[ms_idx];

    // iterate over basins that contribute to this macrostate
    for ( auto &x : pi ) {
	size_t basin_idx  = x.first;  // (original) index of contributing basin
	double basin_frac = x.second; // fraction

	double basin_pf = ppfs.Z(basin_idx); //total partition function of basin

	//std::cerr << "  contributions "<<basin_idx<<" "<<basin_frac<<" "<<basin_pf<<std::endl;

	// iterate over (k1,k2) pairs (that  exist in the sparse
	// representation of the current basin)

	auto &basin_pps = ppfs.pair_probs(basin_idx);

	for ( auto &y : basin_pps ) {
	    size_t k1=y.first.first;
	    size_t k2=y.first.second;
	    double pr=y.second; // Pr[(k1,k2)|basin]

	    // std::cerr << "    pair "<<k1<<" "<<k2<<" "<<pr<<" "<<(basin_frac * basin_pf * pr)<<"\t";

	    ms_ppfs(k1,k2) = ms_ppfs(k1,k2) +  basin_frac * basin_pf * pr;
	}
	//std::cerr << std::endl;
    }
    return true;
}

std::ostream &
BarrierGraph::write_pruning_pps(std::ostream &out,
				const PairPfs &ppfs,
				double min_prob) const {
    assert(track_pruning_);

    //iterate over macrostates
    for (size_t ms_idx = 0; ms_idx<basins_.size(); ms_idx++) {
	// pair partition functions of ms
	LocARNA::SparseMatrix<double> ms_ppfs;
	// total partition function of ms
	double Z;

	if ( !compute_pruning_pps(ms_idx,ppfs,Z,ms_ppfs) ) {continue;}

	// now, ppfs contains the pair partition functions for the current macrostate
	// ==> write them
	PairPfs::write_state(out,ms_idx,Z,ms_ppfs,min_prob);
    }

    return out;
}

gzFile
BarrierGraph::gzwrite_pruning_pps(gzFile fh,
				  const PairPfs &ppfs,
				  double min_prob) const {
    // variant of write_pruning_pps_track() for gz-compressed output

    assert(track_pruning_);

    for (size_t ms_idx = 0; ms_idx<basins_.size(); ms_idx++) {
	LocARNA::SparseMatrix<double> ms_ppfs;
	double Z;
	if ( !compute_pruning_pps(ms_idx,ppfs,Z,ms_ppfs) ) {continue;}

	std::ostringstream out;
	PairPfs::write_state(out,ms_idx,Z,ms_ppfs,min_prob);

	gzputs(fh,out.str().c_str());
    }

    return fh;
}
