#ifndef BASIN_TRANSITION
#define BASIN_TRANSITION

#include <unordered_map>
#include <iostream>
#include <vector>
#include <stack>
#include "basin.hh"

/**
 * Partition function of a single basin transition
*/
using BasinTransition = double;

/**
 * Symmetric matrix of basin transitions partition functions
 */
class BasinTransitions {
private:

    using value_t = double;

    // Note on implementation:
    // The matrix of transition pfs is symmetric, so it would suffice to store values for
    // i<=j.
    // Moreover, it is sparse, representing only non-zero transitions pfs.
    // In the implementation, we need to traverse over all (non-zero transition) neighbors of a
    // given state. This can be implemented much more conveniently by keeping the matrix
    // symmetric.


    //! a row of a transition map "matrix"
    using transitions_map_row_t = std::unordered_map<size_t,BasinTransition>;

    //! map from basin indices to basin transition define as map of
    //! maps, such that we can iterate over first index (equivalently,
    //! second index, in the case of a symmetric matrix; see below!)
    using transitions_map_t = std::unordered_map<size_t,transitions_map_row_t>;
public:

    // read access
    const BasinTransition
    get(size_t i, size_t j) const {
        auto it = transitions_.find(i);
        if (it==transitions_.end()) return 0;
        auto it2 = it->second.find(j);
        if (it2==it->second.end()) return 0;
        return it2->second;
    };

    /** @brief Write access
    */
    void
    set(size_t i, size_t j, const BasinTransition &value) {
        transitions_[i][j] = value;
        transitions_[j][i] = transitions_[i][j];
    };

    /** @brief Add
    */
    void
    add(size_t i, size_t j, value_t value) {
        transitions_[i][j] += value;
        transitions_[j][i] = transitions_[i][j];
    };

    /** @brief Write access
    */
    void
    multiply(size_t i, size_t j, value_t factor) {
        transitions_[i][j] *= factor;
        transitions_[j][i] = transitions_[i][j];
    };

    /** @brief Transitions to neighbors
     *
     * @param i basin index
     * @returns transitions to neighbors as iterable map
     */
    const transitions_map_row_t &
    neighbors(size_t i) const {
        auto row = transitions_.find(i);
        if (row == transitions_.end()) return empty_row_;
        return row->second;
    }

    /** Erase all transitions from and to a basin
        @param i index of basin
     */
    void
    erase_basin(size_t i) {
        // 1) go through neighbors and delete their transition to i
        auto &trs = neighbors(i);
        for (auto it=trs.begin(); trs.end()!=it; ++it) {
    	    transitions_[it->first].erase(i);
        }
        // 2) delete all transitions from i
        transitions_.erase(i);
    }

    template<class Predicate>
    size_t
    filter(size_t i, Predicate filter_fun) {
        size_t removed_transitions = 0;
        auto &trs = transitions_[i];
        for (transitions_map_row_t::iterator it=trs.begin(); trs.end()!=it; ) {
    	    //assert(!basins_[i].merged());

    	    if ( ! filter_fun(it->first, it->second) ) {
    	        //remove transition
    	        removed_transitions++;

                transitions_[it->first].erase(i);
                it=trs.erase(it);
    	    } else {
    	        ++it;
    	    }
        }
        return removed_transitions;
    }

    void
    connected_components(const std::vector<Basin> &basins,
        std::vector<size_t> &components,
        std::vector<size_t> &component_sizes,
        std::vector<double> &component_pfs
        ) const
    {
        // components[i] is the number of the component of basin i
        // or 0 if the basin is still unassigned
        components.clear();
        components.assign(basins.size(),0);

        component_sizes.clear();
        component_pfs.clear();

        // current component number
        size_t c=0;

        typedef std::pair<size_t,transitions_map_row_t::const_iterator> stack_entry;
        std::stack<stack_entry> stack;

        for (size_t i=0; i<basins.size(); i++) {
            if (components[i] > 0) continue;
            c++;

            stack.push(stack_entry(i,transitions_.find(i)->second.begin()));
            components[i]=c;
            component_sizes.push_back(1);
            component_pfs.push_back(0);

            while (!stack.empty()) {
               stack_entry &top = stack.top();

               while (top.second != transitions_.find(top.first)->second.end()
           	   && components[top.second->first]>0) {
           	   top.second++;
                }
                if (top.second == transitions_.find(top.first)->second.end()) {
                    stack.pop();
                    continue;
                }

                size_t j=top.second->first;
                stack.push(stack_entry(j,transitions_.find(j)->second.begin()));

                components[j]=c;
                component_sizes[c-1]++;
                component_pfs[c-1] += basins[j].Z();
            }
        }
    }

    void
    swap(size_t x, size_t y) {
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

private:
    transitions_map_t transitions_;

    transitions_map_row_t empty_row_;
};

#endif
