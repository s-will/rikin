import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import gzip as gz
import re
import scipy
import RNA
import os
import copy

suffix = [".svg", ".pdf"]
output_directory = "."

seqfontname = "FreeMono"
#import matplotlib.font_manager
#matplotlib.font_manager.findSystemFonts(fontpaths=None)

# ## Plotting code

def boltzmann_weight(e, temperature=37):
    R = RNA.GASCONST/1000
    T = temperature + 273.15
    return np.exp( - e /R/T )
def inv_boltzmann_weight(Z, temperature=37):
    R = RNA.GASCONST/1000
    T = temperature + 273.15
    return -R*T*np.log( Z )

        


# -

# ### A color map for probabilities
#
# Suitable color maps can be generated algorithmically e.g. with cubehelix; this provides much flexibility, e.g. to generate a different map for intra-molecular bp...  

# +
# Define a color palette (same colors as in the R scripts)
# to be used as colors of states 

state_colors = [(0,0,0.4,1), (0.8,0,0,1)]
for cname in ["Dark2", "Set2", "Pastel2"]:
    state_colors.extend( mpl.colormaps[cname].colors )

state_palette = sns.color_palette(state_colors)

# set as seaborn default palette
sns.set_palette(state_palette)


#cmap = sns.color_palette("YlOrBr", as_cmap=True)
#cmap.set_gamma(0.5)
#cmap.set_extremes(under='white')

#inter_probability_cmap = sns.cubehelix_palette(start=.5, rot=-.75, dark=0, light=1.0, as_cmap=True)
#inter_probability_cmap = sns.cubehelix_palette(start=.5, rot=-.5, dark=0, light=1.0, as_cmap=True)
inter_probability_cmap = sns.cubehelix_palette(start=1.5, rot=-0.75, dark=0.15, light=1.0, hue=1, gamma=1.5, as_cmap=True)
#inter_probability_cmap = sns.cubehelix_palette(start=0, rot=1, dark=0.15, light=1.0, as_cmap=True)
inter_probability_cmap = sns.color_palette("Blues", as_cmap=True)

intra_probability_cmap = sns.cubehelix_palette(start=-0.25, rot=-0.2, dark=0.15, light=1.0, hue=1, gamma=1.5, as_cmap=True)
intra_probability_cmap = sns.color_palette("Greens", as_cmap=True)

dotplot_probability_cmap = sns.cubehelix_palette(
    start=0.3,      # tuned toward magenta on cubehelix's hue wheel
    rot=0,          # no hue sweep - stays a fixed hue
    dark=0.25,      # dark end (not pure black)
    light=1.0,      # pure white at the low end
    hue=1,          # max saturation
    gamma=1.0,
    as_cmap=True
)
dotplot_probability_cmap = sns.color_palette("Greens", as_cmap=True)

# ### Plotting of state probabilities over time 

# +
def max_probabilities_over_time(state_probabilities, maxp_threshold=None):
    """
    Maximal probabilities over time
    
    Args:
        probabilities of states over time
        maxp_threshold filter by probabilty threshold
    Returns:
        dictionary of maximum state probabilities
    """
    
    return { state:max(state_probabilities[state]) for state in state_probabilities.columns
           if maxp_threshold and max(state_probabilities[state])>=maxp_threshold }
    
def prominent_states(state_probabilities, maxp_threshold=0.02):
    # filter states by their maximum probabilty maxp_threshold at any time
    return max_probabilities_over_time(state_probabilities, maxp_threshold).keys()


# -

def _plot_state_probability_kinetics(state_probabilities, *, shown_states=None, ax=None):
    """
    Classical kinetics plot of state probabilities over time
    """    
    
    important_states = shown_states if shown_states!=None else state_probabilities.columns()
    
    ax = sns.lineplot(state_probabilities[important_states],
                      #palette = mypalette,
                      dashes=False,
                      ax=ax)
    
    ax.set_ylim(-0.02,1.02)
    ax.set_xscale("log")
    ax.set_xlabel("Time")
    ax.set_ylabel("State Probability")
    ax.grid(alpha=0.5)

    return ax


# +
def _seqpos_ticks(length, units_per_tick = 10, addone=True):
    while length < units_per_tick:
        units_per_tick = int(units_per_tick/2)
    units_per_xtick=max(1,units_per_tick)
    posticks = list(range(units_per_tick,length+1,units_per_tick))
    if addone:
        posticks = [1] + posticks
    return posticks
    
def _plot_size(ax):    
    # dimensions of plot area in inch
    figsize = ax.get_figure().get_size_inches()
    bounds = ax.get_position().bounds
    return bounds[2]*figsize[0], bounds[3]*figsize[1]

def _unit_csize(ax, fontsize, xdim, ydim):
    width, height = _plot_size(ax)
    return (fontsize/72) / width * xdim, (fontsize/72) / height * ydim


# -

# ### Plotting states

def _plot_states0(seq, paired_probabilities, states, *, seqfs=8, revseq=False,
                 vmin=0, vmax=1, cmap = inter_probability_cmap, no_label=False, state_names=None, **kwargs):
    
    """
    Show paired probabilities per state projected to one or the other dimension (dim)
    no_label: if True, do not show state labels on y-axis
    """
    
    seqlen = len(seq)
    matrix = np.zeros((len(states), seqlen+1))
    for sidx,state in enumerate(states):
        if state not in paired_probabilities: continue
        for i,prob in paired_probabilities[state].items():
            if revseq: ## revert sequence
                i = len(seq)-i+1
            matrix[sidx, i] = prob
        
    ax = sns.heatmap(matrix, linewidth=0.5, cmap=cmap, vmin=vmin, vmax=vmax, **kwargs)
    
    ax.set_xlim(1,seqlen+1)
    seqposticks = _seqpos_ticks(seqlen, 10)
    ax.set_xticks([x+0.5 for x in seqposticks],seqposticks,rotation=0)

    the_state_names = [(state_names[i] if state_names and i in state_names else i) for i in states]
    ax.set_yticks([x+0.5 for x in range(len(states))], the_state_names, rotation=0)

    if no_label:
        ax.set_yticklabels([])

    ax.grid(alpha=0.5)
        
    
    # Add sequence
    unit_cwidth, unit_cheight = _unit_csize(ax,seqfs,seqlen,len(states))
    
    for i,s in enumerate(seq):
        ax.text(i+1.5,-unit_cheight,s,fontsize=seqfs,name=seqfontname, ha="center")
    
    return ax


def _plot_states(seq, paired_probabilities, states, *, seqfs=8, revseq=False,
                  vmin=0, vmax=1, cmap=inter_probability_cmap, no_label=False,
                  state_names=None, row_gap=0.5, ax=None,
                  cbar=True, cbar_ax=None, cbar_kws=None,
                  state_colors=None, label_pad=1, **kwargs):


    if state_colors is None:
        palette = sns.color_palette(n_colors=len(states))
        state_colors = dict(zip(states, palette))

    seqlen = len(seq)
    n_states = len(states)

    matrix = np.zeros((n_states, seqlen + 1))
    for sidx, state in enumerate(states):
        if state not in paired_probabilities:
            continue
        for i, prob in paired_probabilities[state].items():
            if revseq:
                i = len(seq) - i + 1
            matrix[sidx, i] = prob

    if ax is None:
        ax = plt.gca()

    x_edges = np.arange(seqlen + 2)

    y_edges = [0.0]
    row_centers = []
    for i in range(n_states):
        bottom = y_edges[-1]
        top = bottom + 1
        row_centers.append((bottom + top) / 2)
        y_edges.append(top)
        if i < n_states - 1:
            y_edges.append(top + row_gap)
    y_edges = np.array(y_edges)

    n_ycells = len(y_edges) - 1
    full_matrix = np.full((n_ycells, seqlen + 1), np.nan)
    yi = 0
    for i in range(n_states):
        full_matrix[yi, :] = matrix[i, :]
        yi += 1
        if i < n_states - 1:
            yi += 1
    masked = np.ma.masked_invalid(full_matrix)

    plot_cmap = copy.copy(plt.get_cmap(cmap) if isinstance(cmap, str) else cmap)
    plot_cmap.set_bad(alpha=0)

    mesh = ax.pcolormesh(x_edges, y_edges, masked, cmap=plot_cmap,
                          vmin=vmin, vmax=vmax, edgecolors='white',
                          linewidth=0.5, **kwargs)

    ax.grid(alpha=0.5)

    # Remove frame around the heatmap axes
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Remove frame around the colorbar
    if cbar:
        cbar_kws = cbar_kws or {}
        if cbar_ax is not None:
            cb = ax.figure.colorbar(mesh, cax=cbar_ax, **cbar_kws)
        else:
            cb = ax.figure.colorbar(mesh, ax=ax, **cbar_kws)
        cb.outline.set_visible(False)

    ax.invert_yaxis()
    ax.set_xlim(1, seqlen + 1)

    seqposticks = _seqpos_ticks(seqlen, 10)
    ax.set_xticks([x + 0.5 for x in seqposticks], seqposticks, rotation=0)

    the_state_names = [(state_names[i] if state_names and i in state_names else i) for i in states]
    ax.set_yticks(row_centers, the_state_names, rotation=0)

    if no_label:
        ax.set_yticklabels([])
    else:
        ax.set_yticklabels([])
        from matplotlib.transforms import offset_copy
        trans = offset_copy(ax.transData, fig=ax.figure, x=-label_pad, y=0, units='points')
        the_state_names = [(state_names[i] if state_names and i in state_names else i) for i in states]
        for state, name, y in zip(states, the_state_names, row_centers):
            color = (state_colors.get(state, "grey") if state_colors else "grey")
            ax.text(x_edges[0], y, str(name), fontsize=seqfs, color="white",
                    fontweight="bold", ha="right", va="center",
                    transform=trans,
                    bbox=dict(facecolor=color, edgecolor="none",
                              boxstyle="square,pad=0.4"))

    from matplotlib.transforms import offset_copy
    trans = offset_copy(ax.transData, fig=ax.figure, x=0, y=seqfs, units='points')
    for i, s in enumerate(seq):
        ax.text(i + 1.5, y_edges[0], s, fontsize=seqfs, name=seqfontname,
                ha="center", va="bottom", transform=trans)

    return ax


# ### Plotting of positional probabilities

# +
def _total_event_probs(event_probs,state_probabilities):
    """
    Compute total pair probabilities given state-wise pair probabilities and state probabilities
    
    Args:
        event_probs: conditional probabilities of events for all states
        state_probabilities: probabilities of all states
    """
        
    total = dict()
    for state in event_probs:
        for event in event_probs[state]:
            total[event] = total.get(event,0) + event_probs[state][event] * state_probabilities[state]
    return total

def _total_event_probs_over_time(event_probs, state_probabilities):
    """Compute total event probabilities at all time points
    """
    return { time:_total_event_probs(event_probs, row.values)    
        for time,row in state_probabilities.iterrows() }


# -

def _select_pair_events_for_seq(bpp_events, seqno):
    """
    Filter pair events of all states for one sequence
    
    Args:
        bpp_events pair events for both sequences
        seqno sequence number
    Returns:
        pair events for single sequence
    """
    return { state:{ (i,j):p for (k,i,j),p in xs.items() if k==seqno }
             for state,xs in bpp_events.items() }


def _marginalize_pair_probabilities(pair_events, seqno=None):
    """
    Positionwise-marginalized pair probabilities for all states
    
    Args:
        pair_events: dict of dict of probabilities of base pairing events
        seqno: if not None, marginalize to only given sequence number
        
    Returns:
        positional probabilities per state
    """
    
    assert(seqno is None or seqno in [0,1])
    if seqno == None:
        seqnos = [0,1]
    else:
        seqnos = [seqno]
    
    def _marginalize_per_state(xs):
        ys = dict()
        for pair,p in xs.items():
            for k in seqnos:
                ys[pair[k]] = ys.get(pair[k],0) + p
        return ys

    return { state:_marginalize_per_state(xs) for state,xs in pair_events.items() }


def _plot_positional_probability_kinetics(seq, state_probabilities, positional_probabilities, *,
                                            seqfs=8, revseq=False, cmap=inter_probability_cmap, **kwargs):
    """
    Plot interaction probabilities over time
    
    Args:
        state_probabilities: data frame of probabilities of states at all time points
        positional_probabilities: positional probabilities in all states
    """
    
    probs_time = _total_event_probs_over_time(positional_probabilities, state_probabilities)
    
    seqlen=len(seq)
    
    timepoints=list(probs_time.keys())
    probs_time_matrix = np.zeros(shape=(seqlen+1,len(timepoints)), dtype=float)

    for i,t in enumerate(probs_time.keys()):
        for j in probs_time[t]:
            x = probs_time[t][j]
            jj=j-1
            if revseq: ## revert sequence
                jj=seqlen-j
            probs_time_matrix[jj,i] = x

        
    ax = sns.heatmap(probs_time_matrix, cmap=cmap, vmin=0, vmax=1, **kwargs)

    seqposticks = _seqpos_ticks(seqlen, 10)
    ax.set_yticks([x-0.5 for x in seqposticks],seqposticks)

    mintime=timepoints[0]
    maxtime=timepoints[-1]

    xtick_factor=100
    nxticks=int(np.log(maxtime)/np.log(xtick_factor))+1

    time_ticks = [round(np.log10(xtick_factor**i)) for i in range(nxticks)]
    def linearmap(x,frs,fre,tos,toe):
        return (x-frs)/(fre-frs)*(toe-tos)+tos

    lm = lambda x: linearmap(np.log10(x),np.log10(mintime),np.log10(maxtime),
                        0,len(timepoints))

    #ax.set_xlim(0,len(timepoints))
    ax.set_ylim(seqlen,0)

    ax.set_xticks([lm(10**i) for i in time_ticks],
       [f'$10^{{{i}}}$'
        for i in time_ticks],
        rotation=0
        )

    # Add sequence
    unit_cwidth, unit_cheight = _unit_csize(ax,seqfs,len(timepoints),seqlen)

    for i,s in enumerate(seq):
        ax.text(len(timepoints)+1*unit_cwidth,i+0.5,s,
                fontsize=seqfs, ha='center', va='center', name=seqfontname)

    ax.grid(alpha=0.5)
    
    return ax


# ### Dotplots

def _dotplot(seqA, seqB, state_probabilities, pair_probabilities, *,
             start_time=0, end_time=None, kind="maxbpp", states=None, reverseB=False, swap_axes=False, seqfs=8,
             cmap=dotplot_probability_cmap, vmin=0, vmax=None, ax_labels="xy", seq_labels="ab", **kwargs):

    """
    Dot plot at maximum state probabilities or by integrating probablities over log time from start_time to end_time (maximally until convergence)
    
    Note: the current implementation is resticted to start and end time in the range represented in state_probabilities
    
    maxbpp:     plot maximum probability of base pair over time
    maxsp:      maximum state probabilities
    integrate:  integrate over state probabilities
    """
    
    assert(kind in ["integrate","maxsp", "maxbpp"])
    
    lenA, lenB = len(seqA), len(seqB)
    matrix = np.zeros((lenA+1, lenB+1))
    
    if states is None:
        states = state_probabilities.keys()
        
    if kind in ["integrate", "maxsp"]:
        for state in states:
            stprobs = state_probabilities[state][start_time:end_time]

            if kind=="integrate":
                x = scipy.integrate.trapezoid(stprobs)
                x=x/(len(stprobs)-1)
            elif kind=="maxsp":
                x = max(stprobs)

            if state not in pair_probabilities: continue

            for pair,prob in pair_probabilities[state].items():
                i,j = pair
                if reverseB:
                    j = len(seqB)-j+1
                matrix[i,j] += x * prob

    elif kind=="maxbpp":
        for time,row in state_probabilities.iterrows():
            if not ( start_time <= time and (end_time is None or time < end_time) ): continue
            
            bpp_at_time = dict()
            
            for state in states:
                state_prob = row[state] # probability of state at time

                if state not in pair_probabilities: continue
                for pair,bp_prob in pair_probabilities[state].items():
                    i,j = pair
                    bpp_at_time[(i,j)] = bpp_at_time.get((i,j), 0) + state_prob * bp_prob
                    
            ## keep the maximax in matrix
            
            for (i,j),p in bpp_at_time.items():
                if reverseB:
                    j = len(seqB)-j+1
                matrix[i,j] = max(matrix[i,j], p)            
            
    if swap_axes:
        matrix=matrix.transpose()
    
    ax = sns.heatmap(matrix, cmap=cmap, vmin=vmin, vmax=vmax, **kwargs)
    
    if swap_axes:
        lenA,lenB,seqA,seqB = lenB,lenA,seqB,seqA
    
    # Coordinate system
    #        x/B
    #     +------->
    #     |
    # y/A |
    #     |
    #     v
    #
    
    
    ax.set_xlim(1,lenB+1)
    ax.set_ylim(1,lenA+1)
    
    range_x = _seqpos_ticks(lenB)
    range_y = _seqpos_ticks(lenA)
        
    ax.set_xticks([x+0.5 for x in range_x], range_x)
    ax.set_yticks([x+0.5 for x in range_y], range_y)
    ax.grid(alpha=0.5)

    ax.invert_yaxis()

    if not "x" in ax_labels:
        ax.set_xticklabels([])

    if not "y" in ax_labels:
        ax.set_yticklabels([])

    # Add sequences
        
    unit_cwidth, unit_cheight = _unit_csize(ax,seqfs,lenB,lenA)

    if "a" in seq_labels:
        for i,s in enumerate(seqA):
            ax.text(lenB + 1 + 0.5*unit_cwidth, i+1.5, s, fontsize=seqfs, name=seqfontname, ha="left", va="center")

    if "b" in seq_labels:
        for i,s in enumerate(seqB):
            ax.text(i+1.5, 1 - 0.5*unit_cheight, s, fontsize=seqfs,name=seqfontname, ha="center")
    
    return ax


# ## Code to handle intramolecular base pairs as events over coarsegraining
#
# Here, our objective (and motivation for some back-and-forth computation over coarse-grainin) is the computation of the probabilities of single base pairs in the final states for the purpose of visualization.
#
# Using McCaskill's algo implmented in the Vienna RNA package, one efficiently computes intramolecular base pair probabilities for (interaction regions of) the hybridization states.
#
# From the coarsegraining, we moreover memorize the coarsegraining matrices. Together with hybridization state energies, under the assumption of Boltzmnann distrbibution, this information defines the conditional the base pair probabilities in the final states.
#
# To limit the computation of base pair probabilities, we first determine the hybridization states that contribute substantially to the kinetics. For this purpose, we derive the probabilities of hybridization states, assuming that all final states occur with their maximum probability in the kinetics.
#
# In a second step, we approximate the probabilities of base pairs in the final states from the base pairs that occur with sufficient maximum probability in hybridization states.
#

# #### Reading of coarse graining matrices and state weights from files
#
# The rikin pipeline can write full tracking information of the entire coarse graining procedure to files.
#
# We read in these files in order to calculate the association between initial and final states. 

# +
def _read_prune_track(filename):
    """
    Read compressed prune track file
    
    Returns:
        sparse coarse graining matrix dict(final state index,dict(basin index, float)),
        where entries are the proportion of the basin in the final state
        """
    
    with gz.open(filename) as fh:
        cgmatrix = dict()
        i=0
        while True:
            line = fh.readline().decode()
            if not line: break
            line = line.strip().split('\t')
            line = [x.split(':') for x in line]
            line = {int(k):float(v) for k,v in line}
            cgmatrix[i] = line
            i += 1
        return cgmatrix


def _reorient_state(state,lenB):
    """
    Reorient state
    
    Reorients the position coordinates for sequence B, such that they are oriented
    in the opposite orientation as before.
    This is useful, since states are written to files using in 3'->5' indexation for sequence B.
    """
    i1,i2,j1,j2 = state
    return i1,lenB-j2+1,j1,lenB-i2+1

    
def _read_barriers_track(filename,lenA,lenB):
    """
    Read compressed barrier track file
    
    Returns:
        discrete sparse coarse graining matrix dict(basin index,dict(hybrid state, 1))
        only "1" entries are represented, others are 0
        AND
        sparse vector of partition functions of each basin
    
    Hybridization states are represented as tuples (i1,i2,j1,j2)
    where  i1..j1 is the interaction region in seqA 5'->3'
    and    i2..j2 is the interaction region in seqB 5'->3'
    
    NOTE: does not support double site states
    """
        
    def parse_state(s):
        if s=='{}':
            return ()
        else:
            m = re.match(r'\{\((\d+),(\d+)\)-\((\d+),(\d+)\)\}',s)    
            state = tuple([int(m[i]) for i in range(1,5)])
            return _reorient_state(state,lenB)
    with gz.open(filename) as fh:
        cgmatrix = dict()
        pfvec = dict()
        while True:
            line = fh.readline().decode()
            if not line: break
            line = line.strip().split('\t')
            states=[]
            idx = int(line[0])
            statesnum = int(line[1])
            representative = parse_state(line[3])
            
            # ignore these entries
            #representative_energy = float(line[4])
            energy = float(line[2])
            
            states = [representative]
            if statesnum > 1:
                 states.extend(parse_state(x) for x in line[5].split())
            
            assert(idx not in cgmatrix)
            cgmatrix[idx] = {k:1 for k in states}
            
            pfvec[idx] = boltzmann_weight(energy)
            
        return cgmatrix,pfvec

def _read_states(filename,lenA,lenB):
    """
    Read list of states from file
    
    Args:
        filename: name of input file
        lenB: length of sequence A
        lenB: length of sequence B
    
    Returns: dictionary of state energies indexed by tuples (i1,i2,j1,j2)
    where  i1..j1 is the interaction region in seqA 5'->3'
    and    i2..j2 is the interaction region in seqB 5'->3'
    """
    
    with gz.open(filename) as fh:
        states = dict()
        while True:
            line = fh.readline().decode()
            if not line: break
            line = line.strip().split('\t')
            if len(line)==5:
                state = tuple([int(x) for x in line[1:]])
                state = _reorient_state(state,lenB)
            else:
                assert(len(line)==1)
                state = ()
            states[state] = float(line[0])
        if () not in states:
            states[()] = -4.1 ## default hybridization energy
        return states


# -

# ### Calculation of partition functions after and probabilities before coarse graining

# +
def _coarsegrained_pfs( cgmatrix, zvec ):
    """
    Partition functions after coarse graining step
    
    Args:
        cgmatrix sparse coarse graining matrix
        zvec sparse vector of partitions functions of input states
    Result:
        vector of partition functions
    """
    zvec_p = dict()
    
    for alpha in cgmatrix:
        for x in cgmatrix[alpha]:
            if x in zvec:
                zvec_p[alpha] = zvec_p.get(alpha,0) + cgmatrix[alpha][x] * zvec[x]

    return zvec_p


def _backpropagate(cgmatrix, zvec, pvec_p, threshold=None):
    """
    Backpropagate coarse graining step
    
    Args:
        cgmatrix sparse coarse graining matrix
        zvec sparse vector of partitions functions of input states
        pvec sparse vector / dictionary of probabilities output states
        threshold: probability threshold for optional sparsification
    Result:
        probability vector of out states
    """
    pvec = dict()
    
    zvec_p = _coarsegrained_pfs( cgmatrix, zvec )
    
    for alpha in pvec_p:
        if alpha not in cgmatrix: continue
        for x in cgmatrix[alpha]:
            if x in zvec:
                z = cgmatrix[alpha][x] * zvec[x] /zvec_p[alpha]
                pvec[x] = pvec.get(x,0) + z * pvec_p[alpha]
        
    pvec = { alpha:z for alpha,z in pvec.items()}

    if threshold is not None:
        pvec = { alpha:p for alpha,p in pvec.items() if p>=threshold}

    return pvec


# -

# ### Calculation of the relevant intra-molecular base pair events in hybridization states and their probabilities

# +
def _relevant_hybstates(state_probabilities, prune_cgm, barriers_cgm, basin_pfs, state_pfs):
    max_probabilities = {i:max(row) for i,row in state_probabilities.items() if max(row)>=0.01}
    #print(max_probabilities)
    #print({k:p for k,p in prune_cgm[1].items() if p>=0.8})

    basin_pvec = _backpropagate(prune_cgm, basin_pfs, max_probabilities, 0.01)
    #print('basins', basin_pvec)

    hybstate_pvec = _backpropagate(barriers_cgm, state_pfs, basin_pvec, 0.01)
    #print('hybrid states', hybstate_pvec)

    return hybstate_pvec.keys()

def _hybstates_regions(hybstates):
    """
    The regions in a list of hybridization states
    
    Args:
        hybstates list of hybridization states
    Returns:
        Pair of sets of regions in seqA and seqB
        
    Note: includes regions () for the open state
    """
    
    def hybregions(hybstates, seqno):
        return set( [()]+[(x[seqno],x[seqno+2]) for x in hybstates if len(x)==4] )
    
    a, b = hybregions(hybstates, 0), hybregions(hybstates, 1)
    return (a,b)

def _conditional_basepair_probabilities(seq,region):
    fc = RNA.fold_compound(seq)
    
    if region != ():
        x,y = region
        for i in range(x,y+1):
            fc.hc_add_up(i)
    
    fc.pf()
    bpp = fc.bpp()
    return bpp

def _calculate_region_bpps(seqA, seqB, rhybstates):
    rregions = _hybstates_regions(rhybstates)

    seqs = (seqA,seqB)
    cbpps = dict()
    for i,regions in enumerate(rregions):
        cbpps[i] = dict()
        for region in regions:
            cbpps[i][region] = _conditional_basepair_probabilities(seqs[i],region)
    return cbpps

# generate a dictionary of 'events' of every relevant hybstate
# the propery vector contains one dictionary per state,
# the dict is another sparse vector, indexed by a event identifier and a probability of this event

def _bp_events_of_hs(hs, cbpps, threshold):    
    props = dict()
    
    for seqno in [0,1]:
        if len(hs)==4:
            s,e = (hs[seqno], hs[seqno+2])
            bpp = cbpps[seqno][(s,e)]
        elif len(hs)==0:
            bpp = cbpps[seqno][()]
        props.update( {(seqno,i,j):p for i,_ in enumerate(bpp) for j,p in enumerate(bpp[i]) if p>=threshold} )
            
    return props

def _calculate_bp_events(seqA, seqB, rhybstates):
    cbpps = _calculate_region_bpps(seqA, seqB, rhybstates)
    return { hs:_bp_events_of_hs(hs, cbpps, 0.01) for hs in rhybstates }


# -

# ### Probabilities of events in final states and total probabilities

def _cg_event_probabilities(events, *, cgmatrix, zvec, zvec_p):
    """
    Compute probabilities of events in coarse grainined states 
    Args:
        events: dictionary of dictionaries of events and their probability for every state 
        cgmatrix: sparse coarse graining matrix
        zvec: weights of the input states
        zvec_p: weight of the output states (after coarse graining)
    """
    events_p = dict()
    for alpha in cgmatrix:
        assert(alpha in zvec_p)
        for x in events:
            if x not in cgmatrix[alpha]: continue
            assert(x in zvec)
            if events[x]:
                events_p.setdefault(alpha,dict())
            for evt,p in events[x].items():
                events_p[alpha][evt] = events_p[alpha].get(evt,0) + cgmatrix[alpha][x] * p * zvec[x] / zvec_p[alpha]
                
    return events_p


# ## Class to represent the results of one Rikin pipeline run 

class RikinRun:
    """
    A run of the rikin pipeline and its results
    
    * accesses the data written to disk by rikin_pipeline
    * processes the data
    * plots the data in various ways
    
    Note: We infer intra-molecular base pair probabilities for the final states. From this, we determine the accumulated maximal probabilities of states over time and plot them together with the corresponding info for inter-molecular base pairs.
    """
    
    def __init__(self, seqA, seqB, input_directory, *, shown_state_threshold=0.01, autosave=None):
        self.seqA = seqA
        self.seqB = seqB
        self.input_directory = input_directory
        self.shown_state_threshold = shown_state_threshold
        self.autosave = autosave

    @staticmethod
    def save(name, suffix=suffix):
        
        if not isinstance(suffix, list):
            suffix = [suffix]
        
        if output_directory:
            for suf in suffix:
                print(f"Saving figure to {output_directory}/{name}{suf}")
                plt.savefig(f"{output_directory}/{name}{suf}")

        
    def optsave(self,tag):
        if self.autosave:
            RikinRun.save(f'{self.autosave}_{tag}')
        
    @property
    def seqs(self):
        return [self.seqA,self.seqB]

    @property
    def lenA(self):
        return len(self.seqA)

    @property
    def lenB(self):
        return len(self.seqB)
    
    @property
    def shown_states(self):
        return prominent_states(self.state_probabilities, self.shown_state_threshold)
    
    @property
    def state_probabilities(self):
        return self._cache("state_proabilities", 
                           lambda: self._read_state_probabilities(f"{self.input_directory}/kin"))
 
    @property
    def ipps(self):
        return self._cache("ipps",
                           lambda: self._read_ipps(f"{self.input_directory}/track-ipps-prune.gz"))
    
    @property
    def interaction_probabilities(self):
        return self.ipps[0]

    @property
    def states(self):
        return self._cache("states",
            lambda: _read_states(f'{self.input_directory}/sorted_states.gz', self.lenA, self.lenB))

    def _cache(self, name, fun):
        if not hasattr(self,"_data"):
            self._data=dict()

        if name not in self._data:
            self._data[name] = fun()
        
        return self._data[name]
    
    @staticmethod
    def _read_state_probabilities(filename):
    
        df = pd.read_csv(filename, sep=r"\s+", index_col=0, header=None)
        # 0-based state indices
        df.rename(mapper=lambda x:x-1, axis=1, inplace=True)
        
        return df
    
    @staticmethod
    def _read_ipps(filename):
        """
        Read compressed interaction pp track file

        The result is a dictionary of 0-based states;
         and entries for all (relevant) interaction base pairs.
        Each entry is the (conditional) probability of the base pair in
        this state.

        Base pairs have 1-based indices.
        """

        with gz.open(filename) as fh:
            ppdata = dict()
            pfs=dict()
            while True:
                line = fh.readline().decode()
                if not line: break
                line = line.strip().split('\t')

                pfs[int(line[0])] = float(line[1])

                ld = dict()
                for x in line[2:]:
                    x = x.split(' ')
                    ld[(int(x[0]),int(x[1]))]=float(x[2])
                ppdata[int(line[0])] = ld
            return ppdata, pfs
    
    @property
    def state_pfs(self):
        return { state:boltzmann_weight(e) for state,e in self.states.items() }
    
    @property
    def barriers_track(self):
        return self._cache("barriers_track",
                           lambda: _read_barriers_track(f'{self.input_directory}/barriers_track.gz',
                                                       self.lenA, self.lenB))  
    @property 
    def barriers_cgm(self):
        return self.barriers_track[0]
    
    @property
    def basin_pfs(self):
        return self.barriers_track[1]
    
    @property
    def prune_cgm(self):
        return self._cache("prune_cgm",
                          lambda: _read_prune_track(f'{self.input_directory}/prune_track.gz'))

    @property
    def final_pfs(self):
        return _coarsegrained_pfs( self.prune_cgm, self.basin_pfs )
        
        
    @property
    def bpp_final_events(self):
        def compute_them():
            # ... and calculate the probabilities of intra-mol bp events in these states
            rhybstates =  _relevant_hybstates(self.state_probabilities, self.prune_cgm, self.barriers_cgm,
                                              self.basin_pfs, self.state_pfs)

            bpp_hybstate_events = _calculate_bp_events(self.seqA, self.seqB, rhybstates)
            
            # calculate the probabilities of events in the final states
            bpp_basin_events = _cg_event_probabilities(bpp_hybstate_events, cgmatrix=self.barriers_cgm, zvec=self.state_pfs, zvec_p=self.basin_pfs)

            return _cg_event_probabilities(bpp_basin_events, cgmatrix=self.prune_cgm,
                                                              zvec=self.basin_pfs, zvec_p=self.final_pfs)

        return self._cache("bpp_final_events", compute_them)
    
    def plot_state_probability_kinetics(self, **kwargs):
        return _plot_state_probability_kinetics(self.state_probabilities, shown_states=self.shown_states, **kwargs)
        
    #def plot_interaction_probability_kinetics(self, dim, **kwargs):
    #    _plot_interaction_probability_kinetics(self.seqs[dim], self.state_probabilities, self.interaction_probabilities, dim=dim, **kwargs)

    def plot_states(self, figsize, seqfs, hspace=1, wspace=0.2, cbar_ratio=0.025, state_names=None, state_colors=None, row_gap=0.5, **kwargs):
        
        layout=[['i0','i1','ibar'],
                ['s0','s1','sbar']]
        
        fig,ax = plt.subplot_mosaic(layout,figsize=figsize,
                                    width_ratios=[self.lenA,self.lenB,cbar_ratio*(self.lenA+self.lenB)],
                                    gridspec_kw = dict(hspace=hspace,wspace=wspace),
                                    **kwargs)
        
        vmin, vmax = 0, 1
        
        for dim in [0,1]:
            interpaired = _marginalize_pair_probabilities(self.interaction_probabilities, dim)
            _plot_states(self.seqs[dim], interpaired, states=self.shown_states, revseq=(dim==1),
                         vmin = vmin, vmax = vmax,
                         state_names=state_names, state_colors=state_colors,  row_gap=row_gap,
                         ax=ax[f'i{dim}'], cbar_ax=ax['ibar'], no_label=dim==1, **kwargs)

            pair_events = _select_pair_events_for_seq(self.bpp_final_events, dim)
            intrapaired = _marginalize_pair_probabilities(pair_events)
            _plot_states(self.seqs[dim], intrapaired, states=self.shown_states,
                         vmin = vmin, vmax = vmax,
                         cmap=intra_probability_cmap,
                         state_names=state_names, state_colors=state_colors,  row_gap=row_gap,
                         ax=ax[f's{dim}'], cbar_ax=ax['sbar'], no_label=dim==1, **kwargs)
        
        self.optsave("states")
        return fig,ax
        
    def interaction_dotplot(self, swap_axes=False, **kwargs):
        return _dotplot(self.seqA, self.seqB, self.state_probabilities, self.interaction_probabilities,
                        cmap = inter_probability_cmap,
                        states=self.shown_states, reverseB=True, swap_axes=swap_axes, **kwargs)
       
    def intramolecular_dotplot(self, dim, **kwargs):
        pair_probabilities = _select_pair_events_for_seq(self.bpp_final_events, dim)
        return _dotplot(self.seqs[dim], self.seqs[dim], self.state_probabilities, pair_probabilities,
                        cmap = intra_probability_cmap, **kwargs)

    def plot_inaccessibility_kinetics(self, dim, **kwargs):
        pair_events = _select_pair_events_for_seq(self.bpp_final_events, dim)
        mar_pair_events = _marginalize_pair_probabilities(pair_events)
        return _plot_positional_probability_kinetics(self.seqs[dim], self.state_probabilities, mar_pair_events, **kwargs)
        
    def plot_interaction_probability_kinetics(self, dim, **kwargs):
        mar_pair_events = _marginalize_pair_probabilities(self.interaction_probabilities, dim)
        return _plot_positional_probability_kinetics(self.seqs[dim], self.state_probabilities,
                                                 mar_pair_events, **kwargs)
        
    def dotplot(self, *, figsize, barsize=0.5, seqfs=8, vmax=1, space=0.5, hspace=None, wspace=None, **kwargs):
        lenA,lenB = self.lenA, self.lenB
        bar_ratio = (lenA+lenB)*barsize/figsize[0]

        if space is not None:
            hspace = hspace if hspace is not None else space
            wspace = wspace if wspace is not None else space
        
        layout='''
        AI.b.c
        .B.b.c
        '''
        fig,ax = plt.subplot_mosaic(layout, figsize=figsize, 
                                    width_ratios=[lenA,lenB,bar_ratio/2,bar_ratio,bar_ratio/2,bar_ratio], height_ratios=[lenA,lenB],
                                    gridspec_kw = dict(hspace=hspace,wspace=wspace),
                                    **kwargs)

        self.intramolecular_dotplot(0, seqfs=seqfs,cbar=True, cbar_ax=ax['b'], vmax=vmax, ax=ax['A'], seq_labels="b", ax_labels="xy")
        self.intramolecular_dotplot(1, seqfs=seqfs, cbar=False, vmax=vmax,ax=ax['B'], seq_labels="a", ax_labels="xy")
        self.interaction_dotplot(seqfs=seqfs, vmax=vmax, cbar_ax=ax['c'], ax=ax['I'], seq_labels="ab", ax_labels="")

        self.optsave("dotplot")
        return fig,ax
    
    def plot_paired_probability_kinetics(self, figsize, seqfs, space=0.5, hspace=None, wspace=None, **kwargs):
        layout=[['s0','i0','sbar','ibar',],
                ['s1','i1','sbar','ibar',]]
        
        if space is not None:
            hspace = hspace if hspace is not None else space
            wspace = wspace if wspace is not None else space
        
        fig,ax = plt.subplot_mosaic(layout,figsize=figsize,
                                    width_ratios=[8,8,0.5,0.5],
                                    height_ratios=[self.lenA,self.lenB],
                                    gridspec_kw=dict(hspace=hspace,wspace=wspace),
                                    **kwargs)

        for k in [0,1]:
            self.plot_inaccessibility_kinetics(k,
                seqfs=seqfs, cmap=intra_probability_cmap, ax=ax[f's{k}'], cbar_ax=ax['sbar'])

        for k in [0,1]:
            self.plot_interaction_probability_kinetics(k,
                ax=ax[f'i{k}'],
                revseq=(k==1),
                seqfs=seqfs,
                cmap=inter_probability_cmap,
                cbar_ax=ax['ibar'])
        
        self.optsave("paired_kinetics")
        return fig,ax


