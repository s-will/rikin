# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import gzip as gz
from collections import defaultdict

# +
fontname = "FreeMono"

#import matplotlib.font_manager
#matplotlib.font_manager.findSystemFonts(fontpaths=None)
# -

# # Plotting code

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

#probability_cmap = sns.cubehelix_palette(start=.5, rot=-.75, dark=0, light=1.0, as_cmap=True)
#probability_cmap = sns.cubehelix_palette(start=.5, rot=-.5, dark=0, light=1.0, as_cmap=True)
probability_cmap = sns.cubehelix_palette(start=1.5, rot=-0.75, dark=0.15, light=1.0, hue=1, gamma=1.5, as_cmap=True)
#probability_cmap = sns.cubehelix_palette(start=0, rot=1, dark=0.15, light=1.0, as_cmap=True)
probability_cmap


# -

# ### Reading Rikin pipeline results and plotting over time 

# +
def read_state_probabilities(filename):
    df = pd.read_csv(filename, sep="\s+", index_col=0, header=None)
    # 0-based state indices
    df.rename(mapper=lambda x:x-1, axis=1, inplace=True)
    return df

def highest_states(state_probabilities, maxp_threshold=0.02):
    # filter states by their maximum probabilty maxp_threshold at any time
    return [ state for state in state_probabilities.columns 
        if max(state_probabilities[state])>=maxp_threshold]

def plot_state_probabilities_over_time(state_probabilities, *, show_states=None, ax=None):
    
    important_states = show_states if show_states!=None else state_probabilities.columns()
    
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


# -

def read_ipps(filename):
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


# +
def total_pps(pps,state_probabilities):
    """
    Compute total pair probabilities given state-wise pair probabilities and state probabilities
    
    Args:
        pps: conditional pair probababilities for all states
        state_probabilities: probabilities of all states
    """
        
    total = defaultdict(float)
    for state in pps:
        for pair in pps[state]:
            total[pair] += pps[state][pair] * state_probabilities[state]
    return total

def total_pps_over_time(pps, state_probabilities):
    """Compute total probabilities at all time points"""
    return { time:total_pps(pps, row.values)    
        for time,row in state_probabilities.iterrows() }

def marginalize_pps(pps,dim=0):
    """
    Marginalize pair probabilities to position-wise paired
    probabilities in the given dimension
    """
    xs = defaultdict(float)
    for pair in pps:
        xs[pair[dim]] += pps[pair]
    return xs


# -

def plot_interaction_probabilities_over_time(seq, state_probabilities, interaction_probabilities, *, seqfs=8, ax=None, dim=0, revseq=False):
    """
    Plot interaction probabilities over time
    
    Args:
        state_probabilities: data frame of probabilities of states at all time points
        interaction_probabilities: interaction base pair probabilities in all states
    """
    
    pps_time = total_pps_over_time(interaction_probabilities, state_probabilities)
    mpps_time = {k:marginalize_pps(pps,dim) for k,pps in pps_time.items()}
    
    seqlen=len(seq)
    
    timepoints=list(mpps_time.keys())
    mpps_time_matrix = np.zeros(shape=(seqlen+1,len(timepoints)), dtype=float)

    for i,t in enumerate(pps_time.keys()):
        for j in mpps_time[t]:
            x = mpps_time[t][j]
            jj=j-1
            if revseq: ## revert sequence
                jj=seqlen-j
            mpps_time_matrix[jj,i] = x

        
    ax = sns.heatmap(mpps_time_matrix, cmap=probability_cmap, ax=ax)

    ax.set_yticks([0.5]+[x-0.5 for x in range(10,seqlen+1,10)],[1]+list(range(10,seqlen+1,10)))

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
    for i,s in enumerate(seq):
        ax.text(len(timepoints)+2.5,i+0.5,s,
                fontsize=seqfs, ha='center', va='center', name=fontname)

    ax.grid(alpha=0.5)
    
    return ax


# ### Plotting states

def plot_states(seq, dim, state_probabilities, interaction_probabilities, states, *, seqfs=8, revseq=False, cbar=True, ax=None):
    
    """
    Show paired probabilities per state projected to one or the other dimension (dim)
    """
    
    seqlen = len(seq)
    matrix = np.zeros((len(states), seqlen+1))
    for sidx,state in enumerate(states):
        for pair,prob in interaction_probabilities[state].items():
            i = pair[dim]
            if revseq: ## revert sequence
                i = len(seq)-i+1
            matrix[sidx, i] = prob
    
        
    cmap = probability_cmap
    
    ax = sns.heatmap(matrix, cmap=cmap, linewidth=0.5, cbar=cbar, ax=ax)
    ax.set_xlim(1,seqlen+1)
    ax.set_xticks([1.5]+[x+0.5 for x in range(10,seqlen+1,10)],[1]+list(range(10,seqlen+1,10)),rotation=0)
    ax.set_yticks([x+0.5 for x in range(len(states))], states, rotation=0)
    ax.grid(alpha=0.5)
    
    # Add sequence
    for i,s in enumerate(seq):
        ax.text(i+1.5,-0.5,s,fontsize=seqfs,name=fontname, ha="center")


def plot_accmaxbpprobs(seqA, seqB, state_probabilities, interaction_probabilities, states, *, seqfs=8, cbar=True, ax=None):
    
    """
    Experimentell: akkumuliere pair probs at max prob of state
    
    Idea for extensions: we could integrate over time instead of using maximum probability
    """
    
    lenA, lenB = len(seqA), len(seqB)
    matrix = np.zeros((lenA+1, lenB+1))
    for state in states:
        max_state_prob = max(state_probabilities[state])
        for pair,prob in interaction_probabilities[state].items():
            i,j = pair
            matrix[i,len(seqB)-j+1] += prob * max_state_prob
    
    
    matrix=matrix.transpose()
    
    ax = sns.heatmap(matrix, cmap=probability_cmap, ax=ax, cbar=cbar)
    ax.set_xlim(1,lenA+1)
    ax.set_ylim(1,lenB+1)
    ax.set_xticks([x+0.5 for x in range(10,lenA+1,10)],list(range(10,lenA+1,10)))
    ax.set_yticks([x+0.5 for x in range(10,lenB,10)],list(range(10,lenB+1,10)))
    ax.grid(alpha=0.5)

    ax.invert_yaxis()    
        
    # Add sequences
    for i,s in enumerate(seqA):
        ax.text(i+1.5,-0.25 ,s,fontsize=seqfs,name=fontname, ha="center")
    for i,s in enumerate(seqB):
        ax.text(lenA+1, i+1.5,s,fontsize=seqfs,name=fontname, ha="left", va="center")

# # Result plots

# ## MicA -- MicA

# +
input_directory = "/home/will/Research/Projects/RIKinetics/Experiments/NEW/MicA/MicA-MicA"

seqA="GAAAGACGCGCAUUUGUUAUCAUCAUCCCUGAAUUCAGAGAUGAAAUUUUGGCCACUCACGAGUGGCCUUUU"
seqB=seqA

# +
state_probabilities = read_state_probabilities(f"{input_directory}/kin")
shown_states = highest_states(state_probabilities, 0.02)

#print(state_probabilities)
fig,axs = plt.subplots(1,3,figsize=(24,8))

plot_state_probabilities_over_time(state_probabilities, show_states=shown_states, ax=axs[0])
#ax.set_xlim(0.1,1e10)

interaction_probabilities, eqpfs = read_ipps(f"{input_directory}/track-ipps-prune.gz")
plot_interaction_probabilities_over_time(seqA, state_probabilities, interaction_probabilities, dim=0, ax=axs[1])

plot_interaction_probabilities_over_time(seqB, state_probabilities, interaction_probabilities, dim=1,
                                         revseq=True, ax=axs[2])


#plt.savefig(f"{input_directory}/kinetics-plot.svg")
plt.show()
# -

fig,ax = plt.subplots(1,2,figsize=(16,4),width_ratios=[len(seqA),len(seqB)+2])
plot_states(seqA, 0, state_probabilities, interaction_probabilities, shown_states,ax=ax[0], cbar=False)
plot_states(seqB, 1, state_probabilities, interaction_probabilities, shown_states, revseq=True,ax=ax[1])
plt.show()

fig = plt.figure(figsize=(7,6))
plot_accmaxbpprobs(seqA, seqB, state_probabilities, interaction_probabilities, shown_states, seqfs=6)
plt.show()

# ## Example

seqA="AAAGGGGGGAAAAAAAGGGUGGGAAAAAAAGGGCGGGAAA"
seqB="CCCGCCC"
input_directory = "/home/will/Research/Projects/RIKinetics/RIkin/Doc/example"

# +
state_probabilities = read_state_probabilities(f"{input_directory}/kin")
shown_states = highest_states(state_probabilities, 0.02)

#print(state_probabilities)
fig,axs = plt.subplots(1,3,figsize=(20,6))

plot_state_probabilities_over_time(state_probabilities, show_states=shown_states, ax=axs[0])
#ax.set_xlim(0.1,1e10)

interaction_probabilities, eqpfs = read_ipps(f"{input_directory}/track-ipps-prune.gz")
plot_interaction_probabilities_over_time(seqA, state_probabilities, interaction_probabilities, dim=0, seqfs=10, ax=axs[1])

plot_interaction_probabilities_over_time(seqB, state_probabilities, interaction_probabilities, dim=1,
                                         revseq=True, ax=axs[2], seqfs=10)

plt.show()
# -

fig,ax = plt.subplots(1,2,figsize=(10,1.5),width_ratios=[len(seqA),len(seqB)+2])
plot_states(seqA, 0, state_probabilities, interaction_probabilities, shown_states,ax=ax[0], cbar=False)
plot_states(seqB, 1, state_probabilities, interaction_probabilities, shown_states, revseq=True,ax=ax[1])
plt.show()

fig = plt.figure(figsize=(8,1.25))
plot_accmaxbpprobs(seqA, seqB, state_probabilities, interaction_probabilities, shown_states)
plt.show()

# ## OmpA -- MicA

# +
input_directory = "/home/will/Research/Projects/RIKinetics/Experiments/NEW/MicA/ompA-MicA"

seqA="CUUUUUUUUCAUAUGCCUGACGGAGUUCACACUUGUAAGUUUUCAACUACGUUGUAGACUUUACAUCGCCAGGGGUGCUCGGCAUAAGCCGAAGAUAUCGGUAGAGUUAAUAUUGAGCAGAUCCCCCGGUGAAGGAUUUAACCGUGUUAUCUCGUUGGAGAUAUUCAUGGCGUAUUUUGGAUGAUAACGAGGCGCAAAAAAUGAAAAAGACAGCUAUCGCGAUUGCAGUGGCACUGGCUGGUUUCGCUACCGUAGCGCAGGCCGCUCCGAAAGAUAACACCUGGUACACUGGUGCUAAAC"
seqB="GAAAGACGCGCAUUUGUUAUCAUCAUCCCUGAAUUCAGAGAUGAAAUUUUGGCCACUCACGAGUGGCCUUUU"

# +
state_probabilities = read_state_probabilities(f"{input_directory}/kin")
shown_states = highest_states(state_probabilities, 0.02)

#print(state_probabilities)
fig,axs = plt.subplots(1,3,figsize=(20,6))

plot_state_probabilities_over_time(state_probabilities, show_states=shown_states, ax=axs[0])
#ax.set_xlim(0.1,1e10)

interaction_probabilities, eqpfs = read_ipps(f"{input_directory}/track-ipps-prune.gz")
plot_interaction_probabilities_over_time(seqA, state_probabilities, interaction_probabilities, dim=0, seqfs=4, ax=axs[1])

plot_interaction_probabilities_over_time(seqB, state_probabilities, interaction_probabilities, dim=1,
                                         revseq=True, ax=axs[2], seqfs=6)

plt.show()
# -

fig,ax = plt.subplots(1,2,figsize=(30,3),width_ratios=[len(seqA),len(seqB)+2])
plot_states(seqA, 0, state_probabilities, interaction_probabilities, shown_states,ax=ax[0], cbar=False)
plot_states(seqB, 1, state_probabilities, interaction_probabilities, shown_states, revseq=True,ax=ax[1])
plt.show()

fig = plt.figure(figsize=(20,4))
plot_accmaxbpprobs(seqA, seqB, state_probabilities, interaction_probabilities, shown_states, seqfs=5)
plt.show()

# ## KHP plots

# +
input_dir = "/home/will/Research/Projects/RIKinetics/Experiments/NEW/KHP-Kinetics"

seqs = {"hp1": "GGACGAGGCAUUUCCCCUUGU",
    "hp2": "GGACAAGGGGAAAUGCCUUGU",
    "hp3": "GGACGAUCAGCAUUUCCCUGAUGU"}
# -

for a,b in [(1,2),(2,3)]:
    input_directory = f"{input_dir}/HP{a}-HP{b}"
    seqA = seqs[f"hp{a}"]
    seqB = seqs[f"hp{b}"]

    print(f"HP{a} - HP{b}")
    
    state_probabilities = read_state_probabilities(f"{input_directory}/kin")
    shown_states = highest_states(state_probabilities, 0.02)

    #print(state_probabilities)
    fig,axs = plt.subplots(1,3,figsize=(20,6))

    plot_state_probabilities_over_time(state_probabilities, show_states=shown_states, ax=axs[0])
    #ax.set_xlim(0.1,1e10)

    interaction_probabilities, eqpfs = read_ipps(f"{input_directory}/track-ipps-prune.gz")
    plot_interaction_probabilities_over_time(seqA, state_probabilities, interaction_probabilities, dim=0, seqfs=10, ax=axs[1])

    plot_interaction_probabilities_over_time(seqB, state_probabilities, interaction_probabilities, dim=1,
                                             revseq=True, ax=axs[2], seqfs=10)

    plt.show()
    
    fig,ax = plt.subplots(1,2,figsize=(8,1),width_ratios=[len(seqA),len(seqB)+2])
    plot_states(seqA, 0, state_probabilities, interaction_probabilities, shown_states,ax=ax[0], cbar=False)
    plot_states(seqB, 1, state_probabilities, interaction_probabilities, shown_states, revseq=True,ax=ax[1])
    plt.show()
    
    fig = plt.figure(figsize=(2.25,2))
    plot_accmaxbpprobs(seqA, seqB, state_probabilities, interaction_probabilities, shown_states, seqfs=5)
    plt.show()


# ## APPENDIX / Unused code

def read_prune_track(filename):
    """
    Read compressed prune track file
    
    Returns list of dictionaries, indices(!?) and keys are 0-based state indices
    The values prune_track[i][j] are the proportion of state j in coarse grained state i.
    """
    
    with gz.open(filename) as fh:
        data = list()
        while True:
            line = fh.readline().decode()
            if not line: break
            line = line.strip().split('\t')
            line = [x.split(':') for x in line]
            line = {int(k):v for k,v in line}
            data.append(line)
        return data
