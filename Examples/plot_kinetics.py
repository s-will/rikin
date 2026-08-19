# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: infrared
#     language: python
#     name: python3
# ---

# # Kinetics plots

# +

#### autoreload trick from https://switowski.com/blog/ipython-autoreload
##
# autoreload is even more elegant than reload via importlib:
# %load_ext autoreload
# %autoreload 2

import rikinplotlib as rkplt
from rikinplotlib import *



# +
rkplt.output_directory = "Figs"
#suffix = ".svg"

if rkplt.output_directory and not os.path.isdir(rkplt.output_directory):
    os.mkdir(rkplt.output_directory)
        
# # Results --- generate plots

# ## Example

example = RikinRun(
    "AAAGGGGGGAAAAAAAGGGUGGGAAAAAAAGGGCGGGAAA",
    "CCCGCCC",
    "/home/will/Research/Projects/RIKinetics/RIkin/Doc/example",
    shown_state_threshold = 0.02,
    autosave="example"
)


# +
#-----

fig,axs = plt.subplots(1,1,figsize=(5,5))
example.plot_state_probability_kinetics(ax=axs)
example.optsave("kinetics")

plt.show()

#-----

example.plot_states(figsize=(10,2.5),seqfs=6)

plt.show()

#-----

fig,ax = plt.subplots(1,2,figsize=(8,1.75),width_ratios=[30,1])
example.interaction_dotplot(swap_axes=True,vmax=1,ax=ax[0],cbar_ax=ax[1])
fig.tight_layout()

example.optsave("interaction")
plt.show()

#-----

example.dotplot(figsize = (11,10), seqfs = 8, barsize=0.25, space=0.25)
plt.show()

#-----

example.plot_paired_probability_kinetics(figsize = (12,6), seqfs = 7, hspace=0.25, wspace=0.25)
plt.show()
# -


# ## MicA -- MicA
#
# *Note*: For homodimers, the plots for seqA and seqB are necessarily symmetric; we could ommit plotting the redundant info for both sequences for a more terse presentation. This is not implemented yet (**TODO**).

# +
seq="GAAAGACGCGCAUUUGUUAUCAUCAUCCCUGAAUUCAGAGAUGAAAUUUUGGCCACUCACGAGUGGCCUUUU" 
micA_micA = RikinRun(
    seq, seq, input_directory = "/home/will/Research/Projects/RIKinetics/Experiments/NEW/MicA/MicA-MicA",
    shown_state_threshold = 0.05,
    autosave="micA_micA"
)

seqlen = len(seq)

fig,axs = plt.subplots(1,1,figsize=(6,6))
micA_micA.plot_state_probability_kinetics(ax=axs)

micA_micA.optsave("kinetics")
plt.show()

# -----

#-----

micA_micA.plot_states(figsize=(12,4), seqfs=6, wspace=0.15, hspace=0.6,row_gap=1.2)
plt.show()

# -----

#fig = plt.figure(figsize=(7,6))
#micA_micA.interaction_dotplot(seqfs=6)
#plt.show()

#-----

micA_micA.plot_paired_probability_kinetics(figsize = (12,4), seqfs = 6, hspace=0.25, wspace=0.25, show_rnas=[1], suffix="paired_kinetics")
plt.show()

# -

#-----
micA_micA.dotplot(figsize = (13.5,12), seqfs = 6, barsize=0.5, wspace=0.125, hspace=0.075)
plt.show()

# ## OmpA -- MicA

# +
ompA_micA = RikinRun(
    input_directory = "/home/will/Research/Projects/RIKinetics/Experiments/NEW/MicA/ompA-MicA",
    seqA="CUUUUUUUUCAUAUGCCUGACGGAGUUCACACUUGUAAGUUUUCAACUACGUUGUAGACUUUACAUCGCCAGGGGUGCUCGGCAUAAGCCGAAGAUAUCGGUAGAGUUAAUAUUGAGCAGAUCCCCCGGUGAAGGAUUUAACCGUGUUAUCUCGUUGGAGAUAUUCAUGGCGUAUUUUGGAUGAUAACGAGGCGCAAAAAAUGAAAAAGACAGCUAUCGCGAUUGCAGUGGCACUGGCUGGUUUCGCUACCGUAGCGCAGGCCGCUCCGAAAGAUAACACCUGGUACACUGGUGCUAAAC",
    seqB="GAAAGACGCGCAUUUGUUAUCAUCAUCCCUGAAUUCAGAGAUGAAAUUUUGGCCACUCACGAGUGGCCUUUU",
    shown_state_threshold = 0.05,
    autosave="ompA_micA"
)

#print(state_probabilities)
fig,axs = plt.subplots(1,1,figsize=(6,6))
ompA_micA.plot_state_probability_kinetics(ax=axs)

ompA_micA.optsave("kinetics")
plt.show()

# -----


ompA_micA.plot_states(figsize=(14,8), seqfs=6, wspace=0.15, hspace=0.3 )
plt.show()

# -----

#fig = plt.figure(figsize=(20,4))
#ompA_micA.interaction_dotplot(seqfs=5, swap_axes=True)
#plt.show()

# -----

#fig = plt.figure(figsize=(20,4))
#ompA_micA.interaction_dotplot(kind="integrate", start_time=0, end_time=None, seqfs=5, swap_axes=True)
#plt.show()

ompA_micA.plot_paired_probability_kinetics(figsize=(26,24), seqfs=5, wspace=0.2, hspace=0.07)
plt.show()

ompA_micA.dotplot(figsize = (26,24), seqfs = 8, barsize=0.5, wspace=0.12, hspace=0.08)
plt.show()
# -

ompA_micA.plot_paired_probability_kinetics(figsize=(26,24), seqfs=5, wspace=0.2, hspace=0.07, show_pairing=False, suffix="paired_kinetics_interaction_only")
plt.show()

# +
ompA_micA.plot_paired_probability_kinetics(figsize = (12,4), seqfs = 6, hspace=0.25, wspace=0.25, show_rnas=[1], suffix="paired_kinetics_only_srna")

plt.show()
# -

# ## KHP plots

# +
input_dir = "/home/will/Research/Projects/RIKinetics/Experiments/NEW/KHP-Kinetics"

seqs = {
    "hp1": "GGACGAGGCAUUUCCCCUUGU",
    "hp2": "GGACAAGGGGAAAUGCCUUGU",
    "hp3": "GGACGAUCAGCAUUUCCCUGAUGU"}

# HP1 GGACGAG  GCAUUUCCC CUUGU
# HP3 GGACGAUCAGCAUUUCCCUGAUGU



for a,b in [(1,2),(2,3)]:
    khp = RikinRun(
        input_directory = f"{input_dir}/HP{a}-HP{b}",
        seqA = seqs[f"hp{a}"],
        seqB = seqs[f"hp{b}"],
        shown_state_threshold = 0.02,
        autosave=f"HP{a}-HP{b}"
    )
    
    print("------------------")
    print(f"HP{a} - HP{b}")
    
    fig,ax = plt.subplots(1,1,figsize=(6,6))
    khp.plot_state_probability_kinetics(ax=ax)
    khp.optsave("kinetics")
    plt.show()
    
    khp.plot_states(figsize=(8,3),seqfs=8, wspace=0.6, row_gap=1)
    plt.show()
    
    khp.plot_paired_probability_kinetics(figsize=(8,8), seqfs=8,wspace=0.5,hspace=0.2)
    plt.show()
    
    khp.dotplot(figsize=(7.5,6),seqfs=8,space=0.35)
    plt.show()

# +
state_names = {0:"ED", 1:"Diss", 2:"KC1", 3:"KC2", 9:"OFF1", 15:"OFF2"}

a,b = 1,2
khp = RikinRun(
    input_directory = f"{input_dir}/HP{a}-HP{b}",
    seqA = seqs[f"hp{a}"],
    seqB = seqs[f"hp{b}"],
    shown_state_threshold = 0.02,
    autosave=f"HP{a}-HP{b}"
)

rna_names = [f"HP{a}", f"HP{b}"]

print("------------------")
print(f"HP{a} - HP{b}")

khp.plot_states(figsize=(8,4.5),seqfs=8, rna_names=rna_names, state_names=state_names,row_gap=1,hspace=0.5,wspace=0.25)
plt.show()


# +
state_names = {0:"ED", 1:"Diss", 2:"KC1", 3:"KC2", 14:"OFF1", 29:"OFF2"}

a,b = 2,3
khp = RikinRun(
    input_directory = f"{input_dir}/HP{a}-HP{b}",
    seqA = seqs[f"hp{a}"],
    seqB = seqs[f"hp{b}"],
    shown_state_threshold = 0.02,
    autosave=f"HP{a}-HP{b}"
)

rna_names = [f"HP{a}", f"HP{b}"]

print("------------------")
print(f"HP{a} - HP{b}")

khp.plot_states(figsize=(8,4.5),seqfs=8, rna_names=rna_names, state_names=state_names,row_gap=1,hspace=0.5,wspace=0.25)
plt.show()

# -

xs= np.linspace(0.1, 4, 50)
sns.lineplot(x=xs,y=[boltzmann_weight(x) for x in xs])

boltzmann_weight(3)


