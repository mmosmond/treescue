# ------------ python modules ----------

import pyslim, msprime, tskit
import numpy as np
import os

# ------------- paths ------------------

DATADIR = 'data/' #where to put outputs
PROGRAMDIR = 'programs/'
SLIM = PROGRAMDIR + 'SLiM_3.7.1/slim' #command to run SLiM
RELATE = PROGRAMDIR + 'relate'
CLUES = PROGRAMDIR + 'clues'

#  ----------- parameters -----------------

# for simulations 
Ks = [int(1e4)] #carrying capacity
ds = [0] #decline rate of ancestral homozygote: absolute fitness 1-d
ss = [0.01] #selection coefficient of derived allele: absolute fitness of derived homozygote (1-d)(1+s)
hs = [0.5] #dominance coefficient of derived allele: absolute fitness of the heterozygote (1-d)(1+hs)
Bs = [2] #number of offspring per parent (note we divide viability by this so that it does not affect the expected number of offspring (absolute fitness), but it will affect the variance in the number of offspring (Ne))
us = [0] #mutation rate at selected site (will add neutral mutations to other sites after with msprime)
qs = [1] #number of beneficial mutations at time 0
Ls = [int(1e5)] #number of sites minus 1 (note SLiM 0 indexes sites; will put selected site in middle) 
rs = [1.25e-8] #recombination rate per site
fs = [0.75] #beneficial allele frequency to stop at
ns = range(1) #replicates
t = int(1e4) #max number of generations (to prevent simulations from running forever) 
# for vcfs
ks = [25] #number of individuals to sample (get twice as many genomes)
Us = [2.5e-8] #mutation rate at neutral sites

# ------- simulations -------

sim_output = DATADIR + "sim_{K}K_{d}d_{s}s_{h}h_{B}B_{u}u_{q}q_{L}L_{r}r_{f}f_{n}n.{END}"
sim_ends = ['dynamics','trees']

rule simulate_all:
  input:
    expand(sim_output, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, q=qs, L=Ls, r=rs, f=fs, n=ns, END=sim_ends)

rule simulate:
  input:
#    "scripts/sim.slim"
    "scripts/simWF.slim" #testing with WF
  output:
    expand(sim_output, END=sim_ends, allow_missing=True)
  shell:
    """
    module load gcc/8.3.0 #needed for slim
    mkdir -p {DATADIR}
    {SLIM} \
      -d K={wildcards.K} \
      -d d={wildcards.d} \
      -d s={wildcards.s} \
      -d h={wildcards.h} \
      -d B={wildcards.B} \
      -d u={wildcards.u} \
      -d q={wildcards.q} \
      -d L={wildcards.L} \
      -d r={wildcards.r} \
      -d f={wildcards.f} \
      -d t={t} \
      -d "output_dynamics='{output[0]}'" \
      -d "output_trees='{output[1]}'" \
      {input}
    """

# ---------- trees and vcfs --------

sim_trees = sim_output.replace('.{END}', '_{k}k_{U}U.{END}')
TREES = ['trees', 'vcf']

rule trees_all:
  input:
    expand(sim_trees, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, q=qs, L=Ls, r=rs, f=fs, n=ns, k=ks, U=Us, END=TREES)
 
rule trees:
  input:
    expand(sim_output, END=['trees'], allow_missing=True)
  output:
    expand(sim_trees, END=TREES, allow_missing=True)
  params:
#    Ne = lambda wildcards: int(wildcards.K)*4/(2+4*int(wildcards.B)-3)
    Ne = lambda wildcards: int(wildcards.K) #for WF test
  run:
    ts = tskit.load(input[0]) #load tree sequence
    inds = np.random.choice(range(ts.num_individuals), int(wildcards.k), replace=False) #sample k individuals at random
    ind_nodes = [] #empty vector for nodes
    for ind in inds: #for each individual
      for node in ts.individual(ind).nodes:
        ind_nodes.append(node) #append each genome
    ts = ts.simplify(ind_nodes, keep_input_roots=True) #simplify down to sample while keeping roots for recapitation
    ts = pyslim.update(ts) #using beta pyslim for slim 4.0 but trees made with slim 3.7.1, so update
    ts = pyslim.recapitate(ts, recombination_rate=float(wildcards.r), ancestral_Ne=params.Ne) #recapitate trees
    #ts = msprime.sim_mutations(ts, rate=float(wildcards.U), model=msprime.SLiMMutationModel(type=0), keep=True) #layer on neutral mutations while keeping the selected mutation
    #ts = msprime.sim_mutations(ts, rate=float(wildcards.U), model=msprime.BinaryMutationModel(), keep=True) #layer on neutral mutations while keeping the selected mutation
    ts = pyslim.generate_nucleotides(ts) #generate random nucleotides for slim mutations, https://github.com/tskit-dev/pyslim/pull/174
    ts = pyslim.convert_alleles(ts) #convert slim alleles (0,1) to nucleotides
    ts = msprime.sim_mutations(ts, rate=float(wildcards.U), model=msprime.JC69(), keep=True) #layer on neutral mutations while keeping mutation at selected site (set keep=False if selected allele fixed since then not biallelic)
    ts.dump(output[0]) #save true treesequence
    with open(output[1], 'w') as vcffile: 
      ts.write_vcf(vcffile, individuals = range(int(wildcards.k))) #write out as VCF with arbitrary ID names for Relate
    #ts.write_fasta(output[2]) #write as fasta for argweaver, nevermind, will just convert vcf to sites

# ------------ vcfs to haps -------------

sim_hapsample = sim_trees
HAPSAMPLE = ['haps','sample']

rule haps_all:
  input:
    expand(sim_hapsample, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, q=qs, L=Ls, r=rs, f=fs, n=ns, k=ks, U=Us, END=HAPSAMPLE)

rule haps:
  input:
    expand(sim_trees, END=['vcf'], allow_missing=True)
  output:
    expand(sim_hapsample, END=HAPSAMPLE, allow_missing=True)
  params:
    prefix=sim_trees.replace('.{END}','')
  shell:
    '''
    module load gcc/8.3.0 #needed for relate
    {RELATE}/bin/RelateFileFormats \
                 --mode ConvertFromVcf \
                 --haps {output[0]} \
                 --sample {output[1]} \
                 -i {params.prefix} \
                 --chr 1
    '''

# --------------- recombination map ---------

sim_map = DATADIR + 'sim_{L}L_{r}r.map'

rule map_all:
  input:
    expand(sim_map, L=Ls, r=rs)

rule map:
  input:
  output:
    sim_map
  run:
    L = int(wildcards.L)
    R = (1 - (1 - 2 * float(wildcards.r))**L)/2 #recombination distance from one end of chromosome to other
    cm = 50 * np.log(1/(1-2*R)) #length in centiMorgans
    cr = cm/L * 1e6 #cM per million bases
    script = "pos COMBINED_rate Genetic_Map \n0 %f 0 \n%d %f %f" %(cr, L, cr, cm)
    os.system("echo '" + script + "' >"  + output[0])

# --------------- haps to anc/mut ----------------

inf_trees = sim_trees
ANCMUT = ['anc','mut'] 

rule infer_trees_all:
  input:
    expand(inf_trees, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, q=qs, L=Ls, r=rs, f=fs, n=ns, k=ks, U=Us, END=ANCMUT)

rule infer_trees:
  input:
    expand(sim_hapsample, END=HAPSAMPLE, allow_missing=True),
    sim_map
  output:
    expand(inf_trees, END=ANCMUT, allow_missing=True)
  params:
    prefix=inf_trees.replace(DATADIR,'').replace('.{END}',''),
#    twoNe=lambda wildcards: 2*int(wildcards.K)*4/(2+4*int(wildcards.B)-3)
    twoNe = lambda wildcards: 2*int(wildcards.K) #for WF test
  shell:
    '''
    module load gcc/8.3.0 #needed for relate
    {RELATE}/bin/Relate \
      --mode All \
      -m {wildcards.U} \
      -N {params.twoNe} \
      --haps {input[0]} \
      --sample {input[1]} \
      --map {input[2]} \
      --seed 1 \
      -o {params.prefix}
    mv {params.prefix}.* {DATADIR} 
    '''

# ------------- poplabels ------------

poplabels = DATADIR + 'sim_{k}k.poplabels'

rule poplabels_all:
  input:
    expand(poplabels, k=ks)

rule poplabels:
  input:
  output:
    poplabels
  shell:
    '''
    echo "sample population group sex" > {output[0]}
    for i in {{1..{wildcards.k}}}; do echo "$i 1 1 NA" >> {output[0]}; done
    '''
# ------------- coalescence rates ---------

coal_rates = inf_trees.replace('.{END}','_popsize.{END}')
ANCMUTCOAL = ['anc','mut','coal']

rule coal_rates_all:
  input:
    expand(coal_rates, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, q=qs, L=Ls, r=rs, f=fs, n=ns, k=ks, U=Us, END=ANCMUTCOAL)

rule coal_rates:
  input:
    expand(inf_trees, END=ANCMUT, allow_missing=True),
    poplabels
  output:
    expand(coal_rates, END=ANCMUTCOAL, allow_missing=True)
  params:
    prefix=inf_trees.replace('.{END}',''),
    prefix_out=coal_rates.replace('.{END}','')
  shell:
    '''
    #module load intel/2019u4 #needed for r (gcc also has r but doesnt have needed libraries)
    #module load r/4.1.2 #needed for plots
    # i edited the relate script to avoid plotting
    module load gcc/8.3.0 #needed for relate
    {RELATE}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
              -i {params.prefix} \
              -m {wildcards.U} \
              --years_per_gen 1 \
              --poplabels {input[2]} \
              --seed 1 \
              --num_iter 5 \
              --threshold 0 \
              -o {params.prefix_out} \
              --bins 1,5,0.5
    '''

# -------------- plot tree at selected site -------------

plot_tree = coal_rates.replace('.{END}','.pdf')

rule plot_tree_all:
  input:
    expand(plot_tree, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, q=qs, L=Ls, r=rs, f=fs, n=ns, k=ks, U=Us)

rule plot_tree:
  input:
    expand(sim_hapsample, END=HAPSAMPLE, allow_missing=True),
    expand(coal_rates, END=ANCMUT, allow_missing=True),
    poplabels
  output:
    plot_tree
  params:
    prefix_out = plot_tree.replace('.pdf',''),
    L0=lambda wildcards: round(int(wildcards.L)/2)
  shell:
    '''
    module load gcc/8.3.0 #needed for relate
    module load intel/2019u4 #needed for r (gcc also has r but doesnt have needed libraries)
    module load r/4.1.2 #needed for plots
    {RELATE}/scripts/TreeView/TreeViewMutation.sh \
                   --haps {input[0]} \
                   --sample {input[1]} \
                   --anc {input[2]} \
                   --mut {input[3]} \
                   --poplabels {input[4]} \
                   --bp_of_interest {params.L0} \
                   --years_per_gen 1 \
                   -o {params.prefix_out} 
    '''

# -------------- detect selection -------------

selection = coal_rates
LINFREQSELE = ['lin','freq','sele']

rule detect_selection_all:
  input:
    expand(selection, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, q=qs, L=Ls, r=rs, f=fs, n=ns, k=ks, U=Us, END=LINFREQSELE)

rule detect_selection:
  input:
    expand(coal_rates, END=ANCMUTCOAL, allow_missing=True)
  output:
    expand(selection, END=LINFREQSELE, allow_missing=True) 
  params:
    prefix_in=coal_rates.replace('.{END}',''),
    prefix_out=selection.replace('.{END}','')
  shell:
    '''
    module load gcc/8.3.0
    {RELATE}/scripts/DetectSelection/DetectSelection.sh \
                 -i {params.prefix_in} \
                 -o {params.prefix_out} \
                 -m {wildcards.U} \
                 --years_per_gen 1 \
    '''

# ----------- sample branch lengths ----------------

branch_lengths = coal_rates.replace('.{END}','_sub.timeb')

rule sample_branch_lengths_all:
  input:
    expand(branch_lengths, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, q=qs, L=Ls, r=rs, f=fs, n=ns, k=ks, U=Us)

rule sample_branch_lengths:
  input:
    expand(coal_rates, END=ANCMUTCOAL, allow_missing=True)
  output:
    branch_lengths
  params:
    prefix=coal_rates.replace('.{END}',''),
    prefix_out=branch_lengths.replace('.timeb',''),
    L0=lambda wildcards: round(int(wildcards.L)/2)-1 #selected site - 1
  shell:
    '''
    # find SNP close to selected site (equal to or larger)
    S=$(awk -F';' '(NR>1 && $2>{params.L0}){{print $2; exit}}' {input[1]})
    # then sample branch lengths at that site
    module load gcc/8.3.0
    {RELATE}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
                 -i {params.prefix}\
                 -o {params.prefix_out} \
                 -m {wildcards.U} \
                 --coal {input[2]} \
                 --format b \
                 --num_samples 20 \
                 --first_bp $S \
                 --last_bp $S \
                 --seed 1 
    '''

# ------------- selection coefficient and allele frequency inference --------------

allele_freq = coal_rates.replace('.{END}', '_clues.{END}')
CLUESENDS=['timeBins', 'epochs.npy','freqs.npy','post.npy','png']

rule clues_all:
  input:
    expand(allele_freq, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, q=qs, L=Ls, r=rs, f=fs, n=ns, k=ks, U=Us, END=CLUESENDS)

rule clues:
  input:
    expand(coal_rates, END=['coal'], allow_missing=True),
    branch_lengths
  output:
    expand(allele_freq, END=CLUESENDS, allow_missing=True)
  params:
    prefix_in=branch_lengths.replace('.timeb',''),
    prefix_out=allele_freq.replace('.{END}','')
  shell:
    '''
    for i in {{0,1000}}; do echo "$i" >> {output[0]}; done
    cd {CLUES} #needs to be run in clues directory
    python inference.py \
      --coal ../../{input[0]} \
      --times ../../{params.prefix_in} \
      --out ../../{params.prefix_out} 
    #  --sMax 1 
    #  --tCutoff {t}
    # --timeBins ../../{output[0]}
    # --popFreq {wildcards.f} 
    python plot_traj.py --ext {CLUESENDS[4]} ../../{params.prefix_out} ../../{params.prefix_out} 
    '''

# --------------- true coalescence times ---------------

coal_times = sim_trees
ANCDERTXT = ['anc.txt','der.txt']

rule coal_times_all:
  input:
    expand(coal_times, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, q=qs, L=Ls, r=rs, f=fs, n=ns, k=ks, U=Us, END=ANCDERTXT)

rule coal_times:
  input:
    expand(sim_trees, END=['trees'], allow_missing=True)
  output:
    expand(coal_times, END=ANCDERTXT, allow_missing=True)
  params:
    L0=lambda wildcards: round(int(wildcards.L)/2)
  run:
    ts = tskit.load(input[0]) #load true tree sequence
    tree = ts.at(params.L0) #get tree at selected site
    mutnode = [i.mutations[0].node for i in ts.sites() if i.position==params.L0][0] #get node where mutation appeared (note this assumes a single mutation (update to allow for multiple))
    muttime = [i.mutations[0].time for i in ts.sites() if i.position==params.L0][0] #get time mutation appeared
    dertimes = [] #coalescence times for derived subpop
    anctimes = [] #coalescence times for ancestral subpop
    for node in tree.nodes(): #iterate through all nodes
      if tree.num_children(node) == 2: #check that this is a coalescence node in this tree
        time = ts.node(node).time #time of coalescence
        if time < muttime: # only go back to time of mutation (everything beyond that is neutral)
          if tree.is_descendant(node, mutnode): #if it is below the mutation
            dertimes.append(time) #add coalescence time
          else: # if not below the mutation
            anctimes.append(time) #add coalescence time
    np.savetxt(output[0], anctimes) #output
    np.savetxt(output[1], dertimes)

# ------------- selection coefficient and allele frequency inference with true data --------------

allele_freq_true = coal_rates.replace('.{END}', '_clues_true.{END}')

rule clues_true_all:
  input:
    expand(allele_freq_true, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, q=qs, L=Ls, r=rs, f=fs, n=ns, k=ks, U=Us, END=CLUESENDS)

rule clues_true:
  input:
    expand(coal_rates, END=['coal'], allow_missing=True),
    expand(coal_times, END=ANCDERTXT, allow_missing=True)
  output:
    expand(allele_freq_true, END=CLUESENDS, allow_missing=True)
  params:
    prefix_in=coal_times.replace('.{END}',''),
    prefix_out=allele_freq_true.replace('.{END}','')
  shell:
    '''
    for i in {{0,1000}}; do echo "$i" >> {output[0]}; done
    cd {CLUES} #needs to be run in clues directory
    python inference.py \
      --coal ../../{input[0]} \
      --times ../../{params.prefix_in} \
      --out ../../{params.prefix_out}
    #  --popFreq {wildcards.f}
    #  --sMax 1 #this didnt help 
    #  --tCutoff {t}
    # --timeBins ../../{output[0]}
    python plot_traj.py --ext {CLUESENDS[4]} ../../{params.prefix_out} ../../{params.prefix_out} 
    '''

# -------------- argweaver -----------------

# first we convert the vcf to argweaver's sites format (must do this to include phased samples, see chapter 10 in https://library.oapen.org/handle/20.500.12657/23339)
# we do this with a function used in https://www.biorxiv.org/content/10.1101/2021.11.15.468686v4.abstract
rule vcf_to_sites_all:
  input:
    expand(sim_trees, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, q=qs, L=Ls, r=rs, f=fs, n=ns, k=ks, U=Us, END=['sites'])

rule vcf_to_sites:
  input:
    expand(sim_trees, END=['vcf'], allow_missing=True)
  output:
    expand(sim_trees, END=['sites'], allow_missing=True)
  shell:
    '''
    #git clone https://github.com/deboraycb/ARGsims.git programs/ARGsims
    programs/ARGsims/scripts/argweaver/2_vcf2sites.py \
      {input} \
      {wildcards.k} \
      {wildcards.L} \
      {DATADIR}
    ''' 

# now we sample args as explained in https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008384
args = sim_trees.replace('.{END}','_argweaver/out.{END}')
niters = 3000 #number of MCMC iterations (save every 10th by default)
ARGS = [str(i) + '.smc.gz' for i in range(0,niters+1,10)]

# let's use the time discretization from the CLUES paper (to make sure we can recreate their results)
# wget https://github.com/standard-aaron/clues-v0/blob/master/misc/N_10000_timesfile.txt
timesfile = DATADIR + 'N_10000_timesfile.txt'

rule argweaver_all:
  input:
    expand(args, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, q=qs, L=Ls, r=rs, f=fs, n=ns, k=ks, U=Us, END=ARGS)

rule argweaver:
  input:
    expand(sim_trees, END=['sites'], allow_missing=True),
    timesfile
  output:
    expand(args, END=ARGS, allow_missing=True)
  params:
    prefix = args.replace('.{END}',''),
#    Ne=lambda wildcards: int(wildcards.K)*4/(2+4*int(wildcards.B)-3)
    Ne = lambda wildcards: int(wildcards.K) #for WF test
  shell:
    '''
    $HOME/local/bin/arg-sample \
      -s {input[0]} \
      -o {params.prefix} \
      -N {params.Ne} \
      -m {wildcards.U} \
      -r {wildcards.r} \
      --overwrite \
      --times-file {input[1]} \
      -c 25 \
      -n {niters} \
      --resample-window 40000 \
      --resample-window-iters 8 \
      --infsites
    '''
# above takes 9m

# and now we follow https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008384 to extract the newick trees at the selected site
arg_trees_ends = ['bed.gz','trees']
rule argweaver_trees_all:
  input:
    expand(args, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, q=qs, L=Ls, r=rs, f=fs, n=ns, k=ks, U=Us, END=arg_trees_ends)

rule argweaver_trees:
  input:
    timesfile,
    expand(args, END=ARGS, allow_missing=True)
  output:
    expand(args, END=arg_trees_ends, allow_missing=True)
  params:
    prefix = args.replace('.{END}',''),
    L0=lambda wildcards: round(int(wildcards.L)/2)
  shell:
    '''
    module load samtools # we need tabix and bgzip
    export PATH=$SCRATCH/projects/tsrescue/programs/bedops/bin:$PATH #this gives us bedops' sort-bed
    export PATH=$HOME/local/bin:$PATH #argweaver tools
    smc2bed-all {params.prefix} #this makes output[0] 
    arg-summarize \
      -a {output[0]} \
      -r 1:{params.L0}-{params.L0} \
      -E > {output[1]}
   '''

# to run CLUES-v0 (discretized, like argweaver) we need to precompute some likelihoods
# if we use a different model (otherthan constant N=1e4 and f=0.75) do this
#rule clues_prelim:
#  input:
#  output:
#  shell:

# and now we run CLUES-v0, again following the CLUES paper supplementary material
rule argweaver_clues_all:
  input:
    expand(args, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, q=qs, L=Ls, r=rs, f=fs, n=ns, k=ks, U=Us, END=['clues.h5'])

rule argweaver_clues:
  input:
    expand(args, END=['trees'], allow_missing=True),
    expand(sim_trees, END=['sites'], allow_missing=True)
  output:
    expand(args, END=['clues.h5'], allow_missing=True)
  params:
    prefix = args.replace('{END}','clues'),
    L0=lambda wildcards: round(int(wildcards.L)/2)
  shell:
    '''
    # git clone https://github.com/standard-aaron/clues-v0.git programs/clues-v0
    cd programs/clues-v0
    python \clues.py \
      ../../{input[0]} \
      example.f_75.hdf5 \
      ../../{input[1]} \
      {wildcards.f} \
      --output ../../{params.prefix} \
      --posn {params.L0} \
      --noAncientHap \
      --thin 10 \
      --burnin 100 \
    '''
 
