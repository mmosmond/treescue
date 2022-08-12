# ------------ python modules ----------

import pyslim, msprime, tskit
import numpy as np
import os

# ------------- paths ------------------

DATADIR = 'data/' #where to put outputs
PROGRAMDIR = 'programs/'
SLIM = PROGRAMDIR + 'SLiM_3.7.1/slim' #command to run SLiM
RELATE = PROGRAMDIR + 'relate'

#  ----------- parameters -----------------

# for simulations 
Ns = [int(1e4)] #intial population size
Ks = [int(1e4)] #carrying capacity
ds = [0.05] #decline rate of ancestral homozygote: absolute fitness 1-d
ss = [0.13] #selection coefficient of derived allele: absolute fitness of derived homozygote (1-d)(1+s)
hs = [0.5] #dominance coefficient of derived allele: absolute fitness of the heterozygote (1-d)(1+hs)
Bs = [2] #number of offspring per parent (note we divide viability by this so that it does not affect the expected number of offspring (absolute fitness), but it will affect the variance in the number of offspring (Ne))
us = [1e-5] #mutation rate at selected site (will add neutral mutations to other sites after with msprime)
Ls = [int(2e7)] #number of sites minus 1 (since SLiM 0 indexes sites; will put selected site in middle) 
rs = [2e-8] #recombination rate per site
ts = [int(1e3)] #max number of generations (to prevent simulations from running forever) 
ns = range(1) #replicate
# for vcfs
ks = [100] #number of individuals to sample
Nes = [int(1e4)] #effective population size during recapitation
Us = [7e-9] #mutation rate at neutral sites

# ------- simulations -------

sim_trees = DATADIR + "sim_{N}N_{K}K_{d}d_{s}s_{h}h_{B}B_{u}u_{L}L_{r}r_{t}t_{n}n.trees"

rule simulate_all:
  input:
    expand(sim_trees, N=Ns, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, L=Ls, r=rs, t=ts, n=ns)

rule simulate:
  input:
    "scripts/sim.slim"
  output:
    sim_trees
  shell:
    """
    mkdir -p {DATADIR}
    {SLIM} \
      -d N={wildcards.N} \
      -d K={wildcards.K} \
      -d d={wildcards.d} \
      -d s={wildcards.s} \
      -d h={wildcards.h} \
      -d B={wildcards.B} \
      -d u={wildcards.u} \
      -d L={wildcards.L} \
      -d r={wildcards.r} \
      -d t={wildcards.t} \
      -d "output='{output}'" \
      {input}
    """

# ---------- trees to vcfs --------

sim_vcf = sim_trees + '_{k}k_{Ne}Ne_{U}U.vcf'

rule vcf_all:
  input:
    expand(sim_vcf, N=Ns, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, L=Ls, r=rs, t=ts, n=ns, k=ks, Ne=Nes, U=Us)
 
rule vcf:
  input:
    sim_trees
  output:
    sim_vcf
  run:
    ts = tskit.load(input[0]) #load tree sequence
    inds = np.random.choice(range(ts.num_individuals), int(wildcards.k), replace=False) #sample k individuals at random
    ind_nodes = [] #empty vector for nodes
    for ind in inds: #for each individual
      for node in ts.individual(ind).nodes:
        ind_nodes.append(node) #append each genome
    ts = ts.simplify(ind_nodes, keep_input_roots=True) #simplify down to sample while keeping roots for recapitation
    ts = pyslim.update(ts) #using beta pyslim for slim 4.0 but trees made with slim 3.7.1, so update
    ts = pyslim.recapitate(ts, recombination_rate=float(wildcards.r), ancestral_Ne=float(wildcards.Ne)) #recapitate trees
    #ts = msprime.sim_mutations(ts, rate=float(wildcards.U), model=msprime.SLiMMutationModel(type=0), keep=True) #layer on neutral mutations while keeping the selected mutation
    #ts = msprime.sim_mutations(ts, rate=float(wildcards.U), model=msprime.BinaryMutationModel(), keep=True) #layer on neutral mutations while keeping the selected mutation
    #ts = pyslim.generate_nucleotides(ts) #https://github.com/tskit-dev/pyslim/pull/174
    #ts = pyslim.convert_alleles(ts)
    ts = msprime.sim_mutations(ts, rate=float(wildcards.U), model=msprime.JC69()) #layer on neutral mutations (drop fixed rescue mutation because not biallelic)
    with open(output[0], "w") as vcf_file:
      ts.write_vcf(vcf_file, individuals = range(int(wildcards.k))) #write out as VCF with arbitrary ID names

# ------------ vcfs to haps -------------

sim_hapsample = sim_vcf + '.{END}'
HAPSAMPLE = ['haps','sample']

rule haps_all:
  input:
    expand(sim_hapsample, N=Ns, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, L=Ls, r=rs, t=ts, n=ns, k=ks, Ne=Nes, U=Us, END=HAPSAMPLE)

rule haps:
  input:
    sim_vcf
  output:
    expand(sim_hapsample, END=HAPSAMPLE, allow_missing=True)
  params:
    prefix=sim_vcf.replace('.vcf','')
  shell:
    '''
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

inf_trees = sim_vcf.replace('vcf','{END}')
ANCMUT = ['anc','mut'] 

rule infer_trees_all:
  input:
    expand(inf_trees, N=Ns, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, L=Ls, r=rs, t=ts, n=ns, k=ks, Ne=Nes, U=Us, END=ANCMUT)

rule infer_trees:
  input:
    expand(sim_hapsample, END=HAPSAMPLE, allow_missing=True),
    sim_map
  output:
    expand(inf_trees, END=ANCMUT, allow_missing=True)
  params:
    prefix=inf_trees.replace(DATADIR,'').replace('.{END}',''),
    twoNe=lambda wildcards: 2*int(wildcards.Ne)
  shell:
    '''
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
    expand(coal_rates, N=Ns, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, L=Ls, r=rs, t=ts, n=ns, k=ks, Ne=Nes, U=Us, END=ANCMUTCOAL)

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
    module load intel/2019u4
    module load r/4.1.2
    {RELATE}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
              -i {params.prefix} \
              -m {wildcards.U} \
              --years_per_gen 1 \
              --poplabels {input[2]} \
              --seed 1 \
              --num_iter 5 \
              -o {params.prefix_out}
    '''

# ----------- branch lengths ----------------

branch_lengths = coal_rates.replace('.{END}','_sub.{END}')
TIMEB = ['timeb'] 

rule sample_branch_lengths_all:
  input:
    expand(branch_lengths, N=Ns, K=Ks, d=ds, s=ss, h=hs, B=Bs, u=us, L=Ls, r=rs, t=ts, n=ns, k=ks, Ne=Nes, U=Us, END=TIMEB)

rule sample_branch_lengths:
  input:
    expand(coal_rates, END=ANCMUTCOAL, allow_missing=True)
  output:
    expand(branch_lengths, END=TIMEB, allow_missing=True)
  params:
    prefix=coal_rates.replace('.{END}',''),
    prefix_out=branch_lengths.replace('.{END}',''),
    L0=lambda wildcards: round(int(wildcards.L)/2)
  shell:
    '''
    # find SNP close to selected site (closest that is larger)
    S=$(awk -F';' '(NR>1 && $2>{params.L0}){{print $2; exit}}' {input[1]})
    {RELATE}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
                 -i {params.prefix}\
                 -o {params.prefix_out} \
                 -m {wildcards.U} \
                 --coal {input[2]} \
                 --format b \
                 --num_samples 5 \
                 --first_bp $S \
                 --last_bp $S \
                 --seed 1 
    '''
# this fails since we dont have a SNP at the selected site, L0. so we need to either end the simulation before fixation or track a neighboring SNP.
