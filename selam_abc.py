#!/usr/bin/env python3

''' Description: ABC simulations of admixture between populations with
                 a mitonuclear incompatibility (mni) - i.e.
                 one nuclear locus and mtDNA). Used to estimate selection 
                 coef against MNI (both directions), and estimate positive
                 selection coef for the foreign mito that would allow for
                 introgression.

    Grammar notes:
      AAM = minor (foreign) parent nuclear (AA), minor parent mitochondria (M)
          = 111
      aam = major (native) parent nuclear (aa), major parent mitochondria (m)
          = 000

    Outputs:
      selam-log.err:     lists parameter combos that resulted in faulty 
                         SELAM run i.e. for which the sim was skipped
      accepted_sims.tsv: accepted simulations, with relevant parameters

    Why r u using py2? tl;dr issues I am too lazy 2 resolve 

    usage:  ./selam_abc.py > selam_abc.out

    #XXX shooting for 1000 accepted sims, go for 500000 sims

    cyp X-2020
'''

## Want to debug? Setting DEBUG = True is useful
DEBUG = False 

import sys, subprocess
import numpy as np
import pandas as pd
# custom module, calc_selam_freq.py
from calc_selam_freq import calc_freqs, output_freqs

#### XXX HYPERPARAMETERS XXX ####

numsims              = 1000  # number of simulations to run (start w 10000)
                              # depends on acceptance rate
nuc_loc              = 0.2    # location of nuclear locus on chrom (out of 1)
foreign_infreq_mean  = 0.35   # foreign input freq for initial subpop
foreign_infreq_sd    = 0.02   # foreign input freq std dev
sel_group            = "A"    # XXX sex to select on ("A": all, 
                              # "F": females, "M": males)
tolerance            = 0.05   # threshold for accepting sim data
                              # ie how close sims should be to observed

## Choose whether to vary the following variables
#  1: yes, 0: no
vary_n_gen           = 1      # if 0, num_gen  = 50
vary_pop_size        = 1      # if 0, pop_size = 1000
vary_foreign_inf     = 1      # if 0, foreign_infreq = 0.35
vary_fit_mito        = 0      # vary positive sel coef for foreign mito?

## Define observed incompatibility frequencies
obs_aaM              = 0.01   # observed aaM: major nuclear / minor mito, 2/193
obs_AAm              = 0.005  # observed AAm: minor nuclear / major mito, <1/193 

## Bounds of varied parameters (chosen randomly)
n_gen_min            = 20 # 20
n_gen_max            = 80 # 80
pop_size_min         = 200 # 200
pop_size_max         = 2000 # 2000
incompat1_min        = 0 # 0 sel coef against native nuclear, foreign mito (AaM, aaM)
incompat1_max        = 1 # 1
incompat2_min        = 0 # 0 sel coef against foreign nuclear, native mito (Aam, AAm)
incompat2_max        = 1 # 1
fit_mito_min         = 0 # 0
fit_mito_max         = 0.05 # 0.05

## Auto varies dominance coefficients, h1_m & h2_m, 0 to 1

#################################

# essential files
demography_file   = "selam-demography.txt" 
selection_file    = "selam-selection.txt"
output_file       = "selam-outputfile.txt"

# dictionary for storing parameters per run
store_params      = {}

print("\nTotal simulations to run: " + str(numsims) + "\n")
for sim in range(numsims):

  if DEBUG: print("\n~~~ Simulation #" + str(sim) + " ~~~")

  ### draw simulation parameters ###

  ### DEMOGRAPHY.TXT
  ## Choose random population size (200 to 10000)
  if vary_pop_size: 
    pop_size = np.random.randint(pop_size_min, pop_size_max)
    pop_size = int((pop_size / 2) * 2) # to get even random pop size
  else: pop_size = 1000
  ## Choose frequency of foreign input to initial subpop (AAM) 
  # random sample from gaussian with set mean and std dev 
  random_admix   = np.random.normal(foreign_infreq_mean, \
                   foreign_infreq_sd,1)  
  if vary_foreign_inf: foreign_infreq = round(random_admix[0],2) 
  else: foreign_infreq = 0.35 
  native_infreq  = 1 - foreign_infreq
  # write to demography file
  with open(demography_file,'w') as df:
    # header
    df.write("Pop0\tPop1\tSex\t0\t1\n")
    # specify subpopulation size at gen0 (initial) and gen1->N
    df.write("0\t0\tA\t" + str(pop_size) + "\t" + str(pop_size) + "\n")
    # subpop is initially foreign_infreq  major ancestral pop a0 (aam), 
    # no additional input from a0 after gen0
    df.write("0\ta0\tA\t" + str(native_infreq) + "\t0\n")
    # subpop is initially native_infreq minor ancestral pop a1 (AAM)
    # no additional input from a1 after gen0
    df.write("0\ta1\tA\t" + str(foreign_infreq) + "\t0\n")


  ### OUTPUT.TXT
  ## Choose random number of generations (CHAF)
  if vary_n_gen: num_gen = np.random.randint(n_gen_min, n_gen_max)
  else: num_gen = 50 
  out_num_females = int(pop_size / 2) # output all females (1/2 pop) 
  out_num_males   = out_num_females   # output all males (1/2 pop)

  results_file = "sim." + str(sim) + ".selam-mitonuc.out"
  updated_of_line = str(num_gen) + "\t0\t" + str(out_num_females) + \
                    "\t" + str(out_num_males) + "\t" + results_file
  # write num gens to SELAM output parameter file
  with open(output_file,'w') as of:
    of.write(updated_of_line)


  ### SELECTION.TXT
  ## Set selection on mitochondria
  if vary_fit_mito: fit_mito = np.random.uniform(fit_mito_min,fit_mito_max)
  else: fit_mito = 0

  ## Choose random sel coef against native nuclear, foreign mito (AaM, aaM)
  incompat1 = np.random.uniform(incompat1_min, incompat1_max) 
  # Choose random sel coef against foreign nuclear, native mito (AAm, Aam)
  incompat2 = np.random.uniform(incompat2_min, incompat2_max) 


  # Choose random dominance coefficients
  # incompatibility seems to be  recessive 
  # ie. many indvs het at nuc (Aa), with mito from either sp
  #     probably should be close to zero
  h1_m = np.random.uniform(0,1)  #A[aM]
  h2_m = np.random.uniform(0,1)  #[A]a[m]

  ## Establish genotype fitnesses
  fit_AAM = (1 + fit_mito)
  fit_AAm = (1 - incompat2) 
  fit_AaM = (1 + fit_mito)*(1 - (incompat1*h1_m))
  fit_Aam = (1 - (incompat2*h2_m)) 
  fit_aaM = (1 + fit_mito)*(1 - incompat1)
  fit_aam = 1 
  # list fitnesses in the following order:
  # aam, Aam, AAm, aaM, AaM, AAM (000, 100, 001, 101, 111)
  fitnesses = [fit_aam, fit_Aam, fit_AAm, fit_aaM, fit_AaM, fit_AAM]
  # update selection.txt file
  # structure: N A 0 0.2 fit_aam fit_Aam fit_AAm fit_aaM fit_AaM fit_AAM 
  # XXX 0.2 is pos on nuc chromosome in cM --> we could draw a random
  # interval between 0 and 1
  # if you did this, you'd have to find the right position 
  #   should maybe limit position space to front or end of chrom
  str_fitnesses = '\t'.join(str(x) for x in fitnesses) 

  if DEBUG: 
    print("Results in " + results_file)
    print("Fitness of [aam, Aam, AAm, aaM, AaM, AAM]:\n" + str_fitnesses) 
  
  updated_sel_info = "N\t" + sel_group + "\t0\t" + str(nuc_loc) + \
                     "\t" + str_fitnesses
  with open(selection_file, 'w') as sf:
    sf.write(updated_sel_info)


  print("SELAM run#: " + str(sim))
  ## run SELAM
  subprocess.call(["SELAM", "--seed", "1234", "-d", demography_file, "-o", output_file, \
                   "-c", "2", "1", "0", "-s", selection_file])
  print("End SELAM run#: " + str(sim))

  ## calculate and output genoytpe frequencies for each generation file
  # header and simulation info to output 
  extra_header = ["sim#","pop_size","num_gen","initial_admix","incompat1_aM",
                  "incompat2_Am","mito_fit","h1_AaM","h1_Aam"]
  extra_out_info = [str(sim),str(pop_size),str(num_gen),str(foreign_infreq), \
                    str(incompat1),str(incompat2),str(fit_mito), str(h1_m), \
                    str(h2_m)]

  # check that results_file exists (if not, skip this sim)
  try:
    with open(results_file,'r') as f:
      pass
  except IOError as e:
    print("WARNING: " + results_file + " wasn't generated, skipping sim# " + str(sim))
    with open("./selam-log.err", "a") as ef:
      if ef.tell() == 0: ef.write("\t".join(extra_header) + "\n")
      ef.write("\t".join(extra_out_info) + "\n")
    continue

  # calculate genotype frequencies (function from calc_selam_freq.py)
  F_dict, M_dict, total_dict  = calc_freqs(results_file, nuc_loc) 

  if DEBUG:
    print("#gens=" + str(num_gen) + ", #indvs=" + \
          str(pop_size) + "\nFrequencies: " + str(total_dict)) 
  #print(results_file + "\nFemales (" + str(out_num_females) + \
  #      "): " + str(F_dict) + "\nMales (" + str(out_num_males) + \
  #      "): " + str(M_dict))

  # write out frequencies in nice format to output file
  outfile_prefix = "all-mni-sim.selam-freqs"

  # output female and male freqs in separate files
  output_freqs(results_file,outfile_prefix,F_dict,dict2=M_dict,setting="sep", \
               extra_header=extra_header,extra_output_info=extra_out_info)
  # output total freqs in another file
  output_freqs(results_file,outfile_prefix,total_dict,setting="all", \
               extra_header=extra_header,extra_output_info=extra_out_info)

  aaM_freq = total_dict["001"] # freq of homozygous incompat1
  AAm_freq = total_dict["110"] # freq of homozygous incompat2
  store_params[sim] = [pop_size,num_gen,foreign_infreq,incompat1,incompat2,fit_mito, \
                       h1_m,h2_m,aaM_freq,AAm_freq]
  if DEBUG:
    print(["pop_size","gens","admix","i1","i2","fit_mito","h1","h2","aaM","AAm"])
    print(store_params[sim])

# first check that store_params is not empty
# if it is, exit program
if not store_params: 
  print("WARNING: No simulations were completed, exiting...")
  exit()

# choose simulations that satisfy threshold
# sims are kept if: freq of AAm < observed & freq of aaM is within tolerated range of observed
sims_df = pd.DataFrame.from_dict(store_params, orient='index')
sims_df.columns = ["pop_size","gens","admix","i1","i2","fit_mito","h1","h2","aaM","AAm"]
accepted_sims = sims_df.query('AAm < @obs_AAm & aaM > (1-@tolerance)*@obs_aaM \
                              & aaM < (1+@tolerance)*@obs_aaM')

if accepted_sims.empty:
  print("\nNo simulations yielded observed frequencies, no outfile generated.")
else:
  # write results to outfile
  accepted_sims.to_csv("./accepted_sims.tsv", sep="\t", index=True)

  # STDOUT results
  print("\nSome simulations yielded observed frequencies!")
  pd.set_option('display.max_columns', None)
  pd.set_option('display.width', None)
  print(accepted_sims)
