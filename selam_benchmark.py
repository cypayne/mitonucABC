#!/usr/bin/env python2

''' Description: For benchmarking SELAM

genotypes are population dependent:
  AAM = minor parent nuclear (AA), minor parent mitochondria (M)
      = 111
  aam = major parent nuclear (aa), major parent mitochondria (m)
      = 000

usage:  ./selam_benchmark.py > selam_benchmark.out

'''
import sys, subprocess
import numpy as np
# custom module, calc_selam_freq.py
from calc_selam_freq import calc_freqs, output_freqs


#### XXX HYPERPARAMETERS XXX ####

numsims              = 1    # number of simulations to run (start w 10000)
                            # depends on acceptance rate
foreign_infreq       = 0.35 # foreign input freq for initial subpop
pop_size             = 10000 # subpopulation size
out_num_females      = 5000 # num females to output from subpop
out_num_males        = 5000 # num males to output from subpop
num_gen              = 1000 # number of generations (1000)
gen_step             = 10   # output generation info at each step
nuc_loc              = 0.2  # location of nuclear locus on chrom (out of 1)
sel_group            = "A"  # sex to select on ("A": all, "F": females, "M": males)
incompat1            = 0.05 # AaM, aaM
incompat2            = 0.05 # AAm, Aam
h1_m                 = 0    # dominance coef for AaM
h2_m                 = 0    # dominance coef for Aam
fit_mito             = 0 # sel coef if foreign mito M has advantage

#################################

# essential files
demography_file   = "selam-demography.txt" 
selection_file    = "selam-selection.txt"
output_file       = "selam-outputfile.txt"

# add 1 to num_gen to output last generation
num_gen += 1

# loop through simulations
for k in range(numsims):
  
  print("\n~~~ Simulation #" + str(k) + " ~~~")

  ### draw simulation parameters ###

  ### OUTPUT.TXT
  ## Want to write out every 10 generations / 1000 gens 
  with open(output_file,'w') as of:
    for g in range(0,num_gen,gen_step): 
      results_file = "sim." + str(k) + ".gen." + str(g) + ".selam-MNI.out"
      updated_of_line = str(g) + "\t0\t" + str(out_num_females) + \
                        "\t" + str(out_num_males) + "\t" + results_file + "\n"
      of.write(updated_of_line)

  ### DEMOGRAPHY.TXT
  ## Choose frequency of foreign input to initial subpop (AAM) 
  native_infreq  = 1 - foreign_infreq
  # put into demography file
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

  ### SELECTION.TXT
  ## Establish genotype fitnesses
  fit_AAM = (1 + fit_mito)
  fit_AAm = (1 - incompat2) 
  fit_AaM = (1 + fit_mito)*(1 - (incompat1*h1_m))
  fit_Aam = (1 - (incompat2*h2_m)) 
  fit_aaM = (1 + fit_mito)*(1 - incompat1)
  fit_aam = 1 
  # list fitnesses in the following order:
  #   aam, Aam, AAm, aaM, AaM, AAM
  fitnesses = [fit_aam, fit_Aam, fit_AAm, fit_aaM, fit_AaM, fit_AAM]
  # update selection.txt file
  # format: N A 0 0.2 fit_aam fit_Aam fit_AAm fit_aaM fit_AaM fit_AAM 
  #         where nuclear locus is on chrom 0 at position 0.2
  print("aam, Aam, AAm, aaM, AaM, AAM fitnesses: " + str(fitnesses))
  str_fitnesses = '\t'.join(str(x) for x in fitnesses) 
  updated_sel_info = "N\t" + sel_group + "\t0\t" + str(nuc_loc) + \
                     "\t" + str_fitnesses 
  with open(selection_file, 'w') as sf:
    sf.write(updated_sel_info)

  ## run SELAM
  subprocess.call(["SELAM", "-d", demography_file, "-o", output_file, \
                   "-c", "2", "1", "0", "-s", selection_file])

  ## calculate and output genoytpe frequencies for each generation file
  # extra header info to output with output_freqs
  extra_header = ["gen","initial_admix","incompat1_aM","incompat2_Am",\
                  "mito_fit","h1_AaM","h1_Aam"]
  for g in range(0,num_gen,gen_step):
    extra_out_info = [str(g),str(foreign_infreq),str(incompat1),\
                      str(incompat2),str(fit_mito), str(h1_m),str(h2_m)]
    results_file = "sim." + str(k) + ".gen." + str(g) + ".selam-MNI.out"
    # calculate genotype frequencies
    F_dict, M_dict, total_dict  = calc_freqs(results_file, nuc_loc) 
    # print results to STDOUT
    print(results_file + "\nFemales (" + str(out_num_females) + \
          "): " + str(F_dict) + "\nMales (" + str(out_num_males) + \
          "): " + str(M_dict))
    # write out frequencies in nice format to output files
    outfile_prefix = "sim." + str(k) + ".selam-freqs"
    # output female and male freqs in separate files
    output_freqs(results_file,outfile_prefix,F_dict,dict2=M_dict,setting="sep", \
                 extra_header=extra_header,extra_output_info=extra_out_info)
    # output total freqs in another file
    output_freqs(results_file,outfile_prefix,total_dict,setting="all", \
                 extra_header=extra_header,extra_output_info=extra_out_info)

