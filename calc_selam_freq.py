#!/usr/bin/env python3
''' Description: extract mitonuc genotype frequencies
                 including aaM (001) and AAm (110)
    
    Usage: ./calc_selam_freq.py infile nuc_loc 

    infile=  SELAM results file for interaction between
             one nuc chrom pair and one mtDNA, with format:
           [0] gen, [1] subpop, [2] [indv:0=fem,1=male],
           [3] indv#, [4] chr#, [5] [chrom:0=mother,1=father],
           [6] tract_ancestry, [7] tract_start, [8] tract_end
    nuc_loc= location of nuclear interacting locus on chrom
             (ex. 0.2)
'''

import pandas as pd

def calc_freqs(infile, nuc_loc):
  # read in results file
  res=pd.read_table(infile,header=None,engine='python')
  res.columns = ['gen','subpop','sex','indv','chr_num', \
                 'chr_parent','tract_an','tract_start','tract_end']

  ## collect ancestry tracts containing the nuc locus
  nuc=res[(res['tract_start']<=nuc_loc) & (res['tract_end']>=nuc_loc)]
  # make dictionary of { individual_#: ancestry@locus }
  # note this collects four vals per key (2 for female, 2 for male)
  nuc_dict  = nuc.groupby('indv')['tract_an'].apply(lambda x: x.tolist()) \
              .to_dict()

  ## collect mito ancestry 
  # selam results file mtdna lines don't follow nuc chrom output
  # format, so the names we use here are wonky
  # 'chr_parent' is filled with the 'mtdna' value
  # 'chr_num' is the individual number
  # 'tract_an' is still the mtdna ancestry
  mito = res[res['chr_parent'] == 'mtdna']
  mito_dict = mito.groupby('chr_num')['tract_an'].apply(lambda x: x.tolist()) \
              .to_dict()

  ## now get frequencies of each nuc-mito combo
  # make variables for frequency counts, set them to 0
  aam, aaM, Aam, AaM, AAm, AAM = 0, 0, 0, 0, 0, 0
  # counts for mitonuc genotypes:
  # key is (nuc_ancestry, mtdna_ancestry), 
  # value is count of individuals with this mitonuc genotype
  # nuc_ancestry encoding: 0 = 00, 1 = 01/10, 2 = 11 
  # {(0,0): aam, (0,1): aaM, (1,0): Aam, (1,1): AaM,
  #  (2,0): AAm, (2,1): AAM }
  F_mitonuc_counts = {(0,0): 0, (0,1): 0, (1,0): 0, (1,1): 0,
                      (2,0): 0, (2,1): 0} # females
  M_mitonuc_counts = {(0,0): 0, (0,1): 0, (1,0): 0, (1,1): 0,
                      (2,0): 0, (2,1): 0} # males

  for i,nuc_list in nuc_dict.items():
    # first two nuc ancestries are for female nuc
    # first mito ancestry is for female mtdna
    Fn = nuc_list[0]+nuc_list[1] # add for nuc ancestry coding
    Fm = mito_dict[i][0]         # mito coding is the same 
    # add to total count in mitonuc_counts dict
    F_mitonuc_counts[(Fn,Fm)] += 1

    # second two ancestries are for male nuc
    # second mito ancestry is for male mtdna
    Mn = nuc_list[2]+nuc_list[3]
    Mm = mito_dict[i][1]
    M_mitonuc_counts[(Mn,Mm)] += 1

  # convert counts to frequencies
  total_Findvs    = sum(F_mitonuc_counts.values(), 0.0)
  print("total_Findvs")
  print(total_Findvs)
  print(F_mitonuc_counts)
  F_mitonuc_freqs = {k: v / total_Findvs 
                     for k, v in F_mitonuc_counts.items()}
  print(F_mitonuc_freqs)
  total_Mindvs    = sum(M_mitonuc_counts.values(), 0.0)
  print("total_Mindvs")
  print(total_Mindvs)
  print(M_mitonuc_counts)
  M_mitonuc_freqs = {k: v / total_Mindvs
                     for k, v in M_mitonuc_counts.items()}
  print(M_mitonuc_freqs)

  ## once each individual (female and male) genotype has been
  ## calculated, return frequencies
  # make encoding dictionary for mitonuc geno names
  geno_names = {(0,0): "000", (0,1): "001", 
                (1,0): "100", (1,1): "101",
                (2,0): "110", (2,1): "111" }

  F_output_dict, M_output_dict, total_output_dict = {}, {}, {}
  for g in geno_names.keys():
    F_output_dict[geno_names[g]] = F_mitonuc_freqs[g]
    M_output_dict[geno_names[g]] = M_mitonuc_freqs[g]
    total_output_dict[geno_names[g]] = \
            (F_mitonuc_counts[g] + M_mitonuc_counts[g])/ \
            (total_Findvs + total_Mindvs)

  ## return both sex and total output dicts
  return F_output_dict, M_output_dict, total_output_dict
  # uncomment to return counts dicts too, to debug
  #return F_output_dict, M_output_dict, F_mitonuc_counts, M_mitonuc_counts

def output_freqs(results_file_name, outfile_prefix, dict1, dict2=None, \
                 setting="sep", extra_header=[], extra_output_info=[]):
  ''' Outputs genotype frequency information from a SELAM run
      Arguments:  setting="sep" keep sex info separate, "all" to combine 
                  dict1=dictionary with female frequencies, or total if
                        setting="all"
                  dict2=dicitonary with male frequencies
                  outfile_prefix=desired prefix for output file
  '''
  
  header = "results_file\t" + "\t".join(extra_header) + "\t" + \
           "000_freq\t001_freq\t100_freq\t101_freq\t110_freq\t111_freq\tA_freq\tM_freq\n"
  out_line = results_file_name + "\t" + "\t".join(extra_output_info) + "\t"
  freq_order = ["000","001","100","101","110","111"]
  if setting == "all":
    # output one freq file, with female and male freqs combined
    total_dict = dict1
    out_file = outfile_prefix+ "_all.tsv"
    with open(out_file, "a") as of:
      if of.tell() == 0: of.write(header)
      for k in freq_order:
        out_line += str(total_dict[k]) + "\t"
      # calculate allele frequencies of A and M
      # A = 0.5*AaM + 0.5*Aam + AAM + AAm
      # M = AAM + AaM + aaM
      A_freq = 0.5*total_dict["101"] + 0.5*total_dict["100"] + \
               total_dict["111"] + total_dict["110"]
      M_freq = total_dict["111"] + total_dict["101"] + total_dict["001"]
      out_line += str(A_freq) + "\t" + str(M_freq)
      of.write(out_line + "\n")
  else:
    # output two freq files, one per sex
    F_dict = dict1
    M_dict = dict2
    # write female freqs to female file
    out_file = outfile_prefix + "_females.tsv"
    with open(out_file, "a") as of:
      if of.tell() == 0: of.write(header)
      out_line_f = out_line
      for k in freq_order:
        out_line_f += str(F_dict[k]) + "\t"
      # calculate foreign allele freqs
      A_freq = 0.5*F_dict["101"] + 0.5*F_dict["100"] + \
               F_dict["111"] + F_dict["110"]
      M_freq = F_dict["111"] + F_dict["101"] + F_dict["001"]
      out_line_f += str(A_freq) + "\t" + str(M_freq)
      of.write(out_line_f + "\n")

    # write male freqs to male file
    out_file = outfile_prefix + "_males.tsv"
    with open(out_file, "a") as of:
      if of.tell() == 0: of.write(header)
      out_line_m = out_line
      for k in freq_order:
        out_line_m += str(M_dict[k]) + "\t"
      # calculate foreign allele freqs
      A_freq = 0.5*M_dict["101"] + 0.5*M_dict["100"] + \
               M_dict["111"] + M_dict["110"]
      M_freq = M_dict["111"] + M_dict["101"] + M_dict["001"]
      out_line_m += str(A_freq) + "\t" + str(M_freq)
      of.write(out_line_m + "\n")


# execute
if __name__ == "__main__":
  import sys
  results_file = sys.argv[1]
  nuc_loc      = float(sys.argv[2])
  F_freqs, M_freqs = calc_freqs(results_file,nuc_loc)
  print(F_freqs)
  print(M_freqs)
  output_freqs(F_freqs,M_freqs,results_file,"t3_outfile",extra_header=["gen"],extra_output_info=["0"])
  #F_freqs,M_freqs,F_counts,M_counts = calc_freqs(sys.argv[1], \
  #                                    float(sys.argv[2]))
  #print(F_counts)
  #print(M_counts)
