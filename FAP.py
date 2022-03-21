from olexFunctions import OlexFunctions
OV = OlexFunctions()

import os
import htmlTools
import olex
import olx
import gui
import shutil
import sys
from return_WL_energy import ret_wl
from time import sleep
import logging


import time
debug = bool(OV.GetParam("olex2.debug", False))

instance_path = OV.DataDir()

try:
  from_outside = False
  p_path = os.path.dirname(os.path.abspath(__file__))
except:
  from_outside = True
  p_path = os.path.dirname(os.path.abspath("__file__"))

l = open(os.sep.join([p_path, 'def.txt'])).readlines()
d = {}
for line in l:
  line = line.strip()
  if not line or line.startswith("#"):
    continue
  d[line.split("=")[0].strip()] = line.split("=")[1].strip()

p_name = d['p_name']
p_htm = d['p_htm']
p_img = eval(d['p_img'])
p_scope = d['p_scope']

OV.SetVar('FAP_plugin_path', p_path)

from PluginTools import PluginTools as PT

class FAP(PT):

  def __init__(self):
    super(FAP, self).__init__()
    self.p_name = p_name
    self.p_path = p_path
    self.p_scope = p_scope
    self.p_htm = p_htm
    self.p_img = p_img
    self.deal_with_phil(operation='read')
    self.print_version_date()
    self.base_path = ""
    self.solution_path = ""
    if not from_outside:
      self.setup_gui()
    OV.registerFunction(self.print_formula,True,"FAP")
    OV.registerFunction(self.setBasePath,True,"FAP")
    OV.registerFunction(self.extract_disps,True,"FAP")
    OV.registerFunction(self.setSolutionPath,True,"FAP")

    # END Generated =======================================
    
  def setBasePath(self):
    out = olex.f('fileOpen("Please choose a file inside the folder where your data is stored", "*", filepath())')
    buffer = out.split("\\")
    buffer = buffer[:-1]
    FAP_instance.base_path = "\\".join(buffer)
    print(f"Your data lies at:\t{buffer}")

  def getBasePath(self):
     return self.base_path

  def setSolutionPath(self):
    out = olex.f('fileOpen("Choose Your solution .ins file", "*", filepath())')
    buffer = out.split("\\")
    FAP_instance.solution_path = "\\".join(buffer)
    print(f"Your solution lies at:\t{buffer}")

  def getSolutionPath(self):
    return self.solution_path

  def print_formula(self):
    use_nosphera2 = OV.GetParam("fap.use_nosphera2")
    elements=OV.GetParam("fap.element_string")
    adjustment = float(OV.GetParam("fap.adjustment_eV"))
    resolution = OV.GetParam("fap.resolution")
    
    base_path = FAP_instance.getBasePath()
    solution_path = FAP_instance.getSolutionPath()

    with open(f"{base_path}/log.txt", "w+") as log:
        log.write(str(elements))
        
    if base_path == "" or solution_path == "" :
        raise ImportError("Please provide a base and solution path!")

    #define the string to refine dispersion values for all elements
    disp_ref_string = f'REM RefineDisp="{elements}"\n'
    
    #Search for all sub dirs in base dir and storing their names in sub_dirs
    out = [x for x in os.walk(base_path)]
    sub_dirs = out[0][1]
    name = None    

    #import the .ins file provided in the solutionpath
    copy_ins = []
    state = False
    with open(fr"{solution_path}", "r") as sol:
        for line in sol:
            if line.startswith("LATT"):
                state = True
            if state:
                copy_ins.append(line)
    
    if use_nosphera2:
        #Define output buffers as empty lists. This has to be redone drastically
        energies = []
        wls = []
        f_primes = []
        f_double_primes = []
        r1 = []
        wr2 = []
        f_primes_henke = []
        f_double_primes_henke = []
        r1_henke = []
        wr2_henke = []
        f_primes_refined = []
        f_double_primes_refined = []
        r1_refined = []
        wr2_refined = []
        f_primes_nosphera2 = []
        f_double_primes_nosphera2 = []
        r1_nosphera2 = []
        wr2_nosphera2 = []
        f_primes_henke_nosphera2 = []
        f_double_primes_henke_nosphera2 = []
        r1_henke_nosphera2 = []
        wr2_henke_nosphera2 = []
        f_primes_refined_nosphera2 = []
        f_double_primes_refined_nosphera2 = []
        r1_refined_nosphera2 = []
        wr2_refined_nosphera2 = []

        with open(f"{base_path}/log.txt", "w+") as log:
            
            
            for dir in sub_dirs:
                if name == None:
                    name = dir.split("_")[0]
                else:
                    if name != dir.split("_")[0] and (name != "tempFAP" or name != "solution"):
                        if name == "tempFAP" or name != "solution":
                            continue
                        raise ValueError("There are folders in your working dir belonging to a different structure!")
                
                energy = int(dir.split("_")[1])
                cor_energy = energy + adjustment
                wl = ret_wl(cor_energy)
                energies.append(cor_energy)
                wls.append(wl)

                print(f"Currently working at:\nName:\t{name}\t\t(Corrected) Energy:\t{cor_energy}")

                with open(fr"{base_path}\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution_refined_nosphera2.ins", "w+") as no2_refined_out:
                    with open(fr"{base_path}\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution_henke_nosphera2.ins", "w+") as no2_out_henke:
                        with open(fr"{base_path}\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution_nosphera2.ins", "w+") as no2_out:
                            with open(fr"{base_path}\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution_refined.ins", "w+") as output_ref:
                                with open(fr"{base_path}\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution_henke.ins", "w+") as output_henke:
                                    with open(fr"{base_path}\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution.ins", "w+") as output:
                                        with open(fr"{base_path}\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}.ins", "r") as input:
                                            curr_in = [line for line in input]
                                            state = True
                                            for line in curr_in:
                                                buffer = []
                                                if state:
                                                    if line.startswith("CELL"):
                                                        buffer = line[:5] + str(round(wl,6))[:7] + line[12:]
                                                    else:
                                                        buffer = line
                                                    if line.startswith("LATT"):
                                                        state = False
                                                        continue
                                                    no2_out.write(buffer)
                                                    no2_out_henke.write(buffer)
                                                    no2_refined_out.write(buffer)
                                                    output.write(buffer)
                                                    output_henke.write(buffer)
                                                    output_ref.write(buffer)
                                                else:
                                                    break
                                            for line in copy_ins:
                                                if line == "REM <olex2.extras>\n":
                                                    output_ref.write(line)
                                                    no2_refined_out.write(line)
                                                    no2_out.write(line)
                                                    no2_out_henke.write(line)
                                                    output.write(line)
                                                    output_henke.write(line)
                                                else:
                                                    no2_out.write(line)
                                                    no2_refined_out.write(line)
                                                    output_ref.write(line)
                                                    output.write(line)
                                                    output_henke.write(line)
                                                    no2_out_henke.write(line)
                global nosphera2_switch
                nosphera2_switch = True
                
            
                log.write(f"Currently working at:\nName:\t{name}\t\t(Corrected) Energy:\t{cor_energy}\with elements:\t{elements}\n")
                    
                for nosphera2 in [False, True]:
                    if nosphera2:
                        ins_path = base_path + fr"\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution_nosphera2.ins"
                        ins_path_henke = base_path + fr"\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution_henke_nosphera2.ins"
                        ins_disp_ref_path = base_path + fr"\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution_refined_nosphera2.ins"

                    else: 
                        ins_path = base_path + fr"\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution.ins"
                        ins_path_henke = base_path + fr"\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution_henke.ins"
                        ins_disp_ref_path = base_path + fr"\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution_refined.ins"
                        
                    #Refinement fpr Sasaki Table / Sasaki NoSpherA2
                    refine_please(ins_path, nosphera2, resolution, henke = False)

                    # if nosphera2_switch & nosphera2:
                    #     template_tsc = base_path + fr"\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution_nosphera2.tsc"
                    #     nosphera2_switch = False
                    #     for dir in sub_dirs:
                    #         shutil.copy2(template_tsc, dir)
                    # else:
                    #     continue

                    if nosphera2:
                        r1e, wr2e = float(OV.GetParam("snum.refinement.last_R1")),float(OV.GetParam("snum.refinement.last_wR2"))
                        r1_nosphera2.append(r1e)
                        wr2_nosphera2.append(wr2e)
                        log.write(f"{energy}\tR1sasnos:\t{str(r1e)}\twR2sasnos:\t{str(wr2e)}\n")
                        
                    else:
                        r1e, wr2e = float(OV.GetParam("snum.refinement.last_R1")),float(OV.GetParam("snum.refinement.last_wR2"))
                        r1.append(r1e)
                        wr2.append(wr2e)
                        log.write(f"{energy}\tR1sas:\t{str(r1e)}\twR2sas:\t{str(wr2e)}\n")

                    with open(ins_path, "r") as out_ins:
                        for line in out_ins:
                            if line.startswith("DISP") and "Mo" in line:
                                l = line.split()
                                if nosphera2:
                                    fp, fdp = float(l[2]),float(l[3])
                                    f_primes_nosphera2.append(fp)
                                    f_double_primes_nosphera2.append(fdp)
                                    log.write(f"{energy}\tfpsasnos:\t{str(fp)}\tfdpsas_nos:\t{str(fdp)}\n")
                                else:
                                    fp, fdp = float(l[2]),float(l[3])
                                    f_primes.append(fp)
                                    f_double_primes.append(fdp)
                                    log.write(f"{energy}\tfpsas:\t{str(fp)}\tfdpsas:\t{str(fdp)}\n")
                    
                    #Refinement for Henke table / Henke NospherA2            
                    refine_please(ins_path_henke, nosphera2, resolution, henke = True)            

                    if nosphera2:
                        r1e, wr2e = float(OV.GetParam("snum.refinement.last_R1")),float(OV.GetParam("snum.refinement.last_wR2"))
                        r1_henke_nosphera2.append(r1e)
                        wr2_henke_nosphera2.append(wr2e)
                        log.write(f"{energy}\tR1hennos:\t{str(r1e)}\twR2hennos:\t{str(wr2e)}\n")

                    else:
                        r1e, wr2e = float(OV.GetParam("snum.refinement.last_R1")),float(OV.GetParam("snum.refinement.last_wR2"))
                        r1_henke.append(r1e)
                        wr2_henke.append(wr2e)
                        log.write(f"{energy}\tR1hen:\t{str(r1e)}\twR2hen:\t{str(wr2e)}\n")

                    with open(ins_path_henke, "r") as out_ins:
                        for line in out_ins:
                            if line.startswith("DISP") and "Mo" in line:
                                l = line.split()
                                if nosphera2:
                                    fp, fdp = float(l[2]),float(l[3])
                                    f_primes_henke_nosphera2.append(fp)
                                    f_double_primes_henke_nosphera2.append(fdp)
                                    log.write(f"{energy}\tfphennos:\t{str(fp)}\tfdphen_nos:\t{str(fdp)}\n")
                                else:
                                    fp, fdp = float(l[2]),float(l[3])
                                    f_primes_henke.append(fp)
                                    f_double_primes_henke.append(fdp)
                                    log.write(f"{energy}\tfphen:\t{str(fp)}\tfdphen:\t{str(fdp)}\n")

                    #Refinement for refined / refined NospherA2

                    refine_please(ins_disp_ref_path, nosphera2, resolution, elements,henke = False)

                    if nosphera2:
                        r1e, wr2e = float(OV.GetParam("snum.refinement.last_R1")),float(OV.GetParam("snum.refinement.last_wR2"))
                        r1_refined_nosphera2.append(r1e)
                        wr2_refined_nosphera2.append(wr2e)
                        log.write(f"{energy}\tR1refnos:\t{str(r1e)}\twR2refnos:\t{str(wr2e)}\n")
                    else:
                        r1e, wr2e = float(OV.GetParam("snum.refinement.last_R1")), float(OV.GetParam("snum.refinement.last_wR2"))
                        r1_refined.append(r1e)
                        wr2_refined.append(wr2e)
                        log.write(f"{energy}\tR1ref:\t{str(r1e)}\twR2ref:\t{str(wr2e)}\n")                  

                    with open(ins_disp_ref_path, "r") as out_ins:           # This searches the INS for the information on Mo DISPS. For other disps, this section yet has to be edited
                        for line in out_ins:
                            if line.startswith("REM  <Mo"):
                                l = line.split()
                                if nosphera2:
                                    fp, fdp = float(l[2].replace('"', '')),float(l[3].replace('">', ''))
                                    f_primes_refined_nosphera2.append(fp)
                                    f_double_primes_refined_nosphera2.append(fdp)
                                    log.write(f"{energy}\tfprefnos:\t{str(fp)}\tfdpref_nos:\t{str(fdp)}\n")
                                else:
                                    fp, fdp = float(l[2].replace('"', '')),float(l[3].replace('">', ''))
                                    f_primes_refined.append(fp)
                                    f_double_primes_refined.append(fdp)
                                    log.write(f"{energy}\tfpref:\t{str(fp)}\tfdpref:\t{str(fdp)}\n")
                

            with open(fr"{base_path}\output_{resolution}.csv", "w+") as outfile:
                outfile.write("Energy,Wavelength,f_prime_sasaki,f_double_prime_sasaki,r1_sasaki,wr2_sasaki,f_prime_henke,f_double_prime_henke,r1_henke,wr2_henke,f_prime_refined,f_double_prime_refined,r1_refined,wr2_refined,f_prime_sasaki_nosphera2,f_double_prime_sasaki_nosphera2,r1_sasaki_nosphera2,wr2_sasaki_nosphera2,f_prime_henke_nosphera2,f_double_prime_henke_nosphera2,r1_henke_nosphera2,wr2_henke_nosphera2,f_prime_refined_nosphera2,f_double_prime_refined_nosphera2,r1_refined_nosphera2,wr2_refined_nosphera2"+"\n")
                for i,_ in enumerate(energies): 
                    print(energies,wls,f_primes,f_double_primes,r1, \
                        wr2,f_primes_henke,f_double_primes_henke,r1_refined, \
                        wr2_henke,f_primes_refined, f_double_primes_refined,r1_refined,wr2_refined, \
                        f_primes_nosphera2,f_double_primes_nosphera2,r1_nosphera2, \
                        wr2_nosphera2,f_primes_henke_nosphera2,f_double_primes_henke_nosphera2,r1_henke_nosphera2, \
                        wr2_henke_nosphera2,f_primes_refined_nosphera2,f_double_primes_refined_nosphera2,r1_refined_nosphera2,wr2_refined_nosphera2)
      
                    outfile.write(
                        str(energies[i]) + "," + str(wls[i]) + "," + str(f_primes[i]) + "," + str(f_double_primes[i]) + "," + str(r1[i]) + "," \
                        + str(wr2[i]) + ","  + str(f_primes_henke[i]) + "," + str(f_double_primes_henke[i]) + "," + str(r1_refined[i]) + "," \
                        + str(wr2_henke[i]) + ","  + str(f_primes_refined[i]) + "," + str(f_double_primes_refined[i]) + "," + str(r1_refined[i]) + "," + str(wr2_refined[i]) \
                        + "," + str(f_primes_nosphera2[i]) + "," + str(f_double_primes_nosphera2[i]) + "," + str(r1_nosphera2[i]) + "," \
                        + str(wr2_nosphera2[i]) + ","  + str(f_primes_henke_nosphera2[i]) + "," + str(f_double_primes_henke_nosphera2[i]) + "," + str(r1_henke_nosphera2[i]) + "," \
                        + str(wr2_henke_nosphera2[i]) + ","  + str(f_primes_refined_nosphera2[i]) + "," + str(f_double_primes_refined_nosphera2[i]) + "," + str(r1_refined_nosphera2[i]) \
                        + "," + str(wr2_refined_nosphera2[i]) + "\n")
            print(f"Wrote results to {base_path}\output_{resolution}.csv")
        
    else:
        energies = []
        wls = []
        f_primes = []
        f_double_primes = []
        r1 = []
        wr2 = []
        f_primes_henke = []
        f_double_primes_henke = []
        r1_henke = []
        wr2_henke = []
        f_primes_refined = []
        f_double_primes_refined = []
        r1_refined = []
        wr2_refined = []
        
        for dir in sub_dirs:
            if name == None:
                name = dir.split("_")[0]
            else:
                if name != dir.split("_")[0] and (name != "tempFAP" or name != "solution"):
                    if name == "tempFAP" or name != "solution":
                        continue
                    raise ValueError("There are folders in your working dir belonging to a different structure!")
            
            energy = int(dir.split("_")[1])
            cor_energy = energy + adjustment
            wl = ret_wl(cor_energy)
            energies.append(cor_energy)
            wls.append(wl)
            
            print(f"Currently working at:\nName:\t{name}\t\t(Corrected) Energy:\t{cor_energy}")
        
            with open(fr"{base_path}\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution_refined.ins", "w+") as output_ref:
                with open(fr"{base_path}\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution_henke.ins", "w+") as output_henke:
                    with open(fr"{base_path}\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution.ins", "w+") as output:
                        with open(fr"{base_path}\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}.ins", "r") as input:
                            curr_in = [line for line in input]
                            state = True
                            for line in curr_in:
                                buffer = []
                                if state:
                                    if line.startswith("CELL"):
                                        buffer = line[:5] + str(round(wl,6))[:7] + line[12:]
                                    else:
                                        buffer = line
                                    if line.startswith("LATT"):
                                        state = False
                                        continue
                                    output.write(buffer)
                                    output_henke.write(buffer)
                                    output_ref.write(buffer)
                                else:
                                    break
                            for line in copy_ins:
                                if line == "REM <olex2.extras>\n":
                                    output_ref.write(line)
                                    output.write(line)
                                    output_henke.write(line)
                                else:
                                    output_ref.write(line)
                                    output.write(line)
                                    output_henke.write(line)
                                    
            nosphera2 = False
            
            ins_path = base_path + fr"\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution.ins"
            ins_path_henke = base_path + fr"\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution_henke.ins"
            ins_disp_ref_path = base_path + fr"\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution_refined.ins"
                
            #Refinement fpr Sasaki Table / Sasaki NoSpherA2
            refine_please(ins_path, nosphera2, resolution, henke = False)
            
            r1.append(float(OV.GetParam("snum.refinement.last_R1")))
            wr2.append(float(OV.GetParam("snum.refinement.last_wR2")))
            
            with open(ins_path, "r") as out_ins:
                for line in out_ins:
                    if line.startswith("DISP") and "Mo" in line:
                        l = line.split()
                        f_primes.append(float(l[2]))
                        f_double_primes.append(float(l[3]))
                    
            refine_please(ins_path_henke, nosphera2, resolution, henke = True)
            
            r1_henke.append(float(OV.GetParam("snum.refinement.last_R1")))
            wr2_henke.append(float(OV.GetParam("snum.refinement.last_wR2")))
            
            with open(ins_path_henke, "r") as out_ins:
                for line in out_ins:
                    if line.startswith("DISP") and "Mo" in line:
                        l = line.split()
                        f_primes_henke.append(float(l[2]))
                        f_double_primes_henke.append(float(l[3]))
        
            refine_please(ins_disp_ref_path, nosphera2, resolution, elements, henke = False)
            
            r1_refined.append(float(OV.GetParam("snum.refinement.last_R1")))
            wr2_refined.append(float(OV.GetParam("snum.refinement.last_wR2"))) 
            
            with open(ins_disp_ref_path, "r") as out_ins:
                for line in out_ins:
                    if line.startswith("REM  <Mo"):
                        l = line.split()
                        fp, fdp = float(l[2].replace('"', '')),float(l[3].replace('">', ''))
                        f_primes_refined.append(fp)
                        f_double_primes_refined.append(fdp)


            
            # with open(ins_disp_ref_path, "r") as out_ins:
            #     for line in out_ins:
            #         if line.startswith("DISP") and "Mo" in line:
            #             l = line.split()
            #             f_primes_refined.append(float(l[2]))
            #             f_double_primes_refined.append(float(l[3]))
            
        with open(fr"{base_path}\output_{resolution}.csv", "w+") as outfile:
            outfile.write("Energy,Wavelength,f_prime_sasaki,f_double_prime_sasaki,r1_sasaki,wr2_sasaki,f_prime_henke,f_double_prime_henke,r1_henke,wr2_henke,f_prime_refined,f_double_prime_refined,r1_refined,wr2_refined\n")
            for i,_ in enumerate(energies): 
                print(energies,wls,f_primes,f_double_primes,r1, \
                    wr2,f_primes_henke,f_double_primes_henke,r1_refined, \
                    wr2_henke,f_primes_refined, f_double_primes_refined,r1_refined,wr2_refined)          
                outfile.write(
                    str(energies[i]) + "," + str(wls[i]) + "," + str(f_primes[i]) + "," + str(f_double_primes[i]) + "," + str(r1[i]) + "," \
                    + str(wr2[i]) + ","  + str(f_primes_henke[i]) + "," + str(f_double_primes_henke[i]) + "," + str(r1_refined[i]) + "," \
                    + str(wr2_henke[i]) + ","  + str(f_primes_refined[i]) + "," + str(f_double_primes_refined[i]) + "," + str(r1_refined[i]) + "," + str(wr2_refined[i]) + "\n")
        print(f"Wrote results to {base_path}\output_keinnosphera2_{resolution}.csv")
            
            
  def extract_disps(self):
      print("Hi")
    # elements=OV.GetParam("fap.element_string")
    # adjustment = OV.GetParam("fap.adjustment_eV")
    # resolution = OV.GetParam("fap.resolution")

    # base_path = FAP_instance.getBasePath()
    # solution_path = FAP_instance.getSolutionPath()

    # if base_path == "" or solution_path == "" :
    #     raise ImportError("Please provide a base and solution path!")

    # #define the string to refine dispersion values for all elements
    # disp_ref_string = f'REM RefineDisp="{elements}"\n'
    
    # #Search for all sub dirs in base dir and storing their names in sub_dirs
    # out = [x for x in os.walk(base_path)]
    # sub_dirs = out[0][1]
    # name = None
    
    # #Define output buffers as empty lists. 
    # energies = []
    # wls = []
    # f_primes = []
    # f_double_primes = []
    # r1 = []
    # wr2 = []
    # f_primes_henke = []
    # f_double_primes_henke = []
    # r1_henke = []
    # wr2_henke = []
    # f_primes_refined = []
    # f_double_primes_refined = []
    # r1_refined = []
    # wr2_refined = []
    # f_primes_nosphera2 = []
    # f_double_primes_nosphera2 = []
    # r1_nosphera2 = []
    # wr2_nosphera2 = []
    # f_primes_henke_nosphera2 = []
    # f_double_primes_henke_nosphera2 = []
    # r1_henke_nosphera2 = []
    # wr2_henke_nosphera2 = []
    # f_primes_refined_nosphera2 = []
    # f_double_primes_refined_nosphera2 = []
    # r1_refined_nosphera2 = []
    # wr2_refined_nosphera2 = []
    
    # for dir in sub_dirs:
    #     energy = int(dir.split("_")[1])
    #     cor_energy = energy + adjustment
    #     wl = ret_wl(cor_energy)
    #     energies.append(cor_energy)
    #     wls.append(wl)
        
    #     global nosphera2_switch
    #     nosphera2_switch = True
    #     for nosphera2 in [False, True]:
          
    #         if nosphera2:
    #             ins_path = base_path + fr"\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution_nosphera2.ins"
    #             ins_path_henke = base_path + fr"\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution_henke_nosphera2.ins"
    #             ins_disp_ref_path = base_path + fr"\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution_refined_nosphera2.ins"

    #         else: 
    #             ins_path = base_path + fr"\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_solution.ins"
    #             ins_path_henke = base_path + fr"\{name}_{energy}\struct\olex2_{name}_{energy}\{name}_{energy}_henke_solution.ins"
    #             ins_disp_ref_path = base_path + fr"\{name}_{energy}\struct\olex2_{name}_{energy}_solution_refined_nosphera2.ins"
                
    #         refine_please2(ins_path, nosphera2, resolution, henke = False)
            
    #         if nosphera2:
    #             r1_nosphera2.append(float(OV.GetParam("snum.refinement.last_R1")))
    #             wr2_nosphera2.append(float(OV.GetParam("snum.refinement.last_wR2")))
    #         else:
    #             r1.append(float(OV.GetParam("snum.refinement.last_R1")))
    #             wr2.append(float(OV.GetParam("snum.refinement.last_wR2")))

    #         with open(ins_path, "r") as out_ins:
    #             for line in out_ins:
    #                 if line.startswith("DISP") and "Mo" in line:
    #                   l = line.split()
    #                   if nosphera2:
    #                     f_primes_nosphera2.append(float(l[2]))
    #                     f_double_primes_nosphera2.append(float(l[3]))
    #                   else:
    #                     f_primes.append(float(l[2]))
    #                     f_double_primes.append(float(l[3]))
                        
    #         refine_please2(ins_path_henke, nosphera2, resolution, henke = True)
            
    #         if nosphera2:
    #             r1_henke_nosphera2.append(float(OV.GetParam("snum.refinement.last_R1")))
    #             wr2_henke_nosphera2.append(float(OV.GetParam("snum.refinement.last_wR2")))
    #         else:
    #             r1_henke.append(float(OV.GetParam("snum.refinement.last_R1")))
    #             wr2_henke.append(float(OV.GetParam("snum.refinement.last_wR2")))

    #         with open(ins_path_henke, "r") as out_ins:
    #             for line in out_ins:
    #                 if line.startswith("DISP") and "Mo" in line:
    #                   l = line.split()
    #                   if nosphera2:
    #                     f_primes_henke_nosphera2.append(float(l[2]))
    #                     f_double_primes_henke_nosphera2.append(float(l[3]))
    #                   else:
    #                     f_primes_henke.append(float(l[2]))
    #                     f_double_primes_henke.append(float(l[3]))
                        
    #         refine_please2(ins_disp_ref_path, nosphera2, resolution,henke = False)
            
    #         with open(ins_disp_ref_path, "r") as out_ins:
    #             for line in out_ins:
    #                 if line.startswith("DISP") and "Mo" in line:
    #                     l = line.split()
    #                     if nosphera2:
    #                         f_primes_refined_nosphera2.append(float(l[2]))
    #                         f_double_primes_refined_nosphera2.append(float(l[3]))
    #                     else:
    #                         f_primes_refined.append(float(l[2]))
    #                         f_double_primes_refined.append(float(l[3]))
                            
    #         with open(fr"{base_path}\output.csv", "w+") as outfile:
    #           outfile.write("Energy,Wavelength,f_prime_sasaki,f_double_prime_sasaki,r1_sasaki,wr2_sasaki,f_prime_henke,f_double_prime_henke,r1_henke,wr2_henke,f_prime_refined,f_double_prime_refined,r1_refined,wr2_refined,f_prime_sasaki_nosphera2,f_double_prime_sasaki_nosphera2,r1_sasaki_nosphera2,wr2_sasaki_nosphera2,f_prime_henke_nosphera2,f_double_prime_henke_nosphera2,r1_henke_nosphera2,wr2_henke_nosphera2,f_prime_refined_nosphera2,f_double_prime_refined_nosphera2,r1_refined_nosphera2,wr2_refined_nosphera2"+"\n")
    #           for i,_ in enumerate(energies): 
    #               print(energies,wls,f_primes,f_double_primes,r1, \
    #                   wr2,f_primes_henke,f_double_primes_henke,r1_refined, \
    #                   wr2_henke,f_primes_refined, f_double_primes_refined,r1_refined,wr2_refined, \
    #                   f_primes_nosphera2,f_double_primes_nosphera2,r1_nosphera2, \
    #                   wr2_nosphera2,f_primes_henke_nosphera2,f_double_primes_henke_nosphera2,r1_henke_nosphera2, \
    #                   wr2_henke_nosphera2,f_primes_refined_nosphera2,f_double_primes_refined_nosphera2,r1_refined_nosphera2,wr2_refined_nosphera2)           
    #               outfile.write(
    #                   str(energies[i]) + "," + str(wls[i]) + "," + str(f_primes[i]) + "," + str(f_double_primes[i]) + "," + str(r1[i]) + "," \
    #                   + str(wr2[i]) + ","  + str(f_primes_henke[i]) + "," + str(f_double_primes_henke[i]) + "," + str(r1_refined[i]) + "," \
    #                   + str(wr2_henke[i]) + ","  + str(f_primes_refined[i]) + "," + str(f_double_primes_refined[i]) + "," + str(r1_refined[i]) + "," + str(wr2_refined[i]) \
    #                   + "," + str(f_primes_nosphera2[i]) + "," + str(f_double_primes_nosphera2[i]) + "," + str(r1_nosphera2[i]) + "," \
    #                   + str(wr2_nosphera2[i]) + ","  + str(f_primes_henke_nosphera2[i]) + "," + str(f_double_primes_henke_nosphera2[i]) + "," + str(r1_henke_nosphera2[i]) + "," \
    #                   + str(wr2_henke_nosphera2[i]) + ","  + str(f_primes_refined_nosphera2[i]) + "," + str(f_double_primes_refined_nosphera2[i]) + "," + str(r1_refined_nosphera2[i]) \
    #                   + "," + str(wr2_refined_nosphera2[i]) + "\n")
                
    # print(f"Wrote results to {solution_path}\output.csv")     
        
def configure_ORCA():
  import time
  olx.xf.EndUpdate()
  if OV.HasGUI():
    olx.Refresh()
  olex.m("spy.NoSpherA2.set_default_cpu_and_mem()")
  OV.SetParam('snum.NoSpherA2.basis_name',"x2c-TZVPP")
  if not nosphera2_switch:
    OV.SetParam('snum.NoSpherA2.file', 'template.tsc')
  OV.SetParam('snum.NoSpherA2.multiplicity', "1")
  OV.SetParam('snum.NoSpherA2.charge', "0")
  OV.SetParam('snum.NoSpherA2.method',"PBE0")
  OV.SetParam('snum.NoSpherA2.h_afix', 'True')
  OV.SetParam('snum.NoSpherA2.h_aniso', 'False')  
  OV.SetParam('snum.NoSpherA2.becke_accuracy',"Normal")
  OV.SetParam('snum.NoSpherA2.ORCA_SCF_Conv',"NormalSCF")
  OV.SetParam('snum.NoSpherA2.ORCA_SCF_Strategy',"NormalConv")
  OV.SetParam('snum.NoSpherA2. Max_HAR_Cycles',"5")
  OV.SetParam('snum.NoSpherA2.ORCA_Solvation',"Vacuum")
  OV.SetParam('snum.NoSpherA2.Relativistic',True)
  OV.SetParam('snum.NoSpherA2.full_HAR',True)
  OV.SetParam('snum.NoSpherA2.use_aspherical',True)
  OV.SetParam('snum.NoSpherA2.Calculate',True)
  OV.SetParam('snum.NoSpherA2.source','ORCA')

def refine_please2(file, nosphera2, resolution, henke = False):     
    olex.m(fr"reap {file}")
    olex.m("refine 20")
    
def refine_please(file, nosphera2, resolution, elements = "", henke = False):     
    olex.m(fr"reap {file}")
    if elements != "":
        free_DISPS(elements)
    olx.AddIns("EXTI")
    olx.AddIns("ACTA")
    olex.m(f"SHEL 999 {resolution}")
    if henke:
      olex.m("gendisp -force -source=Henke")
    else:
      olex.m("gendisp -force -source=Sasaki")
    OV.SetParam('snum.NoSpherA2.use_aspherical',False)
    OV.SetParam('snum.refinement.update_weight',True)
    olex.m("grow")
    olex.m("refine")
    if nosphera2:
      exti = OV.GetExtinction()
      if exti < 0.0001:
        olex.m("delins EXTI")
        olex.m("spy.set_refinement_program(olex2.refine, Gaussian-Newton)")
      else: 
        olex.m("spy.set_refinement_program(olex2.refine, Levenberg-Marquardt)")
      olex.m("refine 20")
      configure_ORCA()
      olex.m("refine")
    else:
      olex.m("refine 30")
      exti = OV.GetExtinction()
      if exti < 0.0001:
          olex.m("delins EXTI")
          olex.m("spy.set_refinement_program(olex2.refine, Gaussian-Newton)")
      else: 
          olex.m("spy.set_refinement_program(olex2.refine, Levenberg-Marquardt)")
      for i in range(9):
          olex.m("refine 99")
      olex.m("refine 20")
      olex.m("refine 20")

def free_DISPS(elements):
    for element in elements:
        olex.m(f"free Disp $Mo")

FAP_instance = FAP()
