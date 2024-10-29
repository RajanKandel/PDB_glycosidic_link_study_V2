import pandas as pd
import numpy as np 

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

import wget

import os 
import re
import glob
import multiprocessing

import subprocess

import json 

################
print ("===============================================================================")
print ("Script: Find glycosidic dihedral, average B-factor, and ring shape of glycans from PDB")
print ("Author:Rajan Kandel <rajan.kandel@uga.edu>")
print ("      *Woods Group, CCRC UGA")
print ("  ")
print ("Currently under development.......")
print ("Last Update: Oct 2024")
print ("")
print('Citations:')
print('1. https://glycam.org/cb/')
print('2. https://glycam.org/portal/gf_home/')
print (" ")
print ("===============================================================================")

################
class parse_gf_result:
    def __init__(self, gf_file, phi_dihedral, psi_dihedral, di_sugar, patterns, chimerax_path, bfmp_path):
        self.gf_file = gf_file

        self.phi_dihedral = phi_dihedral
        self.psi_dihedral = psi_dihedral
        self.di_sugar=di_sugar
        self.patterns1 = patterns[0]
        self.patterns2 = patterns[1]

        self.input_df= pd.DataFrame()
        self.pdbs=[]
        self.pdb_glycan_dict={}
        self.selected_pdbs={}
        self.selected_pdbs2={}

        self.torsions =pd.DataFrame()

        self.chimerax_path= chimerax_path
        self.bfmp_path = bfmp_path
        self.dnl_pdbs_dir= '../../pdbs/'

    def read_gf_file(self):
        input_df = pd.read_csv(self.gf_file, sep=',')

        print('========== Reading_gf_results =====================================')
        pdbs = input_df['pdb'].tolist()
        print('pdbs:', pdbs, '\n # pdbs:',len(pdbs))

        #remove the duplicate pdbID from the list:
        pdbs = list(set(pdbs))
        print('unique pdbs: \n',pdbs, '\n Number of unique pdbs:',len(pdbs))
        self.pdbs=pdbs
        self.input_df= input_df

        print('===================================================================')
    
    def dnl_pdbs(self):
        print('\n============= Downloading pdd file from RCSB pdb ===================')
        download_dir = self.dnl_pdbs_dir

        # Define the base URL for PDB files
        base_url = "https://files.rcsb.org/download/"

        if not os.path.exists(download_dir):
            try:
                os.mkdir(download_dir)

                # Loop through the list of PDB IDs and download the corresponding files
                failed_to_dl_count = 0
                for pdb_id in self.pdbs:
                    pdb_url = base_url + pdb_id + ".pdb"
                    try:
                        wget.download(pdb_url, download_dir)  
                    except :
                        print(pdb_id, "failed to download")
                        failed_to_dl_count += 1
                print("failed to download count", failed_to_dl_count)
            except OSError as error:
                print("Error creating directory:", error)
        else:
            print("Dir: ", download_dir, "- already exists")

        # Loop through the list of PDB IDs and download the not downloaded pdb files
        failed_to_dl_count = 0
        not_dl_pdbs = []
        dl_pdbs = [pdb[-8:-4].upper() for pdb in os.listdir(download_dir)]
        for pdb in self.pdbs:
            if(pdb not in dl_pdbs):
                not_dl_pdbs.append(pdb)

        for pdb_id in not_dl_pdbs:
            pdb_url = base_url + pdb_id + ".pdb"
            try:
                wget.download(pdb_url, download_dir)  
            except :
                print(pdb_id, "failed to download")
                failed_to_dl_count += 1
        print("failed to download count", failed_to_dl_count)

        dnl_pdbs= os.listdir(download_dir)
        pdbs_filtered=[]
        for pdb in self.pdbs:
            pdb_file = (pdb.lower()+'.pdb')
            if pdb_file  in dnl_pdbs:
                pdbs_filtered.append(pdb)

        print(f'pdbs_filtered: {pdbs_filtered}')
        self.pdbs = pdbs_filtered                                                     
        print('===================================================================')

    def create_pdb_glycan_dict(self):
        input_df = self.input_df
        pdb_glycan_dict = {}

        for i in range (input_df.shape[0]):  

            pdb =  input_df['pdb'][i]
            oligo_sequence = input_df['oligo_sequence'][i]
            residue_links = input_df['residue_links'][i]

            if pdb not in pdb_glycan_dict.keys(): 
                glycan_info = {}
                glycan_name = []
                glycan_residueID = []

                glycan_name.append(oligo_sequence)
                glycan_residueID.append(residue_links)

                glycan_info['glycan_name'] = glycan_name
                glycan_info['glycan_residueID'] = glycan_residueID
            
                pdb_glycan_dict[pdb] = glycan_info

            else:
                pdb_glycan_dict[pdb]['glycan_name'].append(oligo_sequence)
                pdb_glycan_dict[pdb]['glycan_residueID'].append(residue_links)

        print(pdb_glycan_dict) 
        print(f'total_pdb_count: {len(pdb_glycan_dict.keys())}')

        self.pdb_glycan_dict = pdb_glycan_dict

    def glycan_count(self):
        #glycan count

        glycosite_count = 0
        for pdb in self.pdb_glycan_dict:
            glycosite_count += len(self.pdb_glycan_dict[pdb]['glycan_name'])

        print(f'glycosite_count: {glycosite_count}')   

    def check_atom_exist(self, pdb, chain, ResName,  ResID, atom):
        #line to search in the pdb HETATM 7294  C2  NAG C   2 
        with open(pdb, "r") as pdb:
            for line in pdb:
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    pattern = rf"{atom}\s*{ResName}\s*{chain}\s*{ResID}\s"
                    matches = re.findall(pattern, line)
                    if len(matches) >1:
                        print(line.rstrip())
                        print(matches)
                    
                    if (len(matches) == 1):
                        # print(line.rstrip())
                        return(1)
        return(0)
        
    def check_pdb_model_count(self, pdb):
        model_count =0
        with open(pdb, "r") as pdb:
            for line in pdb:
                if line.startswith('Model') or line.startswith('MODEL'):
                    # print(line)
                    model_count +=1
        return((model_count>1))
    
    #for given glycan_residueID, identify the glycan_residueName
    def find_glycan_residueName(self, pdb, glycan_residueID):
        for i, resID in enumerate(self.pdb_glycan_dict[pdb]['glycan_residueID']):
            if resID == glycan_residueID:
                return (self.pdb_glycan_dict[pdb]['glycan_name'][i])
        
        return(None)
    

    # Look for specific pattern in the resID and based on that identify the residue that needed to be evaluated further 
    # to handle the multiple pattern found, return a list 
    def required_residues(self, resID, pattern):
        regex_pattern = pattern.replace('(', r'\(').replace(')', r'\)')
        regex_pattern = regex_pattern.replace('[', r'\[').replace(']', r'\]')
        regex_pattern = regex_pattern.replace('*', r'[^-\[\]]+')
        
        # Replace [*] with a pattern that matches anything inside square brackets
        regex_pattern = regex_pattern.replace(r'\[[^-\[\]]+\]', r'\[.*?\]')
        
        # Find all matches of the pattern in the resID, including overlapping ones
        matches = []
        start = 0
        while True:
            match = re.search(regex_pattern, resID[start:])
            if not match:
                break
            matches.append(resID[start:][match.start():match.end()])
            start += match.start() + 1
        
        return matches if matches else None
    
    def check_glycosidic_bond_exist(self, pdb, glycan_residueID, patterns):
        glycan_name_idx=None
        for i, residueID in enumerate(self.pdb_glycan_dict[pdb]['glycan_residueID']):
            if residueID == glycan_residueID:
                glycan_name_idx = i
                break
        if glycan_name_idx !=None:    
            glycan_name= self.pdb_glycan_dict[pdb]['glycan_name'][glycan_name_idx] 
            match=[] 
            for pattern in patterns:
                temp=self.required_residues(glycan_name, pattern)
                if temp!=None:
                    match.extend(temp)
            
            return(len(match))
        else:
            return(0)
        
    #write selected pdbs to a json file
    def write_json_file(self, dict, filename):
        with open(f'{filename}.txt', 'w') as file:
            json.dump(dict, file, indent=4)

    #read_json_file
    def read_json_file(self, json_file):
        with open(json_file, 'r') as file:
            return(json.load(file))

    # from a dict of pdbs with their N-glycosylated ASN, find the pdbs with N-glycans that pass <some filtration criteria>
    def select_pdbs(self):
        selected_pdbs={}

        for pdb in self.pdbs:
            if pdb in self.pdb_glycan_dict:
                temp_glycan_residueID_list = []
                for i, glycan_residueID in enumerate(self.pdb_glycan_dict[pdb]['glycan_residueID']):
                    temp_glycan_residueID_list.append(glycan_residueID)

                if len(temp_glycan_residueID_list) >0:
                    selected_pdbs[pdb] = temp_glycan_residueID_list                              
        
        self.selected_pdbs = selected_pdbs
        self.write_json_file(selected_pdbs, 'selected_pdbs')
        # json_file = 'selected_pdbs.txt'
        # self.selected_pdbs = self.read_json_file(json_file)

    
        
    def generate_chimerax_script_to_compute_dihedral(self):
        selected_pdbs2={}
        temp_oligo=[]
        temp_pattern=[]

        #from a dict of pbs with their N-glycosylated ASN, create chimerax commands for NAG dihedral calc check also if the atom exist 
        warning_pdb =[]
        warning_pdb2=[]
        warning_pdb3=[]
        warning_pdb4=[]
            
        for pdb in self.pdbs:           
            if pdb in self.selected_pdbs:        
                pdbpath = self.dnl_pdbs_dir + pdb.lower() + '.pdb'
                if not (self.check_pdb_model_count(pdbpath)):
                    for resID in self.selected_pdbs[pdb]:
                        results=[]                
                        if(self.check_glycosidic_bond_exist(pdb, resID, self.patterns2)):                    
                            for pattern in self.patterns1:
                                if self.required_residues(resID, pattern) != None:
                                    results.extend(self.required_residues(resID, pattern))  
                        
                        if len(results) >0:
                            is_in_temp_oligo=0
                            for result in results:
                                NAG = [res for res in result.split(sep = '-') ]
                                if (len(NAG) > 1):
                                    i = 0 #position of nag in the list

                                    str2 = re.search(r"\(.*\)", NAG[i]).group(0)[1:-1].split(sep='_')
                                    nag2 = str2[1] + ':' + str2[0]

                                    str1 = re.search(r"\(.*\)", NAG[-1]).group(0)[1:-1].split(sep='_')
                                    nag1 = str1[1] + ':' + str1[0]    

                                    if ((self.check_atom_exist(pdbpath, str2[1], self.di_sugar[0], str2[0], self.phi_dihedral[0])) and
                                        (self.check_atom_exist(pdbpath, str2[1], self.di_sugar[0], str2[0], self.phi_dihedral[1]) ) and
                                        (self.check_atom_exist(pdbpath, str1[1], self.di_sugar[1], str1[0], self.phi_dihedral[2]) ) and
                                        (self.check_atom_exist(pdbpath, str1[1], self.di_sugar[1], str1[0], self.phi_dihedral[3]) ) and
                                        (self.check_atom_exist(pdbpath, str1[1], self.di_sugar[1], str1[0], self.psi_dihedral[3]) )):
     
                                        print ('close #1\n')
                                        print(f"open ../{pdbpath}")

                                        print(f"torsion /{nag2}@{self.phi_dihedral[0]} /{nag2}@{self.phi_dihedral[1]} /{nag1}@{self.phi_dihedral[2]} /{nag1}@{self.phi_dihedral[3]}")
                                        print(f"torsion /{nag2}@{self.psi_dihedral[0]} /{nag1}@{self.psi_dihedral[1]} /{nag1}@{self.psi_dihedral[2]} /{nag1}@{self.psi_dihedral[3]}")

                                        temp_pattern.append(result)

                                        if not is_in_temp_oligo:
                                            temp_oligo.append(resID)
                                            is_in_temp_oligo=1

                                        # just for record
                                        if (len(NAG) > 2):
                                            warning_pdb.append(pdbpath)                         

                                    else:
                                        warning_pdb2.append((pdbpath, results))
                            
                    if len(temp_oligo)>0:
                        selected_pdbs2[pdb.lower()] = temp_oligo
                        temp_oligo=[]
                else:
                    warning_pdb3.append(pdbpath)
            else:
                warning_pdb4.append(pdbpath)           
        print('close #1\n\nexit')
        self.selected_pdbs2 = selected_pdbs2
        self.write_json_file(selected_pdbs2, 'selected_pdbs2')
        ########## write to json file 
        ########## read json file as a dictionary


        # print(warning_pdb2)
        # print(f'pdbs with two or more model: {warning_pdb3}')
        # print(warning_pdb4)

    def custom_split_chimerax_script(self, input_file, output_prefix):
        line_count = 0
        file_count = 0
        current_file = None

        try:
            path = './calc_dihedral_chimerax'
            os.makedirs(path, exist_ok=True)
            print(f'chimerax script is run parallely in {path}')

        except:
            print(f'ERROR: not able to create {path}')

        with open(input_file, 'r') as infile:
            for line in infile:
                if line_count == 0 or (line_count >= 500 and line.startswith('open ')):
                    if current_file:
                        current_file.write('\nexit\n')  # Add 'exit' at the end of the last file
                        current_file.close()
                    file_count += 1
                    current_file = open(f"{path}/{output_prefix}{file_count:02d}.cxc", 'w')
                    current_file.write('close #1\n') 
                    line_count = 0
                
                current_file.write(line)
                line_count += 1

        if current_file:        
            current_file.close()

        print(f"Split into {file_count} files.")

    ##run the chimerax script
    def run_chimerax_script(self, chimerax_exec_file, output):
        command = f"{self.chimerax_path} --nogui {chimerax_exec_file} > {output}"
        subprocess.run(command, shell=True, check=True)


    def run_chimerax_script_parallel(self, script_pattern):
        script_files = sorted(glob.glob(script_pattern))
        output_files = [f.replace('.cxc', '.out') for f in script_files]

        # Create a pool of workers with limited concurrency
        max_concurrent_jobs=6
        with multiprocessing.Pool(processes=max_concurrent_jobs, maxtasksperchild=1) as pool:
            # Run scripts in parallel
            pool.starmap(self.run_chimerax_script, zip(script_files, output_files))
            pool.close()
            pool.join()
        
        print(f"Executed {len(script_files)} ChimeraX scripts in parallel.")
        
        # Stitch outputs sequentially
        combined_output = 'combined_chimerax_output.out'
        with open(combined_output, 'w') as outfile:
            for output_file in output_files:
                with open(output_file, 'r') as infile:
                    outfile.write(infile.read())
        
        print(f"Combined outputs into:  {combined_output}")


    #clean up the output from chimerax
    def clean_up_chimerax_output(self, chimerax_output_path):
        temp_line =[]
        with open(chimerax_output_path, "r") as file:
            for line in file:
                if line.startswith('Executing: open ../') or line.startswith('Torsion angle for atoms'):
                    temp_line.append(line)
                    # print(line)

        with open('temp_file1.txt', "w") as file:
            file.writelines(temp_line) 
        temp_line.clear()
                
        pdb=[]
        pdb2=[]
        phi_site=[]
        psi_site=[]
        torsion1=[]
        torsion2=[]
        with open('temp_file1.txt', 'r') as file:
            prev_pdb=1
            for line in file:
                if line.startswith('Executing: open ../') and prev_pdb==1:
                    pdbID=(line.split(sep='/')[-1])
                    pdb.append(pdbID)
                    prev_pdb=0
                if (f' {self.phi_dihedral[0]} {self.phi_dihedral[1]}') in line:
                    pdb2.append(pdbID[:-1])
                    torsion1.append(float((line.split()[-1]).split(sep='°')[0]))
                    prev_pdb=1
                    phi_site.append(line[25:50])
                if (f' {self.psi_dihedral[3]} is') in line:
                    torsion2.append(float((line.split()[-1]).split(sep='°')[0]))
                    prev_pdb=1
                    psi_site.append(line[25:50])
                if prev_pdb ==0:
                    pdbID=(line.split(sep='/')[-1])
                    pdb[-1]= pdbID
                    continue

        print(f'number of pdb read: {len(pdb2)} \nvalid torsion1 and torsion2: {len(torsion1), len(torsion2)}')

        #convert angle into scale of 0-360
        def convert360(angle):
            conversion=[]
            for measure in angle:

                try:                
                    if(measure>0):
                        conversion.append(measure)
                    elif(measure<0):
                        conversion.append(360+measure)
                except:
                    print(measure)
                
            return(conversion)

        torsions =pd.DataFrame()
        torsions['pdb'] = pdb2
        torsions['phi site'] =phi_site
        torsions['psi site'] =psi_site
        torsions['phi'] = convert360(torsion1)
        torsions['psi'] = convert360(torsion2)

        self.torsions =torsions
        
        return(torsions)

    # find B-factor   
    def find_bfactor_residue(self, pdb, residue):
        with open(f'{self.dnl_pdbs_dir}/{pdb.lower()}.pdb', 'r') as file:
            lines = file.readlines()

            sum_Bfactor = 0
            count_atom = 0

            chainID = residue[-3]
            resID = residue[4:-4]
            resName = residue[0:3]
            # print(chainID, resID, resName)
            pattern1 = rf'^ATOM\s*\d+\s+\S+\s+{resName}\s+{chainID}\s*{resID}'
            pattern2 = rf'^HETATM\s*\d+\s+\S+\s+{resName}\s+{chainID}\s*{resID}'
            for line in lines:
                if re.match(pattern1, line) or re.match(pattern2, line) :
                    sum_Bfactor +=(float(line[60:67]))
                    count_atom += 1
                    
            if count_atom >0:
                average_Bfactor = round(sum_Bfactor / count_atom, 2)
                # print(f'#atom = {count_atom},  sum_Bfactor = {sum_Bfactor}, average Bfactor = {average_Bfactor}')
                return(average_Bfactor)

            else:
                average_Bfactor='nan'
                return(average_Bfactor)
            

    def calc_BFactor(self):

        self.selected_pdbs2 = self.read_json_file('selected_pdbs2.txt')

        temp_bfactor = []
        temp_oligo=[]
        temp_glycam_oligo=[]
        for pdb in self.selected_pdbs2:
            
            for n_oligo in self.selected_pdbs2[pdb]:

                #store the sequences for future references
                glycam_oligo = self.find_glycan_residueName(pdb.upper(), n_oligo)

                results=[]
                if(self.check_glycosidic_bond_exist(pdb.upper(), n_oligo, self.patterns2)):                    
                            for pattern in self.patterns1:
                                if self.required_residues(n_oligo, pattern) != None:
                                    results.extend(self.required_residues(n_oligo, pattern))               

                if len(results) >0:
                    for result in results:
                        NAG = [res for res in result.split(sep = '-') ]
                        if (len(NAG) > 1):
                            i = 0 #position of nag in the list

                            str2 = re.search(r"\(.*\)", NAG[i]).group(0)[1:-1].split(sep='_')
                            nag2 = self.di_sugar[0]+'(' + str2[0] + '_' + str2[1] + '_)'

                            str1 = re.search(r"\(.*\)", NAG[-1]).group(0)[1:-1].split(sep='_')
                            nag1   = self.di_sugar[1]+'(' + str1[0] + '_' + str1[1] + '_)'
                            nag1_1 = self.di_sugar[1]+'(' + str1[0] + '_' + str1[1] + '_)'


                            pdbpath = self.dnl_pdbs_dir + '/' + pdb.lower() + '.pdb'
                            if ((self.check_atom_exist(pdbpath, str2[1], self.di_sugar[0], str2[0], self.phi_dihedral[0])) and
                                (self.check_atom_exist(pdbpath, str2[1], self.di_sugar[0], str2[0], self.phi_dihedral[1])) and
                                (self.check_atom_exist(pdbpath, str1[1], self.di_sugar[1], str1[0], self.phi_dihedral[2])) and
                                (self.check_atom_exist(pdbpath, str1[1], self.di_sugar[1], str1[0], self.phi_dihedral[3])) and
                                (self.check_atom_exist(pdbpath, str1[1], self.di_sugar[1], str1[0], self.psi_dihedral[3]))):

                                residues =[[nag2,],[nag1_1, nag1]]
                                temp1 = []
                                temp_monosacch = []

                                for residue in residues:
                                    Bfactor = self.find_bfactor_residue(pdb,residue[0]) 
                                    if Bfactor == 'nan':
                                        try:
                                            Bfactor = self.find_bfactor_residue(pdb,residue[1])
                                            monosaccharide = residue[1]
                                        except:
                                            monosaccharide = residue[0]
                                    else:
                                        monosaccharide = residue[0]

                                    temp1.append(Bfactor)
                                    temp_monosacch.append(monosaccharide)

                                    temp_oligo.append(n_oligo)
                                    temp_glycam_oligo.append(glycam_oligo)

                                temp_bfactor.append((temp_monosacch, temp1))           

        self.torsions['B factor(Avg)'] = temp_bfactor
        self.torsions['glycan']=temp_oligo[::2]
        self.torsions['glycam name']= temp_glycam_oligo[::2]


            
    # copy the specific ring into the seperate file
    def create_ring_file(self, pdb, chain, ResName,  ResID):
        ring_content = []
        try:
            with open(pdb, "r") as pdb:
                for line in pdb:
                    if line.startswith("HETATM") and (ResID == line[22:26].lstrip()):
                        if ((ResName in line[17:20]) and chain == line[21]):
                            # print(line.rstrip())
                            ring_content.append(line)
                pdb.close()
            
            with open('./ring_file.pdb', 'w') as file:
                for line in ring_content:
                    file.write(line)  
                file.close()

        except FileNotFoundError:
            print(f"Error: PDB file '{pdb}' not found.")

    def clear_file_content(self, file_path):
        with open(file_path, 'w') as file:
            pass  # Do nothing, opening in 'w' mode truncates the file

    #analyze the output file from bfmp run 
    def confirm_ring_shape(self, bfmp_ring_conformation_file, pdbpath):    
        pdb = []
        residue = []
        resName = []
        ring_conformation = []

        with open(bfmp_ring_conformation_file, "r") as file:
            for line in file:           
                if line.startswith(pdbpath):
                    pdb.append(line.split(sep='/')[-1])
                    residue.append(line.split(sep='/')[-1].split(sep=' ')[1][:-1])
                    resName.append('NAG')
                if line.startswith('1'):
                    ring_conformation.append(line.split()[1])
                
        self.clear_file_content(bfmp_ring_conformation_file)
        return((residue,ring_conformation))
    
    #run bfmp script to check for ring puckering
    def check_ring_puckering(self, pdbpath, resIDs):

        bfmp_ring_conformation_file = './monosaccharide_ring_conformations.txt'

        for resID in resIDs:
            if (self.check_atom_exist(pdbpath, resID[-3], resID[0:3], resID[4:-4], 'C1')):
                self.create_ring_file(pdbpath, resID[-3], resID[0:3], resID[4:-4])
                
            content = []
            with open(f'./monosaccharide_bfmp_run.sh', 'r') as file:                        
                for line in file:
                    content.append(line)
                file.close()

            if len(content) <8:
                print('something went wrong here')

            content[2] = 'pdbpath=' + '\'' + 'ring_file.pdb' + '\'' + '\n'
            content[3] = 'original_pdbpath=' + '\'' + pdbpath + '\'' + '\n'
            content[4] = 'sugar_ring_ResID=' + '\''+ resID[4:-4] + '.' + resID[-3] +  '\''+ '\n'
            content[5] = 'bfmp_path=' + '\''+ self.bfmp_path + '\''+ '\n'

            with open('./monosaccharide_bfmp_run.sh', 'w') as file:
                for line in content:
                    file.write(line)
                file.close()

            os.system('./monosaccharide_bfmp_run.sh > bfmp_log.txt 2>&1')
            
        
        #analyze the ring_puckering_output_file  
        ring_shape = self.confirm_ring_shape(bfmp_ring_conformation_file, pdbpath)
        if len(ring_shape[1]) == 2:
            return(ring_shape) 
        else:
            return None   
        
    def find_bfmp_ring_shape(self):
        self.selected_pdbs2 = self.read_json_file('selected_pdbs2.txt')
        #from a dict of pdbs with their N-glycosylated ASN, run BFMP script to find the ring_conformation of GlcNAc
        ring_shape=[]

        with open('monosaccharide_ring_conformations.txt', 'w') as file:
            file.write('')

        for pdb in self.selected_pdbs2:
            
            for n_oligo in self.selected_pdbs2[pdb]:
                results=[]
                if(self.check_glycosidic_bond_exist(pdb.upper(), n_oligo, self.patterns2)):                    
                    for pattern in self.patterns1:
                        if self.required_residues(n_oligo, pattern) != None:
                            results.extend(self.required_residues(n_oligo, pattern))
            
                if len(results) >0:
                    for result in results:
                        NAG = [res for res in result.split(sep = '-') ]

                        if (len(NAG) > 1 ): 
                            i = 0 #position of glc in the list

                            str2 = re.search(r"\(.*\)", NAG[i]).group(0)[1:-1].split(sep='_')
                            nag2 = self.di_sugar[0]+'(' + str2[0] + '_' + str2[1] + '_)'

                            str1 = re.search(r"\(.*\)", NAG[-1]).group(0)[1:-1].split(sep='_')
                            nag1   = self.di_sugar[1]+'(' + str1[0] + '_' + str1[1] + '_)'
                            nag1_1 = self.di_sugar[1]+'(' + str1[0] + '_' + str1[1] + '_)'
                        
                            pdbpath = self.dnl_pdbs_dir+'/' + pdb.lower() + '.pdb'
                            if ((self.check_atom_exist(pdbpath, str2[1], self.di_sugar[0], str2[0], self.phi_dihedral[0])) and
                                (self.check_atom_exist(pdbpath, str2[1], self.di_sugar[0], str2[0], self.phi_dihedral[1])) and
                                (self.check_atom_exist(pdbpath, str1[1], self.di_sugar[1], str1[0], self.phi_dihedral[2])) and
                                (self.check_atom_exist(pdbpath, str1[1], self.di_sugar[1], str1[0], self.phi_dihedral[3])) and
                                (self.check_atom_exist(pdbpath, str1[1], self.di_sugar[1], str1[0], self.psi_dihedral[3]))):

                                residues =[[nag2,],[nag1, nag1_1 ]]
                                
                                ring_conformation_check = self.check_ring_puckering(pdbpath,[residues[0][0],residues[1][0]])
                                if ring_conformation_check == None:
                                    ring_conformation_check = self.check_ring_puckering(pdbpath,[residues[0][0],residues[1][1]])
                                    print(ring_conformation_check)

                                print(ring_conformation_check)  
                                ring_shape.append(ring_conformation_check)    

        self.torsions['BFMP ring shape'] = ring_shape

    # find the atom number of atom
    def find_atom_number(self, pdb, chain, ResName,  ResID, atom):
        #line to search in the pdb HETATM 7294  C2  NAG C   2 
        with open(pdb, "r") as pdb:
            for line in pdb:
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    pattern = rf"{atom}\s*{ResName}\s*{chain}\s*{ResID}\s"
                    matches = re.findall(pattern, line)
                    if len(matches) >1:
                        print(line.rstrip())
                        print(matches)
                    
                    if (len(matches) == 1):
                        # print(line.rstrip())
                        atom_number = int(line[6:11])
                        return(atom_number)
        return(False)
    
    
    def check_connect_record(self, pdb_file, atom1, atom2):
        try:
            with open(pdb_file, 'r') as pdb:
                conect_lines_found = False
                for line in pdb:
                    if line.startswith('CONECT'):
                        conect_lines_found = True
                        connect_record = line.split()[1:]
                        if atom1 in connect_record and atom2 in connect_record:
                            print(f"Connection found between {atom1} and {atom2}")
                            return 1
                        if atom1 in line and atom2 in line:
                            print(f"Connection found between(2) {atom1} and {atom2}")
                            return 1
                
                if not conect_lines_found:
                    print("No CONECT lines found in the PDB file.")
                else:
                    print(f"No connection found between {atom1} and {atom2}")            
                return 0
        except FileNotFoundError:
            print(f"Error: File {pdb_file} not found.")
            return 0
        except Exception as e:
            print(f"An error occurred: {e}")
            return 0
        

    def match_glycan_tree(self,):

        rows_to_drop = []

        for index, row in self.torsions.iterrows():
            # Access row data
            # Example: print(row['column_name'])
            pdb  = self.dnl_pdbs_dir + '/' +row['pdb']            
            sugar1 = row['B factor(Avg)'][0][0]
            sugar2 = row['B factor(Avg)'][0][1]

            
            chain1 = sugar1[-3]
            ResName1 = sugar1[0:3]
            ResID1 = sugar1[4:-4]
            atom1= self.phi_dihedral[1]
            
            chain2 = sugar2[-3]
            ResName2 = sugar2[0:3]
            ResID2 = sugar2[4:-4]
            atom2=self.phi_dihedral[2]

            # print(pdb, sugar1, sugar2)
            # print(pdb, chain1, ResName1,  ResID1, atom1 )

            if self.check_atom_exist(pdb, chain1, ResName1,  ResID1, atom1 ) and self.check_atom_exist(pdb, chain2, ResName2,  ResID2, atom2 ):
                atomID1 = str(self.find_atom_number(pdb, chain1, ResName1,  ResID1, atom1 ))
                atomID2 = str(self.find_atom_number(pdb, chain2, ResName2,  ResID2, atom2 ))
                # print(atom1, atom2)
                print(pdb, chain1, ResName1,  ResID1, atom1)
                print(pdb, chain2, ResName2,  ResID2, atom2)
                print(atomID1, atomID2)

                if not self.check_connect_record(pdb, atomID1, atomID2):
                    print(index)
                    rows_to_drop.append(index)
                # else:
                    # print(index)

        # Drop the rows outside the loop
        df = self.torsions.drop(rows_to_drop)

        # Reset the index if needed
        df = df.reset_index(drop=True)

        print(f'droped rows: {rows_to_drop}')
        return (df)
    


    # on the scale of 0-360
    def plot(self, torsions, title, max_histogram_scale, step_histogram_tick):

        fig, ((ax3, ax4), (ax1, ax2)) = plt.subplots(2, 2, figsize=(8,8))

        x =torsions['psi']
        y =torsions['phi']

        ax1.scatter(x, y, linewidth=1, linestyle='-', color ='black',marker='.', s=15)
        ax1.set_xlim(xmin=(0), xmax=360)
        ax1.set_ylim(ymin=(0), ymax=360)

        # set the x and y axis labels
        ax1.set_xlabel('ψ', fontsize=12)
        ax1.set_ylabel('Φ', fontsize=12)

        ax1.set_xticks(range(0,361,60))
        ax1.set_yticks(range(0,361,60))
        ax1.tick_params(axis='both', which='major', labelsize=10)


        #####################
        # Plot the histogram of the frequency on the right subplot (ax2)
        sns.histplot(y=y, ax=ax2, bins=20, color ='black', kde=True)
        ax2.set_xlabel('Frequency', fontsize =10)
        ax2.set_ylabel('') 
        ax2.set_xlim(xmin = 0, xmax = max_histogram_scale)
        ax2.set_ylim(ymin = 0, ymax =360)
        ax2.set_yticks(range(0,361,60))
        ax2.tick_params(axis='both', which='both', labelleft =False)
        ax2.xaxis.set_major_locator(ticker.MultipleLocator(base=step_histogram_tick))  # Set x-ticks at a distance of 50,000
            

        #####################
        # Plot the histogram of the frequency on the top subplot (ax3)
        sns.histplot(x=x, ax=ax3, bins=20, color ='black', kde=True)
        ax3.set_ylabel('Frequency', fontsize =10)
        ax3.set_xlabel('') 
        ax3.set_ylim(ymin = 0, ymax = max_histogram_scale)
        ax3.set_xlim(xmin = 0, xmax =360)
        ax3.set_xticks(range(0,361,60))
        ax3.tick_params(axis='both', which='both', labelbottom = False)  # Set tick font size
        ax3.yaxis.set_major_locator(ticker.MultipleLocator(base=step_histogram_tick))  # Set x-ticks at a distance of 50,000

        # Hide ax4 
        ax4.axis('off')

        plt.subplots_adjust(wspace=0.06, hspace=0.06)
        fig.suptitle(title, fontsize=20)

        # Show grid lines on the plots
        ax1.grid(True)
        ax2.grid(True)
        ax3.grid(True)
        ax4.grid(True)

        # show the plot
        plt.show()