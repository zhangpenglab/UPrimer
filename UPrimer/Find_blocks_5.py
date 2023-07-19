#!/usr/bin/python
# -*- coding: UTF-8 -*-

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
import glob
import os, shutil
import argparse
linesep = os.linesep

def Get_species_list(species_list):
        reference_list = []
        ingroup_list = []
        outgroup_list = []

        input_species_file = os.path.join(os.getcwd(),species_list)
        lines = open(input_species_file,'r').readlines()

        in_number = 0
        out_number = 0
        for i in range(0,len(lines)):
                line_info = lines[i].strip()
                if line_info == '##Ingroup species##':
                        in_number = i
                if line_info == '##Outgroup species##':
                        out_number  = i

        reference_list.append(lines[1].strip())

        for line in lines[in_number+1:out_number]:
                ingroup_list.append(line.strip() + '_genome')

        for line in lines[out_number + 1:]:
                outgroup_list.append(line.strip() + '_cds')

        taxa_list = ingroup_list + outgroup_list

        return reference_list, ingroup_list, outgroup_list, taxa_list
    
def Get_nucleotide_MSA_file(file_pep):
    file_list = os.path.split(file_pep)
    file_nucl = os.path.join(file_list[0].replace('3.5.Refined_peptide_aligned_MSAs', '3.5.Refined_nucleotide_aligned_MSAs'), file_list[1].replace('_pep', '_nucl'))
    return file_nucl

def Calculate_block_identity(s):
    s = s.upper()
    length = len(s)
    majority = 0
    gap_num = s.count('-')
    occ_num = s.count('X')

    
    all_per = float(gap_num + occ_num) / len(s) * 100
    if all_per > 50:
            return 0
    else:
            if float(occ_num/length) < 0.2 and float(gap_num/length) < 0.2:
                    
                    s = s.replace('-','').replace('X','')
                    temp_s = s
                    m = len(set(s))
                    
                    for i in range(0, len(s)):
                            n = temp_s.count(s[i])			
                            temp_s = temp_s.replace(s[i], '') 
                            
                            if majority < n: 
                                    AA = s[i]
                                    majority = n
                                    
                            if temp_s == '':
                                    return float(majority) * 1.0 / length
            else:
                   return 0 

def Calculate_taxa_number(fh):
    ingroup = 0
    outgroup = 0
    for records in SeqIO.parse(fh, 'fasta'):
        if 'Ingroup_' in str(records.id):
            ingroup += 1
        elif 'Outgroup_' in str(records.id):
            outgroup += 1
        else:
            pass
    return ingroup, outgroup
    
def Select_primer_design_region(filename, out_pep, out_nucl, Expect_identity, sp_list):
    file_nucl = Get_nucleotide_MSA_file(filename)
    records_pep = AlignIO.read(filename, 'fasta')
    records_nucl = AlignIO.read(file_nucl, 'fasta')
    
    cut_start = 0
    cut_end = 1
    MinLen = 100
    WindowWidth = 7
    for i in range(0, len(records_pep[0].seq) - WindowWidth):
        total_identity = 0
        for j in range(0, WindowWidth):
            total_identity = total_identity + Calculate_block_identity(records_pep[:, i + j])
        if total_identity/WindowWidth >= float(Expect_identity):
            cut_start = i
            break
        
    for i in range(len(records_pep[0].seq) - 1, WindowWidth, -1):
        total_identity = 0
        for j in range(0, WindowWidth):
            total_identity = total_identity + Calculate_block_identity(records_pep[:, i - j])
        if total_identity/WindowWidth >= float(Expect_identity):
            cut_end = i + 1
            break

    length = cut_end - cut_start + 1
    in_out = Calculate_taxa_number(filename)
    taxa_num_in_msa = len(records_pep)
    
    if length >= MinLen:
        pep_out = os.path.join(out_pep, os.path.split(filename)[1])
        nucl_out = os.path.join(out_nucl, os.path.split(filename)[1].replace('_pep', '_nucl'))
        
        my_alignments_pep = records_pep[:, cut_start:cut_end]
        my_alignments_nucl = records_nucl[:, cut_start*3:cut_end*3]
        
        AlignIO.write(my_alignments_pep, pep_out, "fasta")
        AlignIO.write(my_alignments_nucl, nucl_out, "fasta")
        
        line = os.path.split(filename)[1].replace('Remove_rogue_taxa_', '').replace('.fasta', '') + '\t' + str(len(sp_list) - 1)\
               + '\t' + str(taxa_num_in_msa) + '\t' + str(in_out[0]) + '\t' + str(in_out[1]) + '\tYes\t' + str(cut_start)+'\t'\
               + str(cut_end) + '\t'+str(length) + linesep
        return line
    else:
        line = os.path.split(filename)[1].replace('Remove_rogue_taxa_', '').replace('.fasta', '') + '\t' + str(len(sp_list) -1 )\
               + '\t' + str(taxa_num_in_msa) + '\t' + str(in_out[0]) + '\t' + str(in_out[1]) + '\tNo\t' + str(cut_start)+'\t'\
               + str(cut_end) + '\t'+str(length) + linesep
        return line
    
def Pipeline(species_list, Expect_identity):
    candidate_xls = '4.Information_for_candidate_iden' + str(Expect_identity) + 'MSAs_for_primer_design.xls'
    filename_lists = glob.glob(os.path.join(os.getcwd(), '3.5.Refined_peptide_aligned_MSAs', '*.fasta'))
    
    candidate_pep = '4.Candidate_peptide_iden' + str(Expect_identity) + '_MSAs_for_primer_design'
    candidate_nucl = '4.Candidate_nucleotide_iden' + str(Expect_identity) + '_MSAs_for_primer_design'
    
    path_pep = os.path.join(os.getcwd(), candidate_pep)
    path_nucl = os.path.join(os.getcwd(), candidate_nucl)
    if not os.path.exists(path_pep):
        os.mkdir(path_pep)
    if not os.path.exists(path_nucl):
        os.mkdir(path_nucl)

    sp_list = os.path.join(os.getcwd(),species_list)
    reference_list = Get_species_list(sp_list)[0]
    ingroup_list =  Get_species_list(sp_list)[1]
    outgroup_list =  Get_species_list(sp_list)[2]
    sp_list = Get_species_list(sp_list)[3]
    

##    with open(species_list, 'r') as fh_sp:
##        for line in fh_sp.readlines()[1:]:
##            sp = line.strip().replace('.fasta', '')
##            sp_list.append(sp)
            
    with open(candidate_xls, 'a') as f1:
        li = ['MSA name',
              'No. of species used for developing primer sets',
              'No. of total species in MSA',
              'No. of ingroups in MSA',
              'No. of outgroups in MSA',
              'Whether finding forward and reverse conserved 7AA blocks in MSA (Yes or No)',
              'Primer start position in trimmed MSA',
              'Primer end position in trimmed MSA',
              'Final MSA length (from start to end)']
        line_1 = '\t'.join(li) + linesep
        f1.write(line_1)
        for file in filename_lists:
            line_2 = str(Select_primer_design_region(file, path_pep, path_nucl, Expect_identity, sp_list))
            f1.write(line_2)
    return candidate_pep
            
if __name__ == '__main__':
    dirpath_o = os.path.abspath(os.path.join(os.path.dirname('__file__'), os.path.pardir))
    os.chdir(dirpath_o)
    parser = argparse.ArgumentParser(description='''python Find_blocks_5.py -F species_list.txt''')

    parser.add_argument('--file', '-F', help='A text table containing information on the reference species, ingroup species, and outgroup species (filename: species_list.txt)', required=True)
    parser.add_argument('--exp', '-E', help='Value of the expect identity of all species in 7AA/8AA blocks', default=float(0.5))
    args = parser.parse_args()
    print args
    try:
        print Pipeline(args.file, args.exp)
    except Exception as e:
        print(e)
    print 'The work of generating candidate MSAs for primer designing is complete!'
