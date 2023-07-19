#!/usr/bin/python
# -*- coding: UTF-8 -*-

import glob, os
import shutil
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import re
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

def Creat_ingroup_weight_ratio_dic(species_list):
    
	
	sp_list = os.path.join(os.getcwd(),species_list)
	ingroup_list =  Get_species_list(sp_list)[1]
	outgroup_list =  Get_species_list(sp_list)[2]
	
	weight_ratio_min = 0.5
	weight_ratio_max = 1         
	ingroup_number = len(ingroup_list)
	       
	ingroup_max = ingroup_number + 1 #inculuding reference
	ingroup_min = 2
	Tolerances = (weight_ratio_max - weight_ratio_min)/(ingroup_max - 1)

	count = 0
	ingroup_num = 2
	ingroup_weight_ratio_dic = {}
	while ingroup_num < ingroup_max:

		ingroup_weight_ratio_dic[ingroup_num] = 0.5 + count * Tolerances
		ingroup_num += 1
		count += 1
	
	ingroup_weight_ratio_dic[ingroup_max] = 1
	return ingroup_weight_ratio_dic


def Amplicon_score_weight_dic(max_score, min_score):

        max_weight = 1.0
        min_weight = 0.5

        number = int(max_score * 100 - min_score * 100)  
        Tolerances = float(max_weight - min_weight)/float(number)

        count = 0
        amplicon_score_weight_dic = {}
        score = round(float(min_score),2)

        for i in range(0,int(number) + 2):
                amplicon_score_weight_dic[round(score,2)] = round(float(0.5 + count * Tolerances) * 100,2)
                score += 0.01
                count += 1

        return amplicon_score_weight_dic

def Information_score_weight_dic(max_score, min_score):
        max_weight = 1.0
        min_weight = 0.5

        number = int(max_score * 100 - min_score * 100) 
        Tolerances = float(max_weight - min_weight)/float(number)
        
        count = 0
        information_score_weight_dic = {}
        score = round(float(min_score),8)

        for i in range(0,int(number) + 2):
                information_score_weight_dic[round(score,8)] = round(float((0.5 + count * Tolerances) * 100), 8)

                score += 0.01
                count += 1
        return information_score_weight_dic
    	    

def Simplify_txt_for_nested_primers(txt_file, Max_length, Min_length, path):

    print txt_file, 'start to simplify!'

   
    suffix = '_simp.txt'
    out_file = os.path.split(txt_file)[1].replace('.txt', suffix)
    path_out = os.path.join(path, out_file)
      
    if not os.path.exists(path_out):

	with open(txt_file) as txt_in:
                
	    line1 = txt_in.readline().strip()#title
	    count = 0

	    fh_out = open(path_out,'a')
	    fh_out.write(line1 +  linesep)

	    for lines in txt_in.readlines():
		line = lines.split('\t')
		Var_length = float(line[11].strip()) / 3
                if count <= 10:
                    if Min_length <= Var_length and Var_length <= Max_length:
                        fh_out.write(lines)
                        count += 1
                    else:
                        pass
                else:
                    break
            fh_out.close()
                    
    else:
        print out_file, ' have exists!!!'                              

def Stat_primer_set_feature(simplify_file_dir):
    target_length_list = []
    files = glob.glob(os.path.join(os.getcwd(), simplify_file_dir, '*.txt'))#Folder: 5.5.Simplifed_information_for_nested-PCR_primers
    for File in files:
        line = open(File,'r').readlines()[1]
        length = line.strip().split('\t')[11]
        target_length_list.append(length)

    total = 0
    for length in target_length_list:
        total += int(length)

    mean_length = int(total)/int(len(target_length_list))

    return str(len(target_length_list)), str(mean_length)

	    
def Sorted_Nested_primers(out_file_1, out_file_2, dir_file, PId, File, taxa_group):# summrazy primer information

    dic = {}
    Score_sorted = {}
    with open(out_file_1, 'w') as fh_out_1:
	line1 = '\t'.join(['NPCL ID',
                           'MSA ID',
			   'Predicted Amplicon Length (bp)',
			   'No. of Ingroups in MSA',		
			   'F2 Primer Start Position in the Peptide MSA',
			   'R2 Primer Start Position in the Peptide MSA',
			   'Peptide Sequence of F2 Primer (5\' to 3\')',	
			   'Peptide Sequence of R2 Primer (5\' to 3\')',						   
			   'Nucleotide Sequence of F2 Primer (5\' to 3\')',
			   'Nucleotide Sequence of R2 primer (5\' to 3\',after reverse completement)',
			   'Degeneracy of F2 Primer',
			   'Degeneracy of R2 Primer',
			   'Ratio of Variable Colunms for F2R2',
			   'Information Score for F2R2',
			   'PCR Score for F2R2',
                           'Identify Score for F2R2',
			   'F1 Primer Start Position in the Peptide MSA',
			   'R1 Primer Start Position in the Peptide MSA',
			   'Peptide Sequence of F1 Primer (5\' to 3\')',
			   'Peptide Sequence of R1 Primer (5\' to 3\')',
			   'Nucleotide Sequence of F1 Primer (5\' to 3\')',
			   'Nucleotide Sequence of R1 Primer (5\' to 3\',after reverse completement)',
			   'Degeneracy of F1 Primer',
			   'Degeneracy of R1 Primer',
			   'PCR Score for F1R1',
                           'Identify Score for F1R1',
			   'Final Score for nested primer pairs among different MSAs'])
	
	fh_out_1.write(line1 + '\n')
	files = glob.glob(os.path.join(os.getcwd(), dir_file, '*.txt'))
	Amplicon_score_list = []
	Information_score_list = []
	sorted_primer_list = []
	
	for file in files:
            with open(file, 'r') as fh_in:
		all_line = fh_in.readlines()
		lines = all_line[1].strip()#for the highest-scoring primer pair
                line = str(lines).split('\t')
                Amplicon_score_for_F2R2 = float(line[17])
                Amplicon_score_for_F1R1 = float(line[30])
                Amplicon_score = Amplicon_score_for_F2R2 + Amplicon_score_for_F1R1                
                Information_score = float(line[14])

                Amplicon_score_list.append(round(Amplicon_score,2))
                Information_score_list.append(round(Information_score,2))

        
	for file in files:
	    file_name = str(os.path.split(file)[1])
	    
	    AA_num = file_name.split('_')[-5].replace('AA', '')
	    name = file_name
	    align_name = file_name.split('_8AA_primer_F1R1_simp.txt')[0]

	    with open(file, 'r') as fh_in:
		all_line = fh_in.readlines()

		lines = all_line[1].strip() + '\t' + align_name
		line = str(lines).split('\t')

                Amplicon_score_for_F2R2 = float(line[17])
                Amplicon_score_for_F1R1 = float(line[30])

                Amplicon_score = round(Amplicon_score_for_F2R2 + Amplicon_score_for_F1R1, 2)              
                Information_score = round(float(line[14]),2)

		Weighted_Amplicon_score = Amplicon_score_weight_dic(max(Amplicon_score_list), min(Amplicon_score_list))[Amplicon_score]
		Weighted_Information_score = Information_score_weight_dic(max(Information_score_list),min(Information_score_list))[Information_score]         
		
		Creat_ingroup_weight_ratio_dic(File)
		ingoup_num = float(line[12])
		ingroup_weight_ratio = Creat_ingroup_weight_ratio_dic(File)[ingoup_num]
		Final_Score = (float(PId) * Weighted_Amplicon_score + Weighted_Information_score)/(float(PId) + 1)  * ingroup_weight_ratio
		lines = all_line[1].strip() + '\t' + align_name + '\t' + str(Final_Score)
		dic[lines] = Final_Score

        Score_sorted = sorted(dic.items(), key=lambda x: x[1], reverse=True)
        count = 0
        for i in range(0, len(Score_sorted)):
            records = Score_sorted[i][0].split('\t')
            count += 1            
            NPCL_name = str(taxa_group)[0:3].capitalize() + '_' + str(count)
            re_record = '\t'.join([str(NPCL_name),str(records[-2]),records[11],records[12], records[0], records[1], records[3],
                                   records[7],records[2], records[6],records[5], records[9],
                                   records[13],records[14],records[15],records[16],
                                   records[18], records[23],
                                   records[20], records[25], records[19], records[24], records[21], records[26], records[28], records[29], records[-1]])
            fh_out_1.write(re_record + linesep)
            sorted_primer_list.append(str(records[-2]))



    with open(out_file_2, 'w') as fh_out_2:
	line2 = '\t'.join(['MSA ID',
			   'Predicted Amplicon Length (bp)',
			   'No. of Ingroups in MSA',		
			   'F2 Primer Start Position in the Peptide MSA',
			   'R2 Primer Start Position in the Peptide MSA',
			   'Peptide Sequence of F2 Primer (5\' to 3\')',	
			   'Peptide Sequence of R2 Primer (5\' to 3\')',						   
			   'Nucleotide Sequence of F2 Primer (5\' to 3\')',
			   'Nucleotide Sequence of R2 primer (5\' to 3\',after reverse completement)',
			   'Degeneracy of F2 Primer',
			   'Degeneracy of R2 Primer',
			   'Ratio of Variable Colunms for F2R2',
			   'Information Score for F2R2',
			   'PCR Score for F2R2',
                           'Identify Score for F2R2',
			   'F1 Primer Start Position in the Peptide MSA',
			   'R1 Primer Start Position in the Peptide MSA',
			   'Peptide Sequence of F1 Primer (5\' to 3\')',
			   'Peptide Sequence of R1 Primer (5\' to 3\')',
			   'Nucleotide Sequence of F1 Primer (5\' to 3\')',
			   'Nucleotide Sequence of R1 Primer (5\' to 3\',after reverse completement)',
			   'Degeneracy of F1 Primer',
			   'Degeneracy of R1 Primer',
			   'PCR Score for F1R1',
                           'Identify Score for F1R1'])
	fh_out_2.write(line2 + '\n')

	for name in sorted_primer_list:
            file_basename = name + '_8AA_primer_F1R1_simp.txt'#
            file_in = os.path.join(os.getcwd(), dir_file, file_basename)
            lines = open(file_in,'r').readlines()
            count = 0
            for line in lines[1:]:
                records = line.strip().split('\t')

                if count == 0:                         
                    re_record = '\t'.join([str(name),records[11],records[12], records[0], records[1], records[3],
                                           records[7],records[2], records[6],records[5], records[9],
                                           records[13],records[14],records[15],records[16],
                                           records[18], records[23],
                                           records[20], records[25], records[19], records[24], records[21], records[26], records[28], records[29]])                
                
                    fh_out_2.write(re_record + linesep)
                    count += 1
                else:
                    new_name = name + '_' + str(count)
                    re_record = '\t'.join([str(new_name),records[11],records[12], records[0], records[1], records[3],
                                           records[7],records[2], records[6],records[5], records[9],
                                           records[13],records[14],records[15],records[16],
                                           records[18], records[23],
                                           records[20], records[25], records[19], records[24], records[21], records[26], records[28], records[29]])                
                
                    fh_out_2.write(re_record + linesep)
                    count += 1



def Synthezise_nested_primer_result(file_in,file_out,taxa_group):
    with open(file_out, 'w') as fh_out:
	with open(file_in, 'r') as fh_in:
	    line2 = '\t'.join(['Primer name',
			       'Nucleotide Sequence of Primer (5\' to 3\')',
                               'Number of base pair'])
	    
	    fh_out.write(line2 + '\n')
	    record_num = 0
	    for lines in fh_in.readlines()[1:]:
                record_num += 1
		line = lines.strip().split('\t')

                original_primer_name_F1 = str(line[0]) + '_F1'
		synthezise_primer_name_F1 = '_'.join([str(taxa_group)[0:3].capitalize(),str(record_num),'F1'])#1_Het_1_XP12312.1_F1
		seq_F1 = str(line[19])
		numBp_F1 = str(len(line[19]))
		write_record_F1 = '\t'.join([synthezise_primer_name_F1,seq_F1,numBp_F1])

                original_primer_name_R1 = str(line[0]) + '_R1'
		synthezise_primer_name_R1 = '_'.join([str(taxa_group)[0:3].capitalize(),str(record_num),'R1'])#1_Het_1_XP12312.1_R1
		seq_R1 = str(line[20])
		numBp_R1 = str(len(line[20]))
		write_record_R1 = '\t'.join([synthezise_primer_name_R1,seq_R1,numBp_R1])

                original_primer_name_F2 = str(line[0]) + '_F2'		
		synthezise_primer_name_F2 = '_'.join([str(taxa_group)[0:3].capitalize(),str(record_num),'F2'])#1_Het_1_XP12312.1_F2
		seq_F2 = str(line[7])
		numBp_F2 = str(len(line[7]))
		write_record_F2 = '\t'.join([synthezise_primer_name_F2,seq_F2,numBp_F2])

                original_primer_name_R2 = str(line[0]) + '_R2'		
		synthezise_primer_name_R2 = '_'.join([str(taxa_group)[0:3].capitalize(),str(record_num),'R2'])#1_Het_1_XP12312.1_F2
		seq_R2 = str(line[8])
		numBp_R2 = str(len(line[8]))
		write_record_R2 = '\t'.join([synthezise_primer_name_R2,seq_R2,numBp_R2])		

		fh_out.write(write_record_F1 + '\n' + write_record_R1 + '\n'+ write_record_F2 + '\n' + write_record_R2 + '\n')


		
def Extract_reference_sequences(file_in, fasta_out, dir, num, taxa_group):
    pro_files = glob.glob(os.path.join(os.getcwd(), dir, '*.fasta'))
    nucl_files = glob.glob(os.path.join(os.getcwd(), dir.replace('_peptide', '_nucleotide'), '*.fasta'))

    ref_cds_file = glob.glob(os.path.join(os.getcwd(), 'Reference','single_copy_*_cds.fasta'))
    ref_pep_file = glob.glob(os.path.join(os.getcwd(), 'Reference','single_copy_*_pep.fasta'))


    nucl_fasta_out_full = fasta_out + ' nucleotide sequences of NPCLs for ' + str(taxa_group) +'.fasta'
    pro_fasta_out_full = fasta_out + ' peptide sequences of NPCLs for ' + str(taxa_group) +'.fasta'

    
    ref_cds_dic = {}
    ref_pep_dic = {}
    #creat ref CDS dic
    for record in SeqIO.parse(ref_cds_file[0],'fasta'):
        if'_unknown' not in record.description:
            pre_msa_name = str(record.description)
            mas_name = pre_msa_name.split('_')[0] + '_' + pre_msa_name.split('_')[1] + '_' + \
                       pre_msa_name.split('_')[6].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','') + '_' +\
                       pre_msa_name.split('_')[7].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','')
            ref_cds_dic[mas_name] = str(record.seq)

            
        elif'_unknown' in record.description:
            pre_msa_name = str(record.description).replace('_unknown','')
            mas_name = pre_msa_name.split('_')[0] + '_' + \
                       pre_msa_name.split('_')[2].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','') +\
                       '_' + pre_msa_name.split('_')[3].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','')
            ref_cds_dic[mas_name] = str(record.seq)

            
        else:
            seq_name_info = str(record.description).split('|')
            if len(seq_name_info) == 3:                
                for i in seq_name_info:
                    if '-R' and '-E' in i:#AGLA000415|AGLA000415-RA|AGLA000415-RA-E1, i = AGLA000415-RA-E1
                        ref_cds_dic[i] = str(record.seq)
                    elif '-E' in i and '-R' not in i:#FBgn0036133|FBpp0075940|FBtr0076210-E2|Q9VTF5, i = FBtr0076210-E2
                        ref_cds_dic[i] = str(record.seq)
            else:
                ref_cds_dic[seq_name_info[-1]] = str(record.seq)

                    

    #creat ref pep dic
    for record in SeqIO.parse(ref_pep_file[0],'fasta'):
        if '_unknown' not in record.description:
            pre_msa_name = str(record.description)
            mas_name = pre_msa_name.split('_')[0] + '_' + pre_msa_name.split('_')[1] + '_' + \
                       pre_msa_name.split('_')[6].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','') + '_' +\
                       pre_msa_name.split('_')[7].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','')
            ref_pep_dic[mas_name] = str(record.seq)
        elif '_unknown' in record.description:
            pre_msa_name = str(record.description).replace('_unknown','')
            mas_name = pre_msa_name.split('_')[0] + '_' + \
                       pre_msa_name.split('_')[2].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','') +\
                       '_' + pre_msa_name.split('_')[3].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','')
            ref_pep_dic[mas_name] = str(record.seq)
        else:
            seq_name_info = str(record.description).split('|')
            if len(seq_name_info) == 3:
                for i in seq_name_info:
                    if '-R' and '-E' in i:#AGLA000415|AGLA000415-RA|AGLA000415-RA-E1, i = AGLA000415-RA-E1
                        ref_pep_dic[i] = str(record.seq)
                    elif '-E' in i and '-R' not in i:#FBgn0036133|FBpp0075940|FBtr0076210-E2|Q9VTF5, i = FBtr0076210-E2
                        ref_pep_dic[i] = str(record.seq)
            else:
                ref_pep_dic[seq_name_info[-1]] = str(record.seq)

     
    with open(nucl_fasta_out_full,'w') as nucl_fh_out_full,open(pro_fasta_out_full,'w') as pro_fh_out_full:
        name = {}
        new_name = {}
        with open(file_in, 'r') as fh_in:
            count = 1
            for lines in fh_in.readlines()[1:]:
                line = lines.split('\t')
                msa_name = str(line[1]).strip()
                name[msa_name] = '\t'.join([msa_name, str(line[3]), str(line[4])])
                new_name[msa_name] = str(taxa_group)[0:3].capitalize() + str(count)+ ' '+ str(line[0]).strip()
                count += 1

        with open(file_in, 'r') as fh_in:
            #for nucl files
            nucl_out_dic = {}
            count_nucl = 1
            for nucl_file in nucl_files:
                nucl_name = str(os.path.split(nucl_file)[1]).replace('Re_align_', '').replace('_nucl.fasta', '')
                if nucl_name in name.keys():
                    n = 0
                    for records in SeqIO.parse(open(nucl_file), 'fasta'):
                        n += 1
                        if n == 1:
                            names = name[nucl_name].split('\t')
                            new_names = new_name[nucl_name]
                            start = (int(names[1]) - 1) * 3
                            end = (int(names[2]) - 1 + (num - 1) + 1) * 3
                            seq = str(records.seq[start:end]).replace('-', '')
                            id = str(new_names).split(' ')[0]
                            record_full = '>' + id +  '\n' + ref_cds_dic[names[0]] + '\n'
                            nucl_out_dic[count_nucl] = record_full
                            count_nucl += 1
            for i in range(1,count_nucl):
                    nucl_fh_out_full.write('>' + str(taxa_group)[0:3].capitalize() + str(i) +  '\n' + nucl_out_dic[i].split('\n')[-2] + '\n')                       

            #for pep files
            pro_out_dic = {}
            count = 1
            for pro_file in  pro_files:
                
                pro_name = str(os.path.split(pro_file)[1]).replace('Re_align_', '').replace('_pep.fasta', '')                
                if pro_name in name.keys():
                    n = 0                    
                    for records in SeqIO.parse(open(pro_file), 'fasta'):
                        n += 1
                        if n == 1:
                            names = name[pro_name].split('\t')
                            new_names = new_name[pro_name]
                            start = (int(names[1]) - 1)
                            end = (int(names[2]) - 1 + (num - 1) + 1)
                            seq = str(records.seq[start:end]).replace('-', '')
                            id = str(new_names).split(' ')[0]                            
                            record_full = '>' + id +  '\n' + ref_pep_dic[names[0]] + '\n'
                            pro_out_dic[count] = record_full
                            count += 1
            for i in range(1,count):
                    pro_fh_out_full.write('>' + str(taxa_group)[0:3].capitalize() + str(i) +  '\n' + pro_out_dic[i].split('\n')[-2] + '\n')

      

def Pipeline(PId, max, min, File, dir, taxa_group):
    Simp_dir_8 = '5.5.Simplifed_information_for_nested-PCR_primers'      
    Result_path = os.path.join(os.getcwd(), '6.Designed_nested-PCR_primer_set_of_NPCLs_PId-' + str(PId))
    Refer_fasta = os.path.join(os.getcwd(), '6.Reference_sequences_for_target_regions')

    Simp_path_8 = os.path.join(os.getcwd(), Simp_dir_8)
    if not os.path.exists(Simp_path_8):
        os.mkdir(Simp_path_8)        

    if not os.path.exists(Result_path):
        os.mkdir(Result_path)

    if not os.path.exists(Refer_fasta):
        os.mkdir(Refer_fasta)

    #output 
    txt_files_8AA = glob.glob(os.path.join(os.getcwd(), '5.4.8AA_F1R1_Primers', '*.txt'))
    for txt_file_8AA in txt_files_8AA:
        Simplify_txt_for_nested_primers(txt_file_8AA, int(max), int(min), Simp_path_8)

    #Creat output folder
    primer_pair_num = Stat_primer_set_feature(Simp_dir_8)[0]
    mean_length = Stat_primer_set_feature(Simp_dir_8)[1]

    Result_file_8AA_1 = os.path.join(Result_path, 'Highest-scoring_nested-PCR_primer_set_PId' + str(PId) + '_num' + str(primer_pair_num) + '_mlen' + str(mean_length) +'.xls')
    Result_file_8AA_candidate_top10 = os.path.join(Result_path, 'Candidate_nested-PCR_primer_set_PId' + str(PId) + '_num' + str(primer_pair_num) + '_mlen' + str(mean_length) +'.xls')
    Result_file_8AA_1_syn = os.path.join(Result_path, 'Synthezised_highest-scoring_nested-PCR_primer_PId' + str(PId) + '_num' + str(primer_pair_num) + '_mlen' + str(mean_length) +'.xls')
    Result_fasta_8AA_1 = os.path.join(Refer_fasta, 'Reference')

    #output    
    Sorted_Nested_primers(Result_file_8AA_1, Result_file_8AA_candidate_top10, Simp_path_8, float(PId), File, taxa_group)
    Synthezise_nested_primer_result(Result_file_8AA_1, Result_file_8AA_1_syn, taxa_group)            
    Extract_reference_sequences(Result_file_8AA_1, Result_fasta_8AA_1, dir, 8,taxa_group)

    return 'The work of designing NPC universal primer sets is complete! '   


if __name__ == '__main__':
    dirpath_o = os.path.abspath(os.path.join(os.path.dirname('__file__'), os.path.pardir))
    os.chdir(dirpath_o)
    parser = argparse.ArgumentParser(description='''python _select_primer_pairs_7.py -F species_list.txt -D 4.Candidate_peptide_iden0.5_MSAs_for_primer_design -TG [The name of target group]''')

    parser.add_argument('--dir', '-D', help='Input a directory that contains all candidate peptide MSAs (for example: 4.Candidate_peptide_iden0.5_MSAs_for_primer_design)', required=True)
    parser.add_argument('--max', '-M', help='The maximal length of target amplicons (amino acid length, from F2 to R2)', default=int(700))
    parser.add_argument('--min', '-N', help='The minimal length of target amplicons (amino acid length, from F2 to R2)', default=int(100))
    parser.add_argument('--PId', '-PId', help='The weighting ratio of PCR score and Information score when scoring highest-scoring nested-PCR primers derived from different MSAs', default=float(1.0))
    parser.add_argument('--File', '-F', help='A text table containing information on the reference species, ingroup species, and outgroup species (filename: species_list.txt)', required=True)
    parser.add_argument('--TaxaGroup', '-TG', help='Name of the target group (for example: Vertebrata, Bivalvia and so on)', required=True)
    args = parser.parse_args()
    print args
    print Pipeline(args.PId, args.max, args.min, args.File, args.dir, args.TaxaGroup)
