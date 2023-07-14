#!/usr/bin/python
# -*- coding: UTF-8 -*-
import glob, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from collections import Counter
import multiprocessing
import shutil
import argparse
import re
linesep = os.linesep


Amino_acid_degeneracy_dict = dict(zip('GAVIPFYWTCMNQDEKH-XLSRBJZ', [4, 4, 4, 3, 4, 2, 2, 1, 4, 2, 1, 2, 2, 2, 2, 2, 2, 0, 0, 8, 16, 8, 2, 2, 8]))
Amino_acid_codon_dict = dict(zip('GAVIPFYWTCMNQDEKH-XLSRBJZ', ['GGN', 'GCN', 'GTN', 'ATH', 'CCN', 'TTY', 'TAY', 'TGG', 'ACN', 'TGY',
										    'ATG', 'AAY', 'CAR', 'GAY', 'GAR', 'AAR', 'CAY', '---', 'XXX',
										    'YTN', 'WSN', 'MGN', 'AAY','GAR','YTN']))
Base_degeneracy_dict = dict(zip('-XATCGRYMKSWBVHDN', [0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4]))

def Calculate_column_identity(s):
	m = 0
	s = s.upper()
	length = len(s)
	majority = 0
	AA = ''
	gap_num = s.count('-')
	occ_num = s.count('X')

	
	all_per = float(gap_num + occ_num) / len(s) * 100
	if all_per > 50:
		AA = '-'
		return 0, AA, m
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
					return float(majority) * 1.0 / length, AA, m
		else:
		       return 0, AA, m 


def Identify_two_conserved_columns_for_F2R2(n, align, num_sp, ConsiderOutgroups):
        
	point = 0
	out_num = 0
	m = n
	aa = ''
	rm_aa_list_count = 0
	
	for records in align:
		if "Outgroup_" in str(records.id):
			out_num += 1
		else:
			pass			
	ins_num = num_sp - out_num


        if ConsiderOutgroups == 1:	
                for j in range(0, 2):
                        out_col = str(align[num_sp - out_num:, m + j])
                        ins_col = str(align[:num_sp - out_num, m + j])


                        if '-' not in ins_col:
                                first_AA = ins_col[0].upper()
                                AA_num = ins_col.count(first_AA)
                                all_length = int(ins_col.count('X')) + int(AA_num)

                                if all_length == ins_num:
                                        AA_out = out_col.count(first_AA)
                                        if AA_out >= int(0.50 * out_num):                         
                                                n += 1
                                                aa += first_AA
                                                point += 0.5
                                        else:
                                                n += 1
                                                point = 0
                                                break
                                else:
                                        n += 1
                                        point = 0
                                        break
                        else:
                                n += 1
                                point = 0
                                break
        else:
                for j in range(0, 2):
                        ins_col = str(align[:ins_num, m + j])

                        if '-' not in ins_col:
                                first_AA = ins_col[0].upper()
                                AA_num = ins_col.count(first_AA)
                                all_length = int(ins_col.count('X')) + int(AA_num)

                                if all_length == ins_num:                       
                                        n += 1
                                        aa += first_AA
                                        point += 0.5

                                else:
                                        n += 1
                                        point = 0
                                        break
                        else:
                                n += 1
                                point = 0
                                break                

	if point == 1: 
		LSR_list = ['L', 'S', 'R']
		WM_list = ['W', 'M']
		all_list = ['L', 'S', 'R', 'W', 'M']
		L = aa.count('L')
		S = aa.count('S')
		R = aa.count('R')
		M = aa.count('M')
		W = aa.count('W')
		lsr = L + S + R
		mw = M + W
		
		if lsr == 0 and mw <= 1:
			return m, point, aa

		else:
			return m, point - 1, None
	else:
		return m, point, None


def Identify_two_conserved_columns_for_F1R1(n, align, num_sp, ConsiderOutgroups):
	point = 0
	out_num = 0
	m = n
	aa = ''
	rm_aa_list_count = 0
	
	for records in align:
		if "Outgroup_" in str(records.id):
			out_num += 1
		else:
			pass			
	ins_num = num_sp - out_num

        if ConsiderOutgroups == 1:                
                for j in range(0, 2):

                        out_col = str(align[num_sp - out_num:, m + j])
                        ins_col = str(align[:num_sp - out_num, m + j])


                        if '-' not in ins_col:
                                first_AA = ins_col[0].upper()
                                AA_num = ins_col.count(first_AA)
                                all_length = int(ins_col.count('X')) + int(AA_num)

                                if all_length == ins_num:
                                        AA_out = out_col.count(first_AA)
                                        if AA_out >= int(0.50 * out_num):                              
                                                n += 1
                                                aa += first_AA
                                                point += 0.5
                                        else:
                                                n += 1
                                                point = 0
                                                break
                                else:
                                        n += 1
                                        point = 0
                                        break
                        else:
                                n += 1
                                point = 0
                                break
        else:

                for j in range(0, 2):

                        ins_col = str(align[:num_sp - out_num, m + j])

                        if '-' not in ins_col:
                                first_AA = ins_col[0].upper()
                                AA_num = ins_col.count(first_AA)
                                all_length = int(ins_col.count('X')) + int(AA_num)

                                if all_length == ins_num:                            
                                        n += 1
                                        aa += first_AA
                                        point += 0.5

                                else:
                                        n += 1
                                        point = 0
                                        break
                        else:
                                n += 1
                                point = 0
                                break                

	if point == 1: 
		LSR_list = ['L', 'S', 'R']
		WM_list = ['W', 'M']
		all_list = ['L', 'S', 'R', 'W', 'M']
		L = aa.count('L')
		S = aa.count('S')
		R = aa.count('R')
		M = aa.count('M')
		W = aa.count('W')
		lsr = L + S + R
		mw = M + W
			  
		if aa[1] not in LSR_list and lsr + mw <= 1:
			return m, point, aa

		else:
			return m, point - 1, None
	else:
		return m, point, None

def Get_primer_block_consensus_seq(align, start, thr_s, AA, ConsiderOutgroups):
	out_num = 0
	in_seqs_f = []
	out_seqs_f = []
	in_seqs_r = []	
	out_seqs_r = []
	con_seq_f = ''
	con_seq_r = ''
	f = 0
	out_f = 0
	r = 0
	out_r = 0
	num_sp = len(align)

	for records in align:
		if "Outgroup_" in str(records.id):
			out_num += 1
		else:
			pass


	if ConsiderOutgroups == 1:		
                for i in range(AA - 2, 0, -1):
                        in_seq_f = str(align[:num_sp - out_num, start - i])
                        in_seqs_f.append(in_seq_f)
                        
                        out_seq_f = str(align[num_sp - out_num:, start - i])
                        out_seqs_f.append(out_seq_f)

                for j in range(2, AA):

                        in_seq_r = str(align[:num_sp - out_num, start + j])	
                        in_seqs_r.append(in_seq_r)
                        
                        out_seq_r = str(align[num_sp - out_num:, start + j])		
                        out_seqs_r.append(out_seq_r)		
                
                for k in range(0, AA - 2):
                        
                        iden_ins_f = Calculate_column_identity(in_seqs_f[k])
                        iden_out_f = Calculate_column_identity(out_seqs_f[k])
                        iden_ins_r = Calculate_column_identity(in_seqs_r[k])
                        iden_out_r = Calculate_column_identity(out_seqs_r[k])
                        con_seq_f += iden_ins_f[1]
                        con_seq_r += iden_ins_r[1]
                        
                con_seq_f = con_seq_f + thr_s
                con_seq_r = thr_s + con_seq_r

                ingroup_colunm_iden_f = 0
                outgroup_colunm_iden_f = 0
                all_colunm_iden_f = 0
                
                ingroup_colunm_iden_r = 0
                outgroup_colunm_iden_r = 0
                all_colunm_iden_r = 0
                
                in_seqs_all_f = []
                in_seqs_all_r = []
                
                all_seqs_all_f = []
                all_seqs_all_r = []
                
                out_seqs_all_f = []
                out_seqs_all_r = []

                #extract consencus primers	
                for m in range(AA, 0, -1):
                        in_seqs_a_f = str(align[:num_sp - out_num, start - m + 2])
                        in_seqs_all_f.append(in_seqs_a_f)

                        all_seqs_a_f = str(align[:, start - m + 2])
                        all_seqs_all_f.append(all_seqs_a_f)

                        out_seqs_a_f = str(align[num_sp - out_num:, start - m + 2])
                        out_seqs_all_f.append(out_seqs_a_f)

                for p in range(0, AA):
                        in_seqs_a_r = str(align[:num_sp - out_num, start + p])	
                        in_seqs_all_r.append(in_seqs_a_r)

                        all_seqs_a_r = str(align[:, start + p])
                        all_seqs_all_r.append(all_seqs_a_r)		

                        out_seqs_a_r = str(align[num_sp - out_num:, start + p])	
                        out_seqs_all_r.append(out_seqs_a_r)		


                #calculate block identity (ingroup,outgroup,and all species)		
                for l in range(0, AA):

                        ingroup_colunm_iden_f += Calculate_column_identity(in_seqs_all_f[l])[0]
                        outgroup_colunm_iden_f += Calculate_column_identity(out_seqs_all_f[l])[0]
                        all_colunm_iden_f += Calculate_column_identity(all_seqs_all_f[l])[0]

                        ingroup_colunm_iden_r += Calculate_column_identity(in_seqs_all_r[l])[0]                
                        outgroup_colunm_iden_r += Calculate_column_identity(out_seqs_all_r[l])[0]
                        all_colunm_iden_r += Calculate_column_identity(all_seqs_all_r[l])[0]    
                
                f = round(float(ingroup_colunm_iden_f)/float(AA),2)
                out_f = round(float(outgroup_colunm_iden_f)/float(AA),2)
                all_f = round(float(all_colunm_iden_f)/float(AA),2)
                
                r = round(float(ingroup_colunm_iden_r)/float(AA),2)
                out_r = round(float(outgroup_colunm_iden_r)/float(AA),2)
                all_r = round(float(all_colunm_iden_r)/float(AA),2)
                
                return con_seq_f, str(f), str(out_f), str(all_f), con_seq_r, str(r), str(out_r), str(all_r), out_num

        else:
                for i in range(AA - 2, 0, -1):
                        in_seq_f = str(align[:num_sp - out_num, start - i])
                        in_seqs_f.append(in_seq_f)

                for j in range(2, AA):

                        in_seq_r = str(align[:num_sp - out_num, start + j])	
                        in_seqs_r.append(in_seq_r)
		                
                for k in range(0, AA - 2):
                        
                        iden_ins_f = Calculate_column_identity(in_seqs_f[k])
                        iden_ins_r = Calculate_column_identity(in_seqs_r[k])

                        con_seq_f += iden_ins_f[1]
                        con_seq_r += iden_ins_r[1]
                        
                con_seq_f = con_seq_f + thr_s
                con_seq_r = thr_s + con_seq_r

                ingroup_colunm_iden_f = 0
                ingroup_colunm_iden_r = 0
 
                
                in_seqs_all_f = []
                in_seqs_all_r = []

                #extract consencus primers	
                for m in range(AA, 0, -1):
                        in_seqs_a_f = str(align[:num_sp - out_num, start - m + 2])
                        in_seqs_all_f.append(in_seqs_a_f)

                for p in range(0, AA):
                        in_seqs_a_r = str(align[:num_sp - out_num, start + p])	
                        in_seqs_all_r.append(in_seqs_a_r)

                #calculate block identity (ingroup,outgroup,and all species)		
                for l in range(0, AA):

                        ingroup_colunm_iden_f += Calculate_column_identity(in_seqs_all_f[l])[0]
                        ingroup_colunm_iden_r += Calculate_column_identity(in_seqs_all_r[l])[0]                
   
                
                f = round(float(ingroup_colunm_iden_f)/float(AA),2)                
                r = round(float(ingroup_colunm_iden_r)/float(AA),2)
                
                return con_seq_f, str(f), str(0), str(f), con_seq_r, str(r), str(0), str(r), out_num
                

def Calculate_degeneracy_for_two_conserved_columns(seq_AA, pos):
	de_AA = 1
	de_nucl_f = 1
	de_nucl_r = 1
	nucl = ''
	out_num = 0
	for i in range(0, len(seq_AA)):
		aa = seq_AA[i]
		de_nucl_f *= Amino_acid_degeneracy_dict[aa]
		nucl += Amino_acid_codon_dict[aa]
		
	nucl_r = Seq(nucl).reverse_complement()
	de_nucl_r = de_nucl_f

	return nucl[:-1], str(nucl_r), de_nucl_f, de_nucl_r #'nucl[:-1]' means that removing the 3'end of base from forward nucleotide primers


	

def Calculate_degeneracy_for_consensus_primer(sim_list, pos, degen_list, AA_num):
	out_num = sim_list[-1] # ('TESLWRAN', '0.97', '1.0', '0.99', 'ANATRNAS', '1.0', '0.95', '0.99', 7)
	AA3_nucl_seq_f = degen_list[0] #('GCNAA', 'TTNGC', 8, 8)
	AA3_nucl_seq_r = degen_list[1]
	AA3_nucl_f = degen_list[2]
	AA3_nucl_r = degen_list[3]

	#calculate forward primer degeneracy
	de_nucl_f = 1
	degen_f = 1
	degen_r = 1
	nucl_f = ''
	con_seq_f = sim_list[0][:AA_num - 2] #con_seq_f = TESLWR
	for j in range(0, len(con_seq_f)):
		aa = con_seq_f[j]
		nucl_f += Amino_acid_codon_dict[aa]
				
	seq_f = str(nucl_f + AA3_nucl_seq_f).replace('=','')
	for base in list(seq_f):
			degen_f *= Base_degeneracy_dict[base]
				
	#calculate reverse primer degeneracy
	de_nucl_r = 1
	nucl_r = ''
	con_seq_r = sim_list[4][2:] #con_seq_r='ATRNAS'		
	for j in range(0, len(con_seq_r)):                
		aa = con_seq_r[j]
		nucl_r += Amino_acid_codon_dict[aa]
		
	nucl_r_all =  Seq(AA3_nucl_seq_r).reverse_complement().upper() +  nucl_r	
	for base in list(nucl_r_all[:-1]):
		degen_r *= Base_degeneracy_dict[base]
			
	nucl_s_r = Seq(nucl_r).reverse_complement().upper()	
	seq_r = str(nucl_s_r[1:] + AA3_nucl_seq_r)

	return seq_f, str(degen_f), seq_r, str(degen_r)


def Calculate_MSA_identity(fasta, num, ConsiderOutgroups):
	out_num = 0
	sim_all_cal = 0
	sim_in_cal = 0
	sim_out_cal = 0
	length_all = 0
	length_ins = 0
	length_out = 0
	align_pep = AlignIO.read(fasta, 'fasta')
	species_num = len(align_pep)

	for records in align_pep:
		if "Outgroup_" in str(records.id):
			out_num += 1
		else:
			pass
	if ConsiderOutgroups == 1:		
                for i in range(0, num):
                        all = align_pep[:, i]
                        out = align_pep[species_num - out_num:, i]
                        ins = str(align_pep[:species_num - out_num, i])

                        sim_1 = Calculate_column_identity(all)[0]
                        if sim_1 != 0:
                                length_all += 1
                                sim_all_cal += sim_1
                        sim_2 = Calculate_column_identity(ins)[0]
                        if sim_2 != 0:
                                length_ins += 1
                                sim_in_cal += sim_2

                        sim_3 = Calculate_column_identity(out)[0]
                        if sim_3 != 0:
                                length_out += 1
                                sim_out_cal += sim_3

                return round(float(sim_all_cal) / length_all, 4), round(float(sim_in_cal) / length_ins, 4), round(
                        float(sim_out_cal) / length_out, 4)
        else:
                for i in range(0, num):
                        all = align_pep[:, i]
                        ins = str(align_pep[:species_num - out_num, i])

                        sim_1 = Calculate_column_identity(all)[0]
                        if sim_1 != 0:
                                length_all += 1
                                sim_all_cal += sim_1
                        sim_2 = Calculate_column_identity(ins)[0]
                        if sim_2 != 0:
                                length_ins += 1
                                sim_in_cal += sim_2

                return round(float(sim_all_cal) / length_all, 4), round(float(sim_in_cal) / length_ins, 4), 0
                


def Calculate_variable_columns(start, end, pep_fasta):
	fasta = AlignIO.read(pep_fasta, 'fasta')
	num_sp = len(fasta)
	out_num = 0
	n = 0
	m = 0
	for records in fasta:
		if "Outgroup_" in str(records.id):
			out_num += 1
		else:
			pass
		
	for i in range(start, end):
		ins = str(fasta[:num_sp - out_num,i])
		first_AA = ins[0]
		all_len = len(ins)
		len_1 = all_len - 1
		num_all = ins.count(first_AA) + ins.count('X')

		if first_AA != '-':
			m += 1
			if num_all != all_len:
				n += 1
			else:
				pass
			
		elif num_all != len_1:
			m += 1
			n += 1
		else:
			pass
	return n, m

def Calculate_information_score(start, end, pep_fasta, ingroup_num):
	align_pep = AlignIO.read(pep_fasta, 'fasta')
	information_score = 0
	num_seq = len(align_pep)	

	for i in range(start, end):
		colunm = str(align_pep[0: ingroup_num + 1, i])#for ingroups only
		X_num = colunm.count('X')
		gap_num = colunm.count('-')

		column_identity = float(Calculate_column_identity(colunm)[0])
		if column_identity == 0:
			colunm_diversity = 0
			information_score += colunm_diversity
		else:
			colunm_diversity = 1 - column_identity
			information_score += colunm_diversity
	return information_score

	
def Assess_the_consecutive_of_amino_acid_in_primer(primer_seq):
	primer_seq_list = list(primer_seq)
	repeat_count = 0
	for AA in primer_seq_list:
		
		if re.search(r'(.)\1\1', primer_seq):
		    repeat_count += 1
		    break
		else:
		    pass
	return repeat_count


def Assess_the_diversity_of_amino_acid_in_primer(F,R):#Assessing amino acid Diversity
##    F = EDDKLTFH
##    R = HCEITIKY
        F_counter = Counter(str(F))
        R_counter = Counter(str(R))

        F_bio_score = 0
        R_bio_score = 0

        for aa in F_counter:
                aa_bio_score = int(F_counter[aa]) * (int(F_counter[aa]) - 1)
                F_bio_score += aa_bio_score

        for aa in R_counter:
                aa_bio_score = int(R_counter[aa]) * (int(R_counter[aa]) - 1)
                R_bio_score += aa_bio_score

        FR_bio_score = F_bio_score + R_bio_score

        return FR_bio_score

def Calculate_primer_degeneracy_score(F,R):
##    F = EDDKLTFH
##    R = HCEITIKY

    primer_length = len(list(F))
    if primer_length == 8:            
            weight_ratio_min = 0.3
            weight_ratio_max = 1
            Tolerances = (weight_ratio_max - weight_ratio_min)/(primer_length - 1)

            # generate a weight ratio list for 8aa primers
            weight_ratio_list = []
            for i in range(1,primer_length + 1):
                    weight_ratio_for_aa = 0.3 + Tolerances * (int(i)-1)
                    weight_ratio_list.append(round(weight_ratio_for_aa,4))
            weight_ratio_list_for_R = weight_ratio_list[::-1]

            # calculate PCR score for 8aa F and R primers            
            degen_score_for_F_primer = 0
            degen_score_for_R_primer = 0
            for i in range (0, primer_length):
                    F_aa = F[i]
                    degen_score_for_F_primer += weight_ratio_list[i] * Amino_acid_degeneracy_dict[F_aa]
                    R_aa = R[i]
                    degen_score_for_R_primer += weight_ratio_list_for_R[i] * Amino_acid_degeneracy_dict[R_aa]
            Primer_degeneracy_score =  degen_score_for_F_primer + degen_score_for_R_primer

    if primer_length == 7:            
            weight_ratio_min = 0.4
            weight_ratio_max = 1
            Tolerances = (weight_ratio_max - weight_ratio_min)/(primer_length - 1)
            weight_ratio_list = []

            # generate a weight ratio list for 7aa primers
            for i in range(1,primer_length + 1):
                    weight_ratio_for_aa = 0.4 + Tolerances * (int(i)-1)
                    weight_ratio_list.append(round(weight_ratio_for_aa,4))
            weight_ratio_list_for_R = weight_ratio_list[::-1]

            # calculate PCR score for 7aa F and R primers
            degen_score_for_F_primer = 0
            degen_score_for_R_primer = 0
            for i in range (0, primer_length):
                    F_aa = F[i]
                    degen_score_for_F_primer += weight_ratio_list[i] * Amino_acid_degeneracy_dict[F_aa]

                    R_aa = R[i]
                    degen_score_for_R_primer += weight_ratio_list_for_R[i] * Amino_acid_degeneracy_dict[R_aa]
            Primer_degeneracy_score =  degen_score_for_F_primer + degen_score_for_R_primer
    return Primer_degeneracy_score

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

def Generate_primer_information(pep_fasta, AA_num, conserved_identity, Dir, ConsiderOutgroups):
	              
        path_out = os.path.join(os.getcwd(), '5.2.Information_for_' + str(AA_num) + 'AA_Primers')
        if not os.path.exists(path_out):
                os.mkdir(path_out)

	txt_AA = str(AA_num) + 'AA.txt'
	align_name = str(os.path.split(pep_fasta)[1]).replace('Re_align_', '').replace('pep.fasta', txt_AA)
	out_file = os.path.join(path_out, align_name)

	if not os.path.exists(out_file):
		list_line = ['2AA Start Position in the Peptide Alignment',
			     '2AA Peptide Sequence in the Forward Primer (5\' to 3\')',
			     '2AA Nucleotide Sequence in the Forward Primer (5\' to 3\')',
			     '2AA Nucleotide Sequence in the Reverse Primer (5\' to 3\', reverse complement)',
			     'Consensus Peptide Sequence of the Forward Primer (5\' to 3\')',
			     'Ingroup mean identity in Forward Primer Blocks',
			     'Outgroup mean identity in Forward Primer Blocks',
			     'All species mean identity in Forward Primer Blocks',
			     'Consensus Peptide Sequence of the Reverse Primer (5\' to 3\', no reverse complement)',
			     'Ingroup mean identity in Reverse Primer Blocks',
			     'Outgroup mean identity in Reverse Primer Blocks',
			     'All species mean identity in Forward Primer Blocks',
			     'Degeneracy of 2AA in the Forward Primer',
			     'Degeneracy of 2AA in the Reverse Primer',
			     'Nucleotide Sequence of the Forward Primer (5\' to 3\', after removing the last base)',
			     'Degeneracy of the Forward Primer',
			     'Nucleotide Sequence of the Reverse Primer (5\' to 3\',reverse complement, after removing the first base)',
			     'Degeneracy of the Reverse Primer',
			     'No. of Ingroups in Alignment']
		line_1 = '\t'.join(list_line)
		

		align_pep = AlignIO.read(pep_fasta, 'fasta')
		seq_num = len(align_pep)##species number
		seq_len = len(align_pep[0].seq)
		sim_list = Calculate_MSA_identity(pep_fasta, seq_len, ConsiderOutgroups)#MSA similarity
		
		if sim_list[0] >= conserved_identity:
			return 'Filter', None #do not design primers for highly conserved MSAs
		else:
			thr_AA_n = []
			thr_AA_s = []
			i = AA_num - 2
			j = seq_len - AA_num + 1
			
			while (i < j): 
				re_li = Identify_two_conserved_columns_for_F2R2(i, align_pep, seq_num, ConsiderOutgroups)
				i = re_li[0]
				if re_li[1] == 1:
					thr_AA_n.append(i)
					thr_AA_s.append(re_li[2].upper())
				i += 1

			thr_AA_length = len(thr_AA_n)			
			if thr_AA_length != 0:
				with open(out_file, 'w') as fh:
					
					fh.write(line_1 + linesep)					
					for n in range(0, thr_AA_length):						
						t_AA = thr_AA_s[n]						
						position = thr_AA_n[n]

						if position - AA_num > 0 and (position + AA_num) < seq_len:

							aa3_degen = Calculate_degeneracy_for_two_conserved_columns(t_AA, position)				
							aa_consensus_list = Get_primer_block_consensus_seq(align_pep, position, t_AA, AA_num, ConsiderOutgroups)
							Ingroup_num = seq_num - int(aa_consensus_list[-1])
							nucl_degen = Calculate_degeneracy_for_consensus_primer(aa_consensus_list, position, aa3_degen, AA_num)

							
							if '-' not in aa_consensus_list[0] and 'X' not in aa_consensus_list[0] and \
							   '-' not in aa_consensus_list[4] and 'X' not in aa_consensus_list[4]: 
								line_2 = str(int(position) + 1)  + '\t' + t_AA + '\t' + aa3_degen[0] + '\t' + aa3_degen[1] + '\t' \
										 + '\t'.join(aa_consensus_list[:-1]) + '\t' + str(aa3_degen[2]) + '\t' \
										 + str(aa3_degen[3]) + '\t' + '\t'.join(nucl_degen) + '\t' + str(Ingroup_num)					
								fh.write(str(line_2) + linesep)
								
				return pep_fasta, out_file
			else:
				return 'Filter', None
	else:
		return pep_fasta, out_file



def Select_F2R2(pep_file, txt, AA_num, Dir, all_deg, inthred2, allthred2, max, min):
	              
        pick_AA = '5.1.Selected_Peptide_MSAs_for_' + str(AA_num) + 'AA_Primer_design'
        picked_AA = os.path.join(os.getcwd(), pick_AA)
        if not os.path.exists(picked_AA):
                os.mkdir(picked_AA)

        path_out = os.path.join(os.getcwd(), '5.3.' + str(AA_num) + 'AA_F2R2_Primers')
        if not os.path.exists(path_out):
                os.mkdir(path_out)


	print 'Searching ' + str(AA_num) +'AA F2R2 primer pairs for ' + str(os.path.basename(txt).replace('.txt','')) + ' peptide MSA!'
	
	
	primer_file = os.path.join(path_out, os.path.split(txt)[1].replace('.txt', '_primer.txt'))

	if not os.path.exists(primer_file):                        
		line_1 = ('\t').join(['F2 Primer Start Position in the Peptide Alignment',
				      'R2 Primer Start Position in the Peptide Alignment',
				      'Nucleotide Sequence of F2 Primer(5\' to 3\')',
				      'Peptide Sequence of F2 Primer(5\' to 3\')',
				      'Degeneracy of 2AA in F2 primer',
				      'Degeneracy of F2 Primer',
				      'Nucleotide Sequence of R2 primer(5\' to 3\',reverse complement)',
				      'Peptide Sequence of R2 Primer(5\' to 3\')',
				      'Degeneracy of 2AA in R2 Primer',
				      'Degeneracy of R2 Primer',
				      'No. of Variable colunms in the Peptide Alignment for Target Regions (from F2 to R2)',
				      'Reference Length of PCR',
				      'No. of Ingroup Taxa in Alignment(include reference species)',
				      '% Variable Colunms for F2R2',
				      'Information Score for F2R2',
				      'PCR_score_for_F2R2',
                                      'Identity_score_for_F2R2',
				      'Amplify_score_for_F2R2'])
		
		line_2 = []
		with open(txt, 'r') as fh_in:
			lines = fh_in.readlines()
			length = len(lines)

			n = 0
			F2 = []
			R2 = []
			for i in range(1, length):
				lisp_f = lines[i].split('\t')

				F2_2aa_start = str(lisp_f[0])
				F2_2aa_pep = str(lisp_f[1])
				R2_2aa_nucl = str(lisp_f[2])
				F2_primer_pep_con = str(lisp_f[4])
				F2_primer_nucl_con = str(lisp_f[14])
				F2_2aa_degen = str(lisp_f[12])              
				F2_primer_degen = str(lisp_f[15]) 
				F2_in_mean_iden = str(lisp_f[5]) 
				F2_out_mean_iden = str(lisp_f[6])
				F2_all_mean_iden = str(lisp_f[7]) 
				F2_Ingroup_num = str(lisp_f[-1].strip())
				F2_line = '\t'.join([F2_2aa_start,F2_2aa_pep, R2_2aa_nucl,F2_primer_pep_con,F2_primer_nucl_con,F2_2aa_degen,F2_primer_degen,F2_in_mean_iden,F2_out_mean_iden,F2_all_mean_iden,F2_Ingroup_num])

				if '-' not in F2_primer_pep_con and 'X' not in F2_primer_pep_con and Assess_the_consecutive_of_amino_acid_in_primer(F2_primer_pep_con) == 0:# primers do not contain gap, X and repeat AA                                     
					F2.append(F2_line)


				R2_2aa_start = str(lisp_f[0])
				R2_2aa_pep = str(lisp_f[1])
				R2_2aa_nucl = str(lisp_f[3])
				R2_primer_pep_con = str(lisp_f[8])
				R2_primer_nucl_con = str(lisp_f[16])
				R2_2aa_degen = str(lisp_f[13])
				R2_primer_degen = str(lisp_f[17])
				R2_in_mean_iden = str(lisp_f[9])
				R2_out_mean_iden = str(lisp_f[10])
				R2_all_mean_iden = str(lisp_f[11])
				R2_Ingroup_num = str(lisp_f[-1].strip())
				R2_line = '\t'.join([R2_2aa_start,R2_2aa_pep, R2_2aa_nucl,R2_primer_pep_con,R2_primer_nucl_con,R2_2aa_degen,R2_primer_degen,R2_in_mean_iden,R2_out_mean_iden,R2_all_mean_iden,R2_Ingroup_num])                                

				if '-' not in R2_primer_pep_con and 'X' not in R2_primer_pep_con and Assess_the_consecutive_of_amino_acid_in_primer(R2_primer_pep_con) == 0:# primers do not contain gap, X and repeat AA  
					R2.append(R2_line)
					

			F2_len = len(F2)
			R2_len = len(R2)

			if F2_len > 0 and R2_len > 0:
				for i in range(0,F2_len):                                        
					F2_line = F2[i].strip().split('\t')
					
					F2_2aa_start = int(F2_line[0]) - 1                                        
					F2_2aa_pep = str(F2_line[1])
					R2_2aa_nucl = str(F2_line[2])
					F2_primer_pep_con = str(F2_line[3])
					F2_primer_nucl_con = str(F2_line[4])                                        
					F2_2aa_degen = int(F2_line[5])
					F2_primer_degen = int(F2_line[6])
					F2_in_mean_iden = float(F2_line[7])
					F2_all_mean_iden = float(F2_line[9])
					F2_Ingroup_num = int(F2_line[10])

					for j in range(0, R2_len):
						R2_line = R2[j].strip().split('\t')
						
						R2_2aa_start = int(R2_line[0]) - 1                                                
						R2_2aa_pep = str(R2_line[1])
						R2_2aa_nucl = str(R2_line[2])
						R2_primer_pep_con = str(R2_line[3])
						R2_primer_nucl_con = str(R2_line[4])                                                
						R2_2aa_degen = int(R2_line[5])
						R2_primer_degen = int(R2_line[6])
						R2_in_mean_iden = float(R2_line[7])
						R2_all_mean_iden = float(R2_line[9])
						R2_Ingroup_num = int(R2_line[10])

																								
						if F2_2aa_start != R2_2aa_start: #F and R primers do not derive from the same 2aa block
							if F2_2aa_degen <= 8 and F2_primer_degen <= int(all_deg) and F2_in_mean_iden >= float(inthred2) and F2_all_mean_iden >= float(allthred2)\
							   and R2_2aa_degen <= 8 and R2_primer_degen <= int(all_deg) and R2_in_mean_iden >= float(inthred2) and R2_all_mean_iden >= float(allthred2):

								var_site = Calculate_variable_columns(F2_2aa_start, R2_2aa_start, pep_file)                                                                
								F2R2_var_col_num = str(var_site[0])
								information_score = Calculate_information_score(F2_2aa_start, R2_2aa_start, pep_file, F2_Ingroup_num)
								#count primer degen score
								F2R2_Primer_degen_score = Calculate_primer_degeneracy_score(F2_primer_pep_con, R2_primer_pep_con)
								#count primer biodiversity score
								F2R2_Primer_bio_score = Assess_the_diversity_of_amino_acid_in_primer(F2_primer_pep_con, R2_primer_pep_con)
								#count PCR score for F2R2
								PCR_score_for_F2R2 = 100 - (F2R2_Primer_degen_score + F2R2_Primer_bio_score)
								#count identity score for F2R2
                                                                F2_identity_score = F2_in_mean_iden * 100
                                                                R2_identity_score = R2_in_mean_iden * 100
                                                                Identity_score_for_F2R2 = F2_identity_score * 0.5 + R2_identity_score * 0.5

                                                                
								Amplify_score = PCR_score_for_F2R2 + Identity_score_for_F2R2                                                         

								if var_site[1] >= int(min) and var_site[1] <= int(max):

									F2_primer_start = F2_2aa_start + 1 - (AA_num - 2)
									R2_primer_end = R2_2aa_start + 1
									F2R2_length = ((R2_primer_end + AA_num) - F2_primer_start + 1) * 3#Total length of target regions (including primers)
									Ave_info = float(var_site[0]) / float(var_site[1])
									
									F2R2_line = '\t'.join([str(F2_primer_start), str(R2_primer_end), F2_primer_nucl_con, F2_primer_pep_con,str(F2_2aa_degen), str(F2_primer_degen),\
											      R2_primer_nucl_con, R2_primer_pep_con, str(R2_2aa_degen),str(R2_primer_degen),\
											      F2R2_var_col_num, str(F2R2_length), str(F2_Ingroup_num), str(Ave_info), str(information_score), str(PCR_score_for_F2R2), \
                                                                                              str(Identity_score_for_F2R2), str(Amplify_score)])
									line_2.append(F2R2_line)
									
								else:
									pass

							else:
								pass
						else:
							pass
			else:
				pass

		#sort all F2R2 pairs 
		F2R2 = {}                
		if len(line_2) > 0:
                        shutil.copy(pep_file, picked_AA)
                        with open(primer_file, 'w') as fh_out:
                                fh_out.write(line_1 + linesep)
                                
                                for i in range(0, len(line_2)):
                                        key = i
                                        F2R2[i] = [float(line_2[i].split('\t')[-1]),line_2[i]]##Sorted by Score_for_F2R2(PCR score + Identity score)
                                        
                                if F2R2:
                                        Sort_F2R2 = sorted(F2R2.items(), key=lambda x: x[1], reverse =True)

                                for i in range(0, len(line_2)):
                                        fh_out.write(Sort_F2R2[i][1][1] + linesep)
		else:
			pass                        

	else:
		print 'There exists ' + str(primer_file) + ' file.'


def Select_F1R1(primer_info, directory, F, R, AA_num, all_deg, inthred1, allthred1, ConsiderOutgroups, PIs):

	primer_txt = os.path.split(primer_info)[1]

	print 'Searching ' + str(AA_num) +'AA F1R1 primer pairs for ' + str(os.path.basename(primer_txt).replace('_primer.txt','')) + ' peptide MSA!'
	               
        path_out = os.path.join(os.getcwd(), '5.4.' + str(AA_num) +'AA_F1R1_Primers')
        if not os.path.exists(path_out):
                os.mkdir(path_out)
              
		
	path_file = primer_txt.replace('.txt', '_F1R1.txt')
	out_file = os.path.join(path_out, path_file)
	
	if not os.path.exists(out_file):
		pep = '_'.join(primer_txt.split('_')[:-2])
		pep_fasta = ''
		files_pep = [os.path.split(x)[1] for x in glob.glob(os.path.join(os.getcwd(), directory, '*.fasta'))]
		for fh_pep in files_pep:
			if pep in fh_pep:
				pep_fasta = os.path.join(os.getcwd(), directory, fh_pep)
				
		align_pep = AlignIO.read(pep_fasta, 'fasta')
		seq_num = len(align_pep)
		seq_len = align_pep.get_alignment_length()
		LSR_list = ['L','S','R']

		with open(primer_info) as fh_in:
			line_1 = fh_in.readline()
			line_2 = line_1.strip() + '\t' + '\t'.join(['F1 Primer Start Postion in the Peptide Alignment',
								    'Nucleotide Sequence of F1 Primer(5\' to 3\')',
								    'Peptide Sequence of F1 Primer(5\' to 3\')',
								    'Degeneracy of F1 Primer',
                                                                    'Identity_score_for_F1',
								    'R1 Primer Start Position in the Peptide Alignment',
								    'Nucleotide Sequence of R1 Primer(5\' to 3\')',
								    'Peptide Sequence of R1 Primer(5\' to 3\')',
								    'Degeneracy of R1 Primer',
                                                                    'Identity_score_for_R1',
								    'PCR_score_for_F1R1',
                                                                    'Identity_score_for_F1R1',
                                                                    'Amplify_score_for_F1R1',
                                                                    'Score_for_nested_PCR_primers_within_MSA'])

			
			with open(out_file, 'w') as out:
				out.write(str(line_2) + linesep)

                                Amplicon_score_list = []
                                Information_score_list = []			
                                all_nested_primers_list = []
                                all_nested_primers_dic = {}
                                nested_list = []

				for lines in fh_in.readlines():

					F1R1_list = []					
					
					line = lines.split('\t')
					thr_R_n = []
					thr_R_s = []
					R2_pos = int(line[1]) - 1					
					R1_pos = R2_pos + AA_num 									
					R1_max = seq_len - 7 
					
					j = 4					
					while (j < R and R1_pos < R1_max):
						re_li = Identify_two_conserved_columns_for_F1R1(R1_pos, align_pep, seq_num, ConsiderOutgroups)
						if re_li[1] == 1:
							thr_R_n.append(R1_pos)
							thr_R_s.append(str(re_li[2].upper()))
							
						R1_pos = re_li[0]
						j = R1_pos - R2_pos
						j += 1
						R1_pos += 1

						
					thr_F_n = []
					thr_F_s = []

					F2_pos = int(line[0]) - 1  
					F1_end = F2_pos - 2
					F1_pos = int(F2_pos - F + 5)
					
					if F1_pos < 0:                                                
						F1_pos = 5
						while (F1_pos < F1_end):
							re_li = Identify_two_conserved_columns_for_F1R1(F1_pos, align_pep, seq_num, ConsiderOutgroups)
							if re_li[1] == 1:
								thr_F_n.append(F1_pos)
								thr_F_s.append(str(re_li[2].upper()))
							F1_pos = re_li[0]
							F1_pos += 1							
															
					else:
						while (F1_pos < F1_end):                                                        
							re_li = Identify_two_conserved_columns_for_F1R1(F1_pos, align_pep, seq_num, ConsiderOutgroups)
							if re_li[1] == 1:
								thr_F_n.append(F1_pos) # - 1)
								thr_F_s.append(str(re_li[2].upper()))							
							F1_pos = re_li[0]
							F1_pos += 1
															
					F_len = len(thr_F_n)					
					R_len = len(thr_R_n)
					
					if F_len > 0 and R_len > 0:
						F1 = []
						R1 = []
						for i in range(0, F_len):
							AA_pos_F = thr_F_n[i]							
							AA_seq_F = thr_F_s[i]
							aa_consensus_list_F = Get_primer_block_consensus_seq(align_pep, AA_pos_F, AA_seq_F, 7, ConsiderOutgroups)##F1R1 primer length is 7aa

							if '-' not in aa_consensus_list_F[0] and 'X' not in aa_consensus_list_F[0] and Assess_the_consecutive_of_amino_acid_in_primer(aa_consensus_list_F[0]) == 0:####
															
								aa3_degen_F = Calculate_degeneracy_for_two_conserved_columns(AA_seq_F, AA_pos_F)#, nucl_file)
								nucl_degen_F = Calculate_degeneracy_for_consensus_primer(aa_consensus_list_F, AA_pos_F, aa3_degen_F, 7)#nucl_file,
								
								F1_degen = int(nucl_degen_F[1])
								F1_tAA_degen = int(aa3_degen_F[2])
								F1_in_mean_iden = float(aa_consensus_list_F[1])
								F1_all_mean_iden = float(aa_consensus_list_F[3])

								if F1_tAA_degen <= 32 and F1_degen <= int(all_deg) and F1_degen != 0 and F1_in_mean_iden >= float(inthred1) and F1_all_mean_iden >= float(allthred1):
									F1_record = '\t'.join([str(AA_pos_F + 1  - (7 - 2)), str(nucl_degen_F[0]), str(aa_consensus_list_F[0]), str(nucl_degen_F[1]), str(F1_in_mean_iden * 100)])                                                                        
									F1.append(str(F1_record))

									
						for j in range(0, R_len):

							AA_pos_R = thr_R_n[j]							
							AA_seq_R = str(thr_R_s[j])							
							aa_consensus_list_R = Get_primer_block_consensus_seq(align_pep, AA_pos_R, AA_seq_R, 7, ConsiderOutgroups)

														
							if '-' not in aa_consensus_list_R[4] and 'X' not in aa_consensus_list_R[4] and aa_consensus_list_R[4][0] not in LSR_list and Assess_the_consecutive_of_amino_acid_in_primer(aa_consensus_list_R[4]) == 0:

																
								aa3_degen_R = Calculate_degeneracy_for_two_conserved_columns(AA_seq_R, AA_pos_R)#, nucl_file)
								nucl_degen_R = Calculate_degeneracy_for_consensus_primer(aa_consensus_list_R, AA_pos_R, aa3_degen_R, 7)#nucl_file, 
								R1_degen = int(nucl_degen_R[3])
								R1_tAA_degen = int(aa3_degen_R[3])
								R1_in_mean_iden = float(aa_consensus_list_R[5])
								R1_out_mean_iden = float(aa_consensus_list_R[7])


								if R1_tAA_degen <= 32 and R1_degen <= int(all_deg) and R1_degen != 0 and R1_in_mean_iden >= float(inthred1) and R1_out_mean_iden >= float(allthred1):
								
									R1_record = '\t'.join([str(AA_pos_R + 1), str(nucl_degen_R[2]), str(aa_consensus_list_R[4]), str(nucl_degen_R[3]), str(R1_in_mean_iden * 100)])
									R1.append(str(R1_record))
						

						#sort all nested-pcr primer combinations                                                                
						len_F = len(F1)
						len_R = len(R1)
						
						if len_F > 0 and len_R > 0:
							for n in range(0, len_F):
								record_f = F1[n].split('\t')
##								F_degen = int(record_f[3].strip())
								for m in range(0, len_R):
									record_r = R1[m].split('\t')
									#count primer degen score
									F1R1_Primer_degen_score = Calculate_primer_degeneracy_score(record_f[2],record_r[2])
									#count primer biodiversity score
									F1R1_Primer_bio_score = Assess_the_diversity_of_amino_acid_in_primer(record_f[2],record_r[2])
									#count PCR score for F1R1
									PCR_score_for_F1R1 = 100 - (F1R1_Primer_degen_score + F1R1_Primer_bio_score)
									#count identity score for F1R1
									Identity_score_for_F1R1 = float(record_f[-1]) * 0.5  + float(record_r[-1]) * 0.5
									
									Amplify_score = PCR_score_for_F1R1 + Identity_score_for_F1R1

									line_new = lines.strip() + '\t' + F1[n] + '\t' + R1[m] + '\t' + str(PCR_score_for_F1R1) + '\t' + str(Identity_score_for_F1R1) + '\t' + str(Amplify_score)
									F1R1_list.append(str(line_new))
                                                                                                        									
						else:
                                                        pass


					
					

                                                F1R1 = {}
                                                Sort_F1R1 = {}
                                                if len(F1R1_list) > 0:
                                                        for i in range (0,len(F1R1_list)):
                                                                key = i
                                                                F1R1[i] = [float(F1R1_list[i].split('\t')[-1]),F1R1_list[i]]
                                                                                                                        
                                                        if F1R1:
                                                                Sort_F1R1 = sorted(F1R1.items(), key=lambda x: x[1], reverse =True)

                                                        nested_list.append(Sort_F1R1[0][1][1]+ linesep)##select the best F1R1 for each F2R2 primer pair                                                        
                                                
                               

                                if len(nested_list) > 0:
                                        #get score list
                                        for i in range(0,len(nested_list)):
                                        
                                                info_list = nested_list[i].split('\t')                                  
                                                Amplicon_score_for_F2R2 = float(info_list[17])
                                                Amplicon_score_for_F1R1 = float(info_list[30])
                                                Amplicon_score = Amplicon_score_for_F2R2 + Amplicon_score_for_F1R1
                                                Information_score = float(info_list[14])

                                                Amplicon_score_list.append(round(Amplicon_score,2))
                                                Information_score_list.append(round(Information_score,2))                                        

                                        #count a total score for nested-PCR primers within the same MSA                               
                                        for i in range(0,len(nested_list)):
                                                info_list = nested_list[i].strip().split('\t')                                   
                                                Amplicon_score_for_F2R2 = float(info_list[17])
                                                Amplicon_score_for_F1R1 = float(info_list[30])

                                                Amplicon_score = round(Amplicon_score_for_F2R2 + Amplicon_score_for_F1R1, 2)
                                                Information_score = round(float(info_list[14]),2)

                                                if max(Amplicon_score_list) - min(Amplicon_score_list) == 0 :
                                                        Weighted_Amplicon_score = 100
                                                else:
                                                        Weighted_Amplicon_score = Amplicon_score_weight_dic(max(Amplicon_score_list), min(Amplicon_score_list))[Amplicon_score]

                                                if max(Information_score_list) - min(Information_score_list) == 0:
                                                        Weighted_Information_score = 100
                                                else:
                                                        Weighted_Information_score = Information_score_weight_dic(max(Information_score_list),min(Information_score_list))[Information_score]
                                                
                                                Score_for_nested_PCR_primers = (float(PIs)* Weighted_Amplicon_score + Weighted_Information_score)/(float(PIs) + 1)#
                                                all_nested_primers_list.append(nested_list[i].strip() + '\t' + str(Score_for_nested_PCR_primers))

                                else:
                                        pass

                                #Sort possible nested-PCR primer pairs by PCR score and information score.
                                if len(all_nested_primers_list) > 0:
                                        for i in range (0,len(all_nested_primers_list)):
                                                key = i
                                                all_nested_primers_dic[i] = [float(all_nested_primers_list[i].split('\t')[-1]),all_nested_primers_list[i]]

                                        if all_nested_primers_dic:
                                                Sorted_nested_primers = sorted(all_nested_primers_dic.items(),key=lambda x: x[1], reverse =True)

                                        for i in range(0,len(Sorted_nested_primers)):
                                                out.write(str(Sorted_nested_primers[i][1][1]) + linesep)

	else:
		print path_file, 'file have exist!!!'

	if len(open(out_file).readlines()) == 1:
                os.system('rm ' + str(out_file))
		
def Pipeline(dir, con, fwd, rvs, all, inthred2, allthred2, inthred1, allthred1, max, min, ConsiderOutgroups, PIs):
        		
    # output 8AA F2R2 + 7AA F1R1 nested-primers
    align_files = glob.glob(os.path.join(os.getcwd(), str(dir), '*.fasta'))
    for align_file in align_files:
            print 'Generating primer information for ' + str(os.path.split(align_file)[1]) + ' file!'

            txt_8 = Generate_primer_information(align_file, 8, con, dir, ConsiderOutgroups)
            
            if txt_8[0] != 'Filter':
                    Select_F2R2(txt_8[0], txt_8[1], 8, dir, all, inthred2, allthred2, max, min)                                
            else:
                    pass
            
    primer_txts_8 = glob.glob(os.path.join(os.getcwd(), '5.3.8AA_F2R2_Primers', '*.txt'))
    for primer_txt_8 in primer_txts_8:
            Select_F1R1(primer_txt_8, dir, int(fwd), int(rvs), 8, all, inthred1, allthred1, ConsiderOutgroups, PIs)

    return 'The work of searching primers is complete!' 

if __name__ == '__main__':
	dirpath_o = os.path.abspath(os.path.join(os.path.dirname('__file__'), os.path.pardir))
	os.chdir(dirpath_o)
	parser = argparse.ArgumentParser(description='''python _design_primer_pairs_6.py -D 4.Candidate_peptide_iden0.5MSAs_for_primer_design -CO [1: Consider outgroups;0: Do not consider outgroups]''')
	
	parser.add_argument('--con', '-C', help='The conservation/identitfy of MSAs. If similarity of one MSA is larger than this value, UPrimer will not design primer for it', default=float(0.9))
	parser.add_argument('--dir', '-D', help='Input a directory that contains all candidate peptide MSAs (for example: 4.Candidate_peptide_iden0.5_MSAs_for_primer_design)', required=True)
	parser.add_argument('--fwd', '-W', help='The interval length between F1 and F2', default=int(150))
	parser.add_argument('--rvs', '-R', help='The interval length between R1 and R2', default=int(150))
	parser.add_argument('--all', '-A', help='The maximal degeneracy of primer', default=int(8192))
	parser.add_argument('--inthred2', '-IT2', help='The mean identity of ingroups (including reference species) in 7AA/8AA F2R2 primer blocks',default=float(0.85))
	parser.add_argument('--allthred2', '-AT2', help='The mean identity of all species in 7AA/8AA F2R2 primer blocks',default=float(0.75))
	parser.add_argument('--inthred1', '-IT1', help='The mean identity of ingroups (including reference species) in 7AA/8AA F1R1 primer blocks',default=float(0.80))
	parser.add_argument('--allthred1', '-AT1', help='The mean identity of all species in 7AA/8AA F1R1 primer blocks',default=float(0.70))
	parser.add_argument('--max', '-M', help='The maximal length of target amplicons (from F2 to R2)', default=int(700))
	parser.add_argument('--min', '-N', help='The minimal length of target amplicons (from F2 to R2)', default=int(100))
	parser.add_argument('--ConsiderOutgroups', '-CO', help='Whether outgroups are considered when designing universal primers. 0: Do not consider; 1: Consider.', required=True)
	parser.add_argument('--PIs','-PIs', help=' The weighting ratio of PCR score and Information score when scoring all possible nested-PCR primers within the same MSA. Increasing this value will reduce the length of the amplified regions. Please set this parameter with caution',default=int(1))

	args = parser.parse_args()
	print args
	try:
		print Pipeline(args.dir, args.con, args.fwd, args.rvs, args.all, args.inthred2, args.allthred2, args.inthred1, args.allthred1, args.max, args.min, args.ConsiderOutgroups, args.PIs)
	except Exception as e:
		print(e)
