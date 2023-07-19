#!/usr/bin/python
# -*- coding: UTF-8 -*-
from Bio.Blast import NCBIXML
from Bio import SeqIO, SearchIO
import os, subprocess, sys, getopt, glob
import pandas as pd
import numpy as np
from Make_codon_alignments_3 import Megacc
import Filter_alignments_4 as Tb
import Make_codon_alignments_3 as Mega
from Bio.Alphabet import IUPAC
import argparse
linesep = os.linesep


def Analyze_forward_blast_results(infile):
	RATIO = 1.0
	dic = {}
	with open(infile) as in_file:
		for records in NCBIXML.parse(in_file):
			query = str(records.query)
			first_score = 1
			second_score = 1
			pointer_1 = 0
			pointer_2 = 0
			f_hit = None
			
			for alignment in records.alignments:
				pointer_1 += 1
				for HSP in alignment.hsps:
					pointer_2 += 1
					if pointer_1 == 1:
						if pointer_2 == 1:
							first_score = HSP.score
							f_hit = str(alignment.hit_def)
							break
					elif pointer_1 == 2:
						if pointer_2 == 1:
							second_score = HSP.score
							break
			if first_score / second_score >= RATIO:
				dic[query] = f_hit
	return dic

def Analyze_reverse_blast_results(infile):
	RATIO = 1.0
	dic = {}
	r_local = {}
	with open(infile) as in_file:
		for records in NCBIXML.parse(in_file):
			query = str(records.query)
			first_score = 1
			second_score = 1
			pointer_1 = 0
			pointer_2 = 0
			r_hit = None
			start = 1
			end = 1
			frame = 1
			for alignment in records.alignments:
				pointer_1 += 1
				for HSP in alignment.hsps:
					pointer_2 += 1
					if pointer_1 == 1:
						if pointer_2 == 1:
							start = int(HSP.query_start)
							end = int(HSP.query_end)
							first_score = HSP.score
							r_hit = str(alignment.hit_def)
							frame = int(HSP.frame[0])
							break
					elif pointer_1 == 2:
						if pointer_2 == 1:
							second_score = HSP.score
							break
			if first_score / second_score >= RATIO:
				dic[query] = r_hit
				START = str(start - 1)
				END = str(end)
				if frame < 0:
					r_local.update({query: (START + ',' + END + '_reverse_complement')})
				else:
					r_local.update({query: (START + ',' + END)})
	return dic, r_local
	
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

def Summary_MHB_results(species_list):

	sp_list = os.path.join(os.getcwd(),species_list)
	reference_list = Get_species_list(sp_list)[0]
	taxa_list = Get_species_list(sp_list)[3]
	
	XML_filename_lists = []
	total_r_local = []
	total_r_list = []
	total_f_list = []

	refer_name_pep = str(reference_list[0]) + '_pep'
	refer_name_pre = str(reference_list[0])

	taxa_names = []
	for i in range(0, len(taxa_list)):
                taxa_name = str(taxa_list[i])
                taxa_names.append(taxa_name)
                xml_1 = str(refer_name_pep) + '_tblastn_' + str(taxa_name) + '_hitseq.xml'
                xml_2 = str(taxa_name) + '_hitseq_blastx_' + str(refer_name_pep) + '.xml'
                XML_filename_lists.append(xml_1)
                XML_filename_lists.append(xml_2)                        

   # print len(XML_filename_lists)
	for i in range(0, len(XML_filename_lists) / 2):
		f_name = os.path.join(os.getcwd(), '1.HitSeq_xml', XML_filename_lists[i * 2])
		r_name = os.path.join(os.getcwd(), '1.HitSeq_xml', XML_filename_lists[i * 2 + 1])

		print os.path.basename(r_name)
		r_dic = Analyze_reverse_blast_results(r_name)
		r_dic_1 = r_dic[0]
		r_dic_2 = r_dic[1]
		total_r_list.append(r_dic_1)
		total_r_local.append(r_dic_2)

		print os.path.basename(f_name)
		f_dic = Analyze_forward_blast_results(f_name)
		total_f_list.append(f_dic)
		print i+1, '  xml file-pair analyzed!'
	return total_f_list, total_r_list, total_r_local, taxa_names


def Create_sequence_dictionary(species_list,cut_length):
	dic_list = []
        sp_list = os.path.join(os.getcwd(),species_list)
	reference_list = Get_species_list(sp_list)[0]
	taxa_list = Get_species_list(sp_list)[3]

	refer_sp = os.path.join(os.getcwd(), 'Reference', 'single_copy_Long_' + str(cut_length) + '_' + str(reference_list[0]) +  '_cds.fasta')
	dic = SeqIO.to_dict(SeqIO.parse(refer_sp, "fasta"), key_function=lambda r: r.description)

	## XP_038957404.1_GeneID120100082_LOC120100082_NC_051336.1_1335317_1335628_-1_retinoic acid early-inducible protein 1-gamma-like
	print refer_sp, ' dictionary done!'
	dic_list.append(dic)

	ref_pre = str(reference_list[0])

	for taxa in taxa_list:
                filename = str(taxa) + '_hitseq.fasta'
                print filename, ' begin creating dictionary !'

                species = os.path.join(os.getcwd(), '1.HitSeq', filename)
                dic = SeqIO.to_dict(SeqIO.parse(species, "fasta"), key_function=lambda r: r.description) ## r.description = XP_047753575.1|1232131|312324
                dic_list.append(dic)
        return dic_list
                


def Extract_orthologous_sequence_groups(line, dic_list, title, min_number):

	m = str(line.strip()).count('None')
	n = len(title)
	number = n - m
	if number >= min_number:
		genes = line.strip().split('\t')

                
		#Format: >XP_047753575.1_GeneID71966332_GeneID71966332_NC_062999.1_100866_101500_1_ISWI chromatin-remodeling complex ATPase ISW2
		if '_unknown' not in genes[0]:
                        pre_msa_name = genes[0]
                        mas_name = pre_msa_name.split('_')[0] + '_' + pre_msa_name.split('_')[1] + '_' + pre_msa_name.split('_')[6].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','') + '_' +\
                                   pre_msa_name.split('_')[7].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','')

                        path_out_nucl = os.path.join(os.getcwd(), '2.Nucleotide_orthologous_groups', mas_name + '_nucl.fasta')
                        path_out_pep = os.path.join(os.getcwd(), '2.Peptide_orthologous_groups', mas_name +  '_pep.fasta')
                                   
                #Format: >EUB53855.1_unknown_unknown_APAU02000760.1_556_895_-1_hypothetical protein
                elif '_unknown' in genes[0]:
                        pre_msa_name = genes[0].replace('_unknown','') #>EUB53855.1_APAU02000760.1_556_895_-1_hypothetical protein
                        mas_name = pre_msa_name.split('_')[0] + '_' + pre_msa_name.split('_')[2].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','') +\
                                   '_' + pre_msa_name.split('_')[3].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','')

                        path_out_nucl = os.path.join(os.getcwd(), '2.Nucleotide_orthologous_groups', mas_name + '_nucl.fasta')
                        path_out_pep = os.path.join(os.getcwd(), '2.Peptide_orthologous_groups', mas_name +  '_pep.fasta')
                     
                #Format: RPRC010500|RPRC010500-RA|RPRC010500-RA-E3|T1I2I1
                else:
##                        print genes[0]
                        path_out_nucl = os.path.join(os.getcwd(), '2.Nucleotide_orthologous_groups', str(genes[0]).split('|')[-2] + '_nucl.fasta')
                        path_out_pep = os.path.join(os.getcwd(), '2.Peptide_orthologous_groups', str(genes[0]).split('|')[-2] +  '_pep.fasta')
                        

		check = 'T'
		if not os.path.exists(path_out_nucl) and not os.path.exists(path_out_pep):
			with open(path_out_nucl, 'w') as fh_out_nucl, open(path_out_pep, 'w') as fh_out_pep:
				for i in range(0, len(genes)):
					if str(genes[i]) != 'None':						
						### for referecne exon                        
						if i == 0:
                                                        if 'GeneID ' in genes[0] or 'XP_' in genes[0] or 'NC_' in genes[0]:
                                                                ## Format: >XP_047753575.1_GeneID71966332_GeneID71966332_NC_062999.1_100866_101500_1_ISWI chromatin-remodeling complex ATPase ISW2                                                                                                      
                                                                stop_seq = str(dic_list[i][str(genes[i])].seq[-3:]).upper()
                                                                all_length_seq = dic_list[i][str(genes[i])].seq[:-3]
                                           
                                                                if '*' in str(all_length_seq.translate()):
                                                                        check = 'F'
                                                                elif 'TAG' == stop_seq or 'TGA' == stop_seq or 'TAA' == stop_seq: #new_format: >XP_047753575.1_100866_101500                                                                                                                               
                                                                        
                                                                        record_nucl = '>' + title[i].strip() + '_' + str(genes[0]).split('_')[0] + '_' +\
                                                                                                  str(genes[0]).split('_')[1] + '_' +\
                                                                                                  str(genes[0]).split('_')[6].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','') + '_' +\
                                                                                                  str(genes[0]).split('_')[7].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','') + linesep + str(all_length_seq)
                                                                        fh_out_nucl.write(record_nucl + linesep)
                                                                        record_pep = '>' + title[i].strip() + '_' + str(genes[0]).split('_')[0] + '_' +\
                                                                                                  str(genes[0]).split('_')[1] + '_' +\
                                                                                                  str(genes[0]).split('_')[6].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','') + '_' +\
                                                                                                  str(genes[0]).split('_')[7].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','') +  linesep + str(all_length_seq.translate())
                                                                        fh_out_pep.write(record_pep + linesep)
                                                                else:
                                                                        all_length_seq = dic_list[i][str(genes[i])].seq
                                                                        record_nucl = '>' + title[i].strip() + '_' + str(genes[0]).split('_')[0] + '_' +\
                                                                                                  str(genes[0]).split('_')[1] + '_' +\
                                                                                                  str(genes[0]).split('_')[6].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','') + '_' +\
                                                                                                  str(genes[0]).split('_')[7].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','') + linesep + str(all_length_seq)
                                                                        fh_out_nucl.write(record_nucl + linesep)
                                                                        record_pep = '>' + title[i].strip() + '_' + str(genes[0]).split('_')[0] + '_' +\
                                                                                                  str(genes[0]).split('_')[1] + '_' +\
                                                                                                  str(genes[0]).split('_')[6].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','') + '_' +\
                                                                                                  str(genes[0]).split('_')[7].replace('<','').replace('>','').replace('/','_').replace('(','').replace(')','') + linesep + str(all_length_seq.translate())
                                                                        fh_out_pep.write(record_pep + linesep)

                                                        else:#Format: RPRC010500|RPRC010500-RA|RPRC010500-RA-E3|T1I2I1
                                                                stop_seq = str(dic_list[i][str(genes[i])].seq[-3:]).upper()
                                                                all_length_seq = dic_list[i][str(genes[i])].seq[:-3]
                                                                
                                                                if '*' in str(all_length_seq.translate()):
                                                                        check = 'F'
                                                                elif 'TAG' == stop_seq or 'TGA' == stop_seq or 'TAA' == stop_seq:
                                                                        record_nucl = '>' + title[i].strip() + '_' + str(genes[0]).split('|')[-2] + '|' + str(genes[0]).split('|')[-1] + linesep + str(all_length_seq)
                                                                        fh_out_nucl.write(record_nucl + linesep)

                                                                        record_pep = '>' + title[i].strip() + '_' + str(genes[0]).split('|')[-2] + '|' + str(genes[0]).split('|')[-1] + linesep + str(all_length_seq.translate())
                                                                        fh_out_pep.write(record_pep + linesep)
                                                                else:
                                                                        all_length_seq = dic_list[i][str(genes[i])].seq
                                                                        record_nucl = '>' + title[i].strip() + '_' + str(genes[0]).split('|')[-2] + '|' + str(genes[0]).split('|')[-1] + linesep + str(all_length_seq)
                                                                        fh_out_nucl.write(record_nucl + linesep)

                                                                        record_pep = '>' + title[i].strip() + '_' + str(genes[0]).split('|')[-2] + '|' + str(genes[0]).split('|')[-1] + linesep + str(all_length_seq.translate())
                                                                        fh_out_pep.write(record_pep + linesep)                                                  
                                                                        

						### for ingroup and outgroup species
                                                else:
                                                        
                                                        ## Format: XP_047753575.1|1232131|312324 8814,9171_reverse_complement
							if 'reverse_complement' not in genes[i]:
								name = str(genes[i].split(' ')[0])
								start = int(genes[i].split(' ')[1].split(',')[0])
                                                                                                        
								end = int(genes[i].split(' ')[1].split(',')[1].strip())
								stop_seq_genome = str(dic_list[i][name].seq[end-3:end]).upper()


								
								
								if 'TAG' == stop_seq_genome or 'TGA' == stop_seq_genome or 'TAA' == stop_seq_genome:
									nucl_seq = str(dic_list[i][name].seq[start:end - 3])
									pep_seq = str(dic_list[i][name].seq[start:end - 3].translate())
									if '*' not in pep_seq and len(pep_seq) >= 100:
										record_nucl = '>' + title[i].strip() + '_' + str(name) + linesep + nucl_seq
										record_pep = '>' + title[i].strip() + '_' + str(name) + linesep + pep_seq
										fh_out_nucl.write(record_nucl + linesep)
										fh_out_pep.write(record_pep + linesep)
									else:
										pass
								else:
									nucl_seq = str(dic_list[i][name].seq[start:end])
									pep_seq = str(dic_list[i][name].seq[start:end].translate())
									if '*' not in pep_seq and len(pep_seq) >= 100:
										record_nucl = '>' + title[i].strip() + '_' + str(name) + linesep + nucl_seq
										record_pep = '>' + title[i].strip() + '_' + str(name) + linesep + pep_seq
										fh_out_nucl.write(record_nucl + linesep)
										fh_out_pep.write(record_pep + linesep)
									else:
										pass
							else:
								name = str(genes[i].split(' ')[0])
								start = int(genes[i].split(' ')[1].split(',')[0])
								end = int(genes[i].split(' ')[1].split(',')[1].replace('_reverse_complement', '').strip())
								stop_seq_genome = str(dic_list[i][name].seq[:3]).upper()



								
								if 'TCA' == stop_seq_genome or 'TTA' == stop_seq_genome or 'CTA' == stop_seq_genome:
									nucl_seq = str(dic_list[i][name].seq[start + 3:end].reverse_complement())
									pep_seq = str(dic_list[i][name].seq[start + 3:end].reverse_complement().translate())
									if '*' not in pep_seq and len(pep_seq) >= 100:
										record_nucl = '>' + title[i].strip() + '_' + str(name) + linesep + nucl_seq
										record_pep = '>' + title[i].strip() + '_' + str(name) + linesep + pep_seq
										fh_out_nucl.write(record_nucl + linesep)
										fh_out_pep.write(record_pep + linesep)
									else:
										pass
								else:
									nucl_seq = str(dic_list[i][name].seq[start:end].reverse_complement())
									pep_seq = str(dic_list[i][name].seq[start:end].reverse_complement().translate())
									if '*' not in pep_seq and len(pep_seq) >= 100:
										record_nucl = '>' + title[i].strip() + '_' + str(name) + linesep + nucl_seq
										record_pep = '>' + title[i].strip() + '_' + str(name) + linesep + pep_seq
										fh_out_nucl.write(record_nucl + linesep)
										fh_out_pep.write(record_pep + linesep)
									else:
										pass

					
			records = list(SeqIO.parse(path_out_pep, 'fasta'))
			if len(records) >= min_number:
				if check == 'T':
					return str(path_out_pep)
				else:
					return 'None'
			else:
				return 'None'
		else: 
			print 'Extracted unaligned '+ str(genes[0]) + '.fasta file have existed !'
			records = list(SeqIO.parse(path_out_pep, 'fasta'))

			if len(records) >= min_number:
				return str(path_out_pep)
			else:
				return 'None'
	else:
		return 'None'

	
def Generate_MHB_txt_result(single_cds, species_list):

	sp_list = os.path.join(os.getcwd(), species_list)
	reference_list = Get_species_list(sp_list)[0]
	taxa_list = Get_species_list(sp_list)[3]
	
	
	dic_list = Summary_MHB_results(sp_list)
	dic_f_list = dic_list[0]
	dic_r_list = dic_list[1]
	dic_r_local = dic_list[2]
	species_names = dic_list[3]

	
	with open("1.Homology_orthology_species(1.0-E10).txt", "a") as out_file:
                temp = str(reference_list[0]) + '_exon'
                for name in taxa_list:
                        temp = temp + '\t' + name
                temp = temp + linesep
                out_file.write(temp)
			
		for seq_record in SeqIO.parse(single_cds, "fasta"):
			
			temp1 = str(seq_record.description)
			seq_name = str(seq_record.description)

			for i in range(0, len(species_names)):
				kknd = 'None'
				if temp1 in dic_f_list[i].keys():
					query = dic_f_list[i][str(seq_record.description)]

					if query in dic_r_list[i].keys():
						hit = dic_r_list[i][query]

						if hit == str(seq_record.description):
							if dic_r_local[i] != 'None':
								kknd = query + ' ' + dic_r_local[i][query]
							else:
								kknd = query
				else:
					pass
				seq_name = seq_name + '\t' + kknd
			seq_name = seq_name + linesep
			out_file.write(seq_name)    



def Get_OGs_from_the_MHB_result_and_make_alignments(txt,species_list, cut_length, min, select_scale, remove_rigor):
	path_nucl = os.path.join(os.getcwd(), '2.Nucleotide_orthologous_groups')
	path_pep = os.path.join(os.getcwd(), '2.Peptide_orthologous_groups')
	
	if not os.path.exists(path_nucl):
		os.mkdir(path_nucl)
	if not os.path.exists(path_pep):
		os.mkdir(path_pep)
	with open(txt, 'r') as fh_in:
		
		dic_list = Create_sequence_dictionary(species_list, cut_length)
		print 'Dictionary have created!'
		
		title = fh_in.readline().strip().split('\t')
		for line in fh_in.readlines():
			align_prot = Extract_orthologous_sequence_groups(line, dic_list, title, min)
			
			if align_prot != 'None':
				fa_prot = Megacc(align_prot, 'muscle_align_protein.mao')
				meg_prot = fa_prot.mega_align_prot()
				meg_file_prot = meg_prot + '.meg'
				after_align_prot = fa_prot.tra_meg_fasta(meg_file_prot, error_seq='error_seq_1.txt')
				cut_prot = Mega.Cut_two_sides_of_MSAs(after_align_prot, '1_Taxa_information.txt', 'Cut_Alignment')
				Mega.Remove_rogue_taxon(cut_prot, species_list, '1_Rogue_taxa_remove.txt',select_scale, remove_rigor)
			else:
				pass

def Calculate_taxa_number(file):
    outgroup = 0
    ingroup = 0
    for record in SeqIO.parse(open(file), 'fasta'):
        if '_cds_' in str(record.id):
            outgroup += 1
        elif '_genome_' in str(record.id):
            ingroup += 1
    return outgroup, ingroup


def Pipeline(single_cds, species_list, cut_length, select_scale, remove_rigor, trim_length, trim_times, outgroupNum, ingroupNum, totaltaxanum, ConsiderOutgroups):
	if not os.path.exists("1.Homology_orthology_species(1.0-E10).txt"):

		Generate_MHB_txt_result(single_cds, species_list)		
		Get_OGs_from_the_MHB_result_and_make_alignments("1.Homology_orthology_species(1.0-E10).txt", species_list, cut_length, 3, select_scale, remove_rigor)
	else:
		print 'The MHB analysis result is exist!'
		Get_OGs_from_the_MHB_result_and_make_alignments("1.Homology_orthology_species(1.0-E10).txt", species_list, cut_length, 3, select_scale, remove_rigor)


	re_files_pep = Mega.Realign_MSAs_after_removing_rogue_taxon()    
	for re_file in re_files_pep:
		num_list = Calculate_taxa_number(re_file)
 
	        if int(ConsiderOutgroups) == 1:                        
                        if num_list[0] >= 1 and num_list[1] >= 1:#at least one ingroup and one outgroup
                                re_fa_prot = Megacc(re_file, 'muscle_align_protein.mao')
                                re_meg_prot = re_fa_prot.mega_align_prot()
                                re_meg_file_prot = re_meg_prot + '.meg'
                                re_after_align = re_fa_prot.tra_meg_fasta(re_meg_file_prot, error_seq='error_seq_2.txt')
                                re_after_cut = Mega.Cut_two_sides_of_MSAs(re_after_align, '2_Taxa_information_pep.txt', 'Cut_Alignment_2')	
                                after_trim = Tb.Trim_pipetide_MSAs(re_after_cut, '3.5.Refined_peptide_aligned_MSAs',trim_length, trim_times, outgroupNum, ingroupNum, totaltaxanum, ConsiderOutgroups)[0]			
                        else:
                                pass
                else:
                        if num_list[1] >= 1:#at least one ingroup and one outgroup
                                re_fa_prot = Megacc(re_file, 'muscle_align_protein.mao')
                                re_meg_prot = re_fa_prot.mega_align_prot()
                                re_meg_file_prot = re_meg_prot + '.meg'
                                re_after_align = re_fa_prot.tra_meg_fasta(re_meg_file_prot, error_seq='error_seq_2.txt')
                                re_after_cut = Mega.Cut_two_sides_of_MSAs(re_after_align, '2_Taxa_information_pep.txt', 'Cut_Alignment_2')		
                                after_trim = Tb.Trim_pipetide_MSAs(re_after_cut, '3.5.Refined_peptide_aligned_MSAs',trim_length, trim_times, outgroupNum, ingroupNum, totaltaxanum, ConsiderOutgroups)[0]			
                        else:
                                pass                        
                      
		
	Mega.Realign_nucleotide_MSAs(outgroupNum, ingroupNum, totaltaxanum, ConsiderOutgroups)

	

if __name__ == '__main__':
	dirpath_o = os.path.abspath(os.path.join(os.path.dirname('__file__'), os.path.pardir))
	os.chdir(dirpath_o)

	parser = argparse.ArgumentParser(description='''python Append_orthologous_sequences_2.py -F species_list.txt -D Coding_exon -CO [1: Consider outgroups;0: Do not consider outgroups]''')
                                         
	parser.add_argument('--len', '-L', help='The length cutoff for reference exons', default=int(300))
	parser.add_argument('--file', '-F', help='A text table containing information on the reference species, ingroup species, and outgroup species (filename: species_list.txt)', required=True)
	parser.add_argument("--SelectScale", "-SS", help="Value of a criteria used in removing rogue sequences, larger and more strict. This parameter usually does not need to be adjusted", default=float(1.1))
	parser.add_argument("--RemoveRigor", "-RR", help="Value of a sequence identity cutoff used in removing rogue sequences from MSAs. Set this parameter based on the quality of genome data. Larger value is suitbale for high quality genome data (50 or 55), lower value is suitbale for low quality genome data (30 or 35)", default=float(55))
	parser.add_argument('--outgroupNum', '-ON', help='The minimal number of outgroups', default=int(1))
	parser.add_argument('--ingroupNum', '-IN', help='The minimal number of ingroups', default=int(1))
	parser.add_argument('--totaltaxanum', '-TN', help='The minimal number of all species', default=int(3))
        parser.add_argument('--tril', '-l', help='The permissible length cutoff of amino acid during trimming sequences. This parameter usually does not need to be adjusted', default=int(10))
        parser.add_argument('--trit', '-t', help='The ratio of gap length and AA length during trimming sequences. This parameter usually does not need to be adjusted', default=float(1.5))
	parser.add_argument('--dir', '-D', help='Input a folder contains exon and peptide sequences of the reference species (Folder name: Reference)', required=True)
	parser.add_argument('--ConsiderOutgroups', '-CO', help='Whether considering outgroup sequences when designing primers. 0: Not considering. 1: Considering', required=True)
	args = parser.parse_args()
	print args

        dir_path = os.path.join(os.getcwd(), args.dir)                       
        cds = os.path.join(os.getcwd(), args.dir, '*_cds.fasta')
        
        try:
                print Pipeline(cds, args.file, args.len, args.SelectScale, args.RemoveRigor, args.tril, args.trit, \
                               args.outgroupNum, args.ingroupNum, args.totaltaxanum, args.ConsiderOutgroups)
        except Exception as e:
                print(e)
