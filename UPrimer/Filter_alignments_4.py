#!/usr/bin/python
# -*- coding: UTF-8 -*-

import glob, os, subprocess
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import AlignIO
import argparse
linesep = os.linesep

def Get_gap_and_seq_position(start, record):
	i = start
	if record[start] != '-':
		for n in range(start, len(record)):
			if record[n] == '-':
				break
			else:
				i = i + 1
	else:
		for n in range(start, len(record)):
			if record[n] != '-':
				break
			else:
				i = i + 1
	return i

def Get_gap_and_seq_distribution_from_one_seq(seq):
	range_1 = []
	range_2 = []
	all = []
	n = 0
	all.append(n)
	
	while (n < len(seq)):
		n = Get_gap_and_seq_position(n, seq)
		all.append(n)
	if len(all) % 2 == 0:
		for i in range(0, len(all), 2):
			first = all[i]
			second = all[i + 1]
			range_1.append([first, second])
		for j in range(1, len(all) - 1, 2):
			first = all[j]
			second = all[j + 1]
			range_2.append([first, second])
	else:
		for i in range(0, len(all) - 1, 2):
			first = all[i]
			second = all[i + 1]
			range_1.append([first, second])
		for j in range(1, len(all) - 1, 2):
			first = all[j]
			second = all[j + 1]
			range_2.append([first, second])

	return range_1, range_2, all

def Get_trim_forward_position_1(range_A, range_G, trim_length, trim_times):
	length_A_f = []
	length_G_f = []
	length_A = 0
	length_G = 0
	start = range_A[0][0]
	point = 0
	for n in range(0, len(range_A) - 1):
		length = range_A[n][1] - range_A[n][0]
		if length <= trim_length:
			length_A_f.append(range_A[n])
			length_G_f.append(range_G[n])
		else:
			break
		
	if len(length_A_f) > 0:
		for k in range(0, len(length_A_f)):
			length_A = length_A + length_A_f[k][1] - length_A_f[k][0] # 0 + 
			length_G = length_G + length_G_f[k][1] - length_G_f[k][0]
			if float(length_G) / length_A >= trim_times:
				start = length_G_f[k][1]
				point = k + 1
			else:
				pass
		if point != 0:
			return start, length_A_f[:point]
		else:
			return start, []
	else:
		return start, []
	
def Get_trim_forward_position_2(range_A, range_G,trim_length, trim_times):
	length_A_f = []
	length_G_f = []
	length_A = 0
	length_G = 0
	start = range_A[0][0]
	point = 0
	length_G_f.append(range_G[0])
	for n in range(0, len(range_A) - 1):
		length = range_A[n][1] - range_A[n][0]
		if length <= trim_length:
			length_A_f.append(range_A[n])
			length_G_f.append(range_G[n + 1])
		else:
			break
	if len(length_A_f) > 0:
		length_G = length_G + length_G_f[0][1] - length_G_f[0][0]
		for k in range(0, len(length_A_f)):
			length_A = length_A + length_A_f[k][1] - length_A_f[k][0]
			length_G = length_G + length_G_f[k + 1][1] - length_G_f[k + 1][0]
			if float(length_G) / length_A >= trim_times:
				start = length_G_f[k + 1][1]
				point = k + 1
			else:
				pass
		if point != 0:
			return start, length_A_f[:point]
		else:
			return start, []
	else:
		return start, []
	
def Get_trim_reverse_position_1(range_A, range_G, trim_length, trim_times):
	length_A_r = []
	length_G_r = []
	length_A = 0
	length_G = 0
	end = range_A[-1][1]
	point = 0
	for n in range(-1, -len(range_A)-1, -1):
		length = range_A[n][1] - range_A[n][0]
		if length <= trim_length:
			if abs(n) > len(range_G):
				pass
			else:
				length_A_r.append(range_A[n])
				length_G_r.append(range_G[n])
		else:
			break
	if len(length_A_r) > 0:
		for k in range(0, len(length_A_r)):
			length_A = length_A + length_A_r[k][1] - length_A_r[k][0]
			length_G = length_G + length_G_r[k][1] - length_G_r[k][0]
			if float(length_G) / length_A >= trim_times:
				end = length_G_r[k][0]
				point = -k - 1
			else:
				pass
		if point != 0:
			return end, length_A_r[point:]
		else:
			return end, []
	else:
		return end, []
	
def Get_trim_reverse_position_2(range_A, range_G, trim_length, trim_times):
	length_A_r = []
	length_G_r = []
	end = range_G[-1][0]
	point = 0
	for n in range(-1, -len(range_A) - 1, -1):
		length = range_A[n][1] - range_A[n][0]
		if length <= trim_length:
			length_A_r.append(range_A[n])
			length_G_r.append(range_G[n])
		else:
			length_G_r.append(range_G[n])
			break
	length_A = 0
	length_G = 0
	
	if len(length_A_r) > 0:
		for k in range(0, len(length_A_r)):
			length_A = length_A + length_A_r[k][1] - length_A_r[k][0]
			length_G = length_G + length_G_r[k][1] - length_G_r[k][0]
			if float(length_G) / length_A >= trim_times:
                                end = length_G_r[k][0]
				point = k + 1
			else:
				pass
		if point != 0:
			return end, length_A_r[:point]
		else:
			return end, []
	else:
		return end, []

def Calculate_taxa_number(seq_list):
    outgroup = 0
    ingroup = 0
	    
    for seq in seq_list:
	if 'Outgroup' in seq:
	    outgroup += 1
	elif 'Ingroup' in seq:
	    ingroup += 1
    return outgroup, ingroup


def Trim_pipetide_MSAs(files, trim_dir, trim_length, trim_times,outgroup_num,ingroup_num,totaltaxanum,ConsiderOutgroups):
	
	path_out = os.path.join(os.getcwd(), trim_dir)
	if not os.path.exists(path_out):
		os.mkdir(path_out)
		
	path_file = os.path.join(path_out, os.path.split(files)[1])
	if not os.path.exists(path_file):
		records = path_file + linesep
		outlist = []


		align_pep = AlignIO.read(files, 'fasta')
		MSA_length = len(align_pep[0].seq)
                
		for record in SeqIO.parse(files, 'fasta'):			
			seq = str(record.seq)			
			range_all = Get_gap_and_seq_distribution_from_one_seq(seq)
			
			delete_list = []
			char = 'X'
			
			ID = str(record.id)		
			if '_exon_' in ID:
				record.id = 'Refer_' + ID.split('_exon_')[0]
				record.description = 'Refer_' + ID.split('_exon_')[0]
			elif '_cds_' in ID:
				record.id = 'Outgroup_' + ID.split('_cds_')[0]
				record.description = 'Outgroup_' + ID.split('_cds_')[0]
			elif '_genome_' in ID:
				record.id = 'Ingroup_' + ID.split('_genome_')[0]
				record.description = 'Ingroup_' + ID.split('_genome_')[0]

				
			if len(range_all[2]) < 4:
				outlist.append('>' + str(record.description) + '\n' + str(record.seq) + '\n')
				
			elif seq[0] != '-': 
					range_AA = range_all[0]
					range_GAP = range_all[1]

					forward = Get_trim_forward_position_1(range_AA, range_GAP, trim_length, trim_times)

					start = forward[0]
					delete_list = delete_list + forward[1]
					
					if seq[-1] != '-':
						reverse = Get_trim_reverse_position_1(range_AA, range_GAP, trim_length, trim_times)
						end = reverse[0]
						if end > start:
							delete_list = delete_list + reverse[1]
							reverse_length = range_AA[-1][1] - end
							if float(start + reverse_length)/MSA_length <0.4:
                                                                record.seq = Seq(str(char * start + seq[start:end] + char * reverse_length))
                                                                outlist.append('>' + str(record.description) + '\n' + str(record.seq) + '\n')

					else:
                                                
						reverse = Get_trim_reverse_position_2(range_AA, range_GAP, trim_length, trim_times)
						
						end = reverse[0]
						if end > start:
							delete_list = delete_list + reverse[1]
							reverse_length = range_GAP[-1][1] - end
							if float(start + reverse_length)/MSA_length <0.4:
                                                                record.seq = Seq(str(char * start + seq[start:end] + char * reverse_length))
                                                                outlist.append('>' + str(record.description) + '\n' + str(record.seq) + '\n')
					
					
			else:                                
				range_AA = range_all[1]				
				range_GAP = range_all[0]
				forward = Get_trim_forward_position_2(range_AA, range_GAP, trim_length, trim_times)
				start = forward[0]				
				delete_list = delete_list + forward[1]

				if seq[-1] != '-':
					reverse = Get_trim_reverse_position_1(range_AA, range_GAP, trim_length, trim_times)
					end = reverse[0]
					if end > start:
						delete_list = delete_list + reverse[1]
						reverse_length = range_AA[-1][1] - end
						if float(start + reverse_length)/MSA_length <0.4:
                                                        record.seq = Seq(str(char * start + seq[start:end] + char * reverse_length))
                                                        outlist.append('>' + str(record.description) + '\n' + str(record.seq) + '\n')
				else:
					reverse = Get_trim_reverse_position_2(range_AA, range_GAP, trim_length, trim_times)
					end = reverse[0]
					delete_list = delete_list + reverse[1]
					reverse_length = range_GAP[-1][1] - end
					if float(start + reverse_length)/MSA_length <0.4:
                                                record.seq = Seq(str(char * start + seq[start:end] + char * reverse_length))
                                                outlist.append('>' + str(record.description) + '\n' + str(record.seq) + '\n')

			
			ID = record.id
			if len(delete_list) > 0:
				records = records + ID + ':' + str(delete_list) + linesep
			else:
				pass

                if len(outlist) != 0:
                        if int(ConsiderOutgroups) == 1:
                                
                                if len(outlist) >= int(totaltaxanum) and Calculate_taxa_number(outlist)[0] >= outgroup_num and Calculate_taxa_number(outlist)[1] >= ingroup_num:

                                        outfile = open(path_file,'w')
                                        for seq in outlist:
                                                outfile.write(seq)
                                        outfile.close()
                                        print path_file, ' trim done!'
                        else:
                                if len(outlist) >= int(totaltaxanum) and Calculate_taxa_number(outlist)[1] >= int(ingroup_num):

                                        outfile = open(path_file,'w')
                                        for seq in outlist:
                                                outfile.write(seq)
                                        outfile.close()
                                        print path_file, ' trim done!' 

			
		return path_file, records

	else:
		print path_file, 'have been trimed !'
		return path_file, None

def Trim_nucleotide_MSAs(files, trim_dir, trim_length, trim_times, outgroup_num, ingroup_num, totaltaxanum, ConsiderOutgroups):
	
	path_out = os.path.join(os.getcwd(), trim_dir)
	if not os.path.exists(path_out):
		os.mkdir(path_out)
		
	path_file = os.path.join(path_out, os.path.split(files)[1])
	if not os.path.exists(path_file):
		records = path_file + linesep
		outlist = []


		align_pep = AlignIO.read(files, 'fasta')
		MSA_length = len(align_pep[0].seq)
                
		for record in SeqIO.parse(files, 'fasta'):			
			seq = str(record.seq)
			range_all = Get_gap_and_seq_distribution_from_one_seq(seq)
			
			delete_list = []
			char = 'N'
			
			ID = str(record.id)		
			if '_exon_' in ID:
				record.id = 'Refer_' + ID.split('_exon_')[0]
				record.description = 'Refer_' + ID.split('_exon_')[0]
			elif '_cds_' in ID:
				record.id = 'Outgroup_' + ID.split('_cds_')[0]
				record.description = 'Outgroup_' + ID.split('_cds_')[0]
			elif '_genome_' in ID:
				record.id = 'Ingroup_' + ID.split('_genome_')[0]
				record.description = 'Ingroup_' + ID.split('_genome_')[0]

				
			if len(range_all[2]) < 4:
				outlist.append('>' + str(record.description) + '\n' + str(record.seq) + '\n')
				
			elif seq[0] != '-': 
					range_AA = range_all[0]
					range_GAP = range_all[1]

					forward = Get_trim_forward_position_1(range_AA, range_GAP, trim_length, trim_times)

					start = forward[0]
					delete_list = delete_list + forward[1]
					
					if seq[-1] != '-':
						reverse = Get_trim_reverse_position_1(range_AA, range_GAP, trim_length, trim_times)
						end = reverse[0]
						if end > start:
							delete_list = delete_list + reverse[1]
							reverse_length = range_AA[-1][1] - end
							if float(start + reverse_length)/MSA_length <0.4:
                                                                record.seq = Seq(str(char * start + seq[start:end] + char * reverse_length))
                                                                outlist.append('>' + str(record.description) + '\n' + str(record.seq) + '\n')

					else:
                                                
						reverse = Get_trim_reverse_position_2(range_AA, range_GAP, trim_length, trim_times)
						
						end = reverse[0]
						if end > start:
							delete_list = delete_list + reverse[1]
							reverse_length = range_GAP[-1][1] - end
							if float(start + reverse_length)/MSA_length <0.4:
                                                                record.seq = Seq(str(char * start + seq[start:end] + char * reverse_length))
                                                                outlist.append('>' + str(record.description) + '\n' + str(record.seq) + '\n')
					
					
			else:                                
				range_AA = range_all[1]				
				range_GAP = range_all[0]
				forward = Get_trim_forward_position_2(range_AA, range_GAP, trim_length, trim_times)
				start = forward[0]				
				delete_list = delete_list + forward[1]

				if seq[-1] != '-':
					reverse = Get_trim_reverse_position_1(range_AA, range_GAP, trim_length, trim_times)
					end = reverse[0]
					if end > start:
##                                                print 'end > start'
						delete_list = delete_list + reverse[1]
						reverse_length = range_AA[-1][1] - end
						if float(start + reverse_length)/MSA_length <0.4:
                                                        record.seq = Seq(str(char * start + seq[start:end] + char * reverse_length))
                                                        outlist.append('>' + str(record.description) + '\n' + str(record.seq) + '\n')
				else:
					reverse = Get_trim_reverse_position_2(range_AA, range_GAP, trim_length, trim_times)
					end = reverse[0]
					delete_list = delete_list + reverse[1]
					reverse_length = range_GAP[-1][1] - end
					if float(start + reverse_length)/MSA_length <0.4:
                                                record.seq = Seq(str(char * start + seq[start:end] + char * reverse_length))
                                                outlist.append('>' + str(record.description) + '\n' + str(record.seq) + '\n')

			
			ID = record.id
			if len(delete_list) > 0:
				records = records + ID + ':' + str(delete_list) + linesep
			else:
				pass

                if len(outlist) != 0:
                        if int(ConsiderOutgroups) == 1:
                                
                                if len(outlist) >= int(totaltaxanum) and Calculate_taxa_number(outlist)[0] >= outgroup_num and Calculate_taxa_number(outlist)[1] >= ingroup_num:

                                        outfile = open(path_file,'w')
                                        for seq in outlist:
                                                outfile.write(seq)
                                        outfile.close()
                                        print path_file, ' trim done!'
                        else:
                                if len(outlist) >= int(totaltaxanum) and Calculate_taxa_number(outlist)[1] >= int(ingroup_num):

                                        outfile = open(path_file,'w')
                                        for seq in outlist:
                                                outfile.write(seq)
                                        outfile.close()
                                        print path_file, ' trim done!' 

			
		return path_file, records

	else:
		print path_file, 'have been trimed!'
		return path_file, None

		
if __name__ == '__main__':
	dirpath_o = os.path.abspath(os.path.join(os.path.dirname('__file__'), os.path.pardir))
	os.chdir(dirpath_o)
	parser = argparse.ArgumentParser(description='''python Filter_alignments_4.py -O 3.5.Trim_peptide_alignments -CO [1: Consider outgroups;0: Do not consider outgroups]''')

	
	parser.add_argument('--Out', '-O', help='Output directory (Folder name: 3.5.Trim_peptide_alignments)', required=True)
	parser.add_argument('--trml', '-l', help='The permissible length cutoff of AA during trimming sequences', default=int(10))
	parser.add_argument('--trmt', '-t', help='The ratio of gap length and AA length during trimming sequences', default=float(1.5))	
	parser.add_argument('--outgroupNum', '-ON', help='The minimal number of outgroups in trimed MSAs', default=int(1))
	parser.add_argument('--ingroupNum', '-IN', help='The minimal number of ingroups in trimed MSAs', default=int(1))
	parser.add_argument('--totaltaxanum', '-TN', help='The minimal number of total species in trimed MSAs', default=int(3))
	parser.add_argument('--ConsiderOutgroups', '-CO', help='Whether outgroups are considered when designing universal primers. 0: Do not consider; 1: Consider.', required=True)
	args = parser.parse_args()
	print args

	Input_dir = os.path.join(os.getcwd(), '3.1.Megacc_aligned_peptide_MSAs', 'Tra_fasta_peptide','Cut_Alignment_2')
	Files = glob.glob(os.path.join(os.getcwd(), Input_dir, '*.fasta'))
	for File in Files:
		Trim_pipetide_MSAs(File, args.Out, args.trml, args.trmt, args.outgroupNum, args.ingroupNum, args.totaltaxanum, args.ConsiderOutgroups)

