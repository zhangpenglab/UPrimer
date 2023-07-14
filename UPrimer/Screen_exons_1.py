#!/usr/bin/python
# -*- coding: UTF-8 -*-

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbitblastnCommandline
import gzip, glob
import os, subprocess, sys, getopt
import argparse
import pandas as pd
import numpy as np
linesep = os.linesep

class Blast(object):
	def __init__(self, query, db):
		self.query = query
		self.db = db

	def makeblastdb(self, dbtype):
		print "Make BLAST database..."
		outDB = str(self.db).replace('.fasta', '')
		outDBabs = os.path.abspath(outDB)
		if not os.path.isdir(os.path.dirname(outDBabs)):
			os.mkdir(os.path.dirname(outDBabs))
		if dbtype == 'nucl':
			if not os.path.exists(outDBabs + '.nhr') or not os.path.exists(outDBabs + '.nin') or not os.path.exists(outDBabs + '.nsq'):
				command1 = "makeblastdb -in " + str(self.db) + " -out " + str(outDBabs) + " -dbtype " + str(dbtype)
				os.system(command1)
			else:
				print "A BLAST database has existed."
		elif dbtype == 'prot':
			if not os.path.exists(outDBabs + ".phr") or not os.path.exists(outDBabs + ".pin") or not os.path.exists(
					outDBabs + ".psq"):
				command1 = "makeblastdb -in " + str(self.db) + " -out " + str(outDBabs) + " -dbtype " + str(dbtype)
				os.system(command1)
			else:
				print "A BLAST database has existed."
		else:
			print "There is something wrong with this dataset."

	def blastn(self, outXML, num_thread):
		print 'Run BLASTN...'
		path = os.path.join(os.getcwd(), 'pre_xml')
		if not os.path.exists(path):
			os.mkdir(path)
		outXMLabs = os.path.join(path, outXML)
		if not os.path.isdir(os.path.dirname(outXMLabs)):
			os.makedirs(os.path.dirname(outXMLabs))
		if not os.path.exists(outXMLabs):
			commandn = 'blastn -task blastn -max_hsps 5 -evalue 1E-10 -outfmt 5 -query ' + str(self.query) + ' -num_threads ' + str(
				num_thread) + ' -out ' + str(outXMLabs) + ' -db ' + str(self.db).replace('.fasta', '')
			os.system(commandn)

			print 'BLASTN OK'
		else:
			print 'A BLAST XML has existed.'

	def blastx(self, outXML, num_thread):
		print 'Run BLASTX...'
		path = os.path.join(os.getcwd(), 'pre_xml')
		if not os.path.exists(path):
			os.mkdir(path)
		outXMLabs = os.path.join(path, outXML)
		if not os.path.isdir(os.path.dirname(outXMLabs)):
			os.makedirs(os.path.dirname(outXMLabs))
		if not os.path.exists(outXMLabs):
			commandx = 'blastx -evalue 1E-10 -outfmt 5 -query ' + str(self.query) + ' -max_target_seqs 5 -num_threads ' + str(
				num_thread) + ' -out ' + str(outXMLabs) + ' -db ' + str(self.db).replace('.fasta', '')
			os.system(commandx)
			print "BLASTX OK"
		else:
			print "A BLAST XML has existed."

	def tblastn(self, outXML, num_thread):
		print 'Run TBLASTN...'
		path = os.path.join(os.getcwd(), 'pre_xml')
		if not os.path.exists(path):
			os.mkdir(path)
		outXMLabs = os.path.join(path, outXML)
		if not os.path.exists(outXMLabs):
			tcommandn = 'tblastn -evalue 1E-10 -outfmt 5 -query ' + str(self.query) + ' -max_target_seqs 5 -num_threads ' + str(
				num_thread) + ' -out ' + str(outXMLabs) + ' -db ' + str(self.db).replace('.fasta', '')
			os.system(tcommandn)

			print "TBLASTN OK."
		else:
			print "A BLAST XML has existed."

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
                        
        

def Select_long_exons(infile, out, n):
	exon_file = os.path.join(os.getcwd(), 'Reference', infile)
	out_file = os.path.join(os.getcwd(), 'Reference', out)
	with open(out_file, 'w') as fh:
		for records in SeqIO.parse(exon_file, 'fasta'):
			m = len(records.seq)
			if m >= int(n):
				seqs = '>' + str(records.description) + linesep + str(records.seq) + linesep
				fh.write(seqs)
	print infile + ' cut done!'
	return out

def Select_single_copy_exons(infile, in_xml, outfile, Cover, Similarity):
	path = os.getcwd()
	in_xmlabs = os.path.join(path, 'pre_xml', in_xml)
	out_path = os.path.join(path, 'Reference')
	if not os.path.exists(out_path):
		os.mkdir(out_path)
	out_file = os.path.join(out_path, outfile)
	if not os.path.exists(out_file):
		with open(out_file, 'w') as out:
			Dict = SeqIO.to_dict(SeqIO.parse(os.path.join(path, 'Reference', infile), "fasta"), key_function=lambda r: r.description)
			with open(in_xmlabs) as inxml:
				print 'Dictionary Done!'
				
				blast_records = NCBIXML.parse(inxml)
				for blast_record in blast_records:
					query = str(blast_record.query)
					query_len = str(blast_record.query_length)
					n = 0
					q = 0
					h = 0
					two_cover = 1
					two_identify = 1
					record = None
					for alignment in blast_record.alignments:
						n += 1
						for Hsp in alignment.hsps:
							h += 1
							if n == 2 and h == 1:
								two_identify = float(Hsp.identities) / float(Hsp.align_length)
								two_cover = float(Hsp.align_length) / float(query_len)
								break
					if n == 1:
						record = '>' + str(query) + linesep + str(Dict[query].seq) + linesep
						q = 1
					elif n != 1:
						if two_identify > float(Similarity) and two_cover > float(Cover):
							pass
						else:
							record = '>' + str(query) + linesep + str(Dict[query].seq) + linesep
							q = 1
					else:
						pass
					if q == 1:
						out.write(record)
						q += 1
					else:
						pass
				print 'single-copy filter done!'
	else:
		print 'There is a single_copy file!'
	return out_file

def Get_coding_exons(single_copy, xml): 
	queryDict = SeqIO.to_dict(SeqIO.parse(single_copy, "fasta"), key_function=lambda r: r.description)
	print "Single Copy Dictionary Done!"
	path = os.path.join(os.getcwd(), 'Reference')
	if not os.path.exists(path):
		os.mkdir(path)
	fh_1 = os.path.split(single_copy)[-1].replace('_exon', '_cds')
	fh_2 = os.path.split(single_copy)[-1].replace('_exon', '_pep')
	out_cds = os.path.join(path, fh_1)
	out_pep = os.path.join(path, fh_2)
	if not os.path.exists(out_cds):
		with open(out_cds, 'w') as cds_out, open(out_pep, 'w') as pep_out, open(xml, 'r') as fh_xml:
			blast_records = NCBIXML.parse(fh_xml)
			for blast_record in blast_records:
				query = str(blast_record.query)
				pointer = 0
				for alignment in blast_record.alignments:
					pointer = pointer + 1
					for HSP in alignment.hsps:
						if pointer == 1:
							frame = HSP.frame[0]
							start = int(HSP.query_start)
							end = int(HSP.query_end)
							if frame < 0:
								seq = str(queryDict[query].seq[start - 1:end].reverse_complement())
							else:
								seq = str(queryDict[query].seq[start - 1:end])
							seq_pep = Seq(seq).translate()
							if str(seq_pep).count('*') > 1:
								with open('stop_codon_log.txt', 'a') as log:
									lo = str(query) + '\t' + str(str(seq_pep).count('*'))
									log.write(lo + linesep)
							else:
								record_cds = '>' + str(query) + linesep + seq + linesep
								record_pep = '>' + str(query) + linesep + str(seq_pep) + linesep
								cds_out.write(record_cds)
								pep_out.write(record_pep)
								break
	else:
		print 'Coding exon file have existed!'
	return out_cds, out_pep


def Convert_blast_xml_to_txt(species_list, xml_dic):
    
    sp_list = os.path.join(os.getcwd(),species_list)
    reference_list = Get_species_list(sp_list)[0]
    ingroup_list =  Get_species_list(sp_list)[1]
    outgroup_list =  Get_species_list(sp_list)[2]
    taxa_list = Get_species_list(sp_list)[3]


    
    out_dir = '1.HitSeq'
    path = os.path.join(os.getcwd(), out_dir)
    if not os.path.exists(path):
	    os.mkdir(path)
    
    in_path = os.path.join(os.getcwd(),xml_dic)

    for taxa in taxa_list:
        in_file_name = reference_list[0] + '_pep_tblastn_' + str(taxa) + '.xml'
        in_file = os.path.join(in_path, in_file_name)
        out_file = os.path.basename(in_file)
        outname = os.path.join(path, out_file.replace('.xml','_xml_info.txt'))

	if not os.path.exists(outname):
	    print "Start to summary BLAST result for " + out_file + " file."
	    
	    out = open(outname,'w')       
	    out.write('Query\tSubject\tQuery_start\tQuery_end\tSubject_start\tSubject_end\tQuery_seq\tSubject_seq\tIdentity\tQuery_length\tSubject_length\tHit_length\tScore\tFrame\n')
		    
	    result_handle = open(in_file)
	    blast_records = NCBIXML.parse(result_handle)            
	    for blast_record in blast_records:        
		query = str(blast_record.query)
		query_length = str(blast_record.query_letters)


		for alignment in blast_record.alignments: 
		    subject = str(alignment.hit_def)
		    if '<unknown id> ' in subject:                            
                        subject = str(subject.split('<unknown id> ')[-1].split('|')[0])
                    else:
                        subject = str(subject.split(' ')[0])
                            
		    subject_length = int(alignment.length)

		    for HSP in alignment.hsps:                            
                        score = str(HSP.score)
                        query_sta = str(HSP.query_start)
                        query_end = str(HSP.query_end)
                        subject_sta = int(HSP.sbjct_start)
                        subject_end = int(HSP.sbjct_end)
                        identity = HSP.identities
                        query_seq = str(HSP.query)
                        subject_seq = str(HSP.sbjct)
                        hit_length = HSP.align_length
                        identity_percent = float(identity)/hit_length
                        frame = str(HSP.frame[1])
                        if identity_percent >= 0.30:
                                out.write(query+'\t'+subject+'\t'+query_sta+'\t'+query_end+'\t'+str(subject_sta)+'\t'+str(subject_end)+'\t'+query_seq+'\t'+subject_seq+'\t')
                                out.write(str(identity_percent)+'\t'+query_length+'\t'+str(subject_length)+'\t'+str(hit_length)+'\t'+score+'\t'+frame+'\n')
						
	    out.close()
	    result_handle.close()

	else:
	    print "There has exsisted " + outname + " file."


    for taxa in taxa_list:
        txt_in_file_name = reference_list[0]+ '_pep_tblastn_' + str(taxa) + '_xml_info.txt'
        txt_in_file = os.path.join(path, txt_in_file_name)
	filename_sort = os.path.basename(txt_in_file).replace('_info.txt','_sorted.txt')
	filename_simplify = os.path.basename(txt_in_file).replace('_info.txt','_simplified.txt')

	
	out_filename_sort = os.path.join(path, filename_sort)
	out_filename_simplify = os.path.join(path, filename_simplify)

	if not os.path.exists(out_filename_simplify):
		print 'Start to simplify BLAST result for ' + txt_in_file + ' file.'
	
		open_txt_File = open(txt_in_file,'r')
		lines = open_txt_File.readlines()
		LinesDic = {}
		if len(lines) != 1:
                        for i in range(1,len(lines)):
                            infoList = lines[i].split('\n')[0].split('\t')
                            query = infoList[0]
                            subject = infoList[1]
                            subject_start = infoList[4]
                            subject_end = infoList[5]
                            score = infoList[-2]
                            LinesDic[i] = [query,subject,subject_start,subject_end,score]
                        ResultTable = pd.DataFrame(LinesDic).T
                        ResultTable.columns = ['Query','Subject','Subject_start','Subject_end','Score']
                        ResultTable['Score'] = ResultTable['Score'].astype('float')

                        ResultTable_group = ResultTable.groupby(ResultTable['Query'])
                        ResultTable_new = pd.DataFrame(columns=['Query','Subject','Subject_start','Subject_end','Score'])
                        

                        ResultTable_new = ResultTable_group.apply(lambda ResultTable_group: ResultTable_group.sort_values(by=['Score'],ascending = False))
                        ResultTable_new.to_csv(out_filename_sort,sep='\t')
                        
                        ResultTable_drop = ResultTable_new.drop_duplicates(subset=['Query'], keep='first') #
                        ResultTable_drop.set_index(['Query'], inplace=True)
                        ResultTable_drop.to_csv(out_filename_simplify,sep='\t')
                else:
                        pass
	else:
		print "There has exsisted " + out_filename_simplify + " file."

    print 'Simplify information done.'

	    
	
def Extract_hit_seqs(species_list,ingroup_fasta, outgroup_fasta):

    sp_list = os.path.join(os.getcwd(),species_list)
    reference_list = Get_species_list(sp_list)[0]
    ingroup_list =  Get_species_list(sp_list)[1]
    outgroup_list =  Get_species_list(sp_list)[2]
    taxa_list = Get_species_list(sp_list)[3]


    txt_in_path = os.path.join(os.getcwd(),'1.HitSeq')

    All_seq_dic = {}
    ingroup_fasta_file = ingroup_fasta + '/' + '*genome.fasta'
    outgroup_fasta_file = outgroup_fasta + '/' + '*cds.fasta'

        
#Creat seq dict fro ingroups and outgroups
    for ingroup in ingroup_list:
	in_filename = os.path.join(txt_in_path, ingroup + '_hitseq.fasta')

	if not os.path.exists(in_filename):
		ingroup_seq_file = ingroup_fasta + '/' + ingroup + '.fasta'
		print 'Create sequence dictionary for ' + ingroup_seq_file + '.'
		ingroup_seq_records = SeqIO.parse(ingroup_seq_file,'fasta')
		
		dictOfAIngroup = {}        
		for ingroup_seq_record in ingroup_seq_records:
		    dictOfAIngroup[ingroup_seq_record.description.split(' ')[0]] = ingroup_seq_record            
		All_seq_dic[ingroup] = dictOfAIngroup
	else:
		print ingroup + ' hitseq file exsit.'
    print "Ingroup fasta dict done."

    for outgroup in outgroup_list:
	out_filename = os.path.join(txt_in_path, outgroup + '_hitseq.fasta')

	if not os.path.exists(out_filename):        
		outgroup_seq_file = outgroup_fasta + '/' + outgroup + '.fasta'
		print 'Create sequence dictionary for ' + outgroup_seq_file + '.'
		outgroup_seq_records = SeqIO.parse(outgroup_seq_file,'fasta')

		dictOfAOutgroup = {}
		for outgroup_seq_record in outgroup_seq_records:
                        if '<unknown id> ' in outgroup_seq_record.description:
                                dictOfAOutgroup[outgroup_seq_record.description.split('<unknown id> ')[-1].split('|')[0]] = outgroup_seq_record
                        else:
                                dictOfAOutgroup[outgroup_seq_record.description.split(' ')[0]] = outgroup_seq_record
                        
		All_seq_dic[outgroup] = dictOfAOutgroup
	else:
		print outgroup + ' hitseq file exsit.'
    print "Outgroup fasta dict done."



    for taxa in taxa_list:
        txt_file_name = reference_list[0] + '_pep_tblastn_' + str(taxa) + '_xml_simplified.txt'
        txt_file = os.path.join(txt_in_path, txt_file_name)

	taxaname = os.path.basename(txt_file).split('_tblastn_')[-1].split('_xml_simplified.txt')[0] #Rattus_norvegicus_genomic
	filename = taxaname + '_hitseq.fasta'#Rattus_norvegicus_genomic_hitseq.fasta
	outname = os.path.join(txt_in_path, filename)

	if not os.path.exists(outname):

	    print "Start to extract hit sequences for " + taxaname + "."
	    out = open(outname,'w')
	    
	    ID_list = []
	    txt_file_open = open(txt_file,'r')
	    lines = txt_file_open.readlines()
	    for line in lines[1:]: 
		info = line.strip().split('\t')
		subject_name = info[1]
		subject_start = int(info[2])
		subject_end = int(info[3])
		if '|' not in subject_name:
		    hitseq_id = str(subject_name) + '|' + str(subject_start) + '|' +  str(subject_end)#for genome format
		else:
		    hitseq_id = str(subject_name).split('|')[-1]+ '|' + str(subject_start) + '|' +  str(subject_end)# for cds format


		if hitseq_id not in ID_list:                    
		    out.write('>' + str(hitseq_id) + '\n' + str(All_seq_dic[taxaname][subject_name].seq[subject_start - 1 :subject_end].upper()) + '\n')#
		    ID_list.append(hitseq_id)
		else:
		    continue
			  
	    out.close()

	else:
	    print "There has exsisted " + outname + " file."
	    

def Run_MBH_analysis_for_hit_seqs(species_list, reference_dic, hit_seq_dic, num_thread):
        

    sp_list = os.path.join(os.getcwd(),species_list)
    reference_list = Get_species_list(sp_list)[0]
    ingroup_list =  Get_species_list(sp_list)[1]
    outgroup_list =  Get_species_list(sp_list)[2]
    taxa_list = Get_species_list(sp_list)[3]

    out_dir = '1.HitSeq_xml'
    path = os.path.join(os.getcwd(), out_dir)
    if not os.path.exists(path):
	    os.mkdir(path)

    base_path = os.path.join(os.getcwd(), '1.HitSeq') 
    
    #reference pep database
    for reference_seq in glob.glob(reference_dic + '/' + 'single_copy_Long*_pep.fasta'):
	    ref_database_out = os.path.join(reference_dic, os.path.basename(reference_seq).split('.fasta')[0])
	    if not os.path.exists(ref_database_out + '.nhr') or not os.path.exists(ref_database_out + '.nin') or not os.path.exists(ref_database_out + '.nsq'):
		    command1 = "makeblastdb -in " + reference_seq + " -out " + ref_database_out + " -dbtype prot"
		    os.system(command1)

    #hitseq pep database
    for taxa in taxa_list:
        hit_seq_file_name = str(taxa) + '_hitseq.fasta'
        hit_seq_file = os.path.join(hit_seq_dic, hit_seq_file_name)
	hit_database_out = os.path.join(base_path, os.path.basename(hit_seq_file).split('.fasta')[0])

	if not os.path.exists(hit_database_out + '.nhr') or not os.path.exists(hit_database_out + '.nin') or not os.path.exists(hit_database_out + '.nsq'):
		command2 = "makeblastdb -in " + hit_seq_file + " -out " + hit_database_out + " -dbtype nucl"        
		os.system(command2)

    
    for taxa in taxa_list:
        hit_seq_file_name = str(taxa) + '_hitseq.fasta'
        hit_seq_file = os.path.join(hit_seq_dic, hit_seq_file_name)

	for ref_file in glob.glob(reference_dic + '/' + 'single_copy_Long*_pep.fasta'):
		reference_seq = ref_file
		
	tblastn_out_name = str(os.path.basename(reference_seq).split('single_copy_Long_300_')[-1].split('.fasta')[0]) + '_tblastn_' + str(os.path.basename(hit_seq_file).split('.fasta')[0]) + '.xml'
	blastx_out_name = str(os.path.basename(hit_seq_file).split('.fasta')[0]) + '_blastx_' + str(os.path.basename(reference_seq).split('single_copy_Long_300_')[-1].split('.fasta')[0]) + '.xml'
	tblastn_out = os.path.join(path, tblastn_out_name)
	blastx_out = os.path.join(path, blastx_out_name)        

	query_tblastn = os.path.join(os.getcwd(),reference_seq)
	query_blastx = os.path.join(base_path,os.path.basename(hit_seq_file))

	db_tblastn = os.path.join(base_path,os.path.basename(hit_seq_file).split('.fasta')[0])
	db_blastx = os.path.join(os.getcwd(),reference_seq).split('.fasta')[0]
	

	if not os.path.exists(tblastn_out):
	    print 'Run TBLASTN for ' + os.path.basename(hit_seq_file).split('.fasta')[0]
	    command3 = 'tblastn -evalue 1E-10 -outfmt 5 -query '+ query_tblastn + ' -max_target_seqs 5 -num_threads ' + str(num_thread) + ' -out ' + tblastn_out + ' -db ' + db_tblastn
	    os.system(command3)
	else:
	    print 'There has exsisted '+ os.path.basename(hit_seq_file).split('.fasta')[0] + ' tblastn result.'
	

	if  not os.path.exists(blastx_out):
	    print 'Run BLASTX for ' + os.path.basename(hit_seq_file).split('.fasta')[0]
	    command4 = 'blastx -evalue 1E-10 -outfmt 5 -query '+ query_blastx + ' -max_target_seqs 5 -num_threads ' + str(num_thread) + ' -out ' + blastx_out + ' -db ' + db_blastx
	    os.system(command4)
	else:
	    print 'There has exsisted '+ os.path.basename(hit_seq_file).split('.fasta')[0] + ' blastx result.'


def Pipeline(cut_length, number_thread, Cover, Similarity, species_list):
	
	sp_list = os.path.join(os.getcwd(), species_list)
        reference_list = Get_species_list(sp_list)[0]
        ingroup_list =  Get_species_list(sp_list)[1]
        outgroup_list =  Get_species_list(sp_list)[2]
        taxa_list = Get_species_list(sp_list)[3]	

	outfile_cutlength = 'Long_' + str(cut_length) + '_' + reference_list[0] + '_exon.fasta' #Long_300_Rattus_norvegicus_exon.fasta
	Select_long_exons(reference_list[0] + '_exon.fasta', outfile_cutlength, cut_length)
	long = os.path.join(os.getcwd(), 'Reference', outfile_cutlength) #Long_300_Rattus_norvegicus_exon.fasta

	genome = os.path.join(os.getcwd(), 'Reference', str(reference_list[0]) + '_genome.fasta')#Rattus_norvegicus_genome.fasta
	outXML_n = 'Long_' + str(cut_length) + '_blastn_' + str(reference_list[0]) + '_genome.xml'
	file_n = Blast(long, genome)
	file_n.makeblastdb('nucl')
	file_n.blastn(outXML_n, number_thread)

	single = 'single_copy_' + outfile_cutlength #single_copy_Long_300_Rattus_norvegicus_exon.fasta
	single_long = Select_single_copy_exons(reference_list[0] + '_exon.fasta', outXML_n, single, Cover, Similarity)#######################

	peptide = os.path.join(os.getcwd(), 'Reference', str(reference_list[0]) +'_pep.fasta')
	outXML_x = 'single_copy_Long_' + str(cut_length) + '_blastx_' + str(reference_list[0]) + '_pep.xml'
	file_x = Blast(single_long, peptide)
	file_x.makeblastdb('prot')
	file_x.blastx(outXML_x, number_thread)
	xml_x = os.path.join(os.getcwd(), 'pre_xml', outXML_x)
	single_list = Get_coding_exons(single_long, xml_x)

	single_cds = single_list[0]
	single_pep = single_list[1]

        count = 0
	for taxa in taxa_list:
                count += 1		
                if '_cds' in taxa:
                        exon_cds = os.path.join(os.getcwd(), 'Outgroups', str(taxa) + '.fasta')
                        outXML_1 = str(reference_list[0]) + '_pep_tblastn_' + str(taxa) + '.xml'
                        blast_1 = Blast(single_pep, exon_cds)
                        blast_1.makeblastdb('nucl')
                        blast_1.tblastn(outXML_1, number_thread)
                        
                elif '_genome' in taxa:

                        genome_s = os.path.join(os.getcwd(), 'Ingroups', str(taxa) + '.fasta')
                        outXML_1 = str(reference_list[0]) + '_pep_tblastn_' + str(taxa) + '.xml'
                        blast_1 = Blast(single_pep, genome_s)
                        blast_1.makeblastdb('nucl')
                        blast_1.tblastn(outXML_1, number_thread)
                        
        return single_cds
                


def Run_extract_hit_seqs(species_list, num_thread):

	xml_dic = 'pre_xml'
	print 'Start to convert tblastn xml file to txt file.'
	Convert_blast_xml_to_txt(species_list, xml_dic)

	ingroup_fasta = 'Ingroups'
	outgroup_fasta = 'Outgroups'
	print 'Start to extract hit sequences.'
	Extract_hit_seqs(species_list, ingroup_fasta, outgroup_fasta)

	reference_dic = 'Reference'
	hit_seq_dic = '1.HitSeq'
	print 'Start to run MBH analysis between the single copy exons and hit sequences of ingroups and outgroups.'
	Run_MBH_analysis_for_hit_seqs(species_list, reference_dic, hit_seq_dic, num_thread)
	
if __name__ == '__main__':
	dirpath_o = os.path.abspath(os.path.join(os.path.dirname('__file__'), os.path.pardir))
	os.chdir(dirpath_o)
	parser = argparse.ArgumentParser(description='''python Screen_exons_1.py -F species list''')
	
	parser.add_argument('--len', '-L', help='The length cutoff for reference exons', default=int(300))
	parser.add_argument('--thd', '-T', help='Number of threads when running BLAST', default=int(10))
	parser.add_argument('--cov', '-COV', help='Single-copy exon selection criteria: The length coverage between exon and its hit seq', default=float(0.3))
	parser.add_argument('--sim', '-S', help='Single-copy exon selection criteria: The similarity between exon and its hit seq', default=float(0.5))
	parser.add_argument('--file', '-F', help='A text table containing information on the reference species, ingroup species, and outgroup species (filename: species_list.txt)', required=True)
	args = parser.parse_args()
	print args
	try:
		print Pipeline(args.len, args.thd, args.cov, args.sim, args.file)
		Run_extract_hit_seqs(args.file, args.thd)
	except Exception as e:
		print(e)
