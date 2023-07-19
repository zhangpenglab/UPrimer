#!/usr/bin/python
# -*- coding: UTF-8 -*-

import os, subprocess, sys
import glob
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import Filter_alignments_4 as Tb
import argparse
linesep = os.linesep

class Megacc(object):
    def __init__(self, file, mao):
        self.file = file
        self.path_out_1 = os.path.join(os.getcwd(), '3.1.Megacc_aligned_peptide_MSAs')
        if not os.path.exists(self.path_out_1):
            os.mkdir(self.path_out_1)
        self.mao = os.path.join(os.getcwd(), 'UPrimer', mao)
        self.outfile = str(os.path.split(file)[-1]).replace('.fasta', '')
    
    def mega_align_prot(self):
        self.out_file = os.path.join(self.path_out_1, self.outfile)
        if not os.path.exists(self.out_file + '.meg'):
            command = "megacc -a " + str(self.mao) + " -d " + str(self.file) + " -o " + str(self.out_file) + '.meg'
            os.system(command)

        else:
            print self.outfile + ' protein megacc align have done !'
        return str(self.out_file)
        
    def tra_meg_fasta(self, meg, error_seq):
        out_path = os.path.split(meg)
        if not os.path.exists(meg):
            out_error = os.path.join(out_path[0], error_seq)
            with open(out_error, 'a') as out_error_seq:
                out_error_seq.write(out_path[1].replace('.meg', '.fasta') + linesep)
        elif 'pep' in str(out_path[1]):
            out_file_1 = os.path.join(out_path[0], 'Tra_fasta_peptide')
            if not os.path.exists(out_file_1):
                os.mkdir(out_file_1)
            out_file_2 = os.path.join(out_file_1, out_path[1].replace('meg', 'fasta'))
            if not os.path.exists(out_file_2):
                with open(out_file_2, 'w') as fh_out:
                    with open(meg, 'r') as fh:
                        for lines in fh.readlines()[4:]:
                            if '#' in lines:
                                gene = str(lines).replace('#', ">") + linesep
                                fh_out.write(gene)
                            else:
                                fh_out.write(str(lines).replace('?', 'X'))
                return out_file_2
            else:
                print str(out_path[1]).replace('.fasta', '') + ' protein have been translated !'
                return out_file_2

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
    
def Realign_MSAs_after_removing_rogue_taxon():
    remove_pep = glob.glob(os.path.join(os.getcwd(), '3.2.After_remove_rogue_taxa_peptide_MSAs', '*.fasta'))
    origin_pep = glob.glob(os.path.join(os.getcwd(), '2.Peptide_orthologous_groups', '*.fasta'))
    path_out = os.path.join(os.getcwd(), '3.3.Re_align_peptide_MSAs')
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    dic = {}
    for file_1 in remove_pep:
        records = SeqIO.parse(file_1, 'fasta')
        taxa_list = []
        gene_name = os.path.split(file_1)[-1].replace('Remove_rogue_taxa_', '').replace('.fasta', '')
        for record in records:
            taxa_list.append('_'.join(str(record.description).split('_')[0:2]))
        dic.update({gene_name: taxa_list})
    out_list = []
    for file_2 in origin_pep:
        Gene = os.path.split(file_2)[-1].replace('.fasta', '')
        if Gene in dic.keys():
            outname = os.path.join(path_out, 'Re_align_' + str(os.path.split(file_2)[-1]))
            if not os.path.exists(outname):
                with open(outname, 'w') as out:
                    for record in SeqIO.parse(file_2, 'fasta'):
                        ID = str('_'.join(str(record.description).split('_')[0:2]))
                        if ID in dic[Gene]:
                            SeqIO.write(record, out, 'fasta')
                out_list.append(outname)
            else:
                out_list.append(outname)
        else:
            pass
    return out_list

def Cut_two_sides_of_MSAs(fh, taxa_information, out_dir):#
    path_list = os.path.split(fh)
    out_file_1 = path_list[1]
    path_out = os.path.join(path_list[0], out_dir)
    if not os.path.exists(path_out):
        os.mkdir(path_out)
        
    out_file = os.path.join(path_out, out_file_1)
    if not os.path.exists(out_file):
        start = 0
        end = 1
        alignment = AlignIO.read(fh, "fasta")
        line_1 = str(out_file_1) + "_Taxa_Number: %i" % len(alignment)
        with open(fh, 'r') as fh_in:
            refer = fh_in.readline().replace('>', '').strip()

        for record in alignment:
            if refer == record.description:
                for i in range(0, len(record.seq)):
                    if record.seq[i] != '-' and record.seq[i] != 'X':
                        start = i
                        break
                for j in range(len(record.seq)-1, -1, -1):
                    if record.seq[j] != '-' and record.seq[j] != 'X':
                        end = j + 1
                        break
            else:
                pass
            
        line_2 = 'Cut alignment range:' + str(start) + '\t' + str(end)
        path_log = os.path.join(path_out, taxa_information)
        with open(path_log, 'a') as log:
            log.write(str(line_1) + linesep + str(line_2) + linesep)
        my_alignments = alignment[:, start:end]
        AlignIO.write(my_alignments, out_file, "fasta")
        return out_file
    else:
        print 'Cutted (two sides of MSA) '+ str(out_file_1) + ' have existed !'
        return out_file

def Calculate_identity_of_two_seqs(record1, record2):
    count = 0

    case1 = 0
    case2 = 0
    case3 = 0
    
    list1 = list(str(record1.seq))
    list2 = list(str(record2.seq))
    bad_list = list('X-')
    aa_list = list('GAVIPFYWTCMNQDEKHLSRBZJ')
    
    for i in zip(list1, list2):
        if i[1] == i[0] and i[1] in aa_list and i[0] in aa_list:# AA and  --
            count += 1

        if i[1] in aa_list and i[0] in aa_list:
            case1 += 1
        if i[1] in aa_list and i[0] == '-':
            case2 += 1
        if i[1] == '-' and i[0] in aa_list:
            case3 += 1
        times = case1 +  case2 + case3

    if times == 0: ###avoid 0/0 error
       times = 1
    Identity = float(count) / times
    return Identity*100

def Calculate_Simpsonindex_and_GapPercentage(seq):
    seq = str(list(seq))
    length = 1
    length_list = list()
    converted_seq = ''
    for i in range(0,len(seq)):
        if seq[i]=='-':
            converted_seq = converted_seq + '0'
        else:
            converted_seq = converted_seq + '1'
    gap_percentage = float(converted_seq.count('0'))/float(len(seq))#calculate gap percentage

    for i in range(0,len(seq)-1): 
        if converted_seq[i] == converted_seq[i+1]:
            length += 1
        else:
            length_list.append(length)
            length = 1
    length_list.append(length)
    Simpsonindex_score = 0
    for block_len in length_list:
        Simpsonindex_score = Simpsonindex_score + float(block_len**2)/float(len(seq)**2) #calculate Simpson value
    return Simpsonindex_score,gap_percentage


def Remove_rogue_taxon(file, species_list, Rogue_taxa_remove, select_scale, remove_rigor):

    sp_list = os.path.join(os.getcwd(),species_list)
    reference_list = Get_species_list(sp_list)[0]
    ingroup_list =  Get_species_list(sp_list)[1]
    outgroup_list =  Get_species_list(sp_list)[2]
    taxa_list = Get_species_list(sp_list)[3]
    
##    taxa_list = []
    path_out = os.path.join(os.getcwd(), '3.2.After_remove_rogue_taxa_peptide_MSAs')
    out_name = 'Remove_rogue_taxa_' + os.path.split(file)[1]
    out_file = os.path.join(path_out, out_name)

    if not os.path.exists(path_out):
        os.mkdir(path_out)

    if not os.path.exists(out_file):
##        with open(species_list, 'r') as fh:
##            for taxa in fh.readlines()[1:]:
##                taxa_list.append(taxa.strip().replace('.fasta', '').replace('exon', 'cds'))

        records = list(SeqIO.parse(file, 'fasta'))# Description: Euschistus_heros_genomic_RCWM01000406.1|376006|376476
        remove_taxa = []

        Identity_dic = {}
        Identity_list = []

        Simpsonindex_dic = {}
        Simpsonindex_list = []
        
        GapPercentage_dic = {}
        GapPercentage_list = []
        Gene_sequence = os.path.split(file)[1].replace('.fasta', '')#msa name

        for i in range(0, len(records)):
            Simpsonindex_value = Calculate_Simpsonindex_and_GapPercentage(records[i])[0]
            GapPercentage_value = Calculate_Simpsonindex_and_GapPercentage(records[i])[1]

            Simpsonindex_dic[records[i].description] = Simpsonindex_value
            GapPercentage_dic[records[i].description] = GapPercentage_value

            Simpsonindex_list.append(Simpsonindex_value)
            GapPercentage_list.append(GapPercentage_value)

        
        average_gap = np.mean(GapPercentage_list)
        gap_std = np.std(GapPercentage_list)

        average_Simpsonindex = np.mean(Simpsonindex_list)
        Simpsonindex_std = np.std(Simpsonindex_list)

        #calculate taxa mean identity
        for i in range(0, len(records)):
            I = 0
            times = len(records) - 1
            for j in range(0, len(records)):
                if i != j:
                    I = I + Calculate_identity_of_two_seqs(records[i], records[j])
            Taxa_Ave_Identity = I / times
            Identity_dic.update({records[i].description: Taxa_Ave_Identity})

        for taxa in Identity_dic:
            taxa_identity_value = Identity_dic[taxa]
            Identity_list.append(taxa_identity_value)

        average_taxa_identity = np.mean(Identity_list)
        taxa_identity_std = np.std(Identity_list)

        #remove rogue taxa seqs
        if min(Identity_list) <= float(remove_rigor):
            for taxa in Identity_dic:
                taxa_identity_value = Identity_dic[taxa]
                if taxa_identity_value < 30:
                    remove_taxa.append(taxa)
                elif taxa_identity_value < average_taxa_identity - float(select_scale) * taxa_identity_std:
                    remove_taxa.append(taxa)
                else:
                    pass
        else:
            for taxa in Identity_dic:
                taxa_identity_value = Identity_dic[taxa]
                if taxa_identity_value < 30:
                    remove_taxa.append(taxa)
                

        #remove seqs with high gap prtcentage
        for taxa in GapPercentage_dic:
            taxa_GapPercentage_value = GapPercentage_dic[taxa]
            if taxa_GapPercentage_value >= 0.6:
                remove_taxa.append(taxa)

            
        #Get unique remove_taxa list
        remove_taxa = np.unique(remove_taxa)

        #write rogue taxa record
        if len(remove_taxa) != 0:
            txt_out = open(os.path.join(os.path.split(file)[0], Rogue_taxa_remove), 'a')
            for i in range(0, len(remove_taxa)):
                line_1 = ('\t').join([str(Gene_sequence), str(remove_taxa[i]), str(Identity_dic[remove_taxa[i]]),\
                                     str(GapPercentage_dic[remove_taxa[i]]),str(Simpsonindex_dic[remove_taxa[i]])])
                txt_out.write(line_1 + linesep)
            
        #write non-rogue taxa sequences
        new_records = []
        new_records.append(records[0])
        for record in records:                
            if record.description not in remove_taxa and 'exon' not in record.description:###########
                new_records.append(record)

        if len(new_records) >= 1:
            with open(out_file, 'w') as out:
                SeqIO.write(new_records, out, 'fasta')
                
                
    else:
        print 'There are removed rogue taxa files!'  

def Run_pal2nal_aligning(pep, nucl, out, outgroupNum, ingroupNum, totaltaxanum, ConsiderOutgroups):
    print os.path.split(nucl)[1], 'begin re_align!'
    arg = "pal2nal.pl " + pep + ' ' + nucl + ' -output fasta ' + '>' + out
    os.system(arg)
    after_cut = Cut_two_sides_of_MSAs(out, '2_Taxa_information_nucl.txt', 'Cut_Alignment_2')
    after_trim = Tb.Trim_nucleotide_MSAs(after_cut, '3.5.Refined_nucleotide_aligned_MSAs', 30, 1.5, outgroupNum, ingroupNum, totaltaxanum, ConsiderOutgroups)[0]

    
    
def Calculate_taxa_number(file):#calculate_num
    outgroup = 0
    ingroup = 0
    for record in SeqIO.parse(open(file), 'fasta'):
        if '_cds_' in str(record.id):
            outgroup += 1
        elif '_genomic_' in str(record.id):
            ingroup += 1
    return outgroup, ingroup

def Realign_nucleotide_MSAs(outgroupNum, ingroupNum, totaltaxanum, ConsiderOutgroups):
    remove_pep = glob.glob(os.path.join(os.getcwd(), '3.1.Megacc_aligned_peptide_MSAs', 'Tra_fasta_peptide', 'Re_align_*.fasta'))
    origin_nucl = glob.glob(os.path.join(os.getcwd(), '2.Nucleotide_orthologous_groups', '*.fasta'))
    path_out = os.path.join(os.getcwd(), '3.4.Pal2nal_nucleotide_aligned_MSAs')
    path_out_pal2nal = os.path.join(path_out, 'Re_align_nucleotide')
    
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    if not os.path.exists(path_out_pal2nal):
        os.mkdir(path_out_pal2nal)
    
    dic = {}
    for file_1 in remove_pep:
        records = SeqIO.parse(file_1, 'fasta')
        taxa_list = []
        gene_name = os.path.split(file_1)[-1].replace('Re_align_', '').replace('pep.fasta', 'nucl')
        for record in records:
            taxa_list.append('_'.join(str(record.description).split('_')[0:2]))
        dic.update({gene_name: taxa_list})
    
    for file_2 in origin_nucl:
        Gene = os.path.split(file_2)[-1].replace('.fasta', '')
        if Gene in dic.keys():
            outname = os.path.join(path_out, 'Re_align_' + str(os.path.split(file_2)[1]))
            if not os.path.exists(outname):
                with open(outname, 'w') as out:
                    for record in SeqIO.parse(file_2, 'fasta'):
                        ID = ('_'.join(str(record.description).split('_')[0:2]))
                        if ID in dic[Gene]:
                            SeqIO.write(record, out, 'fasta')
            else:
                pass
            out_align = os.path.join(path_out_pal2nal, 'Re_align_' + str(os.path.split(file_2)[1]))
            align_pep = os.path.join(os.getcwd(), '3.1.Megacc_aligned_peptide_MSAs', 'Tra_fasta_peptide',
                                     'Re_align_' + str(os.path.split(file_2)[1]).replace('_nucl', '_pep'))
            Run_pal2nal_aligning(align_pep, outname, out_align,outgroupNum, ingroupNum, totaltaxanum, ConsiderOutgroups)
        else:
            pass

if __name__ == '__main__':
    dirpath_o = os.path.abspath(os.path.join(os.path.dirname('__file__'), os.path.pardir))
    os.chdir(dirpath_o)
    parser = argparse.ArgumentParser(description='''python Make_codon_alignment_3.py -F species_list.txt -D 2.Peptide_orthologous_groups -CO [1: Consider outgroups;0: Do not consider outgroups]''')
                                     
    parser.add_argument('--outgroupNum', '-ON', help='The minimal number of outgroups', default=int(1))
    parser.add_argument('--ingroupNum', '-IN', help='the minimal number of ingroups', default=int(1))
    parser.add_argument('--totaltaxanum', '-TN', help='The minimal number of all species', default=int(3))
    parser.add_argument('--dir', '-D', help='Input a folder including all unaligned peptide MSAs (Folder name: 2.Peptide_orthologous_groups)', required=True)
    parser.add_argument('--file', '-F', help='Name of a txt file that contains the name of all used genome resource files (for example: species_list.txt)', required=True)
    parser.add_argument('--tril', '-l', help='The permissible length cutoff of amino acid during trimming sequences. This parameter usually does not need to be adjusted', default=int(10))
    parser.add_argument('--trit', '-t', help='The ratio of gap length and AA length during trimming sequences. This parameter usually does not need to be adjusted', default=float(1.5))
    parser.add_argument('--ConsiderOutgroups', '-CO', help='Whether outgroups are considered when designing universal primers. 0: Do not consider; 1: Consider.', required=True)
    args = parser.parse_args()
    print args
