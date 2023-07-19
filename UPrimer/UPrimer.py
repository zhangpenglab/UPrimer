#!/usr/bin/python
# -*- coding: UTF-8 -*-

import os
import Screen_exons_1 as Screen_exons
import Append_orthologous_sequences_2 as Append_orthologous_sequences
import Find_blocks_5 as Find_blocks
import _design_primer_pairs_6 as Design_primer_pairs
import _select_primer_pairs_7 as Select_primer_pairs
import argparse

linesep = os.linesep
dirpath_o = os.path.abspath(os.path.join(os.path.dirname('__file__'), os.path.pardir))
os.chdir(dirpath_o)

def Pipeline(cut_length, number_thread, Cover, Similarity, species_list, select_scale, remove_rigor, Expect_identity, outgroupNum, ingroupNum, totaltaxanum, trim_length,
             trim_times, con, PIs, fwd, rvs, max, min, all, inthred2, allthred2, inthred1, allthred1, ConsiderOutgroups,taxa_group, PId):

	print 'Identify single copy and long exons for the reference species...'        
	single_cds = Screen_exons.Pipeline(cut_length, number_thread, Cover, Similarity, species_list)
	print 'Extract BLAST hit sequences from ingroups and outgroups and proceed the MHB analysis...'        
	Screen_exons.Run_extract_hit_seqs(species_list, number_thread)
	print 'Generate multiple sequence alignments (MSAs) and refine them...'
	Append_orthologous_sequences.Pipeline(single_cds, species_list, cut_length, select_scale, remove_rigor, trim_length, trim_times, outgroupNum, ingroupNum, totaltaxanum, ConsiderOutgroups)
	print 'Select alignments suitable for primer design...'
	pep_dir = Find_blocks.Pipeline(species_list, Expect_identity)
	print 'Design universal primers...'
	Design_primer_pairs.Pipeline(pep_dir, con, fwd, rvs, all, inthred2, allthred2, inthred1, allthred1, max, min, ConsiderOutgroups, PIs)
	Select_primer_pairs.Pipeline(PId, max, min, species_list, pep_dir, taxa_group)

	return 'The work of designing NPC universal primer sets is complete!'

parser = argparse.ArgumentParser(description ='''python UPrimer.py -F species_list.txt -CO [1: Consider outgroups;0: Do not consider outgroups] -TG [Name of the target group]''')


parser.add_argument('--len', '-L', help='The length cutoff for reference exons', default=int(300))
parser.add_argument('--thd', '-T', help='Number of threads for running BLAST', default=int(10))
parser.add_argument('--cov', '-COV', help='Single-copy exon criteria: The length coverage between exon and its hit seq', default=float(0.3))
parser.add_argument('--sim', '-S', help='Single-copy exon criteria: The similarity between exon and its hit seq', default=float(0.5))
parser.add_argument('--file', '-F', help='A text table containing information on the reference species, ingroup species, and outgroup species (filename: species_list.txt)', required=True)
parser.add_argument("--SelectScale", "-SS", help="A criteria used in removing the rogue seqs, larger and more strict", default=float(1.1))
parser.add_argument("--RemoveRigor", "-RR", help="A sequence identity cutoff used in removing the rogue seqs", default=float(55))
parser.add_argument('--exp', '-E', help='The expect identity of all species in 7AA/8AA blocks in MSAs', default=float(0.5))
parser.add_argument('--outgroupNum', '-ON', help='The minimal number of outgroups', default=int(1))
parser.add_argument('--ingroupNum', '-IN', help='The minimal number of ingroups', default=int(1))
parser.add_argument('--totaltaxanum', '-TN', help='The minimal number of all species', default=int(3))
parser.add_argument('--tril', '-l', help='The permissible length cutoff of AA during trimming sequences', default=int(10))
parser.add_argument('--trit', '-t', help='The ratio of gap length and AA length during trimming sequences', default=float(1.5))
parser.add_argument('--con', '-C', help='The conservation of MSAs', default=float(0.9))
parser.add_argument('--PIs','-PIs', help=' The weighting ratio of PCR score and Information score when scoring all possible nested-PCR primers within the same MSA. Increasing this value will reduce the length of the amplified regions. Please set this parameter with caution',default=int(1))
parser.add_argument('--fwd', '-W', help='The length between F1 and F2', default=int(150))
parser.add_argument('--rvs', '-R', help='The length between R1 and R2', default=int(150))
parser.add_argument('--max', '-M', help='The maximal length of target amplicons (from F2 to R2)', default=int(700))
parser.add_argument('--min', '-N', help='The minimal length of target amplicons (from F2 to R2)', default=int(100))
parser.add_argument('--all', '-A', help='The maximal degeneracy level of primers', default=int(8192))
parser.add_argument('--inthred2', '-IT2', help='The mean identity of ingroups (including reference species) in 7AA/8AA F2R2 primer blocks',default=float(0.85))
parser.add_argument('--allthred2', '-AT2', help='The mean identity of all species in 7AA/8AA F2R2 primer blocks',default=float(0.75))
parser.add_argument('--inthred1', '-IT1', help='The mean identity of ingroups (including reference species) in 7AA/8AA F1R1 primer blocks',default=float(0.80))
parser.add_argument('--allthred1', '-AT1', help='The mean identity of all species in 7AA/8AA F1R1 primer blocks',default=float(0.70))
parser.add_argument('--ConsiderOutgroups', '-CO', help='Whether considering outgroup sequences when designing primers. 0: Not considering. 1: Considering', required=True)
parser.add_argument('--TaxaGroup', '-TG', help='Input the taxa group name used for outputing primers, such as vertebrate, mammalian and so on.', required=True)
parser.add_argument('--PId', '-PId', help='The weighting ratio of PCR score and Information score when scoring highest-scoring nested-PCR primers derived from different MSAs', default=float(1))
args = parser.parse_args()
print args

if __name__ == '__main__':
        print Pipeline(args.len, args.thd, args.cov, args.sim, args.file, args.SelectScale, args.RemoveRigor, args.exp, \
                       args.outgroupNum, args.ingroupNum, args.totaltaxanum, args.tril,
                       args.trit, args.con, args.PIs, args.fwd, args.rvs, args.max, args.min, args.all, \
                       args.inthred2, args.allthred2, args.inthred1, args.allthred1, args.ConsiderOutgroups,args.TaxaGroup, args.PId)
