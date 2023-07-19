#!/usr/bin/python
# created by HGC

import os
import Screen_exons_1 as Screen_exons
import Append_orthologous_sequences_2 as Append_orthologous_sequences
import Find_blocks_5 as Find_blocks
import argparse
linesep = os.linesep
dirpath_o = os.path.abspath(os.path.join(os.path.dirname('__file__'), os.path.pardir))
os.chdir(dirpath_o)

def Pipeline(cut_length, thread, Cover, Similarity, species_list, select_scale, remove_rigor, Expect_identity, outgroupNum, ingroupNum, totaltaxanum, trim_length, trim_times):

	print 'Identify single copy and long exons for reference species...'
	single_cds = Screen_exons.Pipeline(cut_length, number_thread, Cover, Similarity, species_list)
	print 'Extract hit sequences from ingroup and outgroups...'    
	Screen_exons.Run_extract_hit_seqs(species_list, number_thread)
	print 'Generate multiple sequence alignments (MSAs) and refine them...'    
	Append_orthologous_sequences.Pipeline(single_cds, species_list, cut_length, select_scale, remove_rigor, outgroupNum, ingroupNum, totaltaxanum, trim_length, trim_times)
	print 'Select alignments suitable for primer design...'
	Find_blocks.Pipeline(species_list, Expect_identity)
	return 'The work of making MSAs suitbale for primer design is complete!'

parser = argparse.ArgumentParser(description = 'python Make_MSAs_suitable_for_primer_design_Part_I.py '
                                                '-L The length cutoff for reference exons '                                 
                                                '-T Number of threads for running BLAST '                                 
                                                '-COV The length coverage between exon and its hit seq '
                                                '-S The similarity between exon and its hit seq '
                                                '-F Input a txt file containing all species name '                                 
                                                '-SS The criteria used in removing the rogue seqs '
                                                '-RR A sequence identity cutoff used in removing the rogue seqs '
                                                '-E The expect identity of all species in 7AA/8AA blocks in MSAs '
                                                '-ON The minimal number of outgroups '
                                                '-IN The minimal number of ingroups '
                                                '-TN The minimal number of all species '
                                                '-l The permissible length cutoff of AA during trimming sequences '
                                                '-t The ratio of gap length and AA length during trimming sequences ')

                                 
parser.add_argument('--len', '-L', help='The length cutoff for reference exons', default=int(300))
parser.add_argument('--thd', '-T', help='Number of threads for running BLAST', default=int(10))
parser.add_argument('--cov', '-COV', help='Single-copy exon criteria: The length coverage between exon and its hit seq', default=float(0.3))
parser.add_argument('--sim', '-S', help='Single-copy exon criteria: The similarity between exon and its hit seq', default=float(0.5))
parser.add_argument('--file', '-F', help='A text table containing information on the reference species, ingroup species, and outgroup species (filename: species_list.txt)', required=True)
parser.add_argument("--SelectScale", "-SS", help="A criteria used in removing the rogue seqs, larger and more strict", default=float(1.1))
parser.add_argument("--RemoveRigor", "-RR", help="A sequence identity cutoff used in removing the rogue seqs", default=float(55))
parser.add_argument('--exp', '-E', help='The expect identity of all species in 7AA/8AA blocks', default=float(0.5))
parser.add_argument('--outgroupNum', '-ON', help='The minimal number of outgroups', default=int(1))
parser.add_argument('--ingroupNum', '-IN', help='The minimal number of ingroups', default=int(1))
parser.add_argument('--totaltaxanum', '-TN', help='The minimal number of all species', default=int(3))
parser.add_argument('--trml', '-l', help='The permissible length cutoff of AA during trimming sequences', default=int(10))
parser.add_argument('--trmt', '-t', help='The ratio of gap length and AA length during trimming sequences', default=float(1.5))
args = parser.parse_args()
print args

if __name__ == '__main__':
	try:
		print Pipeline(args.len, args.thd, args.cov, args.sim, args.file, args.SelectScale, args.RemoveRigor, args.exp, args.outgroupNum, args.ingroupNum, args.totaltaxanum, args.trml, args.trmt)
	except Exception as e:
		print(e)
