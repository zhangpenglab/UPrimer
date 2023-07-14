#!/usr/bin/python
# created by HGC
import os
import _design_primer_pairs_6 as Design_primer_pairs
import _select_primer_pairs_7 as Select_primer_pairs
import argparse
dirpath_o = os.path.abspath(os.path.join(os.path.dirname('__file__'), os.path.pardir))
os.chdir(dirpath_o)
linesep = os.linesep

def Pipeline(dir, con, PIs, fwd, rvs, max, min, species_list, all, inthred2, allthred2, inthred1, allthred1, ConsiderOutgroups, taxa_group, PId):

	Design_primer_pairs.Pipeline(dir, con, fwd, rvs, all, inthred2, allthred2, inthred1, allthred1, max, min, ConsiderOutgroups, PIs)
	Select_primer_pairs.Pipeline(PId, piratio, max, min, species_list, dir, taxa_group)
	return 'The work of designing NPC universal primer sets is complete!'

parser = argparse.ArgumentParser(description = '''python Design_universal_primer_sets_Part_II.py -F species_list.txt -D 4.Candidate_peptide_iden0.5_MSAs_for_primer_design -TG [Name of the target group]''')

                                 
parser.add_argument('--con', '-C', help='The conservation of MSAs', default=float(0.9))
parser.add_argument('--dir', '-D', help='Input a folder containing all candidate peptide MSAs that used for designing primers (Name = 4.Candidate_peptide_iden0.5_MSAs_for_primer_design)', required=True)
parser.add_argument('--max', '-M', help='The maximal length of target amplicons (from F2 to R2)', default=int(700))
parser.add_argument('--min', '-N', help='The minimal length of target amplicons (from F2 to R2)', default=int(100))
parser.add_argument('--PIs','-PIs', help=' The weighting ratio of PCR score and Information score when scoring all possible nested-PCR primers within the same MSA. Increasing this value will reduce the length of the amplified regions. Please set this parameter with caution',default=int(1))
parser.add_argument('--fwd', '-W', help='The length between F1 and F2', default=int(150))
parser.add_argument('--rvs', '-R', help='The length between R1 and R2', default=int(150))
parser.add_argument('--file', '-F', help='Input a txt file containing all species name (filename = species_list.txt)', required=True)
parser.add_argument('--all', '-A', help='The maximal degeneracy level of primers', default=int(8192))
parser.add_argument('--inthred2', '-IT2', help='The mean identity of ingroups (including reference species) in 7AA/8AA F2R2 primer blocks',default=float(0.85))
parser.add_argument('--allthred2', '-AT2', help='The mean identity of all species in 7AA/8AA F2R2 primer blocks',default=float(0.75))
parser.add_argument('--inthred1', '-IT1', help='The mean identity of ingroups (including reference species) in 7AA/8AA F1R1 primer blocks',default=float(0.80))
parser.add_argument('--allthred1', '-AT1', help='The mean identity of all species in 7AA/8AA F1R1 primer blocks',default=float(0.70))
parser.add_argument('--ConsiderOutgroups', '-CO', help='Whether considering outgroup sequences when designing primers. 0: Not considering. 1: Considering', required=True)
parser.add_argument('--TaxaGroup', '-TG', help='Name of the target group (for example: Vertebrata, Bivalvia and so on)', required=True)
parser.add_argument('--PId', '-PId', help='The weighting ratio of PCR score and Information score when scoring highest-scoring nested-PCR primers derived from different MSAs', default=float(1))

args = parser.parse_args()
print args

if __name__ == '__main__':
        print Pipeline(args.dir, args.con, args.PIs, args.fwd, args.rvs, args.max, args.min, args.file, \
                       args.all, args.inthred2, args.allthred2, args.inthred1, args.allthred1, args.ConsiderOutgroups,args.TaxaGroup, args.PId)
