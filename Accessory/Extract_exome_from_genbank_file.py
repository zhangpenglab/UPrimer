import os, glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def extractexon(genome):
        startList = []
        endList = []
        descList = []
	for chromosome in SeqIO.parse(genome, "genbank"):
		for feature in chromosome.features:
			if feature.type == "CDS":
				exonLocList = feature.location.parts
				try:
					geneID = feature.qualifiers["db_xref"][0].replace(":", "")
				except KeyError:
					geneID = 'unknown'
				try:
					proteinID = feature.qualifiers["protein_id"][0]
				except KeyError:
					proteinID = ""
				if not proteinID:
					proteinID = "Prot_"+geneID
				product = feature.qualifiers["product"][0]
				if feature.qualifiers.has_key("gene"):
					gene = feature.qualifiers["gene"][0]
				else:
					gene = geneID
				for exonLoc in exonLocList:
					exonSeq = chromosome.seq[exonLoc.start:exonLoc.end]
					exonStrand = exonLoc.strand
					if exonStrand > 0:
						cdsSeq = exonSeq
					else:
						cdsSeq = exonSeq.reverse_complement()
					startList.append(int(exonLoc.start))
					endList.append(int(exonLoc.end))
					cdsRecord = SeqRecord(cdsSeq)
					descList = [proteinID, geneID, gene, chromosome.id, str(exonLoc.start), str(exonLoc.end), str(exonStrand), product]
					cdsRecord.description = "_".join(descList)
					yield cdsRecord

def extractExonInBatch(inDir, outDir):
	if not os.path.isdir(outDir):
		os.makedirs(outDir)
	genomeList = glob.glob(os.path.join(inDir, "*.gbff"))
	for genome in genomeList:
		outBase = os.path.basename(genome).replace('_genomic.gbff', '_exon.fasta')
		print outBase
		outFile = os.path.join(outDir, outBase)
		with open(outFile, "w") as outHandle:
			for exon in extractexon(genome):
				SeqIO.write(exon, outHandle, 'fasta')

if __name__ == '__main__':
	
	inDir = os.path.join(os.getcwd(), 'gbff')
	outDir = os.path.join(os.getcwd(), 'exon')
	extractExonInBatch(inDir, outDir)
	
