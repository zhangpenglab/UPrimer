#!/usr/bin/python
# -*- coding: UTF-8 -*-

'''

Before running this script, the user needs to prepare 3 folders:

1)A folder contains assembled contigs of each sample in '.fasta' format. 
Then  the user needs to offer the ABSOLUTE PATH of this folder to '-c' argument when running this script.
Files should be named in 'SampleName.fasta' format, for example, 'LepSample001.fasta','LepSample002.fasta'.

2)A folder named 'reference' that contains two sequence files.
One is the reference nucleotide sequence file and it needs to be named 'ref_nuc.fasta'.
Another is the reference peptide sequence file. It needs to be named 'ref_pro.fasta'.
The format of reference sequence id is '>Lep1'.

3)A folder named 'pep_for_exonerate' that contains the file 'ref_pro.fasta'


In brief, the whole process procedure can be divided into 6 steps: 
1) TBLASTN to search target gene from assembled contigs
2) Get best hit from blast result
3) Sort sequences by gene
4) Run exonerate to separate exon and intron
5) Add reference sequences
6) Filter by gene length and translate nucleotide sequences

This script has been tested with Python 2.7.15, Biopython 1.74, BLAST 2.12.0+ and exonerate v2.4.0 on the liunx operating system.

'''

import os,shutil,subprocess,re,glob,sys,getopt,time,argparse
from Bio import SeqIO,AlignIO
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def makeblastdb(contigDir):
    pathout_blastdb=os.path.join(contigDir,os.path.pardir,'1.Blast_result','blast_db_file')
    if not os.path.exists(pathout_blastdb):
        os.makedirs(pathout_blastdb)
    list_contig = os.listdir(contigDir)
    for i in range(0,len(list_contig)):
        contig_file = os.path.join(contigDir,list_contig[i])
        if os.path.isfile(contig_file) and contig_file.split('.')[-1]=='fasta':
            sample=os.path.basename(contig_file).split('.')[0]
            if not os.path.exists(os.path.join(pathout_blastdb,sample)+".nhr") or \
                   not os.path.exists(os.path.join(pathout_blastdb,sample)+".nin") or \
                   not os.path.exists(os.path.join(pathout_blastdb,sample)+".nsq"):
                makeblastdbCommand='makeblastdb -in '+contig_file+' -dbtype nucl '+\
                                    '-out '+os.path.join(pathout_blastdb,sample)
                os.system(makeblastdbCommand)
            else:
                print sample+' blastdb existed already.'
                
def tblastn(contigDir):
    pathout_blastdb=os.path.join(contigDir,os.path.pardir,'1.Blast_result','blast_db_file')
    pathout_blastresult=os.path.join(contigDir,os.path.pardir,'1.Blast_result')
    if not os.path.exists(pathout_blastresult):
        os.makedirs(pathout_blastresult)  
    list_contig = os.listdir(contigDir)
    for i in range(0,len(list_contig)):
        contig_file = os.path.join(contigDir,list_contig[i])
        if os.path.isfile(contig_file) and contig_file.split('.')[-1]=='fasta':
            sample=os.path.basename(contig_file).split('.')[0]
            reference=os.path.join(contigDir,os.path.pardir,'reference','ref_pro.fasta')
            blast_xml=os.path.join(pathout_blastresult,sample)+'.xml'
            if not os.path.exists(blast_xml):
                tblastnCommand='tblastn -db '+os.path.join(pathout_blastdb,sample)+' -query '+\
                        reference+' -outfmt 5 -max_target_seqs 1 -num_threads 40  -evalue 1e-5'+\
                        ' -out '+blast_xml
                os.system(tblastnCommand)
            else:
                print sample+' tblastn done before.'

def blast_result_xml_to_txt(contigDir):
    pathout_blastresult=os.path.join(contigDir,os.path.pardir,'1.Blast_result')
    xml_filenames=glob.glob(os.path.join(pathout_blastresult,'*.xml'))
    for filename in xml_filenames:
        outname=filename.replace('.xml','_xml_info.txt')
        out=open(outname,'w')
        out.write('Subject\tQuery\tQuery_start\tQuery_end\tSubject_start\tSubject_end\tQuery_seq\tSubject_seq\tSimilarity\tIdentity\tQuery_length\tSubject_length\tHit_length\tScore\tframe\n')
        result_handle=open(filename)
        blast_records=NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            query=str(blast_record.query)
            query_length=str(blast_record.query_letters)
            for alignment in blast_record.alignments:
                subject=str(alignment.title.split(' ')[1])
                subject_length=str(alignment.length)
                for HSP in alignment.hsps:
                    score=str(HSP.score)
                    query_sta=str(HSP.query_start)
                    query_end=str(HSP.query_end)
                    subject_sta=str(HSP.sbjct_start)
                    subject_end=str(HSP.sbjct_end)
                    identity=HSP.identities
                    query_seq=str(HSP.query)
                    subject_seq=str(HSP.sbjct)
                    hit_length=int(HSP.align_length)
                    identity_percent=float(identity)/hit_length
                    frame=str(HSP.frame[1])
                    out.write(subject+'\t'+query+'\t'+query_sta+'\t'+query_end+'\t'+subject_sta+'\t'+subject_end+'\t'+query_seq+'\t'+subject_seq+'\t')
                    out.write(str(identity_percent)+'\t'+str(identity)+'\t'+query_length+'\t'+subject_length+'\t'+str(hit_length)+'\t'+score+'\t'+frame+'\n')
        out.close()
        result_handle.close()

def get_best_Hit_info(contigDir):
    pathout=os.path.join(contigDir,os.path.pardir,'2.Filter_result')
    if not os.path.exists(pathout):
        os.makedirs(pathout)
    filenames=glob.glob(os.path.join(contigDir,os.path.pardir,'1.Blast_result','*.txt'))
    for filename in filenames:
        taxon=os.path.basename(filename).split('_')[0]
        refername=os.path.join(contigDir,os.path.pardir,'reference','ref_pro.fasta')
        referdic={}
        genedic={}   
        for record in SeqIO.parse(refername,'fasta'):
            referdic[record.id]=[]
            genedic[record.id]={}
        for line in open(filename).readlines()[1:]:
            temp=line.split('\t')
            Subject=temp[0]
            Qgene=temp[1]
            Score=float(temp[-2])
            Identity=float(temp[-7])
            if Subject not in genedic[Qgene]:
                genedic[Qgene][Subject]=Score
                referdic[Qgene].append(Subject)
            else:
                if Identity > 0.6:
                    genedic[Qgene][Subject]+=Score
        newdic={}
        for key in genedic:
            if genedic[key]!={}:
                if len(genedic[key])==1:
                    newdic[key]=genedic[key].keys()[0]
                else:
                    BestHSP=referdic[key][0]
                    BestScore=genedic[key][BestHSP]
                    newdic[key]=BestHSP
                    contiglist=sorted(genedic[key].items(),key=lambda item:item[1],reverse=True)
        outname=os.path.basename(filename).replace('.txt','_filter.txt')
        out=open(os.path.join(pathout,outname),'w')
        count=0
        GenLis=[]
        ContigLis=[]
        for line in open(filename):
            if count==0:
                out.write(line)
                count=1
            else:
                Subjectcontig=line.split('\t')[0]
                Query=line.split('\t')[1]
                identity=float(line.split('\t')[-7])
                GenCover=float(line.split('\t')[-3])/float(line.split('\t')[-5])
                if Query in newdic and Subjectcontig==newdic[Query]:
                    if identity >= 0.6:
                        out.write(line)
                        GenLis.append(Query)
                        ContigLis.append(Subjectcontig)
        out.close()

def get_hit_seqs(contigDir):
    pathout=os.path.join(contigDir,os.path.pardir,'2.Filter_result','Hit_seqs')
    if not os.path.exists(pathout):
        os.makedirs(pathout)
    count=0
    refernames=glob.glob(os.path.join(contigDir,os.path.pardir,'2.Filter_result','*.txt'))
    for refername in refernames:
        taxon=os.path.basename(refername).split('_')[0]
        genedic={}
        framedic={}
        for line in open(refername).readlines()[1:]:
            temp=line.split('\t')
            Subject=temp[0]
            Qgene=temp[1]
            start=float(temp[4])
            end=float(temp[5])
            frame=int(temp[-1])
            newname=Qgene+'_'+Subject
            if newname not in genedic:
                genedic[newname]=[]
                genedic[newname].append(start)
                genedic[newname].append(end)
            else:
                genedic[newname].append(start)
                genedic[newname].append(end)
            if newname not in framedic:
                framedic[newname]=frame
            else:
                if frame*framedic[newname]<0:
                    print refername,newname,line,framedic[newname],frame,'Wrong frame, please check!'
        
        filename=os.path.join(contigDir,taxon+'.fasta')
        seqdic=SeqIO.to_dict(SeqIO.parse(filename,'fasta'))
        outname=os.path.join(pathout,taxon+'_filtered_exon.fasta')
        out=open(outname,'w')
        for key in genedic:
            contigname=key.replace(key.split('_')[0]+'_','')
            genedic[key].sort()
            sta=genedic[key][0]
            en=genedic[key][-1]
            rename=key.replace(':','_')+'_'+taxon
            newseq=str(seqdic[contigname].seq)[int(sta)-1:int(en)]###get hit length,not whole sequence
            if framedic[key]>0:
                record=SeqRecord(Seq(newseq),id=rename,description='')
            else:
                record=SeqRecord(Seq(newseq).reverse_complement(),id=rename,description='')
            SeqIO.write(record,out,'fasta')
        out.close()

def sort_by_gene(contigDir):
    pathout=os.path.join(contigDir,os.path.pardir,'3.Sorted_by_gene')
    if not os.path.exists(pathout):
        os.makedirs(pathout)
    refername=os.path.join(contigDir,os.path.pardir,'reference','ref_pro.fasta')
    genelist=[]
    for record in SeqIO.parse(refername,'fasta'):
        genelist.append(record.id.split('_')[0])
    filenames=glob.glob(os.path.join(contigDir,os.path.pardir,'2.Filter_result','Hit_seqs','*.fasta'))
    for gene in genelist:
        outname=os.path.join(pathout,gene+'.fasta')
        out=open(outname,'w')
        for filename in filenames:
            judge=False
            seqdic=SeqIO.to_dict(SeqIO.parse(filename,'fasta'),lambda rec:rec.id.split('_')[0])
            if gene in seqdic.keys():
                recName=os.path.basename(filename).split('_')[0]+'_'+gene
                record=SeqRecord(seqdic[gene].seq,id=recName,description='')
                SeqIO.write(record,out,'fasta')
                judge=True
            if judge==False:
                print filename,gene,'missing gene!'
        out.close()

def separate_pep_to_single_fasta(contigDir):
    pathout=os.path.join(contigDir,os.path.pardir,'pep_for_exonerate')
    files=glob.glob(os.path.join(pathout,'ref_pro.fasta'))
    for fi in files:
        seqs=SeqIO.parse(fi,'fasta')
        for se in seqs:
            outname=os.path.join(pathout,se.id.split('_')[-1]+'_query.fasta')
            out=open(outname,'w')
            commend2='>'+str(se.id)+'\n'+str(se.seq)+'\n' 
            out.write(commend2)
            out.close()

def exonerate(contigDir):
    pathout=os.path.join(contigDir,os.path.pardir,'4.exonerate_outfile')
    if not os.path.exists(pathout):
        os.makedirs(pathout)
    filenames=glob.glob(os.path.join(contigDir,os.path.pardir,'3.Sorted_by_gene','*.fasta'))
    for filename in filenames:
        genename=os.path.basename(filename).replace('.fasta','')
        exomerate_command='exonerate --model p2g --query '\
                  +os.path.join(contigDir,os.path.pardir,'pep_for_exonerate',genename)+\
                  '_query.fasta --target '+os.path.join(contigDir,os.path.pardir,'3.Sorted_by_gene',genename)+\
				  '.fasta --showtargetgff T'+\
                  ' >'+os.path.join(pathout,genename+'_exonerate.out')+'\n'
        os.system(exomerate_command)

###get exon from exonerate def
def Query(lis):
    for i in lis:
        if 'Query: ' in i:
            gene=i.replace('Query: ','').strip().split('_')[-1]
    return gene
def Target(lis):
    for i in lis:
        if 'Target: ' in i:
            taxon=i.split(' ')[1].split('_')[0]
    return taxon
def sequence(lis):
    for i in lis:
        if 'Target: ' in i:
            sequence=i.replace('Target: ','').strip()
    return sequence
def score(lis):
    for i in lis:
        if 'Raw score:' in i:
            score=int(i.split(':')[-1].strip())
    return score
def start_end(lis):
    for i in lis:
        if 'Target range:' in i:
            start=str(i.split(':')[-1].split('->')[0].strip())
            end=str(i.split(':')[-1].split('->')[1].strip())
            exon_range=start+'\t'+end
    return exon_range
def frame(lis):
    for i in lis:
        if 'Target range' in i:
            Sta=i.split('->')[0].split(':')[-1].strip()
            End=i.split('->')[-1].strip()
            if int(End)-int(Sta)>0:
                frame=1
            else:
                frame=-1
    return frame
def takesecond(elem):
    return elem[0]

###get exon from exonerate
def get_exon_and_intron(contigDir):
    Lib=[]
    check=False
    
    path_in=os.path.join(contigDir,os.path.pardir,'4.exonerate_outfile')
    for filename in glob.glob(os.path.join(path_in,'*.out')):
        for line in open(filename):
            if 'C4 Alignment:' in line:
                check=True
                lis=[]
            if check==True:
                lis.append(line.strip())
            if 'END OF GFF DUMP' in line:
                check=False 
                Lib.append(lis)
    taxa_d={}
    for element in Lib:
        gene=Query(element)
        scor=score(element)
        taxa=Target(element)
        if taxa not in taxa_d:
            taxa_d[taxa]={}
            taxa_d[taxa][gene]=[]
        else:
            if gene not in taxa_d[taxa]:
                taxa_d[taxa][gene]=[]
        taxa_d[taxa][gene].append(scor)       
        taxa_d[taxa][gene]=sorted(taxa_d[taxa][gene])
    target_id={}
    target_seq={}
    for element in Lib:
        gene=Query(element)
        scor=score(element)
        taxa=Target(element)
        if scor == max(taxa_d[taxa][gene]):
            sequence_id=sequence(element)
            frame_id=frame(element)
            target_id[taxa+'_'+gene]=sequence_id
            query_range=start_end(element)
            target_seq[taxa+'_'+gene]={}
            target_seq[taxa+'_'+gene][query_range]=''
            for r in range(len(element)):
                if 'vulgar:' in element[r]:
                    end=int(r)
            for r in range(9, end):
                if (r-4)%5 == 0:
                    sequence_short=str(element[r+3]).split(':')[1].strip().upper()
                    if '*' in element[r+2] and '*' not in element[r]:
                        if element[r+3].split(':')[1].strip()[0] == '.':
                            black=len(element[r+3].split(':')[1].strip())-len(element[r+2].strip())
                        else:
                            black=0
                        l=element[r+2].strip().count('*')
                        new=element[r+2].strip()
                        for m in range(0,l):
                            loca=new.index('*')+black
                            sequence_short=sequence_short[:loca]+'N'+sequence_short[loca+1:]
                            lo=new.index('*')
                            new=new[:lo]+'N'+new[lo+1:]
                    if '#' in element[r+2]:
                        if element[r+3].split(':')[1].strip()[0] == '.':
                            black=len(element[r+3].split(':')[1].strip())-len(element[r+2].strip())
                        else:
                            black=0
                        h=element[r+2].strip().count('#')
                        new=element[r+2].strip()
                        for k in range(0,h):
                            loca=new.index('#')+black
                            sequence_short=sequence_short[:loca]+sequence_short[loca+1:]
                            lo=new.index('#')
                            new=new[:lo]+new[lo+1:]    
                    sequence_short=sequence_short.replace('-','').replace('{','').replace('}','')
                    target_seq[taxa+'_'+gene][query_range]+=sequence_short
            target_seq[taxa+'_'+gene][query_range]=re.sub('[a-zA-Z]{1,2}\.+[a-zA-Z]{1,2}','',target_seq[taxa+'_'+gene][query_range])
    for element in Lib:
        gene=Query(element)
        scor=score(element)
        taxa=Target(element)
        sequence_id=sequence(element)
        if len(taxa_d[taxa][gene])>1:
            if scor in taxa_d[taxa][gene][0:-1]:
                if sequence_id == target_id[taxa+'_'+gene]:
                    query_range=start_end(element)
                    target_seq[taxa+'_'+gene][query_range]=''
                    for r in range(len(element)):
                        if 'vulgar:' in element[r]:
                            end=int(r)
                    for r in range(9,end):
                        if (r-4)%5 == 0:
                            sequence_short=str(element[r+3]).split(':')[1].strip().upper()
                            if '*' in element[r+2] and '*' not in element[r]:
                                if element[r+3].split(':')[1].strip()[0] == '.':
                                    black=len(element[r+3].split(':')[1].strip())-len(element[r+2].strip())
                                else:
                                    black=0
                                l=element[r+2].strip().count('*')
                                new=element[r+2].strip()
                                for q in range(0,l):
                                    loca=new.index('*')+black
                                    sequence_short=sequence_short[:loca]+'N'+sequence_short[loca+1:]
                                    lo=new.index('*')
                                    new=new[:lo]+'N'+new[lo+1:]
                            if '#' in element[r+2]:
                                if element[r+3].split(':')[1].strip()[0] == '.':
                                    black=len(element[r+3].split(':')[1].strip())-len(element[r+2].strip())
                                else:
                                    black=0
                                h=element[r+2].strip().count('#')
                                new=element[r+2].strip()
                                for z in range(0,h):
                                    loca=new.index('#')+black
                                    sequence_short=sequence_short[:loca]+sequence_short[loca+1:]
                                    lo=new.index('#')
                                    new=new[:lo]+new[lo+1:]
                            sequence_short=sequence_short.replace('-','').replace('{','').replace('}','')
                            target_seq[taxa+'_'+gene][query_range]+=sequence_short
                    target_seq[taxa+'_'+gene][query_range]=re.sub('[a-zA-Z]{1,2}\.+[a-zA-Z]{1,2}','',target_seq[taxa+'_'+gene][query_range])
    gene_lis=[]
    path_out=os.path.join(contigDir,os.path.pardir,'4.exonerate_outfile','exon_out_fasta')
    if not os.path.exists(path_out):
        os.makedirs(path_out)    
    for key1 in target_seq:
        gene=str(key1).split('_')[-1]
        if gene not in gene_lis:
            gene_lis.append(gene)
    for gen in gene_lis:
        fp=open(os.path.join(path_out,gen+'.fasta'),'w')
        for key1 in target_seq:
            gene=str(key1).split('_')[-1]
            if gene == gen:
                if len(target_seq[key1])==1:
                    for key2 in target_seq[key1]:
                        fp.write('>'+str(key1)+'\n'+str(target_seq[key1][key2])+'\n')
                else:
                    fp.write('>'+str(key1)+'\n')
                    location_sta=[]
                    location_end=[]
                    dic_loca={}
                    for key2 in target_seq[key1]:
                        start=int(key2.split('\t')[0])
                        end=int(key2.split('\t')[1])
                        dic_loca[start]=end
                        location_sta.append(start)
                    location_sta=sorted(location_sta)
                    for ram in location_sta:
                        location_end.append(dic_loca[ram])
                    fp.write(str(target_seq[key1][str(location_sta[0])+'\t'+str(location_end[0])]))
                    for y in range(len(location_sta)):
                        if y > 0:
                            if int(location_sta[y])-int(location_end[y-1])>0:
                                fp.write(str(target_seq[key1][str(location_sta[y])+'\t'+str(location_end[y])]))
                            else:
                                chazhi=int(location_end[y])-int(location_end[y-1])
                                if chazhi > 0 :
                                    f=int(location_end[y-1])- int(location_sta[y])
                                    fp.write(str(target_seq[key1][str(location_sta[y])+'\t'+str(location_end[y])])[f:])
                    fp.write('\n')   
        fp.close()

###get intron from exonerate
    path_out_2 = os.path.join(contigDir,os.path.pardir, 'intron_test')
    path_out_1 = os.path.join(contigDir,os.path.pardir, 'intron')
    if not os.path.exists(path_out_1):
        os.mkdir(path_out_1)
    if not os.path.exists(path_out_2):
        os.mkdir(path_out_2)
        
    outs = glob.glob(os.path.join(contigDir,os.path.pardir,'4.exonerate_outfile','*.out'))
    for out in outs:
        query = os.path.split(out)[1].split('_')[0]
        fh_1 = os.path.join(path_out_1,query + '_intron.txt')
        test = os.path.join(path_out_2,query + '_test.txt')
        with open(test,'w') as fh_out:
            with open(out,'r') as fh_in:
                for lines in fh_in.readlines():
                    if 'Query:' in lines:
                        query = str(lines).replace('Query: ','').strip().split('_')[-1]
                    if 'exonerate:protein2genome:local	intron' in lines:
                        fh_out.write(str(query) + '|' + lines.replace(' ',''))
        query = os.path.split(out)[1].split('_')[0]
        with open(fh_1,'w') as fh_2:
            with open(test,'r') as fh_3:
                for lines in fh_3.readlines():
                    if query in lines:
                        line = lines.strip().split('\t')
                        target = line[0].split('|')[1]
                        start = line[3]
                        end = line[4]
                        record = target + '\t' + str(start) + '\t' + str(end)
                        fh_2.write(record + '\n')
        out_fasta=os.path.split(out)[1].replace('_exonerate.out','.fasta')
        fasta = os.path.join(contigDir,os.path.pardir,'3.Sorted_by_gene',out_fasta)
        out_intron = os.path.join(contigDir,os.path.pardir,'4.exonerate_outfile','intron_out_fasta')
        if not os.path.exists(out_intron):
            os.mkdir(out_intron)
        intron_out=os.path.split(out)[1].replace('_exonerate.out','_intron.fasta')
        fasta_exon = os.path.join(out_intron,intron_out)
        queryDict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"), key_function=lambda r: r.description)
        with open(fasta_exon,'w') as fas_out:
            with open(fh_1,'r') as f_in:
                dic_all = {}
                for lines in f_in.readlines():
                    line = lines.split('\t')
                    gene = line[0]
                    start = int(line[1])
                    end = int(line[2])
                    li = [start,end]
                    dic_all.setdefault(gene,[]).append(li)
                for key in dic_all.keys():
                    lis = sorted(dic_all[key], key=lambda x: x[0])
                    if len(lis) == 1:
                        st = lis[0][0]
                        en = lis[0][1]
                        record = '>' + key + '\n' + str(queryDict[key].seq[st-1:en])
                        fas_out.write(record +'\n')
                    else:
                        lis_2 = []
                        lis_2.append(lis[0])
                        for item in lis:
                            if lis_2[-1][1] < item[0]:
                                lis_2.append(item)
                            elif lis_2[-1][1] < item[1]:
                                lis_2[-1][1] = item[1]
                            else:
                                pass
                        seq = ''
                        for i in range(0,len(lis_2)):
                            st_2 = lis_2[i][0]
                            en_2 = lis_2[i][1]
                            seq = seq + str(queryDict[key].seq[st_2-1:en_2])
                        record = '>' + key + '\n' + str(seq)
                        fas_out.write(record + '\n')
    shutil.rmtree(path_out_1)
    shutil.rmtree(path_out_2)

def add_reference(contigDir):
    pathout=os.path.join(contigDir,os.path.pardir,'5.1.add_reference')
    if not os.path.exists(pathout):
        os.mkdir(pathout)
    refSeq=SeqIO.parse(os.path.join(contigDir,os.path.pardir,'reference','ref_nuc.fasta'),'fasta')
    ref_dic={}
    for ref in refSeq:
        ref_dic[ref.id]=ref.seq
    fasta_files=glob.glob(os.path.join(contigDir,os.path.pardir,'4.exonerate_outfile','exon_out_fasta','*.fasta'))
    for fasta_file in fasta_files:
        gene=os.path.basename(fasta_file).split('.')[0]
        out=open(os.path.join(contigDir,os.path.pardir,'5.1.add_reference',gene+'.fasta'),'w')
        handle=str('>'+'reference_'+gene+'\n'+ref_dic[gene]+'\n')
        out.write(handle)
        records=SeqIO.parse(fasta_file,'fasta')
        for record in records:
            handle=str('>'+record.id+'\n'+record.seq+'\n')
            out.write(handle)
        out.close()

def check_seq(contigDir):
    pathin = os.path.join(os.path.dirname(contigDir),'5.1.add_reference')
    pathout = os.path.join(os.path.dirname(contigDir),'5.2.for_translate_seq')

    if not os.path.exists(pathout):
        os.makedirs(pathout)
    else:
        pass
    fasta_files = os.listdir(pathin)
    handle = open(os.path.join(contigDir,os.path.pardir,'5.2.for_translate_seq','error_cds_seq.txt'),'w')
    handle.write('Error sequences: '+ '\n')
    for fasta in fasta_files:
        all_rec= []
        genename = fasta.split('.fasta')[0]
        fasta_name = os.path.join(pathin,fasta)
        for seq_record in SeqIO.parse(fasta_name, "fasta"):
            choise = 0
            end0 = int(len(seq_record.seq)/3)*3
            end1 = int((len(seq_record.seq)-1)/3)*3+1
            end2 = int((len(seq_record.seq)-2)/3)*3+2
            pro_seq0 = list(seq_record.seq[:end0].translate(stop_symbol="."))
            pro_seq1 = list(seq_record.seq[1:end1].translate(stop_symbol="."))
            pro_seq2 = list(seq_record.seq[2:end2].translate(stop_symbol="."))

            if '.' not in pro_seq0[0:-1]:
                start = 0
                end = end0
            elif '.' not in pro_seq1[0:-1]:
                start = 1
                end = end1 
            elif '.' not in pro_seq2[0:-1]:
                start = 2
                end = end2
            else:
                print'blastn error seq:'+genename+'_'+seq_record.id
                handle.write(genename+': '+seq_record.id+'\n')
                choise = 1
            if choise == 0:
                rec = SeqRecord(seq_record.seq[start:end],id=seq_record.id,description="")   
                all_rec.append(rec)
        outname = os.path.join(pathout,fasta)
        SeqIO.write(all_rec,outname,"fasta")
    handle.close()
    
def filt_gene_length_and_translate(contigDir):###
    cutoff=0.4#####
    pathout1=os.path.join(contigDir,os.path.pardir,'6.Final_seq_for_align','nuc')
    pathout2=os.path.join(contigDir,os.path.pardir,'6.Final_seq_for_align','pro')
    if not os.path.exists(pathout1):
        os.makedirs(pathout1)
    if not os.path.exists(pathout2):
        os.makedirs(pathout2)
    filenames=glob.glob(os.path.join(contigDir,os.path.pardir,'5.2.for_translate_seq','*.fasta'))
    DelLis=[]
    for filename in filenames:
        gene=os.path.basename(filename).replace('.fasta','')
        TotLen=0
        records=list(SeqIO.parse(filename,'fasta'))
        taxa=len(records)
        if taxa>1:
            out1=open(os.path.join(pathout1,gene+'_nuc.fasta'),'w')
            out2=open(os.path.join(pathout2,gene+'_pro.fasta'),'w')
            for record in records:
                if record.id.split('_')[0]!='reference':
                    TotLen+=len(record.seq)
            MeanLen=float(TotLen)/(taxa-1)
            for record in SeqIO.parse(filename,'fasta'):
                if record.id.split('_')[0]!='reference':##
                    proRec=SeqRecord(record.seq.translate(),id=record.id.split('_')[0],description='')####
                    nucRec=SeqRecord(record.seq,id=record.id.split('_')[0],description='')####
                    length=len(record.seq)
                    Cover=length/MeanLen
                    if Cover >= cutoff:
                        SeqIO.write(nucRec,out1,'fasta')
                        SeqIO.write(proRec,out2,'fasta')
                    else:
                        candidate=gene+'|'+record.id+'|'+str(Cover)
                        DelLis.append(candidate)
                else:#####
                    proRec=SeqRecord(record.seq.translate(),id=record.id.split('_')[0],description='')
                    nucRec=SeqRecord(record.seq,id=record.id.split('_')[0],description='')####
                    SeqIO.write(nucRec,out1,'fasta')
                    SeqIO.write(proRec,out2,'fasta')
        out1.close()
        out2.close()
        
    files=glob.glob(os.path.join(pathout1,'*.fasta'))
    for file1 in files:
        recos=list(SeqIO.parse(file1,'fasta'))
        taxa1=len(recos)
##        print file1,taxa1
        if taxa1<2:
            os.remove(file1)
    fils=glob.glob(os.path.join(pathout2,'*.fasta'))
    for file2 in fils:
        res=list(SeqIO.parse(file2,'fasta'))
        taxa2=len(res)
##        print file2,taxa2
        if taxa2<2:
            os.remove(file2)
    print 'Filtered sequence:',DelLis,len(DelLis)
    

def run(contigDir):
    time1=time.time()
    
    separate_pep_to_single_fasta(contigDir)
    print "Data analyse for amplicon capture data..."
    
    print "Step:1 tblastn to search target gene from contigs"
    print "makeblastdb, result in 1.Blast_result/blast_db_file/*"
    makeblastdb(contigDir)
    print "tblastn, result in 1.Blast_result/*.xml"
    tblastn(contigDir)
    print "Step:1 tblastn to search target gene from contigs complete!\n"

    print "Step:2 get best hit from blast result"
    print 'blast result xml to txt, result in 1.Blast_result/*.txt'
    blast_result_xml_to_txt(contigDir)
    print 'get best Hit info, result in 2.Filter_result/*.txt'
    get_best_Hit_info(contigDir)
    print 'get hit seqs from best Hit info, result in 2.Filter_result/Hit_seqs/*.fasta'
    get_hit_seqs(contigDir)
    print "Step:2 get best hit from blast result complete!\n"
    
    print "Step:3 sort sequences by gene, result in 3.Sorted_by_gene/*.fasta"
    sort_by_gene(contigDir)
    print "Step:3 sort sequences by gene complete!\n"
    
    print "Step:4 run exonerate to separate exon and intron"
    print "exonerate, result in 4.exonerate_outfile/*.out"
    exonerate(contigDir)
    print "separate exon and intron, result in 4.exonerate_outfile/exon_out_fasta/*.fasta & "+\
          "4.exonerate_outfile/intron_out_fasta/*.fasta"
    get_exon_and_intron(contigDir)
    print "Step:4 run exonerate to separate exon and intron complete!\n"
    
    print "Step:5 add reference,result in 5.add_reference"
    add_reference(contigDir)
    print "Step:5 add reference complete!\n"
    check_seq(contigDir)
    print"check seq complete!"
    
    print "Step:6 filter gene length and translate, result in 6.Final_seq_for_align/nuc/*.fasta & "+\
          "6.Final_seq_for_align/pro/*.fasta"
    filt_gene_length_and_translate(contigDir)
    print "Step:6 filter gene length and translate complete!\n"
    

    time2=time.time()
    print 'Analyse time: '+str(time2-time1)
    print "The work of extrcting orthologous sequence groups from assembled contigs is complete!!!"
	
if __name__ == '__main__':
    dirpath_o = os.path.abspath(os.path.join(os.path.dirname('__file__'), os.path.pardir))
    os.chdir(dirpath_o)
    parser = argparse.ArgumentParser(description='python Extract_target_coding_sequences_from_amplicon_capture_data.py'
	' -c contigs directory')
    parser.add_argument('--contigDir', '-c', help='contigs directory', required=True)
    args = parser.parse_args()
    run(args.contigDir)



