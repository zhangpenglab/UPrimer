<p align="center">
<img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/7a09be4f-3b89-47a8-a12b-281bc77ccb5a" alt="Drawing" width="200"/>
</p>


UPrimer: A Clade-Specific Primer Design Program Based on Nested PCR Strategy and Its Applications in Amplicon Capture Phylogenomics
=================
Jiaxuan Li, Guangcheng Han, Xiao Tian, Dan Liang*, Peng Zhang*<br>
*State Key Laboratory of Biocontrol, School of Life Sciences, Sun Yat-Sen University, Guangzhou, China*


## Citation
Jiaxuan Li, Guangcheng Han, Xiao Tian, Dan Liang*, Peng Zhang*. 2023. UPrimer: A Clade-Specific Primer Design Program Based on Nested PCR Strategy and Its Applications in Amplicon Capture Phylogenomics. Molecular Biology and Evolution. DOI: 10.1093/molbev/msad230

<div class="page-break"></div>


# Table of content

1. [Installation and Requirements](#installation-and-requirements)<br>
3. [Design architecture for UPrimer](#design-architecture-for-uprimer)<br>
 * [Module I: Obtain candidate multiple sequence alignments based on genome data for the target taxon](#module-i-obtain-candidate-multiple-sequence-alignments-based-on-genome-data-for-the-target-taxon)<br>
 * [Module II: Design nested-PCR primers of NPCLs](#module-ii-design-nested-pcr-primers-of-npcls)<br>
   [##Scoring algorithms for primers](#scoring-algorithms-for-primers)<br>
3. [Preparation before running UPrimer](#preparation-before-running-uprimer)<br>
 * [Scan available genome data for the target taxon](#scan-available-data)<br>
 * [Select suitable data for developing NPCL primers (highly important)](#select-suitable-data)<br>
 * [Download data](#download-data)<br>
4. [Usage](#usage)<br>
5. [Input](#input)<br>
6. [Output](#output)<br>
7. [A step-by-step instruction for using UPrimer to develop NPCL primers-taking Lepidoptera as an example](#a-step-by-step-instruction-for-using-uprimer-to-develop-npcl-primers---taking-lepidoptera-as-an-example)
   
  [##A part of analyzing amplicon capture data: Extract orthologous sequence groups from assembled contigs](#extract-orthologous-sequence-groups-from-assembled-contigs)


<div class="page-break"></div>

# Installation and Requirements

**As UPrimer is a Python script, it does not require installation. However, UPrimer relies on several other software and packages that need to be installed on your system and available in your PATH：**

- Python == 2.7.15
- Biopython == 1.74
- Pandas == 0.24.2
- Numpy == 1.16.6
- BLAST+ version 2.12.0 
- MEGA-CC version 10.18 (https://github.com/zhangpenglab/UPrimer/tree/main/Accessory/megacc_10.1.8_amd64.tar.gz)
- PAL2NAL version 14
- Exonerate version 2.4.0 #This software will be used for subsequent capture data analysis.


## A step-by-step installation guide (For Linux)

Before installing, you should to download and install the software Miniconda (https://docs.conda.io/en/latest/miniconda.html) on your system. After that, you can follow the steps below:

(1) Download the UPrimer package.
~~~
git clone https://github.com/zhangpenglab/UPrimer.git 
~~~

(2) Install python, biopython, BLAST+, pandas, pal2nal and exonerate through conda.
~~~
unzip /path/to/UPrimer-main.zip
conda env create -f /path/to/UPrimer-main/UPrimer-conda-env.yml
~~~

(3) Install megacc and put it on your PATH.

~~~
cd /path/to/UPrimer-main/Accessory/
tar xvfz megacc_10.1.8_amd64.tar.gz

vim ~/.bashrc                                           
export PATH=$PATH:/path/to/UPrimer-main/Accessory
source ~/.bashrc
~~~

(4) After installation, check if the software and packages required for UPrimer are functioning properly.
~~~
bash
conda activate UPrimer
~~~

~~~
megacc
------------------
MEGA-CC has logged the following error:
  When         = Wednesday, June 28, 2023 PM03:32:28 PM
  Data file    =   When         
  AnalysisFile = 
  Message      = The analysis (MAO) file you specified does not exist.  Please check your spelling and try again.

  Please see the summary file for warnings/messages
  Please validate all your inputs and try again. If you think this is a
  bug, please report it at http://www.megasoftware.net

~~~

~~~
pal2nal.pl
------------------
pal2nal.pl  (v14)

Usage:  pal2nal.pl  pep.aln  nuc.fasta  [nuc.fasta...]  [options]
...
~~~

~~~
blastn -version
------------------
blastn: 2.12.0+
~~~

~~~
python -V
------------------
Python 2.7.15
~~~

~~~
pip list
------------------

biopython  1.74
numpy      1.16.6
pandas     0.24.2
...
~~~~

## FAQ
If you encounter difficulties installing UPrimer, it is possible that the included packages are incompatible with your system. Here are some clues to help you troubleshoot the problems:

1. Please make sure that **Python 2.7.15** and properly installed and are in your PATH. Type ```python``` to check. 
2. UPrimer is **NOT compatible to Python 3+**. We are sorry for the inconvenience.


If you have any further questions, please contact us.

---

# Design architecture for UPrimer
## Module I Obtain candidate multiple sequence alignments based on genome data for the target taxon

The whole workflow of UPrimer comprises two main modules. **The first module (Make_MSAs_suitable_for_primer_design_Part_I.py) aims to obtain candidate MSAs based on the genome data of the target taxon.** The module contains five main steps, and each step corresponds to a script:
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/8358bad4-f5cf-4069-8d92-13e7801dbe4e" width="770" height="600"/>
</div>
<br /><br />

- Step 1 —— **Identify long and single-copy exons from the genome of a reference species** (Screen_exons_1.py)：The reference species can be any species of the target taxon, but should have well-annotated genome data available. **The input data of this step are exome, proteome, and genome sequences of the reference species.** UPrimer first uses *BLASTX* to trim each exon in the exome to the correct translation frame using the proteome as a guide. Subsequently, it discards exons shorter than a predefined value (default: 300 bp). The program then uses *BLASTN* to search the remaining exons against the genome to remove exons that are not single-copy. The criterion is as followed: if an exon has a second BLAST hit with similarity > 50% and coverage > 30%, this exon is considered to have a similar copy in the genome and is not single-copy.
<br /><br />

- Step 2 —— **Obtain orthologous sequences of the exons of the reference species from ingroup and outgroup species** (Append_orthologous_sequences_2.py)：**The ingroup species (required)** belong to the same target taxon as the reference species, and it is better to use more ingroup species to cover the whole phylogenetic span of the target taxon. **Outgroup species (optional)** do not belong to the target taxon, and having more outgroup species in the analysis can ensure finding conserved regions for primer design. The input data of this step are **genome sequences of the ingroup species** and **coding sequences (CDS or transcriptome) of outgroup species**. Based on the exons of the reference species, UPrimer employs a mutual best-hit blast strategy *(MBH BLAST)* to extract orthologous sequences from the genomes or coding sequences of both ingroup and outgroup species. Among the identified orthologous sequences, only those with a length greater than 300 bp and without stop codons are retained. For each exon of the reference species, the program combines all its filtered orthologous sequences from both ingroup and outgroup species, constructing orthologous sequence groups (OGs) at both the DNA and protein levels.
<br /><br />

- Step 3 —— **Construct MSAs for each orthologous sequence group** (Make_codon_alignments_3.py)： UPrimer first aligns the OGs' protein sequences using *MEGA-CC*. Subsequently, *PAL2NAL* is utilized to construct codon alignments by incorporating the DNA sequences of the OGs and the resulting protein alignments. The program then trims the protein and DNA alignments on both ends, guided by the exons of the reference species. 
<br /><br />

- Step 4 —— **Remove problematic sequences and trim alignments to increase the quality of MSAs** (Filter_alignments_4.py)：UPrimer removes sequences from each multiple sequence alignment (MSA) if they exhibit high levels of missing data (> 60%) or a high gap ratio (> 60%). Furthermore, if a sequence's average similarity to all other sequences within the alignments falls below 30%, it will also be discarded. Although the selection of the 30% cut-off value is somewhat subjective, it is a reasonable threshold to eliminate problematic sequences in a MSA resulting from incorrect orthology assignment or sequence errors.
<br /><br />

- Step 5 —— **Pick out suitable MSAs for subsequent primer design** (Find_blocks_5.py):  After refining the alignments, UPrimer proceeds to select suitable candidate alignments for primer design. It first searches each alignment from both the left and right ends to ensure the presence of two conserved primer blocks, each consisting of eight amino acids and exhibiting a sequence similarity greater than 50%. Subsequently, the alignment is trimmed by removing the regions outside the leftmost and rightmost primer blocks. The trimmed alignments must have a length exceeding 300 bp; otherwise, they are discarded. Additionally, UPrimer will discard highly conserved alignments (with a similarity greater than 90%) that contain too few informative sites.
<br /><br />  

## Module II Design nested-PCR primers of NPCLs

**The second module (Design_universal_primer_sets_Part_II.py) aims to design universal nested-PCR primer sets of NPCLs based on candidate MSAs**. The workflow of this module is as follows:
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/7b904c8e-81e8-47cf-8a47-c28d163a2caa" alt="Drawing" />
</div>
<br /><br />

- Step 1 —— **Search for primer blocks**: For each alignment, UPrimer searches for all conserved primer blocks that are 7 or 8 amino acids in length. The sequence similarity between the conserved primer block within the ingroup and the reference species is required to be at least 85%. The sequence similarity across all species should be at least 75%.
<br /><br />

- Step 2 —— **Design forward and reverse primers**: Firstly, UPrimer infers the consensus amino acid for each column within the ingroup species and the reference species using the amino acid primer block. Then, nucleotide alignments are used to obtain the consensus nucleotide sequence. Finally, based on this consensus nucleotide sequence, UPrimer designs both forward and reverse primers.
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/23957157-9f5f-47be-8467-0d40afb800c4" alt="Drawing" width="500" height="300"/>
</div>
<br /><br />

- Step 3 —— **Match and filter primer pairs**: UPrimer matches all forward and reverse primers to list all possible primer pairs and filters them by primer degeneracy (the overall degeneracy of the nucleotide primers is less than or equal to 8192) and amplification length (300~2100 bp).
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/ff37784e-12c4-4be7-98f8-9c595d7f3dd9" alt="Drawing" width="500" height="500"/>
</div>
<br /><br />

- Step 4 —— **Search for outer forward and reverse primers**: For every retained primer pair, UPrimer searches for their outer forward and reverse primers (if they exist) within a flanking region of 450 bp.
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/9e393997-9874-48ae-8bdb-94fcbd213ae9" alt="Drawing" width="500" height="300"/>
</div>
<br /><br />


- Step 5 —— **Obtain a list of all possible nested-PCR primer pairs**
<br /><br />

- Step 6 —— **Score all nested-PCR primer pairs based on their PCR performances and phylogenetic information**: UPrimer scores each nested-PCR primer pair based on its potential PCR performance (named “ScorePCR”), taking into account the conservativeness of primer blocks, primer degeneracy, and primer complexity, as well as the phylogenetic informativeness of the locus (named “ScoreINFOR”), which depends on the variability of the amplification regions (*detailed calculation algorithm see [Scoring algorithms for primers](#scoring-algorithms-for-primers)*). 

- Step 7 —— **Calculate a total score for each primer pair and sort them**: UPrimer then calculates a total score for each primer pair, based on a weighting parameter “PIs” (default = 1) of PCR performance score and phylogenetic information score. The formula for calculation is as follows: 
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/d64734b9-890c-427e-bb0e-e42339732045" alt="Drawing" width="450" height="60"/>
</div>
<br /><br />

- Step 8 —— **Output NPCL primer development results**: UPrimer selects *the highest-scoring* nested-PCR primer from the list of primers for each alignment. Then, UPrimer outputs ***a nested-PCR primer table***, as well as ***the reference nucleotide and peptide sequences*** of the NPCL regions, which will be used for subsequent capture data analysis.
<br /><br />


### Scoring algorithms for primers
**@PCR amplification performance: ScorePCR = (PrimerCON + PrimerDEG + PrimerCOM)/3, maximum score: 100**

(1) The conservation score of amino acid primer blocks (PrimerCON, maximum score: 100)

The conservation score of amino acid primer blocks is determined by assessing **the similarity of these primer blocks among all ingroup species and the reference species**. A higher similarity corresponds to a higher score, with a maximum score of 100, indicating 100% similarity for that specific primer block. UPrimer calculates the conservation score for each of the four primer blocks and then computes the average score to obtain the overall 'PrimerCON' score for a pair of nested PCR primers.

##For example, calculating the conservation score of a 8 amino acid primer block:
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/26394d84-5bc4-404c-b678-2fe65fb894a5" alt="Drawing" width="450" height="400"/>
</div>  
<br /><br />

(2) The degeneracy score of amino acid primer sequences (PrimerDEG, maximum score: 100)

Higher degeneracy of amino acid primer sequences can potentially decrease the binding efficiency between the primer and the template, as the wider range of choices for degenerate bases may lead to unstable binding. Therefore, selecting primers with lower degeneracy can improve binding efficiency and amplification specificity. Considering the 5' end of the amino acid primer sequence has a smaller impact on amplification efficiency compared to the 3' end, UPrimer assigns weight ratios to each amino acid position based on its distance from the primer's 3' end. The amino acid site at the 3' end is assigned the highest weight ratio of 1, the amino acid site at the second position from the 3' end has a weight ratio of 0.9, and the amino acid site at the third position from the 3' end has a weight ratio of 0.8. The weight ratios decrease progressively until the last amino acid site. The penalty score for each amino acid site is calculated by multiplying its degeneracy by the weight ratio. The penalty score for each primer is obtained by subtracting the sum of penalty scores for all amino acid sites from 100. UPrimer calculates the degeneracy score for each of the four amino acid primers separately, and then calculates the average score to obtain the ‘PrimerDEG’ score for a pair of nested PCR primers.

##For example, calculating the degeneracy score of a 8 amino acid primer sequence:
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/aea4c39d-b6ec-49f4-a9c6-cb38076413a5" alt="Drawing" width="610" height="300"/>
</div>
<br /><br />

(3) The complexity score of amino acid primer sequences (PrimerCOM, maximum score: 100)

The complexity of amino acid primer sequences is determined by evaluating the amino acid composition of the primer sequences. UPrimer employs a counting method to tally the occurrence of each type of amino acid within the primer sequence. UPrimer first applies the formula *N × (N-1)* to calculate the penalty score for each type of amino acid, where N represents the number of occurrences of that particular amino acid. The complexity score for each primer is obtained by subtracting the sum of penalty scores for all types of amino acid from 100. UPrimer calculates the complexity score for each of the four amino acid primers separately, and then calculates the average score to obtain the ‘PrimerCOM’ score for a pair of nested PCR primers.

##For example, calculating the complexity score of a 8 amino acid primer sequence:
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/1c47e223-334b-463e-b981-d31e74d9d88e" alt="Drawing" width="500" height="300"/>
</div>
<br /><br />


**@Phylogenetic informativeness: ScoreINFOR, maximum score: 100**

The phylogenetic information of the amplification region is determined by the sequence variability. To convert this value into a percentage scale, UPrimer utilizes the following formula: 
ScoreINFOR = (Infor - InforMin)/(InForMax-InForMin) * (Max_Infor_score - Min_Infor_score) + Min_Infor_score.


Where:
- **Infor** represents the phylogenetic information value of a specific amplification region, which is calculated by summing the sequence diveregence in each column of all ingroup species and the reference species. The formula for calculation is as follows: InforSelf = D1 + D2 + D3 + ... + Dn. Here, n represents the length of the amplification region (i.e., the number of column positions). Dn represents the amino acid divergence of all ingroup species and the reference species in each column position.
- **InforMin** represents the **lowest value** of phylogenetic information within different primer pairs' amplification regions in the same MSA.
- **InforMax** represents the **highest value** of phylogenetic information within different primer pairs' amplification regions in the same MSA.
- **Max_Infor_score** represents the maximum score limit, which is 100 points.
- **Min_Infor_score** represents the minimum score limit, which is 50 points.

##For example, calualting the phylogenetic information score of all nested primer pairs in the same MSA.
<div>
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/264461d6-eba6-4149-af76-c413eda6b2d4" alt="Drawing" width="1000" height="440"/>
</div>
<br /><br />


# Preparation before running UPrimer

## Scan available data
What available genomic resources can be used for the target group? Such information can be obtained from the genome list in the NCBI Genome database (https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/).

It is worth noting that there is no taxonomy information in the downloaded genome list. To do this, the user need first collect the species name in the downloaded genome list and use these information to extract their corresponding NCBI taxonomy information by using **TaxonKit** (https://github.com/shenwei356/taxonkit) and its related package **taxdump** (https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz). After successfully installing, run:
<br /><br /> 

(1) Extract **Taxonomy ID** based on **Species name**:
~~~
taxonkit name2taxid Genome_list_only_contain_species_name.txt > Genome_list_Name_TaxonomyID.txt
~~~
For example:
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/bc09a565-1ffb-4f05-9663-5bb8259ceac8" alt="Drawing" width="650" height="300"/>
</div>
<br /><br />


(2) Extract **Taxonomy information** based on **Taxonomy ID**:
~~~
taxonkit lineage Genome_list_Name_TaxonomyID.txt -i 2 > Genome_list_Name_TaxonomyID_Lineage.txt
~~~
For example:
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/88c6eafc-5fa1-4872-855d-7f7ce1bed7f9" alt="Drawing" width="600" height="300"/>
</div>
<br /><br />

(3) Collate lineage information of the target group as the user wants:
~~~
taxonkit reformat Genome_list_Name_TaxonomyID_Lineage.txt -i 3 > Genome_list_Name_TaxonomyID_Lineage_Collated.txt
~~~

(4) Finally, the user needs to obtain a genome list like this:
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/75a41f93-78ba-4c6a-86d4-929c3f830aae" alt="Drawing" width="800" height="250"/>
</div>
<br /><br />


## Select suitable data
- To successfully develop universal primers of NPCLs for a target taxonomic group by using UPrimer, the **input genomic data** plays a crucial role. The input genomic data of UPrimer included **the exome, proteome and genome sequences of the reference**, **the genome sequences of ingroup species**, and **the coding sequences (CDS) of outgroup species**. The recommended suggestions for the input genomic data are as follows (carefully reading): 

    - **REFERENCE**: The selection of an appropriate reference species is crucial for obtaining more NPCL primers. The reference species should represent the targeted clade and possess a **high-quality genome and annotations** to facilitate the identification of a greater number of long and single-copy exonic loci. The criteria for selecting the reference species can include indicators such as **Genome assembly level (priority level: chromosome > scaffold > contig), Scaffold N50 (Higher better), BUSCO complete ratio (Higher better)**. If these indicators are similar among different candidate reference species, **the number of exons** will be compared to determine the most suitable reference species (higher better). All these information can be found at the NCBI genome database (https://www.ncbi.nlm.nih.gov/data-hub/genome/);
    
    - **INGROUP**：There should be adequate ingroup species representation. To ensure that the designed NPCL primers are universal across the targeted clade, **each major linage** of the target clade should have at least one representative ingroup species included in the analysis and the total number of ingroups is **relatively moderate (10-15)**. This will help minimize the bias in locating conserved primer blocks and increase the accuracy of NPCL primer design;
    
    - **OUTGROUP**: The inclusion of outgroups can enhance the credibility of discovering conserved primer blocks and increase the PCR success rate of the designed primers. However, improper utilization of outgroups may impact the number of the ultimate NPCL outputs. **When no closely related taxa are available as outgroups** for the target clade, introducing distantly related ones may pose difficulties in identifying orthologous sequences, obtaining a sufficient number of candidate MSAs, and identifying conserved primer blocks. When developing universal primers, it is an optional step to analyze the outgroup. The user need to consider this based on the taxonomic level of the target group. If developing universal primer sets **at lower taxonomic levels, such as class or order**, we suggest using closely related sister taxa to the target group as outgroups. When developing universal primers **at higher taxonomic levels, such as phylum or subphylum**, we no longer recommend using outgroups, because ingroups already cover most of the major clades (as required by developed), and the role of outgroups becomes less important.

## Download data
If the genomic resources of the target group meet the requirements of UPrimer, we strongly recommend using **the ncbi-genome-download tool** (https://github.com/kblin/ncbi-genome-download) to download all used genome data (**FASTA-formatted**), because this tool allows for efficient bulk downloading of genome data at a high speed.


For example, if we want to download genome data for the order Lepidoptera:

(1) Download genome, proteome, and genome annotation GB files based on the accession numbers of the reference species *Bombyx mori* (GCF_014905235.1):

~~~
ncbi-genome-download all --section refseq --formats fasta -A GCF_014905235.1 --flat-output -o /path/to/Reference -r 100
ncbi-genome-download all --section refseq --formats protein-fasta -A GCF_014905235.1 --flat-output -o /path/to/Reference -r 100
ncbi-genome-download all --section refseq --formats genbank -A GCF_014905235.1 --flat-output -o /path/to/Reference -r 100
~~~

(2) Extract exome data from the genome annotation GenBank file for the referecne sepecies *Bombyx mori*:
<br /><br />
**CAUSION**: The exome of the reference species cannot be directly downloaded through the ncbi-genome-download tool. To resolve this issue, we provide a custom script named **"Extract_exome_from_genbank_file.py"** (https://github.com/zhangpenglab/UPrimer/tree/main/Accessory) that can extract exome data from the annotation 'gbff' file. Commands as below: 

~~~
cp /path/to/UPrimer-main/Accessory/Extract_exome_from_genbank_file.py /path/to/Reference
cd /path/to/Reference
mkdir gbff
cp GCF_014905235.1_Bmori_2016v1.0_genomic.gbff ./gbff
python Extract_exome_from_genbank_file.py
~~~

The extracted exome FASTA file of the reference species *Bombyx mori* is like this:
  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/cdeb5083-9eb6-47e3-9539-4eab30926356" alt="Drawing" width="780" height="400"/>
  </div>    
<br /><br />


(3) Download genome sequence data based on the accession numbers of the ingroup species *Plutella xylostella* (GCA_932276165.1) and *Aricia agestis* (GCA_905147365.1):

~~~
ncbi-genome-download all --section refseq --formats fasta -A GCA_932276165.1 --flat-output -o /path/to/Ingroups -r 100
ncbi-genome-download all --section refseq --formats fasta -A GCA_905147365.1 --flat-output -o /path/to/Ingroups -r 100
~~~


(4) Download CDS data based on the accession numbers of the related outgroup species *Bombus terrestris* (Hymenoptera, Apidae) (GCF_014905235.1) and *Tribolium castaneum* (Coleoptera, Tenebrionidae) (GCF_000002335.3):
~~~
ncbi-genome-download all --section genbank --formats cds-fasta -A GCF_014905235.1 --flat-output -o /path/to/Outgroups -r 100
ncbi-genome-download all --section genbank --formats cds-fasta -A GCF_000002335.3 --flat-output -o /path/to/Outgroups -r 100
...
~~~

***###NOTE：After completing the download of genome data through the NCBI genome database, the user only needs to unzip all files and rename them according to UPrimer's requirements (see [Input](#input)). There is no need to modify the sequence ID or description information. Please pay attention to this point.***

---

# Usage
Here is the help message you get by running the script without parameter:

~~~
Usage :
python UPrimer.py --help

Required parameters：
  --file SRT :Input a text file containing a list of species names for the reference, ingroup, and outgroup species (filename: 'species_list.txt') [short: -F]
  --ConsiderOutgroups 0|1 :Whether to consider outgroup species when designing NPCL primers (0: Do not consider. 1: Consider) [short: -CO]
  --TaxaGroup SRT :Input the name of the target taxon for NPCL primer development (for example: Lepidoptera, Vertebrata, Bivalvia) [short: -TG]  

Optioinal parameters：
  --len INT : The length cutoff for exons of the reference species (DEFAULT: 300) [short: -L]
  --thd INT : The number of threads to use when running BLAST (DEFAULT: 10) [short: -T]
  --cov FLOAT : Single-copy exon selection criteria: The length coverage between exons and their second-hit sequences (DEFAULT: 0.3) [short: -C]
  --sim FLOAT : Single-copy exon selection criteria: The similarity between exons and their second-hit sequences (DEFAULT: 0.5) [short: -S]
  --SelectScale INT : The value of a criteria used in removing rogue sequences, which is typically larger and more strict. This parameter usually does not need to be adjusted (DEFAULT: 1.1)[short: -SS]
  --RemoveRigor INT : The value of a sequence identity cutoff used in removing rogue sequences from MSAs. Set this parameter based on the quality of genome data. A larger value (e.g., 50 or 55) is suitable for high-quality genome data, while a lower value (e.g., 30 or 35) is suitable for low-quality genome data (DEFAULT: 55) [short: -RR]
  --exp FLOAT : The value of the expected identity of all species in 7AA/8AA primer blocks (DEFAULT: 0.5) [short: -E]
  --outgroupNum INT : The minimum number of outgroup species (DEFAULT: 1) [short: -ON]
  --ingroupNum INT : The minimum number of ingroup species (DEFAULT: 1) [short: -IN]
  --totaltaxanum INT : The minimum number of all species (DEFAULT: 3) [short: -TN]
  --tril INT : The permissible length cutoff of amino acid during trimming sequences. This parameter usually does not need to be adjusted (DEFAULT: 10) [short: -l]
  --trit FLOAT : The ratio of gap length to AA length during trimming sequences. This parameter usually does not need to be adjusted (DEFAULT: 1.5) [short: -t]
  --con FLOAT : The overall conservation/identity of MSAs. If the overall similarity of one MSA is larger than this value, UPrimer will not design primers for it because the MSA is too conserved (DEFAULT: 0.9) [short: -C]
  --PIs INT : The weighting ratio of PCR performance score and phylogenetic information score when scoring all possible nested-PCR primers within the same MSA. Increasing this value will reduce the length of the amplified regions. Please set this parameter with caution (DEFAULT: 1) [short: -PIs]
  --fwd INT : The length of amino acid between the outer primer F1 and the inner primer F2. This parameter usually does not need to be adjusted (DEFAULT: 150) [short: -W]
  --rvs INT : The length of amino acid between the outer primer R1 and the inner primer R2. This parameter usually does not need to be adjusted (DEFAULT: 150) [short: -R]
  --max INT : The maximum length of target NPCLs (amino acid length, from F2 to R2) (DEFAULT: 700) [short: -M]
  --min INT : The minimum length of target NPCLs (amino acid length, from F2 to R2) (DEFAULT: 100) [short: -N]
  --all INT : The upper threshold for the degeneracy of the entire nucleotide primer sequences. This parameter usually does not need to be adjusted (DEFAULT: 8192) [short: -A]
  --inthred2 FLOAT : The mean identity of the reference species and all ingroups in 7AA/8AA inner (F2 and R2) primer blocks (DEFAULT: 0.85) [short: -IT2]
  --allthred2 FLOAT : The mean identity of all species in 7AA/8AA inner (F2 and R2) primer blocks (DEFAULT: 0.75) [short: -AT2]
  --inthred1 FLOAT : The mean identity of the reference species and all ingroups in 7AA/8AA outer (F1 and R1) primer blocks (DEFAULT: 0.80) [short: -IT1]
  --allthred1 FLOAT : The mean identity of all species in 7AA/8AA outer (F1 and R1) primer blocks (DEFAULT: 0.70) [short: -AT1]
  --PId INT : The weighting ratio of PCR performance score and phylogenetic information score when scoring highest-scoring nested-PCR primers derived from different MSAs (DEFAULT: 1) [short: -PId]

Here are some examples of how to develop nested-PCR primers for NPCLs using UPrimer:

(1) Developing nested-PCR primers for NPCLs for the target taxon at an order or class level:
python UPrimer.py -F species_list.txt -CO 1 -TG TargetGroupName -IT2 0.85 -AT2 0.75 -IT1 0.80 -AT1 0.70 -PId 3

(2) Developing nested-PCR primers for NPCLs for the target taxon at a subphylum or phylum level:
python UPrimer.py -F species_list.txt -CO 0 -TG TargetGroupName -IT2 0.90 -AT2 0.90 -IT1 0.90 -AT1 0.90 -PId 5

(3) Developing nested-PCR primers for NPCLs for the target taxon with relatively low-quality genome data:
python UPrimer.py -F species_list.txt -CO 1 -TG TargetGroupName -RR 30 -IT2 0.85 -AT2 0.75 -IT1 0.80 -AT1 0.70 -PId 3

(4) Only running the module of obtaining candidate MSAs from genome data of the target taxon:
python Make_MSAs_suitable_for_primer_design_Part_I.py -F species_list.txt 

(5) Only running the module of designing NPCL primers based on candidate MSAs:
python Design_universal_primer_sets_Part_II.py -F species_list.txt -D 4.Candidate_peptide_iden0.5_MSAs_for_primer_design -TG TargetGroupName -IT2 0.85 -AT2 0.75 -IT1 0.80 -AT1 0.70 -PId 3
~~~
<br /><br />


# Input
The user should begin by creating a new folder (e.g., XXX_universal_primer_development) and placing the UPrimer program package inside it. Additionally, the user need to prepare a text file and three data folders to store the respective data files. The details are as follows:


(1) A text file (***"species_list.txt"***): This file should contain a list of species names, including the names of the reference species, ingroup species and outgroup species.
  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/0a02bcaa-1c27-4fa8-b257-4421036ebdcb" alt="Drawing" width="480" height="400"/>
  </div>
<br /><br />


(2) The folder ***"Reference"***: This folder contains exonome data (reference_name_exon.fasta), genome data (reference_name_genome.fasta), and proteome data (reference_name_pep.fasta) of the reference species. The FASTA file format is illustrated in the figure below.
  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/dcc011e8-277e-4b48-a85b-4f85ee4e5c1e" alt="Drawing" width="600" height="500"/>
  </div>    
<br /><br />

(3) The folder ***"Ingroups"***: This folder contains the genome data of all ingroup species. The naming convention for these files is as follows: species_name + _genome.fasta. The FASTA file format is illustrated in the figure below.
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/efd58e5d-12c3-4873-afeb-b36e7664b4a1" alt="Drawing" width="460" height="570"/>
  </div>
<br /><br />    

(4) The folder ***"Outgroups"***: This folder contains the CDS data of all outgroup species. The naming convention for these files is as follows: species_name + _cds.fasta. The FASTA file format is illustrated in the figure below.
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/e1ab4ec2-d13d-4673-991b-30dae87c9b9c" alt="Drawing" width="800" height="600"/>
  </div>
<br /><br />



# Output
(1) The candidate nucleotide and peptide MSAs used for primer design can be found in the folder "4.Candidate_nucleotide_iden0.5_MSAs_for_primer_design"*** and ***"4.Candidate_peptide_iden0.5_MSAs_for_primer_design".
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/9c020d28-ee83-44d3-865a-0e83685a72ea" alt="Drawing" width="900" height="310"/>
  </div>
<br /><br />

(2) The final designed NPCL primer sets can be found in the folders "6.Designed_nested-PCR_primer_set_of_NPCLs_PIratio_xxx". In this folder, the user could find three tables:
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/3648085b-3d06-4c6f-a650-ce373c7ee849" alt="Drawing" width="750" height="320"/>
  </div>
<br /><br />

- **Highest-scoring_nested-PCR_primer_set_PIdXX_numXXX_mlenXXX.xls**: This table presents detailed information for the highest-scoring nested primers obtained from various MSAs, including such as NPCL ID, Target region length, primer sequence, primer position, primer block identity. It is important to note that the highest-scoring primers for different MSAs have already be sorted based on their total scores (from high to low). Researchers should prioritize using the primer pairs that are ranked higher on the list.

  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/f40be05c-fdea-4ce0-ba55-6ae69a74b81d" alt="Drawing" width="950" height="360"/>
  </div>
  <br /><br />

- **Synthezised_highest-scoring_nested-PCR_primer_PIdXXX_numXXX_mlenXXX.xls**: This table is a simplified version of the "Highest-scoring_nested-PCR_primer_set_PIdXX_numXXX_mlenXXX.xls" table, containing only the NPCL ID, nested PCR primer sequences, and primer lengths. It can be directly used for primer synthesis.
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/1b8d525a-69d0-466c-8f42-d4a26fa0a8e3" alt="Drawing" width="650" height="550"/>
  </div>
  <br /><br />

- **Candidate_nested-PCR_primer_set_PIdXXX_numXXX_mlenXXX.xls**: If PCR amplification using the highest-scoring primer pair fails, the user can explore other candidate primers listed in the table. Each MSA provides ten candidate nested PCR primer pairs. These candidate primers are sorted based on their total scores, and the author can try them sequentially, one by one, to identify a successful amplification.
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/62d3b994-e110-4d30-b766-8e4e05fb57a4" alt="Drawing" width="920" height="350"/>
  </div>
  <br /><br />

(3) The reference NPCL nucletide and peptide sequences can be found in the folder **"6.Reference_sequences_for_target_regions"**. These two seqeunce files will be used in the extraction of orthologous sequence groups from assembled contigs. For detailed information, please refer to [##A part of analyzing amplicon capture data: Extract orthologous sequence groups from assembled contigs](#extract-orthologous-sequence-groups-from-assembled-contigs)
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/be26d86c-516c-45c3-a4c1-1b44901222cc" alt="Drawing" width="900" height="700"/>
  </div>
  <br /><br />



# A step-by-step instruction for using UPrimer to develop NPCL primers - taking Lepidoptera as an example

Assuming we are now developing a set of NPCL primers for the order Lepidoptera, the steps are as follows:

**Step 1**: Collect available genomic data information for Lepidoptera from the NCBI database (https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/).
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/a0023a7b-bb26-43e4-9837-438046df1458" alt="Drawing" width="800" height="520"/>
  </div>
  <br /><br />


**Step 2**: Select the reference species, ingroup species, and outgroup species for the development of NPCL primers.
- Reference: When selecting the reference species, it is advisable to prioritize species with high-quality genome assembly and annotation. In this case, we have chosen the *Bombyx mori* as the reference species.
<br /><br />

- Ingroups: The selection of ingroup species can be based on the **taxonomy information** of available genome resources. Since our goal is to develop NPCL primer pairs for Lepidoptera, we aim to cover a broad representation of the entire order. Therefore, we can choose representative species at the family or superfamily level to achieve this. In this case, we selected one representative species from each of the 11 families, considering both the quality of genome assembly and annotation (the larger the scaffold N50 value, the genome assemblies with chromosome-level assembly are prioritized for selection).
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/713d5727-ad84-4fa2-8653-160ff2fb213c" alt="Drawing" width="800" height="250"/>
  </div>
  <br /><br />

- Outgroups: When selecting the outgroup species, we need to review the existing literature to identify the outgroup taxa related to our target group. Then, we search for available genome resources of these outgroup taxa in the NCBI database. Finally, we select species with available CDS data as our outgroup representatives. It is generally sufficient to include **4 to 5 outgroup taxa** for our analysis.
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/f620bfdf-525a-4a8a-9193-c2ef721368df" alt="Drawing" width="800" height="130"/>
  </div>
  <br /><br />

**Step 3**: Download the genome data in FASTA format for NPCL primer development (see [Download data](#download-data)). For the reference species, we need to download the proteome, exonome, and genome data for *Bombyx mori*. For the 11 ingroup species, we need to download the genome data. And for the outgroup species, we need to download the CDS data. ***It should be noted that the downloaded genome data from the NCBI database does not require any additional processing, except for modifying the file names (see step 4)***. 
<br /><br />

**Step 4**: Prepare UPrimer input files (see [Input](#input)). We need to start by creating a new folder, which can be named "Lepidoptera_NPCL_primer_development." Inside this new folder, we should first place the UPrimer program package. Then, we need to prepare a text file. The text file should be named "species_list.txt", and it should contain the species names of the reference, ingroup, and outgroup species. Its purpose is to ask UPrimer about the selected reference, ingroup, and outgroup species. In this text table, only species names are displayed, and spaces in the names are replaced with underscores. Finally, we need to prepare three folders, named "Reference," "Ingroups," and "Outgroups". The "Reference" folder contains exonome, genome, and proteome data of the reference species. The "Ingroups" folder contains genome data of all ingroup species and the "Outgroups" folder contains CDS data of all outgroup species. The FASTA data of all species in these three folders need to be labeled with the suffixes "_exon," "_pep," "_genome," and "_cds" to indicate their respective types. 
<br /><br />

**To help the user better understand the workflow of UPrimer, we have provided test data for this case in the *"Test_data_of_developing_NPCL_primers_for_Lepidoptera.zip"* file, which can be accessed at https://github.com/zhangpenglab/UPrimer/. Please note that the test genome data for Lepidoptera has been reduced in size for the purpose of guide the user on how to use UPrimer for designing NPCL primers effectively. The test data files are as shown in the following figure:**
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/f522bfca-39a8-426a-b019-eca403737242" alt="Drawing" width="950" height="700"/>
  </div>
  <br /><br />

**Step 5**: Run UPrimer. The file "Test_data_of_developing_NPCL_primers_for_Lepidoptera.zip" can be found in the UPrimer program package. Please run the following command：

~~~
conda activate UPrimer
cd /path/to/UPrimer-main
unzip Test_data_of_developing_NPCL_primers_for_Lepidoptera.zip
cd /path/to/UPrimer-main/Test_data_of_developing_NPCL_primers_for_Lepidoptera/Lepidoptera_NPCL_primer_development/UPrimer
python UPrimer.py -F species_list.txt -CO 1 -TG Lepidoptera -IT2 0.85 -AT2 0.75 -IT1 0.80 -AT1 0.70 -PId 3
~~~
<br /><br />

**Step 6**: Review and collect the results of the NPCL nested-PCR primer design conducted by UPrimer. The final NPCL primers of the order Lepidoptera can be found in the folder "6.Designed_nested-PCR_primer_set_of_NPCLs_PIratio_3". The reference NPCL nucleotide and peptide sequences can be found in the folder "6.Reference_sequences_for_target_regions".
<br /><br />

**####After the NPCL primer development work is completed, the user can use the primer information to collect amplicons, prepare probes, and conduct sequence capture experiments. The detailed experimental process can be found in the *"Protocol for Amplicon capture.pdf"* (https://github.com/zhangpenglab/UPrimer/tree/main/Accessory)**
<br /><br />
<div class="page-break"></div>
<br /><br />



# Extract orthologous sequence groups from assembled contigs 

The workflow for amplicon capture data analysis is illustrated in the following figure, consisting of four main steps: 1) Capture data assembly; 2) Identification of orthologous sequence groups; 3) Multiple sequence alignment and quality control; 4) Phylogenetic analysis. Among these four steps, the extraction of orthologous sequence groups is crucial and lacks a universal program or software. Therefore, we provide a data analysis script called **'Extract_orthologous_sequence_groups_from_assembled_contigs.py'** (available at https://github.com/zhangpenglab/UPrimer/tree/main/Accessory) specifically designed to extract orthologous sequence groups from assembled contigs. 

To extract target NPCL sequences from the assembled contigs, a specific number of **reference nucleotide and peptide sequences from NPCLs provided by UPrimer** (depending on which NPCLs and how many were captured by the user) can be used as guide sequences. The working principle of this script is roughly as follows: First, **TBLASTN** (e-value < 1e-5, identity > 50%) was performed to identify orthologous contigs based on the reference peptide sequences. Then, a reversed **BLASTN** (e-value < 1e-5, identity > 50%) was performed on the identified orthologous contigs against the reference nucleotide sequence to detect potential chimeras. As the orthologous contigs contained flanking sequences of the target regions, **EXONERATE**  was employed to identify potential intron-exon boundaries based on the reference protein sequence of each target NPCL. Finally, this script combines all identified orthologous exons from all samples, constructing orthologous sequence groups (OGs) at both the DNA and protein levels.

  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/277956e3-6a36-4544-a2bf-d2fc69814a9a" alt="Drawing" width="720" height="400"/>
  </div>
  <br /><br />


## Inputs 
The user should begin by creating a new folder (e.g., Amplicon_capture_data_analysis) and placing the python script inside it. Additionally, the user needs to prepare three data folders to store the respective data files. The details are as follows:
  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/d0fe4ff0-1894-489b-b79d-5ad5f9419047" alt="Drawing" width="500" height="130"/>
  </div>
<br /><br />

(1) **The folder "contigs"**: This folder contains the assembled contig files for all captured samples. These contig files should be named according to the following format: 'sampleID.fasta', for example: 'LepSample1.fasta', 'LepSample2.fasta'. Please ensure that the filenames only consist of the sample ID without any additional characters such as "_" or "|".
  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/d12d26bb-396b-4b13-85a2-091e3c9d05b8" alt="Drawing" width="500" height="400"/>
  </div>
<br /><br />

(2) **The folder "reference"**: This folder contains the reference nucleotide (named ***"ref_nuc.fasta"***) and peptide sequence (named ***"ref_pro.fasta"***) files for the target capture regions. These two sequence files are subsets of "Reference nucleotide sequences.fasta" and "Reference peptide sequences.fasta" generated by UPrimer. The user can select the target NPCLs from the above two files based on the actual captured region's NPCL ID.
  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/816bb167-98a7-4345-8e75-9a42712d9191" alt="Drawing" width="800" height="350"/>
  </div>
<br /><br />

  
(3) **The folder "pep_for_exonerate"**: This folder contains the reference peptide sequences file for the target capture regions. This file is identical to the reference peptide sequence file in the folder "reference".
 
  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/bcf0ef45-3317-4c0c-bc06-417fbcd9ee1d" alt="Drawing" width="500" height="380"/>
  </div>
<br /><br />  

## Outputs
The final orthologous sequence (OG) files can be found in the folder **'6.Final_seq_for_align'**. These files can be used directly for subsequent **sequence alignment and phylogenetic analysis**.
  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/6af59b1f-88c6-4d7b-bb90-4aaef9306659" alt="Drawing" width="800" height="250"/>
  </div>
<br /><br />

## A step-by-step instruction for extracting OGs from contigs - taking Lepidoptera as an example

######To facilitate the user in quickly understanding the usage and processing flow of the script, we provide test data for the following example, which can be found in https://github.com/zhangpenglab/UPrimer/tree/main/Accessory/Test_data_for_extracting_OGs_from_assembled_contigs.zip
<br /><br />

Assuming we captured ***100 NPCLs*** (NPCL ID: Lep1-Lep100) from ***six Lepidoptera samples*** (sample ID: LepSample1-LepSample6), and we have obtained the contigs (in our case, these contigs were assembled using **metaSPAdes**, detailed methods can be found in our paper) for these six samples. The process of extracting orthologous sequence groups from these assembled contigs is as follows:

**Step 1——Prepare reference sequence files** 

From the *"Reference nucleotide sequences.fasta" and "Reference peptide sequences.fasta"* files outputted by UPrimer, extract the nucleotide sequences and protein sequences corresponding to Lep1~Lep100. Rename the extracted nucleotide sequence file as ***"ref_nuc.fasta"*** and the extracted protein sequence file as ***"ref_pro.fasta"***.

  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/edb5c905-34ab-4a00-b5de-3d5961b18922" alt="Drawing" width="800" height="350"/>
  </div>
<br /><br />

**Step 2——Prepare input files** 

Created a new folder named "Lep_amplicon_data_analysis". The script 'Extract_orthologous_sequence_groups_from_assembled_contigs.py' for extracting OGs should be placed within this new folder. Then, three new folders have been created: "contigs", "reference", and "pep_for_exonerate".

- In the "contigs" folder, place the contigs FASTA files for the six Lepidoptera samples. Ensure that the file names follow the format: SampleID.fasta, without the use of "_" and "|".
- In the "reference" folder, place the prepared "ref_nuc.fasta" and "ref_pro.fasta" files from Step 1.
- In the "pep_for_exonerate" folder, place the "ref_pro.fasta" file from Step 1.


  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/ef5ea457-7507-4c04-a7e2-0cefca876da0" alt="Drawing" width="500" height="320"/>
  </div>
<br /><br />

**Step 3——Run the script** 

~~~
conda activate UPrimer
cd /path/to/UPrimer-main/Accessory
unzip Test_data_for_extracting_OGs_from_assembled_contigs.zip
cd /path/to/Lep_amplicon_data_analysis
python Extract_orthologous_sequence_groups_from_assembled_contigs.py -c /path/to/Lep_amplicon_data_analysis/contigs
~~~

**Step 4——Collect the extracted OGs from the output folder '6.Final_seq_for_align' and proceed with the subsequent bioinformatic steps** 

### FAQ
Contact us if you have any questions
