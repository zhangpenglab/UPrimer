<p align="center">
<img src="https://user-images.githubusercontent.com/73507779/223490314-5084f3e0-e695-4acc-a513-b36d35d8e630.png" alt="Drawing" width="200"/>
</p>



UPrimer : a program to automatically develop nuclear protein-coding locus (NPCL) amplification primers from genome data for amplicon capture phylogenomics
=================
Jiaxuan Li, Guangcheng Han, Xiao Tian, Dan Liang*, Peng Zhang*<br>
*State Key Laboratory of Biocontrol, School of Life Sciences, Sun Yat-Sen University, Guangzhou, China*


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
- MEGA-CC version 10.18 (https://github.com/Lijiaxuan420/UPrimer/tree/main/Accessory/megacc_10.1.8_amd64.tar.gz)
- PAL2NAL version 14
- Exonerate version 2.4.0 #This software will be used for subsequent capture data analysis.


## A step-by-step installation guide (For Linux)

Before installing, you should to download and install the software Miniconda (https://docs.conda.io/en/latest/miniconda.html) on your system. After that, you can follow the steps below:

(1) Download the UPrimer package.
~~~
git clone https://github.com/Lijiaxuan420/UPrimer.git 
~~~

(2) Install python, biopython, BLAST+, pandas, pal2nal and exonerate through conda.
~~~
unzip /path/to/UPrimer.zip
conda env create -f /path/to/UPrimer/UPrimer-conda-env.yml
~~~

(3) Install megacc and put it on your PATH.

~~~
cd /path/to/UPrimer/Accessory/
tar xvfz megacc_10.1.8_amd64.tar.gz

vim ~/.bashrc                                           
export PATH=$PATH:/path/to/UPrimer/Accessory
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
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/b2d8fb66-49dc-423d-9a5c-32f97d98baef" width="770" height="600"/>
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
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/1240fd1a-52ab-4f72-88ae-68fa093e6a86" alt="Drawing" />
</div>
<br /><br />

- Step 1 —— **Search for primer blocks**: For each alignment, UPrimer searches for all conserved primer blocks that are 7 or 8 amino acids in length. The sequence similarity between the conserved primer block within the ingroup and the reference species is required to be at least 85%. The sequence similarity across all species should be at least 75%.
<br /><br />

- Step 2 —— **Design forward and reverse primers**: Firstly, UPrimer infers the consensus amino acid for each column within the ingroup species and the reference species using the amino acid primer block. Then, nucleotide alignments are used to obtain the consensus nucleotide sequence. Finally, based on this consensus nucleotide sequence, UPrimer designs both forward and reverse primers.
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/6b09c898-c765-4968-93a5-4712126da6d2" alt="Drawing" width="500" height="300"/>
</div>
<br /><br />

- Step 3 —— **Match and filter primer pairs**: UPrimer matches all forward and reverse primers to list all possible primer pairs and filters them by primer degeneracy (the overall degeneracy of the nucleotide primers is less than or equal to 8192) and amplification length (300~2100 bp).
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/c4c990be-4b7f-448b-a4ca-6dbfe8a9a154" alt="Drawing" width="500" height="500"/>
</div>
<br /><br />

- Step 4 —— **Search for outer forward and reverse primers**: For every retained primer pair, UPrimer searches for their outer forward and reverse primers (if they exist) within a flanking region of 450 bp.
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/fb40e784-b2b4-4650-a6a5-d9c38ad51d1a" alt="Drawing" width="500" height="300"/>
</div>
<br /><br />


- Step 5 —— **Obtain a list of all possible nested-PCR primer pairs**
<br /><br />

- Step 6 —— **Score all nested-PCR primer pairs based on their PCR performances and phylogenetic information**: UPrimer scores each nested-PCR primer pair based on its potential PCR performance (named “ScorePCR”), taking into account the conservativeness of primer blocks, primer degeneracy, and primer complexity, as well as the phylogenetic informativeness of the locus (named “ScoreINFOR”), which depends on the variability of the amplification regions (*detailed calculation algorithm see [Scoring algorithms for primers](#scoring-algorithms-for-primers)*). 

- Step 7 —— **Calculate a total score for each primer pair and sort them**: UPrimer then calculates a total score for each primer pair, based on a weighting parameter “PIs” (default = 1) of PCR performance score and phylogenetic information score. The formula for calculation is as follows: 
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/14b33f5d-323a-40f5-b2a9-9cfad8a8b95f" alt="Drawing" width="450" height="60"/>
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
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/2727f1c3-9f2f-4227-9daf-fade2c91debe" alt="Drawing" width="450" height="400"/>
</div>  
<br /><br />

(2) The degeneracy score of amino acid primer sequences (PrimerDEG, maximum score: 100)

Higher degeneracy of amino acid primer sequences can potentially decrease the binding efficiency between the primer and the template, as the wider range of choices for degenerate bases may lead to unstable binding. Therefore, selecting primers with lower degeneracy can improve binding efficiency and amplification specificity. Considering the 5' end of the amino acid primer sequence has a smaller impact on amplification efficiency compared to the 3' end, UPrimer assigns weight ratios to each amino acid position based on its distance from the primer's 3' end. The amino acid site at the 3' end is assigned the highest weight ratio of 1, the amino acid site at the second position from the 3' end has a weight ratio of 0.9, and the amino acid site at the third position from the 3' end has a weight ratio of 0.8. The weight ratios decrease progressively until the last amino acid site. The penalty score for each amino acid site is calculated by multiplying its degeneracy by the weight ratio. The penalty score for each primer is obtained by subtracting the sum of penalty scores for all amino acid sites from 100. UPrimer calculates the degeneracy score for each of the four amino acid primers separately, and then calculates the average score to obtain the ‘PrimerDEG’ score for a pair of nested PCR primers.

##For example, calculating the degeneracy score of a 8 amino acid primer sequence:
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/34bc5cd5-1bf0-4758-8bec-579d250e70d5" alt="Drawing" width="610" height="300"/>
</div>
<br /><br />

(3) The complexity score of amino acid primer sequences (PrimerCOM, maximum score: 100)

The complexity of amino acid primer sequences is determined by evaluating the amino acid composition of the primer sequences. UPrimer employs a counting method to tally the occurrence of each type of amino acid within the primer sequence. UPrimer first applies the formula *N × (N-1)* to calculate the penalty score for each type of amino acid, where N represents the number of occurrences of that particular amino acid. The complexity score for each primer is obtained by subtracting the sum of penalty scores for all types of amino acid from 100. UPrimer calculates the complexity score for each of the four amino acid primers separately, and then calculates the average score to obtain the ‘PrimerCOM’ score for a pair of nested PCR primers.

##For example, calculating the complexity score of a 8 amino acid primer sequence:
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/462c9eec-eeae-42f2-9e0a-b488a08674a2" alt="Drawing" width="500" height="300"/>
</div>
<br /><br />


**@Phylogenetic informativeness: ScoreINFOR, maximum score: 100**

The phylogenetic information of the amplification region is determined by the sequence variability. To convert this value into a percentage scale, UPrimer utilizes the following formula: 
ScoreINFOR = (Infor - InforMin)/(InForMax-InForMin) * (Max_Infor_score - Min_Infor_score) + Min_Infor_score
Where:
- **Infor** represents the phylogenetic information value of a specific amplification region, which is calculated by summing the sequence diveregence in each column of all ingroup species and the reference species. The formula for calculation is as follows: InforSelf = D1 + D2 + D3 + ... + Dn. Here, n represents the length of the amplification region (i.e., the number of column positions). Dn represents the amino acid divergence of all ingroup species and the reference species in each column position.
- **InforMin** represents the **lowest value** of phylogenetic information within different primer pairs' amplification regions in the same MSA.
- **InforMax** represents the **highest value** of phylogenetic information within different primer pairs' amplification regions in the same MSA.
- **Max_Infor_score** represents the maximum score limit, which is 100 points.
- **Min_Infor_score** represents the minimum score limit, which is 50 points.

##For example, calualting the phylogenetic information score of all nested primer pairs in the same MSA.
<div>
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/a1b0ea79-300c-4dd7-be02-65970c7c3934" alt="Drawing" width="1000" height="440"/>
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
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/d2b7287a-c284-47ee-9938-9050eefe9636" alt="Drawing" width="650" height="300"/>
</div>
<br /><br />


(2) Extract **Taxonomy information** based on **Taxonomy ID**:
~~~
taxonkit lineage Genome_list_Name_TaxonomyID.txt -i 2 > Genome_list_Name_TaxonomyID_Lineage.txt
~~~
For example:
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/40d8a9a9-c4df-4cd2-b5eb-06253037d408" alt="Drawing" width="600" height="300"/>
</div>
<br /><br />

(3) Collate lineage information of the target group as the user wants:
~~~
taxonkit reformat Genome_list_Name_TaxonomyID_Lineage.txt -i 3 > Genome_list_Name_TaxonomyID_Lineage_Collated.txt
~~~


(4) Finally, the user needs to obtain a genome list like this:
<div align="center">
  <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/7f322122-5358-4fe2-99d4-597e791c08e5" alt="Drawing" width="800" height="250"/>
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

(1) Download genome, proteome, and genome annotation GB files based on the accession numbers of the reference species *Bombyx mori*:

~~~
ncbi-genome-download all --section refseq --formats fasta -A GCF_014905235.1 --flat-output -o /path/to/Reference -r 100
ncbi-genome-download all --section refseq --formats protein-fasta -A GCF_014905235.1 --flat-output -o /path/to/Reference -r 100
ncbi-genome-download all --section refseq --formats genbank -A GCF_014905235.1 --flat-output -o /path/to/Reference -r 100
~~~

(2) Extract exome data from the genome annotation GB file for the referecne sepecies *Bombyx mori*:
**CAUSION**: the exome of the reference can not be directly downloaded through the ncbi-genome-download tool, so we need to use a custome script named **Extract_exome_from_genbank_file.py** (https://github.com/Lijiaxuan420/UPrimer/tree/main/Accessory) to extract exome data from the annotation GB file. Commands as below: 

~~~
cp /path/to/UPrimer/Accessory/Extract_exome_from_genbank_file.py /path/to/Reference
cd /path/to/Reference
mkdir gbff
cp GCF_014905235.1_Bmori_2016v1.0_genomic.gbff ./gbff
python Extract_exome_from_genbank_file.py
~~~

(3) Download genome sequence data based on the accession numbers of the ingroup species *Plutella xylostella* and *Aricia agestis*:

~~~
ncbi-genome-download all --section refseq --formats fasta -A GCA_932276165.1 --flat-output -o /path/to/Ingroups -r 100
ncbi-genome-download all --section refseq --formats fasta -A GCA_905147365.1 --flat-output -o /path/to/Ingroups -r 100
~~~


(4) Download CDS data based on the accession numbers of the related outgroup species *Bombus terrestris* (Hymenoptera, Apidae) and *Tribolium castaneum* (Coleoptera, Tenebrionidae):
~~~
ncbi-genome-download all --section genbank --formats cds-fasta -A GCF_014905235.1 --flat-output -o /path/to/Outgroups -r 100
ncbi-genome-download all --section genbank --formats cds-fasta -A GCF_000002335.3 --flat-output -o /path/to/Outgroups -r 100
...
~~~

---

# Usage
Here is the help message you get by running the script without parameter:

~~~
Usage :
python UPrimer.py --help

Required parameters：
  --file SRT :A text table containing a list of species names for the reference, ingroups and outgroups (filename: 'species_list.txt') [short: -F]
  --ConsiderOutgroups 0|1 :Whether outgroups are considered when designing universal primers. 0: Do not consider. 1: Consider [short: -CO]
  --TaxaGroup SRT :Name of the target group (for example: Vertebrata, Bivalvia) [short: -TG]  

Optioinal parameters：
  --len INT : The length cutoff for reference exons (DEFAULT: 300) [short: -L]
  --thd INT : Number of threads when running BLAST (DEFAULT: 10) [short: -T]
  --cov FLOAT : Single-copy exon selection criteria: The length coverage between exon and its hit seq (DEFAULT: 0.3) [short: -C]
  --sim FLOAT : Single-copy exon selection criteria: The similarity between exon and its hit seq (DEFAULT: 0.5) [short: -S]
  --SelectScale INT : Value of a criteria used in removing rogue sequences, larger and more strict. This parameter usually does not need to be adjusted (DEFAULT: 1.1)[short: -SS]
  --RemoveRigor INT : Value of a sequence identity cutoff used in removing rogue sequences from MSAs. Set this parameter based on the quality of genome data. Larger value is suitbale for high quality genome data (50 or 55), lower value is suitbale for low quality genome data (30 or 35) (DEFAULT: 55) [short: -RR]
  --exp FLOAT : Value of the expect identity of all species in 7AA/8AA blocks (DEFAULT: 0.5) [short: -E]
  --outgroupNum INT : The minimal number of outgroups (DEFAULT: 1) [short: -ON]
  --ingroupNum INT : The minimal number of ingroups (DEFAULT: 1) [short: -IN]
  --totaltaxanum INT : The minimal number of all species (DEFAULT: 3) [short: -TN]
  --tril INT : The permissible length cutoff of amino acid during trimming sequences. This parameter usually does not need to be adjusted (DEFAULT: 10) [short: -l]
  --trit FLOAT : The ratio of gap length and AA length during trimming sequences. This parameter usually does not need to be adjusted (DEFAULT: 1.5) [short: -t]
  --con FLOAT : The conservation/identitfy of MSAs. If similarity of one MSA is larger than this value, UPrimer will not design primer for it (DEFAULT: 0.9) [short: -C]
  --PIs INT : The weighting ratio of PCR score and Information score when scoring all possible nested-PCR primers within the same MSA. Increasing this value will reduce the length of the amplified regions. Please set this parameter with caution. (DEFAULT: 1) [short: -PIs]
  --fwd INT : The length of amino acid between outer primer F1 and inner primer F2. This parameter usually does not need to be adjusted (DEFAULT: 150) [short: -W]
  --rvs INT : The length of amino acid between outer primer R1 and inner primer R2. This parameter usually does not need to be adjusted (DEFAULT: 150) [short: -R]
  --max INT : The maximal length of target amplicons (amino acid length, from F2 to R2) (DEFAULT: 700) [short: -M]
  --min INT : The minimal length of target amplicons (amino acid length, from F2 to R2) (DEFAULT: 100) [short: -N]
  --all INT : The maximal degeneracy level of primers. This parameter usually does not need to be adjusted (DEFAULT: 8192) [short: -A]
  --inthred2 FLOAT : The mean identity of ingroups (including reference species) in 7AA/8AA F2R2 primer blocks (DEFAULT: 0.85) [short: -IT2]
  --allthred2 FLOAT : The mean identity of all species in 7AA/8AA F2R2 primer blocks (DEFAULT: 0.75) [short: -AT2]
  --inthred1 FLOAT : The mean identity of ingroups (including reference species) in 7AA/8AA F1R1 primer blocks (DEFAULT: 0.80) [short: -IT1]
  --allthred1 FLOAT : The mean identity of all species in 7AA/8AA F1R1 primer blocks (DEFAULT: 0.70) [short: -AT1]
  --PId INT : The weighting ratio of PCR score and Information score when scoring highest-scoring nested-PCR primers derived from different MSAs (DEFAULT: 1) [short: -PId]

Example for developing universal primers for the target group at a order or class level:
python UPrimer.py -F species_list.txt -CO 1 -TG TargetGroupName -IT2 0.85 -AT2 0.75 -IT1 0.80 -AT1 0.70 -PId 3

Example for developing universal primers for the target group at a subphylum or phylum level:
python UPrimer.py -F species_list.txt -CO 0 -TG TargetGroupName -IT2 0.90 -AT2 0.90 -IT1 0.90 -AT1 0.90 -PId 5

Example for developing universal primers for the target group whose quality of genome data is relative low:
python UPrimer.py -F species_list.txt -CO 1 -TG TargetGroupName -RR 30 -IT2 0.85 -AT2 0.75 -IT1 0.80 -AT1 0.70 -PId 3

Example for running the module of obtaining candidate MSAs from genome data:
python Make_MSAs_suitable_for_primer_design_Part_I.py -F species_list.txt 

Example for running the module of designing universal primers based on candidate MSAs:
python Design_universal_primer_sets_Part_II.py -F species_list.txt -D 4.Candidate_peptide_iden0.5_MSAs_for_primer_design -TG TargetGroupName -IT2 0.85 -AT2 0.75 -IT1 0.80 -AT1 0.70 -PId 3
~~~
<br /><br />


# Input
The user should begin by creating a new folder (e.g., XXX_universal_primer_development) and placing the UPrimer program package inside it. Additionally, the user need to prepare a text file and three data folders to store the respective data files. The details are as follows:


(1) Text file: This file should contain a list of species names, including the names of the reference species, ingroup species and outgroup species.
  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/b6201808-791f-408c-b8d7-7fc3b425126a" alt="Drawing" width="460" height="400"/>
  </div>
<br /><br />


(2) The folder "Reference": it contains the genomic resources of the reference species, consisting of three files: exonome data (reference_name_exon.fasta), genome data (reference_name_genome.fasta), and proteome data (reference_name_pep.fasta).
  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/efd73dd3-803d-4b6d-af54-7d1c69bfa85a" alt="Drawing" width="350" height="110"/>
  </div>    
<br /><br />

(3) The folder "Ingroups": it contains the genomic data of all ingroups species. The naming convention for these files is as follows: species_name + _genome.fasta.
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/28f8bb28-c120-4c43-8f52-4fdad0de9316" alt="Drawing" width="420" height="240"/>
  </div>
<br /><br />    

(4) The folder "Outgroups": it contains the CDS data of all outgroup species. The naming convention for these files is as follows: species_name + _cds.fasta.
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/efb4020c-1ff3-450d-8f6a-a5dcb6560696" alt="Drawing" width="370" height="110"/>
  </div>
<br /><br />


# Output
(1) The candidate nucleotide and peptide MSAs used for primer design can be found in the folder "4.Candidate_nucleotide_iden0.5_MSAs_for_primer_design" and "4.Candidate_peptide_iden0.5_MSAs_for_primer_design".
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/004236aa-d2f4-429a-b864-96b66971583f" alt="Drawing" width="800" height="280"/>
  </div>
<br /><br />

(2) The final designed NPCL primer sets can be found in the folders "6.Designed_nested-PCR_primer_set_of_NPCLs_PIratio_xxx". In this folder, the user could find three tables:
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/aa776f66-7acf-4e05-a8cd-fdff4882fe04" alt="Drawing" width="750" height="320"/>
  </div>
<br /><br />

- **Highest-scoring_nested-PCR_primer_set_PIdXX_numXXX_mlenXXX.xls**: This table presents detailed information for the highest-scoring nested primers obtained from various MSAs, including such as NPCL ID, Target region length, primer sequence, primer position, block identity. It is important to note that the highest-scoring primers for different MSAs have already be sorted based on their total scores (from high to low). Researchers should prioritize using the primer pairs that are ranked higher on the list.

  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/9b96714d-860c-4a8b-bb32-857ec5a701af" alt="Drawing" width="800" height="310"/>
  </div>
  <br /><br />

- **Synthezised_highest-scoring_nested-PCR_primer_PIdXXX_numXXX_mlenXXX.xls**: This table is a simplified version of the "Highest-scoring_nested-PCR_primer_set_PIdXX_numXXX_mlenXXX.xls" table, containing only the NPCL ID, nested PCR primer sequences, and primer lengths. It can be directly used for primer synthesis.
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/9b352fca-dc33-4d64-a9ee-a33904f009b2" alt="Drawing" width="650" height="550"/>
  </div>
  <br /><br />

- **Candidate_nested-PCR_primer_set_PIdXXX_numXXX_mlenXXX.xls**: If PCR amplification using the highest-scoring primer pair fails, users can explore other candidate primers listed in the table. Each MSA provides ten candidate nested PCR primer pairs. These candidate primers are sorted based on their total scores, and the author can try them sequentially, one by one, to identify a successful amplification.
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/3f4fe1b9-8fc9-4595-a6d6-cfe6e69dfb78" alt="Drawing" width="800" height="310"/>
  </div>
  <br /><br />

(3) The reference NPCL nucletide and peptide sequences can be found at the folder **"6.Reference_sequences_for_target_regions"**. These two seqeunce sets will be used in the extraction of orthologous sequence groups from assembled contigs. For detailed information, please refer to [Additional section: Extract orthologous sequence groups from amplicon capture data](#extract-orthologous-sequence-groups-from-amplicon-capture-data)
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/7cd18bd9-731a-4af4-b3a2-ca58a9d0f409" alt="Drawing" width="750" height="400"/>
  </div>
  <br /><br />


# A step-by-step instruction for using UPrimer to develop NPCL primers - taking Lepidoptera as an example

Assuming we are now developing a set of NPCL primers for the order Lepidoptera, the steps are as follows:

**Step 1**: Collect available genomic data information for Lepidoptera from the NCBI database (https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/).
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/600ae008-4bcd-40f9-99fa-f7b0f565e9d9" alt="Drawing" width="800" height="520"/>
  </div>
  <br /><br />


**Step 2**: Select the reference species, ingroup species, and outgroup species for the development of NPCL primers.
- Reference: When selecting the reference species, it is advisable to prioritize species with high-quality genome assembly and annotation. In this case, we have chosen the *Bombyx mori* as the reference species.
- Ingroups: The selection of ingroup species can be based on the **taxonomy information** of available genome resources. Since our goal is to develop NPCL primer pairs for Lepidoptera, we aim to cover a broad representation of the entire order. Therefore, we can choose representative species at the family or superfamily level to achieve this. In this case, we selected one representative species from each of the 11 families, considering both the quality of genome assembly and annotation.
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/e9f994cc-1db0-4169-88ed-adda89e345d9" alt="Drawing" width="800" height="250"/>
  </div>
  <br /><br />

- Outgroups: When selecting the outgroup species, we need to review the existing literature to identify the outgroup taxa related to our target group. Then, we search for available genome resources of these outgroup taxa in the NCBI database. Finally, we select species with available CDS data as our outgroup representatives. It is generally sufficient to include 4 to 5 outgroup taxa for our analysis.
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/1cfb22af-caea-46fc-94a7-0e507a08841c" alt="Drawing" width="800" height="130"/>
  </div>
  <br /><br />

**Step 3**: Download the genome data in FASTA format for NPCL primer development (see [Download data](#download-data)). For the reference species, we need to download the proteome, exonome, and genome data for *Bombyx mori*. For the 11 ingroup species, we need to download the genome data. And for the outgroup species, we need to download the CDS data. It should be noted that the downloaded genome data from the NCBI database does not require any additional processing, except for modifying the file names (see step 4). 


**Step 4**: Prepare UPrimer input files (see [Input](#input)). We need to start by creating a new folder, which can be named "Lepidoptera_NPCL_primer_development." Inside this new folder, we should place the UPrimer program package. Then, we need to prepare a text file and three data folders. The text file contains the species names of the reference, ingroup, and outgroup species. Its purpose is to inform UPrimer about the selected reference, ingroup, and outgroup species. The three data folders correspond to "Reference," "Ingroups," and "Outgroups," respectively. We need to place the corresponding resources in these folders. Pay attention to the naming conventions of the files. In the text file, only the species names are displayed, with spaces replaced by underscores. In the folders, file names include additional data attributes such as "_exon", "_pep", "_genome", and "_cds". 

**To help users better understand the workflow of UPrimer, we have provided test data for this case in the *"Example.zip"* file, which can be accessed at https://github.com/Lijiaxuan420/UPrimer/. Please note that the test genome data for Lepidoptera has been reduced in size for the purpose of guide users on how to use UPrimer for designing NPCL primers effectively. The data files in the Example folder are as shown in the following figure:**
  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/5a5e95c2-e303-4a4f-9e77-2280bd0163b7" alt="Drawing" width="900" height="700"/>
  </div>
  <br /><br />

**Step 5**: Run UPrimer. The user can download the Example.zip file and run these command lines are as follows:
~~~
conda activate UPrimer
cd /path/to/UPrimer
unzip Example.zip
cd /path/to/Example/Lepidoptera_NPCL_primer_development/UPrimer
python UPrimer.py -F species_list.txt -CO 1 -TG Lepidoptera -IT2 0.85 -AT2 0.75 -IT1 0.80 -AT1 0.70 -PId 3
~~~

**Step 6**: Review and collect the results of the NPCL nested PCR primer design conducted by UPrimer. The final NPCL primer sets can be found at the folder 6.Designed_nested-PCR_primer_set_of_NPCLs_PIratio_3". The reference nucleotide and peptide sequences cam be found at the folder "6.Reference_sequences_for_target_regions".
<br /><br />

**####After the NPCL primer development work is completed, the user can use the primer information to collect amplicons, prepare probes, and conduct sequence capture experiments. The detailed experimental process can be found in the *"Protocol for Amplicon capture.pdf"* (https://github.com/Lijiaxuan420/UPrimer/tree/main/Accessory)**
<br /><br />
<div class="page-break"></div>
<br /><br />



# Extract orthologous sequence groups from assembled contigs 

The workflow for amplicon capture data analysis is illustrated in the following diagram, consisting of four main steps: 1) Capture data assembly, 2) Identification of orthologous sequence groups, 3) Multiple sequence alignment and quality control, and 4) Phylogenetic analysis. Among these four steps, the extraction of orthologous sequence groups is crucial and lacks a universal program or software. Therefore, we provide a data analysis script called **'Extract_orthologous_sequence_groups_from_assembled_contigs.py'** (available at https://github.com/Lijiaxuan420/UPrimer/tree/main/Accessory) specifically designed to extract orthologous sequence groups from assembled contigs. 

To extract target NPCL sequences from the assembled contigs, a specific number of **reference nucleotide and peptide sequences from NPCLs provided by UPrimer** (depending on which NPCLs and how many were captured by the user) can be used as guide sequences. The working principle of this script is roughly as follows: First, **TBLASTN** (e-value < 1e-5, identity > 50%) was performed to identify orthologous contigs based on the reference peptide sequences. Then, a reversed **BLASTN** (e-value < 1e-5, identity > 50%) was performed on the identified orthologous contigs against the reference nucleotide sequence to detect potential chimeras. As the orthologous contigs contained flanking sequences of the target regions, **EXONERATE version 2.4.0**  was employed to identify potential intron-exon boundaries based on the reference protein sequence of each target NPCL. Finally, this script combines all identified orthologous exons from all samples, constructing orthologous sequence groups (OGs) at both the DNA and protein levels.

  <div align="center">
     <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/1a4ce2c6-ad7e-46b1-8b72-4fa95c4279e1" alt="Drawing" width="700" height="400"/>
  </div>
  <br /><br />


## Inputs 
The user should begin by creating a new folder (e.g., Amplicon_capture_data_analysis) and placing the python script inside it. Additionally, the user need to prepare three data folders to store the respective data files. The details are as follows:
  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/962516c4-e8ed-4517-a72d-911460e8eff3" alt="Drawing" width="500" height="130"/>
  </div>
<br /><br />

(1) **The folder "contigs"**: it contains the assembled contigs files for all captured samples. The contigs files should be named according to the following format: 'sampleID.fasta', for example: 'LepSample1.fasta', 'LepSample2.fasta'. Please ensure that the filenames only consist of the sample ID without any additional characters such as "_" or "|".
  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/687ff7b7-567a-4297-9287-885de2d6d167" alt="Drawing" width="500" height="400"/>
  </div>
<br /><br />

(2) **The folder "reference"**: it contains the reference nucleotide (named "ref_nuc.fasta") and peptide sequence (named "ref_pro.fasta") files for the target capture regions. These two sequence files are subsets of "Reference nucleotide sequences.fasta" and "Reference peptide sequences.fasta" generated by UPrimer. The user can select the target NPCLs from the above two files based on the actual captured region's NPCL ID.
  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/d7556af7-0450-4c73-b8ae-bb92518482e3" alt="Drawing" width="800" height="350"/>
  </div>
<br /><br />

  
(3) **The folder "pep_for_exonerate"**: it contains the reference peptide sequences file for the target capture regions. This file is identical to the reference peptide sequence file in the "reference" folder.
 
  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/72c3ed10-b0dc-4fb6-923e-fdd4d23b630c" alt="Drawing" width="500" height="380"/>
  </div>
<br /><br />  

## Outputs
The final orthologous sequence (OG) files can be found in the **'6.Final_seq_for_align'** folder. These files can be used directly for subsequent **sequence alignment and phylogenetic analysis**.
  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/c2596d94-375b-4f51-9f06-69fd5b421f9f" alt="Drawing" width="800" height="250"/>
  </div>
<br /><br />

## A step-by-step instruction for extracting OGs from contigs - taking Lepidoptera as an example

######To facilitate the user in quickly understanding the usage and processing flow of the script, we provide test data for the following example, which can be found in https://github.com/Lijiaxuan420/UPrimer/tree/main/Accessory/Test_data_for_extracting_OGs_from_assembled_contigs.zip
<br /><br />

Assuming we captured ***100 NPCLs*** (NPCL ID: Lep1-Lep100) from ***six Lepidoptera samples*** (sample ID: LepSample1-LepSample6), and we have obtained the contigs (in our case, these contigs were assembled using **metaSPAdes**, detailed methods can be found in our paper) for these six samples. The process of extracting orthologous sequence groups from these assembled contigs is as follows:

**Step 1-Prepare reference sequence files** 

From the *"Reference nucleotide sequences.fasta" and "Reference peptide sequences.fasta"* files outputted by UPrimer, extract the nucleotide sequences and protein sequences corresponding to Lep1~Lep100. Rename the extracted nucleotide sequence file as ***"ref_nuc.fasta"*** and the extracted protein sequence file as ***"ref_pro.fasta"***.

  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/c1e22c0b-ca88-4e1b-b103-5548979537cb" alt="Drawing" width="800" height="350"/>
  </div>
<br /><br />

**Step 2-Prepare input files** 

Create a new folder named "Lep_amplicon_data_analysis". Within this folder, create three new folders named "contigs", "reference", and "pep_for_exonerate".

- In the "contigs" folder, place the contigs FASTA files for the six Lepidoptera samples. Ensure that the file names follow the format: SampleID.fasta, without the use of "_" and "|".
- In the "reference" folder, place the prepared "ref_nuc.fasta" and "ref_pro.fasta" files from Step 1.
- In the "pep_for_exonerate" folder, place the "ref_pro.fasta" file from Step 1.
- Additionally, the script 'Extract_orthologous_sequence_groups_from_assembled_contigs.py' for extracting OGs should also be placed within the "Lep_amplicon_data_analysis" folder.

  <div align="center">
    <img src="https://github.com/zhangpenglab/UPrimer/assets/139540726/90d31593-fbbf-4a27-b9c1-0384215e37bd" alt="Drawing" width="500" height="320"/>
  </div>
<br /><br />

**Step 3-Run the script** 

~~~
conda activate UPrimer
cd /my/complete/path/to/Lep_amplicon_data_analysis
python Extract_orthologous_sequence_groups_from_assembled_contigs.py -c /my/complete/path/to/Lep_amplicon_data_analysis/contigs
~~~

**Step 4-Collect the extracted OGs from the output folder '6.Final_seq_for_align' and proceed with the subsequent bioinformatic steps** 

### FAQ
Contact us if you have any questions
