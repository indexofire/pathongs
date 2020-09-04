# Other Format

## Genomatix annotation syntax

Some Genomatix tools, e.g. Gene2Promoter or GPD allow the extraction of sequences. Genomatix uses the following syntax to annotate sequence information: each information item is denoted by a keyword, followed by a "=" and the value. These information items are separated by a pipe symbol "|".
The keywords are the following:

|keyword|meaning|
|-|-|
|loc|The Genomatix Locus Id, consisting of the string "GXL_" followed by a number.|
|sym|The gene symbol. This can be a (comma-separated) list.
|geneid|The NCBI Gene Id. This can be a (comma-separated) list.
|acc|A unique identifier for the sequence. E.g. for Genomatix promoter regions, the Genomatix Promoter Id is listed in this field.|
|taxid|The organism's Taxon Id|
|spec|The organism name|
|chr|The chromosome within the organism.|
|ctg|The NCBI contig within the chromosome.|
|str|Strand, (+) for sense, (-) for antisense strand.|
|start|Start position of the sequence (relative to the contig).|
|end|End position of the sequence (relative to the contig).|
|len|Length of the sequence in basepairs.|
|tss|A (comma-separated list of) UTR-start/TSS position(s). If there are several TSS/UTR-starts, this means that several transcripts share the same promoter (e.g. when they are splice variants). The positions are relative to the promoter region.|
|probe|A (comma-separated list of) Affymetrix Probe Id(s).|
|unigene|A (comma-separated list of) UniGene Cluster Id(s).|
|homgroup|An identifier (a number) for the homology group (available for promoter sequences only). Orthologously related sequences have the same value in this field.|
|promset|If the sequence is a promoter region, the promoter set is denoted here.|
|descr|The gene description. If several genes (i.e. NCBI gene ids) are associated with the sequence, the descriptions for all of the genes are note, separated by ";"|
|comm|A comment field, used for additional annotation. For promoter sequences, this field contains information about the transcripts associated with the promoter. For each transcript the Genomatix Transcript Id, accession number, TSS position and quality is listed, separated by "/". For Genomatix CompGen promoters no transcripts are assigned, in this case the string "CompGen promoter" is denoted.|

This syntax is currently used only for sequences in the FASTA and GenBank formats.

Example (a promoter sequence in GenBank format):

```
LOCUS       GXP_170357    743 bp    DNA
DEFINITION  loc=GXL_141619|sym=TPH2|geneid=121278|acc=GXP_170357|
            taxid=9606|spec=Homo sapiens|chr=12|ctg=NC_000012|str=(+)|
            start=70618393|end=70619135|len=743|tss=501,632|
            homgroup=4612|promset=1|descr=tryptophan hydroxylase 2|
            comm=GXT_2756574/AK094614/632/gold;
            GXT_2799672/NM_173353/501/bronze
ACCESSION   GXP_170357
BASE COUNT    216 a  180 c  147 g  200 t
ORIGIN
        1 TTGATTACCT TATTTGATCA TTACACATTG TACGCTTGTG TCAAAATATC ACATGTGCCT
       61 TATAAATGTG TACAACTATT AGTTATCCAT AAAAATTAAA AATTAAAAAA TCCGTAAAAT
      121 GGTTTAAGCA TTCAGCAGTG CTGATCTTTC TTAAATTATT TTTCTAATTT TGGAAAGAAA
      181 GCACAAAATC TTTGAATTCA CAATTGCTTA AAGACTGAGG TTAACTTGCC AGTGGCAGGC
      241 TTGAGAGATG AGAGAACTAA CGTCAGAGGA TAGATGGTTT CTTGTACAAA TAACACCCCC
      301 TTATGTATTG TTCTCCACCA CCCCCGCCCA AAAAGCTACT CGACCTATGA AACAAATCAC
      361 ACTATGAGCA CAGATAACCC CAGGCTTCAG GTCTGTAATC TGACTGTGGC CATCGGCAAC
      421 CAGAAATGAG TTTCTTTCTA ATCAGTCTTG CATCAGTCTC CAGTCATTCA TATAAAGGAG
      481 CCCGGGGATG GGAGGATTCG CATTGCTCTT CAGCACCAGG GTTCTGGACA GCGCCCCAAG
      541 CAGGCAGCTG ATCGCACGCC CCTTCCTCTC AATCTCCGCC AGCGCTGCTA CTGCCCCTCT
      601 AGTACCCCCT GCTGCAGAGA AAGAATATTA CACCGGGATC CATGCAGCCA GCAATGATGA
      661 TGTTTTCCAG TAAATACTGG GCACGGAGAG GGTTTTCCCT GGATTCAGCA GTGCCCGAAG
      721 AGCATCAGCT ACTTGGCAGC TCA
//
```

#### IUPAC nucleic acid codes

To represent ambiguity in DNA sequences the following letters can be used (following the rules of the International Union of Pure and Applied Chemistry (IUPAC)):

```
A = adenine
C = cytosine
G = guanine
T = thymine
U = uracil
R = G A (purine)
Y = T C (pyrimidine)
K = G T (keto)
M = A C (amino)
S = G C
W = A T
B = G T C
D = G A T
H = A C T
V = G C A
N = A G C T (any)
```
