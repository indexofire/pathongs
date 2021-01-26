# Mafft 使用

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! Abstract "内容简介"
    Mafft 是常用的大规模序列比对的工具，在速度和准确度上有很好的平衡。

## 安装

```bash
# conda 安装
(msa)$ conda install mafft
```

## 使用

```bash
# 最基本的序列比对用法
(msa)$ mafft seq.fas > seq.maf
```

PubMLST 现有物种数据库：

- Achromobacter spp.
- Acinetobacter baumannii#1
- Acinetobacter baumannii#2
- Aeromonas spp.
- Aggregatibacter actinomycetemcomitans
- Anaplasma phagocytophilum
- Arcobacter spp.
- Aspergillus fumigatus
- Bacillus cereus
- Bacillus licheniformis
- Bacillus subtilis
- Bartonella bacilliformis
- Bartonella henselae
- Bartonella washoensis
- Bordetella spp.
- Borrelia spp.
- Brachyspira hampsonii
- Brachyspira hyodysenteriae
- Brachyspira intermedia
- Brachyspira pilosicoli
- Brachyspira spp.
- Brucella spp.
- Burkholderia cepacia complex
- Burkholderia pseudomallei
- Campylobacter concisus/curvus
- Campylobacter fetus
- Campylobacter helveticus
- Campylobacter hyointestinalis
- Campylobacter insulaenigrae
- Campylobacter jejuni
- Campylobacter lanienae
- Campylobacter lari
- Campylobacter sputorum
- Campylobacter upsaliensis
- Candida albicans
- Candida glabrata
- Candida krusei
- Candida tropicalis
- Candidatus Liberibacter solanacearum
- Carnobacterium maltaromaticum
- Chlamydiales spp.
- Citrobacter freundii
- Clonorchis sinensis
- Clostridioides difficile
- Clostridium botulinum
- Clostridium perfringens
- Clostridium septicum
- Corynebacterium diphtheriae
- Cronobacter spp.
- Dichelobacter nodosus
- Edwardsiella spp.
- Enterobacter cloacae
- Enterococcus faecalis
- Enterococcus faecium
- Escherichia coli#1
- Escherichia coli#2
- Flavobacterium psychrophilum
- Gallibacterium anatis
- Geotrichum spp.
- Glaesserella parasuis
- Haemophilus influenzae
- Helicobacter cinaedi
- Helicobacter pylori
- Helicobacter suis
- Kingella kingae
- Klebsiella aerogenes
- Klebsiella oxytoca
- Klebsiella pneumoniae
- Kudoa septempunctata
- Lactobacillus salivarius
- Leptospira spp.
- Leptospira spp.#2
- Leptospira spp.#3
- Listeria monocytogenes
- Macrococcus canis
- Macrococcus caseolyticus
- Mannheimia haemolytica
- Melissococcus plutonius
- Moraxella catarrhalis
- Mycobacteria spp.
- Mycobacterium abscessus
- Mycobacterium massiliense
- Mycoplasma agalactiae
- Mycoplasma anserisalpingitidis
- Mycoplasma bovis
- Mycoplasma flocculare
- Mycoplasma gallisepticum#1
- Mycoplasma gallisepticum#2
- Mycoplasma hominis
- Mycoplasma hyopneumoniae
- Mycoplasma hyorhinis
- Mycoplasma iowae
- Mycoplasma pneumoniae
- Mycoplasma synoviae
- Neisseria spp.
- Orientia tsutsugamushi
- Ornithobacterium rhinotracheale
- Paenibacillus larvae
- Pasteurella multocida#1
- Pasteurella multocida#2
- Pediococcus pentosaceus
- Photobacterium damselae
- Piscirickettsia salmonis
- Porphyromonas gingivalis
- Propionibacterium acnes
- Pseudomonas aeruginosa
- Pseudomonas fluorescens
- Pseudomonas putida
- Rhodococcus spp.
- Riemerella anatipestifer
- Salmonella enterica
- Saprolegnia parasitica
- Shewanella spp.
- Sinorhizobium spp.
- Staphylococcus aureus
- Staphylococcus chromogenes
- Staphylococcus epidermidis
- Staphylococcus haemolyticus
- Staphylococcus hominis
- Staphylococcus lugdunensis
- Staphylococcus pseudintermedius
- Stenotrophomonas maltophilia
- Streptococcus agalactiae
- Streptococcus bovis/equinus complex (SBSEC)
- Streptococcus canis
- Streptococcus dysgalactiae equisimilis
- Streptococcus gallolyticus
- Streptococcus oralis
- Streptococcus pneumoniae
- Streptococcus pyogenes
- Streptococcus suis
- Streptococcus thermophilus
- Streptococcus thermophilus#2
- Streptococcus uberis
- Streptococcus zooepidemicus
- Streptomyces spp
- Taylorella spp.
- Tenacibaculum spp.
- Treponema pallidum
- Trichomonas vaginalis
- Ureaplasma spp.
- Vibrio cholerae
- Vibrio cholerae#2
- Vibrio parahaemolyticus
- Vibrio spp.
- Vibrio tapetis
- Vibrio vulnificus
- Wolbachia
- Xylella fastidiosa
- Yersinia pseudotuberculosis
- Yersinia ruckeri

```bash
(ariba)$ ariba pubmlstget "Salmonella enterica" salmonella
(ariba)$ ariba run salmonella/ref_db S1_1.fastq.gz S1_2.fastq.gz output
```


## Reference
