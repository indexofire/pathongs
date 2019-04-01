# Biopython

## 学习材料

- Biopython官方文档
- Biopython官方文档中文

## 使用

```python
from Bio import SeqIO
```

## 用 Biopython 处理数据

```python
from Bio import Seq

```


## 使用 Entrez

```python
# example from dbSNP
import time
from Bio import Entrez

# email address
Entrez.email = "your@email.com"

# API key from NCBI (https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/).
# 10 queries per second a valid API key, otherwise 3 queries per seconds are allowed for 'None'
Entrez.api_key = None

# entrez query (term) can be build and test online using web query builder (https://www.ncbi.nlm.nih.gov/snp/advanced)
# esearch handle
# db: dbsnp
# term: gene LPL
# usehistory: cache result on server for download in batches
# retmax: return max RSID number
eShandle = Entrez.esearch(db="snp",
  term='LPL[All Fields] AND (pathogenic[Clinical_Significance] AND missense[Function_Class])',
  usehistory="y",
  retmax=20)

# get esearch result
eSresult = Entrez.read(eShandle)
# review results
for k in eSresult:
    print (k, ":", eSresult[k])
#Output: Web environment (&WebEnv) and query key (&query_key) parameters specifying the location on the Entrez history server of the list of UIDs matching the Entrez query
#https://www.ncbi.nlm.nih.gov/books/NBK25500/#chapter1.Storing_Search_Results

# get result RSIDs list 'Idlist'
# total rs count
rslist = (eSresult['IdList'])

# retmax = 20 so print only 20 RSIDs
# additional results can be retrieved by batches
# download in batches example http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc139 or see below.
for rs in rslist:
    print(rs)

# get the WebEnv session cookie, and the QueryKey:
webenv = eSresult["WebEnv"]
query_key = eSresult["QueryKey"]
total_count = int(eSresult["Count"])
query_key = eSresult["QueryKey"]
retmax = 5

# sample codes adopted from http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc139.
fetch_count = 0
for start in range(0, total_count, retmax):
    end = min(total_count, start+retmax)
    print("Going to download record %i to %i" % (start+1, end))
    attempt = 0
    #fetch_count += 1
    while (attempt < 3): # & (fetch_count < 2):
        attempt += 1
        try:
            fetch_handle = Entrez.efetch(db="snp",
                                         rettype="uilist", #available types [uilist | docsum (use retmode=xml))
                                         #retmode="xml",
                                         retstart=start,
                                         retmax=retmax,
                                         webenv=webenv,
                                         query_key=query_key )
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print("Received error from server %s" % err)
                print("Attempt %i of 3" % attempt)
                time.sleep(15)
            else:
                raise
    if (fetch_handle):        
        data = fetch_handle.read()
        print(data)
        fetch_handle.close()

```
