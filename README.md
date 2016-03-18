# Eutils-Tools
Use of the Eutils from NCBI to retrieve genomes using Biopython

##Requirements
- Biopython
- A valid email address (required by NCBI)
- Enough disk space for the genomes to download to
- Knowledge on the length of the required genome


##Purpose
To retrieve all genomes (WGS and complete) from refseq

#### Example
`WGS-fetch.py -q Escherichia+coli -e you@example.com -o /some/directory -l 4-7`


#### Usage
```
usage: `WGS-fetch.py [-h] [--version] -q QUERY -e EMAIL -o OUTPUT -l LENGTH [-c] [-s START]``

Download genomes for organism

optional arguments:
  -h, --help            show this help message and exit
  
  --version             show program's version number and exit
  
  -q QUERY, --query QUERY
  
                        Query for genome database separated by plus sign(s)
                        
  -e EMAIL, --email EMAIL
  
                        A valid email address is required
                        
  -o OUTPUT, --output OUTPUT
  
                        Specify output directory
  -l LENGTH, --length LENGTH
  
                        The range of length for the full genome, the default is 4-7 Mb for E.coli. The default a range in megabases
                        
  -c, --chromosome      Download only complete genomes
  
  -s START, --start START
  
                        Specify start location if downloaded is interrupted


#### For a more stringent retrieval

usage: `WGS-fetchv2.py [-h] [--version] -q QUERY -e EMAIL -o OUTPUT -l LENGTH [-c CONTIGS] [-f COVERAGE]`

Download genomes for organism

optional arguments:
  -h, --help            show this help message and exit
  
  --version             show program's version number and exit
  
  -q QUERY, --query QUERY
  
                        Query for genome database separated by plus sign(s)
                        
  -e EMAIL, --email EMAIL
  
                        A valid email address is required
                        
  -o OUTPUT, --output OUTPUT
  
                        Specify output directory
                        
  -l LENGTH, --length LENGTH
  
                        The range of length for the full genome, the default is 4-7 Mb for E.coli. The default a range in megabases
  -c CONTIGS, --contigs CONTIGS
  
                        Upper limit of contig quatity (default = 250)
  -f COVERAGE, --coverage COVERAGE
  
                        Lower limit of coverage for genome (default = 10.0)
```
