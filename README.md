# RStudio project describing the data processing for the Alpalga Zonation manuscript

An overview of the data, scripts, and results can be consulted at the following address

<center>https://alpalga.github.io/Zonation/</center>

## Directory contents

### Data

#### Sites metadata

- `orchamp_env.csv` : Environmental measures describing the sampling sites
- `label_metabar_trie.csv` : Metadata describing the sampling location.

#### Bioclimatic data

- `climatic_only2016.Rdata` : the estimated bioclimatic variable used for the ecological 
                  analysis

#### EMBL release 140 related data  

Every files are grouped into the `embl-140` directory
In that directory every files are compressed using gzip and must be uncompressed 
to be usable

```bash
pushd embl-140
gunzip *
popd
```


#### The curated data sets

The curated data sets (see the filtering repports below)
are stored into the `cleaned_datasets` directory

They are all coma separated values files (.csv)

- For the Chlo01 marker
  + `Chlo01.cleanned.motus.csv`   : the MOTU descriptions table
  + `Chlo01.cleanned.samples.csv` : the sample descriptions table
  + `Chlo01.cleanned.reads.csv`   : the MOTUs x sample contengency table
- For the Chlo02 marker
  + `Chlo02.cleanned.motus.csv`   : the MOTU descriptions table
  + `Chlo02.cleanned.samples.csv` : the sample descriptions table
  + `Chlo02.cleanned.reads.csv`   : the MOTUs x sample contengency table
- For the Euka03 marker
  + `euka03.cleanned.motus.csv`   : the MOTU descriptions table
  + `euka03.cleanned.samples.csv` : the sample descriptions table
  + `euka03.cleanned.reads.csv`   : the MOTUs x sample contengency table

#### Raw data

The complete raw data are too large to be placed on GitHub, and cannot be easily deposited on the SRA servers because many samples are multiplexed in each sequencing library. Files containing the most possible raw data were produced from the raw sequencing files and uploaded to the GitHub server. They consist of three fasta files annotated following the *OBITools* format. The three steps of processing that allowed their construction are : 

- The initial forward and reverse files were assembled by `illuminapairedends`
- The resulting file was dereplicated using `obiuniq`
- The sequences appearing only once in the whole experiment (singletons) were discarded using `obigrep`.
  
The three *almost* raw data sequence files compressed by bzip2 are stored in the directory corresponding to each marker:

- `chlo01.no_singleton.fasta.bz`
- `chlo02.no_singleton.fasta.bz`
- `euka03.no_singleton.fasta.bz`

Raw data statistics are available for each markers.

- `chlo01`
  + `chlo01.length.stat` : stats on the amplicon lenght
  + `chlo01.max_per_sample.stat` : stats on the maximum occurrence of a MOTUs in a PCR
- chlo02
  + `chlo02.length.stat` : stats on the amplicon lenght
  + `chlo02.max_per_sample.stat` : stats on the maximum occurrence of a MOTUs in a PCR
- euka03
  + `euka03.length.stat` : stats on the amplicon lenght
  + `euka03.max_per_sample.stat` : stats on the maximum occurrence of a MOTUs in a PCR

### Data analysis

#### Preprocessing of the raw data by OBITools

`obitools_processing.sh` resumes the obitools commands used to process to raw fastq files
to a preliminary MOTU table per PCR.

#### Filtering of the raw data 

Resumes the filtering of the preliminary MOTU table per PCR produced
by obitools to build the MOTU table used for ecological analysis

A RMarkdown resumes the filtering for each marker

- `filtering_chlo01.Rmd`
- `filtering_chlo02.Rmd`
- `filtering_euka03.Rmd`

And the corresponding HTML files (the ones to look over)

- [`filtering_chlo01.html`](filtering_chlo01.html)
- [`filtering_chlo02.html`](filtering_chlo02.html)
- [`filtering_euka03.html`](filtering_euka03.html)

The gererated curated data tables are stored in the
directory `cleaned_datasets`

