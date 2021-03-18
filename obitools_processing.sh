##
## Treatement of the raw sequencing files
## Using OBITools : http://metabarcoding.org/obitools
##

#########
#
# Files renaming
#
#

mv chlo01/GWM-1383.chlo1.filtered.uniq.fasta chlo01/chlo01.fasta
mv chlo02/GWM-1384.filtered.uniq.fasta chlo02/chlo02.fasta
mv euka03/GWM-1107.filtered.uniq.fasta euka03/euka03.fasta


#########
#
# Makes stats about sequence variants length
#
#

obistat -c seq_length  chlo01/chlo01.fasta >  chlo01/chlo01.length.stat
obistat -c seq_length  chlo02/chlo02.fasta >  chlo02/chlo02.length.stat
obistat -c seq_length  euka03/euka03.fasta >  euka03/euka03.length.stat


#########
#
# Keeps only good length sequences
#
#

obigrep -l 65 -L 190 chlo01/chlo01.fasta > chlo01/chlo01.good_length.fasta
obigrep -l 65 -L 130 chlo02/chlo02.fasta > chlo02/chlo02.good_length.fasta
obigrep -l 65 -L 200 euka03/euka03.fasta > euka03/euka03.good_length.fasta


#########
#
# Removes singletons (variants occuring as a single read)
#
#

obigrep -p "count > 1" chlo01/chlo01.good_length.fasta > chlo01/chlo01.no_singleton.fasta
### number of variants :    130646
### number of reads    :   7658050

obigrep -p "count > 1" chlo02/chlo02.good_length.fasta > chlo02/chlo02.no_singleton.fasta
###  number of variants :   172078
### number of reads     :  8316740

obigrep -p "count > 1" euka03/euka03.good_length.fasta > euka03/euka03.no_singleton.fasta
###  number of variants :   402999
### number of reads     : 11399428

#########
#
# Removes sequence variants with non-ACGT nucleotides
#
#

obigrep -s '^[acgt]+$' chlo01/chlo01.no_singleton.fasta > chlo01/chlo01.onlyACGT.fasta
###  number of variants :   129811
### number of reads     :  7655639

obigrep -s '^[acgt]+$' chlo02/chlo02.no_singleton.fasta > chlo02/chlo02.onlyACGT.fasta
### number of variants :    171933
### number of reads    :   8315854

obigrep -s '^[acgt]+$' euka03/euka03.no_singleton.fasta > euka03/euka03.onlyACGT.fasta
### number of variants :    390802
### number of reads    :  11246979


#########
#
# Makes stats the maximum occurrency of a sequence variant
# in a PCR
#
#

obiannotate -S 'max_per_sample:max(merged_sample.values())' \
            chlo01/chlo01.onlyACGT.fasta  \
    | obistat -c max_per_sample            \
    | sort -n  \
    > chlo01/chlo01.max_per_sample.stat

obiannotate -S 'max_per_sample:max(merged_sample.values())' \
            chlo02/chlo02.onlyACGT.fasta  \
    | obistat -c max_per_sample             \
    | sort -n  \
    > chlo02/chlo02.max_per_sample.stat

obiannotate -S 'max_per_sample:max(merged_sample.values())' \
            euka03/euka03.onlyACGT.fasta  \
    | obistat -c max_per_sample           \
    | sort -n \
    > euka03/euka03.max_per_sample.stat


#########
#
# Keeps only variants occuring at least 10 time in at least one PCR
#
#


obiannotate -S 'max_per_sample:max(merged_sample.values())' \
            chlo01/chlo01.onlyACGT.fasta  \
    | obigrep -p 'max_per_sample >= 10'   \
    | sort -n  \
    > chlo01/chlo01.min_per_sample_10.fasta

### number of variants :     16187
### number of reads    :   6910562

obiannotate -S 'max_per_sample:max(merged_sample.values())' \
            chlo02/chlo02.onlyACGT.fasta  \
    | obigrep -p 'max_per_sample >= 10'   \
    | sort -n  \
    > chlo02/chlo02.min_per_sample_10.fasta

### number of variants :     33421
### number of reads    :   7698821


obiannotate -S 'max_per_sample:max(merged_sample.values())' \
            euka03/euka03.onlyACGT.fasta  \
    | obigrep -p 'max_per_sample >= 10'   \
    | sort -n  \
    > euka03/euka03.min_per_sample_10.fasta

### number of variants :     17490
### number of reads    :   8798319

#########
#
# Runs OBIClean to tag PCR point mutations.
#
#

obiclean -C -r 0.2 -s merged_sample chlo01/chlo01.onlyACGT.fasta \
         > chlo01/chlo01.obiclean.fasta

obiclean -C -r 0.2 -s merged_sample chlo02/chlo02.onlyACGT.fasta \
         > chlo02/chlo02.obiclean.fasta

obiclean -C -r 0.2 -s merged_sample euka03/euka03.onlyACGT.fasta \
         > euka03/euka03.obiclean.fasta

#########
#
# Runs taxonomic assignment
#
#

ecotag -d embl-140/ncbi20190930 \
       -R embl-140/chlo01.uniq.by_taxid.clean.uid.fasta \
       -m 0.75 \
       chlo01/chlo01.obiclean.fasta > chlo01/chlo01.ecotag.fasta

ecotag -d embl-140/ncbi20190930 \
       -R embl-140/chlo02.uniq.by_taxid.clean.uid.fasta \
       -m 0.75 \
       chlo02/chlo02.obiclean.fasta > chlo02/chlo01.ecotag.fasta

ecotag -d embl-140/ncbi20190930 \
       -R embl-140/euka03.uniq.by_taxid.clean.uid.fasta \
       -m 0.75 \
       euka03/euka03.obiclean.fasta > euka03/euka03.ecotag.fasta

#########
#
# Converts fasta files to a tabular format
#
#

obitab -o -d chlo01/chlo01.ecotag.fasta > chlo01/chlo01.ecotag.tab
obitab -o -d chlo02/chlo02.ecotag.fasta > chlo02/chlo02.ecotag.tab
obitab -o -d euka03/euka03.ecotag.fasta > euka03/euka03.ecotag.tab

