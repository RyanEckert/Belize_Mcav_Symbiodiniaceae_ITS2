Symbiodiniaceae *ITS2* amplicon sequencing
==========================================

### Ryan Eckert -- <ryan.j.eckert@gmail.com>

### version: February 13, 2020

------------------------------------------------------------------------

This repository contains scripts and data associated with the
publication **Eckert RJ, Reaume A, Sturm AB, Studivan MS, and Voss
JD (in review) Symbiodiniaceae associations among *Montastraea
cavernosa* populations along a depth gradient on the Belize Barrier
Reef. Front Micro**

Protocols adapted from [Meiog et
al. 2009](https://doi.org/10.1111/j.1755-0998.2008.02222.x); [Klepac et
al. 2013](https://doi.org/10.3354/meps11369).\
<br><br>

------------------------------------------------------------------------

#### Protocols and walkthroughs accompanying this manuscript:

1.  [Protocol for DNA extraction and ITS2
    amplification](https://ryaneckert.github.io/Belize_Mcav_Symbiodiniaceae_ITS2/lab_protocol/)
2.  [Statistical analysis of sequencing
    reads](https://ryaneckert.github.io/Belize_Mcav_Symbiodiniaceae_ITS2/stats/)

------------------------------------------------------------------------

#### Repsitory contents:

- lab_protocol
    *barcodeMM.csv* -- Barcoding PCR mastermix recipe
    *bcPCR.csv* -- Barcoding PCR cycling profile
    *ctab.csv* -- CTAB extraction buffer recipe
    *etbrGel.csv* -- Ethidium bromide gel recipe
    *index.html* -- Symbiodiniaceae lab prorocol webpage
    *its2MM.csv* -- Symbiodiniaceae *ITS2* mastermix recipe
    *its2PCR.csv* -- Symbiodiniaceae *ITS2* PCR cycling profile
    *its2Primers.csv* -- Symbiodiniaceae specific *ITS2* primers
    *reagents.csv* -- Reagent stock recipes
    *supplies.csv* -- Supplies for CTAB extraction
    *sybrGel.csv* -- SYBR gel recipe

- scripts
    *sampleRename.py* -- rename files based on a .csv table
    
- stats
    *62_20190310_DBV_2019-03-11_01-11-25.167036.profiles.absolute.clean.txt* -- *SymPortal* *ITS2* type profile absolute abundance output file cleaned for importing into R
    *62_20190310_DBV_2019-03-11_01-11-25.167036.profiles.absolute.txt* -- *SymPortal* *ITS2* type profile absolute abundance output file
    *62_20190310_DBV_2019-03-11_01-11-25.167036.profiles.relative.txt* -- *SymPortal* *ITS2* type profile relative abundance output file
    *62_20190310_DBV_2019-03-11_01-11-25.167036.seqs.absolute.clean.txt* -- *SymPortal* *ITS2* sequence absolute abundance output file cleaned for importing into R
    *62_20190310_DBV_2019-03-11_01-11-25.167036.seqs.absolute.txt* -- *SymPortal* *ITS2* sequence absolute abundance output file
    *62_20190310_DBV_2019-03-11_01-11-25.167036.seqs.relative.txt* -- *SymPortal* *ITS2* sequence relative abundance output file
    *CBC_MCAV_sampling_metadata.txt* -- *M. cavernosa* sample metadata
    *Simbiodiniaceae_ITS2_statistical_analyses.Rmd* -- Symbiodiniaceae statistical analysis Rmarkdown document
    *index.html* -- Symbiodiniaceae statistical analysis webpage
    *its2Primer.pwk* -- PRIMER7 workbook
    
