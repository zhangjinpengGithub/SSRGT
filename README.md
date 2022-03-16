# **SSRGT**：An integrated tool for SSR **genotype calling** **in a hybrid population with resequencing data**
# Introduction
SSRGT is a software used to genotype and SSR calling across a hybrid population with genome resequencing data. It will automatically complete the following five steps.

1. Identification of a large amount of SSR loci in the reference genome.

2. Mapping paired reads to reference genome and extracting high quality comparison results.

3. Genotyping of both parents and progeny SSR genotypes.

4. Filtering SSR loci by the Mendelian segregation ratio and the percent of missing genotypes.


# Requirements
To compile SSRGT, the following packages are required:
- [BWA-men2](https://github.com/bwa-mem2/bwa-mem2)  or  [BWA](https://github.com/lh3/bwa)
- [SAMtools](http://samtools.sourceforge.net/)
- [bcftools](https://github.com/samtools/bcftools)
- [pandas](https://github.com/pandas-dev/pandas)
- scipy
- perl   (version >=5.26)
- python (version >=3.6)
```
You can do this with the following command：
conda install bwa
conda install samtools
pip install pandas
pip install scipy (version =1.3.1)
```
# Usage
To run SSRGT, users should install two prerequisite packages: BWA and SAMtools.  Furthermore, an additional setting file parameter is required, namely `parameters.ini`. The parameter file contains three parts: folders, parameters and files. As the first part, ‘folders’ gives the software paths of the BWA、 SAMtools and BCFtools. In addition, a script storage path and resequencing data path need to be provided. The second part ‘parameter’ consists of the number of threads used for parallel computing, the depth of allele coverage,the percent of the maximum missing genotypes at an SSR locus and the minimum p-value allowed for testing the segregation ratio of an SSR locus, the minimum GQ scores. Finally, the ‘data files’ includes reference genome file, the names of parents and progeny and the first read files and the second read files.  A typical parameter file looks as following:

        [folders]
        BWA_FOLD:~/bwa-0.7.15
        SAMTOOLS_FOLD:~/software/samtools-1.9
        BCFTOOLS_FOLD:/mnt/sda/tong/yuyingxuan/software/bcftools-1.9
        SSRGT_FOLD:/mnt/sda/tong/yuyingxuan/zjp/SSRGT
        RADDATA_FOLD:/mnt/sda/tong/yuyingxuan/zjp/exampledata
        [parameters]
        THREADS:30
        ALLELE_DEPTH:5
        MISS_GENOTYPES:0.3
        PVALUE: 0.01
        GQ:20
        [data files]
        REFERENCE_FILE:reference.fasta
        MALE:male.R1.fq male.R2.fq
        FEMALE:female.R1.fq female.R2.fq
        SAMPLE01:sample01.R1.fq sample01.R2.fq
        SAMPLE02:sample02.R1.fq sample02.R2.fq
        SAMPLE03:sample03.R1.fq sample03.R2.fq
        SAMPLE04:sample04.R1.fq sample04.R2.fq
        SAMPLE05:sample05.R1.fq sample05.R2.fq
        SAMPLE06:sample06.R1.fq sample06.R2.fq
        SAMPLE07:sample07.R1.fq sample07.R2.fq
        SAMPLE08:sample08.R1.fq sample08.R2.fq
        SAMPLE09:sample09.R1.fq sample09.R2.fq
        SAMPLE10:sample10.R1.fq sample10.R2.fq
        SAMPLE11:sample11.R1.fq sample11.R2.fq
        SAMPLE12:sample12.R1.fq sample12.R2.fq
        SAMPLE13:sample13.R1.fq sample13.R2.fq
        SAMPLE14:sample14.R1.fq sample14.R2.fq
        SAMPLE15:sample15.R1.fq sample15.R2.fq
        SAMPLE16:sample16.R1.fq sample16.R2.fq
        SAMPLE17:sample17.R1.fq sample17.R2.fq
        SAMPLE18:sample18.R1.fq sample18.R2.fq
        SAMPLE19:sample19.R1.fq sample19.R2.fq
        SAMPLE20:sample20.R1.fq sample20.R2.fq
 

Additionally, users need to install python modules [pandas](https://github.com/pandas-dev/pandas) in python 3 environment.

  When the required software packages are installed and the parameter file is parepared and saved in a work directory, you can go to the work directory and get started with the command:  
 ```
 python SSRGT.py
 
 If you want to change the reference sequence motif threshold for SSR genotype calling, you can run the command in background as :
  
 python SSRGT.py -mo 1=10,2=6,3=5,4=5,5=5,6=4  `
 ```
 In addition, you can select the population type (CP,F2,BC) with the -p parameter, If your population type  is CP or F2,
 you can run the command in background as :
 ```
   python SSRGT -p CP
   python SSRGT -p F2
```
 If your population type  is BC, and the maternal parent is recurrent parent, you can run the command in background as :
 ```
   python SSRGT -p BC:female
 ```
 If it's BC population type and the paternal parent is recurrent parent, you can run the command in background as :
 ```
   python SSRGT -p BC:male
 ```
 Since it will take quite a few hours or even several days to finish a pratical computing, we usually run the command in background as  
```
  nohup python SSRGT &
```
 It can be seen that users can perform the analytical steps independently by adding some options described as above if some prerequisite files are avaiable. For example, if the user modifies the PVALUE and MISS_GENOTYPES parameters in the parameter file, the following command can be used to perform only the fifth step.
 ```
 python SSRGT.py --nosearch  --nocatalog --nomap --nocall
 ```
 

  You can run with the 'help' option (`pthon SSRGT.py -h`) to show the usage of SSRGT:

    usage: python SSRGT.py [-h] [-mo MOTIF]
     ----------------------------------------------------------------------------------------
    SSRGT provides a usable and accurate SSR genotyping platform for distant hybridization populations.
    It has the advantage that a large number of SSR loci segregated by Mendelian segregation ratio in a
    population can be discovered with a simple command and an input format file can be generated for the
    genetic mapping software. It will facilitate the construction of SSR genetic linkage maps,
    locating quantitative trait loci, marker-assisted selection breeding and various genetic studies.
    Contact: Chunfa Tong  <tongchf@njfu.edu.cn>
    Version: 1.0
    Usage: python SSRGT.py  [Options]
    
    -----------------------------------------------------------------------------------------
    optional arguments:
      -h, --help            show this help message and exit
      -mo MOTIF, --motif MOTIF
                            set the threshold of motif, default: 1=10,2=5,3=4,4=4,5=4,6=4
                            left  of equal : length of motif
                            right of equal : the minimum number of repeat
      -p POPULATION, --population POPULATION
                        set the population type [CP,F2,BC:female,BC:male], default : CP
                        If your population type  is CP or F2, you should choose 'CP' or 'F2', if it's BC population type and the maternal parent is recurrent parent, choose 'BC:female', if it's BC population type and the paternal parent is recurrent parent, choose 'BC:male'
      --nosearch            skip the step for searching SSR in reference sequence.
      --nocatalog           skip the step for generating parental SSR catalogs.
      --nomap               skip the step for mapping the progeny reads to the reference sequence.
      --nocall              skip the step for calling  SSR genotypes for each individual.
      --nofilter            skip the step for filtering  SSR genotype data.

# Test Data

We provide a test data for users to quickly grasp the use of SSRGT. All data files can be downloaded in a compressed file as [SSRGT_ExampleData.tar.gz](https://figshare.com/articles/dataset/SSRGT_fq_tar_gz/15094539). With the parameter file named `parameters.ini` is configured, we can perform the analysis process by inputting the command: `python SSRGT.py`. When the computing finishes, we can obtain 7 segregation types of SSR genotype files, including aaxab.txt,  abxaa.txt, abxab.txt, abxcd.txt, abxac.txt,  aaxbc.txt,and abxcc.txt, where the first two letters of each type denote the maternal genotype and the last two the paternal genotype. Finally, the user can use the above results to easily construct a genetic linkage map.(e.g., [JoinMap](https://www.kyazma.nl/index.php/JoinMap/), [FsLinkageMap](https://link.springer.com/article/10.1007/s11295-010-0281-2) ). 






