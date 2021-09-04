# **SSRGT**：An integrated tool for SSR **genotype calling** **in a hybrid population with resequencing data**
# Introduction
SSRGT is a software used to genotype and SSR calling across a hybrid population with genome resequencing data. It will automatically complete the following five steps.

1. Identification of a large amount of SSR loci in the reference genome.

2. Mapping paired reads to reference genome and extracting high quality comparison results.

3. Genotyping of both parents and progeny SSR genotypes.

4. Filtering SSR loci by the Mendelian segregation ratio and the percent of missing genotypes.

5. Generating input format files for the genetic mapping software JoinMap and FsLinkageMap.
# Requirements
To compile SSRGT, the following packages are required:
- [BWA-men2](https://github.com/bwa-mem2/bwa-mem2)  or  [BWA](https://github.com/lh3/bwa)
- [SAMtools](http://samtools.sourceforge.net/)
- [pandas](https://github.com/pandas-dev/pandas)
- [xlwt](https://github.com/python-excel/xlwt)
- perl   (version >=5.26)
- python (version >=3.6)
- java   (version >=1.8)
You can do this with the following command：
```
conda install bwa
conda install samtools
pip install pandas
pip install xlwt
```
# Usage
To run SSRGT, users should install two prerequisite packages: BWA and SAMtools.  Furthermore, an additional setting file parameter is required, namely `parameters.ini`. The parameter file contains three parts: folders, parameter and files. As the first part, ‘folders’ gives the software paths of the BWA and SAMtools. In addition, a script storage path and resequencing data path need to be provided. The second part ‘parameter’ consists of the specified mapping quality (MAPQ) threshold, the number of threads used for parallel computing, the depth of allele coverage, the minimum percentage frequency of a homozygotes, the minimum percentage frequency of a heterozygotes, the percent of the maximum missing genotypes at an SSR locus and the minimum p-value allowed for testing the segregation ratio of an SSR locus. Finally, the ‘ files’ includes reference genome file, the names of parents and progeny and the first read files and the second read files.  A typical parameter file looks as following:

        [folders]  
        BWA:~/bwa-0.7.15
        SAMTOOLS:~/software/samtools-1.9
        SSRGT FOLD:~/zjp/SSRGT
        RADDATA FOLD:/mnt/sda/tong/yuyingxuan/zjp/exampledata
        [parameter]
        MAPQ:30
        THREADS:30
        Depth Of Coverage:5
        Frequency Of Homozygotes:0.7
        Frequency Of Heterozygotes:0.35
        PVALUE: 0.01
        MISS GENOTYPES:0.3
        [files]
        reference genome:p.ref
        male:male.R1.fq male.R2.fq
        female:female.R1.fq female.R2.fq
        sample01:sample01.R1.fq sample01.R2.fq
        sample02:sample02.R1.fq sample02.R2.fq
        sample03:sample03.R1.fq sample03.R2.fq
        '''
        sample18:sample18.R1.fq sample18.R2.fq
        sample19:sample19.R1.fq sample19.R2.fq
        sample20:sample20.R1.fq sample20.R2.fq
 

Additionally, users need to install python modules in python 3 environment, including [pandas](https://github.com/pandas-dev/pandas) and [xlwt](https://github.com/python-excel/xlwt).

  When the required software packages are installed and the parameter file is parepared and saved in a work directory, you can go to the work directory and get started with the command:  
 ```
 python SSRGT.py
 
 If you want to change the reference sequence motif threshold for SSR genotype calling, you can run the command in background as :
  
 python SSRGT.py -mo 1=10,2=6,3=5,4=5,5=5,6=4  `
 ```

 Since it will take quite a few hours or even several days to finish a pratical computing, we usually run the command in background as  
```
  nohup python SSRGT &
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
                            set the threshold of motif, default: 1=11,2=5,3=4,4=4,5=4,6=4
                            left  of equal : length of motif
                            right of equal : the minimum number of repeat
      -p POPULATION, --population POPULATION
                            set the population type [CP,F2,BC], default : CP
                            If your population type  is CP or F2, you should choose 'CP' or 'F2', if it's BC population type and the    maternal parent is recurrent parent, choose 'BC:female'
After the program has finished running, if the user wants to locate the SSR loci by reference genome，they can use the following command  in the 'script' folder. This step will facilitate functional annotation (e.g., GO, KEGG, Pfam) studies of genes containing SSRs.

```
usage: python SSRlocation.py [-h] -g GFF -s SSR_TYPES -t THREADS  -o OUT

Example: python SSRlocation.py -g ../p.gff3 -s ../WorkingDirectory/femaleallSSR_type.txt 
-t 2 -o female.location

optional arguments:
  -h, --help            show this help message and exit
  -g GFF, --gff GFF     Please input GFF file
  -s SSR_TYPES, --SSR_types SSR_TYPES
                        Please input SSR_types.txt file
  -t THREADS, --threads THREADS
                        Please input threads
  -o OUT, --out OUT     Please input out_put file_name

```

# Test Data

We provide a test data for users to quickly grasp the use of SSRGT. All data files can be downloaded in a compressed file as [SSRGTfq.tar.gz](https://figshare.com/articles/dataset/SSRGT_fq_tar_gz/15094539). With the parameter file named `parameters.ini` is configured, we can perform the analysis process by inputting the command: `python SSRGT.py`. When the computing finishes, we can obtain an excel file including all SSR genotype.  In addition, we can obtain 7 segregation types of SSR genotype files, including aaxab.txt,  abxaa.txt, abxab.txt, abxcd.txt, abxac.txt,  aaxbc.txt,and abxcc.txt, where the first two letters of each type denote the maternal genotype and the last two the paternal genotype. Finally, SSRGT can generate an input format for the construction of genetic linkage mapping software.(e.g., [JoinMap](https://www.kyazma.nl/index.php/JoinMap/), [FsLinkageMap](https://link.springer.com/article/10.1007/s11295-010-0281-2) ). 






