# **SSRGT**：An integrated tool for SSR **genotype calling** **in a hybrid population with resequencing data**
# Introduction
SSRGT is a software used to genotype and SSR calling across a hybrid population with genome resequencing data. It will automatically complete the following five steps.

1. Identification of a large amount of SSR loci in the reference genome.

2. Mapping paired reads to reference genome and extracting high quality comparison results.

3. Genotyping of both parents and progeny SSR genotypes.

4. Filtering SSR loci by the Mendelian segregation ratio and the percent of missing genotypes.

5. Generating input format files for the genetic mapping software JoinMap and FsLinkageMap.

# Usage
To run SSRGT, users should install the two prerequisite packages: [BWA-men2](https://github.com/bwa-mem2/bwa-mem2) and [SAMtools](http://samtools.sourceforge.net/).  Furthermore, an additional setting file parameter is required, namely `parameters.ini`. The parameter file contains three parts, i.e., folders, data files and parameter. As the first part, 'folders' gives the BWA-men2 and SAMtools paths to use to run SSRGT software. In addition, 'SCRIPT' provides a built-in script storage path and, 'PROGENY' provides a path to the progeny resequencing data. The second part 'data files' includes the reference genome file path and the two parents sequencing data files. In addition, ' PROGENY_ID ' provides the first read files and the second read files for all progeny. The third part 'parameter' includes the specified mapping quality (MAPQ) threshold, the number of threads used for parallel computing, the threshold for the frequency of a genotype considered to be a homozygotes genotype or a heterozygous genotype, the threshold for the sum of all reads for an allele, the percent of the maximum missing genotypes at an SSR locus and the minimum p-value allowed for testing the segregation ratio of an SSR locus. A typical parameter file looks as following:

        [folders]  
    	BWA:/mnt/sda/tong/yuyingxuan/software/bwa-mem2/bwa-mem2
    	SAMTOOLS:/mnt/sda/tong/yuyingxuan/software/samtools-1.9/samtools
    	SCRIPT:/mnt/sda/tong/yuyingxuan/zjp/SSRGT/script
    	PROGENY:/mnt/sda/tong/yuyingxuan/zjp/SSRGT
    
        [data files]  
    	REFERENCE_GENOME:/mnt/sda/tong/yuyingxuan/zjp/SSRGT/p.ref
    	PROGENY_ID:/mnt/sda/tong/yuyingxuan/zjp/SSRGT/Progeny.id
    	MALEPARENT_1:/mnt/sda/tong/yuyingxuan/zjp/SSRGT/male.R1.fq
    	MALEPARENT_2:/mnt/sda/tong/yuyingxuan/zjp/SSRGT/male.R2.fq
    	FEMALEPARENT_1:/mnt/sda/tong/yuyingxuan/zjp/SSRGT/female.R1.fq
    	FEMALEPARENT_2:/mnt/sda/tong/yuyingxuan/zjp/SSRGT/female.R2.fq
    
        [parameter]  
    	MAPQ:40
    	THREADS:10
    	Frequency_Of_Homozygotes:0.7
    	Frequency_Of_Heterozygotes:0.3
    	Depth_Of_Coverage:10
    	PVALUE: 0.01
    	MISS_GENOTYPES:0.3


​    In addition, the user is required to provide the progeny information file (PROGENY_ID), each line of which includes the name of the progeny and the first read files and the second read files. A typical PROGENY_ID file looks as following:

```
sample01        sample01.R1.fq  sample01.R2.fq
sample02        sample02.R1.fq  sample02.R2.fq
sample03        sample03.R1.fq  sample03.R2.fq
...... ...... ...... ...... ...... ...... ...... 
sample18        sample18.R1.fq  sample18.R2.fq
sample19        sample19.R1.fq  sample19.R2.fq
sample20        sample20.R1.fq  sample20.R2.fq
```

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
    SSRGT is a software used to genotype and SSR calling across a hybrid population with genome resequencing data.
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

We provide a test data for users to quickly grasp the use of SSRGT. All data files can be downloaded in a compressed file as SSRGT_fq.tar.gz. With the parameter file named `parameters.ini` is configured, we can perform the analysis process by inputting the command: `python SSRGT.py`. When the computing finishes, we can obtain an excel file including all SSR genotype.  In addition, we can obtain 7 segregation types of SSR genotype files, including aaxab.txt,  abxaa.txt, abxab.txt, abxcd.txt, abxac.txt,  aaxbc.txt,and abxcc.txt, where the first two letters of each type denote the maternal genotype and the last two the paternal genotype. Finally, SSRGT can generate an input format for the construction of genetic linkage mapping software.(e.g., [JoinMap](https://www.kyazma.nl/index.php/JoinMap/), [FsLinkageMap](https://link.springer.com/article/10.1007/s11295-010-0281-2) ). 






