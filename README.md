# Admixslug

[![Tests](https://github.com/BenjaminPeter/admixslug/actions/workflows/tests.yaml/badge.svg)](https://github.com/BenjaminPeter/admixslug/actions/workflows/tests.yaml)
![PyPI Wheel](https://img.shields.io/pypi/wheel/admixslug?label=pypi)



Admixslug is a genotype likelihood method for contaminated low-coverage nuclear DNA data from
archaic humans. It works by computing a conditional site-frequency spectrum. Outputs of admixslug include contamination estimates and f-statistics.

## Hardware and Software requirements

No special hardware is required.
The package development version is tested on Linux operating systems. The initial version has been tested on Linux - Ubuntu 24.04.4 LTS machines.

All package and version dependencies are listed in pyproject.toml and poetry.lock, and automatically installed with admixslug.

## Installation
Requires `python3.10+`

You can create a conda environment using the .yml file provided (admixslug_env.yml).

```
conda env create --name admixslug_env --file=admixslug_env.yml
conda activate admixslug_env
```
This should take about 3-4 minutes on a desktop computer.

Install latest release of `admixslug` (from pypi):
```
pip install admixslug
```

Or install the latest development version from github (potentially unstable):
```
pip install git+https://github.com/benjaminpeter/admixslug
```
This should take less than 10 seconds on a desktop computer.

## Input files
admixslug requires two input files;

  - a reference file with information from high-quality samples, including the contaminant panel (e.g. tests/data/ref_bigsteffi.csv.xz)
  - a sample file that stores read information for a sample in compact format (e.g. tests/data/FQ.L30MQ25_bigsteffi52.bam, which contains nuclear DNA data from the Forbes' Quarry Neandertal from Gibraltar).


```
admixslug-bam --bam tests/data/FQ.L30MQ25_bigsteffi52.bam --out tests/data/FQ_bigsteffi52.in.xz  --ref tests/data/ref_bigsteffi.csv.xz
```
The reference file is created exactly the same way as in admixfrog (https://github.com/benjaminpeter/admixfrog). The bamfile
contains the reads to be analyzed, and the `--out` flag designates where the
input file will be stored. see `admixslug-bam --help` for details.

On a desktop computer, this step takes less than a minute.

## Quickstart

The following command runs admixslug on a single sample stored in
`tests/data/FQ_bigsteffi52.in.xz` using the sites from `tests/data/ref_bigsteffi.csv.xz`
and saving the output files in `outputs/FQ_bigsteffi52`

```
admixslug --infile tests/data/FQ_bigsteffi52.in.xz \
        --ref tests/data/ref_bigsteffi.csv.xz \
        -o outputs/FQ_bigsteffi52  \
        --states ALT VIN CHA DEN  \
        --cont-id EUR  \
        --ancestral PAN \
        --ll-tol 0.01  \
        --ptol 0.001   \
        --max-iter 100 \
        --filter-pos 50 \
        --filter-ancestral  \
        --len-bin-size 2000\
        --jk-resamples 10\
        --output-jk-sfs\
        --output-fstats
```

On a desktop computer, this should take a little over 1 minute.

The remaining arguments are

  - `--states` : The reference samples or populations to condition the SFS on
  - `--cont-id` : The putative contaminant panel
  - `--ancestral` The ancestral state (these three need to be defined in the
    reference file)
  - `--ll-tol, -ptol`: Convergence criteria in terms of log-likelihood and
    changes in parameter values, respectively
  - `--max-iter` : The maximum number of iterations
  - `--filter-pos` : filter position to be at least x bases apart
  - `--filter-ancestral` : only retain sites with ancestral allele info
  - `--len-bin k `: attempt to bin reads into bins with around k sites. Higher
    numbers of k will result in fewer length-bins for contamination estimation,
    and lower numbers will result in many uncertain estimates
  - `--jk-resamples` : the number of jackknife resamples for standard error
    estimation
  - `--output-jk-sfs`: write a SFS file for each JK resample
  - `--output-fstats`: write all F-stats involving target


## Output
The main outputs are

#### Contamination file
This file, named {out}.cont.xz contains contamination information for different sequence groups
and different categories (sequences with no deamination, sequences with terminal
deamination etc.).

```
rg            cont      n_exact  n_sites  se_cont   l_cont    h_cont
library1_255_0  0.972531  53135    669098   0.006565  0.959664  0.985398
library2_255_0  0.889277  36158    669098   0.009955  0.869765  0.908789
library3_255_0  0.947004  76550    669098   0.009063  0.929241  0.964766
```

#### SFS file
This file, named {out}.sfs.xz contains info on the estimated SFS

```
sex_chrom  VIN_anc  VIN_der  CHA_anc  CHA_der  ALT_anc  ALT_der  DEN_anc  DEN_der  PAN_anc  PAN_der  F         tau       n_snps  n_reads  n_endo        read_ratio  cont_est  psi        se_tau    l_tau     h_tau     se_F      l_F       h_F
autosome   0        2        0        2        0        2        0        2        1        0        0.289165  0.970069  47982   80951    9278.798637   0.617151    0.885378  0.571462   0.024662  0.921731  1.000000  0.074776  0.142604  0.435727
autosome   0        2        0        2        0        2        2        0        1        0        0.232955  0.893832  35125   57736    6577.306616   0.344759    0.886080  0.274166   0.068626  0.759325  1.000000  0.080282  0.075602  0.390308
autosome   2        0        2        0        2        0        0        2        1        0        0.367088  0.010731  42742   69756    7865.868251   0.169577    0.887237  0.189765   0.010181  0.000000  0.030685  0.239400  0.000000  0.836312
autosome   1        1        0        2        1        1        2        0        1        0        0.423531  0.519948  665     1068     128.319719    0.183521    0.879850  0.137579   0.093208  0.337260  0.702636  0.237326  0.000000  0.888690
```

Admixslug outputs another file named {out}.jksfs.xz which contains the same
information but for each JK resample.

#### vcf-file
This file, named {out}.vcf contains a vcf file with i) random read samples, ii)
genotype likelihoods and iii) genotype probabilities for all sites with coverage.

#### snp-file
This file, named {out}.snp.xz contains similar info as the VCF file, but more
easily readable in R.

```
chrom  pos        map         ref  alt  tref  talt  G0          G1          G2          p         random_read  sfs
1      834832     0.000000    G    C    1     1     -2.209128   -1.459541   -0.018131   0.976466  1            0
1      839495     0.000000    G    T    1     0     -1.454317   -0.826969   -0.088351   0.890396  1            1
1      846864     0.000000    G    C    3     0     -0.004820   -2.028945   -2.774279   0.006359  0            2
1      851204     0.000000    G    C    4     1     -1.674939   -0.910296   -0.067567   0.917391  1            1
1      853267     0.000000    G    T    1     0     -2.084439   -1.408942   -0.021013   0.972267  1            0
```

#### pi-file

This file, named {out}.pi.xz contains the pairwise differences used for the
calculation of the f-statistics.
```
sex_chrom  pop1              pop2              is_between  pi        sd        sterr
autosome   ALT               ALT               False       0.046972  0.002154  0.046409
autosome   ALT               CHA               True        0.079208  0.001168  0.034179
autosome   ALT               DEN               True        0.308555  0.001155  0.033984
autosome   ALT               PAN               True        0.305752  0.001218  0.034899
```
#### f-files
Currently there are seperate output files for f2, f3 and f4 statistics, names
{out}.f2.xz, {out}.f3.xz and {out}.f4.xz.
These files contain the names of the individuals the statistics is calculated
for, type of data (autosomal vs. sex chromosome), value of the statistics and
the uncertainity.

Admixslug outputs additional files named {out}.f2.jk.xz, {out}.f3.jk.xz and
{out}.f4.jk.xz. These files contain the names of the individuals the statistics
is calculated for, type of data (autosomal vs. sex chromosome) and the value of
the statistic for each JK resample.

## Contact
Benjamin Peter [benjamin_peter@eva.mpg.de](benjamin_peter@eva.mpg.de), 
Arev Sümer [arev_suemer@eva.mpg.de](arev_suemer@eva.mpg.de)

## License
This project is covered under the BSD 3-Clause License.


## Documentation
Full description of the algorithm is available in the supplementary information of our manuscript which will be soon on bioarxiv.
Some less up to date files can be found currently in [docs/admixslug.pdf](docs/admixslug.pdf).
Changes are that `admixslug --help` will give more up-to-date information.

Full command and parameters can be found by typing `admixslug --help`:
```
usage: admixslug [-h] [-v] [--target-file TARGET_FILE] [--ref REF_FILES] [--seed SEED]                
                 [--sex-chroms SEX_CHROMS] [--bamfile BAMFILE] [--force-target-file]                  
                 [--deam-cutoff DEAM_CUTOFF] [--minmapq MINMAPQ] [--min-length MIN_LENGTH]            
                 [--length-bin-size LENGTH_BIN_SIZE] [--report-alleles] [--vcfgt VCFGT]               
                 [--target TARGET] [--dont-est-contamination] [--dont-est-error] [--est-bias]         
                 [--dont-est-F] [--est-tau] [--F0 [F0 ...]] [--tau0 [TAU0 ...]] [--e0 E0] [--b0 B0]   
                 [--c0 C0] [--max-iter MAX_ITER] [--ll-tol LL_TOL] [--ptol PTOL] [--dont-split-lib]   
                 [--autosomes-only] [--downsample DOWNSAMPLE]                                         
                 [--fake-contamination FAKE_CONTAMINATION] [--deam-bin-size DEAM_BIN_SIZE]            
                 [--len-bin-size LEN_BIN_SIZE] [--jk-resamples JK_RESAMPLES] [--male] [--female]      
                 [--chroms CHROMS] [--outname OUTNAME] [--no-snp] [--no-cont] [--no-pars]             
                 [--no-sfs] [--output-vcf] [--output-jk-sfs] [--output-fstats]                        
                 [--states [STATES ...]] [--het-states [HET_STATES ...]]                              
                 [--homo-states [HOMO_STATES ...]] [--state-file STATE_FILE]                          
                 [--random-read-samples [RANDOM_READ_SAMPLES ...]] [--cont-id CONT_ID]                
                 [--ancestral ANCESTRAL] [--filter-delta FILTER_DELTA] [--filter-pos FILTER_POS]      
                 [--filter-map FILTER_MAP] [--filter-high-cov FILTER_HIGH_COV] [--filter-ancestral]   

Infer sfs and contamination from low-coverage and contaminated genomes                                

options:                                                                                              
  -h, --help            show this help message and exit                                               
  -v, --version         show program's version number and exit                                        
  --target-file TARGET_FILE, --infile TARGET_FILE, --in TARGET_FILE                                   
                        Sample input file (csv). Contains individual specific data, obtained from a   
                        bam file. - Fields are chrom, pos, map, lib, tref, talt" - chrom:             
                        chromosome - pos : physical position (int) - map : rec position (float) -     
                        lib : read group. Any string, same string assumes same contamination - tref   
                        : number of reference reads observed - talt: number of alt reads observed     
  --ref REF_FILES, --ref-file REF_FILES                                                               
                        refernce input file (csv). - Fields are chrom, pos, ref, alt, map, X_alt,     
                        X_ref - chrom: chromosome - pos : physical position (int) - ref : refrence    
                        allele - alt : alternative allele - map : rec position (float) - X_alt,       
                        X_ref : alt/ref alleles from any number of sources / contaminant              
                        populations. these are used later in --cont-id and --state-id flags           
  --seed SEED           random number generator seed for resampling                                   
  --sex-chroms SEX_CHROMS                                                                             
                        The chromosomes to be used as sex chromosomes. If not set, - chromsomes       
                        starting wth any of [XYZW] are sex chromosomes                                
  --vcfgt VCFGT, --vcf-gt VCFGT, --vcf-target_file VCFGT                                              
                        VCF input file. To generate input format for admixfrog in genotype mode,      
                        use this.                                                                     
  --target TARGET, --name TARGET, --sample-id TARGET                                                  
                        sample name if target is read from vcf or geno file. written in output of     
                        f-stats                                                                       
  --no-sfs              Disable output of sfs                                                         
  --output-vcf          Enable output of vcf                                                          
  --output-jk-sfs       write a SFS file for each JK resample                                         
  --output-fstats       write all F-stats involving target                                            
  --states [STATES ...], --state-ids [STATES ...]                                                     
                        the allowed sources. The target will be made of a mix of all homozygous and   
                        heterozygous combinations of states. More than 4 or 5 sources have not been   
                        tested and are not recommended. Must be present in the ref file               
  --het-states [HET_STATES ...]                                                                       
                        Exact het states to be given. If missing or empty, will use all possible      
                        het states                                                                    
  --homo-states [HOMO_STATES ...]                                                                     
                        Which homozygous states to include. If missing or empty, use all homozygous   
                        states                                                                        
  --state-file STATE_FILE, --pop-file STATE_FILE                                                      
                        Population assignments (yaml format). Doesn't currently support het/homo      
                        states                                                                        
  --random-read-samples [RANDOM_READ_SAMPLES ...], --pseudo-haploid [RANDOM_READ_SAMPLES ...]         
                        Set a sample as a pseudo-haploid random-read sample for the reference. This   
                        means when creating a reference, only one allele is taken.                    
  --cont-id CONT_ID, --cont CONT_ID                                                                   
                        the source of contamination. Must be specified in ref file                    
  --ancestral ANCESTRAL, -a ANCESTRAL                                                                 
                        Outgroup population with the ancestral allele. By default, assume ancestral   
                        allele is unknown                                                             
  --filter-delta FILTER_DELTA                                                                         
                        only use sites with allele frequency difference bigger than DELTA (default    
                        off)                                                                          
  --filter-pos FILTER_POS                                                                             
                        greedily prune sites to be at least POS positions apart                       
  --filter-map FILTER_MAP                                                                             
                        greedily prune sites to be at least MAP recombination distance apart          
  --filter-high-cov FILTER_HIGH_COV, --filter-highcov FILTER_HIGH_COV                                 
                        remove SNP with highest coverage (default 0.001, i.e. 0.1% of SNP are         
                        removed)                                                                      
  --filter-ancestral    remove sites with no ancestral allele information                             

bam parsing:                                                                                              
  --bamfile BAMFILE, --bam BAMFILE                                                                        
                        Bam File to process. Choose this or target_file. The resulting input file         
                        will be writen in {out}.in.xz, so it doesn't need to be regenerated. If the       
                        input file exists, an error is generated unless --force-target-file is set        
  --force-target-file, --force-bam, --force-infile                                                        
  --deam-cutoff DEAM_CUTOFF                                                                               
                        reads with deamination in positions < deam-cutoff are considered separately       
  --minmapq MINMAPQ     reads with mapq < MINMAPQ are removed                                             
  --min-length MIN_LENGTH                                                                                 
                        reads with length < MIN_LENGTH are removed                                        
  --length-bin-size LENGTH_BIN_SIZE                                                                       
                        if set, reads are binned by length for contamination estimation                   
  --report-alleles      whether contamination/error rates should be conditioned on alleles present        
                        at locus                                                                          

options that control estimation of model                                                                  
                                  parameters:                                                             
  --dont-est-contamination                                                                                
                        Don't estimate contamination (default do)                                         
  --dont-est-error      estimate sequencing error per rg                                                  
  --est-bias            estimate reference bias independent from error                                    
  --dont-est-F          Estimate F (distance from ref, default False)                                     
  --est-tau, -tau       Estimate tau (population structure in references)                                 
  --F0 [F0 ...]         initial F (should be in [0;1]) (default 0)                                        
  --tau0 [TAU0 ...]     initial log-tau (default 0), at most 1 per source                                 
  --e0 E0, -e E0        initial error rate                                                                
  --b0 B0, -b B0        initial ref bias rate                                                             
  --c0 C0, -c C0        initial contamination rate                                                        

options that control the algorithm behavior:                                                              
  --max-iter MAX_ITER, -m MAX_ITER                                                                        
                        maximum number of iterations                                                      
  --ll-tol LL_TOL       stop EM when DeltaLL < ll-tol                                                     
  --ptol PTOL           stop EM when parameters change by less than ptol                                  
  --dont-split-lib      estimate one global contamination parameter (default: one per read group)         
  --autosomes-only      Only run autosomes                                                                
  --downsample DOWNSAMPLE                                                                                 
                        downsample coverage to a proportion of reads                                      
  --fake-contamination FAKE_CONTAMINATION                                                                 
                        Adds fake-contamination from the contamination panel                              
  --deam-bin-size DEAM_BIN_SIZE, --deam-bin DEAM_BIN_SIZE                                                 
                        bin size for deamination                                                          
  --len-bin-size LEN_BIN_SIZE, --len-bin LEN_BIN_SIZE                                                     
                        bin size for read length                                                          
  --jk-resamples JK_RESAMPLES, --n-resamples JK_RESAMPLES                                                 
                        number of resamples for Jackknife standard error estimation                       
  --male                Assumes haploid X chromosome. Default is guess from coverage. currently           
                        broken                                                                            
  --female              Assumes diploid X chromosome. Default is guess from coverage                      
  --chroms CHROMS, --chromosome-files CHROMS                                                              
                        The chromosomes to be used in vcf-mode.                                           

output name and files to be generated:                                                                    
  By default, all files are generated. However, if any of the --no-* options are used to disable          
  specific files                                                                                          

  --outname OUTNAME, --out OUTNAME, -o OUTNAME                                                            
                        Output file path (without extensions)                                             
  --no-snp              Disable writing posterior genotype likelihood to file with extension              
                        .snp.xz                                                                           
  --no-cont             Disable writing contamination estimates to file with extension .bin.xz            
  --no-pars             Disable writing parameters to file with extension .pars.yaml                      
```
