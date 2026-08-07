"""option groups"""

# set up a list of option groups to keep things organized
# population comparison options
POP_OPTIONS = [
    "cont_id",
    "ref_files",
    "sex",
    "states",
    "het_steates",
    "homo_states",
    "state_file",
    "random_read_samples" "ancestral",
]


# generate target_file from vcf, bam or geno
INFILE_OPTIONS = [
    "bamfile",
    "deam_cutoff",
    "minmapq",
    "min_length",
    "force_target_file",
    "length_bin_size",
    "report_alleles",
    "tstv",
    "alleles",
    "vcfgt",
    "target",
    "target_file",
]

# generate reference from vcf or geno
REFFILE_OPTIONS = [
    "vcf_ref",
    "pop_file",
    "rec_file",
    "rec_rate",
    "pos_id",
    "map_id",
    "default_map",
    "chroms",
    "force_ref",
    "haplo_chroms",
]


ALGORITHM_OPTIONS_SLUG = [
    "autosomes_only",
    "downsample",
    "fake_contamination",
    "ll_tol",
    "max_iter",
    "split_lib",
    "jk_resamples" "deam_bin_size",
    "len_bin_size",
    "bin_reads",
    "gt_mode" "rrgt",
]


def add_pop_options(parser, states_only=False):
    parser.add_argument(
        "--states",
        "--state-ids",
        nargs="*",
        default=["AFR", "NEA", "DEN"],
        help="""the allowed sources. The target will be made of a mix of all homozygous
        and heterozygous combinations of states. More than 4 or 5 sources have not been
        tested and are not recommended. Must be present in the ref file
        """,
    )
    parser.add_argument(
        "--het-states",
        nargs="*",
        default=None,
        help="""Exact het states to be given. If missing or empty, will use all possible het states
        """,
    )
    parser.add_argument(
        "--homo-states",
        nargs="*",
        default=None,
        help="""Which homozygous states to include. If missing or empty, use all homozygous states
        """,
    )
    parser.add_argument(
        "--state-file",
        "--pop-file",
        default=None,
        help="""Population assignments (yaml format). Doesn't currently support het/homo states""",
    )
    parser.add_argument(
        "--random-read-samples",
        "--pseudo-haploid",
        nargs="*",
        default=[],
        help="""Set a sample as a pseudo-haploid random-read sample for the reference. This means when creating a reference,
        only one allele is taken.
        """,
    )
    if not states_only:
        parser.add_argument(
            "--cont-id",
            "--cont",
            default="AFR",
            help="""the source of contamination. Must be specified in ref file""",
        )
        parser.add_argument(
            "--ancestral",
            "-a",
            type=str,
            default=None,
            help="""Outgroup population with the ancestral allele. By default, assume
            ancestral allele is unknown
            """,
        )


def add_output_options_slug(parser):
    g = parser.add_argument_group(
        "output name and files to be generated",
        """By default, all files are generated. However, if any of 
                                  the --no-* options are used to disable specific files
                                  """,
    )
    g.add_argument(
        "--outname",
        "--out",
        "-o",
        default="admixfrog",
        help="""Output file path (without extensions)""",
    )
    g.add_argument(
        "--no-snp",
        action="store_false",
        default=True,
        dest="output_snp",
        help="""Disable writing posterior genotype likelihood to file with extension .snp.xz""",
    )

    g.add_argument(
        "--no-cont",
        action="store_false",
        default=True,
        dest="output_cont",
        help="""Disable writing contamination estimates to file with extension .bin.xz""",
    )
    g.add_argument(
        "--no-pars",
        action="store_false",
        default=True,
        dest="output_pars",
        help="""Disable writing parameters 
            to file with extension .pars.yaml""",
    )

    parser.add_argument(
        "--no-sfs",
        action="store_false",
        default=True,
        dest="output_sfs",
        help="""Disable output of sfs""",
    )

    parser.add_argument(
        "--output-vcf",
        action="store_true",
        default=False,
        help="""Enable output of vcf""",
    )

    parser.add_argument(
        "--output-jk-sfs",
        default=False,
        action="store_true",
        help="""write a SFS file for each JK resample""",
    )

    parser.add_argument(
        "--output-fstats",
        default=False,
        action="store_true",
        help="""write all F-stats involving target""",
    )

    # parser.add_argument(
    #    "--output-pi",
    #    default=False,
    #    action="store_true",
    #    help="""write all pw differences"""
    # )


def add_target_file_options(parser):
    g = parser.add_argument_group("bam parsing")
    g.add_argument(
        "--bamfile",
        "--bam",
        help="""Bam File to process. Choose this or target_file. 
                   The resulting input file will be writen in {out}.in.xz,
                   so it doesn't need to be regenerated.
                   If the input file exists, an error is generated unless
                   --force-target-file is set
                   """,
    )
    g.add_argument(
        "--force-target-file",
        "--force-bam",
        "--force-infile",
        default=False,
        action="store_true",
    )
    # g.add_argument("--bedfile", "--bed",
    #               help="Bed file with anc/der allele to restrict to")
    g.add_argument(
        "--deam-cutoff",
        type=int,
        default=3,
        help="""reads with deamination in positions < deam-cutoff are
                   considered separately""",
    )
    g.add_argument(
        "--minmapq",
        type=int,
        default=25,
        help="""reads with mapq < MINMAPQ are removed""",
    )
    g.add_argument(
        "--min-length",
        type=int,
        default=35,
        help="""reads with length < MIN_LENGTH are removed""",
    )
    g.add_argument(
        "--length-bin-size",
        type=int,
        default=None,
        help="""if set, reads are binned by length for contamination estimation""",
    )
    g.add_argument(
        "--report-alleles",
        dest="report_alleles",
        action="store_true",
        default=False,
        help="""whether contamination/error rates should be conditioned on alleles present at locus""",
    )
    parser.add_argument(
        "--vcfgt",
        "--vcf-gt",
        "--vcf-target_file",
        help="""VCF input file. To generate input format for admixfrog in genotype mode, use this.
        """,
    )
    parser.add_argument(
        "--target",
        "--name",
        "--sample-id",
        default=None,
        help="""sample name if target is read from vcf or geno file. 
        written in output of f-stats
        """,
    )


def add_ref_options(parser):
    g = parser.add_argument_group("creating reference file")
    g.add_argument(
        "--vcf-ref",
        "--vcf",
        help="""VCF File to process. Choose this or reffile. 
                   The resulting ref file will be writen as {out}.ref.xz,
                   so it doesn't need to be regenerated.
                   If the input file exists, an error is generated unless
                   --force-ref is set
                   """,
    )
    g.add_argument(
        "--rec-file",
        "--rec",
        help="""Recombination rate file. Modelled after 
        https://www.well.ox.ac.uk/~anjali/AAmap/
        If file is split by chromosome, use {CHROM} as 
        wildcards where the chromosome id will be included
        """,
        default=None,
    )
    g.add_argument(
        "--rec-rate",
        default=1e-8,
        type=float,
        help="""Constant recombination rate (per generation per base-pair)""",
    )
    g.add_argument(
        "--pos-id",
        default="Physical_Pos",
        help="""column name for position (default: Physical_Pos)
        """,
    )
    g.add_argument(
        "--map-id",
        default=["AA_Map", "deCODE", "YRI_LD", "CEU_LD"],
        nargs="*",
        help="""column name for genetic map (default: AA_Map)
        """,
    )
    parser.add_argument(
        "--default-map", default="AA_Map", help="""default recombination map column"""
    )
    parser.add_argument(
        "--chroms",
        "--chromosome-files",
        default="1-22,X",
        help="""The chromosomes to be used in vcf-mode.
        """,
    )
    parser.add_argument(
        "--haplo-chroms",
        "--haploid-chroms",
        default=None,
        help="""The chromosomes to be used as haploid. If not set, the following rules apply:
         - chromsomes starting wth one of X, x, Y, y are haploid for males
         - chromsomes starting wth one of Z, z, W, w are haploid for females
         - chromosomes starting with the string "hap" are haploid
         - everything else is diploid
        """,
    )
    parser.add_argument(
        "--sex-chroms",
        default="X,Y,Z,W",
        help="""The chromosomes to be used as sex chromosomes. If not set, 
         - chromsomes starting wth any of [XYZW] are sex chromosomes
        """,
    )
    g.add_argument("--force-ref", "--force-vcf", default=False, action="store_true")


def add_estimation_options_slug(P):
    parser = P.add_argument_group("""options that control estimation of model
                                  parameters""")
    parser.add_argument(
        "--dont-est-contamination",
        action="store_false",
        dest="est_contamination",
        default=True,
        help="""Don't estimate contamination (default do)""",
    )
    parser.add_argument(
        "--dont-est-error",
        action="store_false",
        default=True,
        dest="est_error",
        help="""estimate sequencing error per rg""",
    )
    parser.add_argument(
        "--est-bias",
        action="store_true",
        default=False,
        dest="est_bias",
        help="""estimate reference bias independent from error""",
    )
    parser.add_argument(
        "--dont-est-F",
        action="store_false",
        dest="est_F",
        default=True,
        help="""Estimate F (prob of coalescence in SFS category)""",
    )
    parser.add_argument(
        "--dont-est-tau",
        "-tau",
        action="store_false",
        dest="est_tau",
        default=True,
        help="""Don't estimate tau (prop derived alleles in cSFS)""",
    )
    parser.add_argument(
        "--F0",
        nargs="*",
        type=float,
        default=0.5,
        help="initial F (should be in [0;1]) (default 0)",
    )
    parser.add_argument(
        "--tau0",
        nargs="*",
        type=float,
        default=0.8,
        # help="initial tau (should be in [0;1]) (default 1), at most 1 per source",
        help="initial log-tau (default 0), at most 1 per source",
    )
    parser.add_argument(
        "--e0", "-e", type=float, default=1e-2, help="initial error rate"
    )
    parser.add_argument(
        "--b0", "-b", type=float, default=1e-2, help="initial ref bias rate"
    )
    parser.add_argument(
        "--c0", "-c", type=float, default=1e-2, help="initial contamination rate"
    )


def add_base_options_slug(P):
    parser = P.add_argument_group("options that control the algorithm behavior")
    parser.add_argument(
        "--gt-mode",
        "--gt",
        help="""Assume genotypes are known.
        """,
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--rrgt",
        "--random-read-gt",
        help="""Assume target has known genotypes and random read samples
        """,
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--max-iter",
        "-m",
        type=int,
        default=1000,
        help="""maximum number of iterations""",
    )
    parser.add_argument(
        "--ll-tol", type=float, default=1e-2, help="""stop EM when DeltaLL < ll-tol"""
    )
    parser.add_argument(
        "--ptol",
        type=float,
        default=1e-4,
        help="""stop EM when parameters change by less than ptol""",
    )
    parser.add_argument(
        "--dont-split-lib",
        action="store_false",
        dest="split_lib",
        default=True,
        help="""estimate one global contamination parameter (default: one per read
        group)""",
    )
    parser.add_argument(
        "--autosomes-only",
        action="store_true",
        default=False,
        help="Only run autosomes",
    )
    parser.add_argument(
        "--downsample",
        type=float,
        default=1.0,
        help="downsample coverage to a proportion of reads",
    )
    parser.add_argument(
        "--fake-contamination",
        type=float,
        default=0.0,
        help="Adds fake-contamination from the contamination panel",
    )

    parser.add_argument(
        "--deam-bin-size",
        "--deam-bin",
        type=int,
        default=-1,
        help="bin size for deamination",
    )

    parser.add_argument(
        "--len-bin-size",
        "--len-bin",
        type=int,
        default=50000,
        help="bin size for read length",
    )

    parser.add_argument(
        "--jk-resamples",
        "--n-resamples",
        type=int,
        default=0,
        help="number of resamples for Jackknife standard error estimation",
    )
    parser.add_argument(
        "--male",
        dest="sex",
        action="store_const",
        const="m",
        default=None,
        help="Assumes haploid X chromosome. Default is guess from coverage. currently broken",
    )
    parser.add_argument(
        "--female",
        dest="sex",
        action="store_const",
        const="f",
        default=None,
        help="Assumes diploid X chromosome. Default is guess from coverage",
    )
    parser.add_argument(
        "--chroms",
        "--chromosome-files",
        default="1-22,X",
        help="""The chromosomes to be used in vcf-mode.
        """,
    )


def add_filter_options(parser):
    parser.add_argument(
        "--filter-delta",
        type=float,
        help="""only use sites with allele frequency difference bigger than DELTA (default off)""",
    )
    parser.add_argument(
        "--filter-pos",
        type=int,
        help="""greedily prune sites to be at least POS positions apart""",
    )
    parser.add_argument(
        "--filter-map",
        type=float,
        help="""greedily prune sites to be at least MAP recombination distance apart""",
    )
    parser.add_argument(
        "--filter-high-cov",
        "--filter-highcov",
        type=float,
        default=0.001,
        help="""remove SNP with highest coverage (default 0.001, i.e. 0.1%% of SNP are removed)""",
    )
    parser.add_argument(
        "--filter-ancestral",
        action="store_true",
        default=False,
        help="""remove sites with no ancestral allele information""",
    )
