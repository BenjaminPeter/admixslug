import logging
from collections import defaultdict, namedtuple

import numpy as np
import pandas as pd
from numba import njit
from scipy.linalg import expm

Probs = namedtuple(
    "Probs", ("O", "N", "P_cont", "alpha", "beta", "rg", "alpha_hap", "beta_hap", "S")
)
Pars = namedtuple(
    "Pars",
    ("alpha0", "alpha0_hap", "trans", "trans_hap", "cont", "error", "F", "tau", "sex"),
)  # for returning from varios functions


def get_haploid_stuff(snp, chroms, sex):
    haplo_chroms, diplo_chroms = [], []

    if sex is None:
        diplo_chroms = chroms
    else:
        for c in chroms:
            if c[0] in "YyWw":
                haplo_chroms.append(c)
            elif c[0] in "Xx" and sex == "m":
                haplo_chroms.append(c)
            elif c[0] in "Zz" and sex == "f":
                haplo_chroms.append(c)
            elif c.startswith("hap"):
                haplo_chroms.append(c)
            else:
                diplo_chroms.append(c)

    haploid_snps = snp.snp_id[snp.chrom.isin(haplo_chroms)]
    if len(haploid_snps) > 0:
        haploid_snps = slice(min(haploid_snps), max(haploid_snps) + 1)
    else:
        haploid_snps = slice(0, 0)

    return haplo_chroms, haploid_snps


def make_flipped(snp, anc_ref, anc_alt):
    FLIPPED = (snp[anc_ref] == 0) & (snp[anc_alt] > 0)
    return FLIPPED.to_numpy()


def polarize(snp, pop, flipped):
    sfs = pd.DataFrame()
    sfs[f"{pop}_anc"] = snp[f"{pop}_ref"]
    sfs[f"{pop}_der"] = snp[f"{pop}_alt"]
    sfs.loc[flipped, f"{pop}_anc"] = snp.loc[flipped, f"{pop}_alt"]
    sfs.loc[flipped, f"{pop}_der"] = snp.loc[flipped, f"{pop}_ref"]
    return sfs


def obs2sfs(snp, flipped, states, sex_chroms=["Z", "W", "X", "Y"]):
    """create sfs data structure taking ancestral allele into account

    basic strat
    1. make dict[state] : index for all possible indices
    2. use dict to create SNP2SFS


    """

    snp.reset_index(drop=True, inplace=True)
    sfs = pd.DataFrame()
    if sex_chroms is not None:
        sfs["sex_chrom"] = np.where(snp["chrom"].isin(sex_chroms), "sex", "autosome")

    """polarize all input data"""
    for s in states:
        pol1 = polarize(snp, s, flipped)
        sfs = pd.concat((sfs, pol1), axis=1)

    sfs_rows = sfs.drop_duplicates().reset_index(drop=True)
    sfs_dict = dict(
        (tuple(v.values()), k) for (k, v) in sfs_rows.to_dict("index").items()
    )
    """use dicts to create SNP2SFS"""
    SNP2SFS = np.array([sfs_dict[tuple(i)] for i in sfs.to_numpy()], dtype=np.uint16)

    return sfs_rows, SNP2SFS


@njit
def make_full_read_df(df, n_reads):
    """generate vectors relating read groups to allele they carry, read group and snp

    returns vectors of length n_reads
    """
    READS = np.empty(n_reads, np.uint8)
    READ2RG = np.empty(n_reads, np.uint32)
    READ2SNP = np.empty(n_reads, np.uint32)

    i = 0
    for snp_id, ref, alt, rg in df:
        for _ in range(ref):
            READS[i] = 0
            READ2RG[i] = rg
            READ2SNP[i] = snp_id
            i += 1
        for _ in range(alt):
            READS[i] = 1
            READ2RG[i] = rg
            READ2SNP[i] = snp_id
            i += 1

    return READS, READ2RG, READ2SNP


def init_ftau(n_states, F0=0.5, tau0=0):
    """initializes F and tau, which exist for each homozygous state"""
    try:
        if len(F0) == n_states:
            F = F0
        elif len(F0) == 1:
            F = F0 * n_states
        else:
            F = [F0]
    except TypeError:
        F = [F0] * n_states
    try:
        if len(tau0) == n_states:
            tau = tau0
        elif len(F0) == 1:
            tau = tau0 * n_states
        else:
            tau = [tau0]
    except TypeError:
        tau = [tau0] * n_states

    return np.array(F), np.array(tau)


def init_ce(c0=0.01, e0=0.001):
    cont = defaultdict(lambda: c0)
    error = defaultdict(lambda: e0)
    return cont, error


def trans_mat_hap_to_dip(tmat):
    """given a haploid transition rate matrix tmat, returns a diploid transition
    matrix

    assumptions:
     - only one transition at a time
     - "canonical" order of heterozygous states
     - independence between haplotypes
    """
    n = tmat.shape[0]
    n_homo = n
    n_het = int(n * (n - 1) / 2)
    tmat2 = np.zeros((n_homo + n_het, n_homo + n_het))

    # homo -> homo transition
    # for i in range(n_homo):
    #    for j in range(n_homo):
    #        tmat2[i, j] = tmat[i, j] ** 2 #prob both change at once

    # homo -> het transition
    for i in range(n_homo):
        c = n_homo  # state
        for h1 in range(n_homo):  # first het
            for h2 in range(h1 + 1, n_homo):  # second het
                if i == h1:  # only transition second haplotype
                    tmat2[i, c] = tmat[h1, h2]
                    tmat2[c, i] = tmat[h2, h1]
                elif i == h2:  # only transition first haplotype
                    tmat2[i, c] = tmat[h2, h1]
                    tmat2[c, i] = tmat[h1, h2]
                else:  # transition both
                    pass
                    # tmat2[i, c] = tmat[i, h1] * tmat[i, h2] + tmat[i, h2] * tmat[j, h1]
                    # tmat2[c, i] = tmat[h1,i] * tmat[h2, i] + tmat[h2, i] * tmat[h1, j]
                c += 1

    # het -> het transition
    c1 = n_homo
    for i in range(n_homo):  # first het from
        for j in range(i + 1, n_homo):  # second het from
            c2 = n_homo
            for h1 in range(n_homo):  # first het of target
                for h2 in range(h1 + 1, n_homo):  # second het of target
                    if i == h1 and j == h2:  # no transition
                        continue
                    elif i == h1:  # transition second haplotype from j to h2
                        tmat2[c1, c2] = tmat[j, h2]
                    elif j == h2:  # transition first haplotype from i to h1
                        tmat2[c1, c2] = tmat[i, h1]
                    else:  # transition both
                        pass
                    c2 += 1
            c1 += 1

    s = np.sum(tmat2, 1)
    for i in range(tmat2.shape[0]):
        tmat2[i, i] -= s[i]

    return tmat2


def init_pars(
    states,
    homo_ids=None,
    het_ids=None,
    sex=None,
    F0=0.001,
    tau0=1,
    e0=1e-2,
    c0=1e-2,
    init_guess=None,
    transition_matrix=None,
    bin_size=1.0,
    **kwargs,
):
    """initialize parameters

    returns a pars object
    """

    n_states, n_hap = states.n_states, states.n_hap

    alpha0 = np.array([1 / n_states] * n_states)
    alpha0_hap = np.array([1 / n_hap] * n_hap)

    if transition_matrix is None:
        trans_mat = np.zeros((n_states, n_states)) + 2e-2
        trans_mat_hap = np.zeros((n_hap, n_hap)) + 2e-2

        np.fill_diagonal(trans_mat, 1 - (n_states - 1) * 2e-2)
        np.fill_diagonal(trans_mat_hap, 1 - (n_hap - 1) * 2e-2)

        if init_guess is not None:
            guess = [i for i, n in enumerate(states.state_names) if n in init_guess]
            logging.info("starting with guess %s " % guess)
            trans_mat[:, guess] = trans_mat[:, guess] + 1
            trans_mat /= np.sum(trans_mat, 1)[:, np.newaxis]
    else:
        trans_mat_hap = pd.read_csv(transition_matrix, header=None).to_numpy()
        trans_mat = trans_mat_hap_to_dip(trans_mat_hap)
        trans_mat_hap = expm(trans_mat_hap * bin_size)
        trans_mat = expm(trans_mat * bin_size)

    cont, error = init_ce(c0, e0)
    F, tau = init_ftau(states.n_homo, F0, tau0)

    return Pars(
        alpha0,
        alpha0_hap,
        trans_mat,
        trans_mat_hap,
        cont,
        error,
        F,
        tau,
        sex=sex,
    )


def posterior_table_slug(pg, data, gtll=None):
    mu = np.sum(pg * np.arange(3) / 2.0, 1)
    random = np.random.binomial(1, np.clip(mu, 0, 1))
    log_g = np.log10(pg + 1e-40)
    log_g = np.minimum(0.0, log_g)
    df = np.hstack((log_g, mu[:, np.newaxis], random[:, np.newaxis]))
    df = pd.DataFrame(df, columns=["G0", "G1", "G2", "p", "random_read"])
    df.random_read = df.random_read.astype(np.uint8)
    if gtll is not None:
        log_ll = np.log10(gtll + 1e-40)
        df_ll = pd.DataFrame(log_ll, columns=["L0", "L1", "L2"])
        df = pd.concat((df, df_ll), axis=1)
    return df


def guess_sex(ref, data, sex_ratio_threshold=0.75):
    """
    guessing the sex of individuals by comparing heterogametic chromosomes.
    By convention, all chromosomes are assumed to be diploid unless they start
    with an `X` or `Z` or `h`
    """
    ref["heterogametic"] = [
        v[0] in "XZxzh" for v in ref.index.get_level_values("chrom")
    ]
    data["heterogametic"] = [
        v[0] in "XZxzh" for v in data.index.get_level_values("chrom")
    ]

    n_sites = ref.heterogametic.value_counts()
    n_reads = data.groupby(data.heterogametic)[["heterogametic", "tref", "talt"]].apply(
        lambda df: np.sum(df.tref + df.talt)
    )
    cov = n_reads / n_sites
    del data["heterogametic"]
    del ref["heterogametic"]

    # no heteogametic data
    if True not in cov:
        return "f"

    if cov[True] / cov[False] < sex_ratio_threshold:
        sex = "m"
        logging.info("guessing sex is male, X/A = %.4f/%.4f" % (cov[True], cov[False]))
    else:
        sex = "f"
        logging.info(
            "guessing sex is female, X/A = %.4f/%.4f" % (cov[True], cov[False])
        )
    return sex


def parse_chroms(arg):
    if arg is None:
        return None
    chroms = []
    for s in arg.split(","):
        if "-" in s:
            a, b = s.split("-")
            chroms.extend([str(s) for s in range(int(a), int(b) + 1)])
        else:
            chroms.append(s)
    return chroms
