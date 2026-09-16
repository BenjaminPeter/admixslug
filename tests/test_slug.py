from admixslug.slug.classes import *
from admixslug.slug.emissions import *
from admixslug.slug.em import *
import numpy as np
import pytest


def test_slug_p_gt_diploid():
    tau = np.array([0.4, 0.2, 0, 1])
    F = np.array([0, 0.1, 1, 1])
    res = np.empty((4, 3))

    SNP2SFS = np.arange(4, dtype=int)

    res = p_gt_diploid(tau, F, SNP2SFS=SNP2SFS)

    pred0 = F * (1 - tau) + (1 - F) * (1 - tau) ** 2
    pred1 = 2 * (1 - F) * tau * (1 - tau)
    pred2 = F * tau + (1 - F) * tau**2

    pred = np.vstack((pred0, pred1, pred2)).T

    assert np.allclose(res - pred, 0)


def test_slug_p_gt_haploid():
    tau0 = np.arange(10) / 10.0
    res = np.empty((10, 3))

    res = p_gt_haploid(tau=tau0, SNP2SFS=np.arange(10, dtype=int))
    assert np.allclose(res[:, 1], 0)
    assert np.allclose(res[:, 2], tau0)
    assert np.allclose(res[:, 0], 1 - tau0)


def test_slug_fwd_p_x():
    """test example to calculate P(X|.)
    - 3 SNPs,
    - SNP3 is flipped (alt = anc)
    - psi is 100% der for all SNPs
    - 4 Reads (2, 1, 1) with alleles (ref, alt, alt, ref)
    - 2 Read groups (0, 1, 0, 1)
    - contamination probs are (0.3, 0)


    """

    pg = np.array([[0, 0.5, 0.5], [1, 0, 0], [1 / 3, 1 / 3, 1 / 3]])
    pa = np.array([1, 1, 1])
    pc = np.array([0.3, 0])

    READS = np.array((0, 1, 1, 0))
    FLIPPED_READS = np.array((0, 0, 0, 1), "bool")

    n_reads = 4
    n_rgs = 2
    n_snps = 3
    READ2RG = np.array([0, 1, 0, 1])
    READ2SNP = np.array([0, 0, 1, 2])

    bwd_x = bwd_p_o_given_x(READS, FLIPPED_READS, 0, 0)
    # Pr(O|X, e=b=0); 1 if X is the same as READS, except for the flipped SNP
    assert np.allclose(bwd_x[:, 1], (0, 1, 1, 1))

    bwd_g1 = bwd_p_one_o_given_g(bwd_x, pa, pc, READ2SNP, READ2RG, n_reads)
    """Pr(O|G)
        - Read1: we observe ref=anc allele, 30% chance cont; Pr(O|G=0) =  .3 *0 + .7 * 1 = 0.7
        - Read2: we observe alt=der allele, 0% chance cont; Pr(O|G=0) =  .0 *0 + 1 * 0 = 0
        - Read3: we observe alt=der allele, 30% chance cont; Pr(O|G=0) =  .3 *1 + .7 * 0 = 0.3
        - Read4: we observe ref=der allele, 0% chance cont; Pr(O|G=0 (anc)) = 0
    """
    assert np.allclose(bwd_g1[:, 0], (0.7, 0, 0.3, 0))

    bwd_g = bwd_p_all_o_given_g(bwd_g1, READ2SNP, n_snps)
    assert np.allclose(bwd_g[0], bwd_g1[0] * bwd_g1[1])
    assert np.allclose(bwd_g[1:], bwd_g1[2:])

    px_cont = fwd_p_x_cont(pa, READ2SNP)  # all 1 in this case
    assert np.allclose(px_cont[:, 1], 1)
    px_nocont = fwd_p_x_nocont(pg, READ2SNP)
    res = fwd_p_x(px_cont, px_nocont, pc, READ2RG)

    """fwd probabilities:
     - reads1, 2: 
        - case1: both are endogenous (0.7 * 1), 
                P(X_1=0 X_2=0| .) = pg[0,1] * .5 * .5 = 1/8
                P(X_1=0 X_2=1| .) = pg[0,1] * .5 * .5 + pg[1,1] * 0  = 1/8
                P(X_1=1 X_2=0| .) = pg[0,1] * .5 * .5 + pg[1,1] * 0  = 1/8
                P(X_1=1 X_2=1| .) = pg[0,1] * .5 * .5 + pg[1,1] * 1  = 5/8
        - case2: only 2 endogenous (0.3 * 1)
                P(X_1=0 X_2=0| .) = 0
                P(X_1=0 X_2=1| .) = 0
                P(X_1=1 X_2=0| .) = pg[0,1] * 0.5 = 0.25
                P(X_1=1 X_2=1| .) = pg[0,1] * 0.5 + pg[1,1] * 1  = 3/4
        total:
                P(X_1=0 X_2=0| .) = 1/8 * 7/10 = 7/80
                P(X_1=0 X_2=1| .) = 7/80
                P(X_1=1 X_2=0| .) = 13/80
                P(X_1=1 X_2=1| .) = 53/80 
        ergo
                p(X_1 = 0) = 14/80 = 0.175
                p(X_2 = 0) = 20/80 = 0.25

      - read3: endo is always 0, cont always 1
            P(X_3=0) = 1 * .7 = 0.7
      - read4: p(cont) = 0, pg is symmetric, P(X_4 = 0.5)
    """
    assert np.allclose(res[:, 0], [14 / 80, 20 / 80, 0.7, 0.5])


def test_slug_fwd_p_x2():
    fwd_g = np.array([[0.5, 0.5, 0]])
    pa = np.array([0.3, 1])
    pc = np.array([0.0])

    READS = np.array((0, 0, 0, 0, 1))
    FLIPPED_READS = np.array([False] * 5)

    n_reads = len(READS)
    n_rgs = 1
    n_snps = 1
    READ2RG = np.array([0, 0, 0, 0, 0])
    READ2SNP = np.array([0, 0, 0, 0, 0])

    bwd_x = bwd_p_o_given_x(READS, FLIPPED_READS, 0.01, 0.01)
    bwd_g1 = bwd_p_one_o_given_g(bwd_x, pa, pc, READ2SNP, READ2RG, n_reads)
    bwd_g = bwd_p_all_o_given_g(bwd_g1, READ2SNP, n_snps)

    px_cont = fwd_p_x_cont(pa, READ2SNP)
    px_nocont = fwd_p_x_nocont(fwd_g, READ2SNP)
    res = fwd_p_x(px_cont, px_nocont, pc, READ2RG)

    assert np.allclose(res[:, 0], 3 / 4)


def test_error_est():
    """simple test dataset for ensuring algorithm is correct"""
    data = SlugReads(
        READS=[0, 0, 0, 0, 1, 1, 1, 1],
        psi=[1.0],
        READ2RG=[0, 0, 0, 0, 0, 0, 1, 1],
        READ2SNP=np.zeros(8, dtype="i"),
        FLIPPED=np.zeros(1, dtype="bool"),
        SNP2SFS=[0],
    )

    pars = SlugPars(cont=[0, 0], tau=[0], F=[0], e=0.01, b=0.001)
    controller = SlugController(update_eb=True, update_ftau=False, update_cont=False)
    update_pars(pars, data, controller)
    print(f"e : {pars.prev_e} -> {pars.e}")
    print(f"b : {pars.prev_b} -> {pars.b}")
    assert pars.e == 1 / 2.0
    assert pars.b == 1 / 2.0

    update_pars(pars, data, controller)
    assert pars.prev_ll <= pars.ll
    assert pars.e == pars.prev_e
    assert pars.b == pars.prev_b
    print(f"ll : {pars.prev_ll} -> {pars.ll}")


def test_slug_bwd_p_o_given_g_gt_sums_to_one():
    """Pr(R | G) must be a valid distribution over the observed genotype call
    R for every fixed true genotype G, in both flip states."""
    R = np.array([0, 1, 2])
    for e, b in [(0.01, 0.02), (0.1, 0.05), (0.2, 0.3)]:
        for flip in (False, True):
            F = np.full(3, flip)
            res = bwd_p_o_given_g_gt(R, F, e, b)
            assert np.allclose(np.sum(res, axis=0), 1)


def test_slug_bwd_p_o_given_g_gt_flip_symmetry():
    """Flipping ref/alt only changes which allele is ancestral/derived, so it
    is equivalent to reversing the genotype (G <-> 2-G) columns of the
    unflipped table; the observed read-counts R are unaffected."""
    R = np.array([0, 1, 2])
    F_noflip = np.zeros(3, dtype=bool)
    F_flip = np.ones(3, dtype=bool)
    e, b = 0.05, 0.1
    res_noflip = bwd_p_o_given_g_gt(R, F_noflip, e, b)
    res_flip = bwd_p_o_given_g_gt(R, F_flip, e, b)
    assert np.allclose(res_flip, res_noflip[:, ::-1])


def test_slug_bwd_p_o_given_g_gt_rr_sums_to_one():
    """Unlike bwd_p_o_given_g_gt, this variant models a single randomly
    sampled read, so R only takes the two values 0 (ref) and 1 (alt) - any
    nonzero R is treated as "observed alt" by the implementation."""
    R = np.array([0, 1])
    for e, b in [(0.01, 0.02), (0.1, 0.05)]:
        for flip in (False, True):
            F = np.full(2, flip)
            res = bwd_p_o_given_g_gt_rr(R, F, e, b)
            assert np.allclose(np.sum(res, axis=0), 1)


def test_slug_bwd_p_o_given_g_gt_rr_flip_symmetry():
    R = np.array([0, 1, 2])
    F_noflip = np.zeros(3, dtype=bool)
    F_flip = np.ones(3, dtype=bool)
    e, b = 0.02, 0.07
    res_noflip = bwd_p_o_given_g_gt_rr(R, F_noflip, e, b)
    res_flip = bwd_p_o_given_g_gt_rr(R, F_flip, e, b)
    assert np.allclose(res_flip, res_noflip[:, ::-1])


def test_update_eb_gt_counts_heterozygote_evidence():
    """Regression test for update_eb_gt: evidence of errors/bias occurring at
    heterozygous-genotype SNPs (calls R=2 or R=0 at a true het, G=1), and
    evidence of correctly-called heterozygotes (R=1, G=1), must be counted in
    the e/b M-step. An earlier version of this function had the relevant
    lines commented out, so this evidence was silently dropped.
    """
    # 10 correctly-called homozygous-ref SNPs (unambiguous "no error" cells)
    R = [0] * 10
    post_g = [[1.0, 0.0, 0.0]] * 10
    # 5 homozygous-ALT calls at true HET sites: a pure "error" cell that only
    # ever shows up in the G=1 (het) column, which the buggy code never read
    R += [2] * 5
    post_g += [[0.0, 1.0, 0.0]] * 5

    R = np.array(R)
    F = np.zeros(len(R), dtype=bool)
    post_g = np.array(post_g)

    e, b = update_eb_gt(post_g, R, F, two_errors=False)
    assert np.isclose(e, 1 / 6)
    assert np.isclose(b, e)

    # Now add correctly-called heterozygotes. They carry no error/bias signal
    # themselves, but they are real evidence that the data is mostly
    # accurate, so they must lower the estimated error rate through the
    # denominator (which the buggy code also never incremented for them).
    R2 = np.concatenate([R, np.ones(20, dtype=R.dtype)])
    F2 = np.zeros(len(R2), dtype=bool)
    post_g2 = np.concatenate([post_g, np.tile([0.0, 1.0, 0.0], (20, 1))])

    e2, b2 = update_eb_gt(post_g2, R2, F2, two_errors=False)
    assert np.isclose(e2, 1 / 14)
    assert e2 < e


def test_update_eb_gt_joint_optimization_symmetric():
    """With symmetric error and bias evidence, the jointly (numerically)
    optimized e and b should come out equal, and match the closed-form
    result of the equivalent single-parameter (two_errors=False) problem."""
    R = [1] * 5 + [1] * 5 + [0] * 10 + [2] * 10
    post_g = (
        [[1.0, 0.0, 0.0]] * 5  # R=1, G=0: pure error cell
        + [[0.0, 0.0, 1.0]] * 5  # R=1, G=2: pure bias cell
        + [[1.0, 0.0, 0.0]] * 10  # R=0, G=0: pure no-error cell
        + [[0.0, 0.0, 1.0]] * 10  # R=2, G=2: pure no-bias cell
    )
    R = np.array(R)
    F = np.zeros(len(R), dtype=bool)
    post_g = np.array(post_g)

    e, b = update_eb_gt(post_g, R, F, two_errors=True)
    assert np.isclose(e, 1 / 6, atol=1e-4)
    assert np.isclose(b, 1 / 6, atol=1e-4)
    assert np.isclose(e, b, atol=1e-6)


def test_update_tau_only_matches_closed_form_when_F_zero():
    """When F is fixed at zero, the numeric solver in update_tau_only should
    reproduce the (exact, in this case) closed-form allele-frequency update
    used previously."""
    F = np.array([0.0])
    old_tau = np.array([0.2])
    post_g = np.array([[3.0, 4.0, 5.0]])
    data = SlugReads(
        READS=[0], psi=[1.0], READ2RG=[0], READ2SNP=[0], SNP2SFS=[0], FLIPPED=[False]
    )

    tau = update_tau_only(F, old_tau, data, post_g)
    expected = (4 / 2.0 + 5) / (3 + 4 + 5)
    assert np.allclose(tau, expected)


def test_update_tau_only_respects_fixed_nonzero_F():
    """Regression test: when F is fixed at a nonzero value, the naive
    (G1/2+G2)/(G0+G1+G2) formula is *not* the maximizer of the Q-function, so
    update_tau_only must not simply fall back to it. Verify the returned tau
    attains a higher Q-function value than the naive formula, and is close to
    a fine grid-search optimum.
    """
    F = np.array([0.9])
    old_tau = np.array([0.5])
    G0, G1, G2 = 2.0, 0.0, 1.0
    post_g = np.array([[G0, G1, G2]])

    data = SlugReads(
        READS=[0], psi=[1.0], READ2RG=[0], READ2SNP=[0], SNP2SFS=[0], FLIPPED=[False]
    )

    def Q(tau, F_):
        c0 = F_ * (1 - tau) + (1 - F_) * (1 - tau) ** 2
        c1 = (1 - F_) * 2 * tau * (1 - tau)
        c2 = F_ * tau + (1 - F_) * tau * tau
        return G0 * np.log(c0) + G1 * np.log(c1) + G2 * np.log(c2)

    naive_tau = (G1 / 2.0 + G2) / (G0 + G1 + G2)  # 1/3, the old formula
    new_tau = update_tau_only(F, old_tau, data, post_g)[0]

    assert Q(new_tau, F[0]) > Q(naive_tau, F[0])

    grid = np.linspace(1e-6, 1 - 1e-6, 200001)
    best_grid_tau = grid[np.argmax(Q(grid, F[0]))]
    assert abs(new_tau - best_grid_tau) < 1e-3


def test_gt_mode_em_runs_and_improves_ll():
    """Smoke/integration test for the --gtmode EM update (update_pars_gt),
    including a heterozygous genotype call, checking the log-likelihood does
    not decrease across iterations."""
    n = 7
    data = SlugReads(
        READS=[0, 0, 0, 1, 1, 2, 2],
        psi=np.ones(n),
        READ2RG=np.zeros(n, dtype=int),
        READ2SNP=np.arange(n),
        SNP2SFS=np.zeros(n, dtype=int),
        FLIPPED=np.zeros(n, dtype=bool),
    )
    pars = SlugPars(cont=[0], tau=[0.4], F=[0], e=0.05, b=0.05)
    controller = SlugController(
        gt_mode=True,
        update_eb=True,
        update_ftau=True,
        update_F=False,
        update_cont=False,
        do_ll=True,
    )

    pars1 = update_pars_gt(pars, data, controller)
    assert pars1.ll >= pars1.prev_ll

    pars2 = update_pars_gt(pars1, data, controller)
    assert pars2.ll >= pars1.ll - 1e-8


def test_cont_est():
    """simple test dataset for ensuring algorithm is correct

    one RG with only endo (all ref) and one RG with only cont ( all 1)
    """
    reads = SlugReads(
        READS=[0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1],
        psi=[1],
        READ2RG=[0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2],
        READ2SNP=[0] * 12,
        SNP2SFS=[0],
        FLIPPED=[False],
    )

    pars = SlugPars(cont=[0.8, 0.5, 0.5], tau=[0], F=[0], e=0.0, b=0.0)
    controller = SlugController(update_eb=False, update_ftau=False, update_cont=True)
    pars = update_pars(pars, reads, controller)
    assert pars.cont[0] == 0
    assert pars.cont[2] == 1
    assert pars.cont[1] == 0.5
    print(f"C = {pars.cont}")
    ll = calc_full_ll(reads, pars)
    assert pars.ll <= ll
    print(f"ll : {pars.ll} -> {ll}")
