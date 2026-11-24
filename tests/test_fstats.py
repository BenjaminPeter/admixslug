from itertools import permutations

import numpy as np
import pandas as pd

from admixslug.slug.fstats import calc_fstats, summarize_f4


def test_f4_counts_each_deleted_block_once():
    # One SFS category per replicate makes each signed f4 a single frequency product.
    frequencies = pd.DataFrame(
        {
            "P": [0.0, 0.5, 1.0],
            "Q": [1.0, 0.0, 0.5],
            "R": [0.5, 1.0, 0.0],
            "S": [1.0, 1.0, 0.5],
            "target": [0.25, 0.75, 0.5],
        }
    )
    sfs = pd.DataFrame(
        {"rep": frequencies.index, "sex_chrom": "autosome", "tau": frequencies["target"]}
    )
    for pop in ["P", "Q", "R", "S"]:
        sfs[f"{pop}_der"] = 2 * frequencies[pop]
        sfs[f"{pop}_anc"] = 2 - sfs[f"{pop}_der"]

    _, _, f4, _ = calc_fstats(sfs, name="target")
    keys = ["A", "B", "C", "D"]
    assert not f4.duplicated(["sex_chrom", *keys, "rep"]).any()
    replicates = f4.set_index(keys + ["rep"]).f4.unstack("rep")
    summary = summarize_f4(f4).set_index(keys)
    quartets = list(permutations(frequencies.columns, 4))
    assert set(replicates.index) == set(quartets)
    for a, b, c, d in quartets:
        expected = (frequencies[a] - frequencies[b]) * (frequencies[c] - frequencies[d])
        # Three equal blocks give jackknife variance (2/3) * sum of squared deviations.
        se = np.sqrt(2 / 3 * ((expected - expected.mean()) ** 2).sum())
        quartet = (a, b, c, d)
        np.testing.assert_allclose(replicates.loc[quartet].astype(float), expected)
        np.testing.assert_allclose(
            summary.loc[quartet, ["f4", "sterr"]].astype(float), [expected.mean(), se]
        )
