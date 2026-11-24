"""Factorial simulated copy-state experiment with disjoint development and test cohorts."""

from itertools import product

import numpy as np
from ranking import save


def cohort(seed):
    rng = np.random.default_rng(seed)
    gc = rng.uniform(0.28, 0.72, 100)
    eff = rng.lognormal(0, 0.35, 100)

    def read(cn, batch):
        coverage = rng.uniform(60, 150)
        bias = (batch * 3 + rng.normal(0, 0.6)) * (gc - 0.5)
        return np.maximum(
            1,
            rng.poisson(
                coverage * eff * cn / 2 * np.exp(bias + rng.normal(0, 0.18, 100))
            ),
        )

    neutral = np.full(100, 2)
    refs = [[read(neutral, b) for _ in range(12)] for b in (-1, 1)]
    for i in range(60):
        cn = neutral.copy()
        w = [2, 4, 7, 10][i // 3 % 4]
        start = rng.integers(5, 85)
        if i % 3:
            cn[start : start + w] = i % 3 * 2 - 1
        batch = i % 2
        yield gc, refs, batch, read(cn, 2 * batch - 1), cn


def experiment(seed):
    counts = {}
    for gc, refs, batch, d, y in cohort(seed):
        for norm, matched, correct in product([False, True], repeat=3):
            ref = np.array(refs[batch] if matched else refs[0] + refs[1], dtype=float)
            x = d.astype(float)
            if norm:
                ref /= np.median(ref, axis=1)[:, None]
                x /= np.median(x)
            z = np.log2(x / np.median(ref, axis=0))
            if correct:
                z -= np.polyval(np.polyfit(gc, z, 1), gc)
            smooth = np.median(
                np.lib.stride_tricks.sliding_window_view(np.pad(z, 2, mode="edge"), 5),
                axis=1,
            )
            for minimum, gain in product([1, 3], [0.35, 0.5, 0.65]):
                calls = np.where(smooth < -0.5, 1, np.where(smooth > gain, 3, 2))
                bounds = np.r_[0, np.flatnonzero(np.diff(calls)) + 1, len(calls)]
                for l, r in zip(bounds, bounds[1:]):
                    if r - l < minimum:
                        calls[l:r] = 2
                k = (norm, matched, correct, minimum, gain)
                c = counts.setdefault(k, np.zeros(3, dtype=int))
                c += np.array(
                    [
                        np.sum((y != 2) & (calls == y)),
                        np.sum((calls != 2) & (calls != y)),
                        np.sum((y != 2) & (calls != y)),
                    ]
                )
    return counts


def metrics(k, c):
    tp, fp, fn = map(int, c)
    return dict(
        normalize=k[0],
        matched_controls=k[1],
        gc_correction=k[2],
        minimum_targets=k[3],
        gain_threshold=k[4],
        tp=tp,
        fp=fp,
        fn=fn,
        precision=tp / (tp + fp),
        recall=tp / (tp + fn),
        f1=2 * tp / (2 * tp + fp + fn),
    )


if __name__ == "__main__":
    dev = experiment(504)
    winner = max(dev, key=lambda k: metrics(k, dev[k])["f1"])
    save("sequencing/cnv_development.tsv", [metrics(k, v) for k, v in dev.items()])
    test = experiment(8841)
    save(
        "sequencing/cnv_holdout.tsv",
        [
            dict(model=label, **metrics(k, test[k]))
            for label, k in [
                ("fixed comparator", (True, False, False, 1, 0.35)),
                ("development selection", winner),
            ]
        ],
    )
