"""External clustering scores, independent of the team's distance algorithms.

ARI uses Hubert-Arabie adjustment. NMI uses arithmetic-mean entropy.
"""
from collections import Counter
from math import fsum, log


def clustering_scores(truth, prediction):
    truth, prediction = list(truth), list(prediction)
    if not truth or len(truth) != len(prediction):
        raise ValueError("Truth and prediction must have the same positive length")
    n = len(truth)
    rows, cols = Counter(truth), Counter(prediction)
    cells = Counter(zip(truth, prediction))
    choose2 = lambda size: size * (size - 1) // 2
    if n < 2:
        return {"ari": 1.0, "nmi": 1.0}
    observed = sum(choose2(size) for size in cells.values())
    sum_rows = sum(choose2(size) for size in rows.values())
    sum_cols = sum(choose2(size) for size in cols.values())
    # An integer form avoids cancellation when the expected agreement is large.
    numerator = 2 * (observed * choose2(n) - sum_rows * sum_cols)
    denominator = (sum_rows + sum_cols) * choose2(n) - 2 * sum_rows * sum_cols
    ari = numerator / denominator if denominator else 1.0
    h_true = -fsum(count / n * log(count / n) for count in rows.values())
    h_prediction = -fsum(count / n * log(count / n) for count in cols.values())
    mutual = fsum(count / n * log(count * n / (rows[a] * cols[b])) for (a, b), count in cells.items())
    normalizer = (h_true + h_prediction) / 2
    nmi = mutual / normalizer if normalizer else 1.0
    # Remove harmless negative or above-one roundoff; never hide a structural error.
    if -1e-12 < nmi < 0:
        nmi = 0.0
    if 1 < nmi < 1 + 1e-12:
        nmi = 1.0
    return {"ari": ari, "nmi": nmi}
