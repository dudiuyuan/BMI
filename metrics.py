from collections import Counter
from math import fsum, log


def choose_two(number):
    return number * (number - 1) // 2


def entropy(counts, total):
    terms = []
    for count in counts:
        probability = count / total
        terms.append(-probability * log(probability))
    return fsum(terms)


def clustering_scores(truth, prediction):
    truth = list(truth)
    prediction = list(prediction)
    if not truth or len(truth) != len(prediction):
        raise ValueError("Truth and prediction must have the same positive length")
    n = len(truth)
    if n == 1:
        return {"ari": 1.0, "nmi": 1.0}

    true_counts = Counter(truth)
    predicted_counts = Counter(prediction)
    joint_counts = Counter(zip(truth, prediction))

    same_in_both = 0
    for count in joint_counts.values():
        same_in_both += choose_two(count)
    same_in_truth = 0
    for count in true_counts.values():
        same_in_truth += choose_two(count)
    same_in_prediction = 0
    for count in predicted_counts.values():
        same_in_prediction += choose_two(count)

    all_pairs = choose_two(n)
    numerator = 2 * (same_in_both * all_pairs - same_in_truth * same_in_prediction)
    denominator = (same_in_truth + same_in_prediction) * all_pairs - 2 * same_in_truth * same_in_prediction
    if denominator == 0:
        ari = 1.0
    else:
        ari = numerator / denominator

    terms = []
    for (true_label, predicted_label), count in joint_counts.items():
        ratio = count * n / (true_counts[true_label] * predicted_counts[predicted_label])
        terms.append(count / n * log(ratio))
    mutual_information = fsum(terms)
    true_entropy = entropy(true_counts.values(), n)
    predicted_entropy = entropy(predicted_counts.values(), n)
    normalizer = (true_entropy + predicted_entropy) / 2
    if normalizer == 0:
        nmi = 1.0
    else:
        nmi = mutual_information / normalizer
    if -1e-12 < nmi < 0:
        nmi = 0.0
    if 1 < nmi < 1 + 1e-12:
        nmi = 1.0
    return {"ari": ari, "nmi": nmi}
