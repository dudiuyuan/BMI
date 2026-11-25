"""Silhouette scores and plots for K-medoids."""

from k_medoids import k_medoids


def silhouette_score(distances, labels):
    """Return the mean silhouette score for one clustering."""
    n = len(distances)
    for row in distances:
        if len(row) != n:
            raise ValueError("The distance matrix must be square.")
    if len(labels) != n:
        raise ValueError("There must be one cluster label per gene.")

    groups = {}
    for gene, label in enumerate(labels):
        if label not in groups:
            groups[label] = []
        groups[label].append(gene)
    if not 2 <= len(groups) < n:
        raise ValueError("Silhouette requires between 2 and n - 1 clusters.")

    score_sum = 0.0
    for gene, label in enumerate(labels):
        own_group = groups[label]
        if len(own_group) == 1:
            continue  # Silhouette is zero for a single gene in a cluster.

        # Average distance to the other genes in this cluster.
        total = 0
        for other in own_group:
            if other != gene:
                total += distances[gene][other]
        a = total / (len(own_group) - 1)

        # Find the closest other cluster by average distance.
        b = float("inf")
        for other_label, members in groups.items():
            if other_label == label:
                continue
            total = 0
            for other in members:
                total += distances[gene][other]
            average = total / len(members)
            if average < b:
                b = average

        scale = max(a, b)
        if scale > 0:
            score_sum += (b - a) / scale

    return score_sum / n


def compare_k(distances, k_values, seed=42, n_init=5, max_iter=100):
    """Return the silhouette score for each k."""
    scores = {}
    for k in k_values:
        if not 2 <= k < len(distances):
            raise ValueError("Candidate k values must be between 2 and n - 1.")
        labels, _ = k_medoids(
            distances, k=k, seed=seed, n_init=n_init, max_iter=max_iter
        )
        scores[k] = silhouette_score(distances, labels)
    return scores


def plot_medoid_profiles(expression, time_points, medoids, gene_ids=None):
    """Plot the medoid rows of expression, ordered like the DTW matrix."""
    import matplotlib.pyplot as plt

    if len(medoids) == 0 or len(time_points) == 0:
        raise ValueError("Medoids and sampling times are required.")
    if gene_ids is not None and len(gene_ids) != len(expression):
        raise ValueError("There must be one gene ID per expression row.")
    for gene in medoids:
        if not 0 <= gene < len(expression):
            raise ValueError("Each medoid must index a row of the expression matrix.")
        if len(expression[gene]) != len(time_points):
            raise ValueError("Each profile needs one expression value per time point.")

    fig, ax = plt.subplots(figsize=(8, 4), layout="constrained")
    for cluster, gene in enumerate(medoids):
        label = f"Cluster {cluster}"
        if gene_ids is not None:
            label += f" ({gene_ids[gene]})"
        ax.plot(time_points, expression[gene], marker="o", label=label)
    ax.set_xlabel("Sampling time")
    ax.set_ylabel("Expression")
    ax.set_title("Medoid profiles")
    ax.legend()
    return fig


def plot_k_sensitivity(scores, fixed_k=5):
    """Plot the silhouette scores against k."""
    import matplotlib.pyplot as plt

    if not scores:
        raise ValueError("At least one silhouette score is required.")
    k_values = sorted(scores)
    fig, ax = plt.subplots(figsize=(6, 4), layout="constrained")
    ax.plot(k_values, [scores[k] for k in k_values], marker="o")
    if fixed_k in scores:
        ax.axvline(fixed_k, linestyle="--", color="grey", label=f"Reference k={fixed_k}")
        ax.legend()
    ax.set_xticks(k_values)
    ax.set_xlabel("Number of clusters (k)")
    ax.set_ylabel("Mean silhouette score")
    ax.set_title("Sensitivity to the number of clusters")
    return fig
