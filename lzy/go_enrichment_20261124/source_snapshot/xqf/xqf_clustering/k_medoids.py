"""K-medoids clustering with a DTW distance matrix."""

import math
import random


def assign_clusters(distances, medoids):
    labels = [-1] * len(distances)

    # Keep each medoid in its own cluster, even if two medoids are identical.
    for cluster, gene in enumerate(medoids):
        labels[gene] = cluster

    for gene, row in enumerate(distances):
        if labels[gene] != -1:
            continue

        closest = 0
        for cluster in range(1, len(medoids)):
            candidate = medoids[cluster]
            current_best = medoids[closest]

            if row[candidate] < row[current_best]:
                closest = cluster

        labels[gene] = closest
    return labels


def distance_to_members(distances, candidate, members):
    total = 0
    for gene in members:
        total += distances[gene][candidate]
    return total


def k_medoids(distances, k, seed=42, n_init=5, max_iter=100):
    """Return the cluster label of each gene and the medoid of each cluster.

    distances[i][j] is the distance between genes i and j.
    Labels and medoid indices start at 0.
    """
    n = len(distances)
    if not 1 <= k <= n:
        raise ValueError("k must be between 1 and the number of genes.")
    if n_init < 1 or max_iter < 1:
        raise ValueError("n_init and max_iter must be positive.")
    for row in distances:
        if len(row) != n:
            raise ValueError("The distance matrix must be square.")
    for i in range(n):
        if distances[i][i] != 0:
            raise ValueError("Diagonal distances must be zero.")
        for j in range(i):
            distance = distances[i][j]
            if not math.isfinite(distance) or distance < 0:
                raise ValueError("Distances must be finite and non-negative.")
            if distance != distances[j][i]:
                raise ValueError("The distance matrix must be symmetric.")

    rng = random.Random(seed)
    best_cost = None

    for _ in range(n_init):
        # 1. Pick k different genes at random.
        medoids = rng.sample(range(n), k)

        for _ in range(max_iter):
            # 2. Assign each gene to the closest medoid.
            labels = assign_clusters(distances, medoids)
            groups = [[] for _ in range(k)]
            for gene, label in enumerate(labels):
                groups[label].append(gene)

            # 3. Pick the gene with the smallest total distance in each cluster.
            updated_medoids = []
            for cluster, members in enumerate(groups):
                current_medoid = medoids[cluster]
                best_gene = current_medoid
                lowest_cost = distance_to_members(distances, current_medoid, members)

                for candidate in members:
                    if candidate == current_medoid:
                        continue
                    cost = distance_to_members(distances, candidate, members)
                    if cost < lowest_cost:
                        best_gene = candidate
                        lowest_cost = cost
                updated_medoids.append(best_gene)

            if updated_medoids == medoids:
                break
            medoids = updated_medoids

        # Assign again in case the last update reached the iteration limit.
        labels = assign_clusters(distances, medoids)
        cost = 0
        for gene, label in enumerate(labels):
            cost += distances[gene][medoids[label]]

        if best_cost is None or cost < best_cost:
            best_cost = cost
            best_labels = labels
            best_medoids = medoids

    return best_labels, best_medoids
