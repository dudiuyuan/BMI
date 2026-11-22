"""Small, deterministic clustering/evaluation helpers. Standard library only.

Clustering accepts distances only. Labels enter evaluation after clustering.
PAM k-medoids and average linkage both operate on saved distances only.
"""
from collections import Counter
from math import isfinite


def _validate_matrix(matrix, k):
    n = len(matrix)
    if not 1 <= k <= n:
        raise ValueError('k must be between 1 and number of genes')
    if any(len(row) != n for row in matrix):
        raise ValueError('Distance matrix must be square')
    if any(not isfinite(v) or v < 0 for row in matrix for v in row):
        raise ValueError('Distances must be finite and nonnegative')
    if any(abs(matrix[i][j]-matrix[j][i]) > 1e-10 for i in range(n) for j in range(i)):
        raise ValueError('Distance matrix must be symmetric')
    if any(abs(matrix[i][i]) > 1e-10 for i in range(n)):
        raise ValueError('Distance diagonal must be zero')


def average_linkage(matrix, k):
    _validate_matrix(matrix, k)
    n = len(matrix)
    groups = {i: [i] for i in range(n)}
    distances = {(i,j): matrix[i][j] for i in range(n) for j in range(i+1,n)}
    next_id = n
    while len(groups) > k:
        a,b = min(distances, key=lambda pair: (distances[pair],pair))
        na,nb = len(groups[a]),len(groups[b])
        updates = {}
        for c in groups:
            if c not in (a,b):
                da=distances[tuple(sorted((a,c)))]; db=distances[tuple(sorted((b,c)))]
                updates[c,next_id]=(na*da+nb*db)/(na+nb)
        distances = {p:v for p,v in distances.items() if a not in p and b not in p}
        distances.update(updates)
        groups[next_id]=groups.pop(a)+groups.pop(b)
        next_id+=1
    clusters = sorted((sorted(g) for g in groups.values()),key=lambda g:g[0])
    assignment=[0]*n
    for c,group in enumerate(clusters):
        for i in group: assignment[i]=c
    return assignment,clusters


def k_medoids(matrix, k):
    """Deterministic PAM: greedy BUILD followed by best improving SWAPs.

    Minimize sum of distances to nearest medoids, without labels or a size
    balancing constraint. Each medoid owns itself on zero-distance ties.
    SWAP terminates at a one-swap local optimum, not a guaranteed global one.
    """
    _validate_matrix(matrix, k)
    n = len(matrix)
    first = min(range(n), key=lambda j: (sum(matrix[i][j] for i in range(n)), j))
    medoids = [first]
    nearest = [matrix[i][first] for i in range(n)]
    while len(medoids) < k:
        candidate = min((j for j in range(n) if j not in medoids),
                        key=lambda j: (sum(min(nearest[i], matrix[i][j]) for i in range(n)), j))
        medoids.append(candidate)
        nearest = [min(nearest[i], matrix[i][candidate]) for i in range(n)]
    medoids.sort()
    while True:
        owner, closest, second = [], [], []
        for i in range(n):
            ranked = sorted(range(k), key=lambda pos: (matrix[i][medoids[pos]],
                                                       i != medoids[pos], medoids[pos]))
            owner.append(ranked[0])
            closest.append(matrix[i][medoids[ranked[0]]])
            second.append(matrix[i][medoids[ranked[1]]] if k > 1 else float('inf'))
        cost = sum(closest)
        best_cost, best_swap = cost, None
        tolerance = 1e-12 * max(1.0, cost)
        for candidate in range(n):
            if candidate in medoids:
                continue
            common = 0.0
            adjustments = [0.0] * k
            for i in range(n):
                distance = matrix[i][candidate]
                base = min(closest[i], distance)
                common += base
                adjustments[owner[i]] += min(second[i], distance) - base
            for pos in range(k):
                swapped_cost = common + adjustments[pos]
                if swapped_cost < best_cost - tolerance:
                    best_cost, best_swap = swapped_cost, (pos, candidate)
        if best_swap is None:
            break
        pos, candidate = best_swap
        medoids[pos] = candidate
        medoids.sort()
    groups = [[] for _ in medoids]
    for i in range(n):
        pos = min(range(k), key=lambda p: (matrix[i][medoids[p]], i != medoids[p], medoids[p]))
        groups[pos].append(i)
    groups.sort(key=lambda group: group[0])
    assignment = [0] * n
    for cluster, group in enumerate(groups):
        for i in group:
            assignment[i] = cluster
    return assignment, groups


def evaluate(labels, assignment, matrix):
    n=len(labels)
    if n != len(assignment) or n < 2:
        raise ValueError('Evaluation needs at least two labeled genes')
    cell=Counter(zip(labels,assignment)); rows=Counter(labels); cols=Counter(assignment)
    choose2=lambda m:m*(m-1)/2
    observed=sum(map(choose2,cell.values()))
    r=sum(map(choose2,rows.values())); c=sum(map(choose2,cols.values()))
    expected=r*c/choose2(n); maximum=(r+c)/2
    ari=(observed-expected)/(maximum-expected) if maximum != expected else 1.0
    purity=sum(max(cell[t,g] for t in rows) for g in cols)/n
    nearest=[min((j for j in range(n) if j!=i),key=lambda j:(matrix[i][j],j)) for i in range(n)]
    return {'ari':ari,'purity':purity,
            'nearest_same_type':sum(labels[i]==labels[j] for i,j in enumerate(nearest))/n}
