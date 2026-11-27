"""Array implementation of frozen xqf alternating k-medoids; no fastmath."""
import random
import numpy as np
from numba import njit


@njit(cache=False)
def assign(D,medoids):
    labels=np.full(len(D),-1,np.int64)
    for c,g in enumerate(medoids): labels[g]=c
    for g in range(len(D)):
        if labels[g]!=-1: continue
        closest=0
        for c in range(1,len(medoids)):
            if D[g,medoids[c]]<D[g,medoids[closest]]: closest=c
        labels[g]=closest
    return labels


@njit(cache=False)
def cost(D,candidate,members):
    total=0.
    for g in members: total+=D[g,candidate]
    return total


@njit(cache=False)
def core(D,starts,max_iter):
    best_cost=np.inf
    best_labels=np.empty(len(D),np.int64)
    best_medoids=starts[0].copy()
    for initial in starts:
        medoids=initial.copy()
        for iteration in range(max_iter):
            labels=assign(D,medoids)
            updated=medoids.copy()
            for c in range(len(medoids)):
                members=np.where(labels==c)[0]
                best_gene=medoids[c];lowest=cost(D,best_gene,members)
                for candidate in members:
                    if candidate==medoids[c]: continue
                    value=cost(D,candidate,members)
                    if value<lowest:
                        best_gene=candidate;lowest=value
                updated[c]=best_gene
            if np.array_equal(updated,medoids): break
            medoids=updated
        labels=assign(D,medoids)
        objective=0.
        for i,c in enumerate(labels): objective+=D[i,medoids[c]]
        if objective<best_cost:
            best_cost=objective;best_labels=labels.copy();best_medoids=medoids.copy()
    return best_labels,best_medoids,best_cost


def fit(D,k=5,seed=42,n_init=5,max_iter=100):
    rng=random.Random(seed)
    starts=np.array([rng.sample(range(len(D)),k) for _ in range(n_init)],dtype=np.int64)
    return core(D,starts,max_iter)


@njit(cache=False)
def pair_score(labels,edges):
    numerator=0
    for pair in edges:
        if labels[pair[0]]==labels[pair[1]]:numerator+=1
    counts=np.bincount(labels,minlength=5)
    denominator=0
    for n in counts:denominator+=n*(n-1)//2
    return numerator/denominator if denominator else np.nan


@njit(cache=False)
def null_scores(labels,edges,n_permutations,seed):
    np.random.seed(seed)
    scores=np.empty(n_permutations)
    for b in range(n_permutations):
        scores[b]=pair_score(np.random.permutation(labels),edges)
    return scores
