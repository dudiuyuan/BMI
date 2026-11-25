"""JIT the team's original arithmetic; leave source files untouched."""
import importlib.util
import sys
import typing
import numpy as np
from numba import njit, prange

sys.dont_write_bytecode = True


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@njit(cache=False)
def left(a, x):
    lo, hi = 0, len(a)
    while lo < hi:
        mid = (lo + hi) // 2
        if a[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo


@njit(cache=False)
def right(a, x):
    lo, hi = 0, len(a)
    while lo < hi:
        mid = (lo + hi) // 2
        if x < a[mid]:
            hi = mid
        else:
            lo = mid + 1
    return lo


@njit(cache=False)
def array_distance(x,y,t,vx,vy,cfg):
    """Same three-state recurrence as the original; use arrays instead of lists."""
    n=len(t)
    previous=np.full((n+1,3),np.inf)
    current=np.empty((n+1,3))
    previous[0,0]=0.
    a=cfg.step_penalty+cfg.open_penalty
    b=cfg.step_penalty+cfg.extend_penalty
    for i in range(1,n+1):
        current[:,:]=np.inf
        for j in range(1,n+1):
            delta=abs(t[i-1]-t[j-1])
            if cfg.window>0 and delta>cfg.window+1e-12:
                continue
            local=((x[i-1]-y[j-1])**2 + cfg.time_penalty*delta**cfg.power*(1+cfg.dynamics*(vx[i-1]+vy[j-1])/2))
            current[j,0]=local+min(previous[j-1,0],previous[j-1,1],previous[j-1,2])
            current[j,1]=local+min(current[j-1,0]+a,current[j-1,1]+b,current[j-1,2]+a)
            current[j,2]=local+min(previous[j,0]+a,previous[j,1]+a,previous[j,2]+b)
        swap=previous;previous=current;current=swap
    return np.sqrt(min(previous[n,0],previous[n,1],previous[n,2]))


def compile_modules(root):
    original = load('go_original_dtw', root / 'zrz/dtw/dtw.py')
    compiled = load('go_compiled_dtw', root / 'zrz/dtw/dtw.py')
    compiled._speed = njit(cache=False)(compiled._speed)
    compiled.bisect_left, compiled.bisect_right = left, right
    function = njit(cache=False)(compiled.dtw_distance)
    ConfigTuple = typing.NamedTuple('ConfigTuple', [(name, float) for name in original.Config.__dataclass_fields__])

    @njit(parallel=True, cache=False)
    def matrix(X, times, cfg):
        span=(times[-1]-times[0]) or 1.0
        t=(times-times[0])/span
        speeds=np.empty_like(X)
        for i in range(len(X)):
            for j in range(len(t)):
                lo,hi=max(0,j-1),min(len(t)-1,j+1)
                speeds[i,j]=abs((X[i,hi]-X[i,lo])/(t[hi]-t[lo])) if len(t)>1 else 0.
        D = np.zeros((len(X), len(X)), dtype=np.float64)
        for i in prange(len(X)):
            for j in range(i + 1, len(X)):
                D[i, j] = D[j, i] = array_distance(X[i], X[j], t,speeds[i],speeds[j],cfg)
        return D

    raw_cluster = load('go_original_cluster', root / 'xqf/xqf_clustering/k_medoids.py')
    fast_cluster = load('go_fast_cluster', root / 'xqf/xqf_clustering/k_medoids.py')
    # Main algorithm and Python random.Random seeding remain the original code.
    fast_cluster.assign_clusters = njit(cache=False)(fast_cluster.assign_clusters)
    member_cost = njit(cache=False)(fast_cluster.distance_to_members)
    # Convert a group's indices once, rather than reflecting the same Python list
    # again for every candidate. The original summation order stays unchanged.
    last = [None, None]

    def distance_to_members(distances, candidate, members):
        if last[0] is not members:
            last[0], last[1] = members, np.asarray(members,dtype=np.int64)
        return member_cost(distances,candidate,last[1])

    fast_cluster.distance_to_members = distance_to_members
    return original, function, matrix, ConfigTuple, raw_cluster, fast_cluster
