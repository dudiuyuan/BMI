"""Exact, endpoint-constrained DTW with composable penalties (Python 3.10+).

No third-party packages. Horizontal/vertical states repeat a matched sample;
they are not sequence-deletion gaps. All penalties are nonnegative.
"""
from dataclasses import dataclass
from math import isfinite, sqrt, inf


@dataclass(frozen=True)
class Config:
    window: float | None = 0.0  # 0/None disables the elapsed-time constraint
    step_penalty: float = 0.0
    time_penalty: float = 0.0
    power: float = 2.0
    dynamics: float = 0.0
    open_penalty: float | None = 0.0
    extend_penalty: float | None = 0.0

    def __post_init__(self):
        # Accept nulls in old experiment records; all switches now use numeric zero.
        for key in ('window', 'open_penalty', 'extend_penalty'):
            if getattr(self, key) is None:
                object.__setattr__(self, key, 0.0)
        for key in ('step_penalty', 'time_penalty', 'dynamics', 'open_penalty', 'extend_penalty'):
            value = getattr(self, key)
            if not isfinite(value) or value < 0:
                raise ValueError(f'{key} must be finite and nonnegative')
        if not isfinite(self.power) or self.power <= 0:
            raise ValueError('power is an exponent and must be finite and positive')
        if not isfinite(self.window) or not 0 <= self.window <= 1:
            raise ValueError('window must be in [0, 1]; zero disables the constraint')



@dataclass
class Result:
    cost: float
    distance: float
    path: list[tuple[int, int]] | None
    diagnostics: dict | None


def _values(values):
    values = [float(v) for v in values]
    if not values or not all(isfinite(v) for v in values):
        raise ValueError('Sequences must be nonempty and finite; impute missing values first')
    return values


def _times(times, n):
    if times is None:
        return [i / max(1, n - 1) for i in range(n)]
    t = _values(times)
    if len(t) != n or any(b <= a for a, b in zip(t, t[1:])):
        raise ValueError('Time coordinates must match sequence length and strictly increase')
    return t


def _speed(x, t):
    if len(x) == 1:
        return [0.0]
    result = []
    for i in range(len(x)):
        lo, hi = max(0, i - 1), min(len(x) - 1, i + 1)
        result.append(abs((x[hi] - x[lo]) / (t[hi] - t[lo])))
    return result


def align(x, y, config=None, *, times_x=None, times_y=None, return_path=True):
    """Minimize sum of local costs plus transition penalties, then take sqrt.

    Time coordinates are normalized using the shared time span. Supply both
    coordinate arrays or neither. No post-hoc path-length normalization.
    Ties prefer diagonal, horizontal, vertical states in that order.
    Cost-only mode uses O(len(y)) working memory; traceback uses O(m*n).
    """
    cfg = config or Config()
    x, y = _values(x), _values(y)
    if (times_x is None) != (times_y is None):
        raise ValueError('Supply both time arrays or neither')
    tx, ty = _times(times_x, len(x)), _times(times_y, len(y))
    start, end = min(tx[0], ty[0]), max(tx[-1], ty[-1])
    span = end - start or 1.0
    tx, ty = [(v-start)/span for v in tx], [(v-start)/span for v in ty]
    vx, vy = _speed(x, tx), _speed(y, ty)
    # Each non-diagonal move pays step_penalty plus its run opening/extension.
    a = cfg.step_penalty + cfg.open_penalty
    b = cfg.step_penalty + cfg.extend_penalty
    n = len(y)
    previous = [[inf] * 3 for _ in range(n + 1)]
    previous[0][0] = 0.0
    pointers = {} if return_path else None
    for i in range(1, len(x) + 1):
        current = [[inf] * 3 for _ in range(n + 1)]
        for j in range(1, n + 1):
            delta = abs(tx[i-1] - ty[j-1])
            if cfg.window > 0 and delta > cfg.window + 1e-12:
                continue
            local = (x[i-1]-y[j-1])**2 + cfg.time_penalty * delta**cfg.power * (1 + cfg.dynamics * (vx[i-1]+vy[j-1])/2)
            candidates = [previous[j-1],
                          [current[j-1][s] + (b if s == 1 else a) for s in range(3)],
                          [previous[j][s] + (b if s == 2 else a) for s in range(3)]]
            for state, options in enumerate(candidates):
                source = min(range(3), key=options.__getitem__)
                current[j][state] = local + options[source]
                if return_path:
                    pointers[i, j, state] = source
        previous = current
    state = min(range(3), key=previous[n].__getitem__)
    cost = previous[n][state]
    if not isfinite(cost):
        raise ValueError('No finite endpoint-to-endpoint alignment; check window or numeric scale')
    if not return_path:
        return Result(cost, sqrt(cost), None, None)
    i, j, path = len(x), n, []
    while i and j:
        path.append((i-1, j-1))
        source = pointers[i, j, state]
        if state == 0:
            i, j = i-1, j-1
        elif state == 1:
            j -= 1
        else:
            i -= 1
        state = source
    path.reverse()
    moves = [(q[0]-p[0], q[1]-p[1]) for p, q in zip(path, path[1:])]
    runs, length, last = [], 0, None
    for move in moves:
        if move != (1, 1):
            if move == last:
                length += 1
            else:
                if length:
                    runs.append(length)
                length = 1
        else:
            if length:
                runs.append(length)
            length = 0
        last = move
    if length:
        runs.append(length)
    diagnostics = {'path_length': len(path), 'non_diagonal_steps': sum(runs),
                   'stretch_runs': len(runs), 'stretch_run_lengths': runs,
                   'max_time_displacement': max(abs(tx[i]-ty[j]) for i, j in path)}
    return Result(cost, sqrt(cost), path, diagnostics)
