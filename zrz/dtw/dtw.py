import argparse
import csv
import json
from array import array
from bisect import bisect_left, bisect_right
from dataclasses import asdict, dataclass
from math import inf, sqrt
from pathlib import Path

# parameters
MODES = {
    'standard': {},
    'time': {'time_penalty': 2},
    'time-dynamics': {'extend_penalty': 0.05, 'time_penalty': 0.5, 'dynamics': 0.3},
    'diagonal-time-dynamics': {'window': 0.25, 'extend_penalty': 0.05, 'time_penalty': 0.5, 'dynamics': 0.3},
    'diagonal-time': {'window': 0.1667, 'extend_penalty': 0.05, 'time_penalty': 2},
}


@dataclass
class Config:
    window: float = 0.0         
    step_penalty: float = 0.0    
    extend_penalty: float = 0.0  
    time_penalty: float = 0.0
    power: float = 2.0
    dynamics: float = 0.0       


def get_speed(values, times):
    n = len(values)
    speed = []
    for i in range(n):
        left = max(0, i - 1)
        right = min(n - 1, i + 1)
        slope = (values[right] - values[left]) / (times[right] - times[left])
        speed.append(abs(slope))
    return speed


def dtw_distance(x, y, times, config):
    n = len(times)

    # rescale time to [0, 1]
    span = times[-1] - times[0]
    t = []
    for time in times:
        t.append((time - times[0]) / span)

    speed_x = get_speed(x, t)
    speed_y = get_speed(y, t)
    open_cost = config.step_penalty + config.open_penalty
    extend_cost = config.step_penalty + config.extend_penalty

    previous_row = []
    for j in range(n + 1):
        previous_row.append([inf, inf, inf])
    previous_row[0][0] = 0.0

    for i in range(1, n + 1):
        current_row = []
        for j in range(n + 1):
            current_row.append([inf, inf, inf])

        # columns inside the window 
        if config.window > 0:
            start = bisect_left(t, t[i - 1] - config.window - 1e-9)
            end = bisect_right(t, t[i - 1] + config.window + 1e-9)
        else:
            start = 0
            end = n

        for j in range(start + 1, end + 1):
            time_gap = abs(t[i - 1] - t[j - 1])
            if config.window > 0 and time_gap > config.window + 1e-12:
                continue

            expression_cost = (x[i - 1] - y[j - 1]) ** 2
            speed_factor = 1 + config.dynamics * (speed_x[i - 1] + speed_y[j - 1]) / 2
            time_cost = config.time_penalty * time_gap ** config.power * speed_factor
            cost = expression_cost + time_cost

            diagonal = previous_row[j - 1]
            left = current_row[j - 1]
            up = previous_row[j]
            current_row[j][0] = cost + min(diagonal)
        
            current_row[j][1] = cost + min(left[0] + open_cost, left[1] + extend_cost, left[2] + open_cost)
            current_row[j][2] = cost + min(up[0] + open_cost, up[1] + open_cost, up[2] + extend_cost)

        previous_row = current_row

    return sqrt(min(previous_row[n]))


def read_matrix(path):
    gene_ids = []
    values = []
    with open(path, encoding='utf-8-sig', newline='') as file:
        reader = csv.reader(file)
        header = next(reader)
        times = []
        for name in header[1:]:
            times.append(float(name[1:]))   
        for row in reader:
            gene_ids.append(row[0])
            values.append([float(value) for value in row[1:]])
    return gene_ids, times, values


def get_args(with_mode=False):
    parser = argparse.ArgumentParser()
    if with_mode:
        parser.add_argument('--mode', choices=list(MODES), required=True)
    parser.add_argument('--input', type=Path, required=True, help='CSV: gene_id,t0,t10,...')
    parser.add_argument('--output', type=Path, required=True, help='output folder')
    return parser.parse_args()


def run(method, distance, input_path, output_path, info=None):
    gene_ids, times, values = read_matrix(input_path)
    n = len(gene_ids)

  
    distances = array('d', [0.0]) * (n * n)
    for i in range(n):
        for j in range(i + 1, n):
            value = distance(values[i], values[j], times)
            distances[i * n + j] = value
            distances[j * n + i] = value    
        if i % 500 == 0:
            print(f'row {i + 1}/{n}')

    output_path.mkdir(parents=True, exist_ok=True)
    with open(output_path / 'distance_matrix.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['gene_id'] + gene_ids)
        for i in range(n):
            row = [gene_ids[i]]
            for j in range(n):
                row.append(distances[i * n + j])
            writer.writerow(row)

    
    record = {'method': method}
    if info is not None:
        record.update(info)
    record['input'] = str(input_path)
    record['genes'] = n
    record['timepoints'] = times
    with open(output_path / 'run.json', 'w', encoding='utf-8') as file:
        json.dump(record, file, indent=2, ensure_ascii=False)
    print('done:', output_path)


def main():
    args = get_args(with_mode=True)
    config = Config(**MODES[args.mode])

    def distance(x, y, times):
        return dtw_distance(x, y, times, config)

    run('dtw', distance, args.input, args.output, {'mode': args.mode, 'params': asdict(config)})


if __name__ == '__main__':
    main()
