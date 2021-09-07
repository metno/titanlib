import time
import numpy as np
import argparse
import sys
import titanlib

def main():
    parser = argparse.ArgumentParser(description='Runs titanlib benchmarks for processing performance')
    parser.add_argument('-j', type=int, help='Do a scaling test, by running on multiple cores (>= 2)', dest='num_cores')
    parser.add_argument('-s', type=float, default=1, help='Enlarge the inputs by this scaling factor to run a bigger test', dest='scaling')
    parser.add_argument('-n', type=int, default=1, help='Number of iterations to average over', dest='iterations')
    parser.add_argument('-t', help='Run only this function', dest="function")

    args = parser.parse_args()

    if args.num_cores is not None and args.num_cores < 2:
        print("Error: Number of cores must be 2 or more")
        sys.exit(1)

    input = dict()
    points = dict()
    np.random.seed(1000)

    N = 10000
    for i in [1000, N]:
        input[i] = np.random.rand(int(i * args.scaling))*3
    for i in [1000, N]:
        # points[i] = titanlib.Points(np.linspace(0, 1, i * args.scaling), np.linspace(0, 1, i * args.scaling), np.zeros(i * args.scaling))
        points[i] = titanlib.Points(np.random.rand(int(i * args.scaling)) * args.scaling,
                np.random.rand(int(i * args.scaling)) * args.scaling, np.random.rand(int(i * args.scaling)))
    run = dict()
    radius = 10000
    run[(titanlib.Points, "1e6")] = {"expected": 0.82, "args":(np.linspace(0, 1, int(1e6 * args.scaling)), np.linspace(0, 1, int(1e6 * args.scaling)), np.zeros(int(1e6 * args.scaling)))}
    run[(titanlib.buddy_check, "1e4")] = {"expected": 0.71, "args":(points[N], input[N],
        np.full(int(N * args.scaling), radius, float), np.ones(int(N * args.scaling), int) * 10, 0.3, 100.0, 0.0, 1.0, 1)}
    run[(titanlib.isolation_check, "1e4")] = {"expected": 0.71, "args":(points[N], 15, 3000)}
    num_min = 10
    num_max = 50
    inner_radius = 5000
    outer_radius = 50000
    num_iterations = 1
    num_min_prof = 50
    min_elev_diff = 100
    min_horizontal_scale = 10000
    vertical_scale = 200
    Nsct = 1000
    run[(titanlib.sct, "1000 in 100 km2")] = {"expected": 3.37, "args":(points[Nsct], input[Nsct], num_min, num_max,
        inner_radius, outer_radius, num_iterations, num_min_prof, min_elev_diff,
        min_horizontal_scale, vertical_scale, np.full(int(Nsct * args.scaling), 4), np.full(int(Nsct *
            args.scaling), 4), np.full(int(Nsct * args.scaling), 0.5))}

    print("Titanlib benchmark (expected results from Intel i7 3.40 Ghz)")
    print("Titanlib version %s" % titanlib.version())
    if args.num_cores is not None and args.num_cores != 1:
        print("Function                             Expected     Time     Diff    Scaling")
    else:
        print("Function                             Expected     Time     Diff")
    num_cores = [1]
    if args.num_cores is not None:
        num_cores += [args.num_cores]
    for key in run.keys()       :
        timings = dict()
        for num_core in num_cores:
            timings[num_core] = 0

        if isinstance(key, tuple):
            name = key[0].__name__ + " " + str(key[1])
            func = key[0]
        else:
            name = key.__name__
            func = key
        if args.function is not None:
            if func.__name__ != args.function:
                continue
        for num_core in num_cores:
            titanlib.set_omp_threads(num_core)
            for it in range(args.iterations):
                s_time = time.time()
                results = func(*run[key]["args"])
                e_time = time.time()
                # print(np.mean(results))
                curr_time = e_time - s_time
                timings[num_core] += curr_time
                # print("%s() Expected: %.2f s Time: %.2f s" % (func.__name__, run[key]["expected"], e_time - s_time))
        for num_core in num_cores:
            timings[num_core] /= args.iterations

        diff = (timings[1] - run[key]["expected"] * args.scaling) / (run[key]["expected"]  * args.scaling) * 100
        string = "%-36s %8.2f %8.2f %8.2f %%" % (name, run[key]["expected"] * args.scaling, timings[1], diff)
        if args.num_cores is not None:
            scaling = timings[1] / timings[args.num_cores] / args.num_cores
            expected = timings[1] / args.num_cores
            scaling = 1 - (timings[args.num_cores] - expected) / (timings[1] - expected)
            # scaling = (1 - timings[args.num_cores] / timings[1]) * (args.num_cores + 1)

            string += " %8.2f %%" % (scaling * 100)
        print(string)
    # gridpp.neighbourhood(input, radius, gridpp.Mean)


if __name__ == "__main__":
    main()
