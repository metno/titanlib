import time
import numpy as np
import argparse
import collections
import titanlib

def main():
    parser = argparse.ArgumentParser(description='Runs titanlib benchmarks for processing performance')
    parser.add_argument('-j', type=int, help='Do a scaling test, by running on multiple cores (>= 2)', dest='num_cores')
    parser.add_argument('-s', type=float, default=1, help='Enlarge the inputs by this scaling factor to run a bigger test', dest='scaling')
    parser.add_argument('-n', type=int, default=1, help='Number of iterations to average over', dest='iterations')
    parser.add_argument('-t', help='Run only this function', dest="function")

    args = parser.parse_args()

    if args.num_cores is not None and args.num_cores < 2:
        raise Exception("Error: Number of cores must be 2 or more")

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
    radius = 10000
    run = collections.OrderedDict()
    run[("Points", "1e6")] = {"expected": 0.82, "args":(np.linspace(0, 1, int(1e6 * args.scaling)), np.linspace(0, 1, int(1e6 * args.scaling)), np.zeros(int(1e6 * args.scaling)))}
    run[("buddy_check", "1e4")] = {"expected": 0.71, "args":(points[N], input[N],
        np.full(int(N * args.scaling), radius, float), np.ones(int(N * args.scaling), int) * 10, 0.3, 100.0, 0.0, 1.0, 1)}
    run[("isolation_check", "1e4")] = {"expected": 0.71, "args":(points[N], 15, 3000)}
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
    run[("sct", "1000 in 100 km2")] = {"expected": 3.37, "args":(points[Nsct], input[Nsct], num_min, num_max,
        inner_radius, outer_radius, num_iterations, num_min_prof, min_elev_diff,
        min_horizontal_scale, vertical_scale, np.full(int(Nsct * args.scaling), 4), np.full(int(Nsct *
            args.scaling), 4), np.full(int(Nsct * args.scaling), 0.5))}

    if args.num_cores is not None:
        print("Titanlib parallelization test (titanlib version %s)" % titanlib.version())
    else:
        print("Titanlib benchmark (titanlib version %s)" % titanlib.version())
        print("Expected results from Intel i7 3.40 Ghz")
    print("-----------------------------------------------------------------")
    if args.num_cores is not None:
        print("Function                               1 core %2d cores  Scaling" % args.num_cores)
    else:
        print("Function                             Expected     Time     Diff")
    num_cores = [1]
    if args.num_cores is not None:
        num_cores += [args.num_cores]
    for key in run.keys()       :
        try:
            timings = dict()
            for num_core in num_cores:
                timings[num_core] = 0

            if isinstance(key, tuple):
                name = key[0] + " " + str(key[1])
                func = eval("titanlib." + key[0])
            else:
                name = key
                func = eval("titanlib." + key)
            if args.function is not None:
                if func.__name__ != args.function:
                    continue

            # Allow functions to fail (useful when benchmarking older versions of gridpp
            # where functions may not be defined).
            for num_core in num_cores:
                titanlib.set_omp_threads(num_core)
                for it in range(args.iterations):
                    s_time = time.time()
                    func(*run[key]["args"])
                    e_time = time.time()
                    curr_time = e_time - s_time
                    timings[num_core] += curr_time
        except Exception as e:
            print("Could not run", key, e)
            continue

        for num_core in num_cores:
            timings[num_core] /= args.iterations

        if args.num_cores is None:
            diff = (timings[1] - run[key]["expected"] * args.scaling) / (run[key]["expected"]  * args.scaling) * 100
            string = "%-36s %8.2f %8.2f %8.2f %%" % (name, run[key]["expected"] * args.scaling, timings[1], diff)
        else:
            scaling = timings[1] / timings[args.num_cores] / args.num_cores
            expected = timings[1] / args.num_cores
            scaling = 1 - (timings[args.num_cores] - expected) / (timings[1] - expected)
            # scaling = (1 - timings[args.num_cores] / timings[1]) * (args.num_cores + 1)

            string = "%-36s %8.2f %8.2f %8.2f %%" % (name, timings[1], timings[args.num_cores], scaling * 100)
        print(string)


if __name__ == "__main__":
    main()
