import time
import numpy as np
import titanlib

def main():

    N = 30000
    np.random.seed(1000)
    lats = np.random.randn(N) * 1;
    lons = np.random.randn(N) * 1;
    elevs = np.random.rand(N) * 100;
    points = titanlib.Points(lats, lons, elevs)
    values = np.random.randn(N) * 6;
    pos = np.ones(N) * 4;
    neg = np.ones(N) * 4;
    eps2 = np.ones(N) * 0.5;
    iterations = 1
    s_time = time.time()
    flags, sct, rep = titanlib.sct(points, values, 5, 100, 4000, 10000, iterations, 30, 100, 1000, 100, pos, neg, eps2)
    print("SCT time: %.2f s" % (time.time() - s_time))


if __name__ == "__main__":
    main()
