import numpy as np


def summer_temperature_example():
    directory = '/'.join(__file__.split('/')[0:-1])
    np.random.seed(1)
    N = 10000
    lats = np.random.randn(N) + 60
    lons = np.random.randn(N) + 10.7
    elevs = np.random.rand(N) * 100
    values = np.random.rand(N) * 2 + 20

    return lats, lons, elevs, values


def parse(filename):
    with open(filename, 'r') as file:
        header = file.readline().strip().split(';')
        lats = list()
        lons = list()
        elevs = list()
        values = list()
        Ilat = header.index('lat')
        Ilon = header.index('lon')
        Ielev = header.index('elev')
        Ivalue = header.index('value')
        for line in file:
            words = line.strip().split(';')
            lats += [float(words[Ilat])]
            lons += [float(words[Ilon])]
            elevs += [float(words[Ielev])]
            values += [float(words[Ivalue])]
        return lats, lons, elevs, values
