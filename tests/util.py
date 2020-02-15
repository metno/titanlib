def summer_temperature_example():
    directory = '/'.join(__file__.split('/')[0:-1])
    return parse('%s/data/obs_ta_20190601T12Z.txt' % directory)

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
