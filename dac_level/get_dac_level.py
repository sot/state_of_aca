
from mica.archive import aca_hdr3
for year in range(2002, 2014):
    print year
    dac = aca_hdr3.MSID('dac', '{}:001'.format(year), '{}:001'.format(year+1))
    np.save('{}_vals'.format(year), dac.vals.data)
    np.save('{}_mask'.format(year), dac.vals.mask)
    np.save('{}_times'.format(year), dac.times)

