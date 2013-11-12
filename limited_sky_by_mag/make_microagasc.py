import numpy as np

from astropy.table import Table
import tables

print 'Reading miniagasc'
h5 = tables.openFile('miniagasc.h5')
dat = h5.root.data[:]
h5.close()

dat = Table(dat, copy=False)
dat = np.array(dat['AGASC_ID', 'RA', 'DEC', 'PM_RA', 'PM_DEC', 'MAG_ACA', 'MAG_ACA_ERR',
                   'CLASS', 'ASPQ1', 'COLOR1'])

table_desc, bo = tables.descr_from_dtype(dat.dtype)

micro_h5 = tables.openFile('microagasc.h5', mode='w')
micro_tbl = micro_h5.createTable('/', 'data', table_desc, title='AGASC 1.6')
print 'Appending stars to microagasc.h5 file'
micro_tbl.append(dat)
micro_tbl.flush()

print 'Creating indexes in micro_agasc.h5 file'
micro_tbl.cols.RA.createCSIndex()
micro_tbl.cols.DEC.createCSIndex()
micro_tbl.cols.AGASC_ID.createCSIndex()

print 'Flush and close microagasc.h5 file'
micro_tbl.flush()
micro_h5.flush()
micro_h5.close()
