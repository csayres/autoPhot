import glob
import flow

objlist = glob.glob('/Users/csayres/data/arcsat/dcSBS1310_561/gFilter/*.FIT')
x = flow.Driver(objlist)
x.chooseTarget()
x.chooseComparisons()