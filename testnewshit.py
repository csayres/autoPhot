import glob
import flow
import config
import viz
import matplotlib.pyplot as plt

objlist = glob.glob('/Users/csayres/data/arcsat/GJ2069A/GJ1243/gFilter/oneNight/*.FIT')
flareCamConfig = config.Config(camera=config.flareCam, phot=config.PhotConfig(), objFileList=objlist)
photOut = flow.Driver(flareCamConfig)
photOut.chooseTarget()
photOut.chooseComparisons()