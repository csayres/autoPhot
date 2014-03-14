#!/usr/bin/env python
import glob
import os

import matplotlib.pyplot as plt

import autoPhot

def simpleExample(imageDir="/Users/csayres/arcsatData/arcsat/GJ1243/gFilter/oneNight", imageExtension=".FIT"):
    """A simple example.  Will produce output in the same directory
    @param[in] imageDir: string, directory where images are stored
    @param[in] imageExtension: string, extension of images.  All images of this type will be grabbed
    """
    objlist = glob.glob(os.path.join(imageDir, "*"+imageExtension))
    flareCamConfig = autoPhot.config.Config(camera=autoPhot.config.flareCam, phot=autoPhot.config.PhotConfig(), objFileList=objlist)
    photOut = autoPhot.flow.Driver(flareCamConfig)
    photOut.chooseTarget() # figure 1 will pop up, click a target star and close figure
    plt.show(block=True)
    photOut.chooseComparisons() # figure 2 will pop up, click as many comparison stars as you like then close
    plt.show(block=True)
    photOut.crunch()
    outputter = autoPhot.output.Dump(imageDir, photOut)
    outputter.saveDatFile()
    outputter.saveField()

if __name__ == "__main__":
    simpleExample()