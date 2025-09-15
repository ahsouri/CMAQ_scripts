from pykml import parser
from GOESVisualizer import GSVis
import codecs
import numpy as np
import matplotlib.pyplot as plt

coordinates = np.loadtxt('./geolocation/MethaneAIR21rf03.kml', delimiter=",")
#print(coords)
print(np.shape(coordinates))

for h in range(19,25):
    try:
        GSobj = GSVis('east', 2021, 7, 28, h, -112, -99, 36, 42, coordinates, gamma = 3)
        GSobj.plotGS(True,"RF03_" + str(h) +".png")
    except:
        print('no hour')






