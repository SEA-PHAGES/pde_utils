from primer3 import *

def filter_thermostable_oligomers(oligomers, tmMin=52, tmMax=58):
    thermostable_oligomers = []
    for oligomer in oligomers:
        if calcTm(oligomer) < tmMin or calcTm(oligomer) > tmMax:
            continue



        thermostable_oligomers.append(oligomer)

        

    return thermostable_oligomers

