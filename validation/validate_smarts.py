smarts_filename = 'smarts.txt'
smiles_filename = '../data/1K.smi'
helium_filename = '../data/1K.hel'

babel_binary = 'babel'
helium_binary = './../build/tools/helium'

import os
import json

def smarts_iterator(filename):
    f = open(filename)
    lines = f.readlines()
    for line in lines:
        yield line[:-1]

def run_openbabel(smarts):
    os.system(babel_binary + ' -s \'' + smarts + '\' ' + smiles_filename + ' -osmi > openbabel_tmp')
    f = open('openbabel_tmp')
    return len(f.readlines())

def run_helium(smarts):
    os.system(helium_binary + ' smarts \'' + smarts + '\' ' + helium_filename + ' > helium_tmp')
    f = open('helium_tmp')
    data = json.load(f)
    return data['num_hits']



for smarts in smarts_iterator(smarts_filename):
    print 'Testing SMARTS: ' + smarts

    openbabel_hits = run_openbabel(smarts)
    helium_hits = run_helium(smarts)

    if openbabel_hits != helium_hits:
        print 'SMARTS: ' + smarts + '  # openbabel hits: ' + str(openbabel_hits) + ', # helium hits: ' + str(helium_hits) + ' (FAIL)'
    

