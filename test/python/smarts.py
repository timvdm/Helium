import helium


smarts = helium.Smarts()
smarts.init('C')
if smarts.error():
    print 'Could not init smarts:', smarts.error()

mol = helium.Molecule()

#smarts.search(mol)

