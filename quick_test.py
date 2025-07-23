from bfh_python.dstructure import infTypeD, zeroTypeD
from bfh_python.pmc import splitPMC
from bfh_python.arcslide import Arcslide
from bfh_python.arcslideda import ArcslideDA

Z = splitPMC(1)
solid_torus = infTypeD(1)
tau_m = ArcslideDA(Arcslide(Z,1,2))
tau_l = ArcslideDA(Arcslide(Z,2,1))
slides = [tau_m, tau_m, tau_l, tau_l]
for s in slides:
    solid_torus = s.tensorD(solid_torus)
oth_side = zeroTypeD(1)
cf = oth_side.morToD(solid_torus)
cf.simplify()
print(len(cf))

