import sys
sys.path.append('/home/seshu/dev/punctures')

from punctures import *

if __name__ == '__main__':
#    setup_view()
    pdisk2 = Pdisk.split_matching(2)
    bc = pdisk2.zero_type_D()
    uslide = PUnderslide(pdisk2, 6, 7).dd_mod()
    pdisk2p = uslide.pmc_2.opposite()
    oslide = POverslide(pdisk2p, 1, 2).dd_mod()
    pdisk2pp = oslide.pmc_2.opposite()
    oslide2 = POverslide(pdisk2pp, 4, 5).dd_mod()
    g1 = uslide.mor_to_d(bc)
    g2 = oslide.mor_to_d(g1)
    g3 = oslide2.mor_to_d(g2)

    print bc
    print len(g1.basis)
    print len(g2.basis)
    print len(g3.basis)
