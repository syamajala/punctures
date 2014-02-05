import sys
import os
sys.path.append(os.path.expanduser('~/dev/punctures'))

from punctures import *

if __name__ == '__main__':
    setup_view()
    pdisk2 = Pdisk.split_matching(2)
    print pdisk2.alg_basis()

