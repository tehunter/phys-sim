#!/usr/bin/env python3

import os
import sys
sys.path.insert(0, os.getcwd())

from pypsim.boid import *

def main():
    vis = BoidVisualizer()
    vis.run()

if __name__ == '__main__':
    main()
