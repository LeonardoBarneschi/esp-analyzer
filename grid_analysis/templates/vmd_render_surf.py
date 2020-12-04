#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
from jinja2 import Environment, FileSystemLoader

def vmdsurfrender(template, mymol, dxpot, expath):
    file_loader = FileSystemLoader(os.path.dirname(template))
    env = Environment(loader=file_loader)
    template = env.get_template(os.path.basename(template))
    vmdtemp = template.render( mymol = mymol,
                               dxpot = dxpot)

    with open(expath, "w") as v:
        v.write(vmdtemp)

    return expath


def vmdchrmrender(template, mymol, expath):
    file_loader = FileSystemLoader(os.path.dirname(template))
    env = Environment(loader=file_loader)
    template = env.get_template(os.path.basename(template))
    vmdtemp = template.render( mymol = mymol )

    with open(expath, "w") as v:
        v.write(vmdtemp)

    return expath

if __name__ == "__main__":

    template = sys.argv[1]
    mymol    = sys.argv[2]
    dxpot    = sys.argv[3]
    vmdsurfrender(mymol, dxpot, template, expath)
