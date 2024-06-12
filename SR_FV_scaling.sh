#!/bin/bash

N=10
system="long_hopping"
state="alt"
alpha=1.5

python3 correlators_direct.py "$N" "$system" "$state" "$alpha"
python3 SR_frm_corr.py "$N" "$system" "$state" "$alpha"
