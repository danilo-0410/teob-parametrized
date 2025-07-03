#!/bin/bash

# Generate figures for paper

python shift_wfs.py --par delta_a6c        --min -8   --max 8    --q 1 --chi 0.6 --np 50 --real --align start --tlim 200. --save --hide
python shift_wfs.py --par delta_a6c        --min -100 --max 30   --q 1 --chi 0.6 --np 50 --real --align start --tlim 225. --save --nqc --hide --dyn
python shift_wfs.py --par delta_cN3LO      --min -10  --max 10   --q 1 --chi 0.6 --np 50 --real --align start --tlim 200. --save --hide
python shift_wfs.py --par delta_alpha220   --min -0.8 --max 1    --q 1 --chi 0.6 --np 50 --save --hide --tlim 175.
python shift_wfs.py --par delta_omg220     --min -0.8 --max 1    --q 1 --chi 0.6 --np 50 --save --hide --tlim 175.
python shift_wfs.py --par delta_Mbhf       --min -0.5 --max 1    --q 1 --chi 0.6 --np 50 --save --hide --tlim 175.
python shift_wfs.py --par delta_abhf       --min -0.9 --max 0.16 --q 1 --chi 0.6 --np 50 --save --hide --tlim 175.
python shift_wfs.py --par delta_A22_mrg    --min -0.8 --max 1    --q 1 --chi 0.6 --np 50 --save --hide --tlim 175.
python shift_wfs.py --par delta_Omg22_mrg  --min -0.8 --max 1    --q 1 --chi 0.6 --np 50 --save --hide --tlim 175.
python shift_wfs.py --par delta_A21_mrg    --min -0.8 --max 1    --q 2 --chi 0.6 --np 50 --save --hide --tlim 175.
python shift_wfs.py --par delta_A22_nqc    --min -0.8 --max 1    --q 1 --chi 0.6 --np 50 --save --hide --tlim 175.