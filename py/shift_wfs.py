import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import EOBRun_module as EOB
from PyART.models.teob import CreateDict
from PyART.utils.wf_utils import mode_to_k
import seaborn as sns
import argparse
import re

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--par',  type=str,   help='Parameter to shift', required=True)
    parser.add_argument('--min',  type=float, help='Minimum value',      required=True)
    parser.add_argument('--max',  type=float, help='Maximum value',      required=True)
    parser.add_argument('--dp',   type=float, help='Step',               default=None)
    parser.add_argument('--np',   type=int,   help='Number of values',   default=None)
    parser.add_argument('--q',    type=float, help='Mass ratio',         default=1.)
    parser.add_argument('--chi',  type=float, help='Spin',               default=0.)
    parser.add_argument('--chi1',             help='Spin 1',             default=None)
    parser.add_argument('--chi2',             help='Spin 2',             default=None)
    parser.add_argument('--f0',   type=float, help='Initial frequency',  default=0.004)
    parser.add_argument('--ecc',  type=float, help='Eccentricity',       default=0.)

    args = parser.parse_args()

    parlabel = {'delta_a6c'       : r'$\delta a_6^c$',
                'delta_cN3LO'     : r'$\delta c_{\rm{N}^3 \rm{LO}}$',
                'delta_abhf'      : r'$\delta a_{\rm BH}^{\rm f}$',
                'delta_Mbhf'      : r'$\delta M_{\rm BH}^{\rm f}$',
                'delta_Alm_mrg'   : r'$\delta A_{\ell m}^{\rm mrg}$',
                'delta_Omglm_mrg' : r'$\delta \omega_{\ell m}^{\rm mrg}$',
                'delta_Alm_nqc'   : r'$\delta A_{\ell m}^{\rm NQC}$',
                'delta_dAlm_nqc'  : r'$\delta \dot{A}_{\ell m}^{\rm NQC}$',
                'delta_Omglm_nqc' : r'$\delta \omega_{\ell m}^{\rm NQC}$',
                'delta_dOmglm_nqc': r'$\delta \dot{\omega}_{\ell m}^{\rm NQC}$',
                'delta_alphalm0'  : r'$\delta \alpha_{\ell m 0}$',
                'delta_omglm0'    : r'$\delta \omega_{\ell m 0}$'}

    if 'mrg' in args.par or 'nqc' in args.par or 'alpha' in args.par or 'omg' in args.par:
        mode    = re.findall(r'\d+', args.par)[0]
        if len(mode) == 3:
            mode = mode[:2]
        klm     = mode_to_k(int(mode[0]), int(mode[1]))
        parname = args.par.replace(mode, 'lm')
    else:
        parname = args.par
        klm     = None

    if args.chi1 is None and args.chi2 is None:
        chi1 = chi2 = [0., 0., args.chi]
    else:
        if isinstance(eval(args.chi1), float):
            chi1 = [0., 0., eval(args.chi1)]
        elif isinstance(eval(args.chi1), list):
            chi1 = eval(args.chi1)
        
        if isinstance(eval(args.chi2), float):
            chi2 = [0., 0., eval(args.chi2)]
        elif isinstance(eval(args.chi2), list):
            chi2 = eval(args.chi2)
    
    pardic = CreateDict(q=args.q, 
                        chi1x=chi1[0], chi1y=chi1[1], chi1z=chi1[2],
                        chi2x=chi2[0], chi2y=chi2[1], chi2z=chi2[2],
                        f0=args.f0, ecc=args.ecc)
    
    if args.dp is not None and args.np is not None:
        raise ValueError("Specify only one of dp and np!")
    elif args.dp is not None:
        vals = np.arange(args.min, args.max, args.dp, dtype=float)
    elif args.np is not None:
        vals = np.linspace(args.min, args.max, args.np)
    else:
        raise ValueError("Specify one of dp and np!")
    
    
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['text.usetex'] = True
    matplotlib.rcParams['font.size'] = 15

    fig, ax = plt.subplots(2, 2, layout='constrained', sharex='col', figsize=(10, 7), width_ratios=(2, 1))
    cmap = sns.color_palette('coolwarm', as_cmap=True)
    cols = cmap((vals - args.min)/(args.max - args.min))
    
    for jj, parval in enumerate(vals):
        if klm is None:
            pardic[parname] = parval
        else:
            pardic[parname] = {klm: parval}
        
        try:
            t, hp, hc, hlm, dyn = EOB.EOBRunPy(pardic)
            omg22 = np.gradient(hlm['1'][1], t)
    
            ax[0, 0].plot(t, hlm['1'][0], color=cols[jj])
            ax[0, 1].plot(t, hlm['1'][0], color=cols[jj])
            ax[1, 0].plot(t, omg22,       color=cols[jj])
            ax[1, 1].plot(t, omg22,       color=cols[jj])
        except ValueError:
            print("Not this value!")
    
    for any_ax in ax[:, 0]:
        any_ax.set_xlim([1.05*t[0], -75])
    for any_ax in ax[:, 1]:
        any_ax.set_xlim([-75, 75])
    for any_ax in ax[1, :]:
        any_ax.set_xlabel(r'$t/M$')
    ax[0, 0].set_ylabel(r'$|h_{22}|$')
    ax[1, 0].set_ylabel(r'$\omega_{22}$')
    
    norm_cbar = matplotlib.colors.Normalize(vmin=args.min, vmax=args.max)
    cbar      = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm_cbar, cmap=cmap), ax=ax[:, 1])
    cbar.ax.set_ylabel(parlabel[parname])

    plt.show()