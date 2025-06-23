import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import EOBRun_module as EOB
from PyART.models.teob import CreateDict
from PyART.utils.wf_utils import mode_to_k, k_to_ell, k_to_emm
import seaborn as sns
import argparse
import re

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--par',               type=str,   help='Parameter to shift', required=True)
    parser.add_argument('--min',               type=float, help='Minimum value',      required=True)
    parser.add_argument('--max',               type=float, help='Maximum value',      required=True)
    parser.add_argument('--dp',                type=float, help='Step',               default=None)
    parser.add_argument('--np',                type=int,   help='Number of values',   default=None)
    parser.add_argument('--q',                 type=float, help='Mass ratio',         default=1.)
    parser.add_argument('--chi',               type=float, help='Spin',               default=0.)
    parser.add_argument('--chi1',                          help='Spin 1',             default=None)
    parser.add_argument('--chi2',                          help='Spin 2',             default=None)
    parser.add_argument('--f0',                type=float, help='Initial frequency',  default=0.004)
    parser.add_argument('--ecc',               type=float, help='Eccentricity',       default=0.)
    parser.add_argument('-lb2', '--LambdaBl2', type=float, help='LambdaBl2',          default=0.)
    parser.add_argument('--real',                          help='Plot real part',     action='store_true')
    parser.add_argument('--dyn',                           help='Plot dynamics',      action='store_true')
    parser.add_argument('--align',             type=str,   help='Where to align',     default='peak')
    parser.add_argument('--save',                          help='Save figure?',       action='store_true')

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
                'delta_taulm0'    : r'$\delta \tau_{\ell m 0}$',
                'delta_omglm0'    : r'$\delta \omega_{\ell m 0}$',
                'd_delta_t_nqc'   : r'$\delta \Delta t_{\rm NQC}$'}

    if args.par != 'd_delta_t_nqc' and ('mrg' in args.par or 'nqc' in args.par or 'alpha' in args.par or 'tau' in args.par or 'omg' in args.par):
        mode    = re.findall(r'\d+', args.par)[0]
        if len(mode) == 3:
            mode = mode[:2]
        klm     = mode_to_k(int(mode[0]), int(mode[1]))
        parname = args.par.replace(mode, 'lm')
        dictq   = True
    else:
        parname = args.par
        klm     = 1
        dictq   = False

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

    l = k_to_ell(klm)
    m = k_to_emm(klm)
    modesvec = [jj for jj in range(max(klm + 1, 2))]
    
    pardic = CreateDict(q=args.q, 
                        chi1x=chi1[0], chi1y=chi1[1], chi1z=chi1[2],
                        chi2x=chi2[0], chi2y=chi2[1], chi2z=chi2[2],
                        f0=args.f0, ecc=args.ecc, use_mode_lm=modesvec,
                        LambdaBl2=args.LambdaBl2, use_nqc=True)
    
    if args.dp is not None and args.np is not None:
        raise ValueError("Specify only one of dp and np!")
    elif args.dp is not None:
        vals = np.arange(args.min, args.max, args.dp, dtype=float)
    elif args.np is not None:
        vals = np.linspace(args.min, args.max, args.np)
    else:
        raise ValueError("Specify one of dp and np!")
    
    # Make sure vanilla wf also drawn, and last
    if 0. in vals:
        j0 = np.where(vals == 0.)[0][0]
        vals = np.delete(vals, j0)
    vals = np.append(vals, 0.)
    
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['text.usetex'] = True
    matplotlib.rcParams['font.size'] = 15

    fig, ax = plt.subplots(2, 2, layout='constrained', sharex='col', figsize=(10, 7), width_ratios=(1, 1))
    if args.dyn:
        fid, ad = plt.subplots(2, 2, layout='constrained', sharex='col', figsize=(10, 7), width_ratios=(1, 1))
    cmap = sns.color_palette('coolwarm', as_cmap=True)
    cols = cmap((vals - args.min)/(args.max - args.min))
    
    for jj, parval in enumerate(vals):
        if not dictq:
            pardic[parname] = parval
            if parval == 0.:
                col = 'k'
            else:
                col = cols[jj]
        else:
            if parval == 0.:
                col = 'k'
            else:
                col = cols[jj]
            pardic[parname] = {klm: parval}
        
        try:
            t, hp, hc, hlm, dyn = EOB.EOBRunPy(pardic)
            omglm = np.gradient(hlm[f'{klm}'][1], t)
    
            if args.align == 'peak':
                tp = t
            elif args.align == 'start':
                tp = t - t[0]
            dyn['t'] = dyn['t'] - dyn['t'][0] + tp[0]
            if args.real:
                ax[0, 0].plot(tp, hlm['1'][0]*np.cos(hlm['1'][1]), color=col)
                ax[0, 1].plot(tp, hlm['1'][0]*np.cos(hlm['1'][1]), color=col)
            else:
                ax[0, 0].plot(tp, hlm[f'{klm}'][0], color=col)
                ax[0, 1].plot(tp, hlm[f'{klm}'][0], color=col)
            ax[1, 0].plot(tp, omglm, color=col)
            ax[1, 1].plot(tp, omglm, color=col)

            if args.dyn:
                ad[0, 0].plot(dyn['t'], dyn['r'], color=col)
                ad[0, 1].plot(dyn['t'], dyn['r'], color=col)
                ad[1, 0].plot(dyn['t'], dyn['MOmega'], color=col)
                ad[1, 1].plot(dyn['t'], dyn['MOmega'], color=col)
        except ValueError:
            print("Not this value!")
    
    for any_ax in ax[:, 0]:
        any_ax.set_xlim([1.05*t[0] + (tp[0] - t[0]), -125 + (tp[0] - t[0])])
    for any_ax in ax[:, 1]:
        any_ax.set_xlim([-125 + (tp[0] - t[0]), 75 + (tp[0] - t[0])])
    for any_ax in ax[1, :]:
        any_ax.set_xlabel(r'$t/M$')
    if args.real:
        ax[0, 0].set_ylabel(r'$\Re h_{{{}{}}}$'.format(l, m))
    else:
        ax[0, 0].set_ylabel(r'$|h_{{{}{}}}|$'.format(l, m))
    ax[1, 0].set_ylabel(r'$\omega_{{{}{}}}$'.format(l, m))

    if args.dyn:
        for any_ax in ad[:, 0]:
            any_ax.set_xlim([1.05*t[0] + (tp[0] - t[0]), -125 + (tp[0] - t[0])])
        for any_ax in ad[:, 1]:
            any_ax.set_xlim([-125 + (tp[0] - t[0]), 125 + (tp[0] - t[0])])
        for any_ax in ad[1, :]:
            any_ax.set_xlabel(r'$t/M$')
        ad[0, 1].set_ylim([0., 7.])
        ad[0, 0].set_ylabel(r'$r/M$')
        ad[1, 0].set_ylabel(r'$\Omega$')
    
    norm_cbar = matplotlib.colors.Normalize(vmin=args.min, vmax=args.max)
    cbar      = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm_cbar, cmap=cmap), ax=ax[:, 1])
    cbar.ax.set_ylabel(parlabel[parname])
    if args.dyn:
        cbar      = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm_cbar, cmap=cmap), ax=ad)
        cbar.ax.set_ylabel(parlabel[parname])

    if args.save:
        fig.savefig(f'figs/{args.par}_{args.min}_{args.max}.png')
        if args.dyn:
            fid.savefig(f'figs/{args.par}_{args.min}_{args.max}_dyn.png')
    plt.show()