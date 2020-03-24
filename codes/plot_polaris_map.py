import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as clr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import Slider
import astropy.io.fits as fits
from astropy.convolution import convolve, Gaussian2DKernel
import astropy.units as u
from astropy.visualization import quantity_support
quantity_support()
import argparse
import os
import copy

import salamander as sm


#=========================================================================================================================#
#                                                                                                                         #
#=========================================================================================================================#


prim_map_r = cm.get_cmap('inferno_r')
prim_map = cm.get_cmap('inferno')
color_scale_r = np.zeros((255,4))
color_scale = np.zeros((255,4))
for i in range(255):
    color_scale[i] = prim_map(i)
    color_scale_r[i] = prim_map_r(i)
# set max value to white
color_scale[-1] = (1,1,1,1)
# set min value to white
color_scale_r[0] = (1,1,1,1)
my_cmap = clr.ListedColormap(color_scale)
my_cmap_r = clr.ListedColormap(color_scale_r)

# map bad pixels (e.g. 0 in LogNorm) with the lowest value of cmap (default for bad is transparent)
my_cmap.set_bad(my_cmap(0))
my_cmap_r.set_bad(my_cmap_r(0))


aspect_ratio = 1.1
height_inch  = 7

left, width = .06, .85
bottom, height = .06, .9
right = left + width
top = bottom + height


#=========================================================================================================================#
#                                                                                                                         #
#=========================================================================================================================#


HOME = os.path.expanduser("~")

WORK_DIR = HOME + '/polaris/'
POLARIS_DIR = WORK_DIR
POLARIS_BIN_DIR = POLARIS_DIR + 'bin/'


#=========================================================================================================================#
#                                                                                                                         #
#=========================================================================================================================#


parser = argparse.ArgumentParser(
    description=('Plot results of POLARIS dust_mc runs'),
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )


parser.add_argument('--model_name', dest='model_name', type=str, required=True,
#                    choices=['disk','modified_disk','yang_disk','sphere','globule','kata_disk',
 #                            'multi_LB-P_disk'],
                    help='name of the PolarisTools model.')

parser.add_argument('--project_name', dest='pname', type=str, required=True,
                    help='name of the polaris project.')

parser.add_argument('--rad_prof',dest='radial_profile', action='store_true')
parser.set_defaults(radial_profile=False)

parser.add_argument('--simulation',dest='simulation', nargs='*', type=str,
                    help='POLARIS simulation type.', default=['dust'],
                    choices=['dust','dust_mc'])

parser.add_argument('--fwhm', type=float, default=0,
                    help='FWHM of the 2D Gaussian. [px]')

parser.add_argument('--no_pol_vec',dest='pol_vectors', action='store_false')
parser.set_defaults(pol_vectors=True)

parser.add_argument('--bin_size', dest='bin_size', type=int, default=5,
                    help='Bin size of radial bins in pixel for radial profiles.')

parser.add_argument('--vec_stretch', dest='vec_stretch', type=float,
                    default=1, help='Stretch vector arrows by this factor.')

parser.add_argument('--vmin_I', dest='vmin_I', type=float,
                    default=1e-6, help='Minimum relative intensity to compute polarisation degree.')

parser.add_argument('--unit', dest='unit', type=str, default='au',
                    choices=['au','as'], help='name of the polaris project.')

args = parser.parse_args()


#=========================================================================================================================#
#                                                                                                                         #
#=========================================================================================================================#


matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'font.family': 'serif'})
matplotlib.rcParams.update({'axes.labelsize': 'large'})
matplotlib.rcParams.update({'image.interpolation': 'none'})
matplotlib.rcParams.update({'figure.figsize': [aspect_ratio*height_inch,height_inch]})
matplotlib.rcParams.update({'savefig.dpi': 400})
matplotlib.rcParams.update({'savefig.bbox': 'tight'})
matplotlib.rcParams.update({'savefig.format': 'pdf'})

matplotlib.rcParams.update({'text.usetex': True})
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage[T1]{fontenc}'
                                            , r'\usepackage[utf8]{inputenc}'
                                            , r'\usepackage{newtxtext}'
                                            , r'\usepackage{newtxmath}'
                                            , r'\usepackage[tight-spacing=true]{siunitx}'
                                            , r'\DeclareMathAlphabet\mathcal{OMS}{txsy}{m}{n}'
                                            , r'\DeclareSIUnit\jy{Jy}'
                                            , r'\DeclareSIUnit\au{au}'
                                            , r'\DeclareSIUnit\msun{M\ensuremath{_{\odot}}}'
                                            ]


#=========================================================================================================================#
#                                                                                                                         #
#=========================================================================================================================#


polaris_project = sm.polaris.Polaris(args.pname, args.model_name, POLARIS_DIR=POLARIS_DIR)

first_simu = True
stokes = {}
R_map = {}
VOV_map = {}
N_photon_map = {}

for simu_type in args.simulation:
    stokes_tmp, wave, boundary_map, pixel, further_maps, as_per_pix = polaris_project.get_continuum_maps(simu_type=simu_type)
    if first_simu:
        stokes['I'] = {'data': stokes_tmp['I'] / as_per_pix**2}
        stokes['Q'] = {'data': stokes_tmp['Q'] / as_per_pix**2}
        stokes['U'] = {'data': stokes_tmp['U'] / as_per_pix**2}
        stokes['V'] = {'data': stokes_tmp['V'] / as_per_pix**2}
        if simu_type == 'dust_mc':
            stokes['I_scatter'] = {'data': further_maps['I_scatter'] / as_per_pix**2}
        else:
            stokes['I_scatter'] = {'data': np.zeros_like(stokes_tmp['I']) / as_per_pix**2}
        first_simu = False
    else:
        stokes['I']['data'] = stokes['I']['data'] + stokes_tmp['I'] / as_per_pix**2
        stokes['Q']['data'] = stokes['Q']['data'] + stokes_tmp['Q'] / as_per_pix**2
        stokes['U']['data'] = stokes['U']['data'] + stokes_tmp['U'] / as_per_pix**2
        stokes['V']['data'] = stokes['V']['data'] + stokes_tmp['V'] / as_per_pix**2
        if simu_type == 'dust_mc':
            stokes['I_scatter']['data'] = stokes['I_scatter']['data'] + further_maps['I_scatter'] / as_per_pix**2

    if simu_type == 'dust_mc' and len(further_maps) > 2:
        N_photon_map = {'data': further_maps['N_photon']}
        R_map['I'] = {'data': further_maps['R_I']}
        R_map['Q'] = {'data': further_maps['R_Q']}
        R_map['U'] = {'data': further_maps['R_U']}
        R_map['V'] = {'data': further_maps['R_V']}
        VOV_map['I'] = {'data': further_maps['VOV_I']}
        VOV_map['Q'] = {'data': further_maps['VOV_Q']}
        VOV_map['U'] = {'data': further_maps['VOV_U']}
        VOV_map['V'] = {'data': further_maps['VOV_V']}
    else:
        N_photon_map = {'data': np.ones_like(stokes['I']) * u.dimensionless_unscaled}
        R_map['I'] = {'data': np.zeros_like(stokes['I']) * u.dimensionless_unscaled}
        R_map['Q'] = {'data': np.zeros_like(stokes['I']) * u.dimensionless_unscaled}
        R_map['U'] = {'data': np.zeros_like(stokes['I']) * u.dimensionless_unscaled}
        R_map['V'] = {'data': np.zeros_like(stokes['I']) * u.dimensionless_unscaled}
        VOV_map['I'] = {'data': np.zeros_like(stokes['I']) * u.dimensionless_unscaled}
        VOV_map['Q'] = {'data': np.zeros_like(stokes['I']) * u.dimensionless_unscaled}
        VOV_map['U'] = {'data': np.zeros_like(stokes['I']) * u.dimensionless_unscaled}
        VOV_map['V'] = {'data': np.zeros_like(stokes['I']) * u.dimensionless_unscaled}

n_wave = len(wave)

# imshow uses (y,x) indexing, so mrange must also be (y,x)
mrange = np.array([-boundary_map['y'][args.unit].value, boundary_map['y'][args.unit].value,
                   -boundary_map['x'][args.unit].value, boundary_map['x'][args.unit].value]) * boundary_map['y'][args.unit].unit
# if unit is as, then x axis (right ascension) is inverted to match on-sky view
if args.unit == 'as':
    mrange[2:] *= -1

if args.fwhm > 0:
    stdev = args.fwhm / ( 2*np.sqrt(2*np.log(2)) )
    for i_wave in range(n_wave):
        stokes['I']['data'][i_wave] = convolve(stokes['I']['data'][i_wave],Gaussian2DKernel(stdev,stdev))
        stokes['Q']['data'][i_wave] = convolve(stokes['Q']['data'][i_wave],Gaussian2DKernel(stdev,stdev))
        stokes['U']['data'][i_wave] = convolve(stokes['U']['data'][i_wave],Gaussian2DKernel(stdev,stdev))
        stokes['V']['data'][i_wave] = convolve(stokes['V']['data'][i_wave],Gaussian2DKernel(stdev,stdev))
        stokes['I_scatter']['data'][i_wave] = convolve(stokes['I_scatter']['data'][i_wave],Gaussian2DKernel(stdev,stdev))


stokes['I_norm'] = {'data': np.empty((n_wave,pixel[0].value,pixel[1].value))}
for i_wave in range(n_wave):
    stokes['I_norm']['data'][i_wave] = stokes['I']['data'][i_wave]/np.nanmax(stokes['I']['data'][i_wave])

stokes['PI'] = {'data': []}
stokes['P'] = {'data': []}
stokes['P_angle'] = {'data': []}

stokes['PI']['data'], stokes['P']['data'], stokes['P_angle']['data'] = sm.get_polarization(stokes['I']['data'],
                                                                                           stokes['Q']['data'],
                                                                                           stokes['U']['data'],
                                                                                           stokes['V']['data'],
                                                                                           I_norm=stokes['I_norm']['data'],
                                                                                           vmin=args.vmin_I)

pol_vectors = sm.get_pol_vectors(stokes['P']['data'].value,stokes['P_angle']['data'].value, mrange.value)

vec_scale = 0.8 * pixel[0].value/ args.vec_stretch

# if args.radial_profile:
#     radial_profile = {'I': []*n_wave,
#                       'PI': []*n_wave,
#                      }
#     for i_wave in range(n_wave):
#         x_radial_profile, radial_profile['I'][i_wave], std = sm.radial_profile(stokes['I'][i_wave,:,:],
#                                                                                 center=((pixel[0]-1)/2,(pixel[1]-1)/2),
#                                                                                 bin_size=args.bin_size)
#         x_radial_profile, radial_profile['PI'][i_wave], std = sm.radial_profile(stokes['PI'][i_wave,:,:],
#                                                                                 center=((pixel[0]-1)/2,(pixel[1]-1)/2),
#                                                                                 bin_size=args.bin_size)


#=========================================================================================================================#
#                                                                                                                         #
#=========================================================================================================================#


stokes['I']['string'] = 'I'
stokes['Q']['string'] = 'Q'
stokes['U']['string'] = 'U'
stokes['V']['string'] = 'V'
stokes['I_scatter']['string'] = r'I$_{\text{scatter}}$'
stokes['PI']['string'] = r'P$\cdot$I'
stokes['P']['string'] = 'P'

R_map['I']['string'] = r'R$_{\text{I}}$'
R_map['Q']['string'] = r'R$_{\text{Q}}$'
R_map['U']['string'] = r'R$_{\text{U}}$'
R_map['V']['string'] = r'R$_{\text{V}}$'

VOV_map['I']['string'] = r'VOV$_{\text{I}}$'
VOV_map['Q']['string'] = r'VOV$_{\text{Q}}$'
VOV_map['U']['string'] = r'VOV$_{\text{U}}$'
VOV_map['V']['string'] = r'VOV$_{\text{V}}$'

N_photon_map['string'] = r'N$_{\text{photon}}$'



stokes['I']['clabel'] = r'Flux density [{}]'.format(stokes['I']['data'].unit)
stokes['Q']['clabel'] = r'Flux density [{}]'.format(stokes['Q']['data'].unit)
stokes['U']['clabel'] = r'Flux density [{}]'.format(stokes['U']['data'].unit)
stokes['V']['clabel'] = r'Flux density [{}]'.format(stokes['V']['data'].unit)
stokes['I_scatter']['clabel'] = r'Flux density [{}]'.format(stokes['I_scatter']['data'].unit)
stokes['PI']['clabel'] = r'lin. pol. Flux density [{}]'.format(stokes['PI']['data'].unit)
stokes['P']['clabel'] = r'Degree of polarisation [{}]'

R_map['I']['clabel'] = 'R value for I'
R_map['Q']['clabel'] = 'R value for Q'
R_map['U']['clabel'] = 'R value for U'
R_map['V']['clabel'] = 'R value for V'

VOV_map['I']['clabel'] = 'VOV value for I'
VOV_map['Q']['clabel'] = 'VOV value for Q'
VOV_map['U']['clabel'] = 'VOV value for U'
VOV_map['V']['clabel'] = 'VOV value for V'

N_photon_map['clabel'] = 'Number of photons'



stokes['I']['vmin'] = lambda max_val,min_val : max_val*1e-4
stokes['Q']['vmin'] = lambda max_val,min_val : min_val
stokes['U']['vmin'] = lambda max_val,min_val : min_val
stokes['V']['vmin'] = lambda max_val,min_val : min_val
stokes['I_scatter']['vmin'] = lambda max_val,min_val : max_val*1e-4
stokes['PI']['vmin'] = lambda max_val,min_val : max_val*1e-4
stokes['P']['vmin'] = lambda max_val,min_val : 0*u.percent

R_map['I']['vmin'] = lambda max_val,min_val : 0*u.dimensionless_unscaled
R_map['Q']['vmin'] = lambda max_val,min_val : 0*u.dimensionless_unscaled
R_map['U']['vmin'] = lambda max_val,min_val : 0*u.dimensionless_unscaled
R_map['V']['vmin'] = lambda max_val,min_val : 0*u.dimensionless_unscaled

VOV_map['I']['vmin'] = lambda max_val,min_val : 0*u.dimensionless_unscaled
VOV_map['Q']['vmin'] = lambda max_val,min_val : 0*u.dimensionless_unscaled
VOV_map['U']['vmin'] = lambda max_val,min_val : 0*u.dimensionless_unscaled
VOV_map['V']['vmin'] = lambda max_val,min_val : 0*u.dimensionless_unscaled

N_photon_map['vmin'] = lambda max_val,min_val : 1*u.dimensionless_unscaled



stokes['I']['norm'] = lambda vmin,vmax,linthresh : clr.LogNorm(vmin=vmin,
                                                               vmax=vmax, clip=True)
stokes['Q']['norm'] = lambda vmin,vmax,linthresh : clr.SymLogNorm(vmin=vmin, linthresh=linthresh,
                                                                  vmax=vmax, clip=True)
stokes['U']['norm'] = lambda vmin,vmax,linthresh : clr.SymLogNorm(vmin=vmin, linthresh=linthresh,
                                                                  vmax=vmax, clip=True)
stokes['V']['norm'] = lambda vmin,vmax,linthresh : clr.SymLogNorm(vmin=vmin, linthresh=linthresh,
                                                                  vmax=vmax, clip=True)
stokes['I_scatter']['norm'] = lambda vmin,vmax,linthresh : clr.LogNorm(vmin=vmin,
                                                                       vmax=vmax, clip=True)
stokes['PI']['norm'] = lambda vmin,vmax,linthresh : clr.LogNorm(vmin=vmin,
                                                                vmax=vmax, clip=True)
stokes['P']['norm'] = lambda vmin,vmax,linthresh : clr.Normalize(vmin=vmin,
                                                                 vmax=vmax, clip=True)

R_map['I']['norm'] = lambda vmin,vmax,linthresh : clr.Normalize(vmin=vmin,
                                                                vmax=0.3, clip=True)
R_map['Q']['norm'] = lambda vmin,vmax,linthresh : clr.Normalize(vmin=vmin,
                                                                vmax=0.3, clip=True)
R_map['U']['norm'] = lambda vmin,vmax,linthresh : clr.Normalize(vmin=vmin,
                                                                vmax=0.3, clip=True)
R_map['V']['norm'] = lambda vmin,vmax,linthresh : clr.Normalize(vmin=vmin,
                                                                vmax=0.3, clip=True)

VOV_map['I']['norm'] = lambda vmin,vmax,linthresh : clr.Normalize(vmin=vmin,
                                                                  vmax=0.3, clip=True)
VOV_map['Q']['norm'] = lambda vmin,vmax,linthresh : clr.Normalize(vmin=vmin,
                                                                  vmax=0.3, clip=True)
VOV_map['U']['norm'] = lambda vmin,vmax,linthresh : clr.Normalize(vmin=vmin,
                                                                  vmax=0.3, clip=True)
VOV_map['V']['norm'] = lambda vmin,vmax,linthresh : clr.Normalize(vmin=vmin,
                                                                  vmax=0.3, clip=True)

N_photon_map['norm'] = lambda vmin,vmax,linthresh : clr.LogNorm(vmin=vmin,
                                                                vmax=vmax, clip=True)


#=========================================================================================================================#
#                                                                                                                         #
#=========================================================================================================================#


plot_quantities = [stokes['I'], stokes['Q'], stokes['U'], stokes['V'],
                   #stokes['I_scatter'],
                   stokes['PI'], stokes['P'], N_photon_map,
                   R_map['I'], R_map['Q'], R_map['U'],
                   #R_map['V'],
                   #VOV_map['I'], VOV_map['Q'], VOV_map['U'],
                   #VOV_map['V'],
                   ]

# init slider values
i_stokes = 0
wave_ind = 0

x_max = 200 * u.au

x_ref_vector_box = 0.02
y_ref_vector_box = 0.985
width_ref_vector_box = 0.132
height_ref_vector_box = -0.089


fig = plt.figure()
ax = plt.Subplot(fig, 111)
fig.add_axes(ax)
plt.subplots_adjust(bottom=0.2)


max_val = np.amax(plot_quantities[i_stokes]['data'][wave_ind,:,:])
min_val = np.amin(plot_quantities[i_stokes]['data'][wave_ind,:,:])

if i_stokes > 6:
    max_val = 0.2

vmax = max_val
vmin = plot_quantities[i_stokes]['vmin'](max_val,min_val)
linthresh = 10**np.floor(np.log10(vmax.value*1e-3)) * max_val.unit

map_ax = ax.imshow(plot_quantities[i_stokes]['data'][wave_ind,:,:].value,
                   extent=mrange.value,
                   interpolation='none',
                   origin='lower',
                   aspect='equal',
                   cmap=my_cmap_r)
map_ax.set_norm(plot_quantities[i_stokes]['norm'](vmin.value,vmax.value,linthresh.value))

ax.set_xlabel(r'$x$ [{}]'.format(mrange.unit))
ax.set_ylabel(r'$y$ [{}]'.format(mrange.unit))

ax.yaxis.labelpad = -10

ax.set_xlim(mrange[:2])
if mrange[:2][-1] > x_max:
    ax.set_xlim(-x_max,x_max)
ax.set_ylim(mrange[2:])
if mrange[2:][-1] > x_max:
    ax.set_ylim(-x_max,x_max)

text_wave_ax = ax.text(left-0.02, bottom-0.03
                       , r'$\lambda = \SI{%.1f}{\um}$' %(wave[wave_ind].to(u.micron).value)
                       , ha='left'
                       , va='bottom'
                       , transform=ax.transAxes
                       , color='k'
                       , backgroundcolor='w'
                       )

text_stokes_ax = ax.text(right, top, plot_quantities[i_stokes]['string'],
                         ha='right',
                         va='top',
                         transform=ax.transAxes,
                         color='k',
                         backgroundcolor='w')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.07)
clb = fig.colorbar(map_ax,cax=cax)
clb.set_label(plot_quantities[i_stokes]['clabel'])

if n_wave > 1:
    ax_wave = plt.axes([0.3, 0.09, 0.4, 0.02])
    slider_wave = Slider(ax_wave, r'wave$\_$ID', 1 , n_wave, valinit=1, valstep=1, valfmt='%i')

ax_stokes = plt.axes([0.3, 0.05, 0.4, 0.02])
slider_stokes = Slider(ax_stokes, r'Stokes', 1 , len(plot_quantities), valinit=1, valstep=1, valfmt='%i')


if args.radial_profile:
    fig_pol, ax_pol = plt.subplots()

    x_radial_profile, y_radial_profile, std = sm.radial_profile(plot_quantities[i_stokes]['data'][wave_ind,:,:].value,
                                                            center=((pixel[0].value-1)/2,(pixel[1].value-1)/2),
                                                            bin_size=args.bin_size)

    # convert x from pixel to physical units
    x_radial_profile *= (2*boundary_map['x'][args.unit].value / pixel[0].value)
    x_radial_profile = x_radial_profile * mrange.unit

    plot_p, _, err_bar = ax_pol.errorbar(x_radial_profile,y_radial_profile, yerr=std)

    ax_pol.set_ylabel(plot_quantities[i_stokes]['clabel'])
    ax_pol.set_xlim(0,min(x_max,boundary_map['x'][args.unit]))

    if i_stokes < 5:
        ax_pol.set_yscale('log')
        ax_pol.set_ylim(1e-4*np.nanmax(y_radial_profile),
                        1.1*np.nanmax(y_radial_profile))
    elif i_stokes == 6:
        ax_pol.set_yscale('log')
        ax_pol.set_ylim(1,
                        1.1*np.nanmax(y_radial_profile))
    else:
        ax_pol.set_yscale('linear')
        ax_pol.set_ylim(bottom=0)

    ax_pol.set_xlabel(r'$r$ [\si{\au}]')

    text_wave_pol_ax = ax_pol.text(0.45, top, r'$\lambda = \SI{%.1f}{\um}$' %(wave[wave_ind].to(u.micron).value), ha='left', va='top', transform=ax.transAxes, color='k')


if args.pol_vectors:
    # parameters for the vector arrow visualisation
    units='width'
    scale_units='width'
    width=0.004
    headwidth=0.
    headlength=0.
    headaxislength=0.
    pivot='middle'

    # plot the vectors
    ax.quiver(pol_vectors['x_pos'][wave_ind,:,:],
            pol_vectors['y_pos'][wave_ind,:,:],
            pol_vectors['x_length'][wave_ind,:,:],
            pol_vectors['y_length'][wave_ind,:,:],
            scale=vec_scale, units=units, width=width,
            scale_units=scale_units, headwidth=headwidth,
            headlength=headlength, headaxislength=headaxislength,
            pivot=pivot,color='b')

    # plot the reference vector with maximum degree of polarization
    ax.add_patch(matplotlib.patches.Rectangle((x_ref_vector_box,y_ref_vector_box)
                                            ,width_ref_vector_box
                                            ,height_ref_vector_box
                                            ,color='lightgrey'
                                            ,zorder=2
                                            ,transform=ax.transAxes
                                            ))

    ax.quiver(x_ref_vector_box+0.5*width_ref_vector_box,
            y_ref_vector_box+0.75*height_ref_vector_box,
            1/args.vec_stretch,
            0
            ,scale=vec_scale, units=units, width=width*2
            ,scale_units=scale_units, headwidth=headwidth
            ,headlength=headlength, headaxislength=headaxislength
            ,pivot=pivot
            ,transform=ax.transAxes
            ,color='b'
            ,zorder=4
            )

    ax.text(x_ref_vector_box+0.5*width_ref_vector_box
            ,y_ref_vector_box+0.6*height_ref_vector_box
            ,r'\SI{%.1f}{\percent}' %(np.amax(stokes['P']['data'][wave_ind,:,:])/args.vec_stretch)
            ,ha='center'
            ,va='bottom'
            ,color='k'
            ,backgroundcolor='lightgrey'
            ,transform=ax.transAxes
            ,zorder=3
            )

def update(val):
    if n_wave > 1:
        wave_ind = int(slider_wave.val) - 1
    else:
        wave_ind = 0
    i_stokes = int(slider_stokes.val) - 1

    if args.radial_profile:
        # update radial profile
        x_radial_profile, y_radial_profile, std = sm.radial_profile(plot_quantities[i_stokes]['data'][wave_ind,:,:].value,
                                                                center=((pixel[0].value-1)/2,(pixel[1].value-1)/2),
                                                                bin_size=args.bin_size)

        # convert x from pixel to physical units
        x_radial_profile *= (2*boundary_map['x'][args.unit].value / pixel[0].value)
        x_radial_profile = x_radial_profile * mrange.unit

        y_error = [np.array([[x, yt], [x,yb]]) for x, yt, yb in zip(x_radial_profile.value, y_radial_profile - std, y_radial_profile + std)]

        plot_p.set_ydata(y_radial_profile)
        err_bar[0].set_segments(y_error)

        # recompute the ax.dataLim
        ax_pol.relim()
        # update ax.viewLim using the new dataLim
        ax_pol.autoscale()
        ax_pol.set_xlim(0,min(x_max,boundary_map['x'][args.unit]))
        ax_pol.set_ylabel(plot_quantities[i_stokes]['clabel'])

        if i_stokes < 5:
            ax_pol.set_yscale('log')
            ax_pol.set_ylim(1e-4*np.nanmax(y_radial_profile),
                            1.1*np.nanmax(y_radial_profile))
        elif i_stokes == 6:
            ax_pol.set_yscale('log')
            ax_pol.set_ylim(1,
                            1.1*np.nanmax(y_radial_profile))
        else:
            ax_pol.set_yscale('linear')
            ax_pol.set_ylim(bottom=0)

    max_val = np.amax([np.abs(np.amax(plot_quantities[i_stokes]['data'][wave_ind,:,:].value)),
                      np.abs(np.amin(plot_quantities[i_stokes]['data'][wave_ind,:,:].value))])
    max_val = max_val * plot_quantities[i_stokes]['data'][wave_ind,:,:].unit

    if i_stokes > 6:
        max_val = 0.2 * u.dimensionless_unscaled
    min_val = -max_val

    vmax = max_val
    vmin = plot_quantities[i_stokes]['vmin'](max_val,min_val)
    linthresh = 10**np.floor(np.log10(vmax.value*1e-3)) * max_val.unit

    # clear the figure
    global ax, cax
    cax.remove()
    ax.remove()

    # create new axis on figure
    ax = plt.Subplot(fig, 111)
    fig.add_axes(ax)
    plt.subplots_adjust(bottom=0.2)

    # update map
    map_ax = ax.imshow(plot_quantities[i_stokes]['data'][wave_ind,:,:].value,
                   extent=mrange.value,
                   interpolation='none',
                   origin='lower',
                   aspect='equal',
                   cmap=my_cmap_r)
    map_ax.set_norm(plot_quantities[i_stokes]['norm'](vmin.value,vmax.value,linthresh.value))

    ax.set_xlabel(r'$x$ [{}]'.format(mrange.unit))
    ax.set_ylabel(r'$y$ [{}]'.format(mrange.unit))

    ax.yaxis.labelpad = -10

    ax.set_xlim(mrange[:2])
    if mrange[:2][-1] > x_max:
        ax.set_xlim(-x_max,x_max)
    ax.set_ylim(mrange[2:])
    if mrange[2:][-1] > x_max:
        ax.set_ylim(-x_max,x_max)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.07)
    clb = fig.colorbar(map_ax,cax=cax)
    clb.set_label(plot_quantities[i_stokes]['clabel'])


    text_wave_ax = ax.text(left-0.02, bottom-0.03
                            , r'$\lambda = \SI{%.1f}{\um}$' %(wave[wave_ind].to(u.micron).value)
                            , ha='left'
                            , va='bottom'
                            , transform=ax.transAxes
                            , color='k'
                            , backgroundcolor='w'
                            )

    text_stokes_ax = ax.text(right, top, plot_quantities[i_stokes]['string'],
                            ha='right',
                            va='top',
                            transform=ax.transAxes,
                            color='k',
                            backgroundcolor='w')

    if args.radial_profile:
        text_wave_pol_ax.set_text(r'$\lambda = \SI{%.1f}{\um}$' %(wave[wave_ind].to(u.micron).value))


    if args.pol_vectors:
        # update pol. vectors
        ax.quiver(pol_vectors['x_pos'][wave_ind,:,:],
                pol_vectors['y_pos'][wave_ind,:,:],
                pol_vectors['x_length'][wave_ind,:,:],
                pol_vectors['y_length'][wave_ind,:,:],
                scale=vec_scale, units=units, width=width,
                scale_units=scale_units, headwidth=headwidth,
                headlength=headlength, headaxislength=headaxislength,
                pivot=pivot,color='b')

        # plot the reference vector with maximum degree of polarization
        ax.add_patch(matplotlib.patches.Rectangle((x_ref_vector_box,y_ref_vector_box)
                                                ,width_ref_vector_box
                                                ,height_ref_vector_box
                                                ,color='lightgrey'
                                                ,zorder=2
                                                ,transform=ax.transAxes
                                                ))

        ax.quiver(x_ref_vector_box+0.5*width_ref_vector_box,
                y_ref_vector_box+0.75*height_ref_vector_box,
                vec_scale/pixel[0],
                0
                ,scale=vec_scale, units=units, width=width*2
                ,scale_units=scale_units, headwidth=headwidth
                ,headlength=headlength, headaxislength=headaxislength
                ,pivot=pivot
                ,transform=ax.transAxes
                ,color='b'
                ,zorder=4
                )

        ax.text(x_ref_vector_box+0.5*width_ref_vector_box
                ,y_ref_vector_box+0.6*height_ref_vector_box
                ,r'\SI{%.1f}{\percent}' %(np.amax(stokes['P'][wave_ind,:,:])*args.vec_stretch)
                ,ha='center'
                ,va='bottom'
                ,color='k'
                ,backgroundcolor='lightgrey'
                ,transform=ax.transAxes
                ,zorder=3
                )

    # update the figure to make changes visible
    fig.canvas.draw_idle()
    if args.radial_profile:
        fig_pol.canvas.draw_idle()

if n_wave > 1:
    slider_wave.on_changed(update)
slider_stokes.on_changed(update)

plt.show()
