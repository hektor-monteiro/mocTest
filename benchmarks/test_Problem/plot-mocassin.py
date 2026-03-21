import numpy as np
import matplotlib.pyplot as plt
import sys, os

from matplotlib.colors import LogNorm
from astropy.io import fits
from astropy import units as u
from astropy import wcs
from astropy.visualization import make_lupton_rgb
import scipy.ndimage as ndimage
from scipy.signal import savgol_filter, medfilt
import warnings
warnings.filterwarnings('ignore')

plt.close('all')
############################################################
# Read input file from mocassin

print ('Reading in input parameters from input.in ...')
moc_in = open("input/input.in", "r")
lines = moc_in.readlines()
moc_pars = []
for line in lines:
    if (line.split() != []):
        moc_pars.append(line.split())
moc_in.close()    

############################################################
# Reading Grid Structure
        
print ('Reading Grid Structure ...')
grid0 = open("output/grid0.out", "r")
lines = grid0.readlines()

nx,ny,nz, aux = [int(x) for x in (lines[1].split())[0:4]]
ncell = nx*ny*nz

aux,Rout = [float(x) for x in (lines[1].split())[4:]]
X = np.array([float(x) for x in lines[2:2+nx]])
Y = np.array([float(x) for x in lines[2+nx:2+nx+ny]])
Z = np.array([float(x) for x in lines[2+nx+ny:2+nx+ny+nz]])

aux = np.array([[int(y) for y in x.split()] for x in lines[2+nx+ny+nz:]])
converged = np.reshape(aux[:,1],(nx,ny,nz))

plt.figure()
plt.imshow(converged[:,int(nx/2),:], origin='lower')
plt.title('Cell convergente')

if (nx == ny == nz or Rout == 0):
    Rout = np.abs(Z).max() 
############################################################
# reading  ionic fractions
Nelements =9

H_ionic_frac = np.zeros((2,ncell))
He_ionic_frac = np.zeros((3,ncell))
Heavy_ionic_frac = np.zeros((Nelements-2,7,ncell))

grid2 = open("output/grid2.out", "r")
aux = grid2.readlines()

for i in range(ncell):
    H_ionic_frac[:,i] = aux[i*Nelements].split()
    He_ionic_frac[:,i] = aux[i*Nelements+1].split()
    Heavy_ionic_frac[:,:,i] = np.array([x.split() for x in aux[i*Nelements+2:i*Nelements+2+(Nelements-2)]])
    
# np.save('OI_ionization_strut_51',o16300cube)
# np.save('HI_frac_strut_51',np.reshape(H_ionic_frac[0,:],(nx,ny,nz)))

############################################################
# reading  electron temperatures, electron densities and hydrogen densities

print ('Reading Ne and Te Structure ...')
grid1 = np.loadtxt('output/grid1.out')

Te = np.reshape(grid1[:,0],(nx,ny,nz))
Ne = np.reshape(grid1[:,1],(nx,ny,nz))
H0 = np.reshape(grid1[:,2],(nx,ny,nz))

print('Peak density in H0: ',H0.max())

f,ax = plt.subplots(1,3,figsize=(11,3))
p1 = ax[0].imshow(H0.mean(axis=0))
f.colorbar(p1, ax=ax[0])
ax[0].set_title('H density ($cm^{-3}$)')
ax[0].axis('off')
p2 = ax[1].imshow(Ne.mean(axis=0))
f.colorbar(p2, ax=ax[1])
ax[1].set_title('Ne ($cm^{-3}$)')
ax[1].axis('off')
p3 = ax[2].imshow(Te.mean(axis=0))
f.colorbar(p3, ax=ax[2])
ax[2].set_title('Te (K)')
ax[2].axis('off')
plt.tight_layout()
plt.savefig('figs/projected_H0_Ne_Te_mod.png', dpi=300)


f,ax = plt.subplots(1,3,figsize=(11,3))
p1 = ax[0].imshow(H0[:,int(nx/2),:], origin='lower',norm=LogNorm(vmin=50))
f.colorbar(p1, ax=ax[0])
ax[0].set_title('H density ($cm^{-3}$)')
ax[0].axis('off')
p2 = ax[1].imshow(Ne[:,int(nx/2),:], origin='lower',norm=LogNorm(vmin=50))
f.colorbar(p2, ax=ax[1])
#ax[1].contour(H0,colors='w')
ax[1].set_title('Ne ($cm^{-3}$)')
ax[1].axis('off')
p3 = ax[2].imshow(Te[:,int(nx/2),:], origin='lower',norm=LogNorm(vmin=5000))
f.colorbar(p3, ax=ax[2])
ax[2].set_title('Te (K)')
ax[2].axis('off')
plt.tight_layout()
plt.savefig('figs/Zcut_H0_Ne_Te_mod.png', dpi=300)

############################################################
# Read input file from mocassin

print ('Reading in lineflux.out file ...')
lineflux = open("output/lineFlux.out", "r")
lines = lineflux.readlines()
line_flux_out = []
for line in lines:
    if (line.split() != []):
        line_flux_out.append(line.split())
lineflux.close()    

hbeta = float(line_flux_out[2][3])

############################################################
# Other data needed
dist = 1170.

mod_hbeta = (hbeta*1.e36)/4./np.pi/(3.086e18*dist)**2
print ('Using distance of d =%5.0f parsecs.'%dist)
print ('Model Hbeta at earth: %6.2e   Lit.: 7.76e-10 '%mod_hbeta)
 
############################################################
# observed size for given distance
arcsec_sz = 3600.*180./np.pi*np.arctan(Rout/(dist*3.086e18))

############################################################
# Read SED file
sed = np.genfromtxt("output/SED.out",skip_header=4,skip_footer=4,dtype=None)
sed[:,3:] = np.where(sed[:,3:] >= 1e-8, sed[:,3:], np.nan)

sed_ori = np.genfromtxt("output/SED_MIE1000.out",skip_header=4,skip_footer=4,dtype=None)
sed_ori[:,3:] = np.where(sed_ori[:,3:] >= 1e-8, sed_ori[:,3:], np.nan)

sed_cde = np.genfromtxt("output/SED_CDE1000.out",skip_header=4,skip_footer=4,dtype=None)
sed_cde[:,3:] = np.where(sed_cde[:,3:] >= 1e-8, sed_cde[:,3:], np.nan)

############################################################
# Read ionization source file
#%%

from scipy import integrate
from astropy.constants import sigma_sb, M_sun, L_sun, M_earth
from astropy.modeling.models import BlackBody

sigma = sigma_sb.to(u.erg / u.s / u.cm**2 / u.K**4)          # Stefan-Boltzmann constant

ion_file_ind = [i for i, val in enumerate(moc_pars) if 'contShape' in val][0]

Lstar = float(moc_pars[[i for i, val in enumerate(moc_pars) if 'LStar' in val][0]][1]) * 1.e36 * u.erg / u.s 
Tstar = float(moc_pars[[i for i, val in enumerate(moc_pars) if 'TStellar' in val][0]][1]) * u.K
Rstar = np.sqrt(Lstar / (4*np.pi*sigma*Tstar**4))
Fstar = Lstar / u.cm**2 #(4*np.pi*Rstar**2) / (4*np.pi)

print("Central source with Teff= %8.3f K and L= %8.3f Lsun"%(Tstar.value,Lstar/L_sun.to(u.erg/u.s)))
bb = BlackBody(temperature=Tstar)
bb_flux = (bb(sed[:,1]*u.micron)).to(u.erg/u.cm**2/u.s/u.sr/u.AA, equivalencies=u.spectral_density(sed[:,1]*u.micron)) * 4 * np.pi * u.sr

bb_CS = BlackBody(temperature=100*u.K)
bb_CS_flux = (bb_CS(sed[:,1]*u.micron)).to(u.erg/u.cm**2/u.s/u.sr/u.AA, equivalencies=u.spectral_density(sed[:,1]*u.micron)) * 4 * np.pi * u.sr
           
if (moc_pars[ion_file_ind][1].replace('"','') == 'blackbody'):
    bb = BlackBody(temperature=Tstar)
    ion_source_flux = (bb(sed[:,1]*u.micron)).to(u.erg/u.cm**2/u.s/u.sr/u.AA, equivalencies=u.spectral_density(sed[:,1]*u.micron)) * 4 * np.pi * u.sr
    ion_source_wave = sed[:,1]*u.micron
else:
    ion_source = np.genfromtxt(moc_pars[ion_file_ind][1].replace('"',''),skip_header=29,skip_footer=0,dtype=None)
    ion_source_wave = (ion_source[:,0] * u.AA).to(u.micron)
    ion_source_flux = ion_source[:,1] * u.erg / u.s / u.cm**2 / u.AA

# scale to luminosity
ion_source_flux = ion_source_flux / np.abs(integrate.trapezoid(ion_source_flux, ion_source_wave.to(u.AA))) * Fstar
bb_flux = bb_flux / np.abs(integrate.trapezoid(bb_flux, (sed[:,1] * u.micron).to(u.AA))) * Fstar
bb_CS_flux = bb_CS_flux / np.abs(integrate.trapezoid(bb_CS_flux, (sed[:,1] * u.micron).to(u.AA))) * Fstar

# Dilute to same units in Mocassin
# user must still divide by D^2 in pc to get observed flux
ion_source_flux = ion_source_flux / (4*np.pi*(3.08e18)**2)  / dist**2 
bb_flux = bb_flux / (4*np.pi*(3.08e18)**2)  / dist**2 
bb_CS_flux = bb_CS_flux / (4*np.pi*(3.08e18)**2)  / dist**2 


# convert to Jy
ion_source_flux = ion_source_flux.to(u.Jy, equivalencies=u.spectral_density(ion_source_wave))
bb_flux = bb_flux.to(u.Jy, equivalencies=u.spectral_density(sed[:,1]*u.micron))
bb_CS_flux = bb_CS_flux.to(u.Jy, equivalencies=u.spectral_density(sed[:,1]*u.micron))

ion_source_flux = 0.8*ion_source_flux # correction for limited wavelength range and resolution

###############################################################################
# Plot escaped SED
#%%

rmax = (Rout*u.cm).to(u.pc)

cor_fac = np.nansum(sed[:200,3])/np.nansum(sed[:200,2])

plt.figure()
plt.xscale('log')
plt.yscale('log')
plt.ylim(1.0e-2,1.5e3)
plt.xlim(0.015,500000)
plt.xlabel('$\lambda$ (micron)')
plt.ylabel(r'$F_{\nu}$ (Jy)')

plt.plot(ion_source_wave, ion_source_flux, '--',label='Central source',color='C1',zorder=0)
#plt.plot(sed[:,1],sed[:,3]/(dist)**2,'C1',label='x axis SED',zorder=0,lw=2, alpha=0.5)

plt.plot(sed[:,1]*u.micron, bb_flux, '--',label='Blackbody',color='C2',zorder=0)

plt.plot(sed[:,1]*u.micron,sed[:,2]/ dist**2,'C0',label='avg. all angles',zorder=0,lw=2)
plt.plot(sed_ori[:,1]*u.micron,sed_ori[:,2]/ dist**2,'C0:',label='Moc. Original',zorder=0,lw=1)
plt.plot(sed_cde[:,1]*u.micron,sed_cde[:,2]/ dist**2,'C1--',label='Moc. CDE',zorder=1,lw=1)

plt.plot(sed[:,1],sed[:,3]*(1/dist)**2 / (0.154783085*4*np.pi),'C1',label='line of sight SED',zorder=5,lw=2)

plt.plot(sed[:,1],sed[:,4]*(1/dist)**2   ,'C3',label='CS SED',zorder=5,lw=2)



plt.legend()
# plt.scatter(iso_wave_s,iso_flux_s,label='ISO spectrum',alpha=0.25,s=1,c='k')

plt.tight_layout()
plt.savefig('figs/model_escaped_SED.png', dpi=300)

############################################################
# reading  line emissivities
#%%

# determine arcsec per cell
arcsec_cell = arcsec_sz/(nz/2)

print ('Reading line emissivities ...')
lines = np.loadtxt('input/plot.in', skiprows=1, usecols=2)
emiss = np.loadtxt('output/plot.out')

# generate line maps  and save fits files
plot_lines = [4861., 5008., 6585., 4686., 6718.,1906.]

fig=plt.figure()
cont = 1
line_maps = []

for i in range(len(lines)):
#    print ('Doing line %5.0f '%lines[i], i)
    cube = np.reshape(emiss[:,i+1],(nx,ny,nz))
    if(lines[i] == 6563.):
        hacube = cube
    if(lines[i] == 4861.):
        hbcube = cube
    if(lines[i] == 6585.):
        n26583cube = cube
    if(lines[i] == 5756.):
        n25756cube = cube
    if(lines[i] == 9072.):
        s39072cube = cube
    if(lines[i] == 4364.):
        o34363cube = cube
    if(lines[i] == 6302.):
        o16300cube = cube
    line_maps.append(np.nansum(cube,axis=0))

    hdu = fits.PrimaryHDU(np.flip(np.nansum(cube,axis=0),axis=0))
    hdulist = fits.HDUList([hdu])
        
    hdulist[0].header['Lambda'] = str(lines[i])
    hdulist[0].header['CTYPE1'] = 'RA---TAN'
    hdulist[0].header['CTYPE2'] = 'DEC--TAN'
    hdulist[0].header['CDELT1'] = (arcsec_cell*u.arcsec).to(u.degree).value 
    hdulist[0].header['CDELT2'] = (arcsec_cell*u.arcsec).to(u.degree).value
    hdulist[0].header['CUNIT1'] = 'deg'
    hdulist[0].header['CUNIT2'] = 'deg'
    hdulist[0].header['CRVAL1'] =   258.43500          
    hdulist[0].header['CRVAL2'] =   -37.10325        
    hdulist[0].header['CRPIX1'] = int(nz/2)+1
    hdulist[0].header['CRPIX2'] = int(ny/2)+1
    
    hdulist.writeto('maps/model_'+str(lines[i])+'.fits',overwrite=True, output_verify='fix')


###############################################################################
# Plot line maps
plot_lines = [4861., 5008., 6585., 4686., 6718.,1906.]

plt.figure()
cont = 1
for line in plot_lines:
    ax = plt.subplot(2,3,cont)
    ind = np.where(lines == line)
    ax.imshow((line_maps[ind[0][0]]+1.0e-16)/line_maps[ind[0][0]].max(),norm=LogNorm(vmin=1.e-6,vmax=0.99), origin='lower')
    ax.axis('off')
    #ax.imshow((line_maps[ind[0][0]]+1.0e-16)/line_maps[ind[0][0]].max(), vmin=1.e-8,vmax=0.6)
    ax.set_title(str(line))
    cont += 1
plt.tight_layout()
plt.savefig('figs/model_linemaps-log.png', dpi=300)


plt.figure()
cont = 1
for line in plot_lines:
    ax = plt.subplot(2,3,cont)
    ind = np.where(lines == line)
    ax.axis('off')
    ax.imshow((line_maps[ind[0][0]]+1.0e-16)/line_maps[ind[0][0]].max(), vmin=0.,vmax=0.99, origin='lower')
    ax.set_title(str(line))
    cont += 1
plt.tight_layout()
plt.savefig('figs/model_linemaps.png', dpi=300)


#########################################################################
# Plot tau values for each axis
#########################################################################
#%%

# File name (replace with your actual file path)
file_name = 'output/tauNu.out'

# Initialize data storage
data = {'1,0,0': {'lambda': [], 'tau': []},
        '0,1,0': {'lambda': [], 'tau': []},
        '0,0,1': {'lambda': [], 'tau': []}}

# Read the file
with open(file_name, 'r') as f:
    current_dir = None
    for line in f:
        line = line.strip()

        # Detect the direction
        if line.startswith('direction:'):
            current_dir = line.split(':')[-1].strip()
        
        # Read the data lines
        elif current_dir and line and not line.startswith('lambda'):
            parts = line.split()
            if len(parts) == 2:
                wavelength = float(parts[0])
                tau = float(parts[1])
                data[current_dir]['lambda'].append(wavelength)
                data[current_dir]['tau'].append(tau)

# Plotting
fig, ax = plt.subplots()
for direction, values in data.items():
    ax.plot(values['lambda'][1:], values['tau'][1:], label=f'Direction {direction}')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Wavelength [μm]')
ax.set_ylabel('Tau (τ)')
ax.set_title('Optical Depth from the Center to the Edge of the Nebula')
ax.legend()
plt.tight_layout()
plt.savefig('figs/model_axisTau.png', dpi=300)





plt.show()
