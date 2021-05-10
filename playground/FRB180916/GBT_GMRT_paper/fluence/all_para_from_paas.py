import numpy as np
from glob import glob
import psrchive
import sys

def von_mises(x, p):
  #To open .m files from paas
  #mu, k and A are the three parameters in the fil
  #x is defined in phase
  mu, k, A = p
  vm = np.exp(k*np.cos((x-mu)*2.*np.pi))/2./np.pi/np.i0(k)
  return vm / vm.max() * A

def von_mises_approx(x, p):
  #Approximate von Mises function with normal distribution
  #x is defined in phase
  mu, k, A = p
  s2 = 1./k
  vm = 1./(2*s2*np.pi)**0.5*np.exp(-((x-mu)*2*np.pi)**2/2/s2) 
  return vm / vm.max() * A

def curve(x, m_name):
  m = np.loadtxt(m_name)
  if len(m.shape) == 1:  y = von_mises_approx(x, m)
  else:
    y = np.zeros_like(x)
    for mi in m:
      y += von_mises_approx(x, mi)
  return y

def get_period(ar):
  std = ar[:-1] + 'std'
  load_archive = psrchive.Archive_load(std)  
  return (load_archive.end_time() - load_archive.start_time()).in_seconds() * 1000  #ms


x = np.linspace(0,1,1e6)

ar_list = glob('*.m')

for ar in ar_list:
  y = curve(x, ar)
  HM = x[y >= y.max()/2.]
  period = get_period(ar)
  FWHM = (HM.max() - HM.min()) * period
  F = (FWHM*y.max())/(2.35*0.3989*1000)
  print "{:.30}: w = {:.3f} ms F = {:.3f} Jy ms".format(ar, FWHM, F)

