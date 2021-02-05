
# ##############################################################
# 1D case for tke model
# ############################################################## 

import sys; sys.path.append('../py_src')
from pyOM_cdf import pyOM_cdf as pyOM
import numpy as np

N_0     = 0.02

class test_tke1(pyOM):

 def set_parameter(self):
	
    M=self.fortran.main_module
    M.nx = 1
    M.ny = 1
    M.nz = 50 
    M.snapint = 86400 /6
    M.dt_mom = 360 # time step in second for momentum
    M.dt_tracer = 360 # time step for tracer

    M.enable_momentum_equation = 0
    M.enable_thermodynamic_equation = 1
    M.enable_conserve_energy = 0
    M.coord_degree  = 0
    M.enable_cyclic_x = 1
    M.enable_cyclic_y = 1
    M.eq_of_state_type = 5 # 1: linear, 3: non-linear, 5:TEOS ? 
    M.enable_diag_snapshots  = 1

    M.enable_implicit_vert_friction = 1
    T=self.fortran.tke_module 
    T.enable_tke = 1
    T.alpha_tke = 30.0
    T.tke_mxl_choice = 2
    return     

 def set_grid(self):
	
  M=self.fortran.main_module
  M.dxt[:] = 5e3
  M.dyt[:] = 5e3
  M.dzt[:] = 1.0
  return
 
 def set_coriolis(self):
  return

 def set_topography(self):
   M=self.fortran.main_module 
   M.kbot[:]=0
   return  

 def set_initial_conditions(self):
	
  M=self.fortran.main_module
  M.salt[:]= 35.

  for k in range(M.nz):
     alpha = self.fortran.gsw_eq_of_state.gsw_drhodt(35.,20.,0.)
     M.temp[:,:,k,:]  = 32+M.rho_0/M.grav/alpha*(-N_0**2*M.zt[k])  

  T = self.fortran.tke_module
  T.tke[:] = 1e-12
  T.forc_tke_surface[:] = np.sqrt( (0.1e-3)**2 )**(3./2.) 
  return

 def set_forcing(self):
  M=self.fortran.main_module
  return	
 
 def set_diagnostics(self):
    return

if __name__ == "__main__": 
	model = test_tke1()  
	model.run( snapint = 8640.0, runlen = 365*86400.)


 


