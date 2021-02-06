
import sys; sys.path.append('../py_src')
from pyOM_gui import pyOM_gui as pyOM
from numpy import *

BETA=2e-11

class enso1(pyOM):
   """ Enso response
   """
   def set_parameter(self):
     """set main parameter
     """
     M=self.fortran.main_module   
     M.nx    = 64*2
     M.nz    = 1
     M.ny    = 64*2
     M.dt_mom    = 3600.0 /2.
     M.dt_tracer = 3600.0 /2.

     M.congr_epsilon = 1e-12
     M.enable_streamfunction = 0
     M.enable_free_surface   = 1
     M.eq_of_state_type      = 1
     
     M.enable_conserve_energy = 0
     M.coord_degree           = 0
     
     M.enable_biharmonic_friction  = 1
     M.a_hbi  = .05e11
     
     M.enable_hydrostatic          = 1
     M.enable_cyclic_x             = 0
     M.enable_superbee_advection   = 1
     #M.enable_quicker_mom_advection= 0
     #M.enable_no_mom_advection     = 1
     
     M.enable_momentum_sources = 1
     return

   def set_grid(self):
     M=self.fortran.main_module   
     M.dxt[:]    = 50e3
     M.dyt[:]    = 50e3
     M.dzt[:] = (500e3**2*BETA) **2 /9.81
     return
 

   def set_coriolis(self):
     """ vertical and horizontal Coriolis parameter on yt grid
         routine is called after initialization of grid
     """
     M=self.fortran.main_module   
     y0=M.ny*M.dxt[0]*0.5
     for j in range( M.yt.shape[0] ): M.coriolis_t[:,j]   =  BETA*(M.yt[j]-y0)
     return

   def set_forcing(self):
     M=self.fortran.main_module  
     t = M.dt_tracer*M.itt 
     if t/86400.>3.0: M.surface_taux[:] = 0.
     
     y1=10*M.dyt[0]
     y2=(M.ny-10)*M.dyt[0]
     dy = 2*y1
     lam = 1/(3600.*10)
     for j in range(M.js_pe,M.je_pe+1):
      #print j,M.yt[j],y1,exp( -(M.yt[j]-y1)**2/dy**2)
      M.u_source[:,j,:]= -lam*M.u[:,j,:,M.tau-1]*(exp(-(M.yt[j]-y1)**2/dy**2) +
                                                  exp(-(M.yt[j]-y2)**2/dy**2) )
      M.v_source[:,j,:]= -lam*M.v[:,j,:,M.tau-1]*(exp(-(M.yt[j]-y1)**2/dy**2)  +
                                                  exp(-(M.yt[j]-y2)**2/dy**2) )
     #stop                                                                                                  
     return
     
   def set_initial_conditions(self):
     """ setup all initial conditions
     """
     M=self.fortran.main_module  
     M.u[:]=0.
     M.v[:]=0.
     M.w[:]=0. 
     M.psi[:] = 0.
     M.itt = 0
     cn =  (M.dzt[0]*9.81)**0.5  
     hn=cn**2/9.81
     Re = (cn /BETA)**0.5  
     y0=M.ny*M.dxt[0]*0.5
     g=9.81
     print("Eq. Rossby radius = ",Re/1e3," km")
     print("c_n = ",cn," m/s")
     print("h_n = ",hn," m" )
     
     #for i in range(M.xt.shape[0]):
     #  for j in range(M.yt.shape[0]):
     #    M.psi[i,j,:]=0.1*exp( -(M.xt[i]-y0*0.5)**2/(0.5*Re)**2 -(M.yt[j]-y0)**2/(0.5*Re)**2 )
     
     # wind stress forcing
     x0 = 1.5e6
     for i in range(M.is_pe,M.ie_pe+1):
          ii = self.if2py(i)
          taux=0.0
          #if  M.xt[ii]>1e6 and M.xt[ii]<2e6: 
          M.surface_taux[ii,:] =  .1e-3*exp( -(M.yt-y0)**2/(3*Re)**2 )*exp(-(M.xt[ii]-x0)**2/0.5e6**2 )
             
     return

   def make_plot(self):
     """ make a plot using methods of self.figure
     """
     if hasattr(self,'figure'):
       M=self.fortran.main_module         # fortran module with model variables
       x=M.xt[2:-2]/1e3
       y=(M.yt[2:-2]-M.ny/2*M.dyt[0])/1e3

       self.figure.clf()
       ax=self.figure.add_subplot(111)
       a=M.psi[2:-2,2:-2,M.tau-1] 
       co=ax.contourf(x,y,minimum(10,a.transpose()),arange(-10,12,1),cmap='jet')
       
       a=M.u[2:-2:4,2:-2:4,0,M.tau-1] 
       b=M.v[2:-2:4,2:-2:4,0,M.tau-1] 
       ax.quiver(x[::4],y[::4],a.transpose(),b.transpose(),angles='xy', scale_units='xy', scale=2e-3 )
       
       t = M.dt_tracer*M.itt 
       if t/86400.<3.0: 
         ax.contour(x,y,M.surface_taux[2:-2,2:-2].transpose(),colors='k')
       
       self.figure.colorbar(co)
       ax.set_title('thermocline depth at t='+str(M.dt_tracer*M.itt /86400.)+" days")
       ax.set_xlabel('x [km]')
       ax.set_ylabel('y [km]')
     return

if __name__ == "__main__":
   model= enso1()
   model.run(snapint=0.5*86400.0,runlen=365*86400.)
   
