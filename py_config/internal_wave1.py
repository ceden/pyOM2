
import sys; sys.path.append('../py_src')
from pyOM_gui import pyOM_gui as pyOM
from numpy import *

fac=1
N_0 = 2*pi/10.
OM0 = 1./(1.5*10)

class internal_wave1(pyOM):
   
   """ internal wave maker
   """
   
   def set_parameter(self):
     """set main parameter
     """
     M=self.fortran.main_module
     (M.nx,M.ny,M.nz) = (64*fac,1,64*fac)
     
     M.dt_mom   =20*0.025/fac
     M.dt_tracer=20*0.025/fac
     
     M.enable_conserve_energy = 0
     M.coord_degree           = 0
     M.enable_cyclic_x        = 1
     M.enable_hydrostatic     = 0
     M.eq_of_state_type       = 1

     M.congr_epsilon = 1e-12
     M.congr_max_iterations = 5000
     M.congr_epsilon_non_hydro=   1e-9
     M.congr_max_itts_non_hydro = 5000    

     M.enable_explicit_vert_friction = 1
     M.kappam_0 = 5e-3/fac**2   
     M.enable_hor_friction = 1
     M.a_h      = 5e-3/fac**2
     M.enable_superbee_advection = 1

     M.enable_tempsalt_sources = 1
     M.enable_momentum_sources = 1
     return

   def set_grid(self):
       M=self.fortran.main_module   
       M.dxt[:]= 0.25/fac
       M.dyt[:]= 0.25/fac
       M.dzt[:]= 0.25/fac
       return

   def set_topography(self):
       M=self.fortran.main_module   
       M.kbot[:]=0
       M.kbot[:,2:-2]=1
       return


   def set_initial_conditions(self):
     """ setup all initial conditions
     """
     M=self.fortran.main_module
     alpha = self.fortran.linear_eq_of_state.linear_eq_of_state_drhodt()
     grav = 9.81; rho0 = 1024.0

     self.t0  = zeros( (M.i_blk+2*M.onx,M.j_blk+2*M.onx,M.nz) , 'd', order='F')
     self.dt0 = zeros( (M.i_blk+2*M.onx,M.j_blk+2*M.onx,M.nz,3) , 'd', order='F')
     self.u0  = zeros( (M.i_blk+2*M.onx,M.j_blk+2*M.onx,M.nz) , 'd', order='F')

     # prescribed background stratification
     for k in range(M.nz):
        self.t0[:,:,k]=-N_0**2*M.zt[k]/grav/alpha*rho0*M.maskt[:,:,k]

     # wave maker
     for k in range(M.nz):
        self.u0[:,:,k]= M.masku[:,:,k]*1./(100*60.*M.dt_tracer)*exp( -(M.zt[k]-M.zt[M.nz/2-1])**2/(M.dzt[0]*1)**2 )
     # find x for nx/2
     x0  = zeros( (1,) , 'd', order='F')
     if M.nx/2 >= M.is_pe and M.nx/2<= M.ie_pe:
         x0 = M.xu[self.if2py(M.nx/2)]
     self.fortran.global_max(x0)
     for i in range( self.u0.shape[0]):
        self.u0[i,:,:]= self.u0[i,:,:]*exp( -(M.xu[i]-x0)**2/(M.dxu[0]*1)**2 )
     return

   def set_forcing(self):
      M=self.fortran.main_module   
      # implement effect of background state 
      # update density, etc of last time step
      M.temp[:,:,:,M.tau-1] = M.temp[:,:,:,M.tau-1] + self.t0
      self.fortran.calc_eq_of_state(M.tau)
      M.temp[:,:,:,M.tau-1] = M.temp[:,:,:,M.tau-1] - self.t0
      
      # advection of background temperature
      self.fortran.advect_tracer(M.is_pe-M.onx,M.ie_pe+M.onx,M.js_pe-M.onx,M.je_pe+M.onx,self.t0,self.dt0[...,M.tau-1],M.nz) 
      M.temp_source[:] = (1.5+M.ab_eps)*self.dt0[...,M.tau-1] - ( 0.5+M.ab_eps)*self.dt0[...,M.taum1-1]

      # wave maker   
      M.u_source[:,:,:]= self.u0[:,:,:]*sin(2*pi*OM0*M.itt*M.dt_tracer)
      return

   def user_defined_signal(self):
       """ this routine must be called by all processors
       """
       M=self.fortran.main_module  
       a = zeros( (M.nx,M.ny), 'd', order = 'F')
       a[M.is_pe-1:M.ie_pe,0] = M.xt[2:-2]
       self.fortran.pe0_recv_2d(a)
       self.xt_gl = a[:,0].copy()
       
       self.temp_gl = zeros( (M.nx,M.ny,M.nz), 'd', order = 'F')
       self.t0_gl   = zeros( (M.nx,M.ny,M.nz), 'd', order = 'F')
       self.u_gl    = zeros( (M.nx,M.ny,M.nz), 'd', order = 'F')
       self.w_gl    = zeros( (M.nx,M.ny,M.nz), 'd', order = 'F')
       for k in range(M.nz):
         a[M.is_pe-1:M.ie_pe,M.js_pe-1:M.je_pe] = where( M.maskt[2:-2,2:-2,k] >0,  M.temp[2:-2,2:-2,k,M.tau-1] , NaN) 
         self.fortran.pe0_recv_2d(a)
         self.temp_gl[:,:,k]=a.copy()

         a[M.is_pe-1:M.ie_pe,M.js_pe-1:M.je_pe] = where( M.maskt[2:-2,2:-2,k] >0,  self.t0[2:-2,2:-2,k] , NaN) 
         self.fortran.pe0_recv_2d(a)
         self.t0_gl[:,:,k]=a.copy()

         a[M.is_pe-1:M.ie_pe,M.js_pe-1:M.je_pe] = where( M.masku[2:-2,2:-2,k] >0,  M.u[2:-2,2:-2,k,M.tau-1] , NaN) 
         self.fortran.pe0_recv_2d(a)
         self.u_gl[:,:,k]=a.copy()
        
         a[M.is_pe-1:M.ie_pe,M.js_pe-1:M.je_pe] = where( M.maskw[2:-2,2:-2,k] >0,  M.w[2:-2,2:-2,k,M.tau-1] , NaN) 
         self.fortran.pe0_recv_2d(a)
         self.w_gl[:,:,k]=a.copy()
       
  
   def make_plot(self):
       M=self.fortran.main_module
       
       self.set_signal('user_defined') # following routine is called by all PEs
       self.user_defined_signal()
       
       self.figure.clf()
       ax=self.figure.add_subplot(111)
       t = self.temp_gl[:,0,:]
       u = self.u_gl[:,0,:]
       w = self.w_gl[:,0,:]
       t[M.nx/2-3*fac:M.nx/2+2*fac,M.nz/2-3*fac:M.nz/2+2*fac]=0
       u[M.nx/2-3*fac:M.nx/2+2*fac,M.nz/2-3*fac:M.nz/2+2*fac]=0
       w[M.nx/2-3*fac:M.nx/2+2*fac,M.nz/2-3*fac:M.nz/2+2*fac]=0
       co=ax.contourf(self.xt_gl,M.zt,t.transpose() )
       self.figure.colorbar(co)
       ax.contour(self.xt_gl,M.zt,(t+self.t0_gl[:,0,:]).transpose(),10,colors='red')
       ax.quiver(self.xt_gl[::2*fac],M.zt[::2*fac],u[::2*fac,::2*fac].transpose(),w[::2*fac,::2*fac].transpose() )
       
       ax.set_title('Temperature')
       ax.set_xlabel('x [m]')
       ax.set_ylabel('z [m]')
       ax.axis('tight')
       return
  
if __name__ == "__main__": 
     internal_wave1().run(snapint = 0.5 ,runlen = 50.0)

