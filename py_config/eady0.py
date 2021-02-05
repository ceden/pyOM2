import sys; sys.path.append('../py_src')
from pyOM_ave import pyOM_ave as pyOM
from numpy import *

EDDY_LENGTH_SCALES = 2 # number of eddies in domain


Ro  = 0.05
Ek  = 0.001
CFL = 0.02    # CFL number
delta = 0.01 # aspect ratio

if 1:
  Ro  = sqrt( 1./20)
  CFL = 0.05    # CFL number  
  delta = 0.02 # aspect ratio
  Ek  = 0.010

H0 = 200.0
beta= 0e-11  # beta 
f0  = 1e-4   # coriolis freq.
N0=f0/delta  # stability freq.
U0 = Ro*N0*H0



class eady0(pyOM):
   """ Eady (1941) solution
   """
   def set_parameter(self):
     """set main parameter
     """
     M=self.fortran.main_module   
     
     #(M.nx,M.ny,M.nz)    = (120,120,40)
     (M.nx,M.ny,M.nz)    = (64,64,20)
     
     M.congr_epsilon = 1e-12
     M.congr_max_iterations = 5000
     #M.enable_streamfunction = 1
     M.enable_free_surface = 1
     
     M.enable_hydrostatic      = 1
     M.enable_cyclic_x         = 1
     
     M.enable_superbee_advection  = 1
     M.enable_explicit_vert_friction = 1
     M.enable_biharmonic_friction = 1 

     #M.enable_ray_friction = 1
     #M.r_ray = 1/(15*86400.)

     M.enable_conserve_energy = 0
     M.coord_degree           = 0
     M.eq_of_state_type       = 1 
     M.enable_tempsalt_sources = 1 
     M.enable_momentum_sources = 1 
     return



   def set_grid(self):
     import lin_stab
     M=self.fortran.main_module   
     # print some numbers
     Lr    = N0*H0/f0          # Rossby radius
     Ro    = U0/(f0*Lr)       # Rossby number
     Ri    = N0**2*H0**2/U0**2 # Richardson number
     M0    = sqrt(f0*U0/H0)    #    
     if M.my_pe == 0:
       print
       print  ' L  = %f km'%(Lr/1e3)
       print  ' Ro = %f '%Ro
       print  ' Ri = %f '%Ri
       print  ' delta = %f '%delta
       print  ' ell = %f '%(Lr/6400e3)

     # solve linear stability problem first
     ksx=linspace(0,3.2,50);    kx=ksx/Lr
     ky = array( [0./Lr] )
     #ksy=linspace(-3.2,3.2,30); ky=ksy/Lr
     M.dzt[:]    = H0/M.nz
     zw=arange(M.nz)*M.dzt[0]+ M.dzt[0]
     zt=zw-M.dzt[0]/2.0
     zt = zt - zw[-1]
     zw = zw - zw[-1]
     U=U0/2+U0*zt/H0
     V=U*0
     B=N0**2*zt
     om_max,om,kmax,lmax,u,v,w,b,p=lin_stab.pe(U,V,B,M.dzt[0],kx,ky,0.,beta,0.,f0,0.,0)
     if M.my_pe == 0:
       print ' Max. growth rate %f 1/days ' % (-imag(om)*86400)
       print ' k_max = %f Lr , l_max = %f Lr' % (kmax*Lr,lmax*Lr)
       #print ' theoretical om_max = ',sqrt(5./54)/sqrt(1+Ri)*f0 *86400
       #L0 = 2*pi*sqrt(2./5.)*sqrt(1+Ri)*M0**2/f0**2*H0
       #print ' theoretical k_max = ',2*pi/L0 * Lr

     self.lin_stab_om = om
     self.lin_stab_kmax = kmax
     self.lin_stab_b, self.lin_stab_u, self.lin_stab_v, self.lin_stab_w, self.lin_stab_p = b,u,v,w,p
     L = EDDY_LENGTH_SCALES*2*pi/kmax   
     M.dxt[:]  =   L/M.nx 
     M.dyt[:]  =   L/M.ny 
     M.dt_mom    = CFL/U0*M.dxt[0]   # CFL=U*dt/dx
     M.dt_tracer = CFL/U0*M.dxt[0]
     if M.my_pe == 0:
       print " dx=%f km, dt= %f s "%(M.dxt[0]/1e3,M.dt_mom)
       print " CFL  = ",U0*M.dt_mom/M.dxt[0]
       print " CFL  = ",real(om)/kmax*M.dt_mom/M.dxt[0]
     M.congr_epsilon = 1e-12 *(M.dxt[0]/20e3)**2

     #M.a_h      = Ek*f0*M.dxt[0]**2 
     M.a_hbi    = Ek*f0*M.dxt[0]**4 
     M.kappam_0 = Ek*f0*M.dzt[0]**2
     
     #M.k_h = Ek*f0*M.dxt[0]**2 
     #M.k_v = Ek*f0*M.dzt[0]**2 
     if M.my_pe == 0:
       print " A_h = %f m^2/s  Ek = %f"%(M.a_h,Ek)
       print " A_v = %f m^2/s  Ek = %f"%(M.kappam_0,Ek)
     return



   
   def set_topography(self):
       M=self.fortran.main_module   
       M.kbot[:]=0
       M.kbot[:,2:-2]=1
       return

   
   def set_coriolis(self):
     M=self.fortran.main_module   
     for j in range( M.yt.shape[0] ): M.coriolis_t[:,j] = f0+beta*M.yt[j]
     return

   def set_forcing(self):
     M=self.fortran.main_module   
     t_rest = abs(imag(self.lin_stab_om))/5.
     u_rest = abs(imag(self.lin_stab_om))/5.

     if (M.n_pes_i >1): 
         print 'Error: numbers of PES in i direction >1'
         stop
     if M.enable_tempsalt_sources: 
         mm = mean( M.temp[2:-2,:,:,M.tau-1] , axis = 0 )
         #mm = zeros( (M.temp.shape[1], M.temp.shape[2] )  )
         #a = zeros( (M.nx,M.ny), 'd', order = 'F')
         #for k in range( M.temp.shape[2]):
         #  a[M.is_pe-1:M.ie_pe,M.js_pe-1:M.je_pe] = M.temp[2:-2,2:-2,k,M.tau-1]
         #  self.fortran.pe0_recv_2d(a)
         #  m_ = mean(a, axis = 0 )
         #  self.fortran.pe0_bcast(m_,M.ny)
         #  mm[2:-2,k] = m_[M.js_pe-1:M.je_pe]
         #
         for i in range(M.temp.shape[0]): 
              M.temp_source[i,:,:]=t_rest*(self.t0[i,:,:]-mm)*M.maskt[i,:,:]
     if M.enable_momentum_sources: 
         mm = mean( M.u[2:-2,:,:,M.tau-1] , axis = 0 )
         #mm = zeros( (M.u.shape[1], M.u.shape[2] )  )
         #a = zeros( (M.nx,M.ny), 'd', order = 'F')
         #for k in range( M.temp.shape[2]):
         #  a[M.is_pe-1:M.ie_pe,M.js_pe-1:M.je_pe] = M.u[2:-2,2:-2,k,M.tau-1]
         #  self.fortran.pe0_recv_2d(a)
         #  m_ = mean(a, axis = 0 )
         #  self.fortran.pe0_bcast(m_,M.ny)
         #  mm[2:-2,k] = m_[M.js_pe-1:M.je_pe]
         for i in range(M.u.shape[0]): 
              M.u_source[i,:,:]=u_rest*(self.u0[i,:,:]-mm)*M.masku[i,:,:]

     # to diagnose pressure
     for k in range(M.u.shape[2]): 
          self.diag_press[:,:,k,M.tau-1] = M.p_hydro[:,:,k] + M.psi[:,:,M.tau-1]
     return

   def set_initial_conditions(self):
     """ setup all initial conditions
     """
     M=self.fortran.main_module   
     kx=1.6*f0/(N0*H0)
     ky=pi/((M.ny-2)*M.dxt[0])
     d=f0/N0/(kx**2+ky**2)**0.5

     fxa=(exp(H0/d)+exp(-H0/d))/(exp(H0/d)-exp(-H0/d))
     c1= (1+0.25*(H0/d)**2-H0/d*fxa )*complex(1,0)
     c1=(sqrt(c1)*d/H0+0.5)*U0
     A=(U0-c1)/U0*H0/d
     
     alpha = self.fortran.linear_eq_of_state.linear_eq_of_state_drhodt()
     grav = 9.81; rho0 = 1024.0
        
     self.u0  = zeros( (M.i_blk+2*M.onx,M.j_blk+2*M.onx,M.nz) , 'd', order='F')
     self.t0  = zeros( (M.i_blk+2*M.onx,M.j_blk+2*M.onx,M.nz) , 'd', order='F')
     self.diag_press  = zeros( (M.i_blk+2*M.onx,M.j_blk+2*M.onx,M.nz,3) , 'd', order='F')
     
     # zonal velocity 
     for k in range(M.nz):
       self.u0[:,:,k]= (U0/2+U0*M.zt[k]/(M.nz*M.dzt[0]))*M.masku[:,:,k] 
       M.u[:,:,k,M.tau-1]= self.u0[:,:,k]
     M.u[...,M.taum1-1] = M.u[...,M.tau-1]
     

     # rho = alpha T ,  N^2 = b_z = - g/rho0 rho_z = - g/rho0 alpha T_z,  T = - N^2 z rho0/(g alpha)
     for k in range(M.nz):
        self.t0[:,:,k]=-N0**2*M.zt[k]/grav/alpha*rho0*M.maskt[:,:,k]
     
     # fu = -p_y, p_z = -g rho,  f u_z = -g rho_y,  rho_y = - f u_z/g = alpha T_y
     
     a = zeros( (M.nx,M.ny), 'd', order = 'F')
     a[0,M.js_pe-1:M.je_pe] = M.dyt[2:-2]
     self.fortran.pe0_recv_2d(a)
     dyt_gl = a[0,:].copy()

     a = zeros( (M.nx,M.ny), 'd', order = 'F')
     for k in range(M.nz):
        for j in range(M.ny-1):
            a[:,j+1]=a[:,j]+dyt_gl[j]*U0/H0*f0/grav/alpha*rho0
            self.fortran.pe0_bcast(a[:,j+1],M.nx)
        self.t0[2:-2,2:-2,k] += a[M.is_pe-1:M.ie_pe,M.js_pe-1:M.je_pe]

         
     #for k in range(M.nz):
     # for j in range(1,M.ny):
     #   jj = self.jf2py(j)
     #   uz = U0/H0#M.ht[:,jj]
     #   self.t0[:,jj+1,k]=(self.t0[:,jj,k]+M.dyt[jj]*uz*f0/grav/alpha*rho0)*M.maskt[:,jj,k]

      
     # perturbation buoyancy
     for k in range(M.nz):
        #for j in range(1,M.ny+1):
        #jj = self.jf2py(j)
        #phiz=A/d*sinh(M.zt[k]/d)+cosh(M.zt[k]/d)/d
        #M.temp[:,jj,k,M.tau-1] =self.t0[:,jj,k] + \
        #      0.1*sin(kx*M.xt)*sin(ky*M.yt[jj])*abs(phiz)*M.maskt[:,jj,k]*rho0/grav/alpha
        phiz=A/d*sinh(M.zt[k]/d)+cosh(M.zt[k]/d)/d
        for j in range(M.temp.shape[1]):
          M.temp[:,j,k,M.tau-1] =self.t0[:,j,k] + \
              0.1*sin(kx*M.xt)*sin(ky*M.yt[j])*abs(phiz)*M.maskt[:,j,k]*rho0/grav/alpha
     M.temp[...,M.taum1-1] = M.temp[...,M.tau-1]
     return

   
   def set_diagnostics(self):
     """ register variables to be averaged here.
     """
     M=self.fortran.main_module   
     self.register_average(name='u',   long_name='Zonal velocity',      units = 'm/s' ,   grid = 'UTT', var = M.u)
     self.register_average(name='v',   long_name='Meridional velocity', units = 'm/s' ,   grid = 'TUT', var = M.v)
     self.register_average(name='w',   long_name='Vertical velocity',   units = 'm/s' ,   grid = 'TTU', var = M.w)
     self.register_average(name='temp',long_name='Temperature',         units = 'deg C' , grid = 'TTT', var = M.temp)
     self.register_average(name='press',long_name='Pressure',           units = 'm^2/s^2' , grid = 'TTT', var = self.diag_press)
     return

if __name__ == "__main__": 
    model = eady0()
    model.run( snapint = 86400.0, runlen = 365*86400.)




