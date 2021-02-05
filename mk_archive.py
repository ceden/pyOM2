
import os,glob,tarfile

def always_true(file):
   if file == 'bin': return False
   return True

archive_name = 'tmp.tar.gz'
print 'creating archive ',archive_name 
tar = tarfile.open(archive_name,'w:gz')
tar.add('bin',exclude=always_true)
tar.add('doc/pyOM2_2.pdf')
for file in glob.glob('for_config/*.f90'): tar.add(file)
tar.add( 'for_config/Makefile' )

tar.add( 'for_src/Makefile' )
for dir in ['density','diagnostics','eke','etc','external','idemix','isoneutral','qg_filter',
            'rossmix','main','non_hydrostatic','parallel','tke','obc','advection','tracer']:
  tar.add( 'for_src/'+dir+'/Makefile' )
  for file in glob.glob('for_src/'+dir+'/*.f90'): tar.add(file)

for file in glob.glob('py_src/*.py'):   tar.add(file)
tar.add( 'py_src/Makefile' )
for file in glob.glob('py_config/*.py'):  tar.add(file)
for file in glob.glob('site_specific.mk_*'):  tar.add(file)
for file in glob.glob('machenhauer/*.py'):  tar.add(file)
for file in glob.glob('machenhauer/*.f90'):  tar.add(file)
for file in glob.glob('machenhauer/Makefile'):  tar.add(file)
    
for dir in ['flame_high','flame_low','flame_med','global_1deg','global_2deg','global_4deg',
            'global_2deg_45level','global_4deg_45level','pop']: 
   tar.add( 'setups/'+dir+'/Makefile' )
   for file in glob.glob('setups/'+dir+'/setup1.f90' ): tar.add(file)
   file='setups/'+dir+'/setup1.py' 
   if glob.glob(file): tar.add(file)

tar.close()

if 1:
  tar = tarfile.open(archive_name,'r:gz')
  tar.list()
  tar.close()

