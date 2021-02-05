ny=201
nx=201
n_pes_j=15
n_pes_i=12
j_blk = (ny-1)/n_pes_j + 1  
i_blk = (nx-1)/n_pes_i + 1 
for my_pe in range(n_pes_j*n_pes_i):

 my_blk_i = my_pe % n_pes_i +1
 is_pe = (my_blk_i-1)*i_blk + 1 
 ie_pe = min(my_blk_i*i_blk,nx)
 my_blk_j = (my_pe)/n_pes_i + 1 
 js_pe = (my_blk_j-1)*j_blk + 1
 je_pe = min(my_blk_j*j_blk,ny)
 print 'pe=',my_pe,' j=',js_pe,':',je_pe,' i=',is_pe,':',ie_pe
 if ie_pe<= is_pe:
    print "not possible"
    raise NameError
 if je_pe<= js_pe:
    print "not possible"
    raise NameError

print n_pes_i*n_pes_j/36.

