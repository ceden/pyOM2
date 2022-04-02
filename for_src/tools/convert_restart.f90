! converts a restart in different unformated files written
! by each PE to files of a different domain decomposition
! needs at least three 3D fields memory
! reads from restart_PE_?????.dta and writes to new_restart_PE_?????.dta

program convert_restart
 implicit none
 integer :: onx =2,err
 integer :: nx,ny,nz,itt
 integer :: is,ie,js,je,n_pes_i,n_pes_j
 integer :: npes,n,k,io=10
 real*8,allocatable :: a(:,:,:,:)
 character*80 :: filename,arg
 logical :: enable_streamfunction = .false.
 logical :: enable_eke            = .false.
 logical :: enable_tke            = .true.
 logical :: enable_idemix         = .true.

 ! read new domain decomposition from command line input
 if (iargc() < 2) then
    print*,' not enough command line input'
    stop
 endif
 call getarg(1,arg); read(arg,*) n_pes_i
 call getarg(2,arg); read(arg,*) n_pes_j
 print*,'converting to ',n_pes_i,' x ',n_pes_j ,' PEs'

 ! read basic information from first input file
 n=0
 write(filename,'(a,i5,a)')  'restart_PE_',n,'.dta'
 call replace_space_zero(filename)
 print*,'reading from file ',filename(1:len_trim(filename))
 open(io,file=filename,form='unformatted',status='old')
 read(io) nx,ny,nz,itt
 print*,' nx,ny,nz,itt = ',nx,ny,nz,itt 
 close(io)

 ! allocate the memory
 print*,'allocating memory'
 allocate(a(1-onx:nx+onx,1-onx:ny+onx,1:nz,1:3),stat=err ); a=0.     

 ! check how many input files are there
 do n=1,9999999
   write(filename,'(a,i5,a)')  'restart_PE_',n,'.dta'
   call replace_space_zero(filename)
   open(io,file=filename,form='unformatted',status='old',err=10)
   close(io)
 enddo
10 continue
 npes = n-1
 print*,' found ',npes+1,' files'


 ! initialize the output files
 call init_output(itt,nx,ny,nz,onx,n_pes_i,n_pes_j)

 ! read from input and append to output files
 print*,'reading/writing rho' 
 call read_chunk(0,npes,a,nx,ny,nz,onx,.False.,.False.)
 call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
 print*,'reading/writing Nsqr' 
 call read_chunk(1,npes,a,nx,ny,nz,onx,.False.,.False.)
 call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
 print*,'reading/writing HD' 
 call read_chunk(2,npes,a,nx,ny,nz,onx,.False.,.False.)
 call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
 print*,'reading/writing int_drhodT' 
 call read_chunk(3,npes,a,nx,ny,nz,onx,.False.,.False.)
 call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
 print*,'reading/writing int_drhodS' 
 call read_chunk(4,npes,a,nx,ny,nz,onx,.False.,.False.)
 call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
 print*,'reading/writing dHd' 
 call read_chunk(5,npes,a,nx,ny,nz,onx,.False.,.False.)
 call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
 print*,'reading/writing temp' 
 call read_chunk(6,npes,a,nx,ny,nz,onx,.False.,.False.)
 call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
 print*,'reading/writing dtemp' 
 call read_chunk(7,npes,a,nx,ny,nz,onx,.False.,.False.)
 call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
 print*,'reading/writing salt' 
 call read_chunk(8,npes,a,nx,ny,nz,onx,.False.,.False.)
 call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
 print*,'reading/writing dsalt' 
 call read_chunk(9,npes,a,nx,ny,nz,onx,.False.,.False.)
 call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
 print*,'reading/writing u' 
 call read_chunk(10,npes,a,nx,ny,nz,onx,.False.,.False.)
 call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
 print*,'reading/writing du' 
 call read_chunk(11,npes,a,nx,ny,nz,onx,.False.,.False.)
 call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
 print*,'reading/writing v' 
 call read_chunk(12,npes,a,nx,ny,nz,onx,.False.,.False.)
 call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
 print*,'reading/writing dv' 
 call read_chunk(13,npes,a,nx,ny,nz,onx,.False.,.False.)
 call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
 print*,'reading/writing w' 
 call read_chunk(14,npes,a,nx,ny,nz,onx,.False.,.False.)
 call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
 print*,'reading/writing psi' 
 call read_chunk(15,npes,a,nx,ny,nz,onx,.True.,.False.)
 call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.True.,.False.)

 if (enable_streamfunction) then
         print *,' not yet implemented'
         stop
 endif


 n=16
 if (enable_eke) then
   print*,'reading/writing eke'
   call read_chunk(n,npes,a,nx,ny,nz,onx,.False.,.False.)
   call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
   print*,'reading/writing deke'
   call read_chunk(n+1,npes,a,nx,ny,nz,onx,.False.,.False.)
   call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
   n=n+2
 endif  
 
 if (enable_tke) then
   print*,'reading/writing tke'
   call read_chunk(n,npes,a,nx,ny,nz,onx,.False.,.True.)
   call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.True.)
   print*,'reading/writing dtke'
   call read_chunk(n+1,npes,a,nx,ny,nz,onx,.False.,.False.)
   call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
   n=n+2
 endif

 if (enable_idemix) then
   print*,'reading/writing E_iw'
   call read_chunk(n,npes,a,nx,ny,nz,onx,.False.,.False.)
   call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
   print*,'reading/writing dE_iw'
   call read_chunk(n+1,npes,a,nx,ny,nz,onx,.False.,.False.)
   call write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,.False.,.False.)
   n=n+2
 endif


end program

       

subroutine init_output(itt,nx,ny,nz,onx,n_pes_i,n_pes_j)
 implicit none
 integer :: itt,nx,ny,nz,onx,n_pes_i,n_pes_j
 integer :: my_pe,i_blk,j_blk,is_pe,ie_pe,js_pe,je_pe
 integer :: is,ie,js,je,my_blk_i,my_blk_j
 integer :: io=10
 character*80 :: filename

 do my_pe = 0,n_pes_i*n_pes_j-1
   i_blk = (nx-1)/n_pes_i + 1    ! i-extent of each block
   j_blk = (ny-1)/n_pes_j + 1    ! j-extent of each block
   my_blk_i = mod(my_pe,n_pes_i)+1! number of PE in i-dir.
   my_blk_j = (my_pe)/n_pes_i + 1 ! number of PE in j-dir.
   is_pe = (my_blk_i-1)*i_blk + 1 ! start index in i-dir of this PE
   ie_pe = min(my_blk_i*i_blk,nx)
   js_pe = (my_blk_j-1)*j_blk + 1
   je_pe = min(my_blk_j*j_blk,ny)
   is=is_pe-onx; ie=ie_pe+onx; js=js_pe-onx; je=je_pe+onx

   write(filename,'(a,i5,a)')  'new_restart_PE_',my_pe,'.dta'
   call replace_space_zero(filename)
   open(io,file=filename,form='unformatted',status='new')
   write(io) nx,ny,nz,itt
   write(io) is,ie,js,je
   close(io)
 enddo
end subroutine


subroutine read_chunk(k,npes,a,nx,ny,nz,onx,is2D,is3arr)
 implicit none
 logical :: is2D,is3arr
 integer :: kk,k,npes,is,ie,js,je,i
 integer :: nx,ny,nz,onx,n,io=10
 real*8 :: a(1-onx:nx+onx,1-onx:ny+onx,1:nz,1:3)
 real*8,allocatable :: b1(:,:,:),b2(:,:,:),b3(:,:,:)
 character*80 :: filename
 
 do n=0,npes
   write(filename,'(a,i5,a)')  'restart_PE_',n,'.dta'
   call replace_space_zero(filename)
   !print*,'reading from file ',filename(1:len_trim(filename))
   open(io,file=filename,form='unformatted',status='old')
   read(io) 
   read(io) is,ie,js,je
   if (k>0) then 
     do kk=1,k; read (io); enddo
   endif  
   allocate(b1(is:ie,js:je,1:nz),b2(is:ie,js:je,1:nz),b3(is:ie,js:je,1:nz))
   b1=0;b2=0;b3=0
   if (is2D) then 
    if (is3arr) then 
      read (io) b1(is:ie,js:je,1),b2(is:ie,js:je,1),b3(is:ie,js:je,1) 
    else
      read (io) b1(is:ie,js:je,1),b2(is:ie,js:je,1) 
    endif        
   else        
    if (is3arr) then 
     read (io) b1(is:ie,js:je,1:nz),b2(is:ie,js:je,1:nz),b3(is:ie,js:je,1:nz) 
    else 
     read (io) b1(is:ie,js:je,1:nz),b2(is:ie,js:je,1:nz) 
    endif 
   endif
   ! transfer only inner part of subdomain
   a(is+onx:ie-onx,js+onx:je-onx,:,1)=b1(is+onx:ie-onx,js+onx:je-onx,:)
   a(is+onx:ie-onx,js+onx:je-onx,:,2)=b2(is+onx:ie-onx,js+onx:je-onx,:)
   a(is+onx:ie-onx,js+onx:je-onx,:,3)=b3(is+onx:ie-onx,js+onx:je-onx,:)
   deallocate(b1,b2,b3)
   close(io)
 enddo
 ! cyclic boundary conditions
 do i=1,onx
   a(nx+i,:,:,:) = a(i     ,:,:,:)
   a(1-i ,:,:,:) = a(nx-i+1,:,:,:) 
 enddo
end subroutine 



subroutine write_chunk(a,nx,ny,nz,onx,n_pes_i,n_pes_j,is2D,is3arr)
 implicit none
 logical :: is2D,is3arr
 integer :: nx,ny,nz,onx,n_pes_i,n_pes_j
 integer :: my_pe,i_blk,j_blk,is_pe,ie_pe,js_pe,je_pe
 integer :: is,ie,js,je,my_blk_i,my_blk_j
 integer :: io=10
 real*8 :: a(1-onx:nx+onx,1-onx:ny+onx,1:nz,1:3)
 character*80 :: filename

 do my_pe = 0,n_pes_i*n_pes_j-1
   i_blk = (nx-1)/n_pes_i + 1    ! i-extent of each block
   j_blk = (ny-1)/n_pes_j + 1    ! j-extent of each block
   my_blk_i = mod(my_pe,n_pes_i)+1! number of PE in i-dir.
   my_blk_j = (my_pe)/n_pes_i + 1 ! number of PE in j-dir.
   is_pe = (my_blk_i-1)*i_blk + 1 ! start index in i-dir of this PE
   ie_pe = min(my_blk_i*i_blk,nx)
   js_pe = (my_blk_j-1)*j_blk + 1
   je_pe = min(my_blk_j*j_blk,ny)
   is=is_pe-onx; ie=ie_pe+onx; js=js_pe-onx; je=je_pe+onx

   write(filename,'(a,i5,a)')  'new_restart_PE_',my_pe,'.dta'
   call replace_space_zero(filename)
   !print*,'appending to file ',filename(1:len_trim(filename))
   open(io,file=filename,form='unformatted',status='old',access='append')
   if (is2D) then
    if (is3arr) then
      write (io) a(is:ie,js:je,1,1),a(is:ie,js:je,1,2),a(is:ie,js:je,1,3) 
    else  
      write (io) a(is:ie,js:je,1,1),a(is:ie,js:je,1,2) 
    endif  
   else        
    if (is3arr) then
     write (io) a(is:ie,js:je,1:nz,1),a(is:ie,js:je,1:nz,2),a(is:ie,js:je,1:nz,3)  
    else 
     write (io) a(is:ie,js:je,1:nz,1),a(is:ie,js:je,1:nz,2) 
    endif
   endif
   close(io)
 enddo
end subroutine


 subroutine replace_space_zero(name)
      implicit none
      character (len=*) :: name
      integer  :: i
      do i=1,len_trim(name)
          if (name(i:i)==' ')name(i:i)='0'
      enddo
 end subroutine replace_space_zero


