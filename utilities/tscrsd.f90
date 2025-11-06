!! Include RSD effect of particles, and CIC interpolate into grid to get the redshift-space density field
!#define RSD_ELUCID
program tscrsd
  use parameters
  implicit none
  save
  ! nc: coarse grid per node per dim
  ! nf: fine grid per node per dim
  real,parameter :: density_buffer=1.2
  real,parameter :: cen(3)=[30,370,370]
  integer,parameter :: p=2 ! TSC order
  integer,parameter :: ndim=3
  integer(8) i,j,k,l,i_dim,iq(3),nplocal,nplocal_nu,itx,ity,itz
  integer k0,j0,i0,idx(ndim,p+1)
  integer(8) nlast,ip,np,idx1(3),idx2(3)
  integer cur_checkpoint

  ! real(4) rho_grid(0:ng+1,0:ng+1,0:ng+1)[*]
  real,allocatable :: rho_grid(:,:,:)[:,:,:],rho_c(:,:,:),rho_nu(:,:,:),rho_two(:,:,:)
  real(4) sigma_vi,zshift,los(3),sx(3)  !rho_c(ng,ng,ng),
  real(4) mass_p,pos1(3),dx1(3),dx2(3)
  real(8) rho8[*]
  integer(izipx),allocatable :: xp(:,:)
  integer(izipv),allocatable :: vp(:,:)
  ! integer(4) rhoc(nt,nt,nt,nnt,nnt,nnt)
  real(4) vc(3,nt,nt,nt,nnt,nnt,nnt)
  real dx(ndim,p+1),l3(ndim)
  integer(4),allocatable :: rhoc(:,:,:,:,:,:)
  character(20) str_z,str_i
  
  call geometry
  if (head) then
    print*, 'tscrsd on resolution:'
    print*, 'ng=',ng
    print*, 'ng*nn=',ng*nn
    print*, 'cen=',cen
  endif
  sync all

  allocate(rho_grid(0:ng+1,0:ng+1,0:ng+1)[nn,nn,*])
  allocate(rho_c(ng,ng,ng),rho_nu(ng,ng,ng))
  allocate(rhoc(nt,nt,nt,nnt,nnt,nnt))
  allocate(rho_two(ng,ng,ng))

  if (head) then
    print*, 'checkpoint at:'
    open(16,file='../main/z_checkpoint.txt',status='old')
    do i=1,nmax_redshift
      read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
      print*, z_checkpoint(i)
    enddo
    71 n_checkpoint=i-1
    close(16)
    print*,''
  endif

  sync all
  n_checkpoint=n_checkpoint[1]
  z_checkpoint(:)=z_checkpoint(:)[1]
  sync all

  do cur_checkpoint= n_checkpoint,n_checkpoint
    sim%cur_checkpoint=cur_checkpoint
    print*,cur_checkpoint,sim%cur_checkpoint
    if (head) print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))

    !call particle_initialization
    print*,output_name('info')
    open(11,file=output_name('info'),status='old',action='read',access='stream')
    read(11) sim
    close(11)
    sigma_vi=sim%sigma_vi
    !call print_header(sim); stop
    if (sim%izipx/=izipx .or. sim%izipv/=izipv) then
      print*, 'zip format incompatable'
      close(11)
      stop
    endif


    
    !mass_p=sim%mass_p
    mass_p=1.
    nplocal=sim%nplocal
    nplocal_nu=sim%nplocal_nu
    if (head) then
      print*, 'mass_p =',mass_p
      print*, 'nplocal =',nplocal
      print*, 'nplocal_nu =',nplocal_nu
    endif
    !cdm
    allocate(xp(3,nplocal),vp(3,nplocal))
    open(11,file=output_name('xp'),status='old',action='read',access='stream')
    read(11) xp
    close(11)
    open(11,file=output_name('vp'),status='old',action='read',access='stream')
    read(11) vp
    close(11)
    open(11,file=output_name('np'),status='old',action='read',access='stream')
    read(11) rhoc
    close(11)
    open(11,file=output_name('vc'),status='old',action='read',access='stream')
    read(11) vc
    close(11)

    rho_grid=0
    nlast=0
    do itz=1,nnt
    do ity=1,nnt
    do itx=1,nnt
      do k=1,nt
      do j=1,nt
      do i=1,nt
        np=rhoc(i,j,k,itx,ity,itz)
        do l=1,np
          ip=nlast+l
          pos1=nt*((/itx,ity,itz/)-1)+ ((/i,j,k/)-1) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution ! in coarse grid
          pos1=pos1*real(ng)/real(nc)
#ifdef RSD_ELUCID
          sx=vc(:,i,j,k,itx,ity,itz)
          sx=sx+tan((pi*real(vp(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sigma_vi*vrel_boost))
          los=pos1-(cen/500.)*ng ! line of sight vector
          los=los/norm2(los) ! line of sight unit vector
          sx=sum(sx*los)*los ! project velocity to line of sight
          sx=sx*sim%vsim2phys/sim%a/100 ! convert to km/h and multiply 1/aH, in Mpc
          sx=sx/(box/nf_global) ! convert to find grid
          pos1=pos1+sx
          pos1=modulo(pos1,real(ng))
#else
          zshift=vc(zdim,i,j,k,itx,ity,itz) ! coarse grid velocity field
          zshift=zshift+tan((pi*real(vp(zdim,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sigma_vi*vrel_boost))
          zshift=zshift*sim%vsim2phys/sim%a/100 ! convert to km/h and multiply 1/aH, in Mpc
          zshift=zshift*(nf_global/box) ! convert to find grid
          pos1(zdim)=pos1(zdim)+zshift/h0**2 ! add shift field
          pos1(zdim)=modulo(pos1(zdim),real(ng)) ! peridoc over box
#endif 
          idx(:,2)=floor(pos1)+1
          idx(:,1)=idx(:,2)-1
          idx(:,3)=idx(:,2)+1
  
          l3=pos1-floor(pos1)
          dx(:,3)=l3**2/2
          dx(:,1)=(1-l3)**2/2
          dx(:,2)=1-dx(:,1)-dx(:,3)
  
          do k0=1,p+1
          do j0=1,p+1
          do i0=1,p+1
            rho_grid(idx(1,i0),idx(2,j0),idx(3,k0))=rho_grid(idx(1,i0),idx(2,j0),idx(3,k0))+dx(1,i0)*dx(2,j0)*dx(3,k0)*mass_p
          enddo
          enddo
          enddo
  
        enddo
        nlast=nlast+np
      enddo
      enddo
      enddo
    enddo
    enddo
    enddo
    sync all
    deallocate(xp)

    if (head) print*, 'Start sync from buffer regions'
    sync all
    rho_grid(1,:,:)=rho_grid(1,:,:)+rho_grid(ng+1,:,:)[inx,icy,icz]
    rho_grid(ng,:,:)=rho_grid(ng,:,:)+rho_grid(0,:,:)[ipx,icy,icz]; sync all
    rho_grid(:,1,:)=rho_grid(:,1,:)+rho_grid(:,ng+1,:)[icx,iny,icz]
    rho_grid(:,ng,:)=rho_grid(:,ng,:)+rho_grid(:,0,:)[icx,ipy,icz]; sync all
    rho_grid(:,:,1)=rho_grid(:,:,1)+rho_grid(:,:,ng+1)[icx,icy,inz]
    rho_grid(:,:,ng)=rho_grid(:,:,ng)+rho_grid(:,:,0)[icx,icy,ipz]; sync all
    !rho_c=rho_grid(1:ng,1:ng,1:ng)
    do i=1,ng
    do j=1,ng
    do k=1,ng
      rho_c(k,j,i)=rho_grid(k,j,i)
    enddo
    enddo
    enddo

    print*,rho_grid(1:2,1:2,1)
    print*,rho_grid(1:2,1:2,2)
    print*,rho_c(1:2,1:2,1)
    print*,rho_c(1:2,1:2,2)
    print*, 'check: min,max,sum of rho_grid = '
    print*, minval(rho_c),maxval(rho_c),sum(rho_c*1d0)

    rho8=sum(rho_c*1d0); sync all
    ! co_sum
    if (head) then
      do i=2,nn**3
        rho8=rho8+rho8[i]
      enddo
      print*,'rho_global',rho8,ng_global
    endif; sync all
    rho8=rho8[1]; sync all
    ! convert to density contrast
    do i=1,ng
      rho_c(:,:,i)=rho_c(:,:,i)/(rho8/ng_global/ng_global/ng_global)-1
    enddo
    !rho_c=rho_c/(rho8/ng_global/ng_global/ng_global)-1
    ! check normalization
    print*,'min',minval(rho_c),'max',maxval(rho_c),'mean',sum(rho_c*1d0)/ng/ng/ng; sync all

    if (head) print*,'Write delta_tscrsd into',output_name('delta_tscrsd')
    open(11,file=output_name('delta_s'),status='replace',access='stream')
    write(11) rho_c
    close(11); sync all
    
    !open(11,file=output_name('delta_tscrsd_proj'),status='replace',access='stream')
    !write(11) sum(rho_c(:,:,:50),dim=3)/50
    !close(11); sync all

    ! power spectrum
    write(str_i,'(i6)') image
    write(str_z,'(f7.3)') z_checkpoint(cur_checkpoint)


    sync all

  enddo
  sync all
  if (head) print*,'tscrsd done'
  sync all
end
