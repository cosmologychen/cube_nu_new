program nfw_force
  use omp_lib
  use parameters
  use variables
  use,intrinsic :: ISO_C_BINDING
  implicit none
  include 'fftw3.f03'

  integer,parameter :: Interlace=0
  integer,parameter :: nlayer=3
  integer,parameter :: nsample=1 ! undersample particles
  real,parameter :: x_zoom=0.2 ! rescale factor ! grid = x_zoom * kpc/h


  integer i,j,k,l,i1,i2,i3,np,ip,it(3),ix(3),l3(3),idx(3,3),dx(3,3),npsum,ixp,nzero,ilayer,iw,ntotal
  integer,allocatable :: cumsum(:,:,:,:,:,:),rho_local(:,:,:,:,:,:)
  real npr(2,1000),x_center(3),x_sim(3),x_coarse(3),rvec(3),rmag,fmag,rmass,f0mag
  real,allocatable :: rho(:,:,:),v(:,:,:,:),xpos(:,:),xi(:,:)

  x_center=ng/2+0.5 ! center of the NFW halo

  call geometry
  call initialize_pp
  call omp_set_num_threads(ncore)
  if (head) then
    print*,'kernel_initialization on',int(ncore,kind=2),' cores'
    print*,'  nc,ngt,nft =',nc,ngt,nft
  endif

  allocate(rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[nn,nn,*])

  open(11,file='../../visualization/plane/nfw_uniform_particle.bin',access='stream')
  read(11) ntotal
  allocate(xpos(3,ntotal))
  read(11) xpos
  close(11)
  nplocal=ntotal/nsample
  sim%nplocal=nplocal
  print*,'nplocal =',nplocal
  print*,minval(xpos,2),maxval(xpos,2)

  ! convert to zip format
  rhoc=0
  np_image_max=nplocal*8
  np_tile_max=np_image_max
  allocate(xp(3,np_image_max)[nn,nn,*],vp(3,np_image_max)[nn,nn,*],vp_store(3,nplocal,4))
  allocate(xp_new(3,np_tile_max),vp_new(3,np_tile_max))
  xp=0; vp=0; vp_store=0
  do ip=nsample,ntotal,nsample ! loop 1
    x_sim=x_center+xpos(:,ip)*x_zoom
    x_coarse=x_sim/ratio_cs
    it=floor(x_coarse/nt)+1
    ix=floor(x_coarse-nt*(it-1))+1
    rhoc(ix(1),ix(2),ix(3),it(1),it(2),it(3))=rhoc(ix(1),ix(2),ix(3),it(1),it(2),it(3))+1
  enddo

  allocate(cumsum(nt,nt,nt,nnt,nnt,nnt))
  npsum=0
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    do k=1,nt
    do j=1,nt
    do i=1,nt
      cumsum(i,j,k,itx,ity,itz)=npsum
      npsum=npsum+rhoc(i,j,k,itx,ity,itz)
    enddo
    enddo
    enddo
  enddo
  enddo
  enddo
  print*,'npsum',npsum
  print*,'sum(rhoc)',sum(rhoc)

  allocate(rho_local(nt,nt,nt,nnt,nnt,nnt))
  do ip=nsample,ntotal,nsample
    x_sim=x_center+xpos(:,ip)*x_zoom
    x_coarse=x_sim/ratio_cs
    it=floor(x_coarse/nt)+1
    ix=floor(x_coarse-nt*(it-1))+1
    rho_local(ix(1),ix(2),ix(3),it(1),it(2),it(3))=rho_local(ix(1),ix(2),ix(3),it(1),it(2),it(3))+1
    ixp = cumsum(ix(1),ix(2),ix(3),it(1),it(2),it(3))+rho_local(ix(1),ix(2),ix(3),it(1),it(2),it(3))
    xp(:,ixp)=floor( (x_coarse-floor(x_coarse))/x_resolution,kind=8)
  enddo

  deallocate(xpos,cumsum,rho_local)

  print*,'buffer_grid'
  call buffer_grid
  print*,'buffer_x'
  call buffer_x
  print*,'buffer_v'
  call buffer_v

  ! calculate buffered index
  allocate(ip2iw(sum(rhoc)))
  iw=0
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
      do k=1,nt
      do j=1,nt
      do i=1,nt
        np=rhoc(i,j,k,itx,ity,itz)
        nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
        do l=1,np
          ip=nzero+l
          iw=iw+1
          ip2iw(ip)=iw
        enddo
      enddo
      enddo
      enddo
  enddo
  enddo
  enddo
  !print*,ip,iw,minval(ip2iw),maxval(ip2iw); stop
  

  ! initialize Green's functions
  print*,'initialize Green''s functions'
  allocate(Gk1_iso(nc+1,2*nc,2*nc),Gk2(ngt/2+1,ngt,ngt),Gk3(nft/2+1,nft,nft))
  call Green_3D_iso(Gk1_iso,2*nc,nc+1,2*nc,2*nc,apm1c,0.,real(ratio_cs))
  call Green_3D(Gk2,ngt,ngt/2+1,ngt,ngt,apm2,apm1,1.)
  call Green_3D(Gk3,nft,nft/2+1,nft,nft,apm3f,apm2f,1./ratio_sf)

  ! FFT plans
  print*,'create FFT plans'
  call system_clock(t1,t_rate)
  call sfftw_plan_dft_r2c_3d(plan1, 2*nc, 2*nc, 2*nc, rho1_iso,rho1_iso,FFTW_MEASURE)
  call sfftw_plan_dft_c2r_3d(iplan1,2*nc, 2*nc, 2*nc, rho1_iso,rho1_iso,FFTW_MEASURE)
  call sfftw_plan_dft_r2c_3d(plan2,  ngt,  ngt,  ngt, rho2,rho2,FFTW_MEASURE)
  call sfftw_plan_dft_c2r_3d(iplan2, ngt,  ngt,  ngt, rho2,rho2,FFTW_MEASURE)
  call sfftw_plan_dft_r2c_3d(plan3,  nft,  nft,  nft, rho3,rho3,FFTW_MEASURE)
  call sfftw_plan_dft_c2r_3d(iplan3, nft,  nft,  nft, rho3,rho3,FFTW_MEASURE)
  call system_clock(t2,t_rate)
  print*,'    elapsed time =',real(t2-t1)/t_rate,'secs'; print*,''

  call kick

  print*,'force analysis'


  open(11,file='../../visualization/plane/analytical_particle_number.bin',access='stream')
  read(11) npr
  close(11)

  print*, 'interpolation test'
  print*,npr(1,1),npr(2,1)
  print*,npr(1,1000),npr(2,1000)
  print*, interp_npr(1.),interp_npr(450.)

  ! density field
  !allocate(rho(ng,ng,ng),v(3,ng,ng,ng))
  allocate(xi(10,nplocal))
  !rho=0; v=0; 
  iw=0
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    print*,itx,ity,itz
    !do ilayer=0,nlayer-1
      !!$omp paralleldo default(shared) schedule(dynamic)&
      !!$omp& private(k,j,i,np,nzero,l,ip,x_sim,idx,l3,dx,i3,i2,i1)
      !do k=1+ilayer,nt,nlayer
      do k=1,nt
        !print*,k
      do j=1,nt
      do i=1,nt
        np=rhoc(i,j,k,itx,ity,itz)
        nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
        do l=1,np
          ip=nzero+l
          iw=iw+1
          x_sim=ratio_cs*(nt*([itx,ity,itz]-1)+([i,j,k]-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution)
          
          rvec=x_sim-x_center
          rmag=norm2(rvec)
          fmag=norm2(vp(:,ip))
          rmass=interp_npr(rmag/x_zoom)
          f0mag=rmass/rmag**2/nsample

          xi(1,iw)=rmag/x_zoom
          xi(2,iw)=fmag
          xi(3,iw)=f0mag
          
          !idx(:,2)=floor(x_sim)+1
          !idx(:,1)=idx(:,2)-1
          !idx(:,3)=idx(:,2)+1
          !l3=x_sim-floor(x_sim)
          !dx(:,1)=(1-l3)**2/2
          !dx(:,3)=l3**2/2
          !dx(:,2)=1-dx(:,1)-dx(:,3)
          !do i3=1,3
          !do i2=1,3
          !do i1=1,3
          !  rho(idx(1,i1),idx(2,i2),idx(3,i3))=rho(idx(1,i1),idx(2,i2),idx(3,i3))+dx(1,i1)*dx(2,i2)*dx(3,i3)
          !  v(:,idx(1,i1),idx(2,i2),idx(3,i3))=v(:,idx(1,i1),idx(2,i2),idx(3,i3))+vp(:,ip)*dx(1,i1)*dx(2,i2)*dx(3,i3)
          !enddo
          !enddo
          !enddo
        enddo
      enddo
      enddo
      enddo
      !!$omp endparalleldo
    !enddo
  enddo
  enddo
  enddo
  !v(1,:,:,:)=v(1,:,:,:)/(rho+0.000001)
  !v(2,:,:,:)=v(2,:,:,:)/(rho+0.000001)
  !v(3,:,:,:)=v(3,:,:,:)/(rho+0.000001)



  !open(11,file='../../visualization/plane/rho.bin',status='replace',access='stream')
  !write(11) rho
  !close(11)  
  !open(11,file='../../visualization/plane/v.bin',status='replace',access='stream')
  !write(11) v
  !close(11)

  print*,'write',iw,' particles'

  open(11,file='../../visualization/plane/xi.bin',status='replace',access='stream')
  write(11) iw,x_zoom
  write(11) xi
  write(11) vp_store
  close(11)


  contains

  real function interp_npr(r) ! interpolation in log space
    implicit none
    integer ii,i1,i2
    real r,xx,yy,x1,x2,y1,y2
    i1=1
    i2=1000
    do while (i2-i1>1)
      ii=(i1+i2)/2
      if (r>npr(1,ii)) then
        i1=ii
      else
        i2=ii
      endif
    enddo
    x1=log(npr(1,i1))
    y1=log(npr(2,i1))
    x2=log(npr(1,i2))
    y2=log(npr(2,i2))
    xx=log(r)
    yy=y1+(y2-y1)*(xx-x1)/(x2-x1)
    interp_npr=exp(yy)
  endfunction

  subroutine initialize_pp
    implicit none
    print*,'  initialize PP neighbors'
    l=0
    do k=-nrange,-1
    do j=-nrange,nrange
    do i=-nrange,nrange
      l=l+1
      ijk(:,l)=[i,j,k]
    enddo
    enddo
    enddo
    k=0
    do j=-nrange,-1
    do i=-nrange,nrange
      l=l+1
      ijk(:,l)=[i,j,k]
    enddo
    enddo
    j=0
    do i=-nrange,-1
      l=l+1
      ijk(:,l)=[i,j,k]
    enddo

  endsubroutine

end
