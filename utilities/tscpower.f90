!#define write_xv
!define RSD
program tscpower
  use parameters
  use pencil_fft
  use powerspectrum
  implicit none
  save

  real,parameter :: density_buffer=1.2
  integer k0,j0,i0,idx(ndim,p+1)
  integer(8) i,j,k,l,i_dim,iq(3),nplocal,npglobal,itx,ity,itz,nlast,ip,np,idx1(3),idx2(3)
  integer cur_checkpoint, num

  real,allocatable :: rho_grid(:,:,:)[:,:,:],rho_c(:,:,:),rho_etc(:,:,:)
  real mass_p,pos1(3),dx(ndim,p+1),l3(ndim),xi(10,nbin)[*]
  real(8) rho8[*]

  integer(izipx),allocatable :: xp(:,:)
  real(4),allocatable :: xv(:,:)
  integer(4),allocatable :: rhoc(:,:,:,:,:,:)
  character(20) str_z,str_i
  call geometry
  if (head) then
    print*, 'tscpower on resolution:'
    print*, 'ng=',ng
    print*, 'ng*nn=',ng*nn
  endif
  sync all

  allocate(rho_grid(0:ng+1,0:ng+1,0:ng+1)[nn,nn,*])
  allocate(rho_c(ng,ng,ng),rho_etc(ng,ng,ng))
  allocate(rhoc(nt,nt,nt,nnt,nnt,nnt))

  if (head) then
    print*, 'checkpoint at:'
    open(16,file='../z_checkpoint.txt',status='old')
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

  call create_penfft_plan

  do cur_checkpoint= n_checkpoint,1,-1
    sim%cur_checkpoint=cur_checkpoint
    if (head) print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))
    print*,output_name('info')
    open(11,file=output_name('info'),status='old',action='read',access='stream')
    read(11) sim
    close(11)
    if (sim%izipx/=izipx .or. sim%izipv/=izipv) then
      print*, 'zip format incompatable'
      close(11)
      stop
    endif
    mass_p=1.
    nplocal=sim%nplocal
    npglobal=sim%npglobal
    if (head) then
      print*, 'mass_p =',mass_p
      print*, 'nplocal =',nplocal
    endif
    !cdm
    allocate(xp(3,nplocal))
#ifdef write_xv
    allocate(xv(3,nplocal))
#endif
    open(11,file=output_name('xp'),status='old',action='read',access='stream')
    read(11) xp
    close(11)
    open(11,file=output_name('np'),status='old',action='read',access='stream')
    read(11) rhoc
    close(11)

#ifdef write_xv
    open(12,file=output_name('xv'),status='replace',access='stream')
#endif

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

#ifdef write_xv
          xv(:3,ip)=pos1*real(ng)/real(nc)
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

    print*,'a'
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
    print*,'sample'
    print*,rho_c(1:2,1:2,1)
    print*,rho_c(1:2,1:2,2)



    if (head) print*,'Write delta_c into',output_name('delta_c_tsc')
    open(11,file=output_name('delta_c_tsc'),status='replace',access='stream')
    write(11) rho_c
    close(11); sync all
    open(11,file=output_name('delta_c_proj_tsc'),status='replace',access='stream')
    write(11) sum(rho_c(:,:,:),dim=3)/ng
    close(11); sync all

    ! power spectrum
    write(str_i,'(i6)') image
    write(str_z,'(f7.3)') z_checkpoint(cur_checkpoint)

    open(15,file=output_dir()//'delta_L'//output_suffix(),status='old',access='stream')
    read(15) rho_etc
    close(15)
    ! call cross_power(xi,rho_c,rho_etc)
    call auto_power(xi,rho_c,npglobal,3)
    !call density_to_potential(rho_c)
    sync all
    if (head) then
      open(15,file=output_name('tscpower'),status='replace',access='stream')
      write(15) xi
      close(15)
    endif
    sync all

#ifdef write_xv
    write(12) xv
    close(12)
    deallocate(xv)
#endif

  enddo
  deallocate(rho_c,rho_etc,rho_grid,rhoc)
  call destroy_penfft_plan
  sync all
  if (head) print*,'tscpower done'
  sync all
end