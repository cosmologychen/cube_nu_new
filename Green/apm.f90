
subroutine apm(density)
   use omp_lib
   use variables
   implicit none
   save
   real,allocatable :: rho_ex(:,:,:),phi(:,:,:)[:,:,:],force(:,:,:,:)
   integer,parameter :: nlayer=3 ! thread save for TSC intepolation
   integer i,j,k,l,np,idx(ndim,p+1),ilayer,i1,i2,i3,ix,iy,iz,n1(ndim),n2(ndim)
   integer(8) ip,nzero
   real :: density(pm%nstart(1):pm%nend(1),pm%nstart(2):pm%nend(2),pm%nstart(3):pm%nend(3))
   real dx(ndim,p+1),xpos(ndim),l3(ndim),dv(ndim),vreal(ndim)

   !print*,'apm',pm%pm_layer
   !print*,'  allocate rho_ex:',1-pm%nex,pm%nphy+pm%nex
   allocate(rho_ex(1-pm%nex:pm%nphy+pm%nex,1-pm%nex:pm%nphy+pm%nex,1-pm%nex:pm%nphy+pm%nex))
   density=0
   do iz=pm%tile1(3),pm%tile2(3) ! interpolate particles to density field
      do iy=pm%tile1(2),pm%tile2(2)
         do ix=pm%tile1(1),pm%tile2(1)
            rho_ex=0
            do ilayer=0,nlayer-1
               !$omp paralleldo default(shared) schedule(dynamic)&
               !$omp& private(k,j,i,np,nzero,l,ip,xpos,idx,l3,dx,i3,i2,i1)
               do k=pm%nc1(3)-pm%nloop+ilayer,pm%nc2(3)+pm%nloop,nlayer
                  !do k=pm%nc1(3)-pm%nloop,pm%nc2(3)+pm%nloop
                  do j=pm%nc1(2)-pm%nloop,pm%nc2(2)+pm%nloop
                     do i=pm%nc1(1)-pm%nloop,pm%nc2(1)+pm%nloop
                        np=rhoc(i,j,k,ix,iy,iz)
                        nzero=idx_b_r(j,k,ix,iy,iz)-sum(rhoc(i:,j,k,ix,iy,iz))
                        do l=1,np
                           ip=nzero+l
                           xpos=([i,j,k]-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
                           xpos=xpos*ratio_cs/pm%gridsize-pm%utile_shift*([ifx,ify,ifz]-1)*nfp
                           idx(:,2)=floor(xpos)+1
                           idx(:,1)=idx(:,2)-1
                           idx(:,3)=idx(:,2)+1
                           l3=xpos-floor(xpos)
                           dx(:,1)=(1-l3)**2/2
                           dx(:,3)=l3**2/2
                           dx(:,2)=1-dx(:,1)-dx(:,3)
                           if (minval(idx)<1-pm%nex .or. maxval(idx)>pm%nphy+pm%nex) then
                              print*,xpos
                              print*,idx
                              error stop
                           endif
                           do i3=1,p+1
                              do i2=1,p+1
                                 do i1=1,p+1
                                    rho_ex(idx(1,i1),idx(2,i2),idx(3,i3))=rho_ex(idx(1,i1),idx(2,i2),idx(3,i3))+dx(1,i1)*dx(2,i2)*dx(3,i3)
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
               !$omp endparalleldo
            enddo ! ilayer
            !print*,'    from rho_ex',pm%m1,pm%m2
            n1(:)=pm%nstart+pm%tile_shift*([ix,iy,iz]-1)*pm%nphy
            n2(:)=n1(:)+pm%nwork-1
            !print*,'     to density (',n1(1),':',n2(1),',',n1(2),':',n2(2),',',n1(3),':',n2(3),')'
            density(n1(1):n2(1),n1(2):n2(2),n1(3):n2(3))=rho_ex(pm%m1:pm%m2,pm%m1:pm%m2,pm%m1:pm%m2)
            testrho=testrho+sum(1d0*rho_ex(1:pm%nphy,1:pm%nphy,1:pm%nphy))
         enddo
      enddo
   enddo
   deallocate(rho_ex)
   !print*,testrho,sum(density*1d0)
   open(11,file='../../visualization/plane/rhoc.bin',status='replace',access='stream')
   write(11) rho1_iso(:2*nc,:,:)
   close(11)

   selectcase(pm%pm_layer) ! solve Poisson Eq. in Fourier space, diff potential to get force field
    case(0) ! single PM on standard grid
      !call sfftw_execute(plan0)
      !rho0k=rho0k*Gk0
      !call sfftw_execute(iplan0)
      !rho0=rho0/real(ng)/real(ng)/real(ng)
      !allocate(phi(-2:ng+3,-2:ng+3,-2:ng+3)[nn,nn,*])
      !phi(1:ng,1:ng,1:ng)=rho0(:ng,:,:)
      !call buffer_potential(phi,ng)
      !allocate(force(ndim,0:pm%nforce+1,0:pm%nforce+1,0:pm%nforce+1))
      !call grad(force,phi)
      !deallocate(phi)
    case(1) ! coarse grid PM
      call sfftw_execute(plan1)
      rho1k_iso=rho1k_iso*Gk1_iso
      call sfftw_execute(iplan1)
      rho1_iso=rho1_iso/real(2*nc)/real(2*nc)/real(2*nc)
      open(11,file='../../visualization/plane/phi1.bin',status='replace',access='stream')
      write(11) rho1_iso(:2*nc,:,:)
      close(11)
      allocate(force(ndim,3:pm%nforce-2,3:pm%nforce-2,3:pm%nforce-2))
      call grad_iso(force,rho1_iso)
      open(11,file='../../visualization/plane/force1.bin',status='replace',access='stream')
      write(11) force
      close(11)
    case(2) ! standard grid PM on extended tile
      call sfftw_execute(plan2)
      rho2k=rho2k*Gk2
      call sfftw_execute(iplan2)
      rho2=rho2/real(ngt)/real(ngt)/real(ngt)
      allocate(force(ndim,0:pm%nforce+1,0:pm%nforce+1,0:pm%nforce+1))
      call grad(force,rho2)
    case(3) ! fine grid APM on subtile
      call sfftw_execute(plan3)
      rho3k=rho3k*Gk3
      call sfftw_execute(iplan3)
      rho3=rho3/real(nft)/real(nft)/real(nft)
      allocate(force(ndim,0:pm%nforce+1,0:pm%nforce+1,0:pm%nforce+1))
      call grad(force,rho3)
   endselect
   !pm%f2max=max(pm%f2max,maxval(sum(force**2,1)))

   do iz=pm%tile1(3),pm%tile2(3) ! update velocity
      do iy=pm%tile1(2),pm%tile2(2)
         do ix=pm%tile1(1),pm%tile2(1)
            !$omp paralleldo default(shared) schedule(dynamic)&
            !$omp& private(k,j,i,np,nzero,l,ip,xpos,idx,l3,dx,vreal,dv,i3,i2,i1) &
            !$omp& reduction(max:vmax)
            do k=pm%nc1(3),pm%nc2(3)
               do j=pm%nc1(2),pm%nc2(2)
                  do i=pm%nc1(1),pm%nc2(1)
                     np=rhoc(i,j,k,ix,iy,iz)
                     nzero=idx_b_r(j,k,ix,iy,iz)-sum(rhoc(i:,j,k,ix,iy,iz))
                     do l=1,np ! loop over particle
                        ip=nzero+l
                        xpos=pm%tile_shift*([ix,iy,iz]-1)*nt+([i,j,k]-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
                        xpos=xpos*ratio_cs/pm%gridsize-pm%utile_shift*([ifx,ify,ifz]-1)*nfp
                        idx(:,2)=floor(xpos)+1
                        idx(:,1)=idx(:,2)-1
                        idx(:,3)=idx(:,2)+1
                        l3=xpos-floor(xpos)
                        dx(:,1)=(1-l3)**2/2
                        dx(:,3)=l3**2/2
                        dx(:,2)=1-dx(:,1)-dx(:,3)
                        !vreal=tan(pi*real(vp(:,ip))/real(nvbin-1))/(sqrt(pi/2)/(pm%sigv1*vrel_boost))
                        dv=0
                        if (minval(idx)<0 .or. maxval(idx)>pm%nforce+1) then
                           print*,xpos
                           print*,idx
                           error stop
                        endif
                        do i3=1,p+1
                           do i2=1,p+1
                              do i1=1,p+1
                                 dv=dv+force(:,idx(1,i1),idx(2,i2),idx(3,i3))*dx(1,i1)*dx(2,i2)*dx(3,i3)
                              enddo
                           enddo
                        enddo
                        !vreal=vreal+dv*a_mid*dt/6/pi
                        !vmax=max(vmax,abs(vreal+vfield(:,i,j,k,ix,iy,iz)))
                        !vp(:,ip)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(pm%sigv2*vrel_boost)*vreal)/pi,kind=izipv)
                        vp(:,ip)=vp(:,ip)+dv
                        !print*,ip,ip2iw(ip)
                        vp_store(:,ip2iw(ip),pm%pm_layer)=dv
                     enddo
                  enddo
               enddo
            enddo
            !$omp endparalleldo
         enddo
      enddo
   enddo
   deallocate(force)
endsubroutine

subroutine grad(force,phi)
   use omp_lib
   use variables
   implicit none
   save
   integer i_0,j_0,k_0,i_n(4),j_n(4),k_n(4)
   real force(ndim,0:pm%nforce+1,0:pm%nforce+1,0:pm%nforce+1)
   real phi(pm%m1phi(1):pm%m2phi(1),pm%m1phi(2):pm%m2phi(2),pm%m1phi(3):pm%m2phi(3))
   !$omp paralleldo default(shared) schedule(dynamic)&
   !$omp& private(k_0,j_0,i_0,i_n,j_n,k_n)
   do k_0=0,pm%nforce+1
      do j_0=0,pm%nforce+1
         do i_0=0,pm%nforce+1
            i_n=i_0+[-2,-1,1,2]
            j_n=j_0+[-2,-1,1,2]
            k_n=k_0+[-2,-1,1,2]
            force(1,i_0,j_0,k_0)=sum(phi(i_n,j_0,k_0)*weight)
            force(2,i_0,j_0,k_0)=sum(phi(i_0,j_n,k_0)*weight)
            force(3,i_0,j_0,k_0)=sum(phi(i_0,j_0,k_n)*weight)
         enddo
      enddo
   enddo
   !$omp endparalleldo
endsubroutine

subroutine grad_iso(force,phi)
   use omp_lib
   use variables
   implicit none
   save
   integer i_0,j_0,k_0,i_n(4),j_n(4),k_n(4)
   real force(ndim,3:pm%nforce-2,3:pm%nforce-2,3:pm%nforce-2)
   real phi(pm%m1phi(1):pm%m2phi(1),pm%m1phi(2):pm%m2phi(2),pm%m1phi(3):pm%m2phi(3))
   !$omp paralleldo default(shared) schedule(dynamic)&
   !$omp& private(k_0,j_0,i_0,i_n,j_n,k_n)
   do k_0=3,pm%nforce-2
      do j_0=3,pm%nforce-2
         do i_0=3,pm%nforce-2
            i_n=i_0+[-2,-1,1,2]
            j_n=j_0+[-2,-1,1,2]
            k_n=k_0+[-2,-1,1,2]
            force(1,i_0,j_0,k_0)=sum(phi(i_n,j_0,k_0)*weight)
            force(2,i_0,j_0,k_0)=sum(phi(i_0,j_n,k_0)*weight)
            force(3,i_0,j_0,k_0)=sum(phi(i_0,j_0,k_n)*weight)
         enddo
      enddo
   enddo
   !$omp endparalleldo
endsubroutine
