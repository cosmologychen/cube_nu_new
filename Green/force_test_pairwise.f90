use omp_lib
use parameters
use,intrinsic :: ISO_C_BINDING
implicit none
include 'fftw3.f03'

integer,parameter :: nmass=50
integer,parameter :: ntest=2000
integer,parameter :: npair=nmass*ntest
integer,parameter :: Interlace=0

real,allocatable :: Gk1(:,:,:),Gk2(:,:,:),Gk3(:,:,:),rho1(:,:,:),rho2(:,:,:),rho3(:,:,:)
real,allocatable :: xpos(:,:),testpos(:,:),rmag(:),force_correct(:,:)
real,allocatable :: force_PM1(:,:),force_PM2(:,:),force_PM3(:,:),force_PP(:,:)

integer(8) plan1,plan2,plan3,iplan1,iplan2,iplan3

call geometry
call omp_set_num_threads(ncore)
if (head) then
  print*,'kernel_initialization on',int(ncore,kind=2),' cores'
  print*,'  nc,ngt,nft =',nc,ngt,nft
endif

allocate(Gk1(nc/2+1,nc,nc),Gk2(ngt/2+1,ngt,ngt),Gk3(nft/2+1,nft,nft))

call Green_3D(Gk1,nc,nc/2+1,nc,nc,apm1c,0.,real(ratio_cs))
call Green_3D(Gk2,ngt,ngt/2+1,ngt,ngt,apm2,apm1,1.)
call Green_3D(Gk3,nft,nft/2+1,nft,nft,apm3f,apm2f,1./ratio_sf)

allocate(xpos(ndim,nmass),testpos(ndim,ntest))
call generate_particles

allocate(force_PM1(ndim,npair),force_PM2(ndim,npair),force_PM3(ndim,npair),force_PP(ndim,npair))
allocate(rho1(nc+2,nc,nc),rho2(ngt+2,ngt,ngt),rho3(nft+2,nft,nft))
allocate(rmag(npair),force_correct(ndim,npair))

print*,'create FFT plans'
call system_clock(t1,t_rate)
call sfftw_plan_dft_r2c_3d(plan1 , nc, nc, nc,rho1,rho1,FFTW_MEASURE)
call sfftw_plan_dft_c2r_3d(iplan1, nc, nc, nc,rho1,rho1,FFTW_MEASURE)
call sfftw_plan_dft_r2c_3d(plan2, ngt,ngt,ngt,rho2,rho2,FFTW_MEASURE)
call sfftw_plan_dft_c2r_3d(iplan2,ngt,ngt,ngt,rho2,rho2,FFTW_MEASURE)
call sfftw_plan_dft_r2c_3d(plan3, nft,nft,nft,rho3,rho3,FFTW_MEASURE)
call sfftw_plan_dft_c2r_3d(iplan3,nft,nft,nft,rho3,rho3,FFTW_MEASURE)
call system_clock(t2,t_rate)
print*,'    elapsed time =',real(t2-t1)/t_rate,'secs'; print*,''

print*,'PM'
call system_clock(t1,t_rate)
call Particle_Mesh(force_PM1,xpos,testpos,nc ,  0,  real(ratio_cs),Gk1,        0,rho1,plan1,iplan1)
call Particle_Mesh(force_PM2,xpos,testpos,ngp,ngb,             1.0,Gk2,Interlace,rho2,plan2,iplan2)
call Particle_Mesh(force_PM3,xpos,testpos,nfp,nfb,1/real(ratio_sf),Gk3,Interlace,rho3,plan3,iplan3)
call system_clock(t2,t_rate)
print*,'    elapsed time =',real(t2-t1)/t_rate,'secs'; print*,''

print*,'PP'
call system_clock(t1,t_rate)
call PP_force(force_PP,xpos,testpos,app,appr,pp_range)
call system_clock(t2,t_rate)
print*,'    elapsed time =',real(t2-t1)/t_rate,'secs'; print*,''

print*,'write in file'
open(11,file='force_data.bin',access='stream')
write(11) ndim,npair
write(11) rmag,force_PM1,force_PM2,force_PM3,force_PP,force_correct
close(11)

print*,force_PM1(:,1)


deallocate(Gk1,Gk2,Gk3,rho1,rho2,rho3,xpos,testpos)
deallocate(force_PM1,force_PM2,force_PM3,force_PP,rmag,force_correct)
call sfftw_destroy_plan(plan1)
call sfftw_destroy_plan(iplan1)
call sfftw_destroy_plan(plan2)
call sfftw_destroy_plan(iplan2)
call sfftw_destroy_plan(plan3)
call sfftw_destroy_plan(iplan3)




contains

  subroutine PP_force(force_PP,xpos,testpos,app,appr,pp_range)
    integer imass,itest,ipair
    real rvec(ndim),rhat(ndim),app,appr,pp_range
    real,allocatable :: force_PP(:,:),xpos(:,:),testpos(:,:)
    do imass=1,nmass
    do itest=1,ntest
        ipair=(imass-1)*ntest+itest
        rvec=xpos(:,imass)-testpos(:,itest);
        rvec=modulo(rvec+ng/2,real(ng))-ng/2;
        rmag(ipair)=norm2(rvec);
        rhat=rvec/rmag(ipair);
        if (rmag(ipair)<pp_range) then
            force_PP(:,ipair)=rhat*(F_ra(rmag(ipair),app)-F_ra(rmag(ipair),appr));
        endif
        force_correct(:,ipair)=rhat*F_ra(rmag(ipair),app);
    enddo
    enddo
  endsubroutine

  subroutine Particle_Mesh(force_PM,xpos,testpos,nphys,nbuffer,gridsize,Gk,Interlace,rho,plan,iplan)
    integer(8) nphys,plan,iplan
    integer nbuffer,Interlace,ngrid,idx(ndim,p+1),imass,i,j,k,in(4),jn(4),kn(4),itest,i_nnt(ndim),ipair
    real gridsize,dx(ndim,p+1),l3(ndim)
    real,allocatable :: force_PM(:,:),xpos(:,:),testpos(:,:),Gk(:,:,:),rho(:,:,:),force(:,:,:,:)

    print*,'particle mesh'
    ngrid=nphys+2*nbuffer
    print*,'  ngrid =',ngrid
    print*,'  shape(Gk)',shape(Gk)
    rho1=0
    allocate(force(ndim,ngrid,ngrid,ngrid))
    do imass=1,nmass
      idx(:,2)=nbuffer+floor(xpos(:,imass)/gridsize-Interlace*0.5)+1
      idx(:,1)=idx(:,2)-1
      idx(:,3)=idx(:,2)+1
      l3=(xpos(:,imass)/gridsize-Interlace*0.5)-floor(xpos(:,imass)/gridsize-Interlace*0.5)
      dx(:,1)=(1-l3)**2/2
      dx(:,3)=l3**2/2
      dx(:,2)=1-dx(:,1)-dx(:,3)
      idx=modulo(idx-1,ngrid)+1
      rho=0
      do k=1,p+1
      do j=1,p+1
      do i=1,p+1
        rho(idx(1,i),idx(2,j),idx(3,k))=rho(idx(1,i),idx(2,j),idx(3,k))+dx(1,i)*dx(2,j)*dx(3,k)
      enddo
      enddo
      enddo

      call sfftw_execute(plan)
      rho(::2,:,:)=rho(::2,:,:)*Gk
      rho(2::2,:,:)=rho(2::2,:,:)*Gk
      call sfftw_execute(iplan)
      rho=rho/real(ngrid)**3

      !$omp paralleldo default(shared) schedule(dynamic)&
      !$omp& private(k,j,i,in,jn,kn)
      do k=1,ngrid
      do j=1,ngrid
      do i=1,ngrid
        in=modulo(i-1+[-2,-1,1,2],ngrid)+1
        jn=modulo(j-1+[-2,-1,1,2],ngrid)+1
        kn=modulo(k-1+[-2,-1,1,2],ngrid)+1
        force(1,i,j,k)=sum(rho(in,j,k)*weight)
        force(2,i,j,k)=sum(rho(i,jn,k)*weight)
        force(3,i,j,k)=sum(rho(i,j,kn)*weight)
      enddo
      enddo
      enddo
      !$omp endparalleldo

      !$omp paralleldo default(shared) schedule(dynamic)&
      !$omp& private(itest,i_nnt,ipair,idx,l3,dx,k,j,i)
      do itest=1,ntest
        i_nnt=floor(testpos(:,itest) / (nphys*gridsize))+1
        ipair=(imass-1)*ntest+itest
        if (sum(i_nnt)==ndim) then
          idx(:,2) = nbuffer+floor(testpos(:,itest)/gridsize-Interlace*0.5)+1
          idx(:,1)=idx(:,2)-1
          idx(:,3)=idx(:,2)+1
          l3=(testpos(:,itest)/gridsize-Interlace*0.5)-floor(testpos(:,itest)/gridsize-Interlace*0.5)
          dx(:,1)=(1-l3)**2/2
          dx(:,3)=l3**2/2
          dx(:,2)=1-dx(:,1)-dx(:,3)
          idx=modulo(idx-1,ngrid)+1
          do k=1,p+1
          do j=1,p+1
          do i=1,p+1
              force_PM(1,ipair)=force_PM(1,ipair)+force(1,idx(1,i),idx(2,j),idx(3,k))*dx(1,i)*dx(2,j)*dx(3,k)
              force_PM(2,ipair)=force_PM(2,ipair)+force(2,idx(1,i),idx(2,j),idx(3,k))*dx(1,i)*dx(2,j)*dx(3,k)
              force_PM(3,ipair)=force_PM(3,ipair)+force(3,idx(1,i),idx(2,j),idx(3,k))*dx(1,i)*dx(2,j)*dx(3,k)
          enddo
          enddo
          enddo
        endif
      enddo
      !$omp endparalleldo

    enddo ! imass=1,nmass
    deallocate(force)
  endsubroutine

  subroutine generate_particles
    integer seedsize,i
    integer,allocatable :: iseed(:)
    real temp_r(ndim),temp_theta(ndim),x_cen(ndim)

    call system_clock(t1,t_rate)
    print*,'generate_particles'
    call random_seed(size=seedsize)
    print*,'  seedsize =',seedsize
    allocate(iseed(seedsize))
    iseed=777
    call random_seed(put=iseed)
    call random_number(xpos)
    call random_number(testpos)
    deallocate(iseed)
    x_cen=real(nfp)/ratio_sf/2
    xpos=spread(x_cen,2,nmass)+(xpos-0.5)*ratio_cs
    ! Box-Muller transform
    !$omp paralleldo default(shared) schedule(dynamic)&
    !$omp& private(i,temp_theta,temp_r)
    do i=1,ntest,2
      temp_theta=2*pi*testpos(:,i)
      temp_r=sqrt(-2*log(1-testpos(:,i+1)))
      testpos(:,i)=temp_r*cos(temp_theta)
      testpos(:,i+1)=temp_r*sin(temp_theta)
    enddo
    !$omp endparalleldo
    ! set mean and std
    testpos(:,:ntest/2)  = spread(x_cen,2,ntest/2)+apm2*testpos(:,:ntest/2)
    testpos(:,ntest/2+1:)= spread(x_cen,2,ntest/2)+apm1*testpos(:,ntest/2+1:)
    testpos=modulo(testpos,real(ng))
    !print*,shape(spread([1,2,3],2,4))
    !print*,spread([1,2,3],2,4)
    call system_clock(t2,t_rate)
    print*,'    elapsed time =',real(t2-t1)/t_rate,'secs';
    print*,''
  endsubroutine

end