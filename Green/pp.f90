subroutine pp
   use omp_lib
   use variables
   implicit none
   save
   print*, 'PP',nrange,n_neighbor,app,appr
   call system_clock(t1,t_rate)
   f2_max_pp=0
   do itz=1,nnt
      do ity=1,nnt
         do itx=1,nnt
            print*,itx,ity,itz
            call PP_Force()
         enddo
      enddo
   enddo
   !sim%dt_pp = 1. / (sqrt(f2_max_pp)*a_mid*G_grid)
   call system_clock(t2,t_rate)
   print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'
endsubroutine

subroutine PP_Force
   use omp_lib
   use variables
   implicit none
   save
   integer,parameter :: mc=ceiling(2.*nrange/ratio_cs)
   real,parameter :: rsoft=0.3
   integer i,j,k,ilayer,np,l,idxhoc(ndim),idx(ndim),idxf,ipp,ipll1_size,ipll2_size,icellpp,icellpp2
   integer ip,nzero,ip_offset,ipll1,ipll2,nl
   real xvec(ndim),vreal(ndim),rvec(ndim),rhat(ndim),rmag,fpp(ndim)
   integer,allocatable :: ll(:),hoc(:,:,:),ipf(:),ip_list(:)
   real,allocatable :: xf(:,:),vf(:,:),af(:,:),x_tmp(:,:),force_tmp(:,:)

   !print*,'ll'
   !print*,sum(rhoc(0:nt+1,0:nt+1,0:nt+1,itx,ity,itz)),nplocal
   allocate(ll(nplocal),hoc(1-mc*ratio_cs:ngp+mc*ratio_cs,1-mc*ratio_cs:ngp+mc*ratio_cs,1-mc*ratio_cs:ngp+mc*ratio_cs))
   allocate(ipf(nplocal),xf(3,nplocal),vf(3,nplocal),af(3,nplocal))
   nl=maxval(rhoc(0:nt+1,0:nt+1,0:nt+1,itx,ity,itz))*n_neighbor
   !print*,'nl',nl
   allocate(x_tmp(3,nl),force_tmp(3,nl),ip_list(nl))
   ip_offset=idx_b_l(1-mc,1-mc,itx,ity,itz)+sum(rhoc(:-mc,1-mc,1-mc,itx,ity,itz))
   hoc=0; ll=0; idxf=0; af=0
   do k=0,nt+1
      do j=0,nt+1
         do i=0,nt+1
            np=rhoc(i,j,k,itx,ity,itz)
            nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
            do l=1,np
               ip=nzero+l
               xvec=ratio_cs*([i,j,k]-1)+ratio_cs*(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
               idxf=idxf+1 ! index in the xf
               xf(:,idxf)=xvec ! make temp list xf
               vf(:,idxf)=vp(:,ip) ! make temp list vf
               ipf(idxf)=ip-ip_offset ! record index of ip in xf order, not consecutive due to buffer depth
               idx=floor(xvec)+1
               ll(idxf)=hoc(idx(1),idx(2),idx(3))
               hoc(idx(1),idx(2),idx(3))=idxf
            enddo
         enddo
      enddo
   enddo

   !print*,'  calculate force'
   ! update force
   do ilayer=0,nrange
      !$omp paralleldo default(shared) schedule(dynamic)&
      !$omp& private(k,j,i,ipll1_size,ipll1,ip_list)&
      !$omp& private(ipll2_size,icellpp,idxhoc,ipll2)&
      !$omp& private(x_tmp,force_tmp,icellpp2,rvec,rmag,fpp)
      do k=1-nrange+ilayer,ngp+nrange,nrange+1
         do j=1-nrange,ngp+nrange
            do i=1-nrange,ngp+nrange
               ipll1_size=0 ! list 1
               ipll1=hoc(i,j,k)
               do while (ipll1/=0)
                  ipll1_size=ipll1_size+1
                  ip_list(ipll1_size)=ipll1
                  ipll1=ll(ipll1)
               enddo
               if( ipll1_size==0 ) then
                  cycle
               endif
               !print*,'ipll1_size',ipll1_size
               ipll2_size=ipll1_size ! list 2
               do icellpp=1,n_neighbor
                  !print*,ijk(:,icellpp)
                  idxhoc=[i,j,k]+ijk(:,icellpp)
                  if (maxval(idxhoc)>ngp+ngb .or. minval(idxhoc)<1-ngb) then
                     print*,[i,j,k],idxhoc
                     error stop
                  endif
                  !print*,'find',idxhoc(1),idxhoc(2),idxhoc(3)
                  !print*,hoc(idxhoc(1),idxhoc(2),idxhoc(3))
                  ipll2=hoc(idxhoc(1),idxhoc(2),idxhoc(3))
                  !print*,ipll2
                  !print*,'a'
                  do while (ipll2/=0)
                     !print*,'发现'
                     !print*,ipll2,ll(ipll2),ipll2_size,icellpp
                     ipll2_size=ipll2_size+1
                     ip_list(ipll2_size)=ipll2
                     ipll2=ll(ipll2)
                     !print*,ipll2
                  enddo

               enddo
               !print*,'ipll2_size',ipll2_size
               !stop
               x_tmp(:,1:ipll2_size)=xf(:,ip_list(1:ipll2_size))
               !print*, ipll1_size
               !print*, ipll2_size
               !print*, ip_list(1:ipll2_size)
               !print*, x_tmp(:,1:ipll2_size)
               !stop
               force_tmp(:,:ipll2_size)=0
               do icellpp=1,ipll1_size
                  do icellpp2=icellpp+1,ipll2_size
                     rvec=x_tmp(:,icellpp2)-x_tmp(:,icellpp)
                     rmag=norm2(rvec)
                     rhat=rvec/rmag

                     !print*,icellpp,icellpp2
                     !print*,rvec
                     !print*,rmag,app,appr
                     !print*,F_ra(rmag,app),F_ra(rmag,appr)
                     !print*,''
                     !stop

                     if (rmag<pp_range) then
                        fpp=rhat*(F_ra(rmag,app)-F_ra(rmag,appr))
                        force_tmp(:,icellpp)=force_tmp(:,icellpp)+fpp
                        force_tmp(:,icellpp2)=force_tmp(:,icellpp2)-fpp
                     endif
                  enddo
               enddo
               do icellpp=1,ipll2_size
                  af(:,ip_list(icellpp))=af(:,ip_list(icellpp))+force_tmp(:,icellpp)
               enddo
            enddo
         enddo
      enddo
      !$omp endparalleldo
   enddo

   !do k=1,ngp
   !do j=1,ngp
   !do i=1,ngp
   !  ipll1=hoc(i,j,k)
   !  do while (ipll1/=0)
   !print*,i,j,k,ipll1
   !print*,af(:,ipll1); stop
   !    f2_max_pp=max(f2_max_pp,sum(af(:,ipll1)**2))
   !    ipll1=ll(ipll1)
   !  enddo
   !enddo
   !enddo
   !enddo

   !print*,'  update velocity'
   do ipp=1,idxf ! update vp
      ip=ipf(ipp)+ip_offset
      !vreal=vf(:,ipp)+af(:,ipp)*a_mid*dt/6/pi
      !vp(:,ip)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(sigma_vi*vrel_boost)*vreal)/pi,kind=izipv)
      vp(:,ip)=vp(:,ip)+af(:,ipp)
      vp_store(:,ip2iw(ip),4)=af(:,ipp)
   enddo
   deallocate(ll,hoc,ipf,xf,vf,af)
   deallocate(x_tmp,force_tmp,ip_list)
endsubroutine
