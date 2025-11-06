subroutine kick
  use omp_lib
  use variables
  implicit none

  logical,parameter :: PM_12=.false.
  integer i
  ! single layer PM

  print*,'kick'

  if (PM_12) then
    !pm%pm_layer=0
    !pm%nwork=ngp
    !pm%nstart(:)=1
    !pm%nend(:)=ng; pm%nend(1)=pm%nend(1)+2
    !pm%tile1(:)=1
    !pm%tile2(:)=nnt
    !pm%nloop=1
    !pm%nex=5
    !pm%nphy=ngp
    !pm%gridsize=1
    !pm%tile_shift=1
    !pm%utile_shift=0
    !pm%m1=1; pm%m2=ngp
    !pm%nforce=ng
    !pm%m1phi(:)=-2; pm%m2phi(:)=ng+3;
    !pm%sigv1=sigma_vi
    !pm%sigv2=sigma_vi_new
    !pm%nc1(:)=1
    !pm%nc2(:)=nt
    !testrho=0
    !call apm(rho0)
    !sigma_vi=sigma_vi_new
  else
    print*,'PM1'
    call system_clock(t1,t_rate)
    pm%pm_layer=1
    pm%nwork=nt
    pm%nstart(:)=1
    pm%nend(:)=2*nc; pm%nend(1)=pm%nend(1)+2
    pm%tile1(:)=1
    pm%tile2(:)=nnt
    pm%nloop=1
    pm%nex=2
    pm%nphy=nt
    pm%gridsize=ratio_cs
    pm%tile_shift=1
    pm%utile_shift=0
    pm%m1=1; pm%m2=nt
    pm%nforce=nc
    pm%m1phi(:)=1; pm%m2phi(:)=2*nc; pm%m2phi(1)=pm%m2phi(1)+2
    pm%f2max=0
    pm%sigv1=sigma_vi
    pm%sigv2=sigma_vi
    pm%nc1(:)=1
    pm%nc2(:)=nt
    testrho=0
    call apm(rho1_iso)
    sim%dt_pm1=sqrt( 1. / (sqrt(pm%f2max)*a_mid*G_grid) )
    print*,'  pm%f2max =',pm%f2max
    print*,'  testrho =',testrho
    call system_clock(t2,t_rate)
    print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';

    !return

    print*,'PM2'
    vmax=0
    f2_max_pm2=0
    call system_clock(t1,t_rate)
    testrho=0; pm%f2max=0
    do itz=1,nnt
    do ity=1,nnt
    do itx=1,nnt
      pm%pm_layer=2
      pm%nwork=ngt
      pm%nstart(:)=1-ngb
      pm%nend(:)=ngp+ngb; pm%nend(1)=pm%nend(1)+2
      pm%tile1(:)=[itx,ity,itz]
      pm%tile2(:)=[itx,ity,itz]
      pm%nloop=ncb
      pm%nex=ngb+1
      pm%nphy=ngp
      pm%gridsize=1
      pm%tile_shift=0
      pm%utile_shift=0
      pm%m1=1-ngb; pm%m2=ngp+ngb
      pm%nforce=ngp
      pm%m1phi(:)=1-ngb; pm%m2phi(:)=ngp+ngb; pm%m2phi(1)=pm%m2phi(1)+2
      pm%sigv1=sigma_vi
      pm%sigv2=sigma_vi_new
      pm%nc1(:)=1
      pm%nc2(:)=nt
      call apm(rho2)
    enddo
    enddo
    enddo ! itz
    sigma_vi=sigma_vi_new
    sim%vsim2phys=(1.5/sim%a)*box*h0*100.*sqrt(omega_m)/ng_global
    sim%dt_pm2=sqrt( 1. / (sqrt(pm%f2max)*a_mid*G_grid) )
    call system_clock(t2,t_rate)
    print*,'  sum of PM2 density',testrho
    print*,'  elapsed time =',real(t2-t1)/t_rate,'secs'
  endif


  !return



  print*, 'PM3'
  call system_clock(t1,t_rate)
  testrho=0; f2_max_pm3=0
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt  
    do ifz=1,nns
    do ify=1,nns
    do ifx=1,nns
      call system_clock(tt1,t_rate)
      !call PM_ultrafine
      pm%pm_layer=3
      pm%nwork=nft
      pm%nstart(:)=1-nfb
      pm%nend(:)=nfp+nfb; pm%nend(1)=pm%nend(1)+2
      pm%tile1(:)=[itx,ity,itz]
      pm%tile2(:)=[itx,ity,itz]
      pm%nloop=1
      pm%nex=nfb+1
      pm%nphy=nfp
      pm%gridsize=1./ratio_sf
      pm%tile_shift=0
      pm%utile_shift=1
      pm%m1=1-nfb; pm%m2=nfp+nfb
      pm%nforce=nfp
      pm%m1phi(:)=1-nfb; pm%m2phi(:)=nfp+nfb; pm%m2phi(1)=pm%m2phi(1)+2
      pm%sigv1=sigma_vi
      pm%sigv2=sigma_vi
      pm%nc1(:)=([ifx,ify,ifz]-1)*ncfp+1
      pm%nc2(:)=[ifx,ify,ifz]*ncfp
      call apm(rho3)
      call system_clock(tt2,t_rate)
      !print*,'  elapsed time =',real(tt2-tt1)/t_rate,'secs'
    enddo
    enddo
    enddo
  enddo
  enddo
  enddo
  sim%dt_pm3=sqrt( 1. / (sqrt(pm%f2max)*a_mid*G_grid) )
  print*,'  sum of PM3 density',testrho
  print*,'  dt_pm3 =',sim%dt_pm3
  call system_clock(t2,t_rate)
  print*,'  elapsed time =',real(t2-t1)/t_rate,'secs'

  call system_clock(t1,t_rate)
  call pp
  call system_clock(t2,t_rate)
  print*,'    elapsed time =',real(t2-t1)/t_rate,'secs'; print*,''

  !sim%dt_vmax=vbuf*20/maxval(vmax)
  !sim%vz_max=vmax(3)
  !sync all
  !do i=1,nn**3
  !  sim%dt_pm2=min(sim%dt_pm2,sim[i]%dt_pm2)
  !  sim%dt_pm1=min(sim%dt_pm1,sim[i]%dt_pm1)
  !  sim%dt_vmax=min(sim%dt_vmax,sim[i]%dt_vmax)
  !enddo

endsubroutine