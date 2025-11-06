module variables
  use omp_lib
  use parameters
  implicit none
  save

  !integer(8),parameter :: np_image=(nc*np_nc)**3*merge(2,1,body_centered_cubic) ! average number of particles per image
  !integer(8),parameter :: np_image_max=np_image*(nte*1./nt)**3*image_buffer
  integer(8) np_image,np_image_max,np_tile_max
  !integer(8),parameter :: np_tile_max=np_image/nnt**3*(nte*1./nt)**3*tile_buffer
  real,parameter :: dt_max=1
  real,parameter :: G_grid=1.0/6.0/pi

  integer istep,itx,ity,itz,ifx,ify,ifz,nplocal
  integer(8) plan1,plan2,plan3,iplan1,iplan2,iplan3,plan0,iplan0 ! FFT plans
  integer(4) tt1,tt2,ttt1,ttt2,ijk(3,n_neighbor)

  real dt[*],dt_old[*],dt_mid[*],dt_e,da[*],a_mid[*]
  real f2_max_pm1[*],f2_max_pm2[*],f2_max_pm3[*],f2_max_pp[*]
  real vmax(3),overhead_tile[*],overhead_image[*],sigma_vi,sigma_vi_new,svz(500,2),svr(100,2)
  real(8) testrho,std_vsim_c[*],std_vsim_res[*],std_vsim[*]
  
  real,allocatable :: Gk1_iso(:,:,:),Gk2(:,:,:),Gk3(:,:,:),Gk0(:,:,:) ! Green's functions
  integer(izipx),allocatable :: xp(:,:)[:,:,:],xp_new(:,:)
  real,allocatable :: vp(:,:)[:,:,:],vp_new(:,:)
  real,allocatable :: vp_store(:,:,:)
  integer,allocatable :: rhoc(:,:,:,:,:,:)[:,:,:],ip2iw(:)
  real(4) vfield(3,1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[nn,nn,*] ! cannot have >7 dims

  real rho0(ng+2,ng,ng)
  complex rho0k(ng/2+1,ng,ng)
  equivalence(rho0,rho0k)
  real rho1_iso(2*nc+2,2*nc,2*nc)
  complex rho1k_iso(nc+1,2*nc,2*nc)
  equivalence(rho1_iso,rho1k_iso)
  real rho2(1-ngb:ngp+ngb+2,1-ngb:ngp+ngb,1-ngb:ngp+ngb)
  complex rho2k(ngt/2+1,ngt,ngt)
  equivalence(rho2,rho2k)
  real rho3(1-nfb:nfp+nfb+2,1-nfb:nfp+nfb,1-nfb:nfp+nfb)
  complex rho3k(nft/2+1,nft,nft)
  equivalence(rho3,rho3k)

  integer(8),dimension(1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt),codimension[nn,nn,*] :: idx_b_l,idx_b_r
  integer(8),dimension(nt,nt,nnt,nnt,nnt),codimension[nn,nn,*] :: ppl0,pplr,pprl,ppr0,ppl,ppr
  
  type type_pm
    integer pm_layer
    integer nwork
    integer nstart(ndim),nend(ndim),tile1(ndim),tile2(ndim)
    integer nloop,nex,nphy
    real gridsize,f2max,sigv1,sigv2
    integer tile_shift,utile_shift,m1,m2,nforce
    integer m1phi(ndim),m2phi(ndim),nc1(ndim),nc2(ndim)
  endtype
  type(type_pm) pm

contains
  subroutine buffer_potential(pot,n)
    implicit none
    integer(8) n
    real pot(-2:n+3,-2:n+3,-2:n+3)[nn,nn,*]
    sync all
    pot(:0,:,:)=pot(n-2:n,:,:)[inx,icy,icz]
    pot(n+1:,:,:)=pot(1:3,:,:)[ipx,icy,icz]
    sync all
    pot(:,:0,:)=pot(:,n-2:n,:)[icx,iny,icz]
    pot(:,n+1:,:)=pot(:,1:3,:)[icx,ipy,icz]
    sync all
    pot(:,:,:0)=pot(:,:,n-2:n)[icx,icy,inz]
    pot(:,:,n+1:)=pot(:,:,1:3)[icx,icy,ipz]
  endsubroutine

  subroutine spine_tile(rhoce,idx_ex_r,pp_l,pp_r,ppe_l,ppe_r)
    !! make a particle index (cumulative sumation) on tile
    !! used in update_particle, initial_conditions
    !! input:
    !! rhoce -- particle number density on tile, with 2x buffer depth
    !! output:
    !! idx_ex_r -- last extended index on extended right boundary
    !! ppe_r -- last extended index on physical right boundary
    !! pp_l -- first physical index on physical left boundary
    !! pp_r -- last physical index on physical right boundary
    !! ppe_l -- first extended index on physical left boundary
    use omp_lib
    implicit none
    integer(4),intent(in) :: rhoce(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
    integer(8),dimension(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb),intent(out) :: idx_ex_r
    integer(8),dimension(nt,nt),intent(out) :: pp_l,pp_r,ppe_l,ppe_r
    integer(8) nsum,np_phy
    integer igz,igy
    ! spine := yz-plane to record cumulative sum
    nsum=0
    do igz=1-2*ncb,nt+2*ncb
    do igy=1-2*ncb,nt+2*ncb
      nsum=nsum+sum(rhoce(:,igy,igz))
      idx_ex_r(igy,igz)=nsum ! right index
    enddo
    enddo
    !$omp paralleldo default(shared) private(igz,igy)
    do igz=1,nt
    do igy=1,nt
      ppe_r(igy,igz)=idx_ex_r(igy,igz)-sum(rhoce(nt+1:,igy,igz))
    enddo
    enddo
    !$omp endparalleldo
    nsum=0
    do igz=1,nt
    do igy=1,nt
      pp_l(igy,igz)=nsum+1
      np_phy=sum(rhoce(1:nt,igy,igz))
      nsum=nsum+np_phy
      pp_r(igy,igz)=nsum
      ppe_l(igy,igz)=ppe_r(igy,igz)-np_phy+1
    enddo
    enddo
  endsubroutine

  subroutine spine_image(rhoc,idx_b_l,idx_b_r,ppe_l0,ppe_lr,ppe_rl,ppe_r0,ppl,ppr)
    !! make a particle index (cumulative sumation) on image
    !! used in buffer_grid
    !! input:
    !! rhoc -- particle number density on image, with 1x buffer depth
    !! output:
    !! idx_b_l -- zeroth extended index on extended left boundary
    !! idx_b_r -- last extended index on extended right boundary
    !! ppl -- zeroth physical index on physical left boundary
    !! ppr -- last physical index on physical right boundary
    !! ppe_r0 -- last extended index on physical right boundary
    !! ppe_l0 -- zeroth extended index on physical left boundary
    !! ppe_rl -- last extended index on inner right boundary
    !! ppe_lr -- zeroth extended index on inner left boundary
    implicit none
    integer(4),intent(in) :: rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)
    integer(8),dimension(1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt),intent(out) :: idx_b_l,idx_b_r
    integer(8),dimension(nt,nt,nnt,nnt,nnt),intent(out) :: ppe_l0,ppe_lr,ppe_rl,ppe_r0,ppl,ppr
    integer(8) nsum,nsum_p,np_phy,ihz,ihy,ihx,igz,igy,ctile_mass(nnt,nnt,nnt),ctile_mass_p(nnt,nnt,nnt)
    ! spine := yz-plane to record cumulative sum
    nsum=0;nsum_p=0
    do ihz=1,nnt ! sum cumulative tile mass first
    do ihy=1,nnt
    do ihx=1,nnt
      ctile_mass(ihx,ihy,ihz)=nsum
      ctile_mass_p(ihx,ihy,ihz)=nsum_p
      nsum=nsum+sum(rhoc(:,:,:,ihx,ihy,ihz))
      nsum_p=nsum_p+sum(rhoc(1:nt,1:nt,1:nt,ihx,ihy,ihz))
    enddo
    enddo
    enddo
    !$omp paralleldo default(shared) private(ihz,ihy,ihx,nsum,igz,igy)
    do ihz=1,nnt ! calculate extended spine cumulative index on both sides
    do ihy=1,nnt
    do ihx=1,nnt
      nsum=ctile_mass(ihx,ihy,ihz)
      do igz=1-ncb,nt+ncb
      do igy=1-ncb,nt+ncb
        idx_b_l(igy,igz,ihx,ihy,ihz)=nsum
        nsum=nsum+sum(rhoc(:,igy,igz,ihx,ihy,ihz))
        idx_b_r(igy,igz,ihx,ihy,ihz)=nsum
      enddo
      enddo
    enddo
    enddo
    enddo
    !$omp endparalleldo
    !$omp paralleldo default(shared) private(ihz,ihy,ihx,nsum_p,igz,igy,np_phy)
    do ihz=1,nnt ! calculate physical spine
    do ihy=1,nnt
    do ihx=1,nnt
      nsum_p=ctile_mass_p(ihx,ihy,ihz)
      do igz=1,nt
      do igy=1,nt
        ppl(igy,igz,ihx,ihy,ihz)=nsum_p
        np_phy=sum(rhoc(1:nt,igy,igz,ihx,ihy,ihz))
        nsum_p=nsum_p+np_phy
        ppr(igy,igz,ihx,ihy,ihz)=nsum_p

        ppe_r0(igy,igz,ihx,ihy,ihz)=idx_b_r(igy,igz,ihx,ihy,ihz)-sum(rhoc(nt+1:,igy,igz,ihx,ihy,ihz))
        ppe_rl(igy,igz,ihx,ihy,ihz)= ppe_r0(igy,igz,ihx,ihy,ihz)-sum(rhoc(nt-ncb+1:nt,igy,igz,ihx,ihy,ihz))

        ppe_l0(igy,igz,ihx,ihy,ihz)=idx_b_l(igy,igz,ihx,ihy,ihz)+sum(rhoc(:0,igy,igz,ihx,ihy,ihz))
        ppe_lr(igy,igz,ihx,ihy,ihz)=ppe_l0(igy,igz,ihx,ihy,ihz) +sum(rhoc(1:ncb,igy,igz,ihx,ihy,ihz))
      enddo
      enddo
    enddo
    enddo
    enddo
    !$omp endparalleldo
  endsubroutine

endmodule
