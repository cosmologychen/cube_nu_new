!#define gadget

program CUBE_FoF
  use omp_lib
  use variables
  implicit none
  character(len = 4) str_refine
  ! integer(8),parameter:: n_refine = 5
  integer(8),parameter:: fof_buffer = ceiling(10.0/box*nc*nn)
  ! integer(8),parameter:: nfof = nc+2*fof_buffer
  ! integer(8),parameter:: nfof2 = nfof/2
  integer i,j,k,l,cur_checkpoint,np,n1,n2,idx(3),nh[*],im,idx1(3),idx2(3)
  integer iq1,iq2,iq3,i_neighbor,jq(3),ijk_neighbor(3,3),neighbor_b(3,3),nfofs(4),nfof,ifof,nfof1,n_layer
  integer(4) nlast,ip,jp,np_iso,np_head,np_max,np_neighbors(3,3,3)[nn,nn,*],np_need(3,3,3),max_nei[nn,nn,*]
  integer(4),allocatable :: np_halo_all(:),np_halo(:)[:]
  integer(4),allocatable :: hoc(:,:,:),ll(:),llgp(:),hcgp(:),ecgp(:),iph_halo_all(:),iph_halo(:)!,h2(:,:,:,:),l2(:,:)
  real rp2,rsq,pos1(3),dx1(3),dx2(3),shift_xv(3)
  real(8) rho8,n_refine,L_b,L_bL,L_bLb,L_bLb2,dxv(3)
  real,allocatable :: xv(:,:),xv_mean(:,:),xp_neighbors(:,:,:,:,:)[:,:,:],rho_grid(:,:,:)[:,:,:]
  type(type_halo_catalog_header) halo_header
  type(type_halo_catalog_array),allocatable :: hcat(:)[:]
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call geometry

  call omp_set_num_threads(1)

  if (head) then
    open(16,file='z_checkpoint.txt',status='old')
    do i=1,nmax_redshift-1
      read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
    enddo
    71 n_checkpoint=i-1
    close(16)
    if (n_checkpoint==0) stop 'z_checkpoint.txt empty'
  endif
  sync all
  n_checkpoint=n_checkpoint[1]
  z_checkpoint(:)=z_checkpoint(:)[1]
  sync all

  if (head) print*,'  initialize FoF cell neighbors',fof_buffer
  l=0
  do j=-1,1
  do i=-1,1
    l=l+1; ijk(:,l)=[i,j,-1]
  enddo
  enddo
  do i=-1,1
    l=l+1; ijk(:,l)=[i,-1,0]
  enddo
  l=l+1; ijk(:,l)=[-1,0,0]

  if (head) print*,'  initialize tile index'
  ipm2=0;
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    ipm2=ipm2+1
    ixyz2(:,ipm2)=[itx,ity,itz]
  enddo
  enddo
  enddo
  ! stop
  rp2=(b_link/ratio_cs)**2
  L_b    = fof_buffer
  L_bL   = L_b  + nc
  L_bLb  = L_bL + fof_buffer
  L_bLb2 = L_bLb/2
  ! nfofs = [400,800,1200,1600]
  nfofs = [300,400,500,800]
  shift_xv  = [(icx-1)*nc*1d0-fof_buffer,(icy-1)*nc*1d0-fof_buffer,(icz-1)*nc*1d0-fof_buffer]
  if (head) print*,'linking,L _b,L_bL,L_bLb',b_link/ratio_cs,L_b,L_bL,L_bLb

  do cur_checkpoint=n_checkpoint,n_checkpoint
  ! do cur_checkpoint=4,4
    sim%cur_checkpoint=cur_checkpoint
    if (head) print*, ''
    if (head) print*, ''
    if (head) print*, 'FoF at redshift ',z2str(z_checkpoint(cur_checkpoint))
    if (head) print*, '  read checkpoint header',output_name('info')

    ! load particle
    open(11,file=output_name('info'),access='stream'); read(11) sim; close(11)
    allocate(rhoc(nt,nt,nt,nnt,nnt,nnt)[nn,nn,*],xp(3,sim%nplocal)[nn,nn,*])
    sim%cur_checkpoint=cur_checkpoint
    open(11,file=output_name('xp'),access='stream'); read(11) xp; close(11)
    open(11,file=output_name('np'),access='stream'); read(11) rhoc; close(11)

    neighbor_b = reshape([int(0,kind=8),fof_buffer,nc  ,int(0,kind=8),nc,int(0,kind=8) ,nc-fof_buffer,nc,-nc ],[3,3])

    ! count  neighbor particle  -所有buffer在ic是的密度一致
    if (head) print*,' count  neighbor particle '
    ! do im = 1,nn**3
    !   if (image == im ) then
    !     print*,image,sum(rhoc)*1d0/nc/nc/nc
    np_neighbors = 0
    do iq1 = 1, 3
    ijk_neighbor(:, 1) = neighbor_b(:, iq1)
    do iq2 = 1, 3
    ijk_neighbor(:, 2) = neighbor_b(:, iq2)
    do iq3 = 1, 3
      np_neighbors(iq1,iq2,iq3) = 0
      if (iq1 == 2 .and. iq2 == 2 .and. iq3 == 2) cycle
      ijk_neighbor(:, 3) = neighbor_b(:, iq3)
      ! 循环遍历各个维度
      do itz = floor((ijk_neighbor(1, 3))*1d0/nt)+1, floor((ijk_neighbor(2, 3)-1)*1d0/nt)+1
      do ity = floor((ijk_neighbor(1, 2))*1d0/nt)+1, floor((ijk_neighbor(2, 2)-1)*1d0/nt)+1
      do itx = floor((ijk_neighbor(1, 1))*1d0/nt)+1, floor((ijk_neighbor(2, 1)-1)*1d0/nt)+1
        do k = merge(1, mod(ijk_neighbor(1, 3)+1, nt), itz*nt-nt >= ijk_neighbor(1, 3)),merge(nt, mod(ijk_neighbor(2, 3)+1, nt)-1, itz*nt <= ijk_neighbor(2, 3))
        do j = merge(1, mod(ijk_neighbor(1, 2)+1, nt), ity*nt-nt >= ijk_neighbor(1, 2)),merge(nt, mod(ijk_neighbor(2, 2)+1, nt)-1, ity*nt <= ijk_neighbor(2, 2))
        do i = merge(1, mod(ijk_neighbor(1, 1)+1, nt), itx*nt-nt >= ijk_neighbor(1, 1)),merge(nt, mod(ijk_neighbor(2, 1)+1, nt)-1, itx*nt <= ijk_neighbor(2, 1))
          np_neighbors(iq1,iq2,iq3) = np_neighbors(iq1,iq2,iq3) + rhoc(i,j,k,itx,ity,itz)
        enddo
        enddo
        enddo
      enddo
      enddo
      enddo
      ! if (head) print*,iq1-2,iq2-2,iq3-2,np_neighbors(iq1,iq2,iq3)!,np_neighbors(iq1,iq2,iq3)*1d0/(ijk_neighbor(2, 1) -  ijk_neighbor(1, 1))/ (ijk_neighbor(2, 2) -  ijk_neighbor(1, 2))/  (ijk_neighbor(2, 3) -  ijk_neighbor(1, 3))
    enddo
    enddo
    enddo
    max_nei = maxval(np_neighbors)
    
    ! init space for all particles -检查过image——neighbor编号和np_need
    sync all
    if (head) print*,' init space for all particles'
    np_need = 0
    np_max = sum(rhoc)
    do iq1 = 1,3
      i = modulo(icx-3+iq1,nn)+1
    do iq2 = 1,3
      j = modulo(icy-3+iq2,nn)+1
    do iq3 = 1,3
      k = modulo(icz-3+iq3,nn)+1
      if (iq1 == 2 .and. iq2 == 2 .and. iq3 == 2) cycle
      np_need(iq1,iq2,iq3) = np_neighbors(4-iq1,4-iq2,4-iq3)[i,j,k]
      ! print*,iq1-2,iq2-2,iq3-2,np_need(iq1,iq2,iq3)
    enddo
    enddo
    enddo
    np_max = np_max + sum(np_need)
    sync all
    if (head) then
      print*,'sync  max_nei'
      do i = 1,nn
      do j = 1,nn
      do k = 1,nn
        max_nei = max(max_nei,max_nei[i,j,k])
      enddo
      enddo
      enddo
    endif
    sync all
    max_nei = max_nei[1,1,1]
    do im = 1,nn**3
      if (image == im ) then
        print*,image,np_max
      endif
      sync all
    enddo
    stop

    !  initialize particles  -检查过xv范围和sum得到nlast
    if (head) print*,' initialize particles'
    allocate(xv(3,np_max))
    xv = 0
    nlast=0
    !$omp paralleldo  COLLAPSE(3) PRIVATE(itz,ity,itx,k,j,i,np,nlast,l,ip)
    do itz=1,nnt
    do ity=1,nnt
    do itx=1,nnt
      do k=1,nt
      do j=1,nt
      do i=1,nt
        np=rhoc(i,j,k,itx,ity,itz)
        nlast = sum(rhoc(:,:,:,:,:,:itz-1))     &
              + sum(rhoc(:,:,:,:,:ity-1,itz))   &
              + sum(rhoc(:,:,:,:itx-1,ity,itz)) &
              + sum(rhoc(:,:,:k-1,itx,ity,itz)) &
              + sum(rhoc(:,:j-1,k,itx,ity,itz)) &
              + sum(rhoc(:i-1,j,k,itx,ity,itz))
        do l=1,np
          ip=nlast+l
          if (ip < 1 .or.  ip > sim%nplocal) then
            print*, image
            print*, i,j,k,itx,ity,itz
            print*,nlast,l,np,ip
            stop
          endif
#ifdef ZIPX
          xv(:,ip)=nt*((/itx,ity,itz/)-1)+ ((/i,j,k/)-1) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution+[fof_buffer,fof_buffer,fof_buffer] 
#else 
          xv(:,ip)=xp(:,ip)/ratio_cs+[fof_buffer,fof_buffer,fof_buffer] 
#endif
        enddo
      enddo
      enddo
      enddo
      if(head) print*,itz,ity,itx,sum(rhoc(:,:,:,itx,ity,itz))
    enddo
    enddo
    enddo
    !$omp endparalleldo
    deallocate(xp)
    np_max = sum(rhoc)
    ! do im = 1,nn**3
    !   if (image == im ) then
    !     print*,image,sum(xv(1,:np_max))/np_max,minval(xv(1,:np_max)),maxval(xv(1,:np_max))
    !     print*,image,sum(xv(2,:np_max))/np_max,minval(xv(2,:np_max)),maxval(xv(2,:np_max))
    !     print*,image,sum(xv(3,:np_max))/np_max,minval(xv(3,:np_max)),maxval(xv(3,:np_max))
    !   endif
    !   sync all
    ! enddo
    ! stop


    ! init particles for neighbors  -检查粒子的中心位置和最大最小值
    if (head) print*,' init particles for neighbors',max_nei,maxval(np_neighbors)

    allocate(xp_neighbors(3,max_nei,3,3,3)[nn,nn,*])
    xp_neighbors = 0
    ! do im = 1,nn**3
    !   if (image == im ) then
    !     print*,image
    !$omp paralleldo  COLLAPSE(3) PRIVATE(iq1,iq2,iq3,jp,ijk_neighbor,itz,ity,itx,k,j,i,nlast,np,l,ip)
    do iq1 = 1, 3
    do iq2 = 1, 3
    do iq3 = 1, 3
      if (iq1 == 2 .and. iq2 == 2 .and. iq3 == 2) cycle
      jp = 0
      ijk_neighbor(:, 2) = neighbor_b(:, iq2)
      ijk_neighbor(:, 1) = neighbor_b(:, iq1)
      ijk_neighbor(:, 3) = neighbor_b(:, iq3)
      ! 循环遍历各个维度
      do itz = floor((ijk_neighbor(1, 3))*1d0/nt)+1, floor((ijk_neighbor(2, 3)-1)*1d0/nt)+1
      do ity = floor((ijk_neighbor(1, 2))*1d0/nt)+1, floor((ijk_neighbor(2, 2)-1)*1d0/nt)+1
      do itx = floor((ijk_neighbor(1, 1))*1d0/nt)+1, floor((ijk_neighbor(2, 1)-1)*1d0/nt)+1
        do k = merge(1, mod(ijk_neighbor(1, 3)+1, nt), itz*nt-nt >= ijk_neighbor(1, 3)),merge(nt, mod(ijk_neighbor(2, 3)+1, nt)-1, itz*nt <= ijk_neighbor(2, 3))
        do j = merge(1, mod(ijk_neighbor(1, 2)+1, nt), ity*nt-nt >= ijk_neighbor(1, 2)),merge(nt, mod(ijk_neighbor(2, 2)+1, nt)-1, ity*nt <= ijk_neighbor(2, 2))
        do i = merge(1, mod(ijk_neighbor(1, 1)+1, nt), itx*nt-nt >= ijk_neighbor(1, 1)),merge(nt, mod(ijk_neighbor(2, 1)+1, nt)-1, itx*nt <= ijk_neighbor(2, 1))
          nlast = sum(rhoc(:,:,:,:,:,:itz-1))     &
                + sum(rhoc(:,:,:,:,:ity-1,itz))   &
                + sum(rhoc(:,:,:,:itx-1,ity,itz)) &
                + sum(rhoc(:,:,:k-1,itx,ity,itz)) &
                + sum(rhoc(:,:j-1,k,itx,ity,itz)) &
                + sum(rhoc(:i-1,j,k,itx,ity,itz))
          np=rhoc(i,j,k,itx,ity,itz)
          do l=1,np
              ip=nlast+l
              jp = jp + 1
              xp_neighbors(:,jp,iq1,iq2,iq3) = xv(:,ip) + ijk_neighbor(3,:)
          enddo
        enddo
        enddo
        enddo
      enddo
      enddo
      enddo
      if (head)  print*,iq1-2,iq2-2,iq3-2
      ! print*,iq1-2,iq2-2,iq3-2,np_neighbors(iq1,iq2,iq3),jp
      ! print*,sum(xp_neighbors(1,:jp,iq1,iq2,iq3))/jp,minval(xp_neighbors(1,:jp,iq1,iq2,iq3)),maxval(xp_neighbors(1,:jp,iq1,iq2,iq3)),ijk_neighbor(3,1)
      ! print*,sum(xp_neighbors(2,:jp,iq1,iq2,iq3))/jp,minval(xp_neighbors(2,:jp,iq1,iq2,iq3)),maxval(xp_neighbors(2,:jp,iq1,iq2,iq3)),ijk_neighbor(3,2)
      ! print*,sum(xp_neighbors(3,:jp,iq1,iq2,iq3))/jp,minval(xp_neighbors(3,:jp,iq1,iq2,iq3)),maxval(xp_neighbors(3,:jp,iq1,iq2,iq3)),ijk_neighbor(3,3)
    enddo
    enddo
    enddo
    !$omp endparalleldo
    !   endif
    !   sync all
    ! enddo
    ! stop
    


    ! get  particles form neighbors  -检查粒子的中心位置和最大最小值
    sync all
    if (head) print*,' get particles form neighbors',sim%nplocal, sum(rhoc)
    np_max = sim%nplocal
    do iq1 = 1,3
    do iq2 = 1,3
    do iq3 = 1,3
      j = modulo(icy-3+iq2,nn)+1
      i = modulo(icx-3+iq1,nn)+1
      k = modulo(icz-3+iq3,nn)+1
      if (iq1 == 2 .and. iq2 == 2 .and. iq3 == 2) cycle
      xv(:,np_max+1:np_max+np_need(iq1,iq2,iq3)) = xp_neighbors(:,:np_need(iq1,iq2,iq3),4-iq1,4-iq2,4-iq3)[i,j,k]
      np_max = np_max +  np_need(iq1,iq2,iq3)
      if (head)  print*,iq1-2,iq2-2,iq3-2
    enddo
    enddo
    enddo
    ! do im = 1,nn**3
    !   if (image == im ) then
    !     print*,image,np_max
    !     print*,sum(xv(1,:np_max))/np_max,minval(xv(1,:np_max)),maxval(xv(1,:np_max))
    !     print*,sum(xv(2,:np_max))/np_max,minval(xv(2,:np_max)),maxval(xv(2,:np_max))
    !     print*,sum(xv(3,:np_max))/np_max,minval(xv(3,:np_max)),maxval(xv(3,:np_max))
    !   endif
    !   sync all
    ! enddo
    ! stop

    sync all

    do ifof=1,4
      call tic(11)
      nfof=nfofs(ifof)
      ! nfof=nfofs(1)
      nfof1 = nfof+1
      n_refine = nfof*1d0/(nc+fof_buffer*2+b_link/ratio_cs)!-(i-1)*0.1
      n_layer = ceiling(fof_buffer*n_refine)+1
      write(str_refine,'(i4)')  nfof
      if(head) print*,ifof,nfof,n_refine,4/n_refine
      
      ! create hoc ll
      if (head) print*,' create hoc ll'
      allocate(hoc(nfof,nfof,nfof),ll(np_max),llgp(np_max),hcgp(np_max),ecgp(np_max))
      hoc=0; ll=0
      do ip=1,np_max
        idx=floor(xv(1:3,ip)*n_refine)+1 ! index of the grid
        ll(ip)=hoc(idx(1),idx(2),idx(3)) ! linked list 
        hoc(idx(1),idx(2),idx(3))=ip ! head of chain
        hcgp(ip)=ip ! initialize hcgp(ip)=ip for isolated particles
      enddo
      llgp=0; ecgp=0; ! initialize group link list
      if (head) print*, 'hoc ll created'
      
      ! rho
      if (1) then
        allocate(rho_grid(0:nfof+1,0:nfof+1,0:nfof+1)[nn,nn,*])
        rho_grid=0
        do l=0,n_layer-1
        !$omp paralleldo  PRIVATE(iq3,iq2,iq1,ip,pos1,idx1,idx2,dx1,dx2)
        do iq3=1+l,nfof,n_layer
        do iq2=1,nfof
        do iq1=1,nfof
          ip=hoc(iq1,iq2,iq3)
          do while (ip/=0)
            pos1 = xv(1:3,ip)
            if ( pos1(1) == L_bLb ) pos1(1)= 0
            if ( pos1(2) == L_bLb ) pos1(2)= 0
            if ( pos1(3) == L_bLb ) pos1(3)= 0
            pos1 = pos1*n_refine
            idx1=floor(pos1)+1; idx2=idx1+1
            dx1=idx1-pos1;      dx2=1-dx1
            rho_grid(idx1(1),idx1(2),idx1(3))=rho_grid(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)
            rho_grid(idx2(1),idx1(2),idx1(3))=rho_grid(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)
            rho_grid(idx1(1),idx2(2),idx1(3))=rho_grid(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)
            rho_grid(idx1(1),idx1(2),idx2(3))=rho_grid(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)
            rho_grid(idx1(1),idx2(2),idx2(3))=rho_grid(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)
            rho_grid(idx2(1),idx1(2),idx2(3))=rho_grid(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)
            rho_grid(idx2(1),idx2(2),idx1(3))=rho_grid(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)
            rho_grid(idx2(1),idx2(2),idx2(3))=rho_grid(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)
            ip=ll(ip) ! find next particle in the chain
          enddo
        enddo
        enddo
        enddo
        !$omp endparalleldo
        enddo
        rho8 = sum(rho_grid(1:nfof,1:nfof,1:nfof))/(nfof**3)
        rho_grid = rho_grid/rho8-1
        if (head) print*,'min',minval(rho_grid(1:nfof,1:nfof,1:nfof)),'max',maxval(rho_grid(1:nfof,1:nfof,1:nfof)),'mean',sum(rho_grid(1:nfof,1:nfof,1:nfof)*1d0)/nfof/nfof/nfof
        if (head) print*,'Write delta_fof into',output_name('delta_fof_'//trim(adjustl(str_refine)))
        open(11,file=output_name('delta_fof_'//trim(adjustl(str_refine))),status='replace',access='stream')
        write(11) rho_grid(1:nfof,1:nfof,1:nfof)
        close(11)
        deallocate(rho_grid)
      endif
    
      ! loop over fof cells
      do l=0,n_layer-1
      !$omp paralleldo  PRIVATE(iq3,iq2,iq1,ip,jp,rsq,i_neighbor)
      do iq3=1+l,nfof,n_layer
      do iq2=1,nfof
      do iq1=1,nfof
        ip=hoc(iq1,iq2,iq3)
        do while (ip/=0)
          jp=ll(ip)
          do while (jp/=0)
            rsq=sum((xv(1:3,ip)-xv(1:3,jp))**2)
            if (rsq<=rp2) call merge_chain(ip,jp)
            jp=ll(jp)
          enddo
          do i_neighbor=1,13
            ! jq=modulo([iq1,iq2,iq3]+ijk(:,i_neighbor)-1,nfof)+1
            jq = [iq1,iq2,iq3]+ijk(:,i_neighbor)
            if (maxval(jq)>nfof1 .or. minval(jq)<1) cycle

            jp=hoc(jq(1),jq(2),jq(3))
            do while (jp/=0)
              rsq=sum((xv(1:3,ip)-xv(1:3,jp))**2)
              if (rsq<=rp2) call merge_chain(ip,jp)
              jp=ll(jp)
            enddo
          enddo
          ip=ll(ip) ! find next particle in the chain
        enddo
      enddo
      enddo
      enddo
      !$omp endparalleldo
      enddo
      deallocate(hoc,ll,ecgp)

      allocate(np_halo_all(np_max/np_halo_min),xv_mean(3,np_max/np_halo_min),np_halo(np_max/np_halo_min)[*])
      np_iso = 0; np_head = 0; np_halo=0; np_halo_all=0
      do i=1,np_max
        if (hcgp(i)==i) then
          np_iso=np_iso+1
        elseif (hcgp(i)/=0) then
          ip = hcgp(i); jp = ip; np = 0; dxv = 0
          do while (jp /= 0)
            np = np + 1
            dxv = dxv + modulo(xv(:,jp) - xv(:,ip) + L_bLb2, L_bLb*1d0)-L_bLb2
            jp = llgp(jp)
          enddo
          dxv = modulo(dxv/np + xv(:,ip), L_bLb*1d0)
          ! if (np > np_halo_min) then
          if (np > np_halo_min .and. minval(dxv(1:3)) > L_b .and. maxval(dxv(1:3)) < L_bL) then
            np_head=np_head+1
            np_halo_all(np_head) = hcgp(i)
            xv_mean(:,np_head) =  dxv!(dxv + shift_xv)*box/nn/nc
            np_halo(np_head) = np
          endif
        endif
      enddo

      nh=np_head

      halo_header%nhalo_tot=0
      halo_header%nhalo=nh; halo_header%ninfo=ninfo; halo_header%linking_parameter=b_link
      if (head) print*, output_name(trim('halo_'//trim(adjustl(str_refine))))
      open(21,file=output_name(trim('halo_'//trim(adjustl(str_refine)))),status='replace',access='stream')
      write(21) halo_header

      allocate(hcat(nh)[*])
      hcat%hmass = np_halo(1:nh)
      do i=1,3
        hcat%xv(i)=xv_mean(i,1:nh)
      enddo
      write(21) hcat 

      if (head) then
        do i=2,nn**3
          nh=nh+nh[i]
        enddo
      endif
      sync all; nh=nh[1]; halo_header%nhalo_tot=nh; 
      rewind(21); write(21) halo_header; close(21)
      sync all
      if (head) print*,'halo_header',halo_header,maxval(hcat%hmass)
      deallocate(hcgp,llgp,xv_mean,np_halo_all)

      ! nh=np_head
      ! if (head) then
      !   do i=2,nn**3
      !     np_head = nh[i]
      !     print*,i,np_head
      !   enddo
      ! endif
      ! sync all

      ! nh=sum(hcat%hmass)
      ! if (head) then
      !   np_head =  nh[1]
      !   do i=2,nn**3
      !     np_head = np_head + nh[i]
      !   enddo
      !   print*,np_head
      ! endif


      nh=np_head
      sync all
      if (head) then
        allocate(ll((ng_global)**3/np_halo_min))
        ll(1:nh) = hcat%hmass
        do i=2,nn**3
          np_head = nh[i]
          print*,i,np_head
          ll(nh+1:nh+np_head) = np_halo(1:np_head)[i]
          nh=nh+np_head
        enddo
        if (head) print*, output_name(trim('halo_mass_'//trim(adjustl(str_refine))))
        open(21,file=output_name(trim('halo_mass_'//trim(adjustl(str_refine)))),status='replace',access='stream')
        write(21) nh
        write(21) ll(1:nh)
        close(21)
        deallocate(ll)
      endif
      deallocate(np_halo,hcat)
      call toc(11)
      if (head) then
         print*,nfof,'  real time =',tcat(11,istep),'secs'; print*,''
      endif
    enddo
    deallocate(xv)
  enddo ! cur_checkpoint



  contains

  subroutine merge_chain(ii,jj)
    ! llgp is a linked list: llgp(ip) means ip->llgp
    ! ip1->ip2->...->ipn
    ! ip1 is the head of chain (hoc)
    ! ipn is the end of chain (eoc)
    integer(4) ii,jj,ihead,jhead,iend,jend,ipart
    jend=merge(jj,ecgp(jj),ecgp(jj)==0)
    iend=merge(ii,ecgp(ii),ecgp(ii)==0)
    if (iend==jend) return ! same chain
    ihead=max(hcgp(ii),hcgp(iend))
    jhead=max(hcgp(jj),hcgp(jend))
    ipart=jhead ! change eoc of the chain-j
    do while (ipart/=0)
      ecgp(ipart)=iend ! set chain-j's eoc to be iend
      ipart=llgp(ipart) ! next particle
    enddo
    llgp(jend)=ihead ! link j group to i group
    ecgp(ii)=iend
    hcgp(iend)=jhead ! change hoc
    hcgp(jend)=0 ! set jend as a member
  endsubroutine
  
end
