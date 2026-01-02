!#define gadget

program CUBE_FoF
  use variables
  implicit none
  integer(8),parameter:: n_refine = 4
  integer(8),parameter:: fof_buffer = ceiling(30.0/box*nc*nn)*n_refine
  integer(8),parameter:: nfof = nc*n_refine+2*fof_buffer
  integer(8),parameter:: nfof2 = nfof/2
  integer i,j,k,l,cur_checkpoint,np,n1,n2,idx(3),nh[*],im,idx1(3),idx2(3)
  integer iq1,iq2,iq3,i_neighbor,jq(3),ijk_neighbor(3,3),neighbor_b(3,3)
  integer(4) nlast,ip,jp,np_iso,np_head,np_max,np_neighbors(3,3,3)[nn,nn,*],np_need(3,3,3),max_nei[nn,nn,*]
  integer(4),allocatable :: np_halo_all(:),np_halo(:)
  integer(4),allocatable :: hoc(:,:,:),ll(:),llgp(:),hcgp(:),ecgp(:),iph_halo_all(:),iph_halo(:)
  real rp2,rsq,dxv(3),pos1(3),dx1(3),dx2(3)
  real(8) rho8
  real,allocatable :: xv(:,:),xv_mean(:,:),xp_neighbors(:,:,:,:,:)[:,:,:],rho_grid(:,:,:)[:,:,:]
  type(type_halo_catalog_header) halo_header
  type(type_halo_catalog_array),allocatable :: hcat(:)[:]
  
  ! ! 检查循环上下限
  ! neighbor_b = reshape([int(0,kind=8),fof_buffer/n_refine,n_refine+fof_buffer  ,int(0,kind=8),nc,int(0,kind=8) ,nc-fof_buffer/n_refine,nc,-fof_buffer ],[3,3])
  ! do iq1 = 1, 3
  !   ijk_neighbor(:, 1) = neighbor_b(:, iq1)
  !   do iq2 = 1, 3
  !   ijk_neighbor(:, 2) = neighbor_b(:, iq2)
  !   do iq3 = 1, 3
  !     np_neighbors(iq1,iq2,iq3) = 0
  !     if (iq1 == 2 .and. iq2 == 2 .and. iq3 == 2) cycle
  !     ijk_neighbor(:, 3) = neighbor_b(:, iq3)
  !     ! 循环遍历各个维度
  !     do itz = floor((ijk_neighbor(1, 3))*1d0/nt)+1, floor((ijk_neighbor(2, 3)-1)*1d0/nt)+1
  !     do ity = floor((ijk_neighbor(1, 2))*1d0/nt)+1, floor((ijk_neighbor(2, 2)-1)*1d0/nt)+1
  !     do itx = floor((ijk_neighbor(1, 1))*1d0/nt)+1, floor((ijk_neighbor(2, 1)-1)*1d0/nt)+1
  !       do k = merge(1, mod(ijk_neighbor(1, 3)+1, nt), itz*nt-nt >= ijk_neighbor(1, 3)),merge(nt, mod(ijk_neighbor(2, 3)+1, nt)-1, itz*nt <= ijk_neighbor(2, 3))
  !       do j = merge(1, mod(ijk_neighbor(1, 2)+1, nt), ity*nt-nt >= ijk_neighbor(1, 2)),merge(nt, mod(ijk_neighbor(2, 2)+1, nt)-1, ity*nt <= ijk_neighbor(2, 2))
  !       do i = merge(1, mod(ijk_neighbor(1, 1)+1, nt), itx*nt-nt >= ijk_neighbor(1, 1)),merge(nt, mod(ijk_neighbor(2, 1)+1, nt)-1, itx*nt <= ijk_neighbor(2, 1))

  !       if (   itz*nt - nt + k - 1 < ijk_neighbor(1, 3) .or. itz*nt - nt + k > ijk_neighbor(2, 3)&
  !         .or. ity*nt - nt + j - 1 < ijk_neighbor(1, 2) .or. ity*nt - nt + j > ijk_neighbor(2, 2)&
  !         .or. itx*nt - nt + i - 1 < ijk_neighbor(1, 1) .or. itx*nt - nt + i > ijk_neighbor(2, 1)) then

  !         ! 计算边界范围
  !         print*,''
  !         print*,'',itx ,ity,itz,i,j,k
  !         print*,'-----------------------------',itx*nt - nt + i,ity*nt - nt + j,itz*nt - nt + k
  !         write(*, '(A, I1, A, 2I6)') 'Level 3 - iq3:', iq3, ' | ijk_neighbor:', int(ijk_neighbor(:2, 3))  
  !         write(*, '(A, I1, A, 2I6)') 'Level 2 - iq2:', iq2, ' | ijk_neighbor:', int(ijk_neighbor(:2, 2))
  !         write(*, '(A, I1, A, 2I6)') 'Level 1 - iq1:', iq1, ' | ijk_neighbor:', int(ijk_neighbor(:2, 1)) 
  !         write(*, '(A, 2I8)') '      Z-range:', int(floor((ijk_neighbor(1, 3)-1)*1d0/nt)+1)*nt-nt+1, &
  !                                                 int(floor((ijk_neighbor(2, 3)-1)*1d0/nt)+1)*nt
  !         write(*, '(A, 2I8)') '      Y-range:', int(floor((ijk_neighbor(1, 2)-1)*1d0/nt)+1)*nt-nt+1, &
  !                                                 int(floor((ijk_neighbor(2, 2)-1)*1d0/nt)+1)*nt
  !         write(*, '(A, 2I8)') '      X-range:', int(floor((ijk_neighbor(1, 1)-1)*1d0/nt)+1)*nt-nt+1, &
  !                                                 int(floor((ijk_neighbor(2, 1)-1)*1d0/nt)+1)*nt
  !         write(*, '(A, I8)') '        Z-tile:', itz
  !         write(*, '(A, I8)') '        Y-tile:', ity
  !         write(*, '(A, I8)') '        X-tile:', itx
  !         ! 计算每个维度的实际索引范围
  !         write(*, '(A, 2I8)') '          Z-indices:', &
  !             int(itz*nt - nt + merge(1, mod(ijk_neighbor(1, 3)+1, nt)-1, itz*nt-nt >= ijk_neighbor(1, 3))), &
  !             int(itz*nt - nt + merge(nt, mod(ijk_neighbor(2, 3)+1, nt)-1, itz*nt <= ijk_neighbor(2, 3)))
  !         write(*, '(A, 2I8)') '          Y-indices:', &
  !             int(ity*nt - nt + merge(1, mod(ijk_neighbor(1, 2)+1, nt)-1, ity*nt-nt >= ijk_neighbor(1, 2))), &
  !             int(ity*nt - nt + merge(nt, mod(ijk_neighbor(2, 2)+1, nt)-1, ity*nt <= ijk_neighbor(2, 2)))
  !         write(*, '(A, 2I8)') '          X-indices:', &
  !             int(itx*nt - nt + merge(1, mod(ijk_neighbor(1, 1)+1, nt)-1, itx*nt-nt >= ijk_neighbor(1, 1))), &
  !             int(itx*nt - nt + merge(nt, mod(ijk_neighbor(2, 1)+1, nt)-1, itx*nt <= ijk_neighbor(2, 1)))
  !       endif
  !       enddo
  !       enddo
  !       enddo
  !     enddo
  !     enddo
  !     enddo
  !   enddo
  !   enddo
  ! enddo
  ! print* ,'------------------------------------------------------'
  ! stop
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call geometry

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

  if (head) print*,'  initialize FoF cell neighbors',n_checkpoint
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

  n1=1; n2=nfof
  rp2=b_link**2
  if (head) print*,nfof,fof_buffer

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

    ! first demension [upper_boundary, lower_boundary, refine_shift] second demension [-1,0,1] image 
    neighbor_b = reshape([int(0,kind=8),fof_buffer/n_refine,nc*n_refine  ,int(0,kind=8),nc,int(0,kind=8) ,nc-fof_buffer/n_refine,nc,-nc*n_refine ],[3,3])
    ! if (head) print *,neighbor_b
    ! stop

  
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
    !   endif
    !   sync all
    ! enddo
    ! if (head) print*,sum(rhoc),sum(np_neighbors),(sum(rhoc)+sum(np_neighbors))*1d0/nfof/nfof/nfof,(sum(rhoc)+sum(np_neighbors))-nfof**3
    ! stop
    max_nei = maxval(np_neighbors)

    ! print*,maxval(np_neighbors),maxval(np_neighbors)*3*3*3*3*4*1d0/1024/1024/1024
    ! stop

    ! init space for all particles -检查过image——neighbor编号和np_need
    sync all
    if (head) print*,' init space for all particles'
    ! do im = 1,nn**3
    !   if (image == im ) then
    !     print*,image
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
    !   endif
    !   sync all
    ! enddo
    ! stop
    np_max = np_max + sum(np_need)
    ! do im = 1,nn**3
    !   if (image == im ) then
    !     print*,image
    !     print*,np_max,sum(rhoc),sim%nplocal,sum(np_need)
    !     print*,np_max*1d0/n2/n2/n2,sum(rhoc)*1d0/nc/nc/nc/n_refine/n_refine/n_refine
    !     print*,3*4*np_max*1d0/1024/1024/1024
    !     print*,''
    !   endif
    !   sync all
    ! enddo
    ! stop

    !  initialize particles  -检查过xv范围和sum得到nlast
    print*,image,np_max
    sync all
    if (head) print*,' initialize particles'
    allocate(xv(3,np_max))
    ! do im = 1,1!nn**3
    !   if (image == im ) then
    !     print*,image,sum(rhoc)*1d0/nc/nc/nc
    xv = 0
    nlast=0
    do itz=1,nnt
    do ity=1,nnt
    do itx=1,nnt
      ! dxv = 0; np_max = 0
      do k=1,nt
      do j=1,nt
      do i=1,nt
        np=rhoc(i,j,k,itx,ity,itz)
        ! print*,i,j,k,nlast,sum(rhoc(:,:,:,:,:,:itz-1))            &
        !             + sum(rhoc(:,:,:,:,:ity-1,itz))         &
        !             + sum(rhoc(:,:,:,:itx-1,ity,itz))       &
        !             + sum(rhoc(:,:,:k-1,itx,ity,itz))       &
        !             + sum(rhoc(:,:j-1,k,itx,ity,itz))    &
        !             + sum(rhoc(:i-1,j,k,itx,ity,itz))
        ! if (nlast /= sum(rhoc(:,:,:,:,:,:itz-1))            &
        !             + sum(rhoc(:,:,:,:,:ity-1,itz))         &
        !             + sum(rhoc(:,:,:,:itx-1,ity,itz))       &
        !             + sum(rhoc(:,:,:k-1,itx,ity,itz))       &
        !             + sum(rhoc(:,:j-1,k,itx,ity,itz))    &
        !             + sum(rhoc(:i-1,j,k,itx,ity,itz))) then 
        !   print *,image,itz,ity,itx
        !   print*,i,j,k,nlast,sum(rhoc(:,:,:,:,:,:itz-1))            &
        !               + sum(rhoc(:,:,:,:,:ity-1,itz))         &
        !               + sum(rhoc(:,:,:,:itx-1,ity,itz))       &
        !               + sum(rhoc(:,:,:k-1,itx,ity,itz))       &
        !               + sum(rhoc(:,:j-1,k,itx,ity,itz))    &
        !               + sum(rhoc(:i-1,j,k,itx,ity,itz))
        ! endif
        do l=1,np
          ip=nlast+l
          if (ip > np_max  .or. ip < 1 .or.  ip > sim%nplocal) then
            print*, image
            print*, i,j,k,itx,ity,itz
            print*,nlast,l,np,ip
            stop
          endif
#ifdef ZIPX
          xv(:,ip)=(nt*((/itx,ity,itz/)-1)+ ((/i,j,k/)-1) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution)*n_refine+[fof_buffer,fof_buffer,fof_buffer] 
#else 
          xv(:,ip)=xp(:,ip)*n_refine/ratio_cs+[fof_buffer,fof_buffer,fof_buffer] 
#endif
          !  dxv =  modulo(xv(:,nlast) - xv(:,ip)+nfof2,nfof*1d0)-nfof2
          !  dxv =  dxv + xv(:,nlast)
        enddo
        ! np_max =  np_max + np
        ! dxv = modulo(sum(dxv)/np+xv(:,nlast),nfof*1d0)
        ! dxv = sum(dxv)/np+xv(:,nlast)
        ! print*,'  ',i,j,k,dxv
        nlast=nlast+np
      enddo
      enddo
      enddo
      ! dxv = sum(dxv)/np_max
      ! print*,'   ',itx*nt-nt,ity*nt-nt,itz*nt-nt
      ! print*,'   ',dxv
      ! print*,'   ',itx*nt,ity*nt,itz*nt
    enddo
    enddo
    enddo
    !   endif
    !   sync all
    ! enddo
    ! stop
    deallocate(xp)
    ! do im = 1,nn**3
    !   if (image == im ) then
    !     print*,image,nlast,sum(rhoc),sum(xv(1,:))/nlast,sum(xv(2,:))/nlast,sum(xv(3,:))/nlast
    !   endif
    !   sync all
    ! enddo
    ! stop


    ! init particles for neighbors  -检查粒子的中心位置和最大最小值
    if (head) then
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
    if (head) print*,' init particles for neighbors',max_nei,maxval(np_neighbors)

    allocate(xp_neighbors(3,max_nei,3,3,3)[nn,nn,*])
    xp_neighbors = 0
    ! do im = 1,nn**3
    !   if (image == im ) then
    !     print*,image
    do iq1 = 1, 3
    ijk_neighbor(:, 1) = neighbor_b(:, iq1)
    do iq2 = 1, 3
    ijk_neighbor(:, 2) = neighbor_b(:, iq2)
    do iq3 = 1, 3
      jp = 0
      if (iq1 == 2 .and. iq2 == 2 .and. iq3 == 2) cycle
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
      ! print*,iq1-2,iq2-2,iq3-2,np_neighbors(iq1,iq2,iq3),jp
      ! print*,sum(xp_neighbors(1,:jp,iq1,iq2,iq3))/jp,minval(xp_neighbors(1,:jp,iq1,iq2,iq3)),maxval(xp_neighbors(1,:jp,iq1,iq2,iq3)),ijk_neighbor(3,1)
      ! print*,sum(xp_neighbors(2,:jp,iq1,iq2,iq3))/jp,minval(xp_neighbors(2,:jp,iq1,iq2,iq3)),maxval(xp_neighbors(2,:jp,iq1,iq2,iq3)),ijk_neighbor(3,2)
      ! print*,sum(xp_neighbors(3,:jp,iq1,iq2,iq3))/jp,minval(xp_neighbors(3,:jp,iq1,iq2,iq3)),maxval(xp_neighbors(3,:jp,iq1,iq2,iq3)),ijk_neighbor(3,3)
    enddo
    enddo
    enddo
    !   endif
    !   sync all
    ! enddo
    ! stop

    ! do im = 1,nn**3
    !   if (image == im ) then
    !     print*,image,maxval(np_neighbors)
    !   endif
    !   sync all
    ! enddo
    ! sync all
    ! print*,'112'
    ! xv(:,1:557) = xp_neighbors(:,:557,3,3,3)[1,1,2]
    ! if (image == 8) xv(:,sum(rhoc)+1:sum(rhoc)+557) = xp_neighbors(:,:557,3,3,3)
    ! sync all
    ! print*,'222'
    ! xv(:,1:557) = xp_neighbors(:,:557,3,3,3)[2,2,2]
    ! stop

    ! sync all
    ! ! sync images(8)  ! 等同于sync all，但更明确
    ! ! if (image == 1) then
    ! !     sync images([2,2,2])  ! 先与目标image同步
    ! !     xv(:,sum(rhoc)+1:sum(rhoc)+557) = xp_neighbors(:,:557,3,3,3)[2,2,2]
    ! ! endif
    ! ! sync images(1)
    ! ! stop
    ! ! if (head) print*,' 1'
    ! ! if (image == 8) xv(:,sum(rhoc)+1:sum(rhoc)+557) = xp_neighbors(:,:557,3,3,3)[2,2,2]
    ! ! sync all
    ! ! if (head) print*,' 2'
    ! ! if (image == 8) xv(:,sum(rhoc)+1:sum(rhoc)+557) = xp_neighbors(:,:557,3,3,3)
    ! sync all
    ! xp_neighbors = image
    ! xv(:,557*image-557+1:557*image)[1,1,1] = xp_neighbors(:,:557,3,3,3)
    ! print *,image
    ! sync all
    ! if (head) then
    !   do i = 1,nn**3
    !   print*,xv(1,557*i)
    !   enddo
    ! endif
    ! sync all 
    ! stop
    ! if (image == 8) then
    !    xv(:,sum(rhoc)+1:sum(rhoc)+557)[1,1,1] = xp_neighbors(:,:557,3,3,3)
    ! end if
    ! sync all
    ! if (image == 1) xv(:,sum(rhoc)+1:sum(rhoc)+557) = xp_neighbors(:,:557,3,3,3)[2,2,2]
    ! ! if (image == 1) print *,image,xp_neighbors(1,1,3,3,3)[2,2,2]
    ! ! if (head) print*,' 3'
    ! ! critical 
    ! !   if (image == 1) xv(:,sum(rhoc)+1:sum(rhoc)+557) = xp_neighbors(:,:557,3,3,3)[2,2,2]
    ! ! end critical
    ! sync all
    ! stop


    ! get  particles form neighbors  -检查粒子的中心位置和最大最小值
    sync all
    if (head) print*,' get particles form neighbors'
    np_max = sum(rhoc)
    do iq1 = 1,3
      i = modulo(icx-3+iq1,nn)+1
    do iq2 = 1,3
      j = modulo(icy-3+iq2,nn)+1
    do iq3 = 1,3
      k = modulo(icz-3+iq3,nn)+1
      if (iq1 == 2 .and. iq2 == 2 .and. iq3 == 2) cycle
      xv(:,np_max+1:np_max+np_need(iq1,iq2,iq3)) = xp_neighbors(:,:np_need(iq1,iq2,iq3),4-iq1,4-iq2,4-iq3)[i,j,k]
      np_max = np_max +  np_need(iq1,iq2,iq3)
    enddo
    enddo
    enddo
    ! stop


    ! ! get  particles form neighbors  -检查粒子的中心位置和最大最小值
    ! sync all
    ! if (head) print*,' get particles form neighbors'
    ! ! do im = 1,nn**3
    !   if (image == 1 ) then
    !     print*,image
    ! np_max = sum(rhoc)
    ! do iq1 = 1,3
    !   i = modulo(icx-3+iq1,nn)+1
    ! do iq2 = 1,3
    !   j = modulo(icy-3+iq2,nn)+1
    ! do iq3 = 1,3
    !   k = modulo(icz-3+iq3,nn)+1
    !   if (iq1 == 2 .and. iq2 == 2 .and. iq3 == 2) cycle
    !   write(*,'(i0xxx,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0)') image, i, j, k, np_need(iq1,iq2,iq3), max_nei[i,j,k], np_max+np_need(iq1,iq2,iq3),sum(rhoc)+sum(np_need)
    !   ! if ((np_max+np_need(iq1,iq2,iq3) >  sum(rhoc) + sum(np_need)) .or. (np_need(iq1,iq2,iq3) > max_nei[i,j,k]) .or. maxval([i,j,k]) > nn .or. minval([i,j,k]) <1) then
    !   !   print*,image
    !   !   print*,iq1-2,iq2-2,iq3-2
    !   !   print*,i,j,k,np_need(iq1,iq2,iq3)
    !   !   print*,np_max,np_max+np_need(iq1,iq2,iq3)
    !   !   print*,max_nei[i,j,k]
    !   !   print*,'get  particle error'
    !   !   stop
    !   ! endif
    !   print*,4-iq1,4-iq2,4-iq3
    !   xv(:,np_max+1:np_max+np_need(iq1,iq2,iq3)) = xp_neighbors(:,:np_need(iq1,iq2,iq3),4-iq1,4-iq2,4-iq3)[i,j,k]
    !   np_max = np_max +  np_need(iq1,iq2,iq3)
    !   ! print*,iq1-2,iq2-2,iq3-2
    !   ! print*,i,j,k,np_need(iq1,iq2,iq3)
    !   ! print*,sum(xv(1,np_max+1:np_max+np_need(iq1,iq2,iq3)))/np_need(iq1,iq2,iq3),minval(xv(1,np_max+1:np_max+np_need(iq1,iq2,iq3))),maxval(xv(1,np_max+1:np_max+np_need(iq1,iq2,iq3)))
    !   ! print*,sum(xv(2,np_max+1:np_max+np_need(iq1,iq2,iq3)))/np_need(iq1,iq2,iq3),minval(xv(2,np_max+1:np_max+np_need(iq1,iq2,iq3))),maxval(xv(2,np_max+1:np_max+np_need(iq1,iq2,iq3)))
    !   ! print*,sum(xv(3,np_max+1:np_max+np_need(iq1,iq2,iq3)))/np_need(iq1,iq2,iq3),minval(xv(3,np_max+1:np_max+np_need(iq1,iq2,iq3))),maxval(xv(3,np_max+1:np_max+np_need(iq1,iq2,iq3)))
    ! enddo
    ! enddo
    ! enddo
    !   endif
    ! !   sync all
    ! ! enddo
    ! stop
    ! sync all
    ! ! if (head) print*, 'particle initialized, np_max = ', np_max
    

    ! np_max = sum(rhoc)
    ! do im = 1,nn**3
    !   if (image == im ) then
    !     print*,image,sum(xv(1,:np_max))/np_max,minval(xv(1,:np_max)),maxval(xv(1,:np_max))
    !     print*,image,sum(xv(2,:np_max))/np_max,minval(xv(2,:np_max)),maxval(xv(2,:np_max))
    !     print*,image,sum(xv(3,:np_max))/np_max,minval(xv(3,:np_max)),maxval(xv(3,:np_max))
    !   endif
    !   sync all
    ! enddo
    ! sync all
    ! if  (head) print*, ' '

    ! np_max = sum(rhoc) + sum(np_need)
    ! deallocate(xp_neighbors,rhoc)
    ! do im = 1,nn**3
    !   if (image == im ) then
    !     ! print*,image,'min',minval(rho_grid(n1:n2,n1:n2,n1:n2)),'max',maxval(rho_grid(n1:n2,n1:n2,n1:n2)),'mean',sum(rho_grid(n1:n2,n1:n2,n1:n2)*1d0)/nw/nw/nw
    !     print*,image,sum(xv(1,:np_max))/np_max,minval(xv(1,:np_max)),maxval(xv(1,:np_max))
    !     print*,image,sum(xv(2,:np_max))/np_max,minval(xv(2,:np_max)),maxval(xv(2,:np_max))
    !     print*,image,sum(xv(3,:np_max))/np_max,minval(xv(3,:np_max)),maxval(xv(3,:np_max))
    !   endif
    !   sync all
    ! enddo
    ! stop
    ! deallocate(xp_neighbors,rhoc)
    ! if (head) print*, 'particle initialized, np_max = ', np_max



    ! create hoc ll
    sync all
    if (head) print*,' create hoc ll'
    ! do im = 1,nn**3
    !   if (image == im ) then
    !     print*,image
    np_max = sum(rhoc) + sum(np_need)
    allocate(hoc(n1:n2,n1:n2,n1:n2),ll(np_max),llgp(np_max),hcgp(np_max),ecgp(np_max))

    allocate(rho_grid(n1-1:n2+1,n1-1:n2+1,n1-1:n2+1)[nn,nn,*])
    rho_grid=0
    hoc=0; ll=0
    do ip=1,np_max
      idx=floor(mod(xv(1:3,ip),nfof*1d0))+1 ! index of the grid
      if (minval(idx)<n1 .or. maxval(idx)>n2) then
        print*, 'idx out of range'
        print*, idx
        stop
      endif
      rho_grid(idx(1),idx(2),idx(3)) = rho_grid(idx(1),idx(2),idx(3))+1
      ll(ip)=hoc(idx(1),idx(2),idx(3)) ! linked list 
      hoc(idx(1),idx(2),idx(3))=ip ! head of chain
      hcgp(ip)=ip ! initialize hcgp(ip)=ip for isolated particles
    enddo
    ! do im = 1,nn**3
    !   if (image == im ) then
    !     print*,image,'min',minval(rho_grid(n1:n2,n1:n2,n1:n2)),'max',maxval(rho_grid(n1:n2,n1:n2,n1:n2)),'mean',sum(rho_grid(n1:n2,n1:n2,n1:n2)*1d0)/nw/nw/nw
    !   endif
    !   sync all
    ! enddo
    ! stop
    llgp=0; ecgp=0; ! initialize group link list
    if (head) print*, 'hoc ll created'


    ! ! allocate(rho_grid(n1-1:n2+1,n1-1:n2+1,n1-1:n2+1)[nn,nn,*])
    ! rho_grid=0
    ! do iq3=n1,n2
    ! do iq2=n1,n2
    ! do iq1=n1,n2
    !   ip=hoc(iq1,iq2,iq3)
    !   do while (ip/=0)
    !     pos1 = xv(1:3,ip)
    !     if ( pos1(1) == nfof ) pos1(1)= 0
    !     if ( pos1(2) == nfof ) pos1(2)= 0
    !     if ( pos1(3) == nfof ) pos1(3)= 0
    !     idx1=floor(pos1)+1; idx2=idx1+1
    !     dx1=idx1-pos1;      dx2=1-dx1
    !     rho_grid(idx1(1),idx1(2),idx1(3))=rho_grid(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)
    !     rho_grid(idx2(1),idx1(2),idx1(3))=rho_grid(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)
    !     rho_grid(idx1(1),idx2(2),idx1(3))=rho_grid(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)
    !     rho_grid(idx1(1),idx1(2),idx2(3))=rho_grid(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)
    !     rho_grid(idx1(1),idx2(2),idx2(3))=rho_grid(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)
    !     rho_grid(idx2(1),idx1(2),idx2(3))=rho_grid(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)
    !     rho_grid(idx2(1),idx2(2),idx1(3))=rho_grid(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)
    !     rho_grid(idx2(1),idx2(2),idx2(3))=rho_grid(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)
    !     ip=ll(ip) ! find next particle in the chain
    !   enddo
    ! enddo
    ! enddo
    ! enddo
    ! rho8 = sum(rho_grid(n1:n2,n1:n2,n1:n2))/(nfof**3)
    ! do im = 1,nn**3
    !   if (image == im ) then
    !     print*,image,rho8
    !     print*,'min',minval(rho_grid(n1:n2,n1:n2,n1:n2)),'max',maxval(rho_grid(n1:n2,n1:n2,n1:n2)),'mean',sum(rho_grid(n1:n2,n1:n2,n1:n2)*1d0)/n2/n2/n2
    !   endif
    !   sync all
    ! enddo
    ! rho_grid = rho_grid/rho8-1
    ! if (head) print*,'Write delta_fof into',output_name('delta_fof')
    ! open(11,file=output_name('delta_fof'),status='replace',access='stream')
    ! write(11) rho_grid(n1:n2,n1:n2,n1:n2)
    ! close(11)
    ! do im = 1,nn**3
    !   if (image == im ) then
    !     print*,image,rho8
    !     print*,'min',minval(rho_grid(n1:n2,n1:n2,n1:n2)),'max',maxval(rho_grid(n1:n2,n1:n2,n1:n2)),'mean',sum(rho_grid(n1:n2,n1:n2,n1:n2)*1d0)/n2/n2/n2
    !   endif
    !   sync all
    ! enddo
    ! stop
  

    ! loop over fof cells
    do iq3=n1,n2
    do iq2=n1,n2
    do iq1=n1,n2
      ip=hoc(iq1,iq2,iq3)
      do while (ip/=0)
        jp=ll(ip)
        do while (jp/=0)
          rsq=sum((xv(1:3,ip)-xv(1:3,jp))**2)
          if (rsq<=rp2) call merge_chain(ip,jp)
          jp=ll(jp)
        enddo
        do i_neighbor=1,13
          jq=[iq1,iq2,iq3]+ijk(:,i_neighbor)
          if (minval(jq)<n1 .or. maxval(jq)>n2) cycle
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
    deallocate(hoc,ll,ecgp)



    allocate(np_halo_all(np_max/np_halo_min),xv_mean(3,np_max/np_halo_min),np_halo(np_max/np_halo_min))
    np_iso = 0; np_head = 0; np_halo=0; np_halo_all=0
    do i=1,np_max
      if (hcgp(i)==i) then
        np_iso=np_iso+1
      elseif (hcgp(i)/=0) then
        ip = hcgp(i); jp = ip; np = 0; dxv = 0
        do while (jp /= 0)
          np = np + 1
          dxv = dxv + modulo(xv(:,jp) - xv(:,ip)+nfof2,nfof*1d0)-nfof2
          jp = llgp(jp)
        enddo
        dxv = modulo(sum(dxv)/np+xv(:,ip),nfof*1d0)
        if (np > np_halo_min .and. minval(dxv(1:3)) > fof_buffer .and. maxval(dxv(1:3)) < nfof-fof_buffer) then
          np_head=np_head+1
          np_halo_all(np_head) = hcgp(i)
          xv_mean(:,np_head) =  dxv 
          np_halo(np_head) = np
        endif
      endif
    enddo

    nh=np_head

    halo_header%nhalo_tot=0
    halo_header%nhalo=nh; halo_header%ninfo=ninfo; halo_header%linking_parameter=b_link
    if (head) print*, output_name('halo')
    open(21,file=output_name('halo'),status='replace',access='stream')
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
    if (head) print*,'halo_header',halo_header
    deallocate(hcgp,llgp,xv,xv_mean)

    nh=np_head
    if (head) then
      allocate(ll((ng_global)**3/np_halo_min))
      ll(1:nh) = hcat%hmass
      do i=2,nn**3
        np_head = nh[i]
        ll(nh+1:nh+np_head) = hcat[i]%hmass
        nh=nh+np_head
      enddo
      if (head) print*, output_name('halo_mass')
      open(21,file=output_name('halo_mass'),status='replace',access='stream')
      write(21) nh
      write(21) ll(1:nh)
      close(21)
    endif
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
