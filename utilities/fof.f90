!#define gadget

program CUBE_FoF
  use variables
  implicit none
  integer(8),parameter:: n_refine = 4
  integer(8),parameter:: fof_buffer = ceiling(10/box*nc*nn)*n_refine
  integer(8),parameter:: nfof = nc*n_refine+2*fof_buffer
  integer(8),parameter:: nfof2 = nfof/2
  integer i,j,k,l,cur_checkpoint,np,n1,n2,idx(3),nh[*]
  integer iq1,iq2,iq3,i_neighbor,jq(3),ijk_neighbor(3,3),neighbor_b(3,3)
  integer(4) nlast,ip,jp,np_iso,np_head,np_max,np_neighbors(3,3,3)[nn,nn,*],np_need(3,3,3)
  integer(4),allocatable :: np_halo_all(:),np_halo(:)
  integer(4),allocatable :: hoc(:,:,:),ll(:),llgp(:),hcgp(:),ecgp(:),iph_halo_all(:),iph_halo(:)
  real rp2,rsq,dxv(3)
  real,allocatable :: xv(:,:),xv_mean(:,:),xp_neighbors(:,:,:,:,:)[:,:,:]
  type(type_halo_catalog_header) halo_header
  type(type_halo_catalog_array),allocatable :: hcat(:)[:]
  
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

  if (head) print*,'  initialize FoF cell neighbors'
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

  n1=1; n2=nc*n_refine+2*fof_buffer
  rp2=b_link**2

  do cur_checkpoint=100,100
    sim%cur_checkpoint=cur_checkpoint
    if (head) print*, ''
    if (head) print*, ''
    if (head) print*, 'FoF at redshift ',z2str(z_checkpoint(cur_checkpoint))
    if (head) print*, '  read checkpoint header',output_name('info')

    !load particle
    open(11,file=output_name('info'),access='stream'); read(11) sim; close(11)
    allocate(rhoc(nt,nt,nt,nnt,nnt,nnt)[nn,nn,*],xp(3,sim%nplocal)[nn,nn,*])
    open(11,file=output_name('xp'),access='stream'); read(11) xp; close(11)
    open(11,file=output_name('np'),access='stream'); read(11) rhoc; close(11)

    np_neighbors = 0
    ! first demension [upper_boundary, lower_boundary, refine_shift] second demension [-1,0,1] image 
    neighbor_b = reshape([int(0,kind=8),fof_buffer/n_refine,n_refine+fof_buffer  ,int(0,kind=8),nc,int(0,kind=8) ,nc-fof_buffer/n_refine,nc,-fof_buffer ],[3,3])

    ! count  neighbor particle 
    do iq1 = 1,3
    ijk_neighbor(:,1) = neighbor_b(:,iq1)
    do iq2 = 1,3
    ijk_neighbor(:,2) = neighbor_b(:,iq2)
    do iq3 = 1,3
      if (iq1 ==0 .and. iq2 ==0 .and. iq3==0)  cycle
      ijk_neighbor(:,3) = neighbor_b(:,iq3)
      ! do itz=1,nnt
      ! jq(3) = itz * nt
      ! if (jq(3) > ijk_neighbor(2,3) + nt .and.  jq(3) < ijk_neighbor(1,3))  cycle
      ! do ity=1,nnt
      ! jq(2) = ity * nt
      ! if (jq(2) > ijk_neighbor(2,2) + nt .and.  jq(2) < ijk_neighbor(1,2))  cycle
      ! do itx=1,nnt
      ! jq(1) = itx * nt
      ! if (jq(1) > ijk_neighbor(2,1) + nt .and.  jq(1) < ijk_neighbor(1,1))  cycle
      !   do k=1,nt
      !   if (jq(3) + k  > ijk_neighbor(2,3) .and.  jq(3) + k < ijk_neighbor(1,3))  cycle
      !   do j=1,nt
      !   if (jq(2) + j  > ijk_neighbor(2,2) .and.  jq(2) + j < ijk_neighbor(1,2))  cycle
      !   do i=1,nt
      !   if (jq(1) + i  > ijk_neighbor(2,1) .and.  jq(1) + i < ijk_neighbor(1,1))  cycle
      do itz = floor(ijk_neighbor(1,3)*1d0/nt),floor(ijk_neighbor(2,3)*1d0/nt)+1
      do ity = floor(ijk_neighbor(1,2)*1d0/nt),floor(ijk_neighbor(2,2)*1d0/nt)+1
      do itx = floor(ijk_neighbor(1,1)*1d0/nt),floor(ijk_neighbor(2,1)*1d0/nt)+1
        do k = merge(1,mod(ijk_neighbor(1,3)+1,nt),itz*nt<ijk_neighbor(1,3)),merge(nt,mod(ijk_neighbor(2,3)+1,nt),itz*nt+nt>ijk_neighbor(2,3))
        do j = merge(1,mod(ijk_neighbor(1,2)+1,nt),ity*nt<ijk_neighbor(1,2)),merge(nt,mod(ijk_neighbor(2,2)+1,nt),ity*nt+nt>ijk_neighbor(2,2))
        do i = merge(1,mod(ijk_neighbor(1,1)+1,nt),itx*nt<ijk_neighbor(1,1)),merge(nt,mod(ijk_neighbor(2,1)+1,nt),itx*nt+nt>ijk_neighbor(2,1))
          np_neighbors(iq1,iq2,iq3)= np_neighbors(iq1,iq2,iq3) + rhoc(i,j,k,itx,ity,itz)
        enddo
        enddo
        enddo
      enddo
      enddo
      enddo
    enddo
    enddo
    enddo
    sync all

    ! init space for all particles
    np_max = sum(rhoc)
    do iq1 = 1,3
      i = modulo(icx-3+iq1,nn)+1
    do iq2 = 1,3
      j = modulo(icy-3+iq2,nn)+1
    do iq3 = 1,3
      k = modulo(icz-3+iq3,nn)+1
      np_need(iq1,iq2,iq3) = np_neighbors(4-iq1,4-iq2,4-iq3)[i,j,k]
    enddo
    enddo
    enddo
    np_max = np_max + sum(np_need)
    allocate(xv(3,np_max))

    !  initialize particles
    do itz=1,nnt
        do ity=1,nnt
          do itx=1,nnt
              do k=1,nt
                do j=1,nt
                    do i=1,nt
                      np=rhoc(i,j,k,itx,ity,itz)
                      do l=1,np
                          ip=nlast+l
#ifdef ZIPX
                          xv(:,ip)=nt*((/itx,ity,itz/)-1)+ ((/i,j,k/)-1) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution*n_refine+[fof_buffer,fof_buffer,fof_buffer]
#else 
                          xv(:,ip)=xp(:,ip)*n_refine/ratio_cs+[fof_buffer,fof_buffer,fof_buffer]
#endif
                      enddo
                      nlast=nlast+np
                    enddo
                enddo
              enddo
          enddo
        enddo
    enddo
    deallocate(xp)

    ! init particles for neighbors
    allocate(xp_neighbors(3,maxval(np_neighbors),3,3,3)[nn,nn,*])
    do iq1 = 1,3
    ijk_neighbor(:,1) = neighbor_b(:,iq1)
    do iq2 = 1,3
    ijk_neighbor(:,2) = neighbor_b(:,iq2)
    do iq3 = 1,3
      if (iq1 ==0 .and. iq2 ==0 .and. iq3==0)  cycle
      ijk_neighbor(:,3) = neighbor_b(:,iq3)
      ! do itz=1,nnt
      ! jq(3) = itz * nt
      ! if (jq(3) > ijk_neighbor(2,3) + nt .and.  jq(3) < ijk_neighbor(1,3))  cycle
      ! do ity=1,nnt
      ! jq(2) = ity * nt
      ! if (jq(2) > ijk_neighbor(2,2) + nt .and.  jq(2) < ijk_neighbor(1,2))  cycle
      ! do itx=1,nnt
      ! jq(1) = itx * nt
      ! if (jq(1) > ijk_neighbor(2,1) + nt .and.  jq(1) < ijk_neighbor(1,1))  cycle
      !   do k=1,nt
      !   if (jq(3) + k  > ijk_neighbor(2,3) .and.  jq(3) + k < ijk_neighbor(1,3))  cycle
      !   do j=1,nt
      !   if (jq(2) + j  > ijk_neighbor(2,2) .and.  jq(2) + j < ijk_neighbor(1,2))  cycle
      !   do i=1,nt
      !   if (jq(1) + i  > ijk_neighbor(2,1) .and.  jq(1) + i < ijk_neighbor(1,1))  cycle
      do itz = floor(ijk_neighbor(1,3)*1d0/nt),floor(ijk_neighbor(2,3)*1d0/nt)+1
      do ity = floor(ijk_neighbor(1,2)*1d0/nt),floor(ijk_neighbor(2,2)*1d0/nt)+1
      do itx = floor(ijk_neighbor(1,1)*1d0/nt),floor(ijk_neighbor(2,1)*1d0/nt)+1
        do k = merge(1,mod(ijk_neighbor(1,3)+1,nt),itz*nt<ijk_neighbor(1,3)),merge(nt,mod(ijk_neighbor(2,3)+1,nt),itz*nt+nt>ijk_neighbor(2,3))
        do j = merge(1,mod(ijk_neighbor(1,2)+1,nt),ity*nt<ijk_neighbor(1,2)),merge(nt,mod(ijk_neighbor(2,2)+1,nt),ity*nt+nt>ijk_neighbor(2,2))
        do i = merge(1,mod(ijk_neighbor(1,1)+1,nt),itx*nt<ijk_neighbor(1,1)),merge(nt,mod(ijk_neighbor(2,1)+1,nt),itx*nt+nt>ijk_neighbor(2,1))
          nlast = sum(rhoc(:,:,:,:,:,:itz-1))         &
                + sum(rhoc(:,:,:,:,:ity-1,itz))        &
                + sum(rhoc(:,:,:,:itx-1,ity,itz))       &
                + sum(rhoc(:,:,:k-1,itx,ity,itz))       &
                + sum(rhoc(:,:j-1,:k-1,itx,ity,itz))    &
                + sum(rhoc(:i-1,:j-1,:k-1,itx,ity,itz))
          np=rhoc(i,j,k,itx,ity,itz)
          do l=1,np
              ip=nlast+l
              xp_neighbors(:,jp,iq1,iq2,iq3) = xv(:,ip) + ijk_neighbor(3,:)
          enddo
        enddo
        enddo
        enddo
      enddo
      enddo
      enddo
    enddo
    enddo
    enddo
    sync all

    ! get  particles form neighbors
    np_max = sum(rhoc)
    do iq1 = 1,3
      i = modulo(icx-3+iq1,nn)+1
    do iq2 = 1,3
      j = modulo(icy-3+iq2,nn)+1
    do iq3 = 1,3
      k = modulo(icz-3+iq3,nn)+1
      np_max = np_max +  np_need(iq1,iq2,iq3)
      xv(:,np_max+1:np_max+np_need(iq1,iq2,iq3)) = xp_neighbors(:,1:np_need(iq1,iq2,iq3),4-iq1,4-iq2,4-iq3)[i,j,k]
    enddo
    enddo
    enddo
    sync all
    np_max = sum(rhoc) + sum(np_need)
    deallocate(xp_neighbors,rhoc)
    print*, 'particle initialized, np_max = ', np_max
    




    ! create hoc ll
    allocate(hoc(n1:n2,n1:n2,n1:n2),ll(np_max),llgp(np_max),hcgp(np_max),ecgp(np_max))
    hoc=0; ll=0
    do ip=1,np_max
      idx=floor(mod(xv(1:3,ip),nfof*1d0))+1 ! index of the grid
      if (minval(idx)<n1 .or. maxval(idx)>n2) then
        print*, 'idx out of range'
        print*, idx
        stop
      endif
      ll(ip)=hoc(idx(1),idx(2),idx(3)) ! linked list 
      hoc(idx(1),idx(2),idx(3))=ip ! head of chain
      hcgp(ip)=ip ! initialize hcgp(ip)=ip for isolated particles
    enddo
    llgp=0; ecgp=0; ! initialize group link list
    if (head) print*, 'hoc ll created'

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
