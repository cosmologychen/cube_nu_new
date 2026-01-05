!#define gadget

program CUBE_FoF
  use variables
  implicit none
  character(len = 4) str_refine
  ! integer(8),parameter:: n_refine = 5
  integer(8),parameter:: fof_buffer = ceiling(10.0/box*nc*nn)
  integer(8),parameter:: nfof1 = nfof+1
  real(8),parameter:: n_refine = nfof*1d0/(nc+fof_buffer*2+b_link/ratio_cs)
  real(8),parameter:: L_b    = fof_buffer
  real(8),parameter:: L_bL   = L_b  + nc
  real(8),parameter:: L_bLb  = L_bL + fof_buffer
  real(8),parameter:: L_bLb2 = L_bLb/2
  real(8),parameter:: rp2=(b_link/ratio_cs)**2

  integer(8),parameter:: nlayer = ceiling(fof_buffer*n_refine)+1
  integer,parameter:: real_images = 4
  integer,parameter:: layer_image = nn**3/real_images

  integer :: my_id          ! 当前 image 的 ID
  integer :: log_unit       ! 文件单元号
  character(len=64) :: log_filename ! 文件名字符串
  
  integer image_now,i,j,k,l,cur_checkpoint,np,n1,n2,idx(3),nh[*],im,idx1(3),idx2(3),ft(5,nn*3)[*],ftr(nn*3)[*]
  integer(8) iq1,iq2,iq3,i_neighbor,jq(3),ijk_neighbor(3,3),neighbor_b(3,3)
  integer(8) nhalo_all,halo_images(nn**3)[*]
  integer(4) nlast,ip,jp,np_iso,np_head,np_max,np_neighbors(3,3,3)[nn,nn,*],np_need(3,3,3),max_nei[nn,nn,*]
  integer(4),allocatable :: np_halo_all(:),np_halo(:)[:], rhoc_local(:,:,:,:,:,:)
  integer(4),allocatable :: hoc(:,:,:),ll(:),llgp(:),hcgp(:),ecgp(:),iph_halo_all(:),iph_halo(:)!,h2(:,:,:,:),l2(:,:)
  real rsq,pos1(3),dx1(3),dx2(3),shift_xv(3)
  real(8) rho8,dxv(3)
  real,allocatable :: xv(:,:),xv_mean(:,:),xp_neighbors(:,:,:,:,:)[:,:,:],rho_grid(:,:,:)[:,:,:]
  type(type_halo_catalog_header) halo_header
  type(type_halo_catalog_array),allocatable :: hcat(:)[:]
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (real_images*layer_image /= nn**3) then
    print*,real_images,layer_image,nn**3
    stop 'real_images*layer_image /= nn*3'
  endif
  head=(this_image()==1)
  write(str_refine,'(i4)')  nfof
  if(head) print*,nfof,n_refine,4/n_refine

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

  shift_xv  = [(icx-1)*nc*1d0-fof_buffer,(icy-1)*nc*1d0-fof_buffer,(icz-1)*nc*1d0-fof_buffer]
  neighbor_b = reshape([int(0,kind=8),fof_buffer,nc  ,int(0,kind=8),nc,int(0,kind=8) ,nc-fof_buffer,nc,-nc ],[3,3])
  if (head) print*,'linking,L _b,L_bL,L_bLb',b_link/ratio_cs,L_b,L_bL,L_bLb

  do cur_checkpoint=n_checkpoint,n_checkpoint
  ! do cur_checkpoint=4,4
    sim%cur_checkpoint=cur_checkpoint
    if (head) print*, ''
    if (head) print*, ''
    if (head) print*, 'FoF at redshift ',z2str(z_checkpoint(cur_checkpoint))

    nhalo_all = 0
    do image_now = this_image(),nn**3,real_images
      write(log_filename, "('run_output_', I4.4, '.txt')") image_now


      open(newunit=log_unit, file=trim(log_filename), status='replace', action='write')
      write(log_unit, *) '  Fof at image',image_now
      
      ! init image
      call geometry_images

      ! count  neighbor particle 
      call system_clock(ft(1,image_now),ftr(image_now))
      write(log_unit, *)' count  neighbor particle '
      allocate(rhoc_local(nt,nt,nt,nnt,nnt,nnt))
      ! do im = 1,nn**3
      !   if (image == im ) then
      !     print*,image,sum(rhoc_local)*1d0/nc/nc/nc
      np_neighbors = 0
      do iq1 = 1, 3
      do iq2 = 1, 3
      do iq3 = 1, 3
        ijk_neighbor(:, 1) = neighbor_b(:, 4-iq1)
        ijk_neighbor(:, 2) = neighbor_b(:, 4-iq2)
        ijk_neighbor(:, 3) = neighbor_b(:, 4-iq3)

        inx=modulo(icx+1-iq1,nn)+1
        iny=modulo(icy+1-iq2,nn)
        inz=modulo(icz+1-iq3,nn)
        image = inz*nn**2+iny*nn+inx
        ! write(log_unit, *)image,inx,iny,inz

        open(11,file=output_name('np'),access='stream'); read(11) rhoc_local; close(11)
        ! 循环遍历各个维度
        do itz = floor((ijk_neighbor(1, 3))*1d0/nt)+1, floor((ijk_neighbor(2, 3)-1)*1d0/nt)+1
        do ity = floor((ijk_neighbor(1, 2))*1d0/nt)+1, floor((ijk_neighbor(2, 2)-1)*1d0/nt)+1
        do itx = floor((ijk_neighbor(1, 1))*1d0/nt)+1, floor((ijk_neighbor(2, 1)-1)*1d0/nt)+1
          do k = merge(1, mod(ijk_neighbor(1, 3)+1, nt), itz*nt-nt >= ijk_neighbor(1, 3)),merge(nt, mod(ijk_neighbor(2, 3)+1, nt)-1, itz*nt <= ijk_neighbor(2, 3))
          do j = merge(1, mod(ijk_neighbor(1, 2)+1, nt), ity*nt-nt >= ijk_neighbor(1, 2)),merge(nt, mod(ijk_neighbor(2, 2)+1, nt)-1, ity*nt <= ijk_neighbor(2, 2))
          do i = merge(1, mod(ijk_neighbor(1, 1)+1, nt), itx*nt-nt >= ijk_neighbor(1, 1)),merge(nt, mod(ijk_neighbor(2, 1)+1, nt)-1, itx*nt <= ijk_neighbor(2, 1))
            np_neighbors(iq1,iq2,iq3) = np_neighbors(iq1,iq2,iq3) + rhoc_local(i,j,k,itx,ity,itz)
          enddo
          enddo
          enddo
        enddo
        enddo
        enddo
        ! write(log_unit, *)iq1-2,iq2-2,iq3-2,np_neighbors(iq1,iq2,iq3)!,np_neighbors(iq1,iq2,iq3)*1d0/(ijk_neighbor(2, 1) -  ijk_neighbor(1, 1))/ (ijk_neighbor(2, 2) -  ijk_neighbor(1, 2))/  (ijk_neighbor(2, 3) -  ijk_neighbor(1, 3))
      enddo
      enddo
      enddo
      np_max = sum(np_neighbors)





      !  initialize particles
      write(log_unit, *)' initialize particles', np_max
      jp = 0
      allocate(xv(3,np_max))
      do iq1 = 1, 3
      do iq2 = 1, 3
      do iq3 = 1, 3
        ijk_neighbor(:, 1) = neighbor_b(:, 4-iq1)
        ijk_neighbor(:, 2) = neighbor_b(:, 4-iq2)
        ijk_neighbor(:, 3) = neighbor_b(:, 4-iq3)

        inx=modulo(icx+1-iq1,nn)+1
        iny=modulo(icy+1-iq2,nn)
        inz=modulo(icz+1-iq3,nn)
        image = inz*nn**2+iny*nn+inx
        write(log_unit, '(I4,I4,I4,I4)')image,inx,iny,inz
        open(11,file=output_name('np'),access='stream'); read(11) rhoc_local; close(11)
        open(11,file=output_name('info'),access='stream'); read(11) sim; close(11)
        sim%cur_checkpoint=cur_checkpoint
        allocate(xp(3,sim%nplocal)[nn,nn,*])
        open(11,file=output_name('xp'),access='stream'); read(11) xp; close(11)
        ! 循环遍历各个维度
        do itz = floor((ijk_neighbor(1, 3))*1d0/nt)+1, floor((ijk_neighbor(2, 3)-1)*1d0/nt)+1
        do ity = floor((ijk_neighbor(1, 2))*1d0/nt)+1, floor((ijk_neighbor(2, 2)-1)*1d0/nt)+1
        do itx = floor((ijk_neighbor(1, 1))*1d0/nt)+1, floor((ijk_neighbor(2, 1)-1)*1d0/nt)+1
          do k = merge(1, mod(ijk_neighbor(1, 3)+1, nt), itz*nt-nt >= ijk_neighbor(1, 3)),merge(nt, mod(ijk_neighbor(2, 3)+1, nt)-1, itz*nt <= ijk_neighbor(2, 3))
          do j = merge(1, mod(ijk_neighbor(1, 2)+1, nt), ity*nt-nt >= ijk_neighbor(1, 2)),merge(nt, mod(ijk_neighbor(2, 2)+1, nt)-1, ity*nt <= ijk_neighbor(2, 2))
          do i = merge(1, mod(ijk_neighbor(1, 1)+1, nt), itx*nt-nt >= ijk_neighbor(1, 1)),merge(nt, mod(ijk_neighbor(2, 1)+1, nt)-1, itx*nt <= ijk_neighbor(2, 1))
            np=rhoc_local(i,j,k,itx,ity,itz)
            nlast = sum(rhoc_local(:,:,:,:,:,:itz-1))     &
                  + sum(rhoc_local(:,:,:,:,:ity-1,itz))   &
                  + sum(rhoc_local(:,:,:,:itx-1,ity,itz)) &
                  + sum(rhoc_local(:,:,:k-1,itx,ity,itz)) &
                  + sum(rhoc_local(:,:j-1,k,itx,ity,itz)) &
                  + sum(rhoc_local(:i-1,j,k,itx,ity,itz))
            do l=1,np
              ip = nlast+l
              jp = jp+1
              if (ip < 1 .or.  ip > sim%nplocal .or. jp > np_max) then
                print*, image
                print*, i,j,k,itx,ity,itz
                print*,nlast,l,np,ip,jp
                stop 'particle index error'
              endif
#ifdef ZIPX
              xv(:,jp)=nt*((/itx,ity,itz/)-1)+ ((/i,j,k/)-1) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution+[fof_buffer,fof_buffer,fof_buffer]  + ijk_neighbor(3,:)
#else 
              xv(:,jp)=xp(:,ip)/ratio_cs+[fof_buffer,fof_buffer,fof_buffer]  + ijk_neighbor(3,:)
#endif
            enddo
          enddo
          enddo
          enddo
        enddo
        enddo
        enddo
        deallocate(xp)
        ! write(log_unit, *)iq1-2,iq2-2,iq3-2,np_neighbors(iq1,iq2,iq3)!,np_neighbors(iq1,iq2,iq3)*1d0/(ijk_neighbor(2, 1) -  ijk_neighbor(1, 1))/ (ijk_neighbor(2, 2) -  ijk_neighbor(1, 2))/  (ijk_neighbor(2, 3) -  ijk_neighbor(1, 3))
      enddo
      enddo
      enddo
      deallocate(rhoc_local)
      image = image_now


      ! create hoc ll
      call system_clock(ft(2,image_now),ftr(image_now))
      write(log_unit, *) '    real time =',real(ft(2,image_now)-ft(1,image_now))/ftr(image_now),'secs';write(log_unit, *) ''



      write(log_unit, *)' create hoc ll'
      allocate(hoc(nfof,nfof,nfof),ll(np_max),llgp(np_max),hcgp(np_max),ecgp(np_max))
      hoc=0; ll=0
      do ip=1,np_max
        idx=floor(xv(1:3,ip)*n_refine)+1 ! index of the grid
        ll(ip)=hoc(idx(1),idx(2),idx(3)) ! linked list 
        hoc(idx(1),idx(2),idx(3))=ip ! head of chain
        hcgp(ip)=ip ! initialize hcgp(ip)=ip for isolated particles
      enddo
      llgp=0; ecgp=0; ! initialize group link list
      
      ! rho
      if (0) then
        allocate(rho_grid(0:nfof+1,0:nfof+1,0:nfof+1)[nn,nn,*])
        rho_grid=0
        do iq3=1,nfof
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
        rho8 = sum(rho_grid(1:nfof,1:nfof,1:nfof))/(nfof**3)
        rho_grid = rho_grid/rho8-1
        write(log_unit, *)'min',minval(rho_grid(1:nfof,1:nfof,1:nfof)),'max',maxval(rho_grid(1:nfof,1:nfof,1:nfof)),'mean',sum(rho_grid(1:nfof,1:nfof,1:nfof)*1d0)/nfof/nfof/nfof
        write(log_unit, *)'Write delta_fof into',output_name('delta_fof_'//trim(adjustl(str_refine)))
        open(11,file=output_name('delta_fof_'//trim(adjustl(str_refine))),status='replace',access='stream')
        write(11) rho_grid(1:nfof,1:nfof,1:nfof)
        close(11)
        deallocate(rho_grid)
      endif


      call system_clock(ft(3,image_now),ftr(image_now))
      write(log_unit, *) '    real time =',real(ft(3,image_now)-ft(2,image_now))/ftr(image_now),'secs';write(log_unit, *) ''
    
      ! loop over fof cells
      write(log_unit, *) 'Loop over fof cells'
      do iq3=1,nfof
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
      deallocate(hoc,ll,ecgp)

      call system_clock(ft(4,image_now),ftr(image_now))
      write(log_unit, *) '    real time =',real(ft(4,image_now)-ft(3,image_now))/ftr(image_now),'secs';write(log_unit, *) ''


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
      deallocate(xv)
      call system_clock(ft(5,image_now),ftr(image_now))
      write(log_unit, *) '    real time =',real(ft(5,image_now)-ft(4,image_now))/ftr(image_now),'secs';write(log_unit, *) ''


      nh=np_head
      halo_header%nhalo_tot=0
      halo_header%nhalo=nh; halo_header%ninfo=ninfo; halo_header%linking_parameter=b_link
      write(log_unit, *) output_name(trim('halo_'//trim(adjustl(str_refine))))
      open(21,file=output_name(trim('halo_'//trim(adjustl(str_refine)))),status='replace',access='stream')
      write(21) halo_header

      allocate(hcat(nh)[*])
      hcat%hmass = np_halo(1:nh)
      do i=1,3
        hcat%xv(i)=xv_mean(i,1:nh)
      enddo


      write(21) hcat
      close(21)
      write(log_unit, *)'halo_header',halo_header,maxval(hcat%hmass)

      write(log_unit, *) output_name(trim('halo_mass_'//trim(adjustl(str_refine))))
      open(21,file=output_name(trim('halo_mass_'//trim(adjustl(str_refine)))),status='replace',access='stream')
      write(21) nh
      write(21) np_halo(1:nh)
      close(21)
      deallocate(np_halo,hcgp,llgp,xv_mean,np_halo_all)

      deallocate(hcat)
      halo_images(image_now) = np_head


      write(log_unit, *) np_head,'all time =',real(ft(5,image_now)-ft(1,image_now))/ftr(image_now),'secs'

      write(*,'(A, I0, A, A, I0, A, I0, A, F7.3, A)') '  Fof at image',image_now,'->'//trim(log_filename),' |  np_max', np_max,' |  find halo :',np_head,' |  used time =',real(ft(5,image_now)-ft(1,image_now))/ftr(image_now),'secs'
      
      close(log_unit)
    enddo



    ! write halo mass
    sync all
    if (head) then
      print*,''
      print*,''
      print*,''
      do image=2,real_images
        halo_images = halo_images + halo_images(:)[image]
      enddo
      nhalo_all = sum(halo_images)
      print*,nhalo_all
      do image=1,nn**3
        print*,image,halo_images(image)
      enddo
    endif
    sync all
    if (head) then
      image = 1
      write(*, *) output_name(trim('halo_mass_all_'//trim(adjustl(str_refine))))
      open(111,file=output_name(trim('halo_mass_all_'//trim(adjustl(str_refine)))),status='replace',access='stream')
      write(21) nhalo_all
      do image=1,nn**3
        write(*, *) output_name(trim('halo_mass_'//trim(adjustl(str_refine))))
        open(21,file=output_name(trim('halo_mass_'//trim(adjustl(str_refine)))),status='old',access='stream')
        read(21) np_head; allocate(np_halo_all(np_head)); read(21) np_halo_all; close(21)
        write(111) np_halo_all
        deallocate(np_halo_all)
      enddo
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

  subroutine geometry_images
    write(log_unit, *) 'geometry_images'
    rank=image_now-1               ! MPI_rank
    icz=rank/(nn**2)+1             ! image_z
    icy=(rank-nn**2*(icz-1))/nn+1  ! image_y
    icx=mod(rank,nn)+1             ! image_x
  endsubroutine
  
end