program fof_time_slide
   use variables
   implicit none

   character(len = 4) str_refine
   integer,parameter:: fofcore = ncore
   integer(8),parameter:: fof_buffer = ceiling(20.0/box*nc*nn)
   integer(8),parameter:: nfof = (nc+2*fof_buffer)*ratio_cs
   real(8),parameter:: n_refine = nfof*1d0/(nc + fof_buffer*2 + b_link/ratio_cs)
   real(8),parameter:: L_b    = fof_buffer
   real(8),parameter:: L_bL   = L_b  + nc
   real(8),parameter:: L_bLb  = L_bL + fof_buffer
   real(8),parameter:: L_bLb2 = L_bLb/2
   real(8),parameter:: rp2=(b_link/ratio_cs)**2
   integer,parameter:: N_mesh = 1000
   logical,parameter :: bound_check = .true.

   integer :: log_unit
   character(len=64) :: log_filename

   integer total_images,layer_image,image_now,i,j,k,l,nlayer,layer_core
   integer cur_checkpoint,iz_checkpoint,np,n1,n2,idx3(3),nh,im
   integer idx1(3),idx2(3),ft(6,nn*3),ftr(nn*3)
   integer(8) iq1,iq2,iq3,i_neighbor,jq(3),ijk_neighbor(3,3),neighbor_b(3,3)
   integer(4) nlast,ip,jp,np_iso(fofcore),np_head(fofcore),np_max
   integer(4) np_neighbors(3,3,3)[nn,nn,*],np_need(3,3,3)
   integer(4) max_nei[nn,nn,*],offset_nei(3,3,3),offset_team(fofcore)
   integer(4),allocatable :: rhoc_local(:,:,:,:,:,:),offset_map(:,:,:,:,:,:)
   integer(4),allocatable :: hoc(:,:,:),ll(:),llgp(:),hcgp(:),ecgp(:)
   real rsq,pos1(3),dx1(3),dx2(3),shift_xv(3)
   real(8) rho8,dxv(3)
   real,allocatable :: xp_all(:,:)
   integer(8),allocatable :: pid_all(:)
   real,allocatable :: vp_all(:,:)

   ! halo粒子相关变量
   type(type_halo_particle) halo_particles
   integer,allocatable :: i_start_array(:)
   integer(8),allocatable :: PID_tmp(:), idx_to_zin(:)
   real,allocatable :: xv_halo(:,:)
   integer np_all, np_halo_total
   integer nz_current
   integer np_fof, np_temp, np_bound
   real,allocatable :: z_list_tmp(:)
   real,allocatable :: xv_mean_tmp(:,:,:)
   real dx(3), dv(3), r, v_rel_sq
   integer np_found, idx

   ! 局部FOF变量
   integer istart, ihalo
   integer jp_prev

   ! 网格+链表变量
   integer N_x, N_y, N_z
   real L_mesh
   real :: xp_zero(3)
   integer,allocatable :: hoc_local(:,:,:), ll_local(:)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   total_images = num_images()
   layer_image = nn**3/total_images
   if (total_images*layer_image /= nn**3) then
      stop 'total_images*layer_image /= nn*3'
   endif
   ft = 0

   nlayer = ceiling(fof_buffer*n_refine*2)+1
   do while (mod(nfof,nlayer) /= 0)
      nlayer = nlayer + 1
   enddo

   head=(this_image()==1)
   call omp_set_num_threads(fofcore)
   layer_core = min(fofcore,nfof/nlayer)

   if(head) call system("mkdir -p "//opath//'/fof')

   write(str_refine,'(i4)')  nfof
   if(head) print*,nfof,nlayer,layer_core,total_images,layer_image

   ! 读取z_checkpoint
   if (head) then
      open(16,file='z_checkpoint.txt',status='old')
      do i=1,nmax_redshift-1
         read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
      enddo
71    n_checkpoint=i-1
      close(16)
      if (n_checkpoint==0) stop 'z_checkpoint.txt empty'
   endif
   sync all
   n_checkpoint=n_checkpoint[1]
   z_checkpoint(:)=z_checkpoint(:)[1]
   sync all

   ! 每个image独立读取halo_particles文件
   write(*,'(A, I0, A)') '  Reading halo particles file at image', this_image()
   image = this_image()
   open(11,file=output_name(trim('halo_particles_'//trim(adjustl(str_refine)))),status='old',access='stream')
   read(11) halo_particles%box
   read(11) halo_particles%linking_parameter
   read(11) halo_particles%nh
   read(11) halo_particles%nz
   read(11) halo_particles%np_halo_all
   allocate(halo_particles%np_halo(halo_particles%nh))
   allocate(halo_particles%xv_mean(6,halo_particles%nh,halo_particles%nz))
   allocate(halo_particles%iq(halo_particles%np_halo_all))
   allocate(halo_particles%z_list(halo_particles%nz))
   allocate(halo_particles%z_in(halo_particles%np_halo_all))
   read(11) halo_particles%np_halo
   read(11) halo_particles%xv_mean
   read(11) halo_particles%iq
   read(11) halo_particles%z_list
   read(11) halo_particles%z_in
   close(11)

   ! 初始化
   np_all = halo_particles%np_halo_all
   allocate(i_start_array(halo_particles%nh))
   allocate(PID_tmp(np_all), idx_to_zin(np_all))
   allocate(xv_halo(6, np_all))
   xv_halo = 0.0

   ! 构建i_start_array
   i_start_array(1) = 1
   do ihalo = 2, halo_particles%nh
      i_start_array(ihalo) = i_start_array(ihalo-1) + halo_particles%np_halo(ihalo-1)
   enddo

   ! 初始化idx_to_zin
   do i = 1, np_all
      idx_to_zin(i) = i
   enddo

   ! 复制PID到PID_tmp
   do i = 1, np_all
      PID_tmp(i) = int(halo_particles%iq(i), kind=4)
   enddo

   ! 按PID排序
   call sort_by(PID_tmp, idx_to_zin, 1, np_all)

   ! 预分配z_list空间
   deallocate(halo_particles%z_list)
   allocate(halo_particles%z_list(n_checkpoint))
   halo_particles%nz = 0

   if(head) print*, 'Starting time slide from z=', z_checkpoint(n_checkpoint), ' to z=', z_checkpoint(1)

   ! 遍历红移（从低到高）
   do iz_checkpoint = n_checkpoint, 1, -1
      cur_checkpoint = iz_checkpoint
      sim%cur_checkpoint = cur_checkpoint

      ! 在遍历红移处更新z_list
      ! 查找当前红移是否已存在于z_list中
      nz_current = 0
      do i = 1, halo_particles%nz
         if (abs(halo_particles%z_list(i) - z_checkpoint(cur_checkpoint)) < 1d-6) then
            nz_current = i
            exit
         endif
      enddo
      
      ! 如果当前红移不存在于z_list中，则添加
      if (nz_current == 0) then
         halo_particles%nz = halo_particles%nz + 1

         if (halo_particles%nz > size(halo_particles%z_list)) then
            allocate(z_list_tmp(halo_particles%nz))
            z_list_tmp(1:halo_particles%nz-1) = halo_particles%z_list(:)
            z_list_tmp(halo_particles%nz) = z_checkpoint(cur_checkpoint)
            deallocate(halo_particles%z_list)
            allocate(halo_particles%z_list(halo_particles%nz))
            halo_particles%z_list(:) = z_list_tmp(:)
            deallocate(z_list_tmp)
            
            ! 重新分配xv_mean
            allocate(xv_mean_tmp(6, halo_particles%nh, halo_particles%nz))
            xv_mean_tmp(:,:,1:halo_particles%nz-1) = halo_particles%xv_mean(:,:,:)
            xv_mean_tmp(:,:,halo_particles%nz) = 0.0
            deallocate(halo_particles%xv_mean)
            allocate(halo_particles%xv_mean(6, halo_particles%nh, halo_particles%nz))
            halo_particles%xv_mean(:,:,:) = xv_mean_tmp(:,:,:)
            deallocate(xv_mean_tmp)
         else
            halo_particles%z_list(halo_particles%nz) = z_checkpoint(cur_checkpoint)
         endif
         nz_current = halo_particles%nz
      endif

      if (head) print*,  ''
      if (head) print*,  'Processing redshift ', z2str(z_checkpoint(cur_checkpoint)), ' nz=', nz_current

      do image_now = this_image(),nn**3,total_images

         write(log_filename, "(A,'/fof/track_output_', I4.4, '.log')") opath,image_now
         open(newunit=log_unit, file=trim(log_filename), status='replace', action='write')
         write(*,'(A, I0, A)') '  Tracking at image',image_now,'->'//trim(log_filename)
         write(log_unit,'(A, I0)') '  Tracking at image',image_now

         ! init image
         call geometry_images

         ! count neighbor particle
         write(log_unit, *)' count  neighbor particle '
         allocate(rhoc_local(nt,nt,nt,nnt,nnt,nnt))
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

                  open(11,file=output_name('np'),access='stream'); read(11) rhoc_local; close(11)
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
               enddo
            enddo
         enddo
         deallocate(rhoc_local)
         np_max = sum(np_neighbors)

         ! initialize particles
         write(log_unit, *)' initialize particles', np_max
         allocate(xp_all(3,np_max))
         allocate(pid_all(np_max))
         allocate(vp_all(3,np_max))
         xp_all = 0
         pid_all = 0
         vp_all = 0
         nlast = 0
         do iq1 = 1, 3;  do iq2 = 1, 3; do iq3 = 1, 3
            offset_nei(iq1,iq2,iq3) =  nlast
            nlast = nlast + np_neighbors(iq1,iq2,iq3)
         enddo; enddo; enddo

         shift_xv  = [(icx-1)*nc*1d0-fof_buffer,(icy-1)*nc*1d0-fof_buffer,(icz-1)*nc*1d0-fof_buffer]
         neighbor_b = reshape([int(0,kind=8),fof_buffer,nc  ,int(0,kind=8),nc,int(0,kind=8) ,nc-fof_buffer,nc,-nc ],[3,3])

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

            iteam=omp_get_thread_num()+111
            open(iteam,file=output_name('info'),access='stream'); read(iteam) sim; close(iteam)
            sim%cur_checkpoint=cur_checkpoint
            allocate(xp_new(3,sim%nplocal),rhoc_local(nt,nt,nt,nnt,nnt,nnt),offset_map(nt,nt,nt,nnt,nnt,nnt))
#ifdef PID
            allocate(pid_new(sim%nplocal))
#endif
#ifdef ZIPV
            allocate(vp_new(3,sim%nplocal))
#endif
            open(iteam,file=output_name('xp'),access='stream'); read(iteam) xp_new; close(iteam)
            open(iteam,file=output_name('np'),access='stream'); read(iteam) rhoc_local; close(iteam)
#ifdef PID
            open(iteam,file=output_name('pid'),access='stream'); read(iteam) pid_new; close(iteam)
#endif
#ifdef ZIPV
            open(iteam,file=output_name('vp'),access='stream'); read(iteam) vp_new; close(iteam)
#endif

            nlast = 0
            do itz=1,nnt; do ity=1,nnt; do itx=1,nnt; do k=1,nt; do j=1,nt; do i=1,nt
                              np = rhoc_local(i,j,k,itx,ity,itz)
                              if (np < 0 ) then
                                 write(log_unit, '(7I4,I10)')  i,j,k,itx,ity,itz,image,np
                                 error stop 'particle index error'
                              endif
                              offset_map(i,j,k,itx,ity,itz) = nlast
                              nlast = nlast + np
                           enddo; enddo; enddo; enddo; enddo; enddo

            jp = offset_nei(iq1,iq2,iq3)
            n1 = jp+1
            do itz = floor((ijk_neighbor(1, 3))*1d0/nt)+1, floor((ijk_neighbor(2, 3)-1)*1d0/nt)+1
            do ity = floor((ijk_neighbor(1, 2))*1d0/nt)+1, floor((ijk_neighbor(2, 2)-1)*1d0/nt)+1
            do itx = floor((ijk_neighbor(1, 1))*1d0/nt)+1, floor((ijk_neighbor(2, 1)-1)*1d0/nt)+1
               do k = merge(1, mod(ijk_neighbor(1, 3)+1, nt), itz*nt-nt >= ijk_neighbor(1, 3)),merge(nt, mod(ijk_neighbor(2, 3)+1, nt)-1, itz*nt <= ijk_neighbor(2, 3))
               do j = merge(1, mod(ijk_neighbor(1, 2)+1, nt), ity*nt-nt >= ijk_neighbor(1, 2)),merge(nt, mod(ijk_neighbor(2, 2)+1, nt)-1, ity*nt <= ijk_neighbor(2, 2))
               do i = merge(1, mod(ijk_neighbor(1, 1)+1, nt), itx*nt-nt >= ijk_neighbor(1, 1)),merge(nt, mod(ijk_neighbor(2, 1)+1, nt)-1, itx*nt <= ijk_neighbor(2, 1))
                  np=rhoc_local(i,j,k,itx,ity,itz)
                  nlast = offset_map(i,j,k,itx,ity,itz)
#ifdef ZIPX
                  xp_all(:,jp+1:jp+np)=(int(xp_new(:,nlast+1:nlast+np)+ishift,izipx)+rshift)*x_resolution &
                     +spread(nt*((/itx,ity,itz/)-1)+((/i,j,k/)-1)+[fof_buffer,fof_buffer,fof_buffer]+ijk_neighbor(3,:),dim=2,ncopies=np)
#else
                  xp_all(:,jp+1:jp+np)=xp_new(:,nlast+1:nlast+np)+spread([fof_buffer,fof_buffer,fof_buffer]+ijk_neighbor(3,:),dim=2,ncopies=np)
#endif
#ifdef PID
                  pid_all(jp+1:jp+np) = pid_new(nlast+1:nlast+np)
#endif
#ifdef ZIPV
                  vp_all(:,jp+1:jp+np) = (int(vp_new(:,nlast+1:nlast+np)+ishift,izipv)+rshift)*v_resolution
#else
                  vp_all(:,jp+1:jp+np) = vp_new(:,nlast+1:nlast+np)
#endif
                  jp = jp+np
               enddo
               enddo
               enddo
            enddo
            enddo
            enddo
            deallocate(xp_new,rhoc_local,offset_map)
#ifdef PID
            deallocate(pid_new)
#endif
#ifdef ZIPV
            deallocate(vp_new)
#endif
            n2 = jp-n1+1
            if (n2 /= np_neighbors(iq1,iq2,iq3) ) then
               write(log_unit, *)  iq1-2,iq2-2,iq3-2,image
               write(log_unit, *) np_neighbors(iq1,iq2,iq3),n2,jp,n1
               write(log_unit, *) sum(xp_all(1,n1:jp))/n2,minval(xp_all(1,n1:jp)),maxval(xp_all(1,n1:jp)),ijk_neighbor(3,1)
               write(log_unit, *) sum(xp_all(2,n1:jp))/n2,minval(xp_all(2,n1:jp)),maxval(xp_all(2,n1:jp)),ijk_neighbor(3,2)
               write(log_unit, *) sum(xp_all(3,n1:jp))/n2,minval(xp_all(3,n1:jp)),maxval(xp_all(3,n1:jp)),ijk_neighbor(3,3)
               write(log_unit, *)  'particle len error'
               print*,image
               stop 'particle len error'
            endif
         enddo
         enddo
         enddo
         image = image_now

         ! 提取halo粒子到xv_halo
         write(log_unit, *) 'Extracting halo particles'
         allocate(hcgp(np_all))
         allocate(llgp(np_all))
         hcgp = 0
         llgp = 0

         do ip = 1, np_max
            idx = pid_to_idx(pid_all(ip))
            if (idx > 0) then
               ihalo = find_halo(idx, i_start_array, halo_particles%nh)
               xv_halo(1:3, idx) = modulo(xp_all(:,ip) - halo_particles%xv_mean(1:3,ihalo,1) + L_bLb2, L_bLb*1d0) - L_bLb2
               xv_halo(4:6, idx) = vp_all(:,ip) - halo_particles%xv_mean(4:6,ihalo,1)
               hcgp(idx) = idx
            endif
         enddo

         ! 检查是否所有粒子都被找到
         np_found = count(hcgp /= 0)
         if (np_found /= np_all) then
            write(log_unit, *) 'Error: particle count mismatch, found', np_found, 'expected', np_all
            call exit(1)
         endif

         deallocate(xp_all, pid_all, vp_all)

         ! 执行局部FOF
         write(log_unit, *) 'Performing local FOF'
         L_mesh = b_link * 2

         ! 遍历每个halo
         do ihalo = 1, halo_particles%nh
            istart = i_start_array(ihalo)
            np = halo_particles%np_halo(ihalo)

            if (np < N_mesh) then
               ! 直接遍历当前halo的粒子
               do i = istart, istart+np-1
                  do j = i+1, istart+np-1
                     rsq = sum((xv_halo(1:3,i) - xv_halo(1:3,j))**2)
                     if (rsq < rp2) call merge_chain(i, j)
                  enddo
               enddo
            else
               ! 网格+链表，只循环当前halo的粒子
               N_x = ceiling((maxval(xv_halo(1,istart:istart+np-1))-minval(xv_halo(1,istart:istart+np-1))) / L_mesh)
               N_y = ceiling((maxval(xv_halo(2,istart:istart+np-1))-minval(xv_halo(2,istart:istart+np-1))) / L_mesh)
               N_z = ceiling((maxval(xv_halo(3,istart:istart+np-1))-minval(xv_halo(3,istart:istart+np-1))) / L_mesh)
               xp_zero = halo_particles%xv_mean(1:3,ihalo,1) - [minval(xv_halo(1,istart:istart+np-1)), &
                                                           minval(xv_halo(2,istart:istart+np-1)), &
                                                           minval(xv_halo(3,istart:istart+np-1))]

               allocate(hoc_local(N_x,N_y,N_z), ll_local(np))
               hoc_local = 0
               ll_local = 0

               ! 构建局部网格，考虑周期性边界条件
               do i = 1, np
                  idx3 = floor((xv_halo(1:3,istart+i-1) - xp_zero) / L_mesh) + 1
                  ll_local(i) = hoc_local(idx3(1), idx3(2), idx3(3))
                  hoc_local(idx3(1), idx3(2), idx3(3)) = istart+i-1
               enddo

               ! FOF搜索
               do iq3 = 1, N_z
                  do iq2 = 1, N_y
                     do iq1 = 1, N_x
                        ip = hoc_local(iq1, iq2, iq3)
                        do while (ip /= 0)
                           ! 检查当前粒子邻居
                           jp = ll_local(ip-istart+1)
                           do while (jp /= 0)
                              rsq = sum((xv_halo(1:3,ip) - xv_halo(1:3,jp))**2)
                              if (rsq <= rp2) call merge_chain(ip, jp)
                              jp = ll_local(jp-istart+1)
                           enddo
                           ! 检查周围27个格网单元
                           do i_neighbor = 1, 13
                              jq = [iq1,iq2,iq3] + ijk(:,i_neighbor)
                              if (maxval(jq) > N_x .or. minval(jq) < 1) cycle
                              jp = hoc_local(jq(1), jq(2), jq(3))
                              do while (jp /= 0)
                                 rsq = sum((xv_halo(1:3,ip) - xv_halo(1:3,jp))**2)
                                 if (rsq <= rp2) call merge_chain(ip, jp)
                                 jp = ll_local(jp-istart+1)
                              enddo
                           enddo
                           ip = ll_local(ip-istart+1)
                        enddo
                     enddo
                  enddo
               enddo

               deallocate(hoc_local, ll_local)
            endif

            ! 统计halo并更新z_in
            do i = istart, istart+np-1
               if (hcgp(i) /= i .and. hcgp(i) /= 0) then
                  ! 找到FOF链条的头部
                  ip = hcgp(i); np_fof = 0; jp = ip
                  do while (jp /= 0)
                     np_fof = np_fof + 1
                     jp = llgp(jp)
                  enddo
                  
                  if (np_fof >= np_halo_min) then
                     ! 计算该FOF链条的中心和速度
                     dxv = 0
                     jp = ip
                     do while (jp /= 0)
                        dxv = dxv + xv_halo(:,jp)
                        jp = llgp(jp)
                     enddo
                     dxv(1:3) = modulo(dxv(1:3)/np_fof + halo_particles%xv_mean(1:3,ihalo,1), L_bLb*1d0)
                     dxv(4:6) = dxv(4:6)/np_fof
                     
                     if (bound_check) then
                        ! 束缚检查：丢弃不束缚的粒子，更新链表
                        np_bound = 0
                        jp_prev = 0
                        jp = ip
                        do while (jp /= 0)
                           dx = modulo((xv_halo(1:3, jp) - dxv(1:3)) + L_bLb2, L_bLb*1d0) - L_bLb2
                           dv = xv_halo(4:6, jp) - dxv(4:6)
                           r = norm2(dx)
                           v_rel_sq = sum(dv**2)
                           if (0.5 * v_rel_sq < np_fof / r) then
                              ! 束缚粒子
                              np_bound = np_bound + 1
                              jp_prev = jp
                              jp = llgp(jp)
                           else
                              ! 不束缚粒子，更新链表
                              if (jp == ip) then
                                 ! 头部粒子不束缚，更新头部
                                 ip = llgp(jp)
                                 hcgp(i) = ip
                                 hcgp(ip) = ip
                                 jp = ip
                              else
                                 ! 非头部粒子不束缚，更新前一个粒子的链表
                                 llgp(jp_prev) = llgp(jp)
                                 jp = llgp(jp_prev)
                              endif
                           endif
                        enddo
                        
                        ! 如果有粒子被抛弃，重新计算dxv
                        if (np_fof > np_halo_min .and. np_bound < np_fof) then
                           np_fof = np_bound
                           dxv = 0
                           jp = ip
                           do while (jp /= 0)
                              dxv = dxv + xv_halo(:,jp)
                              jp = llgp(jp)
                           enddo
                           dxv(1:3) = modulo(dxv(1:3)/np_fof + halo_particles%xv_mean(1:3,ihalo,1), L_bLb*1d0)
                           dxv(4:6) = dxv(4:6)/np_fof
                        endif
                     endif
                     
                     ! 更新粒子的z_in和xv_mean
                     if (np_fof >= np_halo_min) then
                        halo_particles%xv_mean(:,ihalo,nz_current) = dxv
                        jp = ip
                        do while (jp /= 0)
                           halo_particles%z_in(jp) = nz_current
                           jp = llgp(jp)
                        enddo
                     endif
                  endif

               endif
            enddo
         enddo

         deallocate(hcgp, llgp)

         close(log_unit)
      enddo

      ! 输出更新后的halo_particles文件
      if (head) then
         image = 1
         write(*, *) 'Writing halo_track file at z=', z_checkpoint(cur_checkpoint)
         open(11,file=output_name(trim('halo_track_zin')),status='replace',access='stream')
         write(11) halo_particles%box
         write(11) halo_particles%linking_parameter
         write(11) halo_particles%nh
         write(11) halo_particles%nz
         write(11) halo_particles%np_halo_all
         write(11) halo_particles%np_halo
         write(11) halo_particles%xv_mean
         write(11) halo_particles%iq
         write(11) halo_particles%z_list
         write(11) halo_particles%z_in
         close(11)
      endif

   enddo

   ! 清理
   deallocate(i_start_array, PID_tmp, idx_to_zin, xv_halo)
   if(allocated(halo_particles%np_halo)) deallocate(halo_particles%np_halo)
   if(allocated(halo_particles%xv_mean)) deallocate(halo_particles%xv_mean)
   if(allocated(halo_particles%iq)) deallocate(halo_particles%iq)
   if(allocated(halo_particles%z_list)) deallocate(halo_particles%z_list)
   if(allocated(halo_particles%z_in)) deallocate(halo_particles%z_in)

contains

   recursive subroutine sort_by(arr, brr, left, right)
      integer(8), intent(inout) :: arr(:)
      integer(8), intent(inout) :: brr(:)
      integer, intent(in) :: left, right

      integer(8) :: key_arr,tmp
      integer(8) :: key_brr
      integer :: i, j, mid

      if (left >= right) return

      mid = (left + right) / 2
      key_arr = arr(mid)
      key_brr = brr(mid)

      i = left; j = right
      do while (i <= j)
         do while (arr(i) < key_arr)
            i = i + 1
         enddo
         do while (arr(j) > key_arr)
            j = j - 1
         enddo
         if (i <= j) then
            tmp = arr(i); arr(i) = arr(j); arr(j) = tmp
            tmp = brr(i); brr(i) = brr(j); brr(j) = tmp
            i = i + 1
            j = j - 1
         endif
      enddo

      if (left < j) call sort_by(arr, brr, left, j)
      if (i < right) call sort_by(arr, brr, i, right)
   end subroutine

   function pid_to_idx(pid_val) result(idx)
      integer(8), intent(in) :: pid_val
      integer :: idx

      integer :: lo, mid, hi
      lo = 1; hi = np_all

      do while (lo <= hi)
         mid = (lo + hi) / 2
         if (PID_tmp(mid) == int(pid_val, kind=4)) then
            idx = idx_to_zin(mid)
            return
         else if (PID_tmp(mid) < int(pid_val, kind=4)) then
            lo = mid + 1
         else
            hi = mid - 1
         endif
      enddo

      idx = 0
   end function

   function find_halo(idx, i_start_array, nhalo) result(ihalo)
      integer, intent(in) :: idx, nhalo
      integer, intent(in) :: i_start_array(nhalo)
      integer :: ihalo

      integer :: lo, mid, hi
      lo = 1; hi = nhalo

      do while (lo < hi)
         mid = (lo + hi) / 2
         if (mid+1 <= nhalo .and. idx >= i_start_array(mid+1)) then
            lo = mid + 1
         else
            hi = mid
         endif
      enddo
      ihalo = lo
   end function

   subroutine merge_chain(ii, jj)
      integer(4), intent(in) :: ii, jj
      integer(4) :: ihead, jhead, iend, jend, ipart

      ihead = hcgp(ii)
      jhead = hcgp(jj)
      if (ihead == jhead) return

      ! 找到j链的末尾
      jend = jj
      do while (llgp(jend) /= 0)
         jend = llgp(jend)
      enddo

      ! 找到i链的末尾
      iend = ii
      do while (llgp(iend) /= 0)
         iend = llgp(iend)
      enddo

      ! 链接j链到i链
      llgp(jend) = ihead

      ! 更新j链中所有粒子的hcgp
      ipart = jhead
      do while (ipart /= 0)
         hcgp(ipart) = ihead
         ipart = llgp(ipart)
      enddo
   end subroutine

   subroutine geometry_images
      write(log_unit, *) 'geometry_images'
      rank=image_now-1
      icz=rank/(nn**2)+1
      icy=(rank-nn**2*(icz-1))/nn+1
      icx=mod(rank,nn)+1
   end subroutine

end program fof_time_slide
