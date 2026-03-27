# FOF 时间追踪程序设计方案

***

## 一、总体架构

本方案包含两个独立程序：

| 程序                 | 功能                                       | 输入                      | 输出                      |
| :------------------- | :----------------------------------------- | :------------------------ | :------------------------ |
| `fof.f90`            | 从完整模拟数据中识别所有halo，保存粒子链表 | 粒子位置数据              | halo目录 + 粒子链表       |
| `fof_time_slide.f90` | 追踪halo粒子在不同红移的束缚情况           | halo粒子列表 + 各红移数据 | 粒子束缚历史（z\_in数组） |

***

## 二、数据结构

使用统一的halo数据结构存储粒子信息：

```fortran
type type_halo_particle
   real box ! box size
   real linking_parameter
   real z_fof_ini
   integer nh
   integer nz
   integer np_halo_all ! number of particles in all halos ,sum(np_halo)
   integer, allocatable :: np_halo(:) ! np_halo(nh) number of particles in this halo

   real, allocatable :: xv_mean(:,:,:) ! xv_mean(6,nh,nz) position of particles in this halo
   integer(8), allocatable :: PID_halo(:) ! PID_halo(np_halo_all) PID of particles in this halo
   real, allocatable :: z_list(:) ! z_list(nz) number of z_checkpoint
   integer, allocatable :: z_in(:) ! z_in(np_halo_all) particle fall into halo in z_list(z_in)
   ! optional
   real, allocatable :: xv_z(:,:,:) ! xv_z(6,nz,np_halo_all) position of particles in this halo in each z_checkpoint, nan if particle not in this halo in this z_checkpoint
endtype
```



***

## 三、fof.f90 更新

### 3.1 输出文件

| 文件名                    | 说明                   |
| :------------------------ | :--------------------- |
| `halo_particles_<refine>` | 每个halo包含的粒子链表 |

### 3.2 输出格式

根据 `type_halo_particle` 结构体定义，输出文件包含以下内容：

```
! halo_particles_<refine>
real :: box
real :: linking_parameter
real :: z_fof_ini
integer :: nh                    ! halo数量
integer :: nz = 1                    ! z_checkpoint数量
integer :: np_halo_all           ! 所有halo的粒子总数


! 每个halo的粒子数
integer :: np_halo(nh)

! 每个halo中心位置
real :: xv_mean(6, nh,1)

! 所有halo粒子的PID（连续存储）
integer(8) :: PID_halo(np_halo_all)

! 红移索引
real :: z_list = [z_fof]

! 所有halo粒子落入halo时的红移索引
integer :: z_in(np_halo_all)
```

### 3.3 数据写入逻辑

在FOF统计halo时，遍历链表同时保存各字段：

```fortran
! 分配内存
allocate(halo_particles)
allocate(halo_particles%xv_mean(6,nh,1))

np_halo_all = 0

do ihalo = 1, nh
   jp = iph_halo(ihalo)
   xp_center = xv_mean(1:3, ihalo)

   ! 初始化包围盒
   x_lo = huge(1.0); x_up = -huge(1.0)
   y_lo = huge(1.0); y_up = -huge(1.0)
   z_lo = huge(1.0); z_up = -huge(1.0)

   np = 0
   do while (jp /= 0)
      np = np + 1
      ! 存储PID到halo_particles%PID_halo数组
      i_start = np_halo_all + np
      halo_particles%PID_halo(i_start) = pid(jp)

      jp = llgp(jp)
   enddo

   halo_particles%np_halo_all = np_halo_all + np
   halo_particles%np_halo(ihalo) = np

   ! 存储halo中心位置
   halo_particles%xv_mean(:,ihalo) = xv_mean(:,ihalo)
enddo
halo_particles%z_in = current_z_index
```

> **注意**：`type_halo_particle%xv_z` 字段仅用于文件输出，用于记录每个红移时粒子的位置信息，不作为fof_time_slide程序的运行变量。

***

## 四、fof\_time\_slide.f90（新程序）

### 4.1 程序框架

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         fof_time_slide.f90                              │
├─────────────────────────────────────────────────────────────────────────┤
│  1. 初始化                                                              │
│     ├── 读取所有halo的PID和np_halo                                       │
│     │   └── total_particles = sum(np_halo)                              │
│     ├── 排序PID_tmp，生成Pt2h映射表                                      │
│     └── 预分配xv_halo(6, total_particles)连续内存空间                    │
├─────────────────────────────────────────────────────────────────────────┤
│  2. 读取当前红移粒子数据                                                  │
│     │                                                                     │
│     ├── 2.1 CUBE2 ZIPX模式: 读取rhoc解压缩为实数xp                     │
│     ├── 2.2 读取所有xp, pid, vp                                        │
│     │         初始化xv_mean(3)=xp, xv_mean(6)=[xp,vp]                             │
│     └── 2.3 遍历粒子，利用pid_in_halo找到halo并写入xv_halo              │
│                                                                          │
│  3. 释放不再需要的数据                                                   │
│     └── xp, pid, vp, PID_tmp                                            │
│                                                                          │
│  4. 对每个halo执行局部FOF                                                │
│     ├── 4.1 自适应FOF算法                                               │
│     │         N < N_mesh: 直接遍历O(N²)                                 │
│     │         N ≥ N_mesh: 网格+链表O(N)                                 │
│     ├── 4.2 可选束缚检查 (相对速度判断)                                  │
│     │         bound = 0.5*v_rel² < G*M/r                               │
│     └── 4.3 更新z_in (高红移覆盖低红移)                                  │
│                                                                          │
│  5. 输出halo_track_zin.bin                                              │
└─────────────────────────────────────────────────────────────────────────┘
```

### 4.2 提取halo粒子

通过排序和二分查找高效定位halo粒子：

```fortran
! 初始化
np_all = type_halo_particle%np_halo_all
allocate(PID_tmp(np_all), idx_to_zin(np_all), i_start_array(nhalo))
allocate(xv_halo(6, np_all), hcgp(np_all), llgp(np_all), ecgp(np_all))

! 构建i_start_array
i_start = 1
do ihalo = 1, nhalo
   i_start_array(ihalo) = i_start
   i_start = i_start + halo_particles%np_halo(ihalo)
enddo

! 生成索引数组，记录排序后PID对应的z_in索引
do i = 1, np_all
   idx_to_zin(i) = i
enddo

! 复制所有halo的PID到PID_tmp
PID_tmp = halo_particles%PID_halo

! 按PID排序，同时调整idx_to_zin
call sort_by(PID_tmp, idx_to_zin)
! 排序后：PID_tmp按升序排列，idx_to_zin记录对应的z_in索引

! 初始化链表
llgp=0; ecgp=0; hcgp=0

! 遍历所有粒子，定位并写入xv_halo
! idx_to_zin(idx)就是z_in的正确索引，也是xv_halo的正确索引
do ip = 1, np_max
   idx = pid_to_idx(pid(ip))  ! halo粒子的索引
   if (idx > 0) then
      ! 通过idx确定ihalo
      ihalo = find_halo(idx, i_start_array, nhalo)
      xv_halo(1:3, idx) = modulo(xp(1:3,ip) - halo_particles%xv_mean(1:3,ihalo,1) + L_bLb2, L_bLb*1d0)-L_bLb2 
      xv_halo(4:6, idx) = vp(ip) - halo_particles%xv_mean(4:6,ihalo,1)
      hcgp(idx) = idx
   endif
enddo

! 检查是否所有粒子都被找到（hcgp非零表示已找到）
np_found = count(hcgp /= 0)
if (np_found /= np_all) then
   print*, 'Error: particle count mismatch, found', np_found, 'expected', np_all
   call exit(1)
endif

! 释放不再需要的数据
deallocate(xp, pid, vp, PID_tmp, idx_to_zin)
```

**pid_to_idx 函数**：返回PID在PID_tmp中的索引位置

```fortran
function pid_to_idx(pid_val) result(idx)
   integer(8), intent(in) :: pid_val
   integer :: idx

   integer :: lo, mid, hi
   lo = 1; hi = np_all

   ! 二分查找
   do while (lo <= hi)
      mid = (lo + hi) / 2
      if (PID_tmp(mid) == pid_val) then
         idx = idx_to_zin(mid)
         return
      else if (PID_tmp(mid) < pid_val) then
         lo = mid + 1
      else
         hi = mid - 1
      endif
   enddo

   idx = 0  ! 未找到
end function
```

**find_halo 函数**：通过halo粒子的索引确定所属halo

```fortran
function find_halo(idx, i_start_array, nhalo) result(ihalo)
   integer, intent(in) :: idx, nhalo
   integer, intent(in) :: i_start_array(nhalo)
   integer :: ihalo

   integer :: lo, mid, hi
   lo = 1; hi = nhalo

   ! 二分查找：找到idx落在哪个halo区间
   do while (lo < hi)
      mid = (lo + hi) / 2
      if (idx >= i_start_array(mid+1)) then
         lo = mid + 1
      else
         hi = mid
      endif
   enddo
   ihalo = lo
end function
```

> **优化**：查找复杂度为O(log N)，总提取复杂度为O(np_max * log np_all)

**sort_by 函数**：快速排序，同时调整关联数组

```fortran
recursive subroutine sort_by(arr, brr, left, right)
   integer(8), intent(inout) :: arr(:)
   integer, intent(inout) :: brr(:)
   integer, intent(in) :: left, right

   integer(8) :: key_arr,tmp
   integer :: key_brr
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

```

### 4.3 局部FOF算法

#### 4.3.1 自适应策略

根据粒子数自适应选择算法：

| 条件        | 算法               | 复杂度 |
| :---------- | :----------------- | :----- |
| N < N\_mesh | 直接遍历所有粒子对 | O(N²)  |
| N ≥ N\_mesh | 构建局部网格+链表  | O(N)   |

#### 4.3.2 直接遍历（小规模halo）

```fortran
rp2 = b_link**2
istart = 0

do ihalo = 1, nhalo
   np = halo_particles%np_halo(ihalo)

   if (np < N_mesh) then
      do i = 1, np
         xp1 = xv_halo(1:3, istart+i)
         do j = i+1, np
            xp2 = xv_halo(1:3, istart+j)
            r2 = sum((xp1-xp2)**2)
            if (r2 < rp2) call merge_chain(istart+i, istart+j)
         enddo
      enddo
   endif
   istart = istart + np
enddo
```

#### 4.3.3 网格+链表（大规模halo）

```fortran
L_mesh = b_link * 2

do ihalo = 1, nhalo
   np = halo_particles%np_halo(ihalo)
   istart = i_start_array(ihalo)

   if (np >= N_mesh) then

      N_x = ceiling((maxval(xv_halo(1,istart:istart+np))-minval(xv_halo(1,istart:istart+np))) / L_mesh)
      N_y = ceiling((maxval(xv_halo(2,istart:istart+np))-minval(xv_halo(2,istart:istart+np))) / L_mesh)
      N_z = ceiling((maxval(xv_halo(3,istart:istart+np))-minval(xv_halo(3,istart:istart+np))) / L_mesh)
      xp_zero = halo_particles%xv_mean(1:3,ihalo,1) - [minval(xv_halo(1,istart:istart+np)),...,...]

      allocate(hoc(N_x,N_y,N_z), ll(np))
      hoc = 0; ll = 0

      ! 构建局部网格
      do i = 1, np
         idx = floor((xv_halo(1:3,istart+i) - xp_zero) / L_mesh) + 1
         ll(istart+i) = hoc(idx(1), idx(2), idx(3))
         hoc(idx(1), idx(2), idx(3)) = istart+i
      enddo

      ! FOF搜索 (参考fof.f90第308-344行)
      do l = 0, nlayer-1
         !$omp parallel do schedule(dynamic,1)
         do iq3 = 1+l, N_z, nlayer
            do iq2 = 1, N_y
               do iq1 = 1, N_x
                  ip = hoc(iq1, iq2, iq3)
                  do while (ip /= 0)
                     ! 检查当前粒子邻居
                     jp = ll(ip)
                     do while (jp /= 0)
                        rsq = sum((xv_halo(1:3,ip) - xv_halo(1:3,jp))**2)
                        if (rsq <= rp2) call merge_chain(ip, jp)
                        jp = ll(jp)
                     enddo
                     ! 检查周围27个格网单元
                     do i_neighbor = 1, 13
                        jq = [iq1,iq2,iq3] + ijk(:,i_neighbor)
                        if (maxval(jq) > N_x .or. minval(jq) < 1) cycle
                        jp = hoc(jq(1), jq(2), jq(3))
                        do while (jp /= 0)
                           rsq = sum((xv_halo(1:3,ip) - xv_halo(1:3,jp))**2)
                           if (rsq <= rp2) call merge_chain(ip, jp)
                           jp = ll(jp)
                        enddo
                     enddo
                     ip = ll(ip)
                  enddo
               enddo
            enddo
         enddo
         !$omp end parallel do
      enddo
   endif
enddo
```
### 4.3.4 统计halo

```fortran

do i=1,np_max
   if (hcgp(i)/=0) then
      ip = hcgp(i); jp = ip; np = 0; dxv = 0
      do while (jp /= 0)
         np = np + 1
         dxv = dxv + xv_halo(:,jp) 
         jp = llgp(jp)
      enddo
      dxv(1:3) = modulo(dxv(1:3)/np + halo_particles%xv_mean(1:3,ihalo,1), L_bLb*1d0)
      dxv(4:6) = dxv(4:6)/np
      if (np > np_halo_min) then
         np_head=np_head+1
         np_halo_team(np_head) = np
         xv_mean_team(:,np_head) =  dxv
      endif
   endif
enddo
```

#### 4.3.5 束缚检查（可选）

使用相对速度判断粒子是否被halo束缚：

```fortran
if (check_bound) then
   jp = ip
   do while (jp /= 0)
      dx = modulo((xv_halo(1:3, jp) - dxv(1:3))+ L_bLb2, L_bLb*1d0)-L_bLb2   ! 相对位置 周期性边界条件
      dv = xv_halo(4:6, jp) - dxv(4:6)  ! 相对速度
      r = norm2(dx)
      v_rel_sq = sum(dv**2)
      ! 束缚条件: 相对动能 < 引力势能
      if (0.5 * v_rel_sq < np / r) then
         halo_particles%z_in(jp) = z_check
      endif
      jp = llgp(jp)
   enddo
else
   halo_particles%z_in(istart:istart+np-1) = z_check
endif
```

***

## 五、大规模模拟的多节点优化

### 5.1 并行策略

```
模拟总节点数 = nn³
每个计算节点负责 = nn³ / n_compute_nodes 个模拟节点的数据
```

### 5.2 两阶段处理

**阶段一：粗筛**

- 各计算节点读取负责区域的粒子数据
- 根据halo的xv判断是否有粒子在当前节点
- 汇总需要处理的halo列表

**阶段二：精确定位**

- 每个节点参考fof.f90加载带buffer的所有粒子
- 每个节点只保留负责区域内的halo粒子
- 分布式执行局部FOF

### 5.3 通信优化

- 对边界halo使用边界扩展减少通信次数

***

## 六、设计决策汇总

| 问题     | 决策                                              |
| :------- | :------------------------------------------------ |
| 数据结构 | 统一使用type\_halo\_particle                      |
| 束缚检查 | 可选开关                                          |
| FOF算法  | 自适应: N < N\_mesh直接遍历, N ≥ N\_mesh网格+链表 |
| 束缚判断 | 相对速度: 0.5*v\_rel² < G*M/r                     |
| 输出文件 | halo\_track\_zin.bin                              |
| 红移遍历 | 从低到高，高红移覆盖低红移                        |
| 粒子匹配 | 按pid精确匹配                                     |
| 位置单位 | 以L\_ng为单位，ratio\_refine转换                  |

