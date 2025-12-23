! #define power_test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CUBE™ in Coarray Fortran  !
!   haoran@xmu.edu.cn         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program CUBE
   use omp_lib
use variables
   use power_nu
   implicit none
   save

#ifdef power_test
   integer(4) c1,c2,cn,cpk

   c1=8
   c2=8
   cn=1

   do cpk=1,1
      do cc=c1,c2,cn
         sim%cur_checkpoint = cc
         if (this_image()==1) print*, ''
         if (this_image()==1) print*, ''
         if (this_image()==1) print*, ''
         if (this_image()==1) print*, ''
         if (this_image()==1) print*, ''
         if (this_image()==1) print*, '---------- get_power ----------',cc,cpk
         call initialize
         call particle_initialization
         call buffer_grid
         call buffer_x
         call buffer_v
         call timestep
         call drift
         call buffer_grid
         call buffer_x
         sim%calculate_PK = cpk
         sim%cur_powerpoint = 2
         do while(1/(1+z_powerpoint(sim%cur_powerpoint)) < sim%a .and. sim%cur_powerpoint<=n_powerpoint)
            sim%cur_powerpoint = sim%cur_powerpoint+1
            ! print*, '  power_nu',cc,cpk,sim%cur_powerpoint,z_powerpoint(sim%cur_powerpoint),1/(1+z_powerpoint(sim%cur_powerpoint)),sim%a
         enddo
         call get_power
         call finalize
      enddo
   enddo
   if (head) then
      print*,''
      print*,''
      print*,''
      print*,''
      print*,'+++++++++  power times  ++++++++++'
      print*,' --------------------------------'
      print*,'   check       |  corse       |     std.   |    fine'
      print*,' --------------------------------'
      do cc=c1,c2,cn
         print*,cc,'|',tps(3,cc),'|',tps(4,cc),'|',tps(5,cc)
         print*,' --------------------------------'
      enddo
      print*,''
      print*,''
      print*,''
      print*,''
      print*,'+++++++++  power times  ++++++++++'
      print*,' --------------------------------'
      print*,'  segment      |    global  |    rate'
      print*,' --------------------------------'
      do cc=c1,c2,cn
         print*,tps(1,cc),'|',tps(2,cc),'|',tps(2,cc)/tps(1,cc)
         print*,' --------------------------------'
      enddo
   endif
#else
   sim%cur_checkpoint=1
   call initialize
   call particle_initialization
   call buffer_grid
   call buffer_x
   call buffer_v
   if (sim%cur_powerpoint == 1 .and. Mass_nu > 0 ) then
      power_step = .true.
      call interp_Pk_CDM
      sim%cur_powerpoint = 2
   endif
   sim%cur_checkpoint = sim%cur_checkpoint+1
   sim%cur_halofind = sim%cur_halofind+1
   ! print*,'checkpoint ',sim%cur_checkpoint
   if (head) open(77,file=output_dir()//'vinfo'//output_suffix(),access='stream',status='replace')

   if (head) print*, '---------- starting main loop ----------'
   do istep=sim%timestep,istep_max
   ! do istep=sim%timestep,sim%timestep+3
      call system_clock(t_start,t_rate)
      call tic(100)
      call timestep
      call drift
      call buffer_grid
      call buffer_x
      call kick
      call buffer_v
      if (checkpoint_step .or. halofind_step) then
         dt_old=0
         call drift
         if (checkpoint_step) then
            call checkpoint
            sim%cur_checkpoint = sim%cur_checkpoint+1
            if (sim%cur_checkpoint .eq. 73) then
               final_step=.true.
               print*, '  final checkpoint'
            endif
         endif
         call buffer_grid
         call buffer_x
         call buffer_v
         if (halofind_step) then
            !call halofind_FoF
            sim%cur_halofind = sim%cur_halofind+1
         endif
         call print_header(sim)
         if (final_step) exit
         dt=0
      endif
      call system_clock(t_end,t_rate)
      call toc(100)
      if(head) print*, 'total elapsed time =',tcat(100,istep),real(t_end-t_start)/t_rate,'secs';
   enddo
   if (head) close(77)
   call finalize
#endif

contains

   function Dgrow(scale_factor)
      implicit none
      real, parameter :: om=omega_m
      real, parameter :: ol=omega_l
      real scale_factor
      real Dgrow
      real g,ga,hsq,oma,ola
      hsq=om/scale_factor**3+(1-om-ol)/scale_factor**2+ol
      oma=om/(scale_factor**3*hsq)
      ola=ol/hsq
      g=2.5*om/(om**(4./7)-ol+(1+om/2)*(1+ol/70))
      ga=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))
      Dgrow=scale_factor*ga/g
   endfunction Dgrow

endprogram

