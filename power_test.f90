! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !   CUBE™ in Coarray Fortran  !
! !   haoran@xmu.edu.cn         !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! program CUBE
!    use omp_lib
!    use variables
!    use cicpower_segment
!    implicit none
!    save
!   integer(4) c1,c2,cn,cpk

!   c1=1
!   c2=9
!   cn=1

!    do cpk=1,2
!       do cc=c1,c2,cn
!         sim%cur_checkpoint = cc
!         if (head) print*, ''
!         if (head) print*, ''
!         if (head) print*, ''
!         if (head) print*, ''
!         if (head) print*, ''
!         if (head) print*, '---------- get_power ----------',z_checkpoint(sim%cur_checkpoint)
!             print*,cc
!             call initialize
!             call particle_initialization
!             call buffer_grid
!             call buffer_x
!             call buffer_v
!             call timestep
!             call drift
!             call buffer_grid
!             call buffer_x
!             sim%calculate_PK = cpk
!             sim%cur_powerpoint = cc+2
!             call get_power
!         if (head) close(77)
!         call finalize
!       enddo
!    enddo
!    print*,''
!    print*,''
!    print*,''
!    print*,''
!    print*,'+++++++++  power times  ++++++++++'
!    print*,' --------------------------------'
!    print*,'   check       |  corse       |     std.   |    fine'
!    print*,' --------------------------------'
!    do cc=c1,c2,cn
!       print*,cc,'|',tps(3,cc),'|',tps(4,cc),'|',tps(5,cc)
!       print*,' --------------------------------'
!    enddo
!    print*,''
!    print*,''
!    print*,''
!    print*,''
!    print*,'+++++++++  power times  ++++++++++'
!    print*,' --------------------------------'
!    print*,'  segment      |    global  |    rate'
!    print*,' --------------------------------'
!    do cc=c1,c2,cn
!       print*,tps(1,cc),'|',tps(2,cc),'|',tps(2,cc)/tps(1,cc)
!       print*,' --------------------------------'
!    enddo
! contains

!    function Dgrow(scale_factor)
!       implicit none
!       real, parameter :: om=omega_m
!       real, parameter :: ol=omega_l
!       real scale_factor
!       real Dgrow
!       real g,ga,hsq,oma,ola
!       hsq=om/scale_factor**3+(1-om-ol)/scale_factor**2+ol
!       oma=om/(scale_factor**3*hsq)
!       ola=ol/hsq
!       g=2.5*om/(om**(4./7)-ol+(1+om/2)*(1+ol/70))
!       ga=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))
!       Dgrow=scale_factor*ga/g
!    endfunction Dgrow

! endprogram
