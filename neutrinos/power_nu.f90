module power_nu
   use cicpower_segment
   use parameters
   use ieee_arithmetic
   implicit none

   ! real k_tf1(nw_global/2+1,nw,npen),k_tf2(ngt/2+1,ngt,ngt),k_tf3_2(nft(2)/2+1,nft(2),nft(2)),k_tf3_4(nft(3)/2+1,nft(3),nft(3)),k_tf3_6(nft(4)/2+1,nft(4),nft(4)),k_tf3_8(nft(5)/2+1,nft(5),nft(5))
   ! real X_temp(npbin)
   ! integer(4) t_nu

contains

   real function expansion(a0,dt0)
      use variables
      ! use parameters
      implicit none
      real(8) :: a_x,adot,t_x,tdoa,a8_0
      real(4) :: a0,dt0
      integer i1,i2

      a8_0=a0
      i1 = 1
      do while(s2a(i1+1)>a8_0 .and. i1 < istep_max)
         i1 = i1+1
      enddo
      ! print*,'a i1',i1,dt0

      tdoa = (stime(i1+1)-stime(i1))/(s2a(i1+1)-s2a(i1))
      t_x = stime(i1)+tdoa*(a8_0-s2a(i1))+dt0
      ! print*,'a t_x',t_x,tdoa


      i2 = i1
      do while(stime(i2+1)<t_x .and. i2 > 1)
         i2 = i2-1
      enddo
      ! print*,'a i2',i2

      adot = (s2a(i2+1)-s2a(i2))/(stime(i2+1)-stime(i2))
      a_x = s2a(i2)+adot*(t_x-stime(i2))
      ! print*,'    i1',i1,stime(i1),s2a(i1)
      ! print*,'    i2',i2,stime(i2),s2a(i2)
      ! print*,a_x,a0,a8_0,t_x

      expansion=a_x-a8_0
      ! da1=(a0-a_x)/2
      ! da2=(a0-a_x)/2
   endfunction

   real function H_a(a0)
      use variables
      implicit none
      real a0,Hdot,t_x,tdoa
      integer i1,i2

      i1 = 1
      do while(s2a(i1+1)>a0 .and. i1 < istep_max)
         i1 = i1+1
      enddo
      i2=i1+1
      tdoa = (stime(i2)-stime(i1))/(s2a(i2)-s2a(i1))
      t_x = stime(i1)+tdoa*(a0-s2a(i1))

      Hdot = (s2H(i2)-s2H(i1))/(stime(i2)-stime(i1))
      H_a = s2H(i1)+Hdot*(t_x-stime(i1))
      ! print*,H_a,'H    i1',i1,stime(i1),s2H(i1)
      ! print*,t_x,'H    i2',i2,s2a(i2),s2H(i2),a0-s2a(i1)
   endfunction

   real function dtau_a(da0)
      use variables
      implicit none
      real a0,da0,taudot,t_x,tdoa
      integer i1,i2
      if (abs((sim%a+da0)-1)<1e-6) then
         dtau_a = s2tau(1)-sim%tau
         print*,'tau0'
      else
         a0 = sim%a+da0
         i1 = 1
         do while(s2a(i1+1)>a0 .and. i1 < istep_max)
            i1 = i1+1
         enddo
         i2=i1+1
         tdoa = (stime(i2)-stime(i1))/(s2a(i2)-s2a(i1))
         t_x = stime(i1)+tdoa*(a0-s2a(i1))

         taudot = (s2tau(i2)-s2tau(i1))/(stime(i2)-stime(i1))
         dtau_a = s2tau(i1)+taudot*(t_x-stime(i1))-sim%tau
      endif
   endfunction

   real function sf_a(a0)
      use variables
      implicit none
      real a0,tdoa
      integer i1,i2

      i1 = 1
      do while(s2a(i1+1)>a0 .and. i1 < istep_max)
         i1 = i1+1
      enddo
      i2=i1+1
      tdoa = (stime(i2)-stime(i1))/(s2a(i2)-s2a(i1))
      sf_a = (stime(i1)+tdoa*(a0-s2a(i1)))*omhsq0-s_fi
      ! sf_a = (stime(i1)+tdoa*(a0-s2a(i1)))*omhsq0-s_fi
      ! print*,sf_a,'sf    i1',i1,stime(i1),s2a(i1),a0

   endfunction

   subroutine get_f_nr(a,f_nr_a)
      use omp_lib
      use parameters
      implicit none

      real a, f_nr_a(3)
      real(8) :: x(3),y(3), Fy(3) , dFy(3)
      real(8) ::  k1(3), k2(3), k3(3), k4(3), h
      integer ix,n

      ! print*,'para',Mass_nu,a,sigma_nu,N_eff,k_b,T_gama

      if (a_nu == 0) then
         y = m_nu*a/(sigma_nu*N_eff*k_b*T_gama)
         ! print*,'y',y

         ! 定义步长和循环变量
         h = 1./2
         n = 300./h

         ! 初始化 Fy
         Fy = 0.0

         ! RK4 循环
         do ix=0,n
            x =  h * ix
            k1 = h * (x**2*sqrt(x**2+y**2))/(1+exp(x))
            k3 = h * ((x+h/2)**2*sqrt((x+h/2)**2+y**2))/(1+exp((x+h/2)))
            k4 = h * ((x+h)**2*sqrt((x+h)**2+y**2))/(1+exp((x+h)))
            Fy = Fy + (k1 + 4*k3 + k4)/6
         enddo


         ! 初始化 dFy
         dFy = 0.0

         ! RK4 循环
         do ix=0,n
            x =  h * ix
            k1 = h * (x**2)/((1+exp(x))*sqrt(x**2+y**2))
            k3 = h * ((x+h/2)**2)/((1+exp((x+h/2)))*sqrt((x+h/2)**2+y**2))
            k4 = h * ((x+h)**2)/((1+exp((x+h)))*sqrt((x+h)**2+y**2))
            dFy = dFy + (k1 + 4*k3 + k4)/6
         enddo

         f_nr_a = y**2*dFy/Fy
      elseif (a > a_nu) then
         f_nr_a=1
      else
         f_nr_a=0
      endif
   endsubroutine

   subroutine interp_quad(a, cur_powerpoint, Pk, z_point, result)
      use variables, only: ifs
      use parameters
      implicit none
      integer(8) cur_powerpoint
      real, intent(in) :: a, Pk(npbin, nmax_redshift), z_point(nmax_redshift)
      real, intent(out) :: result(npbin)
      real a1(npbin), a2(npbin), k2(npbin), b2(npbin), c2(npbin)
      real log_a, log_a1, log_a2

      log_a = log(a)
      log_a1 = log(1 / (z_point(1) + 1))
      log_a2 = log(1 / (z_point(2) + 1))

      a1 = (log(Pk(:, 1)) - log(Pk(:, 2))) / (log_a1 - log_a2)
      b2 = log(Pk(:, 1)) - a1 * log_a1
      result = exp(a1 * log_a + b2)

      if ((minval(result(:ifs)) < -1e-5) .or. any((ieee_is_nan(result(:ifs))))) then
         if (head) then
            print *, ''
            print *, ''
            print *, ''
            print *, 'Pk interp err ', minval(result), 1. / a - 1, ifs
            print *, 'cur_powerpoint', cur_powerpoint, z_point(cur_powerpoint)
            print *, 'Pk', result
            print *, 'a1', a1
            print *, 'b', b2
            print *, 'Pk1',Pk(:, 1)
            print *, 'Pk2',Pk(:, 2)
         endif
         error stop 'Pk interp err '
      endif
   endsubroutine

   subroutine exterp_quad(a, cur_powerpoint, Pk, z_point, result)
      use variables
      use parameters
      implicit none
      integer(8) cur_powerpoint
      real, intent(in) :: a, Pk(npbin, nmax_redshift), z_point(nmax_redshift)
      real, intent(out) :: result(npbin)
      real a1(npbin), a2(npbin), k2(npbin), b2(npbin), c2(npbin)
      real log_a, log_ap, log_ap1, log_ap2

      log_a = log(a)
      log_ap = log(1 / (z_point(cur_powerpoint) + 1))
      log_ap1 = log(1 / (z_point(cur_powerpoint - 1) + 1))
      log_ap2 = log(1 / (z_point(cur_powerpoint - 2) + 1))

      a1 = (log(Pk(:, cur_powerpoint)) - log(Pk(:, cur_powerpoint - 1))) / (log_ap - log_ap1)
      a2 = (log(Pk(:, cur_powerpoint - 1)) - log(Pk(:, cur_powerpoint - 2))) / (log_ap1 - log_ap2)
      k2 = (a1 - a2) / (log_ap - log_ap2)
      b2 = a1 - k2 * (log_ap + log_ap1)
      c2 = log(Pk(:, cur_powerpoint)) - k2 * log_ap ** 2 - b2 * log_ap
      result = exp(k2 * log_a ** 2 + b2 * log_a + c2)

      if ((minval(result(:ifs)) < -1e-5) .or. any((ieee_is_nan(result(:ifs))))) then
         if (head) then
            print *, ''
            print *, ''
            print *, ''
            print *, 'Pk exterp err ', minval(result(:ifs)), 1. / a - 1, ifs
            print *, 'cur_powerpoint', cur_powerpoint, z_point(cur_powerpoint)
            print *, 'Pk0', exp(Pk(:, cur_powerpoint))
            print *, 'Pk-1', exp(Pk(:, cur_powerpoint - 1))
            print *, 'Pk-2', exp(Pk(:, cur_powerpoint - 2))
            print *, 'a1', a1
            print *, 'a2', a2
            print *, 'k', k2
            print *, 'b', b2
            print *, 'c', c2
            print *, 'result', result(:ifs)
         endif
         error stop 'Pk exterp err '
      endif
   endsubroutine

   real function interp_tf_F(k_lin_log,k_need)
      use variables
      implicit none
      real, intent(in) :: k_lin_log(npbin),k_need
      real a1,a2,k2,b2,c2
      integer(4) i_mid,i1,i2,i3
      real k_need_log

      k_need_log = log(k_need)

      i1=1; i2=npbin
      do while (i2-i1>1)
         i_mid=(i1+i2)/2
         if (k_need_log>k_lin_log(i_mid)) then
            i1=i_mid
         else
            i2=i_mid
         endif
      enddo

      if (k_need_log-k_lin_log(i1) < k_lin_log(i2)-k_need_log) then
         i3 = i2
         i2 = i1
         i1 = i1-1
      else
         i3 = i2+1
      endif

      if (i3 > npbin) then
         i1 = i1-1
         i2 = i2-1
         i3 = i3-1
      endif

      a1 = (tf_F_log(i1)-tf_F_log(i2))/(k_lin_log(i1)-k_lin_log(i2))
      a2 = (tf_F_log(i2)-tf_F_log(i3))/(k_lin_log(i2)-k_lin_log(i3))
      k2 = (a1-a2)/(k_lin_log(i1)-k_lin_log(i3))
      b2 = a1-k2*(k_lin_log(i1)+k_lin_log(i2))
      c2 = tf_F_log(i1)-k2*k_lin_log(i1)**2-b2*k_lin_log(i1)
      interp_tf_F = exp(k2*k_need_log**2+b2*k_need_log+c2)

      ! print*,'interp:++++++++++++++'
      ! print*,i1,i2,i3
      ! print*,kh_lin(i1),kh_lin(i2),kh_lin(i3),k_need
      ! print*,tf_F(i1),tf_F(i2),tf_F(i3),interp_tf_F
      ! print*,'interp,done'

   endfunction

   subroutine tf_F_correction
      use variables
      implicit none
      real(8) a1,b1,d1,d2
      integer(4) i1,i2

      i1 = ifs-1
      d1 = (tf_F(i1+1)-tf_F(i1))/(kh_lin(i1+1)-kh_lin(i1))
      do while (i1 > 1 .and. d1 < -1)
         i1 = i1-1
         d1 = (tf_F(i1+1)-tf_F(i1))/(kh_lin(i1+1)-kh_lin(i1))
      enddo

      if (head) print*,'  tf_F_correction',i1,kh_lin(i1),k_fs,tf_F(i1)

      a1 = ((1-f_nu)-tf_F(i1))/(10-kh_lin(i1))
      b1 = tf_F(i1)-kh_lin(i1)*a1

      i2 = i1
      do while (i2 <= npbin .and. kh_lin(i2)<10)
         tf_F(i2) = a1*kh_lin(i2)+b1
         i2 = i2+1
      enddo
      if (i2<npbin) tf_F(i2:)=1-f_nu

   endsubroutine

   subroutine get_L(sf,k_lin_2d,Pk,Lx)
      use variables, only: ifs,m_nu_2d
      use parameters
      implicit none
      real, intent(in) ::  k_lin_2d(npbin,3),Pk(npbin),sf
      real(8), intent(out) ::  Lx(npbin,3)
      real(8) X(npbin,3)
      real(8) sq_Pk_2d(npbin,3)

      sq_Pk_2d = spread(sqrt(Pk), dim=2, ncopies=3)  ! reshape Pk to a 2D array
      X = sf*k_lin_2d/(m_nu_2d)*T_nu0*C
      Lx = ((1+0.0168*X**2+0.0407*X**4)/(1+2.1734*X**2+1.6787*X**4.1811+0.1467*X**8))*sq_Pk_2d*sf


      ! print*,'  Pk',Pk(1),Pk(npbin/3),Pk(npbin/3*2),Pk(npbin)
      ! print*,'  sq_Pk_2d1',sq_Pk_2d(1,1),sq_Pk_2d(npbin/3,1),sq_Pk_2d(npbin/3*2,1),sq_Pk_2d(npbin,1)
      ! print*,'  sq_Pk_2d2',sq_Pk_2d(1,2),sq_Pk_2d(npbin/3,2),sq_Pk_2d(npbin/3*2,2),sq_Pk_2d(npbin,2)
      ! print*,'  sq_Pk_2d3',sq_Pk_2d(1,3),sq_Pk_2d(npbin/3,3),sq_Pk_2d(npbin/3*2,3),sq_Pk_2d(npbin,3)
      ! X=(k_lin_2d/m_nu_2d)
      ! print*,'  (k_lin_2d/m_nu_2d)',X(1,1),X(npbin/3,1),X(npbin/3*2,1),X(npbin,1)

      ! print*,'  X1',X(1,1),X(npbin/3,1),X(npbin/3*2,1),X(npbin,1)
      ! print*,'  X2',X(1,2),X(npbin/3,2),X(npbin/3*2,2),X(npbin,2)
      ! print*,'  X3',X(1,3),X(npbin/3,3),X(npbin/3*2,3),X(npbin,3)
      ! print*,'  Lx1',Lx(1,1),Lx(npbin/3,1),Lx(npbin/3*2,1),Lx(npbin,1)
      ! print*,'  Lx2',Lx(1,2),Lx(npbin/3,2),Lx(npbin/3*2,2),Lx(npbin,2)
      ! print*,'  Lx3',Lx(1,3),Lx(npbin/3,3),Lx(npbin/3*2,3),Lx(npbin,3)
      if ((minval(Lx(:ifs,:)) < -1e-5) .or. any((ieee_is_nan(Lx(:ifs,:)))))then
         if(head) then
            print*,''
            print*,''
            print*,''
            print*,'get_L err ',minval(Lx),sf
            print*,'Pk is nan ',any((ieee_is_nan(Pk)))
            print*,'X',X(npbin-40:npbin,:)
            X = ((1+0.0168*X**2+0.0407*X**4)/(1+2.1734*X**2+1.6787*X**4.1811+0.1467*X**8))
            print*,''
            print*,''
            print*,''
            print*,' XX',X(npbin-40:npbin,:)
            print*,''
            print*,''
            print*,''
            print*,' (Pk)',Pk(:ifs)
            print*,' (sqrt(Pk))',sqrt(Pk(:ifs))
            print*,'Lx',Lx(ifs-40:ifs,:)
         endif
         error stop 'get_L err '
      endif
      ! stop
   endsubroutine

   subroutine get_sqrt_pk_nu1(sf,a0,Lx)
      use variables, only: ifs,sq_Pk_nu_ic,Pk_nu_ic,kh_lin_2d,m_nu_2d,H_i,H_0
      use parameters
      implicit none
      real, intent(in) ::  sf,a0
      real, intent(out) ::  Lx(npbin,3)
      real(8) X(npbin,3)
      real f,omega_m_zi,omega_l_zi

      omega_m_zi=(H_0/H_i)**2*sim%omega_m*(1/a0)**3
      omega_l_zi=(H_0)**2/H_i**2*sim%omega_l
      f=omega_m_zi**(4.0/7.0)+omega_l_zi/70.0*(1+omega_m_zi/2.0)
      X = sf*kh_lin_2d/(m_nu_2d)*T_nu0*C
      Lx = ((1+0.0168*X**2+0.0407*X**4)/(1+2.1734*X**2+1.6787*X**4.1811+0.1467*X**8))*sq_Pk_nu_ic*(1+sf*a0**2*H_i*f)

      ! print*,'  Lx0',Lx0(1,1),Lx0(npbin/3,1),Lx0(npbin/3*2,1),Lx0(npbin,1)

      ! print*,'  sq_Pk_nu_ic',sq_Pk_nu_ic(1,2),sq_Pk_nu_ic(npbin/3,2),sq_Pk_nu_ic(npbin/3*2,2),sq_Pk_nu_ic(npbin,2)
      ! print*,'  X',X(1,1),X(npbin/3,1),X(npbin/3*2,1),X(npbin,1)
      if ((minval(Lx(:ifs,:)) < -1e-5) .or. any((ieee_is_nan(Lx(:ifs,:)))))then
         if(head) then
            print*,''
            print*,''
            print*,''
            print*,'get_Lic err ',minval(Lx),sf
            print*,'Pk_ic is nan ',any((ieee_is_nan(Pk_nu_ic)))
            print*,'X',X(npbin-40:npbin,:)
            X = ((1+0.0168*X**2+0.0407*X**4)/(1+2.1734*X**2+1.6787*X**4.1811+0.1467*X**8))
            print*,''
            print*,''
            print*,''
            print*,' XX',X(npbin-40:npbin,:)
            print*,''
            print*,''
            print*,''
            print*,' (Pk)',Pk_nu_ic(npbin-40:npbin,:)
            print*,' (sqrt(Pk))',sq_Pk_nu_ic(npbin-40:npbin,:)
            print*,'Lx',Lx(npbin-40:npbin,:)
         endif
         error stop 'get_Lic err '
      endif
   endsubroutine

   subroutine interp_line_pk2(k_lin_2d,Pk1,Pk2,s1,s2,tau1,tau2,a1,a2,n,Pk_out)
      use variables, only: ifs
      use parameters
      implicit none
      real, intent(in) ::  k_lin_2d(npbin,3),Pk1(npbin),Pk2(npbin),s1,s2,tau1,tau2,a1,a2
      integer, intent(in) ::  n
      real, intent(inout) ::  Pk_out(npbin,3)
      real(8) L0(npbin,3),L1(npbin,3)
      real Pk_a(npbin),Pk_b(npbin),Pk_now(npbin)
      real ai,si,s_b,a_i,s_a,tau_i
      integer ipk

      ! print*,s1,s2,tau1,tau2,a1,a2,n
      if (abs(a1/a2-1)>1e-7) then
         Pk_a = (Pk1-Pk2)/(a1-a2)
         Pk_b = Pk1 - Pk_a*a1
         s_a = (s1-s2)/(1/a1-1/a2)
         s_b = s1 - s_a*(1/a1-1)
         a_i = (a2-a1)/n
         tau_i = (tau2-tau1)/n
         ! print*,s_a,s_b,a_i,tau_i

         call get_L(s1,k_lin_2d,Pk1,L1)
         do ipk=1,n
            ! print*,'       ipk',ipk
            ai = (a1+a_i*ipk)
            si = s_a*(1/ai-1)+s_b
            Pk_now = Pk_a*ai+Pk_b
            if (si<0) then
               if (si<-1e-4) then
                  if (head) print*,'si           =',si,s1,s2
                  if (head) print*,'sa/sb        =',s_a,s_b
                  if (head) print*,'a_i/tau_i/ai =',a_i,tau_i,ai
                  error stop 'si is nagtive in interp_line_pk2'
               endif
               si = 0
            endif
            if (isnan(si)) then
               if (head) print*,'si           =',si
               if (head) print*,'sa/sb        =',s_a,s_b
               if (head) print*,'a_i/tau_i/ai =',a_i,tau_i,ai
               if (head) print*,'a1/a2        =',a1,a2,a1 < a2,a2/a1-1
               error stop 'si is nan in interp_line_pk2'
            endif
            if (minval(Pk_now(:ifs)) < -1e-5 .or. any((ieee_is_nan(Pk_now(:ifs)))))then
               if(head) then
                  print*,''
                  print*,''
                  print*,''
                  print*,'interp_line_pk2 Pk_interp err ',minval(Pk_now),ai
                  print*,'Pk is nan ',any((ieee_is_nan(Pk_now)))
                  print*,' a,tau',a2,a1,tau2,tau1
                  print*,' (Pk)',Pk_now
               endif
               error stop 'interp_line_pk2 Pk_interp err '
            endif

            L0 = L1
            call get_L(si,k_lin_2d,Pk_now,L1)
            Pk_out = Pk_out+(L0+L1)*tau_i
         enddo
      endif
   endsubroutine

   subroutine interp_Pk_CDM
      use omp_lib
      use variables
      use parameters
      implicit none

      if (power_step) then
         if (head) print*,''
         if (head) print*,'Power_step'
         if (calculate_PK > 0)then
            call get_power
            Pk_cb_check(:,sim%cur_powerpoint) = xi_cdm(5,1:)[1]

            ! if ( sim%cur_powerpoint > 3 .and. z_powerpoint(sim%cur_powerpoint)>=0) then
            !    write(str_z,'(f8.4)') z_powerpoint(sim%cur_powerpoint)
            !    open(11,file=nupath//'Pk_cb_'//trim(adjustl(str_z))//'.txt',form='formatted')
            !    read(11,*) Pk_cb_check(:,sim%cur_powerpoint)
            !    close(11)
            ! endif
         else
            if ( sim%cur_powerpoint > 3 .and. z_powerpoint(sim%cur_powerpoint)>=0) then
               write(str_z,'(f8.4)') z_powerpoint(sim%cur_powerpoint)
               open(11,file=nupath//'Pk_cb_'//trim(adjustl(str_z))//'.txt',form='formatted')
               read(11,*) Pk_cb_check(:,sim%cur_powerpoint)
               close(11)
               if (calculate_PK == -1) then
                  open(11,file=nupath//'Pk_nu_'//trim(adjustl(str_z))//'.txt',form='formatted')
                  read(11,*) Pk_nu_check(:,sim%cur_powerpoint)
                  close(11)
               endif
            endif
            if (head) print*,'      Pk_check',sim%cur_powerpoint,z_powerpoint(sim%cur_powerpoint),str_z
            print*,''
         endif
         sim%cur_powerpoint=sim%cur_powerpoint+1
      endif

      if (sim%cur_powerpoint <= 3) then
         if (head) print*,'      interp',1/a_step(istep)-1,sim%cur_powerpoint,istep,a_step(istep),sim%a
         call interp_quad(a_step(istep),sim%cur_powerpoint-1,Pk_cb_check,z_powerpoint,Pk_step(:,istep))
         if (calculate_PK == -1) then
            call interp_quad(a_step(istep),sim%cur_powerpoint-1,Pk_nu_check,z_powerpoint,Pk_nus(:,1))
            call interp_quad(a_step(istep),sim%cur_powerpoint-1,Pk_nu_check,z_powerpoint,Pk_nus(:,2))
            call interp_quad(a_step(istep),sim%cur_powerpoint-1,Pk_nu_check,z_powerpoint,Pk_nus(:,3))
         endif
      else
         if (head) print*,'      exterp',1/a_step(istep)-1,sim%cur_powerpoint,istep
         call exterp_quad(a_step(istep),sim%cur_powerpoint-1,Pk_cb_check,z_powerpoint,Pk_step(:,istep))
         if (calculate_PK == -1) then
            call exterp_quad(a_step(istep),sim%cur_powerpoint-1,Pk_nu_check,z_powerpoint,Pk_nus(:,1))
            call exterp_quad(a_step(istep),sim%cur_powerpoint-1,Pk_nu_check,z_powerpoint,Pk_nus(:,2))
            call exterp_quad(a_step(istep),sim%cur_powerpoint-1,Pk_nu_check,z_powerpoint,Pk_nus(:,3))
         endif
      endif
   endsubroutine

   subroutine get_tf_cb2matter
      use variables
      use parameters
      implicit none

      integer is,n
      real s_pre,s_post
      real sqrt_pk_nu1(npbin,3),sqrt_pk_nu2(npbin,3),sqrt_pk_nu2_step(npbin,3)
      real(8) Lx0(npbin,3),Lx1(npbin,3),Lx2(npbin,3)

      ifs = npbin
      k_fs = max((0.08/sqrt(1/a_step(istep)) * sqrt(sim%omega_m/0.3) * Mass_nu/3/0.1 * h0 )*tf_smooth,10)
      do while (kh_lin(ifs) > k_fs)
         ifs = ifs-1
      enddo
      call get_f_nr(a_step(istep),f_nr)
      s_f = sf_a(a_step(istep))
      sf_step(istep) = s_f
      ! if(head) print*,'      ifs = ',ifs

      ! print*,'istep',istep
      ! print*,Pk_step(:,istep)

      if(head) then
         print*,''
         print*,'Get Tf_nu'
         print*,'  save Pk in z = ',1/a_step(istep)-1,istep
         print*,'  s_f    = ',s_f
         print*,'  1-f_nu = ',(1-f_nu)
         print*,'  f_nus  = ',f_nus
         print*,'  f_nr   = ',f_nr
         print*,'  ifs    = ',ifs,npbin
         print*,'  k_fs   = ',k_fs
      endif


      call tic(96)
      call interp_Pk_CDM
      sync all; call toc(96)

      call tic(97)
      if (istep >0) then
         if(head) then
            if (calculate_PK > -1) then
               call get_sqrt_pk_nu1(s_f,a_step(0),sqrt_pk_nu1)

               s_pre = 0
               s_post = 0
               sqrt_pk_nu2 = 0
               n = 2000/istep+1
               do is=1,istep
                  ! print*,'      is=',is
                  ! print*,Pk_step(:,is-1)
                  ! print*,'=================================='
                  ! print*,Pk_step(:,is)
                  call interp_line_pk2(kh_lin_2d,Pk_step(:,is-1),Pk_step(:,is),(s_f-sf_step(is-1)),(s_f-sf_step(is)),tau_step(is-1),tau_step(is),a_step(is-1),a_step(is),n,sqrt_pk_nu2)
               enddo
               sqrt_pk_nu2 = 1/2.*1.5*sim%omega_m*H_0**2*sqrt_pk_nu2
               ! print*,'  sqrt_pk_nu2',sqrt_pk_nu2(1,2),sqrt_pk_nu2(npbin/3,2),sqrt_pk_nu2(npbin/3*2,2),sqrt_pk_nu1(npbin,2)
               Pk_nus = (sqrt_pk_nu1+sqrt_pk_nu2)

               if ((minval(Pk_nu(:ifs)) < -1e-5) .or. any((ieee_is_nan(Pk_nu(:ifs)))))then
                  if(head) then
                     print*,''
                     print*,''
                     print*,''
                     print*,'Pk_nu err ',minval(Pk_nu),1/a_step(istep)-1
                     print*,'Pk_check(-1 0) is nan ',any((ieee_is_nan(Pk_cb_check(:,sim%cur_powerpoint)))),any((ieee_is_nan(Pk_cb_check(:,sim%cur_powerpoint-1))))
                     print*,'Pk_step(-1 0) is nan ',any((ieee_is_nan(Pk_step(:,istep-1)))),any((ieee_is_nan(Pk_step(:,istep))))
                     print*,' a,tau',a_step(istep),a_step(istep-1),tau_step(istep),tau_step(istep-1)
                     print*,'Pk_m -1',Pk_step(:,istep-1)
                     print*,''
                     print*,''
                     print*,''
                     print*,'Pk_m 0',Pk_step(:,istep)
                     print*,''
                     print*,''
                     print*,''
                     print*,'Pk_nu',Pk_nu(:ifs)
                  endif
                  error stop 'Pk_nu err '
               endif
            endif
         endif
      else
         Pk_nus=sq_Pk_nu_ic
      endif

      ! print*,'  Pk_nu',(f_nus(1)*f_nr(1)*1+f_nus(2)*f_nr(2)*1+f_nus(3)*f_nr(3)*1)**2
      Pk_nu = (f_nus(1)*f_nr(1)*Pk_nus(:,1)+f_nus(2)*f_nr(2)*Pk_nus(:,2)+f_nus(3)*f_nr(3)*Pk_nus(:,3))**2
      Pk_nus = Pk_nus**2
      tf_F = ((1-f_nu)*sqrt(abs(Pk_step(:,istep)))+f_nu*sqrt(abs(Pk_nu)))/sqrt(abs(Pk_step(:,istep)))
      if (calculate_PK > -1) call tf_F_correction

      if (a_nu == 0 .or. a_step(istep) > a_nu) then
         Pk_step(:,istep) = tf_F**2*Pk_step(:,istep) !Pk_step_cdm --> Pk_step_matter
         write(str_z,'(f8.4)') 1/a_step(istep)-1
         if(head) then
            print*,'  Pk_step',Pk_step(1,istep),Pk_step(npbin/3,istep),Pk_step(npbin/3*2,istep),Pk_step(npbin,istep)
            print*,'  Pk_nu',Pk_nu(1),Pk_nu(npbin/3),Pk_nu(npbin/3*2),Pk_nu(npbin)
            print*,'  Pk_nu1',Pk_nus(1,1),Pk_nus(npbin/3,1),Pk_nus(npbin/3*2,1),Pk_nus(npbin,1)
            print*,'  Pk_nu2',Pk_nus(1,2),Pk_nus(npbin/3,2),Pk_nus(npbin/3*2,2),Pk_nus(npbin,2)
            print*,'  Pk_nu3',Pk_nus(1,3),Pk_nus(npbin/3,3),Pk_nus(npbin/3*2,3),Pk_nus(npbin,3)
            print*,'  tf_F',tf_F(1),tf_F(npbin/3),tf_F(npbin/3*2),tf_F(npbin)
            print*,'  path',nupath//'*/*_'//trim(adjustl(str_z))//'.txt'
         endif
         open(211,file=nupath//'a_step.txt',status='replace',access='stream'); write(211) a_step(:istep); close(211)
         open(311,file=nupath//'Pk_nu/Pk_nu_'//trim(adjustl(str_z))//'.txt',status='replace',access='stream'); write(311) Pk_nu; close(311)
         open(311,file=nupath//'Pk_nus/Pk_nus_'//trim(adjustl(str_z))//'.txt',status='replace',access='stream'); write(311) Pk_nus; close(311)
         open(411,file=nupath//'Pk_m/Pk_m_'//trim(adjustl(str_z))//'.txt',status='replace',access='stream'); write(411) Pk_step(:,istep); close(411)
         open(511,file=nupath//'tf/Tf_nu_'//trim(adjustl(str_z))//'.txt',status='replace',access='stream'); write(511) tf_F; close(511)
         if (minval(tf_F) < 1-1e-6-f_nu .or. any((ieee_is_nan(tf_F)))) then
            if(head) print*,'tf_F err',minval(tf_F),1-f_nu
            if(head) print*,'Pk_nu'
            if(head) print*,Pk_nu
            if(head) print*,'Pk_cdm'
            if(head) print*,Pk_step(:,istep)
            if(head) print*,'Tf_nu'
            if(head) print*,tf_F
            error stop 'tf_F err'
         endif
      endif
      call toc(97)

      call tic(98)
      if (sim%a < 1) then
         tf_F_log = log(tf_F(:)[1])
         call nu_correction_coarse()
         pm%nwork = ngt   ;call nu_correction(tf2   ,Gk2   ,tile   )!,k_tf2   )
         pm%nwork = nft(2);call nu_correction(tf3_2 ,Gk3_2 ,subtile)!,k_tf3_2 )
         pm%nwork = nft(3);call nu_correction(tf3_4 ,Gk3_4 ,subtile)!,k_tf3_4 )
         pm%nwork = nft(4);call nu_correction(tf3_6 ,Gk3_6 ,subtile)!,k_tf3_6 )
         pm%nwork = nft(5);call nu_correction(tf3_8 ,Gk3_8 ,subtile)!,k_tf3_8 )
         pm%nwork = nft(6);call nu_correction(tf3_12,Gk3_12,subtile)!,k_tf3_12)
         ! if (head) print*,'  tf3_8 ',minval(tf3_8),maxval(tf3_8),interp_tf_F(kh_lin_log,real(pm%nwork*pi/subtile))
         ! if (head) print*,'  tf3_12',minval(tf3_12),maxval(tf3_12),interp_tf_F(kh_lin_log,real(pm%nwork*pi/subtile))
         ! stop
      endif
      sync all; call toc(98)
   endsubroutine

   subroutine nu_correction(tfk,Gk,size)!,k_tf)
      use omp_lib
      use variables
      use parameters
      implicit none

      integer i,j,k
      real kr,kx(3),size
      real tfk(pm%nwork/2+1,pm%nwork,pm%nwork)
      ! real k_tf(pm%nwork/2+1,pm%nwork,pm%nwork)
      real Gk(pm%nwork/2+1,pm%nwork,pm%nwork)

      tfk = 1
      ! k_tf =0
      if (Mass_nu > 0) then
         !$omp paralleldo default(shared) schedule(dynamic)&
         !$omp& private(i,j,k,kx,kr)
         do k=1,pm%nwork
            do j=1,pm%nwork
               do i=1,pm%nwork/2+1
                  kx=(mod([i,j,k]+pm%nwork/2,pm%nwork)-pm%nwork/2-1)*(2*pi)/size
                  if (i==1 .and. j==1 .and. k==1) cycle ! zero frequency
                  kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
                  tfk(i,j,k) = interp_tf_F(kh_lin_log,kr)*Gk(i,j,k)
               enddo
            enddo
         enddo
         !$omp endparalleldo
      endif
      ! call toc_nu
   endsubroutine

   subroutine nu_correction_coarse()
      use omp_lib
      use variables
      use parameters
      use pencil_fft
      implicit none

      integer i,j,k,ig,jg,kg
      real kr,kx(3)

      tf1 = 1
      ! call tic_nu
      if (Mass_nu > 0) then
         !$omp paralleldo default(shared) schedule(dynamic)&
         !$omp& private(kg,jg,ig,i,j,k,kx,kr)
         do k=1,npen
            do j=1,nw
               do i=1,nyquist+1
                  kg=(nn*(icz-1)+icy-1)*npen+k
                  jg=(icx-1)*nw+j
                  ig=i
                  kx=(mod([ig,jg,kg]+nyquist-1,nw_global)-nyquist)*(2*pi)/box
                  if (ig==1.and.jg==1.and.kg==1) cycle ! zero frequency
                  kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
                  tf1(i,j,k) = interp_tf_F(kh_lin_log,kr)*Gk1(i,j,k)
                  ! k_tf1(i,j,k) = kr
               enddo
            enddo
         enddo
         !$omp endparalleldo
      endif

      ! call toc_nu
   endsubroutine
endmodule
