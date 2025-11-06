program test
    use omp_lib
    use variables
    use power_nu
    implicit none

    real z_step(0:istep_max),a_next,z_next,Pk_now(npbin),Lxn(npbin),X,tau_step_read(0:istep_max),eta_step(0:istep_max),fnu_step(0:istep_max),Hz_step(0:istep_max),fnu,Hz,eta,eta_core(ncore),f_nr_step(0:istep_max),t_step(0:istep_max),dt00
    real ai,da_1,da_2
    integer i,icore,n_z,j,n
    integer,parameter :: calc = 1
    integer,parameter :: rr = 0

    real(16) ::  k1, k2, k3, k4, h,tau0
    integer(4) pm3t1,pm3t2,pm3tt1,pm3tt2,omplist(6,2) ,xx
    integer nnteam,nnnest

    ! do xx=1,6
    !     nnteam = omplist(xx,1)
    !     nnnest = omplist(xx,2)
    !     print*,xx,nnteam,nnnest
    ! enddo
    ! stop
    ! i = 100
    ! ai = 1./(i+1)
    ! print*,'z:',i,'f_nr',f_nr_a(ai)
    ! stop
    ! call system_clock(tt1,t_rate)
    ! do i=1,200
    !     ai = 1./(i+1)
    !     ! print*,'z:',i,'f_nr',f_nr_a(ai) 
    ! enddo
    ! call system_clock(tt2,t_rate);print*,'f_nu',real(tt2-tt1)/t_rate
    ! stop

    ! print*, 'export PYTHONUNBUFFERED=1'
    ! call system('export PYTHONUNBUFFERED=1')
    ! print*, 'python ./neutrinos/Pk.py'
    ! call system('python ./neutrinos/Pk.py')

    open(10,file=nupath//'s_a_tau_H.txt',form='formatted')
    read(10,*) stime
    read(10,*) s2a
    read(10,*) s2tau
    read(10,*) s2H
    close(10)
    s_fi = 0
    s_fi = sf_a(1./(1+200))

    print*,sf_a(1)
    print*,s2tau
    stop

    do while(sim%a<1 .and. sim%timestep<istep_max)
        ! call expansion(a_step(sim%timestep-1),dt00,da_1,da_2)
        da_1 = expansion(a_step(sim%timestep-1),dt00)
        sim%a = da_1+sim%a
        a_step(sim%timestep) = sim%a
        t_step(sim%timestep) = t_step(sim%timestep-1)+dt00
        sim%tau = sim%tau + dtau_a(da_1)
        print*,t_step(sim%timestep),sf_a(sim%a),sim%a,sim%tau,H_a(sim%a)
        ! print*,''
        ! print*,''
        sim%timestep = sim%timestep+1
        ! sim%a=a_step(sim%timestep)
    enddo
    open(18,file='a_step.txt',status='replace',access='stream'); write(18) a_step(:sim%timestep); close(18)
    open(18,file='t_step.txt',status='replace',access='stream'); write(18) t_step(:sim%timestep); close(18)
    ! print*,sim%timestep,a_step(sim%timestep-1),a_step(sim%timestep-2)
    ! print*,''
    ! print*,''
    ! print*,''
    ! print*,''
    ! print*,''
    ! print*,''
    ! print*,''
    ! print*,''
    ! print*,''
    ! print*,''
    ! print*,''
    ! print*,'t'
    ! print*,t_step(:101)
    ! print*,'a'
    ! print*,a_step(:101)

    stop
    ! call system_clock(tt1,t_rate)
    ! print*,f_nu_a(1./203)
    ! call system_clock(tt2,t_rate);print*,'f_nu',real(tt2-tt1)/t_rate
    ! call system('ls ../output1/1200_128_nu_math1/neutrinos/')
    ! stop

    ! do i=10,30
    !     ai = 10.**(-i)
    !     print*,i,H_a(ai),f_nu_a(ai),omega_nu*f_nu_a(ai)*h0**2,'N',1.0_16/(ai*(1.0_16))**4,HUGE(1.0_16)
    ! enddo
    ! stop

    ! call system_clock(tt1,t_rate)
    ! ai = (1./203)
    ! n = ncore*3e2
    ! h = ai/n
    ! 99 FORMAT (E15.7,E15.7)
    ! print 99,1e-10,ai
    ! 100 FORMAT (A,I4,E15.7,E15.7)
    ! 101 FORMAT (A,I4,I4,E15.7,E15.7)
    ! eta_core = 0
    ! !$omp paralleldo num_threads(ncore) default(shared) schedule(dynamic)&
    ! !$omp& private(i,j,k1,k2,k3,k4,icore,ai)
    ! do j=0,ncore-1
    !     icore=omp_get_thread_num()+1
    !     ai = 1.e-10+h*j*n/ncore
    !     print 99,ai,ai+h
    !     k4 = h * 1/H_a(ai)/(ai)**2
    !     ! print 101,'j ',j*n/ncore,(j+1)*n/ncore-1,1.e-30+h*j*n/ncore,1.e-30+h*((j+1)*n/ncore)
    !     do i=j*n/ncore,(j+1)*n/ncore-1
    !         k1 = k4
    !         ai = ai+h/2
    !         k2 = h * 1/H_a(ai)/(ai)**2
    !         ai = ai+h/2
    !         k4 = h * 1/H_a(ai)/(ai)**2
    !         eta_core(icore) = eta_core(icore) + (k1+4*k2+k4)
    !     enddo
    ! enddo
    ! !$omp endparalleldo
    ! ! h = ai/n
    ! eta = sum(eta_core(:))/6
    ! call system_clock(tt2,t_rate);print*,'f_nu',real(tt2-tt1)/t_rate,'s    eta_p',eta,eta/0.0029704539532831135-1
    ! stop

    ! !!!!! test eta
    ! if (0) then
    !     eta_step = -1
    !     open(13,file='./neutrinos/eta200.txt',status='old')
    !     do i=0,istep_max-1
    !     read(13,end=83,fmt='(f20.10)') eta_step(i)
    !     enddo
    !     83 close(13)
    !     print*,'eta',eta_step(0),eta_step(1),eta_step(2)

    !     fnu_step = -1
    !     open(13,file='./neutrinos/fnu200.txt',status='old')
    !     do i=0,istep_max-1
    !     read(13,end=84,fmt='(f20.10)') fnu_step(i)
    !     enddo
    !     84 close(13)
    !     print*,'fnu',fnu_step(0),fnu_step(1),fnu_step(2)

    !     Hz_step = -1
    !     open(13,file='./neutrinos/Hz200.txt',status='old')
    !     do i=0,istep_max-1
    !     read(13,end=85,fmt='(f20.10)') Hz_step(i)
    !     enddo
    !     85 close(13)
    !     print*,'Hz',Hz_step(0),Hz_step(1),Hz_step(2)

    !     z_step = -1
    !     open(11,file='./neutrinos/z200.txt',status='old')
    !     do i=0,istep_max-1
    !     read(11,end=81,fmt='(f6.4)') z_step(i)
    !     n_z = n_z+1
    !     enddo
    !     81 close(11)
    !     print*,'z',z_step(0),z_step(1),z_step(2)

    !     print*,''
    !     print*,''
    !     print*,'',omega_nu*h0**2,omega_nu

    !     call system_clock(tt1,t_rate)
    !     ai = 1./(1+z_step(0))
    !     n = ncore*3e2
    !     h = ai/n
    !     99 FORMAT (E15.7,E15.7)
    !     print 99,1e-10,ai
    !     100 FORMAT (A,I4,E15.7,E15.7)
    !     101 FORMAT (A,I4,I4,E15.7,E15.7)
    !     eta_core = 0
    !     !$omp paralleldo num_threads(ncore) default(shared) schedule(dynamic)&
    !     !$omp& private(i,j,k1,k2,k3,k4,icore,ai)
    !     do j=0,ncore-1
    !         icore=omp_get_thread_num()+1
    !         ai = 1.e-10+h*j*n/ncore
    !         ! print 99,ai,ai+h
    !         k4 = 1/H_a(ai)/(ai)**2
    !         ! print 101,'j ',j*n/ncore,(j+1)*n/ncore-1,1.e-30+h*j*n/ncore,1.e-30+h*((j+1)*n/ncore)
    !         do i=j*n/ncore,(j+1)*n/ncore-1
    !             k1 = k4
    !             ai = ai+h/2
    !             k2 = 1/H_a(ai)/(ai)**2
    !             ai = ai+h/2
    !             k4 = 1/H_a(ai)/(ai)**2
    !             eta_core(icore) = eta_core(icore) + (k1+4*k2+k4)
    !         enddo
    !     enddo
    !     !$omp endparalleldo
    !     eta = h * sum(eta_core(:))/6
    !     call system_clock(tt2,t_rate);print*,'eta',real(tt2-tt1)/t_rate
        
    !     print*,ai,eta
    !     ! stop
    !     do i= 0,n_z-2
    !         sim%a = 1./(1+z_step(i))
    !         da = 1./(1+z_step(i+1))-1./(1+z_step(i))

    !         call system_clock(tt1,t_rate)
    !         fnu = f_nu_a(sim%a )
    !         Hz = H_a(sim%a )
    !         dtau = 0
    !         ai = sim%a
    !         n = 4
    !         h = da/n
    !         k4 = 1/H_a(ai)/(ai)**2
    !         do j= 0,n-1
    !             k1 = k4
    !             ai = ai+h
    !             k4 = 1/H_a(ai)/(ai)**2
    !             dtau = dtau+(k1+k4)!+4*k2)
    !         enddo
    !         dtau = h*dtau/2
    !         call system_clock(tt2,t_rate);print*,'z',z_step(i),ai,i,'time',real(tt2-tt1)/t_rate,'s'

    !         print*,'    fnu',fnu,fnu_step(i),fnu/fnu_step(i)-1
    !         print*,'    Hz',Hz,Hz_step(i),Hz/Hz_step(i)-1
    !         print*,'    deta',299792.458*dtau,eta_step(i+1)-eta_step(i),(299792.458*dtau)/(eta_step(i+1)-eta_step(i))-1
    !         print*,'    eta',299792.458*eta,eta_step(i),(299792.458*eta)/eta_step(i)-1
    !         print*,''
    !         print*,''
    !         print*,''
            
    !         eta = eta+dtau
    !     enddo 

    !     stop
    ! endif

    !!!!! read
    if (1) then 
        head = 1
        ! call system_clock(tt1,t_rate)
        ! print*,'f_un(z=200):',f_nu_a(1./201),npbin
        ! call system_clock(tt2,t_rate);print*,'f_nu',real(tt2-tt1)/t_rate
        ! print*,'H_a(z=200):',H_a(1./201)
        n_z = 0
        z_step = -1
        open(11,file=nupath//'z_values.txt',status='old')
        do i=0,istep_max-1
        read(11,end=91,fmt='(f8.4)') z_step(i)
        n_z = n_z+1
        enddo
        91 close(11)
        print*,'z',z_step(0),z_step(1),z_step(2)

        a_step = -1
        open(12,file=nupath//'a_values.txt',status='old')
        do i=0,istep_max-1
        read(12,end=92,fmt='(f15.12)') a_step(i)
        enddo
        92 close(12)
        print*,'a',a_step(0),a_step(1),a_step(2)

        tau_step_read = -1
        open(13,file=nupath//'tau_values.txt',status='old')
        do i=0,istep_max-1
        read(13,end=93,fmt='(f15.12)') tau_step_read(i)
        enddo
        93 close(13)
        print*,'tau',tau_step_read(0),tau_step_read(1),tau_step_read(2)


        kh_lin = -1
        open(13,file=nupath//'k_values.txt',status='old')
        do i=1,npbin
        read(13,end=95,fmt='(f15.12)') kh_lin(i)
        enddo
        95 close(13)
        kh_lin_log = log(kh_lin)
        print*,'k',kh_lin(1),kh_lin(2),kh_lin(npbin)



        sim%timestep = 0

        z_checkpoint(1)=200
    endif

    






    ! !!!!!!!!! initialize !!!!!!!!!!


    call geometry
    allocate(Gk1(nc_global/2+1,nc,npen))
    allocate(Gk2(ngt/2+1,ngt,ngt))
    allocate(Gk3_2(nft(2)/2+1,nft(2),nft(2)))
    allocate(Gk3_4(nft(3)/2+1,nft(3),nft(3)))
    allocate(Gk3_6(nft(4)/2+1,nft(4),nft(4)))
    allocate(Gk3_8(nft(5)/2+1,nft(5),nft(5)))

    Gk1 = 1
    Gk2 = 1
    Gk3_2 = 1
    Gk3_4 = 1
    Gk3_6 = 1
    Gk3_8 = 1

    open(10,file=nupath//'s_a_tau_H.txt',form='formatted')
    read(10,*) stime
    read(10,*) s2a
    read(10,*) s2tau
    read(10,*) s2H
    close(10)
    s2tau = s2tau/299792.458
    ! print*,''
    ! print*,''
    ! print*,'s'
    ! print*,stime(:101)
    ! print*,''
    ! print*,''
    ! print*,'a'
    ! print*,s2a(:101)
    ! print*,''
    ! print*,''
    ! print*,'tau'
    ! print*,s2tau(:101)
    ! print*,''
    ! print*,''
    ! print*,'H'
    ! print*,s2H(:101)
    ! print*,''
    ! print*,''
    s_fi = 0
    print*,1./(1+z_checkpoint(1))
    print*,sf_a(1.),sf_a(1./(1+z_checkpoint(1)))
    s_fi = sf_a(1./(1+z_checkpoint(1)))
    print*,s_fi

    !nu
    s_f=0
    sim%cur_powerpoint=1
    z_powerpoint=-9999
    power_step=.false.
    ! a_step = 0
    Pk_step = 0
    tau_step = 0
    Pk_cb_check =0
    Pk_nu_check =0
    if (Mass_nu > 0) then
        if (head) then
        print*,'init nu_info'
        print*,'  npb ncb nnb',npbin,ncbin,nnbin
        print*,'  tile subtile',tile,subtile
        print*,'  z_nu_start',1/a_nu-1
        print*,'  nupath:',nupath
        endif
        
        H_i = H_a(1/(z_checkpoint(1)+1))
        H_0 = h0*100

        open(16,file=nupath//'z_powerpoint.txt',status='old')
        do i=1,nmax_redshift-1
        read(16,end=72,fmt='(f15.12)') z_powerpoint(i)
        enddo
        72 n_powerpoint=i-1
        close(16)
        if (n_powerpoint==0) stop 'z_powerpoint.txt empty'  
        n_powerpoint=n_powerpoint[1]
        z_powerpoint(:)=z_powerpoint(:)[1]

        do i=1,3 ! Read the three earliest Pk_nu for linear interpolation
            write(str_z,'(f8.4)') z_powerpoint(i)
            open(10,file=nupath//'Pk_cb_'//trim(adjustl(str_z))//'.txt',form='formatted')
            read(10,*) Pk_cb_check(:,i)
            close(10)
            open(10,file=nupath//'Pk_nu_'//trim(adjustl(str_z))//'.txt',form='formatted')
            read(10,*) Pk_nu_check(:,i)
            close(10)
        enddo
        Pk_step(:,0)=Pk_cb_check(:,1)

        open(10,file=nupath//'Pk_nu_ic.txt',form='formatted')
        read(10,*) Pk_nu_ic
        close(10)

        kh_lin = -1
        open(13,file=nupath//'k_values.txt',status='old')
        do j=1,npbin
        read(13,end=75,fmt='(f15.12)') kh_lin(j)
        enddo
        75 close(13)
        kh_lin_log = log(kh_lin)
        kh_lin_log = log(kh_lin)
        if (head) print*,'k',kh_lin(1),kh_lin(2),kh_lin(npbin)
    endif

    sim%omega_m = omega_m
    sim%h0 = h0
    sim%omega_l = omega_l
    H_0 = 70

    print*,image,head,omhsq0


    sim%a=1./200
    print*,sim%a,H_a(sim%a),dtau_a(0.1),sf_a(sim%a)
    sim%a=1./20
    print*,sim%a,H_a(sim%a),dtau_a(0.1),sf_a(sim%a)
    sim%a=1./10
    print*,sim%a,H_a(sim%a),dtau_a(0.1),sf_a(sim%a)
    sim%a=1./5
    print*,sim%a,H_a(sim%a),dtau_a(0.1),sf_a(sim%a)
    sim%a=1./2
    print*,sim%a,H_a(sim%a),dtau_a(0.1),sf_a(sim%a)
    sim%a=1.
    print*,sim%a,H_a(sim%a),dtau_a(0.1),sf_a(sim%a)
    stop
    if (calc) then
        !init tau(0)
        tau_step = 0
        call system_clock(tt1,t_rate)
        ai = 1./(1+z_checkpoint(1))
        tau0 = s2tau(1)
        call system_clock(tt2,t_rate);print*,'eta',real(tt2-tt1)/t_rate
    
        sim%tau = tau0
        istep = 0
        tau_step(istep) = sim%tau
        print*,'step',istep,tau_step(istep),tau_step_read(istep),tau_step(istep)/tau_step_read(istep)-1

        do istep=sim%timestep,istep_max
            print*,''
            print*,''
            print*,'step',istep,z_step(istep)
            ! write(str_z,'(f8.4)') z_step(istep)
            ! open(10,file=nupath//'Pk_nu_'//trim(adjustl(str_z))//'.dat.txt',form='formatted')
            ! read(10,*) Pk_step(:,istep)
            ! close(10)
            power_step=.false.
            sim%a = a_step(istep)
            ! print*,nupath//'Pk_nu_'//trim(adjustl(str_z))//'.dat.txt'
            ! print*,istep,str_z,'Pk_istep',Pk_step(npbin-4:npbin,istep)


            ! tau_step(istep) = tau_step_read(istep)

            if (Mass_nu > 0) then

                da = a_step(istep+1)-a_step(istep)
                dtau = dtau_a(sim%a+da)
                sim%tau = sim%tau + dtau
                tau_step(istep) = sim%tau
                tau_step(istep) = tau_step_read(istep)

                z_next = z_powerpoint(sim%cur_powerpoint)
                a_next = 1.0/(1+z_next)
                if (a_step(istep)>=a_next) then
                    power_step = .true.
                    print*,'    ************** check point',sim%cur_powerpoint,z_powerpoint(sim%cur_powerpoint)
                endif
            endif
            print*,'    tau',sim%tau,tau_step_read(istep),(tau_step(istep)-tau_step_read(istep))/tau_step_read(istep)
            ! call interp_Pk_CDM

            ! if (istep > 0) then
            !     s_f = s_f+ 0.5*((1/a_step(istep)+1/a_step(istep-1))*(tau_step(istep)-tau_step(istep-1)) ) !superconformal time s=integral(d(z)/a(z)) 
            ! endif
            call get_tf_cb2matter
            f_nr_step(istep) = f_nr

            ! print*,'sf=',s_f,istep,'Pk_step',Pk_step(1:4,istep)-Pk_now(1:4)
            if (z_step(istep) <= 0.0)  exit
        enddo
        ! print*,'        sf=',s_f,istep


        open(10,file=nupath//'f_nr.txt',status='replace',access='stream'); write(10) f_nr_step; close(10)
        ! open(10,file='./tf_test/Tf1.txt',status='replace',access='stream'); write(10) tf1; close(10)
        ! open(10,file='./tf_test/Tf2.txt',status='replace',access='stream'); write(10) tf2; close(10)
        ! open(10,file='./tf_test/Tf3_2.txt',status='replace',access='stream'); write(10) tf3_2; close(10)
        ! open(10,file='./tf_test/Tf3_4.txt',status='replace',access='stream'); write(10) tf3_4; close(10)
        ! open(10,file='./tf_test/Tf3_6.txt',status='replace',access='stream'); write(10) tf3_6; close(10)
        ! open(10,file='./tf_test/Tf3_8.txt',status='replace',access='stream'); write(10) tf3_8; close(10)

        ! open(10,file='./tf_test/k_Tf1.txt',status='replace',access='stream'); write(10) k_tf1; close(10)
        ! open(10,file='./tf_test/k_Tf2.txt',status='replace',access='stream'); write(10) k_tf2; close(10)
        ! open(10,file='./tf_test/k_Tf3_2.txt',status='replace',access='stream'); write(10) k_tf3_2; close(10)
        ! open(10,file='./tf_test/k_Tf3_4.txt',status='replace',access='stream'); write(10) k_tf3_4; close(10)
        ! open(10,file='./tf_test/k_Tf3_6.txt',status='replace',access='stream'); write(10) k_tf3_6; close(10)
        ! open(10,file='./tf_test/k_Tf3_8.txt',status='replace',access='stream'); write(10) k_tf3_8; close(10)
        ! open(10,file='Pk_nu_interp.txt',status='replace',access='stream'); write(10) Pk_nu; close(10)
        ! open(10,file='Tf_nu_interp.txt',status='replace',access='stream'); write(10) tf_F; close(10)

        ! Pk_now = (f_nu*Pk_nu**0.5+(1-f_nu)*Pk_step(:,istep)**0.5)**2
        ! open(10,file='Pk_mater_interp.txt',status='replace',access='stream'); write(10) Pk_now; close(10)
        ! open(10,file='Pk_CB.txt',status='replace',access='stream'); write(10) Pk_step(:,istep); close(10)
    endif

    if (rr) then
        !nu
        s_f=0
        sim%timestep = 0
        sim%cur_powerpoint=1
        power_step=.false.
        s_f=0
        tau_step = tau_step_read
        do istep=sim%timestep,istep_max
            ! write(str_z,'(f8.4)') z_step(istep)
            ! open(10,file=nupath//'Pk_nu_'//trim(adjustl(str_z))//'.dat.txt',form='formatted')
            ! read(10,*) Pk_now
            ! close(10)
            power_step=.false.
            sim%a = a_step(istep)
            ! print*,nupath//'Pk_nu_'//trim(adjustl(str_z))//'.dat.txt'
            ! print*,istep,str_z,'Pk_istep',Pk_step(npbin-4:npbin,istep)



            if (Mass_nu > 0) then
                z_next = z_powerpoint(sim%cur_powerpoint)
                a_next = 1.0/(1+z_next)
                if (a_step(istep)>=a_next) then
                    power_step = .true.
                    print*,'r************** check point',sim%cur_powerpoint,z_powerpoint(sim%cur_powerpoint)
                endif
            endif
            call interp_Pk_CDM
            ! call get_tf_cb2matter

            ! if (istep > 0) then
            !     s_f = s_f+ 0.5*((1/a_step(istep)+1/a_step(istep-1))*(tau_step(istep)-tau_step(istep-1)) ) !superconformal time s=integral(d(z)/a(z)) 
            ! endif
            ! print*,'************************'
            ! print*,istep,z_step(istepx),'Pk_step',(Pk_step(1:4,istep)-Pk_now(1:4))/Pk_now(1:4),(Pk_step(npbin-4:npbin,istep)-Pk_now(npbin-4:npbin))/Pk_now(npbin-4:npbin)

            ! print*,'sf=',s_f,istep,'Pk_step',Pk_step(1:4,istep)-Pk_now(1:4)
            if (z_step(istep) <=20 ) then
                write(str_z,'(f8.4)') z_step(istep)
                call get_tf_cb2matter
                print*,'r************** save point',sim%cur_powerpoint,z_powerpoint(sim%cur_powerpoint),str_z
                open(10,file='./tf_test/Pk_nu_interp_r'//trim(adjustl(str_z))//'.txt',status='replace',access='stream'); write(10) Pk_nu; close(10)
                open(10,file='./tf_test/Tf_nu_interp_r'//trim(adjustl(str_z))//'.txt',status='replace',access='stream'); write(10) tf_F; close(10)
                open(10,file='./tf_test/Tf1'//trim(adjustl(str_z))//'.txt',status='replace',access='stream'); write(10) tf1; close(10)
                open(10,file='./tf_test/Tf2'//trim(adjustl(str_z))//'.txt',status='replace',access='stream'); write(10) tf2; close(10)
                open(10,file='./tf_test/Tf3_2'//trim(adjustl(str_z))//'.txt',status='replace',access='stream'); write(10) tf3_2; close(10)
                open(10,file='./tf_test/Tf3_4'//trim(adjustl(str_z))//'.txt',status='replace',access='stream'); write(10) tf3_4; close(10)
                open(10,file='./tf_test/Tf3_6'//trim(adjustl(str_z))//'.txt',status='replace',access='stream'); write(10) tf3_6; close(10)
                open(10,file='./tf_test/Tf3_8'//trim(adjustl(str_z))//'.txt',status='replace',access='stream'); write(10) tf3_8; close(10)
            endif
            if (z_step(istep) <= 20)  exit
        enddo
        ! call get_tf_cb2matter
        print*,'************************',z_step(istep)
        print*,sim%omega_m,f_nu,'r_sf=',s_f,f_nu,istep,Pk_nu(1:5),npbin,nft(2),nft(3),nft(4),nft(5),'gk1',nc_global/2+1,nc,npen



        ! open(10,file='./tf_test/k_Tf1.txt',status='replace',access='stream'); write(10) k_tf1; close(10)
        ! open(10,file='./tf_test/k_Tf2.txt',status='replace',access='stream'); write(10) k_tf2; close(10)
        ! open(10,file='./tf_test/k_Tf3_2.txt',status='replace',access='stream'); write(10) k_tf3_2; close(10)
        ! open(10,file='./tf_test/k_Tf3_4.txt',status='replace',access='stream'); write(10) k_tf3_4; close(10)
        ! open(10,file='./tf_test/k_Tf3_6.txt',status='replace',access='stream'); write(10) k_tf3_6; close(10)
        ! open(10,file='./tf_test/k_Tf3_8.txt',status='replace',access='stream'); write(10) k_tf3_8; close(10)

        ! open(10,file='Pk_nu_interp_r.txt',status='replace',access='stream'); write(10) Pk_nu; close(10)
        ! Pk_now = (f_nu*Pk_nu**0.5+(1-f_nu)*Pk_step(:,istep)**0.5)**2
        ! open(10,file='Pk_mater_interp_r.txt',status='replace',access='stream'); write(10) Pk_now; close(10)
        ! open(10,file='Pk_CB_r.txt',status='replace',access='stream'); write(10) Pk_step(:,istep); close(10)
        ! open(10,file='Pk_CB_r_200.txt',status='replace',access='stream'); write(10) Pk_step(:,0); close(10)
    endif


    ! open(1,file='./Tf_nu_interp_r.txt',status='old',access='stream')
    ! read(1) tf_F
    ! close(1)
    ! tf_F_log = log(tf_F)

    ! i = 400
    ! h = 1.1
    ! ai = interp_tf_F(kh_lin_log,kh_lin(i)*h)
    ! print*,i,'tf',ai,kh_lin(i)
    
    ! i = 11
    ! h = 1.1
    ! ai = interp_tf_F(kh_lin_log,kh_lin(i)*h)
    ! print*,i,'tf',ai

end program test
