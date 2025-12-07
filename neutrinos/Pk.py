import camb,sys,os,re,time
import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d

test = 0

all_reps = False
all_camb = True
n=50
ni=2

def match_para(para):
    # 读取Fortran文件
    file_path = './parameters.f90'
    with open(file_path, 'r') as file:
        content = file.read()

    # 使用正则表达式匹配模式
    pattern = r'parameter\s*::\s*%s\s*=\s*([^\s]+)'%para
    match = re.search(pattern, content)

    if match:
        variable_value = match.group(1).strip()
        return float(variable_value)
    else:
        print(f"Pattern({para}) not found in the file.")
        sys.exit()

def get_z(n): #get array a
    z = [0]*n
    a = [0]*n
    z[0] = z_max
    a = np.linspace(1/(1+z[0]),1,n+1)
    z=1/a-1
    return z[::-1]
    
def run_reps():
    os.system('rm %s/IC/IC*'%(nupath))
    file_path = './neutrinos/IC.ini'
    with open(file_path, 'r') as f:
        IC_old = f.read()

    # 使用正则表达式匹配模式
    pattern = r'z_initial \s*=(.*)'
    replacement = "z_initial          =  10000"
    IC_new = re.sub(pattern, replacement, IC_old)

    pattern = r'outputfile         \s*=(.*)'
    replacement = "outputfile         =  %s/IC/IC"%(nupath)
    IC_new = re.sub(pattern, replacement, IC_new)

    if Mass_nu>0:
        IC_Neff = 'Neff               =  0.00641'
        IC_N_nu = 'N_nu               =  3.0'
        IC_M_nu = 'M_nu               =  %.2f'%Mass_nu
    else:
        IC_Neff = 'Neff               =  3.046'
        IC_N_nu = 'N_nu               =  0.0'
        IC_M_nu = 'M_nu               =  0.0'
    pattern = r'Neff \s*=(.*)'
    replacement = IC_Neff
    IC_new = re.sub(pattern, replacement, IC_new)
    pattern = r'N_nu \s*=(.*)'
    replacement = IC_N_nu
    IC_new = re.sub(pattern, replacement, IC_new)
    pattern = r'M_nu \s*=(.*)'
    replacement = IC_M_nu
    IC_new = re.sub(pattern, replacement, IC_new)

    pattern = r'As \s*=(.*)'
    replacement = "As                 =  %.4e"%As
    IC_new = re.sub(pattern, replacement, IC_new)
    pattern = r'ns \s*=(.*)'
    replacement = "ns                 =  %.4f"%ns
    IC_new = re.sub(pattern, replacement, IC_new)


    pattern = r'OC0 \s*=(.*)'
    replacement = "OC0                =  %.10f"%(omega_cdm-Mass_nu/93.14/((H0/100)**2))
    IC_new = re.sub(pattern, replacement, IC_new)
    pattern = r'OB0 \s*=(.*)'
    replacement = "OB0                =  %.10f"%omega_bar
    IC_new = re.sub(pattern, replacement, IC_new)

    if all_reps:
        z_need = z_nonlin
    else:
        z_need = np.loadtxt('z_checkpoint.txt')

    n_z = np.array(z_need).shape[0]

    pattern = r'output_number \s*=(.*)'
    replacement = "output_number      =  %d"%n_z
    IC_new = re.sub(pattern, replacement, IC_new)

    z_str = ''
    for i in z_need:
        z_str+=' %.4f'%i
    pattern = r'z_output \s*=(.*)'
    replacement = "z_output           = "+z_str
    IC_new = re.sub(pattern, replacement, IC_new)
    with open(file_path, 'w') as f:
        f.write(IC_new)
    
    if not os.path.exists('./reps-master/reps'):
        os.system('cd ./reps-master/ && make && ./reps ../neutrinos/IC.ini && cd ..')
    else:
        os.system('cd ./reps-master/ && ./reps ../neutrinos/IC.ini && cd ..')
    print('REPS done!')

def interp_pk(z,s,k_new):
    Pk0 = np.loadtxt('%s/IC/IC_P%s_rescaled_z%.4f.txt'%(nupath,s,z))[:,1]
    interp_00 = interp1d(k,Pk0, kind='linear', bounds_error=False, fill_value="extrapolate")
    new_Pk0 =interp_00(k_new)
    new_Pk0[new_Pk0 <= 0] = 1e-10
    return new_Pk0

def get_f_nr(z):
    y = Mass_nu/3/(sigma_nu*N_eff*k_b*T_gama*(1+z))
    Fy=integrate.quad(lambda u: (u**2*np.sqrt(u**2+y**2))/(1+np.exp(u)), 0,300)[0]
    f1=integrate.quad(lambda u: (u**2)/((1+np.exp(u))*np.sqrt(u**2+y**2)), 0,300)[0]
    f_nr = y**2*f1/Fy#*f_nu
    return f_nr
    
def get_Pk_nonlin_CDM(n):
    Pk_nonlin_CDM = [0]*len(z_nonlin)
    Pk_nonlin_nu = [0]*len(z_nonlin)
    nz = int(np.ceil(len(z_nonlin)/n_PK))
    for i in range(0,nz):
        z_n=z_nonlin[i*n_PK:(i+1)*n_PK]
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, omk=omk, neutrino_hierarchy=neutrino_hierarchy, num_massive_neutrinos=num_massive_neutrinos, mnu=mnu, nnu=nnu, standard_neutrino_neff=standard_neutrino_neff)
        pars.omch2=omch2-pars.omnuh2
        pars.InitPower.set_params(As=As,ns=ns)
        print('i = %d/%d;z_max = %.3f;len = %d'%(i+1,nz,z_n[-1],len(z_n)))
        pars.set_matter_power(redshifts=z_n, kmax=k_ic_max, nonlinear=True)
        pars.NonLinear = camb.model.NonLinear_both
        results = camb.get_results(pars)
        kh_nonlin, z_nonlin[i*n_PK:i*n_PK+len(z_n)],Pk_nonlin_CDM[i*n_PK:i*n_PK+len(z_n)]= results.get_matter_power_spectrum(minkh=k_ic_min, maxkh=k_ic_max, npoints = npbin,var1='delta_nonu',var2='delta_nonu')
        _, _                                        ,Pk_nonlin_nu[i*n_PK:i*n_PK+len(z_n)]= results.get_matter_power_spectrum(minkh=k_ic_min, maxkh=k_ic_max, npoints = npbin,var1='delta_nu'  ,var2='delta_nu')
        pars = 0
    return Pk_nonlin_CDM[::-1],Pk_nonlin_nu[::-1],kh_nonlin,z_nonlin[::-1]

#get cube parm
if (1):
    H0=match_para('h0')*100
    omega_bar=match_para('omega_bar')
    omega_cdm=match_para('omega_cdm')
    omega_r=match_para('omega_r')
    omk=0.0

    file_path = './parameters.f90'
    with open(file_path, 'r') as file:
        content = file.read()

    # 使用正则表达式匹配模式
    pattern = r'parameter\s*::\s*%s\s*=*([^\s]+)'%'m_nu'
    match = re.search(pattern, content)

    if match:
        m_nus = match.group(1).strip()
        str_cleaned = m_nus.replace("[", "").replace("]", "").replace(",", " ")
        m_nus = np.array(str_cleaned.split(), dtype=float)
    else:
        print("Pattern (m_nu) not found in the file.")
        sys.exit()

    mnu=m_nus.sum()
    if (mnu>0) :
        if (np.all(m_nus == m_nus[0])):
            neutrino_hierarchy = 'degenerate'
            num_massive_neutrinos=3
        elif (m_nus[0] < m_nus[2]):
            neutrino_hierarchy = 'normal'
            num_massive_neutrinos=3
        elif (m_nus[0] > m_nus[2]):
            neutrino_hierarchy = 'inverted'
            num_massive_neutrinos=2
    else:
        neutrino_hierarchy = 'degenerate'
        num_massive_neutrinos=0
    nnu=3.044
    standard_neutrino_neff=match_para('N_eff')
    ns=match_para('n_s')
    As=match_para('A_s')
    z_max = np.loadtxt('./z_checkpoint.txt')[0]

    pi=np.pi
    ratio_cs=match_para('ratio_cs')
    ratio_sf=np.array([1,2,4,6,8,12,16])
    ng = match_para('ng')
    box = match_para('box')
    nns = match_para('nns')
    nnt = match_para('nnt')
    nn = match_para('nn')
    ngb = match_para('ngb')
    istep_max = match_para('istep_max')
    nc =  ng/nnt
    ng_global = ng*nn
    k_smooth = max(10.0,ng_global*pi/box)
    nfg = min(np.ceil((k_smooth*box/pi/nn) / nc)*nc,ng*1.)
    nfg_global=nfg*nn
    npbin=int(nfg_global/2*np.sqrt(3.))+21
    kmin = 2*pi/box
    kmax = pi/box*nfg_global*np.sqrt(3)
    ombh2 = omega_bar*(H0/100)**2
    omch2 = omega_cdm*(H0/100)**2

    sigma_nu = (4/11)**(1/3)
    N_nu = 3
    N_eff = standard_neutrino_neff
    k_b = 8.617342e-5
    T_gama = 2.7255
    Mass_nu = mnu
    omega_m = omega_bar+omega_cdm
    omega_cb = omega_bar+omega_cdm-Mass_nu/93.14/((H0/100)**2)
    omega_l = 1-omega_bar-omega_cdm-omega_r
    f_nu = Mass_nu/93.14/((H0/100)**2)/omega_m
    omHsq = 2/3*np.sqrt(1/omega_m)/H0


    file_path = './parameters.f90'
    with open(file_path, 'r') as file:
        content = file.read()

    # 使用正则表达式匹配模式
    pattern = r'parameter\s*::\s*opath\s*=\s*([^\s]+)'
    match = re.search(pattern, content)
    opath = os.path.expanduser(match.group(1).strip()[1:-1])
    print(f"Variable value of opath : {opath}")
    nupath =  opath+"neutrinos"
    print(f"Variable value of nupath : {nupath}")
    print('\n'+('+'*40+'\n')*2)
    print('Cosmology  Paras:\n\n   omega_r:   %.6f\n   omega_b:   %.6f\n   omega_c:   %.6f\n   omega_l:   %.6f\n   mass_nu:   %.3f              eV\n  mass_nus:   %.3f %.3f %.3f %s  eV\n      f_nu:   %.6f\n\n'%(omega_r,omega_bar,omega_cdm,omega_l,Mass_nu,m_nus[0],m_nus[1],m_nus[2],neutrino_hierarchy,f_nu))
    print('Simulation Paras:\n\n     nfg:   %d\n     npbin:   %d\n     kmin:   %.3f\n      kmax:   %.3f\n\n\n'%(nfg,npbin,kmin,kmax))



    kh_nonlin = np.concatenate((np.array([1.0, 1.4142135623730951, 1.7320508075688772, 2.0, 
                                        2.23606797749979, 2.449489742783178, 2.8284271247461903, 
                                        3.0, 3.1622776601683795, 3.3166247903554, 3.4641016151377544, 
                                        3.605551275463989, 3.7416573867739413, 4.0, 4.123105625617661, 
                                        4.242640687119285, 4.358898943540674, 4.47213595499958, 
                                        4.58257569495584, 4.69041575982343, 4.898979485566356, 
                                        5.0, 5.0990195135927845, 5.196152422706632, 5.385164807134504, 
                                        5.477225575051661]) * 2 * np.pi / box
                                , np.exp(np.linspace(np.log(6*kmin), np.log(kmax),  npbin - 26))))


# print(kh_nonlin.shape)
# exit()



#mkdir
try:
    os.system('mkdir -p '+nupath+'/TF')
    os.system('mkdir -p '+nupath+'/Pk_m')
    os.system('mkdir -p '+nupath+'/Pk_nu')
    os.system('mkdir -p '+nupath+'/Pk_nus')
    os.system('mkdir -p %s/IC'%(nupath))
except:
    None
    
#set n_z
if test:
    n=n*ni
else:
    n=n
n_PK = 150 # max num of Pk(z) at once
k_ic_min = 1e-4
k_ic_max=max(1e2,kmax*1.5)

z_nonlin=get_z(n)



#get Expansion History
if (all_camb):
    print('all_camb')
    par = camb.CAMBparams()
    par.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, omk=omk, neutrino_hierarchy=neutrino_hierarchy, num_massive_neutrinos= num_massive_neutrinos, mnu= mnu, nnu= nnu, standard_neutrino_neff= standard_neutrino_neff)
    par.omch2=omch2-par.omnuh2
    par.InitPower.set_params(As=As,ns=ns)
    par.set_matter_power(redshifts=[z_max], kmax=k_ic_max, nonlinear=True)
    par.NonLinear = camb.model.NonLinear_both
    result = camb.get_results(par)
    kh_ic, _,Pk_cb_ic= result.get_matter_power_spectrum(minkh=k_ic_min, maxkh=k_ic_max, npoints = npbin,var1='delta_nonu',var2='delta_nonu')
    Pk_cb_ic = Pk_cb_ic[0]
    if (neutrino_hierarchy == 'degenerate'):
        print('neutrino_hierarchy = degenerate')
        kh_ic, _,Pk_nu_ic= result.get_matter_power_spectrum(minkh=k_ic_min, maxkh=k_ic_max, npoints = npbin,var1='delta_nu'  ,var2='delta_nu')
        Pk_nu_ic = [Pk_nu_ic[0],Pk_nu_ic[0],Pk_nu_ic[0]]#list(np.zeros(npbin)),list(np.zeros(npbin))]
    else :
        Pk_nu_ic=[0,0,0]
        i=0
        for m_nu_i in m_nus:
            print(i,'m_nu = ',m_nu_i)
            if (m_nu_i >0):
                par = camb.CAMBparams()
                par.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, omk=omk, neutrino_hierarchy='degenerate', num_massive_neutrinos=1, mnu=m_nu_i, nnu=nnu, standard_neutrino_neff=standard_neutrino_neff)
                par.omch2=omch2-par.omnuh2
                par.InitPower.set_params(As=As,ns=ns)
                par.set_matter_power(redshifts=[z_max], kmax=k_ic_max, nonlinear=True)
                par.NonLinear = camb.model.NonLinear_both
                result1 = camb.get_results(par)
                kh_ic, _,Pk_nu_ic_i= result1.get_matter_power_spectrum(minkh=k_ic_min, maxkh=k_ic_max, npoints = npbin,var1='delta_nu'  ,var2='delta_nu')
                Pk_nu_ic[i]=Pk_nu_ic_i[0]
            else:
                print(i,'m_nu = 0')
                Pk_nu_ic[i]=list(np.zeros(npbin))
            i+=1
    def Ha(a):
        return result.hubble_parameter(1/a-1)
    def Hz(z):
        return result.hubble_parameter(z)
else:
    #run reps
    run_reps()

    k = np.loadtxt('%s/IC/IC_Pcb_rescaled_z%.4f.txt'%(nupath,z_max))[:,0]
    kh_ic = np.exp(np.linspace(np.log(k_ic_min), np.log(k_ic_max),  npbin))
    Pk_cb_ic = interp_pk(z_max,'cb',kh_ic)
    Pk_nu_ic = interp_pk(z_max,'n',kh_nonlin)
    H00 = np.loadtxt('%s/IC/IC_hubble.txt'%nupath)
    Ha = interp1d(1/(H00[:,0]+1),H00[:,1], kind='cubic', bounds_error=False, fill_value="extrapolate")
    def Hz(z):
        return Ha(1/(1+z))
# print('*'*200)
# print(nupath+'/IC/Pnu_ic.txt')
# print(nupath+'/IC/Pcb_ic.txt')
# print('#'*100)
# print(Pk_nu_ic)
    
def taua(a):
    z=1/a-1
    return integrate.quad(lambda z0: 1/Hz(z0), z,10000)[0]
T1 = time.time()

print('calculating Expansion History')
n_a = int(istep_max)
dt0 = 5e-4
t = -np.arange(n_a)*dt0
H_ex = np.ones(n_a)*H0
tau = np.ones(n_a)
a_ex = np.zeros(n_a)
a_ex[0] = 1
tau[0] = taua(1)
for i in range(n_a-1):
    Hai = Ha(a_ex[i])
    H_ex[i] = Hai
    a_ex[i+1] = -omHsq*Hai*a_ex[i]**3 * dt0 +a_ex[i]
    dz = 1/a_ex[i]-1/a_ex[i+1]
    tau[i+1]  = tau[i]+(1/Hz(1/a_ex[i]-1)+2/Hz(1/a_ex[i]-1+dz)+2/Hz(1/a_ex[i+1]-1-dz)+1/Hz(1/a_ex[i+1]-1))*dz/6

    if (a_ex[i+1]<=1./301 or np.isnan(a_ex[i+1])):
        break
i_end = i
T2 = time.time()
tau[-1] = taua(1/201)
# tau = 299792.458*tau
os.system('rm '+nupath+'/*.txt')
np.savetxt(nupath+'/s_a_tau_H.txt',np.array([t,a_ex,tau,H_ex]))
print("EH:\n\n      time:   %.2f seconds\n      step:   %d\n      save:   '%s'\n\n\n"%(T2-T1,i_end,nupath+'/s_a_tau_H.txt'))

print('get Pk')
n=int(z_nonlin.shape[0])
Pk_nl = [0]*n
Pk_nu = [0]*n
if (all_reps):
    #get pk from reps
    for i in range(n):
        # print(z_nonlin[i])
        Pk_nl[i] = interp_pk(z_nonlin[i],'cb',kh_nonlin)
        Pk_nu[i] = interp_pk(z_nonlin[i],'n',kh_nonlin)
else:
    #get pk from camb
    Pk_cdm_cambs,Pk_nu_cambs,kh_nl,z_nonlin = get_Pk_nonlin_CDM(n)
    Pk_nu_ic_k = interp1d(kh_ic,Pk_nu_ic, kind='linear', bounds_error=False, fill_value="extrapolate")
    Pk_nu_ic = Pk_nu_ic_k(kh_nonlin)
    print(kh_nl.max(),kh_nl.min())
    print(kh_nonlin.max(),kh_nonlin.min())
    # exit()

    for i in range(n):
        # 找到Pk是否有小于0的值
        if np.any(Pk_cdm_cambs[i]<0):
            print('Pk_cdm_cambs[i]<0')
            print(Pk_cdm_cambs[i])
            print('z_nonlin[i] = ',z_nonlin[i])
            exit()
        
        Pk_cdm_z = interp1d(kh_nl,Pk_cdm_cambs[i], kind='linear', bounds_error=False, fill_value="extrapolate")
        Pk_nl[i] = Pk_cdm_z(kh_nonlin)
        Pk_nu_z = interp1d(kh_nl,Pk_nu_cambs[i], kind='linear', bounds_error=False, fill_value="extrapolate")
        Pk_nu[i] = Pk_nu_z(kh_nonlin)

#write Pk to nupath
z_powerpoint=open(nupath+'/z_powerpoint.txt', 'w')
z_values=open(nupath+'/z_values.txt', 'w')
k_values=open(nupath+'/k_values.txt', 'w')

print('write Pk to nupath')
# print('save path : ',nupath+'/IC/Pcb_ic.txt')
# print(type(kh_ic), type(Pk_cb_ic))
# print(kh_ic.dtype, Pk_cb_ic.dtype)
np.savetxt(nupath+'/IC/Pcb_ic.txt',np.array([kh_ic,Pk_cb_ic]).T)
np.savetxt(nupath+'/IC/Pnu_ic.txt',np.array(Pk_nu_ic))
for i in range(len(kh_nonlin)):
    k_values.write('%3.12f\n'%kh_nonlin[i])

if test:
    for i in range(n):
        z_values.write('%3.4f\n'%z_nonlin[i])
        
        if i%ni==0:
            z_powerpoint.write('%3.4f\n'%z_nonlin[i])
        np.savetxt(nupath+'/Pk_cb_%3.4f.txt'%z_nonlin[i],Pk_nl[i])
        f_nr = get_f_nr(z_nonlin[i])
        np.savetxt(nupath+'/Pk_nu_%3.4f.txt'%z_nonlin[i],Pk_nu[i])
        np.savetxt(nupath+'/Tf_nu_%3.4f.txt'%z_nonlin[i],((1-f_nu)*np.sqrt(Pk_nl[i])+f_nr*(np.sqrt(Pk_nu[i])))/np.sqrt(Pk_nl[i]))
else:
    z_str = ''
    if (Mass_nu > 0):
        for i in range(n):
            
            np.savetxt(nupath+'/Pk_cb_%3.4f.txt'%z_nonlin[i],Pk_nl[i])

            # 找到Pk是否有小于0的值
            if np.any(Pk_nl[i]<0):
                print('Pk_cb has negative values')
                print(Pk_nl[i])
                print(z_nonlin[i])
                print(i)
                print(n)
                print(nupath)
                print(nupath+'/Pk_cb_%3.4f.txt'%z_nonlin[i])
                exit()
            f_nr = get_f_nr(z_nonlin[i])
            np.savetxt(nupath+'/Pk_nu_%3.4f.txt'%z_nonlin[i],Pk_nu[i][:])
            # 找到Pk是否有小于0的值
            if np.any(Pk_nu[i]<0):
                print('Pk_nu has negative values')
                print(Pk_nu[i])
                print(z_nonlin[i])
                print(i)
                print(n)
                print(nupath)
                print(nupath+'/Pk_nu_%3.4f.txt'%z_nonlin[i])
                exit()
            np.savetxt(nupath+'/Tf_nu_%3.4f.txt'%z_nonlin[i],((1-f_nu)*np.sqrt(Pk_nl[i])+f_nr*(np.sqrt(Pk_nu[i])))/np.sqrt(Pk_nl[i]))
            # print(nupath+'/Tf_nu_%3.4f.txt'%z_nonlin[i])
            z_str+='%3.4f\n'%z_nonlin[i]
        # for i in z_nonlin:
        # print(z_str)
        z_powerpoint.write(z_str)

z_powerpoint.close()
z_values.close()

print('******************** Pk_init done ********************\n\n\n')