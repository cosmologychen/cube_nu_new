import camb,sys,os,math,re,time
from matplotlib import pyplot as plt
import numpy as np
from camb import model
import multiprocessing as mp    
from functools import partial
from matplotlib.gridspec import GridSpec
import scipy.integrate as integrate

test = 0
nocamb = 1

n=50
ni=10



def match_para(para):
    # 读取Fortran文件
    file_path = './parameters.f90'
    with open(file_path, 'r') as file:
        content = file.read()

    # 使用正则表达式匹配模式
    pattern = r'parameter\s*::\s*%s\s*=\s*([^\s]+)'%para
    match = re.search(pattern, content)
    # print(match)

    if match:
        variable_value = match.group(1).strip()
        print(f"Variable value of {para}: {variable_value}")
        # sys.stdout.flush()  # 刷新输出缓冲
        return float(variable_value)
    else:
        print("Pattern not found in the file.")
        sys.exit()

def L(X):
    T_nu0=0.000168
    X = X*T_nu0*C
    # print('X=',X)
    L=(1+0.0168*X**2+0.0407*X**4)/(1+2.1734*X**2+1.6787*X**4.1811+0.1467*X**8)
    return L

def nu(Pk,tau,z,a,kh_lin,step_num):
    global C
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, omk=omk, neutrino_hierarchy=neutrino_hierarchy, num_massive_neutrinos=num_massive_neutrinos, mnu=mnu, nnu=nnu, standard_neutrino_neff=standard_neutrino_neff)
    pars.omch2=omch2-pars.omnuh2
    pars.InitPower.set_params(ns=ns)
    pars.set_matter_power(redshifts=[0,z_max],kmax=kmax)
    pars.NonLinear = model.NonLinear_none
    results = camb.get_results(pars)
    kh_nl, z_nl,Pk_nu_lin= results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints = npbin,var1='delta_nu',var2='delta_nu')
    pars.NonLinear = model.NonLinear_both
    results.calc_power_spectra(pars)
    kh_nl, z_nl,Pk_nu_nonlin= results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints = npbin,var1='delta_nu',var2='delta_nu')


    h0=pars.H0/100.0
    M_nu=pars.omnuh2*93.14
    m_nu=M_nu/3.0

    omega_m=(pars.omch2+pars.ombh2+pars.omnuh2)/(h0**2)
    omega_l=1-omega_m
    C=299792.458*h0

    Pk_nu_new=[0]*npbin
    Pk_nu_new1=[0]*npbin
    Pk_nu_new2=[0]*npbin

    s_pre =0
    s_post=0
    s_f=0

    results = camb.get_results(pars)
    H=results.hubble_parameter(100)
    H0_local=results.hubble_parameter(0)

    omega_m_zi=(H0_local/H)**2*omega_m*(1+z[0])**3
    omega_l_zi=(H0_local)**2/H**2*omega_l
    f=omega_m_zi**(4.0/7.0)+omega_l_zi/70.0*(1+omega_m_zi/2.0)

    for i in range(1,step_num):
        # print('s_f=',s_f)
        s_f=s_f+0.5*(1/a[i-1]+1/a[i])*(tau[i]-tau[i-1]) #superconformal time s=integral(d(z)/a(z)) 
    # print(a)
    # print(tau)

    # print("s_f=",s_f,C,a[0:2],tau[0:2],kh_nl[0:2])
    # for j in range(0,npbin):
    #     print('Lx=',L(s_f*kh_lin[j]/m_nu)*(Pk_nu_lin[1][j])**(0.5))


    print(s_f,m_nu,C)
    np.savetxt(nupath+'/Pk_nu_ic.txt',Pk_nu_lin[1][:])
    # for j in range(0,npbin):
    #     print('L=',L(s_f*kh_lin[j]/m_nu)*(Pk_nu_lin[1][j])**(0.5))
        
        
    # for j in range(0,npbin):
    #     Pk_nu_new1[j]=L(s_f*kh_lin[j]/m_nu)*(Pk_nu_lin[1][j])**(0.5)*(1+s_f*a[0]*a[0]*H*f) #好像不对
    # print('_'*100)
    # for i in range(1,step_num):
    #     s_pre=s_post
    #     s_post=s_pre+0.5*(1/a[i-1]+1/a[i])*(tau[i]-tau[i-1]) #superconformal time s=integral(d(z)/a(z))
    #     for j in range(0,npbin):
    #         Pk_nu_new2[j]=Pk_nu_new2[j]+1.5*omega_m*H0_local**2*(L((s_f-s_pre)*kh_lin[j]/m_nu)*(Pk[i-1][j])**0.5*(s_f-s_pre)+L((s_f-s_post)*kh_lin[j]/m_nu)*(Pk[i][j])**0.5*(s_f-s_post))*(tau[i]-tau[i-1])/2.0
                                                        
    # for j in range(0,npbin):                                                       
    #     Pk_nu_new[j] = Pk_nu_new1[j]+Pk_nu_new2[j]
    #     Pk_nu_new[j] = Pk_nu_new[j]**2


    return Pk_nu_new,Pk_nu_nonlin[0], Pk_nu_lin[0]

def get_line(y1,y2,x1,x2,x):
    a = (y2-y1)/(x2-x1)
    b = y1 - a*x1
    y = a*x+b
    return y

def interp_line(arr,ni,a=None):
    if a is None:
        a = range((len(arr)-1)*ni+1)
    # print('len_in=',len(arr)-1)
    arr_new =[]
    for i in range(0,len(arr)-1):
        arr_new.append(arr[i])
        for j in range(1,ni):
            arr_new.append(get_line(arr[i],arr[i+1],a[i*ni],a[i*ni+ni],a[i*ni+j] ))
    arr_new.append(arr[-1])
    arr_new = np.array(arr_new)
    # print('len_out',len(arr_new))
    return arr_new

def get_a(n,typ): #get array a
    z = [0]*n
    a = [0]*n
    z[0] = z_max
    if typ == 'line':
        a = np.linspace(1/(1+z[0]),1,n+1)
    elif typ =='log':
        log_a = [np.log(1/(1+z[0])),0]
        log_a = interp_line(log_a,n)
        a = np.exp(log_a)
    else:
        print('no such type')
        sys.exit()
    return a

def get_Pk_nonlin_CDM(n,typ):
    a = get_a(n,typ)
    z = 1/a-1
    z.sort()
    nz = int(np.ceil(len(z)/n_PK))
    tau_i=[[0]*n_PK for i in range(0,nz)]
    Pk_nonlin_CDM = [0]*len(z)
    Pk_nonlin_nu = [0]*len(z)
    z_nonlin = [0]*len(z)
    tau = [0]*len(z)
    for i in range(0,nz):
        z_n=z[i*n_PK:(i+1)*n_PK]
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, omk=omk, neutrino_hierarchy=neutrino_hierarchy, num_massive_neutrinos=num_massive_neutrinos, mnu=mnu, nnu=nnu, standard_neutrino_neff=standard_neutrino_neff)
        pars.omch2=omch2-pars.omnuh2
        pars.InitPower.set_params(ns=ns)
        print('i = %d/%d;z_max = %.3f;len = %d'%(i+1,nz,z_n[-1],len(z_n)))
        pars.set_matter_power(redshifts=z_n, kmax=kmax, nonlinear=True)
        pars.NonLinear = model.NonLinear_both
        results = camb.get_results(pars)
        kh_nonlin, z_nonlin[i*n_PK:i*n_PK+len(z_n)],Pk_nonlin_CDM[i*n_PK:i*n_PK+len(z_n)]= results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints = npbin,var1='delta_nonu',var2='delta_nonu')
        kh_nonlin, z_nonlin[i*n_PK:i*n_PK+len(z_n)],Pk_nonlin_nu[i*n_PK:i*n_PK+len(z_n)]= results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints = npbin,var1='delta_nu',var2='delta_nu')
        tau_i=results.conformal_time(z_n, presorted=True, tol=None)
        tau[i*n_PK:i*n_PK+len(z_n)]=tau_i/299792.458
        pars = 0

    tau.sort()
    return Pk_nonlin_CDM[::-1],Pk_nonlin_nu[::-1],kh_nonlin,z_nonlin[::-1],a,tau

def get_f_nr(z):
    y = Mass_nu/(sigma_nu*N_eff*k_b*T_gama*(1+z))
    Fy=integrate.quad(lambda u: (u**2*np.sqrt(u**2+y**2))/(1+np.exp(u)), 0,300)[0]
    f1=integrate.quad(lambda u: (u**2)/((1+np.exp(u))*np.sqrt(u**2+y**2)), 0,300)[0]
    f_nr = y**2*f1/Fy#*f_nu
    return f_nr
    
if nocamb:
    def Ha(a0):
        wde = -1
        omega_r = 5.046734693877551e-05
        omega_l = 1-omega_m-omega_r
        dt_x = 1e-4 / 2
        a_x = a0
        omhsq = 4.0 / 9.0
        a3rlm = a_x**(-3*wde) * omega_l / omega_m
        arkm = a_x * (1.0 - omega_m - omega_l-omega_r) / omega_m

        adot = np.sqrt(omhsq * a_x**3 * (1.0 + arkm + a3rlm))
        addot = a_x**2 * omhsq * (1.5 + 2.0 * arkm + 1.5 * (1.0 - wde) * a3rlm)
        atdot = a_x * adot * omhsq * (3.0 + 6.0 * arkm + 1.5 * (2.0 - 3.0 * wde) * (1.0 - wde) * a3rlm)

        da1 = adot * dt_x + (addot * dt_x**2) / 2.0 + (atdot * dt_x**3) / 6.0

        a_x = a0 + da1
        omhsq = 4.0 / 9.0
        a3rlm = a_x**(-3*wde) * omega_l / omega_m
        arkm = a_x * (1.0 - omega_m - omega_l) / omega_m

        adot = np.sqrt(omhsq * a_x**3 * (1.0 + arkm + a3rlm))
        addot = a_x**2 * omhsq * (1.5 + 2.0 * arkm + 1.5 * (1.0 - wde) * a3rlm)
        atdot = a_x * adot * omhsq * (3.0 + 6.0 * arkm + 1.5 * (2.0 - 3.0 * wde) * (1.0 - wde) * a3rlm)

        da2 = adot * dt_x + (addot * dt_x**2) / 2.0 + (atdot * dt_x**3) / 6.0

        Hz = (da1+da2)/2/dt_x/a0/a0/a0/omHsq

        return Hz
    def Hz(z):
        return Ha(1/(1+z))
    def taua(a):
        z=1/a-1
        return 299792.458*integrate.quad(lambda z0: 1/Hz(z0), z, np.inf)[0]
    def Dgrow_2(a,results0):
        return 5/2*omega_m*results0.hubble_parameter(1/a-1)*H0**2*(integrate.quad(lambda a00: 1/(a00*results0.hubble_parameter(1/a00-1))**3,0,a)[0])

    

else:
    def Ha(a):
        return results.hubble_parameter(1/a-1)
    def taua(a):
        return results.conformal_time(1/a-1)/299792.458
    def Dgrow_2(a,results0):
        return 5/2*omega_m*results0.hubble_parameter(1/a-1)*H0**2*(integrate.quad(lambda a00: 1/(a00*results0.hubble_parameter(1/a00-1))**3,0,a)[0])

#set cosmology parm
H0=match_para('h0')*100
omega_bar=match_para('omega_bar')
omega_cdm=match_para('omega_cdm')
omega_r=match_para('omega_r')
omk=0.0
neutrino_hierarchy='degenerate'
mnu=match_para('Mass_nu')
if (mnu>0) :
    num_massive_neutrinos=3
else:
    num_massive_neutrinos=0
nnu=3.044
standard_neutrino_neff=match_para('N_eff')
ns=match_para('n_s')
As=match_para('A_s')
z_max = 200

pi=np.pi
ratio_cs=match_para('ratio_cs')
ratio_sf=np.array([1,2,4,6,8,12,16])
ng = match_para('ng')
box = match_para('box')
nns = match_para('nns')
nnt = match_para('nnt')
nn = match_para('nn')
ngp = ng/nnt
ngb = match_para('ngb')
istep_max = match_para('istep_max')
ngt=ngp+2*ngb 
nfp=(ngp/nns)*ratio_sf
nfb=ratio_sf*ratio_cs
nft=nfp+2*nfb
cic_iapm = 5-1
tile = box/nn/nnt*ngt/ngp
nw = ng/ratio_cs
nw_global=nw*nn
nyquist=nw_global/2
subtile = box/nn/nnt/nns*(nft[cic_iapm]/nfp[cic_iapm])
ncbin = int(nyquist*np.sqrt(3.))+1
nnbin= int(np.max([nft[cic_iapm],ngt])/2*np.sqrt(3.))+1
npbin  = ncbin+nnbin
kmin = 2*pi/box
kmax = 2*pi/(box/nn/nnt/nns/nfp[cic_iapm]*2/np.sqrt(3))
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


file_path = './parameters.f90'
with open(file_path, 'r') as file:
    content = file.read()

# 使用正则表达式匹配模式
pattern = r'parameter\s*::\s*opath\s*=\s*([^\s]+)'
match = re.search(pattern, content)
opath = match.group(1).strip()[1:-1]
print(f"Variable value of opath : {opath}")
nupath =  opath+"neutrinos"
print(f"Variable value of nupath : {nupath}")
print(omega_r,omega_bar,omega_cdm,omega_l,Mass_nu)



try:
    os.system('mkdir -p '+nupath+'/tf')
    os.system('mkdir -p '+nupath+'/Pk')
except:
    None
    
# os.system('ls '+nupath)
os.system('rm '+nupath+'/*.txt')

omHsq = 2/3*np.sqrt(1/omega_m)/H0
pars = camb.CAMBparams()
pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, omk=omk, neutrino_hierarchy=neutrino_hierarchy, num_massive_neutrinos=num_massive_neutrinos, mnu=mnu, nnu=nnu, standard_neutrino_neff=standard_neutrino_neff)
pars.omch2=pars.omch2-pars.omnuh2
f_nu = pars.omnuh2/(pars.omch2+pars.ombh2+pars.omnuh2)
pars.InitPower.set_params(As=As,ns=ns)
pars.set_matter_power(redshifts=[0,200],kmax=kmax)
pars.NonLinear = camb.model.NonLinear_both
pars.Transfer.high_precision = True
pars.Transfer.accurate_massive_neutrinos = True
pars.Transfer.k_per_logint = 10
pars.Transfer.kmax = 100
pars.Transfer.PK_num_redshifts = 1
pars.Transfer.PK_redshifts = [0]
results = camb.get_results(pars)
print(pars)
print('s8  =',results.get_sigma8_0())
T1 = time.time()
n_a = int(istep_max)
dt0 = 5e-4
t = -np.arange(n_a)*dt0
H_ex = np.ones(n_a)*H0
tau = np.ones(n_a)
a_ex = np.zeros(n_a)
a_ex[0] = 1
for i in range(n_a-1):
    Hai = Ha(a_ex[i])
    H_ex[i] = Hai
    tau[i]  = taua(a_ex[i])
    a_ex[i+1] = -omHsq*Hai*a_ex[i]**3 * dt0 +a_ex[i]

    if (a_ex[i+1]<=1./301 or np.isnan(a_ex[i+1])):
        break
i_end = i
T2 = time.time()
print(f"{T2-T1}"+"seconds")
tau[-1] = taua(1/201)
t[-1] = Dgrow_2(1/201,results)/ Dgrow_2(1,results)
np.savetxt(nupath+'/s_a_tau_H.txt',np.array([t,a_ex,tau,H_ex]))
print(i_end,n_a,nupath+'/s_a_tau_H.txt')


kh_nl, z_nl,Pk_mater_nonlin= results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints = 200,var1='delta_nonu',var2='delta_nonu')
np.savetxt(nupath+'/tf_ic.txt', results.get_matter_transfer_data().transfer_data[[0, 7], :, 0].T)
print(npbin,ncbin,nnbin,kmin,kmax)
print(f_nu)

# os.system('cd ./reps-master/ && ./reps ../neutrinos/IC01.ini && cd ..')

if test:
    n_all=n*ni
else:
    n_all=n
n_PK = 150 # max num of Pk(z) at once
Pk_nl,Pk_nu,kh_nonlin,z_nonlin,a,tau = get_Pk_nonlin_CDM(n_all,'line')






PK_hr_nu,Pk_nu_nonlin, Pk_nu_lin = nu(Pk_nl,tau,z_nonlin,a,kh_nonlin,n_all+1)
np.savetxt(nupath+'/Pk_nu_hr.txt',PK_hr_nu)
np.savetxt(nupath+'/Pk_nu_nl.txt',Pk_nu_nonlin)
np.savetxt('./Pk_mater_nl.txt',Pk_mater_nonlin[0])
Pk_nl = np.array(Pk_nl)
z_powerpoint=open(nupath+'/z_powerpoint.txt', 'w')
z_values=open(nupath+'/z_values.txt', 'w')
a_values=open(nupath+'/a_values.txt', 'w')
tau_values=open(nupath+'/tau_values.txt', 'w')
k_values=open(nupath+'/k_values.txt', 'w')

for i in range(len(kh_nonlin)):
    k_values.write('%3.12f\n'%kh_nonlin[i])

if test:
    for i in range(len(z_nonlin)):
        z_values.write('%3.4f\n'%z_nonlin[i])
        a_values.write('%3.12f\n'%a[i])
        tau_values.write('%3.12f\n'%tau[i])
        
        if i%ni==0:
            z_powerpoint.write('%3.12f\n'%z_nonlin[i])
        np.savetxt(nupath+'/Pk_cb_%3.4f.txt'%z_nonlin[i],Pk_nl[i])
        f_nr = get_f_nr(z_nonlin[i])
        np.savetxt(nupath+'/Pk_nu_%3.4f.txt'%z_nonlin[i],Pk_nu[i])
        np.savetxt(nupath+'/Tf_nu_%3.4f.txt'%z_nonlin[i],((1-f_nu)*np.sqrt(Pk_nl[i])+f_nr*(np.sqrt(Pk_nu[i])))/np.sqrt(Pk_nl[i]))

else:
    for i in range(len(z_nonlin)):
        z_powerpoint.write('%3.12f\n'%z_nonlin[i])
        a_values.write('%3.12f\n'%a[i])
        tau_values.write('%3.12f\n'%tau[i])
        
        np.savetxt(nupath+'/Pk_cb_%3.4f.txt'%z_nonlin[i],Pk_nl[i])
        f_nr = get_f_nr(z_nonlin[i])
        np.savetxt(nupath+'/Pk_nu_%3.4f.txt'%z_nonlin[i],Pk_nu[i])
        np.savetxt(nupath+'/Tf_nu_%3.4f.txt'%z_nonlin[i],((1-f_nu)*np.sqrt(Pk_nl[i])+f_nr*(np.sqrt(Pk_nu[i])))/np.sqrt(Pk_nl[i]))



z_powerpoint.close()
a_values.close()
z_values.close()
tau_values.close()

print('Pk_init done')