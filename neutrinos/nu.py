import camb,sys,os
from matplotlib import pyplot as plt
import numpy as np
from camb import model
import multiprocessing as mp
from functools import partial
from matplotlib.gridspec import GridSpec

import math

def L(X):
    T_nu0=0.000168
    L=(1+0.0168*(X*T_nu0*C)**2+0.0407*(X*T_nu0*C)**4)/(1+2.1734*(X*T_nu0*C)**2+1.6787*(X*T_nu0*C)**4.1811+0.1467*(X*T_nu0*C)**8)
    return L

def get_Pk_nonlin_CDM():
    fid = open('./neutrinos/Pk_nu/a_values.txt')
    a = np.loadtxt(fid, dtype='float32')
    fid = open('./neutrinos/Pk_nu/tau_values.txt')
    tau = np.loadtxt(fid, dtype='float32')
    fid = open('./neutrinos/Pk_nu/z_values.txt')
    z_nonlin = np.loadtxt(fid, dtype='float32')
    fid = open('./neutrinos/Pk_nu/k_values.txt')
    kh_nonlin = np.loadtxt(fid, dtype='float32')
    Pk_nonlin_CDM = [0]*len(a)
    for i in range(len(a)):
        print(z_nonlin[i])
        fid = open('./neutrinos/Pk_nu/Pk_nu_%3.4f.dat.txt'%z_nonlin[i])
        Pk_nonlin_CDM[i] = np.loadtxt(fid, dtype='float32')
        
    return Pk_nonlin_CDM,kh_nonlin,z_nonlin,a,tau
    
def nu_z(z0):
    global C
    z_index=np.where(np.array(z_nonlin)<=z0)[0][0]
    z0 = z_nonlin[z_index]
    print('z_index',z_index,z0)
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, omk=omk, neutrino_hierarchy=neutrino_hierarchy, num_massive_neutrinos=num_massive_neutrinos, mnu=mnu, nnu=nnu, standard_neutrino_neff=standard_neutrino_neff)
    pars.omch2=omch2-pars.omnuh2
    f_nu = pars.omnuh2/(pars.omch2+pars.ombh2+pars.omnuh2)
    pars.InitPower.set_params(ns=ns)
    pars.set_matter_power(redshifts=[z0,z_max],kmax=kmax)
    pars.NonLinear = model.NonLinear_none
    results = camb.get_results(pars)
    kh_nl, z_nl,Pk_nu_lin= results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints = npbin,var1='delta_nu',var2='delta_nu')
    pars.NonLinear = model.NonLinear_both
    results.calc_power_spectra(pars)
    kh_nl, z_nl,Pk_nu_nonlin= results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints = npbin,var1='delta_nu',var2='delta_nu')
    kh_nl, z_nl,Pk_matter_nonlin= results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints = npbin,var1='delta_tot',var2='delta_tot')
    kh_nl, z_nl,Pk_cdm_nonlin= results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints = npbin,var1='delta_cdm',var2='delta_cdm')
    kh_nl, z_nl,Pk_b_nonlin= results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints = npbin,var1='delta_baryon',var2='delta_baryon')

    h0=pars.H0/100.0
    M_nu=pars.omnuh2*93.14
    m_nu=M_nu/3.0

    omega_m=(pars.omch2+pars.ombh2+pars.omnuh2)/(h0**2)
    omega_l=1-omega_m
    C=299792.458*h0

    print(h0,m_nu,omega_m)
    Pk_nu_new=[0]*npbin
    Pk_nu_new1=[0]*npbin
    Pk_nu_new2=[0]*npbin

    s_pre =0
    s_post=0
    s_f=0

    results = camb.get_results(pars)
    H=results.hubble_parameter(z_max)
    H0_local=results.hubble_parameter(z0)

    omega_m_zi=(H0_local/H)**2*omega_m*(1+z_max)**3
    omega_l_zi=(H0_local)**2/H**2*omega_l
    f=omega_m_zi**(4.0/7.0)+omega_l_zi/70.0*(1+omega_m_zi/2.0)
    print(H0_local,H,(H0_local/H)**2,omega_m,H0_local,H)

    for i in range(1,z_index):
        # print(i,'sf=',s_f)
        s_f=s_f+0.5*(1/a[i-1]+1/a[i])*(tau[i]-tau[i-1]) #superconformal time s=integral(d(z)/a(z)) 
    # print(a)
    # print(tau)

    print("s_f=",s_f)
        
    for j in range(0,npbin):
        Pk_nu_new1[j]=L(s_f*kh_nonlin[j]/m_nu)*(Pk_nu_lin[1][j])**(0.5)*(1+s_f*a[0]*a[0]*H*f) #好像不对
    # print(np.array(Pk_nu_new1))
    # print(Pk_nu_lin[1])
    for i in range(1,z_index):
        s_pre=s_post
        s_post=s_pre+0.5*(1/a[i-1]+1/a[i])*(tau[i]-tau[i-1]) #superconformal time s=integral(d(z)/a(z))
        for j in range(0,npbin):
            Pk_nu_new2[j]=Pk_nu_new2[j]+(L((s_f-s_pre)*kh_nonlin[j]/m_nu)*(Pk_nl[i-1][j])**0.5*(s_f-s_pre)+L((s_f-s_post)*kh_nonlin[j]/m_nu)*(Pk_nl[i][j])**0.5*(s_f-s_post))*(tau[i]-tau[i-1])/2.0
        # print(i,'step',(tau[i]-tau[i-1])/2.0,s_pre,s_post) 
        # print(i,s_f,s_pre,'Lx1',[L((s_f-s_pre)*kh_lin[-3]/m_nu)*(Pk[i-1][-3])**0.5,L((s_f-s_pre)*kh_lin[-2]/m_nu)*(Pk[i-1][-2])**0.5,L((s_f-s_pre)*kh_lin[-1]/m_nu)*(Pk[i-1][-1])**0.5])    
        # # print((s_f-s_post),'Lx2',[L((s_f-s_post)*kh_lin[0]/m_nu)*(Pk[i][0])**0.5,L((s_f-s_post)*kh_lin[1]/m_nu)*(Pk[i][1])**0.5,L((s_f-s_post)*kh_lin[2]/m_nu)*(Pk[i][2])**0.5])    
        # print('Pk',Pk[i-1][-4:])  
        # print(i,'pk2',Pk_nu_new2[-3:])       
    # Pk_nu_new2=np.array(Pk_nu_new2)*1.5*omega_m*H0_local**2
    # print('pk2',1.5*omega_m*H0_local**2*np.array(Pk_nu_new2))
    for j in range(0,npbin):                                                       
        Pk_nu_new[j] = Pk_nu_new1[j]+Pk_nu_new2[j]*1.5*omega_m*H0_local**2
        Pk_nu_new[j] = Pk_nu_new[j]**2

    return Pk_nu_new,Pk_nu_nonlin,Pk_matter_nonlin,Pk_cdm_nonlin,z0


def dplt(ax,kh_nonlin, font,Pk1,Pk2,dPk_max=None,dPk_min=None,sq=True):
    fontsize = font['size']
    if sq:
        dPk = np.sqrt(Pk1[0])/np.sqrt(Pk2[0])-1
        ax.plot(kh_nonlin, dPk, label=r'$\sqrt{%s(k)}/\sqrt{%s(k)}-1$'%(Pk1[1],Pk2[1]), color='g', linewidth=1, linestyle='-')
    else:
        dPk = np.array(Pk1[0])/np.array(Pk2[0])-1
        ax.plot(kh_nonlin, dPk, label=r'$%s(k)/%s(k)-1$'%(Pk1[1],Pk2[1]), color='g', linewidth=1, linestyle='-')
    ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
    ax.set_xlabel(r'$k\, [h Mpc^{-1}]$', font)
    ax.set_ylabel('Residuals', font)
    #设置坐标轴的粗细
    ax.tick_params(left=True, right=True, bottom=True, top=True, width=1, labelsize=fontsize*0.5, direction='in')
    ax.grid(True, which='both', linestyle='dotted', linewidth=0.3)
    ax.spines['bottom'].set_linewidth(1);#设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(1);#设置左边坐标轴的粗细
    ax.spines['right'].set_linewidth(1);#设置右边坐标轴的粗细
    ax.spines['top'].set_linewidth(1);#设置上部坐标轴的粗细
    ax.set_yscale('symlog')  # 使用对数坐标轴，同时保留正负

    # 设置y轴的范围
    if dPk_max == None:
        dPk_max = np.max(dPk)
    if dPk_min == None:
        dPk_min = np.min(dPk)
    if abs(dPk_max/dPk_min) > 2:
        dPk_min = -dPk_max/2
    elif abs(dPk_max/dPk_min) < 1/2:
        dPk_max = -dPk_min/2
    ax.set_ylim(dPk_min*1.2,dPk_max*1.2)
    ax.set_yticks([float('%.1e'%(dPk_min*1.2)),float('%.1e'%(dPk_min*0.8)),float('%.1e'%(dPk_min*0.4)),0,float('%.1e'%(dPk_max*0.4)),float('%.1e'%(dPk_max*0.8)),float('%.1e'%(dPk_max*1.2))])
    yticks, ytick_labels = plt.yticks()
    ytick_labels[0].set_text(None)
    ytick_labels[-1].set_text(None)
    ax.set_yticklabels(ytick_labels)
    ax.legend(loc='lower right', ncol=1, fontsize=fontsize)
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    return ax

def plot_Pk2k(n, ni, kh_nonlin,stitle, s, Pk, dPk1 = None, dPk2 = None,sq = True,show = 0):
    fig = plt.figure(figsize=(8*size, 6*size))
    if dPk1 == None:
        gs = GridSpec(1, 1, hspace=0.05)
    elif dPk2 == None:
        gs = GridSpec(2, 1, height_ratios=[3,1], hspace=0.05)
    else:
        gs = GridSpec(3, 1, height_ratios=[10,3,3], hspace=0.05)
    
    # format the axes
    fontsize = 10*size
    font = {
            'weight': 'normal',
            'size': fontsize,
            }

    Pmax = 0
    Pmin = 0
    ax1 = fig.add_subplot(gs[0])
    for Pki in Pk:
        ax1.loglog(kh_nonlin, Pki[0], label=r'$%s(k)$'%Pki[1], color=Pki[2],linewidth=Pki[3], linestyle=Pki[4], alpha=Pki[5])
        Pmax = np.max([np.max(Pki[0]),Pmax])
        Pmin = np.min([np.min(Pki[0]),Pmin])

    ax1.set_xlabel(r'$k\, [h Mpc^{-1}]$', font)
    ax1.set_ylabel(r'$P(k)$', font)
    ax1.set_title(r'${%s}(k)_{(nstep=%d*%d,Z_{max} = %d)}\, \, \,math:%s$'%(stitle,n,ni,z_max,s)+'\n',fontsize=fontsize)
    #设置坐标轴的粗细
    ax1.tick_params(left=True, right=True, bottom=True, top=True, width=1, labelsize=fontsize*0.7, direction='in')
    ax1.grid(True, which='both', linestyle='dotted', linewidth=0.3)
    ax1.spines['bottom'].set_linewidth(1);#设置底部坐标轴的粗细
    ax1.spines['left'].set_linewidth(1);#设置左边坐标轴的粗细
    ax1.spines['right'].set_linewidth(1);#设置右边坐标轴的粗细
    ax1.spines['top'].set_linewidth(1);#设置上部坐标轴的粗细
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    # 设置y轴的范围
    if (Pmax >10 or Pmin<0.1):
        ax1.set_ylim(Pmin*(1/1.5),Pmax*1.5)
    else:
        ax1.set_ylim(0.99,1)
        # ax1.set_yscale('linear')
    # ax1.set_yticks([1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1, 1e2, 1e4])
    yticks, ytick_labels = plt.yticks()
    ytick_labels[0].set_text(None)
    ax1.set_yticklabels(ytick_labels)
    ax1.legend(loc='upper right', ncol=1, fontsize=fontsize)
    legend = ax1.legend()
    legend.get_frame().set_alpha(0.5)
    ax1.set_xlim(kmin*(1/1.5), kmax*1.5)

    if dPk1 == None:
        None
    else :
        ax2 = fig.add_subplot(gs[1], sharex=ax1)
        ax2 = dplt(ax2,kh_nonlin, font,Pk[dPk1[0]],Pk[dPk1[1]],dPk1[2],dPk1[3],sq)

        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_position([0.1, 0.3, 0.8, 0.6])
        ax2.set_position([0.1, 0.1, 0.8, 0.2])
        if dPk2 == None:
            None
        else:
            ax3 = fig.add_subplot(gs[2], sharex=ax1)
            ax3 = dplt(ax3,kh_nonlin, font,Pk[dPk2[0]],Pk[dPk2[1]],dPk2[2],dPk2[3],sq)
            ax1.set_position([0.1, 0.4, 0.8, 0.5])
            ax2.set_position([0.1, 0.25, 0.8, 0.15])
            ax3.set_position([0.1, 0.1, 0.8, 0.15])
            plt.setp(ax2.get_xticklabels(), visible=False)
    if show == 1:
        plt.show()
    else:
        try: 
            os.mkdir(imagepath)
            print('-'*200)
        except: 
            print('*'*200)
        # fig.savefig(imagepath+'%s/%d/Pk_%s_%d_%d_%d_%s.jpg'%(typ,z_max,sP,z_max,n,ni,s.replace('\\','')),dpi=300)
        fig.savefig('%s/%s_%d_%s.jpg'%(imagepath,stitle.replace("\\",''),z_max,s.replace('\\','')),dpi=300)
    # plt.close()


imagepath = './tf_img'
size = 1
mnu = 0.1


pi=np.pi
box = 400
nf = 64
nns = 4
nnt = 4
ng =128
nn = 1
npbin = math.ceil(nf*nns*nnt*nn/2*(3**0.5))
kmin = 2*pi/box
kmax = 2*pi*npbin/box

print(npbin,kmin,kmax)
#set cosmology parm
H0=70.0
ombh2=0.0463*(H0/100)**2
omch2=0.2327*(H0/100)**2
omk=0.0
neutrino_hierarchy='degenerate'
num_massive_neutrinos=3
mnu=0.1
nnu=3.044
standard_neutrino_neff=3.044
ns=0.972
As=2.6025e-09
z_max = 200

n=50
ni=50
n_all=n*ni
n_PK = 150 # max num of Pk(z) at once
Pk_nl,kh_nonlin,z_nonlin,a,tau = get_Pk_nonlin_CDM()


fid = open('./neutrinos/Pk_nu/z_powerpoint.txt')
z_powerpoint = np.loadtxt(fid)
z_str = ['%.4f'%z for z in z_powerpoint[z_powerpoint<5]]
for z in z_str:
    print(z)
    Pk_nu_new,Pk_nu_nonlin,Pk_matter_nonlin,Pk_cdm_nonlin,z0 = nu_z(float(z))
    fid = open('./tf_test/Pk_nu_interp_r%s.txt'%z)
    Pk_iqr = np.fromfile(fid, dtype='float32')
    Pk = [[Pk_nu_nonlin[0]] +['Pk_{camb,nu}' ,'g',6,'-',0.5]
        ,[Pk_iqr]   +['Pk_{fort,nu}','y',2,'-',.5]
        ,[Pk_nu_new]   +['Pk_{simi,nu}','b',4,'-',.5]]

    plot_Pk2k(n,ni,kh_nonlin,'Pk\_neutrinos','z=%.2f'%z0,Pk,[1,0,None,None],[2,1,None,None])