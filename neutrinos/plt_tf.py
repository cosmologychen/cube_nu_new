import glob,re, camb,sys,os,math
from matplotlib import pyplot as plt
import numpy as np
from camb import model
import multiprocessing as mp
from functools import partial
from matplotlib.gridspec import GridSpec

# Pk_nu 

Dir = '/image1/'
Path_nu = '../output/400_256_3_2_2_pp1_new'
z_checkpoint = np.loadtxt('./z_checkpoint.txt')
z = z_checkpoint[70]
z = 0.0006
mnu = 0.1
size = 1
show = 0
kh_r = np.loadtxt(Path_nu+'/neutrinos/k_values.txt')# 读取数据到一维数组

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

def get_pk_nu(z0):

      pars = camb.CAMBparams()
      pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, omk=omk, neutrino_hierarchy=neutrino_hierarchy, num_massive_neutrinos=num_massive_neutrinos, mnu=mnu, nnu=nnu, standard_neutrino_neff=standard_neutrino_neff)
      pars.omch2=omch2-pars.omnuh2
      f_nu = pars.omnuh2/(pars.omch2+pars.ombh2+pars.omnuh2)
      pars.InitPower.set_params(ns=ns)
      pars.set_matter_power(redshifts=[z0,z_max],kmax=kmax)
      pars.NonLinear = model.NonLinear_both
      results = camb.get_results(pars)
      results.calc_power_spectra(pars)
      kh_nl, z_nl,Pk_nu_nonlin= results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints = npbin,var1='delta_nu',var2='delta_nu')
      kh_nl, z_nl,Pk_cdm_nonlin= results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints = npbin,var1='delta_cdm',var2='delta_cdm')
    #   print('f_nu = ',f_nu)
    #   print('z0   = ',z0)

      return Pk_nu_nonlin[0],Pk_cdm_nonlin[0],kh_nl,f_nu

def plot_Pk2k1(stitle, s, Pk,show = 0):
    fig = plt.figure(figsize=(8*size, 6*size))
    gs = GridSpec(1, 1, hspace=0.05)
    
    # format the axes
    fontsize = 10*size
    font = {
            'weight': 'normal',
            'size': fontsize,
            }

    Pmax = 0
    Pmin = 1
    ax1 = fig.add_subplot(gs[0])
    for Pki in Pk:
        ax1.loglog(Pki[0], Pki[1], label=r'$%s(k)$'%Pki[2], color=Pki[3],linewidth=Pki[4], linestyle=Pki[5], alpha=Pki[6])
        Pmax = np.max([np.max(Pki[1]),Pmax])
        Pmin = np.min([np.min(Pki[1]),Pmin])

    ax1.set_xlabel(r'$k\, [h Mpc^{-1}]$', font)
    ax1.set_ylabel(r'$P(k)$', font)
    ax1.set_title(r'${%s}(k)_{Z_{max} = %d}\, \, \,math:%s$'%(stitle,z_max,s)+'\n',fontsize=fontsize)
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
    ax1.legend(loc='upper right', ncol=1, fontsize=fontsize)
    legend = ax1.legend()
    legend.get_frame().set_alpha(0.5)
    ax1.set_xlim(kmin*(1/1.5), kmax*1.5)

    subtile = ((stitle.replace("\\",'')).replace("{",'')).replace("}",'')
    if show == 1:
        plt.show()
    else:
        os.system('mkdir -p neutrinos/image/%s'%subtile)
        # 关闭图形，但仍然保存图像文件
        plt.close()
        fig.savefig('./neutrinos/image/%s/%s.jpg'%(subtile,(subtile+(s.replace('\\,','')).replace(' ',''))),dpi=300)
  
def worker(z):
      print('*'*20)
      print(z)
      fid = open(Path_nu+'/neutrinos/Pk/Pk_nu_{:.4f}.txt'.format(z), 'rb')  # 以二进制模式打开文件
      Pk_nu_r = np.fromfile(fid, dtype='float32').tolist()  # 读取数据到一维数组
      Pk_nu_nonlin,Pk_cdm_nonlin,kh_nl,f_nu = get_pk_nu(z)

      Pk = [
           [kh_r,Pk_nu_r] +['Pk_{fort,nu}' ,'r',2,'-',.5],
            [kh_nl,Pk_nu_nonlin]   +['Pk_{camb,nu}','.5',6,'-',0.5],
            # [(0.90*np.sqrt(Pk_nu_r))**2]   +['0.9*Pk_{fort,nu}','b',2,'-',0.5]
            ]
      plot_Pk2k1('Pk_{nu}','m_{nu}=%.2feV \,\,\, z = %.2f'%(mnu,z),Pk,show=show)

      fid = open(Path_nu+'/neutrinos/Pk/Pk_step_{:.4f}.txt'.format(z), 'rb')  # 以二进制模式打开文件
      Pk_CDM_r = np.fromfile(fid, dtype='float32').tolist()  # 读取数据到一维数组

      Pk = [
           [kh_r,Pk_CDM_r]   +['Pk_{fort,cdm}','r',2,'-',.5],
            [kh_nl,Pk_cdm_nonlin] +['Pk_{camb,cdm}' ,'.5',6,'-',0.5]
            ]
      plot_Pk2k1('Pk_{cdm}','m_{nu}=%.2feV \,\,\, z = %.2f'%(mnu,z),Pk,show=show)


      fid = open(Path_nu+'/neutrinos/tf/Tf_nu_{:.4f}.txt'.format(z), 'rb')  # 以二进制模式打开文件
      TF_iqr = np.fromfile(fid, dtype='float32').tolist()  # 读取数据到一维数组
      TF_camb = ((1-f_nu)*np.sqrt(Pk_cdm_nonlin)+f_nu*(np.sqrt(Pk_nu_nonlin)))/np.sqrt(Pk_cdm_nonlin)

      
      Pk = [[kh_r,TF_iqr]   +['Tf_{fort,dm}','r',2,'-',.5],
            [kh_nl,TF_camb] +['Tf_{camb,dm}' ,'.5',6,'-',0.5]]
            
      plot_Pk2k1('TF','m_{nu}=%.2feV \,\,\, z = %.2f'%(mnu,z),Pk,show=show)
      return z

# init
#set cosmology parm
H0=match_para('h0')*100
ombh2=match_para('omega_bar')*(H0/100)**2
omch2=match_para('omega_cdm')*(H0/100)**2
omk=0.0
neutrino_hierarchy='degenerate'
num_massive_neutrinos=3
mnu=match_para('Mass_nu')
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

print(npbin,kmin,kmax)
n=50
ni=50


# 使用glob.glob获取所有匹配 Tf_nu_*.txt 模式的文件名
file_list = glob.glob(Path_nu+'/neutrinos/tf/Tf_nu*.txt')

# 初始化一个列表来保存解析出的数字
numbers = []

# 使用正则表达式来匹配文件名中的数字部分
pattern = re.compile(Path_nu+r'/neutrinos/tf/Tf_nu_(.*?).txt')

# 遍历文件名列表
for file_name in file_list:
    # 使用正则表达式查找匹配
    match = pattern.search(file_name)
    if match:
        # 将匹配到的数字转换为浮点数并添加到列表中
        number = float(match.group(1))
        numbers.append(number)

z_checkpoint = sorted(numbers, reverse=True)

# 打印结果
print(z_checkpoint)


# 使用multiprocessing.Pool创建进程池
pool = mp.Pool(processes=mp.cpu_count())

# 使用pool.map将worker函数应用到z_checkpoint的每个元素上
results = pool.map(worker, z_checkpoint)

# 关闭进程池
pool.close()
pool.join()

# 打印结果
print(results)

