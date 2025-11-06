import os,re,sys,subprocess,time
import concurrent.futures
import multiprocessing as mp
import os,re


# code para
origin_code_path = '/public/home/acd4q6s2ve/nu_test/cube_nu'
test_path = '/public/home/acd4q6s2ve/nu_test'

# simulation para

paras = [
     {'id': 'x_0'  , 'Mass_nu': 0.1 , 'z_start': 1, 'tf': 1, 'z_math': 1, 'nn':2 , 'box': 400 , 'PP': 0,'i':0}

    ,{'id': '1_x'  , 'Mass_nu': 0.0 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 400 , 'PP': 0,'i':0}
    # ,{'id': '1_x'  , 'Mass_nu': 0.0 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 400 , 'PP': 2,'i':3}
    # ,{'id': '1_x'  , 'Mass_nu': 0.0 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 400 , 'PP': 2,'i':2}
    # ,{'id': '1_x'  , 'Mass_nu': 0.0 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 400 , 'PP': 2,'i':1}
    # ,{'id': '1_1'  , 'Mass_nu': 0.01, 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 400 , 'PP': 2,'i':0}
    # ,{'id': '1_2'  , 'Mass_nu': 0.05, 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 400 , 'PP': 2,'i':0}
    # ,{'id': '1_3'  , 'Mass_nu': 0.2 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 400 , 'PP': 2,'i':0}
    # ,{'id': '1_4'  , 'Mass_nu': 0.5 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 400 , 'PP': 2,'i':0}

    # ,{'id': '2_1_x', 'Mass_nu': 0.0 , 'z_start': 2, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 400 , 'PP': 2,'i':0}
    # ,{'id': '2_1'  , 'Mass_nu': 0.1 , 'z_start': 2, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 400 , 'PP': 2,'i':0}
    # ,{'id': '2_2_x', 'Mass_nu': 0.0 , 'z_start': 3, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 400 , 'PP': 2,'i':0}
    # ,{'id': '2_2'  , 'Mass_nu': 0.1 , 'z_start': 3, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 400 , 'PP': 2,'i':0}

    # ,{'id': '3_1'  , 'Mass_nu': 0.1 , 'z_start': 1, 'tf': 0 , 'z_math': 1, 'nn':2 , 'box': 400 , 'PP': 2,'i':0}
    # ,{'id': '3_2'  , 'Mass_nu': 0.1 , 'z_start': 1, 'tf':-1 , 'z_math': 1, 'nn':2 , 'box': 400 , 'PP': 2,'i':0}

    # ,{'id': '4_1'  , 'Mass_nu': 0.1 , 'z_start': 1, 'tf': 1 , 'z_math': 0, 'nn':2 , 'box': 400 , 'PP': 2,'i':0}

    # ,{'id': '5_1_x', 'Mass_nu': 0.0 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':4 , 'box': 400 , 'PP': 2,'i':0}
    # ,{'id': '5_1'  , 'Mass_nu': 0.1 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':4 , 'box': 400 , 'PP': 2,'i':0}
    # ,{'id': '5_2_x', 'Mass_nu': 0.0 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':1 , 'box': 400 , 'PP': 2,'i':0}
    # ,{'id': '5_2'  , 'Mass_nu': 0.1 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':1 , 'box': 400 , 'PP': 2,'i':0}

    # ,{'id': '6_1_x', 'Mass_nu': 0.0 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 200 , 'PP': 2,'i':0}
    # ,{'id': '6_1'  , 'Mass_nu': 0.1 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 200 , 'PP': 2,'i':0}
    # ,{'id': '6_2_x', 'Mass_nu': 0.0 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 600 , 'PP': 2,'i':0}
    # ,{'id': '6_2'  , 'Mass_nu': 0.1 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 600 , 'PP': 2,'i':0}
    # ,{'id': '6_3_x', 'Mass_nu': 0.0 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 1200, 'PP': 2,'i':0}
    # ,{'id': '6_3'  , 'Mass_nu': 0.1 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 1200, 'PP': 2,'i':0}
    # ,{'id': '6_4_x', 'Mass_nu': 0.0 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 2400, 'PP': 2,'i':0}
    # ,{'id': '6_4'  , 'Mass_nu': 0.1 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 2400, 'PP': 2,'i':0}

    # ,{'id': '7_1'  , 'Mass_nu': 0.1 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 400 , 'PP': 1,'i':0}

    # ,{'id': '8_x'  , 'Mass_nu': 0.0 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 400 , 'PP': 0,'i':0}
    # ,{'id': '8_1'  , 'Mass_nu': 0.1 , 'z_start': 1, 'tf': 1 , 'z_math': 1, 'nn':2 , 'box': 400 , 'PP': 0,'i':0}
]

# nn = 2
# ngs = [1024]
# boxs = [1140]
# Mass_nus = [0,2]
# PPs = [0,1,2]
# z_starts = [1,2]

# paras = [
#         {'ng': ng, 'box': box, 'Mass_nu': Mass_nu, 'PP': PP, 'z_start': z_start,'i':i}
#         for ng in ngs
#         for box in boxs
#         for Mass_nu in Mass_nus
#         for PP in PPs
#         for z_start in z_starts
#         for i in range(times)
#         ]


def match_para(fn,para):
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

def check_id(id,sleeptime):
    check = True
    while (check):
        sq = subprocess.Popen('squeue', stdout=subprocess.PIPE, shell=True, text=True)
        output, error = sq.communicate()
        # print('output',output)
        # print('error',error)

        check = re.search(id, output)
        time.sleep(sleeptime)
    
def make_code(para):
    if (para['Mass_nu'] == 0):
        code_path = "{}/code/{id}/{i}".format(test_path, **para)
        opath =  "{}/output/{id}/{i}".format(test_path, **para)
        
    else:
        code_path = "{}/code/{id}_{Mass_nu:1.2f}/{i}".format(test_path, **para)
        opath = "{}/output/{id}_{Mass_nu:1.2f}/{i}".format(test_path, **para)

    os.system("mkdir -p {}".format(code_path))
    os.system("cp -r {}/* {}".format(origin_code_path,code_path))
    os.chdir(code_path)
    # os.system('pwd')

 #set parameters
    with open('./parameters.f90', 'r') as f:
        parameters_old = f.read()
    
    # set opath
    pattern = r"character\(\*\),parameter :: opath='(.*)'"
    replacement = "character(*),parameter :: opath='%s/'"%opath
    parameters_new = re.sub(pattern, replacement, parameters_old)

    # set Mass_nu
    pattern = r"real\(8\),parameter :: Mass_nu=(.*) ! Mass_nu/eV"
    replacement = "real(8),parameter :: Mass_nu={Mass_nu} ! Mass_nu/eV".format(**para)
    parameters_new = re.sub(pattern, replacement, parameters_new)

    # set tf_math
    pattern = r"integer\(8\),parameter :: calculate_PK = (.*) ! -1:read cb and nu 0:read  cb 1:calculate from pure CDM sim"
    replacement = "integer(8),parameter :: calculate_PK = {tf} ! -1:read cb and nu 0:read  cb 1:calculate from pure CDM sim".format(**para)
    parameters_new = re.sub(pattern, replacement, parameters_new)

    # set z_math
    pattern = r"real\(8\),parameter :: a_nu=(.*) ! nu is matter in a_nu"
    if (para['z_math'] == 1):
        replacement = "real(8),parameter :: a_nu=0 ! 1./(595./5.47*(Mass_nu/3)/0.1+0.01) ! nu is matter in a_nu"
    elif (para['z_math'] == 0):
        replacement = "real(8),parameter :: a_nu=1./(595./5.47*(Mass_nu/3)/0.1+0.01) ! nu is matter in a_nu"
    parameters_new = re.sub(pattern, replacement, parameters_new)

    # set ng
    pattern = r"integer,parameter :: nn=(.*) ! number of imgages (nodes) /dim"
    replacement = "integer,parameter :: nn={nn} ! number of imgages (nodes) /dim".format(**para)
    parameters_new = re.sub(pattern, replacement, parameters_new)

    # set box
    pattern = r"real,parameter :: box=(.*)"
    replacement = "real,parameter :: box={box}  ! simulation scale /dim, in unit of Mpc/h".format(**para)
    parameters_new = re.sub(pattern, replacement, parameters_new)
    
    with open('./parameters.f90', 'w') as f:
        f.write(parameters_new)

 #set kick
    with open('./kick.f90', 'r') as f:
        kick_old = f.read()

    #set PP
    if para['PP'] == 0 :
        pattern = r"logical,parameter :: PP=.(.*)."
        replacement ="logical,parameter :: PP=.false."
        kick_new = re.sub(pattern, replacement, kick_old)

        pattern = r"logical,parameter :: PP_corr=.(.*)."
        replacement ="logical,parameter :: PP_corr=.false."
        kick_new = re.sub(pattern, replacement, kick_new)

    if para['PP'] == 1 :
        pattern = r"logical,parameter :: PP=.(.*)."
        replacement ="logical,parameter :: PP=.true."
        kick_new = re.sub(pattern, replacement, kick_old)

        pattern = r"logical,parameter :: PP_corr=.(.*)."
        replacement ="logical,parameter :: PP_corr=.false."
        kick_new = re.sub(pattern, replacement, kick_new)

    if para['PP'] == 2 :
        pattern = r"logical,parameter :: PP=.(.*)."
        replacement ="logical,parameter :: PP=.true."
        kick_new = re.sub(pattern, replacement, kick_old)

        pattern = r"logical,parameter :: PP_corr=.(.*)."
        replacement ="logical,parameter :: PP_corr=.true."
        kick_new = re.sub(pattern, replacement, kick_new)

    with open('./kick.f90', 'w') as f:
        f.write(kick_new)

 #set ic
    with open('./utilities/ic.f90', 'r') as f:
        ic_old = f.read()
    
    # set start cur_checkpoint
    pattern = r"cur_checkpoint=(.*) !! set current checkpoint"
    replacement = "cur_checkpoint={z_start}!! set current checkpoint".format(**para)
    ic_new = re.sub(pattern, replacement, ic_old)

    # set READ_NOISE
    pattern = r"(.*)define READ_NOISE"
    if (para['Mass_nu'] == 0):
        replacement = "!#define READ_NOISE"
    else:
        replacement = "#define READ_NOISE"
    ic_new = re.sub(pattern, replacement, ic_new)

    with open('./utilities/ic.f90', 'w') as f:
        f.write(ic_new)
        
        
    os.system('make clean')
    os.system('make')
    os.chdir(code_path+'/utilities')
    os.system('make clean')
    os.system('make')


    
    return code_path

def sb_ic(para):

    code_path = make_code(para)

    if (para['i'] == 0 and para['Mass_nu'] == 0):
        os.chdir(code_path)
        os.system('pwd')
        print('*'*20+'\nget ic in path :%s'%code_path)
        ic_done = False
        sb_count = 0
        while(ic_done == False and sb_count < 10):
            sb_count = sb_count + 1
            get_id = False
            while(get_id == False):
                try:
                    sb_ic = subprocess.Popen('sbatch ic%d.sh'%(para['nn']**3), stdout=subprocess.PIPE, shell=True, text=True)
                    output, error = sb_ic.communicate()
                    id = re.search(r"Submitted batch job (\d+)", output).group(1)
                    get_id = True
                except:
                    None
            print('*'*20+id+'\nstart ic in path :%s'%code_path)

            check_id(id,2)
            
            try:
                ta_ic = subprocess.Popen('tail -n 2 ic%d.out'%(para['nn']**3), stdout=subprocess.PIPE, shell=True, text=True)
                ic_out, error = ta_ic.communicate()
                if (error != None):
                    print('*'*20+'ic tail_err in path :%s\n%s'%(code_path,error), file=sys.stderr)
                # print(code_path+id+'\n',ic_out)
                id_done = re.search("initial condition done", ic_out)
                if id_done:
                    print('*'*20+'ic done in path :%s'%code_path)
                    ic_done = True
                    break
            except:
                None
        if sb_count > 10:
            ta_ic = subprocess.Popen('tail -n 300 ic%d.err'%(para['nn']**3), stdout=subprocess.PIPE, shell=True, text=True)
            ic_err, error = ta_ic.communicate()
            print('*'*20+'ic run_err in path :%s\n%s'%(code_path,ic_err), file=sys.stderr)
        
def sb_cicpower(code_path,para):
    os.chdir(code_path)
    # os.system('pwd')
    print('*'*20+'\nget cic in path :%s'%code_path)
    cicpower_done = False
    sb_count = 0
    while(cicpower_done == False and sb_count < 10):
        sb_count = sb_count + 1
        get_id = False
        while(get_id == False):
            try:
                sb_ic = subprocess.Popen('sbatch cicpower%d.sh'%(para['nn']**3), stdout=subprocess.PIPE, shell=True, text=True)
                output, error = sb_ic.communicate()
                id = re.search(r"Submitted batch job (\d+)", output).group(1)
                get_id = True
            except:
                None
        print('*'*20+id+'\nstart cic in path :%s'%code_path)

        check_id(id,5)

        try:
            ta_cicpower = subprocess.Popen('tail -n 2 cicpower.out', stdout=subprocess.PIPE, shell=True, text=True)
            cicpower_out, error = ta_cicpower.communicate()
            # print(code_path+id+'\n',cicpower_out)
            id_done = re.search("cicpower done", cicpower_out)
            if id_done:
                print('*'*20+'cicpower done in path :%s'%code_path)
                cicpower_done = True
                break
        except:
            None
    if sb_count > 10:
        ta_cicpower = subprocess.Popen('tail -n 300 cicpower.err', stdout=subprocess.PIPE, shell=True, text=True)
        cicpower_err, error = ta_cicpower.communicate()
        print('*'*20+'cicpower run_err in path :%s\n%s'%(code_path,cicpower_err), file=sys.stderr)

def sb_main(para):
    print('para:',para,time.time())
    # time.sleep(1)
    # print('para:',para,time.time())
    # return 0
    if (para['Mass_nu'] == 0):
        code_path = "{}/code/{id}/{i}".format(test_path, **para)
        
    else:
        code_path = "{}/code/{id}_{Mass_nu:1.2f}/{i}".format(test_path, **para)
        try:
            id0=[para1['id'] for para1 in paras 
                if (para1['Mass_nu'] == 0.0 
                and para1['z_start'] == para['z_start'] 
                and para1['tf'] == para['tf'] 
                and para1['z_math'] == para['z_math'] 
                and para1['nn'] == para['nn'] 
                and para1['box'] == para['box'] 
                and para1['PP'] == para['PP'] 
                and para1['i'] == para['i'])
                ][0]
        except:
            id0='1-x'
        # print(id0,para['nn']**3)
        for j in range(para['nn']**3):
            # print("cp -vr {}/output/{}/{i}/image{j}/noise*".format(test_path, id0,**para, j=j+1)
            #             +"  {}/output/{id}_{Mass_nu:1.2f}/{i}/image{j}".format(test_path, **para, j=j+1))
            os.system("mkdir -p {}/output/{id}_{Mass_nu:1.2f}/{i}/image{j}".format(test_path, **para, j=j+1))
            os.system("cp -vr {}/output/{}/{i}/image{j}/noise*".format(test_path, id0,**para, j=j+1)
                        +"  {}/output/{id}_{Mass_nu:1.2f}/{i}/image{j}".format(test_path, **para, j=j+1))
        # print('cp done')
        os.chdir(code_path)
        os.system('pwd')
        print('\n\n'+'!'*20+'\nget ic in path :%s'%code_path)
        ic_done = False
        sb_count = 0
        while(ic_done == False and sb_count < 10):
            sb_count = sb_count + 1
            get_id = False
            while(get_id == False):
                try:
                    sb_ic = subprocess.Popen('sbatch ic%d.sh'%(para['nn']**3), stdout=subprocess.PIPE, shell=True, text=True)
                    output, error = sb_ic.communicate()
                    id = re.search(r"Submitted batch job (\d+)", output).group(1)
                    get_id = True
                except:
                    None
            print('*'*20+id+'\nstart ic in path :%s'%code_path)

            check_id(id,2)

            try:
                ta_ic = subprocess.Popen('tail -n 2 ic%d.out'%(para['nn']**3), stdout=subprocess.PIPE, shell=True, text=True)
                ic_out, error = ta_ic.communicate()
                if (error != None):
                    print('*'*20+'ic tail_err in path :%s\n%s'%(code_path,error), file=sys.stderr)
                id_done = re.search("initial condition done", ic_out)
                if id_done:
                    print('*'*20+'ic done in path :%s'%code_path)
                    ic_done = True
                    break
            except:
                None


    os.chdir(code_path)
    id='00000'
    os.system('pwd')
    main_done = False
    sb_count = 0
    get_id = False
    while(get_id == False):
        try:
            sb_ic = subprocess.Popen('sbatch main%d.sh'%(para['nn']**3), stdout=subprocess.PIPE, shell=True, text=True)
            output, error = sb_ic.communicate()
            id = re.search(r"Submitted batch job (\d+)", output).group(1)
            get_id = True
        except:
            None
    print('\n'+'='*20+id+'\nstart main in path :%s'%code_path)
    while(main_done == False and sb_count < 5):

        sb_count = sb_count + 1

        check_id(id,5)
        checkpoint = -1
        try :
            with open('main%d.out'%(para['nn']**3), 'r') as file:
                text = file.read()

            matches = re.finditer(r'cur_checkpoint\s*=\s*(\d+)', text, re.MULTILINE)
            for match in matches:
                if match:
                    # 假设我们想要找到最后一个匹配的整数
                    checkpoint = int(match.group(1))
            
            if checkpoint > 103:
                # print('*'*20+id+'\n main out in path :%s   \n   z_checkpoint = %d'%(code_path,checkpoint))
                ta_main = subprocess.Popen('tail -n 30 main%d.out'%(para['nn']**3), stdout=subprocess.PIPE, shell=True, text=True)
                main_out, error = ta_main.communicate()
                # print('*'*20+'main done in path :%s \n main.out: \n %s'%(code_path,main_out))
                main_done = True


                sb_cicpower(code_path,para)



                break
        except:
            None
                
        if (checkpoint < 1):
            checkpoint = 2
        # print('*'*20+id+'\n main out in path :%s   \n   z_checkpoint = %d'%(code_path,checkpoint))
        checkpoint = checkpoint-1

        print('*'*20+id+code_path,main_done,sb_count,checkpoint)
        if (main_done == False):
            ta_main = subprocess.Popen('tail -n 30 main%d.out'%(para['nn']**3), stdout=subprocess.PIPE, shell=True, text=True)
            main_out, error = ta_main.communicate()
            print('\n\n'*5+'*'*20+id+'\n                    Fail  main in path :%s \n               restart in z_checkpoint = %d %d'%(code_path,checkpoint,sb_count)+main_out+'\n'*5)

            os.system('cp main%d.out log/main%d_%d.log'%(para['nn']**3,para['nn']**3,checkpoint))

            with open('./initialize.f90', 'r') as f:
                init_code = f.read()

            pattern_init = r"\s*sim%cur_checkpoint=(\d+)\s*!\s*change\s*for\s*resuming\s*checkpoints"
            replacement_init =f"\n  sim%cur_checkpoint={checkpoint} ! change for resuming checkpoints"
            modified_init_code = re.sub(pattern_init, replacement_init, init_code)

            pattern_init = r"logical,parameter :: read_Gks=.false."
            replacement_init =f"logical,parameter :: read_Gks=.true."
            modified_init_code = re.sub(pattern_init, replacement_init, modified_init_code)

            with open('./initialize.f90', 'w') as f:
                f.write(modified_init_code)

            os.system('make')
            get_id = False
            while(get_id == False):
                try:
                    sb_ic = subprocess.Popen('sbatch main%d.sh'%(para['nn']**3), stdout=subprocess.PIPE, shell=True, text=True)
                    output, error = sb_ic.communicate()
                    id = re.search(r"Submitted batch job (\d+)", output).group(1)
                    get_id = True
                except:
                    None


    if (sb_count > 10):
        ta_main = subprocess.Popen('tail -n 300 main%d.err'%(para['nn']**3), stdout=subprocess.PIPE, shell=True, text=True)
        main_err, error = ta_main.communicate()
        print('*'*20+'main run_err in path :%s\n%s'%(code_path,main_err), file=sys.stderr)




ncore = min([mp.cpu_count(),len(paras)])
print(ncore)

##### get ic   #####
# 使用multiprocessing.Pool创建进程池
pool = mp.Pool(processes=ncore)

# 使用pool.map将worker函数应用到z_checkpoint的每个元素上
pool.map(sb_ic, paras)

# 关闭进程池
pool.close()
pool.join()


print('\n\n\n\n\n\n\n')
print('\n\n\n\n\n\n\n')
print('!'*10+'run main'+'!'*10)
print('\n\n\n\n\n\n\n')
print('\n\n\n\n\n\n\n')

##### run main   #####
# 使用multiprocessing.Pool创建进程池
pool = mp.Pool(processes=ncore)

# 使用pool.map将worker函数应用到z_checkpoint的每个元素上
pool.map(sb_main, paras)

# 关闭进程池
pool.close()
pool.join()

print('\n\n\n\n\n\n\n')
print('+'*30+'\n\n\n        mutile run done       \n\n\n'+'+'*30)

