#smatrix 定义一个函数analyze_qsh3()，返回一个量子自旋霍尔系统的传输矩阵
def analyze_qsh3():
    params= dict(A=3.65,B=-68.6,D=-51.1, M=-0.01,C=0)
    syst = qsh_system ()

    smatrix=kwant.smatrix(syst, energy=0.0, params=[params])
    return matrix.conductance_matrix(smatrix)
analyze_qsh3()


#通过传输矩阵和电流分布求各端口电压，该段代码只能计算由一端口流入，另一端口流出时的电压分布
def find_voltages(tmatrix,current):                                 #tmatrix是传输矩阵，在此处还未定义，current是电流分布[1, 0, 0, -1, 0, 0] 表示电流从端口1流入，端口4流出。
    try:                                                            #np.linalg.solve()函数用于求解线性方程组，返回值是一个数组，数组中的每个元素是一个未知数的解    current = [1, 0, 0, -1, 0, 0] → current[1:] = [0, 0, -1, 0, 0]
        voltage = np.linalg.solve(tmatrix[1:, 1:], current[1:])     #tmatrix[1:, 1:]表示tmatrix的第一行和第一列被删去 current[1:]表示current的第一个元素被删去  
    except np.linalg.LinAlgError as err:                            #当检测到奇异矩阵时，即该函数无解，避免程序崩溃，强制返回零电压（假设所有端口短路）
        if 'Singular matrix' in str(err):
            print("Hi there, singular matrix here")
            voltage = [0,0,0,0,0]
        else:
            raise
    return [0, *voltage]    ##在电压列表前添加一个 0，表示第一个端口的电压为0（参考接地）。

##
def resistance_xx(tmatrix):
    current = [1, 0, 0, -1, 0, 0]
    voltage = find_voltages(tmatrix,current)
    
    def resistance(lead1, lead2): 
        return voltage[lead1-1]-voltage[lead2-1] ##默认电流为1 ###将端口的序号索引转化为数组索引

    return resistance(2, 3) ##电流从端口2流入，端口3流出，电压差为V2-V3
##
def resistance_xy(tmatrix):
    current = [1, 0, 0, -1, 0, 0]
    voltage = find_voltages(tmatrix,current)
    
    def resistance(lead1, lead2):
        return voltage[lead1-1]-voltage[lead2-1]

    return resistance(2, 6)
###将smatrix转化为tmatrix
tmatrix=analyze_qsh3()



(resistance_xx(tmatrix),
resistance_xy(tmatrix))




