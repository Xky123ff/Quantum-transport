import kwant
import tinyarray
import numpy as np
import matplotlib.pyplot as plt

sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])

I_4 = tinyarray.array(np.kron(sigma_0, sigma_0))       # 4x4 单位矩阵
gamma_1 = tinyarray.array(np.kron(sigma_x, sigma_x))   # Γ1 矩阵
gamma_2 = tinyarray.array(np.kron(sigma_y, sigma_x))   # Γ2 矩阵
gamma_3 = tinyarray.array(np.kron(sigma_z, sigma_x))   # Γ3 矩阵
gamma_4 = tinyarray.array(np.kron(sigma_0, sigma_y))   # Γ4 矩阵
gamma_5 = tinyarray.array(np.kron(sigma_0, sigma_z))   # Γ5 矩阵
gamma_6 = tinyarray.array(np.kron(sigma_z, sigma_z)) 
gamma_7 = tinyarray.array(np.kron(sigma_0, sigma_x))

#定义参数
B,m,a,g,las = -300,20,1,300,1

def ma_sy(width=60, length=300):

    structure=kwant.lattice.square(las)


    onsite_1 = - (4*B)/(a**2) * gamma_6
    onsite_2 = m/(2)*gamma_5
    hop_x = B/(a**2)*gamma_6 + 1j*g/(2*a)*gamma_4
    hop_y = B/(a**2)*gamma_6 - 1j*g/(2*a)*gamma_7

    def onsit(site): #不能略去（site）
        return onsite_1 + onsite_2

    def form(pos):
        x, y = pos
        return 0 < x <= length and 0 < y <= width

    sy= kwant.Builder()
    sy[structure.shape(form,(1,1))] = onsit
    sy[kwant.builder.HoppingKind((1,0), structure)]  = hop_x
    sy[kwant.builder.HoppingKind((0,1), structure)]  = hop_y

    def electronGas_model():
    
        onsite_leads= onsit
        hopping_leads_x= hop_x
        hopping_leads_y= hop_y

        ########### create a translational invariant system
        model_sym= kwant.TranslationalSymmetry(structure.vec((-1,0)),structure.vec((0,-1))) #model_sym是具有两个方向的对称性
        model_1= kwant.Builder(model_sym)

        model_1[structure(0, 0)] = onsite_leads
        model_1[kwant.builder.HoppingKind((1, 0), structure, structure)]  = hopping_leads_x
        model_1[kwant.builder.HoppingKind((0, 1), structure, structure)]  = hopping_leads_y
        return model_1

    width_lead = length//5 

    lead0= kwant.Builder(kwant.TranslationalSymmetry((-1,0)))#lead0是具有(-1,0)方向的对称性
    lead5= kwant.Builder(kwant.TranslationalSymmetry((0,-1)))
    lead4= kwant.Builder(kwant.TranslationalSymmetry((0,-1)))

    def shape_lead0(site):#定义导线的形状
        x, y = site.pos
        return 0 < y <= width

    def shape_lead5(site):
        x, y = site.pos
        return width_lead   < x <= 2*width_lead

    def shape_lead4(site):
        x, y = site.pos
        return 3*width_lead < x <= 4*width_lead

    lead0.fill(electronGas_model(), shape_lead0, (0, 1)) #第一个参数：电极的紧束缚模型（已预先定义）；​第二个参数：电极在空间中的形状定义函数​；第三个参数：电极在系统中的起始坐标（原胞位置）
    lead5.fill(electronGas_model(), shape_lead5, (width_lead+1, 0))
    lead4.fill(electronGas_model(), shape_lead4, (3*width_lead+1, 0))
    lead1=lead5.reversed() #lead1是lead5的逆向
    lead2=lead4.reversed() #lead2是lead4的逆向
    lead3=lead0.reversed() #lead3是lead0的逆向

    for lead in [lead0,lead1,lead2,lead3,lead4,lead5]:   #lead0.reversed()表示lead0的逆向
        sy.attach_lead(lead)  #attach_lead()函数用于将导线连接到系统中
    return sy

kwant.plot(ma_sy(),site_color='blue', hop_color='red', lead_color='green')

sys=ma_sy()
sys=sys.finalized()
########################################计算电阻##############################################################################


#smatrix 定义一个函数smatx()，返回一个量子自旋霍尔系统的电导传输矩阵
def smatx():
    syst = sys
    smatrix=kwant.smatrix(syst, energy=3.0)
    
    print("Is Hermitian:", np.allclose(smatrix.data, smatrix.data.conj().T))
    
    return  smatrix.conductance_matrix()  # 返回电导矩阵
smatx()
tmatrix=smatx()


#通过传输矩阵和电流分布求各端口电压，该段代码只能计算由一端口流入，另一端口流出时的电压分布
def find_voltages(tmatrix,current):  #注意：该函数返回的V是所有V的矩阵                         #tmatrix是传输矩阵，在此处还未定义，current是电流分布[1, 0, 0, -1, 0, 0] 表示电流从端口1流入，端口4流出。
    try:                                                                                    #np.linalg.solve()函数用于求解线性方程组，返回值是一个数组，数组中的每个元素是一个未知数的解   
        voltage = np.linalg.solve(tmatrix[1:, 1:], current[1:])                             #tmatrix[1:, 1:]表示tmatrix的第一行和第一列被删去 current[1:]表示current的第一个元素被删去  
    except np.linalg.LinAlgError as err:                                                    #当检测到奇异矩阵时，即该函数无解，避免程序崩溃，强制返回零电压
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
        return voltage[1]-voltage[2] ##默认电流为1 ###将端口的序号索引转化为数组索引

    return resistance(1, 2) ##电压差为V2-V3
##
def resistance_xy(tmatrix):
    current = [1, 0, 0, -1, 0, 0]
    voltage = find_voltages(tmatrix,current)
    
    def resistance(lead1, lead2):
        return voltage[1]-voltage[2]

    return resistance(1, 5)
###将smatrix转化为tmatrix


(resistance_xx(tmatrix),
resistance_xy(tmatrix))

print(resistance_xx(tmatrix))
print(resistance_xy(tmatrix))