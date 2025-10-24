import kwant
import tinyarray
import numpy as np
import random
import matplotlib.pyplot as plt

sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])
I_4     = tinyarray.array(np.kron(sigma_0, sigma_0))
gamma_1 = tinyarray.array(np.kron(sigma_x, sigma_x))
gamma_2 = tinyarray.array(np.kron(sigma_y, sigma_x))
gamma_3 = tinyarray.array(np.kron(sigma_z, sigma_x))
gamma_5 = tinyarray.array(np.kron(sigma_0, sigma_z))
kron_z0 = tinyarray.array(np.kron(sigma_z, sigma_0))
kron_x0 = tinyarray.array(np.kron(sigma_x, sigma_0))

C_0 = -0.0068
C_1 = 0.  #隐式浮点数
C_2 = -0.196
M_0  = 0.28
M_1  = 0.
M_2  = -0.566
v = 0.41
m = 0.4

lat_spacing = 1.
las = lat_spacing

def ma_sy(width=20, length=30, imp_percent=5, imp_dis='default'):

    structure=kwant.lattice.square(las)

    ########### onsite and hopping energies 
    onsite_1 = (C_0 + (2*C_1 + 4*C_2)/(las**2))*I_4 + (M_0 + (2*M_1 + 4*M_2)/(las**2))*gamma_5
    onsite_2 = m*kron_z0
    hop_x = -C_2/(las**2)*I_4 - 1j*v/(2*las)*gamma_1 - M_2/(las**2)*gamma_5
    hop_y = -C_2/(las**2)*I_4 - 1j*v/(2*las)*gamma_2 - M_2/(las**2)*gamma_5

    def onsit(): #略去（site）
        return onsite_1 + onsite_2

    def form(pos):
        x, y = pos
        return 0 < x <= length and 0 < y <= width

    sy= kwant.Builder()
    sy[structure.shape(form,(1,1))] = onsit
    sy[kwant.builder.HoppingKind((1,0), structure)]  = hop_x
    sy[kwant.builder.HoppingKind((0,1), structure)]  = hop_y

    def electronGas_model():
        ########### crystal structure
        structure=kwant.lattice.square(las)

        onsite_leads= 0.3 *I_4
        hopping_leads= -0.3 *I_4

        ########### create a translational invariant system
        model_sym= kwant.TranslationalSymmetry(structure.vec((-1,0)),structure.vec((0,-1)))
        model_1= kwant.Builder(model_sym)

        model_1[structure(0, 0)] = onsite_leads
        model_1[structure.neighbors()]  = hopping_leads

        return model_1

    width_lead = length//5 

    lead0= kwant.Builder(kwant.TranslationalSymmetry((-1,0)))
    lead1= kwant.Builder(kwant.TranslationalSymmetry((0,1)))
    lead2= kwant.Builder(kwant.TranslationalSymmetry((0,1)))

    def shape_lead0(site):
        x, y = site.pos
        return 0 < y <= width

    def shape_lead1(site):
        x, y = site.pos
        return width_lead < x <= 2*width_lead

    def shape_lead2(site):
        x, y = site.pos
        return 3*width_lead < x <= 4*width_lead

    lead0.fill(electronGas_model(), shape_lead0, (0, 1)) #第一个参数：电极的紧束缚模型（已预先定义）；​第二个参数：电极在空间中的形状定义函数​；第三个参数：电极在系统中的起始坐标（原胞位置）
    lead1.fill(electronGas_model(), shape_lead1, (width_lead+1, 0))
    lead2.fill(electronGas_model(), shape_lead2, (3*width_lead+1, 0))

    for lead in [lead0,lead1,lead2, lead0.reversed(),lead2.reversed(),lead1.reversed()]:   #lead0.reversed()表示lead0的逆向
        sy.attach_lead(lead)  #attach_lead()函数用于将导线连接到系统中
    return sy
kwant.plot(ma_sy(),site_color='blue', hop_color='red', lead_color='green')




