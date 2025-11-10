import kwant
import numpy as np
import matplotlib.pyplot as plt
import tinyarray

# 定义二维正方晶格（norbs=4 对应 4x4 Gamma 矩阵）
lat = kwant.lattice.square(norbs=4)  # 二维晶格，非三维立方

# ================== 定义 Gamma 矩阵 ==================
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
  # Γ6 矩阵
# ================== 模型参数 ==================
B,m, a,g = -300,30,1,300


# ================== 定义在位能和跃迁项 ==================
onsite =  m/(2)*gamma_5 - (4*B)/(a**2) * gamma_6
hopping_x = B/(a**2)*gamma_6 + 1j*g/(2*a)*gamma_4  # x方向跃迁
hopping_y = B/(a**2)*gamma_6 - 1j*g/(2*a)*gamma_7  # y方向跃迁

# ================== 构建沿 x 和 y 方向的导带系统 ==================
def make_syst(app):#该函数用于构造指定方向，例如app（*,&）
    """构造导带系统，沿指定方向无限延伸"""
    ssd = kwant.Builder(kwant.TranslationalSymmetry(app))#构建了ssd系统，ssd具有(-1,0)方向的对称性
    ssd[lat(0, 0)] = onsite
    ssd[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopping_x
    ssd[kwant.builder.HoppingKind((0, 1), lat, lat)] = hopping_y
    return ssd.finalized()

# 沿 x 方向的导带
ssd_x = kwant.Builder(kwant.TranslationalSymmetry((-1, 0)))#构建了ssd系统，ssd具有(-1,0)方向的对称性
##限制y方向的范围为[0,30-1]nm ###问题1：lat.shape(shape, (0, 0))为什么有两个shape？该命名与函数名相同，需更改
def form(pos):
    x, y = pos
    return 0 <= y < 30 #注意限制之前为1，不为无穷，现在为30
ssd_x[lat.shape(form, (0, 0))] = onsite  #函数lat.shape
ssd_x[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopping_x
ssd_x[kwant.builder.HoppingKind((0, 1), lat, lat)] = hopping_y
ssd_x = ssd_x.finalized()

# 沿 y 方向的导带
ssd_y = make_syst(app=(0, -1))  #上一步我们定义了make_syst函数，这里我们直接调用了make_syst函数，app=(-1,0)表示沿x方向

# ================== 计算能带 ==================
def plot_bands(lead, title, k_label):
    momenta = np.linspace(-np.pi, np.pi, 100) ##np.访问numpy库中的linspace函数，生成-π到π的100个数 ,布里渊区标准定义 
    ban = kwant.physics.Bands(lead) ##函数kwant.physics.Bands(),提取导带系统在时空间的hamiltonian并转换到动量空间（通过Bloch`s theorem）
    energies = [ban(k) for k in momenta]

    plt.figure(figsize=(8, 5)) ##设置图形尺寸为 ​8英寸（宽）×5英寸（高）
    for band in range(len(energies[0])): ##？？？#####################################
        plt.plot(momenta, [e[band] for e in energies], lw=1.5)
    plt.xlabel(k_label, fontsize=12) ##字体大小12
    plt.ylabel("Energy", fontsize=12) ##不是变量的确定的名称用""包裹
    plt.title(title)
    plt.ylim(-150,150)
    plt.xlim(-0.5,0.5)
    plt.grid(ls="--", alpha=0.5) ##添加虚线网格，线的透明度为0.5
    plt.show() ##显示绘制好的窗口

# 绘制沿 x 方向的能带（Γ → X）
plot_bands(ssd_x, "Band Structure of 4-site Ribbon along Γ → X", r"$k_x$")

# 绘制沿 y 方向的能带（Γ → Y）
plot_bands(ssd_y, "Band Structure along Γ → Y (ky)", r"$k_y$")