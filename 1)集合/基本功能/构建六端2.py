import kwant
import kwant.continuum
import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt  # 修正：原图中的错误库名

# 定义量子自旋霍尔系统的哈密顿量
def qsh_element(a=20, L=2000, W=1000):
    """构建二维量子自旋霍尔系统的紧束缚哈密顿量"""
    # 原图关键参数修正：
    # 使用正确的张量积符号 kron，修正矩阵符号 sigma_x/y/z
    hamiltonian = (
        "c * identity(4) + M * kron(sigma_0, sigma_z)"      # 带隙项
        "- B * (k_x**2 + k_y**2) * kron(sigma_0, sigma_z)"   # 质量项
        "- D * (k_x**2 + k_y**2) * kron(sigma_0, sigma_0)"   # 动能项（原图疑似笔误）
        "+ A * k_x * kron(sigma_z, sigma_x)"                 # 自旋轨道耦合项
        "- A * k_y * kron(sigma_0, sigma_y)"                 # 自旋轨道耦合项
    )
    return hamiltonian

# 离散化连续模型
template = kwant.continuum.discretize(
    qsh_element(), 
    grid=a  # 晶格常数
)

# 定义系统形状（矩形纳米带）
def shape(site):
    x, y = site.pos
    return (0 <= y < W) and (0 <= x < L)

# 定义三个铅的形状
def lead0_shape(site):
    x, y = site.pos
    return (0 <= y < W)  # 左侧铅

def lead1_shape(site):
    x, y = site.pos
    return (L/5 <= x < 2*L/5)  # 中间铅

def lead2_shape(site):
    x, y = site.pos
    return (3*L/5 <= x < 4*L/5)  # 右侧铅

# 初始化系统
syst = kwant.Builder()

# 填充主系统
syst.fill(template, shape, (0, 0))

# 创建并附加铅（修正原图的拼写错误：Boulder→Builder）
lead0 = kwant.Builder(kwant.TranslationalSymmetry([-a, 0]))
lead0.fill(template, lead0_shape, (0, 0))

lead1 = kwant.Builder(kwant.TranslationalSymmetry([0, -a]))
lead1.fill(template, lead1_shape, (L/5, 0))

lead2 = kwant.Builder(kwant.TranslationalSymmetry([0, -a]))
lead2.fill(template, lead2_shape, (3*L/5, 0))

# 附加铅到系统（修正原图的拼写错误：fr→fill，attach_load→attach_leads）
syst.attach_leads(lead0)
syst.attach_leads(lead1)
syst.attach_leads(lead2)
syst.attach_leads(lead0.reversed())
syst.attach_leads(lead2.reversed())
syst.attach_leads(lead1.reversed())

# 最终化系统并求解
syst.finalize()

# 绘制系统结构（修正原图的拼写错误：split→show）
kwant.plot(syst)
plt.show()