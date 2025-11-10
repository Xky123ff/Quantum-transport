import kwant  
import numpy as np
import matplotlib.pyplot as plt  

# 定义方格lattice
lat = kwant.lattice.square()  

# 定义六端霍尔棒系统
def make_hall_bar(lead_width=2, W=10, L=20):
    syst = kwant.Builder()  
    
    def hall_bar_shape(pos):
        x, y = pos
        return (-L//2 <= x <= L//2) and (-W//2 <= y <= W//2)  # 修复了变量 X 未定义的问题

    # 添加格子点和连接
    syst[lat.shape(hall_bar_shape, (0, 0))] = -4.0
    syst[lat.neighbors()] = -1.0

    # 添加六个引线
    # 左侧引线
    sym_left = kwant.TranslationalSymmetry((-1, 0))
    lead_left = kwant.Builder(sym_left)
    for y in range(-lead_width//2, lead_width//2):    
        lead_left[lat(0, y)] = -4.0
    lead_left[lat.neighbors()] = -1.0
    syst.attach_lead(lead_left, origin=lat(-L//2, 0))

    # 右侧引线
    sym_right = kwant.TranslationalSymmetry((1, 0))
    lead_right = kwant.Builder(sym_right)
    for y in range(-lead_width//2, lead_width//2):    
        lead_right[lat(0, y)] = -4.0
    lead_right[lat.neighbors()] = -1.0
    syst.attach_lead(lead_right, origin=lat(L//2, 0))

    # 上侧引线(两对）
    sym_up = kwant.TranslationalSymmetry((0, 1))
    lead_up = kwant.Builder(sym_up)
    for x in range(-lead_width//2, lead_width//2):    
        lead_up[lat(x, 0)] = -4.0
    lead_up[lat.neighbors()] = -1.0
    syst.attach_lead(lead_up, origin=lat(-L//4, W//2))
    syst.attach_lead(lead_up, origin=lat(L//4, W//2))

    # 下侧引线(两对）
    sym_down = kwant.TranslationalSymmetry((0, -1))
    lead_down = kwant.Builder(sym_down)
    for x in range(-lead_width//2, lead_width//2):    
        lead_down[lat(x, 0)] = -4.0
    lead_down[lat.neighbors()] = -1.0
    syst.attach_lead(lead_down, origin=lat(-L//4, -W//2))
    syst.attach_lead(lead_down, origin=lat(L//4, -W//2))

    return syst.finalized()  # 修复了 return 语句的位置问题

# 构建系统
syst = make_hall_bar()

# 可视化系统
kwant.plot(syst, site_color='blue', hop_color='red', lead_color='green')  # 绘制系统，点用蓝色表示，连接线用红色表示
plt.title("六端霍尔棒系统")  # 给图形加标题
plt.show()  # 显示图形

# 计算传输（transmission）此处仅计算左右输出
energies = np.linspace(-2.0, 2.0, 41)  # 能量范围
transmission = []  # 传输系数
for energy in energies:
    smatrix = kwant.smatrix(syst, energy)  # 计算S矩阵
    transmission.append(smatrix.transmission(0, 1))  # 从左侧到右侧的传输

# 绘制传输系数
plt.figure()
plt.plot(energies, transmission, 'b-')
plt.xlabel("能量(mev)")
plt.ylabel("传输系数")
plt.title("六端霍尔棒系统传输系数")
plt.grid(True)
plt.show()

print("代码运行成功！")