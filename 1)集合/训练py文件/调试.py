import kwant
import matplotlib.pyplot as plt
import numpy as np

# 1. 创建一个 2D 方格晶格
lat = kwant.lattice.square(norbs=1)

# 2. 构建一个 10x10 的系统
syst = kwant.Builder()

# 3. 定义 onsite 能量
def onsite(site):
    return -4.0

# 4. 定义 hopping 能量
def hopping(site1, site2):
    return -1.0

# 5. 添加格子点和连接
for x in range(10):
    for y in range(10):
        syst[lat(x, y)] = onsite(lat(x, y))
for x in range(10):
    for y in range(10):        
        if x < 9:
            syst[lat(x, y), lat(x + 1, y)] = hopping(lat(x, y), lat(x + 1, y))
        if y < 9:
            syst[lat(x, y), lat(x, y + 1)] = hopping(lat(x, y), lat(x, y + 1))

# 6. 添加左侧引线
sym_left = kwant.TranslationalSymmetry((-1, 0))
lead_left = kwant.Builder(sym_left)
for y in range(10):
    lead_left[lat(0, y)] = onsite(lat(0, y))
lead_left[kwant.builder.HoppingKind((0, 1), lat, lat)] = hopping
lead_left[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopping
syst.attach_lead(lead_left)

# 7. 添加右侧引线
sym_right = kwant.TranslationalSymmetry((1, 0))
lead_right = kwant.Builder(sym_right)
for y in range(10):
    lead_right[lat(0, y)] = onsite(lat(0, y))  # 注意这里坐标应为 (0, y)，引线内部定义
lead_right[kwant.builder.HoppingKind((0, 1), lat, lat)] = hopping
lead_right[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopping
syst.attach_lead(lead_right)

# 8. 完成系统
syst = syst.finalized()

# 9. 可视化系统
kwant.plot(syst, site_color='blue', hop_color='red', lead_color='green')
plt.title("带 onsite 和引线的 10x10 方格晶格")
plt.show()

# 10. 计算传输系数
energies = np.linspace(-2.0, 2.0, 41)
transmission = []
for energy in energies:
    smatrix = kwant.smatrix(syst, energy)
    trans = smatrix.transmission(0, 1)
    transmission.append(trans)

# 11. 绘制传输系数随能量的变化
plt.figure()
plt.plot(energies, transmission, 'b-')
plt.xlabel("能量 (meV)")
plt.ylabel("传输系数")
plt.title("传输系数随能量的变化")
plt.grid(True)
plt.show()

print("代码1运行成功！")