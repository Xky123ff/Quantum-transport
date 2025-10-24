import kwant
import numpy as np
import matplotlib.pyplot as plt  
lat = kwant.lattice.square()  
syst = kwant.Builder()  
#定义onsite能量
def onsite(site):
  return -4.0
#定义hopping能量
def hopping(site1,site2,B=0.1):
  return -1.0

#添加格子点
for x in range(10):
    for y in range(10):
        syst[lat(x, y)] = onsite(lat(x,y)) 

#添加连接
for x in range(10):
    for y in range(10):

          if x<9: # 水平方向（右移）
           syst[lat(x, y), lat(x + 1, y)] = hopping(lat(x, y), lat(x + 1, y))
          if y<9: 
           syst[lat(x, y), lat(x, y + 1)] = hopping(lat(x, y), lat(x, y + 1))

#添加左侧引线
sym_left = kwant.TranslationalSymmetry((-1,0))
lead_left = kwant.Builder(sym_left)
for y in range(10):
   lead_left[lat(0,y)] = onsite(lat(0,y))
# 1.添加 y 方向的 hopping
lead_left[kwant.builder.HoppingKind((0,1),lat,lat)] = hopping
# 2.添加 x 方向的 hopping（关键修复）
lead_left[kwant.builder.HoppingKind((1,0),lat,lat)] = hopping
syst.attach_lead(lead_left)

#添加右侧引线
sym_right = kwant.TranslationalSymmetry((1,0))
lead_right = kwant.Builder(sym_right)
for y in range(10):
   lead_right[lat(0,y)] = onsite(lat(9,y))
# 1.添加 y 方向的 hopping
lead_right[kwant.builder.HoppingKind((0,1),lat,lat)] = hopping
# 2.添加 x 方向的 hopping（关键修复）
lead_right[kwant.builder.HoppingKind((1,0),lat,lat)] = hopping
syst.attach_lead(lead_right)

#完成系统
syst = syst.finalized()

#可视化系统
kwant.plot(syst, site_color='blue', hop_color='red',lead_color='green')  # 绘制系统，点用蓝色表示，连接线用红色表示
plt.title("带 onsite, hopping 和终端的 10x10 方格 lattice")  # 给图形加标题
plt.show()  # 显示图形

#计算传输（transmission）在不同能量下的值
energies=np.linspace(-2.0,2.0,41) # 能量范围：-2.0 到 2.0 meV，41 个点
transmission=[]
for energy in energies:
   smatrix = kwant.smatrix(syst,energy) #计算散射矩阵
   trans = smatrix.transmission(0,1)   # 从左侧终端 (0) 到右侧终端 (1) 的传输系数
   transmission.append(trans)

# 10. 绘制传输随能量的变化
plt.figure()
plt.plot(energies,transmission,'b-')
plt.xlabel("能量(mev)")
plt.ylabel("传输系数")
plt.title("传输系数随能量的变化")
plt.grid(True)
plt.show()

print("代码运行成功！")