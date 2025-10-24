import kwant  
import matplotlib.pyplot as plt  
lat = kwant.lattice.square()  
syst = kwant.Builder()  
def onsite(site):
  return -4.0
def hopping(site1,site2):
  return -1.0

for x in range(10):
    for y in range(10):
        syst[lat(x, y)] = onsite(lat(x,y)) 
for x in range(10):
    for y in range(10):

          if x<9: # 水平方向（右移）
           syst[lat(x, y), lat(x + 1, y)] = hopping(lat(x, y), lat(x + 1, y))
          if y<9: 
           syst[lat(x, y), lat(x, y + 1)] = hopping(lat(x, y), lat(x, y + 1))
# 3. 完成系统（让 kwant 准备好进行计算）
syst = syst.finalized()

# 4. 可视化系统
kwant.plot(syst, site_color='blue', hop_color='red')  # 绘制系统，点用蓝色表示，连接线用红色表示
plt.title("带 onsite 和 hopping 能量的 10x10 方格 lattice")  # 给图形加标题
plt.show()  # 显示图形
print("代码运行成功！")