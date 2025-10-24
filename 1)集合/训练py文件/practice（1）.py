import kwant  # 导入 kwant 库，用于构建量子系统
import matplotlib.pyplot as plt  # 导入 matplotlib，用于绘制图形

lat = kwant.lattice.square() 

syst = kwant.Builder()  # 创建一个空的系统

for x in range(1,11):
    for y in range(10):
        syst[lat(x, y)] = 0  # 每个点能量设为 0

syst = syst.finalized()

kwant.plot(syst, site_color='blue', hop_color='red')  # 绘制系统，点用蓝色表示，连接线用红色表示
plt.title("一个 10*10 方格 lattice")  # 给图形加标题
plt.show()  # 显示图形

print("代码运行成功！这是一个基础的 2D 方格，接下来我们可以添加能量和连接。")