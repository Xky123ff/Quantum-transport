import kwant
import numpy as np
import matplotlib.pyplot as plt
import tinyarray

lat = kwant.lattice.cubic(norbs=4)  

# 构建导带（沿x方向无限周期性）
syst = kwant.Builder(kwant.TranslationalSymmetry((-1, 0, 0)))  # 晶格常数设为1，平移对称方向为x轴

# 定义伽马矩阵
sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])
#定义4*4厄米矩阵
I_4 = tinyarray.array(np.kron(sigma_0, sigma_0))
gamma_1 = tinyarray.array(np.kron(sigma_x, sigma_x))
gamma_2 = tinyarray.array(np.kron(sigma_y, sigma_x))
gamma_3 = tinyarray.array(np.kron(sigma_z, sigma_x))
gamma_5 = tinyarray.array(np.kron(sigma_0, sigma_z))
# 参数
C_0, C_1, C_2 = 1, 2, 3
v, v_z, las = 5, 3, 1
M_0, M_1, M_2 = 0, 1, 2

# 导带的在位能和跃迁项
syst[lat(0, 0, 0)] = (C_0 + (2*C_1 + 4*C_2)/las**2) * I_4 + (M_0 + (2*M_1 + 4*M_2)/las**2) * gamma_5
syst[lat.neighbors(1)] = -C_2/(las**2)*I_4 - 1j*v/(2*las)*gamma_1 - M_2/(las**2)*gamma_5  # x方向跃迁

# 其他方向（y/z）的跃迁项（若需要三维周期性，需定义更多导带，此处假设仅沿x方向无限）

# 将导带系统 finalize
syst = syst.finalized()

# 计算沿 Γ → X 方向的能带
momenta = np.linspace(0, np.pi, 100)  # kx 从 0 到 π
bands = kwant.physics.Bands(syst, params=None)  # 无需额外参数
energies = [bands(k) for k in momenta]

# 绘图
plt.figure(figsize=(8, 5))
for band in range(4):  # 4个能带（对应4×4哈密顿量）
    plt.plot(momenta, [energy[band] for energy in energies])
plt.xlabel(r"$k_x$ (Γ → X)", fontsize=12)
plt.ylabel("Energy", fontsize=12)
plt.title("3D Topological Insulator Band Structure ")
plt.grid(ls="--", alpha=0.5)
plt.show()