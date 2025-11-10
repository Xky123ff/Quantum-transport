import kwant  
import matplotlib.pyplot as plt  
import tinyarray
import numpy as np

lat = kwant.lattice.cubic(norbs=4)
syst = kwant.Builder()  
sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])

# 定义伽马矩阵（使用 tinyarray）
I_4 = tinyarray.array(np.kron(sigma_0, sigma_0))
gamma_1 = tinyarray.array(np.kron(sigma_x, sigma_x))
gamma_2 = tinyarray.array(np.kron(sigma_y, sigma_x))
gamma_3 = tinyarray.array(np.kron(sigma_z, sigma_x))
gamma_5 = tinyarray.array(np.kron(sigma_0, sigma_z))

# 参数设置
C_0, C_1, C_2 = 1, 2, 3
v, v_z, las = 5, 3, 1
M_0, M_1, M_2 = 0, 1, 2

# 构造哈密顿量项
onsite = (C_0 + (2*C_1 + 4*C_2)/las**2) * I_4 + (M_0 + (2*M_1 + 4*M_2)/las**2) * gamma_5
hop_x = -C_2/(las**2)*I_4 - 1j*v/(2*las)*gamma_1 - M_2/(las**2)*gamma_5
hop_y = -C_2/(las**2)*I_4 - 1j*v/(2*las)*gamma_2 - M_2/(las**2)*gamma_5
hop_z = -C_1/(las**2)*I_4 - 1j*v_z/(2*las)*gamma_3 - M_1/(las**2)*gamma_5

# 在动量空间​中直接构造三维紧束缚模型的哈密顿量矩阵
def hamiltonian_k(kx, ky, kz):
    H = tinyarray.array(onsite)  # 复制，确保后续操作不会修改原始 onsite 矩阵
    H += hop_x * np.exp(1j*kx) + np.conjugate(np.transpose(hop_x)) * np.exp(-1j*kx)
    H += hop_y * np.exp(1j*ky) + np.conjugate(np.transpose(hop_y)) * np.exp(-1j*ky)
    H += hop_z * np.exp(1j*kz) + np.conjugate(np.transpose(hop_z)) * np.exp(-1j*kz)
    return H

# 计算能带
k_points = np.linspace(0, np.pi, 100)
energies = [np.linalg.eigvalsh(np.array(hamiltonian_k(kx, 0, 0))) for kx in k_points]

# 绘图
plt.figure(figsize=(8, 5))
for band in range(4):
    plt.plot(k_points, [e[band] for e in energies])
plt.xlabel(r"$k_x$ (Γ → X)", fontsize=12)
plt.ylabel("Energy", fontsize=12)
plt.title("3D Topological Insulator Band Structure")
plt.grid(ls="--", alpha=0.5)
plt.show()