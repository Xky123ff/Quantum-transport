import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 设置 Matplotlib 支持中文显示
plt.rcParams['font.sans-serif'] = ['SimHei']  # 指定默认字体为黑体
plt.rcParams['axes.unicode_minus'] = False  # 解决保存图像是负号'-'显示为方块的问题

def calculate_W_mb(N, g_list, n_list):
    """计算玻耳兹曼统计的微观状态数"""
    term = 1.0
    for g, n in zip(g_list, n_list):
        if n < 0: return 0
        term *= (g**n) / math.factorial(n)
    return math.factorial(N) * term

def calculate_W_be(g_list, n_list):
    """计算玻色-爱因斯坦统计的微观状态数"""
    term = 1.0
    for g, n in zip(g_list, n_list):
        if n < 0: return 0
        term *= math.factorial(n + g - 1) / (math.factorial(n) * math.factorial(g - 1))
    return term

def calculate_W_fd(g_list, n_list):
    """计算费米-狄拉克统计的微观状态数"""
    term = 1.0
    for g, n in zip(g_list, n_list):
        # 泡利不相容原理：粒子数不能超过能级简并度
        if n > g or n < 0:
            return 0
        term *= math.factorial(g) / (math.factorial(n) * math.factorial(g - n))
    return term

# --- 系统参数设定 ---
N_particles = 10  # 总粒子数
g_levels = [8, 8]  # 两个能级的简并度 g1, g2

# --- 遍历所有可能的分布并计算 ---
results = []
for n1 in range(N_particles + 1):
    n2 = N_particles - n1
    n_list = [n1, n2]
    
    w_mb = calculate_W_mb(N_particles, g_levels, n_list)
    w_be = calculate_W_be(g_levels, n_list)
    w_fd = calculate_W_fd(g_levels, n_list)
    
    results.append({
        "Distribution (n1, n2)": f"({n1}, {n2})",
        "W_MB": w_mb,
        "W_BE": w_be,
        "W_FD": w_fd,
        "n1": n1 # 用于绘图
    })

# --- 结果展示 ---
# 1. 使用 Pandas 创建表格
df = pd.DataFrame(results)
print("微观状态数随粒子分布的变化:")
print(df[["Distribution (n1, n2)", "W_MB", "W_BE", "W_FD"]].to_string(index=False))

# 2. 使用 Matplotlib 绘图
fig, ax = plt.subplots(figsize=(12, 7))

ax.plot(df['n1'], df['W_MB'], 'o-', label='玻耳兹曼 (MB)')
ax.plot(df['n1'], df['W_BE'], 's--', label='玻色-爱因斯坦 (BE)')
ax.plot(df['n1'], df['W_FD'], '^-.', label='费米-狄拉克 (FD)')

ax.set_xlabel('分布在第一个能级的粒子数 $n_1$')
ax.set_ylabel('微观状态数 $W$ (对数坐标)')
ax.set_title(f'N={N_particles}, g=[{g_levels[0]},{g_levels[1]}] 系统微观状态数随分布的变化')
ax.set_xticks(df['n1'])
ax.set_xticklabels(df["Distribution (n1, n2)"], rotation=45)
ax.set_yscale('log') # 使用对数坐标轴，因为W_MB的值远大于其他
ax.legend()
ax.grid(True, which="both", ls="--", linewidth=0.5)

plt.tight_layout()
plt.show()