import kwant
import numpy as np
import matplotlib.pyplot as plt
import tinyarray

# [新] 导入 Matplotlib 的 Slider 控件
from matplotlib.widgets import Slider

# --- 矩阵定义 (和您原来的一样) ---
sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])

I_4= tinyarray.array(np.kron(sigma_0, sigma_0))
gamma_1 = tinyarray.array(np.kron(sigma_0, sigma_z))
gamma_2 = tinyarray.array(np.kron(sigma_z, sigma_z))
gamma_3 = tinyarray.array(np.kron(sigma_0, sigma_y))
gamma_4 = tinyarray.array(np.kron(sigma_0, sigma_x))

# --- 晶格定义 ---
las_val = 1
structure=kwant.lattice.square(las_val, norbs=4)
model_sym= kwant.TranslationalSymmetry(structure.vec((-1,0)),structure.vec((0,-1)))

# --- 参数化的哈密顿量函数 ---
# (与我们上次修改的一样，以便 kwant.plotter.bands 传递参数)
def onsite(site, m, B, las):
    return m/2*gamma_1 - (4*B)/las**2 * gamma_2

def hopping_x(site1, site2, g, B, las):
    return B/(las**2)*gamma_2 + 1j*g/(2*las)*gamma_3

def hopping_y(site1, site2, g, B, las):
    return B/(las**2)*gamma_2 - 1j*g/(2*las)*gamma_4

# --- model_MTI (不变) ---
def model_MTI():
    model_1= kwant.Builder(model_sym)
    model_1[structure(0, 0)] = onsite
    model_1[kwant.builder.HoppingKind((1,0), structure)]  = hopping_x
    model_1[kwant.builder.HoppingKind((0,1), structure)]  = hopping_y
    return model_1

# --- model_EG (不变) ---
def model_EG():
    onsite_leads= 150 *I_4
    hopping_leads= -150 *I_4
    model_2= kwant.Builder(model_sym)
    model_2[structure(0, 0)] = onsite_leads
    model_2[structure.neighbors()]  = hopping_leads
    return model_2

# --- make_syst (不变) ---
def make_syst(width=50,length=100,leads_model=model_EG()):
    def shape_center(site):
        x, y = site.pos
        return 0 < x <= length and 0 < y <= width
    syst= kwant.Builder()
    syst.fill(model_MTI(), shape_center,(1,1))
    width_lead = length//5
    lead0= kwant.Builder(kwant.TranslationalSymmetry((-1,0)))
    lead1= kwant.Builder(kwant.TranslationalSymmetry((0,1)))
    lead2= kwant.Builder(kwant.TranslationalSymmetry((0,1)))
    def shape_lead0(site):
        x, y = site.pos
        return 0 < y <= width
    def shape_lead1(site):
        x, y = site.pos
        return width_lead < x <= 2*width_lead
    def shape_lead2(site):
        x, y = site.pos
        return 3*width_lead < x <= 4*width_lead
    lead0.fill(leads_model, shape_lead0, (0, 1))
    lead1.fill(leads_model, shape_lead1, (width_lead+1, 0))
    lead2.fill(leads_model, shape_lead2, (3*width_lead+1, 0))
    for lead in [lead0,lead1,lead2, lead0.reversed(),lead2.reversed(),lead1.reversed()]:
        syst.attach_lead(lead)
    return syst.finalized()

# --- [关键修改] 使用 Matplotlib 控件进行交互 ---

# 1. 像之前一样，构建一个以 model_MTI 为引线的系统
syst_lead_template = make_syst(width=60, leads_model=model_MTI()).leads[0]

# 2. 定义绘图的固定参数
xlim = 0.05
ylim = 20

# 3. 创建画布 (Figure) 和主绘图区 (Axes)
#    我们增加 figsize 使窗口更大，并使用 subplots_adjust 在底部留出空间给滑块
fig, ax = plt.subplots(figsize=(7, 6))
plt.subplots_adjust(left=0.1, bottom=0.35)

# 4. 定义初始参数
m_init = 0.0
B_init = -300.0
g_init = 300.0

# 5. [核心] 定义滑块的更新函数
#    每次滑块变动时，此函数会被调用
def update(val):
    # 从滑块对象中获取当前值
    m = slider_m.val
    B = slider_B.val
    g = slider_g.val
    
    # 打包 kwant 需要的参数字典
    params = dict(m=m, B=B, g=g, las=las_val)
    
    # [重要] 清除旧的能带图
    ax.clear()
    
    # 重新绘制能带
    kwant.plotter.bands(syst_lead_template, momenta=np.linspace(-xlim, xlim, 100), 
                      ax=ax, params=params)
                      
    # [重要] 重新设置所有标签和范围 (因为 ax.clear() 会清除它们)
    ax.set_ylabel('Energy (meV)')
    ax.set_xlabel('k (nm^-1)')
    ax.set_ylim(-ylim, ylim)
    ax.set_xlim(-xlim, xlim)
    ax.set_title(f'm={m:.1f}, B={B:.1f}, g={g:.1f}')
    
    # 告诉画布重新绘制
    fig.canvas.draw_idle()

# 6. 为滑块创建新的 "Axes" (绘图区)
#    这些是放置滑块的小矩形区域 [left, bottom, width, height]
ax_m = plt.axes([0.15, 0.20, 0.65, 0.03])
ax_B = plt.axes([0.15, 0.15, 0.65, 0.03])
ax_g = plt.axes([0.15, 0.10, 0.65, 0.03])

# 7. 创建 Slider 对象
slider_m = Slider(
    ax=ax_m,
    label='m',
    valmin=-100,
    valmax=100,
    valinit=m_init, # 初始值
    valstep=1
)

slider_B = Slider(
    ax=ax_B,
    label='B',
    valmin=-500,
    valmax=500,
    valinit=B_init,
    valstep=10
)

slider_g = Slider(
    ax=ax_g,
    label='g',
    valmin=0,
    valmax=500,
    valinit=g_init,
    valstep=10
)

# 8. 将滑块的 "on_changed" 事件 链接到 我们的 "update" 函数
slider_m.on_changed(update)
slider_B.on_changed(update)
slider_g.on_changed(update)

# 9. 手动调用一次 update() 来绘制初始状态的图像
update(None)

# 10. [核心] 显示窗口。程序会在此处暂停，等待用户交互
plt.show()