import numpy as np

# 定义变量
N = 3  # DHN中的N的数量
Np = 3  # 管道的数量
Nd = 2  # 负载的数量
D = np.ones(3) * 150e-3  # 管道直径
ep = np.ones(3) * 1.25e-3  # 管道粗糙度
s = np.pi * D**2 / 4  # 管道横截面积
length = np.array([400, 400, 600])  # 管道长度
mq = np.array([2, 3, -5])  # 节点质量流量
rho = 958.4  # 水的密度 (kg/m^3) 在100度C
g = 9.81  # 重力加速度
viscosity = 0.294e-6  # 温度为100度时的动力粘度 unit:m2/s
A = np.array([[1, -1, 0], [0, 1, 1], [-1, 0, -1]])  # 网络节点-弧关联矩阵
B = np.array([[1, 1, -1]])  # 环路关联矩阵

# 由水力计算得到的管段流量，此处作为给定值
m = np.array([3, 1, 2])

# 热力计算数据
lamda = 0.2  # 管道导热系数 W/(m*K)
cp = 4182  # 水的比热容 J/(kg*K)
Load_num = [0, 1]  # 负荷节点索引列表
Ak = np.array([[-1, 1, 0], [1, -1, 2], [0, 2, -1]])  # 检索两节点间连接的管道编号
Ts = np.zeros(Np)  # 节点供水温度初始化
Ts[-1] = 100  # 给定热源供水温度 100 ℃
To = np.zeros(Np)  # 节点回水温度计算初始化
To[:2] = 50  # 给定节点回水温度 50 ℃
Ta = 10  # 给定环境温度 10 ℃
Ts = Ts - Ta  # 减去环境温度
To = To - Ta  # 减去环境温度

# 负荷节点供水温度计算，算法流程见教材图2.6
Cs = np.zeros((Nd, Nd))  # 初始化
bs = np.zeros(Nd)  # 初始化
for i in range(Nd):  # 节点i
    numberOfOnesInRow2 = np.sum(A[i, :] == 1)  # 确认有几条支路汇入
    if numberOfOnesInRow2 == 1:  # 判断为独立节点
        Cs[i, i] = 1
        k = np.where(A[i, :] == 1)[0][0]  # 确定注入节点i的是哪条管道k
        j = np.where(A[:, k] == -1)[0][0]  # 确定和管道k、节点i连接的是哪个节点j
        bs[i] = Ts[j] * np.exp(-lamda * length[k] / cp / m[k])
    else:  # 混合节点
        K = np.where(A[i, :] == 1)[0]  # 确定注入节点i的所有管道k
        Cs[i, i] = np.sum(m[K])
        for j in range(Np):
            if j == i:
                # i==j跳过
                continue
            if j in Load_num:  # 判断节点j是否为负荷节点
                k = Ak[i, j]  # 确定节点i和j之间连接的管道
                Cs[i, j] = -m[k] * np.exp(-lamda * length[k] / cp / m[k])
            else:
                k = Ak[i, j]
                bs[i] += m[k] * Ts[j] * np.exp(-lamda * length[k] / cp / m[k])

# 求解线性方程组
sol1 = np.linalg.solve(Cs, bs)
sol1 = sol1 + Ta  # 加回环境温度
print('节点1供水温度：{:.3f}\n节点2供水温度：{:.3f}'.format(sol1[0], sol1[1]))

"""
请继续补充回水温度计算代码，算法流程见教材图2.7
需要输出：每个负荷节点的回水温度，最终得到热源的回水温度
"""


