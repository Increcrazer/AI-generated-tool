import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom, poisson

# 设置参数
lambda_val = 5
n_values = [10, 20, 50, 100]
k = np.arange(0, 16)  # 横坐标范围

# 创建子图
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

for idx, n in enumerate(n_values):
    p = lambda_val / n
    
    # 计算二项分布概率
    binom_probs = binom.pmf(k, n, p)
    
    # 计算泊松分布概率（作为参考）
    poisson_probs = poisson.pmf(k, lambda_val)
    
    # 画图
    ax = axes[idx]
    ax.bar(k - 0.15, binom_probs, width=0.3, label=f'Binomial (n={n}, p={p:.3f})', 
           alpha=0.7, color='skyblue', edgecolor='black')
    ax.bar(k + 0.15, poisson_probs, width=0.3, label=f'Poisson (λ={lambda_val})', 
           alpha=0.7, color='salmon', edgecolor='black')
    
    ax.set_xlabel('k')
    ax.set_ylabel('Probability')
    ax.set_title(f'n = {n}, p = λ/n = {lambda_val}/{n} = {p:.3f}')
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_xticks(k[::2])

plt.suptitle(f'Binomial Approximation to Poisson (λ = {lambda_val})', fontsize=14, y=1.02)
plt.tight_layout()
plt.show()
