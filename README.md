# 等离子体鞘套散射矩阵计算程序

## 文件说明
- SMM.m: 散射矩阵法(Scattering Matrix Method, SMM)计算多层介质平板反射系数、透射系数算法实现
- ne2wp.m: 等离子体密度到等离子体角频率的转换函数
- wv2epsr.m: 微波角频率&等离子体参数到等离子体相对介电常数的转换函数
- main.m: 主程序入口，对双高斯电子密度分布、均匀碰撞频率分布等离子体鞘套反射系数进行CPU多核并行计算
- ParWaiter.m: 实现CPU多核并行计算进度统计与可视化

## 使用方法
1. 打开 main.m 文件
2. 根据需要修改以下参数：
   - nemax: 电子密度峰值
   - la: 电子密度峰值位置
   - sigma1, sigma2: 双高斯分布展宽
   - ve: 平均碰撞角频率
   - zmax: 等离子体分布范围(-∞ ~ zmax为空气, 0 ~ zmax为等离子体, zmax ~ +∞为空气)
   - k: 网格细分比例(每层最大厚度 = lambda / k)
3. 运行程序，将得到：
   - 等离子体密度分布图
   - 反射系数复平面图

## 技术特点
- 使用Redheffer星积组合散射矩阵
- 使用parfor循环实现并行计算加速

