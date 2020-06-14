---
title: Standard Form Of Equation/二阶偏微分方程的标准型
author: Huang
date: '2020-06-14'
slug: standard-form-of-equation
markup: mmark
categories:
  - Course
tags:
  - PDE
---

# 前言

在本科学习数学物理方程这门课程时，二阶方程的特征理论与分类是绕不开的内容。

下面将分别对含有两个自变量的二阶线性偏微分方程以及主部具有常系数的多个自变量的二阶线性偏微分方程，就其如何变换为标准型做一下基本方法的总结，并为了提高计算速度，加入了使用矩阵$J^3,J_0$ 以及合同矩阵 $B$ 的方法。

# 两个变量的情形

## 基本方法

一般的含有两个自变量的二阶线性偏微分方程可以写成形式：
$$
au_{xx}+2bu_{xy}+c_{yy}+du_x+eu_y+gu=f
$$
其中 $a,b,c,d,e,f,g$ 都是 $x,y$ 的已知函数。

略过前面的步骤（如何求特征线，取可逆自变量变换等），假设我们已知对上述方程作可逆自变量变换：
$$
\begin{cases}
\xi = \varphi(x,y)\\
\eta = \psi(x,y)
\end{cases}
$$
接下来我们应当如何求解方程关于新自变量的形式？

一种方法：我们可以利用求导的链式法则，不停地求导再合并同类项，将方程转化为关于新变量$\xi,\eta$的形式。

另一种方法，我们可以利用矩阵记忆求导关系。

**当方程为常系数方程的时候**，我们可以将上面的链式法则求导过程机械化，利用雅可比矩阵 $J$ 构造 $J^3$ ，较快得到结果。下面进行阐述。

注意到：
$$
\begin{bmatrix}
u_x\\
u_y
\end{bmatrix}
=
\begin{bmatrix}
\xi_x & \eta_x\\
\xi_y & \eta_y
\end{bmatrix}
\begin{bmatrix}
u_{\xi}\\
u_{\eta}
\end{bmatrix}
= J^{\mathrm T}
\begin{bmatrix}
u_{\xi}\\
u_{\eta}
\end{bmatrix}
$$
其中 $J$ 为变换的雅可比矩阵。并且有
$$
\begin{bmatrix}
u_{xx}\\
u_{xy}\\
u_{yy}
\end{bmatrix}
=
\begin{bmatrix}
\xi_x^2 & 2\xi_x\eta_x & \eta_x^2\\
\xi_x\xi_y & \xi_x\eta_y+\xi_y\eta_x & \eta_x\eta_y\\
\xi_y^2 & 2\xi_y\eta_y & \eta_y^2
\end{bmatrix}
\begin{bmatrix}
u_{\xi\xi}\\
u_{\xi\eta}\\
u_{\eta\eta}
\end{bmatrix}
= J^3
\begin{bmatrix}
u_{\xi\xi}\\
u_{\xi\eta}\\
u_{\eta\eta}
\end{bmatrix}
$$
并且可以看到 $J^3$ 与 $J^{\mathrm T}$ 有很大的联系，如果记住 $J^3$ ，意味着你可以不需要用链式法则繁琐地求二阶导。

然后由于
$$
au_{xx}+2bu_{xy}+c_{yy}=
\begin{bmatrix}
a & 2b & c
\end{bmatrix}
\begin{bmatrix}
u_{xx}\\
u_{xy}\\
u_{yy}
\end{bmatrix}
= 
\begin{bmatrix}
a & 2b & c
\end{bmatrix}
J^3
\begin{bmatrix}
u_{\xi\xi}\\
u_{\xi\eta}\\
u_{\eta\eta}
\end{bmatrix}
$$

$$
du_x+eu_y
=
\begin{bmatrix}
d & e
\end{bmatrix}
\begin{bmatrix}
u_x\\
u_y
\end{bmatrix}
=
\begin{bmatrix}
d & e
\end{bmatrix}
J^{\mathrm T}
\begin{bmatrix}
u_{\xi}\\
u_{\eta}
\end{bmatrix}
$$

即可得到 $u_{\xi\xi},u_{\xi\eta},u_{\eta\eta},u_{\xi},u_{\eta}$ 的系数，其余项不变，标准型即可求得。

**当方程为变系数方程的时候**，对原方程二阶导的转化时可能会出现对新变量的一阶导，具体关系为：
$$
\begin{bmatrix}
u_{xx}\\
u_{xy}\\
u_{yy}
\end{bmatrix}
= J^3
\begin{bmatrix}
u_{\xi\xi}\\
u_{\xi\eta}\\
u_{\eta\eta}
\end{bmatrix}
+
\begin{bmatrix}
\xi_{xx} & \eta_{xx}\\
\xi_{xy} & \eta_{xy}\\
\xi_{yy} & \eta_{yy}
\end{bmatrix}
\begin{bmatrix}
u_{\xi}\\
u_{\eta}
\end{bmatrix}
= J^3
\begin{bmatrix}
u_{\xi\xi}\\
u_{\xi\eta}\\
u_{\eta\eta}
\end{bmatrix}
+
J_0
\begin{bmatrix}
u_{\xi}\\
u_{\eta}
\end{bmatrix}
$$
原变量方程的二阶导数项转化为
$$
au_{xx}+2bu_{xy}+c_{yy}
= 
\begin{bmatrix}
a & 2b & c
\end{bmatrix}
J^3
\begin{bmatrix}
u_{\xi\xi}\\
u_{\xi\eta}\\
u_{\eta\eta}
\end{bmatrix}
+
\begin{bmatrix}
a & 2b & c
\end{bmatrix}
J^0
\begin{bmatrix}
u_{\xi}\\
u_{\eta}
\end{bmatrix}
$$
一阶导数的转化跟前面一样。

## 举例

方程：
$$
4u_{xx}+5u_{xy}+u_{yy}+u_x+u_y+2=0
$$
容易知道，这是一个双曲型方程，并可作变换：
$$
\begin{cases}
\xi = y-x\\
\eta = y-\dfrac x4
\end{cases}
$$
我们由上面的理论，先求
$$
J^{\mathrm T}
=
\begin{bmatrix}
\xi_x & \eta_x\\
\xi_y & \eta_y
\end{bmatrix}
=
\begin{bmatrix}
-1 & -\dfrac14\\
1 & 1
\end{bmatrix}
$$
由 $J^{\mathrm T}$ 与 $J^3$ 之间的关系，可以直接得到
$$
J^3
=
\begin{bmatrix}
1 & \dfrac12 & \dfrac1{16}\\
-1 & -\dfrac54 & -\dfrac14\\
1 & 2 & 1
\end{bmatrix}
$$
进而可以得到 $u_{\xi\xi},u_{\xi\eta},u_{\eta\eta}$ 的系数为
$$
\begin{bmatrix}
4 & 5 & 1
\end{bmatrix}
\begin{bmatrix}
1 & \dfrac12 & \dfrac1{16}\\
-1 & -\dfrac54 & -\dfrac14\\
1 & 2 & 1
\end{bmatrix}
= 
\begin{bmatrix}
0 & -\dfrac94 & 0
\end{bmatrix}
$$
$u_{\xi},u_{\eta}$ 的系数为：
$$
\begin{bmatrix}
1 & 1
\end{bmatrix}
\begin{bmatrix}
-1 & -\dfrac14\\
1 & 1
\end{bmatrix}
=
\begin{bmatrix}
0 & \dfrac 34
\end{bmatrix}
$$
所以方程标准型为：
$$
-\frac94u_{\xi\eta}+\frac34u_\eta+2=0
$$
整理一下：
$$
u_{\xi\eta}-\frac13u_{\eta}-\frac89=0
$$
这就结束了。

实际上这个例子举得并不好，然而我懒得改了。。

# 多个变量的情形

## 基本方法

仅考虑主部具有常系数的多个变量的二阶线性偏微分方程：
$$
\sum_{i,j=1}^na_{ij}\frac{\partial^2u}{\partial x_i\partial x_j}
+\sum_{i=1}^nc_i(x_1,\cdots,x_n)\frac{\partial u}{\partial x_i}
+b(x_1,\cdots,x_n)u
=f(x_1,\cdots,x_n)
$$
其中 $a_{ij}=a_{ji}$ 为常数。

同样地，我们略过前面将特征二次型$\mathcal D$化为标准型的步骤。假设我们求得了矩阵 $A=\{a_{ij}\}$ 的合同矩阵 $B$，令 $\alpha=B\beta$ ，得
$$
\mathcal D=\alpha^{\mathrm T}A\alpha=\beta^{\mathrm T}\Lambda\beta
$$
此时特征二次型$\mathcal D$化为了标准型。其中 $\Lambda=\mathrm{diag}\{\lambda_1\cdots,\lambda_n\},\lambda_i=0\ or\ 1$ 为特征二次型的标准型。

那么此时我们可以确定方程的类型，并且确定了方程二阶导数项的系数就是 $\lambda_1\cdots,\lambda_n$ 。

然后我们做自变量变换：
$$
\mathbf Y = B^{\mathrm T}\mathbf X
$$
这时我们只需要确定标准型的一阶导数项的系数就可以了。

容易知道，矩阵 $B$ 是数值矩阵，也即上面的变换，对于每个 $y_i$ ，它都是关于 $x_1,\cdots,x_n$ 的线性函数。于是容易验证：
$$
B =
\begin{bmatrix}
\dfrac{\partial y_1}{\partial x_1} & \cdots & \dfrac{\partial y_n}{\partial x_1}\\
\vdots & & \vdots\\
\dfrac{\partial y_1}{\partial x_n} & \cdots & \dfrac{\partial y_n}{\partial x_n}
\end{bmatrix}
$$
所以
$$
\begin{bmatrix}
u_{x_1}\\
\vdots\\
u_{x_n}
\end{bmatrix}
=B
\begin{bmatrix}
u_{y_1}\\
\vdots\\
u_{y_n}
\end{bmatrix}
$$
假设原方程一阶导数项的系数为 $C=[c_1,\cdots,c_n]$ ，则标准型的一阶导数项的系数为 $CB$ 。

## 举例

方程：
$$
u_{xx}-4u_{xy}+2u_{xz}+4u_{yy}+u_{zz}+2u_x+u_y-3u=0
$$
可以求得矩阵
$$
B=
\begin{bmatrix}
1 & \dfrac12 & \dfrac32\\
0 & \dfrac12 & \dfrac12\\
0 & \dfrac12 & -\dfrac12
\end{bmatrix}
$$
与 $\Lambda=\mathrm{diag}\{1,1,-1\}$ ，标准型的二阶导数项系数已经确定。

当作自变量变换：
$$
\begin{bmatrix}
\xi\\
\eta\\
\zeta
\end{bmatrix}
= B^{\mathrm T}
\begin{bmatrix}
x\\
y\\
z
\end{bmatrix}
$$
标准型一阶导数项的系数为
$$
CB
=
\begin{bmatrix}
2 & 1 & 0
\end{bmatrix}
\begin{bmatrix}
1 & \dfrac12 & \dfrac32\\
0 & \dfrac12 & \dfrac12\\
0 & \dfrac12 & -\dfrac12
\end{bmatrix}
=
\begin{bmatrix}
2 & \dfrac32 & \dfrac72
\end{bmatrix}
$$
其他项不变。所以立即可得方程的标准型为：
$$
u_{\xi\xi}+u_{\eta\eta}-u_{\zeta\zeta}+2u_{\xi}+\frac32u_{\eta}+\frac72u_{\zeta}-3u=0
$$

# 补充说明

最后求得的标准型的系数中可能含原变量，需要进一步转化，即利用自变量之间的关系将原变量转化为新变量。