# BI-GY-7453-Algorithms-and-Data-Structures
BI-GY 7453 Algorithms and Data Structures course
final project



## Algorithm

### Introduction

$X$ is an $n \times p$ matrix, with $n$ observations and $p$ features. 

$$\mathbf{X}_{n,p} = \left[\begin{array}
{rrrr}
X_{11} & X_{12} & ... & X_{1p}\\
X_{21} & X_{22} & ... & X_{2p}\\
... & ... & ... & ... \\
X_{n1} & X_{n2} & ... & X_{np}
\end{array}\right] = (X_{\cdot 1}, X_{\cdot 2},..., X_{\cdot p})$$

Assume that the $n$ observations belong to $R$ unknown and non-overlape classes. 

Assume that the $p$ features belong to $C$ unknown and non-overlap calsses. 

Let $A = \{a_{ij}\} \in \mathbb{R}^{n \times p}$, where $A$ matrix is the estimated the cluster labels/pattern of matrix $X$

$\mathbf{\mu} = \{\mu_{r,c}\} \in \mathbb{R}^{R \times C}$

$\mu_{r,c} = \frac{1}{|R||C|} \sum_{i \in R, j \in C} x_{ij}$


What the X and A matrix looks like? 

![image](https://github.com/sakuramomo1005/BI-GY-7453-Algorithms-and-Data-Structures/blob/master/File/xa_matrix.gif)
