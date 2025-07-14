---
layout: post
title:  "The Ellipsoidal Nature of GGX distribution - Part 1"
date:   2025-07-09 16:05:00 +0800
use_math: true
---

We first have a quick recap of how GGX distribution is formulated.

According to the reasoning given by the article by Heitz[^a], the formulation of the GGX distribution of normals $D$ and shadowing term $G_1$ can be simplified to these logical steps

1. Assuming the Smith model $\Lambda$
2. $P^{22}$ ...

If we choose to have the microfacets slopes follow Gaussian distribution in slope space, in other words, substitute $P^{22}$ with a two dimentional Gaussian distribution, we arrive at Beckmann distribution for $D$ and $\Lambda$.

Another choice of $P^{22}$ is ....


This distribution will lead us to what Walter referred to as the GGX distribution for $D$ and $G_1$ in his famous GGX paper[^b].

However, even though logically $P^{22}$ determines $D$ and $\Lambda$, $\Lambda$ and $G_1$, which means $P^{22}$ should come first, the GGX paper does not explicitly present equation (). It gives $D$ and $G_1$ directly, call them the GGX distribution and shows GGX more closely fits measured data than Beckmann does.

Up until this point, GGX was treated as yet another distribution of microfacets, albeit it was better at describing real life surfaces under certain circumstances.


Then in 2018, another article by Heitz(reference) proposed a new method for sampling visible normals from GGX distribution, which is commonly referred to as VNDF sampling. In practice, this method samples the 2D projection of an ellipsoid -- a routine fundamentally different from the generic inversion method for sampling a PDF. As the authors put it, "This routine leverages the property that GGX is the distribution of normals of a truncated ellipsoid". This means the GGX distribution is unique compared to other distributions, such as Beckmann, as it has a clear geometric meaning behind it.

This is curious. GGX, which was originally proposed by Walter without much explanation of incentive, has a $P^{22}$ that 


This blog post aims to give a clear understanding of the ellipsoidal nature of GGX distribution. (What does part I cover, what does part II cover?)

# Distribution of Normals on an Ellipsoidal Surface

A valid distribution of normals $D({\mathbf{m}})$ must satisfy these normalization constraint:

$$
\begin{equation}
\int_\Omega D(\mathbf{m})(\mathbf{m}\cdot\mathbf{n})d\mathbf{m} = 1
\end{equation}
$$

$$
\begin{equation}
\int_S (\mathbf{m}\cdot\mathbf{n})dA = A^{\perp}(\mathbf{n})
\end{equation}
$$

$A^{\perp}(\mathbf{n})$ is the projected area of the ellipsoidal surface in direction $\mathbf{n}$. $\Omega$ is the semi-spherical domain over the projection plane. $S$ is the domain of the surface above the projection plane (might need to explain n, m in more details). 

Calling $D({\mathbf{m}})$ a **distribution** is not entirely accurate because the term that is normalized is actually $D(\mathbf{m})(\mathbf{m}\cdot\mathbf{n})$. But in compliance with convention and to differentiate between these two terms, we will refer to $D(\mathbf{m})$ as the NDF and $D(\mathbf{m})(\mathbf{m}\cdot\mathbf{n})$ as the PDF (short for probability distribution function).

The two constraints allows us to map a surface patch to its normal direction. This mapping is referred to as the Gauss map in differential geometry. The Jacobian of this mapping exists and is equal to the Gaussian curvature $K$. Thus we can safely perform a change of variables from domain $S$ to $\Omega$:

$$
\begin{eqnarray}
D(\mathbf{m})(\mathbf{m}\cdot\mathbf{n})d\mathbf{m} &=& \frac{1}{A^{\perp}(\mathbf{n})}(\mathbf{m}\cdot\mathbf{n})dA \\
\Rightarrow D(\mathbf{m}) &=& \frac{1}{A^{\perp}(\mathbf{n})}\left(\frac{d\mathbf{m}}{dA}\right)^{-1} \\
&=& \frac{1}{A^{\perp}(\mathbf{n})K_g(\mathbf{m})} \label{eq:surface-normal-pdf}
\end{eqnarray}
$$

The term $\frac{d\mathbf{m}}{dA}$ is the Gaussian curvature $K$ in the form of the Gauss map Jacobian. $K_g({\mathbf{m}})$ is $K$ written as a function of surface normal $\mathbf{m}$. We arrive at the same conclusion as Eq.(9) in this article(add proper reference).

For an ellipsoid defined by the equation

$$
\begin{equation}
\frac{x^2}{a^2} + \frac{y^2}{b^2} + \frac{z^2}{c^2} = 1
\end{equation}
$$

Write it as a generic quadratic surface that consists of all points $\mathbf{p}$ that satisfies:

$$
\begin{equation}
f = \mathbf{p}^\mathsf{T}\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{p} = 1
\end{equation}
$$

where 

$$
A = \left[\begin{matrix}
a^{-1} & 0 & 0 \\
0 & b^{-1} & 0 \\
0 & 0 & c^{-1}
\end{matrix}\right]
$$

Gaussian curvature via Hessian matrix $H$ and gradient $\nabla f$, known as the Brioschi formula in 3D space:

$$
\begin{eqnarray}
K &=& -{\frac{
    \begin{vmatrix}
    H & \nabla f\\
    \nabla f^{\mathsf{T}} & 0
    \end{vmatrix}}{\left\lVert \nabla f \right\rVert ^{4}}} \\ 
&=& -{\frac{
    \begin{vmatrix}
    f_{xx} & f_{xy} & f_{xz} & f_{x}\\
    f_{xy} & f_{yy} & f_{yz} & f_{y}\\
    f_{xz} & f_{yz} & f_{zz} & f_{z}\\
    f_{x} & f_{y} & f_{z} & 0\\
    \end{vmatrix}
    }{\left\lVert \nabla f \right\rVert^{4}}} \\
&=& \frac{16}{a^2b^2c^2} \left( \frac{x^2}{a^2} + \frac{y^2}{b^2} + \frac{z^2}{c^2} \right) \left\lVert \nabla f \right\rVert ^{-4} \\
&=& \frac{16}{a^2b^2c^2} \left\lVert \nabla f \right\rVert ^{-4} \label{eq:ellipsoid-k}
\end{eqnarray}
$$

Now we just need to figure out the magnitude of $\nabla f$ in terms of surface normal $\mathbf{m}$, which for an implicit surface is defined as:

$$
\begin{equation}
\mathbf{m} = \frac{\nabla f}{\left\lVert \nabla f \right\rVert} = \frac{2\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{p}}{\left\lVert \nabla f \right\rVert}
\end{equation}
$$

Surface point $\mathbf{p}$ can in turn be written as:

$$
\begin{equation}
\mathbf{p} = \frac{1}{2}\left\lVert \nabla f \right\rVert \mathbf{A}^{-1}\mathbf{A}^{-\mathsf{T}} \mathbf{m}
\end{equation}
$$

The implicit ellipsoidal surface in terms of surface normals:

$$
\mathbf{p}^\mathsf{T}\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{p} = \frac{1}{4} \left\lVert \nabla f \right\rVert^2 \mathbf{m}^\mathsf{T} \mathbf{A}^{-1} \mathbf{A}^{-\mathsf{T}}\mathbf{m} = 1
$$

Magnitude of gradient follows:

$$
\begin{equation}
\left\lVert \nabla f \right\rVert^2 = 4 \left(\mathbf{m}^{\mathsf{T}} \mathbf{A}^{-1} \mathbf{A}^{-\mathsf{T}}\mathbf{m}\right)^{-1} \label{eq:mag-gradient}
\end{equation}
$$

Substitute \eqref{eq:mag-gradient} into \eqref{eq:ellipsoid-k}, we have $K$ in terms of $\mathbf{m}$:

$$
\begin{equation}
K_g(\mathbf{m}) = \frac{\left(\mathbf{m}^{\mathsf{T}} \mathbf{A}^{-1} \mathbf{A}^{-\mathsf{T}}\mathbf{m}\right)^2}{a^2b^2c^2} \label{eq:k_g}
\end{equation}
$$

Substitute \eqref{eq:k_g} back into \eqref{eq:surface-normal-pdf}, also note that projected area $A^{\perp}(\mathbf{n}) = \pi ab$ when $\mathbf{n}$ takes the direction of $+z$:

$$
\begin{equation}
D(\mathbf{m}) = \frac{abc^2}{\pi \left(\mathbf{m}^{\mathsf{T}} \mathbf{A}^{-1} \mathbf{A}^{-\mathsf{T}}\mathbf{m}\right)^2}
\end{equation}
$$

Substitute $a$ and $b$ with reciprocal of directional roughness $\alpha_x^{-1}$ and $\alpha_y^{-1}$. Also make substitutions $c=1$ and $\mathbf{m} = (m_x, m_y, m_z)^\mathsf{T}$:

$$
\begin{equation}
D(\mathbf{m}) = \frac{1}{\pi \alpha_x \alpha_y \left( \frac{m_x^2}{\alpha_x^2}+\frac{m_y^2}{\alpha_y^2}+m_z^2 \right)^2} 
\end{equation}
$$

Finally, write surface normal in spherical coordinates $(\sin{\theta_m}\cos{\phi_m}, \sin{\theta_m}\sin{\phi_m}, \cos{\theta_m})$:

$$
\begin{equation}
D(\mathbf{m}) = \frac{1}{\pi \alpha_x \alpha_y \cos^4{\theta_m} \left(1 + \tan^2{\theta_m}\left( \frac{\cos^2{\phi_m}}{\alpha_x^2}+\frac{\sin^2{\phi_m}}{\alpha_y^2} \right) \right)^2} 
\end{equation}
$$


(also, the heavyside function)
We arrive at GGX NDF. The substition of $\alpha_x^{-1}$ and $\alpha_y^{-1}$ implies that the smoother the microfacets are in one direction, the more the equivalent ellipsoid is stretched in that direction. It is also important to note that although $D(\mathbf{m})$ is presented as a function of $\theta_m$ and $\phi_m$, it is a distribution of normal **directions**, not a distribution of the two spherical angles. A change of domain is required to make this conversion:

$$
\begin{equation}
D(\theta_m, \phi_m) = D(\mathbf{m}) \left| \frac{d\mathbf{m}}{d\theta_m d\phi_m} \right| = D(\mathbf{m})\sin\theta_m
\end{equation}
$$

# The Masking Term $G1$

Since we have established that the GGX NDF is equivalent to a truncated ellipsoidal surface NDF, we can build a mental model in which microfacets are viewed as tiny pieces broken off from the upper half of an ellipsoid and then displaced onto the $z=0$ plane. By displacing all of them onto the same plane, masking occurs between microfacets when viewed from an angle. Walter(reference) derived the masking term from a physical standpoint by comparing masking behavior to attenuation in a volumetric medium. Heitz on the other hand juggles around ....

There is also an geometric interpretation of the masking term.



# Sampling Normals

The ellipsoidal geometry also provides another insight into the process of sampling normals. Before VNDF sampling was introduced, the standard approach to sampling normals revolved around the PDF $D(\mathbf{m})(\mathbf{m}\cdot\mathbf{n})$, which did not account for the masking term.

Since we have established the ellipsoidal nature of the NDF, we now show that the sampling equations can also be derived by uniformly choosing points from the projected disk and then generating directions by projecting those points upward onto the ellipsoidal hemisphere above. In a sense, this can be seen as a generalized version of Malley's method(quote).

The equations for uniformly sampling an ellipse disk $\frac{x^2}{a^2} + \frac{y^2}{b^2} = 1$ are:

$$
\begin{eqnarray}
x &=& ar\cos{\theta} \\
y &=& br\sin{\theta}
\end{eqnarray}
$$

where 

$$
\begin{eqnarray}
r &=& \sqrt{\xi_1} \\
\theta &=& 2 \pi \xi_2
\end{eqnarray}
$$

$$
z = c\sqrt{1 - \frac{x^2}{a^2} - \frac{y^2}{b^2}} = c\sqrt{1 - \xi_1}
$$


$$
\mathbf{m} = \frac{\nabla f}{\left\lVert \nabla f \right\rVert} = \frac{2}{\left\lVert \nabla f \right\rVert} \left( \frac{x}{a^2}, \frac{y}{b^2}, \frac{z}{c^2} \right) = \frac{2}{\left\lVert \nabla f \right\rVert} \left( \frac{\sqrt{\xi_1}\cos(2\pi\xi_2)}{a}, \frac{\sqrt{\xi_1}\sin(2\pi\xi_2)}{b}, \frac{\sqrt{1-\xi_1}}{c} \right)
$$

$$
\begin{eqnarray}
\theta_m &=& \arctan\left(\frac{\sqrt{m_x^2 + m_y^2}}{m_z}\right) \\
&=& \arctan\left(c \sqrt{\frac{\xi_1}{1-\xi_1}\left(\frac{\cos^2(2\pi\xi_2)}{a^2}+\frac{\sin^2(2\pi\xi_2)}{b^2}\right)}\right) \\
\phi_m &=& \arctan\left( \frac{y_m}{x_m} \right) = \arctan\left( \frac{a}{b} \tan\left( 2\pi\xi_2 \right) \right)
\end{eqnarray}
$$


This result can be obtained by applying inversion method on the PDF:

$$
p(\theta_m, \phi_m) = D(\theta_m, \phi_m)(\mathbf{m}\cdot\mathbf{n}) = D(\mathbf{m})\sin\theta_m\cos\theta_m
$$

as demonstrated here(quote).



$$
p(x, y) = \frac{1}{\pi \alpha^2}
$$

$$
\begin{equation}
p(r_d, \phi_d) = \frac{p(x, y)}{\left| J_T \right|} = \frac{r_d}{\pi\alpha^2}
\end{equation}
$$

$$
\begin{eqnarray}
p(r_d) &=& \int_{0}^{2\pi}p(r_d, \phi_d)d\phi_d = \frac{2r_d}{\alpha^2} \\
p(\phi_d) &=& \frac{p(r_d, \phi_d)}{p(r_d)} = \frac{1}{2\pi}
\end{eqnarray}
$$

$$
\begin{eqnarray}
\xi_1 = \int p(r_d) dr_d = \frac{r_d^2}{\alpha^2} &\Rightarrow& r_d = \alpha\sqrt{\xi_1} \\
\xi_2 = \int p(\phi_d) d\phi_d = \frac{\phi_d}{2\pi} &\Rightarrow& \phi_d = 2 \pi \xi_2
\end{eqnarray}
$$

(Tangent of ellipsoid is not correct) 

$$
\tan{\theta_m} = \frac{r_d}{z_m} = \frac{r_d}{\sqrt{1 - \frac{x_m^2 + y_m^2}{\alpha^2}}} = \frac{r_d}{\sqrt{1 - \frac{r_d^2}{\alpha^2}}} = \frac{\alpha \sqrt{\xi_1}}{\sqrt{1 - \xi_1}}
$$


# Sampling Normals from a Lune

$$
\begin{eqnarray}
t_1^{\prime} &=& t_1 \\
t_2^{\prime} &=& (1-s)\sqrt{1-t_1^2}+st_2
\end{eqnarray}
$$

The Jacobian of this transformation is:

$$
\begin{equation}
J = 
\begin{vmatrix}
\frac{\partial t_1^{\prime}}{t_1} & \frac{\partial t_1^{\prime}}{t_2} \\
\frac{\partial t_2^{\prime}}{t_1} & \frac{\partial t_2^{\prime}}{t_2}
\end{vmatrix}

= \begin{vmatrix}
1 & 0 \\
0 & s
\end{vmatrix}
\end{equation}
=s
$$

What this means is that, if we sample uniformly from a unit disk, an area-preserving transformation on the samples will not change the uniform nature of the samples' distribution.


# Prove that after area-preserving transform, the new distriubtion is still a uniform distribution, by calculating jacobian of t1 and t2

We take the idea of sampling from ... one step further and sample 

$$

$$

# How to intuitively understand VNDF, comparing between sampling from above and from an angle


# Prove that 

Suface points after being transformed to unit sphere

$$
\begin{equation}
\mathbf{p}^{\prime} = \mathbf{A}\mathbf{p}
\end{equation}
$$

With viewing direction denoted as $\mathbf{v}$, surface normals on silhouette must satisfy:

$$
\mathbf{v}\cdot\mathbf{m} = 0
$$

Substitute $\mathbf{m}$, surface point on the silhouette must satisfy:

$$
0 = \mathbf{v} \cdot \left(\mathbf{A}^{\mathsf{T}}\mathbf{A}\mathbf{p}\right) = \mathbf{v}^\mathsf{T} \mathbf{A}^{\mathsf{T}}\mathbf{A}\mathbf{p} = \left( \mathbf{A}^{\mathsf{T}}\mathbf{A}\mathbf{v} \right)^\mathsf{T} \mathbf{p} = \left( \mathbf{A}^{\mathsf{T}}\mathbf{A}\mathbf{v} \right) \cdot \mathbf{p}
$$

This happens to be the equation of the plane that has $\mathbf{A}^{\mathsf{T}}\mathbf{A}\mathbf{v}$ as normal and also passes through origin (?). The silhouette is the intersection of this plane and the ellipsoid. Notice that the silhouette plane's normal direction is not necessarily aligned with the viewing direction, given an arbitrary $\mathbf{A}$.

Now we 

Transformed surface point:

$$
\mathbf{p}^{\prime} = \mathbf{A}\mathbf{p}
$$

We can prove all $\mathbf{p}^{\prime}$ constitute a unit sphere. Thus surface normal:

$$
\mathbf{m}^{\prime} = \mathbf{p}^{\prime}
$$



# References
[^a]: [Eric Heitz. Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs](https://jcgt.org/published/0003/02/03/paper.pdf)
[^b]: [Bruce Walter. Microfacet Models for Refraction through Rough Surfaces](https://www.graphics.cornell.edu/~bjw/microfacetbsdf.pdf)