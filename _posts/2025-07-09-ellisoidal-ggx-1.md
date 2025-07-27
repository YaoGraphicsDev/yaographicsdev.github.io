---
layout: post
title:  "The Ellipsoidal Nature of GGX distribution - Part 1"
date:   2025-07-09 16:05:00 +0800
use_math: true
---

We first have a quick recap of how GGX distribution is formulated.

According to the reasoning given by Heitz[^Heitz2014], the formulation of the GGX distribution of normals $D$ and shadowing term $G_1$ can be simplified to these logical steps:

1. Introduce a function $P_{22}(p, q)$ to characterize the distribution of slopes of microfacets. $p$ and $q$ are the slopes of a micorfacet in $x$ and $y$ direction. 
2. A change of domain from slope space to micorfacet normal space constructs the expression of distribution of normals $D(\mathbf{m})$:

   $$
   \begin{equation}
   D(\mathbf{m}) = \frac{P_{22}(p, q)}{\cos^4\theta_m}
   \end{equation}
   $$

3. By assuming the Smith model, where the heights and slopes of microfacets are uncorrelated, the effect of inter microfacet masking takes the generalized form of 
   
    $$
    \begin{eqnarray}
    G_1(\mathbf{v}) &=& \frac{1}{1+\Lambda(\mathbf{v})} \\
    \Lambda(\mathbf{v}) &=& \frac{1}{\mu}\int_\mu^\infty(q-\mu)P_2(q)dq \\
    P_2(q) &=& \int_{-\infty}^{\infty}P_{22}(p, q)dp
    \end{eqnarray} 
    $$

    where $\mu = \lvert \cot\theta_v \rvert $.

The choice of $P_{22}$ uniquely determines the expressions for both $D$ and $G_1$. If we choose to have the microfacet slopes to follow Gaussian distribution in slope space, then by substituting $P_{22}$ with a two dimensional Gaussian distribution, we arrive at Beckmann distribution for $D$ and $G_1$.

Another choice of $P_{22}$ is 

$$
\begin{equation}
P_{22} = \frac{1}{\pi \alpha_x \alpha_y \left( 1 + \frac{p^2}{\alpha_x^2} + \frac{q^2}{\alpha_y^2} \right)^2} \label{eq:p_22}
\end{equation}
$$

This formualtion will lead us to what Walter refers to as the GGX distribution for $D$ and $G_1$ in his famous GGX paper[^Walter2007].

Even though logically $P_{22}$ comes before $D$ and $\Lambda$, the GGX paper does not explicitly present \eqref{eq:p_22}. It gives $D$ and $G_1$ directly, calls them the GGX distribution and shows GGX more closely fits measured data than Beckmann does.

Up until this point, GGX was treated as yet another distribution of microfacets, albeit it was better at describing real life surfaces under certain circumstances.

Then in 2018, another article by Heitz[^Heitz2018] proposed a new method for sampling visible normals from GGX distribution, which is commonly referred to as VNDF sampling. In practice, this method samples the 2D projection of an ellipsoid -- a routine fundamentally different from the generic inversion method for sampling a PDF. As the authors put it, "This routine leverages the property that GGX is the distribution of normals of a truncated ellipsoid". This means the GGX distribution is unique compared to other distributions, such as Beckmann, as it has a clear geometric meaning beyond an arithmetic derivation from $P_{22}$.

This blog post aims to give an understanding of the ellipsoidal nature of GGX distribution. Then in part two we will explore the implications of this observation, specifically, the geometric intuition behind $G_1$ term and VNDF sampling method.

# Distribution of Normals on an Ellipsoidal Surface

A valid distribution of normals $D(\mathbf{m}, \mathbf{u})$ along an arbitrary direction $\mathbf{u}$ must satisfy:

$$
\begin{equation}
\int_\Omega D(\mathbf{m}, \mathbf{u})(\mathbf{m}\cdot\mathbf{u})d\mathbf{m} = 1\label{eq:NDF_norm}
\end{equation}
$$

We have to map surface normal to something to get any further. We choose to map surface patch:

$$
\begin{equation}
A^{\perp}(\mathbf{u}) = \int_S (\mathbf{m}\cdot\mathbf{u})dA \label{eq:A_perp}
\end{equation}
$$

$A^{\perp}(\mathbf{u})$ is the projected area of the ellipsoidal surface in direction $\mathbf{u}$. $\Omega$ is the semi-spherical domain over the projection plane. $S$ is the domain of surface above the projection plane.

\eqref{eq:NDF_norm} and \eqref{eq:A_perp} allow us to map a surface patch to its normal direction. This mapping is referred to as the Gauss map in differential geometry. The Jacobian of this mapping exists and is equal to the Gaussian curvature $K$. Thus we can safely perform a change of variables from domain $S$ to $\Omega$:

$$
\begin{equation}
D(\mathbf{m}, \mathbf{u})(\mathbf{m}\cdot\mathbf{u})d\mathbf{m} = \frac{1}{A^{\perp}(\mathbf{u})}(\mathbf{m}\cdot\mathbf{u})dA
\end{equation}
$$

$$
\begin{equation}
D(\mathbf{m}, \mathbf{u}) = \left(\frac{d\mathbf{m}}{dA}\right)^{-1} = \frac{1}{A^{\perp}(\mathbf{u})K_g(\mathbf{m})} \label{eq:d-generic}
\end{equation}
$$

The term $\frac{d\mathbf{m}}{dA}$ is the Gaussian curvature $K$ in the form of the Gauss map Jacobian. $K_g({\mathbf{m}})$ is $K$ written as a function of surface normal $\mathbf{m}$. In the context of microfacet theory, we are interested in the case of $\mathbf{u} = \mathbf{n}$ where $\mathbf{n}$ is the normal of the macro surface. Introduce a specialize $D(\mathbf{m})$:

$$
\begin{equation}
\int_\Omega D(\mathbf{m})(\mathbf{m}\cdot\mathbf{n})d\mathbf{m} = 1
\end{equation}
$$

$$
\begin{equation}
D(\mathbf{m}) = \frac{1}{A^{\perp}(\mathbf{n})K_g(\mathbf{m})} \label{eq:d-normalized}
\end{equation}
$$

We arrive at the same conclusion as Eq.(9) in Walter's 2016 article[^Walter2016].

Calling $D(\mathbf{m})$ a **distribution** is not entirely accurate because the term that is normalized is actually $D(\mathbf{m})(\mathbf{m}\cdot\mathbf{n})$. But in compliance with convention and to differentiate between these two terms, we will refer to $D(\mathbf{m})$ as the NDF and $D(\mathbf{m})(\mathbf{m}\cdot\mathbf{n})$ as the PDF (short for probability distribution function).

# Gaussian Curvature of an Ellipsoidal Surface

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

The projected area of an ellipsoid in an arbitrary direction $\mathbf{u}$ as given by Eq.(26) in the Walter 2016 article[^Walter2016]:

$$
\begin{equation}
A^{\perp}(\mathbf{u}) = \frac{\pi \left\lVert \mathbf{A}\mathbf{u} \right\rVert}{\lvert \mathbf{A} \rvert} = \pi abc \left\lVert \mathbf{A}\mathbf{u} \right\rVert \label{eq:A_perp_expression}
\end{equation}
$$

A Math StackExchange answer[^ProjClosedForm] derives a general closed-form expression for this projected area. Substituting $Q_e = \mathbf{A}^\mathsf{T}\mathbf{A}$, light direction $\mathbf{d}_0 = \mathbf{u}$ and projection plane normal $\mathbf{n} = \mathbf{u}$ into the shadow expression yields a result identical to \eqref{eq:A_perp_expression}.

Substitute \eqref{eq:k_g} and \eqref{eq:A_perp_expression} back into \eqref{eq:d-normalized}:

$$
\begin{equation}
D(\mathbf{m}, \mathbf{u}) = \frac{abc}{\pi \left\lVert \mathbf{A}\mathbf{u} \right\rVert \left(\mathbf{m}^{\mathsf{T}} \mathbf{A}^{-1} \mathbf{A}^{-\mathsf{T}}\mathbf{m}\right)^2}
\end{equation}
$$

Conventionally a macrosurface lies in the $xOy$ plane so we only care about the NDF when $\mathbf{u}$ takes the direction of $+z$ where $\mathbf{u} = \mathbf{n} = (0, 0, 1)^{\mathsf{T}}$. Also substitute $a$ and $b$ with reciprocal of directional roughness $\alpha_x^{-1}$ and $\alpha_y^{-1}$. Assuming there is no scaling in the $z$ direction, we have $c=1$. $\mathbf{m} = (m_x, m_y, m_z)^\mathsf{T}$:

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

We arrive at GGX NDF. The substition of $\alpha_x^{-1}$ and $\alpha_y^{-1}$ implies that the smoother the microfacets are in one direction, the more the equivalent ellipsoid is stretched in that direction. It is also important to note that although $D(\mathbf{m})$ is presented as a function of $\theta_m$ and $\phi_m$, it is a distribution of normal **directions**, not a distribution of the two spherical angles. A change of domain is required to accommodate this conversion:

$$
\begin{equation}
D(\theta_m, \phi_m) = D(\mathbf{m}) \left| \frac{d\mathbf{m}}{d\theta_m d\phi_m} \right| = D(\mathbf{m})\sin\theta_m
\end{equation}
$$


# References
[^Heitz2014]: [Eric Heitz. Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs](https://jcgt.org/published/0003/02/03/paper.pdf)
[^Walter2007]: [Bruce Walter. Microfacet Models for Refraction through Rough Surfaces](https://www.graphics.cornell.edu/~bjw/microfacetbsdf.pdf)
[^Heitz2018]: [Eric Hritz. Sampling the GGX Distribution of Visible Normals](https://jcgt.org/published/0007/04/01/paper.pdf)
[^Walter2016]: [Bruce Walter. The Ellipsoid Normal Distribution Function](https://www.cs.cornell.edu/Projects/metalappearance/SupplementalEllipsoidNDF.pdf)
[^ProjClosedForm]: ["user948761. Closed form expression for shadow area of an ellipsoid subject to uniform direction light"](https://math.stackexchange.com/a/4799289)
