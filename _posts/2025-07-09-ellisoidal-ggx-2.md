---
layout: post
title:  "The Ellipsoidal Nature of GGX distribution - Part 2"
date:   2025-07-09 16:05:00 +0800
use_math: true
---

# Sampling Normals

The ellipsoidal geometry also provides another insight into the process of sampling normals. Before VNDF sampling was introduced, the standard approach to sampling normals revolved around the PDF $D(\mathbf{m})(\mathbf{m}\cdot\mathbf{n})$, which did not account for the masking term.

We now show that the sampling equations can be derived by uniformly choosing points from the projected ellipsoid disk and then generating surface points $p$ by projecting those points upward onto the ellipsoidal hemisphere above. The normal direction $\mathbf{m}$ we sample is the normal direction at $\mathbf{p}$. In a sense, this can be seen as a generalized Malley's method of sampling cosine-weighted hemisphere.

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

Project this point on the disk to ellipsoidal surface:

$$
\begin{equation}
z = c\sqrt{1 - \frac{x^2}{a^2} - \frac{y^2}{b^2}} = c\sqrt{1 - \xi_1}
\end{equation}
$$

Surface normal at $(x, y, z)$ is:

$$
\mathbf{m} = \frac{\nabla f}{\left\lVert \nabla f \right\rVert} = \frac{2}{\left\lVert \nabla f \right\rVert} \left( \frac{x}{a^2}, \frac{y}{b^2}, \frac{z}{c^2} \right) = \frac{2}{\left\lVert \nabla f \right\rVert} \left( \frac{\sqrt{\xi_1}\cos(2\pi\xi_2)}{a}, \frac{\sqrt{\xi_1}\sin(2\pi\xi_2)}{b}, \frac{\sqrt{1-\xi_1}}{c} \right)
$$

Spherical coordinates of $\mathbf{m}$ are as follows. Also substitute $\mathbf{a}$ and $\mathbf{b}$ with $(\alpha_x)^{-1}$ and $(\alpha_y)^{-1}$: 

$$
\begin{eqnarray}
\theta_m &=& \arctan\left(\frac{\sqrt{m_x^2 + m_y^2}}{m_z}\right) \\
&=& \arctan\left(c \sqrt{\frac{\xi_1}{1-\xi_1}\left(\frac{\cos^2(2\pi\xi_2)}{a^2}+\frac{\sin^2(2\pi\xi_2)}{b^2}\right)}\right) \\
&=& \arctan\left(c \sqrt{\frac{\xi_1}{1-\xi_1}\left(\alpha_x^2\cos^2(2\pi\xi_2)+\alpha_y^2\sin^2(2\pi\xi_2)\right)}\right) \label{eq:theta_m_sample} \\
\phi_m &=& \arctan\left( \frac{y_m}{x_m} \right) = \arctan\left( \frac{a}{b} \tan\left( 2\pi\xi_2 \right) \right) = \arctan\left( \frac{\alpha_y}{\alpha_x} \tan\left( 2\pi\xi_2 \right) \right) \label{eq:phi_m_sample}
\end{eqnarray}
$$

The sampling equations for $\theta_m$ and $\phi_m$ match equations obtained via the inversion method listed here[^GraphicsGuysNotes]. One minor discrepancy concerns the parameterization of $\theta_m$: while the referenced blog post expresses $\theta_m$ as a function of $\xi_1$ and $\phi_m$, our formulation presents it as a function of $\xi_1$ and $\xi_2$. This difference can be resolved by susbtituting the following into their $\phi_m$ expression:

$$
\begin{eqnarray}
\cos\phi_m &=& \frac{\alpha_x\cos(2\pi\xi_2)}{\sqrt{\alpha_x^2\cos^2(2\pi\xi_2)+\alpha_y^2\sin^2(2\pi\xi_2)}} \\
\sin\phi_m &=& \frac{\alpha_y\sin(2\pi\xi_2)}{\sqrt{\alpha_x^2\cos^2(2\pi\xi_2)+\alpha_y^2\sin^2(2\pi\xi_2)}}
\end{eqnarray}
$$



# The Masking Term $G_1$

Since we have established that the GGX NDF is equivalent to a truncated ellipsoidal surface NDF, we can build a mental model in which microfacets are viewed as tiny pieces broken off from an upper half of an ellipsoid and then displaced onto the $z=0$ plane. By displacing all of them onto the same plane, masking occurs between microfacets when viewed obliquely. Appendix A of the GGX paper[^Walter2007] derived the masking term from a physical standpoint through analogy between inter microfacet masking behavior and attenuation in a volumetric medium. While Heitz arrived the same result through more rigorous mathematical derivation in Appendix A of his article[^Heitz2014], his approach is less intuitive.

<div style="text-align: center; margin: 40px 0;">
  <figure id="fig1">
    <img src="/images/truncated_ellipsoid.svg" alt="truncated ellipsoid showcasing shadowing term" style="width: 80%; height: auto; display: block; margin: 0 auto;">
    <figcaption>Figure 1: All surface points $\mathbf{p}$ with normals $\mathbf{m} \perp \mathbf{v}$ satisfy the plane equation $(\mathbf{A}^{\mathsf{T}}\mathbf{A}\mathbf{v}) \cdot \mathbf{p} = 0$. $\mathbf{A}^{\mathsf{T}}\mathbf{A}\mathbf{v}$ is not necessarily parallel to $\mathbf{v}$. $S_u$ represents the area of the entire projected shape, while $S_m$ denotes the darker elliptical region. </figcaption>
  </figure>
</div>

There is also a geometric interpretation of the masking term, as illustrated in [Figure 1](#fig1). When microfacets are spread across the ellipsoidal surface, as long as they are not facing away from the viewing direction v (i.e., $\mathbf{m}\cdot\mathbf{v}>0$), they will be visible. The total unmasked projected area of these microfacets is

$$
\begin{equation}
S_u = \frac{A^\perp(\mathbf{v})+A^\perp(\mathbf{n})(\mathbf{v}\cdot\mathbf{n})}{2} \label{eq:unmasked-projection}
\end{equation}
$$


After they break off and are displaced to the $z=0$ plane, the macroscopic geometry of a flat plane limits the maximum total projected area of microfacets to 
 
$$
\begin{equation}
S_m = A^\perp(\mathbf{n})(\mathbf{v}\cdot\mathbf{n}) \label{eq:masked-projection}
\end{equation}
$$
 
We can conclude by looking at the illustration that for any $\mathbf{v}$ not parallel to $\mathbf{n}$, $S_m$ is smaller than $S_u$. This reduction in the projected area of non-facing-away microfacets is the effect of inter-microfacet masking. Thus, $G_1$ when expressed as a ratio of these two projected area is:

$$
\begin{equation}
G_1(\mathbf{v}) = \frac{S_m}{S_u} = \frac{2A^\perp(\mathbf{n})}{A^\perp(\mathbf{v})+A^\perp(\mathbf{n})(\mathbf{v}\cdot\mathbf{n})}(\mathbf{v}\cdot\mathbf{n})
\end{equation}
$$

The projected area of an ellipsoid in an arbitrary direction $\mathbf{u}$ as mentioned in part one is:

$$
\begin{equation}
A^\perp(\mathbf{u}) = \frac{\pi \left\lVert \mathbf{A}\mathbf{u} \right\rVert}{\left|\mathbf{A}\right|} \label{eq:projected-area}
\end{equation}
$$

Thus:

$$
\begin{equation}
G_1(\mathbf{v}) = \frac{2\left\lVert \mathbf{A}\mathbf{n} \right\rVert}{\left\lVert \mathbf{A}\mathbf{v} \right\rVert + \left\lVert \mathbf{A}\mathbf{n} \right\rVert (\mathbf{v}\cdot\mathbf{n})} (\mathbf{v}\cdot\mathbf{n}) \label{eq:g1}
\end{equation}
$$

# VNDF sampling
With an understanding of $D$ and $G_1$ in place, we can finally derive the PDF of visible normals, a crucial step if we intend to perform importance sampling on visible normals.

According to the VNDF sampling article[^Heitz2018], VNDF sampling routine consists of the following steps:

1. Transform equivalent ellipsoid and view direction to sphere space. View direction $\mathbf{v}$ and surface normal $\mathbf{m}$ in ellipsoid space transform to $\mathbf{v}^\prime$ and $\mathbf{m}^\prime$ in sphere space.
2. Draw a sample $(t_1, t_2)$ from a uniform unit disk.
3. Map the sample onto a crescent shape. Mapped sample $(t_1^\prime, t_2^\prime)$.
4. Project that sample back to spherical surface
5. Get normal direction.
6. Transform normal direction back to ellipsoid space.

Starting from the transformed sphere space, surface normal sampled through unit disk from view direction $\mathbf{v}^\prime$ follows PDF:

$$
p_s(\mathbf{m^{\prime}}) = D_s(\mathbf{m^{\prime}}, \mathbf{v^{\prime}})(\mathbf{m^{\prime}}\cdot\mathbf{v^{\prime}})
$$

where $D_s$ denotes the NDF on unit sphere surface and

$$
\begin{eqnarray}
\mathbf{v^{\prime}} &=& \frac{\mathbf{A}\mathbf{v}}{\left\lVert \mathbf{A}\mathbf{v} \right\rVert} \\
\mathbf{m^{\prime}} &=& \frac{\mathbf{A}^{-\mathsf{T}}\mathbf{m}}{\left\lVert \mathbf{A}^{-\mathsf{T}}\mathbf{m} \right\rVert} \label{eq:m-to-m-prime}
\end{eqnarray}
$$

Next, the transformaion between $(t_1, t_2)$ and $(t_1^\prime, t_2^\prime)$ as given by Eq.(11) of the VNDF sampling article[^Heitz2018]:

$$
\begin{eqnarray}
t_1^{\prime} &=& t_1 \\
t_2^{\prime} &=& (1-s)\sqrt{1-t_1^2}+st_2
\end{eqnarray}
$$

Account for the transformation of mapping a sample from unit disk to crescent disk with a Jacobian determinant:  

$$
\begin{equation}
\frac{\partial(t_1^\prime, t_2^\prime)}{\partial(t_1, t_2)} = 
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

where $s$ is the mapping constant in sphere space:

$$
\begin{equation}
s = \frac{1+\left( \frac{\mathbf{A}\mathbf{v}}{\left\lVert \mathbf{A}\mathbf{v} \right\rVert} \cdot \frac{\mathbf{A}^{-\mathsf{T}}\mathbf{n}}{\left\lVert \mathbf{A}^{-\mathsf{T}}\mathbf{n} \right\rVert} \right)}{2} \label{eq:s-mapping}
\end{equation}
$$

Since $\mathbf{n}$ represents the normal vector of a cross-section in ellipsoid space, its transformation to sphere space must be given by the inverse transpose of $\mathbf{A}$. The PDF of sampling $\mathbf{m}^\prime$ on the sphere through the crescent disk projection can be obtained by a change of domain from unit disk:

$$
\begin{equation}
p_c(\mathbf{m^{\prime}}) = p_s(\mathbf{m^{\prime}}) \frac{\partial(t_1, t_2)}{\partial(t_1^\prime, t_2^\prime)} = \frac{1}{s}D_s(\mathbf{m}^\prime, \mathbf{v}^\prime)(\mathbf{m}^\prime\cdot\mathbf{v}^\prime) \label{eq:m-prime-pdf}
\end{equation}
$$

This is the PDF of $\mathbf{m}^\prime$ obtained by following step 2 through 5. Below, we will demonstrate that when we transform $\mathbf{m}^\prime$ sampled in this manner to its ellipsoid space counterpart $\mathbf{m}$, the PDF of $\mathbf{m}$, denoted by $p(\mathbf{m})$, matches the expression given by Heitz[^Heitz2018] in Eq.3.

We can make the following observations:

$$
\begin{eqnarray}
D(\mathbf{m}, \mathbf{v})
&=& \frac{\left\lvert \mathbf{A} \right\rvert}{\pi \left\lVert \mathbf{Av} \right\rVert \left( \left( \mathbf{A}^{-\mathsf{T}}\mathbf{m} \right)^{-\mathsf{T}} \left( \mathbf{A}^{-\mathsf{T}}\mathbf{m} \right) \right)^2 }
= \frac{\left\lvert \mathbf{A} \right\rvert}{\left\lVert \mathbf{Av} \right\rVert} D_s(\mathbf{m}^\prime, \mathbf{v}^\prime) \label{eq:d-sphere-to-ellipsoid} \\

\mathbf{m^{\prime}}\cdot\mathbf{v^{\prime}}
&=& \frac{(\mathbf{A}^{-\mathsf{T}}\mathbf{m})^{\mathsf{T}} \mathbf{A}\mathbf{v}}{\left\lVert \mathbf{A}\mathbf{v} \right\rVert \left\lVert \mathbf{A}^{-\mathsf{T}}\mathbf{m} \right\rVert}
= \frac{\mathbf{m}\cdot\mathbf{v}}{\left\lVert \mathbf{A}\mathbf{v} \right\rVert \left\lVert \mathbf{A}^{-\mathsf{T}}\mathbf{m} \right\rVert} \label{eq:mv-sphere-to-ellipsoid}
\end{eqnarray}
$$

Substitute \eqref{eq:mv-sphere-to-ellipsoid} and \eqref{eq:d-sphere-to-ellipsoid} into \eqref{eq:m-prime-pdf}:

$$
\begin{equation}
p_c(\mathbf{m^{\prime}})= \frac{1}{s} \frac{D(\mathbf{m}, \mathbf{v}) (\mathbf{m}\cdot\mathbf{v})}{\left\lvert \mathbf{A} \right\rvert \left\lVert \mathbf{A}^{-\mathsf{T}}\mathbf{m} \right\rVert}
\end{equation}
$$

$p_c(\mathbf{m}^\prime)$ and $p(\mathbf{m})$ are related by

$$
\begin{equation}
p_c(\mathbf{m}^\prime) d\Omega^\prime = p(\mathbf{m}) d\Omega \Rightarrow p(\mathbf{m}) = p_c(\mathbf{m}^\prime) \left\lvert \frac{d\Omega^\prime}{d\Omega} \right\rvert
\end{equation}
$$

$\left\lvert \frac{d\Omega^\prime}{d\Omega} \right\rvert$ is the Jacobian determinant accounting for the transformation of the solid angle from ellipsoid space to sphere space. Note that it is the reciprocal of $\left\lvert \frac{d\mathbf{m}^\prime}{d\mathbf{m}} \right\rvert$ given by \eqref{eq:m-to-m-prime}:

$$
\begin{equation}
\left\lvert \frac{\partial\Omega^\prime}{\partial\Omega} \right\rvert = \left\lvert \mathbf{A} \right\rvert \left\lVert \mathbf{A}^{-\mathsf{T}}\mathbf{m} \right\rVert \label{eq:solod-angle-jacobian}
\end{equation}
$$

We now arrive at:

$$
\begin{equation}
p(\mathbf{m}) = \frac{1}{s}D(\mathbf{m},\mathbf{v})(\mathbf{m}\cdot\mathbf{v})
\end{equation}
$$

According to part 1 the projected area on $\mathbf{n}$ direction and on $\mathbf{v}$ direction are related by a ratio of projected area in $\mathbf{n}$ and $\mathbf{v}$:

$$
\begin{equation}
D(\mathbf{m}, \mathbf{v}) = \frac{A^\perp(\mathbf{n})}{A^\perp(\mathbf{v})} D(\mathbf{m}) = \frac{\left\lVert \mathbf{A}\mathbf{n} \right\rVert}{\left\lVert \mathbf{A}\mathbf{v} \right\rVert}D(\mathbf{m})
\end{equation}
$$

We have:

$$
\begin{equation}
p(\mathbf{m}) = \frac{1}{s}\frac{\left\lVert \mathbf{A}\mathbf{n} \right\rVert}{\left\lVert \mathbf{A}\mathbf{v} \right\rVert}D(\mathbf{m})(\mathbf{m}\cdot\mathbf{v}) \label{eq:pm}
\end{equation}
$$

Susbstitute $s$ with \eqref{eq:s-mapping}:

$$
\begin{equation}
\frac{1}{s}\frac{\left\lVert \mathbf{A}\mathbf{n} \right\rVert}{\left\lVert \mathbf{A}\mathbf{v} \right\rVert} = \frac{2\left\lVert \mathbf{A}\mathbf{n} \right\rVert \left\lVert \mathbf{A}^{-\mathsf{T}}\mathbf{n} \right\rVert}{\left\lVert \mathbf{A}\mathbf{v} \right\rVert \left\lVert \mathbf{A}^{-\mathsf{T}}\mathbf{n} \right\rVert + (\mathbf{v}\cdot\mathbf{n})} \label{eq:g1-pdf}
\end{equation}
$$

Considering $\mathbf{n}$ points in the $+z$ direction, both $\left\lVert \mathbf{A}\mathbf{n} \right\rVert$ and $\left\lVert \mathbf{A}^{-\mathsf{T}}\mathbf{n} \right\rVert$ take the value of 1. A comparison between \eqref{eq:g1} and \eqref{eq:g1-pdf} reveals:

$$
\begin{equation}
\frac{1}{s}\frac{\left\lVert \mathbf{A}\mathbf{n} \right\rVert}{\left\lVert \mathbf{A}\mathbf{v} \right\rVert} = \frac{G_1(\mathbf{v})}{\mathbf{v}\cdot\mathbf{n}} \label{eq:g1-equality}
\end{equation}
$$

Substitute \eqref{eq:g1-equality} into \eqref{eq:pm}. We obtain the expression of $p(\mathbf{m})$ that exactly matches the PDF given by Heitz[^Heitz2018]: 

$$
\begin{equation}
p(\mathbf{m}) = \frac{G_1(\mathbf{v})D(\mathbf{m})(\mathbf{m}\cdot\mathbf{v})}{\mathbf{v}\cdot\mathbf{n}} \label{eq:vndf-pdf}
\end{equation}
$$ 

# Final Considerations
This blog post is largely inspired by Walter's 2016 ellipsoid NDF article[^Walter2016]. It addresses the same core topic while further generalizing the ellipsoid through the introduction of rotational components to $\mathbf{A}$. The article provides a detailed derivation regarding the ellipsoid NDF and masking terms. However, its treatment of the PDF of VNDF sampling, as demonstrated in the fourth chapter, contains logical discontinuities. This blog post attempts to explain this PDF directly from the sampling routine, thereby demonstrating that the PDF indeed takes the form of \eqref{eq:vndf-pdf}.

# References
[^GraphicsGuysNotes]: [Sampling Anisotropic Microfacet BRDF](https://agraphicsguynotes.com/posts/sample_anisotropic_microfacet_brdf/)
[^Walter2007]: [Bruce Walter. Microfacet Models for Refraction through Rough Surfaces](https://www.graphics.cornell.edu/~bjw/microfacetbsdf.pdf)
[^Heitz2014]: [Eric Heitz. Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs](https://jcgt.org/published/0003/02/03/paper.pdf)
[^Heitz2018]: [Eric Heitz. Sampling the GGX Distribution of Visible Normals](https://jcgt.org/published/0007/04/01/paper.pdf)
[^Walter2016]: [Bruce Walter. The Ellipsoid Normal Distribution Function](https://www.cs.cornell.edu/Projects/metalappearance/SupplementalEllipsoidNDF.pdf)