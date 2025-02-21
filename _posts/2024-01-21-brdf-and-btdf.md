---
layout: post
title:  "Micro-BRDF and BTDF of Ideal Microfacets"
date:   2025-01-22 10:32:50 +0800
use_math: true
---
The article by Eric Heitz[^a] provides a rigorous construction process of the BRDF. By the end of section 3.2, it arrives at equation (26), which is the macro-BRDF expressed as an integral of micro-BRDF over hemisphere.

Then immediately in section 3.3, it dumps the micro-BRDF for mirror-like microfacets with minimum amount of explanation, merely mentioning the Jacobian of the domain transformation from halfway vector $\omega_h$'s domain to its corresponding incident vector $\omega_i$'s domain.

I followed the citation link and looked at the famous GGX article by Bruce Walter[^c]. In it, an elegant geometric derivation of the Jacobian is illustrated in Figure 6. But still, it glosses over the reason why terms like $\left\lvert \mathbf{o}\cdot\mathbf{m} \right\rvert$ and $\mathrm{\rho}$ exist in equation (9) through (11). 

Then I found some clues in the slides of a presentation given by Earl Hammon Jr.[^b] Page 30 gives an idea as how the specular micro-BRDF is normalized. I went down this path and came up with an explanation of micro-BRDF and micro-BTDF that I find sufficiently coherent.

# Revisiting the Integral Form of Macro-BxDF
This section aims to provide some intuition about the integral form of macro-BxDF. Without loss of generality, we will use the term 'BxDF' to cover both BRDF and BTDF. The term 'scatter' covers both reflection and refraction. Equations and notations in this section are mostly taken from section 3.1 and 3.2 of Eric Heitz's article[^a].

The total amount of radiance going out at direction $\omega_o$ is given by:

$$
\begin{equation}
L(\omega_o)=\int_{\Omega}L(\omega_o, \omega_m)D_{\omega_o}(\omega_m)d\omega_m \label{eq:macro-out-radiance}
\end{equation}
$$

+ $\Omega$ -- Hemispherical integral domain of $\omega_m$
+ $L(\omega_o, \omega_m)$ -- exitant radiance in $\omega_o$ scattered by all microfacets facing direction $\omega_m$
+ $D_{\omega_o}(\omega_m)$ -- Distribution of visible normals, viewed from direction $\omega_o$.

It says that the total exitant radiance at certain viewing direction $\omega_o$ equals the sum of exitant radiance in that direction emitted by microfacets facing direction $\omega_m$, weighted by their corresponding 'visibility' in terms of a visible normal distribution $D_{\omega_o}(\omega_m)$ . Notice how this equation says nothing about incident radiance.

We now turn to look at how microfacets facing direction $\omega_m$ interacts with incident irradiance. Micro-BxDF is defined as[^a]:

$$
\begin{equation}
\rho(\omega_o, \omega_i, \omega_m)=\frac{dL(\omega_o, \omega_m)}{(\omega_i \cdot \omega_m) L(\omega_i)d\omega_i} \label{eq:micro-brdf-def}
\end{equation}
$$

It is defined by differentiating exitant **radiance** at $\omega_o$ over differential incident **irradiance** at $\omega_i$. We replace the clamp dot notation of $\left\langle\omega_i, \omega_m\right\rangle$ with a simple dot product by assuming $\omega_i$ and $\omega_m$ are always in the same hemisphere. The exitant radiance of this specific microfacet can then be expressed as:

$$
\begin{equation}
L(\omega_o, \omega_m)=\int_{\Omega_i}\rho(\omega_o, \omega_i, \omega_m)(\omega_i \cdot \omega_m) L(\omega_i)d\omega_i \label{eq:micro-out-radiance}
\end{equation}
$$

It shows that microfacet facing certain direction $\omega_m$ collects irrandiance coming from $\omega_i$ through a patch $d\omega_i$, modulate it with micro-BxDF $\rho(\omega_o, \omega_i, \omega_m)$, integrate over the entirety of incident domain, then scatter radiance towards one specific exitant direction $\omega_o$. Account for all microfacets by substituting \eqref{eq:micro-out-radiance} into \eqref{eq:macro-out-radiance}:

$$
\begin{equation}
L(\omega_o)=\int_{\Omega}\int_{\Omega_i}\rho(\omega_o, \omega_i, \omega_m)(\omega_i \cdot \omega_m) L(\omega_i) d\omega_i D_{\omega_o}(\omega_m)d\omega_m \label{eq:maco-out-radiance-expanded}
\end{equation}
$$

+ $\Omega_i$ -- Hemispherical integral domain of $\omega_i$


We now know,
1. (integral over $\Omega_i$) how microfacets facing $\omega_m$ collects all incident irradiance and scatter it towards $\omega_o$, and, 
2. (integral over $\Omega$) how to take sum of exitant radiance emitted by microfacets over all visible $\omega_m$ s. 

Next, define macro-BxDF using an approach similar to that of micro-BxDF:

$$
\begin{equation}
\rho(\omega_o, \omega_i)=\frac{dL(\omega_o)}{\cos{\theta_i} L(\omega_i)d\omega_i} \label{eq:macro-brdf-def}
\end{equation}
$$

+ $\theta_i$ -- the angle of incidence between $\omega_i$ and macrosurface normal.

Notice how \eqref{eq:macro-brdf-def} relates exitant radiance to incident irradiance and has nothing to do with microfacet normal $\omega_m$ because of its macroscopic nature.

Substitute $\frac{dL(\omega_o)}{d\omega_i}$ by taking derivative of \eqref{eq:maco-out-radiance-expanded}. \eqref{eq:macro-brdf-def} expands to:

$$
\begin{eqnarray}
\rho(\omega_o, \omega_i)
&=& \frac{1}{\cos{\theta_i}} \int_{\Omega} \rho(\omega_o, \omega_i, \omega_m)(\omega_i \cdot \omega_m) D_{\omega_o}(\omega_m) d\omega_m \\
&=& \frac{1}{\cos{\theta_i}\cos{\theta_o}} \int_{\Omega} \rho(\omega_o, \omega_i, \omega_m)(\omega_i \cdot \omega_m) G_2(\omega_o, \omega_i, \omega_m)\left\lvert\omega_o \cdot \omega_m\right\rvert D(\omega_m) d\omega_m
\label{eq:macro-brdf-integral}
\end{eqnarray}
$$

The rationale behind the expansion of $D_{\omega_o}(\omega_m)$ and the masking-shadowing term can be found in section 3.2 of Eric Heitz's article[^a]. Taking absolute value of $\omega_o\cdot\omega_m$ term accommodates both reflection and refraction cases. $G_2$ is also formulated in a way that makes it inherently invariant to the sign of $\omega_o\cdot\omega_m$.

We eventually come to the integral form of macro-BxDF with $d\omega_m$ as the only differential element. This equation tells us how micro-BxDF is related to macro-BxDF by integrating the former over visible normals. At this point, the final puzzle piece needed to instantiate an explicit macro-BxDF expression is a $\rho(\omega_o, \omega_i, \omega_m)$ that could allow \eqref{eq:macro-brdf-integral} to be solved analytically.

# Micro-BxDF of Ideal Microfacets
We first introduce notations:
+ $\mathbf{o}$ -- one specific exitant direction
+ $\mathbf{i}$ -- one specific incident direction defined in the $\Omega_i$ domain
+ $\mathbf{h}$ -- half-direction vector defined in the $\Omega$ domain. A generalized term that covers both cases illustrated in [Figure 1](#fig1).


<div style="text-align: center; margin: 40px 0;">
  <figure id="fig1">
    <img src="/images/reflection.svg" alt="refraction Jacobian geometric derivation" style="width: 80%; height: auto; display: block; margin: 0 auto;">
    <figcaption>Figure 1: Configurations of $\mathbf{i}$, $\mathbf{o}$ and $\mathbf{h}$ in reflection and refraction. </figcaption>
  </figure>
</div>

We are going to focus on the case of ideal microfacets. An ideal microfacet: 
1. Is Perfectly smooth. Given $\mathbf{o}$ and microfacet normal $\omega_m$, the only direction the incident radiance should come from is $\mathbf{i}$, where $\omega_m$ coincides with $\mathbf{h}$.
2. Turns the entirety of the energy carried by incident radiance into either reflected or refracted energy in the form of exitant radiance. Fresnel term will come into play later.

The above two constraints can be formulated as:

$$
\begin{equation}
\rho(\omega_o, \omega_i, \omega_m) = k\delta_\mathbf{i}(\omega_i) \label{eq:micro-brdf-incident-space}
\end{equation}
$$

$$
\begin{equation}
\int_{\Omega_i}\rho(\omega_o, \omega_i, \omega_m)(\omega_i \cdot \omega_m)L(\omega_i)d\omega_i=L(\omega_o) \label{eq:energy-conserve}
\end{equation}
$$

$\delta_\mathbf{i}(\omega_i)$ is a Dirac delta defined in the $\Omega_i$ domain that spikes at $\mathbf{i}$. Introducing a delta function to the integrand allows integral like \eqref{eq:macro-brdf-integral} to be evaluated as its integrand at a fixed configuration of $\mathbf{i}$, $\mathbf{o}$ and $\mathbf{h}$.  \eqref{eq:energy-conserve} is a generalized version of the 'normalized BRDF' equation on page 30 of Earl Hammon Jr.'s lecture[^b].


Substitute \eqref{eq:micro-brdf-incident-space} into \eqref{eq:energy-conserve}, we can solve for $k$. But there's a catch. Equation \eqref{eq:macro-brdf-integral} integrates over microfacts' normal domain $\Omega$ while \eqref{eq:micro-brdf-incident-space} is defined in the incident domain $\Omega_i$. We resolve this mismatch of domains by rewriting $\rho(\omega_o, \omega_i, \omega_m)$ in the normal domain $\Omega$.

$$
\begin{equation}
\rho(\omega_o, \omega_i, \omega_m) = k\delta_{\mathbf{h}}(\omega_m) 
\label{eq:micro-brdf-normal-space}
\end{equation}
$$

An additional change of domain operation should also be applied to \eqref{eq:energy-conserve} to accomodate \eqref{eq:micro-brdf-normal-space}. It takes the form of a Jacobian $\left\lVert\frac{\partial \omega_i}{\partial \omega_m}\right\rVert$:

$$
\begin{equation}
\int_{\Omega}k\delta_{\mathbf{h}}(\omega_m)(\omega_i \cdot \omega_m)L(\omega_i)\left\lVert \frac{\partial \omega_i}{\partial \omega_m} \right\rVert d\omega_m=L(\omega_o) \label{eq:norm-constraint-normal-space}
\end{equation}
$$

$\delta_{\mathbf{h}}$ evaluates the integral:

$$
\begin{equation}
k(\mathbf{i} \cdot \mathbf{h})L(\mathbf{i})\left\lVert \frac{\partial \omega_\mathbf{i}}{\partial \omega_\mathbf{h}} \right\rVert= L(\mathbf{o}) \label{eq:k-solve}
\end{equation}
$$

$$
\begin{equation}
k=\frac{L(\mathbf{o})}{L(\mathbf{i})}\frac{1}{\mathbf{i} \cdot \mathbf{h}}\left\lVert \frac{\partial \omega_\mathbf{h}}{\partial \omega_\mathbf{i}} \right\rVert \label{eq:k-definition}
\end{equation}
$$

Account for the Fresnel term:

$$
\begin{equation}
\rho(\omega_o, \omega_i, \omega_m) = F(\mathbf{o},\mathbf{h})\frac{L(\mathbf{o})}{L(\mathbf{i})}\frac{1}{\mathbf{i} \cdot \mathbf{h}}\left\lVert \frac{\partial \omega_\mathbf{h}}{\partial \omega_\mathbf{i}} \right\rVert\delta_{\mathbf{h}}(\omega_m) 
\label{eq:micro-brdf-full}
\end{equation}
$$

$F(\mathbf{o},\mathbf{h})$ is a blanket term that includes both reflectance and transmittance. It should be replaced with $F_r(\mathbf{o},\mathbf{h})$ in BRDF and $1-F_r(\mathbf{o},\mathbf{h})$ in BTDF.

Substitute \eqref{eq:micro-brdf-full} into \eqref{eq:macro-brdf-integral} and do the same $\delta_\mathbf{h}$ evaluation trick:

$$
\begin{equation}
\rho(\omega_o, \omega_i) 
= \frac{F(\mathbf{o},\mathbf{h})G_2(\mathbf{o},\mathbf{i},\mathbf{h})D(\mathbf{h})}{\cos{\theta_i}\cos{\theta_o}} \frac{L(\mathbf{o})}{L(\mathbf{i})}\left\lVert \frac{\partial \omega_\mathbf{h}}{\partial \omega_\mathbf{i}} \right\rVert \left\lvert \mathbf{o} \cdot \mathbf{h} \right\rvert
\label{eq:macro-bxdf-general}
\end{equation}
$$


At this end, the only two factors that need to finalize are $\frac{L(\mathbf{o})}{L(\mathbf{i})}$ and $\left\lVert \frac{\partial \omega_\mathbf{h}}{\partial \omega_\mathbf{i}} \right\rVert$.

# Microfacets as Ideal Reflectors

There is really nothing to it in this case. $L(\mathbf{o})$ will always be equivalent to $L(\mathbf{i})$. Figure 6 in the GGX paper proves $\left\lVert \frac{\partial \omega_\mathbf{h}}{\partial \omega_\mathbf{i}} \right\rVert=\frac{1}{4\left\lvert \mathbf{o}\cdot\mathbf{h} \right\rvert}$.[^c] Thus \eqref{eq:macro-bxdf-general} is instantiated as:

$$
\begin{equation}
\rho(\omega_o, \omega_i) 
= \frac{F_r(\mathbf{o},\mathbf{h})G_2(\mathbf{o},\mathbf{i},\mathbf{h})D(\mathbf{h})}{4\cos{\theta_i}\cos{\theta_o}}
\label{eq:brdf}
\end{equation}
$$

We arrive at the most commonly known BRDF equation.

# Microfacets as Ideal Refractors

The relationship between incident and exitant radiance is:

$$
\begin{equation}
\frac{L(\mathbf{o})}{L(\mathbf{i})}=\frac{\eta_o^2}{\eta_i^2} \label{eq:radiance-compression}
\end{equation}
$$

$\eta$ denotes the index of refraction. An explanation of this equation based on energy conservation can be found in 9.5.2 of the PBRT book fourth edition.[^d] Intuitively it can be explained as incident irradiance getting compressed to a smaller solid angle at the interface when light propagates from a optically rarer medium to a denser one, resulting in a stronger exitant radiance. 


As for the Jacobian, there is also a geometric derivation, shown in Figure 7 of the GGX paper[^c]. The Jacobian is formulated as $\left\lVert\frac{\partial\omega_{\mathbf{h}}}{\partial\omega_{\mathbf{o}}}\right\rVert$.  While in our case, we need $\left\lVert\frac{\partial\omega_{\mathbf{h}}}{\partial\omega_{\mathbf{i}}}\right\rVert$, the rate of change of $\omega_\mathbf{h}$ with respect to $\omega_\mathbf{i}$. [Figure 2](#fig2) illustrates how $\left\lVert\frac{\partial\omega_{\mathbf{h}}}{\partial\omega_{\mathbf{i}}}\right\rVert$ is derived. The idea is similar to that of $\left\lVert\frac{\partial\omega_{\mathbf{h}}}{\partial\omega_{\mathbf{o}}}\right\rVert$ presented in the GGX paper.


<div style="text-align: center; margin: 40px 0;">
  <figure id="fig2">
    <img src="/images/transmission_jacobian.svg" alt="refraction Jacobian geometric derivation" style="width: 40%; height: auto; display: block; margin: 0 auto;">
    <figcaption>Figure 2: Snell's Law implies $\eta_i\mathbf{i}+\eta_o\mathbf{o}$ and $\mathbf{h}$ are colinear. Thus $c\mathbf{h}=-(\eta_i\mathbf{i}+\eta_o\mathbf{o})$. The negative sign is added to ensure $c$ is always positive. $c=c\mathbf{h}\cdot\mathbf{h}=\eta_i(\mathbf{i}\cdot\mathbf{h})+\eta_o(\mathbf{o}\cdot\mathbf{h})$. </figcaption>
  </figure>
</div>

The rate of change of $\omega_\mathbf{h}$ with respect to $\omega_\mathbf{i}$ is the ratio of the infinitesimal surface area on the $\mathbf{h}$ unit sphere to that on the $\mathbf{i}$ unit sphere: 


$$
\begin{equation}
\left\lVert\frac{\partial\omega_{\mathbf{h}}}{\partial\omega_{\mathbf{i}}}\right\rVert
=\left\lvert\frac{d\omega_\mathbf{h}}{d\omega_\mathbf{i}}\right\rvert
=\frac{\eta_i^2\left\lvert\mathbf{i}\cdot\mathbf{h}\right\rvert}{c^2}
=\frac{\eta_i^2\left\lvert\mathbf{i}\cdot\mathbf{h}\right\rvert}{(\eta_i(\mathbf{i}\cdot\mathbf{h})+\eta_o(\mathbf{o}\cdot\mathbf{h}))^2}
\label{eq:jacobian-h-i}
\end{equation}
$$

Substitiute \eqref{eq:radiance-compression} and \eqref{eq:jacobian-h-i} into \eqref{eq:macro-bxdf-general}:

$$
\begin{equation}
\rho(\omega_o, \omega_i) 
= \frac{(1-F_r(\mathbf{o},\mathbf{h}))G_2(\mathbf{o},\mathbf{i},\mathbf{h})D(\mathbf{h})}{\cos{\theta_i}\cos{\theta_o}} \frac{\eta_o^2\left\lvert\mathbf{i}\cdot\mathbf{h}\right\rvert \left\lvert \mathbf{o} \cdot \mathbf{h} \right\rvert}{(\eta_i(\mathbf{i}\cdot\mathbf{h})+\eta_o(\mathbf{o}\cdot\mathbf{h}))^2}
\end{equation}
$$

We have the rough dielectric BTDF.


# References
[^a]: [Eric Heitz. Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs](https://jcgt.org/published/0003/02/03/paper.pdf)
[^b]: [Earl Hammon Jr. PBR Diffuse Lighting for GGX+Smith Microsurfaces](https://ubm-twvideo01.s3.amazonaws.com/o1/vault/gdc2017/Presentations/Hammon_Earl_PBR_Diffuse_Lighting.pdf)
[^c]: [Bruce Walter. Microfacet Models for Refraction through Rough Surfaces](https://www.graphics.cornell.edu/~bjw/microfacetbsdf.pdf)
[^d]: [Physically Based Rendering: From Theory to Implementation. 9.5.2 Non-Symmetric Scattering and Refraction](https://pbr-book.org/4ed/Reflection_Models/Dielectric_BSDF#Non-SymmetricScatteringandRefraction)