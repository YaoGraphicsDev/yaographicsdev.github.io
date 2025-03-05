---
layout: post
title:  "Solving the Radiative Transfer Equation"
date:   2025-02-28 16:05:00 +0800
use_math: true
---

One of the fundamental equations in the field of volumetric rendering is the the Radiative Transfer Equation(RTE)[^a].

$$
\begin{equation}
\left(\omega\cdot \nabla \right)L(\mathbf{x}, \omega) = -\sigma_t(\mathbf{x})L(\mathbf{x}, \omega) + \sigma_a(\mathbf{x})L_e(\mathbf{x}, \omega) + \sigma_s(\mathbf{x})L_s(\mathbf{x}, \omega) \label{eq:RTE}
\end{equation}
$$

We focus on $x$ direction only by substituting $\omega = 1i + 0j + 0k$ into \eqref{eq:RTE} and drop the three-dimensional notation of position $\mathbf{x}$

$$
\begin{equation}
L'(x) = -\sigma_t(x)L(x) + \sigma_a(x)L_e(x) + \sigma_s(x)L_s(x) \label{eq:RTE-x}
\end{equation}
$$

We can further define a source term as contribution to radiance from emission and in-scattering at position $x$:

$$
\begin{equation}
S(x) = \sigma_a(x)L_e(x) + \sigma_s(x)L_s(x)
\end{equation}
$$

The Volume Rendering Equation(VRE) is given immediately afterwards

$$
\begin{equation}
L(x) = \int_0^x e^{-\int_u^x\sigma_t(v)dv}S(u)du + e^{-\int^x_0{\sigma_t(v)}dv} L(0) \label{eq:VRE}
\end{equation}
$$

This integral makes perfect intuitive sense. The first term indicates that, along a light path, the radiance contributed by a source $S(u)$ located at $u$ is attenuated exponentially by the integral of extinction coefficient $\sigma_t$ from position $u$ to the end of light path within medium, $x$. The second term accounts for the incoming radiance $L(0)$ at entry point and attenuates it along the entire light path from $0$ to $x$.


However, how VRE is derived from RTE is usually not covered by the articles I can find. The derivation process should not be difficult as RTE, when rearranged and written as: 

$$
\begin{equation}
L'(x)+\sigma_t(x)L(x)=S(x) \label{eq:RTE-linear-first-order}
\end{equation}
$$

falls into the well-studied category of linear first order ODE.

# Solution of Linear First Order ODE

A linear first order ODE takes the form of[^b],

$$
\begin{equation}
y'(x)+p(x)y(x)=f(x) \label{eq:LFO}
\end{equation}
$$

It can be solved with the integrating factor method[^c]. Multiply both sides of the equation with $I(x)$, we have:

$$
\begin{equation}
I(x)y'(x)+I(x)p(x)y(x)=I(x)f(x) \label{eq:LFO-multiply-factor}
\end{equation}
$$

According to the product rule, the left side of \eqref{eq:LFO-multiply-factor} can be seen as a derivative where 

$$
\begin{equation}
I'(x)=I(x)p(x) \label{eq:I-ODE}
\end{equation}
$$

Equation \eqref{eq:LFO-multiply-factor} then becomes:

$$
\begin{equation}
(I(x)y(x))'=I(x)f(x) \label{eq:LFO-multiply-factor-derivative}
\end{equation}
$$

Integrate both sides and divide by $I(x)$, we have:

$$
\begin{equation}
y(x)=\frac{1}{I(x)}\left(\int^x{I(u)f(u)}du+C\right) \label{eq:y-solution}
\end{equation}
$$

\eqref{eq:I-ODE} is also a linear first order ODE with a straightforward particular solution of:

$$
\begin{equation}
I(x)=e^{\int^x{p(v)}dv}
\end{equation}
$$

Substitute $I(x)$ into \eqref{eq:y-solution}:

$$
\begin{equation}
y(x)=e^{-\int^x{p(v)}dv}\left(\int^x{e^{\int^u{p(v)}dv}f(u)}du+C\right) \label{eq:y-solution-expanded}
\end{equation}
$$

Variable $u$ and $v$ are introduced to differentiate integral variable from function variable $x$.

For a function $f(x)$ that is analytic at $x=0$, the antiderivative defined by the definite integral $\int_0^xf(u)du$ is equivalent to an indefinite integral of $f(x)$. This can be shown by integrating the Taylor expansion of $f(x)$ term-by-term.

Assume all integrands in \eqref{eq:y-solution-expanded} are analytic. Instantiate the indefinite integrals as their definite counterpart:

$$
\begin{equation}
y(x)=e^{-\int^x_0{p(v)}dv}\left(\int^x_0{e^{\int^u_0{p(v)}dv}f(u)}du+C\right) \label{eq:y-solution-expanded-definite}
\end{equation}
$$

Apply boundary constraint $y(0)$, we can solve for constant $C$.

$$
\begin{equation}
y(x)=e^{-\int^x_0{p(v)}dv}\left(\int^x_0{e^{\int^u_0{p(v)}dv}f(u)}du+y(0)\right) \label{eq:y-solution-expanded-with-boundary}
\end{equation}
$$

# The Volume Rendering Equation

By comparing \eqref{eq:RTE-linear-first-order} and \eqref{eq:LFO}, we make the following substitutions to \eqref{eq:y-solution-expanded-with-boundary}:

$$
\begin{eqnarray}
y(x) &=& L(x) \\
p(x) &=& \sigma_t(x) \\
f(x) &=& S(x)
\label{eq:macro-brdf-integral}
\end{eqnarray}
$$

We have:

$$
\begin{equation}
L(x)=e^{-\int^x_0{\sigma_t(v)}dv}\left(\int^x_0{e^{\int^u_0{\sigma_t(v)}dv}S(u)}du+L(0)\right) \label{eq:RTE-1}
\end{equation}
$$

At this point, our solution does not exactly look like the VRE \eqref{eq:VRE}. We are going to take one final step of manipulation.

Integration variable $u$ goes from $0$ to $x$. We can expand the exponential terms:

$$
\begin{equation}
e^{\int^u_0{\sigma_t(v)}dv} = e^{\int^x_0{\sigma_t(v)}dv-\int^x_u{\sigma_t(v)}dv} = e^{\int^x_0{\sigma_t(v)}dv}e^{-\int^x_u{\sigma_t(v)}dv} \label{eq:e-expansion}
\end{equation}
$$

Equation \eqref{eq:e-expansion} is a function of $u$, so the term $e^{\int^x_0{\sigma_t(v)}dv}$ can be treated as a constant and be moved out of the $\int_0^xdu$ integral in \eqref{eq:RTE-1}. Substitute \eqref{eq:e-expansion} into \eqref{eq:RTE-1}:

$$
\begin{eqnarray}
L(x) &=& e^{-\int^x_0{\sigma_t(v)}dv}\left(\int^x_0{e^{\int^u_0{\sigma_t(v)}dv}S(u)}du+L(0)\right) \\
&=& e^{-\int^x_0{\sigma_t(v)}dv}\left(\int^x_0{e^{\int^x_0{\sigma_t(v)}dv}e^{-\int^x_u{\sigma_t(v)}dv}S(u)}du+L(0)\right) \\
&=& e^{-\int^x_0{\sigma_t(v)}dv}\left(e^{\int^x_0{\sigma_t(v)}dv}\int^x_0{e^{-\int^x_u{\sigma_t(v)}dv}S(u)}du+L(0)\right) \\
&=& \int^x_0{e^{-\int^x_u{\sigma_t(v)}dv}S(u)}du+e^{-\int^x_0{\sigma_t(v)}dv}L(0)
\end{eqnarray}
$$


# References
[^a]: [Julian Fong. Production Volume Rendering SIGGRAPH 2017 Course](https://graphics.pixar.com/library/ProductionVolumeRendering/paper.pdf)
[^b]: [Linear differential equation](https://en.wikipedia.org/wiki/Linear_differential_equation)
[^c]: [Solving Differential Equations with Integrating Factors](https://www.mathcentre.ac.uk/resources/uploaded/mathcentre-ode.pdf)
