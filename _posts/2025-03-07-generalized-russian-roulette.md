---
layout: post
title:  "Generalized Russian Roulette for Monte Carlo Estimator"
date:   2025-03-07 16:05:00 +0800
use_math: true
---

When tracing backwards along a light path, we may encounter situations where multiple interactions occur between the light and an object at a point in space. One example is a dielectric surface. The interactions are reflection and refraction. Exitant radiance going out in direction $\mathbf{o}$ is computed as:

$$
\begin{eqnarray}
L(\mathbf{\omega_o}) &=& \int_{\Omega_i} \left(f_r + f_t \right) \lvert\cos\theta_i\rvert L(\omega_i) d\omega_i \\
&=& \int_{\Omega_i} \left(FR + (1-F)T \right) \lvert\cos\theta_i\rvert L(\omega_i) d\omega_i
\end{eqnarray}
$$

$F$ is the Fresnel reflectance. $R$ is the product of all the other terms in $f_r$ except for $F$. Separate transmittance $1-F$ from $f_t$, product of remaining terms is denoted as $T$.


A simplified path tracing solution may look like

```plaintext
function scatter(wo)
    wm <- sample_visible_normal(wo)
    p <- visible_normal_PDF(wo, wm)
    F <- fresnel(wo, wm)
    q <- sample_random()
    if q < F then
        wi <- reflect(wo, wm)
        brdf <- compute_BRDF(wo, wm, wi, F)
        pr <- reflection_jacobian(wi, wm) * p
        return (brdf, pr, wi)
    else 
        if total_intenal_reflection(wo, wm) then
            return invalid
        wi <- refract(wo, wm)
        btdf <- compute_BTDF(wo, wm, wi, 1-F)
        pt <- refraction_jacobian(wi, wm) * p
        return (btdf, pt, wi)
end function
```

Whether to sample reflection or refraction is determined by sampling a random viariable and checking if it is larger than reflectance $R$. Whichever one of $R$ or $1-R$ makes a greater contribution to outgoing radiance in $\omega_o$ direction is more likely to get selected and have $\omega_i$ sampled against.

Another example is when light interacts with a medium. The interactions are absorption and scattering. The recursive VRE when null-scattering is accounted for can be found in section 2.3.3 of the Pixar course presented in SIGGRAPH 2017[^a]. In it, algorithm 2 shows a piece of pseudocode of the delta tracking process. For simplicity's sake, We will assume the sum of absoption coefficient $\sigma_a$ and scattering coefficient $\sigma_s$ is a constant at every point in a medium. The null-scattering coefficient $\sigma_n$ is always $0$ so that light interactions with null-particles never gets sampled.  

We can see a pattern in the above two examples. Each one of them implements an estimator that introduces a random $Y$ to stochastically select between different types of light interactions and then proceed to perform Monte Carlo sampling on the selected interaction. $Y$ follows Bernoulli distribution:

$$
P(Y = y) =
\begin{cases} 
q, & \text{if } y = 1 \\
1 - q, & \text{if } y = 0
\end{cases}
$$

In the case of sampling surface BSDF, the two different types of interactions are reflection and refraction. Selection probability $q$ takes the value of the Fresnel term. In the case of volume tracking, the interactions are absorption and scattering. $q$ is defined as $\sigma_a/\sigma_t=\sigma_a/(\sigma_a + \sigma_s)$.

# Estimator formulation

If we are trying to estimate the integral of a linear interpolation of functions $f(x)$ and $g(x)$:  

$$
\begin{equation}
\int_{-\infty}^{\infty} qf(x) + (1-q)g(x) dx
\end{equation}
$$

The estimator is formulated as:

$$
\begin{equation}
F_n=\frac{1}{n}\sum_{i=1}^{n}\left( \frac{Y_if(X_i)}{p(X_i)} + \frac{(1-Y_i)g(Z_i)}{p(Z_i)} \right) \label{eq:estimator-general}
\end{equation}
$$

In the case of surface BSDF, random variables in \eqref{eq:estimator-general} corresponds to:

- $X_i$ -- incoming direction of reflection 
- $Y_i$ -- selection variable
- $Z_i$ -- incoming direction of refraction
- $f(X_i)$ -- BRDF
- $g(Z_i)$ -- BTDF
- $p(X_i)$ -- probability of $X_i$ being sampled. Visible normal probability transformed into incoming domain of reflection.
- $p(Z_i)$ -- probability of $Z_i$ being sampled. Visible normal probability transformed into incoming domain of refraction.

In the case of delta tracking, only one random variable that follows exponential distribution is sampled every round of Monte Carlo evaluation. Thus $Z_i$ and $X_i$ in \eqref{eq:estimator-general} can be combined:

$$
\begin{equation}
F_{delta} = \frac{1}{n}\sum_{i=1}^{n}\left( \frac{Y_if(X_i) + (1-Y_i)g(X_i)}{p(X_i)} \right) 
\end{equation}
$$

- $X_i$ -- incoming direction 
- $Y_i$ -- selection variable
- $f(X_i)$ -- emittance
- $g(Z_i)$ -- in-scattering
- $p(X_i)$ -- probability of $X_i$ being sampled. Follows an exponential distribution that corresponds to volumn transmittance

# Expectation


Is this estimator unbiased? To answer this question we need to take a look at its expectation. Expectation of $F_n$ is:

$$
\begin{eqnarray}
E(F_n) &=& E\left(\frac{1}{n}\sum_{i=1}^{n}\left( \frac{Y_if(X_i)}{p(X_i)} + \frac{(1-Y_i)g(Z_i)}{p(Z_i)} \right)\right) \label{eq:expectation-1} \\
&=& \frac{1}{n}\sum_{i=1}^{n}E\left( \frac{Y_if(X_i)}{p(X_i)} + \frac{(1-Y_i)g(Z_i)}{p(Z_i)} \right) \label{eq:expectation-2} \\
&=& \frac{1}{n}\sum_{i=1}^{n} \left( E(Y_i) E\left(\frac{f(X_i)}{p(X_i)}\right) + E(1-Y_i) E\left(\frac{g(Z_i)}{p(Z_i)}\right) \right) \label{eq:expectation-3} \\
&=& \frac{1}{n}\sum_{i=1}^{n} \left( q \int_{-\infty}^{\infty}f(x)dx + (1-q) \int_{-\infty}^{\infty}g(x)dx\right) \label{eq:expectation-4} \\
&=& \int_{-\infty}^{\infty} qf(x) + (1-q)g(x) dx \label{eq:expectation-final}
\end{eqnarray}
$$

The condition that allows the expansion from \eqref{eq:expectation-2} to \eqref{eq:expectation-3} is that $X$ and $Y$ are independent. The expectation is exactly the integral we are trying to estimate, thus the estimator is unbiased. This estimator is not limited to estimating the integral of two linearly combined functions. The way this Pixar course[^a] carries out Delta tracking when null-scattering coefficient $\sigma_n$ is non-zero involves an interpolation between three functions. By extending Bernoulli distribution to three outcomes each assigned a probability, with similar reasoning we can formulate an estimator similar to \eqref{eq:expectation-2}. It is also unbiased.  

# Another Specialized Case: Russian Roulette

Russian Roulette is a technique that aims to improve the efficiency of Monte Carlo estimator. It is characterized by
1. stochastically skipping evaluation of samples with a probability $q$
2. replacing skipped samples with a constant $c$
3. subtracting $qc$ from non-skipped samples $f(X_i)$ 
4. weighting non-skipped samples by a factor of $1/(1-q)$ 

Following this description, the estimation given by Russian Roulette is:

$$
\begin{equation}
F_{roulette}=\frac{1}{n}\sum_{i=1}^{n}\left(\frac{\frac{(1-Y_i)(f(X_i)-qc)}{1-q}}{p(X_i)} +Y_ic\right)
\end{equation}
$$

We can tell $F_{roulette}$ is also a specialized version of $F_n$ if we make the following substitution in $\eqref{eq:estimator-general}$:

$$
\begin{eqnarray}
f(X_i) &\leftarrow& \frac{f(X_i)-qc}{1-q} \\
g(Z_i) &\leftarrow& c \\
p(Z_i) &\leftarrow& 1 \\
Y_i &\leftarrow& 1 - Y_i \\
1 - Y_i &\leftarrow& Y_i \\
\end{eqnarray}
$$

The above substitutions when translated to expectation are:

$$
\begin{eqnarray}
f(x) &\leftarrow& \frac{f(x)-qc}{1-q} \\
g(x) &\leftarrow& c \\
q &\leftarrow& 1 - q \\
1 - q &\leftarrow& q \\
\end{eqnarray}
$$

When applied to \eqref{eq:expectation-final} gives:

$$
\begin{eqnarray}
E(F_{roulette}) &=& \int_{-\infty}^{\infty}  \left((1-q)\frac{f(x)-qc}{1-q} + qc \right) dx \\
&=& \int_{-\infty}^{\infty} f(x) dx
\end{eqnarray}
$$

# References
[^a]: [Julian Fong. Production Volume Rendering SIGGRAPH 2017 Course](https://graphics.pixar.com/library/ProductionVolumeRendering/paper.pdf)
