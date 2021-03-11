Algorithm Comparison
========================================================
author: Matt Weldy
date: 03/10/2021
autosize: true



Goal for today 1
========================================================

Present the fit for two commonly used models in ecology using: 

- JAGS
- NIMBLE
- Stan
- greta

MCMC Samplers
========================================================

Quick refresher:
- Gibbs 
- No-U-Turn sampling
- Hamiltonian Monte Carlo

https://chi-feng.github.io/mcmc-demo/

Key points:
- Choice of sampler is really important for ecologists
- Latent Discreet Parameters
- Are there any benefits to marginalization?

JAGS
========================================================

Just Another Gibbs Sampler
Martyn Plummer

- BUGS based
- Slice sampling (Neal 2003) Gibbs with another name
- Lots... of example code out in the wild (never trust it blindly)
- Samples latent discreet variables
- Limited flow control use of the step function
- Support for functions in manual is limited


NIMBLE
========================================================

Perry de Valpine, UC Berkley
Daniel Turek gets quite a bit of credit

- Perry is an ecologist
- BUGS based
- NimbleEcology
- Daniel spent time optimizing spatial models
- Modular
- Extensible
- step function
- DAG support: pump$plotGraph()

Stan
========================================================

Stanislaw Ulam, a mathematician who was one of the developers of the MonteCarlo method in the 1940s
34 primary contributors currently listed
Andrew Gelman, Bob Carpenter, Matt Hoffman, and Daniel Lee...

- No-U-Turn Sampling
- Hamiltonian Monte Carlo
- C++
- Well documented
- Supported by a community
- Stan ecology website
- UBMS Ken Kellner
- If/else control flow

Stan
========================================================

````
functions {
  // ... function declarations and definitions ...
}
data {
  // ... declarations ...
}
transformed data {
   // ... declarations ... statements ...
}
parameters {
   // ... declarations ...
}
transformed parameters {
   // ... declarations ... statements ...
}
model {
   // ... declarations ... statements ...
}
generated quantities {
   // ... declarations ... statements ...
}
````

greta
========================================================

Grete Hermann wasn't a statistician. She wrote the first algorithms for computer algebra; in the 1920s, well before the first electronic computer was built.

Nick Golding, active developer

- Nick is an ecologist
- Very unfamiliar for those used to JAGS
- One person show (Although I don't think Nick sleeps)
- Seems to be well connected with Hadley Wickam and other key R/RStudio developers R6
- Make Bayes analyses more contained within R
- DAG support
- Hamiltonian Monte Carlo
- Wants to add support for discreet parameters
- R control flow

greta
========================================================


```r
library(greta)
(a <- variable(dim = c(3, 3)))
```

```
greta array (variable)

     [,1] [,2] [,3]
[1,]  ?    ?    ?  
[2,]  ?    ?    ?  
[3,]  ?    ?    ?  
attr(,"class")
[1] matrix array 
```


```r
(z <- normal(0, 1, truncation = c(-1, 1)))
```

```
greta array (variable following a normal distribution)

     [,1]
[1,]  ?  
attr(,"class")
[1] matrix array 
```

greta
========================================================

```r
c <- variable(dim = c(40, 1))
dim(c)
```

```
[1] 40  1
```

```r
head(c)
```

```
greta array (operation)

     [,1]
[1,]  ?  
[2,]  ?  
[3,]  ?  
[4,]  ?  
[5,]  ?  
[6,]  ?  
attr(,"class")
[1] matrix array 
```

Huggins
========================================================

Derivative from the full likelihood abundance estimator. The population size N is conditioned out of the likelihood. 

The capture history $y_{i,t}$ is used to estimate the capture probability of individual $i$ as a Bernoulli trial,
$$y_{i,t} \sim Bernoulli(p_{i,t}) $$
$$\mathcal{L}(p| y) = \prod_{i=1}^n \prod_{t=1}^t p_{i,t}^{z_{i,t}}(1-p_{i,t})^{1-z_{i,t}}$$

Abundance $\hat{N}$ is derived conditional on the count of known individuals ($C$), sometimes referred to the minimum number of known alive ($MNKA$).
$$\hat{N} = \frac{C}{1-\prod^{t}(1-p_t)}$$

CJS
========================================================

State-space parametrization 
$$\mathcal{L}(\phi, p, z| y) = f(z_1|\phi) \prod_{t=2}^T f(z_t| z_{t-1}, \phi) \prod_{t=1}^T f(y_t| z_t, p)$$
State Process
$$z_{i_f} = 1$$
$$z_{i,t+1}|z_{i,t} \sim Bernoulli(z_{i,t} \phi_{i,t}) $$
Observation Process
$$y_{i,t}|z_{i,t} \sim Bernoulli(z_{i,t},p_{i,t})$$ 

