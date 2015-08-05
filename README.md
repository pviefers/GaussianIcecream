# GaussianIcecream
Adds the fit of a Gaussian process with logistic link function to the discussion by [Markus Gesmann](http://www.magesblog.com/2015/08/generalised-linear-models-in-r.html).

The model that is implemented in addition to those in the original blog post is a Gaussian process with a logistic link function. That is, the Gaussian process is wrapped into the inverse logit function
to yield the probability of a binomial trial of size K.

The model is implemented in Stan and also predicts the value of the process out-of-sample to generate counterfactual outcomes.
These simulations are from the posterior predictive distribution of the process, i.e. they already incorporate any parameter and estimation uncertainty in the model.
