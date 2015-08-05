// Fit a Gaussian process's hyperparameters

data {
int<lower=1> N1;    
int<lower=1> N2;    
int<lower=1> N3;    
int K;
vector[N1 + N2 + N3] x;
int<lower=0, upper=K> y2[N2];
}
transformed data {
    int<lower=1> N;             // Total length of the feature vector
    vector[N1 + N2 + N3] mu;         // Prior mean of the process

    N <- N1 + N2 + N3;
    mu <- rep_vector(0, N);
}
parameters {
real<lower=0> beta1;
real<lower=0> beta2;
real<lower=0> beta3;
real<lower=0> beta4;

vector[N1] z1;
vector[N2] z2;
vector[N3] z3;
}
model {
    vector[N] z;
    matrix[N, N] Sigma;
    
    for (a in 1:N1) 
        z[a] <- z1[a];
        
    for (b in 1:N2)
        z[N1 + b] <- z2[b];
    
    for (c in 1:N3)
        z[N1 + N2 + c] <- z3[c];
    
    for(i in 1:N)
        for(j in i:N){
            Sigma[i, j] <- beta1 * exp(-beta2 * pow(x[i] - x[j], 2)) + beta3 + beta4 * x[i] * x[j];
        }
    
    for(i in 1:N)
        for(j in (i+1):N){
            Sigma[j, i] <- Sigma[i, j];
        }

beta1 ~ cauchy(0,2.5);
beta2 ~ cauchy(0,2.5);
beta3 ~ cauchy(0,2.5);
beta4 ~ cauchy(0,2.5);

z ~ multi_normal(mu, Sigma);

for(i in 1:N2)
    y2[i] ~ binomial_logit(K, z2[i]);
} 
generated quantities {
    vector[N] n_sold;
    
    for(i in 1:N1)
        n_sold[i] <- binomial_rng(K, inv_logit(z1[i]));
    for(i in 1:N2)
        n_sold[N1 + i] <- binomial_rng(K, inv_logit(z2[i]));
    for(i in 1:N3)
        n_sold[N1 + N2 + i] <- binomial_rng(K, inv_logit(z3[i]));
}