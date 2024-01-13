#Posterior probability calculator

# we assume a vector of observations as the input
#If the parameter of the distribution is not bounded above or below, and the number of points of the prior distribution
#the prior probability initially will be assumed to have equal probability (discrete uniform)

posterior_calc=function(params_vector, vector_obs, prob_function,disc_or_cont="discrete",cont_cdf=NULL){
  params_vec_length=length(params_vector)
  params_prob_prior=rep(c(1/params_vec_length),params_vec_length)
  mat_obs_probs = matrix(nrow = params_vec_length, ncol = length(vector_obs))
  # if the observations from discrete distribution, (the default).
  if (disc_or_cont=="discrete"){
    # Loop through each params in params_vector and calculate observation probabilities
    mat_obs_probs = sapply(params_vector, function(lambda) prob_function(vector_obs, lambda), simplify = FALSE)
    # Convert the list of vectors into a matrix
    mat_obs_probs = do.call(cbind, mat_obs_probs)
  }else# if the observations from continuous distribution, 
  {
    # Sort observations
    sorted_obs <- sort(vector_obs)
    # Calculate the differences between consecutive elements
    diff <- diff(sorted_obs)
    # Find the minimum difference
    min_diff = min(diff)
    #make delta so that observation +never overlaps, so at least divide min_diff by 2
    # the smaller it is made, the better, so divided by 3,(not too large though)
    delta=min_diff/3
    
    upper_cdf=sapply(params_vector, function(lambda) cont_cdf(vector_obs+delta, lambda), simplify = FALSE)
    
    lower_cdf=sapply(params_vector, function(lambda) cont_cdf(vector_obs-delta, lambda), simplify = FALSE)
    #They need to be turned into a matrix for the subtraction operator to work
    upper_cdf_mat=do.call(cbind, upper_cdf)
    lower_cdf_mat=do.call(cbind, lower_cdf)
    # subtract the cdf to get the "discretized" probability
    mat_obs_probs=upper_cdf_mat-lower_cdf_mat
  }
  
  # Calculate the vector of products
  #(product of the probabilities of observations equaling the realization conditioned on the parameter i is the ith)
  row_prod = apply(mat_obs_probs, 2, prod)
  # Calculate the denominator (sum of products multiplied by the prior probabilities)
  post_prob_denominator = sum(row_prod * params_prob_prior)
  
  # Calculate the numerator (prior probabilities multiplied by the products)
  post_prob_numerator = params_prob_prior * row_prod
  
  # Calculate the posterior probabilities
  post_prob = post_prob_numerator / post_prob_denominator
  return(post_prob)
}


###############################################################
# TWO EXAMPLES FOR DISCRETE bayesian and normal(frequentist) point estimates
# POISSON FOR SINGLE PARAMETER.
# The parameter is not bounded above, so there must be some bound placed based on prior information


# BINOMIAL FOR TWO PARAMETTER
#The method above is implemented in code in a way that only allows a single 
# parameter to be unknown, so in case of two parameter like binomial,
# one of the parameters must be fixed.
###############################################################


#######################
#With poisson (lambda), comparison with the frequentist method
# note sum(vector_obs)/length(vector_obs), the sample mean, is the UMVUE of lambda (frequentist)
#######################

# suppose the parameter true value is smaller than 50 (this is the artificial bound, in practice based on prior information known)
params_vector=seq(1,50)

#15 observation from poisson (25) is the data
vector_obs=rpois(30,25)

#probability function (this is the pmf not CDF) for poisson 
prob_function = function(x, lambda) dpois(x, lambda)


#MAIN CALL TO THE GET POSTERIOR
post_prob=posterior_calc(params_vector, vector_obs, prob_function)

#plotting posterior distribution of parameter
plot(params_vector, post_prob, main = "Scatter Plot", xlab = "value of parameter", ylab = "posterior prob")
#Expected value of the parameter based on posterior probabilities
Bay_est_poisson=sum(post_prob*params_vector)
Bay_est_poisson

#UMVUE is just the sample mean for lambda in poisson(lambda)
UMVUE_poisson=sum(vector_obs)/(length(vector_obs))
UMVUE_poisson

#IN addition, the following are for the point estimates for lambda square
Bay_est_pois_lambda_sqr=sum(post_prob*params_vector^2)
Bay_est_pois_lambda_sqr
#It can be shown that the following is UMVUE for lambda square. ((sumxi)^2-sum(xi))/(n^2)
UMVUE_pois_lamda_sqr=(sum(vector_obs)^2-sum(vector_obs))/(length(vector_obs)^2)
UMVUE_pois_lamda_sqr


#######################
#With Binomial (5,lambda), comparison with the frequentist method
# note sum(vector_obs)/(length(vector_obs)*5), sample mean divided by 5, is the UMVUE of lambda (frequentist)
#######################

#parameter support at each 0.1 increment.
# If there is more information then the prior support can be modified directly accordingly
params_vector=seq(0.1,0.9,by=0.1)

#15 observations of binomial (5,0.25) as the  data
vector_obs=rbinom(15,5,0.25)

#note that for two parameter distribution, one of the parameter should be known and fixed
#binoimial distribution here has fixed n as 5.
prob_function = function(x, lambda) dbinom(x, 5,lambda)

#MAIN CALL TO THE GET POSTERIOR
post_prob=posterior_calc(params_vector, vector_obs, prob_function)

#plotting the posterior distribution of the parameter
plot(params_vector, post_prob, main = "Scatter Plot", xlab = "value of parameter", ylab = "posterior prob")

#Expected value of the parameter based on posterior probabilities
Bay_est_binom=sum(post_prob*params_vector)
Bay_est_binom

#UMVUE
UMVUE_binom=sum(vector_obs)/(length(vector_obs)*5)
UMVUE_binom



###############################################################
#Two continuous examples
# Nrmal (lambda, 1)
#Uniform (0, theta))
#parameters need to have bound, used -50 to 50 for normal, and 50 for uniform
###############################################################




#######################
#With Normal(lambda,1), comparison with the frequentist method
# note sum(vector_obs)/length(vector_obs), the sample mean, is the UMVUE of lambda (frequentist)
#######################


#parameter support from -50 to 50 at 0.5 increment. Increment is based on knowing sigma is 2
# If there is more information then the prior support can be modified directly accordingly
params_vector=seq(-50,50,by=0.5)

#15 observations of normal (29,2) as the  data
vector_obs=rnorm(15,29,2)

#note that for two parameter distribution, one of the parameter should be known and fixed
#binoimial distribution here has fixed n as 5.
prob_function = function(x, lambda) dnorm(x, lambda,2)
cont_cdf= function(x,lambda) pnorm(x,lambda,2)



#MAIN CALL TO THE GET POSTERIOR
post_prob=posterior_calc(params_vector, vector_obs, prob_function,"continuous",cont_cdf)

#plotting the posterior distribution of the parameter
plot(params_vector, post_prob, main = "Scatter Plot", xlab = "value of parameter", ylab = "posterior prob")

#Expected value of the parameter based on posterior probabilities
Bay_est_norm=sum(post_prob*params_vector)
Bay_est_norm

#UMVUE
UMVUE_norm=sum(vector_obs)/(length(vector_obs))
UMVUE_norm

#Again as with the poisson, to get the bayes point estimate
#for the lambda^2, we can simply just square the params_vector as below
Bay_est_norm_lambda_sqr=sum(post_prob*params_vector^2)
Bay_est_norm_lambda_sqr

#It can be shown the UMVUE for the lambda^2 is (sum(xi))^2-n)/n^2
UMVUE_norm_lambda_sqr=(sum(vector_obs)^2-length(vector_obs))/(length(vector_obs)^2)
UMVUE_norm_lambda_sqr



#######################
#Uniform (0,theta)
# (1+1/n)*max{observations} is UMVUE (frequentist)
#######################

#believe that the parameter theta should be lower than 50.
params_vector=seq(1,50,by=1)

#15 observations of normal (29,2) as the  data
vector_obs=runif(10,0,32)

#note that for two parameter distribution, one of the parameter should be known and fixed
#binoimial distribution here has fixed n as 5.
prob_function = function(x, lambda) dunif(x,0,lambda)
cont_cdf= function(x,lambda) punif(x,0,lambda)


#MAIN CALL TO THE GET POSTERIOR
post_prob=posterior_calc(params_vector, vector_obs, prob_function,"continuous",cont_cdf)

#plotting the posterior distribution of the parameter
plot(params_vector, post_prob, main = "Scatter Plot", xlab = "value of parameter", ylab = "posterior prob")

#Expected value of the parameter based on posterior probabilities
Bay_est_unif=sum(post_prob*params_vector)
Bay_est_unif

#UMVUE
UMVUE_unif=(1+1/length(vector_obs))*(max(vector_obs))
UMVUE_unif




