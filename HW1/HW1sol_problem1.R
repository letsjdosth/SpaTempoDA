# R code created by Choi Seokjun
# UTF-8 CRLF
# windows 10 64-bit.
# R 3.5.2

# vscode run as source : ctrl+shift+s
# vscode run selected block only : ctrl+enter


# 1. Suppose we want to simulate a random vector Y ∼ N(µ,Σ). If Σ (Matern) is symmetric and positive definite, 
# it can be represented using the Cholesky decomposition Σ = LL', where L is a lower triangular matrix. 
# Consider the following algorithm for simulating Y : 
# (1) Calculate the matrix L. 
# (2) Sample Z ∼ N(0,I), where I is the n×n identity matrix. 
# (3) Let Y = µ + LZ. 

# (a) Show that Y generated in this way has the correct distribution
# (b) Write a function or a few lines of code in R to implement this method for arguments mu and Sigma. 
# You may use the built-in function chol for the Cholesky decomposition and rnorm to generate Z.
# (c) For a mean and covariance function of your choosing, 
# use your code from (b) and make a few plots illustrating realizations of a Gaussian process on [0;1], 
# but changing the diﬀerent parameters in the model. These diﬀerences will be easier to see if you keep the same Z sample but just change mu and Sigma.


rnorm_chol_algorithm = function(mu, Sigma) {
    # mu is vector
    # Sigma should be symmetric, positive definite matrix
    dimension = dim(Sigma)
    vecZ = rnorm(dimension, 0, 1)
    chol_upperL = chol(Sigma)
    return(mu + t(chol_upperL) %*% vecZ)    
}

test_suite_for_2d = function(mu_2dim, cov_2by2dim, sim_num=10000){
    test_result_x1 = rep(0, sim_num)
    test_result_x2 = rep(0, sim_num)
    for(i in 1:sim_num){
        iter_sample = rnorm_chol_algorithm(mu_2dim, cov_2by2dim)
        test_result_x1[i] = iter_sample[1]
        test_result_x2[i] = iter_sample[2]
    }
    par(mfrow=c(3,1))
    plot(test_result_x1, test_result_x2, main='scatterplot', xlab='1st dim', ylab='2nd dim')
    plot(seq(0,1,length.out=sim_num), test_result_x1, type='l', xlab='iteration(normalized)', ylab='value')
    plot(seq(0,1,length.out=sim_num), test_result_x2, type='l', xlab='iteration(normalized)', ylab='value')
}

# test 1 : ordinary gaussian
test1_mu = c(1, -1)
test1_Sigma = matrix(c(1,-1,-1,2),2,2)
test_suite_for_2d(test1_mu, test1_Sigma)


# for test 2,3,4, I use 
# fields$Matern
# if you don't have 'fields' packages, install first and run below.
library(fields)


# test 2 : with Matern covariance (= exp cov)
test2_mu = c(1, -1)

test2_matern_param_nu = 0.5
test2_matern_param_range = 1
test2_distance_mat = matrix(c(0,0.2,0.2,0),2,2)
test2_Sigma = Matern(test2_distance_mat, test2_matern_param_range, test2_matern_param_nu)
# nu = smoothness in lec.note, range=rho?? in lec.note (maybe!)
print("set covariance matrix of test 2")
print(test2_Sigma)
test_suite_for_2d(test2_mu, test2_Sigma)



# test 3 : with Matern covariance (small nu)
test3_mu = c(1, -1)

# library(fields)
# fields$Matern
test3_matern_param_nu = 0.1
test3_matern_param_range = 1
test3_distance_mat = matrix(c(0,0.2,0.2,0),2,2)
test3_Sigma = Matern(test3_distance_mat, test3_matern_param_range, test3_matern_param_nu)
print("set covariance matrix of test 3")
print(test3_Sigma)
test_suite_for_2d(test3_mu, test3_Sigma)


# test 4 : with Matern covariance (big nu)
test4_mu = c(1, -1)

# library(fields)
# fields$Matern
test4_matern_param_nu = 1
test4_matern_param_range = 1
test4_distance_mat = matrix(c(0,0.2,0.2,0),2,2)
test4_Sigma = Matern(test4_distance_mat, test4_matern_param_range, test4_matern_param_nu)
print("set covariance matrix of test 4")
print(test4_Sigma)
test_suite_for_2d(test4_mu, test4_Sigma)

