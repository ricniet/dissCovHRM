#############################################################################
# EXAMPLE 1: Simulated data using the immer_hrm_simulate() function
#############################################################################

# define data generating parameters
set.seed(1997)
N <- 500  # number of persons
I <- 4    # number of items
R <- 3    # number of raters
K <- 3    # maximum score 
sigma <- 2  # standard deviation
theta <- stats::rnorm( N , sd=sigma )  # abilities
# item intercepts
b <- outer( seq( -1.5 , 1.5 , len=I) , seq( -2 , 2 , len=K)  , "+" )
# item loadings
a <- rep(1,I)
# rater severity parameters
phi <- matrix( c(-.3 , -.2 , .5) , nrow=I , ncol=R , byrow=TRUE )
phi <- phi + stats::rnorm( phi  , sd = .3 )
phi <- phi - rowMeans(phi)
# rater variability parameters
psi <- matrix( c(.1 , .4 , .8) , nrow=I , ncol=R , byrow=TRUE )
# simulate HRM data
data <- immer::immer_hrm_simulate( theta , a , b , phi  = phi, psi = psi )
pid <- data$pid
rater <- data$rater
dat <- data[ , - c(1:2) ]