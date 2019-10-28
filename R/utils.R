sample_edges <- function(A) {

}

sample_elements <- function(X) {

}

set.seed(27)

n <- 1000
k <- 5

B <- matrix(runif(k * k), nrow = k, ncol = k) # mixing probabilities

theta <- round(rlnorm(n, 2)) # degree parameter
pi <- c(1, 2, 4, 1, 1) / 5 # community membership

A <- dcsbm(theta, pi, B, avg_deg = 50)


library(dplyr)

# sampling proportional to absolute value of edge weight


# get a 5% elementwise sample of a sparse matrix





library(Matrix)

alpha = 0.05

A <- rsparsematrix(1000, 1000, 0.1)

s <- as.data.frame(summary(A))
p <- abs(s$x) / sum(abs(s$x))   # probabilities proportional to abs(entry)

index <- sample.int(nrow(s), size = round(nrow(s) * alpha), prob = p)

sampled <- s[index, ]

sparseMatrix(sampled$i, sampled$j, x = sampled$x)






train_edges <- s[-index, ]
test_edges  <- s[index, ]

train <- sparseMatrix(train_edges$i, train_edges$j, x = train_edges$x)
test <- sparseMatrix(test_edges$i, test_edges$j, x = test_edges$x)

list(train = train, test = test)

nnzero(train)
nnzero(test)
nnzero(A)

Matrix()

s$split <- sample()
sample_frac(s, size = 0.05, weight = x)

