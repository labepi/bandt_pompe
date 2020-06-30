# script to compute the maximum and minimum limits for the CCEP

source('bandt_pompe.R')
source('measures.R')
source('helpers.R')

bp_path = './'

D=7
lenmin=5000
lenmax=1000

# size of probability space
N = factorial(D)
#N = factorial(D)^2

# matrix of barycenters (sub-simplex)
b = matrix(0, ncol=N, nrow=N)

# computing the barycenters
for(i in 1:N)
{
    j = N-(i-1)
    b[i,i:N] = rep(1/j, j)
}

# TODO: evaluate in the intervals between barycenters


# C_min

# the computed metrics (HxC) for C_min curve
d1 = matrix(0, ncol=2, nrow=0)

P = seq(b[1,N], b[N,N], length.out=lenmin)
b2 = b[1,]
for(i in 2:length(P))
{
    # p
    p = P[i]
    b2[N] = p
    # 1-p
    b2[1:(N-1)] = rep((1-p)/(N-1), N-1)
    d1 = rbind(d1,
                c(
                    shannon_entropy(b2, normalized=TRUE),
                    complexity(b2, normalized=TRUE)
                )
            )
}
# adjusting the order of C_min limits
d1 = d1[nrow(d1):1,]

#C_max

# the computed metrics (HxC) for C_max curve
d2 = matrix(0, ncol=2, nrow=0)

for(i in 1:(N-1))
{
    P = seq(b[i,i], b[i+1,i], length.out=lenmax)
    b2 = b[i,]
    #print(P)
    for(j in 2:length(P))
    {
        # p
        p = P[j]
        b2[i] = p
        # 1-p
        b2[(i+1):N] = rep((1-p)/(N-i), N-i)
        d2 = rbind(d2,
                   c(
                     shannon_entropy(b2, normalized=TRUE),
                     complexity(b2, normalized=TRUE)
                     )
                   )
    }
}

# adjusting the order of C_max limits
#d2 = d2[nrow(d2):1,]

# concatenating the min and max curves to save
# NOTE: adding a third column indicating whether the point 
#   is for the C_min or C_max curves
d = rbind(cbind(d1, 1), cbind(d2, 2))

# saving the limits to a file
colnames(d) = c("H", "SC", "Z")
write.table(x=d, file=paste("limits/limits_N",N,".dat",sep=''), row.names=F)


quit()

#pdf('test_limits.pdf')
plot(d[,1], d[,2], col=2, type='l')
#plot.ccep(H=d[,1], SC=d[,2], D=3, ylim=c(0,0.3))
plot.ccep(D=7)
lines(d[,1], d[,2], col=2)
#dev.off()

