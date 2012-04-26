
# entropy function
H <- function(x,sigma) {-sum(dnorm(x,0,sigma)*log2(dnorm(x,0,sigma)))}

nit <- 200
sigma_seq <- seq(0.01,0.1,0.001)
nyr_seq <- 1:20

res <- array(dim=c(length(sigma_seq),length(nyr_seq),nit))

for(s in 1:length(sigma_seq))
  for(y in 1:length(nyr_seq))
    for(i in 1:nit)
      res[s,y,i] <- H(rnorm(y,0,sigma_seq[s]),sigma_seq[s])
      
rmat <- apply(res,1:2,mean)

image(x=sigma_seq,y=nyr_seq,z=rmat,main='Entropy (H)',ylab='n',xlab=expression(sigma),cex=2)
contour(x=sigma_seq,y=nyr_seq,z=rmat,add=T,drawlabels=F,nlevels=20)

savePlot('../res/entropy.pdf',type='pdf')


