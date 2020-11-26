## Modello 4: numero 1 in dim 2
Nt = 20000; dt = 0.001
Nsample = 1
#
lam=1; d=0.1; beta=0.1; a=0.2; p=0.3; c=0.3; b=0.1
nat = (beta*lam*c)/(beta*b+c*d)
mor = a
{
  R0=(lam*beta)/(d*a)
  Q0=c*(beta*lam-a*d)-a*beta*d
  #
  # Non c'è una risposta immunitaria:
  if (Q0<0) {
    if (R0>1) {
      sc = lam/d
      ic = (lam*beta - d*a)/a*beta
      ctl = 0
    }
    if (R0<1) {
      sc = lam/d
      ic = 0
      ctl = 0
    }
  }
  # C'è una risposta immunitaria:
  if (Q0>0) {
    sc = (lam*c)/(beta*b+c*d)
    ic = b/c
    ctl = (c*(beta*lam-a*d)-a*beta*d)/(p*(beta*b+c*d))
  }
}
epsilon = 1-p*ctl
p1 = (p*ctl)/(epsilon+p*ctl)
#
sigma = 0.1
#
V = matrix(0,Nsample,Nt)
I = V
V[,1] = 0.1
I[,1] = 0.1
#
#
# Metodo di Eulero
for (t in (1:(Nt-1))) {
  V[,t+1] = V[,t] + dt*(nat-mor-p1)*V[,t] + atan((t*(V[,t])))*sqrt(dt)*sigma*rnorm(Nsample)
  I[,t+1] = I[,t] + dt*I[,t]*(c*V[,t]-b)
  }
#tmp.neg.traj = rep(0,Nsample)
plot(c(0,Nt),c(0,1.5),type="n")
legend("topleft",legend = c("Virus","Immune"), col=c("green","pink"),lty=rep(1,2),horiz=FALSE, bty='n', cex=0.8)
for (i in (1:Nsample)) {
  lines(V[i,], col = "green")
  lines(I[i,], col = "pink")
  #tmp.neg.traj[i] = sum(V[i,]<0)
}





