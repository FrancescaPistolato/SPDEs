# Modello 1: noise limitato
#
Nt = 200; dt = 0.1
Nsample = 10
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
#
V = matrix(0,Nsample,Nt)
V[,1] = 0.1
#
sigma = 0.1
#
# Metodo di Eulero
for (t in (1:(Nt-1))) {
  V[,t+1] = V[,t] + dt*(nat-mor-p1)*V[,t] + 1/(1+V[,t])*sqrt(dt)*sigma*rnorm(Nsample)
}
plot(c(0,Nt),c(-1,1),type="n")
for (i in (1:Nsample)) {
  lines(V[i,])
}