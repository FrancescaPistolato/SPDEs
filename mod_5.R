## Modello 5: numero 2 in dim 2
# Parametri
{
  Nt = 80000; dt = 0.000001; Nsample = 6
  lam=1; d=0.1; beta=0.1; a=0.2; 
  p=0.3; c=0.3; b=0.1
  nat = (beta*lam*c)/(beta*b+c*d); mor = a
}
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
#
sigma = 1
# Valori iniziali
{
  V = matrix(0,Nsample,Nt)
  I = V
  Vmax = 5
  V[,1] = 0.1
  I[,1] = 0.1
}
# Metodo di Eulero
{
  max_rep = 1
  neq.traj = rep(0,max_rep)
  for (j in (1:max_rep)) {
  # Metodo di Eulero
    for (t in (1:(Nt-1))) {
      V[,t+1] = V[,t] + dt*(nat-mor-p*I[,t])*V[,t] + V[,t]*(Vmax-V[,t])*sqrt(dt)*sigma*rnorm(Nsample)
      I[,t+1] = I[,t] + dt*I[,t]*(c*V[,t]-b) 
      }
    tmp.neg.traj = rep(0,Nsample)
    for (i in (1:Nsample)) {
      tmp.neg.traj[i] = sum(V[i,]<0)+sum(I[i,]<0)
    }
    neq.traj[j] = sum(tmp.neg.traj>0)
    }
  print('In media le traiettorie negative sono')
  print(mean(neq.traj))
}
# Plot
{
  plot(c(0,Nt),c(0,0.5),type="n")
  library(RColorBrewer)
  colV <- brewer.pal(Nsample, "Blues")
  colI <- brewer.pal(Nsample, "Oranges")
  legend("topright",legend = c("Virus","Immune"), col=c(colV[Nsample],colI[Nsample]),lty=rep(1,2),horiz=FALSE, bty='n', cex=0.8)
  for (i in (1:Nsample)) {
    lines(I[i,], col = colI[i])
    lines(V[i,], col = colV[i])
  }
}


