# Modello 6 : aleatorio di questo
#
# 3. CTL replication with susceptible and infected cells
#    based on the simplified and pessimistic model for virus replication
#
# Parametri 
{
  N=10; Nt=80000; dt=0.009
  lam=1; d=0.1; beta=0.1; a=0.2
  p=0.3; c=0.3; b=0.1
}
# Scenari 
{
  R0=(lam*beta)/(d*a)
  Q0=c*(beta*lam-a*d)-a*beta*d
  # 
  # Non c'è una risposta immunitaria:
  if (Q0<=0) {
    if (R0>1) {
      s = lam/d
      i = (lam*beta - d*a)/a*beta
      ctl = 0
    }
    if (R0<1) {
      s = lam/d
      i = 0
      ctl = 0
    }
  }
  # C'è una risposta immunitaria:
  if (Q0>0) {
    s = (lam*c)/(beta*b+c*d)
    i = b/c
    ctl = (c*(beta*lam-a*d)-a*beta*d)/(p*(beta*b+c*d))
  }
}
# Valori iniziali
{
  S=1:Nt; I=S; CTL=S
  S[1]=5; I[1]=0.1; CTL[1]=0.1
}
# Metodo di Eulero
{
  for (t in 1:(Nt-1)){
    S[t+1] = S[t] + dt*(lam - d*S[t] - beta*S[t]*I[t])
    I[t+1] = I[t] + dt*(beta*S[t]*I[t] - a*I[t] - p*I[t]*CTL[t]) 
    CTL[t+1] = CTL[t] + dt*(c*I[t]*CTL[t] - b*CTL[t])  
  }
}
# Plots 
{
  plot(c(0,Nt),c(0,3.5), type="n")
  # lines(SC,col="brown")
  lines(I,col="steelblue")
  lines(CTL,col="indianred")
  legend("topright",legend = c("Virus","Immune"), col=c("steelblue","indianred"),lty=rep(1,2),horiz=FALSE, bty='n', cex=0.8)
}
#
#