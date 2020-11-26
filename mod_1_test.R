# Modello 1: noise limitato
#
Nt = 200; dt = 0.001
Nsample = 10
#
lam=1; d=0.1; beta=0.1; a=0.2; p=0.3; c=0.3; b=0.1
nat = (beta*lam*c)/(beta*b+c*d)
mor = a
#
sc = (lam*c)/(beta*b+c*d)
ic = b/c
ctl = (c*(beta*lam-a*d)-a*beta*d)/(p*(beta*b+c*d))
#
epsilon = 1-p*ctl ## l’ho introdotta per controllare # il fattore p, quello associato all’effetto antivirale # delle ctl. cosa succede se fossero sempre un po’ di
# piu` o un po’ di meno dell’equilibrio?
p1 = (p*ctl)/(epsilon+p*ctl)
#
sigma = 0.1
#
V = matrix(0,Nsample,Nt)
V[,1] = 0.1
#
#
max_rep = 100
neq.traj = rep(0,max_rep)
for (j in (1:max_rep)) {
  # Metodo di Eulero
  for (t in (1:(Nt-1))) {
    V[,t+1] = V[,t] + dt*(nat-mor-p1)*V[,t] + atan((t*(V[,t])))*sqrt(dt)*sigma*rnorm(Nsample)
  }
  tmp.neg.traj = rep(0,Nsample)
  #plot(c(0,Nt),c(-0.2,0),type="n")
  for (i in (1:Nsample)) {
    #lines(V[i,])
    tmp.neg.traj[i] = sum(V[i,]<0)
  }
  neq.traj = sum(tmp.neg.traj>0)
}
print('In media le traiettorie negative sono')
print(mean(neq.traj))



