library(GoFKernel)
library(numDeriv)
library(rootSolve)
library(stats19)
library(MASS)
library(bbmle)
library(maxLik)
library(cubature)
library(ggplot2)



a=0.1;b=0.5;c=1
w<-runif(1)
u<-runif(1,0,1-w)
op <- par(mar = c(4,4,4,2)- 0.1)
ff1<-function(t){
if(t<a)return(0)
if(t>=a&t<=b) return(w*(t-a)/(b-a))
if(t>b&t<=c) return(w*(c-t)/(c-b))
if(t>c)return(0)}
functions<-Vectorize(ff1)
curve(functions,0.01,1.2,xlim=c(0.01,1.1),
ylim=c(0,1.2),xlab="x",xaxt="n",
ylab="functions",frame.plot=FALSE,
axes =FALSE,lwd=2)
axis(side=1, pos=0,at=c(0.1,.5,1),
labels=letters[1:3],lwd.ticks=1)
axis(side = 2, pos=0,label=FALSE,
lwd.ticks=0)
axis(side=2, pos=0,at=c(0,u,w,1),
labels =c("",expression("u"[A]),
expression("w"[A]),1),lwd.ticks=1)
segments(0,w,b,w,col="gray",lty="dashed"
,lwd=2)
ff2<-function(t){
if(t<a)return(1)
if(t>=a&t<=b) return((b-t+u*(t-a))/(b-a))
if(t>b&t<=c) return((t-b+u*(c-t))/(c-b))
if(t>c)return(1)}
functions<-Vectorize(ff2)
curve(functions,0.01,1.2,xlim=c(0.01,1.1),
col=2,add=T,lwd=2)
legend(0.15,1.15,
legend=c("Membership function",
"Non-membership function"),
lty=c(1,1),col=c(1,2),bty=’n’)
segments(0,u,0.5,u,col="gray",
lty="dashed",lwd=2)
segments(-0.005,0,0.1,0,lwd=2)
rm(list=ls(all=TRUE))
library(ggplot2)
gamma<-5;beta<-10;M<-50;n<-30;a<-c()
b<-c();w<-c();u<-c();H<-c()
y<-rweibull(n,gamma,beta)
L<-matrix(ncol=M,nrow=M)
y<-rweibull(n, shape=gamma, scale=beta)
for(i in 1:n){
a[i]=runif(1,y[i]-1,y[i])
b[i]=runif(1,y[i],y[i]+1)
w[i]<-runif(1)
u[i]<-runif(1,0,1-w[i])
H[i]<-(1+w[i]-u[i])/2}
A <- seq(2,10, len = M);
B <- seq(5,15, len = M)
for(m in 1:M){for(k in 1:M){
gamma<-A[m]
beta<-B[k]
f<-function(t){
f<-dweibull(t,shape=gamma,
scale=beta)}
mu<-function(t,i){
if(t>=a[i] & t<=y[i])
return(H[i]*(t-a[i])/(y[i]-a[i]))
if(t>y[i] & t<=b[i])
return(H[i]*(b[i]-t)/(b[i]-y[i]))
else
return(0)}
mu<-Vectorize(mu)
h1<-function(t,i) {
h1<-mu(t,i)*f(t)}
I<-function(z){
i<-z[1]
I<-log(integrate(h1,a[i],b[i],i)$value)}
ss<-0
for(j in 1:n){
c<-c(beta,gamma,j)
ss<-ss+I(c)}
L[m,k]<-ss}}
df <- data.frame(A,B,L)
gamma<-A;beta<-B
ggplot(df, aes(x = A, y = B)) +
geom_density_2d(col="blue2")
+theme(axis.title=element_text(size=18,
face="bold"), plot.background =
element_rect(fill = "white"),
axis.text=element_text(size=18),
panel.background = element_rect(fill =
"white"), axis.line.x = element_line(
color = "grey2"),axis.line.y =
element_line(color = "grey2"))+
labs(x=expression(paste(gamma)))+
labs(y=expression(paste(beta)))+
geom_segment(data=df,
mapping=aes(x=5.32, y=4,
xend=5.32, yend=9.47), size=1,
col="red",lty="dashed")
+geom_segment(data=df,
mapping=aes(x=2, y=9.47,
xend=5.31, yend=9.47), size=1,col="red",lty="dashed")+
geom_text(aes(x=5.32, label=paste(5.32),y=3.5),col="red",cex=7)+geom_text(aes(y=10, label=paste(9.47), x=2.1),col="red",cex=7)



a=0.1;b=0.5;c=1
w<-runif(1)
u<-runif(1,0,1-w)
ff1<-function(t){
if(t<a)return(0)
if(t>=a&t<=b) return(w*(t-a)/(b-a))
if(t>b&t<=c) return(w*(c-t)/(c-b))
if(t>c)return(0)}
functions<-(ff1)

ff2<-function(t){
if(t<a)return(1)
if(t>=a&t<=b) return((b-t+u*(t-a))/(b-a))
if(t>b&t<=c) return((t-b+u*(c-t))/(c-b))
if(t>c)return(1)}
functions<-(ff2)
rm(list=ls(all=TRUE))
gamma<-5;beta<-10;M<-50;n<-30;a<-c()
b<-c();w<-c();u<-c();H<-c()
y<-rweibull(n,gamma,beta)
L<-matrix(ncol=M,nrow=M)
y<-rweibull(n, shape=gamma, scale=beta)
for(i in 1:n){
a[i]=runif(1,y[i]-1,y[i])
b[i]=runif(1,y[i],y[i]+1)
w[i]<-runif(1)
u[i]<-runif(1,0,1-w[i])
H[i]<-(1+w[i]-u[i])/2}
A <- seq(2,10, len = M);
B <- seq(5,15, len = M)
for(m in 1:M){for(k in 1:M){
gamma<-A[m]
beta<-B[k]
f<-function(t){
f<-dweibull(t,shape=gamma,
scale=beta)}
mu<-function(t,i){
if(t>=a[i] & t<=y[i])
return(H[i]*(t-a[i])/(y[i]-a[i]))
if(t>y[i] & t<=b[i])
return(H[i]*(b[i]-t)/(b[i]-y[i]))
else
return(0)}
mu<-Vectorize(mu)
h1<-function(t,i) {
h1<-mu(t,i)*f(t)}
I<-function(z){
i<-z[1]
I<-log(integrate(h1,a[i],b[i],i)$value)}
ss<-0
for(j in 1:n){
c<-c(beta,gamma,j)
ss<-ss+I(c)}
L[m,k]<-ss}}
df <- data.frame(A,B,L)
gamma<-A;beta<-B

n=30
gamma=3;beta<-10
x<-rweibull(n, shape=gamma, scale =beta)
a<-c();b<-c();w<-c();u<-c()
for(i in 1:n){
a[i]<-runif(1,x[i]-3,x[i])
b[i]<-runif(1,x[i],x[i]+3)}
h<-1
while(h<n){
w[h]<-runif(1);u[h]<-runif(1)
if(w[h]+u[h]<=1){h<-h+1}
else {h<-h}}

#EM Algo

rm(list=ls(all=TRUE))
n=30;gamma<-3;beta<-10
y<-c();a<-c();b<-c();w<-c()
u<-c();H<-c()
y<-rweibull(n,gamma,beta)
bethat<-c();gamhat<-c()
bett<-c();gamm<-c()
bett[1]<-beta-1.5
gamm[1]<-gamma-0.5
eps1<-c();eps2<-c()
eps1<-c();eps2<-c()
eps1[1]<-1;eps2[1]<-1
for(it in 1:1000){
hh<-1
while(eps1[hh]>0.3| eps2[hh]>0.3){
for(i in 1:n){
a[i]=runif(1,y[i]-1,y[i])
b[i]=runif(1,y[i],y[i]+1)
w[i]<-runif(1)
u[i]<-runif(1,0,1-w[i])
H[i]<-(1+w[i]-u[i])/2}
f<-function(t,beta,gamma){
f<-dweibull(t,shape=gamma,
scale=beta)}
f1<-function(t,beta,gamma){
f1<-log(t/beta)}

f2<-function(t,beta,gamma){
f2<-((t/beta)^(gamma)*log(t/beta))}
f3<-function(t,beta,gamma){
f3<-(t^gamma)/(beta^(gamma+1))}
mu<-function(t,i){
if(t>=a[i] & t<=y[i])
return(H[i]*(t-a[i])/(y[i]-a[i]))
if(t>y[i] & t<=b[i])
return(H[i]*(b[i]-t)/(b[i]-y[i]))
else
return(0)}
mu<-Vectorize(mu)
h1<-function(t,i,beta,gamma){
mu(t,i)*f(t,beta,gamma)}
I1<-function(beta,gamma,i){
I1<-integrate(h1,max(0,a[i]),b[i],
i,beta,gamma)$value}
e1<-function(t,i,beta,gamma){
f1(t,beta,gamma)*h1(t,i,beta,
gamma)/I1(beta,gamma,i)}
E1<-function(beta,gamma,i){
E1<-integrate(e1,max(0,a[i]),
b[i],i,beta,gamma)$value}
e2 <- function(t, i, beta, gamma) {
  value <- f2(t, beta, gamma) * h1(t, i, beta, gamma) / I1(beta, gamma, i)
  return(value)
}
E2<-function(beta,gamma,i){
E2<-integrate(e2,max(0,a[i]),
b[i],i,beta,gamma)$value}
e3 <- function(t, i, beta, gamma) {
  value <- f3(t, beta, gamma) * h1(t, i, beta, gamma) / I1(beta, gamma, i)
  return(value)
}
E3<-function(beta,gamma,i){
E3<-integrate(e3,max(0,a[i]),b[i],i,
beta,gamma)$value}
ss1<-0;ss2<-0;ss3<-0
for(j in 1:n){
ss1<-suppressWarnings(ss1+E1(beta,gamma,j))
ss2<-suppressWarnings(ss2+E2(beta,gamma,j))
ss3<-suppressWarnings(ss3+E3(beta,gamma,j))}
E.1<-ss1;E.2<-ss2;E.3<-ss3
hh<-hh+1
gamm[hh]<-n/(E.2-E.1)
bett[hh]<-n/(E.3)
eps1[hh]<-abs(bett[hh]-bett[hh-1])
eps2[hh]<-abs(gamm[hh]-gamm[hh-1])}
gamhat[it]<-gamm[hh]
bethat[it]<-bett[hh]}
mean(gamhat);mean(bethat)
mean(gamhat)-gamma;mean(bethat)-beta
mean((gamhat-gamma)^2)
mean((bethat-beta)^2)
x<-5
R<-exp(-((x/beta)^gamma))
print(R)
RR<-exp(-((x/bethat)^gamhat))
print(mean(RR))
print(mean(RR-R))
print(mean((RR-R)^2))



#NR Method

rm(list=ls(all=TRUE))
n=100;gamma<-5;beta<-10;y<-c()
a<-c();b<-c();w<-c();u<-c();H<-c()
y<-rweibull(n,gamma,beta)
bethat<-c();gamhat<-c()
for(it in 1:1000){
for(i in 1:n){
a[i]=runif(1,y[i]-1,y[i])
b[i]=runif(1,y[i],y[i]+1)
w[i]<-runif(1)
u[i]<-runif(1,0,1-w[i])
H[i]<-(1+w[i]-u[i])/2}
f<-function(t,beta,gamma){
f<-dweibull(t,shape=gamma,scale=beta)}

mu<-function(t,i){
if(t>=a[i] & t<=y[i])
return(H[i]*(t-a[i])/(y[i]-a[i]))
if(t>y[i] & t<=b[i])
return(H[i]*(b[i]-t)/(b[i]-y[i]))
else
return(0)}
mu<-Vectorize(mu)

h1<-function(t,i,beta,gamma){
h1<-mu(t,i)*f(t,beta,gamma)}
I<-function(z){
beta<-z[1]
gamma<-z[2]
i<-z[3]
I<-log(integrate(h1,max(0,a[i]),
b[i],i,beta,gamma)$value)}
logl<-function(z){
beta<-z[1]
gamma<-z[2]
ss<-0
for(j in 1:n){
c<-c(beta,gamma,j)
ss<-ss+I(c)}
logl<--ss}
c<-c(beta,gamma)
out<-suppressWarnings(nlminb(c,logl))
bethat[it]<-out$par[1]
gamhat[it]<-out$par[2]}
mean(gamhat);mean(bethat)
mean(gamhat)-gamma;mean(bethat)-beta
mean((gamhat-gamma)^2)
mean((bethat-beta)^2)
x<-5
R<-exp(-((x/beta)^gamma))
print(R)
RR<-exp(-((x/bethat)^gamhat))
print(mean(RR))
print(mean(RR-R))
print(mean((RR-R)^2))





