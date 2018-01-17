#H时变函数，all_theta表示时变参数,nu自由度
Hfunc1=function(family,u,v,all_theta,nu){
  UMAX=1-1e-10
  UMIN=1e-10
  n=length(u)
  XEPS=1e-4
  h =rep(NA,length=n)
  for(j in 1:n){
    theta=all_theta[j]
    if(family==0) {h[j]=u[j]}#independent
    pow=function(a,b){return(a^b)}
    if(family==1) 
    {
      x = (qnorm(u[j],0.0,1.0,1,0) - theta*qnorm(v[j],0.0,1.0,1,0))/sqrt(1.0-pow(theta,2.0))
      h[j] = pnorm(x,0.0,1.0,1,0)
    } #gaussian
    if(family==2)#student
    {
      t1 = qt(u[j],nu,1,0)
      t2 = qt(v[j],nu,1,0)
      mu =theta*t2
      sigma2 = ((nu+t2*t2)*(1.0-theta*(theta)))/(nu+1.0)
      h[j] = pt((t1-mu)/sqrt(sigma2),nu+1.0,1,0)
    }
    if(family==3) { #clayton
      x = pow(u[j],-theta)+pow(v[j],-theta)-1.0 
      h[j] =pow(v[j],-theta-1.0)*pow(x,-1.0-1.0/(theta))
      if(theta < 0)
      {
        if(x < 0) h[j] = 0
      }
    }
    if(family==4) {#gumbel
      h[j] = -(exp(-pow(pow(-log(v[j]),theta)+pow(-log(u[j]),theta),1.0/(theta)))*
                 pow(pow(-log(v[j]),theta)+pow(-log(u[j]),theta),1.0/(theta)-1.0)*
                 pow(-log(v[j]),theta))/(v[j]*log(v[j]))
    }
    if(family==5){ #frank
      h[j] = -(exp(theta)*(exp(theta*u[j])-1.0))/
        (exp(theta*v[j]+theta*u[j])-exp(theta*v[j]+theta)-exp(theta*u[j]+theta)+exp(theta))
    }
    
    if(family==6){#joe
      h[j] = pow(pow(1.0-u[j],theta) + pow(1.0-v[j],theta) - pow(1.0-u[j],theta)*pow(1.0-v[j],theta),1.0/(theta)-1) * 
        pow(1.0-v[j],theta-1.0)*(1-pow(1-u[j],theta))
    }
    if(family==13){# rotated clayton
      u[j]=1-u[j]; 
      v[j]=1-v[j];
      x = pow(u[j],-theta)+pow(v[j],-theta)-1.0 ;
      h[j] =   pow(v[j],-theta-1.0)*pow(x,-1.0-1.0/(theta)); 
      h[j]= 1-h[j];
      u[j]=1-u[j];
      v[j]=1-v[j];
    }
    if(family==14){ #rotated gumbel (180°)
      v[j]= 1-v[j]
      u[j]= 1-u[j]
      h[j]= -(exp(-pow(pow(-log(v[j]),theta)+pow(-log(u[j]),theta),1.0/(theta)))*
                pow(pow(-log(v[j]),theta)+pow(-log(u[j]),theta),1.0/(theta)-1.0)*
                pow(-log(v[j]),theta))/(v[j]*	log(v[j]))
      h[j]= 1-h[j]
      u[j]=1-u[j]
      v[j]=1-v[j]
    }
    if(family==16){
      v[j]= 1-v[j]
      u[j]= 1-u[j]
      h[j] = pow(pow(1.0-u[j],theta) + pow(1.0-v[j],theta) - pow(1.0-u[j],theta)*
                   pow(1.0-v[j],theta),1.0/(theta)-1) * pow(1.0-v[j],theta-1.0)*
        (1-pow(1-u[j],theta))
      h[j]= 1-h[j]
      u[j]=1-u[j]
      v[j]=1-v[j]
    }
  }
  return(h)
}
#因为只有一个数据跟静态一样
HNumInv=function(family,u,v,theta,nu)
{
  br=0
  ans=0.0
  tol=0.000001
  x0=1e-10
  x1=1-1e-10
  fl=0.0 
  fh=0.0
  val=0.0
  fl=Hfunc1(family,x0,v,theta,nu)
  fl=fl-u 
  fh=Hfunc1(family,x1,v,theta,nu)
  fh=fh-u
  if(abs(fl)<=tol) { ans=x0; br=1; }
  if(abs(fh)<=tol) { ans=x1; br=1; }
  while(br!=1){
    ans = (x0+x1)/2.0;
    val=Hfunc1(family,ans,v,theta,nu)
    val=val-u;
    if(abs(val)<=tol) br=1;
    if(abs(x0-x1)<=1e-10) br=1; #stop if values become too close (avoid infinite loop)
    if(val > 0.0) {x1 = ans; 
    fh=val;}else{
      x0=ans; fl = val;
    }
    
  }
  return(out=ans);
} 

#H时变逆函数，all_theta表示时变参数,nu自由度
Hinv=function(family,u,v,all_theta,nu){
  n=length(u)
  hinv =rep(NA,n);
  UMAX=1-1e-10;
  UMIN=1e-10;
  XEPS=1e-4;
  out=rep(NA,n)
  pow=function(a,b){return(a^b)}
  for(i in 1:n){
    if(u[i]<UMIN) u[i]=UMIN;
    if(u[i]>UMAX) u[i]=UMAX;
    if(v[i]<UMIN) v[i]=UMIN;
    if(v[i]>UMAX) v[i]=UMAX;
  }
  for(j in 1:n){
    theta=all_theta[j]
    if(family==0){
      hinv[j]=u[j]
    }
    if(family==1)#gaussian
    {
      hinv[j] = pnorm(qnorm(u[j],0.0,1.0,1,0)*sqrt(1.0-pow(theta,2.0))+
                        theta*qnorm(v[j],0.0,1.0,1,0),0.0,1.0,1,0)
    }
    if(family==2)#student
    {
      temp1 = qt(u[j],nu+1.0,1,0)
      temp2 = qt(v[j],nu,1,0)
      mu = theta*temp2
      var=((nu+(temp2*temp2))*(1.0-(theta*(theta))))/(nu+1.0)
      hinv[j] = pt((sqrt(var)*temp1)+mu,nu,1,0)
    }
    else if(family==3)#clayton
    {
      if(theta <XEPS){hinv[j]=u[j]}else{
        hinv[j] = pow(pow(u[j]*pow(v[j],theta+1.0),-theta/
                            (theta+1.0))+1.0-pow(v[j],-theta),-1.0/(theta))
        }
    }
    if(family==4)#gumbel - must turn to numerical inversion
    {
      nu=0
      hinv[j]=HNumInv(family,u[j],v[j],theta,nu)
    }
    if(family==5)#frank
    {
      hinv[j]=-1/(theta)*log(1-(1-exp(-theta)) / ((1/u[j]-1)*exp(-theta*v[j])+1))
    }
    if(family==6)#joe - numerical inversion
    {
      nu=0.0
      hinv[j]=HNumInv(family,u[j],v[j],theta,nu)
    }
    if(family==13){
      u[j]=1-u[j]
      v[j]=1-v[j]
      hinv[j] = pow(pow(u[j]*pow(v[j],theta+1.0),-theta/(theta+1.0))+
                      1.0-pow(v[j],-theta),-1.0/(theta))
      hinv[j]=1-hinv[j]
      u[j]=1-u[j]
      v[j]=1-v[j]
    }
    if(family==14) #rotated gumbel (180°) 
    {
      u[j]=1-u[j]
      v[j]=1-v[j]
      hinv[j]=HNumInv(4,u[j],v[j],theta,nu)
      hinv[j]=1-hinv[j]
      u[j]=1-u[j]
      v[j]=1-v[j]
    }
    if(family==16){
      u[j]=1-u[j]
      v[j]=1-v[j]			
      hinv[j]=HNumInv(6,u[j],v[j],theta,nu)			
      hinv[j]=1-hinv[j]
      u[j]=1-u[j]
      v[j]=1-v[j]
    }
    out[j] = max(min(hinv[j],UMAX),UMIN) 
  }
  return(out)
}  

indep_test <- function(u1, u2) {
  tau <- cor(u1, u2,method = "kendall")
  N <- length(u1)
  f <- sqrt((9 * N * (N - 1))/(2 * (2 * N + 5))) * abs(tau)
  return( p.value = 2 * (1 - pnorm(f)))
} 
