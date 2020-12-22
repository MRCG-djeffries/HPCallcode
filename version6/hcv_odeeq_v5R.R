hcv_odeeq_v5R=function (){
nt=81
M=readMat('F:/HCV/pomp/M.mat')$M
Mdash=readMat('F:/HCV/pomp/Mdash.mat')$Mdash
extra_parms_vals=c(0.133,-log(1-0.027),1/17); # was 5.6, nu was 1/17
phiminusvals1=1;phiminusvals2=22;phiminusvals3=43;phiminusvals4=64;phiminusvals5=85
phiplussvals1=6;phiplussvals2=28;phiplussvals3=49;phiplussvals4=70;phiplussvals5=91
phiminusvals=c(phiminusvals1,phiminusvals2,phiminusvals3,phiminusvals4,phiminusvals5)
phiplussvals=c(phiplussvals1,phiplussvals2,phiplussvals3,phiplussvals4,phiplussvals5)

age_matrix=matrix(0,nrow=9,ncol=9)
age_matrix[1,1] = -1/5
age_matrix[2,1] = 1/5;age_matrix[2,2] = -1/5
age_matrix[3,2] = 1/5;age_matrix[3,3] = -1/5
age_matrix[4,3] = 1/5;age_matrix[4,4] = -1/10
age_matrix[5,4] = 1/10;age_matrix[5,5] = -1/10
age_matrix[6,5] = 1/10;age_matrix[6,6] = -1/10
age_matrix[7,6] = 1/10;age_matrix[7,7] = -1/10
age_matrix[8,7] = 1/10;age_matrix[8,8] = -1/10
age_matrix[9,8] = 1/10
curr_mort_pwid=c(0.96 ,0.96 ,1.12 ,0.18 ,0.22 ,0.53 ,1.38 ,4.28 ,14.96)/1000 # mortality per year for PWID
curr_mort_former=c(0.044 ,0.051 ,0.062 ,0.1 ,0.222 ,0.534 ,1.376 ,4.282 ,14.956 )/1000; # mortality per year for PWID

mort_current = t(matrix(rep(curr_mort_pwid,20),9,20))
mort_former = t(matrix(rep(curr_mort_former,20),9,20))
Irows = c(6,7,8,9,10,11,12,13,14,15) # rows of infectious PWID
death_vec= c(-log(1-0.138), -log(1-0.605), -log(1-0.169), -log(1-0.034)) #DC,HCC,LT1,LT2



dum=matrix(t(rep(seq(360,520,20),10)),ncol=length(Irows))+matrix(rep(Irows,9),ncol=length(Irows),byrow=TRUE)
top10_index = matrix(t(dum),ncol=1)#reshape((repmat(Irows,9,1)+(360:20:520)')',length(Irows)*9,1);
dum=matrix(t(rep(seq(540,700,20),10)),ncol=length(Irows))+matrix(rep(Irows,9),ncol=length(Irows),byrow=TRUE)
top11_index = matrix(t(dum),ncol=1)#reshape((repmat(Irows,9,1)+(540:20:700)')',length(Irows)*9,1);

                                            
#**********************************************************
# Note phi is assumed to be zero for former
#**********************************************************
                                              
# former, failed, ages grps 1 to 9, i=0, j = 0
# current, failed, ages grps 1 to 9, i=1, j = 0
# former, never failed, ages grps 1 to 9, i=0, j = 1
# current, never failed, ages grps 1 to 9, i=1, j = 1


NT=matrix(rep(rep(0,720),nt),ncol=nt)
# starting values
NT[1,1]=560;NT[361,1]=435.6;NT[367,1]=4.4
for (t in 1 : (nt-1)){
  N=NT[,t]
  Mdash=fillval(N,top10_index,top11_index,t,extra_parms_vals,Mdash,phiminusvals,phiplussvals)
  X00 = matrix(N[1:(20*9)],nrow=20,ncol=9)#reshape(N(1:20*9),20,9);
  X01 = matrix(N[181:(180+20*9)],nrow=20,ncol=9)#reshape(N(181:(180+20*9)),20,9);
  X10 = matrix(N[361:(360+20*9)],nrow=20,ncol=9)#reshape(N(361:(360+20*9)),20,9);
  X11 = matrix(N[541:(540+20*9)],nrow=20,ncol=9)#reshape(N(541:(540+20*9)),20,9);
  NT[1:(20*9),(t+1)]=NT[1:(20*9),t] + matrix(M%*%X00-extra_parms_vals[2]*X00+extra_parms_vals[3]*X10-mort_former*X00+t(age_matrix%*%t(X00)),ncol=1)
  NT[181:(180+20*9),(t+1)]=NT[181:(180+20*9),t] + matrix(M%*%X01-extra_parms_vals[2]*X01+extra_parms_vals[3]*X11-mort_former*X01+t(age_matrix%*%t(X01)),ncol=1)
  NT[361:(360+20*9),(t+1)]=NT[361:(360+20*9),t] + matrix(Mdash%*%X10+extra_parms_vals[2]*X00-extra_parms_vals[3]*X10-mort_current*X10+t(age_matrix%*%t(X10)),ncol=1)
  NT[541:(540+20*9),(t+1)]=NT[541:(540+20*9),t] + matrix(Mdash%*%X11+extra_parms_vals[2]*X01-extra_parms_vals[3]*X11-mort_current*X11+t(age_matrix%*%t(X11)),ncol=1)
  NT[361,(t+1)] = NT[361,(t+1)] + sum(mort_former*X00) + sum(mort_former*X01) + sum(mort_current*X10) + sum(mort_current*X11) +
                              death_vec[1]*sum(X00[12,])  + death_vec[2]*sum(X00[13,])  + death_vec[3]*sum(X00[14,])  + death_vec[4]*sum(X00[15,]) +
                              death_vec[1]*sum(X01[12,])  + death_vec[2]*sum(X01[13,])  + death_vec[3]*sum(X01[14,])  + death_vec[4]*sum(X01[15,]) +
                              death_vec[1]*sum(X10[12,])  + death_vec[2]*sum(X10[13,])  + death_vec[3]*sum(X10[14,])  + death_vec[4]*sum(X10[15,]) +
                              death_vec[1]*sum(X11[12,])  + death_vec[2]*sum(X11[13,])  + death_vec[3]*sum(X11[14,])  + death_vec[4]*sum(X11[15,])
}
return(NT)
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
fillval=function(N,top10_index,top11_index,t,extra_parms_vals,Mdash,phiminusvals,phiplussvals){
  top = sum(N[top10_index])+sum(N[top11_index])
  bot = sum(N[361:720])
  I = top/bot;
  
  if (t>=0 && t<51){
    scaleI = 1+1.5*t/51
  }else if (t>51 && t<=56){
    scaleI = 2.5-1.5*(t-51)/5;
  }else{
    scaleI=1;
  }
  
  phi = scaleI*extra_parms_vals[1]*I
  dum = Mdash[5,5]
  endy = length(phiminusvals)
  Mdash[phiminusvals[1:(endy-1)]]=-phi 
  Mdash[phiminusvals[endy]]=-phi+dum 
  Mdash[phiplussvals]=phi
  return(Mdash)
}