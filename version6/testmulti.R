testmulti=function(){
  library(caRamel)

  nvar = 3 # number of variables
  bounds <- matrix(data = 1, nrow = nvar, ncol = 2) # upper and lower bounds
  bounds[1,]=c(0.01,0.5)
  bounds[2,]=c(0.001,0.1)
  bounds[3,]=c(1/20,1/2)
  nobj <- 2 # number of objectives
  minmax <- c(FALSE, FALSE) # min and min
  popsize <- 100 # size of the genetic population
  archsize <- 5 # size of the archive for the Pareto front
  maxrun <- 10000 # maximum number of calls
  prec <- matrix(1.e-3, nrow = 1, ncol = nobj) # accuracy for the convergence phase
  
  results <- 
    caRamel(nobj,
            nvar,
            minmax,
            bounds,
            asode,
            popsize,
            archsize,
            maxrun,
            prec,
            carallel=FALSE,
            numcores = 2,
            graph = TRUE) # no parallelism
  saveRDS(list(R=results),'F:/HCV/version6/multi3.rds')
  
}
asode=function(i){
  library(deSolve)
  library(ggplot2)
  library(R.matlab)
  # get the parameters
  M=readMat('F:/HCV/pomp/M.mat')$M
  Mdash=readMat('F:/HCV/pomp/Mdash.mat')$Mdash
  Mvec=as.vector(t(M))
  Mdashvec=as.vector(t(Mdash))
  #L=coeffsof_v5()
  #saveRDS(L,file="F:/HCV/version6/Lcomps.RDS")
  L=readRDS(file="F:/HCV/version6/Lcomps.RDS")
  #extra_parms_nams={'piv','relapse','nu'}; % infection rate (piv instead of pi), relapse to IDU, 1/duration of injecting span
  #extra_parms_vals=c(0.133,-log(1-0.027),1/17); # was 5.6, nu was 1/17
  #extra_parms_vals=c(0.233,-log(1-0.007),1/10); # was 5.6, nu was 1/17
  phiminusvals1=1;phiminusvals2=22;phiminusvals3=43;phiminusvals4=64;phiminusvals5=85
  phiplussvals1=6;phiplussvals2=28;phiplussvals3=49;phiplussvals4=70;phiplussvals5=91
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
  Magevec=as.vector(t(age_matrix))
  curr_mort_pwid=c(0.96 ,0.96 ,1.12 ,0.18 ,0.22 ,0.53 ,1.38 ,4.28 ,14.96)/1000 # mortality per year for PWID
  curr_mort_former=c(0.044 ,0.051 ,0.062 ,0.1 ,0.222 ,0.534 ,1.376 ,4.282 ,14.956 )/1000; # mortality per year for PWID
  death_rate_dc= 0.138
  death_rate_hcc=0.605
  death_rate_lt1=0.169
  death_rate_lt2=0.034
  extra_parms_vals=c(x[i,1],x[i,2],x[i,3])# place holder, these are the calibration parameters
  #cat(paste0(x[1],",",x[1],",",x[1],"\n"))
  param_vals=c(Mvec,Mdashvec,extra_parms_vals,
               phiminusvals1,phiminusvals2,phiminusvals3,phiminusvals4,phiminusvals5,
               phiplussvals1,phiplussvals2,phiplussvals3,phiplussvals4,phiplussvals5,
               Magevec,
               curr_mort_pwid,curr_mort_former,
               death_rate_dc,death_rate_hcc,death_rate_lt1,death_rate_lt2)
  N=rep(0,720)
  # starting values
  N[1]=560;N[361]=435.6;N[367]=4.4
  
  init       = N
  parameters = param_vals
  times      = seq(0, 81, by = 1)
  lt=length(times)
  X = ode(y=init, times=times, func=MODEL, parms=parameters)
  X=t(X[,2:721])
  
  f1 = (sum(X[L$chronic10,66])/sum(X[L$chronic00,66])-40/144)^2
  f2 = (sum(X[L$chronic10,66])/sum(X[361:540,66])-0.5)^2
  # from matlab
  # f1 = (sum(X(chronic10,66))/sum(X(chronic00,66))-)^2;
  # f2 = (sum(X(chronic10,66))/sum(X(361:540,66))-0.5)^2;
  cat(paste0(f1,",",f2,"\n"))
  return(c(f1,f2))
}




MODEL <- function(time, state, parameters) {
  
  Mvec=parameters[1:400]
  M=t(matrix(Mvec,nrow=20,ncol=20))
  Mdashvec=parameters[401:800]
  Mdash=t(matrix(Mdashvec,nrow=20,ncol=20))
  extra_parms_vals=parameters[801:803]
  phiminusvals=parameters[804:808]
  phiplussvals=parameters[809:813]
  Magevec=parameters[814:894]
  age_matrix=t(matrix(Magevec,nrow=9,ncol=9))
  curr_mort_pwid=parameters[895:903]
  curr_mort_former=parameters[904:912]
  death_rate_dc=parameters[913]
  death_rate_hcc=parameters[914]
  death_rate_lt1=parameters[915]
  death_rate_lt2=parameters[916]
  
  mort_current = t(matrix(rep(curr_mort_pwid,20),9,20))
  mort_former = t(matrix(rep(curr_mort_former,20),9,20))
  
  Irows = c(6,7,8,9,10,11,12,13,14,15) # rows of infectious PWID
  death_vec= c(-log(1-death_rate_dc), -log(1-death_rate_hcc), -log(1-death_rate_lt1), -log(1- death_rate_lt2)); #DC,HCC,LT1,LT2
  t=time
  N=state
  # percentage PWID infected
  dum=matrix(rep(Irows,each=9),nrow=9)+matrix(seq(360,520,20),nrow=9,ncol=10)
  top10_index = matrix(t(dum),nrow=90,ncol=1)
  dum=matrix(rep(Irows,each=9),nrow=9)+matrix(seq(540,700,20),nrow=9,ncol=10)
  top11_index =  matrix(t(dum),nrow=90,ncol=1)
  top = sum(N[top10_index])+sum(N[top11_index]) # 40 is j = 0, 60 is j = 1 and i = 1 for both
  bot = sum(N[361:length(N)]) # 40 is j = 0, 60 is j = 1 and i = 1 for both
  I = top/bot;
  
  # incidence function , note put in =51
  if (t>=0 && t<=51){
    scaleI = 1+1.5*t/51;
  }else if (t>51 && t<=56){
    scaleI = 2.5-1.5*(t-51)/5;
  }else{
    scaleI=1;
  }
  endy=length(phiminusvals)
  phi = scaleI*extra_parms_vals[1]*I
  dum = Mdash[5,5]
  Mdash[phiminusvals[1:(endy-1)]]=-phi
  Mdash[phiminusvals[endy]]=-phi+dum; #dum = -r_svr4DC-rsvr4HCC
  Mdash[phiplussvals]=phi;
  
  X00 = matrix(N[1:(20*9)],20,9)
  X01 = matrix(N[181:(180+20*9)],20,9)
  X10 = matrix(N[361:(360+20*9)],20,9)
  X11 = matrix(N[541:(540+20*9)],20,9)
  
  # with deaths set as zero - this should be closed
  # age movement also set as zero
  # there are treatment so no failed treatment
  # so only D10 and D00
  # only for the youngest age group since no transistion
  d00=M%*%X00-extra_parms_vals[2]*X00+extra_parms_vals[3]*X10-mort_former*X00+t(age_matrix%*%t(X00))
  d01=M%*%X01-extra_parms_vals[2]*X01+extra_parms_vals[3]*X11-mort_former*X01+t(age_matrix%*%t(X01))
  d10=Mdash%*%X10+extra_parms_vals[2]*X00-extra_parms_vals[3]*X10-mort_current*X10+t(age_matrix%*%t(X10))
  d11=Mdash%*%X11+extra_parms_vals[2]*X01-extra_parms_vals[3]*X11-mort_current*X11+t(age_matrix%*%t(X11))
  d10[1,1] =d10[1,1] + sum(mort_former*X00) + sum(mort_former*X01) + sum(mort_current*X10) + sum(mort_current*X11) +
    death_vec[1]*sum(X00[12,])  + death_vec[2]*sum(X00[13,])  + death_vec[3]*sum(X00[14,])  + death_vec[4]*sum(X00[15,]) +
    death_vec[1]*sum(X01[12,])  + death_vec[2]*sum(X01[13,])  + death_vec[3]*sum(X01[14,])  + death_vec[4]*sum(X01[15,]) +
    death_vec[1]*sum(X10[12,])  + death_vec[2]*sum(X10[13,])  + death_vec[3]*sum(X10[14,])  + death_vec[4]*sum(X10[15,]) +
    death_vec[1]*sum(X11[12,])  + death_vec[2]*sum(X11[13,])  + death_vec[3]*sum(X11[14,])  + death_vec[4]*sum(X11[15,]) + 
    0*0.005*sum(state) # no growth 
  # note this is exponential growth exp(percentage*t)*P0
  # cat(paste0("sum state=",sum(state)," t=",t,"\n"))
  return(list(c(as.vector(d00), as.vector(d01), as.vector(d10), as.vector(d11))))
  
}

coeffsof_v5=function (){
  # compnames - root name sof the 20 compartments
  # chronic_nams - names of ij ordering
  # chronic_numsij - compartment index for all chronically infected with HCV
  
  
  # Gives the index numbers for categories
  # There are 20*9*4 compartments
  # 20 represents the X notation
  # element 1 to 5 S0,S1,S2,S3,S4
  # infection naïve or previously achieving spontaneous clearance or SVR through treatment from liver fibrosis stage F0 to F4 respectively
  # elemment 6 is A acute stage
  # element 7 to 11 F0,F1,F2,F3,F4 - the fibrosis stages, chronic
  # element 12 DC stage, chronic
  # element 13 HCC stage, chronic
  # element 14 and 15 LT1 and LT2 - liver transplant stages, chronic
  # element 16 to 20 T0 to T4 
  # chronically infected and in treatment achieving sustained viral response (SVR) (T0 to T4-treated from liver fibrosis stage F0 to F4 respectively)
  compnames=c('S0','S1','S2','S3','S4','A','F0','F1','F2','F3','F4','DC','HCC','LT1','LT2', 'T0','T1','T2','T3','T4')
  chronic_nams=c('00=former,never','01=former,failed','10=current,never','11=current,failed')
  
  # each elemet has three subscripts i,j,k
  # i = 0 formwerPWID, i = 1 current PWID
  # j = 0 never failed treatment, j = 1 had treatment failure
  # k = 1 to 9 the 9 age groups
  
  #for the 20*9*4 compartments
  # ordering is 
  # (i=0,j=0,k = 1);(i=0,j=0,k = 2);...;(i=0,j=0,k = 9) former, never failed
  # (i=0,j=1,k = 1);(i=0,j=0,k = 2);...;(i=0,j=0,k = 9) former, failed
  # (i=1,j=0,k = 1);(i=0,j=0,k = 2);...;(i=0,j=0,k = 9) current, never failed
  # (i=1,j=1,k = 1);(i=0,j=0,k = 2);...;(i=0,j=0,k = 9) current, failed
  
  chronic_numsA =7:15
  chronic_nums00 = chronic_numsA
  chronic_nums01 = 180+chronic_numsA
  chronic_nums10 = 360+chronic_numsA
  chronic_nums11 = 540+chronic_numsA
  
  for (i in 2 : 9){
    chronic_nums00 = c(chronic_nums00 ,20*(i-1)+chronic_numsA)
    chronic_nums01 = c(chronic_nums01 ,180+20*(i-1)+chronic_numsA)
    chronic_nums10 = c(chronic_nums10 ,360+20*(i-1)+chronic_numsA)
    chronic_nums11 = c(chronic_nums11 ,540+20*(i-1)+chronic_numsA)
  }
  
  A_comps_00=seq(6,180,20)
  A_comps_01=A_comps_00+180
  A_comps_10=A_comps_00+360
  A_comps_11=A_comps_00+540
  A_comps=c(A_comps_00,A_comps_01,A_comps_10,A_comps_11)
  
  F0_comps_00 = seq(7,180,20)
  F0_comps_01 = F0_comps_00+180
  F0_comps_10 = F0_comps_00+360
  F0_comps_11 = F0_comps_00+540
  F0_comps=c(F0_comps_00,F0_comps_01,F0_comps_10,F0_comps_11)
  
  F1_comps_00 = seq(8,180,20)
  F1_comps_01 = F1_comps_00+180
  F1_comps_10 = F1_comps_00+360
  F1_comps_11 = F1_comps_00+540
  F1_comps=c(F1_comps_00,F1_comps_01,F1_comps_10,F1_comps_11)
  
  F2_comps_00 = seq(9,180,20)
  F2_comps_01 = F2_comps_00+180
  F2_comps_10 = F2_comps_00+360
  F2_comps_11 = F2_comps_00+540
  F2_comps=c(F2_comps_00,F2_comps_01,F2_comps_10,F2_comps_11)
  
  F3_comps_00 = seq(10,180,20)
  F3_comps_01 = F3_comps_00+180
  F3_comps_10 = F3_comps_00+360
  F3_comps_11 = F3_comps_00+540
  F3_comps=c(F3_comps_00,F3_comps_01,F3_comps_10,F3_comps_11)
  
  F4_comps_00 = seq(11,180,20)
  F4_comps_01 = F4_comps_00+180
  F4_comps_10 = F4_comps_00+360
  F4_comps_11 = F4_comps_00+540
  F4_comps=c(F4_comps_00,F4_comps_01,F4_comps_10,F4_comps_11)
  
  DC_comps_00 = seq(12,180,20)
  DC_comps_01 = DC_comps_00+180
  DC_comps_10 = DC_comps_00+360
  DC_comps_11 = DC_comps_00+540
  DC_comps=c(DC_comps_00,DC_comps_01,DC_comps_10,DC_comps_11)
  
  HCC_comps_00 = seq(13,180,20)
  HCC_comps_01 = HCC_comps_00+180
  HCC_comps_10 = HCC_comps_00+360
  HCC_comps_11 = HCC_comps_00+540
  HCC_comps=c(HCC_comps_00,HCC_comps_01,HCC_comps_10,HCC_comps_11)
  
  LT1_comps_00 = seq(14,180,20)
  LT1_comps_01 = LT1_comps_00+180
  LT1_comps_10 = LT1_comps_00+360
  LT1_comps_11 = LT1_comps_00+540
  LT1_comps=c(LT1_comps_00,LT1_comps_01,LT1_comps_10,LT1_comps_11)
  
  LT2_comps_00 = seq(15,180,20)
  LT2_comps_01 = LT2_comps_00+180
  LT2_comps_10 = LT2_comps_00+360
  LT2_comps_11 = LT2_comps_00+540
  LT2_comps=c(LT2_comps_00,LT2_comps_01,LT2_comps_10,LT2_comps_11)
  
  
  S0_comps_00 = seq(1,180,20)
  S0_comps_01 = S0_comps_00+180
  S0_comps_10 = S0_comps_00+360
  S0_comps_11 = S0_comps_00+540
  S0_comps=c(S0_comps_00,S0_comps_01,S0_comps_10,S0_comps_11)
  
  S1_comps_00 = seq(2,180,20)
  S1_comps_01 = S1_comps_00+180
  S1_comps_10 = S1_comps_00+360
  S1_comps_11 = S1_comps_00+540
  S1_comps=c(S1_comps_00,S1_comps_01,S1_comps_10,S1_comps_11)
  
  S2_comps_00 = seq(3,180,20)
  S2_comps_01 = S2_comps_00+180
  S2_comps_10 = S2_comps_00+360
  S2_comps_11 = S2_comps_00+540
  S2_comps=c(S2_comps_00,S2_comps_01,S2_comps_10,S2_comps_11)
  
  S3_comps_00 = seq(4,180,20)
  S3_comps_01 = S3_comps_00+180
  S3_comps_10 = S3_comps_00+360
  S3_comps_11 = S3_comps_00+540
  S3_comps=c(S3_comps_00,S3_comps_01,S3_comps_10,S3_comps_11)
  
  S4_comps_00 = seq(5,180,20)
  S4_comps_01 = S4_comps_00+180
  S4_comps_10 = S4_comps_00+360
  S4_comps_11 = S4_comps_00+540
  S4_comps=c(S4_comps_00,S4_comps_01,S4_comps_10,S4_comps_11)
  
  T0_comps_00 = seq(16,180,20)
  T0_comps_01 = T0_comps_00+180
  T0_comps_10 = T0_comps_00+360
  T0_comps_11 = T0_comps_00+540
  T0_comps=c(T0_comps_00,T0_comps_01,T0_comps_10,T0_comps_11)
  
  T1_comps_00 = seq(17,180,20)
  T1_comps_01 = T1_comps_00+180
  T1_comps_10 = T1_comps_00+360
  T1_comps_11 = T1_comps_00+540
  T1_comps=c(T1_comps_00,T1_comps_01,T1_comps_10,T1_comps_11)
  
  T2_comps_00 = seq(18,180,20)
  T2_comps_01 = T2_comps_00+180
  T2_comps_10 = T2_comps_00+360
  T2_comps_11 = T2_comps_00+540
  T2_comps=c(T2_comps_00,T2_comps_01,T2_comps_10,T2_comps_11)
  
  T3_comps_00 = seq(19,180,20)
  T3_comps_01 = T3_comps_00+180
  T3_comps_10 = T3_comps_00+360
  T3_comps_11 = T3_comps_00+540
  T3_comps=c(T3_comps_00,T3_comps_01,T3_comps_10,T3_comps_11)
  
  T4_comps_00 = seq(20,180,20)
  T4_comps_01 = T4_comps_00+180
  T4_comps_10 = T4_comps_00+360
  T4_comps_11 = T4_comps_00+540
  T4_comps=c(T4_comps_00,T4_comps_01,T4_comps_10,T4_comps_11)
  # S0, A,F0,F1,F2,F3,F4,DC,HCC,LT1,LT2 - all other comps are empty
  agemat_current=rbind(  # current never failed
    c(361, 366:375),
    c(381 ,386:395),
    c(401 ,406:415),
    c(421 ,426:435),
    c(441 ,446:455),
    c(461 ,466:475),
    c(481 ,486:495),
    c(501 ,506:515),
    c(521 ,526:535))
  
  agemat_former=rbind( # former never failed
    c(361, 366:375),
    c(381 ,386:395),
    c(401 ,406:415),
    c(421 ,426:435),
    c(441 ,446:455),
    c(461 ,466:475),
    c(481 ,486:495),
    c(501 ,506:515),
    c(521 ,526:535))-360
  
  
  return(list(comps=compnames,states=chronic_nams,
              chronic00=chronic_nums00,chronic01=chronic_nums01,chronic10=chronic_nums10,chronic11=chronic_nums11,
              A_comps=A_comps,A_comps_00=A_comps_00,A_comps_01=A_comps_01,A_comps_10=A_comps_10,A_comps_11=A_comps_11,
              F0_comps=F0_comps,F0_comps_00=F0_comps_00,F0_comps_01=F0_comps_01,F0_comps_10=F0_comps_10,F0_comps_11=F0_comps_11,
              F1_comps=F1_comps,F1_comps_00=F1_comps_00,F1_comps_01=F1_comps_01,F1_comps_10=F1_comps_10,F1_comps_11=F1_comps_11,
              F2_comps=F2_comps,F2_comps_00=F2_comps_00,F2_comps_01=F2_comps_01,F2_comps_10=F2_comps_10,F2_comps_11=F2_comps_11,
              F3_comps=F3_comps,F3_comps_00=F3_comps_00,F3_comps_01=F3_comps_01,F3_comps_10=F3_comps_10,F3_comps_11=F3_comps_11,
              F4_comps=F4_comps,F4_comps_00=F4_comps_00,F4_comps_01=F4_comps_01,F4_comps_10=F4_comps_10,F4_comps_11=F4_comps_11,
              DC_comps=DC_comps,DC_comps_00=DC_comps_00,DC_comps_01=DC_comps_01,DC_comps_10=DC_comps_10,DC_comps_11=DC_comps_11,
              HCC_comps=HCC_comps,HCC_comps_00=HCC_comps_00,HCC_comps_01=HCC_comps_01,HCC_comps_10=HCC_comps_10,HCC_comps_11=HCC_comps_11,
              LT1_comps=LT1_comps,LT1_comps_00=LT1_comps_00,LT1_comps_01=LT1_comps_01,LT1_comps_10=LT1_comps_10,LT1_comps_11=LT1_comps_11,
              LT2_comps=LT2_comps,LT2_comps_00=LT2_comps_00,LT2_comps_01=LT2_comps_01,LT2_comps_10=LT2_comps_10,LT2_comps_11=LT2_comps_11,
              S0_comps=S0_comps,S0_comps_00=S0_comps_00,S0_comps_01=S0_comps_01,S0_comps_10=S0_comps_10,S0_comps_11=S0_comps_11,
              S1_comps=S1_comps,S1_comps_00=S1_comps_00,S1_comps_01=S1_comps_01,S1_comps_10=S1_comps_10,S1_comps_11=S1_comps_11,
              S2_comps=S2_comps,S2_comps_00=S2_comps_00,S2_comps_01=S2_comps_01,S2_comps_10=S2_comps_10,S2_comps_11=S2_comps_11,
              S3_comps=S3_comps,S3_comps_00=S3_comps_00,S3_comps_01=S3_comps_01,S3_comps_10=S3_comps_10,S3_comps_11=S3_comps_11,
              S4_comps=S4_comps,S4_comps_00=S4_comps_00,S4_comps_01=S4_comps_01,S4_comps_10=S4_comps_10,S4_comps_11=S4_comps_11,
              T0_comps=T0_comps,T0_comps_00=T0_comps_00,T0_comps_01=T0_comps_01,T0_comps_10=T0_comps_10,T0_comps_11=T0_comps_11,
              T1_comps=T1_comps,T1_comps_00=T1_comps_00,T1_comps_01=T1_comps_01,T1_comps_10=T1_comps_10,T1_comps_11=T1_comps_11,
              T2_comps=T2_comps,T2_comps_00=T2_comps_00,T2_comps_01=T2_comps_01,T2_comps_10=T2_comps_10,T2_comps_11=T2_comps_11,
              T3_comps=T3_comps,T3_comps_00=T3_comps_00,T3_comps_01=T3_comps_01,T3_comps_10=T3_comps_10,T3_comps_11=T3_comps_11,
              T4_comps=T4_comps,T4_comps_00=T4_comps_00,T4_comps_01=T4_comps_01,T4_comps_10=T4_comps_10,T4_comps_11=T4_comps_11,
              agemat_current=agemat_current,agemat_former=agemat_former))
}

age_weighted=function(X,agemat_current,agemat_former){
  # X is the compartments matrix 720 *times
  # agemat_current and agemat_former
  # these are 9 by 11 matrices
  # rows are age group compartment numbers
  # cols are S0,A,F0,F1,F2,F3,F4,DC,HCC,LT1,LT2
  # Assumes 25% of current are in first age group in year 66 (i.e. 2015 assuming 1950:2031 for model)
  # current
  bot = colSums(X[361:540,]); # currently infected
  age_weights_current=rep(9,1)
  for (i in 1 : 9){
    agetot=colSums(X[agemat_current[i,],])
    if (i == 1){
      age_weights_current[i] = 0.25*bot[66]/agetot[66]
    }else{
      age_weights_current[i] = (0.75/8)*bot[66]/agetot[66]
    }
  }
  # former
  bot = colSums(X[1:360,]); # currently infected
  age_weights_former=rep(9,1)
  for (i in 1 : 9){
    agetot=colSums(X[agemat_former[i,],])
    if (i == 1){
      age_weights_former[i] = 0.125*bot[66]/agetot[66]
    }else{
      age_weights_former[i] = (0.875/8)*bot[66]/agetot[66]
    }
  }
  return(list(age_weights_current=age_weights_current,age_weights_former=age_weights_former))
  
  
}
getallpop=function(X,L,age_weights_current,age_weights_former,nt,age_scale00,age_scale10){
  #X is the 720 * 82 matrix
  #L is the list of all the compartment index
  Aformer = matrix(rep(age_weights_former,nt),nrow=9,ncol=nt)
  Apwid = matrix(rep(age_weights_current,nt),nrow=9,ncol=nt)
  with(L,{
    lt2=age_scale00*colSums(Aformer*X[LT2_comps_00,])+age_scale10*colSums(Apwid*X[LT2_comps_10,])+
      age_scale00*colSums(Aformer*X[LT2_comps_01,])+age_scale10*colSums(Apwid*X[LT2_comps_11,])
    lt1=age_scale00*colSums(Aformer*X[LT1_comps_00,])+age_scale10*colSums(Apwid*X[LT1_comps_10,])+
      age_scale00*colSums(Aformer*X[LT1_comps_01,])+age_scale10*colSums(Apwid*X[LT1_comps_11,])
    lt=lt1+lt2;
    dc=age_scale00*colSums(Aformer*X[DC_comps_00,])+age_scale10*colSums(Apwid*X[DC_comps_10,])+
      age_scale00*colSums(Aformer*X[DC_comps_01,])+age_scale10*colSums(Apwid*X[DC_comps_11,])
    hcc=age_scale00*colSums(Aformer*X[HCC_comps_00,])+age_scale10*colSums(Apwid*X[HCC_comps_10,])+
      age_scale00*colSums(Aformer*X[HCC_comps_01,])+age_scale10*colSums(Apwid*X[HCC_comps_11,])
    f4=age_scale00*colSums(Aformer*X[F4_comps_00,])+age_scale10*colSums(Apwid*X[F4_comps_10,])+
      age_scale00*colSums(Aformer*X[F4_comps_01,])+age_scale10*colSums(Apwid*X[F4_comps_11,])
    f3=age_scale00*colSums(Aformer*X[F3_comps_00,])+age_scale10*colSums(Apwid*X[F3_comps_10,])+
      age_scale00*colSums(Aformer*X[F3_comps_01,])+age_scale10*colSums(Apwid*X[F3_comps_11,])
    f2=age_scale00*colSums(Aformer*X[F2_comps_00,])+age_scale10*colSums(Apwid*X[F2_comps_10,])+
      age_scale00*colSums(Aformer*X[F2_comps_01,])+age_scale10*colSums(Apwid*X[F2_comps_11,])
    f1=age_scale00*colSums(Aformer*X[F1_comps_00,])+age_scale10*colSums(Apwid*X[F1_comps_10,])+
      age_scale00*colSums(Aformer*X[F1_comps_01,])+age_scale10*colSums(Apwid*X[F1_comps_11,])
    f0=age_scale00*colSums(Aformer*X[F0_comps_00,])+age_scale10*colSums(Apwid*X[F0_comps_10,])+
      age_scale00*colSums(Aformer*X[F0_comps_01,])+age_scale10*colSums(Apwid*X[F0_comps_11,])
    
    # s compartments
    s4=age_scale00*colSums(Aformer*X[S4_comps_00,])+age_scale10*colSums(Apwid*X[S4_comps_10,])+
      age_scale00*colSums(Aformer*X[S4_comps_01,])+age_scale10*colSums(Apwid*X[S4_comps_11,])
    s3=age_scale00*colSums(Aformer*X[S3_comps_00,])+age_scale10*colSums(Apwid*X[S3_comps_10,])+
      age_scale00*colSums(Aformer*X[S3_comps_01,])+age_scale10*colSums(Apwid*X[S3_comps_11,])
    s2=age_scale00*colSums(Aformer*X[S2_comps_00,])+age_scale10*colSums(Apwid*X[S2_comps_10,])+
      age_scale00*colSums(Aformer*X[S2_comps_01,])+age_scale10*colSums(Apwid*X[S2_comps_11,])
    s1=age_scale00*colSums(Aformer*X[S1_comps_00,])+age_scale10*colSums(Apwid*X[S1_comps_10,])+
      age_scale00*colSums(Aformer*X[S1_comps_01,])+age_scale10*colSums(Apwid*X[S1_comps_11,])
    s0=age_scale00*colSums(Aformer*X[S0_comps_00,])+age_scale10*colSums(Apwid*X[S0_comps_10,])+
      age_scale00*colSums(Aformer*X[S0_comps_01,])+age_scale10*colSums(Apwid*X[S0_comps_11,])
    a0=age_scale00*colSums(Aformer*X[A_comps_00,])+age_scale10*colSums(Apwid*X[A_comps_10,])+
      age_scale00*colSums(Aformer*X[A_comps_01,])+age_scale10*colSums(Apwid*X[A_comps_11,])
    
    t4=age_scale00*colSums(Aformer*X[T4_comps_00,])+age_scale10*colSums(Apwid*X[T4_comps_10,])+
      age_scale00*colSums(Aformer*X[T4_comps_01,])+age_scale10*colSums(Apwid*X[T4_comps_11,])
    t3=age_scale00*colSums(Aformer*X[T3_comps_00,])+age_scale10*colSums(Apwid*X[T3_comps_10,])+
      age_scale00*colSums(Aformer*X[T3_comps_01,])+age_scale10*colSums(Apwid*X[T3_comps_11,])
    t2=age_scale00*colSums(Aformer*X[T2_comps_00,])+age_scale10*colSums(Apwid*X[T2_comps_10,])+
      age_scale00*colSums(Aformer*X[T2_comps_01,])+age_scale10*colSums(Apwid*X[T2_comps_11,])
    t1=age_scale00*colSums(Aformer*X[T1_comps_00,])+age_scale10*colSums(Apwid*X[T1_comps_10,])+
      age_scale00*colSums(Aformer*X[T1_comps_01,])+age_scale10*colSums(Apwid*X[T1_comps_11,])
    t0=age_scale00*colSums(Aformer*X[T0_comps_00,])+age_scale10*colSums(Apwid*X[T0_comps_10,])+
      age_scale00*colSums(Aformer*X[T0_comps_01,])+age_scale10*colSums(Apwid*X[T0_comps_11,])
    
    Tformer=age_scale00*(colSums(Aformer*X[LT2_comps_00,])+colSums(Aformer*X[LT2_comps_01,])+
                           colSums(Aformer*X[LT1_comps_00,])+colSums(Aformer*X[LT1_comps_01,])+
                           colSums(Aformer*X[DC_comps_00,]) +colSums(Aformer*X[DC_comps_01,]) + 
                           colSums(Aformer*X[HCC_comps_00,])+colSums(Aformer*X[HCC_comps_01,])+
                           colSums(Aformer*X[F4_comps_00,]) +colSums(Aformer*X[F4_comps_01,])+
                           colSums(Aformer*X[F3_comps_00,]) +colSums(Aformer*X[F3_comps_01,])+
                           colSums(Aformer*X[F2_comps_00,]) +colSums(Aformer*X[F2_comps_01,])+
                           colSums(Aformer*X[F1_comps_00,]) +colSums(Aformer*X[F1_comps_01,])+
                           colSums(Aformer*X[F0_comps_00,]) +colSums(Aformer*X[F0_comps_01,])+ 
                           colSums(Aformer*X[S4_comps_00,]) +colSums(Aformer*X[S4_comps_01,])+
                           colSums(Aformer*X[S3_comps_00,]) +colSums(Aformer*X[S3_comps_01,])+
                           colSums(Aformer*X[S2_comps_00,]) +colSums(Aformer*X[S2_comps_01,])+
                           colSums(Aformer*X[S1_comps_00,]) +colSums(Aformer*X[S1_comps_01,])+
                           colSums(Aformer*X[S0_comps_00,]) +colSums(Aformer*X[S0_comps_01,])+ 
                           colSums(Aformer*X[T4_comps_00,]) +colSums(Aformer*X[T4_comps_01,])+
                           colSums(Aformer*X[T3_comps_00,]) +colSums(Aformer*X[T3_comps_01,])+
                           colSums(Aformer*X[T2_comps_00,]) +colSums(Aformer*X[T2_comps_01,])+
                           colSums(Aformer*X[T1_comps_00,]) +colSums(Aformer*X[T1_comps_01,])+
                           colSums(Aformer*X[T0_comps_00,]) +colSums(Aformer*X[T0_comps_01,])+  
                           colSums(Aformer*X[A_comps_00,]) +colSums(Aformer*X[A_comps_01,])                         
    )   
    
    Tcurrent=age_scale10*(colSums(Apwid*X[LT2_comps_10,])+colSums(Apwid*X[LT2_comps_11,])+
                            colSums(Apwid*X[LT1_comps_10,])+colSums(Apwid*X[LT1_comps_11,])+
                            colSums(Apwid*X[DC_comps_10,]) +colSums(Apwid*X[DC_comps_11,]) + 
                            colSums(Apwid*X[HCC_comps_10,])+colSums(Apwid*X[HCC_comps_11,])+
                            colSums(Apwid*X[F4_comps_10,]) +colSums(Apwid*X[F4_comps_11,])+
                            colSums(Apwid*X[F3_comps_10,]) +colSums(Apwid*X[F3_comps_11,])+
                            colSums(Apwid*X[F2_comps_10,]) +colSums(Apwid*X[F2_comps_11,])+
                            colSums(Apwid*X[F1_comps_10,]) +colSums(Apwid*X[F1_comps_11,])+
                            colSums(Apwid*X[F0_comps_10,]) +colSums(Apwid*X[F0_comps_11,])+ 
                            colSums(Apwid*X[S4_comps_10,]) +colSums(Apwid*X[S4_comps_11,])+
                            colSums(Apwid*X[S3_comps_10,]) +colSums(Apwid*X[S3_comps_11,])+
                            colSums(Apwid*X[S2_comps_10,]) +colSums(Apwid*X[S2_comps_11,])+
                            colSums(Apwid*X[S1_comps_10,]) +colSums(Apwid*X[S1_comps_11,])+
                            colSums(Apwid*X[S0_comps_10,]) +colSums(Apwid*X[S0_comps_11,])+ 
                            colSums(Apwid*X[T4_comps_10,]) +colSums(Apwid*X[T4_comps_11,])+
                            colSums(Apwid*X[T3_comps_10,]) +colSums(Apwid*X[T3_comps_11,])+
                            colSums(Apwid*X[T2_comps_10,]) +colSums(Apwid*X[T2_comps_11,])+
                            colSums(Apwid*X[T1_comps_10,]) +colSums(Apwid*X[T1_comps_11,])+
                            colSums(Apwid*X[T0_comps_10,]) +colSums(Apwid*X[T0_comps_11,])+  
                            colSums(Apwid*X[A_comps_10,]) +colSums(Apwid*X[A_comps_11,])                         
    )   
    
    
    
    return(list(lt=lt,dc=dc,hcc=hcc,s0=s0,s1=s1,s2=s2,s3=s3,s4=s4,a0=a0,f0=f0,f1=f1,f2=f2,f3=f3,f4=f4,
                t0=t0,t1=t1,t2=t2,t3=t3,t4=t4,Tformer=Tformer,Tcurrent=Tcurrent))
  })
  
}

plot_allcomps=function(times,C){
  y=C$lt+C$dc+C$hcc+C$s0+C$s1+C$s2+C$s3+C$s4+C$a0+C$f0+C$f1+C$f2+C$f3+C$f4+C$t0+C$t1+C$t2+C$t3+C$t4
  df=data.frame(times,y)
  p1=ggplot(data=subset(df,times>=66), aes(x=times+1950, y=y)) + geom_line() +
    xlab("Year") + ylab("Population") + ggtitle("Base model - Number in all compartments")+
    scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))+
    scale_x_continuous(breaks = c(seq(2015,2031,2)))+
    theme(plot.background = element_rect(fill = "#BFD5E3"),axis.line = element_line(colour = "black"))
  #scale_x_discrete(labels=c("2015","2017","2019","2021","2023","2025","2027","2029","2031"))
  return(p1)  
}


plot_areacomps=function(times,C){
  y=c(C$f0,C$f1,C$f2,C$f3,C$f4,C$dc,C$hcc,C$lt)
  Stage=c(rep("F0",82),rep("F1",82),rep("F2",82),rep("F3",82),rep("F4",82),rep("DC",82),rep("HCC",82),rep("LT",82))
  t=rep(0:81,8)
  df=data.frame(y,t,Stage)
  df$Stage=factor(df$Stage,levels=c("F0","F1","F2","F3","F4","DC","HCC","LT"))
  p1=ggplot(data=subset(df,t>=66), aes(x=t+1950, y=y,fill=Stage)) + geom_area() +
    xlab("Year") + ylab("Population") + ggtitle("Base model - Number in chronic compartments")+
    scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))+
    scale_x_continuous(breaks = c(seq(2015,2031,2)))+
    theme(plot.background = element_rect(fill = "#BFD5E3"),axis.line = element_line(colour = "black"))
  #scale_x_discrete(labels=c("2015","2017","2019","2021","2023","2025","2027","2029","2031"))
  return(p1)  
}