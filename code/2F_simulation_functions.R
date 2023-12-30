###############################################################################
#  X: the independent variable
#  M: the measurement of the mediator
#  Dm: the manipulation of the mediator
#  Y: the dependent variable
#  V: the confounder representing the potentially omitted mediators

gen_data <- function(N,a,b,cp,d,g,h,NT,p=.5){
  # preset variances
  v.X = v.Dm = p*(1-p); v.M = v.Y = 1
  
  # computing variances of error terms and interaction terms
  v.eM = v.M - a^2*v.X - d^2*v.Dm - h^2*v.X*v.Dm
  
  v.XM = d^2*v.X*v.Dm + h^2*v.X^2*v.Dm + v.X*v.eM
  
  cov.XM = a^2*v.X
  
  v.eY = v.Y - cp^2*v.X - b^2*v.M - g^2*v.XM - 2*cp*b*cov.XM 
  
  # covariance matrix of the error terms
  Ve = matrix(c(v.eM, 0,
                0, v.eY), 2,2,byrow = T); colnames(Ve)=rownames(Ve) = c('eM','eY')
  
  errors = mvrnorm(N, c(0,0), Ve)
  
  if(!NT){
    errors = apply(errors, 2, tras_ADF, skew=-2, kurtosis = 6)
  }
  #X = rbinom(N, 1, 0.5); Dm = rbinom(N, 1, 0.5)
  X = rbinom(N, 1, 0.5)-0.5; Dm = rbinom(N, 1, 0.5)-0.5
  
  M = a*X+d*Dm+h*X*Dm+errors[,'eM']
  Y = cp*X+b*M+g*X*M+errors[,'eY']
  
  dat = data.frame(X=X, Dm=Dm, M=M, Y=Y)  
  
  return(dat)
}

fn = function(pars, skew, kurtosis){
  b = pars[1]; c = pars[2]; d = pars[3]
  eq1 = b^2 + 6*b*d + 2*c^2 + 15*d^2 -1
  eq2 = 2*c* (b^2 + 24*b*d + 105*d^2 + 2) - skew
  eq3 = 24*(b*d + c^2*(1+b^2+28*b*d) + d^2*(12 + 48*b*d + 141*c^2 + 225*d^2)) - kurtosis
  v = c(eq1, eq2, eq3)
  Q = sum(v^2)
  return(Q)
}

tras_ADF <- function(error, skew, kurtosis){
  bcd = nlm(fn, c(1,1,1), skew, kurtosis)$estimate
  b = bcd[1]; c=bcd[2]; d=bcd[3]; a=-c
  e.ADF  = a + b*error+ c*error^2 + d*error^3
  return(e.ADF)
}


# proposed method
IVSEM_main <- function(data, g, h, colnum){
  model = 'M~a*X+d*Dm+h*X:Dm
         Y~cp*X+b*M+g*X:M'
  fit = try(sem(model=model, data=data), silent=T)
  if(inherits(fit,"try-error")){
    res = rep(NA,colnum)
  }else{
    if(fit@optim[["converged"]]==F){
      res = rep(NA,colnum)
    }else{
      para.est = c(coef(fit)['a'], coef(fit)['b'], coef(fit)['cp'],
                   coef(fit)['d'],coef(fit)['g'],coef(fit)['h'])
      sum = summary(fit)
      pValues = sum$pe$pvalue[1:6]
      resV.M = sum$pe$est[7]; resV.Y = sum$pe$est[8]
      covs = vcov(fit)
      if(det(covs)<0){
        sigma.est=rep(NA,6);infoDef=-1;
        UPNIE.SEest = APNIE.SEest = UTNIE.SEest = ATNIE.SEest = NA
      }else{
        sigma.est = c(covs['a','a'],covs['b','b'],covs['cp','cp'],covs['d','d'],
                      covs['g','g'],covs['h','h'])
        
        IE.vcov = matrix(c(covs['a','a'], covs['a','b'], covs['a','g'], covs['a','h'],
                           covs['a','b'], covs['b','b'], covs['b','g'], covs['b','h'],
                           covs['a','g'], covs['b','g'], covs['g','g'], covs['g','h'],
                           covs['a','h'], covs['b','h'], covs['g','h'], covs['h','h']),4,4)
        
        delta.f.XnDn = matrix(c((coef(fit)['b']-0.5*coef(fit)['g']), (coef(fit)['a']-0.5*coef(fit)['h']),
                                (-.5*coef(fit)['a']+.25*coef(fit)['h']), (-.5*coef(fit)['b']+.25*coef(fit)['g'])), 4, 1)
        delta.f.XnDp = matrix(c((coef(fit)['b']-0.5*coef(fit)['g']), (coef(fit)['a']+0.5*coef(fit)['h']),
                                (-.5*coef(fit)['a']-.25*coef(fit)['h']), (.5*coef(fit)['b']-.25*coef(fit)['g'])), 4, 1)
        
        delta.f.XpDn = matrix(c((coef(fit)['b']+0.5*coef(fit)['g']), (coef(fit)['a']-0.5*coef(fit)['h']),
                                (.5*coef(fit)['a']-.25*coef(fit)['h']), (-.5*coef(fit)['b']-.25*coef(fit)['g'])), 4, 1)
        delta.f.XpDp = matrix(c((coef(fit)['b']+0.5*coef(fit)['g']), (coef(fit)['a']+0.5*coef(fit)['h']),
                                (.5*coef(fit)['a']+.25*coef(fit)['h']), (.5*coef(fit)['b']+.25*coef(fit)['g'])), 4, 1)
        
        UPNIE.SEest = sqrt(t(delta.f.XnDn) %*% IE.vcov %*% delta.f.XnDn)
        APNIE.SEest = sqrt(t(delta.f.XnDp) %*% IE.vcov %*% delta.f.XnDp)
        UTNIE.SEest = sqrt(t(delta.f.XpDn) %*% IE.vcov %*% delta.f.XpDn)
        ATNIE.SEest = sqrt(t(delta.f.XpDp) %*% IE.vcov %*% delta.f.XpDp)
        
        infoDef=1
      }
      res = c(para.est, sigma.est, pValues,UPNIE.SEest, APNIE.SEest,
              UTNIE.SEest, ATNIE.SEest, resV.M, resV.Y, infoDef)
    }
  }
  return(res)
}

TradMed <- function(data, g){
  models = list(m1 = lm(M~X*Dm,data), m2 = lm(Y~X*M,data))
  fit = try(mediate(models$m1,models$m2,treat='X',mediator='M', sims=1000, treat.value = 0.5, control.value = -0.5), silent=T)
  if (inherits(fit, 'try-error')){
    res = rep(NA,8)
  }else{
    sum=summary(fit)
    IE = c(sum$d0, sum$d1)
    IE.CI = c(sum$d0.ci, sum$d1.ci)
    IE.p = c(sum$d0.p, sum$d1.p)
    res = c(IE, IE.CI, IE.p)
  }
  return(res)
}

AnaRES <- function(nrep, N,a,b,cp,d,g,h,NT=TRUE,skew=0,kurtosis=0,p=.5){
  res.IVSEM = matrix(NA,nrep,25); colnames(res.IVSEM) <- c('a.est','b.est','cp.est','d.est','g.est','h.est',
                                                             'a.sigma','b.sigma','cp.sigma','d.sigma','g.sigma','h.sigma',
                                                             'a.p','d.p','h.p','cp.p','b.p','g.p',
                                                             'UPNIE.SEest','APNIE.SEest','UTNIE.SEest','ATNIE.SEest','resV.M','resV.Y','InfoDef')
  res.TradMed = matrix(NA,nrep,8); colnames(res.TradMed) <- c('PNIE','TNIE','PNIE.lower','PNIE.upper',
                                                                'TNIE.lower','TNIE.upper','PNIE.p','TNIE.p')


  
  sd.X = sd.Dm = sqrt(p*(1-p))
  a.u = a/sd.X;  cp.u = cp/sd.X; d.u = d/sd.Dm; g.u = g/sd.X; h.u = h/(sd.X*sd.Dm)
  
  for(i in 1:nrep){
    data = try(gen_data(N,a.u,b,cp.u,d.u,g.u,h.u,NT,p=.5))
    # the proposed modeling approach
    res.IVSEM[i,] = IVSEM_main(data, g, h, colnum.IVSEM)
    
    # Imai2013
    res.TradMed[i,] = TradMed(data, g.u)
  }
  res = list(IVSEM=res.IVSEM, TradMed=res.TradMed)
  return(res)
}

## performance measures
# estimation bias
cal.eBias <- function(theta.hat, theta){
  if(theta==0){
    Ebias <- mean(theta.hat, na.rm=T)-theta
  }else{Ebias <- (mean(theta.hat, na.rm=T)-theta)/theta}
  return(Ebias)
}

# coverage rate
cal.CR <- function(SE, theta, theta.hat){
  cv.error = qnorm(1-0.05/2) * SE
  l = theta.hat-cv.error  # lower bound
  u = theta.hat+cv.error  # upper bound
  lcriterion = (theta >= l)
  ucriterion = (theta <= u)
  CR = mean(lcriterion*ucriterion, na.rm=T)
  return(CR)
}
# rejection rate
cal.RR <- function(Est, seEst, para){
  Z.est = abs(Est)/seEst
  Z.cri = qnorm(1-.025)
  RR = mean(Z.est>=Z.cri, na.rm=T)
  
  if(para == 0){
    power = NA
    typeIerror = RR
  }else{
    power = RR
    typeIerror = NA
  }
  return(c(power,typeIerror))
}

sep.RR <- function(RR, para){
  if(para == 0){
    power = NA
    typeIerror = RR
  }else{
    power = RR
    typeIerror = NA
  }
  return(c(power,typeIerror))
}

# deleting results with erroneous variance estimation 
del_vcov <- function(resi, g, h){
  if(g==0&h==0){
    range = 5:8
  }
  if((g!=0&h==0)|(g==0&h!=0)){
    range = 6:10
  }
  if(g!=0&h!=0){
    range = 7:12
  }
  for(i in 1:nrow(resi)){
    if(is.na(resi[i,'resV.M'])){
      resi[i,] = rep(NA, ncol(resi))
    }else if(resi[i,'resV.M']<0|resi[i,'resV.Y']<0){
      resi[i,] = rep(NA, ncol(resi))
    }
  }
  resi = remove_empty(resi, which='rows')
  for(i in 1:nrow(resi)){
    for(j in range){
      if(is.na(resi[i,j])){
        resi[i,]=rep(NA, ncol(resi))
      }
    }
  }
  resi = remove_empty(resi, which='rows')
  return(resi)
}

IVSEM.main.reg <- function(res, a, b, cp, d, g, h){
  
  ConvRate = 1-mean(is.na(res[,'a.est']))
  res = remove_empty(res, which='rows')
  
  res = del_vcov(res, g, h)
  infoDef = mean(res[,'InfoDef'], na.rm = T)
  
  ## performance measures for the proposed method
  # point Estimate Bias
  a.ebias = cal.eBias(res[,'a.est'], a)
  b.ebias = cal.eBias(res[,'b.est'], b)
  cp.ebias = cal.eBias(res[,'cp.est'], cp)
  d.ebias = cal.eBias(res[,'d.est'], d)
  g.ebias = cal.eBias(res[,'g.est'], g)
  h.ebias = cal.eBias(res[,'h.est'], h)

  # se para
  a.se = sd(res[,'a.est'])
  b.se = sd(res[,'b.est'])
  cp.se = sd(res[,'cp.est'])
  d.se = sd(res[,'d.est'])
  g.se = sd(res[,'g.est'])
  h.se = sd(res[,'h.est'])

  # se bias
  a.sebias = cal.eBias(sqrt(res[,'a.sigma']), a.se)
  b.sebias = cal.eBias(sqrt(res[,'b.sigma']), b.se)
  cp.sebias = cal.eBias(sqrt(res[,'cp.sigma']), cp.se)
  d.sebias = cal.eBias(sqrt(res[,'d.sigma']), d.se)
  g.sebias = cal.eBias(sqrt(res[,'g.sigma']), g.se)
  h.sebias = cal.eBias(sqrt(res[,'h.sigma']), h.se)

  # coverage rate
  a.CR = cal.CR(sqrt(res[,'a.sigma']), a, res[,'a.est'])
  b.CR = cal.CR(sqrt(res[,'b.sigma']), b, res[,'b.est'])
  cp.CR = cal.CR(sqrt(res[,'cp.sigma']), cp, res[,'cp.est'])
  d.CR = cal.CR(sqrt(res[,'d.sigma']), d, res[,'d.est'])
  g.CR = cal.CR(sqrt(res[,'g.sigma']), g, res[,'g.est'])
  h.CR = cal.CR(sqrt(res[,'h.sigma']), h, res[,'h.est'])

  # rejection rate
  a.RR = sep.RR(mean(res[,'a.p']<.05), a)
  b.RR = sep.RR(mean(res[,'b.p']<.05), b)
  cp.RR = sep.RR(mean(res[,'cp.p']<.05), cp)
  d.RR = sep.RR(mean(res[,'d.p']<.05), d)
  g.RR = sep.RR(mean(res[,'g.p']<.05), g)
  h.RR = sep.RR(mean(res[,'h.p']<.05), h)
  
  # combining results for regression coefficients
  reg.res = matrix(c(a.ebias, b.ebias, cp.ebias, d.ebias, g.ebias, h.ebias,
                     a.se, b.se, cp.se, d.se, g.se, h.se,
                     a.sebias,b.sebias,cp.sebias,d.sebias,g.sebias,h.sebias,
                     a.CR, b.CR, cp.CR, d.CR, g.CR, h.CR,
                     a.RR, b.RR, cp.RR, d.RR, g.RR, h.RR,ConvRate,infoDef), 1,38)
  
  colnames(reg.res) = c('a.ebias','b.ebias','cp.ebias','d.ebias','g.ebias','h.ebias',
                        'a.empSE','b.empSE','cp.empSE','d.empSE','g.empSE','h.empSE',
                        'a.sebias','b.sebias','cp.sebias','d.sebias','g.sebias','h.sebias',
                        'a.CR','b.CR','cp.CR','d.CR','g.CR','h.CR',
                        'a.power','a.typeIerror','b.power','b.typeIerror',
                        'cp.power','cp.typeIerror','d.power','d.typeIerror',
                        'g.power','g.typeIerror','h.power','h.typeIerror','ConvRate','NegDef')
  
  return(reg.res)
}

IVSEM.main.IE <- function(res, a, b, g, h){
  ConvRate = 1-mean(is.na(res[,'a.est']))
  res = remove_empty(res, which='rows')
  
  res = del_vcov(res, g, h)
  infoDef = mean(res[,'InfoDef'], na.rm = T)
  
  UPNIE.tv = (a-0.5*h)*(b-0.5*g)
  APNIE.tv = (a+0.5*h)*(b-0.5*g)
  UTNIE.tv = (a-0.5*h)*(b+0.5*g)
  ATNIE.tv = (a+0.5*h)*(b+0.5*g)
  
  UPNIE.est = (res[,'a.est']-0.5*res[,'h.est'])*(res[,'b.est']-0.5*res[,'g.est'])
  APNIE.est = (res[,'a.est']+0.5*res[,'h.est'])*(res[,'b.est']-0.5*res[,'g.est'])
  UTNIE.est = (res[,'a.est']-0.5*res[,'h.est'])*(res[,'b.est']+0.5*res[,'g.est'])
  ATNIE.est = (res[,'a.est']+0.5*res[,'h.est'])*(res[,'b.est']+0.5*res[,'g.est'])
  
  UPNIE.ebias = cal.eBias(UPNIE.est, UPNIE.tv)
  APNIE.ebias = cal.eBias(APNIE.est, APNIE.tv)
  UTNIE.ebias = cal.eBias(UTNIE.est, UTNIE.tv)
  ATNIE.ebias = cal.eBias(ATNIE.est, ATNIE.tv)
  
  UPNIE.empSE = sd(UPNIE.est)
  APNIE.empSE = sd(APNIE.est)
  UTNIE.empSE = sd(UTNIE.est)
  ATNIE.empSE = sd(ATNIE.est)
  
  UPNIE.sebias <- cal.eBias(res[,'UPNIE.SEest'], UPNIE.empSE)
  APNIE.sebias <- cal.eBias(res[,'APNIE.SEest'], APNIE.empSE)
  UTNIE.sebias <- cal.eBias(res[,'UTNIE.SEest'], UTNIE.empSE)
  ATNIE.sebias <- cal.eBias(res[,'ATNIE.SEest'], ATNIE.empSE)
  
  UPNIE.CR <- cal.CR(res[,'UPNIE.SEest'], UPNIE.tv, UPNIE.est)
  APNIE.CR <- cal.CR(res[,'APNIE.SEest'], APNIE.tv, APNIE.est)
  UTNIE.CR <- cal.CR(res[,'UTNIE.SEest'], UTNIE.tv, UTNIE.est)
  ATNIE.CR <- cal.CR(res[,'ATNIE.SEest'], ATNIE.tv, ATNIE.est)
  
  UPNIE.RR <- cal.RR(UPNIE.est, res[,'UPNIE.SEest'], UPNIE.tv)
  APNIE.RR <- cal.RR(APNIE.est, res[,'APNIE.SEest'], APNIE.tv)
  UTNIE.RR <- cal.RR(UTNIE.est, res[,'UTNIE.SEest'], UTNIE.tv)
  ATNIE.RR <- cal.RR(ATNIE.est, res[,'ATNIE.SEest'], ATNIE.tv)
  
  LMIE.res = matrix(c(UPNIE.tv, APNIE.tv, UTNIE.tv, ATNIE.tv, 
                      UPNIE.ebias, APNIE.ebias, UTNIE.ebias, ATNIE.ebias, 
                      UPNIE.empSE, APNIE.empSE, UTNIE.empSE, ATNIE.empSE, 
                      UPNIE.sebias, APNIE.sebias, UTNIE.sebias, ATNIE.sebias,
                      UPNIE.CR, APNIE.CR, UTNIE.CR, ATNIE.CR, 
                      UPNIE.RR, APNIE.RR, UTNIE.RR, ATNIE.RR),1,28)
  
  colnames(LMIE.res) = c('UPNIE.tv','APNIE.tv','UTNIE.tv','ATNIE.tv',
                         'UPNIE.ebias','APNIE.ebias','UTNIE.ebias','ATNIE.ebias',
                         'UPNIE.empSE','APNIE.empSE','UTNIE.empSE','ATNIE.empSE',
                         'UPNIE.sebias','APNIE.sebias','UTNIE.sebias','ATNIE.sebias',
                         'UPNIE.CR','APNIE.CR','UTNIE.CR','ATNIE.CR',
                         'UPNIE.power','UPNIE.typeIerror','APNIE.power','APNIE.typeIerror',
                         'UTNIE.power','UTNIE.typeIerror','ATNIE.power','ATNIE.typeIerror')
  
  return(LMIE.res)
}

TradMed.res.IE <- function(res, a,b,cp,d,g,h){
  ConvRate = 1-mean(is.na(res[,1]))
  res = remove_empty(res, which='rows')
  
  UPNIE.tv = (a-0.5*h)*(b-0.5*g)
  APNIE.tv = (a+0.5*h)*(b-0.5*g)
  UTNIE.tv = (a-0.5*h)*(b+0.5*g)
  ATNIE.tv = (a+0.5*h)*(b+0.5*g) 
  
  TM.UPNIE.ebias = cal.eBias(res[,'PNIE'], UPNIE.tv)
  TM.APNIE.ebias = cal.eBias(res[,'PNIE'], APNIE.tv)
  TM.UTNIE.ebias = cal.eBias(res[,'TNIE'], UTNIE.tv)
  TM.ATNIE.ebias = cal.eBias(res[,'TNIE'], ATNIE.tv)
  
  TM.UPNIE.CR = mean((res[,'PNIE.lower']<=UPNIE.tv)*(res[,'PNIE.upper']>=UPNIE.tv), na.rm=T)
  TM.APNIE.CR = mean((res[,'PNIE.lower']<=APNIE.tv)*(res[,'PNIE.upper']>=APNIE.tv), na.rm=T)
  TM.UTNIE.CR = mean((res[,'TNIE.lower']<=UTNIE.tv)*(res[,'TNIE.upper']>=UTNIE.tv), na.rm=T)
  TM.ATNIE.CR = mean((res[,'TNIE.lower']<=ATNIE.tv)*(res[,'TNIE.upper']>=ATNIE.tv), na.rm=T)
  
  TM.UPNIE.RR = sep.RR(mean(res[,'PNIE.p']<=0.05, na.rm=T), UPNIE.tv)
  TM.APNIE.RR = sep.RR(mean(res[,'PNIE.p']<=0.05, na.rm=T), APNIE.tv)
  TM.UTNIE.RR = sep.RR(mean(res[,'TNIE.p']<=0.05, na.rm=T), UTNIE.tv)
  TM.ATNIE.RR = sep.RR(mean(res[,'TNIE.p']<=0.05, na.rm=T), ATNIE.tv)
  
  TradMed.res <- matrix(c(UPNIE.tv, APNIE.tv, UTNIE.tv, ATNIE.tv,
                          TM.UPNIE.ebias, TM.APNIE.ebias, TM.UTNIE.ebias, TM.ATNIE.ebias,
                          TM.UPNIE.CR, TM.APNIE.CR, TM.UTNIE.CR, TM.ATNIE.CR,
                          TM.UPNIE.RR, TM.APNIE.RR, TM.UTNIE.RR, TM.ATNIE.RR,ConvRate), 1, 21)
  colnames(TradMed.res) <- c('UPNIE.tv','APNIE.tv','UTNIE.tv','ATNIE.tv',
                             'UPNIE.ebias','APNIE.ebias','UTNIE.ebias','ATNIE.ebias',
                             'UPNIE.CR','APNIE.CR','UTNIE.CR','ATNIE.CR',
                             'UPNIE.power','UPNIE.typeIerror','APNIE.power','APNIE.typeIerror',
                             'UTNIE.power','UTNIE.typeIerror','ATNIE.power','ATNIE.typeIerror','ConvRate')
  
  return(TradMed.res)
}


simRES <- function(res, a, b, cp, d, g, h){
  IVSEM.res = IVSEM.main.reg(res[[1]], a, b, cp, d, g, h)
  IVSEM.IE.res = IVSEM.main.IE(res[[1]], a, b, g, h)
  
  TradMed.res = TradMed.res.IE(res[[2]],a,b,cp,d,g,h)
  
  final_res = list(IVSEM.res = IVSEM.res, IVSEM.IE.res = IVSEM.IE.res, 
                   TradMed.res=TradMed.res)
  return(final_res)
}


