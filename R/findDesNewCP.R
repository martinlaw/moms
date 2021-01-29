createCovMat <- function(J.,
                         K.,
                         rho.vec.){
  stage.row <- matrix(rep(1:J., each=K.), J.*K., J.*K.)
  stage.col <- t(stage.row)
  Lambda <- sqrt(pmin(stage.row, stage.col)/pmax(stage.row, stage.col))
  rho_submatrix <- matrix(1, K., K.)
  rho_submatrix[which(lower.tri(rho_submatrix))] <- rho.vec.
  rho_submatrix <- t(rho_submatrix)
  rho_submatrix[which(lower.tri(rho_submatrix))] <- rho.vec.
  rho_matrix <- matrix(NA, J.*K., J.*K.)
  for(j1 in 1:J.){
    for(j2 in 1:J.){
      rho_matrix[(1+(j1-1)*K.):(j1*K.), (1+(j2-1)*K.):(j2*K.)] <- rho_submatrix
    }
  }
  Lambda <- Lambda*rho_matrix
  return(Lambda)
}



findR <- function(bounds,
                  typeI.power="typeI",
                  ts,
                  n.stage,
                  return.optimisation=FALSE,
                  drop.outcomes=TRUE,
                  alpha.combine=TRUE,
                  nsims.,
                  K.,
                  m.,
                  max.outs.s2.,
                  working.outs.,
                  vars.,
                  delta0.,
                  delta1.,
                  delta.true.=NULL,
                  alpha.k.,
                  cp.l.,
                  cp.u.
                  ){
  if(length(bounds)==1){
    bounds <- rep(bounds, K.)
  }
  J <- 2
  denom <-  vars.
  numer <- rep(1*n.stage, each=K.)
  information.j <- numer/denom
  numer.J <- rep(J*n.stage, each=K.)
  information.final <- numer.J/denom
  if(typeI.power=="power"){
    lfc.effects <- delta0.
    lfc.effects[working.outs.] <- delta1.[working.outs.]
    tau <- rep(lfc.effects, times=J) * c(sqrt(information.j), sqrt(information.final))
    # Note that here, the LFC is used, but for calculating CP, all deltas are taken to be equal to their delta1
    ts <- sweep(ts, 2,  tau, "+") # TS for obtaining power
  }
  if(typeI.power=="truedelta"){
 #   browser()
    tau <- rep(delta.true., times=J) * c(sqrt(information.j), sqrt(information.final))
    # Note that here, the LFC is used, but for calculating CP, all deltas are taken to be equal to their delta1
    ts <- sweep(ts, 2,  tau, "+") # TS for obtaining power
  }
  ts.s1 <- ts[,1:K.]
  ts.s2 <- ts[,(K.+1):(2*K.)]
  # An important change here is that the number of outcomes dropped/retained is no longer fixed, but depends on the CP bounds:
  # drop all outcomes that fall below cp.l at the interim. Futility stopping remains a possibility if all outcomes drop below cp.l.
  # Calculate CP:
  z.alpha <- bounds
  cp.component1 <- sweep(ts.s1, 2, sqrt(information.j), "*") 
  cp.component23 <- -z.alpha*sqrt(information.final) + (information.final-information.j)*delta1. # Note: always delta1 here, whether typeIerror or power
  cp.numer <- sweep(cp.component1, 2, cp.component23, "+") # Add components 1 and 23 to create numerator
  cp.denom <- sqrt(information.final-information.j)
  cp <- pnorm(sweep(cp.numer, 2, cp.denom, "/"))
  #stop.futility <- apply(cp, 1, function(x) all(x<cp.l.)) # Too slow. Use line below
  below.cp.bounds <- cp<cp.l.
  stop.futility <- rowSums(below.cp.bounds)>=(K.-m.+1) # Stop for futility if K-m+1 outcomes are below CP_L at the interim.
  #stop.efficacy <- rep(FALSE, nsims.) # If no efficacy stopping permitted.
  #stop.efficacy <- apply(cp, 1, function(x) any(x>cp.u)) # only correct if m==1.
  stop.efficacy <- rowSums(cp>=cp.u.)>=m.
  # XXX IMPORTANT QUESTION: what is the purpose of cp.u if the null is only rejected at the end?
  # If there is interim rejection, this must be based on cp.u, as there are no interim boundaries.
  stop.any <- stop.futility | stop.efficacy
  continue <- !stop.any # Trials that carry on to full N.
  exceed.cpl <- !below.cp.bounds  # Binary: 1: CP>=CP_L
  
  
  ### minus.cp <- -cp
  ### greatest.cps <- t(apply(minus.cp, 1, function(x) rank(x) <= max.outs.s2.))
  ### retained.outcomes <- greatest.cps*exceed.cpl[continue] # should this subset be removed?
  ### retain.continue <- retained.outcomes*continue 
  ### retain.continue=1 : trial continues to stage 2 and outcome is retained. 
  ### retain.continue=0 : trial stopped at stage 1 or trial continues to stage 2 and outcome is dropped
  ### exceed.bounds <- t(t(ts.s2) > z.alpha)
  ### exceed.bounds.binary <- exceed.bounds*retain.continue # outcomes that exceed final boundary AND were retained AND only for trials that continued to stage 2.
  ### This seems fine, but I think there is still some amibguity regarding what is contributing to the type I error:
  ### If m=2, i.e. we reject H0 only if two outcomes exceed their boundary, then a single outcome exceeding its boundary is not a type I error,
  ### and so one could argue that such an instance should not be counted as contributing to that outcome's share of the type I error.
  
  # The section above, using apply, is clearer but approx. half the speed of the conditional loop below:
  minus.cp <- -cp
  # j<-1
  # greatest.cps.list <- vector("list", nsims.)
  # for(i in (1:nsims.)[continue]){
  #   greatest.cps.list[[j]] <- rank(minus.cp[i,]) <= max.outs.s2.
  #   j <- j+1
  # }
  # greatest.cps <- do.call(rbind, greatest.cps.list)
  # ^ Superceded by code below:
  minus.cp.continue.subset <- minus.cp[continue,]
  ranked.rows <- t(apply(minus.cp.continue.subset, 1, rank))
  greatest.cps <- ranked.rows <= max.outs.s2.
  
  retained.outcomes <- greatest.cps*exceed.cpl[continue,]
  n.obs.s2 <- sum(retained.outcomes)
  exceed.bounds.binary <- t(t(ts.s2[continue,]) > z.alpha)*retained.outcomes
  reject.s1 <- sum(stop.efficacy)
  reject.s2 <- sum(rowSums(exceed.bounds.binary)>=m.)
  p.reject <- (reject.s1+reject.s2)/nsims.
  PET <- sum(stop.any)/nsims.

  ENM.pp <- (K*nsims. + n.obs.s2)/nsims. # Expected no. observations is K*(no. times stop at S1) + sum of all observations in S2 ("n.obs.s2")
  
  # drop: number of outcomes to drop (one approach. The other is to drop all outcomes below a certain CP (cp.l))
  # retain.binary <- t(apply(cp.rank, 1, function(x) x > drop)) # Assuming we drop a fixed number of outcomes, set as "drop".
  # if(drop.outcomes==TRUE){
  #   if(always.drop.==FALSE){
  #     retained.outcomes <- cp > cp.l.
  #     }else{
  #   retained.outcomes <- t(apply(cp, 1, function(x) x==max(x)))
  #     }
  #   }else{
  #     retained.outcomes <- matrix(TRUE, ncol=K., nrow=nsims.)
  #     }
  # retain.continue <- retained.outcomes*continue 
  # retain.continue=1 : trial continues to stage 2 and outcome is retained. 
  # retain.continue=0 : trial stopped at stage 1 or trial continues to stage 2 and outcome is dropped
  # exceed.bounds <- t(apply(ts.s2, 1, function(x) x>z.alpha)) # Too slow -- use line below:
  # exceed.bounds <- t(t(ts.s2) > z.alpha)
  # exceed.bounds.binary <- exceed.bounds*retain.continue # outcomes that exceed final boundary AND were retained AND only for trials that continued to stage 2.
  # This seems fine, but I think there is still some amibguity regarding what is contributing to the type I error:
  # If m=2, i.e. we reject H0 only if two outcomes exceed their boundary, then a single outcome exceeding its boundary is not a type I error,
  # and so one could argue that such an instance should not be counted as contributing to that outcome's share of the type I error.
  
#  reject.s1 <- sum(stop.efficacy)
#  reject.s2 <- sum(rowSums(exceed.bounds.binary)>=m.)
#  prob.reject <- (reject.s1+reject.s2)/nsims.
  
  # Output:
  PET <- sum(stop.any)/nsims.
  
  if(return.optimisation==TRUE){
    typeIerror.diff2 <- sum((alpha.k. - p.reject)^2)
    return(typeIerror.diff2)
  }else{
    return(list(prob.reject=p.reject,
                pet=PET,
                enm.pp=ENM.pp))
  }
}



findCPloserDes <- function(nsims=default.nsims.dtl,
                           max.outcomes=default.max.outcomes.dtl,
                           vars=default.vars.dtl,
                           delta0=default.delta0.dtl,
                           delta1=default.delta1.dtl,
                           delta.true=NULL,
                           reuse.deltas=TRUE,
                           alpha.k=default.alpha.dtl,
                           alpha.combine=TRUE,
                           cp.l=default.cp.l.dtl,
                           cp.u=default.cp.u.dtl,
                           n.min=default.nmin.dtl,
                           n.max=default.nmax.dtl,
                           power=default.power.dtl,
                           rho.vec=default.cor.dtl,
                           working.outs=NULL,
                           fix.n=FALSE)
{
  n.init <- ceiling(n.min+(n.max-n.min)/2)
  
  K <- max.outcomes[1]
  max.outs.s2 <- max.outcomes[2]
  m <- max.outcomes[3]
  
  recycleDeltas <- function(vec, working.outs., K.){
    full.delta.vec <- rep(vec[2], K.)
    full.delta.vec[working.outs.] <- vec[1]
    return(full.delta.vec)
  }
  
  #### Warnings, checks ####
  if(is.null(delta0)){
    warning("No uninteresting treatment effects delta0 supplied. Using delta0=-1000, for all outcomes.", call. = FALSE)
    delta0 <- rep(-1000, K)
  }
  if(is.null(rho.vec)){
    warning("No correlations supplied. Using rho=0.1 for all correlations.", call. = FALSE)
    rho.vec <- rep(0.1, times=sum(1:(K-1)))
  } 
  if(is.null(vars) | length(vars)==1){
    warning("Either zero or one outcome variance supplied. Using var=1 for all outcomes.", call. = FALSE)
    vars <- rep(1, K)
  }
  if(is.null(working.outs)){
    warning("Indices of working outcomes not supplied. Taking indices of working outcomes as outcomes 1 to m.", call. = FALSE)
    working.outs <- 1:m
  }
  if(is.null(max.outs.s2)){
    warning("Maxmimum number of outcomes allowed in stage 2 not supplied. Using the number of outcomes required to show promise, m")
    max.outs.s2 <- m
  }
  if(max.outs.s2 < m){
    stop("Maxmimum number of outcomes allowed in stage 2, max.outcomes[2], must be greater than or equal to the number of outcomes required to show promise, m", call=F)
  }
  if(length(rho.vec)!=1 & length(rho.vec)!=sum(1:(K-1))){
    stop("The number of correlation terms must be equal to 1 (in which case it is equal across all pairs of outcomes) or equal to the number of pairs of outcomes",
         call.=FALSE)
  }
  if(length(rho.vec)==1 & K>2){
    warning("Single value supplied for correlation rho.vec. Using supplied value for all correlations", call=F)
    rho.vec <- rep(rho.vec, sum(1:(K-1)))
  }
  if(alpha.combine==TRUE & length(alpha.k)>1){
    stop("When alpha.combine is set to TRUE, a single overall alpha.k is required, not a vector.", call=F)
  }
  
  if(reuse.deltas==TRUE){
      # !!! IMPORTANT: Currently, the only delta values used are delta1[1] and delta0[2].
      # !!! They are used to obtain the power, which is found given outcome effects equal to  delta1[1] for the first m outcomes and equal to delta0[2] for the remaining K-m outcomes.
      if(length(delta0)==1){
        delta0 <- rep(delta0, 2)
      }
      if(length(delta1)==1){
        delta1 <- rep(delta1, 2)
      }
      return.delta0 <- delta0 
      return.delta1 <- delta1
      rm(delta0, delta1)
      delta0 <- rep(return.delta0[2], K)
      delta1 <- rep(return.delta1[2], K)
      delta0[working.outs] <- return.delta0[1] 
      delta1[working.outs] <- return.delta1[1]
      if(K>2){
        warning("reuse.deltas set to TRUE: As K>2, will take delta1[1] as anticipated effect size for outcomes 1 to m, and \n delta0[2] as anticipated effect size for outcomes m+1 to K")
      }
    if(!is.null(delta.true)){
      means.true <- t(apply(delta.true, 1, recycleDeltas, working.outs.=working.outs, K.=K))
      if(K>2){
        warning("reuse.deltas set to TRUE: As K>2, will take delta.true[1] as true delta for all working outcomes (i.e. 1 to m) and \n delta.true[2] as true delta for all non-working outcomes (i.e. m+1 to K.")
      }
    }
  }  
  if(length(vars)!=K | length(delta0)!=K | length(delta1)!=K){
    stop("The arguments vars, delta0 and delta1 must all have length equal to the number of outcomes, max.outcomes[1]", call.=FALSE)
  }
  

  library(minqa)
  J <- 2
  cov.mat <- createCovMat(J.=J, K.=K, rho.vec.=rho.vec)
  #set.seed(seed)
  ts.global.null <- mvtnorm::rmvnorm(nsims, mean=rep(0, J*K), sigma = cov.mat)
  # Find optimal final bounds for an initial n under the global null, using drop the loser design, and find the type I error at these bounds:

#### Find optimal DtL design ##### 
#  (i.e. find final boundary and sample size that gives correct type I error and power) 
# Use the bisection method to find the n that gives the appropriate power (using drop the loser design):
  n.all <- n.min:n.max
  a <- 1
  b <- length(n.all)
  d <- which(n.all==n.init)
  while(b-a>1){
    r.k <- bobyqa(par=c(2),
                  fn = findR,
                  lower=0.01,
                  upper=10,
                  typeI.power="typeI",
                  ts=ts.global.null, 
                  n.stage=n.all[d],
                  return.optimisation=TRUE,
                  nsims.=nsims,  
                  K.=K,
                  m.=m,
                  max.outs.s2.=max.outs.s2,
                  working.outs.=working.outs,
                  vars.=vars,
                  delta0.=delta0,
                  delta1.=delta1,
                  alpha.k.=alpha.k,
                  cp.l.=cp.l,
                  cp.u.=cp.u
                  # always.drop.=always.drop
    )$par

    pwr.output <- findR(bounds=r.k,
                        typeI.power="power",
                        ts=ts.global.null,
                        n.stage=n.all[d],
                        nsims.=nsims,
                        K.=K,
                        m.=m,
                        max.outs.s2.=max.outs.s2,
                        working.outs.=working.outs,
                        vars.=vars,
                        delta0.=delta0,
                        delta1.=delta1,
                        alpha.k.=alpha.k,
                        cp.l.=cp.l,
                        cp.u.=cp.u
                        # always.drop.=always.drop
    )
    print(paste("Power is ", format(pwr.output$prob.reject, digits=4), " when n per stage is ", n.all[d]), q=F)
    if(pwr.output$prob.reject < power){
      a <- d
      d <- ceiling(a+(b-a)/2)
    } else {
      b <- d
      d <- ceiling(a+(b-a)/2)
    }
  } # end of while
  final.n.stage <- n.all[d]
  print(paste("Final n per stage: ", final.n.stage), q=F)
  final.r.k <- bobyqa(par=c(2),
                               fn = findR,
                               lower=0.01,
                               upper=10,
                               typeI.power="typeI",
                               ts=ts.global.null, 
                               n.stage=final.n.stage,
                               return.optimisation=TRUE,
                               nsims.=nsims,  
                               K.=K,
                               m.=m,
                               max.outs.s2.=max.outs.s2,
                               working.outs.=working.outs,
                               vars.=vars,
                               delta0.=delta0,
                               delta1.=delta1,
                               alpha.k.=alpha.k,
                               cp.l.=cp.l,
                               cp.u.=cp.u
                               # always.drop.=always.drop
                      )$par
  final.pwr <- findR(bounds=final.r.k,
                     typeI.power="power",
                     ts=ts.global.null,
                     n.stage=final.n.stage,
                     nsims.=nsims,
                     K.=K,
                     m.=m,
                     max.outs.s2.=max.outs.s2,
                     working.outs.=working.outs,
                     vars.=vars,
                     delta0.=delta0,
                     delta1.=delta1,
                     alpha.k.=alpha.k,
                     cp.l.=cp.l,
                     cp.u.=cp.u)
  t1.final.n <- findR(bounds=final.r.k,
                      typeI.power="typeI",
                      ts=ts.global.null,
                      n.stage=final.n.stage,
                      nsims.=nsims,
                      K.=K,
                      m.=m,
                      max.outs.s2.=max.outs.s2,
                      working.outs.=working.outs,
                      vars.=vars,
                      delta0.=delta0,
                      delta1.=delta1,
                      alpha.k.=alpha.k,
                      cp.l.=cp.l,
                      cp.u.=cp.u)
  ess.h0.dtl <- t1.final.n$pet*final.n.stage + (1-t1.final.n$pet)*2*final.n.stage
  ess.h1.dtl <- final.pwr$pet*final.n.stage + (1-final.pwr$pet)*2*final.n.stage
  
  
  
  
  

  p.reject.single.stage <- function(bounds,
                                    ts,
                                    alpha,
                                    m.,
                                    K.,
                                    opt){
    ts.single.stage <- ts[,(K.+1):(2*K.)]
    rejections <- rowSums(ts.single.stage > bounds) >= m.
    t1.err <- sum(rejections)/nrow(ts)
    #t1.err.single.stage <-  sum(ts.single.stage[,1] > bounds)/nrow(ts) # first outcome only. 
    # When first outcome only, boundary is ~1.645. Power unchanged.
    if(opt==TRUE){
      alpha.diff <- (t1.err-alpha)^2
      return(alpha.diff)
    }else{
      return(t1.err)
    }
  }
  
  
  #####  Find optimal single-stage design #####
  # (i.e. find boundary and sample size that gives correct type I error and power
  r.k.single.stage <- bobyqa(par=2,
                             fn=p.reject.single.stage,
                             lower=0.01,
                             upper=10,
                             ts=ts.global.null,
                             alpha=alpha.k,
                             m.=m,
                             K.=K,
                             opt=TRUE)$par 
  t1.err.single.stage <- p.reject.single.stage(bounds=r.k.single.stage,
                                               ts=ts.global.null,
                                               alpha=alpha.k,
                                               m.=m,
                                               K.=K,
                                               opt=FALSE)
  # Find sample size with correct power:
  denom <-  vars
  lfc.effects <- delta0
  lfc.effects[working.outs] <- delta1[working.outs]
  current.power.single.stage <- 0
  n.all.single.stage <- n.min:(2*n.max)
  i <- 1
  while(current.power.single.stage < power & i<=length(n.all.single.stage)){
  #  numer.J <- rep(J*n.all[i], each=K)
    numer.J <- rep(n.all.single.stage[i], each=K)
    information.final <- numer.J/denom
    tau <- lfc.effects * sqrt(information.final)
    ts.power.single.stage <- sweep(ts.global.null[, (K+1):(2*K)], 2, tau, "+")
    reject.outcome.binary <- ts.power.single.stage > r.k.single.stage
    current.power.single.stage <- sum(rowSums(reject.outcome.binary)>=m)/nsims
    i <- i + 1 
  }
  power.single.stage <- current.power.single.stage
  n.single.stage <- n.all.single.stage[i-1]
  ess.h0.single.stage <- n.single.stage
  ess.h1.single.stage <- n.single.stage
  

#### Fix power etc. for N: ####
if(fix.n==TRUE){
  # Increase sample size for design with the lower sample size of the two designs and calculate the power:
  N.dtl <- J*final.n.stage
  if(N.dtl==n.single.stage){
    # If sample size is already identical for both drop the loser and single-stage designs:
    final.n.stage.max.n <- final.n.stage
    pwr.dtl.max.n <- final.pwr$prob.reject
    t1.dtl.max.n <- t1.final.n$prob.reject
    ess.h0.dtl.max.n <- ess.h0.dtl
    ess.h1.dtl.max.n <- ess.h1.dtl
    
    n.single.stage.max.n <- n.single.stage 
    power.single.stage.max.n <- power.single.stage 
    t1.err.single.stage.max.n <- t1.err.single.stage
    ess.h0.single.stage.max.n <- ess.h0.single.stage
    ess.h1.single.stage.max.n <- ess.h1.single.stage 
  }
  if(N.dtl > n.single.stage){
    # If sample size for single stage design is lower:
    final.n.stage.max.n <- final.n.stage
    pwr.dtl.max.n <- final.pwr$prob.reject
    t1.dtl.max.n <- t1.final.n$prob.reject
    ess.h0.dtl.max.n <- ess.h0.dtl
    ess.h1.dtl.max.n <- ess.h1.dtl
    
    n.single.stage.max.n <- N.dtl
    t1.err.single.stage.max.n <- t1.err.single.stage # type I error doesn't change for single-stage when n changes.
    # Recalculate power and ESS for single stage design:
    numer.J <- rep(n.single.stage.max.n, each=K)
    information.max.n <- numer.J/denom
    tau <- lfc.effects * sqrt(information.max.n)
    ts.power.single.stage <- sweep(ts.global.null[, (K+1):(2*K)], 2, tau, "+")
    reject.outcome.binary <- ts.power.single.stage > r.k.single.stage
    power.single.stage.max.n <- sum(rowSums(reject.outcome.binary)>=m)/nsims
    ess.h0.single.stage.max.n <- n.single.stage.max.n
    ess.h1.single.stage.max.n <- n.single.stage.max.n
  }
  if(N.dtl < n.single.stage){
    # If sample size for drop the loser design is lower:
    final.n.stage.max.n <- n.single.stage/J # divide by J b/c final.n.stage.max.n is n per stage, while n.single.stage is for both stages.
    # Recalculate power and ESS for drop the loser design:
    power.pet.dtl.max.n <- findR(bounds=final.r.k,
                                 typeI.power="power",
                                 ts=ts.global.null,
                                 n.stage=final.n.stage.max.n,
                                 nsims.=nsims,
                                 K.=K,
                                 m.=m,
                                 max.outs.s2.=max.outs.s2,
                                 working.outs.=working.outs,
                                 vars.=vars,
                                 delta0.=delta0,
                                 delta1.=delta1,
                                 alpha.k.=alpha.k,
                                 cp.l.=cp.l,
                                 cp.u.=cp.u)
    
    t1.pet.dtl.max.n <- findR(bounds=final.r.k,
                              typeI.power="typeI",
                              ts=ts.global.null,
                              n.stage=final.n.stage.max.n,
                              nsims.=nsims,
                              K.=K,
                              m.=m,
                              max.outs.s2.=max.outs.s2,
                              working.outs.=working.outs,
                              vars.=vars,
                              delta0.=delta0,
                              delta1.=delta1,
                              alpha.k.=alpha.k,
                              cp.l.=cp.l,
                              cp.u.=cp.u)
    pwr.dtl.max.n <- power.pet.dtl.max.n$prob.reject
    t1.dtl.max.n <- t1.pet.dtl.max.n$prob.reject
    ess.h0.dtl.max.n <- t1.pet.dtl.max.n$pet*final.n.stage.max.n + (1-t1.pet.dtl.max.n$pet)*2*final.n.stage.max.n
    ess.h1.dtl.max.n <- power.pet.dtl.max.n$pet*final.n.stage.max.n + (1-power.pet.dtl.max.n$pet)*2*final.n.stage.max.n
    
    n.single.stage.max.n <- n.single.stage 
    power.single.stage.max.n <- power.single.stage 
    t1.err.single.stage.max.n <- t1.err.single.stage
    ess.h0.single.stage.max.n <- ess.h0.single.stage
    ess.h1.single.stage.max.n <- ess.h1.single.stage 
  }
  
  N.max.n <- c(J*final.n.stage.max.n, n.single.stage.max.n)
  ess0.max.n <- c(ess.h0.dtl.max.n, ess.h0.single.stage.max.n)
  ess1.max.n <- c(ess.h1.dtl.max.n, ess.h1.single.stage.max.n)
  power.max.n <- c(pwr.dtl.max.n, power.single.stage.max.n)
  t1.max.n <- c(t1.dtl.max.n, t1.err.single.stage.max.n)

}else{
  N.max.n <- rep(NA, 2)
  ess0.max.n <- rep(NA, 2)
  ess1.max.n <- rep(NA, 2)
  power.max.n <- rep(NA, 2)
  t1.max.n <- rep(NA, 2)
}
#### True delta ####
if(!is.null(delta.true)){
  nrows <- nrow(means.true)
  true.results.list <- vector("list", nrows)
  power.true.single.stage <- numeric(nrows)
  for(i in 1:nrows){
    true.results.list[[i]] <- unlist(findR(bounds=final.r.k,
                     typeI.power="truedelta",
                     ts=ts.global.null,
                     n.stage=final.n.stage,
                     nsims.=nsims,
                     K.=K,
                     m.=m,
                     max.outs.s2.=max.outs.s2,
                     working.outs.=working.outs,
                     vars.=vars,
                     delta1.=delta1,
                     delta.true.=means.true[i,],
                     alpha.k.=alpha.k,
                     cp.l.=cp.l,
                     cp.u.=cp.u))
    tau.true.current <- means.true[i, ] * sqrt(information.final)
    ts.true.current <-  sweep(ts.global.null[, (K+1):(2*K)], 2, tau.true.current, "+")
    true.reject.power.single.stage <- ts.true.current > r.k.single.stage
    power.true.single.stage[i] <- sum(rowSums(true.reject.power.single.stage)>=m)/nsims
    }
  true.results.mat <- do.call(rbind, true.results.list)
 # browser()
  dtl.true <- data.frame(true.results.mat, delta.true, means.true)
  dtl.true$ess <- final.n.stage*dtl.true$pet + J*final.n.stage*(1-dtl.true$pet)
  dtl.true$enm.total <- dtl.true$enm.pp * final.n.stage
  dtl.true$design <- "DtL"
  single.stage.pet <- rep(0, nrows)
  single.stage.enm.pp <- rep(K, nrows)
  single.stage.true <- data.frame(prob.reject=power.true.single.stage, pet=single.stage.pet, enm.pp=single.stage.enm.pp, delta.true, means.true)
  single.stage.true$ess <- n.single.stage
  single.stage.true$enm.total <- single.stage.true$enm.pp * n.single.stage
  single.stage.true$design <- "Single stage"
  #output.true <- cbind(delta.true, means.true, true.results.mat, power.true.single.stage, true.results.mat[,1]/power.true.single.stage)
  #colnames(output.true) <- c("mu.working", "mu.nonworking", paste("mu.", 1:ncol(means.true), sep=""), "p.reject.dtl", "pet.dtl", "p.reject.single.s", "p.reject.ratio")
  output.true <- rbind(dtl.true, single.stage.true)
  colnames(output.true) <- c("prob.reject", "pet", "enm.pp", "mu.working", "mu.nonworking", paste("mu.", 1:ncol(means.true), sep=""), "ess", "enm", "design")
  ratios <- data.frame(ess.ratio=dtl.true$ess/single.stage.true$ess,
                       enm.ratio=dtl.true$enm.total/single.stage.true$enm.total)
  output.true.ratios <- data.frame(dtl.true$prob.reject, power.true.single.stage, delta.true, means.true, ratios)
  names(output.true.ratios) <- c("p.reject.dtl", "p.reject.ss", "mu.working", "mu.nonworking", paste("mu.", 1:ncol(means.true), sep=""), "ess.ratio", "enm.ratio")
}
  
#### output  ####
  # Collate results for output:
  #typeIerr.k.c <- rbind(typeIerr.k, t1.err.single.stage)
  #colnames(typeIerr.k.c) <- paste("typeIerr.k", 1:length(alpha.k), sep="")
  typeIerr.total.c <- c(sum(t1.final.n$prob.reject), t1.err.single.stage)
  power.c <- c(final.pwr$prob.reject, power.single.stage)
  final.bounds <- rbind(final.r.k, r.k.single.stage)
  ess.h0 <- c(ess.h0.dtl, ess.h0.single.stage)
  ess.h1 <- c(ess.h1.dtl, ess.h1.single.stage)
  enm.pp.h0 <- c(t1.final.n$enm.pp, K)
  enm.pp.h1 <- c(final.pwr$enm.pp, K)
  enm.tot.h0 <- enm.pp.h0*c(final.n.stage, n.single.stage) # Total ENM is (ENM per person)*(n per stage) for DtL, and K*N for single stage.
  enm.tot.h1 <- enm.pp.h1*c(final.n.stage, n.single.stage)
  
  #colnames(final.bounds) <- paste("r.k", 1:K, sep="") # Only needed when bounds differ
  final.n.vec <- c(final.n.stage, NA)
  final.N.vec <- c(J*final.n.stage, n.single.stage)
  
  design <- c("Drop the loser", "Single stage")
  
#  browser()
  design.results <- cbind(final.bounds, final.n.vec, final.N.vec, ess.h0, ess.h1, enm.pp.h0, enm.pp.h1, enm.tot.h0, enm.tot.h1, typeIerr.total.c, power.c, N.max.n, ess0.max.n, ess1.max.n, power.max.n, t1.max.n)
  colnames(design.results) <- c("r.k", "n", "N", "ESS0", "ESS1", "ENM.pp.0", "ENM.pp.1", "ENM0", "ENM1", "typeIerr", "power", "N.max.n", "ESS0.max.n", "ESS1.max.n", "power.max.n", "typeIerr.max.n")
  rownames(design.results) <- c("Drop", "Single stage")
  design.results <- as.data.frame(design.results)
  design.results$design <- design
  # Shared results:
  if(reuse.deltas==TRUE){
    des.chars <- data.frame(K, max.outs.s2, m, cp.l, cp.u, t(return.delta0), t(return.delta1), sum(alpha.k), power,  rho.vec[1]) 
    colnames(des.chars) <- c("K", "max.outs.s2", "m", "cp.l", "cp.u", "delta0.1", "delta0.2", "delta1.1", "delta1.2", "alpha", "req.power", "cor")
  }else{
    des.chars <- data.frame(K, max.outs.s2, m, cp.l, cp.u, t(delta0), t(delta1), sum(alpha.k), power,  rho.vec[1]) # This includes all delta0 and delta1 values (use this if specifying separate delta0/1 values for each outcome)
    colnames(des.chars) <- c("K", "max.outs.s2", "m", "cp.l", "cp.u", paste("delta0.k", 1:K, sep=""), paste("delta1.k", 1:K, sep=""), "alpha", "req.power", "cor")
  }
  des.chars$ess0.ratio <- ess.h0.dtl/ess.h0.single.stage
  des.chars$ess1.ratio <- ess.h1.dtl/ess.h1.single.stage
  des.chars$enm0.pp.ratio <- t1.final.n$enm.pp/K
  des.chars$enm1.pp.ratio <- final.pwr$enm.pp/K
  des.chars$enm0.ratio <- enm.tot.h0[1]/enm.tot.h0[2]
  des.chars$enm1.ratio <- enm.tot.h1[1]/enm.tot.h1[2]
  output <- list(input=des.chars,
                 results=design.results)
                 #paths=rbind(final.pwr$paths, pwr.nodrop.output$paths)
                 #cp=pwr.output$cp
  if(!is.null(delta.true)){  
    output$true.results=output.true
    output$true.ratios=output.true.ratios
    }
  return(output)
}




# Function to tidy up list of outputs from main function :
tidyOutput <- function(output.list){
  input.list <- lapply(output.list, function(x) x$input)
  input.mat <- do.call(rbind, input.list)
  input.mat <- data.frame(input.mat)
  input.mat$km <- paste("K=", input.mat$K, ", m=", input.mat$m, sep="")
  input.mat$max.s2.prop <- "K"
  input.mat$max.s2.prop[input.mat$K-input.mat$max.outs.s2==1] <- "K-1"
  input.mat$max.s2.prop[input.mat$max.outs.s2/input.mat$K==0.5] <- "K/2"
  input.mat$outs.names <- paste("K=", input.mat$K, ", Kmax=", input.mat$max.outs.s2, ", m=", input.mat$m, sep="")
  input.repeated <- input.mat[rep(1:nrow(input.mat), each=2),]
  results.list <- lapply(output.list, function(x) x$results)
  results.mat <- do.call(rbind, results.list)
  results.df <- data.frame(results.mat)
  #results.df <- data.frame(results.mat, type=rownames(results.mat))
  #powerdiff <- unlist(lapply(delta0.par1, function(x) x$results["Drop", "power"] - x$results["Single stage", "power"]))
  return(list(input=input.mat,
              results=cbind(input.repeated, results.df)))
}


createSearchMatrix <- function(K.Kmax.m.data.frame=K.Kmax.m.df, parameter.values){
  search.matrix <- as.matrix(K.Kmax.m.data.frame) %x% rep(1, length(parameter.values))
  search.matrix <- cbind(search.matrix, rep(parameter.values, nrow(K.Kmax.m.data.frame)))
  search.matrix
}


#### Plot/tabulate output ####
# Plot ENM ratio using tidied output, as parameters vary:
plotENM <- function(tidied.output, param, xaxis){
  param <- ensym(param)
  pl <- ggplot(data=tidied.output$input,  mapping=aes(x=!!param, y=enm1.ratio, col=km, linetype=max.s2.prop))+
    geom_line(size=1)+
    geom_hline(yintercept=1,
               linetype="dashed")+
    labs(title="ENM ratio under LFC",
         col="K, m",
         linetype="Kmax",
         y=expression(paste("(", ENM[DtL],"/", ENM[single], ")",  " | LFC")),
         x=xaxis)
  pl
}

# Plot ESS ratio using tidied output, for any parameter:
plotESS <- function(tidied.output, param, xaxis){
  param <- ensym(param)
  pl <- ggplot(data=tidied.output$input,  mapping=aes(x=!!param, y=ess1.ratio, col=km, linetype=max.s2.prop))+
    geom_line(size=1)+
    geom_hline(yintercept=1,
               linetype="dashed")+
    labs(title="ESS ratio under LFC",
         col="K, m",
         linetype="Kmax",
         y=expression(paste("(", ESS[DtL],"/", ESS[single], ")",  " | LFC")),
         x=xaxis)
  pl
}

plotESSENM <-  function(tidied.output, param, xaxis){
  library(gridExtra)
  param <- ensym(param)
  plot.ess <- ggplot(data=tidied.output$input,  mapping=aes(x=!!param, y=ess1.ratio, col=km, linetype=max.s2.prop))+
    geom_line(size=1)+
    geom_hline(yintercept=1,
               linetype="dashed")+
    labs(title="ESS ratio under LFC",
         col="K, m",
         linetype="Kmax",
         y=expression(paste("(", ESS[DtL],"/", ESS[single], ")",  " | LFC")),
         x=xaxis)+
  theme(legend.position = "none")
  
 plot.enm <- ggplot(data=tidied.output$input,  mapping=aes(x=!!param, y=enm1.ratio, col=km, linetype=max.s2.prop))+
   geom_line(size=1)+
   geom_hline(yintercept=1,
              linetype="dashed")+
   labs(title="ENM ratio under LFC",
        col="K, m",
        linetype=expression(paste(K[max])),
        y=expression(paste("(", ENM[DtL],"/", ENM[single], ")",  " | LFC")),
        x=xaxis)#+
   #theme(legend.position = "bottom")
 both.plots <- grid.arrange(plot.ess, plot.enm, ncol=2, widths=c(1.9,2.5))
 both.plots
}
# plotESSENM(tidyOutput(output.cor), param=cor, xaxis="correlation")

#### Plots for changing true effects ####
plotTrueENMdtl <- function(raw.output){
  ggplot(raw.output$true.results[raw.output$true.results$design=="DtL", ], aes(x=mu.1, y=mu.2))+
    geom_raster(aes(fill = enm))+
    scale_x_continuous(breaks=sort(unique(raw.output$true.results$mu.1)))+
    scale_y_continuous(breaks=sort(unique(raw.output$true.results$mu.2)))+
    labs(title=expression(paste(ENM[DtL], ", powered for ")),
         subtitle=bquote( mu[1] == .(raw.output$input$delta1.1) ~~~ mu[2] == .(raw.output$input$delta0.2)),
         y=expression(paste(mu[2])),
         x=expression(paste(mu[1]))
    )+
    geom_text(aes(label = round(enm, 2)), size=5) +
    scale_fill_gradient(low = "white", high = "darkred") 
}


###### plot ratios of ENM_dtl/ENM_ss and ESS_dtl/ESS_ss ####
plotTrueENMratio <- function(raw.output){
  ggplot(raw.output$true.ratios, aes(x=mu.1, y=mu.2))+
    geom_raster(aes(fill = enm.ratio))+
    scale_x_continuous(breaks=sort(unique(raw.output$true.ratios$mu.1)))+
    scale_y_continuous(breaks=sort(unique(raw.output$true.ratios$mu.2)))+
    labs(title=bquote("ENM"[DtL]*"/"*"ENM"[single]*", powered for "*mu[1] == .(raw.output$input$delta1.1)*","~ mu[2] == .(raw.output$input$delta0.2)),
         #subtitle=bquote( mu[1] == .(raw.output$input$delta1.1)*"," ~ mu[2] == .(raw.output$input$delta0.2)),
         y=expression(paste(mu[2])),
         x=expression(paste(mu[1])),
         fill="ENM ratio"
    )+
    geom_text(aes(label = round(enm.ratio, 2)), size=5) +
    scale_fill_gradient2(midpoint=1, low = "darkred", mid="white", high = "darkblue")+
    coord_cartesian(expand=0)
}


plotTrueESSratio <- function(raw.output, method){
  library(ggplot2)
  plot1 <- ggplot(raw.output$true.ratios, aes(x=mu.1, y=mu.2))+
    geom_raster(aes(fill = ess.ratio))+
    scale_x_continuous(breaks=sort(unique(raw.output$true.ratios$mu.1)))+
    scale_y_continuous(breaks=sort(unique(raw.output$true.ratios$mu.2)))+
    labs(#subtitle=bquote( mu[1] == .(raw.output$input$delta1.1)*"," ~ mu[2] == .(raw.output$input$delta0.2)),
         y=expression(paste(mu[2])),
         x=expression(paste(mu[1])))+
    geom_text(aes(label = round(ess.ratio, 2)), size=5) +
    scale_fill_gradient2(midpoint=1, low = "darkred", mid="white", high = "darkblue")+
    coord_cartesian(expand=0)
  if(method=="mo"){
    plot1 <- plot1 +
      labs(title=expression(paste(ESS[MO],"/",ESS[comp], ", powered for ")),
           fill="ESS ratio"
           )
  }
  if(method=="dtl"){
    plot1 <- plot1 +
      labs(title=bquote("ESS"[DtL]*"/"*"ESS"[single]*", powered for "*mu[1] == .(raw.output$input$delta1.1)*","~ mu[2] == .(raw.output$input$delta0.2)),
           fill="ESS ratio"
           )
  }
  plot1
}


plotTrueESSENMratios <- function(raw.out, cols=1){
  essplot <- plotTrueESSratio(raw.out, method="dtl")
  enmplot <- plotTrueENMratio(raw.out)
  both.plots <- grid.arrange(essplot, enmplot, ncol=cols)
  both.plots
}



plotTruePreject <- function(raw.output, method){
  library(ggplot2)
  library(gridExtra)
  if(method=="mo"){
    new.approach.df <- raw.output$true.results[raw.output$true.results$design=="MO", ]
    old.approach.df <- raw.output$true.results[raw.output$true.results$design=="Composite", ]
  }
  if(method=="dtl"){ # XXX CHECK these labels ####
    new.approach.df <- raw.output$true.results[raw.output$true.results$design=="DtL", ]
    old.approach.df <- raw.output$true.results[raw.output$true.results$design=="Single stage", ]
  }
  plot1 <- ggplot(new.approach.df, aes(x=mu.1, y=mu.2))+
    geom_raster(aes(fill = prob.reject))+
    scale_x_continuous(breaks=sort(unique(new.approach.df$mu.1)))+
    scale_y_continuous(breaks=sort(unique(new.approach.df$mu.2)))+
    labs(#subtitle=bquote( mu[1] == .(raw.output$input$delta1.1) ~~~ mu[2] == .(raw.output$input$delta0.2)),
         y=expression(paste(mu[2])),
         x=expression(paste(mu[1])))+
    geom_text(aes(label = round(prob.reject, 2)), size=5) +
    scale_fill_gradient(low = "white", high = "grey40")+
    coord_cartesian(expand = 0) +
    theme(legend.position = "none") 
  
  plot2 <- ggplot(old.approach.df, aes(x=mu.1, y=mu.2))+
    geom_raster(aes(fill = prob.reject))+
    scale_x_continuous(breaks=sort(unique(old.approach.df$mu.1)))+
    scale_y_continuous(breaks=sort(unique(old.approach.df$mu.2)))+
    labs(#subtitle=bquote( mu[1] == .(raw.output$input$delta1.1) ~~~ mu[2] == .(raw.output$input$delta0.2)),
         y=expression(paste(mu[2])),
         x=expression(paste(mu[1])))+
    geom_text(aes(label = round(prob.reject, 2)), size=5) +
    scale_fill_gradient(low = "white", high = "grey40") +
    coord_cartesian(expand = 0) 
  if(method=="mo"){
    plot1 <- plot1 +
      # labs(title=expression(paste(P, "(reject ", H[0], ")", [MO], ", powered for ")),
      #   fill=expression(paste(P, "(reject ", H[0], ")")))
      labs(title="P(reject H0) MO, powered for ",
           fill=expression(paste("P(reject ", H[0])))
    plot2 <- plot2 +
      labs(title="P(reject H0) composite,  powered for ",
           fill=expression(paste("P(reject ", H[0], ")")))
  }
  if(method=="dtl"){
    plot1 <- plot1 +
      # labs(title=expression(paste(P, "(reject ", H[0], ")", [MO], ", powered for ")),
      #   fill=expression(paste(P, "(reject ", H[0], ")")))
      # labs(title=expression(paste("P(reject ", H[0],")"[DtL], ", powered for")),
      #      fill=expression(paste("P(reject ", H[0], ")")))
      labs(title=bquote("R("*mu[1] == .(raw.output$input$delta1.1)*","~ mu[2] == .(raw.output$input$delta0.2)*")"[DtL]),
           fill=expression(paste("R(", mu, ")")))
    plot2 <- plot2 +
      # labs(title=expression(paste("P(reject ", H[0],")"[single], ", powered for")),
      #      fill=expression(paste("P(reject ", H[0], ")")))
      labs(title=bquote("R("*mu[1] == .(raw.output$input$delta1.1)*","~ mu[2] == .(raw.output$input$delta0.2)*")"[single]),
           fill=expression(paste("R(", mu, ")")))
  }
  both.plots <- grid.arrange(plot1, plot2, widths=c(2.4,3))
  both.plots
}


########### Create subset of true results  #########
# subset to just the true values of interest:
createSubset <- function(raw.output, delta.matrix){
  index <- apply(delta.matrix, 1, function(x) which(abs(raw.output$true.ratios$mu.working-x[1])<0.0001 & abs(raw.output$true.ratios$mu.nonworking-x[2])<0.0001))
  df.subset <-   raw.output$true.ratios[index, ]
  df.subset
}


# Find interim stopping boundaries ####
findZ <- function(cp,
                  j=1,
                  J=2,
                  n.stage,
                  vars,
                  z.alpha,
                  d.1){
  numer <- 1*n.stage
  denom <-  vars
  information.j <- numer/denom
  numer.J <- J*n.stage
  information.final <- numer.J/denom
  Z <- (sqrt(information.final-information.j)*qnorm(cp) +  z.alpha*sqrt(information.final) - (information.final-information.j)*d.1) / sqrt(information.j)
  Z
}

findZ(cp=0.1,
      n.stage = 10,
      vars=1,
      z.alpha = 2,
      d.1=0.4)
