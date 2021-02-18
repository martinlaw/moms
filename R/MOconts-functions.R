# Original code used to find results in manuscript.  Only difference compared
# to findDesis that this function also obtains composite results are obtained.
  findDesSubmission <- function(K=default.K,
                      m=default.m,
                      J=default.J,
                      rho.vec=default.cor,
                      nsims=default.nsims,
                      wang.delta=0,
                      alpha=default.alpha,
                      power=default.power,
                      delta0=default.delta0,
                      delta1=default.delta1,
                      reuse.deltas=TRUE,
                      delta.true=NULL,
                      reuse.true.deltas=TRUE,
                      vars.true=NULL,
                      vars=NULL,
                      working.outs=NULL,
                      #delta.matrix=NULL,
                      #fix.n=FALSE,
                      maxn.stage=200,
                      return.boundaries=FALSE,
                      return.ts=FALSE
                      )
  {
    recycleDeltas <- function(vec, working.outs., K.){
      full.delta.vec <- rep(vec[2], K.)
      full.delta.vec[working.outs.] <- vec[1]
      return(full.delta.vec)
    }
  #### Warnings, checks: ####
    if(is.null(delta0)){
      warning("No uninteresting treatment effects delta0 supplied. Using delta0=0, for all outcomes.", call. = FALSE)
      delta0 <- rep(0, K)
    }
    if(is.null(rho.vec)){
      warning("No correlations supplied. Using rho=0.5 for all correlations.", call. = FALSE)
      rho.vec <- rep(0.5, times=sum(1:(K-1)))
    }
    if(length(rho.vec)==1 & K>2){
      warning("Single value supplied for correlations supplied. Using this value for all correlations.", call. = FALSE)
      rho.vec <- rep(rho.vec, times=sum(1:(K-1)))
    }
    if(is.null(vars)){
      warning("No outcome variances supplied. Using vars=1 for all outcomes.", call. = FALSE)
      vars <- rep(1, K)
    }
    if(is.null(vars.true) & !is.null(delta.true)){
      warning("No TRUE outcome variances supplied. Using anticipated vars as true vars (default is 1).", call. = FALSE)
      vars.true <- vars
    }
    if(is.null(working.outs)){
      warning("Indices of working treatments not supplied. Taking indices of working treatments as treatments 1 to m.", call. = FALSE)
      working.outs <- 1:m
    }
    if(is.null(m)){
      warning("Number of outcomes required to show promise, m, not given. Using number of working treatments as treatments.", call. = FALSE)
      m <- length(working.outs)
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
        warning("reuse.deltas set to TRUE: If K>2, will take delta0[1] as delta0 for all working outcomes and delta0[2] as delta0 for all non-working outcomes. Ditto delta1")
      }

    }else{
      if(!is.null(delta.true)){
        warning("CAUTION: reuse.deltas set to FALSE, but values supplied for true delta. This could cause problems re. delta.true and means.true")
      }
    }

    if(!is.null(delta.true)){
      if(reuse.true.deltas==TRUE){
        means.true <- t(apply(delta.true, 1, recycleDeltas, working.outs.=working.outs, K.=K))
          if(K>2){
          warning("reuse.deltas set to TRUE: If K>2, will take delta.true[1] as true delta for all working outcomes and delta.true[2] as true delta for all non-working outcomes.")
        }
      }else{
        means.true <- delta.true
        if(ncol(delta.true)!=K){
          stop("The number of columns (outcomes) in delta.true does not equal K and reuse.delta.true==FALSE.")
        }
      }
    }

    # Checks:
    if(length(rho.vec)!=sum(1:(K-1))) stop("Number of rho values, i.e. length of rho.vec, should be equal to sum(1:(K-1)) ", call. = FALSE)
    if(length(delta0)!=K & !is.null(delta0)) stop ("Number of supplied uninteresting treatment effects, i.e. delta0, should be equal to K, the number of outcomes",
                                                   call. = FALSE)
    if(length(delta1)!=K) stop ("Number of supplied treatment effects, i.e. delta1, should be equal to K, the number of outcomes", call. = FALSE)
    if(length(vars)!=K) stop ("Number of outcome variances, i.e. vars, should be equal to K (the number of outcomes)", call. = FALSE)
  #### P(rejection) ####

    ##### new pRejectFast, with shared constant C:
    pRejectFastCommonC <- function(const,
                                   J=J,
                                   K=K,
                                   m=m,
                                   wang.delta,
                                   ts,
                                   nsims,
                                   prob.only=TRUE,
                                   alpha=alpha,
                                   composite=FALSE
    ){
      e.vec <-  const*((1:J)/J)^(wang.delta-0.5)
      f.vec <- -e.vec
      f.vec[length(f.vec)] <- e.vec[length(e.vec)]
      if(J==1){
        if(composite==FALSE) {
          go.overall <- colSums(t(ts) > e.vec)>=m
        }else{
          go.overall <- ts > e.vec
        }
        prob.reject <- sum(go.overall)/nsims
        minimise.prob <- abs(prob.reject - alpha)
        expd.no.stages <- 1
      } else { # ie if J!=1:
        if(composite==FALSE){
          nogo.overall <- vector("list", J)
          go.overall <- vector("list", J)
          first.m.outs.exceed.e <- vector("list", J)
          other.combn.exceed.e <- vector("list", J)
          for(j in 1:J){
            nogo.overall[[j]] <- colSums(t(ts[,(1+(j-1)*K):(j*K)]) < f.vec[j])>=(K-m+1)
            go.overall[[j]] <- colSums(t(ts[,(1+(j-1)*K):(j*K)]) > e.vec[j])>=m
            first.m.outs.exceed.e[[j]] <- colSums(t(ts[,(1+(j-1)*K):(m+(j-1)*K)]) > e.vec[j])==m
            other.combn.exceed.e[[j]] <- go.overall[[j]]==TRUE & first.m.outs.exceed.e[[j]]==FALSE
          }
          # Are m or more boundaries crossed? Row=simulation, col=stage
          nogo.trial.binary <- t(do.call(rbind, nogo.overall))
          go.trial.binary <- t(do.call(rbind, go.overall))
          first.m.outs.exceed.e <- t(do.call(rbind, first.m.outs.exceed.e))
          other.combn.exceed.e <-  t(do.call(rbind, other.combn.exceed.e))
        } else{ # if composite==TRUE
          # Only J boundaries for composite:
          nogo.trial.binary <-  t(t(ts) < f.vec)
          go.trial.binary <- t(t(ts) > e.vec)
         # browser()
        }
        # The first stage at which a nogo decision is made (and analogous for go):
        # Add extra column of 1's so that there is always some maximum even if NOGO boundary is never crossed:
        nogo.trial.binary.plus <- cbind(nogo.trial.binary, 1)
        # Add extra column of 1's so that there is always some maximum even if GO boundary is never crossed:
        go.trial.binary.plus <- cbind(go.trial.binary, 1)
        first.nogo.stage <- Rfast::rowMaxs(nogo.trial.binary.plus, value=FALSE) # Rfast
        first.go.stage <- Rfast::rowMaxs(go.trial.binary.plus, value=FALSE) # Rfast
        first.stop.stage <- cbind(first.nogo.stage, first.go.stage)
        mode(first.stop.stage) <- "numeric"
        # Does the trial make a nogo or a go decision first?
        # Final decision: 1=nogo, 2=go
        final.decision <- Rfast::rowMins(first.stop.stage, value=FALSE) # value=FALSE returns indices
        prob.reject <- sum(final.decision==2)/nsims
        minimise.prob <- abs(prob.reject - alpha)^2 # 29th Jul: added square.
        stop.stage <- Rfast::rowMins(first.stop.stage, value=TRUE) #Rfast
        expd.no.stages <- sum(stop.stage)/nsims
      } # end of if J==1 else
      if(prob.only==TRUE){
        return(minimise.prob)
      } else{
        if(composite==FALSE & J!=1){
          go.decision.index <- which(final.decision==2)
          correct.go <- rep(NA, length(go.decision.index))
          incorrect.go <- rep(NA, length(go.decision.index))
          for(i in 1:length(go.decision.index)){
            correct.go[i] <- first.m.outs.exceed.e[go.decision.index[i], stop.stage[go.decision.index[i]]]
            incorrect.go[i] <- other.combn.exceed.e[go.decision.index[i], stop.stage[go.decision.index[i]]]
          }
          prob.correct.go <- sum(correct.go)/nsims
          prob.incorrect.go <- sum(incorrect.go)/nsims
        }else{
          prob.correct.go <- NA
          prob.incorrect.go <- NA
        }
        return(list(prob.reject=prob.reject,
                    expd.no.stages=expd.no.stages,
                    f.vec=f.vec,
                    e.vec=e.vec,
                    prob.correct.go=prob.correct.go,
                    prob.incorrect.go=prob.incorrect.go)
        )
      }
    } # end of function

  trueReject <- function(
      J=J,
      K=K,
      m=m,
      ts=ts,
      nsims,
      prob.only=FALSE,
      alpha=alpha,
      means.true=means.true,
      vars.true=vars.true,
      f.vec=pwr.ess$f.vec,
      e.vec=pwr.ess$e.vec,
      n.final=n.final,
      composite=FALSE
    ){
      # Need test statistic for true effects:
      means.true <- rep(means.true, times=J)
      # Need to multiply each outcome's true trt effect by sqrt(j*n/vars), where vars is the variance of outcome K. Vector should end up having length J*K
      denom <-  rep(vars.true, times=J)
      numer <- rep((1:J)*n.final, each=K)
      information <- numer/denom
      tau.true <- means.true*sqrt(information)
      ts.true <- sweep(ts, 2, tau.true, "+")   # Add the above tau vector to every row in the matrix ts
      if(J==1){
        if(composite==FALSE) {
          go <- t(apply(ts.true, 1, function(x) x > e.vec))
          go.overall <- apply(go, 1, function(x) sum(x)>=m)
        } else{
          ts.composite <- rowSums(ts.true)
          go.overall <- ts.composite > e.vec
        }
        prob.reject <- sum(go.overall)/nsims
        ess <- n.final
      } else {
      if(composite==FALSE){
        nogo <- t(apply(ts.true, 1, function(x) x < f.vec))
        go <- t(apply(ts.true, 1, function(x) x > e.vec))
        nogo.overall <- vector("list", J)
        go.overall <- vector("list", J)
        for(j in 1:J){
          # Subset to the K outcomes for stage j:
          current.stage.nogo <- nogo[, (1+(j-1)*K):(j*K)]
          current.stage.go <- go[, (1+(j-1)*K):(j*K)]
          # Would the trial stop for either go or nogo at stage j?
          nogo.overall[[j]] <- apply(current.stage.nogo, 1, function(x) sum(x)>=(K-m+1))
          go.overall[[j]] <- apply(current.stage.go, 1, function(x) sum(x)>=m)
        }
        # Is a boundary crossed? Row=simulation, col=stage
        nogo.trial.binary <- t(do.call(rbind, nogo.overall))
        go.trial.binary <- t(do.call(rbind, go.overall))
      } else {
        # Sum the test statistics at each stage to form the composite test statistics
        ts.composite <- matrix(NA, nrow=nsims, ncol=J)
        for(i in 1:J){
          ts.composite[,i] <- rowSums(ts.true[,(1+(i-1)*K):(i*K)])
        }
        nogo.trial.binary <- t(apply(ts.composite, 1, function(x) x < f.vec)) + 0
        go.trial.binary <- t(apply(ts.composite, 1, function(x) x > e.vec)) + 0
      }
      # The first stage at which a nogo decision is made (and analogous for go):
      first.nogo.stage <- apply(nogo.trial.binary, 1, function(x) min(which(x==1), Inf))
      first.go.stage <- apply(go.trial.binary, 1, function(x) min(which(x==1), Inf))
      first.stop.stage <- cbind(first.nogo.stage, first.go.stage)
      # Does the trial make a nogo or a go decision first?
      # Final decision: 1=nogo, 2=go
      final.decision <- apply(first.stop.stage, 1, which.min)
      prob.reject <- sum(final.decision==2)/nsims
      stop.stage <- apply(first.stop.stage, 1, min)
      expd.no.stages <- sum(stop.stage)/nsims
      ess <- expd.no.stages*n.final
      #minimise.prob <- abs(prob.reject - alpha)
      } # end of J==1 else
      return(c(prob.reject, ess))
    } # end of function

    constToBounds <- function(const, J., wang.d){
      e.vec <- const*((1:J.)/J.)^(wang.d-0.5)
      f.vec <- -e.vec
      f.vec[length(f.vec)] <-e.vec[length(e.vec)]
      bounds <- list(e=e.vec, f=f.vec)
      return(bounds)
    }

    createTrueTS <- function(ts., mu., vars., J., K., n., composite=FALSE){
      mu.full.vec <- rep(mu., times=J.)
      # Need to multiply each outcome's true trt effect by sqrt(j*n/vars), where vars is the variance of outcome K. Vector should end up having length J*K
      denom <-  rep(vars., times=J.)
      numer <- rep((1:J.)*n., each=K.)
      info <- numer/denom
      tau <- mu.full.vec*sqrt(info)
      ts.true <- sweep(ts., 2, tau, "+")   # Add the above tau vector to every row in the matrix ts
      if(composite){
        ts.composite <- matrix(NA, nrow=nrow(ts.), ncol=J.)
        for(i in 1:J.){ ts.composite[,i] <- rowSums(ts.true[, (1+(i-1)*K.):(i*K.)])
        }
        return(ts.composite)
      }
      return(ts.true)
    }

  ############### Covariance matrix ######################
    #outcome_covars <- rep(rho.scalar, sum(1:(K-1)))
    # ^^^ All covariances in one vector. Begin with covariances of the first outcome,
    # i.e. p12, p13,...,p1k, then second outcome, i.e. p23, p24,...p2k, etc.
    # Currently all equal, but code allows different value.
    stage.row <- matrix(rep(1:J, each=K), J*K, J*K)
    stage.col <- t(stage.row)
    Lambda <- sqrt(pmin(stage.row, stage.col)/pmax(stage.row, stage.col))
    rho_submatrix <- matrix(1, K, K)
    rho_submatrix[which(lower.tri(rho_submatrix))] <- rho.vec
    rho_submatrix <- t(rho_submatrix)
    rho_submatrix[which(lower.tri(rho_submatrix))] <- rho.vec
    rho_matrix <- matrix(NA, J*K, J*K)
    for(j1 in 1:J){
      for(j2 in 1:J){
        rho_matrix[(1+(j1-1)*K):(j1*K), (1+(j2-1)*K):(j2*K)] <- rho_submatrix
      }
    }
    Lambda <- Lambda*rho_matrix

    # The means for the K test statistics at stage 1, 2, ..., J.
    means.typeI <- rep(0, times=J*K)
    ts <- mvtnorm::rmvnorm(nsims, mean=means.typeI, sigma = Lambda)
  # Sum the test statistics at each stage to create the composite test statistics
    ts.composite <- matrix(NA, nrow=nsims, ncol=J)
    for(i in 1:J){
      ts.composite[,i] <- rowSums(ts[,(1+(i-1)*K):(i*K)])
    }

    ############### Optimisation ###################
  # FIND THE SET OF BOUNDARIES THAT MINIMISES (probability of rejection - alpha)^2
  # NOTE: In composite function, we use method="Brent", which allows limits and one-dimensional optimisation.
    #browser()
    # No longer using bobyqa because optimisation is now 1-D.
    final.const <- minqa::bobyqa(par=1,
                        fn = pRejectFastCommonC,
                        lower=0.01,
                        upper=30,
                        J=J,
                        K=K,
                        m=m,
                        wang.delta=wang.delta,
                        ts=ts,
                        nsims=nsims,
                        alpha=alpha,
                        prob.only=TRUE
                        )$par
    # final.const <- optim(par=1,
    #                      fn=pRejectFastCommonC,
    #                      J=J,
    #                      K=K,
    #                      m=m,
    #                      wang.delta=wang.delta,
    #                      ts=ts,
    #                      nsims=nsims,
    #                      alpha=alpha,
    #                      composite=FALSE,
    #                      prob.only=TRUE,
    #                      method="Brent",
    #                      lower=0.01,
    #                      upper=10)$par
  # Use boundaries to obtain the type one error:
    typeIerr <- pRejectFastCommonC(const=final.const,
                        J=J,
                        K=K,
                        m=m,
                        wang.delta=wang.delta,
                        ts=ts,
                        nsims=nsims,
                        alpha=alpha,
                        prob.only = FALSE)

  # Now optimise and find type I error for composite outcome.
  # NOTE: Here, using Brent to optimise, because R recommends this or optimize() when optimisation is in one dimension.
    final.const.composite <- minqa::bobyqa(par=1,
                                   fn=pRejectFastCommonC,
                                   lower=0.01,
                                   upper=30,
                                   J=J,
                                   K=K,
                                   m=m,
                                   wang.delta=wang.delta,
                                   ts=ts.composite,
                                   nsims=nsims,
                                   alpha=alpha,
                                   composite=TRUE,
                                   prob.only=TRUE)$par
    # final.const.composite <- optim(par=1,
    #                                fn=pRejectFastCommonC,
    #                                J=J,
    #                                K=K,
    #                                m=m,
    #                                wang.delta=wang.delta,
    #                                ts=ts.composite,
    #                                nsims=nsims,
    #                                alpha=alpha,
    #                                composite=TRUE,
    #                                prob.only=TRUE,
    #                                method="Brent",
    #                                lower=0.01,
    #                                upper=10)$par
    typeIerr.composite <- pRejectFastCommonC(const=final.const.composite,
                                  J=J,
                                  K=K,
                                  m=m,
                                  wang.delta=wang.delta,
                                  ts=ts.composite,
                                  nsims=nsims,
                                  alpha=alpha,
                                  composite=TRUE,
                                  prob.only = FALSE)



  ############### Obtain power ############
  # Define this as the probability of rejecting the null, even by incorrectly concluding that the treatment has an effect on some outcome.
  # Can use LFC or set the number and index of working treatments in the initial arguments.
    pwr.ess <- list(0, NA, NA, NA)
    means.power <- delta0
    if(working.outs[1]=="lfc"){
      std.trt.effects <- delta1/vars
      # Choose the smallest standarised treatment effects to be the working treatment effects:
      working.outs.idx <- rank(std.trt.effects)<=m
    }else{
      working.outs.idx <- working.outs
    }
    means.power[working.outs.idx] <- delta1[working.outs.idx]
    means.power <- rep(means.power, times=J)
  # Need to multiply each outcome's trt effect by sqrt(j*n/vars), where vars is the variance of outcome K. Vector should end up having length J*K
    denom <-  rep(vars, times=J)
    i <- 0
    n.vec <- 1:maxn.stage
    while(pwr.ess[[1]]<power & i<length(n.vec)){
     # browser()
      i <- i+1
      numer <- rep((1:J)*n.vec[i], each=K)
      information <- numer/denom
      tau <- means.power*sqrt(information)
      ts.power <- sweep(ts, 2, tau, "+")   # Add the above tau vector to every row in the matrix ts
      pwr.ess <- pRejectFastCommonC(J=J,
                         K=K,
                         m=m,
                         const=final.const,
                         wang.delta=wang.delta,
                         ts=ts.power,
                         nsims=nsims,
                         alpha = alpha,
                         prob.only = FALSE)
     print(paste("Power is ", format(pwr.ess[1], digits=3), " when n per stage is ", n.vec[i], ". Expected no. of stages: ", pwr.ess$expd.no.stages, sep=""), q=F)

    }
    if(n.vec[i]==maxn.stage & pwr.ess[[1]]<power){
      warning("For multi outcome approach, desired power not achieved with current max N", call. = FALSE)
    }
    n.final <- n.vec[i]

  ############### Obtain composite power  ################
  # Now find the composite power:
    pwr.ess.composite <- list(0, NA, NA, NA)
    i <- 0
    while(pwr.ess.composite[[1]]<power & i<length(n.vec)){
      i <- i+1
      numer <- rep((1:J)*n.vec[i], each=K)
      information <- numer/denom
      tau <- means.power*sqrt(information)
      ts.power <- sweep(ts, 2, tau, "+")   # Add the above tau vector to every row in the matrix ts
      ts.power.composite <- matrix(NA, nrow=nsims, ncol=J)
      for(j in 1:J){
        ts.power.composite[,j] <- rowSums(ts.power[,(1+(j-1)*K):(j*K)])
      }
      pwr.ess.composite <- pRejectFastCommonC(J=J,
                         K=K,
                         m=m,
                         const=final.const.composite,
                         wang.delta=wang.delta,
                         ts=ts.power.composite,
                         nsims=nsims,
                         alpha = alpha,
                         composite = TRUE,
                         prob.only = FALSE)
      print(paste("Composite power is ", format(pwr.ess.composite[[1]], digits=3), " when n per stage is ", n.vec[i]), q=F)

    }
    if(n.vec[i]==maxn.stage & pwr.ess.composite[[1]]<power){
      warning("For composite approach, desired power not achieved with current max N", call. = FALSE)
    }
    n.final.composite <- n.vec[i]



  # In the case that we want both approaches to have the same N, for comparison purposes:
  # set n for both approaches to be the greater of the two values:
    n.final.max.n <- max(n.final, n.final.composite)
  # Update power and expected no. of stages:
    numer <- rep((1:J)*n.final.max.n, each=K)
    information <- numer/denom
    tau <- means.power*sqrt(information)
    ts.power.max.n <- sweep(ts, 2, tau, "+")   # Add the above tau vector to every row in the matrix ts
    ts.power.composite.max.n <- matrix(NA, nrow=nsims, ncol=J)
    for(j in 1:J){
      ts.power.composite.max.n[,j] <- rowSums(ts.power.max.n[,(1+(j-1)*K):(j*K)])
    }
    pwr.ess.max.n <- pRejectFastCommonC(J=J,
                       K=K,
                       m=m,
                       const=final.const,
                       wang.delta=wang.delta,
                       ts=ts.power.max.n,
                       nsims=nsims,
                       alpha = alpha,
                       prob.only = FALSE)
    # Composite approach:
    pwr.ess.composite.max.n <- pRejectFastCommonC(J=J,
                                 K=K,
                                 m=m,
                                 const=final.const.composite,
                                 wang.delta=wang.delta,
                                 ts=ts.power.composite.max.n,
                                 nsims=nsims,
                                 alpha = alpha,
                                 composite = TRUE,
                                 prob.only = FALSE)

    ess <- n.final*pwr.ess$expd.no.stages
    ess.composite <- n.final.composite*pwr.ess.composite$expd.no.stages

    # enm.mo <- K*n.final*pwr.ess$expd.no.stages
    # enm.composite <- K*n.final.composite*pwr.ess.composite$expd.no.stages

    ess.max.n <- n.final.max.n*pwr.ess.max.n$expd.no.stages
    ess.composite.max.n <- n.final.max.n*pwr.ess.composite.max.n$expd.no.stages
    # NOTE: ESS0 is independent of n:
    ess0 <- n.final*typeIerr$expd.no.stages
    ess0.c <- n.final.composite*typeIerr.composite$expd.no.stages


  ############### True treatment effects =/= delta0 OR delta1 #################

    createTrueTSFindpReject <- function(mu.matrix,
                                        ts.,
                                        vars,
                                        J.,
                                        K.,
                                        n.,
                                        const,
                                        m,
                                        wang.delta,
                                        nsims,
                                        alpha,
                                        prob.only=FALSE,
                                        composite=FALSE){
      results.matrix <-  matrix(NA, nrow(mu.matrix), ncol=2)
      for(i in 1:nrow(mu.matrix)){
        ts.true.current <- createTrueTS(mu.=unlist(mu.matrix[i,]),
                                        ts.=ts.,
                                        vars=vars,
                                        J.=J.,
                                        K.=K.,
                                        n.=n.,
                                        composite=composite)
      true.results <- pRejectFastCommonC(const=const,
                                         J=J.,
                                         K=K.,
                                         m=m,
                                         wang.delta=wang.delta,
                                         ts=ts.true.current,
                                         nsims=nsims,
                                         prob.only=FALSE,
                                         alpha=alpha,
                                         composite=composite)
      results.matrix[i,] <- c(true.results$prob.reject, n.*true.results$expd.no.stages)
      }
      #results.matrix <- data.frame()
      return(as.data.frame(results.matrix))
    }

    if(!is.null(delta.true)){
      true.oc <- createTrueTSFindpReject(mu.matrix=means.true,
                                       ts.=ts,
                                       vars=vars.true,
                                       J.=J,
                                       K.=K,
                                       n.=n.final,
                                       const=final.const,
                                       m=m,
                                       wang.delta=wang.delta,
                                       nsims=nsims,
                                       alpha=alpha,
                                       prob.only=FALSE,
                                       composite=FALSE)
    true.oc.composite <- createTrueTSFindpReject(mu.matrix=means.true,
                                       ts.=ts,
                                       vars=vars.true,
                                       J.=J,
                                       K.=K,
                                       n.=n.final.composite,
                                       const=final.const.composite,
                                       m=m,
                                       wang.delta=wang.delta,
                                       nsims=nsims,
                                       alpha=alpha,
                                       prob.only=FALSE,
                                       composite=TRUE)
  # Now get results for true effects for fixed N:
    true.oc.fix.n <- createTrueTSFindpReject(mu.matrix=means.true,
                                       ts.=ts,
                                       vars=vars.true,
                                       J.=J,
                                       K.=K,
                                       n.=n.final.max.n,
                                       const=final.const,
                                       m=m,
                                       wang.delta=wang.delta,
                                       nsims=nsims,
                                       alpha=alpha,
                                       prob.only=FALSE,
                                       composite=FALSE)
    true.oc.composite.fix.n <- createTrueTSFindpReject(mu.matrix=means.true,
                                                 ts.=ts,
                                                 vars=vars.true,
                                                 J.=J,
                                                 K.=K,
                                                 n.=n.final.max.n,
                                                 const=final.const.composite,
                                                 m=m,
                                                 wang.delta=wang.delta,
                                                 nsims=nsims,
                                                 alpha=alpha,
                                                 prob.only=FALSE,
                                                 composite=TRUE)

    # mo.true <- cbind(true.oc, true.oc.fix.n, delta.true, means.true)
    # mo.true$design <- "MO"
    # composite.true <- cbind(true.oc.composite, true.oc.composite.fix.n, delta.true, means.true)
    # composite.true$design <- "Composite"
    # true.results <- rbind(mo.true, composite.true)
    # names(true.results) <- c("prob.reject", "ess", "prob.reject.fix.n", "ess.fix.n", "pi.working", "pi.nonworking", paste("pi.", 1:ncol(means.true), sep=""), "design")

    mo.true <- cbind(true.oc, delta.true, means.true)
    mo.true$design <- "MO"
    mo.true$fixed.n <- "No"
    if(ncol(delta.true)==2){
      names(mo.true) <- c("prob.reject", "ess", "mu.working", "mu.nonworking", paste("mu.", 1:ncol(means.true), sep=""), "design", "fixed.n")
    }else{
      names(mo.true) <- c("prob.reject", "ess",  paste("delta.", 1:ncol(means.true), sep=""), paste("mu.", 1:ncol(means.true), sep=""), "design", "fixed.n")
    }

    composite.true <- cbind(true.oc.composite, delta.true, means.true)
    composite.true$design <- "Composite"
    composite.true$fixed.n <- "No"
    names(composite.true) <- names(mo.true)

    mo.true.fixed <- cbind(true.oc.fix.n, delta.true, means.true)
    mo.true.fixed$design <- "MO"
    mo.true.fixed$fixed.n <- "Yes"
    names(mo.true.fixed) <- names(mo.true)

    composite.true.fixed <- cbind(true.oc.composite.fix.n, delta.true, means.true)
    composite.true.fixed$design <- "Composite"
    composite.true.fixed$fixed.n <- "Yes"
    names(composite.true.fixed) <- names(mo.true)

    true.results <- rbind(mo.true, composite.true, mo.true.fixed, composite.true.fixed)
    #names(true.results) <- c("prob.reject", "ess", "pi.working", "pi.nonworking", paste("pi.", 1:ncol(means.true), sep=""), "design", "fixed.n")
    # We currently (sep 2020) only care about non-fixed results:
    #true.results <- true.results[true.results$fixed.n=="No", ]


    ratios <- mo.true$ess/composite.true$ess
    output.true.ratios <- data.frame(mo.true$prob.reject, composite.true$prob.reject, delta.true, means.true, ratios)
    if(ncol(delta.true==2)){
      names(output.true.ratios) <- c("p.reject.mo", "p.reject.comp", "mu.working", "mu.nonworking", paste("mu.", 1:ncol(means.true), sep=""), "ess.ratio")
    }else{
      names(output.true.ratios) <- c("p.reject.mo", "p.reject.comp", paste("delta.", 1:ncol(means.true), sep=""), paste("mu.", 1:ncol(means.true), sep=""), "ess.ratio")
    }

  }  # end of if statement

  #  browser()

    # Shared design characteristics:
    if(reuse.deltas==TRUE){
      des.chars <- data.frame(K, m, J, t(return.delta0), t(return.delta1), alpha, power, rho.vec[1], wang.delta)
      colnames(des.chars) <- c("K", "m", "J", "delta0.1", "delta0.2", "delta1.1", "delta1.2", "alpha", "req.power", "cor", "WangDelta")
    }else{
      des.chars <- data.frame(K, m, J, t(delta0), t(delta1), alpha, power, rho.vec[1], wang.delta) # This includes all delta0 and delta1 values (use this if specifying separate delta0/1 values for each outcome)
      colnames(des.chars) <- c("K", "m", "J", paste("delta0.k", 1:K, sep=""), paste("delta1.k", 1:K, sep=""), "alpha", "req.power", "cor", "WangDelta")
    }
    des.chars$ess0.ratio <- ess0/ess0.c
    des.chars$ess1.ratio <- ess/ess.composite

    final.n.stage <- c(n.final, n.final.composite)
    N <- J*final.n.stage
    ess0.vec <- c(ess0, ess0.c)
    ess1.vec <- c(ess, ess.composite)
    type1.vec <- c(typeIerr$prob.reject, typeIerr.composite$prob.reject)
    power.vec <- c(pwr.ess$prob.reject, pwr.ess.composite$prob.reject)
    bounds.const <- c(final.const, final.const.composite)
# for fixed n:
    N.fixed <- rep(n.final.max.n, 2)
    ess1.vec.max.n <- c(ess.max.n, ess.composite.max.n)
    power.vec.max.n <- c(pwr.ess.max.n$prob.reject, pwr.ess.composite.max.n$prob.reject)
    design.results <- cbind(final.n.stage, N, ess0.vec, ess1.vec, type1.vec, power.vec, bounds.const, N.fixed, ess1.vec.max.n, power.vec.max.n)
    colnames(design.results) <- c("n.stage", "N", "ESS0", "ESS1", "typeIerr", "power", "C", "N.fix.n", "ESS1.fix.n", "power.fix.n")
    rownames(design.results) <- c("MO", "composite")
    design.results <- as.data.frame(design.results)
    design.results$p.correct.go <- c(pwr.ess$prob.correct.go, NA)
    design.results$p.incorrect.go <- c(pwr.ess$prob.incorrect.go, NA)
    design.results$design <- c("MO", "Composite")
    to.return <- list(input=des.chars,
                      results=design.results)
    if(!is.null(delta.true)){
      to.return$true.results <- true.results
      to.return$true.ratios=output.true.ratios
    }
    # if(return.boundaries==TRUE){
    #   to.return$bounds <- boundaries
    #   to.return$bounds.c <- boundaries.composite
    # }
    if(return.ts==TRUE){
      to.return$ts <- ts
      to.return$ts.pwr <- ts.power
      to.return$ts.comp <- ts.power.composite
    }
    class(to.return) <- "moms"
    return(to.return)
  } # end of overall function
  #### End of main function ####


  #' Find Multi-Outcome Multi-Stage Trials That Allow a General Number of Efficacious Outcomes
  #'
  #' This function allows users to find a single-arm multi-outcome multi-stage trial that
  #' allows ending for trial success if promising effects are observed on a general number
  #' of outcomes, specified by the user.
  #' @import ggplot2
  #' @import gridExtra
  #' @import minqa
  #' @import Rfast
  #' @export
  findDes <- function(K=default.K,
                      m=default.m,
                      J=default.J,
                      rho.vec=default.cor,
                      nsims=default.nsims,
                      wang.delta=0,
                      alpha=default.alpha,
                      power=default.power,
                      delta0=default.delta0,
                      delta1=default.delta1,
                      reuse.deltas=TRUE,
                      delta.true=NULL,
                      reuse.true.deltas=TRUE,
                      vars.true=NULL,
                      vars=NULL,
                      working.outs=NULL,
                      maxn.stage=200,
                      return.boundaries=FALSE,
                      return.ts=FALSE
  )
  {
    recycleDeltas <- function(vec, working.outs., K.){
      full.delta.vec <- rep(vec[2], K.)
      full.delta.vec[working.outs.] <- vec[1]
      return(full.delta.vec)
    }
    #### Warnings, checks: ####
    if(is.null(delta0)){
      warning("No uninteresting treatment effects delta0 supplied. Using delta0=0, for all outcomes.", call. = FALSE)
      delta0 <- rep(0, K)
    }
    if(is.null(rho.vec)){
      warning("No correlations supplied. Using rho=0.5 for all correlations.", call. = FALSE)
      rho.vec <- rep(0.5, times=sum(1:(K-1)))
    }
    if(length(rho.vec)==1 & K>2){
      warning("Single value supplied for correlations supplied. Using this value for all correlations.", call. = FALSE)
      rho.vec <- rep(rho.vec, times=sum(1:(K-1)))
    }
    if(is.null(vars)){
      warning("No outcome variances supplied. Using vars=1 for all outcomes.", call. = FALSE)
      vars <- rep(1, K)
    }
    if(is.null(vars.true) & !is.null(delta.true)){
      warning("No TRUE outcome variances supplied. Using anticipated vars as true vars (default is 1).", call. = FALSE)
      vars.true <- vars
    }
    if(is.null(working.outs)){
      warning("Indices of working treatments not supplied. Taking indices of working treatments as treatments 1 to m.", call. = FALSE)
      working.outs <- 1:m
    }
    if(is.null(m)){
      warning("Number of outcomes required to show promise, m, not given. Using number of working treatments as treatments.", call. = FALSE)
      m <- length(working.outs)
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
        warning("reuse.deltas set to TRUE: If K>2, will take delta0[1] as delta0 for all working outcomes and delta0[2] as delta0 for all non-working outcomes. Ditto delta1")
      }

    }else{
      if(!is.null(delta.true)){
        warning("CAUTION: reuse.deltas set to FALSE, but values supplied for true delta. This could cause problems re. delta.true and means.true")
      }
    }

    if(!is.null(delta.true)){
      if(reuse.true.deltas==TRUE){
        means.true <- t(apply(delta.true, 1, recycleDeltas, working.outs.=working.outs, K.=K))
        if(K>2){
          warning("reuse.deltas set to TRUE: If K>2, will take delta.true[1] as true delta for all working outcomes and delta.true[2] as true delta for all non-working outcomes.")
        }
      }else{
        means.true <- delta.true
        if(ncol(delta.true)!=K){
          stop("The number of columns (outcomes) in delta.true does not equal K and reuse.delta.true==FALSE.")
        }
      }
    }

    # Checks:
    if(length(rho.vec)!=sum(1:(K-1))) stop("Number of rho values, i.e. length of rho.vec, should be equal to sum(1:(K-1)) ", call. = FALSE)
    if(length(delta0)!=K & !is.null(delta0)) stop ("Number of supplied uninteresting treatment effects, i.e. delta0, should be equal to K, the number of outcomes",
                                                   call. = FALSE)
    if(length(delta1)!=K) stop ("Number of supplied treatment effects, i.e. delta1, should be equal to K, the number of outcomes", call. = FALSE)
    if(length(vars)!=K) stop ("Number of outcome variances, i.e. vars, should be equal to K (the number of outcomes)", call. = FALSE)
    #### P(rejection) ####

    ##### new pRejectFast, with shared constant C:
    pRejectFastCommonC <- function(const,
                                   J=J,
                                   K=K,
                                   m=m,
                                   wang.delta,
                                   ts,
                                   nsims,
                                   prob.only=TRUE,
                                   alpha=alpha,
                                   composite=FALSE
    ){
      e.vec <-  const*((1:J)/J)^(wang.delta-0.5)
      f.vec <- -e.vec
      f.vec[length(f.vec)] <- e.vec[length(e.vec)]
      if(J==1){
        if(composite==FALSE) {
          go.overall <- colSums(t(ts) > e.vec)>=m
        }else{
          go.overall <- ts > e.vec
        }
        prob.reject <- sum(go.overall)/nsims
        minimise.prob <- abs(prob.reject - alpha)
        expd.no.stages <- 1
      } else { # ie if J!=1:
        if(composite==FALSE){
          nogo.overall <- vector("list", J)
          go.overall <- vector("list", J)
          first.m.outs.exceed.e <- vector("list", J)
          other.combn.exceed.e <- vector("list", J)
          for(j in 1:J){
            nogo.overall[[j]] <- colSums(t(ts[,(1+(j-1)*K):(j*K)]) < f.vec[j])>=(K-m+1)
            go.overall[[j]] <- colSums(t(ts[,(1+(j-1)*K):(j*K)]) > e.vec[j])>=m
            first.m.outs.exceed.e[[j]] <- colSums(t(ts[,(1+(j-1)*K):(m+(j-1)*K)]) > e.vec[j])==m
            other.combn.exceed.e[[j]] <- go.overall[[j]]==TRUE & first.m.outs.exceed.e[[j]]==FALSE
          }
          # Are m or more boundaries crossed? Row=simulation, col=stage
          nogo.trial.binary <- t(do.call(rbind, nogo.overall))
          go.trial.binary <- t(do.call(rbind, go.overall))
          first.m.outs.exceed.e <- t(do.call(rbind, first.m.outs.exceed.e))
          other.combn.exceed.e <-  t(do.call(rbind, other.combn.exceed.e))
        } else{ # if composite==TRUE
          # Only J boundaries for composite:
          nogo.trial.binary <-  t(t(ts) < f.vec)
          go.trial.binary <- t(t(ts) > e.vec)
          # browser()
        }
        # The first stage at which a nogo decision is made (and analogous for go):
        # Add extra column of 1's so that there is always some maximum even if NOGO boundary is never crossed:
        nogo.trial.binary.plus <- cbind(nogo.trial.binary, 1)
        # Add extra column of 1's so that there is always some maximum even if GO boundary is never crossed:
        go.trial.binary.plus <- cbind(go.trial.binary, 1)
        first.nogo.stage <- Rfast::rowMaxs(nogo.trial.binary.plus, value=FALSE) # Rfast
        first.go.stage <- Rfast::rowMaxs(go.trial.binary.plus, value=FALSE) # Rfast
        first.stop.stage <- cbind(first.nogo.stage, first.go.stage)
        mode(first.stop.stage) <- "numeric"
        # Does the trial make a nogo or a go decision first?
        # Final decision: 1=nogo, 2=go
        final.decision <- Rfast::rowMins(first.stop.stage, value=FALSE) # value=FALSE returns indices
        prob.reject <- sum(final.decision==2)/nsims
        minimise.prob <- abs(prob.reject - alpha)^2 # 29th Jul: added square.
        stop.stage <- Rfast::rowMins(first.stop.stage, value=TRUE) #Rfast
        expd.no.stages <- sum(stop.stage)/nsims
      } # end of if J==1 else
      if(prob.only==TRUE){
        return(minimise.prob)
      } else{
        if(composite==FALSE & J!=1){
          go.decision.index <- which(final.decision==2)
          correct.go <- rep(NA, length(go.decision.index))
          incorrect.go <- rep(NA, length(go.decision.index))
          for(i in 1:length(go.decision.index)){
            correct.go[i] <- first.m.outs.exceed.e[go.decision.index[i], stop.stage[go.decision.index[i]]]
            incorrect.go[i] <- other.combn.exceed.e[go.decision.index[i], stop.stage[go.decision.index[i]]]
          }
          prob.correct.go <- sum(correct.go)/nsims
          prob.incorrect.go <- sum(incorrect.go)/nsims
        }else{
          prob.correct.go <- NA
          prob.incorrect.go <- NA
        }
        return(list(prob.reject=prob.reject,
                    expd.no.stages=expd.no.stages,
                    f.vec=f.vec,
                    e.vec=e.vec,
                    prob.correct.go=prob.correct.go,
                    prob.incorrect.go=prob.incorrect.go)
        )
      }
    } # end of function

    trueReject <- function(
      J=J,
      K=K,
      m=m,
      ts=ts,
      nsims,
      prob.only=FALSE,
      alpha=alpha,
      means.true=means.true,
      vars.true=vars.true,
      f.vec=pwr.ess$f.vec,
      e.vec=pwr.ess$e.vec,
      n.final=n.final,
      composite=FALSE
    ){
      # Need test statistic for true effects:
      means.true <- rep(means.true, times=J)
      # Need to multiply each outcome's true trt effect by sqrt(j*n/vars), where vars is the variance of outcome K. Vector should end up having length J*K
      denom <-  rep(vars.true, times=J)
      numer <- rep((1:J)*n.final, each=K)
      information <- numer/denom
      tau.true <- means.true*sqrt(information)
      ts.true <- sweep(ts, 2, tau.true, "+")   # Add the above tau vector to every row in the matrix ts
      if(J==1){
        if(composite==FALSE) {
          go <- t(apply(ts.true, 1, function(x) x > e.vec))
          go.overall <- apply(go, 1, function(x) sum(x)>=m)
        } else{
          ts.composite <- rowSums(ts.true)
          go.overall <- ts.composite > e.vec
        }
        prob.reject <- sum(go.overall)/nsims
        ess <- n.final
      } else {
        if(composite==FALSE){
          nogo <- t(apply(ts.true, 1, function(x) x < f.vec))
          go <- t(apply(ts.true, 1, function(x) x > e.vec))
          nogo.overall <- vector("list", J)
          go.overall <- vector("list", J)
          for(j in 1:J){
            # Subset to the K outcomes for stage j:
            current.stage.nogo <- nogo[, (1+(j-1)*K):(j*K)]
            current.stage.go <- go[, (1+(j-1)*K):(j*K)]
            # Would the trial stop for either go or nogo at stage j?
            nogo.overall[[j]] <- apply(current.stage.nogo, 1, function(x) sum(x)>=(K-m+1))
            go.overall[[j]] <- apply(current.stage.go, 1, function(x) sum(x)>=m)
          }
          # Is a boundary crossed? Row=simulation, col=stage
          nogo.trial.binary <- t(do.call(rbind, nogo.overall))
          go.trial.binary <- t(do.call(rbind, go.overall))
        } else {
          # Sum the test statistics at each stage to form the composite test statistics
          ts.composite <- matrix(NA, nrow=nsims, ncol=J)
          for(i in 1:J){
            ts.composite[,i] <- rowSums(ts.true[,(1+(i-1)*K):(i*K)])
          }
          nogo.trial.binary <- t(apply(ts.composite, 1, function(x) x < f.vec)) + 0
          go.trial.binary <- t(apply(ts.composite, 1, function(x) x > e.vec)) + 0
        }
        # The first stage at which a nogo decision is made (and analogous for go):
        first.nogo.stage <- apply(nogo.trial.binary, 1, function(x) min(which(x==1), Inf))
        first.go.stage <- apply(go.trial.binary, 1, function(x) min(which(x==1), Inf))
        first.stop.stage <- cbind(first.nogo.stage, first.go.stage)
        # Does the trial make a nogo or a go decision first?
        # Final decision: 1=nogo, 2=go
        final.decision <- apply(first.stop.stage, 1, which.min)
        prob.reject <- sum(final.decision==2)/nsims
        stop.stage <- apply(first.stop.stage, 1, min)
        expd.no.stages <- sum(stop.stage)/nsims
        ess <- expd.no.stages*n.final
        #minimise.prob <- abs(prob.reject - alpha)
      } # end of J==1 else
      return(c(prob.reject, ess))
    } # end of function

    constToBounds <- function(const, J., wang.d){
      e.vec <- const*((1:J.)/J.)^(wang.d-0.5)
      f.vec <- -e.vec
      f.vec[length(f.vec)] <-e.vec[length(e.vec)]
      bounds <- list(e=e.vec, f=f.vec)
      return(bounds)
    }

    createTrueTS <- function(ts., mu., vars., J., K., n., composite=FALSE){
      mu.full.vec <- rep(mu., times=J.)
      # Need to multiply each outcome's true trt effect by sqrt(j*n/vars), where vars is the variance of outcome K. Vector should end up having length J*K
      denom <-  rep(vars., times=J.)
      numer <- rep((1:J.)*n., each=K.)
      info <- numer/denom
      tau <- mu.full.vec*sqrt(info)
      ts.true <- sweep(ts., 2, tau, "+")   # Add the above tau vector to every row in the matrix ts
      if(composite){
        ts.composite <- matrix(NA, nrow=nrow(ts.), ncol=J.)
        for(i in 1:J.){ ts.composite[,i] <- rowSums(ts.true[, (1+(i-1)*K.):(i*K.)])
        }
        return(ts.composite)
      }
      return(ts.true)
    }

    ############### Covariance matrix ######################
    #outcome_covars <- rep(rho.scalar, sum(1:(K-1)))
    # ^^^ All covariances in one vector. Begin with covariances of the first outcome,
    # i.e. p12, p13,...,p1k, then second outcome, i.e. p23, p24,...p2k, etc.
    # Currently all equal, but code allows different value.
    stage.row <- matrix(rep(1:J, each=K), J*K, J*K)
    stage.col <- t(stage.row)
    Lambda <- sqrt(pmin(stage.row, stage.col)/pmax(stage.row, stage.col))
    rho_submatrix <- matrix(1, K, K)
    rho_submatrix[which(lower.tri(rho_submatrix))] <- rho.vec
    rho_submatrix <- t(rho_submatrix)
    rho_submatrix[which(lower.tri(rho_submatrix))] <- rho.vec
    rho_matrix <- matrix(NA, J*K, J*K)
    for(j1 in 1:J){
      for(j2 in 1:J){
        rho_matrix[(1+(j1-1)*K):(j1*K), (1+(j2-1)*K):(j2*K)] <- rho_submatrix
      }
    }
    Lambda <- Lambda*rho_matrix

    # The means for the K test statistics at stage 1, 2, ..., J.
    means.typeI <- rep(0, times=J*K)
    ts <- mvtnorm::rmvnorm(nsims, mean=means.typeI, sigma = Lambda)


    ############### Optimisation ###################
    # FIND THE SET OF BOUNDARIES THAT MINIMISES (probability of rejection - alpha)^2
    # NOTE: In composite function, we use method="Brent", which allows limits and one-dimensional optimisation.
    #browser()
    # No longer using bobyqa because optimisation is now 1-D.
    final.const <- minqa::bobyqa(par=1,
                                 fn = pRejectFastCommonC,
                                 lower=0.01,
                                 upper=30,
                                 J=J,
                                 K=K,
                                 m=m,
                                 wang.delta=wang.delta,
                                 ts=ts,
                                 nsims=nsims,
                                 alpha=alpha,
                                 prob.only=TRUE
    )$par
    # Use boundaries to obtain the type one error:
    typeIerr <- pRejectFastCommonC(const=final.const,
                                   J=J,
                                   K=K,
                                   m=m,
                                   wang.delta=wang.delta,
                                   ts=ts,
                                   nsims=nsims,
                                   alpha=alpha,
                                   prob.only = FALSE)


    ############### Obtain power ############
    # Define this as the probability of rejecting the null, even by incorrectly concluding that the treatment has an effect on some outcome.
    # Can use LFC or set the number and index of working treatments in the initial arguments.
    pwr.ess <- list(0, NA, NA, NA)
    means.power <- delta0
    if(working.outs[1]=="lfc"){
      std.trt.effects <- delta1/vars
      # Choose the smallest standarised treatment effects to be the working treatment effects:
      working.outs.idx <- rank(std.trt.effects)<=m
    }else{
      working.outs.idx <- working.outs
    }
    means.power[working.outs.idx] <- delta1[working.outs.idx]
    means.power <- rep(means.power, times=J)
    # Need to multiply each outcome's trt effect by sqrt(j*n/vars), where vars is the variance of outcome K. Vector should end up having length J*K
    denom <-  rep(vars, times=J)
    i <- 0
    n.vec <- 1:maxn.stage
    while(pwr.ess[[1]]<power & i<length(n.vec)){
      # browser()
      i <- i+1
      numer <- rep((1:J)*n.vec[i], each=K)
      information <- numer/denom
      tau <- means.power*sqrt(information)
      ts.power <- sweep(ts, 2, tau, "+")   # Add the above tau vector to every row in the matrix ts
      pwr.ess <- pRejectFastCommonC(J=J,
                                    K=K,
                                    m=m,
                                    const=final.const,
                                    wang.delta=wang.delta,
                                    ts=ts.power,
                                    nsims=nsims,
                                    alpha = alpha,
                                    prob.only = FALSE)
      print(paste("Power is ", format(pwr.ess[1], digits=3), " when n per stage is ", n.vec[i], ". Expected no. of stages: ", pwr.ess$expd.no.stages, sep=""), q=F)

    }
    if(n.vec[i]==maxn.stage & pwr.ess[[1]]<power){
      warning("For multi outcome approach, desired power not achieved with current max N", call. = FALSE)
    }
    n.final <- n.vec[i]

    ess <- n.final*pwr.ess$expd.no.stages
    # enm.mo <- K*n.final*pwr.ess$expd.no.stages
    # NOTE: ESS0 is independent of n:
    ess0 <- n.final*typeIerr$expd.no.stages

    ############### True treatment effects =/= delta0 OR delta1 #################
    createTrueTSFindpReject <- function(mu.matrix,
                                        ts.,
                                        vars,
                                        J.,
                                        K.,
                                        n.,
                                        const,
                                        m,
                                        wang.delta,
                                        nsims,
                                        alpha,
                                        prob.only=FALSE,
                                        composite=FALSE){
      results.matrix <-  matrix(NA, nrow(mu.matrix), ncol=2)
      for(i in 1:nrow(mu.matrix)){
        ts.true.current <- createTrueTS(mu.=unlist(mu.matrix[i,]),
                                        ts.=ts.,
                                        vars=vars,
                                        J.=J.,
                                        K.=K.,
                                        n.=n.,
                                        composite=composite)
        true.results <- pRejectFastCommonC(const=const,
                                           J=J.,
                                           K=K.,
                                           m=m,
                                           wang.delta=wang.delta,
                                           ts=ts.true.current,
                                           nsims=nsims,
                                           prob.only=FALSE,
                                           alpha=alpha,
                                           composite=composite)
        results.matrix[i,] <- c(true.results$prob.reject, n.*true.results$expd.no.stages)
      }
      return(as.data.frame(results.matrix))
    }

    if(!is.null(delta.true)){
      true.oc <- createTrueTSFindpReject(mu.matrix=means.true,
                                         ts.=ts,
                                         vars=vars.true,
                                         J.=J,
                                         K.=K,
                                         n.=n.final,
                                         const=final.const,
                                         m=m,
                                         wang.delta=wang.delta,
                                         nsims=nsims,
                                         alpha=alpha,
                                         prob.only=FALSE,
                                         composite=FALSE)

      mo.true <- cbind(true.oc, delta.true, means.true)
      if(ncol(delta.true)==2){
        names(mo.true) <- c("prob.reject", "ess", "mu.working", "mu.nonworking", paste("mu.", 1:ncol(means.true), sep=""))
      }else{
        names(mo.true) <- c("prob.reject", "ess",  paste("delta.", 1:ncol(means.true), sep=""), paste("mu.", 1:ncol(means.true), sep=""))
      }
      true.results <- mo.true
      if(ncol(delta.true==2)){
        names(output.true.ratios) <- c("p.reject.mo", "mu.working", "mu.nonworking", paste("mu.", 1:ncol(means.true), sep=""))
      }else{
        names(output.true.ratios) <- c("p.reject.mo", paste("delta.", 1:ncol(means.true), sep=""), paste("mu.", 1:ncol(means.true), sep=""))
      }
    }  # end of if statement

    # Shared design characteristics:
    if(reuse.deltas==TRUE){
      des.chars <- data.frame(K, m, J, t(return.delta0), t(return.delta1), alpha, power, rho.vec[1], wang.delta)
      colnames(des.chars) <- c("K", "m", "J", "delta0.1", "delta0.2", "delta1.1", "delta1.2", "alpha", "req.power", "cor", "WangDelta")
    }else{
      des.chars <- data.frame(K, m, J, t(delta0), t(delta1), alpha, power, rho.vec[1], wang.delta) # This includes all delta0 and delta1 values (use this if specifying separate delta0/1 values for each outcome)
      colnames(des.chars) <- c("K", "m", "J", paste("delta0.k", 1:K, sep=""), paste("delta1.k", 1:K, sep=""), "alpha", "req.power", "cor", "WangDelta")
    }

    final.n.stage <- n.final
    N <- J*final.n.stage
    ess0.vec <- ess0
    ess1.vec <- ess
    type1.vec <- typeIerr$prob.reject
    power.vec <- pwr.ess$prob.reject
    bounds.const <- final.const
    design.results <- data.frame(final.n.stage, N, ess0.vec, ess1.vec, type1.vec, power.vec, bounds.const)
    names(design.results) <- c("n.stage", "N", "ESS0", "ESS1", "typeIerr", "power", "C")
    #design.results <- as.data.frame(design.results)
    design.results$p.correct.go <- c(pwr.ess$prob.correct.go)
    design.results$p.incorrect.go <- c(pwr.ess$prob.incorrect.go)
    to.return <- list(input=des.chars,
                      results=design.results)
    if(!is.null(delta.true)){
      to.return$true.results <- true.results
    }
    # if(return.boundaries==TRUE){
    #   to.return$bounds <- boundaries
    #   to.return$bounds.c <- boundaries.composite
    # }
    if(return.ts==TRUE){
      to.return$ts <- ts
      to.return$ts.pwr <- ts.power
    }
    class(to.return) <- "moms"
    return(to.return)
  } # end of overall function
  #### End of main function ####


  #### Wrap output obtained from foreach/parallelisation: ####
  combine.input.results <- function(output.list) {cbind(rbind(output.list$input, output.list$input),
                                                        output.list$results)}

  wrapOutput <- function(full.output, combine=TRUE){
      if(combine==TRUE){
        all.output.list <- lapply(full.output, combine.input.results)
        return.output <- do.call(rbind, all.output.list)
        return.output$km <- paste("K=", return.output$K, ", m=", return.output$m, sep="")
      }else{
        input.unbound <- lapply(full.output, "[[", "input")
        input.bind <- do.call(rbind, input.unbound)
        input.bind$km <- paste("K=", input.bind$K, ", m=", input.bind$m, sep="")
        results.unbound <- lapply(full.output, "[[", "results")
        results.bind <- do.call(rbind, results.unbound)
        return.output <- list(input=input.bind,
                              results=results.bind)
      }
    return(return.output)
  }




  ####### Plotting #######
  # Plot power: Aug 2020:
  plotPowerOld <- function(plot.df, input.df, varying.param, fix.n=TRUE){
    # plot.df is the result of wrapOutput(combine=TRUE)
    # input.df is the result of wrapOutput(combine=FALSE)$input
    # varying.param is the parameter that was varied in the results
    # For title: show all parameters:
    non.varying.params <- setdiff(names(input.df), varying.param)
    #   browser()
    plot.title.pt1 <- paste(non.varying.params[1:6], "=", input.df[1, non.varying.params[1:6]], sep="", collapse=", ")
    plot.title.pt2 <- paste(non.varying.params[7:length(non.varying.params)], "=", input.df[1, non.varying.params[7:length(non.varying.params)]] , sep="", collapse=", ")
    plot.title <- paste(plot.title.pt1, "\n", plot.title.pt2, sep="")
    #print(plot.title)
    #browser()
    pow <- if(fix.n){plot.df$power.fix.n} else {no=plot.df$power}
    p0 <- ggplot2::ggplot(data=plot.df, mapping = aes(x=plot.df[, varying.param], y=pow)) +
      # geom_hline(yintercept = plot.df$req.power[1], linetype=2)+
      geom_line(aes(col=design), size=1.2)+
      geom_point(aes(col=design), size=2)+
      labs(x=varying.param, y = 'P(reject null)' )+
      labs(title=plot.title)
    p0
  }

  # Plot power: Aug 2020:
  plotPower <- function(plot.df, varying.param, fix.n=TRUE){
    # plot.df is the result of wrapOutput(combine=TRUE)
    plot.title <- paste("Power as ", varying.param, " varies", sep="")
    #print(plot.title)
    #browser()
    pow <- if(fix.n){plot.df$power.fix.n} else {no=plot.df$power}
    p0 <- ggplot2::ggplot(data=plot.df, mapping = aes(x=plot.df[, varying.param], y=pow)) +
      # geom_hline(yintercept = plot.df$req.power[1], linetype=2)+
      geom_line(aes(col=design), size=1.2)+
      geom_point(aes(col=design), size=2)+
      labs(x=varying.param, y = 'Power' )+
      labs(title=plot.title)+
      theme_grey(base_size=22)
    print(p0)
    p0
  }



  plot.moms <- function(output){
    J <- output$chars$J
    K <- ifelse(test=output$composite=="Yes", yes=1, no=output$chars$K)
    lower <- NULL
    upper <- NULL
    for(i in 1:K){
      lower <- c(lower, c(output$bounds[[i]][1,], rep(-Inf, J)))
      upper <- c(upper, c(output$bounds[[i]][2,], rep(Inf, J)))
    }
    plot.df <- data.frame(x=rep(c(1:J, J:1), K),
                          lower=lower,
                          upper=upper,
                          outcome=rep(1:K, each=2*J)
    )
    details.df <- plot.df[!is.infinite(plot.df$lower),]
    p <- ggplot2::ggplot() +
      geom_polygon(data=plot.df, mapping=aes(x=x, y=lower), fill = "red", alpha=0.3)+
      geom_polygon(data=plot.df, mapping=aes(x=x, y=upper), fill = "lightgreen", alpha=0.5)+
      geom_line(data=details.df, mapping=aes(x=x, y=lower), size=0.5)+
      geom_line(data=details.df, mapping=aes(x=x, y=upper), size=0.5)+
      geom_point(data=details.df, mapping=aes(x=x, y=lower), size=1) +
      geom_point(data=details.df, mapping=aes(x=x, y=upper), size=1) +
      geom_text(data=details.df, mapping=aes(x=x, y=upper, label=round(upper,3)),nudge_y=0.5,show.legend = FALSE) +
      geom_text(data=details.df, mapping=aes(x=x, y=lower, label=round(lower,3)),nudge_y=-0.5,show.legend = FALSE) +
      labs(x="Stage j", y = 'Z score' ) +
      labs(title=paste("m/K: ", output$chars$m, "/", output$chars$K, "    J: ", output$chars$J, "    ESS: ", round(output$chars$ess,1),"   ",
                       "Max SS: ", output$chars$max.ss,"   ","alpha: ", round(output$chars$typeIerr,4),"   ","power: ", round(output$chars$pwr,4),sep=""))+
      scale_x_continuous(breaks=1:J)+
      facet_grid(rows=vars(outcome), labeller=label_both)
    p
  }

  plot.true <- function(output){
    K <- output$chars$K
    expd.trt.effects <- paste(output$means.power[1:K], collapse = ", ")
    plot.df <- data.frame(rbind(output$true.oc, output$true.oc.c))
    plot.df$trt.eff <- rep(output$means.true[,K], times=2)
    plot.df$design <- c(rep("multi", nrow(output$means.true)), rep("composite", nrow(output$means.true)))

    p1 <- ggplot2::ggplot(data=plot.df, mapping = aes(x=trt.eff, y=prob.reject)) +
      geom_hline(yintercept = output$chars$power, linetype=2)+
      geom_vline(xintercept = output$means.power[K], linetype=2)+
      geom_line(aes(col=design), size=2)+
      labs(x="True treatment effect of 'non-effective' treatment", y = 'P(reject null)' ) +
      labs(title=paste("m/K: ", output$chars$m, "/", output$chars$K, "    J: ", output$chars$J, "     Type I error (multi): ",
                       round(output$chars$typeIerr,4),"   ","Power (multi): ", round(output$chars$pwr,4),
                       "\n  Powered for treatment effects (", expd.trt.effects, ")", sep=""))

    p2 <-  ggplot2::ggplot(data=plot.df, mapping = aes(x=trt.eff, y=ess))+
      geom_vline(xintercept = output$means.power[K], linetype=2)+
      geom_line(aes(col=design), size=2)+
      labs(x="True treatment effect of 'non-effective' treatment", y = 'ESS' ) +
      labs(title=paste("Max N (multi): ", output$chars$max.ss, "   " ,"Max N (composite): ", output$chars$max.ss.c))
    #plots <- gridExtra::grid.arrange(p1, p2)
    #plots
    return(list(p1, p2))
  }


  plotTrue <- function(true.out, fixn){
    # Argument fixn must equal "both", "No" or "Yes".
    # Plot p(reject) as true values of delta vary:
    anticipated.delta0 <- true.out$input$delta0.2
    anticipated.delta1 <- true.out$input$delta1.1
    m <- true.out$input$m
    K <- true.out$input$K
    J <- true.out$input$J
    powered.for <- paste(c(rep(true.out$input$delta1.1, m), rep(true.out$input$delta0.2, K-m)), collapse=", ")
    title.overall <- paste("m/K: ", m, "/", K, ".  J: ", J,
                           ".\n Powered for outcome effects (", powered.for, ")", sep="", collapse=", ")
    # browser()
    if(fixn=="both"){
      plotting.df <- true.out$true.results[true.out$true.results$mu.working==anticipated.delta1, ]
      p.preject <- ggplot2::ggplot(data=plotting.df, mapping = aes(x=mu.nonworking, y=prob.reject)) +
        geom_line(aes(col=design, linetype=fixed.n), alpha=0.3, size=2)+
        scale_linetype_manual(values=c("dashed", "dotted"))+
        geom_point(aes(col=design, shape=fixed.n), size=2)+
        scale_shape_manual(values=c(2, 3))+
        geom_hline(yintercept = 0.8, linetype=2)+
        geom_vline(xintercept = true.out$input$delta0.2, linetype=2)+
        labs(x="True effect of 'non-effective' outcome", y = 'P(reject null)' )+
        labs(title=title.overall)+
        theme_grey(base_size=22)
      p.ess <- ggplot2::ggplot(data=plotting.df, mapping = aes(x=mu.nonworking, y=ess)) +
        geom_line(aes(col=design, linetype=fixed.n), alpha=0.3, size=2)+
        scale_linetype_manual(values=c("dashed", "dotted"))+
        geom_point(aes(col=design, shape=fixed.n), size=2)+
        scale_shape_manual(values=c(2, 3))+
        geom_vline(xintercept = true.out$input$delta0.2, linetype=2)+
        labs(x="True effect of 'non-effective' outcome", y = 'ESS' )+
        theme_grey(base_size=22)
    }else{
      plotting.df <- true.out$true.results[true.out$true.results$mu.working==anticipated.delta1 & true.out$true.results$fixed.n==fixn, ]
      title.overall <- paste("m/K: ", m, "/", K, ".  J: ", J,
                             ".\n Powered for treatment effects (", powered.for, ")", sep="", collapse=", ")
      p.preject <- ggplot2::ggplot(data=plotting.df, mapping = aes(x=mu.nonworking, y=prob.reject)) +
        geom_line(aes(col=design), alpha=1, size=2)+
        geom_hline(yintercept = 0.8, linetype=2)+
        geom_vline(xintercept = true.out$input$delta0.2, linetype=2)+
        labs(x="True effect of 'non-effective' outcome", y = 'P(reject null)' )+
        labs(title=title.overall)+
        theme_grey(base_size=22)
      p.ess <- ggplot2::ggplot(data=plotting.df, mapping = aes(x=mu.nonworking, y=ess)) +
        geom_line(aes(col=design), alpha=1, size=2)+
        geom_vline(xintercept = true.out$input$delta0.2, linetype=2)+
        labs(x="True effect of 'non-effective' outcome", y = 'ESS' )+
        theme_grey(base_size=22)
    }
    gridExtra::grid.arrange(p.preject, p.ess)
  }

  ####### Show how designs change as all delta0 and delta1 are varied #####
  varyDelta <- function(J=1,
                        delta0.lower,
                        delta0.upper,
                        delta1.lower,
                        delta1.upper,
                        increments,
                        delta0.leq.delta1=TRUE,
                        delta0.k,
                        delta1.k,
                        plotting=FALSE,
                        fix.n=TRUE
                        )
  {
    varied.delta1 <- seq(from=delta1.upper, to=delta1.lower, by=-increments)
    varied.delta0 <- seq(from=delta0.upper, to=delta0.lower, by=-increments)
    delta.combs <- expand.grid(varied.delta1, varied.delta0)
    names(delta.combs) <- c("d1", "d0")
    if(delta0.leq.delta1==TRUE){
      delta.combs <- delta.combs[delta.combs$d0 <= delta.combs$d1, ]
    }

    nrows <- nrow(delta.combs)
    m <- delta1.k
    K <- delta1.k + delta0.k

    output.delta.change <- vector("list", nrows)
    for(i in 1:nrows){
      current.delta0 <- rep(delta.combs$d0[i], K)
      current.delta1 <- rep(delta.combs$d1[i], K)
      output.delta.change[[i]] <- findDes(m=m,
                                          J=J,
                                          K=K,
                                          delta0=current.delta0,
                                          delta1=current.delta1,
                                          vars=rep(1, K),
                                          rho.vec =  rep(0.5, times=sum(1:(K-1))),
                                          vars.true = NULL,
                                          nsims,
                                          wang.delta,
                                          working.outs = 1:m
                                          )
    }
    output.delta.change.chars <- matrix(NA, nrow=nrows, ncol=ncol(output.delta.change[[1]]$chars))
    for(i in 1:nrows){
      output.delta.change.chars[i,] <- as.numeric(output.delta.change[[i]]$chars)
    }
    output.delta.change.chars <- as.data.frame(output.delta.change.chars)
    names(output.delta.change.chars) <- names(output.delta.change[[1]]$chars)
    output.delta.change.chars$pwr.diff <- output.delta.change.chars$pwr - output.delta.change.chars$pwr.c
    #output.delta.change.chars$varied.delta1 <- varied.delta1

    if(plotting==TRUE){
      plot1 <- plotVaryDelta(output.delta.change.chars)
      return(list(output.delta.change.chars, plot1))
    } else { # end of if statement
      return(output.delta.change.chars)
    }
  }
  # Plot the output of varyDelta:
  plotVaryDelta <- function(output){
    p6 <- ggplot2::ggplot(data=output, mapping = aes(x=delta0, y=delta1))+
      geom_tile(aes(fill = pwr.diff)) +
      scale_fill_gradient2(midpoint=0, low = 'darkred', mid="white", high = 'darkblue',
                           breaks=seq(from=-1, to=1, by=0.5),
                           labels=seq(from=-1, to=1, by=0.5),
                           limits=c(-1,1))+
      labs(title=paste("Power (multi) - power (composite) \n m/K: ", output$m[1], "/", output$K[1], ";    J: ", output$J[1], ";     alpha: ",
                       output$alpha[1], ".\nNumber of effective, ineffective outcomes: ", output$m[1], ", ",  output$K[1]- output$m[1], sep=""))+
      scale_x_continuous(breaks=sort(unique(output$delta0)))+
      scale_y_continuous(breaks=sort(unique(output$delta1)))+
      theme_tufte()
    p6
  }
  #dev.print(pdf, file="varying_delta_pt2.pdf")

  # Plot ESS ratio using tidied output, for any parameter:
  plotESS <- function(tidied.output, param, xaxis, method, no.legend=FALSE){
    param <- ensym(param)
    if(method=="mo")  {
      yaxis.text <- expression(paste("(", ESS[MO],"/", ESS[comp], ")",  " | LFC"))
      pl <- ggplot2::ggplot(data=tidied.output$input,  mapping=aes(x=!!param, y=ess1.ratio, col=km))
    }
    if(method=="dtl")  {
      yaxis.text <- expression(paste("(", ESS[DtL],"/", ESS[single], ")",  " | LFC"))
      pl <- ggplot2::ggplot(data=tidied.output$input,  mapping=aes(x=!!param, y=ess1.ratio, col=km, linetype=max.s2.prop))+
        labs(linetype="Kmax")
    }
    pl <- pl+
      geom_line(size=1)+
      geom_hline(yintercept=1,
                 linetype="dashed")+
      labs(title="ESS ratio under LFC",
           col="K, m",
           y=yaxis.text,
           x=xaxis)
    if(no.legend==TRUE){
      pl <- pl + theme(legend.position = "none")
    }
    pl
  }

  #### Plots for changing true effects ####
  plotTrueESSratio <- function(raw.output, method){
   plot1 <- ggplot2::ggplot(raw.output$true.ratios, aes(x=mu.1, y=mu.2))+
     geom_raster(aes(fill = ess.ratio))+
     scale_x_continuous(breaks=sort(unique(raw.output$true.ratios$mu.1)))+
     scale_y_continuous(breaks=sort(unique(raw.output$true.ratios$mu.2)))+
     labs(#subtitle=bquote( mu[1] == .(raw.output$input$delta1.1)*","~ mu[2] == .(raw.output$input$delta0.2)),
          y=expression(paste(mu[2])),
          x=expression(paste(mu[1])))+
     geom_text(aes(label = round(ess.ratio, 2)), size=5) +
     scale_fill_gradient2(midpoint=1, low = "darkred", mid="white", high = "darkblue")
   if(method=="mo"){
     plot1 <- plot1 +
       labs(title=bquote("ESS"[MO]*"/"*"ESS"[comp]*", powered for "*mu[1] == .(raw.output$input$delta1.1)*","~ mu[2] == .(raw.output$input$delta0.2)),
            fill="ESS ratio")
   }
   if(method=="dtl"){
     plot1 <- plot1 +
       labs(title=expression(paste(ESS[DtL],"/",ESS[single], ", powered for ")),
            fill="ESS ratio")
   }
   plot1
  }


  plotTruePrejectOLD <- function(raw.output, method){
    if(method=="mo"){
      new.approach.df <- raw.output$true.results[raw.output$true.results$design=="MO" & raw.output$true.results$fixed.n=="No", ]
      old.approach.df <- raw.output$true.results[raw.output$true.results$design=="Composite" & raw.output$true.results$fixed.n=="No", ]
    }
    if(method=="dtl"){ # XXX CHECK these labels ####
      new.approach.df <- raw.output$true.results[raw.output$true.results$design=="DtL" & raw.output$true.results$fixed.n=="No", ]
      old.approach.df <- raw.output$true.results[raw.output$true.results$design=="Single stage" & raw.output$true.results$fixed.n=="No", ]
    }
    plot1 <- ggplot2::ggplot(new.approach.df, aes(x=mu.1, y=mu.2))+
      geom_raster(aes(fill = prob.reject))+
      scale_x_continuous(breaks=sort(unique(new.approach.df$mu.1)))+
      scale_y_continuous(breaks=sort(unique(new.approach.df$mu.2)))+
      labs(subtitle=bquote( mu[1] == .(raw.output$input$delta1.1) ~~~ mu[2] == .(raw.output$input$delta0.2)),
           y=expression(paste(mu[2])),
           x=expression(paste(mu[1])))+
      geom_text(aes(label = round(prob.reject, 2)), size=5) +
      scale_fill_gradient(low = "white", high = "darkred")

    plot2 <- ggplot2::ggplot(old.approach.df, aes(x=mu.1, y=mu.2))+
      geom_raster(aes(fill = prob.reject))+
      scale_x_continuous(breaks=sort(unique(old.approach.df$mu.1)))+
      scale_y_continuous(breaks=sort(unique(old.approach.df$mu.2)))+
      labs(subtitle=bquote( mu[1] == .(raw.output$input$delta1.1) ~~~ mu[2] == .(raw.output$input$delta0.2)),
           y=expression(paste(mu[2])),
           x=expression(paste(mu[1])))+
      geom_text(aes(label = round(prob.reject, 2)), size=5) +
      scale_fill_gradient(low = "white", high = "darkred")
    if(method=="mo"){
      plot1 <- plot1 +
        # labs(title=expression(paste(P, "(reject ", H[0], ")", [MO], ", powered for ")),
        #   fill=expression(paste(P, "(reject ", H[0], ")")))
        labs(title=expression(paste("R(", mu, ")"[MO], " powered for")),
             fill=expression(paste("R(", mu, ")")))

      plot2 <- plot2 +
        labs(title=expression(paste("R(", mu, ")"[comp], " powered for")),
             fill=expression(paste("R(", mu, ")")))
    }
    if(method=="dtl"){
      plot1 <- plot1 +
        # labs(title=expression(paste(P, "(reject ", H[0], ")", [MO], ", powered for ")),
        #   fill=expression(paste(P, "(reject ", H[0], ")")))
        labs(title=expression(paste("R(", mu, ")"[DtL], " powered for")),
             fill=expression(paste("R(", mu, ")")))
      plot2 <- plot2 +
        labs(title=expression(paste("R(", H[0], ")"[single], " powered for")),
             fill=expression(paste("R(", H[0], ")")))
    }
    both.plots <- gridExtra::grid.arrange(plot1, plot2, ncol=1)
    both.plots
  }


  plotTruePreject <- function(raw.output, method){
    if(method=="mo"){
      new.approach.df <- raw.output$true.results[raw.output$true.results$design=="MO" & raw.output$true.results$fixed.n=="No", ]
      old.approach.df <- raw.output$true.results[raw.output$true.results$design=="Composite" & raw.output$true.results$fixed.n=="No", ]
    }
    if(method=="dtl"){ # XXX CHECK these labels ####
      new.approach.df <- raw.output$true.results[raw.output$true.results$design=="DtL" & raw.output$true.results$fixed.n=="No", ]
      old.approach.df <- raw.output$true.results[raw.output$true.results$design=="Single stage" & raw.output$true.results$fixed.n=="No", ]
    }
    plot1 <- ggplot2::ggplot(new.approach.df, aes(x=mu.1, y=mu.2))+
      geom_raster(aes(fill = prob.reject))+
      scale_x_continuous(breaks=sort(unique(new.approach.df$mu.1)))+
      scale_y_continuous(breaks=sort(unique(new.approach.df$mu.2)))+
      labs(y=expression(paste(mu[2])),
           x=expression(paste(mu[1])))+
      geom_text(aes(label = round(prob.reject, 2)), size=5) +
      scale_fill_gradient(low = "white", high = "grey40")+
      coord_cartesian(expand = 0) +
      theme(legend.position = "none")

    plot2 <- ggplot2::ggplot(old.approach.df, aes(x=mu.1, y=mu.2))+
      geom_raster(aes(fill = prob.reject))+
      scale_x_continuous(breaks=sort(unique(old.approach.df$mu.1)))+
      scale_y_continuous(breaks=sort(unique(old.approach.df$mu.2)))+
      labs(y=expression(paste(mu[2])),
           x=expression(paste(mu[1])))+
      geom_text(aes(label = round(prob.reject, 2)), size=5) +
      scale_fill_gradient(low = "white", high = "grey40")+
      coord_cartesian(expand = 0)
    if(method=="mo"){
      plot1 <- plot1 +
        # labs(title=expression(paste(P, "(reject ", H[0], ")", [MO], ", powered for ")),
        #   fill=expression(paste(P, "(reject ", H[0], ")")))
        labs(title=bquote("R("*mu[1] == .(raw.output$input$delta1.1)*","~ mu[2] == .(raw.output$input$delta0.2)*")"[MO]),
             fill=expression(paste("R(", mu, ")")))

      plot2 <- plot2 +
        labs(title=bquote("R("*mu[1] == .(raw.output$input$delta1.1)*","~ mu[2] == .(raw.output$input$delta0.2)*")"[comp]),
             fill=expression(paste("R(", mu, ")")))
    }
    if(method=="dtl"){
      plot1 <- plot1 +
        # labs(title=expression(paste(P, "(reject ", H[0], ")", [MO], ", powered for ")),
        #   fill=expression(paste(P, "(reject ", H[0], ")")))
        labs(title=bquote("R("*mu[1] == .(raw.output$input$delta1.1)*","~ mu[2] == .(raw.output$input$delta0.2)*")"[DtL]),
             fill=expression(paste("R(", mu, ")")))
      plot2 <- plot2 +
        labs(title=bquote("R("*mu[1] == .(raw.output$input$delta1.1)*","~ mu[2] == .(raw.output$input$delta0.2)*")"[single]),
             fill=expression(paste("R(", mu, ")")))
    }
    #browser()
    #leg <- get_legend(plot1)
    both.plots <- gridExtra::grid.arrange(plot1, plot2, widths=c(2.4,3))
    both.plots
  }
########### Create subset of true results  #########
  # subset to just the true values of interest:
  createSubset <- function(raw.output, delta.matrix){
    index <- apply(delta.matrix, 1, function(x) which(abs(raw.output$true.ratios$mu.working-x[1])<0.0001 & abs(raw.output$true.ratios$mu.nonworking-x[2])<0.0001))
    df.subset <-   raw.output$true.ratios[index, ]
    df.subset
  }


  createSearchMatrix <- function(input.df=K.Kmax.m.df, parameter.values){
    expanded.input.df <- as.matrix(input.df) %x% rep(1, length(parameter.values))
    search.matrix <- data.frame(expanded.input.df, rep(parameter.values, nrow(input.df)))
    names(search.matrix) <- c(names(input.df), deparse(substitute(parameter.values)))
    search.matrix
  }




#### Plot diagrams: ####
makeOutcomeCoords <- function(...){
  # input: any number of vectors. Each vector: y coords for one outcome.
  y.vec <- c(...)
  y.list <- list(...)
  lengths <- sapply(y.list, length)
  x.list <- lapply(y.list, function(x) 1:length(x))
  outcome <- factor(rep(1:length(y.list), times=lengths))
  line.df <- data.frame(x.line=unlist(x.list),
                        y.line=y.vec,
                        Outcome=outcome)
  line.df
}

#' Finds Wang & Tsiatis Stopping Boundaries
#'
#' This function finds stopping boundaries using the formula of Wang & Tsiatis
#' @param C Constant
#' @param delta Parameter determining shape of stopping regions.
#' @param J Number of stages
#' @export
#' @return Returns a data frame of J rows and two columns, representing the lower and upper stopping boundaries for each of the J stages.
#'
findWangTsiatisBounds <- function(C, delta, J){
    j <- 1:J
    upper <- C*(j^(delta-0.5))
    lower <- -upper
    lower[J] <- upper[J]
    bounds <- data.frame(lower=lower, upper=upper)
    bounds
  }

createWTplottingBounds <- function(C, J, delta, ymin, ymax){
  wang <- findWangTsiatisBounds(C=C, J=J, delta=delta)
  x <- c(1:J, J:1, 1:J, J:1)
  y.lower <- c(wang$lower, rep(ymin, J))
  y.upper <- c(wang$upper, rep(ymax, J))
  y <- c(y.lower, y.upper)
  label <- rep(c("lower", "upper"), each=2*J)
  polygon.df <- data.frame(x.coords=x,
                          y.coords=y,
                          label=label)
  bounds.only.y <- unlist(wang)
  bounds.only.x <- rep(1:J, times=2)
  bounds.labels <- as.character(round(bounds.only.y,2))
  bounds.df <- data.frame(x.coords=bounds.only.x,
                          y.coords=bounds.only.y,
                          labels.rounded=bounds.labels,
                          stringsAsFactors = FALSE)
  output.lists <- list(polygon.df=polygon.df,
                       bounds.df=bounds.df)
  output.lists
}

#' Plots Stopping regions and Shows Stopping Boundaries
#'
#' This function plots the upper and lower stopping regions and labels the stopping boundaries,
#' for a design realisation found using the function findDes.
#'
#' @param find.des.output The output from a call to findDes.
#' @export
#'
plotBounds <- function(find.des.output,
                       xlabel="Stage",
                       ylabel="Test Statistic",
                       title.main=NULL,
                       line.df=NULL,
                       ymin=-5,
                       ymax=5){
  WT.plotting.bounds <- createWTplottingBounds(C=find.des.output$results$C,
                                               J=find.des.output$input$J,
                                               delta=find.des.output$input$WangDelta,
                                               ymin=ymin,
                                               ymax=ymax)
  bounds.plot <- ggplot2::ggplot()+
    geom_blank()+
    geom_polygon(data=WT.plotting.bounds$polygon.df,
                 mapping=aes(x=x.coords, y=y.coords, fill=label),
                 alpha=0.2)+
    geom_point(data=WT.plotting.bounds$bounds.df,
               mapping = aes(x=x.coords, y=y.coords))+
    geom_text(data=WT.plotting.bounds$bounds.df, aes(x=x.coords, y=y.coords, label=labels.rounded),
              hjust=1.25,
              vjust=1.25)+
    ylab(ylabel)+
    scale_x_continuous(name=xlabel,
                       breaks=sort(unique(WT.plotting.bounds$polygon.df$x.coords)),
                       labels=waiver())+
    guides(fill=FALSE)
  if(!is.null(line.df)){
    bounds.plot <- bounds.plot+
      geom_line(data=line.df,
                mapping=aes(x=x.line, y=y.line, col=Outcome),
                size=1)+
      geom_point(data=line.df,
                 mapping=aes(x=x.line, y=y.line, col=Outcome),
                 size=3)
  }
  if(!is.null(title.main)){
    bounds.plot <- bounds.plot+
      labs(title = title.main)
  }
  bounds.plot
}

# Plot rejection coordinates: ####
# (comparing MO and composite approaches)

# Create rejection region coords:
createRejectionCoords <- function(C, xymax.prop=1.1){
  xymax <- xymax.prop*C[2]
  reject.mo <- rbind(c(xymax, 0),
                     c(C[1], 0),
                     c(C[1], C[1]),
                     c(0, C[1]),
                     c(0, xymax),
                     c(xymax, xymax))
  reject.mo <- data.frame(reject.mo)
  reject.mo$Design <- "MO"

  reject.comp <- rbind(c(xymax, 0),
                       c(C[2], 0),
                       c(0, C[2]),
                       c(0, xymax),
                       c(xymax, xymax))
  reject.comp <- data.frame(reject.comp)
  reject.comp$Design <- "Comp"
  reject.region.df <- rbind(reject.mo, reject.comp)
  names(reject.region.df)[1:2] <- c("x", "y")

  comp.line.df <- data.frame(x=c(C[2], 0),
                             y=c(0, C[2]))

  output <- list(reject.region.df=reject.region.df,
                 comp.line.df=comp.line.df,
                 xymax=xymax,
                 C=C)
  output
}


plotRejectionCoords <- function(rejection.list){
  # rejection.list should be a list object returned by the function createRejectionCoords()
  theme_set(theme_bw(base_size = 11)) # Increase font size and set theme for plots.
  p <- ggplot2::ggplot()+
    geom_hline(yintercept = rejection.list$C[1], linetype=2)+
    geom_vline(xintercept = rejection.list$C[1], linetype=2)+
    geom_path(data=rejection.list$comp.line.df, aes(x=x, y=y), linetype="dashed")+
    geom_polygon(data=rejection.list$reject.region.df, aes(x=x, y=y, fill=Design), alpha=0.2)+
    coord_cartesian(xlim=c(0, rejection.list$xymax),
                    ylim=c(0, rejection.list$xymax),
                    expand=0)+
    labs(title="Rejection regions",
         x=expression(paste(Z[11])),
         y=expression(paste(Z[12])))+
    scale_x_continuous(breaks=sort(c(0:3, rejection.list$C)),
                       labels=c("0", "1", "2", expression(paste(e[1])), "3", expression(paste(e[1]^(c)))),
                       minor_breaks=1:4)+
    scale_y_continuous(breaks=sort(c(0:3, rejection.list$C)),
                       labels=c("0", "1", "2", expression(paste(e[1])), "3", expression(paste(e[1]^(c)))),
                       minor_breaks=1:4)
  p
}

# Not an original fn:
get_legend<-function(myggplot){
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}




# Quickly compare bounds for MO design vs composite design
findBounds <- function(bounds, sims=100000, K=3, m=2, return.preject=FALSE){
  mat <- matrix(rnorm(K*sims), ncol=K)
  # p(reject) MO design:
  preject.mo <- sum(rowSums(mat>bounds[1])>=m)/sims
  minimised.mo <- (preject.mo-0.05)^2
  # p(reject) composite design:
  preject.comp <- sum(rowSums(mat) > bounds[2])/sims
  minimised.comp <- (preject.comp-0.05)^2
  minimised.both <- minimised.mo+minimised.comp
  if(return.preject==FALSE){
    return(minimised.both)
  }else{
    return(list(mo=preject.mo,
                comp=preject.comp))
  }
}
# # find bounds:
# bounds <- optim(par=c(1,1), fn=findBounds)$par
# findBounds(bounds=bounds, return.preject = TRUE)
# bounds[1]
# bounds[1]*sqrt(3)
# bounds[2]
