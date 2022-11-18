# © Copyright University of Exeter. This code is open source and released under the MIT Licence

pca = function(Y,W = NULL){
  # computes pca and gives weights and basis vectors of all design points
  # INPUTS:
  # Y: noncentred ensemble from which to form basis and weights
  # W: weights of each output. If NULL then assumed weight is uniform
  
  
  if(is.null(W)){
    # conduct pca using specified number of principal components. this 
    # approximates the full output as a linear combination of no.pcs basis
    # vectors
    pca_ = prcomp(Y,center = TRUE)
    
    basis.vectors = pca_$rotation
    weights = Y%*%pca_$rotation
    VarExpl = pca_$sdev^2/sum(pca_$sdev^2)
  } else{
    Winv = GetInverse(W)
    
    # weighted svd
    pca_ = MakeDataBasis(data = t(Y),weightinv = Winv,RemoveMean = TRUE)
    basis.vectors = pca_$tBasis
    weights = CalcScores(data = t(Y),basis = basis.vectors,
                         weightinv = Winv)
    VarExpl = rep(0,ncol(basis.vectors))
    for(i in 1:length(VarExpl)){
      VarExpl[i] = VarExplained(basis.vectors[,i],t(Y),weightinv = Winv)
    }
  }
  return(list(pca = pca_,weights = weights,basis.vectors = basis.vectors,
              Y = Y,VarExpl = sort(VarExpl,decreasing = TRUE)))
}

salter_rotation_basis = function(basis.vectors,Y,z,W,kmax=5,
                                 v = c(rep(0.1,5)),MaxTime=60,
                                 no.change = FALSE){
  # Takes basis vectors and uses the ensemble data Y to rotate the space 
  # ensuring that they span the observations
  
  # INPUTS:
  # basis.vectors: min(n,q) x min(n,q) matrix.
  # Y: n x q matrix.  Output data from ensemble
  # z: q-length vector.
  # W: q x q diag. matrix. Each element is addition of variances of model_disc 
  # and obs_err for each output
  # no.change: if reconstruction error is small on original basis then by 
  # setting no.change = TRUE, this puts it in a 'salter_RB' format to be used 
  # later
  
  # OUTPUTS:
  # rot.basis.vectors: min(n,q) x min(n,q) matrix. Completed rotation of 'basis.vectors'.
  # rot.weights: n x min(n,q). Weights associated with 'rot.basis.vectors'.
  # recon_err: min(n,q)-length vector. Reconstruction error of the observation 
  # from the rotated basis vectors. This is so one can choose number of BFs to 
  # get desired accuracy
  
  # reformat my basis vectors so it's compatible with framework below
  EnsembleMean = apply(Y, 2, mean)
  CentredField = t(Y) - EnsembleMean
  eigW = eigen(W)
  Q = eigW$vectors
  Lambda = 1/eigW$values
  
  DataBasis = list(tBasis = basis.vectors,CentredField = CentredField,
                   EnsembleMean = EnsembleMean,Q = Q,Lambda = Lambda)
  
  Winv = GetInverse(W)
  
  if(no.change == TRUE){
    RotatedCoeffs <- CalcScores(data = DataBasis$CentredField, 
                                basis = DataBasis$tBasis)#, weightinv = Winv)
    
    list1 = list(rot.weights = RotatedCoeffs)
    
    DataBasis = append(DataBasis,list1)
    
    return(DataBasis)
  }
  
  obsc = z - DataBasis$EnsembleMean
  
  # Rotate!
  # Can set MaxTime low here as have small ensemble and l = 100, although MaxTime=30 or 60 usually fine
  RotatedBasis <- RotateBasis(DataBasis = DataBasis, obs = obsc, 
                              weightinv = Winv, kmax = kmax, v = v, 
                              MaxTime = MaxTime, random.init = random.init)
  
  # how close is reconstruction of observation to observation as basis vectors are added (MSE)
  recon_err = vector('numeric',dim(RotatedBasis$tBasis)[2])
  for(q2rot in 1:length(recon_err)){
    ObsReconRot = DataBasis$EnsembleMean + ReconObs(obsc, RotatedBasis$tBasis[,1:q2rot], weightinv = Winv)
    recon_err[q2rot] = sum((z - ObsReconRot)^2)
  }
  
  # Then emulate each set of coeffs (each column of RotatedCoeffs)
  RotatedCoeffs <- CalcScores(data = RotatedBasis$CentredField, basis = RotatedBasis$tBasis, weightinv = Winv)
  
  return(list(tBasis = RotatedBasis$tBasis,
              CentredField = RotatedBasis$CentredField,
              EnsembleMean = RotatedBasis$EnsembleMean,
              scaling = RotatedBasis$scaling,RW = RotatedBasis$RW,
              VarExp = RotatedBasis$VarExp,Opt = RotatedBasis$Opt,
              rot.weights = RotatedCoeffs,recon_err = recon_err))
}

emulation_weights = function(X,Y,setseed = 1,HM = FALSE,
                        prev.wave = NULL,size_NROY = 10^6,var.names=var.names,
                        upper.range=upper.range,lower.range=lower.range){
  # function to build emulator, calculate mean and variance of points still in NROY
  
  # q: no. of princ. comps; n: no. design points; p = no. inputs
  
  # INPUTS:
  # X: n x p matrix
  # Y: n x q' matrix
  # setseed: to ensure reproducibility. can be changed if required
  # HM: if TRUE then allows user to take previous wave of history matching and 
  # evaluate all of NROY space in emulators built
  # var.names: vector of length p with names of all variables
  # upper.range: vector of length p with upper range of all input variables. 
  # Used to discritise input space if prev.wave == NULL
  # lower.range: vector of length p with upper range of all input variables. 
  # Used to discritise input space if prev.wave == NULL
  # prev.wave: default 'NULL'. Output of previous wave of history matching. If 
  # prev. wave has been conducted, input it so NROY can be updated from previous
  # wave
  # size_NROY: if HM is TRUE and prev.wave = NULL, how many discritised points
  # should the input space be split into?
  
  # OUTPUTS:
  # emulators
  # mean and var of points in NROY if HM = TRUE
  
  # bookkeeping. If only one output then need to convert to a matrix
  if(is.null(dim(Y)) | is.null(dim(X))){
    stop('X or Y must not be in vector form')
  }
  
  # checks that things are the sizes they should be
  stopifnot(nrow(X)==nrow(Y))
  
  q = ncol(Y)
  n = nrow(Y)
  p = ncol(X)
  
  # build emulators for each output
  emulators = lapply(1:q,function(i){
    set.seed(setseed)
    emulator = km(design = X,response = Y[,i],covtype = 'matern5_2',
                  multistart = 3,
                  control = list(trace=FALSE, pop.size=200, max.generations=50))
    return(emulator)
  })
  
  
  if(HM){
    size_input_space = size_NROY # number of points used to discritise input space
    
    # now to find mean and variance of points in NROY space
    if(is.null(prev.wave)){
      stopifnot(length(var.names)==length(lower.range),
                length(lower.range)==length(upper.range),
                length(upper.range)==ncol(X))
      
      # generate huge lhs over entire input space to estimate NROY space % if this is first wave
      X = matrix(0,size_input_space,p + 2*q)
      big_sample = randomLHS(n = size_input_space,k = p)
      X[,1:p] = t(t(big_sample)*(upper.range-lower.range) + lower.range)
      colnames(X) = c(var.names,c(rbind(paste0('mean',1:q),paste0('var',1:q))))
      
      NROY = 1:size_input_space # the whole space is in NROY
    } else{
      # if not then assume old NROY space from previous wave
      X = matrix(0,size_input_space,p + 2*q)
      X[,1:p] = prev.wave$X[,1:p]
      colnames(X) = c(var.names,c(rbind(paste0('mean',1:q),paste0('var',1:q))))
      NROY = prev.wave$NROY # NROY space determined from previous wave
    }
    
    # evaluate mean and variance of each weight for each of the inputs in X
    for(i in 1:q){
      em_pred = predict.km(emulators[[i]],newdata = X[NROY,1:p],
                           type = 'UK',checkNames = FALSE,light.return = TRUE)
      X[NROY,p + 2*i - 1] = em_pred$mean
      X[NROY,p + 2*i] = em_pred$sd^2
    }
    return(list(X = X,emulators = emulators,NROY.size = size_input_space))
  }
  return(emulators = emulators)
}

predict_outputs = function(pred.inputs,emulators,full_output = FALSE,
                           basis.vectors = NULL){
  # function to predict weights and/or full outputs given input values and the
  # emulators
  
  # INPUTS:
  # pred.inputs: n x p matrix
  # emulators: output from vector emulation_weights
  # full_output: want to emulate full output or just weights? FALSE means just 
  # the latter
  # basis.vectors: if full_output = TRUE, need basis vectors to reconstruct
  # full output
  
  # OUTPUTS:
  # output.pred: predicted outputs, whether it be the weights or full outputs
  
  no.emulators = nrow(summary(emulators))
  
  # initialise matrix for the weights
  pred.weights = matrix(0,nrow = nrow(pred.inputs),ncol = no.emulators*2)
  colnames(pred.weights) = c(rbind(paste0('PC.mean',1:no.emulators),
                                              paste0('PC.sd',1:no.emulators)))
  
  # for each emulator, get predictive mean and variance for weights
  for(i in 1:no.emulators){
    prediction = predict.km(object = emulators[[i]],newdata = pred.inputs,
                            type = 'UK',se.compute = TRUE,light.return = TRUE,
                            checkNames = FALSE)
    
    pred.weights[,2*i - 1] = prediction$mean
    pred.weights[,2*i] = prediction$sd
  }
  
  # return weights if specified
  if(full_output == FALSE){
    return(output = pred.weights)
  } else{ # reconstruct full output from weights and basis vectors
    full.output.mean = pred.weights[,seq(1,no.emulators*2,by=2)] %*% t(basis.vectors)
    full.output.sd = pred.weights[,seq(2,no.emulators*2,by=2)] %*% t(basis.vectors)
    
    return(list(output.mean = full.output.mean,output.sd = full.output.sd))
  }
}

prediction_curves = function(pred.input,emulators,basis.vectors,no.draws){
  # function to produce no.draws curves based on mean and variance from predicting pred.input
  
  # INPUTS:
  # pred.input: 1 x p matrix from which to produce output curves
  # emulators: emulators to evaluate pred.input and produce weights
  # basis.vectors: basis vectors in which to construct full output
  # no.draws: number of draws from each emulator
  
  # OUTPUTS:
  # output.draws: matrix of size no.draws x q with each draw containing one 
  # realisation of the output
  
  no.emulators = nrow(summary(emulators))
  
  # make multiple draws of gp to obtain the weights
  prediction_weights = matrix(0,no.draws,no.emulators)
  for(j in 1:no.emulators){
    set.seed(1)
    single_draw_data = predict.km(emulators[[j]],newdata = pred.input,
                                  type = 'UK',cov.compute = TRUE,
                                  checkNames = FALSE)
    
    prediction_weights[,j] = mvrnorm(no.draws,mu = single_draw_data$mean,
                                     Sigma = single_draw_data$cov)
  }
  
  # evaluate full output from drawn weights
  output.draws = prediction_weights %*% t(basis.vectors)
  
  return(output.draws)
}

hmwave = function(wave.no,X,Y,new.pca = FALSE,no.pcs = NULL,
                  pca = NULL,z,obs.err,model.disc,
                  var.names = var.names,upper.range = upper.range,
                  lower.range = lower.range,prev.wave = NULL,size_NROY = 10^6,
                  setseed = 1,no.new.samples = 10*length(var.names)){
  # completes one wave of history matching.
  # INPUTS:
  # wave.no: what number wave. if >1 need to give prev.wave
  # X: inputs of the design points. matrix of size n x p
  # Y: outputs of design points. matrix of size n x q
  # pca: 'salter_RB' object where weights are predicted from. Required 
  # if new.pca == FALSE.
  # no.pcs: number of principal components being used
  # new.pca: conduct PCA again for this wave? Recommended TRUE if dealing with 
  # model with high variability
  # no.pcs: number of weights/basis vectors to use
  # z: observation obvs
  # obs.err: same length as 'z'. sd of error of observation
  # model.disc: same length as 'z'. sd of model discrepency
  # var.names: names of input variables
  # upper/lower.range: same length of 'var.names'. upper and lower range of each
  # input variable
  # prev.wave: 'hmwave' object.
  # size_NROY: if wave.no = 1, how many discritised points should the input 
  # space be split into?
  # setseed: default 1 for reproducibility. after hm wave, sample from it and 
  # put in output of func. therefore it can easily be sampled externally
  # no.new.samples: how many new samples to generate from NROY space after hm
  # wave? default = 10*p
  #
  # OUTPUTS:
  # HM: 'HistoryMatch' object
  # emulation: 'emulation_weights' object
  # new.samples: new samples generated from NROY space to run externally and 
  # input into next wave
  # NROY: which points in NROY remain?
  # NROY.prop: proportion of NROY space remaining
  
  # safety check: if new pca isn't calculated then need pca and no.pcs
  stopifnot((new.pca == TRUE) | ((!is.null(pca) & !is.null(no.pcs))))
  
  # dimension of input
  p = ncol(X)
  
  # dimension of output
  k = length(z)
  
  # no.design.points
  n = nrow(X)
  
  W = diag(obs.err + model.disc,k)
  
  # do we want to conduct another PCA given the new outputs?
  if(new.pca == TRUE){
    pca = pca(Y)
    rotated.pca = salter_rotation_basis(basis.vectors = pca$basis.vectors,
                                        Y = Y,z = z,W = W)
    pca_weights = rotated.pca$rot.weights
    print(rotated.pca$VarExp)
    print(rotated.pca$recon_err)
    no.pcs = as.numeric(unlist(strsplit(readline('How many pcs? '),',')))
  } else{
    # calculate weights using the basis given in 'rotated.pca' and 'Y'
    pca_weights = Reduce('rbind',
                          lapply(1:n,function(i){
                           CalcScores(data = Y[i,],
                                      basis = rotated.pca$tBasis[,1:no.pcs],
                                      weightinv = GetInverse(W))
                         }))
    
  }
  
  # use the given training data to build emulators for each pc, generate NROY 
  # space and give pred. mean and var. for each point in NROY
  emulated.pca = emulation_weights(X = X,
                                   Y = matrix(pca_weights[,1:no.pcs],ncol = no.pcs),
                                   var.names = var.names,
                                   upper.range = upper.range,
                                   lower.range = lower.range,
                                   prev.wave = prev.wave,size_NROY = size_NROY,
                                   setseed = setseed,HM = TRUE)
  
  # define NROY in case of prev.wave or wave 1
  if(wave.no == 1){
    NROY = 1:emulated.pca$NROY.size
  } else{
    NROY = prev.wave$NROY
  }
  
  # use the predictions for mean and var. to HM to z
  hm.pca = HistoryMatch(DataBasis = rotated.pca,
                        Obs = z,
                        Expectation = emulated.pca$X[NROY,seq(p+1,p+2*no.pcs,by=2)],
                        Variance = emulated.pca$X[NROY,seq(p+2,p+2*no.pcs,by=2)],
                        Error = diag(obs.err,k),Disc = diag(model.disc,k),
                        BasisUncertainty = TRUE)
  
  # update NROY space
  NROY = NROY[which(hm.pca$inNROY==TRUE)]
  
  # prop. of NROY space remaining
  NROY.prop = length(NROY)/emulated.pca$NROY.size
  
  samples = NULL
  if(length(NROY)<no.new.samples){
    set.seed(setseed)
    samples = emulated.pca$X[sample(NROY,length(NROY)),]
  } else if(NROY.prop!=0){
    # sample from nroy space using seed to ensure reproducibility
    set.seed(setseed)
    samples = emulated.pca$X[sample(NROY,no.new.samples),]
  }
  
  setwd(paste0(base_folder,'/models/',model.name))
  
  return(list(HM = hm.pca,emulation = emulated.pca,X = emulated.pca$X[,1:p],
              output = emulated.pca$X[,-c(1:p)],pca = rotated.pca,
              new.samples = samples,NROY = NROY,NROY.prop = NROY.prop))
}