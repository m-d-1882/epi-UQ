# © Copyright University of Exeter. This code is open source and released under the MIT Licence

# This script is the hub for Uncertainty Quantification (UQ) for this project.
# We cover emulation and history matching for large outputs.
# 
# Prior to using this hub, the user must have supplied a script in 
# '~/models/"model-name".R' that gives the following when sourced:
#   X: n x p size matrix. These are design point inputs.
#   Y: n x q size matrix. These are design point outputs
#   X.valid: n' x p size matrix. These are validation point inputs.
#   Y.valid: n' x q size matrix. These are validation point outputs.
#   var.names: p length vector. Names of all input variables.
#   lower.bounds: p length vector. Lower bound for each input variable.
#   upper.bounds: p length vector. Upper bound for each input variable.
#   z: q length vector. Observation for purposes of history matching.
#   obs_err: q length vector. Observation error for purposes of history matching and constraining basis vectors.
#   model_disc: q length vector. Model discrepancy for purposes of history matching and constraining basis vectors.
# Alternatively one can load these in via '~/models/"model-name".Rdata'

# This script covers emulation and history matching for large numbers of 
# outputs. We use principal component analysis (PCA) to reduce the number of 
# outputs for emulation purposes and use rotation of the bases' to ensure the 
# reproducibility of the observation for history matching. We also conduct 
# history matching when data is supplied for each wave.

# install packages needed and load them in
install.packages(c('tidyverse','GenSA','lhs','tensor','far','parallel',
                   'DiceKriging','MASS'),quiet = TRUE)

library(tidyverse)
library(GenSA)
library(lhs)
library(tensor)
library(far)
library(parallel)
library(DiceKriging)
library(MASS)

#define folder path where repository was downloaded
base_folder = ''

setwd(base_folder)

# load in R functions
source(paste0(base_folder,'/functions/UQ_functions.R'))
source(paste0(base_folder,'/functions/rotation_functions.R'))
source(paste0(base_folder,'/functions/FastHM.R'))

# model name
model.name = 'WSS'

setwd(paste0(base_folder,'/models/',model.name))

# load in model data to perform UQ
load(paste0(base_folder,'/models/',model.name,'/',model.name,'.Rdata'))

### EMULATION

# conduct PCA on outputs to reduce dimensionality
pca1 = pca(Y = Y,W = diag(obs_err + model_disc))

# plot the variance explained by each principal component
plot(x = 1:length(pca1$VarExpl),y = pca1$VarExpl,type = 'p',
     xlab = 'lambda',ylab = 'Prop. variance expl.')

# based off graph, choose number of principal components
no.pcs = 7

# use emulation_RB to emulate weights
pca_emulation = emulation_weights(X = X,
                             Y = matrix(pca1$weights[,1:no.pcs],ncol=no.pcs),
                             setseed = 1)

# predict weights/full outputs given set of inputs
pca_prediction = predict_outputs(pred.inputs = matrix(X.valid[1:30,],ncol = ncol(X.valid)),
                                 emulators = pca_emulation,full_output = TRUE,
                                 basis.vectors = pca1$basis.vectors[,1:no.pcs])

# draw curves from emulator given input
pred.curves = prediction_curves(pred.input = matrix(X.valid[1,],nrow = 1),
                                emulators = pca_emulation,
                                basis.vectors = pca1$basis.vectors[,1:no.pcs],
                                no.draws = 200)

### HISTORY MATCHING
# check size of reconstruction error of observation. wmse
ReconError(obs = z,basis = pca1$basis.vectors[,1:5],
           weightinv = GetInverse(diag(obs_err+model_disc)))

# rotate the basis towards the observation if needed
rotated.pca = salter_rotation_basis(basis.vectors = pca1$basis.vectors,Y = Y,
                                    z = z,W=diag(obs_err + model_disc),
                                    no.change = FALSE)

rotated.pca$VarExp
# if necessary reduce/increase number of principal components (no.pcs) if 
# required. looking for lowest number of pcs that explain desired amount of 
# variance.
rotated.pca$recon_err
# also look at the reconstruction error for each principal component added
# (MSE).

# set no.pcs in 'hmwave' to be the maximum of the choices for rotated.pca$VarExp
# and rotated.pca$recon_err

# History match: am using initial model runs as the inputs/outputs for this wave
wave1 = hmwave(wave.no = 1,X = X,Y = Y,
               pca = rotated.pca,no.pcs = 1,z = z,
               obs.err = obs_err,model.disc = model_disc,var.names = var.names,
               upper.range = upper.range,lower.range = lower.range,setseed = 2)
save('wave1',file = 'historymatch/wave1.Rdata')
#load(file = 'historymatch/wave1.Rdata')

# load in data for next wave. Need:
# X.wave2: n' x p
# Y.wave2: n' x q

wave2 = hmwave(wave.no = 2,X = X.wave2,Y = Y.wave2,z = z,obs.err = obs_err,
               model.disc = model_disc,var.names = var.names,
               upper.range = upper.range,lower.range = lower.range,
               prev.wave = wave1,pca = wave1$pca,no.pcs = 1)
save('wave2',file = 'historymatch/wave2.Rdata')
load(file = 'historymatch/wave2.Rdata')
