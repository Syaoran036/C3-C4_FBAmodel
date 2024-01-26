# flux variability analysis of a minimization of total flux (MTF) result
#09/04/14

# to do: 
# - get correct namespace for .generateWT and so on

#  sysBiolAlg_mtffvClass.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2013 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of sybil.
#
#  Sybil is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Sybil is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybil.  If not, see <http://www.gnu.org/licenses/>.


#------------------------------------------------------------------------------#
#                   definition of the class sysBiolAlg_mtf                     #
#------------------------------------------------------------------------------#

setClass(Class = "sysBiolAlg_mtffv",#dh
         representation(
           maxobj = "numeric"
         ),
         contains = "sysBiolAlg"
)


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

# contructor for class sysBiolAlg_mtf
setMethod(f = "initialize",
          signature = "sysBiolAlg_mtffv",# dh
          definition = function(.Object,
                                model,
                                sumAbsFlux = NULL, # dh: A single numeric value giving the optimal value of the weighted sum of absolute fluxes
                                wtobj = NULL,
                                react = NULL, lb = NULL, ub = NULL,
                                costcoefSumAbsFlux = NULL,#dh: the weights that were used to calculate sumAbsFlux
                                #dh:# costcoefbw = NULL,
                                #dh:#costcoeffw = NULL,
                                useNames = SYBIL_SETTINGS("USE_NAMES"),
                                cnames = NULL,
                                rnames = NULL,
                                pname = NULL,
                                scaling = NULL,
                                writeProbToFileName = NULL, ...) {
            
            if ( ! missing(model) ) {
              #browser()
              if(is.null(sumAbsFlux)){  #dh
                    stop("sumAbsFlux has to be specified.") #dh
              } #dh
              
              if (is.null(wtobj)) {
                  tmp <- .generateWT(model, react, lb, ub, ...)
                  wtobj  <- tmp[["obj"]]
              }
              
              
              
              stopifnot(is(model, "modelorg"), is(sumAbsFlux, "numeric"), is(wtobj, "numeric")) #dh
              
              # If wtobj is longer than 1, mtf algorithm has to run several
              # times. In that case, wtobj is not written in the problem
              # object, it is written separately (maxobj) and used for
              # each iteration.
              
              if (length(sumAbsFlux) > 1) {
                stop("sumAbsFlux should be of length 1.")
                #maxobj <- wtobj
                #currmo <- 0
              }
              else {
                maxobj <- NULL
                currmo <- wtobj[1]
              }
              
              if (length(wtobj) > 1) { # dh
                stop("wtobj should be of length 1.") #dh
              } # dh
              
              #  the problem: maximize:
              #
              #            |      |      |
              #         S  |  0   |  0   |  = b
              #            |      |      |
              #       -------------------------
              #            |      |      |
              #         1  |  1   |  0   | >= 0
              #            |      |      |
              #       -------------------------
              #            |      |      |
              #         -1 |  0   |  1   | >= 0
              #            |      |      |
              #       -------------------------
              #       c_wt |  0   |  0   | >= c^T * v_wt
              #         0  |  costcoefsumAbsFlux   |  costcoefsumAbsFlux   | <= sumAbsFlux
              #            |      |      |
              #  lb   wt_lb|  0   |  0   |
              #  ub   wt_ub|10000 |10000 |
              #            |      |      |
              #  obj obj_coef(model)  |  0   |  0   |   # this will be set by optimizer() in fluxVar()
              
              
              # ---------------------------------------------
              # problem dimensions
              # ---------------------------------------------
              
              nc     <- react_num(model)
              nr     <- met_num(model)
              
              nCols  <- 3*nc
              #nRows  <- nr + 2*nc + 1   
              nRows  <- nr + 2*nc + 2   #dh, further row to constrain sum of absolute fluxes
              
              absMAX <- SYBIL_SETTINGS("MAXIMUM")
              
              
              # ---------------------------------------------
              # constraint matrix
              # ---------------------------------------------
              
              # the initial matrix dimensions
              LHS <- Matrix::Matrix(0, 
                                    nrow = nRows,
                                    ncol = nCols,
                                    sparse = TRUE)
              
              # rows for the mutant strain
              LHS[1:nr,1:nc] <- S(model)
              
              # location of the mutant strain
              fi <- c(1:nc)
              
              # rows for the delta match matrix
              diag(LHS[(nr+1)   :(nr+nc)  ,1       :nc    ]) <-  1
              diag(LHS[(nr+1)   :(nr+nc)  ,(nc+1)  :(2*nc)]) <-  1
              diag(LHS[(nr+nc+1):(nr+2*nc),1       :nc    ]) <- -1
              diag(LHS[(nr+nc+1):(nr+2*nc),(2*nc+1):(3*nc)]) <-  1
              
              # fix the value of the objective function
              LHS[(nr+2*nc+1),1:nc] <- obj_coef(model)    
              
              #dh:
              if (is.null(costcoefSumAbsFlux)) {
                costcoefSumAbsFlux<- rep(1, nc)
              }
              stopifnot(is(costcoefSumAbsFlux,"numeric"))
              
              # dh: fix the value of the weighted mtf objective (ie sum of absolute flux)
              LHS[(nr+ 2*nc +2),(nc+1):(3*nc)]<-rep(costcoefSumAbsFlux,2)
              
              # ---------------------------------------------
              # columns
              # ---------------------------------------------
              
              
              
              # boundaries:
              lower  <- c(lowbnd(model), rep(0, 2*nc))
              upper  <- c(uppbnd(model), rep(absMAX, 2*nc))
              
              
              # ---------------------------------------------
              # rows
              # ---------------------------------------------
              
              # dh: ?loadLPprob explains rtype
              
              #rlower <- c(rhs(model), rep(0, 2*nc), currmo)
              #rupper <- c(rhs(model), rep(absMAX, 2*nc + 1))
              
              
              # dh: # rlower <- c(rep(0, nr), rep(0, 2*nc), currmo)
              rlower <- c(rep(0, nr), rep(0, 2*nc), currmo, 0)
              #dh:# rupper <- c(rep(0, nr), rep(absMAX, 2*nc + 1))
              rupper <- c(rep(0, nr), rep(absMAX, 2*nc + 1), sumAbsFlux)
              # dh: # rtype  <- c(rep("E", nr), rep("L", 2*nc + 1))
              rtype  <- c(rep("E", nr), rep("L", 2*nc + 1),"U")
              
              # ---------------------------------------------
              # objective function
              # ---------------------------------------------
              
              # dh:
#               if (is.null(costcoeffw)) {
#                 fw <- rep(1, nc)
#               }
#               else {
#                 stopifnot(is(costcoeffw, "numeric"),
#                           (length(costcoeffw) == nc))
#                 fw <- costcoeffw
#               }
#               
#               if (is.null(costcoefbw)) {
#                 bw <- fw
#               }
#               else {
#                 stopifnot(is(costcoefbw, "numeric"),
#                           (length(costcoefbw) == nc))
#                 bw <- costcoefbw
#               }
              
             

              
              # dh: # cobj <- c(rep(0, nc), bw, fw)
              cobj<- c(obj_coef(model),rep(0,2*nc)) # use the objective function of the model
              
              
              # ---------------------------------------------
              # row and column names for the problem object
              # ---------------------------------------------
              
              # dh: leave untouched for now
              
              if (isTRUE(useNames)) {
                if (is.null(cnames)) {
                  cn <- c(react_id(model),
                          paste("bw", react_id(model), sep = "_"),
                          paste("fw", react_id(model), sep = "_")
                  )
                  colNames <- .makeLPcompatible(cn, prefix = "x")
                }
                else {
                  stopifnot(is(cnames, "character"),
                            length(cnames) == nCols)
                  colNames <- cnames
                }
                
                if (is.null(rnames)) {
                  rn <- c(met_id(model),
                          paste("bw", 1:nc, sep = "_"),
                          paste("fw", 1:nc, sep = "_"),
                          "obj_wt"
                  )
                  rowNames <- .makeLPcompatible(rn, prefix = "r")
                }
                else {
                  stopifnot(is(rnames, "character"),
                            length(rnames) == nRows)
                  rowNames <- rnames
                }
                
                if (is.null(pname)) {
                  probName <- .makeLPcompatible(
                    paste("MTF", mod_id(model), sep = "_"),
                    prefix = "")
                }
                else {
                  stopifnot(is(pname, "character"),
                            length(pname) == 1)
                  probName <- pname
                }
              }
              else {
                colNames <- NULL
                rowNames <- NULL
                probName <- NULL
              }
              
              
              # ---------------------------------------------
              # build problem object
              # ---------------------------------------------
              
              # dh: callNextMethod will call first Method of contained class, in this case sysBiolAlg (?)

              .Object <- callNextMethod(.Object, 
                                        sbalg      = "mtffv",#dh
                                        pType      = "lp",
                                        scaling    = scaling, 
                                        fi         = fi,  # position von S
                                        nCols      = nCols,
                                        nRows      = nRows,
                                        mat        = LHS,
                                        ub         = upper,
                                        lb         = lower,
                                        obj        = cobj,
                                        rlb        = rlower, # row lower bounds
                                        rub        = rupper, # row upper bounds
                                        rtype      = rtype,  # upper or lower or equal
                                        lpdir      = "min",
                                        ctype      = NULL,
                                        cnames     = colNames,
                                        rnames     = rowNames,
                                        pname      = probName,
                                        algPar     = list("wtobj" = wtobj,
                                                          "sumAbsFlux"=sumAbsFlux,
                                                          "costcoefSumAbsFlux"=costcoefSumAbsFlux
                                                          ),#dh 
                                                          #dh:#"costcoefbw" = bw,
                                                          #dh:#"costcoeffw" = fw),
                                        ...)
              
              .Object@maxobj <- as.numeric(maxobj)
              
              if (!is.null(writeProbToFileName)) {
                writeProb(problem(.Object),
                          fname = as.character(writeProbToFileName))
              }
              #
              #                  # ---------------------------------------------
              #                  # build problem object
              #                  # ---------------------------------------------
              #
              #                  lp <- optObj(solver = solver, method = method)
              #                  lp <- initProb(lp, nrows = nRows, ncols = nCols)
              #
              #                  # ---------------------------------------------
              #                  # set control parameters
              #                  # ---------------------------------------------
              #
              #                  if (!any(is.na(solverParm))) {
              #                      setSolverParm(lp, solverParm)
              #                  }
              #    
              #
              #                  loadLPprob(lp,
              #                             nCols = nCols,
              #                             nRows = nRows,
              #                             mat   = LHS,
              #                             ub    = upper,
              #                             lb    = lower,
              #                             obj   = cobj,
              #                             rlb   = rlower,
              #                             rub   = rupper,
              #                             rtype = rtype,
              #                             lpdir = "min"
              #                  )
              #                  
              #                  if (!is.null(scaling)) {
              #                      scaleProb(lp, scaling)
              #                  }
              #                  
              #                  .Object@problem   <- lp
              #                  .Object@algorithm <- "mtf"
              #                  .Object@nr        <- as.integer(nRows)
              #                  .Object@nc        <- as.integer(nCols)
              #                  .Object@fldind    <- as.integer(fi)
              #                  validObject(.Object)
              
            }
            return(.Object)
          }
)


#------------------------------------------------------------------------------#
#                                other methods                                 #
#------------------------------------------------------------------------------#

setMethod("changeMaxObj", signature(object = "sysBiolAlg_mtffv"), #dh
          function(object, j) {
            
            if (!is.null(object@maxobj)) {
              changeRowsBnds(problem(object), i = nr(object),
                             lb = object@maxobj[j], ub = SYBIL_SETTINGS("MAXIMUM"))
            }
            
            return(invisible(TRUE))
          }
)
