

#' @title Oaxaca Decomposition for linear models
#' 
#' @author Christopher Lee
#' 
#' @description The function creates a seperate linear regression for each level in the \code{group}
#'  and apply the Blinder-Oaxaca Decomposition technique.
#'  
#'  @inheritParams lm
#'  
#'  @param group The string character for the column in the \code{data} which identifies different groups.
#'  The group variable must have at least 2 unique values.
#'  @param decomposition Which decomposition to use. Possible values are \code{oaxaca_0}, \code{oaxaca_1},
#'  \code{Cotton}, \code{Reimer} and \code{Neumark}. For more information on the meaning of these decompositions
#'  see '\code{Details}'.
#'  @param detailed If \code{TRUE}, prints the summaries of all models on top of decomposition summaries
#'  
#'  @details The gap between the predicted means of different groups will always be calculated so that 
#'  the gap is positive. For the moment, the function cannot show the interactions between the
#'  endowment gap (explained) and the coefficient gap (unexplained). The function can only handle up to
#'  decomposition for 3 groups.
#'  @export
#'  
oaxaca <- function(formula, data, group, detailed = FALSE, decomposition = 'oaxaca_0'){
  #Code to allow for the use of formulae in function
   call <- match.call(expand.dots = FALSE)
   f <- match(c("formula", "data"), names(call), 0L)
   call <- call[c(1L, f)]
   call[[1L]] <- quote(stats::model.frame)
   call <- eval(call, parent.frame())
   eqn <- formula(formula)
   call$group <- as.factor(data[,group])
   #Preparing empty lists for multiple for-loops
   out <- list()
   coef <- list()
   mean <- list()
   if(length(levels(call$group)) > 1) {
     for(i in 1:length(levels(call$group)))
       out$Models[[paste(group, "=",
                         levels(call$group)[i],sep="")]] <- lm(eqn, data=subset(call, group == levels(call$group)[i]))       
     out$Models[['Pooled']] <- lm(eqn, data)
    gap <- data.frame(matrix(ncol=1,nrow=length(levels(call$group))))
    for(i in 1 : length(levels(call$group)))
       gap[i,] <- mean(out$Models[[i]]$fitted.values)
    ord <- order(gap, decreasing = TRUE)
    A <- ord[1]
    B <- ord[2]
    C <- ord[3]
    gap <- mean(out$Models[[A]]$fitted.values) - mean(out$Models[[B]]$fitted.values)
   
    for(i in 1 : length(levels(call$group))){
      for(j in 1 : out$Model[[1]]$rank){
        coef[[paste('immigrant',levels(call$group)[i],sep="")]][[names(coef(out$Model[[1]]))[j]]] <- 
          coef(out$Model[[i]])[j]
      }
      for(j in 2 : out$Model[[1]]$rank){
        mean[[paste('immigrant',levels(call$group)[i],sep="")]][[names(coef(out$Model[[1]]))[j]]] <- 
          mean(call[j][call$group == levels(call$group)[i],])
      }
     }
    intercept <- function(x){
      x <- c(1,x)
      attr(x, "names")[1] <- "(Intercept)"
      return(x)
    }
    mean <- llply(mean, intercept)
    tab <- data.frame(Variables = names(coef(out$Model[[1]])))
    D <- 0
    if(decomposition == "oaxaca_0")
      D <- 0
    if(decomposition == "oaxaca_1")
      D <- 1
    if(decomposition == "Reimers")
      D <- 0.5
    if(decomposition == "Cotton")
      D <- length(out$Models[[A]]$resid) / length(out$Models[[B]]$resid)
   for(i in 1 : out$Model[[1]]$rank){
     tab$Endowment.Gap[i] <- (D * coef[[A]][i] + (1-D)* coef[[B]][i]) * (mean[[A]][i] - mean[[B]][i])
     tab$Coefficient.Gap[i] <- ((1-D) * mean[[A]][i] + D * mean[[B]][i]) * (coef[[A]][i] - coef[[B]][i])
   }
   if(decomposition == "Neumark")
     for (i in 1 : out$Model[[1]]$rank){
       tab$Endowment.Gap[i] <- coef(out$Models[['Pooled']])[i] * (mean[[A]][i] - mean[[B]][i])
       tab$Coefficient.Gap[i] <- mean[[A]][i] * (coef[[A]][i] - coef(out$Models[['Pooled']])[i]) + 
         mean[[B]][i] * (coef(out$Models[['Pooled']])[i] - coef[[B]][i])
     }
   tab2 <- data.frame(c(mean(out$Model[[A]]$fitted.value), 
                        mean(out$Model[[B]]$fitted.value), gap, 
                        colSums(tab[,-1])[1], colSums(tab[,-1])[2],
                        colSums(tab[,-1])[1] / gap, colSums(tab[,-1])[2] / gap))
   rownames(tab2) <- c(paste("Mean prediction for", ls(out$Model[A])),
                       paste("Mean prediction for", ls(out$Model[B])),
                       "Total Gap","- Endowment Gap","- Coefficient Gap",
                       "Explained Gap %", "Unexplained Gap %")
   colnames(tab2) <- attr(call,'names')[1]
   out$summary1 <- tab2
   out$summary2 <- tab

   if(detailed == TRUE){
     out$Models <- llply(out$Models,summary)
     return(out)
   }
   return(out)
   } else {
     stop("There are fewer than 2 groups")
   }
}
