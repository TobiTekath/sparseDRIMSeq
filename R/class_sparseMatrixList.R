#' @include class_show_utils.R
#' @import Matrix
NULL

#' @name matrixORsparseMatrix-class
setClassUnion("matrixORsparseMatrix", c("matrix", "sparseMatrix"))

###############################################################################
### sparseMatrixList class
###############################################################################

#' sparseMatrixList object
#' 
#' A sparseMatrixList object is a container for a list of matrices which have the same
#' number of columns but can have varying number of rows. Additionally, one can
#' store an extra information corresponding to each of the matrices in
#' \code{metadata} matrix.
#' 
#' @return
#' 
#' \itemize{ 
#' \item \code{names(x)}, \code{names(x) <- value}: Get or set names
#' of matrices. 
#' \item \code{rownames(x)}, \code{rownames(x) <- value},
#' \code{colnames(x)}, \code{colnames(x) <- value}: Get or set row names or
#' column names of unlistData slot. 
#' \item \code{length(x)}: Get the number of
#' matrices in a list. 
#' \item \code{elementNROWS(x)}: Get the number of rows of each of
#' the matrices. 
#' \item \code{dim(x)}, \code{nrow(x)}, \code{ncol(x)}: Get the
#' dimensions, number of rows or number of columns of unlistData slot. 
#' \item
#' \code{x[[i]]}, \code{x[[i, j]]}: Get the matrix i, and optionally, get only
#' columns j of this matrix. 
#' \item \code{x$name}: Shortcut for
#' \code{x[["name"]]}. 
#' \item \code{x[i, j]}: Get a subset of sparseMatrixList that
#' consists of matrices i with columns j. }
#' 
#' @param x sparseMatrixList object.
#' @param value,i,j,name Parameters used for subsetting and assigning new
#'   attributes to x.
#'   
#' @slot unlistData Matrix which is a row binding of all the matrices in a list.
#' @slot partitioning List of indexes which defines the row partitioning of
#'   unlistData matrix into the original matrices.
#'   
#' @slot metadata Matrix of additional information where each row corresponds to
#'   one of the matrices in a list.
#' @exportClass sparseMatrixList
setClass(Class = "sparseMatrixList",
         slots = c(unlistData = "matrixORsparseMatrix",
                   partitioning = "list",
                   metadata = "matrix"), package = "DTUrtle")


###################################

setValidity("sparseMatrixList", function(object){
  # has to return TRUE when valid object!
  
  partitioning_unlist <- unlist(object@partitioning)
  
  if(length(partitioning_unlist) == nrow(object@unlistData))
    out <- TRUE
  else
    return(paste0("Unequal lengths of partitioning indexes and 
      rows in unlistData: ", length(partitioning_unlist), " and ", 
      nrow(object@unlistData)))
  
  if(nrow(object@metadata) > 0){
    
    if(nrow(object@metadata) == length(object@partitioning))
      out <- TRUE
    else
      return(paste0("Unequal lengths of partitioning and metadata: ", 
        length(object@partitioning), " and ", nrow(object@metadata)))
    
  }
  
  return(out)
  
})


###############################################################################
### sparseMatrixList
###############################################################################

sparseMatrixList <- function(..., metadata){
  
  listData <- list(...)
  
  if (length(listData) == 1L && is.list(listData[[1L]]))
    listData <- listData[[1L]]
  
  if (length(listData) == 0L) {
    return(new("sparseMatrixList"))
    
  } else {
    
    if (!all(sapply(listData, is, "matrix"))&!all(sapply(listData, is, "sparseMatrix")))
      stop("all elements in '...' must be matrices!")
    
    unlistData <- do.call(rbind, listData)
    
    w <- sapply(listData, nrow)
    
    partitioning <- vector("list", length(w))
    
    inds <- 1:nrow(unlistData)
    names(inds) <- rownames(unlistData)
    
    partitioning[w != 0] <- split(inds, rep(1:length(w), w))
    
    if(!is.null(names(listData)))
      names(partitioning) <- names(listData)
    
    if(!missing(metadata))
      return(new("sparseMatrixList", unlistData = unlistData, 
                 partitioning = partitioning, 
                 metadata = metadata))
    else
      return(new("sparseMatrixList", unlistData = unlistData, 
                 partitioning = partitioning))
    
  }
  
}



################################################################################
### show method
################################################################################


setMethod("show", "sparseMatrixList", function(object){
  
  nhead <- 2
  
  nl <- length(object)  
  cat(mode(object@unlistData),"MatrixList of length", nl,"\n")
  
  if(nl > 0){
    np <- min(nl, nhead)
    
    object_sub <- object[1:np]
    
    if(is.null(names(object_sub)))
      print_names <- paste0("[[", 1:np, "]]\n")
    else 
      print_names <- paste0("$", names(object_sub), "\n")
    
    # for(i in 1:np){
    #   # i = 1
    #   cat(print_names[i])
    #   show_matrix(object_sub[[i]])
    #   cat("\n")
    #   
    # }
    
    if(np < nl){
      
      if(is.null(names(object_sub)))
        cat(paste0("[[...]]\n"))
      else 
        cat(paste0("$...\n"))
      
    }
    
  }
  
  # if(nrow(object@metadata) != 0){
  #   cat("\nwith metadata slot\n")
  #   show_matrix(object@metadata)
  # }
  
  
})


###############################################################################
### accessing methods
###############################################################################

#' @rdname sparseMatrixList-class
#' @export
setMethod("names", "sparseMatrixList", function(x){
  
  names(x@partitioning)
  
})


#' @rdname sparseMatrixList-class
#' @export
setMethod("names<-", "sparseMatrixList", function(x, value){
  
  names(x@partitioning) <- value
  x
  
})


#' @rdname sparseMatrixList-class
#' @export
setMethod("rownames", "sparseMatrixList", function(x){
  
  rownames(x@unlistData)
  
})


#' @rdname sparseMatrixList-class
#' @export
setMethod("rownames<-", "sparseMatrixList", function(x, value){
  
  rownames(x@unlistData) <- value
  x
})

#' @rdname sparseMatrixList-class
#' @export
setMethod("colnames", "sparseMatrixList", function(x){
  
  colnames(x@unlistData)
  
})

#' @rdname sparseMatrixList-class
#' @export
setMethod("colnames<-", "sparseMatrixList", function(x, value){
  
  colnames(x@unlistData) <- value
  x
})


#' @rdname sparseMatrixList-class
#' @export
setMethod("length", "sparseMatrixList", function(x){
  
  length(x@partitioning)
  
})


#' @rdname sparseMatrixList-class
#' @export
#' @importFrom S4Vectors elementNROWS
setMethod("elementNROWS", "sparseMatrixList", function(x){
  
  sapply(x@partitioning, length)
  
})


#' @rdname sparseMatrixList-class
#' @export
setMethod("dim", "sparseMatrixList", function(x){
  
  dim(x@unlistData)
  
})


#' @rdname sparseMatrixList-class
#' @export
setMethod("nrow", "sparseMatrixList", function(x){
  
  nrow(x@unlistData)
  
})

#' @rdname sparseMatrixList-class
#' @export
setMethod("ncol", "sparseMatrixList", function(x){
  
  ncol(x@unlistData)
  
})


###############################################################################
### subsetting methods
###############################################################################


#' @aliases [[,sparseMatrixList-method
#' @rdname sparseMatrixList-class
#' @export
setMethod("[[", signature(x = "sparseMatrixList"), function(x, i, j){
  
  if(!missing(j))
    return(x@unlistData[x@partitioning[[i]], j , drop = FALSE])
  else
    return(x@unlistData[x@partitioning[[i]], , drop = FALSE])
  
})


#' @rdname sparseMatrixList-class
#' @export
setMethod("$", "sparseMatrixList", function(x, name){
  
  x[[name]]
  
})


################################

#' @aliases [,sparseMatrixList-method [,sparseMatrixList,ANY-method
#' @rdname sparseMatrixList-class
#' @export
setMethod("[", signature(x = "sparseMatrixList"), function(x, i, j){
  
  if(!missing(i)){
    
    if(!missing(j)){
      
      if(nrow(x@metadata) != 0)
        return(new("sparseMatrixList", 
          unlistData = x@unlistData[unlist(x@partitioning[i]), j, drop = FALSE], 
          partitioning = relist(1:nrow(x@unlistData), x@partitioning[i]), 
          metadata = x@metadata[i, , drop = FALSE]))
      else
        return(new("sparseMatrixList", 
          unlistData = x@unlistData[unlist(x@partitioning[i]), j, drop = FALSE], 
          partitioning = relist(1:nrow(x@unlistData), x@partitioning[i]), 
          metadata = x@metadata))
      
      
    }else{
      
      if(nrow(x@metadata) != 0)
        return(new("sparseMatrixList", 
          unlistData = x@unlistData[unlist(x@partitioning[i]), , drop = FALSE], 
          partitioning = relist(1:nrow(x@unlistData), x@partitioning[i]), 
          metadata = x@metadata[i, , drop = FALSE]))
      else
        return(new("sparseMatrixList", 
          unlistData = x@unlistData[unlist(x@partitioning[i]), , drop = FALSE], 
          partitioning = relist(1:nrow(x@unlistData), x@partitioning[i]), 
          metadata = x@metadata))   
      
    } 
    
  }else{
    
    if(!missing(j)){
      
      return(new("sparseMatrixList", 
        unlistData = x@unlistData[, j, drop = FALSE], 
        partitioning = x@partitioning, 
        metadata = x@metadata))
      
    }else{
      
      return(x)
      
    } 
    
  }
  
})


################################

#' @aliases [,matrixORsparseMatrix-method [,matrixORsparseMatrix,ANY-method
#' @rdname matrixORsparseMatrix-class
#' @export
setMethod("[", signature(x = "matrixORsparseMatrix"), function(x, i, j, drop = F){
  if(is(x, "matrix")){
    x <- as(x, "matrix")
  } else{
    x <- as(x, "sparseMatrix")
  }
  
  if(!missing(i)){
    
    if(!missing(j)){
      return(x[i, j, drop = drop])
    }else{
      return(x[i, , drop = drop])
    }
  }else{
    if(!missing(j)){
      
      return(x[, j, drop = drop])
      
    }else{
      return(x)
    }
  }
})


#' @aliases [[,matrixORsparseMatrix-method
#' @rdname matrixORsparseMatrix-class
#' @export
setMethod("[[", signature(x = "matrixORsparseMatrix"), function(x, i, j, drop = F){
  if(!missing(i)){
    if(!missing(j)){
      return(x[i, j, drop = drop])
    }else{
      return(x[i, , drop = drop])
    }
  }else{
    if(!missing(j)){
      
      return(x[, j, drop = drop])
      
    }else{
      return(x)
    }
  }
})







