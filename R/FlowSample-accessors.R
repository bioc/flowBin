#' @export
#' @docType methods
setGeneric("name", function(object) standardGeneric("name"))
#' @export
#' @docType methods
setGeneric("name<-", function(object,value) standardGeneric("name<-"))


setMethod("name", signature="FlowSample", function(object) {object@name})
setReplaceMethod("name", signature=c("FlowSample", "character"), function(object, value) {object@name <- value; object})

#' @export
setGeneric("tube.set", function(object) standardGeneric("tube.set"))
#' @export
setGeneric("tube.set<-", function(object,value) standardGeneric("tube.set<-"))
setMethod("tube.set", signature="FlowSample", function(object) {object@tube.set})
setReplaceMethod("tube.set", signature=c("FlowSample", "list"), function(object, value) {object@tube.set <- value; object})


#' @export
setGeneric("control.tubes", function(object) standardGeneric("control.tubes"))
#' @export
setGeneric("control.tubes<-", function(object,value) standardGeneric("control.tubes<-"))
setMethod("control.tubes", signature="FlowSample", function(object) {object@control.tubes})
setReplaceMethod("control.tubes", signature=c("FlowSample", "vector"), function(object, value) {object@control.tubes <- value; object})

#' @export
setGeneric("bin.pars", function(object) standardGeneric("bin.pars"))
#' @export
setGeneric("bin.pars<-", function(object,value) standardGeneric("bin.pars<-"))
setMethod("bin.pars", signature="FlowSample", function(object) {object@bin.pars})
setReplaceMethod("bin.pars", signature=c("FlowSample", "vector"), function(object, value) {object@bin.pars <- value; object})
#' @export
setGeneric("measure.pars", function(object) standardGeneric("measure.pars"))
#' @export
setGeneric("measure.pars<-", function(object,value) standardGeneric("measure.pars<-"))
setMethod("measure.pars", signature="FlowSample", function(object) {object@measure.pars})
setReplaceMethod("measure.pars", signature=c("FlowSample", "list"), function(object, value) {object@measure.pars <- value; object})
