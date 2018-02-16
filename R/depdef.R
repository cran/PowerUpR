# optimal functions are removed
optimal.cra4r4 <- function(...){
  .Defunct(new="cosa.crd4r4", package = "cosa")
}
optimal.bcra4r3 <- function(...){
  .Defunct(new="cosa.bcrd4r3", package = "cosa")
}
optimal.bcra4f3 <- function(...){
  .Defunct(new="cosa.crd3r3", package = "cosa")
}
optimal.bcra4r2 <- function(...){
  .Defunct(new="cosa.bcrd4r2", package = "cosa")
}
optimal.bira4r1 <- function(...){
  .Defunct(new="cosa.bird4r1", package = "cosa")
}
optimal.cra3r3 <- function(...){
  .Defunct(new="cosa.crd3r3", package = "cosa")
}
optimal.bcra3r2 <- function(...){
  .Defunct(new="cosa.bcrd3r2", package = "cosa")
}
optimal.bcra3f2 <- function(...){
  .Defunct(new="cosa.crd2r2", package = "cosa")
}
optimal.bira3r1 <- function(...){
  .Defunct(new="cosa.bird3r1", package = "cosa")
}
optimal.cra2r2 <- function(...){
  .Defunct(new="cosa.crd2r2", package = "cosa")
}
optimal.bira2r1 <- function(...){
  .Defunct(new="cosa.bird2r1", package = "cosa")
}
optimal.bira2f1 <- function(...){
  .Defunct(new="cosa.ird1r1", package = "cosa")
}
optimal.bira2c1 <- function(...){
  .Defunct(new="cosa.ird1r1", package = "cosa")
}
optimal.ira1r1 <- function(...){
  .Defunct(new="cosa.ird1r1", package = "cosa")
}
optimal.to.mdes <- function(...){
  .Defunct(msg="This function is no longer in use")
}
optimal.to.power <- function(...){
  .Defunct(msg="This function is no longer in use")
}


# deprecated and defunct arguments
.depdef <- function(x, envir = parent.frame()) {
  names.x <- names(x)
  if("R12" %in% names.x) {
    envir$r21 <- x$R12
    warning("'R12' is renamed and will be removed in the next version, use 'r21' instead", call. = FALSE)
  }
  if("R22" %in% names.x) {
    envir$r22 <- x$R22
    warning("'R22' is renamed and will be removed in the next version, use 'r22' instead", call. = FALSE)
  }
  if("R32" %in% names.x) {
    envir$r23 <- x$R32
    warning("'R32' is renamed and will be removed in the next version, use 'r23' instead", call. = FALSE)
  }
  if("R42" %in% names.x) {
    envir$r24 <- x$R42
    warning("'R42' is renamed and will be removed in the next version, use 'r24' instead", call. = FALSE)
  }

  if("RT22" %in% names.x) {
    envir$r2t2 <- x$RT22
    warning("'RT22' is renamed and will be removed in the next version, use 'r2t2' instead", call. = FALSE)
  }
  if("RT32" %in% names.x) {
    envir$r2t3 <- x$RT32
    warning("'RT32' is renamed and will be removed in the next version, use 'r2t3' instead", call. = FALSE)
  }
  if("RT42" %in% names.x) {
    envir$r2t4 <- x$RT42
    warning("'RT42' is renamed and will be removed in the next version, use 'r2t4' instead", call. = FALSE)
  }
  if("P" %in% names.x) {
    envir$p <- x$P
    warning("'P' is renamed and will be removed in the next version, use 'p' instead", call. = FALSE)
  }
  if("mdes" %in% names.x) {
    envir$es <- x$mdes
    warning("'mdes' is renamed and will be removed in the next version, use 'es' instead", call. = FALSE)
  }
}
