# modified based on function 'genExprs' in R Bioconductor package 'iCheck'

genEset = function (ex, pDat, fDat = NULL, annotation = "") 
{
  cn.dat <- colnames(ex)
  rn.pDat <- rownames(pDat)
  aa <- match(cn.dat, rn.pDat)
  if (length(cn.dat) != length(rn.pDat)) {
    stop("No. of columns of ex is not equalt to that of pDat!\n")
  }
  if(!identical(cn.dat, rn.pDat))
  {
    stop("Column names of ex != row names of pDat!\n")
  }
  pDat2 <- as(pDat, "data.frame")
  aa <- new("AnnotatedDataFrame", data = pDat2)
  exprs <- as(ex, "matrix")
  es <- new("ExpressionSet", exprs = exprs, phenoData = aa, 
    annotation = annotation)

  if(!is.null(fDat))
  {
    if(identical(rownames(ex), rownames(fDat)))
    {
      Biobase::fData(es) = fDat
    } else {
      stop("Row names of ex != row names of fDat!\n")
    }
  }

  invisible(es)
}
