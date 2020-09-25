# scatter plot with fitted lines and p-values

# pval - p-value of testing for differential correlation
plotDiffCor = function(x1,
                       z1,
                       x0,
                       z0,
                       pval = NULL,
                       xlab="gene1",
                       ylab="gene2",
                       title = "scatter plots"
                       )
{
  if(is.null(pval))
  {
    title2 = title
  } else {
    title2=paste(title, "\n(pval=", sprintf("%3.2e", pval), ")", sep="")
  }


  dat=data.frame(x=c(x1, x0), z=c(z1, z0), grp=c(rep("case", length(x1)), rep("control", length(x0)) ) )

  grp = NULL
  x = NULL
  z = NULL

  g = ggplot(dat, aes(x=x, y=z, shape=grp, color=grp)) + geom_point(size=2L)
  g = g + geom_smooth(method=lm, se=FALSE, fullrange=TRUE)
  g = g + facet_grid(. ~ grp)
  g = g + xlab(xlab) + ylab(ylab) + ggtitle(title2)
  g = g + theme(
    axis.line = element_line(arrow = arrow()),
    panel.background=element_rect(fill = "white", colour = "white"),
    plot.background=element_rect(fill = "white", colour = "white"))
  g = g + theme(legend.position = "none")

  print(g)

  dat1=data.frame(x1=x1, z1=z1)
  res1=lm(z1~x1, data=dat1)
  coef1=res1$coefficients

  dat0=data.frame(x0=x0, z0=z0)
  res0=lm(z0~x0, data=dat0)
  coef0=res0$coefficients

  res = list(g = g, dat=dat, coef1=coef1, coef0=coef0)
  invisible(res)

}
