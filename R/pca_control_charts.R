
#' Multivariate Control Charts based on Principal Components
#'
#' Draws a confidence ellipse based on the first two principal components and a multivariate quality control chart using the others principals components
#'
#' @param x dataset.
#' @param alpha confidence level of the quality control charts.
#' @param scaled a logical value indicating whether the variables should be standardized.
#'
#' @author Gualberto Agamez Montalvo; Bruno Marinho Estevam de Oliveira

#'
#'
#' @export
#'
#' @examples
#' data(Overtime)
#' pca_control_charts(Overtime,alpha = .05,scaled = TRUE)
#'
pca_control_charts <- function(x,alpha = 0.05,scaled = FALSE){
  if(REdaS::bart_spher(x)$p.value<0.05){
    cp = stats::prcomp(x,scale. = scaled);
    graphics::par(mfrow = c(1, 2))
    e = MSQC::ellip(type = "chi", cp$x[,1:2], alpha = alpha)
    ylimit = max(cp$x[,2],e[,2])
    xlimit = max(cp$x[,1],e[,1])
    plot(e,
         type = 'l',
         main = expression(Confidence~Ellipse),
         xlab = "PC1",
         ylab = "PC2",
         lty = 2,col = 'red',
         xlim = c(-xlimit-1,xlimit+1),
         ylim = c(-ylimit-1,ylimit+1));
    graphics::points(cp$x[,1],cp$x[,2],pch=20)

    k = ncol(cp$x)
    autovalor = cp$sdev^2
    T2 = 0
    for (i in 3:k) {
      T2 = T2 + cp$x[,i]^2/autovalor[i]
    }

    LSC = stats::qchisq(1-alpha,k-2);LSC
    LIC = 0
    plot(T2,pch=20,type = 'o',
         main = expression(Hotelling~T^2~Control~Chart),
         xlab = "Observation",
         ylab = expression(T^2~statistic),
         ylim = c(0,max(T2+1,LSC)))
    graphics::abline(b=0,LSC,lty = 2,col = 'red')
    graphics::abline(b=0,LIC,lty = 2,col = 'red')


  }else{warning("The correlation matrix not diverges significantly from the identity matrix")}
}


