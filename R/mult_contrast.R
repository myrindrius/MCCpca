#' Multivariate Mean Comparison Test
#'
#' This function applies a multivariate hypothesis test to compare means between groups
#'  according to a given vector or matrix of contrast.
#'
#' @param data a data matrix
#' @param groups a vector or factor object giving the group for the
#' corresponding elements of data matrix.
#' @param contraste a contrast vector or contrast matrix where the sum of each
#'  row is zero.
#'
#' @returns
#' A list of class "mult_contrast" containing
#'
#' @return xbar.ip: a matrix of mean vectors of wich group; delta: a linear
#' combination of the means; var.delta: a variance's matrix of delta; T2: the
#' Hotelling t^2 value; wilks: the Wilk's Lambda value.
#'
#'
#' @author Gualberto Agamez Montalvo; Bruno Marinho Estevam de Oliveira
#'
#'
#' @export
mult_contrast <- function(data,groups,contraste){
  if(ncol(data)>=2){
    dados = cbind(groups,data)
    p=ncol(dados)
    N = nrow(dados);N
    n=nrow(dados)/length(unique(groups));n
    k=length(unique(groups));k
    p=ncol(dados)-1;p

    ybar.ip=apply(dados[groups==levels(groups)[1],2:(p+1)], 2, mean)
    for (i in 2:k) {
      ybar=apply(dados[groups==levels(groups)[i],2:(p+1)], 2, mean)
      ybar.ip=cbind(ybar.ip,ybar)
      ybar=0
    }
    colnames(ybar.ip)=levels(groups)

    delta = (contraste%*%t(ybar.ip));delta

    model = stats::manova(data~groups)
    saida = summary(model)
    E = saida$SS$Residuals;E
    ve = model$df.residual;ve

    var.delta = (E/ve*n)*sum(contraste^2);var.delta

    # T2 de Hotelling
    T2 = delta%*%var.delta%*%t(delta);T2

    #Lambda de Wilks
    H1 = (n/sum(contraste^2))*t(delta)%*%t(t(delta))
    wilks = det(E)/det(H1+E);wilks

    #Aproximação pela F
    approx.F =((1 - wilks)/wilks)*(ve-p+1)/p;approx.F

    #P-valor
    p.value = 1 - stats::pf(approx.F,p,ve-p+1)

    #Resposta
    cat("Multivariate Mean Comparison Test\n\n Hotelling T2 =", T2 ,
        "\n Wilk's Lambda =", wilks,
        "\n Approximation by F =",approx.F,
        "\n p-value =",p.value)

    res <- structure(list(ybar.ip=ybar.ip,
                          delta=delta,
                          var.delta=var.delta,
                          T2=as.numeric(T2),
                          wilks=as.numeric(wilks),
                          approx.F=as.numeric(approx.F),
                          p.value = as.numeric(p.value)), class = "mctest")

  }else{warning("The data matrix must have one column for each variable")}
}
