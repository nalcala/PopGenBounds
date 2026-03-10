#' Computes differentiation statistics
#'
#' @param freqs An allele frequency matrix, where rows are populations and columns alleles
#'
#' @return A numerical vector
#' @export
#'
#' @examples
#' freqs <- matrix(c(1,0.5,0,0.5),nrow=2)
#' Diff(freqs)
Diff = function(freqs){
  if( !all(round(rowSums(freqs*100))==100) ){
    stop("Rows of freqs matrix do not sum to 1. Do rows correspond to populations and columns to alleles?")
  }

  K = nrow(freqs)
  M = max(colMeans(freqs))
  HS = mean(1-rowSums(freqs**2))
  HT = 1-sum(colMeans(freqs)**2)
  FST = 1-HS/HT
  GpST= FST*(K-1+HS)/(K-1)/(1-HS)
  D = K/(K-1)*(HT-HS)/(1-HS)

  stats = tibble::tibble(statistic=c("FST","G'ST","D"), value = c(FST,GpST,D) )

  # compute bounds
  bounds = tibble::tibble(statistic=c("FST","G'ST","D") , bound=c(Fup(K=K,M=M), Gpup(K=K,M=M), Dup(K=K,M=M) ) )
  # compute closeness to upper bound
  closeness = tibble::tibble( statistic=c("FST","G'ST","D"), bound_closeness = stats$value/bounds$bound)

  return(dplyr::bind_rows(tibble::tibble(statistic="M",value=M),
                          dplyr::left_join(dplyr::left_join(stats,bounds,by="statistic"),closeness,by="statistic")))
}


#' Computes the maximal FST in terms of the most frequent allele
#'
#' @param K Number of subpopulations
#' @param M Frequency of the most frequent allele
#'
#' @return A numerical value
#' @export
#'
#' @examples
#' #' Fup(2,0.5)
Fup = function(K,M){
  s1 = K*M
  J = ceiling(1/s1)
  res = (K-1)*(1-s1*(J-1)*(2-J*s1))/(K-(1-s1*(J-1)*(2-J*s1)))*(s1<=1) +
    (K*(K-1)-s1**2+floor(s1)-2*(K-1)*(s1%%1)+(2*K-1)*(s1%%1)**2)/(K*(K-1) -s1**2-floor(s1)+2*s1-(s1%%1)**2)*(s1>1)
  return(res)
}


#' Computes the maximal D in terms of the most frequent allele
#'
#' @param K Number of subpopulations
#' @param M Frequency of the most frequent allele
#'
#' @return A numerical value
#' @export
#'
#' @examples
#' Dup(2,0.5)
Dup = function(K,M){
  s1 = K*M
  res = 1*(s1<=1) + (1-(s1**2-(s1%%1)**2-floor(s1))/(K-1)/(K-1+(s1%%1)**2+(1-s1%%1)**2) )*(s1>1)
  return(res)
}

#' Computes the maximal G'ST in terms of the most frequent allele
#'
#' @param K Number of subpopulations
#' @param M Frequency of the most frequent allele
#'
#' @return A numerical value
#' @export
#'
#' @examples
#' Gpup(2,0.5)
Gpup = function(K,M){
  s1 = K*M
  J = ceiling(1/s1)
  res = 1*(s1<=1) +
    (K*(K-1)-s1**2+floor(s1)-2*(K-1)*(s1%%1)+(2*K-1)*(s1%%1)**2)/(K*(K-1) -s1**2-floor(s1)+2*s1-(s1%%1)**2)*(s1>1)/(K-1)/(K-1+(s1%%1)^2+(1-s1%%1)^2 )*(K^2 - K+1-(s1%%1)^2-(1-s1%%1)^2 )
  return(res)
}


#' Plots differentiation statistics along with the maximal and minimal values of the statistic
#'
#' @param FST Numerical vector of FST values for each locus
#' @param GpST Numerical vector of G prime ST values for each locus
#' @param D Numerical vector of D values for each locus
#' @param M Numerical vector of frequency of the most frequent allele for each locus
#' @param K Number of subpopulations
#'
#' @return A ggplot object
#' @import tidyverse
#' @export
#'
#' @examples
#' library(tidyverse)
#' freqs_locus1 <- matrix(c(1,0.5,0,0.5),nrow=2)
#' freqs_locus2 <- matrix(c(1,0.8,0,0.2),nrow=2)
#' freqs_locus3 <- matrix(c(1,0.2,0,0.8),nrow=2)
#' data  <- rbind(Diff(freqs_locus1),Diff(freqs_locus2),Diff(freqs_locus2))
#' ggbounds(M=data %>% filter(statistic=="M") %>% pull(value),FST=data %>%
#' filter(statistic=="FST") %>% pull(value))
ggbounds = function(M,FST,GpST=NULL,D=NULL,K=2){
  nudge = (mean(M)<0.5)*0.16-(mean(M)>=0.5)*0.25
  MFtmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                        FST= Fup(K,seq(0.001,1-0.001,0.001)) )
  MGptmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                        GpST= Gpup(K,seq(0.001,1-0.001,0.001)) )
  MDtmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                        D= Dup(K,seq(0.001,1-0.001,0.001)) )
  plotFST <- ggplot2::ggplot(dplyr::tibble(M=M,FST=FST),ggplot2::aes(x=M,y=FST)) + ggpointdensity::geom_pointdensity() +
    ggplot2::geom_line(data=MFtmp,ggplot2::aes(x=M,y=FST)) + ggplot2::xlab(expression(italic(M))) +
    ggplot2::geom_point(data=dplyr::tibble(M=mean(M),FST=mean(FST)),col="red",
               pch=16,size=3,stroke=2) +
    ggplot2::geom_segment(data=dplyr::tibble(M=mean(M),FST=mean(FST)),
                          ggplot2::aes(x=M,xend=M,y=0,yend=Fup(K,M)), col="red",size=1 ) +
    ggplot2::geom_label(data=dplyr::tibble(M=mean(M),FST=mean(FST)),
                        ggplot2::aes(x=M,y=FST,
                   label= paste0("mean FST=", format(FST,digits=2)," (",
                                 format(mean(FST)/Fup(K,mean(M))*100,digits=2),"% of range)" )),
               nudge_x = nudge,nudge_y=0.07,col="red") +
    ggplot2::ylab(expression(italic(F[ST]))) +
    ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand = F) +
    ggplot2::theme_bw()

  if(!is.null(GpST)){
    plotGpST <- ggplot2::ggplot(dplyr::tibble(M=M,GpST=GpST),ggplot2::aes(x=M,y=GpST)) + ggpointdensity::geom_pointdensity() +
      ggplot2::geom_line(data=MGptmp,ggplot2::aes(x=M,y=GpST)) + ggplot2::xlab(expression(italic(M))) +
    ggplot2::geom_point(data=dplyr::tibble(M=mean(M),GpST=mean(GpST)),col="red",
                        pch=16,size=3,stroke=2) +
    ggplot2::geom_segment(data=dplyr::tibble(M=mean(M),GpST=mean(GpST)),
                          ggplot2::aes(x=M,xend=M,y=0,yend=Gpup(K,M)), col="red",size=1 ) +
    ggplot2::geom_label(data=dplyr::tibble(M=mean(M),GpST=mean(GpST)),
                        ggplot2::aes(x=M,y=GpST,
                                     label= paste0("mean GpST=", format(GpST,digits=2)," (",
                                                   format(mean(GpST)/Gpup(K,mean(M))*100,digits=2),"% of range)" )),
                        nudge_x = nudge,nudge_y=0.07,col="red") +
    ggplot2::ylab(expression(italic(G[ST]))) +
    ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand = F) +
    ggplot2::theme_bw()
  }else{plotGpST = NULL}

  if(!is.null(D)){
  plotD <- ggplot2::ggplot(dplyr::tibble(M=M,D=D),ggplot2::aes(x=M,y=D)) + ggpointdensity::geom_pointdensity() +
    ggplot2::geom_line(data=MDtmp,ggplot2::aes(x=M,y=D)) + ggplot2::xlab(expression(italic(M))) +
    ggplot2::geom_point(data=dplyr::tibble(M=mean(M),D=mean(D)),col="red",
                        pch=16,size=3,stroke=2) +
    ggplot2::geom_segment(data=dplyr::tibble(M=mean(M),D=mean(D)),
                          ggplot2::aes(x=M,xend=M,y=0,yend=Dup(K,M)), col="red",size=1 ) +
    ggplot2::geom_label(data=dplyr::tibble(M=mean(M),D=mean(D)),
                        ggplot2::aes(x=M,y=D,
                                     label= paste0("mean D=", format(D,digits=2)," (",
                                                   format(mean(D)/Dup(K,mean(M))*100,digits=2),"% of range)" )),
                        nudge_x = nudge,nudge_y=0.07,col="red") +
    ggplot2::ylab(expression(italic(D))) +
    ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand = F) +
    ggplot2::theme_bw()
  }else{plotD=NULL}

  if(length(M)>2) plot <- plot + ggplot2::scale_color_viridis_b()
  return(list(plotFST,plotGpST,plotD) )
}

#' Plots differentiation statistics along with the maximal and minimal values of the statistic
#'
#' @param freqs a matrix of allele frequency
#'
#' @return A ggplot object
#' @import tidyverse
#' @importFrom magrittr "%>%"
#' @export
#'
#' @examples
#' popa = matrix(c(0.3,0.3,0.3,0.1,rep(0,4*5),
#'   rep(0,4),0.3,0.3,0.3,0.1,rep(0,4*4),
#'   rep(0,4*2),0.3,0.3,0.3,0.1,rep(0,4*3),
#'   rep(0,4*3),0.3,0.3,0.3,0.1,rep(0,4*2),
#'   rep(0,4*4),0.3,0.3,0.3,0.1,rep(0,4),
#'   rep(0,4*5),0.3,0.3,0.3,0.1),ncol=24,byrow = TRUE)
#' ggfreqtable(popa)
#'
#' @importFrom rlang .data
ggfreqtable = function(freqs){
  freqs.tmp = freqs
  if(is.null(colnames(freqs)) ) colnames(freqs.tmp) = paste0(1:ncol(freqs.tmp))
  if(is.null(rownames(freqs)) ) rownames(freqs.tmp) = paste0(1:nrow(freqs.tmp))

  freqs.long = dplyr::bind_cols("Subpopulation"=rownames(freqs.tmp),dplyr::as_tibble(freqs.tmp)) %>%
    tidyr::pivot_longer(-"Subpopulation",names_to = "Allele",values_to = "Frequency")

  freqs.long$Allele = factor(freqs.long$Allele,levels=unique(freqs.long$Allele))
  freqs.long$Subpopulation = factor(freqs.long$Subpopulation,levels=rev(unique(freqs.long$Subpopulation)) )

  freqs.long$Frequency = as.numeric(format(freqs.long$Frequency,digits=2))

ggres = ggplot2::ggplot(freqs.long,ggplot2::aes(x=.data$Allele,y=.data$Subpopulation,fill=.data$Frequency)) + ggplot2::geom_tile(color="black",size=0.25) +
  ggplot2::geom_text(ggplot2::aes(label=.data$Frequency))+
  ggplot2::theme_bw() +
  ggplot2::scale_fill_gradientn(limits=c(0,1),colors=c("white","darkorange","darkmagenta"))+
  ggplot2::theme(
    #bold font for legend text
    legend.text=ggplot2::element_text(face="bold"),
    #set thickness of axis ticks
    axis.ticks=ggplot2::element_line(size=0.4),
    #remove plot background
    plot.background=ggplot2::element_blank(),
    #remove plot border
    panel.border=ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 90,vjust = 0.5,hjust=1)
  )

 return(ggres)
}
