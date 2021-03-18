

#
# Loads an eventually installs the required libraries
#

if (!require(ROBITools)) {

  # ROBITools are not available on CRAN and have to be installed
  # from http://git.metabarcoding.org using devtools

  if (!require(devtools)) {
    install.packages("devtools")
  }

  devtools::install_git("https://git.metabarcoding.org/obitools/ROBIUtils.git")
  devtools::install_git("https://git.metabarcoding.org/obitools/ROBITaxonomy.git")
  devtools::install_git("https://git.metabarcoding.org/obitools/ROBITools.git")

  library(ROBITools)
}

library(ROBITaxonomy)

if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}


plot_raw_length_dist <- function(stats,
                                 primer_label,
                                 min_cut,
                                 max_cut) {
  ys <- max(stats$total / stats$count)
  stats %>%
    ggplot(aes(y = total / count, x = seq_length)) +
    geom_point(
      col = "blue",
      cex = 0.5
    ) +
    scale_y_log10() +
    labs(title = primer_label) +
    annotate("rect",
      xmin = min_cut,
      xmax = max_cut,
      ymin = 1,
      ymax = ys,
      alpha = .5
    ) +
    geom_vline(xintercept = min_cut, col = "red") +
    geom_vline(xintercept = max_cut, col = "red")
}

#'
best_rank <- function(taxonomy,taxids) {
    sapply(mapply(taxonomicrank,
                  c(taxonomy),
                  path(taxonomy,taxids),
                  SIMPLIFY = FALSE),
                  function(x) rev(x[x != "no rank"])[1])

}

star_cluster <- function(sumatra,motu,threshold) {
  unique(c(motu,as.character(unlist(sumatra[(sumatra[,1] == motu |
                                             sumatra[,2] == motu) &
                                             sumatra[,3] <= threshold,1:2]))))
}


lcafast <- function(taxonomy,taxid,name=FALSE,check.taxid=TRUE) {


            #
            # Remove not valid taxid
            #

            if (check.taxid){
            	taxid = validate(taxonomy,taxid)
                if (any(is.na(taxid)))
                   return(NA)
            }

            ntaxid=length(taxid)
            ttaxid <- table(taxid)
            utaxid <- as.integer(names(ttaxid))
            allpath  = path(taxonomy,utaxid)
            maxlength= max(vapply(allpath,length,0))

            mp <- matrix(NA,nrow = length(utaxid),ncol=maxlength)

            for (i in seq_along(allpath)) {
              p <- allpath[[i]]
              mp[i,1:length(p)] <- p
            }

            lca    <- rep(NA_integer_,maxlength)
            lerror <- rep(NA_integer_,maxlength)
            selection <- rep(TRUE,length(utaxid))

            for (i in seq_len(maxlength)) {
               tt <- tapply(ttaxid[selection],mp[selection,i], sum, na.rm=TRUE)
               if (length(tt) == 0) break
               maxtaxid_pos <- which.max(tt)[1]
               lca[i] <- as.integer(names(maxtaxid_pos))
               selection <- mp[,i]==lca[i]
               lerror[i] <- ntaxid - tt[maxtaxid_pos]
            }

            ok_lca_na <- !is.na(lca)
            lerror <- lerror[ok_lca_na]
            names(lerror) <- if (name) scientificname(taxonomy,lca[ok_lca_na])
                             else lca[ok_lca_na]

            lerror[c((lerror > 0)[-1],TRUE)]
}

#' Title
#'
#' @param taxonomy a ROBITaxonomy object
#' @param motu     the motu to analyse
#' @param refdb    the complete reference DB provided as a metabarcoding object
#' @param sumatra  the sumatra result file
#' @param cluster_threshold Similarity distance used to cluster around the target motu
#' @param keep_threshold minimum limit bad taxa/total taxa
#' @param excluded_motus
#'
#' @return
#' @export
#'
#' @examples
progressive_lca <- function(taxonomy,motu,refdb,sumatra,
                            cluster_threshold = 0.1,
                            keep_threshold = 0.05,
                            excluded_motus = character(0),
                            limited = FALSE) {

  taxids  <- as.integer(rownames(refdb@reads))
  cluster <- intersect(unique(c(motu,unlist(sumatra[(sumatra[,1] == motu |
                                                     sumatra[,2] == motu) &
                                                     sumatra[,3] <= cluster_threshold,1:2]))),
                       colnames(refdb@reads))
  #cluster <- intersect(star_cluster(sumatra,motu,cluster_threshold),colnames(refdb@reads))
  weight <- rowSums(refdb@reads[,setdiff(cluster,excluded_motus),drop = FALSE])

  taxids <- rep(as.integer(taxids),weight)
  nt <- length(taxids)
  err_max = ceiling(nt * keep_threshold)

  uerror <- lcafast(taxonomy,taxids)

  ulevels <- names(uerror)
  ulcas   <- as.numeric(ulevels)


  ubad    <- lapply(ulcas, function(x) {t <- as.numeric(unique(taxids[which(!(is.subcladeof(taxonomy,taxids,x) | taxids == x))]))
                                        names(t) <- scientificname(taxonomy,t)
                                        t
                                       })
  ubadlca <- c(NA,sapply(seq_along(ubad[-1]) + 1,
                         function(x) lowest.common.ancestor(taxonomy,setdiff(ubad[[x]], ubad[[x - 1]]))))
  ubad <- ubad
  names(ubad) <- ulcas

  uscname <- scientificname(taxonomy,ulcas)
  urank <- taxonomicrank(taxonomy,ulcas)
  ubrank <- best_rank(taxonomy,ulcas)
  ulp    <- sapply(path(taxonomy,ulcas),length)


  rep <- data.frame(error = uerror,
                    good  = nt - uerror,
                    ratio = uerror/nt,
                    down    = ulp - min(ulp),
                    lca   = ulcas,
                    lca_name = uscname,
                    rank = urank,
                    best_rank = ubrank,
                    bad_lca = ubadlca,
                    bad_name = scientificname(taxonomy,ubadlca),
                    keep_ratio = c((uerror/nt)[-1] > keep_threshold,TRUE),
                    keep_lca   = c(FALSE,(c(NA,ulcas[-length(uerror)]) != ubadlca)[-1])
                    )

  rep$keep <- apply(rep[, grep("^keep_", colnames(rep))],
                    MARGIN = 1,
                    any)

  attr(rep,"bad") <- ubad
  attr(rep,"to_remove") <- intersect(as.character(ubad[[as.character(rep[rep$keep,][1,"lca"])]]),names(taxids(refdb,motu)))

  rep

}



best_match <- function(sumatra,motu) {
  matches <- sumatra[(sumatra[,1] == motu | sumatra[,2] == motu),]
  best <- matches[which.min(matches[,3]),]
  d <- best[,3]
  names(d) <- as.character(apply(best,MARGIN = 1, function(x) if (x[1] == motu) x[2] else x[1]))
  d
}


taxids <- function(refdb,motu) {
  t <- refdb@reads[,motu]
  t[t > 0]
}
##################
#
# Hill's number related functions
#
#################

log_q = function(x,q=1) {
  if (q == 1)
    log(x)
  else (x^(1 - q) - 1)/(1 - q)
}

exp_q = function(x, q = 1, base = exp(1)) {
  if (q == 1)
    exp(x)
  else
    (1 + (1 - q) * x)^(1 / (1 - q))
}

H_q = function(x,q=1,clades = NULL) {
  x = x / sum(x)
  x = x * log_q(1/x,q)
  if (is.null(clades)) sum(x,na.rm = TRUE)
  else tapply(x, clades, sum,na.rm = TRUE)
}

D_q = function(x,q=1,clades = NULL) {
  exp_q(H_q(x,q,clades = clades),q)
}

H_spectrum = function(x,q=1) {
  sapply(q,function(Q) H_q(x,Q))
}

D_spectrum = function(x,q=1) {
  sapply(q,function(Q) D_q(x,Q))
}

min_max_filter <- function(metabar,
                           threshold = 0,
                           rm_empty_sample = TRUE) {
  freq <- decostand(metabar@reads,method = "total")
  toplot <- metabar[,apply(freq,2,max) >= threshold]
  if (rm_empty_sample)
    toplot <- toplot[rowSums(toplot@reads) > 0,]
  toplot
}


plot_reads_x_motus <- function(metabar,threshold=0,q=0,limit=NULL) {
  toplot = min_max_filter(metabar,threshold)
  read_per_pcr <- rowSums(toplot@reads)
  hill_per_pcr <- apply(decostand(toplot@reads,
                                  method = "total"),
                        MARGIN = 1,D_q, q=q)

  data = data.frame(reads = read_per_pcr,
                    hills = hill_per_pcr,
                    category = toplot@samples$category)

  data %>% ggplot(aes(y = hills,
                      x = reads,
                      col = category)) +
         geom_point() +
         scale_x_log10() +
         scale_y_log10() -> p

  if (! is.null(limit))
    p <- p + geom_vline(xintercept = limit,
                        col = 2,
                        lty = 2)

  p
}



mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

tag_bad_pcr = function(samples,counts,plot = TRUE,threshold = NA) {
  counts = decostand(counts,method = "hellinger")

  bc = aggregate(counts,
                 by=list(factor(as.character(samples))),
                 mean)
  bc.name = as.character(bc[,1])
  bc = bc[-1]
  rownames(bc)=bc.name
  bc = bc[as.character(samples),]

  d = sqrt(rowSums((counts - bc)^2))
  names(d) = as.character(samples)

  d.m  = mode(d)
  d.sd = sqrt(sum((d[d <= d.m] - d.m)^2)/sum(d <= d.m))

  d.max = aggregate(d,
                    by = list(factor(as.character(samples))),
                    max)

  d.max.names = d.max[,1]
  d.max = d.max[,2]
  names(d.max) = d.max.names
  d.max = d.max[as.character(samples)]

  d.len = aggregate(d,
                    by = list(factor(as.character(samples))),
                    length)

  d.len.names = d.len[,1]
  d.len = d.len[,2]
  names(d.len) = d.len.names
  d.len = d.len[as.character(samples)]

  if (is.na(threshold))
    threshold = d.m + (d.sd*2)

  keep = ((d < threshold) | d!=d.max) & d.len > 1

  selection = data.frame(samples = as.character(samples),
                         distance= d,
                         maximum = d.max,
                         repeats = d.len,
                         keep    = keep,
                         stringsAsFactors = FALSE)

  rownames(selection)=rownames(counts)
  attributes(selection)$dist.mode = d.m
  attributes(selection)$dist.sd = d.sd

  if (plot) {
    hist(d, breaks = 20)
    abline(v=d.m,lty=2,col="green")
    abline(v=threshold,lty=2,col="red")
  }

  return(selection)
}


pcr_saturation <- function(reads, nsamples=100, nrep=10, q=0,from = 0) {
  n_tot  <- sum(reads)
  s_size <- floor(seq(from=1, to=n_tot,length.out = nsamples))
  rep    <- matrix(0.0,nrow = nsamples, ncol = nrep)
  i_from <- floor(nsamples * from) + 1
  for (i in i_from:nsamples)
    for (j in seq_len(nrep)) {
      rep[i,j] <- D_q(ROBITools::rarefy(reads,s_size[i])/s_size[i],q=q)
    }
  n_lm = max(floor((1-0.1) * nsamples),i_from)
  a = lm(rowMeans(rep)[n_lm:nsamples]~s_size[n_lm:nsamples])$coef[2]
  if (is.na(a))
    a=0
  attr(rep,"slope") = a / mean(rep[nsamples,])
  rep
}

