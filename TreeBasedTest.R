# Title: Multiple Sample Metagenomic Hypothesis Testing on a Phylogenetic Tree
# Version: 0.0.1
# Authors: Jun Chen (chen.jun2@mayo.edu)

# Description: Two-sample test for metagenomic sequencing samples, where the sequence
# reads are distributed on the phylogenetic tree. Reads could be distributed on the leaf nodes
# only or on both leaf nodes and internal nodes.  The latter could happen when the reads
# could not be assigned to the leaf nodes due to the incompleteness of the reference tree.

# The input of the data are the edges of the tree, the counts on the edges (equivalently, the 
# counts assigned to the node below), the lengths of the edges, and the group indicator variable

# Date: 2016/06/12

#############################################################################################
# Main function
require(CompQuadForm)

KRtest <- function (edge, edge.ct, edge.len, grp.ind) {
	# Args:
	#       edge: a two-column matrix of mode numeric where each row represents an edge of the tree; 
	#              the nodes and the tips are symbolized with numbers; the tips are numbered 1, 2, ...,
	#              and the nodes are numbered after the tips. For each row, the first column gives 
	#              the ancestor. See the class "phylo" from "ape" package for more details.
	#       edge.ct: a numeric matrix of the counts on the edges, row - edges, column - samples
	#       edge.len: a numeric vector giving the lengths of the branches given by edge. 
	#       grp.ind: a character/numeric vector or a factor indicating the group assignment.
	#               Currently only two groups are allowed.
	#
	# Returns:
	# 	   A numeric value (p value)
	#
	#
	nSeq <- colSums(edge.ct)
	edge.prop <- t(t(edge.ct) / nSeq)
	
	Fstat.obj <- comput.F.stat(edge, edge.prop, edge.len, nSeq, grp.ind, alpha=2)
	Fstat <- Fstat.obj$Fstat
	cum.prop <- Fstat.obj$cum.prop
	
	# Estimate the population mean
	P0 <- rowMeans(cum.prop)
	
	P0[P0 == 0] <- min(P0[P0 != 0]) / 2
	
	cov.mat <- comput.cov.mat(P0, edge)
	lambdas <- comput.lambdas(cov.mat, nSeq, grp.ind)
	
	num.obj <- comput.var.ZPQ(lambdas[['lambda1']], edge.len, grp.ind)
	den.obj <- comput.var.PQ(lambdas[['lambda2']], edge.len, P0, grp.ind)
	
#	Fstat.mod <- ((Fstat.obj$distPQ - num.obj$b) / num.obj$a) /  ((Fstat.obj$P0.v + Fstat.obj$Q0.v - den.obj$b) / den.obj$a) / num.obj$d * den.obj$k
	q <- Fstat * den.obj$b - num.obj$b
	Pval <- davies(q, lambda=c(num.obj$a, - Fstat * den.obj$a), h=c(num.obj$d, den.obj$k))$Qq
	
	return(Pval)
}

#############################################################################################
# Internal function
comput.prop.cum <- function (edge, tab) {
	n <- ncol(tab)
	nbr <- nrow(edge)	
	edge2 <- edge[, 2]
	cum <- matrix(0, nbr, n)
	for (i in 1:nbr) {
		cum[i, ] <- cum[i, ] + tab[i, ]
		node <- edge[i, 1]
		node.loc <- which(edge2 == node)
		while (length(node.loc)) {
			cum[node.loc, ] <- cum[node.loc, ] + tab[i, ]
			node <- edge[node.loc, 1]
			node.loc <- which(edge2 == node)
		}
	}
	cum
}

comput.KR.dist0 <- function (cum1, cum2, br.len, alpha=2, normalize=FALSE) {
	# Comput the distance between two samples
	ind <- (cum1 != 0) | (cum2 != 0)
	cum1 <- cum1[ind]
	cum2 <- cum2[ind]
	br.len2 <- br.len[ind]
	if (!normalize) {
		res <- sum(abs(cum1 - cum2)^alpha * br.len2)^(1/alpha)
	} else {
		res <- sum(abs(cum1 - cum2)^alpha * br.len2)^(1/alpha) / sum(abs(cum1 + cum2)^alpha * br.len2)^(1/alpha)
	}
    return(res)
}

comput.KR.dist <- function (edge, edge.prop, edge.len, is.prop=TRUE, alpha=2, verbose=T) {
	# Denote n = number of sample, 
	#        m = number of leaves/nodes (except the root node), equal to edge number
	# edge: m * 2 matrix containing the node IDs for each edge, first column higher level node (close to root)
	# edge.prop:  m * n matrix containing the number or proportion of mapped reads for each leaf/node 
    # edge.len: n vector containing edge lengths
	# alpha: L-alpha norm 
	# Reads mapped to the root are not considered
	tab <- edge.prop
	br.len <- edge.len
	col.sum <- colSums(tab)
	if (nrow(tab) != nrow(edge)) {
		stop("Number of columns in 'tab' should be equal to number of rows in 'edge'!\n")
	}

	if (sum(col.sum == 0)) {
		if (verbose == TRUE) {
		    cat("Some samples have no reads! They will be dropped!\n")
		}
		tab <- tab[, col.sum != 0, ]
	}
	
	if (is.prop != TRUE) {
		if (verbose == TRUE) {
			cat("Counts are converted into proportions!\n")
		}
		tab <- t(t(tab) / colSums(tab))  
	}

    cum <- comput.prop.cum(edge, tab)
	

	KR.dist <- KR.dist.n <- matrix(0, n, n)
	for (i in 2:n) {
		for (j in 1:(i - 1)) {
			cum1 <- cum[, i]
			cum2 <- cum[, j]
			ind <- (cum1 != 0) | (cum2 != 0)
			cum1 <- cum1[ind]
			cum2 <- cum2[ind]
			br.len2 <- br.len[ind]
			
			KR.dist[i, j] <- KR.dist[j, i] <-  sum(abs(cum1 - cum2)^alpha * br.len2)^(1/alpha)
			KR.dist.n[i, j] <- KR.dist.n[j, i] <- KR.dist[i, j] / sum(abs(cum1 + cum2)^alpha * br.len2)^(1/alpha)
		}
	}

	return(list(KR.dist=KR.dist, KR.dist.n=KR.dist.n, cum.prop=cum))
}


comput.F.stat <- function (edge, edge.prop, edge.len, nSeq, grp.ind, alpha=2) {
	# Start with 2 groups with equal numbers
	# Use unnormalized KR.dist
	ind1 <- as.numeric(factor(grp.ind)) == 1
	ind2 <- as.numeric(factor(grp.ind)) == 2
	m <- sum(ind1)
	n <- sum(ind2)
    
	nSeq1 <- nSeq[ind1]
	nSeq2 <- nSeq[ind2]
	cum.prop <- comput.prop.cum(edge, edge.prop)
	
	P0.cum <- cum.prop[, ind1, drop=FALSE]
	Q0.cum <- cum.prop[, ind2, drop=FALSE]
	
	P0 <- rowMeans(cum.prop)	
#	P0 <- rowSums(t(t(cum.prop) * nSeq)) / sum(nSeq)   # squash

	P0.cum <- t(t(P0.cum - P0) * sqrt(nSeq1))
	Q0.cum <- t(t(Q0.cum - P0) * sqrt(nSeq2))
	
	P0.m <- rowMeans(P0.cum)
	Q0.m <- rowMeans(Q0.cum)
	
	P0.v <- sum(apply(P0.cum, 2, comput.KR.dist0, cum2=P0.m, br.len=edge.len)^2) 
	Q0.v <- sum(apply(Q0.cum, 2, comput.KR.dist0, cum2=Q0.m, br.len=edge.len)^2) 
	
	distPQ <- comput.KR.dist0(P0.m, Q0.m, edge.len)^2
	
	F <- distPQ / (P0.v + Q0.v)
	
	return(list(Fstat=F, cum.prop=cum.prop, distPQ=distPQ, P0.v=P0.v, Q0.v=Q0.v))	
}


comput.cov.mat <- function (P0, edge) {
	nbr <- nrow(edge)	
	edge2 <- edge[, 2]
	subset.mat <- matrix(0, nbr, nbr)
	diag(subset.mat) <- 1:nbr
	for (i in 1:nbr) {
		node <- edge[i, 1]
		node.loc <- which(edge2 == node)
		while (length(node.loc)) {
			subset.mat[i, node.loc] <- subset.mat[node.loc, i] <- i
			node <- edge[node.loc, 1]
			node.loc <- which(edge2 == node)
		}
	}
	subset.mat[subset.mat != 0] <- P0[subset.mat]
	return(subset.mat - tcrossprod(P0))
}

comput.lambdas <- function (cov.mat, nSeq, grp.ind) {
	ind1 <- as.numeric(factor(grp.ind)) == 1
	ind2 <- as.numeric(factor(grp.ind)) == 2
	m <- sum(ind1)
	n <- sum(ind2)
	
    mat1 <- cov.mat  * (m + n) / m / n
	lambda1 <- eigen(mat1, symmetric=TRUE, only.values=TRUE)$values
	
	mat2 <- cov2cor(cov.mat)
	lambda2 <- eigen(mat2, symmetric=TRUE, only.values=TRUE)$values
	
	return(list(lambda1=lambda1, lambda2=lambda2))
}

comput.var.ZPQ <- function (lambdas, edge.len, grp.ind=NULL) {
#	cov.mat <- comput.cov.mat(P0, edge, nSeq, grp.ind)
#	lambdas <- eigen(cov.mat, symmetric=TRUE, only.values=TRUE)$values
	ll <- edge.len * lambdas
	ll2 <- sum(ll^2)
	ll3 <- sum(ll^3)
	ll1 <- sum(ll)
	a <- ll3 / ll2
	b <- ll1 - ll2^2 / ll3
	d <- ll2^3 / ll3^2
	
	return(list(a=a, b=b, d=d))
}

comput.var.PQ <- function (lambdas, edge.len, P0, grp.ind) {
	
#	cor.mat <- comput.cor.mat(P0, edge)
#	lambdas <- eigen(cor.mat, symmetric=TRUE, only.values=TRUE)$values
	
	PQ0 <- P0 * (1 - P0)
	ind1 <- as.numeric(factor(grp.ind)) == 1
	ind2 <- as.numeric(factor(grp.ind)) == 2
	m <- sum(ind1)
	n <- sum(ind2)
	c1 <- 1 + 1 / (1:(m-1) + 1)
	c2 <- 1 + 1 / (1:(n-1) + 1)
    c3 <- c(c1, c2)

	
	ld1 <- ld2 <- ld3 <- 0
	for (i in 1:length(c3)) {
		lambdas2 <- c3[i] * lambdas - c3[i] + 1
		ld <- PQ0 *edge.len * lambdas2
		ld1 <- ld1 + sum(ld)
		ld2 <- ld2 + sum(ld^2)
		ld3 <- ld3 + sum(ld^3)
	}
	a <- ld3 / ld2
	b <- (ld1 -  ld2^2 / ld3)
	k <- ld2^3 / ld3^2
	
	return(list(a=a, b=b, k=k))
}

#############################################################################################





