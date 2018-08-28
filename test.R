ls <- sapply(seq_along(population), function(k) sapply(seq_along(population), function(p) sum(population[[p]])) %*% contactMatrix[,k]*sum(population[[k]])/sum(unlist(population)))
sum(ls)

sapply(seq_along(population), function(p) sum(population[[p]])) %*% contactMatrix[,k]/
