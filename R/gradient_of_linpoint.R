#' Numerically approximate gradient
#' @export

gradient_of_linpoint <- function(linpoint,smesh, tmesh){
  ns = smesh$n; nt = tmesh$n
  coords <- smesh$loc[,c(1,2)]
  distances <- as.matrix(dist(coords, upper = T))
  near.neighbours <- apply(distances, 2, order)[2:(min(10,ns)),]
  grad <- matrix(nrow = ns*nt, ncol = 2)
  for(i in 1:ns){
    diffmat <- matrix(c(coords[near.neighbours[1,i],1]- coords[i,1], 
                        coords[near.neighbours[1,i],2]- coords[i,2],
                        coords[near.neighbours[2,i],1]- coords[i,1], 
                        coords[near.neighbours[2,i],2]- coords[i,2]),
                      byrow = T, nrow = 2)
    diffmat[which(abs(diffmat) < .Machine$double.eps, arr.ind = T)] <- 0
    if(abs(Matrix::det(diffmat))<=.Machine$double.eps){ # if both nearest neighbours are exactly horizontal or both vertical from point, then go to 
      #1st and 3rd near neighbours
      j <-3
      diffmat2 <- matrix(c(coords[near.neighbours[1,i],1]- coords[i,1], 
                           coords[near.neighbours[1,i],2]- coords[i,2],
                           coords[near.neighbours[j,i],1]- coords[i,1], 
                           coords[near.neighbours[j,i],2]- coords[i,2]),
                         byrow = T, nrow = 2)
      while(abs(det(diffmat2)) <= .Machine$double.eps){
        j <- j+1
        diffmat2 <- matrix(c(coords[near.neighbours[1,i],1]- coords[i,1], 
                             coords[near.neighbours[1,i],2]- coords[i,2],
                             coords[near.neighbours[j,i],1]- coords[i,1], 
                             coords[near.neighbours[j,i],2]- coords[i,2]),
                           byrow = T, nrow = 2)
        if(j == 9){
          warning(paste("Mesh behaving strangely. All nearest points to point", i, "lie on a straight line."))
          #browser()
          break
        }
      }
      diffmat2 <- matrix(c(coords[near.neighbours[1,i],1]- coords[i,1], 
                           coords[near.neighbours[1,i],2]- coords[i,2],
                           coords[near.neighbours[j,i],1]- coords[i,1], 
                           coords[near.neighbours[j,i],2]- coords[i,2]),
                         byrow = T, nrow = 2)
      diffmat2[which(abs(diffmat2) < .Machine$double.eps, arr.ind = T)] <- 0
      #print(det(diffmat2))
      for(t in 0:(nt-1)){
        grad[t*ns+i,] <- solve(diffmat2,
                               c(linpoint[near.neighbours[1,i]+t*ns] - linpoint[i + t*ns],
                                 linpoint[near.neighbours[j,i]+t*ns]- linpoint[i + t*ns]))
      }
    }else{
      for(t in 0:(nt-1)){
        grad[t*ns+i,] <- solve(diffmat,
                               c(linpoint[near.neighbours[1,i]+t*ns] - linpoint[i + t*ns],
                                 linpoint[near.neighbours[2,i]+t*ns]- linpoint[i + t*ns]))
      }
    }
    }
  return(grad)
}
