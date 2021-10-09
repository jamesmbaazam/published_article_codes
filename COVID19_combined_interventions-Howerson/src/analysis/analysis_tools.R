#### ANALYSIS TOOLS ####

# out_sub: data.frame, columns containing x y and z vars
# xvar: character, colname of x variable
# yvar: character, colname of y variable
# zvar: character, colname of z variable
#lvsl: z values to create contours fsor
# levels = c(500,1000,2000)
contour = function(out_sub,xvar, yvar, zvar, lvls){
  x = sort(unique(out_sub[,xvar]))
  y = sort(unique(out_sub[,yvar]))
  z = matrix(0,length(x),length(y))
  for(i in 1:length(x)){
    for(j in 1:length(y)){
      z[i,j] = out_sub[out_sub[,xvar] == x[i] & out_sub[,yvar] == y[j], zvar] 
    }
  }
  cntr = contourLines(x,y,z, levels = lvls)
  cntr = lapply(cntr, function(x){do.call(cbind, x)})
  cntr = do.call(rbind, cntr)
  return(cntr)
}