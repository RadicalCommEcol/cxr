#functional response type 2
function2 <-function(a,b,x){
  vector <-vector("numeric",length(x))
  for(i in 1:length(x)){
    if (x[i]!=0){vector[i]<- a*(1-exp(-b*x[i]))
             }
  }
  return(vector)
}

#functional response type 3
function3 <-function(a,b,x){
  vector <-vector("numeric",length(x))
  for(i in 1:length(x)){
    if (x[i]!=0){vector[i]<- (a*x[i]^2/(b+x[i]^2))}
  }
  return(vector)
}
