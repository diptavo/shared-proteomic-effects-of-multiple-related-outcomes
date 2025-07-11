#function to harmonize the matrices (Step1)
harmonize_matrices <- function(matrices) {
  # Find the common column names
  com_cols <- Reduce(intersect, lapply(matrices, colnames))
  # If no common columns are found
  if (length(com_cols) == 0) {
    message("No common columns found.")
    return(NULL)
  }
  # Subset matrices to include only the common columns
  updated_matrices <- lapply(matrices, function(mat) {
    mat[, com_cols, drop = FALSE]
  })
  return(updated_matrices)
}


#function to run JIVE (Step2)
TARIF.D=function(matrices,mi,status)
{
  library(r.jive)
  rr=jive(matrices,maxiter=mi,showProgress=status)
  jr=rr$rankJ
  jm=do.call(rbind, rr$joint)#rbind the list of matrices.
  result=list(jr,jm)
  return(result)
}


#Extracting the proteins (Step3)
TARIF.S=function(pro_name,jefs,lower,upper,lambda.v,gamma)
{
  library(PMA)
  R=jefs[[2]] #rbind the list of matrices.
  
  pro.list <- vector("list", length = jefs[[1]])
  n1=pro_name
  c1=lower
  c2=upper
  for(i in 1:(jefs[[1]]))
  {
    
    l=c2+1
    z=1
    while (l<c1 || l>c2 && z< lambda.v) {
      tt1=SPC(h_sub2,sumabsv=z, K=1)
      l=length(tt1$v[tt1$v!=0])
      z=z+gamma
    }
    ind=which(tt1$v!=0) # indices of protein previously selected
    pro.list[[i]]=n1[ind]
    x=(tt1$d)*(tt1$u)%*%t((tt1$v))
    R.star=R-x
    R=R.star[,-ind]
    n1=n1[-ind]
  }
  return(pro.list)
}

#overall procedure
TARIF=function(matrices,mi=20,status="False",lowl=30,upl=100,lambda.v,gamma)
{
  # Check if the input is a list of matrices
  if (!all(sapply(matrices, is.matrix))) {
    stop("All inputs must be matrices.")
  }
  
  hm=harmonize_matrices(matrices)
  
  if (is.null(hm)) {
    stop("No common columns found.")
  }
  #pro_name=colnames(hm[[1]])
  
  jx=TARIF.D(hm,mi,status)
  
  protein_list=TARIF.S(colnames(hm[[1]]),jx,lowl,upl,,lambda.v,gamma)
  return(protein_list)