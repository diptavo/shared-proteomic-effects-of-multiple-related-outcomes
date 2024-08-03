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
jk=function(matrices)
{
  library(r.jive)
  rr=jive(matrices,maxiter=20,showProgress="False")
  jr=rr$rankJ
  jm=do.call(rbind, rr$joint)#rbind the list of matrices.
  result=list(jr,jm)
  return(result)
}


#Extracting the proteins (Step3)
extract=function(pro_name,jefs,lower,upper)
{
  library(PMA)
  h2=jefs[[2]] #rbind the list of matrices.
  h_sub2=h2
  pro.list <- vector("list", length = jefs[[1]])
  n1=pro_name
  c1=lower
  c2=upper
  for(i in 1:(jefs[[1]]))
  {
    pp=length(n1)
    l=104
    z=1
    while (l<c1 || l>c2 && z< pp^(0.25)) {
      tt1=SPC(h_sub2,sumabsv=z, K=1)
      l=length(tt1$v[tt1$v!=0])
      z=z+0.1
    }
    ind=which(tt1$v!=0) # indices of protein previously selected
    pro.list[[i]]=n1[ind]
    x=(tt1$d)*(tt1$u)%*%t((tt1$v))
    hs1=h_sub2-x
    h_sub2=hs1[,-ind]
    n1=n1[-ind]
  }
  return(pro.list)
}

#overall procedure
effProtein=function(matrices,lowl=30,upl=100)
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

  jx=jk(hm)

  protein_list=extract(colnames(hm[[1]]),jx,lowl,upl)
  return(protein_list)
}
