soma.control <- function(
  seed = NULL, pathLength = 1.5,step = 0.51,
  PRT = 0.1,
  popSize = 20,nDirections = 1,
  migrations = 10,minDiv = 1e-4
)
{
  
  if(is.null(seed))
  {
    seed <- runif(1,min=0,max=.Machine$integer.max)
  }
  
  out <- list(
    seed = seed,
    pathLength = pathLength,
    step=step,
    PRT = PRT,
    popSize = popSize,
    nDirections = nDirections,
    migrations = migrations,
    minDiv = minDiv
  )
  
  return(out)
}