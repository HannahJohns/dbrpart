library(dbrpart)


unit_test_soma(function_id = 0,seed = 1,
                     pathLength = 2,step = 0.11,
                     PRT = 0.5,
                     popSize = 10,nDirections = 10,
                     migrations = 30,minDiv = 1e-4,dims = 0)


x1 <- unit_test_soma(function_id = 0,seed = 1,
                     pathLength = 2,step = 0.11,
                     PRT = 0.5,
                     popSize = 10,nDirections = 10,
                     migrations = 30,minDiv = 1e-4,dims = 6)

x1 <- unit_test_soma(function_id = 0,seed = 1,
                     pathLength = 2,step = 0.11,
                     PRT = 0.5,
                     popSize = 10,nDirections = 10,
                     migrations = 30,minDiv = 1e-4,dims = 6)

x2 <- unit_test_soma(function_id = 0,seed = 1,
                     pathLength = 2,step = 0.11,
                     PRT = 0.5,
                     popSize = 10,nDirections = 10,
                     migrations = 30,minDiv = 1e-4,dims = 6)


x3 <- unit_test_soma(function_id = 0,seed = 2,
                     pathLength = 2,step = 0.11,
                     PRT = 0.5,
                     popSize = 30,nDirections = 30,
                     migrations = 30,minDiv = 1e-4,dims = 6)


sapply(1:100,function(x){
               min(unit_test_soma(function_id = 0,seed = x,
               pathLength = 2,step = 0.11,
               PRT = 0.5,
               popSize = 30,nDirections = 30,
               migrations = 30,minDiv = 1e-4,dims = 6)$utility)
})

x3$utility

yx1$pars-x2$pars

x1$pars-x3$pars



{
sink(file="log.txt")
unit_test_soma(function_id = 0,seed = 1,
               pathLength = 2,step = 0.11,
               PRT = 0.5,
               popSize = 10,nDirections = 10,
               migrations = 30,minDiv = 1e-4,dims = 2)
sink()
}
