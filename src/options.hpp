// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace RcppParallel;
using namespace Rcpp;

// Global options and parameters all get thrown in here

/*
 // Relevant OPTIONS
 const int mode,
 const int nObs,
 const int nCategories,
 const std::vector<double> priors_altered,
 const std::vector<int>  classFreq,
 const std::vector<std::vector<double>> loss
*/


struct SomaOptionsStructure{

  double pathLength;
  double step;
  double PRT;
  double popSize;
  double nDirections;
  double migrations;
  double minDiv;
  
};

SomaOptionsStructure* SOMA_OPTIONS;



struct OptionsStructure{
  
  
  bool FOLDING; // internal flag to indicate if we're in the main tree or not.
                // Used to skip BFS search during cross-validation because we don't care about
                // node labels then.
  
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
   *      UI AND TRACING CONTROLS
   */
  
  
  int trace; // Level of tracing to provide
  bool DEBUG; // Is low-level debugging enabled?

  int mode; // Mode of operation for DB-CART
            // 0 for classification tree
            // 1 for regression tree
            // 2 for distance-based regression tree
  
  
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
   *      TARGET VARIABLE CONTROLS
   */
            
            
    int nObs; // length of target vector, etc, used for navigating stuff
    int model_nObs; // number of observations in the model itself
    
            
            
  // For classification:
    
    
  int nCategories;
  std::vector<int> classFreq;
  std::vector<double> priors;
  std::vector<double> priors_altered;
  
  std::vector< std::vector<double> > loss;
  int classMethod; // How classification occurs.
                   // 0 for most likely (mode)
                   // 1 for minimum class which gives 50% of cumulative probability (median)


 /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
  *      PREDICTIVE VARIABLE CONTROLS
  */              
  

  // For unordered factor predictors
  int nFactorVars;
  CharacterVector factorNames;
  std::vector<int> levelCount;
  
  // For numeric and ordered factors
  int nNumericVars;
  CharacterVector numericNames;
  std::vector<double> numeric_acceptable_error; 
  
  // For distance based predictors
  int nDistanceVars;
  int nDistanceCombinations;
  CharacterVector distanceNames;
  std::vector<double> distance_acceptable_error; 

  double combinationPower;
  
  bool allCenters; // Do we use all possible centers or just within the node?
  double pCenters; // What proportion of centers should we use? (Will round up)
  
  
  
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
   *    STOPPING METHOD AND OTHER TREE PROPERTIES
   */
  

  int stopMethod; // How do we stop the split process?
                  // 0 = Cross-validation
                  // 1 = Hypothesis Testing

  int nSims; // Permutations for hypothesis test
  double signifLevel; // Significance level for hypothesis test
  
  
  // Standard CART controls
  int maxDepth;
  int minBucket;
  int minSplit;
  
  OptionsStructure()
  {
    
  };

  ~OptionsStructure()
  {
    
  };
  
};

OptionsStructure* OPTIONS;
