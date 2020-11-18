#include "rcknntree.hpp"

/*
 *  UNIT TESTS
 */

// [[Rcpp::export]]
List unit_test_soma(int function_id, int seed,
                    double pathLength, double step, double PRT,
                    int popSize, int nDirections,
                    double migrations, double minDiv,
                    int dims)
{
  
  /* Function list:
   
   See the following: https://en.wikipedia.org/wiki/Test_functions_for_optimization
   
   Implemented benchmarks functions are:
   
   0 Rangstrigin function
   1 Sphere function
   2 Rosenbrock function
   
   */

  // std::srand(seed);
  
  // Set some constraints
  std::vector<double> parLb;
  std::vector<double> parUb;
  for(int i=0; i<dims; i++)
  {
    switch(function_id)
    {
    case 0: parLb.push_back(-5.12);
      parUb.push_back(5.12);
      break;
    case 1: parLb.push_back(-10);
      parUb.push_back(10);
      break;
    case 2: parLb.push_back(-10);
      parUb.push_back(10);
      break;
    }
  }
  
  Soma soma(pathLength,
            step,
            PRT,
            popSize,
            nDirections,
            migrations,
            minDiv,
            dims,
            0,
            parLb,
            parUb 
  );
  

  while(soma.exitCode <0)
  {
    
    soma.move();
    //soma.printStatus();
    
    for(int i=0; i< soma.nTest; i++)
    {
     
      double value=0;
    
      if(dims==0)
      {
       value = -99; 
      }
      else
      {
        if(function_id==0)
        { 
          value = 10 * dims;
          for(int j=0; j<dims; j++)
          {
            value = value + (soma.pars[i][j]*soma.pars[i][j] - 10 * std::cos(2* 3.141592653589793238462643383280 * soma.pars[i][j]) );
          }
        }
        else if(function_id==1)
        {
          for(int j=0; j<dims; j++)
          {
            value = value + (soma.pars[i][j]*soma.pars[i][j]);
          }
        }
        else if(function_id==2)
        {
          for(int j=0; j<dims; j++)
          {
            value = value + (soma.pars[i][j]*soma.pars[i][j]);
          }
        }
      }
    
      soma.utility[i] = value;
    
    }
    
    soma.update();
    
  }
  
  // Wrap everything up
  
  Rcpp::NumericVector utility(popSize);
  NumericMatrix parameters(popSize,dims);
    
  if(dims>0)
  {
  for(int i=0; i<popSize; i++)
  {
    utility[i] = soma.utility_current[i];
    for(int j=0;j<dims;j++)
    {
      parameters[i + j*popSize] = soma.pars_current[i][j];
    }
  } 
  }
  else
  {
    utility[0] = soma.utility_current[0];
  }
    
  
  List out = List::create(Named("statusCode") = soma.exitCode,
                          Named("pars") = parameters,
                          Named("utility") = utility);
    
  return out;
  
  
  
}


// [[Rcpp::export]]
List fit_dbrpart(IntegerVector in_target_discrete,
                NumericVector in_target_continuous,
                NumericMatrix in_target_distance,
                List in_factorList,
                List in_numericList,
                List in_distanceList,
                List in_distanceCombinations,
                double in_combinationPower,
                NumericVector in_acceptable_error,
                int in_mode,
                int in_stopMethod,
                int in_nSims,
                double in_signifLevel,
                IntegerVector in_fold,
                int in_nCategories,
                NumericVector in_priors,
                NumericMatrix in_loss,
                int in_minSplit,
                int in_minBucket,
                int in_maxDepth,
                bool in_alteredPriors,
                int in_trace,
                int in_classMethod,
                int in_pCenters,
                double soma_pathLength,
                double soma_step,
                double soma_PRT,
                double soma_popSize,
                int soma_nDirections,
                int soma_migrations,
                double soma_minDiv,
                int soma_seed
                )
{
  
  Rcout << "a";
  
  std::srand(soma_seed);
  
  SOMA_OPTIONS = new SomaOptionsStructure();
  
  SOMA_OPTIONS->pathLength = soma_pathLength;
  SOMA_OPTIONS->step = soma_step;
  SOMA_OPTIONS->PRT = soma_PRT;
  SOMA_OPTIONS->popSize = soma_popSize;
  SOMA_OPTIONS->nDirections = soma_nDirections;
  SOMA_OPTIONS->migrations = soma_migrations;
  SOMA_OPTIONS->minDiv = soma_minDiv;
  
  Rcout << "b";
  
  OPTIONS = new OptionsStructure();
  
  // Next step: Allow for conventional split options too!
  OPTIONS->nDistanceVars=in_distanceList.length();
  OPTIONS->nDistanceCombinations=in_distanceCombinations.length();
  
  CharacterVector distanceNames = in_distanceList.names();
  for (int i=0; i< OPTIONS->nDistanceVars;i++)
  {
    OPTIONS->distanceNames.push_back(distanceNames[i]);
  }
  
  switch (in_mode)
  {
  case 0: OPTIONS->nObs = in_target_discrete.length();
    break;
  case 1: OPTIONS->nObs = in_target_continuous.length();
    break;
  case 2: OPTIONS->nObs = in_target_distance.nrow();
    break;
  default: Rcout << std::endl << "Mode not recognised! Aborting!" << std::endl;
  }
  
  std::vector<int> idList;
  for(int i=0;i<OPTIONS->nObs;i++)
  {
    idList.push_back(i);
  }
  
  Rcout << "c";
  
  ///////////////////////////////////////////////////////////////////////
  //
  //                    Target variable
  //
  ///////////////////////////////////////////////////////////////////////
  
  
  OPTIONS->mode = in_mode;
  
  std::vector<int> target_discrete;
  if(OPTIONS->mode==0){
    for(int i=0;i<OPTIONS->nObs;i++)
    {
      target_discrete.push_back(in_target_discrete[i]);
    }
  }
  
  std::vector<double> target_continuous;
  if(OPTIONS->mode==1){
    for(int i=0;i<OPTIONS->nObs;i++)
    {
      target_continuous.push_back(in_target_continuous[i]);
    }
  }
  
  std::vector< std::vector<double> > target_distance;
  if(OPTIONS->mode==2){
    
    for(int i=0;i<OPTIONS->nObs;i++)
    {
      std::vector<double> thisRow;
      for(int j=0;j<OPTIONS->nObs;j++)
      {
        thisRow.push_back(in_target_distance[i+j*OPTIONS->nObs]);
      }
      target_distance.push_back(thisRow);
    }
  }
  
  Rcout << "d";
  
  
  ///////////////////////////////////////////////////////////////////////
  //
  //                    Predictors
  //
  ///////////////////////////////////////////////////////////////////////
  
  
  std::vector< std::vector< std::vector<double> > > predictor_distance;
  for(int i=0;i<OPTIONS->nDistanceVars;i++)
  {
    NumericMatrix in_thisD = in_distanceList[i];
    
    std::vector< std::vector<double> > thisD;
    for(int j=0;j<OPTIONS->nObs;j++)
    {
      std::vector<double> thisRow;
      for(int k=0;k<OPTIONS->nObs;k++)
      {
        thisRow.push_back(in_thisD[k+j*OPTIONS->nObs]);
      }
      thisD.push_back(thisRow);
    }
    predictor_distance.push_back(thisD);
  }

  OPTIONS->combinationPower = in_combinationPower;
  
  std::vector< std::vector<int> > predictor_distance_combinations;  
  for(int i=0; i<OPTIONS->nDistanceCombinations;i++)
  {
    IntegerVector in_thisCombination = in_distanceCombinations[i];
    
    std::vector<int> thisCombination;
    for(int j=0; j<in_thisCombination.size(); j++)
    {
      thisCombination.push_back(in_thisCombination[j]);
    }
    predictor_distance_combinations.push_back(thisCombination);
  }  
  
  
  Rcout << "e";
  
  // Pass options to a global variable to keep them all in one place
  // Matrices are dropped to std::vector<> format to play nice with
  // multithreadded code later on.
  // OPTIONS should not be written to after initial setup.
  
  OPTIONS->trace = in_trace;
  OPTIONS->classMethod = in_classMethod;
  OPTIONS->maxDepth = in_maxDepth;
  OPTIONS->minBucket = in_minBucket;
  OPTIONS->minSplit = in_minSplit,
  OPTIONS->nCategories = in_nCategories;
  
  OPTIONS->mode=in_mode;
  OPTIONS->stopMethod=in_stopMethod;
  
  // These should be modifiable in R
  OPTIONS->nSims = in_nSims;
  OPTIONS->signifLevel = in_signifLevel;
  
  OPTIONS->pCenters = in_pCenters;
  OPTIONS->allCenters = false; // FLAG
  
  if(OPTIONS->trace == 999)
  {
    // This will cause race conditions in RcppParallel unless
    // RcppParallel is set to 1 thread
    OPTIONS->DEBUG=true;
  }
  else
  {
    OPTIONS->DEBUG=false;
  }
  
  for(int i=0;i<OPTIONS->nDistanceVars;i++)
  {
    OPTIONS->distance_acceptable_error.push_back(in_acceptable_error[i]);
  }
  
  Rcout << "f";
  
  if(OPTIONS->mode==0)
  {
    // Construct priors vectors and loss matrix
    // with zeroes for the time being
    for(int i=0;i<OPTIONS->nCategories;i++)
    {
      OPTIONS->classFreq.push_back(0);
      OPTIONS->priors.push_back(0);
      OPTIONS->priors_altered.push_back(0);
    }
    
    
    Rcout << "g";
    
    // There's something weird going on where
    // R is reusing the values stored in these vectors
    // between runs, even though they should have been initialised
    // to zero above. I have no idea what the hell the problem is here
    // other than it's some bullshit memory thing or Rcpp not doing garbage
    // collection properly, but we're just gonna force the classfreq scores to zero
    // for real this time
    
    for(int i=0;i<OPTIONS->nCategories;i++)
    {
      OPTIONS->classFreq[i] = 0;
    }
    
    Rcout << "i";
    
    for(int i=0;i<OPTIONS->nCategories;i++)
    {
      OPTIONS->loss.push_back(OPTIONS->priors);
    }
    
    Rcout << "j";
    
    for(int i=0;i<OPTIONS->nObs; i++)
    {
      Rcout << target_discrete[i];
      OPTIONS->classFreq[target_discrete[i]]++;
    }
    
    Rcout << "k";
    
    
    for(int i=0;i<OPTIONS->nCategories;i++)
    {
      OPTIONS->priors[i] = in_priors[i];
      OPTIONS->loss[i][i] = in_loss[i+i*OPTIONS->nCategories];
      for(int j=i+1;j<OPTIONS->nCategories;j++)
      {
        OPTIONS->loss[i][j] = in_loss[i+j*OPTIONS->nCategories];
        OPTIONS->loss[j][i] = in_loss[j+i*OPTIONS->nCategories];
      }
    }
    
    Rcout << "h";
    
    // Calculate altered priors.
    if(in_alteredPriors)
    {
      // Get loss for each category
      NumericVector L(OPTIONS->nCategories);
      double sum_piL = 0;
      
      for(int i=0;i<OPTIONS->nCategories;i++)
      {
        L[i] = 0;
        for(int j=0;j<OPTIONS->nCategories;j++)
        {
          L[i] = L[i]+in_loss(i,j);
        }
        sum_piL = sum_piL + in_priors[i]*L[i];
      }
      
      Rcout << sum_piL;
      
      for(int i=0;i<OPTIONS->nCategories;i++)
      {
        OPTIONS->priors_altered[i] = in_priors[i] * L[i] / sum_piL;
      }
    }
    else
    {
      // If altered priors aren't actually used,
      // then just set altered_pi=pi
      for(int i=0;i<OPTIONS->nCategories;i++)
      {
        OPTIONS->priors_altered[i] = in_priors[i];
      }
    }
    
    Rcout << "i";
    
  }
  else if(OPTIONS->mode==1)
  {
    
    Rcout << "nObs " << OPTIONS->nObs << " ySize " << target_continuous.size() << std::endl;
    
  }
  
  if(OPTIONS->trace>0)
  {
    Rcout << "All values passed in!" << std::endl;
  } 
  
  
  
  
    
    OPTIONS->model_nObs = OPTIONS->nObs;
    OPTIONS->FOLDING = false;
    
    RcknnTree_bfs tree(target_discrete,
                       target_continuous,
                       target_distance,
                       predictor_distance,
                       predictor_distance_combinations,
                       idList);
    
    
    
    // What we return to R varies with the method
    
    List out;
    
    // Are we using the classical method of pruning?
    if(OPTIONS->stopMethod == 0)
    {
      // Before we launch into cross-validation,
      // Get the threshold values for pruning
      // as well as things like tree risk
      // and number of terminals at each cut point
      std::vector<double> alpha = tree.alpha;
      std::vector<int> nTerminal = tree.nTerms;
      std::vector<double> risk = tree.risk;
      
      Rcout << std::endl << "a: ";
      for (int i=0;i<alpha.size();i++)
      {
        Rcout << alpha[i] << " ";
      }
      Rcout << std::endl;
      
      Rcout << std::endl << "n: ";
      for (int i=0;i<nTerminal.size();i++)
      {
        Rcout << nTerminal[i] << " ";
      }
      Rcout << std::endl;
      
      Rcout << std::endl << "r: ";
      for (int i=0;i<risk.size();i++)
      {
        Rcout << risk[i] << " ";
      }
      Rcout << std::endl;
  
      std::vector<double> beta;
      // If there's no alpha value for "keep all sensible splits
      // then force it to be included as an option
      if(alpha.size()>1)
      {
        if(alpha[0] > 0 & alpha[1] > 0) // Allow for -1 i.e. Tmax > T1 
        {
          beta.push_back(0);
        }
      }
      
      for(int i=0;i<(alpha.size()-1);i++)
      {
        // alpha < 0 corresponds to Tmax and will cause nan issues when we
        // take the geometric mean. Skip for now and patch in the next step
        if(alpha[i] >= 0)
        {
          beta.push_back(sqrt(alpha[i] * alpha[i+1]));
        }
      }
      beta.push_back(std::numeric_limits<double>::infinity());
      
      // If 0 wasn't in the alpha list (happens when Tmax == T1 i.e. no meaningless partitions)
      // then we would have skipped the "full tree" option above.
      // If this happened, just stick a 0 at the start of the beta vector.
      if(beta[0] > 0)
      {
        beta.insert(beta.begin(),0);  
      }
     
      Rcout << std::endl << std::endl << "b: ";
      for (int i=0;i<beta.size();i++)
      {
        Rcout << beta[i] << " ";
      }
      Rcout << std::endl;
  
      std::vector< double > predictionRisk;
      std::vector< std::vector<int> > termNode;
      std::vector< std::vector<double> > termPrediction;
      std::vector< std::vector< std::vector<double> > > termPredictionDetails;
  
      // Terminal node label and resulting prediction vary with beta.
      // Initialise the vectors that will store this data, then pass them by reference.
      // (This mode of operation is more important later when cross-validating)
      for(int i=0;i<beta.size();i++)
      {
        
        std::vector<int> thisTermNode(OPTIONS->nObs,-1);
        std::vector<double> thisPrediction(OPTIONS->nObs,-1);
        
        // Details Vector varies in length depending on the mode
        std::vector<double> thisDetailsVector;
        if(OPTIONS->mode==0)
        {
          for(int j=0;j<OPTIONS->nCategories;j++)
          {
            thisDetailsVector.push_back(-1);
          }
        }
        else
        {
          thisDetailsVector.push_back(-1);
        }
        
        std::vector< std::vector<double> > thisDetails;
        for(int j=0; j<OPTIONS->nObs; j++)
        {
          thisDetails.push_back(thisDetailsVector);
        }
        
        predictionRisk.push_back(tree.predictionRisk(target_discrete,
                                           target_continuous,
                                           target_distance, 
                                           predictor_distance,
                                           predictor_distance_combinations,
                                           idList,
                                           i));
      
        tree.classTree(predictor_distance,
                       predictor_distance_combinations,
                       0,
                       idList,
                       i,
                       thisTermNode
                       );
          
        tree.prediction(predictor_distance,
                        predictor_distance_combinations,
                        0,
                        idList,
                        i,
                        thisPrediction
                        );
        
        tree.predictionDetails(predictor_distance,
                               predictor_distance_combinations,
                               0,
                               idList,
                               i,
                               thisDetails,
                               target_distance);
        
        
        termNode.push_back(thisTermNode);
        termPrediction.push_back(thisPrediction);
        termPredictionDetails.push_back(thisDetails);
      }
      
      
      if(OPTIONS->trace>0)
      {
        Rcout << std::endl << "Starting Cross-validation" << std::endl;
      }
    
      // Drop trace by 1 level
      OPTIONS->trace--;
      
      OPTIONS->FOLDING = true;
      
      // Get risk under cross validation, along with the actual predictions made
      // just in case the user wants to apply some other metric instead
      
      
      // Initialise outputs
      
      std::vector<double> cv_risk(beta.size(),0);
      std::vector< std::vector<double> > cv_termPrediction;
      std::vector< std::vector< std::vector<double> > > cv_termPredictionDetails;
      
      for(int i=0;i<beta.size();i++)
      {
        std::vector<double> thisDoubleVector(OPTIONS->nObs,-1);
        
        // Details Vector varies in length depending on the mode
        std::vector<double> detailsVector;
        if(OPTIONS->mode==0)
        {
          for(int j=0;j<OPTIONS->nCategories;j++)
          {
            detailsVector.push_back(-1);
          }
        }
        else
        {
          detailsVector.push_back(-1);
        }
        
        std::vector< std::vector<double> > detailsVectorVector;
        for(int j=0; j<OPTIONS->nObs; j++)
        {
          detailsVectorVector.push_back(detailsVector);
        }
        
        
        cv_termPrediction.push_back(thisDoubleVector);
        cv_termPredictionDetails.push_back(detailsVectorVector);
      }
      
      
      // Build crossvalidated trees
      for( int iFold=0; iFold<10;iFold++)
      {
        OPTIONS->model_nObs = 0;
        
        if(OPTIONS->trace>0)
        {
          Rcout << std::endl << "Fold " <<iFold ;
        }
        
        // Separate data into training and testing subsets
        
        std::vector<int> trainList;
        std::vector<int> testList;
        
        // Set nObs to refer to the training subset
        for(int i=0;i<OPTIONS->nObs;i++)
        {
          if(in_fold[i] == iFold)
          {
            testList.push_back(idList[i]);
          }
          else
          {
            trainList.push_back(idList[i]);
            OPTIONS->model_nObs++;
          }
        }
        
        // Build tree on training subset
        
        if(OPTIONS->trace>0)
        {
          Rcout << std::endl << "Training";
        }
        
        
        
        if(OPTIONS->mode==0)
        {
          // Set class frequencies to reflect the fold.
          // Keep priors the same.
          
          for(int j=0;j<OPTIONS->nCategories; j++)
          {
            OPTIONS->classFreq[j]=0;
          }
          for(int j=0;j<OPTIONS->model_nObs; j++)
          {
            OPTIONS->classFreq[target_discrete[trainList[j]]]++;
          }
        }
        
        
        RcknnTree_bfs cv_tree(target_discrete,
                              target_continuous,
                              target_distance,
                              predictor_distance,
                              predictor_distance_combinations,
                              trainList);
        
        
        
        if(OPTIONS->trace>0)
        {
          Rcout << std::endl << "Testing";
        }
        
        for(int i=0;i<beta.size();i++)
        {
          
          // Find which Tv0, Tv1, ... best corresponds to beta
          // This is whichever is the largest tree
          // such that a >= beta 
          int nSkipped = 0;
          int pruneIndex = 0;
          bool largestTreeFound = false;
          do{
            if(cv_tree.alpha[pruneIndex] < -0.1)
            {
              // If alpha==-1 then skip it, and track the fact that we skip it
              // so we can get the correct T0, T1, ... index
              pruneIndex++;
              nSkipped++;
            }
            if(pruneIndex == (cv_tree.alpha.size()-1))
            {
              largestTreeFound=true;
            }
            else
            {
              if(cv_tree.alpha[pruneIndex+1] < beta[i])
              {
                pruneIndex++;
              }
              else
              {
                largestTreeFound=true;
              }
            }
            
          } while (!largestTreeFound);
          
          
          pruneIndex = pruneIndex-nSkipped;
          
          
          cv_risk[i] = cv_risk[i] + cv_tree.predictionRisk(target_discrete,
                                                           target_continuous,
                                                           target_distance, 
                                                           predictor_distance,
                                                           predictor_distance_combinations,
                                                           testList,
                                                           pruneIndex);
          
          cv_tree.prediction(predictor_distance,
                             predictor_distance_combinations,
                          0,
                          testList,
                          pruneIndex,
                          cv_termPrediction[i]
          );
          
          cv_tree.predictionDetails(predictor_distance,
                                    predictor_distance_combinations,
                                 0,
                                 testList,
                                 pruneIndex,
                                 cv_termPredictionDetails[i],
                                 target_distance);
          
        }
      }
      
      
      
      // Restore trace
      OPTIONS->trace++;
      
      if(OPTIONS->trace>0)
      {
        Rcout << std::endl << "Preparing exports to R" << std::endl;
      }  
      
      // Export the tree and return to R
      List treeExport = tree.exportTree();
  
      
      // Flatten out nested vectors into R-friendly export formats
      IntegerMatrix out_termNode(OPTIONS->nObs,beta.size());
      
      NumericMatrix out_termPrediction(OPTIONS->nObs,beta.size());
      List out_predictionDetails;
      
      NumericMatrix out_cv_termPrediction(OPTIONS->nObs,beta.size());
      List out_cv_predictionDetails;
      
      int nDetails = OPTIONS->mode==0 ? OPTIONS->nCategories : 1;
      for(int i=0;i<beta.size();i++)
      {
        NumericMatrix this_predictionDetails(OPTIONS->nObs, nDetails);
       NumericMatrix this_cv_predictionDetails(OPTIONS->nObs, nDetails);
        
        for(int j=0;j<OPTIONS->nObs;j++)
        {
          out_termNode[j+i*OPTIONS->nObs] = termNode[i][j];
          out_termPrediction[j+i*OPTIONS->nObs] = termPrediction[i][j];
          out_cv_termPrediction[j+i*OPTIONS->nObs] = cv_termPrediction[i][j];
          
          for(int k=0;k<nDetails; k++)
          {
            this_predictionDetails[j+k*OPTIONS->nObs] = termPredictionDetails[i][j][k];
            this_cv_predictionDetails[j+k*OPTIONS->nObs] = cv_termPredictionDetails[i][j][k];
          }
        }
        
        out_predictionDetails.push_back(this_predictionDetails);
       out_cv_predictionDetails.push_back(this_cv_predictionDetails);
      }
      
      out = List::create(Named("a") = wrap(tree.alpha),
                         Named("b") = wrap(beta),
                         Named("nTerm") = wrap(nTerminal),
                         Named("model_risk") = wrap(risk),
                         Named("risk") = wrap(predictionRisk),
                         Named("cv_risk") = wrap(cv_risk),
                         Named("classTree") = out_termNode,
                         Named("prediction") = out_termPrediction,
                         Named("details") = out_predictionDetails,
                         Named("cv_prediction") = out_cv_termPrediction,
                         Named("cv_details") = out_cv_predictionDetails,
                         Named("tree") = treeExport
                        );
      
    }
    else
    {
      List treeExport = tree.exportTree();
      out = List::create(Named("tree") = treeExport);
    }
    
    if(OPTIONS->trace>0)
    {
      Rcout << std::endl << "Method complete" << std::endl;
    }  
    
    delete OPTIONS;
    
    if(OPTIONS->trace>0)
    {
      Rcout << std::endl << "Returning to R" << std::endl;
    }  
    return out;
}

