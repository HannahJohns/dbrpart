#include "bubble.hpp"
#include <queue>

using namespace RcppParallel;
using namespace Rcpp;


/*
 * DFS definition of rcknnTree appears to be causing issues at the end
 * of the function, so instead we're going to run this in a BFS approach.
 */

struct RcknnTree_bfs{
  
  // Store each node separately.
  // RcknnTree_bfs contains a std::vector<*Node> variable
  // which is used to access and navigate these.
  
  struct Node{
    
    // Node label.
    // NOTE: This assumes no pruning takes place,
    // i.e. label is equivalent to the poisition in
    // the vector given the tree structure.
    // The actual label will need to be calculated in R
    // and is dependent on the level of pruning used.
    int label;

    // Depth of the node in the tree
    int depth;
    
    double probability; // Probability of being in this node
    
    // Members of this node
    std::vector<int> groupMembers;
    int groupSize;
    
    
    // If this node were to be terminal, what's the decision we'd make?
    
    double risk;     // Risk of misclassification of the node
    double impurity; // Gini or SST, depending on mode (=risk for regression)
    
    // If mode==0 i.e. classification:
    std::vector<int> counts; //Count number of each classification in node
    std::vector<double> pClassGivNode; // Probability of classification given this node
    int classif;
    
    // If mode==1 i.e. regression:
    double mean;
    
    // Split variables
    
    bool terminal; // Do we stop here?
    
    int splitMethod; // Right now does nothing but will be either 0=factor, 1=numeric or 2=distance
    
    int bestSplitIndex; // Which split was the best?
    std::vector<double> weights; // If there's a combination of splits, what's the relative weight? (Sums to 1)
    
    std::string splitName; // Name of the split variable
    double splitUtility; // Split quality (Gini or SST). Note: impurity - splitutility is goodness of fit
    
    
    std::vector<double> goodnessOfFit; // Store goodness of fit for each possible split.
                                       // Allows us to see what variables were "runners up", etc
                                       
    
    // For factor split (splitMethod==0)
    // Not implemented yet
    
    // For numeric split (splitMethod==1)
    // Not implemented yet
    
    // For distance based split (splitMethod==2)
    int center;
    double radius;
    
    // If we're running with perutation tests, what was the p-value for the split?
    double pVal;

   
    
    
    // Labels for in/out children
    int inChild;
    int outChild;
    
        
    // Pruning variables
    // NOTE: all of these need to be calculated recursively.
    // We can't get them until BFS is finished.
   
    // Flags the order in which prunes happen in the tree.
    // Could just be a bool about if this node has been pruned or not,
    // but the order in which the tree was pruned is probably useful to know
    // In terms of the tree's growth.
    int pruneAtStep;
    
    
    // NOTE: These values are largely meaningless
    // outside of the tree. They are used for
    // identifying the weakest link to prune
    // in the tree, but are overridden whenever we make a new cut.
    
    int nTerminal; // Count number of terminal/pruned nodes below this node
    double subTreeRisk; // If terminal, prob*risk, otherwise sum of prob*risk below
    double criticalAlpha; // Critical Value for pruning tree at this node (Breiman, Friedman and Olshen 1984, Classification and Regression Trees, p. 69)
    
    // Constructor for node
    Node(std::vector<int> &target_discrete,
         std::vector<double> &target_continuous,
         std::vector< std::vector<double> > &target_distance,
         std::vector< std::vector< std::vector<double> > > &D_list,
         std::vector< std::vector<int> > &predictor_distance_combinations,
         std::vector<int> idList,
         int label,
         int depth
         )
    {
      
      // Placeholder value to be calculated later.      
      this -> pruneAtStep = -1; 
      
      
      
      this->label = label;
      
      
      // No children upon construction
      this->inChild=-1;
      this->outChild=-1;
      
      if(OPTIONS->trace>1)
      {
        Rcout << std::endl << "Constructing Node " << this->label << " ";
      }
      
  
      /*****************************************************
       *
       *  Node value input and calculations
       *
       */
      
      // Initialise the node
      this->depth = depth;
      this->terminal=true; // Default to terminal node
      
      // Feed in the id list
      this->groupSize = idList.size();
      for(int i=0; i<(this->groupSize);i++)
      {
        this->groupMembers.push_back(idList[i]);
      }
      
      // Feed in tallies, means, etc.
      if(OPTIONS->mode==0)
      {
        
        // Set up a classification tree node
        for(int i=0;i<OPTIONS->nCategories;i++)
        {
          this->counts.push_back(0);
        }
        
        for(int i=0; i<(this->groupSize);i++)
        {
          this->counts[target_discrete[idList[i]]]++;
        }
        
        // With tally done, we can calculate probabilities,
        // risk, etc for this node
        
        // With no prior split to tell us the result,
        // we need to calculate gini/risk/etc
        
        // Probability of visiting this node
        this->probability = 0;
        for(int i=0;i<OPTIONS->nCategories;i++)
        {
          this->probability += OPTIONS->priors[i]*(this->counts[i])/OPTIONS->classFreq[i];
        }
        
        //Probability of each outcome given in this node
        for(int i=0;i<OPTIONS->nCategories;i++)
        {
          this->pClassGivNode.push_back(
              (OPTIONS->priors[i]*(this->counts[i])/OPTIONS->classFreq[i]) /
                (this->probability)
          );
        }
        
        // With probability calculated, we can make a classification
        
        double classProb = 0;
        this->classif = -1;
        for(int i=0;i<OPTIONS->nCategories;i++)
        {
          if(OPTIONS->classMethod==2) // Pick mode
          {
            if((this->pClassGivNode[i]) > classProb)
            {
              this->classif = i;
              classProb = (this->pClassGivNode[i]);
            }
          }
          else if(OPTIONS->classMethod==1) // Pick median
          {
            classProb += (this->pClassGivNode[i]);
            if((this->classif)==-1 & classProb>=0.5)
            {
              this->classif = i;
            }
          }
          else
          {
            Rcout << "Something's gone wrong!" << std::endl;
          }
        }
        
        // With classification and counts, we can get impurity and risk
        this->risk = 0;
        this->impurity = 0;
        for (int iCat=0;iCat<OPTIONS->nCategories;iCat++)
        {
          
          this->risk += (this->pClassGivNode[iCat]) * OPTIONS->loss[iCat][classif];
          
          for(int jCat=0;jCat<OPTIONS->nCategories;jCat++)
          {
            this->impurity += OPTIONS->loss[iCat][jCat] * this->pClassGivNode[iCat]*this->pClassGivNode[jCat];
          }
        }
      }
      else if (OPTIONS->mode == 1)
      {
        
        
        // Set up regression model
        
        // Probability of visiting this node
        // We should integrate across some prior distribution
        // to match what's done for regression trees, but for the moment
        // we'll just assume that P(A)=groupSize/totalSize.
        // This *should* be the same as assuming that the
        // sample distribution is exactly representative of the
        // underlying distribution but I haven't proved this
        
        this->probability = groupSize/((double)OPTIONS->model_nObs);
        
        this->mean = 0;
        for(int i=0; i<groupSize;i++)
        {
          this->mean += target_continuous[idList[i]];
        }
        this->mean=(this->mean)/(this->groupSize);
        
        
        // With mean calculated, we can calculate SST
        this->impurity = 0;
        for(int i=0; i<(this->groupSize);i++)
        {
          this->impurity += pow(this->mean - target_continuous[idList[i]],2);
        }
        // Risk is the variance within the node
        if(this->groupSize>1)
        {
          this->risk = this->impurity/((double) this->groupSize-1);
        }
        else
        { // Hopefully not needed but never underestimate the
          // stupidity of users
          this->risk=0;
        }
      }
      else if (OPTIONS->mode == 2)
      {
        // See regression tree comments
        this->probability = groupSize/((double)OPTIONS->model_nObs);
        // Set up distance-based regression model
        
        // For db-regression tree, we skip straight to
        // calculating the sum of squares within each group
        
        // This method pulls from the following:
        
        // Le Meur, N., Vigneau, C., Lefort, M., Lebbah, S., Jais, J. P.,
        // Daugas, E., & Bayat, S. (2019).
        // Categorical state sequence analysis and regression tree
        // to identify determinants of care trajectory in chronic disease:
        // Example of end-stage renal disease.
        // Statistical methods in medical research, 28(6), 1731-1740.
        
        // This paper uses permutation tests for stopping criterion.
        // To encapsulate it in conventional CART methods, we need
        // a risk measure. Fortunately, the paper gives us the variance
        // as
        //
        //    1/(2W^2) sum_i=1^n sum_j=(i+1)^n w_i * w_j * d_(i,j)
        //
        //  Or more plainly:
        //
        //   risk = 1/(2*(groupSize)^2) * impurity
        //
        
        this->impurity = 0;
        for(int i=0; i<(this->groupSize);i++)
        {
          for(int j=i+1; j<(this->groupSize);j++)
          {
            this->impurity = (this->impurity) + target_distance[idList[i]][idList[j]];
          }
        }
        this->risk = this->impurity/(2.0*pow(this->groupSize,2));
        this->impurity;
      }
      
      
      // Now that node calculations have been made, we can find the best split
      
      
      if(!OPTIONS->FOLDING)
      {
        Rcout << this->label << ":  i " << this->impurity << " n " << this->groupSize << " r " << this->risk;
      }
      
      this->split(target_discrete,
                  target_continuous,
                  target_distance,
                  D_list,
                  predictor_distance_combinations);
      
      
      if(OPTIONS->DEBUG)
      {
        Rcout << "Split complete" << std::endl;
      }
      
      // If we've selected permutation-based stopping criterion, we can look at the quality of the split here
      if(OPTIONS->stopMethod==1 & !this->terminal)
      {
        if(OPTIONS->DEBUG)
        {
          Rcout << "Getting in and out groups for permutation test" << std::endl;
        }
        
        // First,  we need to calculate the test statistic using the split we just found.
        // Regardless of the mode we're running in, we'll need the in group and out group
        std::vector<int> inGroup;
        std::vector<int> outGroup;
        
        this->feed( D_list,
                    predictor_distance_combinations,
                    this->groupMembers,
                    inGroup,
                    outGroup
                    );
        
        
        if(OPTIONS->DEBUG)
        {
          Rcout << "Calculating Test Statistic" << std::endl;
        }
        
        double testStatistic = 0;
        if(OPTIONS->mode == 0)
        {
          // For classification tasks, we use a Chisq independence test:
          // X2 = sum_groups( sum_categories ( ( O - E / E )^2 ))
          // Where O is the tallies in each group and E is the expected tally assuming group and outcome are independent
          
          // First, get the tally
          
          std::vector< std::vector<int> > tallyFinal;
  
          std::vector<int> tallyRow(OPTIONS->nCategories);
          
          for(int i=0; i<OPTIONS->nCategories; i++) tallyRow[i] = 0;
          for(int i=0; i<inGroup.size();i++)
          {
            tallyRow[target_discrete[inGroup[i]]]++;
          }
          tallyFinal.push_back(tallyRow);
          
          for(int i=0; i<OPTIONS->nCategories; i++) tallyRow[i] = 0;
          for(int i=0; i<outGroup.size();i++)
          {
            tallyRow[target_discrete[outGroup[i]]]++;
          }
          tallyFinal.push_back(tallyRow);
          
          
          
          if(OPTIONS->DEBUG)
          {
          for(int iGroup=0;iGroup<2;iGroup++)
          {
            for(int iCat=0;iCat<OPTIONS->nCategories;iCat++)
            {
              Rcout << tallyFinal[iGroup][iCat] << " ";
            }
            Rcout << std::endl;
          }
          Rcout << std::endl << std::endl;
          }
          
          
          // Now that the tally is calculated we can get the test statistic
          
          for(int iGroup=0;iGroup<2;iGroup++)
          {
            for(int iCat=0;iCat<OPTIONS->nCategories;iCat++)
            {
              if(tallyFinal[0][iCat] + tallyFinal[1][iCat] > 0)
              {
                
                double E = (
                  tallyFinal[0][iCat] + tallyFinal[1][iCat]
                )*(
                    iGroup==0 ? inGroup.size() : outGroup.size()
                )/((double) this->groupSize);
                
              
                if(OPTIONS->DEBUG)
                {
                 Rcout << iGroup << " " << iCat << " " << E << std::endl; 
                }
              
                testStatistic += (tallyFinal[iGroup][iCat] - E) * (tallyFinal[iGroup][iCat] - E)/ E;
              }
            }
          }
          
          
          
        }
        else if(OPTIONS->mode == 1)
        {
          Rcout << "Regression tree permutation not implemented yet";
        }
        else if(OPTIONS->mode == 2)
        {
          // this->impurity gives us SST. We need to calculate SSW and SSB in order to get the test statistic
          double SSB = 0;
          
          for(int i=0; i<inGroup.size();i++)
          {
            for(int j=0; j<outGroup.size();j++)
            {
              SSB += target_distance[inGroup[i]][outGroup[j]];
            }
          }
          
          // SST = SSW + SSB
          double SSW = this->impurity - SSB;

          Rcout << " SSB " << SSB << " SSW "  << SSW <<std::endl;
          
          // F =  (SSB/(m-1)) / (SSW/(W-m))
          // where m is the number of groups (in our case, 2)
          // and W is summed weights (in our case, sample size)
          
          testStatistic = SSB / SSW * (inGroup.size() + outGroup.size() -1);
        }
        else
        {
          Rcout << "This shouldn't be seen" << std::endl;
        }
        
        
        if(OPTIONS->DEBUG)
        {
          Rcout << "Test Statistic is " << testStatistic << std::endl;
          Rcout << "Starting permutation test" << std::endl;
        }
        
        double pVal = permutationTest(target_discrete,
                                      target_continuous,
                                      target_distance,
                                      this->groupMembers,
                                      inGroup.size(),
                                      testStatistic);
        
        if(OPTIONS->DEBUG)
        {
          Rcout << "Done, P-value is " << pVal << std::endl;
        }
        
        
        this->pVal=pVal;
        
        // // If p-value is above significance level,
        // // then don't make the split
        // if(pVal>OPTIONS->signifLevel)
        // {
        //   best.bubbleValid=0;
        // }
      }
      else
      {
        this->pVal=-1;
      }
      
      
      if(!OPTIONS->FOLDING)
      {
        Rcout << std::endl;
      }
      
    }
    
    ~Node()
    {
      
    }
    
    
    // Feeder function
    // Take a list of records and feed them through the node
    void feed(std::vector< std::vector<std::vector<double> > > &D_list,
              std::vector< std::vector<int> > &predictor_distance_combinations,
              std::vector<int> &input,
              std::vector<int> &inGroup,
              std::vector<int> &outGroup)
    {
      int inputSize = input.size();
      
      for(int i=0;i<inputSize;i++)
      {
        // SplitMethod is *always* 2 right now
        // but check is in here for future extensions
        if(this->splitMethod==2)
        {
          double thisDistance =0;
          for(int j=0;j<predictor_distance_combinations[this->bestSplitIndex].size();j++)
          {
            int thisFeature = predictor_distance_combinations[this->bestSplitIndex][j];
            thisDistance = thisDistance + std::pow(
              this->weights[j] * 
                D_list[thisFeature][this->center][input[i]],
                OPTIONS->combinationPower);
          }
          thisDistance = std::pow(thisDistance,1.0/OPTIONS->combinationPower);
          
          if( thisDistance < this->radius)
          {
            inGroup.push_back(input[i]);
          }
          else
          {
            outGroup.push_back(input[i]);
          }
        }
      } 
    }
                
    
    
    
    
    // Splitter function
    // Note that we don't actually process any splits here.
    // All this function will do is calculate the instructions
    // for how a split should be calculated, which is then read
    // by the constructor for the overall tree,
    // using a BFS approach.
    
    void split(std::vector<int> &target_discrete,
          std::vector<double> &target_continuous,
          std::vector< std::vector<double> > &target_distance,
          std::vector< std::vector<std::vector<double> > > &D_list,
          std::vector< std::vector<int> > &predictor_distance_combinations
    )
    {
      
      if(OPTIONS->trace>2)
      {
        Rcout << "| Split:";
      }
      
      
      // Find the best split
      double bestImpurity = std::numeric_limits<double>::infinity(); // Track best impurity
      this->bestSplitIndex=-1;
      
      // Split method is always 2 right now
      this->splitMethod = 2;
      
      // Default to no split unless we find a valid one
      this->terminal = true;
      this->splitName = "";
      this->center=-1;
      this->radius = -1;
      
      // Should we attempt a split?
      if((this->groupSize>=OPTIONS->minSplit) &      // Are the enough to attempt a split
         (this->groupSize>=2*OPTIONS->minBucket) &   // Can a split of sufficient size even exist
         (this->depth < OPTIONS->maxDepth) &        // Are we not too deep
         (this->impurity > 0)                      // Is there any point to testing
      )
      {

        // Not implemented yet, but:
        // We treat every class of split type separately.
        // I.e. we get the best factor split, the best numeric
        // split and the best distance split,
        // then we choose the best of the best out of these.
        
        // For every distance measure
        for(int i=0;i<OPTIONS->nDistanceCombinations;i++)
        {
          
          if(OPTIONS->trace>2)
          {
            Rcout << std::endl << "   " << i << ": ";
          }
          
          ////////////////////////////////////////////////////////////////////////////////////
          //
          //                                SOMA START
          //
          //
          //      Search weights in hyperspherical coordinates.
          //
          //      https://en.wikipedia.org/wiki/N-sphere#Spherical_coordinates
          //
          //      w1 = sin(t1)
          //      w2 = sin(t1)cos(t2)
          //      w3 = sin(t1)sin(t2)cos(t3)
          //         ...
          //      w(n-1) = sin(t1)sin(t2)sin(t3)...cos(t(n-1))
          //      wn = sin(t1)sin(t2)sin(t3)...sin(t(n-1))
          //
          //      Total number of angles is equal to number of components - 1
          // 
          //      As sum_i((w_i)^2) = 1, set weights as pow(w_i,2) for each evaluation
          //
          ////////////////////////////////////////////////////////////////////////////////////
          
          int dims = predictor_distance_combinations[i].size();
          
          // Set some constraints
          std::vector<double> parLb;
          std::vector<double> parUb;
          for(int i=0; i<(dims-1); i++)
          {
            // Don't actually test from 0 - 1 in weights
            // Combining weights in this area is ridiculous,
            // just use the simpler weight instead.
            // Adjusting the angles in by 0.1 radians
            // means the actual boundaries for sin(t)^2 and cos(t)^2 are from
            // 0.009966711 to 0.9900333 which
            // is about as close to 0.01-0.99 as we can get without
            // being overly pedantic about the size of this adjustment
            
            parLb.push_back(0 + 0.1);
            parUb.push_back(1.57079632679 - 0.1); // pi/2
          }

          Soma soma(SOMA_OPTIONS->pathLength,
                    SOMA_OPTIONS->step,
                    SOMA_OPTIONS->PRT,
                    SOMA_OPTIONS->popSize,
                    SOMA_OPTIONS->nDirections,
                    SOMA_OPTIONS->migrations,
                    SOMA_OPTIONS->minDiv,
                    dims-1,
                    3,
                    parLb,
                    parUb 
          );
        
          
          
          //With Soma optimiser initialised, find optimal weights
          while(soma.exitCode <0)
          {
            
            if(OPTIONS->DEBUG & soma.nTest>1)
            {
              Rcout << std::endl << "======================" <<std::endl;
              Rcout << "SOMA migration step " << soma.migrationCount << std::endl;
              Rcout << "======================" <<std::endl;
              
            }
          
            // Set up next round of parameters to check
            soma.move();
            
          
            if(OPTIONS->DEBUG & soma.nTest>1)
            {
              Rcout << "There are " << soma.nTest << " weight combinations to test" << std::endl;
            }
          
            for(int popi=0;popi < soma.nTest; popi++)
            {
           
              
              if(OPTIONS->DEBUG & soma.nTest>1)
              {
                Rcout << std::endl << "Testing weights: ";
              }
              
              std::vector<double> test_weights;  
              // Awww yiss hyperspherical coordinates time
              if(dims > 1)
              {
                for(int j=0;j<dims;j++)
                {
                  double thisWeight = 1;
                  for(int k=0; k<j; k++)
                  {
                    thisWeight = thisWeight * sin(soma.pars[popi][k]);
                  }
                  
                  if(j<(dims-1))
                  {
                    thisWeight = thisWeight * cos(soma.pars[popi][j]);
                  }
                  
                  // Square weights should sum to 1
                  thisWeight = thisWeight * thisWeight;
                  
                  if(OPTIONS->DEBUG & soma.nTest>1)
                  {
                    Rcout << thisWeight << " ";
                  }
                  
                  test_weights.push_back(thisWeight);
                }
                
              }
              else
              {
                // If only one feature to optimise then we don't need to actually run SOMA
                // Set the feature weight to 1
                test_weights.push_back(1);
              }
              
              ///////////////////////////////////////////////////////////////////////////////////
              //
              //                         Begin bubble blowing
              //
              ///////////////////////////////////////////////////////////////////////////////////
              
              
              // We don't need to test bubbles centered on every record.
              // If the distance between records A and B is 0 then
              // by the identity of indiscernibles the results will be equivalent.
              //
              // To this end, we can just use a unique list of centers instead of the full thing
              
              std::vector<int> groupMembers_unique;
              
              IntegerVector isUnique(this->groupSize);
              UniqueFinder uniqueObs(this->groupMembers,
                                     D_list,
                                     predictor_distance_combinations[i],
                                     OPTIONS->distance_acceptable_error,
                                     isUnique);
              parallelFor(0, this->groupSize, uniqueObs);
               
              for(int j=0;j<(this->groupSize);j++)
              {
                if(isUnique[j]==1)
                {
                groupMembers_unique.push_back(this->groupMembers[j]);
                }
              }
              
              int nUnique = groupMembers_unique.size();
              
              if(nUnique==0)
              {
                // This should never occur and will be the result of a bug
                Rcout << "WARNING: No unique records found" << std::endl; 
              }
       
              // Once we have a list of unique records, we can look through them to
              // work out what the threshold for distinctiveness is
              
              
              
              double tolerance = 0;
              if(dims > 1)
              {
                ThresholdFinder threshold(groupMembers_unique,
                                          D_list,
                                          predictor_distance_combinations[i],
                                          test_weights,
                                          OPTIONS->combinationPower,
                                          OPTIONS->distance_acceptable_error);
                parallelReduce(0, nUnique, threshold);
                
                tolerance = threshold.minVal/2;
              }
              else
              {
                tolerance = OPTIONS->distance_acceptable_error[predictor_distance_combinations[i][0]];
              }
              
              int nTest = 0;
              
              // Now that we have the unique records,
              // we select the centers we're actually going to
              // evaluate.
              
              // Get the number of centers that we should test
              // as a proportion of the number of unique records.
              // We should always test at least one centre though
              // so just check the most outliery if pCenters <=0.                                                  
              int nCenters = ceil(OPTIONS->pCenters * nUnique);
              if(nCenters<1)
              {
                nCenters = 1;
              }
              std::vector<int> groupMembers_toTest;
              
              // If the number of centers to check is less than literally all
              // possible checks, then we prioritise the "outer" layer of results
              
              if(nCenters < nUnique)
              {
                
                // This algorithm as it is prioritises centering on "extreme" points.
                // This provides a harder solution to the disjointedness problem by preventing
                // splits based in the center of the data.
                // Including variance as a criteria for the next step means that it also captures
                // the center mass of variables as well
                
                
                // Step 1: Select initial record
                
                // How outliery is each potential center? 
                NumericVector sumD(nUnique);
                RowSummer sumDistances(groupMembers_unique,
                                       groupMembers_unique,
                                       D_list,
                                       predictor_distance_combinations[i],
                                                                      test_weights,
                                                                      sumD
                );
                parallelFor(0, nUnique, sumDistances);
                
                
                // Which is the outlieriest?
                MaxFinder maxD(sumD);
                parallelReduce(0, nUnique, maxD);
                
                // Add the chosen record to the list to test and remove it from unique list
                groupMembers_toTest.push_back(groupMembers_unique[maxD.currentLocation]);
                groupMembers_unique.erase(groupMembers_unique.begin()+maxD.currentLocation);
                nUnique--;
                nTest++;
                
                
                // Step 2+:
                
                // Sum distances to all currently selected records and pick
                // the one which is farthest away from all existing records
                while((nUnique>1) & (groupMembers_toTest.size()<nCenters))
                {
                  // Get sum of distances from current selected observations
                  
                  NumericVector sumD(nUnique);
                  RowSummer sumDistances(groupMembers_unique,
                                         groupMembers_toTest,
                                         D_list,
                                         predictor_distance_combinations[i],
                                                                        test_weights, 
                                                                        sumD
                  );
                  parallelFor(0, nUnique, sumDistances);
                  
                  
                  MaxFinder maxD(sumD);
                  parallelReduce(0, nUnique, maxD);
                  
                  groupMembers_toTest.push_back(groupMembers_unique[maxD.currentLocation]);
                  groupMembers_unique.erase(groupMembers_unique.begin()+maxD.currentLocation);
                  nUnique--;
                  nTest++;
                  
                }
                
              }
              else
              {
                // If we'd select everyone anyway skip the algorithm and just copy everything across
                for(int j=0;j<nUnique;j++)
                {
                  groupMembers_toTest.push_back(groupMembers_unique[j]);
                }
                nTest = nUnique;
              }
              
              //Find optimal split for this node
              Bubble bubble; // toil and touble
              bubble.blow(target_discrete,
                          target_continuous,
                          target_distance,
                          D_list,
                          predictor_distance_combinations[i],
                          test_weights, 
                          tolerance,
                          groupMembers,
                          groupMembers_toTest
              );
              
              
              ///////////////////////////////////////////////////////////////////////////////////
              //
              //                         End bubble blowing
              //
              ///////////////////////////////////////////////////////////////////////////////////
              
              // The best possible bubble has been blown using these weights.
              // Pass this result back to SOMA and evaluate the next set of weights
              
              soma.utility[popi] = bubble.bubbleValid>0.5? bubble.impurity : 99999; // Mark utility as arbitrarily large if invalid
              soma.other_values[popi][0] = bubble.center+0.1; // floor(int val + epsilon) to get center index later
              soma.other_values[popi][1] = bubble.radius; // floor(int val + epsilon) to get center index later
              soma.other_values[popi][2] = bubble.bubbleValid + 0.1;
              
            }
            
            if(OPTIONS->DEBUG & soma.nTest>1)
            {
              Rcout << std::endl << "All points calculated" << std::endl;
              
              if(dims>1)
              {
                
                for(int popi=0;popi < soma.nTest; popi++)
                {
                
                  Rcout << soma.popId[popi] << " | (";
                  
                  for(int j=0;j<dims-1;j++)
                  {
                    Rcout << soma.pars[popi][j] << " ";
                  }
                  Rcout << ") ";
                  
                  for(int j=0;j<dims;j++)
                  {
                    double thisWeight = 1;
                    for(int k=0; k<j; k++)
                    {
                      thisWeight = thisWeight * sin(soma.pars[popi][k]);
                    }
                    
                    if(j<(dims-1))
                    {
                      thisWeight = thisWeight * cos(soma.pars[popi][j]);
                    }
                    
                    // Square weights should sum to 1
                    thisWeight = thisWeight * thisWeight;
                    
                    
                    
                    Rcout << thisWeight << " ";
                  }
                  Rcout << "| ";
                  
                  Rcout << std::floor(soma.other_values[popi][0]) << " | "; 
                  Rcout << soma.other_values[popi][1] << " | "; 
                  Rcout << std::floor(soma.other_values[popi][2]) << " | "; 
                  Rcout << soma.utility[popi] << std::endl; 
                }
              }
              else
              {
                Rcout << "N/A | ";
                Rcout << "1 | ";
                Rcout << std::floor(soma.other_values[0][0]) << " | "; 
                Rcout << soma.other_values[0][1] << " | "; 
                Rcout << std::floor(soma.other_values[0][2]) << " | "; 
                Rcout << soma.utility[0] << std::endl; 
              }
            }
            
            
            
            // All current weight configurations have been evaluated, time to update agents and start again
            soma.update();
            
            
            if(OPTIONS->DEBUG & soma.nTest>1)
            {
              Rcout << std::endl << "Agents moved!" << std::endl;
              Rcout << "SOMA status is " << soma.exitCode << std::endl;
              
              if(dims>1)
              {
                
                for(int popi=0;popi < soma.popSize; popi++)
                {
                  
                  Rcout << popi << " | ";
                  
                  for(int j=0;j<dims;j++)
                  {
                    double thisWeight = 1;
                    for(int k=0; k<j; k++)
                    {
                      thisWeight = thisWeight * sin(soma.pars_current[popi][k]);
                    }
                    
                    if(j<(dims-1))
                    {
                      thisWeight = thisWeight * cos(soma.pars_current[popi][j]);
                    }
                    
                    // Square weights should sum to 1
                    thisWeight = thisWeight * thisWeight;
                    
                    Rcout << thisWeight << " ";
                  }
                  Rcout << "| ";
                  
                  Rcout << std::floor(soma.other_values_current[popi][0]) << " | "; 
                  Rcout << soma.other_values_current[popi][1] << " | "; 
                  Rcout << std::floor(soma.other_values_current[popi][2]) << " | "; 
                  Rcout << soma.utility_current[popi] << std::endl; 
                }
              }
              else
              {
                Rcout << "N/A | ";
                Rcout << "1 | ";
                Rcout << std::floor(soma.other_values[0][0]) << " | "; 
                Rcout << soma.other_values[0][1] << " | "; 
                Rcout << std::floor(soma.other_values[0][2]) << " | "; 
                Rcout << soma.utility[0] << std::endl; 
              }
            }
            
          
            
          }
          
          
         
          // SOMA should have terminated, we can extract the best split now
          
          std::vector<double> selectedWeights;  
          // Awww yiss hyperspherical coordinates time
          if(dims > 1)
          {
            for(int j=0;j<dims;j++)
            {
              double thisWeight = 1;
              for(int k=0; k<j; k++)
              {
                thisWeight = thisWeight * sin(soma.pars_current[soma.rank[0]][k]);
              }
              
              if(j<(dims-1))
              {
                thisWeight = thisWeight * cos(soma.pars_current[soma.rank[0]][j]);
              }
              
              // Square weights should sum to 1
              thisWeight = thisWeight * thisWeight;
             
              selectedWeights.push_back(thisWeight);
              
            }
          }
          else
          {
            // If only one feature to optimise then we don't need to actually run SOMA
            // Set the feature weight to 1
            selectedWeights.push_back(1);
          }
          
          if(OPTIONS->DEBUG & soma.nTest>1)
          {
            Rcout << "SOMA is complete!" << std::endl << std::endl;
            
            Rcout << "Leader is " << soma.rank[0] << ": U=" << soma.utility_current[soma.rank[0]] << std::endl;
            Rcout << "Pars are ";
            for(int iPar = 0; iPar < dims-1; iPar++)
            {
              Rcout << soma.pars_current[soma.rank[0]][iPar] << " ";
            }
            Rcout << std::endl;
            
            
            
          }
          
          
          int selectedCenter = std::floor(soma.other_values_current[soma.rank[0]][0]);
          double selectedRadius = soma.other_values_current[soma.rank[0]][1];
          int selectedValid = std::floor(soma.other_values_current[soma.rank[0]][2]);
          double selectedImpurity = soma.utility_current[soma.rank[0]];
          
          
          ////////////////////////////////////////////////////////////////////////////////////
          //
          //                                SOMA END
          //
          ////////////////////////////////////////////////////////////////////////////////////
          
          
          if(OPTIONS->trace>2)
          {
            if(dims>1)
            {
              Rcout << "c " << selectedCenter << "/ w " <<  "( " << selectedWeights[0];
              for(int iDims = 1; iDims < dims; iDims++)
              {
                Rcout << ", " << selectedWeights[iDims];
              }
              Rcout <<")" << "/ r " << selectedRadius << ":" << selectedImpurity;
            }
            else
            {
              Rcout << "c " << selectedCenter << "/ r " << selectedRadius << ":" << selectedImpurity;
            }          
            
            
            
          }
          
          
          
          goodnessOfFit.push_back(this->impurity - selectedImpurity);
          
          // Is this bubble the best we've seen thus far?
          bool update=false;
          if(selectedValid>0)
          {
            if(selectedImpurity<bestImpurity)
            {
              update=true;
            }
          }
          
          
          
          if(update)
          {
            // This split option is better than the last candidate!
            
            // Stick a star next to it in the trace if it caused an update
            if(OPTIONS->trace>2)
            {
              Rcout << "*";
            }
            
            
            bestImpurity = selectedImpurity;
              
            this->terminal = false;
            this->bestSplitIndex=i;
            
            if(dims==1)
            {
              this->splitName = OPTIONS->distanceNames[i];
            }

            //  this->splitName = OPTIONS->distanceNames[this->bestSplitIndex];
            this->center=selectedCenter;
            this->radius = selectedRadius;
            
            
            this->weights.clear();
            for(int j=0;j<dims;j++)
            {
              this->weights.push_back(selectedWeights[j]);
            }
            
          }
        }
        
        
        if(OPTIONS->DEBUG)
        {
          Rcout << "Distance splits checked" << std::endl;
        }
        
      
        if(OPTIONS->DEBUG)
        {
          Rcout << "Split is complete, updating results" << std::endl;
        }
        
        
        if(this->terminal & OPTIONS->trace>2)
        {
          Rcout << " No valid split found.";
        }
       
        
      }
      else // No point in checking split
      {
        
        if(OPTIONS->trace>2)
        {
          Rcout << " Not allowed.";
        }
        // No split was attempted, so fill in instructions with dummy information. 
        
        this->splitMethod = 2;
        
        this->terminal = true;
        this->splitName = "";
        
        this->center=-1;
        this->radius = -1;
        
      }
    } // end split() definition
    
    

  }; // end struct Node
  
  // Store nodes in a std::vector
  std::vector<Node*> nodeList;
  
  
  // Store tree properties at each prune
  std::vector<int> nTerms;
  std::vector<double> risk;
  std::vector<double> alpha;
  
  
  // The meat of this approach is here
  RcknnTree_bfs(std::vector<int> &target_discrete,
                std::vector<double> &target_continuous,
                std::vector< std::vector<double> > &target_distance,
                std::vector< std::vector< std::vector<double> > > &predictor_distance,
                std::vector< std::vector<int> > &predictor_distance_combinations,
                std::vector<int> idList)
  {
    
    //////////////////////////////////////////////////////////////////
    //
    //
    //                    TREE CONSTRUCTION 
    //
    //
    ///////////////////////////////////////////////////////////////////
  
    // Construct level 0 node
    
    int label=0;

    Node *rootNode = new Node(target_discrete,
                              target_continuous,
                              target_distance,
                              predictor_distance,
                              predictor_distance_combinations,
                              idList,
                              label,
                              0
    );
    
    // Start splitting from root node
    
    std::queue<Node*> nodeQueue;
    nodeQueue.push(rootNode);
    
    while(!nodeQueue.empty())
    {
      // Get the next node in the queue
      Node *thisNode = nodeQueue.front();
      
      if(!thisNode->terminal)
      {
        
        // If we need to make any splits, generate the new nodes
        // and add them to the end of the queue  
       
        // First, get in/out groups and pass to the new nodes
        std::vector<int> ingroupMembers;
        std::vector<int> outgroupMembers;
        thisNode->feed( predictor_distance,
                        predictor_distance_combinations,
                        thisNode->groupMembers,
                        ingroupMembers,
                        outgroupMembers);

        if(OPTIONS->DEBUG)
        {
          Rcout << "Creating In Child" << std::endl;
        }
        
        // Now generate the in group as well as the out group.
        // Add them both to the end of the queue.
        Node *newNode;
        
        label++;
        newNode = new Node(target_discrete,
                           target_continuous,
                           target_distance,
                           predictor_distance,
                           predictor_distance_combinations,
                           ingroupMembers,
                           label,
                           thisNode->depth+1
        );
        nodeQueue.push(newNode);
        thisNode->inChild = label;
          
          
          
        if(OPTIONS->DEBUG)
        {
          Rcout << "Creating Out Child" << std::endl;
        }
          
          
        label++;
        newNode = new Node(target_discrete,
                           target_continuous,
                           target_distance,
                           predictor_distance,
                           predictor_distance_combinations,
                           outgroupMembers,
                           label,
                           thisNode->depth+1
        );
        nodeQueue.push(newNode);
        thisNode->outChild = label;
        
        if(OPTIONS->DEBUG)
        {
          Rcout << "Children created and added to queue" << std::endl;
        }
        
        
      }
      
      // Finally, add this node to the master list and delete it from the queue
      nodeList.push_back(thisNode);
      nodeQueue.pop();
    }
    
    //////////////////////////////////////////////////////////////////////////
    //
    //
    //            TREE CLEANUP AND LAST-MINUTE CALCULATIONS
    //
    //
    /////////////////////////////////////////////////////////////////////////
    if(OPTIONS->DEBUG)
    {
      printTree(0);  
    }
    
    // Prune the tree
    prune();
    
    
  } // End constructor
  
  ~RcknnTree_bfs()
  { // Delete all of the nodes when the tree is destroyed
    if(OPTIONS->DEBUG)
    {
      Rcout << "Destroying the tree!" << std::endl;
    }

    // While this *should* be deleted on the stack as far as I'm aware,
    // Rcpp is holding onto values which it really shouldn't bother keeping around.
    while(this->nTerms.size()>0){this->nTerms.pop_back();}
    while(this->risk.size()>0){this->risk.pop_back();}
    while(this->alpha.size()>0){this->alpha.pop_back();}
    
    while(nodeList.size()>0)
    {
      if(OPTIONS->DEBUG)
      {
        Rcout << "Deleting node " << nodeList[0]->label << std::endl;
      }
      delete nodeList[0];
      nodeList.erase(nodeList.begin());
    }
  }
  
  ///////////////////////////////
  //
  // Recursive functions
  //
  // Even though we've constructed
  // the tree in a BFS manner,
  // it's navigated in a DFS manner
  
  // Prune tree as set out in Chapter 3 of
  // Breiman et al (1984).
  
  // Critical alpha levels start with the largetst tree
  // then iteratively remove one record at a time. 
  // This function doesn't perform pruning, it
  // just sets up critical values for alpha, at each subtree
  // and returns the label for the weakest identified link.
  // Pruning occurs in getBetas()
  
  int getCriticalAlpha(int label, int pruneStep)
  {
    // Get the weakest link in the tree, corresponding to
    // the minimum alpha value
    
    Node *thisNode = this->nodeList[label];
    
    // If this node is terminal or was pruned off
    if(thisNode->terminal | ( thisNode->pruneAtStep >=0 & thisNode->pruneAtStep < pruneStep) )
    {
      thisNode->nTerminal = 1;
      thisNode->subTreeRisk = thisNode->probability * thisNode->risk;
      thisNode->criticalAlpha = -1;
      
    
      return label;
    }
    else
    {
      // We can't calculate anything until we know what
      // the children of this node look like
      int inLabel = getCriticalAlpha(thisNode->inChild,pruneStep);
      int outLabel = getCriticalAlpha(thisNode->outChild,pruneStep);
      
      // Before we calculate the critical value of alpha for this node,
      // work out which of the child nodes was the weaker link
      Node *min_in = this->nodeList[inLabel];
      Node *min_out = this->nodeList[outLabel];
       
      int weakestLink = inLabel;
      double minAlpha = min_in->criticalAlpha;
      
      if( (min_out->criticalAlpha >-0.1 & min_out->criticalAlpha < minAlpha) |  minAlpha < -0.1 )
      {
        weakestLink = outLabel;
        minAlpha = min_out->criticalAlpha;
      }
      
      // Calculate subtree risk and critical alpha for this node
      
      Node *inChild = this->nodeList[thisNode->inChild];
      Node *outChild = this->nodeList[thisNode->outChild];
      
      thisNode->nTerminal = inChild->nTerminal + outChild->nTerminal;
      thisNode->subTreeRisk = inChild->subTreeRisk + outChild->subTreeRisk;
      
      
      // We're getting floating point errors here. Just look at difference greater than epsilon
      // for if it's a real split
      if((thisNode->probability)*(thisNode->risk) - thisNode->subTreeRisk > 1e-6)
      {

        thisNode->criticalAlpha = (
                                    (thisNode->probability)*(thisNode->risk)-
                                    (thisNode->subTreeRisk)
                                   )/(thisNode->nTerminal-1);
      }
      else
      {
        thisNode->criticalAlpha = 0;
      }
      
      // Correct floating point error in critical Alpha value if neccesary
      if(thisNode->criticalAlpha < 0)
      {
        thisNode->criticalAlpha = 0;
      }
      
      
      // If this node's critical alpha is smaller than the smallest subtree,
      // mark it as the smallest alpha value
      if( ( thisNode->criticalAlpha > -0.1 & thisNode->criticalAlpha < minAlpha) | minAlpha < -0.1)
      {
        weakestLink = label;
        minAlpha = thisNode->criticalAlpha;
      }

      
      return weakestLink;
    }
  }
  
  
  
  void prune()
  {

    // Hoo boy okay here we go
    //
    // So whenever we calculate criticalAlpha for getting Tn, we calculate the size and risk of T(n-1)
    // As such, the way to populate vectors of alpha values is to lag them. We store alpha
    // at each step to stick into the results vector on the next iteration
    //
    
    int weakestLinkIndex = -1;
    double thisAlpha = -1; // Tmax corresponds to no pruning whatsoever even if it doesn't reduce the risk
    int pruneStep = 0;
    
    bool pointlessSplitFound = false; // See bottom of function for explanation
    
    
    // We evaluate the risk of each pruned tree from the root
    Node* rootNode = this->nodeList[0]; 
  
    do{
    
    
      if(OPTIONS->DEBUG)
      {
        Rcout << "Prune step " << pruneStep << std::endl;
      }
      
      // Calculate subtree risk and subtree size.
      // At the root, this is just tree risk and size.
      
      // Also calculate critical alpha for each node in the tree.
      // This tells us what node is target for the next prune step.
      // Returns index corresponding to the smallest critical alpha value we see.
      weakestLinkIndex = getCriticalAlpha(0, pruneStep);
      
      if(OPTIONS->DEBUG)
      {
        // Print the alphas for tree Ti
        Rcout << "Node:   PruneStep:   Alpha:" << std::endl;
        for(int i=0;i<this->nodeList.size(); i++)
        {
          
          Node* thisNode = this->nodeList[i];
          
          Rcout << thisNode->label << " ";
          Rcout << thisNode->pruneAtStep << " ";
          Rcout << thisNode->criticalAlpha << std::endl;
        }
        
        Rcout << "Weakest link at: " << weakestLinkIndex << std::endl;
      }
      
      
      
      this->alpha.push_back(thisAlpha);
      this->nTerms.push_back(rootNode->nTerminal);
      this->risk.push_back(rootNode->subTreeRisk);
      
      Node* weakestLink  = this->nodeList[weakestLinkIndex];

      thisAlpha = weakestLink->criticalAlpha;
      
      

      if(thisAlpha == 0)
      {
        pointlessSplitFound = true;
      }

      // Loop over all nodes and prune all that have an alpha which match this one.
      // This should be unneseary if Proposition 3.6 and 3.7 of Friedman (1980) holds
      // but if at the first step alpha = 0 then we'll need to trim a bunch of nodes
      // did not result in any reduced complexity.
      

      for(int i=0;i<this->nodeList.size();i++)
      {
        
        Node* thisNode = this->nodeList[i];
        
        if(thisNode->criticalAlpha == thisAlpha & thisNode->pruneAtStep == -1)
        {
          thisNode->pruneAtStep = pruneStep;
        }
      }
      
            
      pruneStep++;
    } while (weakestLinkIndex!=0); // Until we've pruned at the root node
    
    // Final prune is at the root node, add it to the alpha/beta/nTerm list
    
    this->alpha.push_back(thisAlpha);
    this->nTerms.push_back(1);            // There is one node 
    this->risk.push_back(rootNode->risk); // Risk is just the risk of the root node
    
    // Before we end this function, we need to hack around a problem with how to
    // interpret these numbers based on if there's pointless splits (Tmax > T0) or
    // if there's no pointless splits Tmax = T0.
    // Where Tmax>T0, pruneAtStep==0 corresponds to the initial prune described on page
    // 68 of Friedman (1980) where all pointless splits are removed.
    // Where Tmax == T0, pruneAtStep==1 corresponds to the first real prune step
    
    // To correct this (i.e. pruneAtStep==0 should return T0, pruneAtStep==1 should return T1, etc)
    // we need to add 1 to all values of pruneAtStep where no futile split was found.
    
    if(!pointlessSplitFound)
    {
      
      for(int i=0; i<this->nodeList.size();i++)
      {
        Node* thisNode = this->nodeList[i];
        thisNode->pruneAtStep++;
      }
      
    }
  }
  
  
  
  // Print the current state of the tree
  void printTree(int label)
  {
    Node *thisNode = this->nodeList[label];
    
    Rcout << std::endl;
    for(int i=0; i<thisNode->depth;i++)
    {
      Rcout << " ";
    }
    
    Rcout << thisNode-> label << ": " <<  thisNode->probability << " ";
    Rcout << thisNode->risk;
    
    if(!thisNode->terminal)
    {
      printTree(thisNode->inChild);
      printTree(thisNode->outChild);
    }
    
    if(label==0)
    {
      Rcout << std::endl;
    }
  }
  
  ///////////////////////////////////////////////////////////
  //
  //
  //          Clustering, classification, etc
  //
  //
  
  // All functions in this section feature a pass-by-reference output vector.
  // This is potentially sparse and should be of length OPTIONS->nObs.
  //
  // This allows us to construct a full vector by running classTree on
  // different CV models using the same output vector.
  
  // Get the terminal node label corresponding to each record in a list
  void classTree(std::vector< std::vector< std::vector<double> > > &predictor_distance,
                 std::vector< std::vector< int > > &combinations,
                 int label, 
                 std::vector<int> predList,
                 int pruneIndex, // Do we want T0, T1, ..., or Tn?
                 std::vector<int> & output
                 )
  {
    Node* thisNode = this->nodeList[label];
    
    int nPreds = predList.size();
    if( (thisNode->terminal) | ( thisNode->pruneAtStep>=0 & (thisNode->pruneAtStep <= pruneIndex) )  )
    {
      // Stopping criteria are met, get the label that these elements belong to
      
      for(int i=0;i<nPreds;i++)
      {
        output[predList[i]] = thisNode->label; 
      }
      
    }
    else
    {
      std::vector<int> inGroup;
      std::vector<int> outGroup;
      
      int nFeatures = combinations[thisNode->bestSplitIndex].size();
      
      for(int i=0;i<nPreds;i++)
      {
        if(thisNode->splitMethod==2)
        {
          double thisDistance = 0;
          for(int j=0; j<nFeatures; j++)
          {
            double addD = thisNode->weights[j] * predictor_distance[combinations[thisNode->bestSplitIndex][j]][thisNode->center][predList[i]];
            thisDistance = thisDistance + std::pow(addD,OPTIONS->combinationPower) ;
          }
          thisDistance = std::pow(thisDistance,1.0/OPTIONS->combinationPower);
          
          if(thisDistance < thisNode->radius)
          {
            inGroup.push_back(predList[i]);
          }
          else
          {
            outGroup.push_back(predList[i]);
          }
        }
      }
      
      classTree(predictor_distance,
                combinations,
                thisNode->inChild, 
                inGroup,
                pruneIndex,
                output);
      
      classTree(predictor_distance,
                combinations,
                thisNode->outChild, 
                outGroup,
                pruneIndex,
                output);
      
    }
  }
  
  // Get prediction for each record at specified beta
  void prediction(std::vector< std::vector< std::vector<double> > > &predictor_distance,
                  std::vector< std::vector< int > > &combinations,
                  int label, 
                  std::vector<int> predList,
                  int pruneIndex,
                  std::vector<double> & output
                  )
  {
    // Get labels corresponding to each patient
    
    int nPreds = predList.size();
    std::vector<int> labels(OPTIONS->nObs,-1);
    this->classTree(predictor_distance,
                    combinations,
                    0,
                    predList,
                    pruneIndex,
                    labels);
    
    // Then look up the corresponding prediction
    
    for(int i=0;i<nPreds;i++)
    {
      Node *thisNode = this->nodeList[labels[predList[i]]];
      
      if(OPTIONS->mode==0)
      {
        output[predList[i]] = thisNode->classif;
      }
      else if(OPTIONS->mode==1)
      {
        output[predList[i]] = thisNode->mean;
      }
      else if(OPTIONS->mode==2)
      {
        // Prediction doesn't make sense for distance-based regression,
        // but we can still get a measure of variability in each record.
        
        // For the sake of consistency this is handled by predictionDetails() below
        
      }
    }
    
  }
  
  
  // For each prediction, get more detailed info
  void predictionDetails(std::vector< std::vector< std::vector<double> > > &predictor_distance,
                         std::vector< std::vector< int > > &combinations,
                  int label, 
                  std::vector<int> predList,
                  int pruneIndex,
                  std::vector< std::vector<double> > &output,
                  std::vector< std::vector<double> > &target_distance
  )
  {
    
    // Get labels corresponding to each patient
    
    int nPreds = predList.size();
    std::vector<int> labels(OPTIONS->nObs,-1);
    this->classTree(predictor_distance,
                    combinations,
                    0,
                    predList,
                    pruneIndex,
                    labels);
    
    
    // Then look up the corresponding prediction
    for(int i=0;i<nPreds;i++)
    {
      Node *thisNode = this->nodeList[labels[predList[i]]];
      
      if(OPTIONS->mode==0)
      {
        // For classification, get the probability of each classification in the node
        for(int j=0;j<OPTIONS->nCategories;j++)
        {
          output[predList[i]][j] = thisNode->pClassGivNode[j];
        }
      }
      else if(OPTIONS->mode==1)
      {
        // For regression, get the variance in the node (this is the risk)
        output[predList[i]][0] = thisNode->risk;
      }
      else if(OPTIONS->mode==2)
      {
        // Prediction doesn't make sense for distance-based regression,
        // but we can instead return the average squared distance between
        // each record and whatever the group member was.
        // This is roughly equivalent to the variance within the node,
        // but will be on a different scale and is probabily biased
        
        double SS = 0;
        for(int j=0; j<thisNode->groupSize;j++)
        {
          SS = SS + pow(target_distance[predList[i]][thisNode->groupMembers[j]],2);
        }
        output[predList[i]][0] = SS/(thisNode->groupSize);
      }
    }
  }
  
  // Get risk associated with a node, either
  // using the full dataset or as part of training/testing
  // cross-validation
  double predictionRisk(std::vector<int> &target_discrete,
                        std::vector<double> &target_continuous,
                        std::vector< std::vector<double> > &target_distance, 
                        std::vector< std::vector< std::vector<double> > > &predictor_distance,
                        std::vector< std::vector< int > > &combinations,
                        std::vector<int> predList,
                        int pruneIndex
                        )
  {
    int nPreds = predList.size();
    
    // Get terminal node for each record
    std::vector<int> labels(OPTIONS->nObs,-1);
    
    this->classTree(predictor_distance,
                    combinations,
                    0,
                    predList,
                    pruneIndex,
                    labels);
    
    // With that done it's fairly easy to get the risk
    double risk = 0;
    if(OPTIONS->mode==0)
    {
      // For classification trees, risk is sum of misclassification costs
      
      for(int i=0;i<nPreds;i++)
      {
        Node *thisNode = this->nodeList[labels[predList[i]]];
        risk = risk + OPTIONS->loss[target_discrete[predList[i]]][thisNode->classif];
      }
    }
    else if(OPTIONS->mode==1)
    {
      // For regression trees, risk is the sum of square error from the mean
      for(int i=0;i<nPreds;i++)
      {
        Node *thisNode = this->nodeList[labels[predList[i]]];
        risk = risk + pow(target_continuous[predList[i]] - thisNode->mean,2);
      }
    }
    else if(OPTIONS->mode==2)
    {
      // For distance-based regression trees, rather than summing the
      // square error from the mean, we calculate the average squared distance
      // for each record from records in the node.
      
      for(int i=0;i<nPreds;i++)
      {
        Node *thisNode = this->nodeList[labels[predList[i]]];
        
        double nodeRisk = 0;
        for(int j=0; j<thisNode->groupSize;j++)
        {
          nodeRisk = nodeRisk + pow(target_distance[predList[i]][thisNode->groupMembers[j]],2);
        }
        nodeRisk = nodeRisk/(thisNode->groupSize);
        risk = risk + nodeRisk;
      }
    }
    
    return risk;
    
  }
  
  
  List exportTree()
  {
    
    IntegerVector    out_labels;
    IntegerVector    out_parents;
    IntegerVector    out_depths;
    NumericVector    out_probability;
    List             out_groupMembers;
    IntegerVector    out_groupSizes;
    NumericVector    out_risks;
    IntegerVector    out_pruneAtStep;

    NumericVector    out_impurities;
    IntegerVector    out_classifs;
    NumericVector    out_means;
    LogicalVector    out_terminals;
    NumericVector    out_alpha;
    IntegerVector    out_splitMethods;
    IntegerVector    out_bestSplitIndexes;
    CharacterVector  out_splitNames;
    List             out_goodnessOfSplit;
    NumericVector    out_splitUtilities;
    IntegerVector    out_centers;
    NumericVector    out_radii;
    NumericVector    out_pVals;
    IntegerVector    out_inChildren;
    IntegerVector    out_outChildren;
    
    for(int i=0; i<this->nodeList.size(); i++)
    {
      Node* thisNode = this->nodeList[i];
      
      out_labels.push_back(thisNode -> label);
      out_depths.push_back(thisNode -> depth);
      out_probability.push_back(thisNode -> probability);
      
      out_pruneAtStep.push_back(thisNode->pruneAtStep);
        
      NumericVector theseGoodnesses(OPTIONS->nDistanceCombinations);
      if(thisNode->terminal)
      {
        for(int j=0;j<OPTIONS->nDistanceCombinations;j++)
        {
          theseGoodnesses[j] = -1;
        }
      }
      else
      {
        for(int j=0;j<OPTIONS->nDistanceCombinations;j++)
        {
          theseGoodnesses[j] = thisNode->goodnessOfFit[j];
        }
      }
      out_goodnessOfSplit.push_back(theseGoodnesses);
      
      IntegerVector theseMembers(thisNode->groupSize);
      for(int j=0; j< thisNode->groupSize;j++)
      {
        theseMembers[j] = thisNode->groupMembers[j];
      }
      out_groupMembers.push_back(theseMembers);
      
      out_groupSizes.push_back(thisNode -> groupSize);
      out_risks.push_back(thisNode ->risk );
      out_impurities.push_back(thisNode -> impurity );
      out_classifs.push_back(thisNode -> classif );
      out_means.push_back(thisNode ->mean );
      out_terminals.push_back(thisNode ->terminal );
      out_alpha.push_back(thisNode->criticalAlpha);
      out_splitMethods.push_back(thisNode -> splitMethod);
      out_bestSplitIndexes.push_back(thisNode -> bestSplitIndex);
      out_splitNames.push_back(thisNode -> splitName);
      out_splitUtilities.push_back(thisNode -> splitUtility );
      out_centers.push_back(thisNode ->center );
      out_radii.push_back(thisNode ->radius );
      out_pVals.push_back(thisNode ->pVal );
      out_inChildren.push_back(thisNode ->inChild );
      out_outChildren.push_back(thisNode ->outChild );
    }
    
    
    // List::create() only allows up to 20 elements, see https://stackoverflow.com/questions/27371543/how-many-vectors-can-be-added-in-dataframecreate-vec1-vec2
    // We have 21
    List out = List::create(
     
      Named("groupMembers") = out_groupMembers,
      Named("groupSize") = out_groupSizes,
      Named("risk") = out_risks,
      Named("pruneAtStep") = out_pruneAtStep,
      Named("impurity") = out_impurities,
      Named("classif") = out_classifs,
      Named("mean") = out_means,
      Named("terminal") = out_terminals,
      Named("alpha") = out_alpha,
      Named("splitMethod") = out_splitMethods,
      Named("bestSplitIndex") = out_bestSplitIndexes,
      Named("splitName") = out_splitNames,
      Named("splitUtility") = out_splitUtilities,
      Named("goodnessOfFit") = out_goodnessOfSplit,
      Named("center") = out_centers,
      Named("radius") = out_radii,
      Named("pVal") = out_pVals,
      Named("inChild") = out_inChildren,
      Named("outChild") = out_outChildren
    );

    // Mad haxx
    out.push_front(out_probability,"probability");
    out.push_front(out_depths,"depths");
    out.push_front(out_parents,"parent");
    out.push_front(out_labels,"label");

    return out;
    
  }
  
}; // end struct Tree
