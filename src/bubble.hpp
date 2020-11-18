#include "optim.hpp"

/*
 *    Defines bubble blowing/popping parallel workers.
 *
 *    Note: The current method isn't working
 *    as well as it could because it's can't detect outliers
 *    or where the "out" group corresponds to two disjoint groups.
 *
 *    Checking for disjointness when making a cut is far too expensive and
 *    can't feasibly be done.
 *
 *    How to fix this is an open question - right now we're just ignoring
 *    the method under the justification of not letting perfect be
 *    the enemy of the good
 */


// BubbleBlower will, in parallel, find an optimal "bubble" of data
// centered on each record.

struct BubbleBlower : public Worker
{

  // Version 3 of BubbleBlower
  // Pre-sorts bubble sizes in order to re-use information
  // between size checks

  //For parallelised idList[i], check each serial idList[j] and calculate
  //gini index, minimum size, etc.
  //As we iterate across j, keep track of best bubble size

  // Finally, return  bubble center idList[j], radius,
  // gini index and minimum group size if split

  // Inputs:

  // Target variable
  const std::vector<int>  target_discrete;
  const std::vector<double>  target_continuous;
  const std::vector< std::vector<double> > target_distance;

  const std::vector< std::vector< std::vector<double> > > D; // Distance matrix
  
  const std::vector<int> matrixSelection;
  const std::vector<double> weights;
  double combinationP;
  
  const double tolerance; // Half minimum distance in the above matrix

  const std::vector<int> idList; // List of records in this node
  const std::vector<int> idList_unique; // List of unique records in this node

  // Outputs

  RVector<int> validBubble; // If a valid split was found centered on this bubble
                            //
                            // 0 = not valid
                            // 1 = valid
                            

  RVector<double> radius; // Best radius for center i

  // What's the impurity of this split?
  RVector<double> impurity;
  
  const int nObs;
  const int mode;
  const int nCategories;
  const int minBucket;
  const std::vector<double> priors_altered;
  const std::vector<int> classFreq;
  const std::vector< std::vector<double> > loss;
  
  BubbleBlower(const std::vector<int> &target_discrete,
               const std::vector<double> &target_continuous,
               const std::vector< std::vector<double> > &target_distance,
               const std::vector< std::vector< std::vector<double> > > &D,
               
               const std::vector<int> &matrixSelection,
               const std::vector<double> &weights,
               const double combinationP,
               const double tolerance,
               const std::vector<int> &idList,
               const std::vector<int> &idList_unique,
               
               IntegerVector validBubble,
               NumericVector radius,
               NumericVector impurity,

               //relevant OPTIONS
               const int nObs,
               const int mode,
               const int nCategories,
               const int minBucket,
               const std::vector<double> priors_altered,
               const std::vector<int> classFreq,
               const std::vector< std::vector<double> > loss

  )
  :
  target_discrete(target_discrete),
  target_continuous(target_continuous),
  target_distance(target_distance),
  D(D),
  matrixSelection(matrixSelection),
  weights(weights),
  combinationP(combinationP),
  
  tolerance(tolerance),
  idList(idList),
  idList_unique(idList_unique),
  validBubble(validBubble),
  radius(radius),
  impurity(impurity),
  nObs(nObs),
  mode(mode),
  nCategories(nCategories),
  minBucket(minBucket),
  priors_altered(priors_altered),
  classFreq(classFreq),
  loss(loss)
  {}

  void operator()(size_t begin, size_t end)
  {
    
  unsigned int groupSize = idList.size();
  unsigned int uniqueSize = idList_unique.size();
  unsigned int nMatrices = matrixSelection.size();

  // Initialise tallies once for each worker. Reset these
  // for every bubble center.
  //
  // This method is slow and we need every corner cut that we can possibly cut
  
  int base_countIn = 0;
  int base_countOut = groupSize;

  // Track mode-specific criteria
  // Always initialise everything, but don't count actually
  // calculate things unless we need it
  
  
  // Classification: in/out group tallies.
  std::vector<int> base_tallyInClass(nCategories,0);
  std::vector<int> base_tallyOutClass(nCategories,0);
  if(mode==0)
  {
    // OutClass starts with all observations
    for(unsigned int j=0;j<groupSize;j++)
    {
      base_tallyOutClass[target_discrete[idList[j]]]++;
    }
  }

  
  // Regression or DB-Regression: Sum of target variables or sum square distances
  // If Regression: Convert into into mean whenever we need it
  double base_sumIn = 0;
  double base_sumOut = 0;
  double grandMean;
  if(mode==1)
  {
    for(unsigned int j=0;j<groupSize;j++)
    {
      base_sumOut += target_continuous[idList[j]];
    }
    grandMean = base_sumOut/base_countOut;
  }
  else if(mode==2)
  {
    for(unsigned int j=0;j<groupSize;j++)
    {
      for(unsigned int k=0;k<j;k++)
      {
        base_sumOut += target_distance[idList[j]][idList[k]];
      }
    }
  }

  for(unsigned int i=begin; i<end; i++) // for parallelised bubble center i
  {
  
    
  int thisCenter = idList_unique[i];

  // Initialise variables

  // Get some default values.
  int bestValidBubble=-1;
  double bestRadius = -1;
  double bestImpurity = std::numeric_limits<double>::infinity();


  // Restore initial tallies
  int countIn = base_countIn;
  int countOut = base_countOut;

  std::vector<int> tallyInClass(nCategories,0);
  std::vector<int> tallyOutClass(nCategories,0);
  if(mode==0)
  {
    // OutClass starts with all observations
    for(unsigned int j=0;j<nCategories;j++)
    {
      tallyOutClass[j] = base_tallyOutClass[j];
    }
  }
  double sumIn = base_sumIn;
  double sumOut = base_sumOut;

  /////////////////////////////////////////////////////////////////
  //
  // Copy every element of idList into an ordered version
  // and do the same with the individual distances.
  // This lets us combine multipledistances into a weighted average
  // (see: feature-weighted composite k nearest neighbours)
  
  
  std::vector<int> idList_ordered;
  std::vector<double> d_ordered;
  
  for(unsigned int j=0;j<groupSize;j++)
  {
    int id = idList[j];
    
    double d = 0;
    for(int k=0;k<nMatrices;k++)
    {
      double addD = weights[k]*D[matrixSelection[k]][thisCenter][id];
      d=d+std::pow(addD,combinationP);
    }
    d = std::pow(d,1.0/combinationP);

    if(idList_ordered.size()==0) // If list is empty, just add it
    {
      idList_ordered.push_back(id);
      d_ordered.push_back(d);
    }
    else
    {
      unsigned int lower = 0;
      unsigned int upper = idList_ordered.size()-1;

      if( d < d_ordered[lower]+tolerance) // We're below or equal to the current low
      {
        // Insert at the front
        idList_ordered.insert(idList_ordered.begin(),id);
        d_ordered.insert(d_ordered.begin(),d);
      }
      else if ( d > d_ordered[upper]-tolerance)  // We're above or equal to the current upper
      {
        // Insert at the back
        idList_ordered.push_back(id);
        d_ordered.push_back(d);
      }
      else
      {
        // Bisection method madness goes here
        unsigned int position = ceil(0.5*upper);
        while(lower!=upper)
        {
          if((d_ordered[position-1]-tolerance < d) &  (d < d_ordered[position]+tolerance) )
          {    // Are we where we're meant to be?
            lower=position;
            upper=position;
          }
          else if ( (d_ordered[position-1]-tolerance < d) &  (d_ordered[position]-tolerance < d))
          {   // Is our current position too low?
              // If so, shift upward
              lower=position;
              position=ceil(0.5*(lower+upper));
          }
          else if ( (d < d_ordered[position-1]+tolerance) &  (d < d_ordered[position]+tolerance))
          { // Is our current position too high?
            // If so, shift downward
              upper=position;
              position=floor(0.5*(lower+upper));
          }
        } // End while

        idList_ordered.insert(idList_ordered.begin()+position,id);
        d_ordered.insert(d_ordered.begin()+position,d);

      } // End if
    } // End if
  } // End for
 
 /* 
  Rcout << std::endl;
  Rcout << "==============" << std::endl ;
  Rcout << "CENTER" << thisCenter << ":" << std::endl;
  Rcout << "==============" << std::endl ;
  
  for(int j=0;j<groupSize;j++)
  {
    Rcout << idList_ordered[j] << ": ";
    Rcout << d_ordered[j] << " | ";
    Rcout << target_discrete[idList_ordered[j]] << std::endl;
  }
  */
 
  // Test if the above failed
  bool sorted = true;
  for(int j=0; j<(groupSize-1); j++)
  {
    // Test that ordering of idList is correct
    if(d_ordered[j] > d_ordered[j+1]+tolerance)
    {
      Rcout << idList_ordered[j] << " -> " << idList_ordered[j+1] << ": ";
      Rcout << d_ordered[j] << " > " << d_ordered[j+1] << std::endl;
      sorted = false;
    }
    
    // Test that d_ordered contains the correct values 
    // if(abs(D[thisCenter][idList_ordered[j]]-d_ordered[j])>tolerance)
    // {
    //   Rcout << "Mismatch between idList_ordered and d_ordered!" << std::endl;
    //   sorted = false;
    // }
    // 
  }
  if(!sorted)
  {
    Rcout << "Sort failed!";
  }

  
  /////////////////////////////////////////////////////////////////
  //
  // Now that idList is sorted, we can inflate a bubble around it
  //

  for(unsigned int j=0;j<(groupSize-1);j++) // for serial j.
  {                            // Note we have to start from 0 to properly track
                               // the tallies, we can't just skip to OPTIONS->minBucket
                               
    // Blow a larger and larger bubble around i, successively capturing each
    // idList_ordered[j].

    // Keep track of summing variables as we inflate the bubble
    countIn++;
    countOut--;

    if(mode==0)
    {
      tallyInClass[target_discrete[idList_ordered[j]]]++;
      tallyOutClass[target_discrete[idList_ordered[j]]]--;
    }
    else if(mode==1)
    {
      sumIn += target_continuous[idList_ordered[j]];
      sumOut -= target_continuous[idList_ordered[j]];

      // Calculate SST later. It costs us an additional O(n) each time
      // so we'll only do it when the bubble is valid and needs testing
    }
    else if(mode==2)
    {
      // Here's an explanation as to how this bit works because I'm
      // bad at my job and will forget this later.
      //
      // Imagine we're at step j=2, i.e. records 0 and 1 are already in the
      // "in" group. The response matrix looks like this:
      //
      //    (2)(3)
      //  0 1|2|3 4 5
      //    0|1|4 2 3
      //  ---+-+-----
      //     |0|2 1 3
      //       +-----
      //       |0 1 3
      //          0 4
      //            0
      //
      // At the end of step j=2, we have:
      //
      // SST = 39
      // SST_in = 1
      // SST_out = 14
      //
      // At the end of step j=3, we should have:
      //
      // SST=39
      // SST_in = 1+(2+1) = 4
      // SST_out = 14 - (2+1+3) = 8
      //
      // i.e. the trick is to add the upper diagonal portion of column j to
      // SST_in, and subtract the upper diagonal portion of row j from SST_out
      //
      //  Alternatively, if we've added element j to the in group,
      //  we need to add comparisons between it and all 0,1,...,j-1
      //  elements of the in group.
      //
      //  Likewise, if we've removed element j from the out group,
      //  we need to remove comparisons between it and all j+1,j+2,...,n
      //  elements of the out group.
      //
      for(int k=0; k<j;k++)
      {
        sumIn += target_distance[idList_ordered[j]][idList_ordered[k]];
      }
      for(int k=j; k<groupSize; k++)
      {
        sumOut -= target_distance[idList_ordered[k]][idList_ordered[j]];
      }
    }
    
    
    // Only test this radius if:
    int test=0;
    if(j<(groupSize-1)) // It's not the last bubble (never valid), AND
    {
      // The next bubble is bigger than this one (i.e. skip redundant bubble sizes)
      if(d_ordered[j+1] >d_ordered[j] + tolerance)
      {
        // AND one of the following holds:

        // 1) The group sizes are acceptably large (this bubble is valid)
        if( (countIn >= minBucket) & (countOut >= minBucket) )
        {
          test=1;
        }
      }
      
    }

    
    if(test!=0) // We care about this bubble, evaluate it
    {
      // Temprorary variables to store the results of this bubble size
      // Put the radius in the middle of this bubble size and the next
      double thisRadius= 0.5*(d_ordered[j]+d_ordered[j+1]);
      int thisBubbleValid = test;

      double thisImpurity = 0;

      // Actual evaluation methods for testing depend
      // on the type of tree we're constructing

      if(mode==0)
      {
        // Classification Tree
        std::vector< std::vector<int> > tallyFinal; // Final tallies for each group
        tallyFinal.push_back(tallyInClass);
        tallyFinal.push_back(tallyOutClass);

        // Calculate gini impurity for this split.
        // Sum across groups
        for(int iGroup=0;iGroup<2;iGroup++)
        {
          int thisCount = 0;
          double thisP = 0;

          // Get probability of visiting this node
          for(int iCat=0;iCat<nCategories;iCat++)
          {
            thisCount=thisCount+tallyFinal[iGroup][iCat];
            thisP = thisP + (priors_altered[iCat] * tallyFinal[iGroup][iCat] )/(double)classFreq[iCat];
          }

          // Get class probabilities
          std::vector<double> pGivenA;
          for(int iCat=0;iCat<nCategories;iCat++)
          {
            double p_i = (priors_altered[iCat] * tallyFinal[iGroup][iCat] )/(double)classFreq[iCat];
            p_i = p_i/thisP;
            pGivenA.push_back(p_i);
          }

          // Get Gini impurity
          for(int iCat=0;iCat<nCategories;iCat++)
          {
            for(int jCat=0;jCat<nCategories;jCat++)
            {
              double giniContribution = thisP*( loss[iCat][jCat] * pGivenA[iCat] * pGivenA[jCat] );
              thisImpurity=thisImpurity+giniContribution;
            }
          }
        } // End for each group
      }
      else if(mode==1)
      {
        // Regression Tree

        // Calculate means and sum of squares from means
        // for in and out groups.
        // Impurity is the sum of these.

        // See cran/rpart anova.c for details
        
        double meanIn = sumIn/countIn;
        double meanOut = sumOut/countOut;
        
        thisImpurity = countIn*(meanIn  - grandMean)*(meanIn  - grandMean) +
                      countOut*(meanOut - grandMean)*(meanOut - grandMean);
        
        
        // The above measure is for improvement not impurity
        thisImpurity = -thisImpurity;
        
      }
      else if(mode==2)
      {
        // DB-Regression Tree

        // This is calculated in the exact same way as the
        // conventional regression tree, except we've already
        // run calculated the sum of squares within each group
        thisImpurity = sumIn/((double) countIn)+sumOut/((double) countOut);
      }


      // With all that done, is this split better than the one we've seen thus far?
      // The rule here depends on if it's a valid bubble or one that's been flagged
      // as disjoint
      bool update=false;
      if(
          (thisBubbleValid==2) & // Bubble was flagged for disjointness
          (
            (bestValidBubble!=2) | // Current bubble is not flagged, OR
            (
                (thisImpurity<bestImpurity) | // If we improve impurity, OR
                (
                  (thisImpurity==bestImpurity) & // If we're tied for impurity but have a larger radius
                  (thisRadius>=bestRadius)
                )
            )
          )
      )
      {
        // Rcout << "a";
        update = true;
      }
      else if ( (thisBubbleValid==1) & (bestValidBubble!=2) & // Bubble is valid and best bubble is not flagged
                (
                  // If we're tied for impurity but have a larger radius
                  ((thisImpurity==bestImpurity) & (thisRadius>=bestRadius)) | // OR
                  (thisImpurity<bestImpurity) //  If we improve impurity
                )
      )
      {
        update = true;
      }


      if(update)
      {
        // If so, update our results
        bestValidBubble=thisBubbleValid;
        bestRadius=thisRadius;
        bestImpurity=thisImpurity;

        // Rcout << bestValidBubble << " " << bestRadius << " " << bestGini << std::endl;

      }

    } // End if we tested this bubble

  } // End for serial j

  // Finally, store the best results

  // Rcout << "valid " << bestValidBubble << " r: " << bestRadius << " I: " << bestImpurity;
  // Rcout << std::endl;
  
  validBubble[i] = bestValidBubble;
  radius[i] = bestRadius;
  impurity[i] = bestImpurity;

  } // End parallelised i
  } // End operator()

}; // End struct BubbleBlower


// BubblePopper will, in parallel, find the best bubble returned by
// BubbleBlower

struct BubblePopper : public Worker
{
  // Select rercord i for final bubble

  // Get inputs from BubbleBlower
  const RVector<int> valid;
  const RVector<double> radius;
  const RVector<double> impurity;

  // Track best results seen thus far
  int    finalValid;
  double finalRadius;
  double finalImpurity;
  double finalUtility;

  // index for selected bubble
  int finalBubble;

  BubblePopper(
    const IntegerVector valid,
    const NumericVector radius,
    const NumericVector impurity
  )
  :
  valid(valid),
  radius(radius),
  impurity(impurity),
  finalValid(-1),
  finalRadius(-1),
  finalImpurity(std::numeric_limits<double>::infinity()),
  finalUtility(std::numeric_limits<double>::infinity()),
  finalBubble(-1)
  {};

  BubblePopper(const BubblePopper& bubble, Split)
  :
  valid(bubble.valid),
  radius(bubble.radius),
  impurity(bubble.impurity),
  finalValid(bubble.finalValid),
  finalRadius(bubble.finalRadius),
  finalImpurity(bubble.finalImpurity),
  finalUtility(bubble.finalUtility),
  finalBubble(bubble.finalBubble)
  {};

  void operator()(size_t begin, size_t end)
  {
    for (unsigned int i=begin;i<end;i++)
    {

      // We talk about utility here rather than impurity mainly
      // for future developments. We can state a preference for specific
      // node centers and use a weighted average of this preference and
      // impurity to select the best result.
      
      double thisUtility;
      thisUtility = impurity[i];

      bool update = false;

      // Valid code of 2 indicates a priority code for bubble center
      // Currently doesn't do anything

      if((valid[i]==1) & (finalValid!=2))
      {
        if( (thisUtility < finalUtility) |
            ( (thisUtility == finalUtility) & (radius[i]>finalRadius) )
        )
        {
          // Rcout << i ;
          // if(thisUtility < finalUtility)
          // {
          //   Rcout << " Check1 passed : " << thisUtility << "<" << finalUtility;
          // }
          // if((thisUtility == finalUtility) & (radius[i]>finalRadius))
          // {
          //   Rcout << " Check2 passed ";
          // }
          // Rcout << std::endl;
          // 
          update = true;
        }
      }
      else if(valid[i]==2)
      {
        if( (finalValid !=2) |
            (
              (thisUtility < finalUtility) |
              ( (thisUtility == finalUtility) & (radius[i]>finalRadius) )
            )
        )
        update = true;
      }

      if(update)
      {
        finalValid = valid[i];
        finalRadius = radius[i];
        finalImpurity = impurity[i];
        finalUtility = thisUtility;
        finalBubble = i;
      }
    }
  }

  void join(const BubblePopper& rhs)
  {

    bool update = false;

    if(rhs.finalValid==1)
    {
      if( (finalValid!=2) &
          (
          (rhs.finalImpurity < finalImpurity) |
          ( (rhs.finalImpurity == finalImpurity) & (rhs.finalRadius>finalRadius) )
          )
      )
      {
        update=true;
      }
    }
    else if(rhs.finalValid==2)
    {
      if( (finalValid!=2) |
          (rhs.finalImpurity < finalImpurity) |
          ( (rhs.finalImpurity == finalImpurity) & (rhs.finalRadius>finalRadius) )
      )
      {
        update=true;
      }

    }

    if(update)
    {
      finalValid = rhs.finalValid;
      finalRadius = rhs.finalRadius;
      finalImpurity = rhs.finalImpurity;
      finalBubble = rhs.finalBubble;
    }
  }
};



/*
 *  Class for generating a new bubble
 *
 */

struct Bubble {

  int bubbleValid;
  int center;
  double radius;
  int inSize;
  double impurity;

  Bubble():
  bubbleValid(0),
  center(0),
  radius(-1),
  impurity(std::numeric_limits<double>::infinity())
  {}

  void blow(const std::vector<int> &target_discrete,
            const std::vector<double> &target_continuous,
            const std::vector< std::vector<double> > &target_distance,
            const std::vector< std::vector< std::vector<double> > > &D,
            const std::vector<int> &matrixSelection,
            const std::vector<double> &weights,
            const double tolerance,
            const std::vector<int> &idList,
            const std::vector<int> &bubbleCenters
  )
  {
    
    
    // Rcout << "=====================" <<std::endl;
    // Initialise everything needed to blow bubble.
    // Note: we'll iterate over idList_unique, for
    // the purposes of blowing bubble. All non-unique bubbles
    // are treated as if they are invalid
    // To this end, we need to initialise our
    // vectors to contain garbage results

    int groupSize = idList.size();
    int uniqueSize = bubbleCenters.size();

    IntegerVector valid_vec(uniqueSize);
    NumericVector radius_vec(uniqueSize);
    NumericVector impurity_vec(uniqueSize);

    //Blow bubbles
    BubbleBlower bubbles(target_discrete,
                         target_continuous,
                         target_distance,
                         D,
                         matrixSelection,
                         weights,
                         OPTIONS->combinationPower,
                         tolerance,
                         idList,
                         bubbleCenters,
                         valid_vec, radius_vec, impurity_vec,
                         OPTIONS->nObs,
                         OPTIONS->mode,
                         OPTIONS->nCategories,
                         OPTIONS->minBucket,
                         OPTIONS->priors_altered,
                         OPTIONS->classFreq,
                         OPTIONS->loss);
    parallelFor(0, uniqueSize, bubbles);
    
    // Rcout << "AAAAAAAAAAAAAAAAAAAAAA" <<std::endl;
    // 
    // for(int i=0; i< uniqueSize; i++)
    // {
    //   Rcout << bubbleCenters[i] << " ";
    //   Rcout <<valid_vec[i] << " ";
    //   Rcout <<radius_vec[i]  << " ";
    //   Rcout <<impurity_vec[i]  << std::endl;
    // 
    // }
    // 
    // 
    
    
    // Select the best bubble
    BubblePopper finalBubble(valid_vec,
                             radius_vec,
                             impurity_vec
                             );
    parallelReduce(0, uniqueSize, finalBubble);
    

    // Store the results
    this->bubbleValid=finalBubble.finalValid;
    this->center=bubbleCenters[finalBubble.finalBubble];
    this->radius=finalBubble.finalRadius;
    this->impurity=finalBubble.finalImpurity;
    
    
  }
  
};


