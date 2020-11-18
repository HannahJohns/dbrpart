#include "distanceScorer.hpp"


// Function taken from
// https://gsamaras.wordpress.com/code/random-numbers-%E2%88%88min-max/
//
// Returns random number on [0,n]
int RandomUniform(int n) {
  int top = ((((RAND_MAX - n) + 1) / n) * n - 1) + n;
  int r;
  do {
    r = rand();
  } while (r > top);
  return (r % n);
}



struct NullCalculator: public Worker
{

  // Takes an ID list, a response variable, and a split size.
  // Computes a vector of impurities when selecting randomly (with replacement)
  // from the id List

  const std::vector<int>  target_discrete;
  const std::vector<double> target_continuous;
  const std::vector< std::vector<double> > target_distance;

  const std::vector<int> idList;
  const int splitSize;

  RVector<double> null_statistic;

  const int mode;
  const int nObs;
  const int nCategories;

  const std::vector<double> priors_altered;
  const std::vector<int>  classFreq;

  const std::vector< std::vector<double> > loss;

  NullCalculator(const std::vector<int> target_discrete,
                 const std::vector<double> target_continuous,
                 const std::vector< std::vector<double> > target_distance,

                 const std::vector<int> idList,
                 const int splitSize,
                 NumericVector null_statistic,

                 // Relevant OPTIONS
                 const int mode,
                 const int nObs,
                 const int nCategories,
                 const std::vector<double> priors_altered,
                 const std::vector<int>  classFreq,
                 const std::vector< std::vector<double> > loss
                )
    :
    target_discrete(target_discrete),
    target_continuous(target_continuous),
    target_distance(target_distance),
    idList(idList),
    splitSize(splitSize),
    null_statistic(null_statistic),
    mode(mode),
    nObs(nObs),
    nCategories(nCategories),
    priors_altered(priors_altered),
    classFreq(classFreq),
    loss(loss)
    {};


  void operator()(size_t begin, size_t end)
  {
    unsigned int groupSize = idList.size();

    // Calculate running tallies at baseline, which we'll
    // restore each time we run a permutation.

    int countIn = 0;
    int countOut = groupSize;

    // Classification: in/out group tallies.
    std::vector<int> baseline_tallyInClass(nCategories,0);
    std::vector<int> baseline_tallyOutClass(nCategories,0);
    if(mode==0)
    {
      // OutClass starts with all observations
      for(unsigned int j=0;j<groupSize;j++)
      {
        baseline_tallyOutClass[target_discrete[idList[j]]]++;
      }
    }

    // Regression or DB-Regression: Sum of target variables or sum square distances
    // If Regression: Convert into into mean whenever we need it
    // NOTE: for DB-Regression, baseline_sumOut is SST
    double baseline_sumIn = 0;
    double baseline_sumOut = 0;
    if(mode==1)
    {
      for(unsigned int j=0;j<groupSize;j++)
      {
        baseline_sumOut += target_continuous[idList[j]];
      }
    }
    else if(mode==2)
    {
      for(unsigned int j=0;j<groupSize;j++)
      {
        for(unsigned int k=0;k<j;k++)
        {
          baseline_sumOut += target_distance[idList[j]][idList[k]];
        }
      }
    }


    // Run null tests!
    for(unsigned int i=begin; i<end; i++) // for parallelised null calculation i
    {

      // We're not testing anything, we're just after the value at the specified
      // count. All we need is the risk
      double thisStatistic = 0;

      // First, restore iterative summands

      // Classification: in/out group tallies.
      std::vector<int> tallyInClass(nCategories,0);
      std::vector<int> tallyOutClass(nCategories,0);
      if(mode==0)
      {
        // OutClass starts with all observations
        for(unsigned int j=0;j<nCategories;j++)
        {
          tallyOutClass[j] = baseline_tallyOutClass[j];
        }
      }

      double sumIn = baseline_sumIn;
      double sumOut = baseline_sumOut;

      // Now that everything is initialised, we can get a random null risk

      // Get a list from 0 groupSize-1.
      // This tells us what elements of idList we should select.
      std::vector<int> pickList(groupSize);
      for(int j=0;j<groupSize;j++) pickList[j] = j;

      // Randomly pick elements and update our tallies
      // until we've hit the appropriate number of elements
      for(int j=0;j<splitSize;j++)
      {
        // Sample from list of options without replacement
        unsigned int pick = RandomUniform(pickList.size()-1);
        pickList.erase(pickList.begin() + pick);

        // Keep track of summing variables as we inflate the bubble
        countIn++;
        countOut--;

        if(mode==0)
        {
          tallyInClass[target_discrete[idList[pick]]]++;
          tallyOutClass[target_discrete[idList[pick]]]--;
        }
        else if(mode==1)
        {
          sumIn += target_continuous[idList[pick]];
          sumOut -= target_continuous[idList[pick]];
        }
        else if(mode==2)
        {
          for(int k=0; k<j;k++)
          {
            sumIn += target_distance[idList[pick]][idList[k]];
          }
          for(int k=j; k<groupSize; k++)
          {
            sumOut -= target_distance[idList[pick]][idList[k]];
          }
        }
      }

      // Now that we have a random selection of the same size, we can
      // calculate the test statistic for it
      if(mode==0)
      {
        
        // For classification tasks, we use a Chisq independence test:
        
        // X2 = sum_groups( sum_categories ( ( O - E / E )^2 ))
        // Where O is the tallies in each group and E is the expected tally assuming group and outcome are independent
        
        // Final tallies for each group
        std::vector< std::vector<int> > tallyFinal; 
        tallyFinal.push_back(tallyInClass);
        tallyFinal.push_back(tallyOutClass);
        
        thisStatistic = 0;
        for(int iGroup=0;iGroup<2;iGroup++)
        {
          for(int iCat=0;iCat<nCategories;iCat++)
          {
            if(tallyFinal[0][iCat] + tallyFinal[1][iCat] > 0)
            {
              double E = (
                            tallyFinal[0][iCat] + tallyFinal[1][iCat]
                          )*(
                            iGroup==0 ? splitSize : groupSize - splitSize
                          )/((double) groupSize);
              thisStatistic += (tallyFinal[iGroup][iCat] - E)*(tallyFinal[iGroup][iCat] - E)/E;
            }
          }
        }

      }
      else if(mode==1)
      {
        // Regression Tree

        // Calculate means and sum of squares from means
        // for in and out groups.
        // Impurity is the sum of these.

        double meanIn = sumIn/countIn;
        double SST_in = 0;
        for(int k=0;k<splitSize;k++)
        {
          SST_in += pow(meanIn - target_continuous[idList[k]],2);
        }

        double meanOut = sumOut/countOut;
        double SST_out = 0;
        for(int k=splitSize;k<groupSize;k++)
        {
          SST_out += pow(meanOut - target_continuous[idList[k]],2);
        }
        thisStatistic = SST_in+SST_out;
      }
      else if(mode==2)
      {
        // DB-Regression Tree

        double SSW = sumIn + sumOut;
        double SSB = baseline_sumOut - SSW;
        
        thisStatistic = SSB / SSW * (countIn + countOut - 1);
      }
      null_statistic[i] = thisStatistic;
    } // End parallelised i
  } // End operator
}; // Emd Struct NullCalculator


struct NullTester : public Worker
{
  // Counts the number of times a vector of numbers falls below a threshold
  const double threshold;
  const RVector<double> values;

  int count;

  NullTester(const double threshold,
             const NumericVector values)
    :
    threshold(threshold),
    values(values),
    count(0)
  {};

  NullTester(NullTester& null, Split)
    :
    threshold(null.threshold),
    values(null.values),
    count(0)
  {};

  void operator()(size_t begin, size_t end)
  {
    for (unsigned int i=begin;i<end;i++)
    {
      if(values[i]>threshold)
      {
        count++;
      }
    }
  }

  void join(const NullTester& rhs)
  {
    this->count = (this->count) + rhs.count;
  }
};

double permutationTest(std::vector<int> target_discrete,
                       std::vector<double> target_continuous,
                       std::vector< std::vector<double> > target_distance,
                       std::vector<int> idList,
                       int splitSize,
                       double threshold)
{

  NumericVector null_impurity(OPTIONS->nSims);

  NullCalculator nullDist(target_discrete,
                          target_continuous,
                          target_distance,
                          idList,
                          splitSize,
                          null_impurity,
                          OPTIONS->mode,
                          OPTIONS->nObs,
                          OPTIONS->nCategories,
                          OPTIONS->priors_altered,
                          OPTIONS->classFreq,
                          OPTIONS->loss
                          );
  parallelFor(0, OPTIONS->nSims, nullDist);


  NullTester nullTest(threshold,
                      null_impurity);
  parallelReduce(0, OPTIONS->nSims, nullTest);


  return nullTest.count/(double)OPTIONS->nSims;
}



