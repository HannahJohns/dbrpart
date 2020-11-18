#include "options.hpp"

/*
 *  Summarises distance information by observation
 */

// Flags unique rows of a set of distance matrices, assuming identity of indiscernibles holds
// The first occurrence of a given observation
// is flagged as unique

struct UniqueFinder: public Worker
{
  
  const std::vector<int>  rows;
  const std::vector< std::vector< std::vector<double> > > matrix; // Input distance matrix
  const std::vector<int> matrixSelection;
  const std::vector<double> tolerance;

  RVector<int> out;
    
  UniqueFinder(const std::vector<int> &rows,
              const std::vector< std::vector< std::vector<double> > > &matrix, // Input distance matrix
              const std::vector<int> &matrixSelection,
              const std::vector<double> &tolerance,
              IntegerVector out
  ) :
    rows(rows),
    matrix(matrix),
    matrixSelection(matrixSelection),
    tolerance(tolerance),
    out(out)
  {};
  
  void operator()(size_t begin, size_t end)
  {
    for(unsigned int i=begin; i<end; i++) // for parallelised row i
    {
      int nRows = rows.size();
      int nMatrices = matrixSelection.size();
      
      // Assume that this record is unique
      out[i] = 1;

      // Check previous records in the list
      for(int j=(i-1); j>=0;j--)
      {
        
        int countMatches = 0;
        for(int k=0;k<nMatrices;k++)
        {
          if(matrix[matrixSelection[k]][rows[i]][rows[j]]<=tolerance[matrixSelection[k]])
          {
            // We found a match on feature matrixSelection[k] 
            countMatches++;
          }
        }
        if(countMatches==nMatrices)
        {
          // If all variables match, mark as not unique and break this for loop
          out[i] = 0;
          break;
        }
      }
    }
  }
};



// Sums rows of distance matrix with specific rows/columns.
// Used to determine how outliery every element is
struct RowSummer: public Worker
{
  std::vector<int>  rows;
  std::vector<int>  columns;
  const std::vector< std::vector< std::vector<double> > > matrix; // Input distance matrix
  const std::vector<int> matrixSelection;
  const std::vector<double> matrixWeights;
  
  RVector<double> out;
  
  RowSummer(std::vector<int> &rows,
            std::vector<int> &columns,
            const std::vector< std::vector< std::vector<double> > > &matrix,
            const std::vector<int> &matrixSelection,
            const std::vector<double> &matrixWeights,
            NumericVector out
  ) :
    rows(rows),
    columns(columns),
    matrix(matrix),
    matrixSelection(matrixSelection),
    matrixWeights(matrixWeights),
    out(out)
  {};
  
  void operator()(size_t begin, size_t end)
  {
    
    int nCols = columns.size();
    int nMatrices = matrixSelection.size();
    
    for(unsigned int i=begin; i<end; i++) // for parallelised row i
    {
      out[i] = 0;
      
      for(unsigned int j=0; j<nCols;j++)
      {
        
        for(unsigned int k = 0; k<nMatrices; k++)
        {
          out[i] = out[i] + matrixWeights[k]*matrix[matrixSelection[k]][rows[i]][columns[j]];
        }
      }
    }
  }
};

// Sums rows of distance matrix with specific rows/columns.
// Used to determine how outliery every element is
struct ThresholdFinder: public Worker
{
  const std::vector<int>  rows;
  const std::vector< std::vector< std::vector<double> > > matrix; // Input distance matrix
  const std::vector<int> matrixSelection;
  const std::vector<double> matrixWeights;
  const double combinationP;
  
  const std::vector<double> tolerance;
  
  double minVal;
  int row;
  int column;
  
  ThresholdFinder( const std::vector<int>  &rows,
            const std::vector< std::vector< std::vector<double> > > &matrix,
            const std::vector<int> &matrixSelection,
            const std::vector<double> &matrixWeights,
            const double combinationP,
            const std::vector<double> &tolerance
  ) :
    rows(rows),
    matrix(matrix),
    matrixSelection(matrixSelection),
    matrixWeights(matrixWeights),
    combinationP(combinationP),
    tolerance(tolerance),
    minVal(std::numeric_limits<double>::infinity()),
    row(-1),
    column(-1)
  {};
  
  
  ThresholdFinder(ThresholdFinder& tmpThresholdFinder, Split):
    rows(tmpThresholdFinder.rows),
    matrix(tmpThresholdFinder.matrix),
    matrixSelection(tmpThresholdFinder.matrixSelection),
    matrixWeights(tmpThresholdFinder.matrixWeights),
    combinationP(tmpThresholdFinder.combinationP),
    tolerance(tmpThresholdFinder.tolerance),
    minVal(std::numeric_limits<double>::infinity()),
    row(-1),
    column(-1)
  {};
  
  void operator()(size_t begin, size_t end)
  {
    int nObs = rows.size();
    int nMatrices = matrixSelection.size();
    
    for(unsigned int i=begin; i<end; i++) // for parallelised row i
    {
      for(unsigned int j=i+1; j<nObs;j++)
      {
        // Check if records are distinct
        
        int countDifferent = 0;
        for(int k=0;k<nMatrices;k++)
        {
          if(matrix[matrixSelection[k]][rows[i]][rows[j]]>tolerance[matrixSelection[k]])
          {
            // We found a difference on feature matrixSelection[k] 
            countDifferent++;
          }
        }
        
        if(countDifferent>0)
        {
          // At least one relevant feature was different
          double d = 0;
          for(int k=0;k<nMatrices;k++)
          {
            d=d+matrixWeights[k]*std::pow(matrix[matrixSelection[k]][rows[i]][rows[j]],combinationP);
          }
          d = std::pow(d,1.0/combinationP);
          
          if(d < minVal)
          {
            minVal = d;
            row = i;
            column = j;
          }
          
        }
        
         
        
      }
    }
  }
  
  void join(const ThresholdFinder& rhs)
  {
    if(rhs.minVal < this->minVal)
    {
      this->minVal = rhs.minVal;
      this->row = rhs.row;
      this->column = rhs.column;
    }
  }
  
};




// Finds the Maximum value in a list of values
struct MaxFinder: public Worker
{
  
  int currentLocation = -1;
  double currentMax = -1;
  
  RVector<double> y;
  
  MaxFinder(const NumericVector y):
    y(y)
  {};
  
  
  MaxFinder(MaxFinder& tmpMaxFinder, Split):
    y(tmpMaxFinder.y)
  {};
  
  void operator()(size_t begin, size_t end)
  {
    for (unsigned int i=begin;i<end;i++)
    {
      if(this->currentLocation==-1 | this->currentMax <= this->y[i])
      {
        this->currentLocation = i;
        this->currentMax = this->y[i];
      }
    }
  }
  
  void join(const MaxFinder& rhs)
  {
    if(this->currentLocation==-1 | this->currentMax <= rhs.currentMax)
    {
      this->currentLocation = rhs.currentLocation;
      this->currentMax = rhs.currentMax;
    }
  }
};








