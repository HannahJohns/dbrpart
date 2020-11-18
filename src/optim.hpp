#include "permute.hpp"

/*
 *  Optimiser for dbrpart
 * 
 *  Each of these structs should operate as follows:
 *  
 *  1. Ask for values for a series of parameters
 *  2. The rest of the code provides these values
 *  3. From this, the method generates a new request
 *  4. Repeat this until end terminating conditions have been met
 */

struct Soma{
  
  // Control parameters
  
  double pathLength;
  double step;
  double PRT;
  int popSize;
  int nDirections;
  int migrations;
  double minDiv;
  int dims;
  
  int otherValuesCount;
  
  // Parameter bounds
  std::vector<double> parLb;
  std::vector<double> parUb;
  
  // Current Status
  
  int migrationCount; // How many migrations have we done?
  int exitCode;

  // Current position and utility
  std::vector< std::vector<double> > pars_current; // Parameter list
  std::vector< double > utility_current; // The values
  std::vector< std::vector< double > > other_values_current;
  
  // Relatvive fitness of the population
  std::vector<int> rank;
  
  // Where are we going next?
  int nTest; // How many tests do we need to do now?
  std::vector<int> popId; // List of new values to query
  std::vector< std::vector<double> > pars;
  std::vector< double > utility;
  std::vector< std::vector< double > > other_values;
  
  
  Soma(double pathLength,
       double step,
       double PRT,
       int popSize,
       int nDirections,
       int migrations,
       double minDiv,
       int dims,
       int otherValuesCount,
       std::vector<double> parLb,
       std::vector<double> parUb
  )
  {
    
    // Feed parameters in
    
    this->pathLength = pathLength;
    this->step = step;
    this->PRT = PRT;
    this->popSize = popSize;
    this->nDirections = nDirections;
    this->migrations = migrations;
    this->minDiv = minDiv;
    this->dims = dims;
    this->otherValuesCount = otherValuesCount;
    
    for(int i=0;i<dims;i++)
    {
      this->parLb.push_back(parLb[i]);
      this->parUb.push_back(parUb[i]);
    }
    
    this->exitCode = -1; // -1 means still running
    this->migrationCount = 0;
    
    // Generate a random starting population
    
    // Only want to do this if there's a dimension to optimise
    
    if(this->dims==1)
    {
      // When there is only a single dimension to optimise,
      // override some of the basic settings
      
      this->popSize = 3;
      this->PRT = 1;
      this->nDirections = 1;
      
      
      // Empty vals vector for initialisation
      std::vector<double> emptyVals;
      for(int j=0;j<otherValuesCount;j++)
      {
        emptyVals.push_back(-1);
      }
      
      // Don't initialise randomly either
      pars_current = {{0},{2},{2}};
      
      pars_current[0][0] = parLb[0];
      pars_current[1][0] = parUb[0];
      pars_current[2][0] = 0.5*(parLb[0]+parUb[0]);
      
      pars = {{0},{2},{2}};
      pars[0][0] = parLb[0];
      pars[1][0] = parUb[0];
      pars[2][0] = 0.5*(parLb[0]+parUb[0]);
      
      
      for(int i=0; i<this->popSize; i++)
      {
        
        // current pars can be whatever, they'll be overwritten
        popId.push_back(i);
        
        // Set at infinity if not initialised
        utility_current.push_back(std::numeric_limits<double>::infinity());
        utility.push_back(std::numeric_limits<double>::infinity());
        
        other_values_current.push_back(emptyVals);
        other_values.push_back(emptyVals);
        
        // Rank is arbitrary right now
        rank.push_back(i);
      }
      
    }
    else if(this->dims>1)
    {
    
    // Empty vals vector for initialisation
    std::vector<double> emptyVals;
    for(int j=0;j<otherValuesCount;j++)
    {
      emptyVals.push_back(-1);
    }
      
    for(int i=0; i<this->popSize; i++)
    {
      std::vector<double> thesePars;
      for(int j=0; j<dims; j++)
      {
          double roll = (double) std::rand()/ (double) RAND_MAX;
          thesePars.push_back( (parUb[j]-parLb[j])*roll+parLb[j] );
      }
      
      // current pars can be whatever, they'll be overwritten
      pars_current.push_back(thesePars);
      
      popId.push_back(i);
      pars.push_back(thesePars);
      
      // Set at infinity if not initialised
      utility_current.push_back(std::numeric_limits<double>::infinity());
      utility.push_back(std::numeric_limits<double>::infinity());
      
      other_values_current.push_back(emptyVals);
      other_values.push_back(emptyVals);
      
      // Rank is arbitrary right now
      rank.push_back(i);
    }
    }
    else
    {
      // If the function is constant, then initialise return values and force pop size to 1.
      
      this->popSize =1;
      
      // Empty vals vector for initialisation
      std::vector<double> emptyVals;
      for(int j=0;j<otherValuesCount;j++)
      {
        emptyVals.push_back(-1);
      }
    
      utility_current.push_back(std::numeric_limits<double>::infinity());
      utility.push_back(std::numeric_limits<double>::infinity());
      
      other_values_current.push_back(emptyVals);
      other_values.push_back(emptyVals);
      
      rank.push_back(0);
      
    }
    
    this->nTest = this->popSize;
  }
  
  
  // Find the new leader and population locations
  void update() 
  {
    if(this->dims>0)
    {
      
    // Get agent locations
    for(int i=0; i< this->nTest; i++)
    {
      if(this->utility[i] < this->utility_current[this->popId[i]])
      {
        // Update this agent
        this->utility_current[this->popId[i]] = this->utility[i];
        for(int j=0; j < (this->dims); j++)
        {
          this->pars_current[this->popId[i]][j] = pars[i][j];
        }
        for(int j=0; j < (this->otherValuesCount); j++)
        {
          this->other_values_current[this->popId[i]][j] = other_values[i][j];
        }
      }
    }
      
    // Clear the update list now that we're done with it
    this->popId.clear();
    this->pars.clear();
    this->utility.clear();
    this->other_values.clear();
    
    // With agent locations finalised for this step, rank the relative fitness of all
    // agents

    // Bubble sort is hot garbage but in practice this is always O(small) so either sue me
    // or submit a pull request
    bool flag = true;
    while (flag)
    {
      flag = false;
      for(int i=0;i<(this->popSize-1);i++)
      {
        if(utility_current[rank[i]] > utility_current[rank[i+1]])
        {
          int tmp = rank[i+1];
          rank[i+1] = rank[i];
          rank[i] = tmp;
          
          flag=true;
        }
      }
    }
    
    // Update migration count and check for terminating conditions
    
    this->migrationCount++;
    if(std::abs(this->utility_current[this->rank[0]] - this->utility_current[this->rank[this->popSize-1]]) < this->minDiv)
    {
      this->exitCode = 0;
    }
    
    if(migrationCount>this->migrations)
    {
      this->exitCode = 1;
    }
    
    
    }
    else
    {
      // Return success code if no dims to optimise

      this->utility_current[0] = this->utility[0];
      for(int i=0;i<this->otherValuesCount;i++)
      {
        this->other_values_current[0][i] = this->other_values[0][i];
      }
      
      this->exitCode = 0;
    }
  }
  
  
  
  // Get the new locations
  void move()
  {
    // Only move if a leader has been set (no need to move until initial location is set)
    // Also only move if there is a real problem so solve
    if(this->migrationCount>0 & this->dims>0)
    {
      
      // Empty vals vector for initialisation
      std::vector<double> emptyVals;
      for(int j=0;j<otherValuesCount;j++)
      {
        emptyVals.push_back(-1);
      }
      
     // Move all agents
     for(int i=0; i<this->popSize;i++)
     {
       
       //To all agents. Change nDirections to control how many we move towards
       for(int j=0; j<  this->nDirections ; j++)
       {
         if(i!=rank[j]) // We can't move an agent to itself
         {
           
           // Set which directions this agent can move in
           std::vector<int> prtVector;
           int countDims = 0; // Count how many dimensions we move in
           for(int k = 0; k<this->dims;k++)
           {
             double roll = (double) std::rand()/ (double) RAND_MAX;
             if(roll < this->PRT)
             {
               prtVector.push_back(1);
               countDims++;
             }
             else
             {
               prtVector.push_back(0);
             }
           }
           
           
           if(countDims>0)
           { 
             double t= this->step;
             while(t <= this->pathLength)
             {
               
               bool validPars = true; // Constraint check
               std::vector<double> newPars;
               // Rcout << i << " " << rank[j] << "|" << " t=" << t << " | ";
               for(int k=0; k < this->dims; k++)
               {
                 
                 double newPar = this->pars_current[i][k] + (this->pars_current[rank[j]][k] - this->pars_current[i][k])*prtVector[k]*t ;
                 
                 // Rcout << newPar << " ";
                 
                 newPars.push_back( newPar);
                 if( newPars[k] > this->parUb[k] | newPars[k] < this->parLb[k] )
                 {
                   validPars = false;
                 }
               }
               // Rcout << std::endl;
               
               if(validPars)
               {
                 // Only consider this move if it obeys our constraints
                 this->popId.push_back(i);
                 this->utility.push_back(std::numeric_limits<double>::infinity());
                 this->other_values.push_back(emptyVals);
                 this->pars.push_back(newPars);
               }
               else
               {
                 break;
               }
               t = t + this->step;
             }
           } // End if move
         }
       } // End move to these agents
     } // End move agents
    } // End if past migration 0
    
    if(this->dims>0)
    {
      this-> nTest = popId.size();
    }
    else
    {
      this->nTest = 1;
    }
  }
  
  
  void printStatus()
  {
    
    Rcout << std::endl << std::endl;
    Rcout << "Migration count: " << this->migrationCount << std::endl;
    
    if(this->migrationCount>0)
    {
      Rcout << "Ranks: " << " | ";
      for(int i=0;i<this->popSize;i++)
      {
        Rcout << this->rank[i] << " ";
      }
      Rcout << std::endl << std::endl;
      
      
      Rcout << "Current locations:" << std::endl;
      
      for(int i=0; i<this->popSize; i++)
      {
        Rcout << i << " | " << this->utility_current[i] << " | "; 
        for(int j=0;j<this->dims;j++)
        {
          Rcout << this->pars_current[i][j] << " ";
        }
        Rcout << std::endl;
      }
      
    }

    // Rcout << std::endl << std::endl << "About to test:" << std::endl;
    // for(int i=0; i<this->nTest; i++)
    // {
    //   Rcout << this->popId[i] << "| ";
    //   for(int j=0; j<this->dims; j++)
    //   {
    //     Rcout << this->pars[i][j] << " ";
    //   }
    //   Rcout << std::endl;
    // }
    
  }
  
};
  
  
  
  
  