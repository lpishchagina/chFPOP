#ifndef ALGOS_H
#define ALGOS_H

#include "Candidate.h"
#include "pCost.h"
#include "pRectangle.h"


//qhull library
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullVertexSet.h>
#include <libqhullcpp/PointCoordinates.h>

#include <Rcpp.h>
#include "math.h"

/*Class Algos
 * -------------------------------------------------------------------------------
 * OP, PELT, GeomFPOP(rectangle:random/random), chFPOP
 * p - dimension in {2,...,4}
 * -------------------------------------------------------------------------------
 */

using namespace Rcpp;
using namespace std;

//qhull library
//using namespace orgQhull;

template <size_t p>
  
class Algos {
public:
  unsigned int N;
  double Penalty;
  double* VectOfCosts;          //UnpenalizedCost = VectOfCosts[n] - Changes.size()*Penalty
  unsigned int* LastChpt;       //vector of the best last changepoints
  std::vector <unsigned int> Changes;
  std::vector <std::vector <double>> SegmentMeans;
  std::vector <unsigned int> NbOfCandidats;
  double UnpenalizedCost;
  std::array<double,p>* Data;
  
  Algos<p>(Rcpp::NumericMatrix data, double penalty) {
    N = (unsigned int)data.ncol();
    Penalty = penalty;
    VectOfCosts = new double[N + 1];
    LastChpt = new unsigned int[N];
    Data = new std::array<double, p>[N];
    for (size_t  i = 0; i < N; i++) {
      for(size_t k = 0; k < p; k++) { Data[i][k] = data(k,i);}
    }
  }
  
  ~Algos<p>() {
    delete [] VectOfCosts;
    delete [] LastChpt;
    VectOfCosts = NULL;
    LastChpt = NULL;
    delete [] Data;
    Data = NULL;
  }
  
  Algos<p>(const Algos &cand) {
    N = cand.N;
    Penalty = cand.Penalty;
    VectOfCosts = new double[N + 1];
    LastChpt = new unsigned int[N];
    Data = new std::array<double, p>[N];
    for (size_t i = 0; i < N + 1; i++) {VectOfCosts[i] = cand.VectOfCosts[i];}
    for (size_t i = 0; i < N; i++) {
      LastChpt[i] = cand.LastChpt[i];
      for (size_t k = 0; k < p; k++) {Data[i][k] = cand.Data[i][k];}
    }
    Changes = cand.Changes;
    SegmentMeans = cand.SegmentMeans;
    NbOfCandidats = cand.NbOfCandidats;
    UnpenalizedCost = cand.UnpenalizedCost;
  }
  
  //PRUNING-------------------------------------------------------------------//
  void pruning_t(std::list<unsigned int> &changes_t, candidate<p> ** &candidates, unsigned int nbCds) {
    std::list<unsigned int>::iterator iter = changes_t.begin();
    for (size_t i = 0; i < nbCds; i++) {
      if (candidates[*iter]->isPruning) {iter = changes_t.erase(iter);} 
      else { iter++;}
    }
  }
  //RETURN OPTIMAL SEGMENTATION-----------------------------------------------//
  void backtracking() {
    unsigned int chp = N;
    while (chp > 0) {
      Changes.push_back(chp);
      chp = LastChpt[chp-1];
    }
    Changes.push_back(0);
    unsigned int j = 1;
    chp = N - 1;
    while (chp > 0) {
      chp = Changes[j];
      j++;
    }
    reverse(Changes.begin(), Changes.end());//{0,...,N}
    std::vector<double> MeanOneSegment;
    for (unsigned int s = 0; s < Changes.size()-1; s++){
      double coef = Changes[s+1] - Changes[s];
      MeanOneSegment.clear();
      for (unsigned int k = 0; k < p; k++){
        double mean = 0.;
        for (unsigned int count = Changes[s]; count < Changes[s+1]; count++) {mean = mean + Data[count][k];}
        mean = mean/coef;
        MeanOneSegment.push_back(mean);
      }
      SegmentMeans.push_back(MeanOneSegment);
    }
    Changes.pop_back();//remove N
    reverse(Changes.begin(), Changes.end());
    Changes.pop_back();//remove 0
    reverse(Changes.begin(), Changes.end());
    UnpenalizedCost = VectOfCosts[N] - Penalty * (Changes.size());
  }
  //SHOW RESULT---------------------------------------------------------------//
  List ResAlgo(bool showNbCands) {
    List res;
    res["changes"] = Changes;
    res["means"] =  SegmentMeans;
    res["UnpenalizedCost"] = UnpenalizedCost;
    if (showNbCands) {res["NumberOfCandidats"] = NbOfCandidats;}
    return res;
  }
  
  //ALGORITHM-----------------------------------------------------------------//
  List algosOP(unsigned int type_algo, bool showNbCands) {
    srand(time(0));
    //INITIALIZATION------------------------------------//
    candidate<p> ** candidates;
    candidates = new candidate<p> * [N+1];// for all methods
    for (size_t i = 0; i < (N+1); i++) {
      candidates[i] = new candidate<p>(i);
    }
    
    pCost<p> *costDif = new pCost<p>();
    VectOfCosts[0] = 0;
    
    //Initialize list with potential candidates at t
    std::list<unsigned int> changes_t;
    changes_t.push_front(0);
    
    //For GeomFPOP
    pRectangle<p> ** blocks;
    
    // For chFPOP : maxSize, alpha, beta (By default), PointsOnHull and CumSumData
    int maxSize = (p + 2);
    unsigned int alpha = 2; 
    unsigned int beta = 1; 
    double** CumSumData = new double*[N + 1];
    for (unsigned int i = 0; i <= N; i++) 
      CumSumData[i] = new double[p + 1];
    std::vector<unsigned int>* PointsOnHull = new std::vector<unsigned int>[N + 1];
    
    unsigned int index;
    
    
    
    if (type_algo == 2) {//for GeomFPOP
      blocks = new pRectangle<p> * [N+1];
      for (size_t i = 0; i < (N+1); i++) {
        blocks[i] = new pRectangle<p>(i);
      }
    }
    if(type_algo == 3) {//for chFPOP
      //get CumSumData
      for (unsigned int k = 0; k <= p; k++) 
        CumSumData[0][k] = 0;
      for (unsigned int j = 1; j <= N; j++) {
        CumSumData[j][0] = j;
        for (unsigned int k = 1; k <= p; k++) 
          CumSumData[j][k] = CumSumData[j - 1][k] + Data[j-1][k-1];
        //clear points on Hull vector
        for (unsigned int j = 0; j <= N; j++)
          PointsOnHull[j].clear();
      }
    }
    //INITIALIZATION END//
    
    //STEP 1. OPTIMAL PARTITIONING: MIN and ARGMIN//
    pCost<p> * funCtt = new pCost<p>();
    double q_bestTau_t;
    double q_Tau_t;
    //ALGO
    for (size_t t = 0; t < N; t++) {
      // add new data-point in cost functions 
      funCtt->initialize(Data[t]);
      for (std::list<unsigned int>::iterator iter = changes_t.begin(); iter != changes_t.end(); iter++) {
        (candidates[*iter]->cost_tau_t)->addCost(funCtt);
      }
      
      // min
      q_bestTau_t = (candidates[changes_t.front()]->cost_tau_t)->getMin();
      LastChpt[t] = 0;//best solution of last change at time t
      for (std::list<unsigned int>::iterator iter = changes_t.begin(); iter != changes_t.end(); iter++) {
        q_Tau_t = (candidates[*iter]->cost_tau_t)->getMin();
        if(q_Tau_t <= q_bestTau_t) {
          q_bestTau_t = q_Tau_t;
          LastChpt[t] = candidates[*iter]->tau;
        }
      }
      // store min
      VectOfCosts[t+1] = q_bestTau_t + Penalty;//m_t+\beta
      NbOfCandidats.push_back(changes_t.size());
      (candidates[t+1]->cost_tau_t)->initialize(VectOfCosts[t+1]); //add m_t+\beta in cost at time t+1
      
    //STEP 1: END//
    
    
    // STEP 2: PRUNING//
    //PELT
    if (type_algo == 1) oneIterPruningPELT(VectOfCosts[t+1], t, changes_t, candidates);
    //GEOMFPOP
    if (type_algo == 2) oneIterPruningGeomFPOP(t, changes_t,candidates, blocks); 
    //chFPOP
    if ((type_algo == 3) && changes_t.size() > maxSize) {
      //CH pruning
      oneIterPruningChFPOP(VectOfCosts[t+1], t, changes_t, candidates, LastChpt, PointsOnHull, CumSumData); 
      maxSize = floor(alpha * changes_t.size()) + beta;
    }
     
    changes_t.push_back(t+1);//ADD NEW CHANGE
    //STEP 2: END//
    }
    //delete memory
    for (unsigned int i = 0; i <= N; i++)
      delete(CumSumData[i]);
    delete [] CumSumData;
    CumSumData = NULL;
    delete [] PointsOnHull;
    PointsOnHull = NULL;
    //RETURN OPTIMAL SEGMENTATION//
    backtracking();
    return ResAlgo(showNbCands);
  }
  
  //PRUNING-------------------------------------------------------------------//
  void oneIterPruningChFPOP(double boundary, unsigned int t, std::list<unsigned int> &changes_t, candidate<p> ** & candidates, unsigned int* & LastChpt, std::vector<unsigned int>* &PointsOnHull, double** &CumSumData){
    std::list<unsigned int>::iterator iter = changes_t.begin();
    unsigned int index1;
    unsigned int index2;
    //PELT
    for (unsigned int i = 0; i < NbOfCandidats[t]; i++) {
      if ( (candidates[*iter]->cost_tau_t)->getMin() >= boundary) {
        iter = changes_t.erase(iter);
      } else {
        index1 = LastChpt[candidates[(*iter)]->tau];
        if (! PointsOnHull[index1].empty()) {
          PointsOnHull[index1].clear();
        }
        iter++;
      }
    }
    //QHull
    for (auto it = changes_t.begin(); it != changes_t.end(); it++) {
      //do convex hull on (LastChpt[candidates[*it]]],t) and return the points on this hull
      index1 = LastChpt[candidates[(*it)] -> tau];
      if ( PointsOnHull[index1].empty()) {
        std::list<unsigned int>::iterator ch_it = changes_t.begin();
        std::vector<coordT> points;       // coordinate vector of all candidates on interval
        //find left boundary
        while((candidates[(*ch_it)] -> tau) <= index1) {
          ch_it ++;
        }
        //get points
        auto ch_it2 = ch_it;
        while(ch_it2 != changes_t.end()) {
          index2 = candidates[(*ch_it2)] -> tau;
          for (unsigned int k = 0; k <= p; k++)
            points.push_back(CumSumData[index2][k] - CumSumData[index1][k]);
          ch_it2++;
        }
        //get points on hull
        //number of points <= p+1
        if(points.size()/(p+1) <= (p+1)) {
          ch_it2 = ch_it;
          while(ch_it2 != changes_t.end()) {
            PointsOnHull[index1].push_back(candidates[(*ch_it2)] -> tau);
            ch_it2++;
          }
        } else{
          //number of points >= p+2 => do Qhull
          orgQhull::Qhull qhull;
          qhull.runQhull("", p+1, points.size()/(p+1), points.data(), "s");
          //get points on hull
          const orgQhull::QhullVertexList& vertices = qhull.vertexList();
          for (const orgQhull::QhullVertex& vertex : vertices) {
            const double* coords = vertex.point().coordinates();
            int index = int(vertex.point().coordinates()[0] + index1);
            PointsOnHull[index1].push_back(index);
          }
        }
      } else{// Convex Hull exists: check point in PointsonHull
        auto find_it = std::find(PointsOnHull[index1].begin(), PointsOnHull[index1].end(), candidates[(*it)] -> tau);
        if (find_it == PointsOnHull[index1].end()) {
          candidates[(*it)] -> isPruning = true;
        } else{
          candidates[(*it)] -> isPruning = false;
        }
      }
    }
    //clear points on Hull vector
    for (unsigned int j = 0; j <= N; j++)
      PointsOnHull[j].clear();
    iter = changes_t.begin();
    while (iter != changes_t.end()) {
      if (candidates[(*iter)]->isPruning) {
        iter = changes_t.erase(iter);
      } else {
        ++iter;
      }
    }
  }
  //PELT----------------------------------------------------------------------//
  void oneIterPruningPELT(double boundary, unsigned int t, std::list<unsigned int> &changes_t, candidate<p> ** &candidates) {
    
    std::list<unsigned int>::iterator iter = changes_t.begin();
    for (unsigned int i = 0; i < NbOfCandidats[t]; i++) {
      if ( (candidates[*iter]->cost_tau_t)->getMin() >= boundary) {
        iter = changes_t.erase(iter);
      } else { 
        iter++;
      }
    }
  }
  //GEOMFPOP------------------------------------------------------------------//
  void RintersectionGeomFPOP (unsigned int k, bool &isPruning, unsigned int nbCds, unsigned int t, unsigned int * &orderedchanges, candidate<p> ** &candidates, pRectangle<p> ** &blocks, pCost<p> *&costDif, pSphere<p> * &sphere) {
    //last intersection
    costDif->initialize(candidates[orderedchanges[k]]->cost_tau_t, candidates[t+1]->cost_tau_t);
    costDif->getSphere(sphere);
    blocks[orderedchanges[k]]->IntersectionSphere(sphere);
    isPruning = blocks[orderedchanges[k]]->IsEmptyRect();
    
    //others intersections (+random)
    if(!isPruning && (k < (nbCds-1))) {
      unsigned int j = k + (std::rand() % (nbCds - k - 1));//{k,..,t-1}
      costDif->initialize(candidates[orderedchanges[k]]->cost_tau_t, candidates[orderedchanges[j+1]]->cost_tau_t);
      costDif->getSphere(sphere);
      blocks[orderedchanges[k]]->IntersectionSphere(sphere);
      isPruning = blocks[orderedchanges[k]]->IsEmptyRect();
    }
  }
  
  //GEOMFPOP(R-type:last+random/random)----------------------------------------------------------//
  void oneIterPruningGeomFPOP(unsigned int t, std::list<unsigned int> &changes_t, candidate<p> ** &candidates, pRectangle<p> ** &blocks) {
    //INITIALIZATION
    unsigned int nbCds = NbOfCandidats[t];
    pCost<p> *costDif = new pCost<p>();
    pSphere<p> * sphere = new pSphere<p>();
    bool isPruning;
    std::list<unsigned int>::iterator iter = changes_t.begin();
    //copy potential candidates at time t
    unsigned int * orderedchanges = new unsigned int[nbCds]; 
    for (size_t i = 0; i < NbOfCandidats[t]; i++) {
      orderedchanges[i] = *iter;
      iter++;
    }
    //Update candidates
    for (size_t k = 0; k < nbCds; k++) {
      isPruning = false;//candidate is not empty
      RintersectionGeomFPOP (k, isPruning, nbCds,  t,orderedchanges, candidates, blocks, costDif, sphere);
      //exclusion(random)
      if(!isPruning && (k != 0)) {
        unsigned int j =  std::rand() % k; //{0,...,k-1}
        costDif->initialize( candidates[orderedchanges[j]]->cost_tau_t, candidates[orderedchanges[k]]->cost_tau_t);
        costDif->getSphere(sphere);
        blocks[orderedchanges[k]]->ExclusionSphere(sphere);
        isPruning = blocks[orderedchanges[k]]->IsEmptyRect();
      }
      if (isPruning) {candidates[orderedchanges[k]]->isPruning = true;}
    }
    //candidate pruning
    pruning_t(changes_t,candidates, NbOfCandidats[t]);
    delete [] orderedchanges;
    orderedchanges = NULL;
  }
  
  

    
};

#endif //ALGOS_H

 
