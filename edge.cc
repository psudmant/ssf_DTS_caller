#include "edge.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

edge::edge(){
}

bool edge::operator<(edge &other){
    if (delta <other.delta){
        return true;
    }else{
        return false;
    }
}

void edge::recalculate_delta(){
    delta = std::abs(right_median-left_median);
    data = delta;
}

edge::edge(long int _pos, 
           long int _left_epos, 
           long int _right_epos, 
           float _delta, 
           float _left_median, 
           float _right_median){
  pos = _pos;
  left_epos = _left_epos;
  right_epos = _right_epos;
  delta = _delta;
  left_median = _left_median;
  right_median = _right_median;
  data = delta;
}

edge::~edge(){
}

void edge::edge_init(long int _pos, 
                     long int _left_epos, 
                     long int _right_epos, 
                     float _delta, 
                     float _left_median, 
                     float _right_median){
  pos = _pos;
  left_epos = _left_epos;
  right_epos = _right_epos;
  delta = _delta;
  left_median = _left_median;
  right_median = _right_median;
  data = delta;

}
       
void edge::getValStr(int width, char * c){
    sprintf(c,"%*.*f",width,width,delta);
}



