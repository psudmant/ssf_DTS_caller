#include "./heap/heap.h"

class edge: public heap_element<float> { 
    
    public:

        edge();
        edge(long int, long int, long int, float, float, float);
        ~edge();  
        void edge_init(long int, long int, long int, float, float, float);
        void recalculate_delta();
        bool operator<(edge &);
        void getValStr(int, char *);
        long int pos;
        long int left_epos;
        long int right_epos;
        float delta;
        float left_median;
        float right_median;
        edge *right_edge; 
        edge *left_edge; 
};


