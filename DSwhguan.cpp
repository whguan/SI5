/***************************************************************** 
*COMP 550 SI5- Breaking Segments containing endpoints

Input : List of all blue Segments and red flag
Output: List of all blue segments, broken so none contains any flag
	points

Date: Oct 2th, 2014 by  Wenhua Guan

******************************************************************/
#include <stdio.h>
#include <iostream>     // std::cout
#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include <limits>
#include <string>
#include <set>

using namespace std;

//Define 1st class: Point
class Point {

	public:	int x, y;
	
	Point(int x0, int y0){
		x = x0;
		y = y0;
	}
	
	~Point(){}
	//Comparison by x, breaking ties by y
	bool LessThan(Point *q){
		return (x < q->x)||((x == q->x) && (y < q->y));
	}
	//Define two points equal
	bool Equal(Point *q){
		return((x == q->x) && (y == q->y));
	}
	// To check whether P0P2 is on the left of P0P1
	static int Orientation(Point *P0, Point *P1, Point *P2){
		long x1 = P2->x - P0->x;
		long y1 = P2->y - P0->y;
		long x2 = P1->x - P0->x;
		long y2 = P1->y - P0->y;

		long area = x1*y2 - x2*y1;

		if (area > 0){
			return 1; 	// left orientation
		}else if (area < 0){
			return -1; 	// right orientation
		}else{
			return 0; 	//colinear
		}
	
	}

};


class Segment{
 
	//Members
	public:	Point  *start, *end;
	
	//constructor
	Segment(Point  *start0, Point *end0){
		start = start0;
		end = end0;
	}
	~Segment(){}

	//Comparision :P is on which Side of the seg
	static int Side(Segment *Seg, Point *P){
		int side = Point :: Orientation(Seg->start, Seg->end, P);

		if(side > 0){
			return 1; 	// below 
		}else if( side < 0){
			return -1; 	// above
		}else{
			return 0;  	// colinear
		}

	}

	//Check intersection
	static bool SEGEMENT_INTERSECT(Segment *Seg1, Segment *Seg2){
		int d1 = Side(Seg1, Seg2->start);
		int d2 = Side(Seg1, Seg2->end);
		int d3 = Side(Seg2, Seg1->start);
		int d4 = Side(Seg2, Seg1->end);

		if(d1*d2<=0 && d3*d4<=0){
			return true;
		}else{
			return false;
		}
	}

	static double MiniusSlope(Point *start, Point *end){
		if(start->x != end->x){
			return -(start->y-end->y)/(start->x-end->x);
		
		}else{
			return numeric_limits<int>::max();
		
		}
	}

};

class flag {
	public: Point *P;
		Segment *Seg;
		int NoSeg;
		int typetag;
		string color;
			
	flag(Point *P0, Segment *Seg0,int NoSeg0, int typetag0, string color0){
		P = P0;
		Seg = Seg0;
		NoSeg = NoSeg0;
		typetag = typetag0;
		color = color0;
	}
	~flag(){}
};

	// Compare flags, return true if flag f0 < f1, vise versa
bool compare(flag *f0,flag *f1){
	if (!(f0->P->Equal(f1->P))){	
		return f0->P->LessThan(f1->P);// Point: Compare by x, breaking ties by y
	}else if (f0->typetag != f1->typetag){
		return f0->typetag < f1->typetag;            // Tyep : terminal < start
	}else if(Segment::MiniusSlope(f0->Seg->start, f0->Seg->end) !=
		     Segment::MiniusSlope(f1->Seg->start, f1->Seg->end) ){
		return Segment :: MiniusSlope(f0->Seg->start, f0->Seg->end)<
		       Segment::MiniusSlope(f1->Seg->start, f1->Seg->end); // -Slope
	}else if(f0->typetag == 1){ //start blue<  start red  
		return f0->color <  f1->color; 
	}else if(f0->typetag == 0){ // red term < blue term
		return f0->color > f1->color;
	}else{
		 cout <<"Overlapping lines of same color error!" <<endl;
		 return false;
	}
}

//EventLess is the operator of the tree structure Set
/*==============================================================================
   seg1 < seg2
   1.if seg1->start->x < seg2->start->x, Side(Seg1,Seg2->start)<0 => seg1 < seg2
   2.if seg1->start->x > seg2->start->x, Side(Seg1,Seg2->start)<0 => seg1 > seg2
   3 same x for start, then compare y value
   4 same start:  slope1 < slope 2 => seg1 < seg2 
==============================================================================*/   
struct EventLess{
	bool operator()(Segment *Seg1, Segment *Seg2){
		//If Seg1 < Seg2, return true; Vise versa.
		if( Seg1->start->x < Seg2->start->x ){
			return (Segment::Side(Seg1,Seg2->start) < 0);
		}else if(Seg1->start->x > Seg2->start->x ){
			return (Segment::Side(Seg2,Seg1->start) > 0);
		}else if(Seg1->start->y != Seg2->start->y){
			return Seg1->start->y < Seg2->start->y;
		}else{
			double Minius_slope1 = Segment:: MiniusSlope(Seg1->start,Seg1->end);
			double Minius_slope2 = Segment:: MiniusSlope(Seg2->start,Seg2->end);			
			return Minius_slope1 > Minius_slope2;
		}
	}
};

int main(int argc, char *argv[]){
	int m, n, k;
	
	//Read coordinates of segments into variables
	if(argc !=3){
		printf("ERROR for running command!");
		exit(EXIT_FAILURE);
	}
	
	FILE *fp;
	fp = fopen(argv[1],"r");
	if(fp ==NULL){
		printf("ERROR: Can't open input file!");
		exit(EXIT_FAILURE);
	}
    //Input Iteration number
	int runtimes = atoi(argv[2]);
	int jointflag = 0; // jointflag = 1 if there is joint segments 
	
	//Begin to read data
	fscanf(fp,"%d %d %d\n",&m, &n, &k );
      	
	vector<Segment *> RedSeg;
	vector<Segment *> BlueSeg;
	vector<flag *> FlagList;
	
	//For tight Sentinel line
	int xmin, ymin, xmax, ymax;
	xmin = numeric_limits<int>::max();
	ymin = numeric_limits<int>::max();
	xmax = numeric_limits<int>::min();
	ymax = numeric_limits<int>::min(); 
	
	for (int i = 0; i < m; i++){
		int x1,y1,x2,y2;
		fscanf(fp,"%d %d %d %d\n", &x1,&y1, &x2, &y2);
	/*	Point *P = new Point(x1, y1), *Q = new Point(x2, y2);
		Segment *s = new Segment(P,Q); 
		flag *start = new flag(P,s,i,1,"R"), *end = new flag(Q,s,i,0,"R");//start: typetag =1; end: typetag=0
		
		RedSeg.push_back(s); 		//  red: color = 0; blue : color = 1; 
		FlagList.push_back(start);   
		FlagList.push_back(end);*/
	}
	for(int j = 0; j < n; j++){
		int x1,y1,x2,y2;
		fscanf(fp,"%d %d %d %d\n", &x1,&y1, &x2, &y2);	
		xmin = min(xmin, x1); xmin = min(xmin, x2);
		xmax = max(xmax, x1); xmax = max(xmax, x2);
		ymin = min(ymin, y1); ymin = min(ymin, y2);
		ymax = max(ymax, y1); ymax = max(ymax, y2);
		
		Point *P = new Point(x1,y1), *Q = new Point(x2,y2);
		Segment *s = new Segment(P,Q);
		flag *start = new flag(P,s,j,1,"B"), *end = new flag(Q,s,j,0,"B");
		BlueSeg.push_back(s);
		FlagList.push_back(start);
		FlagList.push_back(end);

	}

	fclose(fp);
         
	//Run the code below for many times
	for(int k = 0; k < runtimes; k++){

	    // Sort FlagList
	    sort(FlagList.begin(), FlagList.end(),compare);
	    
	    // print FlagList
	    for(vector<int>:: size_type i=0; i != FlagList.size(); ++i ){
	    	cout<<FlagList[i]->NoSeg+1<<FlagList[i]->color<<
			(FlagList[i]->typetag ==0 ? "T":"S")<<endl;
	    }
	    //Define the tree structure for storing  segements as the sweepline moves along x axis
	    set<Segment*, EventLess> Event;
	    set<Segment*, EventLess> :: iterator Above, Below;
        
	    // Tight Sentinel line
	
	    Point P(xmin-1, ymin-1), Q(xmax+1, ymin-1);
	    Point R(xmin-1, ymax+1), S(xmax+1, ymax+1);
	    Segment segbotom(&P,&Q), segtop(&R,&S);
	    Event.insert(&segbotom);
	    Event.insert(&segtop);
	    
	    //Sweep the line through flags
	    for(int i=0; i < FlagList.size(); i++){
			//start
			if(FlagList[i]->typetag == 1){
				Segment seg(FlagList[i]->Seg->start,FlagList[i]->Seg->end);
				Event.insert(&seg);
			
				Below = Event.lower_bound(FlagList[i]->Seg);
				Above = Event.upper_bound(FlagList[i]->Seg);
			
				Below--;
				Segment *segbelow = *(Below);
				Segment *segabove = *(Above);
				if(Segment::Side(segbelow, FlagList[i]->P) == 0){
					Segment *s = new Segment(FlagList[i]->P,segbelow->end);
					BlueSeg.push_back(s);
					cout<<FlagList[i]->NoSeg<< "check"<<endl;
					BlueSeg[FlagList[i]->NoSeg]->end->x = FlagList[i]->P->x;
					BlueSeg[FlagList[i]->NoSeg]->end->y = FlagList[i]->P->y;
				
				} 
				if(Segment::Side(segabove, FlagList[i]->P) == 0) {
					Segment *s = new Segment(FlagList[i]->P,segabove->end);
					BlueSeg.push_back(s);
					cout<<FlagList[i]->NoSeg<< "check"<<endl;
					BlueSeg[FlagList[i]->NoSeg]->end->x = FlagList[i]->P->x;
					BlueSeg[FlagList[i]->NoSeg]->end->y = FlagList[i]->P->y;
					//BlueSeg[FlagList[i]->NoSeg]->end = FlagList[i]->P;

		                }

			}

			//ternimal
			if(FlagList[i]->typetag==0){
				Below = Event.lower_bound(FlagList[i]->Seg);
				Above = Event.upper_bound(FlagList[i]->Seg);
				Below--;
			
				Segment *segbelow = *(Below);
				Segment *segabove = *(Above);

				if(Segment::Side(segbelow, FlagList[i]->P)==0) {
					Segment *s = new Segment(FlagList[i]->P,segbelow->end);
					BlueSeg.push_back(s);
				//	BlueSeg[FlagList[i]->NoSeg]->end = FlagList[i]->P;
					cout<<"1"<<endl;
					BlueSeg[FlagList[i]->NoSeg]->end->x = FlagList[i]->P->x;
					BlueSeg[FlagList[i]->NoSeg]->end->y = FlagList[i]->P->y;
				}

				if(Segment::Side(segabove, FlagList[i]->P)==0){
			     		Segment *s = new Segment(FlagList[i]->P,segabove->end);
					BlueSeg.push_back(s);
				//	BlueSeg[FlagList[i]->NoSeg]->end = FlagList[i]->P;
			                cout<<"2---"<<FlagList[i]->NoSeg<<endl;	
					Segment *s2 = new Segment(FlagList[i]->Seg->start,
							FlagList[i]->P);
				        segabove = s2;
				//	BlueSeg[FlagList[i]->NoSeg]->end->x = FlagList[i]->P->x;
				//	BlueSeg[FlagList[i]->NoSeg]->end->y = FlagList[i]->P->y;
				}

				Event.erase(FlagList[i]->Seg);
			}
	
	    }
	    
	    //No disjoint
	    if(jointflag == 0){
	         cout<<"VERIFIED!"<<endl;
	    }
	}

	for(vector<int>::size_type i = 0; i != BlueSeg.size(); ++i ){
		cout<< BlueSeg[i]->start->x<< " " << BlueSeg[i]->start->y
		<<" "<< BlueSeg[i]->end->x<<" "<< BlueSeg[i]->end->y<<endl; 
	}

	//Delete
	for(vector<int>::size_type i=0; i < RedSeg.size(); i++){ 
		if(RedSeg[i]){	
			delete RedSeg[i]->start;
			delete RedSeg[i]->end;
			delete RedSeg[i];
		}

	}

	for(vector<int>::size_type j=0; j < BlueSeg.size(); j++){
		if(BlueSeg[j]){
			delete BlueSeg[j]->start;
			delete BlueSeg[j]->end;
			delete BlueSeg[j];
		}
	}

	for(vector<int>::size_type i =0; i < FlagList.size(); i++){
		if(FlagList[i]){
			delete FlagList[i];
		}
	}
	
        // Free memory
	vector<Segment *>().swap(RedSeg);
	vector<Segment *>().swap(BlueSeg);
	vector<flag *>().swap(FlagList);


	return 0;

}
