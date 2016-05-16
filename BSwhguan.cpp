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
	Point(Point *P0){
		x = P0->x;
		y = P0->y;
	}
	
	~Point(){}

	//Overload operator =
	Point& operator=(const Point& P0){
		x = P0.x;
		y = P0.y;

		return *this;
	
	}
	
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
/*		long x1 = P2->x - P0->x;
		long y1 = P2->y - P0->y;
		long x2 = P1->x - P0->x;
		long y2 = P1->y - P0->y;

		long area = x1*y2 - x2*y1;*/

		long area = (P2->x - P0->x)*(P1->y - P0->y) - 
			    (P1->x - P0->x)*(P2->y - P0->y); 

		if (area > 0){
			return 1; 	// left orientation
		}else if (area < 0){
			return -1; 	// right orientation
		}else{
			return 0; 	//colinear
		}
	
	}

};

class flag;
class Segment{
 
	//Members
	public:	flag  *start, *end;
	
	//constructor
	Segment(flag *start0, flag *end0){
		start = start0;
		end = end0;
	}
	~Segment(){}

	//Comparision :P is on which Side of the seg
	static int Side(Segment *, Point *);

	//Check intersection
	static bool SEGEMENT_INTERSECT(Segment *, Segment *);

	static double MiniusSlope(flag *, flag *);

};

class flag {
	public: Point *P;
		Segment *Seg;
		int NoSeg;
		int typetag;
		string color;
			
	flag(Point *P0,int NoSeg0, int typetag0, string color0){
		P = P0;
		NoSeg = NoSeg0;
		typetag = typetag0;
		color = color0;
	}
	// set value to segment
	void SetSegment(Segment *Seg0){
		Seg = Seg0;	
	}
	~flag(){}
};
int Segment::Side(Segment *Seg, Point *P){
	if(!P->Equal(Seg->start->P)&&(!P->Equal(Seg->end->P))){
		int side = Point :: Orientation(Seg->start->P, Seg->end->P, P);

		if(side > 0){
			return 1; 	// below 
		}else if( side < 0){
			return -1; 	// above
		}else{
			return 0;  	// colinear
		}
	}

}
bool Segment::SEGEMENT_INTERSECT(Segment *Seg1, Segment *Seg2){
		int d1 = Side(Seg1, Seg2->start->P);
		int d2 = Side(Seg1, Seg2->end->P);
		int d3 = Side(Seg2, Seg1->start->P);
		int d4 = Side(Seg2, Seg1->end->P);

		if(d1*d2<=0 && d3*d4<=0){
			return true;
		}else{
			return false;
		}
	}	

double  Segment:: MiniusSlope(flag *start, flag *end){
			 if(start->P->x != end->P->x){
				return -(start->P->y - end->P->y)/(start->P->x - end->P->x);
			 }else{
				return numeric_limits<int>::max();	 
		 	 }

		}

// Compare flags, return true if flag f0 < f1, vise versa
 bool compareflags(const flag *f0, const flag *f1){
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
		if( Seg1->start->P->x < Seg2->start->P->x ){
			return (Segment::Side(Seg1,Seg2->start->P) < 0);
		}else if(Seg1->start->P->x > Seg2->start->P->x ){
			return (Segment::Side(Seg2,Seg1->start->P) > 0);
		}else if(Seg1->start->P->y != Seg2->start->P->y){
			return Seg1->start->P->y < Seg2->start->P->y;
		}else{
			double Minius_slope1 = Segment:: MiniusSlope(Seg1->start,Seg1->end);
			double Minius_slope2 = Segment:: MiniusSlope(Seg2->start,Seg2->end);			
			return Minius_slope1 > Minius_slope2;
		}
	}
};

// Insert new segment if there is end point flagtag->P in the middle of segtag.
void InsertNewSeg(vector<Segment *>& Seg, Segment* segtag, flag* flagtag){
	Point *Pstart = new Point(flagtag->P);
	Point *Pend = new Point(segtag->end->P);
				
	flag *sflag = new flag(Pstart, Seg.size(),1,segtag->start->color);
	flag *eflag = new flag(Pend, Seg.size(),0,segtag->end->color);
				
	Segment *s = new Segment(sflag,eflag);
	sflag->SetSegment(s);
	eflag->SetSegment(s);
				
	Seg.push_back(s);
			
	*(segtag->end->P) = *(flagtag->P);

}

void SweepLine( vector<Segment *>& Seg, const vector<flag *> flaglist, 
	        Segment* SegBelow, Segment* SegTop){
			  
    //Define the tree structure for storing segments as the sweep line moves along x axis
	set<Segment*, EventLess> Event;
	set<Segment*, EventLess> :: iterator Above, Below;
        
    //Insert the sentinel line
	Event.insert(SegBelow);
	Event.insert(SegTop);
	    
	//Sweep the line through flags
	for(int i=0; i < flaglist.size(); i++){
	    //start
		if(flaglist[i]->typetag == 1){
			Event.insert(flaglist[i]->Seg);
			
			Below = Event.lower_bound(flaglist[i]->Seg);
			Above = Event.upper_bound(flaglist[i]->Seg);
			
			Below--;
			Segment *segbelow = *(Below);
			Segment *segabove = *(Above);

			if(Segment::Side(segbelow, flaglist[i]->P) == 0){
				
				InsertNewSeg(Seg, segbelow, flaglist[i]);
			} 
			
			if(Segment::Side(segabove, flaglist[i]->P) == 0) {
				
				InsertNewSeg(Seg, segabove, flaglist[i]);
                        }

		}

	    //ternimal
		if(flaglist[i]->typetag==0){
			Below = Event.lower_bound(flaglist[i]->Seg);
			Above = Event.upper_bound(flaglist[i]->Seg);
			Below--;
			
			Segment *segbelow = *(Below);
			Segment *segabove = *(Above);

			if(Segment::Side(segbelow, flaglist[i]->P)==0) {
				
				InsertNewSeg(Seg, segbelow, flaglist[i]);
			}

			if(Segment::Side(segabove, flaglist[i]->P)==0){
				InsertNewSeg(Seg, segabove, flaglist[i]);
			}

			Event.erase(flaglist[i]->Seg);
		}
	
    }
	Event.clear();

}

int main(int argc, char *argv[]){
	int m, n, k;
	
	//Read coordinates of segments into variables
	if(argc !=3){
		cout<<"ERROR for running command!"<<endl;
		exit(EXIT_FAILURE);
	}
	
	FILE *fp;
	fp = fopen(argv[1],"r");
	if(fp ==NULL){
		cout<<"ERROR: Can't open input file!"<<endl;
		exit(EXIT_FAILURE);
	}
    //Input Iteration number
	int runtimes = atoi(argv[2]);
	int jointflag = 0;    // jointflag = 1 if there is joint segments 
	
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
	ymax = numeric_limits<int>::min(); // */
	
	//1. Read red segment
	for (int i = 0; i < m; i++){
		int x1,y1,x2,y2;
		fscanf(fp,"%d %d %d %d\n", &x1,&y1, &x2, &y2);
		
 		xmin = min(xmin, x1); xmin = min(xmin, x2);
		xmax = max(xmax, x1); xmax = max(xmax, x2);
		ymin = min(ymin, y1); ymin = min(ymin, y2);
		ymax = max(ymax, y1); ymax = max(ymax, y2);
		Point *P = new Point(x1, y1), *Q = new Point(x2, y2);
		//start: typetag =1; end: typetag=0
		flag *start = new flag(P,i,1,"R"), *end = new flag(Q,i,0,"R");
		Segment *s = new Segment(start, end);
		start->SetSegment(s);		end->SetSegment(s);
		
		RedSeg.push_back(s); 		
	//	FlagList.push_back(start);   
	//	FlagList.push_back(end); 
	}
	// Tight Sentinel line for Red Segments: Red below and Red above 
	
 	Point RBS(xmin-1, ymin-1), RBT(xmax+1, ymin-1);  // Red below start/terminal 
	Point RAS(xmin-1, ymax+1), RAT(xmax+1, ymax+1);  // Red above start/terminal
	flag RBSF(&RBS, NULL, 1, "R"), RBTF(&RBT, NULL, 0, "R");  //Red Below start/terminal flag
	flag RASF(&RAS, NULL, 1, "R"), RATF(&RAT, NULL, 0, "R");   // Red above start/terminal flag
	Segment Rsegbotom(&RBSF,&RBTF), Rsegtop(&RASF,&RATF);   //Red bottom/top sentinel line */
	
	// 2. Read blue segments
	xmin = numeric_limits<int>::max();
	ymin = numeric_limits<int>::max();
	xmax = numeric_limits<int>::min();
	ymax = numeric_limits<int>::min(); 
	
	for(int j = 0; j < n; j++){
		int x1,y1,x2,y2;
		fscanf(fp,"%d %d %d %d\n", &x1,&y1, &x2, &y2);	
		xmin = min(xmin, x1); xmin = min(xmin, x2);
		xmax = max(xmax, x1); xmax = max(xmax, x2);
		ymin = min(ymin, y1); ymin = min(ymin, y2);
		ymax = max(ymax, y1); ymax = max(ymax, y2);
		
		Point *P = new Point(x1,y1);
		Point *Q = new Point(x2,y2);
		flag *start = new flag(P,j,1,"B"), *end = new flag(Q, j, 0,"B");  
	        Segment *s = new Segment(start, end);
		start->SetSegment(s); 	end->SetSegment(s);

		BlueSeg.push_back(s);
		FlagList.push_back(start);
		FlagList.push_back(end);

	}

	fclose(fp);
	
	Point BBS(xmin-1, ymin-1), BBT(xmax+1, ymin-1);  // Blue below start/terminal 
	Point BAS(xmin-1, ymax+1), BAT(xmax+1, ymax+1);  // Blue above start/terminal
	flag BBSF(&BBS, m, 1, "B"), BBTF(&BBT, m, 0, "B");  //Blue Below start/terminal flag
	flag BASF(&BAS, m+1, 1, "B"), BATF(&BAT, m+1, 0, "B");   //Blue above start/terminal flag
	Segment Bsegbotom(&BBSF,&BBTF), Bsegtop(&BASF,&BATF);   //Blue bottom/top sentinel line
         
	//Run the code below for many times
	for(int k = 0; k < runtimes; k++){

	    // Sort FlagList
	    sort(FlagList.begin(), FlagList.end(),compareflags);
	    
	    // print FlagList
	    for(vector<int>:: size_type i=0; i != FlagList.size(); ++i ){
	    	cout<<FlagList[i]->NoSeg+1<<FlagList[i]->color<<
			(FlagList[i]->typetag ==0 ? "T":"S")<<endl;
	    }
        //Use sweep line algorithm to break segments
	    SweepLine(BlueSeg, FlagList, &Bsegbotom, &Bsegtop);
	   // SweepLine(RedSeg, FlagList, &Rsegbotom, &Rsegtop);
	}

	for(vector<int>::size_type i = 0; i != BlueSeg.size(); ++i ){
		cout<< BlueSeg[i]->start->P->x<< " " << BlueSeg[i]->start->P->y
		<<" "<< BlueSeg[i]->end->P->x<<" "<< BlueSeg[i]->end->P->y<<endl; 
	}

	//Delete elements in vector
	for(vector<int>::size_type i=0; i < RedSeg.size(); i++){ 
		if(RedSeg[i]){	
			delete RedSeg[i]->start->P;
			delete RedSeg[i]->end->P;
			delete RedSeg[i]->start;
			delete RedSeg[i]->end;
			delete RedSeg[i];
		}

	}

	for(vector<int>::size_type j=0; j < BlueSeg.size(); j++){
		if(BlueSeg[j]){
			delete BlueSeg[j]->start->P;
			delete BlueSeg[j]->end->P;
			delete BlueSeg[j]->start;
			delete BlueSeg[j]->end;
			delete BlueSeg[j];
		}
	}
	
    // Free memory
	vector<Segment *>().swap(RedSeg);
	vector<Segment *>().swap(BlueSeg);
	vector<flag *>().swap(FlagList);

	return 0;

}

	  
