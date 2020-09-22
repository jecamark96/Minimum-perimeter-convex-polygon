#include <iostream>
#include <vector>
#include <cmath>


#define epsilon 0.000001
using namespace std;

typedef struct{
	double x;
	double y;
}Point;

typedef struct{
	Point A;
	Point B;
}Segment;

typedef struct {
	double slope;
	Segment s;
	Segment t;
	Point p_s;
	Point p_t;
}CriticalLine;

//using CriticalLines = std::vector<CriticalLine>;

const int RIGHT = 1, LEFT = -1, ZERO = 0; 

bool differentPoint(Point& M, Point& N){
	return !((abs(M.x-N.x)<epsilon) && (abs(M.y-N.y)<epsilon));
}

int directionOfPoint(Point A, Point B, Point P) 
{ 
    // subtracting co-ordinates of point A from 
    // B and P, to make A as origin 
    B.x -= A.x; 
    B.y -= A.y; 
    P.x -= A.x; 
    P.y -= A.y; 
  
    // Determining cross Product 
    double cross_product = B.x * P.y - B.y * P.x; 
  
    // return RIGHT if cross product is positive 
    if (cross_product > 0) 
        return LEFT; 
  
    // return LEFT if cross product is negative 
    if (cross_product < 0) {
	//cout<<"CROSS PRODUCT "<<cross_product<<"hurrary"<<endl;
        return RIGHT; 
    }
  
    // return ZERO if cross product is zero.  
    return ZERO;  
	//cout<<"NEVERRRRRRRRRR CROSS PRODUCT "<<cross_product<<"hurrary"<<endl;
}

void printPoint(const Point& p){
	cout << "Point ("<<p.x<<","<<p.y<<") ";
}

void printSegment(const Segment& s){
	cout << "SEG["<<endl;
	printPoint(s.A);
	cout << "----------";
	printPoint(s.B);
	cout <<"]"<<endl;
}
void translate_to_origin(const Point &rotate_center, Point &rand){
	rand.x -= rotate_center.x;
	rand.y -= rotate_center.y;
}

void translate(double x, double y, Point &rand){
	rand.x +=x;
	rand.y +=y;
}

double quadrant(double x, double y) 
{ 
  
    if (x > 0 and y > 0){ 
        cout << "lies in First quadrant"; 
	return 1;
    }
  
    else if (x < 0 and y > 0) {
        cout << "lies in Second quadrant"; 
  	return 2;
   }
    else if (x < 0 and y < 0){ 
        cout << "lies in Third quadrant"; 
  	return 3;
   }
    else if (x > 0 and y < 0){ 
        cout << "lies in Fourth quadrant"; 
  	return 4;
   }

} 
  
double findSlope(Point A, Point B){
	cout <<"NEW SESSION"<<endl;
	printPoint(A);
	printPoint(B);
	double dx = A.x;
	double dy = A.y;

	translate(-dx,-dy, A);
	translate(-dx,-dy, B);
	printPoint(A);
	printPoint(B);
	//double solution=((A.x == B.x)*M_PI/2) + (A.x != B.x)*atan((B.y - A.y)/(B.x - A.x));
	if (abs(B.x - A.x)<epsilon){
		cout << "B.x-A.x= "<<B.x-A.x<<endl;
		if (B.x > 0){
			translate(dx,dy,A);
			translate(dx,dy,B);
			cout << "SOLUTION 0"<<endl;
			return 0;
		}
		if (B.x < 0){
			translate(dx,dy,A);
			translate(dx,dy,B);
			cout << "SOLUTION 180"<<endl;
			return 180;
		}
	}	
	if (abs(B.y - A.y)<epsilon){
		if (B.y > 0){
			translate(dx,dy,A);
			translate(dx,dy,B);
			cout << "SOLUTION 90"<<endl;
			return 90;
		}
		if (B.y < 0){
			translate(dx,dy,A);
			translate(dx,dy,B);
			cout << "SOLUTION 270"<<endl;
			return 270;
		}
	}	
	

	double sol= (((A.x == B.x)*M_PI/2) + (A.x != B.x)*atan((B.y - A.y)/(B.x - A.x)))*180/M_PI;
	if (quadrant(B.x, B.y) == 1){
		translate(dx,dy,A);
		translate(dx,dy,B);
		return sol;	
	} else if ((quadrant(B.x, B.y) == 2)||(quadrant(B.x, B.y) == 3)){
			translate(dx,dy,A);
			translate(dx,dy,B);
			return sol + 180;
	} else if (quadrant(B.x, B.y) == 4){
			translate(dx,dy,A);
			translate(dx,dy,B);
			return sol + 360;
	}
	cout <<"___________________________________________"<<endl;
}

double distancePoints(Point M, Point N){
	return abs(sqrt((M.x-N.x)*(M.x-N.x)+ (M.y-N.y)*(M.y-N.y)));
}
int betweenPoint(Point M, Point N, Point P){
	
	if (abs(distancePoints(P,M) + distancePoints(P,N) - distancePoints(M,N))<epsilon)
			return 0;
	if (abs(distancePoints(P,N) + distancePoints(M,N) - distancePoints(M,P))<epsilon){
		return -1;
	}

	if (abs(distancePoints(M,N) + distancePoints(P,M) - distancePoints(P,N))<epsilon){
		return 1;
	}
	

	
}
void setPoint(Point &p1,  const Point &p2){
	p1.x = p2.x;
	p1.y = p2.y;
}
void setSegment(Segment &s1, const Segment &s2){
	setPoint(s1.A, s2.A);
	setPoint(s1.B, s2.B);
}

bool differentSegments(Segment &s1, Segment &s2){
	return (differentPoint(s1.A, s2.A) || differentPoint(s1.A, s2.B)) && 
		(differentPoint(s1.B, s2.A) || differentPoint(s1.A, s2.B));
}
void sort(std::vector<CriticalLine>& cll){
	CriticalLine helpp;
	int i,j,pozmin;	
	for (i=0;i<cll.size()-1;i++){
		pozmin =i;
		for (j=i+1;j<cll.size();j++){
			if ((cll[pozmin].slope) > (cll[j].slope)){
				pozmin =j;
			} 
		}
		helpp = cll[i];
		cll[i] = cll[pozmin];
		cll[pozmin]=helpp;
	}
	std::vector<CriticalLine> cl;
	cl.reserve(100);
	int pos;
	std::copy(cll.begin(), cll.end(), std::back_inserter(cl));

	cll.clear();
	//cout << cl.size();
	i=0; 
	CriticalLine help;
	help = cl[i];
	while((i+1) < cl.size()){
		if (abs(cl[i].slope - cl[i+1].slope)<epsilon){
			pos = betweenPoint(help.p_s, help.p_t, cl[i+1].p_s);
			cout << "points";
			printPoint(help.p_s);
			printPoint(help.p_t);
			printPoint(cl[i+1].p_s);
			cout <<endl<<endl;
			if ((pos == 1)){
				if (differentSegments(help.s, cl[i+1].s)){
				cout<< "treca skroz ispred"<<endl;
				help.p_s = cl[i+1].p_s;
				setSegment(help.s,cl[i+1].s);
				printPoint(help.p_s);
				printPoint(cl[i+1].p_s);
				} else {
					help.p_t = help.p_s;
					setSegment(help.t,help.s);
					help.p_s = cl[i+1].p_s;
					setSegment(help.s,cl[i+1].s);	
					
				}
			}		
			if ((pos == -1)&& differentSegments(help.t, cl[i+1].s)){
				help.p_t = cl[i+1].p_s;
				setSegment(help.t,cl[i+1].s);

			}
			pos = betweenPoint(help.p_s, help.p_t, cl[i+1].p_t);
			
			if ((pos == 1) && differentSegments(help.s, cl[i+1].t)){
				cout<<"HERERRERERER AM AM I"<<endl<<endl;
				help.p_s = cl[i+1].p_t;
				setSegment(help.s,cl[i+1].t);
			}
			if ((pos == -1)&& differentSegments(help.t, cl[i+1].t)){
				cout<<"HERERRERERER AM AM I HELLLO"<<endl<<endl;
				help.p_t = cl[i+1].p_t;
				setSegment(help.t,cl[i+1].t);

			}
			
		} else{
		   cll.push_back(help);
		   help = cl[i+1];
		}
		i++;
	}
	cll.push_back(help);
	/*for (int i=0; i < cll.size();i++)
	{
		cout <<"Line number: "<<i<<endl;
		printPoint(cll[i].p_s);
		printPoint(cll[i].p_t);
		cout << "Segment s"<<endl;
		printPoint(cll[i].s.A);
		printPoint(cll[i].s.B);
		cout << "Segment t"<<endl;
		printPoint(cll[i].t.A);
		printPoint(cll[i].t.B);
		cout << endl;
		cout << cll[i].slope<<endl;
	}*/

}


void setCriticalLine(CriticalLine& cl, const Segment& s, const Segment& t, const Point& p_s, const Point& p_t){
	cl.s = s;
	cl.t = t;
        cl.p_s = p_s;
	cl.p_t = p_t;
        cl.slope=findSlope(p_s, p_t);

}


	

void fixCriticalLine(const Segment& s1, const Segment& s2,const Point& s1a, const Point& s2a, const Point& s1b, const Point& s2b, std::vector<CriticalLine>& criticalLiness){
	CriticalLine cl;
	
	if ((directionOfPoint(s1a,s2a,s1b)==LEFT && directionOfPoint(s1a, s2a, s2b)==RIGHT)
	|| (directionOfPoint(s1a,s2a,s1b)==RIGHT && directionOfPoint(s1a, s2a, s2b)==LEFT))
		{
			cout<<"NO ITS NOT"<<endl;
			cout<<"*************************************************"<<endl;
			printPoint(s1a);
			printPoint(s2a);
			printPoint(s1b);
			return;
			
	
	}	
	if (directionOfPoint(s1a,s2a,s1b) == ZERO){
		if (directionOfPoint(s1a, s2a, s2b) == RIGHT){
				setCriticalLine(cl, s1, s2, s1a, s2a);
				criticalLiness.push_back(cl);
		} else{
				setCriticalLine(cl, s2, s1, s2a, s1a);
				criticalLiness.push_back(cl);
		}
	}
	if (directionOfPoint(s1a, s2a, s1b)==LEFT){
			
			       setCriticalLine(cl, s2, s1, s2a, s1a);
			       criticalLiness.push_back(cl);
			
	} 
	if (directionOfPoint(s1a, s2a, s1b)==RIGHT){	
		
				setCriticalLine(cl, s1, s2, s1a, s2a);
				criticalLiness.push_back(cl);
			
	}
}
void criticalLineTwoSegments(const Segment& s1, const Segment& s2, std::vector<CriticalLine>& criticalLiness){
	fixCriticalLine(s1,s2,s1.A,s2.A, s1.B, s2.B, criticalLiness);
	fixCriticalLine(s1,s2,s1.A,s2.B, s1.B, s2.A, criticalLiness);
	fixCriticalLine(s1,s2,s1.B,s2.B, s1.A, s2.A, criticalLiness);
	fixCriticalLine(s1,s2,s1.B,s2.A, s1.A, s2.B, criticalLiness);
	
	//sort(criticalLiness.begin(), criticalLiness.end(), compareCriticalLines);
	sort(criticalLiness);

	
/*	cout << "TWO BY TWO"<<endl;
	for (int i=0; i < criticalLiness.size();i++)
	{
/*
		cout <<"Line number: "<<i<<endl;
		printPoint(lines[i].p_s);
		printPoint(lines[i].p_t);
		cout << "Segment s"<<endl;
		printPoint(criticalLiness[i].s.A);
		printPoint(criticalLiness[i].s.B);
		cout << "Segment t"<<endl;
		printPoint(criticalLiness[i].t.A);
		printPoint(criticalLiness[i].t.B);
		cout << endl;
		cout << criticalLiness[i].slope<<endl;
	}*/
} 

void criticalLineThreeSegments(const Segment &s1, const Segment &s2, const Segment &s3,
std::vector<CriticalLine>& criticalLines){
	
	std::vector<CriticalLine> pomocni;
	pomocni.reserve(100);
	int pos;
	criticalLineTwoSegments(s1, s2, pomocni);
	for (int i=0;i< pomocni.size();i++){
		if ((directionOfPoint(pomocni[i].p_s, pomocni[i].p_t, s3.A)==LEFT) ||
	     (directionOfPoint(pomocni[i].p_s, pomocni[i].p_t, s3.B)==LEFT)){
		if (directionOfPoint(pomocni[i].p_s, pomocni[i].p_t, s3.A)==ZERO){
			pos = betweenPoint(pomocni[i].p_s, pomocni[i].p_t, s3.A);
			if (pos == 1){
				pomocni[i].p_s=s3.A;
				pomocni[i].s = s3;
			}
			if (pos == -1){
				pomocni[i].p_t= s3.A;
				pomocni[i].t = s3;
			}		
			cout <<"HERE AM I AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
			printPoint(s3.A);
			printPoint(pomocni[i].p_s);	
			printPoint(pomocni[i].p_t);
			cout<<"HERE AM I AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
		}
		if (directionOfPoint(pomocni[i].p_s, pomocni[i].p_t, s3.B)==ZERO){
			pos = betweenPoint(pomocni[i].p_s, pomocni[i].p_t, s3.B);
			if (pos == 1){
				pomocni[i].p_s=s3.B;
				pomocni[i].s = s3;
			}
			if (pos == -1){
				pomocni[i].p_t= s3.B;
				pomocni[i].t = s3;
			}	
			cout <<"HERE AM I AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
			printPoint(s3.B);
			printPoint(pomocni[i].p_s);	
			printPoint(pomocni[i].p_t);
			cout<<"HERE AM I AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";	
		}
		   
		criticalLines.push_back(pomocni[i]);
	}}
	pomocni.clear();

	criticalLineTwoSegments(s1, s3, pomocni);
	for (int i=0;i< pomocni.size();i++){
		if ((directionOfPoint(pomocni[i].p_s, pomocni[i].p_t, s2.A)!=RIGHT) ||
	     (directionOfPoint(pomocni[i].p_s, pomocni[i].p_t, s2.B)!=RIGHT)){
			cout <<"!!!!";
			printPoint(pomocni[i].p_s);
			printPoint(pomocni[i].p_t);
			cout <<endl<<endl<<endl;
			cout << directionOfPoint(pomocni[i].p_s, pomocni[i].p_t, s2.A)<<endl;
			cout << directionOfPoint(pomocni[i].p_s, pomocni[i].p_t, s2.B)<<endl;
			if (directionOfPoint(pomocni[i].p_s, pomocni[i].p_t, s2.A)==ZERO){
			pos = betweenPoint(pomocni[i].p_s, pomocni[i].p_t, s2.A);
			if (pos == 1){
				pomocni[i].p_s=s2.A;
				pomocni[i].s = s2;
			}
			if (pos == -1){
				pomocni[i].p_t= s2.A;
				pomocni[i].t = s2;
			}		

			cout <<"HERE AM I AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
			printPoint(s2.A);
			printPoint(pomocni[i].p_s);	
			printPoint(pomocni[i].p_t);
			cout<<"HERE AM I AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
		}
		if (directionOfPoint(pomocni[i].p_s, pomocni[i].p_t, s2.B)==ZERO){
			pos = betweenPoint(pomocni[i].p_s, pomocni[i].p_t, s2.B);
			if (pos == 1){
				pomocni[i].p_s=s2.B;
				pomocni[i].s = s2;
			}
			if (pos == -1){
				pomocni[i].p_t= s2.B;
				pomocni[i].t = s2;
			}		
			cout <<"HERE AM I AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
			printPoint(s2.B);
			printPoint(pomocni[i].p_s);	
			printPoint(pomocni[i].p_t);
			cout<<"HERE AM I AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
		}
		   
		   criticalLines.push_back(pomocni[i]);
	}}
	pomocni.clear();

	criticalLineTwoSegments(s2, s3, pomocni);
	for (int i=0;i< pomocni.size();i++){
		if ((directionOfPoint(pomocni[i].p_s, pomocni[i].p_t, s1.A)==LEFT) ||
	     (directionOfPoint(pomocni[i].p_s, pomocni[i].p_t, s1.B)==LEFT)){
			if (directionOfPoint(pomocni[i].p_s, pomocni[i].p_t, s1.A)==ZERO){
			pos = betweenPoint(pomocni[i].p_s, pomocni[i].p_t, s1.A);
			if (pos == 1){
				pomocni[i].p_s=s1.A;
				pomocni[i].s = s1;
			}
			if (pos == -1){
				pomocni[i].p_t= s1.A;
				pomocni[i].t = s1;
			}		
			cout <<"HERE AM I AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
			printPoint(s1.A);
			printPoint(pomocni[i].p_s);	
			printPoint(pomocni[i].p_t);
			cout<<"HERE AM I AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
		}
		if (directionOfPoint(pomocni[i].p_s, pomocni[i].p_t, s1.B)==ZERO){
			pos = betweenPoint(pomocni[i].p_s, pomocni[i].p_t, s1.B);
			if (pos == 1){
				pomocni[i].p_s=s1.B;
				pomocni[i].s = s1;
			}
			if (pos == -1){
				pomocni[i].p_t= s1.B;
				pomocni[i].t = s1;
			}		
			cout <<"HERE AM I AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
			printPoint(s1.B);
			printPoint(pomocni[i].p_s);	
			printPoint(pomocni[i].p_t);
			cout<<"HERE AM I AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
		}
		   
		   criticalLines.push_back(pomocni[i]);
	}}
	sort(criticalLines);

	cout << "THREE BY THREE"<<endl;
	for (int i=0; i < criticalLines.size();i++)
	{

		//cout <<"Line number: "<<i<<endl;
		//printPoint(lines[i].p_s);
		//printPoint(lines[i].p_t);
		cout << "Segment s"<<endl;
		printPoint(criticalLines[i].s.A);
		printPoint(criticalLines[i].s.B);
		cout << "Segment t"<<endl;
		printPoint(criticalLines[i].t.A);
		printPoint(criticalLines[i].t.B);
		cout << endl;
		cout << criticalLines[i].slope<<endl;
	}


}
double min(double x, double y){
	return x<y ? x:y;
}

double max(double x, double y){
	return x>y ? x:y;
}

void Case34(CriticalLine &a1, CriticalLine &a2, CriticalLine &b1, CriticalLine &b2){
	if (((a1.slope >= b1.slope) && (a1.slope <= b2.slope)) ||
		((b1.slope >=a1.slope) && (b1.slope <=a2.slope))){
		double phi_lo = max(a1.slope, b1.slope);
		double phi_hi = min(a2.slope, b2.slope);
		
	}
}


bool isCriticalLine(Segment& e,Point &p_s, Point& p_t){
		Point first;
		Point second;
	
		Point vector_st;
		
		vector_st.x = p_t.x-p_s.x;
		vector_st.y = p_t.y-p_s.y;

		first = e.A;
		second.x = first.x+vector_st.x;
		second.y = first.y+vector_st.y;
	
		if (directionOfPoint(first, second, e.B)==LEFT){
			first = e.B;
			second.x = first.x+vector_st.x;
			second.y = first.y+vector_st.y;
		}

		if (directionOfPoint(first, second, p_s) == LEFT){
			return false;
		}

		return true;
}




void MakeTwoOfThem(Segment& a, Segment &b, std::vector<CriticalLine>& result){
	std::vector<CriticalLine> lines;
	lines.reserve(100); 
	result.reserve(100);
	criticalLineTwoSegments(a,b,lines);
	for (int i=0;i<lines.size();i++){
		if (!((differentPoint(lines[i].s.A, a.A)) && (differentPoint(lines[i].s.A, a.B)))){	
		result.push_back(lines[i]);
	
}
	}		
}

void mergeSlope(std::vector<Segment>& segments, std::vector<CriticalLine>& result){
     if (segments.size()==2){
		criticalLineTwoSegments(segments[0], segments[1], result);
	} else if (segments.size()==3){
		criticalLineThreeSegments(segments[0], segments[1], segments[2], result);
	} else{
		std::vector<Segment> A_part;
		std::vector<Segment> B_part;
		A_part.reserve(100);
		B_part.reserve(100);
		int n = segments.size();
		std::copy(segments.begin(), segments.begin() + (n/2), std::back_inserter(A_part));
		std::copy(segments.begin()+(n/2)+1, segments.end(), std::back_inserter(B_part));
		std::vector<CriticalLine> a;
		std::vector<CriticalLine> b;
		std::vector<CriticalLine> candidates;
		candidates.reserve(100);
		a.reserve(100);
		b.reserve(100);
		mergeSlope(A_part, a);
		mergeSlope(B_part, b);
		//cout << "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSIZE "<<a.size()<<" "<<b.size()<<endl;
		cout << "A part"<<endl;
		for (int i=0;i<a.size();i++){
			cout << i<< endl;
			printPoint(a[i].p_s);
			printPoint(a[i].p_t);
			cout << a[i].slope<<endl;
			cout <<"*************"<<endl;
		}
		cout << "B part"<<endl;
		for (int i=0;i<b.size();i++){
			cout << i<< endl;
			printPoint(b[i].p_s);
			printPoint(b[i].p_t);
			cout << b[i].slope<<endl;
			cout <<"*************"<<endl;
		}
		int i=0, j=0;
	
		int n_b = b.size()-1;
		int n_a = a.size()-1;
	
	int current_pair_alphas  = n_a-1;
	int current_pair_betas = n_b-1; 
	int betas = 0;
	int alphas = 0;
	while(i< a.size() && j<b.size()){
		if (b[j].slope < a[(current_pair_alphas+1) % n_a].slope){
			alphas = 0;			
			if (betas == 0){
			     current_pair_betas = (current_pair_betas +1) % n_b;
			     betas = 1;
			}
			if (isCriticalLine(a[(current_pair_alphas+1) % n_a].s, b[j].p_s, b[j].p_t)){
				result.push_back(b[j]);
			}
			MakeTwoOfThem(b[j].t, a[(current_pair_alphas+1) % n_a].s, candidates);
			for (int k=0;k<candidates.size();k++){
				if ((candidates[k].slope > b[j].slope) && 
				    (candidates[k].slope < a[(current_pair_alphas+1) % n_a].slope))
					result.push_back(candidates[k]);
			}
			candidates.clear();
			j++;
			//current_pair_betas = (current_pair_betas +1) % n_b;
		} else {
			betas =0;
			if (alphas==0){
				current_pair_alphas = (current_pair_alphas+1)%n_a;
				alphas =1;
			}
			if (isCriticalLine(b[(current_pair_betas+1) % n_b].s, a[i].p_s, a[i].p_t)){
				result.push_back(a[i]);
			}
			MakeTwoOfThem(a[i].t, b[(current_pair_betas+1) % n_b].s, candidates);
			for (int k=0;k<candidates.size();k++){
				if ((candidates[k].slope > a[i].slope) && 
				    (candidates[k].slope < b[(current_pair_betas+1) % n_b].slope))
					result.push_back(candidates[k]);
			}
			candidates.clear();
			i++;
			//current_pair_alphas = (current_pair_alphas +1) % n_a;
		}
		//return result;
	}
	while (i<a.size()){
			if (isCriticalLine(b[0].s, a[i].p_s, a[i].p_t)){
				result.push_back(a[i]);
			}
			MakeTwoOfThem(b[n_b].t, a[i].s, candidates);
			for (int k=0;k<candidates.size();k++){
				if ((candidates[k].slope > b[n_b].slope) && 
				    ((candidates[k].slope < a[i].slope)|| (candidates[k].slope > a[(i+n_a-1)%n_a].slope)))
					result.push_back(candidates[k]);
			}
			candidates.clear();
			i++;
	}
	
	while (j<b.size()){
			if (isCriticalLine(a[0].s, b[j].p_s, b[j].p_t)){
				result.push_back(b[j]);
			}
			MakeTwoOfThem(a[n_a].t, b[j].s, candidates);
			for (int k=0;k<candidates.size();k++){
				if ((candidates[k].slope > a[n_a].slope) && 
				    (candidates[k].slope < b[j].slope))
					result.push_back(candidates[k]);
			}
			candidates.clear();
			j++;
	}
	/*
	while (i<a.size()){
			if (isCriticalLine(b[0].s, a[i].p_s, a[i].p_t)){
				result.push_back(a[i]);
			}
			MakeTwoOfThem(a[i].t, b[0].s, candidates);
			for (int k=0;k<candidates.size();k++){
				if ((candidates[k].slope > a[i].slope) && 
				    ((candidates[k].slope < b[0].slope)|| (candidates[k].slope > b[n_b].slope)))
					result.push_back(candidates[k]);
			}
			candidates.clear();
			i++;
	}
			

	
	while (j<b.size()){
			if (isCriticalLine(a[0].s, b[j].p_s, b[j].p_t)){
				result.push_back(b[j]);
			}
			MakeTwoOfThem(b[j].t, a[0].s, candidates);
			for (int k=0;k<candidates.size();k++){
				if ((candidates[k].slope > b[j].slope) && 
				    ((candidates[k].slope < a[0].slope)|| (candidates[k].slope > a[n_a].slope)))
					result.push_back(candidates[k]);
			}
			candidates.clear();
			j++;
	}*/
	}
}

int main(){/*
	Segment s1, s2, s3;
	CriticalLine c1, c2;
	s1.A.x = -5;
	s1.A.y = 5;
	s1.B.x = -1;
	s1.B.y = 7;
	
	s2.A.x =-9;
	s2.A.y =3;
	s2.B.x = -1;
	s2.B.y =1;

	s3.A.x = -2;
	s3.A.y = -2;
	s3.B.x = -3;
	s3.B.y = -3;
	
	std::vector<CriticalLine> lines;
	lines.reserve(100);	
	/*criticalLineTwoSegments(s1,s2,lines);
	for (int i=0; i < lines.size();i++)
	{
		cout <<"Line number: "<<i<<endl;
		printPoint(lines[i].p_s);
		printPoint(lines[i].p_t);
		cout << endl;
	}	
	lines.clear();
	criticalLineThreeSegments(s1,s2,s3,lines);
	//sort(lines.begin(), lines.end(), compareCriticalLines);
	for (int i=0; i < lines.size();i++)
	{
		cout <<"Line number: "<<i<<endl;
		printPoint(lines[i].p_s);
		printPoint(lines[i].p_t);
		cout << "Segment s"<<endl;
		printPoint(lines[i].s.A);
		printPoint(lines[i].s.B);
		cout << "Segment t"<<endl;
		printPoint(lines[i].t.A);
		printPoint(lines[i].t.B);
		cout << endl;
		cout << lines[i].slope<<endl;
	}

*/
	std::vector<CriticalLine> critical;
	critical.reserve(100);

	std::vector<Segment> segments;
	segments.reserve(100);

	Segment pomocni;
	pomocni.A.x=-5;
	pomocni.A.y=7;
	pomocni.B.x=-6;
	pomocni.B.y=4;
	segments.push_back(pomocni);

	pomocni.A.x=-5;
	pomocni.A.y=5;
	pomocni.B.x=-1;
	pomocni.B.y=7;
	segments.push_back(pomocni);

	pomocni.A.x=-1;
	pomocni.A.y=6;
	pomocni.B.x=3;
	pomocni.B.y=3;
	segments.push_back(pomocni);

	pomocni.A.x=5;
	pomocni.A.y=5;
	pomocni.B.x=5;
	pomocni.B.y=1;
	segments.push_back(pomocni);

	pomocni.A.x=4;
	pomocni.A.y=2;
	pomocni.B.x=6;
	pomocni.B.y=2;
	segments.push_back(pomocni);

	pomocni.A.x=-3;
	pomocni.A.y=4;
	pomocni.B.x=-2;
	pomocni.B.y=3;
	segments.push_back(pomocni);

	pomocni.A.x=-1;
	pomocni.A.y=4;
	pomocni.B.x=-3;
	pomocni.B.y=3;
	segments.push_back(pomocni);

	pomocni.A.x=-1;
	pomocni.A.y=5;
	pomocni.B.x=4;
	pomocni.B.y=-2;
	segments.push_back(pomocni);

	pomocni.A.x=-9;
	pomocni.A.y=3;
	pomocni.B.x=-1;
	pomocni.B.y=1;
	segments.push_back(pomocni);
	
	pomocni.A.x=-8;
	pomocni.A.y=4;
	pomocni.B.x=-8;
	pomocni.B.y=2;
	segments.push_back(pomocni);
	
	pomocni.A.x=-7.8;
	pomocni.A.y=2.2;
	pomocni.B.x=-5;
	pomocni.B.y=-1;
	segments.push_back(pomocni);	

	pomocni.A.x=-7;
	pomocni.A.y=-1;
	pomocni.B.x=-3;
	pomocni.B.y=-1;
	segments.push_back(pomocni);
	
	pomocni.A.x=-6.3;
	pomocni.A.y=0.9;
	pomocni.B.x=-5.7;
	pomocni.B.y=1.2;
	segments.push_back(pomocni);

	pomocni.A.x=-3;
	pomocni.A.y=-3;
	pomocni.B.x=-2;
	pomocni.B.y=-2;
	segments.push_back(pomocni);

	
	mergeSlope(segments, critical);
	cout << "CRITICAL LINES"<<endl;
	for (int i=0;i<critical.size();i++){
		printPoint(critical[i].p_s);
		printPoint(critical[i].p_t);
		cout << endl;
	}

	Point A,B,C;
	A.x= -3;
	A.y= -5;
	C.x= -4.5;
	C.y= -2;	
	B.x= -6;
	B.y= 1;

	int pos = betweenPoint(A,B,C);
	if (pos == 1){
		cout << "A is between B and C"<<endl;
	} else if (pos == -1){
		cout <<"B is between A and C"<<endl;
	} else if (pos == 0){
		cout <<"C is between A and B"<<endl;

	}/*
	Segment s1, s2, s3;
	s2.A.x =-9;
	s2.A.y =3;
	s2.B.x = -1;
	s2.B.y =1;

	s1.A.x =-8;
	s1.A.y =4;
	s1.B.x = -8;
	s1.B.y =2;

	s3.A.x =-7.8;
	s3.A.y =2.2;
	s3.B.x = -5;
	s3.B.y =-1;
	std::vector<CriticalLine> lines;
	lines.reserve(100);
	criticalLineThreeSegments(s1,s2,s3,lines);
	for (int i=0; i < lines.size();i++)
	{
		cout <<"Line number: "<<i<<endl;
		printPoint(lines[i].p_s);
		printPoint(lines[i].p_t);
		cout << "Segment s"<<endl;
		printPoint(lines[i].s.A);
		printPoint(lines[i].s.B);
		cout << "Segment t"<<endl;
		printPoint(lines[i].t.A);
		printPoint(lines[i].t.B);
		cout << endl;
		cout << lines[i].slope<<endl;
	}
cout<<"NEW LINE"<<endl<<endl<<endl;
	lines.clear();
	criticalLineTwoSegments(s1,s3,lines);
	for (int i=0; i < lines.size();i++)
	{
		cout <<"Line number: "<<i<<endl;
		printPoint(lines[i].p_s);
		printPoint(lines[i].p_t);
		cout << "Segment s"<<endl;
		printPoint(lines[i].s.A);
		printPoint(lines[i].s.B);
		cout << "Segment t"<<endl;
		printPoint(lines[i].t.A);
		printPoint(lines[i].t.B);
		cout << endl;
		cout << lines[i].slope<<endl;
	}

	Point A1;
	A1.x=-7.8;
	A1.y = 2.2;
	Point B1;
	B1.x=-8;
	B1.y = 4;
	Point C1;
	C1.x=-8;
	C1.y = 2;
	Point D1;
	D1.x=-5;
	D1.y = -1;
	cout << directionOfPoint(A1,B1,C1);
	cout << directionOfPoint(A1,B1,D1);
	cout<<endl;*/

	return 0;
}
