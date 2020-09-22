#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

typedef struct{
	double x;
	double y;
}Point;

typedef struct{
	Point A;
	Point B;
}Segment;

typedef struct{
	Point P;
	Point P_prim;
}Coresponding;

typedef struct{
	std::vector<Point> upper_chain;
	std::vector<Point> lower_chain;
}Feasible;

void printPoint(const Point& p){
	cout << "Point ("<<p.x<<","<<p.y<<") ";
}
void printPolygon(const Feasible &f){
	cout <<"Poligon [";
	for (int i=0;i<f.upper_chain.size();i++)
		printPoint(f.upper_chain[i]);
	for (int i=f.lower_chain.size()-1;i>=0;i--)
		printPoint(f.lower_chain[i]);
	cout <<"]"<<endl;
}

const int RIGHT = 1, LEFT = -1, ZERO = 0; 
int directionOfPoint(Point A, Point B, Point P) 
{ 
    // subtracting co-ordinates of point A from 
    // B and P, to make A as origin 
    B.x -= A.x; 
    B.y -= A.y; 
    P.x -= A.x; 
    P.y -= A.y; 
  
    // Determining cross Product 
    int cross_product = B.x * P.y - B.y * P.x; 
  
    // return RIGHT if cross product is positive 
    if (cross_product > 0) 
        return RIGHT; 
  
    // return LEFT if cross product is negative 
    if (cross_product < 0) 
        return LEFT; 
  
    // return ZERO if cross product is zero.  
    return ZERO;  
}

int orientation(Point p1, Point p2, Point p3) 
{ 
    // See 10th slides from following link for derivation 
    // of the formula 
    int val = (p2.y - p1.y) * (p3.x - p2.x) - 
              (p2.x - p1.x) * (p3.y - p2.y); 
  
    if (val == 0) return 0;  // colinear 
  
    return (val > 0)? 1: 2; // clock or counterclock wise 
} 

double distancePoints(Point& A,  Point& B){
	return sqrt((A.x-B.x)*(A.x-B.x) + (A.y-B.y)*(A.y-B.y));
}

void visibilityDetermine(const Segment &e1, const Segment &e1_prim, const Point &concaveV, bool lower, std::vector<Coresponding>& first, std::vector<Coresponding>& second ){
	double x1, y1, x2, y2, x1_prim, y1_prim, x2_prim, y2_prim;
	Coresponding cc, cc1, cc2;
	if (!lower){
	        x1= e1.A.x;
		y1 = e1.A.y;

		x2 = e1.B.x;
		y2 = e1.B.y;

		x1_prim = e1_prim.A.x;
		y1_prim = e1_prim.A.y;

		x2_prim = e1_prim.B.x;
		y2_prim = e1_prim.B.y;
} else{
		x2 = e1.A.x;
		y2 = e1.A.y;

		x1 = e1.B.x;
		y1 = e1.B.y;

		x2_prim = e1_prim.A.x;
		y2_prim = e1_prim.A.y;

		x1_prim = e1_prim.B.x;
		y1_prim = e1_prim.B.y;
}
	double x = concaveV.x;
	double y = concaveV.y;

	cout << "(x1,y1)"<<x1<<" "<<y1<<endl;
cout << "(x2,y2)"<<x2<< " "<<y2<<endl;
cout << "(x1_prim,y1_prim)"<<x1_prim<<" "<<y1_prim<<endl;
cout << "(x2_prim,y2_prim)"<<x2_prim<<" "<<y2_prim<<endl;
cout << "(x,y)"<<x<<" "<<y<<endl;
	double A = (y2 - y1)*(x2_prim - x1_prim) - (y2_prim - y1_prim)*(x2 - x1);
	double B = -1*(x2_prim - x1_prim)*(y - y1) -(y2 - y1)*(x - x1_prim) + (x2 - x1)*(y - y1_prim ) + (y2_prim - y1_prim)*(x - x1);
	double C = (y - y1)*(x - x1_prim) - (y - y1_prim)*(x -x1);
	double diskriminanta = B*B - 4*A*C;
	cout <<"A"<<A<<" B"<< B<<" C"<<C <<"diks"<<diskriminanta<<endl;
	if (diskriminanta < 0)
		cout << "Nema resenja"<<endl;
	if (diskriminanta == 0){
		double alpha = (-1)*B / (2*A);
		cout << "Alfa je "<<alpha<<endl;
		if (alpha < 0 || alpha > 1){
			cout <<"Nije resenje"<<endl;
		} else {
			Point P,P_prim;
			P.x = x1 + alpha*(x2-x1);
			P.y = y1 + alpha*(y2-y1);
			cout <<"P"<<endl;
			printPoint(P);
			P_prim.x = x1_prim + alpha*(x2_prim-x1_prim);
			P_prim.y = y1_prim + alpha*(y2_prim-y1_prim);
			cout <<"P_prim"<<endl;
			printPoint(P_prim);
			cc.P = P;
			cc.P_prim = P_prim;
			if (lower){
				first.push_back(cc);		
			} else {
				second.push_back(cc);
			}
		}
	}
	if (diskriminanta > 0){
		int chosen = 0; //0 nije niko, 1 prvi, 2 drugi, 3 oba
		double alpha1 = ((-1)*B + sqrt(diskriminanta)) / (2*A);
		cout << "Alfa1 je"<<alpha1<<endl;
		if (alpha1 < 0 || alpha1 > 1){
			cout <<"Nije resenje"<<endl;
		} else{
			Point P1,P1_prim;
			P1.x = x1 + alpha1*(x2-x1);
			P1.y = y1 + alpha1*(y2-y1);
			cout <<"P1"<<endl;
			printPoint(P1);
			P1_prim.x = x1_prim + alpha1*(x2_prim-x1_prim);
			P1_prim.y = y1_prim + alpha1*(y2_prim-y1_prim);
			cout <<"P1_prim"<<endl;
			printPoint(P1_prim);	
			cc1.P=P1;
			cc1.P_prim = P1_prim;
			chosen = 1;

		}
		double alpha2 = ((-1)*B - sqrt(diskriminanta)) / (2*A);
		cout << "Alfa2 je"<<alpha2<<endl;
		if (alpha2 < 0 || alpha2 > 1){
			cout <<"Nije resenje"<<endl;
		} else{
			Point P2,P2_prim;
			P2.x = x1 + alpha2*(x2-x1);
			P2.y = y1 + alpha2*(y2-y1);
			cout <<"P2"<<endl;
			printPoint(P2);
			P2_prim.x = x1_prim + alpha2*(x2_prim-x1_prim);
			P2_prim.y = y1_prim + alpha2*(y2_prim-y1_prim);
			cout <<"P2_prim"<<endl;
			printPoint(P2_prim);
			cc2.P=P2;
			cc2.P_prim = P2_prim;
			if (chosen ==0){
				chosen = 2;
			} else{
				chosen = 3;
			}

		}
		if (lower){
			if (chosen == 1){
				first.push_back(cc1);
			} else if (chosen == 2){
				first.push_back(cc2);
			} else{
				second.push_back(cc1);
				first.push_back(cc2);
			}

		} else {
			if (chosen == 1){
				second.push_back(cc1);
			} else if (chosen == 2){
				second.push_back(cc2);
			} else{
				first.push_back(cc1);
				second.push_back(cc2);
			}

		}
	}
	
}
void concaveVertex(const Feasible &f){
	Segment d1,d2;
	d1.A = f.upper_chain[0];
	d1.B = f.lower_chain[0];
	printPoint(d1.A);
	printPoint(d1.B);
	d2.A = f.upper_chain.back();
	d2.B = f.lower_chain.back();
	printPoint(d2.A);
	printPoint(d2.B);

	std::vector<Coresponding> first,second;
	first.reserve(10);
	second.reserve(10);
	cout << "concave vertex upper chain"<<endl;
	for (int i=1; i<f.upper_chain.size()-1; i++){
		if (orientation(f.upper_chain[i-1], f.upper_chain[i], f.upper_chain[i+1])==2){
			printPoint(f.upper_chain[i]);
			visibilityDetermine(d1,d2,f.upper_chain[i], false,first,second);
		}
	}
	cout << "concave vertex lower chain"<<endl;
	for (int i=1; i<f.lower_chain.size()-1; i++){
		if (orientation(f.lower_chain[i-1], f.lower_chain[i], f.lower_chain[i+1])==1){
			printPoint(f.lower_chain[i]);
			visibilityDetermine(d1,d2,f.lower_chain[i], true,first,second);
		}
	}
	cout << endl;
	cout <<"First"<<endl;
	int pmax=0;
	for (int i=0;i<first.size();i++){
		printPoint(first[i].P);
		if (distancePoints(d1.B, first[i].P) > distancePoints(d1.B, first[pmax].P))
			pmax=i;
	}
	cout << "Najudaljenija od first"<<endl;
	printPoint(first[pmax].P);

	int pmax2=0;
	cout <<"Second"<<endl;
	for (int i=0;i<second.size();i++){
		printPoint(second[i].P);
		if (distancePoints(d1.A, second[i].P) > distancePoints(d1.A, second[pmax2].P))
			pmax2=i;
	}
	cout << "Najudaljenija od second"<<endl;
	printPoint(second[pmax].P);
	cout <<endl;
	
}

int main(){
	Feasible f;
	f.upper_chain.reserve(10);
	f.lower_chain.reserve(10);

	Point pomocna;
	pomocna.x=-11;
	pomocna.y=7;
	f.upper_chain.push_back(pomocna);

	
	pomocna.x=-3;
	pomocna.y=6;
	f.upper_chain.push_back(pomocna);


	pomocna.x=4;
	pomocna.y=8;
	f.upper_chain.push_back(pomocna);

	
	pomocna.x=11;
	pomocna.y=6;
	f.upper_chain.push_back(pomocna);

	

	pomocna.x=-8;
	pomocna.y=0;
	f.lower_chain.push_back(pomocna);

	pomocna.x=-5;
	pomocna.y=0;
	f.lower_chain.push_back(pomocna);

	
	pomocna.x=0;
	pomocna.y=1.5;
	f.lower_chain.push_back(pomocna);
	
	pomocna.x=5;
	pomocna.y=1;
	f.lower_chain.push_back(pomocna);
	
	pomocna.x=8;
	pomocna.y=-1;
	f.lower_chain.push_back(pomocna);


	concaveVertex(f);
	return 0;
}

