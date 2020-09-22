#include <iostream>
#include <cstdlib>
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
	std::vector<Point> upper_chain;
	std::vector<Point> lower_chain;
}Feasible;

double norm(const Segment &d1){
	return sqrt((d1.A.x-d1.B.x)*(d1.A.x-d1.B.x) +
		(d1.A.y-d1.B.y)*(d1.A.y-d1.B.y));
}

double cross_product(const Segment &d1, const Segment &d2){
	return ((d1.B.x - d1.A.x)*(d2.B.y-d2.A.y) -	
		(d1.B.y - d1.A.y)*(d2.B.x-d2.A.x));
}

double dot_product(const Segment &d1, const Segment &d2){
	return ((d1.B.x - d1.A.x)*(d2.B.x-d2.A.x)+
		(d2.B.y-d2.A.y)*(d1.B.y - d1.A.y));
}

double cosinus(const Segment &d1, const Segment &d2){
	return dot_product(d1, d2) / (norm(d1) * norm(d2));
}

double sinus(const Segment &d1, const Segment &d2){
	return abs(cross_product(d1, d2)) / (norm(d1) * norm(d2));
}

void translate_to_origin(const Point &rotate_center, Point &rand){
	rand.x -= rotate_center.x;
	rand.y -= rotate_center.y;
}

void translate(double x, double y, Point &rand){
	rand.x +=x;
	rand.y +=y;
}

void translatePolygon(double x, double y, Feasible &f){
	for (int i=0;i<f.upper_chain.size();i++)
		translate(x,y, f.upper_chain[i]);
	for (int i=0;i<f.lower_chain.size();i++)
		translate(x,y, f.lower_chain[i]);
}
void reflection(Point &P, const Point& A1, const Point& B1){
	//we need a line in form AX+BY+C=0 
	double A = A1.y - B1.y;
	double B = B1.x - A1.x;
	double C = A1.x*B1.y - A1.y*B1.x;
	
	double x = P.x;
	double y = P.y;
		
	P.x = ((B*B-A*A)*x - 2*A*B*y - 2*A*C)/(A*A+B*B);
	P.y = ((A*A-B*B)*y - 2*A*B*x - 2*B*C)/(A*A+B*B);
}

void write_segments(const Segment &d1, const Segment &d2){
	cout<< d1.A.x<<" "<<d1.A.y <<" "<< d1.B.x<<" "<<d1.B.y<<endl;
	cout<< d2.A.x<<" "<<d2.A.y <<" "<< d2.B.x<<" "<<d2.B.y<<endl;
}
void rotation(Segment &d1, Segment &d2){
	double dx =  d1.A.x - d2.A.x;
	double dy =  d1.A.y - d2.A.y;

	d2.A.x = d1.A.x;
	d2.A.y = d1.A.y;
	
	d2.B.x = d2.B.x + dx;
	d2.B.y = d2.B.y + dy;

	write_segments(d1,d2);

	Point d1_a_help = d1.A;
	
	translate_to_origin(d1_a_help, d1.A);
	translate_to_origin(d1_a_help, d1.B);
	translate_to_origin(d1_a_help, d2.A);
	translate_to_origin(d1_a_help, d2.B);
	
	write_segments(d1,d2);
	double x = d2.B.x;
	double y = d2.B.y;

	double cos = cosinus(d1,d2);
	double sin = sinus(d1,d2);
	cout<<"cos "<<cos<<"sin "<<sin<<endl;
	d2.B.x = x * cos - y * sin;
	d2.B.y = y * cos + x * sin;

	write_segments(d1,d2);

	translate(d1_a_help.x, d1_a_help.y, d1.A);
	translate(d1_a_help.x, d1_a_help.y, d1.B);
	translate(d1_a_help.x, d1_a_help.y, d2.A);
	translate(d1_a_help.x, d1_a_help.y, d2.B);
	
	write_segments(d1,d2);
}



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
void TwoByTwoReflect(Feasible &f2){
	Segment reflect_edge;
	reflect_edge.A = f2.upper_chain[0];
	reflect_edge.B = f2.lower_chain[0];

	for (int i=1;i<f2.upper_chain.size();i++){
		reflection(f2.upper_chain[i], reflect_edge.A, reflect_edge.B);
	}
	for (int i=1;i<f2.lower_chain.size();i++){
		reflection(f2.lower_chain[i], reflect_edge.A, reflect_edge.B);
	}
	cout <<"Reflected polygon"<<endl;
	printPolygon(f2);	
	cout <<"-----------------------------------------------------------------"<<endl;
}


void TwoByTwo(Feasible &f1, Feasible &f2){
	//TwoByTwoReflect(f2);

	Point lastA = f1.upper_chain.back();
	Point lastB = f1.lower_chain.back();
	
	Point firstA = f2.upper_chain[0];
	Point firstB = f2.lower_chain[0];
	translatePolygon(lastA.x-firstA.x, lastA.y-firstA.y,f2);
	
	cout<<"Middle"<<endl;
	printPolygon(f2);
	cout <<"*************************************************************"<<endl;

	Segment d1;
	d1.A = lastA;
	d1.B = lastB; 
	Segment d2;
	d2.A = firstA;
	d2.B = firstB;

	double dx =  d1.A.x - d2.A.x;
	double dy =  d1.A.y - d2.A.y;

	translatePolygon(d1.A.x - d2.A.x, d1.A.y - d2.A.y, f2);

	write_segments(d1,d2);

	Point d1_a_help = d1.A;

	translatePolygon(-d1_a_help.x, -d1_a_help.y, f1);
	translatePolygon(-d1_a_help.x, -d1_a_help.y, f2);
	
	double cos = cosinus(d1,d2);
	double sin = sinus(d1,d2);

	double x,y;
	cout<<"cos "<<cos<<"sin "<<sin<<endl;

	for (int i=0;i<f2.upper_chain.size();i++){
		x= f2.upper_chain[i].x;
		y= f2.upper_chain[i].y;
		f2.upper_chain[i].x = x * cos - y * sin;
		f2.upper_chain[i].y = y * cos + x * sin;	
	}
	
	for (int i=0;i<f2.lower_chain.size();i++){
		x= f2.lower_chain[i].x;
		y= f2.lower_chain[i].y;
		f2.lower_chain[i].x = x * cos - y * sin;
		f2.lower_chain[i].y = y * cos + x * sin;	
	}
	

	translatePolygon(d1_a_help.x, d1_a_help.y, f1);
	translatePolygon(d1_a_help.x, d1_a_help.y, f2);
	
        //write_segments(d1,d2);

	cout<<"First polygon"<<endl;
	printPolygon(f1);
	cout<<"Second polygon"<<endl;
	printPolygon(f2);
	
}
void TwoByTwoMove(Feasible &f1, Feasible &f2){
	Point lastA = f1.upper_chain.back();
	Point lastB = f1.lower_chain.back();
	
	Point firstA = f2.upper_chain[0];
	Point firstB = f2.lower_chain[0];
	translatePolygon(lastA.x-firstA.x, lastA.y-firstA.y,f2);
	printPolygon(f2);
}

void TwoByTwoRotation(Feasible &f1, Feasible &f2){
	
}
int main(){
	Segment d1, d2;
	d1.A.x = 2;
	d1.A.y = 2;
	d2.A.x = -1;
	d2.A.y = -1;
	d1.B.x = 2;
	d1.B.y = 5;
	d2.B.x = -3;
	d2.B.y = -1;

	rotation(d1, d2);
	write_segments(d1,d2);

	Feasible f1;
	Feasible f2;

	f1.upper_chain.reserve(10);
	f1.lower_chain.reserve(10);

	Point pomocna;
	//H
	pomocna.x =11;
	pomocna.y =5;
	f1.upper_chain.push_back(pomocna);
	//A
	pomocna.x =6;
	pomocna.y =7;
	f1.upper_chain.push_back(pomocna);
	//G
	pomocna.x =15;
	pomocna.y =5;
	f1.lower_chain.push_back(pomocna);
   	//J
	pomocna.x =9;
	pomocna.y =7;
	f1.lower_chain.push_back(pomocna);

	//B
	pomocna.x =7;
	pomocna.y =10;
	f1.lower_chain.push_back(pomocna);
	
	f2.upper_chain.reserve(10);
	f2.lower_chain.reserve(10);
	
	//C
	pomocna.x =16;
	pomocna.y =11;
	f2.upper_chain.push_back(pomocna);

	//I
	pomocna.x =15;
	pomocna.y =7;
	f2.upper_chain.push_back(pomocna);
	
	//H
	pomocna.x =11;
	pomocna.y =5;
	f2.upper_chain.push_back(pomocna);

	//D
	pomocna.x =18;
	pomocna.y =9;
	f2.lower_chain.push_back(pomocna);
	
	//G
	pomocna.x =15;
	pomocna.y =5;
	f2.lower_chain.push_back(pomocna);
//	printPolygon(f1);
//	printPolygon(f2);

	Feasible f3;
	f3.upper_chain.reserve(10);
	f3.lower_chain.reserve(10);

	//B
	pomocna.x =7;
	pomocna.y =10;
	
	f3.upper_chain.push_back(pomocna);
	//C	
	pomocna.x =16;
	pomocna.y =11;

	f3.upper_chain.push_back(pomocna);

	//A
	pomocna.x =6;
	pomocna.y =7;
	f3.lower_chain.push_back(pomocna);

	//F
	pomocna.x =9;
	pomocna.y =9;
	f3.lower_chain.push_back(pomocna);
	
	//E
	pomocna.x =14;
	pomocna.y =10;
	f3.lower_chain.push_back(pomocna);

	//D
	pomocna.x =18;
	pomocna.y =9;
	f3.lower_chain.push_back(pomocna);
	printPolygon(f3);
	//TwoByTwo(f1,f2);
	//TwoByTwo(f1,f3);
	cout<<"Drugi"<<endl;
	printPolygon(f2);
	cout<<"Treci"<<endl;
	printPolygon(f3);
	cout<<"Prvi"<<endl;
	TwoByTwoReflect(f2);
	TwoByTwo(f3,f2);
	TwoByTwo(f2,f1);
	return 0;


}
