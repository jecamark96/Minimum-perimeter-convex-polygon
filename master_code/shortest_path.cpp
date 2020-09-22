#include <iostream>
#include <cmath>
#include <vector>
#include <limits>
#include <cfloat>
using namespace std;
#define epsilon 0.000001
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

typedef struct{
	Point P;
	Point P_as;	
}AnchorSeg;

//C++ program to find upper tangent of two polygons. 
#include<bits/stdc++.h> 

  
// stores the center of polygon (It is made 
// global becuase it is used in compare function) 
Point mid; 
  
// determines the quadrant of a point 
// (used in compare()) 
int quad(Point p) 
{ 
    if (p.x >= 0 && p.y >= 0) 
        return 1; 
    if (p.x <= 0 && p.y >= 0) 
        return 2; 
    if (p.x <= 0 && p.y <= 0) 
        return 3; 
    return 4; 
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

// Checks whether the line is crossing the polygon 

// compare function for sorting 
bool compare(Point p1, Point q1) 
{ 
    Point p;
    p.x = p1.x - mid.x; 
    p.y = p1.y - mid.y;
    Point q; 
    q.x = q1.x - mid.x; 
    q.y = q1.y - mid.y; 
  
    int one = quad(p); 
    int two = quad(q); 
  
    if (one != two) 
        return (one < two); 
    return (p.y*q.x < q.y*p.x); 
} 
  
// Finds upper tangent of two polygons 'a' and 'b' 
// represented as two vectors. 
void findLowerTangent(vector<Point> a, 
                      vector<Point> b, Point& a1, Point& b1) 
{ 
    // n1 -> number of points in polygon a 
    // n2 -> number of points in polygon b 
    int n1 = a.size(), n2 = b.size(); 
    cout <<"n1"<<n1<<"n2"<<n2<<endl;
    
    
    // ia -> rightmost point of a 
    int ia = 0, ib = 0; 
    for (int i=1; i<n1; i++) 
        if (a[i].y < a[ia].y) 
            ia = i; 
  
    // ib -> leftmost point of b 
    for (int i=1; i<n2; i++) 
        if (b[i].y > b[ib].y) 
            ib=i; 
  
    // finding the upper tangent 
    int inda = ia, indb = ib; 
    bool done = 0; 
    cout<<"inda"<<inda<<"indb"<<indb<<endl;
    cout << a[inda].x<<" "<<a[inda].y<<endl;
    cout << b[indb].x<<" "<<b[indb].y<<endl;
    while (!done) 
    { 
        done = 1; 
        while (orientation(b[indb], a[inda], a[(n1+inda-1)%n1]) ==2) {
	    cout <<b[indb].x<<"-"<<b[indb].y<<" "<< a[inda].x<<"-"<<a[inda].x<<" "<< a[(n1+inda-1)%n1].x<<"-"<<a[(n1+inda-1)%n1].y<<endl; 
	    cout <<"stack here1 "<<endl;
            inda = (n1+inda-1)%n1;
	    
	}
  
        while (orientation(a[inda], b[indb], b[(indb+1)%n2]) ==2) 
        { 
	    cout <<"stack here2 "<<endl;
	    cout <<a[inda].x<<"-"<<a[inda].y<<" "<< b[indb].x<<"-"<<b[indb].y<<" "<< b[(indb+1)%n2].x<<"-"<<b[(indb+1)%n2].y<<endl;
            indb = (indb+1)%n2; 
            done = 0; 
        } 
    } 
  
    cout << "lower tangent (" << a[inda].x << ","
        << a[inda].y << ") (" << b[indb].x 
        << "," << b[indb].y << ")\n"; 
    a1= a[inda];
    b1= b[indb];
} 


void findUpperTangent(vector<Point> a, 
                      vector<Point> b, Point& a1, Point& b1) 
{ 
    // n1 -> number of points in polygon a 
    // n2 -> number of points in polygon b 
    int n1 = a.size(), n2 = b.size(); 
    cout <<"n1"<<n1<<"n2"<<n2<<endl;
    
    
    // ia -> rightmost point of a 
    int ia = 0, ib = 0; 
    for (int i=1; i<n1; i++) 
        if (a[i].y < a[ia].y) 
            ia = i; 
  
    // ib -> leftmost point of b 
    for (int i=1; i<n2; i++) 
        if (b[i].y > b[ib].y) 
            ib=i; 
  
    // finding the upper tangent 
    int inda = ia, indb = ib; 
    bool done = 0; 
    cout<<"inda"<<inda<<"indb"<<indb<<endl;
    cout << a[inda].x<<" "<<a[inda].y<<endl;
    cout << b[indb].x<<" "<<b[indb].y<<endl;
    while (!done) 
    { 
        done = 1; 
        while (orientation(b[indb], a[inda], a[(inda+1)%n1]) ==2) {
	    cout <<b[indb].x<<"-"<<b[indb].y<<" "<< a[inda].x<<"-"<<a[inda].x<<" "<< a[(inda+1)%n1].x<<"-"<<a[(inda+1)%n1].y<<endl; 
	    cout <<"stack here1 "<<endl;
            inda = (inda + 1) % n1;
	    
	}
  
        while (orientation(a[inda], b[indb], b[(n2+indb-1)%n2]) ==2) 
        { 
	    cout <<"stack here2 "<<endl;
	    cout <<a[inda].x<<"-"<<a[inda].y<<" "<< b[indb].x<<"-"<<b[indb].y<<" "<< b[(n2+indb-1)%n2].x<<"-"<<b[(n2+indb-1)%n2].y<<endl;
            indb = (n2+indb-1)%n2; 
            done = 0; 
        } 
    } 
  
    cout << "upper tangent (" << a[inda].x << ","
        << a[inda].y << ") (" << b[indb].x 
        << "," << b[indb].y << ")\n"; 
    a1= a[inda];
    b1= b[indb];
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


bool lineIntersectionSegment(Point &L1, Point& L2, Point& P1, Point& P2, Point& solution){
	cout<<"Function intersection"<<endl;
	printPoint(L1);
	printPoint(L2);
	printPoint(P1);
	printPoint(P2);
	double A = L1.y - L2.y;
	double B = L2.x - L1.x;
	double C = L1.x*L2.y - L1.y*L2.x;
	double alpha = ((-1)*C-B*P1.y - A*P1.x)/(A*(P2.x-P1.x)+B*(P2.y-P1.y));
	cout <<"A "<< A <<" B "<<B<<" C "<<C<<" alpha "<<alpha;
	if (alpha <0 || alpha > 1)
		return false;
	solution.x = P1.x + alpha*(P2.x-P1.x);
	solution.y = P1.y + alpha*(P2.y-P1.y);
	cout<<"Solution from function"<<endl;
	printPoint(solution);
	return true;
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


double distancePoints(Point& A,  Point& B){
	return sqrt((A.x-B.x)*(A.x-B.x) + (A.y-B.y)*(A.y-B.y));
}


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

void rotation(Segment &d1, Segment &d2, Point& apex2){
	double dx =  d1.A.x - d2.A.x;
	double dy =  d1.A.y - d2.A.y;

	d2.A.x = d1.A.x;
	d2.A.y = d1.A.y;
	
	d2.B.x = d2.B.x + dx;
	d2.B.y = d2.B.y + dy;

	apex2.x=apex2.x + dx;
	apex2.y=apex2.y + dy;
	//write_segments(d1,d2);

	Point d1_a_help = d1.A;
	
	translate_to_origin(d1_a_help, d1.A);
	translate_to_origin(d1_a_help, d1.B);
	translate_to_origin(d1_a_help, d2.A);
	translate_to_origin(d1_a_help, d2.B);
	translate_to_origin(d1_a_help, apex2);
	
	//write_segments(d1,d2);
	double x = d2.B.x;
	double y = d2.B.y;

	double cos = cosinus(d1,d2);
	double sin = sinus(d1,d2);
	cout<<"cos "<<cos<<"sin "<<sin<<endl;
	d2.B.x = x * cos - y * sin;
	d2.B.y = y * cos + x * sin;

	double apex2_x = apex2.x;
	double apex2_y = apex2.y;
	apex2.x = apex2_x * cos - apex2_y * sin;
	apex2.y = apex2_y * cos + apex2_x * sin;
	//write_segments(d1,d2);

	translate(d1_a_help.x, d1_a_help.y, d1.A);
	translate(d1_a_help.x, d1_a_help.y, d1.B);
	translate(d1_a_help.x, d1_a_help.y, d2.A);
	translate(d1_a_help.x, d1_a_help.y, d2.B);
	translate(d1_a_help.x, d1_a_help.y, apex2);
	//write_segments(d1,d2);
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

void solveTriangle(Segment& s1, Segment& s2, Point& apex1, Point& apex2, Point& solution){
	Segment s2_pomocno = s2;
	Point apex2_pomocno = apex2;	
	Point solution_pomocno;
	rotation(s1,s2_pomocno,apex2_pomocno);
	if (directionOfPoint(s1.A, s1.B, apex1) == directionOfPoint(s1.A, s1.B, apex2_pomocno)){
		reflection(apex2_pomocno, s2_pomocno.A, s2_pomocno.B);
	}
	cout <<"SOLVE TRIANGLE"<<endl;
	printPoint(s1.A);
	printPoint(s1.B);
	printPoint(s2_pomocno.A);
	printPoint(s2_pomocno.B);
	printPoint(apex1);
	printPoint(apex2_pomocno);
	cout << endl;
	if (directionOfPoint(apex1, apex2_pomocno, s1.A) == directionOfPoint(apex1,apex2_pomocno, s1.B)){
		if (distancePoints(apex1, s1.A) < distancePoints(apex1,s1.B))
			solution_pomocno = s1.A;
		else 
			solution_pomocno = s1.B;

	} else{
		lineIntersectionSegment(s1.A, s1.B, s2.A, s2.B, solution_pomocno);
	}
	
	rotation(s2_pomocno,s2,solution_pomocno);
	solution = solution_pomocno;
	cout<<"SOLVE TRIANGLE SOLUTION"<<endl;
	printPoint(solution);
	cout <<"++++++++++++++++++++++++++++++++++++++++++++"<<endl;
}

void convexHull(std::vector<Point> points, std::vector<Point>& hull) 
{ 
    // There must be at least 3 points 
    if (points.size() < 3) return; 
  
    // Initialize Result 
   // vector<Point> hull; 
  
    // Find the leftmost point 
    int l = 0; 
    for (int i = 1; i < points.size(); i++) 
        if (points[i].x < points[l].x) 
            l = i; 
  
    // Start from leftmost point, keep moving counterclockwise 
    // until reach the start point again.  This loop runs O(h) 
    // times where h is number of points in result or output. 
    int p = l, q; 
    do
    { 
        // Add current point to result 
        hull.push_back(points[p]); 
  
        // Search for a point 'q' such that orientation(p, x, 
        // q) is counterclockwise for all points 'x'. The idea 
        // is to keep track of last visited most counterclock- 
        // wise point in q. If any point 'i' is more counterclock- 
        // wise than q, then update q. 
        q = (p+1)%points.size(); 
        for (int i = 0; i < points.size(); i++) 
        { 
           // If i is more counterclockwise than current q, then 
           // update q 
           if (orientation(points[p], points[i], points[q]) == 2) 
               q = i; 
        } 
  
        // Now q is the most counterclockwise with respect to p 
        // Set p as q for next iteration, so that q is added to 
        // result 'hull' 
        p = q; 
  
    } while (p != l);  // While we don't come to first point 
  
    // Print Result 
    for (int i = 0; i < hull.size(); i++) 
        cout << "(" << hull[i].x << ", "
              << hull[i].y << ")\n"; 
} 


bool isTangentBetween(std::vector<Point>& Up, std::vector<Point>& Low, Point& ATry, Point& BTry){
	bool up_left = false;
	bool up_right = false;
	
	bool low_left = false;
	bool low_right = false;

	for (int i=0;i<Up.size();i++){
			if (directionOfPoint(ATry, BTry, Up[i])== -1)
				up_left = true;
			if (directionOfPoint(ATry, BTry, Up[i])== 1)
				up_right = true;
	}

	for (int i=0;i<Low.size();i++){
			if (directionOfPoint(ATry, BTry, Low[i])== -1)
				low_left = true;
			if (directionOfPoint(ATry, BTry, Low[i])== 1)
				low_right = true;
	}
	return !((up_left && up_right) || (low_left && low_right));	
}


bool differentPoint(Point& M, Point& N){
	return !((abs(M.x-N.x)<epsilon) && (abs(M.y-N.y)<epsilon));
}

void intersectionTangentPoint(Point& a1, Point& b1, Point& M, Point &N, std::vector<Point>& intersections){
	Point solution;
	if ((differentPoint(a1, M) && differentPoint(a1,N))
		&&  (differentPoint(b1, M) && differentPoint(b1,N)))
		{	if (lineIntersectionSegment(a1, b1, M, N ,solution)){
			printPoint(solution); 			
                        intersections.push_back(solution);
		        //cout <<"ITS SOLUTION";
		        cout << endl;
			}
		}

}
void intersectionEdgePoint(std::vector<Point>& points, Point& M, Point &N, std::vector<Point>& intersections)
{
	Point solution;
	cout << "POINT INTERSECTION"<<endl;
	for (int i=0;i<points.size()-1;i++)
	{
		if ((differentPoint(points[i], M) && differentPoint(points[i],N))
			&&  (differentPoint(points[i+1], M) && differentPoint(points[i+1],N)))
		{	if (lineIntersectionSegment(points[i], points[i+1], M, N ,solution)){
			printPoint(solution); 			
                        intersections.push_back(solution);
		        //cout <<"ITS SOLUTION";
		        cout << endl;
			}
		}
	}
	
}

void invisibleShortestPath(Feasible &f){
	std::vector<Point> concaveUp;
	std::vector<Point> concaveLow;

	cout<<"CONVEX HULL UPPER"<<endl;
	std::vector<Point> UpHull;
	convexHull(f.upper_chain, UpHull);
	std::vector<Point> LowHull;
	cout<<"CONVEX HULL LOWER"<<endl;
	convexHull(f.lower_chain, LowHull);
	cout<<endl;
	sort(UpHull.begin(), UpHull.end(), compare);
	sort(LowHull.begin(), LowHull.end(), compare);
	
	cout << "Upper hull"<<endl;
	for (int i=0;i<UpHull.size();i++){
		printPoint(UpHull[i]);
		cout << endl;
	}

	cout << "Lower hull"<<endl;
	for (int i=0;i<LowHull.size();i++){
		printPoint(LowHull[i]);
		cout << endl;
	}

	Point Up_a, Up_b;
	Point Low_a, Low_b;
        findUpperTangent(UpHull, LowHull,Up_a,Up_b);
	findLowerTangent(UpHull, LowHull,Low_a, Low_b);
	
	concaveUp.reserve(10);
	concaveLow.reserve(10);

	Point A1 = f.upper_chain[0];
	int i;
	for (i=0;i<UpHull.size();i++){
		if ((abs(A1.x-UpHull[i].x)<epsilon) && (abs(A1.y-UpHull[i].y)<epsilon))
			break;
	}
        concaveUp.push_back(UpHull[i]);
	Point B1 = f.upper_chain.back();
	if (orientation(UpHull[i], UpHull[(i+1)%UpHull.size()], UpHull[(i+2)%UpHull.size()])==2){
		cout <<"i am here up"<<endl;
		while (!((abs(B1.x-UpHull[i].x)<epsilon) && (abs(B1.y-UpHull[i].y)<epsilon))){
		concaveUp.push_back(UpHull[(i+1)%UpHull.size()]);
		i = (i+1)%UpHull.size();
		}
	} else if (orientation(UpHull[i], UpHull[(UpHull.size()+i-1)%UpHull.size()], UpHull[(UpHull.size()+i-2)%UpHull.size()])==2){
		cout <<"i am here2 up"<<endl;
		while (!((abs(B1.x-UpHull[i].x)<epsilon) && (abs(B1.y-UpHull[i].y)<epsilon))){
		concaveUp.push_back(UpHull[(UpHull.size()+i-1)%UpHull.size()]);
		i = (UpHull.size()+i-1)%UpHull.size();
	}
	}
	cout << "concaveUp"<<endl;
	for (int i=0;i<concaveUp.size();i++)
		printPoint(concaveUp[i]);
	cout <<endl;

	Point A2 = f.lower_chain[0];
	int k;
	for (k=0;k<LowHull.size();k++){
		if ((abs(A2.x-LowHull[k].x)<epsilon) && (abs(A2.y-LowHull[k].y)<epsilon))
			break;
	}
        concaveLow.push_back(LowHull[k]);
	Point B2 = f.lower_chain.back();
	
	 if (orientation(LowHull[k], LowHull[(k+1)%LowHull.size()], LowHull[(k+2)%LowHull.size()])==1){
		cout <<"i am here"<<endl;
		while (!((abs(B2.x-LowHull[k].x)<epsilon) && (abs(B2.y-LowHull[k].y)<epsilon))){
		concaveLow.push_back(LowHull[(k+1)%LowHull.size()]);
		k = (k+1)%LowHull.size();
		}
	} else if (orientation(LowHull[k], LowHull[(LowHull.size()+k-1)%LowHull.size()], LowHull[(LowHull.size()+k-2)%LowHull.size()])==1){
		cout <<"i am here2"<<endl;
		while (!((abs(B2.x-LowHull[k].x)<epsilon) && (abs(B2.y-LowHull[k].y)<epsilon))){
		concaveLow.push_back(LowHull[(LowHull.size()+k-1)%LowHull.size()]);
		k = (LowHull.size()+k-1)%LowHull.size();
	}
	}	
	
	cout << "concaveLow"<<endl;
	for (int i=0;i<concaveLow.size();i++)
		printPoint(concaveLow[i]);
	cout <<endl;
	
	std::vector<Point> intersection_left;
	std::vector<Point> intersection_right;
	intersection_left.reserve(10);
	intersection_right.reserve(10);
	intersectionEdgePoint(concaveUp, A1, A2, intersection_left);
	intersectionEdgePoint(concaveLow, A1, A2, intersection_left);
	intersectionEdgePoint(concaveUp, B1, B2, intersection_right);
	intersectionEdgePoint(concaveLow, B1, B2, intersection_right);
	intersectionTangentPoint(Up_a, Up_b, A1, A2, intersection_left);
	intersectionTangentPoint(Up_a, Up_b, B1, B2, intersection_right);
	if (differentPoint(Up_a, Low_a) || differentPoint(Up_b, Low_b)){
		intersectionTangentPoint(Low_a, Low_b, A1, A2, intersection_left);
		intersectionTangentPoint(Low_a, Low_b, B1, B2, intersection_right);
	}
	
	cout << "INTERSECTIONS"<<endl;
	for (int i=0;i<intersection_left.size();i++)
	{
		printPoint(intersection_left[i]);
		cout << endl;
	}
	cout <<"BETWEEN THEM"<<endl;
	for (int i=0;i<intersection_right.size();i++)
	{
		printPoint(intersection_right[i]);
		cout << endl;
	}	
	
	
	
}
void visibleShortestPath(Segment &d1, Segment &d2){
	double x1,y1,x2,y2,x1_prim, y1_prim, x2_prim, y2_prim;
	x1 = d1.B.x;
	y1 = d1.B.y;

	x2 = d1.A.x;
	y2 = d1.A.y;

	x1_prim = d2.B.x;
	y1_prim = d2.B.y;

	x2_prim = d2.A.x;
	y2_prim = d2.A.y;

	double c = (-1)*(
(x1 - x1_prim)*(-1*x1 + x2 + x1_prim - x2_prim) + (y1- y1_prim)*(-1*y1 + y2 + y1_prim -y2_prim)
)
/
((-1*x1 + x2 + x1_prim - x2_prim)*(-1*x1 + x2 + x1_prim - x2_prim)+
		(-1*y1 + y2 + y1_prim - y2_prim)*(-1*y1 + y2 + y1_prim - y2_prim));

	cout <<"c "<<c<<endl;
	Point first;
	Point second;
	if (c < 0){
		first.x=x1;
		first.y=y1;

		second.x = x1_prim;
		second.y = y1_prim;
	} else if (c>1){
		first.x=x2;
		first.y=y2;

		second.x = x2_prim;
		second.y = y2_prim;

	} else{
		first.x=x1+c*(x2-x1);
		first.y=y1+c*(y2-y1);

		second.x = x1_prim+c*(x2_prim-x1_prim);
		second.y = y1_prim+c*(y2_prim-y1_prim);
	}
	printPoint(first);
	printPoint(second);
	cout <<"Length "<<distancePoints(first, second)<<endl;

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
	printPoint(second[pmax2].P);
	cout <<endl;

	Segment v1;
	v1.A = first[pmax].P;
	v1.B = second[pmax2].P;
	Segment v1_prim;
	v1_prim.A = first[pmax].P_prim;
	v1_prim.B = second[pmax2].P_prim;
	visibleShortestPath(v1,v1_prim);

	
}


int main(){
	Feasible f;
	f.upper_chain.reserve(10);
	f.lower_chain.reserve(10);

	Point pomocna;
	pomocna.x=-7;
	pomocna.y=5;
	f.upper_chain.push_back(pomocna);

	
	pomocna.x=-4;
	pomocna.y=2;
	f.upper_chain.push_back(pomocna);


	pomocna.x=-1;
	pomocna.y=6;
	f.upper_chain.push_back(pomocna);

	
	pomocna.x=3;
	pomocna.y=5;
	f.upper_chain.push_back(pomocna);

	

	pomocna.x=4;
	pomocna.y=4;
	f.upper_chain.push_back(pomocna);

	pomocna.x=-10;
	pomocna.y=0;
	f.lower_chain.push_back(pomocna);

	
	pomocna.x=-6;
	pomocna.y=-3;
	f.lower_chain.push_back(pomocna);
	
	pomocna.x=-2;
	pomocna.y=1;
	f.lower_chain.push_back(pomocna);
	
	pomocna.x=3;
	pomocna.y=1;
	f.lower_chain.push_back(pomocna);
	
	pomocna.x=7;
	pomocna.y=-1;
	f.lower_chain.push_back(pomocna);

	//concaveVertex(f);
	invisibleShortestPath(f);

	Segment s1,s2;
	Point apex1, apex2, sol;
	s1.A.x=-5;
	s1.A.y=5;
	s1.B.x=-4;
	s1.B.y=3;

	s2.A.x=5;
	s2.A.y=5;
	s2.B.x=4;
	s2.B.y=3;

	apex1.x=-7;
	apex1.y=3;
	apex2.x=3;
	apex2.y=4;

	//solveTriangle(s1,s2,apex1,apex2,sol);
	return 0;
}

