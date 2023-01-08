
/*
C++ functions to support fast compution in R main file "Dead_Reckoning_2018-LI.R"
*/

#include <cmath>
#include <RcppArmadillo.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]


// Max and min macros
// https://stackoverflow.com/a/3437484
#define Max(a,b) \
 ({ __typeof__ (a) _a = (a); \
		 __typeof__ (b) _b = (b); \
	 _a > _b ? _a : _b; })
#define Min(a,b) \
 ({ __typeof__ (a) _a = (a); \
		 __typeof__ (b) _b = (b); \
	 _a < _b ? _a : _b; })
	 

// Max integer value in R minus 1
#define MAX_INT_R (int)(2147483646) // (2147483647-1)

// Struct "Point" consits of 2 coordinates x and y
struct Point{double x; double y; };
#define dist(v,w) sqrt( pow(v.x-w.x,2.0)+pow(v.y-w.y,2.0) )
#define Dot(v,w) v.x*w.x + v.y*w.y

bool pnpoly_test(Rcpp::NumericVector&, Rcpp::NumericVector&, double, double);


// Return minimum distance between line segment uv and Point p
// [[Rcpp::export]]
double min_distance(arma::vec& P, Rcpp::NumericVector& U, Rcpp::NumericVector& V) {
	
	Point p = {P(0),P(1)}, u, v;
	u = {U(0),U(1)}; v = {V(0),V(1)};
	
	// Return minimum distance between line segment uv and Point p
	const double l2 = dist(v,u);  // i.e. |u-v|^2
	if (l2 < 10e-14) return dist(p,u);   // v == u case
	
	Rcpp::Rcout << "dist uv: " << l2 << "\n";
	Rcpp::Rcout << "dist to u: " << dist(p,u) << "\n";
	Rcpp::Rcout << "dist to v: " << dist(p,v) << "\n";
	
	// Shoelace formula 
	// [https://en.wikipedia.org/wiki/Shoelace_formula#Proof_for_a_triangle]
	double A2 = std::abs((p.x-u.x)*(p.y-v.y)-(p.x-v.x)*(p.y-u.y)) ;	
	// double s, a1, a2, b1, b2;
	// a1 = p.x-u.x; b1 = p.y-u.y;
	// a2 = p.x-v.x; b2 = p.y-v.y;
	// s = std::abs(a1*b2-a2*b1);
	
	Rcpp::Rcout << "Area: " << A2 << "\n";	
	
	return Max(A2/l2,Min(dist(p,u),dist(p,v)) ) ;
}



// Calculate the distances between ONE whale position (P) to ONE polygon (coordinates X, Y)
// [[Rcpp::export]]
Rcpp::NumericVector distCoast_test (Rcpp::NumericVector& P, 
									Rcpp::NumericVector& X, Rcpp::NumericVector& Y, 
									int coast_ID = 2147483646 ) 
{
	// Here coast_ID is the one in R minus 1, 
	// because C/C++ array/vector start from 0, and R start from 1
	// (avoid index out of bound)
	if (coast_ID < 0)
		Rcpp::stop( " coast_ID < 0 " ) ;

	Point p = {P(0),P(1)}, u, v;
	int i, N = X.size();
	int k1, k2, K ;
	
	Rcpp::NumericVector d1(1);
	if (N==1) {
		u = {X(0),Y(0)};
		d1(0) = dist(p,u);
		return d1;
	}
	
	// Compute for only TWO coasts of ID "coast_ID"
	// (k1, k1+1) and (k1, k1-1)
	if ( coast_ID < MAX_INT_R ) {
		K = 1 ;
		k1 = coast_ID ;
		if ( coast_ID < (N-1) )
			k2 = k1+1 ;
		else // If it's the "last" point (largest order), connect it to the first point to make a coast, i.e. "circular" - coast_ID = N-1
			k2 = 0 ;
	}
	else { // Compute for all coast of the zone
		K = N ;
	}

	Rcpp::NumericVector d(K), e(K) ;
	double l2, a1, a2, a3 ;
	Point pu, pv, uv;
	
	// "Goto" can't skip over initializations of variables
	// https://stackoverflow.com/a/14274292
	if (K == N)
		goto WHOLE_SHORE_ZONE ;


	u = {X(k1),Y(k1)}; v = {X(k2),Y(k2)};
	pu.x = u.x-p.x; pu.y = u.y-p.y;
	pv.x = v.x-p.x; pv.y = v.y-p.y;
	uv.x = v.x-u.x; uv.y = v.y-u.y;
		
	a1 = Dot(uv,uv); a2 = Dot(pu,pu);	a3 = Dot(pv,pv);
	l2 = std::sqrt(a1);  // i.e. |u-v|_2
	if (l2 < 1e-14) // if u is too close to v
		d(0) = std::sqrt(a2); 
	else { // Collorary from Pythagorean theorem
		if (  ( a1+a3<a2 ) || ( a1+a2<a3 ) ) // altitude outside uv
			d(0) = Min(std::sqrt(a2),std::sqrt(a3));
		else // Shoelace formula
			d(0) = std::abs((pu.x)*(pv.y)-(pv.x)*(pu.y))/l2 ;
	}
	
	return d ;


  WHOLE_SHORE_ZONE:
	for (i=0; i<=(N-1); i++) {
		u = {X(i),Y(i)};
		if ( i < (N-1) )
			v = {X(i+1),Y(i+1)} ;
		else // last point of Zone (i.e. polygon)
			v = {X(0),Y(0)} ;
		
		pu.x = u.x-p.x; pu.y = u.y-p.y;
		pv.x = v.x-p.x; pv.y = v.y-p.y;
		uv.x = v.x-u.x; uv.y = v.y-u.y;
				
		a1 = Dot(uv,uv); a2 = Dot(pu,pu);	a3 = Dot(pv,pv);
		l2 = std::sqrt(a1);  // i.e. |u-v|_2
		if (l2 < 1e-14) // if u is too close to v
			d(i) = std::sqrt(a2); 
		else { // Collorary from Pythagorean theorem, or Law of Cosine
			if (  ( a1+a3<a2 ) || ( a1+a2<a3 ) ) // altitude outside uv
				d(i) = Min(std::sqrt(a2),std::sqrt(a3));
			else // Shoelace formula
				d(i) = std::abs((pu.x)*(pv.y)-(pv.x)*(pu.y))/l2 ;
		}
	}
	
	return d;
}


// Calculate the distances between whale positions to coast modeled by list of polygons (Shores)
// Limited to some polygons pre-selected by "nearest_zone_ID"
// Max_Speed is the assumed the maximum horizontal speed of whales
// [[Rcpp::export]]
Rcpp::DataFrame path_to_coast(
			Rcpp::DataFrame& path, Rcpp::List& Shores, 
			Rcpp::StringVector& nearest_zone_ID, double Max_Speed = 10.0  )
{
	using namespace Rcpp;
	std::string s ;
	int N = path.nrow() ;	
	int n = nearest_zone_ID.size() ;	
	int k, i, j;
		
	NumericVector whale_pos(2), beach_x, beach_y;
	DataFrame beach;
	NumericVector d(n), p(n), X, Y ;
	X = path["X"]; Y = path["Y"];
	NumericVector y(N), P(N), D(N)  ;
	StringVector ID(N) ;
	double d_min = 0.0 ; // nearest distance of last position
	int on_land_ID = -1 ;

	for (i = 0; i < N; i++) {
		whale_pos(0) = X(i);	whale_pos(1) = Y(i);	
		on_land_ID = -1 ;
		for (k = 0; k < n; k++) {
			// https://stackoverflow.com/a/8422324
			s = Rcpp::as<std::string>( nearest_zone_ID(k) ) ;
			beach = Shores[s] ;
			beach_x = beach["X"] ; beach_y = beach["Y"] ;
			// If inside this zone, it's on land so save this ID			
			if ( pnpoly_test(beach_x, beach_y, whale_pos(0), whale_pos(1)) ) {
				on_land_ID = k ;
			}
			// Very "raw" estimation from last nearest distance
			// Useful when there are big differences between distances to zones
			if ( d(k) < (d_min+Max_Speed) ) {
				y = distCoast_test(whale_pos, beach_x, beach_y) ;
				if ( on_land_ID < 0 )
					d(k) = min(y) ; 
				else // On land !!
					d(k) = 0 ; 
			}
			else { // If not, estimate a simple upper bound 
				d(k) -= Max_Speed ;
			}
		}
		
		if ( on_land_ID < 0 ) { // Not on land
			j = which_min(d);
			D(i) = min(d) ;
		}
		else { // On land
			j = on_land_ID ;
			D(i) = 0 ; 			
		}
		
		ID(i) = nearest_zone_ID(j) ;
		d_min = D(i) ;
	}
	
	return DataFrame::create( Named("DateTime") = path["DateTime"],
							  Named("Zone_ID")= ID,
							  Named("Distance")= D ) ;
}


// Test if whale position are on land (position), based on list of polygons (Shores) used to model the coast
// Limited to some polygons pre-selected by "nearest_zone_ID"
// [[Rcpp::export]]
Rcpp::DataFrame find_wrong_position(
					Rcpp::DataFrame& position, Rcpp::List& Shores, 
					Rcpp::StringVector& nearest_zone_ID  )
{
	using namespace Rcpp;

	int i, k, N, n;
	double minn ;
	N = position.nrow() ;	
	NumericVector d(N, 0.0), X, Y ;
	IntegerVector id(N, 0) ;
	X = position["X"]; Y = position["Y"];
	
	NumericVector whale_pos(2), beach_x, beach_y, beach_id;
	NumericVector y(1);
	DataFrame beach;
	
	n = nearest_zone_ID.size() ;
	Rcout << " n = " << n << "\n" ;
	for (k = 0; k < n ; k++) {
		beach = Shores[k] ;
		beach_x = beach["X"] ; beach_y = beach["Y"] ;
		beach_id = beach["Zone_ID"] ;
		
		for (i = 0; i < N; i++) {
			if ( pnpoly_test(beach_x, beach_y, X(i), Y(i)) ) {
				whale_pos(0) = X(i);	whale_pos(1) = Y(i);
				y = distCoast_test( whale_pos, beach_x, beach_y ) ;
				
				minn = -min(y) ;				
				if ( (d(i)==0.0) || (d(i) < minn) ) {
					d(i) = minn ; id(i) = beach_id(0) ;					
				}
			}
		}
	}
	
	return DataFrame::create( Named("d")= d, Named("ID") = id );
}

// Test if point (Point_X,Point_Y) is inside a polygon whose vertices' coordinates are X,Y
// [https://stackoverflow.com/a/2922778]
// [[Rcpp::export]]
bool pnpoly_test(Rcpp::NumericVector& X, Rcpp::NumericVector& Y, double Point_X, double Point_Y)
{
	int i, j;
	bool c = false;
	int nvert = X.size() ;
	
	for (i = 0, j = nvert-1; i < nvert; j = i++) {
		if ( ((Y[i]>Point_Y) != (Y[j]>Point_Y)) &&
			(Point_X < (X[j]-X[i]) * (Point_Y-Y[i]) / (Y[j]-Y[i]) + X[i]) )
			c = !c;
	}
	
	return c;
}
