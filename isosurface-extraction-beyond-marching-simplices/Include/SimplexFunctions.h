#ifndef SIMPLEX_FUNCTIONS_INCLUDED
#define SIMPLEX_FUNCTIONS_INCLUDED

#include "Misha/Geometry.h"

// A class for representing linear functions on a right Dim-dimensinonal simplex
template< unsigned int Dim >
struct SimplexFunction
{
	// Set the linear function's coefficients from its values
	SimplexFunction( const double cornerValues[Dim+1] );

	// Set the linear function's coefficients from its values
	template< typename ... Doubles >
	SimplexFunction( const Doubles ... cornerValues );

	// Evaluate the linear function at a prescribed point
	double operator()( Point< double , Dim > p ) const;

	// Returns the position of the a corner on the right simplex
	static Point< double , Dim > Corner( unsigned int c );

	// Returns the position at which the (Dim+1) linear functions have the same value
	static Point< double , Dim > Intersect( const SimplexFunction< Dim > functions[Dim+1] );

	// Returns the position at which the (Dim+1) linear functions have the same value
	template< typename ... SimplexFunctions >
	static Point< double , Dim > Intersect( SimplexFunctions ... functions );

	// Stream insertion operator
	template< unsigned int _Dim >
	friend std::ostream &operator << ( std::ostream &os , const SimplexFunction< _Dim > &f );

protected:
	// Represent the function as F(p) =  c + < v , p >
	double _c;
	Point< double , Dim > _v;

	void _set( const double cornerValues[Dim+1] );
};

/////////////////
// Definitions //
/////////////////

template< unsigned int Dim >
SimplexFunction< Dim >::SimplexFunction( const double cornerValues[Dim+1] ) { _set( cornerValues ); }

template< unsigned int Dim >
template< typename ... Doubles >
SimplexFunction< Dim >::SimplexFunction( const Doubles ... cornerValues )
{
	static_assert( sizeof...( Doubles )==Dim+1 , "[ERROR] Wrong number of corner values" );
	double cValues[] = { cornerValues... };
	_set( cValues );
}

template< unsigned int Dim >
double SimplexFunction< Dim >::operator()( Point< double , Dim > p ) const { return _c + Point< double , Dim >::Dot( _v , p ); }

template< unsigned int Dim >
Point< double , Dim > SimplexFunction< Dim >::Corner( unsigned int c )
{
	if( c>Dim ) ERROR_OUT( "Corner out of range: 0 <= " , c , " <= " , Dim );
	Point< double , Dim > p;
	if( c ) p[c-1] = 1;
	return p;
}

template< unsigned int Dim >
Point< double , Dim > SimplexFunction< Dim >::Intersect( const SimplexFunction< Dim > functions[Dim+1] )
{
	// The functions are equal at the point p satisfying:
	// functions[0](p) = functions[d](p) for all 1 <= d <= Dim
	// _c[0] + < _v[0] , p > = _c[d] + < _v[d] , p >
	// _c[0] - _c[d] = < _v[d] - _v[0] , p >
	SquareMatrix< double , Dim > A;
	Point< double , Dim > b;
	for( unsigned int d=1 ; d<=Dim ; d++ )
	{
		b[d-1] = functions[0]._c - functions[d]._c;
		Point< double , Dim > v = functions[d]._v - functions[0]._v;
		for( unsigned int dd=0 ; dd<Dim ; dd++ ) A( dd , d-1 ) = v[dd];
	}
	bool success;
	b = A.inverse( success ) * b;
	if( !success ) THROW( "No intersection found" );
	return b;
}

template< unsigned int Dim >
template< typename ... SimplexFunctions >
Point< double , Dim > SimplexFunction< Dim >::Intersect( SimplexFunctions ... functions )
{
	static_assert( sizeof...( SimplexFunctions )==Dim+1 , "[ERROR] Wrong number of simplex functions" );
	const SimplexFunction _functions[] = { functions... };
	return Intersect( _functions );
}

template< unsigned int Dim >
void SimplexFunction< Dim >::_set( const double cornerValues[Dim+1] )
{
	_c = cornerValues[0];
	for( unsigned int d=0 ; d<Dim ; d++ ) _v[d] = cornerValues[d+1] - cornerValues[0];
}

template< unsigned int Dim >
std::ostream &operator << ( std::ostream &os , const SimplexFunction< Dim > &f )
{
	os << " {";
	for( unsigned int d=0 ; d<=Dim ; d++ )
	{
		if( d ) os << " ,";
		os << " " << f( SimplexFunction< Dim >::Corner(d) );
	}
	return os << " }";
}
#endif // SIMPLEX_FUNCTIONS_INCLUDED