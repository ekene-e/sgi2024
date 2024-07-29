#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <array>
#include <iostream>
#include <random>
#include <type_traits>
#include "Misha/Miscellany.h"
#include "Misha/ProgressBar.h"
#include "Misha/CmdLineParser.h"
#include "Misha/RegularGrid.h"
#include "Misha/Ply.h"
#include "Misha/PlyVertexData.h"
#include "Include/GridReader.h"
#include "Include/MultiIndex.h"
#include "Include/CellSimplices.h"
#include "Include/SimplexFunctions.h"
#include "Include/ConvexHull.h"

static const unsigned int Dim = 2;

#define USE_CONVEX_HULL


Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" );
Misha::CmdLineReadable Verbose( "verbose" ) , Progress( "progress" ) , Performance( "performance" ) , ASCII( "ascii" ) , Progess( "progress" ) , NoCulling( "noCulling" ) , NoConvexHull( "noHull" );
Misha::CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&NoCulling ,
	&NoConvexHull ,
	&Progress ,
	&Verbose ,
	&Performance ,
	&Progress ,
	&ASCII ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input grid>\n" , In.name.c_str() );
	printf( "\t[--%s <output curve>]\n" , Out.name.c_str() );
	printf( "\t[--%s]\n" , NoCulling.name.c_str() );
	printf( "\t[--%s]\n" , NoConvexHull.name.c_str() );
	printf( "\t[--%s]\n" , Progress.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
	printf( "\t[--%s]\n" , Performance.name.c_str() );
	printf( "\t[--%s]\n" , ASCII.name.c_str() );
}

template< unsigned int N >
void Process( void )
{
	using Factory = VertexFactory::PositionFactory< double , Dim >;
	using TriangleVertexData =  std::vector< std::pair< MultiIndex< Dim+1 > , unsigned int > >;
	using EdgeVertexData = std::vector< std::pair< MultiIndex< Dim > , unsigned int > >;
	using TriangleVertexMap = std::map< MultiIndex< Dim+1 > , TriangleVertexData >;
	using EdgeVertexMap = std::map< MultiIndex< Dim > , EdgeVertexData >;

	// A function linearizing a grid's index
	auto Linearize = [&]( RegularGrid< Dim >::Index I , RegularGrid< Dim >::Range range )
		{
			unsigned int idx = I[0] - range.first[0];
			for( unsigned int d=1 ; d<Dim ; d++ ) idx = idx * ( range.second[d-1] - range.first[d-1] ) + ( I[d] - range.first[d] );
			return idx;
		};

	Miscellany::Timer timer , subTimer;

	// The output level-set vertices
	std::vector< Factory::VertexType > levelSetVertices;
	// The output level-set edges
	std::vector< SimplexIndex< Dim-1 > > levelSetEdges;
	// An ordered map to track the level-set vertices associated with edges
	EdgeVertexMap edgeVertexMap;
	// An ordered map to track the level-set vertices associated with triangles
	TriangleVertexMap triangleVertexMap;


	// The transformations from grid coordinates to world coordinates
	XForm< double , Dim+1 > gridToWorld;
	// The regular grid of input values
	RegularGrid< Dim , Point< double , N > > grid;
	// Range of grid cells and grid corners
	// [NOTE] Grid values are associated with corners
	RegularGrid< Dim >::Range cellRange , cornerRange;


	// Read in the input grid
	grid = GridReader< Dim , N >::Read( In.value , gridToWorld );
	for( unsigned int d=0 ; d<Dim ; d++ ) cellRange.first[d] = cornerRange.first[d] = 0 , cellRange.second[d] = grid.res(d)-1 , cornerRange.second[d] = grid.res(d);
	if( Verbose.set )
	{
		std::cout << "Grid resolution:";
		for( unsigned int d=0 ; d<Dim ; d++ ) std::cout << " " << grid.res(d);
		std::cout << std::endl;
	}

	// Normalize the values to be weights in the range [0,1]
	{
		subTimer.reset();
		auto Normalize = []( Point< double , N > & v )
			{
				double sum = 0;
				for( unsigned int n=0 ; n<N ; n++ ) if( v[n]>0 ) sum += v[n];
				if( !sum ) ERROR_OUT( "Could not normalize value: " , v );
				for( unsigned int n=0 ; n<N ; n++ )
					if( v[n]<0 ) v[n] = 0;
					else v[n] /= sum;
			};
		cornerRange.process( [&]( RegularGrid< Dim >::Index I ){ Normalize( grid(I) ); } );

		if( Verbose.set ) std::cout << "Normalized: " << subTimer() << std::endl;
	}

	// Functionality for adding the level-set vertices associated with an edge

	auto AddEdgeVertices = [&]( SimplexIndex< Dim-1 , RegularGrid< Dim >::Index > e )
		{
			MultiIndex< Dim > mi( Linearize( e[0] , cornerRange ) , Linearize( e[1] , cornerRange ) );

			// Check if the edge's vertices have already been computed
			if( edgeVertexMap.find( mi )!=edgeVertexMap.end() ) return mi;

			// If they have not already been added, add them now
			std::vector< std::pair< MultiIndex< Dim > , unsigned int > > vertices;

			// Fit functions to the corner values
			SimplexFunction< Dim-1 > f[N];
			for( unsigned int n=0 ; n<N ; n++ ) f[n] = SimplexFunction< Dim-1 >( grid( e[0] )[n] , grid( e[1] )[n] );

			if( !NoConvexHull.set )
			{
				std::vector< Point< double , Dim > > duals( N );
				for( unsigned int n=0 ; n<N ; n++ ) duals[n] = f[n].dual();

				std::vector< SimplexIndex< Dim-1 > > hull = ConvexHull::ConvexHull( duals , false );
				for( unsigned int i=0 ; i<hull.size() ; i++ )
				{
					SimplexIndex< Dim-1 > si = hull[i];
					Simplex< double , Dim , Dim-1 > s;
					for( unsigned int d=0 ; d<Dim ; d++ ) s[d] = duals[ si[d] ];
					if( s.normal()[0]<0 )
					{
						// Find the point of intersection of the three functions, dual to the corners of a hull triangle
						Point< double , Dim-1 > x;
						try{ x = SimplexFunction< Dim-1 >::Intersect( f[ si[0] ] , f[ si[1] ] ); }
						catch( Misha::Exception ){ ERROR_OUT( "Expected intersection" ); }

						// Check that the position is on the edge
						if( x[0]>=0 && x[0]<=1 )
						{
							// Compute the world coordinates of the position
							Point< double , 2 > p = Point< double , 2 >( e[0] ) * ( 1. - x[0] ) + Point< double , 2 >( e[1] ) * x[0] ;

							// Add to the triangle-to-vertex-index map, and add to the list of vertices
							vertices.push_back( std::pair< MultiIndex< Dim > , unsigned int >( MultiIndex< Dim >( si[0] , si[1] ) , (unsigned int)levelSetVertices.size() ) );
							levelSetVertices.push_back( p );
						}
					}
				}
			}
			else
			{
				// For every pair of functions, find the point where the functions are equal, check if that is the maximal value, and add the point if it is
				for( unsigned int i=0 ; i<N ; i++ ) for( unsigned int j=0 ; j<i ; j++ )
				{
					Point< double , 1 > x;
					bool foundIntersection = true;
					try{ x = SimplexFunction< 1 >::Intersect( f[i] , f[j] ); }
					catch( Misha::Exception ){ foundIntersection = false; }
					if( !foundIntersection ) continue;

					// Check that the position is on the edge
					if( x[0]>=0 && x[0]<=1 )
					{
						// Check that the value is maximized by the pair (i,j)
						bool isMax = true;
						for( unsigned int k=0 ; k<N ; k++ ) if( k!=i && k!=j ) if( f[k](x)>f[i](x) ) isMax = false;
						if( isMax )
						{
							// Compute the world coordinates of the position
							Point< double , 2 > p = Point< double , 2 >( e[0] ) * ( 1. - x[0] ) + Point< double , 2 >( e[1] ) * x[0];

							// Add to the edge-to-vertex-index map, and add to the list of vertices
							vertices.push_back( std::pair< MultiIndex< Dim > , unsigned int >( MultiIndex< Dim >(i,j) , (unsigned int)levelSetVertices.size() ) );
							levelSetVertices.push_back( p );
						}
					}
				}
			}

			edgeVertexMap[ mi ] = vertices;
			return mi;
		};

	// Functionality for adding the level-set vertices associated with a triangle
	auto AddTriangleVertices = [&]( SimplexIndex< Dim , RegularGrid< Dim >::Index > t )
		{
			MultiIndex< Dim+1 > mi( Linearize( t[0] , cornerRange ) , Linearize( t[1] , cornerRange ) , Linearize( t[2] , cornerRange ) );
			// Check if the triangle's vertices have already been computed
			if( triangleVertexMap.find( mi )!=triangleVertexMap.end() ) return mi;

			// If they have not already been added, add them now
			std::vector< std::pair< MultiIndex< Dim+1 > , unsigned int > > vertices;

			// Fit functions to the corner values
			SimplexFunction< Dim > f[N];
			for( unsigned int n=0 ; n<N ; n++ ) f[n] = SimplexFunction< Dim >( grid( t[0] )[n] , grid( t[1] )[n] , grid( t[2] )[n] );

			if( !NoConvexHull.set )
			{
				std::vector< Point< double , Dim+1 > > duals( N );
				for( unsigned int n=0 ; n<N ; n++ ) duals[n] = f[n].dual();

				std::vector< SimplexIndex< Dim > > hull = ConvexHull::ConvexHull( duals , false );
				for( unsigned int i=0 ; i<hull.size() ; i++ )
				{
					SimplexIndex< Dim > si = hull[i];
					Simplex< double , Dim+1 , Dim > s;
					for( unsigned int d=0 ; d<=Dim ; d++ ) s[d] = duals[ si[d] ];
					if( s.normal()[0]<0 )
					{
						// Find the point of intersection of the three functions, dual to the corners of a hull triangle
						Point< double , Dim > xy;
						try{ xy = SimplexFunction< Dim >::Intersect( f[ si[0] ] , f[ si[1] ] , f[ si[2] ] ); }
						catch( Misha::Exception ){ ERROR_OUT( "Expected intersection" ); }

						// Check that the position is on the triangle
						if( xy[0]>=0 && xy[0]<=1 && xy[1]>=0 && xy[1]<=1 && (xy[0] + xy[1])<=1 )
						{
							// Compute the world coordinates of the position
							Point< double , 2 > p = Point< double , 2 >( t[0] ) * ( 1. - xy[0] - xy[1] ) + Point< double , 2 >( t[1] ) * xy[0] + Point< double , 2 >( t[2] ) * xy[1];

							// Add to the triangle-to-vertex-index map, and add to the list of vertices
							vertices.push_back( std::pair< MultiIndex< Dim+1 > , unsigned int >( MultiIndex< Dim+1 >( si[0] , si[1] , si[2] ) , (unsigned int)levelSetVertices.size() ) );
							levelSetVertices.push_back( p );
						}
					}
				}
			}
			else
			{
				// For every triplet of functions, find the point where the functions are equal, check if that is the maximal value, and add the point if it is
				for( unsigned int i=0 ; i<N ; i++ ) for( unsigned int j=0 ; j<i ; j++ ) for( unsigned int k=0 ; k<j ; k++ )
				{
					Point< double , 2 > xy;
					bool foundIntersection = true;
					try{ xy = SimplexFunction< 2 >::Intersect( f[i] , f[j] , f[k] ); }
					catch( Misha::Exception ){ foundIntersection = false; }
					if( !foundIntersection ) continue;

					// Check that the position is on the triangle
					if( xy[0]>=0 && xy[0]<=1 && xy[1]>=0 && xy[1]<=1 && (xy[0] + xy[1])<=1 )
					{
						// Check that the value is maximized by the triplet (i,j,k)
						bool isMax = true;
						for( unsigned int l=0 ; l<N ; l++ ) if( l!=i && l!=j && l!=k ) if( f[l](xy)>f[i](xy) ) isMax = false;
						if( isMax )
						{
							// Compute the world coordinates of the position
							Point< double , 2 > p = Point< double , 2 >( t[0] ) * ( 1. - xy[0] - xy[1] ) + Point< double , 2 >( t[1] ) * xy[0] + Point< double , 2 >( t[2] ) * xy[1];

							// Add to the triangle-to-vertex-index map, and add to the list of vertices
							vertices.push_back( std::pair< MultiIndex< Dim+1 > , unsigned int >( MultiIndex< Dim+1 >(i,j,k) , (unsigned int)levelSetVertices.size() ) );
							levelSetVertices.push_back( p );
						}
					}
				}
			}

			triangleVertexMap[ mi ] = vertices;
			return mi;
		};

	// Functionality for adding the level-set associated with a simplex
	auto AddLevelSetGeometry = [&]( SimplexIndex< Dim , RegularGrid< Dim >::Index > s )
		{
			if( !NoCulling.set )
			{
				bool hasDominatingLabel = false;
				// Try if the i-th label dominates all other functions at all the corners
				for( unsigned int i=0 ; i<N ; i++ )
				{
					bool isDominant = true;

					// Try all other functions
					for( unsigned int j=0 ; j<N ; j++ ) if( j!=i )
						// At all other corners
						for( unsigned int d=0 ; d<=Dim ; d++ )
							// If the j-th function is larger at any corner, the i-th function cannot dominate
							if( grid( s[d] )[j] > grid( s[d] )[i] ) isDominant = false;
					if( isDominant ) hasDominatingLabel = true;
				}
				if( hasDominatingLabel ) return;
			}

			//  Add multi-level-set vertices along the edged and in the interior of the triangle
			TriangleVertexMap::iterator triangleVertices;
			EdgeVertexMap::iterator edgeVertices[Dim+1];

			triangleVertices = triangleVertexMap.find( AddTriangleVertices( s ) );
			for( unsigned int d=0 ; d<=Dim ; d++ )
			{
				SimplexIndex< Dim-1 , RegularGrid< Dim >::Index > e;
				e[0] = s[(d+1)%(Dim+1)] , e[1] = s[(d+2)%(Dim+1)];
				edgeVertices[d] = edgeVertexMap.find( AddEdgeVertices( e ) );
			}

			for( unsigned int i=0 ; i<N ; i++ ) for( unsigned int j=0 ; j<i ; j++ )
			{
				SimplexIndex< Dim-1 > levelSetEdge;
				levelSetEdge[0] = levelSetEdge[1] = -1;

				// Look at the vertices generated inside the triangle
				{
					const std::vector< std::pair< MultiIndex< Dim+1 > , unsigned int > > &vertices = triangleVertices->second;
					for( unsigned int v=0 ; v<vertices.size() ; v++ )
						if( ( vertices[v].first[0]==i || vertices[v].first[1]==i || vertices[v].first[2]==i ) && ( vertices[v].first[0]==j || vertices[v].first[1]==j || vertices[v].first[2]==j ) )
						{
							if     ( levelSetEdge[0]==-1 ) levelSetEdge[0] = vertices[v].second;
							else if( levelSetEdge[1]==-1 ) levelSetEdge[1] = vertices[v].second;
							else ERROR_OUT( "Edge is full" );
						}
				}

				// Look at the vertices generated inside the edges
				for( unsigned int d=0 ; d<=Dim ; d++ )
				{
					const std::vector< std::pair< MultiIndex< Dim > , unsigned int > > &vertices = edgeVertices[d]->second;
					for( unsigned int v=0 ; v<vertices.size() ; v++ )
						if( ( vertices[v].first[0]==i || vertices[v].first[1]==i ) && ( vertices[v].first[0]==j || vertices[v].first[1]==j ) )
						{
							if     ( levelSetEdge[0]==-1 ) levelSetEdge[0] = vertices[v].second;
							else if( levelSetEdge[1]==-1 ) levelSetEdge[1] = vertices[v].second;
							else ERROR_OUT( "Edge is full" );
						}
				}

				if     ( levelSetEdge[1]!=-1 ) levelSetEdges.push_back( levelSetEdge );
				else if( levelSetEdge[0]!=-1 ) ERROR_OUT( "Could not complete edge" );
			}
		};

	// Functionality for adding the level-set associated with a cell
	char progressText[1024];
	ProgressBar progressBar( 10 , cellRange.second[0] - cellRange.first[0] , progressText , false );
	auto GetCellLevelSet = [&]( RegularGrid< Dim >::Index I )
		{
			if( Progress.set )
			{
				bool show = true;
				for( unsigned int d=1 ; d<Dim ; d++ ) if( I[d] ) show = false;
				if( show )
				{
					sprintf( progressText , "Processing cells" );
					progressBar.update();
				}
			}

			CellSimplices< Dim > cellSimplices( I );
			AddLevelSetGeometry( cellSimplices[0] );
			AddLevelSetGeometry( cellSimplices[1] );
		};
	if( Progress.set ) std::cout << std::endl;
	// Iterate over the cells and add the level sets
	subTimer.reset();
	cellRange.process( GetCellLevelSet );

	// Transform vertices into world coordinates
	for( unsigned int i=0 ; i<levelSetVertices.size() ; i++ ) levelSetVertices[i] = gridToWorld * levelSetVertices[i];
	if( Verbose.set )
	{
		std::cout << "Got level-set: " << subTimer() << std::endl;
		std::cout << "Vertices/edges: " << levelSetVertices.size() << " / " << levelSetEdges.size() << std::endl;
	}

	if( Out.set )
	{
		Factory vertexFactory;
		PLY::WriteSimplices( Out.value , vertexFactory , levelSetVertices , levelSetEdges , ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE );
	}

	if( Performance.set ) std::cout << "Performance: " << timer() << ", " << Miscellany::MemoryInfo::PeakMemoryUsageMB() << " (MB)" << std::endl;
}

int main( int argc , char* argv[] )
{
	Misha::CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_SUCCESS;
	}

	unsigned int dataDim;
	std::string dataName;
	RegularGrid< Dim >::ReadHeader( In.value , dataDim , dataName );

	switch( dataDim )
	{
	case  2: Process<  2 >() ; break;
	case  3: Process<  3 >() ; break;
	case  4: Process<  4 >() ; break;
	case  5: Process<  5 >() ; break;
	case  6: Process<  6 >() ; break;
	case  7: Process<  7 >() ; break;
	case  8: Process<  8 >() ; break;
	case  9: Process<  9 >() ; break;
	case 10: Process< 10 >() ; break;
	default: ERROR_OUT( "Only grid values of dimension 2..10 supported: " , dataDim );
	}

	return EXIT_SUCCESS;
}
