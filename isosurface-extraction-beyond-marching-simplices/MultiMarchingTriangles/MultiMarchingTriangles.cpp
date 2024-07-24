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

static const unsigned int Dim = 2;


Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" );
Misha::CmdLineReadable Verbose( "verbose" ) , Progress( "progress" ) , Performance( "performance" ) , ASCII( "ascii" ) , Progess( "progress" );
Misha::CmdLineReadable* params[] =
{
	&In ,
	&Out ,
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
	printf( "\t[--%s]\n" , Progress.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
	printf( "\t[--%s]\n" , Progress.name.c_str() );
	printf( "\t[--%s]\n" , Performance.name.c_str() );
	printf( "\t[--%s]\n" , ASCII.name.c_str() );
}

template< unsigned int N >
void Process( void )
{
	using Factory = VertexFactory::PositionFactory< double , Dim >;
	using TriangleVertexMapData =  std::vector< std::pair< MultiIndex< Dim+1 > , unsigned int > >;
	using EdgeVertexMapData = std::vector< std::pair< MultiIndex< Dim > , unsigned int > >;
	using TriangleVertexMap = std::map< MultiIndex< Dim+1 > , TriangleVertexMapData >;
	using EdgeVertexMap = std::map< MultiIndex< Dim > , EdgeVertexMapData >;

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

			///////////////////////
			// Fill in code here //
			// Iterate over all pairs of functions
			// .. Compute the point at which all the functions are equal
			// .. Check that the point is inside the unit-interval
			// .. If it is, add the point (in world coordinates) along with the indices of the two functions generating it to "vertices"
			///////////////////////
			for (unsigned int i = 0; i < N; i++)
				{
					for (unsigned int j = i + 1; j < N; j++)
					{
						double t = (grid(e[0])[i] - grid(e[1])[i]) / ((grid(e[0])[i] - grid(e[1])[i]) + (grid(e[1])[j] - grid(e[0])[j]));
						if (t >= 0 && t <= 1)
						{
							Point<double, Dim> p;
							for (unsigned int d = 0; d < Dim; d++) p[d] = e[0][d] * (1. - t) + e[1][d] * t;
							vertices.emplace_back(mi, (unsigned int)levelSetVertices.size());
							levelSetVertices.push_back(p);
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

			///////////////////////
			// Fill in code here //
			// Iterate over all triplets of functions
			// .. Compute the point at which all the functions are equal
			// .. Check that the point is inside the unit-right-triangle
			// .. If it is, add the point (in world coordinates) along with the indices of the three functions generating it to "vertices"
			///////////////////////
			for (unsigned int i = 0; i < N; i++)
			{
				for (unsigned int j = i + 1; j < N; j++)
				{
					for (unsigned int k = j + 1; k < N; k++)
					{
						Eigen::Matrix<double, 2, 2> A;
						A << grid(t[1])[i] - grid(t[0])[i], grid(t[2])[i] - grid(t[0])[i],
							 grid(t[1])[j] - grid(t[0])[j], grid(t[2])[j] - grid(t[0])[j];
						Eigen::Vector2d b(grid(t[0])[i] - grid(t[0])[j], grid(t[0])[k] - grid(t[0])[i]);

						Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);

						if (x[0] >= 0 && x[1] >= 0 && x[0] + x[1] <= 1)
						{
							Point<double, Dim> p;
							for (unsigned int d = 0; d < Dim; d++) p[d] = t[0][d] * (1. - x[0] - x[1]) + t[1][d] * x[0] + t[2][d] * x[1];
							vertices.emplace_back(mi, (unsigned int)levelSetVertices.size());
							levelSetVertices.push_back(p);
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
			// If there is one function which dominates all the others at every vertex, there will not be a level-set passing through
			{
				unsigned int dominatingIndex = -1;
				for( unsigned int n=0 ; n<N && dominatingIndex==-1; n++ )
				{
					bool dominates = true;
					for( unsigned int _n=0 ; _n<N && dominates ; _n++ ) if( _n!=n )
						for( unsigned int d=0 ; d<=Dim ; d++ ) if( grid( s[d] )[_n]>grid( s[d] )[n] )
						{
							dominates = false;
							break;
						}
					if( dominates ) dominatingIndex = n;
				}
				if( dominatingIndex!=-1 ) return;
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

			///////////////////////
			// Fill in code here //
			// Given the iterators for the vertices generated along the triangle's edges and in the triangle's interior
			// .. Iterate over all pairs of functions
			// .. Connect the vertices labeled with the two function's indices into an edge
			// .. Add the edge "levelSetEdges"
			///////////////////////
			for (unsigned int i = 0; i < N; i++)
			{
				for (unsigned int j = i + 1; j < N; j++)
				{
					std::vector<std::pair<MultiIndex<Dim+1>, unsigned int>> combinedVertices;

					for (auto &edge : edgeVertices)
					{
						for (auto &vertex : edge->second)
						{
							if (vertex.first.contains(i) && vertex.first.contains(j))
								combinedVertices.push_back(vertex);
						}
					}

					for (auto &vertex : triangleVertices->second)
					{
						if (vertex.first.contains(i) && vertex.first.contains(j))
							combinedVertices.push_back(vertex);
					}

					for (unsigned int k = 0; k < combinedVertices.size(); k++)
					{
						for (unsigned int l = k + 1; l < combinedVertices.size(); l++)
						{
							SimplexIndex<Dim-1> edge;
							edge[0] = combinedVertices[k].second;
							edge[1] = combinedVertices[l].second;
							levelSetEdges.push_back(edge);
						}
					}
				}
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
