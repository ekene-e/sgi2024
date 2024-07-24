#include <map>
#include <iostream>
#include <random>
#include <type_traits>
#include "Misha/Miscellany.h"
#include "Misha/ProgressBar.h"
#include "Misha/CmdLineParser.h"
#include "Misha/RegularGrid.h"
#include "Misha/Ply.h"
#include "Misha/PlyVertexData.h"
#include "Include/MultiIndex.h"
#include "Include/GridReader.h"
#include "Include/CellSimplices.h"

static const unsigned int Dim = 2;

Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" );
Misha::CmdLineParameter< double > IsoValue( "iso" , 0. );
Misha::CmdLineReadable Verbose( "verbose" ) , Performance( "performance" ) , ASCII( "ascii" ) , Progress( "progress" );

Misha::CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&IsoValue ,
	&Verbose ,
	&Performance ,
	&Progress ,
	&ASCII ,
	NULL
};

void ShowUsage( const char* ex )
{
	std::cout << "Usage " << std::string( ex ) << ":" << std::endl;
	std::cout << "\t --" << In.name << " <input grid>" << std::endl;
	std::cout << "\t[--" << Out.name << " <output curve>]" << std::endl;
	std::cout << "\t[--" << IsoValue.name << " <iso-value>=" << IsoValue.value << "]" << std::endl;
	std::cout << "\t[--" << Performance.name << "]" << std::endl;
	std::cout << "\t[--" << Progress.name << "]" << std::endl;
	std::cout << "\t[--" << ASCII.name << "]" << std::endl;
	std::cout << "\t[--" << Verbose.name << "]" << std::endl;
}

int main( int argc , char *argv[] )
{
	using Factory = VertexFactory::PositionFactory< double , Dim >;

	// A function linearizing a grid's index
	auto Linearize = [&]( RegularGrid< Dim >::Index I , RegularGrid< Dim >::Range range )
		{
			unsigned int idx = I[0] - range.first[0];
			for( unsigned int d=1 ; d<Dim ; d++ ) idx = idx * ( range.second[d-1] - range.first[d-1] ) + ( I[d] - range.first[d] );
			return idx;
		};

	Misha::CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_SUCCESS;
	}

	Miscellany::Timer timer , subTimer;

	// The output level-set vertices
	std::vector< Factory::VertexType > levelSetVertices;
	// The output level-set edges
	std::vector< SimplexIndex< Dim-1 > > levelSetEdges;

	// The transformations from grid coordinates to world coordinates
	XForm< double , Dim+1 > gridToWorld;
	// The regular grid of input values
	RegularGrid< Dim , double > grid;
	// Range of grid cells and grid corners
	// [NOTE] Grid values are associated with corners
	RegularGrid< Dim >::Range cellRange , cornerRange;


	// Read in the input grid
	grid = GridReader< Dim >::Read( In.value , gridToWorld );

	for( unsigned int d=0 ; d<Dim ; d++ ) cellRange.first[d] = cornerRange.first[d] = 0 , cellRange.second[d] = grid.res(d)-1 , cornerRange.second[d] = grid.res(d);
	if( Verbose.set )
	{
		std::cout << "Grid resolution:";
		for( unsigned int d=0 ; d<Dim ; d++ ) std::cout << " " << grid.res(d);
		std::cout << std::endl;
		double min , max;
		min = max = grid[0];
		for( size_t i=0 ; i<grid.resolution() ; i++ ) min = std::min< double >( min , grid[i] ) , max = std::max< double >( max , grid[i] );
		std::cout << "Min/max: " << min << " / " << max << std::endl;
	}

	// An ordered map to track the level-set vertices associated with edges
	std::map< MultiIndex< Dim > , unsigned int > levelSetVertexMap;

	// Functionality for adding the level-set associated with a simplex
	auto AddLevelSetGeometry = [&]( SimplexIndex< Dim , RegularGrid< Dim >::Index > s )
		{
			// Given a triangle T = { (i_0,j_0) , (i_1,j_1) , (i_2,j_2) }
			double values[] = { grid( s[0] ) , grid( s[1] ) , grid( s[2] ) };

			unsigned int lCount=0 , gCount=0;
			for( unsigned int d=0 ; d<=Dim ; d++ )
				if     ( values[d]<IsoValue.value ) lCount++;
				else if( values[d]>IsoValue.value ) gCount++;

			if( lCount+gCount!=Dim+1 ) ERROR_OUT( "Not in general position" );
			if( lCount==0 || lCount==Dim+1 ) return;

			SimplexIndex< Dim-1 > edge;

			if( lCount==1 )
			{
				unsigned int ltIdx = -1;
				for( unsigned int d=0 ; d<=Dim ; d++ ) if( values[d]<IsoValue.value ) ltIdx = d;

				if( ltIdx==-1 ) ERROR_OUT( "Could not find less than vertex" );

				// The 2D index of the less-than corner is: ( s[ltIdx][0] , s[ltIdx][1] )
				// 20 - 21 - 22 - 23  - 24
				// 15 - 16 - 17 - 18  - 19
				// 10 - 11 - 12 - 13  - 14
				//  5 -  6 -  7 -  8  -  9
				//  0 -  1 -  2 -  3  -  4
				unsigned int ltLinearIndex = Linearize( s[ltIdx] , cornerRange );

				for( unsigned int i=0 ; i<2 ; i++ )
				{
					unsigned int gtIdx = ( ltIdx+1+i ) % ( Dim+1 );
					// The two values are values[ltIdx] and values[gtIdx]
					// alpha = values[ltIdx] * ( 1-t ) + values[gtIdx] * t
					// alpha = values[ltIdx] - values[ltIdx] * t + values[gtIdx] * t
					// alpha - values[ltIdx] = t * ( values[gtIdx] - values[ltIdx] )
					// ( alpha - values[ltIdx] ) / ( values[gtIdx] - values[ltIdx] ) = t
					double t = ( IsoValue.value - values[ltIdx] ) / ( values[gtIdx] - values[ltIdx] );
					Point< double , Dim > p;
					for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = s[ltIdx][d] * ( 1.-t ) + s[gtIdx][d] * t;

					unsigned int gtLinearIndex = Linearize( s[gtIdx] , cornerRange );

					MultiIndex< 2 > mi( ltLinearIndex , gtLinearIndex );
					unsigned int vIdx;

					if( levelSetVertexMap.find( mi )==levelSetVertexMap.end() )
					{
						vIdx = (unsigned int)levelSetVertices.size();
						levelSetVertices.push_back( p );
						levelSetVertexMap[ mi ] = vIdx;
					}
					else vIdx = levelSetVertexMap[ mi ];

					edge[i] = vIdx;
				}
				levelSetEdges.push_back( edge );
			}
			else if( lCount==2 )
			{
				unsigned int gtIdx = -1;
				for( unsigned int d=0 ; d<=Dim ; d++ ) if( values[d]>IsoValue.value ) gtIdx = d;

				if( gtIdx==-1 ) ERROR_OUT( "Could not find less than vertex" );
				unsigned int gtLinearIndex = Linearize( s[gtIdx] , cornerRange );

				for( unsigned int i=0 ; i<2 ; i++ )
				{
					unsigned int ltIdx = ( gtIdx+1+i ) % ( Dim+1 );
					// The two values are values[ltIdx] and values[gtIdx]
					// alpha = values[ltIdx] * ( 1-t ) + values[gtIdx] * t
					// alpha = values[ltIdx] - values[ltIdx] * t + values[gtIdx] * t
					// alpha - values[ltIdx] = t * ( values[gtIdx] - values[ltIdx] )
					// ( alpha - values[ltIdx] ) / ( values[gtIdx] - values[ltIdx] ) = t
					double t = ( IsoValue.value - values[ltIdx] ) / ( values[gtIdx] - values[ltIdx] );
					Point< double , Dim > p;
					for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = s[ltIdx][d] * ( 1.-t ) + s[gtIdx][d] * t;

					unsigned int ltLinearIndex = Linearize( s[ltIdx] , cornerRange );

					MultiIndex< 2 > mi( ltLinearIndex , gtLinearIndex );
					unsigned int vIdx;

					if( levelSetVertexMap.find( mi )==levelSetVertexMap.end() )
					{
						vIdx = (unsigned int)levelSetVertices.size();
						levelSetVertices.push_back( p );
						levelSetVertexMap[ mi ] = vIdx;
					}
					else vIdx = levelSetVertexMap[ mi ];

					edge[i] = vIdx;
				}
				std::swap< unsigned int >( edge[0] , edge[1] );
				levelSetEdges.push_back( edge );
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


	return EXIT_SUCCESS;
}
