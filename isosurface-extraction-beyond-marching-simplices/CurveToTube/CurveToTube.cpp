#include <iostream>
#include "Misha/Miscellany.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Ply.h"
#include "Misha/PlyVertexData.h"

static const unsigned int Dim = 2;

Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" );
Misha::CmdLineParameter< unsigned int > AngularResolution( "res" , 8 );
Misha::CmdLineParameter< double > TubularRadius( "radius" , 0. );

Misha::CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&AngularResolution ,
	&TubularRadius ,
	NULL
};

void ShowUsage( const char* ex )
{
	std::cout << "Usage " << std::string( ex ) << ":" << std::endl;
	std::cout << "\t --" << In.name << " <input curve>" << std::endl;
	std::cout << "\t[--" << Out.name << " <output mesh>]" << std::endl;
	std::cout << "\t[--" << AngularResolution.name << " <angular resolution>=" << AngularResolution.value << "]" << std::endl;
	std::cout << "\t[--" << TubularRadius.name << " <tubular radius>=" << TubularRadius.value << "]" << std::endl;
}

template< unsigned int Dimension >
using Factory = VertexFactory::PositionFactory< double , Dimension >;

int main( int argc , char *argv[] )
{

	Misha::CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_SUCCESS;
	}
	if( AngularResolution.value<3 ) ERROR_OUT( "Angular resolution must be at least three: " , AngularResolution.value );

	// The input/output vertices
	std::vector< Factory< Dim >::VertexType > curveVertices;
	std::vector< Factory< Dim+1 >::VertexType > tubeVertices;
	// The input edges and output quads
	std::vector< SimplexIndex< Dim-1 > > curveEdges;
	std::vector< std::vector< unsigned int > > tubeQuads;

	{
		Factory< Dim > vFactory;
		int file_type;
		PLY::ReadSimplices( In.value , vFactory , curveVertices , curveEdges , NULL , file_type );
	}
	std::cout << "Input vertices/edges: " << curveVertices.size() << " / " << curveEdges.size() << std::endl;

	{
		std::vector< Point< double , Dim > > normals( curveVertices.size() );
		for( unsigned int i=0 ; i<curveEdges.size() ; i++ )
		{
			Point< double , Dim > dir = curveVertices[ curveEdges[i][1] ] - curveVertices[ curveEdges[i][0] ];

			Point< double , Dim > n = Point< double , Dim >::CrossProduct( dir );
			normals[ curveEdges[i][0] ] += n;
			normals[ curveEdges[i][1] ] += n;
		}
		for( unsigned int i=0 ; i<normals.size() ; i++ ) normals[i] /= sqrt( normals[i].squareNorm() );

		std::vector< unsigned int > quad( 4 );

		tubeVertices.resize( curveVertices.size() * AngularResolution.value );
		for( unsigned int i=0 ; i<curveVertices.size() ; i++ )
		{
			Point< double , 3 > v( curveVertices[i][0] , curveVertices[i][1] , 0 );
			Point< double , 3 > n1(0,0,1);
			Point< double , 3 > n2( normals[i][0] , normals[i][1] , 0 );
			for( unsigned int j=0 ; j<AngularResolution.value ; j++ )
			{
				double theta = ( 2. * M_PI * j ) / AngularResolution.value;
				tubeVertices[ i*AngularResolution.value + j ] = v + n1 * sin( theta ) * TubularRadius.value + n2 * cos( theta ) * TubularRadius.value;
			}
		}

		for( unsigned int i=0 ; i<curveEdges.size() ; i++ )
		{
			unsigned int i1 = curveEdges[i][0] , i2 = curveEdges[i][1];

			for( unsigned int j=0 ; j<AngularResolution.value ; j++ )
			{
				quad[0] = i1 * AngularResolution.value + ( ( j + 0 ) % AngularResolution.value );
				quad[1] = i2 * AngularResolution.value + ( ( j + 0 ) % AngularResolution.value );
				quad[2] = i2 * AngularResolution.value + ( ( j + 1 ) % AngularResolution.value );
				quad[3] = i1 * AngularResolution.value + ( ( j + 1 ) % AngularResolution.value );

				tubeQuads.push_back( quad );
			}
		}
	}

	std::cout << "Output vertices/quads: " << tubeVertices.size() << " / " << tubeQuads.size() << std::endl;

	if( Out.set )
	{
		Factory< Dim+1 > vertexFactory;
		PLY::WritePolygons( Out.value , vertexFactory , tubeVertices , tubeQuads , PLY_ASCII );
	}

	return EXIT_SUCCESS;
}
