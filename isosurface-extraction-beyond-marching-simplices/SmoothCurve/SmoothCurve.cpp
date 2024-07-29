#include <iostream>
#include <Eigen/Sparse>
#include "Misha/Miscellany.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Ply.h"
#include "Misha/PlyVertexData.h"

static const unsigned int Dim = 2;

Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" );
Misha::CmdLineParameter< double > DiffusionTime( "diffusion" , 1e-4 );
Misha::CmdLineReadable Verbose( "verbose" ) , Performance( "performance" );

Misha::CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&DiffusionTime ,
	&Verbose ,
	&Performance ,
	NULL
};

void ShowUsage( const char* ex )
{
	std::cout << "Usage " << std::string( ex ) << ":" << std::endl;
	std::cout << "\t --" << In.name << " <input curve>" << std::endl;
	std::cout << "\t[--" << Out.name << " <output curve>]" << std::endl;
	{
		Miscellany::StreamFloatPrecision sfw( std::cout , 1 , true );
		std::cout << "\t[--" << DiffusionTime.name << " <diffusion time>=" << DiffusionTime.value << "]" << std::endl;
	}
	std::cout << "\t[--" << Performance.name << "]" << std::endl;
	std::cout << "\t[--" << Verbose.name << "]" << std::endl;
}

int main( int argc , char* argv[] )
{
	using Factory = VertexFactory::PositionFactory< double , Dim >;

	Misha::CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_SUCCESS;
	}

	Miscellany::Timer timer , subTimer;

	// The input/output vertices
	std::vector< Factory::VertexType > vertices;
	// The input/output edges
	std::vector< SimplexIndex< Dim-1 > > edges;

	Factory vFactory;
	int file_type;
	PLY::ReadSimplices( In.value , vFactory , vertices , edges , NULL , file_type );

	std::cout << "Input vertices/edges: " << vertices.size() << " / " << edges.size() << std::endl;


	// To make the diffusion-time parameter consistent across meshes, we re-scale the geometry to have unit measure
	double scale = 0;
	for( unsigned int i=0 ; i<edges.size() ; i++ ) scale += sqrt( ( vertices[ edges[i][0] ] - vertices[ edges[i][1] ] ).squareNorm() );
	for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i] /= scale;

	subTimer.reset();
	Eigen::SparseMatrix< double > M( vertices.size() , vertices.size() ) , S( vertices.size() , vertices.size() );
	{
		std::vector< Eigen::Triplet< double > > M_triplets , S_triplets;

		for( unsigned int i=0 ; i<edges.size() ; i++ )
		{
			double len = sqrt( ( vertices[ edges[i][0] ] - vertices[ edges[i][1] ] ).squareNorm() );
			// Construct the mass matrix coefficients associated to this edge
			{
				M_triplets.emplace_back( edges[i][0] , edges[i][0] , len/3. );
				M_triplets.emplace_back( edges[i][1] , edges[i][1] , len/3. );
				M_triplets.emplace_back( edges[i][0] , edges[i][1] , len/6. );
				M_triplets.emplace_back( edges[i][1] , edges[i][0] , len/6. );
			}

			// Construct the stiffness matrix coefficients associated to this edge
			{
				S_triplets.emplace_back( edges[i][0] , edges[i][0] ,  1/len );
				S_triplets.emplace_back( edges[i][1] , edges[i][1] ,  1/len );
				S_triplets.emplace_back( edges[i][0] , edges[i][1] , -1/len );
				S_triplets.emplace_back( edges[i][1] , edges[i][0] , -1/len );
			}
		}

		M.setFromTriplets( M_triplets.begin() , M_triplets.end() );
		S.setFromTriplets( S_triplets.begin() , S_triplets.end() );
	}
	if( Verbose.set ) std::cout << "Set mass and stiffness: " << subTimer() << std::endl;

	subTimer.reset();
	Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > > solver( M + S*DiffusionTime.value );
	if( Verbose.set ) std::cout << "Factored system matrix: " << subTimer() << std::endl;

	subTimer.reset();
	Eigen::VectorXd b( vertices.size() );
	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		for( unsigned int i=0 ; i<vertices.size() ; i++ ) b[i] = vertices[i][d];
		b = solver.solve( M * b );
		for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i][d] = b[i];
	}
	if( Verbose.set ) std::cout << "Solved system: " << subTimer() << std::endl;


	for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i] *= scale;
	if( Out.set ) PLY::WriteSimplices( Out.value , vFactory , vertices , edges , file_type );

	if( Performance.set ) std::cout << "Performance: " << timer() << ", " << Miscellany::MemoryInfo::PeakMemoryUsageMB() << " (MB)" << std::endl;

	return EXIT_SUCCESS;
}
