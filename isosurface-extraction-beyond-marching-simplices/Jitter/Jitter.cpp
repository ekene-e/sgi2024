#include <map>
#include <iostream>
#include <random>
#include <type_traits>
#include "Misha/Miscellany.h"
#include "Misha/CmdLineParser.h"
#include "Misha/RegularGrid.h"
#include "Include/GridReader.h"

static const unsigned int Dim = 2;

Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" );
Misha::CmdLineParameter< double > Jitter( "jitter" , 0. );

Misha::CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&Jitter ,
	NULL
};

void ShowUsage( const char* ex )
{
	std::cout << "Usage " << std::string( ex ) << ":" << std::endl;
	std::cout << "\t --" << In.name << " <input grid>" << std::endl;
	std::cout << "\t[--" << Out.name << " <output curve>]" << std::endl;
	std::cout << "\t[--" << Jitter.name << " <jitter magnitude>=" << Jitter.value << "]" << std::endl;
}

template< unsigned int N >
void Execute( void )
{
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

	// If the jitter magnitude is non-zero, jitter the input
	if( Jitter.value!=0 )
	{
		std::random_device rand_dev;
		std::mt19937 generator( rand_dev() );
		std::uniform_real_distribution< double > distr( -fabs(Jitter.value) , fabs(Jitter.value) );

		cornerRange.process( [&]( RegularGrid< Dim >::Index I ){ for( unsigned int n=0 ; n<N ; n++ ) grid( I )[n] += distr( generator ); } );

	}

	if( Out.set ) grid.write( Out.value , gridToWorld );
}


int main( int argc , char *argv[] )
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
	case  1: Execute<  1 >() ; break;
	case  2: Execute<  2 >() ; break;
	case  3: Execute<  3 >() ; break;
	case  4: Execute<  4 >() ; break;
	case  5: Execute<  5 >() ; break;
	case  6: Execute<  6 >() ; break;
	case  7: Execute<  7 >() ; break;
	case  8: Execute<  8 >() ; break;
	case  9: Execute<  9 >() ; break;
	case 10: Execute< 10 >() ; break;
	default: ERROR_OUT( "Only grid values of dimension 1..10 supported: " , dataDim );
	}

	return EXIT_SUCCESS;
}
