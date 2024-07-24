#include <map>
#include <iostream>
#include <random>
#include <type_traits>
#include "Misha/Miscellany.h"
#include "Misha/CmdLineParser.h"
#include "Misha/RegularGrid.h"
#include "Misha/Ply.h"
#include "Misha/PlyVertexData.h"
#include "Include/GridReader.h"

Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" );
Misha::CmdLineParameter< unsigned int > SmoothingIterations( "iters" , 0 ) , Extract( "extract" , -1 );
Misha::CmdLineReadable Normalize( "normalize" ) , Discretize( "discretize" );

Misha::CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&SmoothingIterations ,
	&Extract ,
	&Normalize ,
	&Discretize ,
	NULL
};

void ShowUsage( const char* ex )
{
	std::cout << "Usage " << std::string( ex ) << ":" << std::endl;
	std::cout << "\t --" << In.name << " <input grid>" << std::endl;
	std::cout << "\t[--" << Out.name << " <output curve>]" << std::endl;
	std::cout << "\t[--" << SmoothingIterations.name << " <smoothing iterations>=" << SmoothingIterations.value << "]" << std::endl;
	std::cout << "\t[--" << Extract.name << " <extraction coordinate>]" << std::endl;
	std::cout << "\t[--" << Normalize.name << "]" << std::endl;
	std::cout << "\t[--" << Discretize.name << "]" << std::endl;
}

template< unsigned int Dim , unsigned int N >
RegularGrid< Dim , Point< double , N > > Smooth( const RegularGrid< Dim , Point< double , N > > &in )
{
	RegularGrid< Dim , Point< double , N > > out;
	out.resize( in.res() );
	typename RegularGrid< Dim >::Range range , nRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.first[d] = 0 , range.second[d] = in.res( d );

	const double Weights[] = { 0.25 , 0.5 , 0.25 };
	auto SetWeightedAverage = [&]( typename RegularGrid< Dim >::Index outI )
		{
			Point< double , N > sum = {};
			double wSum = 0;
			auto SumNeighborValues = [&]( typename RegularGrid< Dim >::Index inI )
				{
					double w = 1;
					for( unsigned int d=0 ; d<Dim ; d++ ) w *= Weights[ 1 + outI[d] - inI[d] ];
					sum += in( inI ) * w;
					wSum += w;
				};
			RegularGrid< Dim >::Range::Intersect( ( typename RegularGrid< Dim >::Range( outI ) ).dilate(1) , range ).process( SumNeighborValues );
			out( outI ) = sum / wSum;
		};
	range.process( SetWeightedAverage );
	return out;
}

template< unsigned int Dim , unsigned int N >
void Execute( void )
{
	XForm< double , Dim+1 > xForm;
	RegularGrid< Dim , Point< double , N > > grid = GridReader< Dim , N >::Read( In.value , xForm );
	if( Normalize.set )
	{
		if( N==1 ) WARN( "Normalizing one-dimensional data" );
		for( unsigned int i=0 ; i<grid.resolution() ; i++ )
		{
			double sum = 0;
			for( unsigned int n=0 ; n<N ; n++ ) sum += grid[i][n];
			for( unsigned int n=0 ; n<N ; n++ ) grid[i][n] /= sum;
		}
	}
	for( unsigned int i=0 ; i<SmoothingIterations.value ; i++ ) grid = Smooth( grid );
	if( Out.set )
	{
		if( Extract.set )
		{
			if( Extract.value>N ) ERROR_OUT( "Extraction coordinate out of bounds: 0 <= " , Extract.value , " < " , N );
			RegularGrid< Dim , double > _grid;
			_grid.resize( grid.res() );
			for( unsigned int i=0 ; i<grid.resolution() ; i++ ) _grid[i] = grid[i][ Extract.value ];
			_grid.write( Out.value , xForm );
		}
		else if( Discretize.set )
		{
			RegularGrid< Dim , unsigned int > _grid;
			_grid.resize( grid.res() );
			for( unsigned int i=0 ; i<grid.resolution() ; i++ )
			{
				unsigned int n=0;
				for( unsigned int _n=0 ; _n<N ; _n++ ) if( grid[i][_n]>grid[i][n] ) n = _n;
				_grid[i] = n;
			}
			XForm< unsigned int , Dim+1 > _xForm;
			for( unsigned int i=0 ; i<=Dim ; i++ ) for( unsigned int j=0 ; j<=Dim ; j++ ) _xForm(i,j) = (unsigned int)xForm(i,j);
			_grid.write( Out.value , _xForm );
		}
		else grid.write( Out.value , xForm );
	}
}

template< unsigned int Dim >
void Execute( void )
{
	unsigned int dataDim;
	std::string dataName;
	RegularGrid< Dim >::ReadHeader( In.value , dataDim , dataName );
	switch( dataDim )
	{
	case  1: Execute< Dim ,  1 >() ; break;
	case  2: Execute< Dim ,  2 >() ; break;
	case  3: Execute< Dim ,  3 >() ; break;
	case  4: Execute< Dim ,  4 >() ; break;
	case  5: Execute< Dim ,  5 >() ; break;
	case  6: Execute< Dim ,  6 >() ; break;
	case  7: Execute< Dim ,  7 >() ; break;
	case  8: Execute< Dim ,  8 >() ; break;
	case  9: Execute< Dim ,  9 >() ; break;
	case 10: Execute< Dim , 10 >() ; break;
	default: ERROR_OUT( "Only data dimensions 1..10 are supported: " , dataDim );
	}
}

int main( int argc , char *argv[] )
{
	Misha::CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_SUCCESS;
	}

	unsigned int dim;
	RegularGrid< 0 >::ReadDimension( In.value , dim );
	switch( dim )
	{
		case 1: Execute< 1 >() ; break;
		case 2: Execute< 2 >() ; break;
		case 3: Execute< 3 >() ; break;
		default: ERROR_OUT( "Only dimensions 1, 2, and 3 supported: " , dim );
	}
	return EXIT_SUCCESS;
}
