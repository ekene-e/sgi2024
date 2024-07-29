#ifndef GRID_READER_INCLUDED
#define GRID_READER_INCLUDED

#include "Misha/RegularGrid.h"

template< unsigned int Dim , unsigned int ... Ns > struct GridReader;

template< unsigned int Dim >
struct GridReader< Dim >
{
	static RegularGrid< Dim , double > Read( std::string fileName , XForm< double , Dim+1 > &xForm )
	{
		std::string dataName;
		unsigned int dataDim;
		RegularGrid< Dim >::ReadHeader( fileName , dataDim , dataName );
		if( dataDim!=1 ) ERROR_OUT( "Only one-dimensional values per cell supported: " , dataDim );
		RegularGrid< Dim , double > grid;

		auto ReadAndConvertGrid = [&]< typename InType >( void )
		{
			XForm< InType , Dim+1 > _xForm;
			RegularGrid< Dim , InType > _grid;
			_grid.read( fileName , _xForm );
			grid.resize( _grid.res() );
			typename RegularGrid< Dim >::Range cornerRange;
			for( unsigned int d=0 ; d<Dim ; d++ ) cornerRange.first[d] = 0 , grid.res(d) , cornerRange.second[d] = grid.res(d)+1;

			cornerRange.process( [&]( typename RegularGrid< Dim >::Index I ){ grid(I) = (double)_grid(I); } );
			for( unsigned int i=0 ; i<=Dim ; i++ ) for( unsigned int j=0 ; j<=Dim ; j++ ) xForm(i,j) = (double)_xForm(i,j);
		};

		if     ( dataName==RegularGridDataType< double >::Name ) grid.read( fileName , xForm );
		else if( dataName==RegularGridDataType< float  >::Name ) ReadAndConvertGrid.template operator()< float >();
		else if( dataName==RegularGridDataType< int    >::Name ) ReadAndConvertGrid.template operator()< int   >();
		else ERROR_OUT( "Only float, double, and int type grids supported: " , dataName );
		return grid;
	}
};

template< unsigned int Dim , unsigned int N >
struct GridReader< Dim , N >
{
	static RegularGrid< Dim , Point< double , N > > Read( std::string fileName , XForm< double , Dim+1 > &xForm )
	{
		std::string dataName;
		unsigned int dataDim;
		RegularGrid< Dim >::ReadHeader( fileName , dataDim , dataName );
		if( dataDim!=N ) ERROR_OUT( "Only one-dimensional values per cell supported: " , dataDim );
		RegularGrid< Dim , Point< double , N > > grid;

		auto ReadAndConvertGrid = [&]< typename InType >( void )
		{
			XForm< InType , Dim+1 > _xForm;
			RegularGrid< Dim , Point< InType , N > > _grid;
			_grid.read( fileName , _xForm );
			grid.resize( _grid.res() );
			typename RegularGrid< Dim >::Range cornerRange;
			for( unsigned int d=0 ; d<Dim ; d++ ) cornerRange.first[d] = 0 , grid.res(d) , cornerRange.second[d] = grid.res(d)+1;

			cornerRange.process( [&]( typename RegularGrid< Dim >::Index I ){ grid(I) = Point< double , N >( _grid(I) ); } );
			for( unsigned int i=0 ; i<=Dim ; i++ ) for( unsigned int j=0 ; j<=Dim ; j++ ) xForm(i,j) = (double)_xForm(i,j);
		};

		if     ( dataName==RegularGridDataType< double >::Name ) grid.read( fileName , xForm );
		else if( dataName==RegularGridDataType< float  >::Name ) ReadAndConvertGrid.template operator()< float >();
		else if( dataName==RegularGridDataType< int    >::Name ) ReadAndConvertGrid.template operator()< int   >();
		else ERROR_OUT( "Only float, double, and int type grids supported: " , dataName );
		return grid;
	}
};

#endif // GRID_READER_INCLUDED