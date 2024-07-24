#ifndef CELL_SIMPLICES_INCLUDED
#define CELL_SIMPLICES_INCLUDED

#include "Misha/RegularGrid.h"
#include "Misha/Geometry.h"

// A structure representing the partition of a Dim-dimensional cube/cell into simplices
template< unsigned int Dim > struct CellSimplices;

template<>
struct CellSimplices< 2 >
{
	// The dimension of the cell
	static const unsigned int Dim = 2;
	// The number of simplices it partittions into
	static const unsigned int Num = 2;
	// The simplices (with corners giving the Dim-dimensional indices of the cube)
	SimplexIndex< Dim , RegularGrid< Dim >::Index > simplexIndices[Num];

	// Constructor initializing the simplices' indices, given the index of the cell within the grid
	CellSimplices( RegularGrid< Dim >::Index I )
	{
		simplexIndices[0][0] = I + Point< int , 2 >(0,0);
		simplexIndices[0][1] = I + Point< int , 2 >(1,0);
		simplexIndices[0][2] = I + Point< int , 2 >(0,1);

		simplexIndices[1][0] = I + Point< int , 2 >(1,1);
		simplexIndices[1][1] = I + Point< int , 2 >(0,1);
		simplexIndices[1][2] = I + Point< int , 2 >(1,0);
	}

	SimplexIndex< Dim , RegularGrid< Dim >::Index > &operator[]( unsigned int i ){ return simplexIndices[i]; }
	const SimplexIndex< Dim , RegularGrid< Dim >::Index > &operator[]( unsigned int i ) const { return simplexIndices[i]; }
};

#endif // CELL_SIMPLICES_INCLUDED