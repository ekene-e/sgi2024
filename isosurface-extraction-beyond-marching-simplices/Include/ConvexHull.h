#ifndef CONVEX_HULL_INCLUDED
#define CONVEX_HULL_INCLUDED

#define NEW_CONVEX_HULL
#define FAST_SIMPLE_INCREMENTAL	// For the case where the number of points is Dim+2
#undef USE_CPP_QHULL
#ifdef USE_CPP_QHULL
#define DISABLE_QHULL_TRY_CATCH
#endif // USE_CPP_QHULL

#ifdef USE_CPP_QHULL
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullVertexSet.h"
#else // !USE_CPP_QHULL
extern "C"
{
	#include "libqhull/libqhull.h"
}
#endif // USE_CPP_QHULL
#include "Misha/Geometry.h"
#include "MultiIndex.h"

#ifdef USE_CPP_QHULL
#pragma comment( lib , "qhullcpp.lib" )
#pragma comment( lib , "qhullstatic_r.lib" )
#else // !USE_CPP_QHULL
#pragma comment( lib , "qhullstatic.lib" )
#endif // USE_CPP_QHULL

namespace ConvexHull
{
	//////////////////
	// Declarations //
	//////////////////

	template< unsigned int Dim > constexpr unsigned int MaxSimpleHullSize( void );
	template< unsigned int Dim > constexpr unsigned int DefaultMaxIncrementalHullSize( void ){ if constexpr( Dim<5 ) return 40 ; else return (unsigned int)-1; }

	///////////////////////////////////
	// These methods are thread-safe //
	///////////////////////////////////
	template< unsigned int Dim , unsigned int MaxIncrementalHullSize=DefaultMaxIncrementalHullSize< Dim >() > struct      ConvexHullScratch;
	template< unsigned int Dim >                                                                              struct      SimpleHullScratch;
	template< unsigned int Dim , unsigned int MaxN=-1 >                                                       struct IncrementalHullScratch;

	// MaxIncrementalHullSize: The incremental approach will be used if the number of points is less than MaxIncrementalHullSize
	template< unsigned int Dim , unsigned int MaxIncrementalHullSize >
	std::vector< SimplexIndex< Dim-1 > > ConvexHull( const std::vector< Point< double , Dim > > &points , ConvexHullScratch< Dim , MaxIncrementalHullSize > &scratch , bool generalPosition=true );

	template< unsigned int Dim >
	std::vector< SimplexIndex< Dim-1 > > SimpleHull( const std::vector< Point< double , Dim > > &point , SimpleHullScratch< Dim > &scratch );

	// MaxN: The maximum number of points for which template specialization should be used
	template< unsigned int Dim , unsigned int MaxN >
	std::vector< SimplexIndex< Dim-1 > > IncrementalHull( const std::vector< Point< double , Dim > > &points , IncrementalHullScratch< Dim , MaxN > &scratch );

	template< unsigned int Dim >
	void Orient( const std::vector< Point< double , Dim > > &points , std::vector< SimplexIndex< Dim-1 > > &hull );

	///////////////////////////////////////
	// These methods are not thread-safe //
	///////////////////////////////////////

	// MaxIncrementalHullSize: The incremental approach will be used if the number of points is less than MaxIncrementalHullSize
	template< unsigned int Dim , unsigned int MaxIncrementalHullSize=( Dim<5 ? 40 : Dim ) >
	std::vector< SimplexIndex< Dim-1 > > ConvexHull( const std::vector< Point< double , Dim > > &points , bool generalPosition=true );

	template< unsigned int Dim >
	std::vector< SimplexIndex< Dim-1 > > SimpleHull( const std::vector< Point< double , Dim > > &points );

	// MaxN: The maximum number of points for which template specialization should be used
	template< unsigned int Dim , unsigned int MaxN=(unsigned int)-1 >
	std::vector< SimplexIndex< Dim-1 > > IncrementalHull( const std::vector< Point< double , Dim > > &points );

	// This method is not thread-safe
	template< unsigned int Dim >
	std::vector< SimplexIndex< Dim-1 > > QHull( const std::vector< Point< double , Dim > > &points );

	/////////////////
	// Definitions //
	/////////////////

	template< unsigned int Dim , unsigned int MaxIncrementalHullSize >
	struct ConvexHullScratch
	{
		SimpleHullScratch< Dim > simpleHullScratch;
		IncrementalHullScratch< Dim , MaxIncrementalHullSize > incrementalHullScratch;
	};

	template< unsigned int Dim , unsigned int MaxIncrementalHullSize >
	std::vector< SimplexIndex< Dim-1 > > ConvexHull( const std::vector< Point< double , Dim > > &points, bool generalPosition )
	{
		if     ( points.size()<=MaxSimpleHullSize< Dim >() )                                              return      SimpleHull( points );
		else if( points.size()<=MaxIncrementalHullSize && MaxIncrementalHullSize!=-1 && generalPosition ) return IncrementalHull< Dim , MaxIncrementalHullSize >( points );
		else                                                                                              return           QHull( points );
	}

	template< unsigned int Dim , unsigned int MaxIncrementalHullSize >
	std::vector< SimplexIndex< Dim-1 > > ConvexHull( const std::vector< Point< double , Dim > > &points , ConvexHullScratch< Dim , MaxIncrementalHullSize > &scratch , bool generalPosition )
	{
		if     ( points.size()<=MaxSimpleHullSize< Dim >() )                                              return      SimpleHull( points , scratch.simpleHullScratch );
		else if( points.size()<=MaxIncrementalHullSize && MaxIncrementalHullSize!=-1 && generalPosition ) return IncrementalHull( points , scratch.incrementalHullScratch );
		else
		{
#ifdef USE_CPP_QHULL
			return QHull( points );
#else // !USE_CPP_QHULL
			std::vector< SimplexIndex< Dim-1 > > hull;
#pragma omp critical
			{
				hull = QHull( points );
			}
			return hull;
#endif // USE_CPP_QHULL
		}
	}

#ifdef FAST_SIMPLE_INCREMENTAL
	template< unsigned int Dim >
	struct BoundaryIncident
	{
		unsigned int incident[Dim+1][Dim+1][2];
		BoundaryIncident( void )
		{
			auto Contains = []( SimplexIndex< Dim-1 > si , SimplexIndex< Dim-2 > _si )
				{
					unsigned int count = 0;
					for( unsigned int d=0 ; d<=Dim-1 ; d++ )
						for( unsigned int _d=0 ; _d<=Dim-2 ; _d++ )
							if( si[d]==_si[_d] ) count++;
					return count==(Dim-1);
				};

			// Iterate over the set of (Dim-2)-dimensional sub-simplices
			for( unsigned int v1=0 ; v1<=Dim ; v1++ ) for( unsigned int v2=0 ; v2<v1 ; v2++ )
			{
				SimplexIndex< Dim-2 > _si = SimplexIndex< Dim >::Face( v1 , v2 );

				unsigned int idx = 0;
				// Iterate over the set of (Dim-1)-dimensional sub-simplices
				for( unsigned int v=0 ; v<=Dim ; v++ )
				{
					SimplexIndex< Dim-1 > si = SimplexIndex< Dim >::Face( v );
					// If the (Dim-1)-dimensional sub-simplex contains the (Dim-2)-dimensional sub-simplex, add it
					if( Contains( si , _si ) ) incident[v1][v2][idx++] = v;
				}
			}
		}
	};
	template<>
	struct BoundaryIncident< 1 >{};
	template< unsigned int Dim >
	constexpr unsigned int MaxSimpleHullSize( void ){ return Dim+2; }

	template< unsigned int Dim >
	struct SimpleHullScratch : BoundaryIncident< Dim >{};
#else // !FAST_SIMPLE_INCREMENTAL
	template< unsigned int Dim >
	constexpr unsigned int MaxSimpleHullSize( void ){ return Dim+1; }
#endif // FAST_SIMPLE_INCREMENTAL

#ifdef FAST_SIMPLE_INCREMENTAL
	template< unsigned int Dim >
	std::vector< SimplexIndex< Dim-1 > > SimpleHull( const std::vector< Point< double , Dim > > &points )
	{
		static SimpleHullScratch< Dim > scratch;
		return SimpleHull( points , scratch );
	}

	template< unsigned int Dim >
	std::vector< SimplexIndex< Dim-1 > > SimpleHull( const std::vector< Point< double , Dim > > &points , SimpleHullScratch< Dim > &scratch )
#else // !FAST_SIMPLE_INCREMENTAL
	template< unsigned int Dim >
	std::vector< SimplexIndex< Dim-1 > > SimpleHull( const std::vector< Point< double , Dim > > &points )
#endif // FAST_SIMPLE_INCREMENTAL
	{
		std::vector< SimplexIndex< Dim-1 > > hull;
		if( points.size()<Dim ) ERROR_OUT( "Insufficient points: " , points.size() , " >= " , Dim );
		else if( points.size()==Dim )
		{
			hull.resize( 2 );
			for( unsigned int d=0 ; d<Dim ; d++ ) hull[0][d] = d;
			hull[1] = hull[0];
			std::swap( hull[1][0] , hull[1][1] );
		}
		else if( points.size()==Dim+1 )
		{
			hull.reserve( Dim+1 );
			SimplexIndex< Dim >::template ProcessFaces< Dim-1 >( [&]( SimplexIndex< Dim-1 > si ){ hull.push_back(si); } );
			Orient( points , hull );
		}
#ifdef FAST_SIMPLE_INCREMENTAL
		else if( points.size()==Dim+2 )
		{
			// Incremental algorithm:
			// -- Use the first Dim+1 points to create a simplex with Dim+1 candidate (Dim-1)-dimensional faces
			// -- Dim of these can be extended by joining the boundary of the faces to the last point to create Dim faces

			Point< double , Dim > center;
			for( unsigned int d=0 ; d<=Dim ; d++ ) center += points[d];
			center /= (double)(Dim+1);

			// Identify which faces contain the remaining point in their shadow
			bool faceFlags[ Dim+1 ];
			for( unsigned int d=0 ; d<=Dim ; d++ )
			{
				Simplex< double , Dim , Dim-1 > face;
				for( unsigned int i=0 , _i=0 ; i<=Dim ; i++ ) if( i!=d ) face[ _i++ ] = points[i];
				Point< double , Dim > n = face.normal();
				faceFlags[d] = Point< double , Dim >::Dot( n , center - face[0] ) * Point< double , Dim >::Dot( n , points[Dim+1] - face[0] )>0;

				// Add the back-facing ...
				if( faceFlags[d] ) hull.push_back( SimplexIndex< Dim >::Face( d ) );
			}

			// Iterate over all edges and check if they are are shadow crossing
			SimplexIndex< Dim-1 > si;
			si[ Dim-1 ] = Dim+1;
			for( unsigned int d=0 ; d<=Dim ; d++ ) for( unsigned int dd=0 ; dd<d ; dd++ ) if( faceFlags[ scratch.incident[d][dd][0] ]!=faceFlags[ scratch.incident[d][dd][1] ] )
			{
				SimplexIndex< Dim-2 > _si = SimplexIndex< Dim >::Face( d , dd );
				for( unsigned int i=0 ; i<=Dim-2 ; i++ ) si[i] = _si[i];
				hull.push_back( si );
			}

			Orient( points , hull );
		}
#endif // FAST_SIMPLE_INCREMENTAL
		else ERROR_OUT( "Number of points exceeds number of points supported by simple hull: " , points.size() , " > " , MaxSimpleHullSize< Dim >() );
		return hull;
	}

	template< unsigned int Dim >
	void Orient( const std::vector< Point< double , Dim > > &points , std::vector< SimplexIndex< Dim-1> > &hull )
	{
		auto GetSimplex = []< unsigned int D , unsigned int K >( const std::vector< Point< double , D > > &points , SimplexIndex< K > si )
		{
			Simplex< double , D , K > s;
			for( unsigned int j=0 ; j<=K ; j++ ) s[j] = points[ si[j] ];
			return s;
		};

		Point< double , Dim > center;
		for( unsigned int i=0 ; i<points.size() ; i++ ) center += points[i];
		center /= (double)points.size();
		for( unsigned int i=0 ; i<hull.size() ; i++ )
		{
			Simplex< double , Dim , Dim-1 > s = GetSimplex( points , hull[i] );
			Point< double , Dim > c = s.center() , n = s.normal();
			if( Point< double , Dim >::Dot( center - c , n )>0 ) std::swap( hull[i][0] , hull[i][1] );
		}
	}


	template< unsigned int Dim >
	std::vector< SimplexIndex< Dim-1 > > QHull( const std::vector< Point< double , Dim > > &points )
	{
		std::vector< SimplexIndex< Dim-1 > > hull;
#ifdef USE_CPP_QHULL
		// Code adapted from: https://stackoverflow.com/questions/19530731/qhull-library-c-interface
		std::string comment = ""; // rbox commands, see http://www.qhull.org/html/rbox.htm
		std::string qhull_command = ""; // For qhull commands, see http://www.qhull.org/html/qhull.htm

#ifdef DISABLE_QHULL_TRY_CATCH
#else // !DISABLE_QHULL_TRY_CATCH
		try
#endif // DISABLE_QHULL_TRY_CATCH
		{
			orgQhull::Qhull qhull = orgQhull::Qhull( comment.c_str() , (int)Dim , (int)points.size() , &points[0][0] , qhull_command.c_str() );
			hull.resize( qhull.facetList().size() );

			unsigned int fCount = 0;
			for( const orgQhull::QhullFacet &f : qhull.facetList() )
			{
				unsigned int d=0;
				for( const orgQhull::QhullVertex &v : f.vertices() ) hull[fCount][d++] = ( ( v.getVertexT()->point - &points[0][0] )/Dim );
				fCount++;
			}
		}
#ifdef DISABLE_QHULL_TRY_CATCH
#else // !DISABLE_QHULL_TRY_CATCH
		catch( orgQhull::QhullError &e ){ ERROR_OUT( "qhull failure: " , std::string(e.what() ) ); }
#endif // DISABLE_QHULL_TRY_CATCH
#else // !USE_CPP_QHULL
		qh_init_A( NULL , NULL , stderr , 0 , NULL );  /* sets qh qhull_command */
		int exitcode = setjmp(qh errexit);	/* simple statement for CRAY J916 */
		qh JOGGLEmax = 0.0;					/* 'QJ'   */
		qh PRINTprecision = False;

#if 0
		unsigned int ismalloc = False;
		qh_init_B( &points[0][0] , (int)points.size() , Dim , ismalloc);
#else
		unsigned int ismalloc = True;
		double *pCoordinates = (double*)qh_malloc( sizeof(double)*points.size()*Dim );
		memcpy( pCoordinates , &points[0][0] , sizeof(double)*points.size()*Dim );
		qh_init_B( pCoordinates , (int)points.size() , Dim , ismalloc);
#endif
		qh_qhull();
		qh_check_output();
		if ( qh facet_list )
		{
			unsigned int fCount = 0;
			for ( facetT *facet=( qh facet_list ); facet && facet->next; facet=facet->next ) fCount++;
			hull.resize( fCount );
			fCount = 0;
			for ( facetT *facet=( qh facet_list ); facet && facet->next; facet=facet->next )
			{
				setT* verts = facet->vertices;
				vertexT *vertex, **vertexp;

#define FOREACHsetelement_(type, set, variable) \
		if (((variable= NULL), set)) for (\
		variable##p= (type **)&((set)->e[0].p); \
		(variable= *variable##p++);)

				unsigned int sz=0;
				SimplexIndex< Dim-1 > simplex;
				FOREACHsetelement_( vertexT, verts , vertex ) simplex[sz++] = qh_pointid( vertex->point );
				if( sz!=Dim ) ERROR_OUT( "Degenerate facet: " , sz );
				hull[ fCount++ ] = simplex;
#undef FOREACHsetelement_
			}
		}

		qh NOerrexit= True;  /* no more setjmp */
#ifdef qh_NOmem
		qh_freeqhull( True);
#else
		qh_freeqhull( False);
		int curlong , totlong;
		qh_memfreeshort(&curlong, &totlong);
		if (curlong || totlong)
			fprintf(stderr, "qhull internal warning (main): did not free %d bytes of long memory(%d pieces)\n",
				totlong, curlong);
#endif
#endif // USE_CPP_QHULL

		Orient( points , hull );
		return hull;
	}

	template< unsigned int Dim , unsigned int MaxN >
	struct IncrementalHullScratch
	{
		struct BoundaryInfo
		{
			unsigned int incidentFaces[2];
			SimplexIndex< Dim-2 > fi;

			BoundaryInfo( void ){ incidentFaces[0] = incidentFaces[1] = -1; }

			unsigned int &operator[]( unsigned int i ){ return incidentFaces[i]; }
			const unsigned int &operator[]( unsigned int i ) const { return incidentFaces[i]; }
		};

		struct FaceInfo
		{
			SimplexIndex< Dim-1 > si;
			Point< double , Dim > normal;
			char inShadow;

			FaceInfo( void ){}
			FaceInfo( SimplexIndex< Dim-1 > si , Point< double , Dim > normal , char inShadow=false ) : si(si) , normal(normal) , inShadow(inShadow){}
		};

		// [WARNING] The storage size of the map is Choose( MaxN , Dim )
		typename MultiIndex< Dim-1 , MaxN >::template Map< BoundaryInfo > boundaryInfo;
		std::vector< FaceInfo > faces , _faces;
	};

	template< unsigned int Dim , unsigned int MaxN >
	std::vector< SimplexIndex< Dim-1 > > IncrementalHull( const std::vector< Point< double , Dim > > &points )
	{
		static IncrementalHullScratch< Dim , MaxN > scratch;
		return IncrementalHull( points , scratch );
	}

	template< unsigned int Dim , unsigned int MaxN >
	std::vector< SimplexIndex< Dim-1 > > IncrementalHull( const std::vector< Point< double , Dim > > &points , IncrementalHullScratch< Dim , MaxN > &scratch )
	{
		// [DEFINITON] A face is "in shadow" w.r.t. to a point if it is back-facing to the point
		// Need to remove the faces that are not "in shadow" and replace with new faces

		if constexpr( MaxN!=-1 ) if( points.size()>MaxN ) ERROR_OUT( "Point size mismatch: " , points.size() , " <= " , MaxN );

		if( points.size()<=Dim ) ERROR_OUT( "not enough points" );

		struct FaceBoundaries
		{
			SimplexIndex< Dim-2 > boundaries[ Dim ];
			SimplexIndex< Dim-2 > &operator[]( unsigned int idx ){ return boundaries[idx]; }
			const SimplexIndex< Dim-2 > &operator[]( unsigned int idx ) const { return boundaries[idx]; }
			FaceBoundaries( void ){ for( unsigned int d=0 ; d<=Dim-1 ; d++ ) boundaries[d] = SimplexIndex< Dim-1 >::Face( d ); }
		};

		auto FaceNormal = [&]( SimplexIndex< Dim-1 > fi )
			{
				Simplex< double , Dim , Dim-1 > face;
				for( unsigned int d=0 ; d<=Dim-1 ; d++ ) face[d] = points[ fi[d] ];
				return face.normal();
			};


		static const FaceBoundaries _FaceBoundaries;
		std::vector< typename IncrementalHullScratch< Dim , MaxN >::FaceInfo > &faces = scratch.faces , &_faces = scratch._faces;
		typename MultiIndex< Dim-1 , MaxN >::template Map< typename IncrementalHullScratch< Dim , MaxN >::BoundaryInfo > &boundaryInfo = scratch.boundaryInfo;

		Point< double , Dim > center;

		// Set-up the initial simplices
		{
			faces.resize(0);
			SimplexIndex< Dim > si;
			for( unsigned int d=0 ; d<=Dim ; d++ ) si[d] = d , center += points[d];
			for( unsigned int d=0 ; d<=Dim ; d++ )
			{
				SimplexIndex< Dim-1 > fi = si.face(d);
				faces.emplace_back( fi , FaceNormal( fi ) );
			}
			center /= (double)(Dim+1);
		}

		boundaryInfo.resize( (unsigned int)points.size() );
		// Add the remaining points
		for( unsigned int v=Dim+1 ; v<points.size() ; v++ )
		{
			boundaryInfo.clear();
			_faces.resize(0);

			// Start by marking each faces as in/out of shadow
			for( unsigned int f=0 ; f<faces.size() ; f++ )
			{
				faces[f].inShadow = Point< double , Dim >::Dot( faces[f].normal , center - points[ faces[f].si[0] ] ) * Point< double , Dim >::Dot( faces[f].normal , points[v] - points[ faces[f].si[0] ] )>0;
				if( faces[f].inShadow ) _faces.emplace_back( faces[f] );
			}

			// Set the boundary information

			for( unsigned int f=0 ; f<faces.size() ; f++ ) for( unsigned int d=0 ; d<=Dim-1 ; d++ )
			{
				SimplexIndex< Dim-2 > fi;
				const SimplexIndex< Dim-2 > &_fi = _FaceBoundaries[d];
				for( unsigned int _d=0 ; _d<=Dim-2 ; _d++ ) fi[_d] = faces[f].si[ _fi[_d] ];
				MultiIndex< Dim-1 , MaxN > mi( &fi[0] );
				typename IncrementalHullScratch< Dim , MaxN >::BoundaryInfo &bi = boundaryInfo[ mi ];
				if( bi[0]==-1 ) bi[0] = f , bi.fi = fi;
				else            bi[1] = f;
			}

			SimplexIndex< Dim-1 > fi;
			fi[Dim-1] = v;
			for( unsigned int i=0 ; i<boundaryInfo.size() ; i++ )
			{
				const typename IncrementalHullScratch< Dim , MaxN >::BoundaryInfo &bi = boundaryInfo[i];
				if( faces[ bi[0] ].inShadow!=faces[ bi[1] ].inShadow )
				{
					for( unsigned int d=0 ; d<=Dim-2 ; d++ ) fi[d] = bi.fi[d];
					_faces.emplace_back( fi , FaceNormal(fi) );
				}
			}

			std::swap( faces , _faces );
		}

		std::vector< SimplexIndex< Dim-1 > > hull( faces.size() );
		for( unsigned int i=0 ; i<faces.size() ; i++ ) hull[i] = faces[i].si;
		Orient( points , hull );
		return hull;
	}
}
#endif // CONVEX_HULL_INCLUDED