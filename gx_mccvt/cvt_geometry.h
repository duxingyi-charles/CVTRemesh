/*
 *  _____ _____ ________  _
 * /  __//  __//  __/\  \//
 * | |  _|  \  |  \   \  / 
 * | |_//|  /_ |  /_  /  \ 
 * \____\\____\\____\/__/\\
 *
 * Graphics Environment for EXperimentations.
 *  Copyright (C) 2006 INRIA - Project ALICE
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: 
 *
 *     ALICE Project - INRIA
 *     INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs. 
 */

#ifndef __GEOMETRY_APD_MESH__ 	
#define __GEOMETRY_APD_MESH__

#include <Geex/symbolics/expression.h>
#include <Geex/mathematics/glsl_linear.h>

namespace Geex {


    //template <class T> inline T det2x21(
    //    T a00, T a01,
    //    T a10, T a11
    //) {
    //    return a00*a11 - a10 * a01 ;
    //}

	template <class T> void invert_mat(T m[16], T i[16]) {
		vec3g<T> t(m[12],m[13],m[14]);
		vec3g<T> o=vec3g<T>(vec3g<T>(m[0],m[1],m[2])*t,vec3g<T>(m[4],m[5],m[6])*t,vec3g<T>(m[8],m[9],m[10])*t); 
		i[0] =m[0]; i[1] =m[4]; i[2] =m[8]; i[3] =0;
		i[4] =m[1]; i[5] =m[5]; i[6] =m[9]; i[7] =0;
		i[8] =m[2]; i[9] =m[6]; i[10]=m[10];i[11]=0;
		i[12]=-o.x; i[13]=-o.y; i[14]=-o.z; i[15]=1;
	}

    //template <class T> inline T triangle_area(
    //    const vec2g<T>& p, const vec2g<T>& q, const vec2g<T>& r
    //) {
    //    vec2g<T> U = q-p ;
    //    vec2g<T> V = r-p ;
    //    return T(1) / T(2) * det2x21(U.x,V.x,U.y,V.y) ;
    //}

	template <class T> inline vec3g<T> gx_min(vec3g<T>& v1, vec3g<T>& v2) {
		return vec3g<T>(gx_min(v1.x, v2.x), gx_min(v1.y, v2.y), gx_min(v1.z, v2.z)) ;
	}

	template <class T> inline vec3g<T> gx_max(vec3g<T>& v1, vec3g<T>& v2) {
		return vec3g<T>(gx_max(v1.x, v2.x), gx_max(v1.y, v2.y), gx_max(v1.z, v2.z)) ;
	}

    inline double triangle_area(
        const vec3& v1, const vec3& v2, const vec3& v3
    ) {
        return 0.5f * length(cross(v2-v1, v3-v1)) ;
    }

    template <class T> inline T tetra_radius_ratio(
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r, const vec3g<T>& s
    ) {
		return T(3)*tetra_incirbedradius(p, q, r, s)/tetra_circumradius(p, q, r, s);
	}

    template <class T> inline T tetra_incirbedradius(
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r, const vec3g<T>& s
    ) {
		return T(3)*tetra_volume(p, q, r, s)/(tri_area(q, r, s)+tri_area(r, s, p)+tri_area(s, p, q)+tri_area(p, r, q)) ;
	}

    template <class T> inline T tetra_circumradius(
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r, const vec3g<T>& s
    ) {
		return distance(p, tetra_circumcenter(p, q, r, s) ) ;
	}

	template <class T> inline bool on_segment(const vec3g<T>& p, const vec3g<T>& a, const vec3g<T>& b) {
		if(distance(p, a) + distance(p, b) == distance(a, b))
			return true ;
		return false ; 
	}

	template <class T> inline bool inside_triangle(const vec3g<T>& p, const vec3g<T>& a, const vec3g<T>& b, const vec3g<T>& c) {
		vec3g<T> n = cross(b-a, c-a);
		//T a0 = dot(cross(a-p, b-p), n); 
		//T a1 = dot(cross(b-p, c-p), n); 
		//T a2 = dot(cross(c-p, a-p), n);
		if(dot(cross(a-p, b-p), n)<-1e-10 || dot(cross(b-p, c-p), n)<-1e-10 || dot(cross(c-p, a-p), n)<-1e-10) 
			return false;
		else 
			return true;
	}

	template <class T> inline vec3g<T> braycenter_coord(const vec3g<T>& p, const vec3g<T>& a, const vec3g<T>& b, const vec3g<T>& c) {
		gx_assert(inside_triangle(p, a, b, c)) ;
		T area = cross(b-a, c-a).length() ;
		T a0 = cross(b-p, c-p).length() ; 
		T a1 = cross(c-p, a-p).length() ;
		T a2 = cross(a-p, b-p).length() ;
		return vec3g<T>(a0, a1, a2)/area ;
	}	
		
	template <class T> inline bool ray_triangle_intersection(vec3g<T>& p0, vec3g<T>& dir, vec3g<T>& a, vec3g<T>& b, vec3g<T>& c, vec3g<T>& n, vec3g<T>& intp) {
		Plane<T> P(a, n);
		if(ray_plane_intersects(p0, dir, P)) {
			intp = ray_plane_intersection(p0, dir, P);
			if(inside_triangle(intp, a, b, c))
				return true;
		}
		return false;
	}

	//template <class T> inline bool ray_plane_intersects(vec3g<T>& p0, vec3g<T>& dir, const Plane<T>& plane) {
	//	bool side = plane.side(p0) > 0.0;
	//	bool sign = dot(dir, plane.normal()) > 0.0;

	//	if(side && sign || !side && !sign)
	//		return false;
	//	return true;
	//}

	//template <class T> inline bool seg_plane_intersects(vec3g<T>& p0, vec3g<T>& p1, const Plane<T>& plane) {
	//	if(plane.side(p0)>0.0 != plane.side(p1)>0.0)
	//		return true;
	//	return false;
	//}

	//template <class T> inline bool seg_triangle_intersection(vec3g<T>& p0, vec3g<T>& p1, vec3g<T>& a, vec3g<T>& b, vec3g<T>& c, vec3g<T>& n, vec3g<T>& intp) {
	//	Plane<T> P(a, n);
	//	if(seg_plane_intersects(p0, p1, P)) {
	//		intp = seg_plane_intersection(p0, p1, P);
	//		if(inside_triangle(intp, a, b, c))
	//			return true;
	//	}
	//	return false;
	//}

	//template <class T> inline vec3g<T> project_to_plane(const vec3g<T>& p, const Plane<T>& plane) {
	//	return p-plane.side(p)*plane.normal() ;
	//}

	template <class T> inline bool project_to_triangle(const vec3g<T>& p, const vec3g<T>& a, const vec3g<T>& b, const vec3g<T>& c, vec3g<T>& proj) {
		vec3g<T> pc = p - c;
		vec3g<T> ac = a - c;
		vec3g<T> bc = b - c;
		T a00 = ac.length2();
		T a01 = dot(ac, bc);
		T a11 = bc.length2();
		T b0 = dot(pc, ac);
		T b1 = dot(pc, bc);
		T mdet = a00 * a11 - a01 * a01;
		T s = (a11 * b0 - a01 * b1) / mdet;
		T t = (a00 * b1 - a01 * b0) / mdet;
		proj = s * a + t * b + (1 - s - t ) * c;
		if ( s < 0 || s > 1 || t < 0 || t > 1) 
		{
			return false;
		}
		else
			return true;
		//vec3g<T>& n = cross(b-a, c-a);
		//n = normalize(n);
		//proj = project_to_plane(p, Plane<T>(a, n));
		//if(inside_triangle(proj, a, b, c))
		//	return true;
		//return false;
	}

	template <class T> inline void project_to_linesegment(const vec3g<T>& p, const vec3g<T>& a, const vec3g<T>& b, vec3g<T> & fp, T & sqdist)
	{
		//double t;

		vec3g<T> ab = a - b;
		vec3g<T> pb = p - b;
		T t = dot(pb, ab) / ab.length2();
		//t = (-a.x*b.x + b.x*b.x - a.y*b.y + b.y*b.y - a.z*b.z + b.z*b.z + a.x*p.x
		//	-b.x*p.x +a.y*p.y-b.y*p.y+a.z*p.z-b.z*p.z) 
		//	/ 
		//	(a.x*a.x+a.y*a.y+a.z*a.z-2.*a.x*b.x-2.*a.y*b.y-2.*a.z*b.z+b.x*b.x+b.y*b.y+b.z*b.z);

		if(t < 0)
		{
			fp = b;
			sqdist = pb.length2();
			return;
		}
		else if(t > 1)
		{
			fp = a;
			sqdist = (p-a).length2();
			return;
		}
		else
		{
			fp = t*a+(1-t)*b;
			sqdist = (p-fp).length2();
			return;
		}
	}

    //template <class T> inline T tri_area(
    //    const vec3g<T>& p, const vec3g<T>& q, 
    //    const vec3g<T>& r
    //) {
    //    return length(cross(q-p, r-p)) / T(2);
    //}

    //template <class T> inline vec3g<T> tri_centroid(
    //    const vec3g<T>& p, const vec3g<T>& q, 
    //    const vec3g<T>& r
    //) {
    //    return T(1)/T(3)*(p+q+r);
    //}

	//template <class T> inline T tri_mass(
 //       const vec3g<T>& p, const vec3g<T>& q, const vec3g<T>& r,
	//	T a, T b, T c // a, b, c are density at p, q, r
 //   ) {
	//	return tri_area(p, q, r)/T(3)*(sqrt(fabs(a))+sqrt(fabs(b))+sqrt(fabs(c))) ;
	//}

	//template <class T> inline void tri_centroid(
	//	const vec3g<T>& p, const vec3g<T>& q, const vec3g<T>& r,
	//	vec3g<T>& g, T& V
	//	) {
	//		V = tri_area(p, q, r) ;
	//		g = T(1.)/T(3.) * ( p+q+r ) ;
	//}

 //   template <class T> inline void tri_centroid(
 //       const vec3g<T>& p, const vec3g<T>& q, const vec3g<T>& r,
	//	T a, T b, T c, vec3g<T>& Vg, T& V // a, b, c are density at p, q, r
 //   ) {
	//	T abc = a+b+c ; 
	//	T area = tri_area(p, q, r) ;
	//	V = area/T(3)*abc ;
	//	Vg = (area/T(12))*((a+abc)*p+(b+abc)*q+(c+abc)*r) ;
	//}

    template <class T> inline void tetra_centroid(
        const vec3g<T>& p, const vec3g<T>& q, const vec3g<T>& r, const vec3g<T>& s,
		T a, T b, T c, T d, vec3g<T>& g, T& V // a, b, c d are density at p, q, r, s
    ) {
		T abcd = a+b+c+d ; 
		V = tetra_volume(p, q, r, s) ;
		g =  T(1) / (5*abcd) *  ( (a+abcd)*p+(b+abcd)*q+(c+abcd)*r+(d+abcd)*s );
	}

    template <class T> inline bool inside_tetra(const vec3g<T>& p0, 
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r, const vec3g<T>& s
    ) {
		vec3 dir[4];
		dir[0] = p - p0;
		dir[1] = q - p0;
		dir[2] = r - p0;
		dir[3] = s - p0;

		real msignvolume[4] ;
		msignvolume[0] = Geex::dot(dir[3], Geex::cross(dir[1], dir[2]));
		msignvolume[1] = -Geex::dot(dir[3], Geex::cross(dir[0], dir[2]));
		msignvolume[2] = Geex::dot(dir[3], Geex::cross(dir[0], dir[1]));
		msignvolume[3] = -Geex::dot(dir[2], Geex::cross(dir[0], dir[1]));

		if (msignvolume[0]>= 0 && 
			msignvolume[1]>= 0 && 
			msignvolume[2]>= 0 && 
			msignvolume[3]>= 0 )
		{
			return true;
		}
		return false;
	}

    template <class T> inline bool inside_tetra(const vec3g<T>& p0, 
        const vec3g<T>& p, const vec3g<T>& q, 
        const vec3g<T>& r, const vec3g<T>& s, real msignvolume[] 
    ) {
		vec3 dir[4];
		dir[0] = p - p0;
		dir[1] = q - p0;
		dir[2] = r - p0;
		dir[3] = s - p0;

		msignvolume[0] = Geex::dot(dir[3], Geex::cross(dir[1], dir[2]));
		msignvolume[1] = -Geex::dot(dir[3], Geex::cross(dir[0], dir[2]));
		msignvolume[2] = Geex::dot(dir[3], Geex::cross(dir[0], dir[1]));
		msignvolume[3] = -Geex::dot(dir[2], Geex::cross(dir[0], dir[1]));

		if (msignvolume[0]>= 0 && 
			msignvolume[1]>= 0 && 
			msignvolume[2]>= 0 && 
			msignvolume[3]>= 0 )
		{
			return true;
		}
		return false;
	}

	template<class T> inline bool inside_sphere(const vec3g<T>& p, const vec3g<T>& o, const T r2) {
		return distance2(p, o) >= r2 ;
	}

	template<class T> inline bool segment_sphere_intersects(
		const vec3g<T>& a, 
		const vec3g<T>& b, 
		const vec3g<T>& o,
		const T rr) {

		vec3g<T> n = b-a ;
		double A = dot(n, n) ;
		double B = 2*dot(a-o, n) ;
		double C = dot(a-o, a-o) - rr ;
		double det = B*B-4*A*C ;

		if(det<0) return false ;
		double t0 = (-B-sqrt(det))/(2*A) ;
		double t1 = (-B+sqrt(det))/(2*A) ;
		if(t0>0 && t0<1 || t1>0 && t1<1) 
			return true ;
		return false ;
	}

	template<class T> inline bool triangle_sphere_intersects(
		const vec3g<T>& a, 
		const vec3g<T>& b, 
		const vec3g<T>& c, 
		const vec3g<T>& o,
		const T rr) {

		vec3g<T> A = a - o ;
		vec3g<T> B = b - o ;
		vec3g<T> C = c - o ;
		vec3g<T> V = normalize(cross(B - A, C - A)) ;
		double   d = dot(A, V) ;
		double   e = dot(V, V) ;
		bool     sep1 = d * d > rr * e ;
		if(sep1) 
			return false ;

		double aa = dot(A, A) ;
		double ab = dot(A, B) ;
		double ac = dot(A, C) ;
		double bb = dot(B, B) ;
		double bc = dot(B, C) ;
		double cc = dot(C, C) ;
		bool   sep2 = (aa > rr) && (ab > aa) && (ac > aa) ;
		bool   sep3 = (bb > rr) && (ab > bb) && (bc > bb) ;
		bool   sep4 = (cc > rr) && (ac > cc) && (bc > cc) ;
		if(sep2 && sep3 && sep4) 
			return false ;

		vec3g<T> AB = B - A ;
		vec3g<T> BC = C - B ;
		vec3g<T> CA = A - C ;
		double   d1 = ab - aa ;
		double   d2 = bc - bb ;
		double   d3 = ac - cc ;
		double   e1 = dot(AB, AB) ;
		double   e2 = dot(BC, BC) ;
		double   e3 = dot(CA, CA) ;
		vec3g<T> Q1 = A * e1 - d1 * AB ;
		vec3g<T> Q2 = B * e2 - d2 * BC ;
		vec3g<T> Q3 = C * e3 - d3 * CA ;
		vec3g<T> QC = C * e1 - Q1 ;
		vec3g<T> QA = A * e2 - Q2 ;
		vec3g<T> QB = B * e3 - Q3 ;
		bool     sep5 = (dot(Q1, Q1) > rr * e1 * e1) && (dot(Q1, QC) > 0) ;
		bool     sep6 = (dot(Q2, Q2) > rr * e2 * e2) && (dot(Q2, QA) > 0) ;
		bool     sep7 = (dot(Q3, Q3) > rr * e3 * e3) && (dot(Q3, QB) > 0) ;
		bool     separated = sep1 | sep2 | sep3 | sep4 | sep5 | sep6 | sep7 ;

		if(sep5 && sep6 && sep7) 
			return false ;
		return true ;
	}

    template<class T> inline bool inside_box(const vec3g<T>& p, vec3g<T>& bmin, vec3g<T>& bmax) {
	if(p[0]<bmin[0] || p[0]>bmax[0] || p[1]<bmin[1] || p[1]>bmax[1] || p[1]<bmin[1] || p[1]>bmax[1])
	    return false ; 
	return true ; 
    }
    
    template<class T> inline bool is_nan(vec3g<T>& p) {
	if(Numeric::is_nan(p.x) || Numeric::is_nan(p.y) || Numeric::is_nan(p.z)) {
	    std::cerr << "Nan !" << std::endl ;
	    return true;
	}
	return false ;
    }

	inline double to_degree(const double angle) { return angle/M_PI*180; }
	inline double to_arc(const double angle) { return angle/180*M_PI; }

	template<class T> inline vec3g<T> random_point_tri(vec3g<T>& p0, vec3g<T>& p1, vec3g<T>& p2) {
		real s = Numeric::random_float64() ;
		real t = Numeric::random_float64() ;
		if (s + t > 1)
		{
			s = 1 - s;
			t = 1 - t;
		}
		return (1-s-t)*p0 + s*p1 + t*p2 ;
	}

	template<class T> inline vec3g<T> random_point_quad(
		vec3g<T>& p0, vec3g<T>& p1, vec3g<T>& p2, vec3g<T>& p3
		) {
        real s = Numeric::random_float64() ;
        real t = Numeric::random_float64() ;
        real u = Numeric::random_float64() ;
        if (s+t>1.0) { // cut'n fold the cube into a prism
            s = 1.0 - s ;
            t = 1.0 - t ;
        }
        if (t+u>1.0) { // cut'n fold the prism into a tetrahedron
            real tmp = u ;
            u = 1.0 - s - t ;
            t = 1.0 - tmp ;
        } else if (s+t+u>1.0) {
            real tmp = u ;
            u = s + t + u - 1.0 ;
            s = 1 - t - tmp ;
        }
        real a=1-s-t-u ; // a,s,t,u are the barycentric coordinates of the random point.
		return a*p0 + s*p1 + t*p2 + u*p3 ;
	}

	template<class T> inline vec3g<T> random_point_poly(std::vector<vec3g<T> >& poly) {
		unsigned int n = poly.size() ;
		std::vector<double> area ;
		std::vector<double> cpdf ;
		double total = 0 ;
		double t = Numeric::random_float64() ;
		//double t = Random::UniformRandom() ;

		for(unsigned int i=0; i<n-2; ++i) {
			area.push_back( total+triangle_area(poly[0], poly[i+1], poly[i+2]) );
			total = area[i] ;
		}
		
		for(unsigned int i=0; i<n-2; ++i) {
			cpdf.push_back( area[i]/total );
		}

		for(unsigned int i=0; i<n-2; ++i) {
			if(t<=cpdf[i]) {
				return random_point_tri(poly[0], poly[i+1], poly[i+2]) ;
			}
		}
		gx_assert(false) ;
		return vec3g<T>(T(0), T(0), T(0)) ;
	}

	/**
		randomly select a poly by area: cpdf is normalized cumulative area
	*/
	template<class T> inline int binary_search(std::vector<T>& value, int s, int e, T a)  {
		int mid = (s+e)/2 ;
		if(s==e) {
			return s ;
		}
		if(a>value[mid] && a<=value[mid+1]) {
			return mid+1 ;
		}
		else if(a<value[mid]) {
			return binary_search(value, s, mid, a) ;
		}
		else{
			return binary_search(value, mid+1, e, a) ;
		}
	}
	
	template<class T> inline int random_polygon_by_area(std::vector<T>& value, T a) {
		return binary_search(value, T(0), value.size()-1, a) ;
	}

	template<class T> inline int seg_circle_intersection(vec2g<T>& v0, vec2g<T>& v1, vec2g<T>& o, T r, vec2g<T>& p) {
		vec2g<T> d = normalize(v1 - v0);
		T dov0 = (o-v0).length() ;
		T dpro = dot(o-v0, d) ;
		T h = sqrt(dov0*dov0 - dpro*dpro) ;
		T d0 = (v0-o).length() ;
		T d1 = (v1-o).length() ;

		// v0 out, v1 in
		if(d0>r && d1<r) {
			T t = dpro - sqrt(r*r-h*h) ;
			p = v0 + t*d ;
		}
		// v0 in v1 out
		else if(d1>r && d0<r) {
			T t = dpro + sqrt(r*r-h*h) ;
			p = v0 + t*d ;
		}
		else {
			return 0 ;
		}
	} 

	template<class T> inline int arc_circle_intersection(vec2g<T>& v0, vec2g<T>& v1, vec2g<T>& o, T r, vec2g<T>& p) {
		vec2g<T> d01 = normalize(v1 - v0) ;
		T dist = (v1-v0).length() ;
		T hc  = sqrt(r*r - dist*dist/4) ;
		vec2g<T>   c01 = 0.5*(v0+v1) + hc*vec2(d01.y, -d01.x) ;
		vec2g<T>   dir = normalize(o - c01) ;
		T doc = (o-c01).length() ;
		T h = sqrt(r*r - doc*doc/4) ;
		T d0 = (v0-o).length() ;
		T d1 = (v1-o).length() ;

		// v0 out v1 in
		if(d0>r && d1<r) {
			p = 0.5*(c01 + o) + h * vec2(-dir.y, dir.x) ;
		}
		// v0 in v1 out
		else if (d1>r && d0<r) {
			p = 0.5*(c01 + o) - h * vec2(-dir.y, dir.x) ;
		}
		else {
			return 0 ;
		}
	}

	template<class T> inline T dist_point_segment(vec2g<T>& p, vec2g<T>& v0, vec2g<T>& v1) {
		return dot(p-v0, normalize(v1-v0)) ;
	}

	template<class T> inline T power_distance(const vec3g<T>& p1, const T w1, const vec3g<T>& p2, const T w2) {		
		return distance2(p1, p2)-w1-w2 ;
	}

	//template<class T> inline bool is_obtuse(const vec3g<T>& p, const vec3g<T>& q, const vec3g<T>& r) {
	//	T dpq = distance2(p, q) ;
	//	T dqr = distance2(q, r) ;
	//	T drp = distance2(r, p) ;

	//	return dpq+dqr<drp || dpq+drp<dqr || drp+dqr<dpq ;
	//}

	template<class T> inline unsigned is_obtuse(const vec3g<T>& p, const vec3g<T>& q, const vec3g<T>& r) {
		T dpq = distance2(p, q) ;
		T dqr = distance2(q, r) ;
		T drp = distance2(r, p) ;
		if(drp+dqr<dpq)
			return 3 ;
		else if(dpq+drp<dqr)
			return 1 ;
		else if(dpq+dqr<drp)
			return 2 ;
		else return 0 ;
		//return dpq+dqr<drp || dpq+drp<dqr || drp+dqr<dpq ;
	}

	template<class T> inline void min_max_angle(const vec3g<T>& p, const vec3g<T>& q, const vec3g<T>& r,
		T& amin, T& amax) {
		double d01sqr = distance2(p, q) ;
		double d12sqr = distance2(q, r) ;
		double d20sqr = distance2(r, p) ;

		double d01 = sqrt(d01sqr) ;
		double d12 = sqrt(d12sqr) ;
		double d20 = sqrt(d20sqr) ; 

		double a0 = acos((d12sqr+d20sqr-d01sqr)/(2*d12*d20))*180/M_PI ;
		double a1 = acos((d01sqr+d20sqr-d12sqr)/(2*d01*d20))*180/M_PI ;
		double a2 = acos((d01sqr+d12sqr-d20sqr)/(2*d01*d12))*180/M_PI ;

		amin = gx_min(a0, gx_min(a1, a2)) ;
		amax = gx_max(a0, gx_max(a1, a2)) ;
	}
}

#endif
