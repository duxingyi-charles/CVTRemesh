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

#ifndef __DELAUNAY_IO_H__
#define __DELAUNAY_IO_H__

//#include "delaunay.h"
#include <Geex/graphics/opengl.h>

namespace Geex {

class DelaunayIO {
public:
	//DelaunayIO(Delaunay *delaunay) ;
	DelaunayIO();
	~DelaunayIO() ;

	//void generate_pointset(double radius, int n, const std::string& fname) ;
	//void save_points(const std::string& fname) ;
	//void save_edge_histogram(const std::string& fname) ;
	//void save_area_histogram(const std::string& fname) ;
	//void save_angle_histogram(const std::string& fname) ;
	//void save_vertex_degree(const std::string& name) ;
	//void spetral_analysis(const std::string& fname) ;
	bool snapshot(const std::string& filename) ;
	//std::vector<double>& edge_histogram() { return edge_hist_ ; }
	//std::map<int, std::vector<double> >& val_edge_histogram() { return val_edge_hist_ ;} 
	//GLboolean& show_edge_hist() { return show_edge_hist_ ; }

	//void analyze_pointsets(const std::string& path) ;

private:
	//Delaunay* delaunay_ ;
// 	std::map<int, double> ave_vdeg_ ;
// 	std::map<std::pair<int, int>, double> ave_edeg_ ;
// 	std::vector<double> edge_hist_ ;
// 	std::map<int, std::vector<double> > val_edge_hist_ ;
// 	std::vector<double> angle_hist_ ;
// 	std::vector<double> area_hist_ ;
// 	GLboolean show_edge_hist_ ;

} ;
    
}

#endif
