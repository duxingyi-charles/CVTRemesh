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
 *  You should have received a copy of the GNU General Public License<
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

#include "delaunay_io.h"
#include <Geex/third_party/png/png.h>
#include <glut_viewer/glut_viewer.h>
//#include <psa-src/fileutil.h>
//#include <psa-src/fileio.h>

#include <fstream>

namespace Geex {
// 	DelaunayIO::DelaunayIO(Delaunay *delaunay) 
// 	:delaunay_(delaunay) {
// 		show_edge_hist_ = GL_FALSE ;
// 	}

	DelaunayIO::DelaunayIO() {}
	DelaunayIO::~DelaunayIO() {
	}

#ifndef SQR(x)
#define SQR(x) ((x)*(x))
#endif 
// // // // 	void DelaunayIO::generate_pointset(double radius, int n, const std::string& fname) {
// // // // 		char str[1024] ;
// // // // 		int  npoints = 0 ;
// // // // 		std::string method ;
// // // // 		switch(delaunay_->sampling_mode()) {
// // // // 		case PD_DARTTHROW:
// // // // 			method = "dart" ;
// // // // 			break ;
// // // // 		case PD_VORONOI:
// // // // 			method = "voro" ;
// // // // 			break ;
// // // // 		default:
// // // // 			std::cerr << "algorithm is not implemented..." << std::endl ;
// // // // 			return ;
// // // // 			break ;
// // // // 		}
// // // // 
// // // // 		ave_vdeg_.clear() ;
// // // // 		ave_edeg_.clear() ;
// // // // 		edge_hist_.assign(256, 0) ;
// // // // 		show_edge_hist_ = FALSE ;
// // // // 		delaunay_->sample_radius() = radius ;
// // // // 		
// // // // 		val_edge_hist_.clear() ;
// // // // 		for(int i=3; i<=10; ++i) {
// // // // 			val_edge_hist_[i].assign(256, 0) ;
// // // // 		}
// // // // 
// // // // 		for(int i=0; i<n; ++i) {
// // // // 			delaunay_->generate_poisson_disk() ;
// // // // 
// // // // 			// save points
// // // // 			sprintf(str, "%s%s_r%f_%d.txt", fname.c_str(), method.c_str(), radius, i) ;
// // // // 			save_points(str) ;
// // // // 
// // // // 			// save vertex degree 
// // // // 			sprintf(str, "%s%s_r%f_deg_%d.txt", fname.c_str(), method.c_str(), radius, i) ;
// // // // 			save_vertex_degree(str) ;
// // // // 
// // // // 			// save edge histogram
// // // // 			sprintf(str, "%s%s_r%f_hist_%d.txt", fname.c_str(), method.c_str(), radius, i) ;
// // // // 			save_edge_histogram(str) ;
// // // // 
// // // // 			npoints += delaunay_->nb_vertices() ;
// // // // 		}
// // // // 		// average number of points 
// // // // 		npoints /= n ; 
// // // // 
// // // // 		double density = npoints*M_PI*SQR(radius/2.0) ;
// // // // 		std::cerr << "packing density " << density << std::endl ;
// // // // 
// // // // 		sprintf(str, "%s%s_r%f_deg_ave.txt", fname.c_str(), method.c_str(), radius) ;
// // // // 		std::ofstream out(str) ;
// // // // 		out << "average number of points" << std::endl ;
// // // // 		out << npoints << std::endl ;
// // // // 		out << "average vertex degrees" << std::endl ;
// // // // 		for(std::map<int, double>::iterator it=ave_vdeg_.begin(); it!=ave_vdeg_.end(); ++it) {
// // // // 			it->second/=n ;
// // // // 			out << it->first << "\t" << it->second << std::endl ;
// // // // 		}
// // // // 
// // // // 		out << std::endl << "average edge degrees" << std::endl ;
// // // // 		for(std::map<std::pair<int, int>, double>::iterator it=ave_edeg_.begin(); it!=ave_edeg_.end(); ++it) {
// // // // 			it->second/=n ;
// // // // 			out << it->first.first <<"-" << it->first.second << "\t" << it->second << std::endl ;	
// // // // 		}
// // // // 		out.close() ;
// // // // 
// // // // 		sprintf(str, "%s%s_r%f_hist_ave.txt", fname.c_str(), method.c_str(), radius) ;
// // // // 		out.open(str) ;
// // // // 		for(int i=0; i<edge_hist_.size(); ++i) {
// // // // 			out << i << "\t" << edge_hist_[i]/n << std::endl ;
// // // // 		}
// // // // 		out.close() ;
// // // // 
// // // // 		show_edge_hist_ = TRUE ;
// // // // 		glut_viewer_redraw() ;
// // // // 		sprintf(str, "%s%s_r%f_hist_ave.png", fname.c_str(), method.c_str(), radius) ;
// // // // 		snapshot(str) ;
// // // // 	}
// // // // 
// // // 	void DelaunayIO::save_points(const std::string& fname) {
// // // 		std::vector<Delaunay::Vertex_handle>& vertices = delaunay_->vertices() ;
// // // 		std::ofstream out(fname.c_str()) ;
// // // 
// // // //		out << delaunay_->sample_radius() << std::endl ;
// // // 		out.precision(20) ;
// // // 		out << vertices.size() << std::endl ;
// // // 		for(unsigned int i=0; i<vertices.size(); ++i) {
// // // 			out << vertices[i]->point() << std::endl ;
// // // 		}
// // // 		out.close() ;
// // // 	}
// // // 
// // 	void DelaunayIO::save_vertex_degree(const std::string& fname) {
// // 		std::vector<Delaunay::Vertex_handle>& vertices = delaunay_->vertices() ;
// // 		std::map<int, double> vdeg ;
// // 		std::map<std::pair<int, int>, double> edeg ;
// // 		int vtotal = 0 ;
// // 		int etotal = 0 ;
// // 
// // 		for(unsigned int i=0; i<vertices.size(); ++i) {
// // 			//if(on_convex_hull(all_vertices_[i]) )
// // 			//	continue ;
// // 
// // 			if(delaunay_->is_primary(vertices[i])) {
// // 				int deg = delaunay_->degree(vertices[i]) ;
// // 				if(vdeg.find(deg) != vdeg.end()) {
// // 					vdeg[deg] += 1 ;
// // 				}
// // 				else {
// // 					vdeg[deg] = 1 ;
// // 				}
// // 
// // 				Delaunay::Vertex_circulator cir = delaunay_->incident_vertices(vertices[i]) ;
// // 				do {
// // 					int idx2 = cir->index ;// delaunay_->is_primary(cir) ? cir->index : cir->domain ;
// // 					int d2 = delaunay_->degree(vertices[idx2]) ;
// // 					std::pair<int, int> eij(deg, d2) ;
// // 					if(d2 < deg) eij = std::pair<int, int>(d2, deg) ;
// // 					if(edeg.find(eij) != edeg.end()) {
// // 						edeg[eij] += 1 ;
// // 					}
// // 					else {
// // 						edeg[eij] = 1 ;
// // 					}
// // 					etotal ++ ;
// // 					++cir ;
// // 				} while(cir!=delaunay_->incident_vertices(vertices[i])) ;
// // 				vtotal ++ ;
// // 			}
// // 		}
// // 
// // //		std::ofstream out(fname.c_str()) ;
// // //		out << "vertex degrees" << std::endl ;
// // 		for(std::map<int, double>::iterator it=vdeg.begin(); it!=vdeg.end(); ++it) {
// // 			it->second/=vtotal ;
// // //			out << it->first << "\t" << it->second << std::endl ;
// // 
// // 			if(ave_vdeg_.find(it->first)!=ave_vdeg_.end()) {
// // 				ave_vdeg_[it->first] += it->second ;
// // 			}
// // 			else {
// // 				ave_vdeg_[it->first] = it->second ;
// // 			}
// // 		}
// // 
// // //		out << std::endl << "edge degrees" << std::endl ;
// // 		for(std::map<std::pair<int, int>, double>::iterator it=edeg.begin(); it!=edeg.end(); ++it) {
// // 			it->second/=etotal ;
// // //			out << it->first.first <<"-" << it->first.second << "\t" << it->second << std::endl ;	
// // 
// // 			if(ave_edeg_.find(it->first)!=ave_edeg_.end()) {
// // 				ave_edeg_[it->first] += it->second ;
// // 			}
// // 			else {
// // 				ave_edeg_[it->first] = it->second ;
// // 			}
// // 		}
// // 
// // //		out.close() ;
// // 	}
// // 
// 	static void min_max_edge_length(Delaunay* dt, double& minlen, double& maxlen) {
// 		std::vector<Delaunay::Vertex_handle>& vertices = dt->vertices() ;
// 		minlen = 1e10 ;
// 		maxlen = -1e10 ;
// 
// 		for(unsigned int i=0; i<vertices.size(); ++i) {
// 			if(dt->is_primary(vertices[i])) {
// 				Delaunay::Vertex_circulator cir = dt->incident_vertices(vertices[i]) ;
// 				do {
// 					double curlen = distance(to_geex(vertices[i]->point()), to_geex(cir->point())) ;
// 					if(curlen < minlen) { 
// 						minlen = curlen ;
// 					}
// 					if(curlen > maxlen) {
// 						maxlen = curlen ;
// 					}
// 					++cir  ;
// 				} while(cir!=dt->incident_vertices(vertices[i])) ;
// 			}
// 		}	
// 	}
// 
// // // // 	void DelaunayIO::save_edge_histogram(const std::string& fname) {
// // // // 		std::vector<Delaunay::Vertex_handle>& vertices = delaunay_->vertices() ;
// // // // 		double binmin, binmax ;
// // // // 		min_max_edge_length(delaunay_, binmin, binmax) ;
// // // // 		int ntotal = 0 ;
// // // // 		int bin_size = edge_hist_.size() ;
// // // // 		double step = (binmax-binmin)/(bin_size-1) ;
// // // // 		std::vector<double> edge_lens ;
// // // // 		edge_lens.assign(bin_size, 0) ;
// // // // 		std::map<int, std::vector<double> > val_edge_hist ;
// // // // 		for(int i=3; i<=10; ++i) {
// // // // 			val_edge_hist[i].assign(bin_size, 0) ;
// // // // 		}
// // // // 
// // // // 		for(unsigned int i=0; i<vertices.size(); ++i) {
// // // // 			if(delaunay_->is_primary(vertices[i])) {
// // // // 				Delaunay::Vertex_circulator cir = delaunay_->incident_vertices(vertices[i]) ;
// // // // 				int deg = delaunay_->degree(cir) ;
// // // // 				do {
// // // // 					double curlen = distance(to_geex(vertices[i]->point()), to_geex(cir->point())) ;
// // // // 					if(curlen < binmax) {  // exclude edges on the vonvex hull
// // // // 						int bin_id = (int)((curlen-binmin)/step) ;
// // // // //						edge_hist_[bin_id]++ ;
// // // // 						edge_lens[bin_id]++ ;
// // // // 						if(deg>=3 && deg<11)
// // // // 							val_edge_hist[deg][bin_id]++ ;
// // // // 					}
// // // // 					++cir  ;
// // // // 					ntotal ++ ;
// // // // 				} while(cir!=delaunay_->incident_vertices(vertices[i])) ;
// // // // 			}
// // // // 		}
// // // // 
// // // // 		std::cerr << "min edge " << binmin << ", max edge " << binmax << std::endl ;
// // // // 		//std::ofstream out(fname.c_str()) ;
// // // // 		for(unsigned int i=0; i<bin_size; ++i) {
// // // // 			edge_hist_[i] += edge_lens[i]/ntotal ;
// // // // 		//	out << i << "\t" << edge_lens[i]/ntotal << std::endl ;
// // // // 		}
// // // // 		for(unsigned int i=4; i<=10; ++i) {
// // // // 			for(unsigned int j=0; j<val_edge_hist_[i].size(); ++j) {
// // // // 				val_edge_hist_[i][j] += val_edge_hist[i][j]/ntotal ;
// // // // 			}
// // // // 		}
// // // // 		//out.close() ;
// // // // 
// // // // 	}
// // // // 
// // // 	void DelaunayIO::save_angle_histogram(const std::string& fname) {
// // // 		std::vector<Delaunay::Vertex_handle>& vertices = delaunay_->vertices() ;
// // // 		int ntotal = 0 ;
// // // 		std::vector<double> angle_hist ;
// // // 		angle_hist.assign(120, 0.0) ;
// // // 
// // // 		for(unsigned int i=0; i<vertices.size(); ++i) {
// // // 			if(delaunay_->is_primary(vertices[i])) {
// // // 				vec2 p = to_geex(vertices[i]->point()) ;
// // // 				Delaunay::Vertex_circulator cir = delaunay_->incident_vertices(vertices[i]) ;
// // // 				Delaunay::Vertex_circulator cir2 = cir++ ;
// // // 				do {
// // // 					vec2 e1 = to_geex(cir->point())-p ;
// // // 					vec2 e2 = to_geex(cir2->point())-p ;
// // // 					double angle = to_degree(acos(dot(e1, e2)/(e1.length()*e2.length()))) ;
// // // 					if(angle < 30) {
// // // 						std::cerr << "angle less than 30: " << angle << std::endl ;
// // // //						gx_assert(angle >= 30) ;
// // // 					}
// // // 					angle_hist[int(angle+0.5)] ++ ;
// // // 					++cir  ;
// // // 					++cir2 ;
// // // 					ntotal ++ ;
// // // 				} while(cir!=delaunay_->incident_vertices(vertices[i])) ;
// // // 			}
// // // 		}
// // // 
// // // //		std::ofstream out(fname.c_str()) ;
// // // 		for(unsigned i=0; i<angle_hist.size(); ++i) {
// // // 			angle_hist[i]/=ntotal ;
// // // //			out << i << "\t" << angle_hist[i] << std::endl ;
// // // 		}
// // // //		out.close() ;
// // // 	}
// // // 
// // 	void DelaunayIO::save_area_histogram(const std::string& fname) {
// // 		int ntri = 0 ;
// // 		std::vector<double> area_hist ;
// // 		std::vector<double> areas ;
// // 		area_hist.assign(301, 0) ;
// // 		double min_area=1e10, max_area=-1e10 ;
// // 		double ave_area, total_area=0 ;
// // 
// // 		FOR_EACH_FINITE_FACE_DT(Delaunay, delaunay_, f) {
// // 			if(!f->dual_outside) {
// // 				double area = triangle_area(to_geex(f->vertex(0)->point()), 
// // 											to_geex(f->vertex(1)->point()), 
// // 											to_geex(f->vertex(2)->point())) ;
// // 				if(area < min_area) {
// // 					min_area = area ;
// // 				}
// // 				if(area > max_area) {
// // 					max_area = area ;
// // 				}
// // 
// // 				areas.push_back(area) ;
// // 				total_area += area ;
// // 				ntri ++ ;
// // 			}
// // 		}
// // 
// // 		ave_area = total_area/ntri ;
// // 		//for(unsigned int i=0; i<ntri; ++i) {
// // 		//	int bin = 300*(areas[i]-min_area)/(max_area-min_area) ;
// // 		//	area_hist[bin] ++ ;
// // 		//}
// // 
// // 		std::cerr << "total area " << total_area << ", average area " << ave_area 
// // 			      << ", min area " << min_area << ", max_area " << max_area << std::endl ;
// // 		//std::ofstream out(fname.c_str()) ;
// // 		//for(unsigned i=0; i<area_hist.size(); ++i) {
// // 		//	area_hist[i]/=ntri ;
// // 		//	out << i << "\t" << area_hist[i] << std::endl ;
// // 		//}
// // 		//out.close() ;
// // 	}
// // 	
// 	void DelaunayIO::spetral_analysis(const std::string& fname) {
// 	}
// 
	bool DelaunayIO::snapshot(const std::string& filename) {
		FILE *out = fopen(filename.c_str(), "wb") ;
		Memory::pointer base_mem ;
		int bytes_per_pixel=3, image_size ;
		int w, h ;

		glut_viewer_get_screen_size(&w, &h) ;
		image_size = bytes_per_pixel * w * h ;
		base_mem = new Memory::byte[image_size] ;
		Memory::clear(base_mem, image_size) ;

		glut_viewer_redraw() ; // test

		glPixelStorei(GL_PACK_ALIGNMENT, 1) ; 
		glPixelStorei(GL_PACK_ROW_LENGTH, w) ;
		glReadPixels( 0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, base_mem ) ;
		// Restore default
		glPixelStorei(GL_PACK_ROW_LENGTH, 0) ;

		// save png
		png_structp png_ptr = png_create_write_struct(
			PNG_LIBPNG_VER_STRING,
			(png_voidp) NULL, (png_error_ptr) NULL, (png_error_ptr) NULL 
		) ;
		if ( png_ptr == nil ) {
			fclose(out) ;
			return false ;
		} 

		png_infop info_ptr = png_create_info_struct ( png_ptr );
		png_init_io( png_ptr, out );


		png_byte png_color_encoding ;
		png_color_encoding = PNG_COLOR_TYPE_RGB ;

		png_set_IHDR ( 
			png_ptr, info_ptr, w, h,
			8, png_color_encoding, PNG_INTERLACE_NONE,
			PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE 
		);

		png_text comment;
		comment.compression = PNG_TEXT_COMPRESSION_NONE;
		comment.key = (char*)"Comment" ;
		comment.text = (char*)"Made with Graphite" ;
		png_set_text ( png_ptr, info_ptr, &comment, 1 );

		png_write_info ( png_ptr, info_ptr );
		for (int row = h-1; row >= 0 ; row-- ) {
			png_write_row ( 
				png_ptr, 
				base_mem + w * row * bytes_per_pixel
			);
		}
		png_write_end ( png_ptr, info_ptr );
		fflush ( out );
		png_destroy_info_struct (png_ptr, &info_ptr) ; //(png_infopp)NULL) ;
		png_destroy_write_struct ( &png_ptr, (png_infopp)NULL );

		fclose(out) ;
		delete [] base_mem ;

		return true ;
	}

// 	void DelaunayIO::analyze_pointsets(const std::string& path) {
// //		char str[1024] ;
// //		int  npoints = 0 ;
// //		std::vector<vec2> points ;
// //		std::list<std::string> fnames ;
// //
// //		if(FolderExists(path)) {
// //			FileList(path, fnames) ;
// //		} else if(FileExists(path)) {
// //			fnames.push_back(path) ;
// //		} else { // 
// //			std::cerr << path.c_str() << " is not a folder or file... " << std::endl ;
// //			return ;
// //		}
// //
// //		ave_vdeg_.clear() ;
// //		ave_edeg_.clear() ;
// //		edge_hist_.assign(256, 0) ;
// //		area_hist_.assign(300, 0) ;
// //		angle_hist_.assign(181, 0) ;
// //		show_edge_hist_ = FALSE ;
// //
// //		val_edge_hist_.clear() ;
// //		for(int i=3; i<=10; ++i) {
// //			val_edge_hist_[i].assign(256, 0) ;
// //		}
// //
// //		int nsets = 0 ;
// //		std::list<std::string>::iterator file;
// //		for (file = fnames.begin(); file != fnames.end(); ++file, ++nsets) {
// ////			sprintf(str, "%s_%d.txt", path.c_str(), i) ;
// //			std::cerr << "processing " << file->c_str() << std::endl ;
// //			points.clear() ;
// //			int npts ;
// //
// //			std::ifstream in(file->c_str()) ;
// //			//in >> delaunay_->sample_radius() ;
// //			in >> npts ;
// //			for(int j=0; j<npts; ++j) {
// //				vec2 p ;
// //				in >>p.x >> p.y;
// //				points.push_back(p) ;
// //			}
// //			in.close() ;
// //
// //			delaunay_->set_vertices(points) ;	
// //
// //			// save vertex degree 
// //			sprintf(str, "%s_deg_%d.txt", path.c_str(), nsets) ;
// //			save_vertex_degree(str) ;
// //
// //			// save edge histogram
// //			sprintf(str, "%s_edge_%d.txt", path.c_str(), nsets) ;
// //			save_edge_histogram(str) ;
// //
// //			//// save angle histogram
// //			//sprintf(str, "%s_angle_%d.txt", path.c_str(), nsets) ;
// //			//save_angle_histogram(str) ;
// //
// //			// save area histogram
// //			sprintf(str, "%s_area_%d.txt", path.c_str(), nsets) ;
// //			save_area_histogram(str) ;
// //
// //			npoints += npts ;
// //			double density = npts*M_PI*SQR(delaunay_->sample_radius()/2.0) ;
// //			std::cerr << "packing density " << density << std::endl ;
// //		}
// //		double density = (double)npoints/nsets*M_PI*SQR(delaunay_->sample_radius()/2.0) ;
// //		std::cerr << "average packing density " << density << std::endl ;
// //
// //		sprintf(str, "%s_deg_ave.txt", path.c_str()) ;
// //		std::ofstream out(str) ;
// //		out << "average number of points" << std::endl ;
// //		out << npoints / nsets << std::endl ;
// //		out << "average vertex degrees" << std::endl ;
// //		for(std::map<int, double>::iterator it=ave_vdeg_.begin(); it!=ave_vdeg_.end(); ++it) {
// //			it->second/=nsets ;
// //			out << it->first << "\t" << it->second << std::endl ;
// //		}
// //
// //		out << std::endl << "average edge degrees" << std::endl ;
// //		for(std::map<std::pair<int, int>, double>::iterator it=ave_edeg_.begin(); it!=ave_edeg_.end(); ++it) {
// //			it->second/=nsets ;
// //			out << it->first.first <<"-" << it->first.second << "\t" << it->second << std::endl ;	
// //		}
// //		out.close() ;
// //
// //		sprintf(str, "%s_hist_ave.txt", path.c_str()) ;
// //		out.open(str) ;
// //		for(int i=0; i<edge_hist_.size(); ++i) {
// //			out << i << "\t" << edge_hist_[i]/nsets << std::endl ;
// //		}
// //		out.close() ;
// //
// //		for(int i=4; i<=10; ++i) {
// //			sprintf(str, "%s_vertex_edge_hist_ave_%d.txt", path.c_str(), i) ;
// //			out.open(str) ;
// //			for(int j=0; j<val_edge_hist_[i].size(); ++j) {
// //				out << j << "\t" <<val_edge_hist_[i][j]/nsets << std::endl ;
// //			}
// //			out.close() ;
// //		}
// //
// //		show_edge_hist_ = TRUE ;
// //		glut_viewer_redraw() ;
// //		sprintf(str, "%s_hist_ave.png", path.c_str()) ;
// //		snapshot(str) ;
// 	}
}
