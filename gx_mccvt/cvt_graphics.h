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


#ifndef __GPURVD_GRAPHICS__
#define __GPURVD_GRAPHICS__

#include "delaunay_cvt.h"
#include <Geex/image/colormap.h>
#include <Geex/graphics/opengl.h>
//dxy add: for snapshot
#include "delaunay_io.h"
//

namespace Geex {
    
	enum SurfaceMode { None, Plain, Chart, Cells, CellType, Density, Energy, GeoDist } ;
	enum PrimalMode {PRI_PLAIN, PRI_GAP, PRI_OBTUSE} ;
	//dxy
	enum DrawEnergyMode {CVT_True, CVT_Approx, CVT_True_Minus_Approx, DIR, Cell_Area}; 

	//

    class CVTGraphics {
    public:
		CVTGraphics(DelaunayCVT* arvd) ;
  //      CVTGraphics(DelaunayPD* sampler, SurfacePD*  surfpd) ;
		//Delaunay* delaunay() { return sampler_->delaunay() ; }
		void save_viewpoint(const std::string& filename) ;
		void load_viewpoint(const std::string& filename) ;

        void draw() ;

        GLboolean& show_domain_mesh() { return show_domain_mesh_ ; }
        GLboolean& show_domain() { return show_domain_ ; }
        GLfloat&   vertices_size() { return vertices_size_ ; }
		GLboolean& show_feature() { return show_feature_ ; }
		//GLboolean& show_power_sphere() { return show_power_sphere_ ; }
		//GLboolean& show_boundary_grid() { return show_boundary_grid_ ; }		
  //      GLboolean& show_mesh() { return show_mesh_ ; }
        GLboolean& colorize() { return colorize_ ; }
		GLboolean& shiny() { return shiny_ ; }
  //      GLboolean& show_primal() { return show_primal_ ; }
		//GLboolean& show_primal_mesh() { return show_primal_mesh_ ; }
		GLboolean& show_rvd() { return show_rvd_ ; }
		GLboolean& show_rvd_mesh() { return show_rvd_mesh_ ; }
		GLboolean& show_rvd_vertices() { return show_rvd_vertices_ ; }
		//GLboolean& show_cluster() { return show_cluster_ ; }
		//GLboolean& show_active_triangles() { return show_active_triangles_ ; }

		///dxy add
		GLboolean& show_domain_facet_field()  { return show_domain_facet_field_; }
		GLboolean& show_domain_field() { return show_domain_field_ ; }
		GLboolean& show_seed_field() { return show_seed_field_; }
		GLboolean& show_seed_match_field() { return show_seed_match_field_; }
		GLuint&	   max_match_dir_order() { return max_match_dir_order_; }
		GLint&	   max_match_number()    { return max_match_number_; }
		GLboolean& show_normal()     { return show_normal_; }
		GLboolean& show_short_edges() { return show_short_edges_; }
		GLboolean& show_long_edges() { return show_long_edges_; }
		//
		GLboolean& show_QD_topo_operations() { return show_QD_topo_operations_; }
		//
		GLboolean& show_seed_project_triangle() { return show_seed_project_triangle_; }
		//
		GLboolean& show_lloyd_grad() { return show_lloyd_grad_; }
		GLboolean& show_direction_grad() { return show_direction_grad_; }
		GLfloat& grad_magnification() { return grad_magnification_; }
		//
		GLboolean& show_stretched_long_edges() { return show_stretched_long_edges_; }
		GLboolean& show_stretched_short_edges() { return show_stretched_short_edges_; }
		GLfloat& show_stretched_edges_rate() { return show_stretched_edges_rate_; }
		//
		GLboolean& show_valid_topo_operation() { return show_valid_topo_operation_; }
		GLboolean& show_edge_flip() { return show_edge_flip_; }
		GLboolean& show_edge_collapse() { return show_edge_collapse_; }
		GLboolean& show_vertex_split() { return show_vertex_split_; }
		int& show_topo_operation_deltaR() { return topo_operation_deltaR_; }
		//
		GLboolean& show_singular_rvd_vertices() { return show_singular_rvd_vertices_; }
		//
		GLboolean& show_sharp_edges() { return show_sharp_edges_; }
		//
		GLenum	&draw_energy_mode()  { return draw_energy_mode_; }
		GLboolean &show_cell_energy()  { return show_cell_energy_; }
		double &cell_energy_truncate_value() { return cell_energy_truncate_value_; }
		///

		//GLboolean& show_vertex_gaps() { return show_vertex_gaps_ ; }
		//GLboolean& show_gap_domain() { return show_gap_domain_ ; }
		GLenum& surface_mode()       { return surface_mode_ ; }
		//GLenum& primal_mode()       { return primal_mode_ ; }
		//
  //      float& slice_x() { return slice_x_ ; }
  //      float& slice_y() { return slice_y_ ; }
  //      float& slice_z() { return slice_z_ ; }
  //      float& selection() { return selection_ ; }
		//GLboolean& set_cull_back() {return show_cull_back_; }
		//bool slice_visible(const vec3& p) {	return p.x<=slice_x_ && p.y<=slice_y_ && p.z<=slice_z_ ; }

		//GLboolean& show_slivers() { return show_slivers_ ; } 
		//GLboolean& show_primal_cells() { return show_primal_rvd_facets_ ; } 
		//GLboolean& show_gap_cells()    { return show_gap_rvd_facets_ ; }

		//GLboolean& show_union_balls()  { return show_union_balls_ ; }

		//GLboolean& show_subd()      { return show_subd_ ; }
		//GLboolean& show_subd_mesh() { return show_subd_mesh_ ; }
		//GLboolean& show_approx_rvd()    { return show_approx_rvd_ ; }
		//GLboolean& show_clipped_rvd()   { return show_clipped_rvd_ ; }
		GLboolean& show_primal()        { return show_primal_ ; }
		GLboolean& show_obtuse()        { return show_obtuse_ ; }
		GLboolean& show_primal_features() { return show_primal_features_; }
		//GLboolean& show_edge_color()    { return show_edge_color_ ; }
		//GLboolean& show_facet_gaps(){ return show_facet_gaps_ ; }
		//GLboolean& show_facet_gaps2(){ return show_facet_gaps2_ ; }
		//GLboolean& show_facet_arcgons(){ return show_facet_arcgons_ ; }
		//GLboolean& show_fragments()  { return show_fragments_ ; }
		//GLboolean& show_gap_spheres() { return show_gap_spheres_ ; }

		//dxy add
		void save_snapshot(const std::string& filename);
		//

    protected:
        void draw_domain_mesh() ;
        void draw_domain() ;
		void draw_domain_density() ;
		//void draw_domain_geodist() ;
		//void draw_iso_lines() ;
        void draw_vertices() ;
		void draw_domain_vertices() ;
		//void draw_cluster() ;
		//void draw_active_triangles() ;
		//void draw_features() ;
		//void draw_boundary_grid() ;
		//void draw_power_sphere() ;
		void draw_primal() ;
		void draw_obtuse_triangles() ;
        void draw_rvd() ;
        void draw_rvd_mesh() ;
		void draw_feature() ;
		void draw_chart() ;
  //      void draw_primal() ;
  //      void draw_primal_mesh() ;
        void create_colormap_texture() ;
		//void draw_vertex_gaps() ;
		//void draw_gap_domain() ; // draw domain that contains gaps
		//void draw_bad_triangles() ;
		//void draw_union_balls() ;
		//void draw_subd() ;
		//void draw_subd_mesh() ;
		//void draw_approx_rvd_facets() ;
		//void draw_approx_rvd() ;
		//void draw_clipped_rvd() ;
		//void draw_clipped_facets() ;
		//void draw_facet_gaps() ;
		//void draw_facet_gaps2() ;
		//void draw_facet_arcgons() ;
		//void draw_fragments() ;
		//void draw_gap_spheres() ;
		//void draw_corners() ;
		//void draw_edge_color() ;

		///dxy add
		void draw_domain_facet_field();
		void draw_domain_field();
		void draw_seed_field();
		void draw_seed_match_field();
		//
		void draw_QD_topo_operations(); //Quad Dominant mesh topo operations
		//
		void draw_seed_project_triangle();
		//
		void draw_short_edges();
		void draw_long_edges();
		//
		void draw_lloyd_grad();
		void draw_direction_grad();
		//
		void draw_cell_energy(GLenum energy_mode);
		//
		void draw_primal_map();
		void draw_sharp_edges();
		//
		void draw_stretched_long_edges();
		void draw_stretched_short_edges();
		//
		void draw_primal_map_vertices(bool show_singularity = false);
		//
		void draw_valid_topo_operation(int operator_code, int deltaR_threshold = 0, bool isolation = true);
		//
		void draw_primal_features();
		///

		//// volumetric
		//void draw_slivers() ;
		//void draw_primal_cells() ;
		//void draw_primal_cell(const vec3& p0, const vec3& p2, const vec3& p3, const vec3& p4, bool mesh=false) ;
		//void draw_gap_cells() ;
		//void draw_gap_power_sphere() ;
		
    private:
//		LocalRVD*  lrvd_ ;
		DelaunayCVT* cvt_ ;
  //      DelaunayPD* sampler_ ;
		//SurfacePD*  surfpd_ ;
        GLboolean show_domain_mesh_ ;
        GLboolean show_domain_ ;
		GLboolean show_feature_ ;
        GLfloat   vertices_size_ ;
		GLboolean show_rvd_vertices_ ;
	//	GLboolean show_active_triangles_ ;
	//	GLboolean show_cluster_ ;
	//	GLboolean show_power_sphere_ ;
	//	GLboolean show_boundary_grid_ ;
		GLboolean show_rvd_ ;
        GLboolean show_rvd_mesh_ ;
        GLboolean show_rvd_submesh_ ;
        GLboolean show_primal_ ;
        GLboolean show_primal_mesh_ ;
		GLboolean show_primal_features_;

		///dxy add
		GLboolean show_domain_facet_field_;
		GLboolean show_domain_field_;
		GLboolean show_seed_field_;
		GLboolean show_seed_match_field_;
		GLuint    max_match_dir_order_;
		GLint	  max_match_number_;
		GLboolean show_normal_;
		GLboolean show_short_edges_;
		GLboolean show_long_edges_;
		//
		GLboolean show_QD_topo_operations_;
		//
		GLboolean show_seed_project_triangle_;
		//
		GLboolean show_lloyd_grad_;
		GLboolean show_direction_grad_;
		GLfloat grad_magnification_;
		//
		GLboolean show_stretched_long_edges_;
		GLboolean show_stretched_short_edges_;
		GLfloat show_stretched_edges_rate_;
		//
		GLboolean show_valid_topo_operation_;
		GLboolean show_edge_flip_;
		GLboolean show_edge_collapse_;
		GLboolean show_vertex_split_;
		int topo_operation_deltaR_;
		//
		GLboolean show_singular_rvd_vertices_;
		//
		GLboolean show_sharp_edges_;
		//
		GLenum draw_energy_mode_;
		GLboolean show_cell_energy_;
		double cell_energy_truncate_value_;
		// snapshot
		DelaunayIO IO_;
		///

    //    GLfloat   shrink_ ;
    //    GLboolean show_non_manifold_ ;
        GLuint    colormap_texture_ ;
        GLboolean shiny_ ;
    //    GLboolean show_rvd_facets_ ;
	//	GLboolean show_vertex_gaps_ ;
	//	GLboolean show_gap_domain_ ;
		GLenum    surface_mode_ ;
	//	GLenum    primal_mode_ ;
	//	GLboolean show_subd_ ;
	//	GLboolean show_subd_mesh_ ;
	//	GLboolean show_approx_rvd_ ;
	//	GLboolean show_clipped_rvd_ ;
	//	GLboolean show_facet_gaps_ ;
	//	GLboolean show_facet_gaps2_ ;
	//	GLboolean show_facet_arcgons_ ;
	//	GLboolean show_fragments_ ;
	//	GLboolean show_gap_spheres_ ;
		GLboolean show_obtuse_ ;
	//	GLboolean show_edge_color_ ;

		////---------------------------------- volumetric -----------------------------------------
		//GLboolean show_primal_rvd_facets_ ;
		//GLboolean show_slivers_ ;
		//GLboolean show_gap_tets_ ;
		//GLboolean show_gap_vol_ ;
		//GLboolean show_clipped_vd_ ;
		//GLboolean show_gap_rvd_facets_ ;

        GLboolean show_mesh_ ;
        GLboolean colorize_ ;
  //      float slice_x_ ;
  //      float slice_y_ ;
  //      float slice_z_ ;
  //      float selection_ ;
  //      bool non_convex_ ;

		//GLboolean show_union_balls_ ;
		Colormap  colormap_ ;

        //---------------------------------------------------------------------------------------
  //      GLboolean show_surface_clipping_ ;
		//GLboolean show_volume_clipping_ ;
		//GLboolean show_remesh_ ;
		//GLboolean show_volume_ ;
		//GLboolean show_retet_ ;
		//GLboolean show_cull_back_;
    } ;

} 

#endif
