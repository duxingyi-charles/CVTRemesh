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

#include "cvt_graphics.h"
#include "cvt_geometry.h"
#include <Geex/combinatorics/map_geometry.h>
#include <glut_viewer/glut_viewer.h>
#include <Geex/basics/stopwatch.h>
#include <Geex/graphics/geexob.h>
//dxy add
//#include "delaunay_io.h"
//


namespace Geex {

	static const double c1 = 0.35 ;
	static const double c2 = 0.5 ;
	static const double c3 = 1.0 ;

	static double color_table[12][3] = 
	{
		{c3, c2, c2},
		{c2, c3, c2},
		{c2, c2, c3},
		{c2, c3, c3},
		{c3, c2, c3},
		{c3, c3, c2},

		{c1, c2, c2},
		{c2, c1, c2},
		{c2, c2, c1},
		{c2, c1, c1},
		{c1, c2, c1},
		{c1, c1, c2}

	} ;

	static int random_color_index_ = 0 ;

	static float white[4] = {
		0.8f, 0.8f, 0.3f, 1.0f
	} ;

	static void gl_random_color() {
		glColor3dv(color_table[random_color_index_]);
		random_color_index_ = (random_color_index_ + 1) % 12 ;
	}

	static void gl_random_color(int index) {
		random_color_index_ = index % 12 ;
		gl_random_color() ;
	}

	static void gl_randomize_colors(int index = 0) {
		random_color_index_ = (index % 12) ;
	}

	static void gl_vertex_color(int degree) {
		if(degree==6) {
			//glColor3f(0.85, 0.95, 0.5) ;
			//glColor3f(213/255.0, 222/255.0, 157/255.0) ;
			glColor3f(1.0, 1.0, 1.0);  ///dxy change
		}
		else if(degree==7) {
			//glColor3f(1.0, 0.5, 0.25) ;
			glColor3f(125/255.0, 176/255.0, 223/255.0) ;
		}
		else if(degree==5) {
			//glColor3f(0.25, 0.5, 1.0) ;
			glColor3f(192/255.0, 121/255.0, 165/255.0) ;
		}
		else if(degree>7) {
			//glColor3f(0.5, 0, 0.25) ;
			glColor3f(98/255.0, 57/255.0, 115/255.0) ;
		}
		else {
			glColor3f(0, 0.25, 0.5) ;
			glColor3f(41/255.0, 118/255.0, 102/255.0) ;
		}
	}

	//dxy add: hsv to rgb
	//hsv ([0,360], [0,1], [0,1]),  rgb ([0,1]...)
	static void HSVtoRGB(real& fR, real& fG, real& fB, real& fH, real& fS, real& fV) {
		real fC = fV * fS; // Chroma
		real fHPrime = fmod(fH / 60.0, 6);
		real fX = fC * (1 - fabs(fmod(fHPrime, 2) - 1));
		real fM = fV - fC;

		if(0 <= fHPrime && fHPrime < 1) {
			fR = fC;
			fG = fX;
			fB = 0;
		} else if(1 <= fHPrime && fHPrime < 2) {
			fR = fX;
			fG = fC;
			fB = 0;
		} else if(2 <= fHPrime && fHPrime < 3) {
			fR = 0;
			fG = fC;
			fB = fX;
		} else if(3 <= fHPrime && fHPrime < 4) {
			fR = 0;
			fG = fX;
			fB = fC;
		} else if(4 <= fHPrime && fHPrime < 5) {
			fR = fX;
			fG = 0;
			fB = fC;
		} else if(5 <= fHPrime && fHPrime < 6) {
			fR = fC;
			fG = 0;
			fB = fX;
		} else {
			fR = 0;
			fG = 0;
			fB = 0;
		}

		fR += fM;
		fG += fM;
		fB += fM;
	}
	static void hsv2rgb(vec3& hsv, vec3& rgb) 
	{
		real& fR = rgb[0];
		real& fG = rgb[1];
		real& fB = rgb[2];
		real& fH = hsv[0];
		real& fS = hsv[1];
		real& fV = hsv[2];
		HSVtoRGB(fR, fG, fB, fH, fS, fV);
	}
	///dxy add end
	//dxy add end

	//CVTGraphics::CVTGraphics(DelaunayPD* sampler,  SurfacePD*  surfpd) 
	//: sampler_(sampler), cvt_(surfpd) {
	CVTGraphics::CVTGraphics(DelaunayCVT* arvd) 
		:cvt_(arvd) {
			show_domain_ = GL_TRUE ;
			show_domain_mesh_ = GL_FALSE ;
			show_feature_ = GL_FALSE ;
			vertices_size_ = 0.05 ;
			show_rvd_vertices_ = GL_FALSE ;
			show_mesh_ = GL_TRUE ;
			colorize_ = GL_FALSE ;
			show_primal_ = GL_FALSE ;
			show_primal_mesh_ = GL_FALSE ;
			show_primal_features_ = GL_FALSE;
			show_rvd_mesh_ = GL_FALSE ;
			show_rvd_ = GL_FALSE ;
			shiny_ = GL_FALSE ;
			double xmin, ymin, zmin, xmax, ymax, zmax ;
			colormap_texture_ = -1 ;	
			create_colormap_texture() ;
			show_obtuse_ = GL_FALSE ;
			surface_mode_ = Plain ;


			///dxy add
			show_domain_facet_field_ = GL_FALSE;
			show_domain_field_ = GL_FALSE;
			show_seed_field_ = GL_FALSE;
			show_seed_match_field_ = GL_FALSE;
			show_QD_topo_operations_ = GL_FALSE;
			max_match_dir_order_ = 0;
			max_match_number_ = 6;
			show_normal_ = GL_FALSE;
			show_seed_project_triangle_ = GL_FALSE;
			show_short_edges_ = GL_FALSE;
			show_long_edges_ = GL_FALSE;
			show_lloyd_grad_ = GL_FALSE;
			show_direction_grad_ = GL_FALSE;
			grad_magnification_ = 1.0;
			show_stretched_long_edges_ = GL_FALSE;
			show_stretched_short_edges_ = GL_FALSE;
			show_stretched_edges_rate_ = 0.0;
			show_valid_topo_operation_ = GL_FALSE;
			topo_operation_deltaR_ = 0;
			show_singular_rvd_vertices_ = GL_FALSE;
			show_sharp_edges_ = GL_FALSE;
			draw_energy_mode_ = CVT_True;
			show_cell_energy_ = GL_FALSE;
			cell_energy_truncate_value_ = 0.01;
			IO_ = DelaunayIO();
			///
	}

	inline void draw_triangle(const vec3& v1, const vec3& v2, const vec3& v3) {
		glVertex(v1) ;
		glVertex(v2) ;
		glVertex(v3) ;
	}

	inline void draw_triangle_with_normal(
		const vec3& v1, const vec3& v2, const vec3& v3
		) {
			glNormal(cross(v2-v1, v3-v1)) ;
			glVertex(v1) ;
			glVertex(v2) ;
			glVertex(v3) ;
	}

	class draw_triangle_with_normal_and_colors {
	public:
		draw_triangle_with_normal_and_colors(
			double wmin, double wmax, CVTGraphics *graphics
			) : bias_(wmin), scale_(1.0/(wmax - wmin)), graphics_(graphics) {
		}
		void operator() (
			const TopoPolyVertexEdge& v1, 
			const TopoPolyVertexEdge& v2, 
			const TopoPolyVertexEdge& v3
			) const {
				//	if(graphics_->slice_visible(1.0/3.0*(v1+v2+v3))) {
				glNormal(cross(v2-v1, v3-v1)) ;
				send_vertex(v1) ;
				send_vertex(v2) ;
				send_vertex(v3) ;
				//	}
		}
	protected:
		void send_vertex(const TopoPolyVertexEdge& v) const {
			double w = scale_ * (v.w - bias_) ;
			glTexCoord1d(w) ;
			glVertex(v) ;
		}
	private:
		double bias_ ;
		double scale_ ;
		CVTGraphics *graphics_ ; 
	} ;

	static float spec[4] = {
		0.7f, 0.7f, 0.7f, 1.0f
	} ;

	static float shininess = 20.0f ;

	static void color_ramp(
		GLfloat* cmap, unsigned int i1, vec3 c1, unsigned int i2, vec3 c2
		) {
			for(unsigned int i=i1; i<=i2; i++) {
				double s = double(i-i1) / double(i2-i1) ;
				vec3 c = s*c2 + (1.0 - s)*c1 ;
				cmap[3*i]   = GLfloat(c.x) ;
				cmap[3*i+1] = GLfloat(c.y) ;
				cmap[3*i+2] = GLfloat(c.z) ;
			}
	}

	void CVTGraphics::create_colormap_texture() {
		GLfloat RGB[256*3] ;

		color_ramp(RGB, 0, vec3(0.0, 0.0, 1.0), 127, vec3(0.7, 0.0, 0.0)) ;
		color_ramp(RGB, 127, vec3(0.7, 0.0, 0.0), 255, vec3(1.0, 1.0, 0.0)) ;

		glGenTextures(1, &colormap_texture_) ;
		glEnable(GL_TEXTURE_1D) ;
		glBindTexture(GL_TEXTURE_1D, colormap_texture_) ;
		gluBuild1DMipmaps(GL_TEXTURE_1D, GL_RGB, 256, GL_RGB, GL_FLOAT, RGB) ; 
		glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);            
		glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR) ;
		glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR) ;
		glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP) ;
		glDisable(GL_TEXTURE_1D) ;
	}

	void CVTGraphics::draw() {
		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE) ;
		//
		if (Geexob::program_objects_supported()) {
			glUseProgramObjectARB(0) ;
		}

		if(shiny_) {
			white[0] = 0.8f ;
			white[1] = 0.8f ;
			white[2] = 0.3f ;
			white[3] = 1.0f ;
			spec[0] = 0.7f ;
			spec[1] = 0.7f ;
			spec[2] = 0.7f ;
			spec[3] = 1.0f ;
			shininess = 20.0f ;
		} else {
			white[0] = 1.0f ;
			white[1] = 1.0f ;
			white[2] = 1.0f ;
			white[3] = 1.0f ;
			spec[0] = 0.0f ;
			spec[1] = 0.0f ;
			spec[2] = 0.0f ;
			spec[3] = 1.0f ;
			shininess = 0.0f ;
		}

		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,   spec) ;
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &shininess) ;

		glEnable(GL_COLOR_MATERIAL) ;
		glDisable(GL_CULL_FACE) ;

		if(show_domain_mesh_) {
			draw_domain_mesh() ;
		}

		if(show_rvd_) {
			draw_rvd() ;
			if (!colorize_)
			{
				cvt_->update_valences(true);
			}
		}
		if(show_rvd_mesh_) {
			draw_rvd_mesh() ;
		}
		///dxy add
		if (show_cell_energy_)
		{
			draw_cell_energy(draw_energy_mode_);
		}

		///

		if(show_domain_) {
			glColor4fv(white) ;
			draw_domain() ;
		}
		//dxy comment
		// 		if(show_obtuse_) {
		// 			draw_obtuse_triangles() ;
		// 		}

		if(show_primal_) {
			//draw_primal() ;
			draw_primal_map(); //dxy
		}
		if (show_primal_features_) {
			draw_primal_features();
		}
		if (show_sharp_edges_)
		{
			draw_sharp_edges();
		}
		if(vertices_size_ != 0.0) {
			//draw_vertices() ;	
			draw_primal_map_vertices(show_singular_rvd_vertices_); //dxy
		}
		if(show_rvd_vertices_) {
			draw_domain_vertices() ;
		}
		if(show_feature_) {
			draw_feature() ;
		}

		///dxy add
		if (cvt_->field_loaded())
		{
			if (show_domain_facet_field_)
			{
				draw_domain_facet_field();
			}

			if (show_domain_field_) {
				draw_domain_field();
			}

			if (show_seed_field_) {
				if (show_seed_match_field_) {
					draw_seed_match_field();
				} else {
					draw_seed_field();
				}
			}

			if (show_QD_topo_operations_) {
				draw_QD_topo_operations();
			}
		}

		if (show_seed_project_triangle_)
		{
			draw_seed_project_triangle();
		}

		if (show_short_edges_)
		{
			draw_short_edges();
		}

		if (show_long_edges_)
		{
			draw_long_edges();
		}

		if (show_lloyd_grad_)
		{
			draw_lloyd_grad();
		}

		if (show_direction_grad_)
		{
			draw_direction_grad();
		}

		if (show_stretched_long_edges_ && show_stretched_edges_rate_ > 0)
		{
			draw_stretched_long_edges();
		}

		if (show_stretched_short_edges_ && show_stretched_edges_rate_ > 0)
		{
			draw_stretched_short_edges();
		}

		if (show_valid_topo_operation_)
		{
			int c = (show_edge_flip_?Edge_Flip:Null_Operator)
				| (show_edge_collapse_?Edge_Collapse:Null_Operator)
				| (show_vertex_split_?Vertex_Split:Null_Operator);
			draw_valid_topo_operation(c, -topo_operation_deltaR_);
		}
		///
	}



	void CVTGraphics::draw_domain_mesh() {
		//glDisable(GL_LIGHTING) ;  //dxy change
		glDepthRange(0.0, 0.999) ;
		glLineWidth(1) ;
		glColor3f(/*0.0, 0.5, 1.0*/0,0,0) ; //dxy change
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) ;
		glBegin(GL_TRIANGLES) ;
		Map* m = cvt_->map();
		FOR_EACH_FACET(Map, m, fit) {
			glNormal(fit->normal()) ;
			glVertex(fit->halfedge()->vertex()->point()) ;
			glVertex(fit->halfedge()->next()->vertex()->point()) ;
			glVertex(fit->halfedge()->prev()->vertex()->point()) ;
		}
		//sampler_->boundary()->for_each_triangle(draw_triangle) ;
		glEnd() ;
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;
		glDepthRange(0.0, 1.0) ;
	}

	void CVTGraphics::draw_domain_vertices() {
		std::vector<int>& vcmap = cvt_->vert_cell_map() ;
		Map* m = cvt_->map() ;
		MapVertexProperty<int>& vindex = cvt_->vindex() ;


		glDisable(GL_LIGHTING) ;
		glDepthRange(0.0, 0.998) ;
		glHint(GL_POINT_SMOOTH_HINT, GL_NICEST) ;
		glEnable(GL_POINT_SMOOTH) ;
		glPointSize(vertices_size_ * 20.0f) ;
		glColor3f(0.0, 0.0, 0.0) ;
		glBegin(GL_POINTS) ;
		FOR_EACH_VERTEX(Map, m, v) {
			//gl_random_color(vcmap[vindex[v]]) ; //dxy remove: vcmap is not calculated in this program
			glVertex(v->point()) ;
		}
		glEnd() ;
		glDepthRange(0.0, 1.0) ;	
	}

	void CVTGraphics::draw_domain() {
		TopoPolyMesh* boundary = cvt_->boundary() ;

		if(surface_mode_ == None)
			return ;

		glEnable(GL_LIGHTING) ;

		if(surface_mode_ == Density) {
			draw_domain_density() ;
		}
		else if(surface_mode_ == Chart) {
			std::vector<std::vector<Map::Facet*>>& charts = cvt_->charts() ;
			glEnable(GL_LIGHTING) ;

			for(int i=0; i<charts.size(); ++i) {
				gl_random_color(i) ;
				glBegin(GL_TRIANGLES) ;
				for(int j=0; j<charts[i].size(); ++j) {
					Map::Facet* f = charts[i][j] ;
					Map::Halfedge* h = charts[i][j]->halfedge() ;
					glNormal(f->normal() ) ;
					do {
						glVertex(h->vertex()->point()) ;
						h = h->next() ;
					}
					while (h!=charts[i][j]->halfedge()) ;
				}
				glEnd() ;
			}

		}
		///dxy change
		else {
			Map *m = cvt_->map();
			glEnable(GL_LIGHTING);
			glBegin(GL_TRIANGLES);
			FOR_EACH_FACET(Map, m, fit) {
				const vec3& v1= fit->halfedge()->vertex()->point() ;
				const vec3& v2= fit->halfedge()->next()->vertex()->point() ;
				const vec3& v3= fit->halfedge()->prev()->vertex()->point() ;
				double amin, amax;
				min_max_angle(v1, v2, v3, amin, amax);

				if(show_obtuse_) {
					if(amax > cvt_->angle_max() && amin < cvt_->angle_max()) {
						glColor3f(0.949f, 0.577f, 0.629f) ;
					}
					else if(amin < cvt_->angle_min()) {
						glColor3f(0.609f, 0.838f, 0.962f) ;
					} 
					else if (amax > cvt_->angle_max()) {
						glColor3f(0.609f, 0.962, 0.609f) ;
					}
					else 
						glColor3f(1, 1, 1) ;
				}else
					glColor4fv(white) ;

				glNormal(cross(v2-v1, v3-v1)) ;
				glVertex(v1) ;
				glVertex(v2) ;
				glVertex(v3) ;
			}
			glEnd();
			glDisable(GL_LIGHTING);
		}
		///
		// 		else {
		// 
		// 			glBegin(GL_TRIANGLES) ;
		// 			glColor4fv(white) ;
		// 			cvt_->boundary()->for_each_triangle(draw_triangle_with_normal) ;
		// 			glEnd() ;
		// 		}
	}

	void CVTGraphics::draw_domain_density() {
		glDisable(GL_LIGHTING) ;
		glEnable(GL_TEXTURE_1D) ;
		glBindTexture(GL_TEXTURE_1D, colormap_texture_) ;
		glBegin(GL_TRIANGLES) ;
		double wmin = cvt_->boundary()->vertex(0).w ;
		double wmax = wmin ;
		for(unsigned int i=0; i<cvt_->boundary()->nb_vertices(); i++) {
			wmin = gx_min(wmin, cvt_->boundary()->vertex(i).w) ;
			wmax = gx_max(wmax, cvt_->boundary()->vertex(i).w) ;
		}
		glColor3f(1.0, 1.0, 1.0) ;
		cvt_->boundary()->for_each_triangle(
			draw_triangle_with_normal_and_colors(wmin, wmax, this)
			) ;
		glEnd() ;
		glDisable(GL_TEXTURE_1D) ;
	}

	void CVTGraphics::draw_vertices() {
		//		cvt_->update_primal() ;

		glDisable(GL_LIGHTING) ;
		glDepthRange(0.0, 0.998) ;
		glHint(GL_POINT_SMOOTH_HINT, GL_NICEST) ;
		glEnable(GL_POINT_SMOOTH) ;
		glPointSize(vertices_size_ * 20.0f) ;
		glColor3f(0.0, 0.0, 0.0) ;
		glBegin(GL_POINTS) ;
		//std::vector<SamplePoint>& samples = cvt_->samples() ;
		std::vector<vec3>& vertices = cvt_->vertices() ;
		std::vector<bool>& locked = cvt_->locked() ;
		for(unsigned int i=0; i<vertices.size(); i++) {
			if(locked[i])
				glColor3f(1.0, 0.0, 0.0) ;
			else if(colorize_) {
				//	gl_vertex_color(cvt_->valences()[i]) ;
			}
			else
				glColor3f(0.0, 0.0, 0.0) ;
			glVertex(vertices[i]) ;
		}
		glEnd() ;
		glDepthRange(0.0, 1.0) ;
	}

	///dxy add
	void CVTGraphics::draw_domain_facet_field() {
		glDisable(GL_LIGHTING);
		glDepthRange(0.0, 0.998);
		glLineWidth(2);
		glColor3f(1.0, 0.7, 0.2);
		double scale = 0.3 * sqrt(cvt_->total_area() / cvt_->map()->nb_vertices());
		glBegin(GL_LINES);
		for (int f=0; f<cvt_->boundary()->nb_facets(); ++f)
		{
			vec3 fcenter = cvt_->boundary()->facet_center(f);
			vec3 ffield = cvt_->boundary()->facet_field(f);
			glVertex(fcenter);
			glVertex(fcenter + scale * ffield);
		}

		if (show_normal_)
		{
			const std::vector<vec3> &fnorm = cvt_->boundary()->facet_normal();
			glColor3f(0.3, 0.3, 0.6);
			for (int f=0; f<cvt_->boundary()->nb_facets(); ++f)
			{
				vec3 fcenter = cvt_->boundary()->facet_center(f);
				glVertex(fcenter);
				glVertex(fcenter + scale * fnorm[f]);
			}
		}

		glEnd();
		glDepthRange(0.0, 1.0);
	}

	void CVTGraphics::draw_domain_field() {
		glDisable(GL_LIGHTING);
		glDepthRange(0.0, 0.998);
		glLineWidth(2);
		glColor3f(0.8, 0.5, 0.0);
		double scale = 0.3 * sqrt(cvt_->total_area() / cvt_->map()->nb_vertices());
		glBegin(GL_LINES);
		const std::vector<vec3>& vertices = cvt_->boundary()->original_vertices();
		const std::vector<vec3>& field = cvt_->boundary()->vertex_field();

		for (int i=0; i<vertices.size(); ++i)
		{
			glVertex(vertices[i]);
			glVertex(vertices[i] +  scale * field[i]);

		}
		if (show_normal_)
		{
			const std::vector<vec3> &vnorm = cvt_->boundary()->vertex_normal();
			glColor3f(0.3, 0.3, 0.7);
			for (int i=0; i<vertices.size(); ++i) {
				glVertex(vertices[i]);
				glVertex(vertices[i] + scale * vnorm[i]);
			}
		}

		glEnd();
		glDepthRange(0.0, 1.0);
	}	

	void CVTGraphics::draw_seed_field() {
		glDisable(GL_LIGHTING);
		glDepthRange(0.0, 0.998);
		glLineWidth(2);
		glColor3f(0.0, 0.8, 0.5);
		double scale = 0.3 * sqrt(cvt_->average_area());
		int rosy = cvt_->rot_symmetry();
		glBegin(GL_LINES);
		const std::vector<vec3>& vertices = cvt_->vertices();
		const std::vector<vec3>& field = cvt_->seed_field();
		const std::vector<vec3> &vnorm = cvt_->seed_normal();
		for (int i=0; i<vertices.size(); ++i)
		{
			for (int r=0; r<rosy; ++r)
			{
				glVertex(vertices[i]);
				glVertex(vertices[i] + scale * cvt_->vec_rotate(field[i], vnorm[i], r*2*M_PI/rosy));
			}
		}
		if (show_normal_)
		{
			glColor3f(0.3, 0.3, 0.8);
			for (int i=0; i<vertices.size(); ++i)
			{
				glVertex(vertices[i]);
				glVertex(vertices[i] + scale * vnorm[i]);
			}
		}
		glEnd();
		glDepthRange(0.0, 1.0);	
	}

	void CVTGraphics::draw_seed_project_triangle() {
		const std::vector<vec3> &vertices = cvt_->vertices();
		int nvert = vertices.size();
		const std::vector<vec3> &vw = cvt_->seed_interpolation_weight();
		const std::vector<std::vector<vec3>> &vtri = cvt_->seed_project_triangle();
		if ((vw.size()!=nvert) || (vtri.size()!=nvert)) return;
		//
		glDisable(GL_LIGHTING);
		glDepthRange(0.0, 0.998);
		glLineWidth(2);
		glBegin(GL_LINES);
		for (int c=0; c<nvert; ++c)
		{
			bool abnormal = false;
			const vec3 &w = vw[c];
			for (int i=0; i<3; ++i) if (w[i]<0 || w[i]>1) {
				abnormal = true;
			}
			const std::vector<vec3> &tri = vtri[c];
			vec3 tri_center = (tri[0]+tri[1]+tri[2]) /3.0;
			if (abnormal) {
				glColor3f(1, 0, 0); //red
			} else {
				glColor3f(0, 1, 0); //green
			}
			glVertex(tri_center);
			glVertex(vertices[c]);
		}
		glEnd();
		glDepthRange(0.0, 1.0);	
	}

	void CVTGraphics::draw_seed_match_field() {
		typedef std::pair<int, int> Edge;
		typedef std::pair<double, int> DirAngleOrder;  //<angle, order> 
		const std::vector<vec3> &vertices = cvt_->vertices();
		const std::map<Edge, DirAngleOrder> &matches = cvt_->primal_edge_field_match();
		const std::map<Edge, int>			&matches_history = cvt_->primal_edge_match_history();
		const std::vector<vec3> &field = cvt_->seed_field();
		const std::vector<vec3> &normal = cvt_->seed_normal();
		glDisable(GL_LIGHTING);
		glDepthRange(0.0, 0.998);
		glLineWidth(2);
		//glColor3f(0.0, 0.8, 0.5);
		double scale = 0.3 * sqrt(cvt_->average_area());
		glBegin(GL_LINES);
		for (auto mit=matches.begin(); mit!=matches.end(); ++mit)
		{
			const Edge &e = mit->first;
			const DirAngleOrder &d = mit->second;
			int match_number = matches_history.at(e);
			if (d.second > max_match_dir_order_) continue;
			if (match_number > max_match_number_) continue;
			vec3 hsv(d.first*180.0/M_PI, 1, 1.0/(1+d.second));  //color according to dir angle and order
			vec3 rgb;
			hsv2rgb(hsv, rgb);
			if (!colorize_) rgb=vec3(0,0,0);
			glColor(rgb);

			//edge
			vec3 v0 = vertices[e.first];
			vec3 v1 = vertices[e.second];
			vec3 mid = 0.5 * (v0+v1);
			glVertex(v0);
			glVertex(mid);
			//dir
			//vec3 dir = cvt_->vec_rotate(field[e.first], normal[e.first], d.first);
			//glVertex(v0);
			//glVertex(v0 + scale * dir);
		}
		glEnd();
		glDepthRange(0.0, 1.0);	
	}

	void CVTGraphics::draw_QD_topo_operations() {
		const std::vector<vec3> &vertices = cvt_->vertices();
		const auto &edge_dir_match = cvt_->primal_edge_field_match();

		QDmeshEditor qEd(edge_dir_match, vertices);

		std::vector<vec3> to_insert;
		std::set<int> to_remove;
		qEd.collect_topo_operations(to_insert, to_remove);

		typedef std::pair<int, int> Edge;
		const std::map<Edge, Edge> &bE = qEd.brokenEdges;
		const std::map<Edge, Edge> &lE = qEd.lackEdges;

		//draw broken/lack edges
		glDisable(GL_LIGHTING);
		glDepthRange(0.0, 0.996);
		glLineWidth(2);
		double scale = 0.3 * sqrt(cvt_->average_area());
		glBegin(GL_LINES);
		glColor3f(1, 0, 0);
		for (auto it=bE.begin(); it!=bE.end(); ++it)
		{
			Edge e = it->first;
			glVertex(vertices[e.first]);
			glVertex(0.5 * (vertices[e.first] + vertices[e.second]));
		}
		glColor3f(0, 0.6, 0.9);
		for (auto it=lE.begin(); it!=lE.end(); ++it)
		{
			Edge e = it->first;
			glVertex(vertices[e.first]);
			glVertex(vertices[e.second]);
		}
		glEnd();

		//draw insert/remove vertices
		glPointSize(3);
		glBegin(GL_POINTS);
		glColor3f(1, 0.5, 0.2);
		for (auto it=to_remove.begin(); it!=to_remove.end(); ++it) {
			glVertex(vertices[*it]);
		}
		glColor3f(0, 1, 0.5);
		for (auto it=to_insert.begin(); it!=to_insert.end(); ++it) {
			glVertex(*it);
		}
		glEnd();
		glDepthRange(0.0, 1.0);	
	}


	void CVTGraphics::draw_lloyd_grad() {
		glDisable(GL_LIGHTING);
		glDepthRange(0.0, 0.998);
		glLineWidth(2);
		glColor3f(1.0, 0.0, 0.0); //red
		if (cvt_->average_area() <= 0) cvt_->update_area();
		double area = cvt_->average_area();
		double scale = - sqrt(area) / (area * area);
		scale *= grad_magnification_;
		const std::vector<vec3>& vertices = cvt_->vertices();
		const std::vector<vec3>& grad = cvt_->lloyd_grad();
		glBegin(GL_LINES);
		for (int i=0; i<vertices.size(); ++i)
		{
			glVertex(vertices[i]);
			glVertex(vertices[i] + scale * grad[i]);
		}
		glEnd();
		glDepthRange(0.0, 1.0);	
	}

	void CVTGraphics::draw_direction_grad() {
		glDisable(GL_LIGHTING);
		glDepthRange(0.0, 0.998);
		glLineWidth(2);
		glColor3f(0.0, 1.0, 0.0); //green
		if (cvt_->average_area() <= 0) cvt_->update_area();
		double area = cvt_->average_area();
		double scale = - sqrt(area) / (area * area);
		scale *= grad_magnification_;
		const std::vector<vec3>& vertices = cvt_->vertices();
		const std::vector<vec3>& grad = cvt_->direction_grad();
		glBegin(GL_LINES);
		for (int i=0; i<vertices.size(); ++i)
		{
			glVertex(vertices[i]);
			glVertex(vertices[i] + scale * grad[i]);
		}
		glEnd();
		glDepthRange(0.0, 1.0);	
	}

	//class DrawShortEdges {
	//public:
	//	DrawShortEdges(const std::vector<vec3>& seed_in,
	//		const std::vector<unsigned>& degree_in,
	//		bool isolate_in = false)
	//		: seed_(seed_in), degree_(degree_in), isolate(isolate_in) {
	//			incident_vertices = new std::set<unsigned>();
	//			count = new unsigned[1];
	//			count[0] = 0;
	//	}

	//	~DrawShortEdges() {
	//		delete incident_vertices;
	//		std::cout << count[0] <<  " short edges." << std::endl;
	//		delete[] count;

	//	}

	//	void operator()(unsigned int i, unsigned int j, unsigned int k) const {
	//		unsigned int index[3] = {i, j, k};
	//		for (int id=0; id<3; ++id)
	//		{
	//			if (degree_[index[id]] + degree_[index[(id+1)%3]] <= 10)  //short edge condition
	//			{
	//				if (isolate)
	//				{
	//					if (incident_vertices->find(index[id]) != incident_vertices->end() ||
	//						incident_vertices->find(index[(id+1)%3]) != incident_vertices->end()) continue;
	//					incident_vertices->insert(index[id]);
	//					incident_vertices->insert(index[(id+1)%3]);
	//				}
	//				glVertex(seed_[index[id]]);
	//				glVertex(seed_[index[(id+1)%3]]);
	//				count[0]++;
	//			}
	//		}
	//	}
	//private:
	//	const std::vector<vec3>& seed_ ;
	//	const std::vector<unsigned>& degree_ ;
	//	bool isolate;
	//	std::set<unsigned> *incident_vertices;
	//	unsigned *count;
	//} ;


	//class DrawLongEdges {
	//public:
	//	DrawLongEdges(const std::vector<vec3>& seed_in,
	//		const std::vector<unsigned>& degree_in,
	//		bool isolate_in = false)
	//		: seed_(seed_in), degree_(degree_in), isolate(isolate_in) {
	//			incident_vertices = new std::set<unsigned>();
	//			count = new unsigned[1];
	//			count[0] = 0;
	//	}

	//	~DrawLongEdges() {
	//		delete incident_vertices;
	//		std::cout << count[0] <<  " long edges." << std::endl;
	//		delete[] count;
	//	}

	//	void operator()(unsigned int i, unsigned int j, unsigned int k) const {
	//		unsigned int index[3] = {i, j, k};
	//		for (int id=0; id<3; ++id)
	//		{
	//			if (degree_[index[id]] + degree_[index[(id+1)%3]] >= 14)  //long edge condition
	//			{
	//				if (isolate)
	//				{
	//					if (incident_vertices->find(index[id]) != incident_vertices->end() ||
	//						incident_vertices->find(index[(id+1)%3]) != incident_vertices->end()) continue;
	//					incident_vertices->insert(index[id]);
	//					incident_vertices->insert(index[(id+1)%3]);
	//				}
	//				glVertex(seed_[index[id]]);
	//				glVertex(seed_[index[(id+1)%3]]);
	//				count[0]++;
	//			}
	//		}
	//	}
	//private:
	//	const std::vector<vec3>& seed_ ;
	//	const std::vector<unsigned>& degree_ ;
	//	bool isolate;
	//	std::set<unsigned> *incident_vertices;
	//	unsigned *count;
	//} ;

	// 	void CVTGraphics::draw_short_edges() {
	// 		glDisable(GL_LIGHTING);
	// 		glDepthRange(0.0, 0.997);
	// 		glLineWidth(3);
	// 		glColor3f(0.5, 0.0, 0.0);  //short edge color
	// 		glBegin(GL_LINES);
	// 		cvt_->update_valences();
	// 		cvt_->RVD().for_each_primal_triangle(
	// 			DrawShortEdges(cvt_->vertices(), cvt_->valences(), true)
	// 			) ;
	// 		glEnd();
	// 		glDepthRange(0.0, 1.0);	
	// 	}

	void CVTGraphics::draw_short_edges() {
		Map* m = cvt_->delaunay_map();

		glDisable(GL_LIGHTING);
		glDepthRange(0.0, 0.997);
		glLineWidth(3);
		glColor3f(0.5, 0.0, 0.0);  //short edge color
		glBegin(GL_LINES);
		FOR_EACH_HALFEDGE(Map, m, eit) {
			if (eit->is_edge_key())
			{
				Map::Vertex* v1 = eit->prev()->vertex();
				Map::Vertex* v2 = eit->vertex();
				if (v1->degree() + v2->degree() <= 10)  //short edge condition
				{
					glVertex(v1->point());
					glVertex(v2->point());
				}
			}
		}
		glEnd();
		glDepthRange(0.0, 1.0);	
	}

	// 	void CVTGraphics::draw_long_edges() {
	// 		glDisable(GL_LIGHTING);
	// 		glDepthRange(0.0, 0.997);
	// 		glLineWidth(3);
	// 		glColor3f(0.0, 0.0, 0.5);  //long edge color
	// 		glBegin(GL_LINES);
	// 		cvt_->update_valences();
	// 		cvt_->RVD().for_each_primal_triangle(
	// 			DrawLongEdges(cvt_->vertices(), cvt_->valences(), true)
	// 			) ;
	// 		glEnd();
	// 		glDepthRange(0.0, 1.0);	
	// 	}

	void CVTGraphics::draw_long_edges() {
		Map* m = cvt_->delaunay_map();

		glDisable(GL_LIGHTING);
		glDepthRange(0.0, 0.997);
		glLineWidth(3);
		glColor3f(0.0, 0.0, 0.5);  //long edge color
		glBegin(GL_LINES);
		FOR_EACH_HALFEDGE(Map, m, eit) {
			if (eit->is_edge_key())
			{
				Map::Vertex* v1 = eit->prev()->vertex();
				Map::Vertex* v2 = eit->vertex();
				if (v1->degree() + v2->degree() >= 14)  //long edge condition
				{
					glVertex(v1->point());
					glVertex(v2->point());
				}
			}
		}
		glEnd();
		glDepthRange(0.0, 1.0);	
	}

	typedef std::pair<Map::Halfedge_iterator, double> EdgeValue;
	// 	struct EdgeValueGreater
	// 	{
	// 		bool operator() (const EdgeValue& a, const EdgeValue& b) { return a.second > b.second; }
	// 	};
	// 	struct EdgeValueLess
	// 	{
	// 		bool operator() (const EdgeValue& a, const EdgeValue& b) { return a.second < b.second; }
	// 	};

	void CVTGraphics::draw_stretched_long_edges() {
		cvt_->update_map_edge_stretch();
		Map* m = cvt_->delaunay_map();
		const auto& edge_stretch = cvt_->delaunay_map_edge_stretch();

		std::vector<EdgeValue> long_edges;

		FOR_EACH_HALFEDGE(Map, m, eit) {
			if (eit->is_edge_key())
			{
				Map::Vertex* v1 = eit->prev()->vertex();
				Map::Vertex* v2 = eit->vertex();
				if (v1->degree() != 6 || v2->degree() != 6)
				{
					double stretch = edge_stretch.at(eit);
					if (stretch > 0)
					{
						long_edges.push_back(std::make_pair(eit, stretch));
					}
				}
			}
		}

		std::sort(long_edges.begin(), long_edges.end(), DelaunayCVT::EdgeValueGreater());
		double num_edges = long_edges.size() * show_stretched_edges_rate_;
		glDepthRange(0.0, 0.997);
		glLineWidth(3);
		glColor3f(0.0, 0.0, 0.5);
		glBegin(GL_LINES);
		for (int i=0; i<num_edges; ++i)
		{
			auto& eit = long_edges[i].first;
			glVertex(eit->prev()->vertex()->point());
			glVertex(eit->vertex()->point());
		}
		glEnd();
		glDepthRange(0.0, 1.0);	
	}

	void CVTGraphics::draw_stretched_short_edges() {
		cvt_->update_map_edge_stretch();
		Map* m = cvt_->delaunay_map();
		const auto& edge_stretch = cvt_->delaunay_map_edge_stretch();

		std::vector<EdgeValue> short_edges;

		FOR_EACH_HALFEDGE(Map, m, eit) {
			if (eit->is_edge_key())
			{
				Map::Vertex* v1 = eit->prev()->vertex();
				Map::Vertex* v2 = eit->vertex();
				if (v1->degree() != 6 || v2->degree() != 6)
				{
					double stretch = edge_stretch.at(eit);
					if (stretch < 0)
					{
						short_edges.push_back(std::make_pair(eit, stretch));
					}
				}
			}
		}

		std::sort(short_edges.begin(), short_edges.end(), DelaunayCVT::EdgeValueLess());
		double num_edges = short_edges.size() * show_stretched_edges_rate_;
		glDepthRange(0.0, 0.997);
		glLineWidth(3);
		glColor3f(0.5, 0.0, 0.0);
		glBegin(GL_LINES);
		for (int i=0; i<num_edges; ++i)
		{
			auto& eit = short_edges[i].first;
			glVertex(eit->prev()->vertex()->point());
			glVertex(eit->vertex()->point());
		}
		glEnd();
		glDepthRange(0.0, 1.0);	
	}


	/* draw edge flip, edge collapse and vertex split operations that will decrease regularity energy
	* operator_code : which operators to draw
	* isolation: whether connected operations are allowed
	*/
	void CVTGraphics::draw_valid_topo_operation(int operator_code, int deltaR_threshold /*=0*/, bool isolation /*= true*/)
	{		
		cvt_->update_delaunay_map();
		Map* m = cvt_->delaunay_map();

		// collect valid topo operations
		std::set<Map::Vertex*> visited; //visited vertices
		std::vector<TopoOperation> S;
		double dihedral = cos(cvt_->dihedral_angle_threshold()*M_PI/180);
		//
		//const std::map<Map::Vertex*, int> &map_vert_idx = cvt_->delaunay_map_vert_index();
		FOR_EACH_VERTEX(Map, m, vit) {
			int dv = vit->degree();  //TODO: use virtual degree
			if (dv != 6) {
				TopoOperation o;
				//test
				bool found = (cvt_->lock_feature()) ? cvt_->find_topo_operation_preserve_fea(vit, o) : cvt_->find_topo_operation(vit, o, dihedral);
				if (!found) continue;
				//if(!cvt_->find_topo_operation(vit, o, dihedral)) continue;
				if (o.delta_R() >= deltaR_threshold) continue; //threshold filter
				TopoOperator t = o.type();
				if(t != Geex::Null_Operator) {
					Map::Halfedge* e = o.edge();
					Map::Vertex* v1 = e->vertex();
					Map::Vertex* v2 = e->prev()->vertex();
					if (visited.find(v1) == visited.end() && visited.find(v2) == visited.end()) {
						if (t == Vertex_Split) {
							e = o.counter_edge();
							Map::Vertex* v3 = e->prev()->vertex();
							if (visited.find(v3) == visited.end()) {
								visited.insert(v1);
								visited.insert(v2);
								visited.insert(v3);
								S.push_back(o);
							}
						}
						else { //t = EF, EC, ES
							visited.insert(v1);
							visited.insert(v2);
							S.push_back(o);
						}
					}
				}
			}
		}

		// draw topo operations
		bool show_edge_flip = operator_code & Geex::Edge_Flip;
		bool show_edge_collapse = operator_code & Geex::Edge_Collapse;
		bool show_vertex_split = operator_code & Geex::Vertex_Split;
		glDisable(GL_LIGHTING);
		glDepthRange(0.0, 0.998);
		glLineWidth(2);

		if (show_edge_flip)
		{
			glColor3f(1.0, 0.0, 0.0); //red
			glBegin(GL_LINES);
			for (auto oit = S.begin(); oit != S.end(); ++oit)
			{
				if (oit->type() == Edge_Flip) {
					glVertex(oit->edge()->prev()->vertex()->point());
					glVertex(oit->edge()->vertex()->point());
				}
			}
			glEnd();
		}
		if (show_edge_collapse)
		{
			glColor3f(0.0, 1.0, 0.0); //green
			glBegin(GL_LINES);
			for (auto oit = S.begin(); oit != S.end(); ++oit)
			{
				if (oit->type() == Edge_Collapse) {
					glVertex(oit->edge()->prev()->vertex()->point());
					glVertex(oit->edge()->vertex()->point());
				}
			}
			glEnd();
		}
		if (show_vertex_split)
		{
			glColor3f(0.0, 0.0, 1.0); //blue
			glBegin(GL_LINES);
			for (auto oit = S.begin(); oit != S.end(); ++oit)
			{
				if (oit->type() == Vertex_Split) {

					Map::Halfedge* e1 = oit->edge();
					Map::Halfedge* e2 = oit->counter_edge();
					glVertex(e1->prev()->vertex()->point());
					glVertex(e1->vertex()->point());
					glVertex(e2->prev()->vertex()->point());
					glVertex(e2->vertex()->point());
				}
			}
			glEnd();
		}
		glDepthRange(0.0, 1.0);	
	}


	///

	class DrawPrimal {
	public:
		DrawPrimal( const std::vector<vec3>& seed_in,
			double min_angle_in, 
			double max_angle_in,
			GLboolean show_bad) 
			: seed_(seed_in), 
			min_angle_(min_angle_in), 
			max_angle_(max_angle_in),
			show_obtuse_(show_bad){ }
		void operator()(unsigned int i, unsigned int j, unsigned int k) const {
			//const vec3& v1 = seed_[i] ;
			//const vec3& v2 = seed_[j] ;
			//const vec3& v3 = seed_[k] ;
			//glNormal(cross(v2-v1, v3-v1)) ;
			//glVertex(v1) ;
			//glVertex(v2) ;
			//glVertex(v3) ;
			const vec3& v1=seed_[i] ;
			const vec3& v2=seed_[j] ;
			const vec3& v3=seed_[k] ;
			vec3 ab = v2 - v1 ; //F.vertex[1] - F.vertex[0] ;
			vec3 bc = v3 - v2 ; //F.vertex[2] - F.vertex[1] ;
			vec3 ca = v1 - v3 ; //F.vertex[0] - F.vertex[2] ;
			real lab = length(ab) ;
			real lbc = length(bc) ;
			real lca = length(ca) ;
			real anglea = acos(-dot(ca, ab)/(lca*lab))*180. / M_PI ;
			real angleb = acos(-dot(ab, bc)/(lab*lbc))*180. / M_PI ;
			real anglec = acos(-dot(bc, ca)/(lbc*lca))*180. / M_PI ;
			real amin=180, amax=0 ;
			amin = gx_min(anglea, gx_min(angleb, anglec)) ;
			amax = gx_max(anglea, gx_max(angleb, anglec)) ;

			if(show_obtuse_) {
				if(amax > max_angle_ && amin < min_angle_) {
					glColor3f(0.5f, 0.5f, 0.5f) ;
				}
				else if(amin < min_angle_) {
					glColor3f(0.609f, 0.838f, 0.962f) ;
				} 
				else if (amax > max_angle_) {
					glColor3f(0.949f, 0.577f, 0.629f) ;
					//glColor3f(0.609f, 0.962, 0.609f) ;
				}
				else 
					glColor3f(1, 1, 1) ;
			}else
				glColor4fv(white) ;

			glNormal(cross(v2-v1, v3-v1)) ;
			glVertex(v1) ;
			glVertex(v2) ;
			glVertex(v3) ;
		}
	private:
		const std::vector<vec3>& seed_ ;
		double min_angle_, max_angle_ ;
		GLboolean show_obtuse_ ;
	} ;

	class DrawPrimalMesh {
	public:
		DrawPrimalMesh(const std::vector<vec3>& seed_in)
			: seed_(seed_in) {}

		void operator()(unsigned int i, unsigned int j, unsigned int k) const {
			const vec3& v1 = seed_[i] ;
			const vec3& v2 = seed_[j] ;
			const vec3& v3 = seed_[k] ;
			glNormal(cross(v2-v1, v3-v1)) ;
			glVertex(v1) ;
			glVertex(v2) ;
			glVertex(v3) ;
		}
	private:
		const std::vector<vec3>& seed_ ;
	} ;

	///dxy add
	void CVTGraphics::draw_primal_map() {
		cvt_->update_delaunay_map();
		Map *m = cvt_->delaunay_map();

		glEnable(GL_LIGHTING) ;
		glColor3f(1, 1, 1) ;
		glColor4fv(white) ;
		glBegin(GL_TRIANGLES) ;
		FOR_EACH_FACET(Map, m, fit) {
			const vec3& v1= fit->halfedge()->vertex()->point() ;
			const vec3& v2= fit->halfedge()->next()->vertex()->point() ;
			const vec3& v3= fit->halfedge()->prev()->vertex()->point() ;
			double amin, amax;
			min_max_angle(v1, v2, v3, amin, amax);

			if(show_obtuse_) {
				if(amax > cvt_->angle_max() && amin < cvt_->angle_max()) {
					glColor3f(0.949f, 0.577f, 0.629f) ;
				}
				else if(amin < cvt_->angle_min()) {
					glColor3f(0.609f, 0.838f, 0.962f) ;
				} 
				else if (amax > cvt_->angle_max()) {
					glColor3f(0.609f, 0.962, 0.609f) ;
				}
				else 
					glColor3f(1, 1, 1) ;
			}else
				glColor4fv(white) ;

			glNormal(fit->normal()) ;
			glVertex(v1) ;
			glVertex(v2) ;
			glVertex(v3) ;
		}
		glEnd() ;

		// draw wireframe
		glDepthRange(0.0, 0.999) ;
		glDisable(GL_LIGHTING) ;
		glLineWidth(1.0) ;
		glColor3f(0.0, 0.0, 0.0) ;
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) ;
		glBegin(GL_LINES);
		FOR_EACH_HALFEDGE(Map, m, eit) {
			/*if (eit->is_edge_key())*/
			{
				glVertex(eit->prev()->vertex()->point());
				glVertex(eit->vertex()->point());
			}
		}
		glEnd() ;
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;
		glDepthRange(0.0, 1.0) ;
	}



	void CVTGraphics::draw_sharp_edges() {
		cvt_->update_delaunay_map();
		Map *m = cvt_->delaunay_map();

		double dihedral_dot = cos(cvt_->dihedral_angle_threshold()*M_PI/180);
		//
		glDepthRange(0.0, 0.999) ;
		//glDisable(GL_LIGHTING) ;
		glLineWidth(3.5) ;
		glColor3f(1, 0, 0) ;  //purple
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) ;
		glBegin(GL_LINES);
		FOR_EACH_HALFEDGE(Map, m, eit) {			
			if (show_sharp_edges_ && eit->is_sharp(dihedral_dot))
			{
				glVertex(eit->prev()->vertex()->point());
				glVertex(eit->vertex()->point());
			}
		}
		glEnd() ;
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;
		glDepthRange(0.0, 1.0) ;
	}



	void CVTGraphics::draw_primal_features() {
		const auto &vidx = cvt_->delaunay_map_vert_index();
		const auto &locked = cvt_->locked();
		const auto &on_feature = cvt_->on_feature();
		//cvt_->mark_feature_on_delaunay_map();
		Map *m = cvt_->delaunay_map();
		// feature points
		int fp_cnt = 0;
		glDisable(GL_LIGHTING) ;
		glDepthRange(0.0, 0.998) ;
		glHint(GL_POINT_SMOOTH_HINT, GL_NICEST) ;
		glEnable(GL_POINT_SMOOTH) ;
		glPointSize(vertices_size_ * 20.0f);
		glColor3f(1, 0, 0);
		glBegin(GL_POINTS);
		FOR_EACH_VERTEX(Map, m, it) {
			//if (cvt_->delaunay_is_on_feature(it)) 
			if (vidx.find(it) != vidx.end())
			{
				int id = vidx.at(it);
				if (locked[id] || on_feature[id]) {
					glVertex(it->point());
					++fp_cnt;
				}
			}
		}
		glEnd();
		//feature lines
		// 		int fe_cnt = 0;
		// 		glDisable(GL_LIGHTING) ;
		// 		glLineWidth(1.0) ;
		// 		glColor3f(0.0, 1.0, 0.0) ;
		// 		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) ;
		// 		glBegin(GL_LINES);
		// 		FOR_EACH_EDGE(Map, m, eit) {
		// 			Map::Vertex *v0 = eit->prev()->vertex();
		// 			Map::Vertex *v1 = eit->vertex();
		// 			if (cvt_->delaunay_is_on_feature(v0) && cvt_->delaunay_is_on_feature(v1) 
		// 				&& (cvt_->delaunay_feature_line(v0)==cvt_->delaunay_feature_line(v1) 
		// 				|| cvt_->delaunay_is_corner(v0) || cvt_->delaunay_is_corner(v1)))
		// 			{
		// 				glVertex(v0->point());
		// 				glVertex(v1->point());
		// 				++fe_cnt;
		// 			}
		// 
		// 		}
		// 		glEnd() ;
		// 		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;
		glDepthRange(0.0, 1.0) ;
		//test
		std::cout << fp_cnt << " feature pts." << std::endl;
		//std::cout << fe_cnt << " feature edges." << std::endl;

	}

	void CVTGraphics::draw_primal_map_vertices(bool show_singularity /* = false */) {
		cvt_->update_delaunay_map();
		Map *m = cvt_->delaunay_map();

		glDisable(GL_LIGHTING) ;
		glDepthRange(0.0, 0.998) ;
		glHint(GL_POINT_SMOOTH_HINT, GL_NICEST) ;
		glEnable(GL_POINT_SMOOTH) ;
		glPointSize(vertices_size_ * 20.0f);
		glColor3f(0, 0, 0);
		glBegin(GL_POINTS);
		if (show_singularity)
		{
			int num_singular[13] = {0};
			FOR_EACH_VERTEX(Map, m, vit) {
				int d = vit->degree();
				if (d > 11) num_singular[12]++;
				else num_singular[d]++;

				if (d < 6) {
					glColor3f(0.2, 0.33, 0.65);
					glVertex(vit->point());
				}
				else if (d > 6) {
					glColor3f(0.89, 0.27, 0.51);
					glVertex(vit->point());
				}
				else {
					//glColor3f(0, 0, 0);
				}
				
			}
			//print singularity num info
			std::cout << "Singularity: [";
			for(int i=1; i<13; ++i) std::cout << num_singular[i] << ", ";
			std::cout << "]" << std::endl;
		}
		else
		{
			FOR_EACH_VERTEX(Map, m, vit) {
				glVertex(vit->point());
			}
		}
		glEnd();
		glDepthRange(0.0, 1.0);
	}
	///

	void CVTGraphics::draw_primal() {
		glEnable(GL_LIGHTING) ;
		glColor3f(1, 1, 1) ;
		glColor4fv(white) ;

		glBegin(GL_TRIANGLES) ;
		cvt_->RVD().for_each_primal_triangle(
			DrawPrimal(cvt_->vertices(),  cvt_->angle_min(), cvt_->angle_max(), show_obtuse_)
			) ;
		glEnd() ;

		// draw wireframe
		glDepthRange(0.0, 0.999) ;
		glDisable(GL_LIGHTING) ;
		glLineWidth(1.0) ;
		glColor3f(0.0, 0.0, 0.0) ;

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) ;
		glBegin(GL_TRIANGLES) ;
		cvt_->RVD().for_each_primal_triangle(
			DrawPrimalMesh(cvt_->vertices())
			) ;
		glEnd() ;
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;
		glDepthRange(0.0, 1.0) ;
	}

	void CVTGraphics::draw_obtuse_triangles() {
		/*		cvt_->update_primal() ;
		std::vector<vec3>& vertices = cvt_->vertices() ;
		std::vector<PrimalTri>& tris = cvt_->primal_facets() ;

		glDepthRange(0.0, 0.999) ;

		glEnable(GL_LIGHTING) ;
		glColor3f(1, 1, 1) ;
		glBegin(GL_TRIANGLES) ;
		for(unsigned int i=0; i<tris.size(); ++i) {
		PrimalTri& t = tris[i] ;
		vec3 v0 = vertices[t[0]] ;
		vec3 v1 = vertices[t[1]] ;
		vec3 v2 = vertices[t[2]] ;
		double amin, amax ;
		min_max_angle(v0, v1, v2, amin, amax) ;

		if(amax > cvt_->angle_max() && amin < cvt_->angle_max()) {
		glColor3f(0.949f, 0.577f, 0.629f) ;
		}
		else if(amin < cvt_->angle_min()) {
		glColor3f(0.609f, 0.838f, 0.962f) ;
		} 
		else if (amax > cvt_->angle_max()) {
		glColor3f(0.609f, 0.962, 0.609f) ;
		}
		else 
		glColor3f(1, 1, 1) ;
		glNormal(cross(vertices[t[1]]-vertices[t[0]], vertices[t[2]]-vertices[t[0]])) ;
		for(int j=0; j<3; ++j) {
		glVertex(vertices[t[j]]) ;
		}
		}
		glEnd() ;

		glDisable(GL_LIGHTING) ;	

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) ;
		glLineWidth(1.0) ;
		glColor3f(0.1, 0.1, 0.1) ;
		glBegin(GL_TRIANGLES) ;
		for(unsigned int i=0; i<tris.size(); ++i) {
		PrimalTri& t = tris[i] ;
		glNormal(cross(vertices[t[1]]-vertices[t[0]], vertices[t[2]]-vertices[t[0]])) ;
		for(int j=0; j<3; ++j) {
		glVertex(vertices[t[j]]) ;
		}
		}
		glEnd() ;
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;

		glDepthRange(0.0, 1.0) ;*/	
	}


	inline void c_draw_triangle(
		unsigned int c, const vec3& v1, const vec3& v2, const vec3& v3
		) {
			glVertex(v1) ;
			glVertex(v2) ;
			glVertex(v3) ;
	} 

	inline void c_draw_triangle_with_normal(
		unsigned int c, const vec3& v1, const vec3& v2, const vec3& v3
		) {
			glNormal(cross(v2-v1, v3-v1)) ;
			glVertex(v1) ;
			glVertex(v2) ;
			glVertex(v3) ;
	} 

	inline void c_draw_triangle_with_normal_and_color(
		unsigned int c, const vec3& v1, const vec3& v2, const vec3& v3
		) {
			gl_random_color(c) ;
			glNormal(cross(v2 - v1, v3 - v1)) ;
			glVertex(v1) ;
			glVertex(v2) ;
			glVertex(v3) ;
	}

	class DrawTriangleValence {
	public:
		DrawTriangleValence(std::vector<unsigned>&  degrees)
			: degrees_(degrees) {
		}

		void operator () (unsigned int c, const vec3& v1, const vec3& v2, const vec3& v3) const {
			gl_vertex_color(degrees_[c]) ;
			glNormal(cross(v2 - v1, v3 - v1)) ;
			glVertex(v1) ;
			glVertex(v2) ;
			glVertex(v3) ;
		}

	private: 
		std::vector<unsigned>&  degrees_ ;
	} ;

	void CVTGraphics::draw_rvd() {
		glColor4fv(white) ;
		glEnable(GL_LIGHTING) ;
		glBegin(GL_TRIANGLES) ;
		if(colorize_) {
			cvt_->RVD().for_each_triangle(
				c_draw_triangle_with_normal_and_color
				) ;			
		}
		else {
			///dxy add
			cvt_->update_valences();
			///
			cvt_->RVD().for_each_triangle(
				DrawTriangleValence(cvt_->valences())
				) ;			
		}
		glEnd() ;
	}

	//dxy add
	class DrawCellValue {
	public:
		DrawCellValue(const std::vector<double> &values_in)
			: values(values_in) {}

		virtual ~DrawCellValue() {};

		virtual vec3 colormap(double val) const {
			return vec3(1.0, 1.0, 1.0);
		};

		void operator () (unsigned int c, const vec3& v1, const vec3& v2, const vec3& v3) const {
			glColor(colormap(values[c]));
			glNormal(cross(v2 - v1, v3 - v1)) ;
			glVertex(v1) ;
			glVertex(v2) ;
			glVertex(v3) ;
		}

	protected:
		const std::vector<double> &values;
	};


	/** 
	* Given a truncate value, DrawCellFixedValue maps a value to a unique color, regardless of the min/max value of input values
	*/
	class DrawCellFixedValue : public DrawCellValue {
	public:
		DrawCellFixedValue(const std::vector<double> &values_in, int h_lo_in, int h_hi_in, double truncate_val_in)
			: DrawCellValue(values_in), h_lo(h_lo_in), h_hi(h_hi_in), truncate_val(truncate_val_in) {}

		vec3 colormap(double val) const {
			double s;
			if (truncate_val <= 0) //map R+ to [0, 1]
			{
				s = 1 - exp(-fabs(val));
			}
			else { //map [0, truncate_val] to [0, 1]
				if (fabs(val) > truncate_val) s=1;
				else {
					s = fabs(val)/truncate_val;
				}
			}

			vec3 hsv((val<0 ? h_lo : h_hi), s, 1.0);
			vec3 rgb;
			hsv2rgb(hsv, rgb);
			return rgb;
		}

	protected:
		int h_lo, h_hi;  //hsv hue for value lower/upper bound
		double truncate_val; //if positive, 
	};

	/**
	* DrawCellRelativeValue maps input values so that:
	* min/max of input values is mapped to color_min/color_max (hsv color), 
	* values in between are mapped to hsv colors between color_min and color_max
	*/
	class DrawCellRelativeValue : public DrawCellValue {
	public:
		DrawCellRelativeValue(const std::vector<double> &values_in, vec3 color_min_in, vec3 color_max_in)
			: DrawCellValue(values_in), cmin(color_min_in), cmax(color_max_in)
		{
			auto mm = std::minmax_element(values.begin(), values.end());
			vmin = *(mm.first);
			vmax = *(mm.second);		
		}

		vec3 colormap(double val) const {
			double x = (val-vmin)/(vmax-vmin);

			vec3 hsv = x*cmax + (1-x)*cmin;
			vec3 rgb;
			hsv2rgb(hsv, rgb);
			return rgb;
		}

	protected: 
		double vmin, vmax;
		vec3 cmin, cmax; //hsv color for vmin and vmax
	} ;



	void CVTGraphics::draw_cell_energy(GLenum energy_mode) {
		//hsv color
		vec3 cmin(120, 1, 1); //green
		vec3 cmax(0, 1, 1);   //red

		glEnable(GL_LIGHTING) ;
		glBegin(GL_TRIANGLES) ;
		if (energy_mode == CVT_True_Minus_Approx)
		{
			std::vector<double> true_minus_approx;
			const auto &true_cvt = cvt_->true_cvt_energy();
			const auto &approx_cvt = cvt_->approx_cvt_energy();
			for (int i=0; i<true_cvt.size(); ++i) true_minus_approx.push_back(true_cvt[i] - approx_cvt[i]);
			cvt_->RVD().for_each_triangle(DrawCellFixedValue(true_minus_approx, 240/*blue*/, 0/*red*/, cell_energy_truncate_value_));
		}
		else {
			switch (energy_mode)
			{
			case CVT_True:
				cvt_->RVD().for_each_triangle(DrawCellRelativeValue(cvt_->true_cvt_energy(), cmin, cmax));
				break;
			case CVT_Approx:
				cvt_->RVD().for_each_triangle(DrawCellRelativeValue(cvt_->approx_cvt_energy(), cmin, cmax));
				break;
			case DIR:
				cvt_->RVD().for_each_triangle(DrawCellRelativeValue(cvt_->direction_energy(), cmin, cmax));
				break;
			case Cell_Area:
				cvt_->RVD().for_each_triangle(DrawCellRelativeValue(cvt_->rvd_cell_area(), cmin, cmax));
				break;
			default:
				break;
			}
		}
		glEnd();
	}
	//

	void c_draw_edge(unsigned int c, const TopoPolyVertexEdge& v1, const TopoPolyVertexEdge& v2) {
		glVertex(v1) ;
		glVertex(v2) ;
	}

	void CVTGraphics::draw_rvd_mesh() {
		glDisable(GL_LIGHTING) ;
		glLineWidth(1) ;
		glDepthRange(0.0, 0.998) ;
		glBegin(GL_LINES) ;
		///dxy change
		//glColor3f(0.0, 0.0, 0.0) ;
		glColor3f(0.0, 0.0, 0.8);
		cvt_->RVD().for_each_halfedge(c_draw_edge) ;
		glEnd() ;
		glDepthRange(0.0, 1.0) ;
	}

	///dxy add: to be called by RVD.for_each_facet
	class DrawSeedFeatureLine {
	public:
		DrawSeedFeatureLine (
			const std::map<std::pair<int, int>, int> &featureEdgeLine_in,
			const RestrictedVoronoiDiagram_poly& RVD_in
			) : featureEdgeLine(featureEdgeLine_in), RVD_(RVD_in), B(RVD_in.mesh()) {}

		void operator() (unsigned int c, TopoPolyMesh* M) const {
			int fi = RVD_.current_facet();
			int vidx[3];
			for (int i=0; i<3; ++i)
			{
				vidx[i] = B->facet_begin(fi) + i;
				vidx[i] = B->vertex_index(vidx[i]);
			}
			// detect feature edges of current facet
			std::vector<int> feature_edges;
			std::vector<int> fLine;
			for (int i=0; i<3; ++i) {
				int v0 = vidx[i];
				int v1 = vidx[(i+1)%3];
				int a = (v0<v1) ? v0 : v1;
				int b = (v0<v1) ? v1 : v0;
				auto it = featureEdgeLine.find(std::make_pair(a, b));
				if (it != featureEdgeLine.end()) {
					feature_edges.push_back(i);
					fLine.push_back(it->second);
				}
			}
			// find if feature edges intersect with rvd cell of seed c
			if (!feature_edges.empty()) {
				//compute CVT and dual length
				for(unsigned int f=0; f<M->nb_facets(); ++f) {
					for(unsigned int i = M->facet_begin(f); i<M->facet_end(f); ++i) {
						vec3 v0 = M->vertex(i);
						vec3 v1;
						if (i+1 < M->facet_end(f)) v1 = M->vertex(i+1);
						else v1 = M->vertex(M->facet_begin(f));
						for (int j=0; j<feature_edges.size(); ++j) {
							int i0 = feature_edges[j];
							vec3 s = B->original_vertices()[vidx[i0]];
							vec3 t = B->original_vertices()[vidx[(i0+1)%3]];
							if (partof(v0, v1, s, t)) {								
								glVertex(v0);
								glVertex(v1);
							}
						}
					}
				}
			}				
		}

		bool colinear(const vec3 &e1, const vec3 &e2) const {
			if ((Geex::length2(e1) < 1e-10) || (Geex::length2(e2) < 1e-10)) return true;
			else {
				vec3 n1 = Geex::normalize(e1);
				vec3 n2 = Geex::normalize(e2);
				return fabs(1 - fabs(Geex::dot(n1, n2))) < 1e-10;
			}
		} 

		bool partof(const vec3 &p1, const vec3 & p2, const vec3& q1, const vec3 &q2) const {
			vec3 e1 = p2 - p1;
			vec3 e2 = q2 - q1;
			vec3 e3 = p1 - q1;
			if (colinear(e1, e2) && colinear(e2, e3)) {
				double a = Geex::dot(p1-q1, q2-q1) / Geex::length2(q2-q1);   
				double b = Geex::dot(p2-q1, q2-q1) / Geex::length2(q2-q1);
				if ((fabs(a-0.5)<0.5 + 1e-10) && (fabs(b-0.5)<0.5 + 1e-10)) {
					return true;
				}
			}
			return false;
		}

	private:
		const RestrictedVoronoiDiagram_poly& RVD_;
		const TopoPolyMesh *B;
		const std::map<std::pair<int, int>, int> &featureEdgeLine;
	};

	///

	void CVTGraphics::draw_feature() {
		std::vector<int>& corners = cvt_->corners() ;
		std::vector<std::vector<int> >& features = cvt_->features() ;
		std::vector<Map::Vertex*>& verts = cvt_->verts() ;
		glDisable(GL_LIGHTING) ;
		glDepthRange(0.0, 0.998) ;
		glHint(GL_POINT_SMOOTH_HINT, GL_NICEST) ;
		glEnable(GL_POINT_SMOOTH) ;
		glPointSize(vertices_size_ * 20.0f) ;
		glColor3f(1.0, 0.0, 0.0) ;
		glBegin(GL_POINTS) ;
		for(unsigned int i=0; i<corners.size(); ++i) {
			glVertex(verts[corners[i]]->point()) ;
		}
		glEnd() ;
		//dxy add: seed feature verts
		const std::map<int, std::set<int>> &sfl = cvt_->seed_feature_lines();
		const auto &verices = cvt_->vertices();
		if (sfl.empty()) {
			cvt_->compute_seed_feature_lines();
		}
		glBegin(GL_POINTS) ;
		for (auto it=sfl.begin(); it!=sfl.end(); ++it) {
			int vi = it->first;
			int nline = it->second.size();
			switch (nline)
			{
			case 1:
				glColor3f(1, 0, 0);
				break;
			case 2:
				glColor3f(0, 1, 0);
				break;
			case 3:
				glColor3f(0, 0, 1);
				break;
			default:
				glColor3f(0.5, 0.5, 0.5);
				break;
			}
			glVertex(verices[vi]);
		}
		glEnd() ;
		///


		glLineWidth(3.0f) ;
		glColor3f(0.0, 0.0, 1.0) ;
		glBegin(GL_LINES) ; 
		for(unsigned int i=0; i<features.size(); ++i) {
			std::vector<int>& line = features[i] ;
			for(unsigned int j=0; j+1<line.size(); j++) {
				glVertex(verts[line[j]]->point()) ;
				glVertex(verts[line[j+1]]->point()) ;
			}
		}
		glEnd() ;

		///dxy add: draw feature edges
		// 		std::map<std::pair<int, int>, int> featureEdgeLine;
		// 		for (int i=0; i<features.size(); ++i) {
		// 			const auto &line = features[i];
		// 			for (int j=0; j<line.size()-1; ++j) {
		// 				int v0 = line[j];
		// 				int v1 = line[j+1];
		// 				int a = (v0 < v1) ? v0 : v1;
		// 				int b = (v0 < v1) ? v1 : v0; 
		// 				featureEdgeLine.insert(std::make_pair(std::make_pair(a, b), i));
		// 			}
		// 		}
		// 		glBegin(GL_LINES);
		// 		cvt_->RVD().for_each_facet(DrawSeedFeatureLine(featureEdgeLine, cvt_->RVD()));		
		// 		glEnd() ;
		///

		glLineWidth(1.0f) ;
		glDepthRange(0.0, 1.0) ;
	}

	//dxy add
	void CVTGraphics::save_snapshot(const std::string& filename) {
		IO_.snapshot(filename);
	}
	//

}
