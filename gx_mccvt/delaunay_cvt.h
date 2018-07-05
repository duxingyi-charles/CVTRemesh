#ifndef __DELAUNAY_CVT__
#define __DELAUNAY_CVT__

//#include "local_rvd.h"
#include "cvt_primal.h"
#include <Geex/graphics/opengl.h> // use GLboolean
#include <Geex/numerics/optimizer.h>
#include <Geex/numerics/lbfgs_interface.h>
#include <Geex/combinatorics/map.h>
#include <Geex/combinatorics/map_properties.h>
#include <Geex/CVT/MedialAxis.h>
#include <Geex/CVT/ann_kdtree.h>
#include <Geex/CVT/delaunay_CGAL.h>
#include <Geex/CVT/topo_poly_mesh.h>
#include <Geex/CVT/RVD.h>
#include <Geex/graphics/opengl.h>
#include <Geex/numerics/sparse_matrix.h>
#include <fstream>

///dxy add
#include <Geex/basics/line_stream.h>
#include <Geex/combinatorics/map_editor_ext.h>
#include <functional>
#include "cvt_multithread.h"
///

namespace Geex {

	enum OptimizerMode { LBFGSB, HLBFGS, HM1QN3, HCG, HLBFGS_HESS } ;
	///dxy add
	enum DirectionEnergyMode { Angle_Only, Angle_Edge, Angle_Edge2 } ;
	enum DirectionEdgeWeightMode { Dual_Length, Lloyd_Energy };
	enum InterpolationWeightMode { Barycentric, Cotangent };
	enum TopoOperator {
		Null_Operator = 0,
		Edge_Flip = 1,
		Edge_Collapse = 2,
		Vertex_Split = 4,
		//Edge_Split = 8
	};
	class Experiment;
	class TopoOperation;
	class MatchMaker;
	class CVTMultiThread;
	///
#define OptimizerModeNames  "LBFGSB" | "HLBFGS" | "HM1QN3" | "HCG" | "HESSIAN"

	class AABBTree ;
	//	class PrimalTri ;

	class DelaunayCVT {
	public:
		DelaunayCVT(const std::string& del_algo="CGAL_spatial_sort") ;
		~DelaunayCVT() ;

		void load_map(const std::string& filename) ;      // load Map
		void load_boundary(const std::string& filename) ; // load TopoPolyMesh
		///dxy add
		void load_map(const std::string& filename, bool normalize);
		void model_normalize_map();
		void load_boundary(const std::string& filename, bool normalize) ;
		///
		bool load_feature(const std::string& filename) ;
		void mark_feature() ;
		void save_samples(const std::string& filename) ;
		///dxy add
		//void mark_feature_on_delaunay_map();
		void save_feature(const std::string &filename);
		//
		void save_denormalized_samples(const std::string& filename);
		//Note: reload denormalized samples using load_samples(f, true)
		void load_samples(const std::string& filename, bool normalize);
		void normalize_samples();
		//
		bool load_experiment(const std::string &filename);
		///
		void load_samples(const std::string& filename) ;
		void save_density(const std::string& filename) ;
		void load_density(const std::string& filename) ;
		void save_chart(const std::string& filename) ;
		void load_chart(const std::string& filename) ;
		///dxy add
		GLboolean field_loaded() { return field_loaded_; }
		void load_field(const std::string& filename) ;
		void calcFacetNormalAndDirection();
		///
		//		void load_geodist(const std::string& filename) ;
		void get_bbox(
			real& x_min, real& y_min, real& z_min,
			real& x_max, real& y_max, real& z_max
			) ;
		void save_normalized_boundary(const std::string& filename) ;
		//dxy add
		void save_denormalized_boundary(const std::string& filename);
		void save_denormalized_primal(const std::string& filename);
		//
		//	void save_primal(const std::string& filename) ;
		void save_primal_global(const std::string& filename) ;
		//	void update_primal() ;
		void update_valences(bool message = false) ;
		void count_primal_singular_triangles(std::map<trindex, int>& triangles) ;
		int  split_primal_singular_triangles(std::map<trindex, int>& triangles) ;
		int  split_non_manifold_primal_triangles() ;
		GLboolean& primal_dirty() { return primal_dirty_ ; }
		void extract_features(Map* map, double angle=30) ;
		void extract_chart() ;
		unsigned int& nb_lock_rings() { return nb_lock_rings_ ; }
		void set_locked() ;
		void release_locked() ;

		void init_samples(int nb_pts) ;
		void init_samples_w(int nb_pts) ;
		//dxy add
		void init_samples_barycentric(int& nb_pts);
		//
		void compute_rvd() ;
		bool has_feature(bool fidx) ;

		///dxy add: short/long edges edit
		void split_long_edges();
		void collapse_short_edges();

		void split_all_edges();

		void update_area();
		// optimization info
		void save_optimization_info(const std::string& filename);
		///


		void update_indices() ;
		//void cluster_vertices() ;
		//void init_facet_cells() ;
		//void collect_facet_cells() ;
		//inline bool add_facet_cell(int f, int c) ;
		//void update_facet_cells() ;
		//void compute_facet_rvd(int f, int c, TopoPolyMesh& M) ;
		//std::vector<std::vector<int> >& facet_cells() { return facet_cells_ ; }
		//void update_lrvd_step() ;

		//void compute_approx_rvd() ;
		//void compute_exact_rvd() ;
		//void compute_local_rvd() ;

		// insert a new sample point on the facet f
		void insert_random_sample() ;

		//void insert_sample(int fidx, const vec3& p) ;
		//void insert_sample_obtuse() ;
		//void purtube_obtuse_vertices() ;
		//void remove_low_valence_vertices() ;
		//void purtube_bad_valence_vertices() ;

		void set_vertices(int N, double* x) ;
		void funcgrad(double* f, double* grad_fs, double normal_aniso) ;
		void funcgrad_constrained(int N, double* x_in, double* f, double* g) ;
		GLboolean& lock_feature() { return lock_feature_ ; }
		GLboolean& use_feature() { return use_feature_ ; }
		//dxy add
		GLboolean& lock_corners() { return lock_corners_; }
		void snap_to_corners(int N, double *x);
		void snap_to_feature(int N, double *x);
		void align_feature_grad(int N, double *g);
		//
		std::vector<bool>& on_feature() { return on_feature_ ; }
		vec3 project_to_crease(vec3& pt, vec3& tan) ;

		void clip_by_cell(
			std::vector<int>& neighbors, int c, 
			TopoPolyMesh*& ping, TopoPolyMesh*& pong, 
			bool do_push=false
			) ;
		void clip_by_plane(
			TopoPolyMesh*& ping, TopoPolyMesh*& pong, 
			const Geex::Plane<real>& P, 
			unsigned int cell_in, unsigned int cell_out,
			bool do_push=false
			) ;
		Sign side(
			const Geex::Plane<real>& P, 
			const TopoPolyVertexEdge& v,
			unsigned int cell_in, unsigned int cell_out
			) const ;
		Map*                           map()      { return map_ ; }
		//		std::vector<SamplePoint>&      samples()  { return samples_ ; }
		std::vector<std::vector<int> >& rvd_facets() { return rvd_facets_ ; }
		std::vector<Map::Facet*>&      facets()   { return facets_; }
		std::vector<int>&              facet_cell_map() { return facet_cell_map_ ; }
		MapFacetProperty<int>&         facet_index() { return facet_index_ ; }

		std::vector<Map::Vertex*>&     verts()         { return verts_; }
		std::vector<int>&              vert_cell_map() { return vert_cell_map_ ; }
		MapVertexProperty<int>&        vindex()        { return vindex_ ; }
		std::vector<std::vector<int> > rvd_verts()     { return rvd_verts_ ; }
		//		MapVertexProperty<double>&     vgeodist()      { return vgeodist_ ; }

		bool is_cliping_facet(Map::Facet* f) ;

		// TopoPolyMesh based exact RVD
		Delaunay* delaunay() { return delaunay_ ; }
		TopoPolyMesh* boundary() { return boundary_ ; }
		const TopoPolyMesh* boundary() const { return boundary_ ; }
		const Plane<real>& boundary_facet_plane(unsigned int i) const {
			return boundary_facet_plane_[i] ;
		}
		RestrictedVoronoiDiagram_poly& RVD() { 
			RVD_.set_exact(exact_ == GL_TRUE) ;
			return RVD_ ; 
		}

		const RestrictedVoronoiDiagram_poly& RVD() const { 
			const_cast<DelaunayCVT*>(this)->RVD_.set_exact(exact_ == GL_TRUE) ;
			return RVD_ ; 
		}

		size_t nb_vertices()                { return vertices_.size() ; }
		std::vector<vec3>& vertices()             { return vertices_ ; }
		std::vector<unsigned>&  valences()        { return valences_ ; }
		const std::vector<vec3>& vertices() const { return vertices_ ; }
		std::vector<bool>& locked()               { return locked_ ; }
		GLboolean& exact()                        { return exact_ ; }
		GLboolean& use_density()                  { return use_density_ ; }
		///dxy add
		GLboolean& use_face_center_init()         { return use_face_center_init_; }
		GLboolean& use_direction()				  { return use_direction_; }
		GLboolean& use_cvt_energy()				  { return use_cvt_energy_ ; }
		GLboolean& use_adaptive_mix_factor()      { return use_adaptive_mix_factor_; }
		GLboolean& use_true_gradient()            { return use_true_gradient_; }
		GLboolean& use_special_quad()             { return use_special_quad_; }
		//
		GLboolean& use_auto_save()			  { return use_auto_save_; }
		//edge direction match
		GLboolean& use_edge_dir_match()			  { return use_edge_dir_match_; }
		GLboolean& use_edge_dir_match_topo_opt()  { return use_edge_dir_match_topo_opt_; }
		///
		std::vector<int>& sample_facets()         { return sample_facets_ ; }
		AABBTree* tree()                          { return tree_ ; }
		//		std::vector<PrimalTri>& primal_facets()   { return primal_facets_ ; }
		std::vector<int>&       corners()         { return corners_ ; }
		std::vector<std::vector<int>>& features() { return features_ ;  }
		//		unsigned nb_obtuses() ;
		//		unsigned nb_small() ;
		//		unsigned is_obtuse(PrimalTri& t) ;
		double&  angle_min() { return angle_min_ ; }
		double&  angle_max() { return angle_max_ ; }

		void compute_lfs() ;
		void compute_curvature_tensor(
			bool on_facets = false,
			bool relative = true,
			bool anisotropic = false,
			int aniso_iters  = 3,
			double aniso_factor = 1.5,
			double neighborhood_size = 0.01,
			const std::string& K1_attribute_name = "K1",
			const std::string& K2_attribute_name = "K2",
			const std::string& N_attribute_name  = "N",
			const std::string& k1_attribute_name = "k1",
			const std::string& k2_attribute_name = "k2",
			const std::string& n_attribute_name  = "n",
			bool compute_hk = false,				      
			bool copy_n_to_normal = false
			) ;

		void smooth_density(Map* surface, int nb_iter) ;

		// geodesic isolines
		void extract_iso_lines() ;

		void push_symbolic_mode() {symbolic_mode_backup_ = RVD_.symbolic(); }
		void pop_symbolic_mode() {RVD_.set_symbolic(symbolic_mode_backup_); }

		void export_rvd(std::string& filename) ;

		// CVT
		//		void lloyd_local(int nb_iter = 1) ;
		//		void newton_lloyd_local(int nb_iter = 1) ;
		void lloyd_global(int nb_iter = 1) ;
		void newton_lloyd_global(int nb_iter = 1) ;
		//		void lloyd_feature(int nb_iter = 1) ;

		// Properties
		static DelaunayCVT* instance() { return instance_ ; }
		std::vector<vec3>& all_vertices() { return vertices_ ; }

		void set_vertices(const double* x) ;    
		double& normal_aniso() { return normal_aniso_ ; }
		///dxy add
		double& energy_mixfactor() { return energy_mixfactor_; }
		double& dir_para_a() { return direction_parameter_a_; }
		double& dir_para_quad() { return direction_parameter_quad_; }
		unsigned int& rot_symmetry() { return field_rot_symmetry_; }
		GLboolean& use_facet_field() { return use_facet_field_; }
		double& average_area() { return average_area_; }
		double& total_area() { return total_area_; }
		//
		double calc_map_area();
		///

		////new insertion code
		//void prepare_sort_points(int N, const double *x, std::vector<MyPoint>& mypointvec);
		//void inserte_sort_points(const std::vector<MyPoint>& mypointvec, bool tag);

		//inline GLboolean& weighted()    { return weighted_ ; }
		//inline GLboolean& constrained() { return constrained_ ; }

		//		LocalRVD* lrvd(){return lrvd_ ;}

		GLenum& optimizer_mode() { return optimizer_mode_ ; }
		float&  set_opt_m() {return opt_m_ ; }

		std::vector<std::vector<Map::Facet*>>& charts() { return charts_ ; }

		///dxy add
		std::vector<vec3> &seed_interpolation_weight() {
			return seed_interpolation_weight_;
		}

		std::vector<std::vector<vec3>> &seed_project_triangle() {
			return seed_project_triangle_;
		}

		std::vector<vec3>& seed_field() {
			return seed_field_;
		}
		void set_seed_field(const std::vector<vec3>& f_in) {
			seed_field_ = f_in;
		}

		std::vector<vec3>& seed_normal() {
			return seed_normal_;
		}
		void set_seed_normal(const std::vector<vec3>& n_in) {
			seed_normal_ = n_in;
		}

		std::map<std::pair<int, int>, std::pair<double, int>> &primal_edge_field_match() {
			return primal_edge_field_match_;
		}

		std::map<std::pair<int, int>, int> &primal_edge_match_history() {
			return primal_edge_match_history_;
		}
		//
		const std::map<int, std::set<int>> &seed_feature_lines() const { return seed_feature_lines_; }
		const std::map<int, int> &seed_feature_line() const { return seed_feature_line_; }

		//dxy add
		void save_primal_map(const std::string& filename);

		const std::vector<double>& primal_angles() {
			return primal_angles_;
		}
		void update_primal_angles();
		void save_primal_angles(const std::string& filename); //sort and save

		const std::vector<double>& primal_areas() { return primal_areas_; }
		void update_primal_areas();
		void save_primal_areas(const std::string& filename); //sort and save
		void save_primal_info(const std::string& filename);  //angles, areas, singularities...
		void print_primal_info(bool angle, bool area, bool singularity);

		GLenum& direction_energy_mode() { return direction_energy_mode_; }
		GLenum& direction_edge_weight_mode() { return direction_edge_weight_mode_; }

		//dxy add
		const std::vector<double> &rvd_cell_area() const { return rvd_cell_area_; }
		std::vector<double>	&rvd_cell_area() { return rvd_cell_area_; }

		const std::vector<double> &true_cvt_energy() const { return true_cvt_energy_; }
		std::vector<double> &true_cvt_energy()  { return true_cvt_energy_; }

		const std::vector<double> &approx_cvt_energy() const { return approx_cvt_energy_; }
		std::vector<double> &approx_cvt_energy() { return approx_cvt_energy_; }
		void compute_approx_cvt_energy(const std::vector<std::map<unsigned int, double>> & edge_energy);

		const std::vector<double> &direction_energy() const { return direction_energy_; }
		std::vector<double> &direction_energy() { return direction_energy_; }

		const std::vector<vec3>& lloyd_grad() { return lloyd_grad_; }
		void set_lloyd_grad(const std::vector<vec3>& g_in) { lloyd_grad_ = g_in; }
		void set_lloyd_grad(double* g_in) {
			lloyd_grad_.resize(nb_vertices());
			for (unsigned int i=0; i<nb_vertices(); ++i)
			{
				lloyd_grad_[i] = vec3(g_in[3*i], g_in[3*i + 1], g_in[3*i + 2]);
			}
		}

		const std::vector<vec3>& direction_grad() const { return direction_grad_; }
		std::vector<vec3>& direction_grad() { return direction_grad_; }
		void set_direction_grad(const std::vector<vec3>& g_in) { direction_grad_ = g_in; }
		void set_direction_grad(double* g_in) {
			direction_grad_.resize(nb_vertices());
			for (unsigned int i=0; i<nb_vertices(); ++i)
			{
				direction_grad_[i] = vec3(g_in[3*i], g_in[3*i + 1], g_in[3*i + 2]);
			}
		}

		//dxy add: delaunay_map
		void update_vertices_from_delaunay_map();
		void update_delaunay_map();
		void update_map_edge_stretch();
		Map* delaunay_map() { return delaunay_map_; }
		const std::map<Map::Vertex*, int>& delaunay_map_vert_index() { return delaunay_map_vert_index_; }
		const std::map<Map::Halfedge*, double>& delaunay_map_edge_stretch() { return delaunay_map_edge_stretch_; }
		//bool delaunay_is_corner(Map::Vertex *v) { return delaunay_is_corner_[v]; }
		//bool delaunay_is_on_feature(Map::Vertex *v) { return delaunay_is_on_feature_[v]; }
		//int  delaunay_feature_line(Map::Vertex *v)  { return delaunay_feature_line_[v];  }
		//topo
		void calc_singularities(int& nb);
		bool find_topo_operation(Map::Vertex_iterator v, TopoOperation& op, double dihedral_dot = -1.0);
		bool find_topo_operation_preserve_fea(Map::Vertex_iterator v, TopoOperation& op);
		bool execute_map_topo_operation(MapEditorExt& ed, TopoOperation& op, std::vector<Map::Vertex*>& to_smooth);
		void map_laplace_smooth(std::vector<Map::Vertex*>& to_smooth);
		void map_topo_optimization(int nb_iter, int deltaR_threshold = 0, double dihedral_dot_threshold = -1.0);
		unsigned int& map_topo_optimization_nb_iter() { return map_topo_optimization_nb_iter_; }
		int& map_topo_optimization_deltaR() { return map_topo_optimization_deltaR_; }
		double& dihedral_angle_threshold() { return dihedral_angle_threshold_ ; }
		//stretch opt
		typedef std::pair<Map::Halfedge_iterator, double> EdgeValue;
		struct EdgeValueGreater
		{
			bool operator() (const EdgeValue& a, const EdgeValue& b) { return a.second > b.second; }
		};
		struct EdgeValueLess
		{
			bool operator() (const EdgeValue& a, const EdgeValue& b) { return a.second < b.second; }
		};
		GLboolean& use_stretch_topo_optimize() { return use_stretch_topo_optimize_; }
		GLfloat& stretch_optimize_select_rate() { return stretch_optimize_select_rate_; }
		GLint& stretch_optimize_nb_smooth_iter() { return stretch_optimize_nb_smooth_iter_; }
		void collect_stretched_edges(std::vector<EdgeValue>&);
		void stretch_topo_optimize(unsigned int nb_iter);
		bool quad_dominant_topo_optimize(); //return true if some topo operations are found and executed
		//dxy add: FCVT
		void computeDirectionEnergyAndGradient(const std::vector<vec3> &N, const std::vector<vec3> &D, const std::vector<std::map<unsigned int, double>> &edge_weight, unsigned int RoSy);
		void compute_edge_dir_match(const std::vector<vec3> &N, const std::vector<vec3> &D, const std::vector<std::map<unsigned int, double>> &edge_weight, unsigned int RoSy);
		void computeDirectionEnergyAndGradient_match(const std::vector<vec3> &N, const std::vector<vec3> &D, const std::vector<std::map<unsigned int, double>> &edge_weight, unsigned int RoSy);
		void computeDirectionEnergyAndGradient(double &f, double *g, const std::vector<vec3> &N, const std::vector<vec3> &D, 
			const std::vector<std::map<unsigned int, double>> &edge_weight, unsigned int RoSy);
		void  computeEdgeDirectionEG(double &f, vec3 &g, const vec3& X1, const vec3& X2, const vec3& n, const vec3& d, unsigned int RoSy = 6 );
		void  computeEdgeDirectionEG(double &f, vec3 &g, const vec3& X1X2, const vec3& n, const vec3& d, unsigned int RoSy = 6 );
		//dxy add: interpolation
		GLboolean& use_new_interpolation() { return use_new_interpolation_; }
		void computeNormalAndDirection(std::vector<vec3> &N, std::vector<vec3> &D, std::vector<unsigned int> &Idx, unsigned int RoSy = 6);
		void computeWeight(std::vector<double> &W, const std::vector<unsigned int> &Idx);
		void computeWeight(std::vector<double> &W, const std::vector<std::map<unsigned int, std::vector<std::pair<Geex::TopoPolyVertexEdge, Geex::TopoPolyVertexEdge>>>> &dual_segments);


		void compute_dual_length(std::vector<std::map<unsigned int, double>> &dual_length, const std::vector<std::map<unsigned int, std::vector<std::pair<TopoPolyVertexEdge, TopoPolyVertexEdge>>>> &dual_segments);
		void compute_approximated_lloyd_energy(std::vector<std::map<unsigned int, double>> &energy, const std::vector<std::map<unsigned int, std::vector<std::pair<TopoPolyVertexEdge, TopoPolyVertexEdge>>>> &dual_segments, const std::vector<double> &v_w);


		static void normal_field_interpolate(vec3 &vn, vec3 &vf, const vec3 &v, const vec3 pos[3], const vec3 normal[3], const vec3 field[3], unsigned int RoSy =6);
		//old way of interpolation
		static void simple_normal_field_interpolate(vec3 &vn, vec3 &vf, const vec3 &v, const vec3 pos[3], const vec3 normal[3], const vec3 field[3], unsigned int RoSy =6);
		static void compute_vertex_weight(vec3 &weight, const vec3 &p, const vec3 vert[3], GLenum mode = Barycentric);  //used for interpolation
		static GLenum& interpolationMode() { return interpolationMode_; }
		static vec3 vec_rotate(const vec3 &v, const vec3 &axie, double angle)
		{
			vec3 ez = normalize(axie);
			return cos(angle)*v + sin(angle)*cross(ez, v);
		}

		static void compute_barycentric_weight(vec3 &weight, const vec3 &p, const vec3 vert[3])
		{
			vec3 pc = p - vert[2];
			vec3 ac = vert[0] - vert[2];
			vec3 bc = vert[1] - vert[2];
			double a00 = ac.length2();
			double a01 = dot(ac, bc);
			double a11 = bc.length2();
			double b0 = dot(pc, ac);
			double b1 = dot(pc, bc);
			double mdet = a00 * a11 - a01 * a01;
			double s = (a11 * b0 - a01 * b1) / mdet;
			double t = (a00 * b1 - a01 * b0) / mdet;
			weight = vec3(s, t, 1-s-t);
			//proj = s * a + t * b + (1 - s - t ) * c;
		}

		static void compute_cotangent_weight(vec3 &weight, const vec3 &p, const vec3 vert[3])
		{
			vec3 ab = vert[1] - vert[0];
			vec3 bc = vert[2] - vert[1];
			vec3 ca = vert[0] - vert[2];
			vec3 ap = p - vert[0];
			vec3 bp = p - vert[1];
			vec3 cp = p - vert[2];
			vec3 n = cross(ab, bc);
			normalize(n);
			//refer to OneNote 4.23 for the meaning of notation EA, EB...
			double EA = dot(ab, ap);
			double EB = dot(-ab, bp);
			double ED = dot(cross(ap, bp), n);
			double FA = dot(-ca, ap);
			double FC = dot(ca, cp);
			double FD = dot(cross(cp, ap), n);
			double GB = dot(bc, bp);
			double GC = dot(-bc, cp);
			double GD = dot(cross(bp, cp), n);
			//calc weight
			if (abs(ED) < 1e-10)
			{
				weight = vec3(EB/(EA+EB), EA/(EA+EB), 0);
			}
			else if (abs(GD) < 1e-10) {
				weight = vec3(0, GC/(GB+GC), GB/(GB+GC));
			}
			else if (abs(FD) < 1e-10) {
				weight = vec3(FC/(FA+FC), 0, FA/(FA+FC));
			}
			else {
				double wa = EB/ED + FC/FD;
				double wb = EA/ED + GC/GD;
				double wc = FA/FD + GB/GD;
				double sum = wa + wb + wc;
				weight = vec3(wa/sum, wb/sum, wc/sum);
			}
		}

		//
		double& model_normalize_scale() { return model_normalize_scale_; }
		//
		void compute_feature_edge_line();
		void compute_seed_feature_lines();
		void update_seed_feature_line();
		//
		const std::vector<Experiment> &experiments() { return experiments_; }
		int& experiment_iter() { return experiment_iter_; }
		///
		//dxy test
		static GLenum interpolationMode_;

		void check_map();

		//dxy add: multi-thread
		void init_multithread();
		GLboolean &use_multithread() { return use_multithread_; }
		CVTMultiThread *multithread() { return multithread_; }
		///


	private:
		///dxy add
		double model_normalize_scale_;

		//multi-thread
		GLboolean use_multithread_;
		CVTMultiThread *multithread_;
		//

		GLboolean field_loaded_;
		std::vector<vec3> seed_field_;
		std::vector<vec3> seed_normal_;
		std::vector<vec3> seed_interpolation_weight_;
		std::vector<std::vector<vec3>> seed_project_triangle_;
		std::map<std::pair<int, int>, std::pair<double, int>> primal_edge_field_match_; //map: <vidx, vidx> -> <dir_angle, dir_order>
		std::map<std::pair<int, int>, int> primal_edge_match_history_;
		//
		std::map<std::pair<int, int>, int> feature_edge_line_;
		std::map<int, std::set<int>> seed_feature_lines_;
		std::map<int, int> seed_feature_line_;
		std::map<int, vec3> seed_feature_tan_;
		//
		std::vector<double> primal_angles_;
		std::vector<double> primal_min_angles_; //min angle per triangle
		std::vector<double> primal_max_angles_; //max angle per triangle
		std::vector<double> primal_areas_;
		GLenum direction_energy_mode_;
		GLenum direction_edge_weight_mode_;
		//
		GLboolean use_new_interpolation_;
		//
		double total_area_;
		double average_area_;

		std::vector<vec3> lloyd_grad_;
		std::vector<vec3> direction_grad_;
		std::vector<double> rvd_cell_area_;
		std::vector<double> true_cvt_energy_;
		std::vector<double> approx_cvt_energy_;
		std::vector<double> direction_energy_;
		///

		//		LocalRVD*        lrvd_ ;
		static DelaunayCVT* instance_ ;
		//		GLboolean        weighted_ ;
		GLboolean        constrained_ ;
		double           normal_aniso_ ;
		///dxy add
		double			 energy_mixfactor_;
		double			 direction_parameter_a_;
		double			 direction_parameter_quad_;
		unsigned int	 field_rot_symmetry_;
		GLboolean		 use_facet_field_;
		///

		GLenum           optimizer_mode_ ;
		float            opt_m_;

		// class members of approximate rvd
		std::vector<vec3>             vertices_ ;
		std::vector<int>              sample_facets_ ;  // the triangles that contain each samples
		std::vector<bool>             locked_ ;
		std::vector<unsigned>		  valences_ ;
		std::map<unsigned, std::set<int>> adjtris_ ;
		MedialAxis* lfs_;
		//dxy add
		std::vector<int>			  vertices_corner_index_;
		//

		Map*                          map_ ;
		std::vector<Map::Facet*>      facets_ ;
		std::vector<Map::Vertex*>     verts_ ;
		MapFacetProperty<int>         facet_index_ ;
		MapFacetProperty<int>         fchart_ ;		
		MapVertexProperty<int>        vindex_ ;
		MapVertexProperty<bool>       is_corner_ ;
		MapHalfedgeProperty<bool>     is_feature_ ;

		std::vector<std::vector<int> > rvd_facets_ ;
		std::vector<std::vector<int> > rvd_verts_ ;
		std::vector<int>               facet_cell_map_ ;
		std::vector<int>               vert_cell_map_ ;
		std::vector<std::vector<int> > cliping_facets_ ;
		//std::vector<std::set<int>>    facet_cells_ ;
		std::vector<std::vector<int> > facet_cells_ ;

		//		std::vector<PrimalTri>         primal_facets_ ;
		GLboolean                      primal_dirty_ ;

		// feature handling
		std::vector<int>               corners_ ;
		std::vector<std::vector<int> > features_ ;  
		GLboolean                      lock_feature_ ;
		std::vector<bool>              on_feature_ ;
		GLboolean                      use_feature_ ;
		//std::vector<bool> locked_;
		//dxy add
		GLboolean					   lock_corners_;
		//

		// class members of exact rvd
		TopoPolyMesh*                  boundary_ ;
		std::vector< Plane<real> >     boundary_facet_plane_ ;
		std::string                    del_algo_ ;
		Delaunay*                      delaunay_ ;
		bool						   symbolic_mode_backup_;
		RestrictedVoronoiDiagram_poly  RVD_ ;
		GLboolean                      exact_ ;
		GLboolean                      symbolic_ ;
		GLboolean                      use_density_ ;
		///dxy add
		GLboolean					   use_face_center_init_;
		//
		GLboolean					   use_edge_dir_match_;
		GLboolean					   use_edge_dir_match_topo_opt_;
		//
		GLboolean					   use_direction_ ;
		GLboolean					   use_cvt_energy_;
		GLboolean					   use_adaptive_mix_factor_;
		GLboolean					   use_true_gradient_;
		GLboolean					   use_special_quad_;
		//
		GLboolean					   use_auto_save_;
		///

		TopoPolyStack                  S_ ;
		TopoPolyMesh                   M1, M2 ;
		std::vector<int>               edge_bisectors1_ ;
		std::vector<int>               edge_bisectors2_ ;

		real                           x_min_, y_min_, z_min_, x_max_, y_max_, z_max_ ;
		bool                           bbox_dirty_ ;
		AABBTree*                      tree_ ;
		unsigned int                   lrvd_cur_step_ ;

		// affinity matrix (use vector for test. to be optimized
		//std::vector<std::vector<double>> affinity_ ;
		SparseMatrix                     affmat_ ;
		double                           affinity_weight_ ;
		double							 cvt_normalization_ ;
		double                           affinity_normalization_ ;

		double                          angle_min_, angle_max_ ;
		unsigned int                    nb_lock_rings_ ;
		std::vector<std::vector<Map::Facet*>> charts_ ;
		std::ofstream     obtuse_out ;

		//dxy add: Map of delaunay_
		bool is_delaunay_map_dirty_;
		bool is_map_edge_stretch_dirty_;
		bool is_primal_angles_dirty_;
		bool is_primal_areas_dirty_;

		Map* delaunay_map_;
		std::map<Map::Vertex*, int> delaunay_map_vert_index_;
		std::map<Map::Halfedge*, double> delaunay_map_edge_stretch_;
		//feature preservation
		//MapVertexProperty<bool>		  delaunay_is_corner_;
		//MapVertexProperty<bool>		  delaunay_is_on_feature_;
		//MapVertexProperty<int>		  delaunay_feature_line_;
		//stretch optimize
		GLboolean use_stretch_topo_optimize_;
		GLfloat stretch_optimize_select_rate_;
		GLint stretch_optimize_nb_smooth_iter_;
		//map topo optimize
		unsigned int map_topo_optimization_nb_iter_;
		int map_topo_optimization_deltaR_;
		double dihedral_angle_threshold_;
		//automotive experiment
		std::vector<Experiment> experiments_;
		int experiment_iter_;
		///

	} ;

	//dxy add
	//Experiment : description of an experiment
	class Experiment
	{
	public:
		Experiment() : num(-1) , pts(-1), nb_iter(10000), mix_factor(0), use_density(false), 
			use_topo_opt(false), nb_topo_opt_iter(10), use_quad(false),
			use_all_edge_split(false) {}; //null Experiment
		//Experiment() : num(-1), pts(-1), mixfactor(-1), use_density(false) {};
		//Experiment(int n, int p, double factor, bool density = false, bool quad=false, double pq=0.5)
			//: num(n), pts(p), mixfactor(factor), use_density(density), use_quad(quad), para_quad(pq) {};
		int num;   //order number of the experiment
		int pts;
		//
		int nb_iter;
		double mix_factor;
		bool use_density;
		//
		bool use_topo_opt;
		int nb_topo_opt_iter;
		//quad
		bool use_quad;
		//double para_quad;
		//
		bool use_all_edge_split;

	};

	//Topo Operation {EF, EC, ES, VS} related to a vertex
	class TopoOperation {
	public:
		TopoOperation() : _v(0), _e(0), _type(Null_Operator), _delta_R(0) {}
		TopoOperation(Map::Vertex* v, Map::Halfedge* e, TopoOperator t, int dR) : _v(v), _e(e), _type(t), _delta_R(dR) {}
		TopoOperation(const TopoOperation& rhs) : _v(rhs._v), _e(rhs._e), _type(rhs._type), _delta_R(rhs._delta_R) {}

		Map::Vertex* vertex() const { return _v; }
		Map::Halfedge* edge() const { return _e; }
		TopoOperator type() const { return _type; }
		int delta_R() const { return _delta_R; }

		Map::Halfedge* counter_edge() const {
			int dv = _v->degree();
			int shift = (dv%2) ? ((dv-1)/2) : (dv/2);
			Map::Halfedge* e = _e;
			while (shift>0)
			{
				e = e->next_around_vertex();
				shift--;
			}
			return e;
		}

	private:
		Map::Vertex* _v;
		Map::Halfedge* _e;
		TopoOperator _type;
		int _delta_R;

		friend bool DelaunayCVT::find_topo_operation(Map::Vertex_iterator v, TopoOperation& op, double dihedral_dot);
		friend bool DelaunayCVT::find_topo_operation_preserve_fea(Map::Vertex_iterator v, TopoOperation& op);
	};

	//MatchMaker: match two set of directions(angles)
	class MatchMaker {
	public:
		typedef std::pair<double, int> Direction; //<angle, order> pair

		MatchMaker(unsigned int rosy_in, int hi_order_in) : rosy_(rosy_in), hi_order_(hi_order_in) {
			gx_assert(rosy_ > 0);
			dynamic_order_ = (hi_order_ < 0);
			gx_assert(!dynamic_order_);  //not supported now
			//dirs init
			init();
			//test
// 			std::cout << "MatchMaker:" << std::endl;
// 			std::cout << "number of directions on each order: ";
// 			for (int i=0; i<hi_order_; ++i)
// 			{
// 				std::cout << tot_count_[i] << " ";
// 			}
// 			std::cout << std::endl;
		}

		~MatchMaker() {
			delete [] tot_count_;
		}

		//void match(const std::vector<double> &edge_angle, std::vector<double> &rot_angle) {};


		/**
		* match angles in edge_angle and directions in dirs_
		* Param:
		* edge_angle: sorted angle in range [0, 2pi)
		* matched_dirs: the computed matching directions in the same order as angles in angle_edge
		* match_history: the order in which the edges are matched
		*/
		void match(const std::vector<double> &edge_angle, std::vector<Direction> &matched_dirs, std::vector<int> &match_history) {
			unsigned int nedge = edge_angle.size();
			unsigned int ndir = dirs_.size();
			gx_assert(nedge>0 && nedge<=ndir);

			//find first match between edge_angles and 0-order directions, both in range [0, 2pi)
			double min_diff = 2 * M_PI;
			unsigned int min_eidx, min_didx;
			bool match_forward_edge = false;
			for (int i=0; i<dirs_.size(); ++i)
			{
				if (dirs_[i].second == 0) //0-order direction
				{
					double dir_angle = dirs_[i].first;
					auto lb = std::lower_bound(edge_angle.begin(), edge_angle.end(), dir_angle);
					if (lb==edge_angle.begin()) //the most close angle may be the first or the last one
					{
						double d1 = *lb - dir_angle;
						double d2 = 2 * M_PI + dir_angle - edge_angle.back();
						gx_assert(d1>=0 && d1<2*M_PI && d2>=0 && d2<2*M_PI);
						if (d1<=d2) {
							if (d1 < min_diff) {
								min_diff = d1;
								min_eidx = lb - edge_angle.begin();
								min_didx = i;
								match_forward_edge = true;
							}							
						}
						else if (d2 < min_diff) {
							min_diff = d2;
							min_eidx = edge_angle.size()-1;
							min_didx = i;
						}
					}
					else if (lb==edge_angle.end()) //the most close angle may be the first or the last one
					{
						double d1 = edge_angle.front() + 2*M_PI - dir_angle;
						double d2 = dir_angle - edge_angle.back();
						gx_assert(d1>=0 && d1<2*M_PI && d2>=0 && d2<2*M_PI);
						if (d1<=d2) {
							if (d1 < min_diff) {
								min_diff = d1;
								min_eidx = 0;
								min_didx = i;
								match_forward_edge = true;
							}							
						}
						else if (d2 < min_diff) {
							min_diff = d2;
							min_eidx = edge_angle.size()-1;
							min_didx = i;
						}
					}
					else  //the most close angle may be lower_bound or upper_bound
					{
						auto ub = lb;
						--ub;
						double d1 = *lb - dir_angle;
						double d2 = dir_angle - *ub;
						gx_assert(d1>=0 && d1<2*M_PI && d2>=0 && d2<2*M_PI);
						if (d1<=d2) {
							if (d1 < min_diff) {
								min_diff = d1;
								min_eidx = lb - edge_angle.begin();
								min_didx = i;
								match_forward_edge = true;
							}							
						}
						else if (d2 < min_diff) {
							min_diff = d2;
							min_eidx = ub - edge_angle.begin();
							min_didx = i;
						}
					}
				}
			}

			//move the first matched edge and direction to the front , and translate angle to [0,2pi)
			std::vector<double> edges;
			std::vector<Direction> dirs;
			for (unsigned int i=0; i<nedge; ++i)
			{
				unsigned int cur_eidx = i+min_eidx;
				double cur_ang;
				if (cur_eidx < nedge) cur_ang = edge_angle[cur_eidx];
				else cur_ang = 2*M_PI + edge_angle[cur_eidx % nedge];
				cur_ang -= dirs_[min_didx].first;
				edges.push_back(cur_ang);
				if (i>0) {
					gx_assert(edges[i] >= edges[i-1]);
					gx_assert(edges[i]>0 && (edges[i]-edges[0]<2*M_PI));
				}
			}
			if ((!match_forward_edge) && edges[0]>0) {
				for (int i=0; i<edges.size(); ++i) edges[i] -= 2*M_PI;
			}
			gx_assert(edges[0]>=-M_PI/rosy_);
			gx_assert(edges[0]<=M_PI/rosy_);
			
			for (int unsigned i=0; i<ndir; ++i)  {
				unsigned int cur_didx = i+min_didx;
				Direction cur_d = dirs_[cur_didx % ndir];
				if (cur_didx >= ndir) cur_d.first += 2*M_PI;
				cur_d.first -= dirs_[min_didx].first;
				dirs.push_back(cur_d);
				gx_assert(cur_d.first>=0 && cur_d.first<2*M_PI);
				if (i>0) {
					gx_assert(dirs[i].first >= dirs[i-1].first);
				}
			}

			//define inner_match function
			std::vector<std::pair<unsigned int, unsigned int>> matches;
			int hi_order = hi_order_;

			std::function<bool(unsigned int, unsigned int, unsigned int, unsigned int, int**)> inner_match;
			inner_match = [&edges, &dirs, &matches, hi_order, &inner_match] (unsigned int ebegin, unsigned int eend, unsigned int dbegin, unsigned int dend,
				int **NR /*int *to_match, int *can_match*/) -> bool {
					unsigned int ne = eend - ebegin;
					unsigned int nd = dend - dbegin;
					gx_assert(ne>=0 && nd>=0);
					bool succeed;
					if (ne==0) {
						succeed = true;
					}
					else if (ne > nd) {
						succeed = false;
					}
					else if (ne==nd) { //match one by one
						for (unsigned int i=0; i<ne; ++i) matches.push_back(std::make_pair(ebegin+i, dbegin+i));
						succeed = true;
					}
					else {	//0 < ne < nd
						int *to_match = NR[0];
						int *can_match = NR[1];
						bool can_full_match = true;
						for (int i=0; i<hi_order; ++i)
						{
							gx_assert(to_match[i]>=0 && can_match[i]>=0);
							gx_assert(to_match[i] <= can_match[i]);
							if ((to_match[i]>0) && (to_match[i]<can_match[i])) {
								can_full_match = false;
								break;
							}
						}
						if (can_full_match)  //simple case: no free choice
						{
							int ei = ebegin;
							for (unsigned int i=dbegin; i<dend; ++i) {
								int ord = dirs[i].second;
								if (to_match[ord]>0) {
									matches.push_back(std::make_pair(ei, i));
									++ei;
									to_match[ord] -= 1;
								}
							}
							gx_assert(ei==eend);
							succeed = true;
						}
						else  //complex case: allow free choice
						{
							int cur_order = 0;
							while (cur_order<hi_order && to_match[cur_order]==0) {
								gx_assert(can_match[cur_order]==0);
								++cur_order;
							}
							if (to_match[cur_order] == can_match[cur_order]) //the next cur-order direction must be matched
							{
								//unsigned int d_cnt = 0;
								int *ord_cnt = new int[hi_order];
								memset(ord_cnt, 0, sizeof(int)*hi_order);
								int dnext = dbegin;
								while (dirs[dnext].second!=cur_order)
								{
									//++d_cnt;
									ord_cnt[dirs[dnext].second] += 1;
									++dnext;
								}
								//++d_cnt;
								++ord_cnt[dirs[dnext].second];
								//
								unsigned min_match_cnt = 0;
								unsigned max_match_cnt = 0;
								for (unsigned int ord=cur_order; ord<hi_order; ++ord)
								{
									int ord_min = to_match[ord] + ord_cnt[ord] - can_match[ord];
									int ord_max = (to_match[ord] < ord_cnt[ord]) ? to_match[ord] : ord_cnt[ord];
									min_match_cnt += ((ord_min>0) ? ord_min : 0);
									max_match_cnt += ord_max;
								}
								gx_assert(min_match_cnt<=max_match_cnt);
								//find matching edge
								double cur_d_ang = dirs[dnext].first;
								unsigned int search_begin = ebegin+min_match_cnt-1;
								unsigned int search_end = ebegin+max_match_cnt;
								gx_assert(search_end<=eend);
								auto lb = std::lower_bound(edges.begin()+search_begin, edges.begin()+search_end, cur_d_ang);
								unsigned int enext;
								if (lb-edges.begin() == search_begin) {
									enext = search_begin;
								}
								else if (lb-edges.begin() == search_end) {
									enext = search_end-1;
								}
								else {
									auto ub = lb;
									--ub;
									double d1 = *lb - cur_d_ang;
									double d2 = cur_d_ang - *ub;
									gx_assert(d1>=0 && d2>=0);
									gx_assert(d1<2*M_PI && d2<2*M_PI);
									if (d1<=d2) {
										enext = lb-edges.begin();
									} else {
										enext = ub-edges.begin();
									}
								}
								matches.push_back(std::make_pair(enext, dnext));
								//prepare R (can_match)
								//int *R1 = ord_cnt;
								int *R1 = new int[hi_order];
								memcpy(R1, ord_cnt, sizeof(int)*hi_order);
								delete [] ord_cnt;
								int *R2 = new int[hi_order];
								memset(R2, 0, sizeof(int)*hi_order);
								for (unsigned int ord=cur_order; ord<hi_order; ++ord)
								{
									R2[ord] = can_match[ord] - R1[ord];
								}
								R1[cur_order] -= 1;
								gx_assert(R1[cur_order]==0);
								//prepare N (to_match)
								unsigned int ne1 = enext - ebegin;
								int *N1 = new int[hi_order];
								memset(N1, 0, sizeof(int)*hi_order);
								for (unsigned int ord=cur_order; ord<hi_order; ++ord)
								{
									N1[ord] = (ne1 < R1[ord]) ? ne1 : R1[ord];
									ne1 -= N1[ord];
								}
								gx_assert(ne1==0);
								unsigned int ne2 = eend - (enext+1);
								int *N2 = new int[hi_order];
								memset(N2, 0, sizeof(int)*hi_order);
								for (unsigned int ord=cur_order; ord<hi_order; ++ord)
								{
									N2[ord] = (ne2 < R2[ord]) ? ne2 : R2[ord];
									ne2 -= N2[ord];
									if (ord>cur_order) {
										gx_assert(N1[ord]+N2[ord] == to_match[ord]);
									} else {
										gx_assert(N1[ord]+N2[ord]+1 == to_match[ord]);
									}
								}
								gx_assert(ne2==0);
								//recursion
								int **NR1 = new int*[2];
								NR1[0] = N1;
								NR1[1] = R1;
								int **NR2 = new int*[2];
								NR2[0] = N2;
								NR2[1] = R2;
								succeed = (inner_match(ebegin, enext, dbegin, dnext, NR1) && inner_match(enext+1, eend, dnext+1, dend, NR2));
							}
  							else  //cur-order is the last order remains to match
  							{
								gx_assert(to_match[cur_order]==ne);
 								//assume: the front edge has the highest priority
 								//double e0 = edges[ebegin];
								std::vector<double> dir_angle;
								std::vector<unsigned int> dir_idx;
								for (int i=dbegin; i<dend; ++i) if (dirs[i].second == cur_order)
								{
									dir_angle.push_back(dirs[i].first);
									dir_idx.push_back(i);
								}
								//
								//unsigned int min_match_pos = 0;
								unsigned int max_match_pos = can_match[cur_order] - to_match[cur_order];
								auto lb = std::lower_bound(dir_angle.begin(), dir_angle.begin()+max_match_pos+1, edges[ebegin]);
								unsigned int new_dbegin;
								unsigned int used_dir_cnt;
 								if (lb==dir_angle.begin()+max_match_pos+1)
 								{
									matches.push_back(std::make_pair(ebegin, dir_idx[max_match_pos]));
									new_dbegin = dir_idx[max_match_pos] + 1;	
									used_dir_cnt = max_match_pos + 1;
 								}
								else if (lb==dir_angle.begin()) {
									matches.push_back(std::make_pair(ebegin, dir_idx[0]));
									new_dbegin = dir_idx[0] + 1;
									used_dir_cnt = 1;
								}
								else {
									auto ub=lb;
									--ub;
									double d1 = *lb - edges[ebegin];
									double d2 = edges[ebegin] - *ub;
									gx_assert(d1>=0 && d2>=0);
									gx_assert(d1<2*M_PI && d2<2*M_PI);
									if (d1 <= d2)
									{
										matches.push_back(std::make_pair(ebegin, dir_idx[lb-dir_angle.begin()]));
										new_dbegin = dir_idx[lb-dir_angle.begin()] + 1;
										used_dir_cnt = lb-dir_angle.begin() + 1;
									}
									else
									{
										matches.push_back(std::make_pair(ebegin, dir_idx[ub-dir_angle.begin()]));
										new_dbegin = dir_idx[lb-dir_angle.begin()];
										used_dir_cnt = ub-dir_angle.begin() + 1;
									}
								}
								//recursion
								//int *N1 = new int[hi_order];
								//memcpy(N1, to_match, sizeof(int)*hi_order);
								// R1 (can_match)
								int *R1 = new int[hi_order];
								memset(R1, 0, sizeof(int)*hi_order);
								for (int i=new_dbegin; i<dend; ++i)
								{
									R1[dirs[i].second] += 1;
								}
								for (int ord=0; ord<cur_order; ++ord) gx_assert(R1[ord]==0);
								gx_assert(R1[cur_order]+used_dir_cnt == can_match[cur_order]);
								// N1 (to_match)
								int *N1 = new int[hi_order];
								memset(N1, 0, sizeof(int)*hi_order);
								unsigned int ne1 = eend-ebegin-1;
								for (int ord=cur_order; ord<hi_order; ++ord)
								{
									N1[ord] = (ne1 < R1[ord]) ? ne1 : R1[ord];
									ne1 -= N1[ord];
								}
								gx_assert(ne1==0);
								//memcpy(R1, can_match, sizeof(int)*hi_order);
								//N1[cur_order] -= 1;
								//R1[cur_order] -= used_dir_cnt;
								int **NR1 = new int *[2];
								NR1[0] = N1;
								NR1[1] = R1;
								succeed = inner_match(ebegin+1, eend, new_dbegin, dend, NR1);
  							}
						}
					}

					//release memory
					delete [] NR[0];
					delete [] NR[1];
					delete [] NR;
					//
					return succeed;
			};

			//execute match
			int *R = new int[hi_order_];
			memcpy(R, tot_count_, sizeof(int)*hi_order_);
			int *N = new int[hi_order_];
			memset(N, 0, sizeof(int)*hi_order_);
			int ne = nedge;
			for (unsigned int ord=0; ord<hi_order_; ++ord)
			{
				N[ord] = (ne < R[ord]) ? ne : R[ord];
				ne -= N[ord];
			}
			gx_assert(ne==0);
			int **NR = new int *[2];
			NR[0] = N;
			NR[1] = R;
			//
			bool match_succeed = inner_match(0, nedge, 0, ndir, NR);
			gx_assert(match_succeed);

			//recover matches
			matched_dirs.clear();
			matched_dirs.resize(nedge);
			match_history.clear();
			match_history.resize(nedge);
			gx_assert(matches.size()==nedge); //every edge has a match
			for (auto mit=matches.begin(); mit!=matches.end(); ++mit)
			{
				int eidx = (mit->first + min_eidx) % nedge;  //convert to index of edge_angle
				int didx = (mit->second + min_didx) % ndir;  //convert to index of dirs_
				matched_dirs[eidx] = dirs_[didx];
				match_history[eidx] = mit-matches.begin();  //eidx is i-th to be matched
			}
		}


	private:
		unsigned int rosy_;
		int hi_order_;
		bool dynamic_order_;
		std::vector<Direction> dirs_;
		//std::vector<int> tot_count;
		int *tot_count_;



		void init() {
			//tot_count.assign(hi_order, 0);
			tot_count_ = new int[hi_order_];
			memset(tot_count_, 0, sizeof(int)*(hi_order_));

			std::list<Direction> ldirs;
			//0-order
			double da = 2 * M_PI / rosy_;
			for (unsigned int i=0; i<rosy_; ++i)
			{
				ldirs.push_back(std::make_pair(da*i, 0));
				tot_count_[0] += 1;
			}
			//higher orders
			for (int ord=1; ord<hi_order_; ++ord)
			{
				auto it = ldirs.begin();
				auto jt = it;
				++jt;
				while (jt != ldirs.end())
				{
					ldirs.insert(jt, std::make_pair(0.5*(it->first + jt->first), ord));
					tot_count_[ord] += 1;
					it = jt;
					++jt;
				}
				jt = ldirs.begin();
				ldirs.push_back(std::make_pair(0.5*(it->first + jt->first + 2*M_PI), ord));
				//ldirs.insert(it, std::make_pair(0.5*(it->first + jt->first + 2*M_PI), ord));
				tot_count_[ord] += 1;
			}
			//convert to vector
			for (auto it=ldirs.begin(); it!=ldirs.end(); ++it)
			{
				dirs_.push_back(*it);
				gx_assert(it->first>=0 && it->first<2*M_PI);
				gx_assert(dirs_.size()==1 || (dirs_.back().first >= dirs_[dirs_.size()-2].first));
			}
		}
	};


	class QDmeshEditor {
	public:
		typedef std::pair<int, int> Edge;

		QDmeshEditor(const std::map<std::pair<int, int>, std::pair<double, int>> &edge_dir_match_in,
			const std::vector<vec3> &vertices_in) : edge_dir_match(edge_dir_match_in), vertices(vertices_in), broken_lack_edge_collected(false) {
			init();
		};

		void collect_broken_lack_edges() {
			if (broken_lack_edge_collected) return;
			//brokenEdges
			Edge nullEdge(-1,-1);
			for (int i=0; i<vertices.size(); ++i) {
				for (auto it=OpEdges[i].begin(); it!=OpEdges[i].end(); ++it) {
					int j = it->first;
					if (it->second==-1) brokenEdges.insert(std::make_pair(Edge(i, j), nullEdge));
				}
			}


			//lackEdges
			for (auto eit=brokenEdges.begin(); eit!=brokenEdges.end(); ++eit) {
				//find e1
				int i = eit->first.first;
				int j = eit->first.second;
				int idx = QEdges[i][j];
				for (int k=0; k<3; ++k) {
					//ccw rotate
					idx = (idx+1) % 4;
					j = QD[i][idx];
					//flip
					idx = OpEdges[i][j];
					if (idx==-1) break;
					else {
						gx_assert(QD[j][idx]==i);
						i = j;
						j = QD[i][idx];
					}
				}
				if (idx==-1) continue;;
				Edge e1(i, j);
				//find e2
				i = eit->first.first;
				j = eit->first.second;
				idx = QEdges[i][j];
				for (int k=0; k<3; ++k) {
					//cw rotate
					idx = (idx+4-1) % 4;
					j = QD[i][idx];
					//flip
					idx = OpEdges[i][j];
					if (idx==-1) break;
					else {
						gx_assert(QD[j][idx]==i);
						i = j;
						j = QD[i][idx];
					}
				}
				if (idx==-1) continue;;
				Edge e2(i, j);
				//check and insert
				if (e1.first==e2.second && e1.second==e2.first) {
					eit->second = e1;
					lackEdges.insert(std::make_pair(e1, eit->first));
				}			
			}

			broken_lack_edge_collected = true;
		}

		void collect_topo_operations(std::vector<vec3> &to_insert, std::set<int> &to_remove) {
			if (!broken_lack_edge_collected) collect_broken_lack_edges();
			//
			to_insert.clear();
			to_remove.clear();
			std::set<Edge> usedEdges;

			//vertex remove and edge collapse
			Edge nullEdge(-1,-1);
			for (auto eit=brokenEdges.begin(); eit!=brokenEdges.end(); ++eit)
			{
				if (eit->second==nullEdge) continue;
				Edge eij = eit->first;
				if (usedEdges.find(eij)==usedEdges.end()) {
					int i = eij.first;
					int j = eij.second;
					int idx = QEdges[i][j];
					Edge e1(i, QD[i][(idx+2) % 4]);
					if ((usedEdges.find(e1)==usedEdges.end()) && (brokenEdges.find(e1)!=brokenEdges.end()) && (brokenEdges.at(e1)!=nullEdge)) {
						to_remove.insert(eij.first);
						usedEdges.insert(eij);
						usedEdges.insert(e1);
						//mark corresponding lack edges
						usedEdges.insert(brokenEdges.at(eij));
						usedEdges.insert(brokenEdges.at(e1));
					}
					else { //e1 is not available, check neighbors
						//left neighbor
						bool left_succeed = false;
						i = eij.first;
						j = eij.second;
						idx = QEdges[i][j];
						//ccw rotate
						idx = (idx+1) % 4;
						j = QD[i][idx];
						//flip
						idx = OpEdges[i][j];
						if (idx!=-1) {
							i = j;
							j = QD[i][idx];
							//cw rotate
							idx = (idx+4-1) % 4;
							j = QD[i][idx];
							Edge left_e(i, j);
							if ((usedEdges.find(left_e)==usedEdges.end()) && (brokenEdges.find(left_e)!=brokenEdges.end()) && (brokenEdges.at(left_e)!=nullEdge))
							{
								to_remove.insert(eij.first);
								to_remove.insert(left_e.first);
								to_insert.push_back(0.5 * (vertices[eij.first] + vertices[left_e.first]));
								usedEdges.insert(eij);
								usedEdges.insert(left_e);
								//mark corresponding lack edges
								usedEdges.insert(brokenEdges.at(eij));
								usedEdges.insert(brokenEdges.at(left_e));
								left_succeed = true;
							}
						}
						if (!left_succeed) { //try right neighbor
							i = eij.first;
							j = eij.second;
							idx = QEdges[i][j];
							//cw rotate
							idx = (idx+4-1) % 4;
							j = QD[i][idx];
							//flip
							idx = OpEdges[i][j];
							if (idx!=-1) {
								i = j;
								j = QD[i][idx];
								//ccw rotate
								idx = (idx+1) % 4;
								j = QD[i][idx];
								Edge right_e(i, j);
								if ((usedEdges.find(right_e)==usedEdges.end()) && (brokenEdges.find(right_e)!=brokenEdges.end()) && (brokenEdges.at(right_e)!=nullEdge)) {
									to_remove.insert(eij.first);
									to_remove.insert(right_e.first);
									to_insert.push_back(0.5 * (vertices[eij.first] + vertices[right_e.first]));
									usedEdges.insert(eij);
									usedEdges.insert(right_e);
									usedEdges.insert(brokenEdges.at(eij));
									usedEdges.insert(brokenEdges.at(right_e));
								}
							}
						}
					}
				}
			}

			//edge split and vertex split
			for (auto eit=lackEdges.begin(); eit!=lackEdges.end(); ++eit)
			{
				Edge eij = eit->first;
				if (usedEdges.find(eij)==usedEdges.end()) {
					Edge eji(eij.second, eij.first);
					if ((usedEdges.find(eji)==usedEdges.end()) && (lackEdges.find(eji)!=lackEdges.end())) {
						to_insert.push_back(0.5 * (vertices[eij.first] + vertices[eij.second]));  //edge split
						usedEdges.insert(eij);
						usedEdges.insert(eji);
						//usedEdges.insert(lackEdges.at(eij)); //since brokenedges will not be used, this is not necessary
						//usedEdges.insert(lackEdges.at(eji));
					}
					else { //eji is not available, try neighbors
						//left neighbor
						bool left_succeed = false;
						int i = eij.first;
						int j = eij.second;
						int idx = QEdges[i][j];
						//2 x ccw rotate
						idx = (idx+2) % 4;
						j = QD[i][idx];
						Edge le(i, j);
						if ((usedEdges.find(le)==usedEdges.end()) && (lackEdges.find(le)!=lackEdges.end())) {
							to_remove.insert(eij.first);
							to_insert.push_back(0.5 * (vertices[eij.first] + vertices[eij.second]));
							to_insert.push_back(0.5 * (vertices[le.first] + vertices[le.second]));
							usedEdges.insert(eij);
							usedEdges.insert(le);
							left_succeed = true;
						}
						//
						if (!left_succeed) { //try right neighbor
							i = eij.first;
							j = eij.second;
							idx = QEdges[i][j];
							//flip
							idx = OpEdges[i][j];
							gx_assert(idx!=-1);
							i = j;
							j = QD[i][idx];
							//2 x ccw rotate
							idx = (idx+2) % 4;
							j = QD[i][idx];
							//flip
							idx = OpEdges[i][j];
							if (idx!=-1) {
								i = j;
								j = QD[i][idx];
								Edge re(i, j);
								if ((usedEdges.find(re)==usedEdges.end()) && (lackEdges.find(re)!=lackEdges.end())) {
									to_remove.insert(eij.second);
									to_insert.push_back(0.5 * (vertices[eij.first] + vertices[eij.second]));
									to_insert.push_back(0.5 * (vertices[re.first] + vertices[re.second]));
									usedEdges.insert(eij);
									usedEdges.insert(re);
								}
							}
						}
					}
				}
			}

			//info
			std::cout << "Quad Dominant Topo Optimize:" << std::endl;
			std::cout << "Find " << brokenEdges.size() << " broken edges, " << lackEdges.size() << " lack edges" << std::endl;
			std::cout << to_remove.size() << " verts removed, " << to_insert.size() << " verts inserted." << std::endl;
		}


	//private:
		void init() {
			int nvert = vertices.size();
			//sort incident edges by direction angle
			std::vector<std::map<double, int>> QD0(nvert);  //vec<map<dir_angle, vidx>>
			for (auto it=edge_dir_match.begin(); it!=edge_dir_match.end(); ++it)
			{
				int dOrder = it->second.second;
				if (dOrder==0)
				{
					int vi = it->first.first;
					int vj = it->first.second;
					double dAngle = it->second.first;
					gx_assert(dAngle>=0 && dAngle<2*M_PI);
					QD0[vi].insert(std::make_pair(dAngle, vj));
				}
			}

			//build QD and QEdges
			QD.resize(nvert); 
			QEdges.resize(nvert);
			for (int i=0; i<nvert; ++i)
			{
				auto &map_ang_id = QD0[i];
				for (auto it=map_ang_id.begin(); it!=map_ang_id.end(); ++it) //it: <dir_angle, vidx>
				{
					QD[i].push_back(it->second);
					QEdges[i].insert(std::make_pair(QD[i].back(), QD[i].size()-1));

				}
				gx_assert(QD[i].size()==4); //4 edges per vertex
			}

			//build OpEdges
			OpEdges.resize(nvert);
			for (int i=0; i<nvert; ++i) {
				for (auto it=QD[i].begin(); it!=QD[i].end(); ++it) {
					int j = *it;
					OpEdges[i][j] = (QEdges[j].find(i)==QEdges[j].end()) ? -1 : QEdges[j][i];
				}
			}
		}


		const std::map<std::pair<int, int>, std::pair<double, int>> &edge_dir_match; 
		const std::vector<vec3> &vertices;
		//
		std::vector<std::vector<int>> QD;
		std::vector<std::map<int, int>> QEdges;
		std::vector<std::map<int, int>> OpEdges;
		//
		std::map<Edge, Edge> brokenEdges;
		std::map<Edge, Edge> lackEdges;
		bool broken_lack_edge_collected;
	};
	///dxy add end
}

#endif