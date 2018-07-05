#include "delaunay_cvt.h"
#include "cvt_geometry.h"
#include "cvt_compute.h"

#include <Geex/CVT/geometry.h>
#include <glut_viewer/glut_viewer.h>
#include <Geex/basics/stopwatch.h>
#include <Geex/combinatorics/map_utility.h>
#include <Geex/combinatorics/map_geometry.h>
#include <Geex/combinatorics/map_io.h>
#include <Geex/combinatorics/map_editor.h>
#include <Geex/combinatorics/map_builder.h>
#include <Geex/combinatorics/map_properties.h>
#include <Geex/combinatorics/map_curvature.h>
#include <Geex/CVT/AABB_tree_cgal.h>
#include <Geex/CVT/lloyd_energy.h>
#include <stack>
#include <queue>
#include <fstream>
//#include "rvd_ccvt.h"
//#include <Geex/basics/stopwatch.h>
//#include <fstream>
///dxy add
#include <time.h>
#include <sstream>

bool first = true;
int  curriter = 0 ;
int  call_nb_iter = 0;
int  currfundgraditer = 0 ;
//dxy add
long optimize_cvt_time = 0;
long optimize_dir_time = 0;
long optimize_interpolate_time = 0;
std::vector<double> optimize_cvt_energy;
std::vector<double> optimize_dir_energy;
std::vector<double> optimize_cvt_gnorm;
std::vector<double> optimize_dir_gnorm;
std::vector<double> optimize_tot_gnorm;




////////////////////////////////////////////////////////////////////////// 
void funcgrad_grvd(int N, double* x, double& f, double* g)
{
	//dxy add: info
	long cvt_time = 0;
	long dir_time = 0;
	long interpolate_time = 0;
	bool bad_Edir = (!optimize_dir_energy.empty()) && (optimize_dir_energy.back() != optimize_dir_energy.back());
	if (bad_Edir) std::cout << "Warning: bad direction energy!" << std::endl;
	// log
 	optimize_cvt_energy.push_back(0);
 	optimize_dir_energy.push_back(0);
 	optimize_cvt_gnorm.push_back(0);
 	optimize_dir_gnorm.push_back(0);
 	optimize_tot_gnorm.push_back(0);
	//
	Geex::DelaunayCVT* cvt = Geex::DelaunayCVT::instance() ;
	//dxy add
	if (cvt->lock_corners()) {
		cvt->snap_to_corners(N, x);
	}
	//
	if (cvt->lock_feature()) {
		cvt->snap_to_feature(N, x); //project feature sample points to its feature line, and update feature edge direction
	}
	//
	cvt->set_vertices(N, x) ;
	cvt->primal_dirty() = GL_TRUE ;
	f = 0 ;
	memset(g, 0, sizeof(double)*N) ;
	if(cvt->normal_aniso()>1.0) {
		cvt->RVD().for_each_triangle(Geex::ComputeLloydEnergyAndGradientAnisoN(
			cvt->RVD(), f, g, cvt->vertices(), cvt->normal_aniso())) ;
	}
	else {
		if (bad_Edir || (cvt->use_cvt_energy() && !cvt->use_direction())) {  //Pure CVT
			cvt_time = std::clock();
			if (cvt->use_multithread()) {
				cvt->multithread()->set_target(f, g);
				cvt->multithread()->run(&Geex::CVTThread::run_ComputeLloydEnergyAndGradient);
			}
			else {
				if (cvt->use_density()) {
					cvt->RVD().for_each_triangle(Geex::ComputeLloydEnergyAndGradientw(f, g, cvt->true_cvt_energy(), cvt->rvd_cell_area(), cvt->vertices())) ;
				}
				else {
					cvt->RVD().for_each_triangle(Geex::ComputeLloydEnergyAndGradient(f, g, cvt->true_cvt_energy(), cvt->rvd_cell_area(), cvt->vertices())) ;
				}
			}
			cvt_time = std::clock() - cvt_time;
 			optimize_cvt_energy.back() = f;
 			double gnorm2 = 0;
 			for(int i=0; i<N; ++i) {
 				gnorm2 += g[i]*g[i];
 			}
 			optimize_cvt_gnorm.back() = sqrt(gnorm2);
 			optimize_tot_gnorm.back() = optimize_cvt_gnorm.back();
			cvt->set_lloyd_grad(g);
		}
		else if (cvt->use_cvt_energy() && cvt->use_direction())  //FCVT
		{
			std::vector<unsigned int> Idx(cvt->nb_vertices());
			std::vector<std::map<unsigned int, std::vector<std::pair<Geex::TopoPolyVertexEdge, Geex::TopoPolyVertexEdge>>>> dual_segments(cvt->nb_vertices());
			//
			//calc CVT, Idx and dual_segments
			cvt->RVD().set_symbolic(true); //needed by dual segments computation
			cvt_time = std::clock();
			if (cvt->use_multithread()) {
				cvt->multithread()->set_target(f, g, Idx, dual_segments, cvt->true_cvt_energy(), cvt->rvd_cell_area());
				cvt->multithread()->run(&Geex::CVTThread::run_ComputeCVTPrepareFCVT); //if v.w is not constant, this becomes density fcvt
			}
			else {
				if (cvt->use_density()) cvt->RVD().for_each_facet(Geex::ComputeCVTPrepareFCVTw(f, g, Idx, dual_segments, cvt->true_cvt_energy(), cvt->rvd_cell_area(), cvt->vertices(), cvt->RVD()));
				else {	cvt->RVD().for_each_facet(Geex::ComputeCVTPrepareFCVT(f, g, Idx, dual_segments, cvt->true_cvt_energy(), cvt->rvd_cell_area(), cvt->vertices(), cvt->RVD())); }
			}
			cvt_time = std::clock() - cvt_time;
			cvt->RVD().set_symbolic(false);
			//
			optimize_cvt_energy.back() = f;
			double gnorm2 = 0;
			for(int i=0; i<N; ++i) {
				gnorm2 += g[i]*g[i];
			}
			optimize_cvt_gnorm.back() = sqrt(gnorm2);
			cvt->set_lloyd_grad(g);

			//Normal and Field Interpolation
			std::vector<Geex::vec3> Norm(cvt->nb_vertices(), Geex::vec3(0, 0, 0));
			std::vector<Geex::vec3> D(cvt->nb_vertices(), Geex::vec3(0, 0, 0));
			unsigned int RoSy = cvt->rot_symmetry();
			interpolate_time = std::clock();
			cvt->computeNormalAndDirection(Norm, D, Idx, RoSy);

			//Calc edge weights
			std::vector<std::map<unsigned int, double>> edge_weight;
			if (cvt->direction_edge_weight_mode() == Geex::Lloyd_Energy)
			{
				std::vector<double> v_w(cvt->nb_vertices(), 1.0);
				if (cvt->use_density()) {  //calc seed density
					cvt->computeWeight(v_w, dual_segments);
					//cvt->computeWeight(v_w, Idx);
					// 					for (int i=0; i<v_w.size(); ++i)
					// 					{
					// 						std::cout << "w " << i << " " << v_w[i] << std::endl;
					// 					}
				}
				cvt->compute_approximated_lloyd_energy(edge_weight, dual_segments, v_w); //calc edge_weight
				cvt->compute_approx_cvt_energy(edge_weight); //sum edge_weight to compute approximated cvt energy for each cell
				//test
				double approx_f = 0;
				for (int i=0; i<edge_weight.size(); ++i)
				{
					for (auto jt=edge_weight[i].begin(); jt!=edge_weight[i].end(); ++jt)
					{
						approx_f += jt->second;
					}
				}
				std::cout << "True Lloyd Energy: " << optimize_cvt_energy.back() << " Approx Energy: " << approx_f << std::endl;
				//
			}
			else {  //dual length edge weight
				cvt->compute_dual_length(edge_weight, dual_segments);
			}
			//
			interpolate_time = std::clock() - interpolate_time;

			//calc direction energy and gradient
			//double para_a = cvt->dir_para_a();
			dir_time = std::clock();
			//  			double  tEdir = 0;
			//  			double* tGdir = new double[N];
			//  			memset(tGdir, 0, sizeof(double)*N);
			//  			cvt->computeDirectionEnergyAndGradient(tEdir, tGdir, Norm, D, edge_weight, RoSy);
			if (cvt->use_edge_dir_match_topo_opt()) {
				cvt->compute_edge_dir_match(Norm, D, edge_weight, RoSy);  //only compute match, E and G are set to 0
			}
			else {
				if (cvt->use_edge_dir_match()) {
					cvt->computeDirectionEnergyAndGradient_match(Norm, D, edge_weight, RoSy);
				} else {
					cvt->computeDirectionEnergyAndGradient(Norm, D, edge_weight, RoSy); //results saved in cvt.direction_energy, cvt.direction_grad
				}
			}

			//mix factor
			double mfactor = 1;
			if (cvt->direction_edge_weight_mode() == Geex::Lloyd_Energy)
			{
				mfactor = cvt->energy_mixfactor() * 0.0075 * 30;
			}
			else {
				mfactor = cvt->energy_mixfactor() * 0.005 * cvt->average_area() * sqrt(cvt->average_area());
			}

			double Edir = std::accumulate(cvt->direction_energy().begin(), cvt->direction_energy().end(), 0.0);
			Edir *= mfactor;
			// 			tEdir *= mfactor;
			// 			std::cout << "Edir = " << Edir << " tEdir = " << tEdir << std::endl;

			auto &dir_g = cvt->direction_grad();
			double *Gdir = new double[N];
			for (int i=0; i<dir_g.size(); ++i)
			{
				dir_g[i] *= mfactor;
				for (int j=0; j<3; ++j) Gdir[3*i+j] = dir_g[i][j];		
			}
			//
			dir_time = std::clock() - dir_time;
			optimize_dir_energy.back() = Edir;
			gnorm2 = 0;
			for (int i=0; i<N; ++i)
			{
				gnorm2 += Gdir[i]*Gdir[i];
			}
			optimize_dir_gnorm.back() = sqrt(gnorm2);
			//cvt->set_direction_grad(Gdir);
			if (!cvt->use_stretch_topo_optimize()) {
				std::cout << "direction energy  = " << Edir << std::endl;
			}

			//add to f and g
			f += Edir;
			gnorm2 = 0;
			for (int i=0; i<N; ++i)
			{
				g[i] += Gdir[i];
				gnorm2 += g[i]*g[i];
			}
			optimize_tot_gnorm.back() = sqrt(gnorm2);
			delete[] Gdir;

			//visualize seed field
			cvt->set_seed_field(D);
			cvt->set_seed_normal(Norm);
		}		
	}

	//set grad of corners to 0
	if (cvt->lock_corners())
	{
		const auto& locked = cvt->locked();
		for (int i=0; i<locked.size(); ++i) {
			if (locked[i]) {
				g[3*i+0] = 0.0;
				g[3*i+1] = 0.0;
				g[3*i+2] = 0.0;
			}
		}
	}
	//align grad of feature points to feature edge direction
	if (cvt->lock_feature()) {
		cvt->align_feature_grad(N, g);
	}
	//

	call_nb_iter++;
	//dxy add: time
	optimize_cvt_time += cvt_time;
	optimize_dir_time += dir_time;
	optimize_interpolate_time += interpolate_time;
	//
}

//////////////////////////////////////////////////////////////////////////  
void newiteration_ccvt(int n, const double* x, double f, const double* g, double gnorm)
{
	Geex::DelaunayCVT* cvt = Geex::DelaunayCVT::instance() ;
	if (!cvt->use_stretch_topo_optimize())
	{
		std::cout.precision(16);
		std::cout << curriter <<": " << call_nb_iter <<" " << std::scientific << f <<" " << gnorm  << std::endl << std::endl ;
	}
	//	cvt->log_obtuse() ; 
	curriter ++;
	glut_viewer_redraw() ;
}

void newiteration_ccvt_save_pts(int n, const double* x, double f, const double* g, double gnorm)
{
	Geex::DelaunayCVT* cvt = Geex::DelaunayCVT::instance() ;
	if (!cvt->use_stretch_topo_optimize())
	{
		std::cout.precision(16);
		std::cout << curriter <<": " << call_nb_iter <<" " << std::scientific << f <<" " << gnorm  << std::endl << std::endl ;
	}
	//	cvt->log_obtuse() ; 
	curriter ++;
	//
	//std::string filename = "E:/Research/fcvt_svn/CVTRemesh/data/surfaces/fcvt/svn/simple/ellipsoid/180401-eg fastforward/eg18/autosave/" + std::to_string(std::clock()) + ".pts";
	std::string filename = "E:/Research/fcvt_svn/CVTRemesh/data/surfaces/fcvt/svn/simple/master-defense/autosave/" + std::to_string(std::clock()) + ".pts";
	cvt->save_denormalized_samples(filename);
}
//////////////////////////////////////////////////////////////////////////
namespace Geex {

	//dxy add
	GLenum DelaunayCVT::interpolationMode_ = Barycentric;
	//

	DelaunayCVT* DelaunayCVT::instance_ = nil ;

	DelaunayCVT::DelaunayCVT(const std::string& del_algo
		) : del_algo_(del_algo),
		delaunay_(Delaunay::create(del_algo)),
		boundary_(new TopoPolyMesh),
		RVD_(delaunay_, boundary_),
		exact_(GL_FALSE),
		symbolic_(GL_TRUE),
		use_density_(GL_FALSE), 
		primal_dirty_(GL_TRUE),
		tree_(nil),
		map_(new Map) {
			facet_index_.bind(map_, "fidx") ;
			vindex_.bind(map_, "vindex") ;
			//vgeodist_.bind(map_, "vgeodist") ;
			is_corner_.bind(map_, "vcorner") ;
			is_feature_.bind(map_, "hfeature") ;
			fchart_.bind(map_, "chart_id") ;
			lrvd_cur_step_ = 0 ;
			lfs_ = nil;
			affinity_weight_ = 0.0 ;
			cvt_normalization_ = 1.0 ;
			affinity_normalization_ = 1.0 ;
			//
			lock_corners_ = GL_FALSE/*GL_TRUE*/;
			//
			lock_feature_ = GL_FALSE ;
			use_feature_ = GL_FALSE ;
			angle_min_ = 30 ;
			angle_max_ = 90 ;
			nb_lock_rings_ = 3 ;
			symbolic_mode_backup_ = false;

			gx_assert(instance_ == nil) ;
			instance_ = this ;
			//		weighted_  = GL_FALSE ;
			constrained_ = GL_FALSE ;
			optimizer_mode_ = HLBFGS;
			opt_m_ = 7;
			normal_aniso_ = 1.0 ;

			///dxy add
			field_loaded_ = GL_FALSE;
			use_face_center_init_ = GL_FALSE;
			model_normalize_scale_ = -1; //uninitialized
			use_cvt_energy_ = GL_TRUE;
			use_direction_ = /*GL_TRUE*/ GL_FALSE;
			use_auto_save_ = GL_FALSE;
			use_true_gradient_ = GL_TRUE;
			use_special_quad_ = GL_FALSE;
			use_adaptive_mix_factor_ = GL_FALSE;
			energy_mixfactor_ = /*5.5*/ 0.35;
			direction_parameter_a_ = 1.5;
			direction_parameter_quad_ = 0.5;
			field_rot_symmetry_ = /*4*/ 6;
			use_facet_field_ = GL_FALSE;
			direction_energy_mode_ = Angle_Edge;
			direction_edge_weight_mode_ = Lloyd_Energy;
			use_new_interpolation_ = GL_FALSE;
			use_edge_dir_match_ = GL_FALSE;
			use_edge_dir_match_topo_opt_ = GL_FALSE;

			total_area_ = 0;
			delaunay_map_ = new Map;

			is_delaunay_map_dirty_ = true;
			use_stretch_topo_optimize_ = GL_FALSE;
			stretch_optimize_select_rate_ = 0.01;
			stretch_optimize_nb_smooth_iter_ = 5;
			map_topo_optimization_nb_iter_ = 0;
			map_topo_optimization_deltaR_ = 0;
			dihedral_angle_threshold_ = 180;
			experiment_iter_ = 10;

			//dxy add: multi-thread
			multithread_ = new CVTMultiThread();
			//
	}

	DelaunayCVT::~DelaunayCVT() { 
		instance_ = nil ; 
		facet_index_.unbind() ;
		vindex_.unbind() ;
		is_corner_.unbind() ;
		is_feature_.unbind() ;
		fchart_.unbind() ;
		delete tree_ ;
		delete boundary_ ;
		if(lfs_ != nil) delete lfs_ ;
		//dxy add
		delete delaunay_map_;
		delete multithread_;
		//
	}

	///dxy add: multi-thread
	void DelaunayCVT::init_multithread() {
		multithread_->initialize(instance_);
	}
	///

	// ---------------------------------------------------------------------
	// ----------------------- Lloyd CVT iteration -------------------------
	// ---------------------------------------------------------------------

	void DelaunayCVT::lloyd_global(int nb_iter) {
		unsigned int nv = nb_vertices() ;
		//		std::vector<vec3>& vertices = vertices() ;
		//		std::vector<int>& sample_facets = sample_facets() ;
		SystemStopwatch timer ;
		double time_dt = 0, time_vd = 0 ;

		for(int i=0; i<nb_iter; ++i) {
			std::vector<vec3>   mg_in(nv, vec3(0, 0, 0)) ;
			std::vector<double> m_in(nv, 0) ;

			//compute barycenter: global rvd
			compute_rvd() ;
			RVD().for_each_triangle(ComputeBarycenter(mg_in, m_in)) ;

			for(unsigned int j=0; j<nv; ++j) {
				gx_assert(m_in[j]>0) ;
				vec3 p = mg_in[j]/m_in[j] ;
				vec3 cp ;
				int  fid ;
				tree()->closest_point(p, cp, fid) ;
				vertices_[j] = cp ;
				//				sample_facets[j] = fid ;
			}
			//dxy add:
			is_delaunay_map_dirty_ = true;
			//
			glut_viewer_redraw() ;
		}
		//		std::cout << nb_iter << " global Lloyd iter: " << time_vd+time_dt << "s. dt time " <<time_dt <<", vd time " << time_vd << " s" <<std::endl ;
	}
	//----------------------------------------------------------------------------------------------------

	void DelaunayCVT::newton_lloyd_global(int nb_iter) {
		///dxy add
		if (lock_feature_ ) {
			update_seed_feature_line();
		}
		//

		unsigned int nv = nb_vertices() ;
		int n = nv * 3 ;
		int m = (int)opt_m_ ;
		if (optimizer_mode_ == LBFGSB && m == 0)
		{
			m = 1;
		}

		double *x = new double[n];

		//std::vector<vec3>& vertices = vertices() ;
		for (unsigned int i=0; i< nv; i++) {
			//DelaunayCVT::Vertex_handle it = all_vertices_[i] ;
			//x[3*i  ] = it->point().x();
			//x[3*i+1] = it->point().y();
			//x[3*i+2] = it->point().z();
			x[3*i  ] = vertices_[i][0] ;
			x[3*i+1] = vertices_[i][1] ;
			x[3*i+2] = vertices_[i][2] ;
		}

		curriter = 1;
		call_nb_iter = 0;

		double epsg = 1e-5, epsf=0, epsx=0;
		//dxy add: normalize epsg
		epsg *= nv;
		//

		//dxy test
		//epsg = 0.0;
		//

		Optimizer* opt = nil;
		///dxy add
		if (!use_stretch_topo_optimize_) {
			std::cerr << "Starting Newton (warming up...)" << std::endl ;
		}
		///
		switch(optimizer_mode_) {
		case LBFGSB:
			opt = new LBFGSBOptimizer();
			break;
		case HLBFGS:
			opt = new HLBFGSOptimizer() ;
			break ;
		case HM1QN3:
			opt = new HLBFGSOptimizer() ;
			static_cast<HLBFGSOptimizer*>(opt)->set_m1qn3(true);
			break ;
		case HCG:
			opt = new HLBFGSOptimizer() ;
			static_cast<HLBFGSOptimizer*>(opt)->set_cg(true);
			break ;
		default:
			gx_assert_not_reached ;
		}

		if (optimizer_mode_ == LBFGSB) {
			int *rhs = new int[n];
			memset(rhs, 0, sizeof(int)*n);
			double *lu = new double[n];
			memset(lu, 0, sizeof(double)*n);
			static_cast<LBFGSBOptimizer*>(opt)->set_nbd(n, rhs);
			static_cast<LBFGSBOptimizer*>(opt)->set_l(n, lu);
			static_cast<LBFGSBOptimizer*>(opt)->set_u(n, lu);
			delete[] rhs;
			delete[] lu;
		}

		opt->set_epsg(epsg) ;
		opt->set_epsf(epsf) ;
		opt->set_epsx(epsx) ;

		opt->set_M(m) ;
		opt->set_N(n) ;
		opt->set_max_iter(nb_iter) ;

		if (use_auto_save_) {
			opt->set_newiteration_callback(newiteration_ccvt_save_pts);
		}
		else {
			opt->set_newiteration_callback(newiteration_ccvt) ;
		}
		opt->set_funcgrad_callback(funcgrad_grvd) ;

		long tstart = std::clock();
		optimize_cvt_time = 0;
		optimize_dir_time = 0;
		optimize_interpolate_time = 0;
		optimize_cvt_energy.clear();
		optimize_dir_energy.clear();
		optimize_cvt_gnorm.clear();
		optimize_dir_gnorm.clear();
		optimize_tot_gnorm.clear();
		opt->optimize(x) ;
		long tend = std::clock();
		std::cout << "Outer Run Time(ms): " << tend - tstart << std::endl;
		std::cout << "Inner Run Time(ms): " << optimize_cvt_time+optimize_dir_time+optimize_interpolate_time
			<< " CVT: " << optimize_cvt_time << " Itp: " << optimize_interpolate_time << " Dir: " << optimize_dir_time << std::endl;
		///

		//set_vertices(x) ;
		if (lock_corners_) snap_to_corners(n, x);
		if (lock_feature_) snap_to_feature(n, x);
		set_vertices(n, x) ;
		delete opt;
		delete [] x;
	}

	void DelaunayCVT::load_map(const std::string& filename) {
		Geex::load_map(filename, map_) ;
		normalize_map(map_) ;
	}

	void DelaunayCVT::load_map(const std::string& filename, bool normalize) {
		Geex::load_map(filename, map_) ;
		if (normalize)
		{
			//normalize_map(map_) ;
			normalize_map_scale(map_, model_normalize_scale_);
		}
	}

	void DelaunayCVT::model_normalize_map() {
		normalize_map_scale(map_, model_normalize_scale_);
	}

	//dxy add
	double DelaunayCVT::calc_map_area() {
		double area = 0;
		FOR_EACH_FACET(Map, map_, fit) {
			area += fit->area();
		}
		return area;
	}
	//

	//dxy add
	bool DelaunayCVT::load_experiment(const std::string &filename) {
		std::ifstream in(filename.c_str());
		if (!in.is_open()) {
			std::cout << "cannot load experiment file " << filename << std::endl;
			return false;
		}

		experiments_.clear();

		std::stringstream ss;
		std::string line;

		//std::getline(in, line);
		//ss.clear();
		//ss.str(line);
		///
		bool succeed = true;
		while (std::getline(in, line))
		{
			ss.clear();
			ss.str(line);
			//
			Experiment e;
			std::string keyword;
			while (ss >> keyword) {
				if (keyword == "num") {
					int num;
					ss >> num;
					gx_assert(num>=0);
					e.num = num;
				}
				else if (keyword == "pts") {
					int pts;
					ss >> pts;
					e.pts = pts;
				}
				else if (keyword == "nb_iter") {
					int niter;
					ss >> niter;
					e.nb_iter = niter;
				}
				else if (keyword == "mix_factor") {
					double mf;
					ss >> mf;
					e.mix_factor = mf;
				}
				else if (keyword == "+density") {
					e.use_density = true;
				}
				else if (keyword == "+topo_opt") {
					e.use_topo_opt = true;
				}
				else if (keyword == "nb_topo_opt_iter") {
					int niter;
					ss >> niter;
					e.nb_topo_opt_iter = niter;
				}
				else if (keyword == "+quad") {
					e.use_quad = true;
				}
				else if (keyword == "+all_edge_split") {
					e.use_all_edge_split = true;
				}
				else {
					std::cout << "Error: unknown experiment keyword: " << keyword << std::endl;
					succeed = false;
					break;
				}
			}
			if (succeed)
			{
				experiments_.push_back(e);
			}else {
				break;
			}
		}
		in.close();
		if (succeed)
		{
			std::cout << experiments_.size() << " experiments loaded." << std::endl;
			return true;		
		}
		else {
			std::cout << "Error while loading exp file." << std::endl;
			experiments_.clear();
			return false;
		}
		///	
	}

	//
	void DelaunayCVT::save_feature(const std::string &filename) {
		std::ofstream feafile(filename.c_str());
		if (feafile.is_open()) {
			feafile << corners_.size() << " " << features_.size() << std::endl;
			for (auto it = corners_.begin(); it!=corners_.end(); ++it) {
				feafile << *it << " ";
			}
			feafile << std::endl;
			//
			for (int i=0; i<features_.size(); ++i) {
				const auto &line = features_[i];
				feafile << line.size() << std::endl;
				for (auto it=line.begin(); it!=line.end(); ++it) {
					feafile << *it << " ";
				}
				feafile << std::endl;
			}
			feafile.close();
		} else {
			std::cout << "Can't save fea file:" <<  filename << std::endl;
		}
	}

	bool DelaunayCVT::load_feature(const std::string& filename) {
		std::ifstream in(filename.c_str()) ;
		int nv, /*ne,*/ nl ;

		if( !in.is_open() ) {
			std::cout << "cannot load feature file " << filename << std::endl ;
			return false ;
		}	

		corners_.clear() ;
		features_.clear() ;

		in >> nv >> nl ;

		// load constrained vertex indices
		for(int i=0; i<nv; ++i) {
			int vidx ;
			in >> vidx ;
			corners_.push_back(vidx) ;
		}

		// load feature lines
		for(int i=0; i<nl; ++i) {
			int ni ;
			std::vector<int> line ;
			in >> ni ;
			for(int j=0; j<ni; ++j) {
				int vidx ;
				in >> vidx ;
				line.push_back(vidx) ;
			}
			features_.push_back(line) ;
		}
		in.close();
		//dxy add
		compute_feature_edge_line();
		//
		std::cout << filename << " loaded." << std::endl;
		return true ;
	}

	//dxy add


// 	void DelaunayCVT::mark_feature_on_delaunay_map() {
// 		FOR_EACH_VERTEX(Map, delaunay_map_, it) {
// 			//delaunay_is_corner_[it] = true;
// 			//delaunay_is_on_feature_[it] = true;
// 
// 			//delaunay_is_corner_[it] = false;
// 			//delaunay_is_on_feature_[it] = false;
// 			int vi = delaunay_map_vert_index_[it];
// 			if (locked_[vi]) {
// 				delaunay_is_corner_[it] = true;
// 				delaunay_is_on_feature_[it] = true;
// 			}
// 			auto lit = seed_feature_line_.find(vi);
// 			if (lit != seed_feature_line_.end()) {
// 				delaunay_is_on_feature_[it] = true;
// 				delaunay_feature_line_[it] = lit->second;
// 			}
// 		}
// 	}

	void DelaunayCVT::mark_feature() {
		for(int i=0; i<corners_.size(); ++i) {
			is_corner_[verts_[corners_[i]]] = true ;
		}

		for(int i=0; i<features_.size(); ++i) {
			std::vector<int>& line = features_[i] ;
			for(int j=0; j<line.size()-1; ++j) {
				int ls = line[j] ;
				int le = line[j+1] ;
				Map::Halfedge* h = verts_[ls]->halfedge() ;
				do {
					if(h->opposite()->vertex() == verts_[le]) {
						is_feature_[h] = true ;
						is_feature_[h->opposite()] = true ;
						break ;
					}
					h = h->next_around_vertex() ;
				} while(h!=verts_[ls]->halfedge()) ;
			}
		}

		FOR_EACH_HALFEDGE(Map, map_, h) {
			if(h->is_border() || h->opposite()->is_border()) {
				is_feature_[h] = true ;
			}
		}
	}

	void DelaunayCVT::extract_chart() {
		MapFacetProperty<bool> visited(map_, "visited") ;
		FOR_EACH_FACET(Map, map_, f) {
			visited[f] = false ;
		}

		charts_.clear() ;

		FOR_EACH_FACET(Map, map_, f) {
			if(visited[f]) continue ;
			std::vector<Map::Facet*> chart ;
			std::queue<Map::Facet*> Q ;
			Q.push(f) ;

			while(!Q.empty()) {
				Map::Facet* f = Q.front() ;
				Q.pop() ;
				if(visited[f]) continue ;

				chart.push_back(f) ;
				visited[f] = true ;
				Map::Halfedge* h = f->halfedge() ;
				do {
					Map::Facet* of = h->opposite()->facet() ;
					if(!is_feature_[h] && !visited[of]) {
						Q.push(of) ;
					}
					h = h->next() ;
				} while(h!=f->halfedge()) ;					
			}
			charts_.push_back(chart) ;
		}
		visited.unbind() ;
	}

	void DelaunayCVT::extract_features(Map* map, double angle) {
		MapHalfedgeProperty<bool> visited ;
		//		MapHalfedgeProperty<bool> is_feature ;
		//		MapVertexProperty<bool>   is_corner ;

		//		is_feature.bind(map, "is_feature") ;
		visited.bind(map, "visited") ;
		//		is_corner.bind(map, "is_corner") ;

		corners_.clear() ;
		features_.clear() ;  

		// mark sharp and boundary edges
		FOR_EACH_HALFEDGE(Map, map, h) {
			visited[h] = false ;
			is_feature_[h] = false ;

			if(h->is_border() || h->opposite()->is_border())
				is_feature_[h] = true ;
			else {
				vec3 n0 = Geom::facet_normal(h->facet()) ;
				vec3 n1 = Geom::facet_normal(h->opposite()->facet()) ;
				double d01 = dot(n0, n1) ;
				double a = acos(d01)*180/M_PI ;
				if(d01>=0 && a>=angle || d01<=0 && a<180-angle) {
					is_feature_[h] = true ;
				}
			}			 
		}

		// mark corner vertices
		FOR_EACH_VERTEX(Map, map, v) {
			std::vector<Map::Halfedge*> adjedge ;
			Map::Halfedge* h = v->halfedge() ;

			is_corner_[v] = false ;

			do {
				if(is_feature_[h]) {
					adjedge.push_back(h) ;
				}
				h = h->next_around_vertex() ;
			} while(h!=v->halfedge()) ;

			if(adjedge.size()==1) {
				is_corner_[v] = true ; // dart
				corners_.push_back(vindex_[v]) ;
			}
			else if(adjedge.size()==2) {
				vec3 v0 = v->point() ;
				vec3 v1 = adjedge[0]->opposite()->vertex()->point() ;
				vec3 v2 = adjedge[1]->opposite()->vertex()->point() ;
				double d01 = dot(normalize(v1-v0), normalize(v2-v0)) ;
				double a = acos(d01)*180/M_PI ;
				if(180-a>angle) {
					is_corner_[v] = true ; // crease
					corners_.push_back(vindex_[v]) ;
				}
			}
			else if(adjedge.size()>=3) {
				is_corner_[v] = true ; // corner
				corners_.push_back(vindex_[v]) ;
			}
		}

		std::cout << "number of corners " << corners_.size() << std::endl ;

		// extract feature lines
		FOR_EACH_VERTEX(Map, map, v) {
			// find out going halfedge
			if(is_corner_[v]) {
				Map::Halfedge* h = v->halfedge() ;
				do {
					h = h->next_around_vertex() ;
					if(is_feature_[h] && !visited[h]) {
						std::vector<int> line ;
						line.push_back(vindex_[v]) ;
						Map::Halfedge* h2 = h->opposite() ; 
						do {
							Map::Vertex* v2 = h2->vertex() ;
							visited[h2] = true ;
							visited[h2->opposite()] = true ;
							line.push_back(vindex_[v2]) ;
							if(!is_corner_[v2]) {
								do {
									h2 = h2->next_around_vertex() ;
								} while(!is_feature_[h2]) ;
								h2 = h2->opposite() ;
							} 
						} while(!is_corner_[h2->vertex()]) ;

						if(is_corner_[h2->vertex()]) {
							visited[h2] = true ;
							visited[h2->opposite()] = true ;
							line.push_back(vindex_[h2->vertex()]) ;
							features_.push_back(line) ;
						}
					}
				} while(h!=v->halfedge()) ;
			}
		}

		// extract closed feature loops
		FOR_EACH_HALFEDGE(Map, map, h) {
			if(is_feature_[h] && !visited[h]) {
				Map::Halfedge* hcur = h ;
				std::vector<int> loop ;
				do {
					Map::Halfedge* hnext = hcur ;

					loop.push_back(vindex_[hcur->vertex()]) ;
					visited[hcur] = true ;
					visited[hcur->opposite()] = true ;

					do {
						hnext=hnext->next_around_vertex() ;
					} while(hnext!=hcur && !is_feature_[hnext]) ;

					if(is_feature_[hnext]) {
						hcur = hnext->opposite() ;
					}
				} while(hcur!=h) ;
				loop.push_back(vindex_[hcur->vertex()]) ; // duplicate the starting vertex
				features_.push_back(loop) ;
			}
		}

		std::cout << "number of feature lines " << features_.size() << std::endl ;
		//
		//		is_corner.unbind() ;
		//		is_feature.unbind() ;
		visited.unbind() ;

		//dxy add
		compute_feature_edge_line();
	}

	void DelaunayCVT::update_indices() {
		facets_.clear() ;
		verts_.clear() ;

		int fidx = 0 ;
		FOR_EACH_FACET(Map, map_, f) {
			facet_index_[f] = fidx ;
			facets_.push_back(f) ;			
			fidx ++ ;
		}

		int vidx = 0 ;
		FOR_EACH_VERTEX(Map, map_, v) {
			vindex_[v] = vidx ;
			verts_.push_back(v) ;
			vidx ++ ;
		}
	}

	void DelaunayCVT::load_boundary(const std::string& filename) {
		unsigned int nb_borders = boundary_->load(filename) ;
		double area = 0 ;
		//        use_v_map = (nb_borders > 0) ;
		bbox_dirty_ = true ;
		boundary_facet_plane_.resize(boundary_->nb_facets()) ;
		for(unsigned int f=0; f<boundary_->nb_facets(); f++) {
			boundary_facet_plane_[f] = boundary_->facet_plane(f) ;
		}
		for(unsigned int i=0; i<boundary_->nb_facets(); i++) {
			if(length(boundary_->facet_normal(i)) < 1e-10) {
				std::cerr << " ============ DEGENERATE FACETS !! ========= " << std::endl ;
			}
			area += boundary_->facet_area(i) ;
		}
		std::cout.precision(20) ;
		std::cout << "total area = " << area << std::endl ;

		if(bbox_dirty_) {
			boundary_->get_bbox(x_min_, y_min_, z_min_, x_max_, y_max_, z_max_) ;
			bbox_dirty_ = false ;
		}

		// normalize the area of boundary mesh
		// double scale = sqrt(1.0/area) ;
		//normalize boundary mesh into a bounding box with longest axis equal to 1.0
		double scale = 1.0 / gx_max(x_max_-x_min_, gx_max(y_max_-y_min_, z_max_-z_min_)) ;
		for(unsigned int i=0; i<boundary_->nb_vertices(); ++i) {
			vec3 newp = scale*(boundary_->vertex(i)-vec3(x_min_, y_min_, z_min_)) ;
			boundary_->vertex(i).set_point(newp) ;
		}

		// we have to normalize the original vertices too
		std::vector<vec3>& ovs = boundary_->original_vertices() ;
		for(unsigned int i=0; i<ovs.size(); ++i) {
			ovs[i] = scale*(ovs[i] - vec3(x_min_, y_min_, z_min_)) ;
		}

		area = 0 ;
		for(unsigned int i=0; i<boundary_->nb_facets(); i++) {
			area += boundary_->facet_area(i) ;
		}
		std::cout << "total area = " << area << std::endl ;

		compute_lfs() ;

		//		// update weight of the boundary
		//		update_weight() ;
		update_indices() ;

		// build AABB tree
		tree_ = new AABBTree(boundary_) ;
	}

	void DelaunayCVT::load_boundary(const std::string& filename, bool normalize) {
		unsigned int nb_borders = boundary_->load(filename) ;
		double area = 0 ;
		//        use_v_map = (nb_borders > 0) ;
		bbox_dirty_ = true ;
		boundary_facet_plane_.resize(boundary_->nb_facets()) ;
		for(unsigned int f=0; f<boundary_->nb_facets(); f++) {
			boundary_facet_plane_[f] = boundary_->facet_plane(f) ;
		}
		for(unsigned int i=0; i<boundary_->nb_facets(); i++) {
			if(length(boundary_->facet_normal(i)) < 1e-10) {
				std::cerr << " ============ DEGENERATE FACETS !! ========= " << std::endl ;
			}
			area += boundary_->facet_area(i) ;
		}
		std::cout.precision(20) ;
		std::cout << "total area = " << area << std::endl ;

		if(bbox_dirty_) {
			boundary_->get_bbox(x_min_, y_min_, z_min_, x_max_, y_max_, z_max_) ;
			bbox_dirty_ = false ;
		}

		if (normalize)
		{
			// normalize the area of boundary mesh
			//double scale = sqrt(1.0/area) ;
			// normalize the area per seed
			//double scale = sqrt(nb_points/area);
			//normalize boundary mesh into a bounding box with longest axis equal to 1.0
			//double scale = 1.0 / gx_max(x_max_-x_min_, gx_max(y_max_-y_min_, z_max_-z_min_)) ;
			double scale = model_normalize_scale_;

			for(unsigned int i=0; i<boundary_->nb_vertices(); ++i) {
				vec3 newp = scale*(boundary_->vertex(i)-vec3(x_min_, y_min_, z_min_)) ;
				boundary_->vertex(i).set_point(newp) ;
			}

			// we have to normalize the original vertices too
			std::vector<vec3>& ovs = boundary_->original_vertices() ;
			for(unsigned int i=0; i<ovs.size(); ++i) {
				ovs[i] = scale*(ovs[i] - vec3(x_min_, y_min_, z_min_)) ;
			}
		}

		area = 0 ;
		for(unsigned int i=0; i<boundary_->nb_facets(); i++) {
			area += boundary_->facet_area(i) ;
		}
		std::cout << "total area = " << area << std::endl ;

		//compute_lfs() ; //dxy comment

		//		// update weight of the boundary
		//		update_weight() ;
		update_indices() ;

		// build AABB tree
		tree_ = new AABBTree(boundary_) ;
	}

	void DelaunayCVT::save_samples(const std::string& filename) {
		std::ofstream out(filename.c_str()) ;

		out << vertices_.size() << std::endl ;
		out.precision(20) ;
		for(unsigned int i=0; i<vertices_.size(); ++i) {
			out << vertices_[i] /*<< " " << sample_facets_[i]*/ << std::endl ; ///dxy comment: sample_facets is not updated, thus invalid
		}

		out.close() ;
	}

	///dxy add
	void DelaunayCVT::save_denormalized_samples(const std::string& filename) {
		double scale = 1.0 / model_normalize_scale_;
		std::ofstream out(filename.c_str()) ;
		out << vertices_.size() << std::endl;
		out.precision(20);
		for(unsigned int i=0; i<vertices_.size(); i++) {
			vec3 newp = scale*(vertices_[i]) + vec3(x_min_, y_min_, z_min_) ;
			out << newp << std::endl ;
		}
		out.close();
	}

	void DelaunayCVT::load_samples(const std::string& filename, bool normalize) {
		load_samples(filename);
		if (normalize)
		{
			normalize_samples();
		}
	}

	void DelaunayCVT::normalize_samples() {
		for (int i=0; i<vertices_.size(); ++i)
		{
			vertices_[i] = model_normalize_scale_ * (vertices_[i] - vec3(x_min_, y_min_, z_min_));
		}
	}
	///

	void DelaunayCVT::load_samples(const std::string& filename) {
		std::ifstream in(filename.c_str()) ;
		unsigned int nb = 0 ;

		in >> nb ;
		vertices_.clear() ;
		vertices_.resize(nb) ;
		locked_.clear();
		locked_.resize(nb, false) ;
		//dxy add
		on_feature_.clear();
		on_feature_.resize(nb, false);
		//
		vertices_corner_index_.clear();
		vertices_corner_index_.resize(nb, -1);
		//
		// 		sample_facets_.clear() ;
		// 		sample_facets_.resize(nb) ;

		for(unsigned int i=0; i<nb; ++i) {
			vec3 p ;
			//int  f ;
			in >> p /*>> f*/ ;
			vertices_[i] = p ;
			//sample_facets_[i] = f ;
		}
		// hardcode for reading points of CapCVT
		//for(unsigned int i=0; i<nb; ++i) {
		//	vec3 p ;
		//	int  f ;
		//	//in >> p >> f ;
		//	in >> p ;
		//	vertices_[i] = p ;
		//	vec3 fp ;
		//	tree_->closest_point(p, fp, f) ;
		//	sample_facets_[i] = f ;
		//}

		//// hardcode for reading points of BNOT
		//for(unsigned int i=0; i<nb; ++i) {
		//	vec3 p ;
		//	int  f ;
		//	double w ;
		//	//in >> p >> f ;
		//	in >> p >> w ;
		//	vertices_[i] = p ;
		//	vec3 fp ;
		//	tree_->closest_point(p, fp, f) ;
		//	sample_facets_[i] = f ;
		//}
		in.close() ;

		if(corners_.size()>0) {
			for(unsigned int i=0; i<nb; ++i) {
				for(unsigned int j=0; j<corners_.size(); ++j) {
					if(same_point(vertices_[i], verts_[corners_[j]]->point())) {
						locked_[i] = true ;
						//
						on_feature_[i] = true;
						vertices_corner_index_[i] = j;
					}
				}
			}
		}
		//dxy add:
		is_delaunay_map_dirty_ = true;
	}

	void DelaunayCVT::load_density(const std::string& filename) {
		std::ifstream in(filename.c_str()) ;
		MapVertexProperty<double> density(map_, "vertex_density") ;
		int nv ;
		in >> nv ;
		if(nv!=map_->nb_vertices()) {
			std::cerr << "nb vertices doesn't match..." << std::endl ;
			return ;
		}

		FOR_EACH_VERTEX(Map, map_, v) {
			double w ;
			in >> w ; 
			density[v] = w ;
		}

		in.close() ;
	}

	void DelaunayCVT::save_density(const std::string& filename) {
		MapVertexProperty<double> density ;
		density.bind_if_defined(map_, "vertex_density") ;
		if(!density.is_bound()) {
			std::cerr << "density is not defined..." << std::endl ;
			return ;
		}
		std::ofstream out(filename.c_str()) ;

		out << map_->nb_vertices() << std::endl ;

		FOR_EACH_VERTEX(Map, map_, v) {
			out << density[v] << std::endl ; ;
		}

		out.close() ;
	}

	void DelaunayCVT::get_bbox(
		real& x_min, real& y_min, real& z_min,
		real& x_max, real& y_max, real& z_max
		)  {
			Geex::get_bbox(map_, x_min, y_min, z_min, x_max, y_max, z_max) ;
	}

	///dxy add:
	//init samples at geometric centers
	void DelaunayCVT::init_samples_barycentric(int& nb_pts)
	{
		nb_pts = 0 ;

		// initial sampling
		vertices_.clear() ;
		locked_.clear() ;
		on_feature_.clear();
		vertices_corner_index_.clear();

		for (int i=0; i<verts_.size(); ++i) {
			vertices_.push_back(verts_[i]->point());
			locked_.push_back(false);
			on_feature_.push_back(false);
			vertices_corner_index_.push_back(-1);
		}

		if(corners_.size()>0) {
			for(unsigned int i=0; i<corners_.size(); ++i) {
				locked_[corners_[i]] = true;
				on_feature_[corners_[i]] = true;
				vertices_corner_index_[corners_[i]] = i;
			}
		}

		for(unsigned int i=0; i<boundary_->nb_facets(); ++i) {
			int v0 = boundary_->facet_begin(i) ;
			vertices_.push_back((boundary_->vertex(v0)+boundary_->vertex(v0+1)+boundary_->vertex(v0+2))/3.0);
			locked_.push_back(false);
			on_feature_.push_back(false);
			vertices_corner_index_.push_back(-1);
		}

		//vertices_.insert(vertices_.end(), boundary_->original_vertices().begin(), boundary_->original_vertices().end());

		nb_pts = vertices_.size();

		//dxy add:
		is_delaunay_map_dirty_ = true;
	}
	///

	void DelaunayCVT::init_samples(int nb_pts) {
		// split map for init sampling
		//while(map_->nb_facets()<10*nb_pts) {
		//	MapEditor editor(map_) ;
		//	editor.split_surface(split_cloop) ;
		//}
		//gx_assert(boundary_->nb_facets()>10*nb_pts) ;

		//if(tree_ != nil) delete tree_ ;
		//tree_ = new AABBTree(map_) ;

		int nb_cur = 0 ;
		MapFacetProperty<bool>   valid ;

		valid.bind(map_, "valid") ;
		FOR_EACH_FACET(Map, map_, f) {
			valid[f] = true ;
		}

		// initial sampling
		//samples_.clear() ;
		vertices_.clear() ;
		sample_facets_.clear() ;
		locked_.clear() ;
		on_feature_.clear();
		vertices_corner_index_.clear();

		if(corners_.size()>0) {
			for(unsigned int i=0; i<corners_.size(); ++i) {
				vertices_.push_back(verts_[corners_[i]]->point());
				sample_facets_.push_back(facet_index_[verts_[corners_[i]]->halfedge()->facet()]) ;
				locked_.push_back(true) ;
				on_feature_.push_back(true);
				vertices_corner_index_.push_back(i);
			}
			nb_cur = corners_.size() ;
		}

		do {
			int fidx = rand()%map_->nb_facets() ;
			//	if(!valid[facets_[fidx]]) continue ;
			Map::Facet* f = facets_[fidx] ;
			vec3 v0 = f->halfedge()->vertex()->point() ;
			vec3 v1 = f->halfedge()->next()->vertex()->point() ;
			vec3 v2 = f->halfedge()->next()->next()->vertex()->point() ;
			vec3 p = random_point_tri(v0, v1, v2) ;
			//	samples_.push_back(SamplePoint(p, 0, fidx)) ;
			vertices_.push_back(p) ;
			sample_facets_.push_back(fidx) ;
			valid[f] = false ;
			locked_.push_back(false) ;
			on_feature_.push_back(false);
			vertices_corner_index_.push_back(-1);
			nb_cur++ ;
		} while(nb_cur != nb_pts) ;

		valid.unbind() ;

		//dxy add:
		is_delaunay_map_dirty_ = true;
	}

	void DelaunayCVT::init_samples_w(int nb) {
		double totalw = 0. ; 
		double res = 0 ;
		int    check = 0  ;

		std::vector<double> face_mass(boundary_->nb_facets());
		for(unsigned int i=0; i<boundary_->nb_facets(); ++i) {
			//const Geex::Facet& F = boundary_[i] ;
			//face_mass[i] = Geex::tri_mass(F.vertex[0], F.vertex[1], F.vertex[2], F.vertex_weight[0], F.vertex_weight[1], F.vertex_weight[2]) ;
			//totalw += face_mass[i];
			int v0 = boundary_->facet_begin(i) ;
			if(use_density_)
				face_mass[i] = tri_mass(boundary_->vertex(v0), boundary_->vertex(v0+1), boundary_->vertex(v0+2), boundary_->vertex(v0).w, boundary_->vertex(v0+1).w, boundary_->vertex(v0+2).w) ;
			else
				face_mass[i] = tri_area(boundary_->vertex(v0), boundary_->vertex(v0+1), boundary_->vertex(v0+2)) ;
			totalw += face_mass[i] ;
		}

		vertices_.clear() ;
		locked_.clear() ;
		on_feature_.clear();
		vertices_corner_index_.clear();
		//		sample_facets_.clear() ;

		if(corners_.size()>0) {
			for(unsigned int i=0; i<corners_.size(); ++i) {
				vertices_.push_back(verts_[corners_[i]]->point()) ;
				//				sample_facets_.push_back(facet_index_[verts_[corners_[i]]->halfedge()->facet()]) ;
				locked_.push_back(true) ;
				on_feature_.push_back(true);
				vertices_corner_index_.push_back(i);
			}
		}

		for(unsigned int i=0; i<boundary_->nb_facets(); ++i) {
			//const Geex::Facet& F = boundary_[i] ;
			double nf = (nb-corners_.size())*face_mass[i]/totalw ;
			int n = (int)(nf+res) ;
			if(n >= 1) {
				res = nf+res - n ;
				check += n ;
				//for(int j=0; j<n; ++j) {
				//	double r0 = Numeric::random_float64() ;
				//	double r1 = (1-r0)*Numeric::random_float64() ;
				//	pts.push_back(r0*F.vertex[0]+r1*F.vertex[1]+(1-r0-r1)*F.vertex[2]);
				//	//insert(r0*F.vertex[0]+r1*F.vertex[1]+(1-r0-r1)*F.vertex[2]) ;			
				//}
				int v0 = boundary_->facet_begin(i) ;
				for(int j=0; j<n; ++j) {
					vec3 p = random_point_tri(boundary_->vertex(v0), boundary_->vertex(v0+1), boundary_->vertex(v0+2)) ;
					vertices_.push_back(p) ;
					//					sample_facets_.push_back(i) ;
					locked_.push_back(false) ;
					on_feature_.push_back(false);
					vertices_corner_index_.push_back(-1);
				}
			}
			else
				res += nf ;
		}
		//dxy add
		is_delaunay_map_dirty_ = true;
	}

	///dxy add:
	void DelaunayCVT::update_area() {
		if (total_area_ == 0)
		{
			total_area_ = boundary_->area();
		}
		average_area_ = total_area_ / nb_vertices();

		std::cout << "total area " << total_area_ << std::endl;
		std::cout << "average area " << average_area_ << std::endl;
		std::cout << "nb vertices " << nb_vertices() << std::endl;
	}


	void DelaunayCVT::split_all_edges() {
		update_delaunay_map();
		FOR_EACH_HALFEDGE(Map, delaunay_map_, eit) {
			if (eit->is_edge_key())
			{
				Map::Vertex* v1 = eit->prev()->vertex();
				Map::Vertex* v2 = eit->vertex();
				vec3& p1 = v1->point();
				vec3& p2 = v2->point();
				vertices_.push_back(0.5 * (p1+p2));
				locked_.push_back(false);
				on_feature_.push_back(false);
			}	
		}
		//smoothing
		lloyd_global(1);
		is_delaunay_map_dirty_ = true;
	}

	void DelaunayCVT::split_long_edges() {
		std::vector<std::pair<unsigned, unsigned>> edges;
		RVD_.for_each_primal_triangle(
			CollectLongEdges(edges, valences_)
			);
		for (int i=0; i<edges.size(); ++i)
		{
			unsigned v1 = edges[i].first;
			unsigned v2 = edges[i].second;
			vec3 p = (vertices_[v1] + vertices_[v2]) * 0.5;  //middle point
			vertices_.push_back(p);
			locked_.push_back(false);
			on_feature_.push_back(false);
		}
		//test: smoothing
		lloyd_global(5);
		//
		std::cout << edges.size() << " vertices added." << std::endl;

		is_delaunay_map_dirty_ = true;
	}

	void DelaunayCVT::collapse_short_edges() {
		std::vector<std::pair<unsigned, unsigned>> edges;
		RVD_.for_each_primal_triangle(
			CollectShortEdges(edges, valences_)
			);
		std::vector<bool> to_remove(vertices_.size(), false);
		std::vector<vec3> new_vertices;
		std::vector<bool> new_locked;
		std::vector<bool> new_on_feature;
		for (int i=0; i<edges.size(); ++i)
		{
			unsigned v1 = edges[i].first;
			unsigned v2 = edges[i].second;
			to_remove[v1] = true;
			to_remove[v2] = true;
			vec3 p = (vertices_[v1] + vertices_[v2]) * 0.5;  //middle point
			new_vertices.push_back(p);
			new_locked.push_back(false);
			new_on_feature.push_back(false);
		}
		for (int i=0; i<to_remove.size(); ++i)
		{
			if (!to_remove[i])
			{
				new_vertices.push_back(vertices_[i]);
				new_locked.push_back(locked_[i]);
				new_on_feature.push_back(on_feature_[i]);
			}
		}
		vertices_ = new_vertices;
		locked_ = new_locked;
		on_feature_ = new_on_feature;
		//test: smoothing
		lloyd_global(5);
		//
		std::cout << edges.size() << " vertices removed." << std::endl;

		is_delaunay_map_dirty_ = true;
	}

	void DelaunayCVT::update_vertices_from_delaunay_map() {
		vertices_.clear();
		std::vector<bool> new_locked, new_on_feature;

		FOR_EACH_VERTEX(Map, delaunay_map_, vit) {
			vertices_.push_back(vit->point());
			int vid = (delaunay_map_vert_index_.find(vit)!=delaunay_map_vert_index_.end()) ? delaunay_map_vert_index_[vit] : -1;
			if (vid > -1) {
				new_locked.push_back(locked_[vid]);
				new_on_feature.push_back(on_feature_[vid]);
			}
			else {
				new_locked.push_back(false);
				new_on_feature.push_back(false);
			}
		}
		locked_ = new_locked;
		on_feature_ = new_on_feature;
		//
		compute_rvd();
		//
		is_delaunay_map_dirty_ = true;
		//
		glut_viewer_redraw();

	}

	MapBuilder* BuildMapFacet::builder = 0;  //static var init

	void DelaunayCVT::update_delaunay_map()
	{
		if (!is_delaunay_map_dirty_) return;

		delaunay_map_->clear();
		delaunay_map_vert_index_.clear();
		MapBuilder builder(delaunay_map_);
		builder.begin_surface();
		//vert
		for (int i=0; i<vertices_.size(); ++i)
		{
			builder.add_vertex(vertices_[i]);
			delaunay_map_vert_index_.insert(std::make_pair(builder.vertex(i), i));
		}
		//face
		BuildMapFacet::builder = &builder;
		RVD_.for_each_primal_triangle(BuildMapFacet(/*&builder*/));

		builder.end_surface();
		is_delaunay_map_dirty_ = false;
		is_map_edge_stretch_dirty_ = true;
		is_primal_angles_dirty_ = true;
		is_primal_areas_dirty_ = true;
	}

	void DelaunayCVT::update_map_edge_stretch() {
		update_delaunay_map();
		if (!is_map_edge_stretch_dirty_) return;

		delaunay_map_edge_stretch_.clear();
		FOR_EACH_HALFEDGE(Map, delaunay_map_, eit) {
			if (eit->is_edge_key())
			{
				Map::Vertex* v1 = eit->prev()->vertex();
				Map::Vertex* v2 = eit->vertex();
				vec3& p1 = v1->point();
				vec3& p2 = v2->point();
				vec3 e12 = p2 - p1;
				e12 = normalize(e12);  //test..
				vec3& g1 = direction_grad_[delaunay_map_vert_index_[v1]];
				vec3& g2 = direction_grad_[delaunay_map_vert_index_[v2]];
				double stretch = dot(g1-g2, e12);
				delaunay_map_edge_stretch_.insert(std::make_pair(eit, stretch));
			}
		}
		is_map_edge_stretch_dirty_ = false;
	}

	void DelaunayCVT::calc_singularities(int& nb) {
		nb = 0;
		update_delaunay_map();
		FOR_EACH_VERTEX(Map, delaunay_map_, vit) {
			if (vit->degree() != 6)
			{
				nb++;
			}
		}
	}


	bool DelaunayCVT::find_topo_operation_preserve_fea(Map::Vertex_iterator v, TopoOperation& op) {
		int dv = v->degree();  //TODO: virtual degree
		if (dv == 6 || v->is_on_border()) return false;

		Map::Halfedge* h = v->halfedge();
		op = TopoOperation(v, h, Null_Operator, 0);

		/* feature preservation assumption: 
		*delaunay_map_vert_index_[v] return the right id in vertices_, if v is in delaunay_map_vert_index_
		* in other words, all delaunay map feature Vertex are preserved during map_topo_optimize
		*/

		// feature check for v
		int vid = (delaunay_map_vert_index_.find(v)!=delaunay_map_vert_index_.end()) ? delaunay_map_vert_index_[v] : -1;
		bool v_is_corner = false;
		bool v_is_fea = false;
		int  v_fea_line = -1;
		std::set<int> v_corner_lines;
		if (vid > -1) {
			v_is_corner = locked_[vid];
			v_is_fea = v_is_corner || on_feature_[vid];
			if (v_is_corner) {
				v_corner_lines = seed_feature_lines_[vid];
			}
			else if (v_is_fea) {
				v_fea_line = seed_feature_line_[vid];
			}
		}
		//

		std::vector<int> deg;
		deg.reserve(dv);
		do {
			deg.push_back(h->prev()->vertex()->degree()); //TODO: virtual degree
			h = h->next_around_vertex();
		} while (h != v->halfedge());

		std::vector<Map::Halfedge*> drift_edges;
		drift_edges.reserve(dv);

		int delRmin = 100; //test
		int d2, d3, d4, dk;
		int shift = (dv%2) ? ((dv-1)/2) : (dv/2);
		int xy = (dv%2) ? ((dv-1)*(dv-3)/4) : ((dv/2 -1)*(dv/2 - 1));
		h = v->halfedge();
		for(int i=0; i<dv; ++i) {
			d2 = deg[i];
			d3 = deg[(i+dv-1)%dv];
			d4 = deg[(i+1)%dv];
			dk = deg[(i+shift)%dv];
			// feature check
			Map::Vertex *vh = h->prev()->vertex();
			int vhid = (delaunay_map_vert_index_.find(vh) != delaunay_map_vert_index_.end()) ? delaunay_map_vert_index_[vh] : -1;
			bool vh_is_corner = false;
			bool vh_is_fea = false;
			std::set<int> vh_corner_lines;
			int vh_fea_line = -1;
			if (vhid > -1) {
				vh_is_corner = locked_[vhid];
				vh_is_fea = vh_is_corner || on_feature_[vhid];
				if (vh_is_corner) vh_corner_lines = seed_feature_lines_[vhid];
				else if (vh_is_fea) vh_fea_line = seed_feature_line_[vhid];
			}
			bool h_is_fea = false;
			if (v_is_corner && vh_is_corner) {
				for (auto lit=v_corner_lines.begin(); lit!=v_corner_lines.end(); ++lit) {
					if (vh_corner_lines.find(*lit) != vh_corner_lines.end()) {
						h_is_fea = true;
						break;
					}
				}
			}
			else if (v_is_corner && vh_is_fea) {
				if (v_corner_lines.find(vh_fea_line) != v_corner_lines.end()) {
					h_is_fea = true;
				}
			}
			else if (vh_is_corner && v_is_fea) {
				if (vh_corner_lines.find(v_fea_line) != vh_corner_lines.end()) {
					h_is_fea = true;
				}
			}
			else if (v_is_fea && vh_is_fea) {
				h_is_fea = (v_fea_line == vh_fea_line);
			}

			// Edge Flip
			if (!h_is_fea) {
				int delR_EF = 4 + 2*(d3 + d4 - dv - d2);
				if (delR_EF < delRmin) {
					delRmin = delR_EF;
					op._e = h;
					op._type = Edge_Flip;
					op._delta_R = delR_EF;
				}
				//test: record drift edges
				if (delR_EF == 0)
				{
					drift_edges.push_back(h);
				}
			}

			// Vertex Split
			if (!v_is_fea) {
				int delR_VS = 2*(dv - xy + d2 + dk - 12);
				if (delR_VS < delRmin)
				{
					delRmin = delR_VS;
					op._e = h;
					op._type = Vertex_Split;
					op._delta_R = delR_VS;
				}
			}

			// Edge Collapse
			if ((!v_is_fea) && (!vh_is_fea)) {
				int delR_EC = 4 + 2*(dv-5)*(d2-5) + 2*(dv+d2-d3-d4);
				if (delR_EC < delRmin)
				{
					delRmin = delR_EC;
					op._e = h;
					op._type = Edge_Collapse;
					op._delta_R = delR_EC;
				}
			}
			h = h->next_around_vertex();
		}

		//test: drift edge is always acceptable
		if (op._delta_R >0 || (op._delta_R == 0 && op._type != Edge_Flip))
		{
			if (!drift_edges.empty())
			{
				op._e = drift_edges[0];
				op._type = Edge_Flip;
				op._delta_R = 0;
			}
			else {
				op._type = Null_Operator;
			}
		}
		//

		if (op._type == Null_Operator) return false;
		else {  
// 			//crease check
// 			Map::Vertex *v2 = op._e->opposite()->vertex();
// 			int cnt = 0;
// 			switch (op._type)
// 			{
// 			case Edge_Flip:   
// 				if (op._e->is_sharp(dihedral_dot)) return false;
// 				break;
// 			case Edge_Collapse:  //
// 				if (op._v->is_corner(dihedral_dot) || v2->is_corner(dihedral_dot) ||
// 					(!op._e->is_sharp(dihedral_dot) && (op._v->is_sharp(dihedral_dot) || v2->is_sharp(dihedral_dot))))
// 				{
// 					return false;
// 				}
// 				break;
// 			case Vertex_Split:
// 				if (op._v->is_sharp(dihedral_dot)) return false;
// 				break;
// 			default:
// 				std::cout << "Error: operation type not defined!" << std::endl;
// 				return false;
// 			}
			return true;
		}
	}

	bool DelaunayCVT::find_topo_operation(Map::Vertex_iterator v, TopoOperation& op, double dihedral_dot /*= -1.0*/) {
		int dv = v->degree();  //TODO: virtual degree
		if (dv == 6 || v->is_on_border()) return false;

		Map::Halfedge* h = v->halfedge();
		op = TopoOperation(v, h, Null_Operator, 0);


		std::vector<int> deg;
		deg.reserve(dv);
		do {
			deg.push_back(h->prev()->vertex()->degree()); //TODO: virtual degree
			h = h->next_around_vertex();
		} while (h != v->halfedge());

		std::vector<Map::Halfedge*> drift_edges;
		drift_edges.reserve(dv);

		int delRmin = 100; //test
		int d2, d3, d4, dk;
		int shift = (dv%2) ? ((dv-1)/2) : (dv/2);
		int xy = (dv%2) ? ((dv-1)*(dv-3)/4) : ((dv/2 -1)*(dv/2 - 1));
		h = v->halfedge();
		for(int i=0; i<dv; ++i) {
			d2 = deg[i];
			d3 = deg[(i+dv-1)%dv];
			d4 = deg[(i+1)%dv];
			dk = deg[(i+shift)%dv];

			// Edge Flip
			int delR_EF = 4 + 2*(d3 + d4 - dv - d2);
			if (delR_EF < delRmin) {
				delRmin = delR_EF;
				op._e = h;
				op._type = Edge_Flip;
				op._delta_R = delR_EF;
			}
			//test: record drift edges
			if (delR_EF == 0)
			{
				drift_edges.push_back(h);
			}
			// Vertex Split
			int delR_VS = 2*(dv - xy + d2 + dk - 12);
			if (delR_VS < delRmin)
			{
				delRmin = delR_VS;
				op._e = h;
				op._type = Vertex_Split;
				op._delta_R = delR_VS;
			}
			// Edge Collapse
			int delR_EC = 4 + 2*(dv-5)*(d2-5) + 2*(dv+d2-d3-d4);
			if (delR_EC < delRmin)
			{
				delRmin = delR_EC;
				op._e = h;
				op._type = Edge_Collapse;
				op._delta_R = delR_EC;
			}
			h = h->next_around_vertex();
		}

		//test: drift edge is always acceptable
		if (op._delta_R >0 || (op._delta_R == 0 && op._type != Edge_Flip))
		{
			if (!drift_edges.empty())
			{
				op._e = drift_edges[0];
				op._type = Edge_Flip;
				op._delta_R = 0;
			}
			else {
				op._type = Null_Operator;
			}
		}
		//

		if (op._type == Null_Operator) return false;
		else {  //crease check
			Map::Vertex *v2 = op._e->opposite()->vertex();
			int cnt = 0;
			switch (op._type)
			{
			case Edge_Flip:   
				if (op._e->is_sharp(dihedral_dot)) return false;
				break;
			case Edge_Collapse:  //
				if (op._v->is_corner(dihedral_dot) || v2->is_corner(dihedral_dot) ||
					(!op._e->is_sharp(dihedral_dot) && (op._v->is_sharp(dihedral_dot) || v2->is_sharp(dihedral_dot))))
				{
					return false;
				}
				break;
			case Vertex_Split:
				if (op._v->is_sharp(dihedral_dot)) return false;
				break;
			default:
				std::cout << "Error: operation type not defined!" << std::endl;
				return false;
			}
			return true;

		}
	}

	bool DelaunayCVT::execute_map_topo_operation(MapEditorExt& ed, TopoOperation& op, std::vector<Map::Vertex*>& to_smooth) {
		bool succeed = false;
		to_smooth.clear();

		if (op.type() == Edge_Flip)
		{
			if(ed.edge_flip(op.edge())) {
				//sufficient
				to_smooth.push_back(op.vertex());  //v1
				to_smooth.push_back(op.edge()->vertex());  //v4
				to_smooth.push_back(op.edge()->prev()->vertex()); //v3
				to_smooth.push_back(op.edge()->opposite()->next()->vertex()); //v2
				//return true;
				succeed = true;
			}
		}
		else if (op.type() == Edge_Collapse)
		{
			to_smooth.push_back(op.vertex()); //v1
			to_smooth.push_back(op.edge()->next()->vertex()); //v3
			to_smooth.push_back(op.edge()->opposite()->prev()->prev()->vertex()); //v4
			if(ed.collapse_edge(op.edge())) { //collapse into op.vertex()
				//return true;
				succeed = true;
			}
			else to_smooth.clear();
		}
		else if (op.type() == Vertex_Split) {
			to_smooth.push_back(op.vertex()); //v
			to_smooth.push_back(op.edge()->opposite()->vertex());
			to_smooth.push_back(op.counter_edge()->opposite()->vertex());
			//
			bool success = true;
			Map::Halfedge* next = op.edge()->next_around_vertex()->next_around_vertex();
			Map::Halfedge* hend = op.counter_edge();
			//join faces
			while (next->prev_around_vertex() != hend)
			{
				if (ed.join_facets(next->prev_around_vertex())) {
					next = next->next_around_vertex();
				}
				else {
					success = false;
					break;  //TODO: check feasibility before actually executing operation
				}
			}
			//add center vertex
			if (success)
			{
				to_smooth.push_back(ed.create_center_vertex(hend->facet()));
				//return true;
				succeed = true;
			}
			else {
				to_smooth.clear();
			}
		}

		///test: remove feature vertex from to_smooth
		if (lock_feature_) {
			std::vector<Map::Vertex*> new_to_smooth;
			for (auto vit=to_smooth.begin(); vit!=to_smooth.end(); ++vit) {
				int vid = -1;
				if (delaunay_map_vert_index_.find(*vit) != delaunay_map_vert_index_.end()) vid = delaunay_map_vert_index_[*vit];
				if ((vid==-1) || (!locked_[vid] && !on_feature_[vid])) new_to_smooth.push_back(*vit);
			}
			to_smooth = new_to_smooth;
		}
		///

		//return false;
		return succeed;
	}

	void DelaunayCVT::map_laplace_smooth(std::vector<Map::Vertex*>& to_smooth) {
		std::vector<vec3> newp(to_smooth.size());
		newp.reserve(to_smooth.size());

		int max_iter = 10;
		float w = 0.2;  //step factor
		for (int i=0; i<max_iter; ++i)
		{
			//calc new vertex position
			for (int j=0; j<to_smooth.size(); ++j)
			{
				Map::Vertex* vj = to_smooth[j];
				newp[j] = vec3(0, 0, 0);
				Map::Halfedge* h = vj->halfedge();
				do 
				{
					newp[j] += h->opposite()->vertex()->point();
					h = h->next_around_vertex();
				} while (h != vj->halfedge());
				newp[j] /= vj->degree();  //center
				vec3 N = Geom::vertex_normal(vj);
				vec3 p = vj->point();
				vec3 v = newp[j] - p;
				newp[j] = p + w * (v - Geex::dot(v, N) * N);
			}
			//set new vertex position
			for (int j=0; j<to_smooth.size(); ++j)
			{
				to_smooth[j]->set_point(newp[j]);
			}
		}
	}

	void DelaunayCVT::map_topo_optimization(int nb_iter, int deltaR_threshold /*=0*/, double dihedral_dot_threshold /*= -1.0*/) {
		//calc min deltaR
		int delR_min = 1;
		TopoOperation op;
		if (lock_feature_) {
			gx_assert(lock_corners_);
			FOR_EACH_VERTEX(Map, delaunay_map_, vit) {
				if (find_topo_operation_preserve_fea(vit, op)) {
					delR_min = (op.delta_R() < delR_min) ? op.delta_R() : delR_min;
				}
			}
		} 
		else {
			FOR_EACH_VERTEX(Map, delaunay_map_, vit) {
				if (find_topo_operation(vit, op, dihedral_dot_threshold)) {
					delR_min = (op.delta_R() < delR_min) ? op.delta_R() : delR_min;
				}
			}
		}
		

		//
		MapEditorExt meditor(delaunay_map_);
		std::vector<TopoOperation> S;     //set of operations to execute
		std::set<Map::Vertex*> affected;  //vertices that are affected by some operation
		int iter = 0;
		while (nb_iter>0 && delR_min < deltaR_threshold)
		{
			int current_min = 1;   //min deltaR of current iter
			//collect operations
			FOR_EACH_VERTEX(Map, delaunay_map_, vit) {
				bool found = (lock_feature_) ? (find_topo_operation_preserve_fea(vit, op)) : find_topo_operation(vit, op, dihedral_dot_threshold);
				//if (find_topo_operation(vit, op, dihedral_dot_threshold))
				if (found)
				{
					current_min = (op.delta_R() < current_min) ? op.delta_R() : current_min;
					if (op.delta_R() <= delR_min)
					{
						//check if the define region of op is affected by other operations
						if (affected.find(op.vertex()) == affected.end())
						{
							bool is_define_region_clean = true;
							std::vector<Map::Vertex*> clean_vertices;
							clean_vertices.reserve(op.vertex()->degree() + 1);
							clean_vertices.push_back(op.vertex());

							if (op.type() == Vertex_Split)  //the define region is v and all its neighbors
							{
								Map::Halfedge* h = op.edge();
								do 
								{
									if (affected.find(h->prev()->vertex()) != affected.end())
									{
										is_define_region_clean = false;
										break;
									}
									else {
										clean_vertices.push_back(h->prev()->vertex());
									}
									h = h->next_around_vertex();
								} while (h != op.edge());
							}
							else if ((op.type() == Edge_Flip) || (op.type() == Edge_Collapse)) {  //the define region is 4 vertices of the two adjacent triangles
								Map::Vertex* v2 = op.edge()->prev()->vertex();
								Map::Vertex* v3 = op.edge()->prev_around_vertex()->prev()->vertex();
								Map::Vertex* v4 = op.edge()->next_around_vertex()->prev()->vertex();
								if (affected.find(v2) == affected.end()
									&& affected.find(v3) == affected.end()
									&& affected.find(v4) == affected.end())
								{
									clean_vertices.push_back(v2);
									clean_vertices.push_back(v3);
									clean_vertices.push_back(v4);
								}
								else {
									is_define_region_clean = false;
								}
							}
							//
							if (is_define_region_clean)
							{
								S.push_back(op);
								for (auto vit = clean_vertices.begin(); vit != clean_vertices.end(); ++vit)
								{
									affected.insert(*vit);
								}
							}
						}
					}
				}
			}

			//execute operations
			int op_count = 0;
			std::vector<Map::Vertex*> to_smooth;
			for (auto oit = S.begin(); oit != S.end(); ++oit)
			{
				if(execute_map_topo_operation(meditor, *oit, to_smooth)) {
					map_laplace_smooth(to_smooth);
					op_count++;
				}
			}

			//
			iter++;
			std::cout << "iter " << iter << ": " << S.size() << " found, " << op_count << " executed, delR = " << -current_min << std::endl;
			S.clear();
			affected.clear();
			delR_min = current_min;
			nb_iter--;
		}
	}

	void DelaunayCVT::collect_stretched_edges(std::vector<EdgeValue>& edges) {
		update_map_edge_stretch();
		edges.clear();
		FOR_EACH_HALFEDGE(Map, delaunay_map_, eit) {
			if (eit->is_edge_key())
			{
				Map::Vertex* v1 = eit->prev()->vertex();
				Map::Vertex* v2 = eit->vertex();
				if (v1->degree() != 6 || v2->degree() != 6)
				{
					double stretch = delaunay_map_edge_stretch_.at(eit);
					edges.push_back(std::make_pair(eit, stretch));
				}
			}
		}
	}



	void DelaunayCVT::stretch_topo_optimize(unsigned int nb_iter) {
		std::vector<EdgeValue> edges;

		for (unsigned int it=0; it<nb_iter; ++it)
		{
			collect_stretched_edges(edges);
			std::sort(edges.begin(), edges.end(), EdgeValueLess());
			std::vector<bool> is_visited(vertices_.size(), false);
			int nb_edges = edges.size();
			float nb_select = nb_edges * stretch_optimize_select_rate_;
			nb_select /= 2;
			if (nb_select < 1) return;

			std::vector<bool> to_remove(vertices_.size(), false);
			std::vector<vec3> to_insert;
			//split long
			int nb_selected = 0;
			for (int i=0; i<nb_edges; ++i)
			{
				EdgeValue& ei = edges[nb_edges - 1 - i];
				if (ei.second <= 0) break;
				else {
					Map::Halfedge_iterator eit = ei.first;
					Map::Vertex* v1 = eit->prev()->vertex();
					Map::Vertex* v2 = eit->vertex();
					int id1 = delaunay_map_vert_index_[v1];
					int id2 = delaunay_map_vert_index_[v2];
					if(is_visited[id1] || is_visited[id2]) continue;
					else {
						is_visited[id1] = true;
						is_visited[id2] = true;
						to_insert.push_back(0.5 * (v1->point() + v2->point()));  //insert middle point
						nb_selected++;
						if (nb_selected >= nb_select) break;
					}
				}
			}
			//collapse short
			nb_selected = 0;
			for (int i=0; i<nb_edges; ++i)
			{
				EdgeValue& ei = edges[i];
				if (ei.second >= 0) break;
				else {
					Map::Halfedge_iterator eit = ei.first;
					Map::Vertex* v1 = eit->prev()->vertex();
					Map::Vertex* v2 = eit->vertex();
					int id1 = delaunay_map_vert_index_[v1];
					int id2 = delaunay_map_vert_index_[v2];
					if(is_visited[id1] || is_visited[id2]) continue;
					else {
						is_visited[id1] = true;
						is_visited[id2] = true;
						to_remove[id1] = true;
						to_remove[id2] = true;
						to_insert.push_back(0.5 * (v1->point() + v2->point()));  //insert middle point
						nb_selected++;
						if (nb_selected >= nb_select) break;
					}
				}
			}

			//insert & remove vertices
			std::vector<vec3>& new_vertices = to_insert;
			std::vector<bool> new_locked(to_insert.size(), false);
			std::vector<bool> new_on_feature(to_insert.size(), false);
			for (int i=0; i<vertices_.size(); ++i)
			{
				if (!to_remove[i])
				{
					new_vertices.push_back(vertices_[i]);
					new_locked.push_back(locked_[i]);
					new_on_feature.push_back(on_feature_[i]);
				}
			}
			vertices_ = new_vertices;
			locked_ = new_locked;
			on_feature_ = new_on_feature;

			//smoothing
			newton_lloyd_global(stretch_optimize_nb_smooth_iter_);

			//info
			int nb_singular;
			calc_singularities(nb_singular);
			std::cout << "stretch " << it << ": ";
			std::cout << nb_singular << " singularities, "<< edges.size() << " stretched edges." << std::endl;

			is_delaunay_map_dirty_ = true;
		}
	}




	bool DelaunayCVT::quad_dominant_topo_optimize()
	{
		if (field_rot_symmetry_!=4) return false;
		const std::map<std::pair<int, int>, std::pair<double, int>> &edge_dir_match = primal_edge_field_match_;
		if (edge_dir_match.empty()) return false;
		int nvert = nb_vertices();

		QDmeshEditor qeditor(edge_dir_match, vertices_);


		// 		//sort incident edges by direction angle
		// 		std::vector<std::map<double, int>> QD0(nvert);  //vec<map<dir_angle, vidx>>
		// 		for (auto it=edge_dir_match.begin(); it!=edge_dir_match.end(); ++it)
		// 		{
		// 			int dOrder = it->second.second;
		// 			if (dOrder==0)
		// 			{
		// 				int vi = it->first.first;
		// 				int vj = it->first.second;
		// 				double dAngle = it->second.first;
		// 				gx_assert(dAngle>=0 && dAngle<2*M_PI);
		// 				QD0[vi].insert(std::make_pair(dAngle, vj));
		// 			}
		// 		}
		// 
		// 		//typedef std::pair<int, std::pair<int, int>> MEdge; //matched edge: <vidx, <over, loss>>
		// 		//std::vector<std::vector<MEdge>> QDmesh(nvert); //Quad Dominant Mesh
		// 		std::vector<std::vector<int>> QD1(nvert); //vec<vec<vidx>>
		// 		std::vector<std::map<int, int>> QEdges(nvert);  //
		// 		for (int i=0; i<nvert; ++i)
		// 		{
		// 			auto &map_ang_id = QD0[i];
		// 			for (auto it=map_ang_id.begin(); it!=map_ang_id.end(); ++it) //it: <dir_angle, vidx>
		// 			{
		// 				QD1[i].push_back(it->second);
		// 				QEdges[i].insert(std::make_pair(QD1[i].back(), QD1[i].size()-1));
		// 
		// 			}
		// 			gx_assert(QD1[i].size()==4); //4 edges per vertex
		// 		}
		// 
		// 		//compute opposite edge index and collect broken edges
		// 		typedef std::pair<int, int> Edge;
		// 		Edge nullEdge(-1,-1);
		// 		std::map<Edge, Edge> brokenEdges;
		// 		std::vector<std::map<int, int>> OpEdges(nvert);
		// 		for (int i=0; i<nvert; ++i) {
		// 			for (auto it=QD1[i].begin(); it!=QD1[i].end(); ++it) {
		// 				int j = *it;
		// 				OpEdges[i][j] = (QEdges[j].find(i)==QEdges[j].end()) ? -1 : QEdges[j][i];
		// 				if (OpEdges[i][j]==-1) brokenEdges.insert(std::make_pair(Edge(i, j), nullEdge));
		// 			}
		// 		}
		// 
		// 		//collect lack edges
		// 		std::map<Edge, Edge> lackEdges;
		// 		for (auto eit=brokenEdges.begin(); eit!=brokenEdges.end(); ++eit)
		// 		{
		// 			//find e1
		// 			int i = eit->first.first;
		// 			int j = eit->first.second;
		// 			int idx = QEdges[i][j];
		// 			for (int k=0; k<3; ++k) {
		// 				//ccw rotate
		// 				idx = (idx+1) % 4;
		// 				j = QD1[i][idx];
		// 				//flip
		// 				idx = OpEdges[i][j];
		// 				if (idx==-1) break;
		// 				else {
		// 					i = j;
		// 					j = QD1[i][idx];
		// 				}
		// 			}
		// 			if (idx==-1) break;
		// 			Edge e1(i, j);
		// 			//find e2
		// 			i = eit->first.first;
		// 			j = eit->first.second;
		// 			idx = QEdges[i][j];
		// 			for (int k=0; k<3; ++k) {
		// 				//cw rotate
		// 				idx = (idx+4-1) % 4;
		// 				j = QD1[i][idx];
		// 				//flip
		// 				idx = OpEdges[i][j];
		// 				if (idx==-1) break;
		// 				else {
		// 					i = j;
		// 					j = QD1[i][idx];
		// 				}
		// 			}
		// 			if (idx==-1) break;
		// 			Edge e2(i, j);
		// 			//check and insert
		// 			if (e1.first==e2.second && e1.second==e2.first) {
		// 				eit->second = e1;
		// 				lackEdges.insert(std::make_pair(e1, eit->first));
		// 			}			
		// 		}

		///Collect topo operations
		std::vector<vec3> to_insert;
		std::set<int> to_remove;
		//std::set<Edge> usedEdges;

		qeditor.collect_topo_operations(to_insert, to_remove);
		if (to_insert.empty() && to_remove.empty()) return false;

		// 		//vertex remove and edge collapse
		// 		for (auto eit=brokenEdges.begin(); eit!=brokenEdges.end(); ++eit)
		// 		{
		// 			Edge eij = eit->first;
		// 			if (usedEdges.find(eij)==usedEdges.end()) {
		// 				int i = eij.first;
		// 				int j = eij.second;
		// 				int idx = QEdges[i][j];
		// 				Edge e1(i, QD1[i][(idx+2) % 4]);
		// 				if ((usedEdges.find(e1)==usedEdges.end()) && (brokenEdges.find(e1)!=brokenEdges.end())) {
		// 					to_remove.insert(eij.first);
		// 					usedEdges.insert(eij);
		// 					usedEdges.insert(e1);
		// 					//mark corresponding lack edges
		// 					usedEdges.insert(brokenEdges.at(eij));
		// 					usedEdges.insert(brokenEdges.at(e1));
		// 				}
		// 				else { //e1 is not available, check neighbors
		// 					//left neighbor
		// 					bool left_succeed = false;
		// 					i = eij.first;
		// 					j = eij.second;
		// 					idx = QEdges[i][j];
		// 					//ccw rotate
		// 					idx = (idx+1) % 4;
		// 					j = QD1[i][idx];
		// 					//flip
		// 					idx = OpEdges[i][j];
		// 					if (idx!=-1) {
		// 						i = j;
		// 						j = QD1[i][idx];
		// 						//cw rotate
		// 						idx = (idx+4-1) % 4;
		// 						j = QD1[i][idx];
		// 						Edge left_e(i, j);
		// 						if ((usedEdges.find(left_e)==usedEdges.end()) && (OpEdges[i][j]==-1))
		// 						{
		// 							to_remove.insert(eij.first);
		// 							to_remove.insert(left_e.first);
		// 							to_insert.push_back(0.5 * (vertices_[eij.first] + vertices_[left_e.first]));
		// 							usedEdges.insert(eij);
		// 							usedEdges.insert(left_e);
		// 							//mark corresponding lack edges
		// 							usedEdges.insert(brokenEdges.at(eij));
		// 							usedEdges.insert(brokenEdges.at(left_e));
		// 							left_succeed = true;
		// 						}
		// 					}
		// 					if (!left_succeed) { //try right neighbor
		// 						i = eij.first;
		// 						j = eij.second;
		// 						idx = QEdges[i][j];
		// 						//cw rotate
		// 						idx = (idx+4-1) % 4;
		// 						j = QD1[i][idx];
		// 						//flip
		// 						idx = OpEdges[i][j];
		// 						if (idx!=-1) {
		// 							i = j;
		// 							j = QD1[i][idx];
		// 							//ccw rotate
		// 							idx = (idx+1) % 4;
		// 							j = QD1[i][idx];
		// 							Edge right_e(i, j);
		// 							if ((usedEdges.find(right_e)==usedEdges.end()) && (OpEdges[i][j]==-1)) {
		// 								to_remove.insert(eij.first);
		// 								to_remove.insert(right_e.first);
		// 								to_insert.push_back(0.5 * (vertices_[eij.first] + vertices_[right_e.first]));
		// 								usedEdges.insert(eij);
		// 								usedEdges.insert(right_e);
		// 								usedEdges.insert(brokenEdges.at(eij));
		// 								usedEdges.insert(brokenEdges.at(right_e));
		// 							}
		// 						}
		// 					}
		// 				}
		// 			}
		// 		}
		// 
		// 		//edge split and vertex split
		// 		for (auto eit=lackEdges.begin(); eit!=lackEdges.end(); ++eit)
		// 		{
		// 			Edge eij = eit->first;
		// 			if (usedEdges.find(eij)==usedEdges.end()) {
		// 				Edge eji(eij.second, eij.first);
		// 				if ((usedEdges.find(eji)==usedEdges.end()) && (lackEdges.find(eji)!=lackEdges.end())) {
		// 					to_insert.push_back(0.5 * (vertices_[eij.first] + vertices_[eij.second]));  //edge split
		// 					usedEdges.insert(eij);
		// 					usedEdges.insert(eji);
		// 					//usedEdges.insert(lackEdges.at(eij)); //since brokenedges will not be used, this is not necessary
		// 					//usedEdges.insert(lackEdges.at(eji));
		// 				}
		// 				else { //eji is not available, try neighbors
		// 					//left neighbor
		// 					bool left_succeed = false;
		// 					int i = eij.first;
		// 					int j = eij.second;
		// 					int idx = QEdges[i][j];
		// 					//2 x ccw rotate
		// 					idx = (idx+2) % 4;
		// 					j = QD1[i][idx];
		// 					Edge le(i, j);
		// 					if ((usedEdges.find(le)==usedEdges.end()) && (lackEdges.find(le)!=lackEdges.end())) {
		// 						to_remove.insert(eij.first);
		// 						to_insert.push_back(0.5 * (vertices_[eij.first] + vertices_[eij.second]));
		// 						to_insert.push_back(0.5 * (vertices_[le.first] + vertices_[le.second]));
		// 						usedEdges.insert(eij);
		// 						usedEdges.insert(le);
		// 						left_succeed = true;
		// 					}
		// 					//
		// 					if (!left_succeed) { //try right neighbor
		// 						i = eij.first;
		// 						j = eij.second;
		// 						idx = QEdges[i][j];
		// 						//flip
		// 						idx = OpEdges[i][j];
		// 						gx_assert(idx!=-1);
		// 						i = j;
		// 						j = QD1[i][idx];
		// 						//2 x ccw rotate
		// 						idx = (idx+2) % 4;
		// 						j = QD1[i][idx];
		// 						//flip
		// 						idx = OpEdges[i][j];
		// 						if (idx!=-1) {
		// 							i = j;
		// 							j = QD1[i][idx];
		// 							Edge re(i, j);
		// 							if ((usedEdges.find(re)==usedEdges.end()) && (lackEdges.find(re)!=lackEdges.end())) {
		// 								to_remove.insert(eij.second);
		// 								to_insert.push_back(0.5 * (vertices_[eij.first] + vertices_[eij.second]));
		// 								to_insert.push_back(0.5 * (vertices_[re.first] + vertices_[re.second]));
		// 								usedEdges.insert(eij);
		// 								usedEdges.insert(re);
		// 							}
		// 						}
		// 					}
		// 				}
		// 			}
		// 		}

		//insert & remove vertices
		std::vector<vec3> new_vertices = to_insert;
		std::vector<bool> new_locked(to_insert.size(), false);
		std::vector<bool> new_on_feature(to_insert.size(), false);
		for (int i=0; i<nvert; ++i)
		{
			if (to_remove.find(i)==to_remove.end())
			{
				new_vertices.push_back(vertices_[i]);
				new_locked.push_back(locked_[i]);
				new_on_feature.push_back(on_feature_[i]);
			}
		}
		vertices_ = new_vertices;
		locked_ = new_locked;
		on_feature_ = new_on_feature;

		//smoothing
		newton_lloyd_global(1);

		is_delaunay_map_dirty_ = true;
		return true;
	}

	//////////////////////////////////////////////////////////////////////////
	void DelaunayCVT::computeDirectionEnergyAndGradient(double &f, double *g, 
		const std::vector<Geex::vec3> &N, const std::vector<Geex::vec3> &D, 
		const std::vector<std::map<unsigned int, double>> &edge_weight, unsigned int RoSy)
	{
		int nvert = nb_vertices();

		for (int c=0; c<nvert; ++c)
		{
			const std::map<unsigned int, double> &neighbors = edge_weight[c];
			for (auto it=neighbors.begin(); it!=neighbors.end(); ++it)
			{
				int j = it->first;
				double len = it->second;
				double E1;
				vec3 G1;
				computeEdgeDirectionEG(E1, G1, vertices()[c], vertices()[j], N[c], D[c], RoSy);
				double E2 = E1 * len;
				vec3 G2 = G1 * len;

				//
				f += E2;
				g[3*c  ] += G2[0];
				g[3*c+1] += G2[1];
				g[3*c+2] += G2[2];
			}
		}
	}

	void DelaunayCVT::computeDirectionEnergyAndGradient(const std::vector<Geex::vec3> &N, const std::vector<Geex::vec3> &D, 
		const std::vector<std::map<unsigned int, double>> &edge_weight, unsigned int RoSy)
	{
		int nvert = nb_vertices();
		direction_energy_.assign(nvert, 0);
		direction_grad_.assign(nvert, vec3(0, 0, 0));

		for (int c=0; c<nvert; ++c)
		{
			const std::map<unsigned int, double> &neighbors = edge_weight[c];
			for (auto it=neighbors.begin(); it!=neighbors.end(); ++it)
			{
				int j = it->first;
				double weight = it->second;
				double E1;
				vec3 G1;
				//test
				computeEdgeDirectionEG(E1, G1, vertices()[j]-vertices()[c], N[c], D[c], RoSy);
				//computeEdgeDirectionEG(E1, G1, vertices()[c], vertices()[j], N[c], D[c], RoSy);
				double E2 = E1 * weight;
				vec3 G2 = G1 * weight;

				//
				direction_energy_[c] += E2;
				direction_grad_[c] += G2;
				// 				g[3*c  ] += G2[0];
				// 				g[3*c+1] += G2[1];
				// 				g[3*c+2] += G2[2];
			}
		}
	}


	void DelaunayCVT::compute_edge_dir_match(const std::vector<Geex::vec3> &N, const std::vector<Geex::vec3> &D, 
		const std::vector<std::map<unsigned int, double>> &edge_weight, unsigned int RoSy)
	{
		int nvert = nb_vertices();
		direction_energy_.assign(nvert, 0);
		direction_grad_.assign(nvert, vec3(0, 0, 0));
		primal_edge_field_match_.clear();
		primal_edge_match_history_.clear();

		//match maker
		int max_direction_order = 3;
		MatchMaker mm(RoSy, max_direction_order);
		//


		//compute match

		for (int c=0; c<nvert; ++c)
		{
			const std::map<unsigned int, double> &ew = edge_weight[c];
			std::vector<int> neighbor;
			//std::vector<double> eweight;
			for(auto it=ew.begin(); it!=ew.end(); ++it) {
				neighbor.push_back(it->first);
				//eweight.push_back(it->second);
			}

			//incident edges and their projections 
			vec3 n = N[c];
			vec3 d = D[c];
			std::vector<vec3> edge;
			std::vector<vec3> proj;
			for(int j=0; j<neighbor.size(); ++j) {
				edge.push_back(vertices()[neighbor[j]] - vertices()[c]);
				proj.push_back(edge[j] - dot(edge[j], n) * n);
			}

			//sort edges
			std::map<double, int> angleIdx;
			for (int j=0; j<neighbor.size(); ++j)
			{
				double ang_cos = dot(proj[j], d) / length(proj[j]);
				ang_cos = (ang_cos > 1) ? 1 : ang_cos;
				ang_cos = (ang_cos < -1) ? -1 : ang_cos;
				double ang_cross = dot(n, cross(d, proj[j]));
				double angle = (ang_cross >= 0) ? (acos(ang_cos)) : (2*M_PI - acos(ang_cos));  //angle belongs to [0, 2pi)
				angleIdx.insert(std::make_pair(angle, j));
			}

			std::vector<double> sorted_angle;
			std::vector<int> sorted_idx;
			for (auto it=angleIdx.begin(); it!=angleIdx.end(); ++it)
			{
				sorted_angle.push_back(it->first);
				sorted_idx.push_back(it->second);
			}

			//test
			double last_deg = sorted_angle.front();
			gx_assert(last_deg>=0 && last_deg<2*M_PI);
			for (auto it=sorted_angle.begin()+1; it!=sorted_angle.end(); ++it)
			{
				double deg = *it;
				gx_assert((deg>=last_deg) /*&& (deg>=0) && (deg<2*M_PI)*/);
				gx_assert(deg>=0);
				gx_assert(deg<2*M_PI);
				last_deg = deg;
			}

			//match
			std::vector<MatchMaker::Direction> dir_angle_order;
			std::vector<int> match_history;
			mm.match(sorted_angle, dir_angle_order, match_history);

			//compute energy & grad
			for (int j=0; j<neighbor.size(); ++j)
			{
				//double E1;
				//vec3 G1;
				int idx = sorted_idx[j];
				//computeEdgeDirectionEG(E1, G1, edge[idx], n, vec_rotate(d, n, dir_angle_order[j].first), /*1*/ RoSy);  //last parameter should by 1, but...
				//
				primal_edge_field_match_.insert(std::make_pair(std::make_pair(c, neighbor[idx]), dir_angle_order[j]));
				primal_edge_match_history_.insert(std::make_pair(std::make_pair(c, neighbor[idx]), match_history[j]));
				//
				//double ord_weight = 1.0 / (1+dir_angle_order[j].second);
				//E1 *= ord_weight;
				//G1 *= ord_weight;
				//
				//double weight = eweight[idx];
				//double E2 = E1 * weight;
				//vec3 G2 = G1 * weight;

				//direction_energy_[c] += E2;
				//direction_grad_[c] += G2;
			}
		}
	}

	void DelaunayCVT::computeDirectionEnergyAndGradient_match(const std::vector<Geex::vec3> &N, const std::vector<Geex::vec3> &D, 
		const std::vector<std::map<unsigned int, double>> &edge_weight, unsigned int RoSy)
	{
		int nvert = nb_vertices();
		direction_energy_.assign(nvert, 0);
		direction_grad_.assign(nvert, vec3(0, 0, 0));
		primal_edge_field_match_.clear();
		primal_edge_match_history_.clear();

		//match maker
		int max_direction_order = 3;
		MatchMaker mm(RoSy, max_direction_order);
		//


		//compute match

		for (int c=0; c<nvert; ++c)
		{
			const std::map<unsigned int, double> &ew = edge_weight[c];
			std::vector<int> neighbor;
			std::vector<double> eweight;
			for(auto it=ew.begin(); it!=ew.end(); ++it) {
				neighbor.push_back(it->first);
				eweight.push_back(it->second);
			}

			//incident edges and their projections 
			vec3 n = N[c];
			vec3 d = D[c];
			std::vector<vec3> edge;
			std::vector<vec3> proj;
			for(int j=0; j<neighbor.size(); ++j) {
				edge.push_back(vertices()[neighbor[j]] - vertices()[c]);
				proj.push_back(edge[j] - dot(edge[j], n) * n);
			}

			//sort edges
			std::map<double, int> angleIdx;
			for (int j=0; j<neighbor.size(); ++j)
			{
				double ang_cos = dot(proj[j], d) / length(proj[j]);
				ang_cos = (ang_cos > 1) ? 1 : ang_cos;
				ang_cos = (ang_cos < -1) ? -1 : ang_cos;
				double ang_cross = dot(n, cross(d, proj[j]));
				double angle = (ang_cross >= 0) ? (acos(ang_cos)) : (2*M_PI - acos(ang_cos));  //angle belongs to [0, 2pi)
				angleIdx.insert(std::make_pair(angle, j));
			}

			std::vector<double> sorted_angle;
			std::vector<int> sorted_idx;
			for (auto it=angleIdx.begin(); it!=angleIdx.end(); ++it)
			{
				sorted_angle.push_back(it->first);
				sorted_idx.push_back(it->second);
			}

			//test
			double last_deg = sorted_angle.front();
			gx_assert(last_deg>=0 && last_deg<2*M_PI);
			for (auto it=sorted_angle.begin()+1; it!=sorted_angle.end(); ++it)
			{
				double deg = *it;
				gx_assert((deg>=last_deg) /*&& (deg>=0) && (deg<2*M_PI)*/);
				gx_assert(deg>=0);
				gx_assert(deg<2*M_PI);
				last_deg = deg;
			}

			//match
			std::vector<MatchMaker::Direction> dir_angle_order;
			std::vector<int> match_history;
			mm.match(sorted_angle, dir_angle_order, match_history);

			//compute energy & grad
			for (int j=0; j<neighbor.size(); ++j)
			{
				double E1;
				vec3 G1;
				int idx = sorted_idx[j];
				computeEdgeDirectionEG(E1, G1, edge[idx], n, vec_rotate(d, n, dir_angle_order[j].first), /*1*/ RoSy);  //last parameter should by 1, but...
				//
				primal_edge_field_match_.insert(std::make_pair(std::make_pair(c, neighbor[idx]), dir_angle_order[j]));
				primal_edge_match_history_.insert(std::make_pair(std::make_pair(c, neighbor[idx]), match_history[j]));
				//
				double ord_weight = 1.0 / (1+dir_angle_order[j].second);
				E1 *= ord_weight;
				G1 *= ord_weight;
				//
				double weight = eweight[idx];
				double E2 = E1 * weight;
				vec3 G2 = G1 * weight;

				direction_energy_[c] += E2;
				direction_grad_[c] += G2;
			}
		}
	}

	/**
	* In polar coordinate,
	* E1 = (1 - sin(phi) * cos(RoSy * theta)) /2
	*/
	void  DelaunayCVT::computeEdgeDirectionEG(double &f, vec3 &g,
		const vec3& X1, const vec3& X2, const vec3& n, const vec3& d,
		unsigned int RoSy)
	{
		vec3 X1X2 = X2 - X1;
		vec3 proj = X1X2 - dot(X1X2, n) * n ;

		double phi = acos(dot(X1X2, n) / length(X1X2));
		double theta = 0.0;
		if (proj.length() > 0)
		{
			theta = acos(dot(proj, d) / length(proj));
			if (dot(cross(proj, d), n) > 0)
			{
				theta = 2 * M_PI - theta;
			}
		}

		//
		vec3 g_theta = cross(X1X2, n) / (length2(cross(X1X2, n)));
		vec3 g_phi   = cross(g_theta, X1X2) * length(cross(X1X2, n))/ length2(X1X2);

		//result
		if (RoSy==4 && use_special_quad_)
		{
			double para_q = direction_parameter_quad_;
			double f_quad = para_q * cos(4*theta) + (1-para_q) * cos(8*theta);
			double f_quad_min = para_q *para_q / (8*(para_q - 1)) + (para_q - 1);
			f_quad = (f_quad - f_quad_min) / (1 - f_quad_min);
			vec3 g_quad = (- 4 * para_q * sin(4*theta) - 8 * (1-para_q) * sin(8*theta))* g_theta / (1-f_quad_min); 
			f = 1 - sin(phi) * f_quad;
			g = - cos(phi) * f_quad * g_phi - sin(phi) * g_quad;
		}
		else {
			//f = 0.5 * (1 - sin(phi) * cos(RoSy*theta));
			f = 1 - sin(phi) * (0.5 * (1 + cos(RoSy*theta)));
			//g = -0.5 * cos(phi) * cos(RoSy*theta) * g_phi + 0.5 * RoSy * sin(phi) * sin(RoSy*theta) * g_theta;
			g = 0.5 * (RoSy * sin(RoSy*theta) * sin(phi) * g_theta - (1+cos(RoSy*theta)) * cos(phi) * g_phi);
		}

		//project g_E1 : solve the "shrink" problem
		g = g - dot(g, n) * n;
		//
	}


	/**
	* In polar coordinate,
	* E1 = (1 - sin(phi) * cos(RoSy * theta)) /2
	*/
	void  DelaunayCVT::computeEdgeDirectionEG(double &f, vec3 &g,
		const vec3& X1X2, const vec3& n, const vec3& d,
		unsigned int RoSy)
	{
		gx_assert(fabs(length(n) - 1) < 1e-8);
		vec3 proj = X1X2 - dot(X1X2, n) * n ;
		gx_assert(length2(X1X2) > 0);
		double cos_phi = dot(X1X2, n) / length(X1X2);
		cos_phi = (cos_phi > 1) ? 1: cos_phi;
		cos_phi = (cos_phi < -1) ? -1: cos_phi;
		//
		double phi = acos(cos_phi);
		double theta = 0.0;
		if (proj.length() > 0)
		{
			double cos_theta = dot(proj, d) / length(proj);
			cos_theta = (cos_theta > 1) ? 1 : cos_theta;
			cos_theta = (cos_theta < -1) ? -1 : cos_theta;
			theta = acos(cos_theta);
			if (dot(cross(proj, d), n) > 0)
			{
				theta = 2 * M_PI - theta;
			}
		}

		//
		gx_assert(length2(cross(X1X2, n)) > 0);
		vec3 g_theta = cross(X1X2, n) / (length2(cross(X1X2, n)));
		vec3 g_phi   = cross(g_theta, X1X2) * length(cross(X1X2, n))/ length2(X1X2);

		f = 1 - sin(phi) * (0.5 * (1 + cos(RoSy*theta)));
		//g = -0.5 * cos(phi) * cos(RoSy*theta) * g_phi + 0.5 * RoSy * sin(phi) * sin(RoSy*theta) * g_theta;
		g = 0.5 * (RoSy * sin(RoSy*theta) * sin(phi) * g_theta - (1+cos(RoSy*theta)) * cos(phi) * g_phi);

		//project g_E1 : solve the "shrink" problem
		g = g - dot(g, n) * n;
		//
	}

	void DelaunayCVT::computeNormalAndDirection(
		std::vector<vec3> &N, 
		std::vector<vec3> &D, 
		std::vector<unsigned int> &Idx,
		unsigned int RoSy /*= 6*/) 
	{
		int nvert = nb_vertices();
		const TopoPolyMesh* B = RVD_.mesh();
		void (*interpolator)(vec3 &vn, vec3 &vf, const vec3 &v, const vec3 pos[3],  const vec3 normal[3], const vec3 field[3], unsigned int RoSy);
		interpolator = (use_new_interpolation_) ? (normal_field_interpolate) : (simple_normal_field_interpolate);

		//test
		if (seed_interpolation_weight_.size() != nvert) seed_interpolation_weight_.resize(nvert);
		if (seed_project_triangle_.size() != nvert) seed_project_triangle_.resize(nvert);		
		//

		for (int c=0; c<nvert; ++c)
		{
			vec3 p[3];
			vec3 vnorm[3];
			vec3 field[3];
			unsigned int fidx = Idx[c];
			Geex::vec3g<unsigned int> vidx;
			for (int i=0; i<3; ++i)
			{
				vidx[i] = B->facet_begin(fidx) + i;
				vidx[i] = B->vertex_index(vidx[i]);
			}
			if (use_facet_field_)
			{
				N[c] = B->facet_normal()[fidx];
				D[c] = B->facet_field(fidx);
			}
			else {
				for (int i=0; i<3; ++i)
				{
					p[i] = (B->original_vertices())[vidx[i]];
					vnorm[i] = B->vertex_normal(vidx[i]);
					field[i] = B->vertex_field(vidx[i]);
				}
				interpolator(N[c], D[c], vertices_[c], p, vnorm, field, RoSy);
				//test: record interpolation weight and proj triangle
				vec3 w;
				compute_barycentric_weight(w, vertices_[c], p);
				seed_interpolation_weight_[c] = w;
				seed_project_triangle_[c].clear();
				for (int i=0; i<3; ++i) {
					seed_project_triangle_[c].push_back(p[i]);
				}
				//
			}
		}

	}

	//
	void DelaunayCVT::computeWeight(std::vector<double> &W, 
		const std::vector<std::map<unsigned int, std::vector<std::pair<Geex::TopoPolyVertexEdge, Geex::TopoPolyVertexEdge>>>> &dual_segments) 
	{
		int nvert = nb_vertices();

		for (int c=0; c<nvert; ++c)
		{
			int nseg = 0;
			double wsum = 0.0;
			const auto &neighbor = dual_segments[c];
			for (auto jt=neighbor.begin(); jt!=neighbor.end(); ++jt)
			{
				const auto &jseg = jt->second;
				for(int s=0; s<jseg.size(); ++s) {
					nseg += 2;
					wsum += (jseg[s].first.w + jseg[s].second.w);
				}
			}
			//average weight
			W[c] = wsum/nseg;
		}
	}

	//compute density weight for seed point
	void DelaunayCVT::computeWeight(
		std::vector<double> &W,
		const std::vector<unsigned int> &Idx) 
	{
		int nvert = nb_vertices();
		const TopoPolyMesh* B = RVD_.mesh();

		for (int c=0; c<nvert; ++c)
		{
			vec3 p[3];
			std::vector<double> vweight(3);

			unsigned int fidx = Idx[c];
			Geex::vec3g<unsigned int> vidx;
			for (int i=0; i<3; ++i)
			{
				vidx[i] = B->facet_begin(fidx) + i;
				vidx[i] = B->vertex_index(vidx[i]);
			}

			for (int i=0; i<3; ++i)
			{
				p[i] = (B->original_vertices())[vidx[i]];
				vweight[i] = B->vertex(vidx[i]).w;
			}

			vec3 weight;
			compute_barycentric_weight(weight, vertices_[c], p);

			W[c] = 0;
			//std::cout << "W[" << c << "] = ";
			for (int i=0; i<3; ++i) {
				W[c] += weight[i] * vweight[i];
				//std::cout << weight[i] << "*" << vweight[i] << " + ";
			}
			//std::cout << std::endl;
			//std::cout << "w[" << c << "] " << weight << std::endl;
		}
	}


	void DelaunayCVT::compute_dual_length(
		std::vector<std::map<unsigned int, double>> &res,
		const std::vector<std::map<unsigned int, std::vector<std::pair<TopoPolyVertexEdge, TopoPolyVertexEdge>>>> &dual_segments)
	{
		res.clear();
		res.resize(dual_segments.size());

		for(int i=0; i<dual_segments.size(); ++i) {
			auto &neighbor = dual_segments[i];
			for (auto jt = neighbor.begin(); jt!=neighbor.end(); ++jt)
			{
				int j = jt->first;
				res[i].insert(std::make_pair(j, 0));
				auto &seg = jt->second;
				for (int k=0; k<seg.size(); ++k)
				{
					res[i][j] += length(seg[k].first - seg[k].second);
				}
			}
		}
	}

	void DelaunayCVT::compute_approximated_lloyd_energy(
		std::vector<std::map<unsigned int, double>> &res,
		const std::vector<std::map<unsigned int, std::vector<std::pair<TopoPolyVertexEdge, TopoPolyVertexEdge>>>> &dual_segments,
		const std::vector<double> &v_weight)
	{
		res.clear();
		res.resize(dual_segments.size());

		for(int i=0; i<dual_segments.size(); ++i) {
			auto &neighbor = dual_segments[i];
			for (auto jt = neighbor.begin(); jt!=neighbor.end(); ++jt)
			{
				int j = jt->first;
				res[i].insert(std::make_pair(j, 0));
				auto &seg = jt->second;
				for (int k=0; k<seg.size(); ++k)
				{
					vec3 tgrad;
					double tV;
					if (use_density_)
					{
						res[i][j] += fabs(Lloyd_energy(vertices_[i], vertices_[i], seg[k].first, seg[k].second, 
							v_weight[i], seg[k].first.w, seg[k].second.w, /*1.0, 1.0, 1.0,*/
							tgrad, tV));
						//test: print density
						//std::cout << "w["<< i << "][" << j << "][" << k << "] " << v_weight[i] << ", " << seg[k].first.w << ", " << seg[k].second.w << std::endl;
						//
					}
					else {
						res[i][j] += fabs(Lloyd_energy(vertices_[i], vertices_[i], seg[k].first, seg[k].second, 
							/*v_weight[i], seg[k].first.w, seg[k].second.w,*/ 1.0, 1.0, 1.0,
							tgrad, tV));
					}
				}
			}
		}
	}

	void DelaunayCVT::compute_approx_cvt_energy(const std::vector<std::map<unsigned int, double>> & edge_energy)
	{
		approx_cvt_energy_.assign(edge_energy.size(), 0);
		for (int c=0; c<edge_energy.size(); ++c)
		{
			for (auto it=edge_energy[c].begin(); it!=edge_energy[c].end(); ++it)
			{
				approx_cvt_energy_[c] += it->second;
			}
		}
	}

	/**
	* Old interpolate algorithm: please read paper for details
	*/
	void DelaunayCVT::simple_normal_field_interpolate(
		vec3 &vn, vec3 &vf, const vec3 &v,
		const vec3 pos[3],  const vec3 normal[3], const vec3 field[3],
		unsigned int RoSy /*=6*/)
	{
		vec3 w;  //weight
		compute_barycentric_weight(w, v, pos);
		//compute_vertex_weight(w, v, pos, interpolationMode_);
		//test: print interpolation weight
		//std::cout << "w " << w << std::endl;
		// 		for (int i=0; i<3; ++i)
		// 		{
		// 			gx_assert(w[i]>=0 && w[i]<=1);
		// 		}
		//

		//normal interpolation
		vec3 wn = w[0]*normal[0] + w[1]*normal[1] + w[2]*normal[2];
		if (wn.length2() < 1e-10) {//zero vector wn
			vn = vec3(0, 0, 0);
			vf = vec3(0, 0, 0);
			std::cout << "Warning: Normal Degeneracy!" << std::endl;
			return;
		}
		else {
			vn = normalize(wn);
		}

		//field interpolation: old and simple way
		vec3 proj_field[3];
		for (int i=0; i<3; ++i)
		{
			proj_field[i] = field[i] - dot(field[i], vn) * vn;
			proj_field[i] = normalize(proj_field[i]);
		}
		vec3 e_x = proj_field[0];
		vec3 e_y = cross(vn, e_x);

		double theta[3];
		theta[0] = 0.0;
		for (int i=1; i<3; ++i) {
			double dot_i = dot(proj_field[i], e_x);
			dot_i = (dot_i>1)  ?  1 : dot_i;
			dot_i = (dot_i<-1) ? -1 : dot_i;
			theta[i] = acos(dot_i);
			if (dot(proj_field[i], e_y) < 0) {
				theta[i] = 2 * M_PI - theta[i];
			}
			theta[i] *= RoSy;
		}

		vec2 sym_direction[3];
		for (int i=0; i<3; ++i) {
			sym_direction[i] = vec2(cos(theta[i]), sin(theta[i]));
		}
		vec2 sym_dir = w[0] * sym_direction[0] + w[1] * sym_direction[1] + w[2] * sym_direction[2];
		sym_dir = normalize(sym_dir);
		double sym_angle = acos(sym_dir.x);
		if (sym_dir.y < 0) {
			sym_angle = 2 * M_PI - sym_angle;
		}
		double dir_angle = sym_angle / RoSy;

		vec3 dir = cos(dir_angle) * e_x + sin(dir_angle) * e_y;
		vf = normalize(dir);
	}

	/**
	* New interpolate algorithm: please read paper for details
	*/
	void DelaunayCVT::normal_field_interpolate(
		vec3 &vn, vec3 &vf, const vec3 &v,
		const vec3 pos[3],  const vec3 normal[3], const vec3 field[3],
		unsigned int RoSy /*=6*/)
	{
		vec3 w;  //weight
		compute_barycentric_weight(w, v, pos);
		//compute_vertex_weight(w, v, pos, interpolationMode_);


		//normal interpolation
		vec3 wn = w[0]*normal[0] + w[1]*normal[1] + w[2]*normal[2];
		if (wn.length2() < 1e-10) {//zero vector wn
			vn = vec3(0, 0, 0);
			vf = vec3(0, 0, 0);
			std::cout << "Warning: Normal Degeneracy!" << std::endl;
			return;
		}
		else {
			vn = normalize(wn);
		}

		//field interpolation: transverse all combination of representational field vector
		double maxlen2 = 0;
		vec3 maxf;
		for (unsigned int i=0; i<RoSy; ++i)
		{
			for (unsigned int j=0; j<RoSy; ++j)
			{
				for (unsigned int k=0; k<RoSy; ++k)
				{
					vec3 f0 = vec_rotate(field[0], normal[0], i*2*M_PI/RoSy);
					vec3 f1 = vec_rotate(field[1], normal[1], j*2*M_PI/RoSy);
					vec3 f2 = vec_rotate(field[2], normal[2], k*2*M_PI/RoSy);
					vec3 wf = w[0]*f0 + w[1]*f1 + w[2]*f2;
					vec3 proj_f = wf - dot(wf, vn)*vn;  //project weighted field direction
					if (length2(proj_f) > maxlen2)
					{
						maxlen2 = length2(proj_f);
						maxf = proj_f;
					}
				}
			}
		}

		if (length2(maxf) < 1e-10)
		{
			vf = vec3(0, 0, 0);
			std::cout << "Warning: Field Degeneracy!" << std::endl;
		}
		else {
			vf = normalize(maxf);
		}
		return;
	}

	void DelaunayCVT::compute_vertex_weight(vec3 &weight, const vec3 &p, const vec3 vert[3], GLenum mode /*= Barycentric*/)
	{
		switch (mode)
		{
		case Barycentric:
			compute_barycentric_weight(weight, p, vert);
			break;
		case Cotangent:
			compute_cotangent_weight(weight, p, vert);
			break;
		default:
			weight = vec3(0, 0, 0);
			std::cout << "Error: unrecognized interpolation weight mode!" << std::endl;
			break;
		}
	}

	//feature preservation
	void DelaunayCVT::compute_feature_edge_line() {
		if (features_.empty()) return;
		feature_edge_line_.clear();
		for (int i=0; i<features_.size(); ++i) {
			const auto &line = features_[i];
			for (int j=0; j<line.size()-1; ++j) {
				int v0 = line[j];
				int v1 = line[j+1];
				int a = (v0 < v1) ? v0 : v1;
				int b = (v0 < v1) ? v1 : v0; 
				feature_edge_line_.insert(std::make_pair(std::make_pair(a, b), i));
			}
		}
	}

	void DelaunayCVT::compute_seed_feature_lines(/*std::map<int, std::set<int>> &seed_feature_lines*/) {
		if (feature_edge_line_.empty()) compute_feature_edge_line();
		//compute 
		seed_feature_lines_.clear();
		RVD_.for_each_facet(ComputeSeedFeatureLine(seed_feature_lines_, feature_edge_line_, RVD_));
	}

	void DelaunayCVT::update_seed_feature_line() {
		seed_feature_line_.clear();
		on_feature_.clear();
		on_feature_.resize(vertices_.size(), false);
		compute_seed_feature_lines();
		std::vector<vec3> to_insert;
		std::vector<int> to_insert_line;
		for (auto it=seed_feature_lines_.begin(); it!=seed_feature_lines_.end(); ++it) {
			int vi = it->first;
			if (locked_[vi]) continue; //skip corners
			vec3 p = vertices_[vi];
			const auto &lines = it->second;
			//
			auto lt = lines.begin();
			seed_feature_line_.insert(std::make_pair(vi, *lt));
			on_feature_[vi] = true;
			//to insert
			++lt;
			for (; lt!=lines.end(); ++lt) {
				to_insert.push_back(p);
				to_insert_line.push_back(*lt);
			}
		}
		//insert new vertices
		for (int i=0; i<to_insert.size(); ++i) {
			vertices_.push_back(to_insert[i]);
			locked_.push_back(false);
			on_feature_.push_back(true);
			seed_feature_line_.insert(std::make_pair(vertices_.size()-1, to_insert_line[i]));
		}
	}

	///dxy add end

	void DelaunayCVT::compute_rvd() {
		delaunay_->set_vertices(vertices_) ;
		primal_dirty_ = GL_TRUE ;
		//dxy add
		is_delaunay_map_dirty_ = true;
	}

	//======================= Functors for Lloyd ===================

	void DelaunayCVT::set_vertices(int N, double* x) {
		int nv = N/3 ;
		gx_assert(nv == nb_vertices()) ;
		for(int i=0; i<nv; ++i) {
			vertices_[i][0] = x[i*3+0] ;
			vertices_[i][1] = x[i*3+1] ;
			vertices_[i][2] = x[i*3+2] ;
		}
		compute_rvd() ;
		//dxy add: delaunay_map
		is_delaunay_map_dirty_ = true;
	}

	//dxy add
	void DelaunayCVT::snap_to_corners(int N, double *x) {
		int nv = N/3;
		gx_assert(nv == nb_vertices());
		for (int i=0; i<nv; ++i) {
			if (locked_[i]) {
				gx_assert(vertices_corner_index_[i] != -1);
				vec3 p = verts_[corners_[vertices_corner_index_[i]]]->point();
				x[i*3+0] = p[0];
				x[i*3+1] = p[1];
				x[i*3+2] = p[2];
			}
		}
	}

	//dxy add: project feature points to its feature line, update feature tangent directions
	void DelaunayCVT::snap_to_feature(int N, double *x) {
		int nv = N/3;
		gx_assert(nv == nb_vertices());
		seed_feature_tan_.clear();
		for (auto it = seed_feature_line_.begin(); it !=seed_feature_line_.end(); ++it) {
			int i = it->first;
			int l = it->second;
			gx_assert (l < features_.size());
			vec3 pt(x[3*i], x[3*i+1], x[3*i+2]);
			const std::vector<int> &line = features_[l];
			double mindist2 = 1e10, dist2;
			vec3 fp; // foot point
			vec3 p;
			vec3 tan; //tangent dir
			for (int j=0; j<line.size()-1; ++j) {
				vec3 v0 = verts_[line[j]]->point();
				vec3 v1 = verts_[line[j+1]]->point();
				project_to_linesegment(pt, v0, v1, p, dist2);
				if (mindist2 > dist2) {
					mindist2 = dist2;
					tan = v1 - v0;
					fp = p;
				}
			}
			Geex::normalize(tan);
			x[3*i+0] = fp[0];
			x[3*i+1] = fp[1];
			x[3*i+2] = fp[2];
			seed_feature_tan_.insert(std::make_pair(i, tan));
		}		
	}

	void DelaunayCVT::align_feature_grad(int N, double *g) {
		int nv = N/3;
		gx_assert(nv == nb_vertices());
		for (auto it=seed_feature_tan_.begin(); it!=seed_feature_tan_.end(); ++it) {
			int i = it->first;
			vec3 tan = it->second;
			vec3 gi(g[3*i], g[3*i+1], g[3*i+2]);
			vec3 proj_gi = Geex::dot(gi, tan) * tan;
			g[3*i] = proj_gi[0];
			g[3*i+1] = proj_gi[1];
			g[3*i+2] = proj_gi[2];
		}
	}
	//

	vec3 DelaunayCVT::project_to_crease(vec3& pt, vec3& tan) {
		int idx ;
		double mindist2 = 1e10 , dist2;
		vec3 fp ; // foot point 
		vec3 p ;
		for(unsigned int i=0; i<features_.size(); ++i) {
			std::vector<int>& line = features_[i] ;
			for(unsigned int j=0; j<line.size()-1; ++j) {
				vec3 v0 = verts_[line[j]]->point() ;
				vec3 v1 = verts_[line[j+1]]->point() ;
				project_to_linesegment(pt, v0, v1, p, dist2) ;
				if(mindist2 > dist2) {
					mindist2 = dist2 ;
					tan = v1 - v0 ;
					fp = p ;
				}
			}
		}
		Geex::normalize(tan) ;
		return fp ;
	}

	//unsigned DelaunayCVT::nb_small() {
	//	unsigned nb = 0 ;
	//	//for(unsigned i=0; i<primal_facets_.size(); ++i) {
	//	//	PrimalTri& t = primal_facets_[i] ;
	//	//	double amin, amax ;
	//	//	min_max_angle(vertices_[t[0]], vertices_[t[1]], vertices_[t[2]], amin, amax) ;
	//	//	if(amin < angle_min_)
	//	//		nb ++ ;
	//	//}
	//	return nb ;
	//}

	//unsigned DelaunayCVT::nb_obtuses() {
	//	unsigned nb = 0 ;
	//	for(unsigned i=0; i<primal_facets_.size(); ++i) {
	//		if(is_obtuse(primal_facets_[i])) {
	//			nb++ ;
	//		}
	//	}
	//	return nb ;
	//}

	//unsigned DelaunayCVT::is_obtuse(PrimalTri& t) {
	//	return Geex::is_obtuse(vertices_[t[0]], vertices_[t[1]], vertices_[t[2]]) ;
	//}

	//	void DelaunayCVT::update_valences() {
	//		std::map<unsigned, std::set<unsigned>> neigh ;
	//		adjtris_.clear() ;
	////		update_primal() ;
	//		for(unsigned i=0; i<primal_facets_.size(); ++i) {
	//			PrimalTri& t = primal_facets_[i] ;
	//			neigh[t[0]].insert(t[1]) ;
	//			neigh[t[0]].insert(t[2]) ;
	//			neigh[t[1]].insert(t[0]) ;
	//			neigh[t[1]].insert(t[2]) ;
	//			neigh[t[2]].insert(t[1]) ;
	//			neigh[t[2]].insert(t[0]) ;
	//
	//			adjtris_[t[0]].insert(i) ;
	//			adjtris_[t[1]].insert(i) ;
	//			adjtris_[t[2]].insert(i) ;
	//		}
	//
	//		int nblt3=0, nblt5=0, nbgt7=0, n5=0, n6=0, n7=0 ;
	//		int nv = vertices_.size() ;
	//		valences_.resize(nv) ;
	//		for(unsigned i=0; i<nv; ++i) {
	//			valences_[i] = neigh[i].size() ;
	//			if(valences_[i]<3) nblt3++ ;
	//			else if(valences_[i]<5) nblt5++ ;
	//			else if(valences_[i]>7) nbgt7++ ;
	//			else {
	//				switch(valences_[i]) {
	//				case 5: n5++ ; break ;
	//				case 6: n6++ ; break ;
	//				case 7: n7++ ; break ;
	//				default: break ;
	//				} ;
	//			}
	//		}
	//
	//		std::cout << "nb v<3: " << nblt3 <<", nb v<5: " << nblt5 << ", nb v>7: " << nbgt7 << std::endl ;
	//	}

	///dxy add: rewrite update_valences()

	void DelaunayCVT::update_valences(bool message /* = false */) {
		std::map<unsigned, std::set<unsigned>> neigh ;
		RVD_.for_each_primal_triangle(ComputePrimalNeighborhood(neigh));

		int nblt3=0, nblt5=0, nbgt7=0, n5=0, n6=0, n7=0 ;
		unsigned int nv = vertices_.size() ;
		valences_.resize(nv) ;
		for(unsigned i=0; i<nv; ++i) {
			valences_[i] = neigh[i].size() ;
			if(valences_[i]<3) nblt3++ ;
			else if(valences_[i]<5) nblt5++ ;
			else if(valences_[i]>7) nbgt7++ ;
			else {
				switch(valences_[i]) {
				case 5: n5++ ; break ;
				case 6: n6++ ; break ;
				case 7: n7++ ; break ;
				default: break ;
				} ;
			}
		}

		if (message)
		{
			std::cout << "nb v<3: " << nblt3 <<", nb v<5: " << nblt5 << ", nb v>7: " << nbgt7 << std::endl ;
			std::cout << "nv v5: " << n5 << ", nb v6: " << n6 << ", nb v7: " << n7 << std::endl;
		}
	}

	///dxy add
	// 	void DelaunayCVT::update_primal_angles() {
	// 		primal_angles_.clear();
	// 		primal_min_angles_.clear();
	// 		primal_max_angles_.clear();
	// 		RVD_.for_each_primal_triangle(ComputePrimalAngles(vertices_, primal_angles_, primal_min_angles_, primal_max_angles_));
	// 	}

	void DelaunayCVT::update_primal_angles() {
		update_delaunay_map();
		if (!is_primal_angles_dirty_) return;

		primal_angles_.clear();
		primal_min_angles_.clear();
		primal_max_angles_.clear();
		FOR_EACH_FACET(Map, delaunay_map_, fit) {
			vec3 v1 = fit->halfedge()->prev()->vertex()->point();
			vec3 v2 = fit->halfedge()->vertex()->point();
			vec3 v3 = fit->halfedge()->next()->vertex()->point();
			vec3 ab = v2 - v1 ; //F.vertex[1] - F.vertex[0] ;
			vec3 bc = v3 - v2 ; //F.vertex[2] - F.vertex[1] ;
			vec3 ca = v1 - v3 ; //F.vertex[0] - F.vertex[2] ;
			double lab = length(ab) ;
			double lbc = length(bc) ;
			double lca = length(ca) ;
			double anglea = acos(-dot(ca, ab)/(lca*lab))*180. / M_PI ;
			double angleb = acos(-dot(ab, bc)/(lab*lbc))*180. / M_PI ;
			double anglec = acos(-dot(bc, ca)/(lbc*lca))*180. / M_PI ;

			primal_angles_.push_back(anglea);
			primal_angles_.push_back(angleb);
			primal_angles_.push_back(anglec);

			double min, max;
			if (anglea < angleb && anglea < anglec)
			{
				min = anglea;
				max = (angleb > anglec ? angleb : anglec);
			}
			else {
				if (angleb < anglec)
				{
					min = angleb;
					max = (anglea >  anglec ? anglea : anglec);
				}
				else {
					min = anglec;
					max = (anglea > angleb ? anglea : angleb);
				}
			}
			primal_min_angles_.push_back(min);
			primal_max_angles_.push_back(max);
		}
		is_primal_angles_dirty_ = false;
	}

	void DelaunayCVT::save_primal_angles(const std::string& filename) {
		std::cerr << "Writing primal angles to " << filename << std::endl;
		std::ofstream out(filename.c_str());

		DelaunayCVT::update_primal_angles();
		std::sort(primal_angles_.begin(), primal_angles_.end());
		for (int i=0; i<primal_angles_.size(); ++i) {
			out << "angle " << primal_angles_[i] << std::endl;
		}
		//min
		for (int i=0; i<primal_min_angles_.size(); ++i)
		{
			out << "min " << primal_min_angles_[i] << std::endl;
		}
		//max
		for (int i=0; i<primal_max_angles_.size(); ++i)
		{
			out << "max " << primal_max_angles_[i] << std::endl;
		}
		out.close();

		std::cerr << "min_angle: " << primal_angles_.front() << ", max_angle: " << primal_angles_.back() << std::endl;
		std::cerr << "Done." << std::endl;
	}

	//
	// 	void DelaunayCVT::update_primal_areas() {
	// 		primal_areas_.clear();
	// 		RVD_.for_each_primal_triangle(ComputePrimalAreas(vertices_, primal_areas_));
	// 	}

	void DelaunayCVT::update_primal_areas() {
		update_delaunay_map();
		if (!is_primal_areas_dirty_) return;

		primal_areas_.clear();
		vec3 v1, v2, v3;
		FOR_EACH_FACET(Map, delaunay_map_, fit) {
			v1 = fit->halfedge()->prev()->vertex()->point();
			v2 = fit->halfedge()->vertex()->point();
			v3 = fit->halfedge()->next()->vertex()->point();
			primal_areas_.push_back(tri_area(v1, v2, v3));
		}
		is_primal_areas_dirty_ = false;
	}

	void DelaunayCVT::save_primal_areas(const std::string& filename) {
		std::cerr << "Writing primal areas to " << filename << std::endl;
		std::ofstream out(filename.c_str());

		DelaunayCVT::update_primal_areas();
		std::sort(primal_areas_.begin(), primal_areas_.end());
		for (int i=0; i<primal_areas_.size(); ++i)
		{
			out << "area " << primal_areas_[i] << std::endl;
		}
		out.close();

		//
		double total_area = std::accumulate(primal_areas_.begin(), primal_areas_.end(), 0.0);
		double aver_area = total_area / primal_areas_.size();
		std::cerr << "tot area: " << total_area << ", aver area: " << aver_area << std::endl;
		std::cerr << "min_area: " << primal_areas_.front()/aver_area << ", max_area: " << primal_areas_.back()/aver_area << std::endl;
		std::cerr << "Done." << std::endl;

	}

	void DelaunayCVT::save_primal_map(const std::string& filename) {
		std::cerr << "Writing primal mesh to " << filename << std::endl ;
		update_delaunay_map();
		save_map(filename, delaunay_map_);
		std::cerr << "Done." << std::endl ;
	}

	void DelaunayCVT::save_primal_info(const std::string& filename) {
		std::cerr << "Writing primal info to " << filename << std::endl;
		std::ofstream out(filename.c_str());
		//Angle
		DelaunayCVT::update_primal_angles();
		std::sort(primal_angles_.begin(), primal_angles_.end());
		for (int i=0; i<primal_angles_.size(); ++i) {
			out << "angle " << primal_angles_[i] << std::endl;
		}
		//min
		for (int i=0; i<primal_min_angles_.size(); ++i)
		{
			out << "min_angle " << primal_min_angles_[i] << std::endl;
		}
		//max
		for (int i=0; i<primal_max_angles_.size(); ++i)
		{
			out << "max_angle " << primal_max_angles_[i] << std::endl;
		}
		//Area
		DelaunayCVT::update_primal_areas();
		std::sort(primal_areas_.begin(), primal_areas_.end());
		for (int i=0; i<primal_areas_.size(); ++i)
		{
			out << "area " << primal_areas_[i] << std::endl;
		}
		//Singularity
		int nb_singular;
		DelaunayCVT::calc_singularities(nb_singular);
		out << "nb_singular " << nb_singular << std::endl;
		//
		out << "nb_vert " << delaunay_map_->nb_vertices() << std::endl;
		out << "nb_face " << delaunay_map_->nb_facets() << std::endl;
		//
		out.close();
		std::cout << "Done." << std::endl;
	}

	void DelaunayCVT::print_primal_info(bool angle, bool area, bool singularity) {
		update_delaunay_map();
		if (angle)
		{
			update_primal_angles();
			if (!primal_angles_.empty())
			{
				double sum = 0.0, min = 181, max = -1;
				for (int i=0; i<primal_angles_.size(); ++i)
				{
					sum += primal_angles_[i];
					if (primal_angles_[i] < min)  min = primal_angles_[i];
					if (primal_angles_[i] > max)  max = primal_angles_[i];
				}
				std::cout << "Angles: min " << min << ", aver " << sum/primal_angles_.size() << ", max " << max << std::endl;
			}
		}
		if (area)
		{
			update_primal_areas();
			if (!primal_areas_.empty()) {
				double sum = primal_areas_[0];
				double min = sum, max = sum;
				for (int i=1; i<primal_areas_.size(); ++i)
				{
					sum += primal_areas_[i];
					min = (primal_areas_[i] < min) ? primal_areas_[i] : min;
					max = (primal_areas_[i] > max) ? primal_areas_[i] : max;
				}
				double aver = sum / primal_areas_.size();
				std::cout << "Area: min " << min/aver << ", max " << max/aver << std::endl;
			}
		}
		if (singularity)
		{
			int nb;
			calc_singularities(nb);
			std::cout << "nb singularity: " << nb << std::endl;
		}
	}

	//test
	void DelaunayCVT::check_map() {
		update_delaunay_map();
		int nb_map_vert = delaunay_map_->nb_vertices();
		//int nb_map_face = delaunay_map_->nb_facets();

		//angle
		update_primal_angles();
		std::vector<double> angles, min_angles, max_angles;
		RVD_.for_each_primal_triangle(ComputePrimalAngles(vertices_, angles, min_angles, max_angles));

		auto map_min_max = std::minmax_element(primal_angles_.begin(), primal_angles_.end());
		double map_aver = std::accumulate(primal_angles_.begin(), primal_angles_.end(), 0);
		map_aver /= primal_angles_.size();
		double map_min_aver = std::accumulate(primal_min_angles_.begin(), primal_min_angles_.end(), 0);
		map_min_aver /= primal_min_angles_.size();
		double map_max_aver = std::_Accumulate(primal_max_angles_.begin(), primal_max_angles_.end(), 0);
		map_max_aver /= primal_max_angles_.size();

		auto rvd_min_max = std::minmax_element(angles.begin(), angles.end());
		double rvd_aver = std::accumulate(angles.begin(), angles.end(), 0);
		rvd_aver /= angles.size();
		double rvd_min_aver = std::accumulate(min_angles.begin(), min_angles.end(), 0);
		rvd_min_aver /= min_angles.size();
		double rvd_max_aver = std::_Accumulate(max_angles.begin(), max_angles.end(), 0);
		rvd_max_aver /= max_angles.size();

		//print
		std::cout.precision(5);
		std::cout << "    | vert | face |  min   |aver_min|  max  |aver_max|  aver  |" << std::endl;
		std::cout << "RVD | " << vertices_.size() << " | " << min_angles.size()         << " | " << *rvd_min_max.first << " | " << rvd_min_aver << " | " << *rvd_min_max.second << " | " << rvd_max_aver << " | " << rvd_aver << " | " << std::endl;
		std::cout << "Map | " << nb_map_vert      << " | " << primal_min_angles_.size() << " | " << *map_min_max.first << " | " << map_min_aver << " | " << *map_min_max.second << " | " << map_max_aver << " | " << map_aver << " | " << std::endl;		
	}
	///

	void DelaunayCVT::save_primal_global(const std::string& filename) {
		std::cerr << "Writing primal mesh to " << filename << std::endl ;
		std::ofstream out(filename.c_str()) ;
		for(unsigned int i=0; i<vertices_.size(); i++) {
			out << "v " << vertices_[i] << std::endl ;
		}
		RVD_.for_each_primal_triangle(SavePrimalTriangle(out)) ;
		out.close();
		std::cerr << "Done." << std::endl ;
	}

	///dxy add
	void DelaunayCVT::save_denormalized_boundary(const std::string& filename) {
		std::cerr << "Writing origin mesh to " << filename << std::endl;
		std::ofstream out(filename.c_str()) ;
		double scale = 1.0 / model_normalize_scale_;
		for(unsigned int i=0; i<boundary_->nb_vertices(); i++) {
			vec3 newp = scale*(boundary_->vertex(i)) + vec3(x_min_, y_min_, z_min_) ;
			out << "v " << newp << std::endl ;
		}
		for(unsigned int i=0; i<boundary_->nb_facets(); i++) {
			out << "f" ;
			for(unsigned int j=boundary_->facet_begin(i); j<boundary_->facet_end(i); ++j) {
				out <<" " << j+1 ;
			}
			out << std::endl ;
		}
		out.close();
		std::cerr << "Done." << std::endl ;
	}

	void DelaunayCVT::save_denormalized_primal(const std::string& filename) {
		std::cerr << "Writing primal mesh to " << filename << std::endl ;
		double scale = 1.0 / model_normalize_scale_;
		std::ofstream out(filename.c_str()) ;
		for(unsigned int i=0; i<vertices_.size(); i++) {
			vec3 newp = scale*(vertices_[i]) + vec3(x_min_, y_min_, z_min_) ;
			out << "v " << newp << std::endl ;
		}
		RVD_.for_each_primal_triangle(SavePrimalTriangle(out)) ;
		out.close();
		std::cerr << "Done." << std::endl ;
	}

	void DelaunayCVT::save_optimization_info(const std::string& filename) {
		std::cerr << "Writing optimization information to " << filename << std::endl;
		std::ofstream out(filename.c_str());
		out << "Vert \t " << nb_vertices() << std::endl;
		out << "Time(ms) \t CVT \t " << optimize_cvt_time << " \t Dir \t " << optimize_dir_time+optimize_interpolate_time << std::endl; 
		out << "iter \t Ecvt \t Gcvt \t Edir \t Gdir \t Gtot" << std::endl;
		for(int i=0; i<optimize_cvt_energy.size(); ++i) {
			out << i << " \t " << optimize_cvt_energy[i]
			<< " \t " << optimize_cvt_gnorm[i]
			<< " \t " << optimize_dir_energy[i]
			<< " \t " << optimize_dir_gnorm[i]
			<< " \t " << optimize_tot_gnorm[i] << std::endl;
		}
		out.close();
		std::cerr << "Done." << std::endl ;
	}
	///dxy add end

	void DelaunayCVT::save_normalized_boundary(const std::string& filename) {
		std::cerr << "Writing normalized mesh " << filename << std::endl ;
		save_map(filename, map_) ;
		//std::ofstream out(filename.c_str()) ;
		//for(unsigned int i=0; i<boundary_->nb_vertices(); i++) {
		//	out << "v " << boundary_->vertex(i) << std::endl ;
		//}
		//for(unsigned int i=0; i<boundary_->nb_facets(); i++) {
		//	out << "f" ;
		//	for(unsigned int j=boundary_->facet_begin(i); j<boundary_->facet_end(i); ++j) {
		//		 out <<" " << j+1 ;
		//	}
		//	out << std::endl ;
		//}
		//out.close();
		std::cerr << "Done." << std::endl ;		
	}

	//========================================================

	void DelaunayCVT::compute_curvature_tensor(
		bool on_facets,
		bool relative,
		bool anisotropic,
		int aniso_iters,
		double aniso_factor,
		double neighborhood_size,
		const std::string& K1_attribute_name,
		const std::string& K2_attribute_name,
		const std::string& N_attribute_name,
		const std::string& k1_attribute_name,
		const std::string& k2_attribute_name,
		const std::string& n_attribute_name,
		bool compute_hk,							     
		bool copy_n_to_normal
		) {

			if(relative) {
				Box3d bbox = Geom::get_map_bbox(map_) ;
				neighborhood_size *= (2.0 * bbox.radius()) ;
			}

			MapCurvature curvatures(map_) ;
			curvatures.set_radius(neighborhood_size) ;

			if(on_facets) {
				curvatures.compute_curvature_tensor_on_facets() ;
				//				map_->compute_normals() ; //update() ;
				return ;
			}

			curvatures.set_anisotropic(anisotropic) ;
			curvatures.set_nb_anisotropic_iters(aniso_iters) ;
			curvatures.set_anisotropic_factor(aniso_factor) ;

			if(K1_attribute_name.length() != 0) {
				curvatures.set_Kmin_attribute(K1_attribute_name) ;
			} 

			if(K2_attribute_name.length() != 0) {
				curvatures.set_Kmax_attribute(K2_attribute_name) ;
			} 

			if(N_attribute_name.length() != 0) {
				curvatures.set_N_attribute(N_attribute_name) ;
			} 


			if(k1_attribute_name.length() != 0) {
				curvatures.set_kmin_attribute(k1_attribute_name) ;
			} 

			if(k2_attribute_name.length() != 0) {
				curvatures.set_kmax_attribute(k2_attribute_name) ;
			} 

			if(n_attribute_name.length() != 0) {
				curvatures.set_n_attribute(n_attribute_name) ;
			} 

			curvatures.compute_curvature_tensor() ;

			if(
				copy_n_to_normal && 
				MapVertexProperty<vec3>::is_defined(
				map_, N_attribute_name
				) 
				) {
					MapVertexProperty<vec3> v_normal(
						map_, N_attribute_name
						) ;
					MapTexVertexNormal t_normal(map_) ;
					FOR_EACH_HALFEDGE(Map, map_, it) {
						if(
							dot(Geom::vertex_normal(it->vertex()), v_normal[it->vertex()]) > 0
							) {
								t_normal[it->tex_vertex()] = v_normal[it->vertex()] ;
						} else {
							t_normal[it->tex_vertex()] = -1.0 * v_normal[it->vertex()] ;
						}
						t_normal[it->tex_vertex()] = normalize(t_normal[it->tex_vertex()]) ;
					}
			}

			if(compute_hk) {
				MapVertexProperty<double> k1(map_, k1_attribute_name) ;
				MapVertexProperty<double> k2(map_, k2_attribute_name) ;
				MapVertexProperty<double> h(map_, "h") ;
				MapVertexProperty<double> k(map_, "k") ;
				FOR_EACH_VERTEX(Map, map_, it) {
					h[it] = 0.5 * (k1[it] + k2[it]) ;
					k[it] = k1[it] * k2[it] ;
				}
			}

			MapVertexProperty<double> k1(map_, k1_attribute_name) ;
			MapVertexProperty<double> k2(map_, k2_attribute_name) ;
			MapVertexProperty<double> density(map_, "vertex_density") ;
			//			density_range_.first = 1e30 ; density_range_.second = -1e30 ;
			double average = 0 ;
			FOR_EACH_VERTEX(Map, map_, it) {
				density[it] = gx_max(fabs(k1[it]), fabs(k2[it])) ;/*+1e-5*/ //1.0/(fabs(k1[v])+1e-2) ;
				//density[it] = gx_max(density[it], 1e-3) ;
				//density_range_.first = gx_min(density[it], density_range_.first) ;
				//density_range_.second = gx_max(density[it], density_range_.second) ;
				average += density[it] ;
			}
			average /= map_->nb_vertices() ;
			FOR_EACH_VERTEX(Map, map_, it) {
				density[it] = gx_max(average*1e-3, density[it]) ;
				density[it] = gx_min(density[it], average*1e3) ;
			}


			//for(unsigned int i=0; i<boundary_->nb_vertices(); ++i) {
			//	int orgv = boundary_->vertex_index(i) ;
			//	boundary_->vertex(i).w = density[verts_[orgv]] ;
			//}

			//average/=map_->nb_vertices() ;

			//std::cerr << "Vertex density range = [" 
			//	<< density_range_.first << ", " << density_range_.second << "], average density " << average << std::endl ;

			//// clamp density values
			//density_range_.first = 1e30 ; density_range_.second = -1e30 ;
			//FOR_EACH_VERTEX(Map, map_, it) {
			//	if(density[it]<0.1*average)
			//		density[it] = 0.1*average ;
			//	if(density[it]>10*average)
			//		density[it] = 10.0*average ;
			//	density_range_.first = gx_min(density[it], density_range_.first) ;
			//	density_range_.second = gx_max(density[it], density_range_.second) ;
			//}
			//std::cerr << "Clamped Vertex density range = [" 
			//	<< density_range_.first << ", " << density_range_.second << "], average density " << average << std::endl ;
	}

	void DelaunayCVT::compute_lfs() {
		if(lfs_!=nil) delete lfs_ ;
		std::vector<vec3> vertices(boundary_->original_vertices());
		for(unsigned int i=0; i<boundary_->nb_facets(); i++) {
			vec3 g = boundary_->facet_center(i) ;
			vertices.push_back(g) ;
		}
		lfs_ = new MedialAxis(vertices);
		MapVertexProperty<double> density(map_, "vertex_density") ;
		FOR_EACH_VERTEX(Map, map_, v) {
			density[v] = 1.0/lfs_->squared_lfs(v->point());
		}
		for(unsigned int i=0; i<boundary_->nb_vertices(); i++) {
			double squared_dis = lfs_->squared_lfs(boundary_->vertex(i));
			//    boundary_->vertex(i).w = 1.0 / squared_dis ;
			boundary_->vertex(i).w = 1.0/(squared_dis) ; //1.0 / squared_dis ;
		}
	}

	void DelaunayCVT::smooth_density(Map* surface, int nb_iter) {
		MapVertexProperty<double> vertex_density ;
		MapVertexProperty<double> new_density ;
		vertex_density.bind_if_defined(surface, "vertex_density") ;
		new_density.bind(surface, "new_density") ;

		for(int i=0; i<nb_iter; ++i) {
			std::vector<double> debug_old, debug_new ;

			FOR_EACH_VERTEX(Map, surface, v) {
				Map::Halfedge* h = v->halfedge() ;
				new_density[v] = 0 ;
				double total_w = 0 ;
				do {
					double w = 1.0/distance(v->point(), h->opposite()->vertex()->point()) ;
					total_w += w ;
					new_density[v] += w * vertex_density[h->opposite()->vertex()] ;
					h = h->next_around_vertex() ;
				} while(h!=v->halfedge()) ;	
				new_density[v] /= total_w ;
				new_density[v] = 0.5*(new_density[v]+vertex_density[v]) ;
				debug_old.push_back(vertex_density[v]) ;
				debug_new.push_back(new_density[v]) ;
			}

			FOR_EACH_VERTEX(Map, surface, v) {
				vertex_density[v] = new_density[v] ;
			}
		}



		new_density.unbind() ;
	}

	void DelaunayCVT::save_chart(const std::string& filename) {
		std::ofstream out(filename.c_str()) ;

		for(unsigned int i=0; i<charts_.size(); ++i) {
			for(unsigned int j=0; j<charts_[i].size(); ++j) {
				fchart_[charts_[i][j]] = i ;
			}
		}
		FOR_EACH_FACET(Map, map_, f) {
			out << fchart_[f] << std::endl ;
		}

		out.close() ;
	}

	void DelaunayCVT::load_chart(const std::string& filename) {
		std::ifstream in(filename.c_str()) ;

		if(!fchart_.is_defined(map_, "chart_id")) {
			fchart_.bind(map_, "chart_id") ;
		}

		int fid = 0 ;
		while(!in.eof()) {
			int cid ;
			in >> cid ;

			fchart_[facets_[fid]] = cid ;
			fid ++ ;
		} ;

		in.close() ;
	}

	///dxy add: 
	void DelaunayCVT::load_field(const std::string& filename) {
		std::cout << "Loading field..." << std::endl;

		std::ifstream input(filename.c_str()) ;
		if (!input)
		{
			std::cerr << "could not open file." << std::endl;
			return;
		}
		LineInputStream in(input);

		std::vector<vec3> &field = boundary_->vertex_field();
		std::vector<vec3> &vnorm = boundary_->vertex_normal();

		while (!in.eof())
		{
			in.get_line();
			std::string keyword;

			in >> keyword;

			if (keyword == "vf")
			{
				vec3 f;
				in >> f;
				field.push_back(f);
			}
			else if (keyword == "vn")
			{
				vec3 n;
				in >> n;
				vnorm.push_back(n);
			}
		}

		std::cout << "field size  = " << field.size() << std::endl;
		std::cout << "vnormal size = " << vnorm.size() << std::endl;
		std::cout << "vertex size = " << boundary_->original_vertices().size() << std::endl;
		std::cout << "nb_vertices = " << boundary_->nb_vertices() << std::endl;

		input.close();

		field_loaded_ = GL_TRUE;
		std::cout << "Loading field done." << std::endl;
	}

	void DelaunayCVT::calcFacetNormalAndDirection() {
		if (!field_loaded_) return;

		int RoSy = field_rot_symmetry_;

		//facet normal
		auto &fn = boundary_->facet_normal();
		fn.clear();
		for (unsigned int f=0; f<boundary_->nb_facets(); ++f)
		{
			fn.push_back(normalize(boundary_->facet_normal(f)));
		}

		//facet field
		auto &ff = boundary_->facet_field();
		ff.clear();
		for (unsigned int f=0; f<boundary_->nb_facets(); ++f)
		{
			int v0 = boundary_->facet_begin(f);
			int nvert = boundary_->facet_end(f) - v0;
			if (nvert==0)
			{
				ff.push_back(vec3(0, 0, 0));
				std::cout << f << ": degenerate facet!" << std::endl;
				continue;
			}

			//
			std::vector<std::vector<vec3>> proj_field(nvert);
			for (int i=0; i<nvert; ++i)
			{
				vec3 vf = boundary_->vertex_field(boundary_->vertex_index(v0+i));
				vec3 axie = boundary_->vertex_normal(boundary_->vertex_index(v0+i));
				for (int rot=0; rot<RoSy; ++rot)
				{
					vec3 rotvf = vec_rotate(vf, axie, rot*2*M_PI/RoSy);
					proj_field[i].push_back(rotvf - dot(rotvf, fn[f])*fn[f]);
				}
			}

			//
			vec3  maxvec(0, 0, 0);
			double maxlen2 = 0;

			std::vector<std::pair<int, vec3>> S;
			S.push_back(std::make_pair(0, proj_field[0][0]));
			while (!S.empty())
			{
				if (S.back().first >= RoSy)
				{
					S.pop_back();
					if (!S.empty())
					{
						S.back().first += 1;
					}
				}
				else {
					if (S.size() == nvert)
					{
						double len2 = length2(S.back().second);
						if (len2 > maxlen2)
						{
							maxlen2 = len2;
							maxvec = S.back().second;
						}
						int last_rot = S.back().first;
						if (last_rot + 1 >= RoSy)
						{
							S.back().first = RoSy;
						}
						else {
							S.pop_back();
							S.push_back(std::make_pair(last_rot+1, S.back().second + proj_field[S.size()][last_rot+1]));
						}
					}
					else {
						S.push_back(std::make_pair(0, S.back().second + proj_field[S.size()][0]));
					}
				}
			}

			ff.push_back(normalize(maxvec));
		}
	}
	///

}
