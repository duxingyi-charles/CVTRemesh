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

#include <Geex/graphics/geexapp.h>
#include <Geex/basics/file_system.h>
#include <Geex/basics/processor.h>
#include <ctime>
#include <cstdlib>
#include <AntTweakBar.h>
#include "cvt.h"
#include <stdio.h>

namespace Geex {

	class GPURVDApp : public GeexApp {
	public:
		GPURVDApp(int argc, char** argv) : GeexApp(argc, argv) { 
			hdr_ = false ;
			boundary_filename_ = get_file_arg("obj") ;
			if(boundary_filename_.length() > 0) {
				if(!Geex::FileSystem::is_file(boundary_filename_)) {
					boundary_filename_ = Geex::FileSystem::get_project_root() + "/data/surfaces/" + boundary_filename_ ;
				}
			}
			else {
				std::cerr << "invalid boundary filename" << std::endl ;
				exit(0) ;
			}
			points_filename_ = get_file_arg("pts") ;
			if(points_filename_.length() > 0) {
				if(!Geex::FileSystem::is_file(points_filename_)) {
					points_filename_ = Geex::FileSystem::get_project_root() + "/data/surfaces/" + points_filename_ ;
				}
			}

			nb_points_ = 100 ;
			get_arg("nb_pts", nb_points_) ;
			nb_iters_ = 30 ;
			get_arg("nb_iter", nb_iters_) ;
			//dxy add
			multithread_ = GL_FALSE;
			get_arg("mthread", multithread_);
			//
		}

		~GPURVDApp() {
			std::cout << "deconstruct cvt app. " << std::endl ;
		}

		GPURVD* cvt() { return static_cast<GPURVD*>(scene()) ; }

		virtual void init_scene() {
			scene_ = new GPURVD ;
			//            std::cerr << "Non convex = " << (non_convex_ ? "true" : "false") << "(use +non_convex to set)" << std::endl ;

			if(boundary_filename_.length() > 0) {
				///dxy change
				/*cvt()->load_map(boundary_filename_) ;
				cvt()->load_boundary(boundary_filename_) ;*/
				//
				cvt()->load_map(boundary_filename_, false);
				double area = cvt()->calc_map_area();
				cvt()->model_normalize_scale() = sqrt(nb_points_ / area); //scale to make average area per seed is 1 			
				cvt()->model_normalize_map();
				cvt()->load_boundary(boundary_filename_, true) ;
				///
				std::cout << "boundary loaded..." << std::endl ;
				std::string feature_filename = boundary_filename_ ;
				feature_filename.replace(feature_filename.size()-4, 4, ".fea") ;
				if(cvt()->load_feature(feature_filename)) {
					cvt()->mark_feature() ;  
				} else { //dxy test
					if (cvt()->lock_corners() || cvt()->lock_feature()) {
						cvt()->extract_features(cvt()->map());
					}
				}

				//dxy add: auto experiment
				std::string experiment_filename = boundary_filename_;
				experiment_filename.replace(experiment_filename.size()-4, 4, ".exp");
				cvt()->load_experiment(experiment_filename);
				//
				std::cout << "start insert samplings..." << std::endl ;
				cvt()->init_samples_w(nb_points_) ;
				cvt()->compute_rvd() ;

				//dxy add
				if(multithread_) {
					std::cerr << "Initializing multithread" << std::endl ;
					cvt()->init_multithread() ; 
					cvt()->use_multithread() = GL_TRUE ;
				}
				///
			}
		}

		void reset() {
			if (cvt()->use_face_center_init())
			{
				cvt()->init_samples_barycentric(nb_points_);
			}
			else {
				cvt()->init_samples_w(nb_points_) ;
			}
			cvt()->compute_rvd() ;
		}

		void save_points() {
			std::string filename = boundary_filename_ ;
			filename.replace(filename.size()-4, 4, ".pts") ; ///dxy change: .txt to .pts
			cvt()->save_samples(filename) ;
		}

		///dxy add
		void save_denormalized_points() {
			std::string filename = boundary_filename_ ;
			filename.replace(filename.size()-4, 4, ".pts") ;
			cvt()->save_denormalized_samples(filename) ;
		}
		///

		void save_density() {
			std::string filename = boundary_filename_ ;
			filename.replace(filename.size()-4, 4, ".lfs") ;
			cvt()->save_density(filename) ;
		}

		void load_density() {
			std::string filename = boundary_filename_ ;
			filename.replace(filename.size()-4, 4, ".lfs") ;
			if(FileSystem::is_file(filename)) {
				cvt()->load_density(filename) ;
			}
		}

		void load_points() {
			std::string filename = boundary_filename_ ;
			filename.replace(filename.size()-4, 4, ".pts") ; ///dxy change: .txt to .pts
			//cvt()->load_samples(filename) ;
			cvt()->load_samples(filename, true);
			cvt()->compute_rvd() ;
		}

		//dxy add
		void load_points(const std::string &filename) {
			cvt()->load_samples(filename, true);
			cvt()->compute_rvd();
		}
		//

		void save_primal() {
			std::string filename = boundary_filename_ ;
			filename.insert(filename.size()-4, "_rvd") ;
			//cvt()->save_primal_global(filename) ;
			cvt()->save_denormalized_primal(filename);

			//filename = boundary_filename_;
			//filename.insert(filename.size()-4,  "_org"); //origin
			//cvt()->save_denormalized_boundary(filename);

			///dxy remove
			//filename = boundary_filename_ ;
			//filename.insert(filename.size()-4, "_normalized") ;
			//cvt()->save_normalized_boundary(filename) ;
		}

		///dxy add
		void save_primal_map() {
			std::string filename = boundary_filename_ ;
			filename.insert(filename.size()-4, "_map") ;
			cvt()->save_primal_map(filename) ;
		}

		void save_primal_angle() {
			std::string filename = boundary_filename_;
			filename.insert(filename.size()-4, "_angle");
			filename.replace(filename.size()-4, 4, ".txt");
			cvt()->save_primal_angles(filename);
		}

		void save_primal_area() {
			std::string filename = boundary_filename_;
			filename.insert(filename.size()-4, "_area");
			filename.replace(filename.size()-4, 4, ".txt");
			cvt()->save_primal_areas(filename);
		}

		void save_primal_info() {
			std::string filename = boundary_filename_;
			filename.insert(filename.size()-4, "_info");
			filename.replace(filename.size()-4, 4, ".txt");
			cvt()->save_primal_info(filename);
		}

		void save_opt_info() {
			std::string filename = boundary_filename_;
			filename.replace(filename.size()-4, 4, ".log");
			cvt()->save_optimization_info(filename);
		}
		///

		void lloyd_global() {
			cvt()->lloyd_global(nb_iters_) ;
		}
		void lloyd_global_step() {
			cvt()->lloyd_global(1) ;
		}

		void save_chart() {
			std::string filename = boundary_filename_ ;
			filename.replace(filename.size()-4, 4, ".seg") ;
			cvt()->save_chart(filename) ;
		}
		void load_chart() {
			std::string filename = boundary_filename_ ;
			filename.replace(filename.size()-4, 4, ".seg") ;
			cvt()->load_chart(filename) ;
		}

		///dxy add
		void save_feature() {
			std::string filename = boundary_filename_;
			filename.replace(filename.size()-4, 4, ".fea");
			cvt()->save_feature(filename);
		}

		void load_field() {
			std::string filename = boundary_filename_;
			//filename.replace(filename.size()-4, 4, ".fid") ;
			cvt()->load_field(filename) ;

			//dxy test
			cvt()->calcFacetNormalAndDirection();
		}

		void split_long_edges() {
			cvt()->split_long_edges();
		}

		void collapse_short_edges() {
			cvt()->collapse_short_edges();
		}

		void update_area() {
			cvt()->update_area();
		}

		///

		//void lloyd_local() {
		//	cvt()->lloyd_local(nb_iters_) ;
		//}

		//void lloyd_local_step() {
		//	cvt()->lloyd_local(1) ;
		//}

		//void newton_lloyd_local() {
		//	cvt()->newton_lloyd_local(nb_iters_) ;
		//}
		void newton_lloyd_global() {
			if (cvt()->use_stretch_topo_optimize())
			{
				cvt()->stretch_topo_optimize(nb_iters_);  //dxy add
			}
			else {
				cvt()->newton_lloyd_global(nb_iters_) ;
			}
		}

		void smooth_density() {
			cvt()->smooth_density(cvt()->map(), 1) ;
		}

		//void update_lrvd_one_step() {
		//	cvt()->update_lrvd_step() ;
		//}

		//void split_non_manifold() {
		//	cvt()->split_non_manifold_primal_triangles() ;
		//}

		//void insert_obtuse_centers() {
		//	cvt()->insert_sample_obtuse() ;
		//}

		//dxy add
		void experiment() {
			const auto& experiments = cvt()->experiments();
			for(const auto &e : experiments) {
				std::cout << "Experiment No." << e.num << " begins..." << std::endl;
				cvt()->use_density() = e.use_density;
				if (e.use_quad)
				{
					cvt()->rot_symmetry() = 4;
					cvt()->use_special_quad() = false/*true*/;
				}
				else {
					cvt()->rot_symmetry() = 6;
				}
				// Begin Run
				if (e.pts >= 0) {
					std::string pts_filename = boundary_filename_;
					pts_filename.replace(pts_filename.size()-4, 4, ".pts");
					pts_filename.insert(pts_filename.size()-4, "_"+std::to_string(e.pts));
					cvt()->load_samples(pts_filename, true);
					cvt()->compute_rvd() ;
					//
					if (e.use_all_edge_split) {
						nb_points_ *= 4;
						reset();
						cvt()->load_samples(pts_filename, true);
						cvt()->compute_rvd() ;
						cvt()->split_all_edges();
					}
					else {
						if (fabs(e.mix_factor) < 1e-5) //cvt
						{
							std::cout << "Exp No." << e.num << ": CVT " << e.nb_iter << " iters." << std::endl;
							cvt()->use_cvt_energy() = true;
							cvt()->use_direction() = false;
							cvt()->newton_lloyd_global(e.nb_iter/*nb_iters_*/);
							//cvt()->update_area();

						}
						else if (e.mix_factor > 1e-5) {  //fcvt
							load_field();
							cvt()->use_cvt_energy() = true;
							cvt()->use_direction() = true;
							cvt()->energy_mixfactor() = e.mix_factor;
							cvt()->update_area();
							if (e.use_topo_opt) {
								std::cout << "Exp No." << e.num << ": FCVT with topo opt, " << e.nb_iter << " iters, mf=" << e.mix_factor << std::endl;
								int nfail = 0;
								for (int i=0; i<e.nb_topo_opt_iter; ++i) {
									cvt()->newton_lloyd_global(e.nb_iter);
									//topo optimize
									if (e.use_quad) { //match topo opt
										cvt()->use_edge_dir_match_topo_opt() = true;
										cvt()->newton_lloyd_global(1); //calc edge_dir_match
										cvt()->use_edge_dir_match_topo_opt() = false;
										if (!cvt()->quad_dominant_topo_optimize()) ++nfail;
										else nfail = 0;	
									} else { //traditional topo opt

									}
									//temp save rvd, pts
									std::string rvd_filename = boundary_filename_ ;
									rvd_filename.insert(rvd_filename.size()-4, "_rvd_"+std::to_string(e.num)) ;
									cvt()->save_denormalized_primal(rvd_filename);
									std::string pts_filename = boundary_filename_;
									pts_filename.replace(pts_filename.size()-4, 4, ".pts");
									pts_filename.insert(pts_filename.size()-4, "_"+std::to_string(e.num));
									cvt()->save_denormalized_samples(pts_filename);
									//termination check
									if (nfail > 2) break;
								} 
							} else {
								std::cout << "Exp No." << e.num << ": FCVT, " << e.nb_iter << " iters, mf=" << e.mix_factor << std::endl;
								cvt()->newton_lloyd_global(e.nb_iter/*nb_iters_*/);
							}
						}
					}
					//save opt info
					std::string info_filename = boundary_filename_;
					info_filename.replace(info_filename.size()-4, 4, ".log");
					info_filename.insert(info_filename.size()-4, "_"+std::to_string(e.num));
					cvt()->save_optimization_info(info_filename);
				}
				else { //init
					reset();
				}
				//save rvd, pts
				std::string rvd_filename = boundary_filename_ ;
				rvd_filename.insert(rvd_filename.size()-4, "_rvd_"+std::to_string(e.num)) ;
				cvt()->save_denormalized_primal(rvd_filename);
				std::string pts_filename = boundary_filename_;
				pts_filename.replace(pts_filename.size()-4, 4, ".pts");
				pts_filename.insert(pts_filename.size()-4, "_"+std::to_string(e.num));
				cvt()->save_denormalized_samples(pts_filename);
			}
		}

		void pts_experiment() {
			for (int it=0; it<cvt()->experiment_iter(); ++it) 
			{
				reset();
				std::string filename = boundary_filename_;
				std::string tail = "_";
				tail.append(std::to_string(it));
				tail.append("_1.pts");
				filename.replace(filename.size()-4, 4, tail);
				cvt()->save_denormalized_samples(filename);   //Note:
				newton_lloyd_global();
				filename.replace(filename.size()-5, 5, "2.pts");
				cvt()->save_denormalized_samples(filename);
			}
		}



		void qd_experiment() {
			//CVT
			cvt()->use_cvt_energy() = true;
			cvt()->use_direction() = false;
			cvt()->newton_lloyd_global(100);
			//init fcvt
			load_field();
			cvt()->use_direction() = true;
			cvt()->use_edge_dir_match() = false;
			cvt()->update_area();
			cvt()->newton_lloyd_global(100);
			//Exp iteration
			for (int i=0; i<cvt()->experiment_iter(); ++i)
			{
				std::cout << "Experiment " << i << ":" << std::endl;
				//match_fcvt
				cvt()->use_edge_dir_match() = true;
				cvt()->update_area();
				cvt()->newton_lloyd_global(10);

				//save result
				save_denormalized_points();
				save_primal();
				//topo opt
				cvt()->use_edge_dir_match() = false;
				cvt()->quad_dominant_topo_optimize();
				//fcvt
				cvt()->update_area();
				cvt()->newton_lloyd_global(20);
			}
		}
		///dxy 


		virtual	void init_gui() {
			GeexApp::init_gui() ;

			// New-style GUI =====================================================================================
			TwBar* graphics_bar = TwNewBar("Graphics") ;
			TwDefine(" Graphics position='16 10' size='200 850'") ;    
			TwAddVarRW(graphics_bar, "shiny", TW_TYPE_BOOL8, &cvt()->shiny(), "") ;
			TwAddVarRW(graphics_bar, "colorize", TW_TYPE_BOOL8, &cvt()->colorize(), "") ;
			TwAddVarRW(graphics_bar, "Vertices", TW_TYPE_FLOAT, &cvt()->vertices_size(), "min=0 max=1 step=0.01") ;
			TwAddVarRW(graphics_bar, "RVDVerts", TW_TYPE_BOOL8, &cvt()->show_rvd_vertices(), "") ;
			TwAddVarRW(graphics_bar, "Domain mesh", TW_TYPE_BOOL8, &cvt()->show_domain_mesh(), "") ;
			TwAddVarRW(graphics_bar, "Domain", TW_TYPE_BOOL8, &cvt()->show_domain(), "") ;
			TwAddVarRW(graphics_bar, "Feature", TW_TYPE_BOOL8, &cvt()->show_feature(), "") ;
			//	TwAddVarRW(graphics_bar, "ActiveTri", TW_TYPE_BOOL8, &cvt()->show_active_triangles(), "") ;
			//	TwAddVarRW(graphics_bar, "Cluster", TW_TYPE_BOOL8, &cvt()->show_cluster(), "") ;

			TwAddSeparator(graphics_bar, "RVD", "");
			TwAddVarRW(graphics_bar, "Primal", TW_TYPE_BOOL8, &cvt()->show_primal(), "") ;
			TwAddVarRW(graphics_bar, "Obtuse", TW_TYPE_BOOL8, &cvt()->show_obtuse(), "") ;
			TwAddVarRW(graphics_bar, "Features", TW_TYPE_BOOL8, &cvt()->show_primal_features(), "");
			TwAddVarRW(graphics_bar, "MinAngle", TW_TYPE_DOUBLE, &cvt()->angle_min(), "min=0 max=60 step=0.01") ;
			TwAddVarRW(graphics_bar, "MaxAngle", TW_TYPE_DOUBLE, &cvt()->angle_max(), "min=60 max=180 step=0.01") ;
			TwEnumVal surface_mode_def[] = { {None, "None"}, {Plain, "Plain"}, {Chart, "Chart"}, {Cells, "Cells"}, {CellType, "CellType"}, {Density, "Density"},  {Energy, "Energy"}, {GeoDist, "GeoDist"} } ;
			TwType    tw_surface_mode = TwDefineEnum("SuraceMode", surface_mode_def, 8) ;
			TwAddVarRW(graphics_bar, "Surface Mode", tw_surface_mode, &cvt()->surface_mode(), "") ;
			//		TwAddVarRW(graphics_bar, "LocalRVD", TW_TYPE_BOOL8, &cvt()->show_approx_rvd(), "") ;
			//		TwAddVarRW(graphics_bar, "ClippedRVD", TW_TYPE_BOOL8, &cvt()->show_clipped_rvd(), "") ;
			TwAddVarRW(graphics_bar, "ExactRVD", TW_TYPE_BOOL8, &cvt()->show_rvd(), "") ;
			TwAddVarRW(graphics_bar, "ExactRVD mesh", TW_TYPE_BOOL8, &cvt()->show_rvd_mesh(), "") ;

			TwAddSeparator(graphics_bar, "Optimization", "");
			TwAddVarRW(graphics_bar, "Iteration", TW_TYPE_UINT32, &nb_iters_, "") ;
			TwAddVarRW(graphics_bar, "Density", TW_TYPE_BOOL8, &cvt()->use_density(), "") ;
			if(multithread_) {
				TwAddVarRW(graphics_bar, "Multi-thread", TW_TYPE_BOOL8, &cvt()->use_multithread(), "") ;
			}
			TwAddVarRW(graphics_bar, "Face Center Init", TW_TYPE_BOOL8, &cvt()->use_face_center_init(), "");
			TwAddVarRW(graphics_bar, "Lock Corners", TW_TYPE_BOOL8, &cvt()->lock_corners(), "");
			//TwAddVarRW(graphics_bar, "Use Feature", TW_TYPE_BOOL8, &cvt()->use_feature(), "") ;
			TwAddVarRW(graphics_bar, "Lock Feature", TW_TYPE_BOOL8, &cvt()->lock_feature(), "") ;
			TwAddVarRW(graphics_bar, "NormalAniso", TW_TYPE_DOUBLE, &cvt()->normal_aniso(), "min=0 max=100 step=0.01") ;
			TwAddVarRW(graphics_bar, "LockRing", TW_TYPE_UINT32, &cvt()->nb_lock_rings(), "min=0 max=10 step=1") ;
			toggle_skybox_CB() ;

			///dxy add
			TwAddSeparator(graphics_bar, "Direction", "");
			TwAddVarRW(graphics_bar, "CVT Energy", TW_TYPE_BOOL8, &cvt()->use_cvt_energy(), "");
			TwAddVarRW(graphics_bar, "Direction Field", TW_TYPE_BOOL8, &cvt()->use_direction(), "");
			TwAddVarRW(graphics_bar, "Auto Save pts", TW_TYPE_BOOL8, &cvt()->use_auto_save(), "");
			TwAddVarRW(graphics_bar, "Use True Grad", TW_TYPE_BOOL8, &cvt()->use_true_gradient(), "");
			TwAddVarRW(graphics_bar, "Use Match", TW_TYPE_BOOL8, &cvt()->use_edge_dir_match(), "");
			TwAddVarRW(graphics_bar, "Use Face Field", TW_TYPE_BOOL8, &cvt()->use_facet_field(), "");
			TwAddVarRW(graphics_bar, "Field RoSy", TW_TYPE_UINT32, &cvt()->rot_symmetry(), "min=1 max=100 step=1");
			TwAddVarRW(graphics_bar, "Special 4-RoSy", TW_TYPE_BOOL8, &cvt()->use_special_quad(), "");
			TwAddVarRW(graphics_bar, "Adaptive Mix", TW_TYPE_BOOL8, &cvt()->use_adaptive_mix_factor(), "");
			TwAddVarRW(graphics_bar, "Mix Factor", TW_TYPE_DOUBLE, &cvt()->energy_mixfactor(), "min=0.0 max=1000.0 step=0.0000001") ;
			TwAddVarRW(graphics_bar, "E1 para a", TW_TYPE_DOUBLE, &cvt()->dir_para_a(), "min=1.0 max=5.0 step=0.1");
			TwAddVarRW(graphics_bar, "E1 para quad", TW_TYPE_DOUBLE, &cvt()->dir_para_quad(), "min=0.0 max=0.8 step=0.01");
			//
			// 			TwEnumVal direction_energy_mode_def[] = { {Angle_Only, "Angle"}, {Angle_Edge, "Angle-Edge"}, {Angle_Edge2, "Angle_Edge2"} };
			// 			TwType    tw_direction_energy_mode = TwDefineEnum("EnergyMode", direction_energy_mode_def, 3) ;
			// 			TwAddVarRW(graphics_bar, "Energy Mode", tw_direction_energy_mode, &cvt()->direction_energy_mode(), "") ;

			TwEnumVal direction_edge_weight_mode_def[] = { {Dual_Length, "Dual length"}, {Lloyd_Energy, "Lloyd energy"} };
			TwType	  tw_direction_edge_weight_mode = TwDefineEnum("WeightMode", direction_edge_weight_mode_def, 2);
			TwAddVarRW(graphics_bar, "Weight Mode", tw_direction_edge_weight_mode, &cvt()->direction_edge_weight_mode(), "") ;

			TwAddVarRW(graphics_bar, "New Interpolation", TW_TYPE_BOOL8, &cvt()->use_new_interpolation(), "");
			// 			TwEnumVal interpolation_mode_def[] = { {Barycentric, "Barycentric"}, {Cotangent, "Cot"} };
			// 			TwType    tw_interpolation_mode = TwDefineEnum("InterpolateMode", interpolation_mode_def, 2) ;
			// 			TwAddVarRW(graphics_bar, "Interpolate Mode", tw_interpolation_mode, &cvt()->interpolationMode(), "");
			//
			TwAddSeparator(graphics_bar, "Visualization", "");
			TwAddVarRW(graphics_bar, "Topo Operation", TW_TYPE_BOOL8, &cvt()->show_valid_topo_operation(), "");
			TwAddVarRW(graphics_bar, "deltaR", TW_TYPE_INT32, &cvt()->show_topo_operation_deltaR(), "min=-10 max=50 step=1");
			TwAddVarRW(graphics_bar, "Topo_EF", TW_TYPE_BOOL8, &cvt()->show_edge_flip(), "");
			TwAddVarRW(graphics_bar, "Topo_EC", TW_TYPE_BOOL8, &cvt()->show_edge_collapse(), "");
			TwAddVarRW(graphics_bar, "Topo_VS", TW_TYPE_BOOL8, &cvt()->show_vertex_split(), "");
			TwAddVarRW(graphics_bar, "Sharp Edge", TW_TYPE_BOOL8, &cvt()->show_sharp_edges(), "");
			TwAddVarRW(graphics_bar, "Singularity", TW_TYPE_BOOL8, &cvt()->show_singular_rvd_vertices(), "");
			TwAddVarRW(graphics_bar, "Project Tri", TW_TYPE_BOOL8, &cvt()->show_seed_project_triangle(), "");
			TwAddVarRW(graphics_bar, "Domain Field", TW_TYPE_BOOL8, &cvt()->show_domain_field(), "");
			TwAddVarRW(graphics_bar, "Face Field", TW_TYPE_BOOL8, &cvt()->show_domain_facet_field(), "");
			TwAddVarRW(graphics_bar, "Seed Field", TW_TYPE_BOOL8, &cvt()->show_seed_field(), "");
			TwAddVarRW(graphics_bar, "Match", TW_TYPE_BOOL8, &cvt()->show_se
				ed_match_field(), "");
			TwAddVarRW(graphics_bar, "Max Order", TW_TYPE_UINT32, &cvt()->max_match_dir_order(), "min=0 max=3 step=1");
			TwAddVarRW(graphics_bar, "Max Num", TW_TYPE_INT32, &cvt()->max_match_number(), "min=0 max=20 step=1");
			TwAddVarRW(graphics_bar, "QD Operation", TW_TYPE_BOOL8, &cvt()->show_QD_topo_operations(), "");
			TwAddVarRW(graphics_bar, "Normal", TW_TYPE_BOOL8, &cvt()->show_normal(), "");
			TwAddVarRW(graphics_bar, "Short Edges", TW_TYPE_BOOL8, &cvt()->show_short_edges(), "");
			TwAddVarRW(graphics_bar, "Long Edges", TW_TYPE_BOOL8, &cvt()->show_long_edges(), "");
			TwAddVarRW(graphics_bar, "Cell Energy", TW_TYPE_BOOL8, &cvt()->show_cell_energy(), "");
			TwAddVarRW(graphics_bar, "Truncate Val", TW_TYPE_DOUBLE, &cvt()->cell_energy_truncate_value(), "min=0 max=1 step=0.0001");
			TwEnumVal draw_cell_energy_mode_def[] = { {CVT_True, "true cvt"}, {CVT_Approx, "approx cvt"}, {CVT_True_Minus_Approx, "diff"}, {DIR, "dir"}, {Cell_Area, "area"}};
			TwType	  tw_draw_cell_energy_mode_def = TwDefineEnum("CellEnergyMode", draw_cell_energy_mode_def, 5);
			TwAddVarRW(graphics_bar, "show energy", tw_draw_cell_energy_mode_def, &cvt()->draw_energy_mode(), "") ;
			TwAddVarRW(graphics_bar, "lloyd grad", TW_TYPE_BOOL8, &cvt()->show_lloyd_grad(), "");
			TwAddVarRW(graphics_bar, "direction grad", TW_TYPE_BOOL8, &cvt()->show_direction_grad(), "");
			TwAddVarRW(graphics_bar, "grad mag", TW_TYPE_FLOAT, &cvt()->grad_magnification(), "min=0.0 max=100.0 step=0.1");

			//
			TwBar* optimize_bar = TwNewBar("Optimize") ;
			TwDefine(" Optimize position='240 10' size='200 220'") ;

			TwAddVarRW(optimize_bar, "stretch long", TW_TYPE_BOOL8, &cvt()->show_stretched_long_edges(), "");
			TwAddVarRW(optimize_bar, "stretch short", TW_TYPE_BOOL8, &cvt()->show_stretched_short_edges(), "");
			TwAddVarRW(optimize_bar, "show rate", TW_TYPE_FLOAT, &cvt()->show_stretched_edges_rate(), "min=0.0 max=1.0 step=0.01");
			TwAddSeparator(optimize_bar, "Stretch Opt", "");
			TwAddVarRW(optimize_bar, "stretch opt", TW_TYPE_BOOL8, &cvt()->use_stretch_topo_optimize(), "");
			TwAddVarRW(optimize_bar, "select rate", TW_TYPE_FLOAT, &cvt()->stretch_optimize_select_rate(), "min=0.0 max=0.5 step=0.01");
			TwAddVarRW(optimize_bar, "smooth iter", TW_TYPE_INT32, &cvt()->stretch_optimize_nb_smooth_iter(), "min=0 max=50 step=1");
			TwAddSeparator(optimize_bar, "Traditional Topo Opt", "");
			TwAddVarRW(optimize_bar, "TopoOpt iter", TW_TYPE_UINT32, &cvt()->map_topo_optimization_nb_iter(), "min=0 max=100000 step=1");
			TwAddVarRW(optimize_bar, "deltaR", TW_TYPE_INT32, &cvt()->map_topo_optimization_deltaR(), "min=-10 max=50 step=1");
			TwAddVarRW(optimize_bar, "dihedral", TW_TYPE_DOUBLE, &cvt()->dihedral_angle_threshold(), "min=0 max=180 step=1");
			//TwAddSeparator(optimize_bar, "Match Topo Opt", "");
			//TwAddSeparator(graphics_bar, "Match Topo Opt", TW_TYPE_BOOL8, &cvt()->use_edge_dir_match_topo_opt(), "");

			TwAddSeparator(optimize_bar, "Experiments", "");
			TwAddVarRW(optimize_bar, "Experiment iter", TW_TYPE_INT32, &cvt()->experiment_iter(), "min=0 max=10000 step=1");
			///


			glut_viewer_add_toggle('B', glut_viewer_is_enabled_ptr(GLUT_VIEWER_BACKGROUND), "switch Color/BW") ;
			glut_viewer_add_toggle('T', &viewer_properties_->visible(), "viewer properties") ;
			viewer_properties_->visible() = false ;

		}

	private:
		std::string boundary_filename_ ;
		std::string points_filename_ ; // input points for Delaunay Triangulation
		int nb_points_ ;
		int nb_iters_ ;
		//dxy add
		GLboolean multithread_;
		//
	} ;
}

Geex::GPURVDApp* GPURVD_app() { return static_cast<Geex::GPURVDApp*>(Geex::GeexApp::instance()) ; }

void reset() {
	GPURVD_app()->reset() ;
}

void load() {
	GPURVD_app()->load_points() ;
	//GPURVD_app()->load_density() ; //dxy comment
}

void save() {
	///dxy change: comment save_density()
	//GPURVD_app()->save_points() ;
	GPURVD_app()->save_denormalized_points();
	GPURVD_app()->save_primal() ;
	GPURVD_app()->save_density() ;

	///dxy add
	//GPURVD_app()->save_primal_angle();
	//GPURVD_app()->save_primal_area();

	//GPURVD_app()->save_primal_info(); //temp 
}

//dxy test
void save_Map() {
	GPURVD_app()->save_primal_map();
}
//

void optimize_global() {
	GPURVD_app()->lloyd_global() ;
}
void optimize_global_step() {
	GPURVD_app()->lloyd_global_step() ;
}

void optimize_newton_global() {
	GPURVD_app()->cvt()->update_area();
	GPURVD_app()->newton_lloyd_global() ;
}

void smooth_density() {
	GPURVD_app()->smooth_density() ;
}

void compute_lfs() {
	GPURVD_app()->cvt()->compute_lfs() ;
}
void compute_curv() {
	GPURVD_app()->cvt()->compute_curvature_tensor() ;
}

void extract_features() {
	GPURVD_app()->cvt()->extract_features(
		GPURVD_app()->cvt()->map()
		) ;
	//GPURVD_app()->cvt()->extract_chart() ;
	//GPURVD_app()->save_chart() ;
	GPURVD_app()->save_feature();
}

///dxy add
void batch_snapshot() {
	std::string dir = "E:/Research/fcvt_svn/CVTRemesh/data/surfaces/fcvt/svn/simple/ellipsoid/180401-eg fastforward/eg18";
	char filename[1024] ;
	int N = 0, nSol = 0 ;
	std::cout << "Please specify N and nSol:" << std::endl;
	std::cout << "N:";
	std::cin >> N;
	std::cout << "nSol:";
	std::cin >> nSol;
	for(int i=0; i<nSol; ++i) {
		sprintf(filename, "%s/%d/%d.pts", dir.c_str(), N, i) ;
		std::cout << "load points from " << filename << std::endl;
		GPURVD_app()->load_points(filename);
 		glut_viewer_redraw() ;
		sprintf(filename, "%s/%d/%d.png", dir.c_str(), N, i) ;
		GPURVD_app()->cvt()->save_snapshot(filename) ;
		//
		if (GPURVD_app()->cvt()->use_auto_save()) {
			sprintf(filename, "%s/%d/%d.obj", dir.c_str(), N, i) ;
			GPURVD_app()->cvt()->save_denormalized_primal(filename) ;
		}
	}
}

//dxy add
void load_field() {
	GPURVD_app()->load_field();
}

void edge_split() {
	GPURVD_app()->split_long_edges();
}

void all_edge_split() {
	GPURVD_app()->cvt()->split_all_edges();
}

void edge_collapse() {
	GPURVD_app()->collapse_short_edges();
}

void update_area() {
	GPURVD_app()->update_area();
}

void experiment() {
	GPURVD_app()->experiment();
}

void qd_experiment() {
	GPURVD_app()->qd_experiment();
}

void show_primal_info() {
	GPURVD_app()->cvt()->print_primal_info(true, true, true);
}

void quad_dominant_topo_optimize() {
	GPURVD_app()->cvt()->quad_dominant_topo_optimize();
}

void map_topo_optimize() {
	double d = cos(GPURVD_app()->cvt()->dihedral_angle_threshold()*M_PI/180.0);
	GPURVD_app()->cvt()->map_topo_optimization(GPURVD_app()->cvt()->map_topo_optimization_nb_iter(),
		-GPURVD_app()->cvt()->map_topo_optimization_deltaR(),
		d);
}

void update_boundary() {
	GPURVD_app()->cvt()->update_vertices_from_delaunay_map();
}

void save_opt_info() {
	GPURVD_app()->save_opt_info();
}
//test
void test() {
	//GPURVD_app()->cvt()->check_map();
	GPURVD_app()->pts_experiment();
}
///




#if (defined(WIN32) && defined(_DEBUG))
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif

#include <iostream>
#include <fstream>

int main(int argc, char** argv) {

#if (defined(WIN32) && defined(_DEBUG))
	_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
#endif

	srand((unsigned int)time(0));

	Geex::initialize() ;
	glut_viewer_add_key_func('Z', reset, "reset") ;
	glut_viewer_add_key_func('m', optimize_newton_global, "optimize_g") ;
	glut_viewer_add_key_func('k', optimize_global, "optimize_global") ;
	glut_viewer_add_key_func('K', optimize_global_step, "optimize_global_step") ;
	glut_viewer_add_key_func('d', compute_lfs, "compute_lfs") ;
	glut_viewer_add_key_func('D', compute_curv, "compute_curv") ;
	glut_viewer_add_key_func('s', save, "save") ;
	glut_viewer_add_key_func('o', load, "load") ;
	glut_viewer_add_key_func('c', smooth_density, "smooth_density") ;
	glut_viewer_add_key_func('f', extract_features, "extract_features") ;
	///dxy add
	glut_viewer_add_key_func('F', load_field, "load field") ;
	glut_viewer_add_key_func('S', all_edge_split/*edge_split*/, "edge split");
	glut_viewer_add_key_func('C', edge_collapse, "edge collapse");
	glut_viewer_add_key_func('A', update_area, "update area");
	glut_viewer_add_key_func('e', /*qd_experiment*/experiment, "experiment");
	glut_viewer_add_key_func('i', show_primal_info, "print primal info");
	glut_viewer_add_key_func(',', map_topo_optimize, "map topo optimize");
	glut_viewer_add_key_func('.', update_boundary, "update boundary from map");
	glut_viewer_add_key_func('/', save_opt_info, "save optimization information");
	glut_viewer_add_key_func('b', batch_snapshot, "batch snapshot");
	//test
	glut_viewer_add_key_func('w', quad_dominant_topo_optimize, "test");
	
	//glut_viewer_add_key_func('S', save_Map, "save map");
	///


	Geex::GPURVDApp app(argc, argv) ;
	app.main_loop() ;
	Geex::terminate() ;
}
