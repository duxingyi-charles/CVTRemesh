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

#include "cvt_multithread.h"
#include "cvt_geometry.h"
#include "cvt_compute.h"

#include <Geex/basics/processor.h>

#include <Geex/basics/stopwatch.h>


#ifdef USE_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#else
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#endif

#include <omp.h>




namespace Geex {

    CVTThread::CVTThread(TopoPolyMesh* mesh, DelaunayCVT* geom) : 
		RVD_(geom->delaunay(), mesh), geom_(geom), factor_(1.0) {

    }
   

    void CVTThread::run_ComputeBarycenterWeighted() {
        init_m_and_mg() ;
        RVD_.for_each_triangle(ComputeBarycenterw(mg_, m_)) ;
    }
    
    void CVTThread::run_ComputeLloydEnergyAndGradient() {
        init_f_and_g() ;
        RVD_.for_each_triangle(
            ComputeLloydEnergyAndGradientPtr(f_, &(g_[0]), RVD_.delaunay()->vertices()) 
        ) ;
    }

    void CVTThread::run_ComputeLloydEnergyAndGradientAnisoN() {
        init_f_and_g() ;
        RVD_.for_each_triangle(
            ComputeLloydEnergyAndGradientAnisoNPtr(
                RVD_, f_, &(g_[0]), RVD_.delaunay()->vertices(), 5.0 // TODO: get normal_anisotropy_ from CVTNumerics
            ) 
        ) ;
    }

	//dxy add: FCVT
	void CVTThread::run_ComputeCVTPrepareFCVT() {
		init_fcvt();
		RVD_.set_symbolic(true);
		RVD_.for_each_facet(
			ComputeCVTPrepareFCVT_mthread(
			f_,  &(g_[0]), 
			project_facet_id_, dual_segments_, min_project_distance_,
			cellCVTEnergy_, cell_area_,
			geom_->vertices(), RVD_
			)
		); 	
	}

	///

 //   void CVTThread::run_Compute_Linfty2_LloydEnergyAndGradient() {
 //       init_f_and_g() ;
 //       bool sym_backup = RVD_.symbolic() ;
 //       RVD_.set_symbolic(true);
 //       RVD_.for_each_triangle(
 //           Compute_Linfty2_EnergyAndGradientCentroidcenter(
 //           geom_, data_block_, f_, &(g_[0]), factor_
 //           )
 //       ) ;
 //       RVD_.set_symbolic(sym_backup);    
 //   }

	//void CVTThread::run_Compute_Symbolic_LloydEnergyAndGradient() {
	//	init_f_and_g() ;
	//	bool sym_backup = RVD_.symbolic() ;
	//	RVD_.set_symbolic(true);
	//	RVD_.for_each_triangle(
	//		Compute_Symbolic_EnergyAndGradientCentroidcenter(
	//		geom_, data_block_, f_, &(g_[0]), factor_
	//		)
	//		) ;
	//	RVD_.set_symbolic(sym_backup);    
	//}

    // ---------------------------------------------

    CVTMultiThread::CVTMultiThread() : delaunay_(nil), m_(nil), mg_(nil) {
    }

    CVTMultiThread::~CVTMultiThread() 
    {
        for(unsigned int i=0; i<threads_.size(); i++) {
            delete threads_[i] ;
        }
        threads_.clear();
    }

    void CVTMultiThread::initialize(DelaunayCVT* geom) {
        //geom_ = geom;
        geom_ = geom ; // geom->delaunay();
		delaunay_ = geom->delaunay();
        for(unsigned int i=0; i<threads_.size(); i++) {
            delete threads_[i] ;
        }
        parts_.clear() ;
        unsigned int nb_cores = Geex::Processor::number_of_cores() ;
        
        std::cerr << "CVTMultiThread: initializing for " 
                  << nb_cores << " cores" << std::endl ;
        threads_.resize(nb_cores) ;
        if (nb_cores > 1) {
            //geom->boundary()->partition(nb_cores, parts_) ;
			geom_->boundary()->partition(nb_cores, parts_) ;
            for(unsigned int i=0; i<nb_cores; i++) {
            //    threads_[i] = new CVTThread(&parts_[i], geom) ;
				threads_[i] = new CVTThread(&parts_[i], geom_) ;
            }
        } else {
        //    threads_[0] = new CVTThread(geom->boundary(), geom);
			threads_[0] = new CVTThread(geom_->boundary(), geom_);
        }
    }

#ifdef USE_TBB
    class CVTTask{
        private:
            CVTFunc f_;
        public:
            void operator()( 
                const tbb::blocked_range<std::vector<CVTThread*>::iterator>& r 
            ) const {
                for( 
                    std::vector<CVTThread*>::iterator i=r.begin(); i!=r.end(); ++i 
                ) {
                    ((*i)->*f_)();
                }
            }
            CVTTask(CVTFunc f):f_(f){}
    };
#endif

    void CVTMultiThread::run(CVTFunc f) {

        // Make sure the skeleton is updated before entering
        // multithreading (skeleton creation is not thread-safe !)
        const DelaunaySkeleton* skel = delaunay_->skeleton() ;
        
        //Geex::SystemStopwatch W ;
    
        
#ifdef USE_TBB
        tbb::blocked_range<std::vector<CVTThread*>::iterator> range(
            threads_.begin(),threads_.end()
        ); 
        parallel_for( range, CVTTask(f), tbb::auto_partitioner() );
#else
        // Create a thread group that calls the specified member function
        // on all threads.
        boost::thread_group threads_impl ;
        for(unsigned int i=0; i<nb_threads(); i++) {
            threads_impl.create_thread(boost::bind(f, threads_[i])) ;
        }
        // Wait for termination of all threads.
        threads_impl.join_all() ;
#endif
        //W.print_elapsed_time() ;
        
        // Gather thread local storage into target.
		switch (target_)
		{
		case CVT:
			{
				int nbx = delaunay_->nb_vertices() * 3 ;
				for(unsigned int i=0; i<nb_threads(); i++) {
					CVTThread& th = *(threads_[i]) ;
					int j = 0;
#pragma omp parallel for private(j)
					for(j=0; j<nbx; j++) {
						g_[j] += th.g()[j] ;
					}
					*f_ += th.f() ;
				}
			}
			break;
		case FCVT:
			{
				int nv = delaunay_->nb_vertices();
				int nbx = delaunay_->nb_vertices() * 3;
				std::vector<double> minDist(nv, -1);
				cellCVTEnergy_->assign(nv, 0);
				cell_area_->assign(nv, 0);
				for (unsigned int i=0; i<nb_threads(); ++i) {
					CVTThread& th = *(threads_[i]);
					int j = 0;
#pragma omp parallel for private(j)
					for (j=0; j<nbx; ++j) {
						g_[j] += th.g()[j];
					}
					*f_ += th.f();
					int k = 0;
#pragma omp parallel for private(k)
					for (k=0; k<nv; ++k) {
						(*cell_area_)[k] += th.cell_area()[k];
					}
					for (k=0; k<nv; ++k) {
						(*cellCVTEnergy_)[k] += th.cellCVTEnergy()[k];
					}
					for (k=0; k<nv; ++k) {
						if (th.min_project_distance()[k] >=0 &&
							(minDist[k]==-1 || th.min_project_distance()[k] < minDist[k]))
						{
							(*project_facet_id_)[k] = parts_[i].facet_info(th.project_facet_id()[k]).id_before_partition;
						}
					}
					for (k=0; k<nv; ++k) {
						std::map<unsigned int, std::vector<std::pair<TopoPolyVertexEdge, TopoPolyVertexEdge>>>  segs = th.dual_segments()[k];
						for (auto l=segs.begin(); l!=segs.end(); ++l) {
							auto &cur_edges = (*dual_segments_)[k][l->first];
							cur_edges.insert(cur_edges.end(), l->second.begin(), l->second.end());
						}				
					}
				}
			}
			break;
		case OTHERS:
			{
				for(unsigned int i=0; i<nb_threads(); i++) {
					CVTThread& th = *(threads_[i]) ;
					int j = 0;
					int nb = (int) delaunay_->nb_vertices();
#pragma omp parallel for private(j)
					for(j=0; j<nb; j++) {
						(*m_)[j] += th.m()[j] ;
					}
					for(j=0; j<nb; j++) {
						(*mg_)[j] += th.mg()[j] ;
					}
				}
			}
			break;
		default:
			break;
		}
//         if(funcgrad_mode_) {
// //            int nbx = geom_->has_power_weight()? delaunay_->nb_vertices() * 4 : delaunay_->nb_vertices() * 3 ;
// 			int nbx = delaunay_->nb_vertices() * 3 ;
//             for(unsigned int i=0; i<nb_threads(); i++) {
//                 CVTThread& th = *(threads_[i]) ;
// 				int j = 0;
// #pragma omp parallel for private(j)
//                 for(j=0; j<nbx; j++) {
//                     g_[j] += th.g()[j] ;
//                 }
//                 *f_ += th.f() ;
//             }
//         } else {
//             for(unsigned int i=0; i<nb_threads(); i++) {
//                 CVTThread& th = *(threads_[i]) ;
// 				int j = 0;
// 				int nb = (int) delaunay_->nb_vertices();
// #pragma omp parallel for private(j)
//                 for(j=0; j<nb; j++) {
//                     (*m_)[j] += th.m()[j] ;
//                 }
//                 for(j=0; j<nb; j++) {
//                     (*mg_)[j] += th.mg()[j] ;
//                 }
//             }
//         }
    }
} 

