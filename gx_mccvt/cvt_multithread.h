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

#ifndef CVT_MULTITHREAD_H
#define CVT_MULTITHREAD_H

#include <Geex/basics/arrays.h>
#include <Geex/CVT/topo_poly_mesh.h>
#include <Geex/CVT/RVD.h>
//#include "cvt_geometry.h"
#include "delaunay_cvt.h"
//#include "cvt_differentiation.h"

namespace Geex {

	//dxy add
	class DelaunayCVT;

    /**
     * Represents the CVT computation that can be done on
     * a separate thread and the associated local storage.
     */
    class CVTThread {
    public:
 //       CVTThread(TopoPolyMesh* m,  CVTGeometry* geom) ;
		CVTThread(TopoPolyMesh* m,  DelaunayCVT* geom) ;
        
        void run_ComputeBarycenterWeighted() ;
        void run_ComputeLloydEnergyAndGradient() ;
        void run_ComputeLloydEnergyAndGradientAnisoN() ;
        void run_Compute_Linfty2_LloydEnergyAndGradient() ;
		void run_Compute_Symbolic_LloydEnergyAndGradient() ;
		//dxy add
		void run_ComputeCVTPrepareFCVT();
		//

        std::vector<double>& m() { return m_ ; }
        std::vector<vec3>& mg() { return mg_ ; }
        std::vector<double>& g() { return g_ ; }
        double& f() { return f_ ; }
        RestrictedVoronoiDiagram_poly& RVD() { return RVD_ ; }
		
        void set_factor(double factor) {
            factor_ = factor;
        }

		//dxy add
		std::vector<unsigned int> &project_facet_id() { return project_facet_id_; };
		std::vector<std::map<unsigned int, std::vector<std::pair<TopoPolyVertexEdge, TopoPolyVertexEdge>>>> &dual_segments() { return dual_segments_; }
		std::vector<double> &min_project_distance() { return min_project_distance_; }
		std::vector<double> &cellCVTEnergy() { return cellCVTEnergy_; }
		std::vector<double> &cell_area() { return cell_area_; }
		//


    protected:
        void init_m_and_mg() {
            g_.resize(0) ;
            m_.resize(RVD_.delaunay()->nb_vertices()) ;
            Memory::clear(m_) ;
            mg_.resize(RVD_.delaunay()->nb_vertices()) ;
            Memory::clear(mg_) ;
        }

        void init_f_and_g() {
            m_.resize(0) ;
            mg_.resize(0) ;
  /*          if (geom_->has_power_weight())
                g_.resize(4*RVD_.delaunay()->nb_vertices()) ;
            else*/
                g_.resize(3*RVD_.delaunay()->nb_vertices()) ;
            Memory::clear(g_) ;
            f_ = 0.0 ;
        }

		//dxy add
		void init_fcvt() {
			unsigned int nv = RVD_.delaunay()->nb_vertices();
			g_.resize(3*nv) ;
			Memory::clear(g_) ;
			f_ = 0.0;
			project_facet_id_.assign(nv, -1);
			dual_segments_.clear();
			dual_segments_.resize(nv);
			cellCVTEnergy_.assign(nv, 0);
			cell_area_.assign(nv, 0);
		}
		///

    private:
        RestrictedVoronoiDiagram_poly RVD_ ;
 //       CVTGeometry* geom_;
		DelaunayCVT* geom_;
        CVTThread(const CVTThread& rhs) { }
        CVTThread& operator=(const CVTThread& rhs) { }
        std::vector<double> m_ ;
        std::vector<vec3> mg_ ;
        std::vector<double> g_ ;
        double f_ ;
//		Linfty_Data_Block data_block_;
		double factor_;
		//dxy add: fcvt
		std::vector<unsigned int> project_facet_id_;
		std::vector<std::map<unsigned int, std::vector<std::pair<TopoPolyVertexEdge, TopoPolyVertexEdge>>>> dual_segments_;
		std::vector<double> min_project_distance_;
		std::vector<double> cellCVTEnergy_;
		std::vector<double> cell_area_;
		///
    } ;


    /**
     * A pointer to a member function of CVTThread.
     */
    typedef void (CVTThread::*CVTFunc)() ;


    /**
     * CVTMultiThread partitions the boundary (TopoPolyMesh) into 
     * multiple parts (one per core) and calls a member function
     * of CVTThread in parallel on each core, then gathers the
     * data in the specified target.
     */
    class CVTMultiThread {
    public:
        CVTMultiThread() ;
        ~CVTMultiThread() ;
		//dxy add
		enum Target {CVT, FCVT, OTHERS};
		//
       // void initialize(CVTGeometry* geom) ;
		void initialize(DelaunayCVT* geom) ;
        void run(CVTFunc func) ;
        void set_target(std::vector<vec3>& mg, std::vector<double>& m) {
            m_ = &m ; 
            mg_ = &mg ;
            //funcgrad_mode_ = false ;
			target_ = OTHERS; //dxy
        }
        void set_target(double& f, double* g) {
            f_ = &f ;
            g_ = g ;
            //funcgrad_mode_ = true ;
			target_ = CVT; //dxy 
        }
        unsigned int nb_threads() const { return parts_.size() == 0? 1:parts_.size() ; }

        void set_factor(double factor){
            for(unsigned int i=0; i < threads_.size(); i++)
                threads_[i]->set_factor(factor);
        }

		//dxy add: fcvt
		void set_target(double& f, double* g, std::vector<unsigned int> &Idx_in,
						std::vector<std::map<unsigned int, std::vector<std::pair<TopoPolyVertexEdge, TopoPolyVertexEdge>>>> &dual_segments_in,
						std::vector<double> &cellCVTEnergy_in, std::vector<double> &cell_area_in)
		{
			f_ = &f;
			g_ = g;
			project_facet_id_ = &Idx_in;
			dual_segments_ = &dual_segments_in;
			cellCVTEnergy_ = &cellCVTEnergy_in;
			cell_area_ = &cell_area_in;
			target_ = FCVT;
		}
		//

    private:
        Delaunay* delaunay_ ;
        Array1d<TopoPolyMesh> parts_ ;
        std::vector<CVTThread*> threads_ ;
        std::vector<double>* m_ ;
        std::vector<vec3>* mg_ ;
        //CVTGeometry* geom_;
		DelaunayCVT* geom_ ;
        double* f_ ;
        double* g_ ;
		//dxy add: fcvt
		std::vector<unsigned int> *project_facet_id_;
		std::vector<std::map<unsigned int, std::vector<std::pair<TopoPolyVertexEdge, TopoPolyVertexEdge>>>> *dual_segments_;
		//std::vector<double> *min_project_distance;
		std::vector<double> *cellCVTEnergy_;
		std::vector<double> *cell_area_;
		//
		Target target_;
		///
        bool funcgrad_mode_ ;
    } ;
    
}

#endif
