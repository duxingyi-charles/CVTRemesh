#ifndef __CVT_COMPUTE_H_
#define __CVT_COMPUTE_H_

#include <Geex/CVT/geometry.h>
#include <fstream>
#include <map>
//dxy add
#include <Geex/combinatorics/map_builder.h>  //used by BuildMapFacet class
//
//

namespace Geex {


	/**
	* To be used as a template argument to RVD::for_each_triangle()
	* Computes for each RVD cell:
	*   g (gradient)
	*   f (Lloyd energy)
	*/
	class ComputeLloydEnergyAndGradient {
	public:
		ComputeLloydEnergyAndGradient(
			double& f_in,
			double* g_in,
			std::vector<double> &cell_f_in,
			std::vector<double> &cell_area_in,
			const std::vector<vec3>& seed_in, 
			const double factor = 1.0
			) : f(f_in), g(g_in), cell_f(cell_f_in), cell_area(cell_area_in), seed(seed_in), factor_(factor) {
				cell_f.assign(seed.size(), 0);
				cell_area.assign(seed.size(), 0);
		}


		void operator() (
			unsigned int v, 
			const TopoPolyVertexEdge& v1, 
			const TopoPolyVertexEdge& v2, 
			const TopoPolyVertexEdge& v3
			) const {
				double t_area = factor_ * tri_area(v1, v2, v3);
				vec3 sp[3];
				sp[0] = seed[v] - v1;
				sp[1] = seed[v] - v2;
				sp[2] = seed[v] - v3;
				double cur_f = 0;
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j <= i; j++) {
						cur_f += 5 * dot(sp[i], sp[j]);
					}
				}
				f += t_area * cur_f / 30;
				cell_f[v] += t_area * cur_f / 30;
				cell_area[v] += t_area;
				vec3 tmp_g = (t_area/6) * (12*seed[v] - (4*v1+4*v2+4*v3));
				g[3*v  ] += tmp_g[0];
				g[3*v+1] += tmp_g[1];
				g[3*v+2] += tmp_g[2];
		} 

		double& f ;
		double* g ;
		std::vector<double> &cell_f;
		std::vector<double> &cell_area;
		const std::vector<vec3>& seed ;
		double factor_;
	} ;


	///dxy add 
	/**
	* To be used as a template argument to RVD::for_each_triangle()
	* Computes for each RVD cell:
	*   g (gradient)
	*   f (Lloyd energy)
	*	Idx (index of the origin triangle to which the seed is projected)
	*/
	class ComputeLloydEnergyAndGradientAndIdx {
	public:
		ComputeLloydEnergyAndGradientAndIdx(
			double& f_in,
			double* g_in,
			const std::vector<vec3>& seed_in, 
			std::vector<Geex::vec3g<unsigned int>> &Idx_in,
			const RestrictedVoronoiDiagram_poly& rvd,
			const double factor = 1.0
			) : f(f_in), g(g_in), seed(seed_in), factor_(factor), Idx(Idx_in) , RVD_(rvd) , B(rvd.mesh())
		{
			//init
			last_visited_seed = new std::vector<int>(B->nb_facets(), -1);
			minDist = new std::vector<double>(seed_in.size(), -1);
		}

		~ComputeLloydEnergyAndGradientAndIdx() {
			delete last_visited_seed;
			delete minDist;
		}


		void operator() (
			unsigned int v, 
			const TopoPolyVertexEdge& v1, 
			const TopoPolyVertexEdge& v2, 
			const TopoPolyVertexEdge& v3
			) const {
				///compute Idx
				unsigned int fi = RVD_.current_facet();
				
				if ((*last_visited_seed)[fi] != v)
				{
					(*last_visited_seed)[fi] = v;
					
					//calc distance
					vec3 n = cross(v2-v1, v3-v1);
					n = normalize(n);

					vec3 Xc = seed[v];
					vec3 v1Xc = Xc - v1;
					double dist = abs(dot(v1Xc, n));

					//update Idx
					if ((*minDist)[v] < 0 || dist < (*minDist)[v]) {
						(*minDist)[v] = dist;
						//get triangle vertices
						vec3g<unsigned int> vidx;
						for (int i=0; i<3; ++i)
						{
							vidx[i] = B->facet_begin(fi) + i;
							vidx[i] = B->vertex_index(vidx[i]);
						}
						Idx[v] = vidx;
					}					
				}
				//
				

				//compute CVT
				double t_area = factor_ * tri_area(v1, v2, v3);
				vec3 sp[3];
				sp[0] = seed[v] - v1;
				sp[1] = seed[v] - v2;
				sp[2] = seed[v] - v3;
				double cur_f = 0;
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j <= i; j++) {
						cur_f += 5 * dot(sp[i], sp[j]);
					}
				}
				f += t_area * cur_f / 30;
				vec3 tmp_g = (t_area/6) * (12*seed[v] - (4*v1+4*v2+4*v3));
				g[3*v  ] += tmp_g[0];
				g[3*v+1] += tmp_g[1];
				g[3*v+2] += tmp_g[2];
		} 

		double& f ;
		double* g ;
		const std::vector<vec3>& seed ;
		double factor_;
		//
		std::vector<Geex::vec3g<unsigned int>> &Idx;
		const RestrictedVoronoiDiagram_poly& RVD_;
		std::vector<int> *last_visited_seed;
		std::vector<double> *minDist;
		const TopoPolyMesh* B;
	} ;

	///dxy add end

	class ComputeLloydEnergyAndGradientwAndIdx {
	public:
		ComputeLloydEnergyAndGradientwAndIdx(
			double& f_in,
			double* g_in,
			const std::vector<vec3>& seed_in, 
			std::vector<Geex::vec3g<unsigned int>> &Idx_in,
			const RestrictedVoronoiDiagram_poly& rvd,
			const double factor = 1.0
			) : f(f_in), g(g_in), seed(seed_in), factor_(factor), Idx(Idx_in) , RVD_(rvd) , B(rvd.mesh())
		{
			//init
			last_visited_seed = new std::vector<int>(RVD_.mesh()->nb_facets(), -1);
			minDist = new std::vector<double>(seed_in.size(), -1);
		}

		~ComputeLloydEnergyAndGradientwAndIdx() {
			delete last_visited_seed;
			delete minDist;
		}

		void operator() (
			unsigned int v, 
			const TopoPolyVertexEdge& v1, 
			const TopoPolyVertexEdge& v2, 
			const TopoPolyVertexEdge& v3
			) const {
				///compute Idx
				unsigned int fi = RVD_.current_facet();

				if ((*last_visited_seed)[fi] != v)
				{
					(*last_visited_seed)[fi] = v;

					//calc distance
					vec3 n = cross(v2-v1, v3-v1);
					n = normalize(n);

					vec3 Xc = seed[v];
					vec3 v1Xc = Xc - v1;
					double dist = abs(dot(v1Xc, n));

					//update Idx
					if ((*minDist)[v] < 0 || dist < (*minDist)[v]) {
						(*minDist)[v] = dist;
						//get triangle vertices
						vec3g<unsigned int> vidx;
						for (int i=0; i<3; ++i)
						{
							vidx[i] = B->facet_begin(fi) + i;
							vidx[i] = B->vertex_index(vidx[i]);
						}
						Idx[v] = vidx;
					}					
				}
				//


				//compute CVT
				double t_area = factor_ * tri_area(v1, v2, v3);
				vec3 sp[3];
				sp[0] = seed[v] - v1;
				sp[1] = seed[v] - v2;
				sp[2] = seed[v] - v3;
				double Sp = v1.w + v2.w + v3.w;
				double rho[3], alpha[3];
				rho[0] = v1.w; 
				rho[1] = v2.w; 
				rho[2] = v3.w;
				alpha[0] = Sp + rho[0]; 
				alpha[1] = Sp + rho[1];	
				alpha[2] = Sp + rho[2];

				double cur_f = 0;
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j <= i; j++) {
						cur_f += (alpha[i] + rho[j]) * dot(sp[i], sp[j]);
					}
				}
				f += t_area * cur_f / 30;
				vec3 tmp_g = (t_area/6) * (4*Sp*seed[v] - (alpha[0]*v1+alpha[1]*v2+alpha[2]*v3));
				g[3*v  ] += tmp_g[0];
				g[3*v+1] += tmp_g[1];
				g[3*v+2] += tmp_g[2];

		} 
		double& f ;
		double* g ;
		const std::vector<vec3>& seed ;
		double factor_;
		//
		std::vector<Geex::vec3g<unsigned int>> &Idx;
		const RestrictedVoronoiDiagram_poly& RVD_;
		std::vector<int> *last_visited_seed;
		std::vector<double> *minDist;
		const TopoPolyMesh* B;
	} ;
	///dxy add end

	class ComputeLloydEnergyAndGradientw {
	public:
		ComputeLloydEnergyAndGradientw(
			double& f_in,
			double* g_in,
			std::vector<double> &cell_f_in,
			std::vector<double> &cell_area_in,
			const std::vector<vec3>& seed_in, 
			const double factor = 1.0
			) : f(f_in), g(g_in), cell_f(cell_f_in), cell_area(cell_area_in), seed(seed_in), factor_(factor) {
				cell_f.assign(seed.size(), 0);
				cell_area.assign(seed.size(), 0);
		}

		void operator() (
			unsigned int v, 
			const TopoPolyVertexEdge& v1, 
			const TopoPolyVertexEdge& v2, 
			const TopoPolyVertexEdge& v3
			) const {
				double t_area = factor_ * tri_area(v1, v2, v3);
				vec3 sp[3];
				sp[0] = seed[v] - v1;
				sp[1] = seed[v] - v2;
				sp[2] = seed[v] - v3;
				double Sp = v1.w + v2.w + v3.w;
				double rho[3], alpha[3];
				rho[0] = v1.w; 
				rho[1] = v2.w; 
				rho[2] = v3.w;
				alpha[0] = Sp + rho[0]; 
				alpha[1] = Sp + rho[1];	
				alpha[2] = Sp + rho[2];

				double cur_f = 0;
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j <= i; j++) {
						cur_f += (alpha[i] + rho[j]) * dot(sp[i], sp[j]);
					}
				}
				f += t_area * cur_f / 30;
				cell_f[v] += t_area * cur_f / 30;
				cell_area[v] += t_area;  //use origin area, not weighted area
				vec3 tmp_g = (t_area/6) * (4*Sp*seed[v] - (alpha[0]*v1+alpha[1]*v2+alpha[2]*v3));
				g[3*v  ] += tmp_g[0];
				g[3*v+1] += tmp_g[1];
				g[3*v+2] += tmp_g[2];

		} 
		double& f ;
		double* g ;
		std::vector<double> &cell_f;
		std::vector<double> &cell_area;
		const std::vector<vec3>& seed ;
		double factor_;
	} ;

	// Same as before, but uses an array of doubles
	// for seeds.
	class ComputeLloydEnergyAndGradientPtr {
	public:
		ComputeLloydEnergyAndGradientPtr(
			double& f_in,
			double* g_in,
			const double* seed_in,
			const double factor = 1.0
			) : f(f_in), g(g_in), seed_(seed_in), factor_(factor) {
		}

		void operator() (
			unsigned int v, 
			const TopoPolyVertexEdge& v1, 
			const TopoPolyVertexEdge& v2, 
			const TopoPolyVertexEdge& v3
			) const {
				vec3 seed_v(seed_[3*v], seed_[3*v+1], seed_[3*v+2]) ;
				double t_area = factor_ * tri_area(v1, v2, v3);
				vec3 sp[3];
				sp[0] = seed_v - v1;
				sp[1] = seed_v - v2;
				sp[2] = seed_v - v3;
				double Sp = v1.w + v2.w + v3.w;
				double rho[3], alpha[3];
				rho[0] = v1.w; 
				rho[1] = v2.w; 
				rho[2] = v3.w;
				alpha[0] = Sp + rho[0]; 
				alpha[1] = Sp + rho[1];	
				alpha[2] = Sp + rho[2];

				double cur_f = 0;
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j <= i; j++) {
						cur_f += (alpha[i] + rho[j]) * dot(sp[i], sp[j]);
					}
				}
				f += t_area * cur_f / 30;
				vec3 tmp_g = (t_area/6) * (4*Sp*seed_v - (alpha[0]*v1+alpha[1]*v2+alpha[2]*v3));
				g[3*v  ] += tmp_g[0];
				g[3*v+1] += tmp_g[1];
				g[3*v+2] += tmp_g[2];

		} 
		double& f ;
		double* g ;
		const double* seed_ ;
		double factor_;
	} ;

	class SavePrimalTriangle {
	public:
		SavePrimalTriangle(
			std::ofstream& out
			) : out_(&out) { 
		}
		void operator()(unsigned int i, unsigned int j, unsigned int k) const {
			(*out_) << "f " << i+1 << " " << j+1 << " " << k+1 << std::endl ;
		}

	private:
		std::ofstream* out_ ;
	} ;

	///dxy add:
	/**
	* To be used for a template argument to RVD::for_each_primal_triangle()
	*/
	class ComputePrimalNeighborhood {
	public:
		ComputePrimalNeighborhood(
			std::map<unsigned, std::set<unsigned>>& neigh_in
			) : neigh(neigh_in) {}

		void operator() (unsigned int i, unsigned int j, unsigned int k) const {
			neigh[i].insert(j);
			neigh[i].insert(k);
			neigh[j].insert(i);
			neigh[j].insert(k);
			neigh[k].insert(i);
			neigh[k].insert(j);
		}

		std::map<unsigned, std::set<unsigned>> &neigh ;
	};

	///dxy add:
	/**
	* To be used for a template argument to RVD::for_each_primal_triangle()
	*/
	class BuildMapFacet {
	public:
		BuildMapFacet(/*MapBuilder *builder_in*/) /*: builder(builder_in)*/ {}

		void operator() (unsigned int i, unsigned int j, unsigned int k) const {
			builder->begin_facet();
			builder->add_vertex_to_facet(i);
			builder->add_vertex_to_facet(j);
			builder->add_vertex_to_facet(k);
			builder->end_facet();
		}

		//MapBuilder* builder;
		static MapBuilder* builder;
	};

	///dxy add:
	/**
	* To be used for a template argument to RVD::for_each_primal_triangle()
	*/
	class CollectLongEdges {
	public:
		CollectLongEdges(std::vector<std::pair<unsigned, unsigned>> &edges_in,
			const std::vector<unsigned>& degree_in
			)
			: edges_(edges_in), degree_(degree_in) {
				incident_vertices = new std::set<unsigned>();
		}

		~CollectLongEdges() {
			delete incident_vertices;
		}

		void operator()(unsigned int i, unsigned int j, unsigned int k) const {
			unsigned int index[3] = {i, j, k};
			for (int id=0; id<3; ++id)
			{
				if (degree_[index[id]] + degree_[index[(id+1)%3]] >= 14)  //long edge condition
				{

					if (incident_vertices->find(index[id]) != incident_vertices->end() ||
						incident_vertices->find(index[(id+1)%3]) != incident_vertices->end()) continue;
					incident_vertices->insert(index[id]);
					incident_vertices->insert(index[(id+1)%3]);
					edges_.push_back(std::make_pair(index[id], index[(id+1)%3]));
				}
			}
		}

	private:
		std::vector<std::pair<unsigned, unsigned>> &edges_;
		const std::vector<unsigned> &degree_;
		std::set<unsigned> *incident_vertices;
	};

	/**
	* To be used for a template argument to RVD::for_each_primal_triangle()
	*/
	class CollectShortEdges {
	public:
		CollectShortEdges(std::vector<std::pair<unsigned, unsigned>> &edges_in,
			const std::vector<unsigned>& degree_in
			)
			: edges_(edges_in), degree_(degree_in) {
				incident_vertices = new std::set<unsigned>();
		}

		~CollectShortEdges() {
			delete incident_vertices;
		}

		void operator()(unsigned int i, unsigned int j, unsigned int k) const {
			unsigned int index[3] = {i, j, k};
			for (int id=0; id<3; ++id)
			{
				if (degree_[index[id]] + degree_[index[(id+1)%3]] <= 10)  //short edge condition
				{

					if (incident_vertices->find(index[id]) != incident_vertices->end() ||
						incident_vertices->find(index[(id+1)%3]) != incident_vertices->end()) continue;
					incident_vertices->insert(index[id]);
					incident_vertices->insert(index[(id+1)%3]);
					edges_.push_back(std::make_pair(index[id], index[(id+1)%3]));
				}
			}
		}

	private:
		std::vector<std::pair<unsigned, unsigned>> &edges_;
		const std::vector<unsigned> &degree_;
		std::set<unsigned> *incident_vertices;
	};

	///

	/**
	* To be used for a template argument to RVD::for_each_primal_triangle()
	*/
	class ComputePrimalAngles {
	public:
		ComputePrimalAngles(
			const std::vector<vec3>& seed_in,
			std::vector<double>& angles_in,
			std::vector<double>& min_in,
			std::vector<double>& max_in
			) : seed(seed_in), angles(angles_in), min_angles(min_in), max_angles(max_in) {}

		void operator() (unsigned int i, unsigned int j, unsigned int k) const {
			const vec3& v1=seed[i] ;
			const vec3& v2=seed[j] ;
			const vec3& v3=seed[k] ;
			vec3 ab = v2 - v1 ; //F.vertex[1] - F.vertex[0] ;
			vec3 bc = v3 - v2 ; //F.vertex[2] - F.vertex[1] ;
			vec3 ca = v1 - v3 ; //F.vertex[0] - F.vertex[2] ;
			real lab = length(ab) ;
			real lbc = length(bc) ;
			real lca = length(ca) ;
			real anglea = acos(-dot(ca, ab)/(lca*lab))*180. / M_PI ;
			real angleb = acos(-dot(ab, bc)/(lab*lbc))*180. / M_PI ;
			real anglec = acos(-dot(bc, ca)/(lbc*lca))*180. / M_PI ;

			angles.push_back(anglea);
			angles.push_back(angleb);
			angles.push_back(anglec);

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
			min_angles.push_back(min);
			max_angles.push_back(max);
		}

		const std::vector<vec3>& seed;
		std::vector<double> &angles;
		std::vector<double> &min_angles;
		std::vector<double> &max_angles;
	};


	/**
	* To be used for a template argument to RVD::for_each_primal_triangle()
	*/
	class ComputePrimalAreas {
	public:
		ComputePrimalAreas(
			const std::vector<vec3>& seed_in,
			std::vector<double>& areas_in
			) : seed(seed_in), areas(areas_in) {}

		void operator() (unsigned int i, unsigned int j, unsigned int k) const {
			const vec3& v1=seed[i] ;
			const vec3& v2=seed[j] ;
			const vec3& v3=seed[k] ;
			areas.push_back(tri_area(v1, v2, v3));
		}

		const std::vector<vec3>& seed;
		std::vector<double> &areas;
	};

	///dxy add end

	/**
	* To be used as a template argument to RVD::for_each_triangle()
	* Computes for each RVD cell:
	*   g (gradient)
	*   f (Lloyd energy)
	*/
	class ComputeLloydEnergyAndGradientAnisoN {
	public:
		ComputeLloydEnergyAndGradientAnisoN(
			//const CVTGeometry* geometry_in,
			const RestrictedVoronoiDiagram_poly& rvd,
			double& f_in,
			double* g_in,
			const std::vector<vec3>& seed_in, 
			double normal_aniso = 10.0
			) : RVD_(rvd)/*geometry_(geometry_in)*/, f(f_in), g(g_in), seed(seed_in), aniso_(normal_aniso) {
		}

		void operator() (
			unsigned int v, 
			const TopoPolyVertexEdge& v1, 
			const TopoPolyVertexEdge& v2, 
			const TopoPolyVertexEdge& v3
			) const {
				unsigned int fi = RVD_.current_facet() ;
				const TopoPolyMesh* B = RVD_.mesh() ; 
				//unsigned int fi = geometry_->RVD().current_facet() ;
				//TopoPolyMesh* B = ((CVTGeometry*)geometry_)->boundary() ;
				vec3 N = B->facet_normal(fi) ;
				N = normalize(N) ;
				//vec3 U0 = aniso(v1 - seed[v],N, aniso_);
				//vec3 U1 = aniso(v2 - seed[v],N, aniso_);
				//vec3 U2 = aniso(v3 - seed[v],N, aniso_);

				//double m = tri_area(v1,v2,v3) ;

				//f += m / 6.0 * (
				//    dot(U0,U0)+dot(U0,U1)+dot(U0,U2)+dot(U1,U1)+dot(U1,U2)+dot(U2,U2)
				//) ;

				//vec3 centroid = (1.0 / 3.0) * (v1 + v2 + v3) ;
				//vec3 tmp_g = 2.0 * m * aniso(seed[v] - centroid, N, gx_sqr(aniso_)) ;
				//g[3*v  ] += tmp_g[0];
				//g[3*v+1] += tmp_g[1];
				//g[3*v+2] += tmp_g[2];

				double t_area = tri_area(v1, v2, v3);
				vec3 sp[3];
				sp[0] = aniso(v1 - seed[v], N, aniso_);
				sp[1] = aniso(v2 - seed[v], N, aniso_);
				sp[2] = aniso(v3 - seed[v], N, aniso_);
				double Sp = v1.w + v2.w + v3.w;
				double rho[3], alpha[3];
				rho[0] = v1.w; 
				rho[1] = v2.w; 
				rho[2] = v3.w;
				alpha[0] = Sp + rho[0]; 
				alpha[1] = Sp + rho[1];	
				alpha[2] = Sp + rho[2];

				double cur_f = 0;
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j <= i; j++) {
						cur_f += (alpha[i] + rho[j]) * dot(sp[i], sp[j]);
					}
				}
				f += t_area * cur_f / 30;
				//vec3 tmp_g = (t_area/6) * (4*Sp*p0 - (alpha[0]*v1+alpha[1]*v2+alpha[2]*v3));
				double m = tri_mass(v1, v2, v3, v1.w, v2.w, v3.w) ;
				double V ;
				vec3   Vg ;
				tri_centroid(v1, v2, v3, v1.w, v2.w, v3.w, Vg, V) ;
				vec3 tmp_g = 2*m*aniso(seed[v]-Vg/V, N, gx_sqr(aniso_));
				g[3*v  ] += tmp_g[0];
				g[3*v+1] += tmp_g[1];
				g[3*v+2] += tmp_g[2];
		} 

		// N is supposed to be normalized
		inline vec3 aniso(const vec3& x, const vec3& N, double s) const {
			return x + (s - 1.0) * dot(x,N) * N ;   
		}

		//const CVTGeometry* geometry_ ;
		const RestrictedVoronoiDiagram_poly& RVD_ ;
		double& f ;
		double* g ;
		const std::vector<vec3>& seed ;
		double aniso_ ;
	} ;

	/**
	* To be used as a template argument to RVD::for_each_triangle()
	* Computes for each RVD cell:
	*   g (gradient)
	*   f (Lloyd energy)
	*/
	class ComputeLloydEnergyAndGradientAnisoNPtr {
	public:
		ComputeLloydEnergyAndGradientAnisoNPtr(
			const RestrictedVoronoiDiagram_poly& rvd,
			double& f_in,
			double* g_in,
			const double* seed_in, 
			double normal_aniso = 10.0
			) : RVD_(rvd), f(f_in), g(g_in), seed_(seed_in), aniso_(normal_aniso) {
		}

		void operator() (
			unsigned int v, 
			const TopoPolyVertexEdge& v1, 
			const TopoPolyVertexEdge& v2, 
			const TopoPolyVertexEdge& v3
			) const {
				unsigned int fi = RVD_.current_facet() ;
				const TopoPolyMesh* B = RVD_.mesh() ; 
				vec3 N = B->facet_normal(fi) ;
				N = normalize(N) ;
				vec3 seed_v(seed_[3*v], seed_[3*v+1], seed_[3*v+2]) ;
				vec3 U0 = aniso(v1 - seed_v,N, aniso_);
				vec3 U1 = aniso(v2 - seed_v,N, aniso_);
				vec3 U2 = aniso(v3 - seed_v,N, aniso_);

				double m = tri_area(v1,v2,v3) ;

				f += m / 6.0 * (
					dot(U0,U0)+dot(U0,U1)+dot(U0,U2)+dot(U1,U1)+dot(U1,U2)+dot(U2,U2)
					) ;

				vec3 centroid = (1.0 / 3.0) * (v1 + v2 + v3) ;
				vec3 tmp_g = 2.0 * m * aniso(seed_v - centroid, N, gx_sqr(aniso_)) ;
				g[3*v  ] += tmp_g[0];
				g[3*v+1] += tmp_g[1];
				g[3*v+2] += tmp_g[2];
		} 

		// N is supposed to be normalized
		inline vec3 aniso(const vec3& x, const vec3& N, double s) const {
			return x + (s - 1.0) * dot(x,N) * N ;   
		}

		const RestrictedVoronoiDiagram_poly& RVD_ ;
		double& f ;
		double* g ;
		const double* seed_ ;
		double aniso_ ;
	} ;

	//
	static inline void lloyd_funcgrad(
		unsigned int v, vec3& p0, 
		const TopoPolyVertexEdge& v1, 
		const TopoPolyVertexEdge& v2, 
		const TopoPolyVertexEdge& v3, 
		double& f, double* g) 
	{
		double t_area = tri_area(v1, v2, v3);
		vec3 sp[3];
		sp[0] = p0 - v1;
		sp[1] = p0 - v2;
		sp[2] = p0 - v3;
		double cur_f = 0;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j <= i; j++) {
				cur_f += 5 * dot(sp[i], sp[j]);
			}
		}
		f += t_area * cur_f / 30;
		vec3 tmp_g = (t_area/6) * (12*p0 - (4*v1+4*v2+4*v3));
		g[3*v  ] += tmp_g[0];
		g[3*v+1] += tmp_g[1];
		g[3*v+2] += tmp_g[2];
	}

	static inline void lloyd_funcgrad_w(
		unsigned int v, vec3& p0, 
		const TopoPolyVertexEdge& v1, 
		const TopoPolyVertexEdge& v2, 
		const TopoPolyVertexEdge& v3, 
		double& f, double* g) 
	{
		double t_area = tri_area(v1, v2, v3);
		vec3 sp[3];
		sp[0] = p0 - v1;
		sp[1] = p0 - v2;
		sp[2] = p0 - v3;
		double Sp = v1.w + v2.w + v3.w;
		double rho[3], alpha[3];
		rho[0] = v1.w; 
		rho[1] = v2.w; 
		rho[2] = v3.w;
		alpha[0] = Sp + rho[0]; 
		alpha[1] = Sp + rho[1];	
		alpha[2] = Sp + rho[2];

		double cur_f = 0;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j <= i; j++) {
				cur_f += (alpha[i] + rho[j]) * dot(sp[i], sp[j]);
			}
		}
		f += t_area * cur_f / 30;
		vec3 tmp_g = (t_area/6) * (4*Sp*p0 - (alpha[0]*v1+alpha[1]*v2+alpha[2]*v3));
		g[3*v  ] += tmp_g[0];
		g[3*v+1] += tmp_g[1];
		g[3*v+2] += tmp_g[2];
	}

	// N is supposed to be normalized
	static inline vec3 aniso(const vec3& x, const vec3& N, double s) {
		return x + (s - 1.0) * dot(x,N) * N ;   
	}

	//static inline void lloyd_funcgrad(
	//	unsigned int v, vec3& p0, 
	//	const TopoPolyVertexEdge& v1, 
	//	const TopoPolyVertexEdge& v2, 
	//	const TopoPolyVertexEdge& v3, 
	//	double& f, double* g) {
	//	vec3 seed_v = p0 ;
	//	vec3 U0 = v1 - seed_v,N;
	//	vec3 U1 = v2 - seed_v,N;
	//	vec3 U2 = v3 - seed_v,N;

	//	double m = tri_area(v1,v2,v3) ;

	//	f += m / 6.0 * (
	//		dot(U0,U0)+dot(U0,U1)+dot(U0,U2)+dot(U1,U1)+dot(U1,U2)+dot(U2,U2)
	//		) ;

	//	vec3 centroid = (1.0 / 3.0) * (v1 + v2 + v3) ;
	//	vec3 tmp_g = 2.0 * m * (seed_v - centroid) ;
	//	g[3*v  ] += tmp_g[0];
	//	g[3*v+1] += tmp_g[1];
	//	g[3*v+2] += tmp_g[2];
	//}

	static inline void lloyd_funcgrad_aniso(
		unsigned int v, vec3& p0, vec3& N, double normal_aniso,
		const TopoPolyVertexEdge& v1, 
		const TopoPolyVertexEdge& v2, 
		const TopoPolyVertexEdge& v3, 
		double& f, double* g) {
			vec3 seed_v = p0 ;
			vec3 U0 = aniso(v1 - seed_v,N, normal_aniso);
			vec3 U1 = aniso(v2 - seed_v,N, normal_aniso);
			vec3 U2 = aniso(v3 - seed_v,N, normal_aniso);

			double m = tri_area(v1,v2,v3) ;

			f += m / 6.0 * (
				dot(U0,U0)+dot(U0,U1)+dot(U0,U2)+dot(U1,U1)+dot(U1,U2)+dot(U2,U2)
				) ;

			vec3 centroid = (1.0 / 3.0) * (v1 + v2 + v3) ;
			vec3 tmp_g = 2.0 * m * aniso(seed_v - centroid, N, gx_sqr(normal_aniso)) ;
			g[3*v  ] += tmp_g[0];
			g[3*v+1] += tmp_g[1];
			g[3*v+2] += tmp_g[2];
	}

	static inline void lloyd_funcgrad_aniso_w(
		unsigned int v, vec3& p0, vec3& N, double normal_aniso,
		const TopoPolyVertexEdge& v1, 
		const TopoPolyVertexEdge& v2, 
		const TopoPolyVertexEdge& v3, 
		double& f, double* g) {
			//vec3 seed_v = p0 ;
			//vec3 U0 = aniso(v1 - seed_v,N, normal_aniso);
			//vec3 U1 = aniso(v2 - seed_v,N, normal_aniso);
			//vec3 U2 = aniso(v3 - seed_v,N, normal_aniso);

			//double m = tri_area(v1,v2,v3) ;

			//f += m / 6.0 * (
			//	dot(U0,U0)+dot(U0,U1)+dot(U0,U2)+dot(U1,U1)+dot(U1,U2)+dot(U2,U2)
			//	) ;

			//vec3 centroid = (1.0 / 3.0) * (v1 + v2 + v3) ;
			//vec3 tmp_g = 2.0 * m * aniso(seed_v - centroid, N, gx_sqr(normal_aniso)) ;
			//g[3*v  ] += tmp_g[0];
			//g[3*v+1] += tmp_g[1];
			//g[3*v+2] += tmp_g[2];

			double t_area = tri_area(v1, v2, v3);
			vec3 sp[3];
			sp[0] = aniso(v1 - p0, N, normal_aniso);
			sp[1] = aniso(v2 - p0, N, normal_aniso);
			sp[2] = aniso(v3 - p0, N, normal_aniso);
			double Sp = v1.w + v2.w + v3.w;
			double rho[3], alpha[3];
			rho[0] = v1.w; 
			rho[1] = v2.w; 
			rho[2] = v3.w;
			alpha[0] = Sp + rho[0]; 
			alpha[1] = Sp + rho[1];	
			alpha[2] = Sp + rho[2];

			double cur_f = 0;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j <= i; j++) {
					cur_f += (alpha[i] + rho[j]) * dot(sp[i], sp[j]);
				}
			}
			f += t_area * cur_f / 30;
			//vec3 tmp_g = (t_area/6) * (4*Sp*p0 - (alpha[0]*v1+alpha[1]*v2+alpha[2]*v3));
			double m = tri_mass(v1, v2, v3, v1.w, v2.w, v3.w) ;
			double V ;
			vec3   Vg ;
			tri_centroid(v1, v2, v3, v1.w, v2.w, v3.w, Vg, V) ;
			vec3 tmp_g = 2*m*aniso(p0-Vg/V, N, gx_sqr(normal_aniso));
			g[3*v  ] += tmp_g[0];
			g[3*v+1] += tmp_g[1];
			g[3*v+2] += tmp_g[2];
	}

	/*
	// compute the Voronoi edge length of the RVD
	*/

	inline unsigned int get_opposite_vertex(
		unsigned int v, 
		const TopoPolyVertexEdge& v1, 
		const TopoPolyVertexEdge& v2
		) {
			SymbolicVertex sym ;
			int nb1 = v1.sym.nb_bisectors() ;
			int nb2 = v2.sym.nb_bisectors() ;
			sets_intersect(v1.sym, v2.sym, sym) ;
			if(sym.nb_bisectors()==0) {
				//std::cout << "wrong neighbor" << std::endl ;  ///dxy test comment: so many "wrong neighbor"!
				//gx_debug_assert(sym.nb_bisectors()==1) ;
				return -1 ;
			}
			return sym.bisector(0) ;
	}

	///dxy add
	/**
	* To be used as a template argument to RVD::for_each_facet()
	* Compute for each RVD seed c:
	* func (cvt energy)
	* g (cvt gradient)
	* Idx (index of original facet in which c is projected)
	* dual segments (dual segments of every primal edge)
	*/
	class ComputeCVTPrepareFCVT {
	public:
		ComputeCVTPrepareFCVT(
			double &func_in,
			double *g_in,
			std::vector<unsigned int> &Idx_in,
			std::vector<std::map<unsigned int, std::vector<std::pair<TopoPolyVertexEdge, TopoPolyVertexEdge>>>> &dual_segments_in,
			std::vector<double> &cellCVTEnergy_in,
			std::vector<double> &cell_area_in,
			const std::vector<vec3> &seed_in,
			const RestrictedVoronoiDiagram_poly &rvd
			) : func(func_in), g(g_in), Idx(Idx_in), cellCVTEnergy(cellCVTEnergy_in), cell_area(cell_area_in) ,dualseg(dual_segments_in), seed(seed_in), RVD(rvd), B(rvd.mesh())
		{
			minDist = new std::vector<double>(seed_in.size(), -1);
			cellCVTEnergy.assign(seed.size(), 0);
			cell_area.assign(seed.size(), 0);
		}

		virtual ~ComputeCVTPrepareFCVT() 
		{
			delete minDist;
		}


		void operator() (
			unsigned int c,
			TopoPolyMesh* M
			) const 
		{
			///compute Idx
			updateIdx(c);	

			//compute CVT and dual length
			for(unsigned int f=0; f<M->nb_facets(); ++f) {
				unsigned int i0 = M->facet_begin(f);
				if(M->vertex(i0).check_flag(TopoPolyVertexEdge::INTERSECT)) {
					updateDualseg(c, M->vertex(i0), M->vertex(M->next_around_facet(f, i0)));
				}
				if(i0+1<M->facet_end(f) && M->vertex(i0+1).check_flag(TopoPolyVertexEdge::INTERSECT)) {
					updateDualseg(c, M->vertex(i0+1), M->vertex(M->next_around_facet(f, i0+1)));
				}
				for(unsigned int i = M->facet_begin(f)+1; i+1<M->facet_end(f); i++) {
					if(M->vertex(i+1).check_flag(TopoPolyVertexEdge::INTERSECT)) {
						updateDualseg(c, M->vertex(i+1), M->vertex(M->next_around_facet(f,i+1)));
					}
					calcCVT(c, M->vertex(i0), M->vertex(i), M->vertex(i+1));
				}
			}
		}


		void updateIdx(unsigned int c) const
		{
			unsigned int fi = RVD.current_facet();
			vec3 p  = seed[c];
			vec3 v0 = B->vertex(B->facet_begin(fi));
			vec3 v1 = B->vertex(B->facet_begin(fi)+1);
			vec3 v2 = B->vertex(B->facet_begin(fi)+2);
			
			double dist = 10000;
			vec3 proj;
			if (project_to_triangle(p, v0, v1, v2, proj)) // proj in triangle
			{
				dist = distance(p, proj);
			}
			else { // proj out of triangle
				vec3 d2;
				project_to_linesegment(p, v0, v1, proj, d2[0]);
				project_to_linesegment(p, v1, v2, proj, d2[1]);
				project_to_linesegment(p, v2, v0, proj, d2[2]);
				double mind2 = (d2[1]<d2[0]) ? d2[1] : d2[0];
				mind2 = (d2[2]<mind2) ? d2[2] : mind2;
				dist = sqrt(mind2);
			}

// 			vec3 n = B->facet_normal(fi);
// 			n = normalize(n);
// 
// 			const vec3 &Xc = seed[c];
// 			vec3 v0 = B->vertex(B->facet_begin(fi));
// 			vec3 v0Xc = Xc - v0;
// 			double dist = abs(dot(v0Xc, n));

			if ((*minDist)[c] < 0 || dist < (*minDist)[c]) {
				(*minDist)[c] = dist;
				Idx[c] = fi;
			}
		}


		void updateDualseg(
			unsigned int c, 
			const TopoPolyVertexEdge& v1, 
			const TopoPolyVertexEdge& v2
			) const
		{
			unsigned int j = get_opposite_vertex(c, v1, v2);
			if (j == -1) return;

			//double len = length(v1 - v2);

			std::map<unsigned int, std::vector<std::pair<TopoPolyVertexEdge, TopoPolyVertexEdge>>> &neighbor = dualseg[c];
			if(neighbor.find(j) == neighbor.end()) { //newly found neighbor 
				neighbor.insert(std::make_pair(j, std::vector<std::pair<TopoPolyVertexEdge, TopoPolyVertexEdge>>()));
				neighbor[j].push_back(std::make_pair(v1, v2));
			}
			else {  //exist
				neighbor[j].push_back(std::make_pair(v1, v2));
			}
		}


		virtual void calcCVT (
			unsigned int v, 
			const TopoPolyVertexEdge& v1, 
			const TopoPolyVertexEdge& v2, 
			const TopoPolyVertexEdge& v3
			) const 
		{
			double t_area = tri_area(v1, v2, v3);
			vec3 sp[3];
			sp[0] = seed[v] - v1;
			sp[1] = seed[v] - v2;
			sp[2] = seed[v] - v3;
			double cur_f = 0;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j <= i; j++) {
					cur_f += 5 * dot(sp[i], sp[j]);
				}
			}
			func += t_area * cur_f / 30;
			vec3 tmp_g = (t_area/6) * (12*seed[v] - (4*v1+4*v2+4*v3));
			g[3*v  ] += tmp_g[0];
			g[3*v+1] += tmp_g[1];
			g[3*v+2] += tmp_g[2];
			//
			cellCVTEnergy[v] += t_area * cur_f / 30;
			cell_area[v] += t_area;
			//
		} 


	protected:
		double &func;
		double *g;
		std::vector<unsigned int> &Idx;
		std::vector<std::map<unsigned int, std::vector<std::pair<TopoPolyVertexEdge, TopoPolyVertexEdge>>>> &dualseg;
		const std::vector<vec3> &seed;
		const RestrictedVoronoiDiagram_poly &RVD;
		const TopoPolyMesh* B;
		//
		std::vector<double> *minDist;
		std::vector<double> &cellCVTEnergy;
		std::vector<double> &cell_area;
	};


	///dxy add
	/**
	* To be used as a template argument to RVD::for_each_facet()
	* Compute for each RVD seed c:
	* func (cvt energy)
	* g (cvt gradient)
	* Idx (index of original facet in which c is projected)
	* dual length (dual edge length of every primal edge)
	*/
	class ComputeCVTPrepareFCVTw : public ComputeCVTPrepareFCVT {
	public:
		ComputeCVTPrepareFCVTw(
			double &func_in,
			double *g_in,
			std::vector<unsigned int> &Idx_in,
			std::vector<std::map<unsigned int, std::vector<std::pair<TopoPolyVertexEdge, TopoPolyVertexEdge>>>> &dual_segments_in,
			std::vector<double> &cellCVTEnergy_in,
			std::vector<double> &cell_area_in,
			const std::vector<vec3> &seed_in,
			const RestrictedVoronoiDiagram_poly &rvd
			) : ComputeCVTPrepareFCVT(func_in, g_in, Idx_in, dual_segments_in, cellCVTEnergy_in, cell_area_in, seed_in, rvd)
		{
			minDist = new std::vector<double>(seed_in.size(), -1);
			cellCVTEnergy.assign(seed.size(), 0);
			cell_area.assign(seed.size(), 0);

		}

		virtual void calcCVT (
			unsigned int v, 
			const TopoPolyVertexEdge& v1, 
			const TopoPolyVertexEdge& v2, 
			const TopoPolyVertexEdge& v3
			) const 
		{
			double t_area = tri_area(v1, v2, v3);
			vec3 sp[3];
			sp[0] = seed[v] - v1;
			sp[1] = seed[v] - v2;
			sp[2] = seed[v] - v3;
			double Sp = v1.w + v2.w + v3.w;
			double rho[3], alpha[3];
			rho[0] = v1.w; 
			rho[1] = v2.w; 
			rho[2] = v3.w;
			alpha[0] = Sp + rho[0]; 
			alpha[1] = Sp + rho[1];	
			alpha[2] = Sp + rho[2];

			double cur_f = 0;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j <= i; j++) {
					cur_f += (alpha[i] + rho[j]) * dot(sp[i], sp[j]);
				}
			}
			func += t_area * cur_f / 30;
			vec3 tmp_g = (t_area/6) * (4*Sp*seed[v] - (alpha[0]*v1+alpha[1]*v2+alpha[2]*v3));
			g[3*v  ] += tmp_g[0];
			g[3*v+1] += tmp_g[1];
			g[3*v+2] += tmp_g[2];
			//
			cellCVTEnergy[v] += t_area * cur_f / 30;
			cell_area[v] += t_area;
			//
		} 

	};

	///dxy add
	/**
	* To be used as a template argument to RVD::for_each_facet()
	* Compute for each RVD seed c:
	* func (cvt energy)
	* g (cvt gradient)
	* Idx (index of original facet in which c is projected)
	* dual length (dual edge length of every primal edge)
	*/
	class ComputeCVTPrepareFCVT_mthread : public ComputeCVTPrepareFCVTw {
	public:
		ComputeCVTPrepareFCVT_mthread(
			double &func_in,
			double *g_in,
			std::vector<unsigned int> &Idx_in,
			std::vector<std::map<unsigned int, std::vector<std::pair<TopoPolyVertexEdge, TopoPolyVertexEdge>>>> &dual_segments_in,
			std::vector<double> &minDist_in,
			std::vector<double> &cellCVTEnergy_in,
			std::vector<double> &cell_area_in,
			const std::vector<vec3> &seed_in,
			const RestrictedVoronoiDiagram_poly &rvd
			) : ComputeCVTPrepareFCVTw(func_in, g_in, Idx_in, dual_segments_in, cellCVTEnergy_in, cell_area_in, seed_in, rvd), minDistInfo(minDist_in)
		{ }

		virtual ~ComputeCVTPrepareFCVT_mthread() {
			minDistInfo = *minDist;
		}

	private:
		std::vector<double> &minDistInfo;
	
	};




	//////////////////////////////////////////////////////////////////////////

	class ComputeDualEdgeLength {
	public:
		ComputeDualEdgeLength(std::map<std::pair<int, int>, double>& dual_length) 
			: dual_length_(dual_length) {
		}

		void operator() (
			unsigned int iv1, 
			const TopoPolyVertexEdge& v1, 
			const TopoPolyVertexEdge& v2
			) const {
				unsigned iv2 = get_opposite_vertex(iv1, v1, v2) ;  ///test
				if(iv2!=-1) {
					std::pair<unsigned, unsigned> key ;
					if(iv1<iv2) 		
						key = std::make_pair(iv1, iv2) ;
					else 
						key = std::make_pair(iv2, iv1) ;
					if(dual_length_.find(key)!=dual_length_.end()) {
						dual_length_[key] += distance(v1, v2) ;
					}
					else 
						dual_length_[key] = distance(v1, v2) ;
				}

		} 
	private:
		std::map<std::pair<int, int>, double>& dual_length_ ;
	} ;

	/**
	* To be used as a template argument to RVD::for_each_triangle()
	* Computes for each RVD cell:
	*   mg (total area times centroid)
	*   m  (total area)
	* This version does not take the weights into account.
	*/
	class ComputeBarycenter {
	public:
		ComputeBarycenter(
			std::vector<vec3>& mg_in,
			std::vector<double>& m_in 
			) : mg(mg_in), m(m_in) {
		}

		void operator() (
			unsigned int v,
			const TopoPolyVertexEdge& v1, 
			const TopoPolyVertexEdge& v2, 
			const TopoPolyVertexEdge& v3
			) const {
				double area = 0.5 * length(cross(v2-v1, v3-v1)) ;
				m[v]  += area ;
				mg[v] += (area / 3.0) * (v1 + v2 + v3) ;
		} 

		std::vector<vec3>& mg ;
		std::vector<double>& m ;
	} ;

	class ComputeBarycenterw {
	public:
		ComputeBarycenterw(
			std::vector<vec3>& mg_in,
			std::vector<double>& m_in 
			) : mg(mg_in), m(m_in) {
		}

		void operator() (
			unsigned int v,
			const TopoPolyVertexEdge& v1, 
			const TopoPolyVertexEdge& v2, 
			const TopoPolyVertexEdge& v3
			) const {
				double V ;
				vec3 Vg ;
				tri_centroid(v1, v2, v3, v1.w, v2.w, v3.w, Vg, V) ; 
				m[v]  += V ;
				mg[v] += Vg ;
		} 

		std::vector<vec3>& mg ;
		std::vector<double>& m ;
	} ;

	class SaveRVDTriangles {
	public:
		SaveRVDTriangles(
			std::vector<std::vector<vec3>>& cell_fvs          
			) : cell_fvs_(cell_fvs) {
		}

		void operator()(
			unsigned int v,
			const TopoPolyVertexEdge& v1, 
			const TopoPolyVertexEdge& v2, 
			const TopoPolyVertexEdge& v3
			) const {
				cell_fvs_[v].push_back(v1) ;
				cell_fvs_[v].push_back(v2) ;
				cell_fvs_[v].push_back(v3) ;
		}        
	private:
		std::vector<std::vector<vec3>>& cell_fvs_ ;
	} ;


	///dxy add
	/**
	* To be used as a template argument to RVD::for_each_facet()
	* Compute for each RVD seed c:
	* Idx (index of three original vertex of the triangle in which c is projected)
	*/
	class ComputeProjectedTriangleIndex {
	public:
		ComputeProjectedTriangleIndex (
			std::vector<vec3g<unsigned int>>& Idx_in,
			const std::vector<vec3>& seed_in,
			const RestrictedVoronoiDiagram_poly& rvd
			) : Idx(Idx_in), seed(seed_in), RVD_(rvd) {
				minDist = new std::vector<double>(seed_in.size(), -1);
		}

		~ComputeProjectedTriangleIndex() { delete minDist; }

		void operator() (
			unsigned int c,
			TopoPolyMesh* M
			) const {
				//get normal
				unsigned int fi = RVD_.current_facet();
				const TopoPolyMesh* B = RVD_.mesh();
				vec3 n = B->facet_normal(fi);
				n = normalize(n);

				//get triangle vertices
				if (B->facet_size(fi) < 3) {
					std::cout << "degenerate facet: less than 3 vertices" << std::endl;
					return; //skip
				}
				vec3g<unsigned int> vidx;
				for (int i=0; i<3; ++i)
				{
					vidx[i] = B->facet_begin(fi) + i;
					vidx[i] = B->vertex_index(vidx[i]);
				}

				//calc distance
				vec3 v1 = (B->original_vertices())[vidx[0]];
				vec3 Xc = seed[c];
				vec3 v1Xc = Xc - v1;
				double dist = abs(dot(v1Xc, n));

				//update Idx
				if ((*minDist)[c] < 0 || dist < (*minDist)[c]) {
					(*minDist)[c] = dist;
					Idx[c] = vidx;
				}
		}

	private:
		std::vector<vec3g<unsigned int>> &Idx;
		const std::vector<vec3>& seed;
		const RestrictedVoronoiDiagram_poly& RVD_;
		std::vector<double> *minDist;
	};


	/**
	* To be used as a template argument to RVD::for_each_facet()
	* Compute for each RVD seed:
	* N (normal direction)
	* D (field direction)
	*/
	class ComputeNormalAndDirection {
	public:
		ComputeNormalAndDirection (
			std::vector<vec3>& N_in,
			std::vector<vec3>& D_in,
			const std::vector<vec3>& seed_in,
			const RestrictedVoronoiDiagram_poly& rvd,
			unsigned int rot_symmetry
			) : N(N_in), D(D_in), seed(seed_in), RVD_(rvd), rot_symmetry_(rot_symmetry) {
				minDist = new std::vector<double>(seed_in.size(), -1);
		}

		~ComputeNormalAndDirection() {
			delete minDist;
		}


		void operator() (
			unsigned int c,
			TopoPolyMesh* M
			) const {
				
				//get normal
				unsigned int fi = RVD_.current_facet();
				const TopoPolyMesh* B = RVD_.mesh();
				vec3 n = B->facet_normal(fi);
				n = normalize(n);

				//get triangle vertices
				if (B->facet_size(fi) < 3) {
					std::cout << "degenerate facet: less than 3 vertices" << std::endl;
					return; //skip
				}
				unsigned int vidx[3];
				for (int i=0; i<3; ++i)
				{
					vidx[i] = B->facet_begin(fi) + i;
					vidx[i] = B->vertex_index(vidx[i]);
				}

				vec3 v1 = (B->original_vertices())[vidx[0]];


				//calc distance
				vec3 Xc = seed[c];
				vec3 v1Xc = Xc - v1;
				double dist = abs(dot(v1Xc, n));

				//update N and D
				if ((*minDist)[c] < 0 || dist < (*minDist)[c])
				{
					(*minDist)[c] = dist;
					N[c] = n;  //use facet normal, not interpolation of vertex normal
					vec3 p[3];
					vec3 field[3];
					for (int i=0; i<3; ++i)
					{
						p[i] = (B->original_vertices())[vidx[i]];
						field[i] = B->vertex_field(vidx[i]);
					}
					vec3 d = direction_interpolate(Xc, p, field, n, rot_symmetry_);
					d = d - dot(n, d) * n;
					d = normalize(d);
					D[c] = d;
				}
		}


		/** linear symmetry-representational direction interpolation, inspired by
		* project_to_triangle()
		*/
		vec3 direction_interpolate(
			vec3 p,
			vec3 pos[3],
			vec3 direction[3],
			vec3 normal,
			unsigned int rot_symmetry = 1
			) const
		{
			//get linear factor
			vec3 pc = p - pos[2];
			vec3 ac = pos[0] - pos[2];
			vec3 bc = pos[1] - pos[2];
			double a00 = ac.length2();
			double a01 = dot(ac, bc);
			double a11 = bc.length2();
			double b0 = dot(pc, ac);
			double b1 = dot(pc, bc);
			double mdet = a00 * a11 - a01 * a01;
			double s = (a11 * b0 - a01 * b1) / mdet;
			double t = (a00 * b1 - a01 * b0) / mdet;
			//proj = s * a + t * b + (1 - s - t ) * c;
			//Note: no check for s and t, may be negative
			

			//symmetry direction interpolation
			vec3 proj_direction[3];
			normal = normalize(normal);
			for (int i=0; i<3; ++i)
			{
				proj_direction[i] = direction[i] - dot(direction[i], normal) * normal;
				proj_direction[i] = normalize(proj_direction[i]);
			}
			vec3 e_x = proj_direction[0];
			vec3 e_y = cross(normal, e_x);

			double theta[3];
			theta[0] = 0.0;
			for (int i=1; i<3; ++i) {
				double dot_i = dot(proj_direction[i], e_x);
				dot_i = (dot_i>1)  ?  1 : dot_i;
				dot_i = (dot_i<-1) ? -1 : dot_i;
 				theta[i] = acos(dot_i);
				if (dot(proj_direction[i], e_y) < 0) {
					theta[i] = 2 * M_PI - theta[i];
				}
				theta[i] *= rot_symmetry;
			}

			vec2 sym_direction[3];
			for (int i=0; i<3; ++i) {
				sym_direction[i] = vec2(cos(theta[i]), sin(theta[i]));
			}
			vec2 sym_dir = s * sym_direction[0] + t * sym_direction[1] + (1-s-t) * sym_direction[2];
			sym_dir = normalize(sym_dir);
			double sym_angle = acos(sym_dir.x);
			if (sym_dir.y < 0) {
				sym_angle = 2 * M_PI - sym_angle;
			}
			double dir_angle = sym_angle / rot_symmetry;

			vec3 dir = cos(dir_angle) * e_x + sin(dir_angle) * e_y;

			//debug
			if (dir[0] != dir[0])
			{
				std::cout << "direction_interpolate: return NaN!" << std::endl;
			}
			return dir;
		}

		std::vector<vec3>& N;
		std::vector<vec3>& D;
		const std::vector<vec3>& seed;
		const RestrictedVoronoiDiagram_poly& RVD_;
		std::vector<double> *minDist;
		unsigned int rot_symmetry_;
	};





	/**
	* To be used as a template argument to RVD::for_each_halfedge()
	* Compute for each RVD cell:
	* g (direction gradient)
	* f (direction energy)
	*/
	class ComputeDirectionEnergyAndGradient {
	public:
		typedef std::pair<int, int> Int_pair;
		typedef std::pair<double, vec3> Funcgrad_pair;
		typedef std::map<Int_pair, Funcgrad_pair> Funcgrad_map;

		ComputeDirectionEnergyAndGradient(
			double& f_in,
			double* g_in,
			const std::vector<vec3>& seed_in,
			const std::vector<vec3>& N_in,
			const std::vector<vec3>& D_in,
			const double factor_in,
			const std::vector<double>& area_in,
			GLenum mode_in = Angle_Edge,
			double para_a = 1.5
			) : f(f_in), g(g_in), seed(seed_in), N(N_in), D(D_in), energy_mode_(mode_in), para_a_(para_a), factor_(factor_in), area_(area_in) {
				fgMap = new Funcgrad_map();
		}

		~ComputeDirectionEnergyAndGradient() {
			delete fgMap;
		}

		void operator() (
			unsigned int c,
			const TopoPolyVertexEdge& v1,
			const TopoPolyVertexEdge& v2
			) const {
				vec3 Xc = seed[c];
				vec3 n = N[c];
				vec3 d = D[c];
				unsigned int j = get_opposite_vertex(c, v1, v2);
				if (j == -1)  //why would j be -1 ?
				{
					return;
				}

				vec3 Xj = seed[j];

				//Funcgrad_map::iterator mit = fgMap.find(Int_pair(c, j));
				auto mit = fgMap->find(Int_pair(c, j));
				Funcgrad_pair E1G1;
				if (mit == fgMap->end())
				{
					E1G1 = computeE1G1(Xc, Xj, n, d);
					fgMap->insert(std::make_pair(Int_pair(c, j), E1G1));
				}
				else {
					E1G1 = mit->second;
				}

				//
				double E1 = E1G1.first;
				vec3   G1 = E1G1.second;
				double len = length(v1 - v2);
				double E2 = E1 * len;
				vec3   G2 = G1 * len;
				//test
				double E3 = E1 * len * len;
				vec3   G3 = G1 * len * len;

				//double factor = factor_[c];
				double cell_area = area_[c];
				if (energy_mode_ == Angle_Only)
				{
					double factor = factor_ * cell_area * cell_area;
					f += factor * E1;
					g[3*c  ] += factor * G1[0];
					g[3*c+1] += factor * G1[1];
					g[3*c+2] += factor * G1[2];
				}
				else if (energy_mode_ == Angle_Edge)
				{
					double factor = factor_ * cell_area * sqrt(cell_area);
					f += factor * E2;
					g[3*c  ] += factor * G2[0];
					g[3*c+1] += factor * G2[1];
					g[3*c+2] += factor * G2[2];
				}
				else if (energy_mode_ == Angle_Edge2)
				{
					double factor = factor_ * cell_area;
					f += factor * E3;
					g[3*c  ] += factor * G3[0];
					g[3*c+1] += factor * G3[1];
					g[3*c+2] += factor * G3[2];
				}
				//else: do nothing
		}

		/**
		* In polar coordinate,
		* E1 = a * A(phi) + B(phi) * C(theta)
		*/
		virtual Funcgrad_pair computeE1G1(
			const vec3& X1,
			const vec3& X2,
			const vec3& n,
			const vec3& d
			) const
		{
			vec3 X1X2 = X2 - X1;
			vec3 proj = X1X2 - dot(X1X2, n) * n ;

			double phi = acos(dot(X1X2, n) / length(X1X2));
			//test
			/*if (!(phi == phi))
			{
			std::cout << "phi " << phi << std::endl;
			}*/
			//
			double theta = 0.0;
			if (proj.length() > 0)
			{
				
				theta = acos(dot(proj, d) / length(proj));
				if (!(theta == theta)) //theta is NaN
				{
					if (dot(proj, d) > 0) 
					{
						theta = 0;
					}
					else
					{
						theta = M_PI;
					}
				}
				if (dot(cross(proj, d), n) > 0)
				{
					theta = 2 * M_PI - theta;
				}
			}

			//
			double a = para_a_;
			double A = abs(2 * phi/M_PI - 1);
			double B = 1 - A;
			double C = 0.5 * (1 - cos(6 * theta));

			//
// 			vec3 g_theta = cross(X1X2, n) / length(X1X2) ;
// 			vec3 g_phi   = cross(g_theta, X1X2) / length(X1X2);

			//test
			vec3 g_theta = cross(X1X2, n) / (length2(cross(X1X2, n)));
			vec3 g_phi   = cross(g_theta, X1X2) * length(cross(X1X2, n))/ length2(X1X2);
			//

			vec3 g_A = (2 / M_PI) * g_phi;
			if (phi < M_PI/2)
			{
				g_A *= -1;
			}
			vec3 g_B = -1 * g_A;
			vec3 g_C = (3 * sin(6 * theta)) * g_theta;

			//
			double E1 = a * A + B * C;
			vec3 g_E1 = a * g_A + B * g_C + C * g_B;

			//project g_E1 : solve the "shrink" problem
			g_E1 = g_E1 - dot(g_E1, n) * n;
			//

			return Funcgrad_pair(E1, g_E1);
		}

		double& f ;
		double* g ;
		const std::vector<vec3>& seed ;
		const std::vector<vec3>& N ;
		const std::vector<vec3>& D ;
		const double para_a_;
		double factor_;
		const std::vector<double>& area_;
		Funcgrad_map *fgMap;
		GLenum energy_mode_;
	};


	/**
	* To be used as a template argument to RVD::for_each_halfedge()
	* Compute for each RVD cell:
	* g (direction gradient)
	* f (direction energy)
	*/
	class ComputeDirectionEnergyAndGradient_True : public ComputeDirectionEnergyAndGradient
	{
	public:
		ComputeDirectionEnergyAndGradient_True(
			double& f_in,
			double* g_in,
			const std::vector<vec3>& seed_in,
			const std::vector<vec3>& N_in,
			const std::vector<vec3>& D_in,
			double factor_in,
			std::vector<double>& area_in,
			GLenum mode_in = Angle_Edge
			) : ComputeDirectionEnergyAndGradient(f_in, g_in, seed_in, N_in, D_in, factor_in, area_in, mode_in, 0/*para_a*/) {
		}


		/**
		* In polar coordinate,
		* E1 = (1 - sin(phi) * cos(6 theta)) /2
		*/
		virtual Funcgrad_pair computeE1G1( const vec3& X1, const vec3& X2, const vec3& n, const vec3& d ) const
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
			//
			double E1 = 0.5 * (1 - sin(phi) * cos(6*theta));
			vec3 g_E1 = -0.5 * cos(phi) * cos(6*theta) * g_phi + 3 * sin(phi) * sin(6*theta) * g_theta;

			//project g_E1 : solve the "shrink" problem
			g_E1 = g_E1 - dot(g_E1, n) * n;
			//
			return Funcgrad_pair(E1, g_E1);
		}
	};


	/**
	* To be used as a template argument to RVD::for_each_halfedge()
	* Compute for each RVD cell:
	* g (direction gradient)
	* f (direction energy)
	*/
	class ComputeDirectionEnergyAndGradient_Quad : public ComputeDirectionEnergyAndGradient
	{
	public:
		ComputeDirectionEnergyAndGradient_Quad(
			double& f_in,
			double* g_in,
			const std::vector<vec3>& seed_in,
			const std::vector<vec3>& N_in,
			const std::vector<vec3>& D_in,
			double factor_in,
			std::vector<double>& area_in,
			GLenum mode_in = Angle_Edge,
			double para_quad = 0.5,
			double para_a = 1.5
			) : ComputeDirectionEnergyAndGradient(f_in, g_in, seed_in, N_in, D_in, factor_in, area_in, mode_in, para_a), para_quad_(para_quad) {
		}


		/**
		* In polar coordinate,
		* E1 = a * A(phi) + B(phi) * C(theta, p)
		*/
		virtual Funcgrad_pair computeE1G1( const vec3& X1, const vec3& X2, const vec3& n, const vec3& d ) const
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
			double a = para_a_;
			double A = abs(2 * phi/M_PI - 1);
			double B = 1 - A;
			//
			double p = para_quad_;
			double q = 1-p;
			double C = 0.5 * (1 - p * cos(8 * theta) - q * cos(4 * theta));  //changed for quad meshing

			//
			vec3 g_theta = cross(X1X2, n) / length(X1X2);
			vec3 g_phi   = cross(g_theta, X1X2) / length(X1X2);

			vec3 g_A = (2 / M_PI) * g_phi;
			if (phi < M_PI/2)
			{
				g_A *= -1;
			}
			vec3 g_B = -1 * g_A;
			vec3 g_C = (4 * p * sin(8 * theta) + 2 * q * sin(4 * theta))  * g_theta;  //changed for quad meshing

			//
			double E1 = a * A + B * C;
			vec3 g_E1 = a * g_A + B * g_C + C * g_B;

			//project g_E1 : solve the "shrink" problem
			g_E1 = g_E1 - dot(g_E1, n) * n;
			//
			return Funcgrad_pair(E1, g_E1);
		}


		double para_quad_;  //mix factor of 4-RoSy and 8-RoSy field
	};


	/**
	* To be used for RVD::for_each_triangle()
	* Compute total number of triangles:
	*/
	class CountTriangle {
	public:
		CountTriangle(int& cnt_in) : cnt(cnt_in) { cnt = 0; }

		void operator() (
			unsigned int v,
			const TopoPolyVertexEdge& v1,
			const TopoPolyVertexEdge& v2,
			const TopoPolyVertexEdge& v3
			) const {
				++cnt;
		}
	private:
		int &cnt;
	};


	/** dxy add
	* To be used as a template augument to RVD::for_each_facet()
	* Compute feature line with which seed's rvd cell intersects
	*/
	class ComputeSeedFeatureLine {
	public:
		ComputeSeedFeatureLine (
			std::map<int, std::set<int>> &intersect_feature_line_in,
			const std::map<std::pair<int, int>, int> &featureEdgeLine_in,
			const RestrictedVoronoiDiagram_poly& RVD_in
			) : intersect_feature_line(intersect_feature_line_in), featureEdgeLine(featureEdgeLine_in), RVD_(RVD_in), B(RVD_in.mesh()) {}

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
								intersect_feature_line[c].insert(fLine[j]);
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
		std::map<int, std::set<int>> &intersect_feature_line;
	};


	/**
	* To be used as a template argument to RVD::for_each_facet()
	* Compute total number of facets:
	*/
	class CountFacet {
	public:
		CountFacet (int &cnt_in)  : cnt(cnt_in) { cnt = 0; }

		void operator() (
			unsigned int c,
			TopoPolyMesh* M
			) const {
				++cnt;
		}
	private:
		int &cnt;
	};

	/** 
	* To be used as a template argument to RVD::for_each_triangle()
	* Computes for each RVD cell:
	* area (cell area)
	*/
	class ComputeCellArea {
	public:
		ComputeCellArea(std::vector<double>& area_in) : area(area_in) {}

		void operator() (
			unsigned int v,
			const TopoPolyVertexEdge& v1,
			const TopoPolyVertexEdge& v2,
			const TopoPolyVertexEdge& v3
			) const {
				area[v] += tri_area(v1, v2, v3);
		}

	private:
		std::vector<double>& area;
	};

	///dxy add end
}


#endif
