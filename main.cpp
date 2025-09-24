//
// Created by pengfei on 2025/9/23.
//
#include <iostream>
#include <sstream>
#include <fstream>

#include <vector>
#include <set>
#include <queue>
#include<math.h>
#include <map>
#include <algorithm>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
//#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Aff_transformation_3.h>



#include <Eigen/Dense>
#include <unordered_set>
#include <Eigen/Dense>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/draw_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>


// typedef CGAL::Simple_cartesian<double> K;  // 使用double类型的核
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_3<K>                   DT3;
typedef CGAL::Creator_uniform_3<double, K::Point_3>          Creator;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> Points;


typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>      Vb;
typedef CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3<K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                            Tds;
//typedef CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Parallel_tag>                            Tds;

typedef CGAL::Delaunay_triangulation_3<K, Tds>                    Delaunay;

using namespace std;

typedef CGAL::Polyhedron_3<K>  Polyhedron;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Nef_polyhedron_3<K>  Nef_polyhedron;
typedef K::Point_3 Point_3;
typedef K::Point_2 Point_2;
typedef K::Vector_3 Vector_3;
typedef K::Plane_3 Plane_3;
typedef CGAL::Aff_transformation_3<K> Aff_transformation_3;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;


typedef CGAL::Delaunay_triangulation_3<K>                   DT3;
typedef CGAL::Creator_uniform_3<double, K::Point_3>          Creator;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> Points;


typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>      Vb;
typedef CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3<K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                            Tds;
//typedef CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Parallel_tag>                            Tds;

typedef CGAL::Delaunay_triangulation_3<K, Tds>                    Delaunay;

// typedef double REAL;
typedef K::FT REAL;


struct ConvexPolytope
{
	typedef Eigen::Vector<REAL,2> Point2;
	typedef Eigen::Vector3d Point3;

	const REAL maxHeight=100;

public:
	struct Vertex;
	struct Edge;
	struct Hyperplane;

	map<int, Vertex> vertexs;
	map<int, Edge> edges;
	map<int, Hyperplane> hyperPlanes;

	int hyperPlaneId, vertexId, edgeId;

	//点
	struct Vertex {
		int hplane[3];
		Eigen::Matrix<REAL, 3, 1> g[3];
		REAL w[3];

		Vertex() {}
		Vertex(const Vertex& rhs) {
			for (int i = 0; i < 3; ++i) {
				hplane[i] = rhs.hplane[i];
				g[i]=rhs.g[i];
				w[i]=rhs.w[i];
			}
		};
		Point3 get_point()const {
			Eigen::Matrix<REAL, 3, 1> AinvW;
			Eigen::Matrix<REAL, 3, 3> A;
			Eigen::Matrix<REAL, 3, 1> W;
			for (int i = 0; i < 3; ++i) {
				A(i, 0) = g[i](0);
				A(i, 1) = g[i](1);
				A(i, 2) = g[i](2);
				W(i) = w[i];
			}
			AinvW = A.inverse() * W;
			return Point3(CGAL::to_double(-AinvW(0, 0)),CGAL::to_double(-AinvW(1, 0)),CGAL::to_double(-AinvW(2, 0)));
		}

		void setHyperPlane(const int h1, const int h2, const int h3, map<int, Hyperplane>& hyperPlanes) {
			hplane[0] = h1;
			hplane[1] = h2;
			hplane[2] = h3;
			for (int i=0;i<3;++i) {
				g[i]=hyperPlanes[hplane[i]].g;
				w[i]=hyperPlanes[hplane[i]].w;
			}
		}
		// void setHyperPlane(const int h1, const int h2, const int h3) {
		// 	hplane[0] = h1;
		// 	hplane[1] = h2;
		// 	hplane[2] = h3;
		// }
		Vertex& operator=(const Vertex& rhs) {
			for (int i = 0; i < 3; ++i) {
				hplane[i] = rhs.hplane[i];
				g[i] = rhs.g[i];
				w[i] = rhs.w[i];
			}
			return *this;
		}
		// void initGAndw(map<int, Hyperplane>& hyperPlanes) {
		//
		// }
	};

	//边
	struct Edge {
		int p1, p2;
		int hplanes[2];

		Edge(const Edge& rhs) :p1(rhs.p1), p2(rhs.p2) {
			for (int i = 0; i < 2; ++i)
				hplanes[i] = rhs.hplanes[i];
		}
		Edge(const int _p1, const int _p2, const int _h1, const int _h2) :\
			p1(_p1), p2(_p2) {
			hplanes[0] = _h1;
			hplanes[1] = _h2;
		};
		Edge(const int _p1, const int _p2) :\
			p1(_p1), p2(_p2) {
		};
		Edge() {
		};
	};


struct Hyperplane {
	int sourceIndex;
	Eigen::Matrix<REAL, 3, 1> g;
	REAL w;
	Hyperplane() {}
	Hyperplane& operator=(const Hyperplane& rhs) {
		sourceIndex = rhs.sourceIndex;
		g = rhs.g;
		w = rhs.w;
		return *this;
	}
	Hyperplane(const Hyperplane& rhs) {
		sourceIndex = rhs.sourceIndex;
		g = rhs.g;
		w = rhs.w;
	}
	Hyperplane(REAL d0, REAL d1, REAL d2, REAL d3, int soureceids) : sourceIndex(soureceids) {
		g = Eigen::Matrix<REAL, 3, 1>(d0,d1,d2);
		w = d3;
	}
	Hyperplane(int sourceIds) :sourceIndex(sourceIds) {}

	enum VERTEXOFHYPERPLANE {
		ONPLANE,
		ONUPPERSIDE,
		ONBOTTONSIDE
	};
	REAL zero = REAL();
	VERTEXOFHYPERPLANE onUpperSide(const Vertex& rhs)const {
		Eigen::Matrix<REAL, 4, 4> A_prime;
		A_prime << rhs.g[0](0), rhs.g[0](1), rhs.g[0](2), rhs.w[0],
				   rhs.g[1](0), rhs.g[1](1), rhs.g[1](2), rhs.w[1],
				   rhs.g[2](0), rhs.g[2](1), rhs.g[2](2), rhs.w[2],
				   g(0), g(1), g(2), w;

		Eigen::Matrix<REAL, 3, 3> A;
		A << rhs.g[0](0), rhs.g[0](1), rhs.g[0](2),
			 rhs.g[1](0), rhs.g[1](1), rhs.g[1](2),
			 rhs.g[2](0), rhs.g[2](1), rhs.g[2](2);
		REAL det_A_prime = A_prime.determinant();
		REAL det_A = A.determinant();

		int sign_det_A_prime = (det_A_prime > zero) ? 1 : ((det_A_prime < zero) ? -1 : 0);
		int sign_det_A = (det_A > zero) ? 1 : ((det_A < 0) ? -1 : 0);
		int final_sign = sign_det_A_prime * sign_det_A;

		if (final_sign > 0) { // Example sign adjustment
			return ONBOTTONSIDE;
		} else {
			return ONUPPERSIDE;
		}
	}
};


ConvexPolytope(const ConvexPolytope& rhs) {
	hyperPlanes = rhs.hyperPlanes;
	vertexs = rhs.vertexs;
	edges = rhs.edges;
	hyperPlaneId = rhs.hyperPlaneId;
	vertexId = rhs.vertexId;
	edgeId = rhs.edgeId;
}
ConvexPolytope() {
	//应该是处理到面的程度，切割到。然后利用产生的线，形成新的面。
	for (int i = 1; i <= 8; ++i) {
		vertexs.insert(make_pair(i, Vertex()));
	}
	for (int i = -6; i <= -1; ++i) {
		hyperPlanes.insert(make_pair(i, Hyperplane(i)));
	}
	edges[1] = Edge(1, 2, -1, -5);
	edges[2] = Edge(2, 3, -2, -5);
	edges[3] = Edge(3, 4, -3, -5);
	edges[4] = Edge(1, 4, -5, -4);
	edges[5] = Edge(1, 5, -1, -4);
	edges[6] = Edge(2, 6, -1, -2);
	edges[7] = Edge(3, 7, -2, -3);
	edges[8] = Edge(4, 8, -3, -4);
	edges[9] = Edge(5, 6, -1, -6);
	edges[10] = Edge(6, 7, -2, -6);
	edges[11] = Edge(7, 8, -3, -6);
	edges[12] = Edge(5, 8, -4, -6);


	hyperPlaneId = -1;
	vertexId = 8;
	edgeId = 12;

}
void init() {
	hyperPlanes[-1].g = Eigen::Matrix<REAL, 3, 1>(0,1,0);
	hyperPlanes[-1].w = -maxHeight;

	hyperPlanes[-2].g = Eigen::Matrix<REAL, 3, 1>(1,0,0);
	hyperPlanes[-2].w = -maxHeight;

	hyperPlanes[-3].g = Eigen::Matrix<REAL, 3, 1>(0,-1,0);
	hyperPlanes[-3].w = -maxHeight;

	hyperPlanes[-4].g = Eigen::Matrix<REAL, 3, 1>(-1,0,0);
	hyperPlanes[-4].w = -maxHeight;

	hyperPlanes[-5].g = Eigen::Matrix<REAL, 3, 1>(0,0,-1);
	hyperPlanes[-5].w = -maxHeight;

	hyperPlanes[-6].g = Eigen::Matrix<REAL, 3, 1>(0,0,1);
	hyperPlanes[-6].w = -maxHeight;


	vertexs[1].setHyperPlane(-1, -4, -5,hyperPlanes);
	vertexs[2].setHyperPlane(-1, -2, -5,hyperPlanes);
	vertexs[3].setHyperPlane(-2, -3, -5,hyperPlanes);
	vertexs[4].setHyperPlane(-3, -4, -5,hyperPlanes);
	vertexs[5].setHyperPlane(-1, -4, -6,hyperPlanes);
	vertexs[6].setHyperPlane(-1, -2, -6,hyperPlanes);
	vertexs[7].setHyperPlane(-2, -3, -6,hyperPlanes);
	vertexs[8].setHyperPlane(-3, -4, -6,hyperPlanes);
}

int addHyperPlane(const Hyperplane& plane) {
	set<int> uselessEdges;
	set<int> uselessVertexs;

	map<int, Hyperplane::VERTEXOFHYPERPLANE> vertexsOnPlaneState;
	for (const auto& vertex : vertexs) {
		auto state = plane.onUpperSide(vertex.second);
		if (state == Hyperplane::VERTEXOFHYPERPLANE::ONUPPERSIDE) {
			uselessVertexs.insert(vertex.first);
		}
		vertexsOnPlaneState.insert(make_pair(vertex.first, state));
	}

	int status_code = 0;
	if (uselessVertexs.empty()) {
		return -1;
	}
	hyperPlaneId = plane.sourceIndex;
	hyperPlanes.insert(make_pair(hyperPlaneId, Hyperplane(plane)));

	map<int, set<int>> verticesInFacet;

	for (auto& curEdge:edges) {
		auto& e=curEdge.second;
		auto useEnd1Useless = vertexsOnPlaneState[e.p1];
		auto useEnd2Useless = vertexsOnPlaneState[e.p2];

		if (useEnd1Useless != ConvexPolytope::Hyperplane::VERTEXOFHYPERPLANE::ONBOTTONSIDE && useEnd2Useless != ConvexPolytope::Hyperplane::VERTEXOFHYPERPLANE::ONBOTTONSIDE) {
			uselessEdges.insert(curEdge.first);
			continue;
		}

		if (useEnd1Useless != ConvexPolytope::Hyperplane::VERTEXOFHYPERPLANE::ONBOTTONSIDE) {
			Vertex vertex_new;
			vertex_new.setHyperPlane(e.hplanes[0],e.hplanes[1],hyperPlaneId,hyperPlanes);
			vertexs[++vertexId] = vertex_new;
			e.p1 = vertexId;
			verticesInFacet[e.hplanes[0]].insert(vertexId);
			verticesInFacet[e.hplanes[1]].insert(vertexId);

		}
		else if (useEnd2Useless != ConvexPolytope::Hyperplane::VERTEXOFHYPERPLANE::ONBOTTONSIDE) {
			Vertex vertex_new;
			vertex_new.setHyperPlane(e.hplanes[0],e.hplanes[1],hyperPlaneId,hyperPlanes);

			vertexs[++vertexId] = Vertex(vertex_new);
			e.p2 = vertexId;

			verticesInFacet[e.hplanes[0]].insert(vertexId);
			verticesInFacet[e.hplanes[1]].insert(vertexId);
		}
	}

	for (auto& edg : uselessEdges) {
		edges.erase(edg);
	}

	for (auto& usef : uselessVertexs) {
		vertexs.erase(usef);
	}

	map<int, set<int>> edgesInHyperPlane;

	for (const pair<int, set<int>>& edge_pair : verticesInFacet) {

		if (edge_pair.second.size() != 2) {
			status_code = -1024;
			std::cout << "edge vertex not 2" << endl;
			std::cout << edge_pair.first << endl;
			std::cout << edge_pair.second.size() << endl;
		}

		Edge newEdge(*edge_pair.second.begin(), *edge_pair.second.rbegin());

		newEdge.hplanes[0] = edge_pair.first;
		newEdge.hplanes[1] = hyperPlaneId;

		edges[++edgeId] = newEdge;
	}
	return status_code;
}
vector<int> constructFace(vector<pair<int,int>>& edges) {
	if (edges.empty()) {
		return {};
	}
	unordered_map<int, vector<int>> graph;
	for (const auto& edge : edges) {
		graph[edge.first].push_back(edge.second);
		graph[edge.second].push_back(edge.first);
	}
	vector<int> result;
	int start = edges[0].first;
	int current = start;
	int previous = -1000;
	do {
		result.push_back(current);
		int next = -1000;
		for (int neighbor : graph[current]) {
			if (neighbor != previous) {
				next = neighbor;
				break;
			}
		}
		previous = current;
		current = next;

	} while (current != start && current != -1000);
	return result;
}
	vector<pair<int,vector<Eigen::Vector3d>>> GetConvex()
{
	map<int,vector<pair<int,int>>> regions;
	for (auto cur_edge : edges)
	{
		auto e = cur_edge.second;
		regions[e.hplanes[0]].push_back(make_pair(e.p1,e.p2));
		regions[e.hplanes[1]].push_back(make_pair(e.p2,e.p1));
	}
	vector<pair<int,vector<Eigen::Vector3d>>> result;

	for (auto& region : regions) {
		auto face=constructFace(region.second);
		vector<Eigen::Vector3d> face_with_points;
		for (auto pid:face) {
			face_with_points.push_back(vertexs[pid].get_point());
		}
		result.emplace_back(make_pair(region.first,face_with_points));
	}
	return result;
}
	/**
 * 将多面体数据输出为OBJ文件（简化版，无去重）
 * @param polyhedron 多面体数据结构：vector<pair<面ID, vector<该面的顶点>>>
 * @param outputPath 输出文件路径
 * @return 是否成功输出
 */
	bool exportPolyhedronToOBJ(
		const string& outputPath
	) {
	auto polyhedron = GetConvex();
	// 打开输出文件
	ofstream objFile(outputPath);
	if (!objFile.is_open()) {
		cerr << "错误: 无法打开文件 " << outputPath << endl;
		return false;
	}

	// 写入文件头
	objFile << "# OBJ文件 - 多面体导出" << endl;
	objFile << "# 面数: " << polyhedron.size() << endl;
	objFile << endl;

	// 计算总顶点数
	int totalVertices = 0;
	for (const auto& face : polyhedron) {
		totalVertices += face.second.size();
	}

	// 输出所有顶点
	objFile << "# 顶点坐标 (总数: " << totalVertices << ")" << endl;
	objFile << fixed << setprecision(6);

	for (const auto& face : polyhedron) {
		objFile << "# 面ID: " << face.first << endl;
		for (const auto& vertex : face.second) {
			objFile << "v " << vertex.x() << " " << vertex.y() << " " << vertex.z() << endl;
		}
	}

	objFile << endl;

	// 输出所有面
	objFile << "# 面定义" << endl;
	int currentIndex = 1; // OBJ索引从1开始

	for (const auto& face : polyhedron) {
		int faceId = face.first;
		int vertexCount = face.second.size();

		if (vertexCount < 3) {
			cerr << "警告: 面ID " << faceId << " 顶点数不足3个，跳过" << endl;
			currentIndex += vertexCount;
			continue;
		}

		// 写入面
		objFile << "f";
		for (int i = 0; i < vertexCount; ++i) {
			objFile << " " << (currentIndex + i);
		}
		objFile << endl;

		currentIndex += vertexCount;
	}

	objFile.close();

	cout << "成功导出OBJ文件: " << outputPath << endl;
	cout << "- 总顶点数: " << totalVertices << endl;
	cout << "- 面数: " << polyhedron.size() << endl;

	return true;
}

};

#include <random>

struct PlaneParams {
	double a, b, c, d;

	// 构造函数
	PlaneParams(double a, double b, double c, double d) : a(a), b(b), c(c), d(d) {}

	// 打印平面方程
	void printEquation() const {
		std::cout << "平面方程: " << a << "x + " << b << "y + " << c << "z + " << d << " = 0" << std::endl;
	}

	// 验证点(0,0,0)是否在正侧
	bool isOriginOnPositiveSide() const {
		// 对于点(0,0,0)，有符号距离 = d / ||(a,b,c)||
		double norm = std::sqrt(a*a + b*b + c*c);
		if (norm == 0) return false; // 无效平面
		return d > 0;
	}
};

/**
 * 随机生成平面参数，使得点(0,0,0)位于平面正侧
 * @param seed 随机种子，默认使用当前时间
 * @return PlaneParams 包含a,b,c,d四个参数的结构体
 */
PlaneParams generateRandomPlane(unsigned int seed = 0) {
	// 设置随机数生成器
	std::random_device rd;
	std::mt19937 gen(seed == 0 ? rd() : seed);
	std::uniform_real_distribution<double> dis(-5, 5.0);

	double a, b, c, d;

	do {
		// 随机生成a, b, c
		a = dis(gen);
		b = dis(gen);
		c = dis(gen);

		// 确保至少一个系数不为0（避免法向量为零向量）
	} while (a == 0.0 && b == 0.0 && c == 0.0);

	// 生成d，确保d > 0以使原点在正侧
	// 为了避免d = 0的边界情况，确保d在(0, 5]范围内
	std::uniform_real_distribution<double> d_dis(0.001, 5.0);
	d = d_dis(gen);

	return PlaneParams(a, b, c, d);
}

int main() {
	ConvexPolytope polytope;
	polytope.init();
	for (int i=0;i<100000;++i) {
		auto plane = generateRandomPlane();
		// cout<<"plane "<<plane.a<<" "<<plane.b<<" "<<plane.c<<" "<<plane.d<<endl;
		polytope.addHyperPlane(ConvexPolytope::Hyperplane(plane.a,plane.b,plane.c,plane.d,2*i));
		// polytope.addHyperPlane(ConvexPolytope::Hyperplane(plane.a,plane.b,plane.c,plane.d,2*i+1));

	}

	polytope.exportPolyhedronToOBJ("/Users/pengfei/CLionProjects/ConvexCut/data/convex.obj");
	return 0;
}


