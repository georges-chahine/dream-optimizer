﻿#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include "yaml-cpp/yaml.h"
#include <sstream>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "tf_conversions/tf_eigen.h"
#include <pcl_ros/point_cloud.h>
#include <pcl_ros/transforms.h>
#include "g2o/solvers/cholmod/linear_solver_cholmod.h"
#include "g2o/solvers/dense/linear_solver_dense.h"
#include "g2o/stuff/sampler.h"

#include <pcl/common/transforms.h>
//#include "custom_types/vertex_pose.h"
//#include "custom_types/edge_pose_pose.h"
#include "g2o/types/sba/types_six_dof_expmap.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include <g2o/core/sparse_optimizer.h>
#include <g2o/core/block_solver.h>
#include "g2o/core/solver.h"
#include "edge_se3.h"

using namespace std;
using namespace Eigen;
using namespace g2o;
template<typename M>


M load_csv (const std::string & path) {
    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<double> values;
    uint rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }
    return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, RowMajor>>(values.data(), rows, values.size()/rows);
}


void addEdge(const Eigen::Isometry3d & a_T_b, const int id_a, const int id_b,
             g2o::OptimizableGraph* graph_ptr, Eigen::Matrix<double, 6, 6> Lambda)
{
    //EdgePosePose * e = new EdgePosePose;
    Eigen::Matrix3d rot=a_T_b.matrix().block(0,0,3,3);
    Eigen:: Quaterniond q(rot);
    // q.setIdentity();


    Vector3d trans(a_T_b(0,3),a_T_b(1,3),a_T_b(2,3));

    g2o::SE3Quat pose(q,trans);


    // g2o::VertexSE3Expmap * v_se3= new g2o::VertexSE3Expmap();
    g2o::EdgeSE3* v_se3 = new g2o::EdgeSE3;

    // retrieve vertex pointers from graph with id's
    g2o::OptimizableGraph::Vertex * pose_a_vertex
            = dynamic_cast<g2o::OptimizableGraph::Vertex*>
            (graph_ptr->vertices()[id_a]);

    g2o::OptimizableGraph::Vertex * pose_b_vertex
            = dynamic_cast<g2o::OptimizableGraph::Vertex*>
            (graph_ptr->vertices()[id_b]);

    // error check vertices and associate them with the edge
    assert(pose_a_vertex!=NULL);
    assert(pose_a_vertex->dimension() == 6);
    // e->vertices()[0] = pose_a_vertex;
    v_se3->vertices()[0] = pose_a_vertex;
    assert(pose_b_vertex!=NULL);
    assert(pose_b_vertex->dimension() == 6);
    //e->vertices()[1] = pose_b_vertex;
    v_se3->vertices()[1] = pose_b_vertex;
    // add information matrix
    //Eigen::Matrix<double, 6, 6> Lambda;
    //Lambda.setIdentity();

    // set the observation and imformation matrix

    v_se3->setMeasurement(pose);
    // e->setMeasurement(a_T_b);
    // e->information() = Lambda;
    v_se3->information() = Lambda;
    // finally add the edge to the graph
    if(!graph_ptr->addEdge(v_se3))
    {
        assert(false);
    }
}


void addPosesToGraph(Eigen::Isometry3d pose_estimate, int id, g2o::OptimizableGraph* graph_ptr)
{

    Eigen::Matrix3d rot=pose_estimate.matrix().block(0,0,3,3);
    Eigen:: Quaterniond q(rot);
    // q.setIdentity();


    Vector3d trans(pose_estimate(0,3),pose_estimate(1,3),pose_estimate(2,3));

    g2o::SE3Quat pose(q,trans);
    g2o::VertexSE3 * v_se3= new g2o::VertexSE3();



    // VertexPose * v = new VertexPose();
    v_se3->setEstimate(pose);
    v_se3->setId(id);
    v_se3->setFixed(false);


    if (id==0)
    {
        std::cout<<"id is "<<id<<std::endl;
        v_se3->setFixed(true);
    }


    graph_ptr->addVertex(v_se3);


    /*

    VertexPose * v = new VertexPose();
    v->setEstimate(pose_estimate);
    v->setId(id);
    v->setFixed(false);


    if (id==0)
    {
        std::cout<<"id is "<<id<<std::endl;
        v->setFixed(true);
    }


    graph_ptr->addVertex(v);
*/
}

int main(int argc, char *argv[]){
    // initialize graph

    g2o::SparseOptimizer graph;
    int dir;
    graph.setVerbose(true);
    YAML::Node config = YAML::LoadFile("../config.yaml");

    std::string graphConstraintsStr = config["graphConstraints"].as<std::string>();
    std::string graphInitValuesStr= config["graphInitValues"].as<std::string>();
    std::string pathOut = config["pathOut"].as<std::string>();
    std::string maxIterationsStr = config["maxIterations"].as<std::string>();
    int maxIterations=std::stoi(maxIterationsStr);

    dir=mkdir (pathOut.c_str(),S_IRWXU);
    MatrixXd graphConstraints = load_csv<MatrixXd>(graphConstraintsStr);
    MatrixXd graphInitValues = load_csv<MatrixXd>(graphInitValuesStr);

    for (int i=0; i<graphInitValues.rows(); i++){


        int poseId=graphInitValues(i,0);


        Eigen::Quaterniond q(graphInitValues(i,7),graphInitValues(i,4),graphInitValues(i,5),graphInitValues(i,6));

        Eigen::Matrix3d rot(q);

        Eigen::Isometry3d T=Eigen::Isometry3d::Identity();
        T.matrix().block(0,0,3,3)=rot;

        T(0,3)=graphInitValues(i,1);
        T(1,3)=graphInitValues(i,2);
        T(2,3)=graphInitValues(i,3);

        addPosesToGraph(T, poseId, &graph);

    }

    for (int i=0; i<graphConstraints.rows(); i++){



        const int pose_a_id=graphConstraints(i,0);
        const int pose_b_id=graphConstraints(i,1);


        Eigen::Quaterniond q(graphConstraints(i,8),graphConstraints(i,5),graphConstraints(i,6),graphConstraints(i,7));

        Eigen::Matrix3d rot(q);

        Eigen::Isometry3d T=Eigen::Isometry3d::Identity();
        T.matrix().block(0,0,3,3)=rot;
        //std::cout<<"T was \n"<<T.matrix()<<std::endl;
        T(0,3)=graphConstraints(i,2);
        T(1,3)=graphConstraints(i,3);
        T(2,3)=graphConstraints(i,4);

        //std::cout<<"T is \n"<<T.matrix()<<std::endl;
        if (T.matrix().isIdentity(0.01)){
            std::cout<<"filtering Identity"<<std::endl;
            continue;

        }

        Eigen::Matrix<double, 6, 6> Lambda;
        Lambda.setIdentity();
        Lambda(0,0)=Lambda(1,1)=Lambda(2,2)=Lambda(3,3)=Lambda(4,4)=Lambda(5,5)=1/graphConstraints(i,9);

        addEdge(T, pose_a_id, pose_b_id, &graph, Lambda);



    }



    bool DENSE=false;
    // finished populating graph, export and quit
    /* g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    if (DENSE) {
        linearSolver= new g2o::LinearSolverDense<g2o::BlockSolver_6_3::PoseMatrixType>();
    } else {
        linearSolver= new g2o::LinearSolverCholmod<g2o::BlockSolver_6_3::PoseMatrixType>();
    }
    */
    std::unique_ptr<g2o::BlockSolver_6_3::LinearSolverType> linearSolver (new g2o::LinearSolverCholmod<g2o::BlockSolver_6_3::PoseMatrixType>());


    // linearSolver = new g2o::LinearSolverCholmod<g2o::BlockSolver_6_3::PoseMatrixType>();
    // g2o::BlockSolver_6_3 * solver_ptr= new g2o::BlockSolver_6_3(linearSolver);

    std::unique_ptr<g2o::BlockSolver_6_3> solver_ptr (new g2o::BlockSolver_6_3(std::move(linearSolver)));


    //  g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    g2o::OptimizationAlgorithmLevenberg * solver = new g2o::OptimizationAlgorithmLevenberg(std::move(solver_ptr));


    //  g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(
    //      g2o::make_unique<BlockSolverX>(g2o::make_unique<LinearSolverCholmod<g2o::BlockSolverX::PoseMatrixType>>()));
    graph.setAlgorithm(solver);
    graph.initializeOptimization();

    std::cout << "Saving Graph to example_g2o.g2o...\n";
    graph.save("example_g2o.g2o");
    graph.optimize(maxIterations);
     graph.save("example_g2o_opt.g2o");
    std::cout << "Finished\n";
    std::cout << "rosrun g2o_viewer g2o_viewer example_g2o.g2o\n";


    return 0;
}



