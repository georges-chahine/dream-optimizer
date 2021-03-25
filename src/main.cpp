#include <stdio.h>
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
#include "g2o/core/robust_kernel_impl.h"
#include <pcl/common/transforms.h>
//#include "custom_types/vertex_pose.h"
//#include "custom_types/edge_pose_pose.h"
#include "g2o/types/sba/types_six_dof_expmap.h"
#include "g2o/config.h"
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
            //  std::cout<<cell<<std::endl;
            values.push_back(std::stod(cell));
        }
        ++rows;
    }
    return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, RowMajor>>(values.data(), rows, values.size()/rows);
}


unsigned int returnIndex(unsigned int i, unsigned int j, unsigned int maxKF){
    // maxKF starts at 0
    int idx= (maxKF+1) * j + i;
    return (idx);
}

Eigen::MatrixXd load_g2o (std::string path, MatrixXd transforms0) {

    std::string line;
    std::ifstream inFile(path);

    std::string str;
    double x,y,z,qx,qy,qz,qw;
    std::vector<double> xV,yV,zV,qxV,qyV,qzV,qwV;
    unsigned int idxNbr;
    std::vector<unsigned int> idxNbrV;

    time_t t1;
    float temprature;

    int maxKF=transforms0(transforms0.rows()-1,0);
    std::cout<<"maxKF is " <<maxKF<<std::endl;
    while (inFile >> str >>idxNbr >> x >>y >>z >>qx>>qy>>qz>>qw)
    {
        if (str!="VERTEX_SE3:QUAT"){break;};
        xV.push_back(x);
        yV.push_back(y);
        zV.push_back(z);
        qxV.push_back(qx);
        qyV.push_back(qy);
        qzV.push_back(qz);
        qwV.push_back(qw);
        idxNbrV.push_back(idxNbr);
        //cout << x << " " << y << endl;

    }
    vector<vector<int>> idxNbrV2;

    for (int i=0; i<=idxNbrV.back(); i++){

        // std::cout<<"current KF " <<i%(maxKF+1)<<" current T "<< int (i)/(maxKF+1) <<std::endl;
        //std::cout<<"current KF " <<int(i/maxKF)<<std::endl;
        int KF=i%(maxKF+1);
        int T=int (i/(maxKF+1));
        idxNbrV2.push_back(std::vector<int> {KF, T, i});

    }

    Eigen::MatrixXd transforms(idxNbrV2.size(),10);


    //Eigen::MatrixXd transforms1(1,idxNbrV.size());

    for (int i=0; i<transforms.rows(); i++){


        transforms(i,0)=idxNbrV2[i][0];
        transforms(i,1)=idxNbrV2[i][1];
        bool found=false;
        for (int j=0; j<idxNbrV.size(); j++){
            if (idxNbrV[j]==idxNbrV2[i][2]){
                transforms(i,2)=xV[j];
                transforms(i,3)=yV[j];
                transforms(i,4)=zV[j];
                transforms(i,5)=qxV[j];
                transforms(i,6)=qyV[j];
                transforms(i,7)=qzV[j];
                transforms(i,8)=qwV[j];
                found=true;
            }
        }

        if (!found){
            std::cout<<"couldn't find keyframe/timeperiod matches"<<std::endl;
            transforms(i,2)=0;
            transforms(i,3)=0;
            transforms(i,4)=0;
            transforms(i,5)=0;
            transforms(i,6)=0;
            transforms(i,7)=0;
            transforms(i,8)=1;
        }

        for (int j=0; j<transforms0.rows(); j++){

            if (transforms0(j,0)==idxNbrV2[i][0] && transforms0(j,1)==idxNbrV2[i][1]){

                transforms(i,9)=transforms0(j,9);

            }

        }


    }

    inFile.close();

    std::cout<<transforms<<std::endl;

    return transforms;
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
    //g2o::RobustKernelGemanMcClure* rk = new g2o::RobustKernelGemanMcClure;
    g2o::RobustKernelCauchy* rk = new g2o::RobustKernelCauchy;
    v_se3->setRobustKernel(rk);
    if(!graph_ptr->addEdge(v_se3))
    {
        assert(false);
    }
}


void addOverlappingEdge(int mode, unsigned int skip, int maxKF, int searchKFLimit, std::vector<std::string> autoMatchDir, std::vector<std::vector<std::string>> csvFiles, std::string autoNameTemplate, Eigen::MatrixXd transforms, g2o::OptimizableGraph& graph, Eigen::Matrix<double, 6, 6> Lambda){



    for (unsigned int i=0; i<transforms.rows();i=i+skip){
        int current_KF=transforms(i,0);
        int current_T=transforms(i,1);
        if (transforms(i,2)==0&& transforms(i,3)==0&&transforms(i,4)==0&&transforms(i,5)==0){
            std::cout<<"keyframe not found in graph, skipping KF,T: "<<current_KF<<","<<current_T<<std::endl;
            continue;

        }
        std::string strCsvName=autoMatchDir[current_T]+"/"+autoNameTemplate+std::to_string(current_KF)+".csv";
        MatrixXd current_csv = load_csv<MatrixXd>(strCsvName);
        double currentStamp=current_csv(0,0);
        unsigned int sz=0;
        if(mode==0){
            sz=0;
        }
        if(mode==1){
            sz=current_csv.rows()/2;

        }
        if(mode==2){
            sz=current_csv.rows()-1;
        }


        Eigen::Isometry3d M000=Eigen::Isometry3d::Identity();
        Eigen::Quaterniond q000(current_csv(0,7), current_csv(0,4), current_csv(0,5), current_csv(0,6) );
        M000(0,3)=current_csv(0,1);
        M000(1,3)=current_csv(0,2);
        M000(2,3)=current_csv(0,3);

        currentStamp=current_csv(sz,0);
        Eigen::Isometry3d M00=Eigen::Isometry3d::Identity();
        Eigen::Quaterniond q00(current_csv(sz,7), current_csv(sz,4), current_csv(sz,5), current_csv(sz,6) );
        M00(0,3)=current_csv(sz,1);
        M00(1,3)=current_csv(sz,2);
        M00(2,3)=current_csv(sz,3);


        int currentNode=returnIndex(current_KF,current_T, maxKF);




        for (int j=0; j<csvFiles[current_T].size(); j++){
            // std::cout<<"j, current_KF, size "<<j<<","<<current_T<<","<<csvFiles[current_T].size()<<std::endl;
            if ( abs(current_KF-j)>searchKFLimit || j==current_KF)
            {
                continue;
            }

            else
            {
                MatrixXd temp_csv = load_csv<MatrixXd>(csvFiles[current_T][j]);
                int secondNode=returnIndex(j,current_T, maxKF);
                //std::cout<< "secondnode is "<<secondNode<<std::endl;

                Eigen::Isometry3d M0=Eigen::Isometry3d::Identity();

                unsigned int sz2=0;

                M0.matrix()(0,3)=temp_csv(sz2,1);//-temp_csv(0,1);
                M0.matrix()(1,3)=temp_csv(sz2,2);//-temp_csv(0,2);
                M0.matrix()(2,3)=temp_csv(sz2,3);//-temp_csv(0,3);

                Eigen::Quaterniond qd0(temp_csv(sz2,7),temp_csv(sz2,4),temp_csv(sz2,5),temp_csv(sz2,6));
                Eigen::Matrix3d rot0(qd0);
                M0.matrix().block(0,0,3,3)=rot0;

                Eigen::Isometry3d M0Opt=Eigen::Isometry3d::Identity();
                bool found=false;
                for (int k=0; k<transforms.rows(); k++){

                    if (transforms(k,0)==j && transforms(k,1)==current_T    ){


                        M0Opt.matrix()(0,3)=transforms(k,2);
                        M0Opt.matrix()(1,3)=transforms(k,3);
                        M0Opt.matrix()(2,3)=transforms(k,4);
                        Eigen::Quaterniond q0Opt(transforms(k,8),transforms(k,5),transforms(k,6),transforms(k,7));
                        Eigen::Matrix3d R0Opt(q0Opt);
                        M0Opt.matrix().block(0,0,3,3)=R0Opt;
                        found=true;
                        break;
                    }


                }
                if (!found){

                    std::cout<<"lookup error to debug"<<std::endl;

                }
                if (M0Opt.matrix().isIdentity())
                {
                    std::cout<<"keyframe match not found in graph, skipping KF,T: "<<j<<","<<current_T<<std::endl;
                    break;
                }

                for (int k=1; k<temp_csv.rows(); k++){

                    if (temp_csv(k,0)==currentStamp){

                        //idCounter++;

                        Eigen::Isometry3d M=Eigen::Isometry3d::Identity();
                        M.matrix()(0,3)=temp_csv(k,1);//-temp_csv(0,1);
                        M.matrix()(1,3)=temp_csv(k,2);//-temp_csv(0,2);
                        M.matrix()(2,3)=temp_csv(k,3);//-temp_csv(0,3);

                        Eigen::Quaterniond qd(temp_csv(k,7),temp_csv(k,4),temp_csv(k,5),temp_csv(k,6));
                        Eigen::Matrix3d rot(qd);
                        M.matrix().block(0,0,3,3)=rot;

                        // if (M0Opt.matrix().isIdentity(0.0001)){continue;}
                        //  addPosesToGraph(M0Opt*M0.inverse()*M, idCounter, &graph);
                        //  addEdge(M00, idCounter, currentNode, &graph, Lambda);

                        //  Eigen::Isometry3d Mf=(M000.inverse() * M00)*(M.inverse()*M0);
                        Eigen::Isometry3d Mf=(M000.inverse()*M0);

                        //   M000.inverse()*M
                        addEdge(Mf, currentNode, secondNode, &graph, Lambda);
                    }

                }

            }

        }

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
        v_se3->setFixed(false);
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
    std::string transformsFile= config["transformsCsv"].as<std::string>();
    std::vector<std::string> autoMatchDir = config["autoMatchDir"].as<std::vector<std::string>>();
    int maxIterations=std::stoi(maxIterationsStr);
    std::string autoNameTemplate = config["autoNameTemplate"].as<std::string>();

    dir=mkdir (pathOut.c_str(),S_IRWXU);
    MatrixXd graphConstraints = load_csv<MatrixXd>(graphConstraintsStr);
    MatrixXd graphInitValues = load_csv<MatrixXd>(graphInitValuesStr);
    int maxPoseId=0;
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

        if (poseId>maxPoseId)
        {
            maxPoseId=poseId;
        }
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

        Lambda(0,0)=Lambda(1,1)=Lambda(2,2)=1/graphConstraints(i,9);
        Lambda(3,3)=Lambda(4,4)=Lambda(5,5)=1/graphConstraints(i,9);

        //Lambda(0,0)=Lambda(1,1)=Lambda(2,2)=Lambda(3,3)=Lambda(4,4)=Lambda(5,5)=1/1;

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
    // g2o::OptimizableGraph::VertexContainer V=graph.activeVertices();

    //    VertexSE3* firstRobotPose = dynamic_cast<VertexSE3*>(graph.vertex(0));

    //   double* estimate;
    // firstRobotPose->getEstimateData(estimate);
    //std::cout<<estimate[0]<<" "<<estimate[1]<<" "<<estimate[2]<<" "<<estimate[3]<<" "<<estimate[4]<<" "<<estimate[5]<<" "<<estimate[6]<<std::endl;


    std::string g2oFile="example_g2o_opt.g2o";
    graph.save(g2oFile.c_str());
    std::cout << "Finished\n";
    std::cout << "rosrun g2o_viewer g2o_viewer example_g2o.g2o\n";


    int searchKFLimit=5;
    int dirNumber=autoMatchDir.size();
    MatrixXd transforms;
    MatrixXd transforms0 = load_csv<MatrixXd>(transformsFile);
    int maxKF=transforms0(transforms0.rows()-1,0);
    transforms = load_g2o(g2oFile, transforms0);



    std::vector<std::vector<std::string>> csvFiles;
    for (unsigned int i=0;  i<dirNumber; i++ )
    {
        std::string currentPath=pathOut+"/"+std::to_string(i)+"/";
        dir=mkdir (currentPath.c_str(),S_IRWXU);

        int j=0;

        std::vector<std::string> csvFile;
        while(true){
            std::string strCsvName=autoMatchDir[i]+"/"+autoNameTemplate+std::to_string(j)+".csv";

            std::ifstream fCsv(strCsvName.c_str());


            if (!fCsv.good())
            {
                break;
            }
            //    std::cout<<"string is \n"<<strPcdName<<std::endl;

            csvFile.push_back(strCsvName);
            j++;

        }

        csvFiles.push_back(csvFile);
    }

    // int idCounter=maxPoseId;


    Eigen::Matrix<double, 6, 6> Lambda;
    Lambda.setIdentity();
    Lambda(0,0)=Lambda(1,1)=Lambda(2,2)=Lambda(3,3)=Lambda(4,4)=Lambda(5,5)=1/2;
    //Lambda(3,3)=Lambda(4,4)=Lambda(5,5)=1/2;
    //Lambda(3,3)=Lambda(4,4)=Lambda(5,5)=1/2;

    //          int mode=0;
    //   const unsigned int skip=1;

    //   addOverlappingEdge(mode, skip, maxKF, searchKFLimit, autoMatchDir, csvFiles, autoNameTemplate, transforms, graph, Lambda);
    // mode=0;
    // addOverlappingEdge(mode, skip, maxKF, searchKFLimit, autoMatchDir, csvFiles, autoNameTemplate, transforms, graph, Lambda);
    // mode=2;
    // addOverlappingEdge(mode, skip, maxKF, searchKFLimit, autoMatchDir, csvFiles, autoNameTemplate, transforms, graph, Lambda);
    /*
    graph.initializeOptimization();

    std::cout << "Saving Graph to example_g2o.g2o...\n";
    graph.save("example_g2o2.g2o");
    graph.optimize(maxIterations);
    std::string g2oFile2="example_g2o_opt2.g2o";
    graph.save(g2oFile2.c_str());
    std::cout << "Finished\n";
    // g2o::OptimizableGraph::VertexContainer V=graph.activeVertices();
*/
    return 0;
}
