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
#include <g2o/core/base_unary_edge.h>
#include "g2o/core/solver.h"
#include <g2o/types/slam3d/edge_se3_xyzprior.h>
//#include <g2o/types/slam3d/edge_se3_prior.h>
#include <g2o/types/slam3d/edge_se3.h>
//#include "edge_se3.h"
//#include "edge_se3_xyzprior.h"

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
template<typename N>
N load_xyz (const std::string & path) {
    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<double> values;
    uint rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ' ')) {
            //  std::cout<<cell<<std::endl;
            values.push_back(std::stod(cell));
        }
        ++rows;
    }
    return Map<const Matrix<typename N::Scalar, N::RowsAtCompileTime, N::ColsAtCompileTime, RowMajor>>(values.data(), rows, values.size()/rows);
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

        if (str=="PARAMS_SE3OFFSET"){continue;};
        if (str=="FIX"){continue;}
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

    //std::cout<<transforms<<std::endl;

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


void addXYZEdgeSE3(const Eigen::Vector3d pos_new_z, const int id_a,
                   g2o::OptimizableGraph* graph_ptr, Eigen::Matrix<double, 3, 3> Lambda)
{
    g2o::EdgeSE3XYZPrior* v_se3 = new g2o::EdgeSE3XYZPrior;
    //g2o::EdgeSE3Prior* v_se3_2 = new g2o::EdgeSE3Prior;
    //std::cout<<"id a is " <<id_a<<std::endl;
    // retrieve vertex pointers from graph with id's
    //g2o::OptimizableGraph::Vertex * pose_a_vertex
    //       = dynamic_cast<g2o::OptimizableGraph::Vertex*>
    //      (graph_ptr->vertices()[id_a]);

    g2o::VertexSE3* v1 = (g2o::VertexSE3*)graph_ptr->vertex(id_a);

    // error check vertices and associate them with the edge
    assert(v1!=NULL);
    assert(v1->dimension() == 6);
    // e->vertices()[0] = pose_a_vertex;
    //v_se3->vertices()[0] = pose_a_vertex;
    v_se3->setVertex(0, v1);
    //v_se3_2->setVertex(0, v1);
    //v_se3->setId(id_a);
    v_se3->setParameterId(0, 0);

    //e->vertices()[1] = pose_b_vertex;
    //v_se3->vertices()[1] = pose_b_vertex;
    // add information matrix
    //Eigen::Matrix<double, 6, 6> Lambda;
    //Lambda.setIdentity();

    // set the observation and imformation matrix
    //Eigen::Vector3d pos_new_z2( 10.0, 10.0, 5.0);
    //Eigen::Matrix3d pos_inf_new_z = Eigen::Matrix3d::Zero();
    v_se3->setMeasurement(pos_new_z);
    //Eigen::Isometry3d a;
    //a.setIdentity();
    //v_se3_2->setMeasurement(a);
    //pos_inf_new_z(0, 0) = 0;
    //pos_inf_new_z(1, 1) = 0;
    //pos_inf_new_z(2, 2) = 1E3;
    // e->setMeasurement(a_T_b);
    // e->information() = Lambda;
    //v_se3->information() = Lambda;
    // std::cout<<"lambda is "<<Lambda<<std::endl;
    v_se3->setInformation(Lambda);

    //  Eigen::Matrix<double, 6, 6> Lambda_2;
    //  Lambda_2.setIdentity();
    // v_se3_2->setInformation(Lambda_2);
    //v_se3->setInformation(pos_inf_new_z);
    // v_se3_2->setParameterId(0, 0);
    // finally add the edge to the graph
    //g2o::RobustKernelGemanMcClure* rk = new g2o::RobustKernelGemanMcClure;
    g2o::RobustKernelCauchy* rk = new g2o::RobustKernelCauchy;
    v_se3->setRobustKernel(rk);
    if(!graph_ptr->addEdge(v_se3))
    {
        assert(false);
    }
}



void addDEMconstraints(double DEM_x0, double DEM_y0, MatrixXd& DEM, unsigned int skip, int maxKF, std::vector<std::string> autoMatchDir, std::string autoNameTemplate, Eigen::MatrixXd transforms, g2o::OptimizableGraph& graph, Eigen::Matrix<double, 3, 3> Lambda){

    std::ofstream Zdata;
    Zdata.open("DEM_Z.csv");

    double zInit=0;
    double zRef=0;
    for (unsigned int i=0; i<transforms.rows();i=i+skip){
        int current_KF=transforms(i,0);
        int current_T=transforms(i,1);
        std::cout<<"KF,T: "<<current_KF<<","<<current_T<<std::endl;
        if (transforms(i,2)==0&& transforms(i,3)==0&&transforms(i,4)==0&&transforms(i,5)==0){
            std::cout<<"keyframe not found in graph, skipping KF,T: "<<current_KF<<","<<current_T<<std::endl;
            continue;

        }
        if(current_T!=0){continue;}

        if (current_KF%8!=0){continue;}

        std::string strCsvName=autoMatchDir[current_T]+"/"+autoNameTemplate+std::to_string(current_KF)+".csv";
        MatrixXd current_csv = load_csv<MatrixXd>(strCsvName);


        double x0=current_csv(0,1);
        double y0=current_csv(0,2);
        double z0=current_csv(0,3);


        std::vector<double> distances;
        std::vector<unsigned int> indices;
        //std::cout<<"DEM rows are: "<<DEM.rows()<<std::endl;
        //std::cout<<"DEM cols are: "<<DEM.cols()<<std::endl;
        // std::cout<<"DEM_x0 is: "<<DEM_x0<<std::endl;
        // std::cout<< "x0 is "<<x0<<std::endl;
        for (unsigned int j=0; j<DEM.rows();j=j+1){

            double currentX_DEM=DEM(j,1)+DEM_x0-DEM(0,1);
            double currentY_DEM=DEM(j,0)+DEM_y0-DEM(0,0);

            if ( sqrt(pow(x0-currentX_DEM,2))>0.5 || sqrt(pow(y0-currentY_DEM,2))>0.5 ){

                continue;
            }
            indices.push_back(j);
            double D=sqrt(pow(currentX_DEM-x0,2)+pow(currentY_DEM-y0,2));
            //std::cout<<"D is "<<D<<std::endl;
            distances.push_back(D);

        }
        double lowestD=9999;
        unsigned int chosenIdx;
        if (distances.size()==0){
            std::cout<<"distances size is: "<<distances.size()<<std::endl;
            std::cout<<"FATAL: DEM query failed, results might be unreliable"<<std::endl;
        }
        for (unsigned int j=0; j<distances.size(); j++){

            if (distances[j]<lowestD)
            {
                lowestD=distances[j];
                chosenIdx=indices[j];

            }


        }
        int currentNode=returnIndex(current_KF,current_T, maxKF);
        double zConstraint=-DEM(chosenIdx,2); //-DEM(0,2);

        if (currentNode==0)
        {
            zInit=zConstraint;
            zRef=z0;


        }


        zConstraint=zConstraint-zInit;

        Vector3d T( 0, 0, zConstraint);
        std::cout<<"currentNode is "<<currentNode<<" zConstraint is "<<zConstraint<<" chosenIDX is "<<chosenIdx<<std::endl;
        Zdata<<zConstraint<<","<<z0-zRef<<","<<current_KF<<","<<current_T<<std::endl;
        addXYZEdgeSE3(T, currentNode, &graph, Lambda);
    }
    Zdata.close();

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

}

int main(int argc, char *argv[]){
    // initialize graph

    g2o::SparseOptimizer graph;
    g2o::ParameterSE3Offset* cameraOffset = new g2o::ParameterSE3Offset;
    cameraOffset->setId(0);
    graph.addParameter(cameraOffset);


    int dir;
    graph.setVerbose(true);
    YAML::Node config = YAML::LoadFile("../config.yaml");

    std::string graphConstraintsStr = config["graphConstraints"].as<std::string>();
    std::string useDEMStr = config["useDEM"].as<std::string>();
    bool useDEM=false;

    if (useDEMStr=="True" || useDEMStr=="true")
    {
        useDEM=true;
    }
    MatrixXd DEM;
    double DEM_x0;
    double DEM_y0;
    if (useDEM){

        std::string DEMStr = config["DEM"].as<std::string>();
        std::cout<<"Loading DEM..."<<std::endl;
        DEM= load_xyz<MatrixXd>(DEMStr);
        std::string DEM_x0Str = config["DEM_x0"].as<std::string>();
        std::string DEM_y0Str = config["DEM_y0"].as<std::string>();
        std::string::size_type sz, sz2;
        DEM_y0 = std::stod (DEM_x0Str,&sz);
        DEM_x0 = std::stod (DEM_y0Str,&sz2);
        std::cout<<setprecision(20)<<DEM_x0<<" "<<DEM_y0<<std::endl;
    }

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
    std::vector<int> poseIdVec;

    for (int i=0; i<graphConstraints.rows(); i++){
        poseIdVec.push_back(graphConstraints(i,1));


    }



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


        /*  for (int j=0; j<poseIdVec.size(); j++){

            if (poseId==poseIdVec[j]){
                addPosesToGraph(T, poseId, &graph);
                break;
            }
        }

*/

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
        if (T.matrix().isIdentity()){
            std::cout<<"filtering Identity"<<std::endl;
            continue;

        }

        Eigen::Matrix<double, 6, 6> Lambda;
        Lambda.setIdentity();


        Lambda(0,0)=Lambda(1,1)=Lambda(2,2)=1/graphConstraints(i,9);
        Lambda(3,3)=Lambda(4,4)=Lambda(5,5)=1/graphConstraints(i,9);
        //Lambda(0,0)=Lambda(1,1)=Lambda(2,2)=Lambda(3,3)=Lambda(4,4)=Lambda(5,5)=1/1;

        for (int j=0; j<poseIdVec.size(); j++){

            if (pose_a_id==poseIdVec[j] || pose_a_id==0 ){
                //addEdge(T, pose_a_id, pose_b_id, &graph, Lambda);
                break;
            }
        }


        addEdge(T, pose_a_id, pose_b_id, &graph, Lambda);




    }
    if (useDEM){
        Eigen::Matrix<double, 3, 3> Lambda;
        Lambda.setIdentity();

        Lambda(0,0)=Lambda(1,1)=1/99999999;
        Lambda(2,2)=1/0.1;

        //Lambda(0,0)=Lambda(1,1)=Lambda(2,2)=Lambda(3,3)=Lambda(4,4)=Lambda(5,5)=1/1;
        MatrixXd transforms0 = load_csv<MatrixXd>(transformsFile);
        int maxKF=transforms0(transforms0.rows()-1,0);
        int skip=1;

        MatrixXd transforms;

        transforms = transforms0; //load_g2o(g2oFile, transforms0);

        addDEMconstraints(DEM_x0, DEM_y0, DEM, skip, maxKF, autoMatchDir, autoNameTemplate, transforms, graph, Lambda);

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
