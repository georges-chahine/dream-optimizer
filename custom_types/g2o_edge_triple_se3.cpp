#include "g2o_edge_triple_se3.h"
#include "g2o/core/factory.h"
#include "g2o/stuff/macros.h"
namespace g2o {
using namespace Eigen;

  EdgeTripleSE3::EdgeTripleSE3() :
    BaseMultiEdge<6, Eigen::Isometry3d>()
  {
    resize(4);
  }

  EdgeTripleSE3::EdgeTripleSE3(VertexSE3* data1,VertexSE3* data2, VertexSE3* data3, VertexSE3* data4) :
    BaseMultiEdge<6, Eigen::Isometry3d>()
  {
    resize(4);
    _vertices[0] = data1;
    _vertices[1] = data2;
    _vertices[2] = data3;
    _vertices[3] = data4;
    _measurement = Eigen::Isometry3d::Identity();
  }


  bool EdgeTripleSE3::read(std::istream& is)
  {
    Vector7d meas;
    for (int i=0; i<7; i++) 
      is >> meas[i];
    // normalize the quaternion to recover numerical precision lost by storing as human readable text
    Vector4d::MapType(meas.data()+3).normalize();
    setMeasurement(internal::fromVectorQT(meas));

    if (is.bad()) {
      return false;
    }
    for ( int i=0; i<information().rows() && is.good(); i++)
      for (int j=i; j<information().cols() && is.good(); j++){
        is >> information()(i,j);
        if (i!=j)
          information()(j,i)=information()(i,j);
      }
    if (is.bad()) {
      //  we overwrite the information matrix with the Identity
      information().setIdentity();
    } 
    return true;
  }

  bool EdgeTripleSE3::write(std::ostream& os) const
  {
    //os << _offsetParam->id() <<  " ";
    Vector7d meas=internal::toVectorQT(_measurement);
    for (int i=0; i<6; i++) os  << meas[i] << " ";
    for (int i=0; i<information().rows(); i++)
      for (int j=i; j<information().cols(); j++) {
        os <<  information()(i,j) << " ";
      }
    return os.good();
  }

  G2O_REGISTER_TYPE(EDGE_TRIPLE_SE3, EdgeTripleSE3);
} // end namespace
