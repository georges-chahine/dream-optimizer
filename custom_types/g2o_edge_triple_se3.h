#ifndef EDGE_TRIPLE_SE3_H
#define EDGE_TRIPLE_SE3_H

#include "g2o/types/slam3d/vertex_se3.h"
#include "g2o/core/base_multi_edge.h"

namespace g2o {

  /**
   * \brief Connection between two camera positions and their extrinsic calibration data
   */
  class EdgeTripleSE3: public BaseMultiEdge<6, Eigen::Isometry3d>
  {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      EdgeTripleSE3();
      EdgeTripleSE3(VertexSE3* data1,VertexSE3* data2, VertexSE3* data3, VertexSE3* data4);
      typedef VertexSE3 calibration_vertex_type;

      void computeError()
      {
        //Vertice of Camera 1
        const VertexSE3* data1  = static_cast<const VertexSE3*>(_vertices[0]);
        //Vertice of Camera 2
        const VertexSE3* data2  = static_cast<const VertexSE3*>(_vertices[1]);
        //Transform between Cameras (from 1 to 2)
        const VertexSE3* data3 = static_cast<const VertexSE3*>(_vertices[2]);
        const VertexSE3* data4 = static_cast<const VertexSE3*>(_vertices[3]);


        //Abbreviations
        const Eigen::Isometry3d& data1Tf = data1->estimate();
        const Eigen::Isometry3d& data2Tf = data2->estimate();
        const Eigen::Isometry3d& data3Tf = data3->estimate();
        const Eigen::Isometry3d& data4Tf = data4->estimate();

        //Eigen::Isometry3d delta = (data1Tf * calibTf * _measurement * calibTf.inverse() * data2Tf.inverse());
        Eigen::Isometry3d delta = (data1Tf * data3Tf * _measurement * data3Tf.inverse() * data2Tf.inverse());
        _error=internal::toVectorMQT(delta);
      }

      virtual bool read(std::istream& is);
      virtual bool write(std::ostream& os) const;

  };

}

#endif
