#pragma once


// cant really wrap this
template<typename T, int X=Eigen::Dynamic, int Y=Eigen::Dynamic, int Z=Eigen::RowMajor>
struct VectorwiseRowOp : public Eigen::VectorwiseOp<Eigen::Matrix<T,X,Y,Z>, Eigen::Horizontal>
{
    using ExpressionType = Eigen::VectorwiseOp<Eigen::Matrix<T,X,Y,Z>, Eigen::Horizontal>;

    
    VectorwiseRowOp(const ExpressionType& matrix) : 
        Eigen::VectorwiseOp<Eigen::Matrix<T,X,Y,Z>, Eigen::Horizontal>(matrix)
    {

    }
    


    // these do not go through because you can not instantiate it with swig
    using Eigen::VectorwiseOp<Eigen::Matrix<T,X,Y,Z>, Eigen::Horizontal>::minCoeff;
    using Eigen::VectorwiseOp<Eigen::Matrix<T,X,Y,Z>, Eigen::Horizontal>::maxCoeff;
    using Eigen::VectorwiseOp<Eigen::Matrix<T,X,Y,Z>, Eigen::Horizontal>::squaredNorm;    
    using Eigen::VectorwiseOp<Eigen::Matrix<T,X,Y,Z>, Eigen::Horizontal>::norm;
    using Eigen::VectorwiseOp<Eigen::Matrix<T,X,Y,Z>, Eigen::Horizontal>::blueNorm;
    using Eigen::VectorwiseOp<Eigen::Matrix<T,X,Y,Z>, Eigen::Horizontal>::stableNorm;
    using Eigen::VectorwiseOp<Eigen::Matrix<T,X,Y,Z>, Eigen::Horizontal>::hypotNorm;
    using Eigen::VectorwiseOp<Eigen::Matrix<T,X,Y,Z>, Eigen::Horizontal>::sum;
    using Eigen::VectorwiseOp<Eigen::Matrix<T,X,Y,Z>, Eigen::Horizontal>::mean;
    using Eigen::VectorwiseOp<Eigen::Matrix<T,X,Y,Z>, Eigen::Horizontal>::all;
    using Eigen::VectorwiseOp<Eigen::Matrix<T,X,Y,Z>, Eigen::Horizontal>::any;
    using Eigen::VectorwiseOp<Eigen::Matrix<T,X,Y,Z>, Eigen::Horizontal>::count;
    using Eigen::VectorwiseOp<Eigen::Matrix<T,X,Y,Z>, Eigen::Horizontal>::prod;
    using Eigen::VectorwiseOp<Eigen::Matrix<T,X,Y,Z>, Eigen::Horizontal>::reverse;                        
                
};