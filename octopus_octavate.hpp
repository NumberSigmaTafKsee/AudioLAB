#pragma once


namespace Octopus
{
    Eigen::Matrix<float,1,Eigen::Dynamic,Eigen::RowMajor> Eigenize(const OctopusRowVectorXf & m)
    {
        Eigen::Matrix<float,1,Eigen::Dynamic,Eigen::RowMajor> r(m.cols());
        for(size_t i = 0; i < m.cols(); i++)        
                r(i) = m(i);
        return r;
    }
    Eigen::Matrix<double,1,Eigen::Dynamic,Eigen::RowMajor> Eigenize(const OctopusRowVectorXd & m)
    {
        Eigen::Matrix<double,1,Eigen::Dynamic,Eigen::RowMajor> r(m.cols());
        for(size_t i = 0; i < m.cols(); i++)        
                r(i) = m(i);
        return r;
    }
    Eigen::Matrix<std::complex<float>,1,Eigen::Dynamic,Eigen::RowMajor> Eigenize(const OctopusRowVectorXcf & m)
    {
        Eigen::Matrix<std::complex<float>,1,Eigen::Dynamic,Eigen::RowMajor> r(m.cols());
        for(size_t i = 0; i < m.cols(); i++)        
                r(i) = m(i);
        return r;
    }
    Eigen::Matrix<std::complex<double>,1,Eigen::Dynamic,Eigen::RowMajor> Eigenize(const OctopusRowVectorXcd & m)
    {
        Eigen::Matrix<std::complex<double>,1,Eigen::Dynamic,Eigen::RowMajor> r(m.cols());
        for(size_t i = 0; i < m.rows(); i++)        
                r(i) = m(i);
        return r;
    }

    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> Eigenize(const OctopusMatrixXf & m)
    {
        Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> r(m.rows(),m.cols());
        for(size_t i = 0; i < m.rows(); i++)
            for(size_t j = 0; j < m.cols(); j++)
                r(i,j) = m(i,j);
        return r;
    }
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> Eigenize(const OctopusMatrixXd & m)
    {
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> r(m.rows(),m.cols());
        for(size_t i = 0; i < m.rows(); i++)
            for(size_t j = 0; j < m.cols(); j++)
                r(i,j) = m(i,j);
        return r;
    }
    Eigen::Matrix<std::complex<float>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> Eigenize(const OctopusMatrixXcf & m)
    {
        Eigen::Matrix<std::complex<float>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> r(m.rows(),m.cols());
        for(size_t i = 0; i < m.rows(); i++)
            for(size_t j = 0; j < m.cols(); j++)
                r(i,j) = m(i,j);
        return r;
    }
    Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> Eigenize(const OctopusMatrixXcd & m)
    {
        Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> r(m.rows(),m.cols());
        for(size_t i = 0; i < m.rows(); i++)
            for(size_t j = 0; j < m.cols(); j++)
                r(i,j) = m(i,j);
        return r;
    }

    OctopusRowVectorXf Octavate(const Eigen::Matrix<float,1,Eigen::Dynamic,Eigen::RowMajor> & m)
    {
        OctopusRowVectorXf r(m.rows());
        for(size_t i = 0; i < m.rows(); i++)        
                r(i) = m(i);
        return r;
    }
    OctopusRowVectorXd Octavate(const Eigen::Matrix<double,1,Eigen::Dynamic,Eigen::RowMajor> & m)
    {
        OctopusRowVectorXd r(m.rows());
        for(size_t i = 0; i < m.rows(); i++)        
                r(i) = m(i);
        return r;
    }    
    OctopusRowVectorXcf Octavate(const Eigen::Matrix<std::complex<float>,1,Eigen::Dynamic,Eigen::RowMajor> & m)
    {
        OctopusRowVectorXcf r(m.rows());
        for(size_t i = 0; i < m.rows(); i++)        
                r(i) = m(i);
        return r;
    }
    OctopusRowVectorXcd Octavate(const Eigen::Matrix<std::complex<double>,1,Eigen::Dynamic,Eigen::RowMajor> & m)
    {
        OctopusRowVectorXcd r(m.rows());
        for(size_t i = 0; i < m.rows(); i++)        
                r(i) = m(i);
        return r;
    }  

    OctopusMatrixXf Octavate(const Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> & m)
    {
        OctopusMatrixXf r(m.rows(),m.cols());
        for(size_t i = 0; i < m.rows(); i++)
            for(size_t j = 0; j < m.cols(); j++)
                r(i,j) = m(i,j);
        return r;
    }
    OctopusMatrixXd Octavate(const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> & m)
    {
        OctopusMatrixXd r(m.rows(),m.cols());
        for(size_t i = 0; i < m.rows(); i++)
            for(size_t j = 0; j < m.cols(); j++)
                r(i,j) = m(i,j);
        return r;
    }
    OctopusMatrixXcf Octavate(const Eigen::Matrix<std::complex<float>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> & m)
    {
        OctopusMatrixXcf r(m.rows(),m.cols());
        for(size_t i = 0; i < m.rows(); i++)
            for(size_t j = 0; j < m.cols(); j++)
                r(i,j) = m(i,j);
        return r;
    }
    OctopusMatrixXcd Octavate(const Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> & m)
    {
        OctopusMatrixXcd r(m.rows(),m.cols());
        for(size_t i = 0; i < m.rows(); i++)
            for(size_t j = 0; j < m.cols(); j++)
                r(i,j) = m(i,j);
        return r;
    }

    
    OctopusRowVectorXf Octavate(const AudioDSP::sample_vector<float> & m)
    {
        OctopusRowVectorXf r(m.size());
        for(size_t i = 0; i < m.size(); i++)        
                r(i) = m[i];
        return r;
    }
    AudioDSP::sample_vector<float> Octavate(const OctopusRowVectorXf & m)
    {
        AudioDSP::sample_vector<float> r(m.size(1));
        for(size_t i = 0; i < m.size(1); i++)        
                r[i] = m(i);
        return r;
    }
    OctopusRowVectorXd Octavate(const AudioDSP::sample_vector<double> & m)
    {
        OctopusRowVectorXd r(m.size());
        for(size_t i = 0; i < m.size(); i++)        
                r(i) = m[i];
        return r;
    }
    AudioDSP::sample_vector<double>Octavate(const OctopusRowVectorXd & m)
    {
        AudioDSP::sample_vector<double> r(m.size(1));
        for(size_t i = 0; i < m.size(1); i++)        
                r[i] = m(i);
        return r;
    }
    OctopusRowVectorXcf Octavate(const AudioDSP::complex_vector<float> & m)
    {
        OctopusRowVectorXcf r(m.size());
        for(size_t i = 0; i < m.size(); i++)        
                r(i) = m[i];
        return r;
    }
    AudioDSP::complex_vector<float>Octavate(const OctopusRowVectorXcf & m)
    {
        AudioDSP::complex_vector<float> r(m.size(1));
        for(size_t i = 0; i < m.size(1); i++)        
                r[i] = m(i);
        return r;
    }
    OctopusRowVectorXcd Octavate(const AudioDSP::complex_vector<double> & m)
    {
        OctopusRowVectorXcd r(m.size());
        for(size_t i = 0; i < m.size(); i++)        
                r(i) = m[i];
        return r;
    }  
   

    OctopusMatrixXf Octavate(const AudioDSP::sample_matrix<float> & m)
    {
        OctopusMatrixXf r(m.rows(),m.cols());
        for(size_t i = 0; i < m.rows(); i++)
            for(size_t j = 0; j < m.cols(); j++)
                r(i,j) = m(i,j);
        return r;
    }
    AudioDSP::sample_matrix<float> Octavate(const OctopusMatrixXf & m)
    {
        AudioDSP::sample_matrix<float> r(m.rows(),m.cols());
        for(size_t i = 0; i < m.rows(); i++)
            for(size_t j = 0; j < m.cols(); j++)
                r(i,j) = m(i,j);
        return r;
    }
    OctopusMatrixXd Octavate(const AudioDSP::sample_matrix<double> & m)
    {
        OctopusMatrixXd r(m.rows(),m.cols());
        for(size_t i = 0; i < m.rows(); i++)
            for(size_t j = 0; j < m.cols(); j++)
                r(i,j) = m(i,j);
        return r;
    }
    AudioDSP::sample_matrix<double> Octavate(const OctopusMatrixXd & m)
    {
        AudioDSP::sample_matrix<double> r(m.rows(),m.cols());
        for(size_t i = 0; i < m.rows(); i++)
            for(size_t j = 0; j < m.cols(); j++)
                r(i,j) = m(i,j);
        return r;
    }
    OctopusMatrixXcf Octavate(const AudioDSP::complex_matrix<float> & m)
    {
        OctopusMatrixXcf r(m.rows(),m.cols());
        for(size_t i = 0; i < m.rows(); i++)
            for(size_t j = 0; j < m.cols(); j++)
                r(i,j) = m(i,j);
        return r;
    }
    AudioDSP::complex_matrix<float> Octavate(const OctopusMatrixXcf & m)
    {
        AudioDSP::complex_matrix<float> r(m.rows(),m.cols());
        for(size_t i = 0; i < m.rows(); i++)
            for(size_t j = 0; j < m.cols(); j++)
                r(i,j) = m(i,j);
        return r;
    }
    AudioDSP::complex_matrix<double> Octavate(const OctopusMatrixXcd & m)
    {
        AudioDSP::complex_matrix<double> r(m.rows(),m.cols());
        for(size_t i = 0; i < m.rows(); i++)
            for(size_t j = 0; j < m.cols(); j++)
                r(i,j) = m(i,j);
        return r;
    }           
}