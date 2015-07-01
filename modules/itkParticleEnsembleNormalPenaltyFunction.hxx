/*=========================================================================
Program:   ShapeWorks: Particle-based Shape Correspondence & Visualization
Module:    $RCSfile: itkParticleEnsembleNormalPenaltyFunction.txx,v $
Date:      $Date: 2011/03/24 01:17:33 $
Version:   $Revision: 1.2 $
Author:    $Author: wmartin $

Copyright (c) 2009 Scientific Computing and Imaging Institute.
See ShapeWorksLicense.txt for details.
This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.
=========================================================================*/

#ifndef __itkParticleEnsembleNormalPenaltyFunction_txx
#define __itkParticleEnsembleNormalPenaltyFunction_txx

#include "vnl/vnl_vector_fixed.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "itkParticleEnsembleNormalPenaltyFunction.h"

inline int min(int nX, int nY)
{
    return nX > nY ? nY : nX;
}

namespace itk
{

    template <unsigned int VDimension>
        itk::Matrix<double,VDimension,VDimension>
        ParticleEnsembleNormalPenaltyFunction<VDimension>
        ::RotateNormal(NormalVectorType vec) const

        { 
            RotationMatrixType mat;

            int dimension = vec.GetVectorDimension();
            double norm = vec.GetNorm();
            vec=vec/norm; //They are unit vectors already

            NormalVectorType northPole;
            northPole[0]=0; northPole[1]=0; northPole[2]=1;  

            double alpha,mult;
            alpha = 0; mult = 0;
            for (unsigned int i = 0 ; i < dimension; i ++)
            {
                //std::cout << "position = " << i << " " << vec[i] *northPole[i] << std::endl;
                mult+=vec[i]*northPole[i];
            }

            alpha=acos(mult);

            if (abs(mult) < 1e-15)
            {
                mat.SetIdentity();
            }
            else
            {
                if (abs(mult) < 1e-15)
                {
                    RotationMatrixType mat;
                    mat.SetIdentity();
                    mat*=-1;
                }
            }

            NormalVectorType aux=vec-(northPole*mult);
            norm = aux.GetNorm();
            if (norm != 0)
                aux=aux/norm;

            RotationMatrixType aux1, aux2, aux3;
            for (unsigned int i = 0 ; i < dimension; i ++)
            {
                for (unsigned int j = 0 ; j < dimension; j ++)
                {
                    aux1.GetVnlMatrix().put(i,j,northPole[i]*aux[j]);
                    //std::cout << northPole[i]*aux[j] << std::endl;
                    aux2.GetVnlMatrix().put(i,j,aux[i]*northPole[j]);
                    //std::cout << aux[i]*northPole[j] << std::endl;
                    aux3.GetVnlMatrix().put(i,j,aux1[i][j]-aux2[i][j]);
                    //std::cout << aux1[i][j]-aux2[i][j] << std::endl;
                }
            }

            RotationMatrixType identity;
            identity.SetIdentity();

            for (unsigned int i = 0 ; i < dimension; i ++)
            {
                for (unsigned int j = 0 ; j < dimension; j ++)
                {
                    aux1.GetVnlMatrix().put(i,j,northPole[i]*northPole[j]);
                    aux2.GetVnlMatrix().put(i,j,aux[i]*aux[j]);
                }
            }

            mat = identity.GetVnlMatrix() + (sin(alpha)*aux3.GetVnlMatrix()) + (cos(alpha)-1)*(aux1.GetVnlMatrix()+aux2.GetVnlMatrix());

            return mat;

        }

    template <unsigned int VDimension>
        itk::VariableSizeMatrix<double>
        ParticleEnsembleNormalPenaltyFunction<VDimension>
        ::RepVecMat(itk::VariableLengthVector<double> mat, int repRow, int repCol) const

        { 
            itk::VariableSizeMatrix<double> result(repRow,mat.Size()*repCol);
            int rows, cols;

            rows = result.GetVnlMatrix().rows();
            cols = result.GetVnlMatrix().cols(); 

            //Padding rows 
            for (int i = 0 ; i < rows ; i++)
            {
                for (int j = 0 ; j < cols ; j++)
                {
                    if ((i >= rows) && (j < cols)) //padding rows
                    {
                        result.GetVnlMatrix().put(i,j,result[i-rows][j]);
                    }
                    else
                    {
                        if ((j >= cols) && (i < rows)) //padding columns
                        {
                            result.GetVnlMatrix().put(i,j,result[i][j-cols]);
                        }
                        else
                        {
                            if ((j >= cols) && (i >= rows)) //padding columns and rows at the same time
                            {
                                result.GetVnlMatrix().put(i,j,result[i-rows][j-cols]);
                            }
                            else
                            {
                                result.GetVnlMatrix().put(i,j,mat[j]);
                            }
                        }
                    }
                }   
            }

            return result;
        }

    template <unsigned int VDimension>
        double
        ParticleEnsembleNormalPenaltyFunction<VDimension>
        ::geodesicMeanS1(vnl_vector<double> angles) const

        { 
            int size = angles.size();
            vnl_vector<double> meancandi(size,0), theta(size,0), geodvar(size,0);
            vnl_matrix<double> min_matrix(3,size,0); //there are three different ways to compute distances
            theta=angles;
            double result = -99;


            for (unsigned int i = 0 ; i < size ; i++)
            {
                double aux=(angles.mean()+2*M_PI*i/size);
                meancandi.put(i,aux-floor(aux/(2*M_PI))*(2*M_PI));
                theta.put(i,theta.get(i)-floor(theta.get(i)/(2*M_PI))*(2*M_PI));
            }

            vnl_vector<double> subtract, addition, column_vector(size,0);

            for (unsigned int i = 0 ; i < size ; i++)
            {
                subtract = theta;
                //Calculating min distance 
                subtract -= meancandi[i];

                for (unsigned int j  = 0 ; j < size ; j++)
                {
                    column_vector.put(j,(subtract[j]*subtract[j]));
                }
                min_matrix.set_row(0,column_vector); 

                addition=subtract; addition += 2*M_PI;
                for (unsigned int j  = 0 ; j < size ; j++)
                {
                    column_vector.put(j,(addition[j]*addition[j]));
                }
                min_matrix.set_row(1,column_vector);

                subtract *= -1; addition=subtract; addition += 2*M_PI;
                for (unsigned int j  = 0 ; j < size ; j++)
                {
                    column_vector.put(j,(addition[j]*addition[j]));
                }
                min_matrix.set_row(2,column_vector);

                // Print matrix (debug)
                /*for (unsigned int i = 0 ; i < min_matrix.rows() ; i ++)
                  {
                  for (unsigned int j = 0 ; j < min_matrix.cols() ; j ++)
                  {
                  std::cout << "Position row = " << i << " col = " << j << " " << min_matrix.get(i,j) << std::endl; 
                  }
                  }*/

                vnl_vector<double> min(size,99999);
                for (unsigned int j = 0 ; j < size ; j++)
                {
                    if (min.get(j) >  min_matrix.get_column(j).min_value())
                    {
                        min.put(j,min_matrix.get_column(j).min_value());
                    }
                }

                geodvar[i]+=min.sum();

            }

            double geodvar_final=geodvar.min_value();
            //unsigned int geodvar_final_index = geodvar.arg_min();
            unsigned int geodvar_final_index;
            for (unsigned int j = 0 ; j < size ; j ++)
            {
                if (geodvar[j]==geodvar_final)
                    geodvar_final_index=j;
            }

            result =  meancandi.get(geodvar_final_index)-floor(meancandi.get(geodvar_final_index)/(2*M_PI))*(2*M_PI);

            for (unsigned int j = 0 ; j < size ; j ++)
            {
                geodvar[j]=geodvar[j]/size;
            }

            return (result);
        }

    template <unsigned int VDimension>
        vnl_matrix<double>
        ParticleEnsembleNormalPenaltyFunction<VDimension>
        ::RiemannianLogMap(vnl_matrix<double> mat) const

        { 
            vnl_matrix<double> res(mat.rows()-1,mat.cols(),0);

            int dimension = mat.rows();
            int cols = mat.cols();

            vnl_vector<double> aux = mat.get_row(dimension-1);
            vnl_vector<double> aux_acos = mat.get_row(dimension-1);
            vnl_vector<double> aux_sqrt = mat.get_row(dimension-1);
            //vnl_vector<double> scale = mat.get_row(dimension-1);
            itk::VariableLengthVector<double> scale(cols);

            for (unsigned int i = 0 ; i < aux.size() ; i ++)
            {
                aux_acos[i]=acos(aux[i]);
                aux_sqrt[i]=sqrt(1-(aux[i]*aux[i]));
                scale.SetElement(i,(aux_acos[i]/aux_sqrt[i]));
                if (std::isnan(scale[i]))
                    scale.SetElement(i,1);
            }

            vnl_matrix<double> aux1(dimension-1,cols,0);
            aux1 = RepVecMat(scale,dimension-1,1).GetVnlMatrix();
            vnl_matrix<double> aux2(dimension-1,cols,0);
            aux2 = mat.extract(dimension-1,cols,0,0);

            for (unsigned int i = 0 ; i < aux1.rows() ; i ++)
            {
                for (unsigned int j = 0 ; j < aux1.cols() ; j ++)
                {
                    res[i][j]=aux1[i][j]*aux2[i][j];
                    //std::cout << res[i][j] << " " << i << "-" << j << std::endl;
                }
            }  

            return res.as_ref();
        }

    template <unsigned int VDimension>
        vnl_matrix<double>
        ParticleEnsembleNormalPenaltyFunction<VDimension>
        ::ExponentialLogMap(vnl_matrix<double> mat) const

        { 
            int rows = mat.rows();
            int cols = mat.cols();
            int dimension = rows;

            vnl_matrix<double> mat_square(rows,cols,0);
            double sum=0;
            for (unsigned int i = 0 ; i < rows ; i ++)
            {
                for (unsigned int j = 0 ; j < cols ; j ++)
                {
                    mat_square[i][j]=mat[i][j]*mat[i][j];
                    sum+=mat_square[i][j];
                }
            }

            itk::VariableLengthVector<double> nv(cols);
            nv.SetElement(0,sqrt(sum));
            itk::VariableLengthVector<double> nv_sin(cols);
            itk::VariableLengthVector<double> nv_cos(cols);

            for (unsigned int i = 0 ; i < nv_sin.Size() ; i++)
            {
                nv_sin[i]=sin(nv[i])/nv[i];
                nv_cos[i]=cos(nv[i]);
            }

            vnl_matrix<double> aux(dimension-1,cols,0);
            aux = RepVecMat(nv_sin,dimension,1).GetVnlMatrix();

            vnl_matrix<double> res(rows+1,cols,0);

            for (unsigned int i = 0 ; i < res.rows() ; i ++)
            {
                for (unsigned int j = 0 ; j < res.cols() ; j ++)
                {
                    if (i!=rows)
                    {
                        res[i][j]=mat[i][j]*aux[i][j];
                        if (isnan(res[i][j]))
                            res[i][j]=0;
                    }
                    else
                    {
                        res[i][j]=nv_cos[j];
                        if (isnan(res[i][j]))
                            res[i][j]=0;
                    }
                }
            }

            if (nv[0] < 1e-16) 
            {
                res[0][0]=0; res[1][0]=0; res[2][0]=1;
            }
            //THIS LINE HASNT BEEN CODED Exppx(:,nv < 1e-16) = repmat([zeros(d,1);1],1,sum(nv < 1e-16));
            return res.as_ref();
        }

    template <unsigned int VDimension>
        vnl_matrix<double>
        ParticleEnsembleNormalPenaltyFunction<VDimension>
        ::ComputePrincipalScores (vnl_matrix<double> data) const
        { 
            int dimension = data.rows();
            int nNormals = data.columns();
            itk::Vector<double,VDimension> c0;
            for (unsigned int i = 0 ; i < c0.Size(); i ++)
                c0.SetElement(i, data[i][0]);

            double diff, error, alpha, cnt;
            diff = 1;
            error = 1e-10;
            alpha = 2; //Scalar for gradient works well for alpha in [1,2]
            cnt = 1;

            int nan=1;

            while ((diff>error)&&(cnt<100))
            {
                itk::Matrix<double,VDimension,VDimension> rot = RotateNormal(c0);
                vnl_matrix<double> x = rot.GetVnlMatrix() * data;

                //std::cout << "First log map" << std::endl;
                vnl_matrix<double> logmap = RiemannianLogMap(x);
                itk::VariableLengthVector<double> normlogmap(logmap.cols());  
                itk::VariableLengthVector<double> id(logmap.cols());
                itk::VariableLengthVector<double> norm_division(logmap.cols());
                int size_div = 0;

                for (unsigned int j = 0 ; j < logmap.cols(); j ++)
                {
                    normlogmap.SetElement(j,0);
                    for (unsigned int i = 0 ; i <logmap.rows(); i ++)
                    {
                        double pow_value = logmap.get(i,j)*logmap.get(i,j);
                        double acc_value = normlogmap.GetElement(j)+pow_value;
                        normlogmap.SetElement(j,acc_value);
                    }
                    double sqrt_value = sqrt(normlogmap.GetElement(j));

                    if (sqrt_value > 1e-10)
                    {
                        id.SetElement(j,1);
                        size_div++;
                    }	
                    else
                    {
                        id.SetElement(j,0);
                        norm_division.SetElement(j,0);
                    }
                    normlogmap.SetElement(j,sqrt_value);
                    if (id.GetElement(j))
                    {
                        double division = (normlogmap.GetElement(j) - M_PI/2) / normlogmap.GetElement(j);
                        //std::cout << "division " << division << std::endl;
                        norm_division.SetElement(j,division);
                    }
                }

                vnl_matrix<double> aux1(dimension-1,logmap.cols(),0);
                aux1 = RepVecMat(norm_division,dimension-1,1).GetVnlMatrix();
                //std::cout << "norm logmap " << normlogmap << std::endl;
                vnl_matrix<double> logcnew (logmap.rows(),1,0);

                vnl_matrix<double> final_map = logmap; 
                for (unsigned int j = 0 ; j < logmap.cols(); j ++)
                {
                    if (id.GetElement(j))
                    {
                        for (unsigned int i = 0 ; i < logmap.rows(); i ++)
                        {
                            double value = logmap.get(i,j) * aux1.get(i,j);
                            double sum = value + logcnew.get(i,0);
                            logcnew[i][0]=sum;
                        }
                    }
                }

                double aux = logcnew.get(0,0) ; // removed /size_div
                logcnew[0][0]=aux*2;
                aux = logcnew.get(0,1) ; //removed / size_div
                logcnew[0][1]=aux*2;	

                vnl_matrix<double> expcnew = ExponentialLogMap(logcnew);
                //	std::cout << "expcnew " << expcnew << std::endl;
                vnl_matrix<double> cnew = vnl_inverse(rot.GetVnlMatrix())*expcnew;

                //std::cout << "rot " << rot << std::endl;
                //std::cout << "expnew " << expcnew << std::endl;	

                vnl_vector<double> norm_vector (c0.Size());

                if ((c0[0]==cnew[0][0])&&(c0[1]==cnew[1][0])&&(c0[2]==cnew[2][0]))
                {
                    cnt=101; // To reset the search of this great circle and start searching for a new one	
                }
                //	std::cout << c0 << " " << cnew << std::endl;
                for (unsigned int i = 0 ; i < norm_vector.size(); i ++)
                {
                    norm_vector[i]=cnew[i][0]-c0[i];
                    c0[i]=cnew[i][0];
                }

                diff=norm_vector.two_norm();
                //	std::cout << "diff " << diff << std::endl;
                cnt++;			
            }

            itk::Vector<double,VDimension> c = c0;
            itk::Matrix<double,VDimension,VDimension> rot = RotateNormal(c);
            //  std::cout << "rot " << rot << std::endl;

            vnl_matrix<double> x = rot.GetVnlMatrix() * data;
            // std::cout << "Second log map" << std::endl;
            vnl_matrix<double> logmap = RiemannianLogMap(x);

            vnl_vector<double> angles(logmap.cols(), 0);
            for (unsigned int i = 0 ; i < angles.size() ; i ++)
            {
                double result = atan2(logmap[1][i],logmap[0][i]);
                angles.put(i,result);
            }

            double meantheta = geodesicMeanS1 (angles);
            vnl_matrix<double> Z_scores(2,nNormals,0);
            //vnl_vector<double> Z1_scores (angles.size());
            vnl_vector<double> Z2_scores = c.Get_vnl_vector()*data;
            // double perturbation = M_PI/2000;
            for (unsigned int i = 0 ; i < angles.size() ; i++)
            {
                double residuals1=angles[i] - meantheta + M_PI;
                Z_scores.put(0,i,(residuals1-floor(residuals1/(2*M_PI))*(2*M_PI))-M_PI);
                Z_scores.put(1,i,acos(Z2_scores[i])-M_PI/2); 
            }

            return (Z_scores);
        }

    template <unsigned int VDimension>
        vnl_matrix<double>
        ParticleEnsembleNormalPenaltyFunction<VDimension>
        ::ComputeGeodesicMean (vnl_matrix<double> data) const
        { 
            int dimension = data.rows();
            int nNormals = data.columns();
            itk::Vector<double,VDimension> c0;
            for (unsigned int i = 0 ; i < c0.Size(); i ++)
                c0.SetElement(i, data[i][0]);

            double diff, error, alpha, cnt;
            diff = 1;
            error = 1e-10;
            alpha = 2; //Scalar for gradient works well for alpha in [1,2]
            cnt = 1;

            int nan=1;

            while ((diff>error)&&(cnt<100))
            {
                itk::Matrix<double,VDimension,VDimension> rot = RotateNormal(c0);
                vnl_matrix<double> x = rot.GetVnlMatrix() * data;

                vnl_matrix<double> logmap = RiemannianLogMap(x);
                itk::VariableLengthVector<double> normlogmap(logmap.cols());  
                itk::VariableLengthVector<double> id(logmap.cols());
                itk::VariableLengthVector<double> norm_division(logmap.cols());
                int size_div = 0;

                for (unsigned int j = 0 ; j < logmap.cols(); j ++)
                {
                    normlogmap.SetElement(j,0);
                    for (unsigned int i = 0 ; i <logmap.rows(); i ++)
                    {
                        double pow_value = logmap.get(i,j)*logmap.get(i,j);
                        double acc_value = normlogmap.GetElement(j)+pow_value;
                        normlogmap.SetElement(j,acc_value);
                    }
                    double sqrt_value = sqrt(normlogmap.GetElement(j));

                    if (sqrt_value > 1e-10)
                    {
                        id.SetElement(j,1);
                        size_div++;
                    }	
                    else
                    {
                        id.SetElement(j,0);
                        norm_division.SetElement(j,0);
                    }
                    normlogmap.SetElement(j,sqrt_value);
                    if (id.GetElement(j))
                    {
                        double division = (normlogmap.GetElement(j) - M_PI/2) / normlogmap.GetElement(j);
                        norm_division.SetElement(j,division);
                        //std::cout << norm_division.GetElement(j) << std::endl;
                    }
                }

                vnl_matrix<double> aux1(dimension-1,logmap.cols(),0);
                aux1 = RepVecMat(norm_division,dimension-1,1).GetVnlMatrix();
                vnl_matrix<double> logcnew (logmap.rows(),1,0);

                vnl_matrix<double> final_map = logmap; 
                for (unsigned int j = 0 ; j < logmap.cols(); j ++)
                {
                    if (id.GetElement(j))
                    {
                        for (unsigned int i = 0 ; i < logmap.rows(); i ++)
                        {
                            double value = logmap.get(i,j) * aux1.get(i,j);
                            double sum = value + logcnew.get(i,0);
                            logcnew[i][0]=sum;
                        }
                    }
                }

                double aux = logcnew.get(0,0) / size_div;
                logcnew[0][0]=aux*2;
                aux = logcnew.get(0,1) / size_div;
                logcnew[0][1]=aux*2;

                vnl_matrix<double> expcnew = ExponentialLogMap(logcnew);
                vnl_matrix<double> cnew = vnl_inverse(rot.GetVnlMatrix())*expcnew;


                vnl_vector<double>  norm_vector (c0.Size());

                if ((c0[0]==cnew[0][0])&&(c0[1]==cnew[1][0])&&(c0[2]==cnew[2][0]))
                {
                    cnt=101; // To reset the search of this great circle and start searching for a new one	
                }

                for (unsigned int i = 0 ; i < norm_vector.size(); i ++)
                {
                    norm_vector[i]=cnew[i][0]-c0[i];
                    c0[i]=cnew[i][0];
                }
                diff=c0.Get_vnl_vector().two_norm();
                cnt++;			
            }

            itk::Vector<double,VDimension> c = c0;
            itk::Matrix<double,VDimension,VDimension> rot = RotateNormal(c);
            vnl_matrix<double> x = rot.GetVnlMatrix() * data;
            vnl_matrix<double> logmap = RiemannianLogMap(x);

            vnl_vector<double> angles(logmap.cols(), 0);
            for (unsigned int i = 0 ; i < angles.size() ; i ++)
            {
                double result = atan2(logmap[1][i],logmap[0][i]);
                angles.put(i,result);
            }

            // Compute the geodesicMean that best fits my normal data
            double meantheta = geodesicMeanS1 (angles);

            //Calculate the coordinates of the point in the tangent plane
            vnl_matrix<double> mean_tangent_plane (dimension-1,1,0);
            mean_tangent_plane[0][0]=cos(meantheta);
            mean_tangent_plane[1][0]=sin(meantheta);
            //Calculate the mean vector
            vnl_matrix<double> mean_normal = ExponentialLogMap(mean_tangent_plane);

            return (mean_normal);
        }

    template <unsigned int VDimension>
        vnl_matrix<double>
        ParticleEnsembleNormalPenaltyFunction<VDimension>
        ::ComputePrincipalScoresGradient (unsigned int particleIdx, unsigned int d, const ParticleSystemType *system, double energy) const
        {
            const double epsilon = 1.0e-8;
            const double N = (double)(system->GetNumberOfDomains() / m_DomainsPerShape);

            // Get the position for which we are computing the gradient
            PointType pos = system->GetTransformedPosition(particleIdx, d);

            // get domain information with gradients
            const ParticleImageDomainWithGradients<float, VDimension> * domain
                = static_cast<const ParticleImageDomainWithGradients<float, VDimension> *>(system->GetDomain(d));

            vnl_vector <int> final_level_set (27,1); // We are going to evaluate 27 positions around the voxel we are computing gradient for

            itk::Vector<double> imageSpacing;

            imageSpacing = domain->GetImageDomainSpacing();
            double minImageSpacing = min (min(imageSpacing[0],imageSpacing[1]),imageSpacing[2]);

            std::vector<double> energies;
            std::vector<double> deltaEnergies;
            std::vector<PointType> neighbors;
            int cont = 0;
            int min_energy_neighbor=0;
            double min_energy=99999999999999999;

            for (int x_index = -1; x_index < 2 ; x_index ++)
            {
                for (int y_index = 1; y_index > -2 ; y_index --)
                {
                    for (int z_index = -1; z_index < 2 ; z_index ++)
                    {
                        if ((x_index==0) && (y_index==0) && (z_index==0))
                        {
                            // Do nothing, just print value to really see we have our particle in level set.
                            //std::cout << "CENTRAL VOXEL = " << domain->Sample(pos) << std::endl;		
                        }
                        else
                        {

                            PointType newpos;
                            newpos[0]=pos[0]+x_index*imageSpacing[0]; newpos[1]=pos[1]+y_index*imageSpacing[1]; newpos[2]=pos[2]+z_index*imageSpacing[2];
                            if (( domain->Sample(newpos) < 0.5*minImageSpacing ) && ( domain->Sample(newpos) > - (0.5*minImageSpacing) ))
                            {
                                //POTENTIAL LEVEL SET 

                                //Here compute ppal scores and store the energy
                                //***************** COMPUTE PPAL SCORES
                                vnl_matrix<double> normalDataParticleIdx(VDimension, system->GetNumberOfDomains(), 0);
                                for (unsigned int i = d % m_DomainsPerShape; i < system->GetNumberOfDomains(); i += m_DomainsPerShape)
                                {
                                    PointType neighpos;
                                    if (i != d)
                                        neighpos = system->GetTransformedPosition(particleIdx, i);
                                    else
                                        neighpos = newpos;

                                    domain = static_cast<const ParticleImageDomainWithGradients<float, VDimension> *>(system->GetDomain(i));
                                    typename ParticleImageDomainWithGradients<float,VDimension>::VnlVectorType 
                                        neighnormal = domain->SampleNormalVnl(neighpos);

                                    // Store that normal for entropy computation
                                    for (unsigned int dim = 0 ; dim < VDimension ; dim ++)
                                    {
                                        normalDataParticleIdx.put(dim,i,neighnormal[dim]);
                                    }
                                }
                                //std::cout << "BEFORE PPAL SCORES" << std::endl;
                                vnl_matrix<double> scores = ComputePrincipalScores (normalDataParticleIdx);
                                //***************** COMPUTE ENERGY
                                // Perform eigen-analysis in the ppal scores
                                vnl_matrix_type points_minus_mean(scores.rows(),scores.columns());
                                vnl_vector_type means(PDimensions);       

                                // Compute the mean shape vector. Y
                                for (unsigned int j = 0; j < PDimensions; j++)
                                {
                                    double total = 0.0;
                                    for (unsigned int i = 0; i < scores.columns(); i++)
                                    {
                                        points_minus_mean[j][i]=scores[j][i];
                                        total += points_minus_mean[j][i];
                                    }
                                    means[j] = (double)(total/scores.columns());
                                }

                                for (unsigned int j = 0; j < PDimensions; j++)
                                {
                                    for (unsigned int i = 0; i < scores.columns(); i++)
                                    {
                                        points_minus_mean(j,i) -= means(j);
                                    }
                                }

                                vnl_matrix_type A =  points_minus_mean * points_minus_mean.transpose() ;

                                // Find inverse of covariance matrix
                                vnl_symmetric_eigensystem<double> symEigen(A);

                                double new_energy = 0.0;

                                for (unsigned int i = 0; i < PDimensions; i++)
                                {
                                    //std::cout << "Eigenvalue " << i << " " << symEigen.D(i,i) << std::endl;
                                    if ((isnan(- log(symEigen.D(i,i)))) || (isinf(- log(symEigen.D(i,i)))))
                                    {
                                        new_energy +=0;
                                        //std::cout << "WARNING!!! " << symEigen.D(i,i) << std::endl;
                                    }
                                    else
                                    {
                                        new_energy += - log(symEigen.D(i,i));
                                        //new_energy += log(symEigen.D(i,i));

                                        //std::cout << new_energy << std::endl;
                                    }
                                }

                                //new_energy = new_energy;
                                energies.push_back(new_energy);
                                deltaEnergies.push_back(abs(new_energy-energy));
                                //std::cout << "Energies " << abs(new_energy-energy) << std::endl;
                                neighbors.push_back(newpos);
                                //std::cout << "Neighbors " << neighbors.size() << std::endl;
                                //std::cout <<  "New pos " << newpos << std::endl;
                                //					if ((new_energy < min_energy) && (new_energy > 0)) // For some reason when energy goes to negative the program gets stuck
                                if (new_energy < min_energy) 
                                {
                                    min_energy=new_energy;
                                    min_energy_neighbor=cont;
                                }	
                                cont++;
                            }
                            else
                            {
                                // If it is not in the level set do nothing
                                //std::cout << "x(" << x_index << ") y(" << y_index << ") z(" << z_index << ") neighbour content " << domain->Sample(newpos) << std::endl;
                            }
                        }
                    }
                }
            }

            //******** COMPUTE GRADIENT 
            vnl_matrix<double> gradient(1,3,0);

            vnl_vector<double> h_vector(VDimension,0.0);
            PointType neigh_pos  = neighbors[min_energy_neighbor];
            for ( int dim = 0 ; dim < VDimension ; dim ++ )
            {
                h_vector[dim] = neigh_pos[dim]-pos[dim];
            }

            for ( int dim = 0 ; dim < VDimension ; dim ++ )
            {
                double value = deltaEnergies[min_energy_neighbor] / h_vector.magnitude();
                //std::cout << min_energy_neighbor << " " << deltaEnergies[min_energy_neighbor] <<std::endl;
                gradient.put(0,dim,(h_vector[dim]*value)); 
            }	

            //std::cout << gradient.get(0,0) << " " << gradient.get(0,1) << " " << gradient.get(0,2) << std::endl;
            return (gradient);
        }

    template <unsigned int VDimension>
        typename ParticleEnsembleNormalPenaltyFunction<VDimension>::VectorType
        ParticleEnsembleNormalPenaltyFunction<VDimension>
        ::Evaluate(unsigned int idx, unsigned int d, const ParticleSystemType * system,
                double &maxmove, double &energy) const
        {
            const double epsilon = 1.0e-8;
            const double N = (double)(system->GetNumberOfDomains() / m_DomainsPerShape);

            // Get the position for which we are computing the gradient
            PointType pos = system->GetTransformedPosition(idx, d);

            // get domain information with gradients
            const ParticleImageDomainWithGradients<float, VDimension> * domain
                = static_cast<const ParticleImageDomainWithGradients<float, VDimension> *>(system->GetDomain(d));

            // get normal for current position
            typename ParticleImageDomainWithGradients<float,VDimension>::VnlVectorType posnormal = domain->SampleNormalVnl(pos);

            // find mean normal for sister particles across ensemble 
            vnl_vector<double> mean_normal(VDimension,0.0f);
            vnl_matrix<double> normalDataParticleIdx(VDimension, system->GetNumberOfDomains(), 0);
            for (unsigned int i = d % m_DomainsPerShape; i < system->GetNumberOfDomains(); i += m_DomainsPerShape)
            {
                PointType neighpos = system->GetTransformedPosition(idx, i);

                domain = static_cast<const ParticleImageDomainWithGradients<float, VDimension> *>(system->GetDomain(i));
                typename ParticleImageDomainWithGradients<float,VDimension>::VnlVectorType 
                    neighnormal = domain->SampleNormalVnl(neighpos);

                // Store that normal for entropy computation
                for (unsigned int dim = 0 ; dim < VDimension ; dim ++)
                {
                    normalDataParticleIdx.put(dim,i,neighnormal[dim]);
                }
            }

            //Here is where I have all the normal data for a certain particle, I can compute geodesic mean
            //vnl_matrix<double> geodesic_mean (VDimension,1,0);
            //geodesic_mean= ComputeGeodesicMean (normalDataParticleIdx);

            vnl_matrix<double> scores = ComputePrincipalScores (normalDataParticleIdx);

            //std::cout << "scores: " << scores << std::endl;

            //*********************************************
            // Perform eigen-analysis in the ppal scores
            vnl_matrix_type points_minus_mean(scores.rows(),scores.columns());
            vnl_vector_type means(PDimensions);       

            // Compute the mean shape vector. Y
            for (unsigned int j = 0; j < PDimensions; j++)
            {
                double total = 0.0;
                for (unsigned int i = 0; i < scores.columns(); i++)
                {
                    points_minus_mean[j][i]=scores[j][i];
                    total += points_minus_mean[j][i];
                }
                means[j] = (double)(total/scores.columns());
            }

            for (unsigned int j = 0; j < PDimensions; j++)
            {
                for (unsigned int i = 0; i < scores.columns(); i++)
                {
                    points_minus_mean(j,i) -= means(j);
                }
            }

            // Compute the covariance in the dual space (transposed shape matrix)
            vnl_matrix_type A = points_minus_mean * points_minus_mean.transpose();
            //  std::cout << "HERE A " << A << std::endl;
            // Find inverse of covariance matrix
            vnl_symmetric_eigensystem<double> symEigen(A);

            double m_MinimumEigenValue = symEigen.D(0, 0);
            double m_CurrentEnergy = 0.0;

            for (unsigned int i = 0; i < PDimensions; i++)
            {
                if ((isnan(- log(symEigen.D(i,i)))) || (isinf(- log(symEigen.D(i,i)))))
                {
                    m_CurrentEnergy +=0;
                }
                else
                    m_CurrentEnergy +=  - log(symEigen.D(i,i));
                //m_CurrentEnergy +=  log(symEigen.D(i,i));
            }

            vnl_matrix<double> gradient_ppalScores = ComputePrincipalScoresGradient(idx,d,system, m_CurrentEnergy);
            //*********************************************
            vnl_vector<double> gradE_norm(VDimension,0.0); 
            gradE_norm[0]= gradient_ppalScores.get(0,0); gradE_norm[1]= gradient_ppalScores.get(0,1); gradE_norm[2]= gradient_ppalScores.get(0,2);  
            energy = m_CurrentEnergy ;
            //energy = (-1) * m_CurrentEnergy ;

            maxmove = domain->GetImage()->GetSpacing()[0];
            //maxmove =  fabs(energy) * 0.5 + 0.00001; // to prevent the system getting stuck in the optimizer

            std::cout << "GRADIENT: " << gradient_ppalScores.get(0,0) << " " << gradient_ppalScores.get(0,1) << " " << gradient_ppalScores.get(0,2) << std::endl;
            std::cout << "ENERGY: " << energy << std::endl;

            /* if (energy < 0)
               {
               std::cout << "something happened here" << std::endl;
            //energy = (-1) * energy;
            }*/  

            std::cout << "--" << std::endl;

            //  Transform the gradient according to the transform of the given domain and return.
            return system->TransformVector(gradE_norm, system->GetInversePrefixTransform(d)
                    * system->GetInverseTransform(d));
        }




} // end namespace

#endif

