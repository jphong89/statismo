/*=========================================================================
Program:   ShapeWorks: Particle-based Shape Correspondence & Visualization
Module:    $RCSfile: itkParticleEnsembleNormalPenaltyFunction.h,v $
Date:      $Date: 2011/03/24 01:17:33 $
Version:   $Revision: 1.2 $
Author:    $Author: wmartin $

Copyright (c) 2009 Scientific Computing and Imaging Institute.
See ShapeWorksLicense.txt for details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#ifndef __itkParticleEnsembleNormalPenaltyFunction_h
#define __itkParticleEnsembleNormalPenaltyFunction_h

#include "itkParticleVectorFunction.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include "itkParticleImageDomainWithGradients.h"
#include "itkDerivativeOperator.h"
#include "itkNeighborhoodIterator.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIterator.h"
#include "itkVariableLengthVector.h"
#include "itkVariableSizeMatrix.h"

#define PDimensions 2


namespace itk
{

    /**
     * \class ParticleEnsembleNormalPenaltyFunction
     *
     */
    template <unsigned int VDimension>
        class ParticleEnsembleNormalPenaltyFunction : public ParticleVectorFunction<VDimension>
    {
        public:
            /** Standard class typedefs. */
            typedef ParticleEnsembleNormalPenaltyFunction Self;
            typedef SmartPointer<Self>  Pointer;
            typedef SmartPointer<const Self>  ConstPointer;
            typedef ParticleVectorFunction<VDimension> Superclass;
            itkTypeMacro( ParticleEnsembleNormalPenaltyFunction, ParticleVectorFunction);

            /** Type of particle system. */
            typedef typename Superclass::ParticleSystemType ParticleSystemType;

            /** Vector & Point types. */
            typedef typename Superclass::VectorType VectorType;
            typedef typename ParticleSystemType::PointType PointType;
            typedef vnl_vector<double> vnl_vector_type;
            typedef vnl_matrix<double> vnl_matrix_type;

            /** Vector & Point types for computing normal entropy. */
            typedef std::vector<vnl_matrix_type> NormalDataType;
            typedef itk::Vector<double,VDimension> NormalVectorType;
            typedef itk::Matrix<double,3,3> RotationMatrixType;

            // Normal Penalty related classes
            typedef typename itk::Image<float,VDimension> NormalComponentImageType;
            typedef typename itk::NeighborhoodIterator<NormalComponentImageType> NeighborhoodIteratorType;
            typedef typename itk::ImageRegionIterator<NormalComponentImageType> ImageIteratorType;
            typedef typename itk::DerivativeOperator<float, VDimension> DerivativeOperatorType;

            /** Method for creation through the object factory. */
            itkNewMacro(Self);

            /** Dimensionality of the domain of the particle system. */
            itkStaticConstMacro(Dimension, unsigned int, VDimension);

            /** The first argument is a pointer to the particle system.  The second
              argument is the index of the domain within that particle system.  The
              third argument is the index of the particle location within the given
              domain. */
            virtual VectorType Evaluate(unsigned int a, unsigned int b,
                    const ParticleSystemType *c, double &d, double &e) const;

            virtual VectorType Evaluate(unsigned int a,  unsigned int b, 
                    const ParticleSystemType *c, double &d) const
            {
                double e;
                return this->Evaluate(a,b,c,d,e);
            }

            virtual double Energy(unsigned int a, unsigned int b,
                    const ParticleSystemType *c) const
            {
                double e, d;
                this->Evaluate(a,b,c,d,e);
                return e;
            }

            /** Set/Get the number of domains per shape.  Separate means are used for
              separate shape domains. */
            void SetDomainsPerShape(int i)
            { m_DomainsPerShape = i; }
            int GetDomainsPerShape() const
            { return m_DomainsPerShape; }

        protected:
            ParticleEnsembleNormalPenaltyFunction()
            {
                m_DomainsPerShape = 1;

                imgDuplicator = itk::ImageDuplicator<NormalComponentImageType>::New();

                normComp = NormalComponentImageType::New();
                typename NormalComponentImageType::SizeType tmpSize; tmpSize.Fill(VDimension);
                typename NormalComponentImageType::IndexType tmpIndex; tmpIndex.Fill(0);
                typename NormalComponentImageType::RegionType tmpRegion; tmpRegion.SetSize(tmpSize); tmpRegion.SetIndex(tmpIndex);
                normComp->SetRegions(tmpRegion); normComp->Allocate(); normComp->FillBuffer(0.0f);

                radius.Fill(1);
            }

            virtual ~ParticleEnsembleNormalPenaltyFunction() {}
            //void operator=(const ParticleEnsembleNormalPenaltyFunction &);
            ParticleEnsembleNormalPenaltyFunction(const ParticleEnsembleNormalPenaltyFunction &);

            int m_DomainsPerShape;  


            /** PNS normal entropy computing functions and attributes*/
            itk::Matrix<double,VDimension,VDimension> RotateNormal (NormalVectorType vec) const;
            vnl_matrix<double> RiemannianLogMap (vnl_matrix<double> mat) const;
            vnl_matrix<double> ExponentialLogMap (vnl_matrix<double> mat) const;
            itk::VariableSizeMatrix<double> RepVecMat (itk::VariableLengthVector<double> mat, int repRow, int repCol) const;
            double geodesicMeanS1 (vnl_vector<double> angles) const;
            vnl_matrix<double> ComputePrincipalScores (vnl_matrix<double> mat) const;
            vnl_matrix<double> ComputeGeodesicMean (vnl_matrix<double> mat) const;
            vnl_matrix<double> ComputePrincipalScoresGradient (unsigned int particleIdx, unsigned int d, const ParticleSystemType *system, double energy) const;

            //vnl_matrix<double> m_PrincipalScoresData;
            double m_MinimumVariance;
            //double m_MinimumEigenValue;
            //double m_CurrentEnergy;
            /*std::vector <vnl_matrix<double> > m_PrincipalScoresGradient;
              virtual void ComputeUpdatesNormals();
              virtual void PrintInformationFile(const ParticleSystemType *system);*/
            /** PNS normal entropy computing functions and attributes*/

            typename itk::ImageDuplicator<NormalComponentImageType>::Pointer imgDuplicator;

            typename NormalComponentImageType::Pointer normComp;

            typename NeighborhoodIteratorType::RadiusType radius; 

            typename itk::NeighborhoodInnerProduct<NormalComponentImageType> inprod;
    };

} //end namespace

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkParticleEnsembleNormalPenaltyFunction+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkParticleEnsembleNormalPenaltyFunction.txx"
#endif

#endif
