/*

Copyright (c) 2005-2017, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/*

 Copyright (c) 2005-2015, University of Oxford.
 All rights reserved.

 University of Oxford means the Chancellor, Masters and Scholars of the
 University of Oxford, having an administrative office at Wellington
 Square, Oxford OX1 2JD, UK.

 This file is part of Chaste.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
 contributors may be used to endorse or promote products derived from this
 software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#include "VesselSegment.hpp"
#include "MathsCustomFunctions.hpp"
#include "VesselImpedanceCalculator.hpp"
#include "UnitCollection.hpp"
template<unsigned DIM>
VesselImpedanceCalculator<DIM>::VesselImpedanceCalculator() : AbstractVesselNetworkCalculator<DIM>()
{

}

template<unsigned DIM>
VesselImpedanceCalculator<DIM>::~VesselImpedanceCalculator()
{

}

template <unsigned DIM>
std::shared_ptr<VesselImpedanceCalculator<DIM> > VesselImpedanceCalculator<DIM>::Create()
{
    return std::make_shared<VesselImpedanceCalculator<DIM> >();

}

template<unsigned DIM>
void VesselImpedanceCalculator<DIM>::Calculate()
{
    std::vector<std::shared_ptr<VesselSegment<DIM> > > segments = this->mpNetwork->GetVesselSegments();
    for (unsigned idx = 0; idx < segments.size(); idx++)
    {
        QDynamicViscosity viscosity = segments[idx]->GetFlowProperties()->GetViscosity();
        QFlowImpedance impedance = 8.0 * viscosity * segments[idx]->GetLength() / (M_PI * Qpow4(segments[idx]->GetRadius()));
        segments[idx]->GetFlowProperties()->SetImpedance(impedance);
    }
}

template<unsigned DIM>
void VesselImpedanceCalculator<DIM>::CalculateGetLengthFromMatrix()
{
    std::vector<std::shared_ptr<Vessel<DIM> > > vessels = this->mpNetwork->GetVessels();
    for (unsigned idx = 0; idx < vessels.size(); idx++)
    {
        QDynamicViscosity viscosity = vessels[idx]->GetFlowProperties()->GetViscosity();
        QFlowImpedance impedance = 8.0 * viscosity * vessels[idx]->GetLengthFromMatrix() / (M_PI * Qpow4(vessels[idx]->GetRadius()));
        vessels[idx]->GetFlowProperties()->SetImpedance(impedance);
    }
}

template<unsigned DIM>
void VesselImpedanceCalculator<DIM>::CalculateFromMatrix(std::vector<std::vector<double>> rLengthsMatrix)
{
    	std::vector<std::shared_ptr<VesselSegment<DIM> > > segments = this->mpNetwork->GetVesselSegments();
    	unsigned idx = 0;

    	for (unsigned kdx = 0 ; kdx < rLengthsMatrix.size() ; kdx++)
    	{
		for (unsigned jdx = kdx; jdx < rLengthsMatrix[0].size() ; jdx++)
		{
			if (rLengthsMatrix[kdx][jdx] != 0)
			{
				QDynamicViscosity viscosity = segments[idx]->GetFlowProperties()->GetViscosity();

				QLength LocalLength = 0.000001*rLengthsMatrix[kdx][jdx] * unit::metres;
			       	QFlowImpedance impedance = 8.0 * viscosity * LocalLength / (M_PI*Qpow4(segments[idx]->GetRadius()));
			        segments[idx]->GetFlowProperties()->SetImpedance(impedance);			
				idx++;
			}
		}
    	}

}

template<unsigned DIM>
void VesselImpedanceCalculator<DIM>::CalculateFromMatrixWithPruning(std::vector<std::vector<double>> rLengthsMatrix, std::vector<std::vector<double>> rDiametersMatrix, double DimlessRadiusThreshold)
{
    	std::vector<std::shared_ptr<VesselSegment<DIM> > > segments = this->mpNetwork->GetVesselSegments();
    	unsigned idx = 0;

    	for (unsigned kdx = 0 ; kdx < rLengthsMatrix.size() ; kdx++)
    	{
		for (unsigned jdx = kdx; jdx < rLengthsMatrix[0].size() ; jdx++)
		{
			if (rLengthsMatrix[kdx][jdx] != 0 && rDiametersMatrix[kdx][jdx]>=2.0*DimlessRadiusThreshold)
			{
				QDynamicViscosity viscosity = segments[idx]->GetFlowProperties()->GetViscosity();

				QLength LocalLength = 0.000001*rLengthsMatrix[kdx][jdx] * unit::metres;
			       	QFlowImpedance impedance = 8.0 * viscosity * LocalLength / (M_PI*Qpow4(segments[idx]->GetRadius()));
			        segments[idx]->GetFlowProperties()->SetImpedance(impedance);			
				idx++;
			}
		}
    	}

}

template<unsigned DIM>
void VesselImpedanceCalculator<DIM>::CalculateFromMatrixWithPruningByLength(std::vector<std::vector<double>> rLengthsMatrix, double DimlessLengthThreshold)
{
    	std::vector<std::shared_ptr<VesselSegment<DIM> > > segments = this->mpNetwork->GetVesselSegments();
    	unsigned idx = 0;

    	for (unsigned kdx = 0 ; kdx < rLengthsMatrix.size() ; kdx++)
    	{
		for (unsigned jdx = kdx; jdx < rLengthsMatrix[0].size() ; jdx++)
		{
			if (rLengthsMatrix[kdx][jdx]>=DimlessLengthThreshold)
			{
				QDynamicViscosity viscosity = segments[idx]->GetFlowProperties()->GetViscosity();

				QLength LocalLength = 0.000001*rLengthsMatrix[kdx][jdx] * unit::metres;
			       	QFlowImpedance impedance = 8.0 * viscosity * LocalLength / (M_PI*Qpow4(segments[idx]->GetRadius()));
			        segments[idx]->GetFlowProperties()->SetImpedance(impedance);			
				idx++;
			}
		}
    	}

}

// Explicit instantiation
template class VesselImpedanceCalculator<2> ;
template class VesselImpedanceCalculator<3> ;

