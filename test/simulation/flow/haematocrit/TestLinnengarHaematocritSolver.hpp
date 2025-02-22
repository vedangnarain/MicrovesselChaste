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

#ifndef TESTLINNENGARHAEMATOCRITSOLVER_HPP
#define TESTLINNENGARHAEMATOCRITSOLVER_HPP

#include <cxxtest/TestSuite.h>
#include <boost/lexical_cast.hpp>
#include "VesselImpedanceCalculator.hpp"
#include "StructuralAdaptationSolver.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VesselNetworkGenerator.hpp"
#include "FlowSolver.hpp"
#include "SimulationTime.hpp"
#include "LinnengarHaematocritSolver.hpp"
#include "UnitCollection.hpp"
#include "RegularGrid.hpp"
#include "SimulationTime.hpp"
#include "MicrovesselSolver.hpp"
#include "SimpleLinearEllipticFiniteDifferenceSolver.hpp"
#include "VesselBasedDiscreteSource.hpp"
#include "DiscreteContinuumLinearEllipticPde.hpp"
#include "Owen11Parameters.hpp"
#include "Secomb04Parameters.hpp"
#include "GenericParameters.hpp"
#include "StructuralAdaptationSolver.hpp"
#include "ConstantHaematocritSolver.hpp"
#include "BetteridgeHaematocritSolver.hpp"
#include "VesselNetworkPropertyManager.hpp"
#include "VesselNetworkGeometryCalculator.hpp"

#include "PetscAndVtkSetupAndFinalize.hpp"

class TestLinnengarHaematocritSolver : public CxxTest::TestSuite
{

public:


void TestHexagonalNetworkLinnengarHaematocrit()
{
    // Iterate over vessel lengths
    std::vector<double> lengths;
    std::vector<double> average_oxygen_concentration;
    QLength length_increment = 20_um;
    QLength domain_side_length = 2000_um;
    QLength reference_length = 1_um;
    QLength vessel_radius = 10_um;
    QDynamicViscosity visocity = 1.e-3*unit::poiseuille;
    auto p_file_handler =
                    std::make_shared<OutputFileHandler>("TestLinnengarHaematocritSolver_depl", true);

    double inlet_haematocrit = 0.8;
    unsigned num_samples = 0;

    for(unsigned idx=0; idx<num_samples; idx++)
    {
        std::string extension = boost::lexical_cast<std::string>(idx);
        auto p_internal_file_handler =
                        std::make_shared<OutputFileHandler>("TestLinnengarHaematocritSolver_depl/Length_"+extension, true);

        // Specify the network dimensions
        QLength vessel_length = double(idx+1)*length_increment;

        // Generate the network
        VesselNetworkGenerator<2> generator;
        std::shared_ptr<VesselNetwork<2> > p_network = generator.GenerateHexagonalNetwork(domain_side_length,
                domain_side_length, vessel_length);

        // Assign flow properties
        std::vector<VesselNodePtr<2> > nodes;
        nodes.push_back(VesselNodePtr<2> (VesselNode<2>::Create(0_um,5_um)));
        nodes.push_back(VesselNodePtr<2> (VesselNode<2>::Create(5_um,0_um)));
        std::shared_ptr<VesselSegment<2> > p_segment(VesselSegment<2>::Create(nodes[0], nodes[1]));
        p_segment->SetRadius(vessel_radius);
        p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
        VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);
        std::vector<std::shared_ptr<VesselSegment<2> > > segments = p_network->GetVesselSegments();
        for(unsigned jdx=0; jdx<segments.size(); jdx++)
        {
            segments[jdx]->GetFlowProperties()->SetViscosity(visocity);
        }

        VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
                Vertex<2>(0.0_um));
        VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
                Vertex<2>(domain_side_length, domain_side_length));
        p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
        p_inlet_node->GetFlowProperties()->SetPressure(8000.0*unit::pascals);
        p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
        p_outlet_node->GetFlowProperties()->SetPressure(2000.0*unit::pascals);

        std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
        QLength grid_spacing = 5_um;
        p_grid->SetSpacing(grid_spacing);
        c_vector<unsigned, 3> dimensions;
        dimensions[0] = unsigned(domain_side_length/(grid_spacing)) + 1; // num x
        dimensions[1] = unsigned(domain_side_length/(grid_spacing)) + 1; // num_y
        dimensions[2] = 1;
        p_grid->SetDimensions(dimensions);

        /**
         * Next set up the PDEs for oxygen and VEGF. Cells will act as discrete oxygen sinks and discrete vegf sources.
         */
        auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
        p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
        p_oxygen_pde->SetContinuumLinearInUTerm(-10.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

        /**
        * Vessels release oxygen depending on their haematocrit levels
        */
        auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
        QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                GenericParameters::mpGasConcentrationAtStp->GetValue("User");
        QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                Owen11Parameters::mpReferencePartialPressure->GetValue("User");
        p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
        p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
        p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
        p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

        /*
        * Set up a finite difference solver and pass it the pde and grid.
        */
        auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
        p_oxygen_solver->SetPde(p_oxygen_pde);
        p_oxygen_solver->SetLabel("oxygen");
        p_oxygen_solver->SetGrid(p_grid);

        //std::shared_ptr<LinnengarHaematocritSolver<2> > p_haematocrit_calculator = LinnengarHaematocritSolver<2>::Create();

        auto p_haematocrit_calculator = ConstantHaematocritSolver<2>::Create();
        auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
        auto p_structural_adaptation_solver = StructuralAdaptationSolver<2>::Create();
        p_structural_adaptation_solver->SetTolerance(0.0001);
        p_structural_adaptation_solver->SetMaxIterations(1000);
        p_structural_adaptation_solver->SetTimeIncrement(Owen11Parameters::mpVesselRadiusUpdateTimestep->GetValue("User"));
        p_structural_adaptation_solver->AddPreFlowSolveCalculator(p_impedance_calculator);
        p_structural_adaptation_solver->AddPostFlowSolveCalculator(p_haematocrit_calculator);

        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
        std::shared_ptr<MicrovesselSolver<2> > p_microvessel_solver = MicrovesselSolver<2>::Create();
        p_microvessel_solver->SetVesselNetwork(p_network);
        p_microvessel_solver->SetOutputFileHandler(p_internal_file_handler);
        p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
        p_microvessel_solver->SetStructuralAdaptationSolver(p_structural_adaptation_solver);
        p_microvessel_solver->Run();

        std::vector<double> solution = p_oxygen_solver->GetSolution();
        double average = 0.0;
        for(unsigned jdx=0;jdx<solution.size();jdx++)
        {
            average += solution[jdx];
        }
        average /=double(solution.size());
        average_oxygen_concentration.push_back(average);
        lengths.push_back(double(idx+1)*length_increment/reference_length);

        SimulationTime::Instance()->Destroy();
    }

//    std::shared_ptr<std::ofstream> p_out_file = std::shared_ptr<std::ofstream>(new std::ofstream);
//    p_out_file->open((p_file_handler->GetOutputDirectoryFullPath() + "/av_oxygen.dat").c_str());
//    (*p_out_file) << "Length (micron), Av Oxygen (micro M)"<< std::endl;
//    for(unsigned idx=0;idx<average_oxygen_concentration.size();idx++)
//    {
//        (*p_out_file) << lengths[idx] << "," << average_oxygen_concentration[idx] << std::endl;
//    }
//
//    p_out_file->close();

}
};

#endif // TESTLINNENGARHAEMATOCRITSOLVER_HPP
