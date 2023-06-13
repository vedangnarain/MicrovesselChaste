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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Vedang's Notes (vedang.narain@msdtc.ox.ac.uk)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*

This script contains tests that implement eight different haematocrit solvers in four network architectures. 

Outputs can usually be found in /tmp/<home_directory>/testoutput. The pointwise data for the field can be obtained by opening spreadsheet view in Paraview.

H = Haematocrit
BC = boundary condition
RT = Radiotherapy
PQ = Perfusion Quotient
Fig. = Figure
# = number

TO DO LIST

The Alarcon and Pries solvers are breaking the Voronoi simulation. The solvers are currently replaced by Const. H

In the hetreogeneous case, the Yang and Linninger solvers are breaking the dichotomous simulation. The solvers are currently replaced by Const. H

I should probably change the filename to reflect the current program. TestHaematocritSolvers.

Which folder should I put this script in?

For now, pruning can be turned on by specifying the vessels to remove. I need to work on a univesal pruning rule.

For some reason, the output node extends ~0.4 um beyond the x-extent of the simulation. 
Answer: It's because of the grid spacing resolution.

I've marked unknown functionality with a '???'. Need to read up on it later.

For some reason, you need to have an extra node for an unpruned network to have flow. However, single-path networks don't flow if you leave that node in. WHY? I've added a node at a distance of 100 picometres from the last one to circumvent the problem for now.

I should move my equilateral network build to the Network Generator function, once I figure out why that extra node is needed.

Link the inlet haematocrit (0.45) directly with Owen11 parameters.

The y-extent for the hexagonal network is currently set manually based on the simulations. Find a way to make it appropriately cover the hexagonal network automatically.

13/4/21

*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialisation
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef TESTELEMENTARYNETWORK_HPP_
#define TESTELEMENTARYNETWORK_HPP_
#define _BACKWARD_BACKWARD_WARNING_H 1  // Cut out the VTK deprecated warning

// Essential functionality
#include <boost/lexical_cast.hpp>
#include <cxxtest/TestSuite.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "SimulationTime.hpp"
#include "SmartPointers.hpp"

// Geometry tools
#include "Part.hpp"

// Dimensional analysis
#include "BaseUnits.hpp"
#include "UnitCollection.hpp"
#include "GenericParameters.hpp"
#include "ParameterCollection.hpp"
#include "Owen11Parameters.hpp"
#include "Secomb04Parameters.hpp"
#include "Vertex.hpp"

// Grids and PDEs
#include "CellBasedDiscreteSource.hpp"
#include "CellStateDependentDiscreteSource.hpp"
#include "SimpleLinearEllipticFiniteDifferenceSolver.hpp"
#include "SimpleLinearEllipticFiniteElementSolver.hpp"
#include "DiscreteContinuumBoundaryCondition.hpp"
#include "DiscreteContinuumLinearEllipticPde.hpp"
#include "RegularGrid.hpp"
#include "GridCalculator.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"
#include "VesselBasedDiscreteSource.hpp"

// Vessel networks
#include "Vessel.hpp"
#include "VesselNetwork.hpp"
#include "VesselNetworkGenerator.hpp"
#include "VesselNetworkGeometryCalculator.hpp"
#include "VesselNetworkPropertyManager.hpp"
#include "VesselNode.hpp"
#include "VesselSegment.hpp"

// Flow
#include "NodeFlowProperties.hpp"
#include "VesselImpedanceCalculator.hpp"
#include "FlowSolver.hpp"
#include "WallShearStressCalculator.hpp"
#include "MechanicalStimulusCalculator.hpp"
#include "MetabolicStimulusCalculator.hpp"
#include "ShrinkingStimulusCalculator.hpp"
#include "StructuralAdaptationSolver.hpp"
#include "ViscosityCalculator.hpp"

// Haematocrit
#include "AlarconHaematocritSolver.hpp"
#include "BetteridgeHaematocritSolver.hpp"
#include "ConstantHaematocritSolver.hpp"
#include "LinnengarHaematocritSolver.hpp"
#include "YangHaematocritSolver.hpp"
#include "PriesHaematocritSolver.hpp"
#include "PriesWithMemoryHaematocritSolver.hpp"

// Cells
// #include "CancerCellMutationState.hpp"
// #include "StalkCellMutationState.hpp"
// #include "QuiescentCancerCellMutationState.hpp"
// #include "WildTypeCellMutationState.hpp"
// #include "Owen11CellPopulationGenerator.hpp"
// #include "Owen2011TrackingModifier.hpp"
// #include "CaBasedCellPopulation.hpp"
// #include "ApoptoticCellKiller.hpp"
// #include "SimpleOxygenBasedCellCycleModel.hpp"
// #include "StemCellProliferativeType.hpp"

// Forces
#include "GeneralisedLinearSpringForce.hpp"

// Vessel regression solver
#include "WallShearStressBasedRegressionSolver.hpp"

// General solver to collect all the flows
#include "MicrovesselSolver.hpp"
#include "MicrovesselSimulationModifier.hpp"
#include "OnLatticeSimulation.hpp"
#include "OffLatticeSimulation.hpp"

// Visualisation
#include "MicrovesselVtkScene.hpp"
#include "VtkSceneMicrovesselModifier.hpp"

// Keep this last
#include "PetscAndVtkSetupAndFinalize.hpp"

// ???
using namespace std;

// Make a test class
class TestElementaryNetwork : public CxxTest::TestSuite
{

public:

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tests
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define key parameters here, i.e., things that are constant across architectures
    // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));  // = 40 um
    QLength grid_spacing = 10.0_um;  // the simulation time gets quite long if you reduce the resolution further
    // Note: 10.0_um is a rather low resolution for the tiny equilateral network.

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Equilateral Network
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Change resolution of equilateral network for better visualisation
    // QLength equilateral_grid_spacing = grid_spacing;
    QLength equilateral_grid_spacing = 1.0_um;

    // Make an equilateral network on a PDE grid that acts as a Dirichlet BC in 2D 
    void xTestEquilateralNetworkLineSource2D()
    {
        // Run the simulation with different heterogeneities
        for (unsigned n_alpha=0; n_alpha<5; n_alpha++)
        {
            double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 
        
            // Set up the reference length for the simulation
            QLength reference_length(1.0_um);
            BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

            // Set key vessel parameters
            QLength vessel_length(100.0*unit::microns);
            QLength vessel_height = (pow(2,0.5)*vessel_length)*0.5;

            // Set up the domain parameters
            QLength domain_x = vessel_height+vessel_height+vessel_length+vessel_length+0.0001_um;  // this should extend to the x-position of outlet node
            QLength domain_y = domain_x;  // should be the same as domain_x to make square domain
            QLength mid_domain_y = domain_y*0.5;

            // Set nodes based on an equilateral network
            std::shared_ptr<VesselNode<2> > p_node_1 = VesselNode<2>::Create(0.0_um, mid_domain_y);
            std::shared_ptr<VesselNode<2> > p_node_2 = VesselNode<2>::Create(vessel_length, mid_domain_y);
            std::shared_ptr<VesselNode<2> > p_node_3 = VesselNode<2>::Create(vessel_height+vessel_length, vessel_height+mid_domain_y);
            std::shared_ptr<VesselNode<2> > p_node_4 = VesselNode<2>::Create(vessel_height+vessel_length, -vessel_height+mid_domain_y);
            std::shared_ptr<VesselNode<2> > p_node_5 = VesselNode<2>::Create(vessel_height+vessel_height+vessel_length, mid_domain_y);
            std::shared_ptr<VesselNode<2> > p_node_6 = VesselNode<2>::Create(vessel_height+vessel_height+vessel_length+vessel_length, mid_domain_y);
            std::shared_ptr<VesselNode<2> > p_node_7 = VesselNode<2>::Create(domain_x, mid_domain_y);  // add a node at a distance of 100 picometres from the last node to work around glitch of disconnected flow when setting p_node_6 as the output node

            // Make segments 
            std::shared_ptr<VesselSegment<2> > p_segment_1 = VesselSegment<2>::Create(p_node_1, p_node_2);
            std::shared_ptr<VesselSegment<2> > p_segment_2 = VesselSegment<2>::Create(p_node_2, p_node_3);
            std::shared_ptr<VesselSegment<2> > p_segment_3 = VesselSegment<2>::Create(p_node_2, p_node_4);
            std::shared_ptr<VesselSegment<2> > p_segment_4 = VesselSegment<2>::Create(p_node_3, p_node_5);
            std::shared_ptr<VesselSegment<2> > p_segment_5 = VesselSegment<2>::Create(p_node_4, p_node_5);
            std::shared_ptr<VesselSegment<2> > p_segment_6 = VesselSegment<2>::Create(p_node_5, p_node_6);
            std::shared_ptr<VesselSegment<2> > p_segment_7 = VesselSegment<2>::Create(p_node_6, p_node_7);

            // Make vessels
            std::shared_ptr<Vessel<2> > p_vessel_1 = Vessel<2>::Create(p_segment_1);
            std::shared_ptr<Vessel<2> > p_vessel_2 = Vessel<2>::Create(p_segment_2);
            std::shared_ptr<Vessel<2> > p_vessel_3 = Vessel<2>::Create(p_segment_3);
            std::shared_ptr<Vessel<2> > p_vessel_4 = Vessel<2>::Create(p_segment_4);
            std::shared_ptr<Vessel<2> > p_vessel_5 = Vessel<2>::Create(p_segment_5);
            std::shared_ptr<Vessel<2> > p_vessel_6 = Vessel<2>::Create(p_segment_6);
            std::shared_ptr<Vessel<2> > p_vessel_7 = Vessel<2>::Create(p_segment_7);

            // Add the vessels to a vessel network
            std::shared_ptr<VesselNetwork<2> > p_network = VesselNetwork<2>::Create();
            p_network->AddVessel(p_vessel_1);
            p_network->AddVessel(p_vessel_2);
            p_network->AddVessel(p_vessel_3);
            p_network->AddVessel(p_vessel_4);
            p_network->AddVessel(p_vessel_5);
            p_network->AddVessel(p_vessel_6);
            p_network->AddVessel(p_vessel_7);

            // No pruning
            std::ostringstream strs;
            strs << std::fixed << std::setprecision( 1 );
            strs << "TestEquilateralNetwork/DirichletHaematocrit/Alpha" << alpha << "/NoPruning/";
            std::string str_directory_name = strs.str();
            auto p_output_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);         

            // Prune lower path (literally just removes vessels)
            // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestEquilateralNetwork/LineSource2D/VTK/PruneUpperPath");   
            // p_network->RemoveVessel(p_vessel_2,true);
            // p_network->RemoveVessel(p_vessel_4,true);
            // p_network->RemoveVessel(p_vessel_7,true);

            // Prune upper path (literally just removes vessels)
            // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestEquilateralNetwork/LineSource2D/VTK/PruneLowerPath");  
            // p_network->RemoveVessel(p_vessel_3,true);
            // p_network->RemoveVessel(p_vessel_5,true);
            // p_network->RemoveVessel(p_vessel_7,true);
            
            // // Set up the grid
            // std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
            // p_domain->AddRectangle(domain_x, domain_y);
            // std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
            // // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
            // // QLength grid_spacing(10_um);  // the simulation time gets quite long if you reduce the resolution further
            // p_grid->GenerateFromPart(p_domain, equilateral_grid_spacing);  // set domain and grid spacing

            // Set up the grid for the finite difference solver
            auto p_grid = RegularGrid<2>::Create();
            p_grid->SetSpacing(equilateral_grid_spacing);
            c_vector<unsigned, 3> dimensions;
            dimensions[0] = unsigned((domain_x)/(equilateral_grid_spacing))+1; // num x
            dimensions[1] = unsigned((domain_y)/(equilateral_grid_spacing))+1; // num_y
            dimensions[2] = 1;
            p_grid->SetDimensions(dimensions);

            // Choose the PDE
            std::shared_ptr<DiscreteContinuumLinearEllipticPde<2> > p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();

            // Set the diffusivity and decay terms
            p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
            p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));
    
            // Set the O2 concentration value
            QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") * GenericParameters::mpGasConcentrationAtStp->GetValue("User");
            QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp * Owen11Parameters::mpReferencePartialPressure->GetValue("User");

            // Set the boundary condition to be the network acting as a Dirichlet source
            std::shared_ptr<DiscreteContinuumBoundaryCondition<2> > p_vessel_boundary_condition = DiscreteContinuumBoundaryCondition<2>::Create();
            p_vessel_boundary_condition->SetValue(vessel_oxygen_concentration);  // set the O2 concentration at the source
            p_vessel_boundary_condition->SetType(BoundaryConditionType::VESSEL_LINE);
            p_vessel_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

            // Set up the finite difference solver for oxygen (which handles everything)
            SimpleLinearEllipticFiniteDifferenceSolver<2> solver;
            solver.SetGrid(p_grid);
            solver.SetPde(p_oxygen_pde);
            solver.AddBoundaryCondition(p_vessel_boundary_condition);
            solver.SetVesselNetwork(p_network);
            solver.SetLabel("oxygen");
            solver.SetFileName("oxygen_solution_0");

            // Run the solver and write the output file (visualise with Paraview; set Filters->Alphabetical->Tube)
            solver.SetFileHandler(p_output_file_handler);
            solver.SetWriteSolution(true);
            solver.Solve();
            p_network->Write(p_output_file_handler->GetOutputDirectoryFullPath() + "equilateral_network_results.vtp");
    
            // Dump our parameter collection to an xml file and, importantly, clear it for the next test
            ParameterCollection::Instance()->DumpToFile(p_output_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
            ParameterCollection::Instance()->Destroy();
            BaseUnits::Instance()->Destroy();
        }
    }

    // Make an equilateral network on a PDE grid that has flow with different haematocrit splitting rules
    void xTestEquilateralNetworkWithFlow2D()
    {
        // Run the simulation with different heterogeneities
        for (unsigned n_alpha=0; n_alpha<5; n_alpha++)
        {
            double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 
            
            // Run the simulation with different solvers of interest
            for (unsigned h_solver=1; h_solver<8; h_solver++)
            {    
                // Set file name based on haematocrit solver
                std::ostringstream strs;
                strs << std::fixed << std::setprecision( 1 );
                if (h_solver==1)
                {
                    strs << "TestEquilateralNetwork/ConstantHaematocrit/Alpha" << alpha << "/NoPruning/";
                }
                else if (h_solver==2)
                {
                    strs << "TestEquilateralNetwork/PriesHaematocrit/Alpha" << alpha << "/NoPruning/";
                    
                }
                else if (h_solver==3)
                {
                    strs << "TestEquilateralNetwork/MemoryHaematocrit/Alpha" << alpha << "/NoPruning/";
                }
                else if (h_solver==4)
                {
                    strs << "TestEquilateralNetwork/AlarconHaematocrit/Alpha" << alpha << "/NoPruning/";
                }
                else if (h_solver==5)
                {
                    strs << "TestEquilateralNetwork/BetteridgeHaematocrit/Alpha" << alpha << "/NoPruning/";
                }
                else if (h_solver==6)
                {
                    strs << "TestEquilateralNetwork/LinningerHaematocrit/Alpha" << alpha << "/NoPruning/";
                }
                else if (h_solver==7)
                {
                    strs << "TestEquilateralNetwork/YangHaematocrit/Alpha" << alpha << "/NoPruning/";
                }
                std::string str_directory_name = strs.str();
                auto p_output_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                // Set up the reference length for the simulation
                QLength reference_length(1.0_um);
                BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

                // Set key vessel parameters
                QLength vessel_length(100.0_um);
                QLength vessel_height = (pow(2,0.5)*vessel_length)*0.5;

                // Set up the domain parameters
                QLength domain_x = vessel_height+vessel_height+vessel_length+vessel_length+0.0001_um;  // this should be the x-position of outlet node 
                QLength domain_y = domain_x;  // should be the same as domain_x to make square domain
                QLength mid_domain_y = domain_y*0.5;

                // Set nodes based on an equilateral network
                std::shared_ptr<VesselNode<2> > p_node_1 = VesselNode<2>::Create(0.0_um, mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_2 = VesselNode<2>::Create(vessel_length, mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_3 = VesselNode<2>::Create(vessel_height+vessel_length, vessel_height+mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_4 = VesselNode<2>::Create(vessel_height+vessel_length, -vessel_height+mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_5 = VesselNode<2>::Create(vessel_height+vessel_height+vessel_length, mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_6 = VesselNode<2>::Create(vessel_height+vessel_height+vessel_length+vessel_length, mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_7 = VesselNode<2>::Create(domain_x, mid_domain_y);  // add a node at a distance of 100 picometres from the last node to work around glitch of disconnected flow when setting p_node_6 as the output node

                // Make segments 
                std::shared_ptr<VesselSegment<2> > p_segment_1 = VesselSegment<2>::Create(p_node_1, p_node_2);
                std::shared_ptr<VesselSegment<2> > p_segment_2 = VesselSegment<2>::Create(p_node_2, p_node_3);
                std::shared_ptr<VesselSegment<2> > p_segment_3 = VesselSegment<2>::Create(p_node_2, p_node_4);
                std::shared_ptr<VesselSegment<2> > p_segment_4 = VesselSegment<2>::Create(p_node_3, p_node_5);
                std::shared_ptr<VesselSegment<2> > p_segment_5 = VesselSegment<2>::Create(p_node_4, p_node_5);
                std::shared_ptr<VesselSegment<2> > p_segment_6 = VesselSegment<2>::Create(p_node_5, p_node_6);
                std::shared_ptr<VesselSegment<2> > p_segment_7 = VesselSegment<2>::Create(p_node_6, p_node_7);

                // Make vessels
                std::shared_ptr<Vessel<2> > p_vessel_1 = Vessel<2>::Create(p_segment_1);
                std::shared_ptr<Vessel<2> > p_vessel_2 = Vessel<2>::Create(p_segment_2);
                std::shared_ptr<Vessel<2> > p_vessel_3 = Vessel<2>::Create(p_segment_3);
                std::shared_ptr<Vessel<2> > p_vessel_4 = Vessel<2>::Create(p_segment_4);
                std::shared_ptr<Vessel<2> > p_vessel_5 = Vessel<2>::Create(p_segment_5);
                std::shared_ptr<Vessel<2> > p_vessel_6 = Vessel<2>::Create(p_segment_6);
                std::shared_ptr<Vessel<2> > p_vessel_7 = Vessel<2>::Create(p_segment_7);

                // Add the vessels to a vessel network
                std::shared_ptr<VesselNetwork<2> > p_network = VesselNetwork<2>::Create();
                p_network->AddVessel(p_vessel_1);
                p_network->AddVessel(p_vessel_2);
                p_network->AddVessel(p_vessel_3);
                p_network->AddVessel(p_vessel_4);
                p_network->AddVessel(p_vessel_5);
                p_network->AddVessel(p_vessel_6);
                p_network->AddVessel(p_vessel_7);

                // No pruning
                // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestEquilateralNetwork/VTK/NoPruning/");            

                // Prune lower path (literally just removes vessels)
                // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestEquilateralNetwork/VTK/PruneUpperPath");   
                // p_network->RemoveVessel(p_vessel_2,true);
                // p_network->RemoveVessel(p_vessel_4,true);
                // p_network->RemoveVessel(p_vessel_7,true);

                // Prune upper path (literally just removes vessels)
                // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestEquilateralNetwork/VTK/PruneLowerPath");  
                // p_network->RemoveVessel(p_vessel_3,true);
                // p_network->RemoveVessel(p_vessel_5,true);
                // p_network->RemoveVessel(p_vessel_7,true);
                
                // Specify which nodes are the inlets and outlets
                p_network->GetNode(0)->GetFlowProperties()->SetIsInputNode(true);
                p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));
                p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetIsOutputNode(true);
                p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));

                // Set segment radii values
                QLength vessel_radius(1.0 *GenericParameters::mpCapillaryRadius->GetValue());
                VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

                // Set heterogeneous radii for upper path
                QLength alpha_radius(alpha*vessel_radius);
                p_vessel_2->SetRadius(alpha_radius);
                p_vessel_4->SetRadius(alpha_radius);

                // Set segment viscosity values
                QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();
                auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                p_viscosity_calculator->SetVesselNetwork(p_network);
                p_viscosity_calculator->Calculate();
                // VesselNetworkPropertyManager<2>::SetSegmentViscosity(p_network, viscosity);

                // Set up the impedance calculator
                VesselImpedanceCalculator<2> impedance_calculator = VesselImpedanceCalculator<2>();
                impedance_calculator.SetVesselNetwork(p_network);
                impedance_calculator.Calculate();

                // Set up the flow solver
                FlowSolver<2> flow_solver = FlowSolver<2>();
                flow_solver.SetVesselNetwork(p_network);
                flow_solver.SetUp();
                flow_solver.Solve();

                // Set the haematocrit for all vessels
                // AlarconHaematocritSolver<2> haematocrit_solver = AlarconHaematocritSolver<2>();
                // ConstantHaematocritSolver<2> haematocrit_solver = ConstantHaematocritSolver<2>();
                // haematocrit_solver.SetVesselNetwork(p_network);
                // haematocrit_solver.SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                // haematocrit_solver.Calculate();

                // Set the haematocrit solver
                std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                if (h_solver==1)
                {
                    std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;       
                }
                else if (h_solver==2)
                {
                    std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;                    
                }
                else if (h_solver==3)
                {
                    std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;      
                }
                else if (h_solver==4)
                {
                    std::cout << "Now using AlarconHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = AlarconHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;                    
                }
                else if (h_solver==5)
                {
                    std::cout << "Now using BetteridgeHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;      
                }
                else if (h_solver==6)
                {
                    std::cout << "Now using LinnengarHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = LinnengarHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;                    
                }
                else if (h_solver==7)
                {
                    std::cout << "Now using YangHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;      
                }
                // p_abstract_haematocrit_solver->Calculate();

                // // Set up the grid
                // std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
                // p_domain->AddRectangle(domain_x, domain_y);
                // std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
                // // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                // // QLength grid_spacing(10_um);  // the simulation time gets quite long if you reduce the resolution further
                // p_grid->GenerateFromPart(p_domain, equilateral_grid_spacing);  // set domain and grid spacing

                // Set up the grid for the finite difference solver
                auto p_grid = RegularGrid<2>::Create();
                p_grid->SetSpacing(equilateral_grid_spacing);
                c_vector<unsigned, 3> dimensions;
                dimensions[0] = unsigned((domain_x)/(equilateral_grid_spacing))+1; // num x
                dimensions[1] = unsigned((domain_y)/(equilateral_grid_spacing))+1; // num_y
                dimensions[2] = 1;
                p_grid->SetDimensions(dimensions);
                
                // Choose the PDE
                std::shared_ptr<DiscreteContinuumLinearEllipticPde<2> > p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                
                // Set the diffusivity and decay terms
                p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                // Set up the discrete source
                std::shared_ptr<VesselBasedDiscreteSource<2> > p_vessel_source = VesselBasedDiscreteSource<2>::Create();
                QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                p_vessel_source->SetReferenceConcentration(vessel_oxygen_concentration);
                p_vessel_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_vessel_source->SetVesselPermeability(1.0*Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                p_oxygen_pde->AddDiscreteSource(p_vessel_source);

                // Set up the finite difference solver for oxygen (which handles everything)
                auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                p_oxygen_solver->SetPde(p_oxygen_pde);
                p_oxygen_solver->SetLabel("oxygen");
                p_oxygen_solver->SetGrid(p_grid);

                // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                double initial_haematocrit = 0.45;
                unsigned max_iter = 1000; 
                double tolerance2 = 1.e-10;
                std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                for(unsigned idx=0;idx<max_iter;idx++)
                {
                    // Run the solvers
                    impedance_calculator.Calculate();
                    flow_solver.SetUp();
                    flow_solver.Solve();
                    p_abstract_haematocrit_solver->Calculate();
                    p_viscosity_calculator->Calculate();

                    // Get the residual
                    double max_difference = 0.0;
                    double h_for_max = 0.0;
                    double prev_for_max = 0.0;
                    for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                    {
                        double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                        double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                        if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                        {
                            max_difference = difference;
                            h_for_max = current_haematocrit;
                            prev_for_max = previous_haematocrit[jdx];
                        }
                        previous_haematocrit[jdx] = current_haematocrit;
                    }
                    std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                    // Print the final or intermediary convergence results
                    if(max_difference<=tolerance2)  
                    {
                        std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                        break;
                    }
                    else
                    {
                        if(idx%1==0)
                        {
                            std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                            // std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                            // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                            // p_network->Write(output_file);
                        }
                    }

                    // If there is no convergence after all the iterations, print the error message.
                    if(idx==max_iter-1)
                    {
                        std::cout << "Problem encountered in Equilateral Network with h_solver = " << h_solver << " using alpha = " << n_alpha << std::endl;
                        EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                    }
                }
                
                // Run the simulation 
                SimulationTime::Instance()->SetStartTime(0.0);
                SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                p_microvessel_solver->SetVesselNetwork(p_network);
                p_microvessel_solver->SetOutputFileHandler(p_output_file_handler);
                p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                p_microvessel_solver->Run();

                // Print the average oxygenation
                std::vector<double> solution = p_oxygen_solver->GetSolution();
                double average_oxygen = 0.0;
                for(unsigned jdx=0;jdx<solution.size();jdx++)
                {
                    average_oxygen += solution[jdx];
                }
                average_oxygen /= double(solution.size());
                std::cout << "Average oxygen: " << average_oxygen << std::endl;
                
                // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                std::string output_file = p_output_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                p_network->Write(output_file);

                // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                ParameterCollection::Instance()->DumpToFile(p_output_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                ParameterCollection::Instance()->Destroy();
                BaseUnits::Instance()->Destroy();
                SimulationTime::Instance()->Destroy();

                // // Run the solver and write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                // solver.SetFileHandler(p_output_file_handler);
                // solver.SetWriteSolution(true);
                // solver.Solve();
                // p_network->Write(p_output_file_handler->GetOutputDirectoryFullPath() + "equilateral_network_results.vtp");

                // // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                // ParameterCollection::Instance()->DumpToFile(p_output_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                // ParameterCollection::Instance()->Destroy();
                // BaseUnits::Instance()->Destroy();
            }
        }
    }    

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Forking Network
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Make a multi-generation forking network on a PDE grid as a Dirichlet BC in 2D 
    void xTestDichotomousNetworkLineSource2D()
    {
        // Run the simulation with different heterogeneities
        for (unsigned n_alpha=0; n_alpha<5; n_alpha++)
        {
            double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 
            
            // Set up the reference length for the simulation
            QLength reference_length(1.0_um);
            BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

            // Set key vessel parameters
            // QLength vessel_length(100.0_um);
            // QLength vessel_height = (pow(2,0.5)*vessel_length)*0.5;
            unsigned order = 5;  // number of forking generations
            double dimless_length = 1.0;  

            // ???
            for(unsigned i_aux=1; i_aux<order+1; i_aux++)
            {
                dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
            }

            // Set input radius
            // QLength input_radius = 50_um;
            QLength input_radius(10.0 *GenericParameters::mpCapillaryRadius->GetValue());  // = 50_um
            // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

            // Generate the networks
            VesselNetworkGenerator<2> network_generator;
            double lambda;  // lambda = length/diameter
            double twicelambda;  // used as an input parameter
            // for (unsigned k_aux=1; k_aux<5; k_aux++)  // generate the network for various lambdas
            for (unsigned k_aux=1; k_aux<2; k_aux++)  // generate the network for lambda = 4
            {
                lambda = 2.0+double(k_aux)*2.0;
                twicelambda = 2.0*lambda;

                // Height of of first-order vessels
                QLength main_vert_length = 0.9*twicelambda*input_radius*pow(2.0,-1.0/3.0);
                
                // Height of the domain
                QLength domain_side_length_y = 4.0*main_vert_length;

                // Length of the domain
                QLength domain_side_length_x = dimless_length*2.0*twicelambda*input_radius;

                // Generate the name of the directory based on the value of lambda
                std::ostringstream strs;
                strs << std::fixed << std::setprecision( 1 );
                strs << "TestDichotomousNetwork/DirichletHaematocrit/Alpha" << alpha << "/NoPruning/LambdaEquals" << lambda;
                std::string str_directory_name = strs.str();
                auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                // Generate the network
                std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, main_vert_length, input_radius, alpha, twicelambda);

                // Set up the grid for the finite difference solver
                auto p_grid = RegularGrid<2>::Create();
                // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
                // QLength grid_spacing = 10_um;
                p_grid->SetSpacing(grid_spacing);
                c_vector<unsigned, 3> dimensions;
                dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                dimensions[2] = 1;
                p_grid->SetDimensions(dimensions);
                // std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
                // p_domain->AddRectangle(domain_x, domain_y);
                // std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
                // p_grid->GenerateFromPart(p_domain, grid_spacing);  // set domain and grid spacing

                // Choose the PDE
                auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                
                // Set the diffusivity and decay terms
                p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                // Set the O2 concentration value
                QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") * GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp * Owen11Parameters::mpReferencePartialPressure->GetValue("User");

                // Set the boundary condition to be the network acting as a Dirichlet source
                std::shared_ptr<DiscreteContinuumBoundaryCondition<2> > p_vessel_boundary_condition = DiscreteContinuumBoundaryCondition<2>::Create();
                p_vessel_boundary_condition->SetValue(vessel_oxygen_concentration);  // set the O2 concentration at the source
                p_vessel_boundary_condition->SetType(BoundaryConditionType::VESSEL_LINE);
                p_vessel_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);
                    
                // Set up the finite difference solver for oxygen (which handles everything)
                SimpleLinearEllipticFiniteDifferenceSolver<2> solver;
                solver.SetGrid(p_grid);
                solver.SetPde(p_oxygen_pde);
                solver.AddBoundaryCondition(p_vessel_boundary_condition);
                solver.SetVesselNetwork(p_network);
                solver.SetLabel("oxygen");
                solver.SetFileName("oxygen_solution_0");

                // Run the solver and write the output file (visualise with Paraview; set Filters->Alphabetical->Tube)
                solver.SetFileHandler(p_file_handler);
                solver.SetWriteSolution(true);
                solver.Solve();
                p_network->Write(p_file_handler->GetOutputDirectoryFullPath() + "dichotomous_network_results.vtp");
        
                // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                ParameterCollection::Instance()->Destroy();
                BaseUnits::Instance()->Destroy();
            }
        }
    }

    // Make a multi-generation forking network on a PDE grid with different haematocrit splitting rules
    void xTestDichotomousNetworkWithFlow2D()
    {
        // Run the simulation with different heterogeneities
        for (unsigned n_alpha=0; n_alpha<5; n_alpha++)
        {
            double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 
            
            // Run the simulation with different solvers of interest
            for (unsigned h_solver=1; h_solver<8; h_solver++)
            {      
                // Set up the reference length for the simulation
                QLength reference_length(1.0_um);
                BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

                // Set key vessel parameters
                // QLength vessel_length(100.0_um);
                // QLength vessel_height = (pow(2,0.5)*vessel_length)*0.5;
                unsigned order = 5;  // number of forking generations
                double dimless_length = 1.0;  

                // ???
                for(unsigned i_aux=1; i_aux<order+1; i_aux++)
                {
                    dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
                }

                // Set input radius
                QLength input_radius(10.0 *GenericParameters::mpCapillaryRadius->GetValue());  // = 50_um
                // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

                // Set the viscosity
                QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
                // QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();  // = 0.0012

                // Set the inlet and initial haematocrit
                // p_haematocrit_calculator->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));  // = 0.45
                // double inlet_haematocrit = 0.45;
                double initial_haematocrit = 0.45;

                // Generate the networks
                VesselNetworkGenerator<2> network_generator;
                double lambda;  // lambda = length/diameter
                double twicelambda;  // used as an input parameter
                // for (unsigned k_aux=1; k_aux<5; k_aux++)  // generate the network for various lambdas
                for (unsigned k_aux=1; k_aux<2; k_aux++)  // generate the network for lambda = 4
                {
                    lambda = 2.0+double(k_aux)*2.0;
                    twicelambda = 2.0*lambda;

                    // Height of of first-order vessels
                    QLength main_vert_length = 0.9*twicelambda*input_radius*pow(2.0,-1.0/3.0);
                    
                    // Height of the domain
                    QLength domain_side_length_y = 4.0*main_vert_length;

                    // Length of the domain
                    QLength domain_side_length_x = dimless_length*2.0*twicelambda*input_radius;

                    // // Generate the name of the directory based on the value of lambda
                    // std::ostringstream strs;
                    // if (withMemory)
                    // {
                    //     strs << "TestDichotomousNetwork/WithSplitting2D/VTK/NoPruning/Dichotomous_NoCorners_NewModel_LambdaEquals" << lambda;
                    // }
                    // else
                    // {
                    //     strs << "TestDichotomousNetwork/VTK/NoPruning/Dichotomous_NoCorners_OldModel_LambdaEquals" << lambda;
                    // }
                    // std::string str_directory_name = strs.str();
                    // auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);
                        
                    // Set file name based on haematocrit solver
                    std::ostringstream strs;
                    strs << std::fixed << std::setprecision( 1 );
                    if (h_solver==1)
                    {
                        strs << "TestDichotomousNetwork/ConstantHaematocrit/Alpha" << alpha << "/NoPruning/LambdaEquals" << lambda;
                    }
                    else if (h_solver==2)
                    {
                        strs << "TestDichotomousNetwork/PriesHaematocrit/Alpha" << alpha << "/NoPruning/LambdaEquals" << lambda;
                    }
                    else if (h_solver==3)
                    {
                        strs << "TestDichotomousNetwork/MemoryHaematocrit/Alpha" << alpha << "/NoPruning/LambdaEquals" << lambda;

                    }
                    else if (h_solver==4)
                    {
                        strs << "TestDichotomousNetwork/AlarconHaematocrit/Alpha" << alpha << "/NoPruning/LambdaEquals" << lambda;

                    }
                    else if (h_solver==5)
                    {
                        strs << "TestDichotomousNetwork/BetteridgeHaematocrit/Alpha" << alpha << "/NoPruning/LambdaEquals" << lambda;

                    }
                    else if (h_solver==6)
                    {
                        strs << "TestDichotomousNetwork/LinningerHaematocrit/Alpha" << alpha << "/NoPruning/LambdaEquals" << lambda;

                    }
                    else if (h_solver==7)
                    {
                        strs << "TestDichotomousNetwork/YangHaematocrit/Alpha" << alpha << "/NoPruning/LambdaEquals" << lambda;

                    }
                    std::string str_directory_name = strs.str();
                    auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                    // Generate the network
                    std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, main_vert_length, input_radius, alpha, twicelambda);

                    // Identify input and output nodes and assign them properties
                    VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(0.0_um,2.0*main_vert_length));
                    VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(domain_side_length_x, 2.0*main_vert_length));
                    p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
                    // p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));  // = 3333.06
                    p_inlet_node->GetFlowProperties()->SetPressure(3320.0_Pa);
                    p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
                    // p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));  // = 1999.84
                    p_outlet_node->GetFlowProperties()->SetPressure(2090.0_Pa);

                    // Set up the grid for the finite difference solver
                    auto p_grid = RegularGrid<2>::Create();
                    // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                    // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
                    // QLength grid_spacing = 10_um;
                    p_grid->SetSpacing(grid_spacing);
                    c_vector<unsigned, 3> dimensions;
                    dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    dimensions[2] = 1;
                    p_grid->SetDimensions(dimensions);
                    // std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
                    // p_domain->AddRectangle(domain_x, domain_y);
                    // std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
                    // p_grid->GenerateFromPart(p_domain, grid_spacing);  // set domain and grid spacing

                    // Choose the PDE
                    auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // Set the diffusivity and decay terms
                    p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // Set up the discrete source
                    auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                            GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                            Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // Set up the finite difference solver for oxygen (which handles everything)
                    auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    p_oxygen_solver->SetPde(p_oxygen_pde);
                    p_oxygen_solver->SetLabel("oxygen");
                    p_oxygen_solver->SetGrid(p_grid);

                    // // Switch between solvers for Pries or Pries "with memory"
                    // std::shared_ptr<AbstractHaematocritSolver<2>> p_abs_haematocrit_calculator;
                    // if (withMemory)
                    // {
                    //     auto p_haematocrit_calculator = PriesWithMemoryHaematocritSolver<2>::Create();
                    //     p_haematocrit_calculator->SetVesselNetwork(p_network);
                    //     p_haematocrit_calculator->SetHaematocrit(inlet_haematocrit);
                    //     p_abs_haematocrit_calculator = p_haematocrit_calculator;
                    // }
                    // else
                    // {
                    //     auto  p_haematocrit_calculator = PriesHaematocritSolver<2>::Create();
                    //     p_haematocrit_calculator->SetVesselNetwork(p_network);
                    //     p_haematocrit_calculator->SetHaematocrit(inlet_haematocrit);
                    //     p_abs_haematocrit_calculator = p_haematocrit_calculator;
                    // }

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    if (h_solver==1)
                    {
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;       
                    }
                    else if (h_solver==2)
                    {
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                    
                    }
                    else if (h_solver==3)
                    {
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;      
                    }
                    else if (h_solver==4)
                    {
                        std::cout << "Now using AlarconHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = AlarconHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                    
                    }
                    else if (h_solver==5)
                    {
                        std::cout << "Now using BetteridgeHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;      
                    }
                    else if (h_solver==6)
                    {
                        std::cout << "Now using LinnengarHaematocritSolver..." << std::endl;
                        // auto p_haematocrit_solver = LinnengarHaematocritSolver<2>::Create();
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                    
                    }
                    else if (h_solver==7)
                    {
                        std::cout << "Now using YangHaematocritSolver..." << std::endl;
                        // auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;      
                    }

                    // Set up the viscosity calculator
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance calculator
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    flow_solver.SetUp();

                    // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                    unsigned max_iter = 1000; 
                    double tolerance2 = 1.e-10;
                    std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                    std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                    for(unsigned idx=0;idx<max_iter;idx++)
                    {
                        // Run the solvers
                        p_impedance_calculator->Calculate();
                        flow_solver.SetUp();
                        flow_solver.Solve();
                        p_abstract_haematocrit_solver->Calculate();
                        p_viscosity_calculator->Calculate();

                        // Get the residual
                        double max_difference = 0.0;
                        double h_for_max = 0.0;
                        double prev_for_max = 0.0;
                        for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                        {
                            double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                            double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                            if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                            {
                                max_difference = difference;
                                h_for_max = current_haematocrit;
                                prev_for_max = previous_haematocrit[jdx];
                            }
                            previous_haematocrit[jdx] = current_haematocrit;
                        }
                        std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                        // Print the final or intermediary convergence results
                        if(max_difference<=tolerance2)  
                        {
                            std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                            break;
                        }
                        else
                        {
                            if(idx%1==0)
                            {
                                std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                                // std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                                // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                                // p_network->Write(output_file);
                            }
                        }

                        // If there is no convergence after all the iterations, print the error message.
                        if(idx==max_iter-1)
                        {
                            std::cout << "Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using alpha = " << n_alpha << " and lambda = " << lambda << std::endl;
                            EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                        }
                    }
                
                    // Run the simulation 
                    SimulationTime::Instance()->SetStartTime(0.0);
                    SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                    auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                    p_microvessel_solver->SetVesselNetwork(p_network);
                    p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                    p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                    p_microvessel_solver->Run();

                    // Print the average oxygenation
                    std::vector<double> solution = p_oxygen_solver->GetSolution();
                    double average_oxygen = 0.0;
                    for(unsigned jdx=0;jdx<solution.size();jdx++)
                    {
                        average_oxygen += solution[jdx];
                    }
                    average_oxygen /= double(solution.size());
                    std::cout << "Average oxygen: " << average_oxygen << std::endl;
                    
                    // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                    std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                    p_network->Write(output_file);

                    // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                    ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                    ParameterCollection::Instance()->Destroy();
                    BaseUnits::Instance()->Destroy();
                    SimulationTime::Instance()->Destroy();
                }
            }
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Voronoi Network
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Choose the vessel density for the Voronoi network
    unsigned NumberOfSeedPoints = 25;  // change this to select which Voronoi architecture to use: 25, 100, 400

    // // Make a Voronoi network on a PDE grid as a Dirichlet BC in 2D 
    void xTestVoronoiNetworkLineSource2D()
    {
        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Set key parameters
        // unsigned NumberOfSeedPoints = 400;  // change this to select which Voronoi architecture to use: 25, 100, 400
        // QLength input_radius(10.0 *GenericParameters::mpCapillaryRadius->GetValue());  // = 50_um
        // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

        // Initialise the simulation space
        VesselNetworkGenerator<2> network_generator;
        double dimless_domain_size_x = 2000.0; 
        QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
        
        // Initialise the .vtp outputs
        auto p_file_handler = std::make_shared<OutputFileHandler>("TestVoronoiNetwork/DirichletHaematocrit/NoPruning/SeedPoints"+to_string(NumberOfSeedPoints), true);  // create a folder
        // std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");
        std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetwork.vtp");
        // std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("PostRTNetwork.vtp");

        // Set the simulation paramaters
        std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
        // QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;  // minimum flow
        unsigned NumberOfLayouts = 1;  // number of different point layouts to run simulations with (max. = 100)
        std::ofstream outfileMean;
        // unsigned NumberOfCutOffs = 10;  // upper limit for varying pruning length threshold
        // std::vector<std::vector<double> > QuotientMatrixMean(1,std::vector<double>(NumberOfCutOffs+1));
        // std::vector<std::vector<double> > QuotientMatrixAggregate(1,std::vector<double>(NumberOfCutOffs+1,0.0));

        // Run simulations with different layouts (could be used to compute some average later)
        for(unsigned layout=1;layout < NumberOfLayouts+1; layout++)   
        {
            // Read the network layout from a file
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(NumberOfSeedPoints)+"SeedPoints/EdgesMatrixSampleNumber"+to_string(layout)+".txt");
            std::vector<std::vector<double> > rEdgesMatrix;
            string line;

            // Break down the rows in the file into column values
            while (std::getline(in, line)) 
            {
                rEdgesMatrix.push_back(std::vector<double>());
                std::stringstream split(line);
                double value;
                while (split >> value)
                {
                        rEdgesMatrix.back().push_back(value);
                }
            }

            // Generate the network
            std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateVoronoiNetwork(rEdgesMatrix);

            // Set up the grid for the finite difference solver
            auto p_grid = RegularGrid<2>::Create();
            // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
            // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
            // QLength grid_spacing = 10_um;
            p_grid->SetSpacing(grid_spacing);
            c_vector<unsigned, 3> dimensions;
            dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
            dimensions[1] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num_y
            dimensions[2] = 1;
            p_grid->SetDimensions(dimensions);
            // std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
            // p_domain->AddRectangle(domain_x, domain_y);
            // std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
            // p_grid->GenerateFromPart(p_domain, grid_spacing);  // set domain and grid spacing

            // Choose the PDE
            auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
            
            // Set the diffusivity and decay terms
            p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
            p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

            // Set the O2 concentration value
            QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") * GenericParameters::mpGasConcentrationAtStp->GetValue("User");
            QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp * Owen11Parameters::mpReferencePartialPressure->GetValue("User");

            // Set the boundary condition to be the network acting as a Dirichlet source
            std::shared_ptr<DiscreteContinuumBoundaryCondition<2> > p_vessel_boundary_condition = DiscreteContinuumBoundaryCondition<2>::Create();
            p_vessel_boundary_condition->SetValue(vessel_oxygen_concentration);  // set the O2 concentration at the source
            p_vessel_boundary_condition->SetType(BoundaryConditionType::VESSEL_LINE);
            p_vessel_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

            // Set up the finite difference solver for oxygen (which handles everything)
            SimpleLinearEllipticFiniteDifferenceSolver<2> solver;
            solver.SetGrid(p_grid);
            solver.SetPde(p_oxygen_pde);
            solver.AddBoundaryCondition(p_vessel_boundary_condition);
            solver.SetVesselNetwork(p_network);
            solver.SetLabel("oxygen");
            solver.SetFileName("oxygen_solution_0");

            // Run the solver and write the output file (visualise with Paraview; set Filters->Alphabetical->Tube)
            solver.SetFileHandler(p_file_handler);
            solver.SetWriteSolution(true);
            solver.Solve();
            p_network->Write(p_file_handler->GetOutputDirectoryFullPath() + "voronoi_network_results.vtp");
            
            // Dump our parameter collection to an xml file and, importantly, clear it for the next test
            ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
            ParameterCollection::Instance()->Destroy();
            BaseUnits::Instance()->Destroy();
        }
    }

    // Make a 2D Voronoi network on a PDE grid with flow and H-splitting
    void xTestVoronoiNetworkWithFlow2D() 
    {
        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<8; h_solver++)
        {
            // Set up the reference length for the simulation
            QLength reference_length(1.0_um);
            BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

            // Set key parameters
            // unsigned NumberOfSeedPoints = 400;  // change this to select which Voronoi architecture to use: 25, 100, 400
            QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
            // double inlet_haematocrit = 0.45;
            double initial_haematocrit = 0.45;
            // QLength input_radius(10.0 *GenericParameters::mpCapillaryRadius->GetValue());  // = 50_um
            // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);
            
            // Set file name
            std::ostringstream strs;
            if (h_solver==1)
            {
                strs << "TestVoronoiNetwork/ConstantHaematocrit/NoPruning/SeedPoints" << NumberOfSeedPoints;
            }
            else if (h_solver==2)
            {
                strs << "TestVoronoiNetwork/PriesHaematocrit/NoPruning/SeedPoints" << NumberOfSeedPoints;
            }
            else if (h_solver==3)
            {
                strs << "TestVoronoiNetwork/MemoryHaematocrit/NoPruning/SeedPoints" << NumberOfSeedPoints;
            }
            else if (h_solver==4)
            {
                strs << "TestVoronoiNetwork/AlarconHaematocrit/NoPruning/SeedPoints" << NumberOfSeedPoints;
            }
            else if (h_solver==5)
            {
                strs << "TestVoronoiNetwork/BetteridgeHaematocrit/NoPruning/SeedPoints" << NumberOfSeedPoints;
            }
            else if (h_solver==6)
            {
                strs << "TestVoronoiNetwork/LinningerHaematocrit/NoPruning/SeedPoints" << NumberOfSeedPoints;
            }
            else if (h_solver==7)
            {
                strs << "TestVoronoiNetwork/YangHaematocrit/NoPruning/SeedPoints" << NumberOfSeedPoints;
            }
            std::string str_directory_name = strs.str();
            auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

            // Initialise the simulation space
            VesselNetworkGenerator<2> network_generator;
            double dimless_domain_size_x = 2000.0; 
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            
            // Initialise the .vtp outputs
            // auto p_file_handler = std::make_shared<OutputFileHandler>("TestVoronoiNetwork/VTK/NoPruning/Voronoi_WithMethods_ConstPressure_SeedPoints"+to_string(NumberOfSeedPoints), true);  // create a folder
            // std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");
            std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetwork.vtp");
            // std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("PostRTNetwork.vtp");

            // Set the simulation paramaters
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            double tolerance = 0.001;  // tolerance for solution convergence ???
            // QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;  // minimum flow
            unsigned NumberOfLayouts = 1;  // number of different point layouts to run simulations with (max. = 100)
            std::ofstream outfileMean;
            // unsigned NumberOfCutOffs = 10;  // upper limit for varying pruning length threshold
            // std::vector<std::vector<double> > QuotientMatrixMean(1,std::vector<double>(NumberOfCutOffs+1));
            // std::vector<std::vector<double> > QuotientMatrixAggregate(1,std::vector<double>(NumberOfCutOffs+1,0.0));

            // Run simulations with different layouts (could be used to compute some average later)
            for(unsigned layout=1;layout < NumberOfLayouts+1; layout++)   
            {
                // Read the network layout from a file
                std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(NumberOfSeedPoints)+"SeedPoints/EdgesMatrixSampleNumber"+to_string(layout)+".txt");
                std::vector<std::vector<double> > rEdgesMatrix;
                string line;

                // Break down the rows in the file into column values
                while (std::getline(in, line)) 
                {
                    rEdgesMatrix.push_back(std::vector<double>());
                    std::stringstream split(line);
                    double value;
                    while (split >> value)
                    {
                            rEdgesMatrix.back().push_back(value);
                    }
                }

                // Generate the network
                std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateVoronoiNetwork(rEdgesMatrix);

                // Set inlet and outlet nodes
                std::vector<std::shared_ptr<Vessel<2> > > vessels = p_network->GetVessels();
                for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                {
                    if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
                    {
                        if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] <  tolerance)
                        {
                            (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
                            (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
                        }
                        //if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                        else
                        {
                            (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
                            (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                        }
                    }
                    if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
                    {
                        if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] <  tolerance)
                        {
                            (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
                            (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
                        }
                        //if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                        else
                        {
                            (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
                            (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                        }
                    }
                }

                // Set the haematocrit solver
                std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                if (h_solver==1)
                {
                    std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;       
                }
                else if (h_solver==2)
                {
                    std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                    // auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                    auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;                    
                }
                else if (h_solver==3)
                {
                    std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;      
                }
                else if (h_solver==4)
                {
                    std::cout << "Now using AlarconHaematocritSolver..." << std::endl;
                    // auto p_haematocrit_solver = AlarconHaematocritSolver<2>::Create();
                    auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;                    
                }
                else if (h_solver==5)
                {
                    std::cout << "Now using BetteridgeHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;      
                }
                else if (h_solver==6)
                {
                    std::cout << "Now using LinnengarHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = LinnengarHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;                    
                }
                else if (h_solver==7)
                {
                    std::cout << "Now using YangHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;      
                }

                // Set up the viscosity solver
                auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                p_viscosity_calculator->SetVesselNetwork(p_network);
                p_viscosity_calculator->Calculate();
                
                // Set up the impedance solver
                auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                p_impedance_calculator->SetVesselNetwork(p_network);
                p_impedance_calculator->Calculate();
                // if(layout==1)
                // {
                //     p_network->Write(output_file_initial);
                // }
                
                // Set up the flow solver 
                FlowSolver<2> flow_solver;
                flow_solver.SetVesselNetwork(p_network);
                // flow_solver.SetUp();
                flow_solver.SetUseDirectSolver(true);
                // flow_solver.Solve();

                // Set up the grid for the finite difference solver
                auto p_grid = RegularGrid<2>::Create();
                // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
                // QLength grid_spacing = 10_um;
                p_grid->SetSpacing(grid_spacing);
                c_vector<unsigned, 3> dimensions;
                dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                dimensions[1] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num_y
                dimensions[2] = 1;
                p_grid->SetDimensions(dimensions);
                // std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
                // p_domain->AddRectangle(domain_x, domain_y);
                // std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
                // p_grid->GenerateFromPart(p_domain, grid_spacing);  // set domain and grid spacing

                // Choose the PDE
                auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                
                // Set the diffusivity and decay terms
                p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                // Set up the discrete source
                auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                        GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                        Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                // Set up the finite difference solver for oxygen (which handles everything)
                auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                p_oxygen_solver->SetPde(p_oxygen_pde);
                p_oxygen_solver->SetLabel("oxygen");
                p_oxygen_solver->SetGrid(p_grid);
                // solver.SetFileName("oxygen_solution_0");

                // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                unsigned max_iter = 1000;  // 1000 
                double tolerance2 = 1.e-10;
                std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                for(unsigned idx=0;idx<max_iter;idx++)
                {
                    // Run the solvers
                    p_impedance_calculator->Calculate();
                    flow_solver.SetUp();
                    flow_solver.Solve();
                    p_abstract_haematocrit_solver->Calculate();
                    p_viscosity_calculator->Calculate();

                    // Get the residual
                    double max_difference = 0.0;
                    double h_for_max = 0.0;
                    double prev_for_max = 0.0;
                    for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                    {
                        double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                        double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                        if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                        {
                            max_difference = difference;
                            h_for_max = current_haematocrit;
                            prev_for_max = previous_haematocrit[jdx];
                        }
                        previous_haematocrit[jdx] = current_haematocrit;
                    }
                    std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                    // Print the final or intermediary convergence results
                    if(max_difference<=tolerance2)  
                    {
                        std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                        break;
                    }
                    else
                    {
                        if(idx%1==0)
                        {
                            std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                            std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                            std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                            p_network->Write(output_file);
                        }
                    }

                    // If there is no convergence after all the iterations, print the error message.
                    if(idx==max_iter-1)
                    {
                        EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                    }
                }

                // Run the simulation 
                SimulationTime::Instance()->SetStartTime(0.0);
                SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                p_microvessel_solver->SetVesselNetwork(p_network);
                p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                p_microvessel_solver->Run();

                // Print the average oxygenation
                std::vector<double> solution = p_oxygen_solver->GetSolution();
                double average_oxygen = 0.0;
                for(unsigned jdx=0;jdx<solution.size();jdx++)
                {
                    average_oxygen += solution[jdx];
                }
                average_oxygen /= double(solution.size());
                std::cout << "Average oxygen: " << average_oxygen << std::endl;
                
                // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                p_network->Write(output_file);

                // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                ParameterCollection::Instance()->Destroy();
                BaseUnits::Instance()->Destroy();
                SimulationTime::Instance()->Destroy();
            }
        }    
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hexagonal Network 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Make a hexagonal network on a PDE grid as a Dirichlet BC in 2D 
    void xTestHexagonalNetworkLineSource2D()
    {
        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);
        
        // Initialise the simulation space
        double dimless_domain_size_x = 2000.0; 
        double dimless_domain_size_y = 2000.0 - 87.0; 
        QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
        QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
        std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
        std::ofstream outfileMean;

        // Define the key parameters
        unsigned dimless_vessel_length = 100.0;

	    // Read the network from a file
        VesselNetworkGenerator<2> network_generator;
        std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix.txt");
        std::vector<std::vector<double> > rEdgesMatrix;
        string line;
	
        // // Initialise the .vtp outputs
        // auto p_file_handler = std::make_shared<OutputFileHandler>("TestVoronoiNetwork/LineSource2D/NoPruning/SeedPoints"+to_string(NumberOfSeedPoints), true);  // create a folder
        // // std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");
        // std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetwork.vtp");
        // // std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("PostRTNetwork.vtp");

        // Initialise the outputs
        auto p_file_handler = std::make_shared<OutputFileHandler>("TestHexagonalNetwork/DirichletHaematocrit/NoPruning/", true);  // create a folder
        
        // Break down the rows in the file into column values
        while (std::getline(in, line)) 
        {
            rEdgesMatrix.push_back(std::vector<double>());
            std::stringstream split(line);
            double value;
            while (split >> value)
            {
                rEdgesMatrix.back().push_back(value);
            }
        }

        // Set the simulation paramaters
        std::shared_ptr<VesselNetwork<2> > p_network;
        std::vector<std::shared_ptr<Vessel<2> > > vessels;
        // QLength large_vessel_radius = 10_um;
        string line2;
        std::vector<std::vector<unsigned> > Order;
        VesselPtr<2> p_thin_vessel;

        // Generate the network
        p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);

        // Set up the grid for the finite difference solver
        auto p_grid = RegularGrid<2>::Create();
        // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
        // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
        // QLength grid_spacing = 10_um;
        p_grid->SetSpacing(grid_spacing);
        c_vector<unsigned, 3> dimensions;
        dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
        dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
        dimensions[2] = 1;
        p_grid->SetDimensions(dimensions);

        // Choose the PDE
        auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();

        // Set the diffusivity and decay terms
        p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
        p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

        // Set the O2 concentration value
        QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") * GenericParameters::mpGasConcentrationAtStp->GetValue("User");
        QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp * Owen11Parameters::mpReferencePartialPressure->GetValue("User");

        // Set the boundary condition to be the network acting as a Dirichlet source
        std::shared_ptr<DiscreteContinuumBoundaryCondition<2> > p_vessel_boundary_condition = DiscreteContinuumBoundaryCondition<2>::Create();
        p_vessel_boundary_condition->SetValue(vessel_oxygen_concentration);  // set the O2 concentration at the source
        p_vessel_boundary_condition->SetType(BoundaryConditionType::VESSEL_LINE);
        p_vessel_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

        // Set up the finite difference solver for oxygen (which handles everything)
        SimpleLinearEllipticFiniteDifferenceSolver<2> solver;
        solver.SetGrid(p_grid);
        solver.SetPde(p_oxygen_pde);
        solver.AddBoundaryCondition(p_vessel_boundary_condition);
        solver.SetVesselNetwork(p_network);
        solver.SetLabel("oxygen");
        solver.SetFileName("oxygen_solution_0");
        
        // Run the solver and write the output file (visualise with Paraview; set Filters->Alphabetical->Tube)
        solver.SetFileHandler(p_file_handler);
        solver.SetWriteSolution(true);
        solver.Solve();
        p_network->Write(p_file_handler->GetOutputDirectoryFullPath() + "hexagonal_network_results.vtp");
        
        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
        ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
        ParameterCollection::Instance()->Destroy();
        BaseUnits::Instance()->Destroy();
    }

    // Make a 2D Voronoi network on a PDE grid with flow and H-splitting
    void xTestHexagonalNetworkWithFlow2D()
    {
        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<8; h_solver++)
        {
            // Set up the reference length for the simulation
            QLength reference_length(1.0_um);
            BaseUnits::Instance()->SetReferenceLengthScale(reference_length);
            
            // Initialise the simulation space
            double dimless_domain_size_x = 2000.0; 
            double dimless_domain_size_y = 2000.0 - 87.0; 
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            std::ofstream outfileMean;

            // Define the key parameters
            unsigned dimless_vessel_length = 100.0;
            QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
            double initial_haematocrit = 0.45;
	        double tolerance = 0.001;

            // Read the network from a file
            VesselNetworkGenerator<2> network_generator;
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix.txt");
            std::vector<std::vector<double> > rEdgesMatrix;
            string line;

            // Set file name
            std::ostringstream strs;
            if (h_solver==1)
            {
                strs << "TestHexagonalNetwork/ConstantHaematocrit/NoPruning/";
            }
            else if (h_solver==2)
            {
                strs << "TestHexagonalNetwork/PriesHaematocrit/NoPruning/";
            }
            else if (h_solver==3)
            {
                strs << "TestHexagonalNetwork/MemoryHaematocrit/NoPruning/";
            }
            else if (h_solver==4)
            {
                strs << "TestHexagonalNetwork/AlarconHaematocrit/NoPruning/";
            }
            else if (h_solver==5)
            {
                strs << "TestHexagonalNetwork/BetteridgeHaematocrit/NoPruning/";
            }
            else if (h_solver==6)
            {
                strs << "TestHexagonalNetwork/LinningerHaematocrit/NoPruning/";
            }
            else if (h_solver==7)
            {
                strs << "TestHexagonalNetwork/YangHaematocrit/NoPruning/";
            }
            std::string str_directory_name = strs.str();
            auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

            // Break down the rows in the file into column values
            while (std::getline(in, line)) 
            {
                rEdgesMatrix.push_back(std::vector<double>());
                std::stringstream split(line);
                double value;
                while (split >> value)
                {
                    rEdgesMatrix.back().push_back(value);
                }
            }

            // Set the simulation paramaters
            std::shared_ptr<VesselNetwork<2> > p_network;
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            QLength large_vessel_radius = 10_um;
            string line2;
            std::vector<std::vector<unsigned> > Order;
            VesselPtr<2> p_thin_vessel;

            // Generate the network
            p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);
            
            // Set inlet and outlet nodes
            auto p_segment = p_network->GetVesselSegments()[0];
			vessels = p_network->GetVessels();
			for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
			{
				if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
				{
					if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < -0.5*dimless_vessel_length+tolerance)
					{
						(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
						(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
					}
					if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
					{
						(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
						(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
					}
				}
				if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)   	
				//std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
				{
					if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < -0.5*dimless_vessel_length+tolerance)
					{
						(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
						(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
					}
					if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
					{
						(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
						(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
					}
				}
			}
			p_segment = p_network->GetVesselSegments()[0];
			p_segment->SetRadius(large_vessel_radius);
            VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   

            // Set the haematocrit solver
            std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
            if (h_solver==1)
            {
                std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;       
            }
            else if (h_solver==2)
            {
                std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;                    
            }
            else if (h_solver==3)
            {
                std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;      
            }
            else if (h_solver==4)
            {
                std::cout << "Now using AlarconHaematocritSolver..." << std::endl;
                // auto p_haematocrit_solver = AlarconHaematocritSolver<2>::Create();
                auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;                    
            }
            else if (h_solver==5)
            {
                std::cout << "Now using BetteridgeHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;      
            }
            else if (h_solver==6)
            {
                std::cout << "Now using LinnengarHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = LinnengarHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;                    
            }
            else if (h_solver==7)
            {
                std::cout << "Now using YangHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;      
            }

            // Set up the viscosity solver
            auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
            p_viscosity_calculator->SetPlasmaViscosity(viscosity);
            p_viscosity_calculator->SetVesselNetwork(p_network);
            p_viscosity_calculator->Calculate();
            
            // Set up the impedance solver
            auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
            p_impedance_calculator->SetVesselNetwork(p_network);
            p_impedance_calculator->Calculate();

            // Set up the flow solver 
            FlowSolver<2> flow_solver;
            flow_solver.SetVesselNetwork(p_network);
            // flow_solver.SetUp();
            flow_solver.SetUseDirectSolver(true);
            // flow_solver.Solve();

            // Set up the grid for the finite difference solver
            auto p_grid = RegularGrid<2>::Create();
            p_grid->SetSpacing(grid_spacing);
            c_vector<unsigned, 3> dimensions;
            dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
            dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
            dimensions[2] = 1;
            p_grid->SetDimensions(dimensions);

            // Choose the PDE
            auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
            
            // Set the diffusivity and decay terms
            p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
            p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

            // Set up the discrete source
            auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
            QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    GenericParameters::mpGasConcentrationAtStp->GetValue("User");
            QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    Owen11Parameters::mpReferencePartialPressure->GetValue("User");
            p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
            p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
            p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
            p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

            // Set up the finite difference solver for oxygen (which handles everything)
            auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
            p_oxygen_solver->SetPde(p_oxygen_pde);
            p_oxygen_solver->SetLabel("oxygen");
            p_oxygen_solver->SetGrid(p_grid);
            // solver.SetFileName("oxygen_solution_0");

            // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
            unsigned max_iter = 1000;  // 1000 
            double tolerance2 = 1.e-10;
            std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
            std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
            for(unsigned idx=0;idx<max_iter;idx++)
            {
                // Run the solvers
                p_impedance_calculator->Calculate();
                flow_solver.SetUp();
                flow_solver.Solve();
                p_abstract_haematocrit_solver->Calculate();
                p_viscosity_calculator->Calculate();

                // Get the residual
                double max_difference = 0.0;
                double h_for_max = 0.0;
                double prev_for_max = 0.0;
                for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                {
                    double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                    double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                    if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                    {
                        max_difference = difference;
                        h_for_max = current_haematocrit;
                        prev_for_max = previous_haematocrit[jdx];
                    }
                    previous_haematocrit[jdx] = current_haematocrit;
                }
                std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                // Print the final or intermediary convergence results
                if(max_difference<=tolerance2)  
                {
                    std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                    break;
                }
                else
                {
                    if(idx%1==0)
                    {
                        std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                        std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                        p_network->Write(output_file);
                    }
                }

                // If there is no convergence after all the iterations, print the error message.
                if(idx==max_iter-1)
                {
                    EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                }
            }

            // Run the simulation 
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
            auto p_microvessel_solver = MicrovesselSolver<2>::Create();
            p_microvessel_solver->SetVesselNetwork(p_network);
            p_microvessel_solver->SetOutputFileHandler(p_file_handler);
            p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
            p_microvessel_solver->Run();

            // Print the average oxygenation
            std::vector<double> solution = p_oxygen_solver->GetSolution();
            double average_oxygen = 0.0;
            for(unsigned jdx=0;jdx<solution.size();jdx++)
            {
                average_oxygen += solution[jdx];
            }
            average_oxygen /= double(solution.size());
            std::cout << "Average oxygen: " << average_oxygen << std::endl;
            
            // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
            std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
            p_network->Write(output_file);

            // Dump our parameter collection to an xml file and, importantly, clear it for the next test
            ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
            ParameterCollection::Instance()->Destroy();
            BaseUnits::Instance()->Destroy();
            SimulationTime::Instance()->Destroy();
        }
    }

};

#endif /*TESTELEMENTARYNETWORK_HPP_*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Construction Site: Keep Out
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // // Make an equilateral network on a PDE grid that has flow with constant haematocrit and a cell population
    // void xTestEquilateralNetworkWithFlowAndCells2D()
    // {
    //     // Seed the random number generator
    //     RandomNumberGenerator::Instance()->Reseed(12345);

    //     // Set up the reference length and time for the simulation
    //     BaseUnits::Instance()->SetReferenceLengthScale(1.0_um);
    //     BaseUnits::Instance()->SetReferenceTimeScale(1_h);

    //     // Set key vessel parameters
    //     QLength vessel_length(100.0_um);
    //     QLength vessel_height = (pow(2,0.5)*vessel_length)*0.5;

    //     // Set up the domain parameters
    //     QLength domain_x = vessel_height+vessel_height+vessel_length+vessel_length+0.000000000001_um;  // this should be the x-position of outlet node
    //     QLength domain_y = domain_x;  // should be the same as domain_x to make square domain
    //     QLength mid_domain_y = domain_y*0.5;

    //     // Set nodes based on an equilateral network
    //     std::shared_ptr<VesselNode<2> > p_node_1 = VesselNode<2>::Create(0.0_um, mid_domain_y);
    //     std::shared_ptr<VesselNode<2> > p_node_2 = VesselNode<2>::Create(vessel_length, mid_domain_y);
    //     std::shared_ptr<VesselNode<2> > p_node_3 = VesselNode<2>::Create(vessel_height+vessel_length, vessel_height+mid_domain_y);
    //     std::shared_ptr<VesselNode<2> > p_node_4 = VesselNode<2>::Create(vessel_height+vessel_length, -vessel_height+mid_domain_y);
    //     std::shared_ptr<VesselNode<2> > p_node_5 = VesselNode<2>::Create(vessel_height+vessel_height+vessel_length, mid_domain_y);
    //     std::shared_ptr<VesselNode<2> > p_node_6 = VesselNode<2>::Create(vessel_height+vessel_height+vessel_length+vessel_length, mid_domain_y);
    //     std::shared_ptr<VesselNode<2> > p_node_7 = VesselNode<2>::Create(domain_x, mid_domain_y);  // add a node at a distance of one picometre from the last node to work around glitch of disconnected flow when setting p_node_6 as the output node

    //     // Make vessels and vessel segments 
    //     std::shared_ptr<VesselSegment<2> > p_segment_1 = VesselSegment<2>::Create(p_node_1, p_node_2);
    //     std::shared_ptr<VesselSegment<2> > p_segment_2 = VesselSegment<2>::Create(p_node_2, p_node_3);
    //     std::shared_ptr<VesselSegment<2> > p_segment_3 = VesselSegment<2>::Create(p_node_2, p_node_4);
    //     std::shared_ptr<VesselSegment<2> > p_segment_4 = VesselSegment<2>::Create(p_node_3, p_node_5);
    //     std::shared_ptr<VesselSegment<2> > p_segment_5 = VesselSegment<2>::Create(p_node_4, p_node_5);
    //     std::shared_ptr<VesselSegment<2> > p_segment_6 = VesselSegment<2>::Create(p_node_5, p_node_6);
    //     std::shared_ptr<VesselSegment<2> > p_segment_7 = VesselSegment<2>::Create(p_node_6, p_node_7);

    //     // Add vessels to network
    //     std::shared_ptr<Vessel<2> > p_vessel_1 = Vessel<2>::Create(p_segment_1);
    //     std::shared_ptr<Vessel<2> > p_vessel_2 = Vessel<2>::Create(p_segment_2);
    //     std::shared_ptr<Vessel<2> > p_vessel_3 = Vessel<2>::Create(p_segment_3);
    //     std::shared_ptr<Vessel<2> > p_vessel_4 = Vessel<2>::Create(p_segment_4);
    //     std::shared_ptr<Vessel<2> > p_vessel_5 = Vessel<2>::Create(p_segment_5);
    //     std::shared_ptr<Vessel<2> > p_vessel_6 = Vessel<2>::Create(p_segment_6);
    //     std::shared_ptr<Vessel<2> > p_vessel_7 = Vessel<2>::Create(p_segment_7);

    //     // Add the vessels to a vessel network
    //     std::shared_ptr<VesselNetwork<2> > p_network = VesselNetwork<2>::Create();
    //     p_network->AddVessel(p_vessel_1);
    //     p_network->AddVessel(p_vessel_2);
    //     p_network->AddVessel(p_vessel_3);
    //     p_network->AddVessel(p_vessel_4);
    //     p_network->AddVessel(p_vessel_5);
    //     p_network->AddVessel(p_vessel_6);
    //     p_network->AddVessel(p_vessel_7);

    //     // No pruning
    //     // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestElementaryNetwork/WithFlowAndCells2D/VTK/NoPruning/");            

    //     // Prune lower path (literally just removes vessels)
    //     // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestElementaryNetwork/WithFlowAndCells2D/VTK/PruneUpperPath");   
    //     // p_network->RemoveVessel(p_vessel_2,true);
    //     // p_network->RemoveVessel(p_vessel_4,true);
    //     // p_network->RemoveVessel(p_vessel_7,true);

    //     // Prune upper path (literally just removes vessels)
    //     auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestElementaryNetwork/WithFlowAndCells2D/VTK/PruneLowerPath");  
    //     p_network->RemoveVessel(p_vessel_3,true);
    //     p_network->RemoveVessel(p_vessel_5,true);
    //     p_network->RemoveVessel(p_vessel_7,true);
        
    //     // Specify which nodes are the inlets and outlets
    //     p_network->GetNode(0)->GetFlowProperties()->SetIsInputNode(true);
    //     p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));
    //     p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetIsOutputNode(true);
    //     p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));

    //     // Set segment radii and viscosity values
    //     QLength vessel_radius(1.0 *GenericParameters::mpCapillaryRadius->GetValue());
    //     VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);
    //     QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();
    //     VesselNetworkPropertyManager<2>::SetSegmentViscosity(p_network, viscosity);

    //     // Set up the impedance calculator
    //     VesselImpedanceCalculator<2> impedance_calculator = VesselImpedanceCalculator<2>();
    //     impedance_calculator.SetVesselNetwork(p_network);
    //     impedance_calculator.Calculate();

    //     // Set up the flow solver
    //     FlowSolver<2> flow_solver = FlowSolver<2>();
    //     flow_solver.SetVesselNetwork(p_network);
    //     flow_solver.SetUp();
    //     flow_solver.Solve();

    //     // Set the haematocrit for all vessels
    //     // AlarconHaematocritSolver<2> haematocrit_solver = AlarconHaematocritSolver<2>();
    //     ConstantHaematocritSolver<2> haematocrit_solver = ConstantHaematocritSolver<2>();
    //     haematocrit_solver.SetVesselNetwork(p_network);
    //     haematocrit_solver.SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
    //     haematocrit_solver.Calculate();

    //     // Set up the grid
    //     std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
    //     p_domain->AddRectangle(domain_x, domain_y);
    //     std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
    //     // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
    //     QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
    //     p_grid->GenerateFromPart(p_domain, grid_spacing);  // set domain and grid spacing

    //     // We can write the lattice to file for quick viz. with Paraview
    //     p_grid->Write(p_output_file_handler);
    //     std::shared_ptr<MicrovesselVtkScene<2> > p_scene = std::shared_ptr<MicrovesselVtkScene<2> >(new MicrovesselVtkScene<2> );
    //     p_scene->SetRegularGrid(p_grid);
    //     p_scene->GetRegularGridActorGenerator()->SetVolumeOpacity(0.1);
    //     p_scene->SetIsInteractive(true);
    //     p_scene->Start();

    //     // We can write the vessels for quick viz. with Paraview
    //     p_network->Write(p_output_file_handler->GetOutputDirectoryFullPath() + "initial_network.vtp");
    //     p_scene->SetVesselNetwork(p_network);
    //     p_scene->GetVesselNetworkActorGenerator()->SetEdgeSize(20.0);
    //     p_scene->Start();

    //     // Let's set up a grid filled with normal cells 
    //     std::shared_ptr<GridCalculator<2> > p_grid_calc = GridCalculator<2>::Create();
    //     p_grid_calc->SetGrid(p_grid);
    //     std::shared_ptr<Owen11CellPopulationGenerator<2> > p_cell_population_genenerator = Owen11CellPopulationGenerator<2>::Create();
    //     p_cell_population_genenerator->SetGridCalculator(p_grid_calc);
    //     p_cell_population_genenerator->SetVesselNetwork(p_network);
    //     QLength tumour_radius(30.0 * unit::microns);
    //     p_cell_population_genenerator->SetTumourRadius(tumour_radius);
    //     std::shared_ptr<CaBasedCellPopulation<2> > p_cell_population = p_cell_population_genenerator->Update();

    //     // Set up Paraview viz. for the cells
    //     p_scene->GetRegularGridActorGenerator()->SetShowEdges(false);
    //     p_scene->GetRegularGridActorGenerator()->SetVolumeOpacity(0.0);
    //     p_scene->SetCellPopulation(p_cell_population);
    //     p_scene->GetCellPopulationActorGenerator()->GetDiscreteColorTransferFunction()->AddRGBPoint(1.0, 0.0, 0.0, 0.6);
    //     p_scene->GetCellPopulationActorGenerator()->SetPointSize(20.0);
    //     p_scene->GetCellPopulationActorGenerator()->SetColorByCellMutationState(true);
    //     p_scene->ResetRenderer();
    //     p_scene->Start();

    //     // Choose the PDE
    //     std::shared_ptr<DiscreteContinuumLinearEllipticPde<2> > p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
        
    //     // Set the diffusivity term for the oxygen equation
    //     p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
    //     // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

    //     // Set up the cell-based sink
    //     auto p_cell_oxygen_sink = CellBasedDiscreteSource<2>::Create();
    //     p_cell_oxygen_sink->SetLinearInUConsumptionRatePerCell(Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));
    //     p_oxygen_pde->AddDiscreteSource(p_cell_oxygen_sink);

    //     // Set up the discrete source
    //     std::shared_ptr<VesselBasedDiscreteSource<2> > p_vessel_source = VesselBasedDiscreteSource<2>::Create();
    //     QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
    //            GenericParameters::mpGasConcentrationAtStp->GetValue("User");
    //     QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
    //            Owen11Parameters::mpReferencePartialPressure->GetValue("User");
    //     p_vessel_source->SetReferenceConcentration(vessel_oxygen_concentration);
    //     p_vessel_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
    //     p_vessel_source->SetVesselPermeability(1.0*Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
    //     p_oxygen_pde->AddDiscreteSource(p_vessel_source);

    //     // Set up the finite difference solver for oxygen (which handles everything)
    //     // SimpleLinearEllipticFiniteDifferenceSolver<2> solver;
    //     // solver.SetGrid(p_grid);
    //     // solver.SetPde(p_oxygen_pde);
    //     // solver.SetVesselNetwork(p_network);
    //     // //solver.SetUseDirectSolver(false);  // ???
    //     // solver.SetLabel("oxygen");
    //     // solver.SetFileName("oxygen_solution_0");
    //     auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
    //     p_oxygen_solver->SetPde(p_oxygen_pde);
    //     p_oxygen_solver->SetLabel("oxygen");
    //     p_oxygen_solver->SetGrid(p_grid);

    //     // auto p_normal_and_quiescent_cell_source = CellStateDependentDiscreteSource<2>::Create(); //set up source
    //     // std::map<unsigned, QConcentrationFlowRate > normal_and_quiescent_cell_rates; //?
    //     // std::map<unsigned, QConcentration > normal_and_quiescent_cell_rate_thresholds; //?
    //     // MAKE_PTR(QuiescentCancerCellMutationState, p_quiescent_cancer_state);
    //     // MAKE_PTR(WildTypeCellMutationState, p_normal_cell_state);
    //     // normal_and_quiescent_cell_rates[p_normal_cell_state->GetColour()] = Owen11Parameters::mpCellVegfSecretionRate->GetValue("User"); // set vegf secretipn rate for normal 
    //     // normal_and_quiescent_cell_rate_thresholds[p_normal_cell_state->GetColour()] = 0.27*unit::mole_per_metre_cubed; // set threshold for normal cell
    //     // normal_and_quiescent_cell_rates[p_quiescent_cancer_state->GetColour()] = Owen11Parameters::mpCellVegfSecretionRate->GetValue("User"); // set vegf secretion rate for q
    //     // normal_and_quiescent_cell_rate_thresholds[p_quiescent_cancer_state->GetColour()] = 0_M; // set threshold ffor q cell
    //     // p_normal_and_quiescent_cell_source->SetStateRateMap(normal_and_quiescent_cell_rates); // set source rate based on state
    //     // p_normal_and_quiescent_cell_source->SetLabelName("VEGF");
    //     // p_normal_and_quiescent_cell_source->SetStateRateThresholdMap(normal_and_quiescent_cell_rate_thresholds); // set threshp;d
    //     // p_vegf_pde->AddDiscreteSource(p_normal_and_quiescent_cell_source);  // add ousrce

    //     // Set the solver
    //     auto p_microvessel_solver = MicrovesselSolver<2>::Create();
    //     p_microvessel_solver->SetVesselNetwork(p_network);
    //     p_microvessel_solver->SetOutputFrequency(5);
    //     p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
    //     // p_microvessel_solver->SetStructuralAdaptationSolver(p_structural_adaptation_solver);
    //     // p_microvessel_solver->SetRegressionSolver(p_regression_solver);
    //     // p_microvessel_solver->SetAngiogenesisSolver(p_angiogenesis_solver);

    //     // Set up real-time plotting
    //     auto p_scene_modifier = std::make_shared<VtkSceneMicrovesselModifier<2> >();
    //     p_scene_modifier->SetVtkScene(p_scene);
    //     p_scene_modifier->SetUpdateFrequency(2);
    //     p_microvessel_solver->AddMicrovesselModifier(p_scene_modifier);

    //     // Link the vessel and cell solvers
    //     boost::shared_ptr<MicrovesselSimulationModifier<2> > p_microvessel_modifier =
    //             boost::shared_ptr<MicrovesselSimulationModifier<2> >(new MicrovesselSimulationModifier<2> ());
    //     p_microvessel_modifier->SetMicrovesselSolver(p_microvessel_solver);
    //     std::vector<std::string> update_labels;
    //     update_labels.push_back("oxygen");
    //     // update_labels.push_back("VEGF_Extracellular");
    //     p_microvessel_modifier->SetCellDataUpdateLabels(update_labels);

    //     // Set up the cell-based simulation
    //     OnLatticeSimulation<2> simulator(*p_cell_population);
    //     simulator.AddSimulationModifier(p_microvessel_modifier);

    //     // Add killer to remove apoptotic cells
    //     boost::shared_ptr<ApoptoticCellKiller<2> > p_apoptotic_cell_killer(new ApoptoticCellKiller<2>(p_cell_population.get()));
    //     simulator.AddCellKiller(p_apoptotic_cell_killer);

    //     // Add modifier to update cell cycle quantities
    //     boost::shared_ptr<Owen2011TrackingModifier<2> > p_owen11_tracking_modifier(new Owen2011TrackingModifier<2>);
    //     simulator.AddSimulationModifier(p_owen11_tracking_modifier);

    //     // Set up remainder of simulation
    //     simulator.SetOutputDirectory("TestElementaryFlowingCells");
    //     simulator.SetSamplingTimestepMultiple(5);
    //     simulator.SetDt(1.e-6);
    //     simulator.SetEndTime(2.e-6);
    //     simulator.Solve();



    //     // Run the solver and write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
    //     // solver.SetFileHandler(p_output_file_handler);
    //     // solver.SetWriteSolution(true);
    //     // solver.Solve();
    //     p_network->Write(p_output_file_handler->GetOutputDirectoryFullPath() + "equilateral_network_results.vtp");

    //     // Dump our parameter collection to an xml file and, importantly, clear it for the next test
    //     ParameterCollection::Instance()->DumpToFile(p_output_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
    //     ParameterCollection::Instance()->Destroy();
    //     BaseUnits::Instance()->Destroy();
    // }   

    // // Make an equilateral network on a PDE grid that has flow with constant haematocrit and a different sort of cell population
    // void xTestEquilateralNetworkWithUptake2D()
    // {
    //     // // Set up the reference length for the simulation
    //     // QLength reference_length(1.0_um);
    //     // BaseUnits::Instance()->SetReferenceLengthScale(reference_length);
        
    //     //  // Set up the grid
    //     // std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
    //     // p_domain->AddRectangle(domain_x, domain_y);
    //     // std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
    //     // // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
    //     // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
    //     // p_grid->GenerateFromPart(p_domain, grid_spacing);  // set domain and grid spacing

    //     // // Choose the PDE
    //     // std::shared_ptr<DiscreteContinuumLinearEllipticPde<2> > p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
        
    //     // // Set the diffusivity and decay terms
    //     // p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
    //     // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

    //     // // Set up the discrete source
    //     // std::shared_ptr<VesselBasedDiscreteSource<2> > p_vessel_source = VesselBasedDiscreteSource<2>::Create();
    //     // QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
    //     //        GenericParameters::mpGasConcentrationAtStp->GetValue("User");
    //     // QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
    //     //        Owen11Parameters::mpReferencePartialPressure->GetValue("User");
    //     // p_vessel_source->SetReferenceConcentration(vessel_oxygen_concentration);
    //     // p_vessel_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
    //     // p_vessel_source->SetVesselPermeability(1.0*Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
    //     // p_oxygen_pde->AddDiscreteSource(p_vessel_source);
    //     // /*
    //     // // Set up the cell population
    //     // std::vector<CellPtr> cells;

    //     // // Define a mutation state to be used by cells
    //     // MAKE_PTR(WildTypeCellMutationState, p_state);  // healthy cells
    //     // MAKE_PTR(StemCellProliferativeType, p_stem_type);

    //     // for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
    //     // {
    //     //     /*
    //     //      * ...then create a cell, giving it a {{{SimpleOxygenBasedCellCycleModel}}}.
    //     //      * The spatial dimension (1, 2 or 3) needs to be set on the cell-cycle model before it is passed to the cell.
    //     //      */
    //     //     SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
    //     //     p_model->SetDimension(2);
    //     //     CellPtr p_cell(new Cell(p_state, p_model));
    //     //     p_cell->SetCellProliferativeType(p_stem_type);

    //     //     /*
    //     //      * We also alter the default cell-cycle times.
    //     //      */
    //     //     p_model->SetStemCellG1Duration(8.0);
    //     //     p_model->SetTransitCellG1Duration(8.0);

    //     //     /*
    //     //      * We now define a random birth time, chosen from [-T,0], where
    //     //      * T = t,,1,, + t,,2,,, where t,,1,, is a parameter representing the G,,1,, duration
    //     //      * of a 'stem' cell, and t,,2,, is the basic S+G,,2,,+M phases duration...
    //     //      */
    //     //     double birth_time = - RandomNumberGenerator::Instance()->ranf() *
    //     //                          (  p_model->GetStemCellG1Duration()
    //     //                           + p_model->GetSG2MDuration() );
    //     //     /*
    //     //      * ...then we set the birth time and push the cell back into the vector
    //     //      * of cells.
    //     //      */
    //     //     p_cell->SetBirthTime(birth_time);
    //     //     cells.push_back(p_cell);
    //     // }

    //     // // Create cell population
    //     // MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        
        
    //     // // Set key vessel parameters
    //     // QLength vessel_length(100.0_um);
    //     // QLength vessel_height = (pow(2,0.5)*vessel_length)*0.5;

    //     // // Set up the domain parameters
    //     // QLength domain_x = vessel_height+vessel_height+vessel_length+vessel_length+0.000000000001_um;  // this should be the x-position of outlet node
    //     // QLength domain_y = domain_x;  // should be the same as domain_x to make square domain
    //     // QLength mid_domain_y = domain_y*0.5;

    //     // // Set nodes based on an equilateral network
    //     // std::shared_ptr<VesselNode<2> > p_node_1 = VesselNode<2>::Create(0.0_um, mid_domain_y);
    //     // std::shared_ptr<VesselNode<2> > p_node_2 = VesselNode<2>::Create(vessel_length, mid_domain_y);
    //     // std::shared_ptr<VesselNode<2> > p_node_3 = VesselNode<2>::Create(vessel_height+vessel_length, vessel_height+mid_domain_y);
    //     // std::shared_ptr<VesselNode<2> > p_node_4 = VesselNode<2>::Create(vessel_height+vessel_length, -vessel_height+mid_domain_y);
    //     // std::shared_ptr<VesselNode<2> > p_node_5 = VesselNode<2>::Create(vessel_height+vessel_height+vessel_length, mid_domain_y);
    //     // std::shared_ptr<VesselNode<2> > p_node_6 = VesselNode<2>::Create(vessel_height+vessel_height+vessel_length+vessel_length, mid_domain_y);
    //     // std::shared_ptr<VesselNode<2> > p_node_7 = VesselNode<2>::Create(domain_x, mid_domain_y);  // add a node at a distance of one picometre from the last node to work around glitch of disconnected flow when setting p_node_6 as the output node

    //     // // Make vessels and vessel segments 
    //     // std::shared_ptr<VesselSegment<2> > p_segment_1 = VesselSegment<2>::Create(p_node_1, p_node_2);
    //     // std::shared_ptr<VesselSegment<2> > p_segment_2 = VesselSegment<2>::Create(p_node_2, p_node_3);
    //     // std::shared_ptr<VesselSegment<2> > p_segment_3 = VesselSegment<2>::Create(p_node_2, p_node_4);
    //     // std::shared_ptr<VesselSegment<2> > p_segment_4 = VesselSegment<2>::Create(p_node_3, p_node_5);
    //     // std::shared_ptr<VesselSegment<2> > p_segment_5 = VesselSegment<2>::Create(p_node_4, p_node_5);
    //     // std::shared_ptr<VesselSegment<2> > p_segment_6 = VesselSegment<2>::Create(p_node_5, p_node_6);
    //     // std::shared_ptr<VesselSegment<2> > p_segment_7 = VesselSegment<2>::Create(p_node_6, p_node_7);

    //     // // Add vessels to network
    //     // std::shared_ptr<Vessel<2> > p_vessel_1 = Vessel<2>::Create(p_segment_1);
    //     // std::shared_ptr<Vessel<2> > p_vessel_2 = Vessel<2>::Create(p_segment_2);
    //     // std::shared_ptr<Vessel<2> > p_vessel_3 = Vessel<2>::Create(p_segment_3);
    //     // std::shared_ptr<Vessel<2> > p_vessel_4 = Vessel<2>::Create(p_segment_4);
    //     // std::shared_ptr<Vessel<2> > p_vessel_5 = Vessel<2>::Create(p_segment_5);
    //     // std::shared_ptr<Vessel<2> > p_vessel_6 = Vessel<2>::Create(p_segment_6);
    //     // std::shared_ptr<Vessel<2> > p_vessel_7 = Vessel<2>::Create(p_segment_7);

    //     // // Add the vessels to a vessel network
    //     // std::shared_ptr<VesselNetwork<2> > p_network = VesselNetwork<2>::Create();
    //     // p_network->AddVessel(p_vessel_1);
    //     // p_network->AddVessel(p_vessel_2);
    //     // p_network->AddVessel(p_vessel_3);
    //     // p_network->AddVessel(p_vessel_4);
    //     // p_network->AddVessel(p_vessel_5);
    //     // p_network->AddVessel(p_vessel_6);
    //     // p_network->AddVessel(p_vessel_7);

    //     // // No pruning
    //     // // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestElementaryNetwork/VTK/NoPruning/");            

    //     // // Prune lower path (literally just removes vessels)
    //     // // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestElementaryNetwork/VTK/PruneUpperPath");   
    //     // // p_network->RemoveVessel(p_vessel_2,true);
    //     // // p_network->RemoveVessel(p_vessel_4,true);
    //     // // p_network->RemoveVessel(p_vessel_7,true);

    //     // // Prune upper path (literally just removes vessels)
    //     // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestElementaryNetwork/VTK/PruneLowerPath");  
    //     // p_network->RemoveVessel(p_vessel_3,true);
    //     // p_network->RemoveVessel(p_vessel_5,true);
    //     // p_network->RemoveVessel(p_vessel_7,true);
        
    //     // // Specify which nodes are the inlets and outlets
    //     // p_network->GetNode(0)->GetFlowProperties()->SetIsInputNode(true);
    //     // p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));
    //     // p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetIsOutputNode(true);
    //     // p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));

    //     // // Set segment radii and viscosity values
    //     // QLength vessel_radius(1.0 *GenericParameters::mpCapillaryRadius->GetValue());
    //     // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);
    //     // QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();
    //     // VesselNetworkPropertyManager<2>::SetSegmentViscosity(p_network, viscosity);

    //     // // Set up the impedance calculator
    //     // VesselImpedanceCalculator<2> impedance_calculator = VesselImpedanceCalculator<2>();
    //     // impedance_calculator.SetVesselNetwork(p_network);
    //     // impedance_calculator.Calculate();

    //     // // Set up the flow solver
    //     // FlowSolver<2> flow_solver = FlowSolver<2>();
    //     // flow_solver.SetVesselNetwork(p_network);
    //     // flow_solver.SetUp();
    //     // flow_solver.Solve();

    //     // // Set the haematocrit for all vessels
    //     // // AlarconHaematocritSolver<2> haematocrit_solver = AlarconHaematocritSolver<2>();
    //     // ConstantHaematocritSolver<2> haematocrit_solver = ConstantHaematocritSolver<2>();
    //     // haematocrit_solver.SetVesselNetwork(p_network);
    //     // haematocrit_solver.SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
    //     // haematocrit_solver.Calculate();

    //     // // Set up the finite difference solver for oxygen (which handles everything)
    //     // SimpleLinearEllipticFiniteDifferenceSolver<2> solver;
    //     // solver.SetGrid(p_grid);
    //     // solver.SetPde(p_oxygen_pde);
    //     // solver.SetVesselNetwork(p_network);
    //     // //solver.SetUseDirectSolver(false);  // ???
    //     // solver.SetLabel("oxygen");
    //     // solver.SetFileName("oxygen_solution_0");

    //     // // Run the solver and write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
    //     // solver.SetFileHandler(p_output_file_handler);
    //     // solver.SetWriteSolution(true);
    //     // solver.Solve();
    //     // p_network->Write(p_output_file_handler->GetOutputDirectoryFullPath() + "equilateral_network_results.vtp");

    //     // // Dump our parameter collection to an xml file and, importantly, clear it for the next test
    //     // ParameterCollection::Instance()->DumpToFile(p_output_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
    //     // ParameterCollection::Instance()->Destroy();
    //     // BaseUnits::Instance()->Destroy();
    // }   

    // // Make an equilateral network on a PDE grid that acts as a Dirichlet BC in 3D 
    // void xTestEquilateralNetworkLineSource3D()
    // {
    //     // Set up the reference length for the simulation
    //     QLength reference_length(1_um);
    //     BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

    //     // Set key vessel parameters
    //     QLength vessel_length(100.0*unit::microns);
    //     QLength vessel_height = (pow(2,0.5)*vessel_length)*0.5;

    //     // Set up the domain parameters
    //     QLength domain_x = vessel_height+vessel_height+vessel_length+vessel_length;  // this should extend to the x-position of outlet node
    //     QLength domain_y = domain_x;  // should be the same as domain_x to make square domain
    //     QLength domain_z(0_um);  // since we're working in 2D
    //     QLength mid_domain_y = domain_y*0.5;

    //     // Set nodes based on an equilateral network
    //     std::shared_ptr<VesselNode<3> > p_node_1 = VesselNode<3>::Create(0.0_um, mid_domain_y, domain_z);
    //     std::shared_ptr<VesselNode<3> > p_node_2 = VesselNode<3>::Create(vessel_length, mid_domain_y, domain_z);
    //     std::shared_ptr<VesselNode<3> > p_node_3 = VesselNode<3>::Create(vessel_height+vessel_length, vessel_height+mid_domain_y, domain_z);
    //     std::shared_ptr<VesselNode<3> > p_node_4 = VesselNode<3>::Create(vessel_height+vessel_length, -vessel_height+mid_domain_y, domain_z);
    //     std::shared_ptr<VesselNode<3> > p_node_5 = VesselNode<3>::Create(vessel_height+vessel_height+vessel_length, mid_domain_y, domain_z);
    //     std::shared_ptr<VesselNode<3> > p_node_6 = VesselNode<3>::Create(vessel_height+vessel_height+vessel_length+vessel_length, mid_domain_y, domain_z);

    //     // Make vessels and vessel segments 
    //     std::shared_ptr<VesselSegment<3> > p_segment_1 = VesselSegment<3>::Create(p_node_1, p_node_2);
    //     std::shared_ptr<Vessel<3> > p_vessel_1 = Vessel<3>::Create(p_segment_1);
    //     std::shared_ptr<VesselSegment<3> > p_segment_2 = VesselSegment<3>::Create(p_node_2, p_node_3);
    //     std::shared_ptr<Vessel<3> > p_vessel_2 = Vessel<3>::Create(p_segment_2);
    //     std::shared_ptr<VesselSegment<3> > p_segment_3 = VesselSegment<3>::Create(p_node_2, p_node_4);
    //     std::shared_ptr<Vessel<3> > p_vessel_3 = Vessel<3>::Create(p_segment_3);
    //     std::shared_ptr<VesselSegment<3> > p_segment_4 = VesselSegment<3>::Create(p_node_3, p_node_5);
    //     std::shared_ptr<Vessel<3> > p_vessel_4 = Vessel<3>::Create(p_segment_4);
    //     std::shared_ptr<VesselSegment<3> > p_segment_5 = VesselSegment<3>::Create(p_node_4, p_node_5);
    //     std::shared_ptr<Vessel<3> > p_vessel_5 = Vessel<3>::Create(p_segment_5);
    //     std::shared_ptr<VesselSegment<3> > p_segment_6 = VesselSegment<3>::Create(p_node_5, p_node_6);
    //     std::shared_ptr<Vessel<3> > p_vessel_6 = Vessel<3>::Create(p_segment_6);

    //     // Add the vessels to a vessel network
    //     std::shared_ptr<VesselNetwork<3> > p_network = VesselNetwork<3>::Create();
    //     p_network->AddVessel(p_vessel_1);
    //     p_network->AddVessel(p_vessel_2);
    //     p_network->AddVessel(p_vessel_3);
    //     p_network->AddVessel(p_vessel_4);
    //     p_network->AddVessel(p_vessel_5);
    //     p_network->AddVessel(p_vessel_6);

    //     // Simple pruning (literally just removes vessels)
    //     p_network->RemoveVessel(p_vessel_2,true);
    //     p_network->RemoveVessel(p_vessel_4,true);

    //     // Set up the grid
    //     std::shared_ptr<Part<3> > p_domain = Part<3>::Create();
    //     p_domain->AddCuboid(domain_x, domain_y, domain_z, Vertex<3>(0.0, 0.0, 0.0));
    //     std::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
    //     p_grid->GenerateFromPart(p_domain, 1.0*1_um);  // set domain and grid spacing

    //     // Set up the PDE expression for oxygen
    //     std::shared_ptr<DiscreteContinuumLinearEllipticPde<3> > p_oxygen_pde = DiscreteContinuumLinearEllipticPde<3>::Create();
    //     // QDiffusivity diffusivity(0.0033 * unit::metre_squared_per_second);     
    //     // QRate consumption_rate(2.e6 * unit::per_second);  
    //     // p_oxygen_pde->SetIsotropicDiffusionConstant(diffusivity);
    //     // p_oxygen_pde->SetContinuumLinearInUTerm(-consumption_rate);
    //     p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
    //     p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));
 
    //     // Set the O2 concentration value
    //     // QConcentration boundary_concentration(1.e-3 * unit::mole_per_metre_cubed);  // set the value directly
    //     QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") * GenericParameters::mpGasConcentrationAtStp->GetValue("User");
    //     QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp * Owen11Parameters::mpReferencePartialPressure->GetValue("User");

    //     // Set the boundary condition to be the network acting as a Dirichlet source
    //     std::shared_ptr<DiscreteContinuumBoundaryCondition<3> > p_vessel_boundary_condition = DiscreteContinuumBoundaryCondition<3>::Create();
    //     p_vessel_boundary_condition->SetValue(vessel_oxygen_concentration);  // set the O2 concentration at the source
    //     p_vessel_boundary_condition->SetType(BoundaryConditionType::VESSEL_LINE);
    //     p_vessel_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

    //     // Set up the solver
    //     SimpleLinearEllipticFiniteDifferenceSolver<3> solver;
    //     solver.SetGrid(p_grid);
    //     solver.SetPde(p_oxygen_pde);
    //     solver.AddBoundaryCondition(p_vessel_boundary_condition);
    //     solver.SetVesselNetwork(p_network);
    //     solver.SetLabel("oxygen");
    //     solver.SetFileName("oxygen_solution_0");

    //     // Set up output file management and seed the random number generator.
    //     auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestElementaryNetwork/WithoutFlow", true);        
    //     solver.SetFileHandler(p_output_file_handler);
        
    //     // Run the solver and write the output file (visualise with Paraview; set Filters->Alphabetical->Tube)
    //     solver.SetWriteSolution(true);
    //     solver.Solve();
    //     p_network->Write(p_output_file_handler->GetOutputDirectoryFullPath() + "equilateral_network_results.vtp");
 
    //     // Dump our parameter collection to an xml file and, importantly, clear it for the next test
    //     ParameterCollection::Instance()->DumpToFile(p_output_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
    //     ParameterCollection::Instance()->Destroy();
    //     BaseUnits::Instance()->Destroy();
    // }
    
    // // Make an equilateral network on a PDE grid that has flow with constant haematocrit (broken, doesn't work in pseudo-3D for some reason. Delete after turning previous test to 3D)
    // void xTestEquilateralNetworkLineSourceWithFlowBroken3D()
    // {
    //     // Set up the reference length for the simulation
    //     QLength reference_length(1.0 * unit::microns);
    //     BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

    //     // Set key vessel parameters
    //     QLength vessel_length(100.0*unit::microns);
    //     QLength vessel_height = (pow(2,0.5)*vessel_length)*0.5;

    //     // Set up the domain parameters
    //     QLength domain_x = vessel_height+vessel_height+vessel_length+vessel_length;  // this should extend to the x-position of outlet node
    //     QLength domain_y = domain_x;  // should be the same as domain_x to make square domain
    //     QLength domain_z(1_um);  // since we're working in 2D
    //     QLength mid_domain_y = domain_y*0.5;

    //     // Set nodes based on an equilateral network
    //     std::shared_ptr<VesselNode<3> > p_node_1 = VesselNode<3>::Create(0.0_um, mid_domain_y, domain_z);
    //     std::shared_ptr<VesselNode<3> > p_node_2 = VesselNode<3>::Create(vessel_length, mid_domain_y, domain_z);
    //     std::shared_ptr<VesselNode<3> > p_node_3 = VesselNode<3>::Create(vessel_height+vessel_length, vessel_height+mid_domain_y, domain_z);
    //     std::shared_ptr<VesselNode<3> > p_node_4 = VesselNode<3>::Create(vessel_height+vessel_length, -vessel_height+mid_domain_y, domain_z);
    //     std::shared_ptr<VesselNode<3> > p_node_5 = VesselNode<3>::Create(vessel_height+vessel_height+vessel_length, mid_domain_y, domain_z);
    //     std::shared_ptr<VesselNode<3> > p_node_6 = VesselNode<3>::Create(vessel_height+vessel_height+vessel_length+vessel_length, mid_domain_y, domain_z);

    //     // Make vessels and vessel segments 
    //     std::shared_ptr<VesselSegment<3> > p_segment_1 = VesselSegment<3>::Create(p_node_1, p_node_2);
    //     std::shared_ptr<Vessel<3> > p_vessel_1 = Vessel<3>::Create(p_segment_1);
    //     std::shared_ptr<VesselSegment<3> > p_segment_2 = VesselSegment<3>::Create(p_node_2, p_node_3);
    //     std::shared_ptr<Vessel<3> > p_vessel_2 = Vessel<3>::Create(p_segment_2);
    //     std::shared_ptr<VesselSegment<3> > p_segment_3 = VesselSegment<3>::Create(p_node_2, p_node_4);
    //     std::shared_ptr<Vessel<3> > p_vessel_3 = Vessel<3>::Create(p_segment_3);
    //     std::shared_ptr<VesselSegment<3> > p_segment_4 = VesselSegment<3>::Create(p_node_3, p_node_5);
    //     std::shared_ptr<Vessel<3> > p_vessel_4 = Vessel<3>::Create(p_segment_4);
    //     std::shared_ptr<VesselSegment<3> > p_segment_5 = VesselSegment<3>::Create(p_node_4, p_node_5);
    //     std::shared_ptr<Vessel<3> > p_vessel_5 = Vessel<3>::Create(p_segment_5);
    //     std::shared_ptr<VesselSegment<3> > p_segment_6 = VesselSegment<3>::Create(p_node_5, p_node_6);
    //     std::shared_ptr<Vessel<3> > p_vessel_6 = Vessel<3>::Create(p_segment_6);

    //     // Add the vessels to a vessel network
    //     std::shared_ptr<VesselNetwork<3> > p_network = VesselNetwork<3>::Create();
    //     p_network->AddVessel(p_vessel_1);
    //     p_network->AddVessel(p_vessel_2);
    //     p_network->AddVessel(p_vessel_3);
    //     p_network->AddVessel(p_vessel_4);
    //     p_network->AddVessel(p_vessel_5);
    //     p_network->AddVessel(p_vessel_6);

    //     // Simple pruning (literally just removes vessels)
    //     p_network->RemoveVessel(p_vessel_2,true);
    //     p_network->RemoveVessel(p_vessel_4,true);

    //     ///////////////////////////////////////////////////////////

    //     // Specify which nodes are the inlets and outlets
    //     p_network->GetNode(0)->GetFlowProperties()->SetIsInputNode(true);
    //     p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));
    //     p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetIsOutputNode(true);
    //     p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));

    //     // Set segment radii and viscosity values
    //     QLength vessel_radius(GenericParameters::mpCapillaryRadius->GetValue());
    //     VesselNetworkPropertyManager<3>::SetSegmentRadii(p_network, vessel_radius);
    //     // p_network->SetSegmentRadii(vessel_radius);
    //     QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();
    //     // p_network->SetSegmentViscosity(viscosity);
    //     VesselNetworkPropertyManager<3>::SetSegmentViscosity(p_network, viscosity);

    //     // Set up the impedance calculator
    //     VesselImpedanceCalculator<3> impedance_calculator = VesselImpedanceCalculator<3>();
    //     impedance_calculator.SetVesselNetwork(p_network);
    //     impedance_calculator.Calculate();

    //     // Check the impedance
    //     // units::quantity<unit::flow_impedance> expected_impedance = 8.0 * viscosity* vessel_length/(M_PI*units::pow<4>(vessel_radius));
    //     // TS_ASSERT_DELTA(p_network->GetVessel(0)->GetSegment(0)->GetFlowProperties()->GetImpedance().value(), expected_impedance.value(), 1.e-6);

    //     // Set up the flow solver
    //     FlowSolver<3> flow_solver = FlowSolver<3>();
    //     flow_solver.SetVesselNetwork(p_network);
    //     flow_solver.Solve();

    //     // Check the pressure: it is expected to drop linearly so should be the average of the input and output half way along the network
    //     // units::quantity<unit::pressure> expected_pressure = (inlet_pressure + Owen11Parameters::mpOutletPressure->GetValue())/2.0;
    //     // TS_ASSERT_DELTA(p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->GetPressure().value(), expected_pressure.value(), 1.e-6);

    //     // Set the haematocrit for all vessels
    //     // AlarconHaematocritSolver<3> haematocrit_solver = AlarconHaematocritSolver<3>();
    //     ConstantHaematocritSolver<3> haematocrit_solver = ConstantHaematocritSolver<3>();
    //     haematocrit_solver.SetVesselNetwork(p_network);
    //     haematocrit_solver.SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
    //     haematocrit_solver.Calculate();
    //     // p_network->GetVessels()[0]->GetFlowProperties()->SetHaematocrit(0.45);



    //     ///////////////////////////////////////////////////////////

    //     // Set up the grid
    //     std::shared_ptr<Part<3> > p_domain = Part<3>::Create();
    //     p_domain->AddCuboid(domain_x, domain_y, domain_z, Vertex<3>(0.0, 0.0, 0.0));
    //     std::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
    //     p_grid->GenerateFromPart(p_domain, 1.0*1_um);  // set domain and grid spacing

    //     // Alternate way to set up grid
    //     // auto p_grid = RegularGrid<3>::Create();
    //     // QLength grid_spacing = Owen11Parameters::mpLatticeSpacing->GetValue("User");
    //     // p_grid->SetSpacing(grid_spacing);
    //     // c_vector<unsigned, 3> dimensions;
    //     // dimensions[0] = 51; // num x
    //     // dimensions[1] = 51; // num_y
    //     // dimensions[2] = 1;
    //     // p_grid->SetDimensions(dimensions);

    //     // Set up the PDE expression for oxygen (diffusivity and decay term)
    //     std::shared_ptr<DiscreteContinuumLinearEllipticPde<3> > p_oxygen_pde = DiscreteContinuumLinearEllipticPde<3>::Create();
    //     p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
    //     p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));
    //     // QConcentrationFlowRate consumption_rate(-2.e0 * unit::mole_per_metre_cubed_per_second);
    //     // p_oxygen_pde->SetContinuumConstantInUTerm(consumption_rate);
 
    //     // Set the O2 concentration value
    //     QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") * GenericParameters::mpGasConcentrationAtStp->GetValue("User");
    //     QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp * Owen11Parameters::mpReferencePartialPressure->GetValue("User");
        
    //     // Set up the haematocrit-based vessel source for the PDE  //
    //     auto p_vessel_oxygen_source = VesselBasedDiscreteSource<3>::Create();
    //     p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration); 
    //     p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
    //     p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
    //     p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);   // !!! only problematic line in whole section

    //     // Set the boundary condition to be the network acting as a Dirichlet source
    //     // std::shared_ptr<DiscreteContinuumBoundaryCondition<3> > p_vessel_boundary_condition = DiscreteContinuumBoundaryCondition<3>::Create();
    //     // p_vessel_boundary_condition->SetValue(vessel_oxygen_concentration);  // set the O2 concentration at the source
    //     // p_vessel_boundary_condition->SetType(BoundaryConditionType::VESSEL_LINE);
    //     // p_vessel_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

    //     // Set up the finite difference solver for oxygen (which handles everything)
    //     SimpleLinearEllipticFiniteDifferenceSolver<3> solver;
    //     solver.SetGrid(p_grid);  // add the oxygen solver to the grid 
    //     solver.SetPde(p_oxygen_pde);
    //     // solver.AddBoundaryCondition(p_vessel_boundary_condition);
    //     solver.SetVesselNetwork(p_network);
    //     solver.SetLabel("oxygen");
    //     solver.SetFileName("oxygen_solution_0");

    //     // Set up output file management and seed the random number generator.
    //     auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestElementaryNetwork/WithFlow", true);        
    //     solver.SetFileHandler(p_output_file_handler);
        
    //     // Run the solver and write the output file (visualise with Paraview; set Filters->Alphabetical->Tube)
    //     solver.SetWriteSolution(true);
    //     solver.Solve();
    //     p_network->Write(p_output_file_handler->GetOutputDirectoryFullPath() + "equilateral_network_results.vtp");
 
    //     // Dump our parameter collection to an xml file and, importantly, clear it for the next test
    //     ParameterCollection::Instance()->DumpToFile(p_output_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
    //     ParameterCollection::Instance()->Destroy();
    //     BaseUnits::Instance()->Destroy();
    // }   

    // // Make an equilateral network on a PDE grid that acts as a Dirichlet BC with a Universal Pruning Rule (not functional yet)
    // void xTestEquilateralNetworkLineSourceUniversalPruning()
    // {
    //     // Set up the reference length for the simulation
    //     QLength reference_length(1_um);
    //     BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

    //     // Set key vessel parameters
    //     QLength vessel_length(100.0*unit::microns);
    //     QLength vessel_height = (pow(2,0.5)*vessel_length)*0.5;

    //     // Set up the domain parameters
    //     QLength domain_x = vessel_height+vessel_height+vessel_length+vessel_length;  // this should extend to the x-position of outlet node
    //     QLength domain_y = domain_x;  // should be the same as domain_x to make square domain
    //     QLength domain_z(1_um);  // since we're working in 2D
    //     QLength mid_domain_y = domain_y*0.5;

    //     // Set nodes based on an equilateral network
    //     std::shared_ptr<VesselNode<3> > p_node_1 = VesselNode<3>::Create(0.0_um, mid_domain_y, domain_z);
    //     std::shared_ptr<VesselNode<3> > p_node_2 = VesselNode<3>::Create(vessel_length, mid_domain_y, domain_z);
    //     std::shared_ptr<VesselNode<3> > p_node_3 = VesselNode<3>::Create(vessel_height+vessel_length, vessel_height+mid_domain_y, domain_z);
    //     std::shared_ptr<VesselNode<3> > p_node_4 = VesselNode<3>::Create(vessel_height+vessel_length, -vessel_height+mid_domain_y, domain_z);
    //     std::shared_ptr<VesselNode<3> > p_node_5 = VesselNode<3>::Create(vessel_height+vessel_height+vessel_length, mid_domain_y, domain_z);
    //     std::shared_ptr<VesselNode<3> > p_node_6 = VesselNode<3>::Create(vessel_height+vessel_height+vessel_length+vessel_length, mid_domain_y, domain_z);

    //     // Make vessels and vessel segments 
    //     std::shared_ptr<VesselSegment<3> > p_segment_1 = VesselSegment<3>::Create(p_node_1, p_node_2);
    //     std::shared_ptr<Vessel<3> > p_vessel_1 = Vessel<3>::Create(p_segment_1);
    //     std::shared_ptr<VesselSegment<3> > p_segment_2 = VesselSegment<3>::Create(p_node_2, p_node_3);
    //     std::shared_ptr<Vessel<3> > p_vessel_2 = Vessel<3>::Create(p_segment_2);
    //     std::shared_ptr<VesselSegment<3> > p_segment_3 = VesselSegment<3>::Create(p_node_2, p_node_4);
    //     std::shared_ptr<Vessel<3> > p_vessel_3 = Vessel<3>::Create(p_segment_3);
    //     std::shared_ptr<VesselSegment<3> > p_segment_4 = VesselSegment<3>::Create(p_node_3, p_node_5);
    //     std::shared_ptr<Vessel<3> > p_vessel_4 = Vessel<3>::Create(p_segment_4);
    //     std::shared_ptr<VesselSegment<3> > p_segment_5 = VesselSegment<3>::Create(p_node_4, p_node_5);
    //     std::shared_ptr<Vessel<3> > p_vessel_5 = Vessel<3>::Create(p_segment_5);
    //     std::shared_ptr<VesselSegment<3> > p_segment_6 = VesselSegment<3>::Create(p_node_5, p_node_6);
    //     std::shared_ptr<Vessel<3> > p_vessel_6 = Vessel<3>::Create(p_segment_6);

    //     // Add the vessels to a vessel network
    //     std::shared_ptr<VesselNetwork<3> > p_network = VesselNetwork<3>::Create();
    //     p_network->AddVessel(p_vessel_1);
    //     p_network->AddVessel(p_vessel_2);
    //     p_network->AddVessel(p_vessel_3);
    //     p_network->AddVessel(p_vessel_4);
    //     p_network->AddVessel(p_vessel_5);
    //     p_network->AddVessel(p_vessel_6);

    //     // Simple pruning (literally just removes vessels)
    //     // p_network->RemoveVessel(p_vessel_2,true);
    //     // p_network->RemoveVessel(p_vessel_4,true);

    //     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //     // PRUNING TEST AREA STARTS
    //     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          

    //     // void VesselNetwork<DIM>::RemoveThinVessels(QLength cutoffRadius, bool endsOnly)
    //     // {
    //     //     std::vector<std::shared_ptr<Vessel<DIM> > > vessels_to_remove;  

    //     //     for(unsigned idx=0; idx<mVessels.size(); idx++)
    //     //     {
    //     //         if(mVessels[idx]->GetRadius() < cutoffRadius)
    //     //         {
    //     //             if(endsOnly && (mVessels[idx]->GetStartNode()->GetNumberOfSegments() == 1 || mVessels[idx]->GetEndNode()->GetNumberOfSegments() == 1 ))
    //     //             {
    //     //                 vessels_to_remove.push_back(mVessels[idx]);
    //     //             }
    //     //             else if(!endsOnly)
    //     //             {
    //     //                 vessels_to_remove.push_back(mVessels[idx]);
    //     //             }
    //     //         }
    //     //     }

    //     //     // remove all vessels in list of condemned
    //     //     for(unsigned idx=0; idx<vessels_to_remove.size(); idx++)
    //     //     {
    //     //         RemoveVessel(vessels_to_remove[idx],true);
    //     //     }
    //     // }

    //               //auto minRadius = (*vessels.begin())->GetRadius();
    //       //unsigned IdToKill = 0;
    //     /*
    //       for(it = vessels.begin(); it != vessels.end(); it++)
    //       {	
    //     if  ((minRadius > (*it)->GetRadius())&&((*it)->GetRadius()>0.00000001_um))
    //     {
    //       minRadius = (*it)->GetRadius();
    //       IdToKill = (*it)->GetId();
    //     }
    //       }
    //     vessels[IdToKill]->SetRadius(0.000000000001_um);
    //       vessels[IdToKill]->SetToDie();
    //     */ /*
    //       for(it = vessels.begin(); it != vessels.end(); it++)
    //       {	
    //     if  (minRadius > (*it)->GetRadius())
    //     {
    //       minRadius = (*it)->GetRadius();
    //       IdToKill = (*it)->GetId();
    //     }
    //       }
    //     */

    //     // pick a segment to prune 
    //     //p_network->UniversalPruning(p_vessel_2,true);  

    //     // obtain the input node

    //     // are there other outgoing vessels at that node?

    //     // if no, make a list of the outgoing vessels, and pass each outgoing vessel to the pruning function

    //     // obtain the output node

    //     // are there other incoming vessels at that node?

    //     // if no, make a list of the incoming vessels, and pass each incoming vessel to the pruning function

    //     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //     // PRUNING TEST AREA ENDS
    //     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //     // Set up the grid
    //     std::shared_ptr<Part<3> > p_domain = Part<3>::Create();
    //     p_domain->AddCuboid(domain_x, domain_y, domain_z, Vertex<3>(0.0, 0.0, 0.0));
    //     std::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
    //     p_grid->GenerateFromPart(p_domain, 1.0*1_um);  // set domain and grid spacing

    //     // Set up the PDE expression for oxygen
    //     std::shared_ptr<DiscreteContinuumLinearEllipticPde<3> > p_oxygen_pde = DiscreteContinuumLinearEllipticPde<3>::Create();
    //     // QDiffusivity diffusivity(0.0033 * unit::metre_squared_per_second);     
    //     // QRate consumption_rate(2.e6 * unit::per_second);  
    //     // p_oxygen_pde->SetIsotropicDiffusionConstant(diffusivity);
    //     // p_oxygen_pde->SetContinuumLinearInUTerm(-consumption_rate);
    //     p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
    //     p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));
 
    //     // Set the concentration of O2 at the source 
    //     // QConcentration boundary_concentration(1.e-3 * unit::mole_per_metre_cubed);  // set the value directly
    //     QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") * GenericParameters::mpGasConcentrationAtStp->GetValue("User");
    //     QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp * Owen11Parameters::mpReferencePartialPressure->GetValue("User");

    //     // Set the boundary condition to be the network acting as a Dirichlet source
    //     std::shared_ptr<DiscreteContinuumBoundaryCondition<3> > p_vessel_boundary_condition = DiscreteContinuumBoundaryCondition<3>::Create();
    //     p_vessel_boundary_condition->SetValue(vessel_oxygen_concentration);  // set the O2 concentration at the source
    //     p_vessel_boundary_condition->SetType(BoundaryConditionType::VESSEL_LINE);
    //     p_vessel_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

    //     // Set up the solver
    //     SimpleLinearEllipticFiniteDifferenceSolver<3> solver;
    //     solver.SetGrid(p_grid);
    //     solver.SetPde(p_oxygen_pde);
    //     solver.AddBoundaryCondition(p_vessel_boundary_condition);
    //     solver.SetVesselNetwork(p_network);
    //     solver.SetLabel("oxygen");
    //     solver.SetFileName("oxygen_solution_0");

    //     // Set up output file management and seed the random number generator.
    //     auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestElementaryNetwork", true);        
    //     solver.SetFileHandler(p_output_file_handler);
        
    //     // Run the solver and write the output file (visualise with Paraview; set Filters->Alphabetical->Tube)
    //     solver.SetWriteSolution(true);
    //     solver.Solve();
    //     p_network->Write(p_output_file_handler->GetOutputDirectoryFullPath() + "equilateral_network_results.vtp");
 
    //     // Dump our parameter collection to an xml file and, importantly, clear it for the next test
    //     ParameterCollection::Instance()->DumpToFile(p_output_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
    //     ParameterCollection::Instance()->Destroy();
    //     BaseUnits::Instance()->Destroy();
    // }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Possibly Useful Code
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // // Set up loop for pruning
    	// std::vector<std::shared_ptr<Vessel<3> > > vessels = p_network->GetVessels();  // store all the vessels
        // // vessels = p_network->GetVessels();  // store all the vessels
        // std::vector<std::shared_ptr<Vessel<3> > >::iterator vessel_iterator;
        // int count=0;
        // for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)  // for all the vessels
        // { 
        //         p_network->RemoveVessel(*vessel_iterator,true);

        //         // Set up the solver
        //         SimpleLinearEllipticFiniteDifferenceSolver<3> solver;
        //         solver.SetGrid(p_grid);
        //         solver.SetPde(p_oxygen_pde);
        //         solver.AddBoundaryCondition(p_vessel_boundary_condition);
        //         solver.SetVesselNetwork(p_network);
        //         solver.SetLabel("oxygen");
        //         solver.SetFileName("oxygen_solution_0");

        //         // Set up output file management and seed the random number generator.
        //         auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestElementaryNetwork", true);        
        //         solver.SetFileHandler(p_output_file_handler);
                
        //         // Run the solver and write the output file (visualise with Paraview; set Filters->Alphabetical->Tube)
        //         solver.SetWriteSolution(true);
        //         solver.Solve();
        //         p_network->Write(p_output_file_handler->GetOutputDirectoryFullPath() + "equilateral_network_results.vtp");


        //         // p_network->Write(output_file_initial);
    	//         // std::string output_file = p_output_file_handler->GetOutputDirectoryFullPath().append("FinalNetwork"+to_string(alpha)+".vtp");  // saves the network after adding flow




        //         // Dump our parameter collection to an xml file and, importantly, clear it for the next test
        //         ParameterCollection::Instance()->DumpToFile(p_output_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
        //         ParameterCollection::Instance()->Destroy();
        //         BaseUnits::Instance()->Destroy();
                
        //         // Return the network to its original state
        //         p_network->AddVessel(*vessel_iterator);  // add the vessel that was pruned
        //         count++;  // increase the count
        // }

                // // Set nodes based on an equilateral network
        // std::vector<VesselNodePtr<2> > nodes;
        // nodes.push_back(VesselNode<2>::Create(0.0_um, mid_domain_y));
        // nodes.push_back(VesselNode<2>::Create(vessel_length, mid_domain_y));
        // nodes.push_back(VesselNode<2>::Create(vessel_height+vessel_length, vessel_height+mid_domain_y));
        // nodes.push_back(VesselNode<2>::Create(vessel_height+vessel_length, -vessel_height+mid_domain_y));
        // nodes.push_back(VesselNode<2>::Create(vessel_height+vessel_height+vessel_length, mid_domain_y));
        // nodes.push_back(VesselNode<2>::Create(vessel_height+vessel_height+vessel_length+vessel_length, mid_domain_y));

        // // Make vessels and vessel segments 
        // std::vector<std::shared_ptr<Vessel<2> > > vessels;
        // vessels.push_back(Vessel<2>::Create(VesselSegment<2>::Create(nodes[0], nodes[1])));
        // vessels.push_back(Vessel<2>::Create(VesselSegment<2>::Create(nodes[1], nodes[2])));
        // vessels.push_back(Vessel<2>::Create(VesselSegment<2>::Create(nodes[1], nodes[3])));
        // vessels.push_back(Vessel<2>::Create(VesselSegment<2>::Create(nodes[2], nodes[4])));
        // vessels.push_back(Vessel<2>::Create(VesselSegment<2>::Create(nodes[3], nodes[4])));
        // vessels.push_back(Vessel<2>::Create(VesselSegment<2>::Create(nodes[4], nodes[5])));


        // p_network->MergeCoincidentNodes();
        // p_network->UpdateAll();
