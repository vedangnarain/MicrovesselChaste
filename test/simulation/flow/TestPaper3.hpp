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

This script contains tests that are used to simulate blood flow, oxygenation, and radiotherapy in different vascular architectures. 

The default location for outputs is /tmp/<home_directory>/testoutput. 

H = Haematocrit
BC = boundary condition
RT = Radiotherapy
PQ = Perfusion Quotient
Fig. = Figure
# = number

TO DO LIST
 
Currently, we manually set no-flow vessels to have zero H, but maybe we should either put that bit of code in the solvers. Or set all relevant parameters (like viscosity) to zero too.

The y-extent for a non-hierarchical network is currently set manually based on the simulations. Find a way to make it appropriately cover the hexagonal network automatically.

22/8/24

*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialisation
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef TESTPAPER3_HPP_
#define TESTPAPER3_HPP_
#define _BACKWARD_BACKWARD_WARNING_H 1  // Cut out the VTK deprecated warning

// Essential functionality
#include <boost/lexical_cast.hpp>
#include <cxxtest/TestSuite.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>  // needed for exp function
#include <sstream>
#include <string>
// #include <time.h>
// #include <vector>
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
// #include "Part.hpp"

// Dimensional analysis
#include "BaseUnits.hpp"
#include "UnitCollection.hpp"
#include "GenericParameters.hpp"
#include "ParameterCollection.hpp"
#include "Owen11Parameters.hpp"
#include "Secomb04Parameters.hpp"
// #include "Vertex.hpp"

// Grids and PDEs
#include "CellBasedDiscreteSource.hpp"
#include "CellStateDependentDiscreteSource.hpp"
#include "SimpleLinearEllipticFiniteDifferenceSolver.hpp"
// #include "SimpleLinearEllipticFiniteElementSolver.hpp"
// #include "DiscreteContinuumBoundaryCondition.hpp"
#include "DiscreteContinuumLinearEllipticPde.hpp"
#include "RegularGrid.hpp"
#include "GridCalculator.hpp"
// #include "CellwiseSourceEllipticPde.hpp"
// #include "ConstBoundaryCondition.hpp"
// #include "EllipticGrowingDomainPdeModifier.hpp"
#include "VesselBasedDiscreteSource.hpp"

// Vessel networks
#include "Vessel.hpp"
#include "VesselNetwork.hpp"
// #include "VesselNetworkReader.hpp"
#include "VesselNetworkGenerator.hpp"
#include "VesselNetworkGeometryCalculator.hpp"
#include "VesselNetworkPropertyManager.hpp"
#include "VesselNode.hpp"
#include "VesselSegment.hpp"

// Flow
#include "NodeFlowProperties.hpp"
#include "VesselImpedanceCalculator.hpp"
#include "FlowSolver.hpp"
// #include "WallShearStressCalculator.hpp"
// #include "MechanicalStimulusCalculator.hpp"
// #include "MetabolicStimulusCalculator.hpp"
// #include "ShrinkingStimulusCalculator.hpp"
// #include "StructuralAdaptationSolver.hpp"
#include "ViscosityCalculator.hpp"

// Haematocrit
// #include "AlarconHaematocritSolver.hpp"
#include "BetteridgeHaematocritSolver.hpp"
#include "ConstantHaematocritSolver.hpp"
// #include "LinnengarHaematocritSolver.hpp"
// #include "YangHaematocritSolver.hpp"
#include "PriesHaematocritSolver.hpp"
#include "PriesWithMemoryHaematocritSolver.hpp"

// Cells
#include "CellsGenerator.hpp"
// #include "CellBasedSimulationArchiver.hpp"
#include "CellLabel.hpp"
#include "Owen2011OxygenBasedCellCycleModel.hpp"
// #include "UniformCellCycleModel.hpp"
// #include "HoneycombMeshGenerator.hpp"
// #include "CellProliferativeTypesCountWriter.hpp"
// #include "CellPopulationAreaWriter.hpp"
// #include "CellPopulationElementWriter.hpp"
// #include "RandomCellKiller.hpp"
#include "LQRadiotherapyCellKiller.hpp"
#include "CancerCellMutationState.hpp"
// #include "StalkCellMutationState.hpp"
#include "QuiescentCancerCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "Owen11CellPopulationGenerator.hpp"
#include "Owen2011TrackingModifier.hpp"
#include "CaBasedCellPopulation.hpp"
#include "ApoptoticCellKiller.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
// #include "TransitCellProliferativeType.hpp"
#include "NoCellCycleModel.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"

// Forces
#include "GeneralisedLinearSpringForce.hpp"

// Angiogenesis
#include "Owen2011SproutingRule.hpp"
#include "Owen2011MigrationRule.hpp"
#include "AngiogenesisSolver.hpp"

// Vessel regression solver
// #include "WallShearStressBasedRegressionSolver.hpp"

// General solver to collect all the flows
#include "MicrovesselSolver.hpp"
#include "MicrovesselSimulationModifier.hpp"
#include "OnLatticeSimulation.hpp"
#include "OffLatticeSimulation.hpp"

// Visualisation
// #include "MicrovesselVtkScene.hpp"
// #include "VtkSceneMicrovesselModifier.hpp"

// Keep this last
#include "PetscAndVtkSetupAndFinalize.hpp"

// Includes standard library
using namespace std;

// Tests with cells
class TestPaper3 : public AbstractCellBasedWithTimingsTestSuite
{

public:

    // Simulate a 2D hexagonal network on a PDE grid with flow and H-splitting (single inlet and outlet) and cells
    void xTestSingleFeedHexagonalUnit2D()
    {
        // Set file name
        std::ostringstream strs;
        strs << std::fixed << std::setprecision( 2 );
        strs << "TestDebuggingNetworks/CellularSingleFeedHexagonalUnit/";        

        // Seed the random number generator
        RandomNumberGenerator::Instance()->Reseed(12345);

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);
        // BaseUnits::Instance()->SetReferenceTimeScale(1_h);
        QTime reference_time = 3600.0*unit::seconds;
        BaseUnits::Instance()->SetReferenceTimeScale(reference_time);

        // Match this based on cell size
        QLength cell_grid_spacing = 10.0_um;  // the simulation time gets quite long if you reduce the resolution further
    
        // Define the key parameters
        QLength homogeneous_vessel_radius = 10_um;
        double dimless_domain_size_x = 400.0;  // x-coordinate of output node
        double dimless_domain_size_y = 346.41 + 173.205;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        unsigned dimless_vessel_length = 100.0;
        // std::vector<std::vector<unsigned> > Order;
        std::shared_ptr<VesselNetwork<2> > p_network;
        std::vector<std::shared_ptr<Vessel<2> > > vessels;
        string line2;

        // Initialise the simulation space
        QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
        QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
        std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
        std::ofstream outfileMean;

        // Set up the grid for the finite difference solver
        auto p_grid = RegularGrid<2>::Create();
        p_grid->SetSpacing(cell_grid_spacing);
        c_vector<unsigned, 3> dimensions;
        dimensions[0] = unsigned((domain_side_length_x)/(cell_grid_spacing))+1; // num x
        dimensions[1] = unsigned((domain_side_length_y)/(cell_grid_spacing))+1; // num_y
        dimensions[2] = 1;
        p_grid->SetDimensions(dimensions);

        // Define the key blood flow parameters
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Read the network from a file
        VesselNetworkGenerator<2> network_generator;
        std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/edges/honeycomb/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_SingleFeed.txt");
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

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {
            // std::ostringstream strs;
            if (h_solver==1)
            {
                strs << "ConstantHaematocrit";
            }
            else if (h_solver==2)
            {
                strs << "PriesHaematocrit";
            }
            else if (h_solver==3)
            {
                strs << "MemoryHaematocrit";
            }
            else if (h_solver==4)
            {
                strs << "FungHaematocrit";
            }
            std::string str_directory_name = strs.str();
            auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);            

            // Generate the network
            p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);
            
            // Set inlet and outlet nodes
            auto p_segment = p_network->GetVesselSegments()[0];
            vessels = p_network->GetVessels();
            for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
            {
                if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
                {
                    if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)
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
                {
                    if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)
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
            p_segment->SetRadius(homogeneous_vessel_radius);
            VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);

            // Save just the network architecture without flow
            // p_network->Write(p_file_handler->GetOutputDirectoryFullPath() + "network_structure.vtp");

            // Save just the grid
            // p_grid->Write(p_file_handler);

            // Set up the cell populations
            std::shared_ptr<GridCalculator<2> > p_grid_calc = GridCalculator<2>::Create();
            p_grid_calc->SetGrid(p_grid);
         
            std::shared_ptr<Owen11CellPopulationGenerator<2> > p_cell_population_genenerator = Owen11CellPopulationGenerator<2>::Create();
            p_cell_population_genenerator->SetGridCalculator(p_grid_calc);
            p_cell_population_genenerator->SetVesselNetwork(p_network);
            QLength tumour_radius(50.0 * unit::microns);
            p_cell_population_genenerator->SetTumourRadius(tumour_radius);
            std::shared_ptr<CaBasedCellPopulation<2> > p_cell_population = p_cell_population_genenerator->Update();

            // Set up the PDE
            auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
            p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
            // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));
            auto p_cell_oxygen_sink = CellBasedDiscreteSource<2>::Create();  // Set the cell terms
            // p_cell_oxygen_sink->SetLinearInUConsumptionRatePerCell(13.0);
            p_cell_oxygen_sink->SetLinearInUConsumptionRatePerCell(Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));
            p_oxygen_pde->AddDiscreteSource(p_cell_oxygen_sink);

            // Set up the discrete vessel source
            auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
            QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    GenericParameters::mpGasConcentrationAtStp->GetValue("User");
            QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    Owen11Parameters::mpReferencePartialPressure->GetValue("User");
            p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
            p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
            p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
            p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

            // Set up the finite difference solver for oxyge
            auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
            p_oxygen_solver->SetPde(p_oxygen_pde);
            p_oxygen_solver->SetLabel("oxygen");
            p_oxygen_solver->SetGrid(p_grid);
            // solver.SetFileName("oxygen_solution_0");
            // p_oxygen_solver->SetWriteSolution(true);

            // Set up the VEGF PDE as the O2 one
            auto p_vegf_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
            p_vegf_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpVegfDiffusivity->GetValue("User"));
            // p_vegf_pde->SetIsotropicDiffusionConstant(0.0);
            p_vegf_pde->SetContinuumLinearInUTerm(-Owen11Parameters::mpVegfDecayRate->GetValue("User"));
        
            // Set up a map for different release rates depending on cell type. Also include a threshold intracellular VEGF below which there is no release.
            auto p_normal_and_quiescent_cell_source = CellStateDependentDiscreteSource<2>::Create();  // initialise the source (which combines normal and quiescent cells)
            std::map<unsigned, QConcentrationFlowRate > normal_and_quiescent_cell_rates;  // map to store secretion rates
            std::map<unsigned, QConcentration > normal_and_quiescent_cell_rate_thresholds;  // map to store secretion thresholds
            MAKE_PTR(QuiescentCancerCellMutationState, p_quiescent_cancer_state);
            MAKE_PTR(WildTypeCellMutationState, p_normal_cell_state);
            normal_and_quiescent_cell_rates[p_normal_cell_state->GetColour()] = Owen11Parameters::mpCellVegfSecretionRate->GetValue("User");  // set the rate for normal cells
            normal_and_quiescent_cell_rate_thresholds[p_normal_cell_state->GetColour()] = 0.27*unit::mole_per_metre_cubed;  // set the threshold for normal cells
            // normal_and_quiescent_cell_rates[p_normal_cell_state->GetColour()] = 0.00000000000000000001;  // set the rate for normal cells
            // normal_and_quiescent_cell_rate_thresholds[p_normal_cell_state->GetColour()] = 1000000_M;  // set the threshold for normal cells
            normal_and_quiescent_cell_rates[p_quiescent_cancer_state->GetColour()] = Owen11Parameters::mpCellVegfSecretionRate->GetValue("User");  // set the rate for quiescent cells
            normal_and_quiescent_cell_rate_thresholds[p_quiescent_cancer_state->GetColour()] = 0_M;  // set the threshold for quiescent cells
            // normal_and_quiescent_cell_rates[p_quiescent_cancer_state->GetColour()] = 0.00000000000000000001;  // set the rate for quiescent cells
            // normal_and_quiescent_cell_rate_thresholds[p_quiescent_cancer_state->GetColour()] = 1000000_M;  // set the threshold for quiescent cells
            p_normal_and_quiescent_cell_source->SetStateRateMap(normal_and_quiescent_cell_rates);  // configure the map
            p_normal_and_quiescent_cell_source->SetLabelName("VEGF");  // set the label for the source
            p_normal_and_quiescent_cell_source->SetStateRateThresholdMap(normal_and_quiescent_cell_rate_thresholds);
            // p_vegf_pde->AddDiscreteSource(p_normal_and_quiescent_cell_source);

            // Add a vessel related VEGF sink
            auto p_vessel_vegf_sink = VesselBasedDiscreteSource<2>::Create();
            p_vessel_vegf_sink->SetReferenceConcentration(0.0_M);
            p_vessel_vegf_sink->SetVesselPermeability(Owen11Parameters::mpVesselVegfPermeability->GetValue("User"));
            // p_vessel_vegf_sink->SetVesselPermeability(0.0);
            // p_vegf_pde->AddDiscreteSource(p_vessel_vegf_sink);

            // Set up a finite difference solver as before for the VEGF
            auto p_vegf_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
            p_vegf_solver->SetPde(p_vegf_pde);
            p_vegf_solver->SetLabel("VEGF_Extracellular");  // set the label for the external field
            p_vegf_solver->SetGrid(p_grid);
            
            // Set the haematocrit solver
            std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
            if (h_solver==1)
            {
                std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
                p_abstract_haematocrit_solver = p_haematocrit_solver;       
            }
            else if (h_solver==2)
            {
                std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
                p_abstract_haematocrit_solver = p_haematocrit_solver;                    
            }
            else if (h_solver==3)
            {
                std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
                p_abstract_haematocrit_solver = p_haematocrit_solver;      
            }
            else if (h_solver==4)
            {
                std::cout << "Now using FungHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
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

            // Set up flag for broken solver
            unsigned broken_solver = 0;

            // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
            unsigned max_iter = 1000;  // 1000 
            double tolerance2 = 1.e-10;
            std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
            std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
            for(unsigned idx=0;idx<max_iter;idx++)
            {
                // Run the solvers
                p_viscosity_calculator->Calculate();
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
                    broken_solver = 0;
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
                    std::cout << "Problem encountered in " << str_directory_name << std::endl;
                    error_log << "\n Problem encountered in " << str_directory_name << std::endl; 
                    broken_solver = 1;
                    break;
                }
            }

            // If solver doesn't converge, move on to next one
            if (broken_solver == 1)
            {
                continue;
            }

            // Set up the angiogenesis solver (runs on deactivated VEGF function but seems to be necessary for main solver)
            std::shared_ptr<AngiogenesisSolver<2> > p_angiogenesis_solver = AngiogenesisSolver<2>::Create();
            std::shared_ptr<Owen2011SproutingRule<2> > p_sprouting_rule = Owen2011SproutingRule<2>::Create();
            std::shared_ptr<Owen2011MigrationRule<2> > p_migration_rule = Owen2011MigrationRule<2>::Create();
            p_angiogenesis_solver->SetMigrationRule(p_migration_rule);
            p_angiogenesis_solver->SetSproutingRule(p_sprouting_rule);
            p_sprouting_rule->SetDiscreteContinuumSolver(p_vegf_solver);
            p_migration_rule->SetDiscreteContinuumSolver(p_vegf_solver);
            p_angiogenesis_solver->SetVesselGridCalculator(p_grid_calc);
            p_angiogenesis_solver->SetVesselNetwork(p_network);
            
            // Initialise the solver
            auto p_microvessel_solver = MicrovesselSolver<2>::Create();
            p_microvessel_solver->SetVesselNetwork(p_network);
            p_microvessel_solver->SetOutputFrequency(5_h);
            // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
            p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
            p_microvessel_solver->AddDiscreteContinuumSolver(p_vegf_solver);
            p_microvessel_solver->SetAngiogenesisSolver(p_angiogenesis_solver);

            // Specifies which extracellular fields to update based on PDE
            boost::shared_ptr<MicrovesselSimulationModifier<2> > p_microvessel_modifier =
            boost::shared_ptr<MicrovesselSimulationModifier<2> >(new MicrovesselSimulationModifier<2> ());
            p_microvessel_modifier->SetMicrovesselSolver(p_microvessel_solver);
            std::vector<std::string> update_labels;
            update_labels.push_back("oxygen");
            p_microvessel_modifier->SetCellDataUpdateLabels(update_labels);

            // The full simulation is run as a typical Cell Based Chaste simulation
            OnLatticeSimulation<2> simulator(*p_cell_population);
            simulator.AddSimulationModifier(p_microvessel_modifier);

            // Add a killer to remove apoptotic cells (leave it commented so RT is the only killer)
            // boost::shared_ptr<ApoptoticCellKiller<2> > p_apoptotic_cell_killer(new ApoptoticCellKiller<2>(p_cell_population.get())); (leave it commented so RT is the only killer)
            // simulator.AddCellKiller(p_apoptotic_cell_killer);
 
            // Add a LQ RT killer
            boost::shared_ptr<LQRadiotherapyCellKiller<2> > p_rt_killer =
            boost::shared_ptr<LQRadiotherapyCellKiller<2> >(new LQRadiotherapyCellKiller<2> (p_cell_population.get()));
            p_rt_killer->SetAlphaMax(0.3*unit::per_gray);
            p_rt_killer->SetBetaMax(0.03*unit::per_gray_squared);
            p_rt_killer->SetDoseInjected(2.0*unit::gray);
            p_rt_killer->SetCancerousRadiosensitivity(0.3 * unit::per_gray, 0.03 * unit::per_gray_squared);
            p_rt_killer->SetNormalRadiosensitivity(0.15 * unit::per_gray, 0.05 * unit::per_gray_squared);
            p_rt_killer->SetOerAlphaMax(1.75);
            p_rt_killer->SetOerAlphaMin(1.0);
            p_rt_killer->SetOerBetaMax(3.25);
            p_rt_killer->SetOerBetaMin(1.0);
            // p_rt_killer->SetOerConstant(oxygen_solubility_at_stp * 3.28_Pa);
            p_rt_killer->SetOerConstant(0.004499758549541243*unit::mole_per_metre_cubed);  // K_OER (3.28 mmHg)
            p_rt_killer->SetAlphaMax(0.3 * unit::per_gray);
            p_rt_killer->SetBetaMax(0.03 * unit::per_gray_squared);
            p_rt_killer->UseOer(true);

            // Set Up the dosage times
            // std::vector<QTime > rt_times;
            // rt_times.push_back(3600.0*24.0*unit::seconds);
            // rt_times.push_back(3600.0*48.0*unit::seconds);
            // rt_times.push_back(3600.0*72.0*unit::seconds);
            // p_rt_killer->SetTimeOfRadiation(rt_times);
            // p_rt_killer->AddTimeOfRadiation(3600.0*96.0*unit::seconds);
            // simulator.AddCellKiller(p_rt_killer);

            // If needed
            // boost::shared_ptr<RandomCellKiller<2> > p_cell_killer(new RandomCellKiller<2>(p_cell_population.get(), 0.01));
            // simulator.AddCellKiller(p_cell_killer);
            // for (typename AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population.get()->Begin(); cell_iter != p_cell_population.get()->End(); ++cell_iter)
            // {
            //     // CheckAndLabelSingleCellForApoptosis(*cell_iter);
            //     std::cout << "hi cell" << (*cell_iter)->GetCellData()->GetItem("oxygen") << std::endl;

            //     //  std::cout <<  BaseUnits::Instance()->GetReferenceConcentrationScale() << std::endl;


            //     //     double oxygen_concentration;

            //     // std::cout << "hi cell" << oxygen_concentration << std::endl;
            // }

            // Add another modifier for updating cell cycle quantities.
            boost::shared_ptr<Owen2011TrackingModifier<2> > p_owen11_tracking_modifier(new Owen2011TrackingModifier<2>);
            simulator.AddSimulationModifier(p_owen11_tracking_modifier);

            // Print the average oxygenation
            // std::vector<double> solution = p_oxygen_solver->GetSolution();
            // double average_oxygen = 0.0;
            // for(unsigned jdx=0;jdx<solution.size();jdx++)
            // {
            //     average_oxygen += solution[jdx];
            // }
            // average_oxygen /= double(solution.size());
            // std::cout << "Average oxygen: " << average_oxygen << std::endl;

            // Set up the simulation time and run it
            simulator.SetOutputDirectory(str_directory_name);
            // simulator.SetSamplingTimestepMultiple(10);  // get sample every x steps (multiple*dt = time between s)
            simulator.SetSamplingTimestepMultiple(1);
            // simulator.SetDt(1);  // simulation time step in hours
            simulator.SetDt(1.0);
            simulator.SetEndTime(24*1.0);  // end time in hours
            // simulator.SetEndTime(24);  // end time in hours
            // simulator.SetEndTime(2.0);
            simulator.Solve();

            // Write the output network file (visualise with Paraview: set Filters->Alphabetical->Tube)
            std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
            p_network->Write(output_file);

            // Dump our parameter collection to an xml file and, importantly, clear it for the next test
            ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
            ParameterCollection::Instance()->Destroy();
            BaseUnits::Instance()->Destroy();
            SimulationTime::Instance()->Destroy();
        }
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Simulate a 2D Voronoi network on a PDE grid with flow and H-splitting and O2 and cells
    void TestVoronoiNetwork2DPaper3() 
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Set network in filename
        std::string network_name = "TestVoronoiNetworkTissue/";

        // Choose the vascular parameters
        unsigned NumberOfSeedPoints = 60;  // change this to select which Voronoi architecture to use: 25, 100, 400
        unsigned NumberOfLayouts = 1000;  // number of different point layouts to run simulations with (add an extra selection for a demo)
        double dimless_domain_size_x = 1105.0; 
        double dimless_domain_size_y = 1105.0 + (2.0*85.0);  // same as Mantegazza network
        // double dimless_domain_size_x = 2050.0;  // x-coordinate of output node
        // double dimless_domain_size_y = 1818.65 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        // QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        QDynamicViscosity viscosity = 1.96*1.e-3*unit::poiseuille;
        // double initial_haematocrit = 0.45;        
        // double initial_haematocrit = 0.1;  // conks out somewhere between 0.35 and 0.4
        double initial_haematocrit = 0.3;  // conks out somewhere between 0.35 and 0.4
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        // unsigned n_vessels = 382;  // number of non-inlet/outlet vessels from which to select ones to make thin
        // double percToKill = 0.2;  // percentage of vessels to kill
        // unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
        // unsigned ToBeKilled = 0.5*((3*NumberOfSeedPoints)-6);  // number to kill
        // unsigned ToBeKilled = NumberOfSeedPoints;  // number to kill
        // unsigned ToBeKilled = 100;  // number to kill
        unsigned MaxKills = (3*NumberOfSeedPoints)-6;  // number to kill

        // Match this based on cell size
        QLength cell_grid_spacing = 20.0_um;

        // Seed the random number generator
        RandomNumberGenerator::Instance()->Reseed(12345);

        // Create the log file for broken simulations
        std::ofstream broken_layouts_file;
        broken_layouts_file.open("/scratch/narain/testoutput/TestVoronoiNetwork/broken_layouts.txt");
        broken_layouts_file.close();        

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);
        // BaseUnits::Instance()->SetReferenceTimeScale(1_h);
        QTime reference_time = 3600.0*unit::seconds;  // 1 hour
        BaseUnits::Instance()->SetReferenceTimeScale(reference_time);

        ////////////////////////////////////////////////////////////////
        // Vasculature
        ////////////////////////////////////////////////////////////////

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=3; h_solver<=3; h_solver++)
        {     
            // Initialise the simulation
            std::shared_ptr<VesselNetwork<2> > p_network;
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            VesselPtr<2> p_current_vessel;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            // std::vector<std::vector<unsigned> > Order;

            // Run simulations with different layouts (could be used to compute some average later)
            // std::vector<unsigned> broken_selections = {2, 9, 11, 49, 69, 80, 93, 99};
            // std::vector<unsigned> broken_selections = {99};
            unsigned seed_number = 42;  // which seed to use for the random number generator for the edge matrix generation 
            // for (unsigned list_number : broken_selections)
            for (unsigned layout=4;layout<=NumberOfLayouts; layout++)   
            { 
                // Set up flag for broken solver
                unsigned broken_solver = 0;
                
                // Read the network layout from a file
                VesselNetworkGenerator<2> network_generator;
                std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/edges/voronoi/"+to_string(NumberOfSeedPoints)+"SeedPoints/random_seed_" + to_string(seed_number) + "/verified/EdgesMatrixSampleNumber"+to_string(layout)+".txt");
                std::vector<std::vector<double> > rEdgesMatrix;
                string line;
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

                // Set inlet and outlet nodes (if node (start/end) is at the edge of the domain space, make it input or output depending on which side), and assign each vessel a radius from the list  
                auto p_segment = p_network->GetVesselSegments()[0];
                vessels = p_network->GetVessels();
                
                // Assign the first vessel as the input and the last as the output and remove all other dead ends (i.e., make the network single inlet and outlet)
                unsigned counter = 0;
                for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                {                                         
                    // std::cout << counter << std::endl;
                    if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
                    {
                        if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] <  tolerance)
                        {
                            if (counter==0)
                            {           
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(5530.0*unit::pascals);
                            }
                            else if (counter != vessels.size() - 1)
                            {
                                p_network->RemoveVessel(*vessel_iterator, true);
                                counter++;
                                // std::cout << "DELETE1" << std::endl;
                                continue;
                            }
                        }
                        //if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                        else  // any dead-end not at x=0 is outlet (the geometry generation process removes vessels at top and bottom so only rightmost remain; nodes are designated at start so pruning makes no difference)
                        {
                            if (counter==vessels.size()-1)
                            {    
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                            }
                            // (*vessel_iterator)->SetRadius(inlet_vessel_radius);
                            else if (counter != vessels.size() - 1)
                            {
                                p_network->RemoveVessel(*vessel_iterator, true);                                        
                                counter++;
                                // std::cout << "DELETE2" << std::endl;
                                    continue;
                            }
                        }
                    }
                    if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
                    {
                        if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] <  tolerance)
                        {
                            if (counter==0)
                            {    
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(5530.0*unit::pascals);
                            }
                            else if (counter != vessels.size() - 1)
                            {
                                p_network->RemoveVessel(*vessel_iterator, true);                                        
                                counter++;
                                // std::cout << "DELETE3" << std::endl;
                                continue;

                            }
                        }
                        //if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                        else
                        {
                            if (counter==vessels.size()-1)
                            {   
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                            }
                            // (*vessel_iterator)->SetRadius(inlet_vessel_radius);
                            else if (counter != vessels.size() - 1)
                            {
                                p_network->RemoveVessel(*vessel_iterator, true);                                        
                                counter++;
                                // std::cout << "DELETE4" << std::endl;
                                continue;
                            }
                        }
                    }
                    counter++;
                }
                vessels = p_network->GetVessels();

                // Remove diameter heterogeneity
                p_segment = p_network->GetVesselSegments()[0];
                p_segment->SetRadius(7.5_um);  // same as Mantegazza
                VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
                
                // Prune all vessels up to specified dose 
                for(unsigned KilledVessels=0; KilledVessels < MaxKills; KilledVessels++)
                { 
                    // Display status message
                    std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    std::string solver_name;
                    if (h_solver==1)
                    {
                        solver_name = "ConstantHaematocrit"; 
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==2)
                    {
                        solver_name = "PriesHaematocrit"; 
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                
                    }
                    else if (h_solver==3)
                    {
                        solver_name = "MemoryHaematocrit"; 
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==4)
                    {
                        solver_name = "FungHaematocrit"; 
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
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

                    // // Prune all vessels up to specified dose 
                    // for(unsigned KilledVessels=0; KilledVessels <= ToBeKilled; KilledVessels++)
                    // { 
                        // // Display status message
                        // std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                    // Set filename
                    std::stringstream selection_stream;
                    std::stringstream kill_stream;
                    selection_stream << std::fixed << std::setprecision(0) << to_string(layout);                    
                    kill_stream << std::fixed << std::setprecision(0) << KilledVessels;
                    std::string selection_string = selection_stream.str();
                    std::string kill_string = kill_stream.str();                    
                    std::string file_name = network_name + solver_name + "/Selection" + selection_string + "/Kills" + kill_string;                        
                    std::string str_directory_name = file_name;
                    std::cout << str_directory_name << std::endl;
                    auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                    // Set up the grid for the finite difference solver
                    auto p_grid = RegularGrid<2>::Create();
                    p_grid->SetSpacing(cell_grid_spacing);
                    c_vector<unsigned, 3> dimensions;
                    dimensions[0] = unsigned((domain_side_length_x)/(cell_grid_spacing))+2; // num x
                    dimensions[1] = unsigned((domain_side_length_y)/(cell_grid_spacing))+2; // num_y
                    dimensions[2] = 1;
                    p_grid->SetDimensions(dimensions);

                    // Set up the cell populations
                    std::shared_ptr<GridCalculator<2> > p_grid_calc = GridCalculator<2>::Create();
                    p_grid_calc->SetGrid(p_grid);
                    std::shared_ptr<Owen11CellPopulationGenerator<2> > p_cell_population_genenerator = Owen11CellPopulationGenerator<2>::Create();
                    p_cell_population_genenerator->SetGridCalculator(p_grid_calc);
                    p_cell_population_genenerator->SetVesselNetwork(p_network);
                    QLength tumour_radius(10000.0 * unit::microns);  // just make it bigger than the domain to ensure all cells are tumour cells 
                    p_cell_population_genenerator->SetTumourRadius(tumour_radius);
                    std::shared_ptr<CaBasedCellPopulation<2> > p_cell_population = p_cell_population_genenerator->Update();

                    // Choose the PDE
                    auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                                
                    // Set the diffusivity and decay terms
                    p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    auto p_cell_oxygen_sink = CellBasedDiscreteSource<2>::Create();  // Set the cell terms
                    p_cell_oxygen_sink->SetLinearInUConsumptionRatePerCell(Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));
                    p_oxygen_pde->AddDiscreteSource(p_cell_oxygen_sink);

                    // Set up the vessel source
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
                    // solver.AddBoundaryCondition(p_vessel_boundary_condition);
                    // solver.SetVesselNetwork(p_network);

                    // Set up the VEGF PDE as we set the O2 one (currently turned off)
                    auto p_vegf_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    p_vegf_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpVegfDiffusivity->GetValue("User"));
                    // p_vegf_pde->SetIsotropicDiffusionConstant(0.0);
                    p_vegf_pde->SetContinuumLinearInUTerm(-Owen11Parameters::mpVegfDecayRate->GetValue("User"));
                
                    // Set up a map for different release rates depending on cell type. Also include a threshold intracellular VEGF below which there is no release (currently turned off)
                    auto p_normal_and_quiescent_cell_source = CellStateDependentDiscreteSource<2>::Create();  // initialise the source (which combines normal and quiescent cells)
                    std::map<unsigned, QConcentrationFlowRate > normal_and_quiescent_cell_rates;  // map to store secretion rates
                    std::map<unsigned, QConcentration > normal_and_quiescent_cell_rate_thresholds;  // map to store secretion thresholds
                    MAKE_PTR(QuiescentCancerCellMutationState, p_quiescent_cancer_state);
                    MAKE_PTR(WildTypeCellMutationState, p_normal_cell_state);
                    normal_and_quiescent_cell_rates[p_normal_cell_state->GetColour()] = Owen11Parameters::mpCellVegfSecretionRate->GetValue("User");  // set the rate for normal cells
                    normal_and_quiescent_cell_rate_thresholds[p_normal_cell_state->GetColour()] = 0.27*unit::mole_per_metre_cubed;  // set the threshold for normal cells
                    // normal_and_quiescent_cell_rates[p_normal_cell_state->GetColour()] = 0.00000000000000000001;  // set the rate for normal cells
                    // normal_and_quiescent_cell_rate_thresholds[p_normal_cell_state->GetColour()] = 1000000_M;  // set the threshold for normal cells
                    normal_and_quiescent_cell_rates[p_quiescent_cancer_state->GetColour()] = Owen11Parameters::mpCellVegfSecretionRate->GetValue("User");  // set the rate for quiescent cells
                    normal_and_quiescent_cell_rate_thresholds[p_quiescent_cancer_state->GetColour()] = 0_M;  // set the threshold for quiescent cells
                    // normal_and_quiescent_cell_rates[p_quiescent_cancer_state->GetColour()] = 0.00000000000000000001;  // set the rate for quiescent cells
                    // normal_and_quiescent_cell_rate_thresholds[p_quiescent_cancer_state->GetColour()] = 1000000_M;  // set the threshold for quiescent cells
                    p_normal_and_quiescent_cell_source->SetStateRateMap(normal_and_quiescent_cell_rates);  // configure the map
                    p_normal_and_quiescent_cell_source->SetLabelName("VEGF");  // set the label for the source
                    p_normal_and_quiescent_cell_source->SetStateRateThresholdMap(normal_and_quiescent_cell_rate_thresholds);
                    // p_vegf_pde->AddDiscreteSource(p_normal_and_quiescent_cell_source);

                    // Print the cell types and their colours
                    MAKE_PTR(CancerCellMutationState, p_cancer_state);
                    // MAKE_PTR(StemCellProliferativeType, p_stem_state);
                    std::cout << "Normal state has Legacy Cell type " <<  p_normal_cell_state->GetColour() << std::endl;
                    std::cout << "Cancer state has Legacy Cell type " <<  p_cancer_state->GetColour() << std::endl;
                    std::cout << "Quiescent cancer state has Legacy Cell type " <<  p_quiescent_cancer_state->GetColour() << std::endl;
                    // std::cout << "Quiescent cancer state has Legacy Cell type " <<  p_stem_state->GetColour() << std::endl;

                    // Add a vessel related VEGF sink (currently turned off)
                    auto p_vessel_vegf_sink = VesselBasedDiscreteSource<2>::Create();
                    p_vessel_vegf_sink->SetReferenceConcentration(0.0_M);
                    p_vessel_vegf_sink->SetVesselPermeability(Owen11Parameters::mpVesselVegfPermeability->GetValue("User"));
                    // p_vessel_vegf_sink->SetVesselPermeability(0.0);
                    // p_vegf_pde->AddDiscreteSource(p_vessel_vegf_sink);

                    // Set up a finite difference solver as before for the VEGF
                    auto p_vegf_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    p_vegf_solver->SetPde(p_vegf_pde);
                    p_vegf_solver->SetLabel("VEGF_Extracellular");  // set the label for the external field
                    p_vegf_solver->SetGrid(p_grid);
                    
                    // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                    unsigned max_max_iter = 10000; 
                    unsigned max_iter = 1000;  // 1000 
                    double prev_max_difference = 0.0;
                    double tolerance2 = 1.e-10;
                    std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                    std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                    for(unsigned idx=0;idx<max_max_iter;idx++)
                    {
                        // Run the solvers (order of calculators matters!)
                        p_viscosity_calculator->Calculate();
                        p_impedance_calculator->Calculate();
                        flow_solver.SetUp();
                        flow_solver.Solve();
                        p_abstract_haematocrit_solver->Calculate();

                        // Get the residual
                        double max_difference = 0.0;
                        double h_for_max = 0.0;
                        double prev_for_max = 0.0;
                        for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                        {
                            // Set segments with no flow to be dead (this is only problem for inlet really)
                            if (fabs(segments[jdx]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                            {
                                segments[jdx]->GetFlowProperties()->SetViscosity(0.0);
                                segments[jdx]->GetFlowProperties()->SetImpedance(0.0);
                                segments[jdx]->GetFlowProperties()->SetHaematocrit(0.0);
                            }    
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
                            // broken_solver = 0;
                            break;
                        }
                        else
                        {
                            std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                            if(idx%100==0)
                            {
                                // std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                                // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                                // p_network->Write(output_file);
                            }
                        }

                        // If there is no convergence after all the iterations, print the error message.
                        if(idx>=max_iter-1 && max_difference >= prev_max_difference)
                        {
                            std::cout << "Problem encountered in " << str_directory_name << std::endl;
                            error_log << "\n Problem encountered in " << str_directory_name << std::endl; 
                            broken_solver = 1;
                            break;
                        }                                
                        prev_max_difference = max_difference;
                    }

                    // If simulation doesn't converge, move on to next layout and log problem 
                    if (broken_solver == 1)
                    {
                        broken_layouts_file.open("/scratch/narain/testoutput/TestVoronoiNetwork/broken_layouts.txt", std::ios_base::app);
                        broken_layouts_file << str_directory_name << " \n"; 
                        broken_layouts_file.close();

                        // Move onto the next selection
                        break;
                    }

                    // Set up the angiogenesis solver (runs on deactivated VEGF function but seems to be necessary for main solver)
                    std::shared_ptr<AngiogenesisSolver<2> > p_angiogenesis_solver = AngiogenesisSolver<2>::Create();
                    std::shared_ptr<Owen2011SproutingRule<2> > p_sprouting_rule = Owen2011SproutingRule<2>::Create();
                    std::shared_ptr<Owen2011MigrationRule<2> > p_migration_rule = Owen2011MigrationRule<2>::Create();
                    p_angiogenesis_solver->SetMigrationRule(p_migration_rule);
                    p_angiogenesis_solver->SetSproutingRule(p_sprouting_rule);
                    p_sprouting_rule->SetDiscreteContinuumSolver(p_vegf_solver);
                    p_migration_rule->SetDiscreteContinuumSolver(p_vegf_solver);
                    p_angiogenesis_solver->SetVesselGridCalculator(p_grid_calc);
                    p_angiogenesis_solver->SetVesselNetwork(p_network);
                    
                    // Initialise the solver
                    auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                    p_microvessel_solver->SetVesselNetwork(p_network);
                    p_microvessel_solver->SetOutputFrequency(5_h);
                    // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                    p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                    p_microvessel_solver->AddDiscreteContinuumSolver(p_vegf_solver);
                    p_microvessel_solver->SetAngiogenesisSolver(p_angiogenesis_solver);

                    // Specifies which extracellular fields to update based on PDE
                    boost::shared_ptr<MicrovesselSimulationModifier<2> > p_microvessel_modifier =
                    boost::shared_ptr<MicrovesselSimulationModifier<2> >(new MicrovesselSimulationModifier<2> ());
                    p_microvessel_modifier->SetMicrovesselSolver(p_microvessel_solver);
                    std::vector<std::string> update_labels;
                    update_labels.push_back("oxygen");
                    p_microvessel_modifier->SetCellDataUpdateLabels(update_labels);

                    ////////////////////////////////////////////////////////////////
                    // Cells
                    ////////////////////////////////////////////////////////////////

                    // The full simulation is run as a typical Cell Based Chaste simulation
                    OnLatticeSimulation<2> simulator(*p_cell_population);
                    simulator.AddSimulationModifier(p_microvessel_modifier);

                    // Add a killer to remove apoptotic cells from the grid
                    boost::shared_ptr<ApoptoticCellKiller<2> > p_apoptotic_cell_killer(new ApoptoticCellKiller<2>(p_cell_population.get()));
                    simulator.AddCellKiller(p_apoptotic_cell_killer);
        
                    // Add a LQ RT killer
                    boost::shared_ptr<LQRadiotherapyCellKiller<2> > p_rt_killer =
                    boost::shared_ptr<LQRadiotherapyCellKiller<2> >(new LQRadiotherapyCellKiller<2> (p_cell_population.get()));
                    std::vector<QTime > rt_times;
                    rt_times.push_back(3600.0*0.0*unit::seconds);
                    // rt_times.push_back(3600.0*5.0*unit::seconds);
                    // rt_times.push_back(3600.0*10.0*unit::seconds);
                    // rt_times.push_back(3600.0*15.0*unit::seconds);
                    // rt_times.push_back(3600.0*24.0*unit::seconds);
                    // rt_times.push_back(3600.0*48.0*unit::seconds);
                    // rt_times.push_back(3600.0*72.0*unit::seconds);
                    p_rt_killer->SetTimeOfRadiation(rt_times);
                    // p_rt_killer->AddTimeOfRadiation(3600.0*96.0*unit::seconds);
                    p_rt_killer->SetDoseInjected(2.0*unit::gray);  // modify dose here
                    // p_rt_killer->SetDoseInjected(4.0*unit::gray);  // modify dose here
                    p_rt_killer->UseOer(true);  // turn OER on and off here

                    // Set up the radiobiological parameters for the Lewin model
                    // p_rt_killer->UseLewinModel(true);
                    // double alpha_max = 0.3;
                    // double alpha_beta_ratio = 10.0;

                    // Set up the radiobiological parameters for the Scott model
                    // p_rt_killer->UseScottModel(true);
                    // // double alpha_max = 0.2895*3  // Scott code
                    // double alpha_max = 0.3  // Grogan paper
                    // double alpha_beta_ratio = 10.0;
                    // // double alpha_beta_ratio = 8.64;  // alternative value
                    // p_rt_killer->SetOerAlphaMax(1.75);
                    // p_rt_killer->SetOerBetaMax(3.25);

                    // Set up the radiobiological parameters for the Wouters model
                    p_rt_killer->UseWoutersModel(true);
                    double alpha_max = 0.2;
                    double alpha_beta_ratio = 2.2;
                    p_rt_killer->SetOerAlphaMax(2.5);
                    p_rt_killer->SetOerBetaMax(3.0);

                    // If we're using OER                    
                    double beta_max = alpha_max/alpha_beta_ratio;
                    p_rt_killer->SetOerAlphaMin(1.0);
                    p_rt_killer->SetOerBetaMin(1.0);
                    p_rt_killer->SetOerConstant(0.004499758549541243*unit::mole_per_metre_cubed);  // K_OER (3.28 mmHg)
                    p_rt_killer->SetAlphaMax(alpha_max * unit::per_gray);
                    p_rt_killer->SetBetaMax(beta_max * unit::per_gray_squared);

                    // If OER is off, set the radiosensitivity parameters for different cells explicitly
                    p_rt_killer->SetCancerousRadiosensitivity(alpha_max * unit::per_gray, beta_max * unit::per_gray_squared);  // alpha, beta
                    // p_rt_killer->SetNormalRadiosensitivity(0.15 * unit::per_gray, 0.05 * unit::per_gray_squared);  // no normal cells in simulation
                    p_rt_killer->SetNormalRadiosensitivity(0.0 * unit::per_gray, 0.0 * unit::per_gray_squared);  // no normal cells in simulation

                    // Add the RT cell killer
                    simulator.AddCellKiller(p_rt_killer);

                    // If needed
                    // boost::shared_ptr<RandomCellKiller<2> > p_cell_killer(new RandomCellKiller<2>(p_cell_population.get(), 0.01));
                    // simulator.AddCellKiller(p_cell_killer);
                    // for (typename AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population.get()->Begin(); cell_iter != p_cell_population.get()->End(); ++cell_iter)
                    // {
                    //     // CheckAndLabelSingleCellForApoptosis(*cell_iter);
                    //     std::cout << "hi cell" << (*cell_iter)->GetCellData()->GetItem("oxygen") << std::endl;

                    //     //  std::cout <<  BaseUnits::Instance()->GetReferenceConcentrationScale() << std::endl;


                    //     //     double oxygen_concentration;

                    //     // std::cout << "hi cell" << oxygen_concentration << std::endl;
                    // }

                    // Add another modifier for updating cell cycle quantities.
                    boost::shared_ptr<Owen2011TrackingModifier<2> > p_owen11_tracking_modifier(new Owen2011TrackingModifier<2>);
                    simulator.AddSimulationModifier(p_owen11_tracking_modifier);

                    // Print the average oxygenation
                    // std::vector<double> solution = p_oxygen_solver->GetSolution();
                    // double average_oxygen = 0.0;
                    // for(unsigned jdx=0;jdx<solution.size();jdx++)
                    // {
                    //     average_oxygen += solution[jdx];
                    // }
                    // average_oxygen /= double(solution.size());
                    // std::cout << "Average oxygen: " << average_oxygen << std::endl;

                    // Set up the simulation time and run it
                    simulator.SetOutputDirectory(str_directory_name);
                    // simulator.SetSamplingTimestepMultiple(10);  // get sample every x steps of Dt (multiple*dt = time between s)
                    simulator.SetSamplingTimestepMultiple(1);
                    // simulator.SetDt(1);  // simulation time step in hours
                    simulator.SetDt(1.0);
                    // simulator.SetEndTime(24.0*7.0);  // one week
                    simulator.SetEndTime(1.0*1.0);  // five hours
                    simulator.Solve();

                    // Print the average oxygenation
                    // std::vector<double> solution = p_oxygen_solver->GetSolution();
                    // double average_oxygen = 0.0;
                    // for(unsigned jdx=0;jdx<solution.size();jdx++)
                    // {
                    //     average_oxygen += solution[jdx];
                    // }
                    // average_oxygen /= double(solution.size());
                    // std::cout << "Median oxygen: " << average_oxygen << std::endl;

                    // Print the median oxygenation
                    std::vector<double> solution = p_oxygen_solver->GetSolution();
                    std::sort(solution.begin(), solution.end());
                    double median_oxygen;
                    size_t solution_size = solution.size();
                    if (solution_size % 2 == 0) 
                    {
                        // Even number of elements
                        median_oxygen = (solution[solution_size / 2 - 1] + solution[solution_size / 2]) / 2.0;
                    } 
                    else 
                    {
                        // Odd number of elements
                        median_oxygen = solution[solution_size / 2];
                    }
                    std::cout << "Median oxygen (nM): " << median_oxygen << std::endl;

                    // Write the output network file (visualise with Paraview: set Filters->Alphabetical->Tube)
                    // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                    // p_network->Write(output_file);

                    // // Remove the shortest or lowest-flow vessel
                    // vessels = p_network->GetVessels();
                    // // QLength minimum_length = 10000.0_um;
                    // QFlowRate minimum_flow = 1.e-0*unit::metre_cubed_per_second;
                    // unsigned int minimum_index = 0;
                    // for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                    // {
                    //     // Exclude inlets and outlets
                    //     if (!(vessels[vessel_index]->GetStartNode()->GetFlowProperties()->IsInputNode()
                    //     || vessels[vessel_index]->GetStartNode()->GetFlowProperties()->IsOutputNode()
                    //     || vessels[vessel_index]->GetEndNode()->GetFlowProperties()->IsInputNode()
                    //     || vessels[vessel_index]->GetEndNode()->GetFlowProperties()->IsOutputNode()))
                    //     {   
                    //         // Get the current segment's length
                    //         // QLength current_length = vessels[vessel_index]->GetLength();

                    //         // // If the current length is less than the minimum length, record the new minimum
                    //         // if (current_length < minimum_length)
                    //         // {
                    //         //     minimum_length = current_length;
                    //         //     minimum_index = vessel_index;
                    //         // }                  

                    //         // Get the current segment's flow rate
                    //         QFlowRate current_flow = fabs(vessels[vessel_index]->GetFlowProperties()->GetFlowRate());

                    //         // If the current flow rate is non-zero and less than the minimum flow rate, record the new minimum
                    //         if (current_flow < minimum_flow)
                    //         {
                    //             minimum_flow = current_flow;
                    //             minimum_index = vessel_index;
                    //         }
                    //     }
                    // }
                    // p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel

                    // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                    ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                    ParameterCollection::Instance()->Destroy();
                    BaseUnits::Instance()->Destroy();
                    SimulationTime::Instance()->Destroy();
                    SimulationTime::Instance()->SetStartTime(0.0);

                    // If there is no oxygen in the network, don't prune any more vessels
                    double total_oxygen = 0.0;
                    for (unsigned i = 0; i < solution_size; ++i) 
                    {
                        total_oxygen += solution[i];
                    }
                    std::cout << "Total oxygen (nM): " << total_oxygen << std::endl;
                    if (total_oxygen < tolerance2)
                    {
                        break;
                    }

                    // Remove the shortest or lowest-flow vessel
                    vessels = p_network->GetVessels();
                    QFlowRate minimum_flow = 1.e-0*unit::metre_cubed_per_second;
                    std::vector<unsigned int> minimum_flow_indices;

                    // Step 1: Find the minimum flow rate and collect all vessels with this flow rate
                    for (unsigned vessel_index = 0; vessel_index < vessels.size(); vessel_index++) 
                    {
                        // Exclude inlets and outlets
                        if (!(vessels[vessel_index]->GetStartNode()->GetFlowProperties()->IsInputNode() ||
                            vessels[vessel_index]->GetStartNode()->GetFlowProperties()->IsOutputNode() ||
                            vessels[vessel_index]->GetEndNode()->GetFlowProperties()->IsInputNode() ||
                            vessels[vessel_index]->GetEndNode()->GetFlowProperties()->IsOutputNode())) 
                        {
                            QFlowRate current_flow = fabs(vessels[vessel_index]->GetFlowProperties()->GetFlowRate());
                            if (current_flow < minimum_flow) 
                            {
                                minimum_flow = current_flow;
                                minimum_flow_indices.clear();
                                minimum_flow_indices.push_back(vessel_index);
                            } 
                            else if (current_flow == minimum_flow) 
                            {
                                minimum_flow_indices.push_back(vessel_index);
                            }
                        }
                    }
                    // Step 2: Check for vessels with a degree one node
                    std::vector<unsigned int> degree_one_indices;
                    for (unsigned int index : minimum_flow_indices) 
                    {
                        if (vessels[index]->GetStartNode()->GetNumberOfSegments() == 1 || vessels[index]->GetEndNode()->GetNumberOfSegments() == 1) 
                        {
                            degree_one_indices.push_back(index);
                        }
                    }
                    // Step 3: If there are multiple vessels with degree one nodes, find the shortest one
                    unsigned int vessel_to_remove = minimum_flow_indices[0];
                    if (!degree_one_indices.empty()) 
                    {
                        vessel_to_remove = degree_one_indices[0];
                        QLength minimum_length = vessels[degree_one_indices[0]]->GetLength();
                        for (unsigned int index : degree_one_indices) 
                        {
                            QLength current_length = vessels[index]->GetLength();
                            if (current_length < minimum_length) {
                                minimum_length = current_length;
                                vessel_to_remove = index;
                            }
                        }
                    } 
                    else if (minimum_flow_indices.size() > 1) 
                    {
                        // If no degree one nodes, find the shortest vessel among those with minimum flow
                        vessel_to_remove = minimum_flow_indices[0];
                        QLength minimum_length = vessels[minimum_flow_indices[0]]->GetLength();
                        for (unsigned int index : minimum_flow_indices) 
                        {
                            QLength current_length = vessels[index]->GetLength();
                            if (current_length < minimum_length) 
                            {
                                minimum_length = current_length;
                                vessel_to_remove = index;
                            }
                        }
                    }
                    // Remove the selected vessel
                    p_network->RemoveVessel(vessels[vessel_to_remove], true);

                    // If simulation doesn't converge, move on to next layout
                    if (broken_solver == 1)
                    {
                        break;
                    }
                }
            }
        }
        
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Simulate an experimentally-acquired healthy network on a PDE grid with flow and O2 and cells
    void xTestHealthyNetwork2D() 
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Set network in filename
        std::string network_name = "TestHealthyNetworkTissue/";

        // Choose the vascular parameters
        // unsigned NumberOfLayouts = 1000;  // number of different point layouts to run simulations with (add an extra selection for a demo)
        double dimless_domain_size_x = 850.22; 
        double dimless_domain_size_y = 850.22;  // based on image dimensions
        // QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        QDynamicViscosity viscosity = 1.96*1.e-3*unit::poiseuille;
        // double initial_haematocrit = 0.45;        
        double initial_haematocrit = 0.01;  // conks out somewhere between 0.35 and 0.4
        // unsigned n_vessels = 382;  // number of non-inlet/outlet vessels from which to select ones to make thin
        // double percToKill = 0.2;  // percentage of vessels to kill
        // unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
        // unsigned ToBeKilled = 0.5*((3*NumberOfSeedPoints)-6);  // number to kill
        // unsigned ToBeKilled = NumberOfSeedPoints;  // number to kill
        unsigned ToBeKilled = 0;  // number to kill
        
        // Match this based on cell size
        QLength cell_grid_spacing = 100.0_um;

        // Seed the random number generator
        RandomNumberGenerator::Instance()->Reseed(12345);

        // Create the log file for broken simulations
        std::ofstream broken_layouts_file;
        broken_layouts_file.open("/scratch/narain/testoutput/TestHealthyNetworkTissue/broken_layouts.txt");
        broken_layouts_file.close();        

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);
        // BaseUnits::Instance()->SetReferenceTimeScale(1_h);
        QTime reference_time = 3600.0*unit::seconds;  // 1 hour
        BaseUnits::Instance()->SetReferenceTimeScale(reference_time);

        ////////////////////////////////////////////////////////////////
        // Vasculature
        ////////////////////////////////////////////////////////////////

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {     
            // Initialise the simulation
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            VesselPtr<2> p_current_vessel;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            // std::vector<std::vector<unsigned> > Order;

                // Set up flag for broken solver
                unsigned broken_solver = 0;
                
                // Read the network layout from an edge matrix
                std::shared_ptr<VesselNetwork<2> > p_network;
                VesselNetworkGenerator<2> network_generator;
                std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/edges/biological/CombinedSkeletons.txt");
                // std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/edges/biological/Skeleton_66.txt");
                std::vector<std::vector<double> > rEdgesMatrix;
                string line;
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
                p_network = network_generator.GenerateVoronoiNetwork(rEdgesMatrix);

                // Read the network layout from a VTP file (doesn't work yet)
                // std::shared_ptr<VesselNetworkReader<2> > p_vessel_reader = std::shared_ptr<VesselNetworkReader<2> >(new VesselNetworkReader<2> );
                // FileFinder file_finder = FileFinder("/projects/MicrovesselChaste/test/data/interesting/FinalHaematocrit.vtp", RelativeTo::ChasteSourceRoot);
                // p_vessel_reader->SetFileName(file_finder.GetAbsolutePath());
                // p_vessel_reader->SetMergeCoincidentPoints(true);
                // // p_vessel_reader->SetTargetSegmentLength(40.0e-6*unit::metres);
                // std::shared_ptr<VesselNetwork<2> >  p_network = p_vessel_reader->Read();
                // auto p_clean_file_handler = std::make_shared<OutputFileHandler>(network_name, true);
                // p_network->Write(p_clean_file_handler->GetOutputDirectoryFullPath() + "cleaned_network.vtp");
                // std::cout << "This is line number: " << __LINE__ << std::endl;
                // QLength short_vessel_cutoff = 1.0_um;
                // bool remove_end_vessels_only = true;
                // p_network->RemoveShortVessels(short_vessel_cutoff, remove_end_vessels_only);
                // p_network->UpdateAll();
                // p_network->MergeCoincidentNodes();
                // p_network->UpdateAll();


                // Manually set the inlets and outlets
                vessels = p_network->GetVessels();
                p_network->GetNode(187)->GetFlowProperties()->SetIsInputNode(true);
                p_network->GetNode(187)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));
                p_network->GetNode(1084)->GetFlowProperties()->SetIsOutputNode(true);
                p_network->GetNode(1084)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));

                p_network->GetNode(1557)->GetFlowProperties()->SetIsInputNode(true);
                p_network->GetNode(1557)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));
                p_network->GetNode(1091)->GetFlowProperties()->SetIsOutputNode(true);
                p_network->GetNode(1091)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));

                // Assign the first vessel as the input and the last as the output and remove all other dead ends (i.e., make the network single inlet and outlet)
                // double tolerance = 0.001;  // for location of inlet/outlet nodes
                // vessels = p_network->GetVessels();
                // for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                // {
                //     if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
                //     {
                //         if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < tolerance)
                //         {
                //             (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
                //             (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
                //         }
                //         if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                //         {
                //             (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
                //             (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                //         }
                //     }
                //     if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)   	
                //     {
                //         if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < tolerance)
                //         {
                //             (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
                //             (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
                //         }
                //         if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                //         {
                //             (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
                //             (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                //         }
                //     }
                // }

                // Remove diameter heterogeneity
                auto p_segment = p_network->GetVesselSegments()[0];
                p_segment->SetRadius(7.5_um);  // same as Mantegazza
                VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   

                // Prune all vessels up to specified dose 
                // for(unsigned KilledVessels=0; KilledVessels < ToBeKilled; KilledVessels++)
                // { 
                //     // Remove the shortest vessel
                //     vessels = p_network->GetVessels();
                //     QLength minimum_length = 1000.0_um;
                //     unsigned int minimum_index = 0;
                //     for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                //     {
                //         // Exclude inlets and outlets
                //         if (!(vessels[vessel_index]->GetStartNode()->GetFlowProperties()->IsInputNode()
                //         || vessels[vessel_index]->GetStartNode()->GetFlowProperties()->IsOutputNode()
                //         || vessels[vessel_index]->GetEndNode()->GetFlowProperties()->IsInputNode()
                //         || vessels[vessel_index]->GetEndNode()->GetFlowProperties()->IsOutputNode()))
                //         {   
                //             // Get the current segment's length
                //             QLength current_length = vessels[vessel_index]->GetLength();
                            
                //             // If the current length is less than the minimum length, record the new minimum
                //             if (current_length < minimum_length)
                //             {
                //                 minimum_length = current_length;
                //                 minimum_index = vessel_index;
                //             }                  
                //         }
                //     }
                //     p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel
                    
                //     // Display status message
                //     std::cout << "Now killed " << KilledVessels+1 << " vessels." << std::endl;  

                // }                       

                // Set the haematocrit solver
                std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                std::string solver_name;
                if (h_solver==1)
                {
                    solver_name = "ConstantHaematocrit"; 
                    std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
                    p_abstract_haematocrit_solver = p_haematocrit_solver;     
                }
                else if (h_solver==2)
                {
                    solver_name = "PriesHaematocrit"; 
                    std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
                    p_abstract_haematocrit_solver = p_haematocrit_solver;                
                }
                else if (h_solver==3)
                {
                    solver_name = "MemoryHaematocrit"; 
                    std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
                    p_abstract_haematocrit_solver = p_haematocrit_solver;     
                }
                else if (h_solver==4)
                {
                    solver_name = "FungHaematocrit"; 
                    std::cout << "Now using FungHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
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

                // // Prune all vessels up to specified dose 
                // for(unsigned KilledVessels=0; KilledVessels <= ToBeKilled; KilledVessels++)
                // { 
                    // // Display status message
                    // std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                    // Set filename
                    std::stringstream kill_stream;
                    // kill_stream << std::fixed << std::setprecision(0) << KilledVessels;
                    kill_stream << std::fixed << std::setprecision(0) << ToBeKilled;                      
                    std::string kill_string = kill_stream.str();                    
                    std::string file_name = network_name + solver_name + "/Kills" + kill_string;                        
                    std::string str_directory_name = file_name;
                    std::cout << str_directory_name << std::endl;
                    auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                    // Set up the grid for the finite difference solver
                    auto p_grid = RegularGrid<2>::Create();
                    p_grid->SetSpacing(cell_grid_spacing);
                    c_vector<unsigned, 3> dimensions;
                    dimensions[0] = unsigned((domain_side_length_x)/(cell_grid_spacing))+2; // num x
                    dimensions[1] = unsigned((domain_side_length_y)/(cell_grid_spacing))+2; // num_y
                    dimensions[2] = 1;
                    p_grid->SetDimensions(dimensions);

                    // Set up the PDE
                    auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));
                    auto p_cell_oxygen_sink = CellBasedDiscreteSource<2>::Create();  // Set the cell terms
                    // p_cell_oxygen_sink->SetLinearInUConsumptionRatePerCell(13.0);
                    p_cell_oxygen_sink->SetLinearInUConsumptionRatePerCell(Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));
                    p_oxygen_pde->AddDiscreteSource(p_cell_oxygen_sink);

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
                    unsigned max_max_iter = 10000; 
                    unsigned max_iter = 1000;  // 1000 
                    double prev_max_difference = 0.0;
                    double tolerance2 = 1.e-10;
                    std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                    std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                    for(unsigned idx=0;idx<max_max_iter;idx++)
                    {
                        // Run the solvers (order of calculators matters!)
                        p_viscosity_calculator->Calculate();
                        p_impedance_calculator->Calculate();
                        flow_solver.SetUp();
                        flow_solver.Solve(); std::cout << "This is line number: " << __LINE__ << std::endl;
                        p_abstract_haematocrit_solver->Calculate(); std::cout << "This is line number: " << __LINE__ << std::endl;

                        // Get the residual
                        double max_difference = 0.0;
                        double h_for_max = 0.0;
                        double prev_for_max = 0.0;
                        for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                        {
                            // Set segments with no flow to be dead (this is only problem for inlet really)
                            if (fabs(segments[jdx]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                            {
                                segments[jdx]->GetFlowProperties()->SetViscosity(0.0);
                                segments[jdx]->GetFlowProperties()->SetImpedance(0.0);
                                segments[jdx]->GetFlowProperties()->SetHaematocrit(0.0);
                            }    
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
                            // broken_solver = 0;
                            break;
                        }
                        else
                        {
                            std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                            if(idx%100==0)
                            {
                                std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                                std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                                p_network->Write(output_file);
                            }
                        }

                        // If there is no convergence after all the iterations, print the error message.
                        if(idx>=max_iter-1 && max_difference >= prev_max_difference)
                        {
                            std::cout << "Problem encountered in " << str_directory_name << std::endl;
                            error_log << "\n Problem encountered in " << str_directory_name << std::endl; 
                            broken_solver = 1;
                            break;
                        }                                
                        prev_max_difference = max_difference;
                    }

                    // If simulation doesn't converge, move on to next layout and log problem 
                    if (broken_solver == 1)
                    {
                        broken_layouts_file.open("/scratch/narain/testoutput/TestHealthyNetworkTissue/broken_layouts.txt", std::ios_base::app);
                        broken_layouts_file << str_directory_name << " \n"; 
                        broken_layouts_file.close();

                        // Move onto the next selection
                        break;
                    }

                    // Set up the cell populations
                    std::shared_ptr<GridCalculator<2> > p_grid_calc = GridCalculator<2>::Create();
                    p_grid_calc->SetGrid(p_grid);
                    std::shared_ptr<Owen11CellPopulationGenerator<2> > p_cell_population_genenerator = Owen11CellPopulationGenerator<2>::Create();
                    p_cell_population_genenerator->SetGridCalculator(p_grid_calc);
                    p_cell_population_genenerator->SetVesselNetwork(p_network);
                    QLength tumour_radius(10000.0 * unit::microns);  // just make it bigger than the domain to ensure all cells are tumour cells 
                    p_cell_population_genenerator->SetTumourRadius(tumour_radius);
                    std::shared_ptr<CaBasedCellPopulation<2> > p_cell_population = p_cell_population_genenerator->Update();
                    
                    // Set up the VEGF PDE as we set the O2 one (currently turned off)
                    auto p_vegf_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    p_vegf_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpVegfDiffusivity->GetValue("User"));
                    // p_vegf_pde->SetIsotropicDiffusionConstant(0.0);
                    p_vegf_pde->SetContinuumLinearInUTerm(-Owen11Parameters::mpVegfDecayRate->GetValue("User"));
                
                    // Set up a map for different release rates depending on cell type. Also include a threshold intracellular VEGF below which there is no release (currently turned off)
                    auto p_normal_and_quiescent_cell_source = CellStateDependentDiscreteSource<2>::Create();  // initialise the source (which combines normal and quiescent cells)
                    std::map<unsigned, QConcentrationFlowRate > normal_and_quiescent_cell_rates;  // map to store secretion rates
                    std::map<unsigned, QConcentration > normal_and_quiescent_cell_rate_thresholds;  // map to store secretion thresholds
                    MAKE_PTR(QuiescentCancerCellMutationState, p_quiescent_cancer_state);
                    MAKE_PTR(WildTypeCellMutationState, p_normal_cell_state);
                    normal_and_quiescent_cell_rates[p_normal_cell_state->GetColour()] = Owen11Parameters::mpCellVegfSecretionRate->GetValue("User");  // set the rate for normal cells
                    normal_and_quiescent_cell_rate_thresholds[p_normal_cell_state->GetColour()] = 0.27*unit::mole_per_metre_cubed;  // set the threshold for normal cells
                    // normal_and_quiescent_cell_rates[p_normal_cell_state->GetColour()] = 0.00000000000000000001;  // set the rate for normal cells
                    // normal_and_quiescent_cell_rate_thresholds[p_normal_cell_state->GetColour()] = 1000000_M;  // set the threshold for normal cells
                    normal_and_quiescent_cell_rates[p_quiescent_cancer_state->GetColour()] = Owen11Parameters::mpCellVegfSecretionRate->GetValue("User");  // set the rate for quiescent cells
                    normal_and_quiescent_cell_rate_thresholds[p_quiescent_cancer_state->GetColour()] = 0_M;  // set the threshold for quiescent cells
                    // normal_and_quiescent_cell_rates[p_quiescent_cancer_state->GetColour()] = 0.00000000000000000001;  // set the rate for quiescent cells
                    // normal_and_quiescent_cell_rate_thresholds[p_quiescent_cancer_state->GetColour()] = 1000000_M;  // set the threshold for quiescent cells
                    p_normal_and_quiescent_cell_source->SetStateRateMap(normal_and_quiescent_cell_rates);  // configure the map
                    p_normal_and_quiescent_cell_source->SetLabelName("VEGF");  // set the label for the source
                    p_normal_and_quiescent_cell_source->SetStateRateThresholdMap(normal_and_quiescent_cell_rate_thresholds);
                    // p_vegf_pde->AddDiscreteSource(p_normal_and_quiescent_cell_source);

                    // Print the cell types and their colours
                    MAKE_PTR(CancerCellMutationState, p_cancer_state);
                    // MAKE_PTR(StemCellProliferativeType, p_stem_state);
                    std::cout << "WildTypeCellMutationState has Legacy Cell type " <<  p_normal_cell_state->GetColour() << std::endl;
                    std::cout << "CancerCellMutationState has Legacy Cell type " <<  p_cancer_state->GetColour() << std::endl;
                    std::cout << "QuiescentCancerCellMutationState has Legacy Cell type " <<  p_quiescent_cancer_state->GetColour() << std::endl;
                    // std::cout << "Quiescent cancer state has Legacy Cell type " <<  p_stem_state->GetColour() << std::endl;

                    // Add a vessel related VEGF sink (currently turned off)
                    auto p_vessel_vegf_sink = VesselBasedDiscreteSource<2>::Create();
                    p_vessel_vegf_sink->SetReferenceConcentration(0.0_M);
                    p_vessel_vegf_sink->SetVesselPermeability(Owen11Parameters::mpVesselVegfPermeability->GetValue("User"));
                    // p_vessel_vegf_sink->SetVesselPermeability(0.0);
                    // p_vegf_pde->AddDiscreteSource(p_vessel_vegf_sink);

                    // Set up a finite difference solver as before for the VEGF
                    auto p_vegf_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    p_vegf_solver->SetPde(p_vegf_pde);
                    p_vegf_solver->SetLabel("VEGF_Extracellular");  // set the label for the external field
                    p_vegf_solver->SetGrid(p_grid);

                    // Set up the angiogenesis solver (runs on deactivated VEGF function but seems to be necessary for main solver)
                    std::shared_ptr<AngiogenesisSolver<2> > p_angiogenesis_solver = AngiogenesisSolver<2>::Create();
                    std::shared_ptr<Owen2011SproutingRule<2> > p_sprouting_rule = Owen2011SproutingRule<2>::Create();
                    std::shared_ptr<Owen2011MigrationRule<2> > p_migration_rule = Owen2011MigrationRule<2>::Create();
                    p_angiogenesis_solver->SetMigrationRule(p_migration_rule);
                    p_angiogenesis_solver->SetSproutingRule(p_sprouting_rule);
                    p_sprouting_rule->SetDiscreteContinuumSolver(p_vegf_solver);
                    p_migration_rule->SetDiscreteContinuumSolver(p_vegf_solver);
                    p_angiogenesis_solver->SetVesselGridCalculator(p_grid_calc);
                    p_angiogenesis_solver->SetVesselNetwork(p_network);
                    
                    // Initialise the solver
                    auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                    p_microvessel_solver->SetVesselNetwork(p_network);
                    p_microvessel_solver->SetOutputFrequency(5_h);
                    // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                    p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                    p_microvessel_solver->AddDiscreteContinuumSolver(p_vegf_solver);
                    p_microvessel_solver->SetAngiogenesisSolver(p_angiogenesis_solver);

                    // Specifies which extracellular fields to update based on PDE
                    boost::shared_ptr<MicrovesselSimulationModifier<2> > p_microvessel_modifier =
                    boost::shared_ptr<MicrovesselSimulationModifier<2> >(new MicrovesselSimulationModifier<2> ());
                    p_microvessel_modifier->SetMicrovesselSolver(p_microvessel_solver);
                    std::vector<std::string> update_labels;
                    update_labels.push_back("oxygen");
                    p_microvessel_modifier->SetCellDataUpdateLabels(update_labels);

                    ////////////////////////////////////////////////////////////////
                    // Cells
                    ////////////////////////////////////////////////////////////////

                    // The full simulation is run as a typical Cell Based Chaste simulation
                    OnLatticeSimulation<2> simulator(*p_cell_population);
                    simulator.AddSimulationModifier(p_microvessel_modifier);

                    // Add a killer to remove apoptotic cells from the grid
                    boost::shared_ptr<ApoptoticCellKiller<2> > p_apoptotic_cell_killer(new ApoptoticCellKiller<2>(p_cell_population.get()));
                    simulator.AddCellKiller(p_apoptotic_cell_killer);
        
                    // Add a LQ RT killer
                    boost::shared_ptr<LQRadiotherapyCellKiller<2> > p_rt_killer =
                    boost::shared_ptr<LQRadiotherapyCellKiller<2> >(new LQRadiotherapyCellKiller<2> (p_cell_population.get()));
                    p_rt_killer->SetDoseInjected(2.0*unit::gray);  // modify dose here
                    // p_rt_killer->SetDoseInjected(4.0*unit::gray);  // modify dose here
                    // p_rt_killer->SetCancerousRadiosensitivity(0.3 * unit::per_gray, 0.03 * unit::per_gray_squared);  // only used if OER is off; currently same as default parameters
                    // p_rt_killer->SetNormalRadiosensitivity(0.15 * unit::per_gray, 0.05 * unit::per_gray_squared);  // no normal cells in simulation
                    p_rt_killer->SetOerAlphaMax(1.75);
                    p_rt_killer->SetOerAlphaMin(1.0);
                    p_rt_killer->SetOerBetaMax(3.25);
                    p_rt_killer->SetOerBetaMin(1.0);
                    // p_rt_killer->SetOerConstant(oxygen_solubility_at_stp * 3.28_Pa);  // adjust the radiobiological threshold here
                    p_rt_killer->SetOerConstant(0.004499758549541243*unit::mole_per_metre_cubed);  // K_OER (3.28 mmHg)
                    p_rt_killer->SetAlphaMax(0.3 * unit::per_gray);
                    p_rt_killer->SetBetaMax(0.03 * unit::per_gray_squared);
                    p_rt_killer->UseOer(true);  // turn OER on and off here

                    // Set up the dosage times
                    std::vector<QTime > rt_times;
                    rt_times.push_back(3600.0*1.0*unit::seconds);
                    // rt_times.push_back(3600.0*5.0*unit::seconds);
                    // rt_times.push_back(3600.0*10.0*unit::seconds);
                    // rt_times.push_back(3600.0*15.0*unit::seconds);
                    // rt_times.push_back(3600.0*24.0*unit::seconds);
                    // rt_times.push_back(3600.0*48.0*unit::seconds);
                    // rt_times.push_back(3600.0*72.0*unit::seconds);
                    p_rt_killer->SetTimeOfRadiation(rt_times);
                    // p_rt_killer->AddTimeOfRadiation(3600.0*96.0*unit::seconds);
                    simulator.AddCellKiller(p_rt_killer);

                    // If needed
                    // boost::shared_ptr<RandomCellKiller<2> > p_cell_killer(new RandomCellKiller<2>(p_cell_population.get(), 0.01));
                    // simulator.AddCellKiller(p_cell_killer);
                    // for (typename AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population.get()->Begin(); cell_iter != p_cell_population.get()->End(); ++cell_iter)
                    // {
                    //     // CheckAndLabelSingleCellForApoptosis(*cell_iter);
                    //     std::cout << "hi cell" << (*cell_iter)->GetCellData()->GetItem("oxygen") << std::endl;

                    //     //  std::cout <<  BaseUnits::Instance()->GetReferenceConcentrationScale() << std::endl;


                    //     //     double oxygen_concentration;

                    //     // std::cout << "hi cell" << oxygen_concentration << std::endl;
                    // }

                    // Add another modifier for updating cell cycle quantities.
                    boost::shared_ptr<Owen2011TrackingModifier<2> > p_owen11_tracking_modifier(new Owen2011TrackingModifier<2>);
                    simulator.AddSimulationModifier(p_owen11_tracking_modifier);

                    // Print the average oxygenation
                    // std::vector<double> solution = p_oxygen_solver->GetSolution();
                    // double average_oxygen = 0.0;
                    // for(unsigned jdx=0;jdx<solution.size();jdx++)
                    // {
                    //     average_oxygen += solution[jdx];
                    // }
                    // average_oxygen /= double(solution.size());
                    // std::cout << "Average oxygen: " << average_oxygen << std::endl;

                    // Set up the simulation time and run it
                    simulator.SetOutputDirectory(str_directory_name);
                    // simulator.SetSamplingTimestepMultiple(10);  // get sample every x steps of Dt (multiple*dt = time between s)
                    simulator.SetSamplingTimestepMultiple(1);
                    // simulator.SetDt(1);  // simulation time step in hours
                    simulator.SetDt(1.0);
                    // simulator.SetEndTime(24*7);  // end time in hours
                    // simulator.SetEndTime(24);  // end time in hours
                    simulator.SetEndTime(10.0);  // ten hours
                    // simulator.SetEndTime(24.0*7.0);  // one week
                    simulator.Solve();

                    // Write the output network file (visualise with Paraview: set Filters->Alphabetical->Tube)
                    std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                    p_network->Write(output_file);

                    // // Remove the shortest vessel
                    // vessels = p_network->GetVessels();
                    // QLength minimum_length = 1000.0_um;
                    // unsigned int minimum_index = 0;
                    // for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                    // {
                    //     // Exclude inlets and outlets
                    //     if (!(vessels[vessel_index]->GetStartNode()->GetFlowProperties()->IsInputNode()
                    //     || vessels[vessel_index]->GetStartNode()->GetFlowProperties()->IsOutputNode()
                    //     || vessels[vessel_index]->GetEndNode()->GetFlowProperties()->IsInputNode()
                    //     || vessels[vessel_index]->GetEndNode()->GetFlowProperties()->IsOutputNode()))
                    //     {   
                    //         // Get the current segment's length
                    //         QLength current_length = vessels[vessel_index]->GetLength();
                            
                    //         // If the current length is less than the minimum length, record the new minimum
                    //         if (current_length < minimum_length)
                    //         {
                    //             minimum_length = current_length;
                    //             minimum_index = vessel_index;
                    //         }                  
                    //     }
                    // }
                    // p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel

                    // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                    ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                    ParameterCollection::Instance()->Destroy();
                    BaseUnits::Instance()->Destroy();
                    SimulationTime::Instance()->Destroy();
                // }
                // If simulation doesn't converge, move on to next layout
                if (broken_solver == 1)
                {
                    break;
                }
            // }
        }

        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Simulate an experimentally-acquired metastatic network on a PDE grid with flow and O2 and experimentally-acquired cells
    void xTestMetastaticImage2D() 
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Set network in filename
        std::string network_name = "TestMetastaticNetworkTissue/";

        // Choose the domain dimensions (lectin monotile)
        double dimless_domain_size_x = 1214.56; 
        double dimless_domain_size_y = 1214.56;  // based on image dimensions
                
        // Choose the domain dimensions (control)
        // double dimless_domain_size_x = 850.22; 
        // double dimless_domain_size_y = 850.22;  // based on image dimensions

        // Set the initial blood parameters
        QDynamicViscosity viscosity = 1.96*1.e-3*unit::poiseuille;
        double initial_haematocrit = 0.01;  // conks out somewhere between 0.35 and 0.4
        
        // Change the resolution for the PDE fields
        QLength cell_grid_spacing = 1.0_um;

        // Seed the random number generator
        RandomNumberGenerator::Instance()->Reseed(12345);

        // Create the log file for broken simulations
        std::ofstream broken_layouts_file;
        broken_layouts_file.open("/scratch/narain/testoutput/TestMetastaticNetworkTissue/broken_layouts.txt");
        broken_layouts_file.close();        

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);
        QTime reference_time = 3600.0*unit::seconds;  // 1 hour
        BaseUnits::Instance()->SetReferenceTimeScale(reference_time);

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {     
            // Initialise the simulation
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            VesselPtr<2> p_current_vessel;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            // std::vector<std::vector<unsigned> > Order;

            // Set up flag for broken solver
            unsigned broken_solver = 0;
            
            // Read the network layout from an edge matrix
            std::shared_ptr<VesselNetwork<2> > p_network;
            VesselNetworkGenerator<2> network_generator;
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/edges/biological/lectin_monotile/CombinedSkeletons.txt");
            std::vector<std::vector<double> > rEdgesMatrix;
            string line;
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
            p_network = network_generator.GenerateVoronoiNetwork(rEdgesMatrix);

            // Manually set the inlets and outlets (lectin monotile)
            vessels = p_network->GetVessels();
            p_network->GetNode(1711)->GetFlowProperties()->SetIsInputNode(true);
            p_network->GetNode(1711)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));
            p_network->GetNode(1405)->GetFlowProperties()->SetIsOutputNode(true);
            p_network->GetNode(1405)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));

            //  Manually set the inlets and outlets (control)
            // p_network->GetNode(1557)->GetFlowProperties()->SetIsInputNode(true);
            // p_network->GetNode(1557)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));
            // p_network->GetNode(1091)->GetFlowProperties()->SetIsOutputNode(true);
            // p_network->GetNode(1091)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));

            // Remove diameter heterogeneity
            auto p_segment = p_network->GetVesselSegments()[0];
            p_segment->SetRadius(7.5_um);  // same as Mantegazza
            VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);                   

            // Set the haematocrit solver
            std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
            std::string solver_name;
            if (h_solver==1)
            {
                solver_name = "ConstantHaematocrit"; 
                std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
                p_abstract_haematocrit_solver = p_haematocrit_solver;     
            }
            else if (h_solver==2)
            {
                solver_name = "PriesHaematocrit"; 
                std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
                p_abstract_haematocrit_solver = p_haematocrit_solver;                
            }
            else if (h_solver==3)
            {
                solver_name = "MemoryHaematocrit"; 
                std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
                p_abstract_haematocrit_solver = p_haematocrit_solver;     
            }
            else if (h_solver==4)
            {
                solver_name = "FungHaematocrit"; 
                std::cout << "Now using FungHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
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

            // Set filename
            std::string file_name = network_name + solver_name;                      
            std::string str_directory_name = file_name;
            std::cout << str_directory_name << std::endl;
            auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

            // Set up the grid for the finite difference solver
            auto p_grid = RegularGrid<2>::Create();
            p_grid->SetSpacing(cell_grid_spacing);
            c_vector<unsigned, 3> dimensions;
            dimensions[0] = unsigned((domain_side_length_x)/(cell_grid_spacing))+2; // num x
            dimensions[1] = unsigned((domain_side_length_y)/(cell_grid_spacing))+2; // num_y
            dimensions[2] = 1;
            p_grid->SetDimensions(dimensions);

            // Set up the PDE
            auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
            p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
            auto p_cell_oxygen_sink = CellBasedDiscreteSource<2>::Create();  // Set the cell terms
            p_cell_oxygen_sink->SetLinearInUConsumptionRatePerCell(Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));
            p_oxygen_pde->AddDiscreteSource(p_cell_oxygen_sink);

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

            // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
            unsigned max_max_iter = 10000; 
            unsigned max_iter = 1000;  // 1000 
            double prev_max_difference = 0.0;
            double tolerance2 = 1.e-10;
            std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
            std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
            for(unsigned idx=0;idx<max_max_iter;idx++)
            {
                // Run the solvers (order of calculators matters!)
                p_viscosity_calculator->Calculate();
                p_impedance_calculator->Calculate();
                flow_solver.SetUp();
                flow_solver.Solve(); 
                p_abstract_haematocrit_solver->Calculate();

                // Get the residual
                double max_difference = 0.0;
                double h_for_max = 0.0;
                double prev_for_max = 0.0;
                for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                {
                    // Set segments with no flow to be dead (this is only problem for inlet really)
                    if (fabs(segments[jdx]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                    {
                        segments[jdx]->GetFlowProperties()->SetViscosity(0.0);
                        segments[jdx]->GetFlowProperties()->SetImpedance(0.0);
                        segments[jdx]->GetFlowProperties()->SetHaematocrit(0.0);
                    }    
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
                    // broken_solver = 0;
                    break;
                }
                else
                {
                    std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                    if(idx%100==0)
                    {
                        std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                        p_network->Write(output_file);
                    }
                }

                // If there is no convergence after all the iterations, print the error message.
                if(idx>=max_iter-1 && max_difference >= prev_max_difference)
                {
                    std::cout << "Problem encountered in " << str_directory_name << std::endl;
                    error_log << "\n Problem encountered in " << str_directory_name << std::endl; 
                    broken_solver = 1;
                    break;
                }                                
                prev_max_difference = max_difference;
            }

            // If simulation doesn't converge, move on to next layout and log problem 
            if (broken_solver == 1)
            {
                broken_layouts_file.open("/scratch/narain/testoutput/TestMetastaticNetworkTissue/broken_layouts.txt", std::ios_base::app);
                broken_layouts_file << str_directory_name << " \n"; 
                broken_layouts_file.close();

                // Move onto the next selection
                break;
            }

            // Set up the grid calculator
            std::shared_ptr<GridCalculator<2> > p_grid_calc = GridCalculator<2>::Create();
            p_grid_calc->SetGrid(p_grid);

            // Set the simulation dimensions just to save time
            const int DIM = 2;

            // Set up the cell nodes
            std::vector<Node<2>*> nodes;
            unsigned nodeNum=0;
            // unsigned numStromaNodes = 0;

            // // For staggered rows
            // unsigned rownum = 1;
            // double offset;
            // double node_spacing = 1.0;
            
            // // Mark the stroma nodes
            // for (double x=0; x<=dimless_domain_size_x; x=x+node_spacing)//x++)
            // {
            //     offset = 0;
            //     // if (rownum % 2)
            //     // {
            //     //     offset = 0;
            //     // }
            //     for (double y=0+offset; y<=dimless_domain_size_y; y=y+node_spacing)//y++)
            //     {
            //         if (DIM == 2) 
            //         {
            //             nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
            //             nodeNum++;
            //             numStromaNodes++;
            //         }
            //     }
            //     rownum++;
            // }

            // Create a node to mark the domain origin
            std::vector<unsigned> location_indices;
            nodes.push_back(new Node<DIM>(0, false, 0, 0));                    
            nodeNum++;
            // location_indices.push_back(0);

            // Read the tumour cell coordinates from a file
            unsigned numTumourNodes = 0;
            string line_1;
            std::ifstream dapi_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/edges/biological/lectin_monotile/dapi_data.txt");
            while (std::getline(dapi_file, line_1)) 
            {
                std::stringstream split(line_1);
                double x_coor, y_coor, diameter_calc;
                if (split >> x_coor >> y_coor >> diameter_calc) // Read x, y, and d values from the line
                {
                    nodes.push_back(new Node<DIM>(nodeNum, false, x_coor, y_coor));
                    nodeNum++;
                    numTumourNodes++;

                    // Optionally store d for later use
                    // d_values.push_back(diameter_calc);
                }
            }
            std::cout << "nodes.size() = " << nodes.size() << std::endl;

            // Pick the cell mutation and proliferation state
            std::vector<CellPtr> cells;
            MAKE_PTR(WildTypeCellMutationState, p_wild_state);
            MAKE_PTR(StemCellProliferativeType, p_stem_state);
            // MAKE_PTR(CancerCellMutationState, p_cancer_state);
            // MAKE_PTR(QuiescentCancerCellMutationState, p_quiescent_cancer_state);

            // Create the cancer and stroma labels
            unsigned cancerLabelColour = 9;
            unsigned stromaLabelColour = 1;
            MAKE_PTR_ARGS(CellLabel, p_stroma_label, (stromaLabelColour));
            MAKE_PTR_ARGS(CellLabel, p_cancer_label, (cancerLabelColour));

            // Print the cell types and their colours
            std::cout << "WildTypeCellMutationState has Legacy Cell type " <<  p_wild_state->GetColour() << std::endl;
            std::cout << "StemCellProliferativeType has Legacy Cell type " <<  p_stem_state->GetColour() << std::endl;
            std::cout << "p_stroma_label has Legacy Cell type " <<  p_stroma_label->GetColour() << std::endl;
            std::cout << "p_cancer_label has Legacy Cell type " <<  p_cancer_label->GetColour() << std::endl;
            // std::cout << "Cancer state has Legacy Cell type " <<  p_cancer_state->GetColour() << std::endl;

            // // Seed the stromal cells
            // for (unsigned int i=0; i<1; i++)
		    // {
            //     NoCellCycleModel* p_model = new NoCellCycleModel;
            //     p_model->SetDimension(2);
            //     CellPtr p_cell(new Cell(p_wild_state, p_model));
            //     // p_cell->SetCellProliferativeType(p_stem_type);
            //     // p_model->SetStemCellG1Duration(8.0);
            //     // p_model->SetTransitCellG1Duration(8.0);
            //     p_cell->AddCellProperty(p_stroma_label);
            //     // p_cell->SetBirthTime(birth_time);
            //     cells.push_back(p_cell);
            // }

            // Seed the tumour cells
		    for (unsigned int i=0; i<1+numTumourNodes; i++)  // iterate over the number of tumour nodes plus one to add a cell and location index for the origin
		    {
                SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
                p_model->SetDimension(2);
                CellPtr p_cell(new Cell(p_wild_state, p_model));
                p_cell->SetCellProliferativeType(p_stem_state);
                p_cell->AddCellProperty(p_cancer_label);
                double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                 (  p_model->GetStemCellG1Duration()
                                  + p_model->GetSG2MDuration() );
                p_cell->SetBirthTime(birth_time);
                cells.push_back(p_cell);
                location_indices.push_back(i);
            }
            std::cout << "cells.size = " << cells.size() << std::endl;
            std::cout << "location_indices.size = " << location_indices.size() << std::endl;

            // Set up a mesh for the cells 
            NodesOnlyMesh<2> mesh;                       
            mesh.ConstructNodesWithoutMesh(nodes, 1.5);                       
            NodeBasedCellPopulation<2> cell_population(mesh, cells, location_indices);                       
	    	cell_population.SetAbsoluteMovementThreshold(DBL_MAX);  //Set big movement threshold

            // for (std::list<CellPtr>::iterator cell_iter = this->cell_population.begin();
            //      cell_iter != this->mCells.end();
            //      )
            // {
            //     if ((*cell_iter)->IsDead())
            //     {
            //         // Get the location index corresponding to this cell
            //         unsigned location_index = this->GetLocationIndexUsingCell(*cell_iter);

            //         // Use this to remove the cell from the population
            //         RemoveCellUsingLocationIndex(location_index, (*cell_iter));

            //         // Erase cell and update counter
            //         cell_iter = this->mCells.erase(cell_iter);
            //         num_removed++;
            //     }
            //     else
            //     {
            //         ++cell_iter;
            //     }
            // }

            ////////////////////////////////////////////////////////////////

            // Set up the VEGF PDE as we set the O2 one (currently turned off)
            auto p_vegf_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
            p_vegf_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpVegfDiffusivity->GetValue("User"));
            p_vegf_pde->SetContinuumLinearInUTerm(-Owen11Parameters::mpVegfDecayRate->GetValue("User"));
        
            // Set up a map for different release rates depending on cell type. Also include a threshold intracellular VEGF below which there is no release (currently turned off)
            auto p_normal_and_quiescent_cell_source = CellStateDependentDiscreteSource<2>::Create();  // initialise the source (which combines normal and quiescent cells)
            std::map<unsigned, QConcentrationFlowRate > normal_and_quiescent_cell_rates;  // map to store secretion rates
            std::map<unsigned, QConcentration > normal_and_quiescent_cell_rate_thresholds;  // map to store secretion thresholds
            normal_and_quiescent_cell_rates[p_wild_state->GetColour()] = Owen11Parameters::mpCellVegfSecretionRate->GetValue("User");  // set the rate for normal cells
            normal_and_quiescent_cell_rate_thresholds[p_wild_state->GetColour()] = 0.27*unit::mole_per_metre_cubed;  // set the threshold for normal cells
            normal_and_quiescent_cell_rates[p_wild_state->GetColour()] = Owen11Parameters::mpCellVegfSecretionRate->GetValue("User");  // set the rate for quiescent cells
            normal_and_quiescent_cell_rate_thresholds[p_wild_state->GetColour()] = 0_M;  // set the threshold for quiescent cells
            p_normal_and_quiescent_cell_source->SetStateRateMap(normal_and_quiescent_cell_rates);  // configure the map
            p_normal_and_quiescent_cell_source->SetLabelName("VEGF");  // set the label for the source
            p_normal_and_quiescent_cell_source->SetStateRateThresholdMap(normal_and_quiescent_cell_rate_thresholds);

            // Add a vessel related VEGF sink (currently turned off)
            auto p_vessel_vegf_sink = VesselBasedDiscreteSource<2>::Create();
            p_vessel_vegf_sink->SetReferenceConcentration(0.0_M);
            p_vessel_vegf_sink->SetVesselPermeability(Owen11Parameters::mpVesselVegfPermeability->GetValue("User"));
            // p_vessel_vegf_sink->SetVesselPermeability(0.0);
            // p_vegf_pde->AddDiscreteSource(p_vessel_vegf_sink);

            // Set up a finite difference solver as before for the VEGF
            auto p_vegf_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
            p_vegf_solver->SetPde(p_vegf_pde);
            p_vegf_solver->SetLabel("VEGF_Extracellular");  // set the label for the external field
            p_vegf_solver->SetGrid(p_grid);

            // Set up the angiogenesis solver (runs on deactivated VEGF function but seems to be necessary for main solver)
            std::shared_ptr<AngiogenesisSolver<2> > p_angiogenesis_solver = AngiogenesisSolver<2>::Create();
            std::shared_ptr<Owen2011SproutingRule<2> > p_sprouting_rule = Owen2011SproutingRule<2>::Create();
            std::shared_ptr<Owen2011MigrationRule<2> > p_migration_rule = Owen2011MigrationRule<2>::Create();
            p_angiogenesis_solver->SetMigrationRule(p_migration_rule);
            p_angiogenesis_solver->SetSproutingRule(p_sprouting_rule);
            p_sprouting_rule->SetDiscreteContinuumSolver(p_vegf_solver);
            p_migration_rule->SetDiscreteContinuumSolver(p_vegf_solver);
            p_angiogenesis_solver->SetVesselGridCalculator(p_grid_calc);
            p_angiogenesis_solver->SetVesselNetwork(p_network);
            
            // Initialise the solver
            auto p_microvessel_solver = MicrovesselSolver<2>::Create();
            p_microvessel_solver->SetVesselNetwork(p_network);
            p_microvessel_solver->SetOutputFrequency(1_h);
            // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
            p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
            p_microvessel_solver->AddDiscreteContinuumSolver(p_vegf_solver);
            p_microvessel_solver->SetAngiogenesisSolver(p_angiogenesis_solver);

            // Specifies which extracellular fields to update based on PDE
            boost::shared_ptr<MicrovesselSimulationModifier<2> > p_microvessel_modifier =
            boost::shared_ptr<MicrovesselSimulationModifier<2> >(new MicrovesselSimulationModifier<2> ());
            p_microvessel_modifier->SetMicrovesselSolver(p_microvessel_solver);
            std::vector<std::string> update_labels;
            update_labels.push_back("oxygen");
            p_microvessel_modifier->SetCellDataUpdateLabels(update_labels);

            // Set up the off-lattice simulation
            OffLatticeSimulation<2> simulator(cell_population);
            simulator.AddSimulationModifier(p_microvessel_modifier);

            // Add the spring forces between cells
            // MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
            // simulator.AddForce(p_force);

            // Add a cell killer to remove apoptotic cells
            boost::shared_ptr<ApoptoticCellKiller<2>> p_apoptotic_cell_killer(new ApoptoticCellKiller<2>(&cell_population));
            simulator.AddCellKiller(p_apoptotic_cell_killer);

            // Add a radiotherapy cell killer
            boost::shared_ptr<LQRadiotherapyCellKiller<2>> p_rt_killer =
            boost::shared_ptr<LQRadiotherapyCellKiller<2>>(new LQRadiotherapyCellKiller<2>(&cell_population));

            // Set the RT parameters
            p_rt_killer->SetDoseInjected(2.0*unit::gray);  // modify dose here
            p_rt_killer->SetOerAlphaMax(1.75);
            p_rt_killer->SetOerAlphaMin(1.0);
            p_rt_killer->SetOerBetaMax(3.25);
            p_rt_killer->SetOerBetaMin(1.0);
            // p_rt_killer->SetOerConstant(oxygen_solubility_at_stp * 3.28_Pa);  // adjust the radiobiological threshold here
            p_rt_killer->SetOerConstant(0.004499758549541243*unit::mole_per_metre_cubed);  // K_OER (3.28 mmHg)
            p_rt_killer->SetAlphaMax(0.3 * unit::per_gray);
            p_rt_killer->SetBetaMax(0.03 * unit::per_gray_squared);
            p_rt_killer->UseOer(true);  // turn OER on and off here

            // Set up the dosage times
            std::vector<QTime > rt_times;
            rt_times.push_back(3600.0*1.0*unit::seconds);
            p_rt_killer->SetTimeOfRadiation(rt_times);
            // p_rt_killer->AddTimeOfRadiation(3600.0*96.0*unit::seconds);
            simulator.AddCellKiller(p_rt_killer);

            // Add another modifier for updating cell cycle quantities.
            // boost::shared_ptr<Owen2011TrackingModifier<2> > p_owen11_tracking_modifier(new Owen2011TrackingModifier<2>);
            // simulator.AddSimulationModifier(p_owen11_tracking_modifier);

            // Set up the simulation time and run it
            simulator.SetOutputDirectory(str_directory_name);
            simulator.SetSamplingTimestepMultiple(1);
            simulator.SetDt(1.0);
            simulator.SetEndTime(3.0);  // hours
            simulator.Solve();
            // std::cout << "This is line number: " << __LINE__ << std::endl;

            // Write the output network file (visualise with Paraview: set Filters->Alphabetical->Tube)
            std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
            p_network->Write(output_file);

            // Dump our parameter collection to an xml file and, importantly, clear it for the next test
            ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
            ParameterCollection::Instance()->Destroy();
            BaseUnits::Instance()->Destroy();
            SimulationTime::Instance()->Destroy();

            // If simulation doesn't converge, move on to next layout
            if (broken_solver == 1)
            {
                break;
            }
        }

        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

};

#endif /*TESTPAPER3_HPP_*/