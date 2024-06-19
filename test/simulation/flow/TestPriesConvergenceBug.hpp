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

This script contains two tests to generate a honeycomb neighbourhood with the diameters assigned (a) homogenously (b) heterogeneously.

H = Haematocrit
BC = boundary condition
RT = Radiotherapy
PQ = Perfusion Quotient
Fig. = Figure
# = number

28/5/24

*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialisation
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef TESTPRIESCONVERGENCEBUG_HPP_
#define TESTPRIESCONVERGENCEBUG_HPP_
#define _BACKWARD_BACKWARD_WARNING_H 1  // Cut out the VTK deprecated warning

// Essential functionality
#include <boost/lexical_cast.hpp>
#include <cxxtest/TestSuite.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>  // needed for exp function
#include <sstream>
#include <string>
#include <time.h>
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
#include "BetteridgeHaematocritSolver.hpp"
#include "ConstantHaematocritSolver.hpp"
#include "YangHaematocritSolver.hpp"
#include "PriesHaematocritSolver.hpp"
#include "PriesWithMemoryHaematocritSolver.hpp"

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
using namespace std;

// Make a test class
class TestPriesConvergenceBug : public CxxTest::TestSuite
{

public:
 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Test
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define grid size for PDE
    QLength grid_spacing = 10.0_um;  // the simulation time gets quite long if you reduce the resolution further

    // Make a 2D hexagonal neighbourhood on a PDE grid with a homogeneous diameter assignment
    void TestHomogeneousHexagonalNetwork()
    {
        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Set file name
        std::string network_name = "TestDebuggingNetworks/HexagonalNeighbourhood/HomogeneousDiameters/NoPruning/";   

        // Define the key parameters
        QLength large_vessel_radius = 10_um;
        double dimless_domain_size_x = 700.0;  // x-coordinate of output node
        double dimless_domain_size_y = 606.218 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        unsigned dimless_vessel_length = 100.0;
        std::vector<std::vector<unsigned> > Order;
        std::shared_ptr<VesselNetwork<2> > p_network;
        std::vector<std::shared_ptr<Vessel<2> > > vessels;

        // Initialise the simulation space
        VesselPtr<2> p_current_vessel;
        QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
        QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
        std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
        std::ofstream outfileMean;

        // Define the key parameters
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Read the network from a file
        VesselNetworkGenerator<2> network_generator;
        std::ifstream in("../Chaste/projects/MicrovesselChaste/test/simulation/flow/" + to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_Neighbourhood.txt");
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

        // Run the simulation with different solvers
        std::string solver_name;
        for (unsigned h_solver=2; h_solver<=2; h_solver++)
        {
            // Set the h-solver
            if (h_solver==1)
            {
                solver_name = "ConstantHaematocrit";
            }
            else if (h_solver==2)
            {
                solver_name = "PriesHaematocrit";
            }
            else if (h_solver==3)
            {
                solver_name = "MemoryHaematocrit";
            }
            else if (h_solver==4)
            {
                solver_name = "FungHaematocrit";
            }

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
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;                    
            }
            else if (h_solver==3)
            {
                std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;      
            }
            else if (h_solver==4)
            {
                std::cout << "Now using FungHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
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
            flow_solver.SetUseDirectSolver(true);

            // Set filename                  
            std::string file_name = network_name + solver_name;                        
            std::string str_directory_name = file_name;
            std::cout << str_directory_name << std::endl;
            auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

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

        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a 2D hexagonal neighbourhood on a PDE grid with heterogeneous diameter assignments
    void xTestHeterogeneousHexagonalNetwork()
    {
        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);
        
        // Set the number of diameter distribution variations
        unsigned max_n_sigma = 2;
        unsigned max_n_mu = 2;

        // Set file name
        std::string network_name = "TestDebuggingNetworks/HexagonalNeighbourhood/HeterogeneousDiameters/NoPruning/";   

        // Define the key parameters
        QLength large_vessel_radius = 10_um;
        double dimless_domain_size_x = 700.0;  // x-coordinate of output node
        double dimless_domain_size_y = 606.218 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        unsigned dimless_vessel_length = 100.0;
        std::vector<std::vector<unsigned> > Order;
        std::shared_ptr<VesselNetwork<2> > p_network;
        std::vector<std::shared_ptr<Vessel<2> > > vessels;
        unsigned n_vessels = 42;  // number of vessels

        // Initialise the simulation space
        VesselPtr<2> p_current_vessel;
        QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
        QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
        std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
        std::ofstream outfileMean;

        // Define the key parameters
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Read the network from a file
        VesselNetworkGenerator<2> network_generator;
        std::ifstream in("../Chaste/projects/MicrovesselChaste/test/simulation/flow/" + to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_Neighbourhood.txt");
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

        // Run the simulation with different solvers
        std::string solver_name;
        for (unsigned h_solver=2; h_solver<=2; h_solver++)
        {
            // Loop over different diameter layouts
            for(unsigned list_number=1; list_number<=100; list_number++)
            {                   
                // Loop over increasing level of heteregeneity 
                for(unsigned n_sigma=0; n_sigma<=max_n_sigma; n_sigma++)
                {
                    // Loop over increasing level of mean diameter 
                    for(unsigned n_mu=0; n_mu<=max_n_mu; n_mu++)
                    {
                        // Convert to actual mu
                        std::vector<std::string> mus{"22.76", "28.5", "33.64"};
                        std::string mu = mus[n_mu];  // mu determines the mean

                        // Convert to actual sigma
                        std::vector<std::string> sigmas{"8.68", "13.23", "17.49"};
                        std::string sigma = sigmas[n_sigma];  // sigma determines the SD

                        // Choose a radii list file and read the file to an array
                        std::vector<std::vector<double> > radii_array;
                        radii_array.clear();
                        string line_1;
                        std::ifstream radius_list_file("../Chaste/projects/MicrovesselChaste/test/simulation/flow/" + to_string(dimless_vessel_length)+"VesselLength/descending_hexagonal_diameter_log_normal_distribution_sigma_"+sigma+"/mu_"+mu+"/radii_list_"+to_string(list_number)+".txt");
                        while (std::getline(radius_list_file, line_1)) 
                        {
                            radii_array.push_back(std::vector<double>());
                            std::stringstream split(line_1);
                            double value;
                            while (split >> value)
                            {
                                radii_array.back().push_back(value);
                            }
                        }     

                        // Set the h-solver
                        if (h_solver==1)
                        {
                            solver_name = "ConstantHaematocrit";
                        }
                        else if (h_solver==2)
                        {
                            solver_name = "PriesHaematocrit";
                        }
                        else if (h_solver==3)
                        {
                            solver_name = "MemoryHaematocrit";
                        }
                        else if (h_solver==4)
                        {
                            solver_name = "FungHaematocrit";
                        }

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
                        p_segment->SetRadius(large_vessel_radius);
                        VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   

                         // Assign vessel radii based on imported arrays
                        unsigned int vessel_index = 0;
                        while(vessel_index < n_vessels)
                        {
                            p_current_vessel = p_network->GetVessel(vessel_index);
                            // std::cout << "../Chaste/projects/MicrovesselChaste/test/simulation/flow/" << to_string(dimless_vessel_length)<<"VesselLength/descending_hexagonal_diameter_log_normal_distribution_sigma_"<<sigma<<"/mu_"<<mu<<"/radii_list_"<<to_string(list_number)<<".txt" << std::endl;

                            // Get segment radius from array
                            QLength new_vessel_radius(double(radii_array[vessel_index][0])*unit::microns);

                            p_current_vessel->SetRadius(new_vessel_radius);
                            vessel_index++;
                        }
                        vessels = p_network->GetVessels();

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
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;                    
                        }
                        else if (h_solver==3)
                        {
                            std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;      
                        }
                        else if (h_solver==4)
                        {
                            std::cout << "Now using FungHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
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
                        flow_solver.SetUseDirectSolver(true);

                        // Set filename
                        std::stringstream selection_stream;
                        std::stringstream sd_stream;
                        std::stringstream mean_stream;
                        selection_stream << std::fixed << std::setprecision(0) << to_string(list_number);                    
                        sd_stream << std::fixed << std::setprecision(0) << sigma;                     
                        mean_stream << std::fixed << std::setprecision(0) << mu;                          
                        std::string selection_string = selection_stream.str();
                        std::string sd_string = sd_stream.str();                    
                        std::string mean_string = mean_stream.str();                   
                        std::string file_name = network_name + solver_name + "/Selection" + selection_string + "/Sigma" + sd_string + "/Mu" + mean_string;                        
                        std::string str_directory_name = file_name;
                        std::cout << str_directory_name << std::endl;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

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

        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Pick a diameter distribution and replace one diameter at a time with the average value to isolate the errant diameter
    void xTestChangeOneDiameterPerSimulation()
    {
        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Set file name
        std::string network_name = "TestDebuggingNetworks/HexagonalNeighbourhood/DiameterIsolation/NoPruning/";   

        // Define the key parameters
        QLength large_vessel_radius = 10_um;
        double dimless_domain_size_x = 700.0;  // x-coordinate of output node
        double dimless_domain_size_y = 606.218 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        unsigned dimless_vessel_length = 100.0;
        std::vector<std::vector<unsigned> > Order;
        std::shared_ptr<VesselNetwork<2> > p_network;
        std::vector<std::shared_ptr<Vessel<2> > > vessels;
        unsigned n_vessels = 42;  // number of vessels

        // Initialise the simulation space
        VesselPtr<2> p_current_vessel;
        QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
        QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
        std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
        std::ofstream outfileMean;

        // Define the key parameters
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Read the network from a file
        VesselNetworkGenerator<2> network_generator;
        std::ifstream in("../Chaste/projects/MicrovesselChaste/test/simulation/flow/" + to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_Neighbourhood.txt");
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

        // Run the simulation with different solvers
        std::string solver_name;
        for (unsigned h_solver=2; h_solver<=2; h_solver++)
        {
            // Loop over different diameter layouts
            for(unsigned list_number=0; list_number<=n_vessels; list_number++)
            {                   
                // Choose a radii list file and read the file to an array
                std::vector<std::vector<double> > radii_array;
                radii_array.clear();
                string line_1;
                // std::ifstream radius_list_file("../Chaste/projects/MicrovesselChaste/test/simulation/flow/" + to_string(dimless_vessel_length)+"VesselLength/debugging/Selection1Sigma0Mu0/diameter_list_"+to_string(list_number)+".txt");
                std::ifstream radius_list_file("../Chaste/projects/MicrovesselChaste/test/simulation/flow/" + to_string(dimless_vessel_length)+"VesselLength/debugging/diameter_list_"+to_string(list_number)+".txt");
                while (std::getline(radius_list_file, line_1)) 
                {
                    radii_array.push_back(std::vector<double>());
                    std::stringstream split(line_1);
                    double value;
                    while (split >> value)
                    {
                        radii_array.back().push_back(value);
                    } 
                }     

                // Set the h-solver
                if (h_solver==1)
                {
                    solver_name = "ConstantHaematocrit";
                }
                else if (h_solver==2)
                {
                    solver_name = "PriesHaematocrit";
                }
                else if (h_solver==3)
                {
                    solver_name = "MemoryHaematocrit";
                }
                else if (h_solver==4)
                {
                    solver_name = "FungHaematocrit";
                }

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
                p_segment->SetRadius(large_vessel_radius);
                VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   

                    // Assign vessel radii based on imported arrays
                unsigned int vessel_index = 0;
                while(vessel_index < n_vessels)
                {
                    p_current_vessel = p_network->GetVessel(vessel_index);

                    // Get segment radius from array
                    QLength new_vessel_radius(double(radii_array[vessel_index][0])*unit::microns);

                    p_current_vessel->SetRadius(new_vessel_radius);
                    vessel_index++;
                }
                vessels = p_network->GetVessels();

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
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;                    
                }
                else if (h_solver==3)
                {
                    std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;      
                }
                else if (h_solver==4)
                {
                    std::cout << "Now using FungHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
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
                flow_solver.SetUseDirectSolver(true);

                // Set filename
                std::stringstream selection_stream;
                selection_stream << std::fixed << std::setprecision(0) << to_string(list_number);                                         
                std::string selection_string = selection_stream.str();
                std::string file_name = network_name + solver_name + "/Selection" + selection_string;                        
                std::string str_directory_name = file_name;
                std::cout << str_directory_name << std::endl;
                auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

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

        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

};

#endif /*TESTPRIESCONVERGENCEBUG_HPP_*/
