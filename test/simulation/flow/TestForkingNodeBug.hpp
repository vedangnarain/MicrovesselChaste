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

This script contains a test that generates a forking network and prunes it stochastically. Currently, pruning results in patches of low O2 concentration around nodes that lead to no-flow vessels.

H = Haematocrit
BC = boundary condition
RT = Radiotherapy
PQ = Perfusion Quotient
Fig. = Figure
# = number

9/11/21

*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialisation
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef TESTFORKINGNODEBUG_HPP_
#define TESTFORKINGNODEBUG_HPP_
#define _BACKWARD_BACKWARD_WARNING_H 1  // Cut out the VTK deprecated warning

// General functionality
#include <boost/lexical_cast.hpp>
#include <cxxtest/TestSuite.h>
#include <fstream>
#include <iomanip>
#include <iostream>
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
class TestForkingNodeBug : public CxxTest::TestSuite
{

public:

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tests
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Make a multi-generation forking network with stochastic pruning
    void TestLeakyForkingNodeBug()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Set network in filename
        std::string network_name = "TestLeakyDichotomousNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("leaky_forking_stochastic_pruning_perfusion_quotients.txt", std::ios_base::app);
        outfile.close();

        // Run the simulation with different solvers of interest (currently only set to Constant and Pries)
        for (unsigned h_solver=1; h_solver<=2; h_solver++)
        {   
            // Generate the network for various lambdas
            for (unsigned k_aux=1; k_aux<2; k_aux++)  // generate the network for lambda = 4
            {
                // Set key vessel parameters
                double dimless_length = 1.0;  
                unsigned order = 6;  // choose the number of bifurcating generations (inlet and outlet vessels don't count)
                for(unsigned i_aux=1; i_aux<order+1; i_aux++)
                {
                    dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
                }
   
                // Set up the reference length for the simulation
                QLength reference_length(1.0_um);
                BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

                // Set input radius
                QLength input_radius(10.0 *GenericParameters::mpCapillaryRadius->GetValue());  // = 50_um

                // Set grid spacing
                QLength grid_spacing = 10.0_um;  // the simulation time gets quite long if you reduce the resolution further

                // Set the max. heterog. parameter (alpha = 1+(max_alpha*0.1), i.e., 5 is 1.5)
                unsigned max_alpha = 4;   

                // Set the viscosity
                QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;

                // Set the inlet and initial haematocrit
                double initial_haematocrit = 0.45;

                // Set threshold for perfusion quotient
		        QFlowRate threshold = 1.e-12*unit::metre_cubed_per_second;  // set flow threshold for what counts as a perfused vessel

                // Generate the networks
                VesselNetworkGenerator<2> network_generator;

                // Set lambda
                double lambda;  // lambda = length/diameter
                double twicelambda;  // used as an input parameter
                lambda = 2.0+double(k_aux)*2.0;
                twicelambda = 2.0*lambda;

                // Run the simulation with different heterogeneities (currently only set to 1.4)
                for (unsigned n_alpha=4; n_alpha<=max_alpha; n_alpha++)
                { 
                    // Set up flag for broken solver
                    unsigned broken_solver = 0;

                    // Set alpha
                    double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 

                    // Height of of first-order vessels
                    QLength main_vert_length = 0.9*twicelambda*input_radius*pow(2.0,-1.0/3.0);
                    
                    // Height of the domain
                    QLength domain_side_length_y = 4.0*main_vert_length;

                    // Length of the domain
                    QLength domain_side_length_x = dimless_length*2.0*twicelambda*input_radius;

                    // Set the number of trials and the maximum gamma value
                    int max_trials = 10;
                    int max_gamma = 15;
                    int gamma = 10;

                    // Iterate over gamma values (currently starts at 10 and ends at 15)
                    while (gamma <= max_gamma)
                    {
                        // Run the trials                    
                        int n_trial = 1;
                        while (n_trial <= max_trials)
                        {  
                            // Generate the network
                            std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, main_vert_length, input_radius, alpha, twicelambda);

                            // Identify input and output nodes and assign them properties
                            VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(0.0_um,2.0*main_vert_length));
                            VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(domain_side_length_x, 2.0*main_vert_length));
                            p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
                            p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
                            p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
                            p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);

                            // Set the h-solver
                            std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                            std::string solver_name;
                            if (h_solver==1)
                            {
                                solver_name = "ConstantHaematocrit/"; 
                                std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                                auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                                p_haematocrit_solver->SetVesselNetwork(p_network);
                                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                                p_abstract_haematocrit_solver = p_haematocrit_solver;     
                            }
                            else if (h_solver==2)
                            {
                                solver_name = "PriesHaematocrit/"; 
                                std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                                auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                                p_haematocrit_solver->SetVesselNetwork(p_network);
                                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                                p_abstract_haematocrit_solver = p_haematocrit_solver;                
                            }
                            else if (h_solver==3)
                            {
                                solver_name = "MemoryHaematocrit/"; 
                                std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                                auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                                p_haematocrit_solver->SetVesselNetwork(p_network);
                                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                                p_abstract_haematocrit_solver = p_haematocrit_solver;     
                            }
                            else if (h_solver==4)
                            {
                                solver_name = "FungHaematocrit/"; 
                                std::cout << "Now using FungHaematocritSolver..." << std::endl;
                                auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                                p_haematocrit_solver->SetVesselNetwork(p_network);
                                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                                p_abstract_haematocrit_solver = p_haematocrit_solver;     
                            }
                                            
                            // Store the lambda values for file name
                            std::stringstream lambda_stream;
                            lambda_stream << std::fixed << std::setprecision(0) << lambda;
                            std::string lambda_string = lambda_stream.str();
                            std::string lambda_directory = network_name + solver_name + "/Lambda" + lambda_string; 

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
   
                            // Conduct stochastic pruning (thinner vessels are more likely to be pruned)
                            std::vector<std::shared_ptr<Vessel<2> > > vessels_to_remove;
                            std::vector<std::shared_ptr<Vessel<2> > > vessels = p_network->GetVessels();
                            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
                            for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                            {  
                                // srand((unsigned) time(NULL));  // seed random number
                                float random_number = (float) rand() / RAND_MAX;    
                                float vessel_radius = (*vessel_iterator)->GetRadius();
                                float prob_death = exp(-(vessel_radius*1000000)/gamma);  
                                // std::cout << "vessel radius " << vessel_radius*1000000 << " prob_death " << prob_death << " random number " << random_number << std::endl;
                                if(random_number <= prob_death)
                                {
                                    if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() != 1 && (*vessel_iterator)->GetEndNode()->GetNumberOfSegments() != 1)
                                    {
                                        vessels_to_remove.push_back((*vessel_iterator));
                                    }
                                }
                            }
                            for(unsigned idx=0; idx<vessels_to_remove.size(); idx++)
                            {
                                p_network->RemoveVessel(vessels_to_remove[idx],true);
                            }

                            // Set filename
                            std::stringstream alpha_stream;
                            std::stringstream trial_stream;
                            std::stringstream gamma_stream;
                            alpha_stream << std::fixed << std::setprecision(2) << alpha;
                            trial_stream << std::fixed << std::setprecision(0) << n_trial;
                            gamma_stream << std::fixed << std::setprecision(0) << gamma;
                            std::string alpha_string = alpha_stream.str();
                            std::string trial_string = trial_stream.str();
                            std::string gamma_string = gamma_stream.str();
                            std::string file_name = lambda_directory + "/Alpha" + alpha_string + "/Gamma" + gamma_string + "/Trial" + trial_string;
                            std::string str_directory_name = file_name;
                            auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                            // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance. The loop runs the solvers, calculates the max. difference in the network between the previous H of a vessel and the new H in the vessel, and repeats this until all vessels converge on a value within a certain tolerance.
                            unsigned max_iter = 100;  // simulation fails if it doesn't converge in these many iterations
                            double tolerance2 = 0.01; //1.e-10
                            std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                            std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                            for (unsigned idx=0;idx<max_iter;idx++)
                            {
                                // Run the solvers
                                p_impedance_calculator->Calculate();
                                flow_solver.SetUp();
                                flow_solver.Solve();
                                p_abstract_haematocrit_solver->Calculate();
                                p_viscosity_calculator->Calculate();

                                // Check for convergence 
                                double max_difference = 0.0;
                                double h_for_max = 0.0;
                                double prev_for_max = 0.0;
                                for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)  // for all the segments 
                                {
                                    // Set segments with no flow to have no haematocrit
                                    if (fabs(segments[segment_index]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                    {
                                        segments[segment_index]->GetFlowProperties()->SetHaematocrit(0.0);
                                    }                                    
                                    double current_haematocrit = segments[segment_index]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                    double difference = std::abs(current_haematocrit - previous_haematocrit[segment_index]);  // difference in new H in segment and old H
                                    if(difference>max_difference)  // if the difference is greater than previous max.
                                    {
                                        max_difference = difference;
                                        h_for_max = current_haematocrit;
                                        prev_for_max = previous_haematocrit[segment_index];
                                    }
                                    previous_haematocrit[segment_index] = current_haematocrit;
                                }
                                std::cout << "Segment H. at max difference: " << h_for_max << ", Prev. segment H at max difference:" << prev_for_max << std::endl;
                                if(max_difference<=tolerance2)  
                                {
                                    std::cout << "Converged after: " << idx << " iterations." <<  std::endl;
                                    broken_solver = 0;
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
                                if(idx==max_iter-1)  // if there is no convergence after all the iterations, print the error message
                                {
                                    std::cout << "Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using n_alpha = " << n_alpha << " and lambda = " << lambda << " in trial = " << n_trial << std::endl;
                                    // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                                    error_log << "\n Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using n_alpha = " << n_alpha << " and lambda = " << lambda << " in trial = " << n_trial; 
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

                            // Write the PQs
                            outfile.open("leaky_forking_stochastic_pruning_perfusion_quotients.txt", std::ios_base::app);
                            outfile << network_name << " " << solver_name << " " << lambda_string << " " << alpha_string << " " << gamma_string << " " << trial_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                            outfile.close();
                            
                            // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                            std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                            p_network->Write(output_file);
                            
                            // Move onto next trial
                            n_trial = n_trial + 1;  

                            // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                            ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                            ParameterCollection::Instance()->Destroy();
                            BaseUnits::Instance()->Destroy();
                            SimulationTime::Instance()->Destroy();
                        }
                        // Move on to next gamma
                        gamma = gamma + 1;
                    }
                }
            }
        }
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }
};

#endif /*TESTFORKINGNODEBUG_HPP_*/
