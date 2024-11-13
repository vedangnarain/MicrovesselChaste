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
// Notes
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*

This script contains tests that are used to simulate blood flow in different vascular architectures. 

Outputs can usually be found in /tmp/<home_directory>/testoutput. The pointwise data for the field can be obtained by opening spreadsheet view in Paraview.

H = Haematocrit
BC = boundary condition
RT = Radiotherapy
PQ = Perfusion Quotient
Fig. = Figure
# = number

*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialisation
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef TESTENHANCEDPERFUSIONFOLLOWINGEXPOSURETORADIOTHERAPY_HPP_
#define TESTENHANCEDPERFUSIONFOLLOWINGEXPOSURETORADIOTHERAPY_HPP_
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
// #include "AlarconHaematocritSolver.hpp"
#include "BetteridgeHaematocritSolver.hpp"
#include "ConstantHaematocritSolver.hpp"
// #include "LinnengarHaematocritSolver.hpp"
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
class TestEnhancedPerfusionFollowingExposureToRadiotherapy : public CxxTest::TestSuite
{

public:

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Forking Network
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define key parameters here, i.e., things that are constant across architectures
    // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));  // = 40 um
    QLength grid_spacing = 10.0_um;  // the simulation time gets quite long if you reduce the resolution further
    // Note: 10.0_um is a rather low resolution for the tiny equilateral network.
    unsigned max_alpha = 3;  // set the max. heterog. parameter (alpha = 1+(max_alpha*0.1), i.e., 5 is 1.5) 

    // Choose the number of bifurcating generations (inlet and outlet vessels don't count)
    unsigned order = 6;

    // Make a multi-generation forking network with different h-splitting rules and individual pruning
    void TestDichotomousNetworkWithIndividualPruningAndFlow2DAndVaryingMeansPaper1()
        {
            // Initialise error log
            std::ostringstream error_log;
            error_log << "\n The following simulations failed to converge: \n";

            // Set network in filename
            std::string network_name = "TestDichotomousNetwork/";

            // Create the output file for the PQs
            std::ofstream outfile;
            outfile.open("/tmp/narain/testoutput/TestDichotomousNetwork/forking_nonrandom_individual_pruning_perfusion_quotients.txt", std::ios_base::app);
            outfile.close();

            // Define the key pruning parameters
            unsigned n_vessels = 250;  // number of non-inlet/outlet vessels from which to select ones to kill
            // unsigned n_vessels = 508;  // number of non-inlet/outlet vessels from which to select ones to kill
            // unsigned n_vessels = 124;  // number of non-inlet/outlet vessels from which to select ones to kill
            // double percToKill = 0.2;  // percentage of vessels to kill
            // unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
            unsigned ToBeKilled = n_vessels;  // number to kill

            // Run the simulation with different solvers of interest
            for (unsigned h_solver=1; h_solver<=1; h_solver++)
            {   
                // Generate the network for various lambdas
                for (unsigned k_aux=1; k_aux<2; k_aux++)  // generate the network for lambda = 4
                {
                    // Set key vessel parameters
                    double dimless_length = 1.0;  

                    // Non-dimensionalising the length
                    for(unsigned i_aux=1; i_aux<order+1; i_aux++)
                    {
                        dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
                    }
    
                    // Set up the reference length for the simulation
                    QLength reference_length(1.0_um);
                    BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

                    // Set input radius
                    QLength input_radius(7.5 *GenericParameters::mpCapillaryRadius->GetValue());  // = *5_um (*10_um for diameter)
                    // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

                    // Set the viscosity
                    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
                    // QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();  // = 0.0012

                    // Set the inlet and initial haematocrit
                    double initial_haematocrit = 0.45;
                    // double initial_haematocrit = 0.36;
                    // double initial_haematocrit = 0.54;

                    // Set threshold for perfusion quotient
                    QFlowRate threshold = 3.e-12*unit::metre_cubed_per_second;  // set flow threshold for what counts as a perfused vessel

                    // Generate the networks
                    VesselNetworkGenerator<2> network_generator;

                    // Set lambda
                    double lambda;  // lambda = length/diameter
                    double twicelambda;  // used as an input parameter
                    lambda = 2.0+double(k_aux)*2.0;
                    twicelambda = 2.0*lambda;

                    // Set the number of different means
                    unsigned max_means = 3;

                    // Run the simulation with different heterogeneities
                    for (unsigned n_alpha=0; n_alpha<=max_alpha; n_alpha++)
                    { 
                        // Offset the radii depending on the mean and SD
                        for (unsigned n_mean=0; n_mean<max_means; n_mean++)
                        { 
                            // Set up flag for broken solver
                            unsigned broken_solver = 0;

                            // Set alpha
                            double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 

                            // Height of of first-order vessels
                            QLength main_vert_length = 0.9*twicelambda*input_radius*pow(2.0,-1.0/3.0);
                            
                            // Height of the domain
                            // QLength domain_side_length_y = 4.0*main_vert_length;

                            // Length of the domain
                            QLength domain_side_length_x = dimless_length*2.0*twicelambda*input_radius;
                    
                            // Generate the network
                            std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, main_vert_length, input_radius, alpha, twicelambda);

                            // Identify input and output nodes and assign them properties
                            VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(0.0_um,2.0*main_vert_length));
                            VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(domain_side_length_x, 2.0*main_vert_length));
                            p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
                            // p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));  // = 3333.06
                            p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
                            p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
                            // p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));  // = 1999.84
                            p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);

                            // Specify the minimum offsets for all alphas (100 um, 6 gens)
                            // QLength min_array[] = {-9.95_um, -9.61_um, -8.75_um, -7.56_um, -6.21_um};
                            // QLength mean_array[] = {-4.21_um, -3.87_um, -3.01_um, -1.82_um, -0.47_um};
                            // QLength max_array[] = {0.93_um, 1.27_um, 2.13_um, 3.32_um, 4.67_um};

                            // Specify the minimum offsets for all alphas (75 um, 6 gens)
                            QLength min_array[] = {-1.77_um, -1.52_um, -0.87_um, 0.02_um, 1.04_um};
                            QLength mean_array[] = {3.97_um, 4.22_um, 4.87_um, 5.76_um, 6.78_um};
                            QLength max_array[] = {9.11_um, 9.36_um, 10.01_um, 10.9_um, 11.92_um};

                            // Specify the minimum offsets for all alphas (50 um, 6 gens)
                            // QLength min_array[] = {6.41_um, 6.58_um, 7.01_um, 7.6_um, 8.28_um};
                            // QLength mean_array[] = {12.15_um, 12.32_um, 12.75_um, 13.34_um, 14.02_um};
                            // QLength max_array[] = {17.29_um, 17.46_um, 17.89_um, 18.48_um, 19.16_um};
                        
                            // Specify the minimum offsets for all alphas (75 um, 5 gens)
                            // QLength min_array[] = {-7.64_um, -7.39_um, -6.74_um, -5.84_um, -4.8_um};
                            // QLength mean_array[] = {-1.9_um, -1.65_um, -1.0_um, -0.1_um, 0.94_um};
                            // QLength max_array[] = {3.24_um, 3.49_um, 4.14_um, 5.04_um, 6.08_um};

                            // Specify the minimum offsets for all alphas (50 um, 5 gens)
                            // QLength min_array[] = {2.49_um, 2.66_um, 3.09_um, 3.69_um, 4.38_um};
                            // QLength mean_array[] = {8.23_um, 8.4_um, 8.83_um, 9.43_um, 10.12_um};
                            // QLength max_array[] = {13.37_um, 13.54_um, 13.97_um, 14.57_um, 15.26_um};

                            // Initialise the offsets
                            std::vector<QLength> min_list(min_array, min_array+5);
                            std::vector<QLength> mean_list(mean_array, mean_array+5);
                            std::vector<QLength> max_list(max_array, max_array+5);

                            // Set the radius offset based on the alpha
                            std::vector<std::shared_ptr<Vessel<2> > > vessels;
                            vessels = p_network->GetVessels();
                            if (n_mean==0)
                                {
                                    QLength radius_offset = min_list[n_alpha];
                                    for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)
                                    {                            
                                        // Get the current segment's radius
                                        QLength current_radius = vessels[vessel_index]->GetRadius();

                                        // Calculate the new radius
                                        QLength new_vessel_radius = current_radius + (radius_offset*0.5);

                                        // Set the new radius
                                        vessels[vessel_index]->SetRadius(new_vessel_radius);
                                    }
                                }
                            else if (n_mean==1)
                            {
                                QLength radius_offset = mean_list[n_alpha];
                                for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)
                                {                            
                                    // Get the current segment's radius
                                    QLength current_radius = vessels[vessel_index]->GetRadius();

                                    // Calculate the new radius
                                    QLength new_vessel_radius = current_radius + (radius_offset*0.5);
    
                                    // Set the new radius
                                    vessels[vessel_index]->SetRadius(new_vessel_radius);
                                }
                            }
                            else if (n_mean==2)
                            {
                                QLength radius_offset = max_list[n_alpha];
                                for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)
                                {                            
                                    // Get the current segment's radius
                                    QLength current_radius = vessels[vessel_index]->GetRadius();

                                    // Calculate the new radius
                                    QLength new_vessel_radius = current_radius + (radius_offset*0.5);

                                    // Set the new radius
                                    vessels[vessel_index]->SetRadius(new_vessel_radius);
                                }
                            }
                            vessels = p_network->GetVessels();

                            // Set the h-solver
                            std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                            std::string solver_name;
                            if (h_solver==1)
                            {
                                solver_name = "ConstantHaematocrit/"; 
                                std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                                auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                                p_haematocrit_solver->SetVesselNetwork(p_network);
                                p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
                                p_abstract_haematocrit_solver = p_haematocrit_solver;     
                            }
                                                    
                            // Store the lambda values for file name
                            std::stringstream lambda_stream;
                            lambda_stream << std::fixed << std::setprecision(0) << lambda;
                            std::string lambda_string = lambda_stream.str();
                            std::string lambda_directory = network_name + solver_name + "/Lambda" + lambda_string; 

                            // Store the n_mean values for file name
                            std::stringstream mean_stream;
                            mean_stream << std::fixed << std::setprecision(0) << n_mean;
                            std::string mean_string = mean_stream.str();

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
        
                            // Prune all vessels up to specified dose 
                            for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
                            { 
                                // Display status message
                                std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                                // Set filename
                                std::stringstream alpha_stream;
                                std::stringstream kill_stream;
                                alpha_stream << std::fixed << std::setprecision(2) << alpha;
                                kill_stream << std::fixed << std::setprecision(0) << KilledVessels;
                                std::string alpha_string = alpha_stream.str();                            
                                std::string kill_string = kill_stream.str();
                                std::string file_name = lambda_directory + "/Alpha" + alpha_string + "/Mean" + mean_string + "/Kills" + kill_string;
                                std::string str_directory_name = file_name;
                                auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                                // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance. The loop runs the solvers, calculates the max. difference in the network between the previous H of a vessel and the new H in the vessel, and repeats this until all vessels converge on a value within a certain tolerance.
                                std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                                // Run the solvers
                                p_impedance_calculator->Calculate();
                                flow_solver.SetUp();
                                flow_solver.Solve();
                                p_abstract_haematocrit_solver->Calculate();
                                p_viscosity_calculator->Calculate();

                                // Check for convergence 
                                for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)  // for all the segments 
                                {
                                    // Set segments with no flow to have no haematocrit
                                    if (fabs(segments[segment_index]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                    {
                                        segments[segment_index]->GetFlowProperties()->SetHaematocrit(0.0);
                                    }                                    
                                }

                                // If solver doesn't converge, move on to next one
                                if (broken_solver == 1)
                                {
                                    continue;
                                }

                                // Write the PQs
                                outfile.open("/tmp/narain/testoutput/TestDichotomousNetwork/forking_nonrandom_individual_pruning_perfusion_quotients.txt", std::ios_base::app);
                                outfile << network_name << " " << solver_name << " " << lambda_string << " " << alpha_string << " " << mean_string << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                                outfile.close();
                                
                                // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                                std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                                p_network->Write(output_file);

                                // Remove the smallest vessel
                                QLength minimum_radius = input_radius;
                                unsigned int minimum_index = 0;
                                std::vector<std::shared_ptr<Vessel<2> > > vessels = p_network->GetVessels();
                                for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                                {
                                    // Get the current segment's radius
                                    QLength current_radius = vessels[vessel_index]->GetRadius();

                                    // If the current radius is less than the minimum radius, record the new minimum
                                    if (current_radius < minimum_radius)
                                    {
                                        minimum_radius = current_radius;
                                        minimum_index = vessel_index;
                                    }
                                }
                                p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel

                                // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                                ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                                ParameterCollection::Instance()->Destroy();
                                BaseUnits::Instance()->Destroy();
                                SimulationTime::Instance()->Destroy();
                            }
                        }
                    }
                }
            }
            // Print the error log
            std::string error_message = error_log.str();
            std::cout << error_message << std::endl; 
        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hexagonal Network 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Choose heterogeneity parameters
	unsigned thin_selections = 100;  // number of layouts from which to choose thin vessels (upto 100)
    QLength inlet_vessel_radius = 37.5_um;  // the maximum radius of the distribution
    // double dimless_domain_size_x = 2000.0; 
    // double dimless_domain_size_y = 2000.0 - 86.6025; 
    double dimless_domain_size_x = 2050.0;  // x-coordinate of output node
    double dimless_domain_size_y = 1818.65 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
    unsigned dimless_vessel_length = 100.0;
    unsigned max_n_alpha = 2;  // alpha controls the mean

    // Make a full 2D hexagonal network with flow, H-splitting, non-inlet/outlet radii set according to a log normal distribution based on biological networks, and kills-based pruning on the non-inlet/outlet vessels (correcting for the unexplained -50 um offset in the other simulation).
    void xTestConstantOffsetBiologicalHexagonalNeighbourhoodWithIndividualPruning2DPaper1()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";
        
        // Set network in filename
        std::string network_name = "TestHexagonalNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("/tmp/narain/testoutput/TestHexagonalNetwork/hex_lognormal_individual_perfusion_quotients.txt");
        outfile.close();

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        // double dimless_domain_size_x = 2050.0;  // x-coordinate of output node
        // double dimless_domain_size_y = 1818.65 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        // unsigned dimless_vessel_length = 100.0;
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        // double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        // unsigned n_vessels = 407;  // number of vessels from which to select ones to make thin
        unsigned n_vessels = 386;  // number of non-inlet/outlet vessels from which to select ones to make thin
        // double percToKill = 0.2;  // percentage of vessels to kill
        // unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
        unsigned ToBeKilled = 200;  // number to kill
        // unsigned ToBeKilled = 2;  // number to kill

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {
            // Initialise the simulation
            std::shared_ptr<VesselNetwork<2> > p_network;
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            VesselPtr<2> p_current_vessel;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            std::vector<std::vector<double> > QuotientMatrixMean(max_n_alpha,std::vector<double>(ToBeKilled+1));
            std::vector<std::vector<double> > QuotientMatrixAggregate(max_n_alpha,std::vector<double>(ToBeKilled+1,0.0));
            std::vector<std::vector<unsigned> > Order;
            double PerfusionQuotient;

            // Read the network layout from a file and store it in an array
            VesselNetworkGenerator<2> network_generator;
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_Offset.txt");
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

            // Loop over increasing level of heteregeneity (% of thin vessels)
            for(unsigned n_alpha=0; n_alpha<=max_n_alpha; n_alpha++)
            {
                // Convert to actual mu
                std::vector<std::string> alphas{"22.76", "28.5", "33.64"};
                std::string alpha = alphas[n_alpha];  // alpha determines the mean

                // Loop over different thin vessel layouts
                for(unsigned list_number=1; list_number<=thin_selections; list_number++)
                {                
                    // Choose a radii list file and read the file to an array (change this based on SD)
                    std::vector<std::vector<double> > radii_array;
                    radii_array.clear();
                    // int list_index=0;
                    string line_1;
                    std::ifstream radius_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_diameter_log_normal_distribution_sigma_8.68/mu_"+alpha+"/radii_list_"+to_string(list_number)+".txt");
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

                    // Choose the corresponding vessel ID order file and read the file to an array (change this based on SD
                    std::vector<std::vector<double> > id_array;
                    id_array.clear();
                    // int list_index=0;
                    string line_2;
                    std::ifstream id_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_diameter_log_normal_distribution_sigma_8.68/mu_"+alpha+"/id_list_"+to_string(list_number)+".txt");
                    while (std::getline(id_list_file, line_2)) 
                    {
                        id_array.push_back(std::vector<double>());
                        std::stringstream split(line_2);
                        double value;
                        while (split >> value)
                        {
                            id_array.back().push_back(value);
                        }
                    }

                    // Generate the network
                    p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);

                    // Set inlet and outlet nodes (if node (start/end) is at the edge of the domain space, make it input or output depending on which side), and assign each non-inlet/outlet vessel a radius from the list
                    auto p_segment = p_network->GetVesselSegments()[0];
                    p_segment->SetRadius(inlet_vessel_radius);
                    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
                    vessels = p_network->GetVessels();
                    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                    {
                        if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the start node equals one (no other segments at nodes)
                        {
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the start node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);  // set the node as the input node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the start node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the ouput node pressure                                
                            }
                        }
                        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the end node equals one (no other segments at nodes)	
                        {
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the end node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);;  // set the node as the input node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the end node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the output node pressure
                            }
                        }
                    }
                
                    // Assign vessel radii based on imported arrays
                    unsigned int vessel_index = 0;
                    while(vessel_index < n_vessels)
                    {
                        p_current_vessel = p_network->GetVessel(id_array[vessel_index][0]);

                        // Get segment radius from array
                        QLength new_vessel_radius(double(radii_array[vessel_index][0])*unit::microns);

                        p_current_vessel->SetRadius(new_vessel_radius);
                        vessel_index++;
                    }
                    vessels = p_network->GetVessels();

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    std::string solver_name;
                    if (h_solver==1)
                    {
                        solver_name = "ConstantHaematocrit"; 
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==2)
                    {
                        solver_name = "PriesHaematocrit"; 
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                
                    }
                    else if (h_solver==3)
                    {
                        solver_name = "MemoryHaematocrit"; 
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==4)
                    {
                        solver_name = "FungHaematocrit"; 
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==5)
                    {
                        solver_name = "YangHaematocrit"; 
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

                    // Prune all vessels up to specified dose 
                    for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
                    { 
                        // Display status message
                        std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                        // Set filename
                        std::stringstream selection_stream;
                        std::stringstream heterogeneity_stream;
                        std::stringstream kill_stream;
                        selection_stream << std::fixed << std::setprecision(0) << to_string(list_number);                    
                        heterogeneity_stream << std::fixed << std::setprecision(0) << alpha;                              
                        kill_stream << std::fixed << std::setprecision(0) << KilledVessels;           
                        std::string selection_string = selection_stream.str();
                        std::string heterogeneity_string = heterogeneity_stream.str();                    
                        std::string kill_string = kill_stream.str();                    
                        std::string file_name = network_name + solver_name + "/Mu" + heterogeneity_string + "/Selection" + selection_string + "/Kills" + kill_string;
                        std::string str_directory_name = file_name;
                        std::cout << str_directory_name << std::endl;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                        // Set up flag for broken solver
                        unsigned broken_solver = 0;

                        // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                        // unsigned max_iter = 5;  // 1000 
                        // double tolerance2 = 1.e-10;
                        std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                        // std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                        // for(unsigned idx=0; idx<max_iter; idx++)
                        // {
                            // Run the solvers
                            p_impedance_calculator->Calculate();
                            p_abstract_haematocrit_solver->Calculate();
                            p_viscosity_calculator->Calculate();
                            flow_solver.SetUp();
                            flow_solver.Solve();

                            // Get the residual
                            for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                            {
                                // Set segments with no flow to have no haematocrit
                                if (fabs(segments[jdx]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                {
                                    segments[jdx]->GetFlowProperties()->SetHaematocrit(0.0);
                                }   
                            }

                        // If solver doesn't converge, move on to next one
                        if (broken_solver == 1)
                        {
                            continue;
                        }

                        // Run the simulation 
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                        p_network->Write(output_file);

                        // Write the PQs
                        // QFlowRate threshold = 5.e-14*unit::metre_cubed_per_second;
                        // QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;
                        QFlowRate threshold = 3.e-13*unit::metre_cubed_per_second;
                        PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
                        outfile.open("/tmp/narain/testoutput/TestHexagonalNetwork/hex_lognormal_individual_perfusion_quotients.txt", std::ios_base::app);
                        outfile << network_name << " " << solver_name << " " << heterogeneity_string << " " << selection_string << " " << kill_string << " " << PerfusionQuotient << " \n"; 
                        outfile.close();

                        // Remove the smallest vessel
                        vessels = p_network->GetVessels();
                        QLength minimum_radius = inlet_vessel_radius;
                        unsigned int minimum_index = 0;
                        for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                        {                            
                            // // Exclude inlets and outlets
                            // if (vessels[vessel_index]->GetStartNode()->GetFlowProperties()->IsInputNode()==0
                            // && vessels[vessel_index]->GetStartNode()->GetFlowProperties()->IsOutputNode()==0
                            // && vessels[vessel_index]->GetEndNode()->GetFlowProperties()->IsInputNode()==0
                            // && vessels[vessel_index]->GetEndNode()->GetFlowProperties()->IsOutputNode()==0)
                            // {
                                // Get the current segment's radius
                                QLength current_radius = vessels[vessel_index]->GetRadius();

                                // If the current radius is less than the minimum radius, record the new minimum
                                if (current_radius < minimum_radius)
                                {
                                    minimum_radius = current_radius;
                                    minimum_index = vessel_index;
                                }     
                            // }
                        }
                        p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel

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

};

#endif /*TESTENHANCEDPERFUSIONFOLLOWINGEXPOSURETORADIOTHERAPY_HPP_*/
