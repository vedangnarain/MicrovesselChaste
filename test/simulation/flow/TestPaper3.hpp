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
#include "CellsGenerator.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "CellLabel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellPopulationAreaWriter.hpp"
#include "CellPopulationElementWriter.hpp"
#include "RandomCellKiller.hpp"
#include "LQRadiotherapyCellKiller.hpp"

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
// #include "SimulationTime.hpp"
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
#include "VesselNetworkReader.hpp"
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

// Cells
#include "CancerCellMutationState.hpp"
#include "StalkCellMutationState.hpp"
#include "QuiescentCancerCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "Owen11CellPopulationGenerator.hpp"
#include "Owen2011TrackingModifier.hpp"
#include "CaBasedCellPopulation.hpp"
#include "ApoptoticCellKiller.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"

// Forces
#include "GeneralisedLinearSpringForce.hpp"

// Angiogenesis
#include "Owen2011SproutingRule.hpp"
#include "Owen2011MigrationRule.hpp"
#include "AngiogenesisSolver.hpp"

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

// Includes standard library
using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tests with cells
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Tests with cells
class TestPaper3 : public AbstractCellBasedWithTimingsTestSuite
{

public:

    // Simulate an experimentally-acquired network on a PDE grid with flow and O2 and experimentally-acquired cells
    void TestExperimentalImage2DWithFlowAndO2AndCells() 
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Set network in filename
        std::string network_name = "TestBiologicalNetworkTissue/";

        // Choose the vascular parameters
        // double dimless_domain_size_x = 85.22; 
        // double dimless_domain_size_y = 85.22;  // based on image dimensions

        double dimless_domain_size_x = 1214.56; 
        double dimless_domain_size_y = 1214.56;  // based on image dimensions
        // double dimless_domain_size_x = 850.22; 
        // double dimless_domain_size_y = 850.22;  // based on image dimensions
        QDynamicViscosity viscosity = 1.96*1.e-3*unit::poiseuille;
        double initial_haematocrit = 0.01;  // conks out somewhere between 0.35 and 0.4
        
        // Match this based on cell size
        QLength cell_grid_spacing = 1.0_um;

        // Seed the random number generator
        RandomNumberGenerator::Instance()->Reseed(12345);

        // Create the log file for broken simulations
        std::ofstream broken_layouts_file;
        broken_layouts_file.open("/scratch/narain/testoutput/TestBiologicalNetwork/broken_layouts.txt");
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
            std::vector<std::vector<unsigned> > Order;

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

            // Manually set the inlets and outlets
            vessels = p_network->GetVessels();
            p_network->GetNode(1711)->GetFlowProperties()->SetIsInputNode(true);
            p_network->GetNode(1711)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));
            p_network->GetNode(1405)->GetFlowProperties()->SetIsOutputNode(true);
            p_network->GetNode(1405)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));

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
                broken_layouts_file.open("/scratch/narain/testoutput/TestBiologicalNetwork/broken_layouts.txt", std::ios_base::app);
                broken_layouts_file << str_directory_name << " \n"; 
                broken_layouts_file.close();

                // Move onto the next selection
                break;
            }

            // Set up the grid calculator
            std::shared_ptr<GridCalculator<2> > p_grid_calc = GridCalculator<2>::Create();
            p_grid_calc->SetGrid(p_grid);
            
            ////////////////////////////////////////////////////////////////

            // josh
            		const int DIM = 2;

            std::vector<Node<2>*> nodes;
            // Stroma Nodes
            unsigned nodeNum=0;
            unsigned numStromaNodes = 0;
            // For staggering rows
            unsigned rownum = 1;
            double offset;
            double node_spacing = 1.0;
            for (double x=0; x<=dimless_domain_size_x; x=x+node_spacing)//x++)
            {
                offset = 0;
                // if (rownum % 2)
                // {
                //     offset = 0;
                // }
                for (double y=0+offset; y<=dimless_domain_size_y; y=y+node_spacing)//y++)
                {
                    if (DIM == 2) 
                    {
                        nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
                        nodeNum++;
                        numStromaNodes++;
                    }
                }
                rownum++;
            }

            // Choose a radii list file and read the file to an array
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


            // Tumour Nodes
            // unsigned numTumourNodes = 0;
            // Seed a cluster of tumour cells in the centre
            // double x = dimless_domain_size_x*0.5;
            // double y = dimless_domain_size_y*0.5;
            // nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
            // nodeNum++;
            // nodes.push_back(new Node<DIM>(nodeNum, false, x-0.5, y));
            // nodeNum++;
            // nodes.push_back(new Node<DIM>(nodeNum, false, x, y-0.5));
            // nodeNum++;
            // nodes.push_back(new Node<DIM>(nodeNum, false, x-0.5, y-0.5));
            // nodeNum++;
            // numTumourNodes = numTumourNodes+4;


				std::vector<CellPtr> cells;
MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		// MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

		unsigned cancerLabelColour = 9;
		// unsigned macrophageLabelColour = 7;
		// MAKE_PTR_ARGS(CellLabel, p_stroma_label, (1));
		MAKE_PTR_ARGS(CellLabel, p_cancer_label, (cancerLabelColour));

for (unsigned i=0; i<numStromaNodes; i++)
		{
                    
                        SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
            p_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            // p_model->SetStemCellG1Duration(8.0);
            // p_model->SetTransitCellG1Duration(8.0);
            double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                 (  p_model->GetStemCellG1Duration()
                                  + p_model->GetSG2MDuration() );
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
            }


		for (unsigned i=numStromaNodes; i<numStromaNodes+numTumourNodes; i++)
		{
                        SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
            p_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            // p_model->SetStemCellG1Duration(8.0);
            // p_model->SetTransitCellG1Duration(8.0);
            p_cell->AddCellProperty(p_cancer_label);
            double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                 (  p_model->GetStemCellG1Duration()
                                  + p_model->GetSG2MDuration() );
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
            }














        //     for (unsigned i=0; i<p_generating_mesh->GetNumNodes(); i++)
        //             {
        //     if (i>=10)
        //             {
        //                 SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
        //     p_model->SetDimension(2);
        //     CellPtr p_cell(new Cell(p_state, p_model));
        //     p_cell->SetCellProliferativeType(p_stem_type);
        //     // p_model->SetStemCellG1Duration(8.0);
        //     // p_model->SetTransitCellG1Duration(8.0);
        //     double birth_time = - RandomNumberGenerator::Instance()->ranf() *
        //                          (  p_model->GetStemCellG1Duration()
        //                           + p_model->GetSG2MDuration() );
        //     p_cell->SetBirthTime(birth_time);
        //     cells.push_back(p_cell);
        //     }
        // }




            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 1.5);
            NodeBasedCellPopulation<2> cell_population(mesh, cells);
                                                    std::cout << "This is line number: " << __LINE__ << std::endl;
		cell_population.SetAbsoluteMovementThreshold(DBL_MAX);//Set big movement threshold

		// Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
		// ChastePoint<DIM> lower(0, 0, 0);
		// ChastePoint<DIM> upper(dimless_domain_size_x,dimless_domain_size_y, 0);
		// MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));








            // int dimless_domain_size_x_int = static_cast<int>(20);
            //                                         std::cout << "This is line number: " << __LINE__ << std::endl;
            // double scale_f = 1.0;
            // HoneycombMeshGenerator generator(dimless_domain_size_x_int, dimless_domain_size_x_int, 0, scale_f);
            //                                         std::cout << "This is line number: " << __LINE__ << std::endl;

            // boost::shared_ptr<MutableMesh<2,2>> p_generating_mesh = generator.GetMesh();
            // boost::shared_ptr<MutableMesh<2, 2>> p_generating_mesh(generator.GetMesh());




      
            // std::vector<CellPtr> cells;
        // MAKE_PTR(WildTypeCellMutationState, p_state);
        // // MAKE_PTR(StemCellProliferativeType, p_stem_type);
        //     for (unsigned i=0; i<p_generating_mesh->GetNumNodes(); i++)
        //             {
        //     if (i>=10)
        //             {
        //                 SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
        //     p_model->SetDimension(2);
        //     CellPtr p_cell(new Cell(p_state, p_model));
        //     p_cell->SetCellProliferativeType(p_stem_type);
        //     // p_model->SetStemCellG1Duration(8.0);
        //     // p_model->SetTransitCellG1Duration(8.0);
        //     double birth_time = - RandomNumberGenerator::Instance()->ranf() *
        //                          (  p_model->GetStemCellG1Duration()
        //                           + p_model->GetSG2MDuration() );
        //     p_cell->SetBirthTime(birth_time);
        //     cells.push_back(p_cell);
        //     }
        // }

                                                    std::cout << "This is line number: " << __LINE__ << std::endl;

            // CellsGenerator<UniformCellCycleModel, 2> cell_population_genenerator;
            // cells_generator.SetGridCalculator(p_grid_calc);            
            // cell_population_genenerator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_state);


            // OffLatticeSimulation<2> off_simulator(cell_population);
            // // simulator.SetOutputDirectory("NodeBasedMonolayer");
            // // simulator.SetSamplingTimestepMultiple(12);
            // // simulator.SetEndTime(10.0);

            // MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
            // off_simulator.AddForce(p_force);

            ////////////////////////////////////////////////////////////////


            // Set up the cell populations
            // std::shared_ptr<Owen11CellPopulationGenerator<2> > p_cell_population_genenerator = Owen11CellPopulationGenerator<2>::Create();
            // p_cell_population_genenerator->SetGridCalculator(p_grid_calc);
            // p_cell_population_genenerator->SetVesselNetwork(p_network);
            // QLength tumour_radius(10000.0 * unit::microns);  // just make it bigger than the domain to ensure all cells are tumour cells 
            // p_cell_population_genenerator->SetTumourRadius(tumour_radius);
            // std::shared_ptr<CaBasedCellPopulation<2> > p_cell_population = p_cell_population_genenerator->Update();
            
            // Set up the VEGF PDE as we set the O2 one (currently turned off)
            auto p_vegf_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
            p_vegf_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpVegfDiffusivity->GetValue("User"));
            p_vegf_pde->SetContinuumLinearInUTerm(-Owen11Parameters::mpVegfDecayRate->GetValue("User"));
        
            // Set up a map for different release rates depending on cell type. Also include a threshold intracellular VEGF below which there is no release (currently turned off)
            auto p_normal_and_quiescent_cell_source = CellStateDependentDiscreteSource<2>::Create();  // initialise the source (which combines normal and quiescent cells)
            std::map<unsigned, QConcentrationFlowRate > normal_and_quiescent_cell_rates;  // map to store secretion rates
            std::map<unsigned, QConcentration > normal_and_quiescent_cell_rate_thresholds;  // map to store secretion thresholds
            // MAKE_PTR(QuiescentCancerCellMutationState, p_quiescent_cancer_state);
            // MAKE_PTR(WildTypeCellMutationState, p_normal_cell_state);
            normal_and_quiescent_cell_rates[p_state->GetColour()] = Owen11Parameters::mpCellVegfSecretionRate->GetValue("User");  // set the rate for normal cells
            normal_and_quiescent_cell_rate_thresholds[p_state->GetColour()] = 0.27*unit::mole_per_metre_cubed;  // set the threshold for normal cells
            normal_and_quiescent_cell_rates[p_state->GetColour()] = Owen11Parameters::mpCellVegfSecretionRate->GetValue("User");  // set the rate for quiescent cells
            normal_and_quiescent_cell_rate_thresholds[p_state->GetColour()] = 0_M;  // set the threshold for quiescent cells
            p_normal_and_quiescent_cell_source->SetStateRateMap(normal_and_quiescent_cell_rates);  // configure the map
            p_normal_and_quiescent_cell_source->SetLabelName("VEGF");  // set the label for the source
            p_normal_and_quiescent_cell_source->SetStateRateThresholdMap(normal_and_quiescent_cell_rate_thresholds);

            // Print the cell types and their colours
            MAKE_PTR(CancerCellMutationState, p_cancer_state);
            // MAKE_PTR(StemCellProliferativeType, p_stem_state);
            // std::cout << "Normal state has Legacy Cell type " <<  p_normal_cell_state->GetColour() << std::endl;
            std::cout << "Cancer state has Legacy Cell type " <<  p_cancer_state->GetColour() << std::endl;
            // std::cout << "Quiescent cancer state has Legacy Cell type " <<  p_quiescent_cancer_state->GetColour() << std::endl;
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

            // Specifies which extracellular fields to upadate based on PDE
            boost::shared_ptr<MicrovesselSimulationModifier<2> > p_microvessel_modifier =
            boost::shared_ptr<MicrovesselSimulationModifier<2> >(new MicrovesselSimulationModifier<2> ());
            p_microvessel_modifier->SetMicrovesselSolver(p_microvessel_solver);
            std::vector<std::string> update_labels;
            update_labels.push_back("oxygen");
            p_microvessel_modifier->SetCellDataUpdateLabels(update_labels);

            ////////////////////////////////////////////////////////////////
            
            
            OffLatticeSimulation<2> simulator(cell_population);
            simulator.AddSimulationModifier(p_microvessel_modifier);

            // simulator.SetOutputDirectory("NodeBasedMonolayer");
            // simulator.SetSamplingTimestepMultiple(12);
            // simulator.SetEndTime(10.0);
boost::shared_ptr<ApoptoticCellKiller<2>> p_apoptotic_cell_killer(new ApoptoticCellKiller<2>(&cell_population));
            simulator.AddCellKiller(p_apoptotic_cell_killer);
boost::shared_ptr<LQRadiotherapyCellKiller<2>> p_rt_killer =
    boost::shared_ptr<LQRadiotherapyCellKiller<2>>(new LQRadiotherapyCellKiller<2>(&cell_population));





            // MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
            // simulator.AddForce(p_force);


            ////////////////////////////////////////////////////////////////




            // The full simulation is run as a typical Cell Based Chaste simulation
            // OnLatticeSimulation<2> simulator(*p_cell_population);
            // simulator.AddSimulationModifier(p_microvessel_modifier);

            // Add a LQ RT killer
            p_rt_killer->SetDoseInjected(2.0*unit::gray);  // modify dose here
            p_rt_killer->SetOerAlphaMax(1.75);
            p_rt_killer->SetOerAlphaMin(1.0);
            p_rt_killer->SetOerBetaMax(3.25);
            p_rt_killer->SetOerBetaMin(1.0);
            p_rt_killer->SetOerConstant(oxygen_solubility_at_stp * 3.28_Pa);  // adjust the radiobiological threshold here
            // p_rt_killer->SetOerConstant(oxygen_solubility_at_stp * 3.00_Pa);  // adjust the radiobiological threshold here
            p_rt_killer->SetAlphaMax(0.3 * unit::per_gray);
            p_rt_killer->SetBetaMax(0.03 * unit::per_gray_squared);
            p_rt_killer->UseOer(true);  // turn OER on and off here

            // Set up the dosage times
            std::vector<QTime > rt_times;
            rt_times.push_back(3600.0*3.0*unit::seconds);
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
            // boost::shared_ptr<Owen2011TrackingModifier<2> > p_owen11_tracking_modifier(new Owen2011TrackingModifier<2>);
            // simulator.AddSimulationModifier(p_owen11_tracking_modifier);

            // Set up the simulation time and run it
            simulator.SetOutputDirectory(str_directory_name);
            simulator.SetSamplingTimestepMultiple(1);
            simulator.SetDt(1.0);
            simulator.SetEndTime(10.0);  // ten hours
                                       std::cout << "This is line number: " << __LINE__ << std::endl;
 simulator.Solve();
                            std::cout << "This is line number: " << __LINE__ << std::endl;

            // Write the output network file (visualise with Paraview: set Filters->Alphabetical->Tube)
            std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
            p_network->Write(output_file);
                                        std::cout << "This is line number: " << __LINE__ << std::endl;

            // Dump our parameter collection to an xml file and, importantly, clear it for the next test
            ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
            ParameterCollection::Instance()->Destroy();
            BaseUnits::Instance()->Destroy();
            SimulationTime::Instance()->Destroy();
                            std::cout << "This is line number: " << __LINE__ << std::endl;

            // If simulation doesn't converge, move on to next layout
            if (broken_solver == 1)
            {
                break;
            }
                                        std::cout << "This is line number: " << __LINE__ << std::endl;

        }
                            std::cout << "This is line number: " << __LINE__ << std::endl;

        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

};

#endif /*TESTPAPER3_HPP_*/