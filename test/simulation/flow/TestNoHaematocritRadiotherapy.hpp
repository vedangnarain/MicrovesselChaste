/*

Copyright (c) 2005-2017, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,v
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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*

I've tried to add more detailed comments. At the moment, I'm not sure what several bits of the code do and I've marked them with a '???'

I addded line-by-line comments when I was trying to understand Jakub's code. If you don't know what something is doing, try going back to earlier instances of the code block, usually in preceding tests with the same architecture.

The tests work by pruning vessels over a range of radii and computing the optimum radius threshold for max. perfusion

There are only six tests of importance, used to generate the figures in the draft paper. They've been marked in the comments.

Simulation outputs can usually found in /tmp/<home_directory>/testoutput.

RT = Radiotherapy
PQ = Perfusion Quotient
Fig. = Figure
# = number

TO DO LIST

Check tutorials and check whether I need to indent everything after the first {.

Move hexagonal tests from Biological NWs to appropriate section after testing.

Add comments to second Voronoi and second Hex simulation.

12/3/21

*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialisation
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef TESTNOHAEMATOCRITRADIOTHERAPY_HPP
#define TESTNOHAEMATOCRITRADIOTHERAPY_HPP

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
#include "ConstantHaematocritSolver.hpp"
#include "UnitCollection.hpp"
#include "RegularGrid.hpp"
#include "SimulationTime.hpp"
#include "MicrovesselSolver.hpp"
#include "SimpleLinearEllipticFiniteDifferenceSolver.hpp"
#include "SimpleLinearEllipticFiniteElementSolver.hpp"
#include "VesselBasedDiscreteSource.hpp"
#include "DiscreteContinuumLinearEllipticPde.hpp"
#include "Owen11Parameters.hpp"
#include "ViscosityCalculator.hpp"
#include "Secomb04Parameters.hpp"
#include "GenericParameters.hpp"
#include "StructuralAdaptationSolver.hpp"
#include "VesselNetworkPropertyManager.hpp"
#include "VesselNetworkGeometryCalculator.hpp"
#include <sstream>
#include <fstream>
#include<vector>
#include<iostream>
#include <string>
#include "PetscAndVtkSetupAndFinalize.hpp"

using namespace std;

class TestNoHaematocritRadiotherapy : public CxxTest::TestSuite
{

public:

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Forking Networks
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// The following two tests are used to generate Fig. 4, 5, 6 and the data is used for Fig. 7

// Test to run forking simulations with constant pressure difference (used for figures)
void TestRadiotherapyDichotomousRemoveMethodsUnifiedKilling()
{	
	// Define the key parameters
    unsigned order=6; 										// number of vessel generations (the inlet/outlet vessels are Gen 0)
    QLength max_radius = 50_um;  							// radius of inlet/outlet vessels
    double twicelambda = 8.0;  								// lambda = length/diameter = 4 (multiplied by 2)
    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;  	// viscosity
    double inlet_haematocrit = 0.45;  						// inlet haematocrit

    // Generate the network
    VesselNetworkGenerator<2> network_generator;
    QLength quarter_vertical_length = 0.9*twicelambda*max_radius*pow(2.0,-1.0/3.0);  // Gen 1 height  
    double dimless_length = 1.0;
	
	// Initialise the outputs
	std::ofstream outfile;
    outfile.open("testRemoveMethodsForkingUnifiedKilling.txt");//, std::ios_base::app);  // for Fig. 7a
    std::ofstream outfile2;
    outfile2.open("optimumRemoveMethodsForkingUnifiedKilling.txt");//, std::ios_base::app);  // for Fig 7b

	// ??? Non-dimensionalising the length
    for(unsigned i_aux=1; i_aux<order+1; i_aux++)  // for all generations
    {
		dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
    }

	// Create a folder in which to store .vtp outputs
    auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Dichotomous_Heterogeneous_RemoveMethods_Order6_UnifiedKilling", true);

	// Loop over various heterogeneities to generate networks and prune them
	//for(unsigned j_aux=0; j_aux<11; j_aux++)
    for(unsigned j_aux=0; j_aux<6; j_aux++)  // for all generations
    {
		// Generate and save the networks
    	//double alpha = 1.0+(double)j_aux*0.05;
    	double alpha = 1.0+(double)j_aux*0.1;  // alpha determines the heterogeneity
    	std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, quarter_vertical_length, max_radius, alpha, twicelambda);  // generate the network
		std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");  // saves the network without flow 
		p_network->Write(output_file_initial);
    	std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetwork"+to_string(alpha)+".vtp");  // saves the network after adding flow
    	std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkRT"+to_string(alpha)+".vtp");  // saves the post-RT network with flow

		// Set the haematocrit for all vessels
		std::vector<std::shared_ptr<Vessel<2> > > vessels = p_network->GetVessels();
    	typename std::vector<std::shared_ptr<Vessel<2> > >::iterator it;
		for(it = vessels.begin(); it != vessels.end(); it++)  
    	{
        	(*it)->GetFlowProperties()->SetHaematocrit(inlet_haematocrit); 
    	}
		
		// Set input and output nodes and pressures
		QLength domain_side_length_x = dimless_length*2.0*twicelambda*max_radius;
    	VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(0.0_um,2.0*quarter_vertical_length));  // pick inlet node
    	VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(domain_side_length_x, 2.0*quarter_vertical_length));  // pick outlet node
		p_inlet_node->GetFlowProperties()->SetIsInputNode(true);  	 // set input node
    	p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa); 	 // set inlet pressure
    	p_outlet_node->GetFlowProperties()->SetIsOutputNode(true); 	 // set output node
    	p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);  // set outlet pressure

		// Set up a viscosity solver
        auto p_viscosity_calculator = ViscosityCalculator<2>::Create();	 // create a viscosity calculator
        p_viscosity_calculator->SetPlasmaViscosity(viscosity);  		 // set plasma viscosity
        p_viscosity_calculator->SetVesselNetwork(p_network);			 // assign calculator to network
        p_viscosity_calculator->Calculate(); 							 // calculate the viscosity
        
		// Set up the impedance solver
		auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();  // create an impedance calculator
    	p_impedance_calculator->SetVesselNetwork(p_network); 				   // assign calculator to network
     	p_impedance_calculator->Calculate();								   // calculate the impedance

		// Set up the flow solver 
		FlowSolver<2> flow_solver;
    	flow_solver.SetVesselNetwork(p_network);
    	flow_solver.SetUp();
        flow_solver.SetUseDirectSolver(false);
    	flow_solver.Solve();

		// Get pre-RT Perfusion Quotient
		QFlowRate threshold = 1.e-12*unit::metre_cubed_per_second;  // set flow threshold for what counts as a perfused vessel
		double PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);  		// get perfusion quotient
    	//std::cout << "Perfusion quotient before RT is : " << PerfusionQuotient << " \n"; 
		p_network->Write(output_file_final);

		// ???
		/*
		//unsigned NV_Total = pow(2,order+2)-2;
    	//unsigned NV_ToKillFrom = NV_Total;
    	//unsigned KilledVessels = 0;
    	//double dose = 0.5;
    	//unsigned NV_ToKill = (unsigned)(dose*NV_ToKillFrom);
		//  quarter_vertical_length
		*/

		// Set up pruning
        QLength max_radius_to_kill = 50_um;  // set threshold up to which vessels should be pruned
		QLength radius_step = 1_um;
		QLength current_radius_threshold = 0_um;
		//outfile << alpha << " " << current_radius_threshold << " " << PerfusionQuotient << " \n"; 	
    	double optimumBeta =PerfusionQuotient;  // ???
    	QLength optimumThreshold = 0_um;  		// ???

		// Prunes vessels from thinnest up to simulate RT
		while( current_radius_threshold <= max_radius_to_kill )
		{
    		//auto minRadius = (*vessels.begin())->GetRadius();
    		//unsigned IdToKill = 0;
			/*
    		for(it = vessels.begin(); it != vessels.end(); it++)
    		{	
			if  ((minRadius > (*it)->GetRadius())&&((*it)->GetRadius()>0.00000001_um))
			{
				minRadius = (*it)->GetRadius();
				IdToKill = (*it)->GetId();
			}
    		}
 			vessels[IdToKill]->SetRadius(0.000000000001_um);
    		vessels[IdToKill]->SetToDie();
			*/ /*
    		for(it = vessels.begin(); it != vessels.end(); it++)
    		{	
			if  (minRadius > (*it)->GetRadius())
			{
				minRadius = (*it)->GetRadius();
				IdToKill = (*it)->GetId();
			}
    		}
			*/
 			
			// Remove vessels under radius threshold
			p_network->RemoveThinVessels(current_radius_threshold, false);  
			
			// Calculate new flow
			vessels = p_network->GetVessels();								// get leftover vessels
    		p_viscosity_calculator->Calculate(); 							// calculate the viscosity
    		p_impedance_calculator->Calculate();							// calculate the impedance
    		flow_solver.SetUp();
    		flow_solver.Solve();
    		//std::cout << "The perfusion coefficient for the radius threshold of " << current_radius_threshold << ", is  " << p_network->GetPerfusionQuotientBeta(threshold) << "\n";
			outfile << alpha << " " << current_radius_threshold << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
    
			// ??? Figure out the optimum radius threshold for maximum perfusion? 
    		if(p_network->GetPerfusionQuotientBeta(threshold)>optimumBeta)
    		{
    			optimumBeta=p_network->GetPerfusionQuotientBeta(threshold);
    			optimumThreshold = current_radius_threshold;
				//std::cout << "The perfusion coefficient for the radius threshold of blablabla" << current_radius_threshold << ", is blablabla " << p_network->GetPerfusionQuotientBeta(threshold) << "  and current optimum threshold is: " << optimumThreshold << "\n";
    		}

			// ??? Assuming radius doesn't exceed threshold, write the network to the output file
			if(current_radius_threshold>12.9_um && current_radius_threshold<13.1_um)
			{
				p_network->Write(output_file_final_RT);
			}
			current_radius_threshold = current_radius_threshold + radius_step;  // move onto next radius threshold for next step
		}
    	outfile2 << alpha << " " << optimumThreshold << " \n";  // print the optimum threshold
    }
}

// Test to run forking simulations with constant flow rate (used for figures)
void xTestRadiotherapyDichotomousRemoveMethodsUnifiedKillingConstantFlux()
{
    // Define the key parameters
	//unsigned order=7;
    unsigned order=6;
    QLength max_radius = 50_um;
    double twicelambda = 8.0;
    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
    double inlet_haematocrit = 0.45;

    // Generate the network
    VesselNetworkGenerator<2> network_generator;
    QLength quarter_vertical_length = 0.9*twicelambda*max_radius*pow(2.0,-1.0/3.0);
    double dimless_length = 1.0;
    std::ofstream outfile;

	// Initialise the outputs
    outfile.open("testRemoveMethodsForkingUnifiedKillingConstantFlux.txt");//, std::ios_base::app);  // for Fig. 7c
    std::ofstream outfile2;
    outfile2.open("optimumRemoveMethodsForkingUnifiedKillingConstantFlux.txt");//, std::ios_base::app);  // for Fig. 7d
	
	// ??? Non-dimensionalising the length
    for(unsigned i_aux=1; i_aux<order+1; i_aux++)
    {
		dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
    }

	// Create a folder in which to store .vtp outputs  
    auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Dichotomous_Heterogeneous_RemoveMethods_Order6_UnifiedKilling_ConstantFlux", true);
    
	// Loop over various heterogeneities
	//for(unsigned j_aux=0; j_aux<11; j_aux++)
    for(unsigned j_aux=0; j_aux<6; j_aux++)  // for all generations
    {
    	// Generate and save the networks
		//double alpha = 1.0+(double)j_aux*0.05;
    	double alpha = 1.0+(double)j_aux*0.1;
    	std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, quarter_vertical_length, max_radius, alpha, twicelambda);
		std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");  // saves the network without flow 
		p_network->Write(output_file_initial);
    	std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetwork"+to_string(alpha)+".vtp");  // saves the network after adding flow
    	std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkRT"+to_string(alpha)+".vtp");  // saves the post-RT network with flow
    	
		// Set the haematocrit for all vessels
		std::vector<std::shared_ptr<Vessel<2> > > vessels = p_network->GetVessels();
    	typename std::vector<std::shared_ptr<Vessel<2> > >::iterator it;
		for(it = vessels.begin(); it != vessels.end(); it++)
    	{
        	(*it)->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    	}
		
		// Set input and output nodes and pressures
		QLength domain_side_length_x = dimless_length*2.0*twicelambda*max_radius;
    	VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(0.0_um,2.0*quarter_vertical_length));
    	VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(domain_side_length_x, 2.0*quarter_vertical_length));
    	p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
        p_inlet_node->GetFlowProperties()->SetUseVelocityBoundaryCondition(true);  // ??? not sure what this does
        p_inlet_node->GetSegment(0)->GetFlowProperties()->SetFlowRate(2.e-10*unit::metre_cubed_per_second);  // set the flow rate
    	p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
    	p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);
		
		// Set up a viscosity solver
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
    	flow_solver.SetUp();
        flow_solver.SetUseDirectSolver(false);
    	flow_solver.Solve();
		
		// Get pre-RT Perfusion Quotient
		QFlowRate threshold = 1.e-12*unit::metre_cubed_per_second;
		double PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
    	//std::cout << "Perfusion quotient before RT is : " << PerfusionQuotient << " \n"; 
		p_network->Write(output_file_final);
		
		// ??? seems to be used in non-UnifiedKilling tests
		/*
		//unsigned NV_Total = pow(2,order+2)-2;
    	//unsigned NV_ToKillFrom = NV_Total;
    	//unsigned KilledVessels = 0;
    	//double dose = 0.5;
    	//unsigned NV_ToKill = (unsigned)(dose*NV_ToKillFrom);
  		//  quarter_vertical_length
		*/

		// Set up pruning
        QLength max_radius_to_kill = 12.5_um;
		QLength radius_step = 1_um;
		QLength current_radius_threshold = 0_um;
		//outfile << alpha << " " << current_radius_threshold << " " << PerfusionQuotient << " \n"; 
    	double optimumBeta =PerfusionQuotient;
    	QLength optimumThreshold = 0_um;

		// Prunes vessels from thinnest up to simulate RT
		while( current_radius_threshold <= max_radius_to_kill )
		{
    		//auto minRadius = (*vessels.begin())->GetRadius();
    		//unsigned IdToKill = 0;
			/*
    		for(it = vessels.begin(); it != vessels.end(); it++)
    		{	
			if  ((minRadius > (*it)->GetRadius())&&((*it)->GetRadius()>0.00000001_um))
			{
				minRadius = (*it)->GetRadius();
				IdToKill = (*it)->GetId();
			}
    		}
 			vessels[IdToKill]->SetRadius(0.000000000001_um);
    		vessels[IdToKill]->SetToDie();
			*/ /*
    		for(it = vessels.begin(); it != vessels.end(); it++)
    		{	
			if  (minRadius > (*it)->GetRadius())
			{
				minRadius = (*it)->GetRadius();
				IdToKill = (*it)->GetId();
			}
    		}
			*/

			// Remove vessels under radius threshold
		 	p_network->RemoveThinVessels(current_radius_threshold, false);
			
			// Calculate new flow
			vessels = p_network->GetVessels();
    		p_viscosity_calculator->Calculate();
    		p_impedance_calculator->Calculate();
    		flow_solver.SetUp();
    		flow_solver.Solve();
    		//std::cout << "The perfusion coefficient for the radius threshold of " << current_radius_threshold << ", is  " << 					p_network->GetPerfusionQuotientBeta(threshold) << "\n";
			outfile << alpha << " " << current_radius_threshold << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
			
			// ??? Figure out the optimum radius threshold for maximum perfusion? 
    		if(p_network->GetPerfusionQuotientBeta(threshold)>optimumBeta)
    		{
    			optimumBeta=p_network->GetPerfusionQuotientBeta(threshold);
    			optimumThreshold = current_radius_threshold;
				std::cout << "The perfusion coefficient for the radius threshold of blablabla" << current_radius_threshold << ", is blablabla " << p_network->GetPerfusionQuotientBeta(threshold) << "  and current optimum threshold is: " << optimumThreshold << "\n";
    		}
			current_radius_threshold = current_radius_threshold + radius_step;  // move onto next radius threshold for next step
		}
		outfile2 << alpha << " " << optimumThreshold << " \n"; 
		p_network->Write(output_file_final_RT);
    }
}

// Unused Tests

void xTestRadiotherapyDichotomous()
{

	// Define the key parameters
    unsigned order=7;  // number of vessel generations (the inlet/outlet vessels are generation 0)
    double theta = 1.0;
    QLength max_radius = 50_um;  // radius of inlet/outlet vessels
    double lambda = 8.0;  // lambda = length/diameter
    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;  // viscosity
    double inlet_haematocrit = 0.45;  // inlet haematocrit

    // Generate the network
    VesselNetworkGenerator<2> network_generator;

	// Specify the domain
    //QLength domain_side_length_x = 800_um;
    //QLength domain_side_length_y = 800_um;
    
	// Calculate the vessel length using lambda
	QLength vessel_length = 0.9*lambda*max_radius*pow(2.0,-1.0/3.0);
	
	// Calculate the dimensionless length for all the generations
    double dimless_length = 1.0;
    for(unsigned i_aux=1; i_aux<order+1; i_aux++)
	{
	dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));  // ???
    }

	// ???
    auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Dichotomous_Heterogeneous", true);

	// ??? Do something 11 times
    for(unsigned j_aux=0; j_aux<11; j_aux++)
	{
    double alpha = 1.0+(double)j_aux*0.05; 
    std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, vessel_length, max_radius, alpha, theta, lambda);
    std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");
    p_network->Write(output_file_initial);
    std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetwork.vtp");
    std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkRT.vtp");
    std::vector<std::shared_ptr<Vessel<2> > > vessels = p_network->GetVessels();
    typename std::vector<std::shared_ptr<Vessel<2> > >::iterator it;
    
	for(it = vessels.begin(); it != vessels.end(); it++)
    {
        (*it)->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    }

    //auto p_segment = p_network->GetVessels()[0];

    //p_segment->GetFlowProperties()->SetViscosity(viscosity);
    //VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);

    QLength domain_side_length_x = dimless_length*2.0*lambda*max_radius;
    VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
            Vertex<2>(0.0_um,2.0*vessel_length));
    VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
            Vertex<2>(domain_side_length_x, 2.0*vessel_length));
    p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
    p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
    p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
    p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);



//VesselPtr<2> p_thin_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
 //       Vertex<2>(2.8*real_length_x/4.0,3.0*real_length_y/4.0));
//VesselPtr<2> p_thin_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
//        Vertex<2>(200_um,200_um));
//p_radiated_vessel->SetToDie();
//p_thin_vessel->SetRadius(5_um);
//p_thin_vessel->SetRadius(0.000001_um);
   
    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
    p_viscosity_calculator->SetVesselNetwork(p_network);
    p_viscosity_calculator->Calculate();
    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
    p_impedance_calculator->SetVesselNetwork(p_network);
    p_impedance_calculator->Calculate();

    FlowSolver<2> flow_solver;
    flow_solver.SetVesselNetwork(p_network);
    flow_solver.SetUp();
    flow_solver.SetUseDirectSolver(true);
    flow_solver.Solve();

    //QFlowRate threshold = 3.e-12*unit::metre_cubed_per_second;
    QFlowRate threshold = 1.e-12*unit::metre_cubed_per_second;

    double PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
    std::cout << "Perfusion quotient before RT is : " << PerfusionQuotient << " \n"; 



   p_network->Write(output_file_final);


    unsigned NV_Total = pow(2,9)-2;
    unsigned NV_ToKillFrom = NV_Total;
    unsigned KilledVessels = 0;
    double dose = 0.5;
    unsigned NV_ToKill = (unsigned)(dose*NV_ToKillFrom);
    
    std::ofstream outfile;

    outfile.open("test.txt", std::ios_base::app);
    outfile << alpha << " " << KilledVessels << " " << PerfusionQuotient << " \n"; 

    std::ofstream outfile2;

    outfile2.open("optimum.txt", std::ios_base::app);

    double optimumBeta =PerfusionQuotient;
    unsigned optimumKilled = 0;


    while( KilledVessels < NV_ToKill ) {
    
    auto minRadius = (*vessels.begin())->GetRadius();
    //std::cout << " The first radius is  " << minRadius << " \n";
    unsigned IdToKill = 0;

    for(it = vessels.begin(); it != vessels.end(); it++)
    {	
//std::cout << "the id at this stage is " << (*it)->GetId() << " \n";
	if  ((minRadius > (*it)->GetRadius())&&((*it)->GetRadius()>0.00000001_um)) {
	minRadius = (*it)->GetRadius();
	IdToKill = (*it)->GetId();
        //std::cout << "the new id is " << IdToKill << " \n";
	}
    }

    vessels[IdToKill]->SetRadius(0.000000000001_um);
    vessels[IdToKill]->SetToDie();
    KilledVessels++;

    p_viscosity_calculator->Calculate();
    p_impedance_calculator->Calculate();
    flow_solver.SetUp();
    flow_solver.Solve();
    //std::cout << "The minimum radius is: " << minRadius << " and the vessel id to kill is" << IdToKill << " \n"; 
    //std::cout << "So far killed"  << KilledVessels << "\n";
    std::cout << "The perfusion coefficient for " << KilledVessels << " killed vessels is: " << p_network->GetPerfusionQuotientBeta(threshold) << "\n";

    outfile << alpha << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
    
    if(p_network->GetPerfusionQuotientBeta(threshold)<optimumBeta)
    {
    optimumBeta=p_network->GetPerfusionQuotientBeta(threshold);
    optimumKilled = KilledVessels;
    }

}

    outfile2 << alpha << " " << optimumKilled << " \n"; 

}

//p_network->Write(output_file_final_RT);

//std::cout << "Inlet absolute flow rate is: " <<     p_segment->GetFlowProperties()->GetFlowRate() << "\n";


/* Radiotherapy

    p_thin_vessel->SetRadius(0.000000001_um);
   
    auto p_viscosity_calculatorRT = ViscosityCalculator<2>::Create();
    p_viscosity_calculatorRT->SetPlasmaViscosity(viscosity);
    p_viscosity_calculatorRT->SetVesselNetwork(p_network);
    p_viscosity_calculatorRT->Calculate();
    auto p_impedance_calculatorRT = VesselImpedanceCalculator<2>::Create();
    p_impedance_calculatorRT->SetVesselNetwork(p_network);
    p_impedance_calculatorRT->Calculate();

    FlowSolver<2> flow_solverRT;
    flow_solverRT.SetVesselNetwork(p_network);
    flow_solverRT.SetUp();
    //flow_solver.SetUseDirectSolver(false);
    flow_solverRT.Solve();


    p_network->Write(output_file_final_RT);
std::cout << "Inlet absolute flow rate post radiation is: " <<     p_segment->GetFlowProperties()->GetFlowRate() << "\n";
*/

}

void xTestRadiotherapyDichotomousRandom()
{
    srand(unsigned(time(NULL)));
    unsigned order=7;
    //double theta = 1.0;
    QLength max_radius = 50_um;
    double lambda = 8.0;
    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
    double inlet_haematocrit = 0.45;

    //unsigned NV_Total = pow(2,9)-2;
    unsigned NV_ToKillFrom = pow(2,8)*(pow(2,1)-1);
    double dose = 0.5;
    unsigned NV_ToKill = (unsigned)(dose*NV_ToKillFrom);

    // Generate the network

    VesselNetworkGenerator<2> network_generator;
    QLength vessel_length = 0.9*lambda*max_radius*pow(2.0,-1.0/3.0);
    double dimless_length = 1.0;
    double dimless_offset_x = 1.0;

    for(unsigned i_aux=1; i_aux<order+1; i_aux++)
    {
	dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
    }

    for(unsigned j_aux=1; j_aux<order; j_aux++)
    {
	dimless_offset_x += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(j_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(j_aux-1)));
    }

    QLength domain_side_length_x = dimless_length*2.0*lambda*max_radius;
    QLength offset_x = dimless_offset_x*lambda*max_radius;

    auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Dichotomous_Heterogeneous_Random", true);

    std::ofstream outfile;
    outfile.open("TestRandom.txt");//, std::ios_base::app);
    //std::ofstream outfile2;
    //outfile2.open("OptimumRandom.txt");//, std::ios_base::app)

    unsigned attempts = 100;
    unsigned HetVectorSize = 11;
    //std::vector<std::vector<std::vector<double> > > QuotientMatrixSampling(NHet,std::vector<std::vector<double> >(ToBeKilled+1,std::vector <double>(attempts,0)));

    std::vector<std::vector<double> > QuotientMatrixMeanForking(HetVectorSize,std::vector<double>(NV_ToKill+1));
    std::ofstream outfileMeanForking;
    //std::ofstream outfile2Mean;

    
    std::vector<std::vector<double> > QuotientMatrixAggregateForking(HetVectorSize,std::vector<double>(NV_ToKill+1,0.0));

    std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");
    std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetwork.vtp");
    std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkRT.vtp");

    //QFlowRate threshold = 5.e-13*unit::metre_cubed_per_second;
    QFlowRate threshold = 1.e-12*unit::metre_cubed_per_second;
    std::cout << "Blablabla: " << threshold << " \n";
    for(unsigned att=0; att < attempts;att++)
    {

    	for(unsigned j_aux=0; j_aux<HetVectorSize; j_aux++)
    	//for(unsigned j_aux=0; j_aux<1; j_aux++)
	{
			
    		unsigned KilledVessels = 0;
    		double alpha = 1.0+(double)j_aux*0.05;
    		std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, vessel_length, max_radius, alpha, lambda);
    		//auto p_segment = p_network->GetVesselSegments()[0];
    		//p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    		//VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);
    		p_network->Write(output_file_initial);

    		std::vector<std::shared_ptr<Vessel<2> > > vessels = p_network->GetVessels();
    		typename std::vector<std::shared_ptr<Vessel<2> > >::iterator it;
    		for(it = vessels.begin(); it != vessels.end(); it++)
    		{
        		(*it)->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    		}

    		//auto p_segment = p_network->GetVessels()[0];
    		//p_segment->GetFlowProperties()->SetViscosity(viscosity);
    		//VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);

    		VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
            		Vertex<2>(0.0_um,2.0*vessel_length));
   		VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
           		Vertex<2>(domain_side_length_x, 2.0*vessel_length));
    		p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
    		p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
    		p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
    		p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);
   
    		auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
    		p_viscosity_calculator->SetPlasmaViscosity(viscosity);
    		p_viscosity_calculator->SetVesselNetwork(p_network);
    		p_viscosity_calculator->Calculate();
    		auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
    		p_impedance_calculator->SetVesselNetwork(p_network);
    		p_impedance_calculator->Calculate();

    		FlowSolver<2> flow_solver;
    		flow_solver.SetVesselNetwork(p_network);
    		flow_solver.SetUp();
    		flow_solver.SetUseDirectSolver(true);
    		flow_solver.Solve();
		double PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
    		//std::cout << "Perfusion quotient before RT is : " << PerfusionQuotient << " \n";
    		p_network->Write(output_file_final);
		
    		QuotientMatrixAggregateForking[j_aux][KilledVessels] += p_network->GetPerfusionQuotientBeta(threshold);
    
    		outfile << alpha << " " << KilledVessels << " " << PerfusionQuotient << " \n"; 

    		double optimumBeta =PerfusionQuotient;
    		//unsigned optimumKilled = 0;
		while( KilledVessels < NV_ToKill )
		{
    
			double x_norm = (double)rand()/RAND_MAX;
    			double y_norm = (double)rand()/RAND_MAX;
    			QLength x_Kill = offset_x + x_norm*(domain_side_length_x-2.0*offset_x);
    			QLength y_Kill = y_norm*4.0*vessel_length;
    			VesselPtr<2> p_radiated_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        			Vertex<2>(x_Kill,y_Kill));
    			if (p_radiated_vessel->GetRadius()> 0.000000001_um)
			{
    				p_radiated_vessel->SetRadius(0.000000000001_um);
				KilledVessels++;
				p_viscosity_calculator->Calculate();
    				p_impedance_calculator->Calculate();
    				flow_solver.SetUp();
    				flow_solver.Solve();
				//std::cout << "The perfusion coefficient for " << KilledVessels << " killed vessels is: " << p_network->GetPerfusionQuotientBeta(threshold) << "\n";

    				outfile << alpha << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
    				QuotientMatrixAggregateForking[j_aux][KilledVessels] += p_network->GetPerfusionQuotientBeta(threshold);
    
    				if(p_network->GetPerfusionQuotientBeta(threshold)>optimumBeta)
    				{
    					optimumBeta=p_network->GetPerfusionQuotientBeta(threshold);
    					//optimumKilled = KilledVessels;
    				}
        		}

		}
    		
		//outfile2 << alpha << " " << optimumKilled << " \n"; 
    		p_network->Write(output_file_final_RT);
	}
   
    outfileMeanForking.open("testForkingMeanIntermediate100_HigherThreshold.txt");//, std::ios_base::app);
    //std::cout << "Attempts " << att+1 << " \n";
    outfileMeanForking << "Attempts " << att+1 << " \n";
    
    for (unsigned iii = 0; iii < HetVectorSize; iii++)
    {
        for(unsigned jjj = 0; jjj < NV_ToKill+1; jjj++)
        {
	    QuotientMatrixMeanForking[iii][jjj]= QuotientMatrixAggregateForking[iii][jjj]/double(att+1);
	    //std::cout << "For heterogeneity level " << iii << ", and radiotherapy dose " << jjj << ", the aggregate quotient is: " << QuotientMatrixAggregateForking[iii][jjj] << " , and the mean quotient is: " << QuotientMatrixMeanForking[iii][jjj] << endl;
	    outfileMeanForking << iii << " " << jjj << " " << QuotientMatrixMeanForking[iii][jjj] << " \n";
        }
    }    

    outfileMeanForking.close();


   }

}

void xTestRadiotherapyDichotomousRemoveMethods()
{
    //unsigned order=7;
    unsigned order=6;
    QLength max_radius = 50_um;
    double twicelambda = 8.0;
    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
    double inlet_haematocrit = 0.45;

    // Generate the network
    VesselNetworkGenerator<2> network_generator;
   
    QLength quarter_vertical_length = 0.9*twicelambda*max_radius*pow(2.0,-1.0/3.0);
    double dimless_length = 1.0;
    std::ofstream outfile;
    outfile.open("testRemoveMethodsForking.txt");//, std::ios_base::app);
    std::ofstream outfile2;
    outfile2.open("optimumRemoveMethodsForking.txt");//, std::ios_base::app);


    for(unsigned i_aux=1; i_aux<order+1; i_aux++)
    {
	dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
    }


    auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Dichotomous_Heterogeneous_RemoveMethods_Order6", true);
    //for(unsigned j_aux=0; j_aux<11; j_aux++)
    for(unsigned j_aux=0; j_aux<11; j_aux++)
    {
    	double alpha = 1.0+(double)j_aux*0.05;
    	std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, 			quarter_vertical_length, max_radius, alpha, twicelambda);
	std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");
	p_network->Write(output_file_initial);
    	std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetwork"+to_string(alpha)+".vtp");
    	std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkRT"+to_string(alpha)+".vtp");


    	std::vector<std::shared_ptr<Vessel<2> > > vessels = p_network->GetVessels();
    	typename std::vector<std::shared_ptr<Vessel<2> > >::iterator it;
    	for(it = vessels.begin(); it != vessels.end(); it++)
    	{
        	(*it)->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    	}

	QLength domain_side_length_x = dimless_length*2.0*twicelambda*max_radius;
    	VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
           Vertex<2>(0.0_um,2.0*quarter_vertical_length));
    	VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
           Vertex<2>(domain_side_length_x, 2.0*quarter_vertical_length));
    	p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
    	p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
    	p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
    	p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);

        auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
        p_viscosity_calculator->SetPlasmaViscosity(viscosity);
        p_viscosity_calculator->SetVesselNetwork(p_network);
        p_viscosity_calculator->Calculate();
        auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
    	p_impedance_calculator->SetVesselNetwork(p_network);
    	p_impedance_calculator->Calculate();
	FlowSolver<2> flow_solver;
    	flow_solver.SetVesselNetwork(p_network);
    	flow_solver.SetUp();
        flow_solver.SetUseDirectSolver(true);
    	flow_solver.Solve();
	QFlowRate threshold = 1.e-12*unit::metre_cubed_per_second;
	double PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
    	std::cout << "Perfusion quotient before RT is : " << PerfusionQuotient << " \n"; 
	p_network->Write(output_file_final);

	unsigned NV_Total = pow(2,order+2)-2;
    	unsigned NV_ToKillFrom = NV_Total;
    	unsigned KilledVessels = 0;
    	double dose = 0.5;
    	unsigned NV_ToKill = (unsigned)(dose*NV_ToKillFrom);
  //  quarter_vertical_length
    	outfile << alpha << " " << KilledVessels << " " << PerfusionQuotient << " \n"; 

    	double optimumBeta =PerfusionQuotient;
    	unsigned optimumKilled = 0;

	while( KilledVessels < NV_ToKill )
	{
    		auto minRadius = (*vessels.begin())->GetRadius();
    		unsigned IdToKill = 0;
		/*
    		for(it = vessels.begin(); it != vessels.end(); it++)
    		{	
			if  ((minRadius > (*it)->GetRadius())&&((*it)->GetRadius()>0.00000001_um))
			{
				minRadius = (*it)->GetRadius();
				IdToKill = (*it)->GetId();
			}
    		}
 		vessels[IdToKill]->SetRadius(0.000000000001_um);
    		vessels[IdToKill]->SetToDie();
		*/
    		for(it = vessels.begin(); it != vessels.end(); it++)
    		{	
			if  (minRadius > (*it)->GetRadius())
			{
				minRadius = (*it)->GetRadius();
				IdToKill = (*it)->GetId();
			}
    		}
 		p_network->RemoveVessel(vessels[IdToKill],true);
    		KilledVessels++;

    		p_viscosity_calculator->Calculate();
    		p_impedance_calculator->Calculate();
    		flow_solver.SetUp();
    		flow_solver.Solve();
    		std::cout << "The perfusion coefficient for " << KilledVessels << " killed vessels is: " << 					p_network->GetPerfusionQuotientBeta(threshold) << "\n";
		outfile << alpha << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
    
    		if(p_network->GetPerfusionQuotientBeta(threshold)>optimumBeta)
    		{
    			optimumBeta=p_network->GetPerfusionQuotientBeta(threshold);
    			optimumKilled = KilledVessels;
    		}

	}

    	outfile2 << alpha << " " << optimumKilled << " \n"; 

	p_network->Write(output_file_final_RT);
    }

}

// Same as TestRadiotherapyDichotomousRemoveMethodsUnifiedKillingConstantFlux except with non-unified killing
void xTestRadiotherapyDichotomousRemoveMethodsConstantFlux()
{
    // Define the key parameters
	//unsigned order=7;
    unsigned order=6;
    QLength max_radius = 50_um;
    double twicelambda = 8.0;
    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
    double inlet_haematocrit = 0.45;

    // Generate the network
    VesselNetworkGenerator<2> network_generator;
    QLength quarter_vertical_length = 0.9*twicelambda*max_radius*pow(2.0,-1.0/3.0);
    double dimless_length = 1.0;
	
	// Initialise the outputs
    std::ofstream outfile;
    outfile.open("testRemoveMethodsForkingConstFlux.txt");//, std::ios_base::app);
    std::ofstream outfile2;
    outfile2.open("optimumRemoveMethodsForkingConstFlux.txt");//, std::ios_base::app);
	
	// ??? Non-dimensionalising the length
    for(unsigned i_aux=1; i_aux<order+1; i_aux++)
    {
		dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
    }

	// Create a folder in which to store .vtp outputs  
    auto p_file_handler = std::make_shared<OutputFileHandler>("Const_Flux_Radiotherapy_Dichotomous_Heterogeneous_RemoveMethods_Order6", true);
    
	// Loop over various heterogeneities to compute optimal PQ
	//for(unsigned j_aux=0; j_aux<11; j_aux++)
    for(unsigned j_aux=0; j_aux<11; j_aux++)
    {
    	double alpha = 1.0+(double)j_aux*0.05;
    	std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, 			quarter_vertical_length, max_radius, alpha, twicelambda);
		std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");
		p_network->Write(output_file_initial);
    	std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetwork"+to_string(alpha)+".vtp");
    	std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkRT"+to_string(alpha)+".vtp");
    	
		// Set the haematocrit for all vessels
		std::vector<std::shared_ptr<Vessel<2> > > vessels = p_network->GetVessels();
    	typename std::vector<std::shared_ptr<Vessel<2> > >::iterator it;
		for(it = vessels.begin(); it != vessels.end(); it++)
    	{
        	(*it)->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    	}
		
		// Set input and output nodes and pressures
		QLength domain_side_length_x = dimless_length*2.0*twicelambda*max_radius;
    	VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
        Vertex<2>(0.0_um,2.0*quarter_vertical_length));  // ???
    	VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
        Vertex<2>(domain_side_length_x, 2.0*quarter_vertical_length));  // ???
    	p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
    	//p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
        p_inlet_node->GetFlowProperties()->SetUseVelocityBoundaryCondition(true);
        p_inlet_node->GetSegment(0)->GetFlowProperties()->SetFlowRate(2.e-10*unit::metre_cubed_per_second);
        p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
    	p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);
		
		// Set up a viscosity solver
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
    	flow_solver.SetUp();
        flow_solver.SetUseDirectSolver(true);
    	flow_solver.Solve();
		
		// Get pre-RT Perfusion Quotient
		QFlowRate threshold = 1.e-12*unit::metre_cubed_per_second;
		double PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
    	std::cout << "Perfusion quotient before RT is : " << PerfusionQuotient << " \n"; 
		p_network->Write(output_file_final);

		// Non-unified killing
		unsigned NV_Total = pow(2,order+2)-2;
    	unsigned NV_ToKillFrom = NV_Total;
    	unsigned KilledVessels = 0;
    	double dose = 0.45;
    	unsigned NV_ToKill = (unsigned)(dose*NV_ToKillFrom);
  		// quarter_vertical_length
    	outfile << alpha << " " << KilledVessels << " " << PerfusionQuotient << " \n"; 

		// Initialise the optimum parameters
    	double optimumBeta = PerfusionQuotient;
    	unsigned optimumKilled = 0;
		
		// Until a certain number of vessels have been killed, kill the smallest vessels
		while( KilledVessels < NV_ToKill )
		{
    		auto minRadius = (*vessels.begin())->GetRadius();  // ???vessels are in order of increasing size; get the smallest radius
    		unsigned IdToKill = 0;  // initialise the Vessel ID that will be killed
			
			// ???
			/*
    		for(it = vessels.begin(); it != vessels.end(); it++)
    		{	
				if  ((minRadius > (*it)->GetRadius())&&((*it)->GetRadius()>0.00000001_um))
				{
					minRadius = (*it)->GetRadius();
					IdToKill = (*it)->GetId();
				}
    		}
 			vessels[IdToKill]->SetRadius(0.000000000001_um);
    		vessels[IdToKill]->SetToDie();
			*/

			// ??? For all vessels in order of increasing size, find and remove the smallest vessel 
    		for(it = vessels.begin(); it != vessels.end(); it++)  
    		{	
				if  (minRadius > (*it)->GetRadius())  		  // if the vessel radius < radius threshold to kill
				{
					minRadius = (*it)->GetRadius();   		  // set the radius of the current vessel as the new radius threshold
					IdToKill = (*it)->GetId();  	  		  // store ID of vessel for killing
				}
    		}
 			p_network->RemoveVessel(vessels[IdToKill],true);  // remove the smallest vessel as detected previously 
    		KilledVessels++;  								  // increase the count of vessels removed
			
			// Compute the flow
    		p_viscosity_calculator->Calculate();
    		p_impedance_calculator->Calculate();
    		flow_solver.SetUp();
    		flow_solver.Solve();
    		
			// Print the results
			std::cout << "The perfusion coefficient for " << KilledVessels << " killed vessels is: " << 					p_network->GetPerfusionQuotientBeta(threshold) << "\n";
			outfile << alpha << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
    
    		
			// If the new PQ and kill-count is better than the previous optimum, store the new optimums
			if(p_network->GetPerfusionQuotientBeta(threshold)>optimumBeta)
    		{
    			optimumBeta=p_network->GetPerfusionQuotientBeta(threshold);
    			optimumKilled = KilledVessels;
    		}
		}
    	outfile2 << alpha << " " << optimumKilled << " \n"; 
		p_network->Write(output_file_final_RT);
    }
}

void xTestNoRadiotherapyDichotomousRemoveMethodsRandom()
{
    // Define the key parameters
	//unsigned order=7;
    srand(unsigned(time(NULL)));  // ??? pick a random variable
    unsigned order=6;
    unsigned NV_Total = pow(2,order+2)-2;  // total # of vessels in network
    double dose = 0.5;  // fraction to kill
    unsigned NV_ToKill = (unsigned)(dose*NV_Total);  // # of vessels to kill
    std::cout<< NV_ToKill << " \n";
    //unsigned NV_ToKillFrom = pow(2,7)*(pow(2,1)-1);
    //double dose = 0.5;
    //unsigned NV_ToKill = (unsigned)(dose*NV_ToKillFrom);
    //unsigned IdToKill;
    QLength max_radius = 50_um;
    double twicelambda = 8.0;
    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
    double inlet_haematocrit = 0.45;
    QFlowRate threshold = 1.e-12*unit::metre_cubed_per_second;

    // Generate the network
    VesselNetworkGenerator<2> network_generator;
    QLength quarter_vertical_length = 0.9*twicelambda*max_radius*pow(2.0,-1.0/3.0);
    double dimless_length = 1.0;
    double dimless_offset_x = 1.0;  // ???
	
	// ??? Non-dimensionalising the length
    for(unsigned i_aux=1; i_aux<order+1; i_aux++)
    {
		dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
    }

	// ??? not sure what we're doing here for the offset
    for(unsigned j_aux=1; j_aux<order; j_aux++)
    {
		dimless_offset_x += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(j_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(j_aux-1)));
    }
    QLength domain_side_length_x = dimless_length*2.0*twicelambda*max_radius;
    QLength offset_x = dimless_offset_x*twicelambda*max_radius;

    std::ofstream outfile;
    outfile.open("testRemoveMethodsRandomForking.txt");//, std::ios_base::app);
    std::ofstream outfile2;
    outfile2.open("optimumRemoveMethodsRandomForking.txt");//, std::ios_base::app);
    
	auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Dichotomous_Heterogeneous_RemoveMethods_Order6_Random", true);
    //for(unsigned j_aux=0; j_aux<11; j_aux++)
    //unsigned attempts = 100;
    unsigned attempts = 100;
    unsigned HetVectorSize = 11;
    //std::vector<std::vector<std::vector<double> > > QuotientMatrixSampling(NHet,std::vector<std::vector<double> >(ToBeKilled+1,std::vector <double>(attempts,0)));
    std::vector<std::vector<double> > QuotientMatrixMeanForking(HetVectorSize,std::vector<double>(NV_ToKill+1));
    std::ofstream outfileMeanForking;
    //std::ofstream outfile2Mean;
    std::vector<std::vector<double> > QuotientMatrixAggregateForking(HetVectorSize,std::vector<double>(NV_ToKill+1,0.0));

    for(unsigned att=0; att < attempts;att++)
    {
    	for(unsigned j_aux=0; j_aux<HetVectorSize; j_aux++)
    	//for(unsigned j_aux=0; j_aux<1; j_aux++)
		{
			
    		unsigned KilledVessels = 0;
    		double alpha = 1.0+(double)j_aux*0.05;
    		std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, 			quarter_vertical_length, max_radius, alpha, twicelambda);
			std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork"+to_string(alpha)+".vtp");
			p_network->Write(output_file_initial);
    		std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetwork"+to_string(alpha)+".vtp");
    		std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkRT"+to_string(alpha)+".vtp");

    		p_network->Write(output_file_initial);

    		std::vector<std::shared_ptr<Vessel<2> > > vessels = p_network->GetVessels();
    		typename std::vector<std::shared_ptr<Vessel<2> > >::iterator it;
    		for(it = vessels.begin(); it != vessels.end(); it++)
    		{
        		(*it)->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    		}

    		VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
            Vertex<2>(0.0_um,2.0*quarter_vertical_length));
   			VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
           	Vertex<2>(domain_side_length_x, 2.0*quarter_vertical_length));
    		p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
    		p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
    		p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
    		p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);
   
    		auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
    		p_viscosity_calculator->SetPlasmaViscosity(viscosity);
    		p_viscosity_calculator->SetVesselNetwork(p_network);
    		p_viscosity_calculator->Calculate();
    		auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
    		p_impedance_calculator->SetVesselNetwork(p_network);
    		p_impedance_calculator->Calculate();
    		FlowSolver<2> flow_solver;
    		flow_solver.SetVesselNetwork(p_network);
    		flow_solver.SetUp();
    		flow_solver.SetUseDirectSolver(true);
    		flow_solver.Solve();
			double PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
    		//std::cout << "Perfusion quotient before RT is : " << PerfusionQuotient << " \n";
    		p_network->Write(output_file_final);
		
    		QuotientMatrixAggregateForking[j_aux][KilledVessels] += p_network->GetPerfusionQuotientBeta(threshold);
    
    		outfile << alpha << " " << KilledVessels << " " << PerfusionQuotient << " \n"; 

    		double optimumBeta =PerfusionQuotient;
    		//unsigned optimumKilled = 0;
			while( KilledVessels < NV_ToKill )
			{
				double x_norm = (double)rand()/RAND_MAX;
    			double y_norm = (double)rand()/RAND_MAX;
    			QLength x_Kill = offset_x + x_norm*(domain_side_length_x-2.0*offset_x);
    			QLength y_Kill = y_norm*4.0*quarter_vertical_length;
    			VesselPtr<2> p_radiated_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        			Vertex<2>(x_Kill,y_Kill));
    			if (p_radiated_vessel->GetOwnerRank()==6)
				{
    				//p_radiated_vessel->SetRadius(0.000000000001_um);
					//IdToKill = p_radiated_vessel->GetId();
 					p_network->RemoveVessel(p_radiated_vessel,true);
					KilledVessels++;
					p_viscosity_calculator->Calculate();
    				p_impedance_calculator->Calculate();
    				flow_solver.SetUp();
    				flow_solver.Solve();
					//std::cout << "The perfusion coefficient for " << KilledVessels << " killed vessels is: " << p_network->GetPerfusionQuotientBeta(threshold) << "\n";

    				outfile << alpha << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
    				QuotientMatrixAggregateForking[j_aux][KilledVessels] += p_network->GetPerfusionQuotientBeta(threshold);
    
    				if(p_network->GetPerfusionQuotientBeta(threshold)>optimumBeta)
    				{
    					optimumBeta=p_network->GetPerfusionQuotientBeta(threshold);
    					//optimumKilled = KilledVessels;
    				}
        		
				vessels = p_network->GetVessels();
			}
		}
    		
		//outfile2 << alpha << " " << optimumKilled << " \n"; 
    		p_network->Write(output_file_final_RT);
	}
   
    outfileMeanForking.open("ForkingMeanIntermediateRandom100_6Gen.txt");//, std::ios_base::app);
    //std::cout << "Attempts " << att+1 << " \n";
    outfileMeanForking << "Attempts " << att+1 << " \n";
    
    for (unsigned iii = 0; iii < HetVectorSize; iii++)
    {
        for(unsigned jjj = 0; jjj < NV_ToKill+1; jjj++)
        {
	    QuotientMatrixMeanForking[iii][jjj]= QuotientMatrixAggregateForking[iii][jjj]/double(att+1);
	    //std::cout << "For heterogeneity level " << iii << ", and radiotherapy dose " << jjj << ", the aggregate quotient is: " << QuotientMatrixAggregateForking[iii][jjj] << " , and the mean quotient is: " << QuotientMatrixMeanForking[iii][jjj] << endl;
	    outfileMeanForking << iii << " " << jjj << " " << QuotientMatrixMeanForking[iii][jjj] << " \n";
        }
    }    
    outfileMeanForking.close();
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Voronoi Networks
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void xTestVoronoi_RemoveMethods_ConstantFlux()
{

    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
    double inlet_haematocrit = 0.45;
    VesselNetworkGenerator<2> network_generator;
    double dimless_domain_size_x = 2000.0; 
    QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
    auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Voronoi_Methods_ConstantFlux", true);

    std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");
    std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetwork.vtp");
    std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("PostRTNetwork.vtp");

    std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
    double tolerance = 0.001;  // 
    QFlowRate threshold = 1e-11*unit::metre_cubed_per_second;
    QFlowRate InletFlowRate = 1e-10*unit::metre_cubed_per_second;

    // Generate the network

    unsigned NumberOfAttempts = 10;
    unsigned NumberOfCutOffs = 10;
    std::ofstream outfileMean;
    unsigned NumberOfSeedPoints = 25; // This will be the most outer cycle...400 should be maximum


    std::vector<std::vector<double> > QuotientMatrixMean(1,std::vector<double>(NumberOfCutOffs+1));
    std::vector<std::vector<double> > QuotientMatrixAggregate(1,std::vector<double>(NumberOfCutOffs+1,0.0));
    
    for(unsigned attempt=1;attempt < NumberOfAttempts+1; attempt++)    // This will be the second most outer cycle
    {

    std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(NumberOfSeedPoints)+"SeedPoints/EdgesMatrixSampleNumber"+to_string(attempt)+".txt");
    std::vector<std::vector<double> > rEdgesMatrix;
    string line;

    while (std::getline(in, line)) 
    {
            rEdgesMatrix.push_back(std::vector<double>());
            // Break down the row into column values
            std::stringstream split(line);
            double value;
            while (split >> value)
		{
        	        rEdgesMatrix.back().push_back(value);
        	}
    }

    std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateVoronoiNetwork(rEdgesMatrix);

    // set inlet and outlet nodes

    std::vector<std::shared_ptr<Vessel<2> > > vessels = p_network->GetVessels();

    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
    {
	    (*vessel_iterator)->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);

            if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
            {
                    if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] <  tolerance)
                    {
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
			(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetUseVelocityBoundaryCondition(true);
                        //(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
			(*vessel_iterator)->GetFlowProperties()->SetFlowRate(InletFlowRate);
                    }
		    if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                    {
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
                        //(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                    }
	     }

             if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
	     {

		    if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] <  tolerance)
		    {
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
			(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetUseVelocityBoundaryCondition(true);
                        //(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
			(*vessel_iterator)->GetFlowProperties()->SetFlowRate(InletFlowRate);
                    }
		    if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
		    {
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
                        //(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                    }
            }
    }

    //auto p_segment = p_network->GetVesselSegments()[0];
    //p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    //VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);
    p_network->Write(output_file_initial);
    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
    p_viscosity_calculator->SetVesselNetwork(p_network);
    p_viscosity_calculator->Calculate();
    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
    p_impedance_calculator->SetVesselNetwork(p_network);
    p_impedance_calculator->Calculate();
    FlowSolver<2> flow_solver;
    flow_solver.SetVesselNetwork(p_network);
    flow_solver.SetUp();
    flow_solver.SetUseDirectSolver(false);
    flow_solver.Solve();
    //double PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
    p_network->Write(output_file_final);
    QuotientMatrixAggregate[0][0] += p_network->GetPerfusionQuotientBeta(threshold);
    
    //std::ofstream outfile;
    //outfile.open("VoronoiKilling.txt", std::ios_base::app);
    //outfile.open("VoronoiKilling.txt");
    //outfile << NumberOfSeedPoints << " " << 0.0 << " " << PerfusionQuotient << " \n"; 

    QLength length_cut_off;


    for(unsigned het = 1; het<NumberOfCutOffs+1;het++)
    {
    length_cut_off = 2.0*(double)het*unit::microns;
        
    	for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
        {   // if it's not an input or output vessel
            //if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() != 1 && (*vessel_iterator)->GetEndNode()->GetNumberOfSegments() != 1)
            if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() >2 && (*vessel_iterator)->GetEndNode()->GetNumberOfSegments() >2)  // This one only kills vessels with both ends having connectivity 3 (or more)
            {          // if the length is smaller than the cut off, kill it
                    if((*vessel_iterator)->GetLength() <  length_cut_off)
                    {
                        (*vessel_iterator)->SetRadius(0.0000000000001_um);
		        //p_network->RemoveVessel(*vessel_iterator,true);
                    }

	     }

	}




    //This one kills all vessels including inlets and outlets
    //p_network->RemoveShortVessels(length_cut_off, false);
    p_viscosity_calculator->SetVesselNetwork(p_network);
    p_viscosity_calculator->Calculate();
    p_impedance_calculator->SetVesselNetwork(p_network);
    p_impedance_calculator->Calculate();
    flow_solver.SetVesselNetwork(p_network);
    flow_solver.SetUp();
    flow_solver.SetUseDirectSolver(false);
    flow_solver.Solve();
    //PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
    //outfile << NumberOfSeedPoints << " " << 2.0*(double)het << " " << PerfusionQuotient << " \n"; 
    p_network->Write(output_file_final_RT);
    QuotientMatrixAggregate[0][het] += p_network->GetPerfusionQuotientBeta(threshold);
    }


    outfileMean.open("MeanUnperfusionVoronoi.txt");//, std::ios_base::app);
    std::cout<< "Attempts " << attempt << " \n";
    outfileMean << "Attempts " << attempt << " \n";
    
    for (unsigned iii = 0; iii < 1; iii++)
    {
        for(unsigned jjj = 0; jjj < NumberOfCutOffs+1; jjj++)
        {
	    QuotientMatrixMean[iii][jjj]= QuotientMatrixAggregate[iii][jjj]/double(attempt);
	    std::cout << "For density level " << iii << ", and cutoff level " << jjj << ", the aggregate quotient is: " << QuotientMatrixAggregate[iii][jjj] << " , and the mean quotient is: " << QuotientMatrixMean[iii][jjj] << endl;
	    outfileMean << iii << " " << jjj << " " << QuotientMatrixMean[iii][jjj] << " \n";
        }
    }    

    outfileMean.close();
    }
}

void xTestVoronoi_ScaleRadius()
{

    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
    double inlet_haematocrit = 0.45;
    VesselNetworkGenerator<2> network_generator;
    double dimless_domain_size_x = 2000.0; 
    QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
    auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Voronoi_Methods_ScaleRadius", true);

    std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");
    std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetwork.vtp");
    std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("PostRTNetwork.vtp");

    std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
    double tolerance = 0.001;
    QFlowRate threshold = 1e-12*unit::metre_cubed_per_second;

    // Generate the network

    unsigned NumberOfAttempts = 1;
    unsigned NumberOfCutOffs = 10;
    std::ofstream outfileMean;
    unsigned NumberOfSeedPoints = 25; // This will be the most outer cycle...400 should be maximum


    std::vector<std::vector<double> > QuotientMatrixMean(1,std::vector<double>(NumberOfCutOffs+1));
    std::vector<std::vector<double> > QuotientMatrixAggregate(1,std::vector<double>(NumberOfCutOffs+1,0.0));
    
    for(unsigned attempt=1;attempt < NumberOfAttempts+1; attempt++)    // This will be the second most outer cycle
    {

    std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(NumberOfSeedPoints)+"SeedPoints/EdgesMatrixSampleNumber"+to_string(attempt)+".txt");
    std::vector<std::vector<double> > rEdgesMatrix;
    string line;

    while (std::getline(in, line)) 
    {
            rEdgesMatrix.push_back(std::vector<double>());
            // Break down the row into column values
            std::stringstream split(line);
            double value;
            while (split >> value)
		{
        	        rEdgesMatrix.back().push_back(value);
        	}
    }

    std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateVoronoiNetworkScaleRadius(rEdgesMatrix);
    p_network->Write(output_file_initial);

    // set inlet and outlet nodes

    std::vector<std::shared_ptr<Vessel<2> > > vessels = p_network->GetVessels();

    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
    {
	    (*vessel_iterator)->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);

            if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
            {
                    if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] <  tolerance)
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

		    if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] <  tolerance)
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

    //auto p_segment = p_network->GetVesselSegments()[0];
    //p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    //VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);
    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
    p_viscosity_calculator->SetVesselNetwork(p_network);
    p_viscosity_calculator->Calculate();
    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
    p_impedance_calculator->SetVesselNetwork(p_network);
    p_impedance_calculator->Calculate();
    FlowSolver<2> flow_solver;
    flow_solver.SetVesselNetwork(p_network);
    flow_solver.SetUp();
    flow_solver.SetUseDirectSolver(false);
    flow_solver.Solve();
    //double PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
    p_network->Write(output_file_final);
    QuotientMatrixAggregate[0][0] += p_network->GetPerfusionQuotientBeta(threshold);
    
    //std::ofstream outfile;
    //outfile.open("VoronoiKilling.txt", std::ios_base::app);
    //outfile.open("VoronoiKilling.txt");
    //outfile << NumberOfSeedPoints << " " << 0.0 << " " << PerfusionQuotient << " \n"; 

    QLength length_cut_off;


    for(unsigned het = 1; het<NumberOfCutOffs+1;het++)
    {
    length_cut_off = 2.0*(double)het*unit::microns;
        
    	for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
        {   // if it's not an input or output vessel
            if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() != 1 && (*vessel_iterator)->GetEndNode()->GetNumberOfSegments() != 1)
            // if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() >2 && (*vessel_iterator)->GetEndNode()->GetNumberOfSegments() >2)  // This one only kills vessels with both ends having connectivity 3 (or more)
            {          // if the length is smaller than the cut off, kill it
                    if((*vessel_iterator)->GetLength() <  length_cut_off)
                    {
                        (*vessel_iterator)->SetRadius(0.0000000000001_um);
		        //p_network->RemoveVessel(*vessel_iterator,true);
                    }

	     }

	}




    //This one kills all vessels including inlets and outlets
    //p_network->RemoveShortVessels(length_cut_off, false);
    p_viscosity_calculator->SetVesselNetwork(p_network);
    p_viscosity_calculator->Calculate();
    p_impedance_calculator->SetVesselNetwork(p_network);
    p_impedance_calculator->Calculate();
    flow_solver.SetVesselNetwork(p_network);
    flow_solver.SetUp();
    flow_solver.SetUseDirectSolver(false);
    flow_solver.Solve();
    //PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
    //outfile << NumberOfSeedPoints << " " << 2.0*(double)het << " " << PerfusionQuotient << " \n"; 
    p_network->Write(output_file_final_RT);
    QuotientMatrixAggregate[0][het] += p_network->GetPerfusionQuotientBeta(threshold);
    }


    outfileMean.open("MeanUnperfusionVoronoi.txt");//, std::ios_base::app);
    std::cout<< "Attempts " << attempt << " \n";
    outfileMean << "Attempts " << attempt << " \n";
    
    for (unsigned iii = 0; iii < 1; iii++)
    {
        for(unsigned jjj = 0; jjj < NumberOfCutOffs+1; jjj++)
        {
	    QuotientMatrixMean[iii][jjj]= QuotientMatrixAggregate[iii][jjj]/double(attempt);
	    std::cout << "For density level " << iii << ", and cutoff level " << jjj << ", the aggregate quotient is: " << QuotientMatrixAggregate[iii][jjj] << " , and the mean quotient is: " << QuotientMatrixMean[iii][jjj] << endl;
	    outfileMean << iii << " " << jjj << " " << QuotientMatrixMean[iii][jjj] << " \n";
        }
    }    

    outfileMean.close();
    }
}

// ??? not sure what this does
/*
   p_network->RemoveShortVessels(length_cut_off, false);
    p_network->Write(output_file_final);
std::cout << p_network->GetPerfusionQuotientBeta(threshold);
    length_cut_off=40.0*unit::microns;
    p_network->RemoveShortVessels(length_cut_off,false);
    p_network->Write(output_file_final_RT);
std::cout << p_network->GetPerfusionQuotientBeta(threshold);
*/

// The following two tests are used to generate Fig. 10 and the data is used for Fig. 11

// Test to run Voronoi simulations with constant pressure difference (used for figures)
void xTestVoronoi_RemoveMethods_ConstPressure()
{
	// Define the key parameters
    unsigned NumberOfSeedPoints = 400;  // change this to select which Voronoi architecture to use: 25, 100, 400
    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
    double inlet_haematocrit = 0.45;

    // Initialise the simulation space
    VesselNetworkGenerator<2> network_generator;
    double dimless_domain_size_x = 2000.0; 
    QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
    
	// Initialise the .vtp outputs
	auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Voronoi_WithMethods_ConstPressure_SeedPoints"+to_string(NumberOfSeedPoints), true);  // create a folder
    std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");
    std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetwork.vtp");
    std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("PostRTNetwork.vtp");

	// Set the simulation paramaters
    std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
    double tolerance = 0.001;  // tolerance for solution convergence
    QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;  // minimum flow
    unsigned NumberOfAttempts = 100;  // number of times to run simulation (??? but why, is something stochastic?)
    std::ofstream outfileMean;
    unsigned NumberOfCutOffs = 10;  // upper limit for varying pruning length threshold
    std::vector<std::vector<double> > QuotientMatrixMean(1,std::vector<double>(NumberOfCutOffs+1));
    std::vector<std::vector<double> > QuotientMatrixAggregate(1,std::vector<double>(NumberOfCutOffs+1,0.0));
    
	// ??? loop over the number of random selections of seed points generating the Voronoi diagrams <what is attempts? What is stochastic here?>
    for(unsigned attempt=1;attempt < NumberOfAttempts+1; attempt++)    
    {
		// Read the network from a file
		std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(NumberOfSeedPoints)+"SeedPoints/EdgesMatrixSampleNumber"+to_string(attempt)+".txt");
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

		// Set the haematocrit for all vessels
		auto p_segment = p_network->GetVesselSegments()[0];
		p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
		VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);
		
		// Set up the viscosity solver
		auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
		p_viscosity_calculator->SetPlasmaViscosity(viscosity);
		p_viscosity_calculator->SetVesselNetwork(p_network);
		p_viscosity_calculator->Calculate();
		
		// Set up the impedance solver
		auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
		p_impedance_calculator->SetVesselNetwork(p_network);
		p_impedance_calculator->Calculate();
		if(attempt==1)
    	{
    		p_network->Write(output_file_initial);
    	}
		
		// Set up the flow solver 
		FlowSolver<2> flow_solver;
		flow_solver.SetVesselNetwork(p_network);
		flow_solver.SetUp();
		flow_solver.SetUseDirectSolver(true);
		flow_solver.Solve();

		// Get pre-RT Perfusion Quotient
		double PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
		if(attempt==1)
    	{
    		p_network->Write(output_file_final);
    	}
    	QuotientMatrixAggregate[0][0] += p_network->GetPerfusionQuotientBeta(threshold);

		// Initialise the output PQ file
		std::ofstream outfile;
		//outfile.open("VoronoiKilling.txt", std::ios_base::app);
		outfile.open("VoronoiKillingConstPressureDrop.txt");
		outfile << NumberOfSeedPoints << " " << 0.0 << " " << PerfusionQuotient << " \n"; 

		// Loop over increasing pruning length thresholds upto max. cut-off 
   		for(unsigned het = 1; het<NumberOfCutOffs+1;het++)
    	{
    		QLength length_cut_off(3.0*(double)het*unit::microns);  // cut-off length threshold for pruning

			// This would also kill short inlet and outlet vessels
        	//p_network->RemoveShortVessels(length_cut_off, false);

    		// This kills all vessels below certain length
			vessels = p_network->GetVessels();  // ??? This update is needed, if we use RemoveVessel instead of SetRadius
    		for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
        	{   
				// If it's not an input or output vessel
            	if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() != 1 && (*vessel_iterator)->GetEndNode()->GetNumberOfSegments() != 1)
            	{   
					// If the length is smaller than the cut off, kill it
                    if((*vessel_iterator)->GetLength() <  length_cut_off)
                    {
                        //(*vessel_iterator)->SetRadius(0.0000000000001_um);
		        		p_network->RemoveVessel(*vessel_iterator,true);
                    }
	     		}
			}

			// Calculate new flow
			//p_network->Write(output_file_final);
			p_viscosity_calculator->SetVesselNetwork(p_network);
			p_viscosity_calculator->Calculate();
			p_impedance_calculator->SetVesselNetwork(p_network);
			p_impedance_calculator->Calculate();
			flow_solver.SetVesselNetwork(p_network);
			flow_solver.SetUp();
			flow_solver.SetUseDirectSolver(true);
			flow_solver.Solve();
			PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
			
			// Save output for current attempt in VoronoiKillingConstPressureDrop.txt
			outfile << NumberOfSeedPoints << " " << 3.0*(double)het << " " << PerfusionQuotient << " \n"; 
			//p_network->Write(output_file_final_RT);
			QuotientMatrixAggregate[0][het] += p_network->GetPerfusionQuotientBeta(threshold);
    	}

    	///???
		if(attempt==1)
		{
			p_network->Write(output_file_final_RT);
		}

		// Save PQ for average of all attempts (used for Fig. 11) 
		outfileMean.open("MeanPerfusionVoronoi"+to_string(NumberOfSeedPoints)+"SeedPointsConstPressDrop.txt");//, std::ios_base::app);
		std::cout<< "Attempts " << attempt << " \n";
		outfileMean << "Attempts " << attempt << " \n";		
		for (unsigned iii = 0; iii < 1; iii++)
		{
			for(unsigned jjj = 0; jjj < NumberOfCutOffs+1; jjj++)
			{
				QuotientMatrixMean[iii][jjj]= QuotientMatrixAggregate[iii][jjj]/double(attempt);
				std::cout << "For density level " << iii << ", and cutoff level " << jjj << ", the aggregate quotient is: " << QuotientMatrixAggregate[iii][jjj] << " , and the mean quotient is: " << QuotientMatrixMean[iii][jjj] << endl;
				outfileMean << to_string(NumberOfSeedPoints) << " " << 3*jjj << " " << QuotientMatrixMean[iii][jjj] << " \n";
			}
		}    
    	outfileMean.close();
    }
}

// Test to run Voronoi simulations with constant flow rate (used for figures)
void xTestVoronoi_RemoveMethods_ConstFlowRate()
{

    //bool TooCloseToInlet;

    unsigned NumberOfSeedPoints = 400; // This will be the most outer cycle
    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
    double inlet_haematocrit = 0.45;
    VesselNetworkGenerator<2> network_generator;
    double dimless_domain_size_x = 2000.0; 
    QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
    auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Voronoi_WithMethods_ConstFlowRate_SeedPoints"+to_string(NumberOfSeedPoints), true);

    std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");
    std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetwork.vtp");
    std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("PostRTNetwork.vtp");

    std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
  //  std::vector<std::shared_ptr<Vessel<2> > >::iterator neighbours_iterator;
    double tolerance = 0.001;
    QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;

    // Generate the network

    unsigned NumberOfAttempts = 100;

    std::ofstream outfileMean;

    unsigned NumberOfCutOffs = 10;


    std::vector<std::vector<double> > QuotientMatrixMean(1,std::vector<double>(NumberOfCutOffs+1));
    
    std::vector<std::vector<double> > QuotientMatrixAggregate(1,std::vector<double>(NumberOfCutOffs+1,0.0));
    

    for(unsigned attempt=1;attempt < NumberOfAttempts+1; attempt++)    // This will be the second most outer cycle
    {

    std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(NumberOfSeedPoints)+"SeedPoints/EdgesMatrixSampleNumber"+to_string(attempt)+".txt");

    std::vector<std::vector<double> > rEdgesMatrix;
    string line;

    while (std::getline(in, line)) 
    {
            rEdgesMatrix.push_back(std::vector<double>());
            // Break down the row into column values
            std::stringstream split(line);
            double value;
            while (split >> value)
		{
        	        rEdgesMatrix.back().push_back(value);
        	}
    }

    std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateVoronoiNetwork(rEdgesMatrix);

    // set inlet and outlet nodes

    std::vector<std::shared_ptr<Vessel<2> > > vessels = p_network->GetVessels();

//    std::vector<std::shared_ptr<Vessel<2> > > LocalNeighbours;

    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
    {
            if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
            {
                    if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] <  tolerance)
                    {
                        //(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
                        //(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
			(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
        		(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetUseVelocityBoundaryCondition(true);
        		(*vessel_iterator)->GetFlowProperties()->SetFlowRate(8.e-13*unit::metre_cubed_per_second);
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
                        //(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
                        //(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
			(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
        		(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetUseVelocityBoundaryCondition(true);
        		(*vessel_iterator)->GetFlowProperties()->SetFlowRate(8.e-13*unit::metre_cubed_per_second);
                    }
		    //if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
		    else
		    {
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                    }
            }
    }
    p_network->Write(output_file_initial);
    auto p_segment = p_network->GetVesselSegments()[0];
    p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);
    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
    p_viscosity_calculator->SetVesselNetwork(p_network);
    p_viscosity_calculator->Calculate();
    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
    p_impedance_calculator->SetVesselNetwork(p_network);
    p_impedance_calculator->Calculate();
    FlowSolver<2> flow_solver;
    flow_solver.SetVesselNetwork(p_network);
    flow_solver.SetUp();
    flow_solver.SetUseDirectSolver(true);
    flow_solver.Solve();
    double PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
    p_network->Write(output_file_final);

    QuotientMatrixAggregate[0][0] += p_network->GetPerfusionQuotientBeta(threshold);

    
    std::ofstream outfile;
    //outfile.open("VoronoiKilling.txt", std::ios_base::app);
    outfile.open("VoronoiKillingConstFlowRate.txt");
    outfile << NumberOfSeedPoints << " " << 0.0 << " " << PerfusionQuotient << " \n"; 

    for(unsigned het = 1; het<NumberOfCutOffs+1;het++)
    {
    	QLength length_cut_off(3*(double)het*unit::microns);

	//This would also kill short inlet and outlet vessels

        //p_network->RemoveShortVessels(length_cut_off, false);

    	// This kills all vessels below certain length
        
	vessels = p_network->GetVessels();  // This update is needed, if we use RemoveVessel instead of SetRadius
        
    	for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
        {   // if it's not an input or output vessel
            //if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() != 1 && (*vessel_iterator)->GetEndNode()->GetNumberOfSegments() != 1)
	    if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > 1000 && (*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > 1000)
            {          // if the length is smaller than the cut off, kill it
		    //TooCloseToInlet = false;
		    //LocalNeighbours = (*vessel_iterator)->GetConnectedVessels();

                    //for (neighbours_iterator = LocalNeighbours.begin(); neighbours_iterator != LocalNeighbours.end(); neighbours_iterator++)
        	    //{
		//	if((*neighbours_iterator)->GetStartNode()->GetFlowProperties()->IsInputNode()==true || (*neighbours_iterator)->GetEndNode()->GetFlowProperties()->IsInputNode()==true)
		//	{
		//		TooCloseToInlet = true;
		//	}
		  //  }   
                    if((*vessel_iterator)->GetLength() <  length_cut_off)// && TooCloseToInlet==false)
                    {
                        //(*vessel_iterator)->SetRadius(0.0000000000001_um);
		        p_network->RemoveVessel(*vessel_iterator,true);
                    }

	     }

	}

    //p_network->Write(output_file_final);
    p_viscosity_calculator->SetVesselNetwork(p_network);
    p_viscosity_calculator->Calculate();
    p_impedance_calculator->SetVesselNetwork(p_network);
    p_impedance_calculator->Calculate();
    flow_solver.SetVesselNetwork(p_network);
    flow_solver.SetUp();
    flow_solver.SetUseDirectSolver(true);
    flow_solver.Solve();
    PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
    outfile << NumberOfSeedPoints << " " << 3*(double)het << " " << PerfusionQuotient << " \n"; 
    //p_network->Write(output_file_final_RT);
    QuotientMatrixAggregate[0][het] += p_network->GetPerfusionQuotientBeta(threshold);
        p_network->Write(output_file_final_RT);
    }

    outfileMean.open("MeanPerfusionVoronoi"+to_string(NumberOfSeedPoints)+"SeedPointsConstFlowRate.txt");//, std::ios_base::app);
    std::cout<< "Attempts " << attempt << " \n";
    outfileMean << "Attempts " << attempt << " \n";
    
    for (unsigned iii = 0; iii < 1; iii++)
    {
        for(unsigned jjj = 0; jjj < NumberOfCutOffs+1; jjj++)
        {
	    QuotientMatrixMean[iii][jjj]= QuotientMatrixAggregate[iii][jjj]/double(attempt);
	    std::cout << "For density level " << iii << ", and cutoff level " << jjj << ", the aggregate quotient is: " << QuotientMatrixAggregate[iii][jjj] << " , and the mean quotient is: " << QuotientMatrixMean[iii][jjj] << endl;
	    outfileMean << to_string(NumberOfSeedPoints) << " " << 3*jjj << " " << QuotientMatrixMean[iii][jjj] << " \n";
        }
    }    


    outfileMean.close();


    }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hexagonal Networks (generated using Voronoi networks with regular spacing of seed points)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// The following two tests are used to generate Fig. 8 and the data is used for Fig. 9

// Test to run hexagonal simulations with constant pressure difference (used for figures)
void xTestHexagonalFromVoronoi()
{
	// Define the key parameters
	QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
   	double inlet_haematocrit = 0.45;
   	double dimless_domain_size_x = 2000.0;
	unsigned dimless_vessel_length = 100.0;
	double tolerance = 0.001;
	QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;
	std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
	//unsigned NV_Total=6*units_in_x_direction*units_in_y_direction+units_in_y_direction+units_in_x_direction+1;
	//unsigned NV_ToThinFrom = NV_Total - 2*(3*units_in_y_direction + 1);
	//unsigned MaxHeterogeneity = (unsigned)(NV_ToThinFrom*0.05);
	//double percToKill =0.05;
	//unsigned ToBeKilled = (unsigned)(percToKill*NV_ToThinFrom);
	unsigned NV_ToThinFrom = 386;  // number of vessels from which to select ones to make thin (??? where did this number come from?)
	double percToKill = 0.2;  // percentage of vessels to kill
	unsigned ToBeKilled = (unsigned)(percToKill*NV_ToThinFrom);  // number to kill
	//unsigned MaxHeterogeneity = (unsigned)(NV_ToThinFrom*0.05);

	// Read the network from a file
    VesselNetworkGenerator<2> network_generator;
	std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix.txt");
	
	// Initialise the outputs
   	auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Hexagonal_FromVoronoi", true);  // create a folder
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

	// Set the simulation paramaters
    std::shared_ptr<VesselNetwork<2> > p_network;
	std::vector<std::shared_ptr<Vessel<2> > > vessels;
   	QLength large_vessel_radius = 10_um;
    QLength small_vessel_radius = 5_um;
	double PerfusionQuotient;
	unsigned ToBeThin;
	unsigned NThin;
	//unsigned NumberOfInlets = 11;
    unsigned NHet = 5;
	unsigned attempts = 100;  // number of times to run the simulation
	std::vector<std::vector<double> > QuotientMatrixMean(NHet,std::vector<double>(ToBeKilled+1));
	std::ofstream outfileMean;
    std::vector<std::vector<double> > QuotientMatrixAggregate(NHet,std::vector<double>(ToBeKilled+1,0.0));
    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
	auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
	FlowSolver<2> flow_solver;
	flow_solver.SetUseDirectSolver(true);
    string line2;
	std::vector<std::vector<unsigned> > Order;
	VesselPtr<2> p_thin_vessel;

	// loop over 100 random realizations of which vessels should be thin - these selections are read from text files
	for(unsigned att=0; att < attempts;att++)
	{
		outfileMean.open("HexagonalFromVoronoiMean.txt");//, std::ios_base::app);
		std::ifstream in2("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/Selection"+to_string(att+1)+".txt");
		Order.clear();
		while (std::getline(in2, line2)) 
		{
					Order.push_back(std::vector<unsigned>());
					// Break down the row into column values
				std::stringstream split2(line2);
				unsigned value2;
					while (split2 >> value2)
			{
						Order.back().push_back(value2);
				}
		}

		// loop over increasing level of heteregeneity (percentage of thin vessels)

		for(unsigned i=0; i<NHet;i++)
		{
			std::string output_file_two_radii = p_file_handler->GetOutputDirectoryFullPath().append("TwoRadiiNetwork"+to_string(i)+".vtp");
			std::string output_file_RT = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkRT"+to_string(i)+".vtp");

			p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);
			auto p_segment = p_network->GetVesselSegments()[0];
			// set inlet and outlet nodes
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
			p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
			VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
			double percOfThin = 0.05*i;
   			ToBeThin = (unsigned)(percOfThin*NV_ToThinFrom);
   			NThin = 0;

			// make a certain percentage of vessels thin
   			while( NThin < ToBeThin)
			{
				p_thin_vessel = p_network->GetVessel(Order[NThin][0]-1);
				p_thin_vessel->SetRadius(small_vessel_radius);
				NThin++;
   			}
   			//p_network->Write(output_file_two_radii);
			vessels = p_network->GetVessels();

        	// Kill vessels from smallest, up to the specified dose
    		for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
			{
				//vessels = p_network->GetVessels();
				p_viscosity_calculator->SetVesselNetwork(p_network);
    			p_impedance_calculator->SetVesselNetwork(p_network);
    			flow_solver.SetVesselNetwork(p_network);
		    	p_viscosity_calculator->Calculate();
		    	p_impedance_calculator->Calculate();
		    	flow_solver.SetUp();
				flow_solver.Solve();
				unsigned PictureNumber = (unsigned)(0.5*ToBeKilled);
				if(KilledVessels ==PictureNumber)
				{
					p_network->Write(output_file_RT);  // write the output file halfway through pruning
				}
				if(KilledVessels ==0)
				{
					p_network->Write(output_file_two_radii);  // write the output file before of pruning
				}
				PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
        		QuotientMatrixAggregate[i][KilledVessels] += PerfusionQuotient;  // store the PQ with the associated NHet and number of kills
				p_network->RemoveVessel(vessels[Order[KilledVessels][0]-1],true);  // remove the vessel
				//p_network->RemoveVessel(p_network->GetVessel(Order[KilledVessels][0]-1),true);
			}
    	}    
    	std::cout << "Attempts " << att+1 << " \n";
    	outfileMean << "Attempts " << att+1 << " \n";

		// write the resulting perfusion quotients into an output file
    
    	for (unsigned iii = 0; iii < NHet; iii++)
    	{
        	for(unsigned jjj = 0; jjj < ToBeKilled+1; jjj++)
        	{
	    		QuotientMatrixMean[iii][jjj]= QuotientMatrixAggregate[iii][jjj]/double(att+1);
	   	 		//std::cout << "For heterogeneity level " << iii << ", and radiotherapy dose " << jjj << ", the aggregate quotient is: " << QuotientMatrixAggregate[iii][jjj] << " , and the mean quotient is: " << QuotientMatrixMean[iii][jjj] << endl;
	    		outfileMean << iii << " " << jjj << " " << QuotientMatrixMean[iii][jjj] << " \n";
        	}
    	}    
		outfileMean.close();
	}
}

// Test to run hexagonal simulations with constant flow rate (used for figures)
void xTestHexagonalFromVoronoiConstantFlux()
{
	
	QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
   	double inlet_haematocrit = 0.45;
   	double dimless_domain_size_x = 2000.0;
    	unsigned dimless_vessel_length = 100.0;
    	double tolerance = 0.001;
    	QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;
    	std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
    	//unsigned NV_Total=6*units_in_x_direction*units_in_y_direction+units_in_y_direction+units_in_x_direction+1;
    	//unsigned NV_ToThinFrom = NV_Total - 2*(3*units_in_y_direction + 1);
    	//unsigned MaxHeterogeneity = (unsigned)(NV_ToThinFrom*0.05);
    	//double percToKill =0.05;
    	//unsigned ToBeKilled = (unsigned)(percToKill*NV_ToThinFrom);
	//unsigned NV_ToThinFrom =386;
	unsigned NV_ToThinFrom =346;
    	double percToKill =0.2;
    	unsigned ToBeKilled = (unsigned)(percToKill*NV_ToThinFrom);
	//unsigned MaxHeterogeneity = (unsigned)(NV_ToThinFrom*0.05);

    	// Generate the network

    	VesselNetworkGenerator<2> network_generator;
	std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrixConstFlowRate.txt");

   	auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Hexagonal_FromVoronoi_ConstantFlux", true);

    	std::vector<std::vector<double> > rEdgesMatrix;
    	string line;

    	while (std::getline(in, line)) 
   	{
            	rEdgesMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(line);
           	double value;
            	while (split >> value)
		{
        	        rEdgesMatrix.back().push_back(value);
        	}
    	}

    	std::shared_ptr<VesselNetwork<2> > p_network;
	std::vector<std::shared_ptr<Vessel<2> > > vessels;


   	QLength large_vessel_radius = 10_um;
    	QLength small_vessel_radius = 5_um;
	double PerfusionQuotient;



	unsigned ToBeThin;
	unsigned NThin;
	//unsigned NumberOfInlets = 11;
    	unsigned NHet = 5;
	unsigned attempts = 100;
	std::vector<std::vector<double> > QuotientMatrixMean(NHet,std::vector<double>(ToBeKilled+1));
	std::ofstream outfileMean;
    	std::vector<std::vector<double> > QuotientMatrixAggregate(NHet,std::vector<double>(ToBeKilled+1,0.0));

    	auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
    	p_viscosity_calculator->SetPlasmaViscosity(viscosity);
	auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
	FlowSolver<2> flow_solver;
	flow_solver.SetUseDirectSolver(true);
    	string line2;
	std::vector<std::vector<unsigned> > Order;
	VesselPtr<2> p_thin_vessel;
	
    	for(unsigned att=0; att < attempts;att++)
    	{
    		outfileMean.open("HexagonalFromVoronoiMeanConstantFlux.txt");//, std::ios_base::app);
		std::ifstream in2("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/ConstFlowRateSelection"+to_string(att+1)+".txt");
		Order.clear();
		while (std::getline(in2, line2)) 
   		{
            		Order.push_back(std::vector<unsigned>());
            		// Break down the row into column values
         		std::stringstream split2(line2);
           		unsigned value2;
            		while (split2 >> value2)
			{
        	        	Order.back().push_back(value2);
        		}
    		}
		//std::cout << Order[0][0]<<"  " << Order[1][0] << "  " << Order[192][0] <<"\n";// << Order[1][0]<<"\n";
		
    		for(unsigned i=0; i<NHet;i++)
    		{
 
			std::string output_file_two_radii = p_file_handler->GetOutputDirectoryFullPath().append("TwoRadiiNetwork"+to_string(i)+".vtp");
    			std::string output_file_RT = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkRT"+to_string(i)+".vtp");

			p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);
			auto p_segment = p_network->GetVesselSegments()[0];
		    	// set inlet and outlet nodes
	        	vessels = p_network->GetVessels();

	    		for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
	    		{
	        	    	if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
	            		{
	                	    	if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < -0.5*dimless_vessel_length+tolerance)
	                	    	{
						(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
	                	        	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetUseVelocityBoundaryCondition(true);
	                	        	(*vessel_iterator)->GetStartNode()->GetSegment(0)->GetFlowProperties()->SetFlowRate(8.e-13*unit::metre_cubed_per_second);
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
	                	        	(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetUseVelocityBoundaryCondition(true);
	                	        	(*vessel_iterator)->GetEndNode()->GetSegment(0)->GetFlowProperties()->SetFlowRate(8.e-13*unit::metre_cubed_per_second);
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
	    		p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
   			VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
			double percOfThin = 0.05*i;
   			ToBeThin = (unsigned)(percOfThin*NV_ToThinFrom);
   			NThin = 0;

   			while( NThin < ToBeThin)
			{
				p_thin_vessel = p_network->GetVessel(Order[NThin][0]-1);
				p_thin_vessel->SetRadius(small_vessel_radius);
				NThin++;
   			}
   			//p_network->Write(output_file_two_radii);
			vessels = p_network->GetVessels();

        		// Kill vessels from smallest, up to the specified dose
    			for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
			{
				//vessels = p_network->GetVessels();
				p_viscosity_calculator->SetVesselNetwork(p_network);
    				p_impedance_calculator->SetVesselNetwork(p_network);
    				flow_solver.SetVesselNetwork(p_network);
		    		p_viscosity_calculator->Calculate();
		    		p_impedance_calculator->Calculate();
		    		flow_solver.SetUp();
		    		flow_solver.Solve();
				unsigned PictureNumber = (unsigned)(0.5*ToBeKilled);
				if(KilledVessels ==PictureNumber)
				{
				p_network->Write(output_file_RT);
				}
				if(KilledVessels ==0)
				{
				p_network->Write(output_file_two_radii);
				}
				PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
        			QuotientMatrixAggregate[i][KilledVessels] += PerfusionQuotient;
				p_network->RemoveVessel(vessels[Order[KilledVessels][0]-1],true);
				//p_network->RemoveVessel(p_network->GetVessel(Order[KilledVessels][0]-1),true);
    			
			}
			

    		}    


    		std::cout << "Attempts " << att+1 << " \n";
    		outfileMean << "Attempts " << att+1 << " \n";
    
    		for (unsigned iii = 0; iii < NHet; iii++)
    		{
        		for(unsigned jjj = 0; jjj < ToBeKilled+1; jjj++)
        		{
	    		QuotientMatrixMean[iii][jjj]= QuotientMatrixAggregate[iii][jjj]/double(att+1);
	   	 	std::cout << "For heterogeneity level " << iii << ", and radiotherapy dose " << jjj << ", the aggregate quotient is: " << QuotientMatrixAggregate[iii][jjj] << " , and the mean quotient is: " << QuotientMatrixMean[iii][jjj] << endl;
	    		outfileMean << iii << " " << jjj << " " << QuotientMatrixMean[iii][jjj] << " \n";
        		}
    		}    

	outfileMean.close();
	}


}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Biological Networks (+ some hexagonal stuff)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void xTestRealNetworkNoRT()
{
	
	unsigned inletIndex = 79;
	QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
   	double inlet_haematocrit = 0.45;
	std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/RealNetwork/NodeCoordinates.txt");

    	std::vector<std::vector<double> > rCoordinatesMatrix;
    	string line;

    	while (std::getline(in, line)) 
   	{
            	rCoordinatesMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(line);
           	double value;
            	while (split >> value)
		{
        	        rCoordinatesMatrix.back().push_back(value);
        	}
    	}


	std::ifstream inDia("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/RealNetwork/AdjacencyDiameters.txt");

    	std::vector<std::vector<double> > rDiametersMatrix;
    	string lineDia;

    	while (std::getline(inDia, lineDia)) 
   	{
            	rDiametersMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(lineDia);
           	double valueDia;
            	while (split >> valueDia)
		{
        	        rDiametersMatrix.back().push_back(valueDia);
        	}
    	}

	std::ifstream inLen("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/RealNetwork/AdjacencyLengths.txt");

    	std::vector<std::vector<double> > rLengthsMatrix;
    	string lineLen;

    	while (std::getline(inLen, lineLen)) 
   	{
            	rLengthsMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(lineLen);
           	double valueLen;
            	while (split >> valueLen)
		{
        	        rLengthsMatrix.back().push_back(valueLen);
        	}
    	}


    
    	// Generate the network

    	VesselNetworkGenerator<3> network_generator;
    	std::shared_ptr<VesselNetwork<3> > p_network;

   	auto p_file_handler = std::make_shared<OutputFileHandler>("RealNetwork_NoRT", true);
    	auto p_viscosity_calculator = ViscosityCalculator<3>::Create();
    	p_viscosity_calculator->SetPlasmaViscosity(viscosity);
	auto p_impedance_calculator = VesselImpedanceCalculator<3>::Create();
	FlowSolver<3> flow_solver;
	//flow_solver.SetUseDirectSolver(true);

	std::string output_file_real_network = p_file_handler->GetOutputDirectoryFullPath().append("RealNetwork.vtp");
	p_network = network_generator.GenerateRealNetworkFromMatrices3D(rCoordinatesMatrix, rDiametersMatrix);

	std::vector<std::shared_ptr<Vessel<3> > > vessels = p_network->GetVessels();


	// set inlet and outlet nodes
	std::vector<std::shared_ptr<Vessel<3> > >::iterator vessel_iterator;
	//typename std::vector<std::shared_ptr<Vessel<3> > >::iterator it;


	for (vessel_iterator = vessels.begin();vessel_iterator != vessels.end(); vessel_iterator++)
	{
		(*vessel_iterator)->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);

	       	if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
	        {
	               	if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] == rCoordinatesMatrix[inletIndex][0]  && (*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[1] ==  rCoordinatesMatrix[inletIndex][1] && (*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[2] ==  rCoordinatesMatrix[inletIndex][2])
	                {
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
	                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
	                }
			else
	                {
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
	                }
		}
	
	        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)   	
				//std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
		{
			if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] == rCoordinatesMatrix[inletIndex][0]  && (*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[1] ==  rCoordinatesMatrix[inletIndex][1] && (*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[2] ==  rCoordinatesMatrix[inletIndex][2])
			{
	  			(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
	  	        }
			else
			{
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
	  	        }
	  	}
	}

	p_viscosity_calculator->SetVesselNetwork(p_network);
    	p_impedance_calculator->SetVesselNetwork(p_network);
    	flow_solver.SetVesselNetwork(p_network);
	p_viscosity_calculator->Calculate();
	p_impedance_calculator->CalculateFromMatrix(rLengthsMatrix);
	flow_solver.SetUp();
	flow_solver.Solve();
	p_network->Write(output_file_real_network);


}

void xTestRealNetworkRT()
{
	
	unsigned inletIndex = 79;
	QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
   	double inlet_haematocrit = 0.45;
    	QFlowRate threshold = 1.e-15*unit::metre_cubed_per_second;
	std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/RealNetwork/NodeCoordinates.txt");

    	std::vector<std::vector<double> > rCoordinatesMatrix;
    	string line;

    	while (std::getline(in, line)) 
   	{
            	rCoordinatesMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(line);
           	double value;
            	while (split >> value)
		{
        	        rCoordinatesMatrix.back().push_back(value);
        	}
    	}


	std::ifstream inDia("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/RealNetwork/AdjacencyDiameters.txt");

    	std::vector<std::vector<double> > rDiametersMatrix;
    	string lineDia;

    	while (std::getline(inDia, lineDia)) 
   	{
            	rDiametersMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(lineDia);
           	double valueDia;
            	while (split >> valueDia)
		{
        	        rDiametersMatrix.back().push_back(valueDia);
        	}
    	}

   	QLength InletRadius;
	double DimlessInletRadius;

	std::ifstream inLen("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/RealNetwork/AdjacencyLengths.txt");

    	std::vector<std::vector<double> > rLengthsMatrix;
    	string lineLen;

    	while (std::getline(inLen, lineLen)) 
   	{
            	rLengthsMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(lineLen);
           	double valueLen;
            	while (split >> valueLen)
		{
        	        rLengthsMatrix.back().push_back(valueLen);
        	}
    	}


    
    	// Generate the network

    	VesselNetworkGenerator<3> network_generator;
    	std::shared_ptr<VesselNetwork<3> > p_network;

   	auto p_file_handler = std::make_shared<OutputFileHandler>("RealNetwork_RT", true);
	double PerfusionQuotient;

    	auto p_viscosity_calculator = ViscosityCalculator<3>::Create();
    	p_viscosity_calculator->SetPlasmaViscosity(viscosity);
	auto p_impedance_calculator = VesselImpedanceCalculator<3>::Create();
	FlowSolver<3> flow_solver;
	//flow_solver.SetUseDirectSolver(true);


	p_network = network_generator.GenerateRealNetworkFromMatrices3D(rCoordinatesMatrix, rDiametersMatrix);

	std::vector<std::shared_ptr<Vessel<3> > > vessels = p_network->GetVessels();

	// set inlet and outlet nodes
	std::vector<std::shared_ptr<Vessel<3> > >::iterator vessel_iterator;

	for (vessel_iterator = vessels.begin();vessel_iterator != vessels.end(); vessel_iterator++)
	{
		(*vessel_iterator)->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);

	       	if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
	        {
	               	if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] == rCoordinatesMatrix[inletIndex][0]  && (*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[1] ==  rCoordinatesMatrix[inletIndex][1] && (*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[2] ==  rCoordinatesMatrix[inletIndex][2])
	                {
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
	                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
				InletRadius = (*vessel_iterator)->GetRadius();
				DimlessInletRadius = InletRadius/1_um;
	                }
			else
	                {
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
	                }
		}
	
	        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)   	
				//std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
		{
			if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] == rCoordinatesMatrix[inletIndex][0]  && (*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[1] ==  rCoordinatesMatrix[inletIndex][1] && (*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[2] ==  rCoordinatesMatrix[inletIndex][2])
			{
	  			(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
				InletRadius = (*vessel_iterator)->GetRadius();
				DimlessInletRadius = InletRadius/1_um;
	  	        }
			else
			{
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
	  	        }
	  	}
	}

	// Here we kill vessels from thinnest, up to the threshold
	unsigned MaxThreshold = 10;
	double DimlessRadiusThreshold;
	QLength RadiusThreshold;
	std::string output_file_real_RT;

	std::vector<double> QuotientVector(MaxThreshold+2);
	std::ofstream outfile;
    	outfile.open("RealNetworkPerfusion.txt");//, std::ios_base::app);

	//p_impedance_calculator->SetVesselNetwork(p_network);

    	//for(unsigned KillingLevel=0; KillingLevel < MaxThreshold+3; KillingLevel++)
    	for(unsigned KillingLevel=0; KillingLevel < MaxThreshold+2; KillingLevel++)	
	{
   		output_file_real_RT = p_file_handler->GetOutputDirectoryFullPath().append("RealNetworkRT"+to_string(KillingLevel)+".vtp");
		RadiusThreshold = double(KillingLevel)*InletRadius/double(MaxThreshold);
		DimlessRadiusThreshold = double(KillingLevel)*DimlessInletRadius/double(MaxThreshold);
		vessels = p_network->GetVessels();
		p_network->RemoveThinVessels(RadiusThreshold, false);
		p_viscosity_calculator->SetVesselNetwork(p_network);
    		p_impedance_calculator->SetVesselNetwork(p_network);
    		flow_solver.SetVesselNetwork(p_network);
		p_viscosity_calculator->Calculate();
		p_impedance_calculator->CalculateFromMatrixWithPruning(rLengthsMatrix,rDiametersMatrix,DimlessRadiusThreshold);
		flow_solver.SetUp();
		flow_solver.Solve();
		p_network->Write(output_file_real_RT);
	
		PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
        	QuotientVector[KillingLevel] = PerfusionQuotient;
	    	outfile << RadiusThreshold << " " << QuotientVector[KillingLevel] << " \n";
		//p_network->RemoveShortVessels(length_cut_off, false);	
	}

	outfile.close();

}

void xTestRealNetworkRTPickThickInlets()
{
	
	//unsigned inletIndex = 79;
	QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
   	double inlet_haematocrit = 0.45;
    	QFlowRate threshold = 1.e-15*unit::metre_cubed_per_second;
	std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/RealNetwork/33_2B/NodeCoordinates.txt");

    	std::vector<std::vector<double> > rCoordinatesMatrix;
    	string line;

    	while (std::getline(in, line)) 
   	{
            	rCoordinatesMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(line);
           	double value;
            	while (split >> value)
		{
        	        rCoordinatesMatrix.back().push_back(value);
        	}
    	}


	std::ifstream inDia("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/RealNetwork/33_2B/AdjacencyDiameters.txt");

    	std::vector<std::vector<double> > rDiametersMatrix;
    	string lineDia;

    	while (std::getline(inDia, lineDia)) 
   	{
            	rDiametersMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(lineDia);
           	double valueDia;
            	while (split >> valueDia)
		{
        	        rDiametersMatrix.back().push_back(valueDia);
        	}
    	}

   	//QLength InletRadius;
	//double DimlessInletRadius;
	//QLength RadiusThrForInlet = 40_um;
	QLength RadiusThrForInlet = 20_um;
	double DimlessRadiusThrForInlet = RadiusThrForInlet/1_um;

	std::ifstream inLen("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/RealNetwork/33_2B/AdjacencyLengths.txt");

    	std::vector<std::vector<double> > rLengthsMatrix;
    	string lineLen;

    	while (std::getline(inLen, lineLen)) 
   	{
            	rLengthsMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(lineLen);
           	double valueLen;
            	while (split >> valueLen)
		{
        	        rLengthsMatrix.back().push_back(valueLen);
        	}
    	}


    
    	// Generate the network

    	VesselNetworkGenerator<3> network_generator;
    	std::shared_ptr<VesselNetwork<3> > p_network;

   	auto p_file_handler = std::make_shared<OutputFileHandler>("RealNetwork_RT_ThickInlets_33_2B", true);
	double PerfusionQuotient;

    	auto p_viscosity_calculator = ViscosityCalculator<3>::Create();
    	p_viscosity_calculator->SetPlasmaViscosity(viscosity);
	auto p_impedance_calculator = VesselImpedanceCalculator<3>::Create();
	FlowSolver<3> flow_solver;
	flow_solver.SetUseDirectSolver(false);


	p_network = network_generator.GenerateRealNetworkFromMatrices3D(rCoordinatesMatrix, rDiametersMatrix);

	std::vector<std::shared_ptr<Vessel<3> > > vessels = p_network->GetVessels();

	// set inlet and outlet nodes
	std::vector<std::shared_ptr<Vessel<3> > >::iterator vessel_iterator;

	for (vessel_iterator = vessels.begin();vessel_iterator != vessels.end(); vessel_iterator++)
	{
		(*vessel_iterator)->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);

	       	if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
	        {
	               	if((*vessel_iterator)->GetRadius()>RadiusThrForInlet)
	                {
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
	                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
				//InletRadius = (*vessel_iterator)->GetRadius();
				//DimlessInletRadius = InletRadius/1_um;
	                }
			else
	                {
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
	                }
		}
	
	        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)   	
				//std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
		{
			if((*vessel_iterator)->GetRadius()>RadiusThrForInlet)
			{
	  			(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
				//InletRadius = (*vessel_iterator)->GetRadius();
				//DimlessInletRadius = InletRadius/1_um;
	  	        }
			else
			{
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
	  	        }
	  	}
	}

	// Here we kill vessels from thinnest, up to the threshold
	unsigned MaxThreshold = 10;
	double DimlessRadiusThreshold;
	QLength RadiusThreshold;
	std::string output_file_real_RT;

	std::vector<double> QuotientVector(MaxThreshold+1);
	std::ofstream outfile;
	std::ofstream outfileCorrelation;
    	outfile.open("RealNetworkPerfusionThickInlets33_2B_Indirect.txt");//, std::ios_base::app);
	vessels = p_network->GetVessels();
	//p_impedance_calculator->SetVesselNetwork(p_network);

    	//for(unsigned KillingLevel=0; KillingLevel < MaxThreshold+3; KillingLevel++)
    	//for(unsigned KillingLevel=0; KillingLevel < MaxThreshold+1; KillingLevel++)	
    	for(unsigned KillingLevel=0; KillingLevel < 5; KillingLevel++)		
	{
   		output_file_real_RT = p_file_handler->GetOutputDirectoryFullPath().append("RealNetwork33_2B_RTThickInlets"+to_string(KillingLevel)+".vtp");
		RadiusThreshold = double(KillingLevel)*RadiusThrForInlet/double(MaxThreshold);
		DimlessRadiusThreshold = double(KillingLevel)*DimlessRadiusThrForInlet/double(MaxThreshold);
		p_network->RemoveThinVessels(RadiusThreshold, false);
		vessels = p_network->GetVessels();
		p_viscosity_calculator->SetVesselNetwork(p_network);
    		p_impedance_calculator->SetVesselNetwork(p_network);
    		flow_solver.SetVesselNetwork(p_network);
		p_viscosity_calculator->Calculate();
		p_impedance_calculator->CalculateFromMatrixWithPruning(rLengthsMatrix,rDiametersMatrix,DimlessRadiusThreshold);
		flow_solver.SetUp();
		flow_solver.Solve();
		p_network->Write(output_file_real_RT);

	
		PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
        	QuotientVector[KillingLevel] = PerfusionQuotient;
	    	outfile << RadiusThreshold << " " << QuotientVector[KillingLevel] << " \n";
    		outfileCorrelation.open("RealNetwork33_2B_PerfusionThickInletsCorrelations"+to_string(KillingLevel)+"ManyInlets.txt");
		for (std::vector<std::shared_ptr<Vessel<3> > >::iterator itt = vessels.begin();itt != vessels.end(); itt++)
		{

			outfileCorrelation << (*itt)->GetRadius() << " " <<  fabs((*itt)->GetFlowProperties()->GetFlowRate()) << " \n";
			//outfileCorrelation << fabs((*itt)->GetFlowProperties()->GetFlowRate()) << " \n";
			//(*vessel_iterator)->GetRadius();
			//fabs((*vessel_iterator)->GetFlowRate()) / unit::metre_cubed_per_second
		}
		outfileCorrelation.close();
		//p_network->RemoveShortVessels(length_cut_off, false);	
	}

	outfile.close();

}

void xTestFullRealNetworkRTPickThickInlets()
{
	
	//unsigned inletIndex = 79;
	QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
   	double inlet_haematocrit = 0.45;
    	QFlowRate threshold = 3.5e-12*unit::metre_cubed_per_second;
	QLength RadiusThrForInlet = 10_um;
	QLength RadiusThresholdFromBetti = 10_um;

	unsigned NumberOfKilling = 11;
		
	//string CellLine = "24_2B";
	//string CellLine = "18_4C";
	string CellLine = "35_2C";
        //string CellLine = "33_2B";
	std::ifstream in("/home/narain/RT-enhanced perfusionShared/Networks/RealNetworks/IR/"+CellLine+"/Full Network/FullNodeCoordinates"+CellLine+".txt");

    	std::vector<std::vector<double> > rCoordinatesMatrix;
    	string line;

    	while (std::getline(in, line)) 
   	{
            	rCoordinatesMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(line);
           	double value;
            	while (split >> value)
		{
        	        rCoordinatesMatrix.back().push_back(value);
        	}
    	}


	std::ifstream inDia("/home/narain/RT-enhanced perfusionShared/Networks/NetworksRealNetworks/IR/"+CellLine+"/Full Network/FullAdjacencyDiameters"+CellLine+".txt");

    	std::vector<std::vector<double> > rDiametersMatrix;
    	string lineDia;

    	while (std::getline(inDia, lineDia)) 
   	{
            	rDiametersMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(lineDia);
           	double valueDia;
            	while (split >> valueDia)
		{
        	        rDiametersMatrix.back().push_back(valueDia);
        	}
    	}

   	//QLength InletRadius;
	//double DimlessInletRadius;
	//QLength RadiusThrForInlet = 40_um;

	//QLength RadiusThrForInlet = 25_um; this is for cell line 33_2B
	//double DimlessRadiusThrForInlet = RadiusThrForInlet/1_um;

	std::ifstream inLen("/home/narain/RT-enhanced perfusionShared/Networks/RealNetworks/IR/"+CellLine+"/Full Network/FullAdjacencyLengths"+CellLine+".txt");

    	std::vector<std::vector<double> > rLengthsMatrix;
    	string lineLen;

    	while (std::getline(inLen, lineLen)) 
   	{
            	rLengthsMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(lineLen);
           	double valueLen;
            	while (split >> valueLen)
		{
        	        rLengthsMatrix.back().push_back(valueLen);
        	}
    	}


    
    	// Generate the network

    	VesselNetworkGenerator<3> network_generator;
    	std::shared_ptr<VesselNetwork<3> > p_network;

   	auto p_file_handler = std::make_shared<OutputFileHandler>("FullRealNetwork_RT_ThickInlets_"+CellLine, true);
	double PerfusionQuotient;

    	auto p_viscosity_calculator = ViscosityCalculator<3>::Create();
    	p_viscosity_calculator->SetPlasmaViscosity(viscosity);
	auto p_impedance_calculator = VesselImpedanceCalculator<3>::Create();
	FlowSolver<3> flow_solver;
	flow_solver.SetUseDirectSolver(true);


	p_network = network_generator.GenerateRealNetworkFromMatrices3D(rCoordinatesMatrix, rDiametersMatrix);

	std::vector<std::shared_ptr<Vessel<3> > > vessels = p_network->GetVessels();

	// set inlet and outlet nodes
	std::vector<std::shared_ptr<Vessel<3> > >::iterator vessel_iterator;

	for (vessel_iterator = vessels.begin();vessel_iterator != vessels.end(); vessel_iterator++)
	{
		(*vessel_iterator)->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);

	       	if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
	        {
	               	if((*vessel_iterator)->GetRadius()>RadiusThrForInlet)
	                {
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
	                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
				//InletRadius = (*vessel_iterator)->GetRadius();
				//DimlessInletRadius = InletRadius/1_um;
	                }
			else
	                {
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
	                }
		}
	
	        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)   	
				//std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
		{
			if((*vessel_iterator)->GetRadius()>RadiusThrForInlet)
			{
	  			(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
				//InletRadius = (*vessel_iterator)->GetRadius();
				//DimlessInletRadius = InletRadius/1_um;
	  	        }
			else
			{
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
	  	        }
	  	}
	}

	// Here we kill vessels from thinnest, up to the threshold
	//unsigned MaxThreshold = 10;

	double DimlessRadiusThreshold;
	QLength RadiusThreshold;
	std::string output_file_real_RT;

	std::vector<double> QuotientVector(NumberOfKilling);
	std::ofstream outfile;
	std::ofstream outfileCorrelation;
    	outfile.open("FullRealNetworkPerfusionThickInlets"+CellLine+"_Direct.txt");//, std::ios_base::app);
	vessels = p_network->GetVessels();
	//p_impedance_calculator->SetVesselNetwork(p_network);

    	//for(unsigned KillingLevel=0; KillingLevel < MaxThreshold+3; KillingLevel++)
    	//for(unsigned KillingLevel=0; KillingLevel < MaxThreshold+1; KillingLevel++)	
    	for(unsigned KillingLevel=0; KillingLevel < NumberOfKilling; KillingLevel++)		
	{
   		output_file_real_RT = p_file_handler->GetOutputDirectoryFullPath().append("FullRealNetwork"+CellLine+"_RTThickInlets"+to_string(KillingLevel)+".vtp");
		RadiusThreshold = double(KillingLevel)*RadiusThresholdFromBetti/double(NumberOfKilling-1);
		DimlessRadiusThreshold = RadiusThreshold/1_um;
		p_network->RemoveThinVessels(RadiusThreshold, false);
		vessels = p_network->GetVessels();
		p_viscosity_calculator->SetVesselNetwork(p_network);
    		p_impedance_calculator->SetVesselNetwork(p_network);
    		flow_solver.SetVesselNetwork(p_network);
		p_viscosity_calculator->Calculate();
		p_impedance_calculator->CalculateFromMatrixWithPruning(rLengthsMatrix,rDiametersMatrix,DimlessRadiusThreshold);
		flow_solver.SetUp();
		flow_solver.Solve();
		p_network->Write(output_file_real_RT);

	
		PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
        	QuotientVector[KillingLevel] = PerfusionQuotient;
	    	outfile << RadiusThreshold << " " << QuotientVector[KillingLevel] << " \n";
    		outfileCorrelation.open("FullRealNetwork"+CellLine+"_PerfusionThickInletsCorrelations"+to_string(KillingLevel)+"ManyInlets.txt");
		for (std::vector<std::shared_ptr<Vessel<3> > >::iterator itt = vessels.begin();itt != vessels.end(); itt++)
		{

			outfileCorrelation << (*itt)->GetRadius() << " " <<  fabs((*itt)->GetFlowProperties()->GetFlowRate()) << " \n";
			//outfileCorrelation << fabs((*itt)->GetFlowProperties()->GetFlowRate()) << " \n";
			//(*vessel_iterator)->GetRadius();
			//fabs((*vessel_iterator)->GetFlowRate()) / unit::metre_cubed_per_second
		}
		outfileCorrelation.close();
		//p_network->RemoveShortVessels(length_cut_off, false);	
	}

	outfile.close();

}

void xTestFullRealNetworkRTPickThickInletsConstFlowRate()
{
	
	//unsigned inletIndex = 79;
	QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
   	double inlet_haematocrit = 0.45;
    	QFlowRate threshold = 3.5e-12*unit::metre_cubed_per_second;
	QLength RadiusThrForInlet = 25_um;
	QLength RadiusThresholdFromBetti = 15_um;

	unsigned NumberOfKilling = 2;
		
	//string CellLine = "24_2B";
	//string CellLine = "18_4C";
	string CellLine = "35_2C";
        //string CellLine = "33_2B";
	std::ifstream in("/home/narain/RT-enhanced perfusionShared/Networks/RealNetworks/IR/"+CellLine+"/Full Network/FullNodeCoordinates"+CellLine+".txt");

    	std::vector<std::vector<double> > rCoordinatesMatrix;
    	string line;

    	while (std::getline(in, line)) 
   	{
            	rCoordinatesMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(line);
           	double value;
            	while (split >> value)
		{
        	        rCoordinatesMatrix.back().push_back(value);
        	}
    	}


	std::ifstream inDia("/home/narain/RT-enhanced perfusionShared/Networks/RealNetworks/IR/"+CellLine+"/Full Network/FullAdjacencyDiameters"+CellLine+".txt");

    	std::vector<std::vector<double> > rDiametersMatrix;
    	string lineDia;

    	while (std::getline(inDia, lineDia)) 
   	{
            	rDiametersMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(lineDia);
           	double valueDia;
            	while (split >> valueDia)
		{
        	        rDiametersMatrix.back().push_back(valueDia);
        	}
    	}

   	//QLength InletRadius;
	//double DimlessInletRadius;
	//QLength RadiusThrForInlet = 40_um;

	//QLength RadiusThrForInlet = 25_um; this is for cell line 33_2B
	//double DimlessRadiusThrForInlet = RadiusThrForInlet/1_um;

	std::ifstream inLen("/home/narain/RT-enhanced perfusionShared/Networks/RealNetworks/IR/"+CellLine+"/Full Network/FullAdjacencyLengths"+CellLine+".txt");

    	std::vector<std::vector<double> > rLengthsMatrix;
    	string lineLen;

    	while (std::getline(inLen, lineLen)) 
   	{
            	rLengthsMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(lineLen);
           	double valueLen;
            	while (split >> valueLen)
		{
        	        rLengthsMatrix.back().push_back(valueLen);
        	}
    	}


    
    	// Generate the network

    	VesselNetworkGenerator<3> network_generator;
    	std::shared_ptr<VesselNetwork<3> > p_network;

   	auto p_file_handler = std::make_shared<OutputFileHandler>("ConstFlowRateFullRealNetwork_RT_ThickInlets_"+CellLine, true);
	double PerfusionQuotient;

    	auto p_viscosity_calculator = ViscosityCalculator<3>::Create();
    	p_viscosity_calculator->SetPlasmaViscosity(viscosity);
	auto p_impedance_calculator = VesselImpedanceCalculator<3>::Create();
	FlowSolver<3> flow_solver;
	flow_solver.SetUseDirectSolver(true);


	p_network = network_generator.GenerateRealNetworkFromMatrices3D(rCoordinatesMatrix, rDiametersMatrix);

	std::vector<std::shared_ptr<Vessel<3> > > vessels = p_network->GetVessels();

	// set inlet and outlet nodes
	std::vector<std::shared_ptr<Vessel<3> > >::iterator vessel_iterator;

	for (vessel_iterator = vessels.begin();vessel_iterator != vessels.end(); vessel_iterator++)
	{
		(*vessel_iterator)->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);

	       	if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
	        {
	               	if((*vessel_iterator)->GetRadius()>RadiusThrForInlet)
	                {
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetUseVelocityBoundaryCondition(true);
	                	(*vessel_iterator)->GetStartNode()->GetSegment(0)->GetFlowProperties()->SetFlowRate(1.e-10*unit::metre_cubed_per_second);
	                        //(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
				//InletRadius = (*vessel_iterator)->GetRadius();
				//DimlessInletRadius = InletRadius/1_um;
	                }
			else
	                {
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
	                }
		}
	
	        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)   	
				//std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
		{
			if((*vessel_iterator)->GetRadius()>RadiusThrForInlet)
			{
	  			(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
	                	(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetUseVelocityBoundaryCondition(true);
	                	(*vessel_iterator)->GetEndNode()->GetSegment(0)->GetFlowProperties()->SetFlowRate(1.e-10*unit::metre_cubed_per_second);
	  		        //(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
				//InletRadius = (*vessel_iterator)->GetRadius();
				//DimlessInletRadius = InletRadius/1_um;
	  	        }
			else
			{
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
	  	        }
	  	}
	}

	// Here we kill vessels from thinnest, up to the threshold
	//unsigned MaxThreshold = 10;

	double DimlessRadiusThreshold;
	QLength RadiusThreshold;
	std::string output_file_real_RT;
	std::string output_file_real_Initial;
   	output_file_real_Initial = p_file_handler->GetOutputDirectoryFullPath().append("FullRealNetwork"+CellLine+"_RTThickInletsInitial.vtp");
	p_network->Write(output_file_real_Initial);

	std::vector<double> QuotientVector(NumberOfKilling);
	std::ofstream outfile;
	std::ofstream outfileCorrelation;
    	outfile.open("ConstFlowRateFullRealNetworkPerfusionThickInlets"+CellLine+"_Direct.txt");//, std::ios_base::app);
	vessels = p_network->GetVessels();
	//p_impedance_calculator->SetVesselNetwork(p_network);

    	//for(unsigned KillingLevel=0; KillingLevel < MaxThreshold+3; KillingLevel++)
    	//for(unsigned KillingLevel=0; KillingLevel < MaxThreshold+1; KillingLevel++)	
    	for(unsigned KillingLevel=0; KillingLevel < NumberOfKilling-1; KillingLevel++)		
	{
   		output_file_real_RT = p_file_handler->GetOutputDirectoryFullPath().append("ConstFlowRateFullRealNetwork"+CellLine+"_RTThickInlets"+to_string(KillingLevel)+".vtp");
		RadiusThreshold = double(KillingLevel)*RadiusThresholdFromBetti/double(NumberOfKilling-1);
		DimlessRadiusThreshold = RadiusThreshold/1_um;
		p_network->RemoveThinVessels(RadiusThreshold, false);
		vessels = p_network->GetVessels();
		p_viscosity_calculator->SetVesselNetwork(p_network);
    		p_impedance_calculator->SetVesselNetwork(p_network);
    		flow_solver.SetVesselNetwork(p_network);
		p_viscosity_calculator->Calculate();
		p_impedance_calculator->CalculateFromMatrixWithPruning(rLengthsMatrix,rDiametersMatrix,DimlessRadiusThreshold);
		flow_solver.SetUp();
		flow_solver.Solve();
		p_network->Write(output_file_real_RT);

	
		PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
        	QuotientVector[KillingLevel] = PerfusionQuotient;
	    	outfile << RadiusThreshold << " " << QuotientVector[KillingLevel] << " \n";
    		outfileCorrelation.open("ConstFlowRateFullRealNetwork"+CellLine+"_PerfusionThickInletsCorrelations"+to_string(KillingLevel)+"ManyInlets.txt");
		for (std::vector<std::shared_ptr<Vessel<3> > >::iterator itt = vessels.begin();itt != vessels.end(); itt++)
		{

			outfileCorrelation << (*itt)->GetRadius() << " " <<  fabs((*itt)->GetFlowProperties()->GetFlowRate()) << " \n";
			//outfileCorrelation << fabs((*itt)->GetFlowProperties()->GetFlowRate()) << " \n";
			//(*vessel_iterator)->GetRadius();
			//fabs((*vessel_iterator)->GetFlowRate()) / unit::metre_cubed_per_second
		}
		outfileCorrelation.close();
		//p_network->RemoveShortVessels(length_cut_off, false);	
	}

	outfile.close();

}

void xTestLargestComponentRealNetworkRTPickThickInlets()
{
	
	//unsigned inletIndex = 79;
	QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
   	double inlet_haematocrit = 0.45;
    	QFlowRate threshold = 1.e-14*unit::metre_cubed_per_second;
	//string CellLine = "24_2B";
	//string CellLine = "18_4C";
	//string CellLine = "35_2C";
        string CellLine = "33_2B";
	std::ifstream in("/home/narain/RT-enhanced perfusionShared/Networks/RealNetworks/IR/"+CellLine+"/Largest Component/NodeCoordinates"+CellLine+".txt");

    	std::vector<std::vector<double> > rCoordinatesMatrix;
    	string line;

    	while (std::getline(in, line)) 
   	{
            	rCoordinatesMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(line);
           	double value;
            	while (split >> value)
		{
        	        rCoordinatesMatrix.back().push_back(value);
        	}
    	}


	std::ifstream inDia("/home/narain/RT-enhanced perfusionShared/Networks/RealNetworks/IR/"+CellLine+"/Largest Component/AdjacencyDiameters"+CellLine+".txt");

    	std::vector<std::vector<double> > rDiametersMatrix;
    	string lineDia;

    	while (std::getline(inDia, lineDia)) 
   	{
            	rDiametersMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(lineDia);
           	double valueDia;
            	while (split >> valueDia)
		{
        	        rDiametersMatrix.back().push_back(valueDia);
        	}
    	}

   	//QLength InletRadius;
	//double DimlessInletRadius;
	//QLength RadiusThrForInlet = 40_um;
	QLength RadiusThrForInlet = 20_um;		
	//QLength RadiusThrForInlet = 25_um; this is for cell line 33_2B
	//double DimlessRadiusThrForInlet = RadiusThrForInlet/1_um;

	std::ifstream inLen("/home/narain/RT-enhanced perfusionShared/Networks/RealNetworks/IR/"+CellLine+"/Largest Component/AdjacencyLengths"+CellLine+".txt");

    	std::vector<std::vector<double> > rLengthsMatrix;
    	string lineLen;

    	while (std::getline(inLen, lineLen)) 
   	{
            	rLengthsMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(lineLen);
           	double valueLen;
            	while (split >> valueLen)
		{
        	        rLengthsMatrix.back().push_back(valueLen);
        	}
    	}


    
    	// Generate the network

    	VesselNetworkGenerator<3> network_generator;
    	std::shared_ptr<VesselNetwork<3> > p_network;

   	auto p_file_handler = std::make_shared<OutputFileHandler>("LargestComponentRealNetwork_RT_ThickInlets_"+CellLine, true);
	double PerfusionQuotient;

    	auto p_viscosity_calculator = ViscosityCalculator<3>::Create();
    	p_viscosity_calculator->SetPlasmaViscosity(viscosity);
	auto p_impedance_calculator = VesselImpedanceCalculator<3>::Create();
	FlowSolver<3> flow_solver;
	flow_solver.SetUseDirectSolver(true);


	p_network = network_generator.GenerateRealNetworkFromMatrices3D(rCoordinatesMatrix, rDiametersMatrix);

	std::vector<std::shared_ptr<Vessel<3> > > vessels = p_network->GetVessels();

	// set inlet and outlet nodes
	std::vector<std::shared_ptr<Vessel<3> > >::iterator vessel_iterator;

	for (vessel_iterator = vessels.begin();vessel_iterator != vessels.end(); vessel_iterator++)
	{
		(*vessel_iterator)->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);

	       	if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
	        {
	               	if((*vessel_iterator)->GetRadius()>RadiusThrForInlet)
	                {
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
	                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
				//InletRadius = (*vessel_iterator)->GetRadius();
				//DimlessInletRadius = InletRadius/1_um;
	                }
			else
	                {
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
	                }
		}
	
	        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)   	
				//std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
		{
			if((*vessel_iterator)->GetRadius()>RadiusThrForInlet)
			{
	  			(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
				//InletRadius = (*vessel_iterator)->GetRadius();
				//DimlessInletRadius = InletRadius/1_um;
	  	        }
			else
			{
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
	  	        }
	  	}
	}

	// Here we kill vessels from thinnest, up to the threshold
	unsigned NumberOfKilling = 4;
	double DimlessRadiusThreshold;
	QLength RadiusThreshold;
	QLength RadiusThresholdFromBetti = 15_um;
	std::string output_file_real_RT;

	std::vector<double> QuotientVector(NumberOfKilling);
	std::ofstream outfile;
	std::ofstream outfileCorrelation;
    	outfile.open("LargestComponentRealNetworkPerfusionThickInlets"+CellLine+"_Direct.txt");//, std::ios_base::app);
	vessels = p_network->GetVessels();
	//p_impedance_calculator->SetVesselNetwork(p_network);

    	//for(unsigned KillingLevel=0; KillingLevel < NumberOfKilling+3; KillingLevel++)
    	//for(unsigned KillingLevel=0; KillingLevel < NumberOfKilling+1; KillingLevel++)	
    	for(unsigned KillingLevel=0; KillingLevel < NumberOfKilling; KillingLevel++)		
	{
   		output_file_real_RT = p_file_handler->GetOutputDirectoryFullPath().append("LargestComponentRealNetwork"+CellLine+"_RTThickInlets"+to_string(KillingLevel)+".vtp");
		//RadiusThreshold = 1.5*double(KillingLevel)*RadiusThrForInlet/double(NumberOfKilling);
		//DimlessRadiusThreshold = 1.5*double(KillingLevel)*DimlessRadiusThrForInlet/double(NumberOfKilling);
		RadiusThreshold = double(KillingLevel)*RadiusThresholdFromBetti/double(NumberOfKilling-1);
		DimlessRadiusThreshold = RadiusThreshold/1_um;
		p_network->RemoveThinVessels(RadiusThreshold, false);
		vessels = p_network->GetVessels();
		p_viscosity_calculator->SetVesselNetwork(p_network);
    		p_impedance_calculator->SetVesselNetwork(p_network);
    		flow_solver.SetVesselNetwork(p_network);
		p_viscosity_calculator->Calculate();
		p_impedance_calculator->CalculateFromMatrixWithPruning(rLengthsMatrix,rDiametersMatrix,DimlessRadiusThreshold);
		flow_solver.SetUp();
		flow_solver.Solve();
		p_network->Write(output_file_real_RT);

	
		PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
        	QuotientVector[KillingLevel] = PerfusionQuotient;
	    	outfile << RadiusThreshold << " " << QuotientVector[KillingLevel] << " \n";
    		outfileCorrelation.open("LargestComponentRealNetwork"+CellLine+"_PerfusionThickInletsCorrelations"+to_string(KillingLevel)+"ManyInlets.txt");
		for (std::vector<std::shared_ptr<Vessel<3> > >::iterator itt = vessels.begin();itt != vessels.end(); itt++)
		{

			outfileCorrelation << (*itt)->GetRadius() << " " <<  fabs((*itt)->GetFlowProperties()->GetFlowRate()) << " \n";
			//outfileCorrelation << fabs((*itt)->GetFlowProperties()->GetFlowRate()) << " \n";
			//(*vessel_iterator)->GetRadius();
			//fabs((*vessel_iterator)->GetFlowRate()) / unit::metre_cubed_per_second
		}
		outfileCorrelation.close();
		//p_network->RemoveShortVessels(length_cut_off, false);	
	}

	outfile.close();

}

void xTestRealNetworkRTKillByLength()
{
	
	
	unsigned inletIndex = 79;
	QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
   	double inlet_haematocrit = 0.45;
    	QFlowRate threshold = 1.e-15*unit::metre_cubed_per_second;
	std::ifstream in("/home/narain/RT-enhanced perfusionShared/Networks/RealNetworks/700Edges/NodeCoordinates.txt");

    	std::vector<std::vector<double> > rCoordinatesMatrix;
    	string line;

    	while (std::getline(in, line)) 
   	{
            	rCoordinatesMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(line);
           	double value;
            	while (split >> value)
		{
        	        rCoordinatesMatrix.back().push_back(value);
        	}
    	}


	std::ifstream inDia("/home/narain/RT-enhanced perfusionShared/Networks/RealNetworks/AdjacencyDiameters.txt");

    	std::vector<std::vector<double> > rDiametersMatrix;
    	string lineDia;

    	while (std::getline(inDia, lineDia)) 
   	{
            	rDiametersMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(lineDia);
           	double valueDia;
            	while (split >> valueDia)
		{
        	        rDiametersMatrix.back().push_back(valueDia);
        	}
    	}

   	QLength InletLength;
	double DimlessInletLength;

	std::ifstream inLen("/home/narain/RT-enhanced perfusionShared/Networks/RealNetworks/AdjacencyLengths.txt");

    	std::vector<std::vector<double> > rLengthsMatrix;
    	string lineLen;

	unsigned SpIndex = 0;

    	while (std::getline(inLen, lineLen)) 
   	{
            	rLengthsMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(lineLen);
           	double valueLen;
            	while (split >> valueLen)
		{
        	        rLengthsMatrix.back().push_back(valueLen);
			if (valueLen!=0 && SpIndex == inletIndex)
			{
				DimlessInletLength = valueLen;
				InletLength = valueLen*1_um;
			}
        	}
		SpIndex++;
    	}

	std::cout << DimlessInletLength << " and dimensional... " << InletLength << "\n";

    
    	// Generate the network

    	VesselNetworkGenerator<3> network_generator;
    	std::shared_ptr<VesselNetwork<3> > p_network;

   	auto p_file_handler = std::make_shared<OutputFileHandler>("RealNetworkByLength_RT", true);
	double PerfusionQuotient;

    	auto p_viscosity_calculator = ViscosityCalculator<3>::Create();
    	p_viscosity_calculator->SetPlasmaViscosity(viscosity);
	auto p_impedance_calculator = VesselImpedanceCalculator<3>::Create();
	FlowSolver<3> flow_solver;
	//flow_solver.SetUseDirectSolver(true);


	p_network = network_generator.GenerateRealNetworkFromMatrices3DWithLengths(rCoordinatesMatrix, rDiametersMatrix, rLengthsMatrix);

	std::vector<std::shared_ptr<Vessel<3> > > vessels = p_network->GetVessels();

	// set inlet and outlet nodes
	std::vector<std::shared_ptr<Vessel<3> > >::iterator vessel_iterator;

	for (vessel_iterator = vessels.begin();vessel_iterator != vessels.end(); vessel_iterator++)
	{
		(*vessel_iterator)->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);

	       	if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
	        {
	               	if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] == rCoordinatesMatrix[inletIndex][0]  && (*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[1] ==  rCoordinatesMatrix[inletIndex][1] && (*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[2] ==  rCoordinatesMatrix[inletIndex][2])
	                {
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
	                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
				//InletRadius = (*vessel_iterator)->GetRadius();
				//DimlessInletRadius = InletRadius/1_um;
	                }
			else
	                {
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
	                }
		}
	
	        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)   	
				//std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
		{
			if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] == rCoordinatesMatrix[inletIndex][0]  && (*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[1] ==  rCoordinatesMatrix[inletIndex][1] && (*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[2] ==  rCoordinatesMatrix[inletIndex][2])
			{
	  			(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
				//InletRadius = (*vessel_iterator)->GetRadius();
				//DimlessInletRadius = InletRadius/1_um;
	  	        }
			else
			{
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
	  	        }
	  	}
	}

	// Here we kill vessels from thinnest, up to the threshold
	unsigned MaxThreshold = 10;
	double DimlessLengthThreshold;
	//double PrevDimlessLengthThreshold;
	QLength LengthThreshold;
	std::string output_file_real_RT;

	std::vector<double> QuotientVector(MaxThreshold+5);
	std::ofstream outfile;
    	outfile.open("RealNetworkPerfusionByLength.txt");//, std::ios_base::app);

    	std::vector<std::shared_ptr<Vessel<3> > > vessels_to_remove;

	//p_impedance_calculator->SetVesselNetwork(p_network);

    	//for(unsigned KillingLevel=0; KillingLevel < MaxThreshold+3; KillingLevel++)
    	for(unsigned KillingLevel=0; KillingLevel < MaxThreshold+5; KillingLevel++)	
	{
	//	unsigned KillingLevel=0;
   		output_file_real_RT = p_file_handler->GetOutputDirectoryFullPath().append("RealNetworkRTByLength"+to_string(KillingLevel)+".vtp");
		LengthThreshold = 0.2*double(KillingLevel)*InletLength/double(MaxThreshold);
		DimlessLengthThreshold = 0.2*double(KillingLevel)*DimlessInletLength/double(MaxThreshold);

		//std::cout << DimlessLengthThreshold << " and dimensional... " << LengthThreshold << "\n";
		std::cout << "Dimensionless length threshold is: " << DimlessLengthThreshold << " and dimensional... " << LengthThreshold <<  "\n";

		p_network->RemoveShortVesselsLengthFromMatrix(LengthThreshold, false);

		vessels = p_network->GetVessels();
		p_viscosity_calculator->SetVesselNetwork(p_network);
    		p_impedance_calculator->SetVesselNetwork(p_network);
    		flow_solver.SetVesselNetwork(p_network);
		p_viscosity_calculator->Calculate();
		p_impedance_calculator->CalculateGetLengthFromMatrix();//FromMatrixWithPruningByLength(rLengthsMatrix,DimlessLengthThreshold);
		flow_solver.SetUp();
		flow_solver.Solve();
		p_network->Write(output_file_real_RT);
	
		PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
        	QuotientVector[KillingLevel] = PerfusionQuotient;
	    	outfile << LengthThreshold << " " << QuotientVector[KillingLevel] << " \n";
		//p_network->RemoveShortVessels(length_cut_off, false);	
	}

	outfile.close();

}

/* I might use this to determine the order of killing
  	//string line2;
	//std::vector<std::vector<unsigned> > Order;
	//VesselPtr<2> p_thin_vessel;
	std::ifstream in2("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/Selection"+to_string(att+1)+".txt");
	Order.clear();
	while (std::getline(in2, line2)) 
   	{
            	Order.push_back(std::vector<unsigned>());
            		// Break down the row into column values
         	std::stringstream split2(line2);
           	unsigned value2;
            	while (split2 >> value2)
		{
        	        Order.back().push_back(value2);
        	}
    	}*/

	/* Here we write into output files
    
    		for (unsigned iii = 0; iii < NHet; iii++)
    		{
        		for(unsigned jjj = 0; jjj < ToBeKilled+1; jjj++)
        		{
	    		QuotientMatrixMean[iii][jjj]= QuotientMatrixAggregate[iii][jjj]/double(att+1);
	   	 	std::cout << "For heterogeneity level " << iii << ", and radiotherapy dose " << jjj << ", the aggregate quotient is: " << QuotientMatrixAggregate[iii][jjj] << " , and the mean quotient is: " << QuotientMatrixMean[iii][jjj] << endl;

        		}
    		}    



	*/

void xTestRealNetworkRTKillByLengthPickLongInlets()
{
	
	
	//unsigned inletIndex = 79;
	QLength LengthThrForInlet = 380_um;
	//double DimlessLengthThrForInlet = LengthThrForInlet/1_um;
	QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
   	double inlet_haematocrit = 0.45;
    	QFlowRate threshold = 1.e-15*unit::metre_cubed_per_second;
	std::ifstream in("/home/narain/RT-enhanced perfusionShared/Networks/RealNetworks/NodeCoordinates.txt");

    	std::vector<std::vector<double> > rCoordinatesMatrix;
    	string line;

    	while (std::getline(in, line)) 
   	{
            	rCoordinatesMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(line);
           	double value;
            	while (split >> value)
		{
        	        rCoordinatesMatrix.back().push_back(value);
        	}
    	}


	std::ifstream inDia("/home/narain/RT-enhanced perfusionShared/Networks/RealNetworks/AdjacencyDiameters.txt");

    	std::vector<std::vector<double> > rDiametersMatrix;
    	string lineDia;

    	while (std::getline(inDia, lineDia)) 
   	{
            	rDiametersMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(lineDia);
           	double valueDia;
            	while (split >> valueDia)
		{
        	        rDiametersMatrix.back().push_back(valueDia);
        	}
    	}

   	//QLength InletLength;
	//double DimlessInletLength;

	std::ifstream inLen("/home/narain/RT-enhanced perfusionShared/Networks/RealNetworks/AdjacencyLengths.txt");

    	std::vector<std::vector<double> > rLengthsMatrix;
    	string lineLen;

	//unsigned SpIndex = 0;

    	while (std::getline(inLen, lineLen)) 
   	{
            	rLengthsMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(lineLen);
           	double valueLen;
            	while (split >> valueLen)
		{
        	        rLengthsMatrix.back().push_back(valueLen);
			//if (valueLen!=0 && SpIndex == inletIndex)
			//{
			//	DimlessInletLength = valueLen;
			//	InletLength = valueLen*1_um;
			//}
        	}
		//SpIndex++;
    	}

	//std::cout << DimlessInletLength << " and dimensional... " << InletLength << "\n";

    
    	// Generate the network

    	VesselNetworkGenerator<3> network_generator;
    	std::shared_ptr<VesselNetwork<3> > p_network;

   	auto p_file_handler = std::make_shared<OutputFileHandler>("RealNetworkByLength_RT_PickLongInlets", true);
	double PerfusionQuotient;

    	auto p_viscosity_calculator = ViscosityCalculator<3>::Create();
    	p_viscosity_calculator->SetPlasmaViscosity(viscosity);
	auto p_impedance_calculator = VesselImpedanceCalculator<3>::Create();
	FlowSolver<3> flow_solver;
	//flow_solver.SetUseDirectSolver(true);


	p_network = network_generator.GenerateRealNetworkFromMatrices3DWithLengths(rCoordinatesMatrix, rDiametersMatrix, rLengthsMatrix);

	std::vector<std::shared_ptr<Vessel<3> > > vessels = p_network->GetVessels();

	// set inlet and outlet nodes
	std::vector<std::shared_ptr<Vessel<3> > >::iterator vessel_iterator;

	for (vessel_iterator = vessels.begin();vessel_iterator != vessels.end(); vessel_iterator++)
	{
		(*vessel_iterator)->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);

	       	if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
	        {
	               	if((*vessel_iterator)->GetLengthFromMatrix()>LengthThrForInlet)
	                {
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
	                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
				//InletRadius = (*vessel_iterator)->GetRadius();
				//DimlessInletRadius = InletRadius/1_um;
	                }
			else
	                {
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
	                	(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
	                }
		}
	
	        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)   	
				//std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
		{
			if((*vessel_iterator)->GetLengthFromMatrix()>LengthThrForInlet)
			{
	  			(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
				//InletRadius = (*vessel_iterator)->GetRadius();
				//DimlessInletRadius = InletRadius/1_um;
	  	        }
			else
			{
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
	  		        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
	  	        }
	  	}
	}

	// Here we kill vessels from thinnest, up to the threshold
	unsigned MaxThreshold = 30;
	//double DimlessLengthThreshold;
	//double PrevDimlessLengthThreshold;
	QLength LengthThreshold;
	std::string output_file_real_RT;

	std::vector<double> QuotientVector(MaxThreshold+2);
	std::ofstream outfile;
    	outfile.open("RealNetworkPerfusionByLengthPickLongInlets.txt");//, std::ios_base::app);

    	std::vector<std::shared_ptr<Vessel<3> > > vessels_to_remove;

	//p_impedance_calculator->SetVesselNetwork(p_network);

    	//for(unsigned KillingLevel=0; KillingLevel < MaxThreshold+3; KillingLevel++)
    	for(unsigned KillingLevel=0; KillingLevel < MaxThreshold+2; KillingLevel++)	
	{
	//	unsigned KillingLevel=0;
   		output_file_real_RT = p_file_handler->GetOutputDirectoryFullPath().append("RealNetworkRTByLengthPickLongInlets"+to_string(KillingLevel)+".vtp");
		LengthThreshold = double(KillingLevel)*LengthThrForInlet/double(MaxThreshold); //*0.2
		//DimlessLengthThreshold = double(KillingLevel)*DimlessLengthThrForInlet/double(MaxThreshold);

		//std::cout << DimlessLengthThreshold << " and dimensional... " << LengthThreshold << "\n";
		//std::cout << "Dimensionless length threshold is: " << DimlessLengthThreshold << " and dimensional... " << LengthThreshold <<  "\n";

		p_network->RemoveShortVesselsLengthFromMatrix(LengthThreshold, false);

		vessels = p_network->GetVessels();
		p_viscosity_calculator->SetVesselNetwork(p_network);
    		p_impedance_calculator->SetVesselNetwork(p_network);
    		flow_solver.SetVesselNetwork(p_network);
		p_viscosity_calculator->Calculate();
		p_impedance_calculator->CalculateGetLengthFromMatrix();//FromMatrixWithPruningByLength(rLengthsMatrix,DimlessLengthThreshold);
		flow_solver.SetUp();
		flow_solver.Solve();
		p_network->Write(output_file_real_RT);
	
		PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
        	QuotientVector[KillingLevel] = PerfusionQuotient;
	    	outfile << LengthThreshold << " " << QuotientVector[KillingLevel] << " \n";
		//p_network->RemoveShortVessels(length_cut_off, false);	
	}

	outfile.close();

}

void xTestIrregularHexagonsFromVoronoi()
{
	
	QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
   	double inlet_haematocrit = 0.45;
   	double dimless_domain_size_x = 2000.0;
    	unsigned dimless_vessel_length = 100.0;
    	double tolerance = 0.001;
    	QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;
    	std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
    	//unsigned NV_Total=6*units_in_x_direction*units_in_y_direction+units_in_y_direction+units_in_x_direction+1;
    	//unsigned NV_ToThinFrom = NV_Total - 2*(3*units_in_y_direction + 1);
    	//unsigned MaxHeterogeneity = (unsigned)(NV_ToThinFrom*0.05);
    	//double percToKill =0.05;
    	//unsigned ToBeKilled = (unsigned)(percToKill*NV_ToThinFrom);
	unsigned NV_ToThinFrom =405;
    	double percToKill =0.25;
    	unsigned ToBeKilled = (unsigned)(percToKill*NV_ToThinFrom);
	//unsigned MaxHeterogeneity = (unsigned)(NV_ToThinFrom*0.05);
	unsigned irregularity = 15;

    	// Generate the network


    	VesselNetworkGenerator<2> network_generator;
	std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/Irregular0."+to_string(irregularity)+"EdgesMatrix.txt");
	std::cout << "/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/Irregular0."+to_string(irregularity)+"EdgesMatrix.txt" << "\n";

   	auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Irregular0."+to_string(irregularity)+"_Hexagons_FromVoronoi", true);

    	std::vector<std::vector<double> > rEdgesMatrix;
    	string line;

    	while (std::getline(in, line)) 
   	{
            	rEdgesMatrix.push_back(std::vector<double>());
            	// Break down the row into column values
         	std::stringstream split(line);
           	double value;
            	while (split >> value)
		{
        	        rEdgesMatrix.back().push_back(value);
        	}
    	}

    	std::shared_ptr<VesselNetwork<2> > p_network;
	std::vector<std::shared_ptr<Vessel<2> > > vessels;


   	QLength large_vessel_radius = 10_um;
    	//QLength small_vessel_radius = 5_um;
	double PerfusionQuotient;


	//unsigned NumberOfInlets = 11;
	unsigned attempts = 1;
	std::vector<std::vector<double> > QuotientMatrixMean(1,std::vector<double>(ToBeKilled+1));
	std::ofstream outfileMean;
    	std::vector<std::vector<double> > QuotientMatrixAggregate(1,std::vector<double>(ToBeKilled+1,0.0));

    	auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
    	p_viscosity_calculator->SetPlasmaViscosity(viscosity);
	auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
	FlowSolver<2> flow_solver;
	flow_solver.SetUseDirectSolver(false);
    	string line2;
	std::vector<std::vector<unsigned> > Order;
	VesselPtr<2> p_thin_vessel;
	
    	for(unsigned att=0; att < attempts;att++)
    	{
    		outfileMean.open("Irregular0."+to_string(irregularity)+"HexagonsFromVoronoiMean.txt");//, std::ios_base::app);
		std::ifstream in2("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/Irregular0."+to_string(irregularity)+"Selection"+to_string(att+1)+".txt");
		Order.clear();
		while (std::getline(in2, line2)) 
   		{
            		Order.push_back(std::vector<unsigned>());
            		// Break down the row into column values
         		std::stringstream split2(line2);
           		unsigned value2;
            		while (split2 >> value2)
			{
        	        	Order.back().push_back(value2);
        		}
    		}
		//std::cout << Order[0][0]<<"  " << Order[1][0] << "  " << Order[192][0] <<"\n";// << Order[1][0]<<"\n";
		

 
		std::string output_file_two_radii = p_file_handler->GetOutputDirectoryFullPath().append("PreRTNetwork"+to_string(irregularity)+".vtp");
    		std::string output_file_RT = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkRT"+to_string(irregularity)+".vtp");

		p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);
		auto p_segment = p_network->GetVesselSegments()[0];
		// set inlet and outlet nodes
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
	    	p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
   		VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
   		//p_network->Write(output_file_two_radii);
		vessels = p_network->GetVessels();

        	// Kill vessels from smallest, up to the specified dose
    		for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
		{
			//vessels = p_network->GetVessels();
			p_viscosity_calculator->SetVesselNetwork(p_network);
    			p_impedance_calculator->SetVesselNetwork(p_network);
    			flow_solver.SetVesselNetwork(p_network);
		    	p_viscosity_calculator->Calculate();
		    	p_impedance_calculator->Calculate();
		    	flow_solver.SetUp();
		    	flow_solver.Solve();
			unsigned PictureNumber = (unsigned)(0.5*ToBeKilled);
			if(KilledVessels ==PictureNumber)
			{
				p_network->Write(output_file_RT);
			}
			if(KilledVessels ==0)
			{
				p_network->Write(output_file_two_radii);
			}
			PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
        		QuotientMatrixAggregate[0][KilledVessels] += PerfusionQuotient;
			p_network->RemoveVessel(vessels[Order[KilledVessels][0]-1],true);
			//p_network->RemoveVessel(p_network->GetVessel(Order[KilledVessels][0]-1),true);
    			
		}
			
    


    		std::cout << "Attempts " << att+1 << " \n";
    		outfileMean << "Attempts " << att+1 << " \n";
    
    		for (unsigned iii = 0; iii < 1; iii++)
    		{
        		for(unsigned jjj = 0; jjj < ToBeKilled+1; jjj++)
        		{
	    			QuotientMatrixMean[iii][jjj]= QuotientMatrixAggregate[iii][jjj]/double(att+1);
	   	 		std::cout << "For heterogeneity level " << iii << ", and radiotherapy dose " << jjj << ", the aggregate quotient is: " << QuotientMatrixAggregate[iii][jjj] << " , and the mean quotient is: " << QuotientMatrixMean[iii][jjj] << endl;
	    			outfileMean << iii << " " << jjj << " " << QuotientMatrixMean[iii][jjj] << " \n";
        		}
    		}    

	outfileMean.close();
	}


//std::cout<< inlet_haematocrit<< viscosity << dimless_domain_size_x<<dimless_vessel_length<<tolerance<<threshold<<ToBeKilled<<irregularity;
}

void xTestNoRadiotherapyHexagonal()
{

    //unsigned NV_Total = 6*4*9;
  /*
    unsigned NV_ToKillFrom = 6*3*7;
    unsigned KilledVessels = 0;

    double dose = 0.0;
    unsigned NV_ToKill = (unsigned)(dose*NV_ToKillFrom);
*/
    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;

    double inlet_haematocrit = 0.45;

    // Generate the network

    VesselNetworkGenerator<2> network_generator;
   


// Specify the domain
    QLength vessel_length = 100_um;
    QLength domain_side_length_x = 800_um;
    QLength domain_side_length_y = 800_um;

    auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Hexagonal_Full_BottomTop", true);
    std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateHexagonalNetworkEquilateral(domain_side_length_x, domain_side_length_y, vessel_length,false);
    std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");
    p_network->Write(output_file_initial);
    std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetwork.vtp");
    std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkRT.vtp");
    QLength real_length_x = domain_side_length_x;
    QLength real_length_y = domain_side_length_y;
    //QLength offset_x = 200_um;
    //QLength offset_y = 200_um;
    // Assign radii
    QLength large_vessel_radius = 20_um;
    auto p_segment = p_network->GetVesselSegments()[0];
    p_segment->SetRadius(large_vessel_radius);
    p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    //p_segment->GetFlowProperties()->SetViscosity(viscosity);
    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);


   VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
          Vertex<2>(-vessel_length, 0.0_um));
   p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
   p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);


   //VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
    //      Vertex<2>(domain_side_length_x, 0.0_um));
VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
        Vertex<2>(domain_side_length_x, domain_side_length_y));
   p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
   p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);



//VesselPtr<2> p_thin_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
 //       Vertex<2>(2.8*real_length_x/4.0,3.0*real_length_y/4.0));
VesselPtr<2> p_thin_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        Vertex<2>(200_um,200_um));
//p_radiated_vessel->SetToDie();
p_thin_vessel->SetRadius(5_um);
//p_thin_vessel->SetRadius(0.000001_um);



    auto p_grid = RegularGrid<2>::Create();
    QLength grid_spacing = 5_um;
    p_grid->SetSpacing(grid_spacing);
    c_vector<unsigned, 3> dimensions;
    dimensions[0] = unsigned((real_length_x)/(grid_spacing))+1; // num x
    dimensions[1] = unsigned((real_length_y)/(grid_spacing))+1; // num_y
    dimensions[2] = 1;
    p_grid->SetDimensions(dimensions);

   
    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
    p_viscosity_calculator->SetVesselNetwork(p_network);
    p_viscosity_calculator->Calculate();
    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
    p_impedance_calculator->SetVesselNetwork(p_network);
    p_impedance_calculator->Calculate();

    FlowSolver<2> flow_solver;
    flow_solver.SetVesselNetwork(p_network);
    flow_solver.SetUp();
    //flow_solver.SetUseDirectSolver(false);
    flow_solver.Solve();


    p_network->Write(output_file_final);

std::cout << "Inlet absolute flow rate is: " <<     p_segment->GetFlowProperties()->GetFlowRate() << "\n";



    p_thin_vessel->SetRadius(0.000000001_um);
   
    auto p_viscosity_calculatorRT = ViscosityCalculator<2>::Create();
    p_viscosity_calculatorRT->SetPlasmaViscosity(viscosity);
    p_viscosity_calculatorRT->SetVesselNetwork(p_network);
    p_viscosity_calculatorRT->Calculate();
    auto p_impedance_calculatorRT = VesselImpedanceCalculator<2>::Create();
    p_impedance_calculatorRT->SetVesselNetwork(p_network);
    p_impedance_calculatorRT->Calculate();

    FlowSolver<2> flow_solverRT;
    flow_solverRT.SetVesselNetwork(p_network);
    flow_solverRT.SetUp();
    //flow_solver.SetUseDirectSolver(false);
    flow_solverRT.Solve();


    p_network->Write(output_file_final_RT);
std::cout << "Inlet absolute flow rate post radiation is: " <<     p_segment->GetFlowProperties()->GetFlowRate() << "\n";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Unknown (Multiple Inputs)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void xTestNoRadiotherapyHexagonalMultipleInputs()
{


    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;

    double inlet_haematocrit = 0.45;

    // Generate the network

    VesselNetworkGenerator<2> network_generator;
   


// Specify the domain
    QLength vessel_length = 100_um;
    QLength domain_side_length_x = 2000_um;
    QLength domain_side_length_y = 2000_um;

    auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Hexagonal_Full_MultipleInlets_Random", true);
    std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateHexagonalNetworkEquilateralMultipleInlets(domain_side_length_x, domain_side_length_y, vessel_length,false);


    std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");
    std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetwork.vtp");
    std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkRT.vtp");
    //QLength offset_x = 200_um;
    //QLength offset_y = 200_um;
    // Assign radii
    QLength large_vessel_radius = 20_um;
    QLength small_vessel_radius = 10_um;
    auto p_segment = p_network->GetVesselSegments()[0];
    p_segment->SetRadius(large_vessel_radius);
    p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    //p_segment->GetFlowProperties()->SetViscosity(viscosity);
    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);


// Multiple inlets and outlets
      QLength unit_height = sqrt(2.0)*vessel_length;
      unsigned units_in_y_direction = floor(domain_side_length_y/unit_height);
      for(unsigned i_aux=0; i_aux < units_in_y_direction+1;i_aux++)
      {
       VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
            Vertex<2>(-vessel_length, double(i_aux)*unit_height));
    //VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
        //    Vertex<2>(domain_side_length_x, domain_side_length_y));
 	VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
           Vertex<2>(domain_side_length_x, double(i_aux)*unit_height));
      p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
      p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
      p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
      p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);

      }
   
    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
    p_viscosity_calculator->SetVesselNetwork(p_network);
    p_viscosity_calculator->Calculate();
    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
    p_impedance_calculator->SetVesselNetwork(p_network);
    p_impedance_calculator->Calculate();

    FlowSolver<2> flow_solver;
    flow_solver.SetVesselNetwork(p_network);
    flow_solver.SetUp();
    //flow_solver.SetUseDirectSolver(false);
    flow_solver.Solve();


    p_network->Write(output_file_initial);



 // Kill from a given rectangle 


    QLength offset_x = vessel_length;

    QLength real_length_x = vessel_length*(9.0+5.0*sqrt(2));
    QLength real_length_y = vessel_length*(28.0/sqrt(2));


// random vessel killing


    //unsigned NV_Total = 6*4*9;
  /*
    unsigned NV_ToKillFrom = 6*3*7;
   

    double dose = 0.0;
    unsigned NV_ToKill = (unsigned)(dose*NV_ToKillFrom);
*/
   unsigned NV_Total=6*5*14+20;
   unsigned NV_ToThinFrom = NV_Total - 86;
   unsigned MaxHeterogeneity = (unsigned)(NV_ToThinFrom*0.2);
   double percToKill =0.05;
   unsigned ToBeKilled = (unsigned)(percToKill*NV_ToThinFrom);

   

std::vector<double> xOrder (MaxHeterogeneity);
std::vector<double> yOrder (MaxHeterogeneity);

unsigned NThin;
unsigned ToBeThin;

srand(unsigned(time(NULL)));

for(unsigned i=0; i<6;i++)
{

    //p_segment->SetRadius(large_vessel_radius);
    //p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    //p_segment->GetFlowProperties()->SetViscosity(viscosity);


    std::ofstream outfile;

    outfile.open("testHexagonal.txt", std::ios_base::app);

    std::ofstream outfile2;

    outfile2.open("optimumHexagonal.txt", std::ios_base::app);

    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   


   double percOfThin = 0.01*i;
   ToBeThin = (unsigned)(percOfThin*NV_ToThinFrom);
   NThin = 0;

   while( NThin < ToBeThin ) {
      double x_norm = (double)rand()/RAND_MAX;
      double y_norm = (double)rand()/RAND_MAX;
      QLength x_Kill = offset_x + x_norm*(real_length_x-2.0*offset_x);
      QLength y_Kill = y_norm*real_length_y;
      //std::cout << x_Kill << y_Kill;
      //VesselPtr<2> p_radiated_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
       //  Vertex<2>(x_Kill,y_Kill));
    VesselPtr<2> p_radiated_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        Vertex<2>(x_Kill,y_Kill));
      if (p_radiated_vessel->GetRadius()==large_vessel_radius) {
    	//p_radiated_vessel->SetToDie();
    	p_radiated_vessel->SetRadius(small_vessel_radius);
        xOrder[NThin]=x_Kill;
	yOrder[NThin]=y_Kill;
	NThin++;
        }
      std::cout << "So far this many thin vessels: "  << NThin << " \n"; // << xOrder << " and y coordinates " << yOrder <<  " \n";
      //std::cout << xOrder[0] << " " << xOrder[1] << " " << xOrder[2] << " " << yOrder[0] << " " << yOrder[1] << " " << yOrder[2] << " \n";
    }


    p_viscosity_calculator->Calculate();
    p_impedance_calculator->Calculate();
    flow_solver.SetUp();
    flow_solver.Solve();

    QFlowRate threshold = 3.e-12*unit::metre_cubed_per_second;


    double PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
    std::cout << "Perfusion quotient of the heterogeneous network before RT is : " << PerfusionQuotient << " \n"; 

    p_network->Write(output_file_final);


    unsigned KilledVessels = 0;

   

    outfile << ToBeThin << " " << KilledVessels << " " << PerfusionQuotient << " \n"; 


    double optimumBeta =PerfusionQuotient;
    unsigned optimumKilled = 0;


    for(unsigned idx=0; idx < ToBeThin; idx++) {
    
     QLength x_Kill_Order = xOrder[idx];
     QLength y_Kill_Order = yOrder[idx];
    VesselPtr<2> p_radiated_vesselOrder = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        Vertex<2>(x_Kill_Order,y_Kill_Order));
      p_radiated_vesselOrder->SetRadius(0.000000000000001_um);
     KilledVessels++;

    p_viscosity_calculator->Calculate();
    p_impedance_calculator->Calculate();
    flow_solver.SetUp();
    flow_solver.Solve();
    std::cout << "So far killed"  << KilledVessels << "\n";
    std::cout << "The perfusion coefficient for " << KilledVessels << " killed vessels is: " << p_network->GetPerfusionQuotientBeta(threshold) << "\n";

    outfile << ToBeThin << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 

    if(p_network->GetPerfusionQuotientBeta(threshold)<optimumBeta)
    {
    optimumBeta=p_network->GetPerfusionQuotientBeta(threshold);
    optimumKilled = KilledVessels;
    }

    //outfile2 << ToBeThin << " " << optimumKilled << " \n";  
    p_network->Write(output_file_final_RT);
    }
    
while(KilledVessels<ToBeKilled)  
{

      double x_norm = (double)rand()/RAND_MAX;
      double y_norm = (double)rand()/RAND_MAX;
      QLength x_Kill = offset_x + x_norm*(real_length_x-2.0*offset_x);
      QLength y_Kill = y_norm*real_length_y;
 VesselPtr<2> p_radiated_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        Vertex<2>(x_Kill,y_Kill));
     
		if(p_radiated_vessel->GetRadius()>0.00001_um)
		{
		p_radiated_vessel->SetRadius(0.000000000000001_um);
     		KilledVessels++;

    		p_viscosity_calculator->Calculate();
    		p_impedance_calculator->Calculate();
    		flow_solver.SetUp();
    		flow_solver.Solve();
    		std::cout << "So far killed"  << KilledVessels << "\n";
    		std::cout << "The perfusion coefficient for " << KilledVessels << " killed vessels is: " << p_network->GetPerfusionQuotientBeta(threshold) << "\n";

                outfile << ToBeThin << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 

    		if(p_network->GetPerfusionQuotientBeta(threshold)<optimumBeta)
    		{
    		optimumBeta=p_network->GetPerfusionQuotientBeta(threshold);
    		optimumKilled = KilledVessels;
    		}

    		p_network->Write(output_file_final_RT);
		}

}

outfile2 << ToBeThin << " " << optimumKilled << " \n";  

outfile.close();
outfile2.close();
}


}

void xTestNoRadiotherapyHexagonalMultipleInputsPolished()
{
    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
    double inlet_haematocrit = 0.45;

    // Generate the network and specify the domain
    VesselNetworkGenerator<2> network_generator;
    QLength vessel_length = 100_um;
    QLength domain_side_length_x = 2000_um;
    QLength domain_side_length_y = 2000_um;

    auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Hexagonal_Full_MultipleInlets_Random_Polished", true);
    std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateHexagonalNetworkEquilateralMultipleInlets(domain_side_length_x, domain_side_length_y, vessel_length,false);
    std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");
    std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkHeterogeneousPreRT.vtp");
    std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkHeterogeneousPostRT.vtp");
    // Assign two distinct radii and inlet haematocrit everywhere
    QLength large_vessel_radius = 10_um;
    QLength small_vessel_radius = 5_um;
    auto p_segment = p_network->GetVesselSegments()[0];
    p_segment->SetRadius(large_vessel_radius);
    p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);
    // Multiple inlets and outlets
    QLength unit_height = sqrt(2.0)*vessel_length;
    QLength unit_width = (2.0+sqrt(2.0))*vessel_length;
    unsigned units_in_x_direction = floor(domain_side_length_x/unit_width);
    unsigned units_in_y_direction = floor(domain_side_length_y/unit_height);

    // Select input and output nodes
    for(unsigned i_aux=0; i_aux < units_in_y_direction+1;i_aux++)
    {
       		VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
            		Vertex<2>(-vessel_length, double(i_aux)*unit_height));
		VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
           		Vertex<2>(domain_side_length_x, double(i_aux)*unit_height));
      		p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
      		p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
      		p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
      		p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);
    }
   
    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
    p_viscosity_calculator->SetVesselNetwork(p_network);
    p_viscosity_calculator->Calculate();
    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
    p_impedance_calculator->SetVesselNetwork(p_network);
    p_impedance_calculator->Calculate();
    FlowSolver<2> flow_solver;
    flow_solver.SetVesselNetwork(p_network);
    flow_solver.SetUp();
    //flow_solver.SetUseDirectSolver(false);
    flow_solver.Solve();
    p_network->Write(output_file_initial);

    // Kill from a given rectangle 
    QLength offset_x = vessel_length/sqrt(2.0);
    QLength real_length_x = double(units_in_x_direction)*unit_width - vessel_length;
    QLength real_length_y = double(units_in_y_direction)*unit_height;

    unsigned NV_Total=6*units_in_x_direction*units_in_y_direction+units_in_y_direction+units_in_x_direction+1;
    unsigned NV_ToThinFrom = NV_Total - 2*(3*units_in_y_direction + 1);
    unsigned MaxHeterogeneity = (unsigned)(NV_ToThinFrom*0.05);
    double percToKill =0.05;
    unsigned ToBeKilled = (unsigned)(percToKill*NV_ToThinFrom);

    std::vector<double> xOrder (MaxHeterogeneity);
    std::vector<double> yOrder (MaxHeterogeneity);
    unsigned NThin;
    unsigned ToBeThin;
    srand(unsigned(time(NULL)));
    // vary the initial radius heterogeneity
    for(unsigned i=0; i<6;i++)
    {
    	std::ofstream outfile;
    	outfile.open("testHexagonal.txt", std::ios_base::app);
    	std::ofstream outfile2;
    	outfile2.open("optimumHexagonal.txt", std::ios_base::app);
   	VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   

   	double percOfThin = 0.01*i;
   	ToBeThin = (unsigned)(percOfThin*NV_ToThinFrom);
   	NThin = 0;
        // randomly pick the given number of vessels and make them thin
   	while( NThin < ToBeThin )
	{
      		double x_norm = (double)rand()/RAND_MAX;
      		double y_norm = (double)rand()/RAND_MAX;
      		QLength x_Kill = offset_x + x_norm*(real_length_x-2.0*offset_x);
      		QLength y_Kill = y_norm*real_length_y;
    		VesselPtr<2> p_radiated_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        		Vertex<2>(x_Kill,y_Kill));
      		if (p_radiated_vessel->GetRadius()==large_vessel_radius)
		{
    			p_radiated_vessel->SetRadius(small_vessel_radius);
        		xOrder[NThin]=x_Kill;
			yOrder[NThin]=y_Kill;
			NThin++;
        	}
      		//std::cout << "So far this many thin vessels: "  << NThin << " \n"; // << xOrder << " and y coordinates " << yOrder <<  " \n";
   	}

    	p_viscosity_calculator->Calculate();
    	p_impedance_calculator->Calculate();
    	flow_solver.SetUp();
    	flow_solver.Solve();
	QFlowRate threshold = 2.e-13*unit::metre_cubed_per_second;
	double PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
    	//std::cout << "Perfusion quotient of the heterogeneous network before RT is : " << PerfusionQuotient << " \n"; 
	p_network->Write(output_file_final);

	unsigned KilledVessels = 0;
	outfile << ToBeThin << " " << KilledVessels << " " << PerfusionQuotient << " \n"; 
	double optimumBeta =PerfusionQuotient;
   	unsigned optimumKilled = 0;

        // Kill vessels from smallest, up to the number of thin vessels
    	for(unsigned idx=0; idx < ToBeThin; idx++)
	{
    		QLength x_Kill_Order = xOrder[idx];
     		QLength y_Kill_Order = yOrder[idx];
    		VesselPtr<2> p_radiated_vesselOrder = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        		Vertex<2>(x_Kill_Order,y_Kill_Order));
      		p_radiated_vesselOrder->SetRadius(0.000000000000001_um);
     		KilledVessels++;

  		p_viscosity_calculator->Calculate();
    		p_impedance_calculator->Calculate();
    		flow_solver.SetUp();
    		flow_solver.Solve();
    		//std::cout << "So far killed"  << KilledVessels << "\n";
    		//std::cout << "The perfusion coefficient for " << KilledVessels << " killed vessels is: " << p_network->GetPerfusionQuotientBeta(threshold) << "\n";
		outfile << ToBeThin << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
		// find optimum value of beta and killed vessels
    		if(p_network->GetPerfusionQuotientBeta(threshold)<optimumBeta)
    		{
    			optimumBeta=p_network->GetPerfusionQuotientBeta(threshold);
    			optimumKilled = KilledVessels;
    		}
		p_network->Write(output_file_final_RT);
    	}
        // kill further vessels, until specified dose is achieved
	while(KilledVessels<ToBeKilled)  
	{
		double x_norm = (double)rand()/RAND_MAX;
      		double y_norm = (double)rand()/RAND_MAX;
      		QLength x_Kill = offset_x + x_norm*(real_length_x-2.0*offset_x);
      		QLength y_Kill = y_norm*real_length_y;
 		VesselPtr<2> p_radiated_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        		Vertex<2>(x_Kill,y_Kill));
                // only kill if that vessel has not been killed already
		if(p_radiated_vessel->GetRadius()>0.00001_um)
		{
			p_radiated_vessel->SetRadius(0.000000000000001_um);
     			KilledVessels++;
			p_viscosity_calculator->Calculate();
	    		p_impedance_calculator->Calculate();
    			flow_solver.SetUp();
    			flow_solver.Solve();
    			//std::cout << "So far killed"  << KilledVessels << "\n";
    			//std::cout << "The perfusion coefficient for " << KilledVessels << " killed vessels is: " << p_network->GetPerfusionQuotientBeta(threshold) << "\n";
		        outfile << ToBeThin << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n";
			//  find optimum value of beta and killed vessels
			if(p_network->GetPerfusionQuotientBeta(threshold)<optimumBeta)
    			{
	    			optimumBeta=p_network->GetPerfusionQuotientBeta(threshold);
    				optimumKilled = KilledVessels;
    			}

    		p_network->Write(output_file_final_RT);
		}

	}

    outfile2 << ToBeThin << " " << optimumKilled << " \n";  
    outfile.close();
    outfile2.close();
    }

}

void xTestNoRadiotherapyHexagonalMultipleInputsPolishedStatistics()
{
    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
    double inlet_haematocrit = 0.45;

    // Generate the network and specify the domain
    VesselNetworkGenerator<2> network_generator;
    QLength vessel_length = 100_um;
    QLength domain_side_length_x = 2000_um;
    QLength domain_side_length_y = 2000_um;

    auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Hexagonal_Full_MultipleInlets_Random_Polished_Statistics", true);
    std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateHexagonalNetworkEquilateralMultipleInlets(domain_side_length_x, domain_side_length_y, vessel_length,false);
    std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");
    std::string output_file_intermediate = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkHeterogeneousPreRTIntermediate.vtp");
    std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkHeterogeneousPreRT.vtp");
    std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkHeterogeneousPostRT.vtp");
    // Assign two distinct radii and inlet haematocrit everywhere
    QLength large_vessel_radius = 10_um;
    QLength small_vessel_radius = 5_um;
    // Multiple inlets and outlets
    QLength unit_height = sqrt(2.0)*vessel_length;
    QLength unit_width = (2.0+sqrt(2.0))*vessel_length;
    unsigned units_in_x_direction = floor(domain_side_length_x/unit_width);
    unsigned units_in_y_direction = floor(domain_side_length_y/unit_height);

    // Kill from a given rectangle 
    QLength offset_x = vessel_length/sqrt(2.0);
    QLength real_length_x = double(units_in_x_direction)*unit_width - vessel_length;
    QLength real_length_y = double(units_in_y_direction)*unit_height;



    unsigned NV_Total=6*units_in_x_direction*units_in_y_direction+units_in_y_direction+units_in_x_direction+1;
    unsigned NV_ToThinFrom = NV_Total - 2*(3*units_in_y_direction + 1);
    unsigned MaxHeterogeneity = (unsigned)(NV_ToThinFrom*0.05);
    double percToKill =0.05;
    unsigned ToBeKilled = (unsigned)(percToKill*NV_ToThinFrom);

    std::vector<double> xOrder (MaxHeterogeneity);
    std::vector<double> yOrder (MaxHeterogeneity);

    // Select input and output nodes
    for(unsigned i_aux=0; i_aux < units_in_y_direction+1;i_aux++)
    {
       		VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
            		Vertex<2>(-vessel_length, double(i_aux)*unit_height));
		VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
           		Vertex<2>(domain_side_length_x, double(i_aux)*unit_height));
      		p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
      		p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
      		p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
      		p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);
    }


    auto p_segment = p_network->GetVesselSegments()[0];
    p_segment->SetRadius(large_vessel_radius);
    p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);
    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
    p_viscosity_calculator->SetVesselNetwork(p_network);
    p_viscosity_calculator->Calculate();
    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
    p_impedance_calculator->SetVesselNetwork(p_network);
    p_impedance_calculator->Calculate();
    FlowSolver<2> flow_solver;
    flow_solver.SetVesselNetwork(p_network);
    flow_solver.SetUp();
    //flow_solver.SetUseDirectSolver(false);
    flow_solver.Solve();
    p_network->Write(output_file_initial);
   
    srand(unsigned(time(NULL)));
    unsigned NThin;
    unsigned ToBeThin;
    std::ofstream outfile;
    outfile.open("testHexagonal.txt");//, std::ios_base::app);
    std::ofstream outfile2;
    outfile2.open("optimumHexagonal.txt");//, std::ios_base::app);
    // vary the initial radius heterogeneity

    unsigned NHet = 6;

    //double QuotientMatrix[NHet][ToBeKilled+1];
    //double QuotientMatrix[6][18];
    //std::vector<std::vector<double>> QuotientMatrix[6][18];
    
    //std::vector<std::vector<double> > QuotientMatrix(NHet,std::vector<double>(ToBeKilled+1));

    unsigned attempts = 20;
    std::vector<std::vector<std::vector<double> > > QuotientMatrixSampling(NHet,std::vector<std::vector<double> >(ToBeKilled+1,std::vector <double>(attempts,0)));

    for(unsigned att=0; att < attempts;att++)
    {


    auto p_segment = p_network->GetVesselSegments()[0];
    p_segment->SetRadius(large_vessel_radius);
    p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);
 

    for(unsigned i=0; i<NHet;i++)
    {
   	VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   

   	double percOfThin = 0.01*i;
   	ToBeThin = (unsigned)(percOfThin*NV_ToThinFrom);
   	NThin = 0;
        // randomly pick the given number of vessels and make them thin
   	while( NThin < ToBeThin )
	{
      		double x_norm = (double)rand()/RAND_MAX;
      		double y_norm = (double)rand()/RAND_MAX;
      		QLength x_Kill = offset_x + x_norm*(real_length_x-2.0*offset_x);
      		QLength y_Kill = y_norm*real_length_y;
    		VesselPtr<2> p_radiated_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        		Vertex<2>(x_Kill,y_Kill));
      		if (p_radiated_vessel->GetRadius()==large_vessel_radius)
		{
    			p_radiated_vessel->SetRadius(small_vessel_radius);
        		xOrder[NThin]=x_Kill;
			yOrder[NThin]=y_Kill;
			NThin++;
        	}
      		//std::cout << "So far this many thin vessels: "  << NThin << " \n"; // << xOrder << " and y coordinates " << yOrder <<  " \n";
   	}
	p_network->Write(output_file_intermediate);
    	p_viscosity_calculator->Calculate();
    	p_impedance_calculator->Calculate();
    	flow_solver.SetUp();
    	flow_solver.Solve();
	QFlowRate threshold = 2.e-13*unit::metre_cubed_per_second;
	double PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
    	//std::cout << "Perfusion quotient of the heterogeneous network before RT is : " << PerfusionQuotient << " \n"; 
	p_network->Write(output_file_intermediate);

	unsigned KilledVessels = 0;
	outfile << ToBeThin << " " << KilledVessels << " " << PerfusionQuotient << " \n";
        QuotientMatrixSampling[i][KilledVessels][att] = p_network->GetPerfusionQuotientBeta(threshold);
	double optimumBeta =PerfusionQuotient;
   	unsigned optimumKilled = 0;

        // Kill vessels from smallest, up to the number of thin vessels
    	for(unsigned idx=0; idx < ToBeThin; idx++)
	{
    		QLength x_Kill_Order = xOrder[idx];
     		QLength y_Kill_Order = yOrder[idx];
    		VesselPtr<2> p_radiated_vesselOrder = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        		Vertex<2>(x_Kill_Order,y_Kill_Order));
      		p_radiated_vesselOrder->SetRadius(0.000000000000001_um);
     		KilledVessels++;

  		p_viscosity_calculator->Calculate();
    		p_impedance_calculator->Calculate();
    		flow_solver.SetUp();
    		flow_solver.Solve();
    		//std::cout << "So far killed"  << KilledVessels << "\n";
    		//std::cout << "The perfusion coefficient for " << KilledVessels << " killed vessels is: " << p_network->GetPerfusionQuotientBeta(threshold) << "\n";
		outfile << ToBeThin << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n";
                QuotientMatrixSampling[i][KilledVessels][att] = p_network->GetPerfusionQuotientBeta(threshold);
		// find optimum value of beta and killed vessels
    		if(p_network->GetPerfusionQuotientBeta(threshold)<optimumBeta)
    		{
    			optimumBeta=p_network->GetPerfusionQuotientBeta(threshold);
    			optimumKilled = KilledVessels;
    		}
		p_network->Write(output_file_final_RT);
		p_network->Write(output_file_intermediate);
    	}
        // kill further vessels, until specified dose is achieved
	while(KilledVessels<ToBeKilled)  
	{
		double x_norm = (double)rand()/RAND_MAX;
      		double y_norm = (double)rand()/RAND_MAX;
      		QLength x_Kill = offset_x + x_norm*(real_length_x-2.0*offset_x);
      		QLength y_Kill = y_norm*real_length_y;
 		VesselPtr<2> p_radiated_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        		Vertex<2>(x_Kill,y_Kill));
                // only kill if that vessel has not been killed already
		if(p_radiated_vessel->GetRadius()>0.00001_um)
		{
			p_radiated_vessel->SetRadius(0.000000000000001_um);
     			KilledVessels++;
			p_viscosity_calculator->Calculate();
	    		p_impedance_calculator->Calculate();
    			flow_solver.SetUp();
    			flow_solver.Solve();
    			//std::cout << "So far killed"  << KilledVessels << "\n";
    			//std::cout << "The perfusion coefficient for " << KilledVessels << " killed vessels is: " << p_network->GetPerfusionQuotientBeta(threshold) << "\n";
		        outfile << ToBeThin << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n";
        		QuotientMatrixSampling[i][KilledVessels][att] = p_network->GetPerfusionQuotientBeta(threshold);
			//  find optimum value of beta and killed vessels
			if(p_network->GetPerfusionQuotientBeta(threshold)<optimumBeta)
    			{
	    			optimumBeta=p_network->GetPerfusionQuotientBeta(threshold);
    				optimumKilled = KilledVessels;
    			}

    		p_network->Write(output_file_final_RT);
		p_network->Write(output_file_intermediate);
		}

	}

        outfile2 << ToBeThin << " " << optimumKilled << " \n";  
    }    
    }

    outfile.close();
    outfile2.close();

    for (unsigned atttt = 0 ; atttt < attempts; atttt++)
    {
    for (unsigned ii = 0; ii < NHet; ii++)
    {
        for(unsigned jj = 0; jj < ToBeKilled+1; jj++)
        {
            std::cout << "At attempts number " << atttt <<  ", heterogeneity level is " << ii << ", radiotherapy dose is " << jj << " and the quotient matrix is: " << QuotientMatrixSampling[ii][jj][atttt] << endl;
        }
    }
    }

    std::vector<std::vector<double> > QuotientMatrixMean(NHet,std::vector<double>(ToBeKilled+1));
    std::vector<std::vector<double> > QuotientMatrixAggregate(NHet,std::vector<double>(ToBeKilled+1));

   std::ofstream outfileMean;
    outfileMean.open("testHexagonalMean.txt", std::ios_base::app);
    std::ofstream outfile2Mean;
    outfile2Mean.open("optimumHexagonalMean.txt", std::ios_base::app);

    for (unsigned iii = 0; iii < NHet; iii++)
    {
        for(unsigned jjj = 0; jjj < ToBeKilled+1; jjj++)
        {
	    QuotientMatrixAggregate[iii][jjj]=0.0;
            for (unsigned attttt = 0 ; attttt < attempts; attttt++)
    	    {
	    QuotientMatrixAggregate[iii][jjj]+=QuotientMatrixSampling[iii][jjj][attttt];
	    }
	    QuotientMatrixMean[iii][jjj]= QuotientMatrixAggregate[iii][jjj]/double(attempts);
	    std::cout << "For heterogeneity level " << iii << ", and radiotherapy dose " << jjj << ", the mean quotient is: " << QuotientMatrixMean[iii][jjj] << endl;
	    outfileMean << iii << " " << jjj << " " << QuotientMatrixMean[iii][jjj] << " \n";
        }
    }

    outfileMean.close();
    outfile2Mean.close();

}

void xTestNoRadiotherapyHexagonalMultipleInputsPolishedStatistics_ConvergenceBetter_Original()
{
    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
    double inlet_haematocrit = 0.45;

    // Generate the network and specify the domain
    VesselNetworkGenerator<2> network_generator;
    QLength vessel_length = 100_um;
    QLength domain_side_length_x = 2000_um;
    QLength domain_side_length_y = 2000_um;

    auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Hexagonal_Full_MultipleInlets_Random_Polished_Statistics_ConvergeBetter100_HigherThreshold", true);
    std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateHexagonalNetworkEquilateralMultipleInlets(domain_side_length_x, domain_side_length_y, vessel_length,false);
    std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");
    std::string output_file_intermediate = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkHeterogeneousPreRTIntermediate.vtp");
    std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkHeterogeneousPreRT.vtp");
    std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkHeterogeneousPostRT.vtp");
    // Assign two distinct radii and inlet haematocrit everywhere
    QLength large_vessel_radius = 10_um;
    QLength small_vessel_radius = 5_um;
    // Multiple inlets and outlets
    QLength unit_height = sqrt(2.0)*vessel_length;
    QLength unit_width = (2.0+sqrt(2.0))*vessel_length;
    unsigned units_in_x_direction = floor(domain_side_length_x/unit_width);
    unsigned units_in_y_direction = floor(domain_side_length_y/unit_height);

    // Kill from a given rectangle 
    QLength offset_x = vessel_length/sqrt(2.0);
    QLength real_length_x = double(units_in_x_direction)*unit_width - vessel_length;
    QLength real_length_y = double(units_in_y_direction)*unit_height;



    unsigned NV_Total=6*units_in_x_direction*units_in_y_direction+units_in_y_direction+units_in_x_direction+1;
    unsigned NV_ToThinFrom = NV_Total - 2*(3*units_in_y_direction + 1);
    unsigned MaxHeterogeneity = (unsigned)(NV_ToThinFrom*0.05);
    double percToKill =0.05;
    unsigned ToBeKilled = (unsigned)(percToKill*NV_ToThinFrom);

    std::vector<double> xOrder (MaxHeterogeneity);
    std::vector<double> yOrder (MaxHeterogeneity);

    // Select input and output nodes
    for(unsigned i_aux=0; i_aux < units_in_y_direction+1;i_aux++)
    {
       		VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
            		Vertex<2>(-vessel_length, double(i_aux)*unit_height));
		VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
           		Vertex<2>(domain_side_length_x, double(i_aux)*unit_height));
      		p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
      		p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
      		p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
      		p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);
    }


    auto p_segment = p_network->GetVesselSegments()[0];
    p_segment->SetRadius(large_vessel_radius);
    p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);
    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
    p_viscosity_calculator->SetVesselNetwork(p_network);
    p_viscosity_calculator->Calculate();
    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
    p_impedance_calculator->SetVesselNetwork(p_network);
    p_impedance_calculator->Calculate();
    FlowSolver<2> flow_solver;
    flow_solver.SetVesselNetwork(p_network);
    flow_solver.SetUp();
    flow_solver.SetUseDirectSolver(false);
    flow_solver.Solve();
    p_network->Write(output_file_initial);
   
    srand(unsigned(time(NULL)));
    unsigned NThin;
    unsigned ToBeThin;
    std::ofstream outfile;
    outfile.open("testHexagonal.txt");//, std::ios_base::app);
    std::ofstream outfile2;
    outfile2.open("optimumHexagonal.txt");//, std::ios_base::app);
    // vary the initial radius heterogeneity

    unsigned NHet = 6;

    unsigned attempts = 10;
    //std::vector<std::vector<std::vector<double> > > QuotientMatrixSampling(NHet,std::vector<std::vector<double> >(ToBeKilled+1,std::vector <double>(attempts,0)));


    std::vector<std::vector<double> > QuotientMatrixMean(NHet,std::vector<double>(ToBeKilled+1));

    std::ofstream outfileMean;
    std::ofstream outfile2Mean;
    outfile2Mean.open("optimumHexagonalMean.txt", std::ios_base::app);
    
    std::vector<std::vector<double> > QuotientMatrixAggregate(NHet,std::vector<double>(ToBeKilled+1,0.0));

    for(unsigned att=0; att < attempts;att++)
    {
    auto p_segment = p_network->GetVesselSegments()[0];
    p_segment->SetRadius(large_vessel_radius);
    p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);
 

    for(unsigned i=0; i<NHet;i++)
    {
   	VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   

   	double percOfThin = 0.01*i;
   	ToBeThin = (unsigned)(percOfThin*NV_ToThinFrom);
   	NThin = 0;
        // randomly pick the given number of vessels and make them thin
   	while( NThin < ToBeThin )
	{
      		double x_norm = (double)rand()/RAND_MAX;
      		double y_norm = (double)rand()/RAND_MAX;
      		QLength x_Kill = offset_x + x_norm*(real_length_x-2.0*offset_x);
      		QLength y_Kill = y_norm*real_length_y;
    		VesselPtr<2> p_radiated_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        		Vertex<2>(x_Kill,y_Kill));
      		if (p_radiated_vessel->GetRadius()==large_vessel_radius)
		{
    			p_radiated_vessel->SetRadius(small_vessel_radius);
        		xOrder[NThin]=x_Kill;
			yOrder[NThin]=y_Kill;
			NThin++;
        	}
      		//std::cout << "So far this many thin vessels: "  << NThin << " \n"; // << xOrder << " and y coordinates " << yOrder <<  " \n";
   	}
	p_network->Write(output_file_intermediate);
    	p_viscosity_calculator->Calculate();
    	p_impedance_calculator->Calculate();
    	flow_solver.SetUp();
    	flow_solver.Solve();
	QFlowRate threshold = 4.e-13*unit::metre_cubed_per_second;
	//	QFlowRate threshold = 2.e-13*unit::metre_cubed_per_second;	
	double PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
    	//std::cout << "Perfusion quotient of the heterogeneous network before RT is : " << PerfusionQuotient << " \n"; 
	p_network->Write(output_file_intermediate);

	unsigned KilledVessels = 0;
	outfile << ToBeThin << " " << KilledVessels << " " << PerfusionQuotient << " \n";
        QuotientMatrixAggregate[i][KilledVessels] += p_network->GetPerfusionQuotientBeta(threshold);
        //std::cout << p_network->GetPerfusionQuotientBeta(threshold) << " \n";
	double optimumBeta =PerfusionQuotient;
   	unsigned optimumKilled = 0;

        // Kill vessels from smallest, up to the number of thin vessels
    	for(unsigned idx=0; idx < ToBeThin; idx++)
	{
    		QLength x_Kill_Order = xOrder[idx];
     		QLength y_Kill_Order = yOrder[idx];
    		VesselPtr<2> p_radiated_vesselOrder = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        		Vertex<2>(x_Kill_Order,y_Kill_Order));
      		p_radiated_vesselOrder->SetRadius(0.000000000000001_um);
     		KilledVessels++;

  		p_viscosity_calculator->Calculate();
    		p_impedance_calculator->Calculate();
    		flow_solver.SetUp();
    		flow_solver.Solve();
    		//std::cout << "So far killed"  << KilledVessels << "\n";
    		//std::cout << "The perfusion coefficient for " << KilledVessels << " killed vessels is: " << p_network->GetPerfusionQuotientBeta(threshold) << "\n";
		outfile << ToBeThin << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n";
                QuotientMatrixAggregate[i][KilledVessels] += p_network->GetPerfusionQuotientBeta(threshold);
                //std::cout << p_network->GetPerfusionQuotientBeta(threshold) << " \n";
		// find optimum value of beta and killed vessels
    		if(p_network->GetPerfusionQuotientBeta(threshold)<optimumBeta)
    		{
    			optimumBeta=p_network->GetPerfusionQuotientBeta(threshold);
    			optimumKilled = KilledVessels;
    		}
		p_network->Write(output_file_final_RT);
		p_network->Write(output_file_intermediate);
    	}
        // kill further vessels, until specified dose is achieved
	while(KilledVessels<ToBeKilled)  
	{
		double x_norm = (double)rand()/RAND_MAX;
      		double y_norm = (double)rand()/RAND_MAX;
      		QLength x_Kill = offset_x + x_norm*(real_length_x-2.0*offset_x);
      		QLength y_Kill = y_norm*real_length_y;
 		VesselPtr<2> p_radiated_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        		Vertex<2>(x_Kill,y_Kill));
                // only kill if that vessel has not been killed already
		if(p_radiated_vessel->GetRadius()>0.00001_um)
		{
			p_radiated_vessel->SetRadius(0.000000000000001_um);
     			KilledVessels++;
			p_viscosity_calculator->Calculate();
	    		p_impedance_calculator->Calculate();
    			flow_solver.SetUp();
    			flow_solver.Solve();
    			//std::cout << "So far killed"  << KilledVessels << "\n";
    			//std::cout << "The perfusion coefficient for " << KilledVessels << " killed vessels is: " << p_network->GetPerfusionQuotientBeta(threshold) << "\n";
		        outfile << ToBeThin << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n";
        		QuotientMatrixAggregate[i][KilledVessels] += p_network->GetPerfusionQuotientBeta(threshold);
                        //std::cout << p_network->GetPerfusionQuotientBeta(threshold) << " \n";
			//  find optimum value of beta and killed vessels
			if(p_network->GetPerfusionQuotientBeta(threshold)<optimumBeta)
    			{
	    			optimumBeta=p_network->GetPerfusionQuotientBeta(threshold);
    				optimumKilled = KilledVessels;
    			}

    		p_network->Write(output_file_final_RT);
		p_network->Write(output_file_intermediate);
		}

	}

        outfile2 << ToBeThin << " " << optimumKilled << " \n";  
    }    


    outfileMean.open("testHexagonalMeanIntermediate100_HigherThreshold.txt");//, std::ios_base::app);
    std::cout << "Attempts " << att+1 << " \n";
    outfileMean << "Attempts " << att+1 << " \n";
    
    for (unsigned iii = 0; iii < NHet; iii++)
    {
        for(unsigned jjj = 0; jjj < ToBeKilled+1; jjj++)
        {
	    QuotientMatrixMean[iii][jjj]= QuotientMatrixAggregate[iii][jjj]/double(att+1);
	    std::cout << "For heterogeneity level " << iii << ", and radiotherapy dose " << jjj << ", the aggregate quotient is: " << QuotientMatrixAggregate[iii][jjj] << " , and the mean quotient is: " << QuotientMatrixMean[iii][jjj] << endl;
	    outfileMean << iii << " " << jjj << " " << QuotientMatrixMean[iii][jjj] << " \n";
        }
    }    

    outfileMean.close();
    }

    outfile.close();
    outfile2.close();

    outfile2Mean.close();
/*
    for (unsigned atttt = 0 ; atttt < attempts; atttt++)
    {
    for (unsigned ii = 0; ii < NHet; ii++)
    {
        for(unsigned jj = 0; jj < ToBeKilled+1; jj++)
        {
            std::cout << "At attempts number " << atttt <<  ", heterogeneity level is " << ii << ", radiotherapy dose is " << jj << " and the quotient matrix is: " << QuotientMatrixSampling[ii][jj][atttt] << endl;
        }
    }
    }
*/


}

void xTestNoRadiotherapyHexagonalMultipleInputsPolishedStatistics_ConvergenceBetter_Playing()
{
    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
    double inlet_haematocrit = 0.45;

    // Generate the network and specify the domain
    VesselNetworkGenerator<2> network_generator;
    QLength vessel_length = 100_um;
    QLength domain_side_length_x = 1200_um;
    QLength domain_side_length_y = 2000_um;

    auto p_file_handler = std::make_shared<OutputFileHandler>("Radiotherapy_Hexagonal_Full_MultipleInlets_Random_Polished_Statistics_ConvergeBetter100_Playing_HigherThreshold", true);
    std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateHexagonalNetworkEquilateralMultipleInlets(domain_side_length_x, domain_side_length_y, vessel_length,false);
    std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");
    std::string output_file_intermediate = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkHeterogeneousPreRTIntermediate.vtp");
    std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkHeterogeneousPreRT.vtp");
    std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetworkHeterogeneousPostRT.vtp");
    // Assign two distinct radii and inlet haematocrit everywhere
    QLength large_vessel_radius = 10_um;
    QLength small_vessel_radius = 5_um;
    // Multiple inlets and outlets
    QLength unit_height = sqrt(2.0)*vessel_length;
    QLength unit_width = (2.0+sqrt(2.0))*vessel_length;
    unsigned units_in_x_direction = floor(domain_side_length_x/unit_width);
    unsigned units_in_y_direction = floor(domain_side_length_y/unit_height);

    // Kill from a given rectangle 
    QLength offset_x = vessel_length*(1.0+1.0/sqrt(2.0));
    QLength real_length_x = double(units_in_x_direction)*unit_width - vessel_length;
    QLength real_length_y = double(units_in_y_direction)*unit_height;



    unsigned NV_Total=6*units_in_x_direction*units_in_y_direction+units_in_y_direction+units_in_x_direction+1;
    unsigned NV_ToThinFrom = NV_Total - 2*(4*units_in_y_direction + 1);
    unsigned MaxHeterogeneity = (unsigned)(NV_ToThinFrom*0.4);
    double percToKill =0.4;
    unsigned ToBeKilled = (unsigned)(percToKill*NV_ToThinFrom);

    std::vector<double> xOrder (MaxHeterogeneity);
    std::vector<double> yOrder (MaxHeterogeneity);

    // Select input and output nodes
    for(unsigned i_aux=0; i_aux < units_in_y_direction+1;i_aux++)
    {
       		VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
            		Vertex<2>(-vessel_length, double(i_aux)*unit_height));
		VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
           		Vertex<2>(domain_side_length_x, double(i_aux)*unit_height));
      		p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
      		p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
      		p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
      		p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);
    }


    auto p_segment = p_network->GetVesselSegments()[0];
    p_segment->SetRadius(large_vessel_radius);
    p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);
    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
    p_viscosity_calculator->SetVesselNetwork(p_network);
    p_viscosity_calculator->Calculate();
    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
    p_impedance_calculator->SetVesselNetwork(p_network);
    p_impedance_calculator->Calculate();
    FlowSolver<2> flow_solver;
    flow_solver.SetVesselNetwork(p_network);
    flow_solver.SetUp();
    flow_solver.SetUseDirectSolver(false);
    flow_solver.Solve();
    p_network->Write(output_file_initial);
   
    srand(unsigned(time(NULL)));
    unsigned NThin;
    unsigned ToBeThin;
    std::ofstream outfile;
    outfile.open("testHexagonal.txt");//, std::ios_base::app);
    std::ofstream outfile2;
    outfile2.open("optimumHexagonal.txt");//, std::ios_base::app);
    // vary the initial radius heterogeneity

    unsigned NHet = 11;

    //double QuotientMatrix[NHet][ToBeKilled+1];
    //double QuotientMatrix[6][18];
    //std::vector<std::vector<double>> QuotientMatrix[6][18];
    
    //std::vector<std::vector<double> > QuotientMatrix(NHet,std::vector<double>(ToBeKilled+1));

    unsigned attempts = 8;
    //std::vector<std::vector<std::vector<double> > > QuotientMatrixSampling(NHet,std::vector<std::vector<double> >(ToBeKilled+1,std::vector <double>(attempts,0)));


    std::vector<std::vector<double> > QuotientMatrixMean(NHet,std::vector<double>(ToBeKilled+1));

    std::ofstream outfileMean;
    std::ofstream outfile2Mean;
    outfile2Mean.open("optimumHexagonalMean.txt", std::ios_base::app);
    
    std::vector<std::vector<double> > QuotientMatrixAggregate(NHet,std::vector<double>(ToBeKilled+1,0.0));

    for(unsigned att=0; att < attempts;att++)
    {
    auto p_segment = p_network->GetVesselSegments()[0];
    p_segment->SetRadius(large_vessel_radius);
    p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);
 

    for(unsigned i=0; i<NHet;i++)
    {
   	VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   

   	double percOfThin = 0.04*i;
   	ToBeThin = (unsigned)(percOfThin*NV_ToThinFrom);
   	NThin = 0;
        // randomly pick the given number of vessels and make them thin
   	while( NThin < ToBeThin )
	{
      		double x_norm = (double)rand()/RAND_MAX;
      		double y_norm = (double)rand()/RAND_MAX;
      		QLength x_Kill = offset_x + x_norm*(real_length_x-2.0*offset_x);
      		QLength y_Kill = y_norm*real_length_y;
    		VesselPtr<2> p_radiated_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        		Vertex<2>(x_Kill,y_Kill));
      		if (p_radiated_vessel->GetRadius()==large_vessel_radius)
		{
    			p_radiated_vessel->SetRadius(small_vessel_radius);
        		xOrder[NThin]=x_Kill;
			yOrder[NThin]=y_Kill;
			NThin++;
        	}
      		//std::cout << "So far this many thin vessels: "  << NThin << " \n"; // << xOrder << " and y coordinates " << yOrder <<  " \n";
   	}
	p_network->Write(output_file_intermediate);
    	p_viscosity_calculator->Calculate();
    	p_impedance_calculator->Calculate();
    	flow_solver.SetUp();
    	flow_solver.Solve();
	QFlowRate threshold =6.e-13*unit::metre_cubed_per_second;
	//QFlowRate threshold = 4.e-13*unit::metre_cubed_per_second;
	double PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
    	std::cout << "Perfusion quotient of the heterogeneous network before RT is : " << PerfusionQuotient << " \n"; 
	p_network->Write(output_file_intermediate);

	unsigned KilledVessels = 0;
	outfile << ToBeThin << " " << KilledVessels << " " << PerfusionQuotient << " \n";
        QuotientMatrixAggregate[i][KilledVessels] += p_network->GetPerfusionQuotientBeta(threshold);
        //std::cout << p_network->GetPerfusionQuotientBeta(threshold) << " \n";
	double optimumBeta =PerfusionQuotient;
   	unsigned optimumKilled = 0;

        // Kill vessels from smallest, up to the number of thin vessels
    	for(unsigned idx=0; idx < ToBeThin; idx++)
	{
    		QLength x_Kill_Order = xOrder[idx];
     		QLength y_Kill_Order = yOrder[idx];
    		VesselPtr<2> p_radiated_vesselOrder = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        		Vertex<2>(x_Kill_Order,y_Kill_Order));
      		p_radiated_vesselOrder->SetRadius(0.000000000000001_um);
     		KilledVessels++;

  		p_viscosity_calculator->Calculate();
    		p_impedance_calculator->Calculate();
    		flow_solver.SetUp();
    		flow_solver.Solve();
    		//std::cout << "So far killed"  << KilledVessels << "\n";
    		//std::cout << "The perfusion coefficient for " << KilledVessels << " killed vessels is: " << p_network->GetPerfusionQuotientBeta(threshold) << "\n";
		outfile << ToBeThin << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n";
                QuotientMatrixAggregate[i][KilledVessels] += p_network->GetPerfusionQuotientBeta(threshold);
                //std::cout << p_network->GetPerfusionQuotientBeta(threshold) << " \n";
		// find optimum value of beta and killed vessels
    		if(p_network->GetPerfusionQuotientBeta(threshold)<optimumBeta)
    		{
    			optimumBeta=p_network->GetPerfusionQuotientBeta(threshold);
    			optimumKilled = KilledVessels;
    		}
		p_network->Write(output_file_final_RT);
		p_network->Write(output_file_intermediate);
    	}
        // kill further vessels, until specified dose is achieved
	while(KilledVessels<ToBeKilled)  
	{
		double x_norm = (double)rand()/RAND_MAX;
      		double y_norm = (double)rand()/RAND_MAX;
      		QLength x_Kill = offset_x + x_norm*(real_length_x-2.0*offset_x);
      		QLength y_Kill = y_norm*real_length_y;
 		VesselPtr<2> p_radiated_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        		Vertex<2>(x_Kill,y_Kill));
                // only kill if that vessel has not been killed already
		if(p_radiated_vessel->GetRadius()>0.00001_um)
		{
			p_radiated_vessel->SetRadius(0.000000000000001_um);
     			KilledVessels++;
			p_viscosity_calculator->Calculate();
	    		p_impedance_calculator->Calculate();
    			flow_solver.SetUp();
    			flow_solver.Solve();
    			//std::cout << "So far killed"  << KilledVessels << "\n";
    			//std::cout << "The perfusion coefficient for " << KilledVessels << " killed vessels is: " << p_network->GetPerfusionQuotientBeta(threshold) << "\n";
		        outfile << ToBeThin << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n";
        		QuotientMatrixAggregate[i][KilledVessels] += p_network->GetPerfusionQuotientBeta(threshold);
                        //std::cout << p_network->GetPerfusionQuotientBeta(threshold) << " \n";
			//  find optimum value of beta and killed vessels
			if(p_network->GetPerfusionQuotientBeta(threshold)<optimumBeta)
    			{
	    			optimumBeta=p_network->GetPerfusionQuotientBeta(threshold);
    				optimumKilled = KilledVessels;
    			}

    		p_network->Write(output_file_final_RT);
		p_network->Write(output_file_intermediate);
		}

	}

        outfile2 << ToBeThin << " " << optimumKilled << " \n";  
    }    


    outfileMean.open("testHexagonalMeanIntermediate100_HigherThreshold.txt");//, std::ios_base::app);
    std::cout << "Attempts " << att+1 << " \n";
    outfileMean << "Attempts " << att+1 << " \n";
    
    for (unsigned iii = 0; iii < NHet; iii++)
    {
        for(unsigned jjj = 0; jjj < ToBeKilled+1; jjj++)
        {
	    QuotientMatrixMean[iii][jjj]= QuotientMatrixAggregate[iii][jjj]/double(att+1);
	    std::cout << "For heterogeneity level " << iii << ", and radiotherapy dose " << jjj << ", the aggregate quotient is: " << QuotientMatrixAggregate[iii][jjj] << " , and the mean quotient is: " << QuotientMatrixMean[iii][jjj] << endl;
	    outfileMean << iii << " " << jjj << " " << QuotientMatrixMean[iii][jjj] << " \n";
        }
    }    

    outfileMean.close();
    }

    outfile.close();
    outfile2.close();

    outfile2Mean.close();
/*
    for (unsigned atttt = 0 ; atttt < attempts; atttt++)
    {
    for (unsigned ii = 0; ii < NHet; ii++)
    {
        for(unsigned jj = 0; jj < ToBeKilled+1; jj++)
        {
            std::cout << "At attempts number " << atttt <<  ", heterogeneity level is " << ii << ", radiotherapy dose is " << jj << " and the quotient matrix is: " << QuotientMatrixSampling[ii][jj][atttt] << endl;
        }
    }
    }
*/


}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Unknown (seems to be a testbed)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*

    p_thin_vessel->SetRadius(0.000000001_um);
   
    auto p_viscosity_calculatorRT = ViscosityCalculator<2>::Create();
    p_viscosity_calculatorRT->SetPlasmaViscosity(viscosity);
    p_viscosity_calculatorRT->SetVesselNetwork(p_network);
    p_viscosity_calculatorRT->Calculate();
    auto p_impedance_calculatorRT = VesselImpedanceCalculator<2>::Create();
    p_impedance_calculatorRT->SetVesselNetwork(p_network);
    p_impedance_calculatorRT->Calculate();

    FlowSolver<2> flow_solverRT;
    flow_solverRT.SetVesselNetwork(p_network);
    flow_solverRT.SetUp();
    //flow_solver.SetUseDirectSolver(false);
    flow_solverRT.Solve();


    p_network->Write(output_file_final_RT);
std::cout << "Inlet absolute flow rate post radiation is: " <<     p_segment->GetFlowProperties()->GetFlowRate() << "\n";

*/

//Testing random generator
//srand(unsigned(time(NULL)));
//for(unsigned ii=0;ii<10;ii++)
//{
//	int flip=  rand() % 2;// assign random numbers
//	std::cout<< "ii is: " << ii << ", and flip is: " <<flip << "  ";
//}


 // Kill from a given rectangle 
//QLength domain_side_length_x_ToKill = dimless_length_ToKill*lambda*max_radius;
  //  QLength domain_side_length_x_ToKill_SecondEndPoint = domain_side_length_x - domain_side_length_x_ToKill;

   // QLength LeftEnd_x = domain_side_length_x_ToKill;
  //  QLength RightEnd_x = domain_side_length_x_ToKill_SecondEndPoint;
  //  QLength BottomEnd_y = 0.0*domain_side_length_y;
  //  QLength TopEnd_y = domain_side_length_y;
   // std::cout << "Here, Domain length is:" << domain_side_length_x << "meters \n";
 //   std::cout << "Left point is: " << LeftEnd_x << ", Right point is: " << RightEnd_x << " Bottom end is: " << BottomEnd_y << " Top end is: " << TopEnd_y << "\n";



// Multiple inlets and outlets
  //  QLength unit_height = sqrt(2.0)*vessel_length;
//    unsigned units_in_y_direction = floor(domain_side_length_y/unit_height);
  //  for(unsigned i_aux=0; i_aux < units_in_y_direction+1;i_aux++)
//    {
  //  VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
    //        Vertex<2>(-vessel_length, double(i_aux)*unit_height));
    //VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
        //    Vertex<2>(domain_side_length_x, domain_side_length_y));
//	VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
  //        Vertex<2>(domain_side_length_x, double(i_aux)*unit_height));
 //   p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
  //  p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
    //p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
    //p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);

    //}


// random vessel killing
 /*  srand(unsigned(time(NULL)));

    while( KilledVessels < NV_ToKill ) {
      double x_norm = (double)rand()/RAND_MAX;
      double y_norm = (double)rand()/RAND_MAX;
      QLength x_Kill = offset_x + x_norm*(real_length_x-2.0*offset_x);
      QLength y_Kill = offset_y + y_norm*(real_length_y-2.0*offset_y);
      std::cout << x_Kill << y_Kill;
      //VesselPtr<2> p_radiated_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
       //  Vertex<2>(x_Kill,y_Kill));
    VesselPtr<2> p_radiated_vessel = VesselNetworkGeometryCalculator<2>::GetNearestVessel(p_network,
        Vertex<2>(x_Kill,y_Kill));
      if (p_radiated_vessel->GetIsAlive()==true && p_radiated_vessel->GetId() >24) {
    	p_radiated_vessel->SetToDie();
    	p_radiated_vessel->SetRadius(0.00000000000001_um);
	KilledVessels++;
        }
      std::cout << "Tried again and so far killed"  << KilledVessels << "\n";
    }

*/
//    std::cout << " x to kill: " << x_Kill << " and y to kill: " << y_Kill; //<< " X1 normed is: " << x_norm_1 << " and Y1 normed is: " << y_norm_1 << " X2 normed is: " << x_norm_2 << " and Y2 normed is: " << y_norm_2;
//std::cout << "Is it alive? " << p_radiated_vessel->GetIsAlive();

// oxygen stuff
/*


    auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
    p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
    p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

 
    auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
    QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
            GenericParameters::mpGasConcentrationAtStp->GetValue("User");
    QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
            Owen11Parameters::mpReferencePartialPressure->GetValue("User");
    p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
    p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
    p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
    p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);


    auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
    p_oxygen_solver->SetPde(p_oxygen_pde);
    p_oxygen_solver->SetLabel("oxygen");
    p_oxygen_solver->SetGrid(p_grid);
*/

};

#endif // TESTNOHAEMATOCRITRADIOTHERAPY_HPP
