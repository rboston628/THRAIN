//**************************************************************************************
//							THRAIN
//	Main program.  All other functionality is called from here.
//	This program is not necessary to use the stellar models or the mode functions
//		defined in the src files, but provides a single and simple access to them
//  Reece Boston, Mar 24, 2022
//**************************************************************************************

#include "ThrainMain.h"
#include "ThrainIO.h"

namespace mode {
	int create_modes(Calculation::OutputData &data_out);
}

int main(int argc, char* argv[]){
	//Handle command-line information to extract the input file
	FILE* input_file;
	Calculation::InputData calcdataIn;
	Calculation::OutputData calcdataOut;
	if(argc == 2) {
		//the calculation filename is sent as command-line argument
		//call the read_input routine on the given filename (see GRPulseIO.h)
		if(!io::read_input(argv[1], calcdataIn)) printf("file read\n");
		else return 1;
	}
	else {
		printf("Usage is thrain.out <input file>\n");
		return 1;
	}
	
	//remove any information from past calculations, and create directory for the calculation
	system( ( "rm -r ./output/"+calcdataIn.calcname ).c_str() );
	system( ( "mkdir -p ./output/"+calcdataIn.calcname+"/star" ).c_str() );
	system( ( "mkdir -p ./output/"+calcdataIn.calcname+"/modes" ).c_str() );
	//print a copy of the input file for future reference
	io::echo_input(calcdataIn);
	//prepare the output, based on the input
	io::setup_output(calcdataIn, calcdataOut);
	//make the background stellar model
	create_star(calcdataOut);
	//write output on the bakground stellar model
	io::write_stellar_output(calcdataOut);
	
	mode::create_modes(calcdataOut);
	io::write_mode_output(calcdataOut);
	if(calcdataOut.regime==regime::PN0 && calcdataOut.modetype==modetype::nonradial)
		io::write_tidal_overlap(calcdataOut);
	printf("THRAIN done\n");
	return 0;
}