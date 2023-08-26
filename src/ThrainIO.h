// **************************************************************************************
// GRPulseIO.h
//	This defines the protocols for handling use input and program output
//	Will take a filename form the commandline as an input file
//**************************************************************************************

#ifndef GRPULSEIOH
#define GRPULSEIOH

namespace Calculation {
    struct InputData;
    struct OutputData;
}

namespace io {

//will read user input from the specified file and organize the input data
int read_input(const char input_file_name[128], Calculation::InputData&);
// individual parts of input file
int read_calcname(FILE* fp, Calculation::InputData&);
int read_model(FILE* fp, Calculation::InputData&);
int read_units(FILE* fp, Calculation::InputData&);
int read_frequencies(FILE* fp, Calculation::InputData&);

//will re-print the user input in order to re-run the same calculation
int echo_input(Calculation::InputData &);

//prepare to write output
int setup_output(Calculation::InputData&, Calculation::OutputData&);

//just prints a boc around certain parts of the output file
void print_splash(FILE *fp, char*, int WIDTH);

//write output information about the star itself
int write_stellar_output(Calculation::OutputData&);

//write output data about the calculated modes
int write_mode_output(Calculation::OutputData&);

//write both stellar and mode output
int write_output(Calculation::InputData&);

//write tidal overlap coefficients in separate file
int write_tidal_overlap(Calculation::OutputData&);

} // namespace io

#endif